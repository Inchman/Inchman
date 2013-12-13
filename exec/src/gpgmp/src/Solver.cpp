#include "Solver.h"
#include "DiffusionModel.h"
#include "Species.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/log/trivial.hpp>


namespace gpgmp {
Solver::Solver(DiffusionModel &diffusionModel,
               std::vector<int> *runState,
               std::vector<Real> *deterministicState,
               SolverType solverType,
               bool deterministic, bool hasIndividualSpecies)
    : m_diffusionModel(diffusionModel),
      m_runState(runState),
      m_runDeterministicState(deterministicState),
      m_solverType(solverType),
      m_deterministic(deterministic),
      m_hasIndividualSpecies(hasIndividualSpecies),
      m_hasLocalizedReactions(diffusionModel.hasLocalizedReactions()),
      m_globalRange(diffusionModel.gridWidth(), diffusionModel.gridHeight()),
      m_localRange(LOCAL_RANGE_BLOCK,LOCAL_RANGE_BLOCK),
      m_globalRangeReduce(diffusionModel.gridArea()/(LOCAL_RANGE_BLOCK*LOCAL_RANGE_BLOCK)),
      m_localSizeReduce(diffusionModel.gridArea()/(LOCAL_RANGE_BLOCK*LOCAL_RANGE_BLOCK)*sizeof(Real)),
      m_individualSpecies(diffusionModel.numSpecies(), false),
      d_y0(0), d_ys(0), d_y1(0), d_ym1(0), d_ym2(0), d_scrarray(0), d_qs(0), d_rtaus(0),
      m_errors(GMP_NUM_ERRORS*diffusionModel.gridArea(), 0)
{
    // initialize CL device
    oclInit();

    // create state buffer
    if (deterministic)
        d_state = oclCreateBufferFrom(*m_runDeterministicState);
    else
        d_state = oclCreateBufferFrom(*m_runState);

    // create individual solver object if needed
    if (hasIndividualSpecies) {
        // create individual solver
        m_individualSolver = new IndividualSolver(m_diffusionModel, *this, m_context, m_commandQueue, c_contextDevice);

        // create bool array to quickly work out which species are individual-based
        for (unsigned int tk=0; tk<m_diffusionModel.numSpecies(); tk++)
            m_individualSpecies[tk] = m_diffusionModel.species()[tk]->hasIndividualProperties();
    }

    // the AQSS solver needs a bunch of helpers
    if (m_solverType == gpgmp::deterministic_inhomogeneous_aqss) {
        size_t n = m_diffusionModel.gridArea() * m_diffusionModel.numSpecies();
        d_y0       = oclCreateBuffer<Real>(n);
        d_ys       = oclCreateBuffer<Real>(n);
        d_y1       = oclCreateBuffer<Real>(n);
        d_ym1      = oclCreateBuffer<Real>(n);
        d_ym2      = oclCreateBuffer<Real>(n);
        d_scrarray = oclCreateBuffer<Real>(n);
        d_qs       = oclCreateBuffer<Real>(n);
        d_rtaus    = oclCreateBuffer<Real>(n);
    }

    // Create localized reaction mask if there's any localized reactions in the system
    if (m_hasLocalizedReactions) {
        BOOST_LOG_TRIVIAL(debug) <<"Creating localization mask.."<<std::flush;
        d_reactionMask = oclCreateReactionMask();
    }

    // create error checking array
    d_errors = oclCreateBufferFrom(m_errors);
    // we also need to zero out the device memory for the errors
    for (uint i=0; i<GMP_NUM_ERRORS; i++) m_errors[i]=0;
    oclCopy(m_errors, *d_errors);

    // create buffer for field arrays
    const std::map <std::string, std::vector<Real> > fieldParameters = m_diffusionModel.getFieldParameters();
    size_t fpsize = 0;
    m_hasFields = false;
    if (fieldParameters.size() > 0) {
        m_hasFields = true;
        for (std::map<std::string, std::vector<Real> >::const_iterator itt = fieldParameters.begin();
             itt != fieldParameters.end();
             ++itt) {

            // store offset
            m_fieldOffsets[itt->first] = fpsize;
            fpsize += (itt->second).size();

        }
        BOOST_LOG_TRIVIAL(debug) <<"Total elements in field parameter array: " << fpsize <<" corresponding to "<<fpsize*sizeof(Real)/1024 <<" kB of device memory.";

        // and create the buffer
        d_fields = oclCreateBuffer<Real>(fpsize);

        // and copy it
        synchronizeFieldBufferToDevice();
    }

    // We might need the reaction kernels
    if (m_diffusionModel.numReactions()>0 || (m_hasFields && m_diffusionModel.returnSetFieldMethod().size()>0)) {
        if (!deterministic) {
            BOOST_LOG_TRIVIAL(debug) <<"Creating Gillespie kernel."<<std::flush;

            // create Gillespie kernel
            std::ostringstream gss;

            // Generate CL code and build the program
            oclGenerateHeader(gss);
            oclReadAndReplaceTemplateFile(gss, "Gillespie.cl");
            m_gillespieProgram = oclBuildProgram(gss);

            if (m_diffusionModel.numReactions()>0) {
                // and create the Gillespie kernel
                gillespie = Kernel(m_gillespieProgram,
                                   "performGillespie",
                                   cl::EnqueueArgs(*m_commandQueue, m_globalRange));
                gillespie.setLocalSize(m_localRange);
                if (!m_hasLocalizedReactions)
                    if (!m_debugKernel)
                        gillespie.instance().setArg(3, LOCAL_RANGE_BLOCK*LOCAL_RANGE_BLOCK*sizeof(cl_int)*m_diffusionModel.numSpecies(),0);
                    else {
                        gillespie.instance().setArg(4, LOCAL_RANGE_BLOCK*LOCAL_RANGE_BLOCK*sizeof(cl_int)*m_diffusionModel.numSpecies(),0);
                    }
                else
                    if (!m_debugKernel)
                        gillespie.instance().setArg(4, LOCAL_RANGE_BLOCK*LOCAL_RANGE_BLOCK*sizeof(cl_int)*m_diffusionModel.numSpecies(),0);
                    else {
                        gillespie.instance().setArg(5, LOCAL_RANGE_BLOCK*LOCAL_RANGE_BLOCK*sizeof(cl_int)*m_diffusionModel.numSpecies(),0);
                    }
            }

            if (m_hasFields && m_diffusionModel.returnSetFieldMethod().size()>0) {
                // create the field update kernel
                m_updateFields = Kernel(m_gillespieProgram,
                                   "updateFields",
                                   cl::EnqueueArgs(*m_commandQueue, m_globalRange));
                m_updateFields.setLocalSize(m_localRange);
            }
        }
    } else {
        // create deterministic reaction kernels
        // todo:: put these in a common "reactions.cl" file ..
        std::ostringstream dss;

        // Generate CL code and build the program
        oclGenerateHeader(dss);
        oclReadAndReplaceTemplateFile(dss, "Deterministic.cl");
        m_deterministicProgram = oclBuildProgram(dss);

        if (m_solverType==gpgmp::deterministic_homogeneous_RK4) {
            deterministic_performReaction_rk4 = Kernel(m_deterministicProgram,
                                                       "deterministic_performReaction_rk4",
                                                       cl::EnqueueArgs(*m_commandQueue, m_globalRange));
            deterministic_performReaction_rk4.setLocalSize(m_localRange);
            deterministic_performReaction_rk4.setLocalMemory(LOCAL_RANGE_BLOCK*LOCAL_RANGE_BLOCK*sizeof(Real)*5*m_diffusionModel.species().size());
        }
        else if (m_solverType==gpgmp::deterministic_homogeneous_aqss) {
            deterministic_performReaction_aqss = Kernel(m_deterministicProgram,
                                                        "deterministic_performReaction_aqss",
                                                        cl::EnqueueArgs(*m_commandQueue, m_globalRange));
            deterministic_performReaction_aqss.setLocalSize(m_localRange);
            deterministic_performReaction_aqss.setLocalMemory(LOCAL_RANGE_BLOCK*LOCAL_RANGE_BLOCK*sizeof(Real)*(3*m_diffusionModel.species().size()+1));
        }
    } // create reaction kernels
}

Solver::~Solver() {
    // deleting buffers

    if (m_solverType == gpgmp::deterministic_inhomogeneous_aqss) {
        delete d_y0;
        delete d_ys;
        delete d_y1;
        delete d_ym1;
        delete d_ym2;
        delete d_scrarray;
        delete d_qs;
        delete d_rtaus;
    }

    if (m_hasLocalizedReactions)
        clReleaseMemObject(d_reactionMask);

    if (m_hasFields)
        delete d_fields;

    delete d_errors;
    delete d_state;

    // delete individual solver
    std::cout <<"About to destroy individual solver.."<<std::flush;
    if (m_hasIndividualSpecies)
        delete m_individualSolver;
    std::cout <<"done destroying individual solver."<<std::flush;

    // delete command queue etc.
    oclDestroy();
}

void Solver::doReactions(Real deltat, int &seed)
{
    if (!m_deterministic) {
        // todo: debug code - remove me!

        /*
        if (m_hasIndividualSpecies) {
            // go through all individual species and check for consistency
            for (uint ii=0; ii<m_diffusionModel.numSpecies(); ii++) {
                if (m_individualSpecies[ii])
                    m_individualSolver->checkPositionListConsistency(ii, d_state);
            }
        }
        */
        // end remove me

        // do gillespie
        if (m_hasLocalizedReactions)
            if (!m_debugKernel)
                gillespie(*d_state, seed, deltat, d_reactionMask);
            else
                gillespie(*d_state, seed, deltat, *d_errors, d_reactionMask);
        else
            if (!m_debugKernel)
                gillespie(*d_state, seed, deltat);
            else
                gillespie(*d_state, seed, deltat, *d_errors);

        // check for errors
        if (m_debugKernel) checkKernelError();

        // increase RNG seed
        seed++;

        // update individual species
        if (m_hasIndividualSpecies  && !m_diffusionModel.hasConstantIndividuals()) {
            // go through all individual species and update the population
            for (uint ii=0; ii<m_diffusionModel.numSpecies(); ii++) {
                if (m_individualSpecies[ii])
                    m_individualSolver->maintainPopulation(ii, seed++, d_state);
            }
        }

    } else if (m_deterministic) {
        if (m_solverType == gpgmp::deterministic_homogeneous_RK4)
            if (m_hasLocalizedReactions) {
                deterministic_performReaction_rk4(*d_state, deltat, d_reactionMask);
            } else {
                deterministic_performReaction_rk4(*d_state, deltat);
            }
        else if (m_solverType == gpgmp::deterministic_homogeneous_aqss) {
            if (m_hasLocalizedReactions) {
                deterministic_performReaction_aqss(*d_state, deltat, *d_y0, *d_ys, *d_y1, *d_ym1, *d_ym2,
                                                   *d_scrarray, *d_qs, *d_rtaus, d_reactionMask);
            } else {
                deterministic_performReaction_aqss(*d_state, deltat, *d_y0, *d_ys, *d_y1, *d_ym1, *d_ym2,
                                                   *d_scrarray, *d_qs, *d_rtaus);
            }
        }
    }
}

void Solver::updateFields(Real deltat, Real simtime)
{
    m_updateFields(deltat, simtime, *d_fields, *d_state);
}

bool Solver::checkKernelError() {
    // TODO: We could put this one in the Kernel class ..

    // flush command queue and copy error array
    m_commandQueue->flush();
    oclCopy(*d_errors, m_errors);

    // we need to check for all threads .. super heavy-weight!
    bool ret = false;
    for (uint ii=0; ii<m_diffusionModel.gridArea(); ii++) {
        if (m_errors[ii*GMP_NUM_ERRORS]>0) {
        ret = true;
        BOOST_LOG_TRIVIAL(error) <<"Kernel error detected in thread "<<ii<<". Dumping error array ..\n";
        for (uint i=0; i<GMP_NUM_ERRORS; i++)
            std::cerr <<"error["<<i<<"] = "<<m_errors[ii*GMP_NUM_ERRORS+i]<<" ";
        BOOST_LOG_TRIVIAL(error) << std::flush;
        ret = true;
        }
    }

    return ret;
}

void Solver::synchronizeStateBufferFromDevice()
{
    m_commandQueue->flush();
    if (!m_deterministic)
        oclCopy(*d_state, *m_runState);
    else
        oclCopy(*d_state, *m_runDeterministicState);
}

void Solver::synchronizeStateBufferToDevice()
{
    m_commandQueue->flush();

    // copy state array back to device
    if (!m_deterministic)
        oclCopy(*m_runState, *d_state);
    else
        oclCopy(*m_runDeterministicState, *d_state);
}

void Solver::synchronizeFieldBufferFromDevice()
{
    // copy state buffer from device
    for (std::map<std::string, std::vector<Real> >::iterator itt = m_diffusionModel.getFieldParameters().begin();
         itt != m_diffusionModel.getFieldParameters().end();
         ++itt) {
        BOOST_LOG_TRIVIAL(trace) <<"Copying field from device"<<itt->first<<". Offsets is " << m_fieldOffsets[itt->first] <<", size is "<< (itt->second).size();

        m_commandQueue->enqueueReadBuffer(*d_fields, CL_TRUE, m_fieldOffsets[itt->first]*sizeof(Real),
                    (itt->second).size()*sizeof(Real),
                    (void *) &((itt->second)[0]) );
            m_commandQueue->flush();

            const Real *ptr = &((itt->second)[0]);
    }

}

// copy state buffer to device
void Solver::synchronizeFieldBufferToDevice()
{
    for (std::map<std::string, std::vector<Real> >::const_iterator itt = m_diffusionModel.getFieldParameters().begin();
         itt != m_diffusionModel.getFieldParameters().end();
         ++itt) {
        BOOST_LOG_TRIVIAL(trace) <<"Copying field "<<itt->first<<". Offsets is " << m_fieldOffsets[itt->first] <<", size is "<< (itt->second).size();

            m_commandQueue->enqueueWriteBuffer(*d_fields, CL_TRUE, m_fieldOffsets[itt->first]*sizeof(Real),
                    (itt->second).size()*sizeof(Real),
                    (const void *) &((itt->second)[0]) );
            m_commandQueue->flush();

            const Real *ptr = &((itt->second)[0]);
    }
}

cl_mem Solver::oclCreateReactionMask() const
{
    if (m_diffusionModel.numReactions() <= 0)
        return 0;

    // prepare the mask array on host
    std::vector<float> h_reactionMask(m_diffusionModel.gridArea() * m_diffusionModel.numReactions());
    float *maskIt = &h_reactionMask.front();
    for (size_t r=0; r<m_diffusionModel.numReactions(); r++) {
        for (size_t y=0; y<m_diffusionModel.gridHeight(); y++) {
            for (size_t x=0; x<m_diffusionModel.gridWidth(); x++) {
                *(maskIt++) = m_diffusionModel.reactions().at(r)->reactionMask(x, y, (size_t)0);
            }
        }
    }

/*
    maskIt = &h_reactionMask.front();
     for (size_t r=0; r<m_diffusionModel.numReactions(); r++) {
         for (size_t y=0; y<m_diffusionModel.gridHeight(); y++) {
        for (size_t x=0; x<m_diffusionModel.gridWidth(); x++) {
            std::cout <<*(maskIt);
        }
        std::cout <<"\n";
         }
         std::cout <<"\n";
    }
*/
    // copy mask to OpenCL texture
    cl_mem_flags flags = CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR;
    cl::ImageFormat format(CL_R, CL_FLOAT);
    /*cl::Image *image = (m_diffusionModel.numReactions() == 1
                        ? (cl::Image *)new cl::Image2D(*m_context, flags, format, m_diffusionModel.gridWidth(), m_diffusionModel.gridHeight(), 0, &h_reactionMask[0])
                        : (cl::Image *)new cl::Image3D(*m_context, flags, format, m_diffusionModel.gridWidth(), m_diffusionModel.gridHeight(), m_diffusionModel.numReactions(), 0, 0, &h_reactionMask[0])
                          );
                          */
    cl::Image *image = (cl::Image *)new cl::Image3D(*m_context, flags, format, m_diffusionModel.gridWidth(), m_diffusionModel.gridHeight(), m_diffusionModel.numReactions(), 0, 0, &h_reactionMask[0]);

    cl_mem mem = (cl_mem)(*image)();
    clRetainMemObject(mem); // retain, as cl::Image's dtor will call to release it
    return mem; // return opencl's nice dimension independent pointer (you can't setArg with cl::Image)
}

void Solver::oclInit()
{
    // Create context
    m_context = new cl::Context(CL_DEVICE_TYPE_GPU, 0, oclNotification, this);

    // Find devices
    std::vector<cl::Device> devices;
    m_context->getInfo(CL_CONTEXT_DEVICES, &devices);

    // List devices
    BOOST_LOG_TRIVIAL(debug) << "Found OpenCL " << devices.size() << " GPU Device(s)";
    for (size_t i=0; i < devices.size(); i++) {
        std::string name;
        devices[i].getInfo(CL_DEVICE_NAME, &name);
        BOOST_LOG_TRIVIAL(debug) << "  " << i << ": \"" << name << "\"";

        // query the device
        cl_ulong globalMemSize;
        devices[i].getInfo<cl_ulong>(CL_DEVICE_GLOBAL_MEM_SIZE, &globalMemSize);
        BOOST_LOG_TRIVIAL(debug) << " Device " << i << ": \"" << name << "\".globalMemSize:"<<globalMemSize;
        VECTOR_CLASS<std::size_t> work_items;
        devices[i].getInfo< VECTOR_CLASS<std::size_t> >(CL_DEVICE_MAX_WORK_ITEM_SIZES, &work_items);
        BOOST_LOG_TRIVIAL(debug) << " Device " << i << ": \"" << name << "\".maxWorkItemSize:["<<work_items[0]<<", "<<work_items[1]<<", "<<work_items[2]<<"]";
    }

    // Try to use the specified device
    // TODO: this seems to use both available tesla cards on the cluster.. nvidia-smi shows this thread to grab both GPUs .. not sure why!
    if (m_diffusionModel.oclDeviceIndex() > devices.size()-1) {
        BOOST_LOG_TRIVIAL(error) << "Device number out of range: " << m_diffusionModel.oclDeviceIndex() << ", exiting..." ;
        // TODO: Throw exception!
        exit(1);
    }
    c_contextDevice = devices[m_diffusionModel.oclDeviceIndex()];
    std::string deviceName = c_contextDevice.getInfo<CL_DEVICE_NAME>();
    BOOST_LOG_TRIVIAL(debug) << "Starting on OpenCL GPU Device " << m_diffusionModel.oclDeviceIndex() << ": \"" << deviceName << "\"" ;

    // Create command queue
    m_commandQueue = new cl::CommandQueue(*m_context, c_contextDevice);
}

void Solver::oclDestroy()
{
    delete m_commandQueue;
    delete m_context;
}


void CL_CALLBACK Solver::oclNotification(const char *errinfo, const void * /*private_info*/,
                                 size_t /*cb*/, void * /*user_data*/)
{
    // TODO: Throw exception..
    BOOST_LOG_TRIVIAL(error) << errinfo << std::endl;
}

void Solver::oclGenerateHeader(std::ostringstream &ss) const
{
    ss << std::scientific;

    // Preliminaries
    ss << "// Header\n";
    ss << "typedef float Real;\n";
    ss << "\n";

    // Includes for RNG
    ss << "// Includes for the threefry RNG\n";
#if defined(__APPLE__) || defined(__MACOSX)
    ss << "#pragma clang diagnostic push\n";
    ss << "#pragma clang diagnostic ignored \"-Wmissing-prototypes\"\n";
#endif
    ss << "#include <Random123/threefry.h>\n";
#if defined(__APPLE__) || defined(__MACOSX)
    ss << "#pragma clang diagnostic pop\n";
#endif
    ss << "\n";
    ss << "// this is from www.doornik.com/research/randomdouble.pdf\n";
    ss << "// multiply by 2^-32 to get float value\n";
    ss << "#define M_RAN_INVM32 2.32830643653869628906e-010\n";
    ss << "\n";

    // Memory mapping
    ss << "// Mapping from thread ids to global memory\n";
    ss << "#define get_global_area_2d (GridModelWidth * GridModelHeight)\n";
    ss << "#define get_local_area_2d  (get_local_size(0)  * get_local_size(1))\n";
    ss << "#define get_global_id_2d   (get_global_id(0) + (get_global_id(1) * GridModelWidth))\n";
    ss << "#define get_local_id_2d    (get_local_id(0)  + (get_local_id(1)  * get_local_size(0)))\n";
    ss << "#define get_block_id_2d    (get_group_id(0)  + (get_group_id(1)  * get_num_groups(0)))\n";
    ss << "#define get_global_id_for_species_2d(SPECIES_INDEX) (get_global_id_2d + (SPECIES_INDEX*get_global_area_2d))\n";
    ss << "#define get_local_id_for_species_2d(SPECIES_INDEX) (get_local_id_2d + (SPECIES_INDEX*get_local_area_2d))\n";
    ss << "\n";

    // Localization support
    if (m_diffusionModel.hasLocalizedReactions()) {
        ss << "#define GPGMP_HAS_LOCALIZED_REACTIONS\n\n";
        ss << "// Sampler for the reaction mask\n";
        ss << "const sampler_t _reactionMaskSampler = CLK_NORMALIZED_COORDS_FALSE\n"
              "                                       | CLK_ADDRESS_CLAMP_TO_EDGE\n"
              "                                       | CLK_FILTER_NEAREST;\n";
    }

    // Global GPGMP variables
    ss <<
    "#define GPGMP_DIMENSIONALITY    3\n"
    "\n"
    // Hopefully the compiler optimizes these out:
    "__constant size_t GridModelWidth       = " << m_diffusionModel.gridWidth() << ";\n"
    "__constant size_t GridModelHeight      = " << m_diffusionModel.gridHeight()<< ";\n"
    "__constant Real   PhysicalModelWidth   = " << m_diffusionModel.physicalLength() << "f;" << "\n"
    "__constant Real   PhysicalModelHeight  = " << m_diffusionModel.physicalLength() << "f;" << "\n"
    "__constant Real   PhysicalCellWidth    = " << (m_diffusionModel.physicalLength()/ (Real)m_diffusionModel.gridWidth())  << "f;\n"
    "__constant Real   PhysicalCellHeight   = " << (m_diffusionModel.physicalLength()/ (Real)m_diffusionModel.gridHeight()) << "f;\n"
    "__constant Real   StayProbability      = " << (1.-m_diffusionModel.p0()) << "f;\n"
    "__constant size_t NumReactions         = " << m_diffusionModel.numReactions() << ";\n"
    "__constant size_t NumSpecies           = " << m_diffusionModel.numSpecies() << ";\n"
    // Define num reactions and num species also as macros (and make it clear that they are indeed macros),
    // rather than as just variables, so that they can be used in if-defs and array size definitions.
    "#define GPGMP_NUM_REACTIONS              " << m_diffusionModel.numReactions() << "\n"
    "#define GPGMP_NUM_SPECIES                " << m_diffusionModel.numSpecies() << "\n"
    "\n"
    // note Visual C++ 2010 doesn't allow for accessing cl_int4 via x,y,z,w - use numerical subscripts instead
    "__constant int BoundaryMaskX1 = " << m_diffusionModel.boundaryMask().s[0] << ";\n"
    "__constant int BoundaryMaskX2 = " << m_diffusionModel.boundaryMask().s[1] << ";\n"
    "__constant int BoundaryMaskY1 = " << m_diffusionModel.boundaryMask().s[2] << ";\n"
    "__constant int BoundaryMaskY2 = " << m_diffusionModel.boundaryMask().s[3] << ";\n"
    "\n";

    // Solver type
    switch(m_diffusionModel.solver()) {
    case stochastic_homogeneous:
        ss <<"#define GPGMP_SOLVER_STOCHASTIC_HOMOGENEOUS\n";
        break;
    case stochastic_inhomogeneous:
        ss <<"#define GPGMP_SOLVER_STOCHASTIC_INHOMOGENEOUS\n";
        break;
    case stochastic_inhomogeneous_fpe:
        ss <<"#define GPGMP_SOLVER_STOCHASTIC_INHOMOGENEOUS_FPE\n";
        break;
    case deterministic_homogeneous_RK4:
        ss <<"#define GPGMP_SOLVER_DETERMINISTIC_RK4\n";
        break;
    case deterministic_homogeneous_aqss:
        ss <<"#define GPGMP_SOLVER_DETERMINISTIC_AQSS\n";
        break;
    case deterministic_inhomogeneous_RK4:
        ss <<"#define GPGMP_SOLVER_DETERMINISTIC_INHOMOGENEOUS_RK4\n";
        break;
    case deterministic_inhomogeneous_aqss:
        ss <<"#define GPGMP_SOLVER_DETERMINISTIC_INHOMOGENEOUS_AQSS\n";
        break;
    default:
        BOOST_LOG_TRIVIAL(error) <<"Can't find solver.";
        //TODO: Throw exception here!!
        break;
    }// which solver used

    ss <<"\n";

    // ballistic boundary conditions for individuals
    if (m_diffusionModel.hasBallisticBoundaryConditions())
        ss <<"#define GPGMP_BALLISTIC\n";

    // Since OpenCL doesn't support templates we need to define the state type
    ss <<"#define GPGMP_STATE_TYPE ";
    if (m_diffusionModel.solver() & stochastic_mask)
        ss <<" int\n";
    else if (m_diffusionModel.solver() & deterministic_mask)
        ss <<" Real\n";

    // we choose if we need the kernel debug
    if (m_debugKernel){
        ss <<"#define GPGMP_DEBUG_KERNEL\n";
        ss <<"#define GPGMP_NUM_ERRORS "<<GMP_NUM_ERRORS<<"\n";
    }
}

void Solver::oclReadAndReplaceTemplateFile(std::ostringstream &ss, const char *clFileName) const
{
    // TODO: Probably better to first check if tag is present and then generate the code snippet only if needed
    std::string clFullPath = (boost::filesystem::path(m_diffusionModel.oclSourcePath()) / clFileName).string();
    BOOST_LOG_TRIVIAL(debug) <<"Reading in CL Source file "<<clFullPath<<".";
    std::ifstream ifs(clFullPath.c_str());
    if (ifs.fail()) {
        BOOST_LOG_TRIVIAL(error) << "Error: failed to open " << clFileName << " for reading from. Exiting..." << std::endl;
        exit(1);
    }
    std::string program_str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();

    // Modify the openCL kernel template
    // todo: doesn't understand localization yet

    // update fields kernel
    std::ostringstream ossUpdateFields(std::ostringstream::out);
    for (size_t funcSpeciesIndex=0; funcSpeciesIndex < m_diffusionModel.numSpecies(); ++funcSpeciesIndex)
    {
        const std::string &funcSpeciesId = m_diffusionModel.species()[funcSpeciesIndex]->id();
        const size_t funcSpeciesOffset = funcSpeciesIndex * m_diffusionModel.gridArea();
        const std::map<std::string, std::vector<Real> > fields = m_diffusionModel.getFieldParameters();

        ossUpdateFields
        << "    const size_t _cellIndex = get_global_id_2d;\n"
        << "\n"
        << "    const size_t GridX = get_global_id(0);\n"
        << "    const size_t GridY = get_global_id(1);\n"
        << "#pragma unused(GridX)\n"
        << "#pragma unused(GridY)\n"
        << "\n";
        //<< "    const Real PhysicalX = getPhysicalX();\n"
        //<< "    const Real PhysicalY = getPhysicalY();\n"
        //<< "#pragma unused(PhysicalX)\n"
        //<< "#pragma unused(PhysicalY)\n"
        //<< "\n";

        // Define and load tehe species for all species
        // NOTE: given that this code gets generated for each species and only some actually store the value back, most of this should usually be optimized out
        for (size_t defineSpeciesIndex=0; defineSpeciesIndex < m_diffusionModel.numSpecies(); ++defineSpeciesIndex)
        {
            const std::string &defineSpeciesId = m_diffusionModel.species()[defineSpeciesIndex]->id();
            const size_t defineSpeciesOffset = defineSpeciesIndex * m_diffusionModel.gridArea();
            ossUpdateFields
            << "     Real " << defineSpeciesId << " = _state[get_global_id_for_species_2d("<<defineSpeciesIndex<<")];\n";
        }

        ossUpdateFields << "\n";

        // insert macros for fields
        for (std::map<std::string, std::vector<Real> >::const_iterator itt = fields.begin();
             itt != fields.end();
             ++itt) {

            const std::string id = itt->first;
            const size_t offset = m_fieldOffsets.find(id)->second;
            ossUpdateFields << "     #define " << id <<" _field["<<offset<<" + _cellIndex]\n";
        }

        // copy over the user-defined routine
        ossUpdateFields << m_diffusionModel.returnSetFieldMethod().c_str() << "\n";

        // undef macros for fields
        for (std::map<std::string, std::vector<Real> >::const_iterator itt = fields.begin();
             itt != fields.end();
             ++itt) {

            ossUpdateFields << "     #undef " << itt->first << "\n";
        }

        // Copy back the species to device
        for (size_t defineSpeciesIndex=0; defineSpeciesIndex < m_diffusionModel.numSpecies(); ++defineSpeciesIndex)
        {
            const std::string &defineSpeciesId = m_diffusionModel.species()[defineSpeciesIndex]->id();
            const size_t defineSpeciesOffset = defineSpeciesIndex * m_diffusionModel.gridArea();
            ossUpdateFields
            << "      _state[get_global_id_for_species_2d("<<defineSpeciesIndex<<")] = " << defineSpeciesId << ";\n";
        }

        ossUpdateFields << "\n";


        // and replace it
        boost::replace_first(program_str, "// <<<! updateFields !>>>", ossUpdateFields.str());
    }

    /*
     * inhomogeneous_computeDiffusionConstants_forSpecies_?
     */
    std::ostringstream ossComputeDiffusionConstants_forSpecies_ALL(std::ostringstream::out);
    for (size_t funcSpeciesIndex=0; funcSpeciesIndex < m_diffusionModel.numSpecies(); ++funcSpeciesIndex)
    {
        const std::string &funcSpeciesId = m_diffusionModel.species()[funcSpeciesIndex]->id();
        const size_t funcSpeciesOffset = funcSpeciesIndex * m_diffusionModel.gridArea();
        const std::map<std::string, std::vector<Real> > fields = m_diffusionModel.getFieldParameters();

        // compute the drift and diffusivity field
        // TODO: would be nice to have templates..
        // TODO: Maybe we can put the drift/diffusivity into image object (writable)?
        // TODO: But then we can't read it in the same kernel but we usually don't need to
        ossComputeDiffusionConstants_forSpecies_ALL
        << "__kernel void inhomogeneous_computeDiffusionConstants_forSpecies_" << funcSpeciesId << "(\n"
        << "    __global Real *_diffusionConstantsX,\n"
        << "    __global Real *_diffusionConstantsY,\n"
        << "    __global Real *_rx,\n"
        << "    __global Real *_ry,\n"
        << "    __global GPGMP_STATE_TYPE *_state,\n";

        if (m_hasFields) ossComputeDiffusionConstants_forSpecies_ALL << "    __global Real *d_field,\n";

        ossComputeDiffusionConstants_forSpecies_ALL
        << "    const Real PhysicalSimTime /* exposed to user */,\n"
        << "    __global Real *_errors,\n" // TODO: remove this?? otherwise use it!!
        << "    __global Real *_sumMoments,\n"
        << "    __global Real *_sumStates)\n"
        << "{\n"
        << "    const size_t _cellIndex = get_global_id_2d;\n"
        << "\n"
        << "    const size_t GridX = get_global_id(0);\n"
        << "    const size_t GridY = get_global_id(1);\n"
        << "#pragma unused(GridX)\n"
        << "#pragma unused(GridY)\n"
        << "\n"
        << "    const Real PhysicalX = getPhysicalX();\n"
        << "    const Real PhysicalY = getPhysicalY();\n"
        << "#pragma unused(PhysicalX)\n"
        << "#pragma unused(PhysicalY)\n"
        << "\n"
        << "// <<<! defineUserParameters !>>>\n"
        << "\n";

        // Define and load _species_t for all species
        // NOTE: given that this code gets generated for each species and only some actually store the value back, most of this should usually be optimized out
        for (size_t defineSpeciesIndex=0; defineSpeciesIndex < m_diffusionModel.numSpecies(); ++defineSpeciesIndex)
        {
            const std::string &defineSpeciesId = m_diffusionModel.species()[defineSpeciesIndex]->id();
            const size_t defineSpeciesOffset = defineSpeciesIndex * m_diffusionModel.gridArea();
            ossComputeDiffusionConstants_forSpecies_ALL
            << "     _species_t _" << defineSpeciesId << " = {\n"
            << "          _state[get_global_id_for_species_2d("<< defineSpeciesIndex <<")],\n"
            << "          _state+get_global_area_2d*"<<defineSpeciesIndex<<",\n"
            << "          getMeanX(_sumMoments, _sumStates, " << defineSpeciesIndex << "),\n"
            << "          _diffusionConstantsX[_cellIndex + " << defineSpeciesOffset << "],\n"
            << "          _diffusionConstantsY[_cellIndex + " << defineSpeciesOffset << "],\n"
            << "          _rx[_cellIndex + " << defineSpeciesOffset << "],\n"
            << "          _ry[_cellIndex + " << defineSpeciesOffset << "]\n"
            << "     };\n"
            << "     _species_t *" << defineSpeciesId << " = &_"<<defineSpeciesId<<";\n"
            << "     #pragma unused(_" << defineSpeciesId << ")\n"
            << "     #pragma unused(" << defineSpeciesId << ")\n";
        }

        // NOTE: given that this code gets generated for each species and only some actually store the value back, most of this should usually be optimized out
        ossComputeDiffusionConstants_forSpecies_ALL
        << "\n"
        << "    _species_t *All = " << funcSpeciesId <<";\n";
        // use macros, as C doesn't support refs (we want the two always to point to the same value, without subjecting the user to pointers :)
        // TODO: macros won't work - need to think of something else.
/*        << "     #define AllMeanX "        << funcSpeciesId << ".MeanX\n"
       << "     #define AllDiffusivityX " << funcSpeciesId << ".DiffusivityX\n"
        << "     #define AllDiffusivityY " << funcSpeciesId << ".DiffusivityY\n"
        << "     #define AllDriftX "       << funcSpeciesId << ".DriftX\n"
        << "     #define AllDriftY "       << funcSpeciesId << ".DriftY\n"
        << "     #define AllSpecies " << funcSpeciesId <<"\n";*/


        if (m_hasFields) {
            // insert macros for fields
            for (std::map<std::string, std::vector<Real> >::const_iterator itt = fields.begin();
                 itt != fields.end();
                 ++itt) {

                const std::string id = itt->first;
                const size_t offset = m_fieldOffsets.find(id)->second;
                ossComputeDiffusionConstants_forSpecies_ALL
                        << "     #define " << id <<" d_field["<<offset<<" + _cellIndex]\n";
            }
        }

        ossComputeDiffusionConstants_forSpecies_ALL
            << m_diffusionModel.computeDriftDiffusivityMethod();

        ossComputeDiffusionConstants_forSpecies_ALL << "\n";

        if (m_hasFields) {
            // undef macros for fields
            for (std::map<std::string, std::vector<Real> >::const_iterator itt = fields.begin();
                 itt != fields.end();
                 ++itt) {

                ossComputeDiffusionConstants_forSpecies_ALL
                        << "     #undef " << itt->first << "\n";
            }
        }
        ossComputeDiffusionConstants_forSpecies_ALL << "\n";

        // Save _species_t back
        ossComputeDiffusionConstants_forSpecies_ALL
        << "     _diffusionConstantsX[_cellIndex + " << funcSpeciesOffset << "] = " << funcSpeciesId << "->DiffusivityX;\n"
        << "     _diffusionConstantsY[_cellIndex + " << funcSpeciesOffset << "] = " << funcSpeciesId << "->DiffusivityY;\n"
        << "     _rx[_cellIndex + " << funcSpeciesOffset << "] = " << funcSpeciesId << "->DriftX;\n"
        << "     _ry[_cellIndex + " << funcSpeciesOffset << "] = " << funcSpeciesId << "->DriftY;\n"
// DEBUG       << " dx = _rx[_cellIndex + " << funcSpeciesOffset << "]; diffx =  _diffusionConstantsX[_cellIndex + " << funcSpeciesOffset << "];\n"
        << "}\n";


    }
    boost::replace_first(program_str, "// <<<! inhomogeneous_computeDiffusionConstants_forSpecies_ALL !>>>", ossComputeDiffusionConstants_forSpecies_ALL.str());


    /*
     * Copy global states to local memory (and back)
     * May be used by performReactions and computeDriftDiffusivity
     */
    // TODO: Macro substitution won't work!! need to remove!!
    std::ostringstream ossStateToLocal(std::ostringstream::out);
    for (size_t i=0; i<m_diffusionModel.species().size(); i++) {
        ossStateToLocal << "    _localState[_localId + "<<i<<" * _speciesOffsetLocal] = _globalState[_globalId + "<<i<<" * _speciesOffsetGlobal]; // species "<< m_diffusionModel.species().at(i)->id()<<"\n";
        ossStateToLocal << "#define " << m_diffusionModel.species().at(i)->id() << " _localState[_localId + "<<i<<" * _speciesOffsetLocal]\n"; // not the nicest - C++ refs would be awesome...
    }
    boost::replace_first(program_str, "// <<<! copyStateToLocal !>>>", ossStateToLocal.str());

    std::ostringstream ossLocalToState(std::ostringstream::out);
    for (size_t i=0; i<m_diffusionModel.species().size(); i++) {
        ossLocalToState << "    _globalState[_globalId + "<<i<<" * _speciesOffsetGlobal] = _localState[_localId + "<<i<<" * _speciesOffsetLocal]; // species "<< m_diffusionModel.species().at(i)->id()<<"\n";
        ossLocalToState << "#undef " << m_diffusionModel.species().at(i)->id() << "\n";
    }
    boost::replace_first(program_str, "// <<<! copyLocalToState !>>>", ossLocalToState.str());


    /*
     * Define species to local state mapping
     * May be used by performReactions and computeDriftDiffusivity, computeODERHS and computeODERHSaqss
     */
    // TODO: species get defined twice here
    std::ostringstream ossDefineSpeciesFromLocal(std::ostringstream::out);
    for (size_t i=0; i<m_diffusionModel.species().size(); i++) {
        ossDefineSpeciesFromLocal << "#define " << m_diffusionModel.species().at(i)->id() << " _localState[_localId + "<<i<<" * _speciesOffsetLocal]\n"; // not the nicest - C++ refs would be awesome...
    }
    boost::replace_all/*first*/(program_str, "// <<<! defineSpeciesFromLocal !>>>", ossDefineSpeciesFromLocal.str());

    // different definition for RK4 (and alpha-QSS) solver
    std::ostringstream ossDefineSpeciesFromLocalRK4(std::ostringstream::out);
    for (size_t i=0; i<m_diffusionModel.species().size(); i++) {
        ossDefineSpeciesFromLocalRK4 << "#define " << m_diffusionModel.species().at(i)->id() << " _localState[ "<<i<<" ]\n"; // not the nicest - C++ refs would be awesome...
    }
    boost::replace_all/*first*/(program_str, "// <<<! defineSpeciesFromLocalRK4 !>>>", ossDefineSpeciesFromLocalRK4.str());


    std::ostringstream ossUndefSpeciesFromLocal(std::ostringstream::out);
    for (size_t i=0; i<m_diffusionModel.species().size(); i++) {
        ossDefineSpeciesFromLocal << "#undef " << m_diffusionModel.species().at(i)->id() << "\n";
    }
    boost::replace_all/*first*/(program_str, "// <<<! undefSpeciesFromLocal !>>>", ossUndefSpeciesFromLocal.str());

    /*
     * Generate the parameter code
     * May be used by performReactions and computeDriftDiffusivity, computeODERHS and computeODERHSaqss
     *
     * We let the optimizer deal with them as needed- they may or may not be used by the code below.
     * It's best that we define them as variables in this scope here, rather than macros, otherwise we may run into
     * variable name conflicts.
     */
    // FUTURE: use the value from python in the future? so that it can be modified? but that would mean other changes too...
    std::ostringstream ossDefineUserParameters(std::ostringstream::out);
    typedef std::map<std::string, Real> Parameters_t;
    BOOST_FOREACH(const Parameters_t::value_type &i, m_diffusionModel.parameters()) {
        ossDefineUserParameters
        << "    const Real " << i.first << " = " << i.second << ";\n"
        << "    #pragma unused(" << i.first << ")\n"; // supress the "unused variable" warning if the user doesn't end up using it
    }
    boost::replace_all/*_first*/(program_str, "// <<<! defineUserParameters !>>>", ossDefineUserParameters.str());

    // compute reaction hazards
    std::ostringstream ossComputeReactionHazards(std::ostringstream::out);
    ossComputeReactionHazards << "      _h[0] = 0.;\n";
    // loop through reactions
    for (int i=0; i<m_diffusionModel.numReactions(); i++) {
        // get rate law
        std::string rateLaw = m_diffusionModel.reactions().at(i)->kineticLaw();
        ossComputeReactionHazards << "      // Reaction '"<<m_diffusionModel.reactions().at(i)->id()<<"'\n";
        ossComputeReactionHazards << "      _h["<<i+1<<"] = _h[0] + ";
        // if localized we need to fetch the mask
        if ((m_diffusionModel.reactions().at(i))->isLocalized())
            ossComputeReactionHazards << " read_imagef(_reactionMask, _reactionMaskSampler, (int4)(GridX, GridY, "<<i<<", 0)).x * ";
        ossComputeReactionHazards << rateLaw <<";\n";
        ossComputeReactionHazards << "      _h[0] = _h["<< (i+1) <<"];\n";
    }
    boost::replace_first(program_str, "// <<<! computeReactionHazards !>>>", ossComputeReactionHazards.str());

    // perform the reactions according to stoichiometry
    std::ostringstream ossPerformReaction(std::ostringstream::out);
    // pick reaction index according to random number
    ossPerformReaction <<"        // perform reactions\n";
    ossPerformReaction <<"        int _index;\n\n";
    for (int i=0; i<(m_diffusionModel.numReactions()); i++) {
        ossPerformReaction <<"        // Reaction '"<<m_diffusionModel.reactions().at(i)->id()<<"'\n";
        if (i==0) {
            ossPerformReaction <<"        _index = isless(_r2, _h[1]);\n";
        } else if (i==(m_diffusionModel.numReactions()-1)) {
            ossPerformReaction <<"        _index = isgreaterequal(_r2, _h["<< i <<"]) * isless(_r2, _h["<< (i+1) <<"]);\n";
        } else {
            ossPerformReaction <<"        _index = isgreaterequal(_r2, _h["<< i <<"]) * isless(_r2, _h["<< (i+1) <<"]);\n";
        }

        // read out the reactant stoichiometry for this reaction and write it
        const std::map<Species *, StoichiometryEntry> stoichr = m_diffusionModel.reactions().at(i)->reactantStoichiometryMap();


        // if we need kernel debug store reaction index and reactant state
        if (m_debugKernel) {
            std::ostringstream dosssl(std::ostringstream::out);
            std::map<Species *, StoichiometryEntry>::const_iterator dit;
            ossPerformReaction << "        if (_index > 0) {_error[get_global_id_2d * GPGMP_NUM_ERRORS + 2]="<<i
                               <<"; _error[get_global_id_2d * GPGMP_NUM_ERRORS + 3]=_index; _error[get_global_id_2d * GPGMP_NUM_ERRORS + 4]=_r2; _error[get_global_id_2d * GPGMP_NUM_ERRORS + 5]=_h["<<i
                              <<"]; _error[get_global_id_2d * GPGMP_NUM_ERRORS + 6]=_h["<<i+1<<"];";
            uint ri = 7;
            for ( dit=stoichr.begin() ; dit != stoichr.end(); dit++ ) {
                dosssl << " _error[get_global_id_2d * GPGMP_NUM_ERRORS + "<<ri++<<"] = " << ((*dit).first)->id()<<";";
            }
            ossPerformReaction << dosssl.str()<< "} \n\n";
        }

        std::ostringstream osssl(std::ostringstream::out);
        std::map<Species *, StoichiometryEntry>::const_iterator it;
        for ( it=stoichr.begin() ; it != stoichr.end(); it++ ) {
            osssl << "        " << ((*it).first)->id() << " = " << ((*it).first)->id() << " - "<< ((*it).second) << " * _index; \n";
        }
        // read out the product stoichiometry for this reaction and write it
        const std::map<Species *, StoichiometryEntry> stoichp = m_diffusionModel.reactions().at(i)->productStoichiometryMap();
        for ( it=stoichp.begin() ; it != stoichp.end(); it++ ) {
            osssl  << "        " << ((*it).first)->id() << " = " << ((*it).first)->id() << " + "<< ((*it).second) << " * _index; \n";
        }
        ossPerformReaction << osssl.str() <<"\n";
    }
    boost::replace_first(program_str, "// <<<! performReactions !>>>", ossPerformReaction.str());

    // replace the method do compute diffusivity and drift
    boost::replace_first(program_str, "// <<<! computeDriftDiffusivity !>>>", m_diffusionModel.computeDriftDiffusivityMethod());

    // replace method to store diffusion constants in shared memory (for deterministic)
    if(boost::find_first(program_str, "// <<<! deterministicStoreDiffusionConstants !>>>")) {
        std::ostringstream ossStoreDiffusionConstants(std::ostringstream::out);

        // go through all species and store diffusion constants
        for (size_t i=0; i<m_diffusionModel.species().size(); i++) {
            ossStoreDiffusionConstants <<"    _physicalDiffusionConstants["<<i<<"] = "<< m_diffusionModel.species().at(i)->diffusionConstant()
                                       <<"; // species "<< m_diffusionModel.species().at(i)->id()<<"\n";
        }
        boost::replace_first(program_str, "// <<<! deterministicStoreDiffusionConstants !>>>", ossStoreDiffusionConstants.str());
    }

    // replace method to compute the RHS of the deterministic ODE
    if(boost::find_first(program_str, "// <<<! computeODERHS !>>>")) {
        std::ostringstream ossComputeODERHS(std::ostringstream::out);

        // go through all species
        for (size_t i=0; i<m_diffusionModel.species().size(); i++) {
            // start with zero term in case this species is not reacting
            ossComputeODERHS << "    _deriv[" <<i<<"] = 0.";

            // go through all reactions
            for (int j=0; j<m_diffusionModel.numReactions(); j++) {
                // get total stoichiometry for reaction
                std::map<Species *, StoichiometryEntry> stoichr = m_diffusionModel.reactions().at(j)->reactantStoichiometryMap();
                std::map<Species *, StoichiometryEntry> stoichp = m_diffusionModel.reactions().at(j)->productStoichiometryMap();
                int totalStoich = stoichp[m_diffusionModel.species().at(i)] - stoichr[m_diffusionModel.species().at(i)];

                // add term only if it affects our species
                if (totalStoich!=0) {
                    if (totalStoich > 0) {
                        ossComputeODERHS << " + " << abs(totalStoich)<<". * ";
                    } else {
                        ossComputeODERHS << " - " << abs(totalStoich) <<". * ";
                    }

                    // get rate law
                    if (m_diffusionModel.reactions().at(j)->isLocalized()) {
                            ossComputeODERHS << " read_imagef(_reactionMask, _reactionMaskSampler, (int4)(GridX, GridY, "<<j<<", 0)).x * ";
                    }

                    // get rate law
                    std::string rateLaw = m_diffusionModel.reactions().at(j)->deterministicLaw();
                    ossComputeODERHS << rateLaw;
                }
            }
            ossComputeODERHS <<";\n";
        }
        boost::replace_first(program_str, "// <<<! computeODERHS !>>>", ossComputeODERHS.str());
    }// <<<! computeODERHS !>>>

    // replace method to compute the RHS of the deterministic ODE for alpha-QSS
    if(boost::find_first(program_str, "// <<<! computeODERHSaqss !>>>")) {
        std::ostringstream ossODERHSaqss(std::ostringstream::out);
        // go through all species
        for (size_t i=0; i<m_diffusionModel.species().size(); i++) {
            ossODERHSaqss <<"    // Species '"<<m_diffusionModel.species().at(i)->id()<<"'\n";

            // Need to generate production/destruction terms separately
            std::ostringstream ossProduction(std::ostringstream::out);
            std::ostringstream ossDestruction(std::ostringstream::out);
            // start with zero (in case there's no production/destruction term for this species)
            ossProduction << "    _q[" <<i<<"] = 0.";
            ossDestruction<< "    _d[" <<i<<"] = 0.";

            // go through all reactions
            for (int j=0; j<m_diffusionModel.numReactions(); j++) {
                // get total stoichiometry for reaction
                std::map<Species *, StoichiometryEntry> stoichr = m_diffusionModel.reactions().at(j)->reactantStoichiometryMap();
                std::map<Species *, StoichiometryEntry> stoichp = m_diffusionModel.reactions().at(j)->productStoichiometryMap();
                int totalStoich = stoichp[m_diffusionModel.species().at(i)] - stoichr[m_diffusionModel.species().at(i)];

                // add term only if it affects our species
                if (totalStoich!=0) {
                    // get rate law
                    std::string rateLaw = m_diffusionModel.reactions().at(j)->deterministicLaw();
                    // add to production or destruction array
                    if (totalStoich > 0) {
                        ossProduction << " + " << abs(totalStoich)<<". * ";
                        // localization
                        if (m_diffusionModel.reactions().at(j)->isLocalized()) {
                                ossProduction << " read_imagef(_reactionMask, _reactionMaskSampler, (int4)(GridX, GridY, "<<j<<", 0)).x * ";
                        }
                        // and add rate law
                        ossProduction << rateLaw;
                    } else {
                        ossDestruction<< " + " << abs(totalStoich) <<". * ";
                        // localization
                        if (m_diffusionModel.reactions().at(j)->isLocalized()) {
                                ossDestruction<< " read_imagef(_reactionMask, _reactionMaskSampler, (int4)(GridX, GridY, "<<j<<", 0)).x * ";
                        }
                        // add rate law
                        ossDestruction << rateLaw;
                    }
                } // rate law term (if stoich>0)
            }// reaction
            ossODERHSaqss << ossProduction.str() <<";\n"<< ossDestruction.str() <<";\n\n";
        }// loop through species
        boost::replace_first(program_str, "// <<<! computeODERHSaqss !>>>", ossODERHSaqss.str());
    }//<<<! computeODERHSaqss !>>>


    // and add it to the definitions
    ss << program_str;

    // write the code to file (for debugging purpose)
    // todo: make this behaviour optional
    std::ofstream ofs;
    std::string finalname="final_code_";
    finalname += clFileName;
    ofs.open(finalname.c_str());
    ss.flush();
    ofs << ss.str();
    ofs.close();

    /*
    // for debug purpose .. read in code from file
    std::string line;
    std::ifstream ifs("final_code.cl");
    if (ifs.is_open()) {
      // throw first line
      getline(ifs, line);
      while (ifs.good()) {
    getline(ifs, line);
    std::cout << line << std::endl;
    ss << line << std::endl;
      }
      ifs.close();
    }
  */
    ss.flush();
}//oclReadAndReplaceTemplateFile

cl::Program Solver::oclBuildProgram(std::ostringstream &ss) const
{
    // get the include path right
    std::string cl_options;
#if !defined(__APPLE__) && !defined(__MACOSX) // TODO: replace this with a runtime check for an nvidia platform
    cl_options += " -cl-nv-verbose";
#endif
    cl_options += " -I";
    cl_options += (boost::filesystem::path(m_diffusionModel.oclSourcePath()) / "cl_include").string();

    BOOST_LOG_TRIVIAL(debug) <<"Building program with options: "<<cl_options.c_str();

    // Build the program
    cl::Program program(*m_context, ss.str());
    try {
        program.build(cl_options.c_str());
        //BOOST_LOG_TRIVIAL(debug)
        //        << "OpenCL program build successful. Build log:" << std::endl
        //        << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(c_contextDevice);
    } catch (cl::Error err) {
        BOOST_LOG_TRIVIAL(error)
                << "OpenCL program build failed:" << std::endl
                << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(c_contextDevice)
                << "Exiting.";
        exit(1);
    }
    return program;
}

} // namespace gpgmp
