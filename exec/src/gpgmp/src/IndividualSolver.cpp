#include <boost/foreach.hpp>
#include <boost/filesystem.hpp>

#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <limits>

#include "Solver.h"
#include "IndividualSolver.h"
#include "Species.h"

// these are just for debugging..
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/log/trivial.hpp>

namespace gpgmp {

IndividualSolver::IndividualSolver(DiffusionModel &model, const Solver &solver,
                                   cl::Context *context, cl::CommandQueue *queue,
                                   cl::Device &contextDevice)
    : m_diffusionModel(model), m_solver(solver), m_context(context), m_commandQueue(queue), m_contextDevice(contextDevice)
{
    // go through species and check if there is any individual based properties
    int speciesIndex = -1;

    // get species list from model
    const boost::container::stable_vector<Species *> species = model.species();

    BOOST_FOREACH(Species *s, species)
    {
        speciesIndex++;

        if (s->hasIndividualProperties()) {

            // get handle on properties and positions list
            std::map<std::string, std::vector<Real> > &props = s->getIndividualProperties();
            std::vector<int> ipos = s->getIndividualPositions();

            // we assume that the size of the initial position vector gives us the
            // total initial number of individuals of this species
            m_nIndividuals[speciesIndex] = ipos.size();
            std::cout <<"Found "<<m_nIndividuals[speciesIndex]<<" individuals for species "
                     << speciesIndex << " ("<<s->id()<<").\n";

            // how many properties?
            uint nprops = props.size();

            // create combined array for properties
            std::vector<Real> h_iprops(nprops*c_maxNumIndividuals);
            std::vector<cl_ulong> h_ipos(c_maxNumIndividuals);

            // and fill them in
            for (uint i=0; i<m_nIndividuals[speciesIndex]; i++) {
                std::map<std::string, std::vector<Real> >::iterator iter;

                int k=4; // First four properties are system (diffusivity & drift)
                for (iter = props.begin(); iter != props.end(); iter++) {
                    // check if it's one of the system properties
                    if ((iter->first) == "diffx") {
                        h_iprops[0*c_maxNumIndividuals+i] = (iter->second)[i];
                    } else if ((iter->first) == "diffy") {
                        h_iprops[1*c_maxNumIndividuals+i] = (iter->second)[i];
                    } else if ((iter->first) == "rx") {
                        h_iprops[2*c_maxNumIndividuals+i] = (iter->second)[i];
                    } else if ((iter->first) == "ry") {
                        h_iprops[3*c_maxNumIndividuals+i] = (iter->second)[i];
                    } else {
                        // additional properties will be added starting from count 5
                        h_iprops[k*c_maxNumIndividuals+i] = (iter->second)[i];
                        k++;
                    }
                } // fill in properties

                // now combine key and index
                h_ipos[i] = KEY_INDEX_PAIR((cl_uint) ipos[i], (cl_uint) i);
            }

            // now we fill in the rest with default values
            // we will need to get the particle id right though!
            for (uint i=m_nIndividuals[speciesIndex]; i<c_maxNumIndividuals; i++) {
                // properties default
                // default diffusivity is 1
                // todo: change me ..
                for (uint j=0; j<2; j++)
                    h_iprops[j*c_maxNumIndividuals + i] = 1.;

                for (uint j=2; j<nprops; j++)
                    h_iprops[j*c_maxNumIndividuals + i] = 0.;

                // position default is UINT_MAX
                h_ipos[i] = KEY_INDEX_PAIR(UINT_MAX, (cl_uint) i);
            }

            try {
                // and create a buffer from them
                individualProperties[speciesIndex] = oclCreateBufferFrom(h_iprops);
                individualPositions[speciesIndex] = oclCreateBufferFrom(h_ipos);


                // and copy it over to the device
                oclCopy(h_iprops, *individualProperties[speciesIndex]);
                oclCopy(h_ipos, *individualPositions[speciesIndex]);

                // create cell indices buffers
                d_cellIndicesLeft[speciesIndex]  = oclCreateBuffer<cl_uint>(model.gridArea());
                d_cellIndicesRight[speciesIndex] = oclCreateBuffer<cl_uint>(model.gridArea());

                // fill buffers with default values
                // todo: try to work out how to use clEnqueueFillBuffers for this..
                std::vector<cl_uint> def(model.gridArea(), 0);
                oclCopy(def, *d_cellIndicesLeft[speciesIndex]);
                oclCopy(def, *d_cellIndicesRight[speciesIndex]);

            } catch (cl::Error err) {
                std::cerr <<"ERROR in IndividualSolver while creating buffers for species "<<s->id()<<": "
                         << err.what()<<" ("
                         <<oclErrorString(err.err())<<").\n"<<std::flush;;
                throw err;
            }

           /*
            // now print and exit
            std::cout <<"Species "<<s->id()<<" has individual-based properties.\n";
            std::cout <<"size of long is "<<sizeof(cl_ulong) <<".\n";
            for (int i=0; i<c_maxNumIndividuals; i++) {
                std::cout << h_ipos[i] << " " << GET_INDEX(h_ipos[i]) << " " << GET_KEY(h_ipos[i]) << " " << ipos[i]
                          <<" "<<h_iprops[0*c_maxNumIndividuals + i] <<" "<<h_iprops[1*c_maxNumIndividuals + i]<<"\n";
            }*/

        } // create list if species has individual properties
    } // loop over all species

    // Compute topology for radix sort and fill it in
    // we need to set the WG size and cta size to allow small key numbers
    // 4*wgSizeScan <= arrayLength = clDevice.nkeys/2/clDevice.cta_size*16 <= 4*wgSizeScan*wgSizeScan
    // wg size must be multiple of two
    //m_wgSizeScan = c_maxNumIndividuals/c_ctaSize;
    m_wgSizeScan = c_ctaSize;
    std::cout <<"Setting workgroup size for scan of "<<m_wgSizeScan<<".\n"<<std::flush;

    // Compute number of blocks for radix sort
    m_radixNumBlocks = ((c_maxNumIndividuals % (c_ctaSize * 4)) == 0) ?
                (c_maxNumIndividuals / (c_ctaSize * 4)) : (c_maxNumIndividuals / (c_ctaSize* 4) + 1);
    m_radixNumBlocks2 = ((c_maxNumIndividuals % (c_ctaSize * 2)) == 0) ?
      (c_maxNumIndividuals / (c_ctaSize * 2)) : (c_maxNumIndividuals / (c_ctaSize * 2) + 1);

    // create remaining device buffers
    d_tempKeys = oclCreateBuffer<cl_ulong>(c_maxNumIndividuals);
    d_mCounters = oclCreateBuffer<cl_uint>(c_warp_size * m_radixNumBlocks);
    d_mCountersSum = oclCreateBuffer<cl_uint>(c_warp_size * m_radixNumBlocks);
    d_mBlockOffsets = oclCreateBuffer<cl_uint>(c_warp_size * m_radixNumBlocks);
    d_buffer = oclCreateBuffer<cl_uint>(c_maxNumIndividuals/1024); // TODO: where does the number come from??
    d_reductionOffsets = oclCreateBuffer<cl_int>(c_maxNumIndividuals);
    d_reductionOffsetsBlocks = oclCreateBuffer<cl_int>(c_maxNumIndividuals/2/16);
    d_reductionOffsetsAdding = oclCreateBuffer<cl_int>(c_maxNumIndividuals);
    d_reductionOffsetsBlocksAdding = oclCreateBuffer<cl_int>(c_maxNumIndividuals/2/16);
    d_total = oclCreateBuffer<cl_int>(1);
    d_totalAdded = oclCreateBuffer<cl_int>(1);
    d_lostParticles = oclCreateBuffer<cl_ulong>(c_maxNumIndividuals);
    d_particlesBuffer = oclCreateBuffer<cl_ulong>(c_maxNumIndividuals);
    d_newParticles = oclCreateBuffer<cl_uint>(m_diffusionModel.gridArea());

    // device buffers for segmented scan
    uint nscan = m_diffusionModel.gridArea()/2;
    uint nBlockSize = c_ssWgSize*c_ssWgSize;
    uint nBlocks = nscan/nBlockSize;
    std::cout <<"Allocating "<<nBlocks<<" blocks for segmented scan.\n";

    d_ssFlags = oclCreateBuffer<cl_uchar>(m_diffusionModel.gridArea());
    d_domainReductionOffsets = oclCreateBuffer<cl_uint>(m_diffusionModel.gridArea());
    d_ssDataBlocks = oclCreateBuffer<cl_uint>(nBlocks);
    d_ssFlagsBlocks = oclCreateBuffer<cl_uchar>(nBlocks);
    d_ssFirstBlockFlag = oclCreateBuffer<cl_uchar>(nBlocks);
    d_ssPartialFlags = oclCreateBuffer<cl_uchar>(nscan*2);

    // now build the kernels
    buildKernels();

    //writeProperties("properties_initial.dat", 0);
} // constructor

void IndividualSolver::buildKernels() {
    // Read in source code
    std::string clFileName = "IndividualSolver.cl";

    std::string clFullPath = (boost::filesystem::path(m_diffusionModel.oclSourcePath()) / clFileName).string();
    std::cout <<"Reading in CL Source file "<<clFullPath<<". \n";
    std::ifstream ifs(clFullPath.c_str());
    if (ifs.fail()) {
        std::cerr << "Error: failed to open " << clFileName << " for reading from. Exiting..." << std::endl;
        exit(1);
    }
    std::string template_str((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();

    // add header
    std::ostringstream oss;
    m_solver.oclGenerateHeader(oss);

    std::string program_str = oss.str() + template_str;

    // replace tags
    // TODO: we should integrate the individual solver better into the whole solver framework
    // compute diffusivity/drift method for individual solver
    std::ostringstream ossComputeDiffusivityDrift(std::ostringstream::out);
    ossComputeDiffusivityDrift
            << "inline void computeDiffusivityDrift(_species_t *actualSpecies, int speciesIndex) {\n"
            << "  // define the species structs\n"
            << "  _species_t *_allSpecies[" << m_diffusionModel.numSpecies() <<"]; // contains list of all species\n\n";

    ossComputeDiffusivityDrift
            << "  //define dummies for all species\n";
    for (size_t funcSpeciesIndex=0; funcSpeciesIndex < m_diffusionModel.numSpecies(); ++funcSpeciesIndex)
    {
        const std::string &funcSpeciesId = m_diffusionModel.species()[funcSpeciesIndex]->id();
        ossComputeDiffusivityDrift
                << "  _species_t _" << funcSpeciesId <<" = {0,0,0,0,0,0,0}; // will be filled in later for each individual\n"
                << "  _allSpecies[" << funcSpeciesIndex <<"] = &_"<<funcSpeciesId<<";\n";
    }
    ossComputeDiffusivityDrift
            << "\n"
            << "  // fill in values for this particle\n"
            << "  _allSpecies[speciesIndex] = actualSpecies;\n"
            <<"\n"
           << "  // and define the aliases for the species\n"
           << "  _species_t *All = _allSpecies[speciesIndex];\n";

    for (size_t funcSpeciesIndex=0; funcSpeciesIndex < m_diffusionModel.numSpecies(); ++funcSpeciesIndex)
    {
        const std::string &funcSpeciesId = m_diffusionModel.species()[funcSpeciesIndex]->id();
        ossComputeDiffusivityDrift
                << "  _species_t *" << funcSpeciesId <<" = _allSpecies["<<funcSpeciesIndex<<"];\n";
    }

    // define the parameters
    typedef std::map<std::string, Real> Parameters_t;
    BOOST_FOREACH(const Parameters_t::value_type &i, m_diffusionModel.parameters()) {
        ossComputeDiffusivityDrift
        << "  const Real " << i.first << " = " << i.second << ";\n"
        << "  #pragma unused(" << i.first << ")\n"; // supress the "unused variable" warning if the user doesn't end up using it
    }

    // we can't use macros -- Expressions like B->DiffusivityX will evaluate to B->_allSpecies[speciesIndex] etc..
    // instead we define All->MeanX etc..
    ossComputeDiffusivityDrift
            << "\n"
            << "  // compute diffusivity/drift from user-defined method\n"
            << m_diffusionModel.computeDriftDiffusivityMethod()
            << "}\n";

    boost::replace_first(program_str, "// <<<! computeDiffusivityDrift !>>>", ossComputeDiffusivityDrift.str());


    // method to initialize properties of new particles
    std::ostringstream ossNewIndividuals(std::ostringstream::out);
    ossNewIndividuals
            << "  // define the species structs\n"
            << "  _species_t *_allSpecies[" << m_diffusionModel.numSpecies() <<"]; // contains list of all species\n\n";

    ossNewIndividuals
            << "  //define dummies for all species\n";
    for (size_t funcSpeciesIndex=0; funcSpeciesIndex < m_diffusionModel.numSpecies(); ++funcSpeciesIndex)
    {
        const std::string &funcSpeciesId = m_diffusionModel.species()[funcSpeciesIndex]->id();
        ossNewIndividuals
                << "  _species_t _" << funcSpeciesId <<" = {0,0,0,0,0,0,0}; // will be filled in later for each individual\n"
                << "  _allSpecies[" << funcSpeciesIndex <<"] = &_"<<funcSpeciesId<<";\n";
    }
    ossNewIndividuals
            << "\n"
            << "  // fill in values for this particle\n"
            << "  _allSpecies[speciesIndex] = &_actualSpecies;\n"
            <<"\n"
           << "  // and define the aliases for the species\n"
           << "  _species_t *All = _allSpecies[speciesIndex];\n";

    for (size_t funcSpeciesIndex=0; funcSpeciesIndex < m_diffusionModel.numSpecies(); ++funcSpeciesIndex)
    {
        const std::string &funcSpeciesId = m_diffusionModel.species()[funcSpeciesIndex]->id();
        ossNewIndividuals
                << "  _species_t *" << funcSpeciesId <<" = _allSpecies["<<funcSpeciesIndex<<"];\n";
    }

    // define the parameters
    typedef std::map<std::string, Real> Parameters_t;
    BOOST_FOREACH(const Parameters_t::value_type &i, m_diffusionModel.parameters()) {
        ossNewIndividuals
        << "  const Real " << i.first << " = " << i.second << ";\n"
        << "  #pragma unused(" << i.first << ")\n"; // supress the "unused variable" warning if the user doesn't end up using it
    }

    // we can't use macros -- Expressions like B->DiffusivityX will evaluate to B->_allSpecies[speciesIndex] etc..
    // instead we define All->MeanX etc..
    ossNewIndividuals
            << "\n"
            << "  // compute diffusivity/drift from user-defined method\n"
            << m_diffusionModel.getNewIndividualsMethod()
            << "\n";
    boost::replace_first(program_str, "// <<<! newIndividualProperties !>>>", ossNewIndividuals.str());

    // write the program
    std::ofstream myfile;
    myfile.open ("final_code_individualSolver.cl");
    myfile << program_str.c_str();
    myfile.close();

    // Build the program from stream
    std::string cl_options;
    cl_options += " -cl-nv-verbose";
    cl_options += " -I";
    cl_options += (boost::filesystem::path(m_diffusionModel.oclSourcePath()) / "cl_include").string();
    cl_options += " -D WORKGROUP_SIZE_SCAN="+ static_cast<std::ostringstream*>( &(std::ostringstream() << m_wgSizeScan) )->str();

    std::cout <<"Building program with options: "<<cl_options.c_str()<<" .\n";

    // Build the program
    cl::Program program(*m_context, program_str);

    try {
        program.build(cl_options.c_str());
        std::cout
                << "OpenCL program build successful. Build log:" << std::endl
                << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(m_contextDevice) << std::endl;
    } catch (cl::Error err) {
        std::cerr
                << "OpenCL program build failed:" << std::endl
                << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(m_contextDevice) << std::endl
                << "Exiting..." << std::endl;
        exit(1);
    }

    // and create the kernels
    try {
        m_kernel_radixSortBlocksKeysOnly = Kernel(program, "radixSortBlocksKeysOnly");
        m_kernel_findRadixOffsets = Kernel(program, "findRadixOffsets");
        m_kernel_scanNaive = Kernel(program, "scanNaive");
        m_kernel_reorderDataKeysOnly = Kernel(program, "reorderDataKeysOnly");
        m_kernel_scanExclusiveLocal1 = Kernel(program, "scanExclusiveLocal1");
        m_kernel_scanExclusiveLocal2 = Kernel(program, "scanExclusiveLocal2");
        m_kernel_uniformUpdate = Kernel(program, "uniformUpdate");
        m_kernel_computeCellIndices = Kernel(program, "computeCellIndices");
        m_kernel_computeDiffusionConstants= Kernel(program, "computeDiffusionConstants");
        m_kernel_individualDiffuseX = Kernel(program, "individual_diffuseX");
        m_kernel_individualDiffuseY = Kernel(program, "individual_diffuseY");
        m_kernel_padKeys = Kernel(program, "padKeys");

        if (m_diffusionModel.numReactions() > 0) {
            m_kernel_maintainPopulation = Kernel(program, "maintainPopulation");
            m_kernel_scanIndividuals = Kernel(program, "scanIndividuals");
            m_kernel_scanIndividualsBlocks = Kernel(program, "scanIndividualsBlocks");
            m_kernel_scatterIndividuals = Kernel(program, "scatterIndividuals");
            m_kernel_writeLostParticles = Kernel(program, "writeLostParticles");
            m_kernel_computeReductionOffsets = Kernel(program, "computeReductionOffsets");
            m_kernel_clearReductionOffsetsAdding = Kernel(program, "clearReductionOffsetsAdding");
            m_kernel_updateNewParticles = Kernel(program, "updateNewParticles");
            m_kernel_copyTempToKeys = Kernel(program, "copyTempToKeys");
        }

        // for segmented scan
        m_kernel_segmentedScanBlock = Kernel(program, "segmentedScanBlock");
        m_kernel_segmentedScanUpSweep = Kernel(program, "segmentedScanUpSweep");
        m_kernel_segmentedScanDownSweep = Kernel(program, "segmentedScanDownSweep");

        // build single-scan test kernel
        m_kernel_segmentedScanBlock_single = Kernel(program, "segmentedScanBlock_single");

    } catch (cl::Error err) {
        std::cerr <<"ERROR while building kernels for IndividualSolver: "<< err.what()<<" ("
                 <<oclErrorString(err.err())<<").\n";
        exit(1);
    }
} // build kernels

// destructor
IndividualSolver::~IndividualSolver() {
    // free up all device memory
    delete d_tempKeys;
    delete d_mCounters;
    delete d_mCountersSum;
    delete d_mBlockOffsets;
    delete d_buffer;
    delete d_reductionOffsets;
    delete d_reductionOffsetsBlocks;
    delete d_reductionOffsetsAdding;
    delete d_reductionOffsetsBlocksAdding;
    delete d_total;
    delete d_totalAdded;
    delete d_lostParticles;
    delete d_particlesBuffer;
    delete d_domainReductionOffsets;
    delete d_newParticles;

    delete d_ssDataBlocks;
    delete d_ssFlagsBlocks;
    delete d_ssFirstBlockFlag;
    delete d_ssPartialFlags;
    delete d_ssFlags;

    std::map<int, cl::Buffer *>::iterator iter;
    for (iter = individualProperties.begin(); iter != individualProperties.end(); iter++)
        delete (iter->second);
    for (iter = individualPositions.begin(); iter != individualPositions.end(); iter++)
        delete (iter->second);
    for (iter = d_cellIndicesLeft.begin(); iter != d_cellIndicesLeft.end(); iter++)
        delete (iter->second);
    for (iter = d_cellIndicesRight.begin(); iter != d_cellIndicesRight.end(); iter++)
        delete (iter->second);
}// destructor

void IndividualSolver::diffuseX(int seed, int speciesIndex, float dt, cl::Buffer *d_state, cl::Buffer *d_errors) {
    // call diffusion kernel over whole grid
    uint width = m_diffusionModel.gridWidth();
    uint height = m_diffusionModel.gridHeight();

    m_kernel_individualDiffuseX.setEnqueueArgs(cl::EnqueueArgs(*m_commandQueue,
                                                     cl::NDRange(width, height),
                                                     cl::NDRange(16,16)));

    uint maxNumIndividuals = c_maxNumIndividuals;

    m_kernel_individualDiffuseX(*(individualPositions[speciesIndex]),
                                *(d_cellIndicesLeft[speciesIndex]),
                                *(d_cellIndicesRight[speciesIndex]),
                                *(individualProperties[speciesIndex]),
                                *d_state,
                                 maxNumIndividuals, seed, speciesIndex,
                                dt, *d_errors);
    //std::cout <<"Diffusing X with seed "<<seed<<"\n";
} // diffuseX

void IndividualSolver::diffuseY(int seed, int speciesIndex, float dt, cl::Buffer *d_state, cl::Buffer *d_errors) {
    // call diffusion kernel over whole grid
    uint width = m_diffusionModel.gridWidth();
    uint height = m_diffusionModel.gridHeight();

    m_kernel_individualDiffuseY.setEnqueueArgs(cl::EnqueueArgs(*m_commandQueue,
                                                     cl::NDRange(width, height),
                                                     cl::NDRange(16,16)));

    uint maxNumIndividuals = c_maxNumIndividuals;

    m_kernel_individualDiffuseY(*(individualPositions[speciesIndex]),
                                *(d_cellIndicesLeft[speciesIndex]),
                                *(d_cellIndicesRight[speciesIndex]),
                                *(individualProperties[speciesIndex]),
                                *d_state,
                                maxNumIndividuals, seed, speciesIndex,
                                dt, *d_errors);
    //std::cout <<"Diffusing Y with seed "<<seed<<"\n";
}// diffuseY

void IndividualSolver::maintainPopulation(int speciesIndex, int seed, cl::Buffer *d_state)
{
    uint width = m_diffusionModel.gridWidth();
    uint height = m_diffusionModel.gridHeight();
    uint nIndividuals = m_nIndividuals[speciesIndex];

    // todo: debug code .. remove me
    // copy positions back to device
    std::vector<cl_ulong> h_iposOld(c_maxNumIndividuals);
    oclCopy(*(individualPositions[speciesIndex]), h_iposOld);
    // end remove me

    std::cout <<"Maintain population for species " << speciesIndex << " called. \n"<<std::flush;

    // fill in temp list
    oclCopyDeviceBuffers(*(individualPositions[speciesIndex]), *d_tempKeys, c_maxNumIndividuals*sizeof(type_key));

    writeCellIndices(0, "cellIndicesBefore");

    //writeState("stateBefore.dat", d_state);
    //writePositions(0, "positionsBefore.dat");

    // remove particles which were killed in the reaction
    m_kernel_maintainPopulation.setEnqueueArgs(cl::EnqueueArgs(*m_commandQueue,
                                                               cl::NDRange(width, height),
                                                               cl::NDRange(16,16)));
    m_kernel_maintainPopulation(*(individualPositions[speciesIndex]),
                                *(d_cellIndicesLeft[speciesIndex]),
                                *(d_cellIndicesRight[speciesIndex]),
                                *d_reductionOffsets,
                                speciesIndex, *d_state, seed,
                                *d_newParticles, *d_domainReductionOffsets,
                                *d_ssFlags);
    seed++;

    //writePositions(0, "positionsBeforeReduction.dat", 0);
    //writeCellIndices(0, "cellIndicesBeforeReduction.dat");
    //writeSegmentedScanBuffers("ssScanBuffersBefore.dat");

    // do segmented scan to find out how many new particles need to be inserted
    segmentedScan(d_domainReductionOffsets, d_ssFlags, m_diffusionModel.gridWidth(), m_diffusionModel.gridHeight());
    //writeSegmentedScanBuffers("ssScanBuffersAfterScan.dat");

    // clear and compute the combined reduction offset
    std::cout <<"computing reduction offsets ..\n"<<std::flush;

    if (m_nIndividuals[speciesIndex]>0) {
        m_kernel_clearReductionOffsetsAdding.setEnqueueArgs(cl::EnqueueArgs(*m_commandQueue,
                                                                            cl::NDRange(m_nIndividuals[speciesIndex])));
        m_kernel_clearReductionOffsetsAdding(*d_reductionOffsetsAdding);
    }

    m_kernel_computeReductionOffsets.setEnqueueArgs(cl::EnqueueArgs(*m_commandQueue,
                                                                    cl::NDRange(width, height),
                                                                    cl::NDRange(16,16)));
    m_kernel_computeReductionOffsets(*d_newParticles, *d_domainReductionOffsets, *d_ssFlags,
                                     *(d_cellIndicesLeft[speciesIndex]), *d_reductionOffsetsAdding);
    //writeReductionOffsets("reductionOffsetsBefore.dat");

    // reduce individuals stream
    // we only need to do that if there's at least one active individual
   if (nIndividuals > 0) {
        const size_t wgSize = 64;
        uint nScan = 4*wgSize;
        while (nIndividuals > nScan) {
            nScan <<=1;
        }
        nScan/=2;
        std::cout <<"Scanning with "<<nScan<<" threads.\n"<<std::flush;

        m_kernel_scanIndividuals.setEnqueueArgs(cl::EnqueueArgs(*m_commandQueue,
                                                                cl::NDRange(nScan),
                                                                cl::NDRange(wgSize)));
        m_kernel_scanIndividuals.instance().setArg(4, 2*wgSize*sizeof(cl_uint),0);
        m_kernel_scanIndividuals.instance().setArg(5, 2*wgSize*sizeof(cl_uint),0);
        m_kernel_scanIndividuals(*d_reductionOffsets, *d_reductionOffsetsBlocks, *d_reductionOffsetsAdding, *d_reductionOffsetsBlocksAdding);
        uint nBlockScan = nScan/wgSize/2;

        //writeReductionOffsets("offsetsBeforeBlockReduction.dat", seed);

        std::cout <<"Reducing blocks with "<<nBlockScan<<" threads.\n"<<std::flush;

        m_kernel_scanIndividualsBlocks.setEnqueueArgs(cl::EnqueueArgs(*m_commandQueue,
                                                                      cl::NDRange(nBlockScan),
                                                                      cl::NDRange(nBlockScan)));
        m_kernel_scanIndividualsBlocks.instance().setArg(4, 2*nBlockScan*sizeof(cl_uint),0);
        m_kernel_scanIndividualsBlocks.instance().setArg(5, 2*nBlockScan*sizeof(cl_uint),0);
        m_kernel_scanIndividualsBlocks(*d_reductionOffsets, *d_reductionOffsetsBlocks,
                                       *d_reductionOffsetsAdding, *d_reductionOffsetsBlocksAdding);

        //writeReductionOffsets("reductionOffsetsAfter.dat");
        //writeReductionOffsets("offsetsAfterBlockReduction.dat", seed);

        // scatter individuals
        // same topology as the scan kernel
        m_kernel_scatterIndividuals.setEnqueueArgs(cl::EnqueueArgs(*m_commandQueue,
                                                                   cl::NDRange(nScan),
                                                                   cl::NDRange(wgSize)));
        /*
    m_kernel_scatterIndividuals(*(individualPositions[speciesIndex]),
                                *d_reductionOffsets, *d_reductionOffsetsBlocks,
                                *d_total, *d_lostParticles,
                                m_nIndividuals[speciesIndex]);
                                */
        // todo: check where the keys come from..
        m_kernel_scatterIndividuals(*d_tempKeys, *(individualPositions[speciesIndex]),
                                    *d_reductionOffsets, *d_reductionOffsetsBlocks,
                                    *d_reductionOffsetsAdding, *d_reductionOffsetsBlocksAdding,
                                    *d_totalAdded, *d_total, *d_lostParticles,
                                    m_nIndividuals[speciesIndex],
                                    *d_domainReductionOffsets, *d_newParticles);
    } // if at least one individual is active

    // copy back survivor count ..
    std::vector<cl_uint> h_total(1);
    std::vector<cl_uint> h_totalAdded(1);

    // if we had to reduce, particle counts are in d_total, d_totalAdded
    if (nIndividuals>0) {
        oclCopy(*d_total, h_total);
        oclCopy(*d_totalAdded, h_totalAdded);
    } else {
        // only added particles .. how many is in the reduction offsets
        oclCopy(*d_reductionOffsetsAdding, &(h_totalAdded.front()), 1);
        h_total[0]=0;
    }
        const uint lost = m_nIndividuals[speciesIndex]-h_total[0];
    std::cout <<"After reduction we have "<<h_total[0]<<" particles left. We added " << h_totalAdded[0]<<" and lost " << lost<< ". \n"<<std::flush;
    uint eolNew = h_total[0] + h_totalAdded[0];

    if (eolNew > c_maxNumIndividuals) {
        std::cerr <<"ERROR: Required number of individuals (" << eolNew <<") is greater than maximum"
                    " numer of individuals (" << c_maxNumIndividuals<<"). Exiting..\n";
        exit(1);
    }
    //writePositions(0, "positionsAfterScatter.dat");

    if (h_totalAdded[0]>0) {
        // copy lost particles from end of old list
        std::cout <<"copying new particles ..\n"<<std::flush;

        const uint propertiesOffset = c_maxNumIndividuals;
        m_kernel_updateNewParticles.setEnqueueArgs(cl::EnqueueArgs(*m_commandQueue,
                                                   cl::NDRange(m_diffusionModel.gridWidth(), m_diffusionModel.gridHeight()),
                                                   cl::NDRange(16,16)));
        m_kernel_updateNewParticles(*d_newParticles,
                                    *d_domainReductionOffsets, *(d_cellIndicesLeft[speciesIndex]),
                                    *(individualPositions[speciesIndex]),*d_tempKeys,
                                    *d_reductionOffsetsAdding,
                                    *d_reductionOffsets,
                                    m_nIndividuals[speciesIndex],
                                    eolNew,
                                    *(individualProperties[speciesIndex]),
                                    *d_state,
                                    propertiesOffset,
                                    speciesIndex,
                                    seed
                                    );
        seed++;
    }


    // update total particle count
    m_nIndividuals[speciesIndex] = h_total[0] + h_totalAdded[0];

    if (lost > 0) {
        // append lost particles at the end of the list
        m_kernel_writeLostParticles.setEnqueueArgs(cl::EnqueueArgs(*m_commandQueue,
                                                                   cl::NDRange(lost)));
        m_kernel_writeLostParticles(*d_tempKeys,
                                    *d_lostParticles, m_nIndividuals[speciesIndex]);

    }
    //writePositions(0, "positionsAfterLostParticlesWriting.dat");

    // now copy the temp keys over to the proper key list
    m_kernel_copyTempToKeys.setEnqueueArgs(cl::EnqueueArgs(*m_commandQueue,
                                                           cl::NDRange(c_maxNumIndividuals)));
    m_kernel_copyTempToKeys(*d_tempKeys, *(individualPositions[speciesIndex]));

    // we finally need to rebuild the cell indices
    // todo: we can probably do that in the scatter kernel..
    computeCellIndices(speciesIndex, m_nIndividuals[speciesIndex]);

    //writePositions(0, "positionsAfterAll.dat");
    //exit(1);
    //writeReductionOffsets("reductionOffsetsAfter2.dat");
    //writeCellIndices(0, "cellIndicesAfterAll.dat");
    //exit(1);

    // check for consistency
    try {
        checkPositionListConsistency(speciesIndex, d_state);
    }
    catch (gpgmp::Error& e) {
        // write old particle list
        std::ofstream myfile;
        std::string fn;
        std::ostringstream oss;
        oss << "oldPositions_species_" << speciesIndex;
        fn = oss.str() + ".dat";

        uint nwrite=c_maxNumIndividuals;

        myfile.open (fn.c_str());
        for (uint i=0; i<nwrite; i++) {
            myfile << h_iposOld[i] << " " << GET_INDEX(h_iposOld[i]) << " " << GET_KEY(h_iposOld[i])
                   << "\n";
        }
        myfile.close();

        throw e;
    }
}

void IndividualSolver::sortSpecies(int species, uint nIndividuals)
{
    // only sort if there's more than one individual
    if (nIndividuals == 0) return;

    // first compute the sorting topology

    // number of sort keys - needs to be:
    // * bigger than 1024
    // * an integer exponent of 2, i.e. 2*l = nkeys (scan exclusive)
    // * a multiple of 4*c_ctaSize (total blocks for sorting)
    // * nkeys > m_wgSizeScan*c_ctaSize (scan exclusive)

    //uint nkeys = iSnapUp(nIndividuals, 1024);
    uint nkeys = std::max(1024, 4*c_ctaSize);
    nkeys=std::max(nkeys, uint(m_wgSizeScan*c_ctaSize));
    while (nIndividuals > nkeys) {
        nkeys <<=1;
    }
    //std::cout <<"Using "<<nkeys<<" keys for sorting.\n";

    // number of total blocks
    uint totalBlocks = nkeys/4/c_ctaSize;

    // we need to pad ..
    // todo: once we implement reactions for individuals we can do the
    // padding in the stream compaction step..
    if (nkeys > nIndividuals) {
        m_kernel_padKeys.setEnqueueArgs(cl::EnqueueArgs(*m_commandQueue,
                                                        cl::NDRange(nkeys-nIndividuals)));
        m_kernel_padKeys(*(individualPositions[species]), nkeys);
    }

    // sorts individual species using radix sort
    int i=0;
    while (i*c_bitStep < c_keyBits) {
        radixSortStep(i*c_bitStep, species, nkeys, totalBlocks);
        i++;
    }
}

void IndividualSolver::radixSortStep(uint startbit, int species, uint nkeys, uint totalBlocks) {

    // get number of keys to sort
    // todo: we should keep count of the number of species if we include reactions later
    // for now, we just take the max number of keys
    //uint nkeys = c_maxNumIndividuals;

    // todo: we don't need to compute sorting topology in each step..
    //uint totalBlocks = nkeys/4/c_ctaSize;

    uint bitStep = c_bitStep;
    m_kernel_radixSortBlocksKeysOnly.setEnqueueArgs(cl::EnqueueArgs(*m_commandQueue,
                                                                    cl::NDRange(c_ctaSize*totalBlocks),
                                                                    cl::NDRange(c_ctaSize)));
    m_kernel_radixSortBlocksKeysOnly.instance().setArg(6, 4*c_ctaSize*sizeof(type_key), 0);
    m_kernel_radixSortBlocksKeysOnly(*(individualPositions[species]),
                                     *d_tempKeys,
                                     bitStep, startbit, nkeys, totalBlocks);

    totalBlocks = nkeys/2/c_ctaSize;
    m_kernel_findRadixOffsets.setEnqueueArgs(cl::EnqueueArgs(*m_commandQueue,
                                                             cl::NDRange(c_ctaSize*totalBlocks),
                                                             cl::NDRange(c_ctaSize)));
    m_kernel_findRadixOffsets.instance().setArg(6, 2*c_ctaSize * sizeof(type_key), 0);
    m_kernel_findRadixOffsets(*d_tempKeys, *d_mCounters, *d_mBlockOffsets,
                              startbit, nkeys, totalBlocks);



    // scan over counters
    scanExclusive(d_mCountersSum, d_mCounters, 1, nkeys/2/c_ctaSize*16);


    // reorder the keys
    totalBlocks = nkeys/2/c_ctaSize;
    m_kernel_reorderDataKeysOnly.setEnqueueArgs(cl::EnqueueArgs(*m_commandQueue,
                                                                cl::NDRange(c_ctaSize*totalBlocks),
                                                                cl::NDRange(c_ctaSize)));
    m_kernel_reorderDataKeysOnly.instance().setArg(8, 2*c_ctaSize*sizeof(type_key),0);
    m_kernel_reorderDataKeysOnly(*(individualPositions[species]),
                                 *d_tempKeys, *d_mBlockOffsets, *d_mCountersSum,
                                 *d_mCounters, startbit, nkeys, totalBlocks);

} // one radix sort step

void IndividualSolver::scanExclusive(cl::Buffer *d_Dst, cl::Buffer *d_Src, uint batchSize, uint arrayLength) {
    //Check power-of-two factorization
    unsigned int log2L;
    unsigned int factorizationRemainder = factorRadix2(log2L, arrayLength);
    if (factorizationRemainder != 1) {
        // make sure that arrayLength can be written as 2**log2L
      std::cerr <<"ERROR: factorizationRemainder == "<<factorizationRemainder<< " (should be 1)."
               << "arraylength = "<<arrayLength<<", log2L = "<<log2L <<"\n";
      exit(1);
    }

    // check if array is in bounds
    if (!( (arrayLength >= 8 * m_wgSizeScan) && (arrayLength <= 4 * m_wgSizeScan * m_wgSizeScan))) {
      std::cerr <<"ERROR: array length out of bounds. It is "<< arrayLength
            <<" but should be in between "<< 8 * m_wgSizeScan
            <<" and "<< 4 * m_wgSizeScan * m_wgSizeScan <<"\n";
      exit(1);
    }
    if (! ( (batchSize * arrayLength) <= 64 * 1048576 /*TODO: where does this number come from? .. NVIDIA SDK .. */)) {
      std::cerr <<"ERROR: batch size out of bounds.\n";
      exit(1);
    }

    // scan over field
    uint cel1n = (batchSize * arrayLength) / (4 * m_wgSizeScan);
    uint cel1size = 4 * m_wgSizeScan;
    m_kernel_scanExclusiveLocal1.setEnqueueArgs(cl::EnqueueArgs(*m_commandQueue,
                                                                cl::NDRange(cel1n*cel1size/4),
                                                                cl::NDRange(m_wgSizeScan)));
    m_kernel_scanExclusiveLocal1.instance().setArg(3, 2*m_wgSizeScan*sizeof(cl_uint), 0);
    m_kernel_scanExclusiveLocal1(*d_Dst, *d_Src, cel1size);


    // scan over blocks
    uint sel2n = batchSize;
    uint sel2size = arrayLength / (4 * m_wgSizeScan);
    uint elements = sel2n * sel2size;
    m_kernel_scanExclusiveLocal2.setEnqueueArgs(cl::EnqueueArgs(*m_commandQueue,
                                                                cl::NDRange(iSnapUp(elements, m_wgSizeScan)),
                                                                m_wgSizeScan));
    m_kernel_scanExclusiveLocal2.instance().setArg(5, 2*m_wgSizeScan*sizeof(cl_uint), 0);
    m_kernel_scanExclusiveLocal2(*d_buffer, *d_Dst, *d_Src, elements, sel2size);

    // and update
    uint uun = (batchSize * arrayLength) / (4 * m_wgSizeScan);
    m_kernel_uniformUpdate.setEnqueueArgs(cl::EnqueueArgs(*m_commandQueue,
                                                          cl::NDRange(uun*m_wgSizeScan),
                                                          cl::NDRange(m_wgSizeScan)));
    m_kernel_uniformUpdate(*d_Dst, *d_buffer);
} // scan exclusive

unsigned int IndividualSolver::factorRadix2(unsigned int& log2L, unsigned int L)
{
    /* This function finds the first non-zero bit in L, stores numbers of leading
      zeros in log2L and returns L >> log2L
      */
  if(!L)
    {
      log2L = 0;
      return 0;
    } else {
    for(log2L = 0; (L & 1) == 0; L >>= 1, log2L++);
    return L;
  }
}

unsigned int IndividualSolver::iSnapUp(unsigned int dividend, unsigned int divisor)
{
    return ((dividend % divisor) == 0) ? dividend : (dividend - dividend % divisor + divisor);
}

void IndividualSolver::segmentedScan(cl::Buffer *d_data, cl::Buffer *d_flags, uint arrayLengthX, uint arrayLengthY)
{
    // performs a general segmented scan

    // 2D topology
    uint nscan = arrayLengthX*arrayLengthY/2; // total number of items to scan
    uint wgSizeX = c_ssWgSize;
    uint wgSizeY;
    if (arrayLengthY==1) {
        wgSizeY = 1; // 1D
    } else {
        wgSizeY = c_ssWgSize; // 2D
    }

    uint nBlocksX = arrayLengthX/2/wgSizeX;
    uint nBlocksY = arrayLengthY/wgSizeY;

    uint nBlocks = nBlocksX*nBlocksY;
    uint arrayLengthBlock = 2*wgSizeX*wgSizeY;

    // write out topology
    std::cout <<"Segmented scan has topology:\n";
    std::cout <<"nx = " << arrayLengthX/2 << ", ny = "<<arrayLengthY << ", total = "<< nscan <<"\n";
    std::cout <<"nbx = " << nBlocksX << ", nby = " << nBlocksY << ", total = "<< nBlocks << "\n";
    std::cout <<"wgx = " <<wgSizeX <<", wgy = " << wgSizeY <<"\n";

    //Check power-of-two factorization
    unsigned int log2L;
    unsigned int factorizationRemainder = factorRadix2(log2L, 2*nscan);
    if (factorizationRemainder != 1) {
        // make sure that arrayLength can be written as 2**log2L
      std::cerr <<"ERROR: factorizationRemainder == "<<factorizationRemainder<< " (should be 1)."
               << "arraylength = "<<2*nscan<<", log2L = "<<log2L <<"\n";
      exit(1);
    }

    // array length x/y need to be evenly divisible by c_ssWgSize
    if (! (arrayLengthX % wgSizeX == 0) ||
          ! (arrayLengthY % wgSizeY == 0)) {
        std::cout << std::flush;
        std::cerr <<"ERROR: segmented scan 2D topology mismatch.\n"
                 << "We have lx = "<< arrayLengthX <<", wgSizeX = " << wgSizeX
                 << ", ly = " << arrayLengthY <<" , wgSizeY = " << wgSizeY <<"\n";
        exit(1);
    }

    // check
    if (! (nBlocksX*wgSizeX*nBlocksY*wgSizeY == nscan)) {
        std::cout << std::flush;
        std::cout <<"ERROR: segmented scan computes something wrong. \n"
                 << "We have lx = "<< arrayLengthX <<", wgSizeX = " << wgSizeX
                 << ", ly = " << arrayLengthY <<" , wgSizeY = " << wgSizeY
                 << "\nTotal: " << nBlocksX*wgSizeX*nBlocksY*wgSizeY
                 << " should be " << nscan << "\n";
        exit(1);
    }

    /*
    // create buffers to hold block values
    cl::Buffer *d_dataBlocks = oclCreateBuffer<cl_uint>(nBlocks);
    cl::Buffer *d_flagsBlocks = oclCreateBuffer<cl_uchar>(nBlocks);
    cl::Buffer *d_firstBlockFlag= oclCreateBuffer<cl_uchar>(nBlocks);

    // buffer to hold the partial OR data for the array
    cl::Buffer *d_partialFlags = oclCreateBuffer<cl_uchar>(arrayLength);
    */

    // do up-sweep
    cl::NDRange usglobal(arrayLengthX/2, arrayLengthY);
    cl::NDRange uslocal(wgSizeX, wgSizeY);
    std::cout <<"Performing up-sweep with global="; printRange(std::cout, usglobal); std::cout <<"local = "; printRange(std::cout, uslocal); std::cout <<".\n"<<std::flush;
    m_kernel_segmentedScanUpSweep.setEnqueueArgs(cl::EnqueueArgs(*m_commandQueue,
                                                                 usglobal,
                                                                 uslocal));

    m_kernel_segmentedScanUpSweep.instance().setArg(6, arrayLengthBlock*sizeof(cl_uint),0);
    m_kernel_segmentedScanUpSweep.instance().setArg(7, arrayLengthBlock*sizeof(cl_uchar),0);
    m_kernel_segmentedScanUpSweep(*d_data, *d_flags, *d_ssPartialFlags,
                                  *d_ssDataBlocks, *d_ssFlagsBlocks, *d_ssFirstBlockFlag);

    // do second-level scan
    cl::NDRange slglobal(nBlocksX/2, nBlocksY);
    cl::NDRange sllocal(nBlocksX/2, nBlocksY);
    std::cout <<"Performing second-level sweep with global="; printRange(std::cout, slglobal); std::cout <<"local = "; printRange(std::cout, sllocal); std::cout <<".\n"<<std::flush;
    m_kernel_segmentedScanBlock.setEnqueueArgs(cl::EnqueueArgs(*m_commandQueue,
                                                               slglobal,
                                                               sllocal));
    m_kernel_segmentedScanBlock.instance().setArg(3, nBlocks*sizeof(cl_uint), 0);
    m_kernel_segmentedScanBlock.instance().setArg(4, nBlocks*sizeof(cl_uchar), 0);
    m_kernel_segmentedScanBlock.instance().setArg(5, nBlocks*sizeof(cl_uchar), 0);
    m_kernel_segmentedScanBlock(*d_ssDataBlocks, *d_ssFlagsBlocks, *d_ssFirstBlockFlag);

    // do down-sweep
    std::cout <<"Performing down-sweep with global="; printRange(std::cout, usglobal); std::cout <<"local = "; printRange(std::cout, uslocal); std::cout <<".\n"<<std::flush;
    m_kernel_segmentedScanDownSweep.setEnqueueArgs(cl::EnqueueArgs(*m_commandQueue,
                                                                   usglobal,
                                                                   uslocal));
    m_kernel_segmentedScanDownSweep.instance().setArg(5, arrayLengthBlock*sizeof(cl_uint),0);
    m_kernel_segmentedScanDownSweep.instance().setArg(6, arrayLengthBlock*sizeof(cl_uchar),0);
    m_kernel_segmentedScanDownSweep.instance().setArg(7, arrayLengthBlock*sizeof(cl_uchar),0);
    m_kernel_segmentedScanDownSweep(*d_data, *d_flags, *d_ssPartialFlags,
                                    *d_ssDataBlocks, *d_ssFlagsBlocks);

    /*
    // do single-block scan
    m_kernel_segmentedScanBlock.setEnqueueArgs(cl::EnqueueArgs(*m_commandQueue,
                                                               cl::NDRange(nscan),
                                                               cl::NDRange(nscan)));
    m_kernel_segmentedScanBlock.instance().setArg(3, arrayLength*sizeof(cl_uint), 0);
    m_kernel_segmentedScanBlock.instance().setArg(4, arrayLength*sizeof(cl_uchar), 0);
    m_kernel_segmentedScanBlock.instance().setArg(5, arrayLength*sizeof(cl_uchar), 0);
    m_kernel_segmentedScanBlock(*d_data, *d_flags, *d_flags);
    */

    /*
    // delete temp buffers
    delete d_dataBlocks;
    delete d_flagsBlocks;
    delete d_firstBlockFlag;
    delete d_partialFlags;*/
}

void IndividualSolver::testSegmentedScan()
{
    // tests the segmented scan routine
    // we take uint as data buffer and char as flag buffer
    uint nelements = m_diffusionModel.gridArea(); // needs to be 2** ..
    //uint nea = 512;
    //uint nelements = nea*nea; // 1D

    uint nBlocksX = m_diffusionModel.gridWidth()/2/c_ssWgSize;
    uint nBlocksY = m_diffusionModel.gridHeight()/c_ssWgSize;
    //uint nBlocksX = nelements/2/c_ssWgSize; // 1D
    //uint nBlocksY = 1; // 1D
    //uint nBlocksX = nea/2/c_ssWgSize; // 2D
    //uint nBlocksY = nea/c_ssWgSize; // 2D
    uint nBlocks = nBlocksX*nBlocksY;

    std::vector<cl_uint> h_data(nelements);
    std::vector<cl_uchar> h_flags(nelements);
    std::vector<cl_uint> h_datas(nelements); // for single-block scan
    std::vector<cl_uchar> h_flagss(nelements); // for single-block scan

    // use random data
    boost::random::mt19937 rng;
    rng.seed(m_diffusionModel.getRandomSeed());
    //rng.seed(0);

    boost::random::uniform_int_distribution<> dist(0, 25);
    boost::random::uniform_int_distribution<> dist2(0, nelements-1);


    // set data and erase flags
    for (uint i=0; i<nelements; i++) {
        h_flags[i] = 0;
        h_flagss[i] = 0;
        h_data[i]=dist(rng);
        h_datas[i] = h_data[i];
    }

    // randomly set a few flags
    uint nflags = 25;
    for (uint i=0; i<nflags; i++) {
        uint t = dist2(rng);
        h_flags[t] = 1;
        h_flagss[t] = 1;
    }    

    // save both
    std::ofstream myfile;

    myfile.open("segscan_initial.dat");
    for (uint i=0; i<nelements; i++) {
        myfile << i << " " << h_data[i] << " " << (uint) h_flags[i] /* << " " << h_datas[i] << " " << (uint) h_flagss[i]*/  << "\n";
    }
    myfile.close();

    // create buffers from them
    cl::Buffer *d_data = oclCreateBufferFrom(h_data);
    cl::Buffer *d_flags = oclCreateBufferFrom(h_flags);
    cl::Buffer *d_datas = oclCreateBufferFrom(h_datas);
    cl::Buffer *d_flagss = oclCreateBufferFrom(h_flagss);

    // fill all temp buffers with traceable vaules
    std::vector<cl_uint> db(nBlocks);
    std::vector<cl_uchar> fb(nBlocks);
    std::vector<cl_uchar> fbf(nBlocks);
    std::vector<cl_uchar> pf(nelements);
    for (uint i=0; i<nBlocks; i++) {
        db[i] = 666;
        fb[i] = 222;
        fbf[i] = 111;
    }
    for (uint i=0; i<nelements; i++) pf[i] = 212;
    oclCopy(db, *d_ssDataBlocks);
    oclCopy(fb, *d_ssFlagsBlocks);
    oclCopy(fbf, *d_ssFirstBlockFlag);
    oclCopy(pf, *d_ssPartialFlags);

    /*
    // do single-block scan - this works in 2D
    m_kernel_segmentedScanBlock_single.setEnqueueArgs(cl::EnqueueArgs(*m_commandQueue,
                                                               cl::NDRange(32/2, 32),
                                                               cl::NDRange(32/2, 32)));
    m_kernel_segmentedScanBlock_single.instance().setArg(3, nelements*sizeof(cl_uint), 0);
    m_kernel_segmentedScanBlock_single.instance().setArg(4, nelements*sizeof(cl_uchar), 0);
    m_kernel_segmentedScanBlock_single.instance().setArg(5, nelements*sizeof(cl_uchar), 0);
    m_kernel_segmentedScanBlock_single(*d_datas, *d_flagss, *d_flagss);
    */

    // do blocked scan
    segmentedScan(d_data, d_flags, m_diffusionModel.gridWidth(),m_diffusionModel.gridHeight());
    //segmentedScan(d_data, d_flags, nea,nea);
    //segmentedScan(d_data, d_flags, nelements, 1); // 1D

    // write back to host
    oclCopy(*d_data, h_data);
    oclCopy(*d_flags, h_flags);
    oclCopy(*d_datas, h_datas);
    oclCopy(*d_flagss, h_flagss);

    // copy partial results
    std::cout <<" Copying topology is nbx="<<nBlocksX<<", nby="<<nBlocksY<<".\n";
    std::vector<cl_uchar> h_partialFlags(nelements);
    std::vector<cl_uint> h_dataBlocks(nBlocks);
    std::vector<cl_uchar> h_flagsBlocks(nBlocks);
    std::vector<cl_uchar> h_firstBlockFlag(nBlocks);

    oclCopy(*d_ssDataBlocks, h_dataBlocks);
    oclCopy(*d_ssFlagsBlocks, h_flagsBlocks);
    oclCopy(*d_ssFirstBlockFlag, h_firstBlockFlag);
    oclCopy(*d_ssPartialFlags, h_partialFlags);

    // and output to file
    myfile.open("segscan_final.dat");
    for (uint i=0; i<nelements; i++) {
        myfile << i << " "  /* << h_datas[i] << " " << (uint) h_flagss[i] << " " */ << h_data[i] << " " << (uint) h_flags[i] << " " << (uint) h_partialFlags[i] << "\n";
    }
    myfile.close();

    // output block results
    myfile.open("segscan_blocks_final.dat");
    for (uint i=0; i<nBlocks; i++) {
        myfile << i << " " << h_dataBlocks[i] << "  " << (uint) h_flagsBlocks[i] << " " << (uint) h_firstBlockFlag[i] << "\n";
    }
    myfile.close();

    // delete buffer
    delete d_data;
    delete d_flags;
    exit(1);
}

void IndividualSolver::computeCellIndices(int speciesIndex, uint nIndividuals) {

    // only do something if there's any individuals
    if (nIndividuals ==0 ) return;

    uint numCells = m_diffusionModel.gridArea();
    uint nkeys = nIndividuals;

    m_kernel_computeCellIndices.setEnqueueArgs(cl::EnqueueArgs(*m_commandQueue,
                                                               cl::NDRange(nkeys)));

    m_kernel_computeCellIndices(*(individualPositions[speciesIndex]),
                                *(d_cellIndicesLeft[speciesIndex]),
                                *(d_cellIndicesRight[speciesIndex]),
                                nkeys, numCells);
}

void IndividualSolver::computeDiffusionConstants(int speciesIndex,
                                                 cl::Buffer *d_diffusionConstantsX,
                                                 cl::Buffer *d_diffusionConstantsY,
                                                 cl::Buffer *d_rx, cl::Buffer *d_ry,
                                                 cl::Buffer *d_state,
                                                 Real simTime, cl::Buffer *d_errors,
                                                 cl::Buffer *d_sumMoments,
                                                 cl::Buffer *d_sumStates) {

    //std::cout <<"Start computing diffusion constants for species "<<speciesIndex<<"..\n"<<std::flush;

    // get number of individuals
    uint nIndividuals = m_nIndividuals[speciesIndex];

    //std::cout <<"Sorting "<<nIndividuals<<".. "<<std::flush;
    // first we need to sort the species
    sortSpecies(speciesIndex, nIndividuals);
    //std::cout <<"done.\n"<<std::flush;

    // then we need to compute the cell indices
    //std::cout <<"Computing cell indices ..\n"<<std::flush;
    computeCellIndices(speciesIndex, nIndividuals);
    //std::cout <<"done.\n"<<std::flush;

    // and now we can compute the constants
    uint width = m_diffusionModel.gridWidth();
    uint height = m_diffusionModel.gridHeight();
    const uint propertiesOffset = c_maxNumIndividuals;

    m_kernel_computeDiffusionConstants.setEnqueueArgs(cl::EnqueueArgs(*m_commandQueue,
                                                                      cl::NDRange(width, height),
                                                                      cl::NDRange(16,16)));
    m_kernel_computeDiffusionConstants(*(individualPositions[speciesIndex]),
                                       *(d_cellIndicesLeft[speciesIndex]),
                                       *(d_cellIndicesRight[speciesIndex]),
                                       *(individualProperties[speciesIndex]),
                                       propertiesOffset,
                                       *d_diffusionConstantsX, *d_diffusionConstantsY,
                                       *d_rx, *d_ry,
                                       *d_state,
                                       simTime, speciesIndex,
                                       *d_sumMoments, *d_sumStates
                                       );

/*
    m_commandQueue->flush();
    writeAll(speciesIndex);
    writeDiffusionConstants("diffusionConstants.dat",
                            d_diffusionConstantsX, d_diffusionConstantsY,
                            d_rx, d_ry);
    writeProperties("properties.dat", speciesIndex);
    exit(1);
*/
} // compute drift+diffusion constants

void IndividualSolver::writeAll(int speciesIndex)
{
    writePositions(speciesIndex, "positions", 0, c_maxNumIndividuals);
    writeCellIndices(speciesIndex, "cellIndices");
    writeReductionOffsets("reductionOffsets.dat");
    writeSegmentedScanBuffers("scanBuffers.dat");
}

void IndividualSolver::writeProperties(std::string filename, int speciesIndex) {
    // copy over to host
    std::vector<Real> h_props(4*c_maxNumIndividuals); // todo: might be more than 2 ..
    oclCopy(*(individualProperties[speciesIndex]), h_props);
    std::ofstream myfile;
    myfile.open (filename.c_str());
    for (uint i=0; i<c_maxNumIndividuals; i++) {
        myfile << i <<" " << h_props[i+0*c_maxNumIndividuals] << " " << h_props[i+1*c_maxNumIndividuals]
                  << " " <<h_props[i+2*c_maxNumIndividuals] << " " << h_props[i+3*c_maxNumIndividuals] << "\n";
    }
    myfile.close();
}

void IndividualSolver::writeDiffusionConstants(std::string filename,
                                               cl::Buffer *d_diffusionConstantsX, cl::Buffer *d_diffusionConstantsY,
                                               cl::Buffer *d_rx, cl::Buffer *d_ry) {

    // TODO: !! THIS DOES NOT TAKE INTO ACCOUNT THE SPECIES INDEX !!
    // copy over to host
    std::vector<Real> h_diffx(m_diffusionModel.gridArea());
    std::vector<Real> h_diffy(m_diffusionModel.gridArea());
    std::vector<Real> h_rx(m_diffusionModel.gridArea());
    std::vector<Real> h_ry(m_diffusionModel.gridArea());

    oclCopy(*d_diffusionConstantsX, h_diffx);
    oclCopy(*d_diffusionConstantsY, h_diffy);
    oclCopy(*d_rx, h_rx);
    oclCopy(*d_ry, h_ry);

    std::ofstream myfile;
    myfile.open (filename.c_str());

    for (uint i=0; i<m_diffusionModel.gridArea(); i++) {
        myfile << i <<" "<<h_diffx[i] <<" " << h_diffy[i]
               << " " << h_rx[i] << " " << h_ry[i] <<"\n";
    }

    myfile.close();
}

void IndividualSolver::writePositions(int speciesIndex, std::string filename, int index, int num)
{
    // copy back to device and print
    std::vector<cl_ulong> h_ipos(c_maxNumIndividuals);
    oclCopy(*(individualPositions[speciesIndex]), h_ipos);
    std::vector<cl_ulong> h_iposTemp(c_maxNumIndividuals);
    oclCopy(*d_tempKeys, h_iposTemp);

    std::ofstream myfile;

    std::string fn;
    std::ostringstream oss;
    oss << "_species_" << speciesIndex;
    if (index > 0) {
        oss << "_"<< index;
    }
    fn = filename +  oss.str() + ".dat";

    uint nwrite=c_maxNumIndividuals;
    if (num==0)
        nwrite = m_nIndividuals[speciesIndex];

    myfile.open (fn.c_str());
    for (uint i=0; i<nwrite; i++) {
        myfile << h_ipos[i] << " " << GET_INDEX(h_ipos[i]) << " " << GET_KEY(h_ipos[i])
               << " " << GET_INDEX(h_iposTemp[i]) << " " << GET_KEY(h_iposTemp[i])
               << "\n";
    }
    myfile.close();
}

void IndividualSolver::writeCellIndices(int speciesIndex, std::string filename)
{
    std::string fn;
    std::ostringstream oss;
    oss << "_species_" << speciesIndex;
    fn = filename +  oss.str() + ".dat";

    // create host buffer
    std::vector<cl_uint> left(m_diffusionModel.gridArea());
    std::vector<cl_uint> right(m_diffusionModel.gridArea());

    // copy back to device and print
    oclCopy(*(d_cellIndicesLeft[speciesIndex]), left);
    oclCopy(*(d_cellIndicesRight[speciesIndex]), right);

    std::ofstream myfile;
    myfile.open (filename.c_str());
    for (uint i=0; i<m_diffusionModel.gridArea(); i++) {
        myfile << i <<" " << left[i] <<" " <<right[i] << " \n";
    }
    myfile.close();
}

void IndividualSolver::writeReductionOffsets(std::string filename, int index)
{
    // create host buffer
    std::vector<int> h_offsets(c_maxNumIndividuals);
    std::vector<int> h_offsetsBlocks(c_maxNumIndividuals/2/16);
    std::vector<int> h_offsetsAdding(c_maxNumIndividuals);
    std::vector<int> h_offsetsBlocksAdding(c_maxNumIndividuals/2/16);

    // copy back to host and print
    oclCopy(*(d_reductionOffsets), h_offsets);
    oclCopy(*(d_reductionOffsetsBlocks), h_offsetsBlocks);
    oclCopy(*(d_reductionOffsetsAdding), h_offsetsAdding);
    oclCopy(*(d_reductionOffsetsBlocksAdding), h_offsetsBlocksAdding);

    std::string fn;
    if (index > 0) {
        std::ostringstream oss;
        oss << index;
        fn = filename + "." + oss.str();
    } else {
        fn = filename;
    }

    std::ofstream myfile;
    myfile.open (fn.c_str());
    for (uint i=0; i<c_maxNumIndividuals; i++) {
        myfile << i <<" "<< h_offsets[i]<< " " << h_offsetsAdding[i] <<"\n";
    }
    myfile.close();

    myfile.open ("offsetBlocks.dat");
    for (uint i=0; i<c_maxNumIndividuals/2/16; i++) {
        myfile << i <<" "<< h_offsetsBlocks[i] << " " << h_offsetsBlocksAdding[i] <<"\n";
    }
    myfile.close();

}

std::vector<cl_ulong> IndividualSolver::getPositionList(int speciesIndex)
{
    // copy positions back to device
    std::vector<cl_ulong> h_ipos(m_nIndividuals[speciesIndex]);
    if (m_nIndividuals[speciesIndex] > 0)
        oclCopy(*(individualPositions[speciesIndex]), h_ipos);

    return h_ipos;
}

std::map<int, std::vector<Real> > IndividualSolver::getPropertiesList(int speciesIndex) {
    // copy whole list over to host
    std::vector<Real> h_props(4*c_maxNumIndividuals);
    oclCopy(*(individualProperties[speciesIndex]), h_props);

    // now we need to separate the properties
    std::map<int, std::vector<Real> > props;
    props[0] = std::vector<Real>(h_props.begin(), h_props.begin()+m_nIndividuals[speciesIndex]);
    props[1] = std::vector<Real>(h_props.begin()+c_maxNumIndividuals, h_props.begin()+c_maxNumIndividuals+m_nIndividuals[speciesIndex]);
    props[2] = std::vector<Real>(h_props.begin()+2*c_maxNumIndividuals, h_props.begin()+2*c_maxNumIndividuals+m_nIndividuals[speciesIndex]);
    props[3] = std::vector<Real>(h_props.begin()+3*c_maxNumIndividuals, h_props.begin()+3*c_maxNumIndividuals+m_nIndividuals[speciesIndex]);

    return props;
}

void IndividualSolver::checkPositionListConsistency(int speciesIndex, cl::Buffer *d_state)
{
    try {
        // copy positions back to device
        std::vector<cl_ulong> h_ipos(c_maxNumIndividuals);
        oclCopy(*(individualPositions[speciesIndex]), h_ipos);

        // copy state back to device
        std::vector<cl_uint> h_data(m_diffusionModel.gridArea()*m_diffusionModel.numSpecies());
        oclCopy(*d_state, h_data);

        // vector to count the occurrences of each individual
        std::vector<unsigned int> occ(c_maxNumIndividuals,0);

        // vector to count the "is" state to compare to target state
        std::vector<uint> isState(m_diffusionModel.gridArea(), 0);

        std::cerr <<"Performing consistency check for individuals of species "
                 << speciesIndex << std::flush;

        for (uint i=0; i<c_maxNumIndividuals; i++) {
            // check for multiple particles
            occ[GET_INDEX(h_ipos[i])] += 1;
            if (occ[GET_INDEX(h_ipos[i])] > 1) {
                std::cerr <<"ERROR in individual consistency check: "
                            "Particle " << GET_INDEX(h_ipos[i]) << " (" << i << ") "
                            " of species " << speciesIndex;
                std::cerr << "is a duplicate. Exiting.. \n"<<std::flush;
                writeAll(speciesIndex);
                exit(1);
            }
            // check if they are in order
            if (i>0) {
                if (GET_KEY(h_ipos[i]) < GET_KEY(h_ipos[i-1])) {
                    std::cerr <<"ERROR in individual consistency check: "
                                "Particle " << GET_INDEX(h_ipos[i]) << " (" << i << ") "
                                " of species " << speciesIndex
                             << "is not in order. Cell index is " << GET_KEY(h_ipos[i]) <<
                                " which is bigger than previous particle " << GET_KEY(h_ipos[i-1]);
                    throw gpgmp::Error(-1, "Failed consistency check");
                }
            }
            // add to is state
            if (GET_KEY(h_ipos[i]) < UINT_MAX)
                isState[GET_KEY(h_ipos[i])]++;
        }

        std::cout <<"Performing additional state check .."<<std::flush;

        // finally compare is state to target state
        for (uint i=0; i<m_diffusionModel.gridArea(); i++) {
            if (isState[i] != h_data[speciesIndex*(m_diffusionModel.gridArea())+i]) {
                std::cerr << "ERROR in individual consistency check: "
                             "We have " << isState[i] <<
                             " in cell " << i << " but it should be "
                          <<  h_data[speciesIndex*(m_diffusionModel.gridArea())+i] <<". Exiting. \n"<<std::flush;
                writeAll(speciesIndex);
                exit(1);
            }
        }
        std::cerr <<"done.\n" << std::flush;
    } catch (gpgmp::Error& e) {
        writeAll(speciesIndex);
        writeCellIndices(0, "cellIndices.dat");
        m_diffusionModel.writeAllSpecies(std::numeric_limits<float>::max(), this);
        throw e;
    } // catch errors
}

void IndividualSolver::writeLostParticles(std::string filename)
{
    // copy back to device and print
    std::vector<cl_ulong> h_lost(c_maxNumIndividuals);
    oclCopy(*d_lostParticles, h_lost);
    std::ofstream myfile;
    myfile.open (filename.c_str());
    for (uint i=0; i<c_maxNumIndividuals; i++) {
        myfile << h_lost[i] << " " << GET_INDEX(h_lost[i]) << " " << GET_KEY(h_lost[i]) << "\n";
    }
    myfile.close();
}

void IndividualSolver::writeSegmentedScanBuffers(std::string filename)
{
    // create host arrays
    uint nelements = m_diffusionModel.gridArea();
    std::vector<cl_uint> h_data(m_diffusionModel.gridArea());
    std::vector<cl_uint> h_reductionOffsets(m_diffusionModel.gridArea());
    std::vector<cl_uchar> h_flags(m_diffusionModel.gridArea());

    // write back to host
    std::cout <<"Reading back results ..\n"<<std::flush;
    oclCopy(*d_newParticles, h_data);
    oclCopy(*d_domainReductionOffsets, h_reductionOffsets);
    oclCopy(*d_ssFlags, h_flags);

    // and output to file
    std::ofstream myfile;
    myfile.open(filename.c_str());
    for (uint i=0; i<nelements; i++) {
        myfile << i << " "<< h_data[i] << " " << h_reductionOffsets[i] << " " << (uint) h_flags[i] /* << " " << (uint) h_partialFlags[i]*/ << "\n";
    }
    myfile.close();
}

void IndividualSolver::writeState(std::string filename, cl::Buffer *d_state)
{
    // create host arrays
    uint nelements = m_diffusionModel.gridArea();
    std::vector<cl_uint> h_data(m_diffusionModel.gridArea());

    // copy state
    oclCopy(*d_state, h_data);

    // and output to file
    std::ofstream myfile;
    myfile.open(filename.c_str());
    for (uint i=0; i<nelements; i++)
        myfile << i << " " << h_data[i] << "\n";

    myfile.close();
} // write buffers from segmented scan

void IndividualSolver::oclCopyDeviceBuffers(cl::Buffer &from, cl::Buffer &to, size_t count) {
    m_commandQueue->enqueueCopyBuffer(from, to, 0, 0, count);
    m_commandQueue->flush();
}

} // namespace gpgmp
