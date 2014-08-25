/*
 *  DiffusionModel.cpp
 *
 *  Created on: 09/12/2009
 *  Created by: matthias
 */

// TODO: add checks to ensure that gridWidth,gridHeight is valid (e.g. multiple of 16) - and square?
// TODO: add checks to ensure that ids are unique across species and parameters

#include "DiffusionModel.h"

#include "Kernel.h"
#include "Compartment.h"
#include "Species.h"
#include "IndividualSolver.h"

#include "HomogeneousSolver.h"
#include "InhomogeneousSolver.h"

#if defined(_WIN32) || defined(_WIN64)
# include "win_gettimeofday.h"
#endif

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/math/special_functions/round.hpp>
#include <boost/random.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>

#ifndef __APPLE__
#include <boost/regex.hpp>
#endif

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

// Important: must include Python AFTER the basic std and boost headers
#include <Python.h>
#include <boost/python.hpp>
#include <numpy/noprefix.h>

#include <boost/log/trivial.hpp>

// Misc Prototypes
void CL_CALLBACK oclNotification(const char *errinfo, const void *private_info, size_t cb, void  *user_data);

namespace gpgmp {

const Real DiffusionModel::Avogadro = 6.0221415e23f; ///< Avogadro number.

// Constructor
DiffusionModel::DiffusionModel(const std::string &id, bool managePythonInstance, unsigned int randomSeed)
:   Base(id),
    m_managesPythonInstance(managePythonInstance),
    m_nonlinearDiffusivity(false),
    m_computeMoments(false),
    m_hasContinuousParameterSpecies(false),
    m_hasConstantIndividuals(false),
    m_randomSeed(randomSeed)
{
    init();
}

DiffusionModel::DiffusionModel(Real length, int dx, int dy)
	: DiffusionModel()
{

}

// logs the model
void DiffusionModel::logModel()
{
    BOOST_LOG_TRIVIAL(debug) <<"Diffusionmodel "<<id()<<".width: "<<m_gridWidth;
    BOOST_LOG_TRIVIAL(debug) <<"Diffusionmodel "<<id()<<".heigh: "<<m_gridHeight;
    BOOST_LOG_TRIVIAL(debug) <<"Diffusionmodel "<<id()<<".physicalWidth: "<<m_physicalLength;
    unsigned int i=1;
    BOOST_FOREACH(Compartment *comp, m_compartments) {
        BOOST_LOG_TRIVIAL(debug)  <<"Diffusionmodel "<<id()<<".compartment "<<i++<<": "<<*comp;
    }
    i=1;
    BOOST_FOREACH(Species *species, m_species) {
        BOOST_LOG_TRIVIAL(debug)  <<"Diffusionmodel "<<id()<<".species "<<i++<<": "<<*species;
    }
    i=1;
    BOOST_FOREACH(Reaction *reaction, m_reactions) {
        BOOST_LOG_TRIVIAL(debug)  <<"Diffusionmodel "<<id()<<".reaction "<<i++<<": "<<*reaction;
    }

    // Output Solver Type
    switch (m_solver){
    case gpgmp::stochastic_homogeneous:
        BOOST_LOG_TRIVIAL(debug)  <<"Diffusionmodel "<<id()<<".solverType: STOCHASTIC_HOMOGENEOUS";
        break;

    case gpgmp::stochastic_inhomogeneous:
        BOOST_LOG_TRIVIAL(debug)  <<"Diffusionmodel "<<id()<<".solverType: STOCHASTIC_INHOMOGENEOUS";
        break;

    default:
        BOOST_LOG_TRIVIAL(debug)  <<"Diffusionmodel "<<id()<<".solverType: UNKNOWN";
        break;
    }

    // output drift-diffusion method
    if (m_solver == gpgmp::stochastic_inhomogeneous)
        BOOST_LOG_TRIVIAL(debug)  <<"Diffusionmodel "<<id()<<".computeDriftDiffusivityMethod: "<<computeDriftDiffusivityMethod();

    // set field method
    BOOST_LOG_TRIVIAL(debug) <<"Diffusionmodel "<<id()<<".setFieldMethod: "<<returnSetFieldMethod();

    // new individuals method
    BOOST_LOG_TRIVIAL(debug) <<"Diffusionmodel "<<id()<<".newIndividuals: "<<getNewIndividualsMethod();

    // output boundary conditions
    BOOST_LOG_TRIVIAL(debug)  <<"Diffusionmodel "<<id()<<".boundaryConditions: ["
                             << m_boundaryMask.s[0] << ", "
                             << m_boundaryMask.s[1] << ", "
                             << m_boundaryMask.s[2] << ", "
                             << m_boundaryMask.s[3] <<"].";

    // output nonlinear and compute moments flags
    if (m_nonlinearDiffusivity)
        BOOST_LOG_TRIVIAL(debug) <<"Diffusionmodel "<<id()<<".nonlinear: true";
    else
        BOOST_LOG_TRIVIAL(debug) <<"Diffusionmodel "<<id()<<".nonlinear: false";

    if (m_computeMoments)
        BOOST_LOG_TRIVIAL(debug) <<"Diffusionmodel "<<id()<<".computeMoments: true";
    else
        BOOST_LOG_TRIVIAL(debug) <<"Diffusionmodel "<<id()<<".computeMoments: false";

    BOOST_LOG_TRIVIAL(info) <<"Diffusionmodel "<<id()<<".randomSeed: "<<m_randomSeed;
}

// Initialisation
INIT_RETURN_VALUE DiffusionModel::init()
{
    m_physicalLength = 40.0f;
    m_gridWidth = 32;
    m_gridHeight = 32;
    m_solver = stochastic_homogeneous;
    m_p0 = 0.1f;
    m_numRuns = 1;
    m_maxSteps = InfSteps;
    m_maxTime = 100.0f;
    m_runsOffset = 0;
    m_oclSourcePath = boost::filesystem::current_path().string();
    m_oclDeviceIndex = 0;
    m_oclInhibitWarnings = false;
    m_outputInterval = 100.0f;
    m_outputFormat = gpgmp::OutputHdf5;
    m_outputFilename = "output.h5";

    m_regenerateInitStates = false;

    m_runIndex = 0;
    m_runDumpIndex = 0;
    
    m_hdf5File = 0;
    m_hdf5CurrentGroup = 0;
    
    if (m_managesPythonInstance) {
        // initialize main interpreter
        Py_Initialize();
        BOOST_LOG_TRIVIAL(trace) << "Interpreter init done.";
    }

    // init numpy arrays
    import_array();
    BOOST_LOG_TRIVIAL(trace) << "Importing of numpy arrays module.. done.";
    
    
    clearBoundaryMasks();

    m_hasBallisticBoundaryConditions = false;

    // set random seed if none is given
    m_randomSeed = (unsigned int)time(0);

#ifdef HAS_INIT_RETURN_VALUE
    return 0;
#endif
}

// Destructor
DiffusionModel::~DiffusionModel()
{
    if (m_managesPythonInstance) {
        Py_Finalize();
    }
    
    BOOST_FOREACH(Species *s, m_species) {
        delete s;
    }
    
    BOOST_FOREACH(Reaction *r, m_reactions) {
        delete r;
    }
    
    BOOST_FOREACH(Compartment *c, m_compartments) {
        delete c;
    }
}

void DiffusionModel::setPhysicalLength(Real length) {
    m_physicalLength = length;
}    
void DiffusionModel::setGridModelWidth(size_t width) {
    m_gridWidth = width;
}
void DiffusionModel::setGridModelHeight(size_t height) {
    m_gridHeight = height;
}
void DiffusionModel::setGridDims(size_t width, size_t height) {
    setGridModelWidth(width);
    setGridModelHeight(height);
}


/**
 *
 * @param boundaryMask The boundary mask.
 * @param sourceMask The source mask.
 */
void DiffusionModel::setBoundaryMasks(const cl_int4 &boundaryMask, const cl_int4 &sourceMask) {
    m_boundaryMask = boundaryMask;
    m_sourceMask   = sourceMask;
}
void DiffusionModel::clearBoundaryMasks() {
    for (int i=0; i<4; i++) {
        m_boundaryMask.s[i] = 0;
        m_sourceMask.s[i] = 0;
    }
}


void DiffusionModel::setSolver(SolverType solver) {
    m_solver = solver;
}


void DiffusionModel::setSolver(const std::string &solverName)
{
    SolverType solver = stochastic_inhomogeneous;
    
    // FUTURE: use a hash instead?
    if      (solverName == "stochastic_homogeneous")           solver = stochastic_homogeneous;
    else if (solverName == "stochastic_inhomogeneous")         solver = stochastic_inhomogeneous;
    else if (solverName == "stochastic_inhomogeneous_fpe")     solver = stochastic_inhomogeneous_fpe;
    else if (solverName == "deterministic_homogeneous")        solver = deterministic_homogeneous_RK4; // todo: remove this hack? it's mainly for inchman
    else if (solverName == "deterministic_homogeneous_RK4")    solver = deterministic_homogeneous_RK4;
    else if (solverName == "deterministic_homogeneous_aqss")   solver = deterministic_homogeneous_aqss;
    else if (solverName == "deterministic_homogeneous")        solver = deterministic_homogeneous_RK4; // todo: remove this hack? it's mainly for inchman
    else if (solverName == "deterministic_inhomogeneous_RK4")  solver = deterministic_inhomogeneous_RK4;
    else if (solverName == "deterministic_inhomogeneous_aqss") solver = deterministic_inhomogeneous_aqss;
    
    setSolver(solver);
}

void DiffusionModel::setP0(Real p0) {
    m_p0 = p0;
}


void DiffusionModel::setNumRuns(unsigned numRuns) {
    m_numRuns = numRuns;
}

void DiffusionModel::setMaxSteps(int maxSteps) {
    m_maxSteps = maxSteps;
}
void DiffusionModel::setMaxTime(Real maxTime) {
    m_maxTime = maxTime;
}



void DiffusionModel::setOclSourcePath(const std::string &path) {
    m_oclSourcePath = path;
}

void DiffusionModel::setOclDeviceIndex(size_t deviceIndex) {
    m_oclDeviceIndex = deviceIndex;
}
    
void DiffusionModel::setOclInhibitWarnings(bool inhibit) {
    m_oclInhibitWarnings = inhibit;
}

void DiffusionModel::setOutputInterval(Real interval) {
    m_outputInterval = interval;
}
void DiffusionModel::setOutputFormat(OutputFormat format) {
    m_outputFormat = format;
}
void DiffusionModel::setOutputFilename(const std::string &filename) {
    m_outputFilename = filename;
}

void DiffusionModel::setInitScriptContents(const std::string &contents) {
    m_initScriptContents = contents;
}

void DiffusionModel::loadInitScript(const std::string &filename) {
    loadTextFile(m_initScriptContents, filename);
}
void DiffusionModel::clearInitScript() {
    m_initScriptContents.clear();
}

void DiffusionModel::setComputeDriftDiffusivityMethod(const DiffusionType method) {
    std::ostringstream ss;

    switch(method) {
    case DT_HOMOGENEOUS:
        ss <<
            "// sets the diffusivity and drift to the constants given by the host parameters\n"
            "All->DiffusivityX = diffX;\n"
            "All->DiffusivityY = diffY;\n"
            "\n"
            "All->DriftX = driftX;\n"
            "All->DriftY = driftY;\n";
        break;

    case DT_MULTIPLICATIVE_NOISE:
        ss <<	  
            "Real shiftedX = PhysicalX + PhysicalModelWidth;\n"
            "Real shiftedY = PhysicalY + PhysicalModelHeight;\n"
            "\n"	  
            "All->DiffusivityX = (sigmaX * sigmaX) * (shiftedX * shiftedX) / 2.0;\n"
            "All->DiffusivityY = (sigmaY * sigmaY) * (shiftedY * shiftedY) / 2.0;\n"
            "\n"	  
            "All->DriftX = muX * shiftedX;\n"
            "All->DriftY = muY * shiftedY;\n";
	  
	  /*
	  //"Real xc = PhysicalCellWidth*((Real)get_global_id(0)-(Real)GridModelWidth/2.) + PhysicalModelWidth;\n"
	  //"Real yc = PhysicalCellWidth*((Real)get_global_id(1)-(Real)GridModelWidth/2.) + PhysicalModelWidth;\n"
	  "Real xc = getPhysicalX() + PhysicalModelWidth;\n"
	  "Real yc = getPhysicalY() + PhysicalModelWidth;\n"
	  "DiffusivityX = (sigmaX * sigmaX) * (xc * xc) / 2.0;\n"
	  "DiffusivityY = (sigmaY * sigmaY) * (yc * yc) / 2.0;\n"
	  "DriftX = muX * xc;\n"
	  "DriftY = muY * yc;\n";*/
        break;

    case DT_ORNSTEIN_UHLENBECK:
      ss <<
	"// set up an Ornstein-Uhlenbeck process\n"
    "All->DiffusivityX = diffX;\n"
    "All->DiffusivityY = diffY;\n"
	"\n"	
    "All->DriftX = -(gammaXX*PhysicalX + gammaXY*PhysicalY);\n"
    "All->DriftY = -(gammaYX*PhysicalX + gammaYY*PhysicalY);\n";
      break;

    case DT_NONLINEAR:
        ss <<
            "// diffusivity is constant\n"
            "All->DiffusivityX = diffX;\n"
            "All->DiffusivityY = diffY;\n"
            "\n"
            "Real shiftedX = PhysicalX + origin;\n"
            "Real shiftedMeanX = MeanX + origin;\n"
            "\n"
            "All->DriftX = -(omega * shiftedX + theta*shiftedMeanX);\n"
            "All->DriftY = 0.0;\n";
        break;
    }// which method

    m_computeDriftDiffusivityMethod = ss.str();
}

void DiffusionModel::addEventTime(Real time) {
    m_eventTimes.push_back(time);
}
void DiffusionModel::setEventTimes(const std::vector<Real> &times) {
    m_eventTimes = times;
}
void DiffusionModel::clearEventTimes() {
    m_eventTimes.clear();
}

void DiffusionModel::setEventsScriptContents(const std::string &contents) {
    m_eventsScriptContents = contents;
}
void DiffusionModel::loadEventsScript(const std::string &filename) {
    loadTextFile(m_eventsScriptContents, filename);
}
void DiffusionModel::clearEventsScript() {
    m_eventsScriptContents.clear();
}

void DiffusionModel::setFieldParameter(const std::string &name) {
    m_fieldParameters[name] = std::vector<Real>();

    BOOST_LOG_TRIVIAL(debug) << "  Field \"" << name<<"\"";
}

void DiffusionModel::setParameter(const std::string &name, Real value)
{
    m_parameters[name] = value;
    
    // Set in main namespace/dict, so that it's accessible to everything,
    // including DiffusionModel init and event scripts.
    boost::python::object pyMainModule = boost::python::import("__main__");
    boost::python::object pyMainNamespace = pyMainModule.attr("__dict__");
    pyMainNamespace[name.c_str()] = value;
    BOOST_LOG_TRIVIAL(debug) << "  \"" << name << "\" = " << value;
}

void DiffusionModel::removeParameter(const std::string &name) {
    m_parameters.erase(name);

    // todo: would need to remove from python main namespace?
}

void DiffusionModel::addCompartment(Compartment *compartment)
{
    m_compartments.push_back(compartment);
    c_idToCompartment[compartment->id()] = compartment; // TODO: handle conflicts
}

Compartment * DiffusionModel::addCompartment(const std::string &id) {
    Compartment *c = new Compartment(id);
    addCompartment(c);
    return c;
}

Compartment * DiffusionModel::addCompartment(const std::string &id, int x0, int y0, int x1, int y1) {
    Compartment *c = new Compartment(id, x0, y0, x1, y1);
    addCompartment(c);
    return c;
}

void DiffusionModel::removeCompartment(Compartment *compartment) {
    c_idToCompartment.erase(compartment->id());
    m_compartments.erase(std::remove(m_compartments.begin(), m_compartments.end(), compartment), m_compartments.end());
}

void DiffusionModel::addSpecies(Species *species) {
    m_species.push_back(species);
    c_idToSpecies[species->id()] = species; // TODO: handle conflicts
}

Species * DiffusionModel::addSpecies(const std::string &id, Real diffusionConstant) {
    Species *s = new Species(id, diffusionConstant);
    addSpecies(s);
    return s;
}

void DiffusionModel::removeSpecies(Species *species) {
    c_idToSpecies.erase(species->id());
    m_species.erase(std::remove(m_species.begin(), m_species.end(), species), m_species.end());
}

void DiffusionModel::addReaction(Reaction *reaction) {
    m_reactions.push_back(reaction);
    c_idToReaction[reaction->id()] = reaction; // TODO: handle conflicts
}

Reaction * DiffusionModel::addReaction(const std::string &id) {
    Reaction *r = new Reaction(id);
    addReaction(r);
    return r;
}

// todo: add interface that takes a stoichiometry structure for reactants and products
Reaction * DiffusionModel::addReaction(const std::string &id, Real rate,
                                       Species *reactant1, Species *reactant2,
                                       Species *reactant3, Species *reactant4,
                                       Species *reactant5, Species *reactant6,
                                       Species *reactant7, Species *reactant8)
{
    Reaction *r = new Reaction(id, rate,
                               reactant1, reactant2, reactant3, reactant4,
                               reactant5, reactant6, reactant7, reactant8);
    addReaction(r);
    return r;
}

void DiffusionModel::removeReaction(Reaction *reaction) {
    c_idToReaction.erase(reaction->id());
    m_reactions.erase(std::remove(m_reactions.begin(), m_reactions.end(), reaction), m_reactions.end());
}


void DiffusionModel::solve()
{
    try {
        // print out runtime parameters
        BOOST_LOG_TRIVIAL(info) <<"Runtime configuration: time="<<maxTime()
                               <<", number of Runs="<<numRuns()
                              <<", number of Steps="<<maxSteps()
                             <<", output dt="<<outputInterval();
        createOutputFile(); // create HDF5 file

        if (m_solver & stochastic_mask) {
            BOOST_LOG_TRIVIAL(info) << "Starting stochastic simulation.";

            std::vector<int> stateCache;
            generateInitialStates(stateCache);

            // FUTURE: don't rebuild the kernels for each run???s
            for (m_runIndex = m_runsOffset;
                 m_runIndex < (m_runsOffset+m_numRuns);
                 ++m_runIndex)
            {
                BOOST_LOG_TRIVIAL(info) << "Executing run #" << m_runIndex;
                if (regenerateInitStates())                                
                    generateInitialStates(stateCache);

                time_t startTime = time(0);

                m_runState = stateCache;

                prepareOutputFile();

                BOOST_LOG_TRIVIAL(debug) <<"Solver type is "<<m_solver<<std::flush;

                oclSolve();

                BOOST_LOG_TRIVIAL(info) << "Time used for solver: " << time(0)-startTime << "seconds.";
            }
        }
        else if (m_solver & deterministic_mask) {
            BOOST_LOG_TRIVIAL(info) << "Starting deterministic simulation...";

            time_t startTime = time(0);

            // generate initial states as realS
            generateInitialStates(m_runState, &m_runDeterministicState);

            prepareOutputFile();
            writeAllSpecies(0.f); // write out initial configuration

            oclSolve();

            std::cout << "Time used for solver: " << time(0)-startTime << "seconds.";
        }
    }
    catch (cl::Error error) {
        BOOST_LOG_TRIVIAL(error)
                << "OPENCL ERROR: "
                << error.what()
                << "("
                << oclErrorString(error.err())
                << "). Exiting.";
        // TODO: Do we need to further unwind the stack? probably not..
        exit(1);
    }
    catch (gpgmp::Error& error) {
         BOOST_LOG_TRIVIAL(error) <<"GPGMP ERROR: "
                 << error.what()
                 << ". Exiting.";
        exit(1);
    }
}

void DiffusionModel::generateInitialStates(std::vector<int> &states, std::vector<Real> *realStates)
{
    const size_t dxy = m_gridWidth * m_gridHeight;
    
    // Initialize boost random number generator
    boost::mt19937 seed( m_randomSeed );
    
    // Build States and speciesNumParticles
    states = std::vector<int>(numSpecies() * dxy, 0); // resize and zero-out
    BOOST_FOREACH(Compartment *c, m_compartments)
    {
        // iterate through the species prescription for this compartment
        const std::map<Species *, InitialAmount> &initialAmounts =
                c->initialAmounts();
        
        int speciesIndex = -1;
        BOOST_FOREACH(Species *s, m_species)
        {
            ++speciesIndex; // increment now (we start at -1), as we may "continue" at anytime...
            
            std::map<Species *, InitialAmount>::const_iterator initialAmountIt = initialAmounts.find(s);
            if (initialAmountIt == initialAmounts.end()
                    || initialAmountIt->second.amount == 0) {
                continue;
            }
            const int &initialAmount = initialAmountIt->second.amount;
            const Distribution &dist = initialAmountIt->second.dist;
            
            int *species_init = &states[speciesIndex * dxy];
            
            
            // check if it's there
            assert(speciesIndex!=-1);
            
            // distribute them
            // todo: maybe that should be done in Compartment
            if (dist == gpgmp::HomogeneousDistribution)
            {
                BOOST_LOG_TRIVIAL(debug) << "Distributing " << initialAmount << " particles homogeneously within " << c->id() ;
                
                int area = c->height() * c->width();
                int particlesPerCompartment = initialAmount / area;
                int particlesPerCompartmentPlusOne = particlesPerCompartment + 1;
                int moduloParticlesPerCompartment = initialAmount % area;
                
                int x0, y0, x1, y1;
                c->getBounds(x0, y0, x1, y1);
                
                for (int y=y0; y <= y1; ++y)
                {
                    for (int x=x0; x <= x1; ++x)
                    {
                        if (c->isInCompartment(x, y, 0)) {
                            int offset = y*m_gridWidth + x;
                            species_init[offset] += (offset < moduloParticlesPerCompartment
                                                     ? particlesPerCompartmentPlusOne
                                                     : particlesPerCompartment);
                        }
                    }
                }
            }
            else if (dist == gpgmp::RandomDistribution)
            {
                // distribute them randomly
                BOOST_LOG_TRIVIAL(debug) << "Distributing " << initialAmount << " particles randomly within " << c->id();
                
                for (int i=0; i<initialAmount; i++)
                {
                    int x0, y0, x1, y1;
                    c->getBounds(x0, y0, x1, y1);
                    
                    // we need two int RNGs
                    boost::uniform_int<> gm_int_x( 0, x1-x0);
                    boost::uniform_int<> gm_int_y( 0, x1-x0);
                    boost::variate_generator<boost::mt19937&, boost::uniform_int<> > vgx( seed, gm_int_x);
                    boost::variate_generator<boost::mt19937&, boost::uniform_int<> > vgy( seed, gm_int_y);
                    
                    int x = x0 + vgx();
                    int y = y0 + vgy();
                    while (!c->isInCompartment(x, y, 0)) {
                        x = x0 + vgx();
                        y = y0 + vgy();
                    }
                    size_t offset = y*m_gridWidth + x;
                    species_init[offset] += 1;
                }
            } else if (dist==gpgmp::LinearXDistribution) {
                // distribute according to a gradient parallel to x axis
                 BOOST_LOG_TRIVIAL(debug) << "Distributing " << initialAmount << " particles"
                          << " according to linear x in " << c->id() ;
                
                // get boundaries of compartment
                int x0, y0, x1, y1;
                c->getBounds(x0, y0, x1, y1);
                
                int nx, ny;
                nx = (x1-x0);
                ny = (y1-y0);
                
                // compute gradient
                Real coeff = static_cast<Real>(2.*initialAmount/((nx-1)*nx*ny));
                
                // and fill them in
                for (int y=y0; y <= y1; ++y)
                {
                    for (int x=x0; x <= x1; ++x)
                    {
                        if (c->isInCompartment(x, y, 0)) {
                            size_t offset = y*m_gridWidth + x;
                            species_init[offset] += boost::math::iround(coeff*(x-x0));
                        }
                    }
                }
                
            }
        } // loop through species
        assert(speciesIndex == (int)m_species.size()-1);
    }
    
    // fill in the Real states array if given
    if (realStates != 0) {
        *realStates = std::vector<Real>(numSpecies() * m_gridWidth * m_gridHeight, 0);
        
        for (size_t i=0; i<realStates->size(); i++)
            (*realStates)[i] = static_cast<Real>(states[i]);
    }

    if (!m_initScriptContents.empty())
        callInitScript(states, realStates);

    // set the starting positions in the species list if it allows individual properties
    int speciesIndex = -1;
    BOOST_FOREACH(Species *s, m_species)
    {
        speciesIndex++;

        if (s->hasIndividualProperties()) {
            std::vector<int> ipos;
            int *species_init = &states[speciesIndex * dxy];

            // loop over whole grid and assign positions
            for (size_t i=0; i<dxy; i++) {
                if (species_init[i]>0) {
                    for (int k=0; k<species_init[i]; k++) {
                        ipos.push_back(i);
                    }
                }
            }

            // and set it
            s->setIndividualPositions(ipos);
            BOOST_LOG_TRIVIAL(debug) <<"Found "<<ipos.size()<<" individuals for species "<<s->id()<<"\n";
        }
    } // set individual properties

}


// FUTURE: allow params to be changed? if so, regen calc propensities func with new values?
// todo: for some reason that gets called twice even if only one run is done?
void DiffusionModel::callInitScript(std::vector<int> &states,
                                    std::vector<Real> *realStates)
{
    boost::python::object pyMainModule = boost::python::import("__main__");
    boost::python::object pyMainNamespace = pyMainModule.attr("__dict__");
        
    // inject runtime parameters
    boost::python::dict runtimeDict;
    runtimeDict["nx"] = m_gridWidth;
    runtimeDict["ny"] = m_gridHeight;
    runtimeDict["length"] = m_physicalLength;
    runtimeDict["nspecies"] = m_species.size();
    
    switch (m_solver) {
    case (gpgmp::stochastic_homogeneous):
        runtimeDict["solver"] = "stochastic_homogeneous";
        break;
    case (gpgmp::deterministic_homogeneous_RK4):
        runtimeDict["solver"] = "deterministic_homogeneous";
        break;
    case (gpgmp::deterministic_homogeneous_aqss):
        runtimeDict["solver"] = "deterministic_homogeneous";
        break;
    case (gpgmp::stochastic_inhomogeneous):
        runtimeDict["solver"] = "stochastic_inhomogeneous";
        break;
    default:
        runtimeDict["solver"] = "unknown";
        break;
    }
    pyMainNamespace["runtimeInformation"] = runtimeDict;
    
    // and optional parameters (if any)
    boost::python::dict paramsDict;
    for (std::map<std::string, Real>::const_iterator it = m_parameters.begin();
         it != m_parameters.end();
         ++it) {
        paramsDict[it->first] = it->second;
        BOOST_LOG_TRIVIAL(debug) <<"Injected optional parameter "<<it->first <<"="
                 <<it->second;
    }

    // we also add the field parameters into this dictionary
    std::map<std::string, std::vector<Real> >::iterator itt;
    for (itt = m_fieldParameters.begin();
         itt != m_fieldParameters.end();
         ++itt) {
            // create a buffer for the new parameter
            int fieldSize = states.size();
            std::string id = itt->first;
            m_fieldParameters[itt->first] = std::vector<Real>(fieldSize, 0);

            // make a python array out of it
            intp *N = new intp(1);
            N[0] = fieldSize;

            // todo: is that object ever deleted?
            PyArrayObject *retval_param = (PyArrayObject *)
                    PyArray_SimpleNewFromData(1, N, NPY_FLOAT, &((m_fieldParameters[id])[0]));

            // // convert it to a boost object
            boost::python::handle<> h_param((PyObject*)retval_param);
            boost::python::object my_array_param(h_param);

            // and add that array to the individual dictionary
            paramsDict[id] = my_array_param;
    }

    // and expose it to the python namespace
    pyMainNamespace["parameters"] = paramsDict;

    // if there are any continuous parameter (=individual) species we need to pass a dictionary
    boost::python::dict individualDict;

    // we need to pass the species array as well
    boost::python::dict speciesDict;
    int speciesIndex = -1;
    BOOST_FOREACH(Species *s, m_species)
    {
        ++speciesIndex; // increment now (we start at -1), as we may "continue" at anytime...
        speciesDict[s->id()] = speciesIndex;

        // is it an individual species?
        if (s->hasIndividualProperties()) {
            // add a dictionary for the properties
            boost::python::dict ipDict;
            //ipDict["rx"]=0;
            //ipDict["ry"]=0;

            // create the properties for the species
            // TODO: that might better be done in the Species class?
            std::map<std::string, std::vector<Real> > &properties = s->getIndividualProperties();

            // make a python array out of it
            intp *N = new intp(1);
            N[0] = gpgmp::IndividualSolver::c_maxNumIndividuals;

            // todo: is that object ever deleted?
            PyArrayObject *retval_rx = (PyArrayObject *)
                    PyArray_SimpleNewFromData(1, N, NPY_FLOAT, &((properties["rx"])[0]));

            // // convert it to a boost object
            boost::python::handle<> h_rx((PyObject*)retval_rx);
            boost::python::object my_array_rx(h_rx);

            // and add that array to the individual dictionary
            ipDict["rx"] = my_array_rx;

            // todo: is that object ever deleted?
            PyArrayObject *retval_ry = (PyArrayObject *)
                    PyArray_SimpleNewFromData(1, N, NPY_FLOAT, &((properties["ry"])[0]));

            // // convert it to a boost object
            boost::python::handle<> h_ry((PyObject*)retval_ry);
            boost::python::object my_array_ry(h_ry);

            // and add that array to the individual dictionary
            ipDict["ry"] = my_array_ry;

            // todo: is that object ever deleted?
            PyArrayObject *retval_dx = (PyArrayObject *)
                    PyArray_SimpleNewFromData(1, N, NPY_FLOAT, &((properties["diffx"])[0]));

            // // convert it to a boost object
            boost::python::handle<> h_dx((PyObject*)retval_dx);
            boost::python::object my_array_dx(h_dx);

            // and add that array to the individual dictionary
            ipDict["diffx"] = my_array_dx;

            // todo: is that object ever deleted?
            PyArrayObject *retval_dy = (PyArrayObject *)
                    PyArray_SimpleNewFromData(1, N, NPY_FLOAT, &((properties["diffy"])[0]));

            // // convert it to a boost object
            boost::python::handle<> h_dy((PyObject*)retval_dy);
            boost::python::object my_array_dy(h_dy);

            // and add that array to the individual dictionary
            ipDict["diffy"] = my_array_dy;

            // and add it to the individual dictionary
            individualDict[s->id()] = ipDict;
        }
    }
    pyMainNamespace["species"] = speciesDict;    
    BOOST_LOG_TRIVIAL(debug) << "Injected species array." << std::flush;
    
    // inject the individual property dictionary (if there are any..)
    if (hasContinuousParameterSpecies())
            pyMainNamespace["individualProperties"] = individualDict;

    // the numpy array later (which is a bit easier)
    intp *N = new intp(1);
    N[0] = states.size();
    
    // create a Python array (not numpy)
    // todo: is that object ever deleted?
    PyArrayObject *retval = (PyArrayObject *)
            PyArray_SimpleNewFromData(1, N, NPY_INT, &(states[0]));
    
    // // convert it to a boost object
    boost::python::handle<> h((PyObject*)retval);
    boost::python::object my_array(h);
    
    // // make it visible to main module
    pyMainNamespace["state"] = my_array;
    
    BOOST_LOG_TRIVIAL(debug) << "Injected state array." << std::flush;
    
    // if the deterministic solver is used we provide a pointer to the FPE
    // array too
    if (realStates) {
        
        // create a Python array (not numpy)
        PyArrayObject *retvalReal = (PyArrayObject*)
                PyArray_SimpleNewFromData(1, N, NPY_FLOAT, &((*realStates)[0]));
        
        // // convert it to a boost object
        boost::python::handle<> hReal((PyObject*)retvalReal);
        boost::python::object my_arrayReal(hReal);
        
        // // make it visible to main module
        pyMainNamespace["deterministicState"] = my_arrayReal;
    }
    
    BOOST_LOG_TRIVIAL(debug) << "Calling init script.. " << std::flush;
    
    // and execute the script
    boost::python::object ignored
            = boost::python::exec(m_initScriptContents.c_str(), pyMainNamespace);

	BOOST_LOG_TRIVIAL(debug) << "done. " << std::flush;

    // debug .. check if the species individual property array was set properly
    BOOST_FOREACH(Species *s, m_species)
    {
        if (s->hasIndividualProperties()) {
            BOOST_LOG_TRIVIAL(debug) <<"IP rx for species " << s->id() <<" is "<<s->getIndividualProperties()["rx"][1];
            BOOST_LOG_TRIVIAL(debug) <<"IP ry for species " << s->id() <<" is "<<s->getIndividualProperties()["ry"][29];
            BOOST_LOG_TRIVIAL(debug) <<"IP diffx for species " << s->id() <<" is "<<s->getIndividualProperties()["diffx"][3];
            BOOST_LOG_TRIVIAL(debug) <<"IP diffy for species " << s->id() <<" is "<<s->getIndividualProperties()["diffy"][15];
        }
    }
}

void DiffusionModel::createOutputFile()
{
    if (m_outputFormat != gpgmp::OutputHdf5)
        return;
    
    // HDF5
    size_t numSpecies = m_species.size();
    
    //std::cerr << "Creating HDF5 file " << m_outputFilename << "." << std::endl;
    
    m_hdf5File = H5Fcreate(m_outputFilename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    // write model data
    hid_t modelGroup = H5Gcreate2(m_hdf5File, "/Model", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    hsize_t dims[1];
    dims[0] = 1;
    
    // write number of species attribute
    hid_t status;
    hid_t dataspace = H5Screate_simple(1, dims, 0);
    hid_t dataset = H5Acreate2(modelGroup, "numSpecies", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(dataset, H5T_NATIVE_INT, &numSpecies);
    status = H5Aclose(dataset);
    status = H5Sclose(dataspace);
    
    // write list of species
    dims[0] = numSpecies;
    
    // create list of species
    std::vector<const char *> speciesIds;
    speciesIds.reserve(numSpecies);
    BOOST_FOREACH(Species *species, m_species) {
        speciesIds.push_back(species->id().c_str());
    }
    
    // create HDF5 datatypes for it
    hid_t stringtype = H5Tcopy(H5T_C_S1);
    status = H5Tset_size(stringtype, H5T_VARIABLE);
    
    // create dataspace and write it
    dataspace = H5Screate_simple(1, dims, 0);
    dataset = H5Dcreate2(modelGroup, "Species", stringtype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, stringtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(speciesIds.front()));
    status = H5Dclose(dataset);
    status = H5Sclose(dataspace);
    
    // write number of runs attribute -- TODO: how should be handle runsOffset... is this attribute even needed, as you just check the count inside the group... ???
    dims[0]=1;
    dataspace = H5Screate_simple(1, dims, 0);
    dataset = H5Acreate2(modelGroup, "numRuns", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(dataset, H5T_NATIVE_INT, &(m_numRuns));
    status = H5Aclose(dataset);
    status = H5Sclose(dataspace);
    
    // create compartment group
    hid_t compartmentGroup = H5Gcreate2(modelGroup, "Compartments", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    // iterate through compartment list
    dims[0]=4;
    BOOST_FOREACH(Compartment *c, m_compartments) {
        c->saveHdf5(compartmentGroup);
    }
    
    // close compartment group
    status = H5Gclose(compartmentGroup);
    
    // close model group
    status = H5Gclose(modelGroup);
    
    
    hid_t parentOfRuns = m_hdf5File;
    
    // If present, then the HDF5 runs will be written to /Jobs/{m_outputJobId}/Runs/
    if (!m_outputJobId.empty()) {
        //std::cerr <<"Creating HDF5 Jobs group." << std::endl;
        hid_t jobsRoot = H5Gcreate2(m_hdf5File, "Jobs", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        
        //std::cerr <<"Creating HDF5 job group "<<m_outputJobId<<"." << std::endl;
        
        parentOfRuns = H5Gcreate2(jobsRoot, m_outputJobId.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }
    
    // check if Runs group exists already
    hid_t runs;
    if (H5Lexists(parentOfRuns, "Runs", H5P_DEFAULT) != 1) {
        //std::cerr <<"Runs group does not exist. Creating HDF5 Runs group." << std::endl;
        runs = H5Gcreate2(parentOfRuns, "Runs", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }
}

void DiffusionModel::prepareOutputFile()
{
    // set output data set to 0
    m_runDumpIndex = 0;
    
    if (m_outputFormat != gpgmp::OutputHdf5)
        return;
    
    // HDF5
    herr_t status;
        
    // just open the current one
    m_hdf5File = H5Fopen(m_outputFilename.c_str(),
                         H5F_ACC_RDWR, H5P_DEFAULT);
    
    hid_t parentOfRuns = m_hdf5File;
    
    // If present, then the HDF5 runs will be written to /Jobs/{m_outputJobId}/Runs/
    if (!m_outputJobId.empty()) {
        hid_t jobsRoot = H5Gopen2(m_hdf5File, "Jobs", H5P_DEFAULT);        
        parentOfRuns = H5Gopen2(jobsRoot, m_outputJobId.c_str(), H5P_DEFAULT);
    }
    
    hid_t runs = H5Gopen2(parentOfRuns, "Runs", H5P_DEFAULT);
    
    // create group for current run
    m_hdf5CurrentGroup = H5Gcreate2(runs,
                                    boost::lexical_cast<std::string>(m_runIndex).c_str(),
                                    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);        
    
    // close the data group
    status = H5Gclose(m_hdf5CurrentGroup);
    
    // close the file
    status = H5Fclose(m_hdf5File);
}

void DiffusionModel::openOutputFile()
{    
    m_hdf5File = H5Fopen(m_outputFilename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    
    hid_t parentOfRuns = m_hdf5File;
    
    // If present, then the HDF5 runs will be written to /Jobs/{m_outputJobId}/Runs/
    if (!m_outputJobId.empty()) {
        std::cerr <<"Opening HDF5 Jobs group." << std::endl;
        hid_t jobsRoot = H5Gopen2(m_hdf5File, "Jobs", H5P_DEFAULT);
        
        
        parentOfRuns = H5Gopen2(jobsRoot, m_outputJobId.c_str(), H5P_DEFAULT);
    }
    
    hid_t runs = H5Gopen2(parentOfRuns, "Runs", H5P_DEFAULT);
    
    // create group for current run
    m_hdf5CurrentGroup = H5Gopen2(runs, boost::lexical_cast<std::string>(m_runIndex).c_str(), H5P_DEFAULT);
}

void DiffusionModel::closeOutputFile() {
    
    if (m_outputFormat != gpgmp::OutputHdf5)
        return;
    
    herr_t status;
    
    // close the data group
    status = H5Gclose(m_hdf5CurrentGroup);
    
    // close the file
    status = H5Fclose(m_hdf5File);
}

void DiffusionModel::writeAllSpecies(Real time, IndividualSolver *isolver)
{
    // Gather buffer information, so that we can access it in a generic way
    const char *buffer = 0;
    size_t      bufferStep = 0;
    hid_t       bufferType = 0;
    if (m_solver & stochastic_mask) {
        buffer = (const char *)&m_runState.front();
        bufferStep = sizeof(int) * m_gridWidth * m_gridHeight;
        bufferType = H5T_NATIVE_INT;
    }
    else {
        buffer = (const char *)&m_runDeterministicState.front();
        bufferStep = sizeof(Real) * m_gridWidth * m_gridHeight;
        bufferType = HDF_REAL;
    }
    
    if (m_outputFormat == gpgmp::OutputHdf5) {
        BOOST_LOG_TRIVIAL(debug) << "Writing HDF5 output file at t="<<time<<".";

        // open output file
        openOutputFile();
        
        // status variable
        herr_t status;
        
        // dimensions of the dataset
        hsize_t dims[2];
        
        // data set and data space
        hid_t dataspace, dataset;
        
        // create the data group
        std::string dataGroupName = "Dump_" + boost::lexical_cast<std::string>(time);
        hid_t currentDataGroup = H5Gcreate2(m_hdf5CurrentGroup,
                                            dataGroupName.c_str(),
                                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        
        // write time attribute
        dims[0] = 1;
        dataspace = H5Screate_simple(1, dims, 0);
        dataset = H5Acreate2(currentDataGroup, "time", HDF_REAL, dataspace, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Awrite(dataset, HDF_REAL, &time);
        status = H5Aclose(dataset);
        status = H5Sclose(dataspace);
        

        // and write the species
        writeSingleTimeHDF5(currentDataGroup, buffer, bufferStep, bufferType);

        // write fields
        if (getFieldParameters().size()>0) {
            // create data group for fields
            std::string fieldsGroupName = "Fields";
            hid_t fieldsDataGroup = H5Gcreate2(currentDataGroup,
                                                fieldsGroupName.c_str(),
                                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

            // set dimensions for each dataset (they're all the same .. grid dims)
            dims[0] = m_gridWidth;
            dims[1] = m_gridHeight;

            // go through all parameters
            for (std::map<std::string, std::vector<Real> >::const_iterator itt = getFieldParameters().begin();
                 itt != getFieldParameters().end();
                 ++itt) {

                // create data group name
                hid_t dataspaceField, datasetField;

                dataspaceField = H5Screate_simple(2, dims, 0);
                datasetField = H5Dcreate2(fieldsDataGroup, (itt->first).c_str(), HDF_REAL, dataspaceField,
                                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                status = H5Dwrite(datasetField, HDF_REAL, H5S_ALL, H5S_ALL, H5P_DEFAULT,  (const char *)&(itt->second).front());
                status = H5Dclose(datasetField);
                status = H5Sclose(dataspaceField);
            }

            // close the group for the individual list
            status = H5Gclose(fieldsDataGroup);

        }

        // create data group for individuals list
        std::string individualGroupName = "Individuals";
        hid_t individualDataGroup = H5Gcreate2(currentDataGroup,
                                            individualGroupName.c_str(),
                                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        // go through species again and write indidivual lists
        for (size_t i=0; i<m_species.size(); i++)
        {
            Species *species = m_species[i];
            if (species->hasIndividualProperties()) {
                if (isolver == 0) {
                    BOOST_LOG_TRIVIAL(error) <<"ERROR: WriteAllSpecies needs pointer to individual solver but has only NULL. Exiting.";
                    exit(1);
                }

                // get position list
                std::vector<cl_ulong> pos = isolver->getPositionList(i);

                // only store list if there are any individuals
                if (pos.size()>0)
                {
                    dims[0] = pos.size();
                    // write individuals key pairs
                    dataspace = H5Screate_simple(1, dims, 0);
                    dataset = H5Dcreate2(individualDataGroup, species->id().c_str(), H5T_NATIVE_ULONG, dataspace,
                                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                    const char *ibuffer = (const char *)&pos.front();
                    status = H5Dwrite(dataset, H5T_NATIVE_ULONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, ibuffer);
                    status = H5Dclose(dataset);
                    status = H5Sclose(dataspace);

                    // get properties
                    std::map<int, std::vector<Real> > props = isolver->getPropertiesList(i);

                    // write properties
                    std::string dsname;
                    dsname = species->id()+".DiffusivityX";
                    dataspace = H5Screate_simple(1, dims, 0);
                    dataset = H5Dcreate2(individualDataGroup, dsname.c_str(), HDF_REAL, dataspace,
                                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                    ibuffer = (const char *)&(props[0]).front();
                    status = H5Dwrite(dataset, HDF_REAL, H5S_ALL, H5S_ALL, H5P_DEFAULT, ibuffer);
                    status = H5Dclose(dataset);
                    status = H5Sclose(dataspace);
                    dsname = species->id()+".DiffusivityY";
                    dataspace = H5Screate_simple(1, dims, 0);
                    dataset = H5Dcreate2(individualDataGroup, dsname.c_str(), HDF_REAL, dataspace,
                                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                    ibuffer = (const char *)&(props[1]).front();
                    status = H5Dwrite(dataset, HDF_REAL, H5S_ALL, H5S_ALL, H5P_DEFAULT, ibuffer);
                    status = H5Dclose(dataset);
                    status = H5Sclose(dataspace);
                    dsname = species->id()+".DriftX";
                    dataspace = H5Screate_simple(1, dims, 0);
                    dataset = H5Dcreate2(individualDataGroup, dsname.c_str(), HDF_REAL, dataspace,
                                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                    ibuffer = (const char *)&(props[2]).front();
                    status = H5Dwrite(dataset, HDF_REAL, H5S_ALL, H5S_ALL, H5P_DEFAULT, ibuffer);
                    status = H5Dclose(dataset);
                    status = H5Sclose(dataspace);
                    dsname = species->id()+".DriftY";
                    dataspace = H5Screate_simple(1, dims, 0);
                    dataset = H5Dcreate2(individualDataGroup, dsname.c_str(), HDF_REAL, dataspace,
                                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                    ibuffer = (const char *)&(props[3]).front();
                    status = H5Dwrite(dataset, HDF_REAL, H5S_ALL, H5S_ALL, H5P_DEFAULT, ibuffer);
                    status = H5Dclose(dataset);
                    status = H5Sclose(dataspace);
                }
            }
        }

        // close the group for the individual list
        status = H5Gclose(individualDataGroup);

        // close the model group
        status = H5Gclose(currentDataGroup);
        
        // close output file
        closeOutputFile();
    }
    else {
        BOOST_LOG_TRIVIAL(error) <<"Error: Unknown output format !";
        exit(1);
    } // OUTPUT_FORMAT
    
    // increase output data set number
    m_runDumpIndex++;
}

void DiffusionModel::writeSingleTimeHDF5(hid_t currentDataGroup, const char *buffer, size_t  bufferStep, hid_t bufferType)
{
    // dimensions of the dataset
    hsize_t dims[2];
    dims[0] = m_gridWidth;
    dims[1] = m_gridHeight;

    // data set and data space
    hid_t dataspace, dataset;

    // status variable
    herr_t status;

    BOOST_FOREACH(Species *species, m_species)
    {
        // create data group name
        dataspace = H5Screate_simple(2, dims, 0);
        dataset = H5Dcreate2(currentDataGroup, species->id().c_str(), bufferType, dataspace,
                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, bufferType, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);

        buffer += bufferStep;
    }

}

// FUTURE: allow params to be changed? if so, regen calc propensities func with new values?
void DiffusionModel::eventCallback(Real time) {
    
    boost::python::object pyMainModule = boost::python::import("__main__");
    boost::python::object pyMainNamespace = pyMainModule.attr("__dict__");
    
    // injects the runtime info into scope of the main module -- TODO: use local instead?
    pyMainNamespace["dx_"] = m_gridWidth;
    pyMainNamespace["dy_"] = m_gridHeight;
    pyMainNamespace["nSpecies_"] = m_species.size();
    pyMainNamespace["time_"] = time;
    
    // we need to pass the species array as well
    boost::python::dict speciesDict;
    int speciesIndex = -1;
    BOOST_FOREACH(Species *s, m_species)
    {
        ++speciesIndex; // increment now (we start at -1), as we may "continue" at anytime...
        speciesDict[s->id()] = speciesIndex;
    }
    pyMainNamespace["species"] = speciesDict;
    
    // the numpy array later (which is a bit easier)
    intp *N = new intp(1);
    N[0] = numSpecies() * m_gridWidth * m_gridHeight;
    
    // create a Python array (not numpy)
    PyArrayObject *retval = (PyArrayObject*)
            PyArray_SimpleNewFromData(1, N, NPY_INT, &(m_runState[0]));
    
    // convert it to a boost object
    boost::python::handle<> h((PyObject*)retval);
    boost::python::object my_array(h);
    
    // make it visible to main module
    pyMainNamespace["state"] = my_array;
    
    // and execute the script
    boost::python::object ignored
            = boost::python::exec(m_eventsScriptContents.c_str(), pyMainNamespace);
}


Real DiffusionModel::subvolumeSize() const {
    return static_cast<Real>(pow(m_physicalLength/m_gridWidth*1e-6f, 3)*1e3f);
}


int DiffusionModel::particleNumberFromConcentration(Real concentration) const {
    return boost::math::iround(concentration*Avogadro*subvolumeSize());
}


Real DiffusionModel::zerothOrderReactionRate(Real rate) const {
    return rate * subvolumeSize() * Avogadro;
}

Real DiffusionModel::secondOrderReactionRate(Real rate) const {
    return rate / (Avogadro * subvolumeSize());
}


void DiffusionModel::loadTextFile(std::string &contentsOut,
                                  const std::string &filename) const
{
    contentsOut.clear();
    
    // Open the file
    std::ifstream infile;
    infile.open(filename.c_str(), std::ifstream::in);
    if (!infile.good()) {
        std::cerr << "Failed to open \"" << filename << "\" for reading from." << std::endl;
        exit(1);
    }
    
    // Get the length of the file
    infile.seekg(0, std::ios::end);      // go to the end
    contentsOut.reserve(infile.tellg()); // reserve size = position at end
    infile.seekg(0, std::ios::beg);      // go back to the beginning
    
    // Copy the file contents
    char buffer[1024];
    while (!infile.eof()) {
        infile.read(buffer, 1024);
        contentsOut.append(buffer, infile.gcount());
    }
}

cl_mem DiffusionModel::oclCreateReactionMask() const
{
    if (numReactions() <= 0)
        return 0;
    
    // prepare the mask array on host
    std::vector<float> h_reactionMask(m_gridWidth * m_gridHeight * numReactions());
    float *maskIt = &h_reactionMask.front();
    for (size_t r=0; r<numReactions(); r++) {
        for (size_t y=0; y<m_gridHeight; y++) {
            for (size_t x=0; x<m_gridWidth; x++) {
                *(maskIt++) = m_reactions.at(r)->reactionMask(x, y, (size_t)0);
            }
        }
    }

    // copy mask to OpenCL texture
    cl_mem_flags flags = CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR;
    cl::ImageFormat format(CL_R, CL_FLOAT);
    cl::Image *image = (numReactions() == 1
                        ? (cl::Image *)new cl::Image2D(*m_context, flags, format, m_gridWidth, m_gridHeight, 0, &h_reactionMask[0])
                        : (cl::Image *)new cl::Image3D(*m_context, flags, format, m_gridWidth, m_gridHeight, numReactions(), 0, 0, &h_reactionMask[0])
                          );
    cl_mem mem = (cl_mem)(*image)();
    clRetainMemObject(mem); // retain, as cl::Image's dtor will call to release it
    return mem; // return opencl's nice dimension independant pointer (you can't setArg with cl::Image)
}

void DiffusionModel::oclSolve()
{
    // ------------------------------------------------------------------------
    // Initialization
    
    // all possible event types
    enum EventType { DiffusionEvent, OutputEvent, AlertEvent } nextEvent;

    // create solver
    Solver *solver;

    if (m_solver == stochastic_homogeneous
            || m_solver == deterministic_homogeneous_RK4
            || m_solver == deterministic_homogeneous_aqss)
        solver = new HomogeneousSolver(*this, &m_runState, &m_runDeterministicState,
                                       m_solver, m_solver & deterministic_mask, hasContinuousParameterSpecies());
    else if (m_solver == stochastic_inhomogeneous
             || m_solver == stochastic_inhomogeneous_fpe
             || m_solver == deterministic_inhomogeneous_RK4
             || m_solver == deterministic_inhomogeneous_aqss)
        solver = new InhomogeneousSolver(*this, &m_runState, &m_runDeterministicState,
                                         m_solver,
                                         m_solver & deterministic_mask, hasContinuousParameterSpecies(), computeMoments(), nonlinearDiffusivity() );
    else  {
         BOOST_LOG_TRIVIAL(error) <<"Error in oclSolve: Unknown solver type "<<m_solver<<"\n";
        exit(1);
    }

    writeAllSpecies(0.f, solver->getIndividualSolver()); // write out initial configuration
    std::flush(std::cout);

    // simulation time
    Real simTime=0;
            
    // initialize the step counter for each species to 1.0
    //std::vector<Real> stepCount(numSpecies(), 1.0);
    bool updateFields = (m_fieldParameters.size()>0 && returnSetFieldMethod().size()>0) ? true : false;

    // get initial time steps from solver
    std::vector<Real> dt = solver->getNextDiffusionTime();
    
    // Initialize the data output timer
    Real outputLast = 0.f;
    
    // Initialize alert times iterator
    std::vector<Real>::iterator eventTimesIterator = m_eventTimes.begin();
    
    // use more accurate timer
    timeval tim;
    gettimeofday(&tim, NULL);
    double t1=tim.tv_sec+(tim.tv_usec/1000000.0);

    // take this as a random seed to have different results for different runs
    boost::mt19937 rng( m_randomSeed++);
    boost::uniform_int<> dist(0, INT_MAX);

    BOOST_LOG_TRIVIAL(debug) <<"Inhomogeneous solver seeds with "<<m_randomSeed;

    // main loop
    std::flush(std::cout);

    int n=0;
    //int ndiff = 1;

    // todo: we should have some sort of progress information..
    while ((m_maxSteps < 0 || n<m_maxSteps) // allow any -ve number to mean infinity, not just the InfSteps "sugar"
           && simTime<m_maxTime)
    {
        int seed = dist(rng)+n;
        // increase step count
        n++;
        
        // minimum species (only for stochastic solver)
        size_t minspec;
        
        // compute next dt and determine the event type
        // diffusion time step
        Real mindt;
        //mindt=stepCount[0]*dt[0];
        mindt = dt[0];
        minspec = 0;
        for (size_t i=1; i<numSpecies(); i++) {
            if (dt[i] < mindt) {
                mindt=dt[i];
                minspec = i;
            }
        }
        
        if ((n % 1000) == 0)
            BOOST_LOG_TRIVIAL(trace) <<"Executing step "<<n<< "at t="<<mindt<<"."<<std::flush;

        // by default the next event is diffusion
        nextEvent = DiffusionEvent;
        
        // alert time step
        if (eventTimesIterator != m_eventTimes.end()) {
            Real mindtAlert = *eventTimesIterator;
            
            // set it as mindt if smaller
            if (mindtAlert < mindt) {
                mindt = mindtAlert;
                nextEvent = AlertEvent;
            }
        }
        
        // output time step
        Real mindtOutput = outputLast + m_outputInterval;
        
        // set it as mindt if smaller
        if (mindtOutput < mindt) {
            mindt = mindtOutput;
            nextEvent = OutputEvent;
        }
        
        // first, we perform the reaction step
        Real deltat = mindt-simTime;

        if (numReactions() > 0)
            solver->doReactions(deltat, seed);

        // now we need to update the fields if necessary
        if (updateFields && deltat>0)
            solver->updateFields(deltat, simTime);

        // If we have a diffusion event, we need to diffuse each species
        // whose time has come and update the step counter accordingly
        if (nextEvent == DiffusionEvent)
        {
            for (size_t i=0; i<numSpecies(); i++)
            {
                // To avoid float comparison we diffuse the minimum species
                // AND all species which would be smaller
                if (i==minspec || dt[i] <= mindt)
                    if (dt[i]>0) // diffuse if necessary
                        solver->diffuse(i, mindt, seed);
            } // diffuse species

            // update dt
            dt = solver->getNextDiffusionTime();
        } // diffusion event
        else if (nextEvent == OutputEvent)
        {
            // if output time is reached, write data
            if (simTime >= outputLast+m_outputInterval)
            {

                // update host state array
                solver->synchronizeStateBufferFromDevice();

                // update field array (if any)
                if (m_fieldParameters.size()>0)
                    solver->synchronizeFieldBufferFromDevice();

                // and write it all (including fields and individuals)
                writeAllSpecies(mindt, solver->getIndividualSolver());
                
                // set next output interval
                outputLast+=m_outputInterval;
            }
        }
        else if (nextEvent == AlertEvent)
        {
            // perform alert event
             BOOST_LOG_TRIVIAL(info) << "Alert event reached at t="<<mindt<<".\n";

            ++eventTimesIterator;
            
            // update host state array
            solver->synchronizeStateBufferFromDevice();

            // update field array (if any)
            if (m_fieldParameters.size()>0)
                solver->synchronizeFieldBufferFromDevice();

            // execute event
            this->eventCallback(mindt);

            // and copy buffer back to device
            solver->synchronizeStateBufferToDevice();

            // update field array (if any)
            if (m_fieldParameters.size()>0)
                solver->synchronizeFieldBufferToDevice();

        } // events
        
        // set global clock
        simTime=mindt;
        
    } // end main loop
    

    // write timing
    gettimeofday(&tim, NULL);
    double t2=tim.tv_sec+(tim.tv_usec/1000000.0);
    BOOST_LOG_TRIVIAL(info) << "Algorithm took" << t2-t1 <<"seconds.";

    // get result
    solver->synchronizeStateBufferFromDevice();
    solver->synchronizeFieldBufferFromDevice();

    // write final state
    writeAllSpecies(simTime, solver->getIndividualSolver());

    // clean up
    delete solver;
}

bool DiffusionModel::hasLocalizedReactions() const {
    bool localized=false;

    // go through all reactions
    BOOST_FOREACH(Reaction *r, m_reactions) {
        if (r->isLocalized()) {
            localized=true;
            break;
        }
    }

    return localized;
}

} // namespace gpgmp


// TODO: integrate or remove
void printRange(std::ostream &stream, const cl::NDRange &range)
{
    stream << "[";
    for (size_t i=0; i < range.dimensions(); ++i ) {
        stream << (i>0 ? ", " : "") << range[i];
    }
    stream << "]";
}

const char* oclErrorString(cl_int error)
{
    static const char* errorString[] = {
        "CL_SUCCESS",
        "CL_DEVICE_NOT_FOUND",
        "CL_DEVICE_NOT_AVAILABLE",
        "CL_COMPILER_NOT_AVAILABLE",
        "CL_MEM_OBJECT_ALLOCATION_FAILURE",
        "CL_OUT_OF_RESOURCES",
        "CL_OUT_OF_HOST_MEMORY",
        "CL_PROFILING_INFO_NOT_AVAILABLE",
        "CL_MEM_COPY_OVERLAP",
        "CL_IMAGE_FORMAT_MISMATCH",
        "CL_IMAGE_FORMAT_NOT_SUPPORTED",
        "CL_BUILD_PROGRAM_FAILURE",
        "CL_MAP_FAILURE",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "CL_INVALID_VALUE",
        "CL_INVALID_DEVICE_TYPE",
        "CL_INVALID_PLATFORM",
        "CL_INVALID_DEVICE",
        "CL_INVALID_CONTEXT",
        "CL_INVALID_QUEUE_PROPERTIES",
        "CL_INVALID_COMMAND_QUEUE",
        "CL_INVALID_HOST_PTR",
        "CL_INVALID_MEM_OBJECT",
        "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR",
        "CL_INVALID_IMAGE_SIZE",
        "CL_INVALID_SAMPLER",
        "CL_INVALID_BINARY",
        "CL_INVALID_BUILD_OPTIONS",
        "CL_INVALID_PROGRAM",
        "CL_INVALID_PROGRAM_EXECUTABLE",
        "CL_INVALID_KERNEL_NAME",
        "CL_INVALID_KERNEL_DEFINITION",
        "CL_INVALID_KERNEL",
        "CL_INVALID_ARG_INDEX",
        "CL_INVALID_ARG_VALUE",
        "CL_INVALID_ARG_SIZE",
        "CL_INVALID_KERNEL_ARGS",
        "CL_INVALID_WORK_DIMENSION",
        "CL_INVALID_WORK_GROUP_SIZE",
        "CL_INVALID_WORK_ITEM_SIZE",
        "CL_INVALID_GLOBAL_OFFSET",
        "CL_INVALID_EVENT_WAIT_LIST",
        "CL_INVALID_EVENT",
        "CL_INVALID_OPERATION",
        "CL_INVALID_GL_OBJECT",
        "CL_INVALID_BUFFER_SIZE",
        "CL_INVALID_MIP_LEVEL",
        "CL_INVALID_GLOBAL_WORK_SIZE",
    };

    const int errorCount = sizeof(errorString) / sizeof(errorString[0]);

    const int index = -error;

    return (index >= 0 && index < errorCount) ? errorString[index] : "";

}
