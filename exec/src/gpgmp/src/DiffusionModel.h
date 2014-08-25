/*
 * DiffusionModel.h
 *
 *  Created on: 09/12/2009
 *  Created by: matthias
 */

#ifndef __gpgmp_DiffusionModel_h__
#define __gpgmp_DiffusionModel_h__

#include "Base.h"
#include "Reaction.h"
#include "Globals.h"

#include <hdf5.h>

#include <boost/container/stable_vector.hpp>

#include <list>
#include <iomanip>
#include <map>
#include <vector>

#include <Python.h>

namespace gpgmp {


/**
 * Enumeration for output format. Currently only HDF5 is supported.
 *
 * \todo Implement NETCDF
 * \todo use Hdf5 C++ wrapper instead
 */
enum OutputFormat {
    //OutputAscii,
    OutputHdf5 ///< HDF5 output
//    OutputNetCdf
};

/**
 * Enumeration for the solver type.
 *
 */
// Note: keep setSolver(const std::string &solverName) in sync with this
enum SolverType {
    stochastic_homogeneous           = 1,  //!< stochastic_homogeneous
    stochastic_inhomogeneous         = 2,  //!< stochastic_inhomogeneous
    stochastic_inhomogeneous_fpe     = 4,  //!< stochastic inhomogeneous solver using FPE discretization
    deterministic_homogeneous_RK4    = 8,  //!< deterministic_homogeneous_RK4
    deterministic_homogeneous_aqss   = 16, //!< deterministic alpha QSS
    deterministic_inhomogeneous_RK4  = 32,
    deterministic_inhomogeneous_aqss = 64,
    stochastic_mask           = stochastic_homogeneous | stochastic_inhomogeneous | stochastic_inhomogeneous_fpe,
    deterministic_mask        = deterministic_homogeneous_RK4 | deterministic_homogeneous_aqss | deterministic_inhomogeneous_RK4 | deterministic_inhomogeneous_aqss
};
    
    
class Compartment;
class Kernel;
class Species;
class IndividualSolver;

/**
 * Basic exception class.
 *
 */
class Error : public std::exception
{
private:
    int err_; ///< Error number.
    const char * errStr_; ///< Error meaning.
public:

    /**
     * Constructs the error.
     *
     * @param err Error number
     * @param errStr Meaning
     */
    Error(int err, const char * errStr = NULL) : err_(err), errStr_(errStr)
    {}

    ~Error() throw() {}

    /**
     * Returns human-readable form of the error.
     *
     * @return Error message
     */
    virtual const char * what() const throw ()
    {
        if (errStr_ == NULL) {
            return "empty";
        }
        else {
            return errStr_;
        }
    }

    /**
     * Getter method for the err_ parameter.
     */
    cl_int err(void) const { return err_; }
};

/**
 * The main model class. This class fully defines the model, including simulation parameters, compartments, species and reactions as well as
 * the various methods (such as drift-diffusivity method, initialization) and the events. In brief, any feature that is accessible through the
 * Inchman GUI is represented in this class.
 */
class DiffusionModel : public Base {

public:
    // static variables
    static const Real Avogadro; ///< Avogadro number
    static const int InfSteps = -1; ///< If the number of steps is set to this value, the solver assumes the maximum number of steps is infinite.
    
    /**
     * Constructs an empty model.
     *
     * @param id The model name
     * @param managePythonInstance If true (the default), a Python instance will be
     * created straight away and will be destroyed / finialized on destruction.
     * The numpy arrays Python will be imported regardless of the paramater's value
     * @param randomSeed Can be used to initialize the RNG with a different seed to make experiments reproducable. Not fully implemented yet.
     */
    DiffusionModel(const std::string &id = "Model", bool managePythonInstance = true, unsigned int randomSeed = 0);

	/*
	 * Compability constructor. THis is decrapated.
	 */
	DiffusionModel(Real length, int dx, int dy);

    /**
     * Destructor.
     */
    virtual ~DiffusionModel();
    

    /** @name Model options.
     *  Specifies the main model parameters.
     */
    ///@{
    /**
     * Sets the physical width of the grid. Currently only square simulation domains are supported.
     * @param length The width.
     */
    void setPhysicalLength(Real length);

    /**
     * Sets the number of subvolumes in x-direction. Currently only square domains are supported and the number needs to be a multiple of 8.
     * @param width Number of subvolumes in x-direction.
     */
    void setGridModelWidth(size_t width);
    /**
     * Sets the number of subvolumes in y-direction. Currently only square domains are supported and the number needs to be a multiple of 8.
     * @param height Number of subvolumes in y-direction.
     */
    void setGridModelHeight(size_t height);
    /**
     * Sets the number of subvolumes in both directions. Currently only square domains are supported and the number needs to be a multiple of 8.
     * If the width and height are different, the behaviour of the solver is unspecified (most likely it will fail or yield incorrect results).
     * @param width Number of subvolumes in x-direction.
     * @param height Number of subvolumes in y-direction.
     */
    void setGridDims(size_t width, size_t height);
    /**
     * Returns the physical width of the domain.
     * @return Physical width of the domain.
     */
    inline Real physicalLength() const;
    /**
     * Returns the grid area, i.e. total number of subvolumes.
     * @return Total number of subvolumes.
     */
    inline size_t gridArea() const;
    /**
     * Returns the grid width.
     * @return grid width
     */
    inline size_t gridWidth() const;
    /**
     * Returns the grid height.
     * @return grid height
     */
    inline size_t gridHeight() const;
    /**
     * Sets the boundary and source mask for the integration domain boundaries.
     * The boundary mask determines what happens to particles trying to leave
     * the integration domain. Currently, there are two possibilities:
     * - \c reflective The particle bounces back from the boundary.
     * - \c outflow The particle can freely leave the boundary.
     *
     * \c boundaryMask is a 4-vector with components 0,1,2,3 referring to
     * left (x=0), right (x=dx-1), lower (y=0), and upper (y=dy-1) boundary.
     * Possible values are 0 (reflective) and 1 (outflow).
     *
     * In addition, particle sources at the boundary can be implemented with
     * the source mask (see SemiInfiniteSlab for an example). At the according
     * boundary \c i, the particle number is fixed at a value \c sourceMask[i].
     * E.g., \c sourceMask[0]=10 means that \f$N(x=-1)=10\f$. Note that both,
     * \c boundaryMask and \c sourceMask, apply for \e all species.
     * @param boundaryMask The model boundaries
     * @param sourceMask The source boundaries for inflow conditions
     */
    void setBoundaryMasks(const cl_int4 &boundaryMask, const cl_int4 &sourceMask);
    /**
     * Sets all boundaries to reflecting.
     */
    void clearBoundaryMasks();
    /**
     * @return The boundary conditions of the model.
     */
    inline const cl_int4 &boundaryMask() const;
    /**
     * @return The inflow boundary conditions of the model.
     */
    inline const cl_int4 &sourceMask() const;
    ///@}

    /** @name Solver parameters.
     *  Type and basic solver parameters
     */
    ///@{

    /**
     * Sets the solver type by the enumaration.
     * @param solver The solver type.
     */
    void setSolver(SolverType solver);

    /**
     * Sets the solver type by name.
     *
     * @param solverName The solver type.
     */
    void setSolver(const std::string &solverName);

    /**
     * Gets the current solver type.
     * @return The solver type.
     */
    inline SolverType solver() const;

    /**
     * Sets the probability for a particle to stay in its cell for the stochastic homogeneous solver.
     * If it is 0, all particles
     * will leave the cell at every time step.
     * @param p0
     */
    void setP0(Real p0);
    /**
     * Gets the move probability in the stochastic homogeneous solver.
     * @return Move probability
     */
    inline Real p0() const;

    /**
     * Sets the number of independent runs that are executed for this experiment.
     * @param numRuns Number of runs
     */
    void setNumRuns(unsigned numRuns);
    /**
     * Gets the number of runs.
     */
    inline unsigned numRuns() const;
    /**
     * Sets the maximum number of steps before the simulation stops. If set to -1,
     * the simulation will run until the maximum time is reached.
     *
     * @param maxSteps Maximum step number
     */
    void setMaxSteps(int maxSteps);
    /**
     * Returns the maximum number of steps
     * @return Maximum step number
     */
    inline int maxSteps() const;
    /**
     * Sets the maximum simulation time.
     *
     * @param maxTime Maximum simulation time
     */
    void setMaxTime(Real maxTime);

    /**
     * Returns the maximum simulation time.
     * @return The maximum simulation time.
     */
    inline Real maxTime() const;

    /**
     * Sets the interval (in simulation time) for output dumps.
     *
     * @param interval Output dump interval
     */
    void setOutputInterval(Real interval);

    /**
     * Gets the dump interval for simulation output.
     * @return Dump interval
     */
    inline Real outputInterval() const;

    /**
     * Sets the output format. Currently only HDF5 is fully supported.
     * @param format Output format
     */
    void setOutputFormat(OutputFormat format);
    /**
     * Returns the current output format.
     * @return Current output format
     */
    inline OutputFormat outputFormat() const;
    /**
     * Sets the output file name.
     * @param filename The file name
     */
    void setOutputFilename(const std::string &filename);
    /**
     * Gets the name of the output file
     * @return File name
     */
    inline const std::string &outputFilename() const;
    ///@}


    /** @name Model specifications.
     *  Specifies the model.
     */
    ///@{
    //inline const std::string &initScriptContents() const; // todo: This method is not implemented anywhere?
    /**
     * Sets the contents of the init script that can be used to specify the initial conditions of the model.
     * @param script The script content
     */
    void setInitScriptContents(const std::string &script);
    /**
     * Loads an init script from a file.
     * @param filename File name of the init script.
     */
    void loadInitScript(const std::string &filename);
    /**
     * Deletes the init script.
     */
    void clearInitScript();
    /**
     * If set to true, calls the initialization script for each new run.
     */
    void setRegenerateInitStates(bool regenerateInitStates) {m_regenerateInitStates = regenerateInitStates;}
    /**
     * @return If set to true, calls the initialization script for each new run.
     */
    bool regenerateInitStates() const {return m_regenerateInitStates;}

    // compute diffusion/drift fields for inhomogeneous solver
    /**
     * Sets the method that the inhomgeoenous solvers use to compute the diffusivity and drift. This is
     * basically C code that will be inserted into the Open-CL templates.
     *
     * @param script The script
     */
    inline void setComputeDriftDiffusivityMethod(const std::string &script);
    /**
     * Sets one of the standard methods to compute drift/diffusivity.
     *
     * @param method The standard method.
     */
    void setComputeDriftDiffusivityMethod(const DiffusionType method);
    /**
     * Gets the current method to compute drift and diffusivity.
     * @return The drift/diffusivity method.
     */
    inline const std::string &computeDriftDiffusivityMethod() const;
    /**
     * Loads a method to compute the drift/diffusivity from file.
     * @param filename File name.
     */
    void loadComputeDriftDiffusivityMethod(const std::string &filename);
    /**
     * Clears the drift/diffusivity method.
     */
    void clearComputeDriftDiffusivityMethod();
    /**
     * If set to true, the drift/diffusivity method is called after each
     * diffusion sweep to recompute the drift/diffusivity. This is necessary,
     * for example, if the drift/diffusivity depends on the species count.
     *
     * @param nonlinear If true, re-compute drift/diffusivity after each time step.
     */
    inline void setNonlinearDiffusivity(bool nonlinear);
    /**
     * @return True, if the drift/diffusivity is recomputed after each time step.
     */
    inline bool nonlinearDiffusivity() const;
    /**
     * If set to true, the mean x and y position of each species are computed
     * after each time step and can be used to compute the drift/diffusivity.
     *
     * @param computeMoments If true, computes the moments after each time step.
     */
    inline void setComputeMoments(bool computeMoments);
    /**
     * @return True, iof the moments are recomputed after each time step.
     */
    inline bool computeMoments() const;

    // set field method
    /**
     * Sets the method to update field variables after each time step.
     * @param script The update method
     */
    inline void setFieldMethod(const std::string &script);
    /**
     * @return The update method for the field variables.
     */
    inline const std::string &returnSetFieldMethod() const;

    // new individuals method
    /**
     * Sets the method that initializes newly create individuals.
     * @param method The intialization method.
     */
    inline void setNewIndividualsMethod(const std::string &method);
    /**
     * @return The method to initialize new individuals.
     */
    inline const std::string &getNewIndividualsMethod() const;

    // events
    /**
     * sets a script that is called each time the simulation reaches an
     * event (currently events can only be triggered by simulation time)
     *
     * @param script The event script.
     */
    void setEventsScriptContents(const std::string &script);
    /**
     * Loads an event script from file.
     *
     * @param filename File name for the event script.
     */
    void loadEventsScript(const std::string &filename);
    /**
     * Clears the event script.
     */
    void clearEventsScript();
    /**
     * @return Returns the event script.
     */
    inline const std::string &eventsScriptContents() const;
    /**
     * Adds an event time to the event list. Events allow the user to directly
     * modify
     * the internal state array at pre-defined times.
     *
     * @param time Event time.
     *
     * \sa setEventsScriptContents
     */
    void addEventTime(Real time);
    /**
     * Sets a list of event times.
     * @param times The event times.
     */
    void setEventTimes(const std::vector<Real> &times);
    /**
     * Clears the list of event times.
     */
    void clearEventTimes();
    /**
     * @return Returns the list of event times.
     */
    inline const std::vector<Real> &eventTimes() const;

    // compartments
    /**
     * Adds a new compartment to the model and takes ownership of it
     * (that is, the model will delete it on model destruction).
     * @param compartment The new compartemnt
     */
    virtual void addCompartment(Compartment *compartment);
    /**
     * Retrieves a compartment from the model by id.
     *
     * @param id Id of the compartment
     * @return The comaprtment
     */
    Compartment *addCompartment(const std::string &id);
    /**
     * Creates a new compartment with the given dimensions and id and
     * adds it to the model.
     *
     * @param id Id of the new compartment
     * @param x0, y0 Position of compartment upper-left corner
     * @param x1, y1 Width and height of new compartment
     * @return Pointer to the new compartment
     */
    Compartment *addCompartment(const std::string &id, int x0, int y0, int x1, int y1);
    /**
     * Removes a compartment from the model
     *
     * @param compartment The compartment to be removed.
     */
    void removeCompartment(Compartment *compartment);
    /**
     * @return Number of compartments in the model.
     */
    inline size_t numCompartments() const;
    /**
     * @return Compartment list.
     */
    inline const boost::container::stable_vector<Compartment *> &compartments() const;
    /**
     * Gets compartment by id.
     *
     * @param id id
     * @return Compartment
     */
    inline Compartment *compartment(const std::string &id) const;

    // species
    /**
     * Adds a new species to the model and takes ownership of it
     * (that is, the model will delete it on model destruction).
     * @param species The species
     */
    virtual void addSpecies(Species *species);
    /**
     * Creates a new species with a given ID and diffusivity and adds it to the model.
     * @param id ID of the new species.
     * @param diffusionConstant Diffusivity of the new species
     * @return Pointer to the new species
     */
    Species *addSpecies(const std::string &id, Real diffusionConstant=0.0);
    /**
     * Removes the given species.
     * @param species The species to be removed.
     */
    void removeSpecies(Species *species);
    /**
     * @return Number of species in the model.
     */
    inline size_t numSpecies() const;
    /**
     * @return Species list.
     */
    inline const boost::container::stable_vector<Species *> &species() const;
    /**
     * Retrieves a model species by id.
     *
     * @param id Species id
     * @return The species
     */
    inline Species *species(const std::string &id) const;

    // individual species
    /**
     * If the model has any individual species, this flag needs to be set.
     * @param hasCPS
     */
    void setHasContinuousParameterSpecies(bool hasCPS=true);
    /**
     * @return True if the model has individual species.
     */
    bool hasContinuousParameterSpecies() const;
    /**
     * Returns true if the reactions in the model do not change the number of individual species.
     *
     * @return True if the reactions in the model do not change the number of individual species.
     */
    bool hasConstantIndividuals() const {return m_hasConstantIndividuals;}

    /**
     * If the reactions in the model do not change the number of individuals (e.g. reactions of the type \f$A\rightarrow A\f$),
     * users can set this value to true. Since the code does not need to keep track of the total number of individuals
     * after each time step, this can speed up the simulation substantially. An example for using
     * this property is MajorityVoteProblem.
     *
     * @param hasConstantIndividuals Set to true, if the reactions in the model do not change the number of individual species.
     */
    void setHasConstantIndividuals(bool hasConstantIndividuals) {m_hasConstantIndividuals = hasConstantIndividuals;}
    /**
     * This property can be used to have individual species behave ballistically at the domain
     * boundaries (i.e. they are reflected from the boundaries like billard balls). If this is
     * set, whenever an individual hits a reflecting boundary, it's velocity component is inversed.
     *
     * @param hasBCS Set to true, if individuals should behave ballistically.
     */
    void setHasBallisticBoundaryConditions(bool hasBCS) {m_hasBallisticBoundaryConditions = hasBCS;}

    /**
     * Returns true, if the individuals behave ballistically at the boundaries.
     *
     * @return True, if the individuals behave ballistically at the boundaries.
     */
    bool hasBallisticBoundaryConditions() const {return m_hasBallisticBoundaryConditions;}

    // reactions
    /**
     * @return Number of reactions in the model.
     */
    inline size_t numReactions() const;
    /**
     * @return List of reactions in the model.
     */
    inline const boost::container::stable_vector<Reaction *> &reactions() const;
    /**
     * Retrieves a reaction by ID.
     * @param id ID of the requested reaction
     * @return Pointer to the reaction
     */
    inline Reaction *reaction(const std::string &id) const;
    /**
     * Adds a new reaction and takes ownership of it
     * (that is, the model will delete it on model destruction)..
     * @param reaction The reaction to add
     */
    virtual void addReaction(Reaction *reaction);
    /**
     * Creates a new reaction with a given ID and adds it to the model.
     * @param id The ID of the new reaction
     * @return Pointer to the newly created reaction.
     */
    Reaction *addReaction(const std::string &id);
    /**
     * Adds a new reaction with a given ID, a given rate and a number of up to eight reactants to the model.
     * This is just a convenience function.
     *
     * @param id ID of the new reaction
     * @param rate The rate of the new reaction
     * @param reactant1, reactant2, reactant3, reactant4, reactant5, reactant6, reactant7, reactant8 The reactants
     * @return Pointer to the newly created reaction.
     */
    Reaction *addReaction(const std::string &id, Real rate,
                           Species *reactant1=0, Species *reactant2=0,
                           Species *reactant3=0, Species *reactant4=0,
                           Species *reactant5=0, Species *reactant6=0,
                           Species *reactant7=0, Species *reactant8=0);

    /**
     * Removes the given reaction from the model.
     * @param reaction
     */
    void removeReaction(Reaction *reaction);

    /**
     * @return True if any of the reactions are localized to a particular compartment. If none of the reactions
     * are localized, the simulation can be sped up considerably.
     */
    bool hasLocalizedReactions() const;
    ///@}

    /** @name Model parameters.
     *  Specifies all parameters that can be used in the various scripts and methods.
     */
    ///@{
    /**
     * @return Number of parameters in the model.
     */
    inline size_t numParameters() const;
    /**
     * @return List of parameters in the model.
     */
    inline const std::map<std::string, Real> &parameters() const;
    /**
     * Gets value of a parameter by name. Note that parameter sweeping is handled in the
     * the Inchman invocation script. In this class, every parameter has a particular value.
     *
     * @param name Name of the parameter
     * @return Its value
     */
    inline Real parameter(const std::string &name) const;
    /**
     * Sets a parameter to a particular value. If the parameter is not yet present in the model,
     * it is automatically added.
     *
     * @param name Name of the parameter
     * @param value Its value
     */
    virtual void setParameter(const std::string &name, Real value);
    /**
     * Removes a parameter from the model.
     *
     * @param name Name of the parameter to be removed.
     */
    virtual void removeParameter(const std::string &name);
    /**
     * Adds a field parameter to the model. The parameter will be available in the init/event scripts. It
     * is also exposed for read/write access in user-defined diffusivity/drift scripts.
     *
     * @param name The name of the field parameter.
     */
    void setFieldParameter(const std::string &name);
    /**
     * Returns a constant reference to the field parameter array.
     */
    std::map<std::string, std::vector< Real > > &getFieldParameters() {return m_fieldParameters;}
    ///@}

    /** @name Open-CL options.
     *  Options that concern the Open-CL parameters.
     */
    ///@{
    /**
     * Sets the path where the kernel template files are located. Defaults to the
     * current working directory.
     * @param path Path to kernel templates.
     */
    void setOclSourcePath(const std::string &path);
    /**
     * @return Path to kernel templates.
     */
    inline const std::string &oclSourcePath() const;

    /**
     * Sets the index of the Open-CL device to use.
     * @param deviceIndex Device index
     */
    void setOclDeviceIndex(size_t deviceIndex);
    /**
     * @return Device index
     */
    inline size_t oclDeviceIndex() const;

    /**
     * If set to true, warnings of the Open-CL compiler are suppressed.
     */
    void setOclInhibitWarnings(bool inhibit);
    /**
     * @return If true, warnings of the Open-CL compiler are suppressed.
     */
    inline bool oclInhibitWarnings() const;
    ///@}

    /** @name Helper functions.
     *  Some helper functions that can be used to set up physical models.
     */
    ///@{
    /**
     * Computes the subvolume size in l.
     * @return Size of subvolume (in l).
     */
    Real subvolumeSize() const;

    /**
     * Computes the number of particles in a subvolume from a concentration
     * given in \f$\mathrm{M}=\mathrm{mol}\,\mathrm{l}^{-1}\f$.
     * @param concentration Concentration in M.
     * @return Particles per subvolume.
     */
    int particleNumberFromConcentration(Real concentration) const;

    /**
     * Computes the reaction rate per subvolume
     * (in \f$\mathrm{s}^{-1}\f$) for a 0th-order
     * reaction from a reaction rate given in \f$\mathrm{M}\,\mathrm{s}^{-1}\f$.
     * @param rate Reaction rate in \f$\mathrm{M}\,\mathrm{s}^{-1}\f$.
     * @return Reaction rate per subvolume in \f$\mathrm{s}^{-1}\f$.
     */
    Real zerothOrderReactionRate(Real rate) const;

    /**
     * Computes the reaction rate per subvolume
     * (in \f$\mathrm{s}^{-1}\f$) for a 2nd-order
     * reaction from a reaction rate given in
     * \f$\mathrm{M}^{-1}\,\mathrm{s}^{-1}\f$.
     * @param rate Reaction rate in \f$\mathrm{M}^{-1}\,\mathrm{s}^{-1}\f$.
     * @return Reaction rate per subvolume in \f$\mathrm{s}^{-1}\f$.
     */
    Real secondOrderReactionRate(Real rate) const;
    ///@}


    /** @name Execution.
     *  All methods related to the execution and output functions.
     */
    ///@{
    /**
     * Solves the model.
     */
    virtual void solve();

    /**
     * Prints the current model to log file.
     */
    void logModel();

    /**
     * Writes all species to the output file. This method can be overridden to allow
     * selective accounting of species. By default, all species are dumped.
     *
     * @param time Dump time.
     * @param isolver Handle to the individual solver to allow extraction of the individual species.
     */
    virtual void writeAllSpecies(Real time, gpgmp::IndividualSolver *isolver=0);

    /**
     * Writes a single time dump to the current output HDF5 file. Inheriting classes can override these
     * methods to do data post-processing before the data is written on disk (see MajorityVoteProblem for an
     * example)
     *
     * @param currentDataGroup The current HDF5 data group
     * @param buffer Pointer to the data buffer (this might be Real or int, depending on the solver chosen)
     * @param bufferStep Size (in byte) of the data buffer per species
     * @param bufferType Data type of the buffer (either H5T_NATIVE_INT or HDF_REAL)
     */
    virtual void writeSingleTimeHDF5(hid_t currentDataGroup, const char *buffer, size_t bufferStep, hid_t bufferType);
    ///@}

    /** @name Internal functions.
     */
    ///@{
    /**
     * @return If true, the Python instance is managed by this class.
     */
    inline bool managesPythonInstance() const;
    /**
     * @return The current state array for stochastic runs.
     */
    inline std::vector<int> &runState();
    /**
     * @return The current state array for deterministic runs.
     */
    inline std::vector<Real>&runDeterministicState();
    /**
     * @return The random seed for this run.
     */
    unsigned int getRandomSeed() const {return m_randomSeed;}
    ///@}


protected:
    /**
     * Will be called to generate the initial states. This method can be overriden to implemement complex
     * generation of the initial states. The default implementation distributes the particles according
     * to the compartment prescriptions and then calls the init script if it is given.
     *
     * @param states The state vector that needs to be filled for stochastic simulations.
     * @param realStates The state vector that needs to be filled for deterministic simulations.
     */
    virtual void generateInitialStates(std::vector<int> &states, std::vector<Real> *realStates=0);

    /**
     * Calls the init script. In particular, the default implementation injects all parameters etc. into
     * the python namespace and makes the state variables available to the script.
     *
     * @param states The state vector that needs to be filled for stochastic simulations.
     * @param realStates The state vector that needs to be filled for deterministic simulations.
     */
    virtual void callInitScript(std::vector<int> &states, std::vector<Real> *realStates);
    
    /**
     * Creates the output file. Will be called only once before start.
     */
    virtual void createOutputFile();    
    /**
     * Opens an existing output file.
     */
    virtual void openOutputFile();
    /**
     * Prepares the output file for writing. This will be called
     * before a new run is executed.
     */
    virtual void prepareOutputFile();
    /**
     * Closes the output file.
     */
    virtual void closeOutputFile();

    /**
     * Execute an alert event from the core.
     *
     * @param time Simulation time stamp.
     */
    virtual void eventCallback(Real time);
    
    /**
     * Helper function to load a text file into a string array.
     *
     * @param contentsOut Reference to the string class
     * @param filename Name of the text file.
     */
    void loadTextFile(std::string &contentsOut, const std::string &filename) const;
    
    /**
     * Creates the reaction mask for localized reactions as an Open-CL image object.
     *
     * @return Handle to the image
     */
    cl_mem oclCreateReactionMask() const;

    /**
     * Calls the Open-CL solver.
     */
    void oclSolve();


private:
    // in the numpy API for python version > 3.x it returns with a value instead of void..
    // weird.
#if PY_VERSION_HEX >= 0x03000000
#define INIT_RETURN_VALUE int
#define HAS_INIT_RETURN_VALUE
#else
#define INIT_RETURN_VALUE void
#endif
    INIT_RETURN_VALUE init();

    Real    m_physicalLength; ///< Physical length of the integration domain in \f$\mu\mathrm{m}\f$.
    size_t  m_gridWidth; ///< Number of subvolumes in direction X.
    size_t  m_gridHeight; ///< Number of subvolumes in direction Y.
    
    cl_int4 m_boundaryMask; ///< Boundary conditions for the integration domain.
    cl_int4 m_sourceMask;   ///< For particle source BCs, holds the particle number.
    
    SolverType m_solver; ///< solver to use.
    
    const bool m_managesPythonInstance; ///< If true, this class manages the python instance
    bool m_nonlinearDiffusivity; ///< If true, recomputes drift and diffusion fields after each diffusion sweep
    bool m_computeMoments; ///< If true, compute the moments of each species after each diffusion sweep
    bool m_regenerateInitStates; ///< if true, the init states are re-generated after each run

    Real m_p0; ///< Gives the probability for a particle to jump to a neighbouring subvolume in the stochastic homogeneous solver.

    bool m_hasContinuousParameterSpecies; ///< If true, model contains species with continuous parameters
    bool m_hasConstantIndividuals; ///< True if reactions do not change the number of individuals
    bool m_hasBallisticBoundaryConditions; ///< True if the individuals should behave ballistically at the boundary.

    unsigned m_numRuns; ///< Number of runs to do for this experiment.
    int      m_maxSteps; ///< Number of steps for each run.
    Real     m_maxTime; ///< Maximum simulation runtime.
    
    unsigned m_runsOffset; ///< Offset for the run number (deprecated)

    std::string m_oclSourcePath; ///< Path to the kernel template files.
    size_t      m_oclDeviceIndex; ///< Index of the Open-CL device that is used.
    bool        m_oclInhibitWarnings; ///< If true, Open-CL compiler warnings are suppressed.
    
    Real          m_outputInterval; ///< Interval at which the state will be written (in simulation time).    
    OutputFormat  m_outputFormat;   ///< Current output format.
    std::string   m_outputFilename; ///< Output filename.
    std::string   m_outputJobId;    ///< If present, then the HDF5 runs will be written to /{m_jobId}/Run_{
    
    
    std::string m_initScriptContents; ///< The init script.
    std::string m_computeDriftDiffusivityMethod; ///< The drift/diffusivity method
    std::string m_setFieldMethod; ///< The field update method
    std::string m_newIndividualsMethod; ///< The method to initialize newly created individuals

    std::vector<Real> m_eventTimes; ///< Contains the times for the alert events.
    std::string m_eventsScriptContents; ///< The event script
    
    std::map<std::string, Real>                      m_parameters; ///< Model parameters
    boost::container::stable_vector<Compartment *>   m_compartments; ///< Model compartments
    boost::container::stable_vector<Species *>       m_species; ///< Species in the model
    boost::container::stable_vector<Reaction *>      m_reactions; ///< Reactions in the model
    std::map<std::string, std::vector<Real> > m_fieldParameters; ///< Names of all field parameters.

    // Caches
    std::map<std::string, Species *>     c_idToSpecies; ///< Used by the solver to convert integer id to a species
    std::map<std::string, Reaction *>    c_idToReaction; ///< Used by the solver to convert integer id to a reaction
    std::map<std::string, Compartment *> c_idToCompartment; ///< Used by the solver to convert integer id to a compartment
    
    cl::Context      *m_context; ///< Device context for the Open-CL device
    cl::Device        c_contextDevice;
    cl::CommandQueue *m_commandQueue; ///< Open-CL command queue.
    
    // model information
    std::vector<int> m_runState; ///< Pointer to the state array held for the species.
    std::vector<Real> m_runDeterministicState; ///< Pointer to state array from deterministic sim
    
    // runtime parameters
    unsigned m_runIndex; ///< Number of current run.
    unsigned m_runDumpIndex; ///< Number of output data set in current run.
    
    // hdf5 output data
    hid_t m_hdf5File; ///< File id for current hdf5 output.
    hid_t m_hdf5CurrentGroup; ///< Group for current run.

    unsigned int m_randomSeed; ///< RNG seed
};

// inlines
inline void DiffusionModel::setHasContinuousParameterSpecies(bool hasCPS) {m_hasContinuousParameterSpecies = hasCPS;}
inline bool DiffusionModel::hasContinuousParameterSpecies()  const {return m_hasContinuousParameterSpecies;}

inline Real DiffusionModel::physicalLength() const { return m_physicalLength; }
    
inline size_t DiffusionModel::gridArea() const { return m_gridWidth * m_gridHeight; }
inline size_t DiffusionModel::gridWidth() const { return m_gridWidth; }
inline size_t DiffusionModel::gridHeight() const { return m_gridHeight; }
    
inline const cl_int4 & DiffusionModel::boundaryMask() const { return m_boundaryMask; }
inline const cl_int4 & DiffusionModel::sourceMask() const { return m_sourceMask; }

// todo: rename these to solverType ..
inline SolverType DiffusionModel::solver() const { return m_solver; }
    
inline Real DiffusionModel::p0() const { return m_p0; }
    
inline unsigned DiffusionModel::numRuns() const {
    return m_numRuns;
}
inline int DiffusionModel::maxSteps() const {
    return m_maxSteps;
}
inline Real DiffusionModel::maxTime() const {
    return m_maxTime;
}

inline std::vector<int> & DiffusionModel::runState() {return m_runState;}

inline std::vector<Real> & DiffusionModel::runDeterministicState() {return m_runDeterministicState;}

inline const std::string & DiffusionModel::oclSourcePath() const {
    return m_oclSourcePath;
}

inline size_t DiffusionModel::oclDeviceIndex() const {
    return m_oclDeviceIndex;
}
    
inline bool DiffusionModel::oclInhibitWarnings() const {
    return m_oclInhibitWarnings;
}

inline Real DiffusionModel::outputInterval() const {
    return m_outputInterval;
}
inline OutputFormat DiffusionModel::outputFormat() const {
    return m_outputFormat;
}
inline const std::string & DiffusionModel::outputFilename() const {
    return m_outputFilename;
}

inline bool DiffusionModel::managesPythonInstance() const {
    return m_managesPythonInstance;
}    

inline bool DiffusionModel::nonlinearDiffusivity() const {
    return m_nonlinearDiffusivity;
}

inline void DiffusionModel::setNonlinearDiffusivity(bool nonlinear) {
    m_nonlinearDiffusivity = nonlinear;
}

inline bool DiffusionModel::computeMoments() const {
    return m_computeMoments;
}

inline void DiffusionModel::setComputeMoments(bool computeMoments) {
    m_computeMoments = computeMoments;
}

inline const std::vector<Real> & DiffusionModel::eventTimes() const {
    return m_eventTimes;
}
    
inline size_t DiffusionModel::numParameters() const { return m_parameters.size(); }
inline const std::map<std::string, Real> & DiffusionModel::parameters() const { return m_parameters; }
inline Real DiffusionModel::parameter(const std::string &name) const {
    std::map<std::string, Real>::const_iterator i = m_parameters.find(name);
    return (i != m_parameters.end()
            ? i->second
            : 0); // FUTURE: handle this better
}
    
inline size_t DiffusionModel::numCompartments() const { return m_compartments.size(); }
inline const boost::container::stable_vector<Compartment *> & DiffusionModel::compartments() const {
    return m_compartments;
}
inline Compartment * DiffusionModel::compartment(const std::string &id) const {
    std::map<std::string, Compartment *>::const_iterator i = c_idToCompartment.find(id);
    return (i != c_idToCompartment.end()
            ? i->second
            : 0);
}

inline size_t DiffusionModel::numSpecies() const { return m_species.size(); }        
inline const boost::container::stable_vector<Species *> & DiffusionModel::species() const {
    return m_species;
} 
inline Species * DiffusionModel::species(const std::string &id) const {
    std::map<std::string, Species *>::const_iterator i = c_idToSpecies.find(id);
    return (i != c_idToSpecies.end()
            ? i->second
            : 0);
}

inline size_t DiffusionModel::numReactions() const { return m_reactions.size(); }
inline const boost::container::stable_vector<Reaction *> & DiffusionModel::reactions() const {
    return m_reactions;
}
inline Reaction * DiffusionModel::reaction(const std::string &id) const {
    std::map<std::string, Reaction *>::const_iterator i = c_idToReaction.find(id);
    return (i != c_idToReaction.end()
            ? i->second
            : 0);
}

 const std::string &DiffusionModel::computeDriftDiffusivityMethod() const {
   return m_computeDriftDiffusivityMethod;
 }

 const std::string &DiffusionModel::returnSetFieldMethod() const {
     return m_setFieldMethod;
 }

void DiffusionModel::setComputeDriftDiffusivityMethod(const std::string &method) {
    m_computeDriftDiffusivityMethod = method;
}

void DiffusionModel::setFieldMethod(const std::string &script) {
    m_setFieldMethod = script;
}

void DiffusionModel::setNewIndividualsMethod(const std::string &method) {
    m_newIndividualsMethod = method;
}

const std::string &DiffusionModel::getNewIndividualsMethod() const {
    return m_newIndividualsMethod;
}



} // namespace gpgmp


#endif // !__gpgmp_DiffusionModel_h__
