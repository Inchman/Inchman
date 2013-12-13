#ifndef SOLVER_H
#define SOLVER_H

#include "Globals.h"
#include "DiffusionModel.h"
#include "IndividualSolver.h"

namespace gpgmp {
/**
 * Solver is an abstract class that defines the interface for all solvers. Derived classes
 * must implement computeNextDiffusionTime() and diffuse() for both stochastic and deterministic simulations.
 * This base class provides a standard implementation of the Gillespie algorithm and alpha-QSS/RK4 integrators
 * to deal with reactions.
 */
class Solver
{
public:
    /**
     * Initializes the solver. All device buffers, GPU programs, kernels etc. need to be initialized here.
     * We need a handle to the calling diffusion model in order to have access to the parameters.
     *
     * @todo: maybe have a diffusion model parameter structure so we don't need the handle..
     *
     * @param runState Handle to the stochastic state array
     * @param deterministicState Handle to the deterministic state array
     * @param diffusionModel Handle to the calling diffusionModel.
     * @param solverType The type of the solver
     * @param deterministic
     * @param hasIndividualSpecies
     *
     * @todo Need to put in check if the GPU resources are sufficient (shared mem etc..)
     */
    Solver(DiffusionModel &diffusionModel,
           std::vector<int> *runState,
           std::vector<Real> *deterministicState, SolverType solverType,
           bool deterministic=false, bool hasIndividualSpecies=false);

    /**
     * Destructor. Cleans up all the CL stuff.
     */
    virtual ~Solver();

    /**
     * Performs the reactions in each cell for the time step deltat using the
     * given seed. Depending on the chosen algorithm, this implementation uses
     * Gillespie (for stochastic simulations) or alpha-QSS/RK4 (for deterministic).
     * The seed will be updated.
     *
     * @param deltat The time step
     * @param seed Seed
     */
    void doReactions(Real deltat, int &seed);

    /**
     * Updates the field with the user-defined kernel.
     *
     * @param deltat The time step.
     */
    void updateFields(Real deltat, Real simtime);

    /**
     * Returns the absolute times for the next diffusion step indexed by species.
     * This array needs to be maintained by all deriving classes and updated
     * after each diffusion step.
     *
     * @return Absolute times for the next diffusion step indexed by species.
     */
    virtual std::vector<Real> getNextDiffusionTime() = 0;


    /**
     * Deriving classes need to implement this method to diffuse the given species
     * speciesIndex. The current simulation time is provided in currentTime and
     * the random number seed needs to be updated.
     *
     * @param speciesIndex Species to be diffused
     * @param currentTime Current simulation time
     * @param seed Random number seed
     */
    virtual void diffuse(int speciesIndex, Real currentTime, int &seed) = 0;

    /**
     * Copies the device state buffer to host memory (for write to file).
     */
    void synchronizeStateBufferFromDevice();

    /**
     * Copies the state buffer from host memory to device (after events)
     */
    void synchronizeStateBufferToDevice();

    /**
     * Copies the field buffer from the device to host memory.
     */
    void synchronizeFieldBufferFromDevice();

    /**
     * Copies the field buffer from host memory to device.
     */
    void synchronizeFieldBufferToDevice();

    /**
     * Returns a handle to the individual solver object (if any).
     *
     * @todo: get rid of that dependency..
     *
     * @return
     */
    IndividualSolver *getIndividualSolver() {return m_individualSolver;}

    /**
     * Adds the common headers for the CL source codes to the input string stream.
     * In particular, it makes the global and user-defined variables
     * available to the kernels.
     *
     * @todo No, it doesn't make the variables available to the kernels! It should though..
     *
     * @param ss The stringstream the headers are added to
     */
    void oclGenerateHeader(std::ostringstream &ss) const;

protected:
    // members

    // helper functions copied from diffusion model
    // todo: would be nice to have a proper buffer class..
    template <class T> cl::Buffer *oclCreateBuffer(size_t count, cl_mem_flags flags = CL_MEM_READ_WRITE) const;
    template <class ContainerT> cl::Buffer *oclCreateBufferFrom(ContainerT &from, cl_mem_flags flags = CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR) const;
    template <class T> cl::Buffer *oclCreateBufferFrom(T *from, size_t count, cl_mem_flags flags) const;
    template <class ContainerT> void oclCopy(const cl::Buffer &from, ContainerT &to);
    template <class T> void oclCopy(const cl::Buffer &from, T *to, size_t count);
    template <class ContainerT> void oclCopy(const ContainerT &from, cl::Buffer &to);
    template <class T> void oclCopy(const T *from, cl::Buffer &to, size_t count);
    void oclCopyDeviceBuffers(cl::Buffer &from, cl::Buffer &to, size_t count);

    /**
     * Builds a program from the given source.
     *
     * @param ss The source
     * @return The program
     */
    cl::Program oclBuildProgram(std::ostringstream &ss) const;


    /**
     *  Adds the contents of a template file to the input stream and replaces all tags.
     *
     *  @todo: maybe split it up into common part for reactions and solver-dependent ones for diffusion?
     *
     *  @param ss The stringstream
     *  @param clFileName Name of the template file
     **/
    void oclReadAndReplaceTemplateFile(std::ostringstream &ss, const char *clFileName) const;

    /**
     * Initialises the CL device.
     *
     * @todo: we should use our new classes for that..
     */
    void oclInit();

    /**
     * Deletes the CL device.
     *
     * @todo: use our new classes
     */
    void oclDestroy();

    /**
     * Callback function for CL calls. Prints out the CL info to std::error.
     */
    static void CL_CALLBACK oclNotification(const char *errinfo, const void * /*private_info*/,
                                            size_t /*cb*/, void * /*user_data*/);

    /**
     * Creates the localized reaction mask in device memory.
     *
     * @return Pointer to device buffer for the localized reaction mask
     */
    cl_mem oclCreateReactionMask() const;

    /**
     * Checks if a kernel call resulted in an kernel error. Only works if m_debugKernel is true.
     *
     * @return True, if a kernel error occurred.
     */
    bool checkKernelError();

    // attributes
    DiffusionModel &m_diffusionModel;
    const bool m_debugKernel = false; ///< Set to true if you want very heavy-weight kernel debug code

    // cl files
    cl::Context      *m_context; ///< device context
    cl::Device        c_contextDevice; ///< The current device.
    cl::CommandQueue *m_commandQueue; ///< command queue

    // Buffers
    cl::Buffer *d_state; ///< The device buffer for the state
    std::vector<cl_int>  *m_runState; ///< Pointer to state for stochastic runs
    std::vector<Real> *m_runDeterministicState; ///< Pointer to state for deterministic runs
    cl::Buffer *d_fields; ///< The device buffer for the field variables

    bool m_hasFields; ///< Will be set to true if the solver has field parameters
    std::map<std::string, size_t> m_fieldOffsets; ///< During a run, holds the offsets for the field buffers

    SolverType m_solverType; ///< Type of solver

    bool m_deterministic; ///< True, if the solver should perform a deterministic simulation
    bool m_hasIndividualSpecies; ///< True, if the model has individual species
    bool m_hasLocalizedReactions; ///< True if the model has localized reactions

    // topology for GPU
    cl::NDRange m_globalRange; ///< Global range (whole domain)
    cl::NDRange m_localRange; ///< Local range (one block)
    cl::NDRange m_globalRangeReduce; ///< Global range for reduce over blocks operation
    int m_localSizeReduce; ///< Local size for reduce over blocks operation

    // stuff for the individual solver
    IndividualSolver *m_individualSolver; ///< Handle to the individual solver if needed
    std::vector<bool> m_individualSpecies; ///< Elements are set true, if the species is individual

    // for the alpha-QSS solver
    cl::Buffer *d_y0, *d_ys, *d_y1, *d_ym1, *d_ym2, *d_scrarray, *d_qs, *d_rtaus;

    cl_mem d_reactionMask; ///< Buffer for the localized reaction mask

    std::vector<Real> m_errors; ///< Host array for kernel errors
    cl::Buffer *d_errors; ///< Device buffer for kernel errors

    // Programs
    cl::Program m_gillespieProgram; ///< Program for the Gillespie kernel
    cl::Program m_deterministicProgram; ///< Program for the deterministic reaction kernels

    // Kernels
    Kernel gillespie; ///< Gillespie kernel
    Kernel deterministic_performReaction_rk4, deterministic_performReaction_aqss;
    Kernel m_updateFields; ///< Kernel to update the field variables.
};

template <class ContainerT>
void Solver::oclCopy(const cl::Buffer &from, ContainerT &to) {
    assert(!to.empty());
    oclCopy(from, &to.front(), to.size());
}
template <class T>
void Solver::oclCopy(const cl::Buffer &from, T *to, size_t count) {
    m_commandQueue->enqueueReadBuffer(from, CL_TRUE, 0, count * sizeof(T), (void *)to);
    m_commandQueue->flush();
}

template <class ContainerT>
void Solver::oclCopy(const ContainerT &from, cl::Buffer &to) {
    assert(!from.empty());
    oclCopy(&from.front(), to, from.size());
}
template <class T>
void Solver::oclCopy(const T* from, cl::Buffer &to, size_t count) {
    m_commandQueue->enqueueWriteBuffer(to, CL_TRUE, 0, count * sizeof(T), (const void *)from);
    m_commandQueue->flush();
}

template <class T>
cl::Buffer *Solver::oclCreateBuffer(size_t count, cl_mem_flags flags) const {
    return oclCreateBufferFrom((T *)0, count, flags);
}

template <class ContainerT>
cl::Buffer *Solver::oclCreateBufferFrom(ContainerT &from, cl_mem_flags flags) const {
    return oclCreateBufferFrom(&from.front(), from.size(), flags);
}

template <class T>
cl::Buffer *Solver::oclCreateBufferFrom(T *from, size_t count, cl_mem_flags flags) const
{
    return new cl::Buffer(*m_context, flags, count * sizeof(T), (void *)from);
}
} // namespace gpgmp

#endif // SOLVER_H
