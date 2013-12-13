#ifndef INHOMOGENEOUSSOLVER_H
#define INHOMOGENEOUSSOLVER_H

#include "Solver.h"

namespace gpgmp{
/**
 * This class implements an inhomogeneous stochastic solver. It deals with inhomogeneous diffusivity and
 * drift and provides a stochastic and deterministic solver.
 */
class InhomogeneousSolver : public Solver
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
     * @param deterministic True if it's a deterministic simulation
     * @param hasIndividualSpecies True if the model involves individual species
     * @param computeMoments True if the moments need to be computed before each time step
     * @param nonlinearDiffusivivity True if the drift/diffusivity needs to be re-computed after each step
     */
    InhomogeneousSolver(DiffusionModel &diffusionModel,
                        std::vector<int> *runState,
                        std::vector<Real> *deterministicState,
                        SolverType solverType,
                        bool deterministic, bool hasIndividualSpecies, bool computeMoments,
                        bool nonlinearDiffusivity);

    /**
     * Destructor. Deletes all device buffers.
     *
     * @todo: what about d_states and d_lostParticles?
     */
    ~InhomogeneousSolver();

    /**
     * The implementation of Solver::getNextDiffusionTime() (see there for description).
     *
     * @return The absolute time for the next diffusion step indexed by species
     */
    std::vector<Real> getNextDiffusionTime();
    /**
     * The implementation of Solver::diffuse(int, Real, int &) (see there for description).
     *
     */
    void diffuse(int speciesIndex, Real currentTime, int &seed);


private:

    // methods
    /**
     * Computes the X moment for species speciesIndex. The moments are held on the device
     * in d_sumMoments and d_sumStates (for the moment and the particle count, respectively).
     *
     * @param speciesIndex Species for which to compute the moment
     */
    void computeMoments(cl_int speciesIndex);

    /**
     * Computes the diffusivity and drift field for species speciesIndex. Since this might
     * depend on the absolute simulation time, this needs to be given too. The diffusivity
     * and drift are stored in the device arrays d_diffusionConstantsX, d_diffusionConstantsY,
     * d_rx, and r_ry, respectively.
     *
     * @param speciesIndex Index for the species
     * @param simTime Absolute simulation time
     *
     */
    void computeDiffusivityDrift(cl_int speciesIndex, Real simTime);

    /**
     * Computes the timestep in X direction for the species speciesIndex. This method
     * assumes that the diffusivity and drift field is set properly for each individual cell on
     * the device [see computeDiffusivityDrift(cl_int, double)]. Since this might
     * depend on the absolute simulation time, the simulation time needs to be given too.
     * The routine returns the global minimum time step for this species in X direction and
     * updates it in m_smallestDtX.
     *
     * @param speciesIndex Index for the species
     * @param simTime Absolute simulation time
     * @return Global minimum time step for this species in X direction.
     */
    Real computeTimestepX(cl_int speciesIndex, Real simTime);

    /**
     * Computes the timestep in Y direction for the species speciesIndex. This method
     * assumes that the diffusivity and drift field is set properly for each individual cell on
     * the device [see computeDiffusivityDrift(cl_int, double)]. Since this might
     * depend on the absolute simulation time, the simulation time needs to be given too.
     * The routine returns the global minimum time step for this species in X direction and
     * updates it in m_smallestDtY.
     * @param speciesIndex Index for the species
     * @param simTime Absolute simulation time
     * @return Global minimum time step for this species in Y direction.
     */
    Real computeTimestepY(cl_int speciesIndex, Real simTime);


    // attributes
    bool m_computeMoments; ///< True, if the moments need to be computed before each time step
    bool m_nonlinearDiffusivity; ///< True if the diffusivity needs to be recomputed after each time step

    // host array for smallest diffusion time step
    std::vector<Real> m_smallestDtX; ///< Global minimum time step in X direction indexed by species
    std::vector<Real> m_smallestDtY; ///< Global minimum time step in X direction indexed by species

    // host array for next diffusion time for each species in each dir - with init to 0.
    std::vector<Real> m_nextDiffusionTimeX; ///< Host array for next diffusion time for each species in x dir
    std::vector<Real> m_nextDiffusionTimeY;///< Host array for next diffusion time for each species in y dir

    // host array for diffusion time step for each species in each dir
    // we need those for the diffusion sweeps
    std::vector<Real> m_diffusionDtX; ///< Host array for diffusion time step for each species in x dir
    std::vector<Real> m_diffusionDtY; ///< Host array for diffusion time step for each species in y dir

    // host array for mean and total x
    std::vector<Real> m_meanX; ///< Host array for global summed x indexed by species
    std::vector<Real> m_sumState; ///< Host array for global particle count indexed by species

    cl::Buffer *d_lostParticles; ///< Store the number of lost particles per cell

    // programs
    cl::Program m_stochasticProgram; ///< Program for stochastic solver
    cl::Program m_deterministicProgram; ///< Program for deterministic solver

    // Kernels
    std::vector<Kernel> computeDiffusionConstants_forSpecies;
    Kernel computeIndividualTimestepX, computeIndividualTimestepY;
    Kernel reduceBlocks;
    Kernel diffuseX, diffuseY;
    Kernel updateStateX, updateStateY;
    Kernel computeMomentsX, reduceMomentsBlocks;

    // these are the deterministic kernels
    Kernel inhomogeneous_deterministic_diffuseX, inhomogeneous_deterministic_diffuseY;

    // drift/diffusivity fields
    cl::Buffer *d_diffusionConstantsX;
    cl::Buffer *d_diffusionConstantsY;
    cl::Buffer *d_rx;
    cl::Buffer *d_ry;

    cl::Buffer *d_timeSteps; ///< Device array for individual time steps
    cl::Buffer *d_minTimeBlock; ///< Device array for minimum time step per block
    cl::Buffer *d_globalMinTimestep; ///< Device array for global minimum time step
    cl::Buffer *d_sumMomentsBlocks; ///< Device array for summed X per block
    cl::Buffer *d_sumStatesBlocks; ///< Device array for particle count per block
    cl::Buffer *d_sumMoments; ///< Device array for global summed x indexed by species
    cl::Buffer *d_sumStates; ///< Device array for global particle count indexed by species

    cl::Buffer *d_transitionProbabilities; ///< Device array for transition probabilities
};
} // namespace gpgmp

#endif // INHOMOGENEOUSSOLVER_H
