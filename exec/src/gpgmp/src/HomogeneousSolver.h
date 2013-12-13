#ifndef HOMOGENEOUSSOLVER_H
#define HOMOGENEOUSSOLVER_H

#include "Solver.h"

namespace gpgmp{
/**
 * This class implements an inhomogeneous stochastic solver. It provides a stochastic and deterministic solver as described in Vigelius et al. (2010).
 */
class HomogeneousSolver : public Solver
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
     */
    HomogeneousSolver(DiffusionModel &diffusionModel,
                      std::vector<int> *runState,
                      std::vector<Real> *deterministicState,
                      SolverType solverType,
                      bool deterministic, bool hasIndividualSpecies);

    /**
     * Destructor. Destroys all diffusion-related buffers (only d_lostParticles at the moment).
     */
    ~HomogeneousSolver();

    /**
     * Returns the absolute next diffusion time indexed per species. This
     * is an implementation of Solver::getNextDiffusionTime() [see there for more info].
     *
     * @return Next diffusion time.
     */
    std::vector<Real> getNextDiffusionTime() {return m_nextDiffusionTime;}

    /**
     * Diffuses the given species. This is an
     * implementation of Solver::diffuse(int, Real, int &) [see there for info].
     *
     * @param speciesIndex The species to diffuse
     * @param currentTime Simulation time
     * @param seed RNG seed
     */
    void diffuse(int speciesIndex, Real currentTime, int &seed);

private:

    // methods

    // Attributes

    // Programs
    cl::Program m_diffusionProgram; ///< Program containing the diffusion kernels

    // Kernels
    Kernel stochastic_diffuse; ///< Kernel to diffuse particles
    Kernel stochastic_updateState; ///< Updates the current cell with the incoming particles from the neighbouring cell
    Kernel stochastic_source_boundary; ///< Inflow (source) boundary
    Kernel homogeneous_deterministic_diffuse; ///< Deterministic diffusion kernel

    /// array for diffusion time steps
    std::vector<Real> m_dt;

    /// array for absolute next diffusion time
    std::vector<Real> m_nextDiffusionTime;

    // device buffers
    cl::Buffer *d_lostParticles; ///< Store the number of lost particles per cell

};
}

#endif // HOMOGENEOUSSOLVER_H
