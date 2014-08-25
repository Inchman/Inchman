#include "NonlinearProblem.h"
#include "Compartment.h"

namespace gpgmp {

/**
    The constructor sets up one species \f$A\f$ and a source compartment. We
    use the pre-define drift-diffusion method gpgmp::DT_NONLINEAR. This method
    needs to know the first moment (average) of the particle distribution and
    we hence call DiffusionModel::setComputeMoments. We also need to re-compute
    the drift-diffusion field after each time step. This is achieved by
    calling DiffusionModel::setNonlinearDiffusivity.

    \see DiffusionModel::setComputeMoments
    \see DiffusionModel::setNonlinearDiffusivity
  */
NonlinearProblem::NonlinearProblem(Real length, int dx, int dy,
                                   Real diffX, Real diffY,
                                   Real theta, Real omega,
                                   int numMolecules)
    : DiffusionModel(length, dx, dy)
{
    // one species
    Species *speciesA = addSpecies("A", 1.);

    // set boundary mask
    const cl_int4 boundaryMask = {{1, 1, 1, 1}};
    const cl_int4 sourceMask   = {{0, 0, 0, 0}};
    setBoundaryMasks(boundaryMask, sourceMask);

    // add source compartment
    Compartment *compartmentSource = new Compartment("TrueSource", dx/2, dy/2, dx/2, dy/2);
    compartmentSource->setInitialAmount(speciesA, numMolecules, gpgmp::HomogeneousDistribution);
    addCompartment(compartmentSource);

    // add user parameters for kernel
    setParameter("diffX", diffX);
    setParameter("diffY", 0);
    setParameter("theta", theta);
    setParameter("omega", omega);
    setParameter("origin", length/8.);

    // set non-linear diffusivity and drift
     setComputeDriftDiffusivityMethod(gpgmp::DT_NONLINEAR);

     // we also need to re-compute the moments and the diffusivity
     // after each time step as it's a non-linear problem
     setNonlinearDiffusivity(true);
     setComputeMoments(true);
} // constructor

} // namespace gpgmp
