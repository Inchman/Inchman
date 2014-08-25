#include "OrnsteinUhlenbeckProblem.h"
#include "Compartment.h"


namespace gpgmp {
/**
  Sets up the species and the source compartment. We use the pre-defined
  drift-diffusion method gpgmp::DT_ORNSTEIN_UHLENBECK and set the kernel
  parameters accordingly.

  \see setComputeDriftDiffusivityMethod
  \see setParameter
  */
OrnsteinUhlenbeckProblem::OrnsteinUhlenbeckProblem(Real length, int dx, int dy,
                                                   Real diffusivityX, Real diffusivityY,
                                                   Real gxx, Real gxy, Real gyx, Real gyy,
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
    Compartment *compartmentSource = new Compartment("TrueSource", dx/4,dy/4, dx/4, dy/4);
    compartmentSource->setInitialAmount(speciesA, numMolecules, gpgmp::HomogeneousDistribution);
    addCompartment(compartmentSource);

    // add user parameters for kernel
    setParameter("diffX", diffusivityX);
    setParameter("diffY", diffusivityY);
    setParameter("gammaXX", gxx);
    setParameter("gammaXY", gxy);
    setParameter("gammaYX", gyx);
    setParameter("gammaYY", gyy);

    // set OU method
    setComputeDriftDiffusivityMethod(gpgmp::DT_ORNSTEIN_UHLENBECK);
}
} // namespace gpgmp
