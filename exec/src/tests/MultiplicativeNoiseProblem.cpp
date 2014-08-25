#include "MultiplicativeNoiseProblem.h"
#include "Compartment.h"

namespace gpgmp {

/**
  Creates the species and one compartment "Source".
  We use the pre-defined drift-diffusion method gpgmp::DT_MULTIPLICATIVE_NOISE and set the kernel
  parameters accordingly.

  \see setComputeDriftDiffusivityMethod
  \see setParameter
  */
MultiplicativeNoiseProblem::MultiplicativeNoiseProblem(Real length, int dx, int dy, int numMolecules,
                                                       Real cx, Real cy, Real mux, Real muy)
  : DiffusionModel("Multiplicative Noise Problem", true)
{
  // set dimensions
  setGridDims(dx,dy);
  setPhysicalLength(length);

    // one species
    Species *speciesA = addSpecies("A", 1.);

    // set boundary mask
    const cl_int4 boundaryMask = {{1, 1, 1, 1}};
    const cl_int4 sourceMask   = {{0, 0, 0, 0}};
    setBoundaryMasks(boundaryMask, sourceMask);

    // add source compartment
    Compartment *compartmentSource = new Compartment("Source", dx/4,dy/4, dx/4, dy/4);
    compartmentSource->setInitialAmount(speciesA, numMolecules, gpgmp::HomogeneousDistribution);
    addCompartment(compartmentSource);

    // add user parameters for kernel
    setParameter("sigmaX", cx);
    setParameter("sigmaY", cy);
    setParameter("muX", mux);
    setParameter("muY", muy);

    // set homogeneous diffusivity and drift
    setComputeDriftDiffusivityMethod(gpgmp::DT_MULTIPLICATIVE_NOISE);
} // constructor

}// namespace gpgmp
