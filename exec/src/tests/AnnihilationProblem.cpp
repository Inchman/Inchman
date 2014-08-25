/*
 * AnnihilationModel.cpp
 *
 *  Created on: 03/03/2010
 *      Author: matthias
 */

#include "AnnihilationProblem.h"

#include "Compartment.h"


namespace gpgmp {

/**
  If a non-vanishing drift constant is specified, the constructor will set up
  an inhomogeneous simulation. We use the pre-defined constant drift-diffusion
  method gpgmp::DT_HOMOGENEOUS and set up the kernel parameters accordingly.

  \code
    setParameter("diffx", diffusionConstant);
    setParameter("diffy", diffusionConstant);
    setParameter("rx", drift);
    setParameter("ry", drift);

    setComputeDriftDiffusivityMethod(gpgmp::DT_HOMOGENEOUS);
  \endcode
  */
AnnihilationProblem::AnnihilationProblem(Real length, int dx, int dy,
                                         int numMolecules,
                                         Real diffusionConstant,
                                         Real reactionRate,
                                         Real drift)
    :   DiffusionModel("AB Annihilation", true)
{
  std::cout <<"Constructing.."<<std::flush;
  setPhysicalLength(length);
  setGridDims(dx, dy);

	// two species
    Species *speciesA = addSpecies("A", diffusionConstant);
    Species *speciesB = addSpecies("B", diffusionConstant);

    // set boundary conditions to outflow
    const cl_int4 boundaryMask = {{1, 1, 1, 1}};
    const cl_int4 sourceMask   = {{0, 0, 0, 0}};
    setBoundaryMasks(boundaryMask, sourceMask);

    // create World compartment (background)
    Compartment *world = addCompartment("World", 0, 0, dx-1, dy-1);
    world->setInitialAmount(speciesA, numMolecules, HomogeneousDistribution);

    // create Source compartment
	// we want to set up an equivalent problem to MesoRD
	// so we need to make a 2x2 source compartment in the centre
    Compartment *source = addCompartment("Source", dx/2-1, dy/2-1, dx/2, dy/2);
    source->setInitialAmount(speciesB, numMolecules, HomogeneousDistribution);

	// add the annihilation reaction
    addReaction("A plus B annihilation",
                secondOrderReactionRate(reactionRate),
                speciesA, speciesB);

    // set up inhomogeneous simulation if drift field is given
    if (drift != 0) {
        // add user parameters for kernel
        setParameter("diffX", diffusionConstant);
        setParameter("diffY", diffusionConstant);
        setParameter("driftX", drift);
        setParameter("driftY", drift);

        // set homogeneous diffusivity and drift
        setComputeDriftDiffusivityMethod(gpgmp::DT_HOMOGENEOUS);
    }
  std::cout <<"Constructor finished.\n"<<std::flush;
} // constructor


} // namespace gpgmp
