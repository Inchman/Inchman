/*
* FisherProblem.cpp
*
*  Created on: 11/03/2010
*      Author: matthias
*/

#include "ReactionProblem.h"

#include <iostream>

#include "Reaction.h"


namespace gpgmp {

/**
  Sets up the species, reactions and a "World" compartment.
  */
ReactionProblem::ReactionProblem(Real length, int dx, int dy,
                                   Real diffusionConstant,
                                   Real k1, Real k2, Real k3, Real k4)
:   DiffusionModel(length, dx, dy)
{
    // add both species to main model
    Species *speciesA = addSpecies("A", diffusionConstant);
    Species *speciesB = addSpecies("B", diffusionConstant);

    // add A+A -> 0 annihilation reaction
    Reaction *aplusaAnnihilation = addReaction("A plus A annihilation",
                                               k1, speciesA, speciesA);
    std::cout << "Added reaction "<<aplusaAnnihilation->id() << std::endl;

    // add A+B -> 0 annihilation reaction
    Reaction *aplusbAnnihilation = addReaction("A plus B annihilation",
                                               k2, speciesA, speciesB);
    std::cout << "Added reaction "<<aplusbAnnihilation->id() << std::endl;

    // add A creation reaction
    Reaction *aCreation = addReaction("A creation", k3);
    aCreation-> setProductStoichiometry(speciesA, 1);
    std::cout << "Added reaction " << aCreation->id() << std::endl;

    // add B creation reaction
    Reaction *bCreation = addReaction("B creation", k4);
    bCreation-> setProductStoichiometry(speciesB, 1);
    std::cout << "Added reaction " << bCreation->id() << std::endl;

    // set diffusivity/drift - needed if started with inhomogeneous solver
    // if we use the inhomogeneous solver, we don't get the same results ..
    // the reason is that the time step is different. we'd need to set
    // p0=1 .. but that sort of breaks the code (todo: why?)
    setParameter("diffX", diffusionConstant);
    setParameter("diffY", diffusionConstant);
    setParameter("driftX", 0.);
    setParameter("driftY", 0.);

    // set homogeneous diffusivity and drift
    setComputeDriftDiffusivityMethod(gpgmp::DT_HOMOGENEOUS);

    // set boundary mask
    clearBoundaryMasks();

    // set up world compartment
    addCompartment("World", 0,0, dx-1, dy-1);
}

} // namespace gpgmp
