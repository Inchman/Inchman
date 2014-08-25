/*
 * LocalizedAnnhilationProblem.cpp
 *
 *  Created on: 09/03/2010
 *      Author: matthias
 */

#include "LocalizedAnnihilationProblem.h"

#include "Compartment.h"

#include <boost/math/special_functions/round.hpp>

#include <cmath>
#include <iostream>

namespace gpgmp {

/**
  Localization is achieved by creating the corresponding compartments and
  constraining the reaction using Reaction::addCompartment():

  \code
    gpgmp::Reaction *aCreation = addReaction("A creation", zerothOrderReactionRate(k3));
    gpgmp::Compartment *compartment1 = new gpgmp::Compartment("A Creation compartment", iA0, 0, iA1, dx);

    aCreation->addCompartment(compartment1);
  \endcode

  \see Reaction::addCompartment
  \see Compartment
  */
LocalizedAnnihilationProblem::LocalizedAnnihilationProblem(
        Real length, int dx, int dy,
        Real diffusionConstant,
        Real k1, Real k2, Real k3, Real k4,
        Real xminA, Real xmaxA,
        Real xminB, Real xmaxB)
    :   DiffusionModel(length, dx, dy)
{
    // two species
    Species *speciesA = addSpecies("A", diffusionConstant);
    Species *speciesB = addSpecies("B", diffusionConstant);

    // add A+A annihilation reaction
    Reaction *aplusaAnnihilation
            = addReaction("A plus A annihilation", secondOrderReactionRate(k1), speciesA, speciesA);
    std::cout << "Added reaction " << aplusaAnnihilation->id() << std::endl;
    
    // add A+B annihilation reaction
    Reaction *aplusbAnnihilation
            = addReaction("A plus B annihilation", secondOrderReactionRate(k2), speciesA, speciesB);
    std::cout << "Added reaction " << aplusbAnnihilation->id() << std::endl;

    // compute compartment cell boundaries from physical boundaries
    Real deltax = length/dx;
    int iA0 = boost::math::iround(xminA/deltax);
    int iA1 = boost::math::iround(xmaxA/deltax);
    int iB0 = boost::math::iround(xminB/deltax);
    int iB1 = boost::math::iround(xmaxB/deltax);

    // use new compartment API
    Compartment *compartment1 = new Compartment("A Creation compartment", iA0, 0, iA1, dx);
    Compartment *compartment2 = new Compartment("B Creation compartment", iB0, 0, iB1, dx);

    // add A creation reaction
    Reaction *aCreation = addReaction("A creation", zerothOrderReactionRate(k3));
    aCreation->setProductStoichiometry(speciesA, 1);
    aCreation->addCompartment(compartment1);
    std::cout << "Added reaction " << aCreation->id() << std::endl;

    // add B creation reaction
    Reaction *bCreation = addReaction("B creation", zerothOrderReactionRate(k4));
    bCreation->setProductStoichiometry(speciesB, 1);
    bCreation->addCompartment(compartment2);
    std::cout << "Added reaction " << bCreation->id() << std::endl;

    // set boundary mask
    clearBoundaryMasks();
}


} // namespace gpgmp
