/*
* FisherProblem.cpp
*
*  Created on: 11/03/2010
*      Author: matthias
*/

#include "FisherProblem.h"

#include <boost/filesystem.hpp>
#include <boost/math/special_functions/round.hpp>

#include <iostream>
#include <cmath>


namespace gpgmp {

/**
  The initial configuration (a smoothed step function to initiate the travelling wave) is
  set up using a pre-defined initialization script init_fisher.py which is expected to be
  in the current directory.

  \code
    loadInitScript("init_fisher.py");
  \endcode

  \see loadInitScript
  */
FisherProblem::FisherProblem(Real length, int dx, int dy,
                             Real baseConcentration,
                             const std::string &initScriptsPath)
//:   DiffusionModel(length, dx, dy)
  : DiffusionModel("Fisher problem", true)
{

  // set dimensions
  setGridDims(dx,dy);
  setPhysicalLength(length);

    // one species
    Species *speciesA = addSpecies("A", 1);
    
    // compute annihilation rate from base concentration
    Real annihilationRate = secondOrderReactionRate(1./baseConcentration);
    
    // compute particle number per compartment
    int particlesPerCompartment = particleNumberFromConcentration(baseConcentration);
    
    // print out statistics
    std::cout << "Fisher problem starts with concentration of " << baseConcentration << std::endl;
    std::cout << "Subvolume size is " << subvolumeSize()
              << ", Rate "<<annihilationRate
              << "and initial particle count " << particlesPerCompartment << std::endl;

    // and add it to the script
    setParameter("particles", particlesPerCompartment);

    // add A+A -> A annihilation reaction
    Reaction *aAnnihilation = addReaction("A plus A annihilation",
                                          annihilationRate, speciesA, speciesA);
    aAnnihilation->setProductStoichiometry(speciesA, 1);

    // A->A+A reaction
    Reaction *aCreation = addReaction("A Creation", 1.0, speciesA);
    aCreation->setProductStoichiometry(speciesA, 2);
    
    // set boundary mask to reflective
    // todo: why reflective?
    clearBoundaryMasks();

    // load the init script
    loadInitScript((boost::filesystem::path(initScriptsPath) / "init_fisher.py").string());
} // constructor

} // namespace gpgmp
