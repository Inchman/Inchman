#include "GrayScottModel.h"

#include <boost/filesystem.hpp>

#include <iostream>

namespace gpgmp {

/**
  The initial configuration which is used to trigger a spiral wave is provided by the init script models/init_gray_scott.py.
  */
GrayScottModel::GrayScottModel(Real length, int dx, int dy,
                               Real diffU, Real diffV,
                               Real k, Real F, Real omega,
                               const std::string &initScriptsPath)
    : DiffusionModel(length, dx, dy)
{
    Real kf = F;
    Real k2 = F+k;
    Real k1=1./(omega*omega);
    Real u0=1.*omega;

    std::cout <<"Setting up Gray-Scott problem with diffusion constant :"
            <<diffU << " "<< diffV <<"\n";
    std::cout <<"k1 = "<<k1<<", k2 = "<<k2<<", kf = "<<kf<<", u0 = "<<u0<<"\n";


    // add both species to main model
    Species *speciesU = addSpecies("U", diffU);
    Species *speciesV = addSpecies("V", diffV);

    // add U+2 V -> 3 V (autocatalytic)
    Reaction *autocatalytic = addReaction("autocatalytic", k1, speciesU, speciesV, speciesV);
    autocatalytic->setProductStoichiometry(speciesV, 3);

    // U decay
    addReaction("U decay", kf, speciesU);

    // V decay
    addReaction("V decay", k2, speciesV);

    // U creation
    Reaction *ucreation = addReaction("U creation", kf*u0);
    ucreation->setProductStoichiometry(speciesU, 1);

    // set boundary mask
    clearBoundaryMasks();

    // load init script
    loadInitScript((boost::filesystem::path(initScriptsPath) / "init_gray_scott.py").string());

    // add script parameters
    setParameter("omega", omega);
    setParameter("k", k);
    setParameter("F", F);

    // we need to provide unscaled u0 and k1
    setParameter("u0", u0/omega);
    setParameter("k1", k1*omega*omega);
} // constructor

} // namespace gpgmp
