#include "OregonatorModel.h"

#include <boost/filesystem.hpp>

#include <cmath>

namespace gpgmp {

/**
  All parameters to trigger a spiral wave are hard coded. The initial
  configuration is generated in tests/init_oregonator.py
  */
OregonatorModel::OregonatorModel(Real length, int dx, int dy,
                                 Real omega, const std::string &initScriptsPath)
:   DiffusionModel(length, dx, dy)
{
    Real diffx, diffy, diffz;

    diffx = 1.5e-5f;
    diffy= 0.f;
    diffz = 0.6f*diffx;

    Species *x = addSpecies("X", diffx);
    Species *y = addSpecies("Y", diffy);
    Species *z = addSpecies("Z", diffz);
    Species *ys = addSpecies("YS", diffy);
    Species *yss = addSpecies("YSS", diffy);
    Species *r = addSpecies("R", 0);

    // parameters - will want to compute them more general
    Real k1, k2, k3, k4, k5, k6, k7;
    k1 = 843.2f;
    k2 = 1.04472e-9f/pow(omega,3);
    k3 = 4.216f;
    k4 = 5.2236e-15/pow(omega,3);
    k5 = 0.048f;
    //k6 = 4.353e-6*1e7f/pow(omega,3f);
    k6 = 100.f;
    k7 = k6;

    // and inject them for the script
    setParameter("omega", omega);
    setParameter("k1", k1);
    setParameter("k2", k2);
    setParameter("k3", k3);
    setParameter("k4", k4);
    setParameter("k5", k5);
    setParameter("k6", k6);
    setParameter("k7", k7);

    // add reactions
    // add reaction Y -> X + R
    Reaction *reaction1 = addReaction("reaction 1", k1, y);
    reaction1->setProductStoichiometry(x, 1);
    reaction1->setProductStoichiometry(r, 1);

    // add reaction X + Y -> 2 R
    Reaction *reaction2 = addReaction("reaction 2", k2, x, y);
    reaction2->setProductStoichiometry(r, 2);

   // add reaction X -> 2 X + 2 Z
    Reaction *reaction3 = addReaction("reaction 3", k3, x);
    reaction3->setProductStoichiometry(x, 2);
    reaction3->setProductStoichiometry(z, 2);

    // add reaction 2 X -> R
    Reaction *reaction4 = addReaction("reaction 4", k4, x, x);
    reaction4->setProductStoichiometry(r, 1);

    // add reaction Z -> YS + YSS
    Reaction *reaction5 = addReaction("reaction 5", k5, z);
    reaction5->setProductStoichiometry(ys, 1);
    reaction5->setProductStoichiometry(yss, 1);

    // add reaction 2 YS -> Y
    Reaction *reaction6 = addReaction("reaction 6", k6, ys, ys);
    reaction6->setProductStoichiometry(y, 1);

    // add reaction 2 YSS -> YS
    Reaction *reaction7 = addReaction("reaction 7", k7, yss, yss);
    reaction7->setProductStoichiometry(ys, 1);

    loadInitScript((boost::filesystem::path(initScriptsPath) / "init_oregonator.py").string());

}
} // namespace gpgmp
