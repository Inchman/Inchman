#include "CalciumModel.h"

#include <boost/filesystem.hpp>

#include <cmath>

namespace gpgmp {

/**
  We use customized rate laws to model the complicated dynamics. We
  also require the init script models/init_ca.py to set the initial
  configuration which triggers the spiral waves.

  \see Reaction::setKineticLaw
  */
CalciumModel::CalciumModel(Real length, int dx, int dy,
                           Real diffc, Real diffp,
                           Real nmax, Real kp, Real betac,
                           Real omega,
                           const std::string &initScriptsPath)
    : DiffusionModel(length, dx, dy)
{
    float subvolume = pow(length*omega/float(dx),3)*1e3;
    float na = 6.02214e23;
    float muM = 1e-6*na*subvolume;

    // parameters
    float beta = betac*muM;
    float gamma = 2.*muM;
    float kgamma = 0.1*muM;
    float b = 0.111;
    float k1 = 0.7*muM;
    float v1 = 0.889;

    float mu = 0.583;
    float k2 = 0.7*muM;
    float taun = 2.;
    float kflux = 8.1/nmax*muM;

    float mu0=0.567;
    float mu1=0.433;
    float kmu=4.*muM;

    // add species
    Species *sc = addSpecies("Ca", diffc);
    Species *sn = addSpecies("n", 0.);
    Species *sp = addSpecies("p", diffp);


    // add leak current
    Reaction *r1 = addReaction("Leak current", beta);
    r1->setProductStoichiometry(sc, 1);

    // add SERCA pump
    Reaction *r2 = addReaction("SERCA pump");
    r2->setReactantStoichiometry(sc,1);
    r2->setCustomLaw("gamma*Ca/(kgamma+Ca)");
    setParameter("gamma", gamma);
    setParameter("kgamma", kgamma);

    // add channel flux
    Reaction *r3 = addReaction("IPR channel flux");
    r3->setProductStoichiometry(sc,1);
    r3->setCustomLaw("kflux *"
                     "(b + v1 * Ca/(k1+Ca)) "
                     "* n * (mu0 + mu1 * p/(kmu + p))");
    setParameter("b", b);
    setParameter("v1", v1);
    setParameter("k1", k1);
    setParameter("mu0", mu0);
    setParameter("mu1", mu1);
    setParameter("kmu", kmu);
    setParameter("kflux", kflux);

    // creation of open channels
    Reaction *r4 = addReaction("Open channel creation");
    r4->setProductStoichiometry(sn,1);
    r4->setCustomLaw("cc * (1-Ca*(Ca-1)/"
                     "(k2*k2 + Ca*(Ca-1)))");
    setParameter("cc", nmax/taun);
    setParameter("k2", k2);

    // annihilation of open channels
    Reaction *r5 = addReaction("Open channel annihilation", 1./taun);
    r5->setReactantStoichiometry(sn,1);
    //r5->setRateLaw(1./taun);

    // IP3 decay
    Reaction *r6 = addReaction("IP3 decay", kp);
    r6->setReactantStoichiometry(sp,1);
    //r6->setRateLaw(kp);

    // load init script
    loadInitScript((boost::filesystem::path(initScriptsPath) / "init_ca2p_atri.py").string());

    // add script parameters
    setParameter("omega", omega);
    setParameter("muM", muM);
    setParameter("nmax", nmax);

}//Constructor

} // namespace gpgmp
