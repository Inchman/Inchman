#include "RandomDriftProblem.h"
#include "Species.h"
#include "Compartment.h"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/math/special_functions/erf.hpp>

#include <sstream>
#include <cmath>

namespace gpgmp {

RandomDriftProblem::RandomDriftProblem(Real length, int dx, int dy,
                                       Real diffX, Real diffY,
                                       Real mux, Real sigmax, Real muy, Real sigmay,
                                       uint numMolecules,
                                       uint discreteK)
    : DiffusionModel(length, dx, dy)
{
    // add two species
    Species *speciesA = addSpecies("A");

    // set boundary mask to be empty/clear
    clearBoundaryMasks();

    // add World compartment
    addCompartment("World", 0, 0, dx-1, dy-1);

    // add Source compartment
    int x0 = dx/2;
    int y0 = dy/2;
    Compartment *source = addCompartment("Source", x0, y0, x0, y0);
    source->setInitialAmount(speciesA, numMolecules, gpgmp::RandomDistribution);

    // add random individual drifts
    std::map<std::string, std::vector<Real> > &properties = speciesA->getIndividualProperties();

    std::vector<Real> irx;
    std::vector<Real> iry;
    std::vector<Real> idiffx;
    std::vector<Real> idiffy;

    // init RNG
    // todo: seed with time
    boost::random::mt19937 rng;
    rng.seed();

    // distribution of drift vectors
    boost::random::normal_distribution<> ndx(mux, sigmax);
    boost::random::normal_distribution<> ndy(muy, sigmay);

    for (unsigned int i=0; i<numMolecules; i++) {
        idiffx.push_back(diffX);
        idiffy.push_back(diffY);
        irx.push_back(ndx(rng));
        iry.push_back(ndy(rng));
    }

    // and add to property list
    properties["rx"] = irx;
    properties["ry"] = iry;
    properties["diffx"] = idiffx;
    properties["diffy"] = idiffy;

    speciesA->setHasIndividualProperties(true);

    // now we add a couple of "discretized" species
    std::ostringstream ss; // to generate the compute function
    // need to undefine first .. todo: obviously this needs to be changed!
    ss <<"#undef DiffusivityX\n";
    ss <<"#undef DiffusivityY\n";
    ss <<"#undef DriftX\n";
    ss <<"#undef DriftY\n";

    uint tot = 0.;

    for (uint kx=0; kx<discreteK; kx++) {
        for (uint ky=0; ky<discreteK; ky++) {
            const float bl=-10., br=10.;
            const float db=(br-bl)/float(discreteK);

            std::ostringstream snames;
            snames << "B"<<kx<<ky;

            // add species
            Species *speciesBk = addSpecies(snames.str());

            const float bincx = bl+db/2.+kx*db;
            const float bincy = bl+db/2.+ky*db;
            const float normkx = (-boost::math::erf((bl - mux)/(std::sqrt(2.)*sigmax)) +
                                  boost::math::erf((-bl + br + br*discreteK - discreteK*mux)/(std::sqrt(2)*discreteK*sigmax)))/2.;
            const float normky = (-boost::math::erf((bl - muy)/(std::sqrt(2.)*sigmay)) +
                                  boost::math::erf((-bl + br + br*discreteK - discreteK*muy)/(std::sqrt(2)*discreteK*sigmay)))/2.;

            const uint ia = numMolecules/normkx/normky *
                    1./2. *(boost::math::erf((bincx+db/2.-mux)/(std::sqrt(2.)*sigmax)) - boost::math::erf((bincx-db/2.-mux)/(std::sqrt(2.)*sigmax))) *
                    1./2. *(boost::math::erf((bincy+db/2.-muy)/(std::sqrt(2.)*sigmay)) - boost::math::erf((bincy-db/2.-muy)/(std::sqrt(2.)*sigmay)));
            source->setInitialAmount(speciesBk, ia, gpgmp::RandomDistribution);
            // std::cout << "Creating species "<<snames.str()<<" with initial amount "<<ia<<".\n";
            tot+=ia;

            // and set diffusivity/drift
            ss << snames.str() <<".DiffusivityX = "<< diffX <<";\n";
            ss << snames.str() <<".DiffusivityY = "<< diffY <<";\n";
            ss << snames.str() <<".DriftX = " << bincx <<";\n";
            ss << snames.str() <<".DriftY = " << bincy <<";\n";
        }
    }
    std::cout <<"Finished creating species. Total = "<<tot<<".\n";

    // set homogeneous diffusivity and drift
    // std::cout <<"setting method :\n"<<ss.str()<<"\n"<<std::flush;
    setComputeDriftDiffusivityMethod(ss.str());

} // constructor

} // namespace
