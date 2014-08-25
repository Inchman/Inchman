#include "MajorityVoteProblem.h"

#include "Species.h"
#include "Compartment.h"

#include <sstream>
#include <cstdlib>
#include <ctime>
#include <string>

#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>

#include <stdlib.h>
#include <time.h>
#include <cmath>

namespace gpgmp {

bool MajorityVoteProblem::isMarginal(uint pos, uint n)
{
    bool ret=false;
    for (uint i=0; i<n; i++) {
        if (pos == (-2 - i - i*i + 2 * n + 2 * i * n)/2) {
            ret=true;
        }
    }

    return ret;
}

// gets the next red+1 position
unsigned int MajorityVoteProblem::getNextRedPosition(uint pos, uint n)
{
    // we first need to find the i position
    uint i=0;
    while (pos >= (2*n-i)*(1+i)/2) i++;

    // we get then the red+1 position by adding n-i
    return pos+(n-i);
}

/*
MajorityVoteProblem::MajorityVoteProblem(Real length, int dx, int dy, uint numMolecules, uint nReds, uint numEncounters, Real avoidanceRadius)
    : DiffusionModel(length, dx, dy), m_individual(true)
{
    // we only have one species
    Species *speciesA = addSpecies("A");

    // set boundary mask to be empty/clear
    clearBoundaryMasks();

    // add World compartment
    Compartment *world = addCompartment("World", 0, 0, dx-1, dy-1);
    world->setInitialAmount(speciesA, numMolecules, gpgmp::RandomDistribution);

    // we will have no diffusion and a random drift vector
    Real diff = 0.;
    Real speed = 0.01;

    // property vectors
    std::vector<Real> irx;
    std::vector<Real> iry;
    std::vector<Real> idiffx;
    std::vector<Real> idiffy;
    std::vector<Real> icolor;
    std::vector<Real> inumreds;
    std::vector<Real> inumgreens;

    // initialize RNG
    srand(time(NULL));

    // fill in initial property vectors
    for (unsigned int i=0; i<numMolecules; i++) {
        idiffx.push_back(diff);
        idiffy.push_back(diff);

        // choose random drift
        double pi = std::atan(1.)*4.;
        double th =((double)rand())/((double)RAND_MAX)*2*pi;

        irx.push_back(speed*cos(th));
        iry.push_back(speed*sin(th));

        // make random experiment to decide which starting color (0 = green, 1 = red)
        double r = rand();
        if (r>=0.5)
            icolor.push_back(1);
        else
            icolor.push_back(0);

        inumreds.push_back(0);
        inumgreens.push_back(0);
    }

    // and add to property list
    std::map<std::string, std::vector<Real> > &properties = speciesA->getIndividualProperties();
    properties["rx"] = irx;
    properties["ry"] = iry;
    properties["diffx"] = idiffx;
    properties["diffy"] = idiffy;
    properties["color"] = icolor;
    properties["nreds"] = inumreds;
    properties["ngreens"] = inumgreens;

    // and allow individual properties
    speciesA->setHasIndividualProperties(true);

    // we only have one reaction that handles it all
    // todo: change the reaction prob..
    Reaction *r = addReaction ("individual reaction", 0.485675, speciesA, speciesA);
    r->setProductStoichiometry(speciesA, 2);

    // reactions do not change number of individuals
    this->setHasConstantIndividuals(true);

    // individuals should behave ballistically
    this->setHasBallisticBoundaryConditions(true);
}

*/

MajorityVoteProblem::MajorityVoteProblem(Real length, int dx, int dy,
                                         Real diffX, Real diffY,
                                         Real mux, Real muy,
                                         uint numMolecules, uint nReds,
                                         uint numEncounters,
                                         Real reactionProb,const std::string &initScriptsPath)
    : DiffusionModel(length, dx, dy), m_individual(false), m_isNonDiffusive(false)
{
    // this list will hold all red and green species
    // we need (i^2+i)/2 species per red and green
    // where i is the number of encounters after which the species "votes"
    uint maxSpecies = (numEncounters*numEncounters+numEncounters)/2;
    std::vector<gpgmp::Species *> speciesListRed(maxSpecies);
    std::vector<gpgmp::Species *> speciesListGreen(maxSpecies);

    // add species and reactions corresponding to scheme:
    // Rxy - red, has encountered x red and y green robots
    // Gxy - green, has encountered x red and y green robots
    uint c=0;
    for (uint i=0; i<numEncounters; i++) {
        for (uint j=0; j<numEncounters-i; j++) {
            // add species
            std::ostringstream snamered;
            snamered << "R"<<i<<j;
            speciesListRed[c] = addSpecies(snamered.str(), diffX);
            std::ostringstream snamegreen;
            snamegreen << "G"<<i<<j;
            speciesListGreen[c] = addSpecies(snamegreen.str(),diffX);

            // print out
            std::cout <<"Added species "<<snamered.str()<<" and "<<snamegreen.str()<<".\n";
            // increase counter;
            c++;
        }
    } // create species

    // set diffusivity and drift for species
    setParameter("diffX", diffX);
    setParameter("diffY", diffY);
    setParameter("driftX", mux);
    setParameter("driftY", muy);

    // set homogeneous diffusivity and drift
    setComputeDriftDiffusivityMethod(gpgmp::DT_HOMOGENEOUS);

    // add reactions
    // Gxy+Rwz -> Gx+1y + Rwz+1
    uint cg = 0;
    for (uint x=0; x<numEncounters; x++) {
        for (uint y=0; y<numEncounters-x; y++) {
            // to name the reaction
            std::ostringstream snamers;
            std::ostringstream snamerd;

            // find source and target species for green
            Species *Gxy = speciesListGreen[cg];
            snamers << "G"<<x<<y<<"("<<Gxy->id()<<")";

            Species *Txp1y;
            if ((x+y)==(numEncounters-1)) {
                // target species is either G00 or R00
                if (x>=y) {
                    // turns red
                    Txp1y = speciesListRed[0];
                    snamerd <<"R00"<<"("<<Txp1y->id()<<")";
                } else {
                    // stays greenstd::map<std::string, std::vector<Real> > &propertiesA = speciesA->getIndividualProperties();
                    Txp1y = speciesListGreen[0];
                    snamerd <<"G00"<<"("<<Txp1y->id()<<")";
                }
            } else {
                // target species is Gx+1y
                Txp1y = speciesListGreen[cg+(numEncounters-x)];
                snamerd <<"G"<<x+1<<y<<"("<<Txp1y->id()<<")";
            }

            // now pair with all red robots
            uint cr=0;
            for (uint w=0; w<numEncounters; w++) {
                for (uint z=0; z<numEncounters-w; z++) {
                    // find source and target species for red
                    Species *Rxy = speciesListRed[cr];
                    std::ostringstream sss;
                    std::ostringstream ssd;
                    sss << "R"<<w<<z << "("<<Rxy->id()<<")";

                    Species *Twzp1;
                    if ((w+z)==(numEncounters-1)) {
                        // target species is either G00 or R00
                        if (w>z) {
                            // turns red
                            Twzp1 = speciesListRed[0];
                            ssd <<"R00" << "("<<Twzp1->id()<<")";
                        } else {
                            // stays green
                            Twzp1 = speciesListGreen[0];
                            ssd <<"G00" << "("<<Twzp1->id()<<")";
                        }
                    } else {
                        // target species is Rw+1z
                        Twzp1 = speciesListRed[cr+1];
                        ssd <<"R"<<w<<z+1<< "("<<Twzp1->id()<<")";
                    }

                    std::ostringstream rname;
                    rname <<snamers.str()<<" + "<<sss.str()<<" -> "<<snamerd.str()<<" + " << ssd.str();
                    std::cout <<"Adding reaction "<<rname.str()<<".\n";
                    Reaction *r = addReaction (rname.str(), reactionProb, Gxy, Rxy);
                    r->setProductStoichiometry(Txp1y, 1);
                    if (Txp1y == Twzp1)
                        r->setProductStoichiometry(Twzp1, 2);
                    else
                        r->setProductStoichiometry(Twzp1, 1);

                    // increase counter
                    cr++;
                } // pair with red robots z
            } // pair with red robots w

            // increase counter
            cg++;
        } // green robots y
    } // green robots x

    // Gxy+Gwz -> Gxy+1 + Gwz+1
    cg = 0;
    for (uint x=0; x<numEncounters; x++) {
        for (uint y=0; y<numEncounters-x; y++) {
            // to name the reaction
            std::ostringstream snamers;
            std::ostringstream snamerd;

            // find source and target species for green
            Species *Gxy = speciesListGreen[cg];
            snamers << "G"<<x<<y<<"("<<Gxy->id()<<")";

            Species *Txyp1;
            if ((x+y)==(numEncounters-1)) {
                // target species is either G00 or R00
                if (x>y) {
                    // turns red
                    Txyp1 = speciesListRed[0];
                    snamerd <<"R00"<<"("<<Txyp1->id()<<")";
                } else {
                    // stays green
                    Txyp1 = speciesListGreen[0];
                    snamerd <<"G00"<<"("<<Txyp1->id()<<")";
                }
            } else {
                // target species is Gxyp1
                Txyp1 = speciesListGreen[cg+1];
                snamerd <<"G"<<x<<y+1<<"("<<Txyp1->id()<<")";
            }

            // now pair with all other green robots
            uint cgg=0;
            for (uint w=0; w<numEncounters; w++) {
                for (uint z=0; z<numEncounters-w; z++) {
                    // find source and target species for red
                    Species *GGxy = speciesListGreen[cgg];
                    std::ostringstream sss;
                    std::ostringstream ssd;
                    sss << "G"<<w<<z<<"("<<GGxy->id()<<")";

                    Species *Twzp1;
                    if ((w+z)==(numEncounters-1)) {
                        // target species is either G00 or R00
                        if (w>z) {
                            // turns red
                            Twzp1 = speciesListRed[0];
                            ssd <<"R00"<<"("<<Twzp1->id()<<")";
                        } else {
                            // stays green
                            Twzp1 = speciesListGreen[0];
                            ssd <<"G00"<<"("<<Twzp1->id()<<")";
                        }
                    } else {
                        // target species is Gwzp1
                        Twzp1 = speciesListGreen[cgg+1];
                        ssd <<"G"<<w<<z+1<<"("<<Twzp1->id()<<")";
                    }

                    std::stringstream rname;
                    rname << snamers.str()<<" + "<<sss.str()<<" -> "<<snamerd.str();
                    std::cout <<"Adding reaction "<< rname.str() <<" + " << ssd.str()<<".\n";

                    Reaction *r = addReaction(rname.str(), reactionProb, Gxy, GGxy);
                    r->setProductStoichiometry(Txyp1, 1);
                    if (Twzp1 == Txyp1)
                        r->setProductStoichiometry(Twzp1, 2);
                    else
                        r->setProductStoichiometry(Twzp1, 1);

                    // increase counter
                    cgg++;
                } // pair with other green robots z
            } // pair with other green robots w

            // increase counter
            cg++;
        } // green robots y
    } // green robots x

    // Rxy+Rwz -> Rxy+1 + Rwz+1
    uint cr = 0;
    for (uint x=0; x<numEncounters; x++) {
        for (uint y=0; y<numEncounters-x; y++) {
            // to name the reaction
            std::ostringstream snamers;
            std::ostringstream snamerd;

            // find source and target species for green
            Species *Rxy = speciesListRed[cr];
            snamers << "R"<<x<<y<<"("<<Rxy->id()<<")";

            Species *Txp1y;
            if ((x+y)==(numEncounters-1)) {
                // target species is either G00 or R00
                if (x>=y) {
                    // turns red
                    Txp1y = speciesListRed[0];
                    snamerd <<"R00"<<"("<<Txp1y->id()<<")";
                } else {
                    // stays green
                    Txp1y = speciesListGreen[0];
                    snamerd <<"G00"<<"("<<Txp1y->id()<<")";
                }
            } else {
                // target species is Gxyp1
                Txp1y = speciesListRed[cr+(numEncounters-x)];
                snamerd <<"R"<<x+1<<y<<"("<<Txp1y->id()<<")";
            }

            // now pair with all other red robots
            uint crr=0;
            for (uint w=0; w<numEncounters; w++) {
                for (uint z=0; z<numEncounters-w; z++) {
                    // find source and target species for red
                    Species *RRxy = speciesListRed[crr];
                    std::ostringstream sss;
                    std::ostringstream ssd;
                    sss << "R"<<w<<z<<"("<<RRxy->id()<<")";

                    Species *Twp1z;
                    if ((w+z)==(numEncounters-1)) {
                        // target species is either G00 or R00
                        if (w>=z) {
                            // turns red
                            Twp1z = speciesListRed[0];
                            ssd <<"R00"<<"("<<Twp1z->id()<<")";
                        } else {
                            // stays green
                            Twp1z = speciesListGreen[0];
                            ssd <<"G00"<<"("<<Twp1z->id()<<")";
                        }
                    } else {
                        // target species is Rwzp1
                        Twp1z = speciesListRed[crr+(numEncounters-w)];
                        ssd <<"R"<<w+1<<z<<"("<<Twp1z->id()<<")";
                    }

                    std::stringstream rname;
                    rname << snamers.str()<<" + "<<sss.str()<<" -> "<<snamerd.str()<<" + " << ssd.str();
                    std::cout <<"Adding reaction "<< rname.str() <<".\n";
                    Reaction *r = addReaction(rname.str(), reactionProb, Rxy, RRxy);
                    r->setProductStoichiometry(Txp1y, 1);
                    if (Txp1y == Twp1z)
                        r->setProductStoichiometry(Twp1z, 2);
                    else
                        r->setProductStoichiometry(Twp1z, 1);

                    // increase counter
                    crr++;
                } // pair with other red robots z
            } // pair with other red robots w

            // increase counter
            cr++;
        } // green robots y
    } // green robots x

    // distribute species randomly
    /*
    Compartment *source = addCompartment("Source", 0, 0, dx-1, dy-1);
    source->setInitialAmount(speciesListRed[0], nReds, gpgmp::RandomDistribution);
    source->setInitialAmount(speciesListGreen[0], numMolecules-nReds, gpgmp::RandomDistribution);
    */

    // load init script
    loadInitScript((boost::filesystem::path(initScriptsPath) / "init_majority.py").string());

    // add script parameters
    setParameter("numMolecules", numMolecules);
    setParameter("p", static_cast<float>(nReds)/static_cast<float>(numMolecules));

    setRegenerateInitStates(true);
}

// sets up the toy model
//
MajorityVoteProblem::MajorityVoteProblem(Real length, int dx, int dy, uint numMolecules, uint nReds, Real reactionProb, const std::string &initScriptsPath)
    : DiffusionModel(length, dx, dy), m_individual(false), m_isNonDiffusive(true)
{
    // add four species R_0, R_1, G_0, G_1
    Species *R0 = addSpecies("R0", 0.);
    Species *R1 = addSpecies("R1", 0.);
    Species *G0 = addSpecies("G0", 0.);
    Species *G1 = addSpecies("G1", 0.);

    // set diffusivity and drift for species
    setParameter("diffX", 0.);
    setParameter("diffY", 0.);
    setParameter("driftX", 0.);
    setParameter("driftY", 0.);

    // set homogeneous diffusivity and drift
    setComputeDriftDiffusivityMethod(gpgmp::DT_HOMOGENEOUS);

    // add reactions
    // R0 + G0 -> R1 + G1 (i)
    Reaction *r1 = addReaction("R0 + G0 -> R1 + G1", reactionProb, R0, G0);
    r1 -> setProductStoichiometry(R1, 1);
    r1 -> setProductStoichiometry(G1, 1);

    // R0 + G1 -> R1 + R0 (ii)
    Reaction *r2 = addReaction("R0 + G1 -> R1 + R0", reactionProb, R0, G1);
    r2 -> setProductStoichiometry(R1, 1);
    r2 -> setProductStoichiometry(R0, 1);

    // R1 + G0 -> G0 + G1 (iii)
    Reaction *r3 = addReaction("R1 + G0 -> G0 + G1", reactionProb, R1, G0);
    r3 -> setProductStoichiometry(G1, 1);
    r3 -> setProductStoichiometry(G0, 1);

    // R1 + G1 -> G0 + R0 (iv)
    Reaction *r4 = addReaction("R1 + G1 -> G0 + R0", reactionProb, R1, G1);
    r4 -> setProductStoichiometry(G0, 1);
    r4 -> setProductStoichiometry(R0, 1);

    // R0 + R1 -> R0 + R0 (v)
    Reaction *r5 = addReaction("R0 + R1 -> R0 + R0", reactionProb, R0, R1);
    r5 -> setProductStoichiometry(R0, 2);

    // R1 + R1 -> R0 + R0 (vi)
    Reaction *r6 = addReaction("R1 + R1 -> R0 + R0", reactionProb, R1, R1);
    r6 -> setProductStoichiometry(R0, 2);

    // G0 + G1 -> G0 + G0 (vii)
    Reaction *r7 = addReaction("G0 + G1 -> G0 + G0", reactionProb, G0, G1);
    r7 -> setProductStoichiometry(G0, 2);

    // G1 + G1 -> G0 + G0 (viii)
    Reaction *r8 = addReaction("G1 + G1 -> G0 + G0", reactionProb, G1, G1);
    r8 -> setProductStoichiometry(G0, 2);

    // load init script
    loadInitScript((boost::filesystem::path(initScriptsPath) / "init_majority_no_diffusion.py").string());

    // add script parameters
    setParameter("numMolecules", numMolecules);
    setParameter("p", static_cast<float>(nReds)/static_cast<float>(numMolecules));

    setRegenerateInitStates(true);

}

// this one sets up dx*dy majority vote problems without diffusion
MajorityVoteProblem::MajorityVoteProblem(Real length, int dx, int dy, uint numMolecules, uint nReds, uint numEncounters, Real reactionProb, Real recognitionError, const std::string &initScriptsPath)
    : DiffusionModel(length, dx, dy), m_individual(false), m_isNonDiffusive(true)
{
    if (numEncounters < 3) {
        std::cout <<"ERROR: Need at least three encounters!\n";
        exit(1);
    }

    // this list will hold all red and green species
    // we need (i^2+i)/2 species per red and green
    // where i is the number of encounters after which the species "votes"
    uint maxSpecies = (numEncounters*numEncounters+numEncounters)/2;
    std::vector<gpgmp::Species *> speciesListRed(maxSpecies);
    std::vector<gpgmp::Species *> speciesListGreen(maxSpecies);

    // add species and reactions corresponding to scheme:
    // Rxy - red, has encountered x red and y green robots
    // Gxy - green, has encountered x red and y green robots
    uint c=0;
    for (uint i=0; i<numEncounters; i++) {
        for (uint j=0; j<numEncounters-i; j++) {
            // add species
            std::ostringstream snamered;
            snamered << "R"<<i<<j;
            speciesListRed[c] = addSpecies(snamered.str(), 0.);
            std::ostringstream snamegreen;
            snamegreen << "G"<<i<<j;
            speciesListGreen[c] = addSpecies(snamegreen.str(),0.);

            // print out
            std::cout <<"Added species "<<snamered.str()<<" and "<<snamegreen.str()<<".\n";
            // increase counter;
            c++;
        }
    } // create species

    // set diffusivity and drift for species
    setParameter("diffX", 0.);
    setParameter("diffY", 0.);
    setParameter("driftX", 0.);
    setParameter("driftY", 0.);

    // set homogeneous diffusivity and drift
    setComputeDriftDiffusivityMethod(gpgmp::DT_HOMOGENEOUS);

    // these are the actual reaction probabilities for the correct reaction and the error reaction
    double kcorrect = reactionProb*(1.-2.*recognitionError-recognitionError*recognitionError);
    double kwrong1  = reactionProb*recognitionError;
    double kwrongb  = reactionProb*recognitionError*recognitionError;

    // add reactions
    // Gxy+Rwz -> Gx+1y + Rwz+1
    // Gxy+Rwz -> Gxy+1 + Rw+1z (for recognition error)
    uint cg = 0;
    uint totnreacts=0.; // just debug counter to count number of reactions
    for (uint x=0; x<numEncounters; x++) {
        for (uint y=0; y<numEncounters-x; y++) {
            // to name the reaction
            std::ostringstream snamers;
            std::ostringstream snamerd;
            std::ostringstream snamerdw;

            // find source and target species for green
            Species *Gxy = speciesListGreen[cg];
            snamers << "G"<<x<<y<<"("<<Gxy->id()<<")";

            Species *Txp1y;
            Species *Txyp1;
            if ((x+y)==(numEncounters-1)) {
                // target species is either G00 or R00
                if (x>=y) {
                    // turns red
                    Txp1y = speciesListRed[0];
                    snamerd <<"R00"<<"("<<Txp1y->id()<<")";

                    // wrong decision
                    if (x>y) {
                        Txyp1 = speciesListRed[0];
                        snamerdw <<"R00"<<"("<<Txyp1->id()<<")";
                    }
                    else {
                        Txyp1 = speciesListGreen[0];
                        snamerdw <<"G00"<<"("<<Txyp1->id()<<")";
                    }
                } else {
                    // stays green
                    Txp1y = speciesListGreen[0];
                    snamerd <<"G00"<<"("<<Txp1y->id()<<")";
                    // and the wrong decision stays green too
                    Txyp1 = speciesListGreen[0];
                    snamerdw <<"G00"<<"("<<Txyp1->id()<<")";
                }
            } else {
                // target species is Gx+1y
                Txp1y = speciesListGreen[cg+(numEncounters-x)];
                snamerd <<"G"<<x+1<<y<<"("<<Txp1y->id()<<")";
                // or Gxy+1
                Txyp1 = speciesListGreen[cg+1];
                snamerdw <<"G"<<x<<y+1<<"("<<Txyp1->id()<<")";
            }

            // now pair with all red robots
            uint cr=0;
            for (uint w=0; w<numEncounters; w++) {
                for (uint z=0; z<numEncounters-w; z++) {
                    // find source and target species for red
                    Species *Rxy = speciesListRed[cr];
                    std::ostringstream sss;
                    std::ostringstream ssd;
                    std::ostringstream ssdw;
                    sss << "R"<<w<<z << "("<<Rxy->id()<<")";

                    Species *Twzp1;
                    Species *Twp1z;
                    if ((w+z)==(numEncounters-1)) {
                        // target species is either G00 or R00
                        if (w>z) {
                            // stays red
                            Twzp1 = speciesListRed[0];
                            ssd <<"R00" << "("<<Twzp1->id()<<")";
                            // and the wrong decision .. stays red too
                            Twp1z = speciesListRed[0];
                            ssdw <<"R00" << "("<<Twp1z->id()<<")";
                        } else {
                            // stays green
                            Twzp1 = speciesListGreen[0];
                            ssd <<"G00" << "("<<Twzp1->id()<<")";
                            if (w==z) {
                                Twp1z = speciesListRed[0];
                                ssdw <<"R00" << "("<<Twp1z->id()<<")";
                            } else {
                                Twp1z = speciesListGreen[0];
                                ssdw <<"G00" << "("<<Twp1z->id()<<")";
                            }
                        }
                    } else {
                        // target species is Rw+1z
                        Twzp1 = speciesListRed[cr+1];
                        ssd <<"R"<<w<<z+1<< "("<<Twzp1->id()<<")";
                        // and Rwzp1
                        Twp1z = speciesListRed[cr+(numEncounters-w)];
                        ssdw <<"R"<<w+1<<z<< "("<<Twp1z->id()<<")";
                    }

                    std::ostringstream rname, rnamew;
                    rname <<snamers.str()<<" + "<<sss.str()<<" -> "<<snamerd.str()<<" + " << ssd.str();
                    std::cout <<"Adding reaction "<<rname.str()<<".\n";
                    totnreacts++;
                    Reaction *r = addReaction (rname.str(), kcorrect, Gxy, Rxy);
                    r->setProductStoichiometry(Txp1y, 1);
                    if (Txp1y == Twzp1)
                        r->setProductStoichiometry(Twzp1, 2);
                    else
                        r->setProductStoichiometry(Twzp1, 1);

                    // and add wrong reaction if needed
                    if (recognitionError > 0) {
                        // both wrong
                        rnamew <<snamers.str()<<" + "<<sss.str()<<" -> "<<snamerdw.str()<<" + " << ssdw.str();
                        std::cout <<"Adding errenous reaction "<<rnamew.str()<<".\n"<<std::flush;
                        totnreacts++;
                        r = addReaction(rnamew.str(), kwrongb, Gxy, Rxy);
                        r->setProductStoichiometry(Txyp1, 1);
                        if (Txyp1 == Twp1z)
                            r->setProductStoichiometry(Twp1z, 2);
                        else
                            r->setProductStoichiometry(Twp1z, 1);

                        std::ostringstream rnamew1, rnamew2;
                        // first wrong
                        rnamew1 <<snamers.str()<<" + "<<sss.str()<<" -> "<<snamerdw.str()<<" + " << ssd.str();
                        std::cout <<"Adding errenous reaction (first wrong) "<<rnamew1.str()<<".\n"<<std::flush;
                        totnreacts++;
                        r = addReaction(rnamew.str(), kwrong1, Gxy, Rxy);
                        r->setProductStoichiometry(Txyp1, 1);
                        if (Txyp1 == Twzp1)
                            r->setProductStoichiometry(Twzp1, 2);
                        else
                            r->setProductStoichiometry(Twzp1, 1);

                        // second wrong
                        rnamew2 <<snamers.str()<<" + "<<sss.str()<<" -> "<<snamerd.str()<<" + " << ssdw.str();
                        std::cout <<"Adding errenous reaction (second wrong) "<<rnamew2.str()<<".\n"<<std::flush;
                        totnreacts++;
                        r = addReaction(rnamew2.str(), kwrong1, Gxy, Rxy);
                        r->setProductStoichiometry(Txp1y, 1);
                        if (Txp1y == Twp1z)
                            r->setProductStoichiometry(Twp1z, 2);
                        else
                            r->setProductStoichiometry(Twp1z, 1);
                    }
                    // increase counter
                    cr++;
                } // pair with red robots z
            } // pair with red robots w

            // increase counter
            cg++;
        } // green robots y
    } // green robots x

    // Gxy+Gwz -> Gxy+1 + Gwz+1
    for (uint x=0; x<maxSpecies; x++) {
        for (uint y=x; y<maxSpecies; y++) {
            // get source species
            Species *G1 = speciesListGreen[x];
            Species *G2 = speciesListGreen[y];

            // get target species 1
            Species *T1;
            Species *T1wrong;

            // is it at marginal position?
            if (isMarginal(x, numEncounters)) {
                // check if over changeover limit (see notebook index_computations.nb)
                if (x > getChangeoverLimit(x, numEncounters)) {
                    // change to red
                    T1 = speciesListRed[0];
                    T1wrong = speciesListRed[0];
                } else {
                    T1 = speciesListGreen[0];

                    // wrong species
                    if (x==getChangeoverLimit(x, numEncounters))
                        T1wrong = speciesListRed[0];
                    else
                        T1wrong = speciesListGreen[0];
                }
            } else {
                T1 = speciesListGreen[x+1];
                T1wrong = speciesListGreen[getNextRedPosition(x, numEncounters)];
            }

            // get target species 2
            Species *T2;
            Species *T2wrong;

            // is it at marginal position?
            if (isMarginal(y, numEncounters)) {
                // check if over changeover limit (see notebook index_computations.nb)
                if (y > getChangeoverLimit(y, numEncounters)) {
                    // change to red
                    T2 = speciesListRed[0];
                    T2wrong = speciesListRed[0];
                } else {
                    T2 = speciesListGreen[0];

                    // wrong species
                    if (y==getChangeoverLimit(y, numEncounters))
                        T2wrong = speciesListRed[0];
                    else
                        T2wrong = speciesListGreen[0];
                }


            } else {
                T2 = speciesListGreen[y+1];
                T2wrong = speciesListGreen[getNextRedPosition(y, numEncounters)];
            }

            // print reaction
            std::stringstream rname;
            rname <<G1->id()<<" + "<<G2->id()<<" -> "<< T1->id()<<" + "<<T2->id();
            std::cout <<"Adding reaction "<<rname.str() <<"\n"<<std::flush;

            // and add it
            totnreacts++;
            Reaction *r = addReaction(rname.str(), kcorrect, G1, G2);
            r->setProductStoichiometry(T1, 1);
            if (T1 == T2)
                r->setProductStoichiometry(T2, 2);
            else
                r->setProductStoichiometry(T2, 1);

            // wrong reactions
            if (recognitionError > 0) {
                // first wrong
                std::stringstream rnamew1;
                rnamew1 << G1->id() << " + " << G2->id() <<" -> "<< T1wrong->id() << " + "<<T2->id();
                std::cout <<"Adding erreneous (first wrong) reaction "<<rnamew1.str()<<"\n"<<std::flush;
                totnreacts++;
                r = addReaction(rnamew1.str(), kwrong1, G1, G2);
                r->setProductStoichiometry(T1wrong, 1);
                if (T1wrong == T2)
                    r->setProductStoichiometry(T2, 2);
                else
                    r->setProductStoichiometry(T2, 1);

                // second wrong
                std::stringstream rnamew2;
                rnamew2 << G1->id() << " + " << G2->id() <<" -> "<< T1->id() << " + "<<T2wrong->id();
                std::cout <<"Adding erreneous (second wrong) reaction "<<rnamew2.str()<<"\n"<<std::flush;
                totnreacts++;
                r = addReaction(rnamew2.str(), kwrong1, G1, G2);
                r->setProductStoichiometry(T1, 1);
                if (T1== T2wrong)
                    r->setProductStoichiometry(T2wrong, 2);
                else
                    r->setProductStoichiometry(T2wrong, 1);

                // both wrong
                std::stringstream rnamewb;
                rnamewb << G1->id() << " + " << G2->id() <<" -> "<< T1wrong->id() << " + "<<T2wrong->id();
                std::cout <<"Adding erreneous (both wrong) reaction "<<rnamewb.str()<<"\n"<<std::flush;
                totnreacts++;
                r = addReaction(rnamewb.str(), kwrongb, G1, G2);
                r->setProductStoichiometry(T1wrong, 1);
                if (T1wrong == T2wrong)
                    r->setProductStoichiometry(T2wrong, 2);
                else
                    r->setProductStoichiometry(T2wrong, 1);

            }
        }
    }

    // Rxy+Rwz -> Rxy+1 + Rwz+1
    for (uint x=0; x<maxSpecies; x++) {
        for (uint y=x; y<maxSpecies; y++) {
            // get source species
            Species *R1 = speciesListRed[x];
            Species *R2 = speciesListRed[y];

            // get target species 1
            Species *T1;
            Species *T1wrong;

            // is it at marginal position?
            if (isMarginal(x, numEncounters)) {
                // check if over changeover limit (see notebook index_computations.nb)
                if (x >= getChangeoverLimit(x, numEncounters)) {
                    // change to red
                    T1 = speciesListRed[0];

                    // wrong species
                    if (x == getChangeoverLimit(x, numEncounters))
                        T1wrong = speciesListGreen[0];
                    else
                        T1wrong = speciesListRed[0];
                } else {
                    T1 = speciesListGreen[0];
                    T1wrong = speciesListGreen[0];
                }
            } else {
                T1 = speciesListRed[getNextRedPosition(x, numEncounters)];
                T1wrong = speciesListRed[x+1];
            }

            // get target species 2
            Species *T2;
            Species *T2wrong;

            // is it at marginal position?
            if (isMarginal(y, numEncounters)) {
                // check if over changeover limit (see notebook index_computations.nb)
                if (y >= getChangeoverLimit(y, numEncounters)) {
                    // change to red
                    T2 = speciesListRed[0];

                    // wrong species
                    if (y == getChangeoverLimit(y, numEncounters))
                        T2wrong = speciesListGreen[0];
                    else
                        T2wrong = speciesListRed[0];
                } else {
                    T2 = speciesListGreen[0];
                    T2wrong = speciesListGreen[0];
                }
            } else {
                T2 = speciesListRed[getNextRedPosition(y, numEncounters)];
                T2wrong = speciesListRed[y+1];
            }

            // print reaction
            std::stringstream rname;
            rname <<R1->id()<<" + "<<R2->id()<<" -> "<< T1->id()<<" + "<<T2->id();
            std::cout <<"Adding reaction "<<rname.str() <<"\n"<<std::flush;

            // and add it
            totnreacts++;
            Reaction *r = addReaction(rname.str(), kcorrect, R1, R2);
            r->setProductStoichiometry(T1, 1);
            if (T1 == T2)
                r->setProductStoichiometry(T2, 2);
            else
                r->setProductStoichiometry(T2, 1);

            // wrong reactions
            if (recognitionError > 0) {
                // first wrong
                std::stringstream rnamew1;
                rnamew1 << R1->id() <<" + " << R2->id() <<" -> " << T1wrong->id() << " + " << T2->id();
                std::cout <<"Adding erreneous (first wrong) reaction "<<rnamew1.str() << "\n" << std::flush;
                totnreacts++;
                r = addReaction(rnamew1.str(), kwrong1, R1, R2);
                r->setProductStoichiometry(T1wrong, 1);
                if (T1wrong == T2)
                    r->setProductStoichiometry(T2, 2);
                else
                    r->setProductStoichiometry(T2, 1);

                // second wrong
                std::stringstream rnamew2;
                rnamew2 << R1->id() <<" + " << R2->id() <<" -> " << T1->id() << " + " << T2wrong->id();
                std::cout <<"Adding erreneous (second wrong) reaction "<<rnamew2.str() << "\n" << std::flush;
                totnreacts++;
                r = addReaction(rnamew2.str(), kwrong1, R1, R2);
                r->setProductStoichiometry(T1, 1);
                if (T1== T2wrong)
                    r->setProductStoichiometry(T2wrong, 2);
                else
                    r->setProductStoichiometry(T2wrong, 1);

                // both wrong
                std::stringstream rnamewb;
                rnamewb << R1->id() <<" + " << R2->id() <<" -> " << T1wrong->id() << " + " << T2wrong->id();
                std::cout <<"Adding erreneous (both wrong) reaction "<<rnamewb.str() << "\n" << std::flush;
                totnreacts++;
                r = addReaction(rnamewb.str(), kwrongb, R1, R2);
                r->setProductStoichiometry(T1wrong, 1);
                if (T1wrong == T2wrong)
                    r->setProductStoichiometry(T2wrong, 2);
                else
                    r->setProductStoichiometry(T2wrong, 1);

            }
        }
    }
    std::cout <<"Added "<<totnreacts<<" reactions.\n";
    std::cout <<"Reaction propensities: no error = "<<kcorrect<<", one wrong = "<< kwrong1 <<" both wrong = "<<kwrongb
             <<"Sum = "<<(kcorrect + 2. * kwrong1 + kwrongb)<<"\n"<<std::flush;

    /*
    // Gxy+Gwz -> Gxy+1 + Gwz+1
    cg = 0;
    totnreacts=0;
    for (uint x=0; x<numEncounters; x++) {
        for (uint y=0; y<numEncounters-x; y++) {
            // to name the reaction
            std::ostringstream snamers;
            std::ostringstream snamerd;

            // find source and target species for green
            Species *Gxy = speciesListGreen[cg];
            snamers << "G"<<x<<y<<"("<<Gxy->id()<<")";

            Species *Txyp1;
            if ((x+y)==(numEncounters-1)) {
                // target species is either G00 or R00
                if (x>y) {
                    // turns red
                    Txyp1 = speciesListRed[0];
                    snamerd <<"R00"<<"("<<Txyp1->id()<<")";
                } else {
                    // stays green
                    Txyp1 = speciesListGreen[0];
                    snamerd <<"G00"<<"("<<Txyp1->id()<<")";
                }
            } else {
                // target species is Gxyp1
                Txyp1 = speciesListGreen[cg+1];
                snamerd <<"G"<<x<<y+1<<"("<<Txyp1->id()<<")";
            }

            // now pair with all other green robots
            for (uint w=x; w<numEncounters; w++) {
                for (uint z=x; z<numEncounters-w; z++) {
                    // compute index to second source species
                    // TODODODODO: you'll need to continue here ..
                    // come up with a function to compute the speciesindex from w and z
                    // and then assign cgg to it
                    // and then compare to the notebook with nencounters=3 (mesoscopic_reactions.nb)
                    // and then you need to adapt the HDF5 writer too
                    // and then just do a couple of runs with this problem and see how it compares to spatial
                    uint cgg=getSpeciesIndexFrom2D(w, z, numEncounters);

                    // find source and target species for red
                    Species *GGxy = speciesListGreen[cgg];
                    std::ostringstream sss;
                    std::ostringstream ssd;
                    sss << "G"<<w<<z<<"("<<GGxy->id()<<")";

                    Species *Twzp1;
                    if ((w+z)==(numEncounters-1)) {
                        // target species is either G00 or R00
                        if (w>z) {
                            // turns red
                            Twzp1 = speciesListRed[0];
                            ssd <<"R00"<<"("<<Twzp1->id()<<")";
                        } else {
                            // stays green
                            Twzp1 = speciesListGreen[0];
                            ssd <<"G00"<<"("<<Twzp1->id()<<")";
                        }
                    } else {
                        // target species is Gwzp1
                        Twzp1 = speciesListGreen[cgg+1];
                        ssd <<"G"<<w<<z+1<<"("<<Twzp1->id()<<")";
                    }

                    std::stringstream rname;
                    rname << snamers.str()<<" + "<<sss.str()<<" -> "<<snamerd.str();
                    std::cout <<"Adding reaction "<< rname.str() <<" + " << ssd.str()<<".\n"<<std::flush;
                    totnreacts++;

                    Reaction *r = addReaction(rname.str(), reactionProb, Gxy, GGxy);
                    r->setProductStoichiometry(Txyp1, 1);
                    if (Twzp1 == Txyp1)
                        r->setProductStoichiometry(Twzp1, 2);
                    else
                        r->setProductStoichiometry(Twzp1, 1);

                    // increase counter
                    cgg++;
                } // pair with other green robots z
            } // pair with other green robots w

            // increase counter
            cg++;
        } // green robots y
    } // green robots x

    std::cout <<"Added a total of "<<totnreacts<<" reactions in the green self-reaction part.\n";
    exit(1);

    // Rxy+Rwz -> Rxy+1 + Rwz+1
    uint cr = 0;
    for (uint x=0; x<numEncounters; x++) {
        for (uint y=0; y<numEncounters-x; y++) {
            // to name the reaction
            std::ostringstream snamers;
            std::ostringstream snamerd;

            // find source and target species for green
            Species *Rxy = speciesListRed[cr];
            snamers << "R"<<x<<y<<"("<<Rxy->id()<<")";

            Species *Txp1y;
            if ((x+y)==(numEncounters-1)) {
                // target species is either G00 or R00
                if (x>=y) {
                    // turns red
                    Txp1y = speciesListRed[0];
                    snamerd <<"R00"<<"("<<Txp1y->id()<<")";
                } else {
                    // stays green
                    Txp1y = speciesListGreen[0];
                    snamerd <<"G00"<<"("<<Txp1y->id()<<")";
                }
            } else {
                // target species is Gxyp1
                Txp1y = speciesListRed[cr+(numEncounters-x)];
                snamerd <<"R"<<x+1<<y<<"("<<Txp1y->id()<<")";
            }

            // now pair with all other red robots
            uint crr=0;
            for (uint w=0; w<numEncounters; w++) {
                for (uint z=0; z<numEncounters-w; z++) {
                    // find source and target species for red
                    Species *RRxy = speciesListRed[crr];
                    std::ostringstream sss;
                    std::ostringstream ssd;
                    sss << "R"<<w<<z<<"("<<RRxy->id()<<")";

                    Species *Twp1z;
                    if ((w+z)==(numEncounters-1)) {
                        // target species is either G00 or R00
                        if (w>=z) {
                            // turns red
                            Twp1z = speciesListRed[0];
                            ssd <<"R00"<<"("<<Twp1z->id()<<")";
                        } else {
                            // stays green
                            Twp1z = speciesListGreen[0];
                            ssd <<"G00"<<"("<<Twp1z->id()<<")";
                        }
                    } else {
                        // target species is Rwzp1
                        Twp1z = speciesListRed[crr+(numEncounters-w)];
                        ssd <<"R"<<w+1<<z<<"("<<Twp1z->id()<<")";
                    }

                    std::stringstream rname;
                    rname << snamers.str()<<" + "<<sss.str()<<" -> "<<snamerd.str()<<" + " << ssd.str();
                    std::cout <<"Adding reaction "<< rname.str() <<".\n";
                    Reaction *r = addReaction(rname.str(), reactionProb, Rxy, RRxy);
                    r->setProductStoichiometry(Txp1y, 1);
                    if (Txp1y == Twp1z)
                        r->setProductStoichiometry(Twp1z, 2);
                    else
                        r->setProductStoichiometry(Twp1z, 1);

                    // increase counter
                    crr++;
                } // pair with other red robots z
            } // pair with other red robots w

            // increase counter
            cr++;
        } // green robots y
    } // green robots x
    */
    // distribute species randomly
    /*
    Compartment *source = addCompartment("Source", 0, 0, dx-1, dy-1);
    source->setInitialAmount(speciesListRed[0], nReds, gpgmp::RandomDistribution);
    source->setInitialAmount(speciesListGreen[0], numMolecules-nReds, gpgmp::RandomDistribution);
    */

    // load init script
    loadInitScript((boost::filesystem::path(initScriptsPath) / "init_majority_no_diffusion.py").string());

    // add script parameters
    setParameter("numMolecules", numMolecules);
    setParameter("p", static_cast<float>(nReds)/static_cast<float>(numMolecules));

    setRegenerateInitStates(true);

}


void MajorityVoteProblem::writeSingleTimeHDF5(hid_t currentDataGroup, const char *buffer, size_t bufferStep, hid_t bufferType)
{
    if (m_individual) {
        // just call the standard implementation
        DiffusionModel::writeSingleTimeHDF5(currentDataGroup, buffer, bufferStep, bufferType);
    } if (m_isNonDiffusive)  {
      std::cout <<"MajorityVoteProblem::writeSingleTimeHDF5 called for non-diffusive system.\n";

        // we assume that the buffer is int ..
        const int *iBuffer = &(runState().front());

        // we need two new buffers
        int *iReds = new int[gridArea()];
        int *iGreens = new int[gridArea()];

        int scount = 0; // running counter for species

        // initialize buffers
        for (int x=0; x<gridWidth(); x++) {
            for (int y=0; y<gridHeight(); y++) {
                iReds[x+y*gridWidth()]=0;
                iGreens[x+y*gridWidth()]=0;
            }
        }

        for (int x=0; x<gridWidth(); x++) {
            for (int y=0; y<gridHeight(); y++)
            {
                scount=0;
                BOOST_FOREACH(Species *s, species())
                {
                    // is it red or green?
                    if ((s->id().c_str())[0]=='R') {
                        //std::cout <<"Species "<<s->id().c_str()<<" is RED. Found "<<count <<" particles.\n";
                        iReds[x+y*gridWidth()] += iBuffer[x + y*gridWidth() + scount*gridArea()];
                    } else if ((s->id().c_str())[0]=='G') {
                        //std::cout <<"Species "<<s->id().c_str()<<" is GREEN. Found "<<count <<" particles.\n";
                        iGreens[x+y*gridWidth()] += iBuffer[x + y*gridWidth() + scount*gridArea()];
                    } else {
                        //std::cout <<"Species "<<s->id().c_str()<<" is UNKNOWN. Found "<<count <<" particles.\n";
                    }

                    // increase species index
                    scount++;
                }
		// DEBUG: check if total species count agrees for each cell
		//if ((iReds[x+y*gridWidth()] + iGreens[x+y*gridWidth()]) != 150) {
		//  std::cerr <<"ERROR: at x="<<x<<" and y = "<<y<<" we find r="<<iReds[x+y*gridWidth()]<<" and g="<<iGreens[x+y*gridWidth()]<<".\n"<<std::flush;
		//  exit(1);
		//
            }
        }

        // dimensions of the dataset
        std::cout <<"Grid dimensions are: "<<gridWidth()<<" x "<<gridHeight()<<".\n";
        hsize_t dims[2];
        dims[0] = gridWidth();
        dims[1] = gridHeight();

        // data set and data space
        hid_t dataspace, dataset;

        // status variable
        herr_t status;

        dataspace = H5Screate_simple(2, dims, 0);
        dataset = H5Dcreate2(currentDataGroup, "Red", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, iReds);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);

        dataspace = H5Screate_simple(2, dims, 0);
        dataset = H5Dcreate2(currentDataGroup, "Green", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, iGreens);
        status = H5Dclose(dataset);
        status = H5Sclose(dataspace);

        // delete temp buffers
        delete iReds;
        delete iGreens;
    } else {
        std::cout <<"MajorityVoteProblem::writeSingleTimeHDF5 called.\n";



        // count numbers
        int nreds=0;
        int ngreens=0;

        // we assume that the buffer is int ..
        const int *iBuffer = &(runState().front());

        int scount = 0; // running counter for species

        BOOST_FOREACH(Species *s, species())
        {
            // count all particles for this species
            int count=0;
            for (int x=0; x<gridWidth(); x++)
                for (int y=0; y<gridHeight(); y++)
                    count += iBuffer[x + y*gridWidth() + scount*gridArea()];

            // is it red or green?
            if ((s->id().c_str())[0]=='R') {
                //std::cout <<"Species "<<s->id().c_str()<<" is RED. Found "<<count <<" particles.\n";
                nreds += count;
            } else if ((s->id().c_str())[0]=='G') {
                //std::cout <<"Species "<<s->id().c_str()<<" is GREEN. Found "<<count <<" particles.\n";
                ngreens += count;
            } else {
                //std::cout <<"Species "<<s->id().c_str()<<" is UNKNOWN. Found "<<count <<" particles.\n";
            }

            // increase species index
            scount++;
        }

        //std::cout <<"Found "<<nreds<<" reds and "<< ngreens <<" greens. Total is "<<(nreds+ngreens)<<".\n"<<std::flush;

        // Save into two data groups

        // dimensions of the dataset
        hsize_t dims[1];
        dims[0] = 1;

        // data set and data space
        hid_t dataspace, dataset;

        // status variable
        herr_t status;

        dataspace = H5Screate_simple(1, dims, 0);
        dataset = H5Acreate2(currentDataGroup, "Red", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Awrite(dataset, H5T_NATIVE_INT, &nreds);
        status = H5Aclose(dataset);
        status = H5Sclose(dataspace);

        dataspace = H5Screate_simple(1, dims, 0);
        dataset = H5Acreate2(currentDataGroup, "Green", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Awrite(dataset, H5T_NATIVE_INT, &ngreens);
        status = H5Aclose(dataset);
        status = H5Sclose(dataspace);
    }
} // writeSingleTimeHDF5

} // namespace gpgmp
