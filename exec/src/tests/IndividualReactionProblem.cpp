#include "IndividualReactionProblem.h"
#include "Compartment.h"
#include "Species.h"

namespace gpgmp {

IndividualReactionProblem::IndividualReactionProblem(Real length, int dx, int dy,
                                                     Real diffX, Real diffY,
                                                     Real k1, Real k2, Real k3, Real k4)

    :   DiffusionModel("AB Annihilation", true)
{
  std::cout <<"Constructing.."<<std::flush;
  setPhysicalLength(length);
  setGridDims(dx, dy);

    // set boundary mask to be empty/clear
    clearBoundaryMasks();

    // add two species for each individual/non-individual (for comparison)
    Species *speciesA = addSpecies("A");
    Species *speciesB = addSpecies("B");

    // add diffusivity
    std::map<std::string, std::vector<Real> > &propertiesA = speciesA->getIndividualProperties();
    std::map<std::string, std::vector<Real> > &propertiesB = speciesB->getIndividualProperties();

    std::vector<Real> irxA;
    std::vector<Real> iryA;
    std::vector<Real> idiffxA;
    std::vector<Real> idiffyA;
    std::vector<Real> irxB;
    std::vector<Real> iryB;
    std::vector<Real> idiffxB;
    std::vector<Real> idiffyB;

    // we fill in values for all possible species.. since atm we can't specify a rule
    // to intialize the attributes of new particles
    for (unsigned int i=0; i<262144; i++) {
        idiffxA.push_back(diffX);
        idiffyA.push_back(diffY);
        irxA.push_back(0.);
        iryA.push_back(0.);
        idiffxB.push_back(diffX);
        idiffyB.push_back(diffY);
        irxB.push_back(0.);
        iryB.push_back(0.);
    }

    // and add to property list
    propertiesA["rx"] = irxA;
    propertiesA["ry"] = iryA;
    propertiesA["diffx"] = idiffxA;
    propertiesA["diffy"] = idiffyA;
    propertiesB["rx"] = irxB;
    propertiesB["ry"] = iryB;
    propertiesB["diffx"] = idiffxB;
    propertiesB["diffy"] = idiffyB;

    speciesA->setHasIndividualProperties(true);
    speciesB->setHasIndividualProperties(true);

    // add reactions for individuals

    // add A+A -> 0 annihilation reaction
    Reaction *aplusaAnnihilation = addReaction("A plus A annihilation",
                                               k1, speciesA, speciesA);
    std::cout << "Added reaction "<<aplusaAnnihilation->id() << std::endl;

    // add A creation reaction
    Reaction *aCreation = addReaction("A creation", k3);
    aCreation-> setProductStoichiometry(speciesA, 1);
    std::cout << "Added reaction " << aCreation->id() << std::endl;

    // add A+B -> 0 annihilation reaction
    Reaction *aplusbAnnihilation = addReaction("A plus B annihilation",
                                               k2, speciesA, speciesB);
    std::cout << "Added reaction "<<aplusbAnnihilation->id() << std::endl;

    // add B creation reaction
    Reaction *bCreation = addReaction("B creation", k4);
    bCreation-> setProductStoichiometry(speciesB, 1);
    std::cout << "Added reaction " << bCreation->id() << std::endl;

    // collective species (for comparison)
    Species *speciesAc = addSpecies("Ac");
    Species *speciesBc = addSpecies("Bc");

    // add reactions for collective species

    // add A+A -> 0 annihilation reaction
    Reaction *aplusaAnnihilationc = addReaction("A plus A annihilation (collective)",
                                               k1, speciesAc, speciesAc);
    std::cout << "Added reaction "<<aplusaAnnihilationc->id() << std::endl;

    // add A+B -> 0 annihilation reaction
    Reaction *aplusbAnnihilationc = addReaction("A plus B annihilation (collective)",
                                               k2, speciesAc, speciesBc);
    std::cout << "Added reaction "<<aplusbAnnihilationc->id() << std::endl;

    // add A creation reaction
    Reaction *aCreationc = addReaction("A creation (collective)", k3);
    aCreationc-> setProductStoichiometry(speciesAc, 1);
    std::cout << "Added reaction " << aCreationc->id() << std::endl;

    // add B creation reaction
    Reaction *bCreationc = addReaction("B creation (collective)", k4);
    bCreationc-> setProductStoichiometry(speciesBc, 1);
    std::cout << "Added reaction " << bCreationc->id() << std::endl;

    // add creation reactions
    //Real k1 = 0.00305;
    //Real k1 = 3.05; // todo: change me back. But don't know where that number comes from..
    //Reaction *creationA = addReaction("A creation", k1);
    //creationA->setProductStoichiometry(speciesA, 2);
    //Reaction *creationB = addReaction("B creation", k1);
    //creationB->setProductStoichiometry(speciesB, 2);

    // set diffusivity/drift
    setParameter("diffX", diffX);
    setParameter("diffY", diffY);
    setParameter("driftX", 0.);
    setParameter("driftY", 0.);

    // set homogeneous diffusivity and drift
    setComputeDriftDiffusivityMethod(gpgmp::DT_HOMOGENEOUS);
}

} // namespace gpgmp
