/*
 * Reaction.cpp
 *
 *  Created on: 07/12/2009
 *      Author: matthias
 */

#include "Reaction.h"

#include "Compartment.h"
#include "Species.h"

#include <boost/array.hpp>
#include <boost/foreach.hpp>

#include <iostream>
#include <sstream>


namespace gpgmp {


/*!
 * Constructor
 */
Reaction::Reaction(const std::string &id)
:   Base(id),
//  m_reactionRate(nan),
    m_isLocalized(false)
{}


/*!
 * Convenience constructor
 *
 * Sets the kineticLaw() and reactantStoichiometryMap() as needed, as specified by the \param rate and given reactant
 * species.
 */
Reaction::Reaction(const std::string &id, Real rate,
                   Species *reactant1, Species *reactant2,
                   Species *reactant3, Species *reactant4,
                   Species *reactant5, Species *reactant6,
                   Species *reactant7, Species *reactant8)
:   Base(id),
    m_reactionRate(rate),
    m_isLocalized(false)
{
    boost::array<Species *, 8> reactants =
    {{ reactant1, reactant2, reactant3, reactant4, reactant5, reactant6, reactant7, reactant8 }};


    BOOST_FOREACH(Species *r, reactants) {
        if (!r) continue;
        increaseReactantStoichiometry(r);
    }
}
    
    
std::string Reaction::kineticLaw() const {
    return (!m_customLaw.empty()
            ? m_customLaw
            : generateKineticLaw()); // lazy generation - user will mostly likely set far more than get
}

std::string Reaction::deterministicLaw() const {
    return (!m_customLaw.empty()
            ? m_customLaw
            : generateDeterministicLaw()); // lazy generation - user will mostly likely set far more than get
}


// the stream operator
std::ostream& operator<<(std::ostream& os, const gpgmp::Reaction& reaction)
{
    os << reaction.id() <<" [(";
    // loop through reactants
    for (std::map<Species *, StoichiometryEntry>::const_iterator i = reaction.m_reactantStoichiometryMap.begin();
         i != reaction.m_reactantStoichiometryMap.end(); ++i) {
        const Species *s=(*i).first;
        int stoimax = (*i).second;
        if (i != reaction.m_reactantStoichiometryMap.begin()) os <<" + ";
        if (stoimax > 1) os << stoimax <<" ";
        os <<  s->id();
    }
    os <<") -> (";

    for (std::map<Species *, StoichiometryEntry>::const_iterator i = reaction.m_productStoichiometryMap.begin();
         i != reaction.m_productStoichiometryMap.end(); ++i) {
        const Species *s=(*i).first;
        int stoimax = (*i).second;
        if (i != reaction.m_productStoichiometryMap.begin()) os <<" + ";
        if (stoimax > 1) os << stoimax <<" ";
        os <<  s->id();
    }
    os <<"), rate Law: " <<reaction.kineticLaw()<<", in compartments:(";
    BOOST_FOREACH(Compartment *comp, reaction.m_compartments) {
        os << comp->id() << ", ";
    }
    os <<")]";

    return os;
}


/*!
 * The kinetic law, expressed as a textual formula.
 * For example: "a * sin(t)"
 */
void Reaction::setCustomLaw(const std::string &formula) {
    m_customLaw = formula;
}


/*!
    Generates the kinetic law from the reactant stoichiometry. Non-linear kinetic laws
    will be handled assuming integer copy counts, e.g. A*(A-1) etc.. If a deterministic
    non-linear law (A*A) is desired, it needs to be set manually.
   */
std::string Reaction::generateKineticLaw() const {
    // rate comes first
    std::stringstream ss;
    ss << m_reactionRate;
    
    // loop through reactants
    for (std::map<Species *, StoichiometryEntry>::const_iterator i = m_reactantStoichiometryMap.begin();
         i != m_reactantStoichiometryMap.end(); ++i) {
        const Species *s=(*i).first;
        int stoimax = (*i).second;
        ss << " * " <<s->id();
        for (int stoi=1; stoi<stoimax; stoi++) {
            ss << " * ("<< s->id() <<"-"<<stoi<<")";
        }
    }

    std::cout <<"Kinetic law for reaction " << id() <<" is :" <<ss.str()<<"\n"<<std::flush;

    return ss.str();
}

/*!
    Generates the kinetic law from the reactant stoichiometry. Non-linear kinetic laws
    will be handled assuming integer copy counts, e.g. A*(A-1) etc.. If a deterministic
    non-linear law (A*A) is desired, it needs to be set manually.
   */
std::string Reaction::generateDeterministicLaw() const {
    // rate comes first
    std::stringstream ss;
    ss << m_reactionRate;

    // loop through reactants
    for (std::map<Species *, StoichiometryEntry>::const_iterator i = m_reactantStoichiometryMap.begin();
         i != m_reactantStoichiometryMap.end(); ++i) {
        const Species *s=(*i).first;
        int stoimax = (*i).second;
        ss << " * " << s->id();
        for (int stoi=1; stoi<stoimax; stoi++) {
            ss << " * "<< s->id();
        }
    }

    return ss.str();
} // generateDeterministicLaw

const std::map<Species *, StoichiometryEntry> Reaction::productStoichiometryMap() const
{ return m_productStoichiometryMap; }

StoichiometryEntry Reaction::productStoichiometry(Species *species) const {
    std::map<Species *, StoichiometryEntry>::const_iterator i = m_productStoichiometryMap.find(species);
    return (i != m_productStoichiometryMap.end()
            ? i->second
            : 0);
}

void Reaction::increaseProductStoichiometry(Species *species, int stoichiometry) {
    m_productStoichiometryMap[species] += stoichiometry;
}

void Reaction::setProductStoichiometry(Species *species, int stoichiometry) {
    m_productStoichiometryMap[species] = stoichiometry;
}

void Reaction::unsetProductStoichiometry(Species *species) {
    m_productStoichiometryMap.erase(species);
}

// we'll have reactant stoichiometry too
const std::map<Species *, StoichiometryEntry> Reaction::reactantStoichiometryMap() const
{ return m_reactantStoichiometryMap; }

StoichiometryEntry Reaction::reactantStoichiometry(Species *species) const {
    std::map<Species *, StoichiometryEntry>::const_iterator i = m_reactantStoichiometryMap.find(species);
    return (i != m_reactantStoichiometryMap.end()
            ? i->second
            : 0);
}

void Reaction::increaseReactantStoichiometry(Species *species, int stoichiometry) {
    m_reactantStoichiometryMap[species] += stoichiometry;
}

void Reaction::setReactantStoichiometry(Species *species, int stoichiometry) {
    m_reactantStoichiometryMap[species] = stoichiometry;
}

void Reaction::unsetReactantStoichiometry(Species *species) {
    m_reactantStoichiometryMap.erase(species);
}

void Reaction::addCompartment(Compartment *compartment) {
    m_isLocalized = true;
    m_compartments.push_back(compartment);
}


// standard implementation of the reaction mask
Real Reaction::reactionMask(int x, int y, int z) const {
    if (m_isLocalized) {
        // if the reaction is localized, iterates through all
        // compartments and returns 1, if the coordinates are
        // contained in at least one of them
        for (std::list<Compartment *>::const_iterator i=m_compartments.begin();
             i != m_compartments.end();
             ++i)
        {
            if ((*i)->isInCompartment(x, y, z)) {
                return 1.f;
            }
        }
        // and zero otherwise
        return 0.f;
    } else {
        // if it's a global reaction, we return 1
        return 1.f;
    }
}


} // namespace gpgmp
