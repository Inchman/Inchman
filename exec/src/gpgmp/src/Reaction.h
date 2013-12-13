// TODO: add "parameter species" back in again??

/*
 * Reaction.h
 *
 *  Created on: 07/12/2009
 *      Author: matthias
 */

#ifndef __gpgmp_Reaction_h__
#define __gpgmp_Reaction_h__


#include <Base.h>

#include <map>
#include <list>
#include <ostream>

namespace gpgmp {


class Compartment;
class Reaction;
class Species;


/**
 * Class for all reactions.
 */
class Reaction : public Base {

public:
	/**
	 * Standard constructor.
	 * @param id ID of the reaction.
	 */
	Reaction(const std::string &id);

    /*!
     * Convenience constructor
     *
     * Sets the kineticLaw() and reactantStoichiometryMap() as needed, as specified by the \param rate and given reactant
     * species. The \param rate is given in units \f$\mathrm{s}^{-1}\f$.
     */
    Reaction(const std::string &id, Real rate,
              Species *reactant1=0, Species *reactant2=0,
              Species *reactant3=0, Species *reactant4=0,
              Species *reactant5=0, Species *reactant6=0,
              Species *reactant7=0, Species *reactant8=0);
    
    /**
     * Generates the kinetic law of the given reaction. If a custom law was set via setCustomLaw()
     * this will be returned verbatim. If not, a mass action kinetic law is generated from
     * the reactant stoichiometry.
     *
     * @return The kinetic law in C style ready to be inserted into the kernel template.
     */
    std::string kineticLaw() const;

    /**
     * Generates the kinetic law of the given reaction for deterministic simulations.
     * If a custom law was set via setCustomLaw()
     * this will be returned verbatim. If not, a mass action kinetic law is generated from
     * the reactant stoichiometry.
     *
     * @return The kinetic law in C style ready to be inserted into the kernel template.
     */
    std::string deterministicLaw() const;

    /**
     * Can be used to set custom (i.e. not mass action) kinetic laws. An example would be
     * Michaelis-Menten kinetics. Custom kinetic laws can make full use of the parameters
     * defined in the DiffusionModel.
     *
     * @param formula The kinetic law in C format
     */
    void setCustomLaw(const std::string &formula);

    /*!
        Generates the kinetic law from the reactant stoichiometry. Non-linear kinetic laws
        will be handled assuming integer copy counts, e.g. A*(A-1) for second-order reactions.
       */
    void generateKineticLaw();

    /*!
        Generates the kinetic law from the reactant stoichiometry for deterministic models.
        Non-linear kinetic laws
        will be handled assuming real values, e.g. A*A for second-order reactions.
       */
    void generateDeterministicLaw();

    /**
     * @return The stoichiometry map of the product species.
     */
    const std::map<Species *, StoichiometryEntry> productStoichiometryMap() const;

    /**
     * Returns the product stoichiometry for the given species.
     *
     * @param species The species to be queried
     * @return The stoichiometry
     */
    StoichiometryEntry productStoichiometry(Species *species) const;

    /**
     * Increases the product stoichiometry for the given \param species by \param stoichiometry.
     * If the species
     * is not yet in the stoichiometry map, it will be added.
     *
     * @param species The species
     * @param stoichiometry The amount by which the product stoichiometry map will be increased
     */
    void increaseProductStoichiometry(Species *species, int stoichiometry=1);

    /**
     * Sets the product stoichiometry of \param species to \param stoichiometry.
     */
    void setProductStoichiometry(Species *species, int stoichiometry);

    /**
     * Removes the \param species from the product stoichiometry map.
     */
    void unsetProductStoichiometry(Species *species);

    /**
     * @return The stoichiometry map of the reactant species.
     */
    const std::map<Species *, StoichiometryEntry> reactantStoichiometryMap() const;

    /**
     * Returns the reactant stoichiometry for the given species.
     *
     * @param species The species to be queried
     * @return The stoichiometry
     */
    StoichiometryEntry reactantStoichiometry(Species *species) const;

    /**
     * Increases the reactant stoichiometry for the given \param species by \param stoichiometry.
     * If the species
     * is not yet in the stoichiometry map, it will be added.
     *
     * @param species The species
     * @param stoichiometry The amount by which the reactant stoichiometry map will be increased
     */
    void increaseReactantStoichiometry(Species *species, int stoichiometry=1);

    /**
     * Sets the reactant stoichiometry of \param species to \param stoichiometry.
     */
    void setReactantStoichiometry(Species *species, int stoichiometry);

    /**
     * Removes the \param species from the reactant stoichiometry map.
     */
    void unsetReactantStoichiometry(Species *species);

	/**
	 * Returns a list of all compartments this reaction is localized to.
	 * Adding a compartment automatically sets this reaction to localized,
	 * i.e. isLocalized() will be set to true.
	 *
	 * @return Local compartments
	 */
	inline const std::list<Compartment *> &compartments() const;

	/**
	 * Adds a compartment where this reaction is enabled in.
	 *
	 * @param compartment Compartment for localization
	 */
	void addCompartment(Compartment *compartment);
    // TODO: add method to remove a compartment (and possibly invalidate m_isLocalized)
    
    /**
	 * Returns the reaction mask \f$w(x,y,z)\f$. This checks isLocalized() and
     * either returns 1, regardless of the argument
	 * (if the reaction is not localized) or checks its compartment list and
	 * returns 1, if the coordinates are inside one of its compartments or 0
	 * otherwise.
	 *
	 * @param x,y,z Grid coordinates to compute the reaction mask.
	 *
	 * @return Reaction mask \f$w(x,y,z)\f$.
	 *
	 * \sa RateLaw
	 */
	Real reactionMask(int x, int y, int z) const;


    /*!
      Checks if the reaction is localized.
      \return true, if the reaction is localizded
      */
    inline bool isLocalized() { return m_isLocalized; }

    // Allow a stream operator
    friend std::ostream& operator<<(std::ostream& os, const gpgmp::Reaction& reaction);

private:
    std::string generateKineticLaw() const;
    std::string generateDeterministicLaw() const;
    
	Real m_reactionRate;
    std::string m_customLaw;
    std::map<Species *, StoichiometryEntry> m_productStoichiometryMap; ///< Stoichiometry of the products
    std::map<Species *, StoichiometryEntry> m_reactantStoichiometryMap; ///< Stoichiometry of the reactants
    bool m_isLocalized; ///< Flag whether the reaction is localized or not
	std::list<Compartment *> m_compartments; ///< List of compartments this reaction is localized to
};

inline const std::list<Compartment *> & Reaction::compartments() const {
    return m_compartments;
}


} // namespace gpgmp


#endif // !__gpgmp_Reaction_h__
