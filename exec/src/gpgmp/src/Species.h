/*
 * Species.h
 *
 *  Created on: 08/09/2010
 *      Author: aidan
 */

#ifndef __gpgmp_Species_h__
#define __gpgmp_Species_h__


#include <Base.h>
#include <map>
#include <vector>
#include <ostream>


namespace gpgmp {
    
/**
 * This class represents an Inchman species.
 */
class Species : public Base {
    
public:
    /**
     * Constructs a new species object with a given \param id and diffusivity \param diffusionConstant.
     * The diffusivity assigned here is mainly used for the homogeneous solver. In the
     * inhomogeneous solvers, the user needs to compute the diffusivity and drift manually. However,
     * this value is available in the inhomogeneous computeDriftDiffusivity method as a parameter for
     * the species structure.
     *
     * @param id ID of the species.
     * @param diffusionConstant Diffusion constant in
     * \f$\mu\mathrm{m}^2\,\mathrm{s}^{-1}\f$.
     */
    Species(const std::string &id, Real diffusionConstant = 0.0);
    ~Species();

    /**
     * @return Diffusivity of the species.
     */
    inline Real diffusionConstant() const;

    /**
     * Sets the diffusivity of this species to \param diffusionConstant.
     */*
    void setDiffusionConstant(Real diffusionConstant);


    /** @name Individual species
     * The species can be marked as an individual species. Individual species are
     * handled by the IndividualSolver. Individual species currently only have the
     * properties diffusivity in both directions ("diffx" and "diffy") and drift in
     * both directions ("rx" and "ry").
     */
    ///@{
    /**
     * Marks the species as an individual species if \param hasIP is set to true.
     */
    void setHasIndividualProperties(bool hasIP);
    /**
     * @return True, if the species is an individual species.
     */
    bool hasIndividualProperties() const;
    /**
     * @return A list of all individual properties of the species.
     */
    std::map<std::string, std::vector<Real> > &getIndividualProperties();
    /**
     * Initializes the individuals to a particular starting value \param positions.
     */
    void setIndividualPositions(std::vector<int> positions);
    /**
     * @return The current positions of the individual species. This host array will
     * only be updated on output events.
     */
    std::vector<int> getIndividualPositions() const;
    ///@}

    // Allow a stream operator
    friend std::ostream& operator<<(std::ostream& os, const gpgmp::Species& species)
    {
        os << species.id() <<" [Diffusivity: "<<species.diffusionConstant()<<", individual: "<<species.hasIndividualProperties()<<"]";
        return os;
    }
private:
    Real m_diffusionConstant; ///< Diffusivity of this species (for the homogeneous solvers).
    bool m_hasIndividualProperties; ///< True if the species has individual-based properties associated with it

    std::map<std::string, std::vector<Real> > m_individualProperties; ///< Map with individual properties' name and initial values
    std::vector<Real> *m_rx; ///< The individual drift in x direction.
    std::vector<Real> *m_ry; ///< The individual drift in y direction.
    std::vector<Real> *m_dx; ///< The individual diffusion in x direction.
    std::vector<Real> *m_dy; ///< The individual diffusion in y direction.

    std::vector<int> m_individualPositions; ///< The individual grid positions for the species.
};
    
    
inline Real Species::diffusionConstant() const {
    return m_diffusionConstant;
}

inline void Species::setHasIndividualProperties(bool hasIP) {m_hasIndividualProperties = hasIP;}
inline bool Species::hasIndividualProperties() const {return m_hasIndividualProperties;}
inline std::map<std::string, std::vector<Real> > &Species::getIndividualProperties() {return m_individualProperties;}
inline void Species::setIndividualPositions(std::vector<int> positions) {m_individualPositions = positions;}
inline std::vector<int> Species::getIndividualPositions() const {return m_individualPositions;}
} // namespace gpgmp


#endif // !__gpgmp_Species_h__
