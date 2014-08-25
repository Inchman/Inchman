/*
 * Compartment.h
 *
 *  Created on: 28/05/2010
 *      Author: matthias
 */

// TODO: move to just using x,y,width,height ?

#ifndef __gpgmp_Compartment_h__
#define __gpgmp_Compartment_h__


#include "Base.h"

#include <hdf5.h>
#include <map>
#include <ostream>

namespace gpgmp {


class Species;
    
    
/**
 * Enumeration for initial amount specification. Valid values are "homogeneous" and "random"
 */
enum Distribution {
    HomogeneousDistribution,
    RandomDistribution,
    LinearXDistribution
//    LinearYDistribution
};
    
    
/**
 * Structure to hold initial amounts.
 */
struct InitialAmount {
    int amount; ///< initial Amount
    Distribution dist; ///< distribution of that amount in Compartment.
};


/**
 * Class to represent compartments. Compartment are rectangular (currently)
 * regions in the integration domain.
 */
class Compartment : public Base {
    
public:
	/**
	 * Creates an unspecified named compartment.
	 *
	 * @param id ID of the compartment.
	 */
    Compartment(const std::string &id);

    /**
     * Creates a compartment with the given coordinates (in device units).
     *
     * \param id ID of the compartment.
     * \param x0, y0 Upper-left corner of the compartment rectangle
     * \param x1, y1 Lower-right corner of the compartment
     */
	Compartment(const std::string &id, int x0, int y0, int x1, int y1);


	/**
	 * Sets the dimensions of the compartment.
	 *
	 * @param x0,y0,x1,y1 Dimensions.
	 */
    void setBounds(int x0, int y0, int x1, int y1);
    void getBounds(int &x0, int &y0, int &x1, int &y1) const;
    
    int x0() const { return x0_; }
    int x1() const { return x1_; }
    int y0() const { return y0_; }
    int y1() const { return y1_; }
    void setX0(int x0) { x0_ = x0; }
    void setX1(int x1) { x1_ = x1; }
    void setY0(int y0) { y0_ = y0; }
    void setY1(int y1) { y1_ = y1; }
    
    /**
	 * Returns the width of the compartment.
	 *
	 * @return Width.
	 */
	int width() const {return (x1_-x0_)+1;}
    void setWidth(int width) { x1_=x0_+width-1; }
    
	/**
	 * Returns the height of the compartment.
	 * @return Height.
	 */
	int height() const {return (y1_-y0_)+1;}
    void setHeight(int height) { y1_=y0_+height-1; }

    /**
     * Sets an initial amount for a species in that compartment.
     *
     * \param speciesId ID of the species
     * \param amount Initial amount (number of particles)
     * \param dist How the initial number of particles are distributed
     */
	void setInitialAmount(Species *species, int amount,
                          Distribution distribution=HomogeneousDistribution);

	/**
	 * Returns a constant reference to the initial amounts for the species in this
	 * compartment.
	 *
	 * \return initialAmounts The list of initial amount prescription.
	 */
	const std::map<Species *, InitialAmount> &initialAmounts() const;
    
    /**
	 * Checks whether a coordinate tuple \f$(x,y,z)\f$ is located inside the compartment boundaries \f$\Omega\f$.
	 *
	 * \param x,y,z Coordinates to check.
	 * \return 1 if \f$(x,y,z) \in \Omega\f$ and 0 otherwise.
	 */
	bool isInCompartment(int x, int y, int z) const;

	/**
	 * Saves the compartment to an hdf5 file.
	 *
	 * \param group Id of the group to be saved in.
	 */
	herr_t saveHdf5(hid_t group) const;

    // Allow a stream operator
    friend std::ostream& operator<<(std::ostream& os, const gpgmp::Compartment& compartment)
    {
        os << compartment.id() <<"[(" << compartment.x0() <<", "<< compartment.y0() << ") - (" << compartment.x1() << ", " << compartment.y1() << ")]";
        return os;
    }

private:
	int x0_, y0_, x1_, y1_; ///< Compartment boundaries.

	/**
	 * List of initial values for species.
	 */
	std::map<Species *, InitialAmount> initialAmounts_;
};
    
    
inline bool Compartment::isInCompartment(int x, int y, int /*z*/) const
{
    return (   x >= x0_
            && x <= x1_
            && y >= y0_
            && y <= y1_);
}

} // namespace gpgmp


#endif // !__gpgmp_Compartment_h__
