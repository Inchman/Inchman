/*
 * Compartment.cpp
 *
 *  Created on: 28/05/2010
 *      Author: matthias
 */

#include "Compartment.h"


namespace gpgmp {


Compartment::Compartment(const std::string &id)
:   Base(id) {
    setBounds(0, 0, 0, 0);
}


Compartment::Compartment(const std::string &id, int x0, int y0, int x1, int y1)
:   Base(id) {
    setBounds(x0, y0, x1, y1);
}
 

void Compartment::setBounds(int x0, int y0, int x1, int y1) {
    x0_ = x0;
    y0_ = y0;
    x1_ = x1;
    y1_ = y1;
}
    
void Compartment::getBounds(int &x0, int &y0, int &x1, int &y1) const {
    x0 = x0_;
    y0 = y0_;
    x1 = x1_;
    y1 = y1_;
}

    
void Compartment::setInitialAmount(Species *species, int amount, Distribution dist) {
	InitialAmount init;
	init.amount = amount;
	init.dist = dist;
	initialAmounts_[species]=init;
}
    
const std::map<Species *, InitialAmount> & Compartment::initialAmounts() const {
    return initialAmounts_;
}

    
herr_t Compartment::saveHdf5(hid_t group) const
{
	hid_t dataspace, dataset;
	herr_t status;

	hsize_t dims[1];
	dims[0]=4;

	dataspace = H5Screate_simple(1, dims, NULL);
	// dataset = H5Dcreate(group, id().c_str(), H5T_NATIVE_INT, dataspace,
        //                 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	dataset = H5Dcreate2(group, id().c_str(), H5T_NATIVE_INT, dataspace,
                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
	int boundaries[4] = { x0_, y0_, x1_, y1_ };
    status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(boundaries));
	status = H5Dclose(dataset);
	status = H5Sclose(dataspace);

	return status;
}


} // namespace gpgmp
