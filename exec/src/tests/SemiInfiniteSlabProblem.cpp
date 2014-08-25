/*
* SemiInfiniteSlabProblem.cpp
*
*  Created on: 04/03/2010
*      Author: matthias
*/

#include "SemiInfiniteSlabProblem.h"


namespace gpgmp {


SemiInfiniteSlabProblem::SemiInfiniteSlabProblem(Real length, int dx, int dy,
                                                 int sourceNumber, Real diffusionConstant)
:   DiffusionModel(length, dx, dy)
{
    // just one species
    addSpecies("A", diffusionConstant);
    
    // set boundary mask
    cl_int4 boundaryMask = {{1, 1, 1, 1}};
    cl_int4 sourceMask = {{0, 0, sourceNumber, 0}};
    setBoundaryMasks(boundaryMask, sourceMask);
}


} // namespace gpgmp
