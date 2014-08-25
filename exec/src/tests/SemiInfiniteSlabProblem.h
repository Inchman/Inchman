/*
* SemiInfiniteSlabProblem.h
*
*  Created on: 04/03/2010
*      Author: matthias
*/

#ifndef CSEMIINFINITESLABPROBLEM_H_
#define CSEMIINFINITESLABPROBLEM_H_


#include "DiffusionModel.h"

namespace gpgmp {

/**
  Sets up a test problem with a particle source boundary, one species and no reactions.
  The source boundary provides a steady influx of particles into the domain
  by providing a particle reservoir
  with a constant number of particles in it. An analytic solution
  for the stationary state can be constructed which allows us to determine
  the accuracy of our solver.

  The problem can be accessed in the test suit via

  \verbatim
    ./run-tests.py --infinite-slab
  \endverbatim

  \see setBoundaryMasks
  \see \vigelius10
 */
class SemiInfiniteSlabProblem : public DiffusionModel {

public:
    /**
     * Constructor. Constructs the problem with a particle source reservoir
     * at the y=0 boundary. The source boundary artificially fixes the number of
     * particles in the cells right next to the boundary.
     * @param length Length of the simulation domain [\f$\mu\mathrm{m}\f$].
     * @param dx, dy Number of cells in x and y direction
     * @param diffusionConstant Diffusion constant
     * @param sourceNumber Number of particles in the reservoir.
     */
    SemiInfiniteSlabProblem(Real length, int dx, int dy,
                            int sourceNumber,
                            Real diffusionConstant);
};

} // namespace gpgmp


#endif /* CSEMIINFINITESLABPROBLEM_H_ */
