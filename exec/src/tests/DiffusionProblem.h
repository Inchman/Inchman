/*
* SimpleDiffusion.h
*
*  Created on: 05/03/2010
*      Author: matthias
*/

#ifndef CSIMPLEDIFFUSION_H_
#define CSIMPLEDIFFUSION_H_


#include "DiffusionModel.h"

namespace gpgmp {


/**
  Simple diffusion test species with one species and no reactions.
  Initially, all particles
  are contained in the source compartment which is located in the center
  of the integration domain. Optionally, a drift field can be prescribed.

  This problem can be accessed in the test suite via

  \verbatim
     ./run-tests.py --homogeneous-diffusion
  \endverbatim

  The same problem with a constant drift can be run with

  \verbatim
    ./run-tests.py --homogeneous-drift
  \endverbatim
 */
class DiffusionProblem : public DiffusionModel {

public:
    /**
     * Constructor.
     *
     * @param length Length of the simulation domain [\f$\mu\mathrm{m}\f$].
     * @param dx, dy Number of cells in x and y direction
     * @param diffusionConstantX Diffusion constant
     * \f$[\mu\mathrm{m}^2\,\mathrm{s}^{-1}]\f$
     * @param numMolecules Initial number of molecules.
       \param diffusionConstantY If not zero, diffusivity in \f$y\f$ direction
       \param rx, ry Drift field.
     */
    DiffusionProblem(Real length, int dx, int dy, Real diffusionConstantX, int numMolecules,
                     Real diffusionConstantY=0., Real rx=0., Real ry=0.);
};


} // namespace gpgmp


#endif /* CSIMPLEDIFFUSION_H_ */
