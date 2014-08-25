/*
* AnnihilationModel.h
*
*  Created on: 03/03/2010
*      Author: matthias
*/

#ifndef CANNIHILATIONMODEL_H_
#define CANNIHILATIONMODEL_H_


#include "DiffusionModel.h"

namespace gpgmp {


/**
  Sets up the simple 2D \f$A+B \stackrel{k_1} \rightarrow \emptyset \f$
  annihilation model. The problem consists of an initial amount of molecules
  of species \f$B\f$, located in the center of the domain, and a homogeneous
  background concentration of \f$A\f$. Upon diffusing and reacting, the \f$B\f$ molecules
  develop a concentric annihilation front. Optionally, the problem can be set up
  with a constant drift field. Both, diffusivity and drift are homogeneous and isotropic.

  The purely diffusive problem can be accessed in the test suite using

  \verbatim
    ./run-tests.py --annihilation-2d
  \endverbatim

  The same problem with a constant drift field and homogeneous, isotropic diffusivity
  is accessed via

  \verbatim
    ./run-tests.py --annihilation-2d-drift
  \endverbatim

  \see \vigelius10
  \see \vigelius12
 */
class AnnihilationProblem : public DiffusionModel {

public:
    /**
     * Sets up a purely diffusive (no drift) model with the parameters given.
     * The reaction rate per cell is given in units
     * \f$\mathrm{M}^{-1}\,\mathrm{s}^{-1}\f$, the
     * length and diffusivity as usual in \f$\mu\mathrm{m}\f$ and
       \f$\mu\mathrm{m}^2\,\mathrm{s}^{-1}\f$, respectively. If the drift field is 0, a
       homogeneous simulation is set up.

     *
     * @param length Side length of the integration domain.
     * @param dx, dy Number of subvolumes per dimension.
     * @param numMolecules Initial number of molecules of species \f$B\f$.
     * @param diffusionConstant Diffusion constant
     * @param reactionRate Reaction rate \f$k_1\f$.
       \param drift The drift field [\f$\mu m\ \mathrm{s}^{-1}\f$]
     */
    AnnihilationProblem(Real length, int dx, int dy,
                         int numMolecules,
                         Real diffusionConstant,
                         Real reactionRate,
                         Real drift=0.);

    // Note: DiffusionModel takes ownership of its species, rections and compartments.
    //       Therefore, we don't need a destructor to delete them.
};


} // namespace gpgmp


#endif /* CANNIHILATIONMODEL_H_ */
