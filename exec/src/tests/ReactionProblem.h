/*
* FisherProblem.h
*
*  Created on: 11/03/2010
*      Author: matthias
*/

#ifndef CREACTIONPROBLEM_H_
#define CREACTIONPROBLEM_H_


#include "DiffusionModel.h"


namespace gpgmp {


/**
  Sets up a simple (non-localized) reaction problem with two creation and
  two annihilation reactions.

  \f{eqnarray*}{
    A + A \stackrel{k_1} \rightarrow \emptyset & & A+B \stackrel{k_2} \rightarrow \emptyset \\
    \emptyset \stackrel{k_3} \rightarrow A & & \emptyset \stackrel{k_4} \rightarrow B.
  \f}

  This test problem can be accessed by the test suite using

  \verbatim
    ./run-tests.py --aplusb
  \endverbatim

  \see \vigelius10
 */
class ReactionProblem : public DiffusionModel {
    
public:
    /**
     * Standard constructor.
     *
     * @param length Physical length [\f$\mu\mathrm{m}\f$]
     * @param dx,dy	Number of subvolumes in x and y direction.
     * @param diffusionConstant Diffusion constant [\f$\mu\mathrm{m}^2\,\mathrm{s}^{-1}\f$]
     * @param k1, k2 Annihilation reaction rate constants [\f$\mathrm{s}^{-1}\f$]
     * @param k3, k4 Creation reaction rate constants [\f$\mathrm{s}^{-1}\f$]
     */
    ReactionProblem(Real length, int dx, int dy,
                    Real diffusionConstant,
                    Real k1, Real k2, Real k3, Real k4);
    
    // Note: DiffusionModel takes ownership of its species, rections and compartments.
    //       Therefore, we don't need a destructor to delete them.
};


} // namespace gpgmp


#endif /* CREACTIONPROBLEM_H_ */
