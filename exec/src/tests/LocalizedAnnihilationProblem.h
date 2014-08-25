/*
* LocalizedAnnhilationProblem.h
*
*  Created on: 09/03/2010
*      Author: matthias
*/

#ifndef CLOCALIZEDANNHILATIONPROBLEM_H_
#define CLOCALIZEDANNHILATIONPROBLEM_H_

#include "DiffusionModel.h"


namespace gpgmp {


/**
 Sets up a localized annihilation problem. This is essentially the same problem as AnnihilationProblem,
 except that we here allow diffusion and also locally constrain the creation reactions to compartments.

 This problem can be accessed in the test suite using

 \verbatim
    ./run-tests.py --localised-ab
 \endverbatim

 \see AnnihilationProblem
 \see \vigelius10
 */
class LocalizedAnnihilationProblem : public DiffusionModel {
    
public:
    /**
     * Constructor. Creates the underlying diffusion model and creates the
     * compartments to localize the reactions.
     *
     * @param length Side length of integration area.
     * @param dx, dy Number of subvolumes in each dimension.
     * @param diffusionConstant Diffusivity.
     * @param k1, k2 Rate constants for annihilation reactions [in \f$\mathrm{M}^{-1}\,\mathrm{s}^{-1}\f$].
     * @param k3, k4 Rate constants for creation reactions [in \f$\mathrm{M}\,\mathrm{s}^{-1}\f$].
     * @param xminA, xmaxA Boundaries for the \f$A\f$ creation compartment.
     * @param xminB, xmaxB Boundaries for the \f$B\f$ creation compartment.
     */
    LocalizedAnnihilationProblem(Real length, int dx, int dy,
                                  Real diffusionConstant,
                                  Real k1, Real k2, Real k3, Real k4,
                                  Real xminA, Real xmaxA, Real xminB, Real xmaxB);
};
    

} // namespace gpgmp


#endif /* CLOCALIZEDANNHILATIONPROBLEM_H_ */
