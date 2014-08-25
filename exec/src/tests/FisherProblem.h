/*
* FisherProblem.h
*
*  Created on: 11/03/2010
*      Author: matthias
*/

#ifndef CFISHERPROBLEM_H_
#define CFISHERPROBLEM_H_


#include "DiffusionModel.h"

namespace gpgmp {


/**
 * This class is an example implementation of the Fisher problem. The Fisher problem essentially
 consists of a reversible self-annihilation reaction for a species \f$A\f$

 \f[
    A + A \stackrel{k_0} \rightleftharpoons A.
 \f]

 \f$A\f$ is allowed to diffuse and if the diffusivity is chosen correctly, the system exhibits
 travelling wave solutions. Since an analytic expression for the corresponding deterministic
 PDE can be constructed, the Fisher problem is an excellent test problem to determine the
 accuracy of any stochastic solver.

 The Fisher problem can be accessed in the test suit via

 \verbatim
   ./run-tests.py --fisher
 \endverbatim

 \see J Murray (2007): Mathematical Biology. Springer.
 \see \vigelius10
 */
class FisherProblem : public DiffusionModel {
    
public:
    /**
     * Initializes the Fisher problem.
     *
     * @param length Length of the integration domain.
     * @param dx, dy Number of subvolumes for each dimension.
     * @param baseConcentration Base concentration \f$u_0\f$ [\f$M\f$]
     */
    FisherProblem(Real length, int dx, int dy,
                  Real baseConcentration,
                  const std::string &initScriptsPath);
    
};


} // namespace gpgmp


#endif /* CFISHERPROBLEM_H_ */
