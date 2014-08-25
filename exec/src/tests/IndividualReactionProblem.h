#ifndef INDIVIDUALREACTIONPROBLEM_H
#define INDIVIDUALREACTIONPROBLEM_H

#include "DiffusionModel.h"

namespace gpgmp {

class IndividualReactionProblem : public DiffusionModel
{
public:
    IndividualReactionProblem(Real length, int dx, int dy, Real diffX, Real diffY, Real k1, Real k2, Real k3, Real k4);
};

} // namespace gpgmp
#endif // INDIVIDUALREACTIONPROBLEM_H
