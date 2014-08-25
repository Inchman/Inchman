#ifndef RANDOMDRIFTPROBLEM_H
#define RANDOMDRIFTPROBLEM_H

#include "DiffusionModel.h"

namespace gpgmp {


class RandomDriftProblem : public DiffusionModel
{
public:
    RandomDriftProblem(Real length, int dx, int dy,
                       Real diffX, Real diffY,
                       Real mux, Real sigmax,
                       Real muy, Real sigmay,
                       uint numMolecules, uint discreteK);
};

}
#endif // RANDOMDRIFTPROBLEM_H
