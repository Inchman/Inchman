#ifndef ORNSTEINUHLENBECKPROBLEM_H
#define ORNSTEINUHLENBECKPROBLEM_H

#include "DiffusionModel.h"

namespace gpgmp {

/** This class sets up the classical Ornstein-Uhlenbeck process. It is defined
  by the stochastic differential equation

  \f[
    \mathrm{d} \mathbf{X}_t = -\mathbf{\gamma} \mathbf{X}_t \mathrm{d}t + \mathbf{b}\, \mathrm{d}\mathbf{W}_t,
  \f]

  where \f$\mathbf{\gamma}\f$ is a constant matrix (not necessarily diagonal) and \f$\mathbf{b}=\mathrm{diag}\{2 D_x, 2 D_y\}\f$ encodes the diffusivity. Details can
  be found in Vigelius&Meyer (2012). We use one species \f$A\f$ and no reactions.

  The problem can be accessed through the test suite via

  \verbatim
    ./run-tests.py --ornstein-uhlenbeck
  \endverbatim

  \see \vigelius12
  */
class OrnsteinUhlenbeckProblem : public DiffusionModel
{
public:
    /**
      Sets up the model.

      \param length Domain size in \f$\mu\mathrm{m}\f$
      \param dx, dy Number of grid cells in each dimension
      \param diffusivityX, diffusivityY Diffusivity \f$D_x, D_y\f$ in \f$\mu\mathrm{m}^2\ \mathrm{s}^{-1}\f$
      \param gxx, gxy, gyx, gyy Matrix elements of \f$\gamma\f$ in \f$\mathrm{s}^{-1}\f$
      \param numMolecules Number of molecules
      */
    OrnsteinUhlenbeckProblem(Real length, int dx, int dy,
                             Real diffusivityX, Real diffusivityY,
                             Real gxx, Real gxy, Real gyx, Real gyy,
                             int numMolecules);
};
} // namespace gpgmp

#endif // ORNSTEINUHLENBECKPROBLEM_H
