#ifndef NONLINEARPROBLEM_H
#define NONLINEARPROBLEM_H

#include "DiffusionModel.h"

namespace gpgmp {
/**
    This class sets up a genuinely non-linear test problem.  The time evolution of this test problem is governed by the one-dimensional, nonlinear SDE
\f[
  \label{eq:nonlinear:master}
  \mathrm{d} X_t = - \left[ \omega X_t + \theta \left<x(t)\right> \right] \mathrm{d}t + \sqrt{2 D}\, \mathrm{d}W_t,
\f]
where the drift field implicitly depends on the probability distribution \f$f(x,t)\f$ through its first moment
\f[
  \left< x(t) \right> = \int_{-\infty}^\infty x f(x,t)\,\mathrm{d}x.
\f]

We have one species and no reactions.

\see \vigelius12
**/

class NonlinearProblem : public DiffusionModel
{
public:
    /**
      Sets up the model.

      \param length Size of the integration domain in \f$\mu\mathrm{m}\f$
      \param dx, dy Number of grid cells in each dimension.
      \param diffX, diffY Diffusivity in each direction (in \f$\mu\mathrm{m}^2\ \mathrm{s}^{-1}\f$)
      \param theta Parameter \f$\theta\f$ (in \f$\mathrm{s}^{-1}\f$)
      \param omega Parameter \f$\omega\f$ (in \f$\mathrm{s}\f$)
      */
    NonlinearProblem(Real length, int dx, int dy,
                     Real diffX, Real diffY,
                     Real theta, Real omega,
                     int numMolecules);
};

}// namespace gpgmp
#endif // NONLINEARPROBLEM_H
