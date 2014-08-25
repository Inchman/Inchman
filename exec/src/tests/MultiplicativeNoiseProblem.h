#ifndef MULTIPLICATIVENOISEPROBLEM_H
#define MULTIPLICATIVENOISEPROBLEM_H

#include "DiffusionModel.h"

namespace gpgmp {
/**
  Geometric Brownian motion (or multiplicative noise) problem. This problem is
  defined by the stochastic differential equation

  \f[
     \mathrm{d}\mathbf{X}_t = \mathbf{\mu} \mathbf{X}_t \mathrm{d}t + \mathbf{\sigma}\, \mathbf{X}_t\, \mathrm{d}\mathbf{W}_t,
  \f]

  with \f$\mathbf{\mu}=\mathrm{diag} \{\mu_x, \mu_y\}\f$ and \f$\mathbf{\sigma}=\mathrm{diag}\{\sigma_x, \sigma_y\}\f$ diagonal matrices.

  Details about the problem can be found in Vigelius&Meyer (2012). The initial set up consists of one species \f$A\f$ which is initially
  located in the center of the domain. No reactions are defined.

  The problem can be run with the test suite using
  \verbatim
    ./run-tests.py --multiplicative-noise
  \endverbatim

  \see \vigelius12
  */
class MultiplicativeNoiseProblem : public DiffusionModel
{
public:
    /**
      Sets up the problem.

      \param length Size of the integration domain in \f$\mu\mathrm{m}\f$.
      \param dx, dy Number of compartments in \f$x\f$ and \f$y\f$ direction
      \param numMolecules Number of molecules of species \f$A\f$
      \param cx, cy Matrix entries for diffusivity \f$\sigma\f$ (in \f$\mathrm{s}^{-1/2}\f$)
      \param mux, muy Matrix entries for drift \f$\mu\f$ (in \f$\mathrm{s}^{-1}\f$)
      */
    MultiplicativeNoiseProblem(Real length, int dx, int dy, int numMolecules,
                               Real cx, Real cy, Real mux, Real muy);
};
} // namespace gpgmp
#endif // MULTIPLICATIVENOISEPROBLEM_H
