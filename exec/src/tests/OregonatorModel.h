#ifndef OREGONATORMODEL_H
#define OREGONATORMODEL_H

#include "DiffusionModel.h"

namespace gpgmp {
/**
  Sets up the classical Oregonator model. This model is defined by five species: \f$A\f$ and \f$B\f$, which are assumed
  to be constant over time and space, an activator species \f$X\f$, an inhibitor species \f$Z\f$ and an intermediary \f$Y\f$.
  In order to incorporate rational product stoichiometries, we also introduce two pseudo-species \f$Y^\ast\f$ and \f$Y^{\ast\ast}\f$
  [see Vigelius&Meyer (2012b) for details]. The reaction scheme is given by

\f{eqnarray*}{
  A + Y & \stackrel{k_1}{\rightarrow} & X + R \\
  X + Y & \stackrel{k_2(\Omega)}{\rightarrow} & 2 R \\
  A + X & \stackrel{k_3}{\rightarrow} & 2 X + 2 Z \\
  2 X & \stackrel{k_4(\Omega)} {\rightarrow} & A + R \\
  B + Z & \stackrel{k_5}{\rightarrow} & Y^\ast + Y^{\ast\ast} \\
  2 Y^\ast & \stackrel{k_6}{\rightarrow} & Y \\
  2 Y^{\ast\ast} & \stackrel{k_7}{\rightarrow} & Y^\ast.
\f}

 Similar to GrayScottModel, \f$\Omega\f$ is a scaling parameter.

 \see GrayScottModel
 \see \vigelius12b
  **/
class OregonatorModel : public DiffusionModel
{
public:
    /**
      Sets up the model.

      \param length Size of the integration domain in \f$\mu\mathrm{m}\f$
      \param dx, dy Number of grid cells in each dimension
      \param Omega Scaling parameter \f$\Omega\f$
      */
    OregonatorModel(Real length, int dx, int dy,
                    Real omega, const std::string &initScriptsPath);
};

} //namespace gpgmp
#endif // OREGONATORMODEL_H
