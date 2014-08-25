#ifndef GRAYSCOTTMODEL_H
#define GRAYSCOTTMODEL_H

#include "DiffusionModel.h"

namespace gpgmp {

 /**
  This class implements the Gray-Scott model. The model consists of two species, \f$U\f$ and \f$V\f$, and a reaction network
\f{eqnarray*}{
  U + 2V & \stackrel{k_1(\Omega)}{\rightarrow} & 3 V \\
  U & \stackrel{k_f}{\rightarrow} & \emptyset \\
  V & \stackrel{k_2}{\rightarrow} & \emptyset \\
  \emptyset & \stackrel{k_f u_0(\Omega)}{\rightarrow} & U,
\f}
  where \f$\Omega\f$ allows us to scale the number of particles per subvolume while keeping the
  reaction dynamics unchanged [see Vigelius&Meyer (2012b) for details].
  We furthermore introduce the control parameters \f$F=k_f\f$ and \f$k=k_2-F\f$ such that the properties of our system
  are fully determined by its position in the \f$F-k\f$ plane.
  \see \vigelius12b
  */

 class GrayScottModel : public DiffusionModel
 {
 public:
     /**
       Sets up the model.

       \param length Size of the integration domain in \f$\mu\mathrm{m}\f$
       \param dx, dy Number of grid cells in each dimension
       \param diffU, diffV Diffusivity for both species
       \param k Parameter \f$k\f$ in \f$\mathrm{s}^{-1}\f$
       \param F Parameter \f$F\f$ in \f$\mathrm{s}^{-1}\f$
       \param omega Scaling parameter \f$\Omega\f$
       */
     GrayScottModel(Real length, int dx, int dy,
                    Real diffU, Real diffV,
                    Real k, Real F, Real omega,
                    const std::string &initScriptsPath);
 };
}// namespace gpgmp

#endif // GRAYSCOTTMODEL_H
