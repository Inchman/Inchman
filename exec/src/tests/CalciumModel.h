#ifndef CALCIUMMODEL_H
#define CALCIUMMODEL_H

#include "DiffusionModel.h"

namespace gpgmp{

/**
 Implements a simple model for intracellular Calcium dynamics in Xenopus oocyte cells. Details of the model
 can be found in Vigelius&Meyer (2012b).

 \see \vigelius12
  */
class CalciumModel : public DiffusionModel
{
public:
    /**
      Sets up the model.

      \param length Size of the integration domain in \f$\mu\mathrm{m}\f$
      \param dx, dy Number of grid cells per dimension
      \param diffc, diffp Diffusivity of Calcium ions and \f$\mathrm{IP}_3\f$, respectively
      \param nmax Maximum number of open channels
      \param kp \f$\mathrm{IP}_3\f$ break-down rate
      \param betac Leakage of Calcium ions
      \param omega Scaling parameter
      */
    CalciumModel(Real length, int dx, int dy,
                 Real diffc, Real diffp,
                 Real nmax, Real kp, Real betac,
                 Real omega,
                 const std::string &initScriptsPath);
};

} // namespace gpgmp
#endif // CALCIUMMODEL_H
