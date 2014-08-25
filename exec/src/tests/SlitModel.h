#ifndef SLITMODEL_H
#define SLITMODEL_H

#include "DiffusionModel.h"

namespace gpgmp {
/**
  This class sets up a model for cell migrations from neurons in the brain under the
  influence of a signalling molecule, SLIT. Details can be found in Vigelius&Meyer (2012).

  The standard model consists of one species representing the migrating cells. The
  concentration of SLIT after administration is assumed to be constant in time. In
  keeping with the experiment, we prescribe a one-dimensional exponential distribution
  of SLIT with a variable scale length. The interaction of neurons with SLIT is modelled
  as an exponential function while cell-cell interaction is implemented using the heuristic
  contact inhibition model.

  \see Vigelius12
*/
class SlitModel : public DiffusionModel
{
public:
    /**
      Constructs the model.

      \param length Size of the integration domain in \f$\mu\mathrm{m}\f$
      \param dx, dy Number of grid cells in each dimension
      **/
    SlitModel(Real length, int dx, int dy,
              const std::string &initScriptsPath);
};
} // namespace GPGMP

#endif // SLITMODEL_H
