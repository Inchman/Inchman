#include "SlitModel.h"

#include <boost/filesystem.hpp>

#include <sstream>

namespace gpgmp {
/**
  Currently, all parameters are hard-coded into the constructor. We dynamically construct
  the source code for the drift-diffusion method here with the parameters passed in the source
  rather than using the GPGMP_PARAMS_.. logic. Since the local diffusivity depends on the
  cell concentration (contact inhibition model), we need to set it up as a non-linear model.
  We also require the init script init_slit.py to set the initial configuration

  \todo allow different diffusion models etc .. will come with the real models!

  \see DiffusionModel::setComputeDriftDiffusivityMethod
  \see DiffusionModel::setNonlinearDiffusivity
  \see DiffusionModel::loadInitScript
  */
SlitModel::SlitModel(Real length, int dx, int dy,
                     const std::string &initScriptsPath)
    :   DiffusionModel("AB Annihilation", true)
{
  setPhysicalLength(length);
  setGridDims(dx, dy);
    Real diffCell = 0.56;
    addSpecies("Cells", diffCell);

    // set boundary mask
    const cl_int4 boundaryMask = {{1, 1, 1, 1}};
    const cl_int4 sourceMask   = {{0, 0, 0, 0}};
    setBoundaryMasks(boundaryMask, sourceMask);

    // Subvolume size (in mum^2)
    Real subVolume = length/dx*length/dx;

     // set diffusion type
    int diffType = 1;
    setParameter("diffType", diffType);

    // slit parameters
    Real lambda = 0.01;
    //Real beta = 20;
    Real beta = 400;

    // set diffusion method
    // todo: change according to SLIT type..
    std::ostringstream ss;
    ss << "  // Parameters\n"
          "  const Real L = "<<length/2.<<";\n"
          "  const Real lambda = "<<lambda<<";\n"
          "  const Real a = "<<0.02*subVolume<<";\n"
          "  const Real ta = "<<24.*3600.<<";\n"
          "  const Real xi = "<<1e4/3600.<<";\n"
          "  const Real beta = "<<beta<<";\n"
          "  const Real diffCell = "<<diffCell<<";\n"
          "  \n"
          "  // Is slit switched on already?\n"
      //          "  const Real s0 = (SimTime > ta) ? 1 : 0;\n"
          "  const Real s0 = (SimTime > ta) ? 1 : 1;\n"
          "  \n"
          "  // compute slit concentration\n"
          "  // note that (cx-L) should always be negative..\n"
          "  const Real cx = "<<length/Real(dx)<<"*(get_global_id(0)+0.5)-L;\n"
          "  const Real slitConcentration = exp(-lambda*(L-cx));\n"
          "  const Real slitConcentrationPrime = lambda*exp(-lambda*(L-cx));\n"
          "\n"
          "  const Real apu = a + _state[_index];\n"
          "  const Real a0 = a * diffCell * exp(-beta*s0*slitConcentration);\n"
          "  const Real up = getCentralDifferenceX(_state, _speciesIndex)/(2.*Length/(GridModelWidth));\n"
          "  const Real sp = s0*slitConcentrationPrime;\n"
          "\n"
          "  // and set X diffusion/drift\n"
          "  Cells->DriftX = -a0/(apu*apu) * (beta * apu * sp + up);\n"
          "  Cells->DiffusivityX = a0/apu;"
          "\n"
          "  // for y direction we take only population pressure (no slit gradient here)\n"
          "  const Real uprimey = getCentralDifferenceY(_state, _speciesIndex)/(2.*Length/(GridModelHeight));"
          "  Cells->DriftY = -diffCell* a * uprimey/(apu*apu);\n"
          "  Cells->DiffusivityY = diffCell*a/apu;\n"
          ;
    setComputeDriftDiffusivityMethod(ss.str());
    setNonlinearDiffusivity(true);

    // load init script
    loadInitScript((boost::filesystem::path(initScriptsPath) / "init_slit.py").string());
} // constructor
} //namespace gpgmp
