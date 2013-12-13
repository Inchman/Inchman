#ifndef __gpgmp_Globals_h__
#define __gpgmp_Globals_h__

#define NOMINMAX /* force Visual C++ to not generate min and max macros, so that we can use the cross-platform std::min and std::max */

#define __CL_ENABLE_EXCEPTIONS
#define CL_USE_DEPRECATED_OPENCL_1_1_APIS /** Need these to use cl.hpp (which is not yet 1.2).. TODO: Move to our own CL classes! */
//#include "CL/cl.hpp"
#include "cl.hpp"

#define HDF_REAL H5T_NATIVE_FLOAT

#define NA 6.0221415e23

#define GMP_DIMENSIONALITY 3

// number of error values reserved
#define GMP_NUM_ERRORS 20 // TODO: remove this

// number of threads per block dim
// TODO: We really need to make this dynamic somehow..
#ifdef __APPLE__
#define LOCAL_RANGE_BLOCK 8
#else
#define LOCAL_RANGE_BLOCK 16
#endif

#include <cstring>

// define some kernel parameters for the inhomogeneous diffusion algorithm
namespace gpgmp {

// type definitions for stoichiometry values - TODO: we don't really need these anymore ..
typedef int StoichiometryEntry;

// type definitions for a species index - TODO: we don't really need these anymore .. 
typedef int SpeciesIndex;

// use this type for Real values from now on
 typedef cl_float Real;

// These are standard methods to compute the diffusivity and drift field
 enum DiffusionType {
     DT_HOMOGENEOUS,
     DT_MULTIPLICATIVE_NOISE,
     DT_ORNSTEIN_UHLENBECK,
     DT_NONLINEAR
 };

 /*
   DT_LOCAL_AVERAGE,
   DT_STRICTLY_LOCAL,
   DT_NEIGHBOUR_BASED,
   DT_GRADIENT_BASED,
   DT_ORNSTEIN_UHLENBECK,
   DT_NONLINEAR
*/

 // Helper function to get error string
 // this is from the NVIDIA SDK so we probably shouldn't distribute it..
 // *********************************************************************

} // namespace gpgmp

//todo: integrate somewhere..
const char* oclErrorString(cl_int error);
void printRange(std::ostream &stream, const cl::NDRange &range);

#endif // !__gpgmp_Globals_h__
