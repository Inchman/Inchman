# include src dir for headers
INCLUDE_DIRECTORIES(
  ${GPGMP_CL_ROOT_DIR}/src
  ${PYTHON_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS}
  ${OPENCL_INCLUDE_DIRS}
)

# add the correct HDF5 version to include dirs
INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})

# c++ sources
set (GPGMP_SRC
  ${GPGMP_CL_ROOT_DIR}/src/Base.cpp
  ${GPGMP_CL_ROOT_DIR}/src/Compartment.cpp
  ${GPGMP_CL_ROOT_DIR}/src/DiffusionModel.cpp
  ${GPGMP_CL_ROOT_DIR}/src/Kernel.cpp
  ${GPGMP_CL_ROOT_DIR}/src/Reaction.cpp
  ${GPGMP_CL_ROOT_DIR}/src/Species.cpp
  ${GPGMP_CL_ROOT_DIR}/src/IndividualSolver.cpp
  ${GPGMP_CL_ROOT_DIR}/src/Solver.cpp
  ${GPGMP_CL_ROOT_DIR}/src/InhomogeneousSolver.cpp
  ${GPGMP_CL_ROOT_DIR}/src/HomogeneousSolver.cpp
)

# main target
add_library(gpgmp ${GPGMP_SRC})