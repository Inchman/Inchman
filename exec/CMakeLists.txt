project(inchman-exec)

# we need the FindHDF5 script, as found in CMake 2.8
cmake_minimum_required(VERSION 2.8.0)

# debug definitions
set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG") # for assert()
set(CMAKE_C_FLAGS_RELEASE "-DNDEBUG")   # for assert()
add_definitions(-Wall)
add_definitions(-Wextra)

# define gpgmp_cl source dir
set (GPGMP_CL_ROOT_DIR ${inchman-exec_SOURCE_DIR}/src/gpgmp/)

# define python source dir 
set (PYTHON_SOURCE_DIR ${inchman-exec_SOURCE_DIR}/../python/)

# OpenCL
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
include(FindOpenCL REQUIRED)

# Python - need version 2.7 so it's consistent with the boost python
FIND_PACKAGE(PythonInterp REQUIRED)
FIND_PACKAGE(PythonLibs REQUIRED)

# HDF5
include(FindHDF5 REQUIRED)

# SBML
find_library(SBML_LIBRARIES sbml HINTS /usr/lib /usr/lib64/ /usr/local/lib $ENV{HOME}/lib)
if (${SBML_LIBRARIES} STREQUAL "SBML_LIBRARIES-NOTFOUND")
  message(FATAL_ERROR "Error: Could not locate SBML libraries.")
else (${SBML_LIBRARIES} STREQUAL "SBML_LIBRARIES-NOTFOUND")
  message(STATUS "Found SBML: ${SBML_LIBRARIES}")
endif(${SBML_LIBRARIES} STREQUAL "SBML_LIBRARIES-NOTFOUND")
find_path(SBML_INCLUDE_DIR sbml HINTS /usr/include /usr/local/include $ENV{HOME}/include)
if (${SBML_INCLUDE_DIR} STREQUAL "SBML_INCLUDE_DIR-NOTFOUND")
  message(FATAL_ERROR "Error: Could not locate SBML include directory.")
endif(${SBML_INCLUDE_DIR} STREQUAL "SBML_INCLUDE_DIR-NOTFOUND")

# boost libraries
set(Boost_USE_STATIC_LIBS        OFF)
set(Boost_USE_MULTITHREADED      OFF)
set(Boost_USE_STATIC_RUNTIME    OFF)
find_package( Boost 1.54.0 COMPONENTS log chrono date_time filesystem log_setup system thread program_options python REQUIRED)
add_definitions(-DBOOST_ALL_DYN_LINK)

# recursive make
add_subdirectory(src)
