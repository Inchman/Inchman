include_directories(
   ${inchman-exec_SOURCE_DIR}/src
   ${inchman-exec_SOURCE_DIR}/src/gpgmp
   ${inchman-exec_SOURCE_DIR}/src/inchman
   ${GPGMP_CL_ROOT_DIR}/src
   ${Boost_INCLUDE_DIRS}
   ${HDF5_INCLUDE_DIR}
   ${SBML_INCLUDE_DIR}
   ${PYTHON_INCLUDE_DIRS}
   ${OPENCL_INCLUDE_DIRS}
)

link_directories(
   ${inchman-exec_BINARY_DIR}/src/gpgmp
   ${CUDA_TOOLKIT_ROOT_DIR}/lib
   ${HDF5_LIBRARY_DIRS}
)

add_definitions(${HDF5_DEFINITIONS})

set(inchman-exec_SRCS
   main.cpp
#   InchmanDiffusionModel.h
   InchmanDiffusionModel.cpp
)

source_group(
   NeatSources
   FILES
   ${inchman-exec_SRCS}
)

if (UNIX AND NOT APPLE)
   set(EXTRA_LIBS util dl)
endif (UNIX AND NOT APPLE)

add_executable(inchman-exec ${inchman-exec_SRCS})

target_link_libraries(
   inchman-exec
   gpgmp
   ${Boost_LIBRARIES}
   ${SBML_LIBRARIES}
   ${HDF5_LIBRARIES}
   hdf5_hl # didn't work (at least on the Mac it didn't): ${HDF5_HL_LIBRARIES}
   ${PYTHON_LIBRARIES}
   ${OPENCL_LIBRARIES}
 #  ${EXTRA_LIBS}
 util
 dl
)

# we need to set the runtime path to be the link path..
# otherwise "make install" will strip custom locations of boost if they
# are in non-standard locations (e.g. /usr/local/ instead of /usr/)
set_target_properties(inchman-exec PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)

# Install targets

# executable
install (TARGETS inchman-exec DESTINATION bin)

# cl templates
install (FILES
	"${GPGMP_CL_ROOT_DIR}/src/Deterministic.cl"
	"${GPGMP_CL_ROOT_DIR}/src/Gillespie.cl"
	"${GPGMP_CL_ROOT_DIR}/src/IndividualSolver.cl"
	"${GPGMP_CL_ROOT_DIR}/src/StochasticInhomogeneousDriftDiffusion.cl"
	"${GPGMP_CL_ROOT_DIR}/src/StochasticHomogeneousDiffusion.cl"
         DESTINATION share/inchman/cl_source/)
install (DIRECTORY "${GPGMP_CL_ROOT_DIR}/src/cl_include" DESTINATION share/inchman/cl_source/)

# configure the config file and copy it into the python dest-dir
install(CODE "MESSAGE(STATUS ${INCHMAN_EXEC})")
install(CODE "CONFIGURE_FILE(${PYTHON_SOURCE_DIR}/gpgmp/common/options.py.in ${PYTHON_SOURCE_DIR}/gpgmp/common/options.py)")

# fire up the python install script - we let the user do that instead so they can just
# set the pythonpath if they don't have root access
# install(CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE} ${PYTHON_SOURCE_DIR}/setup.py install)")

# install inchman script
install (PROGRAMS "${PYTHON_SOURCE_DIR}/inchman" DESTINATION bin)