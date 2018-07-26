set( CYTHON_SOURCE_FILE EuclidEmulator_BackEnd.pyx )
set( C_SOURCE_FILES EuclidEmulator.c cosmo.c )

execute_process ( COMMAND bash -c "python -c 'import sys; print(sys.version)'"
                  OUTPUT_VARIABLE PYTHON_VERSION
                )

message ( STATUS "The module is compiled for python version ${PYTHON_VERSION}" )

execute_process ( COMMAND bash -c "python -c 'import numpy; print(numpy.get_include())'"
                  OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
                )

include_directories( ${EUCLID_EMULATOR_SOURCE_DIR}/include 
                     ${NUMPY_INCLUDE_DIR}
                   )

cython_add_module( EuclidEmulator_BackEnd ${CYTHON_SOURCE_FILE} ${C_SOURCE_FILES} )

add_executable( ee master.c ${C_SOURCE_FILES} )

#include GSL library
find_package(GSL REQUIRED)
target_link_libraries( ee GSL::gsl GSL::gslcblas)
target_link_libraries( EuclidEmulator_BackEnd GSL::gsl GSL::gslcblas)
