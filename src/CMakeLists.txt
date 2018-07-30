set( CYTHON_SOURCE_FILE EuclidEmulator_BackEnd.pyx )
set( C_SOURCE_FILES EuclidEmulator.c cosmo.c )

# Set message colors
if (NOT WIN32)
    string( ASCII 27 Esc )
    set ( ColourReset "${Esc}[m" )
    set ( Red "${Esc}[31m" )
    set ( Green "${Esc}[32m" )
endif()
####

execute_process ( COMMAND bash -c "python -c 'import sys; print(sys.version)'"
                  OUTPUT_VARIABLE PYTHON_VERSION
                )

message ( STATUS "${Green}The module is compiled for python version ${PYTHON_VERSION}${ColourReset}" )

execute_process ( COMMAND bash -c "python -c 'import numpy; print(numpy.get_include())'"
                  OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
                )

# Notice that the default Numpy C API used by cython is depricated and
# cimport-ing numpy into the python code will hence raise a warning (see
# the cython docs here https://cython.readthedocs.io/en/latest/src/reference/compilation.html)
# For the time this warning can just be ignored. But people are working on this
# issue, see e.g. https://github.com/cython/cython/issues/2498.
#
# In the case of EuclidEmulator we don't have to cimport numpy but only import
# it. This does NOT raise this warning.
include_directories( ${EUCLID_EMULATOR_SOURCE_DIR}/include 
#                     ${NUMPY_INCLUDE_DIR}
                   )

# If necessary at some point, using the full C-API for numpy is easy: just
# uncomment the line including {NUMPY_INCLUDE_DIR} and in EuclidEmulator_BackEnd.pyx
# "cimport numpy" (instead of just "import numpy")

cython_add_module( EuclidEmulator_BackEnd ${CYTHON_SOURCE_FILE} ${C_SOURCE_FILES} )

add_executable( ee master.c ${C_SOURCE_FILES} )

#include GSL library
find_package(GSL REQUIRED)

if ( !GSL_FOUND )
    message ( SEND_ERROR "${Red}GSL library not installed${ColourReset}" )
elseif ( GSL_VERSION VERSION_LESS 1.13 )
    message ( SEND_ERROR "${Red}Required minimum version of gsl library is 1.13${ColourReset}") 
endif()

target_link_libraries( ee GSL::gsl GSL::gslcblas)
target_link_libraries( EuclidEmulator_BackEnd GSL::gsl GSL::gslcblas)










