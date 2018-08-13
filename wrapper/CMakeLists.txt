cmake_minimum_required ( VERSION 3.1 )
project ( EUCEMU_PYMOD )

# Set message colors
if (NOT WIN32)
    string( ASCII 27 Esc )
    set ( ColourReset "${Esc}[m" )
    set ( Red "${Esc}[31m" )
    set ( Green "${Esc}[32m" )
endif()
####

# Set some variables
set( CYTHON_SOURCE_FILE EuclidEmulator_BackEnd.pyx )

# Make the scripts available in the 'cmake_scripts' directory available for the
# 'include()' command, 'find_package()' command.
set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_LIST_DIR}/cmake_scripts )

# Add the euclid emulator library
add_subdirectory( ../libeuc ../../libeuc/build )

# Include the CMake script UseCython.cmake.  This defines add_cython_module().
# Instruction for use can be found at the top of cmake/UseCython.cmake.
include( UseCython )

get_filename_component( PYLIBPATH ${PYTHON_LIBRARIES} DIRECTORY )
set( PYTHON_MODULE_PATH "${PYLIBPATH}/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages" )

# Compile
cython_add_module( EuclidEmulator_BackEnd ${CYTHON_SOURCE_FILE} )

# Link
target_link_libraries ( EuclidEmulator_BackEnd eucemu )

# Install
install( TARGETS EuclidEmulator_BackEnd DESTINATION ${PYTHON_MODULE_PATH} )
install( DIRECTORY e2py DESTINATION ${PYTHON_MODULE_PATH} )
