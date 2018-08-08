cmake_minimum_required (VERSION 2.8.5)
project ( EUCLID_EMULATOR )

# The version number
set ( EUCLID_EMULATOR_VERSION_MAJOR 1 )
set ( EUCLID_EMULATOR_VERSION_MINOR 1 )

set ( TOP_PROJECT_DIR ${CMAKE_CURRENT_LIST_DIR} )

# Configure the header file to pass CMake setting to source code
configure_file ( "${CMAKE_CURRENT_LIST_DIR}/include/EuclidEmulator.h.in"
                 "${CMAKE_CURRENT_LIST_DIR}/include/EuclidEmulator.h"
               )

# add_definitions( "-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION" )

# Make the scripts available in the 'cmake_scripts' directory available for the
# 'include()' command, 'find_package()' command.
set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_LIST_DIR}/cmake_scripts )

# Include the CMake script UseCython.cmake.  This defines add_cython_module().
# Instruction for use can be found at the top of cmake/UseCython.cmake.
include( UseCython )

# Procedd the CMakeLists.txt in the 'src' and 'bin' directory
add_subdirectory( src )
