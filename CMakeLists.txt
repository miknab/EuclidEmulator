cmake_minimum_required (VERSION 2.8.5)
project ( EUCLID_EMULATOR )

# Make the scripts available in the 'cmake_scripts' directory available for the
# 'include()' command, 'find_package()' command.
set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_LIST_DIR}/cmake_scripts )

# Include the CMake script UseCython.cmake.  This defines add_cython_module().
# Instruction for use can be found at the top of cmake/UseCython.cmake.
include( UseCython )

# Procedd the CMakeLists.txt in the 'src' and 'bin' directory
add_subdirectory( src )

