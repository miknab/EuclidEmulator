cmake_minimum_required (VERSION 2.8.5)
project ( EUCLID_EMULATOR )

# Set message colors
if (NOT WIN32)
    string( ASCII 27 Esc )
    set ( ColourReset "${Esc}[m" )
    set ( Red "${Esc}[31m" )
    set ( Green "${Esc}[32m" )
endif()
####


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

# add the compilation of the class patch as an external project
set ( USEPATCH TRUE CACHE BOOL "Compile and use the patched version of classy (called ''classee'')" )

if( ${USEPATCH} )
    message ( "Building patched class as external project" )
    
    file( MAKE_DIRECTORY ClassPatch )

    include( ExternalProject )

    ExternalProject_Add(project_eeClassPatch
        URL https://github.com/miknab/ClassPatch/archive/master.tar.gz
        SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/ClassPatch
        DOWNLOAD_DIR ${CMAKE_CURRENT_LIST_DIR}/ClassPatch
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND ""
        BUILD_COMMAND make
        INSTALL_COMMAND make
    )

    ExternalProject_Get_Property(project_eeClassPatch install_dir)

    add_library ( eeClassPatch MODULE IMPORTED )
else()
    message ( WARNING "${Red}You will be using EuclidEmulator with a non-patched version of class/classy. This might lead to unexpected or even wrong results. Do this at your own risk! We encourage you to use the patched version. We do not take any responsibility for anything that could go wrong because of using EuclidEmulator with the original classy version.${ColourReset}" )
endif()
>>>>>>> AddClassPatch
