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

# The way the PYTHON_MODULE_PATH might seem a bit strange. There are two main reasons
# why it is done this way:

# 1) We go the detour via PYTHON_EXECUTABLE instead of accessing PYTHON_LIBRARIES
#    directly as this can lead to problems when the found python interpreter was
#    installed with Homebrew. Going the detour ensures that the found interpreter
#    and libraries belong to the same python installation

# 2) There is a function "FindPython" in CMake 3.12 (most recent version of CMake in Aug 2018).
#    that has a return value "Python_SITELIB" which is basically just what we need.
#    But we decided not to rely on it as because it is so recent.

#    ---------------------------------- I M P O R T A N T ---------------------------------------     #
#       ...for users without root access to the machine they try to build EuclidEmulator on.          #
#                                                                                                     #
#    You have two options to build EuclidEmulator on a machine without root access:                   #
#                                                                                                     #
#    Option A:  Build EuclidEmulator within a virtual environment.                                    #
#                                                                                                     #
#    Option B:  Modify this cmake file. To do so, please read each of the following comment line      #
#               carefully and follow the step-by-step instructions.                                   #
#                                                                                                     #
#    For further info, please consult the EuclidEmulator wiki page (https://github.com/miknab/        #
#    EuclidEmulator/wiki/II)-From-download-to-installation).                                          #
#                                                                                                     # 
#    ---------------------------------------------------------------------------------------------    #

# Modification for "Option B": comment the following 'execute_process' as well as the 'message' command
execute_process ( COMMAND
                  ${PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_python_lib; print(get_python_lib())"
                  OUTPUT_VARIABLE PYTHON_SITE_PACKAGES
                  OUTPUT_STRIP_TRAILING_WHITESPACE
                )
message( "Path to python site-packages: ${PYTHON_SITE_PACKAGES}" )

# Modifications for "Option B": uncomment the next two lines and edit the path in quotes setting the
# variable 'PYTHON_MODULE_PATH'. 
#get_filename_component( PYEXPATH ${PYTHON_EXECUTABLE} DIRECTORY )
#set( PYTHON_MODULE_PATH "path/to/userdir")

# Compile
cython_add_module( EuclidEmulator_BackEnd ${CYTHON_SOURCE_FILE} )

# Link
target_link_libraries ( EuclidEmulator_BackEnd eucemu )

# Install
# Modifications for "Option B": comment the next two lines
install( TARGETS EuclidEmulator_BackEnd DESTINATION ${PYTHON_SITE_PACKAGES} )
install( DIRECTORY e2py DESTINATION ${PYTHON_SITE_PACKAGES} )

# Modifications for "Option B": uncomment the next two lines
#install( TARGETS EuclidEmulator_BackEnd DESTINATION ${PYTHON_MODULE_PATH} )
#install( DIRECTORY e2py DESTINATION ${PYTHON_MODULE_PATH} )
