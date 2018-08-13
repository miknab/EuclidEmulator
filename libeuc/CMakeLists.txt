cmake_minimum_required(VERSION 2.4)

project ( eucemu )

# The version number
set ( EUCLID_EMULATOR_VERSION_MAJOR 1 )
set ( EUCLID_EMULATOR_VERSION_MINOR 1 )


# Set the path to where the EuclidEmulator data file ee.dat should be installed
set ( PATH_TO_EUCLIDEMULATOR_DATA_FILE "${CMAKE_INSTALL_PREFIX}/share" )

# Configure the header file to pass CMake setting to source code
configure_file ( "EuclidEmulator.h.in"
                 "${CMAKE_CURRENT_BINARY_DIR}/EuclidEmulator.h"
               )

add_library( ${PROJECT_NAME} STATIC EuclidEmulator.c cosmo.c EuclidEmulator.h.in cosmo.h )
set_property( TARGET ${PROJECT_NAME} PROPERTY POSITION_INDEPENDENT_CODE ON )
target_include_directories( ${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR} )

#include GSL library
find_package(GSL REQUIRED )

#include math library
find_library(LIBM REQUIRED)

target_link_libraries( ${PROJECT_NAME} ${GSL_LIBRARIES} )
target_include_directories( ${PROJECT_NAME} PUBLIC ${GSL_INCLUDE_DIRS} )


install( TARGETS ${PROJECT_NAME} DESTINATION "lib" )
install( FILES ${CMAKE_CURRENT_BINARY_DIR}/EuclidEmulator.h DESTINATION "include" )
install( FILES ee.dat DESTINATION ${PATH_TO_EUCLIDEMULATOR_DATA_FILE} )


