cmake_minimum_required (VERSION 2.8.5)
project ( EuclidEmulator_CLI )

# Process the CMakeLists.txt in the 'src' and 'bin' directory
add_subdirectory( ../libeuc ../../libeuc/build )
add_executable( ee master.c )

target_link_libraries( ee eucemu )

install( TARGETS ee DESTINATION "bin" )
