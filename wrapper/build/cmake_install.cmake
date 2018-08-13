# Install script for directory: /home/ics/mischak/code/EuclidEmulator/EmulatorCode/EucEmu_V1.1.1.clean/NewStructure/wrapper

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/home/ics/mischak/anaconda2/lib/python2.7/site-packages/EuclidEmulator_BackEnd.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/ics/mischak/anaconda2/lib/python2.7/site-packages/EuclidEmulator_BackEnd.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/home/ics/mischak/anaconda2/lib/python2.7/site-packages/EuclidEmulator_BackEnd.so"
         RPATH "")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/ics/mischak/anaconda2/lib/python2.7/site-packages/EuclidEmulator_BackEnd.so")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/ics/mischak/anaconda2/lib/python2.7/site-packages" TYPE MODULE FILES "/home/ics/mischak/code/EuclidEmulator/EmulatorCode/EucEmu_V1.1.1.clean/NewStructure/wrapper/build/EuclidEmulator_BackEnd.so")
  if(EXISTS "$ENV{DESTDIR}/home/ics/mischak/anaconda2/lib/python2.7/site-packages/EuclidEmulator_BackEnd.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/ics/mischak/anaconda2/lib/python2.7/site-packages/EuclidEmulator_BackEnd.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/home/ics/mischak/anaconda2/lib/python2.7/site-packages/EuclidEmulator_BackEnd.so"
         OLD_RPATH "/home/ics/mischak/anaconda2/lib:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/home/ics/mischak/anaconda2/lib/python2.7/site-packages/EuclidEmulator_BackEnd.so")
    endif()
  endif()
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/ics/mischak/anaconda2/lib/python2.7/site-packages/e2py")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/ics/mischak/anaconda2/lib/python2.7/site-packages" TYPE DIRECTORY FILES "/home/ics/mischak/code/EuclidEmulator/EmulatorCode/EucEmu_V1.1.1.clean/NewStructure/wrapper/e2py")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/ics/mischak/code/EuclidEmulator/EmulatorCode/EucEmu_V1.1.1.clean/NewStructure/libeuc/build/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/ics/mischak/code/EuclidEmulator/EmulatorCode/EucEmu_V1.1.1.clean/NewStructure/wrapper/build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
