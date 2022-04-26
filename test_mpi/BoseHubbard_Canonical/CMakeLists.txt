# we need at least cmake version 3.1
cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

project(qmc_worm CXX)

add_compile_options(
  "-Wall" "-Wextra"
  "-Wno-unused-parameter"
  "-Wno-reorder"
  "-Wfatal-errors"
  "-Wno-unknown-pragmas"
  "-fdiagnostics-color=always"
  )

SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(DEBUG OFF CACHE BOOL "Build in Debug Mode")
if(DEBUG)
  add_definitions(-DDEBUGMODE)
  set(CMAKE_BUILD_TYPE "Debug")
  message(STATUS "Build type is Debug")
else()
  remove_definitions(-DDEBUGMODE)
  set(CMAKE_BUILD_TYPE "Release")
  message(STATUS "Build type is Release")
endif()


add_definitions(-DUNISYS)
add_definitions(-DDSTRUC_LIST)
add_definitions(-DCAN_WINDOW)


set(LATTICE "chain" CACHE STRING "Selection of the Lattice")
if(LATTICE STREQUAL chain)
  message(STATUS "Using chain lattice")
elseif(LATTICE STREQUAL square)
  message(STATUS "Using square lattice")
elseif(LATTICE STREQUAL cubic)
  message(STATUS "Using cubic lattice")
elseif(LATTICE STREQUAL ladder)
  message(STATUS "Using ladder lattice")
elseif(LATTICE STREQUAL triangular)
  message(STATUS "Using triangular lattice")
elseif(LATTICE STREQUAL honeycomb)
  message(STATUS "Using honeycomb lattice")
else()
  message(FATAL_ERROR "Unknown lattice ${LATTICE}")
endif()
add_compile_definitions(LATTICE=${LATTICE})


add_executable(${PROJECT_NAME}_mpi
  ../../src/worm.run_mpi.cpp
  ../../src/worm.cpp
  ../../src/worm.update.cpp
  ../../src/worm.output.cpp)


# the project relies on the ALPSCore package. If ALPSCore is not found
# automatically, specify its location using:
# export ALPSCore_DIR=/location/to/ALPSCORE/
find_package(ALPSCore 2.0  REQUIRED mc accumulators params)

# check if Eigen is present:
if (ALPSCore_HAS_EIGEN_VERSION)
  message(STATUS "ALPSCore provides Eigen, version ${ALPSCore_HAS_EIGEN_VERSION}")
else()
  message("WARNING: ALPSCore does not provide Eigen, the code won't compile!")
endif()
# Use ALPSCore_LIBRARIES variable to link to ALPSCore 

target_link_libraries(${PROJECT_NAME}_mpi ${ALPSCore_LIBRARIES})
