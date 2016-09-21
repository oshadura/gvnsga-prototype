cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)
project(libga)
#-------------------
# Requirements
#-------------------
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")
#-------------------
# Dependency for MVA
#find_package(Eigen3 REQUIRED)

#-------------------
# Compiler & Language
#-------------------

enable_language(CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
     message(STATUS "Detecting: Clang")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
     message(STATUS "Detecting:  GCC")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
      message(STATUS "Detecting: Intel C++")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
      message(STATUS "Detecting: MSVC")
endif()

#-------------------
# Flags
#-------------------

if (GNU)
  set(GA_CXX_FLAGS "${GA_CXX_FLAGS} -Wall -fPIC -g -z defs")
  if (CMAKE_BUILD_TYPE MATCHES Debug)
    set(GA_CXX_FLAGS "${GA_CXX_FLAGS} -ggdb -O0")
  else()
    set(GA_CXX_FLAGS "${GA_CXX_FLAGS} -finline-limit=10000000 -ffast-math -ftree-vectorize")
  endif()
  if (Vc)
    set(GA_CXX_FLAGS "${GA_CXX_FLAGS} -fabi-version=6")
  endif()
elseif(Intel)
  set(GA_CXX_FLAGS "${GA_CXX_FLAGS} -Wall -fPIC -diag-disable 3438")
  if (CMAKE_BUILD_TYPE MATCHES Debug)
    set(GA_CXX_FLAGS "${GA_CXX_FLAGS} -g -O0")
  else()
    set(GA_CXX_FLAGS "${GA_CXX_FLAGS} -fno-alias")
  endif()
else(Clang)
  set(GA_CXX_FLAGS "${GA_CXX_FLAGS}  -std=c++11 -Wall -fPIC -g")
  if (CMAKE_BUILD_TYPE MATCHES Debug)
    set(GA_CXX_FLAGS "${GA_CXX_FLAGS}  -ggdb -O0")
  else()
    set(GA_CXX_FLAGS "${GA_CXX_FLAGS} -ffast-math -ftree-vectorize")
  endif()
endif()

#-------------------
# Includes
#-------------------

include_directories(inc)
if(Shark_FOUND)
  include_directories(${SHARK_INCLUDE_DIRS})
endif()
if(ENABLE_GEANTV)
  include_directories(${Vc_INCLUDE_DIR})
  include_directories(${VECGEOM_INCLUDE_DIR})
  include_directories(${HEPMC_INCLUDE_DIR})
  include_directories(${GEANTV_INCLUDE_DIR})
endif()

find_package(ROOT REQUIRED)
include(${ROOT_USE_FILE})

if(ENABLE_PERF)
  include_directories(${Perf_INCLUDE_DIR})
endif()
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif()
if(MLPACK_FOUND)
  include_directories(${ARMADILLO_INCLUDE_DIRS})
  include_directories(${MLPACK_INCLUDE_DIRS})
  include_directories(/usr/include/libxml2/)
endif()
if(Eigen_FOUND)
  include_directories(${EIGEN3_INCLUDE_DIR})
endif()
# Temporary solution - cmake MLPACK & Eigen3 doesnt work..
if(LIBCMAES_FOUND)
  include_directories(${LIBCMAES_INCLUDE_DIR})
endif()

#--------------------------------------------------
include_directories(/usr/include/libxml2/)
include_directories(/usr/include/eigen3)
# Mac OS X
include_directories(/usr/local/Cellar/eigen/3.2.6/include/eigen3/)
# Lib Instrument backward
include_directories(${CMAKE_SOURCE_DIR}/libinstrument/inc)

#-------------------
# Files
#-------------------

file(GLOB sources src/addstructures/*.cxx
                  src/algorithms/*.cxx
                  src/gaoperators/*.cxx
                  #src/hpc/*.cxx
                  src/generic/*.cxx
                  src/instrumentation/*.cxx
                  src/mva/*.cxx
                  src/output/*.cxx
                  src/problem/*.cxx
                  src/tools/*cxx
                  )
file(GLOB headers inc/addstructures/*.h
                  inc/algorithms/*.h
                  inc/gaoperators/*.h
                  inc/generic/*.h
                  #inc/hpc/*.h
                  inc/instrumentation/*.h
                  inc/mva/*.h
                  inc/output/*.h
                  inc/problem/*.h
                  inc/statistic/*.h
                  inc/tools/*.h
                  )

#-------------------
# Dictionary generation
#-------------------

ROOT_GENERATE_DICTIONARY(GADict ${headers} MODULE GA LINKDEF LinkDef.h)
ROOT_GENERATE_ROOTMAP(GADict LINKDEF LinkDef.h)

#-------------------
# Library
#-------------------

add_library(GA SHARED ${sources} GADict.cxx)
target_link_libraries(GA ${VC_LIBRARIES} ${ROOT_LIBRARIES} Hist RIO Tree Gui TreePlayer ${VC_LIBRARIES} ${VECGEOM_LIBRARIES} ${GEANTV_LIBRARIES} ${HEPMC_LIBRARIES} ${Boost_LIBRARY_DIRS} ${PerfMon_LIBRARY} ${MLPACK_LIBRARY})
set_target_properties(GA PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
install(TARGETS GA DESTINATION ${CMAKE_BINARY_DIR})
