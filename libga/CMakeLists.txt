cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)
project(libga)

find_package(ROOT REQUIRED)
include(${ROOT_USE_FILE})
#----------------------------------------------------------------------------

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
#-----------------------------------------------------------------------------

if (GNU)
  set(GA_CXX_FLAGS "${GA_CXX_FLAGS} -Wall -fPIC -g")
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

include_directories(inc)
if(ENABLE_GEANTV)
  include_directories(${GEANTV_INCLUDE_DIR})
  include_directories(${VECGEOM_INCLUDE_DIR})
  include_directories(${Vc_INCLUDE_DIR})
endif()
if(ENABLE_PERFMON)
  include_directories(${PerfMon_INCLUDE_DIR})
endif()
#-----------------------------------------------------------------------------

file(GLOB sources src/*.cxx)
file(GLOB headers inc/PFMWatch.h
 inc/Functions.h
 inc/TGenes.h
 inc/Population.h
 inc/AlgorithmNSGA.h
 inc/HistogramManager.h
 inc/ExceptionMessenger.h
 inc/invoke_cpp11.hpp
 inc/type_traits.hpp
 inc/GeantVFitness.h
 )
#----------------------------------------------------------------------------
# Generate the Dictionary

ROOT_GENERATE_DICTIONARY(GADict ${headers} MODULE GA LINKDEF LinkDef.h)
ROOT_GENERATE_ROOTMAP(GADict LINKDEF LinkDef.h)

add_library(GA SHARED ${sources} GADict.cxx)
target_link_libraries(GA ${ROOT_LIBRARIES} Hist RIO Tree Gui)
#set_target_properties(GA PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
#install(TARGETS GA DESTINATION ${CMAKE_BINARY_DIR})
