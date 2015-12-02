cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)
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
#----------------------------------------------------------------------------
############
if (GNU)
  set(GA_CXX_FLAGS "${GA_CXX_FLAGS} -std=c++11 -Wall -fPIC -g")
  if (CMAKE_BUILD_TYPE MATCHES Debug)
    set(GA_CXX_FLAGS "${GA_CXX_FLAGS} -ggdb -O0")
  else()
    set(GA_CXX_FLAGS "${GA_CXX_FLAGS} -finline-limit=10000000 -ffast-math -ftree-vectorize")
  endif()
  if (Vc)
    set(GA_CXX_FLAGS "${GA_CXX_FLAGS} -fabi-version=6")
  endif()
elseif(Intel)
  set(GA_CXX_FLAGS "${GA_CXX_FLAGS} -std=c++11 -Wall -fPIC -diag-disable 3438")
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

#----------------------------------------------------------------------------

option(ENABLE_OPENMP "Adding OpenMP support" ON)
find_package(OpenMP)
if(ENABLE_OPENMP)
  set(GVNSGA_COMPILATION_FLAGS "${GA_COMPILATION_FLAGS} -DOPENMP")
endif()
#----------------------------------------------------------------------------
#option(ENABLE_BOOST "Adding OpenMP support" ON)
#find_package(Boost)
#if(ENABLE_BOOST)
#  set(GVNSGA_COMPILATION_FLAGS "${GA_COMPILATION_FLAGS} -DBOOST")
#endif()
#----------------------------------------------------------------------------
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GA_CXX_FLAGS} ${GA_COMPILATION_FLAGS}")
message(STATUS "Compiling GA with C++ flags: ${CMAKE_CXX_FLAGS}")
#----------------------------------------------------------------------------

include_directories(inc)
#-----------------------------------------------------------------------------

file(GLOB sources src/*.cxx)
file(GLOB headers inc/Functions.h inc/Genes.h inc/Population.h inc/AlgorithmNSGA.h inc/HistogramManager.h)
#----------------------------------------------------------------------------
# Generate the Dictionary

ROOT_GENERATE_DICTIONARY(GADict ${headers} MODULE GA LINKDEF LinkDef.h)
ROOT_GENERATE_ROOTMAP(GADict LINKDEF LinkDef.h)

add_library(GA SHARED ${sources} GADict.cxx)
target_link_libraries(GA ${ROOT_LIBRARIES} Hist RIO Tree Gui)
set_target_properties(GA PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
install(TARGETS GA DESTINATION ${CMAKE_BINARY_DIR})
