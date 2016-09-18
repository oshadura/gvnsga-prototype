#----------------------------------------------------------------------------
# Setup the project

cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
#----------------------------------------------------------------------------
#project(libinstrument)
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
option(ENABLE_NUMA "Adding NUMA support" OFF)
find_package(Numa QUIET)
if(ENABLE_NUMA)
  set(INSTR_COMPILATION_FLAGS "${INSTR_COMPILATION_FLAGS} -DNUMA")
endif()
#----------------------------------------------------------------------------
option(ENABLE_PIN "Adding Pin benchmarking support" OFF)
find_package(Pin QUIET)
if(ENABLE_PIN)
  set(INSTR_COMPILATION_FLAGS "${INSTR_COMPILATION_FLAGS} -DPIN")
endif()
#----------------------------------------------------------------------------
option(ENABLE_GPERF "Adding GPERF benchmarking support" OFF)
find_package(GooglePerftools QUIET)
if(ENABLE_GPERF)
  set(INSTR_COMPILATION_FLAGS "${INSTR_COMPILATION_FLAGS} -DGPERF")
endif()
#----------------------------------------------------------------------------
option(ENABLE_PAPI "Adding PAPI benchmarking support" ON)
find_package(PAPI QUIET)
if(ENABLE_PAPI)
  set(INSTR_COMPILATION_FLAGS "${INSTR_COMPILATION_FLAGS} -DPAPI")
endif()
#----------------------------------------------------------------------------
option(ENABLE_PERF "Adding Perfmon benchmarking support" ON)
find_package(PERF QUIET)
if(ENABLE_PERF)
  set(INSTR_COMPILATION_FLAGS "${INSTR_COMPILATION_FLAGS} -DPERF")
endif()
#----------------------------------------------------------------------------

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${_build_type}} ${INSTR_CXX_FLAGS} ${INSTR_COMPILATION_FLAGS}")
message(STATUS "Compiling with C++ flags: ${CMAKE_CXX_FLAGS}")
#--------------------------------------------------------------------------

include_directories(inc)
file(GLOB sources src/*.cxx)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -W -Wall -ansi -ldl")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -W -Wall -ansi -ldl")
add_library(Instrument SHARED ${sources})
if(PAPI_FOUND)
  target_link_libraries(Instrument ${PAPI_LIBRARY})
else()
  target_link_libraries(Instrument)
endif()
set_target_properties(Instrument PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
install(TARGETS Instrument DESTINATION ${CMAKE_BINARY_DIR})

message(STATUS "Compiling Instruments with C++ flags: ${CMAKE_CXX_FLAGS}")