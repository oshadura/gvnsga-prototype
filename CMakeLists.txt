cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)

enable_language(CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED on)

#----------------------------------------------------------------------------
project(gvnsga-prototype)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "${CMAKE_BUILD_TYPE}: ${CMAKE_C_FLAGS_RELWITHDEBINFO}")

include(CheckCCompilerFlag)
include(CheckCXXCompilerFlag)

find_package(ROOT REQUIRED)
include(${ROOT_USE_FILE})

#----------------------------------------------------------------------------
# export OMP_NUM_THREADS=1
option(ENABLE_OPENMP "Adding OpenMP support" ON)
find_package(OpenMP)
if(OPENMP_FOUND)
  #set(GVNSGA_COMPILATION_FLAGS " ${GVNSGA_COMPILATION_FLAGS} ${GA_COMPILATION_FLAGS} -DENABLE_OPENMP")
  set(GA_CXX_FLAGS "${GA_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

#----------------------------------------------------------------------------
option(ENABLE_MPI "Adding MPI support" ON)
find_package(MPI)
if(MPI_CXX_FOUND)
  set(GVNSGA_COMPILATION_FLAGS "${GVNSGA_COMPILATION_FLAGS} ${GA_COMPILATION_FLAGS} -DENABLE_MPI")
endif()

#----------------------------------------------------------------------------
option(ENABLE_BOOST "Adding Boost support" ON)
find_package(Boost)
if(BOOST_FOUND)
  set(GVNSGA_COMPILATION_FLAGS " ${GVNSGA_COMPILATION_FLAGS} ${GA_COMPILATION_FLAGS} -DENABLE_BOOST")
endif()

#----------------------------------------------------------------------------
option(ENABLE_PERFMON "Adding PerfMon support" ON)

if(ENABLE_PERFMON)
  #find_package(PerfMon)
  #f(PERFMON_FOUND)
    set(GVNSGA_COMPILATION_FLAGS "${GVNSGA_COMPILATION_FLAGS} ${GA_COMPILATION_FLAGS} -DENABLE_PERFMON")
  #endif()
endif()
#----------------------------------------------------------------------------
option(ENABLE_GEANTV "We are using libGA already for GeantV" OFF)

if(ENABLE_GEANTV)
	find_package(GeantV)
	if(GeantV_FOUND)
  		set(GVNSGA_COMPILATION_FLAGS "${GVNSGA_COMPILATION_FLAGS} ${GA_COMPILATION_FLAGS} -DENABLE_GEANTV")
  		set(GA_CXX_FLAGS "${GA_CXX_FLAGS} ${GEANTV_CXX_FLAGS}")
	endif()
endif()
#----------------------------------------------------------------------------
if(ENABLE_GEANTV)
  set(ENABLE_VECGEOM ON)
  set(ENABLE_VC ON)
endif()

#----------------------------------------------------------------------------
option(NUMERIC_LIB "We are using libGA already for simple models" OFF)

if(NUMERIC_LIB)
  set(ENABLE_GEANTV OFF)
  set(ENABLE_VECGEOM OFF)
  set(ENABLE_VC OFF)
  set(GVNSGA_COMPILATION_FLAGS "${GVNSGA_COMPILATION_FLAGS} ${GA_COMPILATION_FLAGS} -DNUMERIC_LIB")
endif()

#----------------------------------------------------------------------------
option(ENABLE_VECGEOM "If we are using libGA already with VecGeom" OFF)

if(ENABLE_VECGEOM)
	find_package(VecGeom)
	if(VecGeom_FOUND)
  		set(GVNSGA_COMPILATION_FLAGS "${GVNSGA_COMPILATION_FLAGS} ${GA_COMPILATION_FLAGS} -DENABLE_VECGEOM")
  		set(GA_CXX_FLAGS "${GA_CXX_FLAGS} ${VECGEOM_CXX_FLAGS}")
	endif()
endif()
#----------------------------------------------------------------------------

option(ENABLE_VC "If we are using libGA already with Vc" OFF)

if(ENABLE_VC)
	find_package(Vc)
	if(Vc_FOUND)
  		set(GVNSGA_COMPILATION_FLAGS "${GVNSGA_COMPILATION_FLAGS} ${GA_COMPILATION_FLAGS} -DENABLE_VC")
  		set(GA_CXX_FLAGS "${GA_CXX_FLAGS} ${VC_CXX_FLAGS}")
		endif()
endif()
#----------------------------------------------------------------------------

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GA_CXX_FLAGS} ${GVNSGA_COMPILATION_FLAGS}")
message(STATUS "Compiling GA with C++ flags: ${CMAKE_CXX_FLAGS}")
#-----------------------------------------------------------------------------

SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
SET(CMAKE_INSTALL_RPATH "")
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)
SET(CMAKE_INSTALL_RPATH ${CMAKE_BINARY_DIR})
#------------------Backend detection from CMake options------------------------
function (getListOfVarsWith _prefix _varResult)
	set (_variableName "${_prefix}")
    string (REGEX MATCH "${_variableName}" _varResult "${VECGEOM_DEFINITIONS}")
endfunction()

#----------------------------Debug macros-------------------------------------
macro(debug msg)
	message(STATUS "Debug ${msg}")
endmacro()

macro(debugValue _variableName)
	debug("${_variableName}=\${${_variableName}}")
endmacro()
#-----------------------------------------------------------------------------
getListOfVarsWith(CUDA BACKEND)
if(BACKEND MATCHES "")
	getListOfVarsWith(mic BACKEND)
	if(BACKEND MATCHES "")
		set(BACKEND "CPU")
		debugValue(BACKEND)
	endif()
endif()

message(STATUS "Check which is backend of VecGeom/GeantV: ${BACKEND} backend")
#----------------------------------------------------------------------------

add_subdirectory(libga)
add_subdirectory(example)
#----------------------------------------------------------------------------
