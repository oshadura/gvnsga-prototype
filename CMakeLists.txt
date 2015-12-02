cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)

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
#-----------------------------------------------------------------------------

SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
SET(CMAKE_INSTALL_RPATH "")
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)
SET(CMAKE_INSTALL_RPATH ${CMAKE_BINARY_DIR})

# When CMake 2.8.12 is required, change this default to 1.
# When CMake 3.0.0 is required, remove this block (see CMP0042).
if (NOT DEFINED CMAKE_MACOSX_RPATH)
  set(CMAKE_MACOSX_RPATH 0)
endif()

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
