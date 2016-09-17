# Try to find DLIB
# Once done this will define
#
# DLIB_FOUND           - system has DLIB
# DLIB_INCLUDE_DIRS    - the DLIB include directories
# DLIB_LIBRARIES       - Link these to use DLIB

FIND_PATH(DLIB_INCLUDE_DIR clustering.h
	  PATHS
    ${CMAKE_MODULE_PATH}/include/dlib
    )

FIND_LIBRARY(DLIB_LIBRARY NAMES dlib PATHS
${CMAKE_MODULE_PATH}/lib/)

INCLUDE(FindPackageHandleStandardArgs)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(DLIB DEFAULT_MSG DLIB_LIBRARY DLIB_INCLUDE_DIR)

message(STATUS "DLIB_LIBRARY is ${DLIB_LIBRARY}")
message(STATUS "DLIB_INCLUDE_DIR is ${PDLIB_INCLUDE_DIR}")
