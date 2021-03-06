# Try to find MLPACK
# Once done this will define
#
# MLPACK_FOUND           - system has MLPACK
# MLPACK_INCLUDE_DIRS    - the MLPACK include directories
# MLPACK_LIBRARIES       - Link these to use MLPACK

if ($ENV{MLPACK_ROOT_DIR})
  set(MLPACK_ROOT_DIR $ENV{MLPACK_ROOT_DIR})
endif()

#IF (MLPACK_INCLUDE_DIRS)
  # Already in cache, be silent
#  SET(MLPACK_FIND_QUIETLY TRUE)
#ENDIF (MLPACK_INCLUDE_DIRS)

FIND_PATH(MLPACK_INCLUDE_DIR core.hpp
	  PATHS
#    /usr/local/include/mlpack
    /usr/include/mlpack
    ${MLPACK_ROOT_DIR}/include/mlpack
    )

FIND_LIBRARY(MLPACK_LIBRARY NAMES mlpack PATHS 
  /usr/lib64/
  ${MLPACK_ROOT_DIR}/lib64/)

INCLUDE(FindPackageHandleStandardArgs)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(mlpack DEFAULT_MSG MLPACK_LIBRARY MLPACK_INCLUDE_DIR)

message(STATUS "MLPACK_LIBRARY is ${MLPACK_LIBRARY}")
message(STATUS "MLPACK_INCLUDE_DIR is ${PMLPACK_INCLUDE_DIR}")
