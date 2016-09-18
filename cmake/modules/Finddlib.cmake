# Try to find DLIB
# Once done this will define
#
# DLIB_FOUND           - system has DLIB
# DLIB_INCLUDE_DIRS    - the DLIB include directories
# DLIB_LIBRARIES       - Link these to use DLIB

if($ENV{DLIB_ROOT})
  set(DLIB_ROOT $ENV{DLIB_ROOT})
endif()

find_package(dlib NO_MODULE PATHS ${DLIB_ROOT} )
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(dlib CONFIG_MODE)
