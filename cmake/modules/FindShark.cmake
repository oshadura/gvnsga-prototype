# Locate the Shark library.
#
#
# This module defines the following variables:
# 
if ($ENV{Shark_DIR})
  set(Shark_DIR $ENV{Shark_DIR})
endif()
find_package(Shark ${SharkV_FIND_VERSION} NO_MODULE PATHS ${Shark_DIR} )
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Shark CONFIG_MODE)
