# Locate the GeantV library. 
# 
#
# This module defines the following variables:
# GEANTV_FOUND
# GEANTV_INCLUDE_DIR
# GEANTV_LIBRARIES
# GEANTV_DEFINITIONS
# GEANTV_VERSION_MAJOR # not yet
# GEANTV_VERSION_MINOR # not yet
# GEANTV_VERSION_PATCH # not yet
# GEANTV_VERSION # not yet
# GEANTV_VERSION_STRING # not yet
# GEANTV_INSTALL_DIR
# GEANTV_LIB_DIR
# GEANTV_CMAKE_MODULES_DIR
#
include(FindPackageHandleStandardArgs)

if ($ENV{GeantV_DIR})
  set(GeantV_DIR $ENV{GeantV_DIR})
endif()
find_package(GeantV ${GeantV_FIND_VERSION} NO_MODULE PATHS ${GeantV_DIR} )

find_package_handle_standard_args(GeantV CONFIG_MODE)
