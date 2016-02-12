# - Locate HepMC library
# Defines:
#
#  HEPMC_FOUND
#  HEPMC_INCLUDE_DIR
#  HEPMC_INCLUDE_DIRS (not cached)
#  HEPMC_LIBRARIES

find_path(HEPMC_INCLUDE_DIR HepMC/GenEvent.h 
          HINTS ${HEPMC_ROOT_DIR}/include $ENV{HEPMC_ROOT_DIR}/include)
find_library(HEPMC_LIBRARY NAMES HepMC 
             HINTS ${HEPMC_ROOT_DIR}/lib64 $ENV{HEPMC_ROOT_DIR}/lib64)

set(HEPMC_INCLUDE_DIRS ${HEPMC_INCLUDE_DIR})
set(HEPMC_LIBRARIES ${HEPMC_LIBRARY})


# handle the QUIETLY and REQUIRED arguments and set HEPMC_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HepMC DEFAULT_MSG HEPMC_INCLUDE_DIR HEPMC_LIBRARY)
mark_as_advanced(HEPMC_FOUND HEPMC_INCLUDE_DIR HEPMC_LIBRARY)