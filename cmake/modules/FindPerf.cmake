# - Locate PerfMon library
# Defines:
#
#  PerfMon_FOUND
#  PerfMon_INCLUDE_DIR
#  PerfMon_INCLUDE_DIRS (not cached)
#  PerfMon_LIBRARIES

if ($ENV{PerfMon_ROOT_DIR})
  set(PerfMon_ROOT_DIR $ENV{PerfMon_ROOT_DIR})
endif()

find_path(PerfMon_INCLUDE_DIR perfmon/pfmlib_perf_event.h 
          PATHS 
#         /usr/local/include/
          ${PerfMon_ROOT_DIR}/include/)
find_library(PerfMon_LIBRARY NAMES pfm 
             PATHS 
#            /usr/local/lib
             ${PerfMon_ROOT_DIR}/lib)

set(PerfMon_INCLUDE_DIR ${PerfMon_INCLUDE_DIR})
set(PerfMon_LIBRARY ${PerfMon_LIBRARY})

# handle the QUIETLY and REQUIRED arguments and set PerfMon_FOUND to TRUE if
# all listed variables are TRUE

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PerfMon DEFAULT_MSG PerfMon_INCLUDE_DIR PerfMon_LIBRARY)
mark_as_advanced(PERFMON_FOUND PerfMon_INCLUDE_DIR PerfMon_LIBRARY)