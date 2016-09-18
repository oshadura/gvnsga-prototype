# - Locate PerfMon library
# Defines:
#
#  PerfMon_FOUND
#  PerfMon_INCLUDE_DIR
#  PerfMon_INCLUDE_DIRS (not cached)
#  PerfMon_LIBRARIES

if($ENV{PERF_ROOT})
  set(PERF_ROOT $ENV{PERF_ROOT})
endif()

find_path(PERF_INCLUDE_DIR perfmon/pfmlib_perf_event.h 
          PATHS 
#         /usr/local/include/
          ${PERF_ROOT}/include/)
find_library(PERF_LIBRARY NAMES pfm 
             PATHS 
#            /usr/local/lib
             ${PERF_ROOT}/lib)

set(PERF_INCLUDE_DIR ${Perf_INCLUDE_DIR})
set(PERF_LIBRARY ${Perf_LIBRARY})

# handle the QUIETLY and REQUIRED arguments and set PerfMon_FOUND to TRUE if
# all listed variables are TRUE

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PERF DEFAULT_MSG PERF_INCLUDE_DIR PERF_LIBRARY)
mark_as_advanced(PERF_FOUND PERF_INCLUDE_DIR PERF_LIBRARY)