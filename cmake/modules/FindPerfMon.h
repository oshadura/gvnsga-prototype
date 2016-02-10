# - Locate PerfMon library
# Defines:
#
#  PerfMon_FOUND
#  PerfMon_INCLUDE_DIR
#  PerfMon_INCLUDE_DIRS (not cached)
#  PerfMon_LIBRARIES

if ($ENV{PerfMon_DIR})
  set(PerfMon_DIR $ENV{PerfMon_DIR})
endif()

find_path(PerfMon_INCLUDE_DIR perfmon/pfmlib_perf_event.h 
          HINTS /usr/local/include ${PerfMon_DIR}/include $ENV{PerfMon_DIR}/include)
find_library(PerfMon_LIBRARY NAMES pfm 
             HINTS /usr/local/lib ${PerfMon_DIR}/lib $ENV{PerfMon_DIR}/lib)
set(PerfMon_INCLUDE_DIRS ${PerfMon_INCLUDE_DIR})
set(PerfMon_LIBRARIES ${PerfMon_LIBRARY})


# handle the QUIETLY and REQUIRED arguments and set PerfMon_FOUND to TRUE if
# all listed variables are TRUE

include(FindPackageHandleStandardArgs)

FindPackageHandleStandardArgs(PerfMon DEFAULT_MSG PerfMon_INCLUDE_DIR PerfMon_LIBRARY)
mark_as_advanced(PerfMon_FOUND PerfMon_INCLUDE_DIR PerfMon_LIBRARY)