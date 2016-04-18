if ($ENV{LIBCMAES_ROOT_DIR})
  set(LIBCMAES_ROOT_DIR $ENV{LIBCMAES_ROOT_DIR})
endif()

if (LIBCMAES_INCLUDE_DIR AND LIBCMAES_LIBRARIES)
   # in cache already
   SET(Libcmaes_FIND_QUIETLY TRUE)
endif (LIBCMAES_INCLUDE_DIR AND LIBCMAES_LIBRARIES)

find_path(LIBCMAES_INCLUDE_DIR NAMES cmaes.h
  PATHS
  ${LIBCMAES_ROOT_DIR}/include/libcmaes
)

find_library(LIBCMAES_LIBRARIES NAMES libcmaes
  PATHS
   ${LIBCMAES_ROOT_DIR}/lib/
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Libcmaes DEFAULT_MSG LIBCMAES_INCLUDE_DIR LIBCMAES_LIBRARIES )

# show the LIBCMAES_INCLUDE_DIR and LIBCMAES_LIBRARIES variables only in the advanced view
mark_as_advanced(LIBCMAES_INCLUDE_DIR LIBCMAES_LIBRARIES)