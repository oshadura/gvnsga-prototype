find_path(LIBCMAES_INCLUDE_DIR libcmaes/cmaes.h 
          HINTS ${LIBCMAES_ROOT_DIR}/include $ENV{LIBCMAES_ROOT_DIR}/include)
find_library(LIBCMAES_LIBRARY NAMES cmaes 
             HINTS ${LIBCMAES_ROOT_DIR}/lib $ENV{LIBCMAES_ROOT_DIR}/lib)

set(LIBCMAES_INCLUDE_DIRS ${LIBCMAES_INCLUDE_DIR})
set(LIBCMAES_LIBRARIES ${LIBCMAES_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Libcmaes DEFAULT_MSG LIBCMAES_INCLUDE_DIR LIBCMAES_LIBRARY)
mark_as_advanced(LIBCMAES_FOUND LIBCMAES_INCLUDE_DIR LIBCMAES_LIBRARY)