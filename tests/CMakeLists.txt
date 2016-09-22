include(ExternalProject)

###############################################################################
# Google Test
#-------------------
# Includes
#-------------------

include_directories(${gtest_INCLUDE_DIR} ${gmock_INCLUDE_DIR})
include_directories(inc)
include_directories(${CMAKE_SOURCE_DIR}/libga/inc)
include_directories(${CMAKE_SOURCE_DIR}/libinstrument/inc)

# Linux
include_directories(/usr/local/include/eigen3)
# Mac OS X
include_directories(/usr/local/Cellar/eigen/3.2.6/include/eigen3/)

# DLIB
find_package(dlib REQUIRED)
include_directories(${dlib_INCLUDE_DIRS})

# STUPUD TEST
include_directories(/home/geant/Install/perfmon/include)
include_directories(/home/geant/Install/dlib/include)


if(ENABLE_GEANTV)
  include_directories(${GEANTV_INCLUDE_DIR})
  include_directories(${VECGEOM_INCLUDE_DIR})
  include_directories(BEFORE SYSTEM ${Vc_INCLUDE_DIR})
  include_directories(${HEPMC_INCLUDE_DIR})
endif()

include_directories(${CMAKE_BINARY_DIR}/third_party/papi-wrap/src)

#-------------------
# Setting up Ctest also
#-------------------
include(CTest)
enable_testing()

#-------------------
# Testing GA library
#-------------------

file(GLOB_RECURSE TEST_FILES "*.cxx")
add_executable(gagtests ${TEST_FILES})

if(APPLE)
  target_link_libraries(gagtests
          /usr/local/lib/libgtest.a
          /usr/local/lib/libgtest_main.a
          /usr/local/lib/libgmock.a
          /usr/local/lib/libgmock_main.a
          GA
          Instrument
          ${ROOT_LIBRARIES}
          Hist
          RIO
          Tree
          Gui
          Geom
          ${GEANTV_LIBRARIES}
          ${VECGEOM_LIBRARIES}
          ${VC_LIBRARIES} ${HEPMC_LIBRARIES}
          ${Boost_LIBRARY_DIRS} ${PERF_LIBRARY}
          ${MLPACK_LIBRARY}
          ${dlib_LIBRARIES}
          papi
          pfm
          lapack
          hwloc
# STUPID TEST
#          /home/geant/Install/perfmon/lib/libpfm.so
          )
elseif(UNIX)
  target_link_libraries(gagtests
          ${gmock_LIBRARY}
          ${gtest_LIBRARY}
          GA
          Instrument
          ${ROOT_LIBRARIES}
          Hist
          RIO
          Tree
          Gui
          Geom
          ${GEANTV_LIBRARIES}
          ${VECGEOM_LIBRARIES}
          ${VC_LIBRARIES} ${HEPMC_LIBRARIES}
          ${Boost_LIBRARY_DIRS} ${PERF_LIBRARY}
          ${MLPACK_LIBRARY}
          ${dlib_LIBRARIES}
          papi
#          pfm
          lapack
          hwloc
# STUPID TEST
          /home/geant/Install/perfmon/lib/libpfm.so
          ${CMAKE_BINARY_DIR}/third_party/papi-wrap/src/papi-wrap-build/lib/libpapi_wrap.a
          )
endif()

add_test(gagtests gagtests)
install(TARGETS gagtests DESTINATION ${CMAKE_BINARY_DIR})
