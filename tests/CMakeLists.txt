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

if(ENABLE_GEANTV)
  include_directories(${GEANTV_INCLUDE_DIR})
  include_directories(${VECGEOM_INCLUDE_DIR})
  include_directories(${Vc_INCLUDE_DIR})
  include_directories(${HEPMC_INCLUDE_DIR})
endif()

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
          ${Boost_LIBRARY_DIRS} ${PerfMon_LIBRARY}
          ${MLPACK_LIBRARY}
          papi
          pfm
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
          ${Boost_LIBRARY_DIRS} ${PerfMon_LIBRARY}
          ${MLPACK_LIBRARY}
          papi
          pfm
          )
endif()

add_test(gagtests gagtests)
install(TARGETS gagtests DESTINATION ${CMAKE_BINARY_DIR})
