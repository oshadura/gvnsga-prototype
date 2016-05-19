project(gagtests)
include(ExternalProject)

ExternalProject_Add(googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    CMAKE_ARGS -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG:PATH=DebugLibs
               -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE:PATH=ReleaseLibs
               -Dgtest_force_shared_crt=ON
               -DBUILD_GTEST=ON
    PREFIX "${CMAKE_CURRENT_BINARY_DIR}"
    INSTALL_COMMAND ""
    LOG_DOWNLOAD ON
    LOG_CONFIGURE ON
    LOG_BUILD ON
)


#-------------------
# Specify includes
#-------------------
ExternalProject_Get_Property(googletest source_dir)
message(STATUS "DEBUG: source_dir for GoogleTest = ${source_dir}")
set(GTEST_INCLUDE_DIRS ${source_dir}/googletest/include)
set(GMOCK_INCLUDE_DIRS ${source_dir}/googlemock/include)

#-------------------
# Specify libraries
#-------------------
ExternalProject_Get_Property(googletest binary_dir)
message(STATUS "DEBUG: binary_dir for GoogleTest = ${binary_dir}")
set(GTEST_LIBS_DIR ${binary_dir}/googlemock/gtest)
set(GMOCK_LIBS_DIR ${binary_dir}/googlemock/)

#-------------------
# Includes
#-------------------

include_directories(${GTEST_INCLUDE_DIRS} ${GMOCK_INCLUDE_DIRS})
include_directories(inc)
include_directories(${CMAKE_SOURCE_DIR}/libga/inc)
include_directories(/usr/include/eigen3)
# Mac OS X
include_directories(/usr/local/Cellar/eigen/3.2.6/include/eigen3/)

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
          # Fix this
          ${MLPACK_LIBRARY}
          )
elseif(UNIX)
  target_link_libraries(gagtests
          ${GTEST_LIBS_DIR}/libgtest.a
          ${GTEST_LIBS_DIR}/libgtest_main.a
          ${GMOCK_LIBS_DIR}/libgmock.a
          ${GMOCK_LIBS_DIR}/libgmock_main.a
          GA
          Instrument
          # Fix this
          ${MLPACK_LIBRARY}
          )
endif()

add_test(gagtests gagtests)
install(TARGETS gagtests DESTINATION ${CMAKE_BINARY_DIR})
