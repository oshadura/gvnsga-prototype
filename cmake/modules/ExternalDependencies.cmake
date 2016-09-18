include(ExternalProject)

###############################################################################
# Google Test
ExternalProject_Add(
    gtest
    GIT_REPOSITORY https://github.com/google/googletest.git
    TIMEOUT 10
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/third_party/gtest
    # no install required, we link the library from the build tree
    INSTALL_COMMAND ""
    # Wrap download, configure and build steps in a script to log output
    LOG_DOWNLOAD ON
    LOG_BUILD ON)

ExternalProject_Get_Property(gtest BINARY_DIR)
ExternalProject_Get_Property(gtest SOURCE_DIR)
set(gtest_BINARY_DIR ${BINARY_DIR})
set(gtest_SOURCE_DIR ${SOURCE_DIR})
set(gtest_INCLUDE_DIR ${gtest_SOURCE_DIR}/googletest/include)
include_directories(${gtest_INCLUDE_DIR})
set(gtest_LIBRARY ${gtest_BINARY_DIR}/googlemock/gtest/libgtest.a)
set(gtest_MAIN_LIBRARY ${gtest_BINARY_DIR}/googlemock/gtest/libgtest_main.a)

set(gmock_INCLUDE_DIR ${gtest_SOURCE_DIR}/googlemock/include)
include_directories(${gmock_INCLUDE_DIR})
set(gmock_LIBRARY ${gtest_BINARY_DIR}/googlemock/libgmock.a)
set(gmock_MAIN_LIBRARY ${gtest_BINARY_DIR}/googlemock/libgmock_main.a)


################################################################################
# PAPI-WRAP
ExternalProject_Add(
	papi-wrap
	GIT_REPOSITORY https://github.com/bcumming/papi-wrap.git
	TIMEOUT 10
	PREFIX ${CMAKE_CURRENT_BINARY_DIR}/third_party/papi-wrap
	CMAKE_ARGS
	INSTALL_COMMAND ""
	LOG_DOWNLOAD ON
	LOG_BUILD ON)

ExternalProject_Get_Property(papi-wrap BINARY_DIR)
ExternalProject_Get_Property(papi-wrap SOURCE_DIR)
set(papi_BINARY_DIR ${BINARY_DIR})
set(papi_SOURCE_DIR ${SOURCE_DIR})
set(papi_INCLUDE_DIR ${papi_SOURCE_DIR})
include_directories(${papi_INCLUDE_DIR})
set(gtest_LIBRARY ${papi_BINARY_DIR}/lib/libpapi_wrap.a)
