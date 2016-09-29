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
set(papiw_BINARY_DIR ${BINARY_DIR})
set(papiw_SOURCE_DIR ${SOURCE_DIR})
set(papiw_INCLUDE_DIR ${papiw_SOURCE_DIR})
include_directories(${papiw_INCLUDE_DIR})
set(papiw_LIBRARY ${papiw_BINARY_DIR}/lib/libpapi_wrap.a)

#################################################################################
# Cereal

ExternalProject_Add(
    cereal
    GIT_REPOSITORY https://github.com/USCiLab/cereal.git
    TIMEOUT 10
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/third_party/cereal
    STEP_TARGETS builds
    EXCLUDE_FROM_ALL TRUE
    )

ExternalProject_Get_Property(cereal BINARY_DIR)
ExternalProject_Get_Property(cereal SOURCE_DIR)
set(CEREAL_INCLUDE_DIR ${BINARY_DIR})
include_directories(${CEREAL_INCLUDE_DIR})

##################################################################################
# PFM

#ExternalProject_Add(
#    perf
#    GIT_REPOSITORY git://git.code.sf.net/p/perfmon2/libpfm4 perfmon2-libpfm4
#    #GIT_TAG v4.7.1
#    CONFIGURE_COMMAND <SOURCE_DIR>/configure
#                        --prefix=<INSTALL_DIR>
#    BUILD_IN_SOURCE 1
#    LOG_DOWNLOAD ON
#    LOG_BUILD ON)

#ExternalProject_Get_Property(perf BINARY_DIR)
#ExternalProject_Get_Property(perf SOURCE_DIR)
#set(PERF_BINARY_DIR ${BINARY_DIR})
#set(PERF_SOURCE_DIR ${SOURCE_DIR})
#set(PERF_INCLUDE_DIR ${PERF_SOURCE_DIR})
#include_directories(${PERF_INCLUDE_DIR})
#set(PERF_LIBRARY ${PERF_BINARY_DIR}/lib/libpfm4.a)

##################################################################################
# PAPI

#ExternalProject_Add(
#    papi
#    GIT_REPOSITORY https://icl.cs.utk.edu/git/papi.git
#    GIT_TAG papi-5-5-0-t
#    CONFIGURE_COMMAND cd <SOURCE_DIR>/src && ./configure --prefix=<INSTALL_DIR>
#    BUILD_COMMAND cd <SOURCE_DIR>/src && make -j4
#    BUILD_IN_SOURCE 1
#    LOG_DOWNLOAD ON
#    LOG_BUILD ON)

#ExternalProject_Get_Property(papi BINARY_DIR)
#ExternalProject_Get_Property(papi SOURCE_DIR)
#set(PAPI_BINARY_DIR ${BINARY_DIR})
#set(PAPI_SOURCE_DIR ${SOURCE_DIR})
#set(PAPI_INCLUDE_DIR ${PAPI_SOURCE_DIR}/include)
#include_directories(${PAPI_INCLUDE_DIR})
#set(PAPI_LIBRARY ${PAPI_BINARY_DIR}/lib/libpapi.a)
