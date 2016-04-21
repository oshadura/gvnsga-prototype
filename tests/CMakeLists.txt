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

# Specify include dir
ExternalProject_Get_Property(googletest source_dir)
message(STATUS "DEBUG: source_dir for GoogleTest = ${source_dir}")
set(GTEST_INCLUDE_DIRS ${source_dir}/googletest/include)

# Specify MainTest's link libraries
ExternalProject_Get_Property(googletest binary_dir)
message(STATUS "DEBUG: binary_dir for GoogleTest = ${binary_dir}")
set(GTEST_LIBS_DIR ${binary_dir}/googlemock/gtest)

#-------------------
# Includes 
#-------------------

include_directories(${GTEST_INCLUDE_DIRS})

#-------------------
# Testing GA library
#-------------------
enable_testing()
set(PROJECT_TEST_NAME ${PROJECT_NAME_STR}_test)
file(GLOB TEST_SRC_FILES ${CMAKE_SOURCE_DIR}/tests/*.cxx)
add_executable(${PROJECT_TEST_NAME} ${TEST_SRC_FILES})
add_dependencies(${PROJECT_TEST_NAME} googletest)

target_link_libraries(${PROJECT_TEST_NAME}
        ${GTEST_LIBS_DIR}/libgtest.a
        ${GTEST_LIBS_DIR}/libgtest_main.a)

add_test(csvtest ${PROJECT_TEST_NAME})