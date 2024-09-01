############################################################
# Project
############################################################

cmake_minimum_required(VERSION 3.29)

project(shockwave)

list(APPEND CMAKE_MODULE_PATH "cmake")
include(macros)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(FetchContent)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

############################################################
# Options
############################################################

option(SHOCKWAVE_BUILD_TESTS "Enables the test targets" OFF)
option(SHOCKWAVE_TEST_COVERAGE "Enables the test coverage targets" OFF)

############################################################
# External libs
############################################################

# add_subdirectory(libs)

############################################################
# Shockwave
############################################################

shockwave_add_library(shockwave "")
add_library(shockwave::shockwave ALIAS shockwave)

target_include_directories(
    shockwave PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
                     $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

add_subdirectory(include/)
add_subdirectory(src/)

generates_compile_commands(shockwave)

############################################################
# Tests
############################################################

if(SHOCKWAVE_BUILD_TESTS)
    FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG 31588bb # V3.7.0
    )

    FetchContent_MakeAvailable(Catch2)
    list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)

    # hack to prevent CTest targets from showing up
    # https://stackoverflow.com/a/57240389
    set_property(GLOBAL PROPERTY CTEST_TARGETS_ADDED 1)
    include(CTest)

    include(Catch)

    enable_testing()

    add_executable(shockwave-test "")
    target_link_libraries(
        shockwave-test PRIVATE Catch2::Catch2WithMain shockwave
    )

    catch_discover_tests(shockwave-test)

    generates_compile_commands(shockwave-test)

    add_subdirectory(test)

    if(SHOCKWAVE_TEST_COVERAGE)
        shockwave_coverage(shockwave-test)
    endif()
endif()

############################################################
# Docs
############################################################

# TODO:

############################################################
# Installation
############################################################

install(
    TARGETS shockwave
    EXPORT shockwave
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} FILE_SET HEADERS
)

install(
    EXPORT shockwave
    NAMESPACE shockwave::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/shockwave
)

configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/config.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/shockwave-config.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/shockwave
)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/shockwave-config.cmake"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/shockwave
)
