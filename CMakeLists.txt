cmake_minimum_required(VERSION 3.25)

project(
    navier_stokes_simulation
    VERSION 0.1.0
    LANGUAGES CXX
)

option(NSSIM_BUILD_TESTS "Build project tests" ON)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_library(nssim
    src/problem/problem2d.cpp
    src/io/vti_writer2d.cpp
)

add_library(nssim::nssim ALIAS nssim)

target_include_directories(
    nssim
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

if(MSVC)
    target_compile_options(
        nssim
        PRIVATE
            /W4
            /permissive-
    )
else()
    target_compile_options(
        nssim
        PRIVATE
            -Wall
            -Wextra
            -Wpedantic
            -Wshadow
            -Wconversion
    )
endif()

add_executable(
    nssim_app
    app/main.cpp
)

target_link_libraries(
    nssim_app
    PRIVATE
        nssim::nssim
)

if(NSSIM_BUILD_TESTS)
    include(CTest)

    add_executable(
        nssim_core_smoke
        tests/core_smoke.cpp
    )

    target_link_libraries(
        nssim_core_smoke
        PRIVATE
            nssim::nssim
    )

    add_test(
        NAME core_smoke
        COMMAND nssim_core_smoke
    )
endif()