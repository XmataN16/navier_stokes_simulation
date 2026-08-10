cmake_minimum_required(VERSION 3.25)

project(
    navier_stokes_simulation
    VERSION 0.4.0
    DESCRIPTION "Navier-Stokes flow simulation framework"
    LANGUAGES CXX
)

# =============================================================================
# Project options
# =============================================================================

option(
    NSSIM_BUILD_TESTS
    "Build project tests"
    ON
)

option(
    NSSIM_BUILD_BENCHMARKS
    "Build project benchmarks"
    ON
)

# =============================================================================
# C++ configuration
# =============================================================================

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# =============================================================================
# Dependencies
# =============================================================================

find_package(
    PkgConfig
    REQUIRED
)

pkg_check_modules(
    FFTW3
    REQUIRED
    IMPORTED_TARGET
    fftw3
)

# =============================================================================
# Main library
# =============================================================================

add_library(
    nssim
    STATIC

    # Problem description
    src/problem/problem2d.cpp

    # I/O
    src/io/vti_writer2d.cpp
    src/io/vtk_time_series_writer2d.cpp

    # Simulation
    src/simulation/simulation2d.cpp

    # Pressure solvers
    src/solver/pressure/jacobi_pressure_solver2d.cpp
    src/solver/pressure/spectral_pressure_solver2d.cpp

    # Flow solvers
    src/solver/projection/cpu_projection_solver2d.cpp
)

add_library(
    nssim::nssim
    ALIAS
    nssim
)

target_include_directories(
    nssim

    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_compile_features(
    nssim

    PUBLIC
        cxx_std_23
)

target_link_libraries(
    nssim

    PRIVATE
        PkgConfig::FFTW3
)

# =============================================================================
# Compiler warnings
# =============================================================================

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

# =============================================================================
# Main executable
# =============================================================================

add_executable(
    nssim_app

    app/main.cpp
)

target_link_libraries(
    nssim_app

    PRIVATE
        nssim::nssim
)

# =============================================================================
# Tests
# =============================================================================

if(NSSIM_BUILD_TESTS)

    include(CTest)

    # -------------------------------------------------------------------------
    # Core smoke test
    # -------------------------------------------------------------------------

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

    # -------------------------------------------------------------------------
    # Projection smoke test
    # -------------------------------------------------------------------------

    add_executable(
        nssim_projection_smoke

        tests/projection_smoke.cpp
    )

    target_link_libraries(
        nssim_projection_smoke

        PRIVATE
            nssim::nssim
    )

    add_test(
        NAME projection_smoke
        COMMAND nssim_projection_smoke
    )

    # -------------------------------------------------------------------------
    # Simulation smoke test
    # -------------------------------------------------------------------------

    add_executable(
        nssim_simulation_smoke

        tests/simulation_smoke.cpp
    )

    target_link_libraries(
        nssim_simulation_smoke

        PRIVATE
            nssim::nssim
    )

    add_test(
        NAME simulation_smoke
        COMMAND nssim_simulation_smoke
    )

    # -------------------------------------------------------------------------
    # Spectral pressure solver test
    # -------------------------------------------------------------------------

    add_executable(
        nssim_spectral_pressure_test

        tests/spectral_pressure_solver_test.cpp
    )

    target_link_libraries(
        nssim_spectral_pressure_test

        PRIVATE
            nssim::nssim
    )

    add_test(
        NAME spectral_pressure
        COMMAND nssim_spectral_pressure_test
    )

endif()

# =============================================================================
# Benchmarks
# =============================================================================

if(NSSIM_BUILD_BENCHMARKS)

    add_executable(
        nssim_pressure_benchmark

        benchmarks/pressure_solver_benchmark.cpp
    )

    target_link_libraries(
        nssim_pressure_benchmark

        PRIVATE
            nssim::nssim
    )

endif()

# =============================================================================
# MinGW runtime
# =============================================================================

if(MINGW)

    target_link_options(
        nssim_app

        PRIVATE
            -static-libgcc
            -static-libstdc++
    )

    if(NSSIM_BUILD_TESTS)

        target_link_options(
            nssim_core_smoke

            PRIVATE
                -static-libgcc
                -static-libstdc++
        )

        target_link_options(
            nssim_projection_smoke

            PRIVATE
                -static-libgcc
                -static-libstdc++
        )

        target_link_options(
            nssim_simulation_smoke

            PRIVATE
                -static-libgcc
                -static-libstdc++
        )

        target_link_options(
            nssim_spectral_pressure_test

            PRIVATE
                -static-libgcc
                -static-libstdc++
        )

    endif()

    if(NSSIM_BUILD_BENCHMARKS)

        target_link_options(
            nssim_pressure_benchmark

            PRIVATE
                -static-libgcc
                -static-libstdc++
        )

    endif()

endif()