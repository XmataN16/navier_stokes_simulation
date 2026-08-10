#include "nssim/core/field2d.hpp"
#include "nssim/core/grid2d.hpp"

#include "nssim/problem/boundary_conditions2d.hpp"

#include "nssim/solver/pressure/jacobi_pressure_solver2d.hpp"
#include "nssim/solver/pressure/spectral_pressure_solver2d.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using namespace nssim;

struct BenchmarkProblem final {
    UniformGrid2D grid;

    Field2D exact_pressure;
    Field2D rhs;

    explicit BenchmarkProblem(
        const std::size_t n
    )
        : grid{
              n,
              n,
              1.0,
              1.0
          },

          exact_pressure{
              n,
              n
          },

          rhs{
              n,
              n
          } {

        initialize();
    }

private:
    void initialize() {

        const auto nx =
            grid.nx_cells();

        const auto ny =
            grid.ny_cells();

        /*
         * Low-frequency mode is deliberately
         * selected because this is where Jacobi
         * converges slowly.
         */
        constexpr std::size_t mode_x =
            1;

        constexpr std::size_t mode_y =
            1;

        const Real pi =
            std::numbers::pi_v<Real>;

        const Real lambda_x =
            -Real{4} /
            (
                grid.dx() *
                grid.dx()
            ) *
            std::pow(
                std::sin(
                    pi *
                    static_cast<Real>(
                        mode_x
                    ) /
                    (
                        Real{2} *
                        static_cast<Real>(
                            nx
                        )
                    )
                ),
                Real{2}
            );

        const Real lambda_y =
            -Real{4} /
            (
                grid.dy() *
                grid.dy()
            ) *
            std::pow(
                std::sin(
                    pi *
                    static_cast<Real>(
                        mode_y
                    ) /
                    (
                        Real{2} *
                        static_cast<Real>(
                            ny
                        )
                    )
                ),
                Real{2}
            );

        const Real lambda =
            lambda_x +
            lambda_y;

        for (
            std::size_t j = 0;
            j < ny;
            ++j
        ) {
            for (
                std::size_t i = 0;
                i < nx;
                ++i
            ) {

                const Real x_mode =
                    std::cos(
                        pi *
                        static_cast<Real>(
                            mode_x
                        ) *
                        (
                            static_cast<Real>(
                                i
                            ) +
                            Real{0.5}
                        ) /
                        static_cast<Real>(
                            nx
                        )
                    );

                const Real y_mode =
                    std::cos(
                        pi *
                        static_cast<Real>(
                            mode_y
                        ) *
                        (
                            static_cast<Real>(
                                j
                            ) +
                            Real{0.5}
                        ) /
                        static_cast<Real>(
                            ny
                        )
                    );

                exact_pressure(i, j) =
                    x_mode *
                    y_mode;

                rhs(i, j) =
                    lambda *
                    exact_pressure(i, j);
            }
        }
    }
};

struct Measurement final {
    double milliseconds{};

    std::size_t iterations{};

    Real residual_l2{};
    Real error_l2{};

    bool converged{};
};

[[nodiscard]]
Real compute_error_l2(
    const Field2D& numerical,
    const Field2D& exact
) {

    Real squared_sum{};

    for (
        std::size_t j = 0;
        j < exact.ny();
        ++j
    ) {
        for (
            std::size_t i = 0;
            i < exact.nx();
            ++i
        ) {

            const Real difference =
                numerical(i, j) -
                exact(i, j);

            squared_sum +=
                difference *
                difference;
        }
    }

    return std::sqrt(
        squared_sum /
        static_cast<Real>(
            exact.size()
        )
    );
}

[[nodiscard]]
Measurement measure_solver(
    IPressureSolver2D& solver,
    const BenchmarkProblem& problem,
    const PressureSolverConfig& config
) {

    Field2D pressure{
        problem.grid.nx_cells(),
        problem.grid.ny_cells()
    };

    const auto start =
        std::chrono::
            steady_clock::now();

    const PressureSolveResult result =
        solver.solve(
            problem.grid,
            BoundarySet2D{},
            problem.rhs,
            pressure,
            config
        );

    const auto finish =
        std::chrono::
            steady_clock::now();

    const std::chrono::duration<double, std::milli>
        elapsed =
            finish -
            start;

    return {
        .milliseconds =
            elapsed.count(),

        .iterations =
            result.iterations,

        .residual_l2 =
            result.residual_l2,

        .error_l2 =
            compute_error_l2(
                pressure,
                problem.exact_pressure
            ),

        .converged =
            result.converged
    };
}

void print_measurement(
    const std::size_t n,
    const std::string& solver_name,
    const Measurement& measurement
) {

    std::cout
        << std::left
        << std::setw(10)
        << (
            std::to_string(n) +
            "x" +
            std::to_string(n)
        )

        << std::setw(24)
        << solver_name

        << std::right
        << std::setw(14)
        << std::fixed
        << std::setprecision(3)
        << measurement.milliseconds

        << std::setw(14)
        << measurement.iterations

        << std::setw(18)
        << std::scientific
        << measurement.residual_l2

        << std::setw(18)
        << measurement.error_l2

        << std::setw(12)
        << (
            measurement.converged
            ? "yes"
            : "no"
        )

        << '\n';
}

[[nodiscard]]
std::vector<std::size_t>
parse_grid_sizes(
    const int argc,
    char** argv
) {

    if (argc <= 1) {
        return {
            32,
            64,
            128
        };
    }

    std::vector<std::size_t>
        result;

    for (
        int argument = 1;
        argument < argc;
        ++argument
    ) {

        const auto value =
            std::stoull(
                argv[argument]
            );

        if (value < 4) {
            throw std::invalid_argument{
                "Benchmark grid size must "
                "be at least 4"
            };
        }

        result.push_back(
            static_cast<std::size_t>(
                value
            )
        );
    }

    return result;
}

} // namespace

int main(
    const int argc,
    char** argv
) {

    try {

        using namespace nssim;

        const auto sizes =
            parse_grid_sizes(
                argc,
                argv
            );

        std::cout
            << "Pressure solver benchmark\n\n"

            << std::left
            << std::setw(10)
            << "Grid"

            << std::setw(24)
            << "Solver"

            << std::right
            << std::setw(14)
            << "Time [ms]"

            << std::setw(14)
            << "Iterations"

            << std::setw(18)
            << "Residual L2"

            << std::setw(18)
            << "Error L2"

            << std::setw(12)
            << "Converged"

            << '\n';

        std::cout
            << std::string(
                110,
                '-'
            )
            << '\n';

        for (
            const std::size_t n :
            sizes
        ) {

            const BenchmarkProblem
                problem{
                    n
                };

            /*
             * Weighted Jacobi.
             */
            PressureSolverConfig
                jacobi_config;

            jacobi_config.kind =
                PressureSolverKind::
                    weighted_jacobi;

            jacobi_config.max_iterations =
                200000;

            jacobi_config.tolerance =
                1.0e-5;

            jacobi_config.relaxation =
                2.0 / 3.0;

            JacobiPressureSolver2D
                jacobi_solver;

            const Measurement jacobi =
                measure_solver(
                    jacobi_solver,
                    problem,
                    jacobi_config
                );

            print_measurement(
                n,
                "weighted_jacobi",
                jacobi
            );

            /*
             * Spectral solver.
             *
             * The first solve includes FFTW
             * planning time.
             */
            PressureSolverConfig
                spectral_config;

            spectral_config.kind =
                PressureSolverKind::
                    spectral_dct;

            spectral_config.tolerance =
                1.0e-10;

            SpectralPressureSolver2D
                spectral_solver;

            const Measurement
                spectral_first =
                    measure_solver(
                        spectral_solver,
                        problem,
                        spectral_config
                    );

            print_measurement(
                n,
                "spectral_dct_first",
                spectral_first
            );

            /*
             * Same solver object:
             *
             * FFTW plans are now cached,
             * therefore this approximates
             * the normal CFD timestep cost.
             */
            const Measurement
                spectral_cached =
                    measure_solver(
                        spectral_solver,
                        problem,
                        spectral_config
                    );

            print_measurement(
                n,
                "spectral_dct_cached",
                spectral_cached
            );

            if (
                jacobi.converged &&
                spectral_cached.converged &&
                spectral_cached.milliseconds >
                    0.0
            ) {

                const double speedup =
                    jacobi.milliseconds /
                    spectral_cached.milliseconds;

                std::cout
                    << "  speedup cached DCT vs Jacobi: "
                    << std::fixed
                    << std::setprecision(2)
                    << speedup
                    << "x\n";
            }

            std::cout
                << '\n';
        }

        return 0;

    } catch (
        const std::exception& exception
    ) {

        std::cerr
            << "Benchmark failed: "
            << exception.what()
            << '\n';

        return 1;
    }
}