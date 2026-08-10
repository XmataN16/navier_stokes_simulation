#include "nssim/core/field2d.hpp"
#include "nssim/core/grid2d.hpp"

#include "nssim/problem/boundary_conditions2d.hpp"

#include "nssim/solver/pressure/spectral_pressure_solver2d.hpp"

#include <cmath>
#include <cstdlib>
#include <numbers>

int main() {

    using namespace nssim;

    /*
     * Intentionally rectangular grid:
     * verifies that dx != dy is handled.
     */
    const UniformGrid2D grid{
        32,
        24,
        2.0,
        1.0
    };

    const auto nx =
        grid.nx_cells();

    const auto ny =
        grid.ny_cells();

    Field2D exact_pressure{
        nx,
        ny
    };

    Field2D rhs{
        nx,
        ny
    };

    Field2D numerical_pressure{
        nx,
        ny
    };

    constexpr std::size_t mode_x =
        3;

    constexpr std::size_t mode_y =
        2;

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
                static_cast<Real>(mode_x) /
                (
                    Real{2} *
                    static_cast<Real>(nx)
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
                static_cast<Real>(mode_y) /
                (
                    Real{2} *
                    static_cast<Real>(ny)
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
                    static_cast<Real>(mode_x) *
                    (
                        static_cast<Real>(i) +
                        Real{0.5}
                    ) /
                    static_cast<Real>(nx)
                );

            const Real y_mode =
                std::cos(
                    pi *
                    static_cast<Real>(mode_y) *
                    (
                        static_cast<Real>(j) +
                        Real{0.5}
                    ) /
                    static_cast<Real>(ny)
                );

            exact_pressure(i, j) =
                x_mode *
                y_mode;

            /*
             * This is an exact eigenfunction
             * of our discrete Neumann Laplacian.
             */
            rhs(i, j) =
                lambda *
                exact_pressure(i, j);
        }
    }

    PressureSolverConfig config;

    config.kind =
        PressureSolverKind::
            spectral_dct;

    config.tolerance =
        1.0e-10;

    SpectralPressureSolver2D solver;

    const PressureSolveResult result =
        solver.solve(
            grid,
            BoundarySet2D{},
            rhs,
            numerical_pressure,
            config
        );

    if (!result.converged) {
        return EXIT_FAILURE;
    }

    if (
        result.residual_l2 >
        Real{1.0e-10}
    ) {
        return EXIT_FAILURE;
    }

    Real squared_error{};
    Real maximum_error{};

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

            const Real error =
                numerical_pressure(i, j) -
                exact_pressure(i, j);

            squared_error +=
                error *
                error;

            maximum_error =
                std::max(
                    maximum_error,
                    std::abs(error)
                );
        }
    }

    const Real error_l2 =
        std::sqrt(
            squared_error /
            static_cast<Real>(
                nx * ny
            )
        );

    if (
        error_l2 >
        Real{1.0e-11}
    ) {
        return EXIT_FAILURE;
    }

    if (
        maximum_error >
        Real{1.0e-10}
    ) {
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}