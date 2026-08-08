#include "nssim/core/grid2d.hpp"

#include "nssim/problem/boundary_conditions2d.hpp"
#include "nssim/problem/initial_condition2d.hpp"
#include "nssim/problem/problem2d.hpp"

#include <cmath>
#include <cstdlib>

int main() {

    using namespace nssim;

    const UniformGrid2D grid{
        8,
        4,
        2.0,
        1.0
    };

    InitialCondition2D initial;

    initial.velocity =
        [](const Vec2 point) {

            return Vec2{
                .x = point.x,
                .y = -point.y
            };
        };

    initial.pressure =
        [](const Vec2 point) {

            return
                point.x +
                point.y;
        };

    const Problem2D problem{

        grid,

        FluidProperties{
            .density =
                1000.0,

            .dynamic_viscosity =
                1.0e-3
        },

        BoundarySet2D{},

        initial
    };

    const auto state =
        problem.make_initial_state();

    /*
     * Pressure:
     * Nx x Ny
     */
    if (
        state.pressure().nx() != 8 ||
        state.pressure().ny() != 4
    ) {
        return EXIT_FAILURE;
    }

    /*
     * U:
     * (Nx + 1) x Ny
     */
    if (
        state.velocity().u().nx() != 9 ||
        state.velocity().u().ny() != 4
    ) {
        return EXIT_FAILURE;
    }

    /*
     * V:
     * Nx x (Ny + 1)
     */
    if (
        state.velocity().v().nx() != 8 ||
        state.velocity().v().ny() != 5
    ) {
        return EXIT_FAILURE;
    }

    const Vec2 first_cell =
        grid.cell_center(
            0,
            0
        );

    const Real expected_pressure =
        first_cell.x +
        first_cell.y;

    if (
        std::abs(
            state.pressure()(0, 0) -
            expected_pressure
        ) > 1.0e-12
    ) {
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}