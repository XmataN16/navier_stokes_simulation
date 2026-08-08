#include "nssim/problem/problem2d.hpp"

#include "nssim/simulation/simulation2d.hpp"

#include "nssim/solver/projection/cpu_projection_solver2d.hpp"

#include <cmath>
#include <cstdlib>

int main() {

    using namespace nssim;

    const UniformGrid2D grid{
        8,
        8,
        1.0,
        1.0
    };

    const BoundarySet2D boundaries{

        .left =
            BoundaryCondition2D::
                no_slip(),

        .right =
            BoundaryCondition2D::
                no_slip(),

        .bottom =
            BoundaryCondition2D::
                no_slip(),

        .top =
            BoundaryCondition2D::
                moving_wall(
                    {1.0, 0.0}
                )
    };

    const Problem2D problem{

        grid,

        FluidProperties{
            .density =
                1.0,

            .dynamic_viscosity =
                0.01
        },

        boundaries,

        InitialCondition2D::rest()
    };

    ProjectionSolverConfig
        solver_config;

    solver_config.pressure.
        max_iterations =
            5000;

    solver_config.pressure.
        tolerance =
            1.0e-5;

    CpuProjectionSolver2D solver{
        solver_config
    };

    const Simulation2D simulation{

        SimulationConfig2D{

            .end_time =
                0.01,

            .time_step =
                0.001,

            .output_interval_steps =
                5,

            .diagnostics_interval_steps =
                5,

            .maximum_cfl =
                0.8,

            .fail_on_pressure_nonconvergence =
                true
        }
    };

    std::size_t frame_count{};

    const auto summary =
        simulation.run(

            problem,
            solver,

            [&](const FlowState2D&) {

                ++frame_count;
            }
        );

    /*
     * 0.01 / 0.001 = 10 steps.
     */
    if (
        summary.steps != 10
    ) {
        return EXIT_FAILURE;
    }

    if (
        std::abs(
            summary.final_time -
            Real{0.01}
        ) >
        Real{1.0e-12}
    ) {

        return EXIT_FAILURE;
    }

    /*
     * Frames:
     *
     * step 0
     * step 5
     * step 10
     */
    if (
        frame_count != 3
    ) {
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}