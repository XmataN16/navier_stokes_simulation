#include "nssim/core/fluid_properties.hpp"
#include "nssim/core/grid2d.hpp"

#include "nssim/io/vtk_time_series_writer2d.hpp"

#include "nssim/problem/boundary_conditions2d.hpp"
#include "nssim/problem/initial_condition2d.hpp"
#include "nssim/problem/problem2d.hpp"

#include "nssim/simulation/simulation2d.hpp"

#include "nssim/solver/pressure/pressure_solver_config.hpp"
#include "nssim/solver/projection/cpu_projection_solver2d.hpp"
#include "nssim/solver/projection_solver_config.hpp"

#include <chrono>
#include <exception>
#include <iomanip>
#include <iostream>

int main() {

    try {

        using namespace nssim;

        /*
         * Lid-driven cavity benchmark.
         *
         * L = 1
         * U_lid = 1
         * Re = 100
         *
         * Therefore:
         *
         * nu = U * L / Re = 0.01
         */
        constexpr Real cavity_length =
            1.0;

        constexpr Real lid_velocity =
            1.0;

        constexpr Real reynolds_number =
            100.0;

        constexpr Real density =
            1.0;

        const Real kinematic_viscosity =
            lid_velocity *
            cavity_length /
            reynolds_number;

        const Real dynamic_viscosity =
            density *
            kinematic_viscosity;

        /*
         * Reference grid.
         *
         * 32x32 is intentionally modest:
         * weighted Jacobi is still a very
         * slow pressure solver.
         */
        const UniformGrid2D grid{
            32,
            32,
            cavity_length,
            cavity_length
        };

        const FluidProperties fluid{
            .density =
                density,

            .dynamic_viscosity =
                dynamic_viscosity
        };

        /*
         * Lid-driven cavity:
         *
         * left   -> no slip
         * right  -> no slip
         * bottom -> no slip
         * top    -> moving wall
         */
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
                        {
                            lid_velocity,
                            0.0
                        }
                    )
        };

        /*
         * Fluid initially at rest.
         */
        const Problem2D problem{
            grid,
            fluid,
            boundaries,
            InitialCondition2D::rest()
        };

        /*
         * Projection method configuration.
         */
        ProjectionSolverConfig
            solver_config;

        solver_config.pressure.kind =
            PressureSolverKind::
                weighted_jacobi;

        /*
         * Pure Neumann pressure problem in
         * a closed cavity converges slowly
         * with Jacobi.
         */
        solver_config.pressure.
            max_iterations =
                30000;

        solver_config.pressure.
            tolerance =
                1.0e-5;

        solver_config.pressure.
            relaxation =
                2.0 / 3.0;

        CpuProjectionSolver2D solver{
            solver_config
        };

        /*
         * Simulation parameters.
         *
         * dt = 0.01
         * T  = 5
         *
         * 500 physical time steps.
         *
         * Output every 10 steps:
         * Delta t_output = 0.1
         */
        const SimulationConfig2D
            simulation_config{

                .end_time =
                    5.0,

                .time_step =
                    0.01,

                .output_interval_steps =
                    10,

                .diagnostics_interval_steps =
                    10,

                .maximum_cfl =
                    0.8,

                .fail_on_pressure_nonconvergence =
                    true
            };

        /*
         * Generates:
         *
         * output/cavity_re100/
         *   cavity_re100.pvd
         *   cavity_re100_000000.vti
         *   cavity_re100_000010.vti
         *   ...
         *   cavity_re100_000500.vti
         */
        VtkTimeSeriesWriter2D writer{
            "output/cavity_re100",
            "cavity_re100"
        };

        /*
         * Remove old frames belonging
         * to this series.
         */
        writer.reset();

        Simulation2D simulation{
            simulation_config
        };

        std::cout
            << std::fixed
            << std::setprecision(6)

            << "Lid-driven cavity simulation\n"

            << "Solver: "
            << solver.name()
            << '\n'

            << "Grid: "
            << grid.nx_cells()
            << 'x'
            << grid.ny_cells()
            << '\n'

            << "Re = "
            << reynolds_number
            << '\n'

            << "nu = "
            << fluid.
                kinematic_viscosity()
            << '\n'

            << "dt = "
            << simulation_config.
                time_step
            << '\n'

            << "T = "
            << simulation_config.
                end_time

            << "\n\n";

        const auto wall_start =
            std::chrono::
                steady_clock::now();

        const SimulationSummary2D summary =
            simulation.run(

                problem,
                solver,

                /*
                 * Frame output callback.
                 */
                [&](const FlowState2D& state) {

                    writer.write_frame(
                        problem,
                        state
                    );
                },

                /*
                 * Diagnostics callback.
                 */
                [](const SimulationProgress2D&
                       progress) {

                    std::cout

                        << "step = "
                        << std::setw(6)
                        << progress.step

                        << "  t = "
                        << std::setw(10)
                        << progress.time

                        << "  CFL = "
                        << std::setw(10)
                        << progress.
                            diagnostics.cfl

                        << "  div_L2 = "
                        << std::scientific
                        << progress.
                            diagnostics.
                                divergence_l2

                        << "  p_res_L2 = "
                        << progress.
                            diagnostics.
                                pressure_residual_l2

                        << std::fixed

                        << "  p_iter = "
                        << progress.
                            diagnostics.
                                pressure_iterations

                        << '\n';
                }
            );

        const auto wall_end =
            std::chrono::
                steady_clock::now();

        const std::chrono::duration<double>
            elapsed =
                wall_end -
                wall_start;

        std::cout
            << "\nSimulation completed.\n"

            << "Steps: "
            << summary.steps
            << '\n'

            << "Final time: "
            << summary.final_time
            << '\n'

            << "Maximum CFL: "
            << summary.
                maximum_observed_cfl
            << '\n'

            << "Pressure non-converged steps: "
            << summary.
                nonconverged_pressure_steps
            << '\n'

            << "Wall time: "
            << elapsed.count()
            << " s\n"

            << "ParaView collection: "
            << writer.
                collection_path().
                    string()
            << '\n';

        return 0;

    } catch (
        const std::exception& exception
    ) {

        std::cerr
            << "Fatal error: "
            << exception.what()
            << '\n';

        return 1;
    }
}