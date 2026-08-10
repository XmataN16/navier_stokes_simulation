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
         * 32x32 is kept as the first
         * validation grid.
         *
         * Later it can be increased to:
         *
         * 64x64
         * 128x128
         * 256x256
         *
         * for grid-convergence studies.
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
         *
         * Pressure uses homogeneous
         * Neumann boundary conditions
         * on all four walls.
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
         * Fluid initially at rest:
         *
         * u = 0
         * v = 0
         * p = 0
         */
        const Problem2D problem{
            grid,
            fluid,
            boundaries,
            InitialCondition2D::rest()
        };

        /*
         * Projection method configuration.
         *
         * Use a direct spectral Poisson
         * solver based on the DCT.
         *
         * For the closed rectangular cavity
         * with homogeneous Neumann pressure
         * boundary conditions, the discrete
         * pressure Laplacian can be diagonalized
         * using a cosine transform.
         */
        ProjectionSolverConfig
            solver_config;

        solver_config.pressure.kind =
            PressureSolverKind::
                spectral_dct;

        /*
         * The spectral solver is direct,
         * therefore:
         *
         * - max_iterations is not used;
         * - relaxation is not used.
         *
         * tolerance is used to verify the
         * residual after the direct solve.
         */
        solver_config.pressure.tolerance =
            1.0e-10;

        /*
         * Sequential CPU implementation of
         * the projection method.
         *
         * The pressure algorithm is selected
         * through solver_config.
         */
        CpuProjectionSolver2D solver{
            solver_config
        };

        /*
         * Simulation parameters.
         *
         * dt = 0.01
         * T  = 5
         *
         * This gives:
         *
         * 500 physical time steps.
         *
         * Output every 10 steps:
         *
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
         * Output structure:
         *
         * output/cavity_re100/
         *
         *   cavity_re100.pvd
         *
         *   cavity_re100_000000.vti
         *   cavity_re100_000010.vti
         *   cavity_re100_000020.vti
         *   ...
         *   cavity_re100_000500.vti
         *
         * Open cavity_re100.pvd in ParaView
         * to load the complete time series.
         */
        VtkTimeSeriesWriter2D writer{
            "output/cavity_re100",
            "cavity_re100"
        };

        /*
         * Remove frames left from a
         * previous execution of this series.
         */
        writer.reset();

        Simulation2D simulation{
            simulation_config
        };

        /*
         * Print simulation configuration.
         */
        std::cout
            << std::fixed
            << std::setprecision(6)

            << "Lid-driven cavity simulation\n"

            << "Flow solver: "
            << solver.name()
            << '\n'

            << "Pressure solver: spectral DCT\n"

            << "Grid: "
            << grid.nx_cells()
            << 'x'
            << grid.ny_cells()
            << '\n'

            << "L = "
            << cavity_length
            << '\n'

            << "U_lid = "
            << lid_velocity
            << '\n'

            << "Re = "
            << reynolds_number
            << '\n'

            << "rho = "
            << fluid.density
            << '\n'

            << "nu = "
            << fluid.
                kinematic_viscosity()
            << '\n'

            << "mu = "
            << fluid.dynamic_viscosity
            << '\n'

            << "dx = "
            << grid.dx()
            << '\n'

            << "dy = "
            << grid.dy()
            << '\n'

            << "dt = "
            << simulation_config.
                time_step
            << '\n'

            << "T = "
            << simulation_config.
                end_time

            << "\n\n";

        /*
         * Measure complete wall-clock
         * simulation time.
         */
        const auto wall_start =
            std::chrono::
                steady_clock::now();

        const SimulationSummary2D summary =
            simulation.run(

                problem,
                solver,

                /*
                 * Frame output callback.
                 *
                 * Simulation2D calls this for:
                 *
                 * step 0,
                 * every output_interval_steps,
                 * and the final state.
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

        /*
         * Final simulation summary.
         */
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