#include "nssim/core/fluid_properties.hpp"
#include "nssim/core/grid2d.hpp"

#include "nssim/io/vti_writer2d.hpp"

#include "nssim/problem/boundary_conditions2d.hpp"
#include "nssim/problem/initial_condition2d.hpp"
#include "nssim/problem/problem2d.hpp"

#include <exception>
#include <filesystem>
#include <iostream>

int main() {

    try {

        using namespace nssim;

        /*
         * Computational domain
         */
        const UniformGrid2D grid{
            128,
            64,
            2.0,
            1.0
        };

        /*
         * Water-like fluid
         */
        const FluidProperties water{
            .density =
                998.2,

            .dynamic_viscosity =
                1.002e-3
        };

        /*
         * Boundary conditions
         */
        const BoundarySet2D boundaries{

            .left =
                BoundaryCondition2D::
                    velocity_inlet(
                        {1.0, 0.0}
                    ),

            .right =
                BoundaryCondition2D::
                    pressure_outlet(
                        0.0
                    ),

            .bottom =
                BoundaryCondition2D::
                    no_slip(),

            .top =
                BoundaryCondition2D::
                    no_slip()
        };

        /*
         * Initial conditions
         */
        InitialCondition2D initial =
            InitialCondition2D::rest();

        const Real length_x =
            grid.length_x();

        /*
         * Temporary non-zero pressure
         * distribution just to test
         * ParaView output.
         */
        initial.pressure =
            [length_x](
                const Vec2 position
            ) {

                return
                    Real{1} -
                    position.x /
                    length_x;
            };

        /*
         * Complete physical problem
         */
        const Problem2D problem{
            grid,
            water,
            boundaries,
            initial
        };

        const FlowState2D state =
            problem.make_initial_state();

        /*
         * Output
         */
        const std::filesystem::path
            output_path{
                "output/initial_state.vti"
            };

        const VtiWriter2D writer;

        writer.write(
            problem,
            state,
            output_path
        );

        std::cout
            << "Project bootstrap is working.\n"

            << "Grid: "
            << grid.nx_cells()
            << 'x'
            << grid.ny_cells()
            << '\n'

            << "dx = "
            << grid.dx()
            << ", dy = "
            << grid.dy()
            << '\n'

            << "rho = "
            << water.density
            << ", nu = "
            << water.kinematic_viscosity()
            << '\n'

            << "Wrote: "
            << output_path.string()
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