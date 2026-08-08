#include "nssim/problem/problem2d.hpp"

#include "nssim/solver/projection/cpu_projection_solver2d.hpp"

#include <cstdlib>

int main() {

    using namespace nssim;

    const UniformGrid2D grid{
        16,
        16,
        1.0,
        1.0
    };

    const BoundarySet2D boundaries{

        .left =
            BoundaryCondition2D::no_slip(),

        .right =
            BoundaryCondition2D::no_slip(),

        .bottom =
            BoundaryCondition2D::no_slip(),

        .top =
            BoundaryCondition2D::
                moving_wall(
                    {1.0, 0.0}
                )
    };

    const Problem2D problem{

        grid,

        FluidProperties{
            .density = 1.0,
            .dynamic_viscosity = 0.01
        },

        boundaries,

        InitialCondition2D::rest()
    };

    ProjectionSolverConfig config;

    config.pressure.kind =
        PressureSolverKind::
            weighted_jacobi;

    config.pressure.max_iterations =
        5000;

    config.pressure.tolerance =
        1.0e-5;

    config.pressure.relaxation =
        2.0 / 3.0;

    CpuProjectionSolver2D solver{
        config
    };

    solver.initialize(
        problem
    );

    const StepDiagnostics diagnostics =
        solver.advance(
            0.001
        );

    if (
        !diagnostics.pressure_converged
    ) {
        return EXIT_FAILURE;
    }

    if (
        diagnostics.divergence_l2 >
        1.0e-7
    ) {
        return EXIT_FAILURE;
    }

    if (solver.step() != 1) {
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}