#include "nssim/simulation/simulation2d.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace nssim {

Simulation2D::Simulation2D(
    SimulationConfig2D config
)
    : config_{
        std::move(config)
      } {

    config_.validate();
}

SimulationSummary2D
Simulation2D::run(
    const Problem2D& problem,
    IFlowSolver2D& solver,
    const FrameCallback& frame_callback,
    const ProgressCallback& progress_callback
) const {

    config_.validate();

    if (!frame_callback) {
        throw std::invalid_argument{
            "Simulation2D requires "
            "a frame callback"
        };
    }

    /*
     * Initialize the selected solver.
     */
    solver.initialize(
        problem
    );

    FlowState2D output_state{
        problem.grid()
    };

    /*
     * Initial state:
     *
     * step = 0
     * t = 0
     */
    solver.copy_state_to(
        output_state
    );

    frame_callback(
        output_state
    );

    std::size_t last_output_step =
        solver.step();

    SimulationSummary2D summary{};

    /*
     * Small tolerance to avoid an additional
     * floating point time step near end_time.
     */
    const Real time_epsilon =
        Real{1.0e-12} *
        std::max(
            Real{1},
            std::abs(
                config_.end_time
            )
        );

    while (
        solver.time() <
        config_.end_time -
        time_epsilon
    ) {

        const Real remaining_time =
            config_.end_time -
            solver.time();

        /*
         * Shorten the final time step if
         * end_time is not exactly divisible by dt.
         */
        const Real dt =
            std::min(
                config_.time_step,
                remaining_time
            );

        const StepDiagnostics diagnostics =
            solver.advance(
                dt
            );

        ++summary.steps;

        summary.final_time =
            solver.time();

        summary.maximum_observed_cfl =
            std::max(
                summary.maximum_observed_cfl,
                diagnostics.cfl
            );

        /*
         * Pressure convergence control.
         */
        if (
            !diagnostics.pressure_converged
        ) {

            ++summary.
                nonconverged_pressure_steps;

            if (
                config_.
                    fail_on_pressure_nonconvergence
            ) {

                std::ostringstream message;

                message
                    << "Pressure solver failed "
                    << "to converge at step "
                    << solver.step()
                    << ", t = "
                    << solver.time()
                    << ", residual L2 = "
                    << diagnostics.
                        pressure_residual_l2;

                throw std::runtime_error{
                    message.str()
                };
            }
        }

        /*
         * Explicit time integration CFL check.
         */
        if (
            !std::isfinite(
                diagnostics.cfl
            ) ||
            diagnostics.cfl >
                config_.maximum_cfl
        ) {

            std::ostringstream message;

            message
                << "CFL limit exceeded at step "
                << solver.step()
                << ", t = "
                << solver.time()
                << ": CFL = "
                << diagnostics.cfl
                << ", limit = "
                << config_.maximum_cfl;

            throw std::runtime_error{
                message.str()
            };
        }

        const bool is_final_step =
            solver.time() >=
            config_.end_time -
            time_epsilon;

        /*
         * Diagnostics callback.
         */
        if (
            progress_callback &&
            (
                solver.step() %
                    config_.
                        diagnostics_interval_steps ==
                    0 ||
                is_final_step
            )
        ) {

            progress_callback(
                SimulationProgress2D{
                    .step =
                        solver.step(),

                    .time =
                        solver.time(),

                    .diagnostics =
                        diagnostics
                }
            );
        }

        /*
         * VTK/output callback.
         */
        if (
            solver.step() %
                config_.
                    output_interval_steps ==
                0 ||
            is_final_step
        ) {

            solver.copy_state_to(
                output_state
            );

            frame_callback(
                output_state
            );

            last_output_step =
                solver.step();
        }
    }

    /*
     * Guarantee that the final solution
     * is available in output.
     */
    if (
        last_output_step !=
        solver.step()
    ) {

        solver.copy_state_to(
            output_state
        );

        frame_callback(
            output_state
        );
    }

    summary.steps =
        solver.step();

    summary.final_time =
        solver.time();

    return summary;
}

} // namespace nssim