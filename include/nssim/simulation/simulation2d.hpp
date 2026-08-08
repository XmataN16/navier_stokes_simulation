#pragma once

#include "nssim/problem/problem2d.hpp"
#include "nssim/simulation/simulation_config2d.hpp"
#include "nssim/solver/flow_solver2d.hpp"

#include <cstddef>
#include <functional>

namespace nssim {

struct SimulationProgress2D final {
    std::size_t step{};
    Real time{};

    StepDiagnostics diagnostics{};
};

struct SimulationSummary2D final {
    std::size_t steps{};

    Real final_time{};
    Real maximum_observed_cfl{};

    std::size_t nonconverged_pressure_steps{};
};

class Simulation2D final {
public:
    using FrameCallback =
        std::function<
            void(const FlowState2D&)
        >;

    using ProgressCallback =
        std::function<
            void(const SimulationProgress2D&)
        >;

    explicit Simulation2D(
        SimulationConfig2D config = {}
    );

    [[nodiscard]]
    SimulationSummary2D run(
        const Problem2D& problem,
        IFlowSolver2D& solver,
        const FrameCallback& frame_callback,
        const ProgressCallback& progress_callback = {}
    ) const;

private:
    SimulationConfig2D config_{};
};

} // namespace nssim