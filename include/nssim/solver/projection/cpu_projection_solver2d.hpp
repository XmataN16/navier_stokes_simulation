#pragma once

#include "nssim/solver/flow_solver2d.hpp"
#include "nssim/solver/pressure/pressure_solver2d.hpp"
#include "nssim/solver/projection_solver_config.hpp"

#include <memory>
#include <optional>

namespace nssim {

class CpuProjectionSolver2D final
    : public IFlowSolver2D {

public:
    explicit CpuProjectionSolver2D(
        ProjectionSolverConfig config = {}
    );

    [[nodiscard]]
    std::string_view
    name() const noexcept override;

    void initialize(
        const Problem2D& problem
    ) override;

    [[nodiscard]]
    StepDiagnostics advance(
        Real dt
    ) override;

    void copy_state_to(
        FlowState2D& destination
    ) const override;

    [[nodiscard]]
    Real time() const noexcept override;

    [[nodiscard]]
    std::size_t
    step() const noexcept override;

private:
    ProjectionSolverConfig config_{};

    const Problem2D* problem_{};

    std::optional<FlowState2D>
        state_{};

    std::optional<MacVelocityField2D>
        predicted_velocity_{};

    Field2D pressure_rhs_{};

    std::unique_ptr<IPressureSolver2D>
        pressure_solver_{};
};

} // namespace nssim