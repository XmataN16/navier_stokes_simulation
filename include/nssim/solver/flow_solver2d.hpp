#pragma once

#include "nssim/core/state2d.hpp"
#include "nssim/core/types.hpp"

#include "nssim/problem/problem2d.hpp"

#include <cstddef>
#include <string_view>

namespace nssim {

struct StepDiagnostics final {
    std::size_t pressure_iterations{};

    Real pressure_residual_l2{};
    Real divergence_l2{};
    Real cfl{};

    bool pressure_converged{};
};

class IFlowSolver2D {
public:
    virtual ~IFlowSolver2D() = default;

    [[nodiscard]]
    virtual std::string_view
    name() const noexcept = 0;

    virtual void initialize(
        const Problem2D& problem
    ) = 0;

    [[nodiscard]]
    virtual StepDiagnostics advance(
        Real dt
    ) = 0;

    virtual void copy_state_to(
        FlowState2D& destination
    ) const = 0;

    [[nodiscard]]
    virtual Real
    time() const noexcept = 0;

    [[nodiscard]]
    virtual std::size_t
    step() const noexcept = 0;
};

} // namespace nssim