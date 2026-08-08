#pragma once

#include "nssim/core/fluid_properties.hpp"
#include "nssim/core/grid2d.hpp"
#include "nssim/core/state2d.hpp"

#include "nssim/problem/boundary_conditions2d.hpp"
#include "nssim/problem/initial_condition2d.hpp"

namespace nssim {

class Problem2D final {
public:
    Problem2D(
        UniformGrid2D grid,
        FluidProperties fluid,
        BoundarySet2D boundaries,
        InitialCondition2D initial_condition
    );

    [[nodiscard]]
    const UniformGrid2D&
    grid() const noexcept {
        return grid_;
    }

    [[nodiscard]]
    const FluidProperties&
    fluid() const noexcept {
        return fluid_;
    }

    [[nodiscard]]
    const BoundarySet2D&
    boundaries() const noexcept {
        return boundaries_;
    }

    [[nodiscard]]
    FlowState2D make_initial_state() const;

    void validate() const;

private:
    UniformGrid2D grid_;

    FluidProperties fluid_;

    BoundarySet2D boundaries_;

    InitialCondition2D initial_condition_;
};

} // namespace nssim