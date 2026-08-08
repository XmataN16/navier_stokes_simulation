#pragma once

#include "nssim/core/field2d.hpp"
#include "nssim/core/grid2d.hpp"

#include <cstdint>

namespace nssim {

class MacVelocityField2D final {
public:
    explicit MacVelocityField2D(
        const UniformGrid2D& grid)
        : u_{
            grid.nx_cells() + 1,
            grid.ny_cells()
        },
          v_{
            grid.nx_cells(),
            grid.ny_cells() + 1
        } {
    }

    [[nodiscard]]
    Field2D& u() noexcept {
        return u_;
    }

    [[nodiscard]]
    const Field2D& u() const noexcept {
        return u_;
    }

    [[nodiscard]]
    Field2D& v() noexcept {
        return v_;
    }

    [[nodiscard]]
    const Field2D& v() const noexcept {
        return v_;
    }

private:
    Field2D u_;
    Field2D v_;
};

class FlowState2D final {
public:
    explicit FlowState2D(
        const UniformGrid2D& grid)
        : velocity_{grid},
          pressure_{
              grid.nx_cells(),
              grid.ny_cells()
          } {
    }

    [[nodiscard]]
    MacVelocityField2D& velocity() noexcept {
        return velocity_;
    }

    [[nodiscard]]
    const MacVelocityField2D&
    velocity() const noexcept {
        return velocity_;
    }

    [[nodiscard]]
    Field2D& pressure() noexcept {
        return pressure_;
    }

    [[nodiscard]]
    const Field2D&
    pressure() const noexcept {
        return pressure_;
    }

    [[nodiscard]]
    Real time() const noexcept {
        return time_;
    }

    [[nodiscard]]
    std::uint64_t step() const noexcept {
        return step_;
    }

    void set_clock(
        const Real time,
        const std::uint64_t step) noexcept {

        time_ = time;
        step_ = step;
    }

private:
    MacVelocityField2D velocity_;
    Field2D pressure_;

    Real time_{};
    std::uint64_t step_{};
};

} // namespace nssim