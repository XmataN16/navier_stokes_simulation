#pragma once

#include "nssim/core/types.hpp"

namespace nssim {

enum class BoundaryKind {
    no_slip_wall,
    free_slip_wall,
    velocity_inlet,
    pressure_outlet,
    periodic
};

struct BoundaryCondition2D final {
    BoundaryKind kind{
        BoundaryKind::no_slip_wall
    };

    Vec2 velocity{};
    Real pressure{};

    [[nodiscard]]
    static constexpr
    BoundaryCondition2D no_slip() noexcept {
        return {};
    }

    [[nodiscard]]
    static constexpr
    BoundaryCondition2D free_slip() noexcept {

        return {
            .kind =
                BoundaryKind::free_slip_wall
        };
    }

    [[nodiscard]]
    static constexpr
    BoundaryCondition2D velocity_inlet(
        const Vec2 value) noexcept {

        return {
            .kind =
                BoundaryKind::velocity_inlet,

            .velocity = value
        };
    }

    [[nodiscard]]
    static constexpr
    BoundaryCondition2D pressure_outlet(
        const Real value) noexcept {

        return {
            .kind =
                BoundaryKind::pressure_outlet,

            .pressure = value
        };
    }

    [[nodiscard]]
    static constexpr
    BoundaryCondition2D periodic() noexcept {

        return {
            .kind =
                BoundaryKind::periodic
        };
    }
};

struct BoundarySet2D final {
    BoundaryCondition2D left{};
    BoundaryCondition2D right{};
    BoundaryCondition2D bottom{};
    BoundaryCondition2D top{};
};

} // namespace nssim