#pragma once

#include "nssim/core/types.hpp"

#include <functional>

namespace nssim {

struct InitialCondition2D final {

    using VelocityFunction =
        std::function<Vec2(Vec2)>;

    using PressureFunction =
        std::function<Real(Vec2)>;

    VelocityFunction velocity =
        [](const Vec2) {
            return Vec2{};
        };

    PressureFunction pressure =
        [](const Vec2) {
            return Real{};
        };

    [[nodiscard]]
    static InitialCondition2D rest() {
        return {};
    }
};

} // namespace nssim