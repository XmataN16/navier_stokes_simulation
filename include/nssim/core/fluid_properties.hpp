#pragma once

#include "nssim/core/types.hpp"

#include <stdexcept>

namespace nssim {

struct FluidProperties final {
    Real density{};
    Real dynamic_viscosity{};

    [[nodiscard]]
    Real kinematic_viscosity() const noexcept {
        return dynamic_viscosity / density;
    }

    void validate() const {
        if (density <= Real{}) {
            throw std::invalid_argument{
                "Fluid density must be positive"
            };
        }

        if (dynamic_viscosity <= Real{}) {
            throw std::invalid_argument{
                "Dynamic viscosity must be positive"
            };
        }
    }
};

} // namespace nssim