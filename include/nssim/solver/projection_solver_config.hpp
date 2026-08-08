#pragma once

#include "nssim/solver/pressure/pressure_solver_config.hpp"

namespace nssim {

struct ProjectionSolverConfig final {
    PressureSolverConfig pressure{};

    void validate() const {
        pressure.validate();
    }
};

} // namespace nssim