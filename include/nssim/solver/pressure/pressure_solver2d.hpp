#pragma once

#include "nssim/core/field2d.hpp"
#include "nssim/core/grid2d.hpp"

#include "nssim/problem/boundary_conditions2d.hpp"

#include "nssim/solver/pressure/pressure_solver_config.hpp"

#include <cstddef>

namespace nssim {

struct PressureSolveResult final {
    std::size_t iterations{};

    Real residual_l2{};

    bool converged{};
};

class IPressureSolver2D {
public:
    virtual ~IPressureSolver2D() = default;

    [[nodiscard]]
    virtual PressureSolveResult solve(
        const UniformGrid2D& grid,
        const BoundarySet2D& boundaries,
        const Field2D& rhs,
        Field2D& pressure,
        const PressureSolverConfig& config
    ) = 0;
};

} // namespace nssim