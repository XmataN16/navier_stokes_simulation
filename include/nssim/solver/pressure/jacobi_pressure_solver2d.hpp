#pragma once

#include "nssim/solver/pressure/pressure_solver2d.hpp"

namespace nssim {

class JacobiPressureSolver2D final
    : public IPressureSolver2D {

public:
    [[nodiscard]]
    PressureSolveResult solve(
        const UniformGrid2D& grid,
        const BoundarySet2D& boundaries,
        const Field2D& rhs,
        Field2D& pressure,
        const PressureSolverConfig& config
    ) override;

private:
    Field2D next_{};
};

} // namespace nssim