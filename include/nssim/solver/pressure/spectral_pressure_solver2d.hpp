#pragma once

#include "nssim/solver/pressure/pressure_solver2d.hpp"

#include <memory>

namespace nssim {

/*
 * Direct spectral Poisson solver for a rectangular,
 * uniform, cell-centered grid with homogeneous
 * Neumann pressure boundary conditions.
 *
 * The discrete Laplacian is diagonalized using
 * a two-dimensional DCT-II.
 *
 * Current limitations:
 *
 * - uniform Cartesian grid only;
 * - homogeneous Neumann pressure BC on all sides;
 * - no periodic boundaries;
 * - no pressure_outlet / Dirichlet pressure BC.
 */
class SpectralPressureSolver2D final
    : public IPressureSolver2D {

public:
    SpectralPressureSolver2D();

    ~SpectralPressureSolver2D() override;

    SpectralPressureSolver2D(
        const SpectralPressureSolver2D&
    ) = delete;

    SpectralPressureSolver2D&
    operator=(
        const SpectralPressureSolver2D&
    ) = delete;

    SpectralPressureSolver2D(
        SpectralPressureSolver2D&&
    ) noexcept;

    SpectralPressureSolver2D&
    operator=(
        SpectralPressureSolver2D&&
    ) noexcept;

    [[nodiscard]]
    PressureSolveResult solve(
        const UniformGrid2D& grid,
        const BoundarySet2D& boundaries,
        const Field2D& rhs,
        Field2D& pressure,
        const PressureSolverConfig& config
    ) override;

private:
    struct Impl;

    std::unique_ptr<Impl> impl_;
};

} // namespace nssim