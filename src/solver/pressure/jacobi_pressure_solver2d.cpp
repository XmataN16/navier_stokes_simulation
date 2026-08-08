#include "nssim/solver/pressure/jacobi_pressure_solver2d.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace nssim {

namespace {

[[nodiscard]]
bool is_dirichlet_pressure(
    const BoundaryCondition2D& boundary
) noexcept {

    return
        boundary.kind ==
        BoundaryKind::pressure_outlet;
}

[[nodiscard]]
bool has_dirichlet_pressure(
    const BoundarySet2D& boundaries
) noexcept {

    return
        is_dirichlet_pressure(boundaries.left) ||
        is_dirichlet_pressure(boundaries.right) ||
        is_dirichlet_pressure(boundaries.bottom) ||
        is_dirichlet_pressure(boundaries.top);
}

void validate_shapes(
    const UniformGrid2D& grid,
    const Field2D& rhs,
    const Field2D& pressure
) {

    if (
        rhs.nx() != grid.nx_cells() ||
        rhs.ny() != grid.ny_cells() ||

        pressure.nx() != grid.nx_cells() ||
        pressure.ny() != grid.ny_cells()
    ) {
        throw std::invalid_argument{
            "Pressure solver field dimensions "
            "do not match grid"
        };
    }
}

struct EquationTerms final {
    Real sum{};
    Real diagonal{};
};

[[nodiscard]]
EquationTerms equation_terms(
    const UniformGrid2D& grid,
    const BoundarySet2D& boundaries,
    const Field2D& pressure,
    const std::size_t i,
    const std::size_t j
) {

    const Real inv_dx2 =
        Real{1} /
        (grid.dx() * grid.dx());

    const Real inv_dy2 =
        Real{1} /
        (grid.dy() * grid.dy());

    const auto nx =
        grid.nx_cells();

    const auto ny =
        grid.ny_cells();

    EquationTerms result{};

    /*
     * West
     */
    if (i > 0) {

        result.sum +=
            pressure(i - 1, j) *
            inv_dx2;

        result.diagonal +=
            inv_dx2;

    } else if (
        is_dirichlet_pressure(
            boundaries.left
        )
    ) {

        result.sum +=
            Real{2} *
            boundaries.left.pressure *
            inv_dx2;

        result.diagonal +=
            Real{2} *
            inv_dx2;
    }

    /*
     * East
     */
    if (i + 1 < nx) {

        result.sum +=
            pressure(i + 1, j) *
            inv_dx2;

        result.diagonal +=
            inv_dx2;

    } else if (
        is_dirichlet_pressure(
            boundaries.right
        )
    ) {

        result.sum +=
            Real{2} *
            boundaries.right.pressure *
            inv_dx2;

        result.diagonal +=
            Real{2} *
            inv_dx2;
    }

    /*
     * South
     */
    if (j > 0) {

        result.sum +=
            pressure(i, j - 1) *
            inv_dy2;

        result.diagonal +=
            inv_dy2;

    } else if (
        is_dirichlet_pressure(
            boundaries.bottom
        )
    ) {

        result.sum +=
            Real{2} *
            boundaries.bottom.pressure *
            inv_dy2;

        result.diagonal +=
            Real{2} *
            inv_dy2;
    }

    /*
     * North
     */
    if (j + 1 < ny) {

        result.sum +=
            pressure(i, j + 1) *
            inv_dy2;

        result.diagonal +=
            inv_dy2;

    } else if (
        is_dirichlet_pressure(
            boundaries.top
        )
    ) {

        result.sum +=
            Real{2} *
            boundaries.top.pressure *
            inv_dy2;

        result.diagonal +=
            Real{2} *
            inv_dy2;
    }

    return result;
}

void remove_mean(
    Field2D& pressure
) {

    Real sum{};

    for (
        const Real value :
        pressure.values()
    ) {
        sum += value;
    }

    const Real mean =
        sum /
        static_cast<Real>(
            pressure.size()
        );

    for (
        Real& value :
        pressure.values()
    ) {
        value -= mean;
    }
}

[[nodiscard]]
Real compute_residual_l2(
    const UniformGrid2D& grid,
    const BoundarySet2D& boundaries,
    const Field2D& rhs,
    const Field2D& pressure
) {

    Real squared_sum{};

    for (
        std::size_t j = 0;
        j < pressure.ny();
        ++j
    ) {
        for (
            std::size_t i = 0;
            i < pressure.nx();
            ++i
        ) {

            const auto terms =
                equation_terms(
                    grid,
                    boundaries,
                    pressure,
                    i,
                    j
                );

            const Real residual =
                terms.sum -
                terms.diagonal *
                    pressure(i, j) -
                rhs(i, j);

            squared_sum +=
                residual * residual;
        }
    }

    return std::sqrt(
        squared_sum /
        static_cast<Real>(
            pressure.size()
        )
    );
}

} // namespace

PressureSolveResult
JacobiPressureSolver2D::solve(
    const UniformGrid2D& grid,
    const BoundarySet2D& boundaries,
    const Field2D& rhs,
    Field2D& pressure,
    const PressureSolverConfig& config
) {

    config.validate();

    validate_shapes(
        grid,
        rhs,
        pressure
    );

    if (
        config.kind !=
            PressureSolverKind::jacobi &&
        config.kind !=
            PressureSolverKind::weighted_jacobi
    ) {
        throw std::invalid_argument{
            "JacobiPressureSolver2D received "
            "unsupported solver kind"
        };
    }

    if (
        next_.nx() != pressure.nx() ||
        next_.ny() != pressure.ny()
    ) {
        next_ =
            Field2D{
                pressure.nx(),
                pressure.ny()
            };
    }

    const Real omega =
        config.kind ==
            PressureSolverKind::jacobi
        ? Real{1}
        : config.relaxation;

    /*
     * With only Neumann conditions pressure
     * is defined up to an arbitrary constant.
     *
     * We choose a zero-mean pressure gauge.
     */
    const bool remove_pressure_nullspace =
        !has_dirichlet_pressure(
            boundaries
        );

    PressureSolveResult result{};

    for (
        std::size_t iteration = 1;
        iteration <= config.max_iterations;
        ++iteration
    ) {

        for (
            std::size_t j = 0;
            j < pressure.ny();
            ++j
        ) {
            for (
                std::size_t i = 0;
                i < pressure.nx();
                ++i
            ) {

                const auto terms =
                    equation_terms(
                        grid,
                        boundaries,
                        pressure,
                        i,
                        j
                    );

                const Real jacobi_value =
                    (
                        terms.sum -
                        rhs(i, j)
                    ) /
                    terms.diagonal;

                next_(i, j) =
                    (
                        Real{1} -
                        omega
                    ) *
                    pressure(i, j)
                    +
                    omega *
                    jacobi_value;
            }
        }

        if (
            remove_pressure_nullspace
        ) {
            remove_mean(next_);
        }

        std::swap(
            pressure,
            next_
        );

        const Real residual =
            compute_residual_l2(
                grid,
                boundaries,
                rhs,
                pressure
            );

        result.iterations =
            iteration;

        result.residual_l2 =
            residual;

        result.converged =
            residual <=
            config.tolerance;

        if (result.converged) {
            return result;
        }
    }

    return result;
}

} // namespace nssim