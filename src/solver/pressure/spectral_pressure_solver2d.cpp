#include "nssim/solver/pressure/spectral_pressure_solver2d.hpp"

#include <fftw3.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <numbers>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace nssim {

namespace {

static_assert(
    std::is_same_v<Real, double>,
    "SpectralPressureSolver2D currently requires Real = double"
);

[[nodiscard]]
bool is_supported_neumann_boundary(
    const BoundaryCondition2D& boundary
) noexcept {

    switch (boundary.kind) {

    case BoundaryKind::no_slip_wall:
    case BoundaryKind::moving_wall:
    case BoundaryKind::free_slip_wall:
    case BoundaryKind::velocity_inlet:
        return true;

    case BoundaryKind::pressure_outlet:
    case BoundaryKind::periodic:
        return false;
    }

    return false;
}

void validate_boundaries(
    const BoundarySet2D& boundaries
) {

    if (
        !is_supported_neumann_boundary(
            boundaries.left
        ) ||
        !is_supported_neumann_boundary(
            boundaries.right
        ) ||
        !is_supported_neumann_boundary(
            boundaries.bottom
        ) ||
        !is_supported_neumann_boundary(
            boundaries.top
        )
    ) {
        throw std::invalid_argument{
            "SpectralPressureSolver2D currently "
            "supports only homogeneous Neumann "
            "pressure boundary conditions"
        };
    }
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
            "Spectral pressure solver field "
            "dimensions do not match grid"
        };
    }
}

[[nodiscard]]
Real field_mean(
    const Field2D& field
) {

    Real sum{};

    for (
        const Real value :
        field.values()
    ) {
        sum += value;
    }

    return
        sum /
        static_cast<Real>(
            field.size()
        );
}

void remove_mean(
    Field2D& field
) {

    const Real mean =
        field_mean(
            field
        );

    for (
        Real& value :
        field.values()
    ) {
        value -= mean;
    }
}

[[nodiscard]]
Real neumann_laplacian_at(
    const UniformGrid2D& grid,
    const Field2D& pressure,
    const std::size_t i,
    const std::size_t j
) {

    const auto nx =
        grid.nx_cells();

    const auto ny =
        grid.ny_cells();

    const Real center =
        pressure(i, j);

    const Real inv_dx2 =
        Real{1} /
        (
            grid.dx() *
            grid.dx()
        );

    const Real inv_dy2 =
        Real{1} /
        (
            grid.dy() *
            grid.dy()
        );

    Real result{};

    /*
     * At a Neumann boundary the missing
     * neighbour contributes zero because
     *
     * dp/dn = 0.
     *
     * This is exactly the same discrete
     * operator as in the current Jacobi
     * implementation.
     */

    if (i > 0) {
        result +=
            (
                pressure(i - 1, j) -
                center
            ) *
            inv_dx2;
    }

    if (i + 1 < nx) {
        result +=
            (
                pressure(i + 1, j) -
                center
            ) *
            inv_dx2;
    }

    if (j > 0) {
        result +=
            (
                pressure(i, j - 1) -
                center
            ) *
            inv_dy2;
    }

    if (j + 1 < ny) {
        result +=
            (
                pressure(i, j + 1) -
                center
            ) *
            inv_dy2;
    }

    return result;
}

[[nodiscard]]
Real compute_residual_l2(
    const UniformGrid2D& grid,
    const Field2D& rhs,
    const Real rhs_mean,
    const Field2D& pressure
) {

    Real squared_sum{};

    for (
        std::size_t j = 0;
        j < grid.ny_cells();
        ++j
    ) {
        for (
            std::size_t i = 0;
            i < grid.nx_cells();
            ++i
        ) {

            /*
             * The Neumann problem requires
             *
             * mean(rhs) = 0.
             *
             * The spectral solver projects the
             * supplied RHS onto that compatible
             * subspace.
             */
            const Real compatible_rhs =
                rhs(i, j) -
                rhs_mean;

            const Real residual =
                neumann_laplacian_at(
                    grid,
                    pressure,
                    i,
                    j
                ) -
                compatible_rhs;

            squared_sum +=
                residual *
                residual;
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

struct SpectralPressureSolver2D::Impl final {

    std::size_t nx{};
    std::size_t ny{};
    std::size_t element_count{};

    double* work{};

    fftw_plan forward_plan{};
    fftw_plan inverse_plan{};

    ~Impl() {
        reset();
    }

    void reset() noexcept {

        if (forward_plan != nullptr) {

            fftw_destroy_plan(
                forward_plan
            );

            forward_plan = nullptr;
        }

        if (inverse_plan != nullptr) {

            fftw_destroy_plan(
                inverse_plan
            );

            inverse_plan = nullptr;
        }

        if (work != nullptr) {

            fftw_free(
                work
            );

            work = nullptr;
        }

        nx = 0;
        ny = 0;
        element_count = 0;
    }

    void prepare(
        const std::size_t new_nx,
        const std::size_t new_ny
    ) {

        if (
            new_nx == nx &&
            new_ny == ny &&
            work != nullptr &&
            forward_plan != nullptr &&
            inverse_plan != nullptr
        ) {
            return;
        }

        reset();

        if (
            new_nx >
                static_cast<std::size_t>(
                    std::numeric_limits<int>::max()
                ) ||
            new_ny >
                static_cast<std::size_t>(
                    std::numeric_limits<int>::max()
                )
        ) {
            throw std::overflow_error{
                "Grid dimensions exceed FFTW "
                "integer plan limits"
            };
        }

        nx =
            new_nx;

        ny =
            new_ny;

        element_count =
            nx * ny;

        work =
            static_cast<double*>(
                fftw_malloc(
                    sizeof(double) *
                    element_count
                )
            );

        if (work == nullptr) {

            reset();

            throw std::bad_alloc{};
        }

        /*
         * Field2D memory layout:
         *
         * index = j * nx + i
         *
         * therefore memory is logically:
         *
         * [ny][nx]
         *
         * and x is the fastest-changing index.
         *
         * FFTW uses row-major order, so:
         *
         * n0 = ny
         * n1 = nx
         */

        forward_plan =
            fftw_plan_r2r_2d(
                static_cast<int>(ny),
                static_cast<int>(nx),

                work,
                work,

                FFTW_REDFT10,
                FFTW_REDFT10,

                FFTW_MEASURE
            );

        if (forward_plan == nullptr) {

            reset();

            throw std::runtime_error{
                "Failed to create FFTW "
                "forward DCT plan"
            };
        }

        inverse_plan =
            fftw_plan_r2r_2d(
                static_cast<int>(ny),
                static_cast<int>(nx),

                work,
                work,

                FFTW_REDFT01,
                FFTW_REDFT01,

                FFTW_MEASURE
            );

        if (inverse_plan == nullptr) {

            reset();

            throw std::runtime_error{
                "Failed to create FFTW "
                "inverse DCT plan"
            };
        }
    }
};

SpectralPressureSolver2D::
SpectralPressureSolver2D()
    : impl_{
        std::make_unique<Impl>()
      } {
}

SpectralPressureSolver2D::
~SpectralPressureSolver2D() = default;

SpectralPressureSolver2D::
SpectralPressureSolver2D(
    SpectralPressureSolver2D&&
) noexcept = default;

SpectralPressureSolver2D&
SpectralPressureSolver2D::
operator=(
    SpectralPressureSolver2D&&
) noexcept = default;

PressureSolveResult
SpectralPressureSolver2D::solve(
    const UniformGrid2D& grid,
    const BoundarySet2D& boundaries,
    const Field2D& rhs,
    Field2D& pressure,
    const PressureSolverConfig& config
) {

    config.validate();

    if (
        config.kind !=
        PressureSolverKind::spectral_dct
    ) {
        throw std::invalid_argument{
            "SpectralPressureSolver2D received "
            "unsupported solver kind"
        };
    }

    validate_shapes(
        grid,
        rhs,
        pressure
    );

    validate_boundaries(
        boundaries
    );

    const auto nx =
        grid.nx_cells();

    const auto ny =
        grid.ny_cells();

    impl_->prepare(
        nx,
        ny
    );

    /*
     * The Poisson equation with only Neumann
     * boundary conditions is solvable only if
     *
     * sum(rhs) = 0.
     *
     * In the projection method this should be
     * satisfied by mass conservation, but small
     * floating-point errors are unavoidable.
     */
    const Real rhs_mean =
        field_mean(
            rhs
        );

    for (
        std::size_t j = 0;
        j < ny;
        ++j
    ) {
        for (
            std::size_t i = 0;
            i < nx;
            ++i
        ) {

            const std::size_t index =
                j * nx + i;

            impl_->work[index] =
                rhs(i, j) -
                rhs_mean;
        }
    }

    /*
     * Forward two-dimensional DCT-II.
     */
    fftw_execute(
        impl_->forward_plan
    );

    const Real pi =
        std::numbers::pi_v<Real>;

    const Real inv_dx2 =
        Real{1} /
        (
            grid.dx() *
            grid.dx()
        );

    const Real inv_dy2 =
        Real{1} /
        (
            grid.dy() *
            grid.dy()
        );

    /*
     * For the cell-centered discrete Neumann
     * Laplacian:
     *
     * lambda(k,l) =
     *
     * -4/dx^2 sin^2(pi*k/(2*Nx))
     * -4/dy^2 sin^2(pi*l/(2*Ny))
     *
     * Therefore every spectral coefficient can
     * be solved independently.
     */
    for (
        std::size_t j = 0;
        j < ny;
        ++j
    ) {

        const Real angle_y =
            pi *
            static_cast<Real>(j) /
            (
                Real{2} *
                static_cast<Real>(ny)
            );

        const Real sin_y =
            std::sin(
                angle_y
            );

        const Real lambda_y =
            -Real{4} *
            inv_dy2 *
            sin_y *
            sin_y;

        for (
            std::size_t i = 0;
            i < nx;
            ++i
        ) {

            const std::size_t index =
                j * nx + i;

            /*
             * Constant mode.
             *
             * lambda(0,0) = 0 because pressure
             * is defined only up to a constant.
             *
             * Setting it to zero chooses:
             *
             * mean(p) = 0.
             */
            if (
                i == 0 &&
                j == 0
            ) {

                impl_->work[index] =
                    Real{};

                continue;
            }

            const Real angle_x =
                pi *
                static_cast<Real>(i) /
                (
                    Real{2} *
                    static_cast<Real>(nx)
                );

            const Real sin_x =
                std::sin(
                    angle_x
                );

            const Real lambda_x =
                -Real{4} *
                inv_dx2 *
                sin_x *
                sin_x;

            const Real lambda =
                lambda_x +
                lambda_y;

            impl_->work[index] /=
                lambda;
        }
    }

    /*
     * Inverse two-dimensional DCT-III.
     */
    fftw_execute(
        impl_->inverse_plan
    );

    /*
     * FFTW real-to-real transforms are
     * unnormalized.
     *
     * DCT-II followed by DCT-III gives 2*N
     * in one dimension, therefore:
     *
     * scale = 1 / (4 * Nx * Ny)
     *
     * in two dimensions.
     */
    const Real normalization =
        Real{1} /
        (
            Real{4} *
            static_cast<Real>(nx) *
            static_cast<Real>(ny)
        );

    for (
        std::size_t j = 0;
        j < ny;
        ++j
    ) {
        for (
            std::size_t i = 0;
            i < nx;
            ++i
        ) {

            const std::size_t index =
                j * nx + i;

            pressure(i, j) =
                impl_->work[index] *
                normalization;
        }
    }

    /*
     * Remove tiny floating-point drift of
     * the constant pressure mode.
     */
    remove_mean(
        pressure
    );

    const Real residual =
        compute_residual_l2(
            grid,
            rhs,
            rhs_mean,
            pressure
        );

    return {
        .iterations =
            1,

        .residual_l2 =
            residual,

        .converged =
            residual <=
            config.tolerance
    };
}

} // namespace nssim