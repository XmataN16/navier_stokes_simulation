#include "nssim/solver/projection/cpu_projection_solver2d.hpp"

#include "nssim/solver/pressure/jacobi_pressure_solver2d.hpp"
#include "nssim/solver/pressure/spectral_pressure_solver2d.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace nssim {

namespace {

[[nodiscard]]
bool is_periodic(
    const BoundaryCondition2D& boundary
) noexcept {

    return
        boundary.kind ==
        BoundaryKind::periodic;
}

void reject_periodic_boundaries(
    const BoundarySet2D& boundaries
) {

    if (
        is_periodic(boundaries.left) ||
        is_periodic(boundaries.right) ||
        is_periodic(boundaries.bottom) ||
        is_periodic(boundaries.top)
    ) {
        throw std::invalid_argument{
            "CpuProjectionSolver2D does not "
            "support periodic boundaries yet"
        };
    }
}

[[nodiscard]]
Real normal_boundary_value(
    const BoundaryCondition2D& boundary,
    const Real prescribed_component,
    const Real adjacent_value
) {

    switch (boundary.kind) {

    case BoundaryKind::no_slip_wall:
    case BoundaryKind::free_slip_wall:
        return Real{};

    case BoundaryKind::moving_wall:
    case BoundaryKind::velocity_inlet:
        return prescribed_component;

    case BoundaryKind::pressure_outlet:
        return adjacent_value;

    case BoundaryKind::periodic:
        break;
    }

    throw std::logic_error{
        "Unsupported periodic boundary"
    };
}

[[nodiscard]]
Real tangential_ghost_value(
    const Real inside_value,
    const BoundaryCondition2D& boundary,
    const Real prescribed_component
) {

    switch (boundary.kind) {

    case BoundaryKind::no_slip_wall:
        return -inside_value;

    case BoundaryKind::moving_wall:
    case BoundaryKind::velocity_inlet:

        return
            Real{2} *
            prescribed_component -
            inside_value;

    case BoundaryKind::free_slip_wall:
    case BoundaryKind::pressure_outlet:

        return inside_value;

    case BoundaryKind::periodic:
        break;
    }

    throw std::logic_error{
        "Unsupported periodic boundary"
    };
}

void apply_velocity_boundaries(
    MacVelocityField2D& velocity,
    const UniformGrid2D& grid,
    const BoundarySet2D& boundaries
) {

    auto& u =
        velocity.u();

    auto& v =
        velocity.v();

    const auto nx =
        grid.nx_cells();

    const auto ny =
        grid.ny_cells();

    /*
     * Left / right:
     * u is normal to the boundary.
     */
    for (
        std::size_t j = 0;
        j < ny;
        ++j
    ) {

        u(0, j) =
            normal_boundary_value(
                boundaries.left,
                boundaries.left.velocity.x,
                u(1, j)
            );

        u(nx, j) =
            normal_boundary_value(
                boundaries.right,
                boundaries.right.velocity.x,
                u(nx - 1, j)
            );
    }

    /*
     * Bottom / top:
     * v is normal to the boundary.
     */
    for (
        std::size_t i = 0;
        i < nx;
        ++i
    ) {

        v(i, 0) =
            normal_boundary_value(
                boundaries.bottom,
                boundaries.bottom.velocity.y,
                v(i, 1)
            );

        v(i, ny) =
            normal_boundary_value(
                boundaries.top,
                boundaries.top.velocity.y,
                v(i, ny - 1)
            );
    }
}

void compute_predictor(
    const Problem2D& problem,
    const MacVelocityField2D& current,
    MacVelocityField2D& predicted,
    const Real dt
) {

    const auto& grid =
        problem.grid();

    const auto& boundaries =
        problem.boundaries();

    const Real dx =
        grid.dx();

    const Real dy =
        grid.dy();

    const Real inv_dx =
        Real{1} / dx;

    const Real inv_dy =
        Real{1} / dy;

    const Real inv_dx2 =
        inv_dx * inv_dx;

    const Real inv_dy2 =
        inv_dy * inv_dy;

    const Real nu =
        problem.fluid()
            .kinematic_viscosity();

    const auto nx =
        grid.nx_cells();

    const auto ny =
        grid.ny_cells();

    const auto& u =
        current.u();

    const auto& v =
        current.v();

    predicted =
        current;

    auto& u_star =
        predicted.u();

    auto& v_star =
        predicted.v();

    /*
     * Predictor for u.
     *
     * u is located at vertical MAC faces.
     */
    for (
        std::size_t j = 0;
        j < ny;
        ++j
    ) {
        for (
            std::size_t i = 1;
            i < nx;
            ++i
        ) {

            const Real u_center =
                u(i, j);

            const Real v_at_u =
                Real{0.25} *
                (
                    v(i - 1, j) +
                    v(i, j) +
                    v(i - 1, j + 1) +
                    v(i, j + 1)
                );

            const Real south =
                j > 0
                ? u(i, j - 1)
                : tangential_ghost_value(
                    u_center,
                    boundaries.bottom,
                    boundaries.bottom.velocity.x
                );

            const Real north =
                j + 1 < ny
                ? u(i, j + 1)
                : tangential_ghost_value(
                    u_center,
                    boundaries.top,
                    boundaries.top.velocity.x
                );

            const Real du_dx =
                u_center >= Real{}
                ?
                (
                    u_center -
                    u(i - 1, j)
                ) *
                inv_dx
                :
                (
                    u(i + 1, j) -
                    u_center
                ) *
                inv_dx;

            const Real du_dy =
                v_at_u >= Real{}
                ?
                (
                    u_center -
                    south
                ) *
                inv_dy
                :
                (
                    north -
                    u_center
                ) *
                inv_dy;

            const Real laplacian_u =
                (
                    u(i - 1, j) -
                    Real{2} * u_center +
                    u(i + 1, j)
                ) *
                inv_dx2
                +
                (
                    south -
                    Real{2} * u_center +
                    north
                ) *
                inv_dy2;

            const Real advection =
                u_center * du_dx +
                v_at_u * du_dy;

            u_star(i, j) =
                u_center +
                dt *
                (
                    -advection +
                    nu * laplacian_u
                );
        }
    }

    /*
     * Predictor for v.
     *
     * v is located at horizontal MAC faces.
     */
    for (
        std::size_t j = 1;
        j < ny;
        ++j
    ) {
        for (
            std::size_t i = 0;
            i < nx;
            ++i
        ) {

            const Real v_center =
                v(i, j);

            const Real u_at_v =
                Real{0.25} *
                (
                    u(i, j - 1) +
                    u(i + 1, j - 1) +
                    u(i, j) +
                    u(i + 1, j)
                );

            const Real west =
                i > 0
                ? v(i - 1, j)
                : tangential_ghost_value(
                    v_center,
                    boundaries.left,
                    boundaries.left.velocity.y
                );

            const Real east =
                i + 1 < nx
                ? v(i + 1, j)
                : tangential_ghost_value(
                    v_center,
                    boundaries.right,
                    boundaries.right.velocity.y
                );

            const Real dv_dx =
                u_at_v >= Real{}
                ?
                (
                    v_center -
                    west
                ) *
                inv_dx
                :
                (
                    east -
                    v_center
                ) *
                inv_dx;

            const Real dv_dy =
                v_center >= Real{}
                ?
                (
                    v_center -
                    v(i, j - 1)
                ) *
                inv_dy
                :
                (
                    v(i, j + 1) -
                    v_center
                ) *
                inv_dy;

            const Real laplacian_v =
                (
                    west -
                    Real{2} * v_center +
                    east
                ) *
                inv_dx2
                +
                (
                    v(i, j - 1) -
                    Real{2} * v_center +
                    v(i, j + 1)
                ) *
                inv_dy2;

            const Real advection =
                u_at_v * dv_dx +
                v_center * dv_dy;

            v_star(i, j) =
                v_center +
                dt *
                (
                    -advection +
                    nu * laplacian_v
                );
        }
    }

    apply_velocity_boundaries(
        predicted,
        grid,
        boundaries
    );
}

void build_pressure_rhs(
    const Problem2D& problem,
    const MacVelocityField2D& predicted,
    Field2D& rhs,
    const Real dt
) {

    const auto& grid =
        problem.grid();

    const auto& u =
        predicted.u();

    const auto& v =
        predicted.v();

    const Real scale =
        problem.fluid().density /
        dt;

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

            const Real divergence =
                (
                    u(i + 1, j) -
                    u(i, j)
                ) /
                grid.dx()
                +
                (
                    v(i, j + 1) -
                    v(i, j)
                ) /
                grid.dy();

            rhs(i, j) =
                scale *
                divergence;
        }
    }
}

void project_velocity(
    const Problem2D& problem,
    const MacVelocityField2D& predicted,
    const Field2D& pressure,
    MacVelocityField2D& result,
    const Real dt
) {

    const auto& grid =
        problem.grid();

    const Real scale =
        dt /
        problem.fluid().density;

    result =
        predicted;

    auto& u =
        result.u();

    auto& v =
        result.v();

    /*
     * Correct u.
     */
    for (
        std::size_t j = 0;
        j < grid.ny_cells();
        ++j
    ) {
        for (
            std::size_t i = 1;
            i < grid.nx_cells();
            ++i
        ) {

            const Real dp_dx =
                (
                    pressure(i, j) -
                    pressure(i - 1, j)
                ) /
                grid.dx();

            u(i, j) -=
                scale *
                dp_dx;
        }
    }

    /*
     * Correct v.
     */
    for (
        std::size_t j = 1;
        j < grid.ny_cells();
        ++j
    ) {
        for (
            std::size_t i = 0;
            i < grid.nx_cells();
            ++i
        ) {

            const Real dp_dy =
                (
                    pressure(i, j) -
                    pressure(i, j - 1)
                ) /
                grid.dy();

            v(i, j) -=
                scale *
                dp_dy;
        }
    }

    apply_velocity_boundaries(
        result,
        grid,
        problem.boundaries()
    );
}

[[nodiscard]]
Real compute_divergence_l2(
    const UniformGrid2D& grid,
    const MacVelocityField2D& velocity
) {

    const auto& u =
        velocity.u();

    const auto& v =
        velocity.v();

    Real squared_sum{};

    const auto cell_count =
        grid.nx_cells() *
        grid.ny_cells();

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

            const Real divergence =
                (
                    u(i + 1, j) -
                    u(i, j)
                ) /
                grid.dx()
                +
                (
                    v(i, j + 1) -
                    v(i, j)
                ) /
                grid.dy();

            squared_sum +=
                divergence *
                divergence;
        }
    }

    return std::sqrt(
        squared_sum /
        static_cast<Real>(
            cell_count
        )
    );
}

[[nodiscard]]
Real compute_cfl(
    const UniformGrid2D& grid,
    const MacVelocityField2D& velocity,
    const Real dt
) {

    Real max_u{};
    Real max_v{};

    for (
        const Real value :
        velocity.u().values()
    ) {
        max_u =
            std::max(
                max_u,
                std::abs(value)
            );
    }

    for (
        const Real value :
        velocity.v().values()
    ) {
        max_v =
            std::max(
                max_v,
                std::abs(value)
            );
    }

    return
        dt *
        (
            max_u / grid.dx() +
            max_v / grid.dy()
        );
}

} // namespace

CpuProjectionSolver2D::
CpuProjectionSolver2D(
    ProjectionSolverConfig config
)
    : config_{
        std::move(config)
      } {

    config_.validate();

switch (
    config_.pressure.kind
) {

case PressureSolverKind::jacobi:
case PressureSolverKind::weighted_jacobi:

    pressure_solver_ =
        std::make_unique<
            JacobiPressureSolver2D
        >();

    break;

case PressureSolverKind::spectral_dct:

    pressure_solver_ =
        std::make_unique<
            SpectralPressureSolver2D
        >();

    break;

case PressureSolverKind::
    conjugate_gradient:

case PressureSolverKind::
    multigrid:

    throw std::invalid_argument{
        "Selected pressure solver "
        "is not implemented yet"
    };
}
}

std::string_view
CpuProjectionSolver2D::
name() const noexcept {

    return "cpu_projection_2d";
}

void CpuProjectionSolver2D::
initialize(
    const Problem2D& problem
) {

    reject_periodic_boundaries(
        problem.boundaries()
    );

    problem_ =
        &problem;

    state_.emplace(
        problem.make_initial_state()
    );

    predicted_velocity_.emplace(
        problem.grid()
    );

    pressure_rhs_ =
        Field2D{
            problem.grid().nx_cells(),
            problem.grid().ny_cells()
        };

    apply_velocity_boundaries(
        state_->velocity(),
        problem.grid(),
        problem.boundaries()
    );
}

StepDiagnostics
CpuProjectionSolver2D::
advance(
    const Real dt
) {

    if (
        problem_ == nullptr ||
        !state_.has_value() ||
        !predicted_velocity_.has_value()
    ) {
        throw std::logic_error{
            "CpuProjectionSolver2D must "
            "be initialized before advance"
        };
    }

    if (dt <= Real{}) {
        throw std::invalid_argument{
            "Time step must be positive"
        };
    }

    /*
     * Step 1:
     * compute intermediate velocity u*.
     */
    compute_predictor(
        *problem_,
        state_->velocity(),
        *predicted_velocity_,
        dt
    );

    /*
     * Step 2:
     * build pressure Poisson RHS.
     */
    build_pressure_rhs(
        *problem_,
        *predicted_velocity_,
        pressure_rhs_,
        dt
    );

    /*
     * Step 3:
     * solve pressure Poisson equation.
     */
    const auto pressure_result =
        pressure_solver_->solve(
            problem_->grid(),
            problem_->boundaries(),
            pressure_rhs_,
            state_->pressure(),
            config_.pressure
        );

    /*
     * Step 4:
     * project velocity onto the
     * divergence-free field.
     */
    project_velocity(
        *problem_,
        *predicted_velocity_,
        state_->pressure(),
        state_->velocity(),
        dt
    );

    state_->set_clock(
        state_->time() + dt,
        state_->step() + 1
    );

    return {
        .pressure_iterations =
            pressure_result.iterations,

        .pressure_residual_l2 =
            pressure_result.residual_l2,

        .divergence_l2 =
            compute_divergence_l2(
                problem_->grid(),
                state_->velocity()
            ),

        .cfl =
            compute_cfl(
                problem_->grid(),
                state_->velocity(),
                dt
            ),

        .pressure_converged =
            pressure_result.converged
    };
}

void CpuProjectionSolver2D::
copy_state_to(
    FlowState2D& destination
) const {

    if (!state_.has_value()) {
        throw std::logic_error{
            "CpuProjectionSolver2D "
            "is not initialized"
        };
    }

    destination =
        *state_;
}

Real CpuProjectionSolver2D::
time() const noexcept {

    return
        state_.has_value()
        ? state_->time()
        : Real{};
}

std::size_t
CpuProjectionSolver2D::
step() const noexcept {

    return
        state_.has_value()
        ?
        static_cast<std::size_t>(
            state_->step()
        )
        :
        std::size_t{};
}

} // namespace nssim