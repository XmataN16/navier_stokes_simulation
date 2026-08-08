#include "nssim/problem/problem2d.hpp"

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

} // namespace

Problem2D::Problem2D(
    UniformGrid2D grid,
    FluidProperties fluid,
    BoundarySet2D boundaries,
    InitialCondition2D initial_condition
)
    : grid_{std::move(grid)},
      fluid_{fluid},
      boundaries_{boundaries},
      initial_condition_{
          std::move(initial_condition)
      } {

    validate();
}

void Problem2D::validate() const {

    fluid_.validate();

    if (
        is_periodic(boundaries_.left) !=
        is_periodic(boundaries_.right)
    ) {
        throw std::invalid_argument{
            "Periodic left/right boundaries "
            "must be specified as a pair"
        };
    }

    if (
        is_periodic(boundaries_.bottom) !=
        is_periodic(boundaries_.top)
    ) {
        throw std::invalid_argument{
            "Periodic bottom/top boundaries "
            "must be specified as a pair"
        };
    }
}

FlowState2D
Problem2D::make_initial_state() const {

    FlowState2D state{grid_};

    /*
     * Initialize pressure at cell centers.
     */
    auto& pressure =
        state.pressure();

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
            pressure(i, j) =
                initial_condition_.pressure(
                    grid_.cell_center(i, j)
                );
        }
    }

    /*
     * Initialize x velocity component
     * at vertical MAC faces.
     */
    auto& u =
        state.velocity().u();

    for (
        std::size_t j = 0;
        j < u.ny();
        ++j
    ) {
        for (
            std::size_t i = 0;
            i < u.nx();
            ++i
        ) {
            u(i, j) =
                initial_condition_.velocity(
                    grid_.u_face_center(i, j)
                ).x;
        }
    }

    /*
     * Initialize y velocity component
     * at horizontal MAC faces.
     */
    auto& v =
        state.velocity().v();

    for (
        std::size_t j = 0;
        j < v.ny();
        ++j
    ) {
        for (
            std::size_t i = 0;
            i < v.nx();
            ++i
        ) {
            v(i, j) =
                initial_condition_.velocity(
                    grid_.v_face_center(i, j)
                ).y;
        }
    }

    return state;
}

} // namespace nssim