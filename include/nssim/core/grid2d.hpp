#pragma once

#include "nssim/core/types.hpp"

#include <cstddef>
#include <stdexcept>

namespace nssim {

class UniformGrid2D final {
public:
    UniformGrid2D(
        const std::size_t nx_cells,
        const std::size_t ny_cells,
        const Real length_x,
        const Real length_y)
        : nx_cells_{nx_cells},
          ny_cells_{ny_cells},
          length_x_{length_x},
          length_y_{length_y} {

        if (nx_cells < 2 || ny_cells < 2) {
            throw std::invalid_argument{
                "UniformGrid2D requires at least 2x2 cells"
            };
        }

        if (length_x <= Real{} ||
            length_y <= Real{}) {

            throw std::invalid_argument{
                "UniformGrid2D lengths must be positive"
            };
        }

        dx_ =
            length_x_ /
            static_cast<Real>(nx_cells_);

        dy_ =
            length_y_ /
            static_cast<Real>(ny_cells_);
    }

    [[nodiscard]]
    std::size_t nx_cells() const noexcept {
        return nx_cells_;
    }

    [[nodiscard]]
    std::size_t ny_cells() const noexcept {
        return ny_cells_;
    }

    [[nodiscard]]
    Real length_x() const noexcept {
        return length_x_;
    }

    [[nodiscard]]
    Real length_y() const noexcept {
        return length_y_;
    }

    [[nodiscard]]
    Real dx() const noexcept {
        return dx_;
    }

    [[nodiscard]]
    Real dy() const noexcept {
        return dy_;
    }

    [[nodiscard]]
    Vec2 cell_center(
        const std::size_t i,
        const std::size_t j) const noexcept {

        return {
            (static_cast<Real>(i) + Real{0.5}) * dx_,
            (static_cast<Real>(j) + Real{0.5}) * dy_
        };
    }

    [[nodiscard]]
    Vec2 u_face_center(
        const std::size_t i,
        const std::size_t j) const noexcept {

        return {
            static_cast<Real>(i) * dx_,
            (static_cast<Real>(j) + Real{0.5}) * dy_
        };
    }

    [[nodiscard]]
    Vec2 v_face_center(
        const std::size_t i,
        const std::size_t j) const noexcept {

        return {
            (static_cast<Real>(i) + Real{0.5}) * dx_,
            static_cast<Real>(j) * dy_
        };
    }

private:
    std::size_t nx_cells_{};
    std::size_t ny_cells_{};

    Real length_x_{};
    Real length_y_{};

    Real dx_{};
    Real dy_{};
};

} // namespace nssim