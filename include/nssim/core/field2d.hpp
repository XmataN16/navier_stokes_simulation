#pragma once

#include "nssim/core/types.hpp"

#include <algorithm>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <vector>

namespace nssim {

class Field2D final {
public:
    Field2D() = default;

    Field2D(
        const std::size_t nx,
        const std::size_t ny,
        const Real initial_value = Real{})
        : nx_{nx},
          ny_{ny},
          values_(nx * ny, initial_value) {

        if (nx == 0 || ny == 0) {
            throw std::invalid_argument{
                "Field2D dimensions must be non-zero"
            };
        }
    }

    [[nodiscard]]
    std::size_t nx() const noexcept {
        return nx_;
    }

    [[nodiscard]]
    std::size_t ny() const noexcept {
        return ny_;
    }

    [[nodiscard]]
    std::size_t size() const noexcept {
        return values_.size();
    }

    [[nodiscard]]
    bool empty() const noexcept {
        return values_.empty();
    }

    Real& operator()(
        const std::size_t i,
        const std::size_t j) noexcept {

        return values_[index(i, j)];
    }

    const Real& operator()(
        const std::size_t i,
        const std::size_t j) const noexcept {

        return values_[index(i, j)];
    }

    Real& at(
        const std::size_t i,
        const std::size_t j) {

        check_index(i, j);
        return values_[index(i, j)];
    }

    const Real& at(
        const std::size_t i,
        const std::size_t j) const {

        check_index(i, j);
        return values_[index(i, j)];
    }

    [[nodiscard]]
    std::span<Real> values() noexcept {
        return values_;
    }

    [[nodiscard]]
    std::span<const Real> values() const noexcept {
        return values_;
    }

    void fill(const Real value) noexcept {
        std::fill(
            values_.begin(),
            values_.end(),
            value
        );
    }

private:
    [[nodiscard]]
    std::size_t index(
        const std::size_t i,
        const std::size_t j) const noexcept {

        return j * nx_ + i;
    }

    void check_index(
        const std::size_t i,
        const std::size_t j) const {

        if (i >= nx_ || j >= ny_) {
            throw std::out_of_range{
                "Field2D index is out of range"
            };
        }
    }

    std::size_t nx_{};
    std::size_t ny_{};

    std::vector<Real> values_{};
};

} // namespace nssim