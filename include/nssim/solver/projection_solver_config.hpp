#pragma once

#include "nssim/core/types.hpp"

#include <cstddef>
#include <stdexcept>

namespace nssim {

enum class PressureSolverKind {
    jacobi,
    weighted_jacobi,
    conjugate_gradient,
    multigrid
};

struct ProjectionSolverConfig final {

    PressureSolverKind pressure_solver{
        PressureSolverKind::weighted_jacobi
    };

    std::size_t max_pressure_iterations{
        500
    };

    Real pressure_tolerance{
        1.0e-8
    };

    Real jacobi_omega{
        2.0 / 3.0
    };

    void validate() const {

        if (max_pressure_iterations == 0) {
            throw std::invalid_argument{
                "max_pressure_iterations "
                "must be greater than zero"
            };
        }

        if (pressure_tolerance <= Real{}) {
            throw std::invalid_argument{
                "pressure_tolerance "
                "must be positive"
            };
        }

        if (
            jacobi_omega <= Real{} ||
            jacobi_omega >= Real{1}
        ) {
            throw std::invalid_argument{
                "jacobi_omega must be "
                "in the interval (0, 1)"
            };
        }
    }
};

} // namespace nssim