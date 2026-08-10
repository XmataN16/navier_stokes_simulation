#pragma once

#include "nssim/core/types.hpp"

#include <cstddef>
#include <stdexcept>

namespace nssim {

enum class PressureSolverKind {
    jacobi,
    weighted_jacobi,
    spectral_dct,
    conjugate_gradient,
    multigrid
};

struct PressureSolverConfig final {
    PressureSolverKind kind{
        PressureSolverKind::weighted_jacobi
    };

    std::size_t max_iterations{
        2000
    };

    Real tolerance{
        1.0e-6
    };

    Real relaxation{
        2.0 / 3.0
    };

    void validate() const {

        /*
         * Spectral DCT is a direct solver,
         * therefore max_iterations does not
         * apply to it.
         */
        if (
            kind != PressureSolverKind::spectral_dct &&
            max_iterations == 0
        ) {
            throw std::invalid_argument{
                "Pressure solver max_iterations "
                "must be greater than zero"
            };
        }

        if (tolerance <= Real{}) {
            throw std::invalid_argument{
                "Pressure solver tolerance "
                "must be positive"
            };
        }

        if (
            kind ==
                PressureSolverKind::weighted_jacobi &&
            (
                relaxation <= Real{} ||
                relaxation >= Real{1}
            )
        ) {
            throw std::invalid_argument{
                "Weighted Jacobi relaxation "
                "must be in (0, 1)"
            };
        }
    }
};

} // namespace nssim