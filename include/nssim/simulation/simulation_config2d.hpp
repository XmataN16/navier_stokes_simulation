#pragma once

#include "nssim/core/types.hpp"

#include <cstddef>
#include <stdexcept>

namespace nssim {

struct SimulationConfig2D final {
    Real end_time{1.0};
    Real time_step{1.0e-3};

    std::size_t output_interval_steps{100};
    std::size_t diagnostics_interval_steps{10};

    Real maximum_cfl{0.8};

    bool fail_on_pressure_nonconvergence{true};

    void validate() const {

        if (end_time <= Real{}) {
            throw std::invalid_argument{
                "Simulation end_time must be positive"
            };
        }

        if (time_step <= Real{}) {
            throw std::invalid_argument{
                "Simulation time_step must be positive"
            };
        }

        if (output_interval_steps == 0) {
            throw std::invalid_argument{
                "Simulation output_interval_steps "
                "must be greater than zero"
            };
        }

        if (diagnostics_interval_steps == 0) {
            throw std::invalid_argument{
                "Simulation diagnostics_interval_steps "
                "must be greater than zero"
            };
        }

        if (maximum_cfl <= Real{}) {
            throw std::invalid_argument{
                "Simulation maximum_cfl must be positive"
            };
        }
    }
};

} // namespace nssim