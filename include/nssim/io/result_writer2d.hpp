#pragma once

#include "nssim/core/state2d.hpp"
#include "nssim/problem/problem2d.hpp"

#include <filesystem>

namespace nssim {

class IResultWriter2D {
public:
    virtual ~IResultWriter2D() = default;

    virtual void write(
        const Problem2D& problem,
        const FlowState2D& state,
        const std::filesystem::path& file_path
    ) const = 0;
};

} // namespace nssim