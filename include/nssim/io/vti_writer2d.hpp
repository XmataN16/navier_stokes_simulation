#pragma once

#include "nssim/io/result_writer2d.hpp"

namespace nssim {

class VtiWriter2D final
    : public IResultWriter2D {

public:
    void write(
        const Problem2D& problem,
        const FlowState2D& state,
        const std::filesystem::path& file_path
    ) const override;
};

} // namespace nssim