#pragma once

#include "nssim/core/state2d.hpp"
#include "nssim/io/vti_writer2d.hpp"
#include "nssim/problem/problem2d.hpp"

#include <filesystem>
#include <string>
#include <vector>

namespace nssim {

class VtkTimeSeriesWriter2D final {
public:
    VtkTimeSeriesWriter2D(
        std::filesystem::path output_directory,
        std::string series_name
    );

    void reset();

    void write_frame(
        const Problem2D& problem,
        const FlowState2D& state
    );

    [[nodiscard]]
    const std::filesystem::path&
    collection_path() const noexcept;

private:
    struct Entry final {
        Real time{};

        std::filesystem::path
            file_name{};
    };

    void write_collection() const;

    [[nodiscard]]
    std::filesystem::path
    make_frame_path(
        std::uint64_t step
    ) const;

    std::filesystem::path
        output_directory_{};

    std::string
        series_name_{};

    std::filesystem::path
        collection_path_{};

    VtiWriter2D
        vti_writer_{};

    std::vector<Entry>
        entries_{};
};

} // namespace nssim