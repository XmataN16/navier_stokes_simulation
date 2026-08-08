#include "nssim/io/vtk_time_series_writer2d.hpp"

#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace nssim {

namespace {

[[nodiscard]]
std::string escape_xml_attribute(
    const std::string& value
) {

    std::string result;

    result.reserve(
        value.size()
    );

    for (
        const char character :
        value
    ) {

        switch (character) {

        case '&':
            result += "&amp;";
            break;

        case '"':
            result += "&quot;";
            break;

        case '<':
            result += "&lt;";
            break;

        case '>':
            result += "&gt;";
            break;

        default:
            result += character;
            break;
        }
    }

    return result;
}

} // namespace

VtkTimeSeriesWriter2D::
VtkTimeSeriesWriter2D(
    std::filesystem::path output_directory,
    std::string series_name
)
    : output_directory_{
          std::move(
              output_directory
          )
      },

      series_name_{
          std::move(
              series_name
          )
      },

      collection_path_{
          output_directory_ /
          (
              series_name_ +
              ".pvd"
          )
      } {

    if (series_name_.empty()) {
        throw std::invalid_argument{
            "VTK series name "
            "must not be empty"
        };
    }
}

void VtkTimeSeriesWriter2D::
reset() {

    entries_.clear();

    std::filesystem::
        create_directories(
            output_directory_
        );

    /*
     * Remove stale frames from previous run,
     * but only frames belonging to this series.
     */
    const std::string frame_prefix =
        series_name_ +
        "_";

    for (
        const auto& entry :
        std::filesystem::
            directory_iterator{
                output_directory_
            }
    ) {

        if (
            !entry.is_regular_file()
        ) {
            continue;
        }

        const auto file_name =
            entry.path()
                .filename()
                .string();

        if (
            file_name.starts_with(
                frame_prefix
            ) &&
            entry.path().extension() ==
                ".vti"
        ) {

            std::filesystem::remove(
                entry.path()
            );
        }
    }

    /*
     * Create empty collection file.
     */
    write_collection();
}

void VtkTimeSeriesWriter2D::
write_frame(
    const Problem2D& problem,
    const FlowState2D& state
) {

    std::filesystem::
        create_directories(
            output_directory_
        );

    const auto frame_path =
        make_frame_path(
            state.step()
        );

    vti_writer_.write(
        problem,
        state,
        frame_path
    );

    entries_.push_back(
        Entry{
            .time =
                state.time(),

            .file_name =
                frame_path.filename()
        }
    );

    /*
     * Update PVD after every frame.
     *
     * This also leaves a usable collection
     * if simulation is interrupted later.
     */
    write_collection();
}

const std::filesystem::path&
VtkTimeSeriesWriter2D::
collection_path() const noexcept {

    return collection_path_;
}

std::filesystem::path
VtkTimeSeriesWriter2D::
make_frame_path(
    const std::uint64_t step
) const {

    std::ostringstream name;

    name
        << series_name_
        << '_'
        << std::setw(6)
        << std::setfill('0')
        << step
        << ".vti";

    return
        output_directory_ /
        name.str();
}

void VtkTimeSeriesWriter2D::
write_collection() const {

    std::ofstream output{
        collection_path_
    };

    if (!output) {
        throw std::runtime_error{
            "Failed to open PVD "
            "collection file: "
            + collection_path_.string()
        };
    }

    output
        << std::setprecision(17)

        << "<?xml version=\"1.0\"?>\n"

        << "<VTKFile "
        << "type=\"Collection\" "
        << "version=\"0.1\" "
        << "byte_order=\"LittleEndian\">\n"

        << "  <Collection>\n";

    for (
        const Entry& entry :
        entries_
    ) {

        output
            << "    <DataSet "

            << "timestep=\""
            << entry.time
            << "\" "

            << "group=\"\" "

            << "part=\"0\" "

            << "file=\""

            << escape_xml_attribute(
                entry.file_name.
                    generic_string()
            )

            << "\"/>\n";
    }

    output
        << "  </Collection>\n"
        << "</VTKFile>\n";
}

} // namespace nssim