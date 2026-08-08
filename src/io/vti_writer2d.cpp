#include "nssim/io/vti_writer2d.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace nssim {

namespace {

void validate_state_shape(
    const UniformGrid2D& grid,
    const FlowState2D& state) {

    const auto nx =
        grid.nx_cells();

    const auto ny =
        grid.ny_cells();

    const auto& pressure =
        state.pressure();

    const auto& u =
        state.velocity().u();

    const auto& v =
        state.velocity().v();

    if (
        pressure.nx() != nx ||
        pressure.ny() != ny ||

        u.nx() != nx + 1 ||
        u.ny() != ny ||

        v.nx() != nx ||
        v.ny() != ny + 1
    ) {
        throw std::invalid_argument{
            "FlowState2D dimensions do not "
            "match UniformGrid2D"
        };
    }
}

} // namespace

void VtiWriter2D::write(
    const Problem2D& problem,
    const FlowState2D& state,
    const std::filesystem::path& file_path
) const {

    const auto& grid =
        problem.grid();

    validate_state_shape(
        grid,
        state
    );

    if (file_path.has_parent_path()) {
        std::filesystem::create_directories(
            file_path.parent_path()
        );
    }

    std::ofstream output{
        file_path
    };

    if (!output) {
        throw std::runtime_error{
            "Failed to open VTI output file: "
            + file_path.string()
        };
    }

    output
        << std::setprecision(17);

    const auto nx =
        grid.nx_cells();

    const auto ny =
        grid.ny_cells();

    const auto x_max =
        nx - 1;

    const auto y_max =
        ny - 1;

    output
        << "<?xml version=\"1.0\"?>\n"

        << "<VTKFile "
        << "type=\"ImageData\" "
        << "version=\"0.1\" "
        << "byte_order=\"LittleEndian\">\n"

        << "  <ImageData "
        << "WholeExtent=\"0 "
        << x_max
        << " 0 "
        << y_max
        << " 0 0\" "

        << "Origin=\""
        << Real{0.5} * grid.dx()
        << ' '
        << Real{0.5} * grid.dy()
        << " 0\" "

        << "Spacing=\""
        << grid.dx()
        << ' '
        << grid.dy()
        << " 1\">\n"

        << "    <Piece Extent=\"0 "
        << x_max
        << " 0 "
        << y_max
        << " 0 0\">\n"

        << "      <PointData "
        << "Scalars=\"pressure\" "
        << "Vectors=\"velocity\">\n";

    /*
     * Pressure
     */
    output
        << "        <DataArray "
        << "type=\"Float64\" "
        << "Name=\"pressure\" "
        << "NumberOfComponents=\"1\" "
        << "format=\"ascii\">\n";

    const auto& pressure =
        state.pressure();

    for (
        std::size_t j = 0;
        j < ny;
        ++j
    ) {
        output << "          ";

        for (
            std::size_t i = 0;
            i < nx;
            ++i
        ) {
            output
                << pressure(i, j)
                << ' ';
        }

        output << '\n';
    }

    output
        << "        </DataArray>\n";

    /*
     * Velocity
     *
     * MAC velocity is converted
     * to cell centers for visualization.
     */
    output
        << "        <DataArray "
        << "type=\"Float64\" "
        << "Name=\"velocity\" "
        << "NumberOfComponents=\"3\" "
        << "format=\"ascii\">\n";

    const auto& u =
        state.velocity().u();

    const auto& v =
        state.velocity().v();

    for (
        std::size_t j = 0;
        j < ny;
        ++j
    ) {
        output << "          ";

        for (
            std::size_t i = 0;
            i < nx;
            ++i
        ) {
            const Real u_center =
                Real{0.5} *
                (
                    u(i, j) +
                    u(i + 1, j)
                );

            const Real v_center =
                Real{0.5} *
                (
                    v(i, j) +
                    v(i, j + 1)
                );

            output
                << u_center
                << ' '
                << v_center
                << " 0 ";
        }

        output << '\n';
    }

    output
        << "        </DataArray>\n"

        << "      </PointData>\n"
        << "      <CellData/>\n"

        << "    </Piece>\n"
        << "  </ImageData>\n"
        << "</VTKFile>\n";
}

} // namespace nssim