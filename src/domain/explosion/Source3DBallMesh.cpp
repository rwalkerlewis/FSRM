/**
 * @file Source3DBallMesh.cpp
 * @brief Pass-13a (axis-1b foundation) implementation of Source3DBallMesh.
 *        See include/domain/explosion/Source3DBallMesh.hpp for the
 *        contract and design rationale.
 *
 * Implementation notes:
 *
 *  - File parsing tolerates `#`-prefixed comment lines and blank lines
 *    anywhere. TetGen 1.6 emits comment-free output but
 *    build_source_ball_mesh.py prepends a header line for traceability,
 *    so the parser must skip them.
 *  - TetGen indexing in .node / .ele files defaults to 1-based. The
 *    first vertex/cell ID is consulted to detect 0-based files; the
 *    parser stores 0-based internally either way.
 *  - DMPlex construction follows the GmshIO.cpp pattern in this
 *    repository: rank 0 owns the cell list and vertex coords, other
 *    ranks pass empty arrays, then DMPlexDistribute partitions.
 *  - DMPlexCreateFromCellListPetsc with interpolate=PETSC_TRUE builds
 *    the full Hasse diagram (faces and edges in addition to cells and
 *    vertices), so pass-13c can label boundary facets via
 *    DMPlexGetDepthStratum.
 *  - Vertex markers are written into a DMLabel before
 *    DMPlexDistribute. PETSc preserves DMLabels across distribution,
 *    so the post-distribute DM still carries the cavity / elastic
 *    surface tagging.
 */

#include "domain/explosion/Source3DBallMesh.hpp"

#include <petscviewer.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace FSRM {

namespace
{

constexpr int kDim = 3;
constexpr int kCornersPerTet = 4;

/// Strip a trailing `\r` if present (.node/.ele files written on
/// Windows lines round-trip with CRLF).
inline void stripCarriageReturn(std::string& s)
{
    if (!s.empty() && s.back() == '\r') {
        s.pop_back();
    }
}

/// Read the next non-comment, non-blank line from `in`. Returns true if
/// a line was read; false at EOF. Lines beginning with `#` are skipped.
/// Trailing comments after `#` on a data line are not stripped because
/// TetGen's own writer never produces them.
bool readNextDataLine(std::istream& in, std::string& line)
{
    while (std::getline(in, line)) {
        stripCarriageReturn(line);
        // Skip blank lines.
        bool blank = true;
        for (char c : line) {
            if (!std::isspace(static_cast<unsigned char>(c))) {
                blank = false;
                break;
            }
        }
        if (blank) continue;
        // Skip comment lines.
        size_t first_non_ws = 0;
        while (first_non_ws < line.size()
               && std::isspace(static_cast<unsigned char>(line[first_non_ws]))) {
            ++first_non_ws;
        }
        if (first_non_ws < line.size() && line[first_non_ws] == '#') {
            continue;
        }
        return true;
    }
    return false;
}

void parseNodeFile(const std::string& path,
                   std::vector<double>& coords_out,
                   std::vector<int>& markers_out,
                   int& base_index_out)
{
    std::ifstream in(path);
    if (!in.good()) {
        throw std::runtime_error(
            "Source3DBallMesh: cannot open TetGen .node file: " + path);
    }

    std::string line;
    if (!readNextDataLine(in, line)) {
        throw std::runtime_error(
            "Source3DBallMesh: empty .node file: " + path);
    }
    std::istringstream header(line);
    int n_nodes = 0;
    int dim = 0;
    int n_attrs = 0;
    int n_markers = 0;
    if (!(header >> n_nodes >> dim >> n_attrs >> n_markers)) {
        throw std::runtime_error(
            "Source3DBallMesh: malformed .node header: " + path);
    }
    if (dim != kDim) {
        throw std::runtime_error(
            "Source3DBallMesh: .node file dimension is " + std::to_string(dim)
            + ", expected 3: " + path);
    }
    if (n_nodes <= 0) {
        throw std::runtime_error(
            "Source3DBallMesh: .node file has zero nodes: " + path);
    }

    coords_out.assign(static_cast<size_t>(n_nodes) * 3, 0.0);
    markers_out.assign(static_cast<size_t>(n_nodes), 0);
    base_index_out = -1;

    for (int read = 0; read < n_nodes; ++read) {
        if (!readNextDataLine(in, line)) {
            throw std::runtime_error(
                "Source3DBallMesh: .node file truncated at node "
                + std::to_string(read) + ": " + path);
        }
        std::istringstream rs(line);
        int id = 0;
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
        if (!(rs >> id >> x >> y >> z)) {
            throw std::runtime_error(
                "Source3DBallMesh: malformed .node line " + std::to_string(read)
                + ": " + path);
        }
        for (int a = 0; a < n_attrs; ++a) {
            double dummy = 0.0;
            rs >> dummy;
        }
        int marker = 0;
        if (n_markers > 0) {
            rs >> marker;
        }
        if (base_index_out < 0) {
            base_index_out = id;
        }
        int idx0 = id - base_index_out;
        if (idx0 < 0 || idx0 >= n_nodes) {
            throw std::runtime_error(
                "Source3DBallMesh: .node id out of range: id=" + std::to_string(id)
                + ", base=" + std::to_string(base_index_out)
                + ", nNodes=" + std::to_string(n_nodes));
        }
        coords_out[static_cast<size_t>(idx0) * 3 + 0] = x;
        coords_out[static_cast<size_t>(idx0) * 3 + 1] = y;
        coords_out[static_cast<size_t>(idx0) * 3 + 2] = z;
        markers_out[static_cast<size_t>(idx0)] = marker;
    }
}

void parseEleFile(const std::string& path,
                  std::vector<int>& cells_out,
                  int base_index_in,
                  int n_nodes_in)
{
    std::ifstream in(path);
    if (!in.good()) {
        throw std::runtime_error(
            "Source3DBallMesh: cannot open TetGen .ele file: " + path);
    }

    std::string line;
    if (!readNextDataLine(in, line)) {
        throw std::runtime_error(
            "Source3DBallMesh: empty .ele file: " + path);
    }
    std::istringstream header(line);
    int n_tets = 0;
    int n_corners = 0;
    int n_attrs = 0;
    if (!(header >> n_tets >> n_corners >> n_attrs)) {
        throw std::runtime_error(
            "Source3DBallMesh: malformed .ele header: " + path);
    }
    if (n_corners != kCornersPerTet) {
        throw std::runtime_error(
            "Source3DBallMesh: .ele cell type has " + std::to_string(n_corners)
            + " corners; pass-13a requires linear tets (4): " + path);
    }
    if (n_tets <= 0) {
        throw std::runtime_error(
            "Source3DBallMesh: .ele file has zero cells: " + path);
    }

    cells_out.assign(static_cast<size_t>(n_tets) * kCornersPerTet, 0);

    int ele_base = -1;
    for (int read = 0; read < n_tets; ++read) {
        if (!readNextDataLine(in, line)) {
            throw std::runtime_error(
                "Source3DBallMesh: .ele file truncated at cell "
                + std::to_string(read) + ": " + path);
        }
        std::istringstream rs(line);
        int id = 0;
        int n0 = 0;
        int n1 = 0;
        int n2 = 0;
        int n3 = 0;
        if (!(rs >> id >> n0 >> n1 >> n2 >> n3)) {
            throw std::runtime_error(
                "Source3DBallMesh: malformed .ele line "
                + std::to_string(read) + ": " + path);
        }
        for (int a = 0; a < n_attrs; ++a) {
            double dummy = 0.0;
            rs >> dummy;
        }
        if (ele_base < 0) {
            ele_base = id;
        }
        int cell0 = id - ele_base;
        if (cell0 < 0 || cell0 >= n_tets) {
            throw std::runtime_error(
                "Source3DBallMesh: .ele cell id out of range: id="
                + std::to_string(id));
        }
        const int verts[4] = { n0, n1, n2, n3 };
        for (int k = 0; k < 4; ++k) {
            int idx0 = verts[k] - base_index_in;
            if (idx0 < 0 || idx0 >= n_nodes_in) {
                throw std::runtime_error(
                    "Source3DBallMesh: .ele references vertex out of range: id="
                    + std::to_string(verts[k]));
            }
            cells_out[static_cast<size_t>(cell0) * kCornersPerTet + k] = idx0;
        }
    }
}

}  // namespace

Source3DBallMesh::Source3DBallMesh() = default;

Source3DBallMesh::~Source3DBallMesh()
{
    if (dm_) {
        DMDestroy(&dm_);
        dm_ = nullptr;
    }
}

void Source3DBallMesh::loadFromTetGen(MPI_Comm comm, const std::string& basename)
{
    if (dm_) {
        throw std::runtime_error(
            "Source3DBallMesh::loadFromTetGen called on already-loaded mesh.");
    }
    comm_ = comm;
    source_basename_ = basename;

    PetscMPIInt rank = 0;
    MPI_Comm_rank(comm, &rank);

    std::vector<double> coords;
    std::vector<int> vertex_markers;
    std::vector<int> cells_flat;
    int n_nodes_global = 0;
    int n_cells_global = 0;

    if (rank == 0) {
        const std::string node_path = basename + ".node";
        const std::string ele_path = basename + ".ele";
        int base_index = 1;
        parseNodeFile(node_path, coords, vertex_markers, base_index);
        n_nodes_global = static_cast<int>(vertex_markers.size());
        parseEleFile(ele_path, cells_flat, base_index, n_nodes_global);
        n_cells_global = static_cast<int>(cells_flat.size() / kCornersPerTet);
    }

    int header_buf[2] = { n_nodes_global, n_cells_global };
    MPI_Bcast(header_buf, 2, MPI_INT, 0, comm);
    n_nodes_global = header_buf[0];
    n_cells_global = header_buf[1];
    num_global_cells_ = static_cast<PetscInt>(n_cells_global);

    PetscInt local_n_cells = (rank == 0) ? n_cells_global : 0;
    PetscInt local_n_vertices = (rank == 0) ? n_nodes_global : 0;
    const PetscInt* cell_ptr = nullptr;

    std::vector<PetscInt> cells_petsc;
    if (rank == 0) {
        cells_petsc.assign(cells_flat.begin(), cells_flat.end());
        cell_ptr = cells_petsc.data();
    }

    PetscErrorCode ierr = 0;
    DM dm_init = nullptr;
    if (rank == 0) {
        ierr = DMPlexCreateFromCellListPetsc(
            comm, kDim,
            local_n_cells, local_n_vertices, kCornersPerTet,
            PETSC_TRUE,
            cell_ptr,
            kDim,
            coords.data(),
            &dm_init);
    } else {
        ierr = DMPlexCreateFromCellListPetsc(
            comm, kDim,
            0, 0, kCornersPerTet,
            PETSC_TRUE,
            nullptr,
            kDim,
            nullptr,
            &dm_init);
    }
    if (ierr) {
        throw std::runtime_error(
            "Source3DBallMesh: DMPlexCreateFromCellListPetsc failed (rank "
            + std::to_string(rank) + ")");
    }

    if (rank == 0 && !vertex_markers.empty()) {
        ierr = DMCreateLabel(dm_init, kSourceBallVertexMarkerLabel);
        if (ierr) {
            throw std::runtime_error(
                "Source3DBallMesh: DMCreateLabel failed");
        }
        DMLabel pre_label = nullptr;
        ierr = DMGetLabel(dm_init, kSourceBallVertexMarkerLabel, &pre_label);
        if (ierr || !pre_label) {
            throw std::runtime_error(
                "Source3DBallMesh: DMGetLabel returned null");
        }
        PetscInt vStart = 0;
        PetscInt vEnd = 0;
        ierr = DMPlexGetDepthStratum(dm_init, 0, &vStart, &vEnd);
        if (ierr) {
            throw std::runtime_error(
                "Source3DBallMesh: DMPlexGetDepthStratum vertices failed");
        }
        const PetscInt expected_local = static_cast<PetscInt>(vertex_markers.size());
        if ((vEnd - vStart) != expected_local) {
            throw std::runtime_error(
                "Source3DBallMesh: vertex count mismatch after DMPlex create");
        }
        for (PetscInt v = vStart; v < vEnd; ++v) {
            const int marker = vertex_markers[static_cast<size_t>(v - vStart)];
            if (marker != 0) {
                ierr = DMLabelSetValue(pre_label, v, marker);
                if (ierr) {
                    throw std::runtime_error(
                        "Source3DBallMesh: DMLabelSetValue failed");
                }
            }
        }
    } else {
        // Other ranks must create the same label so DMPlexDistribute can
        // migrate it. The label is empty on these ranks.
        ierr = DMCreateLabel(dm_init, kSourceBallVertexMarkerLabel);
        if (ierr) {
            throw std::runtime_error(
                "Source3DBallMesh: DMCreateLabel (non-root) failed");
        }
    }

    DM dm_dist = nullptr;
    ierr = DMPlexDistribute(dm_init, 0, nullptr, &dm_dist);
    if (ierr) {
        throw std::runtime_error("Source3DBallMesh: DMPlexDistribute failed");
    }
    if (dm_dist) {
        DMDestroy(&dm_init);
        dm_ = dm_dist;
    } else {
        dm_ = dm_init;
    }

    PetscInt cStart = 0;
    PetscInt cEnd = 0;
    ierr = DMPlexGetHeightStratum(dm_, 0, &cStart, &cEnd);
    if (ierr) {
        throw std::runtime_error(
            "Source3DBallMesh: DMPlexGetHeightStratum cells failed");
    }
    num_local_cells_ = cEnd - cStart;

    PetscInt vStart = 0;
    PetscInt vEnd = 0;
    ierr = DMPlexGetDepthStratum(dm_, 0, &vStart, &vEnd);
    if (ierr) {
        throw std::runtime_error(
            "Source3DBallMesh: DMPlexGetDepthStratum vertices failed (post-distribute)");
    }
    num_local_vertices_ = vEnd - vStart;

    DMLabel post_label = nullptr;
    ierr = DMGetLabel(dm_, kSourceBallVertexMarkerLabel, &post_label);
    if (ierr) {
        post_label = nullptr;
    }
    vertex_marker_label_ = post_label;
    num_cavity_marked_vertices_ = 0;
    num_elastic_marked_vertices_ = 0;
    if (post_label) {
        for (PetscInt v = vStart; v < vEnd; ++v) {
            PetscInt val = 0;
            DMLabelGetValue(post_label, v, &val);
            if (val == SourceBallVertexMarkerValues::INNER_CAVITY_SURFACE) {
                ++num_cavity_marked_vertices_;
            } else if (val == SourceBallVertexMarkerValues::OUTER_ELASTIC_SURFACE) {
                ++num_elastic_marked_vertices_;
            }
        }
    }

    Vec coord_vec = nullptr;
    ierr = DMGetCoordinatesLocal(dm_, &coord_vec);
    if (ierr || coord_vec == nullptr) {
        local_min_radius_ = 0.0;
        local_max_radius_ = 0.0;
    } else {
        const PetscScalar* coord_data = nullptr;
        VecGetArrayRead(coord_vec, &coord_data);
        PetscInt coord_size = 0;
        VecGetLocalSize(coord_vec, &coord_size);
        const PetscInt n_coords = coord_size / kDim;
        double rmin = std::numeric_limits<double>::max();
        double rmax = 0.0;
        for (PetscInt i = 0; i < n_coords; ++i) {
            const double x = static_cast<double>(coord_data[kDim * i + 0]);
            const double y = static_cast<double>(coord_data[kDim * i + 1]);
            const double z = static_cast<double>(coord_data[kDim * i + 2]);
            const double r = std::sqrt(x * x + y * y + z * z);
            if (r < rmin) rmin = r;
            if (r > rmax) rmax = r;
        }
        VecRestoreArrayRead(coord_vec, &coord_data);
        if (n_coords == 0) {
            local_min_radius_ = 0.0;
            local_max_radius_ = 0.0;
        } else {
            local_min_radius_ = rmin;
            local_max_radius_ = rmax;
        }
    }
}

}  // namespace FSRM
