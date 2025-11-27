#include "GmshIO.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cstring>
#include <limits>
#include <set>

namespace FSRM {

// ============================================================================
// MeshElement Implementation
// ============================================================================

int MeshElement::getNumNodes() const {
    switch (type) {
        case GmshElementType::POINT: return 1;
        case GmshElementType::LINE: return 2;
        case GmshElementType::TRIANGLE: return 3;
        case GmshElementType::QUAD: return 4;
        case GmshElementType::TETRAHEDRON: return 4;
        case GmshElementType::HEXAHEDRON: return 8;
        case GmshElementType::PRISM: return 6;
        case GmshElementType::PYRAMID: return 5;
        case GmshElementType::LINE3: return 3;
        case GmshElementType::TRIANGLE6: return 6;
        case GmshElementType::QUAD9: return 9;
        case GmshElementType::TET10: return 10;
        case GmshElementType::HEX27: return 27;
        case GmshElementType::PRISM18: return 18;
        case GmshElementType::PYRAMID14: return 14;
        case GmshElementType::QUAD8: return 8;
        case GmshElementType::HEX20: return 20;
        case GmshElementType::PRISM15: return 15;
        case GmshElementType::PYRAMID13: return 13;
        default: return 0;
    }
}

int MeshElement::getDimension() const {
    switch (type) {
        case GmshElementType::POINT: return 0;
        case GmshElementType::LINE:
        case GmshElementType::LINE3: return 1;
        case GmshElementType::TRIANGLE:
        case GmshElementType::QUAD:
        case GmshElementType::TRIANGLE6:
        case GmshElementType::QUAD8:
        case GmshElementType::QUAD9: return 2;
        case GmshElementType::TETRAHEDRON:
        case GmshElementType::HEXAHEDRON:
        case GmshElementType::PRISM:
        case GmshElementType::PYRAMID:
        case GmshElementType::TET10:
        case GmshElementType::HEX20:
        case GmshElementType::HEX27:
        case GmshElementType::PRISM15:
        case GmshElementType::PRISM18:
        case GmshElementType::PYRAMID13:
        case GmshElementType::PYRAMID14: return 3;
        default: return -1;
    }
}

// ============================================================================
// UnstructuredMesh Implementation
// ============================================================================

int UnstructuredMesh::getNumCells() const {
    int count = 0;
    for (const auto& elem : elements) {
        if (elem.getDimension() == dimension) {
            count++;
        }
    }
    return count;
}

int UnstructuredMesh::getNumFaces() const {
    int count = 0;
    for (const auto& elem : elements) {
        if (elem.getDimension() == dimension - 1) {
            count++;
        }
    }
    return count;
}

void UnstructuredMesh::computeBoundingBox() {
    if (nodes.empty()) return;
    
    bbox_min = {std::numeric_limits<double>::max(),
                std::numeric_limits<double>::max(),
                std::numeric_limits<double>::max()};
    bbox_max = {std::numeric_limits<double>::lowest(),
                std::numeric_limits<double>::lowest(),
                std::numeric_limits<double>::lowest()};
    
    for (const auto& node : nodes) {
        bbox_min[0] = std::min(bbox_min[0], node.x);
        bbox_min[1] = std::min(bbox_min[1], node.y);
        bbox_min[2] = std::min(bbox_min[2], node.z);
        bbox_max[0] = std::max(bbox_max[0], node.x);
        bbox_max[1] = std::max(bbox_max[1], node.y);
        bbox_max[2] = std::max(bbox_max[2], node.z);
    }
}

std::vector<int> UnstructuredMesh::getElementsByPhysicalTag(int tag) const {
    std::vector<int> result;
    for (size_t i = 0; i < elements.size(); i++) {
        if (elements[i].physical_tag == tag) {
            result.push_back(static_cast<int>(i));
        }
    }
    return result;
}

std::vector<int> UnstructuredMesh::getNodesByPhysicalTag(int tag) const {
    std::vector<int> node_tags;
    std::vector<bool> node_found(nodes.size(), false);
    
    for (const auto& elem : elements) {
        if (elem.physical_tag == tag) {
            for (int node_tag : elem.node_tags) {
                // Find node index (tags may not be consecutive)
                for (size_t i = 0; i < nodes.size(); i++) {
                    if (nodes[i].tag == node_tag && !node_found[i]) {
                        node_found[i] = true;
                        node_tags.push_back(node_tag);
                        break;
                    }
                }
            }
        }
    }
    return node_tags;
}

int UnstructuredMesh::getPhysicalTagByName(const std::string& name) const {
    for (const auto& [tag, group] : physical_groups) {
        if (group.name == name) {
            return tag;
        }
    }
    return -1;  // Not found
}

// ============================================================================
// GmshIO Implementation
// ============================================================================

GmshIO::GmshIO() : mesh_(nullptr), format_version_(""), is_binary_(false) {
}

bool GmshIO::readMshFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "GmshIO: Cannot open file: " << filename << std::endl;
        return false;
    }
    
    mesh_ = std::make_shared<UnstructuredMesh>();
    
    std::string line;
    while (std::getline(file, line)) {
        // Remove carriage return if present (Windows line endings)
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        
        if (line == "$MeshFormat") {
            if (!readMeshFormat(file)) return false;
        }
        else if (line == "$PhysicalNames") {
            if (format_version_.substr(0, 1) == "4") {
                if (!readPhysicalNames(file)) return false;
            } else {
                if (!readPhysicalNames2(file)) return false;
            }
        }
        else if (line == "$Entities") {
            if (!readEntities(file)) return false;
        }
        else if (line == "$Nodes") {
            if (format_version_.substr(0, 1) == "4") {
                if (!readNodes(file)) return false;
            } else {
                if (!readNodes2(file)) return false;
            }
        }
        else if (line == "$Elements") {
            if (format_version_.substr(0, 1) == "4") {
                if (!readElements(file)) return false;
            } else {
                if (!readElements2(file)) return false;
            }
        }
    }
    
    file.close();
    
    // Post-processing
    mesh_->computeBoundingBox();
    
    // Determine mesh dimension
    int max_dim = 0;
    for (const auto& elem : mesh_->elements) {
        max_dim = std::max(max_dim, elem.getDimension());
    }
    mesh_->dimension = max_dim;
    
    return true;
}

bool GmshIO::readMeshFormat(std::ifstream& file) {
    std::string line;
    std::getline(file, line);
    
    std::istringstream iss(line);
    double version;
    int file_type, data_size;
    
    iss >> version >> file_type >> data_size;
    
    format_version_ = std::to_string(version);
    is_binary_ = (file_type == 1);
    
    if (is_binary_) {
        std::cerr << "GmshIO: Binary format not yet supported" << std::endl;
        return false;
    }
    
    // Read until end of section
    while (std::getline(file, line)) {
        if (line.find("$EndMeshFormat") != std::string::npos) break;
    }
    
    return true;
}

bool GmshIO::readPhysicalNames(std::ifstream& file) {
    std::string line;
    std::getline(file, line);
    
    int num_names;
    std::istringstream iss(line);
    iss >> num_names;
    
    for (int i = 0; i < num_names; i++) {
        std::getline(file, line);
        std::istringstream line_iss(line);
        
        int dim, tag;
        std::string name;
        line_iss >> dim >> tag >> name;
        
        // Remove quotes from name
        if (!name.empty() && name.front() == '"') {
            name = name.substr(1);
        }
        if (!name.empty() && name.back() == '"') {
            name.pop_back();
        }
        
        mesh_->physical_groups[tag] = PhysicalGroup(tag, dim, name);
    }
    
    // Read until end of section
    while (std::getline(file, line)) {
        if (line.find("$EndPhysicalNames") != std::string::npos) break;
    }
    
    return true;
}

bool GmshIO::readPhysicalNames2(std::ifstream& file) {
    // Same format as v4 for physical names
    return readPhysicalNames(file);
}

bool GmshIO::readEntities(std::ifstream& file) {
    std::string line;
    std::getline(file, line);
    
    std::istringstream iss(line);
    int num_points, num_curves, num_surfaces, num_volumes;
    iss >> num_points >> num_curves >> num_surfaces >> num_volumes;
    
    // Skip point entities
    for (int i = 0; i < num_points; i++) {
        std::getline(file, line);
    }
    
    // Skip curve entities
    for (int i = 0; i < num_curves; i++) {
        std::getline(file, line);
    }
    
    // Skip surface entities
    for (int i = 0; i < num_surfaces; i++) {
        std::getline(file, line);
    }
    
    // Skip volume entities
    for (int i = 0; i < num_volumes; i++) {
        std::getline(file, line);
    }
    
    // Read until end of section
    while (std::getline(file, line)) {
        if (line.find("$EndEntities") != std::string::npos) break;
    }
    
    return true;
}

bool GmshIO::readNodes(std::ifstream& file) {
    std::string line;
    std::getline(file, line);
    
    std::istringstream iss(line);
    int num_entity_blocks, num_nodes, min_tag, max_tag;
    iss >> num_entity_blocks >> num_nodes >> min_tag >> max_tag;
    
    mesh_->nodes.reserve(num_nodes);
    
    for (int block = 0; block < num_entity_blocks; block++) {
        std::getline(file, line);
        std::istringstream block_iss(line);
        
        int entity_dim, entity_tag, parametric, num_nodes_in_block;
        block_iss >> entity_dim >> entity_tag >> parametric >> num_nodes_in_block;
        
        // Read node tags
        std::vector<int> tags(num_nodes_in_block);
        for (int i = 0; i < num_nodes_in_block; i++) {
            std::getline(file, line);
            tags[i] = std::stoi(line);
        }
        
        // Read node coordinates
        for (int i = 0; i < num_nodes_in_block; i++) {
            std::getline(file, line);
            std::istringstream coord_iss(line);
            double x, y, z;
            coord_iss >> x >> y >> z;
            mesh_->nodes.emplace_back(tags[i], x, y, z);
        }
    }
    
    // Read until end of section
    while (std::getline(file, line)) {
        if (line.find("$EndNodes") != std::string::npos) break;
    }
    
    return true;
}

bool GmshIO::readNodes2(std::ifstream& file) {
    std::string line;
    std::getline(file, line);
    
    int num_nodes = std::stoi(line);
    mesh_->nodes.reserve(num_nodes);
    
    for (int i = 0; i < num_nodes; i++) {
        std::getline(file, line);
        std::istringstream iss(line);
        
        int tag;
        double x, y, z;
        iss >> tag >> x >> y >> z;
        
        mesh_->nodes.emplace_back(tag, x, y, z);
    }
    
    // Read until end of section
    while (std::getline(file, line)) {
        if (line.find("$EndNodes") != std::string::npos) break;
    }
    
    return true;
}

bool GmshIO::readElements(std::ifstream& file) {
    std::string line;
    std::getline(file, line);
    
    std::istringstream iss(line);
    int num_entity_blocks, num_elements, min_tag, max_tag;
    iss >> num_entity_blocks >> num_elements >> min_tag >> max_tag;
    
    mesh_->elements.reserve(num_elements);
    
    for (int block = 0; block < num_entity_blocks; block++) {
        std::getline(file, line);
        std::istringstream block_iss(line);
        
        int entity_dim, entity_tag, elem_type, num_elements_in_block;
        block_iss >> entity_dim >> entity_tag >> elem_type >> num_elements_in_block;
        
        GmshElementType type = intToElementType(elem_type);
        int num_nodes_per_elem = getNumNodesForElementType(type);
        
        for (int i = 0; i < num_elements_in_block; i++) {
            std::getline(file, line);
            std::istringstream elem_iss(line);
            
            MeshElement elem;
            elem_iss >> elem.tag;
            elem.type = type;
            elem.geometrical_tag = entity_tag;
            elem.physical_tag = entity_tag;  // May be updated later
            
            elem.node_tags.resize(num_nodes_per_elem);
            for (int n = 0; n < num_nodes_per_elem; n++) {
                elem_iss >> elem.node_tags[n];
            }
            
            mesh_->elements.push_back(elem);
        }
    }
    
    // Read until end of section
    while (std::getline(file, line)) {
        if (line.find("$EndElements") != std::string::npos) break;
    }
    
    return true;
}

bool GmshIO::readElements2(std::ifstream& file) {
    std::string line;
    std::getline(file, line);
    
    int num_elements = std::stoi(line);
    mesh_->elements.reserve(num_elements);
    
    for (int i = 0; i < num_elements; i++) {
        std::getline(file, line);
        std::istringstream iss(line);
        
        int tag, elem_type, num_tags;
        iss >> tag >> elem_type >> num_tags;
        
        MeshElement elem;
        elem.tag = tag;
        elem.type = intToElementType(elem_type);
        
        // Read tags (usually: physical_tag, geometrical_tag, ...)
        if (num_tags >= 1) {
            iss >> elem.physical_tag;
        }
        if (num_tags >= 2) {
            iss >> elem.geometrical_tag;
        }
        // Skip remaining tags
        for (int t = 2; t < num_tags; t++) {
            int dummy;
            iss >> dummy;
        }
        
        // Read node tags
        int num_nodes = getNumNodesForElementType(elem.type);
        elem.node_tags.resize(num_nodes);
        for (int n = 0; n < num_nodes; n++) {
            iss >> elem.node_tags[n];
        }
        
        mesh_->elements.push_back(elem);
    }
    
    // Read until end of section
    while (std::getline(file, line)) {
        if (line.find("$EndElements") != std::string::npos) break;
    }
    
    return true;
}

int GmshIO::getNumNodesForElementType(GmshElementType type) const {
    switch (type) {
        case GmshElementType::POINT: return 1;
        case GmshElementType::LINE: return 2;
        case GmshElementType::TRIANGLE: return 3;
        case GmshElementType::QUAD: return 4;
        case GmshElementType::TETRAHEDRON: return 4;
        case GmshElementType::HEXAHEDRON: return 8;
        case GmshElementType::PRISM: return 6;
        case GmshElementType::PYRAMID: return 5;
        case GmshElementType::LINE3: return 3;
        case GmshElementType::TRIANGLE6: return 6;
        case GmshElementType::QUAD9: return 9;
        case GmshElementType::TET10: return 10;
        case GmshElementType::HEX27: return 27;
        case GmshElementType::PRISM18: return 18;
        case GmshElementType::PYRAMID14: return 14;
        case GmshElementType::QUAD8: return 8;
        case GmshElementType::HEX20: return 20;
        case GmshElementType::PRISM15: return 15;
        case GmshElementType::PYRAMID13: return 13;
        default: return 0;
    }
}

GmshElementType GmshIO::intToElementType(int type_id) const {
    switch (type_id) {
        case 1: return GmshElementType::LINE;
        case 2: return GmshElementType::TRIANGLE;
        case 3: return GmshElementType::QUAD;
        case 4: return GmshElementType::TETRAHEDRON;
        case 5: return GmshElementType::HEXAHEDRON;
        case 6: return GmshElementType::PRISM;
        case 7: return GmshElementType::PYRAMID;
        case 8: return GmshElementType::LINE3;
        case 9: return GmshElementType::TRIANGLE6;
        case 10: return GmshElementType::QUAD9;
        case 11: return GmshElementType::TET10;
        case 12: return GmshElementType::HEX27;
        case 13: return GmshElementType::PRISM18;
        case 14: return GmshElementType::PYRAMID14;
        case 15: return GmshElementType::POINT;
        case 16: return GmshElementType::QUAD8;
        case 17: return GmshElementType::HEX20;
        case 18: return GmshElementType::PRISM15;
        case 19: return GmshElementType::PYRAMID13;
        default: return GmshElementType::UNKNOWN;
    }
}

bool GmshIO::writeMshFile(const std::string& filename) const {
    if (!mesh_ || mesh_->nodes.empty()) {
        std::cerr << "GmshIO: No mesh to write" << std::endl;
        return false;
    }
    
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "GmshIO: Cannot create file: " << filename << std::endl;
        return false;
    }
    
    // Write mesh format
    file << "$MeshFormat\n";
    file << "4.1 0 8\n";
    file << "$EndMeshFormat\n";
    
    // Write physical names
    if (!mesh_->physical_groups.empty()) {
        file << "$PhysicalNames\n";
        file << mesh_->physical_groups.size() << "\n";
        for (const auto& [tag, group] : mesh_->physical_groups) {
            file << group.dimension << " " << tag << " \"" << group.name << "\"\n";
        }
        file << "$EndPhysicalNames\n";
    }
    
    // Write nodes (simplified - single block)
    file << "$Nodes\n";
    file << "1 " << mesh_->nodes.size() << " 1 " << mesh_->nodes.size() << "\n";
    file << mesh_->dimension << " 1 0 " << mesh_->nodes.size() << "\n";
    
    for (const auto& node : mesh_->nodes) {
        file << node.tag << "\n";
    }
    for (const auto& node : mesh_->nodes) {
        file << node.x << " " << node.y << " " << node.z << "\n";
    }
    file << "$EndNodes\n";
    
    // Write elements (simplified - single block per type)
    std::map<GmshElementType, std::vector<const MeshElement*>> elements_by_type;
    for (const auto& elem : mesh_->elements) {
        elements_by_type[elem.type].push_back(&elem);
    }
    
    file << "$Elements\n";
    file << elements_by_type.size() << " " << mesh_->elements.size() 
         << " 1 " << mesh_->elements.size() << "\n";
    
    for (const auto& [type, elems] : elements_by_type) {
        int type_id = static_cast<int>(type);
        file << mesh_->dimension << " 1 " << type_id << " " << elems.size() << "\n";
        
        for (const auto* elem : elems) {
            file << elem->tag;
            for (int node : elem->node_tags) {
                file << " " << node;
            }
            file << "\n";
        }
    }
    file << "$EndElements\n";
    
    file.close();
    return true;
}

bool GmshIO::writeVtkFile(const std::string& filename) const {
    if (!mesh_ || mesh_->nodes.empty()) {
        std::cerr << "GmshIO: No mesh to write" << std::endl;
        return false;
    }
    
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "GmshIO: Cannot create file: " << filename << std::endl;
        return false;
    }
    
    // Create node tag to index mapping
    std::map<int, int> tag_to_index;
    for (size_t i = 0; i < mesh_->nodes.size(); i++) {
        tag_to_index[mesh_->nodes[i].tag] = static_cast<int>(i);
    }
    
    // VTK XML format header
    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
    file << "  <UnstructuredGrid>\n";
    file << "    <Piece NumberOfPoints=\"" << mesh_->nodes.size() 
         << "\" NumberOfCells=\"" << mesh_->elements.size() << "\">\n";
    
    // Points
    file << "      <Points>\n";
    file << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto& node : mesh_->nodes) {
        file << "          " << node.x << " " << node.y << " " << node.z << "\n";
    }
    file << "        </DataArray>\n";
    file << "      </Points>\n";
    
    // Cells
    file << "      <Cells>\n";
    
    // Connectivity
    file << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (const auto& elem : mesh_->elements) {
        file << "         ";
        for (int tag : elem.node_tags) {
            file << " " << tag_to_index[tag];
        }
        file << "\n";
    }
    file << "        </DataArray>\n";
    
    // Offsets
    file << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    int offset = 0;
    for (const auto& elem : mesh_->elements) {
        offset += elem.node_tags.size();
        file << "          " << offset << "\n";
    }
    file << "        </DataArray>\n";
    
    // Cell types (VTK numbering)
    file << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (const auto& elem : mesh_->elements) {
        int vtk_type;
        switch (elem.type) {
            case GmshElementType::POINT: vtk_type = 1; break;
            case GmshElementType::LINE: vtk_type = 3; break;
            case GmshElementType::TRIANGLE: vtk_type = 5; break;
            case GmshElementType::QUAD: vtk_type = 9; break;
            case GmshElementType::TETRAHEDRON: vtk_type = 10; break;
            case GmshElementType::HEXAHEDRON: vtk_type = 12; break;
            case GmshElementType::PRISM: vtk_type = 13; break;
            case GmshElementType::PYRAMID: vtk_type = 14; break;
            case GmshElementType::LINE3: vtk_type = 21; break;
            case GmshElementType::TRIANGLE6: vtk_type = 22; break;
            case GmshElementType::TET10: vtk_type = 24; break;
            case GmshElementType::HEX20: vtk_type = 25; break;
            default: vtk_type = 0;
        }
        file << "          " << vtk_type << "\n";
    }
    file << "        </DataArray>\n";
    file << "      </Cells>\n";
    
    // Cell data (physical tags)
    file << "      <CellData>\n";
    file << "        <DataArray type=\"Int32\" Name=\"PhysicalGroup\" format=\"ascii\">\n";
    for (const auto& elem : mesh_->elements) {
        file << "          " << elem.physical_tag << "\n";
    }
    file << "        </DataArray>\n";
    file << "      </CellData>\n";
    
    file << "    </Piece>\n";
    file << "  </UnstructuredGrid>\n";
    file << "</VTKFile>\n";
    
    file.close();
    return true;
}

PetscErrorCode GmshIO::createDMPlex(MPI_Comm comm, DM* dm) const {
    PetscFunctionBeginUser;
    
    if (!mesh_ || mesh_->nodes.empty()) {
        SETERRQ(comm, PETSC_ERR_ARG_NULL, "No mesh loaded");
    }
    
    PetscErrorCode ierr;
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    // Count cells (highest dimension elements)
    int num_cells = 0;
    int num_vertices = mesh_->nodes.size();
    
    std::vector<int> cell_indices;
    for (size_t i = 0; i < mesh_->elements.size(); i++) {
        if (mesh_->elements[i].getDimension() == mesh_->dimension) {
            cell_indices.push_back(static_cast<int>(i));
            num_cells++;
        }
    }
    
    // Create node tag to index mapping
    std::map<int, int> tag_to_index;
    for (size_t i = 0; i < mesh_->nodes.size(); i++) {
        tag_to_index[mesh_->nodes[i].tag] = static_cast<int>(i);
    }
    
    // Build cell-vertex connectivity
    std::vector<int> cell_types;
    std::vector<int> cell_vertices;
    
    for (int idx : cell_indices) {
        const auto& elem = mesh_->elements[idx];
        
        // Determine DMPlex cell type
        int plex_type;
        switch (elem.type) {
            case GmshElementType::TRIANGLE: plex_type = DM_POLYTOPE_TRIANGLE; break;
            case GmshElementType::QUAD: plex_type = DM_POLYTOPE_QUADRILATERAL; break;
            case GmshElementType::TETRAHEDRON: plex_type = DM_POLYTOPE_TETRAHEDRON; break;
            case GmshElementType::HEXAHEDRON: plex_type = DM_POLYTOPE_HEXAHEDRON; break;
            case GmshElementType::PRISM: plex_type = DM_POLYTOPE_TRI_PRISM; break;
            case GmshElementType::PYRAMID: plex_type = DM_POLYTOPE_PYRAMID; break;
            default: plex_type = DM_POLYTOPE_UNKNOWN;
        }
        cell_types.push_back(plex_type);
        
        for (int node_tag : elem.node_tags) {
            cell_vertices.push_back(tag_to_index[node_tag]);
        }
    }
    
    // Create DMPlex from cells and vertices on rank 0
    if (rank == 0) {
        // Build cone sizes array
        std::vector<int> cone_sizes(num_cells);
        for (size_t i = 0; i < cell_indices.size(); i++) {
            cone_sizes[i] = mesh_->elements[cell_indices[i]].node_tags.size();
        }
        
        ierr = DMPlexCreateFromCellListPetsc(comm, mesh_->dimension, 
                                            num_cells, num_vertices,
                                            cell_vertices.size() / num_cells,
                                            PETSC_TRUE,
                                            cell_vertices.data(),
                                            mesh_->dimension,
                                            nullptr, dm); CHKERRQ(ierr);
    } else {
        ierr = DMPlexCreateFromCellListPetsc(comm, mesh_->dimension,
                                            0, 0, 0, PETSC_TRUE,
                                            nullptr, mesh_->dimension,
                                            nullptr, dm); CHKERRQ(ierr);
    }
    
    // Distribute mesh
    DM dm_dist = nullptr;
    ierr = DMPlexDistribute(*dm, 0, nullptr, &dm_dist); CHKERRQ(ierr);
    if (dm_dist) {
        ierr = DMDestroy(dm); CHKERRQ(ierr);
        *dm = dm_dist;
    }
    
    // Set coordinates
    ierr = DMSetFromOptions(*dm); CHKERRQ(ierr);
    
    // Create coordinate section
    PetscSection coord_section;
    ierr = DMGetCoordinateSection(*dm, &coord_section); CHKERRQ(ierr);
    ierr = PetscSectionSetNumFields(coord_section, 1); CHKERRQ(ierr);
    ierr = PetscSectionSetFieldComponents(coord_section, 0, mesh_->dimension); CHKERRQ(ierr);
    
    PetscInt vStart, vEnd;
    ierr = DMPlexGetDepthStratum(*dm, 0, &vStart, &vEnd); CHKERRQ(ierr);
    
    for (PetscInt v = vStart; v < vEnd; v++) {
        ierr = PetscSectionSetDof(coord_section, v, mesh_->dimension); CHKERRQ(ierr);
        ierr = PetscSectionSetFieldDof(coord_section, v, 0, mesh_->dimension); CHKERRQ(ierr);
    }
    ierr = PetscSectionSetUp(coord_section); CHKERRQ(ierr);
    
    // Create coordinate vector
    Vec coordinates;
    ierr = DMGetCoordinatesLocal(*dm, &coordinates); CHKERRQ(ierr);
    if (!coordinates) {
        ierr = DMCreateLocalVector(*dm, &coordinates); CHKERRQ(ierr);
    }
    
    PetscScalar* coords;
    ierr = VecGetArray(coordinates, &coords); CHKERRQ(ierr);
    
    for (PetscInt v = vStart; v < vEnd; v++) {
        PetscInt idx = v - vStart;
        if (idx < static_cast<PetscInt>(mesh_->nodes.size())) {
            PetscInt offset;
            ierr = PetscSectionGetOffset(coord_section, v, &offset); CHKERRQ(ierr);
            coords[offset] = mesh_->nodes[idx].x;
            coords[offset + 1] = mesh_->nodes[idx].y;
            if (mesh_->dimension == 3) {
                coords[offset + 2] = mesh_->nodes[idx].z;
            }
        }
    }
    
    ierr = VecRestoreArray(coordinates, &coords); CHKERRQ(ierr);
    ierr = DMSetCoordinatesLocal(*dm, coordinates); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

PetscErrorCode GmshIO::createDMPlexWithBoundaries(
    MPI_Comm comm, 
    DM* dm,
    const std::vector<std::string>& boundary_groups
) const {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    
    // First create the basic DMPlex
    ierr = createDMPlex(comm, dm); CHKERRQ(ierr);
    
    // Label boundaries based on physical groups
    DMLabel label;
    ierr = DMCreateLabel(*dm, "boundary"); CHKERRQ(ierr);
    ierr = DMGetLabel(*dm, "boundary", &label); CHKERRQ(ierr);
    
    for (const auto& group_name : boundary_groups) {
        int tag = mesh_->getPhysicalTagByName(group_name);
        if (tag >= 0) {
            auto element_indices = mesh_->getElementsByPhysicalTag(tag);
            for (int idx : element_indices) {
                // Mark faces in DMPlex corresponding to boundary elements
                ierr = DMLabelSetValue(label, idx, tag); CHKERRQ(ierr);
            }
        }
    }
    
    PetscFunctionReturn(0);
}

void GmshIO::printStatistics() const {
    if (!mesh_) {
        std::cout << "No mesh loaded\n";
        return;
    }
    
    std::cout << "Mesh Statistics:\n";
    std::cout << "  Format version: " << format_version_ << "\n";
    std::cout << "  Dimension: " << mesh_->dimension << "D\n";
    std::cout << "  Number of nodes: " << mesh_->getNumNodes() << "\n";
    std::cout << "  Number of elements: " << mesh_->getNumElements() << "\n";
    std::cout << "  Number of cells: " << mesh_->getNumCells() << "\n";
    std::cout << "  Number of faces: " << mesh_->getNumFaces() << "\n";
    
    std::cout << "  Bounding box:\n";
    std::cout << "    X: [" << mesh_->bbox_min[0] << ", " << mesh_->bbox_max[0] << "]\n";
    std::cout << "    Y: [" << mesh_->bbox_min[1] << ", " << mesh_->bbox_max[1] << "]\n";
    std::cout << "    Z: [" << mesh_->bbox_min[2] << ", " << mesh_->bbox_max[2] << "]\n";
    
    if (!mesh_->physical_groups.empty()) {
        std::cout << "  Physical groups:\n";
        for (const auto& [tag, group] : mesh_->physical_groups) {
            std::cout << "    " << group.name << " (dim=" << group.dimension 
                     << ", tag=" << tag << ")\n";
        }
    }
    
    // Count element types
    std::map<GmshElementType, int> type_counts;
    for (const auto& elem : mesh_->elements) {
        type_counts[elem.type]++;
    }
    
    std::cout << "  Element types:\n";
    for (const auto& [type, count] : type_counts) {
        std::string name;
        switch (type) {
            case GmshElementType::LINE: name = "Line"; break;
            case GmshElementType::TRIANGLE: name = "Triangle"; break;
            case GmshElementType::QUAD: name = "Quad"; break;
            case GmshElementType::TETRAHEDRON: name = "Tetrahedron"; break;
            case GmshElementType::HEXAHEDRON: name = "Hexahedron"; break;
            case GmshElementType::PRISM: name = "Prism"; break;
            case GmshElementType::PYRAMID: name = "Pyramid"; break;
            default: name = "Other";
        }
        std::cout << "    " << name << ": " << count << "\n";
    }
}

bool GmshIO::validate() const {
    if (!mesh_) return false;
    if (mesh_->nodes.empty()) return false;
    if (mesh_->elements.empty()) return false;
    
    // Check that all element node tags reference valid nodes
    std::set<int> node_tags;
    for (const auto& node : mesh_->nodes) {
        node_tags.insert(node.tag);
    }
    
    for (const auto& elem : mesh_->elements) {
        for (int tag : elem.node_tags) {
            if (node_tags.find(tag) == node_tags.end()) {
                std::cerr << "GmshIO: Element " << elem.tag 
                         << " references invalid node " << tag << "\n";
                return false;
            }
        }
    }
    
    return true;
}

} // namespace FSRM
