#include "GmshIO.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cstring>
#include <limits>
#include <set>
#include <cmath>
#include <numeric>

namespace FSRM {

// ============================================================================
// UnstructuredMesh - Material Domain Methods
// ============================================================================

MaterialDomain* UnstructuredMesh::getMaterialDomain(const std::string& name) {
    auto it = material_domains.find(name);
    return (it != material_domains.end()) ? &it->second : nullptr;
}

const MaterialDomain* UnstructuredMesh::getMaterialDomain(const std::string& name) const {
    auto it = material_domains.find(name);
    return (it != material_domains.end()) ? &it->second : nullptr;
}

std::vector<std::string> UnstructuredMesh::getMaterialDomainNames() const {
    std::vector<std::string> names;
    for (const auto& [name, domain] : material_domains) {
        names.push_back(name);
    }
    return names;
}

int UnstructuredMesh::getMaterialIdForElement(int element_idx) const {
    if (element_idx >= 0 && element_idx < static_cast<int>(element_material_ids.size())) {
        return element_material_ids[element_idx];
    }
    return -1;  // Unknown material
}

// ============================================================================
// UnstructuredMesh - Fault Surface Methods
// ============================================================================

FaultSurface* UnstructuredMesh::getFaultSurface(const std::string& name) {
    auto it = fault_surfaces.find(name);
    return (it != fault_surfaces.end()) ? &it->second : nullptr;
}

const FaultSurface* UnstructuredMesh::getFaultSurface(const std::string& name) const {
    auto it = fault_surfaces.find(name);
    return (it != fault_surfaces.end()) ? &it->second : nullptr;
}

std::vector<std::string> UnstructuredMesh::getFaultNames() const {
    std::vector<std::string> names;
    for (const auto& [name, fault] : fault_surfaces) {
        names.push_back(name);
    }
    return names;
}

// ============================================================================
// UnstructuredMesh - Boundary Surface Methods
// ============================================================================

BoundarySurface* UnstructuredMesh::getBoundarySurface(const std::string& name) {
    auto it = boundary_surfaces.find(name);
    return (it != boundary_surfaces.end()) ? &it->second : nullptr;
}

const BoundarySurface* UnstructuredMesh::getBoundarySurface(const std::string& name) const {
    auto it = boundary_surfaces.find(name);
    return (it != boundary_surfaces.end()) ? &it->second : nullptr;
}

std::vector<std::string> UnstructuredMesh::getBoundaryNames() const {
    std::vector<std::string> names;
    for (const auto& [name, boundary] : boundary_surfaces) {
        names.push_back(name);
    }
    return names;
}

std::vector<int> UnstructuredMesh::getBoundaryNodes(const std::string& name) const {
    auto it = boundary_surfaces.find(name);
    if (it != boundary_surfaces.end()) {
        return it->second.node_indices;
    }
    return {};
}

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

// ============================================================================
// Domain Mapping Configuration
// ============================================================================

void GmshIO::setMaterialMapping(const std::vector<GmshMaterialMapping>& mappings) {
    material_mappings_ = mappings;
    // Auto-assign material IDs if not specified
    for (size_t i = 0; i < material_mappings_.size(); i++) {
        if (material_mappings_[i].material_id < 0) {
            material_mappings_[i].material_id = static_cast<int>(i);
        }
    }
}

void GmshIO::addMaterialMapping(const std::string& physical_group,
                                const std::string& material_section,
                                int material_id) {
    int id = (material_id >= 0) ? material_id : static_cast<int>(material_mappings_.size());
    material_mappings_.emplace_back(physical_group, material_section, id);
}

void GmshIO::setFaultMapping(const std::vector<GmshFaultMapping>& mappings) {
    fault_mappings_ = mappings;
}

void GmshIO::addFaultMapping(const std::string& physical_group,
                             const std::string& fault_section,
                             bool use_split_nodes) {
    fault_mappings_.emplace_back(physical_group, fault_section, use_split_nodes);
}

void GmshIO::setBoundaryNames(const std::vector<std::string>& boundary_names) {
    boundary_names_ = boundary_names;
}

bool GmshIO::processDomains() {
    if (!mesh_) {
        std::cerr << "GmshIO: No mesh loaded\n";
        return false;
    }
    
    bool success = true;
    
    // Process material domains
    if (!material_mappings_.empty()) {
        success &= processMaterialDomains();
    }
    
    // Process fault surfaces
    if (!fault_mappings_.empty()) {
        success &= processFaultSurfaces();
    }
    
    // Process boundary surfaces
    if (!boundary_names_.empty()) {
        success &= processBoundarySurfaces();
    }
    
    return success;
}

bool GmshIO::processMaterialDomains() {
    if (!mesh_) return false;
    
    mesh_->material_domains.clear();
    mesh_->element_material_ids.resize(mesh_->elements.size(), -1);
    
    for (const auto& mapping : material_mappings_) {
        // Find the physical group tag
        int tag = mesh_->getPhysicalTagByName(mapping.physical_group_name);
        if (tag < 0) {
            std::cerr << "GmshIO: Physical group '" << mapping.physical_group_name 
                     << "' not found in mesh\n";
            continue;
        }
        
        // Create material domain
        MaterialDomain domain;
        domain.physical_tag = tag;
        domain.name = mapping.physical_group_name;
        domain.material_name = mapping.material_section;
        domain.material_id = mapping.material_id;
        
        // Get dimension from physical group
        auto pg_it = mesh_->physical_groups.find(tag);
        if (pg_it != mesh_->physical_groups.end()) {
            domain.dimension = pg_it->second.dimension;
        }
        
        // Find elements belonging to this domain
        domain.num_elements = 0;
        for (size_t i = 0; i < mesh_->elements.size(); i++) {
            if (mesh_->elements[i].physical_tag == tag) {
                mesh_->element_material_ids[i] = domain.material_id;
                domain.num_elements++;
            }
        }
        
        // Compute domain statistics
        computeDomainStatistics(domain);
        
        mesh_->material_domains[domain.name] = domain;
        
        std::cout << "GmshIO: Material domain '" << domain.name 
                 << "' -> " << domain.material_name
                 << " (" << domain.num_elements << " elements, volume=" 
                 << domain.volume << ")\n";
    }
    
    return true;
}

bool GmshIO::processFaultSurfaces() {
    if (!mesh_) return false;
    
    mesh_->fault_surfaces.clear();
    
    for (const auto& mapping : fault_mappings_) {
        // Find the physical group tag
        int tag = mesh_->getPhysicalTagByName(mapping.physical_group_name);
        if (tag < 0) {
            std::cerr << "GmshIO: Fault physical group '" << mapping.physical_group_name 
                     << "' not found in mesh\n";
            continue;
        }
        
        // Create fault surface
        FaultSurface fault;
        fault.physical_tag = tag;
        fault.name = mapping.physical_group_name;
        fault.requires_split_nodes = mapping.use_split_nodes;
        
        // Get dimension from physical group
        auto pg_it = mesh_->physical_groups.find(tag);
        if (pg_it != mesh_->physical_groups.end()) {
            fault.dimension = pg_it->second.dimension;
        }
        
        // Find elements and nodes belonging to this fault
        std::set<int> unique_nodes;
        for (size_t i = 0; i < mesh_->elements.size(); i++) {
            if (mesh_->elements[i].physical_tag == tag) {
                fault.element_indices.push_back(static_cast<int>(i));
                for (int node_tag : mesh_->elements[i].node_tags) {
                    unique_nodes.insert(node_tag);
                }
            }
        }
        
        fault.node_indices.assign(unique_nodes.begin(), unique_nodes.end());
        
        // Compute fault geometry
        computeFaultGeometry(fault);
        
        mesh_->fault_surfaces[fault.name] = fault;
        
        std::cout << "GmshIO: Fault surface '" << fault.name 
                 << "' (" << fault.element_indices.size() << " elements, "
                 << fault.node_indices.size() << " nodes, area=" 
                 << fault.area << " m²)\n";
    }
    
    return true;
}

bool GmshIO::processBoundarySurfaces() {
    if (!mesh_) return false;
    
    mesh_->boundary_surfaces.clear();
    
    for (const auto& name : boundary_names_) {
        // Find the physical group tag
        int tag = mesh_->getPhysicalTagByName(name);
        if (tag < 0) {
            std::cerr << "GmshIO: Boundary physical group '" << name 
                     << "' not found in mesh\n";
            continue;
        }
        
        // Create boundary surface
        BoundarySurface boundary;
        boundary.physical_tag = tag;
        boundary.name = name;
        
        // Get dimension from physical group
        auto pg_it = mesh_->physical_groups.find(tag);
        if (pg_it != mesh_->physical_groups.end()) {
            boundary.dimension = pg_it->second.dimension;
        }
        
        // Find elements and nodes belonging to this boundary
        std::set<int> unique_nodes;
        boundary.num_elements = 0;
        for (size_t i = 0; i < mesh_->elements.size(); i++) {
            if (mesh_->elements[i].physical_tag == tag) {
                boundary.element_indices.push_back(static_cast<int>(i));
                boundary.num_elements++;
                for (int node_tag : mesh_->elements[i].node_tags) {
                    unique_nodes.insert(node_tag);
                }
            }
        }
        
        boundary.node_indices.assign(unique_nodes.begin(), unique_nodes.end());
        
        // Compute boundary statistics
        computeBoundaryStatistics(boundary);
        
        mesh_->boundary_surfaces[boundary.name] = boundary;
        
        std::cout << "GmshIO: Boundary surface '" << boundary.name 
                 << "' (" << boundary.num_elements << " elements, "
                 << boundary.node_indices.size() << " nodes, area=" 
                 << boundary.area << ")\n";
    }
    
    return true;
}

void GmshIO::computeDomainStatistics(MaterialDomain& domain) {
    if (!mesh_) return;
    
    domain.volume = 0.0;
    domain.centroid = {0.0, 0.0, 0.0};
    
    double total_weight = 0.0;
    
    for (size_t i = 0; i < mesh_->elements.size(); i++) {
        if (mesh_->elements[i].physical_tag == domain.physical_tag) {
            double vol = computeElementVolume(static_cast<int>(i));
            domain.volume += vol;
            
            // Compute element centroid
            const auto& elem = mesh_->elements[i];
            double cx = 0, cy = 0, cz = 0;
            for (int node_tag : elem.node_tags) {
                auto coords = getNodeCoords(node_tag);
                cx += coords[0];
                cy += coords[1];
                cz += coords[2];
            }
            int n = static_cast<int>(elem.node_tags.size());
            cx /= n; cy /= n; cz /= n;
            
            domain.centroid[0] += cx * vol;
            domain.centroid[1] += cy * vol;
            domain.centroid[2] += cz * vol;
            total_weight += vol;
        }
    }
    
    if (total_weight > 0) {
        domain.centroid[0] /= total_weight;
        domain.centroid[1] /= total_weight;
        domain.centroid[2] /= total_weight;
    }
}

void GmshIO::computeFaultGeometry(FaultSurface& fault) {
    if (!mesh_ || fault.element_indices.empty()) return;
    
    fault.area = 0.0;
    fault.centroid = {0.0, 0.0, 0.0};
    fault.normal = {0.0, 0.0, 0.0};
    
    double total_area = 0.0;
    
    for (int elem_idx : fault.element_indices) {
        double area = computeElementVolume(elem_idx);
        fault.area += area;
        
        // Compute element centroid
        const auto& elem = mesh_->elements[elem_idx];
        double cx = 0, cy = 0, cz = 0;
        for (int node_tag : elem.node_tags) {
            auto coords = getNodeCoords(node_tag);
            cx += coords[0];
            cy += coords[1];
            cz += coords[2];
        }
        int n = static_cast<int>(elem.node_tags.size());
        cx /= n; cy /= n; cz /= n;
        
        fault.centroid[0] += cx * area;
        fault.centroid[1] += cy * area;
        fault.centroid[2] += cz * area;
        
        // Accumulate normal (weighted by area)
        auto normal = computeElementNormal(elem_idx);
        fault.normal[0] += normal[0] * area;
        fault.normal[1] += normal[1] * area;
        fault.normal[2] += normal[2] * area;
        
        total_area += area;
    }
    
    if (total_area > 0) {
        fault.centroid[0] /= total_area;
        fault.centroid[1] /= total_area;
        fault.centroid[2] /= total_area;
        
        // Normalize the normal vector
        double norm = std::sqrt(fault.normal[0]*fault.normal[0] + 
                                fault.normal[1]*fault.normal[1] + 
                                fault.normal[2]*fault.normal[2]);
        if (norm > 1e-12) {
            fault.normal[0] /= norm;
            fault.normal[1] /= norm;
            fault.normal[2] /= norm;
        }
    }
    
    // Compute approximate strike and dip from normal
    // Strike is perpendicular to horizontal projection of normal
    // Dip is angle from horizontal
    double nx = fault.normal[0];
    double ny = fault.normal[1];
    double nz = fault.normal[2];
    
    double horiz_len = std::sqrt(nx*nx + ny*ny);
    if (horiz_len > 1e-12) {
        fault.strike = std::atan2(ny, nx);  // Strike direction
        fault.dip = std::atan2(horiz_len, std::abs(nz));  // Dip angle
    } else {
        fault.strike = 0.0;
        fault.dip = 0.0;  // Horizontal fault
    }
    
    // Estimate length and width from bounding box of fault nodes
    double min_x = std::numeric_limits<double>::max();
    double max_x = std::numeric_limits<double>::lowest();
    double min_y = std::numeric_limits<double>::max();
    double max_y = std::numeric_limits<double>::lowest();
    double min_z = std::numeric_limits<double>::max();
    double max_z = std::numeric_limits<double>::lowest();
    
    for (int node_idx : fault.node_indices) {
        auto coords = getNodeCoords(node_idx);
        min_x = std::min(min_x, coords[0]);
        max_x = std::max(max_x, coords[0]);
        min_y = std::min(min_y, coords[1]);
        max_y = std::max(max_y, coords[1]);
        min_z = std::min(min_z, coords[2]);
        max_z = std::max(max_z, coords[2]);
    }
    
    double dx = max_x - min_x;
    double dy = max_y - min_y;
    double dz = max_z - min_z;
    
    fault.length = std::sqrt(dx*dx + dy*dy);  // Horizontal extent
    fault.width = std::sqrt(dz*dz + std::min(dx, dy)*std::min(dx, dy));  // Vertical extent
}

void GmshIO::computeBoundaryStatistics(BoundarySurface& boundary) {
    if (!mesh_ || boundary.element_indices.empty()) return;
    
    boundary.area = 0.0;
    boundary.centroid = {0.0, 0.0, 0.0};
    boundary.normal = {0.0, 0.0, 0.0};
    
    double total_area = 0.0;
    
    for (int elem_idx : boundary.element_indices) {
        double area = computeElementVolume(elem_idx);
        boundary.area += area;
        
        // Compute element centroid
        const auto& elem = mesh_->elements[elem_idx];
        double cx = 0, cy = 0, cz = 0;
        for (int node_tag : elem.node_tags) {
            auto coords = getNodeCoords(node_tag);
            cx += coords[0];
            cy += coords[1];
            cz += coords[2];
        }
        int n = static_cast<int>(elem.node_tags.size());
        cx /= n; cy /= n; cz /= n;
        
        boundary.centroid[0] += cx * area;
        boundary.centroid[1] += cy * area;
        boundary.centroid[2] += cz * area;
        
        // Accumulate normal (weighted by area)
        auto normal = computeElementNormal(elem_idx);
        boundary.normal[0] += normal[0] * area;
        boundary.normal[1] += normal[1] * area;
        boundary.normal[2] += normal[2] * area;
        
        total_area += area;
    }
    
    if (total_area > 0) {
        boundary.centroid[0] /= total_area;
        boundary.centroid[1] /= total_area;
        boundary.centroid[2] /= total_area;
        
        // Normalize the normal vector
        double norm = std::sqrt(boundary.normal[0]*boundary.normal[0] + 
                                boundary.normal[1]*boundary.normal[1] + 
                                boundary.normal[2]*boundary.normal[2]);
        if (norm > 1e-12) {
            boundary.normal[0] /= norm;
            boundary.normal[1] /= norm;
            boundary.normal[2] /= norm;
        }
    }
}

// ============================================================================
// Mesh Access Helpers
// ============================================================================

std::vector<std::string> GmshIO::getPhysicalGroupNames() const {
    std::vector<std::string> names;
    if (mesh_) {
        for (const auto& [tag, group] : mesh_->physical_groups) {
            names.push_back(group.name);
        }
    }
    return names;
}

std::vector<std::string> GmshIO::getPhysicalGroupsByDimension(int dim) const {
    std::vector<std::string> names;
    if (mesh_) {
        for (const auto& [tag, group] : mesh_->physical_groups) {
            if (group.dimension == dim) {
                names.push_back(group.name);
            }
        }
    }
    return names;
}

// ============================================================================
// Geometry Computation Helpers
// ============================================================================

std::array<double, 3> GmshIO::getNodeCoords(int node_tag) const {
    if (!mesh_) return {0, 0, 0};
    
    for (const auto& node : mesh_->nodes) {
        if (node.tag == node_tag) {
            return {node.x, node.y, node.z};
        }
    }
    return {0, 0, 0};
}

int GmshIO::findNodeIndex(int node_tag) const {
    if (!mesh_) return -1;
    
    for (size_t i = 0; i < mesh_->nodes.size(); i++) {
        if (mesh_->nodes[i].tag == node_tag) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

double GmshIO::computeElementVolume(int element_idx) const {
    if (!mesh_ || element_idx < 0 || element_idx >= static_cast<int>(mesh_->elements.size())) {
        return 0.0;
    }
    
    const auto& elem = mesh_->elements[element_idx];
    int dim = elem.getDimension();
    
    if (dim == 3) {
        // Volume element (tetrahedron, hexahedron, etc.)
        if (elem.type == GmshElementType::TETRAHEDRON || elem.type == GmshElementType::TET10) {
            // Tetrahedron volume: |((v1-v0) · ((v2-v0) × (v3-v0)))| / 6
            auto p0 = getNodeCoords(elem.node_tags[0]);
            auto p1 = getNodeCoords(elem.node_tags[1]);
            auto p2 = getNodeCoords(elem.node_tags[2]);
            auto p3 = getNodeCoords(elem.node_tags[3]);
            
            double v10[3] = {p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2]};
            double v20[3] = {p2[0]-p0[0], p2[1]-p0[1], p2[2]-p0[2]};
            double v30[3] = {p3[0]-p0[0], p3[1]-p0[1], p3[2]-p0[2]};
            
            // Cross product v20 × v30
            double cross[3] = {
                v20[1]*v30[2] - v20[2]*v30[1],
                v20[2]*v30[0] - v20[0]*v30[2],
                v20[0]*v30[1] - v20[1]*v30[0]
            };
            
            // Dot product v10 · cross
            double vol = std::abs(v10[0]*cross[0] + v10[1]*cross[1] + v10[2]*cross[2]) / 6.0;
            return vol;
        }
        // For other 3D elements, use approximate volume based on bounding box
        // (proper implementation would require more complex formulas)
        double min_x = std::numeric_limits<double>::max();
        double max_x = std::numeric_limits<double>::lowest();
        double min_y = std::numeric_limits<double>::max();
        double max_y = std::numeric_limits<double>::lowest();
        double min_z = std::numeric_limits<double>::max();
        double max_z = std::numeric_limits<double>::lowest();
        
        for (int node_tag : elem.node_tags) {
            auto coords = getNodeCoords(node_tag);
            min_x = std::min(min_x, coords[0]);
            max_x = std::max(max_x, coords[0]);
            min_y = std::min(min_y, coords[1]);
            max_y = std::max(max_y, coords[1]);
            min_z = std::min(min_z, coords[2]);
            max_z = std::max(max_z, coords[2]);
        }
        return (max_x - min_x) * (max_y - min_y) * (max_z - min_z);
    }
    else if (dim == 2) {
        // Surface element (triangle, quad, etc.)
        if (elem.type == GmshElementType::TRIANGLE || elem.type == GmshElementType::TRIANGLE6) {
            // Triangle area: |((v1-v0) × (v2-v0))| / 2
            auto p0 = getNodeCoords(elem.node_tags[0]);
            auto p1 = getNodeCoords(elem.node_tags[1]);
            auto p2 = getNodeCoords(elem.node_tags[2]);
            
            double v10[3] = {p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2]};
            double v20[3] = {p2[0]-p0[0], p2[1]-p0[1], p2[2]-p0[2]};
            
            // Cross product
            double cross[3] = {
                v10[1]*v20[2] - v10[2]*v20[1],
                v10[2]*v20[0] - v10[0]*v20[2],
                v10[0]*v20[1] - v10[1]*v20[0]
            };
            
            return std::sqrt(cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2]) / 2.0;
        }
        else if (elem.type == GmshElementType::QUAD || elem.type == GmshElementType::QUAD8 || 
                 elem.type == GmshElementType::QUAD9) {
            // Quad area: sum of two triangles
            auto p0 = getNodeCoords(elem.node_tags[0]);
            auto p1 = getNodeCoords(elem.node_tags[1]);
            auto p2 = getNodeCoords(elem.node_tags[2]);
            auto p3 = getNodeCoords(elem.node_tags[3]);
            
            // Triangle 0-1-2
            double v10[3] = {p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2]};
            double v20[3] = {p2[0]-p0[0], p2[1]-p0[1], p2[2]-p0[2]};
            double cross1[3] = {
                v10[1]*v20[2] - v10[2]*v20[1],
                v10[2]*v20[0] - v10[0]*v20[2],
                v10[0]*v20[1] - v10[1]*v20[0]
            };
            double area1 = std::sqrt(cross1[0]*cross1[0] + cross1[1]*cross1[1] + cross1[2]*cross1[2]) / 2.0;
            
            // Triangle 0-2-3
            double v30[3] = {p3[0]-p0[0], p3[1]-p0[1], p3[2]-p0[2]};
            double cross2[3] = {
                v20[1]*v30[2] - v20[2]*v30[1],
                v20[2]*v30[0] - v20[0]*v30[2],
                v20[0]*v30[1] - v20[1]*v30[0]
            };
            double area2 = std::sqrt(cross2[0]*cross2[0] + cross2[1]*cross2[1] + cross2[2]*cross2[2]) / 2.0;
            
            return area1 + area2;
        }
    }
    else if (dim == 1) {
        // Line element
        if (elem.node_tags.size() >= 2) {
            auto p0 = getNodeCoords(elem.node_tags[0]);
            auto p1 = getNodeCoords(elem.node_tags[1]);
            return std::sqrt(
                (p1[0]-p0[0])*(p1[0]-p0[0]) + 
                (p1[1]-p0[1])*(p1[1]-p0[1]) + 
                (p1[2]-p0[2])*(p1[2]-p0[2])
            );
        }
    }
    
    return 0.0;
}

std::array<double, 3> GmshIO::computeElementNormal(int element_idx) const {
    if (!mesh_ || element_idx < 0 || element_idx >= static_cast<int>(mesh_->elements.size())) {
        return {0, 0, 1};
    }
    
    const auto& elem = mesh_->elements[element_idx];
    int dim = elem.getDimension();
    
    if (dim == 2 && elem.node_tags.size() >= 3) {
        // Surface element - compute normal from first three nodes
        auto p0 = getNodeCoords(elem.node_tags[0]);
        auto p1 = getNodeCoords(elem.node_tags[1]);
        auto p2 = getNodeCoords(elem.node_tags[2]);
        
        double v10[3] = {p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2]};
        double v20[3] = {p2[0]-p0[0], p2[1]-p0[1], p2[2]-p0[2]};
        
        // Cross product
        double cross[3] = {
            v10[1]*v20[2] - v10[2]*v20[1],
            v10[2]*v20[0] - v10[0]*v20[2],
            v10[0]*v20[1] - v10[1]*v20[0]
        };
        
        double norm = std::sqrt(cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2]);
        if (norm > 1e-12) {
            return {cross[0]/norm, cross[1]/norm, cross[2]/norm};
        }
    }
    
    return {0, 0, 1};  // Default to z-up
}

double GmshIO::computeElementQuality(int element_idx) const {
    if (!mesh_ || element_idx < 0 || element_idx >= static_cast<int>(mesh_->elements.size())) {
        return 0.0;
    }
    
    const auto& elem = mesh_->elements[element_idx];
    
    // For triangles: quality = 4*sqrt(3)*area / (sum of squared edge lengths)
    // For tetrahedra: quality = 216*sqrt(3)*volume / (sum of squared edge lengths)^(3/2)
    
    if (elem.type == GmshElementType::TRIANGLE || elem.type == GmshElementType::TRIANGLE6) {
        auto p0 = getNodeCoords(elem.node_tags[0]);
        auto p1 = getNodeCoords(elem.node_tags[1]);
        auto p2 = getNodeCoords(elem.node_tags[2]);
        
        // Edge lengths squared
        double e01_2 = (p1[0]-p0[0])*(p1[0]-p0[0]) + (p1[1]-p0[1])*(p1[1]-p0[1]) + (p1[2]-p0[2])*(p1[2]-p0[2]);
        double e12_2 = (p2[0]-p1[0])*(p2[0]-p1[0]) + (p2[1]-p1[1])*(p2[1]-p1[1]) + (p2[2]-p1[2])*(p2[2]-p1[2]);
        double e20_2 = (p0[0]-p2[0])*(p0[0]-p2[0]) + (p0[1]-p2[1])*(p0[1]-p2[1]) + (p0[2]-p2[2])*(p0[2]-p2[2]);
        
        double area = computeElementVolume(element_idx);
        double sum_e2 = e01_2 + e12_2 + e20_2;
        
        if (sum_e2 > 1e-20) {
            return 4.0 * std::sqrt(3.0) * area / sum_e2;
        }
    }
    else if (elem.type == GmshElementType::TETRAHEDRON || elem.type == GmshElementType::TET10) {
        auto p0 = getNodeCoords(elem.node_tags[0]);
        auto p1 = getNodeCoords(elem.node_tags[1]);
        auto p2 = getNodeCoords(elem.node_tags[2]);
        auto p3 = getNodeCoords(elem.node_tags[3]);
        
        // Sum of squared edge lengths
        double sum_e2 = 0;
        for (int i = 0; i < 4; i++) {
            for (int j = i+1; j < 4; j++) {
                auto pi = getNodeCoords(elem.node_tags[i]);
                auto pj = getNodeCoords(elem.node_tags[j]);
                sum_e2 += (pj[0]-pi[0])*(pj[0]-pi[0]) + (pj[1]-pi[1])*(pj[1]-pi[1]) + (pj[2]-pi[2])*(pj[2]-pi[2]);
            }
        }
        
        double volume = computeElementVolume(element_idx);
        double denom = std::pow(sum_e2, 1.5);
        
        if (denom > 1e-20) {
            return 216.0 * std::sqrt(3.0) * volume / denom;
        }
    }
    
    return 1.0;  // Default to perfect quality for unsupported element types
}

bool GmshIO::validateQuality(double min_quality, double max_aspect_ratio) const {
    if (!mesh_) return false;
    
    bool all_pass = true;
    
    for (size_t i = 0; i < mesh_->elements.size(); i++) {
        double quality = computeElementQuality(static_cast<int>(i));
        if (quality < min_quality) {
            std::cerr << "GmshIO: Element " << i << " has low quality: " << quality << "\n";
            all_pass = false;
        }
    }
    
    return all_pass;
}

void GmshIO::getQualityStats(double& min_qual, double& max_qual, double& avg_qual) const {
    if (!mesh_ || mesh_->elements.empty()) {
        min_qual = max_qual = avg_qual = 0.0;
        return;
    }
    
    min_qual = std::numeric_limits<double>::max();
    max_qual = std::numeric_limits<double>::lowest();
    double sum = 0.0;
    int count = 0;
    
    for (size_t i = 0; i < mesh_->elements.size(); i++) {
        double quality = computeElementQuality(static_cast<int>(i));
        min_qual = std::min(min_qual, quality);
        max_qual = std::max(max_qual, quality);
        sum += quality;
        count++;
    }
    
    avg_qual = (count > 0) ? sum / count : 0.0;
}

// ============================================================================
// VTK Output with Materials
// ============================================================================

bool GmshIO::writeVtkFileWithMaterials(const std::string& filename) const {
    if (!mesh_ || mesh_->nodes.empty()) {
        std::cerr << "GmshIO: No mesh to write\n";
        return false;
    }
    
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "GmshIO: Cannot create file: " << filename << "\n";
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
    
    // Cell data
    file << "      <CellData>\n";
    
    // Physical group tags
    file << "        <DataArray type=\"Int32\" Name=\"PhysicalGroup\" format=\"ascii\">\n";
    for (const auto& elem : mesh_->elements) {
        file << "          " << elem.physical_tag << "\n";
    }
    file << "        </DataArray>\n";
    
    // Material IDs
    file << "        <DataArray type=\"Int32\" Name=\"MaterialID\" format=\"ascii\">\n";
    for (size_t i = 0; i < mesh_->elements.size(); i++) {
        int mat_id = -1;
        if (i < mesh_->element_material_ids.size()) {
            mat_id = mesh_->element_material_ids[i];
        }
        file << "          " << mat_id << "\n";
    }
    file << "        </DataArray>\n";
    
    file << "      </CellData>\n";
    
    file << "    </Piece>\n";
    file << "  </UnstructuredGrid>\n";
    file << "</VTKFile>\n";
    
    file.close();
    return true;
}

// ============================================================================
// DMPlex Creation with Materials
// ============================================================================

PetscErrorCode GmshIO::createDMPlexWithMaterials(MPI_Comm comm, DM* dm) const {
    PetscFunctionBeginUser;
    
    // First create basic DMPlex
    PetscErrorCode ierr = createDMPlex(comm, dm); CHKERRQ(ierr);
    
    // Add material labels
    if (mesh_ && !mesh_->material_domains.empty()) {
        DMLabel material_label;
        ierr = DMCreateLabel(*dm, "Material"); CHKERRQ(ierr);
        ierr = DMGetLabel(*dm, "Material", &material_label); CHKERRQ(ierr);
        
        for (const auto& [name, domain] : mesh_->material_domains) {
            auto elem_indices = mesh_->getElementsByPhysicalTag(domain.physical_tag);
            for (int idx : elem_indices) {
                ierr = DMLabelSetValue(material_label, idx, domain.material_id); CHKERRQ(ierr);
            }
        }
    }
    
    PetscFunctionReturn(0);
}

PetscErrorCode GmshIO::createDMPlexFull(MPI_Comm comm, DM* dm) const {
    PetscFunctionBeginUser;
    
    // First create DMPlex with materials
    PetscErrorCode ierr = createDMPlexWithMaterials(comm, dm); CHKERRQ(ierr);
    
    // Add fault labels
    if (mesh_ && !mesh_->fault_surfaces.empty()) {
        DMLabel fault_label;
        ierr = DMCreateLabel(*dm, "Fault"); CHKERRQ(ierr);
        ierr = DMGetLabel(*dm, "Fault", &fault_label); CHKERRQ(ierr);
        
        int fault_id = 0;
        for (const auto& [name, fault] : mesh_->fault_surfaces) {
            for (int idx : fault.element_indices) {
                ierr = DMLabelSetValue(fault_label, idx, fault_id); CHKERRQ(ierr);
            }
            fault_id++;
        }
    }
    
    // Add boundary labels
    if (mesh_ && !mesh_->boundary_surfaces.empty()) {
        DMLabel boundary_label;
        ierr = DMCreateLabel(*dm, "Boundary"); CHKERRQ(ierr);
        ierr = DMGetLabel(*dm, "Boundary", &boundary_label); CHKERRQ(ierr);
        
        int boundary_id = 0;
        for (const auto& [name, boundary] : mesh_->boundary_surfaces) {
            for (int idx : boundary.element_indices) {
                ierr = DMLabelSetValue(boundary_label, idx, boundary_id); CHKERRQ(ierr);
            }
            boundary_id++;
        }
    }
    
    PetscFunctionReturn(0);
}

// ============================================================================
// Detailed Output
// ============================================================================

void GmshIO::printDomainInfo() const {
    if (!mesh_) {
        std::cout << "No mesh loaded\n";
        return;
    }
    
    std::cout << "\n=== Domain Information ===\n\n";
    
    // Material domains
    if (!mesh_->material_domains.empty()) {
        std::cout << "Material Domains:\n";
        for (const auto& [name, domain] : mesh_->material_domains) {
            std::cout << "  " << name << ":\n";
            std::cout << "    Material: " << domain.material_name << " (ID=" << domain.material_id << ")\n";
            std::cout << "    Elements: " << domain.num_elements << "\n";
            std::cout << "    Volume: " << domain.volume << " m³\n";
            std::cout << "    Centroid: (" << domain.centroid[0] << ", " 
                     << domain.centroid[1] << ", " << domain.centroid[2] << ")\n";
        }
        std::cout << "\n";
    }
    
    // Fault surfaces
    if (!mesh_->fault_surfaces.empty()) {
        std::cout << "Fault Surfaces:\n";
        for (const auto& [name, fault] : mesh_->fault_surfaces) {
            std::cout << "  " << name << ":\n";
            std::cout << "    Elements: " << fault.element_indices.size() << "\n";
            std::cout << "    Nodes: " << fault.node_indices.size() << "\n";
            std::cout << "    Area: " << fault.area << " m²\n";
            std::cout << "    Length: " << fault.length << " m\n";
            std::cout << "    Width: " << fault.width << " m\n";
            std::cout << "    Strike: " << fault.strike * 180.0 / M_PI << "°\n";
            std::cout << "    Dip: " << fault.dip * 180.0 / M_PI << "°\n";
            std::cout << "    Normal: (" << fault.normal[0] << ", " 
                     << fault.normal[1] << ", " << fault.normal[2] << ")\n";
            std::cout << "    Centroid: (" << fault.centroid[0] << ", " 
                     << fault.centroid[1] << ", " << fault.centroid[2] << ")\n";
            std::cout << "    Split nodes: " << (fault.requires_split_nodes ? "yes" : "no") << "\n";
        }
        std::cout << "\n";
    }
    
    // Boundary surfaces
    if (!mesh_->boundary_surfaces.empty()) {
        std::cout << "Boundary Surfaces:\n";
        for (const auto& [name, boundary] : mesh_->boundary_surfaces) {
            std::cout << "  " << name << ":\n";
            std::cout << "    Elements: " << boundary.num_elements << "\n";
            std::cout << "    Nodes: " << boundary.node_indices.size() << "\n";
            std::cout << "    Area: " << boundary.area << " m²\n";
            std::cout << "    Normal: (" << boundary.normal[0] << ", " 
                     << boundary.normal[1] << ", " << boundary.normal[2] << ")\n";
            std::cout << "    Centroid: (" << boundary.centroid[0] << ", " 
                     << boundary.centroid[1] << ", " << boundary.centroid[2] << ")\n";
        }
        std::cout << "\n";
    }
}

} // namespace FSRM
