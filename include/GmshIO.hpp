#ifndef GMSH_IO_HPP
#define GMSH_IO_HPP

#include "ReservoirSim.hpp"
#include <string>
#include <vector>
#include <map>
#include <memory>
#include <array>

namespace FSRM {

/**
 * @brief Node in the mesh
 */
struct MeshNode {
    int tag;           ///< Node identifier
    double x, y, z;    ///< Coordinates
    
    MeshNode() : tag(0), x(0), y(0), z(0) {}
    MeshNode(int t, double xx, double yy, double zz) : tag(t), x(xx), y(yy), z(zz) {}
};

/**
 * @brief Element types supported by Gmsh
 */
enum class GmshElementType {
    LINE = 1,           ///< 2-node line
    TRIANGLE = 2,       ///< 3-node triangle
    QUAD = 3,           ///< 4-node quadrangle
    TETRAHEDRON = 4,    ///< 4-node tetrahedron
    HEXAHEDRON = 5,     ///< 8-node hexahedron
    PRISM = 6,          ///< 6-node prism
    PYRAMID = 7,        ///< 5-node pyramid
    LINE3 = 8,          ///< 3-node second order line
    TRIANGLE6 = 9,      ///< 6-node second order triangle
    QUAD9 = 10,         ///< 9-node second order quadrangle
    TET10 = 11,         ///< 10-node second order tetrahedron
    HEX27 = 12,         ///< 27-node second order hexahedron
    PRISM18 = 13,       ///< 18-node second order prism
    PYRAMID14 = 14,     ///< 14-node second order pyramid
    POINT = 15,         ///< 1-node point
    QUAD8 = 16,         ///< 8-node second order quadrangle
    HEX20 = 17,         ///< 20-node second order hexahedron
    PRISM15 = 18,       ///< 15-node second order prism
    PYRAMID13 = 19,     ///< 13-node second order pyramid
    UNKNOWN = 0
};

/**
 * @brief Mesh element
 */
struct MeshElement {
    int tag;                       ///< Element identifier
    GmshElementType type;          ///< Element type
    int physical_tag;              ///< Physical group tag
    int geometrical_tag;           ///< Geometrical entity tag
    std::vector<int> node_tags;    ///< Node tags that form this element
    
    MeshElement() : tag(0), type(GmshElementType::UNKNOWN), physical_tag(0), geometrical_tag(0) {}
    
    int getNumNodes() const;
    int getDimension() const;
};

/**
 * @brief Physical group definition
 */
struct PhysicalGroup {
    int tag;             ///< Physical group tag
    int dimension;       ///< Dimension (0=point, 1=line, 2=surface, 3=volume)
    std::string name;    ///< Group name
    
    PhysicalGroup() : tag(0), dimension(0) {}
    PhysicalGroup(int t, int d, const std::string& n) : tag(t), dimension(d), name(n) {}
};

/**
 * @brief Unstructured mesh data container
 */
struct UnstructuredMesh {
    std::vector<MeshNode> nodes;
    std::vector<MeshElement> elements;
    std::map<int, PhysicalGroup> physical_groups;
    
    // Derived data
    int dimension;                      ///< Mesh dimension (2D or 3D)
    std::array<double, 3> bbox_min;     ///< Bounding box minimum
    std::array<double, 3> bbox_max;     ///< Bounding box maximum
    
    UnstructuredMesh() : dimension(3), bbox_min({0,0,0}), bbox_max({0,0,0}) {}
    
    // Query methods
    int getNumNodes() const { return nodes.size(); }
    int getNumElements() const { return elements.size(); }
    int getNumCells() const;  // Volume elements only
    int getNumFaces() const;  // Surface elements only
    
    // Bounding box
    void computeBoundingBox();
    double getLx() const { return bbox_max[0] - bbox_min[0]; }
    double getLy() const { return bbox_max[1] - bbox_min[1]; }
    double getLz() const { return bbox_max[2] - bbox_min[2]; }
    
    // Filter by physical group
    std::vector<int> getElementsByPhysicalTag(int tag) const;
    std::vector<int> getNodesByPhysicalTag(int tag) const;
    
    // Get physical group by name
    int getPhysicalTagByName(const std::string& name) const;
};

/**
 * @brief Gmsh mesh file reader/writer
 * 
 * Supports Gmsh MSH file format version 2.2 and 4.x
 * 
 * Usage:
 * @code
 * GmshIO gmsh;
 * if (gmsh.readMshFile("mesh.msh")) {
 *     auto mesh = gmsh.getMesh();
 *     std::cout << "Loaded " << mesh->getNumNodes() << " nodes\n";
 * }
 * @endcode
 */
class GmshIO {
public:
    GmshIO();
    ~GmshIO() = default;
    
    // =========================================================================
    // File I/O
    // =========================================================================
    
    /**
     * @brief Read a Gmsh .msh file
     * @param filename Path to the .msh file
     * @return true if successful
     */
    bool readMshFile(const std::string& filename);
    
    /**
     * @brief Write mesh to Gmsh .msh file (v4.1 format)
     * @param filename Output file path
     * @return true if successful
     */
    bool writeMshFile(const std::string& filename) const;
    
    /**
     * @brief Write mesh to VTK format for visualization
     * @param filename Output file path (.vtu)
     * @return true if successful
     */
    bool writeVtkFile(const std::string& filename) const;
    
    // =========================================================================
    // Mesh Access
    // =========================================================================
    
    /**
     * @brief Get the loaded mesh
     * @return Shared pointer to the mesh data
     */
    std::shared_ptr<UnstructuredMesh> getMesh() const { return mesh_; }
    
    /**
     * @brief Check if mesh is loaded
     */
    bool isLoaded() const { return mesh_ != nullptr && !mesh_->nodes.empty(); }
    
    /**
     * @brief Get mesh dimension (2 or 3)
     */
    int getDimension() const { return mesh_ ? mesh_->dimension : 0; }
    
    // =========================================================================
    // PETSc Integration
    // =========================================================================
    
    /**
     * @brief Create a PETSc DMPlex from the loaded mesh
     * @param comm MPI communicator
     * @param dm Output DM object (must be destroyed by caller)
     * @return PETSc error code
     */
    PetscErrorCode createDMPlex(MPI_Comm comm, DM* dm) const;
    
    /**
     * @brief Create DMPlex with specified physical groups for boundaries
     * @param comm MPI communicator
     * @param dm Output DM object
     * @param boundary_groups Physical group names to use as boundaries
     * @return PETSc error code
     */
    PetscErrorCode createDMPlexWithBoundaries(
        MPI_Comm comm, 
        DM* dm,
        const std::vector<std::string>& boundary_groups
    ) const;
    
    // =========================================================================
    // Utility Methods
    // =========================================================================
    
    /**
     * @brief Get the Gmsh file format version
     */
    std::string getFormatVersion() const { return format_version_; }
    
    /**
     * @brief Print mesh statistics
     */
    void printStatistics() const;
    
    /**
     * @brief Validate mesh consistency
     * @return true if mesh is valid
     */
    bool validate() const;
    
private:
    std::shared_ptr<UnstructuredMesh> mesh_;
    std::string format_version_;
    bool is_binary_;
    
    // Version-specific readers
    bool readMsh2(std::ifstream& file);
    bool readMsh4(std::ifstream& file);
    bool readMsh41(std::ifstream& file);
    
    // Section readers for MSH 4.x format
    bool readMeshFormat(std::ifstream& file);
    bool readPhysicalNames(std::ifstream& file);
    bool readEntities(std::ifstream& file);
    bool readNodes(std::ifstream& file);
    bool readElements(std::ifstream& file);
    
    // Section readers for MSH 2.x format
    bool readNodes2(std::ifstream& file);
    bool readElements2(std::ifstream& file);
    bool readPhysicalNames2(std::ifstream& file);
    
    // Helper functions
    int getNumNodesForElementType(GmshElementType type) const;
    GmshElementType intToElementType(int type_id) const;
};

} // namespace FSRM

#endif // GMSH_IO_HPP
