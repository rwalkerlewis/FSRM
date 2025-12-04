#ifndef GMSH_IO_HPP
#define GMSH_IO_HPP

#include "FSRM.hpp"
#include <string>
#include <vector>
#include <map>
#include <memory>
#include <array>
#include <set>
#include <functional>

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
 * @brief Material domain definition from physical groups
 * 
 * Maps Gmsh physical groups to material properties. This allows defining
 * different rock/material types in different regions of the mesh.
 */
struct MaterialDomain {
    int physical_tag;           ///< Physical group tag in Gmsh
    std::string name;           ///< Domain name (from physical group name)
    int dimension;              ///< Dimension (3 for volume, 2 for surface)
    std::string material_name;  ///< Reference to material properties
    int material_id;            ///< Material ID for lookup
    
    // Domain statistics
    int num_elements;           ///< Number of elements in domain
    double volume;              ///< Total volume (for 3D) or area (for 2D)
    std::array<double, 3> centroid;  ///< Domain centroid
    
    MaterialDomain() : 
        physical_tag(0), dimension(3), material_id(0),
        num_elements(0), volume(0.0), centroid({0,0,0}) {}
};

/**
 * @brief Fault surface definition from physical groups
 * 
 * Maps Gmsh physical groups representing fault surfaces. Faults are typically
 * defined as 2D surfaces (physical surfaces) embedded in a 3D domain, or as
 * 1D lines in a 2D domain.
 */
struct FaultSurface {
    int physical_tag;           ///< Physical group tag in Gmsh
    std::string name;           ///< Fault name (from physical group name)
    int dimension;              ///< Dimension (2 for surface, 1 for line)
    
    // Fault properties (can be overridden in config)
    double strike;              ///< Strike angle (radians, computed from geometry)
    double dip;                 ///< Dip angle (radians, computed from geometry)
    double length;              ///< Fault length (m)
    double width;               ///< Fault width/depth (m)
    double area;                ///< Total fault surface area (m²)
    
    // Element information
    std::vector<int> element_indices;   ///< Indices of fault elements
    std::vector<int> node_indices;      ///< Indices of fault nodes
    std::array<double, 3> centroid;     ///< Fault centroid
    std::array<double, 3> normal;       ///< Average fault normal vector
    
    // Split node information (for fault discontinuity)
    bool requires_split_nodes;          ///< True if fault needs split nodes
    std::vector<std::pair<int, int>> split_node_pairs;  ///< Plus/minus node pairs
    
    FaultSurface() : 
        physical_tag(0), dimension(2), 
        strike(0), dip(0), length(0), width(0), area(0),
        centroid({0,0,0}), normal({0,0,1}),
        requires_split_nodes(false) {}
};

/**
 * @brief Boundary surface definition from physical groups
 * 
 * Maps Gmsh physical groups to boundary conditions. Supports named boundaries
 * for easy reference in configuration files.
 */
struct BoundarySurface {
    int physical_tag;           ///< Physical group tag in Gmsh
    std::string name;           ///< Boundary name (from physical group name)
    int dimension;              ///< Dimension (2 for surface, 1 for line, 0 for point)
    std::string bc_type;        ///< Boundary condition type (DIRICHLET, NEUMANN, etc.)
    
    // Boundary statistics
    int num_elements;           ///< Number of elements on boundary
    double area;                ///< Total boundary area/length
    std::vector<int> element_indices;  ///< Indices of boundary elements
    std::vector<int> node_indices;     ///< Indices of boundary nodes
    std::array<double, 3> centroid;    ///< Boundary centroid
    std::array<double, 3> normal;      ///< Average outward normal (for surfaces)
    
    BoundarySurface() : 
        physical_tag(0), dimension(2), num_elements(0), area(0),
        centroid({0,0,0}), normal({0,0,1}) {}
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
 * 
 * Contains all mesh data including nodes, elements, physical groups,
 * and derived structures for material domains, faults, and boundaries.
 */
struct UnstructuredMesh {
    std::vector<MeshNode> nodes;
    std::vector<MeshElement> elements;
    std::map<int, PhysicalGroup> physical_groups;
    
    // Derived data
    int dimension;                      ///< Mesh dimension (2D or 3D)
    std::array<double, 3> bbox_min;     ///< Bounding box minimum
    std::array<double, 3> bbox_max;     ///< Bounding box maximum
    
    // Material domains (from physical volumes)
    std::map<std::string, MaterialDomain> material_domains;
    
    // Fault surfaces (from physical surfaces marked as faults)
    std::map<std::string, FaultSurface> fault_surfaces;
    
    // Boundary surfaces (from physical surfaces for BCs)
    std::map<std::string, BoundarySurface> boundary_surfaces;
    
    // Element-to-domain mapping (element index -> material domain name)
    std::vector<int> element_material_ids;
    
    UnstructuredMesh() : dimension(3), bbox_min({0,0,0}), bbox_max({0,0,0}) {}
    
    // Query methods
    int getNumNodes() const { return static_cast<int>(nodes.size()); }
    int getNumElements() const { return static_cast<int>(elements.size()); }
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
    
    // Material domain methods
    MaterialDomain* getMaterialDomain(const std::string& name);
    const MaterialDomain* getMaterialDomain(const std::string& name) const;
    std::vector<std::string> getMaterialDomainNames() const;
    int getMaterialIdForElement(int element_idx) const;
    
    // Fault surface methods
    FaultSurface* getFaultSurface(const std::string& name);
    const FaultSurface* getFaultSurface(const std::string& name) const;
    std::vector<std::string> getFaultNames() const;
    
    // Boundary surface methods
    BoundarySurface* getBoundarySurface(const std::string& name);
    const BoundarySurface* getBoundarySurface(const std::string& name) const;
    std::vector<std::string> getBoundaryNames() const;
    std::vector<int> getBoundaryNodes(const std::string& name) const;
};

/**
 * @brief Configuration for material domain mapping
 */
struct GmshMaterialMapping {
    std::string physical_group_name;  ///< Name of physical group in Gmsh
    std::string material_section;     ///< Config section name (e.g., "ROCK1")
    int material_id;                  ///< Numeric ID for fast lookup
    
    GmshMaterialMapping() : material_id(0) {}
    GmshMaterialMapping(const std::string& pg, const std::string& mat, int id = 0)
        : physical_group_name(pg), material_section(mat), material_id(id) {}
};

/**
 * @brief Configuration for fault mapping
 */
struct GmshFaultMapping {
    std::string physical_group_name;  ///< Name of physical group in Gmsh
    std::string fault_section;        ///< Config section name (e.g., "FAULT1")
    bool use_split_nodes;             ///< Whether to use split nodes for this fault
    
    GmshFaultMapping() : use_split_nodes(false) {}
    GmshFaultMapping(const std::string& pg, const std::string& fault, bool split = false)
        : physical_group_name(pg), fault_section(fault), use_split_nodes(split) {}
};

/**
 * @brief Gmsh mesh file reader/writer
 * 
 * Supports Gmsh MSH file format version 2.2 and 4.x
 * 
 * Features:
 * - Load unstructured meshes from Gmsh .msh files
 * - Map physical groups to material domains
 * - Map physical surfaces to faults
 * - Map physical surfaces to boundary conditions
 * - Create PETSc DMPlex for simulation
 * - Validate mesh quality and consistency
 * 
 * Usage:
 * @code
 * GmshIO gmsh;
 * if (gmsh.readMshFile("mesh.msh")) {
 *     // Map physical groups to materials
 *     gmsh.setMaterialMapping({
 *         {"reservoir", "ROCK1", 0},
 *         {"caprock", "ROCK2", 1},
 *         {"aquifer", "ROCK3", 2}
 *     });
 *     
 *     // Map physical groups to faults
 *     gmsh.setFaultMapping({
 *         {"main_fault", "FAULT1", true},
 *         {"secondary_fault", "FAULT2", false}
 *     });
 *     
 *     // Process domains and create DMPlex
 *     gmsh.processDomains();
 *     
 *     DM dm;
 *     gmsh.createDMPlexWithMaterials(PETSC_COMM_WORLD, &dm);
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
    
    /**
     * @brief Write mesh with material IDs to VTK format
     * @param filename Output file path (.vtu)
     * @return true if successful
     */
    bool writeVtkFileWithMaterials(const std::string& filename) const;
    
    // =========================================================================
    // Domain Mapping Configuration
    // =========================================================================
    
    /**
     * @brief Set material domain mappings
     * 
     * Maps physical group names from the Gmsh file to material sections
     * in the configuration file.
     * 
     * @param mappings Vector of material mappings
     */
    void setMaterialMapping(const std::vector<GmshMaterialMapping>& mappings);
    
    /**
     * @brief Add a single material mapping
     * @param physical_group Name of physical group in Gmsh file
     * @param material_section Config section name (e.g., "ROCK1")
     * @param material_id Numeric ID for fast lookup
     */
    void addMaterialMapping(const std::string& physical_group, 
                           const std::string& material_section,
                           int material_id = -1);
    
    /**
     * @brief Set fault surface mappings
     * 
     * Maps physical group names representing fault surfaces to fault
     * sections in the configuration file.
     * 
     * @param mappings Vector of fault mappings
     */
    void setFaultMapping(const std::vector<GmshFaultMapping>& mappings);
    
    /**
     * @brief Add a single fault mapping
     * @param physical_group Name of physical group in Gmsh file
     * @param fault_section Config section name (e.g., "FAULT1")
     * @param use_split_nodes Whether to create split nodes for this fault
     */
    void addFaultMapping(const std::string& physical_group,
                        const std::string& fault_section,
                        bool use_split_nodes = false);
    
    /**
     * @brief Set boundary surface names
     * 
     * Specifies which physical groups represent boundary surfaces
     * for boundary condition application.
     * 
     * @param boundary_names Names of physical groups to use as boundaries
     */
    void setBoundaryNames(const std::vector<std::string>& boundary_names);
    
    /**
     * @brief Process all domain mappings
     * 
     * Call this after setting material/fault/boundary mappings to
     * populate the derived mesh structures (material_domains, fault_surfaces,
     * boundary_surfaces, element_material_ids).
     * 
     * @return true if successful
     */
    bool processDomains();
    
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
    
    /**
     * @brief Get list of all physical group names
     */
    std::vector<std::string> getPhysicalGroupNames() const;
    
    /**
     * @brief Get physical groups by dimension
     * @param dim Dimension (0=points, 1=lines, 2=surfaces, 3=volumes)
     */
    std::vector<std::string> getPhysicalGroupsByDimension(int dim) const;
    
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
    
    /**
     * @brief Create DMPlex with material labels
     * 
     * Creates a DMPlex with labels for material domains, allowing
     * material properties to be set per-region.
     * 
     * @param comm MPI communicator
     * @param dm Output DM object
     * @return PETSc error code
     */
    PetscErrorCode createDMPlexWithMaterials(MPI_Comm comm, DM* dm) const;
    
    /**
     * @brief Create DMPlex with material and fault labels
     * 
     * Creates a DMPlex with full labeling for materials, faults, and boundaries.
     * This is the most complete DM creation method.
     * 
     * @param comm MPI communicator
     * @param dm Output DM object
     * @return PETSc error code
     */
    PetscErrorCode createDMPlexFull(MPI_Comm comm, DM* dm) const;
    
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
     * @brief Print detailed domain information
     */
    void printDomainInfo() const;
    
    /**
     * @brief Validate mesh consistency
     * @return true if mesh is valid
     */
    bool validate() const;
    
    /**
     * @brief Validate mesh quality
     * @param min_quality Minimum acceptable element quality (0-1)
     * @param max_aspect_ratio Maximum acceptable aspect ratio
     * @return true if all elements pass quality checks
     */
    bool validateQuality(double min_quality = 0.1, double max_aspect_ratio = 10.0) const;
    
    /**
     * @brief Get mesh quality statistics
     * @param min_qual Output: minimum element quality
     * @param max_qual Output: maximum element quality
     * @param avg_qual Output: average element quality
     */
    void getQualityStats(double& min_qual, double& max_qual, double& avg_qual) const;
    
    /**
     * @brief Compute element quality (0 = degenerate, 1 = ideal)
     * @param element_idx Element index
     * @return Quality metric
     */
    double computeElementQuality(int element_idx) const;
    
    /**
     * @brief Compute element volume/area
     * @param element_idx Element index
     * @return Volume (3D) or area (2D)
     */
    double computeElementVolume(int element_idx) const;
    
    /**
     * @brief Compute surface element normal vector
     * @param element_idx Element index (must be a surface element)
     * @return Normal vector (unit length)
     */
    std::array<double, 3> computeElementNormal(int element_idx) const;
    
private:
    std::shared_ptr<UnstructuredMesh> mesh_;
    std::string format_version_;
    bool is_binary_;
    
    // Domain mappings
    std::vector<GmshMaterialMapping> material_mappings_;
    std::vector<GmshFaultMapping> fault_mappings_;
    std::vector<std::string> boundary_names_;
    
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
    
    // Domain processing helpers
    bool processMaterialDomains();
    bool processFaultSurfaces();
    bool processBoundarySurfaces();
    void computeDomainStatistics(MaterialDomain& domain);
    void computeFaultGeometry(FaultSurface& fault);
    void computeBoundaryStatistics(BoundarySurface& boundary);
    
    // Helper functions
    int getNumNodesForElementType(GmshElementType type) const;
    GmshElementType intToElementType(int type_id) const;
    std::array<double, 3> getNodeCoords(int node_tag) const;
    int findNodeIndex(int node_tag) const;
};

} // namespace FSRM

#endif // GMSH_IO_HPP
