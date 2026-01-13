/**
 * @file gmsh_mesh_domains.cpp
 * @brief Example demonstrating Gmsh mesh loading with material domains and fault surfaces
 * 
 * This example shows how to:
 * 1. Load a Gmsh mesh file
 * 2. Configure material domain mappings
 * 3. Configure fault surface mappings
 * 4. Process and validate the mesh
 * 5. Access domain and fault information
 * 6. Create a PETSc DMPlex for simulation
 * 
 * Prerequisites:
 *   - Gmsh mesh file with physical groups defined
 *   - PETSc and MPI installed
 * 
 * Compile:
 *   mpicxx -o gmsh_mesh_domains gmsh_mesh_domains.cpp \
 *     -I../include -L../lib -lfsrm \
 *     $(pkg-config --cflags --libs petsc)
 * 
 * Run:
 *   mpirun -np 4 ./gmsh_mesh_domains meshes/my_mesh.msh
 */

#include "GmshIO.hpp"
#include "ConfigReader.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/**
 * @brief Print a separator line for output formatting
 */
void printSeparator(const std::string& title) {
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << " " << title << "\n";
    std::cout << std::string(70, '=') << "\n\n";
}

/**
 * @brief Print mesh statistics
 */
void printMeshStats(const FSRM::GmshIO& gmsh) {
    auto mesh = gmsh.getMesh();
    if (!mesh) return;
    
    std::cout << "Mesh Format: " << gmsh.getFormatVersion() << "\n";
    std::cout << "Dimension: " << mesh->dimension << "D\n";
    std::cout << "Nodes: " << mesh->getNumNodes() << "\n";
    std::cout << "Elements: " << mesh->getNumElements() << "\n";
    std::cout << "Cells: " << mesh->getNumCells() << "\n";
    std::cout << "Faces: " << mesh->getNumFaces() << "\n";
    
    std::cout << "\nBounding Box:\n";
    std::cout << "  X: [" << mesh->bbox_min[0] << ", " << mesh->bbox_max[0] << "]\n";
    std::cout << "  Y: [" << mesh->bbox_min[1] << ", " << mesh->bbox_max[1] << "]\n";
    std::cout << "  Z: [" << mesh->bbox_min[2] << ", " << mesh->bbox_max[2] << "]\n";
    
    std::cout << "\nPhysical Groups:\n";
    for (const auto& [tag, group] : mesh->physical_groups) {
        std::cout << "  " << group.name 
                  << " (dim=" << group.dimension 
                  << ", tag=" << tag << ")\n";
    }
}

/**
 * @brief Print material domain information
 */
void printMaterialDomains(const FSRM::GmshIO& gmsh) {
    auto mesh = gmsh.getMesh();
    if (!mesh) return;
    
    for (const auto& name : mesh->getMaterialDomainNames()) {
        auto* domain = mesh->getMaterialDomain(name);
        if (!domain) continue;
        
        std::cout << "Domain: " << name << "\n";
        std::cout << "  Physical Tag: " << domain->physical_tag << "\n";
        std::cout << "  Material: " << domain->material_name << " (ID=" << domain->material_id << ")\n";
        std::cout << "  Dimension: " << domain->dimension << "D\n";
        std::cout << "  Elements: " << domain->num_elements << "\n";
        std::cout << "  Volume: " << std::scientific << std::setprecision(4) 
                  << domain->volume << " m³\n";
        std::cout << "  Centroid: (" << std::fixed << std::setprecision(2)
                  << domain->centroid[0] << ", "
                  << domain->centroid[1] << ", "
                  << domain->centroid[2] << ")\n\n";
    }
}

/**
 * @brief Print fault surface information
 */
void printFaultSurfaces(const FSRM::GmshIO& gmsh) {
    auto mesh = gmsh.getMesh();
    if (!mesh) return;
    
    for (const auto& name : mesh->getFaultNames()) {
        auto* fault = mesh->getFaultSurface(name);
        if (!fault) continue;
        
        std::cout << "Fault: " << name << "\n";
        std::cout << "  Physical Tag: " << fault->physical_tag << "\n";
        std::cout << "  Dimension: " << fault->dimension << "D\n";
        std::cout << "  Elements: " << fault->element_indices.size() << "\n";
        std::cout << "  Nodes: " << fault->node_indices.size() << "\n";
        std::cout << "  Area: " << std::scientific << std::setprecision(4) 
                  << fault->area << " m²\n";
        std::cout << "  Length: " << std::fixed << std::setprecision(1) 
                  << fault->length << " m\n";
        std::cout << "  Width: " << fault->width << " m\n";
        std::cout << "  Strike: " << fault->strike * 180.0 / M_PI << "°\n";
        std::cout << "  Dip: " << fault->dip * 180.0 / M_PI << "°\n";
        std::cout << "  Normal: (" << std::fixed << std::setprecision(3)
                  << fault->normal[0] << ", "
                  << fault->normal[1] << ", "
                  << fault->normal[2] << ")\n";
        std::cout << "  Centroid: (" << std::setprecision(1)
                  << fault->centroid[0] << ", "
                  << fault->centroid[1] << ", "
                  << fault->centroid[2] << ")\n";
        std::cout << "  Split Nodes: " << (fault->requires_split_nodes ? "yes" : "no") << "\n\n";
    }
}

/**
 * @brief Print boundary surface information
 */
void printBoundarySurfaces(const FSRM::GmshIO& gmsh) {
    auto mesh = gmsh.getMesh();
    if (!mesh) return;
    
    for (const auto& name : mesh->getBoundaryNames()) {
        auto* boundary = mesh->getBoundarySurface(name);
        if (!boundary) continue;
        
        std::cout << "Boundary: " << name << "\n";
        std::cout << "  Elements: " << boundary->num_elements << "\n";
        std::cout << "  Nodes: " << boundary->node_indices.size() << "\n";
        std::cout << "  Area: " << std::scientific << std::setprecision(4) 
                  << boundary->area << " m²\n";
        std::cout << "  Normal: (" << std::fixed << std::setprecision(3)
                  << boundary->normal[0] << ", "
                  << boundary->normal[1] << ", "
                  << boundary->normal[2] << ")\n\n";
    }
}

/**
 * @brief Print mesh quality statistics
 */
void printQualityStats(const FSRM::GmshIO& gmsh) {
    double min_q, max_q, avg_q;
    gmsh.getQualityStats(min_q, max_q, avg_q);
    
    std::cout << "Quality Statistics:\n";
    std::cout << "  Minimum: " << std::fixed << std::setprecision(4) << min_q << "\n";
    std::cout << "  Maximum: " << max_q << "\n";
    std::cout << "  Average: " << avg_q << "\n";
    
    // Quality assessment
    if (min_q < 0.1) {
        std::cout << "  WARNING: Some elements have poor quality (< 0.1)\n";
    } else if (min_q < 0.3) {
        std::cout << "  NOTICE: Some elements have moderate quality (< 0.3)\n";
    } else {
        std::cout << "  Good mesh quality (minimum > 0.3)\n";
    }
}

int main(int argc, char** argv) {
    // Initialize PETSc and MPI
    PetscErrorCode ierr = PetscInitialize(&argc, &argv, nullptr, nullptr);
    CHKERRQ(ierr);
    
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    
    // Check command line arguments
    std::string mesh_file = "meshes/reservoir.msh";  // Default
    if (argc > 1) {
        mesh_file = argv[1];
    }
    
    if (rank == 0) {
        std::cout << "======================================================\n";
        std::cout << "  FSRM Gmsh Mesh Domain Example\n";
        std::cout << "======================================================\n\n";
        std::cout << "Mesh file: " << mesh_file << "\n";
    }
    
    // Create GmshIO instance
    FSRM::GmshIO gmsh;
    
    // Load mesh file
    printSeparator("Loading Mesh");
    if (!gmsh.readMshFile(mesh_file)) {
        if (rank == 0) {
            std::cerr << "ERROR: Failed to load mesh file: " << mesh_file << "\n";
            std::cerr << "\nUsage: " << argv[0] << " <mesh_file.msh>\n";
            std::cerr << "\nTo create a test mesh, use Gmsh:\n";
            std::cerr << "  gmsh -3 test.geo -o test.msh\n";
        }
        PetscFinalize();
        return 1;
    }
    
    if (rank == 0) {
        std::cout << "Mesh loaded successfully!\n\n";
        printMeshStats(gmsh);
    }
    
    // Configure material domain mappings
    printSeparator("Configuring Material Domains");
    
    // These names should match Physical Volume names in your Gmsh file
    // Adjust based on your actual mesh
    auto physical_groups = gmsh.getPhysicalGroupNames();
    
    if (rank == 0) {
        std::cout << "Available physical groups:\n";
        for (const auto& name : physical_groups) {
            std::cout << "  - " << name << "\n";
        }
        std::cout << "\n";
    }
    
    // Example material mappings (adjust to match your mesh)
    // Format: addMaterialMapping(physical_group_name, material_section, material_id)
    gmsh.addMaterialMapping("reservoir", "ROCK1", 0);
    gmsh.addMaterialMapping("caprock", "ROCK2", 1);
    gmsh.addMaterialMapping("basement", "ROCK3", 2);
    gmsh.addMaterialMapping("aquifer", "ROCK4", 3);
    
    // Alternative: set all mappings at once
    // gmsh.setMaterialMapping({
    //     {"reservoir", "ROCK1", 0},
    //     {"caprock", "ROCK2", 1},
    //     {"basement", "ROCK3", 2}
    // });
    
    // Configure fault surface mappings
    printSeparator("Configuring Fault Surfaces");
    
    // Example fault mappings
    // Format: addFaultMapping(physical_group_name, fault_section, use_split_nodes)
    gmsh.addFaultMapping("main_fault", "FAULT1", true);    // With split nodes
    gmsh.addFaultMapping("secondary_fault", "FAULT2", false);  // Without split nodes
    
    // Configure boundary surfaces
    printSeparator("Configuring Boundary Surfaces");
    
    gmsh.setBoundaryNames({"inlet", "outlet", "top", "bottom", "north", "south"});
    
    // Process all domain mappings
    printSeparator("Processing Domains");
    
    if (!gmsh.processDomains()) {
        if (rank == 0) {
            std::cerr << "WARNING: Some domains could not be processed.\n";
            std::cerr << "This is normal if your mesh doesn't have all the expected physical groups.\n";
        }
    }
    
    // Print detailed domain information
    printSeparator("Material Domains");
    if (rank == 0) {
        printMaterialDomains(gmsh);
    }
    
    printSeparator("Fault Surfaces");
    if (rank == 0) {
        printFaultSurfaces(gmsh);
    }
    
    printSeparator("Boundary Surfaces");
    if (rank == 0) {
        printBoundarySurfaces(gmsh);
    }
    
    // Validate mesh quality
    printSeparator("Mesh Quality Analysis");
    if (rank == 0) {
        printQualityStats(gmsh);
        
        // Check quality threshold
        if (!gmsh.validateQuality(0.1, 10.0)) {
            std::cout << "\nWARNING: Some elements may cause numerical issues.\n";
            std::cout << "Consider re-meshing with better quality settings.\n";
        }
    }
    
    // Create PETSc DMPlex
    printSeparator("Creating PETSc DMPlex");
    
    DM dm;
    ierr = gmsh.createDMPlexFull(PETSC_COMM_WORLD, &dm); CHKERRQ(ierr);
    
    if (rank == 0) {
        std::cout << "DMPlex created with labels:\n";
        std::cout << "  - Material (material domain IDs)\n";
        std::cout << "  - Fault (fault surface IDs)\n";
        std::cout << "  - Boundary (boundary surface IDs)\n";
    }
    
    // Write VTK output for visualization
    printSeparator("Writing VTK Output");
    
    if (rank == 0) {
        std::string vtk_file = mesh_file.substr(0, mesh_file.find_last_of('.')) + "_domains.vtu";
        if (gmsh.writeVtkFileWithMaterials(vtk_file)) {
            std::cout << "VTK file written: " << vtk_file << "\n";
            std::cout << "Open in ParaView to visualize material domains.\n";
        } else {
            std::cerr << "Failed to write VTK file.\n";
        }
    }
    
    // Example: Access mesh data for custom processing
    printSeparator("Custom Data Access Example");
    
    if (rank == 0) {
        auto mesh = gmsh.getMesh();
        
        // Get material ID for each element
        std::cout << "Element-Material ID mapping (first 10 elements):\n";
        for (int i = 0; i < std::min(10, mesh->getNumElements()); i++) {
            int mat_id = mesh->getMaterialIdForElement(i);
            std::cout << "  Element " << i << " -> Material ID " << mat_id << "\n";
        }
        
        // Get boundary nodes
        std::cout << "\nBoundary node counts:\n";
        for (const auto& name : mesh->getBoundaryNames()) {
            auto nodes = mesh->getBoundaryNodes(name);
            std::cout << "  " << name << ": " << nodes.size() << " nodes\n";
        }
    }
    
    // Cleanup
    ierr = DMDestroy(&dm); CHKERRQ(ierr);
    
    if (rank == 0) {
        printSeparator("Done");
        std::cout << "Example completed successfully.\n";
        std::cout << "\nNext steps:\n";
        std::cout << "  1. Create your own Gmsh mesh with physical groups\n";
        std::cout << "  2. Configure material/fault/boundary mappings\n";
        std::cout << "  3. Use with FSRM Simulator class for full simulation\n\n";
    }
    
    PetscFinalize();
    return 0;
}
