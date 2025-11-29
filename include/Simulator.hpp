#ifndef SIMULATOR_HPP
#define SIMULATOR_HPP

#include "ReservoirSim.hpp"
#include "PhysicsKernel.hpp"
#include "EclipseIO.hpp"
#include "WellModel.hpp"
#include "FractureModel.hpp"
#include "FaultModel.hpp"
#include "GmshIO.hpp"
#include "CoordinateSystem.hpp"
#include "ImplicitExplicitTransition.hpp"
#include <memory>
#include <vector>

namespace FSRM {

// Forward declaration
class ConfigReader;

class Simulator {
public:
    Simulator(MPI_Comm comm);
    ~Simulator();
    
    // Initialization
    PetscErrorCode initialize(const SimulationConfig& config);
    PetscErrorCode initializeFromConfigFile(const std::string& config_file);
    PetscErrorCode setupDM();
    PetscErrorCode setupFields();
    PetscErrorCode setupPhysics();
    PetscErrorCode setupTimeStepper();
    PetscErrorCode setupSolvers();
    
    // IMEX transition setup for induced seismicity modeling
    PetscErrorCode setupIMEXTransition(const IMEXConfig& imex_config);
    PetscErrorCode setupFaultNetwork();
    
    // Load input
    PetscErrorCode loadEclipseInput(const std::string& filename);
    PetscErrorCode loadMaterialPropertiesFromFile(const std::string& filename);
    PetscErrorCode setInitialConditions();
    PetscErrorCode setMaterialProperties();
    
    // Run simulation
    PetscErrorCode run();
    PetscErrorCode solve();
    
    // Physics management
    PetscErrorCode addPhysicsKernel(std::shared_ptr<PhysicsKernel> kernel);
    PetscErrorCode enableCoupling(PhysicsType type1, PhysicsType type2);
    
    // Well management
    PetscErrorCode addWell(const std::string& name, WellType type);
    PetscErrorCode setWellControl(const std::string& name, double target);
    
    // Fracture management
    PetscErrorCode addFracture(FractureType type, const std::vector<double>& coords);
    PetscErrorCode updateFractureNetwork();
    
    // Output
    PetscErrorCode writeOutput(int step);
    PetscErrorCode writeCheckpoint(int step);
    PetscErrorCode writeSummary();
    PetscErrorCode writeVisualization(int step);
    
    // Adaptive timestepping
    PetscErrorCode adaptTimeStep(double& dt);
    bool checkConvergence();
    
    // Analysis and diagnostics
    PetscErrorCode computePerformanceMetrics();
    PetscErrorCode generatePlots();
    PetscErrorCode writeStatistics();
    
private:
    MPI_Comm comm;
    int rank, size;
    
    // PETSc objects
    DM dm;
    TS ts;
    Vec solution, solution_old;
    Mat jacobian;
    SNES snes;
    KSP ksp;
    PC pc;
    
    // Fields
    std::vector<PetscFE> fe_fields;
    PetscDS prob;
    
    // Configuration
    SimulationConfig config;
    GridConfig grid_config;
    
    // Physics kernels
    std::vector<std::shared_ptr<PhysicsKernel>> kernels;
    std::map<PhysicsType, int> kernel_map;
    
    // Wells and fractures
    std::vector<std::shared_ptr<WellModel>> wells;
    std::vector<std::shared_ptr<FractureModel>> fractures;
    
    // Eclipse I/O
    std::unique_ptr<EclipseIO> eclipse_io;
    
    // Gmsh mesh I/O
    std::unique_ptr<GmshIO> gmsh_io;
    
    // Coordinate system manager
    std::unique_ptr<CoordinateSystemManager> coord_manager;
    
    // IMEX transition manager for induced seismicity
    std::unique_ptr<ImplicitExplicitTransitionManager> imex_manager;
    IMEXConfig imex_config;
    
    // Fault models for induced seismicity
    std::unique_ptr<FaultNetwork> fault_network;
    
    // Material properties (per cell or per region)
    std::vector<MaterialProperties> material_props;
    std::vector<FluidProperties> fluid_props;
    
    // Simulation state
    double current_time;
    double dt;
    int timestep;
    int output_counter;
    
    // Performance metrics
    std::vector<double> solve_times;
    std::vector<int> nonlinear_iterations;
    std::vector<int> linear_iterations;
    
    // Static callback functions for PETSc
    static PetscErrorCode FormFunction(TS ts, PetscReal t, Vec U, Vec U_t, Vec F, void *ctx);
    static PetscErrorCode FormJacobian(TS ts, PetscReal t, Vec U, Vec U_t, 
                                       PetscReal a, Mat J, Mat P, void *ctx);
    static PetscErrorCode MonitorFunction(TS ts, PetscInt step, PetscReal t, 
                                         Vec U, void *ctx);
    
    // Residual and Jacobian pointwise functions
    static void f0_SinglePhase(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                              const PetscInt uOff[], const PetscInt uOff_x[],
                              const PetscScalar u[], const PetscScalar u_t[],
                              const PetscScalar u_x[], const PetscInt aOff[],
                              const PetscInt aOff_x[], const PetscScalar a[],
                              const PetscScalar a_x[], const PetscScalar a_t[],
                              PetscReal t,
                              const PetscReal x[], PetscInt numConstants,
                              const PetscScalar constants[], PetscScalar f0[]);
    
    static void f1_SinglePhase(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                              const PetscInt uOff[], const PetscInt uOff_x[],
                              const PetscScalar u[], const PetscScalar u_t[],
                              const PetscScalar u_x[], const PetscInt aOff[],
                              const PetscInt aOff_x[], const PetscScalar a[],
                              const PetscScalar a_x[], const PetscScalar a_t[],
                              PetscReal t,
                              const PetscReal x[], PetscInt numConstants,
                              const PetscScalar constants[], PetscScalar f1[]);
    
    // Helper functions
    PetscErrorCode setupStructuredGrid();
    PetscErrorCode setupUnstructuredGrid();
    PetscErrorCode setupGmshGrid();
    PetscErrorCode setupCoordinateSystem();
    PetscErrorCode createFieldsFromConfig();
    PetscErrorCode setupBoundaryConditions();
    PetscErrorCode interpolateMaterialProperties();
    PetscErrorCode applyMaterialPropertiesToKernels();
    
    // Coordinate transformation helpers
    void transformToModelCoordinates(double& x, double& y, double& z) const;
    void transformToInputCoordinates(double& x, double& y, double& z) const;
};

} // namespace FSRM

#endif // SIMULATOR_HPP
