#ifndef SIMULATOR_HPP
#define SIMULATOR_HPP

#include "FSRM.hpp"
#include "PhysicsKernel.hpp"
#include "EclipseIO.hpp"
#include "WellModel.hpp"
#include "FractureModel.hpp"
#include "FaultModel.hpp"
#include "GmshIO.hpp"
#include "CoordinateSystem.hpp"
#include "ImplicitExplicitTransition.hpp"
#include "ConfigReader.hpp"
#include "SeismometerNetwork.hpp"
#include <memory>
#include <vector>

namespace FSRM {

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
    PetscErrorCode setupIMEXTransition(const ConfigReader::IMEXConfig& imex_config);
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

    // Seismometer output (stations/seismograms)
    std::unique_ptr<SeismometerNetwork> seismometers_;
    SeismometerOutputConfig seismo_out_cfg_;
    std::vector<SeismometerSpec> seismo_specs_;
    
    // IMEX transition manager for induced seismicity
    std::unique_ptr<ImplicitExplicitTransitionManager> imex_manager;
    ConfigReader::IMEXConfig imex_config;
    
    // Seismicity / fault network configuration (optional)
    ConfigReader::SeismicityConfig seismicity_config_;
    std::string seismic_catalog_file_;

    // Explosion coupling (used to drive faults and dynamics from explosion physics)
    // Implemented as a PIMPL in Simulator.cpp to avoid pulling large headers here.
    struct ExplosionCoupling;
    std::unique_ptr<ExplosionCoupling> explosion_;

    // Synthetic "nuclear trigger" stress pulse (optional, for historical nuclear test demos)
    bool nuclear_trigger_enabled_ = false;
    double nuclear_trigger_time0_ = 10.0;
    double nuclear_trigger_rise_ = 0.5;
    double nuclear_trigger_decay_ = 60.0;
    double nuclear_trigger_dsigma_ = 0.0;
    double nuclear_trigger_dtau_ = 2.0e6;
    double nuclear_trigger_dpore_ = 0.0;
    double last_fault_update_time_ = -1.0;

    // Estimated velocity for IMEX monitor (finite difference)
    Vec solution_prev_ = nullptr;
    Vec velocity_est_ = nullptr;
    double prev_time_ = 0.0;
    
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

    // If true, use PETSc FEM helpers (DMPlexTSComputeIFunctionFEM).
    // If false, fall back to a safe placeholder (U_t = 0) to avoid crashes when
    // the discretization/physics residuals are not configured.
    bool use_fem_time_residual_ = false;
    
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
