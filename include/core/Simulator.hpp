#ifndef SIMULATOR_HPP
#define SIMULATOR_HPP

#include "core/FSRM.hpp"
#include "physics/PhysicsKernel.hpp"
#include "io/EclipseIO.hpp"
#include "domain/reservoir/WellModel.hpp"
#include "domain/geomechanics/FractureModel.hpp"
#include "domain/geomechanics/FaultModel.hpp"
#include "io/GmshIO.hpp"
#include "util/CoordinateSystem.hpp"
#include "numerics/ImplicitExplicitTransition.hpp"
#include "core/ConfigReader.hpp"
#include "domain/seismic/SeismometerNetwork.hpp"
#include "domain/geomechanics/PyLithFault.hpp"
#include "numerics/FaultMeshManager.hpp"
#include "physics/CohesiveFaultKernel.hpp"
#include "domain/geomechanics/CoulombStressTransfer.hpp"
#include "core/DerivedFieldComputer.hpp"
#include "io/VelocityModelReader.hpp"
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
    PetscErrorCode labelBoundaries();
    PetscErrorCode setupBoundaryConditions();
    
    // Load input
    PetscErrorCode loadEclipseInput(const std::string& filename);
    PetscErrorCode loadMaterialPropertiesFromFile(const std::string& filename);
    PetscErrorCode setInitialConditions();
    PetscErrorCode applyInitialFaultStress();
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

    // Accessors for testing
    DM getDM() const { return dm; }
    TS getTS() const { return ts; }
    Vec getSolution() const { return solution; }
    DM getAuxDM() const { return auxDM_; }
    Vec getAuxVector() const { return auxVec_; }
    const SimulationConfig& getConfig() const { return config; }
    const DerivedFieldComputer& getDerivedFields() const { return derived_fields_; }
    
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
    
    // Cohesive fault for dynamic rupture (Branch 3)
    std::unique_ptr<FaultCohesiveDyn> cohesive_fault_;
    std::unique_ptr<FaultMeshManager> fault_mesh_manager_;
    std::unique_ptr<CohesiveFaultKernel> cohesive_kernel_;
    std::unique_ptr<CoulombStressTransfer> cfs_transfer_;

    // Derived field computation (stress, strain, CFS)
    DerivedFieldComputer derived_fields_;
    
    // Material properties (per cell or per region)
    std::vector<MaterialProperties> material_props;
    std::vector<FluidProperties> fluid_props;

    // Auxiliary DM for heterogeneous material properties
    DM  auxDM_ = nullptr;
    Vec auxVec_ = nullptr;
    PetscErrorCode setupAuxiliaryDM();
    PetscErrorCode populateAuxFieldsByDepth();
    PetscErrorCode populateAuxFieldsByMaterialLabel();
    PetscErrorCode populateAuxFieldsByVelocityModel();
    PetscErrorCode applyExplosionDamageToAuxFields();
    PetscErrorCode createCohesiveCellLabel();
    
    // Simulation state
    double current_time;
    double dt;
    int timestep;
    int output_counter;

    // If true, use PETSc FEM helpers (DMPlexTSComputeIFunctionFEM).
    // If false, fall back to a safe placeholder (U_t = 0) to avoid crashes when
    // the discretization/physics residuals are not configured.
    bool use_fem_time_residual_ = false;

    // Injection source (Phase 2)
    bool injection_enabled_ = false;
    double injection_x_ = 0, injection_y_ = 0, injection_z_ = 0;
    double injection_rate_ = 0;
    double injection_start_ = 0, injection_end_ = 1e30;
    PetscInt injection_cell_ = -1;
    PetscErrorCode locateInjectionCell();
    PetscErrorCode addInjectionToResidual(PetscReal t, Vec locF);
    PetscErrorCode addFaultPressureToResidual(PetscReal t, Vec locF);
    PetscErrorCode addTractionToResidual(PetscReal t, Vec locF);
    PetscErrorCode subtractLagrangeRegularizationOnCohesive(Vec locF, Vec locU);
    PetscErrorCode zeroLagrangeDiagonalOnCohesive(Mat J);
    PetscErrorCode addCohesivePenaltyToJacobian(Mat J, Vec locU = nullptr);

    // Hydraulic fracture model (Phase 3)
    std::unique_ptr<HydraulicFractureModel> hydrofrac_;
    bool hydrofrac_initiated_ = false;
    bool hydrofrac_fem_pressurized_mode_ = false;
    double hydrofrac_uniform_pressure_pa_ = 0.0;

    // Traction boundary conditions (manual assembly)
    bool has_traction_bc_ = false;

    // Cohesive fracture plane (Phase 4)
    bool fracture_plane_enabled_ = false;
    double fracture_plane_strike_ = 0.0;   // radians
    double fracture_plane_dip_ = M_PI / 2.0;
    double fracture_plane_center_[3] = {0, 0, 0};
    double fracture_plane_length_ = 200.0;
    double fracture_plane_width_ = 100.0;
    double fracture_plane_tensile_strength_ = 5.0e6;  // Pa

    // Fault mode, geometry, and prescribed slip
    std::string fault_mode_ = "locked";
    double fault_center_[3] = {-1.0, -1.0, -1.0};  // negative = use grid center default
    double fault_strike_ = 0.0;    // radians; negative = use default (0)
    double fault_dip_ = -1.0;      // radians; negative = use default (pi/2)
    double fault_length_ = -1.0;   // m; negative = use 2*max(Lx,Ly,Lz)
    double fault_width_ = -1.0;    // m; negative = use 2*max(Lx,Ly,Lz)
    bool fault_geometry_from_config_ = false;
    double fault_slip_strike_ = 0.0;
    double fault_slip_dip_ = 0.0;
    double fault_slip_opening_ = 0.0;
    double fault_friction_coefficient_ = 0.6;
    bool use_slip_weakening_ = false;
    double mu_static_ = 0.677;
    double mu_dynamic_ = 0.525;
    double critical_slip_distance_ = 0.40;
    double fault_slip_onset_time_ = 0.0;   // Time when prescribed slip begins (s)
    double fault_slip_rise_time_ = 0.0;    // Duration of slip ramp (s); 0 = instantaneous

    // Initial fault stress (for SCEC TPV-type benchmarks)
    bool fault_has_initial_stress_ = false;
    double fault_initial_normal_stress_ = 0.0;    // Pa (compressive negative)
    double fault_initial_shear_stress_ = 0.0;     // Pa
    double fault_nucleation_center_x_ = 0.0;      // m
    double fault_nucleation_center_z_ = 0.0;      // m
    double fault_nucleation_radius_ = 0.0;         // m
    double fault_nucleation_shear_stress_ = 0.0;   // Pa (overstressed patch)

    // Output directory and format for HDF5/VTK output
    std::string output_directory_ = "output";
    bool output_topology_written_ = false;

    // Explosion source FEM injection (Phase 5)
    PetscInt explosion_cell_ = -1;
    PetscErrorCode locateExplosionCell();
    PetscErrorCode addExplosionSourceToResidual(PetscReal t, Vec locF);
    
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
    static PetscErrorCode ViscoelasticPostStep(TS ts);
    
    // Helper functions
    PetscErrorCode setupStructuredGrid();
    PetscErrorCode setupUnstructuredGrid();
    PetscErrorCode setupGmshGrid();
    PetscErrorCode setupCoordinateSystem();
    PetscErrorCode createFieldsFromConfig();
    PetscErrorCode interpolateMaterialProperties();
    PetscErrorCode applyMaterialPropertiesToKernels();
    
    // Coordinate transformation helpers
    void transformToModelCoordinates(double& x, double& y, double& z) const;
    void transformToInputCoordinates(double& x, double& y, double& z) const;
};

} // namespace FSRM

#endif // SIMULATOR_HPP
