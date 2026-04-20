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

    // Auxiliary DM for heterogeneous material properties. Session 18 keeps
    // this as the material-only aux: three cell-centered fields (lambda,
    // mu, rho) attached at material_label_ values 1 (neg) and 2 (pos) so
    // the hybrid driver's bulk-side displacement callbacks can read them.
    DM  auxDM_ = nullptr;
    Vec auxVec_ = nullptr;

    // Session 18: slip-only aux DM/Vec attached at (interfaces_label_, 1, 0).
    // Built as a topological clone of the main DM with a single (dim-1)
    // surface FE restricted to interfaces_label_ and marked cohesive via
    // PetscDSSetCohesive; this makes DMGetDS(slipAuxDM_) return a DS whose
    // totDim matches the closure at cohesive cells and keeps PETSc 3.25's
    // hybrid-driver closure check at plexfem.c:3982 consistent. Session 18
    // attempted to fold this subfield into auxDM_ per PyLith's composite
    // Field pattern, but PETSc 3.25's DMGetDS returns probs[0] of the DM
    // while the interfaces-region DS is probs[1]; the two-DM split is the
    // practical analogue that keeps the lookups well-defined.
    DM  slipAuxDM_ = nullptr;
    Vec slipAuxVec_ = nullptr;
    PetscInt fault_slip_aux_field_idx_ = -1;

    // Session 21: rigid-body near-null-space on the displacement field,
    // sized for the full assembled Jacobian (displacement DOFs populated,
    // Lagrange DOFs zero). Built once in setupSolvers when aux-slip is
    // active and GAMG is selected; attached to J in FormJacobian after
    // the hybrid assembly bookend so PCGAMG can find it via
    // MatGetNearNullSpace. The sub-DM + PetscObjectCompose path used
    // through Session 20 did not reach PCGAMG on the full monolithic
    // matrix (diagnosed by the Session 21 probe: MatGetNearNullSpace J
    // returned NULL with only the composed-on-field path active).
    MatNullSpace rigidBodyNullSpace_ = nullptr;
    PetscErrorCode updateFaultAuxiliary(PetscReal t);
    PetscErrorCode setupAuxiliaryDM();
    PetscErrorCode populateAuxFieldsByDepth();
    PetscErrorCode populateAuxFieldsByMaterialLabel();
    PetscErrorCode populateAuxFieldsByVelocityModel();
    PetscErrorCode applyExplosionDamageToAuxFields();
    PetscErrorCode createCohesiveCellLabel();
    PetscErrorCode getOrCreateInterfacesLabel(DMLabel *interfacesLabel);
    PetscErrorCode getOrCreateInterfaceFacetsLabel(DMLabel *interfaceFacetsLabel);
    // Session 10: identify cohesive vertices that lie on the global domain
    // boundary (the rim of the fault). These need an essential Dirichlet BC
    // on the Lagrange field to pin the otherwise rank-deficient constraint
    // block in the saddle-point matrix. PyLith pattern, see
    // libsrc/pylith/faults/FaultCohesive.cc createConstraints (buried-edges).
    PetscErrorCode getOrCreateFaultBoundaryLabel(DMLabel *faultBoundaryLabel);

    // Session 12: PyLith "buried_cohesive" label. Contains cohesive
    // edge-prisms (DM_POLYTOPE_SEG_PRISM_TENSOR) and cohesive vertex-prisms
    // (DM_POLYTOPE_POINT_PRISM_TENSOR) whose cone touches the domain bounding
    // box (fault rim). These are the points that carry Lagrange DOFs in
    // PETSc's hybrid numbering (Session 11 found d1=25 for the 4x4x4 test).
    // Pinning them with a zero Dirichlet on the Lagrange field removes the
    // structural rank deficiency in the saddle-point constraint block.
    // PyLith reference: libsrc/pylith/faults/FaultCohesive.cc:444-480.
    PetscErrorCode getOrCreateBuriedCohesiveLabel(DMLabel *buriedCohesiveLabel);

    // Session 4: material label that tags regular cells adjacent to cohesive
    // prisms with value 1 (negative side) or value 2 (positive side). This is
    // the PETSc hybrid-driver pattern (see dm/impls/plex/tests/ex5.c
    // TestAssembly) required by DMPlexComputeResidualHybridByKey /
    // DMPlexComputeJacobianHybridByKey to dispatch the per-side displacement
    // weak forms.
    PetscErrorCode createMaterialLabel();

    // PyLith-style "cohesive interface" label: cohesive cells plus their
    // closure faces, edges, and vertices on the negative side of the fault.
    // Populated by getOrCreateInterfacesLabel during setupFaultNetwork and
    // consumed by setupBoundaryConditions to register the Lagrange
    // fault_constraint Neumann BC on the cohesive geometry.
    DMLabel interfaces_label_ = nullptr;

    // Facets-only ("cohesive interface facets") label: height-1 subset of
    // interfaces_label_. DMAddBoundary dispatches BdResidual on every labeled
    // point, so passing the full-stratum interfaces_label_ feeds cohesive
    // cells and vertices into the kernel and yields NaN. This label contains
    // only the cohesive facets that are the correct target for BdResidual.
    DMLabel interface_facets_label_ = nullptr;

    // Session 10: subset of cohesive vertices that lie on the global domain
    // boundary (the fault rim). Used to pin Lagrange DOFs via post-assembly
    // MatZeroRowsColumnsLocal in FormJacobian and an explicit zero of the
    // residual entries in FormFunction. Removes the structural rank
    // deficiency in the constraint block diagnosed in Session 9 (PCSVD
    // smallest sigma ~3.7e-08, condition number 5.36e+17).
    DMLabel fault_boundary_label_ = nullptr;
    // Session 12: label built by getOrCreateBuriedCohesiveLabel containing
    // cohesive edge/vertex prisms at the fault rim. Used to register a
    // PetscDSAddBoundary essential BC on the region DS that owns the
    // Lagrange field.
    DMLabel buried_cohesive_label_ = nullptr;
    // Cached LOCAL section indices of Lagrange DOFs at fault-boundary
    // vertices, populated at the end of setupFields once the local section
    // exists. Empty when faults are disabled or no rim vertices were found.
    std::vector<PetscInt> fault_boundary_lagrange_local_dofs_;
    PetscErrorCode cacheFaultBoundaryLagrangeDofs();
    PetscErrorCode pinFaultBoundaryLagrangeJacobian(Mat J, Mat P);
    PetscErrorCode pinFaultBoundaryLagrangeResidual(Vec locF);

    // Session 4: material label with value 1 on regular cells on the negative
    // side of the fault and value 2 on regular cells on the positive side.
    // Used as key[0].label / key[1].label in DMPlexComputeResidualHybridByKey.
    DMLabel material_label_ = nullptr;

    // Field indices cached during setupFields so FormFunction / FormJacobian
    // can build PetscFormKey arrays without re-querying the DS on every call.
    // -1 means the field is not present.
    PetscInt displacement_field_idx_ = -1;
    PetscInt lagrange_field_idx_ = -1;
    
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
    [[deprecated("Session 6: hybrid Jacobian (DMPlexComputeJacobianHybridByKey) "
                 "is the sole driver of cohesive coupling entries. Calling this "
                 "augmentation re-introduces the double-write that caused "
                 "DIVERGED_LINEAR_SOLVE iterations 0 in Session 5.")]]
    PetscErrorCode addCohesivePenaltyToJacobian(Mat J, Vec locU = nullptr);

    /**
     * @brief Zero the Lagrange residual at interior (non-cohesive) DOFs.
     *
     * Per docs/PYLITH_COMPATIBILITY.md Section B Option B, interior
     * Lagrange degrees of freedom are not driven by the cohesive
     * constraint and must not retain any residual contribution. PETSc
     * 3.25 requires the weak volume callback (f0_weak_lagrange) to stay
     * registered for the BdResidual on cohesive cells to fire reliably;
     * this routine overwrites the small epsilon = 1e-4 contribution
     * that callback leaves at non-cohesive Lagrange DOFs so the manual
     * diagonal in addCohesivePenaltyToJacobian is the sole driver
     * there. Cohesive vertex contributions are intentionally preserved.
     *
     * Called from FormFunction after DMPlexTSComputeIFunctionFEM and the
     * other manual residual additions.
     *
     * @param locF Local residual vector to update in place.
     * @return PETSC_SUCCESS or a PETSc error code.
     */
    PetscErrorCode addInteriorLagrangeResidual(Vec locF);

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
