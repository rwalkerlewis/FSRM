# Machine Learning Capabilities for FSRM

This document describes the comprehensive AI/ML enhancements implemented in FSRM (Fully-coupled Subsurface Reservoir Model) to accelerate simulations, improve accuracy, and enable new capabilities.

## Overview

The ML system provides **10 major enhancement categories** that can be used independently or in combination:

| Category | Purpose | Expected Speedup |
|----------|---------|------------------|
| Neural Surrogates | Replace expensive physics calculations | 10-100x |
| Graph Neural Operators | Learn on unstructured meshes | 5-50x |
| Reduced Order Models | Accelerate parametric studies | 100-1000x |
| Neural Linear Algebra | Accelerate linear solvers | 2-10x |
| Bayesian UQ | Quantify prediction uncertainty | N/A (new capability) |
| Neural Inversion | Accelerate history matching | 10-100x |
| Multi-Fidelity Learning | Combine cheap/expensive models | 10-50x |
| Neural Time Stepping | Optimal time step selection | 2-5x |
| Neural AMR | Intelligent mesh refinement | 2-5x |
| Physics-Informed FNO | Guarantee physics constraints | 5-20x |

## Model Selection

Each ML model can be selected **in place of** or **in addition to** existing methods. The `MLModelRegistry` provides a unified interface:

```cpp
// Access registry
auto& registry = FSRM::ML::MLModelRegistry::getInstance();

// Register a model
registry.registerFNO("pressure_solver", fno_config);

// Get and use model
auto* fno = registry.getModel<FNOModel>("pressure_solver");
Tensor output = fno->forward(input);

// Or use wrapper for seamless switching
MLSolverWrapper wrapper(ModelCategory::FNO);
wrapper.setActiveModel("pressure_solver");
wrapper.enableML(true);  // Enable ML, or false for classical
Tensor result = wrapper.solve(input);
```

## Component Details

### 1. Neural Surrogates (`NeuralSurrogates.hpp`)

Replace expensive physics computations with fast neural network inference:

#### Neural PVT Model
```cpp
NeuralPVTConfig config;
config.hidden_layers = {128, 256, 256, 128};
config.properties = {PVTProperty::OIL_FVF, PVTProperty::OIL_VISCOSITY};

NeuralPVTModel pvt(config);
pvt.train(inputs, outputs, 200, 1e-3);

double Bo = pvt.getOilFVF(pressure, temperature, Rs);
```

#### Neural Flash Calculator
```cpp
NeuralFlashCalculator flash(config);
double L;
std::vector<double> x, y;
flash.flash(z_composition, P, T, L, x, y);
```

#### Neural Relative Permeability
```cpp
NeuralRelPermModel relperm(config);
double krw, kro, krg;
relperm.getRelPerm(Sw, Sg, krw, kro, krg);
```

### 2. Graph Neural Operators (`GraphNeuralOperator.hpp`)

For unstructured meshes with complex geometries:

```cpp
GNOConfig config;
config.input_dim = 3;
config.output_dim = 3;
config.hidden_dim = 128;
config.num_layers = 6;
config.layer_type = GNOConfig::LayerType::GAT;

GraphNeuralOperator gno(config);

// Build graph from mesh
MeshGraph graph;
graph.fromMesh(dm);
graph.setNodeFeatures(solution);
graph.computeEdgeFeatures();

// Predict
Tensor output = gno.predict(graph);
```

**Available GNO Types:**
- `BASIC_GNO`: Standard message-passing
- `MESH_GRAPH_NET`: For physics simulation
- `MULTIPOLE_GNN`: Multi-scale learning
- `PHYSICS_INFORMED_GNN`: With conservation laws

### 3. Reduced Order Models (`NeuralReducedOrderModel.hpp`)

For fast parametric studies:

```cpp
PODNNConfig config;
config.num_modes = 50;
config.energy_threshold = 0.999;
config.hidden_layers = {256, 256, 256};

PODNeuralROM rom(config);
rom.buildFromSnapshots(snapshots, times);

Tensor solution = rom.solve(initial_condition, t_final, dt);
```

**ROM Types:**
- `POD_NN`: POD + neural dynamics
- `POD_LSTM`: POD + LSTM for temporal
- `AUTOENCODER`: Nonlinear reduction
- `DEEPONET`: Deep operator network

### 4. Neural Linear Algebra (`NeuralLinearAlgebra.hpp`)

Accelerate linear solvers:

```cpp
NeuralCoarseGridConfig config;
config.fine_dim = 10000;
config.coarse_dim = 100;

NeuralCoarseGridCorrection cgc(config);
cgc.setup(A_fine);
cgc.apply(residual, correction);

// Or use full neural solver
NeuralLinearSolver::SolverConfig solver_config;
solver_config.solver_type = NeuralLinearSolver::SolverType::HYBRID;

NeuralLinearSolver solver(solver_config);
solver.solve(A, b, x);
```

### 5. Uncertainty Quantification (`BayesianNeuralOperator.hpp`)

Quantify prediction uncertainty:

```cpp
DeepEnsembleConfig config;
config.num_models = 5;
config.layer_sizes = {128, 256, 128};

DeepEnsemble ensemble(config);
ensemble.train(inputs, targets);

Tensor mean, variance;
ensemble.predict(input, mean, variance);
```

**UQ Methods:**
- Bayesian Neural Networks
- Deep Ensembles
- MC Dropout
- Probabilistic FNO

### 6. Neural Inversion (`NeuralInversion.hpp`)

Accelerate history matching:

```cpp
NeuralEnKF::EnKFConfig config;
config.ensemble_size = 100;
config.use_neural_forward = true;
config.use_neural_localization = true;

NeuralEnKF enkf(config);
enkf.initialize(prior_mean, prior_covariance);

for (auto& obs : observations) {
    enkf.update(obs, obs_covariance, forward_model);
}

Tensor posterior_mean = enkf.getMean();
```

### 7. Multi-Fidelity Learning (`MultiFidelityLearning.hpp`)

Combine cheap and expensive models:

```cpp
AdditiveCorrection::AdditiveConfig config;
config.hidden_layers = {256, 256};
config.include_lf_output = true;

AdditiveCorrection correction(config);
correction.train(inputs, lf_outputs, hf_outputs);

Tensor hf_prediction = correction.correct(lf_output, input);
```

**Correction Methods:**
- Additive: `y_HF = y_LF + Δ`
- Multiplicative: `y_HF = y_LF * (1 + Δ)`
- Residual: Cascaded corrections
- Transfer Learning: Pre-train on LF

### 8. Neural Time Stepping (`NeuralTimestepping.hpp`)

Optimal time step selection:

```cpp
NeuralTimeStepPredictor::TimeStepConfig config;
config.state_history_length = 5;
config.use_physics_features = true;

NeuralTimeStepPredictor predictor(config);
double dt = predictor.predictTimeStep(state, time, dt_prev);

// Or full controller
NeuralTimeSteppingController controller(controller_config);
controller.step(state, residual_func, dt, mode);
```

### 9. Neural AMR (`NeuralAMR.hpp`)

Intelligent mesh refinement:

```cpp
NeuralAMRController::AMRConfig config;
config.use_predictive = true;
config.use_physics_aware = true;
config.use_feature_tracking = true;

NeuralAMRController amr(config);

if (amr.shouldAdapt(step, time)) {
    std::vector<PetscInt> refine, coarsen;
    amr.getRefinementFlags(dm, solution, time, refine, coarsen);
    amr.adapt(dm, solution, time);
}
```

### 10. Physics-Informed FNO (`FourierNeuralOperator.hpp`)

FNO with physics guarantees:

```cpp
PINNConfig config;
config.physics_type = PINNConfig::PhysicsType::WAVE;
config.physics_loss_weight = 0.1;
config.enforce_energy_conservation = true;

PhysicsInformedFNO pinn(config);
pinn.train(labeled_data);

// Verify physics
bool valid = pinn.checkPhysicsConstraints(prediction, 1e-3);
```

## Configuration Files

Configuration files in `config/` directory:

- `ml_neural_surrogates.config`: PVT, Flash, RelPerm surrogates
- `ml_neural_operators.config`: FNO, GNO, ROM configurations
- `ml_advanced.config`: UQ, inversion, time stepping, AMR

Example usage in main config:

```ini
[ml]
enabled = true
model_directory = ml_models
use_gpu = true

# Enable specific models
fno_enabled = true
fno_solver_method = HYBRID_FNO_REFINE

neural_pvt_enabled = true
neural_flash_enabled = true

use_neural_timestep = true
use_neural_amr = true
```

## Best Practices

### Model Selection Guidelines

| Problem Type | Recommended Model |
|--------------|-------------------|
| Regular grid PDEs | FNO |
| Unstructured mesh | GNO |
| Parametric studies | POD-ROM |
| Multi-phase flow | Neural PVT + Flash |
| History matching | Neural EnKF |
| Uncertainty needed | Deep Ensemble |
| Complex geometry | Physics-Informed GNN |

### Training Data Requirements

| Model | Minimum Samples | Recommended |
|-------|----------------|-------------|
| Neural Surrogates | 1,000 | 10,000 |
| FNO | 500 | 5,000 |
| GNO | 100 | 1,000 |
| POD-ROM | 50 | 200 |
| Deep Ensemble | 1,000 | 10,000 |

### Performance Tips

1. **Start with HYBRID mode**: Use `HYBRID_FNO_REFINE` to combine ML speed with numerical accuracy
2. **Use GPU**: Enable GPU acceleration for 5-50x speedup
3. **Pre-train offline**: Train models offline on representative data
4. **Enable online learning**: Adapt during simulation when conditions change
5. **Monitor uncertainty**: Use ensemble methods to detect out-of-distribution inputs

## File Structure

```
include/
├── FourierNeuralOperator.hpp      # Core FNO + PINN extensions
├── NeuralSurrogates.hpp           # PVT, Flash, RelPerm, Preconditioner
├── GraphNeuralOperator.hpp        # GNO, MeshGraphNet, MultipoleGNN
├── NeuralReducedOrderModel.hpp    # POD-ROM, LSTM-ROM, Autoencoder
├── NeuralLinearAlgebra.hpp        # Coarse grid, AMG, Jacobian
├── BayesianNeuralOperator.hpp     # BNN, Ensemble, MC Dropout
├── NeuralInversion.hpp            # EnKF, Physics-guided inversion
├── MultiFidelityLearning.hpp      # Additive/multiplicative correction
├── NeuralTimestepping.hpp         # Time step, IMEX, stability
├── NeuralAMR.hpp                  # Error indicator, predictive refinement
└── MLModelRegistry.hpp            # Unified model management

src/
├── FourierNeuralOperator.cpp
├── NeuralSurrogates.cpp
├── GraphNeuralOperator.cpp
├── NeuralReducedOrderModel.cpp
└── MLModelRegistry.cpp

config/
├── ml_neural_surrogates.config
├── ml_neural_operators.config
└── ml_advanced.config
```

## Future Extensions

The modular architecture supports easy addition of:

- TensorRT/ONNX export for production deployment
- Distributed training across multiple GPUs
- AutoML for hyperparameter optimization
- Neural Architecture Search (NAS)
- Federated learning for multi-site data
- Reinforcement learning for control problems

## References

1. Li et al. (2020). "Fourier Neural Operator for Parametric PDEs"
2. Pfaff et al. (2020). "Learning Mesh-Based Simulation with Graph Networks"
3. Raissi et al. (2019). "Physics-Informed Neural Networks"
4. Lakshminarayanan et al. (2017). "Simple and Scalable Predictive Uncertainty Estimation using Deep Ensembles"
5. Hesthaven et al. (2018). "Non-intrusive reduced order modeling of nonlinear problems using neural networks"
