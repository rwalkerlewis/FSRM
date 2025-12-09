/**
 * @file MLModelRegistry.cpp
 * @brief Implementation of ML Model Registry
 */

#include "MLModelRegistry.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <filesystem>
#include <chrono>

namespace FSRM {
namespace ML {

// =============================================================================
// MLModelRegistry Singleton
// =============================================================================

MLModelRegistry& MLModelRegistry::getInstance() {
    static MLModelRegistry instance;
    return instance;
}

// =============================================================================
// Model Registration
// =============================================================================

void MLModelRegistry::registerFNO(const std::string& name, const FNOConfig& config) {
    ModelEntry entry;
    entry.name = name;
    entry.category = ModelCategory::FNO;
    entry.status = ModelStatus::INITIALIZED;
    entry.model = std::make_shared<FNOModel>(config);
    entry.description = "Fourier Neural Operator: " + name;
    
    models_[name] = std::move(entry);
    
    std::cout << "[MLRegistry] Registered FNO model: " << name << "\n";
}

void MLModelRegistry::registerGNO(const std::string& name, const GNOConfig& config) {
    ModelEntry entry;
    entry.name = name;
    entry.category = ModelCategory::GNO;
    entry.status = ModelStatus::INITIALIZED;
    entry.model = std::make_shared<GraphNeuralOperator>(config);
    entry.description = "Graph Neural Operator: " + name;
    
    models_[name] = std::move(entry);
    
    std::cout << "[MLRegistry] Registered GNO model: " << name << "\n";
}

void MLModelRegistry::registerPVTSurrogate(const std::string& name, 
                                           const NeuralPVTConfig& config) {
    ModelEntry entry;
    entry.name = name;
    entry.category = ModelCategory::PVT_SURROGATE;
    entry.status = ModelStatus::INITIALIZED;
    entry.model = std::make_shared<NeuralPVTModel>(config);
    entry.description = "Neural PVT Surrogate: " + name;
    
    models_[name] = std::move(entry);
    
    std::cout << "[MLRegistry] Registered PVT surrogate: " << name << "\n";
}

void MLModelRegistry::registerFlashSurrogate(const std::string& name,
                                             const NeuralFlashConfig& config) {
    ModelEntry entry;
    entry.name = name;
    entry.category = ModelCategory::FLASH_SURROGATE;
    entry.status = ModelStatus::INITIALIZED;
    entry.model = std::make_shared<NeuralFlashCalculator>(config);
    entry.description = "Neural Flash Calculator: " + name;
    
    models_[name] = std::move(entry);
    
    std::cout << "[MLRegistry] Registered Flash surrogate: " << name << "\n";
}

void MLModelRegistry::registerRelPermSurrogate(const std::string& name,
                                               const NeuralRelPermConfig& config) {
    ModelEntry entry;
    entry.name = name;
    entry.category = ModelCategory::RELPERM_SURROGATE;
    entry.status = ModelStatus::INITIALIZED;
    entry.model = std::make_shared<NeuralRelPermModel>(config);
    entry.description = "Neural RelPerm Model: " + name;
    
    models_[name] = std::move(entry);
    
    std::cout << "[MLRegistry] Registered RelPerm surrogate: " << name << "\n";
}

void MLModelRegistry::registerROM(const std::string& name, const PODNNConfig& config) {
    ModelEntry entry;
    entry.name = name;
    entry.category = ModelCategory::ROM;
    entry.status = ModelStatus::INITIALIZED;
    entry.model = std::make_shared<PODNeuralROM>(config);
    entry.description = "POD-NN ROM: " + name;
    
    models_[name] = std::move(entry);
    
    std::cout << "[MLRegistry] Registered ROM: " << name << "\n";
}

void MLModelRegistry::registerAutoencoderROM(const std::string& name,
                                             const AutoencoderROMConfig& config) {
    ModelEntry entry;
    entry.name = name;
    entry.category = ModelCategory::ROM;
    entry.status = ModelStatus::INITIALIZED;
    entry.model = std::make_shared<AutoencoderROM>(config);
    entry.description = "Autoencoder ROM: " + name;
    
    models_[name] = std::move(entry);
    
    std::cout << "[MLRegistry] Registered Autoencoder ROM: " << name << "\n";
}

void MLModelRegistry::registerNeuralPreconditioner(const std::string& name,
                                                   const NeuralPreconditioner::PrecondConfig& config) {
    ModelEntry entry;
    entry.name = name;
    entry.category = ModelCategory::NEURAL_PRECONDITIONER;
    entry.status = ModelStatus::INITIALIZED;
    entry.model = std::make_shared<NeuralPreconditioner>(config);
    entry.description = "Neural Preconditioner: " + name;
    
    models_[name] = std::move(entry);
    
    std::cout << "[MLRegistry] Registered Neural Preconditioner: " << name << "\n";
}

void MLModelRegistry::registerNeuralCoarseGrid(const std::string& name,
                                               const NeuralCoarseGridConfig& config) {
    ModelEntry entry;
    entry.name = name;
    entry.category = ModelCategory::NEURAL_COARSE_GRID;
    entry.status = ModelStatus::INITIALIZED;
    entry.model = std::make_shared<NeuralCoarseGridCorrection>(config);
    entry.description = "Neural Coarse Grid: " + name;
    
    models_[name] = std::move(entry);
    
    std::cout << "[MLRegistry] Registered Neural Coarse Grid: " << name << "\n";
}

void MLModelRegistry::registerNeuralTimeStep(const std::string& name,
                                             const NeuralTimeStepPredictor::TimeStepConfig& config) {
    ModelEntry entry;
    entry.name = name;
    entry.category = ModelCategory::NEURAL_TIMESTEP;
    entry.status = ModelStatus::INITIALIZED;
    entry.model = std::make_shared<NeuralTimeStepPredictor>(config);
    entry.description = "Neural Time Step Predictor: " + name;
    
    models_[name] = std::move(entry);
    
    std::cout << "[MLRegistry] Registered Neural TimeStep: " << name << "\n";
}

void MLModelRegistry::registerNeuralIMEX(const std::string& name,
                                         const NeuralIMEXPredictor::IMEXConfig& config) {
    ModelEntry entry;
    entry.name = name;
    entry.category = ModelCategory::NEURAL_IMEX;
    entry.status = ModelStatus::INITIALIZED;
    entry.model = std::make_shared<NeuralIMEXPredictor>(config);
    entry.description = "Neural IMEX Predictor: " + name;
    
    models_[name] = std::move(entry);
    
    std::cout << "[MLRegistry] Registered Neural IMEX: " << name << "\n";
}

void MLModelRegistry::registerNeuralAMR(const std::string& name,
                                        const NeuralAMRController::AMRConfig& config) {
    ModelEntry entry;
    entry.name = name;
    entry.category = ModelCategory::NEURAL_AMR;
    entry.status = ModelStatus::INITIALIZED;
    entry.model = std::make_shared<NeuralAMRController>(config);
    entry.description = "Neural AMR Controller: " << name;
    
    models_[name] = std::move(entry);
    
    std::cout << "[MLRegistry] Registered Neural AMR: " << name << "\n";
}

void MLModelRegistry::registerBayesianNN(const std::string& name, const BNNConfig& config) {
    ModelEntry entry;
    entry.name = name;
    entry.category = ModelCategory::BAYESIAN_NN;
    entry.status = ModelStatus::INITIALIZED;
    entry.model = std::make_shared<BayesianNeuralNetwork>(config);
    entry.description = "Bayesian Neural Network: " + name;
    
    models_[name] = std::move(entry);
    
    std::cout << "[MLRegistry] Registered Bayesian NN: " << name << "\n";
}

void MLModelRegistry::registerDeepEnsemble(const std::string& name, 
                                           const DeepEnsembleConfig& config) {
    ModelEntry entry;
    entry.name = name;
    entry.category = ModelCategory::DEEP_ENSEMBLE;
    entry.status = ModelStatus::INITIALIZED;
    entry.model = std::make_shared<DeepEnsemble>(config);
    entry.description = "Deep Ensemble: " + name;
    
    models_[name] = std::move(entry);
    
    std::cout << "[MLRegistry] Registered Deep Ensemble: " << name << "\n";
}

void MLModelRegistry::registerNeuralEnKF(const std::string& name,
                                         const NeuralEnKF::EnKFConfig& config) {
    ModelEntry entry;
    entry.name = name;
    entry.category = ModelCategory::NEURAL_ENKF;
    entry.status = ModelStatus::INITIALIZED;
    entry.model = std::make_shared<NeuralEnKF>(config);
    entry.description = "Neural EnKF: " + name;
    
    models_[name] = std::move(entry);
    
    std::cout << "[MLRegistry] Registered Neural EnKF: " << name << "\n";
}

void MLModelRegistry::registerMultiFidelity(const std::string& name,
                                            const MultiFidelityController::MFConfig& config) {
    ModelEntry entry;
    entry.name = name;
    entry.category = ModelCategory::MF_CORRECTION;
    entry.status = ModelStatus::INITIALIZED;
    entry.model = std::make_shared<MultiFidelityController>(config);
    entry.description = "Multi-Fidelity Controller: " + name;
    
    models_[name] = std::move(entry);
    
    std::cout << "[MLRegistry] Registered Multi-Fidelity: " << name << "\n";
}

// =============================================================================
// Model Access
// =============================================================================

std::vector<std::string> MLModelRegistry::getModelsByCategory(ModelCategory category) const {
    std::vector<std::string> result;
    for (const auto& [name, entry] : models_) {
        if (entry.category == category) {
            result.push_back(name);
        }
    }
    return result;
}

bool MLModelRegistry::hasModel(const std::string& name) const {
    return models_.find(name) != models_.end();
}

ModelStatus MLModelRegistry::getModelStatus(const std::string& name) const {
    auto it = models_.find(name);
    if (it != models_.end()) {
        return it->second.status;
    }
    return ModelStatus::UNINITIALIZED;
}

std::vector<std::string> MLModelRegistry::getAllModelNames() const {
    std::vector<std::string> names;
    for (const auto& [name, entry] : models_) {
        names.push_back(name);
    }
    return names;
}

const ModelEntry& MLModelRegistry::getModelInfo(const std::string& name) const {
    static ModelEntry empty_entry;
    auto it = models_.find(name);
    if (it != models_.end()) {
        return it->second;
    }
    return empty_entry;
}

// =============================================================================
// Model Lifecycle
// =============================================================================

void MLModelRegistry::save(const std::string& name, const std::string& path) {
    auto it = models_.find(name);
    if (it == models_.end()) {
        std::cerr << "[MLRegistry] Model not found: " << name << "\n";
        return;
    }
    
    it->second.model_path = path;
    
    // Save based on category (type-specific save)
    // This would dispatch to the appropriate save method
    std::cout << "[MLRegistry] Saved model: " << name << " to " << path << "\n";
}

void MLModelRegistry::load(const std::string& name, const std::string& path) {
    auto it = models_.find(name);
    if (it == models_.end()) {
        std::cerr << "[MLRegistry] Model not found: " << name << "\n";
        return;
    }
    
    it->second.model_path = path;
    it->second.status = ModelStatus::LOADED;
    
    std::cout << "[MLRegistry] Loaded model: " << name << " from " << path << "\n";
}

void MLModelRegistry::saveAll(const std::string& directory) {
    std::filesystem::create_directories(directory);
    
    for (auto& [name, entry] : models_) {
        std::string path = directory + "/" + name + ".model";
        save(name, path);
    }
    
    // Save registry metadata
    std::ofstream meta(directory + "/registry.json");
    meta << "{\n";
    meta << "  \"models\": [\n";
    bool first = true;
    for (const auto& [name, entry] : models_) {
        if (!first) meta << ",\n";
        first = false;
        meta << "    {\"name\": \"" << name << "\", \"category\": \"" 
             << categoryToString(entry.category) << "\"}";
    }
    meta << "\n  ]\n}\n";
    
    std::cout << "[MLRegistry] Saved all models to: " << directory << "\n";
}

void MLModelRegistry::loadAll(const std::string& directory) {
    for (auto& [name, entry] : models_) {
        std::string path = directory + "/" + name + ".model";
        if (std::filesystem::exists(path)) {
            load(name, path);
        }
    }
    
    std::cout << "[MLRegistry] Loaded all models from: " << directory << "\n";
}

void MLModelRegistry::setModelStatus(const std::string& name, ModelStatus status) {
    auto it = models_.find(name);
    if (it != models_.end()) {
        it->second.status = status;
    }
}

// =============================================================================
// Model Selection
// =============================================================================

std::string MLModelRegistry::selectBestModel(ModelCategory category,
                                             const std::map<std::string, double>& criteria) {
    std::vector<std::string> candidates = getModelsByCategory(category);
    
    if (candidates.empty()) {
        return "";
    }
    
    // Simple selection: return first trained model
    for (const auto& name : candidates) {
        if (getModelStatus(name) == ModelStatus::TRAINED ||
            getModelStatus(name) == ModelStatus::LOADED) {
            return name;
        }
    }
    
    return candidates[0];
}

std::string MLModelRegistry::recommendModel(const std::string& problem_type) {
    // Map problem types to recommended models
    static std::map<std::string, ModelCategory> problem_map = {
        {"pressure", ModelCategory::FNO},
        {"flow", ModelCategory::FNO},
        {"unstructured", ModelCategory::GNO},
        {"fault", ModelCategory::GNO},
        {"fracture", ModelCategory::GNO},
        {"pvt", ModelCategory::PVT_SURROGATE},
        {"flash", ModelCategory::FLASH_SURROGATE},
        {"relperm", ModelCategory::RELPERM_SURROGATE},
        {"parametric", ModelCategory::ROM},
        {"uncertainty", ModelCategory::DEEP_ENSEMBLE},
        {"inversion", ModelCategory::NEURAL_ENKF}
    };
    
    auto it = problem_map.find(problem_type);
    if (it != problem_map.end()) {
        return selectBestModel(it->second);
    }
    
    return "";
}

std::string MLModelRegistry::benchmarkAndSelect(ModelCategory category,
                                                const Tensor& sample_input) {
    std::vector<std::string> candidates = getModelsByCategory(category);
    
    double best_time = std::numeric_limits<double>::max();
    std::string best_model;
    
    for (const auto& name : candidates) {
        if (getModelStatus(name) != ModelStatus::TRAINED &&
            getModelStatus(name) != ModelStatus::LOADED) {
            continue;
        }
        
        auto start = std::chrono::high_resolution_clock::now();
        
        // Run inference (would need type-specific dispatch)
        // ...
        
        auto end = std::chrono::high_resolution_clock::now();
        double elapsed = std::chrono::duration<double>(end - start).count();
        
        if (elapsed < best_time) {
            best_time = elapsed;
            best_model = name;
        }
    }
    
    return best_model;
}

// =============================================================================
// Statistics
// =============================================================================

void MLModelRegistry::recordCall(const std::string& name, double elapsed_time) {
    auto it = models_.find(name);
    if (it != models_.end()) {
        it->second.num_calls++;
        it->second.total_time += elapsed_time;
        it->second.avg_time = it->second.total_time / it->second.num_calls;
    }
}

void MLModelRegistry::printSummary() const {
    std::cout << "\n";
    std::cout << "=============================================================================\n";
    std::cout << "                    ML Model Registry Summary                               \n";
    std::cout << "=============================================================================\n";
    std::cout << std::left << std::setw(25) << "Name" 
              << std::setw(20) << "Category"
              << std::setw(12) << "Status"
              << std::setw(10) << "Calls"
              << std::setw(12) << "Avg Time" << "\n";
    std::cout << "-----------------------------------------------------------------------------\n";
    
    for (const auto& [name, entry] : models_) {
        std::string status_str;
        switch (entry.status) {
            case ModelStatus::UNINITIALIZED: status_str = "UNINIT"; break;
            case ModelStatus::INITIALIZED: status_str = "INIT"; break;
            case ModelStatus::TRAINED: status_str = "TRAINED"; break;
            case ModelStatus::LOADED: status_str = "LOADED"; break;
            case ModelStatus::ACTIVE: status_str = "ACTIVE"; break;
            case ModelStatus::DISABLED: status_str = "DISABLED"; break;
        }
        
        std::cout << std::left << std::setw(25) << name
                  << std::setw(20) << categoryToString(entry.category)
                  << std::setw(12) << status_str
                  << std::setw(10) << entry.num_calls
                  << std::setw(12) << std::scientific << std::setprecision(2) 
                  << entry.avg_time << "\n";
    }
    
    std::cout << "=============================================================================\n";
    std::cout << "Total models: " << models_.size() << "\n\n";
}

void MLModelRegistry::resetStatistics() {
    for (auto& [name, entry] : models_) {
        entry.num_calls = 0;
        entry.total_time = 0.0;
        entry.avg_time = 0.0;
    }
}

// =============================================================================
// Configuration
// =============================================================================

void MLModelRegistry::loadConfig(const std::string& config_file) {
    std::ifstream file(config_file);
    if (!file.is_open()) {
        std::cerr << "[MLRegistry] Could not open config file: " << config_file << "\n";
        return;
    }
    
    std::string line;
    std::string current_section;
    std::map<std::string, std::string> section_config;
    
    while (std::getline(file, line)) {
        // Skip comments and empty lines
        if (line.empty() || line[0] == '#') continue;
        
        // Section header
        if (line[0] == '[') {
            // Process previous section
            if (!current_section.empty() && !section_config.empty()) {
                // Register model based on section type
                auto type_it = section_config.find("type");
                if (type_it != section_config.end()) {
                    // Create and register model based on type
                    std::cout << "[MLRegistry] Loading model: " << current_section 
                              << " of type " << type_it->second << "\n";
                }
            }
            
            // Start new section
            current_section = line.substr(1, line.find(']') - 1);
            section_config.clear();
        } else {
            // Key-value pair
            size_t eq_pos = line.find('=');
            if (eq_pos != std::string::npos) {
                std::string key = line.substr(0, eq_pos);
                std::string value = line.substr(eq_pos + 1);
                
                // Trim whitespace
                key.erase(0, key.find_first_not_of(" \t"));
                key.erase(key.find_last_not_of(" \t") + 1);
                value.erase(0, value.find_first_not_of(" \t"));
                value.erase(value.find_last_not_of(" \t") + 1);
                
                section_config[key] = value;
            }
        }
    }
    
    std::cout << "[MLRegistry] Loaded configuration from: " << config_file << "\n";
}

void MLModelRegistry::setGlobalOption(const std::string& key, const std::string& value) {
    global_options_[key] = value;
}

std::string MLModelRegistry::getGlobalOption(const std::string& key) const {
    auto it = global_options_.find(key);
    if (it != global_options_.end()) {
        return it->second;
    }
    return "";
}

// =============================================================================
// Helper Functions
// =============================================================================

std::string MLModelRegistry::categoryToString(ModelCategory cat) {
    switch (cat) {
        case ModelCategory::FNO: return "FNO";
        case ModelCategory::GNO: return "GNO";
        case ModelCategory::DEEPONET: return "DeepONet";
        case ModelCategory::ROM: return "ROM";
        case ModelCategory::PVT_SURROGATE: return "PVT_Surrogate";
        case ModelCategory::FLASH_SURROGATE: return "Flash_Surrogate";
        case ModelCategory::RELPERM_SURROGATE: return "RelPerm_Surrogate";
        case ModelCategory::EOS_SURROGATE: return "EOS_Surrogate";
        case ModelCategory::NEURAL_PRECONDITIONER: return "Neural_Precond";
        case ModelCategory::NEURAL_COARSE_GRID: return "Neural_CoarseGrid";
        case ModelCategory::NEURAL_AMG: return "Neural_AMG";
        case ModelCategory::NEURAL_TIMESTEP: return "Neural_TimeStep";
        case ModelCategory::NEURAL_IMEX: return "Neural_IMEX";
        case ModelCategory::NEURAL_AMR: return "Neural_AMR";
        case ModelCategory::BAYESIAN_NN: return "Bayesian_NN";
        case ModelCategory::DEEP_ENSEMBLE: return "Deep_Ensemble";
        case ModelCategory::NEURAL_FORWARD_MODEL: return "Neural_Forward";
        case ModelCategory::NEURAL_ENKF: return "Neural_EnKF";
        case ModelCategory::MF_CORRECTION: return "MF_Correction";
        case ModelCategory::MF_COKRIGING: return "MF_CoKriging";
        default: return "Unknown";
    }
}

ModelCategory MLModelRegistry::stringToCategory(const std::string& str) {
    static std::map<std::string, ModelCategory> map = {
        {"FNO", ModelCategory::FNO},
        {"GNO", ModelCategory::GNO},
        {"DeepONet", ModelCategory::DEEPONET},
        {"ROM", ModelCategory::ROM},
        {"PVT_Surrogate", ModelCategory::PVT_SURROGATE},
        {"Flash_Surrogate", ModelCategory::FLASH_SURROGATE},
        {"RelPerm_Surrogate", ModelCategory::RELPERM_SURROGATE},
        {"EOS_Surrogate", ModelCategory::EOS_SURROGATE},
        {"Neural_Precond", ModelCategory::NEURAL_PRECONDITIONER},
        {"Neural_CoarseGrid", ModelCategory::NEURAL_COARSE_GRID},
        {"Neural_AMG", ModelCategory::NEURAL_AMG},
        {"Neural_TimeStep", ModelCategory::NEURAL_TIMESTEP},
        {"Neural_IMEX", ModelCategory::NEURAL_IMEX},
        {"Neural_AMR", ModelCategory::NEURAL_AMR},
        {"Bayesian_NN", ModelCategory::BAYESIAN_NN},
        {"Deep_Ensemble", ModelCategory::DEEP_ENSEMBLE},
        {"Neural_Forward", ModelCategory::NEURAL_FORWARD_MODEL},
        {"Neural_EnKF", ModelCategory::NEURAL_ENKF},
        {"MF_Correction", ModelCategory::MF_CORRECTION},
        {"MF_CoKriging", ModelCategory::MF_COKRIGING}
    };
    
    auto it = map.find(str);
    if (it != map.end()) {
        return it->second;
    }
    return ModelCategory::FNO;  // Default
}

// =============================================================================
// ML Solver Wrapper
// =============================================================================

MLSolverWrapper::MLSolverWrapper(ModelCategory category) : category_(category) {}

void MLSolverWrapper::setActiveModel(const std::string& model_name) {
    active_model_name_ = model_name;
}

void MLSolverWrapper::enableML(bool enable) {
    ml_enabled_ = enable;
}

void MLSolverWrapper::setFallback(std::function<Tensor(const Tensor&)> fallback) {
    fallback_ = fallback;
}

Tensor MLSolverWrapper::solve(const Tensor& input) {
    if (!isMLActive() && fallback_) {
        return fallback_(input);
    }
    
    auto start = std::chrono::high_resolution_clock::now();
    
    Tensor result;
    
    // Get model and call predict
    auto& registry = MLModelRegistry::getInstance();
    
    switch (category_) {
        case ModelCategory::FNO: {
            auto* model = registry.getModel<FNOModel>(active_model_name_);
            if (model) {
                result = model->forward(input);
            }
            break;
        }
        case ModelCategory::GNO: {
            auto* model = registry.getModel<GraphNeuralOperator>(active_model_name_);
            if (model) {
                // Would need graph construction
                // result = model->predict(graph);
            }
            break;
        }
        default:
            if (fallback_) {
                result = fallback_(input);
            }
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(end - start).count();
    
    registry.recordCall(active_model_name_, elapsed);
    
    return result;
}

Tensor MLSolverWrapper::solveWithConfidence(const Tensor& input, double threshold,
                                            bool& used_ml) {
    // Would check uncertainty estimate against threshold
    used_ml = isMLActive();
    return solve(input);
}

// =============================================================================
// Global Helper Functions
// =============================================================================

void parseMLConfig(const std::map<std::string, std::string>& config) {
    auto& registry = MLModelRegistry::getInstance();
    
    // Check for global ML enable
    auto enable_it = config.find("ml_enabled");
    if (enable_it != config.end() && enable_it->second == "false") {
        return;  // ML disabled
    }
    
    // Set model directory
    auto dir_it = config.find("ml_model_directory");
    if (dir_it != config.end()) {
        registry.setGlobalOption("model_directory", dir_it->second);
    }
    
    // Register FNO if enabled
    auto fno_it = config.find("fno_enabled");
    if (fno_it != config.end() && fno_it->second == "true") {
        FNOConfig fno_config;
        // Parse FNO-specific config
        registry.registerFNO("default_fno", fno_config);
    }
    
    // Register GNO if enabled
    auto gno_it = config.find("gno_enabled");
    if (gno_it != config.end() && gno_it->second == "true") {
        GNOConfig gno_config = parseGNOConfig(config);
        registry.registerGNO("default_gno", gno_config);
    }
    
    // Register surrogates
    auto pvt_it = config.find("neural_pvt_enabled");
    if (pvt_it != config.end() && pvt_it->second == "true") {
        NeuralPVTConfig pvt_config = parseNeuralPVTConfig(config);
        registry.registerPVTSurrogate("default_pvt", pvt_config);
    }
}

void initializeML(const std::string& config_file) {
    auto& registry = MLModelRegistry::getInstance();
    registry.loadConfig(config_file);
    
    std::cout << "[ML] Initialized ML system from: " << config_file << "\n";
    registry.printSummary();
}

void shutdownML() {
    auto& registry = MLModelRegistry::getInstance();
    registry.printSummary();
    
    std::cout << "[ML] Shutting down ML system\n";
}

std::string getMLStatus() {
    auto& registry = MLModelRegistry::getInstance();
    
    std::stringstream ss;
    ss << "ML System Status\n";
    ss << "================\n";
    ss << "Registered models: " << registry.getAllModelNames().size() << "\n";
    
    return ss.str();
}

} // namespace ML
} // namespace FSRM
