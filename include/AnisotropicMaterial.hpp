#ifndef ANISOTROPIC_MATERIAL_HPP
#define ANISOTROPIC_MATERIAL_HPP

#include "MaterialModel.hpp"
#include <array>
#include <vector>
#include <map>
#include <string>

namespace FSRM {

/**
 * @brief Anisotropic elastic material model
 * 
 * Implements fully anisotropic elasticity with 21 independent elastic constants.
 * The stress-strain relation is:
 * 
 * σ_ij = c_ijkl ε_kl
 * 
 * Due to symmetry, the 81-component tensor reduces to 21 independent components.
 * 
 * Special cases:
 * - Isotropic: 2 constants (λ, μ or E, ν)
 * - Transversely isotropic: 5 constants
 * - Orthotropic: 9 constants
 * - Triclinic: 21 constants
 * 
 * Applications:
 * - Layered sedimentary rocks
 * - Metamorphic rocks with preferred orientation
 * - Fractured media with oriented fracture sets
 * - Stress-induced anisotropy
 */
class AnisotropicMaterial {
public:
    /**
     * @brief Elastic stiffness tensor in Voigt notation
     * 
     * Voigt notation maps tensor indices (ij) to vector indices:
     * 11→1, 22→2, 33→3, 23→4, 13→5, 12→6
     * 
     * Stiffness matrix (6x6 symmetric):
     * | c11 c12 c13 c14 c15 c16 |
     * |     c22 c23 c24 c25 c26 |
     * |         c33 c34 c35 c36 |
     * |             c44 c45 c46 |
     * |                 c55 c56 |
     * |                     c66 |
     */
    struct StiffnessTensor {
        // 21 independent components
        double c11, c12, c13, c14, c15, c16;
        double     c22, c23, c24, c25, c26;
        double         c33, c34, c35, c36;
        double             c44, c45, c46;
        double                 c55, c56;
        double                     c66;
        
        StiffnessTensor() { setIsotropic(30e9, 0.25); }
        
        // Initialize from Lamé parameters (isotropic)
        void setIsotropic(double lambda, double mu);
        void setIsotropic(double E, double nu, bool use_youngs);
        
        // Initialize transversely isotropic (5 parameters)
        void setTransverselyIsotropic(double E_parallel, double E_perpendicular,
                                      double nu_parallel, double nu_perpendicular,
                                      double G_perpendicular);
        
        // Initialize orthotropic (9 parameters)
        void setOrthotropic(double E1, double E2, double E3,
                           double nu12, double nu23, double nu31,
                           double G12, double G23, double G31);
        
        // Set from array
        void setFromArray(const double* c_array);  // Length 21
        void setFromMatrix(const double* c_matrix); // 6x6 matrix
        
        // Get as arrays
        void toArray(double* c_array) const;
        void toMatrix(double* c_matrix) const;  // 6x6 matrix
        
        // Get specific component
        double get(int i, int j) const;
        void set(int i, int j, double value);
        
        // Compute compliance tensor (inverse)
        StiffnessTensor inverse() const;
        
        // Symmetry check
        bool checkSymmetry(double tol = 1e-10) const;
        
        // Positive definiteness check (stability)
        bool isStable() const;
        
        // Compute wave speeds in direction n
        struct WaveSpeeds {
            double vp;   // P-wave (fastest)
            double vs1;  // Fast S-wave
            double vs2;  // Slow S-wave
            std::array<double, 3> p_direction;   // P-wave polarization
            std::array<double, 3> s1_direction;  // S1 polarization
            std::array<double, 3> s2_direction;  // S2 polarization
        };
        WaveSpeeds computeWaveSpeeds(const std::array<double, 3>& direction, double density) const;
        
        // Rotate stiffness tensor
        StiffnessTensor rotate(const std::array<std::array<double, 3>, 3>& R) const;
    };
    
    /**
     * @brief Symmetry class of anisotropic material
     */
    enum class SymmetryClass {
        ISOTROPIC,              // 2 constants
        CUBIC,                  // 3 constants
        HEXAGONAL,              // 5 constants (transversely isotropic)
        TETRAGONAL_1,           // 6 constants
        TETRAGONAL_2,           // 7 constants
        TRIGONAL_1,             // 6 constants
        TRIGONAL_2,             // 7 constants
        ORTHOTROPIC,            // 9 constants
        MONOCLINIC,             // 13 constants
        TRICLINIC               // 21 constants (fully anisotropic)
    };
    
    AnisotropicMaterial();
    
    // Set stiffness tensor
    void setStiffness(const StiffnessTensor& C);
    const StiffnessTensor& getStiffness() const { return C; }
    
    // Set density
    void setDensity(double rho) { density = rho; }
    double getDensity() const { return density; }
    
    // Compute stress from strain
    void computeStress(const double* strain, double* stress) const;
    void computeStress(const std::array<double, 6>& strain,
                      std::array<double, 6>& stress) const;
    
    // Compute strain from stress (using compliance)
    void computeStrain(const double* stress, double* strain) const;
    
    // Compute tangent stiffness matrix (for FEM)
    void getTangentStiffness(double* K) const;  // 6x6 matrix
    
    // Detect symmetry class
    SymmetryClass detectSymmetry(double tol = 1e-8) const;
    
    // Rotation
    void rotate(const std::array<std::array<double, 3>, 3>& rotation_matrix);
    void rotateToLocal(const std::array<double, 3>& x_local,
                      const std::array<double, 3>& y_local,
                      const std::array<double, 3>& z_local);
    
    // Coordinate system for anisotropy (e.g., bedding planes)
    void setCoordinateSystem(const std::array<double, 3>& strike_dir,
                            const std::array<double, 3>& dip_dir,
                            const std::array<double, 3>& normal_dir);
    
    // Parse from config
    void configure(const std::map<std::string, std::string>& config);
    
private:
    StiffnessTensor C;            // Stiffness tensor
    StiffnessTensor S;            // Compliance tensor (inverse)
    double density;
    bool compliance_computed;
    
    void computeCompliance();
};

/**
 * @brief Transversely isotropic material (special case)
 * 
 * Transverse isotropy has one axis of symmetry (e.g., layering direction).
 * Common in sedimentary rocks and shales.
 * 
 * 5 independent constants: E_parallel, E_perpendicular, ν_parallel, ν_perpendicular, G
 */
class TransverselyIsotropicMaterial {
public:
    struct Parameters {
        double E_parallel;         // Young's modulus parallel to symmetry axis
        double E_perpendicular;    // Young's modulus perpendicular to axis
        double nu_parallel;        // Poisson's ratio (parallel)
        double nu_perpendicular;   // Poisson's ratio (perpendicular)
        double G_perpendicular;    // Shear modulus perpendicular to axis
        double density;            // Mass density
        
        // Symmetry axis direction (default: vertical)
        std::array<double, 3> axis = {0.0, 0.0, 1.0};
        
        Parameters() :
            E_parallel(50e9),
            E_perpendicular(30e9),
            nu_parallel(0.25),
            nu_perpendicular(0.30),
            G_perpendicular(10e9),
            density(2500.0) {}
    };
    
    TransverselyIsotropicMaterial();
    
    void setParameters(const Parameters& params);
    const Parameters& getParameters() const { return params; }
    
    // Convert to full anisotropic material
    AnisotropicMaterial toAnisotropic() const;
    
    // Wave speeds
    double vpParallel() const;      // P-wave parallel to axis
    double vpPerpendicular() const; // P-wave perpendicular to axis
    double vsParallel() const;      // S-wave parallel
    double vsPerpendicular() const; // S-wave perpendicular
    
private:
    Parameters params;
};

/**
 * @brief Orthorhombic material (9 constants)
 * 
 * Orthorhombic symmetry has three mutually perpendicular symmetry planes.
 * Common in fractured rocks with two orthogonal fracture sets.
 */
class OrthorhombicMaterial {
public:
    struct Parameters {
        double E1, E2, E3;              // Young's moduli in three directions
        double nu12, nu23, nu31;        // Poisson's ratios
        double G12, G23, G31;           // Shear moduli
        double density;
        
        // Principal directions (default: aligned with x,y,z)
        std::array<double, 3> dir1 = {1.0, 0.0, 0.0};
        std::array<double, 3> dir2 = {0.0, 1.0, 0.0};
        std::array<double, 3> dir3 = {0.0, 0.0, 1.0};
        
        Parameters() :
            E1(40e9), E2(35e9), E3(30e9),
            nu12(0.25), nu23(0.28), nu31(0.27),
            G12(12e9), G23(11e9), G31(11.5e9),
            density(2500.0) {}
    };
    
    OrthorhombicMaterial();
    
    void setParameters(const Parameters& params);
    const Parameters& getParameters() const { return params; }
    
    // Convert to full anisotropic material
    AnisotropicMaterial toAnisotropic() const;
    
private:
    Parameters params;
};

/**
 * @brief Effective medium theory for fractured rocks
 * 
 * Compute effective anisotropic properties of rock with oriented fractures
 * using various effective medium theories.
 */
class FracturedRockEMT {
public:
    enum class Theory {
        HUDSON,             // Hudson's crack model
        SCHOENBERG,         // Schoenberg's linear slip model
        ESHELBY,            // Eshelby inclusion theory
        SELF_CONSISTENT,    // Self-consistent method
        DIFFERENTIAL        // Differential effective medium
    };
    
    struct FractureSet {
        std::array<double, 3> normal;  // Normal to fracture plane
        double density;                 // Fracture density (cracks per unit volume)
        double aspect_ratio;            // Width/length ratio
        double compliance;              // Fracture compliance (displacement/stress)
        
        FractureSet() :
            normal({0.0, 0.0, 1.0}),
            density(0.05),
            aspect_ratio(0.001),
            compliance(1e-10) {}
    };
    
    FracturedRockEMT(Theory theory = Theory::HUDSON);
    
    // Set background (intact rock) properties
    void setBackgroundMaterial(double lambda, double mu, double density);
    
    // Add fracture set
    void addFractureSet(const FractureSet& fractures);
    void clearFractureSets() { fracture_sets.clear(); }
    
    // Compute effective anisotropic properties
    AnisotropicMaterial computeEffectiveProperties() const;
    
private:
    Theory theory;
    double lambda_bg, mu_bg, density_bg;
    std::vector<FractureSet> fracture_sets;
    
    // Theory-specific implementations
    AnisotropicMaterial hudsonMethod() const;
    AnisotropicMaterial schoenbergMethod() const;
    AnisotropicMaterial eshelbyMethod() const;
};

/**
 * @brief Stress-induced anisotropy
 * 
 * Account for stress-induced changes in elastic properties due to:
 * - Crack closure/opening
 * - Nonlinear elasticity
 * - Preferred orientation development
 */
class StressInducedAnisotropy {
public:
    StressInducedAnisotropy();
    
    // Set initial (unstressed) properties
    void setInitialMaterial(const AnisotropicMaterial& material);
    
    // Update properties based on current stress state
    void updateFromStress(const std::array<double, 6>& stress);
    
    // Get current effective properties
    const AnisotropicMaterial& getCurrentMaterial() const { return current_material; }
    
    // Enable/disable stress dependence
    void enable(bool enable) { enabled = enable; }
    
private:
    bool enabled;
    AnisotropicMaterial initial_material;
    AnisotropicMaterial current_material;
    
    // Stress sensitivity parameters
    double pressure_sensitivity;    // ∂c/∂p
    double shear_sensitivity;       // ∂c/∂τ
};

/**
 * @brief Configuration helper for anisotropic materials
 */
struct AnisotropicConfig {
    bool enabled = false;
    std::string symmetry_class = "isotropic";  // isotropic, transverse, orthotropic, triclinic
    
    // Stiffness components (for direct specification)
    std::map<std::string, double> stiffness;
    
    // Transverse isotropy parameters
    double E_parallel = 50e9;
    double E_perpendicular = 30e9;
    double nu_parallel = 0.25;
    double nu_perpendicular = 0.30;
    double G_perpendicular = 10e9;
    
    // Symmetry axis orientation (strike, dip for transverse isotropy)
    double axis_strike = 0.0;       // degrees from North
    double axis_dip = 90.0;         // degrees from horizontal (90 = vertical)
    
    // Parse from config file
    void parseConfig(const std::map<std::string, std::string>& config);
    
    // Create material
    AnisotropicMaterial createMaterial() const;
};

} // namespace FSRM

#endif // ANISOTROPIC_MATERIAL_HPP
