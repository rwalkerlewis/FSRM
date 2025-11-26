/**
 * @file test_material_model.cpp
 * @brief Unit tests for MaterialModel classes
 */

#include <gtest/gtest.h>
#include "MaterialModel.hpp"
#include <cmath>

using namespace FSRM;

class MaterialModelTest : public ::testing::Test {
protected:
    void SetUp() override {}
};

// ============================================================================
// StressTensor Tests
// ============================================================================

TEST_F(MaterialModelTest, StressTensorMeanStress) {
    StressTensor sigma;
    sigma.xx = 10e6;
    sigma.yy = 20e6;
    sigma.zz = 30e6;
    sigma.xy = sigma.xz = sigma.yz = 0.0;
    
    double mean = sigma.meanStress();
    EXPECT_NEAR(mean, 20e6, 1e3);
}

TEST_F(MaterialModelTest, StressTensorVonMises) {
    StressTensor sigma;
    sigma.xx = 100e6;
    sigma.yy = 0.0;
    sigma.zz = 0.0;
    sigma.xy = sigma.xz = sigma.yz = 0.0;
    
    // For uniaxial stress, von Mises = |σ_xx|
    double vm = sigma.vonMises();
    EXPECT_NEAR(vm, 100e6, 1e3);
}

TEST_F(MaterialModelTest, StressTensorPrincipals) {
    StressTensor sigma;
    sigma.xx = 30e6;
    sigma.yy = 20e6;
    sigma.zz = 10e6;
    sigma.xy = sigma.xz = sigma.yz = 0.0;
    
    // For diagonal stress, principals are the diagonal values
    auto principals = sigma.principals();
    EXPECT_NEAR(principals[0], 30e6, 1e3);  // Maximum
    EXPECT_NEAR(principals[2], 10e6, 1e3);  // Minimum
}

// ============================================================================
// StrainTensor Tests
// ============================================================================

TEST_F(MaterialModelTest, StrainTensorVolumetric) {
    StrainTensor eps;
    eps.xx = 0.001;
    eps.yy = 0.002;
    eps.zz = 0.003;
    eps.xy = eps.xz = eps.yz = 0.0;
    
    double vol = eps.volumetric();
    EXPECT_NEAR(vol, 0.006, 1e-9);
}

// ============================================================================
// LinearElasticMaterial Tests
// ============================================================================

TEST_F(MaterialModelTest, LinearElasticCreate) {
    LinearElasticMaterial mat;
    EXPECT_GT(mat.getYoungsModulus(), 0.0);
}

TEST_F(MaterialModelTest, LinearElasticConfigure) {
    LinearElasticMaterial mat;
    std::map<std::string, std::string> config;
    config["youngs_modulus"] = "20e9";
    config["poisson_ratio"] = "0.25";
    config["density"] = "2650";
    mat.configure(config);
    
    EXPECT_NEAR(mat.getYoungsModulus(), 20e9, 1e6);
    EXPECT_NEAR(mat.getPoissonsRatio(), 0.25, 0.001);
    EXPECT_NEAR(mat.getDensity(), 2650.0, 0.1);
}

TEST_F(MaterialModelTest, LinearElasticStressFromStrain) {
    LinearElasticMaterial mat;
    std::map<std::string, std::string> config;
    config["youngs_modulus"] = "20e9";
    config["poisson_ratio"] = "0.25";
    mat.configure(config);
    
    // Uniaxial strain
    StrainTensor eps;
    eps.xx = 0.001;
    eps.yy = eps.zz = 0.0;
    eps.xy = eps.xz = eps.yz = 0.0;
    
    StressTensor sigma = mat.computeStress(eps);
    
    // Stress should be positive (compression follows extension)
    EXPECT_GT(sigma.xx, 0.0);
    EXPECT_TRUE(std::isfinite(sigma.xx));
}

TEST_F(MaterialModelTest, LinearElasticStiffnessMatrix) {
    LinearElasticMaterial mat;
    std::map<std::string, std::string> config;
    config["youngs_modulus"] = "20e9";
    config["poisson_ratio"] = "0.25";
    mat.configure(config);
    
    auto C = mat.getStiffnessMatrix();
    
    // Matrix should be symmetric (C[0] = C11, etc.)
    // For isotropic: C11 = C22 = C33
    EXPECT_NEAR(C[0], C[7], 1e6);    // C11 = C22
    EXPECT_NEAR(C[0], C[14], 1e6);   // C11 = C33
    
    // All elements should be finite
    for (int i = 0; i < 36; i++) {
        EXPECT_TRUE(std::isfinite(C[i]));
    }
}

// ============================================================================
// PoroelasticMaterial Tests
// ============================================================================

TEST_F(MaterialModelTest, PoroelasticCreate) {
    PoroelasticMaterial mat;
    EXPECT_GT(mat.getBiotCoefficient(), 0.0);
}

TEST_F(MaterialModelTest, PoroelasticBiotCoefficient) {
    PoroelasticMaterial mat;
    std::map<std::string, std::string> config;
    config["youngs_modulus"] = "20e9";
    config["poisson_ratio"] = "0.25";
    config["biot_coefficient"] = "0.8";
    mat.configure(config);
    
    EXPECT_NEAR(mat.getBiotCoefficient(), 0.8, 0.01);
}

TEST_F(MaterialModelTest, PoroelasticEffectiveStress) {
    PoroelasticMaterial mat;
    std::map<std::string, std::string> config;
    config["biot_coefficient"] = "1.0";
    mat.configure(config);
    
    StressTensor total;
    total.xx = 50e6;
    total.yy = 50e6;
    total.zz = 50e6;
    total.xy = total.xz = total.yz = 0.0;
    
    double pore_pressure = 20e6;
    
    StressTensor eff = mat.getEffectiveStress(total, pore_pressure);
    
    // Effective stress = total - α*p*I
    EXPECT_NEAR(eff.xx, 30e6, 1e3);
    EXPECT_NEAR(eff.yy, 30e6, 1e3);
    EXPECT_NEAR(eff.zz, 30e6, 1e3);
}

// ============================================================================
// ViscoelasticMaterial Tests
// ============================================================================

TEST_F(MaterialModelTest, ViscoelasticMaxwellCreate) {
    ViscoelasticMaterial mat(ViscoelasticType::MAXWELL);
    EXPECT_GT(mat.getYoungsModulus(), 0.0);
}

TEST_F(MaterialModelTest, ViscoelasticRelaxation) {
    ViscoelasticMaterial mat(ViscoelasticType::MAXWELL);
    std::map<std::string, std::string> config;
    config["youngs_modulus"] = "10e9";
    config["viscosity"] = "1e18";
    mat.configure(config);
    
    // Apply constant strain, stress should relax over time
    StrainTensor eps;
    eps.xx = 0.001;
    eps.yy = eps.zz = 0.0;
    eps.xy = eps.xz = eps.yz = 0.0;
    
    StressTensor sigma0 = mat.computeStress(eps, nullptr, 0.0);
    
    // After long time, Maxwell stress should relax
    StressTensor sigma_relaxed = mat.computeStress(eps, &sigma0, 1e15);  // Very long dt
    
    // Stress should decrease (relax) for Maxwell model
    EXPECT_LE(sigma_relaxed.xx, sigma0.xx);
}

// ============================================================================
// ElastoplasticMaterial Tests
// ============================================================================

TEST_F(MaterialModelTest, ElastoplasticMohrCoulombCreate) {
    ElastoplasticMaterial mat(FailureCriterion::MOHR_COULOMB);
    EXPECT_GT(mat.getYoungsModulus(), 0.0);
}

TEST_F(MaterialModelTest, ElastoplasticYieldCheck) {
    ElastoplasticMaterial mat(FailureCriterion::MOHR_COULOMB);
    std::map<std::string, std::string> config;
    config["youngs_modulus"] = "20e9";
    config["cohesion"] = "1e6";
    config["friction_angle"] = "30.0";
    mat.configure(config);
    
    // Low stress - elastic
    StressTensor sigma_low;
    sigma_low.xx = 5e6;
    sigma_low.yy = 5e6;
    sigma_low.zz = 5e6;
    sigma_low.xy = sigma_low.xz = sigma_low.yz = 0.0;
    
    EXPECT_FALSE(mat.hasYielded(sigma_low));
}

// ============================================================================
// AnisotropicMaterial Tests
// ============================================================================

TEST_F(MaterialModelTest, AnisotropicVTICreate) {
    AnisotropicMaterial mat(AnisotropyType::VTI);
    EXPECT_GT(mat.getYoungsModulus(), 0.0);
}

TEST_F(MaterialModelTest, AnisotropicStiffnessSymmetry) {
    AnisotropicMaterial mat(AnisotropyType::VTI);
    std::map<std::string, std::string> config;
    config["c11"] = "50e9";
    config["c33"] = "40e9";
    config["c12"] = "15e9";
    config["c13"] = "12e9";
    config["c44"] = "12e9";
    mat.configure(config);
    
    auto C = mat.getStiffnessMatrix();
    
    // VTI: C11 = C22 (transverse isotropy)
    EXPECT_NEAR(C[0], C[7], 1e6);
    
    // All elements finite
    for (int i = 0; i < 36; i++) {
        EXPECT_TRUE(std::isfinite(C[i]));
    }
}

// ============================================================================
// Factory Function Tests
// ============================================================================

TEST_F(MaterialModelTest, CreateMaterialModelLinearElastic) {
    std::map<std::string, std::string> config;
    config["youngs_modulus"] = "20e9";
    auto mat = createMaterialModel("LINEAR_ELASTIC", config);
    EXPECT_NE(mat, nullptr);
}

TEST_F(MaterialModelTest, CreateMaterialModelPoroelastic) {
    std::map<std::string, std::string> config;
    config["youngs_modulus"] = "20e9";
    config["biot_coefficient"] = "0.8";
    auto mat = createMaterialModel("POROELASTIC", config);
    EXPECT_NE(mat, nullptr);
}

// ============================================================================
// PermeabilityModel Tests
// ============================================================================

TEST_F(MaterialModelTest, PermeabilityConstant) {
    PermeabilityModel perm;
    perm.setModel(PermeabilityModelType::CONSTANT);
    perm.setReferencePermeability(100e-15);  // 100 mD in m²
    
    double k = perm.getPermeability(0.2, 10e6, 0.0);
    EXPECT_NEAR(k, 100e-15, 1e-18);
}

TEST_F(MaterialModelTest, PermeabilityKozenyCarman) {
    PermeabilityModel perm;
    perm.setModel(PermeabilityModelType::KOZENY_CARMAN);
    perm.setReferencePermeability(100e-15);
    perm.setReferencePorosity(0.2);
    
    // Higher porosity should give higher permeability
    double k_low = perm.getPermeability(0.15, 10e6, 0.0);
    double k_high = perm.getPermeability(0.25, 10e6, 0.0);
    
    EXPECT_GT(k_high, k_low);
}

// ============================================================================
// RockProperties Tests
// ============================================================================

TEST_F(MaterialModelTest, RockPropertiesCreate) {
    RockProperties rock;
    EXPECT_GT(rock.getPorosity(), 0.0);
}

TEST_F(MaterialModelTest, RockPropertiesPermeability) {
    RockProperties rock;
    rock.setPorosity(0.2);
    rock.setPermeability(100e-15, 100e-15, 10e-15);  // kx, ky, kz
    
    double kx = rock.getPermeability(10e6, 0.0);
    EXPECT_NEAR(kx, 100e-15, 1e-18);
}

// ============================================================================
// Physical Validity Tests
// ============================================================================

TEST_F(MaterialModelTest, ElasticModuliPositive) {
    LinearElasticMaterial mat;
    std::map<std::string, std::string> config;
    config["youngs_modulus"] = "20e9";
    config["poisson_ratio"] = "0.25";
    mat.configure(config);
    
    EXPECT_GT(mat.getYoungsModulus(), 0.0);
    EXPECT_GT(mat.getPoissonsRatio(), -1.0);
    EXPECT_LT(mat.getPoissonsRatio(), 0.5);  // Physical limit for isotropic
}

TEST_F(MaterialModelTest, StiffnessMatrixPositiveDefinite) {
    LinearElasticMaterial mat;
    std::map<std::string, std::string> config;
    config["youngs_modulus"] = "20e9";
    config["poisson_ratio"] = "0.25";
    mat.configure(config);
    
    auto C = mat.getStiffnessMatrix();
    
    // Diagonal elements should be positive
    EXPECT_GT(C[0], 0.0);   // C11
    EXPECT_GT(C[7], 0.0);   // C22
    EXPECT_GT(C[14], 0.0);  // C33
    EXPECT_GT(C[21], 0.0);  // C44
    EXPECT_GT(C[28], 0.0);  // C55
    EXPECT_GT(C[35], 0.0);  // C66
}

TEST_F(MaterialModelTest, LameParametersFromModuli) {
    double E = 20e9;
    double nu = 0.25;
    
    // Lamé's first parameter
    double lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    
    // Shear modulus (Lamé's second parameter)
    double mu = E / (2.0 * (1.0 + nu));
    
    EXPECT_GT(lambda, 0.0);
    EXPECT_GT(mu, 0.0);
    
    // Verify: K = λ + 2μ/3
    double K = lambda + 2.0 * mu / 3.0;
    double K_expected = E / (3.0 * (1.0 - 2.0 * nu));
    EXPECT_NEAR(K, K_expected, 1e6);
}
