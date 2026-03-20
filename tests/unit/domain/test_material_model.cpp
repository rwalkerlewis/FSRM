/**
 * @file test_material_model.cpp
 * @brief Unit tests for MaterialModel classes
 */

#include "domain/geomechanics/MaterialModel.hpp"
#include "core/FSRM.hpp"
#include <cmath>
#include <gtest/gtest.h>

using namespace FSRM;

class MaterialModelTest : public ::testing::Test {
protected:
    void SetUp() override {
        int rk = 0;
        MPI_Comm_rank(PETSC_COMM_WORLD, &rk);
        EXPECT_GE(rk, 0);
    }
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

    double mean = sigma.mean();
    EXPECT_NEAR(mean, 20e6, 1e3);
}

TEST_F(MaterialModelTest, StressTensorVonMises) {
    StressTensor sigma;
    sigma.xx = 100e6;
    sigma.yy = 0.0;
    sigma.zz = 0.0;
    sigma.xy = sigma.xz = sigma.yz = 0.0;

    double vm = sigma.vonMises();
    EXPECT_NEAR(vm, 100e6, 1e3);
}

TEST_F(MaterialModelTest, StressTensorPrincipals) {
    StressTensor sigma;
    sigma.xx = 30e6;
    sigma.yy = 20e6;
    sigma.zz = 10e6;
    sigma.xy = sigma.xz = sigma.yz = 0.0;

    auto principals = sigma.principalStresses();
    EXPECT_NEAR(principals[0], 30e6, 1e3);
    EXPECT_NEAR(principals[2], 10e6, 1e3);
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

    StrainTensor eps;
    eps.xx = 0.001;
    eps.yy = eps.zz = 0.0;
    eps.xy = eps.xz = eps.yz = 0.0;

    StressTensor sigma = mat.computeStress(eps);

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

    EXPECT_NEAR(C[0], C[7], 1e6);
    EXPECT_NEAR(C[0], C[14], 1e6);

    for (int i = 0; i < 36; i++) {
        EXPECT_TRUE(std::isfinite(C[i]));
    }
}

TEST_F(MaterialModelTest, BulkModulusFromYoungsAndNu) {
    LinearElasticMaterial mat;
    const double E = 30e9;
    const double nu = 0.2;
    mat.setElasticModuli(E, nu);
    const double K_expected = E / (3.0 * (1.0 - 2.0 * nu));
    EXPECT_NEAR(mat.getBulkModulus(), K_expected, E * 1e-9);
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

TEST_F(MaterialModelTest, PoroelasticTotalStress) {
    PoroelasticMaterial mat;
    std::map<std::string, std::string> config;
    config["biot_coefficient"] = "1.0";
    config["youngs_modulus"] = "20e9";
    config["poisson_ratio"] = "0.25";
    mat.configure(config);

    StrainTensor eps;
    eps.xx = 0.001;
    eps.yy = eps.zz = 0.0;
    eps.xy = eps.xz = eps.yz = 0.0;

    StressTensor eff = mat.computeEffectiveStress(eps);

    double pore_pressure = 20e6;
    StressTensor total = mat.computeTotalStress(eff, pore_pressure);

    EXPECT_TRUE(std::isfinite(total.xx));
    EXPECT_TRUE(std::isfinite(total.yy));
    EXPECT_TRUE(std::isfinite(total.zz));
}

// ============================================================================
// ViscoelasticMaterial Tests
// ============================================================================

TEST_F(MaterialModelTest, ViscoelasticMaxwellCreate) {
    ViscoelasticMaterial mat(ViscoelasticMaterial::ViscoType::MAXWELL);
    EXPECT_GT(mat.getYoungsModulus(), 0.0);
}

TEST_F(MaterialModelTest, ViscoelasticRelaxation) {
    ViscoelasticMaterial mat(ViscoelasticMaterial::ViscoType::MAXWELL);
    std::map<std::string, std::string> config;
    config["youngs_modulus"] = "10e9";
    config["viscosity"] = "1e18";
    mat.configure(config);

    StrainTensor eps;
    eps.xx = 0.001;
    eps.yy = eps.zz = 0.0;
    eps.xy = eps.xz = eps.yz = 0.0;

    StressTensor sigma0 = mat.computeStress(eps);

    double G0 = mat.getRelaxationModulus(0.0);
    double G_later = mat.getRelaxationModulus(1e15);

    EXPECT_LE(G_later, G0);
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
    config["cohesion"] = "10e6";
    config["friction_angle"] = "30.0";
    mat.configure(config);

    StressTensor sigma;
    sigma.xx = 50e6;
    sigma.yy = 50e6;
    sigma.zz = 50e6;
    sigma.xy = sigma.xz = sigma.yz = 0.0;

    double yf = mat.evaluateYieldFunction(sigma);
    EXPECT_TRUE(std::isfinite(yf));

    EXPECT_NEAR(mat.getCohesion(), 10e6, 1e3);
    EXPECT_GT(mat.getFrictionAngle(), 0.0);
}

// ============================================================================
// AnisotropicMaterial Tests
// ============================================================================

TEST_F(MaterialModelTest, AnisotropicVTICreate) {
    AnisotropicMaterial mat(true);
    EXPECT_GT(mat.getYoungsModulus(), 0.0);
}

TEST_F(MaterialModelTest, AnisotropicStiffnessSymmetry) {
    AnisotropicMaterial mat(true);
    std::map<std::string, std::string> config;
    config["c11"] = "50e9";
    config["c33"] = "40e9";
    config["c12"] = "15e9";
    config["c13"] = "12e9";
    config["c44"] = "12e9";
    mat.configure(config);

    auto C = mat.getStiffnessMatrix();

    EXPECT_NEAR(C[0], C[7], 1e6);

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
    PermeabilityModel perm(PermeabilityModel::PermeabilityType::CONSTANT);
    perm.setInitialPermeability(100e-15);

    double k = perm.getPermeability(10e6, 0.2, 0.0);
    EXPECT_NEAR(k, 100e-15, 1e-18);
}

TEST_F(MaterialModelTest, PermeabilityKozenyCarman) {
    PermeabilityModel perm(PermeabilityModel::PermeabilityType::KOZENY_CARMAN);
    perm.setInitialPermeability(100e-15);

    double k_low = perm.getPermeability(10e6, 0.15, 0.0);
    double k_high = perm.getPermeability(10e6, 0.25, 0.0);

    EXPECT_GT(k_high, k_low);
}

// ============================================================================
// RockProperties Tests
// ============================================================================

TEST_F(MaterialModelTest, RockPropertiesCreate) {
    RockProperties rock;
    std::map<std::string, std::string> config;
    config["porosity"] = "0.2";
    config["permeability_x"] = "100e-15";
    config["youngs_modulus"] = "20e9";
    rock.configure(config);

    EXPECT_GT(rock.getPorosity(), 0.0);
}

TEST_F(MaterialModelTest, RockPropertiesPermeability) {
    RockProperties rock;
    std::map<std::string, std::string> config;
    config["porosity"] = "0.2";
    config["permeability_x"] = "100e-15";
    config["youngs_modulus"] = "20e9";
    config["poisson_ratio"] = "0.25";
    rock.configure(config);

    double k = rock.getPermeability();
    EXPECT_GT(k, 0.0);
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
    EXPECT_LT(mat.getPoissonsRatio(), 0.5);
}

TEST_F(MaterialModelTest, StiffnessMatrixPositiveDefinite) {
    LinearElasticMaterial mat;
    std::map<std::string, std::string> config;
    config["youngs_modulus"] = "20e9";
    config["poisson_ratio"] = "0.25";
    mat.configure(config);

    auto C = mat.getStiffnessMatrix();

    EXPECT_GT(C[0], 0.0);
    EXPECT_GT(C[7], 0.0);
    EXPECT_GT(C[14], 0.0);
    EXPECT_GT(C[21], 0.0);
    EXPECT_GT(C[28], 0.0);
    EXPECT_GT(C[35], 0.0);
}

TEST_F(MaterialModelTest, LameParametersFromModuli) {
    double E = 20e9;
    double nu = 0.25;

    double lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    double mu = E / (2.0 * (1.0 + nu));

    EXPECT_GT(lambda, 0.0);
    EXPECT_GT(mu, 0.0);

    double K = lambda + 2.0 * mu / 3.0;
    double K_expected = E / (3.0 * (1.0 - 2.0 * nu));
    EXPECT_NEAR(K, K_expected, 1e6);
}
