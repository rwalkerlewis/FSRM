#ifndef PHYSICS_KERNEL_HPP
#define PHYSICS_KERNEL_HPP

#include "ReservoirSim.hpp"
#include <petscfe.h>

namespace ResSim {

// Base class for all physics kernels
class PhysicsKernel {
public:
    PhysicsKernel(PhysicsType type) : physics_type(type) {}
    virtual ~PhysicsKernel() = default;
    
    // Setup kernel with PETSc
    virtual PetscErrorCode setup(DM dm, PetscFE fe) = 0;
    
    // Residual evaluation F(u) = 0
    virtual void residual(const PetscScalar u[], const PetscScalar u_t[],
                         const PetscScalar u_x[], const PetscScalar a[],
                         const PetscReal x[], PetscScalar f[]) = 0;
    
    // Jacobian evaluation dF/du
    virtual void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                         const PetscScalar u_x[], const PetscScalar a[],
                         const PetscReal x[], PetscScalar J[]) = 0;
    
    // Get number of fields and components
    virtual int getNumFields() const = 0;
    virtual int getNumComponents(int field) const = 0;
    
    PhysicsType getType() const { return physics_type; }
    
protected:
    PhysicsType physics_type;
};

// Single phase flow kernel (Darcy flow)
class SinglePhaseFlowKernel : public PhysicsKernel {
public:
    SinglePhaseFlowKernel();
    
    PetscErrorCode setup(DM dm, PetscFE fe) override;
    
    void residual(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar f[]) override;
    
    void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar J[]) override;
    
    int getNumFields() const override { return 1; }
    int getNumComponents(int field) const override { return 1; }
    
    void setProperties(double phi, double k, double ct, double mu, double rho);
    
private:
    double porosity;
    double permeability;
    double compressibility;
    double viscosity;
    double density;
};

// Black oil model kernel (oil, gas, water with dissolution)
class BlackOilKernel : public PhysicsKernel {
public:
    BlackOilKernel();
    
    PetscErrorCode setup(DM dm, PetscFE fe) override;
    
    void residual(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar f[]) override;
    
    void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar J[]) override;
    
    int getNumFields() const override { return 3; } // Po, Sw, Sg
    int getNumComponents(int field) const override { return 1; }
    
    void setFluidProperties(const FluidProperties& props);
    void setRockProperties(double phi, double kx, double ky, double kz);
    
    // PVT functions
    double oilDensity(double P, double Rs) const;
    double gasDensity(double P) const;
    double waterDensity(double P) const;
    double oilViscosity(double P, double Rs) const;
    double gasViscosity(double P) const;
    double waterViscosity(double P) const;
    double solutionGOR(double P) const;
    
    // Relative permeability
    double krw(double Sw) const;
    double kro(double So, double Sg) const;
    double krg(double Sg) const;
    
    // Capillary pressure
    double Pcow(double Sw) const;
    double Pcog(double Sg) const;
    
private:
    FluidProperties fluid_props;
    double porosity;
    double perm_x, perm_y, perm_z;
};

// Compositional flow kernel (multi-component with EOS)
class CompositionalKernel : public PhysicsKernel {
public:
    CompositionalKernel(int num_components);
    
    PetscErrorCode setup(DM dm, PetscFE fe) override;
    
    void residual(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar f[]) override;
    
    void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar J[]) override;
    
    int getNumFields() const override { return nc + 1; } // P + nc-1 compositions
    int getNumComponents(int field) const override { return 1; }
    
    void setComponentProperties(const std::vector<double>& mw,
                               const std::vector<double>& Tc,
                               const std::vector<double>& Pc,
                               const std::vector<double>& omega);
    
    // Equation of state (Peng-Robinson)
    void flashCalculation(double P, double T, const std::vector<double>& z,
                         std::vector<double>& x, std::vector<double>& y,
                         double& S_L, double& S_V) const;
    
    double fugacityCoefficient(double P, double T, double z_factor,
                              const std::vector<double>& comp, int i) const;
    
private:
    int nc; // number of components
    std::vector<double> component_mw;
    std::vector<double> component_Tc;
    std::vector<double> component_Pc;
    std::vector<double> component_omega;
    double porosity;
    double perm_x, perm_y, perm_z;
};

// Geomechanics kernel (linear elasticity or viscoelasticity)
class GeomechanicsKernel : public PhysicsKernel {
public:
    GeomechanicsKernel(SolidModelType model_type);
    
    PetscErrorCode setup(DM dm, PetscFE fe) override;
    
    void residual(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar f[]) override;
    
    void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar J[]) override;
    
    int getNumFields() const override { return 1; }
    int getNumComponents(int field) const override { return 3; } // ux, uy, uz
    
    void setMaterialProperties(double E, double nu, double rho);
    void setViscoelasticProperties(double tau, double eta);
    void setPoroelasticCoupling(double alpha, double M);
    
private:
    SolidModelType model_type;
    double youngs_modulus;
    double poisson_ratio;
    double density;
    double relaxation_time;  // for viscoelasticity
    double viscosity;         // for viscoelasticity
    double biot_coefficient;
    double biot_modulus;
};

// Thermal kernel (heat conduction and convection)
class ThermalKernel : public PhysicsKernel {
public:
    ThermalKernel();
    
    PetscErrorCode setup(DM dm, PetscFE fe) override;
    
    void residual(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar f[]) override;
    
    void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar J[]) override;
    
    int getNumFields() const override { return 1; } // Temperature
    int getNumComponents(int field) const override { return 1; }
    
    void setThermalProperties(double k, double rho, double cp);
    void setFluidThermalProperties(double k_f, double rho_f, double cp_f);
    
private:
    double thermal_conductivity_solid;
    double density_solid;
    double heat_capacity_solid;
    double thermal_conductivity_fluid;
    double density_fluid;
    double heat_capacity_fluid;
    double porosity;
};

// Particle transport kernel (proppant, tracers, etc.)
class ParticleTransportKernel : public PhysicsKernel {
public:
    ParticleTransportKernel();
    
    PetscErrorCode setup(DM dm, PetscFE fe) override;
    
    void residual(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar f[]) override;
    
    void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar J[]) override;
    
    int getNumFields() const override { return 1; } // Concentration
    int getNumComponents(int field) const override { return 1; }
    
    void setParticleProperties(double diameter, double density, double diffusivity);
    void enableGravitationalSettling(bool enable);
    void enableBridging(bool enable);
    
private:
    double particle_diameter;
    double particle_density;
    double diffusivity;
    bool gravity_settling;
    bool enable_bridging;
    double porosity;
    double permeability;
};

// Fracture propagation kernel (cohesive zone model)
class FracturePropagationKernel : public PhysicsKernel {
public:
    FracturePropagationKernel();
    
    PetscErrorCode setup(DM dm, PetscFE fe) override;
    
    void residual(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar f[]) override;
    
    void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar J[]) override;
    
    int getNumFields() const override { return 1; } // Fracture width
    int getNumComponents(int field) const override { return 1; }
    
    void setFractureProperties(double Kc, double Gc, double sigma_c);
    
private:
    double fracture_toughness;
    double fracture_energy;
    double critical_stress;
};

// Tidal forces kernel
class TidalForcesKernel : public PhysicsKernel {
public:
    TidalForcesKernel();
    
    PetscErrorCode setup(DM dm, PetscFE fe) override;
    
    void residual(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar f[]) override;
    
    void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar J[]) override;
    
    int getNumFields() const override { return 0; } // Body force addition
    int getNumComponents(int field) const override { return 0; }
    
    void setLocationAndTime(double lat, double lon, double time);
    void computeTidalStress(const PetscReal x[], PetscScalar stress[]);
    
private:
    double latitude;
    double longitude;
    double current_time;
};

// Elastodynamics kernel (wave propagation with inertia)
class ElastodynamicsKernel : public PhysicsKernel {
public:
    ElastodynamicsKernel();
    
    PetscErrorCode setup(DM dm, PetscFE fe) override;
    
    void residual(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar f[]) override;
    
    void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar J[]) override;
    
    int getNumFields() const override { return 1; }
    int getNumComponents(int field) const override { return 3; } // ux, uy, uz
    
    void setMaterialProperties(double E, double nu, double rho);
    void setWaveProperties(double vp, double vs, double Q);
    void setDamping(double alpha, double beta);
    void setStaticTriggeringMode(bool enable, double threshold, double duration);
    
    // Check if static stress exceeds threshold
    bool checkStaticTrigger(const PetscScalar stress[], double current_time);
    bool isInDynamicEvent(double current_time) const;
    
private:
    double youngs_modulus;
    double poisson_ratio;
    double density;
    double p_wave_velocity;
    double s_wave_velocity;
    double quality_factor;
    double damping_alpha;  // Rayleigh mass damping
    double damping_beta;   // Rayleigh stiffness damping
    
    // Static triggering parameters
    bool enable_static_triggering;
    double trigger_threshold;
    double event_duration;
    double trigger_time;
    bool event_active;
};

// Poroelastodynamics kernel (Biot's equations with dynamics)
class PoroelastodynamicsKernel : public PhysicsKernel {
public:
    PoroelastodynamicsKernel();
    
    PetscErrorCode setup(DM dm, PetscFE fe) override;
    
    void residual(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar f[]) override;
    
    void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar J[]) override;
    
    int getNumFields() const override { return 2; } // Displacement + Pressure
    int getNumComponents(int field) const override { return (field == 0 ? 3 : 1); }
    
    void setMaterialProperties(double E, double nu, double rho, double phi);
    void setFluidProperties(double rho_f, double mu, double K_f);
    void setBiotParameters(double alpha, double M);
    void setWaveProperties(double vp_fast, double vs, double vp_slow);
    void setDamping(double alpha, double beta);
    void setStaticTriggeringMode(bool enable, double threshold, double duration);
    
    // Permeability dynamics
    void enableDynamicPermeabilityChange(bool enable, double strain_coeff, 
                                         double stress_coeff, double recovery_time);
    double computePermeabilityChange(const PetscScalar u_x[], const PetscScalar stress[], 
                                     double k_initial, double dt);
    
    bool checkStaticTrigger(const PetscScalar stress[], double current_time);
    bool isInDynamicEvent(double current_time) const;
    
private:
    // Solid properties
    double youngs_modulus;
    double poisson_ratio;
    double density_solid;
    double porosity;
    
    // Fluid properties
    double density_fluid;
    double viscosity;
    double bulk_modulus_fluid;
    
    // Biot parameters
    double biot_coefficient;
    double biot_modulus;
    
    // Wave properties
    double p_wave_fast;   // Fast P-wave (Biot wave I)
    double s_wave;        // Shear wave
    double p_wave_slow;   // Slow P-wave (Biot wave II)
    double damping_alpha;
    double damping_beta;
    
    // Permeability
    double permeability;
    double permeability_initial;
    
    // Permeability dynamics
    bool enable_dynamic_k;
    double k_strain_coeff;
    double k_stress_coeff;
    double k_recovery_time;
    
    // Static triggering
    bool enable_static_triggering;
    double trigger_threshold;
    double event_duration;
    double trigger_time;
    bool event_active;
};

// Dynamic permeability change model (can be used by other kernels)
class DynamicPermeabilityModel {
public:
    DynamicPermeabilityModel();
    
    void setParameters(double k0, double strain_coeff, double stress_coeff,
                      double recovery_time, double k_min, double k_max);
    
    // Compute instantaneous permeability change from wave passage
    double computeInstantaneousChange(double strain_amplitude, double stress_amplitude);
    
    // Compute time-dependent permeability evolution
    double computeTimeEvolution(double k_current, double k0, double dt);
    
    // Mechanistic models for permeability change
    double strainModel(double volumetric_strain);
    double stressModel(double effective_stress);
    double shearDilationModel(double shear_strain);
    double crackOpeningModel(double normal_stress);
    
private:
    double k_initial;
    double strain_coefficient;
    double stress_coefficient;
    double recovery_tau;
    double k_minimum;
    double k_maximum;
};

} // namespace ResSim

#endif // PHYSICS_KERNEL_HPP
