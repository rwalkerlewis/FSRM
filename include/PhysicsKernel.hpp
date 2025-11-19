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

} // namespace ResSim

#endif // PHYSICS_KERNEL_HPP
