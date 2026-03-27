#ifndef MULTIPHASE_FLOW_KERNEL_HPP
#define MULTIPHASE_FLOW_KERNEL_HPP

/**
 * @file MultiphaseFlowKernel.hpp
 * @brief PetscDS-compatible multiphase flow kernel for CO2-brine systems
 *
 * Provides pointwise residual (f0, f1) and Jacobian (g0, g1, g2, g3) callbacks
 * for coupled pressure-saturation equations in the PETSc finite element framework.
 *
 * Field layout (registered via PetscDSSetResidual/PetscDSSetJacobian):
 *   Field 0: Pressure (P)       [Pa]
 *   Field 1: Water saturation (Sw) [dimensionless]
 *
 * The system models two-phase (water/brine + CO2) flow in porous media:
 *   Pressure equation: conservation of total fluid mass
 *   Saturation equation: conservation of water phase
 *
 * All static member functions follow the exact PetscPointFunc / PetscPointJac
 * signatures for use with PetscDSSetResidual() and PetscDSSetJacobian().
 *
 * Constants layout for PetscDSSetConstants():
 *   [0]  porosity
 *   [1]  permeability_x       [m²]
 *   [2]  permeability_y       [m²]
 *   [3]  permeability_z       [m²]
 *   [4]  water_viscosity      [Pa·s]
 *   [5]  co2_viscosity        [Pa·s]
 *   [6]  water_density        [kg/m³]
 *   [7]  co2_density          [kg/m³]
 *   [8]  water_compressibility [1/Pa]
 *   [9]  co2_compressibility   [1/Pa]
 *   [10] Swr (water residual saturation)
 *   [11] Sgr (CO2 residual/critical saturation)
 *   [12] corey_nw (water Corey exponent)
 *   [13] corey_ng (CO2 Corey exponent)
 *   [14] krw_max
 *   [15] krg_max
 *   [16] gravity              [m/s²]
 *   [17] biot_coefficient
 */

#include "physics/PhysicsKernel.hpp"

namespace FSRM {

class MultiphaseFlowKernel : public PhysicsKernel {
public:
    MultiphaseFlowKernel();

    PetscErrorCode setup(DM dm, PetscFE fe) override;

    void residual(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar f[]) override;

    void jacobian(const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscScalar a[],
                 const PetscReal x[], PetscScalar J[]) override;

    int getNumFields() const override { return 2; }
    int getNumComponents(int field) const override { (void)field; return 1; }

    // =========================================================================
    // Property setters
    // =========================================================================
    void setRockProperties(double phi, double kx, double ky, double kz);
    void setWaterProperties(double rho, double mu, double cw);
    void setCO2Properties(double rho, double mu, double cg);
    void setRelPermParams(double Swr, double Sgr, double nw, double ng,
                          double krw_max, double krg_max);

    // =========================================================================
    // Static PetscDS-compatible pointwise RESIDUAL functions
    // =========================================================================

    /**
     * @brief Pressure equation f0: accumulation term
     *
     * f0_p = φ * [Sw * cw + (1-Sw) * cg] * dP/dt
     *       + φ * (ρw - ρg) * dSw/dt     (density-difference coupling)
     */
    static void f0_pressure(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                            const PetscInt uOff[], const PetscInt uOff_x[],
                            const PetscScalar u[], const PetscScalar u_t[],
                            const PetscScalar u_x[],
                            const PetscInt aOff[], const PetscInt aOff_x[],
                            const PetscScalar a[], const PetscScalar a_t[],
                            const PetscScalar a_x[],
                            PetscReal t, const PetscReal x[],
                            PetscInt numConstants, const PetscScalar constants[],
                            PetscScalar f0[]);

    /**
     * @brief Pressure equation f1: total mobility flux
     *
     * f1_p = Σ_α [k * kr_α/μ_α * (∇P - ρ_α * g * ê_z)]
     */
    static void f1_pressure(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                            const PetscInt uOff[], const PetscInt uOff_x[],
                            const PetscScalar u[], const PetscScalar u_t[],
                            const PetscScalar u_x[],
                            const PetscInt aOff[], const PetscInt aOff_x[],
                            const PetscScalar a[], const PetscScalar a_t[],
                            const PetscScalar a_x[],
                            PetscReal t, const PetscReal x[],
                            PetscInt numConstants, const PetscScalar constants[],
                            PetscScalar f1[]);

    /**
     * @brief Saturation equation f0: accumulation
     *
     * f0_s = φ * dSw/dt
     */
    static void f0_saturation(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                              const PetscInt uOff[], const PetscInt uOff_x[],
                              const PetscScalar u[], const PetscScalar u_t[],
                              const PetscScalar u_x[],
                              const PetscInt aOff[], const PetscInt aOff_x[],
                              const PetscScalar a[], const PetscScalar a_t[],
                              const PetscScalar a_x[],
                              PetscReal t, const PetscReal x[],
                              PetscInt numConstants, const PetscScalar constants[],
                              PetscScalar f0[]);

    /**
     * @brief Saturation equation f1: water phase Darcy flux
     *
     * f1_s = k * kr_w/μ_w * (∇P - ρ_w * g * ê_z)
     */
    static void f1_saturation(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                              const PetscInt uOff[], const PetscInt uOff_x[],
                              const PetscScalar u[], const PetscScalar u_t[],
                              const PetscScalar u_x[],
                              const PetscInt aOff[], const PetscInt aOff_x[],
                              const PetscScalar a[], const PetscScalar a_t[],
                              const PetscScalar a_x[],
                              PetscReal t, const PetscReal x[],
                              PetscInt numConstants, const PetscScalar constants[],
                              PetscScalar f1[]);

    // =========================================================================
    // Static PetscDS-compatible pointwise JACOBIAN functions
    // =========================================================================

    /** @brief g0 for (pressure, pressure): ∂f0_p/∂P */
    static void g0_pp(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                     const PetscInt uOff[], const PetscInt uOff_x[],
                     const PetscScalar u[], const PetscScalar u_t[],
                     const PetscScalar u_x[],
                     const PetscInt aOff[], const PetscInt aOff_x[],
                     const PetscScalar a[], const PetscScalar a_t[],
                     const PetscScalar a_x[],
                     PetscReal t, PetscReal u_tShift, const PetscReal x[],
                     PetscInt numConstants, const PetscScalar constants[],
                     PetscScalar g0[]);

    /** @brief g3 for (pressure, pressure): ∂f1_p/∂(∇P) — total mobility tensor */
    static void g3_pp(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                     const PetscInt uOff[], const PetscInt uOff_x[],
                     const PetscScalar u[], const PetscScalar u_t[],
                     const PetscScalar u_x[],
                     const PetscInt aOff[], const PetscInt aOff_x[],
                     const PetscScalar a[], const PetscScalar a_t[],
                     const PetscScalar a_x[],
                     PetscReal t, PetscReal u_tShift, const PetscReal x[],
                     PetscInt numConstants, const PetscScalar constants[],
                     PetscScalar g3[]);

    /** @brief g0 for (pressure, saturation): ∂f0_p/∂Sw — cross-coupling */
    static void g0_ps(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                     const PetscInt uOff[], const PetscInt uOff_x[],
                     const PetscScalar u[], const PetscScalar u_t[],
                     const PetscScalar u_x[],
                     const PetscInt aOff[], const PetscInt aOff_x[],
                     const PetscScalar a[], const PetscScalar a_t[],
                     const PetscScalar a_x[],
                     PetscReal t, PetscReal u_tShift, const PetscReal x[],
                     PetscInt numConstants, const PetscScalar constants[],
                     PetscScalar g0[]);

    /** @brief g0 for (saturation, saturation): ∂f0_s/∂Sw */
    static void g0_ss(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                     const PetscInt uOff[], const PetscInt uOff_x[],
                     const PetscScalar u[], const PetscScalar u_t[],
                     const PetscScalar u_x[],
                     const PetscInt aOff[], const PetscInt aOff_x[],
                     const PetscScalar a[], const PetscScalar a_t[],
                     const PetscScalar a_x[],
                     PetscReal t, PetscReal u_tShift, const PetscReal x[],
                     PetscInt numConstants, const PetscScalar constants[],
                     PetscScalar g0[]);

    /** @brief g3 for (saturation, pressure): ∂f1_s/∂(∇P) — water mobility tensor */
    static void g3_sp(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                     const PetscInt uOff[], const PetscInt uOff_x[],
                     const PetscScalar u[], const PetscScalar u_t[],
                     const PetscScalar u_x[],
                     const PetscInt aOff[], const PetscInt aOff_x[],
                     const PetscScalar a[], const PetscScalar a_t[],
                     const PetscScalar a_x[],
                     PetscReal t, PetscReal u_tShift, const PetscReal x[],
                     PetscInt numConstants, const PetscScalar constants[],
                     PetscScalar g3[]);

private:
    double porosity_;
    double perm_x_, perm_y_, perm_z_;
    double water_density_, water_viscosity_, water_compressibility_;
    double co2_density_, co2_viscosity_, co2_compressibility_;
    double Swr_, Sgr_, nw_, ng_, krw_max_, krg_max_;

    // =========================================================================
    // Helper: Corey relative permeability from constants array
    // =========================================================================
    static inline double computeKrw(PetscScalar Sw, const PetscScalar constants[]) {
        PetscScalar Swr = constants[10];
        PetscScalar Sgr = constants[11];
        PetscScalar nw  = constants[12];
        PetscScalar krw_max = constants[14];
        PetscScalar Se = (Sw - Swr) / (1.0 - Swr - Sgr);
        Se = std::max(0.0, std::min(1.0, PetscRealPart(Se)));
        return PetscRealPart(krw_max) * std::pow(Se, PetscRealPart(nw));
    }

    static inline double computeKrg(PetscScalar Sw, const PetscScalar constants[]) {
        PetscScalar Swr = constants[10];
        PetscScalar Sgr = constants[11];
        PetscScalar ng  = constants[13];
        PetscScalar krg_max = constants[15];
        PetscScalar Sg = 1.0 - Sw;
        PetscScalar Se = (Sg - Sgr) / (1.0 - Swr - Sgr);
        Se = std::max(0.0, std::min(1.0, PetscRealPart(Se)));
        return PetscRealPart(krg_max) * std::pow(Se, PetscRealPart(ng));
    }

    static inline double computeDkrw_dSw(PetscScalar Sw, const PetscScalar constants[]) {
        PetscScalar Swr = constants[10];
        PetscScalar Sgr = constants[11];
        PetscScalar nw  = constants[12];
        PetscScalar krw_max = constants[14];
        double denom = PetscRealPart(1.0 - Swr - Sgr);
        if (denom < 1e-15) return 0.0;
        double Se = (PetscRealPart(Sw) - PetscRealPart(Swr)) / denom;
        Se = std::max(0.0, std::min(1.0, Se));
        if (Se < 1e-15 || Se > 1.0 - 1e-15) return 0.0;
        return PetscRealPart(krw_max) * PetscRealPart(nw)
             * std::pow(Se, PetscRealPart(nw) - 1.0) / denom;
    }

    static inline double computeDkrg_dSw(PetscScalar Sw, const PetscScalar constants[]) {
        // dkrg/dSw = dkrg/dSg * dSg/dSw = -dkrg/dSg
        PetscScalar Swr = constants[10];
        PetscScalar Sgr = constants[11];
        PetscScalar ng  = constants[13];
        PetscScalar krg_max = constants[15];
        double denom = PetscRealPart(1.0 - Swr - Sgr);
        if (denom < 1e-15) return 0.0;
        double Sg = PetscRealPart(1.0 - Sw);
        double Se = (Sg - PetscRealPart(Sgr)) / denom;
        Se = std::max(0.0, std::min(1.0, Se));
        if (Se < 1e-15 || Se > 1.0 - 1e-15) return 0.0;
        double dkrg_dSg = PetscRealPart(krg_max) * PetscRealPart(ng)
                        * std::pow(Se, PetscRealPart(ng) - 1.0) / denom;
        return -dkrg_dSg;  // Chain rule: dSg/dSw = -1
    }
};

} // namespace FSRM

#endif // MULTIPHASE_FLOW_KERNEL_HPP
