#ifndef COULOMB_STRESS_TRANSFER_HPP
#define COULOMB_STRESS_TRANSFER_HPP

/**
 * @file CoulombStressTransfer.hpp
 * @brief Samples stress tensor from FEM solution at fault vertices
 *
 * Bridges the continuum FEM displacement/pressure fields and the
 * discrete fault mechanics by:
 * 1. Computing strain from the displacement gradient at fault locations
 * 2. Computing stress via Hooke's law + Biot effective stress
 * 3. Projecting stress onto fault orientation (normal, shear)
 * 4. Evaluating delta_CFS relative to initial (pre-injection) state
 *
 * This replaces the placeholder computeVonMisesStress() and
 * computeCoulombFailure() in ImplicitExplicitTransitionManager.
 *
 * References:
 *   Segall & Lu (2015): Injection-induced seismicity, Rev. Geophys.
 */

#include <petscdmplex.h>
#include <petscvec.h>
#include <array>
#include <vector>

namespace FSRM {

class FaultCohesiveDyn;

/**
 * @brief Samples stress from FEM solution and computes CFS at fault vertices
 */
class CoulombStressTransfer {
public:
    CoulombStressTransfer(MPI_Comm comm);
    ~CoulombStressTransfer() = default;

    /**
     * @brief Initialize with DM and fault topology
     *
     * Stores references and sets up material properties for stress computation.
     *
     * @param dm The DMPlex mesh (may have cohesive cells)
     * @param fault The FaultCohesiveDyn with vertex topology
     * @param lambda First Lamé parameter (Pa)
     * @param mu Shear modulus (Pa)
     * @param biot_alpha Biot coefficient
     * @return PETSc error code
     */
    PetscErrorCode initialize(DM dm, FaultCohesiveDyn* fault,
                               double lambda, double mu, double biot_alpha);

    /**
     * @brief Sample stress tensor at fault vertices from FEM displacement field
     *
     * For each fault vertex, finds the adjacent volume cell, extracts the
     * local displacement solution, computes the displacement gradient
     * (strain), and then computes stress via Hooke's law.
     *
     * Also samples pressure for effective stress computation.
     *
     * @param solution The global solution vector (displacement + pressure fields)
     * @return PETSc error code
     */
    PetscErrorCode sampleStressAtFaults(Vec solution);

    /**
     * @brief Project stress tensor onto each fault vertex's orientation
     *
     * Computes:
     *   sigma_n = n^T * sigma * n  (normal stress on fault)
     *   tau = sqrt(|sigma*n|^2 - sigma_n^2)  (shear stress on fault)
     *
     * @return PETSc error code
     */
    PetscErrorCode resolveStressOnFault();

    /**
     * @brief Compute delta_CFS for each vertex relative to initial state
     *
     * delta_CFS = delta_tau - mu_s * delta_sigma_n_eff
     * Positive delta_CFS means moving toward failure.
     *
     * @param static_friction Static friction coefficient
     * @return PETSc error code
     */
    PetscErrorCode computeDeltaCFS(double static_friction);

    /**
     * @brief Update FaultCohesiveDyn's DynamicFaultState with current stress values
     *
     * Writes sigma_n_eff, tau, and is_locked/has_ruptured flags.
     *
     * @param fault The FaultCohesiveDyn to update
     * @return PETSc error code
     */
    PetscErrorCode updateFaultState(FaultCohesiveDyn* fault);

    /**
     * @brief Store current stress state as the initial (pre-injection) reference
     */
    PetscErrorCode storeInitialStress();

    /**
     * @brief Get the maximum delta_CFS across all fault vertices
     */
    double getMaxDeltaCFS() const;

    /**
     * @brief Get the vertex index with maximum delta_CFS
     */
    int getMaxCFSVertex() const;

    /**
     * @brief Per-vertex stress state
     */
    struct VertexStress {
        std::array<double, 6> sigma = {0, 0, 0, 0, 0, 0};  // Voigt: xx,yy,zz,xy,xz,yz
        double pressure = 0.0;
        double sigma_n_eff = 0.0;       // Effective normal stress on fault
        double tau = 0.0;                // Shear stress magnitude on fault
        double delta_cfs = 0.0;          // Change in Coulomb failure stress
    };

    const std::vector<VertexStress>& getCurrentStress() const { return current_; }

private:
    MPI_Comm comm_;
    DM dm_ = nullptr;
    FaultCohesiveDyn* fault_ = nullptr;

    // Material properties
    double lambda_ = 0.0;       // First Lamé parameter
    double mu_ = 0.0;           // Shear modulus
    double biot_alpha_ = 1.0;   // Biot coefficient

    // Per-vertex stress state
    std::vector<VertexStress> current_;
    std::vector<VertexStress> initial_;
    bool has_initial_ = false;

    /**
     * @brief Compute stress tensor from strain using Hooke's law
     *
     * sigma_ij = lambda * eps_kk * delta_ij + 2 * mu * eps_ij
     * sigma'_ij = sigma_ij - alpha * P * delta_ij  (effective stress)
     */
    void computeStressFromStrain(const double eps[6], double P,
                                  double sigma[6]) const;
};

} // namespace FSRM

#endif // COULOMB_STRESS_TRANSFER_HPP
