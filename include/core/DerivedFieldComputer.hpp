/**
 * @file DerivedFieldComputer.hpp
 * @brief Compute cell-centered stress, strain, and CFS from the FEM solution
 */

#ifndef FSRM_DERIVED_FIELD_COMPUTER_HPP
#define FSRM_DERIVED_FIELD_COMPUTER_HPP

#include <petscdmplex.h>
#include <petscfe.h>
#include <array>
#include <vector>

namespace FSRM
{

/**
 * @brief Computes cell-centered derived fields (stress, strain, CFS)
 *        from a DMPlex displacement solution using FE basis tabulation.
 *
 * Voigt ordering: [xx, yy, zz, xy, xz, yz]
 */
class DerivedFieldComputer
{
public:
  DerivedFieldComputer() = default;

  /**
   * @brief Initialize with material parameters
   *
   * @param lambda First Lame parameter (Pa)
   * @param mu Shear modulus (Pa)
   * @param biot_alpha Biot coefficient (dimensionless)
   */
  void setMaterialParameters(double lambda, double mu, double biot_alpha);

  /**
   * @brief Set receiver fault orientation for CFS computation
   *
   * @param strike_deg Strike angle in degrees (clockwise from north)
   * @param dip_deg Dip angle in degrees (0 = horizontal, 90 = vertical)
   * @param friction Static friction coefficient
   */
  void setReceiverOrientation(double strike_deg, double dip_deg, double friction);

  /**
   * @brief Compute cell-centered stress and strain for all cells
   *
   * Creates two Vecs on a cell-centered DM: stress (6 components)
   * and strain (6 components) in Voigt notation.
   *
   * @param dm The DMPlex mesh
   * @param solution Global solution vector
   * @return PetscErrorCode
   */
  PetscErrorCode compute(DM dm, Vec solution);

  /**
   * @brief Write stress to an HDF5 file
   *
   * @param comm MPI communicator
   * @param filename Output file path
   * @param step Timestep index
   * @param time Current simulation time
   * @return PetscErrorCode
   */
  PetscErrorCode writeStressHDF5(MPI_Comm comm, const char* filename,
                                 PetscInt step, PetscReal time);

  /**
   * @brief Write CFS to an HDF5 file
   *
   * @param comm MPI communicator
   * @param filename Output file path
   * @param step Timestep index
   * @param time Current simulation time
   * @return PetscErrorCode
   */
  PetscErrorCode writeCfsHDF5(MPI_Comm comm, const char* filename,
                              PetscInt step, PetscReal time);

  /**
   * @brief Get the number of cells
   */
  PetscInt numCells() const { return num_cells_; }

  /**
   * @brief Get cell-centered stress array (6 * numCells values, Voigt order)
   */
  const std::vector<double>& stress() const { return stress_; }

  /**
   * @brief Get cell-centered strain array (6 * numCells values, Voigt order)
   */
  const std::vector<double>& strain() const { return strain_; }

  /**
   * @brief Get cell-centered CFS array (numCells values)
   */
  const std::vector<double>& cfs() const { return cfs_; }

private:
  void computeStressFromStrain(const double eps[6], double P, double sigma[6]) const;
  void computeReceiverNormal();
  double computeCFS(const double sigma[6]) const;

  double lambda_ = 0.0;
  double mu_ = 0.0;
  double biot_alpha_ = 0.0;

  double strike_deg_ = 0.0;
  double dip_deg_ = 90.0;
  double friction_ = 0.4;
  std::array<double, 3> receiver_normal_ = {1.0, 0.0, 0.0};

  PetscInt num_cells_ = 0;
  std::vector<double> stress_;
  std::vector<double> strain_;
  std::vector<double> cfs_;
};

} // namespace FSRM

#endif // FSRM_DERIVED_FIELD_COMPUTER_HPP
