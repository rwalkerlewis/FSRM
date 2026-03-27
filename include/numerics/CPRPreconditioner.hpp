#ifndef CPR_PRECONDITIONER_HPP
#define CPR_PRECONDITIONER_HPP

/**
 * @file CPRPreconditioner.hpp
 * @brief CPR (Constrained Pressure Residual) preconditioner wrapper for PETSc
 */

#include <petscksp.h>
#include <petscdmplex.h>
#include <petscpc.h>

namespace FSRM {

/**
 * @brief CPR (Constrained Pressure Residual) preconditioner
 *
 * Industry-standard two-stage preconditioner for reservoir simulation:
 * - Stage 1: AMG (algebraic multigrid) on pressure block
 * - Stage 2: ILU on full coupled system
 *
 * Uses PETSc's PCFIELDSPLIT with PC_COMPOSITE_MULTIPLICATIVE.
 */
class CPRPreconditioner {
public:
    static PetscErrorCode setup(PC pc, DM dm, PetscInt pressure_field);
    static PetscErrorCode configureFromOptions(KSP ksp, DM dm);
    static PetscErrorCode setupFieldSplit(PC pc, DM dm, PetscInt pressure_field);
};

} // namespace FSRM

#endif // CPR_PRECONDITIONER_HPP
