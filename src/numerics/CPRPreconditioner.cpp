/**
 * @file CPRPreconditioner.cpp
 * @brief Implementation of CPR (Constrained Pressure Residual) PC wrapper
 */

#include "numerics/CPRPreconditioner.hpp"

#include <petscsection.h>
#include <petscvec.h>

namespace FSRM {

PetscErrorCode CPRPreconditioner::setup(PC pc, DM dm, PetscInt pressure_field) {
    PetscFunctionBeginUser;
    PetscCall(PCSetType(pc, PCFIELDSPLIT));
    PetscCall(PCFieldSplitSetType(pc, PC_COMPOSITE_MULTIPLICATIVE));
    PetscCall(setupFieldSplit(pc, dm, pressure_field));
    PetscFunctionReturn(0);
}

PetscErrorCode CPRPreconditioner::setupFieldSplit(PC pc, DM dm, PetscInt pressure_field) {
    PetscFunctionBeginUser;

    PetscSection local_section = NULL;
    PetscSection global_section = NULL;
    IS is_pressure = NULL;
    IS is_rest = NULL;
    PetscInt *indices = NULL;
    PetscInt n_idx = 0;
    PetscInt n_splits = 0;
    KSP *sub_ksp = NULL;

    PetscCall(DMGetLocalSection(dm, &local_section));
    PetscCall(DMGetGlobalSection(dm, &global_section));
    PetscCheck(local_section, PetscObjectComm((PetscObject)dm), PETSC_ERR_ARG_WRONGSTATE,
               "DM must have a local section (e.g. from DMPlex)");
    PetscCheck(global_section, PetscObjectComm((PetscObject)dm), PETSC_ERR_ARG_WRONGSTATE,
               "DM must have a global section");

    PetscInt p_start = 0;
    PetscInt p_end = 0;
    PetscCall(PetscSectionGetChart(global_section, &p_start, &p_end));

    for (PetscInt p = p_start; p < p_end; ++p) {
        PetscInt fdof = 0;
        PetscCall(PetscSectionGetFieldDof(local_section, p, pressure_field, &fdof));
        if (fdof <= 0) {
            continue;
        }
        PetscInt foff = 0;
        PetscCall(PetscSectionGetFieldOffset(global_section, p, pressure_field, &foff));
        if (foff < 0) {
            continue; /* not in the global vector (e.g. ghost / unowned) */
        }
        n_idx += fdof;
    }

    PetscCheck(n_idx > 0, PetscObjectComm((PetscObject)dm), PETSC_ERR_ARG_OUTOFRANGE,
               "Pressure field %" PetscInt_FMT " has no degrees of freedom on this DM", pressure_field);

    PetscCall(PetscMalloc1(n_idx, &indices));
    {
        PetscInt k = 0;
        for (PetscInt p = p_start; p < p_end; ++p) {
            PetscInt fdof = 0;
            PetscCall(PetscSectionGetFieldDof(local_section, p, pressure_field, &fdof));
            if (fdof <= 0) {
                continue;
            }
            PetscInt foff = 0;
            PetscCall(PetscSectionGetFieldOffset(global_section, p, pressure_field, &foff));
            if (foff < 0) {
                continue;
            }
            for (PetscInt j = 0; j < fdof; ++j) {
                indices[k++] = foff + j;
            }
        }
        PetscCheck(k == n_idx, PETSC_COMM_SELF, PETSC_ERR_PLIB, "Index count mismatch building pressure IS");
    }

    PetscCall(ISCreateGeneral(PetscObjectComm((PetscObject)dm), n_idx, indices, PETSC_OWN_POINTER, &is_pressure));
    indices = NULL;
    PetscCall(ISSort(is_pressure));

    Vec glob = NULL;
    PetscCall(DMGetGlobalVector(dm, &glob));
    PetscInt N = 0;
    PetscCall(VecGetSize(glob, &N));
    PetscCall(DMRestoreGlobalVector(dm, &glob));

    PetscCall(ISComplement(is_pressure, 0, N, &is_rest));

    PetscCall(PCFieldSplitSetIS(pc, "pressure", is_pressure));
    PetscCall(PCFieldSplitSetIS(pc, "rest", is_rest));

    PetscCall(PCSetUp(pc));

    PetscCall(PCFieldSplitGetSubKSP(pc, &n_splits, &sub_ksp));
    PetscCheck(n_splits == 2, PetscObjectComm((PetscObject)pc), PETSC_ERR_PLIB,
               "Expected 2 field-split blocks for CPR, got %" PetscInt_FMT, n_splits);

    {
        PC sub_pc = NULL;
        PetscCall(KSPSetType(sub_ksp[0], KSPPREONLY));
        PetscCall(KSPGetPC(sub_ksp[0], &sub_pc));
        PetscCall(PCSetType(sub_pc, PCGAMG));
    }
    {
        PC sub_pc = NULL;
        PetscCall(KSPSetType(sub_ksp[1], KSPPREONLY));
        PetscCall(KSPGetPC(sub_ksp[1], &sub_pc));
        PetscCall(PCSetType(sub_pc, PCILU));
    }

    PetscCall(PCSetUp(pc));

    PetscCall(PetscFree(sub_ksp));
    PetscCall(ISDestroy(&is_pressure));
    PetscCall(ISDestroy(&is_rest));

    PetscFunctionReturn(0);
}

PetscErrorCode CPRPreconditioner::configureFromOptions(KSP ksp, DM dm) {
    PetscFunctionBeginUser;

    PC pc = NULL;
    PetscCall(KSPGetPC(ksp, &pc));
    PetscCall(setup(pc, dm, 0));
    PetscCall(KSPSetFromOptions(ksp));

    PetscFunctionReturn(0);
}

} // namespace FSRM
