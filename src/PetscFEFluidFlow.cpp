#include "PetscFEFluidFlow.hpp"

#include <algorithm>
#include <cmath>

namespace FSRM {
namespace PetscFEFluidFlow {

// -----------------------------------------------------------------------------
// Helpers / shared utilities
// -----------------------------------------------------------------------------

namespace {
// Constants layout used by black-oil PETScFE callbacks:
//   [0]  phi        (porosity)
//   [1]  kx         (m^2)
//   [2]  ky         (m^2)
//   [3]  kz         (m^2)
//   [4]  cw         (1/Pa)
//   [5]  co         (1/Pa)
//   [6]  cg         (1/Pa)
//   [7]  mu_w       (Pa*s)
//   [8]  mu_o       (Pa*s)
//   [9]  mu_g       (Pa*s)
//   [10] Swr        (-)
//   [11] Sor        (-)
//   [12] Sgr        (-)
//   [13] nw         (-)
//   [14] no         (-)
//   [15] ng         (-)
//   [16] krw0       (-)
//   [17] kro0       (-)
//   [18] krg0       (-)

static inline PetscScalar clamp01(PetscScalar v) {
    return v < 0.0 ? 0.0 : (v > 1.0 ? 1.0 : v);
}

static inline PetscScalar safeDenom(PetscScalar v) {
    const PetscScalar eps = 1e-12;
    if (std::abs(v) < eps) return (v >= 0.0) ? eps : -eps;
    return v;
}

static inline PetscScalar getConst(PetscInt numConstants, const PetscScalar c[], PetscInt idx, PetscScalar fallback) {
    return (numConstants > idx) ? c[idx] : fallback;
}

static inline void blackOilRelPermCorey(const PetscScalar Sw, const PetscScalar Sg,
                                        const PetscScalar Swr, const PetscScalar Sor, const PetscScalar Sgr,
                                        const PetscScalar nw, const PetscScalar no, const PetscScalar ng,
                                        const PetscScalar krw0, const PetscScalar kro0, const PetscScalar krg0,
                                        PetscScalar& krw, PetscScalar& kro, PetscScalar& krg) {
    PetscScalar So = 1.0 - Sw - Sg;
    if (So < 0.0) So = 0.0;

    const PetscScalar denom = safeDenom(1.0 - Swr - Sor - Sgr);
    const PetscScalar Swe = clamp01((Sw - Swr) / denom);
    const PetscScalar Sge = clamp01((Sg - Sgr) / denom);
    const PetscScalar Soe = clamp01((So - Sor) / denom);

    krw = krw0 * std::pow(Swe, nw);
    kro = kro0 * std::pow(Soe, no);
    krg = krg0 * std::pow(Sge, ng);
}
} // namespace

// -----------------------------------------------------------------------------
// Single-phase pressure equation
// -----------------------------------------------------------------------------

void f0_SinglePhase(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[],
                    const PetscScalar u[], const PetscScalar u_t[],
                    const PetscScalar u_x[], const PetscInt aOff[],
                    const PetscInt aOff_x[], const PetscScalar a[],
                    const PetscScalar a_x[], const PetscScalar a_t[],
                    PetscReal t,
                    const PetscReal x[], PetscInt numConstants,
                    const PetscScalar constants[], PetscScalar f0[]) {
    // Suppress unused parameter warnings - part of PETSc pointwise function interface
    (void)dim; (void)Nf; (void)NfAux;
    (void)uOff_x; (void)u; (void)u_x;
    (void)aOff; (void)aOff_x;
    (void)a; (void)a_x; (void)a_t;
    (void)t; (void)x;
    (void)numConstants; (void)constants;

    // Accumulation term: phi * ct * dP/dt
    const PetscReal phi = 0.2;   // TODO: plumb from constants/auxiliary fields
    const PetscReal ct  = 1e-9;  // TODO: plumb from constants/auxiliary fields
    f0[0] = phi * ct * u_t[uOff[0]];
}

void f1_SinglePhase(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                    const PetscInt uOff[], const PetscInt uOff_x[],
                    const PetscScalar u[], const PetscScalar u_t[],
                    const PetscScalar u_x[], const PetscInt aOff[],
                    const PetscInt aOff_x[], const PetscScalar a[],
                    const PetscScalar a_x[], const PetscScalar a_t[],
                    PetscReal t,
                    const PetscReal x[], PetscInt numConstants,
                    const PetscScalar constants[], PetscScalar f1[]) {
    // Suppress unused parameter warnings - part of PETSc pointwise function interface
    (void)Nf; (void)NfAux;
    (void)uOff; (void)u; (void)u_t;
    (void)aOff; (void)aOff_x;
    (void)a; (void)a_x; (void)a_t;
    (void)t; (void)x;
    (void)numConstants; (void)constants;

    // Flux term: -div( k/mu * grad(P) ) in strong form -> provide coefficient for weak form.
    const PetscReal k  = 100.0e-15; // m^2
    const PetscReal mu = 1e-3;      // Pa*s
    const PetscReal mobility = k / mu;
    for (PetscInt d = 0; d < dim; ++d) {
        f1[d] = mobility * u_x[uOff_x[0] + d];
    }
}

// -----------------------------------------------------------------------------
// Black-oil pressure equation
// -----------------------------------------------------------------------------

void f0_BlackOilPressure(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                         const PetscInt uOff[], const PetscInt uOff_x[],
                         const PetscScalar u[], const PetscScalar u_t[],
                         const PetscScalar u_x[], const PetscInt aOff[],
                         const PetscInt aOff_x[], const PetscScalar a[],
                         const PetscScalar a_x[], const PetscScalar a_t[],
                         PetscReal t,
                         const PetscReal x[], PetscInt numConstants,
                         const PetscScalar constants[], PetscScalar f0[]) {
    (void)dim; (void)Nf; (void)NfAux;
    (void)uOff_x; (void)u_x;
    (void)aOff; (void)aOff_x; (void)a; (void)a_x; (void)a_t;
    (void)t; (void)x;

    const PetscScalar Sw = u[uOff[1]];
    const PetscScalar Sg = u[uOff[2]];
    PetscScalar So = 1.0 - Sw - Sg;
    if (So < 0.0) So = 0.0;

    const PetscScalar phi = getConst(numConstants, constants, 0, 0.2);
    const PetscScalar cw  = getConst(numConstants, constants, 4, 4.5e-10);
    const PetscScalar co  = getConst(numConstants, constants, 5, 1.5e-9);
    const PetscScalar cg  = getConst(numConstants, constants, 6, 1.0e-8);

    const PetscScalar ct_eff = Sw * cw + So * co + Sg * cg;
    f0[0] = phi * ct_eff * u_t[uOff[0]];
}

void f1_BlackOilPressure(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                         const PetscInt uOff[], const PetscInt uOff_x[],
                         const PetscScalar u[], const PetscScalar u_t[],
                         const PetscScalar u_x[], const PetscInt aOff[],
                         const PetscInt aOff_x[], const PetscScalar a[],
                         const PetscScalar a_x[], const PetscScalar a_t[],
                         PetscReal t,
                         const PetscReal x[], PetscInt numConstants,
                         const PetscScalar constants[], PetscScalar f1[]) {
    (void)Nf; (void)NfAux;
    (void)u_t;
    (void)aOff; (void)aOff_x; (void)a; (void)a_x; (void)a_t;
    (void)t; (void)x;

    const PetscScalar Sw = u[uOff[1]];
    const PetscScalar Sg = u[uOff[2]];

    const PetscScalar kx = getConst(numConstants, constants, 1, 100e-15);
    const PetscScalar ky = getConst(numConstants, constants, 2, 100e-15);
    const PetscScalar kz = getConst(numConstants, constants, 3, 10e-15);

    const PetscScalar mu_w = getConst(numConstants, constants, 7, 1e-3);
    const PetscScalar mu_o = getConst(numConstants, constants, 8, 5e-3);
    const PetscScalar mu_g = getConst(numConstants, constants, 9, 1e-5);

    const PetscScalar Swr = getConst(numConstants, constants, 10, 0.2);
    const PetscScalar Sor = getConst(numConstants, constants, 11, 0.2);
    const PetscScalar Sgr = getConst(numConstants, constants, 12, 0.05);

    const PetscScalar nw = getConst(numConstants, constants, 13, 2.0);
    const PetscScalar no = getConst(numConstants, constants, 14, 2.0);
    const PetscScalar ng = getConst(numConstants, constants, 15, 2.0);

    const PetscScalar krw0 = getConst(numConstants, constants, 16, 0.5);
    const PetscScalar kro0 = getConst(numConstants, constants, 17, 1.0);
    const PetscScalar krg0 = getConst(numConstants, constants, 18, 0.8);

    PetscScalar krw, kro, krg;
    blackOilRelPermCorey(Sw, Sg, Swr, Sor, Sgr, nw, no, ng, krw0, kro0, krg0, krw, kro, krg);

    const PetscScalar lambda_t = (krw / mu_w) + (kro / mu_o) + (krg / mu_g);

    for (PetscInt d = 0; d < dim; ++d) {
        PetscScalar kd = kx;
        if (dim == 2) {
            kd = (d == 0) ? kx : kz; // common 2D x-z convention
        } else if (dim == 3) {
            kd = (d == 0) ? kx : (d == 1 ? ky : kz);
        }
        f1[d] = kd * lambda_t * u_x[uOff_x[0] + d];
    }
}

void g0_BlackOilPressurePressure(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                                 const PetscInt uOff[], const PetscInt uOff_x[],
                                 const PetscScalar u[], const PetscScalar u_t[],
                                 const PetscScalar u_x[], const PetscInt aOff[],
                                 const PetscInt aOff_x[], const PetscScalar a[],
                                 const PetscScalar a_x[], const PetscScalar a_t[],
                                 PetscReal t, PetscReal u_tShift,
                                 const PetscReal x[], PetscInt numConstants,
                                 const PetscScalar constants[], PetscScalar g0[]) {
    (void)dim; (void)Nf; (void)NfAux;
    (void)uOff_x; (void)u_t; (void)u_x;
    (void)aOff; (void)aOff_x; (void)a; (void)a_x; (void)a_t;
    (void)t; (void)x;

    const PetscScalar Sw = u[uOff[1]];
    const PetscScalar Sg = u[uOff[2]];
    PetscScalar So = 1.0 - Sw - Sg;
    if (So < 0.0) So = 0.0;

    const PetscScalar phi = getConst(numConstants, constants, 0, 0.2);
    const PetscScalar cw  = getConst(numConstants, constants, 4, 4.5e-10);
    const PetscScalar co  = getConst(numConstants, constants, 5, 1.5e-9);
    const PetscScalar cg  = getConst(numConstants, constants, 6, 1.0e-8);

    const PetscScalar ct_eff = Sw * cw + So * co + Sg * cg;
    g0[0] = phi * ct_eff * u_tShift;
}

void g3_BlackOilPressurePressure(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                                 const PetscInt uOff[], const PetscInt uOff_x[],
                                 const PetscScalar u[], const PetscScalar u_t[],
                                 const PetscScalar u_x[], const PetscInt aOff[],
                                 const PetscInt aOff_x[], const PetscScalar a[],
                                 const PetscScalar a_x[], const PetscScalar a_t[],
                                 PetscReal t, PetscReal u_tShift,
                                 const PetscReal x[], PetscInt numConstants,
                                 const PetscScalar constants[], PetscScalar g3[]) {
    (void)Nf; (void)NfAux;
    (void)uOff_x; (void)u_t; (void)u_x;
    (void)aOff; (void)aOff_x; (void)a; (void)a_x; (void)a_t;
    (void)t; (void)u_tShift; (void)x;

    const PetscScalar Sw = u[uOff[1]];
    const PetscScalar Sg = u[uOff[2]];

    const PetscScalar kx = getConst(numConstants, constants, 1, 100e-15);
    const PetscScalar ky = getConst(numConstants, constants, 2, 100e-15);
    const PetscScalar kz = getConst(numConstants, constants, 3, 10e-15);

    const PetscScalar mu_w = getConst(numConstants, constants, 7, 1e-3);
    const PetscScalar mu_o = getConst(numConstants, constants, 8, 5e-3);
    const PetscScalar mu_g = getConst(numConstants, constants, 9, 1e-5);

    const PetscScalar Swr = getConst(numConstants, constants, 10, 0.2);
    const PetscScalar Sor = getConst(numConstants, constants, 11, 0.2);
    const PetscScalar Sgr = getConst(numConstants, constants, 12, 0.05);

    const PetscScalar nw = getConst(numConstants, constants, 13, 2.0);
    const PetscScalar no = getConst(numConstants, constants, 14, 2.0);
    const PetscScalar ng = getConst(numConstants, constants, 15, 2.0);

    const PetscScalar krw0 = getConst(numConstants, constants, 16, 0.5);
    const PetscScalar kro0 = getConst(numConstants, constants, 17, 1.0);
    const PetscScalar krg0 = getConst(numConstants, constants, 18, 0.8);

    PetscScalar krw, kro, krg;
    blackOilRelPermCorey(Sw, Sg, Swr, Sor, Sgr, nw, no, ng, krw0, kro0, krg0, krw, kro, krg);
    const PetscScalar lambda_t = (krw / mu_w) + (kro / mu_o) + (krg / mu_g);

    const PetscInt n = dim * dim;
    for (PetscInt i = 0; i < n; ++i) g3[i] = 0.0;

    if (dim == 2) {
        g3[0] = kx * lambda_t; // xx
        g3[3] = kz * lambda_t; // zz (2D x-z convention)
    } else if (dim == 3) {
        g3[0] = kx * lambda_t; // xx
        g3[4] = ky * lambda_t; // yy
        g3[8] = kz * lambda_t; // zz
    } else {
        for (PetscInt d = 0; d < dim; ++d) {
            g3[d * dim + d] = kx * lambda_t;
        }
    }
}

// -----------------------------------------------------------------------------
// Placeholder saturation equations (safe assembly)
// -----------------------------------------------------------------------------

void f0_BlackOilSw(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                   const PetscInt uOff[], const PetscInt uOff_x[],
                   const PetscScalar u[], const PetscScalar u_t[],
                   const PetscScalar u_x[], const PetscInt aOff[],
                   const PetscInt aOff_x[], const PetscScalar a[],
                   const PetscScalar a_x[], const PetscScalar a_t[],
                   PetscReal t,
                   const PetscReal x[], PetscInt numConstants,
                   const PetscScalar constants[], PetscScalar f0[]) {
    (void)dim; (void)Nf; (void)NfAux;
    (void)uOff_x; (void)u; (void)u_x;
    (void)aOff; (void)aOff_x; (void)a; (void)a_x; (void)a_t;
    (void)t; (void)x; (void)numConstants; (void)constants;
    f0[0] = u_t[uOff[1]];
}

void f0_BlackOilSg(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                   const PetscInt uOff[], const PetscInt uOff_x[],
                   const PetscScalar u[], const PetscScalar u_t[],
                   const PetscScalar u_x[], const PetscInt aOff[],
                   const PetscInt aOff_x[], const PetscScalar a[],
                   const PetscScalar a_x[], const PetscScalar a_t[],
                   PetscReal t,
                   const PetscReal x[], PetscInt numConstants,
                   const PetscScalar constants[], PetscScalar f0[]) {
    (void)dim; (void)Nf; (void)NfAux;
    (void)uOff_x; (void)u; (void)u_x;
    (void)aOff; (void)aOff_x; (void)a; (void)a_x; (void)a_t;
    (void)t; (void)x; (void)numConstants; (void)constants;
    f0[0] = u_t[uOff[2]];
}

void f1_Zero(PetscInt dim, PetscInt Nf, PetscInt NfAux,
             const PetscInt uOff[], const PetscInt uOff_x[],
             const PetscScalar u[], const PetscScalar u_t[],
             const PetscScalar u_x[], const PetscInt aOff[],
             const PetscInt aOff_x[], const PetscScalar a[],
             const PetscScalar a_x[], const PetscScalar a_t[],
             PetscReal t,
             const PetscReal x[], PetscInt numConstants,
             const PetscScalar constants[], PetscScalar f1[]) {
    (void)Nf; (void)NfAux;
    (void)uOff; (void)uOff_x; (void)u; (void)u_t; (void)u_x;
    (void)aOff; (void)aOff_x; (void)a; (void)a_x; (void)a_t;
    (void)t; (void)x; (void)numConstants; (void)constants;
    for (PetscInt d = 0; d < dim; ++d) f1[d] = 0.0;
}

void g0_Identity(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                 const PetscInt uOff[], const PetscInt uOff_x[],
                 const PetscScalar u[], const PetscScalar u_t[],
                 const PetscScalar u_x[], const PetscInt aOff[],
                 const PetscInt aOff_x[], const PetscScalar a[],
                 const PetscScalar a_x[], const PetscScalar a_t[],
                 PetscReal t, PetscReal u_tShift,
                 const PetscReal x[], PetscInt numConstants,
                 const PetscScalar constants[], PetscScalar g0[]) {
    (void)dim; (void)Nf; (void)NfAux;
    (void)uOff; (void)uOff_x; (void)u; (void)u_t; (void)u_x;
    (void)aOff; (void)aOff_x; (void)a; (void)a_x; (void)a_t;
    (void)t; (void)x; (void)numConstants; (void)constants;
    g0[0] = u_tShift;
}

} // namespace PetscFEFluidFlow
} // namespace FSRM

