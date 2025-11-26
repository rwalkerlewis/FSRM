# FSRM Numerical Methods

Comprehensive documentation of the numerical approaches used to solve the governing equations in FSRM.

---

## Table of Contents

1. [Overview](#overview)
2. [Spatial Discretization](#spatial-discretization)
3. [Temporal Integration](#temporal-integration)
4. [Nonlinear Solvers](#nonlinear-solvers)
5. [Linear Solvers](#linear-solvers)
6. [Coupling Strategies](#coupling-strategies)
7. [Convergence and Stability](#convergence-and-stability)
8. [Method of Manufactured Solutions](#method-of-manufactured-solutions)
9. [Performance Considerations](#performance-considerations)

---

## Overview

FSRM uses PETSc (Portable, Extensible Toolkit for Scientific Computation) as its numerical foundation. The simulation solves coupled partial differential equations (PDEs) for:

- **Fluid flow** (Darcy's law)
- **Solid mechanics** (momentum balance)
- **Heat transfer** (energy equation)
- **Chemical transport** (advection-diffusion-reaction)

### Solution Strategy

The general approach follows these steps:

1. **Discretize** governing equations using finite element or finite volume methods
2. **Assemble** residual vector \(F(u)\) and Jacobian matrix \(J = \partial F/\partial u\)
3. **Solve** nonlinear system \(F(u) = 0\) using Newton iteration
4. **Advance** in time using implicit methods for stability

---

## Spatial Discretization

### Finite Element Method (FEM)

FSRM uses the Galerkin finite element method for both flow and mechanics.

#### Weak Form

The strong form PDE:
$$\nabla \cdot \boldsymbol{\sigma} + \mathbf{b} = \rho \ddot{\mathbf{u}}$$

Is converted to weak form by multiplying by test function \(v\) and integrating:
$$\int_\Omega \boldsymbol{\sigma} : \nabla v \, d\Omega - \int_{\Gamma_N} \mathbf{t} \cdot v \, d\Gamma = \int_\Omega (\rho \ddot{\mathbf{u}} - \mathbf{b}) \cdot v \, d\Omega$$

#### Shape Functions

FSRM supports:
- **Q1** (trilinear hex): 8 nodes, linear interpolation
- **Q2** (triquadratic hex): 27 nodes, quadratic interpolation
- **P1** (linear tet): 4 nodes, linear interpolation
- **P2** (quadratic tet): 10 nodes, quadratic interpolation

The solution is approximated as:
$$u(\mathbf{x}) = \sum_{i=1}^{n} N_i(\mathbf{x}) u_i$$

where \(N_i\) are shape functions and \(u_i\) are nodal values.

#### Numerical Integration

Gaussian quadrature is used for element integrals:
$$\int_\Omega f \, d\Omega \approx \sum_{q=1}^{n_q} w_q f(\mathbf{x}_q) |J_q|$$

where \(w_q\) are quadrature weights and \(|J_q|\) is the Jacobian determinant.

Standard integration orders:
- Q1 elements: 2×2×2 Gauss points
- Q2 elements: 3×3×3 Gauss points

### Finite Volume Method (FVM)

For conservative fluid transport, FSRM uses cell-centered finite volumes:

$$\frac{d}{dt} \int_{V_i} \phi \rho \, dV + \sum_f \mathbf{F}_f \cdot \mathbf{n}_f A_f = Q_i$$

Face fluxes \(\mathbf{F}_f\) are computed using:
- **Two-point flux approximation** (TPFA) for structured grids
- **Multi-point flux approximation** (MPFA) for unstructured grids

#### Upwind Scheme

For advection-dominated transport, upwind discretization ensures stability:
$$F_f = \begin{cases}
\phi_L u_f & \text{if } u_f > 0 \\
\phi_R u_f & \text{if } u_f < 0
\end{cases}$$

---

## Temporal Integration

### Implicit Methods

FSRM primarily uses implicit time integration for stability with large timesteps.

#### Backward Euler (First-Order Implicit)

$$\frac{u^{n+1} - u^n}{\Delta t} = f(u^{n+1}, t^{n+1})$$

This leads to a nonlinear system:
$$F(u^{n+1}) = u^{n+1} - u^n - \Delta t \, f(u^{n+1}) = 0$$

**Properties:**
- Unconditionally stable (A-stable)
- First-order accurate: \(O(\Delta t)\)
- Strong numerical damping

#### BDF2 (Second-Order Backward Differentiation)

$$\frac{3 u^{n+1} - 4 u^n + u^{n-1}}{2 \Delta t} = f(u^{n+1})$$

**Properties:**
- A-stable
- Second-order accurate: \(O(\Delta t^2)\)
- Requires two startup steps

#### Generalized-α Method

For structural dynamics (second-order in time):

$$M \ddot{u}_{n+1-\alpha_m} + C \dot{u}_{n+1-\alpha_f} + K u_{n+1-\alpha_f} = F_{n+1-\alpha_f}$$

with interpolation:
$$u_{n+1-\alpha_f} = (1 - \alpha_f) u_{n+1} + \alpha_f u_n$$

**Properties:**
- Controllable numerical damping via \(\rho_\infty\) parameter
- Second-order accurate
- Well-suited for wave propagation

### Adaptive Timestepping

Timestep is adjusted based on:

1. **Newton convergence:** Reduce \(\Delta t\) if iterations exceed threshold
2. **Truncation error:** Estimate local error and adjust
3. **CFL condition:** For explicit dynamics, \(\Delta t \leq \Delta x / V_p\)

```
if (newton_iterations < 3) then
    Δt_new = 1.2 * Δt
else if (newton_iterations > 8) then
    Δt_new = 0.5 * Δt
end if

Δt_new = min(Δt_new, Δt_max)
Δt_new = max(Δt_new, Δt_min)
```

---

## Nonlinear Solvers

### Newton-Raphson Method

For nonlinear system \(F(u) = 0\):

1. Compute residual: \(r = -F(u^k)\)
2. Compute Jacobian: \(J = \frac{\partial F}{\partial u}\Big|_{u^k}\)
3. Solve linear system: \(J \delta u = r\)
4. Update: \(u^{k+1} = u^k + \delta u\)
5. Check convergence: \(\|F(u^{k+1})\| < \epsilon\)

#### Convergence Criteria

FSRM uses both absolute and relative tolerances:

$$\|F(u^{k+1})\| < \text{atol}$$
$$\|F(u^{k+1})\| < \text{rtol} \cdot \|F(u^0)\|$$

Typical values:
- `rtol = 1e-6`
- `atol = 1e-8`

#### Line Search

To improve convergence for difficult problems:

$$u^{k+1} = u^k + \lambda \delta u$$

where \(\lambda \in (0, 1]\) is chosen to ensure sufficient decrease in residual.

### Jacobian Computation

#### Analytical Jacobian

For element \(e\), the Jacobian contribution is:
$$J_{ij}^e = \frac{\partial F_i^e}{\partial u_j}$$

For linear elasticity:
$$J^e = \int_{\Omega^e} B^T C B \, d\Omega$$

where \(B\) is the strain-displacement matrix and \(C\) is the stiffness matrix.

#### Jacobian-Free Newton-Krylov (JFNK)

For complex physics, the Jacobian-vector product is approximated:
$$J v \approx \frac{F(u + \epsilon v) - F(u)}{\epsilon}$$

This avoids explicit Jacobian storage but requires a good preconditioner.

---

## Linear Solvers

### Krylov Methods

For the linear system \(J \delta u = r\):

#### GMRES (Generalized Minimal Residual)

- **Best for:** General nonsymmetric systems
- **Memory:** Grows with iteration count (restart after 30-50 iterations)
- **Convergence:** Minimizes residual over Krylov subspace

#### CG (Conjugate Gradient)

- **Best for:** Symmetric positive definite matrices
- **Memory:** Fixed (3 vectors)
- **Convergence:** Optimal for SPD systems

#### BiCGStab (Bi-Conjugate Gradient Stabilized)

- **Best for:** Nonsymmetric systems, limited memory
- **Memory:** Fixed (8 vectors)
- **Convergence:** May have irregular residual history

### Preconditioners

Preconditioning transforms the system to improve convergence:
$$M^{-1} J \delta u = M^{-1} r$$

#### ILU (Incomplete LU Factorization)

Approximate factorization \(J \approx LU\) with controlled fill-in:
- **ILU(0):** No fill-in, cheapest
- **ILU(k):** k levels of fill-in
- **ILUT:** Threshold-based fill-in

#### ASM (Additive Schwarz Method)

Domain decomposition preconditioner:
$$M^{-1} = \sum_{i=1}^{N} R_i^T A_i^{-1} R_i$$

where \(R_i\) restricts to subdomain \(i\) and \(A_i\) is the local problem.

Options:
- `overlap = 0`: No overlap (block Jacobi)
- `overlap = 1-2`: Typical for good convergence

#### AMG (Algebraic Multigrid)

Hierarchy of coarse grids constructed algebraically:
- **Smoothing:** Gauss-Seidel or Jacobi on each level
- **Restriction:** Transfer to coarser grid
- **Coarse solve:** Direct solve or recursive AMG
- **Prolongation:** Interpolate correction to finer grid

Recommended via Hypre/BoomerAMG:
```
-pc_type hypre -pc_hypre_type boomeramg
```

---

## Coupling Strategies

### Full Coupling (Monolithic)

All physics solved simultaneously:

$$\begin{bmatrix} A_{ff} & A_{fm} \\ A_{mf} & A_{mm} \end{bmatrix} \begin{bmatrix} \delta p \\ \delta u \end{bmatrix} = \begin{bmatrix} r_f \\ r_m \end{bmatrix}$$

**Advantages:**
- Optimal convergence for strongly coupled physics
- No iteration between physics

**Disadvantages:**
- Large memory requirement
- Complex implementation

### Sequential Coupling (Operator Splitting)

Physics solved separately with data exchange:

1. Solve flow: \(p^{n+1}\) given \(\sigma^n\)
2. Solve mechanics: \(u^{n+1}, \sigma^{n+1}\) given \(p^{n+1}\)
3. Update coupling terms
4. Iterate if needed

#### Fixed-Stress Split

For poromechanics, fix volumetric stress during flow solve:
$$\phi^{n+1} = \phi^n + \frac{\alpha}{K_d}(\sigma_m^{n+1} - \sigma_m^n) + \frac{\alpha^2}{K_d}(p^{n+1} - p^n)$$

**Properties:**
- Unconditionally stable for appropriate stabilization
- Faster than monolithic per iteration
- May require sub-iterations

### Staggered Coupling

Alternating solves without sub-iteration:
- Efficient for weakly coupled problems
- May have timestep restrictions for stability

---

## Convergence and Stability

### CFL Condition

For explicit dynamics:
$$\Delta t \leq \frac{\Delta x}{V_p}$$

where \(V_p = \sqrt{(K + 4G/3)/\rho}\) is the P-wave velocity.

### Mesh Quality

Element quality affects accuracy and solver performance:
- **Aspect ratio:** Keep close to 1 (ideally < 5)
- **Jacobian:** Must be positive everywhere
- **Angles:** Avoid very acute or obtuse angles

### Numerical Dispersion

For wave propagation, the element size limits resolved wavelengths:
$$\Delta x \leq \frac{\lambda}{10}$$

where \(\lambda\) is the minimum wavelength of interest.

---

## Method of Manufactured Solutions

FSRM includes MMS tests for verification:

### Approach

1. Choose exact solution \(u_{exact}(\mathbf{x}, t)\)
2. Substitute into governing equation to get source term \(f\)
3. Solve with source \(f\) and appropriate BCs
4. Compare numerical solution to \(u_{exact}\)

### Example: Diffusion Equation

Exact solution:
$$p_{exact} = \sin(\pi x) \sin(\pi y) e^{-2\pi^2 \kappa t}$$

Source term:
$$f = 0$$ (satisfies homogeneous equation)

Boundary conditions:
$$p = p_{exact}$$ on \(\partial\Omega$$

### Convergence Rates

Expected convergence for error \(e = u_h - u_{exact}\):

| Order | Q1/P1 (L²) | Q1/P1 (H¹) | Q2/P2 (L²) | Q2/P2 (H¹) |
|-------|-----------|-----------|-----------|-----------|
| Space | O(h²) | O(h) | O(h³) | O(h²) |
| Time (BE) | O(Δt) | O(Δt) | O(Δt) | O(Δt) |
| Time (BDF2) | O(Δt²) | O(Δt²) | O(Δt²) | O(Δt²) |

---

## Performance Considerations

### Parallel Scalability

FSRM uses MPI for distributed memory parallelism:
- Grid partitioned using ParMETIS or PETSc partitioner
- Ghost cells exchanged at partition boundaries
- Strong scaling: Fixed problem size, increase cores
- Weak scaling: Scale problem size with cores

### Optimal Process Count

Rule of thumb: 5,000-20,000 DOF per process

```bash
# Example: 1M cell mesh, 4 DOF/cell = 4M DOF
# 4M / 10000 = 400 processes optimal
mpirun -np 400 fsrm -c config.config
```

### Memory Estimation

Per-process memory requirements:
- Mesh: ~100 bytes/cell
- Solution vectors: ~8 bytes/DOF
- Jacobian (sparse): ~500-2000 bytes/DOF (depending on stencil)
- Preconditioner: varies (ILU: 2-5× Jacobian)

### Profiling

Enable PETSc logging to identify bottlenecks:
```bash
mpirun -np 4 fsrm -c config.config -log_view
```

Key metrics:
- `SNESSolve`: Total nonlinear solve time
- `KSPSolve`: Total linear solve time
- `MatMult`: Matrix-vector product time
- `PCApply`: Preconditioner application time

---

## References

1. Hughes, T.J.R. (2000). The Finite Element Method: Linear Static and Dynamic Finite Element Analysis.
2. Balay et al. PETSc Users Manual. https://petsc.org/
3. Saad, Y. (2003). Iterative Methods for Sparse Linear Systems.
4. Roache, P.J. (2002). Code Verification by the Method of Manufactured Solutions.
5. Kim, J. et al. (2011). Stability and convergence of sequential methods for coupled flow and geomechanics.
