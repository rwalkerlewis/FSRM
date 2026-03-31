Read CLAUDE.md for build/test instructions.

PoroelasticSolver (src/numerics/PoroelasticSolver.cpp) has complete, correct 
PetscDS callbacks for coupled Biot poroelasticity: pressure diffusion + 
elastic displacement with Biot coupling terms. It creates its own DM, FE 
fields, and PetscDS internally.

The problem: PoroelasticSolver is a standalone solver that bypasses the 
Simulator pipeline. The Simulator has a PoroelastodynamicsKernel but never 
registers PetscDS callbacks for it.

Make the Simulator solve Biot poroelasticity using the same physics as 
PoroelasticSolver but through the standard pipeline:

1. Create PetscFEPoroelasticity callbacks (namespace, not class methods) 
   following the PetscFEFluidFlow / PetscFEElasticity pattern. Port the 
   callback math from PoroelasticSolver's static methods.

   Decision: use a single 3-component vector displacement field (like 
   PetscFEElasticity) plus a scalar pressure field. This means 2 fields 
   total, not the 5 separate scalars that PoroelasticSolver uses. The 
   math is identical but the index arithmetic changes.

   Fields: 0 = pressure (1 component), 1 = displacement (3 components)
   
   Residuals:
   - Field 0 (pressure): f0 = (1/M)*dp/dt + alpha*div(du/dt) - source
                          f1 = (k/mu)*grad(p)
   - Field 1 (displacement): f0 = 0 (quasi-static) or rho*u_tt (dynamic)
                              f1 = sigma_ij - alpha*p*delta_ij

   Jacobians: all 4 blocks (p-p, p-u, u-p, u-u) including Biot cross terms.

2. In Simulator::setupFields(), when poroelastic is selected, create both 
   pressure (1-comp) and displacement (3-comp) PetscFE fields.

3. In Simulator::setupPhysics(), register PetscDS callbacks for both fields 
   with the correct Biot coupling Jacobian blocks.

4. Create config/test_poroelasticity.config:
   - Terzaghi consolidation: 1D column, load on top, drainage at top
   - Known analytical solution: p(z,t) = series expansion
   - This is the classic Biot validation problem

5. Add a physics validation test that checks Terzaghi consolidation against 
   the analytical solution. Pressure at mid-height at t=t_cv/4 should match 
   the series solution within 5%.

Do not modify PoroelasticSolver itself. It remains as an independent solver.
Build and test in Docker.