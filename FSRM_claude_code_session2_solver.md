Read CLAUDE.md for build/test instructions.

Session 1 added PetscFEElasticity callbacks for elastostatics with PetscDS 
registration in Simulator::setupPhysics(). The Simulator now solves static 
elasticity via DMPlexTSComputeIFunctionFEM.

Now add dynamics. The elastodynamic wave equation is:
  rho * d2u/dt2 = div(sigma) + f

In PETSc, second-order-in-time systems use TSSetI2Function / TSSetI2Jacobian 
with TSALPHA2 (generalized-alpha, which includes Newmark-beta as a special case).

The I2Function signature is:
  F(t, U, U_t, U_tt) = M*U_tt + C*U_t + K*U - f = 0

Where M = rho*I (mass), C = damping, K = stiffness (from elasticity f1/g3).

Steps:

1. In PetscFEElasticity, add I2-form callbacks:
   - f0_elastodynamics_I2: rho * u_tt[c] (mass times acceleration)
   - f1_elastodynamics_I2: sigma_{ij} (same as elastostatics f1, stress divergence)
   - g0_elastodynamics_I2: rho * u_ttShift * delta_{cd} (mass Jacobian)
   - g3_elastodynamics_I2: C_{ijkl} (same as elastostatics g3, stiffness Jacobian)

   The I2Function/I2Jacobian callback signatures add U_tt and u_ttShift 
   parameters. Check PETSc ts/tutorials/ex76.c for the exact signatures.

2. In Simulator::setupTimeStepper(), when elastodynamics is enabled:
   - Use TSSetType(ts, TSALPHA2) instead of TSBEULER
   - Use TSSetI2Function / TSSetI2Jacobian instead of TSSetIFunction / TSSetIJacobian
   - Remove the warning about falling back to BEULER

3. In Simulator::setupPhysics(), for elastodynamics:
   - Register PetscDSSetResidual with the f0/f1 dynamic callbacks
   - Register PetscDSSetJacobian with the g0/g3 dynamic callbacks
   - Set constants: [0]=lambda, [1]=mu, [2]=rho
   - Set use_fem_time_residual_ = true

4. Simulator::FormFunction needs a second-order variant. Add:
   static PetscErrorCode FormFunction2(TS ts, PetscReal t, Vec U, Vec U_t, 
                                        Vec U_tt, Vec F, void *ctx);
   static PetscErrorCode FormJacobian2(TS ts, PetscReal t, Vec U, Vec U_t, 
                                        Vec U_tt, PetscReal v, PetscReal a,
                                        Mat J, Mat P, void *ctx);
   These delegate to DMPlexTSComputeI2FunctionFEM / DMPlexTSComputeI2JacobianFEM.

5. Create config/test_elastodynamics.config:
   - 3D box, 10x10x10 elements
   - Granite properties (rho=2700, E=50GPa, nu=0.25)
   - Point source impulse at center (or initial displacement perturbation)
   - End time = 0.01s, dt = 1e-4s
   - Fixed boundaries on all faces

6. Verify: run the config, check that displacement is nonzero and propagates 
   outward from the source. The P-wave should reach the boundary in approximately 
   L/(2*vp) = 0.5/5500 ~ 9e-5 seconds.

7. Update tests: if test_mms_wave_propagation.cpp can now exercise the real 
   PetscDS callbacks, remove the GTEST_SKIP guard.

Do not modify PoroelasticSolver, PetscFEFluidFlow, or any fluid flow code.
Build and test in Docker.