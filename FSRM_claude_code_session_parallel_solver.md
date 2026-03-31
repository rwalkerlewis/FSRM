# FSRM: Fix Constants Layout, Real Validation, End-to-End Solve

Read CLAUDE.md first. Build and test everything in Docker.

## Bug 1 (CRITICAL): Constants Layout Mismatch

Three callback files read from the PetscDS constants array at different indices for the same physical quantity. The Simulator fills one layout, PetscFEPoroelasticity expects another.

### What Simulator fills (25 elements):
```
[0]  = lambda           [1]  = mu              [2]  = rho_solid
[3]  = porosity         [4]  = perm_x (m^2)    [5]  = perm_y (m^2)
[6]  = perm_z (m^2)     [7]  = cw              [8]  = co
[9]  = cg               [10] = mu_w            [11] = mu_o
[12] = mu_g             [13] = Swr             [14] = Sor
[15] = Sgr              [16] = nw              [17] = no
[18] = ng               [19] = krw0            [20] = kro0
[21] = krg0             [22] = biot_alpha      [23] = 1/M
[24] = rho_fluid
```

### What PetscFEPoroelasticity EXPECTS to read:
```
[0] = lambda  [1] = mu  [4] = permeability  [5] = fluid_viscosity
[7] = biot_coefficient  [8] = biot_modulus_inv
```

### The mismatch:
- Poroelasticity reads constants[5] as fluid_viscosity -> gets permeability_y
- Poroelasticity reads constants[7] as biot_alpha -> gets water_compressibility
- Poroelasticity reads constants[8] as 1/M -> gets oil_compressibility
- biot_alpha is actually at [22], 1/M at [23], fluid_viscosity at [10]

### What PetscFEFluidFlow BlackOil EXPECTS to read (OLD layout, never updated):
```
[0] = porosity  [1] = kx  [2] = ky  [3] = kz  [4] = cw  [5] = co ...
```
- FluidFlow reads constants[0] as porosity -> gets lambda (Lame parameter!)
- FluidFlow reads constants[1] as kx -> gets mu (shear modulus!)

### Fix: Pick ONE layout. Update ALL callback files to match.

The Simulator layout is the most complete. Keep it. Update the callback files:

**PetscFEPoroelasticity.cpp** -- change these reads:
```cpp
// OLD (wrong):
constants[4]  // permeability -> actually perm_x, OK for isotropic
constants[5]  // fluid_viscosity -> actually perm_y, WRONG
constants[7]  // biot_alpha -> actually cw, WRONG
constants[8]  // biot_modulus_inv -> actually co, WRONG

// NEW (correct per Simulator layout):
constants[4]  // perm_x (isotropic permeability, use this)
constants[10] // mu_w (water viscosity = fluid viscosity for single-phase)
constants[22] // biot_alpha
constants[23] // 1/M (biot modulus inverse)
```

Update EVERY function in PetscFEPoroelasticity.cpp. There are 9 functions that read from constants. Check each one.

Also update the header comments in PetscFEPoroelasticity.hpp to match the actual Simulator layout.

**PetscFEFluidFlow.cpp** -- change the BlackOil callback reads:
```cpp
// OLD (wrong indices):
getConst(numConstants, constants, 0, 0.2);   // phi -> actually lambda
getConst(numConstants, constants, 1, 100e-15); // kx -> actually mu
getConst(numConstants, constants, 2, 100e-15); // ky -> actually rho_solid
// ... etc

// NEW (correct per Simulator layout):
getConst(numConstants, constants, 3, 0.2);     // phi = porosity
getConst(numConstants, constants, 4, 100e-15); // kx = permeability_x
getConst(numConstants, constants, 5, 100e-15); // ky = permeability_y
getConst(numConstants, constants, 6, 10e-15);  // kz = permeability_z
getConst(numConstants, constants, 7, 4.5e-10); // cw = water_compressibility
getConst(numConstants, constants, 8, 1.5e-9);  // co = oil_compressibility
getConst(numConstants, constants, 9, 1e-8);    // cg = gas_compressibility
getConst(numConstants, constants, 10, 1e-3);   // mu_w = water_viscosity
getConst(numConstants, constants, 11, 5e-3);   // mu_o = oil_viscosity
getConst(numConstants, constants, 12, 1e-5);   // mu_g = gas_viscosity
getConst(numConstants, constants, 13, 0.2);    // Swr
getConst(numConstants, constants, 14, 0.2);    // Sor
getConst(numConstants, constants, 15, 0.05);   // Sgr
getConst(numConstants, constants, 16, 2.0);    // nw
getConst(numConstants, constants, 17, 2.0);    // no
getConst(numConstants, constants, 18, 2.0);    // ng
getConst(numConstants, constants, 19, 0.5);    // krw0
getConst(numConstants, constants, 20, 1.0);    // kro0
getConst(numConstants, constants, 21, 0.8);    // krg0
```

Update ALL BlackOil callbacks: f0_BlackOilPressure, f1_BlackOilPressure, g0_BlackOilPressurePressure, g3_BlackOilPressurePressure, f0_BlackOilSw, f0_BlackOilSg.

**PetscFEFluidFlow.cpp** -- fix f0_SinglePhase to read from constants instead of hardcoding:
```cpp
// OLD (hardcoded):
const PetscReal phi = 0.2;
const PetscReal ct  = 1e-9;

// NEW:
const PetscReal phi = (numConstants > 3) ? PetscRealPart(constants[3]) : 0.2;
const PetscReal ct  = (numConstants > 7) ? PetscRealPart(constants[7]) : 1e-9;
```

**PetscFEElasticity.cpp** -- already reads [0]=lambda, [1]=mu, [2]=rho. No changes needed.

### Verify after fixing:
All three callback files must agree on the layout documented in Simulator::setupPhysics().
Create a static_assert or compile-time constant header if you want, but at minimum add this
comment block to the top of EACH callback .cpp file:
```cpp
// UNIFIED PetscDS CONSTANTS LAYOUT (set once in Simulator::setupPhysics):
// [0]=lambda [1]=mu [2]=rho_s [3]=phi [4]=kx [5]=ky [6]=kz
// [7]=cw [8]=co [9]=cg [10]=mu_w [11]=mu_o [12]=mu_g
// [13]=Swr [14]=Sor [15]=Sgr [16]=nw [17]=no [18]=ng
// [19]=krw0 [20]=kro0 [21]=krg0 [22]=biot_alpha [23]=1/M [24]=rho_f
```


## Bug 2: Validation Tests Don't Test Anything

### test_elastostatics_patch.cpp
This test computes analytical Hooke's law, calls the OLD GeomechanicsKernel::residual()
(strong-form placeholder), then GTEST_SKIPs. It never calls PetscFEElasticity functions.

Rewrite to directly call PetscFEElasticity::f1_elastostatics with known inputs and verify
the output stress matches Hooke's law:

```cpp
TEST_F(ElastostaticsPatchTest, F1ProducesCorrectStress) {
    const double E = 10e9, nu = 0.25;
    const double lambda = E*nu/((1+nu)*(1-2*nu));
    const double mu = E/(2*(1+nu));
    const double eps_0 = 0.001;

    PetscScalar constants[25] = {0};
    constants[0] = lambda; constants[1] = mu; constants[2] = 2700;

    // Displacement gradient: uniform uniaxial strain eps_xx = eps_0
    PetscScalar u_x[9] = {0};
    u_x[0] = eps_0;  // du_x/dx

    PetscInt uOff[1] = {0}, uOff_x[1] = {0};
    PetscScalar f1[9] = {0};
    PetscReal x[3] = {0.5, 0.5, 0.5};

    PetscFEElasticity::f1_elastostatics(
        3, 1, 0, uOff, uOff_x, nullptr, nullptr, u_x,
        nullptr, nullptr, nullptr, nullptr, nullptr,
        0.0, x, 25, constants, f1);

    // f1[c*dim+d] = sigma_{cd}
    EXPECT_NEAR(PetscRealPart(f1[0]), (lambda+2*mu)*eps_0, 1e-6); // sigma_xx
    EXPECT_NEAR(PetscRealPart(f1[4]), lambda*eps_0, 1e-6);         // sigma_yy
    EXPECT_NEAR(PetscRealPart(f1[8]), lambda*eps_0, 1e-6);         // sigma_zz
    EXPECT_NEAR(PetscRealPart(f1[1]), 0.0, 1e-12);                 // sigma_xy
}
```

Also test g3_elastostatics to verify the stiffness tensor entries.

### test_terzaghi.cpp
Currently GTEST_SKIPs saying PetscFEPoroelasticity doesn't exist. It does now.

Add a test that calls PetscFEPoroelasticity callbacks directly with known inputs:

```cpp
TEST_F(TerzaghiConsolidationTest, F1DisplacementBiotCoupling) {
    // Verify that f1_displacement includes -alpha*p*delta_{cd} term
    PetscScalar constants[25] = {0};
    constants[0] = 4e9;   // lambda
    constants[1] = 4e9;   // mu
    constants[22] = 0.8;  // biot alpha

    PetscScalar u[4] = {0};  // [pressure, ux, uy, uz]
    u[0] = 1e6;  // 1 MPa pore pressure

    PetscInt uOff[2] = {0, 1};
    PetscInt uOff_x[2] = {0, 3};
    PetscScalar u_x[12] = {0};  // 3 for pressure grad + 9 for displacement grad
    PetscScalar f1[9] = {0};
    PetscReal x[3] = {0.5, 0.5, 0.5};

    PetscFEPoroelasticity::f1_displacement(
        3, 2, 0, uOff, uOff_x, u, nullptr, u_x,
        nullptr, nullptr, nullptr, nullptr, nullptr,
        0.0, x, 25, constants, f1);

    // With zero strain, f1 should be just -alpha*p*delta
    double alpha = 0.8, p = 1e6;
    EXPECT_NEAR(PetscRealPart(f1[0]), -alpha*p, 1e-3);  // sigma_xx - alpha*p
    EXPECT_NEAR(PetscRealPart(f1[4]), -alpha*p, 1e-3);  // sigma_yy - alpha*p
    EXPECT_NEAR(PetscRealPart(f1[8]), -alpha*p, 1e-3);  // sigma_zz - alpha*p
}
```

Keep the existing analytical Terzaghi series tests -- they're fine as reference value checks. Just remove the GTEST_SKIP that claims poroelasticity doesn't exist, and add callback-level tests alongside.

### Also update test_mms_elasticity.cpp and test_mms_wave_propagation.cpp:
If these still GTEST_SKIP because they test GeomechanicsKernel::residual() (the old strong-form), consider adding parallel tests that exercise PetscFEElasticity callbacks directly, like the patch test above.


## Bug 3: End-to-End Verification

After fixing constants, run the Simulator and verify it produces nonzero output:

```bash
docker run --rm -v $(pwd):/workspace -w /workspace fsrm-ci:local bash -c '
  cd build && cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_TESTING=ON -DENABLE_CUDA=OFF -DBUILD_EXAMPLES=ON &&
  make -j$(nproc) &&
  echo "=== ELASTOSTATICS ===" &&
  ./fsrm -c ../config/test_elastostatics.config 2>&1 | tail -20 &&
  echo "=== ELASTODYNAMICS ===" &&
  ./fsrm -c ../config/test_elastodynamics.config 2>&1 | tail -20 &&
  echo "=== TESTS ===" &&
  ctest --output-on-failure 2>&1 | tail -30
'
```

If the Simulator segfaults or produces NaN, debug by checking:
1. Is `prob` (PetscDS) actually created before PetscDSSetResidual? Check `DMGetDS(dm, &prob)`.
2. Is the displacement field actually added to DM? Check setupFields() for the POROELASTIC / ELASTIC cases.
3. Are boundary conditions set? Without Dirichlet BCs, the elastostatic system is singular.
4. Is the solution vector created? Check `DMCreateGlobalVector(dm, &solution)`.

## Parallel Subagent Tasks

### Subagent A: Fix PetscFEPoroelasticity.cpp constants indices
Update all 9 functions to read from the Simulator's actual layout:
- constants[10] for fluid_viscosity (was [5])
- constants[22] for biot_alpha (was [7])
- constants[23] for 1/M (was [8])

Also update the header doc comments in PetscFEPoroelasticity.hpp.

### Subagent B: Fix PetscFEFluidFlow.cpp constants indices
Update all BlackOil callbacks to use unified layout indices [3]-[21].
Fix f0_SinglePhase to read phi and ct from constants instead of hardcoding.

### Subagent C: Rewrite validation tests
Rewrite test_elastostatics_patch.cpp to call PetscFEElasticity callbacks directly.
Rewrite test_terzaghi.cpp to call PetscFEPoroelasticity callbacks directly.
Remove stale GTEST_SKIPs.

### Main agent after subagents:
1. Build and run all tests
2. Run `./fsrm -c config/test_elastostatics.config` and verify nonzero output
3. Fix any runtime issues (segfaults, NaN, singular systems)

## Do NOT
- Modify PetscFEElasticity.hpp/.cpp (indices are already correct)
- Modify PoroelasticSolver (standalone, not affected)
- Change the Simulator constants array itself (the layout is fine, callbacks need to match it)
- Delete or gut any tests