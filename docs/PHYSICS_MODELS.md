# FSRM Physics Models

Mathematical formulations and theory behind FSRM physics models.

---

## Table of Contents

1. [Governing Equations](#governing-equations)
2. [Fluid Flow Models](#fluid-flow-models)
3. [Geomechanics Models](#geomechanics-models)
4. [Thermal Models](#thermal-models)
5. [Fracture Models](#fracture-models)
6. [Fault Mechanics](#fault-mechanics)
7. [Wave Propagation](#wave-propagation)
8. [Coupled Physics](#coupled-physics)

---

## Governing Equations

FSRM solves coupled PDEs for mass, momentum, and energy conservation.

### Mass Conservation (Fluid Flow)

$$\frac{\partial (\phi \rho_f)}{\partial t} + \nabla \cdot (\rho_f \mathbf{v}) = q$$

where:
- φ = porosity
- ρ_f = fluid density
- **v** = Darcy velocity
- q = source/sink term

### Darcy's Law

$$\mathbf{v} = -\frac{\mathbf{k}}{\mu}(\nabla p - \rho_f \mathbf{g})$$

where:
- **k** = permeability tensor
- μ = dynamic viscosity
- p = pressure
- **g** = gravity vector

### Momentum Conservation (Geomechanics)

$$\nabla \cdot \boldsymbol{\sigma} + \rho_b \mathbf{g} = \rho_b \frac{\partial^2 \mathbf{u}}{\partial t^2}$$

where:
- **σ** = Cauchy stress tensor
- ρ_b = bulk density
- **u** = displacement vector

### Energy Conservation (Thermal)

$$(\rho c_p)_{\text{eff}} \frac{\partial T}{\partial t} + \rho_f c_{p,f} \mathbf{v} \cdot \nabla T = \nabla \cdot (k_T \nabla T) + Q$$

where:
- T = temperature
- c_p = specific heat
- k_T = thermal conductivity
- Q = heat source

---

## Fluid Flow Models

### Single-Phase Flow

Compressible single-phase flow with pressure-dependent density:

$$\rho(p) = \rho_0 \exp\left[c_f(p - p_0)\right]$$

where c_f is fluid compressibility.

### Black Oil Model

Three-phase (oil, water, gas) flow with solution gas:

**Oil phase:**
$$\frac{\partial}{\partial t}\left(\phi \frac{S_o}{B_o}\right) + \nabla \cdot \left(\frac{\mathbf{k} k_{ro}}{B_o \mu_o}\nabla \Phi_o\right) = q_o$$

**Water phase:**
$$\frac{\partial}{\partial t}\left(\phi \frac{S_w}{B_w}\right) + \nabla \cdot \left(\frac{\mathbf{k} k_{rw}}{B_w \mu_w}\nabla \Phi_w\right) = q_w$$

**Gas phase:**
$$\frac{\partial}{\partial t}\left(\phi \left[\frac{S_g}{B_g} + \frac{R_s S_o}{B_o}\right]\right) + \nabla \cdot \left(\frac{\mathbf{k} k_{rg}}{B_g \mu_g}\nabla \Phi_g + \frac{R_s \mathbf{k} k_{ro}}{B_o \mu_o}\nabla \Phi_o\right) = q_g$$

where:
- S_α = saturation of phase α
- B_α = formation volume factor
- k_rα = relative permeability
- R_s = solution gas-oil ratio
- Φ_α = phase potential

**PVT Correlations:**

*Standing (1947) solution GOR:*
$$R_s = \gamma_g \left[\left(\frac{p}{18.2} + 1.4\right) 10^{0.0125 \text{API} - 0.00091 T}\right]^{1.2048}$$

*Vazquez-Beggs (1980) oil FVF:*
$$B_o = 1 + C_1 R_s + C_2(T - 60)\left(\frac{\text{API}}{\gamma_{gs}}\right) + C_3 R_s (T-60)\left(\frac{\text{API}}{\gamma_{gs}}\right)$$

### Compositional Model

Multi-component flow with phase equilibrium:

$$\frac{\partial}{\partial t}(\phi z_i \rho) + \nabla \cdot \sum_\alpha (x_{i,\alpha} \rho_\alpha \mathbf{v}_\alpha) = q_i$$

where z_i is overall mole fraction of component i.

**Peng-Robinson EOS:**
$$p = \frac{RT}{V_m - b} - \frac{a \alpha(T)}{V_m(V_m + b) + b(V_m - b)}$$

where:
- a, b = EOS parameters from critical properties
- α(T) = temperature-dependent function using acentric factor

**Flash Calculations:**
K-value approach with Rachford-Rice equation:
$$\sum_i \frac{z_i(K_i - 1)}{1 + V(K_i - 1)} = 0$$

---

## Geomechanics Models

### Linear Elasticity

Hooke's law (isotropic):
$$\sigma_{ij} = \lambda \epsilon_{kk} \delta_{ij} + 2G \epsilon_{ij}$$

where:
- λ = Lamé's first parameter
- G = shear modulus
- ε_ij = strain tensor

**Relations:**
$$\lambda = \frac{E\nu}{(1+\nu)(1-2\nu)}, \quad G = \frac{E}{2(1+\nu)}$$

### Poroelasticity (Biot Theory)

Coupled deformation-fluid flow:

**Total stress:**
$$\sigma_{ij} = \sigma'_{ij} - \alpha p \delta_{ij}$$

**Effective stress-strain:**
$$\sigma'_{ij} = \lambda \epsilon_{kk} \delta_{ij} + 2G \epsilon_{ij}$$

**Fluid content:**
$$\zeta = \alpha \epsilon_{kk} + \frac{p}{M}$$

where:
- α = Biot coefficient (α = 1 - K_d/K_s)
- M = Biot modulus
- K_d = drained bulk modulus
- K_s = solid grain bulk modulus

**Coupled equations:**
$$\nabla \cdot \boldsymbol{\sigma}' - \alpha \nabla p + \rho_b \mathbf{g} = 0$$
$$\frac{\partial \zeta}{\partial t} + \nabla \cdot \mathbf{v} = q$$

### Viscoelasticity

**Maxwell model:**
$$\dot{\epsilon} = \frac{\dot{\sigma}}{E} + \frac{\sigma}{\eta}$$

Relaxation time: τ = η/E

**Kelvin-Voigt model:**
$$\sigma = E\epsilon + \eta\dot{\epsilon}$$

**Standard Linear Solid (SLS):**
$$\sigma + \frac{\eta}{E_1}\dot{\sigma} = E_2\epsilon + \eta\left(1 + \frac{E_2}{E_1}\right)\dot{\epsilon}$$

### Elastoplasticity

**Yield functions:**

*Mohr-Coulomb:*
$$f = \tau - c - \sigma_n \tan\phi$$

or in terms of principal stresses:
$$f = \frac{1}{2}(\sigma_1 - \sigma_3) + \frac{1}{2}(\sigma_1 + \sigma_3)\sin\phi - c\cos\phi$$

*Drucker-Prager:*
$$f = \sqrt{J_2} - a I_1 - k$$

where:
- J_2 = second deviatoric stress invariant
- I_1 = first stress invariant
- a, k = material parameters

**Flow rule:**
$$\dot{\epsilon}^p_{ij} = \dot{\lambda} \frac{\partial g}{\partial \sigma_{ij}}$$

### Anisotropy (VTI/HTI)

Transverse isotropy stiffness matrix:

$$\begin{bmatrix} \sigma_{11} \\ \sigma_{22} \\ \sigma_{33} \\ \sigma_{23} \\ \sigma_{13} \\ \sigma_{12} \end{bmatrix} = \begin{bmatrix} C_{11} & C_{12} & C_{13} & 0 & 0 & 0 \\ C_{12} & C_{11} & C_{13} & 0 & 0 & 0 \\ C_{13} & C_{13} & C_{33} & 0 & 0 & 0 \\ 0 & 0 & 0 & C_{44} & 0 & 0 \\ 0 & 0 & 0 & 0 & C_{44} & 0 \\ 0 & 0 & 0 & 0 & 0 & C_{66} \end{bmatrix} \begin{bmatrix} \epsilon_{11} \\ \epsilon_{22} \\ \epsilon_{33} \\ 2\epsilon_{23} \\ 2\epsilon_{13} \\ 2\epsilon_{12} \end{bmatrix}$$

with C_66 = (C_11 - C_12)/2.

**Thomsen parameters:**
$$\epsilon_T = \frac{C_{11} - C_{33}}{2C_{33}}, \quad \delta_T = \frac{(C_{13} + C_{44})^2 - (C_{33} - C_{44})^2}{2C_{33}(C_{33} - C_{44})}$$

---

## Thermal Models

### Heat Conduction

$$(\rho c_p)_{\text{eff}} \frac{\partial T}{\partial t} = \nabla \cdot (k_T \nabla T) + Q$$

Effective properties:
$$(\rho c_p)_{\text{eff}} = (1-\phi)\rho_s c_{p,s} + \phi \rho_f c_{p,f}$$
$$k_{T,\text{eff}} = (1-\phi) k_{T,s} + \phi k_{T,f}$$

### Thermal Stress

$$\sigma_{ij} = C_{ijkl}(\epsilon_{kl} - \alpha_T \Delta T \delta_{kl})$$

---

## Fracture Models

### Linear Elastic Fracture Mechanics (LEFM)

**Stress intensity factors:**
$$K_I = Y \sigma \sqrt{\pi a}$$

where:
- Y = geometry factor
- σ = remote stress
- a = crack half-length

**Propagation criterion:**
$$K_I \geq K_{IC}$$

**Energy release rate:**
$$G = \frac{K_I^2}{E'}, \quad E' = E \text{ (plane stress)}, \quad E' = \frac{E}{1-\nu^2} \text{ (plane strain)}$$

### Fracture Permeability

**Cubic law:**
$$k_f = \frac{w^2}{12}$$

where w is aperture.

### Proppant Transport

**Advection-settling:**
$$\frac{\partial (\phi_f c)}{\partial t} + \nabla \cdot (c \mathbf{v}_f) - \nabla \cdot (c v_s \hat{\mathbf{z}}) = 0$$

Settling velocity (Stokes):
$$v_s = \frac{(\rho_p - \rho_f) g d_p^2}{18 \mu}$$

---

## Fault Mechanics

### Coulomb Friction

**Failure criterion:**
$$|\tau| \geq \mu_s (\sigma_n - p) + c$$

**Slip condition:**
$$|\tau| = \mu_d (\sigma_n - p)$$

### Rate-State Friction

**Dieterich-Ruina formulation:**
$$\mu = \mu_0 + a \ln\left(\frac{V}{V_0}\right) + b \ln\left(\frac{V_0 \theta}{D_c}\right)$$

**State evolution laws:**

*Aging (Dieterich):*
$$\frac{d\theta}{dt} = 1 - \frac{V\theta}{D_c}$$

*Slip (Ruina):*
$$\frac{d\theta}{dt} = -\frac{V\theta}{D_c} \ln\left(\frac{V\theta}{D_c}\right)$$

where:
- V = slip velocity
- θ = state variable
- D_c = critical slip distance
- a, b = constitutive parameters

**Stability:**
- a - b > 0: velocity strengthening (stable)
- a - b < 0: velocity weakening (potentially unstable)

### Stress Resolution on Faults

For fault with normal **n** and strike direction **s**:
$$\sigma_n = \mathbf{n} \cdot \boldsymbol{\sigma} \cdot \mathbf{n}$$
$$\tau = \mathbf{s} \cdot \boldsymbol{\sigma} \cdot \mathbf{n}$$

### Coulomb Stress Transfer

$$\Delta CFS = \Delta\tau + \mu(\Delta\sigma_n - \Delta p)$$

### Seismicity Analysis

**Gutenberg-Richter relation:**
$$\log_{10} N = a - b M$$

**Seismic moment:**
$$M_0 = G A D$$

where:
- G = shear modulus
- A = rupture area
- D = average slip

**Moment magnitude:**
$$M_w = \frac{2}{3}\log_{10} M_0 - 6.07$$

**Omori law (aftershocks):**
$$n(t) = \frac{K}{(c + t)^p}$$

---

## Wave Propagation

### Elastic Wave Equation

$$\rho \frac{\partial^2 \mathbf{u}}{\partial t^2} = \nabla \cdot \boldsymbol{\sigma} + \mathbf{f}$$

**P-wave velocity:**
$$V_p = \sqrt{\frac{\lambda + 2G}{\rho}} = \sqrt{\frac{K + \frac{4}{3}G}{\rho}}$$

**S-wave velocity:**
$$V_s = \sqrt{\frac{G}{\rho}}$$

### Poroelastic Waves (Biot Theory)

**Fast P-wave, slow P-wave (Biot), and S-wave:**

Dispersion relations from eigenvalue problem of coupled equations.

**High-frequency limit:**
$$V_{Pf}^2 = \frac{H + C^2/M}{\rho}, \quad V_s^2 = \frac{G}{\rho}$$

### Absorbing Boundaries (PML)

Perfectly Matched Layer for wave absorption:
$$\tilde{\partial}_i = \frac{1}{1 + i\omega^{-1}d_i(x)}\partial_i$$

### Source Time Functions

**Ricker wavelet:**
$$f(t) = (1 - 2\pi^2 f_0^2 t^2) \exp(-\pi^2 f_0^2 t^2)$$

**Gaussian:**
$$f(t) = \exp\left(-\frac{(t-t_0)^2}{2\sigma^2}\right)$$

### Dynamic Permeability Enhancement

**Strain-dependent model:**
$$k = k_0 \left(1 + \beta |\epsilon_{vol}|\right)$$

**Frequency-dependent:**
$$k(\omega) = k_0 \left(1 + A \frac{\omega/\omega_c}{1 + (\omega/\omega_c)^2}\right)$$

---

## Coupled Physics

### Flow-Geomechanics Coupling

**Iterative coupling:**
1. Solve flow: p^{n+1} given σ^n
2. Update effective stress: σ'^{n+1} = f(p^{n+1})
3. Solve mechanics: u^{n+1}, σ^{n+1}
4. Update porosity/permeability
5. Iterate to convergence

**Fixed-stress split:**
$$\phi^{n+1} = \phi^n + \frac{\alpha}{K_d}(\sigma_m^{n+1} - \sigma_m^n) + \frac{\alpha^2}{K_d}(p^{n+1} - p^n)$$

### Thermo-Hydro-Mechanical (THM)

Full coupling includes:
- Thermal expansion in mechanics
- Temperature-dependent fluid properties
- Heat advection by fluid flow
- Viscous dissipation

### Static-to-Dynamic Triggering

Quasi-static stress changes from injection → fault stress analysis → rate-state nucleation → dynamic rupture propagation.

---

## Numerical Methods

### Spatial Discretization

- Finite element (FEM) for mechanics
- Finite volume (FVM) for flow
- DMPlex for all mesh types (structured and unstructured)

### Time Integration

- Implicit Euler for stability
- Adaptive timestepping based on:
  - Newton convergence
  - CFL condition (dynamics)
  - Error estimates

### Nonlinear Solvers

- Newton-Raphson (SNES)
- Line search / trust region
- Jacobian-free Newton-Krylov (JFNK)

### Linear Solvers

- Krylov methods: GMRES, CG, BiCGStab
- Preconditioners: ILU, ASM, AMG (Hypre)

---

## References

1. Biot, M.A. (1941). General theory of three-dimensional consolidation.
2. Dieterich, J.H. (1979). Modeling of rock friction.
3. Standing, M.B. (1947). A pressure-volume-temperature correlation for mixtures of California oils and gases.
4. Peng, D.Y. & Robinson, D.B. (1976). A new two-constant equation of state.
5. Thomsen, L. (1986). Weak elastic anisotropy.
