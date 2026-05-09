# Fidelity Ladder Guide

A single-page guide to the LOW / MED / HIGH / HIGHEST tier selection across
the six fidelity ladders that govern source-physics, EOS, opacity, and time
integration in FSRM. For each ladder this document tells you what each tier
does, what it costs, and when to pick it. Cross-reference
[USER_GUIDE.md](USER_GUIDE.md) for end-to-end run examples and
[CONFIGURATION.md](CONFIGURATION.md) for the full config-key reference.

The canonical record of which gate closed at which tier under which pass is
[AXIS_1A_FIDELITY_REPORT.md](AXIS_1A_FIDELITY_REPORT.md). The per-pass
physics description is [EXPLOSION_IMPACT_PHYSICS.md](EXPLOSION_IMPACT_PHYSICS.md).

## At a glance

All six ladders are configured in `[NEAR_FIELD_SOURCE]` sub-keys (except
`time_integrator` and `time_integrator_diffusion` which are top-level
under `[TIME]`). Defaults are backwards-compatible: a config that does
not set a tier reproduces the previous pass byte-for-byte.

| Knob | LOW | MED | HIGH | HIGHEST |
|---|---|---|---|---|
| `radiation_phase` | `ZELDOVICH_RAIZER` | `MARSHAK_GREY` | `MARSHAK_GREY` + `TABULATED_PATCHED` opacity | `MARSHAK_MULTIGROUP` |
| `cavity_eos` | `IDEAL_GAS` | `TILLOTSON` | `TILLOTSON_TABULATED_PATCH` | `TABULATED_FULL` |
| `opacity_model` | `CONSTANT` | `POWER_LAW_ZR` | `TABULATED_PATCHED` | `TABULATED_FULL` |
| `operator_splitting` | `LIE` | `LIE` | `STRANG` | `STRANG_MULTIGROUP` |
| `time_integrator` (hydro) | `EXPLICIT_EULER` | `EXPLICIT_EULER` | `EXPLICIT_EULER` | `RK3_SSP` |
| `time_integrator_diffusion` | `BACKWARD_EULER` | `BACKWARD_EULER` | `CRANK_NICOLSON` | `BDF2` |

## Picking a tier

| Scenario | Recommended tiers |
|---|---|
| Quick smoke test, comparison-only run | All LOW |
| Close-in V&V (peak velocity, cavity radius) | All MED, optionally HIGH `time_integrator_diffusion` |
| Marshak-tier radiation transport gates | HIGH `radiation_phase`, `cavity_eos`, `opacity_model`; HIGH or HIGHEST `time_integrator_diffusion` |
| Multigroup transport / strict spec on radiation | HIGHEST across all six |

The pass-12 example default for the historic-nuclear suite is **MED**. The
pass-12 example default for the Salmon Marshak showcase is **HIGH**. The
HIGHEST tier is opt-in.

## Per-ladder detail

### `radiation_phase`

Governs how the radiation phase of the explosion source is modeled.

| Tier | Value | Cost | When to pick it |
|---|---|---|---|
| LOW | `ZELDOVICH_RAIZER` | Closed-form end-state. Negligible runtime | Pass-7 closed-form path; fastest available |
| MED | `MARSHAK_GREY` | Adds explicit grey radiation-diffusion solve with Newton T^4 coupling | Marshak-tier radiation gates; matches pass-8 default |
| HIGH | `MARSHAK_GREY` + `TABULATED_PATCHED` opacity | MED + Mihalas-Mihalas opacity blend in [1e5, 1.26e5] K patch | When opacity drift in the patch dominates the residual |
| HIGHEST | `MARSHAK_MULTIGROUP` | G coupled tridiagonal solves per Newton iteration (G=16 default) | Strict-spec radiation gates; multigroup transport required |

`SN_TRANSPORT` is a named-only HIGHEST-tier scaffold; calling it raises a
runtime error. Pass-13 is not adding it.

### `cavity_eos`

Governs the equation of state used inside the inner cavity at the
radiation-to-hydrodynamic transition and during cavity expansion.

| Tier | Value | Cost | When to pick it |
|---|---|---|---|
| LOW | `IDEAL_GAS` | Single closed-form expression | Bring-up only |
| MED | `TILLOTSON` | Two-region analytic Tillotson EOS, four-parameter set per medium (granite, tuff, salt, alluvium) | Pass-7 default; production-baseline runs |
| HIGH | `TILLOTSON_TABULATED_PATCH` | Tillotson blended with Z-R partial-ionization plasma table in [5e10, 6e10] Pa pressure patch | When Tillotson over-predicts cavity pressure in the partial-ionization regime |
| HIGHEST | `TABULATED_FULL` | Pure tabulated EOS everywhere, with Tillotson safety net for out-of-table queries | Strict-spec EOS gates; tabulated reference data required |

Tabulated tier requires `tools/tabulated_data/tables/eos/<medium>_aneos.h5`.
Pass-9 ships ANEOS-derived tables for granite, tuff, salt, alluvium; pass-11
extended Hugoniot coverage is at 30 % envelope (Marsh 1980 granite, McQueen
1970 salt). The 5 % spec target requires a Tillotson refit (axis-4).

### `opacity_model`

Governs the Rosseland and Planck mean opacity used inside the radiation
diffusion solve. Tied to `radiation_phase`: a `radiation_phase = ZELDOVICH_RAIZER`
ignores this knob. Combinations are validated at config load.

| Tier | Value | Cost | When to pick it |
|---|---|---|---|
| LOW | `CONSTANT` | Single per-medium constant | Bring-up only |
| MED | `POWER_LAW_ZR` | Z-R power law in T and rho | Pass-8 default; closed-form scaling |
| HIGH | `TABULATED_PATCHED` | Z-R blended with Mihalas-Mihalas-corrected Rosseland/Planck means in a temperature patch | Marshak-tier and above |
| HIGHEST | `TABULATED_FULL` | Pure tabulated grey opacity, no Z-R fallback | Strict-spec gates |

Tabulated tier requires
`tools/tabulated_data/tables/opacity/<medium>_{rosseland,planck}.h5`. Out-of-
table queries fall back to the analytic with one-time stderr warnings.

### `operator_splitting`

Governs the temporal coupling between the hydrodynamic and Marshak
radiation operators.

| Tier | Value | Convergence | When to pick it |
|---|---|---|---|
| LOW / MED | `LIE` | First order | Pass-8 byte-identical default; bring-up |
| HIGH | `STRANG` | Second order (Strang 1968) | Marshak gates; converges Strang order at substep level |
| HIGHEST | `STRANG_MULTIGROUP` | Second order with per-group Strang | Pair with `radiation_phase = MARSHAK_MULTIGROUP` |

Strang convergence-order can be checked at the substep level via the
`operator_splitting_convergence_diagnostic` flag. Set this to `true`
when validating the splitting; leave at the default `false` for
production runs.

### `time_integrator` (hydro)

Governs the explicit integrator used for the source-region hydrodynamic
operator (1D radial Lagrangian solver).

| Tier | Value | Convergence | When to pick it |
|---|---|---|---|
| LOW / MED / HIGH | `EXPLICIT_EULER` | First order | Pass-9 default; CFL-stable; cheapest |
| HIGHEST | `RK3_SSP` | Third order (Shu-Osher 1988 SSP3) | Strict-spec convergence on shocked profiles |

`TVD_RK2` (Heun's method, second order) is also available as an
intermediate; it is not on the default tier ladder but can be selected
explicitly.

### `time_integrator_diffusion`

Governs the implicit integrator used for the radiation diffusion operator.
Pass-11 introduced this knob to close the Marshak self-similar gate to
factor 2 and the radiation-energy conservation gate to 2 %.

| Tier | Value | Convergence | When to pick it |
|---|---|---|---|
| LOW / MED | `BACKWARD_EULER` | First order, L-stable | Pass-9 default; cheapest |
| HIGH | `CRANK_NICOLSON` | Second order, A-stable | When you need second-order accuracy without BDF2 multi-step state |
| HIGHEST | `BDF2` | Second order, A-stable | Closes Marshak gates at strict spec; pass-11 default for HIGH+ tiers |

## Cavity geometry (axis-1b scaffold)

Pass-11 added `cavity_geometry` as a future-axis scaffold. The current
default and only working value is `RADIAL_LAGRANGIAN` (1D radial
Lagrangian solver). `THREE_DIMENSIONAL` raises a runtime error pointing
at pass-13's axis-1b implementation track. See
[AXIS_1B_DESIGN.md](AXIS_1B_DESIGN.md).

## Backwards compatibility

Every historic-nuclear test config in `examples/` pins the tier values
that match its original pass. Pinned configs continue to reproduce
their pass-N gate envelope byte-for-byte.

| Pass | Default tier (when not overridden) |
|---|---|
| Pre-pass-7 | KINEMATIC_RDP scalar moment-rate path (legacy) |
| Pass-7 | LOW radiation, MED EOS, MED opacity, LOW splitting / hydro, LOW diffusion |
| Pass-8 | MED radiation, MED EOS, MED opacity, LOW splitting, LOW hydro, LOW diffusion |
| Pass-9 | MED radiation, MED EOS, MED opacity, LOW splitting (LIE), LOW hydro, LOW diffusion |
| Pass-10 | unchanged from pass-9; HIGHEST tiers opt-in |
| Pass-11 | unchanged; HIGHEST `time_integrator_diffusion = BDF2` is the new default for HIGH+ tier configs |

When in doubt, leave the tier knobs unset and you reproduce the pinned
behaviour for that example.
