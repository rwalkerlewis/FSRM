# FSRM Configuration

The top-level `config/` directory holds schema-anchor configurations
that are not tied to a specific runnable example. Per-event configs
live next to their runnable examples under `examples/N_<event>/`.

## Schema-anchor templates

| File | Purpose |
|---|---|
| `default.config` | Minimal default config used by the test fixtures |
| `complete_template.config` | Annotated reference covering every recognized section and key |
| `test_elastostatics.config` | Schema fixture for the elastostatic test family |
| `test_elastodynamics.config` | Schema fixture for the elastodynamic test family |
| `test_poroelasticity.config` | Schema fixture for the poroelastic test family |
| `test_locked_fault.config` | Schema fixture for the locked-fault test family |

These are referenced by name in `tests/` and must remain in this
directory.

## Templates (`templates/`)

The `templates/` subdirectory holds per-physics starter configs that
are not driven by integration tests. They are anchored on a single
verified physics path each and may be copied into a fresh example
directory as a starting point. Currently:

- `lambs_problem.config` -- Lamb 1904 surface-wave verification
- `lithostatic_column.config` -- Closed-form lithostatic stress
- `lithostatic_equilibrium.config` -- Initial-state lithostatic balance
- `terzaghi_consolidation.config` -- Terzaghi 1925 1D consolidation
- `layered_elastostatics.config` -- Depth-layered elastic stack
- `injection_pressure_buildup.config` -- Poroelastic injection
- `gmsh_box.config` -- Minimal Gmsh per-region material example
- `minimal_explosion.config` -- Smallest viable Mueller-Murphy explosion
- `explosion_seismogram.config` -- Production-resolution seismogram
- `underground_explosion_template.config` -- Annotated underground template
- `slipping_fault_shear.config` -- Slipping fault under shear
- `prescribed_slip_test.config` -- Prescribed-slip Cartesian-vector exercise
- `fault_compression.config` -- Locked fault under compression
- `cohesive_hydraulic_fracture.config` -- Cohesive hydrofrac stub
- `hydraulic_fracture_pkn.config` -- PKN hydrofrac formulas

## Per-example configs

Pass-12 relocated event-specific configs from `config/examples/<event>.config`
into `examples/N_<event>/config.config`. The `config/examples/` subdirectory
no longer holds production configs.

The previous README index of `config/examples/<event>.config` files is
superseded by the example-directory READMEs.

## See also

- [docs/CONFIGURATION.md](../docs/CONFIGURATION.md) -- Per-section
  per-key reference.
- [docs/USER_GUIDE.md](../docs/USER_GUIDE.md) -- End-to-end worked
  example.
- [docs/FIDELITY_LADDER_GUIDE.md](../docs/FIDELITY_LADDER_GUIDE.md) --
  Fidelity-tier selection.
