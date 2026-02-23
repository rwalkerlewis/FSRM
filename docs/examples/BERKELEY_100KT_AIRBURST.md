# Berkeley 100 kT Airburst Analysis

## Event Summary

| Parameter | Value |
|-----------|-------|
| **Type** | Hypothetical atmospheric nuclear detonation |
| **Location** | 37.8716°N, 122.2727°W (UC Berkeley campus) |
| **Yield** | 100 kt |
| **Height of Burst** | 50 m above ground level |
| **Total Energy** | 4.184 × 10¹⁴ J |
| **Fireball Radius** | ~416 m (maximum) |
| **Cloud Top Height** | ~13.9 km (stabilised) |

## Background

This analysis models the effects of a hypothetical 100 kiloton nuclear airburst over Berkeley, California, at a height of burst (HOB) of 50 metres. The scenario provides a comprehensive assessment of coupled multi-physics effects including blast wave propagation, thermal radiation, prompt nuclear radiation, fallout transport, electromagnetic pulse, ground-coupled seismic effects, and combined hazard mapping.

A 50 m HOB places the detonation low enough to create significant ground coupling and a shallow crater, while still producing substantial Mach stem enhancement of the blast wave in the near field. This produces effects intermediate between a true surface burst and an optimised airburst.

## Physics Models

### Blast Wave (Glasstone-Dolan)

Peak overpressure follows the empirical Glasstone-Dolan scaling with Mach stem enhancement:

$$\Delta P(r) \approx \frac{K}{Z^3} + \frac{L}{Z^2} + \frac{M}{Z}$$

where $Z = R / W^{1/3}$ is the scaled distance and $R = \sqrt{r^2 + h^2}$ is the slant range. Mach stem reflection enhances overpressure by a factor of ~1.8 at angles below 40° from horizontal.

| Damage Level | Overpressure | Radius |
|--------------|-------------|--------|
| Total destruction | >140 kPa | ~1.3 km |
| Reinforced concrete destroyed | >35 kPa | ~2.7 km |
| Most buildings destroyed | >14 kPa | ~4.6 km |
| Moderate damage | >7 kPa | ~6.7 km |
| Light damage | >3.5 kPa | ~9.5 km |
| Glass breakage | >1 kPa | ~18 km |

### Thermal Radiation

Thermal fluence at ground range $r$ with atmospheric transmission:

$$Q(r) = \frac{0.35 \cdot E_{total}}{4\pi R^2} \cdot e^{-R/\lambda}$$

where $\lambda \approx 18$ km is the atmospheric absorption length.

| Effect | Fluence | Radius |
|--------|---------|--------|
| Third-degree burns | >670 kJ/m² | ~3.5 km |
| Second-degree burns | >335 kJ/m² | ~5 km |
| First-degree burns | >125 kJ/m² | ~8 km |
| Pain threshold | >50 kJ/m² | ~13 km |

### Prompt Nuclear Radiation

Prompt gamma and neutron dose with exponential atmospheric attenuation:

$$D(r) = D_{ref} \cdot W \cdot \left(\frac{1000}{R}\right)^2 \cdot e^{-R/\lambda_r}$$

where $\lambda_r \approx 2500$ m is the radiation mean free path in air.

| Effect | Dose | Radius |
|--------|------|--------|
| Fatal | >10 Gy | ~2.3 km |
| LD50/60 | >4.5 Gy | ~2.6 km |
| Radiation sickness | >2 Gy | ~3.0 km |
| Mild symptoms | >0.5 Gy | ~3.7 km |

### Fallout (Gaussian Plume)

Fallout is modelled using a Gaussian plume model with prevailing northwesterly winds (315°, 15 m/s at altitude). The cloud stabilisation height for 100 kt is approximately 13.9 km.

Key parameters:
- **Fission fraction**: 50%
- **Wind direction**: From NW (toward SE) — consistent with Bay Area prevailing winds
- **Downwind deposition**: Oakland, Hayward, Fremont, San Jose corridor

Dose rate decay follows the Way-Wigner rule: $\dot{D}(t) = \dot{D}_1 \cdot t^{-1.2}$

### EMP

For a low-altitude burst (50 m), the EMP is significant but localised:
- E1 component: ~25,000 V/m at burst point, decaying exponentially
- E3 (MHD-EMP): Negligible for surface/low-altitude burst

### Ground-Coupled Seismic Effects

Air-blast couples to the ground through:
1. Direct blast loading
2. Air-coupled Rayleigh wave excitation
3. Minor cratering at ground zero

Estimated air-coupled mb ≈ 1.8 (much weaker than underground coupling).

### Combined Hazard

The dominant lethal mechanism varies with distance:
- **< 1 km**: Blast dominant (total destruction)
- **1-5 km**: Blast and thermal compete
- **5-15 km**: Thermal radiation dominant
- **> 15 km**: Fallout becomes the primary concern (downwind)

## Affected Area

### Immediate Impact Zone (< 5 km)
- UC Berkeley campus
- Downtown Berkeley
- Emeryville
- North Oakland
- Albany
- Parts of El Cerrito and Kensington

### Severe Damage Zone (5-15 km)
- Most of Oakland
- San Francisco (across the Bay — ~13 km)
- Richmond
- Piedmont
- Walnut Creek area (partial)

### Light Damage / Fallout Zone (15-80 km)
- San Leandro, Hayward (blast + fallout)
- Fremont, San Jose (fallout plume)
- Southern Contra Costa County

## Generated Figures

### Matplotlib Analysis (8 figures)

| Figure | Description |
|--------|-------------|
| `fig01_scenario_overview.png` | Scenario parameters, energy partition, overpressure/thermal/radiation vs range |
| `fig02_blast_damage_map.png` | 2D blast damage zones with categorical and continuous overpressure |
| `fig03_thermal_and_radiation.png` | Thermal fluence, prompt radiation dose, and burn zone maps |
| `fig04_fallout_pattern.png` | Fallout plume, dose rate, Way-Wigner decay, cumulative dose |
| `fig05_emp_analysis.png` | E1/E2/E3 waveforms, spatial distribution, power line coupling |
| `fig06_fireball_and_cloud.png` | Fireball expansion, temperature, thermal power, mushroom cloud |
| `fig07_ground_coupling.png` | PGV map, air-ground coupling physics, crater profile |
| `fig08_combined_hazard.png` | Combined hazard index, dominant effect, casualty estimation |

### PyGMT Relief Maps (6 maps)

| Map | Description |
|-----|-------------|
| `map01_regional_overview.png` | Bay Area overview with blast contours on topographic relief |
| `map02_blast_damage.png` | Close-up blast zones on terrain with local landmarks |
| `map03_thermal_radiation.png` | Thermal and radiation contours on terrain |
| `map04_fallout_plume.png` | Regional fallout plume toward SE on terrain |
| `map05_combined_effects.png` | Combined effects with infrastructure (BART, bridges, ports) |
| `map06_radiation_zones.png` | Detailed radiation and thermal zones on terrain |

## Running the Analysis

```bash
# Generate all matplotlib analysis figures
python scripts/model_berkeley_100kt_airburst.py

# Generate all PyGMT maps with terrain relief
python scripts/plot_berkeley_airburst_maps.py
```

## Physical References

- Glasstone & Dolan (1977) — *The Effects of Nuclear Weapons* (3rd Edition)
- Brode (1955) — Numerical Solutions of Spherical Blast Waves
- Sedov (1959) — Similarity and Dimensional Methods in Mechanics
- Way & Wigner (1948) — Rate of Decay of Fission Products
- Bridgman (2001) — *Introduction to the Physics of Nuclear Weapons Effects*
- Mueller & Murphy (1971) — Seismic Characteristics of Underground Nuclear Detonations
