# Berkeley 100 kT Airburst — Hypothetical Atmospheric Nuclear Detonation

## Event Summary

| Parameter | Value |
|-----------|-------|
| **Type** | Hypothetical atmospheric nuclear detonation (airburst) |
| **Location** | 37.8716°N, 122.2727°W |
| **Site** | UC Berkeley campus, Berkeley, California |
| **Yield** | 100 kt |
| **Height of Burst** | 50 m above ground level |
| **Total Energy** | 4.184 × 10¹⁴ J |
| **Fireball Radius** | ~416 m (maximum) |
| **Mushroom Cloud Top** | ~13.9 km (stabilised) |
| **Prevailing Wind** | From NW at 15 m/s (altitude); 5 m/s (surface) |

## Background

This analysis models the coupled multi-physics effects of a hypothetical 100 kiloton nuclear airburst over Berkeley, California, at a height of burst (HOB) of 50 metres. The scenario is chosen to illustrate the full spectrum of nuclear weapons effects — blast, thermal, radiation, fallout, EMP, and ground-coupled seismic response — in an urban environment with complex topography and dense population.

### Why 50 m Height of Burst?

The 50 m HOB produces effects that are intermediate between a true surface burst (HOB = 0) and an optimised airburst. This regime is physically important for several reasons:

1. **Mach Stem Formation**: The incident and reflected shock waves merge near the ground to form a Mach stem, producing peak overpressures up to 1.8× higher than the incident wave alone
2. **Ground Coupling**: The burst is close enough to excavate a shallow crater and couple significant energy into the ground as seismic waves
3. **Fallout**: A 50 m burst entrains substantially less ground material than a contact burst, producing less local fallout — but the fission products in the weapon debris itself still produce significant residual radioactivity
4. **Fireball Contact**: The maximum fireball radius (~416 m) is much larger than the 50 m HOB, so the fireball engulfs the surface, producing cratering and debris lofting

### Geographic Setting

Berkeley sits on the eastern shore of San Francisco Bay at the base of the Berkeley Hills:

- **Bay Area Population**: ~7.8 million in the metropolitan area
- **Terrain**: Flat near the waterfront (0–50 m elevation), rising steeply into the Berkeley/Oakland Hills (300–500 m)
- **Water Bodies**: San Francisco Bay to the west, acting as a thermal and blast buffer for San Francisco (~13 km across the Bay)
- **Infrastructure**: UC Berkeley campus, BART transit system, Bay Bridge, Port of Oakland, Lawrence Berkeley National Laboratory

### Prevailing Winds

The San Francisco Bay Area has a well-characterised wind climatology:

- **Surface winds**: Predominantly from the west/northwest at 3–7 m/s, driven by the Pacific high-pressure system and marine layer dynamics
- **Upper-level winds**: Generally from the west-northwest at 10–25 m/s in the free troposphere
- **Fallout transport**: The NW wind at altitude carries the stabilised mushroom cloud debris toward the southeast — along the Oakland–Hayward–Fremont–San Jose corridor

This analysis uses a prevailing wind direction of 315° (from NW) at 15 m/s at the cloud stabilisation height, and 5 m/s at the surface.

---

## Physics of Atmospheric Nuclear Detonations

An atmospheric nuclear explosion involves coupled phenomena across multiple physical domains. Unlike underground tests (which primarily produce seismic signals), an airburst's energy partitions across blast, thermal, radiation, and electromagnetic channels simultaneously.

### Energy Partition

For a fission/thermonuclear weapon detonated in the atmosphere, the total yield partitions as:

| Energy Channel | Fraction | Physical Mechanism | Primary Effects |
|----------------|----------|--------------------|-----------------|
| Blast / shock | 50% | Hydrodynamic shock wave in air | Structural damage, casualties from debris |
| Thermal radiation | 35% | Blackbody emission from fireball surface | Burns, fires, firestorms |
| Prompt nuclear radiation | 5% | Gamma rays and neutrons emitted in first second | Acute radiation syndrome |
| Residual radiation (fallout) | 10% | Fission product decay | Chronic radiation exposure |

For a 100 kt weapon: $E_{total} = 100 \times 4.184 \times 10^{12} = 4.184 \times 10^{14}$ J.

These fractions are approximate and depend on weapon design (fission vs. thermonuclear), HOB, and atmospheric conditions. A fission-only device produces more fallout; a thermonuclear device with a high fusion fraction produces less.

---

## Blast Wave Physics

### Fireball Formation and Shock Launch

The nuclear detonation releases its energy in less than 1 microsecond. The initial temperature exceeds $10^8$ K, and the energy radiates outward as X-rays that are absorbed within a few metres of air, creating a **fireball** — a luminous sphere of ionised gas.

#### Fireball Maximum Radius (Brode Model)

The fireball expands until its internal pressure equilibrates with the surrounding atmosphere:

$$R_{fb} = 66 \times W^{0.4} \text{ metres}$$

For $W = 100$ kt: $R_{fb} \approx 416$ m.

Since the HOB is only 50 m, the fireball engulfs the ground surface, vaporising and melting surface material and excavating a shallow crater.

#### Fireball Formation Time

$$t_{form} = 0.001 \times W^{0.3} \text{ seconds}$$

For 100 kt: $t_{form} \approx 4$ ms.

### Shock Wave Propagation

#### Sedov-Taylor Solution

At early times, the expanding blast wave is well described by the Sedov-Taylor self-similar solution for a strong point explosion:

$$R(t) = \xi_0 \left(\frac{E}{\rho_0}\right)^{1/5} t^{2/5}$$

Where:
- $R(t)$ = shock front radius at time $t$
- $E$ = total energy deposited
- $\rho_0$ = ambient air density (1.225 kg/m³ at sea level)
- $\xi_0 \approx 1.15$ = dimensionless constant for $\gamma = 1.4$

This gives $R \propto t^{2/5}$ (decelerating expansion). For 100 kt, the shock front reaches 1 km in approximately 0.3 s.

#### Peak Overpressure (Kinney-Graham Formula)

At intermediate to far ranges, the peak static overpressure at ground level is computed using the Kinney-Graham (1985) closed-form empirical fit, which is continuous across the entire range:

$$\frac{\Delta P}{P_0} = \frac{808 \left[1 + \left(\bar{Z}/4.5\right)^2\right]}{\sqrt{1 + \left(\bar{Z}/0.048\right)^2} \cdot \sqrt{1 + \left(\bar{Z}/0.32\right)^2} \cdot \sqrt{1 + \left(\bar{Z}/1.35\right)^2}}$$

Where:
- $\bar{Z} = R / m_{charge}^{1/3}$ is the **Hopkinson-Cranz scaled distance** (m/kg$^{1/3}$)
- $R = \sqrt{r^2 + h^2}$ is the **slant range** from the burst point to ground range $r$ at burst height $h$
- $P_0 = 101.325$ kPa is standard atmospheric pressure
- $m_{charge} = W \times 10^6$ kg TNT equivalent for yield $W$ in kt

This formula is a rational polynomial fit to compiled experimental blast data and reproduces Glasstone-Dolan (1977) reference values within 5–10% across the full range.

*Reference: Kinney, G.F. & Graham, K.J., "Explosive Shocks in Air", 2nd ed., Springer-Verlag, 1985.*

#### Mach Stem Enhancement

When the blast wave reflects from the ground, the incident and reflected waves can merge into a single, stronger **Mach stem**. This occurs when the angle of incidence exceeds the critical angle for regular reflection (~40° from horizontal for strong shocks in air):

$$\alpha = \arctan\left(\frac{h}{r}\right) < 40°$$

In the Mach reflection regime, the effective overpressure is enhanced by a factor of approximately 1.8:

$$\Delta P_{Mach} \approx 1.8 \times \Delta P_{incident}$$

This enhancement is critical for low-altitude bursts like the 50 m HOB scenario. The Mach stem forms at a distance of approximately $r \geq h \cdot \tan(50°) \approx 60$ m from ground zero, and persists to large distances. This means that the blast damage from a 50 m airburst is substantially greater than the free-air blast alone.

#### Dynamic Pressure

The **dynamic pressure** (blast wind) is derived from the overpressure:

$$q = \frac{5}{2} \cdot \frac{(\Delta P)^2}{7P_0 + \Delta P}$$

Where $P_0 = 101.325$ kPa is ambient atmospheric pressure. Dynamic pressure drives wind-borne missile damage and drag forces on structures.

### Damage Criteria

Structural damage is characterised by overpressure thresholds:

| Damage Level | Overpressure | Estimated Radius | Physical Effects |
|--------------|-------------|------------------|-----------------|
| Total destruction | >140 kPa (20 psi) | ~1.4 km | All structures obliterated, crater |
| Reinforced concrete destroyed | >35 kPa (5 psi) | ~2.9 km | Heavy concrete structures collapse |
| Most buildings destroyed | >14 kPa (2 psi) | ~5.7 km | Residential/commercial buildings collapse |
| Moderate damage | >7 kPa (1 psi) | ~10.4 km | Walls cracked, roofs damaged |
| Light damage | >3.5 kPa (0.5 psi) | ~20.2 km | Windows broken, light structural damage |
| Glass breakage | >1 kPa (0.15 psi) | ~70 km | Windows shattered across Bay Area |

The light-damage radius (~20 km) means that window breakage would extend to San Francisco, Oakland hills, Richmond, San Leandro, and most of the inner Bay Area. The 1 kPa glass-breakage threshold potentially extends across the entire Bay Area metropolitan region.

---

## Thermal Radiation Physics

### Emission Mechanism

The fireball surface radiates as an approximate blackbody. Thermal energy is emitted in two pulses:

1. **First pulse** (< 10 ms): Initial X-ray flash, rapidly absorbed by air
2. **Second pulse** (0.1–3 s): Main thermal radiation emission as the shock wave becomes transparent; delivers ~99% of thermal energy

The **total thermal energy** is 35% of yield:
$$E_{thermal} = 0.35 \times 4.184 \times 10^{14} = 1.46 \times 10^{14} \text{ J}$$

### Thermal Fluence at Ground Level

The thermal fluence (energy per unit area) at ground range $r$ is:

$$Q(r) = \frac{E_{thermal}}{4\pi R^2} \cdot \tau(R)$$

Where:
- $R = \sqrt{r^2 + h^2}$ is the slant range from burst to target
- $\tau(R) = e^{-R/\lambda}$ is the atmospheric transmissivity
- $\lambda \approx 18$ km is the atmospheric absorption mean free path (depends on visibility, humidity)

### Atmospheric Transmission

Atmospheric absorption reduces thermal fluence through:

1. **Molecular absorption**: H₂O and CO₂ absorption bands
2. **Aerosol scattering**: Mie scattering by particulates
3. **Rayleigh scattering**: Molecular scattering (minor at visible/IR wavelengths)

For the Bay Area (frequently foggy or hazy), the effective visibility can be significantly reduced, which would decrease thermal damage radii. The analysis uses a clear-day transmissivity ($\lambda = 18$ km), representing a worst case for thermal effects.

### Burn and Fire Thresholds

| Effect | Fluence | Radius | Physical Mechanism |
|--------|---------|--------|-------------------|
| Third-degree burns | >670 kJ/m² | ~3.5 km | Full-thickness skin destruction |
| Second-degree burns | >335 kJ/m² | ~5 km | Blistering of exposed skin |
| First-degree burns | >125 kJ/m² | ~8 km | Reddening, similar to severe sunburn |
| Pain threshold | ~50 kJ/m² | ~13 km | Brief pain on exposed skin |

#### Ignition and Firestorm Potential

Thermal radiation ignites combustible materials:
- **Dark fabrics**: Ignite at ~250 kJ/m² (~4 km)
- **Newspaper/leaves**: Ignite at ~170 kJ/m² (~5 km)
- **Dry wood**: Ignite at ~500 kJ/m² (~3.5 km)

In the 2–5 km zone, numerous fires would start simultaneously. If these merge before firefighting response, a **mass fire or firestorm** can develop — a self-sustaining conflagration with hurricane-force inward winds that greatly increases casualties and destruction beyond what blast alone would produce.

---

## Prompt Nuclear Radiation

### Source Mechanisms

Prompt (initial) radiation is emitted within the first second of detonation:

1. **Gamma rays**: High-energy photons from fission reactions and neutron capture
2. **Neutrons**: Fast and thermal neutrons from fission/fusion reactions

For a 100 kt weapon, approximately 5% of the yield ($2.09 \times 10^{13}$ J) is released as prompt radiation.

### Dose Model

The prompt radiation dose at ground range $r$ is modelled as:

$$D(r) = D_{ref} \cdot W \cdot \left(\frac{1000}{R}\right)^2 \cdot e^{-R / \lambda_r}$$

Where:
- $D_{ref} = 10$ Gy at 1 km for 1 kt (reference dose)
- $W$ = yield in kt
- $R = \sqrt{r^2 + h^2}$ = slant range (m)
- $\lambda_r \approx 2500$ m = radiation mean free path in air

The $1/R^2$ term accounts for geometric spreading, while the exponential term represents atmospheric attenuation (pair production, Compton scattering, and photoelectric absorption for gamma rays; elastic scattering and capture for neutrons).

### Biological Effects

| Effect | Dose | Radius | Onset |
|--------|------|--------|-------|
| Fatal (>10 Gy) | >10 Gy | ~2.3 km | Death within hours to days |
| LD50/60 (~4.5 Gy) | ~4.5 Gy | ~2.6 km | 50% mortality within 60 days (untreated) |
| Acute radiation sickness | >2 Gy | ~3.0 km | Nausea, vomiting, immunosuppression |
| Mild symptoms | >0.5 Gy | ~3.7 km | Transient nausea, reduced blood counts |
| Detectable | >10 mGy | ~5.5 km | Measurable by dosimetry, no clinical symptoms |

**Note**: Within the blast destruction zone (~2.7 km for severe), casualties from blast and thermal effects would dominate over radiation. Prompt radiation becomes the *primary* lethal mechanism only at ranges where blast is survivable but radiation dose is high — a narrow annular zone that shrinks with increasing yield (because blast scales as $W^{1/3}$ while radiation attenuation is exponential and yield-independent).

---

## Fallout Physics

### Formation of Radioactive Debris

A 100 kt detonation produces approximately $3 \times 10^{26}$ fission reactions (assuming 50% fission fraction), creating ~100 kg of highly radioactive fission products comprising over 300 different isotopes.

At a 50 m HOB, the fireball engulfs the ground surface, lofting vaporised soil and debris into the rising mushroom cloud. This **activated debris** mixes with fission products, condenses as the cloud rises and cools, and eventually falls out downwind.

### Mushroom Cloud Stabilisation

The cloud rises buoyantly until it reaches an altitude where its density matches the ambient atmosphere:

$$H_{top} \approx 2200 \times W^{0.4} \text{ metres}$$

For 100 kt: $H_{top} \approx 13{,}900$ m (just above the tropopause at ~11 km). The cloud stabilises near the tropopause and spreads laterally.

### Gaussian Plume Transport Model

Fallout transport is modelled using a Gaussian plume in wind-rotated coordinates. The wind vector at cloud altitude is decomposed into downwind and crosswind components:

$$A(x, y) = \exp\left[-\frac{1}{2}\left(\frac{d - d_0}{\sigma_d}\right)^2\right] \cdot \exp\left[-\frac{1}{2}\left(\frac{c}{\sigma_c}\right)^2\right]$$

Where:
- $d$ = downwind distance (along wind vector, rotated by $\theta = 315°$)
- $c$ = crosswind distance (perpendicular to wind)
- $d_0 = v_{wind} \times H_{top} / v_{fall}$ = mean drift distance (wind speed × fall time)
- $\sigma_d$ = downwind dispersion (proportional to distance + constant)
- $\sigma_c$ = crosswind dispersion (proportional to distance + constant)

The prevailing NW wind ($\theta = 315°$) carries fallout toward the **southeast**: through Oakland, San Leandro, Hayward, Fremont, and potentially reaching San Jose for the lightest particles.

### Dose Rate Decay — Way-Wigner Rule

The fallout dose rate decays following the Way-Wigner law (Way & Wigner, 1948) for fission product decay:

$$\dot{D}(t) = \dot{D}_1 \cdot t^{-1.2}$$

Where $\dot{D}_1$ is the dose rate at $t = 1$ hour after detonation. This power law arises from the superposition of thousands of fission product isotopes with different half-lives, producing an aggregate decay that is faster than any single exponential.

Key consequence: **The 7-10 Rule** — for every 7-fold increase in time, the dose rate decreases by a factor of 10. After 49 hours (~2 days), fallout radiation is reduced to 1% of its 1-hour value.

### Cumulative Dose Integration

The cumulative dose from fallout exposure from time $t_1$ to $t_2$ is:

$$D(t_1, t_2) = \frac{\dot{D}_1}{-0.2} \left[ t_2^{-0.2} - t_1^{-0.2} \right]$$

This integral converges, meaning the total dose accumulated from $t = 1$ h to infinity is finite (approximately $5 \times \dot{D}_1$ Sv·h).

### Airburst vs. Surface Burst Fallout

A critical distinction: a 50 m airburst produces **significantly less fallout** than a contact surface burst of the same yield:

| Feature | 50 m Airburst | Contact Surface Burst |
|---------|--------------|----------------------|
| Ground material lofted | Moderate (crater ejecta) | Very large (deep crater) |
| Local fallout | Moderate | Severe |
| Downwind hazard distance | ~50–80 km | ~200+ km |
| Time to first fallout | 30–60 min | 15–30 min |
| Dominant particle size | Fine (< 100 µm) | Coarse + fine (mm to µm) |

---

## Electromagnetic Pulse (EMP)

### Physical Origin

A nuclear detonation produces electromagnetic pulses through several mechanisms, categorised as E1, E2, and E3 components.

### E1 Component (Fast)

- **Mechanism**: Compton scattering of gamma rays creates a coherent current of relativistic electrons in the atmosphere
- **Rise time**: ~2.5 nanoseconds
- **Duration**: ~1 microsecond
- **Peak field**: ~25,000 V/m at burst point for 100 kt (decays exponentially with distance)
- **Frequency**: 1 MHz – 100 MHz

$$E_1(r) = E_0 \cdot \left(\frac{W}{100}\right)^{0.5} \cdot e^{-r / \lambda_{EMP}}$$

Where $\lambda_{EMP} \approx 30$ km for a surface/low-altitude burst.

### E2 Component (Intermediate)

- **Mechanism**: Scattered gamma rays and inelastic neutron reactions
- **Duration**: 1 µs – 1 s
- **Character**: Similar to lightning-induced surges
- **Field strength**: 100–1000 V/m

### E3 Component (MHD-EMP)

- **Mechanism**: Distortion of Earth's magnetic field by expanding plasma
- **Duration**: 1–1000 seconds
- **Relevance**: Primarily significant for high-altitude bursts (>30 km HOB); **negligible for the 50 m burst scenario**

### Infrastructure Vulnerability

For the 50 m airburst, EMP effects are localised but intense:

- **0–10 km**: Severe damage to unshielded electronics, SCADA systems, telecommunications
- **10–30 km**: Moderate disruption to electronics
- **>30 km**: Minimal EMP effects

The Bay Area's dense electronic infrastructure (data centres, BART control systems, telecommunications) would be vulnerable within ~20 km.

---

## Fireball and Mushroom Cloud Evolution

### Phase 1: Initial Flash (< 1 ms)

X-ray emission from the detonation point is absorbed by surrounding air within a few metres, creating the initial fireball. Peak surface temperature exceeds $10^7$ K.

### Phase 2: First Minimum (~10 ms)

The expanding shock wave outruns the luminous front, and the fireball surface temperature drops as the shock becomes opaque. Brightness decreases sharply — the **first minimum** of the characteristic optical **double flash**.

### Phase 3: Second Maximum (~0.3 s)

As the shock wave expands and weakens, it becomes transparent to thermal radiation from the interior. The fireball surface re-brightens to a peak surface temperature of ~8000 K (hotter than the Sun's photosphere at 5778 K). This **second maximum** delivers ~99% of the total thermal energy.

### Phase 4: Buoyant Rise (1–120 s)

The hot fireball gases rise buoyantly, entraining surrounding air and forming the characteristic mushroom cloud. The stem is fed by inrushing surface air carrying dust and debris.

### Phase 5: Cloud Stabilisation (2–10 minutes)

The cloud rises until it reaches neutral buoyancy near the tropopause:

$$H_{top} = 2200 \times W^{0.4} \approx 13{,}900 \text{ m}$$

The cloud then spreads laterally and begins to be transported by upper-level winds.

### Thermal Power

The instantaneous thermal power emitted by the fireball surface:

$$P(t) = 4\pi R_{fb}(t)^2 \cdot \sigma \cdot T(t)^4$$

Where $\sigma = 5.67 \times 10^{-8}$ W/m²/K⁴ is the Stefan-Boltzmann constant. At the second maximum ($R_{fb} \approx 416$ m, $T \approx 8000$ K), the thermal power exceeds $10^{18}$ W (1 exawatt) — briefly outshining the Sun as seen from the Bay Area.

---

## Ground-Coupled Seismic and Acoustic Effects

### Air-Ground Coupling Mechanisms

An atmospheric blast wave couples energy into the ground through three mechanisms:

#### 1. Direct Blast Loading

The overpressure from the blast wave acts as a transient load on the ground surface, generating compressional (P) and shear (S) waves in the subsurface:

$$\sigma_{zz}(r, t) = \Delta P(r, t)$$

The induced ground motion scales with the overpressure and the acoustic impedance of the near-surface material:

$$v_{ground} = \frac{\Delta P}{\rho_{rock} \cdot V_s} \cdot \eta$$

Where $\eta \approx 0.5$ is an empirical coupling coefficient.

#### 2. Air-Coupled Rayleigh Wave

Because the blast wave propagates faster than seismic surface waves near the source (but decelerates), there exists a range beyond which the blast-wave velocity matches the surface-wave phase velocity (~300 m/s for Rayleigh waves). In this regime, the blast wave continuously excites the ground surface, producing an efficient **air-coupled Rayleigh wave** that can propagate to great distances.

#### 3. Crater Coupling

At 50 m HOB with a 100 kt yield, the fireball contacts the ground and excavates a shallow crater (~130 m diameter, ~15 m deep). This direct ground contact provides additional high-frequency seismic energy.

### Seismic Magnitude

The air-coupled seismic signal is much weaker than a buried explosion of the same yield:

$$m_b^{air} \approx 4.0 + 0.75 \log_{10}(W \cdot \eta_{coupling})$$

With $\eta_{coupling} \approx 10^{-3}$ (air-to-ground coupling efficiency):

$$m_b^{air} \approx 4.0 + 0.75 \log_{10}(0.1) = 4.0 - 0.75 = 1.8$$

Compare with an underground test of the same yield: $m_b^{UG} \approx 4.0 + 0.75 \times 2 = 5.5$.

### Peak Ground Velocity (PGV)

The estimated PGV at distance $r$:

$$PGV(r) = \frac{\Delta P(r)}{\rho \cdot V_s} \cdot \eta$$

For $\rho = 2700$ kg/m³, $V_s = 3000$ m/s:
- At 1 km: PGV ~ 5 cm/s (perceptible, minor damage)
- At 5 km: PGV ~ 0.2 cm/s (felt but not damaging)
- At 20 km: PGV ~ 0.01 cm/s (detectable by seismometers)

### Crater Formation

At 50 m HOB, a shallow crater forms where the fireball contacts the ground. Empirical scaling (Glasstone-Dolan):

- **Apparent crater diameter**: ~130 m
- **Apparent crater depth**: ~15 m
- **Lip height**: ~5 m
- **Ejecta volume**: ~100,000 m³

This crater is much smaller than what a contact surface burst would produce, because most of the energy is deposited in the atmosphere rather than the ground.

---

## Combined Hazard Analysis

### Dominant Lethal Mechanism by Range

The dominant cause of casualties varies with distance from ground zero:

| Range | Dominant Effect | Lethality | Notes |
|-------|----------------|-----------|-------|
| 0–1.3 km | Blast | 100% | Total destruction zone |
| 1.3–3 km | Blast + Thermal | >90% | Severe blast, 3rd-degree burns, lethal radiation |
| 3–5 km | Thermal + Blast | 50–90% | Buildings collapse, mass fires |
| 5–10 km | Thermal | 15–50% | Burns on exposed skin, fires start |
| 10–18 km | Blast (light) | 2–15% | Window breakage, flying glass casualties |
| 18–80 km (downwind) | Fallout | Variable | Depends on shelter, evacuation timing |

### Combined Hazard Index

The combined hazard is computed as the maximum of normalised individual hazards:

$$H_{combined} = \max\left(\frac{\Delta P}{35}, \frac{Q}{670}, \frac{D}{6}\right)$$

Where the denominators are approximate lethality thresholds (35 kPa for blast, 670 kJ/m² for thermal, 6 Gy for radiation). This index varies from 0 (safe) to ≥1 (lethal).

---

## Affected Area — Bay Area Geography

### Immediate Destruction Zone (< 3 km)

- UC Berkeley campus (ground zero)
- Downtown Berkeley
- North Oakland (Temescal, Rockridge)
- Berkeley Marina
- Albany
- Parts of Emeryville
- Kensington
- Lawrence Berkeley National Laboratory

All structures levelled. Near-total casualties among exposed population.

### Severe Damage Zone (3–7 km)

- Most of Oakland (downtown, Lake Merritt, Piedmont)
- Emeryville (full)
- El Cerrito
- Parts of Richmond
- Berkeley Hills (partial shielding by terrain)

Most buildings destroyed or severely damaged. Mass fires likely. High casualty rates.

### Moderate Damage Zone (7–18 km)

- San Francisco Financial District (~13 km across Bay)
- Bay Bridge (severe blast damage on eastern span)
- Richmond
- Walnut Creek (partial)
- San Leandro
- Alameda
- Golden Gate Bridge (light damage, ~20 km)

Moderate to light structural damage. Broken windows throughout. Thermal burns on exposed individuals facing the fireball.

### Fallout Zone (20–80 km downwind)

Under prevailing NW winds, the fallout plume extends to the southeast:

- San Leandro → Hayward → Fremont → Milpitas → San Jose
- Crosswind spread: Union City, Newark, Santa Clara

Sheltering for 24–48 hours would dramatically reduce fallout exposure due to Way-Wigner decay.

---

## Generated Figures

All figures are located in `figures/berkeley_100kt/`.

### Matplotlib Analysis Figures (8 figures)

#### Figure 1 — Scenario Overview (`fig01_scenario_overview.png`)

Six-panel summary:
- **(a)** Energy partition pie chart
- **(b)** Fireball and shock cross-section at $t = 0.1$ s
- **(c)** Overpressure vs. ground range (log-linear) with damage thresholds
- **(d)** Thermal fluence vs. range with burn thresholds
- **(e)** Prompt radiation dose vs. range with biological thresholds
- **(f)** Parameter table with all scenario constants and damage radii

#### Figure 2 — Blast Damage Map (`fig02_blast_damage_map.png`)

Two-panel map:
- **(a)** Categorical damage zones (6 levels from glass breakage to total destruction) with Bay Area landmarks
- **(b)** Continuous overpressure field (log-coloured contour fill) with damage contour overlay

#### Figure 3 — Thermal & Radiation (`fig03_thermal_and_radiation.png`)

Three-panel figure:
- **(a)** Thermal fluence map (log-scale, inferno colourmap)
- **(b)** Prompt radiation dose map (log-scale, YlOrRd colourmap)
- **(c)** Categorical thermal burn zone map (pain / 1st° / 2nd° / 3rd°)

#### Figure 4 — Fallout Pattern (`fig04_fallout_pattern.png`)

Six-panel figure:
- **(a)** Fallout deposition map (normalised, wind arrow showing NW direction)
- **(b)** Dose rate at H+1 hour (Sv/h, log-scale) with evacuation contours
- **(c)** Way-Wigner dose rate decay curves for hotspot, 10%, 1% levels
- **(d)** Cumulative fallout dose vs. time with shelter/evacuation thresholds
- **(e)** Downwind fallout cross-section profiles at different crosswind offsets
- **(f)** Fallout parameter summary and protective action guidance

#### Figure 5 — EMP Analysis (`fig05_emp_analysis.png`)

Six-panel figure:
- **(a)** E1 waveform at 5, 10, 20, 50 km
- **(b)** E1 peak field vs. distance
- **(c)** E1 spatial distribution map
- **(d)** Induced voltage on power lines of various lengths
- **(e)** EMP frequency spectrum (E1/E2/E3 components)
- **(f)** Infrastructure vulnerability assessment table

#### Figure 6 — Fireball & Mushroom Cloud (`fig06_fireball_and_cloud.png`)

Five-panel figure:
- **(a)** Fireball and shock radius vs. time (log-log, Sedov-Taylor solution)
- **(b)** Fireball surface temperature vs. time (showing $10^8$ K → 8000 K → cooling)
- **(c)** Thermal radiation power vs. time (PW scale)
- **(d)** Mushroom cloud rise sequence (1 s to 120 s, with tropopause reference)
- **(e)** Optical double flash time history

#### Figure 7 — Ground Coupling (`fig07_ground_coupling.png`)

Six-panel figure:
- **(a)** Peak ground velocity (PGV) vs. range with damage thresholds
- **(b)** Blast pressure time histories at 1, 3, 5, 10, 20 km (Friedlander waveform)
- **(c)** Acoustic/shock arrival times vs. distance
- **(d)** PGV 2D map (viridis colourmap, log-scale)
- **(e)** Air-ground coupling physics summary
- **(f)** Crater profile estimate (depth, lip height, ejecta)

#### Figure 8 — Combined Hazard (`fig08_combined_hazard.png`)

Six-panel figure:
- **(a)** Combined hazard index map (0=safe, 1=lethal, RdYlGn_r colourmap)
- **(b)** Dominant lethal effect map (blast/thermal/radiation categorical)
- **(c)** All effects vs. range on single axis (overpressure, thermal, dose, dynamic P)
- **(d)** Casualty estimation by annular ring (fatalities + injuries bar chart)
- **(e)** Cumulative casualties vs. radius
- **(f)** Casualty and damage summary table

### PyGMT Relief Maps (8 maps)

All maps use SRTM 3-arc-second (`@earth_relief_03s`) topographic relief as background, providing terrain context for the Bay Area's complex geography. Scale bars are placed below each map frame to avoid obscuring terrain detail.

#### Map 1 — Regional Overview (`map01_regional_overview.png`)

Bay Area regional map showing:
- Topographic relief from coast to East Bay hills
- Blast damage contours (1–140 kPa) overlaid on terrain
- Major cities annotated
- Inset map showing California context
- Scale bar below the map frame
- Legend with overpressure zones and damage radii in km

#### Map 2 — Blast Damage Zones (`map02_blast_damage.png`)

Close-up Berkeley/Oakland map:
- High-resolution terrain with shaded relief
- Blast damage circles at 6 overpressure levels (1–140 kPa)
- Local neighbourhood labels (UC Berkeley, Downtown, Rockridge, etc.)
- Legend with peak overpressure zones and radii in km
- Scale bar below the map frame

#### Map 3 — Thermal & Radiation Zones (`map03_thermal_radiation.png`)

Close-up map with:
- Thermal burn contours (dashed lines: 50, 125, 335, 670 kJ/m²)
- Radiation dose contours (solid lines: 0.5, 2.0, 4.5, 10 Gy)
- Terrain shading showing Berkeley Hills shadow effects
- Legend distinguishing thermal fluence (kJ/m²) from absorbed dose (Gy)

#### Map 4 — Fallout Plume (`map04_fallout_plume.png`)

Regional map (Bay Area to South Bay):
- Fallout contour ellipses oriented along NW→SE wind direction
- 6 activity levels from trace to very high
- Affected cities annotated along the plume corridor
- Wind direction arrow
- Scale bar below the map frame

#### Map 5 — Combined Effects with Infrastructure (`map05_combined_effects.png`)

Intermediate-scale map with:
- Blast damage contours (5 overpressure levels)
- Infrastructure markers: BART stations, Bay Bridge, Golden Gate Bridge, Richmond Bridge, Port of Oakland, SFO/OAK airports
- University/laboratory locations
- Legend distinguishing infrastructure types and overpressure zones with radii

#### Map 6 — Radiation Zones Close-up (`map06_radiation_zones.png`)

High-resolution close-up showing:
- Prompt radiation dose contours (10 mGy to 10 Gy)
- Thermal burn contours (1st° through 3rd°)
- Scale bar below the map frame
- Terrain detail showing topographic shielding

#### Map 7 — Terrain-Aware Blast Damage (`map07_blast_terrain.png`)

Close-up Berkeley/Oakland map with DEM-based terrain interaction:
- Blast overpressure contours modified by terrain slope and shielding
- Forward-facing slopes show enhanced overpressure (up to 2×) from blast reflection
- Terrain shadows behind ridges show reduced overpressure (0.3–0.7×)
- Berkeley/Oakland Hills, Claremont Canyon labelled as terrain features
- Overpressure computed on the full 3-arc-second DEM grid using vectorized numpy operations
- Legend with terrain-modified overpressure zones

#### Map 8 — Terrain-Aware Combined Effects (`map08_combined_terrain.png`)

Intermediate-scale map with full terrain interaction:
- Combined hazard index contours (10%–90%) accounting for terrain
- Blast overpressure modified by terrain slope/shielding factors
- Thermal radiation and prompt radiation blocked by line-of-sight terrain occlusion
- LOS visibility computed via vectorized ray-tracing through the DEM (25 samples per ray)
- Purple dashed contours mark areas where hills block thermal and nuclear radiation
- Legend with combined hazard levels and terrain effect descriptions

---

## Running the Analysis

### Prerequisites

System packages (included in Docker dev environment):
```bash
sudo apt-get install gmt gmt-dcw gmt-gshhg ghostscript
```

Python packages (included in Docker dev environment):
```bash
pip install numpy scipy matplotlib pandas obspy pygmt
```

A `libgmt.so` symlink may be required for PyGMT on Ubuntu:
```bash
sudo ln -sf /usr/lib/x86_64-linux-gnu/libgmt.so.6 /usr/lib/x86_64-linux-gnu/libgmt.so
sudo ldconfig
```

### Execution

```bash
# Generate all 8 matplotlib analysis figures
python3 scripts/model_berkeley_100kt_airburst.py

# Generate all 8 PyGMT maps with terrain relief
python3 scripts/plot_berkeley_airburst_maps.py
```

All output is saved to `figures/berkeley_100kt/`.

Maps 1–6 use simple circular contours from the physics functions. Maps 7–8 load the SRTM DEM grid and compute terrain-aware effects (vectorized numpy — no per-pixel Python loops), so they take slightly longer but still complete in under a minute.

---

## Algorithms and Implementation

### Overpressure Calculation (`peak_overpressure_kpa`)

1. Compute slant range: $R = \sqrt{r^2 + h^2}$
2. Convert yield to kg TNT equivalent: $m = W \times 10^6$ kg
3. Compute Hopkinson-Cranz scaled distance: $\bar{Z} = R / m^{1/3}$ (m/kg$^{1/3}$)
4. Apply Kinney-Graham (1985) continuous closed-form fit:
$$\frac{\Delta P}{P_0} = \frac{808 \left[1 + (\bar{Z}/4.5)^2\right]}{\sqrt{1 + (\bar{Z}/0.048)^2} \cdot \sqrt{1 + (\bar{Z}/0.32)^2} \cdot \sqrt{1 + (\bar{Z}/1.35)^2}}$$
5. Determine angle from horizontal: $\alpha = \arctan(h/r)$
6. Apply Mach stem enhancement factor (1.8×) if $\alpha < 40°$

This formula directly produces overpressure in kPa and is calibrated against compiled experimental blast data. It supersedes the earlier piecewise formulation which had discontinuities at segment boundaries.

### Thermal Fluence Calculation (`thermal_fluence_kj`)

1. Compute slant range: $R = \sqrt{r^2 + h^2}$
2. Compute inverse-square fluence: $Q_0 = 0.35 \cdot E / (4\pi R^2)$
3. Apply atmospheric transmission: $Q = Q_0 \cdot e^{-R/18000}$
4. Return in kJ/m²

### Fallout Transport (`fallout_activity`)

1. Rotate coordinate system to align with wind direction ($\theta = 315°$)
2. Decompose position into downwind ($d$) and crosswind ($c$) components
3. Compute mean drift distance: $d_0 = v_{wind} \times H/v_{fall}$
4. Compute along-wind dispersion: $\sigma_d = 0.08 \cdot d + 800$
5. Compute cross-wind dispersion: $\sigma_c = 0.05 \cdot d + 500$
6. Apply bivariate Gaussian: $A = \exp[-\frac{1}{2}((d-d_0)/\sigma_d)^2] \cdot \exp[-\frac{1}{2}(c/\sigma_c)^2]$

### PyGMT Map Generation

Maps 1–6 follow a consistent pipeline:
1. Load `@earth_relief_03s` (SRTM 3-arc-second) as background grid
2. Apply shaded relief (`shading=True`, `cmap="geo"`)
3. Overlay coastlines and borders from GMT databases
4. Compute damage/effect radii using physics functions
5. Convert radii (metres) to lon/lat circles using local scale factors
6. Plot circles as polylines with styled pens
7. Add infrastructure markers, labels, legends
8. Place scale bar below map frame (not on the terrain) for readability
9. Export at 300 DPI

Maps 7–8 extend this with terrain-aware physics:
1. Load `@earth_relief_03s` DEM as an xarray grid
2. Compute vectorized line-of-sight (LOS) visibility from burst point to every grid cell using bilinear-interpolated ray sampling (25 samples per ray) — blocks thermal and nuclear radiation where terrain intercepts the line of sight
3. Compute vectorized blast terrain modification factors using terrain gradient (slope facing the burst direction): forward-facing slopes get enhanced overpressure (up to 2×), areas behind ridges get reduced overpressure (0.3–0.7×)
4. Multiply base physics values by terrain factors to produce terrain-modified grids
5. Write grids to temporary NetCDF files and contour with `grdcontour`
6. Export at 300 DPI

---

## Comparison with Underground Test Analyses

| Feature | DPRK 2017 / Lop Nor 2020 / Divider | Berkeley 100 kT Airburst |
|---------|-------------------------------------|--------------------------|
| **Source type** | Underground (coupled/decoupled) | Atmospheric (50 m HOB) |
| **Primary output** | Seismic waves | Blast + thermal + radiation |
| **Source model** | Mueller-Murphy RDP | Glasstone-Dolan / Sedov-Taylor |
| **Key observations** | Seismograms at regional stations | Damage contours on terrain |
| **Discrimination** | Explosion vs. earthquake | N/A (hypothetical) |
| **Fallout** | None (contained underground) | Gaussian plume, Way-Wigner decay |
| **Maps** | Station locations on relief | Damage zones on relief |
| **Ground motion** | Strong seismic coupling | Weak air-ground coupling |

---

## C++ Implementation — `NuclearAirburstEffects`

All physics calculations from the Python analysis scripts are also available as a generic, config-driven C++ class in the FSRM library. This allows the entire exercise to be reproduced for different parameters without modifying or recompiling code.

### Files

| File | Purpose |
|------|---------|
| `include/NuclearAirburstEffects.hpp` | Header — structs, class, constants |
| `src/NuclearAirburstEffects.cpp` | Full implementation of all physics models |
| `examples/nuclear_airburst_effects.cpp` | Standalone CLI driver executable |
| `config/nuclear_airburst_berkeley_100kt.config` | INI config reproducing the Berkeley 100 kt scenario |

### Quick Start

```bash
# Build
cd build && cmake .. -DCMAKE_BUILD_TYPE=Release && make -j$(nproc) nuclear_airburst_effects

# Run with the Berkeley 100 kt config
./examples/nuclear_airburst_effects -c examples/config/nuclear_airburst_berkeley_100kt.config

# Generate a template config for a new scenario
./examples/nuclear_airburst_effects --generate-template my_scenario.config
```

### Parameterisation

Every physical constant is exposed as a config parameter:

| Section | Key Parameters |
|---------|---------------|
| `[AIRBURST]` | `yield_kt`, `burst_height_m`, `fission_fraction`, `latitude`, `longitude`, `location_name`, `thermal_attenuation_length_m` |
| `[WIND]` | `speed_surface_mps`, `speed_altitude_mps`, `direction_deg` |
| `[GROUND]` | `density_kgm3`, `vs_mps`, `seismic_coupling` |
| `[POPULATION]` | `density_urban`, `density_suburban` |
| `[RADIATION]` | `prompt_dose_ref_gy`, `prompt_atten_length_m` |
| `[EMP]` | `e0_vm`, `decay_length_m`, `tau_rise_s`, `tau_decay_s` |
| `[OUTPUT]` | `radial_profile`, `grid_csv`, `print_summary`, file names, grid extent/resolution |

### Physics Methods

The `NuclearAirburstEffects` class provides these methods:

- `peakOverpressure(ground_range_m)` — Kinney-Graham (1985) with Mach stem
- `dynamicPressure(overpressure_kpa)` — Rankine-Hugoniot
- `blastRadiusForPressure(target_kpa)` — Bisection search
- `thermalFluence(ground_range_m)` — 35% yield partition, inverse-square with atmospheric absorption
- `fireballMaxRadius()` — $R = 66 W^{0.4}$
- `promptRadiationDose(ground_range_m)` — Parametric dose model
- `cloudTopHeight()` — $H = \min(25000, 2200 W^{0.4})$
- `falloutActivity(x_m, y_m)` — Gaussian plume with wind rotation
- `falloutDoseRate1hr(x_m, y_m)` — Way-Wigner decay
- `empE1Peak(ground_range_m)` — E1 peak field strength
- `empE1Waveform(t_s, ground_range_m)` — E1 time-domain signal
- `fireballRadius(t_s)` — Time-dependent fireball growth
- `shockRadius(t_s)` — Sedov-Taylor self-similar solution
- `peakGroundVelocity(ground_range_m)` — Acoustic impedance coupling
- `seismicMagnitude()` — $m_b = 4.0 + 0.75 \log_{10}(W \times 10^{-3})$

### Output Formats

- **Console summary** — Damage radii, thermal/radiation radii, PGV, casualty estimates
- **Radial profile CSV** — All quantities vs. ground range (configurable resolution)
- **2D grid CSV** — All quantities on a regular latitude/longitude grid with fallout plume

### Example: Modelling a Different Scenario

To model a 20 kt airburst at 500 m over a different city:

```ini
[AIRBURST]
yield_kt       = 20.0
burst_height_m = 500.0
latitude       = 40.7128
longitude      = -74.0060
location_name  = Midtown Manhattan

[WIND]
speed_altitude_mps = 20.0
direction_deg      = 270.0

[POPULATION]
density_urban = 10000
```

---

## References

1. **Glasstone, S. & Dolan, P.J. (1977)** — *The Effects of Nuclear Weapons*, 3rd Edition. US DoD / US DoE.

2. **Brode, H.L. (1955)** — Numerical solutions of spherical blast waves. *Journal of Applied Physics* 26(6), 766-775.

3. **Sedov, L.I. (1959)** — *Similarity and Dimensional Methods in Mechanics*. Academic Press.

4. **Way, K. & Wigner, E.P. (1948)** — The rate of decay of fission products. *Physical Review* 73(11), 1318-1330.

5. **Bridgman, C.J. (2001)** — *Introduction to the Physics of Nuclear Weapons Effects*. Defense Threat Reduction Agency.

6. **Brode, H.L. (1968)** — Review of nuclear weapons effects. *Annual Review of Nuclear Science* 18, 153-202.

7. **Longmire, C.L. (1978)** — On the electromagnetic pulse produced by nuclear explosions. *IEEE Transactions on Antennas and Propagation* 26(1), 3-13.

8. **Glasstone, S. (1962)** — *The Effects of Nuclear Weapons* (revised edition). US Atomic Energy Commission.

9. **Dolan, P.J. (1972)** — *Capabilities of Nuclear Weapons*. Defense Nuclear Agency, Effects Manual EM-1.

10. **Mueller, R.A. & Murphy, J.R. (1971)** — Seismic characteristics of underground nuclear detonations. Part I: Seismic source. *BSSA* 61(6), 1675-1692.

11. **Patton, H.J. (1988)** — Source models of the Harzer explosion from regional observations of fundamental-mode and higher-mode surface waves. *BSSA* 78(5), 1133-1157.

12. **FEMA (2010)** — *Planning Guidance for Response to a Nuclear Detonation*, 2nd Edition.

13. **Kinney, G.F. & Graham, K.J. (1985)** — *Explosive Shocks in Air*, 2nd Edition. Springer-Verlag.

13. **Buddemeier, B.R. & Dillon, M.B. (2009)** — *Key Response Planning Factors for the Aftermath of Nuclear Terrorism*. Lawrence Livermore National Laboratory, LLNL-TR-410067.
