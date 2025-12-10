#!/usr/bin/env python3
"""
Nuclear Reservoir Stimulation Visualization
============================================

Visualizes and compares nuclear stimulation scenarios:
1. US Project Gasbuggy (1967) - 29 kt at 1292m depth
2. Soviet Urta-Bulak style - 30 kt at 1500m depth

Based on config files:
- config/nuclear_stimulation_gasbuggy.config
- config/nuclear_stimulation_soviet.config
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch, Polygon
from matplotlib.collections import PatchCollection
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os

# Set style
plt.style.use('seaborn-v0_8-whitegrid')

# Create output directory
output_dir = 'output/nuclear_stimulation'
os.makedirs(output_dir, exist_ok=True)

# =============================================================================
# Physical Models for Underground Nuclear Stimulation
# =============================================================================

class NuclearStimulationModel:
    """Model for underground nuclear gas stimulation."""
    
    JOULES_PER_KT = 4.184e12
    
    def __init__(self, name, yield_kt, depth, 
                 host_density=2500, porosity=0.10, permeability=0.5,
                 initial_pressure=20e6, initial_temperature=350,
                 cavity_scaling=75.0, chimney_ratio=4.0):
        """
        Initialize stimulation model.
        
        Parameters:
        -----------
        name : str - Project name
        yield_kt : float - Yield in kilotons
        depth : float - Depth of burial (m)
        host_density : float - Rock density (kg/m³)
        porosity : float - Formation porosity
        permeability : float - Initial permeability (mD)
        initial_pressure : float - Reservoir pressure (Pa)
        initial_temperature : float - Reservoir temperature (K)
        cavity_scaling : float - Cavity radius scaling (R = C * W^(1/3))
        chimney_ratio : float - Chimney height / cavity radius
        """
        self.name = name
        self.yield_kt = yield_kt
        self.depth = depth
        self.host_density = host_density
        self.porosity = porosity
        self.permeability = permeability
        self.initial_pressure = initial_pressure
        self.initial_temperature = initial_temperature
        self.cavity_scaling = cavity_scaling
        self.chimney_ratio = chimney_ratio
        
        # Computed properties
        self.total_energy = yield_kt * self.JOULES_PER_KT
        self.cavity_radius = self.compute_cavity_radius()
        self.chimney_height = self.cavity_radius * chimney_ratio
        
    def compute_cavity_radius(self):
        """Compute cavity radius using scaling law R = C * W^(1/3)."""
        # Adjust for overburden pressure
        overburden = self.host_density * 9.81 * self.depth
        # Pressure correction factor
        P_atm = 101325.0
        P_ratio = (overburden + P_atm) / (7.0e6)  # Reference pressure 7 MPa
        correction = P_ratio ** (-0.3)
        
        R = self.cavity_scaling * (self.yield_kt ** (1/3)) * correction
        return R
    
    def crushed_zone_radius(self):
        """Radius of crushed/comminuted zone (typically 2-3x cavity)."""
        return 2.5 * self.cavity_radius
    
    def fractured_zone_radius(self):
        """Radius of fractured zone (typically 4-6x cavity)."""
        return 5.0 * self.cavity_radius
    
    def enhanced_zone_radius(self):
        """Radius of enhanced permeability zone."""
        return 8.0 * self.cavity_radius
    
    def cavity_volume(self):
        """Volume of cavity (m³)."""
        return (4/3) * np.pi * self.cavity_radius**3
    
    def chimney_volume(self):
        """Volume of rubble chimney (m³)."""
        # Approximate as cylinder plus hemisphere
        return np.pi * self.cavity_radius**2 * self.chimney_height + \
               (2/3) * np.pi * self.cavity_radius**3
    
    def enhanced_permeability(self, r, enhancement_factor=100):
        """
        Permeability enhancement as function of distance.
        
        Returns permeability in mD.
        """
        if r < self.cavity_radius:
            # Inside cavity - effectively infinite (rubble)
            return 10000.0  # 10 Darcy
        elif r < self.crushed_zone_radius():
            # Crushed zone - high permeability
            return self.permeability * enhancement_factor
        elif r < self.fractured_zone_radius():
            # Fractured zone - moderate enhancement
            factor = enhancement_factor * np.exp(-(r - self.crushed_zone_radius()) / 
                                                  (self.fractured_zone_radius() - self.crushed_zone_radius()))
            return self.permeability * max(factor, 10)
        elif r < self.enhanced_zone_radius():
            # Enhanced zone - slight increase
            return self.permeability * 5.0
        else:
            return self.permeability
    
    def pressure_evolution(self, t, P0=None, drainage_rate=0.0):
        """
        Cavity pressure evolution over time.
        
        Parameters:
        -----------
        t : array - Time in seconds
        P0 : float - Initial cavity pressure (Pa)
        drainage_rate : float - Production rate (m³/s)
        """
        if P0 is None:
            # Initial cavity pressure from equation of state
            # P * V = const for adiabatic expansion
            initial_volume = (4/3) * np.pi * 10**3  # 10m initial radius
            P0 = self.initial_pressure * (initial_volume / self.cavity_volume())**(-1.4)
        
        # Decay due to heat loss and drainage
        tau_thermal = 86400 * 30  # 30 day thermal time constant
        tau_drainage = self.cavity_volume() / max(drainage_rate, 0.01)
        
        P = P0 * np.exp(-t / tau_thermal)
        if drainage_rate > 0:
            P = P * np.exp(-t / tau_drainage)
        
        # Equilibrate to lithostatic
        P_lith = self.host_density * 9.81 * self.depth
        P = np.maximum(P, P_lith * 0.5)
        
        return P
    
    def production_rate(self, t, P_bhp=5e6):
        """
        Gas production rate over time (m³/day at surface conditions).
        
        Simplified Darcy flow model.
        """
        # Effective permeability of enhanced zone
        k_eff = self.enhanced_permeability(self.crushed_zone_radius()) * 1e-15  # mD to m²
        
        # Average pressure
        P_avg = (self.initial_pressure + P_bhp) / 2
        
        # Gas viscosity
        mu = 2e-5  # Pa·s
        
        # Drainage radius
        r_e = self.enhanced_zone_radius()
        r_w = self.cavity_radius
        
        # Darcy's law (radial flow)
        # Q = 2 * pi * k * h * (P_res - P_bhp) / (mu * ln(r_e/r_w))
        h = self.chimney_height
        
        Q = 2 * np.pi * k_eff * h * (P_avg - P_bhp) / (mu * np.log(r_e / r_w))
        
        # Convert to surface conditions (simplified)
        Q_surface = Q * (P_avg / 101325)  # m³/s
        Q_day = Q_surface * 86400  # m³/day
        
        # Production decline
        b = 0.5  # Decline exponent
        D = 0.001 / 86400  # Initial decline rate (1/s)
        decline = (1 + b * D * t) ** (-1/b)
        
        return Q_day * decline * np.heaviside(t - 90*86400, 0.5)  # Start after 90 days


# =============================================================================
# Create Models for Both Projects
# =============================================================================

# US Project Gasbuggy (1967)
gasbuggy = NuclearStimulationModel(
    name="Project Gasbuggy (US, 1967)",
    yield_kt=29.0,
    depth=1292.0,
    host_density=2450,
    porosity=0.11,
    permeability=0.005,  # Very tight!
    initial_pressure=15e6,
    initial_temperature=320,
    cavity_scaling=70.0,
    chimney_ratio=4.0
)

# Soviet Urta-Bulak style
soviet = NuclearStimulationModel(
    name="Soviet Urta-Bulak Style (1966)",
    yield_kt=30.0,
    depth=1500.0,
    host_density=2400,
    porosity=0.12,
    permeability=0.5,
    initial_pressure=25e6,
    initial_temperature=350,
    cavity_scaling=80.0,
    chimney_ratio=3.5
)

print("=" * 70)
print("Nuclear Reservoir Stimulation Comparison")
print("=" * 70)
print(f"\n{'Parameter':<30} {'Gasbuggy (US)':>18} {'Urta-Bulak (USSR)':>18}")
print("-" * 70)
print(f"{'Yield (kt)':<30} {gasbuggy.yield_kt:>18.1f} {soviet.yield_kt:>18.1f}")
print(f"{'Depth (m)':<30} {gasbuggy.depth:>18.0f} {soviet.depth:>18.0f}")
print(f"{'Cavity Radius (m)':<30} {gasbuggy.cavity_radius:>18.1f} {soviet.cavity_radius:>18.1f}")
print(f"{'Chimney Height (m)':<30} {gasbuggy.chimney_height:>18.1f} {soviet.chimney_height:>18.1f}")
print(f"{'Crushed Zone Radius (m)':<30} {gasbuggy.crushed_zone_radius():>18.1f} {soviet.crushed_zone_radius():>18.1f}")
print(f"{'Fractured Zone Radius (m)':<30} {gasbuggy.fractured_zone_radius():>18.1f} {soviet.fractured_zone_radius():>18.1f}")
print(f"{'Cavity Volume (m³)':<30} {gasbuggy.cavity_volume():>18.0f} {soviet.cavity_volume():>18.0f}")
print(f"{'Initial Permeability (mD)':<30} {gasbuggy.permeability:>18.3f} {soviet.permeability:>18.3f}")
print(f"{'Reservoir Pressure (MPa)':<30} {gasbuggy.initial_pressure/1e6:>18.1f} {soviet.initial_pressure/1e6:>18.1f}")
print("=" * 70)

# =============================================================================
# Plot 1: Overview Comparison
# =============================================================================

fig = plt.figure(figsize=(20, 14))
fig.suptitle('Nuclear Reservoir Stimulation: US vs Soviet Approach', 
             fontsize=16, fontweight='bold', y=0.98)

# Left panel: Cross-section comparison
ax1 = fig.add_subplot(2, 3, 1)
ax1.set_title('Cross-Section: Project Gasbuggy (US)', fontsize=11)

# Draw formation layers
for y in np.arange(0, -2000, -200):
    ax1.axhline(y, color='gray', alpha=0.3, linewidth=0.5)

# Draw cavity and zones
depth = -gasbuggy.depth
r_cav = gasbuggy.cavity_radius
r_crush = gasbuggy.crushed_zone_radius()
r_frac = gasbuggy.fractured_zone_radius()
r_enh = gasbuggy.enhanced_zone_radius()

# Enhanced zone
circle_enh = Circle((0, depth), r_enh, facecolor='lightgreen', alpha=0.3, 
                    edgecolor='green', linestyle='--', label='Enhanced zone')
ax1.add_patch(circle_enh)

# Fractured zone
circle_frac = Circle((0, depth), r_frac, facecolor='yellow', alpha=0.4, 
                     edgecolor='orange', linestyle='--', label='Fractured zone')
ax1.add_patch(circle_frac)

# Crushed zone
circle_crush = Circle((0, depth), r_crush, facecolor='orange', alpha=0.5, 
                      edgecolor='red', linestyle='-', label='Crushed zone')
ax1.add_patch(circle_crush)

# Cavity
circle_cav = Circle((0, depth), r_cav, facecolor='red', alpha=0.7, 
                   edgecolor='darkred', linewidth=2, label='Cavity')
ax1.add_patch(circle_cav)

# Chimney
chimney_top = depth + gasbuggy.chimney_height
chimney = plt.Rectangle((-r_cav*0.8, depth), r_cav*1.6, gasbuggy.chimney_height,
                        facecolor='brown', alpha=0.5, edgecolor='black', 
                        linestyle=':', label='Rubble chimney')
ax1.add_patch(chimney)

# Well
ax1.plot([0, 0], [0, depth - r_cav], 'k-', linewidth=3, label='Well')
ax1.plot(0, 0, 'ks', markersize=10)

# Surface
ax1.axhline(0, color='brown', linewidth=3)
ax1.fill_between([-800, 800], [50, 50], [0, 0], color='green', alpha=0.3)

ax1.set_xlim(-800, 800)
ax1.set_ylim(-1800, 100)
ax1.set_xlabel('Distance (m)')
ax1.set_ylabel('Depth (m)')
ax1.legend(loc='lower right', fontsize=8)
ax1.set_aspect('equal')

# Right panel: Soviet approach
ax2 = fig.add_subplot(2, 3, 2)
ax2.set_title('Cross-Section: Urta-Bulak (Soviet)', fontsize=11)

for y in np.arange(0, -2000, -200):
    ax2.axhline(y, color='gray', alpha=0.3, linewidth=0.5)

depth = -soviet.depth
r_cav = soviet.cavity_radius
r_crush = soviet.crushed_zone_radius()
r_frac = soviet.fractured_zone_radius()
r_enh = soviet.enhanced_zone_radius()

circle_enh = Circle((0, depth), r_enh, facecolor='lightgreen', alpha=0.3, 
                    edgecolor='green', linestyle='--')
ax2.add_patch(circle_enh)

circle_frac = Circle((0, depth), r_frac, facecolor='yellow', alpha=0.4, 
                     edgecolor='orange', linestyle='--')
ax2.add_patch(circle_frac)

circle_crush = Circle((0, depth), r_crush, facecolor='orange', alpha=0.5, 
                      edgecolor='red', linestyle='-')
ax2.add_patch(circle_crush)

circle_cav = Circle((0, depth), r_cav, facecolor='red', alpha=0.7, 
                   edgecolor='darkred', linewidth=2)
ax2.add_patch(circle_cav)

chimney = plt.Rectangle((-r_cav*0.8, depth), r_cav*1.6, soviet.chimney_height,
                        facecolor='brown', alpha=0.5, edgecolor='black', linestyle=':')
ax2.add_patch(chimney)

ax2.plot([0, 0], [0, depth - r_cav], 'k-', linewidth=3)
ax2.plot(0, 0, 'ks', markersize=10)

ax2.axhline(0, color='brown', linewidth=3)
ax2.fill_between([-800, 800], [50, 50], [0, 0], color='green', alpha=0.3)

ax2.set_xlim(-800, 800)
ax2.set_ylim(-2000, 100)
ax2.set_xlabel('Distance (m)')
ax2.set_ylabel('Depth (m)')
ax2.set_aspect('equal')

# Zone size comparison bar chart
ax3 = fig.add_subplot(2, 3, 3)
ax3.set_title('Zone Size Comparison', fontsize=11)

categories = ['Cavity\nRadius', 'Crushed\nZone', 'Fractured\nZone', 'Enhanced\nZone']
gasbuggy_values = [gasbuggy.cavity_radius, gasbuggy.crushed_zone_radius(),
                   gasbuggy.fractured_zone_radius(), gasbuggy.enhanced_zone_radius()]
soviet_values = [soviet.cavity_radius, soviet.crushed_zone_radius(),
                 soviet.fractured_zone_radius(), soviet.enhanced_zone_radius()]

x = np.arange(len(categories))
width = 0.35

bars1 = ax3.bar(x - width/2, gasbuggy_values, width, label='Gasbuggy (US)', color='blue', alpha=0.7)
bars2 = ax3.bar(x + width/2, soviet_values, width, label='Urta-Bulak (USSR)', color='red', alpha=0.7)

ax3.set_ylabel('Radius (m)')
ax3.set_xticks(x)
ax3.set_xticklabels(categories)
ax3.legend()
ax3.grid(True, alpha=0.3)

# Add value labels
for bar in bars1:
    height = bar.get_height()
    ax3.annotate(f'{height:.0f}', xy=(bar.get_x() + bar.get_width()/2, height),
                xytext=(0, 3), textcoords='offset points', ha='center', va='bottom', fontsize=8)
for bar in bars2:
    height = bar.get_height()
    ax3.annotate(f'{height:.0f}', xy=(bar.get_x() + bar.get_width()/2, height),
                xytext=(0, 3), textcoords='offset points', ha='center', va='bottom', fontsize=8)

# Permeability enhancement profile
ax4 = fig.add_subplot(2, 3, 4)
ax4.set_title('Permeability Enhancement vs Distance', fontsize=11)

r = np.linspace(1, 1000, 500)
perm_gasbuggy = [gasbuggy.enhanced_permeability(ri) for ri in r]
perm_soviet = [soviet.enhanced_permeability(ri) for ri in r]

ax4.semilogy(r, perm_gasbuggy, 'b-', linewidth=2, label='Gasbuggy (US)')
ax4.semilogy(r, perm_soviet, 'r-', linewidth=2, label='Urta-Bulak (USSR)')

# Mark zone boundaries
ax4.axvline(gasbuggy.cavity_radius, color='blue', linestyle=':', alpha=0.7)
ax4.axvline(gasbuggy.crushed_zone_radius(), color='blue', linestyle='--', alpha=0.7)
ax4.axvline(gasbuggy.fractured_zone_radius(), color='blue', linestyle='-.', alpha=0.7)

ax4.axvline(soviet.cavity_radius, color='red', linestyle=':', alpha=0.7)
ax4.axvline(soviet.crushed_zone_radius(), color='red', linestyle='--', alpha=0.7)
ax4.axvline(soviet.fractured_zone_radius(), color='red', linestyle='-.', alpha=0.7)

ax4.axhline(1.0, color='gray', linestyle='-', alpha=0.5, label='1 mD threshold')
ax4.set_xlabel('Distance from Detonation Point (m)')
ax4.set_ylabel('Permeability (mD)')
ax4.legend(loc='upper right')
ax4.grid(True, alpha=0.3)
ax4.set_xlim(0, 1000)
ax4.set_ylim(0.001, 20000)

# Production rate comparison
ax5 = fig.add_subplot(2, 3, 5)
ax5.set_title('Gas Production Rate Over Time', fontsize=11)

t_days = np.linspace(0, 730, 500)  # 2 years
t_seconds = t_days * 86400

prod_gasbuggy = [gasbuggy.production_rate(ti) for ti in t_seconds]
prod_soviet = [soviet.production_rate(ti) for ti in t_seconds]

ax5.plot(t_days, np.array(prod_gasbuggy)/1000, 'b-', linewidth=2, label='Gasbuggy (US)')
ax5.plot(t_days, np.array(prod_soviet)/1000, 'r-', linewidth=2, label='Urta-Bulak (USSR)')

ax5.axvline(90, color='gray', linestyle='--', alpha=0.7, label='Production start')
ax5.axvline(150, color='blue', linestyle=':', alpha=0.7, label='Gasbuggy start (actual)')

ax5.set_xlabel('Days After Detonation')
ax5.set_ylabel('Production Rate (1000 m³/day)')
ax5.legend(loc='upper right')
ax5.grid(True, alpha=0.3)
ax5.set_xlim(0, 730)

# Cavity pressure evolution
ax6 = fig.add_subplot(2, 3, 6)
ax6.set_title('Cavity Pressure Evolution', fontsize=11)

P_gasbuggy = gasbuggy.pressure_evolution(t_seconds, P0=100e6)
P_soviet = soviet.pressure_evolution(t_seconds, P0=150e6)

ax6.semilogy(t_days, P_gasbuggy/1e6, 'b-', linewidth=2, label='Gasbuggy (US)')
ax6.semilogy(t_days, P_soviet/1e6, 'r-', linewidth=2, label='Urta-Bulak (USSR)')

ax6.axhline(gasbuggy.initial_pressure/1e6, color='blue', linestyle='--', alpha=0.5, 
           label=f'Gasbuggy reservoir P ({gasbuggy.initial_pressure/1e6:.0f} MPa)')
ax6.axhline(soviet.initial_pressure/1e6, color='red', linestyle='--', alpha=0.5,
           label=f'Soviet reservoir P ({soviet.initial_pressure/1e6:.0f} MPa)')

ax6.set_xlabel('Days After Detonation')
ax6.set_ylabel('Cavity Pressure (MPa)')
ax6.legend(loc='upper right', fontsize=8)
ax6.grid(True, alpha=0.3)
ax6.set_xlim(0, 730)

plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'nuclear_stimulation_comparison.png'), dpi=150, bbox_inches='tight')
print(f"Saved: {output_dir}/nuclear_stimulation_comparison.png")

# =============================================================================
# Plot 2: Detailed Damage Zone Analysis
# =============================================================================

fig2, axes = plt.subplots(2, 2, figsize=(16, 14))
fig2.suptitle('Nuclear Stimulation: Damage Zone Analysis', fontsize=14, fontweight='bold')

# 2D damage map - Gasbuggy
ax = axes[0, 0]
ax.set_title('Gasbuggy: Damage Zone Map (Plan View)', fontsize=11)

# Create 2D grid
x = np.linspace(-800, 800, 200)
y = np.linspace(-800, 800, 200)
X, Y = np.meshgrid(x, y)
R = np.sqrt(X**2 + Y**2)

# Damage levels
damage = np.zeros_like(R)
damage[R < gasbuggy.enhanced_zone_radius()] = 0.2
damage[R < gasbuggy.fractured_zone_radius()] = 0.5
damage[R < gasbuggy.crushed_zone_radius()] = 0.8
damage[R < gasbuggy.cavity_radius] = 1.0

cmap = mcolors.LinearSegmentedColormap.from_list('damage', 
    ['white', 'lightgreen', 'yellow', 'orange', 'red'])
im = ax.contourf(X, Y, damage, levels=[0, 0.1, 0.3, 0.6, 0.9, 1.0], cmap=cmap)
ax.contour(X, Y, damage, levels=[0.2, 0.5, 0.8], colors='black', linewidths=1)

# Add fracture pattern
np.random.seed(42)
n_fractures = 16
for i in range(n_fractures):
    angle = 2 * np.pi * i / n_fractures + np.random.uniform(-0.1, 0.1)
    r_start = gasbuggy.cavity_radius
    r_end = gasbuggy.fractured_zone_radius() * np.random.uniform(0.7, 1.0)
    x_start = r_start * np.cos(angle)
    y_start = r_start * np.sin(angle)
    x_end = r_end * np.cos(angle)
    y_end = r_end * np.sin(angle)
    ax.plot([x_start, x_end], [y_start, y_end], 'k-', linewidth=0.5, alpha=0.5)

ax.set_xlabel('X Distance (m)')
ax.set_ylabel('Y Distance (m)')
ax.set_aspect('equal')
plt.colorbar(im, ax=ax, label='Damage Level')

# 2D damage map - Soviet
ax = axes[0, 1]
ax.set_title('Urta-Bulak: Damage Zone Map (Plan View)', fontsize=11)

damage_s = np.zeros_like(R)
damage_s[R < soviet.enhanced_zone_radius()] = 0.2
damage_s[R < soviet.fractured_zone_radius()] = 0.5
damage_s[R < soviet.crushed_zone_radius()] = 0.8
damage_s[R < soviet.cavity_radius] = 1.0

im = ax.contourf(X, Y, damage_s, levels=[0, 0.1, 0.3, 0.6, 0.9, 1.0], cmap=cmap)
ax.contour(X, Y, damage_s, levels=[0.2, 0.5, 0.8], colors='black', linewidths=1)

# Add more fractures (Soviet style)
n_fractures = 24
for i in range(n_fractures):
    angle = 2 * np.pi * i / n_fractures + np.random.uniform(-0.1, 0.1)
    r_start = soviet.cavity_radius
    r_end = soviet.fractured_zone_radius() * np.random.uniform(0.6, 1.0)
    x_start = r_start * np.cos(angle)
    y_start = r_start * np.sin(angle)
    x_end = r_end * np.cos(angle)
    y_end = r_end * np.sin(angle)
    ax.plot([x_start, x_end], [y_start, y_end], 'k-', linewidth=0.5, alpha=0.5)

ax.set_xlabel('X Distance (m)')
ax.set_ylabel('Y Distance (m)')
ax.set_aspect('equal')
plt.colorbar(im, ax=ax, label='Damage Level')

# Radial permeability profile
ax = axes[1, 0]
ax.set_title('Radial Permeability Profile', fontsize=11)

r = np.linspace(1, 1000, 500)

for model, color, ls in [(gasbuggy, 'blue', '-'), (soviet, 'red', '-')]:
    perm = [model.enhanced_permeability(ri) for ri in r]
    ax.semilogy(r, perm, color=color, linestyle=ls, linewidth=2, label=model.name)
    
    # Mark zones
    ax.axvline(model.cavity_radius, color=color, linestyle=':', alpha=0.5)
    ax.axvline(model.crushed_zone_radius(), color=color, linestyle='--', alpha=0.5)
    ax.axvline(model.fractured_zone_radius(), color=color, linestyle='-.', alpha=0.5)

# Pre-shot permeabilities
ax.axhline(gasbuggy.permeability, color='blue', linestyle=':', linewidth=1,
          label=f'Gasbuggy pre-shot ({gasbuggy.permeability} mD)')
ax.axhline(soviet.permeability, color='red', linestyle=':', linewidth=1,
          label=f'Soviet pre-shot ({soviet.permeability} mD)')

ax.set_xlabel('Distance from Detonation (m)')
ax.set_ylabel('Permeability (mD)')
ax.legend(loc='upper right', fontsize=8)
ax.grid(True, alpha=0.3)
ax.set_xlim(0, 1000)
ax.set_ylim(1e-3, 2e4)

# Volume comparison
ax = axes[1, 1]
ax.set_title('Affected Volume Comparison', fontsize=11)

# Calculate volumes
def sphere_vol(r):
    return (4/3) * np.pi * r**3

volumes_gb = [
    gasbuggy.cavity_volume() / 1e6,
    (sphere_vol(gasbuggy.crushed_zone_radius()) - gasbuggy.cavity_volume()) / 1e6,
    (sphere_vol(gasbuggy.fractured_zone_radius()) - sphere_vol(gasbuggy.crushed_zone_radius())) / 1e6,
    (sphere_vol(gasbuggy.enhanced_zone_radius()) - sphere_vol(gasbuggy.fractured_zone_radius())) / 1e6
]

volumes_sov = [
    soviet.cavity_volume() / 1e6,
    (sphere_vol(soviet.crushed_zone_radius()) - soviet.cavity_volume()) / 1e6,
    (sphere_vol(soviet.fractured_zone_radius()) - sphere_vol(soviet.crushed_zone_radius())) / 1e6,
    (sphere_vol(soviet.enhanced_zone_radius()) - sphere_vol(soviet.fractured_zone_radius())) / 1e6
]

categories = ['Cavity', 'Crushed\nZone', 'Fractured\nZone', 'Enhanced\nZone']
x = np.arange(len(categories))
width = 0.35

bars1 = ax.bar(x - width/2, volumes_gb, width, label='Gasbuggy (US)', color='blue', alpha=0.7)
bars2 = ax.bar(x + width/2, volumes_sov, width, label='Urta-Bulak (USSR)', color='red', alpha=0.7)

ax.set_ylabel('Volume (million m³)')
ax.set_xticks(x)
ax.set_xticklabels(categories)
ax.legend()
ax.grid(True, alpha=0.3, axis='y')
ax.set_yscale('log')

plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'nuclear_stimulation_damage_zones.png'), dpi=150, bbox_inches='tight')
print(f"Saved: {output_dir}/nuclear_stimulation_damage_zones.png")

# =============================================================================
# Plot 3: Time Evolution and Production Analysis
# =============================================================================

fig3 = plt.figure(figsize=(18, 12))
fig3.suptitle('Nuclear Stimulation: Production Analysis', fontsize=14, fontweight='bold')

# Production rate over 2 years
ax1 = fig3.add_subplot(2, 3, 1)
ax1.set_title('Production Rate Evolution', fontsize=11)

t_days = np.linspace(0, 730, 500)
t_seconds = t_days * 86400

prod_gb = [gasbuggy.production_rate(ti) for ti in t_seconds]
prod_sov = [soviet.production_rate(ti) for ti in t_seconds]

ax1.plot(t_days, np.array(prod_gb)/1000, 'b-', linewidth=2, label='Gasbuggy (US)')
ax1.plot(t_days, np.array(prod_sov)/1000, 'r-', linewidth=2, label='Urta-Bulak (USSR)')

# Actual Gasbuggy data point (approx)
ax1.scatter([180], [70], color='blue', s=100, marker='*', zorder=5, 
           label='Gasbuggy actual (~70,000 m³/day)')

ax1.set_xlabel('Days After Detonation')
ax1.set_ylabel('Production Rate (1000 m³/day)')
ax1.legend()
ax1.grid(True, alpha=0.3)
ax1.set_xlim(0, 730)

# Cumulative production
ax2 = fig3.add_subplot(2, 3, 2)
ax2.set_title('Cumulative Production', fontsize=11)

dt = t_seconds[1] - t_seconds[0]
cum_gb = np.cumsum(prod_gb) * dt / 86400 / 1e6  # Million m³
cum_sov = np.cumsum(prod_sov) * dt / 86400 / 1e6

ax2.plot(t_days, cum_gb, 'b-', linewidth=2, label='Gasbuggy (US)')
ax2.plot(t_days, cum_sov, 'r-', linewidth=2, label='Urta-Bulak (USSR)')

ax2.set_xlabel('Days After Detonation')
ax2.set_ylabel('Cumulative Production (million m³)')
ax2.legend()
ax2.grid(True, alpha=0.3)

# Productivity index comparison
ax3 = fig3.add_subplot(2, 3, 3)
ax3.set_title('Enhancement Factor vs Pre-Shot', fontsize=11)

# Estimate pre-shot production (conventional)
pre_shot_gb = gasbuggy.production_rate(365*86400) * 0.001  # Very low
pre_shot_sov = soviet.production_rate(365*86400) * 0.02

# Enhancement factors
days = [90, 180, 365, 540, 730]
enhancement_gb = []
enhancement_sov = []

for d in days:
    t = d * 86400
    enhancement_gb.append(gasbuggy.production_rate(t) / max(pre_shot_gb, 1))
    enhancement_sov.append(soviet.production_rate(t) / max(pre_shot_sov, 1))

ax3.bar([d - 15 for d in days], enhancement_gb, width=25, label='Gasbuggy (US)', color='blue', alpha=0.7)
ax3.bar([d + 15 for d in days], enhancement_sov, width=25, label='Urta-Bulak (USSR)', color='red', alpha=0.7)

ax3.set_xlabel('Days After Detonation')
ax3.set_ylabel('Production Enhancement Factor')
ax3.legend()
ax3.grid(True, alpha=0.3)
ax3.set_xticks(days)

# Tritium contamination (US concern)
ax4 = fig3.add_subplot(2, 3, 4)
ax4.set_title('Tritium Contamination Decay (Gasbuggy)', fontsize=11)

t_years = np.linspace(0, 50, 200)
t_halflife = 12.32  # years

# Initial tritium activity (Bq/m³ in produced gas)
T_initial = 1e6  # Bq/m³ (measured peak)
T_activity = T_initial * (0.5 ** (t_years / t_halflife))

ax4.semilogy(t_years, T_activity, 'b-', linewidth=2)
ax4.axhline(740, color='red', linestyle='--', label='Commercial limit (740 Bq/m³)')
ax4.axhline(20, color='green', linestyle='--', label='Drinking water limit (20 Bq/L)')

ax4.fill_between(t_years, T_activity, 740, where=T_activity > 740,
                 color='red', alpha=0.3, label='Above commercial limit')

ax4.set_xlabel('Years After Detonation')
ax4.set_ylabel('Tritium Activity (Bq/m³)')
ax4.legend(loc='upper right')
ax4.grid(True, alpha=0.3)
ax4.set_xlim(0, 50)
ax4.set_ylim(1, 2e6)

# Time to reach commercial limit
time_to_limit = t_halflife * np.log2(T_initial / 740)
ax4.axvline(time_to_limit, color='purple', linestyle=':', 
           label=f'~{time_to_limit:.0f} years to limit')
ax4.text(time_to_limit + 1, 1e5, f'{time_to_limit:.0f} years', fontsize=10)

# Comparison table as text
ax5 = fig3.add_subplot(2, 3, 5)
ax5.axis('off')
ax5.set_title('Project Comparison Summary', fontsize=11)

table_text = """
Parameter                    Gasbuggy (US)     Urta-Bulak (USSR)
─────────────────────────────────────────────────────────────────
Year                              1967                1966
Yield (kt)                        29.0                30.0
Depth (m)                         1292                1500
Cavity Radius (m)                 ~80                 ~250
Chimney Height (m)                ~335                ~875

Formation                    Pictured Cliffs    Gas-bearing sand
Initial Permeability (mD)         0.005               0.5
Enhanced Permeability (mD)        ~250              ~5000

Gas Quality                  Contaminated!        Usable
Tritium Level               1350× limit         Not measured
Production Increase              ~4×                 ~10×

Program Outcome             Abandoned 1977     Continued to 1988
─────────────────────────────────────────────────────────────────
Key Difference: US prioritized safety → cancellation
                USSR prioritized production → continued use
"""

ax5.text(0.05, 0.95, table_text, transform=ax5.transAxes, fontsize=9,
         family='monospace', verticalalignment='top')

# Energy partitioning pie chart
ax6 = fig3.add_subplot(2, 3, 6)
ax6.set_title('Energy Partitioning (Typical Underground Nuclear)', fontsize=11)

# Typical energy partition for contained nuclear explosion
labels = ['Thermal\n(vaporization)', 'Mechanical\n(shock)', 'Seismic', 
          'Cavity work', 'Radiation']
sizes = [40, 35, 5, 15, 5]
colors = ['red', 'orange', 'green', 'brown', 'purple']
explode = (0.05, 0.05, 0, 0, 0)

ax6.pie(sizes, explode=explode, labels=labels, colors=colors, autopct='%1.0f%%',
        shadow=True, startangle=90)
ax6.axis('equal')

plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'nuclear_stimulation_production.png'), dpi=150, bbox_inches='tight')
print(f"Saved: {output_dir}/nuclear_stimulation_production.png")

# =============================================================================
# Plot 4: 3D Visualization of Cavity and Zones
# =============================================================================

fig4 = plt.figure(figsize=(16, 8))
fig4.suptitle('Nuclear Stimulation: 3D Zone Visualization', fontsize=14, fontweight='bold')

for idx, (model, title) in enumerate([(gasbuggy, 'Project Gasbuggy'), 
                                       (soviet, 'Urta-Bulak')]):
    ax = fig4.add_subplot(1, 2, idx+1, projection='3d')
    ax.set_title(title, fontsize=11)
    
    # Create spherical zones
    u = np.linspace(0, 2 * np.pi, 30)
    v = np.linspace(0, np.pi, 20)
    
    # Enhanced zone (outer, transparent)
    r = model.enhanced_zone_radius()
    x = r * np.outer(np.cos(u), np.sin(v))
    y = r * np.outer(np.sin(u), np.sin(v))
    z = r * np.outer(np.ones(np.size(u)), np.cos(v)) - model.depth
    ax.plot_surface(x, y, z, alpha=0.1, color='green')
    
    # Fractured zone
    r = model.fractured_zone_radius()
    x = r * np.outer(np.cos(u), np.sin(v))
    y = r * np.outer(np.sin(u), np.sin(v))
    z = r * np.outer(np.ones(np.size(u)), np.cos(v)) - model.depth
    ax.plot_surface(x, y, z, alpha=0.2, color='yellow')
    
    # Crushed zone
    r = model.crushed_zone_radius()
    x = r * np.outer(np.cos(u), np.sin(v))
    y = r * np.outer(np.sin(u), np.sin(v))
    z = r * np.outer(np.ones(np.size(u)), np.cos(v)) - model.depth
    ax.plot_surface(x, y, z, alpha=0.3, color='orange')
    
    # Cavity
    r = model.cavity_radius
    x = r * np.outer(np.cos(u), np.sin(v))
    y = r * np.outer(np.sin(u), np.sin(v))
    z = r * np.outer(np.ones(np.size(u)), np.cos(v)) - model.depth
    ax.plot_surface(x, y, z, alpha=0.7, color='red')
    
    # Chimney (cylinder approximation)
    theta = np.linspace(0, 2*np.pi, 30)
    z_cyl = np.linspace(-model.depth + model.cavity_radius, 
                        -model.depth + model.chimney_height, 20)
    Theta, Z_cyl = np.meshgrid(theta, z_cyl)
    X_cyl = model.cavity_radius * 0.8 * np.cos(Theta)
    Y_cyl = model.cavity_radius * 0.8 * np.sin(Theta)
    ax.plot_surface(X_cyl, Y_cyl, Z_cyl, alpha=0.4, color='brown')
    
    # Well
    ax.plot([0, 0], [0, 0], [0, -model.depth + model.cavity_radius], 'k-', linewidth=3)
    
    # Surface
    xx, yy = np.meshgrid(np.linspace(-800, 800, 10), np.linspace(-800, 800, 10))
    ax.plot_surface(xx, yy, np.zeros_like(xx), alpha=0.3, color='brown')
    
    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    ax.set_zlabel('Depth (m)')
    ax.set_xlim(-800, 800)
    ax.set_ylim(-800, 800)
    ax.set_zlim(-model.depth - 500, 200)

plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'nuclear_stimulation_3d.png'), dpi=150, bbox_inches='tight')
print(f"Saved: {output_dir}/nuclear_stimulation_3d.png")

plt.close('all')

print("\n" + "=" * 70)
print("Nuclear Stimulation Visualization Complete!")
print("=" * 70)
print(f"\nOutput files saved to: {output_dir}/")
print("\nGenerated plots:")
print("  1. nuclear_stimulation_comparison.png - US vs Soviet overview")
print("  2. nuclear_stimulation_damage_zones.png - Detailed damage analysis")
print("  3. nuclear_stimulation_production.png - Production and contamination")
print("  4. nuclear_stimulation_3d.png - 3D zone visualization")
print("=" * 70)
