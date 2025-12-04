#!/usr/bin/env python3
"""
Comprehensive Visualization for Atmospheric Nuclear Test Simulation
====================================================================

This script visualizes results from the atmospheric_nuclear_test.config simulation,
including blast wave propagation, thermal radiation, EMP fields, radiation dose,
ground motion, and fallout patterns.

Usage:
    python visualize_atmospheric_nuclear.py [output_directory] [--options]
    
Example:
    python visualize_atmospheric_nuclear.py output/atmospheric_nuclear --save-all

Requirements:
    pip install numpy matplotlib h5py scipy pyvista imageio
"""

import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Circle, Wedge
from matplotlib.collections import PatchCollection
from matplotlib.animation import FuncAnimation, PillowWriter
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import RegularGridInterpolator, interp1d
from scipy.ndimage import gaussian_filter
import warnings
warnings.filterwarnings('ignore')

# Try to import optional dependencies
try:
    import h5py
    HAS_HDF5 = True
except ImportError:
    HAS_HDF5 = False
    print("Warning: h5py not installed. Using synthetic data for demonstration.")

try:
    import pyvista as pv
    HAS_PYVISTA = True
except ImportError:
    HAS_PYVISTA = False

try:
    import imageio
    HAS_IMAGEIO = True
except ImportError:
    HAS_IMAGEIO = False


# =============================================================================
# Physical Constants and Scaling Laws
# =============================================================================

class NuclearExplosionPhysics:
    """Physical models for nuclear explosion phenomenology."""
    
    # Energy conversion
    JOULES_PER_KT = 4.184e12
    
    def __init__(self, yield_kt=500.0, burst_height=500.0):
        self.yield_kt = yield_kt
        self.burst_height = burst_height
        self.total_energy = yield_kt * self.JOULES_PER_KT
        
    def fireball_radius(self, t):
        """Maximum fireball radius (meters) using Brode model."""
        # R_max = 66 * W^0.4 meters for W in kt
        R_max = 66.0 * (self.yield_kt ** 0.4)
        
        # Time-dependent expansion
        t_formation = 0.001 * (self.yield_kt ** 0.3)  # Formation time (s)
        
        if t < 0:
            return 0.0
        elif t < t_formation:
            return R_max * (t / t_formation) ** 0.4
        else:
            return R_max
    
    def shock_radius(self, t):
        """Shock front radius using Sedov-Taylor solution."""
        if t <= 0:
            return 0.0
        
        # Ambient density
        rho_0 = 1.225  # kg/m³ at sea level
        
        # Sedov-Taylor: R = (E/rho)^(1/5) * t^(2/5)
        E = self.total_energy
        R = 1.15 * (E / rho_0) ** 0.2 * t ** 0.4
        
        return R
    
    def peak_overpressure(self, r, height=None):
        """Peak overpressure (Pa) at ground range r."""
        if height is None:
            height = self.burst_height
        
        # Handle both scalar and array inputs
        r = np.atleast_1d(np.asarray(r))
        scalar_input = r.shape == (1,)
            
        # Slant range
        R = np.sqrt(r**2 + height**2)
        
        # Scaled distance (m/kt^1/3)
        Z = R / (self.yield_kt ** (1/3))
        
        # Glasstone-Dolan empirical fit (psi)
        P_psi = np.where(Z < 10, 
                        1e6,  # Near source, very high
                        1.772e5 / Z**3 + 1.41e4 / Z**2 + 5.0e2 / Z)
        
        # Convert psi to Pa
        P_Pa = P_psi * 6894.76
        
        # Ground reflection enhancement (Mach stem)
        # Regular to Mach reflection transition
        angle = np.arctan2(height, r + 1e-10)  # Avoid division by zero
        mach_factor = np.where(angle < np.radians(40), 1.8, 1.0)
        P_Pa = P_Pa * mach_factor
        
        if scalar_input:
            return float(P_Pa[0])
        return P_Pa
    
    def dynamic_pressure(self, overpressure):
        """Dynamic pressure from overpressure."""
        P0 = 101325.0  # Pa
        gamma = 1.4
        
        # q = 5/2 * p^2 / (7*P0 + p)
        q = 2.5 * overpressure**2 / (7 * P0 + overpressure)
        return q
    
    def thermal_fluence(self, r, height=None):
        """Thermal fluence (J/m²) at ground range r."""
        if height is None:
            height = self.burst_height
        
        r = np.atleast_1d(np.asarray(r))
        scalar_input = r.shape == (1,)
            
        # Slant range
        R = np.sqrt(r**2 + height**2)
        
        # Total thermal energy (35% of yield)
        E_thermal = 0.35 * self.total_energy
        
        # Fluence at slant range (ignoring atmospheric absorption)
        Q = E_thermal / (4 * np.pi * R**2)
        
        # Atmospheric transmission (simplified)
        tau = np.exp(-R / 20000.0)  # ~20 km visibility
        
        result = Q * tau
        if scalar_input:
            return float(result[0])
        return result
    
    def prompt_radiation_dose(self, r, height=None):
        """Prompt radiation dose (Gy) at ground range r."""
        if height is None:
            height = self.burst_height
        
        r = np.atleast_1d(np.asarray(r))
        scalar_input = r.shape == (1,)
            
        R = np.sqrt(r**2 + height**2)
        
        # Gamma dose (simplified)
        # About 5% of yield as prompt radiation, mostly gamma
        E_rad = 0.05 * self.total_energy
        
        # Reference dose at 1 km for 1 kt
        D_ref = 10.0  # Gy at 1 km for 1 kt
        
        # Scale with yield and distance
        D = D_ref * self.yield_kt * (1000.0 / R)**2
        
        # Atmospheric attenuation
        D = D * np.exp(-R / 3000.0)
        
        if scalar_input:
            return float(D[0])
        return D
    
    def emp_e1_field(self, r, t):
        """E1 component of EMP (V/m)."""
        # Rise time ~2.5 ns, decay ~1 µs
        tau_rise = 2.5e-9
        tau_decay = 1.0e-6
        
        # Peak field (depends on yield and geometry)
        E_peak = 50000.0 * (self.yield_kt / 500.0) ** 0.5
        
        # Distance dependence
        E_peak *= np.exp(-r / 50000.0)
        
        # Time function (double exponential)
        if t <= 0:
            return 0.0
        
        E = E_peak * (np.exp(-t / tau_decay) - np.exp(-t / tau_rise))
        return E
    
    def emp_e3_field(self, r, t):
        """E3 component of EMP (V/m)."""
        # Much slower: rise ~1 s, duration ~100 s
        tau_rise = 1.0
        tau_decay = 100.0
        
        # Lower peak
        E_peak = 10.0 * (self.yield_kt / 500.0) ** 0.3
        
        if t <= 0:
            return 0.0
        
        E = E_peak * (1 - np.exp(-t / tau_rise)) * np.exp(-t / tau_decay)
        return E


# =============================================================================
# Data Loading Functions
# =============================================================================

def generate_synthetic_data(physics, nx=200, ny=200, nz=100, nt=100):
    """Generate synthetic data for visualization demonstration."""
    
    print("Generating synthetic data for demonstration...")
    
    # Domain
    x = np.linspace(0, 100000, nx)  # 100 km
    y = np.linspace(0, 100000, ny)
    z = np.linspace(-50000, 100000, nz)  # -50 km to +100 km
    t = np.linspace(0, 60, nt)  # 60 seconds
    
    # Center at (50000, 50000)
    x_burst, y_burst = 50000, 50000
    
    X, Y = np.meshgrid(x, y)
    r = np.sqrt((X - x_burst)**2 + (Y - y_burst)**2)
    
    data = {}
    
    # Time snapshots
    data['time'] = t
    data['x'] = x
    data['y'] = y
    data['z'] = z
    
    # Overpressure field (time-varying)
    overpressure = np.zeros((nt, ny, nx))
    for i, ti in enumerate(t):
        shock_r = physics.shock_radius(ti)
        if shock_r < 1.0:
            shock_r = 1.0  # Minimum shock radius
        # Create shock front
        shock_width = 500.0  # meters
        front = np.exp(-((r - shock_r) / shock_width)**2)
        
        # Overpressure behind shock (vectorized)
        P = np.zeros_like(r)
        mask = r < shock_r
        if np.any(mask):
            r_masked = r[mask]
            P_masked = physics.peak_overpressure(r_masked.flatten())
            if isinstance(P_masked, np.ndarray):
                P[mask] = P_masked * (1 - r_masked / (shock_r + 1))
            else:
                P[mask] = P_masked * (1 - r_masked / (shock_r + 1))
        
        # Add shock front
        P_shock = physics.peak_overpressure(shock_r)
        P += front * P_shock * 0.5
        
        overpressure[i] = P / 1000.0  # Convert to kPa
    
    data['overpressure'] = overpressure
    
    # Temperature field (fireball)
    temperature = np.zeros((nt, ny, nx))
    for i, ti in enumerate(t):
        R_fb = physics.fireball_radius(ti)
        T_surface = 8000.0 if ti > 0.001 else 1e6
        if ti > 1.0:
            T_surface = 8000.0 * (1.0 / ti)**0.5
        
        # Temperature decreases with distance from fireball
        T = T_surface * np.exp(-((r - R_fb) / 1000.0)**2)
        T = np.maximum(T, 288.15)  # Ambient
        temperature[i] = T
    
    data['temperature'] = temperature
    
    # Thermal fluence (cumulative) - vectorized
    thermal_fluence = physics.thermal_fluence(r.flatten()).reshape(r.shape)
    data['thermal_fluence'] = thermal_fluence / 1000.0  # kJ/m²
    
    # Radiation dose - vectorized
    radiation_dose = physics.prompt_radiation_dose(r.flatten()).reshape(r.shape)
    data['radiation_dose'] = radiation_dose
    
    # EMP fields
    emp_e1 = np.zeros((nt,))
    emp_e3 = np.zeros((nt,))
    for i, ti in enumerate(t):
        emp_e1[i] = physics.emp_e1_field(10000, ti)
        emp_e3[i] = physics.emp_e3_field(10000, ti)
    data['emp_e1'] = emp_e1
    data['emp_e3'] = emp_e3
    
    # EMP spatial distribution at t=1µs
    emp_spatial = np.zeros((ny, nx))
    for i in range(ny):
        for j in range(nx):
            emp_spatial[i, j] = physics.emp_e1_field(float(r[i, j]), 1e-6)
    data['emp_spatial'] = emp_spatial
    
    # Ground velocity (seismic)
    ground_velocity = np.zeros((nt, ny, nx))
    for i, ti in enumerate(t):
        # P-wave arrives first
        vp = 6000.0  # m/s
        vs = 3500.0
        
        # P-wave front
        rp = vp * ti
        front_p = np.exp(-((r - rp) / 500.0)**2)
        
        # S-wave front
        rs = vs * ti
        front_s = np.exp(-((r - rs) / 300.0)**2)
        
        # Combined with attenuation
        v = (front_p * 0.1 + front_s * 0.2) * np.exp(-r / 50000.0)
        ground_velocity[i] = v
    
    data['ground_velocity'] = ground_velocity
    
    # Receiver time series
    receiver_distances = [1000, 5000, 10000, 20000, 50000]  # meters
    receiver_data = {}
    for d in receiver_distances:
        p = np.zeros(nt)
        v = np.zeros(nt)
        for i, ti in enumerate(t):
            shock_r = physics.shock_radius(ti)
            if d < shock_r:
                p[i] = physics.peak_overpressure(d) * (1 - (shock_r - d) / shock_r)
            
            # Ground velocity
            vp, vs = 6000.0, 3500.0
            rp, rs = vp * ti, vs * ti
            v[i] = (np.exp(-((d - rp) / 500)**2) * 0.1 + 
                   np.exp(-((d - rs) / 300)**2) * 0.2) * np.exp(-d / 50000)
        
        receiver_data[d] = {'pressure': p / 1000, 'velocity': v}
    
    data['receivers'] = receiver_data
    
    # Fallout pattern (simplified)
    # Assume wind from west (positive x direction)
    wind_speed = 10.0  # m/s at cloud height
    cloud_height = 15000.0  # m
    settling_velocity = 1.0  # m/s for average particle
    
    fallout = np.zeros((ny, nx))
    for i in range(ny):
        for j in range(nx):
            dx = x[j] - x_burst
            dy = y[i] - y_burst
            
            # Downwind distance
            downwind = dx
            crosswind = np.abs(dy)
            
            if downwind > 0:
                # Time for fallout to arrive
                t_fall = cloud_height / settling_velocity
                # Drift distance
                drift = wind_speed * t_fall
                
                # Gaussian plume
                sigma_y = 0.1 * downwind + 1000
                sigma_x = 0.05 * downwind + 500
                
                fallout[i, j] = np.exp(-0.5 * ((downwind - drift) / sigma_x)**2) * \
                               np.exp(-0.5 * (crosswind / sigma_y)**2)
    
    # Normalize and scale (arbitrary units → Bq/m²)
    fallout = fallout / np.max(fallout + 1e-10) * 1e6
    data['fallout'] = fallout
    
    # Damage zones (overpressure thresholds) - vectorized
    P = physics.peak_overpressure(r.flatten()).reshape(r.shape)
    data['peak_overpressure'] = P / 1000  # kPa
    
    return data


def load_hdf5_data(directory):
    """Load simulation data from HDF5 files."""
    
    if not HAS_HDF5:
        return None
    
    data = {}
    
    # Look for HDF5 files
    h5_files = [f for f in os.listdir(directory) if f.endswith('.h5') or f.endswith('.hdf5')]
    
    if not h5_files:
        print(f"No HDF5 files found in {directory}")
        return None
    
    for h5_file in h5_files:
        filepath = os.path.join(directory, h5_file)
        try:
            with h5py.File(filepath, 'r') as f:
                print(f"Loading {h5_file}...")
                
                # Recursively load all datasets
                def load_group(group, prefix=''):
                    for key in group.keys():
                        item = group[key]
                        full_key = f"{prefix}/{key}" if prefix else key
                        
                        if isinstance(item, h5py.Dataset):
                            data[full_key] = item[:]
                        elif isinstance(item, h5py.Group):
                            load_group(item, full_key)
                
                load_group(f)
                
        except Exception as e:
            print(f"Error loading {h5_file}: {e}")
    
    return data if data else None


# =============================================================================
# Visualization Functions
# =============================================================================

class AtmosphericNuclearVisualizer:
    """Comprehensive visualizer for atmospheric nuclear test results."""
    
    def __init__(self, data, physics, output_dir='figures'):
        self.data = data
        self.physics = physics
        self.output_dir = output_dir
        
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        # Color maps
        self.cmap_pressure = 'hot'
        self.cmap_thermal = 'inferno'
        self.cmap_radiation = 'YlOrRd'
        self.cmap_velocity = 'seismic'
        self.cmap_fallout = 'YlGnBu'
        self.cmap_emp = 'plasma'
        
        # Damage thresholds (kPa)
        self.damage_levels = {
            'Total destruction': 140.0,
            'Severe damage': 35.0,
            'Moderate damage': 14.0,
            'Light damage': 3.5,
            'Glass breakage': 1.0
        }
        
        # Thermal thresholds (kJ/m²)
        self.thermal_levels = {
            'Third-degree burns': 670,
            'Second-degree burns': 335,
            'First-degree burns': 125,
            'No effect': 50
        }
    
    def plot_overview(self, time_index=-1, save=True):
        """Create an overview figure with multiple panels."""
        
        fig = plt.figure(figsize=(20, 16))
        fig.suptitle(f'Atmospheric Nuclear Test - {self.physics.yield_kt:.0f} kt at {self.physics.burst_height:.0f} m HOB',
                    fontsize=16, fontweight='bold')
        
        # Grid layout
        gs = fig.add_gridspec(3, 4, hspace=0.3, wspace=0.3)
        
        x = self.data['x'] / 1000  # Convert to km
        y = self.data['y'] / 1000
        X, Y = np.meshgrid(x, y)
        
        # 1. Peak overpressure
        ax1 = fig.add_subplot(gs[0, 0])
        P = self.data['peak_overpressure']
        im1 = ax1.contourf(X, Y, P, levels=np.logspace(-1, 3, 20), 
                          cmap=self.cmap_pressure, norm=mcolors.LogNorm())
        ax1.set_title('Peak Overpressure (kPa)')
        ax1.set_xlabel('X (km)')
        ax1.set_ylabel('Y (km)')
        plt.colorbar(im1, ax=ax1, label='kPa')
        
        # Add damage contours
        for label, level in self.damage_levels.items():
            ax1.contour(X, Y, P, levels=[level], colors='white', linewidths=1)
        
        # Mark burst point
        ax1.plot(50, 50, 'w*', markersize=15, markeredgecolor='k')
        
        # 2. Thermal fluence
        ax2 = fig.add_subplot(gs[0, 1])
        Q = self.data['thermal_fluence']
        im2 = ax2.contourf(X, Y, Q, levels=np.logspace(0, 4, 20),
                          cmap=self.cmap_thermal, norm=mcolors.LogNorm())
        ax2.set_title('Thermal Fluence (kJ/m²)')
        ax2.set_xlabel('X (km)')
        ax2.set_ylabel('Y (km)')
        plt.colorbar(im2, ax=ax2, label='kJ/m²')
        ax2.plot(50, 50, 'w*', markersize=15, markeredgecolor='k')
        
        # 3. Radiation dose
        ax3 = fig.add_subplot(gs[0, 2])
        D = self.data['radiation_dose']
        D_clipped = np.clip(D, 1e-3, 1e3)
        im3 = ax3.contourf(X, Y, D_clipped, levels=np.logspace(-3, 3, 20),
                          cmap=self.cmap_radiation, norm=mcolors.LogNorm())
        ax3.set_title('Prompt Radiation Dose (Gy)')
        ax3.set_xlabel('X (km)')
        ax3.set_ylabel('Y (km)')
        plt.colorbar(im3, ax=ax3, label='Gy')
        ax3.plot(50, 50, 'w*', markersize=15, markeredgecolor='k')
        
        # 4. EMP spatial distribution
        ax4 = fig.add_subplot(gs[0, 3])
        E = self.data['emp_spatial']
        im4 = ax4.contourf(X, Y, E, levels=20, cmap=self.cmap_emp)
        ax4.set_title('E1 EMP Field (V/m)')
        ax4.set_xlabel('X (km)')
        ax4.set_ylabel('Y (km)')
        plt.colorbar(im4, ax=ax4, label='V/m')
        ax4.plot(50, 50, 'w*', markersize=15, markeredgecolor='k')
        
        # 5. Overpressure time evolution
        ax5 = fig.add_subplot(gs[1, 0:2])
        t = self.data['time']
        for d, rd in self.data['receivers'].items():
            ax5.plot(t, rd['pressure'], label=f'{d/1000:.0f} km')
        ax5.set_xlabel('Time (s)')
        ax5.set_ylabel('Overpressure (kPa)')
        ax5.set_title('Pressure Time History at Various Ranges')
        ax5.legend()
        ax5.set_xlim([0, max(t)])
        ax5.grid(True, alpha=0.3)
        
        # 6. Ground velocity
        ax6 = fig.add_subplot(gs[1, 2:4])
        for d, rd in self.data['receivers'].items():
            ax6.plot(t, rd['velocity'], label=f'{d/1000:.0f} km')
        ax6.set_xlabel('Time (s)')
        ax6.set_ylabel('Ground Velocity (m/s)')
        ax6.set_title('Seismic Ground Motion')
        ax6.legend()
        ax6.grid(True, alpha=0.3)
        
        # 7. EMP time history
        ax7 = fig.add_subplot(gs[2, 0])
        t_emp = np.linspace(0, 0.001, 1000)  # 1 ms for E1
        e1 = [self.physics.emp_e1_field(10000, ti) for ti in t_emp]
        ax7.plot(t_emp * 1e6, e1)  # Convert to microseconds
        ax7.set_xlabel('Time (µs)')
        ax7.set_ylabel('Electric Field (V/m)')
        ax7.set_title('E1 EMP at 10 km')
        ax7.grid(True, alpha=0.3)
        
        # 8. E3 time history
        ax8 = fig.add_subplot(gs[2, 1])
        t_e3 = np.linspace(0, 200, 1000)
        e3 = [self.physics.emp_e3_field(10000, ti) for ti in t_e3]
        ax8.plot(t_e3, e3)
        ax8.set_xlabel('Time (s)')
        ax8.set_ylabel('Electric Field (V/m)')
        ax8.set_title('E3 MHD-EMP at 10 km')
        ax8.grid(True, alpha=0.3)
        
        # 9. Fallout pattern
        ax9 = fig.add_subplot(gs[2, 2:4])
        fallout = self.data['fallout']
        fallout_plot = np.maximum(fallout, 1)
        im9 = ax9.contourf(X, Y, fallout_plot, levels=np.logspace(0, 6, 20),
                          cmap=self.cmap_fallout, norm=mcolors.LogNorm())
        ax9.set_title('Fallout Deposition (Bq/m²)')
        ax9.set_xlabel('X (km)')
        ax9.set_ylabel('Y (km)')
        plt.colorbar(im9, ax=ax9, label='Bq/m²')
        ax9.plot(50, 50, 'r*', markersize=15)
        
        # Add wind direction arrow
        ax9.annotate('', xy=(70, 50), xytext=(60, 50),
                    arrowprops=dict(arrowstyle='->', color='red', lw=2))
        ax9.text(65, 52, 'Wind', color='red', fontsize=10)
        
        plt.tight_layout()
        
        if save:
            plt.savefig(os.path.join(self.output_dir, 'overview.png'), dpi=150, bbox_inches='tight')
            print(f"Saved: {self.output_dir}/overview.png")
        
        return fig
    
    def plot_damage_zones(self, save=True):
        """Plot damage zone map with legends."""
        
        fig, axes = plt.subplots(1, 3, figsize=(18, 6))
        fig.suptitle('Nuclear Effects Damage Zones', fontsize=14, fontweight='bold')
        
        x = self.data['x'] / 1000
        y = self.data['y'] / 1000
        X, Y = np.meshgrid(x, y)
        
        # 1. Blast damage zones
        ax1 = axes[0]
        P = self.data['peak_overpressure']
        
        # Create categorical damage map
        damage_zones = np.zeros_like(P)
        damage_zones[P >= 1.0] = 1    # Glass breakage
        damage_zones[P >= 3.5] = 2    # Light damage
        damage_zones[P >= 14.0] = 3   # Moderate damage
        damage_zones[P >= 35.0] = 4   # Severe damage
        damage_zones[P >= 140.0] = 5  # Total destruction
        
        cmap_damage = mcolors.ListedColormap(['white', 'lightblue', 'yellow', 'orange', 'red', 'darkred'])
        bounds = [0, 1, 2, 3, 4, 5, 6]
        norm = mcolors.BoundaryNorm(bounds, cmap_damage.N)
        
        im1 = ax1.imshow(damage_zones, extent=[x.min(), x.max(), y.min(), y.max()],
                        origin='lower', cmap=cmap_damage, norm=norm)
        ax1.set_title('Blast Damage Zones')
        ax1.set_xlabel('X (km)')
        ax1.set_ylabel('Y (km)')
        ax1.plot(50, 50, 'w*', markersize=15, markeredgecolor='k')
        
        # Legend
        labels = ['No damage', 'Glass breakage\n(>1 kPa)', 'Light damage\n(>3.5 kPa)',
                 'Moderate\n(>14 kPa)', 'Severe\n(>35 kPa)', 'Total\n(>140 kPa)']
        cbar1 = plt.colorbar(im1, ax=ax1, ticks=[0.5, 1.5, 2.5, 3.5, 4.5, 5.5])
        cbar1.set_ticklabels(labels)
        
        # 2. Thermal damage zones
        ax2 = axes[1]
        Q = self.data['thermal_fluence']
        
        thermal_zones = np.zeros_like(Q)
        thermal_zones[Q >= 50] = 1    # No injury
        thermal_zones[Q >= 125] = 2   # First-degree burns
        thermal_zones[Q >= 335] = 3   # Second-degree burns
        thermal_zones[Q >= 670] = 4   # Third-degree burns
        
        cmap_thermal = mcolors.ListedColormap(['white', 'lightyellow', 'yellow', 'orange', 'red'])
        bounds_t = [0, 1, 2, 3, 4, 5]
        norm_t = mcolors.BoundaryNorm(bounds_t, cmap_thermal.N)
        
        im2 = ax2.imshow(thermal_zones, extent=[x.min(), x.max(), y.min(), y.max()],
                        origin='lower', cmap=cmap_thermal, norm=norm_t)
        ax2.set_title('Thermal Burn Zones')
        ax2.set_xlabel('X (km)')
        ax2.set_ylabel('Y (km)')
        ax2.plot(50, 50, 'w*', markersize=15, markeredgecolor='k')
        
        labels_t = ['No effect', 'Mild\n(>50 kJ/m²)', '1st degree\n(>125 kJ/m²)',
                   '2nd degree\n(>335 kJ/m²)', '3rd degree\n(>670 kJ/m²)']
        cbar2 = plt.colorbar(im2, ax=ax2, ticks=[0.5, 1.5, 2.5, 3.5, 4.5])
        cbar2.set_ticklabels(labels_t)
        
        # 3. Radiation zones
        ax3 = axes[2]
        D = self.data['radiation_dose']
        
        rad_zones = np.zeros_like(D)
        rad_zones[D >= 0.01] = 1   # Detectable
        rad_zones[D >= 0.5] = 2    # Slight symptoms
        rad_zones[D >= 2.0] = 3    # Radiation sickness
        rad_zones[D >= 6.0] = 4    # LD50
        rad_zones[D >= 10.0] = 5   # Fatal
        
        cmap_rad = mcolors.ListedColormap(['white', 'lightgreen', 'yellow', 'orange', 'red', 'purple'])
        bounds_r = [0, 1, 2, 3, 4, 5, 6]
        norm_r = mcolors.BoundaryNorm(bounds_r, cmap_rad.N)
        
        im3 = ax3.imshow(rad_zones, extent=[x.min(), x.max(), y.min(), y.max()],
                        origin='lower', cmap=cmap_rad, norm=norm_r)
        ax3.set_title('Radiation Dose Zones')
        ax3.set_xlabel('X (km)')
        ax3.set_ylabel('Y (km)')
        ax3.plot(50, 50, 'w*', markersize=15, markeredgecolor='k')
        
        labels_r = ['Negligible', 'Detectable\n(>10 mGy)', 'Mild\n(>0.5 Gy)',
                   'Sickness\n(>2 Gy)', 'LD50\n(>6 Gy)', 'Fatal\n(>10 Gy)']
        cbar3 = plt.colorbar(im3, ax=ax3, ticks=[0.5, 1.5, 2.5, 3.5, 4.5, 5.5])
        cbar3.set_ticklabels(labels_r)
        
        plt.tight_layout()
        
        if save:
            plt.savefig(os.path.join(self.output_dir, 'damage_zones.png'), dpi=150, bbox_inches='tight')
            print(f"Saved: {self.output_dir}/damage_zones.png")
        
        return fig
    
    def plot_blast_wave_evolution(self, save=True):
        """Animate blast wave propagation."""
        
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        fig.suptitle('Blast Wave Evolution', fontsize=14, fontweight='bold')
        
        x = self.data['x'] / 1000
        y = self.data['y'] / 1000
        X, Y = np.meshgrid(x, y)
        
        t = self.data['time']
        P = self.data['overpressure']
        
        # Select time snapshots
        time_indices = [0, len(t)//10, len(t)//5, len(t)//3, len(t)//2, len(t)-1]
        
        for idx, (ax, ti) in enumerate(zip(axes.flat, time_indices)):
            im = ax.contourf(X, Y, P[ti], levels=np.linspace(0, np.percentile(P, 99), 30),
                           cmap=self.cmap_pressure)
            ax.set_title(f't = {t[ti]:.2f} s')
            ax.set_xlabel('X (km)')
            ax.set_ylabel('Y (km)')
            ax.plot(50, 50, 'w*', markersize=10, markeredgecolor='k')
            
            # Add shock front circle
            shock_r = self.physics.shock_radius(t[ti]) / 1000
            circle = Circle((50, 50), shock_r, fill=False, color='white', linestyle='--', linewidth=2)
            ax.add_patch(circle)
            
            plt.colorbar(im, ax=ax, label='Overpressure (kPa)')
        
        plt.tight_layout()
        
        if save:
            plt.savefig(os.path.join(self.output_dir, 'blast_evolution.png'), dpi=150, bbox_inches='tight')
            print(f"Saved: {self.output_dir}/blast_evolution.png")
        
        return fig
    
    def plot_emp_analysis(self, save=True):
        """Detailed EMP analysis plots."""
        
        fig = plt.figure(figsize=(16, 12))
        fig.suptitle('Electromagnetic Pulse (EMP) Analysis', fontsize=14, fontweight='bold')
        
        gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)
        
        # 1. E1 waveform (nanoseconds)
        ax1 = fig.add_subplot(gs[0, 0])
        t_ns = np.linspace(0, 100, 1000)  # nanoseconds
        e1_waveform = [self.physics.emp_e1_field(10000, t * 1e-9) for t in t_ns]
        ax1.plot(t_ns, e1_waveform)
        ax1.set_xlabel('Time (ns)')
        ax1.set_ylabel('Electric Field (V/m)')
        ax1.set_title('E1 Component (10 km)')
        ax1.grid(True, alpha=0.3)
        ax1.axhline(y=0, color='k', linestyle='-', linewidth=0.5)
        
        # 2. E1 waveform (microseconds)
        ax2 = fig.add_subplot(gs[0, 1])
        t_us = np.linspace(0, 10, 1000)  # microseconds
        e1_us = [self.physics.emp_e1_field(10000, t * 1e-6) for t in t_us]
        ax2.plot(t_us, e1_us)
        ax2.set_xlabel('Time (µs)')
        ax2.set_ylabel('Electric Field (V/m)')
        ax2.set_title('E1 Decay (10 km)')
        ax2.grid(True, alpha=0.3)
        
        # 3. E1 vs distance
        ax3 = fig.add_subplot(gs[0, 2])
        distances = np.linspace(1000, 100000, 100)
        e1_peak = [self.physics.emp_e1_field(d, 1e-8) for d in distances]
        ax3.semilogy(distances / 1000, e1_peak)
        ax3.set_xlabel('Distance (km)')
        ax3.set_ylabel('Peak E1 Field (V/m)')
        ax3.set_title('E1 Peak vs Distance')
        ax3.grid(True, alpha=0.3)
        
        # 4. E3 waveform
        ax4 = fig.add_subplot(gs[1, 0])
        t_s = np.linspace(0, 200, 1000)
        e3_waveform = [self.physics.emp_e3_field(10000, t) for t in t_s]
        ax4.plot(t_s, e3_waveform)
        ax4.set_xlabel('Time (s)')
        ax4.set_ylabel('Electric Field (V/m)')
        ax4.set_title('E3 MHD-EMP (10 km)')
        ax4.grid(True, alpha=0.3)
        
        # 5. E3 vs distance
        ax5 = fig.add_subplot(gs[1, 1])
        e3_peak = [self.physics.emp_e3_field(d, 50) for d in distances]
        ax5.plot(distances / 1000, e3_peak)
        ax5.set_xlabel('Distance (km)')
        ax5.set_ylabel('Peak E3 Field (V/m)')
        ax5.set_title('E3 Peak vs Distance')
        ax5.grid(True, alpha=0.3)
        
        # 6. Combined spectrum (frequency domain)
        ax6 = fig.add_subplot(gs[1, 2])
        # Approximate frequency content
        f = np.logspace(3, 11, 1000)  # Hz
        
        # E1: high frequency (MHz to GHz)
        e1_spec = 1e4 * np.exp(-(np.log10(f) - 8)**2 / 2)
        # E2: intermediate (kHz to MHz)
        e2_spec = 1e3 * np.exp(-(np.log10(f) - 5)**2 / 2)
        # E3: low frequency (mHz to Hz)
        e3_spec = 1e2 * np.exp(-(np.log10(f) - 0)**2 / 2)
        
        ax6.loglog(f, e1_spec, label='E1', color='red')
        ax6.loglog(f, e2_spec, label='E2', color='orange')
        ax6.loglog(f, e3_spec, label='E3', color='blue')
        ax6.set_xlabel('Frequency (Hz)')
        ax6.set_ylabel('Spectral Density (V/m/Hz)')
        ax6.set_title('EMP Frequency Spectrum')
        ax6.legend()
        ax6.grid(True, alpha=0.3)
        
        # 7. Induced voltage in power line
        ax7 = fig.add_subplot(gs[2, 0])
        line_length = 10000  # 10 km line
        t_line = np.linspace(0, 10e-6, 1000)
        
        # Simplified coupling (V = E * L * coupling_factor)
        coupling = 0.3
        v_induced = [self.physics.emp_e1_field(10000, t) * line_length * coupling for t in t_line]
        
        ax7.plot(t_line * 1e6, np.array(v_induced) / 1000)  # kV
        ax7.set_xlabel('Time (µs)')
        ax7.set_ylabel('Induced Voltage (kV)')
        ax7.set_title('Induced Voltage on 10 km Power Line')
        ax7.grid(True, alpha=0.3)
        
        # 8. EMP spatial distribution
        ax8 = fig.add_subplot(gs[2, 1:3])
        x = self.data['x'] / 1000
        y = self.data['y'] / 1000
        X, Y = np.meshgrid(x, y)
        E = self.data['emp_spatial']
        
        im8 = ax8.contourf(X, Y, E, levels=30, cmap=self.cmap_emp)
        ax8.set_title('E1 Peak Field Distribution (V/m)')
        ax8.set_xlabel('X (km)')
        ax8.set_ylabel('Y (km)')
        ax8.plot(50, 50, 'w*', markersize=15, markeredgecolor='k')
        plt.colorbar(im8, ax=ax8, label='V/m')
        
        # Add contours
        ax8.contour(X, Y, E, levels=[1000, 5000, 10000, 25000], colors='white', linewidths=1)
        
        plt.tight_layout()
        
        if save:
            plt.savefig(os.path.join(self.output_dir, 'emp_analysis.png'), dpi=150, bbox_inches='tight')
            print(f"Saved: {self.output_dir}/emp_analysis.png")
        
        return fig
    
    def plot_seismic_analysis(self, save=True):
        """Seismic wave analysis from ground coupling."""
        
        fig, axes = plt.subplots(2, 3, figsize=(18, 10))
        fig.suptitle('Seismic Ground Motion from Air-Ground Coupling', fontsize=14, fontweight='bold')
        
        t = self.data['time']
        
        # 1. Seismograms at different distances
        ax1 = axes[0, 0]
        for d, rd in sorted(self.data['receivers'].items()):
            offset = d / 10000  # Offset for visibility
            ax1.plot(t, rd['velocity'] + offset, label=f'{d/1000:.0f} km')
        ax1.set_xlabel('Time (s)')
        ax1.set_ylabel('Velocity + offset (m/s)')
        ax1.set_title('Ground Velocity Seismograms')
        ax1.legend(loc='upper right')
        ax1.grid(True, alpha=0.3)
        
        # 2. Record section
        ax2 = axes[0, 1]
        distances = sorted(self.data['receivers'].keys())
        velocities = np.array([self.data['receivers'][d]['velocity'] for d in distances])
        
        extent = [t.min(), t.max(), distances[0]/1000, distances[-1]/1000]
        im2 = ax2.imshow(velocities, aspect='auto', origin='lower', extent=extent,
                        cmap='seismic', vmin=-0.5, vmax=0.5)
        ax2.set_xlabel('Time (s)')
        ax2.set_ylabel('Distance (km)')
        ax2.set_title('Velocity Record Section')
        plt.colorbar(im2, ax=ax2, label='m/s')
        
        # Add moveout lines
        vp, vs = 6.0, 3.5  # km/s
        for v, label, color in [(vp, 'P-wave', 'blue'), (vs, 'S-wave', 'green')]:
            t_moveout = np.array(distances) / 1000 / v
            ax2.plot(t_moveout, np.array(distances)/1000, '--', color=color, label=label)
        ax2.legend(loc='lower right')
        
        # 3. Peak velocity vs distance
        ax3 = axes[0, 2]
        peak_v = [np.max(np.abs(self.data['receivers'][d]['velocity'])) for d in distances]
        ax3.loglog(np.array(distances)/1000, peak_v, 'o-')
        ax3.set_xlabel('Distance (km)')
        ax3.set_ylabel('Peak Ground Velocity (m/s)')
        ax3.set_title('PGV vs Distance')
        ax3.grid(True, alpha=0.3)
        
        # Fit power law
        from scipy.optimize import curve_fit
        def power_law(r, a, b):
            return a * r**b
        try:
            popt, _ = curve_fit(power_law, np.array(distances)/1000, peak_v)
            r_fit = np.linspace(1, 50, 100)
            ax3.loglog(r_fit, power_law(r_fit, *popt), 'r--', 
                      label=f'Fit: {popt[0]:.2f}×r^{popt[1]:.2f}')
            ax3.legend()
        except:
            pass
        
        # 4. Ground velocity field snapshot
        ax4 = axes[1, 0]
        x = self.data['x'] / 1000
        y = self.data['y'] / 1000
        X, Y = np.meshgrid(x, y)
        
        # Find time with maximum activity
        gv = self.data['ground_velocity']
        energy_vs_time = np.sum(gv**2, axis=(1, 2))
        t_max = np.argmax(energy_vs_time)
        
        im4 = ax4.contourf(X, Y, gv[t_max], levels=30, cmap='seismic', 
                          vmin=-np.percentile(np.abs(gv), 99),
                          vmax=np.percentile(np.abs(gv), 99))
        ax4.set_title(f'Ground Velocity at t = {t[t_max]:.2f} s')
        ax4.set_xlabel('X (km)')
        ax4.set_ylabel('Y (km)')
        ax4.plot(50, 50, 'w*', markersize=10, markeredgecolor='k')
        plt.colorbar(im4, ax=ax4, label='m/s')
        
        # 5. Spectrum at one station
        ax5 = axes[1, 1]
        d_ref = 10000  # 10 km
        v_ref = self.data['receivers'][d_ref]['velocity']
        
        # Compute FFT
        dt = t[1] - t[0]
        n = len(v_ref)
        freq = np.fft.rfftfreq(n, dt)
        spectrum = np.abs(np.fft.rfft(v_ref))
        
        ax5.loglog(freq[1:], spectrum[1:])
        ax5.set_xlabel('Frequency (Hz)')
        ax5.set_ylabel('Spectral Amplitude')
        ax5.set_title(f'Velocity Spectrum at {d_ref/1000:.0f} km')
        ax5.grid(True, alpha=0.3)
        ax5.set_xlim([0.01, 10])
        
        # 6. Pressure vs ground motion correlation
        ax6 = axes[1, 2]
        for d in [5000, 10000, 20000]:
            if d in self.data['receivers']:
                p = self.data['receivers'][d]['pressure']
                v = self.data['receivers'][d]['velocity']
                ax6.plot(p, v, '.', alpha=0.5, label=f'{d/1000:.0f} km', markersize=2)
        
        ax6.set_xlabel('Overpressure (kPa)')
        ax6.set_ylabel('Ground Velocity (m/s)')
        ax6.set_title('Air-Ground Coupling')
        ax6.legend()
        ax6.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if save:
            plt.savefig(os.path.join(self.output_dir, 'seismic_analysis.png'), dpi=150, bbox_inches='tight')
            print(f"Saved: {self.output_dir}/seismic_analysis.png")
        
        return fig
    
    def plot_radiation_analysis(self, save=True):
        """Detailed radiation analysis."""
        
        fig = plt.figure(figsize=(16, 12))
        fig.suptitle('Radiation Effects Analysis', fontsize=14, fontweight='bold')
        
        gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)
        
        x = self.data['x'] / 1000
        y = self.data['y'] / 1000
        X, Y = np.meshgrid(x, y)
        
        # 1. Prompt gamma dose
        ax1 = fig.add_subplot(gs[0, 0])
        D = self.data['radiation_dose']
        D_clipped = np.clip(D, 1e-4, 1e3)
        im1 = ax1.contourf(X, Y, D_clipped, levels=np.logspace(-4, 3, 20),
                          cmap=self.cmap_radiation, norm=mcolors.LogNorm())
        ax1.set_title('Prompt Gamma Dose (Gy)')
        ax1.set_xlabel('X (km)')
        ax1.set_ylabel('Y (km)')
        ax1.plot(50, 50, 'w*', markersize=12, markeredgecolor='k')
        plt.colorbar(im1, ax=ax1, label='Gy')
        
        # 2. Dose vs distance
        ax2 = fig.add_subplot(gs[0, 1])
        r = np.linspace(100, 50000, 200)
        dose_gamma = [self.physics.prompt_radiation_dose(ri) for ri in r]
        ax2.semilogy(r/1000, dose_gamma, label='Prompt gamma')
        
        # Add thresholds
        thresholds = {'LD50/60': 4.5, 'Sickness': 1.0, 'Detectable': 0.1}
        for label, val in thresholds.items():
            ax2.axhline(y=val, linestyle='--', alpha=0.7, label=label)
        
        ax2.set_xlabel('Distance (km)')
        ax2.set_ylabel('Dose (Gy)')
        ax2.set_title('Radiation Dose vs Distance')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        ax2.set_xlim([0, 50])
        ax2.set_ylim([1e-3, 1e3])
        
        # 3. Fallout pattern
        ax3 = fig.add_subplot(gs[0, 2])
        fallout = self.data['fallout']
        fallout_plot = np.maximum(fallout, 1)
        im3 = ax3.contourf(X, Y, fallout_plot, levels=np.logspace(0, 6, 20),
                          cmap=self.cmap_fallout, norm=mcolors.LogNorm())
        ax3.set_title('Ground Contamination (Bq/m²)')
        ax3.set_xlabel('X (km)')
        ax3.set_ylabel('Y (km)')
        ax3.plot(50, 50, 'r*', markersize=12)
        plt.colorbar(im3, ax=ax3, label='Bq/m²')
        
        # Wind arrow
        ax3.annotate('', xy=(75, 50), xytext=(60, 50),
                    arrowprops=dict(arrowstyle='->', color='red', lw=2))
        
        # 4. Fallout contours with dose rate
        ax4 = fig.add_subplot(gs[1, 0])
        
        # Convert Bq/m² to dose rate (approximate)
        # 1 MBq/m² ≈ 0.1 mSv/h for Cs-137
        dose_rate = fallout * 1e-7  # mSv/h
        dose_rate_plot = np.maximum(dose_rate, 1e-6)
        
        im4 = ax4.contourf(X, Y, dose_rate_plot, levels=np.logspace(-6, 2, 20),
                          cmap='YlOrRd', norm=mcolors.LogNorm())
        ax4.set_title('Fallout Dose Rate (mSv/h)')
        ax4.set_xlabel('X (km)')
        ax4.set_ylabel('Y (km)')
        ax4.plot(50, 50, 'k*', markersize=12)
        plt.colorbar(im4, ax=ax4, label='mSv/h')
        
        # Add significant contours
        ax4.contour(X, Y, dose_rate, levels=[0.001, 0.01, 0.1, 1, 10], 
                   colors='black', linewidths=0.5)
        
        # 5. Fallout arrival time
        ax5 = fig.add_subplot(gs[1, 1])
        
        # Estimate arrival time based on distance and wind
        wind_speed = 10.0  # m/s
        cloud_height = 15000.0
        settling = 1.0
        
        arrival_time = np.zeros_like(X)
        for i in range(len(y)):
            for j in range(len(x)):
                dx = x[j] - 50  # km from GZ
                if dx > 0:
                    # Time = settling time + drift time
                    t_settle = cloud_height / settling / 3600  # hours
                    t_drift = (dx * 1000) / wind_speed / 3600  # hours
                    arrival_time[i, j] = t_settle + t_drift
                else:
                    arrival_time[i, j] = np.nan
        
        im5 = ax5.contourf(X, Y, arrival_time, levels=np.linspace(0, 24, 25), cmap='viridis')
        ax5.set_title('Fallout Arrival Time (hours)')
        ax5.set_xlabel('X (km)')
        ax5.set_ylabel('Y (km)')
        ax5.plot(50, 50, 'r*', markersize=12)
        plt.colorbar(im5, ax=ax5, label='hours')
        
        # 6. Total dose integration
        ax6 = fig.add_subplot(gs[1, 2])
        
        # Cumulative dose over time (simplified)
        times = np.array([1, 4, 24, 48, 168, 720])  # hours
        labels = ['1h', '4h', '1d', '2d', '1w', '1mo']
        
        # Way-Wigner decay: dose rate ∝ t^-1.2
        # Integrated dose up to time t: ∝ t^-0.2
        
        r_sample = [5, 10, 20, 30]  # km downwind
        
        for ri in r_sample:
            # Peak dose rate at this location
            idx_x = int((50 + ri) / 100 * len(x))
            idx_y = len(y) // 2
            if idx_x < len(x):
                dr0 = dose_rate[idx_y, idx_x]
                cumulative = dr0 * times**0.8 / 0.8  # Integrate t^-1.2
                ax6.loglog(times, cumulative, 'o-', label=f'{ri} km')
        
        ax6.axhline(y=1000, linestyle='--', color='red', alpha=0.7, label='Evacuation (1 Sv)')
        ax6.axhline(y=100, linestyle='--', color='orange', alpha=0.7, label='Shelter (100 mSv)')
        ax6.set_xlabel('Time after detonation (hours)')
        ax6.set_ylabel('Cumulative Dose (mSv)')
        ax6.set_title('Fallout Dose Accumulation')
        ax6.legend()
        ax6.grid(True, alpha=0.3)
        
        # 7. Radiation shielding effectiveness
        ax7 = fig.add_subplot(gs[2, 0])
        
        shielding = np.linspace(0, 100, 100)  # cm of concrete
        # Half-value layer for gamma: ~6 cm concrete
        hvl = 6.0
        transmission = 0.5 ** (shielding / hvl)
        
        ax7.semilogy(shielding, transmission * 100)
        ax7.set_xlabel('Concrete Thickness (cm)')
        ax7.set_ylabel('Transmission (%)')
        ax7.set_title('Gamma Shielding Effectiveness')
        ax7.grid(True, alpha=0.3)
        ax7.axhline(y=10, linestyle='--', color='green', alpha=0.7)
        ax7.axhline(y=1, linestyle='--', color='orange', alpha=0.7)
        
        # 8. Combined hazard zones
        ax8 = fig.add_subplot(gs[2, 1:3])
        
        # Create combined hazard index
        # Normalize each effect to 0-1 scale based on lethality
        P = self.data['peak_overpressure']
        Q = self.data['thermal_fluence']
        D = self.data['radiation_dose']
        
        hazard_blast = np.clip(P / 35.0, 0, 1)  # 35 kPa = severe
        hazard_thermal = np.clip(Q / 670.0, 0, 1)  # 670 kJ/m² = 3rd degree
        hazard_rad = np.clip(D / 6.0, 0, 1)  # 6 Gy = LD50
        
        # Combined (maximum of effects)
        hazard_combined = np.maximum(np.maximum(hazard_blast, hazard_thermal), hazard_rad)
        
        im8 = ax8.contourf(X, Y, hazard_combined, levels=np.linspace(0, 1, 20),
                          cmap='RdYlGn_r')
        ax8.set_title('Combined Hazard Index (0=safe, 1=lethal)')
        ax8.set_xlabel('X (km)')
        ax8.set_ylabel('Y (km)')
        ax8.plot(50, 50, 'k*', markersize=15)
        plt.colorbar(im8, ax=ax8, label='Hazard Index')
        
        # Mark dominant effect zones
        blast_dom = (hazard_blast > hazard_thermal) & (hazard_blast > hazard_rad)
        thermal_dom = (hazard_thermal > hazard_blast) & (hazard_thermal > hazard_rad)
        rad_dom = (hazard_rad > hazard_blast) & (hazard_rad > hazard_thermal)
        
        ax8.contour(X, Y, blast_dom.astype(float), levels=[0.5], colors='blue', linewidths=2)
        ax8.contour(X, Y, thermal_dom.astype(float), levels=[0.5], colors='red', linewidths=2)
        ax8.contour(X, Y, rad_dom.astype(float), levels=[0.5], colors='green', linewidths=2)
        
        # Legend
        ax8.plot([], [], 'b-', linewidth=2, label='Blast dominant')
        ax8.plot([], [], 'r-', linewidth=2, label='Thermal dominant')
        ax8.plot([], [], 'g-', linewidth=2, label='Radiation dominant')
        ax8.legend(loc='upper right')
        
        plt.tight_layout()
        
        if save:
            plt.savefig(os.path.join(self.output_dir, 'radiation_analysis.png'), dpi=150, bbox_inches='tight')
            print(f"Saved: {self.output_dir}/radiation_analysis.png")
        
        return fig
    
    def plot_fireball_evolution(self, save=True):
        """Fireball and mushroom cloud evolution."""
        
        fig = plt.figure(figsize=(16, 10))
        fig.suptitle('Fireball and Mushroom Cloud Evolution', fontsize=14, fontweight='bold')
        
        gs = fig.add_gridspec(2, 3, hspace=0.3, wspace=0.3)
        
        # 1. Fireball radius vs time
        ax1 = fig.add_subplot(gs[0, 0])
        t = np.logspace(-6, 2, 1000)  # 1 µs to 100 s
        r_fb = [self.physics.fireball_radius(ti) for ti in t]
        r_shock = [self.physics.shock_radius(ti) for ti in t]
        
        ax1.loglog(t, r_fb, label='Fireball radius')
        ax1.loglog(t, r_shock, '--', label='Shock front')
        ax1.set_xlabel('Time (s)')
        ax1.set_ylabel('Radius (m)')
        ax1.set_title('Fireball & Shock Expansion')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        ax1.set_xlim([1e-6, 100])
        
        # 2. Fireball temperature vs time
        ax2 = fig.add_subplot(gs[0, 1])
        
        # Temperature evolution (simplified)
        t_temp = np.logspace(-6, 1, 1000)
        T = np.zeros_like(t_temp)
        for i, ti in enumerate(t_temp):
            if ti < 1e-6:
                T[i] = 1e8
            elif ti < 0.01:
                T[i] = 1e8 * (1e-6 / ti)**0.5
            elif ti < 0.1:
                T[i] = 8000
            else:
                T[i] = 8000 * (0.1 / ti)**0.5
        
        ax2.loglog(t_temp, T)
        ax2.axhline(y=5778, linestyle='--', color='orange', alpha=0.7, label='Solar surface')
        ax2.set_xlabel('Time (s)')
        ax2.set_ylabel('Temperature (K)')
        ax2.set_title('Fireball Surface Temperature')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        # 3. Thermal power vs time
        ax3 = fig.add_subplot(gs[0, 2])
        
        # P = 4πR² σT⁴
        sigma = 5.67e-8
        P_thermal = 4 * np.pi * np.array(r_fb)**2 * sigma * np.array([
            (1e8 if ti < 1e-6 else 8000 * max(1, (0.1/ti)**0.5))**4 for ti in t
        ])
        
        ax3.loglog(t, P_thermal / 1e15)  # PW
        ax3.set_xlabel('Time (s)')
        ax3.set_ylabel('Thermal Power (PW)')
        ax3.set_title('Thermal Radiation Power')
        ax3.grid(True, alpha=0.3)
        ax3.set_xlim([1e-6, 100])
        
        # 4. Cloud rise (2D visualization)
        ax4 = fig.add_subplot(gs[1, 0:2])
        
        # Draw mushroom cloud at different times
        times_cloud = [1, 5, 10, 30, 60, 120]  # seconds
        colors = plt.cm.Reds(np.linspace(0.3, 1, len(times_cloud)))
        
        for ti, c in zip(times_cloud, colors):
            # Cloud parameters (empirical)
            cloud_top = min(15000, 1000 * ti**0.5)  # m
            cloud_radius = min(5000, 200 * ti**0.5)
            stem_width = cloud_radius * 0.3
            
            # Draw stem
            stem = plt.Rectangle((-stem_width/2/1000, 0), stem_width/1000, cloud_top/1000 * 0.7,
                                 facecolor=c, alpha=0.3, edgecolor=c)
            ax4.add_patch(stem)
            
            # Draw cap (ellipse)
            cap = plt.matplotlib.patches.Ellipse(
                (0, cloud_top/1000), cloud_radius*2/1000, cloud_radius/1000,
                facecolor=c, alpha=0.5, edgecolor=c, label=f't={ti}s')
            ax4.add_patch(cap)
        
        ax4.set_xlim([-20, 20])
        ax4.set_ylim([0, 20])
        ax4.set_xlabel('Horizontal Distance (km)')
        ax4.set_ylabel('Altitude (km)')
        ax4.set_title('Mushroom Cloud Rise')
        ax4.legend(loc='upper right')
        ax4.set_aspect('equal')
        ax4.grid(True, alpha=0.3)
        
        # Add tropopause and stratosphere labels
        ax4.axhline(y=11, linestyle='--', color='blue', alpha=0.5)
        ax4.text(15, 11.5, 'Tropopause', fontsize=9)
        
        # 5. Optical flash time history
        ax5 = fig.add_subplot(gs[1, 2])
        
        t_flash = np.linspace(0, 10, 1000)  # seconds
        
        # Double-flash phenomenon
        flash = np.zeros_like(t_flash)
        for i, ti in enumerate(t_flash):
            if ti < 0.001:
                # First flash (X-ray)
                flash[i] = np.exp(-ti / 0.0001)
            elif ti < 0.01:
                # First minimum
                flash[i] = 0.1 * np.exp(-(ti - 0.01)**2 / 0.0001)
            elif ti < 1:
                # Second maximum (shock breakout)
                flash[i] = 0.8 * np.exp(-(ti - 0.5)**2 / 0.1)
            else:
                # Decay
                flash[i] = 0.5 * np.exp(-(ti - 1) / 2)
        
        ax5.plot(t_flash * 1000, flash)  # ms
        ax5.set_xlabel('Time (ms)')
        ax5.set_ylabel('Relative Brightness')
        ax5.set_title('Optical Flash (Double Flash)')
        ax5.grid(True, alpha=0.3)
        ax5.axvline(x=10, linestyle='--', color='red', alpha=0.5)
        ax5.text(12, 0.9, 'First minimum', fontsize=9)
        ax5.axvline(x=500, linestyle='--', color='green', alpha=0.5)
        ax5.text(520, 0.9, 'Second maximum', fontsize=9)
        
        plt.tight_layout()
        
        if save:
            plt.savefig(os.path.join(self.output_dir, 'fireball_evolution.png'), dpi=150, bbox_inches='tight')
            print(f"Saved: {self.output_dir}/fireball_evolution.png")
        
        return fig
    
    def create_animation(self, field='overpressure', fps=10, save=True):
        """Create animation of time-evolving field."""
        
        if not HAS_IMAGEIO:
            print("imageio not installed. Skipping animation.")
            return None
        
        print(f"Creating {field} animation...")
        
        x = self.data['x'] / 1000
        y = self.data['y'] / 1000
        X, Y = np.meshgrid(x, y)
        
        t = self.data['time']
        data_field = self.data[field]
        
        # Setup figure
        fig, ax = plt.subplots(figsize=(10, 8))
        
        frames = []
        
        # Create frames
        n_frames = min(100, len(t))
        frame_indices = np.linspace(0, len(t)-1, n_frames, dtype=int)
        
        vmax = np.percentile(data_field, 99)
        vmin = 0
        
        for i, ti in enumerate(frame_indices):
            ax.clear()
            
            im = ax.contourf(X, Y, data_field[ti], levels=30,
                           cmap=self.cmap_pressure, vmin=vmin, vmax=vmax)
            ax.set_title(f'{field.capitalize()} at t = {t[ti]:.2f} s')
            ax.set_xlabel('X (km)')
            ax.set_ylabel('Y (km)')
            ax.plot(50, 50, 'w*', markersize=12, markeredgecolor='k')
            
            # Add shock front
            shock_r = self.physics.shock_radius(t[ti]) / 1000
            circle = Circle((50, 50), shock_r, fill=False, color='white', 
                           linestyle='--', linewidth=2)
            ax.add_patch(circle)
            
            # Capture frame
            fig.canvas.draw()
            image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
            image = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))
            frames.append(image)
            
            if (i + 1) % 10 == 0:
                print(f"  Frame {i+1}/{n_frames}")
        
        plt.close(fig)
        
        if save and frames:
            output_path = os.path.join(self.output_dir, f'{field}_animation.gif')
            imageio.mimsave(output_path, frames, fps=fps)
            print(f"Saved: {output_path}")
        
        return frames
    
    def generate_all_plots(self):
        """Generate all visualization plots."""
        
        print("\nGenerating all visualizations...")
        print("="*50)
        
        self.plot_overview()
        self.plot_damage_zones()
        self.plot_blast_wave_evolution()
        self.plot_emp_analysis()
        self.plot_seismic_analysis()
        self.plot_radiation_analysis()
        self.plot_fireball_evolution()
        
        if HAS_IMAGEIO:
            self.create_animation('overpressure')
        
        print("="*50)
        print(f"All visualizations saved to: {self.output_dir}/")


# =============================================================================
# Main Execution
# =============================================================================

def main():
    """Main execution function."""
    
    parser = argparse.ArgumentParser(
        description='Visualize atmospheric nuclear test simulation results',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python visualize_atmospheric_nuclear.py output/atmospheric_nuclear
  python visualize_atmospheric_nuclear.py --synthetic --yield 500 --burst-height 500
  python visualize_atmospheric_nuclear.py output/run1 --save-all --output-dir figures
        """
    )
    
    parser.add_argument('data_dir', nargs='?', default=None,
                       help='Directory containing simulation output')
    parser.add_argument('--synthetic', action='store_true',
                       help='Generate synthetic data for demonstration')
    parser.add_argument('--yield', type=float, default=500.0, dest='yield_kt',
                       help='Weapon yield in kilotons (default: 500)')
    parser.add_argument('--burst-height', type=float, default=500.0,
                       help='Height of burst in meters (default: 500)')
    parser.add_argument('--output-dir', type=str, default='figures',
                       help='Output directory for figures (default: figures)')
    parser.add_argument('--save-all', action='store_true',
                       help='Generate and save all plots')
    parser.add_argument('--animation', action='store_true',
                       help='Generate animations (slower)')
    parser.add_argument('--interactive', action='store_true',
                       help='Show interactive plots')
    
    args = parser.parse_args()
    
    # Initialize physics model
    physics = NuclearExplosionPhysics(
        yield_kt=args.yield_kt,
        burst_height=args.burst_height
    )
    
    print(f"\n{'='*60}")
    print("Atmospheric Nuclear Test Visualization")
    print(f"{'='*60}")
    print(f"Yield: {args.yield_kt:.0f} kt")
    print(f"Burst Height: {args.burst_height:.0f} m")
    print(f"Total Energy: {physics.total_energy:.2e} J")
    print(f"Max Fireball Radius: {physics.fireball_radius(1):.0f} m")
    print(f"{'='*60}\n")
    
    # Load or generate data
    data = None
    
    if args.data_dir and os.path.exists(args.data_dir):
        print(f"Loading data from: {args.data_dir}")
        data = load_hdf5_data(args.data_dir)
    
    if data is None or args.synthetic:
        print("Using synthetic data for demonstration...")
        data = generate_synthetic_data(physics)
    
    # Create visualizer
    viz = AtmosphericNuclearVisualizer(data, physics, output_dir=args.output_dir)
    
    # Generate plots
    if args.save_all:
        viz.generate_all_plots()
    else:
        # Generate key plots
        viz.plot_overview()
        viz.plot_damage_zones()
        viz.plot_emp_analysis()
        viz.plot_radiation_analysis()
    
    if args.animation:
        viz.create_animation('overpressure')
        viz.create_animation('ground_velocity')
    
    if args.interactive:
        plt.show()
    
    print("\nVisualization complete!")


if __name__ == '__main__':
    main()
