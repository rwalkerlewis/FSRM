#!/usr/bin/env python3
"""
Run Underground Explosion CRAM3D Example and Generate Visualizations

This script demonstrates the near-field deformation effects for underground
explosions based on the CRAM3D methodology.

Usage:
    python run_explosion_example.py [--yield YIELD_KT] [--depth DEPTH_M]
    
Example:
    python run_explosion_example.py --yield 150 --depth 1000
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Wedge
from matplotlib.colors import LinearSegmentedColormap
import argparse
import os
import sys

# Physical constants
JOULES_PER_KT = 4.184e12

def compute_cavity_radius(yield_kt, density=2650.0):
    """
    Compute cavity radius using NTS empirical scaling.
    R_c = 55 * W^0.295 * (ρ/2.65)^(-1/3.4)
    """
    rho_ratio = density / 2650.0
    return 55.0 * (yield_kt ** 0.295) * (rho_ratio ** (-1.0/3.4))

def compute_zone_radii(cavity_radius):
    """Compute damage zone radii from cavity radius."""
    return {
        'cavity': cavity_radius,
        'crushed': 2.5 * cavity_radius,
        'fractured': 5.0 * cavity_radius,
        'damaged': 10.0 * cavity_radius
    }

def shock_pressure(r, r0, P0, n_strong=3.0, n_weak=1.87, r_transition=None):
    """
    Compute peak shock pressure at distance r.
    P(r) = P0 * (r0/r)^n with transition from strong to weak shock.
    """
    if r_transition is None:
        r_transition = 5.0 * r0
    
    if r < r0:
        return P0
    
    if r < r_transition:
        n = n_strong
    else:
        # Smooth transition
        f = min(1.0, (r - r_transition) / (2.0 * r_transition))
        n = n_strong + (n_weak - n_strong) * f
    
    return P0 * (r0 / r) ** n

def compute_damage_profile(r, zones):
    """Compute damage at radius r based on zone boundaries."""
    R_c = zones['cavity']
    R_cr = zones['crushed']
    R_fr = zones['fractured']
    R_dm = zones['damaged']
    
    if r < R_c:
        return 1.0  # Cavity
    elif r < R_cr:
        return 0.9 + 0.09 * (1.0 - (r - R_c) / (R_cr - R_c))
    elif r < R_fr:
        return 0.4 + 0.5 * (1.0 - (r - R_cr) / (R_fr - R_cr))
    elif r < R_dm:
        return 0.1 + 0.3 * (1.0 - (r - R_fr) / (R_dm - R_fr))
    else:
        return 0.0

def compute_seismic_moment(yield_kt, density=2650.0, vp=5500.0):
    """
    Compute scalar seismic moment using empirical relation.
    log10(M0) ≈ 17.0 + log10(W_kt) for contained tests.
    """
    return 10.0 ** (17.0 + np.log10(yield_kt))

def moment_magnitude(M0):
    """Convert scalar moment to moment magnitude."""
    return (np.log10(M0) - 9.1) / 1.5

def body_wave_magnitude(yield_kt):
    """Compute body wave magnitude mb from yield."""
    return 4.0 + 0.75 * np.log10(yield_kt)

def rdp_source_time_function(t, corner_freq, overshoot=1.2):
    """
    Compute RDP (Reduced Displacement Potential) source time function.
    Brune omega-squared model.
    """
    tau = 1.0 / (2.0 * np.pi * corner_freq)
    t_norm = t / tau
    
    if isinstance(t, np.ndarray):
        psi = np.zeros_like(t)
        mask = t > 0
        psi[mask] = (1.0 - (1.0 + t_norm[mask]) * np.exp(-t_norm[mask])) * overshoot
        
        psi_dot = np.zeros_like(t)
        psi_dot[mask] = t_norm[mask] * np.exp(-t_norm[mask]) * overshoot / tau
        
        return psi, psi_dot
    else:
        if t <= 0:
            return 0.0, 0.0
        psi = (1.0 - (1.0 + t_norm) * np.exp(-t_norm)) * overshoot
        psi_dot = t_norm * np.exp(-t_norm) * overshoot / tau
        return psi, psi_dot

def create_damage_zone_plot(zones, depth, yield_kt, save_path=None):
    """Create a cross-section visualization of damage zones."""
    fig, ax = plt.subplots(figsize=(12, 10))
    
    R_dm = zones['damaged']
    extent = 1.5 * R_dm
    
    # Colors for zones
    colors = {
        'intact': '#4477AA',
        'damaged': '#CCCC00',
        'fractured': '#FF8800',
        'crushed': '#CC3311',
        'cavity': '#000000'
    }
    
    # Draw zones as filled circles (from outside in)
    ax.add_patch(Circle((0, -depth), zones['damaged'], 
                        fc=colors['damaged'], ec='black', lw=2, zorder=1,
                        label=f"Damaged zone (r={zones['damaged']:.0f}m)"))
    ax.add_patch(Circle((0, -depth), zones['fractured'], 
                        fc=colors['fractured'], ec='black', lw=2, zorder=2,
                        label=f"Fractured zone (r={zones['fractured']:.0f}m)"))
    ax.add_patch(Circle((0, -depth), zones['crushed'], 
                        fc=colors['crushed'], ec='black', lw=2, zorder=3,
                        label=f"Crushed zone (r={zones['crushed']:.0f}m)"))
    ax.add_patch(Circle((0, -depth), zones['cavity'], 
                        fc=colors['cavity'], ec='white', lw=2, zorder=4,
                        label=f"Cavity (r={zones['cavity']:.0f}m)"))
    
    # Surface line
    ax.axhline(y=0, color='brown', lw=3, zorder=5, label='Surface')
    
    # Fill above surface
    ax.fill_between([-extent, extent], [0, 0], [0.2*extent, 0.2*extent], 
                   color='#DDCCAA', zorder=0, alpha=0.5)
    
    # Fill below damaged zone
    ax.fill_between([-extent, extent], 
                   [-depth - R_dm, -depth - R_dm],
                   [-depth - 1.5*R_dm, -depth - 1.5*R_dm],
                   color=colors['intact'], zorder=0, alpha=0.3)
    
    # Ground zero marker
    ax.plot(0, 0, 'v', color='red', markersize=15, zorder=10, label='Ground Zero')
    ax.plot(0, -depth, '*', color='white', markersize=20, zorder=10, 
            markeredgecolor='black', markeredgewidth=2, label='Detonation Point')
    
    # Labels
    ax.set_xlabel('Distance from Ground Zero (m)', fontsize=12)
    ax.set_ylabel('Depth (m)', fontsize=12)
    ax.set_title(f'Underground Explosion Damage Zones\n{yield_kt} kt at {depth} m depth', 
                fontsize=14, fontweight='bold')
    
    ax.set_xlim(-extent, extent)
    ax.set_ylim(-depth - 1.3*R_dm, 0.3*extent)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper right', fontsize=10)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"  Saved: {save_path}")
    
    return fig, ax

def create_pressure_attenuation_plot(zones, yield_kt, density=2650.0, vp=5500.0, 
                                     save_path=None):
    """Create shock pressure attenuation plot."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    R_c = zones['cavity']
    P0 = 100.0e9  # 100 GPa initial
    
    # Distance array
    r = np.logspace(np.log10(R_c), np.log10(50 * zones['damaged']), 500)
    
    # Compute pressure, velocity, damage
    P = np.array([shock_pressure(ri, R_c, P0) for ri in r])
    v = P / (density * vp)  # Particle velocity
    D = np.array([compute_damage_profile(ri, zones) for ri in r])
    
    # Arrival time (approximate)
    t_arr = (r - R_c) / vp
    
    # Plot 1: Pressure attenuation
    ax = axes[0, 0]
    ax.loglog(r, P / 1e9, 'b-', lw=2, label='Peak Pressure')
    
    # Mark zones
    for name, radius, color in [('Cavity', zones['cavity'], 'black'),
                                ('Crushed', zones['crushed'], 'red'),
                                ('Fractured', zones['fractured'], 'orange'),
                                ('Damaged', zones['damaged'], 'yellow')]:
        ax.axvline(radius, color=color, ls='--', lw=1.5, label=f'{name} ({radius:.0f}m)')
    
    # Pressure thresholds
    ax.axhline(5.0, color='gray', ls=':', label='Crushing (5 GPa)')
    ax.axhline(0.5, color='gray', ls=':', alpha=0.7, label='Fracture (0.5 GPa)')
    
    ax.set_xlabel('Distance from Source (m)', fontsize=11)
    ax.set_ylabel('Peak Pressure (GPa)', fontsize=11)
    ax.set_title('Shock Pressure Attenuation', fontsize=12, fontweight='bold')
    ax.legend(fontsize=8, loc='upper right')
    ax.grid(True, alpha=0.3)
    
    # Plot 2: Particle velocity
    ax = axes[0, 1]
    ax.loglog(r, v, 'g-', lw=2, label='Peak Velocity')
    ax.axhline(1.0, color='red', ls='--', label='Spall threshold (1 m/s)')
    
    for name, radius, color in [('Crushed', zones['crushed'], 'red'),
                                ('Fractured', zones['fractured'], 'orange'),
                                ('Damaged', zones['damaged'], 'yellow')]:
        ax.axvline(radius, color=color, ls='--', lw=1.5, alpha=0.7)
    
    ax.set_xlabel('Distance from Source (m)', fontsize=11)
    ax.set_ylabel('Peak Particle Velocity (m/s)', fontsize=11)
    ax.set_title('Ground Motion', fontsize=12, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    
    # Plot 3: Damage distribution
    ax = axes[1, 0]
    ax.semilogx(r, D, 'r-', lw=2)
    ax.fill_between(r, 0, D, alpha=0.3, color='red')
    
    for name, radius, color in [('Cavity', zones['cavity'], 'black'),
                                ('Crushed', zones['crushed'], 'red'),
                                ('Fractured', zones['fractured'], 'orange'),
                                ('Damaged', zones['damaged'], 'yellow')]:
        ax.axvline(radius, color=color, ls='--', lw=1.5)
    
    ax.set_xlabel('Distance from Source (m)', fontsize=11)
    ax.set_ylabel('Damage Parameter D [-]', fontsize=11)
    ax.set_title('Radial Damage Distribution', fontsize=12, fontweight='bold')
    ax.set_ylim(0, 1.05)
    ax.grid(True, alpha=0.3)
    
    # Plot 4: Arrival time
    ax = axes[1, 1]
    ax.loglog(r, t_arr, 'm-', lw=2, label='Shock Arrival')
    ax.set_xlabel('Distance from Source (m)', fontsize=11)
    ax.set_ylabel('Arrival Time (s)', fontsize=11)
    ax.set_title('Shock Wave Propagation', fontsize=12, fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    fig.suptitle(f'Underground Explosion: {yield_kt} kt TNT', 
                fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"  Saved: {save_path}")
    
    return fig, axes

def create_source_time_function_plot(yield_kt, save_path=None):
    """Create seismic source time function plot."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Corner frequency from yield
    corner_freq = 2.5 * (yield_kt ** (-1.0/3.0))
    
    # Time array
    duration = 10.0 / corner_freq
    t = np.linspace(-0.1 * duration, duration, 500)
    
    # Compute source time function
    psi, psi_dot = rdp_source_time_function(t, corner_freq)
    
    # Normalize for plotting
    psi_norm = psi / np.max(psi) if np.max(psi) > 0 else psi
    psi_dot_norm = psi_dot / np.max(np.abs(psi_dot)) if np.max(np.abs(psi_dot)) > 0 else psi_dot
    
    # Plot 1: RDP and its derivative
    ax = axes[0]
    ax.plot(t, psi_norm, 'b-', lw=2, label='Normalized ψ(t)')
    ax.plot(t, psi_dot_norm, 'r--', lw=2, label='Normalized dψ/dt')
    ax.axvline(0, color='gray', ls=':', alpha=0.5)
    ax.axhline(0, color='gray', ls=':', alpha=0.5)
    
    ax.set_xlabel('Time (s)', fontsize=11)
    ax.set_ylabel('Normalized Amplitude', fontsize=11)
    ax.set_title('Reduced Displacement Potential', fontsize=12, fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot 2: Frequency spectrum
    ax = axes[1]
    f = np.logspace(-2, 2, 500)
    omega = 2 * np.pi * f
    omega_c = 2 * np.pi * corner_freq
    
    # Brune spectrum: 1 / (1 + (ω/ωc)²)
    spectrum = 1.0 / (1.0 + (omega / omega_c) ** 2)
    
    ax.loglog(f, spectrum, 'b-', lw=2)
    ax.axvline(corner_freq, color='red', ls='--', lw=2, 
              label=f'Corner freq = {corner_freq:.2f} Hz')
    ax.axvline(1.0, color='gray', ls=':', alpha=0.5, label='1 Hz reference')
    
    ax.set_xlabel('Frequency (Hz)', fontsize=11)
    ax.set_ylabel('Spectral Amplitude', fontsize=11)
    ax.set_title('Source Spectrum (Brune Model)', fontsize=12, fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    fig.suptitle(f'Seismic Source: {yield_kt} kt (f_c = {corner_freq:.2f} Hz)', 
                fontsize=13, fontweight='bold', y=1.02)
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"  Saved: {save_path}")
    
    return fig, axes

def create_summary_plot(yield_kt, depth, zones, save_path=None):
    """Create comprehensive summary visualization."""
    fig = plt.figure(figsize=(16, 12))
    
    # Seismic parameters
    M0 = compute_seismic_moment(yield_kt)
    Mw = moment_magnitude(M0)
    mb = body_wave_magnitude(yield_kt)
    fc = 2.5 * (yield_kt ** (-1.0/3.0))
    energy = yield_kt * JOULES_PER_KT
    
    # Create subplots
    gs = fig.add_gridspec(3, 3, hspace=0.35, wspace=0.3)
    
    # Plot 1: Damage zones (large)
    ax1 = fig.add_subplot(gs[0:2, 0:2])
    
    R_dm = zones['damaged']
    extent = 1.3 * R_dm
    
    # Draw zones
    colors = {'damaged': '#CCCC00', 'fractured': '#FF8800', 
              'crushed': '#CC3311', 'cavity': '#000000'}
    
    ax1.add_patch(Circle((0, -depth), zones['damaged'], 
                         fc=colors['damaged'], ec='black', lw=2))
    ax1.add_patch(Circle((0, -depth), zones['fractured'], 
                         fc=colors['fractured'], ec='black', lw=2))
    ax1.add_patch(Circle((0, -depth), zones['crushed'], 
                         fc=colors['crushed'], ec='black', lw=2))
    ax1.add_patch(Circle((0, -depth), zones['cavity'], 
                         fc=colors['cavity'], ec='white', lw=2))
    
    ax1.axhline(0, color='brown', lw=3)
    ax1.fill_between([-extent, extent], [0, 0], [0.1*extent, 0.1*extent], 
                    color='#DDCCAA', alpha=0.5)
    
    ax1.plot(0, 0, 'rv', markersize=12)
    ax1.plot(0, -depth, 'w*', markersize=15, markeredgecolor='black', markeredgewidth=2)
    
    ax1.set_xlim(-extent, extent)
    ax1.set_ylim(-depth - 1.1*R_dm, 0.15*extent)
    ax1.set_aspect('equal')
    ax1.set_xlabel('Distance from GZ (m)')
    ax1.set_ylabel('Depth (m)')
    ax1.set_title('Damage Zone Cross-Section', fontweight='bold')
    ax1.grid(True, alpha=0.3)
    
    # Legend for zones
    legend_elements = [plt.Rectangle((0,0), 1, 1, fc=c, ec='black') 
                       for c in colors.values()]
    ax1.legend(legend_elements, 
              [f'Damaged ({zones["damaged"]:.0f}m)',
               f'Fractured ({zones["fractured"]:.0f}m)',
               f'Crushed ({zones["crushed"]:.0f}m)',
               f'Cavity ({zones["cavity"]:.0f}m)'],
              loc='upper right', fontsize=9)
    
    # Plot 2: Pressure attenuation
    ax2 = fig.add_subplot(gs[0, 2])
    
    R_c = zones['cavity']
    r = np.logspace(np.log10(R_c), np.log10(20*R_dm), 200)
    P = np.array([shock_pressure(ri, R_c, 100e9) for ri in r])
    
    ax2.loglog(r, P/1e9, 'b-', lw=2)
    ax2.axvline(zones['crushed'], color='red', ls='--', alpha=0.7)
    ax2.axvline(zones['fractured'], color='orange', ls='--', alpha=0.7)
    ax2.set_xlabel('Distance (m)')
    ax2.set_ylabel('Pressure (GPa)')
    ax2.set_title('Pressure Decay', fontweight='bold')
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Damage profile
    ax3 = fig.add_subplot(gs[1, 2])
    
    D = np.array([compute_damage_profile(ri, zones) for ri in r])
    ax3.semilogx(r, D, 'r-', lw=2)
    ax3.fill_between(r, 0, D, alpha=0.3, color='red')
    ax3.set_xlabel('Distance (m)')
    ax3.set_ylabel('Damage')
    ax3.set_ylim(0, 1.05)
    ax3.set_title('Damage Distribution', fontweight='bold')
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Source time function
    ax4 = fig.add_subplot(gs[2, 0])
    
    duration = 10.0 / fc
    t = np.linspace(0, duration, 200)
    psi, psi_dot = rdp_source_time_function(t, fc)
    psi_norm = psi / np.max(psi) if np.max(psi) > 0 else psi
    psi_dot_norm = psi_dot / np.max(np.abs(psi_dot)) if np.max(np.abs(psi_dot)) > 0 else psi_dot
    
    ax4.plot(t, psi_norm, 'b-', lw=2, label='ψ(t)')
    ax4.plot(t, psi_dot_norm, 'r--', lw=1.5, label='dψ/dt')
    ax4.set_xlabel('Time (s)')
    ax4.set_ylabel('Normalized')
    ax4.set_title(f'Source Function (f_c={fc:.2f} Hz)', fontweight='bold')
    ax4.legend(fontsize=8)
    ax4.grid(True, alpha=0.3)
    
    # Plot 5: Energy partitioning (bar chart)
    ax5 = fig.add_subplot(gs[2, 1])
    
    categories = ['Internal\n(thermal)', 'Kinetic', 'Plastic\n(damage)', 'Seismic', 'Fracture']
    # Approximate energy partitioning
    fractions = [75, 10, 10, 0.01, 4.99]
    
    bars = ax5.bar(categories, fractions, color=['red', 'blue', 'orange', 'green', 'purple'])
    ax5.set_ylabel('Energy Fraction (%)')
    ax5.set_title('Energy Partitioning', fontweight='bold')
    ax5.set_ylim(0, 100)
    
    # Annotate small bars
    for bar, frac in zip(bars, fractions):
        if frac < 5:
            ax5.annotate(f'{frac:.2f}%', xy=(bar.get_x() + bar.get_width()/2, bar.get_height()),
                        ha='center', va='bottom', fontsize=8)
    
    # Plot 6: Summary text
    ax6 = fig.add_subplot(gs[2, 2])
    ax6.axis('off')
    
    summary = f"""
    EXPLOSION PARAMETERS
    ━━━━━━━━━━━━━━━━━━━━
    Yield:           {yield_kt:.1f} kt TNT
    Energy:          {energy:.2e} J
    Depth:           {depth:.0f} m
    
    DAMAGE ZONES
    ━━━━━━━━━━━━━━━━━━━━
    Cavity:          {zones['cavity']:.1f} m
    Crushed zone:    {zones['crushed']:.1f} m
    Fractured zone:  {zones['fractured']:.1f} m
    Damaged zone:    {zones['damaged']:.1f} m
    
    SEISMIC SOURCE
    ━━━━━━━━━━━━━━━━━━━━
    Scalar moment:   {M0:.2e} N·m
    M_w:             {Mw:.2f}
    m_b:             {mb:.2f}
    Corner freq:     {fc:.2f} Hz
    """
    
    ax6.text(0.05, 0.95, summary, transform=ax6.transAxes, 
            fontfamily='monospace', fontsize=9,
            verticalalignment='top')
    
    fig.suptitle(f'CRAM3D Underground Explosion Analysis: {yield_kt} kt at {depth} m', 
                fontsize=14, fontweight='bold', y=0.98)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"  Saved: {save_path}")
    
    return fig

def main():
    parser = argparse.ArgumentParser(
        description='Run CRAM3D Underground Explosion Example')
    parser.add_argument('--yield', '-y', type=float, default=150.0,
                       dest='yield_kt', help='Yield in kilotons (default: 150)')
    parser.add_argument('--depth', '-d', type=float, default=1000.0,
                       help='Depth of burial in meters (default: 1000)')
    parser.add_argument('--output-dir', '-o', type=str, default='.',
                       help='Output directory for plots')
    parser.add_argument('--no-show', action='store_true',
                       help='Do not display plots interactively')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Print header
    print("=" * 70)
    print("  CRAM3D Underground Explosion Simulation")
    print("  Near-field Damage and Seismic Wave Generation")
    print("=" * 70)
    print()
    
    # Compute parameters
    print(f"Configuration:")
    print(f"  Yield:  {args.yield_kt} kt TNT")
    print(f"  Depth:  {args.depth} m")
    print()
    
    # Compute cavity and zones
    cavity_radius = compute_cavity_radius(args.yield_kt)
    zones = compute_zone_radii(cavity_radius)
    
    print("Computed Zone Radii (NTS empirical scaling):")
    print(f"  Cavity radius:     {zones['cavity']:.1f} m")
    print(f"  Crushed zone:      {zones['crushed']:.1f} m")
    print(f"  Fractured zone:    {zones['fractured']:.1f} m")
    print(f"  Damaged zone:      {zones['damaged']:.1f} m")
    print()
    
    # Seismic source
    M0 = compute_seismic_moment(args.yield_kt)
    Mw = moment_magnitude(M0)
    mb = body_wave_magnitude(args.yield_kt)
    fc = 2.5 * (args.yield_kt ** (-1.0/3.0))
    
    print("Seismic Source Parameters:")
    print(f"  Scalar moment:     {M0:.2e} N·m")
    print(f"  Moment magnitude:  Mw {Mw:.2f}")
    print(f"  Body wave mag:     mb {mb:.2f}")
    print(f"  Corner frequency:  {fc:.2f} Hz")
    print()
    
    # Generate plots
    print("Generating visualizations...")
    
    base = os.path.join(args.output_dir, f'explosion_{args.yield_kt:.0f}kt')
    
    # Zone plot
    create_damage_zone_plot(zones, args.depth, args.yield_kt,
                           save_path=f'{base}_zones.png')
    
    # Attenuation plot
    create_pressure_attenuation_plot(zones, args.yield_kt,
                                     save_path=f'{base}_attenuation.png')
    
    # Source time function
    create_source_time_function_plot(args.yield_kt,
                                     save_path=f'{base}_source.png')
    
    # Summary
    create_summary_plot(args.yield_kt, args.depth, zones,
                       save_path=f'{base}_summary.png')
    
    print()
    print("=" * 70)
    print("  Complete!")
    print("=" * 70)
    print()
    print(f"Output files in: {args.output_dir}")
    print(f"  {base}_zones.png      - Damage zone cross-section")
    print(f"  {base}_attenuation.png - Shock attenuation curves")
    print(f"  {base}_source.png      - Seismic source function")
    print(f"  {base}_summary.png     - Complete summary")
    print()
    
    if not args.no_show:
        print("Displaying plots... (close windows to exit)")
        plt.show()

if __name__ == '__main__':
    main()
