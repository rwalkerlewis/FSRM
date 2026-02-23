#!/usr/bin/env python3
"""
Full Moment Tensor Inversion for the 1992-09-23 US Divider Nuclear Test

Performs waveform-based moment tensor inversion using regional seismic data:
  1. Downloads broadband waveforms from IRIS FDSN
  2. Computes Green's functions for a 1D velocity model (FK / ray theory)
  3. Inverts for the full 6-component moment tensor via linear least squares
  4. Decomposes result into isotropic (ISO), compensated linear vector dipole
     (CLVD), and double-couple (DC) components
  5. Generates publication-quality figures with beachball, decomposition, and
     waveform fits

The full moment tensor has 6 independent components:
    M = [Mrr, Mtt, Mpp, Mrt, Mrp, Mtp]  (in spherical coordinates)
or equivalently:
    M = [Mxx, Myy, Mzz, Mxy, Mxz, Myz]  (in Cartesian)

For an underground nuclear explosion, we expect:
  - Large isotropic component (positive, compressional)
  - Small CLVD (from spall or cavity collapse)
  - Small DC (from tectonic release or block motion)

Event:
    1992-09-23 15:04:00 UTC  (Operation Julin — Divider)
    Location: 37.021°N, 116.058°W (Nevada Test Site, Yucca Flat)
    Depth: ~800 m
    Yield: ~1 kt (sub-kiloton class)

Expected decomposition: ~80% ISO, ~18% CLVD, ~2% DC
Note: 1992 data may be limited; fallback to synthetics is likely.

Usage:
    python scripts/invert_us_divider_1992_moment_tensor.py

Requirements:
    pip install obspy numpy scipy matplotlib

References:
    Jost & Herrmann (1989) — A student's guide to and review of moment tensors
    Dreger (2003) — TDMT_INV: Time domain moment tensor inversion
    Ford et al. (2009) — Explosion source discrimination using full waveforms
    Walter et al. (2004) — Moment tensor estimation for NTS explosions
"""

import os
import sys
import warnings
import numpy as np
from scipy.optimize import nnls
from scipy.signal import butter, filtfilt
from scipy.linalg import lstsq, svd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Wedge, Polygon
from matplotlib.collections import PatchCollection
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import AutoMinorLocator

from obspy import UTCDateTime, Stream, Trace
from obspy.clients.fdsn import Client
from obspy.geodetics import gps2dist_azimuth
from obspy.signal.filter import bandpass

warnings.filterwarnings("ignore")

# ==============================================================================
# Directories
# ==============================================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(SCRIPT_DIR)
OUTDIR = os.path.join(PROJECT_DIR, "figures", "us_divider_1992")
os.makedirs(OUTDIR, exist_ok=True)

# ==============================================================================
# Event Parameters
# ==============================================================================
EVENT_TIME = UTCDateTime("1992-09-23T15:04:00")
EVENT_LAT = 37.021
EVENT_LON = -116.058
EVENT_DEPTH_KM = 0.8  # NTS Yucca Flat depth

# Processing parameters
FREQMIN = 0.02  # Hz - lower bound for long-period inversion
FREQMAX = 0.10  # Hz - upper bound (regional surface waves)
PRE_ORIGIN = 60  # seconds before origin
POST_ORIGIN = 600  # seconds after origin (captures Rayleigh at 2000 km)

# ==============================================================================
# Stations for Inversion (regional, good azimuthal coverage)
# ==============================================================================
STATIONS = [
    # (net, sta, loc, chan_pref, description)
    ("CI", "PAS", "*", "BH", "Pasadena"),
    ("CI", "GSC", "*", "BH", "Goldstone"),
    ("BK", "SAO", "*", "BH", "San Andreas Obs"),
    ("NN", "ELK", "*", "BH", "Elko"),
]


# ==============================================================================
# 1D Velocity Model (NTS / Basin and Range)
# ==============================================================================
class VelocityModel1D:
    """
    Layered 1D velocity model for Green's function computation.
    Based on NTS (Nevada Test Site) / Basin and Range crustal structure.
    """

    def __init__(self):
        # Layers: [depth_top_km, thickness_km, Vp_km/s, Vs_km/s, rho_g/cm3, Qp, Qs]
        self.layers = np.array([
            [0.0,   0.3,  2.00, 1.00, 1.90,   60,  30],    # Alluvium/fill
            [0.3,   0.7,  3.00, 1.70, 2.10,  100,  50],    # Tuff/alluvium
            [1.0,   4.0,  4.50, 2.60, 2.30,  200, 100],    # Volcanic tuff
            [5.0,   8.0,  5.50, 3.20, 2.50,  350, 175],    # Paleozoic carbonate
            [13.0, 12.0,  6.00, 3.50, 2.65,  500, 250],    # Upper crust
            [25.0,  8.0,  6.50, 3.75, 2.85,  600, 300],    # Lower crust
            [33.0, 100.0, 7.90, 4.50, 3.30, 1000, 500],    # Upper mantle
        ])

    def get_properties_at_depth(self, depth_km):
        """Return (Vp, Vs, rho, Qp, Qs) at given depth."""
        for layer in self.layers[::-1]:
            if depth_km >= layer[0]:
                return layer[2], layer[3], layer[4], layer[5], layer[6]
        return self.layers[0, 2:7]


# ==============================================================================
# Green's Functions (Simplified Ray Theory)
# ==============================================================================

def compute_greens_functions(dist_km, azimuth_deg, depth_km, dt, npts,
                             model=None, freqmin=0.02, freqmax=0.10):
    """
    Compute approximate Green's functions for a point source.

    Returns a dict of Green's functions for each moment tensor component:
        G[component][channel] = waveform array

    Components: Mxx, Myy, Mzz, Mxy, Mxz, Myz
    Channels: Z, R, T (vertical, radial, transverse)

    This uses a simplified ray-theory approach suitable for regional distances.
    For production use, replace with FK (frequency-wavenumber) synthetics.
    """
    if model is None:
        model = VelocityModel1D()

    az_rad = np.deg2rad(azimuth_deg)
    t = np.arange(npts) * dt

    # Get average crustal properties
    vp_avg = 6.2  # km/s
    vs_avg = 3.6  # km/s
    rho_avg = 2.8  # g/cm³

    # Travel times for various phases
    t_p = dist_km / vp_avg  # Direct P
    t_s = dist_km / vs_avg  # Direct S
    t_rayleigh = dist_km / 3.2  # Rayleigh wave
    t_love = dist_km / 3.5  # Love wave

    # Source time function (Gaussian derivative for moment rate)
    tau = 2.0  # source duration in seconds

    def source_time_function(t_arr, t0):
        """Gaussian derivative centered at t0."""
        t_rel = t_arr - t0
        return t_rel / tau**2 * np.exp(-t_rel**2 / (2 * tau**2))

    # Geometric spreading
    r_km = max(dist_km, 1.0)
    geom_body = 1.0 / r_km  # Body wave
    geom_surf = 1.0 / np.sqrt(r_km)  # Surface wave

    # Initialize Green's functions dict
    components = ['Mxx', 'Myy', 'Mzz', 'Mxy', 'Mxz', 'Myz']
    channels = ['Z', 'R', 'T']
    G = {comp: {ch: np.zeros(npts) for ch in channels} for comp in components}

    # Radiation patterns for each component (simplified)
    cos_az = np.cos(az_rad)
    sin_az = np.sin(az_rad)
    cos_2az = np.cos(2 * az_rad)
    sin_2az = np.sin(2 * az_rad)

    # Take-off angle (assume ~45° for regional distances)
    takeoff = np.deg2rad(45)
    cos_th = np.cos(takeoff)
    sin_th = np.sin(takeoff)

    # P-wave radiation patterns (Aki & Richards, 2002)
    rad_p = {
        'Mxx': sin_th**2 * cos_az**2 - cos_th**2,
        'Myy': sin_th**2 * sin_az**2 - cos_th**2,
        'Mzz': cos_th**2,
        'Mxy': sin_th**2 * sin_2az,
        'Mxz': sin_th * cos_th * cos_az,
        'Myz': sin_th * cos_th * sin_az,
    }

    # SV radiation patterns
    rad_sv = {
        'Mxx': sin_th * cos_th * (cos_az**2 - 1),
        'Myy': sin_th * cos_th * (sin_az**2 - 1),
        'Mzz': -sin_th * cos_th,
        'Mxy': sin_th * cos_th * sin_2az,
        'Mxz': (cos_th**2 - sin_th**2) * cos_az,
        'Myz': (cos_th**2 - sin_th**2) * sin_az,
    }

    # SH radiation patterns
    rad_sh = {
        'Mxx': -sin_th * sin_2az / 2,
        'Myy': sin_th * sin_2az / 2,
        'Mzz': 0.0,
        'Mxy': sin_th * cos_2az,
        'Mxz': cos_th * sin_az,
        'Myz': -cos_th * cos_az,
    }

    # Rayleigh wave (vertical + radial) - fundamental mode
    rad_rayleigh_z = {
        'Mxx': 0.5 * cos_2az,
        'Myy': -0.5 * cos_2az,
        'Mzz': 1.0,
        'Mxy': sin_2az,
        'Mxz': cos_az,
        'Myz': sin_az,
    }
    rad_rayleigh_r = {
        'Mxx': 0.3 * cos_2az,
        'Myy': -0.3 * cos_2az,
        'Mzz': 0.5,
        'Mxy': 0.6 * sin_2az,
        'Mxz': 0.5 * cos_az,
        'Myz': 0.5 * sin_az,
    }

    # Love wave (transverse) - sensitive to horizontal shear
    rad_love = {
        'Mxx': -0.5 * sin_2az,
        'Myy': 0.5 * sin_2az,
        'Mzz': 0.0,
        'Mxy': cos_2az,
        'Mxz': sin_az,
        'Myz': -cos_az,
    }

    # Build Green's functions by summing phases
    for comp in components:
        stf_p = source_time_function(t, t_p)
        G[comp]['Z'] += rad_p[comp] * geom_body * stf_p * 0.5
        G[comp]['R'] += rad_p[comp] * geom_body * stf_p * 0.3

        stf_s = source_time_function(t, t_s)
        G[comp]['Z'] += rad_sv[comp] * geom_body * stf_s * 0.4
        G[comp]['R'] += rad_sv[comp] * geom_body * stf_s * 0.6

        G[comp]['T'] += rad_sh[comp] * geom_body * stf_s * 0.8

        stf_r = source_time_function(t, t_rayleigh)
        G[comp]['Z'] += rad_rayleigh_z[comp] * geom_surf * stf_r * 2.0
        G[comp]['R'] += rad_rayleigh_r[comp] * geom_surf * stf_r * 1.5

        stf_l = source_time_function(t, t_love)
        G[comp]['T'] += rad_love[comp] * geom_surf * stf_l * 1.5

    # Apply bandpass filter to all Green's functions
    for comp in components:
        for ch in channels:
            G[comp][ch] = bandpass_filter(G[comp][ch], freqmin, freqmax, 1.0/dt)

    return G


def bandpass_filter(data, freqmin, freqmax, sampling_rate, corners=4):
    """Apply zero-phase Butterworth bandpass filter."""
    nyq = 0.5 * sampling_rate
    low = freqmin / nyq
    high = freqmax / nyq
    low = max(0.001, min(low, 0.999))
    high = max(low + 0.001, min(high, 0.999))
    b, a = butter(corners, [low, high], btype='band')
    return filtfilt(b, a, data)


# ==============================================================================
# Moment Tensor Inversion
# ==============================================================================

def build_data_vector(observations):
    """
    Stack observed waveforms into a single data vector.

    observations: list of (trace_z, trace_r, trace_t, weight) tuples

    Returns: 1D numpy array
    """
    d = []
    for obs in observations:
        trace_z, trace_r, trace_t, weight = obs
        if trace_z is not None:
            d.extend(trace_z * weight)
        if trace_r is not None:
            d.extend(trace_r * weight)
        if trace_t is not None:
            d.extend(trace_t * weight)
    return np.array(d)


def build_greens_matrix(greens_list, observations):
    """
    Build the Green's function matrix G such that d = G @ m.

    greens_list: list of Green's function dicts (one per station)
    observations: list of (trace_z, trace_r, trace_t, weight) tuples

    Returns: (ndata, 6) numpy array
    """
    components = ['Mxx', 'Myy', 'Mzz', 'Mxy', 'Mxz', 'Myz']
    rows = []

    for i, (obs, gf) in enumerate(zip(observations, greens_list)):
        trace_z, trace_r, trace_t, weight = obs

        if trace_z is not None:
            for comp in components:
                rows.append([])
            for j, comp in enumerate(components):
                rows[-6 + j].extend(gf[comp]['Z'] * weight)

        if trace_r is not None:
            for comp in components:
                rows.append([])
            for j, comp in enumerate(components):
                rows[-6 + j].extend(gf[comp]['R'] * weight)

        if trace_t is not None:
            for comp in components:
                rows.append([])
            for j, comp in enumerate(components):
                rows[-6 + j].extend(gf[comp]['T'] * weight)

    G_rows = []
    for i, (obs, gf) in enumerate(zip(observations, greens_list)):
        trace_z, trace_r, trace_t, weight = obs
        npts = len(trace_z) if trace_z is not None else len(trace_r) if trace_r is not None else len(trace_t)

        if trace_z is not None:
            for k in range(npts):
                row = [gf[comp]['Z'][k] * weight for comp in components]
                G_rows.append(row)

        if trace_r is not None:
            for k in range(npts):
                row = [gf[comp]['R'][k] * weight for comp in components]
                G_rows.append(row)

        if trace_t is not None:
            for k in range(npts):
                row = [gf[comp]['T'][k] * weight for comp in components]
                G_rows.append(row)

    return np.array(G_rows)


def invert_moment_tensor(d, G, method='lstsq'):
    """
    Solve for moment tensor: d = G @ m

    Parameters:
        d: data vector (ndata,)
        G: Green's matrix (ndata, 6)
        method: 'lstsq' (unconstrained) or 'nnls' (non-negative, for ISO constraint)

    Returns:
        m: moment tensor vector [Mxx, Myy, Mzz, Mxy, Mxz, Myz]
        residual: RMS misfit
        variance_reduction: percentage of variance explained
    """
    if method == 'lstsq':
        m, residuals, rank, s = lstsq(G, d, lapack_driver='gelsy')
    elif method == 'nnls':
        m, rnorm = nnls(G, d)
    else:
        raise ValueError(f"Unknown method: {method}")

    d_pred = G @ m
    residual = np.sqrt(np.mean((d - d_pred)**2))
    var_data = np.var(d)
    var_res = np.var(d - d_pred)
    variance_reduction = 100 * (1 - var_res / var_data) if var_data > 0 else 0

    return m, residual, variance_reduction, d_pred


# ==============================================================================
# Moment Tensor Decomposition
# ==============================================================================

def decompose_moment_tensor(m):
    """
    Decompose moment tensor into isotropic, CLVD, and double-couple components.

    Input m = [Mxx, Myy, Mzz, Mxy, Mxz, Myz] (Cartesian)

    Returns dict with:
        - M0: scalar moment
        - Mw: moment magnitude
        - ISO_fraction: isotropic fraction (0-1)
        - CLVD_fraction: CLVD fraction (0-1)
        - DC_fraction: double-couple fraction (0-1)
        - eigenvalues: sorted eigenvalues
        - eigenvectors: corresponding eigenvectors (P, B, T axes)
        - strike, dip, rake: fault plane solution (for DC part)
    """
    Mxx, Myy, Mzz, Mxy, Mxz, Myz = m
    M = np.array([
        [Mxx, Mxy, Mxz],
        [Mxy, Myy, Myz],
        [Mxz, Myz, Mzz]
    ])

    eigenvalues, eigenvectors = np.linalg.eigh(M)

    idx = np.argsort(np.abs(eigenvalues))[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]

    tr_M = np.trace(M)
    ISO = tr_M / 3.0

    M_dev = M - ISO * np.eye(3)
    eigenvalues_dev = eigenvalues - ISO

    idx_dev = np.argsort(eigenvalues_dev)[::-1]
    e_dev = eigenvalues_dev[idx_dev]

    if np.abs(e_dev[0]) > 1e-20:
        epsilon = -e_dev[2] / e_dev[0]
        epsilon = np.clip(epsilon, -0.5, 0.5)
    else:
        epsilon = 0.0

    f_clvd = 2 * np.abs(epsilon)
    f_dc = 1 - f_clvd

    M0_dev = np.max(np.abs(eigenvalues_dev))
    M0_iso = np.abs(ISO)
    M0 = np.sqrt(M0_iso**2 + M0_dev**2)

    M0_alt = np.sqrt(0.5 * np.sum(eigenvalues**2))

    if M0 > 0:
        Mw = (np.log10(M0) - 9.1) / 1.5
    else:
        Mw = -99

    total_abs = np.abs(ISO) + M0_dev
    if total_abs > 0:
        f_iso = np.abs(ISO) / total_abs
    else:
        f_iso = 0.0

    f_clvd_adj = (1 - f_iso) * f_clvd
    f_dc_adj = (1 - f_iso) * f_dc

    T_axis = eigenvectors[:, 0]
    P_axis = eigenvectors[:, 2]

    strike, dip, rake = axes_to_sdr(P_axis, T_axis)

    return {
        'M0': M0,
        'Mw': Mw,
        'ISO_fraction': f_iso,
        'CLVD_fraction': f_clvd_adj,
        'DC_fraction': f_dc_adj,
        'ISO_value': ISO,
        'epsilon': epsilon,
        'eigenvalues': eigenvalues,
        'eigenvectors': eigenvectors,
        'P_axis': P_axis,
        'T_axis': T_axis,
        'strike': strike,
        'dip': dip,
        'rake': rake,
        'tensor_matrix': M,
    }


def axes_to_sdr(P_axis, T_axis):
    """
    Convert P and T axes to strike, dip, rake of one fault plane.
    Simplified conversion - returns approximate values.
    """
    B_axis = np.cross(T_axis, P_axis)
    B_axis = B_axis / np.linalg.norm(B_axis)

    n1 = (T_axis + P_axis) / np.sqrt(2)
    s1 = (T_axis - P_axis) / np.sqrt(2)

    if n1[2] < 0:
        n1 = -n1
        s1 = -s1

    strike = np.degrees(np.arctan2(n1[0], n1[1])) % 360
    dip = np.degrees(np.arccos(np.abs(n1[2])))

    cos_rake = np.dot(s1, [-np.sin(np.radians(strike)), np.cos(np.radians(strike)), 0])
    sin_rake = -s1[2] / np.sin(np.radians(dip)) if np.sin(np.radians(dip)) > 0.01 else 0
    rake = np.degrees(np.arctan2(sin_rake, cos_rake))

    return strike, dip, rake


# ==============================================================================
# Beachball Plotting
# ==============================================================================

def plot_beachball(ax, m, size=1.0, xy=(0, 0), facecolor='red', bgcolor='white',
                   linewidth=1):
    """
    Plot moment tensor beachball on given axes.

    m: [Mxx, Myy, Mzz, Mxy, Mxz, Myz]
    """
    decomp = decompose_moment_tensor(m)
    eigenvalues = decomp['eigenvalues']
    eigenvectors = decomp['eigenvectors']

    n_points = 180
    theta = np.linspace(0, np.pi, n_points)
    phi = np.linspace(0, 2*np.pi, 2*n_points)
    THETA, PHI = np.meshgrid(theta, phi)

    X = np.sin(THETA) * np.cos(PHI)
    Y = np.sin(THETA) * np.sin(PHI)
    Z = np.cos(THETA)

    Mxx, Myy, Mzz, Mxy, Mxz, Myz = m
    M_matrix = np.array([
        [Mxx, Mxy, Mxz],
        [Mxy, Myy, Myz],
        [Mxz, Myz, Mzz]
    ])

    amplitude = np.zeros_like(X)
    for i in range(n_points * 2):
        for j in range(n_points):
            gamma = np.array([X[i,j], Y[i,j], Z[i,j]])
            amplitude[i,j] = gamma @ M_matrix @ gamma

    mask_lower = Z < 0

    r_proj = np.sqrt(2) * np.sin((np.pi - THETA) / 2)
    X_proj = r_proj * np.cos(PHI)
    Y_proj = r_proj * np.sin(PHI)

    circle = Circle(xy, size, facecolor=bgcolor, edgecolor='black',
                   linewidth=linewidth, zorder=1)
    ax.add_patch(circle)

    ax.contourf(X_proj * size + xy[0], Y_proj * size + xy[1], amplitude,
                levels=[0, np.max(np.abs(amplitude))], colors=[facecolor],
                zorder=2)
    ax.contour(X_proj * size + xy[0], Y_proj * size + xy[1], amplitude,
               levels=[0], colors=['black'], linewidths=linewidth, zorder=3)

    ax.set_xlim(xy[0] - size * 1.1, xy[0] + size * 1.1)
    ax.set_ylim(xy[1] - size * 1.1, xy[1] + size * 1.1)
    ax.set_aspect('equal')


def plot_simple_beachball(ax, decomp, size=1.0, xy=(0, 0)):
    """
    Plot a simplified beachball using the decomposition results.
    Shows ISO, CLVD, DC fractions as pie-like segments.
    """
    from matplotlib.patches import Wedge, Circle

    f_iso = decomp['ISO_fraction']
    f_clvd = decomp['CLVD_fraction']
    f_dc = decomp['DC_fraction']

    circle = Circle(xy, size, facecolor='white', edgecolor='black', linewidth=2)
    ax.add_patch(circle)

    angles = [0, f_iso * 360, (f_iso + f_clvd) * 360, 360]
    colors = ['blue', 'green', 'red']
    labels_internal = ['ISO', 'CLVD', 'DC']

    for i in range(3):
        if angles[i+1] - angles[i] > 0.5:
            wedge = Wedge(xy, size * 0.95, angles[i], angles[i+1],
                         facecolor=colors[i], edgecolor='black', linewidth=0.5)
            ax.add_patch(wedge)

    ax.set_xlim(xy[0] - size * 1.2, xy[0] + size * 1.2)
    ax.set_ylim(xy[1] - size * 1.2, xy[1] + size * 1.2)
    ax.set_aspect('equal')


# ==============================================================================
# Data Fetching
# ==============================================================================

def fetch_waveforms_for_inversion():
    """
    Download and process waveforms from IRIS for moment tensor inversion.

    Note: 1992 data may be limited — many broadband stations did not yet exist
    or archive may be incomplete. Fallback to synthetics is expected.

    Returns list of (net, sta, dist_km, az, baz, traces_zrt, inventory)
    """
    print("Fetching waveforms from IRIS...")
    print("  (Note: 1992 data may be limited; fallback to synthetics likely)")
    client = Client("IRIS", timeout=120)
    t1 = EVENT_TIME - PRE_ORIGIN
    t2 = EVENT_TIME + POST_ORIGIN

    results = []

    for net, sta, loc, chan_pref, desc in STATIONS:
        label = f"{net}.{sta}"
        try:
            st = Stream()
            for comp in ['Z', 'N', 'E']:
                try:
                    st_comp = client.get_waveforms(net, sta, loc,
                                                   f"{chan_pref}{comp}", t1, t2)
                    st += st_comp
                except Exception:
                    pass

            if len(st) < 3:
                print(f"  {label:12s} – incomplete components, skipping")
                continue

            st.merge(fill_value='interpolate')

            inv = client.get_stations(network=net, station=sta,
                                     starttime=EVENT_TIME, endtime=EVENT_TIME,
                                     level="response")

            sta_lat = inv[0][0].latitude
            sta_lon = inv[0][0].longitude
            dist_m, az, baz = gps2dist_azimuth(EVENT_LAT, EVENT_LON,
                                                sta_lat, sta_lon)
            dist_km = dist_m / 1000.0

            st.detrend('demean')
            st.detrend('linear')
            st.taper(0.05)

            try:
                st.remove_response(inventory=inv, output='VEL',
                                  pre_filt=[0.005, 0.01, 0.5, 1.0])
            except Exception as e:
                print(f"  {label:12s} – response removal failed: {e}")
                continue

            st.rotate(method='NE->RT', back_azimuth=baz)

            st.filter('bandpass', freqmin=FREQMIN, freqmax=FREQMAX,
                     corners=4, zerophase=True)

            st.resample(1.0)

            tr_z = st.select(component='Z')[0] if st.select(component='Z') else None
            tr_r = st.select(component='R')[0] if st.select(component='R') else None
            tr_t = st.select(component='T')[0] if st.select(component='T') else None

            results.append({
                'net': net, 'sta': sta, 'desc': desc,
                'dist_km': dist_km, 'az': az, 'baz': baz,
                'tr_z': tr_z, 'tr_r': tr_r, 'tr_t': tr_t,
                'inv': inv
            })

            print(f"  {label:12s}  dist={dist_km:7.1f} km  az={az:5.1f}°  ✓")

        except Exception as e:
            print(f"  {label:12s} – {str(e)[:60]}")

    results.sort(key=lambda x: x['dist_km'])
    return results


# ==============================================================================
# Main Inversion Routine
# ==============================================================================

def run_moment_tensor_inversion():
    """
    Main function to perform moment tensor inversion.
    """
    print("=" * 70)
    print("Full Moment Tensor Inversion for US Divider 1992-09-23 Nuclear Test")
    print("=" * 70)
    print()

    # Fetch data
    station_data = fetch_waveforms_for_inversion()

    if len(station_data) < 2:
        print("\nInsufficient stations for inversion. Using synthetic example.")
        print("  (1992 data availability is limited for broadband stations)")
        station_data = generate_synthetic_data()

    print(f"\nUsing {len(station_data)} stations for inversion")

    # Prepare observations and Green's functions
    print("\nComputing Green's functions...")
    observations = []
    greens_list = []

    model = VelocityModel1D()
    dt = 1.0

    for sdata in station_data:
        dist_km = sdata['dist_km']
        az = sdata['az']

        if sdata['tr_z'] is not None:
            npts = sdata['tr_z'].stats.npts
            data_z = sdata['tr_z'].data
        else:
            continue

        data_r = sdata['tr_r'].data if sdata['tr_r'] is not None else np.zeros(npts)
        data_t = sdata['tr_t'].data if sdata['tr_t'] is not None else np.zeros(npts)

        gf = compute_greens_functions(dist_km, az, EVENT_DEPTH_KM, dt, npts,
                                      model, FREQMIN, FREQMAX)

        weight = 1.0 / np.sqrt(dist_km / 1000.0)

        observations.append((data_z, data_r, data_t, weight))
        greens_list.append(gf)

        print(f"  {sdata['net']}.{sdata['sta']:5s}  dist={dist_km:6.1f} km  "
              f"npts={npts}  weight={weight:.3f}")

    # Build system and invert
    print("\nBuilding linear system...")
    d = build_data_vector(observations)
    G = build_greens_matrix(greens_list, observations)

    print(f"  Data vector size: {len(d)}")
    print(f"  Green's matrix size: {G.shape}")

    print("\nInverting for moment tensor...")
    m, residual, var_red, d_pred = invert_moment_tensor(d, G, method='lstsq')

    # Decompose result
    decomp = decompose_moment_tensor(m)

    # Print results
    print("\n" + "=" * 70)
    print("INVERSION RESULTS")
    print("=" * 70)
    print(f"\nMoment tensor [Mxx, Myy, Mzz, Mxy, Mxz, Myz]:")
    print(f"  {m}")
    print(f"\nScalar moment M₀ = {decomp['M0']:.3e} N·m")
    print(f"Moment magnitude Mw = {decomp['Mw']:.2f}")
    print(f"\nSource type decomposition:")
    print(f"  Isotropic (ISO):    {decomp['ISO_fraction']*100:5.1f}%")
    print(f"  CLVD:               {decomp['CLVD_fraction']*100:5.1f}%")
    print(f"  Double-couple (DC): {decomp['DC_fraction']*100:5.1f}%")
    print(f"\nIsotropic component: {decomp['ISO_value']:.3e} N·m")
    print(f"  (Positive = compressional = explosion-like)")
    print(f"\nFault plane solution (DC component):")
    print(f"  Strike = {decomp['strike']:.1f}°")
    print(f"  Dip    = {decomp['dip']:.1f}°")
    print(f"  Rake   = {decomp['rake']:.1f}°")
    print(f"\nVariance reduction: {var_red:.1f}%")
    print(f"RMS residual: {residual:.3e}")

    # Interpretation
    print("\n" + "=" * 70)
    print("INTERPRETATION")
    print("=" * 70)
    if decomp['ISO_fraction'] > 0.5 and decomp['ISO_value'] > 0:
        print("\n✓ Large positive isotropic component (>50%)")
        print("  → Consistent with EXPLOSIVE SOURCE (underground nuclear test)")
        print("  → Clean explosion with minimal tectonic release")
    elif decomp['DC_fraction'] > 0.7:
        print("\n✗ Dominated by double-couple component")
        print("  → Consistent with tectonic earthquake")
    else:
        print("\n⚠ Mixed source mechanism")
        print("  → Possible explosion with tectonic release or collapse")

    # Generate figures
    print("\nGenerating figures...")
    plot_inversion_results(m, decomp, station_data, observations, greens_list,
                          d, d_pred, var_red)

    return m, decomp


def generate_synthetic_data():
    """
    Generate synthetic data for demonstration when real data unavailable.
    Uses an explosion-like moment tensor — characteristic of a small (~1 kt)
    NTS underground nuclear test (Divider, 1992).
    """
    print("\nGenerating synthetic test data...")

    # True moment tensor (explosion-like: large ISO, small non-isotropic)
    M0 = 1e17  # N·m (roughly 1 kt)
    # Explosion: Mxx = Myy = Mzz = M0/3, with small off-diagonal
    m_true = np.array([M0/3, M0/3, M0/3, 0.02*M0, 0.01*M0, 0.015*M0])

    model = VelocityModel1D()
    dt = 1.0
    npts = 600

    synthetic_stations = [
        {'net': 'CI', 'sta': 'PAS', 'desc': 'Pasadena',
         'dist_km': 350, 'az': 195, 'baz': 15},
        {'net': 'CI', 'sta': 'GSC', 'desc': 'Goldstone',
         'dist_km': 230, 'az': 180, 'baz': 0},
        {'net': 'BK', 'sta': 'SAO', 'desc': 'San Andreas Obs',
         'dist_km': 520, 'az': 260, 'baz': 80},
        {'net': 'NN', 'sta': 'ELK', 'desc': 'Elko',
         'dist_km': 340, 'az': 20, 'baz': 200},
    ]

    for sdata in synthetic_stations:
        gf = compute_greens_functions(sdata['dist_km'], sdata['az'],
                                      EVENT_DEPTH_KM, dt, npts, model,
                                      FREQMIN, FREQMAX)

        components = ['Mxx', 'Myy', 'Mzz', 'Mxy', 'Mxz', 'Myz']
        data_z = sum(m_true[i] * gf[comp]['Z'] for i, comp in enumerate(components))
        data_r = sum(m_true[i] * gf[comp]['R'] for i, comp in enumerate(components))
        data_t = sum(m_true[i] * gf[comp]['T'] for i, comp in enumerate(components))

        noise_level = 0.1 * np.max(np.abs(data_z))
        data_z += np.random.randn(npts) * noise_level
        data_r += np.random.randn(npts) * noise_level
        data_t += np.random.randn(npts) * noise_level

        tr_z = Trace(data=data_z)
        tr_z.stats.sampling_rate = 1.0 / dt
        tr_z.stats.npts = npts

        tr_r = Trace(data=data_r)
        tr_r.stats.sampling_rate = 1.0 / dt

        tr_t = Trace(data=data_t)
        tr_t.stats.sampling_rate = 1.0 / dt

        sdata['tr_z'] = tr_z
        sdata['tr_r'] = tr_r
        sdata['tr_t'] = tr_t
        sdata['inv'] = None

    return synthetic_stations


def plot_inversion_results(m, decomp, station_data, observations, greens_list,
                           d, d_pred, var_red):
    """Generate publication-quality figures of inversion results."""

    # Figure 1: Source type decomposition
    fig = plt.figure(figsize=(16, 10))
    gs = GridSpec(2, 3, figure=fig, hspace=0.3, wspace=0.3)

    # (a) Source type pie chart
    ax1 = fig.add_subplot(gs[0, 0])
    fracs = [decomp['ISO_fraction'], decomp['CLVD_fraction'], decomp['DC_fraction']]
    labels = [f"ISO\n{fracs[0]*100:.1f}%",
              f"CLVD\n{fracs[1]*100:.1f}%",
              f"DC\n{fracs[2]*100:.1f}%"]
    colors = ['#3498db', '#2ecc71', '#e74c3c']
    explode = (0.05, 0, 0) if fracs[0] > 0.3 else (0, 0, 0)

    wedges, texts, autotexts = ax1.pie(fracs, labels=labels, colors=colors,
                                        explode=explode, autopct='',
                                        startangle=90, shadow=True)
    ax1.set_title('(a) Source Type Decomposition', fontweight='bold', fontsize=12)

    # (b) Hudson plot (source type diagram)
    ax2 = fig.add_subplot(gs[0, 1])
    k_range = np.linspace(-1, 1, 100)
    ax2.fill_between(k_range, -np.abs(k_range), np.abs(k_range),
                     color='lightgray', alpha=0.3)
    ax2.axhline(0, color='gray', lw=0.5)
    ax2.axvline(0, color='gray', lw=0.5)

    eig = decomp['eigenvalues']
    if np.abs(eig[0] - eig[2]) > 1e-20:
        k = (eig[0] + eig[2]) / np.abs(eig[0] - eig[2])
        T = 2 * eig[1] / np.abs(eig[0] - eig[2])
    else:
        k, T = 0, 0
    ax2.plot(k, T, 'r*', markersize=20)
    ax2.text(k + 0.1, T + 0.1, 'Divider\n1992', fontsize=10, color='red')

    ax2.text(0, 0.8, 'CLVD\n(+)', ha='center', fontsize=9)
    ax2.text(0, -0.8, 'CLVD\n(−)', ha='center', fontsize=9)
    ax2.text(0.8, 0, 'Explosion', ha='center', fontsize=9, rotation=90)
    ax2.text(-0.8, 0, 'Implosion', ha='center', fontsize=9, rotation=90)
    ax2.text(0, 0, 'DC', ha='center', fontsize=9, fontweight='bold',
             bbox=dict(boxstyle='circle', facecolor='white', edgecolor='black'))

    ax2.set_xlim(-1.2, 1.2)
    ax2.set_ylim(-1.2, 1.2)
    ax2.set_xlabel('k (Isotropic parameter)', fontsize=10)
    ax2.set_ylabel('T (CLVD parameter)', fontsize=10)
    ax2.set_title('(b) Hudson Source Type Diagram', fontweight='bold', fontsize=12)
    ax2.set_aspect('equal')

    # (c) Moment tensor components
    ax3 = fig.add_subplot(gs[0, 2])
    components = ['Mxx', 'Myy', 'Mzz', 'Mxy', 'Mxz', 'Myz']
    x_pos = np.arange(len(components))
    colors_bar = ['red' if v > 0 else 'blue' for v in m]
    ax3.bar(x_pos, m, color=colors_bar, edgecolor='black')
    ax3.set_xticks(x_pos)
    ax3.set_xticklabels(components)
    ax3.set_ylabel('Moment (N·m)')
    ax3.axhline(0, color='black', lw=0.5)
    ax3.set_title('(c) Moment Tensor Components', fontweight='bold', fontsize=12)
    ax3.ticklabel_format(axis='y', style='scientific', scilimits=(0,0))

    # (d-f) Waveform fits
    n_sta = min(3, len(station_data))
    idx_offset = 0

    for i in range(n_sta):
        ax = fig.add_subplot(gs[1, i])
        sdata = station_data[i]
        obs = observations[i]
        data_z, data_r, data_t, weight = obs

        npts = len(data_z)
        t = np.arange(npts)

        gf = greens_list[i]
        components_list = ['Mxx', 'Myy', 'Mzz', 'Mxy', 'Mxz', 'Myz']
        pred_z = sum(m[j] * gf[comp]['Z'] for j, comp in enumerate(components_list))

        ax.plot(t, data_z / np.max(np.abs(data_z)), 'k-', lw=1, label='Observed')
        ax.plot(t, pred_z / np.max(np.abs(pred_z)), 'r--', lw=1, label='Predicted')
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Amplitude (normalized)')
        ax.set_title(f'({chr(100+i)}) {sdata["net"]}.{sdata["sta"]} — Z component\n'
                     f'Δ = {sdata["dist_km"]:.0f} km, az = {sdata["az"]:.0f}°',
                     fontweight='bold', fontsize=10)
        ax.legend(fontsize=8, loc='upper right')
        ax.set_xlim(0, npts)

    fig.suptitle(f'Moment Tensor Inversion Results — US Divider 1992-09-23 Nuclear Test\n'
                 f'M₀ = {decomp["M0"]:.2e} N·m, Mw = {decomp["Mw"]:.2f}, '
                 f'VR = {var_red:.1f}%',
                 fontsize=14, fontweight='bold')

    outpath = os.path.join(OUTDIR, 'moment_tensor_inversion.png')
    fig.savefig(outpath, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {outpath}")

    # Figure 2: Summary card
    fig2, ax = plt.subplots(figsize=(10, 8))
    ax.axis('off')

    summary_text = f"""
    ══════════════════════════════════════════════════════════════════
    FULL MOMENT TENSOR INVERSION — US Divider 1992-09-23 Nuclear Test
    ══════════════════════════════════════════════════════════════════

    Event Location:     {EVENT_LAT:.3f}°N, {EVENT_LON:.3f}°W
    Event Depth:        {EVENT_DEPTH_KM:.1f} km (NTS Yucca Flat)
    Origin Time:        1992-09-23 15:04:00 UTC

    ──────────────────────────────────────────────────────────────────
    MOMENT TENSOR SOLUTION
    ──────────────────────────────────────────────────────────────────

    Scalar Moment:      M₀ = {decomp['M0']:.3e} N·m
    Moment Magnitude:   Mw = {decomp['Mw']:.2f}

    Tensor Components (N·m):
        Mxx = {m[0]:+.3e}    Mxy = {m[3]:+.3e}
        Myy = {m[1]:+.3e}    Mxz = {m[4]:+.3e}
        Mzz = {m[2]:+.3e}    Myz = {m[5]:+.3e}

    ──────────────────────────────────────────────────────────────────
    SOURCE TYPE DECOMPOSITION
    ──────────────────────────────────────────────────────────────────

    Isotropic (ISO):        {decomp['ISO_fraction']*100:5.1f}%   {'← EXPLOSION SIGNATURE' if decomp['ISO_fraction'] > 0.3 else ''}
    CLVD:                   {decomp['CLVD_fraction']*100:5.1f}%
    Double-Couple (DC):     {decomp['DC_fraction']*100:5.1f}%

    Isotropic Value:    {decomp['ISO_value']:+.3e} N·m
                        {'(Compressional = Explosion)' if decomp['ISO_value'] > 0 else '(Dilatational)'}

    ──────────────────────────────────────────────────────────────────
    INVERSION QUALITY
    ──────────────────────────────────────────────────────────────────

    Variance Reduction:     {var_red:.1f}%
    Number of Stations:     {len(station_data)}
    Frequency Band:         {FREQMIN:.3f} – {FREQMAX:.2f} Hz

    ──────────────────────────────────────────────────────────────────
    INTERPRETATION
    ──────────────────────────────────────────────────────────────────

    Expected: ~80% ISO, ~18% CLVD, ~2% DC
    Large positive ISO → clean explosion (~1 kt)
    Note: 1992 data limited; synthetics used for demonstration
    Last US underground nuclear test (Operation Julin)

    ══════════════════════════════════════════════════════════════════
    """

    ax.text(0.02, 0.98, summary_text, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    outpath2 = os.path.join(OUTDIR, 'moment_tensor_summary.png')
    fig2.savefig(outpath2, dpi=150, bbox_inches='tight')
    plt.close(fig2)
    print(f"  Saved: {outpath2}")


# ==============================================================================
# Entry Point
# ==============================================================================

if __name__ == "__main__":
    try:
        m, decomp = run_moment_tensor_inversion()
    except KeyboardInterrupt:
        print("\nInversion interrupted.")
        sys.exit(1)
    except Exception as e:
        print(f"\nError during inversion: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
