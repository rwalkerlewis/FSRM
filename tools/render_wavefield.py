#!/usr/bin/env python3
"""Render a polished wavefield figure pack from FSRM solution.h5.

Generic over any underground-explosion example. Reads grid extents and source
location from solution.h5 + a config file (INI sections [EXPLOSION_SOURCE] and
[SEISMOMETER_*]) and produces the same six-panel + animated figure set as the
Sedan reference pack.

Usage:
    python3 render_wavefield.py <output_dir> <config_path> [<event_title>]

Outputs under <output_dir>/figures/:
    wavefield_xz_slice.gif        animation of |u| on the y=mid plane
    wavefield_surface.gif         animation of vertical surface displacement
    wavefield_panels.png          6-panel XZ slice snapshots
    surface_panels.png            6-panel surface vertical-displacement snapshots
    seismograms.png               N-station x 3-comp seismogram traces
"""
import sys, os, struct, glob, re, configparser
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation


def load_sac(path):
    with open(path, 'rb') as fh:
        h = fh.read(632)
        delta = struct.unpack('<f', h[0:4])[0]
        npts = struct.unpack('<i', h[316:320])[0]
        data = np.array(struct.unpack(f'<{npts}f', fh.read(4 * npts)))
    return delta, data


def parse_config(cfg_path):
    """Return (source_xyz_km, stations[list of (x_km,y_km,z_km,name)], dt_step)."""
    cp = configparser.ConfigParser(strict=False, inline_comment_prefixes=('#',))
    cp.optionxform = str
    cp.read(cfg_path)
    src = None
    if cp.has_section('EXPLOSION_SOURCE'):
        try:
            sx = float(cp['EXPLOSION_SOURCE']['location_x'])
            sy = float(cp['EXPLOSION_SOURCE']['location_y'])
            sz = float(cp['EXPLOSION_SOURCE']['location_z'])
            src = (sx / 1000.0, sy / 1000.0, sz / 1000.0)
        except (KeyError, ValueError):
            src = None
    stations = []
    for sec in cp.sections():
        if not re.match(r'^SEISMOMETER_\d+$', sec):
            continue
        try:
            name = cp[sec].get('sta', sec).strip()
            coords = cp[sec]['location_xyz'].split(',')
            x, y, z = (float(c) for c in coords[:3])
            stations.append((x / 1000.0, y / 1000.0, z / 1000.0, name))
        except (KeyError, ValueError):
            continue
    dt_step = 0.001
    if cp.has_section('SIMULATION'):
        try:
            dt_step = float(cp['SIMULATION'].get('dt_initial', '0.001'))
        except ValueError:
            pass
    return src, stations, dt_step


def render(out_dir: str, cfg_path: str, title: str = '') -> None:
    fig_dir = os.path.join(out_dir, 'figures')
    os.makedirs(fig_dir, exist_ok=True)
    src, stations, dt_step = parse_config(cfg_path)

    with h5py.File(os.path.join(out_dir, 'solution.h5'), 'r') as h:
        v = h['geometry/vertices'][:]
        u = h['fields/solution'][:]
    T, Nv, _ = u.shape

    pk = np.abs(u).max(axis=(1, 2))
    nz = np.where(pk > 0)[0]
    if len(nz) == 0:
        print('no nonzero snapshots; aborting')
        return
    print(f'{len(nz)} nonzero snapshots, peak |u|={pk.max():.3e} m')

    x_km_min, x_km_max = v[:, 0].min() / 1000.0, v[:, 0].max() / 1000.0
    y_km_min, y_km_max = v[:, 1].min() / 1000.0, v[:, 1].max() / 1000.0
    z_km_min, z_km_max = v[:, 2].min() / 1000.0, v[:, 2].max() / 1000.0
    y_mid_m = 0.5 * (v[:, 1].min() + v[:, 1].max())
    if src is None:
        src = (0.5 * (x_km_min + x_km_max),
               y_mid_m / 1000.0,
               z_km_max)

    # ----- XZ slice through the source y -----
    src_y_m = src[1] * 1000.0
    slab = np.abs(v[:, 1] - src_y_m) < 1500.0
    if slab.sum() < 4:
        slab = np.abs(v[:, 1] - y_mid_m) < 1500.0
    tri = Triangulation(v[slab, 0] / 1000.0, v[slab, 2] / 1000.0)
    umag_slab = np.linalg.norm(u[:, slab, :], axis=2)
    vmax_xz = float(umag_slab[nz].max())
    if vmax_xz <= 0:
        vmax_xz = 1.0e-9
    levels_xz = np.linspace(0.0, vmax_xz, 30)

    station_xz = [(sx, sz, name) for (sx, sy, sz, name) in stations
                  if abs(sy * 1000.0 - src_y_m) < 1500.0]

    fig, ax = plt.subplots(figsize=(11, 4))
    sm_xz = plt.cm.ScalarMappable(
        cmap='inferno', norm=plt.Normalize(vmin=0.0, vmax=vmax_xz))
    cbar_xz = fig.colorbar(sm_xz, ax=ax, extend='max', pad=0.02, fraction=0.046)
    cbar_xz.set_label('|u| [m]')
    frame_paths = []
    for k, i in enumerate(nz):
        ax.clear()
        ax.tricontourf(tri, umag_slab[i], levels=levels_xz, cmap='inferno',
                       extend='max', vmin=0.0, vmax=vmax_xz)
        ax.plot(src[0], src[2], 'c*', ms=18, mec='white')
        for sx, sz, name in station_xz:
            ax.plot(sx, sz, 'wv', ms=10)
            ax.annotate(name, (sx, sz), xytext=(0, 8),
                        textcoords='offset points',
                        color='white', ha='center', fontsize=8)
        ax.set_xlim(x_km_min, x_km_max)
        ax.set_ylim(z_km_min, z_km_max)
        ax.set_aspect('equal')
        ax.set_xlabel('x [km]'); ax.set_ylabel('z [km]')
        ax.set_title(f'{title} |u(x,z)| at y={src[1]:.1f} km   '
                     f't={i*dt_step:.3f} s   '
                     f'peak={umag_slab[i].max():.3f} m')
        fp = os.path.join(fig_dir, f'xz_{k:03d}.png')
        plt.savefig(fp, dpi=110, bbox_inches='tight')
        frame_paths.append(fp)
    plt.close(fig)

    try:
        from PIL import Image
        imgs = [Image.open(p) for p in frame_paths]
        imgs[0].save(os.path.join(fig_dir, 'wavefield_xz_slice.gif'),
                     save_all=True, append_images=imgs[1:], duration=80, loop=0)
    except Exception as e:
        print('GIF skipped:', e)

    sample = [nz[1] if len(nz) > 1 else nz[0],
              nz[max(1, len(nz)//6)],
              nz[max(1, len(nz)//3)],
              nz[max(1, len(nz)//2)],
              nz[max(1, 3*len(nz)//4)],
              nz[-1]]
    fig, axs = plt.subplots(2, 3, figsize=(16, 6))
    for ax, i in zip(axs.flat, sample):
        ax.tricontourf(tri, umag_slab[i], levels=levels_xz, cmap='inferno',
                       extend='max', vmin=0.0, vmax=vmax_xz)
        ax.plot(src[0], src[2], 'c*', ms=10, mec='white')
        ax.set_xlim(x_km_min, x_km_max)
        ax.set_ylim(z_km_min, z_km_max)
        ax.set_aspect('equal')
        ax.set_xlabel('x [km]'); ax.set_ylabel('z [km]')
        ax.set_title(f't={i*dt_step:.2f}s   peak={umag_slab[i].max():.2f}m')
    fig.suptitle(f'{title} wavefield slice |u| at y={src[1]:.1f} km')
    fig.tight_layout(rect=[0, 0, 0.94, 1.0])
    cax = fig.add_axes([0.955, 0.12, 0.012, 0.76])
    sm_p = plt.cm.ScalarMappable(
        cmap='inferno', norm=plt.Normalize(vmin=0.0, vmax=vmax_xz))
    cb = fig.colorbar(sm_p, cax=cax, extend='max')
    cb.set_label('|u| [m]')
    plt.savefig(os.path.join(fig_dir, 'wavefield_panels.png'), dpi=120,
                bbox_inches='tight')
    plt.close(fig)
    print('XZ slice figures done')

    # ----- Surface vertical displacement at z = max -----
    z_top_m = v[:, 2].max()
    surf = v[:, 2] > z_top_m - 1.0
    if surf.sum() < 4:
        surf = v[:, 2] > z_top_m - 100.0
    tri_s = Triangulation(v[surf, 0] / 1000.0, v[surf, 1] / 1000.0)
    uz_surf = u[:, surf, 2]
    vmax_s = float(max(abs(uz_surf[nz].max()), abs(uz_surf[nz].min())))
    if vmax_s <= 0:
        vmax_s = 1.0e-9
    levels_s = np.linspace(-vmax_s, vmax_s, 31)

    fig, ax = plt.subplots(figsize=(8, 6))
    sm_s = plt.cm.ScalarMappable(
        cmap='RdBu_r', norm=plt.Normalize(vmin=-vmax_s, vmax=vmax_s))
    cbar_s = fig.colorbar(sm_s, ax=ax, extend='both', pad=0.02, fraction=0.046)
    cbar_s.set_label('u_z [m]')
    frame_paths = []
    for k, i in enumerate(nz):
        ax.clear()
        ax.tricontourf(tri_s, uz_surf[i], levels=levels_s, cmap='RdBu_r',
                       extend='both', vmin=-vmax_s, vmax=vmax_s)
        ax.plot(src[0], src[1], 'k*', ms=14)
        for sx, sy, sz, name in stations:
            ax.plot(sx, sy, 'kv', ms=8)
            ax.annotate(name, (sx, sy), xytext=(0, 6),
                        textcoords='offset points', ha='center', fontsize=8)
        ax.set_xlim(x_km_min, x_km_max)
        ax.set_ylim(y_km_min, y_km_max)
        ax.set_aspect('equal')
        ax.set_xlabel('x [km]'); ax.set_ylabel('y [km]')
        ax.set_title(f'{title}  u_z surface  t={i*dt_step:.3f}s  '
                     f'pk={np.abs(uz_surf[i]).max():.2f}m')
        fp = os.path.join(fig_dir, f'surf_{k:03d}.png')
        plt.savefig(fp, dpi=100, bbox_inches='tight')
        frame_paths.append(fp)
    plt.close(fig)
    try:
        from PIL import Image
        imgs = [Image.open(p) for p in frame_paths]
        imgs[0].save(os.path.join(fig_dir, 'wavefield_surface.gif'),
                     save_all=True, append_images=imgs[1:], duration=80, loop=0)
    except Exception as e:
        print('surface GIF skipped:', e)

    fig, axs = plt.subplots(2, 3, figsize=(16, 9))
    for ax, i in zip(axs.flat, sample):
        ax.tricontourf(tri_s, uz_surf[i], levels=levels_s, cmap='RdBu_r',
                       extend='both', vmin=-vmax_s, vmax=vmax_s)
        ax.plot(src[0], src[1], 'k*', ms=10)
        ax.set_xlim(x_km_min, x_km_max)
        ax.set_ylim(y_km_min, y_km_max)
        ax.set_aspect('equal')
        ax.set_xlabel('x [km]'); ax.set_ylabel('y [km]')
        ax.set_title(f't={i*dt_step:.2f}s  pk={np.abs(uz_surf[i]).max():.2f}m')
    fig.suptitle(f'{title} surface vertical displacement u_z(x,y,z=top)')
    fig.tight_layout(rect=[0, 0, 0.94, 1.0])
    cax = fig.add_axes([0.955, 0.10, 0.012, 0.80])
    sm_sp = plt.cm.ScalarMappable(
        cmap='RdBu_r', norm=plt.Normalize(vmin=-vmax_s, vmax=vmax_s))
    cb = fig.colorbar(sm_sp, cax=cax, extend='both')
    cb.set_label('u_z [m]')
    plt.savefig(os.path.join(fig_dir, 'surface_panels.png'), dpi=120,
                bbox_inches='tight')
    plt.close(fig)
    print('surface figures done')

    # ----- Seismograms (any SAC files in this dir tree) -----
    sacs = sorted(glob.glob(os.path.join(out_dir, '**', '*.sac'),
                            recursive=True))
    if not sacs:
        return
    by_sta = {}
    for p in sacs:
        m = re.match(r'^[A-Z]{1,2}\.([^.]+)\.\d+\.(BH[ZNE])\.sac$',
                     os.path.basename(p))
        if not m:
            continue
        by_sta.setdefault(m.group(1), {})[m.group(2)] = p
    if not by_sta:
        return
    stations_sorted = sorted(by_sta.keys())
    comps = ['BHZ', 'BHN', 'BHE']
    fig, axs = plt.subplots(len(stations_sorted), 3,
                            figsize=(14, 2.0 * len(stations_sorted)),
                            sharex=True, squeeze=False)
    for i, sta in enumerate(stations_sorted):
        for j, c in enumerate(comps):
            f = by_sta[sta].get(c)
            if not f:
                axs[i][j].set_axis_off()
                continue
            delta, data = load_sac(f)
            t = np.arange(len(data)) * delta
            axs[i][j].plot(t, data, 'k-', lw=0.7)
            pk = float(np.abs(data).max())
            axs[i][j].set_title(f'{sta} {c}  pk={pk:.2e}', fontsize=8)
            axs[i][j].grid(alpha=0.3)
            if i == len(stations_sorted) - 1:
                axs[i][j].set_xlabel('t [s]')
    fig.suptitle(f'{title} seismograms')
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, 'seismograms.png'), dpi=120,
                bbox_inches='tight')
    plt.close(fig)
    print('seismograms figure done')


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)
    title = sys.argv[3] if len(sys.argv) > 3 else ''
    render(sys.argv[1], sys.argv[2], title)
