#!/usr/bin/env python3
"""Render a polished Sedan wavefield figure pack from solution.h5.

Outputs under <output_dir>/figures/:
  wavefield_xz_slice.gif        animation of |u| on the y=15 km plane
  wavefield_surface.gif         animation of vertical surface displacement
  wavefield_panels.png          6-panel XZ slice snapshots
  surface_panels.png            6-panel surface vertical-displacement snapshots
  velocity_seismograms.png      3-station x 3-comp velocity traces (m/s)
  displacement_seismograms.png  3-station x 3-comp displacement traces (m)
"""
import sys, os, struct, glob
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


def render(out_dir: str, dt_step: float = 0.001) -> None:
    fig_dir = os.path.join(out_dir, 'figures')
    os.makedirs(fig_dir, exist_ok=True)

    with h5py.File(os.path.join(out_dir, 'solution.h5'), 'r') as h:
        v = h['geometry/vertices'][:]
        u = h['fields/solution'][:]
    T, Nv, _ = u.shape

    pk = np.abs(u).max(axis=(1, 2))
    nz = np.where(pk > 0)[0]
    print(f'{len(nz)} nonzero snapshots, peak |u|={pk.max():.3e} m')

    # ----- XZ slice at y = 15 km -----
    slab = np.abs(v[:, 1] - 15000.0) < 1500.0
    tri = Triangulation(v[slab, 0] / 1000.0, v[slab, 2] / 1000.0)
    umag_slab = np.linalg.norm(u[:, slab, :], axis=2)
    vmax_xz = umag_slab[nz].max()
    levels_xz = np.linspace(0.0, vmax_xz, 30)

    fig, ax = plt.subplots(figsize=(11, 4))
    sm_xz = plt.cm.ScalarMappable(
        cmap='inferno', norm=plt.Normalize(vmin=0.0, vmax=vmax_xz))
    cbar_xz = fig.colorbar(sm_xz, ax=ax, extend='max', pad=0.02,
                           fraction=0.046)
    cbar_xz.set_label('|u| [m]')
    frame_paths = []
    station_xz = [(15, 5, 'SPALL'), (17, 5, 'STA02'), (20, 5, 'STA05'),
                  (25, 5, 'STA10')]
    for k, i in enumerate(nz):
        ax.clear()
        ax.tricontourf(tri, umag_slab[i], levels=levels_xz, cmap='inferno',
                       extend='max', vmin=0.0, vmax=vmax_xz)
        ax.plot(15.0, 4.806, 'c*', ms=18, mec='white')
        for sx, sz, name in station_xz:
            ax.plot(sx, sz, 'wv', ms=10)
            ax.annotate(name, (sx, sz), xytext=(0, 8),
                        textcoords='offset points',
                        color='white', ha='center', fontsize=8)
        ax.set_xlim(0, 30); ax.set_ylim(0, 5)
        ax.set_aspect('equal')
        ax.set_xlabel('x [km]'); ax.set_ylabel('z [km]')
        ax.set_title(f'Sedan 1962 |u(x,z)| at y=15 km   t={i*dt_step:.3f} s   '
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

    # 6-panel summary
    sample = [nz[1], nz[len(nz)//6], nz[len(nz)//3], nz[len(nz)//2],
              nz[3*len(nz)//4], nz[-1]]
    fig, axs = plt.subplots(2, 3, figsize=(16, 6))
    for ax, i in zip(axs.flat, sample):
        cf = ax.tricontourf(tri, umag_slab[i], levels=levels_xz,
                            cmap='inferno', extend='max',
                            vmin=0.0, vmax=vmax_xz)
        ax.plot(15.0, 4.806, 'c*', ms=10, mec='white')
        ax.set_xlim(0, 30); ax.set_ylim(0, 5); ax.set_aspect('equal')
        ax.set_xlabel('x [km]'); ax.set_ylabel('z [km]')
        ax.set_title(f't={i*dt_step:.2f}s   peak={umag_slab[i].max():.2f}m')
    fig.suptitle('Sedan 1962 wavefield slice |u| at y=15 km')
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

    # ----- Surface vertical displacement at z = 5 km -----
    surf = v[:, 2] > 4999.0
    tri_s = Triangulation(v[surf, 0] / 1000.0, v[surf, 1] / 1000.0)
    uz_surf = u[:, surf, 2]
    vmax_s = max(abs(uz_surf[nz].max()), abs(uz_surf[nz].min()))
    levels_s = np.linspace(-vmax_s, vmax_s, 31)

    fig, ax = plt.subplots(figsize=(8, 6))
    sm_s = plt.cm.ScalarMappable(
        cmap='RdBu_r', norm=plt.Normalize(vmin=-vmax_s, vmax=vmax_s))
    cbar_s = fig.colorbar(sm_s, ax=ax, extend='both', pad=0.02,
                          fraction=0.046)
    cbar_s.set_label('u_z [m]')
    frame_paths = []
    for k, i in enumerate(nz):
        ax.clear()
        ax.tricontourf(tri_s, uz_surf[i], levels=levels_s, cmap='RdBu_r',
                       extend='both', vmin=-vmax_s, vmax=vmax_s)
        ax.plot(15.0, 15.0, 'k*', ms=14)
        for x, y, name in [(15, 15, 'SPALL'), (17, 15, 'STA02'),
                           (20, 15, 'STA05'), (25, 15, 'STA10'),
                           (15, 20, 'N05')]:
            ax.plot(x, y, 'kv', ms=8)
            ax.annotate(name, (x, y), xytext=(0, 6),
                        textcoords='offset points', ha='center', fontsize=8)
        ax.set_xlim(0, 30); ax.set_ylim(0, 30); ax.set_aspect('equal')
        ax.set_xlabel('x [km]'); ax.set_ylabel('y [km]')
        ax.set_title(f'Sedan 1962  u_z surface  t={i*dt_step:.3f}s  '
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
        ax.plot(15.0, 15.0, 'k*', ms=10)
        ax.set_xlim(0, 30); ax.set_ylim(0, 30); ax.set_aspect('equal')
        ax.set_xlabel('x [km]'); ax.set_ylabel('y [km]')
        ax.set_title(f't={i*dt_step:.2f}s  pk={np.abs(uz_surf[i]).max():.2f}m')
    fig.suptitle('Sedan 1962 surface vertical displacement u_z(x,y,z=5km)')
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

    # ----- Seismogram panels (whatever SAC files exist alongside) -----
    sacs = glob.glob(os.path.join(out_dir, '*.sac'))
    if not sacs:
        return
    stations = sorted({os.path.basename(p).split('.')[1] for p in sacs})
    comps = ['BHZ', 'BHN', 'BHE']
    fig, axs = plt.subplots(len(stations), 3, figsize=(14, 2.0 * len(stations)),
                            sharex=True)
    if len(stations) == 1:
        axs = axs[None, :]
    for i, sta in enumerate(stations):
        for j, c in enumerate(comps):
            f = os.path.join(out_dir, f'XX.{sta}.00.{c}.sac')
            if not os.path.exists(f):
                continue
            delta, data = load_sac(f)
            t = np.arange(len(data)) * delta
            axs[i][j].plot(t, data, 'k-', lw=0.7)
            pk = float(np.abs(data).max())
            axs[i][j].set_title(f'{sta} {c}  pk={pk:.2e}', fontsize=8)
            axs[i][j].grid(alpha=0.3)
            if i == len(stations) - 1:
                axs[i][j].set_xlabel('t [s]')
    fig.suptitle('Sedan 1962 seismograms')
    plt.tight_layout()
    plt.savefig(os.path.join(fig_dir, 'seismograms.png'), dpi=120,
                bbox_inches='tight')
    plt.close(fig)
    print('seismograms figure done')


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(1)
    render(sys.argv[1])
