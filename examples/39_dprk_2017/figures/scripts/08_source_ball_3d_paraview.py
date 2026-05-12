#!/usr/bin/env python3
"""
Build source_ball_3d_vis.h5 + source_ball_3d_vis.xdmf from the
source_ball_3d_*.h5.csv snapshots produced by the DPRK 2017 coupled
RADIAL_LAGRANGIAN + THREE_DIMENSIONAL simulation.

Output (written to output_radial/):
  source_ball_3d_vis.h5    -- HDF5 containing mesh + all per-step scalar fields
  source_ball_3d_vis.xdmf  -- XDMF2 temporal collection; open THIS in ParaView

Fields available in ParaView:
  sigma_xx/yy/zz/xy/xz/yz  -- Cauchy stress components (Pa)
  eps_p_eq                   -- Equivalent plastic strain
  von_mises                  -- Von Mises stress (Pa)  [derived]
  mean_stress                -- Mean normal stress (Pa), positive = tensile [derived]
  rho                        -- Density (kg/m^3)
  e_int                      -- Internal energy (J/m^3)
  T_m                        -- Melt temperature (K)
  E_r                        -- Radiation energy density (J/m^3)
  source_forcing_power       -- Deposited power per cell (W)
  inner_cavity_marker        -- 1 inside cavity wall, 0 outside
  dp_pressure_field          -- Drucker-Prager pressure perturbation (Pa)
"""

import os
import sys
import glob
import re

import numpy as np
import h5py
import pandas as pd

# --- paths ---------------------------------------------------------------

SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
EXAMPLE_DIR = os.path.dirname(os.path.dirname(SCRIPT_DIR))   # examples/39_dprk_2017
OUTPUT_DIR  = os.path.join(EXAMPLE_DIR, "output_radial")
MESH_PREFIX = os.path.join(EXAMPLE_DIR, "../../cache/source_ball_meshes/dprk_2017")

OUT_H5      = os.path.join(OUTPUT_DIR, "source_ball_3d_vis.h5")
OUT_XDMF    = os.path.join(OUTPUT_DIR, "source_ball_3d_vis.xdmf")
H5_BASENAME = "source_ball_3d_vis.h5"   # relative path used inside XDMF DataItems

# -------------------------------------------------------------------------
# Load TetGen mesh
# -------------------------------------------------------------------------

print("Loading TetGen mesh ...", flush=True)

verts = []
with open(MESH_PREFIX + ".node") as f:
    n_nodes = int(f.readline().split()[0])
    for _ in range(n_nodes):
        row = f.readline().split()
        verts.append([float(row[1]), float(row[2]), float(row[3])])
verts = np.array(verts, dtype=np.float64)  # (N_verts, 3)  local-frame coords (m)

tets = []
with open(MESH_PREFIX + ".ele") as f:
    n_tets = int(f.readline().split()[0])
    for _ in range(n_tets):
        row = f.readline().split()
        # TetGen is 1-indexed; convert to 0-indexed for XDMF / VTK
        tets.append([int(row[1]) - 1, int(row[2]) - 1,
                     int(row[3]) - 1, int(row[4]) - 1])
tets = np.array(tets, dtype=np.int32)  # (N_tets, 4)

N_VERTS = len(verts)
N_TETS  = len(tets)
print(f"  {N_VERTS} vertices,  {N_TETS} tetrahedra", flush=True)

# -------------------------------------------------------------------------
# Find CSV snapshots
# -------------------------------------------------------------------------

csv_files = sorted(glob.glob(os.path.join(OUTPUT_DIR, "source_ball_3d_*.h5.csv")))
if not csv_files:
    print(f"ERROR: no source_ball_3d_*.h5.csv files found in {OUTPUT_DIR}", file=sys.stderr)
    sys.exit(1)
print(f"Found {len(csv_files)} snapshots", flush=True)

SCALAR_FIELDS = [
    'sigma_xx', 'sigma_yy', 'sigma_zz',
    'sigma_xy', 'sigma_xz', 'sigma_yz',
    'eps_p_eq',
    'rho', 'e_int', 'T_m', 'E_r',
    'source_forcing_power', 'inner_cavity_marker', 'dp_pressure_field',
]

# -------------------------------------------------------------------------
# Build HDF5 + XDMF
# -------------------------------------------------------------------------

print(f"Writing {OUT_H5} ...", flush=True)

times       = []
snap_labels = []

with h5py.File(OUT_H5, 'w') as hf:

    # Store mesh topology and geometry once (shared by all time steps)
    hf.create_dataset("Mesh/vertices",     data=verts, compression='gzip', shuffle=True)
    hf.create_dataset("Mesh/connectivity", data=tets,  compression='gzip', shuffle=True)

    for csv_path in csv_files:
        m = re.search(r'source_ball_3d_(\d+)\.h5\.csv$', csv_path)
        snap_idx = int(m.group(1))
        step     = snap_idx * 10            # cadence_steps = 10
        t        = step / 3000.0 * 3.0     # linear interpolation to t=3.0 at step=3000
        times.append(t)
        label = f"{snap_idx:04d}"
        snap_labels.append(label)

        # Read CSV; sort by cell_idx to guarantee order matches .ele connectivity
        df    = pd.read_csv(csv_path)
        df    = df.sort_values('cell_idx').reset_index(drop=True)

        grp = hf.create_group(f"Snap/{label}")
        grp.attrs['time'] = t
        grp.attrs['step'] = step

        for field in SCALAR_FIELDS:
            arr = df[field].to_numpy(dtype=np.float32)
            grp.create_dataset(field, data=arr, compression='gzip', shuffle=True)

        # Derived: von Mises stress
        sxx = df['sigma_xx'].to_numpy(dtype=np.float64)
        syy = df['sigma_yy'].to_numpy(dtype=np.float64)
        szz = df['sigma_zz'].to_numpy(dtype=np.float64)
        sxy = df['sigma_xy'].to_numpy(dtype=np.float64)
        sxz = df['sigma_xz'].to_numpy(dtype=np.float64)
        syz = df['sigma_yz'].to_numpy(dtype=np.float64)
        vm  = np.sqrt(0.5 * ((sxx - syy)**2 + (syy - szz)**2 + (szz - sxx)**2
                             + 6.0 * (sxy**2 + sxz**2 + syz**2)))
        grp.create_dataset('von_mises',   data=vm.astype(np.float32),
                           compression='gzip', shuffle=True)

        # Derived: mean (hydrostatic) stress -- positive = tension
        p_mean = (sxx + syy + szz) / 3.0
        grp.create_dataset('mean_stress', data=p_mean.astype(np.float32),
                           compression='gzip', shuffle=True)

        vm_max = float(vm.max())
        print(f"  snap {label}  t={t:.3f}s  von_Mises_max={vm_max:.3e} Pa", flush=True)

ALL_FIELDS = SCALAR_FIELDS + ['von_mises', 'mean_stress']

# -------------------------------------------------------------------------
# Write XDMF2
# -------------------------------------------------------------------------

print(f"Writing {OUT_XDMF} ...", flush=True)

lines = [
    '<?xml version="1.0" ?>',
    '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>',
    '<Xdmf Version="2.0">',
    '  <Domain>',
    '    <Grid Name="SourceBall3D" GridType="Collection" CollectionType="Temporal">',
]

for label, t in zip(snap_labels, times):
    lines += [
        f'      <Grid Name="snap_{label}" GridType="Uniform">',
        f'        <Time Value="{t:.6f}"/>',
        # Topology
        f'        <Topology TopologyType="Tetrahedron" NumberOfElements="{N_TETS}">',
        f'          <DataItem Dimensions="{N_TETS} 4" NumberType="Int" Precision="4" Format="HDF">',
        f'            {H5_BASENAME}:/Mesh/connectivity',
        f'          </DataItem>',
        f'        </Topology>',
        # Geometry
        f'        <Geometry GeometryType="XYZ">',
        f'          <DataItem Dimensions="{N_VERTS} 3" NumberType="Float" Precision="8" Format="HDF">',
        f'            {H5_BASENAME}:/Mesh/vertices',
        f'          </DataItem>',
        f'        </Geometry>',
    ]
    for field in ALL_FIELDS:
        lines += [
            f'        <Attribute Name="{field}" AttributeType="Scalar" Center="Cell">',
            f'          <DataItem Dimensions="{N_TETS}" NumberType="Float" Precision="4" Format="HDF">',
            f'            {H5_BASENAME}:/Snap/{label}/{field}',
            f'          </DataItem>',
            f'        </Attribute>',
        ]
    lines.append('      </Grid>')

lines += [
    '    </Grid>',
    '  </Domain>',
    '</Xdmf>',
]

with open(OUT_XDMF, 'w') as f:
    f.write('\n'.join(lines) + '\n')

print(f"\nDone.")
print(f"  HDF5 : {OUT_H5}")
print(f"  XDMF : {OUT_XDMF}")
print("Open source_ball_3d_vis.xdmf in ParaView to animate the 3D source ball.")
