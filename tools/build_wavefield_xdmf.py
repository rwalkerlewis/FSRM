#!/usr/bin/env python3
"""Build a ParaView-loadable 3D temporal XDMF wrapping FSRM solution.h5.

Usage:
    python3 build_wavefield_xdmf.py <output_dir>

Reads `<output_dir>/solution.h5` and emits `<output_dir>/wavefield.xdmf`
referencing it via hyperslab DataItems with shared Domain-level geometry and
connectivity. Only snapshots with nonzero displacement are written.

Also computes a per-snapshot displacement magnitude and writes it next to the
vector field so ParaView shows scalar contour data without a Calculator step.
"""
import sys, os
import h5py
import numpy as np


def build(out_dir: str) -> None:
    h5_path = os.path.join(out_dir, 'solution.h5')
    with h5py.File(h5_path, 'r') as h:
        v = h['geometry/vertices'][:]
        cells = h['viz/topology/cells'][:]
        u = h['fields/solution'][:]
    T, Nv, _ = u.shape
    Nc = cells.shape[0]

    # Identify snapshots with real data (output_frequency strides leave zeros).
    pk = np.abs(u).max(axis=(1, 2))
    nz = np.where(pk > 0)[0]
    if len(nz) == 0:
        print(f'{h5_path}: no nonzero displacement snapshots; '
              'XDMF/aux not written')
        return

    # Compute and store per-snapshot |displacement| as a scalar dataset.
    mag = np.zeros((len(nz), Nv), dtype=np.float64)
    for k, i in enumerate(nz):
        mag[k] = np.linalg.norm(u[i], axis=1)
    aux_path = os.path.join(out_dir, 'wavefield_aux.h5')
    with h5py.File(aux_path, 'w') as g:
        g.create_dataset('umag', data=mag)
        g.create_dataset('step_index', data=nz)

    # Real wall time per FEM step is read from the config. The HDF5 time
    # dataset is broken upstream, so derive from step number and dt_initial=1ms.
    dt_step = 0.001

    xdmf = os.path.join(out_dir, 'wavefield.xdmf')
    nz_count = len(nz)
    with open(xdmf, 'w') as f:
        f.write('<?xml version="1.0" ?>\n')
        f.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n')
        f.write('<Xdmf Version="3.0">\n<Domain>\n')
        f.write(f'  <DataItem Name="verts" Format="HDF" NumberType="Float" '
                f'Precision="8" Dimensions="{Nv} 3">'
                f'solution.h5:/geometry/vertices</DataItem>\n')
        f.write(f'  <DataItem Name="cells" Format="HDF" NumberType="Int" '
                f'Precision="4" Dimensions="{Nc} 8">'
                f'solution.h5:/viz/topology/cells</DataItem>\n')
        f.write('  <Grid Name="TimeSeries" GridType="Collection" '
                'CollectionType="Temporal">\n')
        for k, i in enumerate(nz):
            t = float(i) * dt_step
            f.write(f'    <Grid Name="t{i}">\n')
            f.write(f'      <Time Value="{t}"/>\n')
            f.write(f'      <Topology TopologyType="Hexahedron" '
                    f'NumberOfElements="{Nc}">\n')
            f.write('        <DataItem Reference="XML">'
                    '/Xdmf/Domain/DataItem[@Name="cells"]</DataItem>\n')
            f.write('      </Topology>\n')
            f.write('      <Geometry GeometryType="XYZ">\n')
            f.write('        <DataItem Reference="XML">'
                    '/Xdmf/Domain/DataItem[@Name="verts"]</DataItem>\n')
            f.write('      </Geometry>\n')
            # Vector displacement via hyperslab into the (T, Nv, 3) cube.
            f.write('      <Attribute Name="displacement" Center="Node" '
                    'AttributeType="Vector">\n')
            f.write(f'        <DataItem ItemType="HyperSlab" '
                    f'Dimensions="{Nv} 3">\n')
            f.write(f'          <DataItem Dimensions="3 3" Format="XML">'
                    f'{i} 0 0 1 1 1 1 {Nv} 3</DataItem>\n')
            f.write(f'          <DataItem Format="HDF" NumberType="Float" '
                    f'Precision="8" Dimensions="{T} {Nv} 3">'
                    f'solution.h5:/fields/solution</DataItem>\n')
            f.write('        </DataItem>\n')
            f.write('      </Attribute>\n')
            # Scalar |u| from the auxiliary file.
            f.write('      <Attribute Name="umag" Center="Node" '
                    'AttributeType="Scalar">\n')
            f.write(f'        <DataItem ItemType="HyperSlab" '
                    f'Dimensions="{Nv}">\n')
            f.write(f'          <DataItem Dimensions="3 2" Format="XML">'
                    f'{k} 0 1 1 1 {Nv}</DataItem>\n')
            f.write(f'          <DataItem Format="HDF" NumberType="Float" '
                    f'Precision="8" Dimensions="{nz_count} {Nv}">'
                    f'wavefield_aux.h5:/umag</DataItem>\n')
            f.write('        </DataItem>\n')
            f.write('      </Attribute>\n')
            f.write('    </Grid>\n')
        f.write('  </Grid>\n</Domain>\n</Xdmf>\n')
    print(f'wrote {xdmf}  ({nz_count} snapshots, peak |u|={mag.max():.3e} m)')
    print(f'wrote {aux_path}')


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(1)
    build(sys.argv[1])
