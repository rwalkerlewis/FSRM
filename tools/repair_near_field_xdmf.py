#!/usr/bin/env python3
"""
Repair near_field_profile.xdmf files that replicate the 200-cell
polyline connectivity inline per snapshot. We:

  1. Stream-scan the old XDMF for the (Grid Name, Time) pairs.
  2. Read num_cells from the companion HDF5 file.
  3. Rewrite the XDMF using shared Domain-level connectivity / Y / Z
     DataItems and per-snapshot Reference="XML" pointers.

Drop-in for examples/<event>/output*/near_field_profile.xdmf produced
by writeRadialProfilesXDMF before the pass-14 fix.
"""
import os, sys, re, h5py

GRID_RE = re.compile(rb'<Grid Name="(profile_\d+)" GridType="Uniform">')
TIME_RE = re.compile(rb'<Time Value="([^"]+)"/>')

def repair(xdmf_path):
    h5_basename = "near_field_profile.h5"
    h5_path = os.path.join(os.path.dirname(xdmf_path), h5_basename)
    if not os.path.exists(h5_path):
        print(f"  SKIP: {h5_path} missing")
        return False

    with h5py.File(h5_path, "r") as h:
        N = int(h["/num_cells"][()])
        snap_ids = sorted(k for k in h["/profiles"].keys())
    Nface = N + 1
    print(f"  HDF5: N={N}, snapshots={len(snap_ids)}")

    # Stream-extract (grid_name -> time) from the old XDMF.
    # Snapshots appear in order; we just zip them.
    times = []
    cur_grid = None
    with open(xdmf_path, "rb") as f:
        for line in f:
            m = GRID_RE.search(line)
            if m:
                cur_grid = m.group(1).decode()
                continue
            m = TIME_RE.search(line)
            if m and cur_grid is not None:
                times.append(m.group(1).decode())
                cur_grid = None
    if len(times) != len(snap_ids):
        print(f"  WARN: extracted {len(times)} times from XDMF "
              f"but HDF5 has {len(snap_ids)} snapshots. Truncating.")
    n = min(len(times), len(snap_ids))

    backup = xdmf_path + ".bloated.bak"
    if not os.path.exists(backup):
        os.rename(xdmf_path, backup)
    else:
        # already repaired once; the .bak is the source of truth
        pass

    with open(xdmf_path, "w") as out:
        out.write('<?xml version="1.0" ?>\n'
                  '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n'
                  '<Xdmf Version="3.0">\n'
                  '  <Domain>\n')
        out.write(f'    <DataItem Name="connectivity" Dimensions="{N} 2"'
                  ' Format="XML" DataType="Int">\n')
        for j in range(N):
            out.write(f"      {j} {j+1}\n")
        out.write('    </DataItem>\n')
        out.write(f'    <DataItem Name="y_zero" Dimensions="{Nface}"'
                  ' NumberType="Float" Precision="8" Format="XML">\n')
        for _ in range(Nface): out.write("      0.0\n")
        out.write('    </DataItem>\n')
        out.write(f'    <DataItem Name="z_zero" Dimensions="{Nface}"'
                  ' NumberType="Float" Precision="8" Format="XML">\n')
        for _ in range(Nface): out.write("      0.0\n")
        out.write('    </DataItem>\n')
        out.write('    <Grid Name="NearFieldProfileSeries"'
                  ' GridType="Collection" CollectionType="Temporal">\n')

        for i in range(n):
            g = snap_ids[i]
            t = times[i]
            out.write(f'      <Grid Name="profile_{g}" GridType="Uniform">\n')
            out.write(f'        <Time Value="{t}"/>\n')
            out.write(f'        <Topology TopologyType="Polyline"'
                      f' NumberOfElements="{N}" NodesPerElement="2">\n')
            out.write(f'          <DataItem Dimensions="{N} 2"'
                      ' Format="XML" DataType="Int" Reference="XML">'
                      "/Xdmf/Domain/DataItem[@Name='connectivity']"
                      '</DataItem>\n')
            out.write('        </Topology>\n')
            out.write('        <Geometry GeometryType="X_Y_Z">\n')
            out.write(f'          <DataItem Dimensions="{Nface}"'
                      ' NumberType="Float" Precision="8" Format="HDF">\n')
            out.write(f'            {h5_basename}:/profiles/{g}/r\n')
            out.write('          </DataItem>\n')
            out.write(f'          <DataItem Dimensions="{Nface}"'
                      ' NumberType="Float" Precision="8"'
                      ' Format="XML" Reference="XML">'
                      "/Xdmf/Domain/DataItem[@Name='y_zero']"
                      '</DataItem>\n')
            out.write(f'          <DataItem Dimensions="{Nface}"'
                      ' NumberType="Float" Precision="8"'
                      ' Format="XML" Reference="XML">'
                      "/Xdmf/Domain/DataItem[@Name='z_zero']"
                      '</DataItem>\n')
            out.write('        </Geometry>\n')
            for attr in ("rho", "p", "sigma_rr", "sigma_tt",
                         "eps_p", "damage", "yield_indicator"):
                out.write(f'        <Attribute Name="{attr}"'
                          ' AttributeType="Scalar" Center="Cell">\n')
                out.write(f'          <DataItem Dimensions="{N}"'
                          ' NumberType="Float" Precision="8" Format="HDF">\n')
                out.write(f'            {h5_basename}:/profiles/{g}/{attr}\n')
                out.write('          </DataItem>\n')
                out.write('        </Attribute>\n')
            out.write('      </Grid>\n')
        out.write('    </Grid>\n  </Domain>\n</Xdmf>\n')

    new_size_mb = os.path.getsize(xdmf_path) / 1e6
    old_size_mb = os.path.getsize(backup) / 1e6
    print(f"  rewrote {xdmf_path}: {old_size_mb:.1f} MB -> {new_size_mb:.1f} MB")
    return True

if __name__ == "__main__":
    targets = sys.argv[1:]
    if not targets:
        print("usage: repair_near_field_xdmf.py file1.xdmf [file2.xdmf ...]")
        sys.exit(1)
    for p in targets:
        print(f"Repairing {p}")
        repair(p)
