#!/usr/bin/env python3
"""Generate a synthetic binary velocity model for testing.

Binary format (little-endian):
  Header (44 bytes):
    nx, ny, nz       (3 x int32)     = 12 bytes
    x_min, x_max     (2 x float32)   =  8 bytes
    y_min, y_max     (2 x float32)   =  8 bytes
    z_min, z_max     (2 x float32)   =  8 bytes
    padding           (2 x float32)  =  8 bytes
  Data (nx * ny * nz * 3 * float32):
    For each (ix, iy, iz) in C order (iz fastest):
      Vp, Vs, rho    (3 x float32)

Usage:
  python3 generate_velocity_model.py output.bin [--nx 10] [--ny 10] [--nz 10]
"""

import argparse
import struct
import sys


def generate_layered_model(filename, nx=10, ny=10, nz=10,
                           x_range=(0.0, 10000.0),
                           y_range=(0.0, 10000.0),
                           z_range=(0.0, 5000.0)):
    """Generate a two-layer velocity model.

    Top half:    Vp=3000, Vs=1732, rho=2200 (soft sediment)
    Bottom half: Vp=6000, Vs=3464, rho=2700 (hard basement)
    """
    z_mid = (z_range[0] + z_range[1]) / 2.0

    with open(filename, 'wb') as f:
        # Header
        f.write(struct.pack('<iii', nx, ny, nz))
        f.write(struct.pack('<ffffff',
                            x_range[0], x_range[1],
                            y_range[0], y_range[1],
                            z_range[0], z_range[1]))
        f.write(struct.pack('<ff', 0.0, 0.0))  # padding

        # Data
        for ix in range(nx):
            for iy in range(ny):
                for iz in range(nz):
                    if nz > 1:
                        z = z_range[0] + (z_range[1] - z_range[0]) * iz / (nz - 1)
                    else:
                        z = z_mid
                    if z > z_mid:
                        f.write(struct.pack('<fff', 3000.0, 1732.0, 2200.0))
                    else:
                        f.write(struct.pack('<fff', 6000.0, 3464.0, 2700.0))

    total_bytes = 44 + nx * ny * nz * 12
    print(f"Wrote {filename}: {nx}x{ny}x{nz} grid, {total_bytes} bytes")


def generate_gradient_model(filename, nx=10, ny=10, nz=20,
                            x_range=(0.0, 10000.0),
                            y_range=(0.0, 10000.0),
                            z_range=(0.0, 5000.0)):
    """Generate a velocity model with a linear velocity gradient with depth.

    Vp increases linearly from 2000 m/s at top to 6000 m/s at bottom.
    Vs = Vp / 1.73 (constant Vp/Vs ratio).
    rho from Gardner's relation: rho = 310 * Vp^0.25
    """
    with open(filename, 'wb') as f:
        # Header
        f.write(struct.pack('<iii', nx, ny, nz))
        f.write(struct.pack('<ffffff',
                            x_range[0], x_range[1],
                            y_range[0], y_range[1],
                            z_range[0], z_range[1]))
        f.write(struct.pack('<ff', 0.0, 0.0))

        # Data
        for ix in range(nx):
            for iy in range(ny):
                for iz in range(nz):
                    if nz > 1:
                        frac = iz / (nz - 1)
                    else:
                        frac = 0.5
                    # Velocity increases with depth (lower z = deeper)
                    vp = 6000.0 - frac * 4000.0  # 6000 at z=0, 2000 at z=max
                    vs = vp / 1.73
                    rho = 310.0 * (vp ** 0.25)
                    f.write(struct.pack('<fff', vp, vs, rho))

    total_bytes = 44 + nx * ny * nz * 12
    print(f"Wrote {filename}: {nx}x{ny}x{nz} gradient model, {total_bytes} bytes")


def main():
    parser = argparse.ArgumentParser(
        description='Generate synthetic binary velocity model')
    parser.add_argument('output', help='Output .bin filename')
    parser.add_argument('--nx', type=int, default=10)
    parser.add_argument('--ny', type=int, default=10)
    parser.add_argument('--nz', type=int, default=10)
    parser.add_argument('--model', choices=['layered', 'gradient'],
                        default='layered',
                        help='Model type (default: layered)')
    parser.add_argument('--xmin', type=float, default=0.0)
    parser.add_argument('--xmax', type=float, default=10000.0)
    parser.add_argument('--ymin', type=float, default=0.0)
    parser.add_argument('--ymax', type=float, default=10000.0)
    parser.add_argument('--zmin', type=float, default=0.0)
    parser.add_argument('--zmax', type=float, default=5000.0)
    args = parser.parse_args()

    x_range = (args.xmin, args.xmax)
    y_range = (args.ymin, args.ymax)
    z_range = (args.zmin, args.zmax)

    if args.model == 'layered':
        generate_layered_model(args.output, args.nx, args.ny, args.nz,
                               x_range, y_range, z_range)
    else:
        generate_gradient_model(args.output, args.nx, args.ny, args.nz,
                                x_range, y_range, z_range)


if __name__ == '__main__':
    main()
