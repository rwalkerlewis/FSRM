#!/usr/bin/env bash
# Pass-7 amplitude diagnostic.
#
# Runs the Sedan 1962 historic-nuclear fixture twice: once with
# solver_kind = CLOSED_FORM (the pass-5 RDP-driven baseline) and once
# with solver_kind = RADIAL_LAGRANGIAN under the pass-7 cavity-EOS /
# physics-based-init defaults. Prints the peak |M0_iso_dot| from each
# near_field_history.csv side-by-side with the ratio.
#
# Usage:
#   scripts/pass7_amplitude_diagnostic.sh [radial_cells]
#
# Default radial_cells = 200 (pass-7 production resolution is 800; the
# default keeps the diagnostic fast while still revealing the
# amplitude band).

set -euo pipefail

CELLS="${1:-200}"
WORKDIR="$(mktemp -d -t pass7_diag_XXXXXX)"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BUILD_DIR="${ROOT}/build"
EXE="${BUILD_DIR}/fsrm"

if [[ ! -x "${EXE}" ]]; then
    echo "fsrm binary not found at ${EXE}; build the project first" >&2
    exit 1
fi

write_config() {
    local solver_kind="$1"
    local out_dir="$2"
    cat <<EOC
[SIMULATION]
name = pass7_diag_${solver_kind,,}
start_time = 0.0
end_time = 0.05
dt_initial = 0.001
dt_min = 0.0001
dt_max = 0.005
max_timesteps = 100
output_frequency = 50
output_format = HDF5
fluid_model = NONE
solid_model = ELASTIC
enable_geomechanics = true
enable_faults = false
enable_elastodynamics = true
rtol = 1.0e-6
atol = 1.0e-8
max_nonlinear_iterations = 20

[GRID]
nx = 4
ny = 4
nz = 4
Lx = 4000.0
Ly = 4000.0
Lz = 2000.0

[ROCK]
density = 2650.0
lambda = 1.87e10
shear_modulus = 1.34e10

[EXPLOSION_SOURCE]
type = UNDERGROUND_NUCLEAR
yield_kt = 104.0
depth_of_burial = 194.0
location_x = 2000.0
location_y = 2000.0
location_z = 1806.0
onset_time = 0.0
rise_time = 0.01
cavity_overpressure = 1.0e10
medium_type = ALLUVIUM

[NEAR_FIELD_SOURCE]
mode = DYNAMIC_PLASTIC
solver_kind = ${solver_kind}
radial_cells = ${CELLS}
elastic_radius_factor = 3.0
near_field_dt = 1.0e-5
output_cadence_microseconds = 1000
profile_output_cadence_microseconds = 5000

[BOUNDARY_CONDITIONS]
bottom = free
sides = free
top = free

[ABSORBING_BC]
enabled = true
x_min = true
x_max = true
y_min = true
y_max = true
z_min = true
z_max = false

[SEISMOMETERS]
enabled = true
formats = SAC
output_dir = ${out_dir}
default_quantity = DISPLACEMENT
default_sample_rate_hz = 200.0

[SEISMOMETER_1]
sta = SPALL
location_xyz = 2000.0,2000.0,2000.0
EOC
}

run_one() {
    local solver_kind="$1"
    local out_dir="${WORKDIR}/${solver_kind,,}_out"
    mkdir -p "${out_dir}"
    local cfg="${WORKDIR}/${solver_kind,,}.config"
    write_config "${solver_kind}" "${out_dir}" > "${cfg}"
    (cd "${WORKDIR}" && "${EXE}" -c "${cfg}" > "${WORKDIR}/${solver_kind,,}.log" 2>&1) || {
        echo "fsrm run failed for ${solver_kind} (see ${WORKDIR}/${solver_kind,,}.log)" >&2
        exit 1
    }
    echo "${out_dir}/near_field_history.csv"
}

peak_m0iso_dot() {
    local csv="$1"
    awk -F, '
        /^#/ || /^t,/ { next }
        NF >= 10 {
            v = $10
            if (v < 0) v = -v
            if (v > peak) peak = v
        }
        END { printf "%.6e\n", peak }
    ' "${csv}"
}

CSV_CF="$(run_one CLOSED_FORM)"
CSV_RL="$(run_one RADIAL_LAGRANGIAN)"
P_CF="$(peak_m0iso_dot "${CSV_CF}")"
P_RL="$(peak_m0iso_dot "${CSV_RL}")"
RATIO="$(awk -v a="${P_RL}" -v b="${P_CF}" 'BEGIN { if (b > 0) printf "%.4f\n", a / b; else print "n/a" }')"

echo "Pass-7 amplitude diagnostic (Sedan 1962, radial_cells=${CELLS}):"
echo "  CLOSED_FORM      peak |M0_iso_dot| = ${P_CF} N*m/s"
echo "  RADIAL_LAGRANGIAN peak |M0_iso_dot| = ${P_RL} N*m/s"
echo "  Ratio (RL / CF)                    = ${RATIO}"
echo "Workdir: ${WORKDIR}"
