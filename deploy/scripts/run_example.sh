#!/bin/bash
# Script to run example simulations

set -e

echo "=========================================="
echo "FSRM Example Runner"
echo "=========================================="
echo ""

# Default values
NPROCS=${NPROCS:-$(nproc)}
OUTPUT_DIR=${OUTPUT_DIR:-./output}
CONFIG_DIR=${CONFIG_DIR:-../config}

# Color codes
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m'

print_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_header() {
    echo -e "${BLUE}========================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}========================================${NC}"
}

# Create output directory
mkdir -p $OUTPUT_DIR

# Example menu
echo "Select an example to run:"
echo ""
echo "1) Shale Reservoir with Hydraulic Fracturing"
echo "2) Enhanced Geothermal System (EGS)"
echo "3) Induced Seismicity from Injection"
echo "4) CO2 Geological Storage"
echo "5) LEFM Fracture Growth"
echo "6) Wave Propagation in Poroelastic Media"
echo "7) Static Triggered Seismicity"
echo "8) Default Configuration (quick test)"
echo ""
read -p "Enter choice [1-8]: " CHOICE

case $CHOICE in
    1)
        CONFIG="shale_reservoir.config"
        DESC="Shale Reservoir Simulation"
        ;;
    2)
        CONFIG="geothermal.config"
        DESC="Geothermal System"
        ;;
    3)
        CONFIG="induced_seismicity.config"
        DESC="Induced Seismicity"
        ;;
    4)
        CONFIG="co2_storage.config"
        DESC="CO2 Storage"
        ;;
    5)
        CONFIG="lefm_fracture_growth.config"
        DESC="LEFM Fracture Growth"
        ;;
    6)
        CONFIG="poroelastodynamic_waves.config"
        DESC="Wave Propagation"
        ;;
    7)
        CONFIG="static_triggered_seismicity.config"
        DESC="Static Triggered Seismicity"
        ;;
    8)
        CONFIG="default.config"
        DESC="Default Configuration"
        ;;
    *)
        echo "Invalid choice"
        exit 1
        ;;
esac

# Number of processors
read -p "Number of MPI processes (default: $NPROCS): " INPUT_PROCS
NPROCS=${INPUT_PROCS:-$NPROCS}

# Run simulation
print_header "Running: $DESC"
print_info "Configuration: $CONFIG"
print_info "MPI processes: $NPROCS"
print_info "Output directory: $OUTPUT_DIR"
echo ""

if [ ! -f "$CONFIG_DIR/$CONFIG" ]; then
    echo "Config file not found: $CONFIG_DIR/$CONFIG"
    exit 1
fi

# Run with timing
START_TIME=$(date +%s)

mpirun -np $NPROCS fsrm \
    -c $CONFIG_DIR/$CONFIG \
    -o $OUTPUT_DIR/$(basename $CONFIG .config) \
    -log_view \
    -ts_monitor \
    -snes_monitor

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

print_header "Simulation Complete"
print_info "Total time: $ELAPSED seconds"
print_info "Output files: $OUTPUT_DIR/$(basename $CONFIG .config)*"
echo ""

# Check for output files
if [ -f "$OUTPUT_DIR/$(basename $CONFIG .config).log" ]; then
    print_info "View log: cat $OUTPUT_DIR/$(basename $CONFIG .config).log"
fi

if ls $OUTPUT_DIR/*.vtu 1> /dev/null 2>&1; then
    print_info "Visualize with ParaView: paraview $OUTPUT_DIR/*.vtu"
fi

if ls $OUTPUT_DIR/*.h5 1> /dev/null 2>&1; then
    print_info "HDF5 files available for post-processing"
fi

echo ""
print_info "To run another example, execute: $0"
