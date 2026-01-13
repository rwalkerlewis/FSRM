#!/bin/bash
# ============================================================================
# FSRM Benchmark Runner Script
# ============================================================================
# Comprehensive benchmark suite runner for performance testing
#
# Usage:
#   ./run_benchmarks.sh [options]
#
# Options:
#   -a, --all          Run all benchmarks (default)
#   -k, --kernel       Run kernel benchmarks only
#   -p, --physics      Run physics benchmarks only
#   -g, --gpu          Run GPU benchmarks only
#   -m, --memory       Run memory/IO benchmarks only
#   -s, --scenario     Run scenario benchmarks only
#   -e, --spe          Run SPE benchmarks only
#   --scaling          Run scaling tests only
#   -n, --nproc NUM    Number of MPI processes (default: 4)
#   -o, --output DIR   Output directory for results (default: benchmark_results)
#   -v, --verbose      Verbose output
#   -h, --help         Show this help message
#
# Examples:
#   ./run_benchmarks.sh -a                    # Run all benchmarks with 4 processes
#   ./run_benchmarks.sh -p -n 8               # Run physics benchmarks with 8 processes
#   ./run_benchmarks.sh -g                    # Run GPU benchmarks
#   ./run_benchmarks.sh --spe -n 16           # Run SPE benchmarks with 16 processes
# ============================================================================

set -e  # Exit on error

# Default values
NPROC=4
OUTPUT_DIR="benchmark_results"
VERBOSE=0
RUN_ALL=1
RUN_KERNEL=0
RUN_PHYSICS=0
RUN_GPU=0
RUN_MEMORY=0
RUN_SCENARIO=0
RUN_SPE=0
RUN_SCEC=0
RUN_SCALING=0

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Helper functions
print_header() {
    echo -e "${BLUE}============================================================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}============================================================================${NC}"
}

print_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -a|--all)
            RUN_ALL=1
            shift
            ;;
        -k|--kernel)
            RUN_ALL=0
            RUN_KERNEL=1
            shift
            ;;
        -p|--physics)
            RUN_ALL=0
            RUN_PHYSICS=1
            shift
            ;;
        -g|--gpu)
            RUN_ALL=0
            RUN_GPU=1
            shift
            ;;
        -m|--memory)
            RUN_ALL=0
            RUN_MEMORY=1
            shift
            ;;
        -s|--scenario)
            RUN_ALL=0
            RUN_SCENARIO=1
            shift
            ;;
        -e|--spe)
            RUN_ALL=0
            RUN_SPE=1
            shift
            ;;
        --scec)
            RUN_ALL=0
            RUN_SCEC=1
            shift
            ;;
        --scaling)
            RUN_ALL=0
            RUN_SCALING=1
            shift
            ;;
        -n|--nproc)
            NPROC="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -v|--verbose)
            VERBOSE=1
            shift
            ;;
        -h|--help)
            grep "^#" "$0" | grep -v "#!/bin/bash" | sed 's/^# //'
            exit 0
            ;;
        *)
            print_error "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Check if we're in the correct directory
if [ ! -f "run_benchmarks.sh" ]; then
    print_error "This script must be run from the tests directory"
    exit 1
fi

# Check if executables exist
BUILD_DIR="../build"
if [ ! -d "$BUILD_DIR" ]; then
    print_error "Build directory not found. Please build the project first."
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Get timestamp for results
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

print_header "FSRM Benchmark Suite"
echo "Timestamp:        $TIMESTAMP"
echo "MPI Processes:    $NPROC"
echo "Output Directory: $OUTPUT_DIR"
echo ""

# Check for GPU support
GPU_AVAILABLE=0
if command -v nvidia-smi &> /dev/null; then
    GPU_AVAILABLE=1
    print_info "GPU detected:"
    nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv,noheader | head -1
    echo ""
fi

# ============================================================================
# Run Benchmarks
# ============================================================================

run_benchmark() {
    local name=$1
    local filter=$2
    local output_file="${OUTPUT_DIR}/${name}_${TIMESTAMP}.log"
    
    print_info "Running $name..."
    
    if [ $VERBOSE -eq 1 ]; then
        mpirun -np $NPROC "$BUILD_DIR/tests/run_performance_tests" --gtest_filter="$filter" | tee "$output_file"
    else
        mpirun -np $NPROC "$BUILD_DIR/tests/run_performance_tests" --gtest_filter="$filter" > "$output_file" 2>&1
    fi
    
    if [ $? -eq 0 ]; then
        print_info "$name completed successfully"
    else
        print_warning "$name had some failures (see $output_file)"
    fi
    echo ""
}

# Kernel Benchmarks
if [ $RUN_ALL -eq 1 ] || [ $RUN_KERNEL -eq 1 ]; then
    print_header "Kernel Performance Benchmarks"
    run_benchmark "kernel_benchmarks" "BenchmarkTest.*"
fi

# Physics Benchmarks
if [ $RUN_ALL -eq 1 ] || [ $RUN_PHYSICS -eq 1 ]; then
    print_header "Physics Benchmarks"
    run_benchmark "physics_benchmarks" "PhysicsBenchmark.*"
fi

# GPU Benchmarks
if [ $RUN_ALL -eq 1 ] || [ $RUN_GPU -eq 1 ]; then
    print_header "GPU Benchmarks"
    if [ $GPU_AVAILABLE -eq 1 ]; then
        run_benchmark "gpu_benchmarks" "GPUBenchmark.*"
    else
        print_warning "GPU not available, skipping GPU benchmarks"
    fi
fi

# Memory/IO Benchmarks
if [ $RUN_ALL -eq 1 ] || [ $RUN_MEMORY -eq 1 ]; then
    print_header "Memory and I/O Benchmarks"
    run_benchmark "memory_io_benchmarks" "MemoryIOBenchmark.*"
fi

# Scaling Tests
if [ $RUN_ALL -eq 1 ] || [ $RUN_SCALING -eq 1 ]; then
    print_header "Parallel Scaling Tests"
    run_benchmark "scaling_tests" "ScalingTest.*"
fi

# SCEC Benchmarks
if [ $RUN_ALL -eq 1 ] || [ $RUN_SCEC -eq 1 ]; then
    print_header "SCEC Benchmarks"
    run_benchmark "scec_benchmarks" "SCECBenchmark.*"
fi

# Scenario Benchmarks (long running)
if [ $RUN_SCENARIO -eq 1 ]; then
    print_header "Scenario Benchmarks"
    print_warning "These benchmarks may take a long time to complete"
    run_benchmark "scenario_benchmarks" "ScenarioBenchmark.*"
fi

# SCEC Full Simulation Benchmarks (if --scec flag given)
if [ $RUN_SCEC -eq 1 ]; then
    print_header "SCEC Full Simulation Benchmarks"
    print_warning "These simulations may take several hours"
    
    # TPV5 - Strike-slip rupture
    if [ -f "$BUILD_DIR/examples/scec_tpv5" ]; then
        print_info "Running SCEC TPV5 (strike-slip dynamic rupture)..."
        output_file="${OUTPUT_DIR}/scec_tpv5_${TIMESTAMP}.log"
        if [ $VERBOSE -eq 1 ]; then
            mpirun -np $NPROC "$BUILD_DIR/examples/scec_tpv5" -c config/scec_tpv5.config | tee "$output_file"
        else
            mpirun -np $NPROC "$BUILD_DIR/examples/scec_tpv5" -c config/scec_tpv5.config > "$output_file" 2>&1
        fi
        print_info "SCEC TPV5 completed"
    else
        print_warning "SCEC TPV5 executable not found"
    fi
    echo ""
    
    # TPV10 - Branching fault
    if [ -f "$BUILD_DIR/examples/scec_tpv10" ]; then
        print_info "Running SCEC TPV10 (branching fault)..."
        output_file="${OUTPUT_DIR}/scec_tpv10_${TIMESTAMP}.log"
        if [ $VERBOSE -eq 1 ]; then
            mpirun -np $NPROC "$BUILD_DIR/examples/scec_tpv10" -c config/scec_tpv10.config | tee "$output_file"
        else
            mpirun -np $NPROC "$BUILD_DIR/examples/scec_tpv10" -c config/scec_tpv10.config > "$output_file" 2>&1
        fi
        print_info "SCEC TPV10 completed"
    else
        print_warning "SCEC TPV10 executable not found"
    fi
    echo ""
    
    # LOH.1 - Wave propagation
    if [ -f "$BUILD_DIR/examples/scec_loh1" ]; then
        print_info "Running SCEC LOH.1 (wave propagation)..."
        output_file="${OUTPUT_DIR}/scec_loh1_${TIMESTAMP}.log"
        if [ $VERBOSE -eq 1 ]; then
            mpirun -np $NPROC "$BUILD_DIR/examples/scec_loh1" -c config/scec_loh1.config | tee "$output_file"
        else
            mpirun -np $NPROC "$BUILD_DIR/examples/scec_loh1" -c config/scec_loh1.config > "$output_file" 2>&1
        fi
        print_info "SCEC LOH.1 completed"
    else
        print_warning "SCEC LOH.1 executable not found"
    fi
    echo ""
    
    # TPV16 - Rough fault (most demanding, skip if < 16 processes)
    if [ $NPROC -lt 16 ]; then
        print_warning "SCEC TPV16 (rough fault) recommended with >= 16 MPI processes"
        print_warning "Skipping TPV16 (current: $NPROC processes). Use -n 16 or more to run."
    else
        if [ -f "$BUILD_DIR/examples/scec_tpv16" ]; then
            print_info "Running SCEC TPV16 (rough fault - high resolution)..."
            output_file="${OUTPUT_DIR}/scec_tpv16_${TIMESTAMP}.log"
            if [ $VERBOSE -eq 1 ]; then
                mpirun -np $NPROC "$BUILD_DIR/examples/scec_tpv16" -c config/scec_tpv16.config | tee "$output_file"
            else
                mpirun -np $NPROC "$BUILD_DIR/examples/scec_tpv16" -c config/scec_tpv16.config > "$output_file" 2>&1
            fi
            print_info "SCEC TPV16 completed"
        else
            print_warning "SCEC TPV16 executable not found"
        fi
    fi
    echo ""
fi

# SPE Benchmarks (very long running)
if [ $RUN_SPE -eq 1 ]; then
    print_header "SPE Benchmarks"
    print_warning "These benchmarks may take a very long time to complete"
    
    # SPE1
    if [ -f "$BUILD_DIR/examples/spe1" ]; then
        print_info "Running SPE1 benchmark..."
        output_file="${OUTPUT_DIR}/spe1_${TIMESTAMP}.log"
        if [ $VERBOSE -eq 1 ]; then
            mpirun -np $NPROC "$BUILD_DIR/examples/spe1" -c config/spe1_benchmark.config | tee "$output_file"
        else
            mpirun -np $NPROC "$BUILD_DIR/examples/spe1" -c config/spe1_benchmark.config > "$output_file" 2>&1
        fi
        print_info "SPE1 completed"
    else
        print_warning "SPE1 executable not found"
    fi
    echo ""
    
    # SPE3
    if [ -f "$BUILD_DIR/examples/spe3" ]; then
        print_info "Running SPE3 benchmark..."
        output_file="${OUTPUT_DIR}/spe3_${TIMESTAMP}.log"
        if [ $VERBOSE -eq 1 ]; then
            mpirun -np $NPROC "$BUILD_DIR/examples/spe3" -c config/spe3_benchmark.config | tee "$output_file"
        else
            mpirun -np $NPROC "$BUILD_DIR/examples/spe3" -c config/spe3_benchmark.config > "$output_file" 2>&1
        fi
        print_info "SPE3 completed"
    else
        print_warning "SPE3 executable not found"
    fi
    echo ""
    
    # SPE9
    if [ -f "$BUILD_DIR/examples/spe9" ]; then
        print_info "Running SPE9 benchmark..."
        output_file="${OUTPUT_DIR}/spe9_${TIMESTAMP}.log"
        if [ $VERBOSE -eq 1 ]; then
            mpirun -np $NPROC "$BUILD_DIR/examples/spe9" -c config/spe9_benchmark.config | tee "$output_file"
        else
            mpirun -np $NPROC "$BUILD_DIR/examples/spe9" -c config/spe9_benchmark.config > "$output_file" 2>&1
        fi
        print_info "SPE9 completed"
    else
        print_warning "SPE9 executable not found"
    fi
    echo ""
    
    # SPE10 (very large, recommend more processes)
    if [ $NPROC -lt 16 ]; then
        print_warning "SPE10 benchmark recommended with >= 16 MPI processes"
        print_warning "Skipping SPE10 (current: $NPROC processes). Use -n 16 or more to run."
    else
        if [ -f "$BUILD_DIR/examples/spe10" ]; then
            print_info "Running SPE10 benchmark..."
            output_file="${OUTPUT_DIR}/spe10_${TIMESTAMP}.log"
            if [ $VERBOSE -eq 1 ]; then
                mpirun -np $NPROC "$BUILD_DIR/examples/spe10" -c config/spe10_benchmark.config | tee "$output_file"
            else
                mpirun -np $NPROC "$BUILD_DIR/examples/spe10" -c config/spe10_benchmark.config > "$output_file" 2>&1
            fi
            print_info "SPE10 completed"
        else
            print_warning "SPE10 executable not found"
        fi
    fi
    echo ""
fi

# ============================================================================
# Generate Summary
# ============================================================================

print_header "Benchmark Summary"

SUMMARY_FILE="${OUTPUT_DIR}/summary_${TIMESTAMP}.txt"

{
    echo "FSRM Benchmark Results"
    echo "====================="
    echo ""
    echo "Timestamp:     $TIMESTAMP"
    echo "MPI Processes: $NPROC"
    echo "System:        $(uname -a)"
    echo ""
    
    if [ $GPU_AVAILABLE -eq 1 ]; then
        echo "GPU:           $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
        echo ""
    fi
    
    echo "Results Files:"
    echo ""
    ls -lh "${OUTPUT_DIR}"/*_${TIMESTAMP}.log
    echo ""
    
    echo "Quick Performance Metrics:"
    echo ""
    
    # Extract some key metrics from log files
    if [ -f "${OUTPUT_DIR}/kernel_benchmarks_${TIMESTAMP}.log" ]; then
        echo "=== Kernel Benchmarks ==="
        grep -E "(Single phase kernel|Geomechanics kernel|μs/eval)" "${OUTPUT_DIR}/kernel_benchmarks_${TIMESTAMP}.log" | tail -5
        echo ""
    fi
    
    if [ -f "${OUTPUT_DIR}/physics_benchmarks_${TIMESTAMP}.log" ]; then
        echo "=== Physics Benchmarks ==="
        grep -E "(Poroelasticity|Fracture|Wave|μs/eval)" "${OUTPUT_DIR}/physics_benchmarks_${TIMESTAMP}.log" | tail -5
        echo ""
    fi
    
    if [ $GPU_AVAILABLE -eq 1 ] && [ -f "${OUTPUT_DIR}/gpu_benchmarks_${TIMESTAMP}.log" ]; then
        echo "=== GPU Benchmarks ==="
        grep -E "(Speedup|Bandwidth|GB/s)" "${OUTPUT_DIR}/gpu_benchmarks_${TIMESTAMP}.log" | tail -10
        echo ""
    fi
    
    echo "Complete results available in: $OUTPUT_DIR"
    
} | tee "$SUMMARY_FILE"

print_info "Summary saved to: $SUMMARY_FILE"
echo ""

print_header "Benchmark Suite Complete"
print_info "All results saved to: $OUTPUT_DIR"
print_info "Summary: $SUMMARY_FILE"

exit 0
