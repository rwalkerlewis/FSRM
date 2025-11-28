#!/bin/bash
################################################################################
# SCEC Benchmark Suite Runner for FSRM
################################################################################
# This script runs all SCEC benchmarks and optionally verifies results
#
# Usage:
#   ./run_scec_suite.sh [options]
#
# Options:
#   --all              Run all benchmarks
#   --basic            Run basic benchmarks (TPV5, TPV10)
#   --advanced         Run advanced benchmarks (TPV13, TPV16)
#   --rate-state       Run rate-and-state benchmarks (TPV101, TPV104)
#   --thermal          Run thermal pressurization (TPV34)
#   --verify           Verify results against reference solutions
#   --report           Generate summary report
#   --parallel N       Run N benchmarks in parallel
#   --gpu              Use GPU acceleration
#   --help             Show this help message
################################################################################

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Default settings
RUN_ALL=false
RUN_BASIC=false
RUN_ADVANCED=false
RUN_RATE_STATE=false
RUN_THERMAL=false
VERIFY=false
REPORT=false
PARALLEL=1
USE_GPU=false
FSRM_BIN="./fsrm"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --all)
            RUN_ALL=true
            shift
            ;;
        --basic)
            RUN_BASIC=true
            shift
            ;;
        --advanced)
            RUN_ADVANCED=true
            shift
            ;;
        --rate-state)
            RUN_RATE_STATE=true
            shift
            ;;
        --thermal)
            RUN_THERMAL=true
            shift
            ;;
        --verify)
            VERIFY=true
            shift
            ;;
        --report)
            REPORT=true
            shift
            ;;
        --parallel)
            PARALLEL="$2"
            shift 2
            ;;
        --gpu)
            USE_GPU=true
            shift
            ;;
        --help)
            head -n 20 "$0" | tail -n 18
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# If no specific category selected, show help
if [ "$RUN_ALL" = false ] && [ "$RUN_BASIC" = false ] && \
   [ "$RUN_ADVANCED" = false ] && [ "$RUN_RATE_STATE" = false ] && \
   [ "$RUN_THERMAL" = false ]; then
    echo -e "${YELLOW}No benchmark category selected. Use --help for options.${NC}"
    exit 1
fi

# Benchmark lists
BASIC_BENCHMARKS=(
    "tpv5:scec_tpv5.config:5 min"
)

ADVANCED_BENCHMARKS=(
    "tpv10:scec_tpv10.config:10 min"
    "tpv13:scec_tpv13.config:30 min"
    "tpv16:scec_tpv16.config:15 min"
)

RATE_STATE_BENCHMARKS=(
    "tpv101:scec_tpv101.config:60 min"
    "tpv104:scec_tpv104.config:90 min"
)

THERMAL_BENCHMARKS=(
    "tpv34:scec_tpv34.config:20 min"
)

# Combine selected benchmarks
BENCHMARKS=()
if [ "$RUN_ALL" = true ] || [ "$RUN_BASIC" = true ]; then
    BENCHMARKS+=("${BASIC_BENCHMARKS[@]}")
fi
if [ "$RUN_ALL" = true ] || [ "$RUN_ADVANCED" = true ]; then
    BENCHMARKS+=("${ADVANCED_BENCHMARKS[@]}")
fi
if [ "$RUN_ALL" = true ] || [ "$RUN_RATE_STATE" = true ]; then
    BENCHMARKS+=("${RATE_STATE_BENCHMARKS[@]}")
fi
if [ "$RUN_ALL" = true ] || [ "$RUN_THERMAL" = true ]; then
    BENCHMARKS+=("${THERMAL_BENCHMARKS[@]}")
fi

# Check if FSRM binary exists
if [ ! -f "$FSRM_BIN" ]; then
    echo -e "${RED}Error: FSRM binary not found at $FSRM_BIN${NC}"
    echo "Please build FSRM first or specify path with FSRM_BIN environment variable"
    exit 1
fi

# Create output directory
mkdir -p benchmarks/results
RESULTS_DIR="benchmarks/results/$(date +%Y%m%d_%H%M%S)"
mkdir -p "$RESULTS_DIR"

# Log file
LOG_FILE="$RESULTS_DIR/benchmark_suite.log"
touch "$LOG_FILE"

# Function to run a single benchmark
run_benchmark() {
    local name=$1
    local config=$2
    local expected_time=$3
    
    echo -e "${BLUE}Running $name...${NC}"
    echo "Benchmark: $name" >> "$LOG_FILE"
    echo "Config: benchmarks/$config" >> "$LOG_FILE"
    echo "Start time: $(date)" >> "$LOG_FILE"
    
    local start_time=$(date +%s)
    
    # Build command
    local cmd="$FSRM_BIN -c benchmarks/$config"
    if [ "$USE_GPU" = true ]; then
        cmd="$cmd --use-gpu"
    fi
    
    # Run benchmark
    if $cmd >> "$LOG_FILE" 2>&1; then
        local end_time=$(date +%s)
        local elapsed=$((end_time - start_time))
        local elapsed_min=$((elapsed / 60))
        local elapsed_sec=$((elapsed % 60))
        
        echo -e "${GREEN}✓ $name completed in ${elapsed_min}m ${elapsed_sec}s${NC}"
        echo "Status: SUCCESS" >> "$LOG_FILE"
        echo "Elapsed time: ${elapsed_min}m ${elapsed_sec}s" >> "$LOG_FILE"
        echo "" >> "$LOG_FILE"
        
        # Copy output to results directory
        if [ -d "output/$name" ]; then
            cp -r "output/$name" "$RESULTS_DIR/"
        fi
        
        return 0
    else
        local end_time=$(date +%s)
        local elapsed=$((end_time - start_time))
        
        echo -e "${RED}✗ $name failed after $elapsed seconds${NC}"
        echo "Status: FAILED" >> "$LOG_FILE"
        echo "Elapsed time: $elapsed seconds" >> "$LOG_FILE"
        echo "" >> "$LOG_FILE"
        
        return 1
    fi
}

# Function to verify benchmark results
verify_benchmark() {
    local name=$1
    
    echo -e "${BLUE}Verifying $name...${NC}"
    
    if python3 benchmarks/compare_with_reference.py "$name" >> "$LOG_FILE" 2>&1; then
        echo -e "${GREEN}✓ $name verification passed${NC}"
        return 0
    else
        echo -e "${YELLOW}⚠ $name verification failed or incomplete${NC}"
        return 1
    fi
}

# Main execution
echo "================================"
echo "SCEC Benchmark Suite for FSRM"
echo "================================"
echo ""
echo "Configuration:"
echo "  Benchmarks to run: ${#BENCHMARKS[@]}"
echo "  Parallel jobs: $PARALLEL"
echo "  GPU acceleration: $USE_GPU"
echo "  Verification: $VERIFY"
echo "  Results directory: $RESULTS_DIR"
echo ""

# Track results
PASSED=0
FAILED=0
VERIFIED=0
VERIFICATION_FAILED=0

# Run benchmarks
if [ "$PARALLEL" -eq 1 ]; then
    # Sequential execution
    for benchmark_info in "${BENCHMARKS[@]}"; do
        IFS=':' read -r name config expected_time <<< "$benchmark_info"
        
        if run_benchmark "$name" "$config" "$expected_time"; then
            ((PASSED++))
            
            if [ "$VERIFY" = true ]; then
                if verify_benchmark "$name"; then
                    ((VERIFIED++))
                else
                    ((VERIFICATION_FAILED++))
                fi
            fi
        else
            ((FAILED++))
        fi
        
        echo ""
    done
else
    # Parallel execution
    echo "Running $PARALLEL benchmarks in parallel..."
    
    # Create temporary directory for parallel jobs
    TEMP_DIR=$(mktemp -d)
    
    # Launch benchmarks in background
    declare -a PIDS
    for benchmark_info in "${BENCHMARKS[@]}"; do
        IFS=':' read -r name config expected_time <<< "$benchmark_info"
        
        # Wait if we've reached parallel limit
        while [ $(jobs -r | wc -l) -ge $PARALLEL ]; do
            sleep 1
        done
        
        # Run in background
        (
            if run_benchmark "$name" "$config" "$expected_time"; then
                echo "PASSED" > "$TEMP_DIR/$name.status"
            else
                echo "FAILED" > "$TEMP_DIR/$name.status"
            fi
        ) &
        
        PIDS+=($!)
    done
    
    # Wait for all jobs to complete
    echo "Waiting for all benchmarks to complete..."
    for pid in "${PIDS[@]}"; do
        wait $pid
    done
    
    # Count results
    for status_file in "$TEMP_DIR"/*.status; do
        if [ -f "$status_file" ]; then
            if grep -q "PASSED" "$status_file"; then
                ((PASSED++))
            else
                ((FAILED++))
            fi
        fi
    done
    
    # Verification (sequential after parallel run)
    if [ "$VERIFY" = true ]; then
        echo ""
        echo "Verifying results..."
        for benchmark_info in "${BENCHMARKS[@]}"; do
            IFS=':' read -r name config expected_time <<< "$benchmark_info"
            if [ -f "$TEMP_DIR/$name.status" ] && grep -q "PASSED" "$TEMP_DIR/$name.status"; then
                if verify_benchmark "$name"; then
                    ((VERIFIED++))
                else
                    ((VERIFICATION_FAILED++))
                fi
            fi
        done
    fi
    
    # Cleanup
    rm -rf "$TEMP_DIR"
fi

# Generate summary
echo ""
echo "================================"
echo "Summary"
echo "================================"
echo "Total benchmarks: $((PASSED + FAILED))"
echo -e "${GREEN}Passed: $PASSED${NC}"
if [ $FAILED -gt 0 ]; then
    echo -e "${RED}Failed: $FAILED${NC}"
fi

if [ "$VERIFY" = true ]; then
    echo ""
    echo "Verification Results:"
    echo -e "${GREEN}Verified: $VERIFIED${NC}"
    if [ $VERIFICATION_FAILED -gt 0 ]; then
        echo -e "${YELLOW}Verification failed: $VERIFICATION_FAILED${NC}"
    fi
fi

echo ""
echo "Results saved to: $RESULTS_DIR"
echo "Log file: $LOG_FILE"

# Generate detailed report
if [ "$REPORT" = true ]; then
    REPORT_FILE="$RESULTS_DIR/benchmark_report.html"
    echo "Generating report: $REPORT_FILE"
    python3 benchmarks/generate_report.py "$RESULTS_DIR" > "$REPORT_FILE"
    echo "Report generated: $REPORT_FILE"
fi

# Exit with appropriate code
if [ $FAILED -gt 0 ]; then
    exit 1
else
    exit 0
fi
