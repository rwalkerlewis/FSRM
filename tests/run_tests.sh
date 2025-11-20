#!/bin/bash
# Script to run FSRM tests with various configurations

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Default values
BUILD_DIR="build"
NPROCS=1
TEST_FILTER=""
VERBOSE=0
QUICK_ONLY=0
CONVERGENCE_ONLY=0
INTEGRATION_ONLY=0

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -b|--build-dir)
            BUILD_DIR="$2"
            shift 2
            ;;
        -n|--nprocs)
            NPROCS="$2"
            shift 2
            ;;
        -f|--filter)
            TEST_FILTER="$2"
            shift 2
            ;;
        -v|--verbose)
            VERBOSE=1
            shift
            ;;
        --quick)
            QUICK_ONLY=1
            shift
            ;;
        --convergence)
            CONVERGENCE_ONLY=1
            shift
            ;;
        --integration)
            INTEGRATION_ONLY=1
            shift
            ;;
        -h|--help)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  -b, --build-dir DIR    Build directory (default: build)"
            echo "  -n, --nprocs N         Number of MPI processes (default: 1)"
            echo "  -f, --filter PATTERN   GoogleTest filter pattern"
            echo "  -v, --verbose          Verbose output"
            echo "  --quick                Run only quick unit tests"
            echo "  --convergence          Run only convergence tests"
            echo "  --integration          Run only integration tests"
            echo "  -h, --help             Show this help message"
            echo ""
            echo "Examples:"
            echo "  $0 --quick                    # Quick unit tests"
            echo "  $0 --nprocs 4                 # Run with 4 MPI processes"
            echo "  $0 --filter SinglePhaseFlow*  # Run specific tests"
            echo "  $0 --convergence --verbose    # Convergence tests with output"
            exit 0
            ;;
        *)
            echo -e "${RED}Unknown option: $1${NC}"
            exit 1
            ;;
    esac
done

# Check if build directory exists
if [ ! -d "$BUILD_DIR" ]; then
    echo -e "${RED}Build directory '$BUILD_DIR' not found!${NC}"
    echo "Run: mkdir -p $BUILD_DIR && cd $BUILD_DIR && cmake .. && make"
    exit 1
fi

# Check if test executable exists
if [ ! -f "$BUILD_DIR/tests/run_tests" ]; then
    echo -e "${RED}Test executable not found!${NC}"
    echo "Run: cd $BUILD_DIR && make run_tests"
    exit 1
fi

cd "$BUILD_DIR"

echo -e "${BLUE}================================${NC}"
echo -e "${BLUE}FSRM Test Suite${NC}"
echo -e "${BLUE}================================${NC}"
echo ""
echo "Build directory: $BUILD_DIR"
echo "MPI processes: $NPROCS"

# Set test filter based on options
if [ $QUICK_ONLY -eq 1 ]; then
    TEST_FILTER="UnitTest*:EclipseIO*:PhysicsKernelTestFixture*:WellModel*"
    echo "Test mode: Quick unit tests"
elif [ $CONVERGENCE_ONLY -eq 1 ]; then
    TEST_FILTER="MMSConvergenceTest*"
    echo "Test mode: Convergence tests"
elif [ $INTEGRATION_ONLY -eq 1 ]; then
    TEST_FILTER="IntegrationTest*"
    echo "Test mode: Integration tests"
elif [ -n "$TEST_FILTER" ]; then
    echo "Test filter: $TEST_FILTER"
else
    echo "Test mode: All tests"
fi

echo -e "${BLUE}--------------------------------${NC}"
echo ""

# Build test command
if [ $NPROCS -gt 1 ]; then
    TEST_CMD="mpirun -np $NPROCS ./tests/run_tests"
else
    TEST_CMD="./tests/run_tests"
fi

# Add filter if specified
if [ -n "$TEST_FILTER" ]; then
    TEST_CMD="$TEST_CMD --gtest_filter=$TEST_FILTER"
fi

# Add verbose flag if requested
if [ $VERBOSE -eq 1 ]; then
    TEST_CMD="$TEST_CMD --gtest_print_time=1"
    CTEST_ARGS="-V"
else
    CTEST_ARGS=""
fi

# Run tests
echo -e "${YELLOW}Running tests...${NC}"
echo "Command: $TEST_CMD"
echo ""

START_TIME=$(date +%s)

# Run using CTest for better output formatting
if [ -n "$TEST_FILTER" ]; then
    ctest $CTEST_ARGS -R "${TEST_FILTER//\*/.*}"
    TEST_RESULT=$?
else
    ctest $CTEST_ARGS
    TEST_RESULT=$?
fi

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

echo ""
echo -e "${BLUE}================================${NC}"

if [ $TEST_RESULT -eq 0 ]; then
    echo -e "${GREEN}✓ All tests passed!${NC}"
else
    echo -e "${RED}✗ Some tests failed!${NC}"
fi

echo -e "${BLUE}================================${NC}"
echo "Duration: ${DURATION}s"
echo ""

# Generate summary
if [ -f "Testing/Temporary/LastTest.log" ]; then
    echo -e "${YELLOW}Test Summary:${NC}"
    grep -E "Test #|Passed|Failed|Total Test time" Testing/Temporary/LastTest.log | tail -20
    echo ""
fi

# Check for test report
if [ -f "tests/test_report.txt" ]; then
    echo -e "${YELLOW}Detailed report available at:${NC} $BUILD_DIR/tests/test_report.txt"
fi

# Check for generated plots
PLOT_COUNT=$(find tests -name "*.png" 2>/dev/null | wc -l)
if [ $PLOT_COUNT -gt 0 ]; then
    echo -e "${YELLOW}Generated $PLOT_COUNT plot(s) in:${NC} $BUILD_DIR/tests/"
fi

echo ""

# Return test result code
exit $TEST_RESULT
