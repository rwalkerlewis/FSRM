#!/bin/bash
# tests/diagnostics/run_all_baselines.sh
#
# Pass-12 followup 3 convenience driver: regenerate every parallel-KSP
# baseline / post-fix log in one go. Iterates the twelve "SNES diverges
# step 0" examples and calls capture_parallel_ksp_baseline.sh for each
# with the same FSRM_FINAL_TIME_OVERRIDE the examples_runtime gate uses.
#
# Usage (from the repo root, inside fsrm-ci):
#   bash tests/diagnostics/run_all_baselines.sh

set -u
REPO_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$REPO_DIR"

declare -A FINAL_TIME
FINAL_TIME[06_gmsh_multimaterial]=0.001
FINAL_TIME[09_gasbuggy_1967]=0.05
FINAL_TIME[17_velocity_model]=0.01
FINAL_TIME[23_milrow_1969]=0.05
FINAL_TIME[24_cannikin_1971]=0.05
FINAL_TIME[25_faultless_1968]=0.05
FINAL_TIME[29_rio_blanco_1973]=0.05
FINAL_TIME[33_dprk_2006]=0.05
FINAL_TIME[34_dprk_2009]=0.05
FINAL_TIME[35_dprk_2013]=0.05
FINAL_TIME[36_dprk_2016a]=0.05
FINAL_TIME[37_dprk_2016b]=0.05

EXAMPLES=(
    06_gmsh_multimaterial
    09_gasbuggy_1967
    17_velocity_model
    23_milrow_1969
    24_cannikin_1971
    25_faultless_1968
    29_rio_blanco_1973
    33_dprk_2006
    34_dprk_2009
    35_dprk_2013
    36_dprk_2016a
    37_dprk_2016b
)

for ex in "${EXAMPLES[@]}"; do
    ft="${FINAL_TIME[$ex]:-0.05}"
    echo "###### capture $ex (final_time=$ft) ######"
    bash "$REPO_DIR/tests/diagnostics/capture_parallel_ksp_baseline.sh" \
         "$REPO_DIR/examples/$ex" "$ft" \
        || echo "[run_all_baselines] capture script exited non-zero for $ex"
done
echo "###### run_all_baselines done ######"
