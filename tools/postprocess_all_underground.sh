#!/bin/bash
# Batch-render the Sedan-style figure pack for every underground-test example
# that already has a solution.h5 on disk. Skips examples that have not been run.
#
# Usage:
#   tools/postprocess_all_underground.sh
set -e

REPO_DIR="$(cd "$(dirname "$0")/.." && pwd)"
TOOLS="${REPO_DIR}/tools"

UNDERGROUND=(
    05_punggye_ri_nuclear_test
    09_gasbuggy_1967
    10_gnome_1961
    11_sedan_1962
    12_degelen_mountain
    13_nts_pahute_mesa
    19_rainier_1957
    20_salmon_1964
    21_sterling_1966
    22_long_shot_1965
    23_milrow_1969
    24_cannikin_1971
    25_faultless_1968
    26_baneberry_1970
    27_schooner_1968
    28_rulison_1969
    29_rio_blanco_1973
    30_chagan_1965
    31_azgir_a1_1966
    32_pokhran_i_1974
    33_dprk_2006
    34_dprk_2009
    35_dprk_2013
    36_dprk_2016a
    37_dprk_2016b
    38_lop_nor_1976
    39_dprk_2017
    40_lop_nor_1996
    41_pokhran_ii_1998
)

RENDERED=0
SKIPPED=0
for d in "${UNDERGROUND[@]}"; do
    EX="${REPO_DIR}/examples/${d}"
    # Pick config: prefer config.config, then any config_*.config
    CFG=""
    if [ -f "${EX}/config.config" ]; then
        CFG="${EX}/config.config"
    else
        CFG=$(ls "${EX}"/config*.config 2>/dev/null | head -1)
    fi
    if [ -z "$CFG" ]; then
        echo "SKIP ${d}: no config"; SKIPPED=$((SKIPPED+1)); continue
    fi
    # Pick output dir: output/ then output_lagrangian/ then output_paraview/
    for OUT in output output_lagrangian output_paraview output_dynamic_fast; do
        if [ -f "${EX}/${OUT}/solution.h5" ]; then
            echo ""
            echo "=== ${d}  (out=${OUT}) ==="
            bash "${TOOLS}/postprocess_example.sh" "${EX}" "${CFG}" "${OUT}" "${d}" || true
            RENDERED=$((RENDERED+1))
            break
        fi
    done
    if [ $? -ne 0 ]; then SKIPPED=$((SKIPPED+1)); fi
done

echo ""
echo "Rendered: ${RENDERED}    Skipped: ${SKIPPED}"
