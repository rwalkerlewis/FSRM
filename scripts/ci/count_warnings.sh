#!/usr/bin/env bash
# =============================================================================
# count_warnings.sh — Count and report compiler warnings from a build log
# =============================================================================
# Parses a GCC/G++ build log for compiler warnings and produces a summary
# showing total count and the most common warning types.
#
# Usage:
#   ./scripts/ci/count_warnings.sh <build-log-file>
#
# Output goes to stdout (markdown format).
# If $GITHUB_STEP_SUMMARY is set, output is also appended there.
# Exit code is always 0 (informational only — does not fail the build).
# =============================================================================

set -uo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: $0 <build-log-file>"
    exit 1
fi

BUILD_LOG="$1"

if [ ! -f "$BUILD_LOG" ]; then
    echo "ERROR: Build log file not found: $BUILD_LOG"
    exit 1
fi

# ---------------------------------------------------------------------------
# Count warnings
# ---------------------------------------------------------------------------
# GCC/G++ warning format: "file.cpp:line:col: warning: message [-Wflag]"
TOTAL_WARNINGS=$(grep -c ': warning:' "$BUILD_LOG" || true)
# Ensure it's a clean integer
TOTAL_WARNINGS=${TOTAL_WARNINGS:-0}
TOTAL_WARNINGS=$((TOTAL_WARNINGS + 0))

# Extract warning types (the [-Wflag] part)
WARNING_TYPES=$(grep -oP '\[-W[^\]]+\]' "$BUILD_LOG" | \
    sort | uniq -c | sort -rn | head -10 || true)

# Count unique files with warnings
WARNING_FILES=$(grep ': warning:' "$BUILD_LOG" | \
    grep -oP '^[^:]+' | sort -u | wc -l || true)
WARNING_FILES=${WARNING_FILES:-0}
WARNING_FILES=$((WARNING_FILES + 0))

# ---------------------------------------------------------------------------
# Generate the markdown report
# ---------------------------------------------------------------------------
generate_report() {
    echo "## Compiler Warning Report"
    echo ""

    if [ "$TOTAL_WARNINGS" -eq 0 ]; then
        echo "✅ **0 warnings** — clean build!"
        return
    fi

    echo "⚠️ **$TOTAL_WARNINGS warnings** found across **$WARNING_FILES files**."
    echo ""

    if [ -n "$WARNING_TYPES" ]; then
        echo "### Top Warning Types"
        echo ""
        echo "| Count | Warning Flag |"
        echo "|------:|-------------|"
        while IFS= read -r line; do
            if [ -n "$line" ]; then
                count=$(echo "$line" | awk '{print $1}')
                flag=$(echo "$line" | awk '{print $2}')
                echo "| $count | \`$flag\` |"
            fi
        done <<< "$WARNING_TYPES"
        echo ""
    fi

    echo "*Warnings are informational and do not fail the build.*"
    echo "*To enable warning-as-error mode, pass \`-DCMAKE_CXX_FLAGS=\"-Werror\"\` to cmake.*"
}

REPORT="$(generate_report)"

# Print to stdout
echo "$REPORT"

# Append to GitHub Step Summary if available
if [ -n "${GITHUB_STEP_SUMMARY:-}" ]; then
    echo "$REPORT" >> "$GITHUB_STEP_SUMMARY"
fi

exit 0
