#!/usr/bin/env bash
# count_warnings.sh — Count compiler warnings from a build log.
#
# Usage:
#   ./scripts/ci/count_warnings.sh build.log
#   make -j$(nproc) 2>&1 | tee build.log && ./scripts/ci/count_warnings.sh build.log
#
# Exit code is always 0 (informational, not gating).

set -euo pipefail

LOG_FILE="${1:--}"  # default to stdin if no argument

if [[ "$LOG_FILE" == "-" ]]; then
    INPUT=$(cat)
else
    if [[ ! -f "$LOG_FILE" ]]; then
        echo "Warning count: 0 (log file not found: $LOG_FILE)"
        exit 0
    fi
    INPUT=$(cat "$LOG_FILE")
fi

# Count unique warning lines (GCC/Clang format: "file:line:col: warning: ...")
WARNING_COUNT=$(echo "$INPUT" | grep -c ' warning: ' || true)

echo "============================================"
echo "  Compiler Warning Summary"
echo "============================================"
echo "  Total warnings: ${WARNING_COUNT}"
echo "============================================"

if [[ "$WARNING_COUNT" -gt 0 ]]; then
    echo ""
    echo "Top warning types:"
    echo "$INPUT" \
        | grep ' warning: ' \
        | sed 's/.*warning: //' \
        | sed 's/ \[.*//' \
        | sort \
        | uniq -c \
        | sort -rn \
        | head -20
fi

exit 0
