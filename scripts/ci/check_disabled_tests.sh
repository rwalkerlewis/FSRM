#!/usr/bin/env bash
# =============================================================================
# check_disabled_tests.sh — Report disabled tests from tests/CMakeLists.txt
# =============================================================================
# Parses the test CMakeLists.txt for commented-out test source files and
# reports them grouped by category with their stated reason.
#
# Usage:
#   ./scripts/ci/check_disabled_tests.sh [path/to/tests/CMakeLists.txt]
#
# Output goes to stdout (markdown format).
# If $GITHUB_STEP_SUMMARY is set, output is also appended there.
# =============================================================================

set -euo pipefail

CMAKELISTS="${1:-tests/CMakeLists.txt}"

if [ ! -f "$CMAKELISTS" ]; then
    echo "ERROR: Cannot find $CMAKELISTS"
    exit 1
fi

# Collect disabled test entries grouped by category
total_disabled=0

declare -A categories
declare -A cat_counts

while IFS= read -r line; do
    # Match lines like:
    #   # unit/test_foo.cpp  # Disabled: reason text
    #   # performance/test_bar.cpp  # Disabled: reason
    if [[ "$line" =~ ^[[:space:]]*#[[:space:]]+(([a-z]+)/test_[a-zA-Z0-9_]+\.cpp)[[:space:]]*#[[:space:]]*[Dd]isabled:[[:space:]]*(.*) ]]; then
        filepath="${BASH_REMATCH[1]}"
        category="${BASH_REMATCH[2]}"
        reason="${BASH_REMATCH[3]}"

        if [ -z "${categories[$category]+x}" ]; then
            categories[$category]=""
            cat_counts[$category]=0
        fi

        categories[$category]+="| \`$filepath\` | $reason |"$'\n'
        cat_counts[$category]=$(( ${cat_counts[$category]} + 1 ))
        total_disabled=$(( total_disabled + 1 ))
    fi
done < "$CMAKELISTS"

# ---------------------------------------------------------------------------
# Generate the markdown report
# ---------------------------------------------------------------------------
generate_report() {
    echo "## Disabled Tests Report"
    echo ""

    if [ "$total_disabled" -eq 0 ]; then
        echo "✅ No disabled tests found."
        return
    fi

    echo "⚠️ **$total_disabled** test source files are currently disabled in \`$CMAKELISTS\`."
    echo ""

    for category in $(echo "${!categories[@]}" | tr ' ' '\n' | sort); do
        count="${cat_counts[$category]}"
        # Capitalize first letter of category name
        cap_category="$(echo "${category:0:1}" | tr '[:lower:]' '[:upper:]')${category:1}"
        echo "### $cap_category Tests ($count disabled)"
        echo ""
        echo "| File | Reason |"
        echo "|------|--------|"
        echo -n "${categories[$category]}"
        echo ""
    done

    echo "---"
    echo ""
    echo "*These tests are commented out in the source lists but may still have CTest*"
    echo "*registrations that pass vacuously (0 tests matched = success in GTest).*"
    echo "*To re-enable a test, uncomment its source file and fix the stated issue.*"
}

REPORT="$(generate_report)"

# Print to stdout
echo "$REPORT"

# Append to GitHub Step Summary if available
if [ -n "${GITHUB_STEP_SUMMARY:-}" ]; then
    echo "$REPORT" >> "$GITHUB_STEP_SUMMARY"
fi

exit 0
