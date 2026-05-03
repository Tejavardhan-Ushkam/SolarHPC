#!/bin/bash
# ============================================================
# run_tests.sh -- SolarHPC test runner
# Usage: bash tests/run_tests.sh
# ============================================================

PASS=0; FAIL=0; SKIP=0
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"

run_test() {
  local name="$1"
  local exe="$2"
  local needs_data="${3:-no}"

  if [ "$needs_data" = "yes" ] && [ ! -f "data/initial_conditions.dat" ]; then
    echo "  SKIP  $name (run 'make fetch' first)"
    SKIP=$((SKIP+1))
    return
  fi

  printf "  %-45s " "$name ..."
  if bash -c "$exe" > /tmp/solarhpc_test_out.txt 2>&1; then
    echo "PASS"
    PASS=$((PASS+1))
  else
    echo "FAIL"
    grep -E "FAIL|ERROR|error" /tmp/solarhpc_test_out.txt 2>/dev/null | head -3
    FAIL=$((FAIL+1))
  fi
}

echo ""
echo "========================================================"
echo "  SolarHPC Validation Test Suite"
echo "========================================================"
echo ""

# Build tests if needed (suppress warnings, allow success with warnings)
[ ! -f tests/test_two_body ] && { make tests 2>&1 | grep -v "Warning\|warning" || true; }

run_test "Force law unit test"         tests/test_two_body
run_test "Hohmann + Tsiolkovsky"       tests/test_hohmann
run_test "Eclipse geometry + JD"       tests/test_eclipse
run_test "Energy conservation (10yr)"  tests/test_energy    "yes"

if [ -f "./solarhpc" ]; then
  run_test "Simulation dry-run" "cd .. && ./solarhpc --dry-run" "yes"
else
  echo "  SKIP  Simulation dry-run (run 'make all' first)"
  SKIP=$((SKIP+1))
fi

echo ""
echo "========================================================"
printf "  Results: %d passed, %d failed, %d skipped\n" $PASS $FAIL $SKIP
echo "========================================================"
echo ""

if [ "$FAIL" -gt 0 ]; then
  echo "  ACTION: Fix failing tests before running full simulation."
  exit 1
fi
if [ "$SKIP" -gt 0 ]; then
  echo "  NOTE: Skipped tests require data. Run: make fetch"
fi
echo "  All runnable tests passed. Safe to proceed."
echo ""
exit 0
