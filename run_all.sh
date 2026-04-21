#!/bin/bash
# ============================================================
# run_all.sh -- SolarHPC
# Complete end-to-end pipeline from scratch.
#
# Usage:
#   bash run_all.sh [options]
#
# Options:
#   --skip-fetch    Skip JPL data fetch (use existing data/)
#   --np=N          MPI ranks for simulation (default: 4)
#   --years=N       Simulation years (default: 2)
#   --target=BODY   MOON or MARS (default: MOON)
#   --dry-run       Only validate, no full simulation
#
# Git stages committed automatically after each success.
# ============================================================

set -uo pipefail

# ── Defaults ─────────────────────────────────────────────────
SKIP_FETCH=0
NP=4
YEARS=2
TARGET=MOON
DRY_RUN=0

for arg in "$@"; do
  case "$arg" in
    --skip-fetch)    SKIP_FETCH=1 ;;
    --np=*)          NP="${arg#*=}" ;;
    --years=*)       YEARS="${arg#*=}" ;;
    --target=*)      TARGET="${arg#*=}" ;;
    --dry-run)       DRY_RUN=1 ;;
  esac
done

# ── Helpers ───────────────────────────────────────────────────
START_TIME=$(date +%s)
LOG_FILE="logs/pipeline_$(date +%Y%m%d_%H%M%S).log"
mkdir -p logs

log() {
  local level="$1"; shift
  local msg="$*"
  local ts
  ts=$(date '+%Y-%m-%d %H:%M:%S')
  echo "[$ts] [$level] $msg" | tee -a "$LOG_FILE"
}

stage() { echo ""; log INFO "======== STAGE $1: $2 ========"; }
ok()    { log INFO "[OK] $1"; }
fail()  { log ERROR "[FAIL] $1"; exit 1; }

# ── Set OMP threads ───────────────────────────────────────────
if [ -z "${OMP_NUM_THREADS:-}" ]; then
  export OMP_NUM_THREADS=4
  log WARN "OMP_NUM_THREADS not set. Using default: 4"
fi

stage 0 "Create directories"
make dirs -s && ok "Directories ready"

stage 1 "Data acquisition"
if [ "$SKIP_FETCH" -eq 0 ]; then
  python3 fetch_jpl.py             || fail "JPL fetch failed"
  python3 fetch_reference_ephemeris.py 2>/dev/null || \
    log WARN "Reference ephemeris fetch skipped (optional)"
  ok "Initial conditions written to data/"
  git add data/ 2>/dev/null && \
    git commit -m "Stage 0: initial conditions from JPL Horizons" 2>/dev/null || true
else
  [ -f data/initial_conditions.dat ] || fail "data/initial_conditions.dat missing"
  ok "Skipping fetch (--skip-fetch)"
fi

stage 2 "Build"
make all 2>&1 | tee -a "$LOG_FILE" | tail -3
[ -f solarhpc ] || fail "Build failed -- see $LOG_FILE"
ok "solarhpc binary built"
git add src/ include/ Makefile 2>/dev/null && \
  git commit -m "Stage 1: Fortran + C++ compiled" 2>/dev/null || true

stage 3 "Validation tests"
bash tests/run_tests.sh 2>&1 | tee -a "$LOG_FILE"
ok "Tests passed"
git add results/validation/ tests/ 2>/dev/null && \
  git commit -m "Stage 2: all validation tests passed" 2>/dev/null || true

if [ "$DRY_RUN" -eq 1 ]; then
  log INFO "Dry-run complete. Exiting."
  exit 0
fi

stage 4 "Full simulation (${YEARS} yr, target=${TARGET}, ${NP} MPI ranks)"
mpirun -np "$NP" ./solarhpc \
  --years="$YEARS" \
  --site=ISRO_SRIHARIKOTA \
  --target="$TARGET" \
  --perturb \
  2>&1 | tee -a "$LOG_FILE"
ok "Simulation complete"
git add results/orbits/ results/eclipses/ results/missions/ results/validation/ \
  2>/dev/null && \
  git commit -m "Stage 3: simulation results (${YEARS}yr, ${TARGET})" 2>/dev/null || true

stage 5 "MPI scaling sweep"
for NP_SWEEP in 1 2 4 8; do
  if [ "$NP_SWEEP" -le "$NP" ] || [ "$NP_SWEEP" -eq 1 ]; then
    log INFO "  mpirun -np $NP_SWEEP ./solarhpc --benchmark-only"
    mpirun -np "$NP_SWEEP" ./solarhpc --benchmark-only \
      2>&1 | tee -a "$LOG_FILE" | grep -E "SWEEP|TIMER|MPI" || true
  fi
done
ok "MPI scaling data in results/benchmarks/mpi_scaling.csv"
git add results/benchmarks/ 2>/dev/null && \
  git commit -m "Stage 4: benchmark data (OpenMP + MPI)" 2>/dev/null || true

stage 6 "Package results for download"
bash package_for_download.sh 2>&1 | tee -a "$LOG_FILE"

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
echo ""
log INFO "Pipeline complete in ${ELAPSED}s"
log INFO "Log saved: $LOG_FILE"
log INFO "Next steps (run LOCALLY):"
log INFO "  1. Download: scp <cluster>:<path>/results_*.tar.gz ."
log INFO "  2. Extract:  tar -xzf results_*.tar.gz"
log INFO "  3. Plot:     make plots"
log INFO "  4. Tables:   python3 plots/generate_tables.py"
echo ""
