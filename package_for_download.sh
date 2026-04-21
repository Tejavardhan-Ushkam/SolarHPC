#!/bin/bash
# ============================================================
# package_for_download.sh -- SolarHPC
# Creates a download archive of results and logs.
# Large orbit trajectory files are excluded if > 50 MB.
# ============================================================

JOB_ID="${SLURM_JOB_ID:-$(date +%Y%m%d_%H%M%S)}"
PKG="results_${JOB_ID}.tar.gz"

EXCLUDE_ARGS=""
if [ -d results/orbits ]; then
  ORBIT_MB=$(du -sm results/orbits 2>/dev/null | cut -f1 || echo "0")
  if [ "$ORBIT_MB" -gt 50 ]; then
    echo "  [WARN] results/orbits/ is ${ORBIT_MB}MB -- excluding from package."
    echo "         Download separately: scp cluster:path/results/orbits/*.csv ."
    EXCLUDE_ARGS="--exclude=results/orbits"
  fi
fi

tar -czf "$PKG" \
  $EXCLUDE_ARGS \
  results/eclipses/ \
  results/missions/ \
  results/benchmarks/ \
  results/validation/ \
  results/run_metadata.txt \
  logs/ \
  2>/dev/null || true

SIZE=$(du -sh "$PKG" 2>/dev/null | cut -f1 || echo "?")
echo ""
echo "  [OK] Package: $PKG  ($SIZE)"
echo ""
echo "  Download with:"
echo "    scp <user>@<cluster>:<full_path>/$PKG ."
echo ""
echo "  Then locally:"
echo "    tar -xzf $PKG"
echo "    make plots          # generates all figures"
echo "    python3 plots/generate_tables.py  # LaTeX tables"
echo ""
