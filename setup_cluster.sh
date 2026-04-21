#!/bin/bash
# ============================================================
# setup_cluster.sh -- SolarHPC
# Source this file before building or running:
#   source setup_cluster.sh
# ============================================================

echo "  [SETUP] Loading SolarHPC environment..."

# ── Load HPC modules (edit for your cluster) ─────────────────
# Try common module names — one of these will work on most clusters
module load gcc/12.2.0     2>/dev/null || \
  module load gcc/11.3.0   2>/dev/null || \
  module load gcc           2>/dev/null || \
  echo "  [INFO]  Using system gcc: $(gcc --version 2>&1 | head -1)"

module load openmpi/4.1.4  2>/dev/null || \
  module load openmpi/4.1   2>/dev/null || \
  module load openmpi        2>/dev/null || \
  echo "  [INFO]  Using system MPI: $(mpicxx --version 2>&1 | head -1)"

module load python/3.10.0  2>/dev/null || \
  module load python/3.9    2>/dev/null || \
  module load python3        2>/dev/null || \
  echo "  [INFO]  Using system Python: $(python3 --version 2>&1)"

# ── Environment variables ─────────────────────────────────────
# Set OMP_NUM_THREADS to number of cores per node
# Check: nproc --all
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-8}"
export OMP_PROC_BIND=close
export OMP_PLACES=cores

# Compilers
export FC=gfortran
export CXX=g++
export MPICC=mpicxx

# ── Verify tools ──────────────────────────────────────────────
echo ""
for cmd in gfortran g++ mpicxx mpirun python3; do
  if command -v "$cmd" &>/dev/null; then
    printf "  [OK]   %-10s  %s\n" "$cmd" "$(which $cmd)"
  else
    printf "  [MISS] %-10s  NOT FOUND -- load the correct module\n" "$cmd"
  fi
done

# ── Python packages ───────────────────────────────────────────
echo ""
python3 -c "import numpy, matplotlib, pandas, requests, scipy; print('  [OK]   Python: numpy matplotlib pandas requests scipy')" 2>/dev/null || \
  echo "  [WARN] Missing Python packages. Run:"                                && \
  echo "         pip install numpy matplotlib pandas requests scipy --user"

# ── Summary ───────────────────────────────────────────────────
echo ""
echo "  OMP_NUM_THREADS = $OMP_NUM_THREADS"
echo ""
echo "  Quick reference:"
echo "  ─────────────────────────────────────────────"
echo "  make dirs           Create output directories"
echo "  make fetch          Fetch JPL data (internet)"
echo "  make all            Compile everything"
echo "  make run-tests      Run validation tests"
echo "  make run NP=4       Run with 4 MPI ranks"
echo "  make benchmark      OpenMP thread sweep"
echo "  make sweep-mpi      MPI scaling experiment"
echo "  bash run_all.sh     Full pipeline"
echo "  ─────────────────────────────────────────────"
