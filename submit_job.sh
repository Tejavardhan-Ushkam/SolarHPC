#!/bin/bash
#SBATCH --job-name=SolarHPC
#SBATCH --partition=day
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --output=logs/slurm_%j.log
#SBATCH --error=logs/slurm_%j.err

mkdir -p logs

source setup_cluster.sh
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo "Job $SLURM_JOB_ID started: $(date)"
echo "Nodes: $SLURM_JOB_NODELIST"
echo "MPI tasks: $SLURM_NTASKS  OMP threads: $OMP_NUM_THREADS"

# Build if needed
[ ! -f solarhpc ] && make all 2>&1 | tee logs/build_${SLURM_JOB_ID}.log

# Validate data exists
if [ ! -f ./data/initial_conditions.dat ]; then
    echo "[ERROR] Missing data file"
    ls -l ./data/
    exit 1
fi

cd /wd/users/b23233/HPSC_Solar/SolarHPC
# Run tests
bash tests/run_tests.sh || { echo "[ERROR] Tests failed"; exit 2; }

# Full simulation
srun ./solarhpc \
    --years=2 --target=MOON \
    2>&1 | tee logs/sim_${SLURM_JOB_ID}.log

# Benchmark sweep
for NP in 1 2; do
    srun -n $NP ./solarhpc --benchmark-only \
        2>&1 | tee -a logs/bench_${SLURM_JOB_ID}.log
done

# Package results
bash package_for_download.sh

echo "Job completed: $(date)"
