# SolarHPC — Parallel N-body Solar System Simulator

A high-performance computing project implementing parallel
gravitational N-body simulation for eclipse prediction and
spacecraft trajectory optimisation.
Uses **Fortran 90 + OpenMP** for compute kernels,
**C++17 + MPI** for orchestration, and **Python 3** for
publication-quality figures.

---

## What it does

- Integrates orbits of all 10 major solar system bodies for
  arbitrary duration using a symplectic leapfrog integrator
- Predicts solar and lunar eclipses from shadow-cone geometry
- Plans spacecraft missions to the Moon and Mars (Hohmann
  transfer and gravity-assist), computing Δv, fuel mass, and
  landing coordinates for any chosen launch site
- Benchmarks 5 parallel kernels across 1–16 OpenMP threads
  and 1–8 MPI ranks, producing speedup/efficiency data
  suitable for a research paper

---

## Requirements

| Tool | Minimum version |
|---|---|
| gfortran | 10 |
| g++ | 10 |
| OpenMPI | 4.0 |
| Python | 3.8 |
| numpy, matplotlib, pandas, requests, scipy | any recent |

Install Python packages:
```bash
pip install numpy matplotlib pandas requests scipy --user
```

---

## Quick start

```bash
git clone <repo-url>
cd SolarHPC

# On a cluster: load modules first
source setup_cluster.sh

# Create output directories
make dirs

# Fetch initial conditions from NASA JPL (needs internet, ~30s)
make fetch
make fetch-ref     # reference ephemeris for validation (optional)

# Compile everything
make all

# Run validation tests
make run-tests     # 3 physics tests pass immediately
                   # 2 more pass after make fetch

# Run the full simulation (4 MPI ranks, 2-year integration, Moon mission)
mpirun -np 4 ./solarhpc --years=2 --target=MOON

# Run benchmark sweep
make benchmark     # OpenMP thread sweep
make sweep-mpi     # MPI scaling (1,2,4,8 ranks)

# Download results and generate plots LOCALLY
tar -czf results.tar.gz results/ logs/
scp results.tar.gz localPC:~/SolarHPC/
# then on local machine:
tar -xzf results.tar.gz
make plots
```

---

## All make targets

```
make all           Build solarhpc binary
make tests         Build all test executables
make run-tests     Build + run validation tests
make fetch         Fetch initial conditions from JPL Horizons API
make fetch-ref     Fetch reference ephemeris for validation
make validate      Run simulation with --dry-run
make run NP=4      Run simulation (4 MPI ranks)
make benchmark     OpenMP thread sweep (all 5 kernels)
make sweep-mpi     MPI scaling sweep (1,2,4,8 ranks)
make plots         Generate all figures (run locally)
make dirs          Create output directories
make clean         Remove build artefacts
make distclean     Remove everything including results
make help          Show all targets with descriptions
```

---

## Command-line flags

```
./solarhpc [options]

--site=NAME         Launch site (default: ISRO_SRIHARIKOTA)
                    Options: ISRO_SRIHARIKOTA, NASA_KENNEDY,
                             ROSCOSMOS_BAIKONUR, ESA_KOUROU,
                             JAXA_TANEGASHIMA
--target=BODY       MOON or MARS (default: MOON)
--mode=MODE         HOHMANN or GRAVITY_ASSIST (default: HOHMANN)
--mass=KG           Spacecraft wet mass (default: 5000)
--payload=KG        Payload / dry mass (default: 1000)
--years=N           Simulation duration in years (default: 2)
--dt=DAYS           Integration timestep in days (default: 1)
--perturb           Enable J2, solar radiation, GR corrections
--benchmark-only    Skip simulation, run thread sweep only
--dry-run           Validate only, no orbital integration
```

Examples:
```bash
# Moon mission from Kennedy Space Center
mpirun -np 4 ./solarhpc --site=NASA_KENNEDY --target=MOON --years=1

# Mars mission with gravity assist from ISRO
mpirun -np 8 ./solarhpc --site=ISRO_SRIHARIKOTA --target=MARS \
  --mode=GRAVITY_ASSIST --years=3 --perturb

# Benchmark only (no simulation output)
mpirun -np 1 ./solarhpc --benchmark-only
```

---

## Output files

| Path | Contents |
|---|---|
| `results/orbits/orbit_trajectories.csv` | Position/velocity per body per timestep |
| `results/eclipses/eclipse_predictions.csv` | Eclipse events with date, type, umbra fraction |
| `results/missions/mission_report.csv` | Δv, fuel, landing coordinates |
| `results/missions/launch_window_sweep.csv` | Δv vs launch date (MPI output) |
| `results/benchmarks/thread_sweep_results.csv` | Timing per kernel per thread count |
| `results/benchmarks/mpi_scaling.csv` | Wall time and speedup per rank count |
| `results/benchmarks/integrator_comparison.csv` | Energy drift per year (leapfrog vs RK4) |
| `results/validation/test_report.txt` | PASS/FAIL for all unit tests |
| `results/run_metadata.txt` | Cluster info, compiler flags, Git hash |

---

## Reproducing paper figures

After downloading `results/` from the cluster:

```bash
# All figures at once
make plots

# Individual scripts
python3 plots/plot_orbits.py           # fig_orbits_2d, fig_orbits_3d
python3 plots/plot_benchmarks.py       # fig_speedup_all_kernels, fig_efficiency, ...
python3 plots/plot_eclipse_timeline.py # fig_eclipse_timeline, fig_saros_pattern
python3 plots/plot_validation.py       # fig_energy_drift, fig_ephemeris_error
python3 plots/plot_mission.py          # fig_launch_window, fig_landing_zone, ...
python3 plots/generate_tables.py       # LaTeX tables -> docs/paper/tables/

# Compile paper (requires LaTeX)
cd docs/paper
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```

---

## Git workflow (commit after each successful stage)

```bash
# Stage 0: after make fetch
git add data/
git commit -m "Stage 0: initial conditions from JPL Horizons"

# Stage 1: after make all
git add src/ include/ Makefile
git commit -m "Stage 1: Fortran + C++ compiled"

# Stage 2: after make run-tests
git add results/validation/ tests/
git commit -m "Stage 2: all validation tests passed"

# Stage 3: after make run
git add results/orbits/ results/eclipses/ results/missions/ results/validation/
git commit -m "Stage 3: 2-year simulation results"

# Stage 4: after make benchmark && make sweep-mpi
git add results/benchmarks/
git commit -m "Stage 4: OpenMP + MPI benchmark data"

# Stage 5: after make plots
git add plots/output/ docs/paper/tables/
git commit -m "Stage 5: publication figures and LaTeX tables"

# Stage 6: after paper draft
git add docs/paper/
git commit -m "Stage 6: paper draft complete"
```

---

## Project structure

```
SolarHPC/
├── src/
│   ├── fortran/           Compute kernels (Fortran 90 + OpenMP)
│   │   ├── solar_constants.f90
│   │   ├── omp_timer.f90
│   │   ├── nbody_force.f90
│   │   ├── integrators.f90
│   │   ├── coordinate_transform.f90
│   │   ├── eclipse_geom.f90
│   │   ├── trajectory_prop.f90
│   │   └── validation_suite.f90
│   └── cpp/               Orchestration (C++17 + MPI)
│       ├── main.cpp
│       ├── jpl_fetch.cpp
│       ├── mission_config.cpp
│       ├── results_writer.cpp
│       └── mpi_launcher.cpp
├── include/               Headers
│   ├── solar_types.h
│   └── fortran_iface.h
├── tests/                 Standalone test programs
├── plots/                 Python plotting scripts
├── results/               Simulation output (gitignored mostly)
├── docs/paper/            IEEE paper (LaTeX)
├── data/                  JPL initial conditions (gitignored)
├── Makefile
├── fetch_jpl.py
├── fetch_reference_ephemeris.py
├── setup_cluster.sh
├── run_all.sh
└── package_for_download.sh
```

---

## Citation

If you use this code, please cite:

```bibtex
@inproceedings{solarhpc2025,
  author    = {[Author Name]},
  title     = {Parallel {N}-body Solar System Simulation for
               Eclipse Prediction and Spacecraft Trajectory
               Optimisation: an {OpenMP/MPI} Performance Study},
  booktitle = {[Conference Name]},
  year      = {2025}
}
```

---

## Physics references

- Meeus, J. *Astronomical Algorithms*, Willmann-Bell, 1998
- Hairer et al. *Geometric Numerical Integration*, Springer, 2006
- Bate, Mueller, White. *Fundamentals of Astrodynamics*, Dover, 1971
- Standish. *JPL Planetary Ephemerides DE405*, NASA TDB, 1998
