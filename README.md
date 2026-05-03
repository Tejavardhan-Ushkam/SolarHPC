# SolarHPC — Parallel N-body Solar System Simulator

[![Build](https://img.shields.io/badge/build-passing-brightgreen)]()
[![Language](https://img.shields.io/badge/language-Fortran%20%7C%20C%2B%2B%20%7C%20Python-blue)]()
[![Parallel](https://img.shields.io/badge/parallel-OpenMP%20%7C%20MPI-orange)]()

A high-performance parallel simulator for solar system orbital
dynamics, eclipse prediction, and spacecraft trajectory
optimisation. 

---

## What it does

- **Integrates 10-body solar system** (Sun + 8 planets + Moon)
  using a symplectic leapfrog integrator over multi-year horizons
- **Predicts solar and lunar eclipses** from first-principles
  shadow-cone geometry, validated against the NASA eclipse catalogue
- **Optimises spacecraft trajectories** (Hohmann transfer or
  gravity assist) to the Moon or Mars from any of 5 real launch
  sites, returning Δv budget, fuel mass, transfer time, and
  predicted landing coordinates
- **Benchmarks parallel performance** — 5 OpenMP kernels across
  {1,2,4,8,16} threads and 1 MPI kernel across {1,2,4,8} ranks,
  with Amdahl-law fitting and publication-quality figures

---

## Quick start

```bash
git clone <repo-url>
cd SolarHPC

# 1. Load cluster environment (edit module names for your cluster)
source setup_cluster.sh

# 2. Create output directories
make dirs

# 3. Fetch initial conditions from NASA JPL Horizons (~30 seconds, needs internet)
make fetch          # writes data/initial_conditions.dat
make fetch-ref      # writes data/reference_ephemeris.csv (for validation)

# 4. Compile everything
make all            # builds solarhpc binary

# 5. Run validation tests
make run-tests      # 5 tests: 3 physics + 2 simulation

# 6. Run 2-year Moon simulation with 4 MPI ranks
make run NP=4 ARGS="--years=2 --target=MOON --site=ISRO_SRIHARIKOTA"

# 7. Run benchmark sweeps
make benchmark      # OpenMP thread sweep (all 5 kernels)
make sweep-mpi      # MPI scaling (1,2,4,8 ranks)
```

---

## Run commands reference

| Command | Purpose | Expected output | Time |
|---------|---------|----------------|------|
| `make dirs` | Create result folders | Directories | instant |
| `make fetch` | Query JPL Horizons API | `data/initial_conditions.dat` | ~30s |
| `make all` | Compile Fortran + C++ | `solarhpc` binary | ~20s |
| `make run-tests` | Physics unit tests | 5 PASS | ~5 min |
| `make run NP=4` | Full 2-year simulation | results/ CSVs | ~10 min |
| `make benchmark` | OpenMP thread sweep | benchmarks CSV | ~15 min |
| `make sweep-mpi` | MPI scaling study | mpi_scaling CSV | ~20 min |
| `bash run_all.sh` | Complete pipeline | everything | ~1 hr |
| `make plots` | Generate all figures | 14 PDFs | ~2 min |

---

## Project structure

```
SolarHPC/
├── src/
│   ├── fortran/              Numerical kernels (OpenMP)
│   │   ├── solar_constants.f90    Physical constants
│   │   ├── omp_timer.f90          Timing + statistics
│   │   ├── nbody_force.f90        O(N²) force computation
│   │   ├── integrators.f90        Leapfrog + RK4
│   │   ├── coordinate_transform.f90
│   │   ├── eclipse_geom.f90       Shadow-cone eclipse detection
│   │   ├── trajectory_prop.f90    Spacecraft trajectories
│   │   └── validation_suite.f90   Unit tests
│   └── cpp/                  Orchestration (MPI)
│       ├── main.cpp               Pipeline driver
│       ├── jpl_fetch.cpp          Load initial conditions
│       ├── mission_config.cpp     Parse command-line config
│       ├── results_writer.cpp     Centralised CSV output
│       └── mpi_launcher.cpp       MPI launch-window search
├── include/
│   ├── solar_types.h              C++ structs and enums
│   └── fortran_iface.h            C→Fortran bindings
├── tests/                    Standalone validation programs
├── plots/                    Python plotting scripts
├── results/                  Simulation output (gitignored if large)
├── docs/paper/               LaTeX research paper
├── data/                     JPL initial conditions (gitignored)
├── Makefile
├── fetch_jpl.py
├── fetch_reference_ephemeris.py
├── setup_cluster.sh
├── run_all.sh
└── package_for_download.sh
```

---

## Mission configuration

All flags optional — sensible defaults apply:

```bash
mpirun -np 4 ./solarhpc \
  --site=ISRO_SRIHARIKOTA    # or NASA_KENNEDY, ROSCOSMOS_BAIKONUR,
                              #    ESA_KOUROU, JAXA_TANEGASHIMA
  --target=MOON              # or MARS
  --mode=HOHMANN             # or GRAVITY_ASSIST
  --mass=5000                # spacecraft net mass (kg)
  --payload=1000             # payload mass (kg)
  --years=2                  # simulation duration
  --perturb                  # enable J2 + GR corrections
```

```bash



## Physical validation results

| Test | Expected | Tolerance | Status |
|------|----------|-----------|--------|
| Two-body circular orbit (1 yr) | Return to start | < 200,000 km | PASS |
| Hohmann Earth→Mars Δv | 5591 m/s | < 1% | PASS |
| Tsiolkovsky fuel (dv=5590, Isp=450) | 2549 kg | < 0.01% | PASS |
| Energy conservation (1 yr leapfrog) | Bounded drift | < 0.01% | PASS |
| JD conversion round-trip | J2000 = 2451545.0 | < 1e-6 days | PASS |
| Eclipse geometry (aligned) | event_flag ≥ 1 | detected | PASS |
| Eclipse 2024-Apr-08 | JD 2460409.5 | ± 2 days | PASS* |

*Requires `--years=24` to cover 2024.

---

