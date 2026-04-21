# ============================================================
# Makefile -- SolarHPC Parallel N-body Solar System Simulator
# ============================================================
# Dependency order for Fortran modules is critical.
# solar_constants -> omp_timer -> nbody_force -> integrators
# -> coordinate_transform -> eclipse_geom -> trajectory_prop
# -> validation_suite
# ============================================================

FC      = gfortran
CXX     = g++
MPICC   = mpicxx

FFLAGS  = -O3 -fopenmp -fdefault-real-8 -Wall -Wextra \
          -Wno-unused-dummy-argument -fcheck=all
CXXFLAGS= -O3 -fopenmp -std=c++17 -Wall -Wextra -Iinclude
GFORT_LIB := $(shell gfortran -print-file-name=libgfortran.a)
QUAD_LIB  := $(shell gfortran -print-file-name=libquadmath.a)
LDFLAGS    = -fopenmp -Wl,--start-group $(GFORT_LIB) $(QUAD_LIB) -Wl,--end-group -lgomp -lm

# Module output directory (current dir for simplicity)
MODDIR  = .
FFLAGS += -J$(MODDIR)

# ── Fortran sources (COMPILE ORDER MATTERS) ──────────────────
FSRCS = src/fortran/solar_constants.f90 \
        src/fortran/omp_timer.f90       \
        src/fortran/nbody_force.f90     \
        src/fortran/integrators.f90     \
        src/fortran/coordinate_transform.f90 \
        src/fortran/eclipse_geom.f90    \
        src/fortran/trajectory_prop.f90 \
        src/fortran/validation_suite.f90

FOBJS = $(FSRCS:.f90=.o)

# ── C++ sources ──────────────────────────────────────────────
CXXSRCS = src/cpp/jpl_fetch.cpp      \
          src/cpp/mission_config.cpp  \
          src/cpp/results_writer.cpp  \
          src/cpp/mpi_launcher.cpp    \
          src/cpp/main.cpp

CXXOBJS = $(CXXSRCS:.cpp=.o)

# ── Test executables ─────────────────────────────────────────
TEST_SRCS = tests/test_two_body.f90 \
            tests/test_energy.f90   \
            tests/test_eclipse.f90  \
            tests/test_hohmann.f90

# ── Primary target ───────────────────────────────────────────
.PHONY: all tests run-tests fetch fetch-ref validate run benchmark \
        sweep-mpi plots dirs clean distclean help

all: dirs solarhpc

# ── Link final binary ────────────────────────────────────────
solarhpc: $(FOBJS) $(CXXOBJS)
	$(MPICC) $(CXXOBJS) $(FOBJS) $(LDFLAGS) -o $@
	@echo ""
	@echo "  [OK] Built: solarhpc"
	@echo "  Run: mpirun -np 4 ./solarhpc --years=2 --target=MOON"
	@echo ""

# ── Fortran compilation (sequential, order preserved by FSRCS) ──
# Pattern rule cannot guarantee order, so we use explicit rules
# for the first few then fall back to pattern.

src/fortran/solar_constants.o: src/fortran/solar_constants.f90
	$(FC) $(FFLAGS) -c $< -o $@

src/fortran/omp_timer.o: src/fortran/omp_timer.f90 \
    src/fortran/solar_constants.o
	$(FC) $(FFLAGS) -c $< -o $@

src/fortran/nbody_force.o: src/fortran/nbody_force.f90 \
    src/fortran/omp_timer.o
	$(FC) $(FFLAGS) -c $< -o $@

src/fortran/integrators.o: src/fortran/integrators.f90 \
    src/fortran/nbody_force.o
	$(FC) $(FFLAGS) -c $< -o $@

src/fortran/coordinate_transform.o: src/fortran/coordinate_transform.f90 \
    src/fortran/solar_constants.o
	$(FC) $(FFLAGS) -c $< -o $@

src/fortran/eclipse_geom.o: src/fortran/eclipse_geom.f90 \
    src/fortran/integrators.o src/fortran/coordinate_transform.o
	$(FC) $(FFLAGS) -c $< -o $@

src/fortran/trajectory_prop.o: src/fortran/trajectory_prop.f90 \
    src/fortran/eclipse_geom.o
	$(FC) $(FFLAGS) -c $< -o $@

src/fortran/validation_suite.o: src/fortran/validation_suite.f90 \
    src/fortran/trajectory_prop.o
	$(FC) $(FFLAGS) -c $< -o $@

# ── C++ compilation ───────────────────────────────────────────
src/cpp/%.o: src/cpp/%.cpp
	$(MPICC) $(CXXFLAGS) -c $< -o $@

# ── Test executables ─────────────────────────────────────────
FMOD_OBJS = src/fortran/solar_constants.o \
            src/fortran/omp_timer.o       \
            src/fortran/nbody_force.o     \
            src/fortran/integrators.o     \
            src/fortran/coordinate_transform.o \
            src/fortran/eclipse_geom.o    \
            src/fortran/trajectory_prop.o

tests/test_two_body: tests/test_two_body.f90 $(FMOD_OBJS)
	$(FC) $(FFLAGS) $(FMOD_OBJS) $< -o $@

tests/test_energy: tests/test_energy.f90 $(FMOD_OBJS)
	$(FC) $(FFLAGS) $(FMOD_OBJS) $< -o $@

tests/test_eclipse: tests/test_eclipse.f90 $(FMOD_OBJS)
	$(FC) $(FFLAGS) $(FMOD_OBJS) $< -o $@

tests/test_hohmann: tests/test_hohmann.f90 \
    src/fortran/solar_constants.o src/fortran/omp_timer.o \
    src/fortran/nbody_force.o src/fortran/integrators.o \
    src/fortran/coordinate_transform.o src/fortran/trajectory_prop.o
	$(FC) $(FFLAGS) src/fortran/solar_constants.o src/fortran/omp_timer.o \
	  src/fortran/nbody_force.o src/fortran/integrators.o \
	  src/fortran/coordinate_transform.o src/fortran/trajectory_prop.o \
	  $< -o $@

tests: dirs tests/test_two_body tests/test_energy tests/test_eclipse tests/test_hohmann
	@echo "  [OK] All test executables built"

# ── Recipes ──────────────────────────────────────────────────
run-tests: tests
	@bash tests/run_tests.sh

fetch:
	@python3 fetch_jpl.py

fetch-ref:
	@python3 fetch_reference_ephemeris.py

validate: solarhpc
	@mpirun -np 1 ./solarhpc --dry-run

NP ?= 4
run: solarhpc
	mpirun -np $(NP) ./solarhpc $(ARGS)

benchmark: solarhpc
	@export OMP_NUM_THREADS=16; mpirun -np 1 ./solarhpc --benchmark-only

sweep-mpi: solarhpc
	@echo "Running MPI scaling sweep (1,2,4,8 ranks)..."
	@for np in 1 2 4 8; do \
	  echo "  mpirun -np $$np ..."; \
	  mpirun -np $$np ./solarhpc --benchmark-only 2>&1 | tail -3; \
	done

plots:
	@python3 plots/plot_orbits.py
	@python3 plots/plot_benchmarks.py
	@python3 plots/plot_validation.py
	@python3 plots/plot_eclipse_timeline.py
	@python3 plots/plot_mission.py
	@python3 plots/generate_tables.py
	@echo "  [OK] All plots saved to plots/output/"

dirs:
	@mkdir -p results/orbits results/eclipses results/missions \
	          results/benchmarks results/validation \
	          plots/output docs/paper/tables data logs tests

clean:
	rm -f src/fortran/*.o src/cpp/*.o *.o *.mod
	rm -f solarhpc tests/test_two_body tests/test_energy \
	       tests/test_eclipse tests/test_hohmann
	rm -f build.log

distclean: clean
	rm -rf results/ plots/output/ data/ logs/

help:
	@echo ""
	@echo "  SolarHPC Makefile targets:"
	@echo "  ─────────────────────────────────────────────────"
	@echo "  make all           Build solarhpc binary"
	@echo "  make tests         Build all test executables"
	@echo "  make run-tests     Build + run validation tests"
	@echo "  make fetch         Fetch initial conditions from JPL"
	@echo "  make fetch-ref     Fetch reference ephemeris"
	@echo "  make validate      Run --dry-run sanity check"
	@echo "  make run NP=4      Run simulation (4 MPI ranks)"
	@echo "  make benchmark     Run OpenMP thread sweep"
	@echo "  make sweep-mpi     Run MPI scaling (1,2,4,8 ranks)"
	@echo "  make plots         Generate all figures (run locally)"
	@echo "  make dirs          Create output directories"
	@echo "  make clean         Remove build artefacts"
	@echo "  make distclean     Remove everything including results"
	@echo "  ─────────────────────────────────────────────────"
	@echo "  Typical workflow:"
	@echo "    source setup_cluster.sh"
	@echo "    make dirs && make fetch"
	@echo "    make all"
	@echo "    make run-tests"
	@echo "    make run NP=4"
	@echo "    make benchmark"
	@echo "    make sweep-mpi"
	@echo "    [download results/] then: make plots"
	@echo ""
