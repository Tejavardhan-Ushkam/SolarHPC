/* ============================================================
   main.cpp  --  SolarHPC
   Top-level orchestrator. Runs the complete simulation pipeline.
   Stages: init -> load -> validate -> integrate -> eclipse
        -> mission -> benchmark -> finalise
   ============================================================ */
#include "solar_types.h"
#include "fortran_iface.h"
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <unistd.h>
#include <iomanip>

/* Forward declarations from other .cpp files */
bool         load_initial_conditions(const std::string&, std::vector<BodyState>&);
void         flatten_to_fortran(const std::vector<BodyState>&,
                                 double*, double*, double*);
void         print_body_table(const std::vector<BodyState>&);
double       bodies_initial_energy(const std::vector<BodyState>&);
MissionConfig parse_config(int, char**);
void         print_config(const MissionConfig&);
bool         init_output_dirs();
void         write_orbit_header(std::ofstream&);
void         write_orbit_step(std::ofstream&, int, double,
                               const double*, const double*,
                               int, const char[][32]);
void         write_benchmark_header();
void         write_mission_report(const MissionResult&, const MissionConfig&);
void         write_run_metadata(const MissionConfig&, int, int);
void         broadcast_initial_conditions(double*, double*, double*, int);
MissionResult mpi_launch_window_search(const MissionConfig&,
                                        double*, double*, double*,
                                        int, int);

/* ── Body name array for write_orbit_step ─────────────────── */
static char BNAMES[11][32] = {
    "Sun","Mercury","Venus","Earth","Moon",
    "Mars","Jupiter","Saturn","Uranus","Neptune",""
};

/* ── main ─────────────────────────────────────────────────── */
int main(int argc, char** argv)
{
    /* ── Stage 0: MPI init ────────────────────────────────── */
    MPI_Init(&argc, &argv);
    int rank, n_ranks;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &n_ranks);

    const char* omp_env = getenv("OMP_NUM_THREADS");
    int max_threads     = omp_env ? atoi(omp_env) : 1;

    if (rank == 0) {
        std::cout << "\n"
          << "  ================================================================\n"
          << "  SolarHPC v1.0 -- Parallel N-body Solar System Simulator\n"
          << "  MPI ranks: " << n_ranks
          << "  |  OMP_NUM_THREADS: " << max_threads << "\n"
          << "  Compiled: gfortran/g++ -O3 -fopenmp\n"
          << "  ================================================================\n\n";
    }

    /* ── Stage 1: Parse config ────────────────────────────── */
    MissionConfig cfg;
    try {
        cfg = parse_config(argc, argv);
    } catch (const std::exception& e) {
        if (rank == 0) std::cerr << "  [ERROR] " << e.what() << "\n";
        MPI_Finalize(); return 1;
    }
    if (rank == 0) print_config(cfg);

    /* ── Create output directories ────────────────────────── */
    if (rank == 0) init_output_dirs();
    MPI_Barrier(MPI_COMM_WORLD);

    /* ── Stage 2: Load initial conditions (all ranks) ──────── */
    int N = N_BODIES_CPP;
    std::vector<double> pos(3*N), vel(3*N), mass(N);

    if (rank == 0) {
        std::vector<BodyState> bodies;
        if (!load_initial_conditions("data/initial_conditions.dat", bodies)) {
            std::cerr << "  [ERROR] Failed to load initial conditions.\n"
                      << "  Run: python3 fetch_jpl.py\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        flatten_to_fortran(bodies, pos.data(), vel.data(), mass.data());
        print_body_table(bodies);
        double E0 = bodies_initial_energy(bodies);
        std::cout << std::scientific << std::setprecision(6)
                  << "  Initial energy E0 = " << E0 << " J\n\n";
    }
    broadcast_initial_conditions(pos.data(), vel.data(), mass.data(), rank);

    /* Save copies of initial state for benchmarking */
    std::vector<double> pos0(pos), vel0(vel);

    /* ── Stage 3: Validation tests (rank 0) ──────────────── */
    if (!cfg.benchmark_only && !cfg.dry_run && rank == 0) {
        int all_passed=0, n_passed=0, n_total=0;
        validation_suite__run_all_tests_(pos0.data(), vel0.data(), mass.data(),
                                          &all_passed, &n_passed, &n_total);
        std::cout << "  [VALIDATION] " << n_passed << "/" << n_total
                  << " tests passed\n\n";
        write_run_metadata(cfg, n_ranks, max_threads);
    }

    /* ── dry-run exits here ────────────────────────────────── */
    if (cfg.dry_run) {
        if (rank == 0)
            std::cout << "  [DRY-RUN] Complete. All systems nominal.\n\n";
        MPI_Finalize(); return 0;
    }

    /* ── benchmark-only: skip straight to sweep ────────────── */
    if (!cfg.benchmark_only) {

        /* ── Stage 4: Orbital integration ──────────────────── */
        if (rank == 0) {
            std::cout << "  [STAGE 4] Orbital integration ("
                      << cfg.n_simulation_years << " years, dt="
                      << cfg.timestep_days << " day)...\n";
        }

        int n_steps  = cfg.n_simulation_years * 365 / cfg.timestep_days;
        int n_stored = n_steps / cfg.output_decimation + 1;
        double dt    = (double)cfg.timestep_days;

        /* Allocate position history for eclipse scan (rank 0 only) */
        std::vector<double> pos_history;   /* 3 * N * n_stored */
        std::ofstream orbit_file;

        if (rank == 0) {
            pos_history.resize(3 * N * n_stored, 0.0);
            orbit_file.open("results/orbits/orbit_trajectories.csv");
            write_orbit_header(orbit_file);
        }

        double E0_saved = 0.0, E_final = 0.0;
        if (rank == 0)
            integrators__compute_energy_(&N, pos.data(), vel.data(),
                                          mass.data(), &E0_saved);

        /* Main integration loop (rank 0 drives; all ranks
           could share but orbit output is serial) */
        if (rank == 0) {
            int stored_idx = 0;
            for (int s = 0; s < n_steps; ++s) {
                integrators__leapfrog_step_(&N, pos.data(), vel.data(),
                                             mass.data(), &dt);

                if (s % cfg.output_decimation == 0) {
                    double jd = J2000_JD + (double)s * dt;
                    write_orbit_step(orbit_file, s, jd,
                                     pos.data(), vel.data(), N, BNAMES);
                    /* Store for eclipse scan */
                    if (stored_idx < n_stored) {
                        for (int b = 0; b < N; ++b) {
                            pos_history[0 + b*3 + stored_idx*3*N] = pos[0+b*3];
                            pos_history[1 + b*3 + stored_idx*3*N] = pos[1+b*3];
                            pos_history[2 + b*3 + stored_idx*3*N] = pos[2+b*3];
                        }
                        ++stored_idx;
                    }
                }
            }
            orbit_file.close();
            integrators__compute_energy_(&N, pos.data(), vel.data(),
                                          mass.data(), &E_final);
            double drift = std::abs(E_final - E0_saved) / std::abs(E0_saved) * 100.0;
            std::cout << std::scientific << std::setprecision(3)
                      << "  Energy drift after " << cfg.n_simulation_years
                      << " yr: " << drift << " %\n\n";
        }

        /* Broadcast final positions to all ranks for mission planning */
        MPI_Bcast(pos.data(), 3*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(vel.data(), 3*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        /* ── Stage 5: Eclipse prediction (rank 0) ─────────── */
        if (rank == 0) {
            std::cout << "  [STAGE 5] Scanning for eclipses...\n";
            int n_solar = 0, n_lunar = 0;
            double dt_hist = (double)(cfg.timestep_days * cfg.output_decimation);
            double t_start = J2000_JD;
            eclipse_geom__scan_eclipses_omp_(&n_stored, pos_history.data(),
                                              &dt_hist, &t_start,
                                              &n_solar, &n_lunar);
            std::cout << "  Found " << n_solar << " solar eclipses, "
                      << n_lunar << " lunar eclipses\n\n";

            /* Verify 2024-04-08 eclipse if simulation covers it */
            if (cfg.n_simulation_years >= 24) {
                int passed = 0;
                const char* ef = "results/eclipses/eclipse_predictions.csv";
                int eflen = (int)strlen(ef);
                eclipse_geom__verify_eclipse_2024_(
                    const_cast<char*>(ef), &eflen, &passed);
            }
        }

        /* ── Stage 6: Mission planning (MPI distributed) ─── */
        if (rank == 0)
            std::cout << "  [STAGE 6] Launch window search ("
                      << n_ranks << " MPI ranks)...\n";

        MissionResult best = mpi_launch_window_search(cfg,
                                                       pos0.data(), vel0.data(),
                                                       mass.data(),
                                                       rank, n_ranks);

        if (rank == 0) write_mission_report(best, cfg);

    } /* end !benchmark_only */

    /* ── Stage 7: Thread sweep benchmarks (rank 0) ─────────── */
    if (rank == 0) {
        std::cout << "  [STAGE 7] OpenMP thread sweep benchmarks...\n\n";
        write_benchmark_header();

        int n_sweep = N;
        nbody_force__run_force_sweep_(&n_sweep, pos0.data(), mass.data());

        int n_steps_bench = 50;
        double dt_bench   = 1.0;
        integrators__run_integrator_sweep_(&n_sweep, pos0.data(), vel0.data(),
                                            mass.data(), &dt_bench, &n_steps_bench);

        /* Eclipse sweep needs a small position history */
        int n_eph = 200;
        std::vector<double> ph(3 * N * n_eph, 0.0);
        /* Fill with the current trajectory positions (first n_eph steps) */
        {
            std::vector<double> ptmp(pos0), vtmp(vel0);
            double dt1 = 1.0;
            for (int s = 0; s < n_eph; ++s) {
                integrators__leapfrog_step_(&N, ptmp.data(), vtmp.data(),
                                             mass.data(), &dt1);
                for (int b = 0; b < N; ++b) {
                    ph[0+b*3+s*3*N] = ptmp[0+b*3];
                    ph[1+b*3+s*3*N] = ptmp[1+b*3];
                    ph[2+b*3+s*3*N] = ptmp[2+b*3];
                }
            }
        }
        double dt_eph = 1.0, t0_eph = J2000_JD;
        eclipse_geom__run_eclipse_sweep_(&n_eph, ph.data(), &dt_eph, &t0_eph);

        trajectory_prop__run_trajectory_sweep_(&n_sweep, pos0.data(),
                                                vel0.data(), mass.data());
    }

    /* ── Stage 8: Finalise ──────────────────────────────────── */
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        std::cout << "  ================================================================\n"
                  << "  Simulation complete. All results in results/\n"
                  << "  Download results/ then run locally: make plots\n"
                  << "  ================================================================\n\n";
    }

    MPI_Finalize();
    return 0;
}

