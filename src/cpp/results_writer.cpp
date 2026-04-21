/* ============================================================
   results_writer.cpp  --  SolarHPC
   Centralised output writer. Every CSV written by C++ goes
   through here ensuring consistent headers and formats.
   ============================================================ */
#include "solar_types.h"
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <cstdio>
#include <cstring>
#include <unistd.h>
#include <sys/stat.h>

/* ── init_output_dirs ─────────────────────────────────────── */
bool init_output_dirs()
{
    const char* dirs[] = {
        "results", "results/orbits", "results/eclipses",
        "results/missions", "results/benchmarks", "results/validation",
        "plots/output", "logs", nullptr
    };
    for (int i = 0; dirs[i]; ++i) {
        std::string cmd = std::string("mkdir -p ") + dirs[i];
        if (system(cmd.c_str()) != 0) {
            std::cerr << "  [ERROR] Could not create: " << dirs[i] << "\n";
            return false;
        }
    }
    return true;
}

/* ── write_orbit_header ───────────────────────────────────── */
void write_orbit_header(std::ofstream& f)
{
    f << "step,julian_day,body_index,body_name,"
      << "x_au,y_au,z_au,vx_au_day,vy_au_day,vz_au_day\n";
}

/* ── write_orbit_step ─────────────────────────────────────── */
/* pos and vel are flat Fortran column-major arrays: pos[k + i*3]
   where k=0,1,2 and i=0..n_bodies-1                            */
void write_orbit_step(std::ofstream& f, int step, double jd,
                       const double* pos, const double* vel,
                       int n_bodies, const char body_names[][32])
{
    for (int i = 0; i < n_bodies; ++i) {
        f << std::fixed << std::setprecision(6)
          << step << ","
          << std::setprecision(4) << jd << ","
          << i+1 << ","                           /* 1-based for paper */
          << body_names[i] << ","
          << std::setprecision(10)
          << pos[0+i*3] << "," << pos[1+i*3] << "," << pos[2+i*3] << ","
          << vel[0+i*3] << "," << vel[1+i*3] << "," << vel[2+i*3] << "\n";
    }
}

/* ── write_spacecraft_path_header ─────────────────────────── */
void write_spacecraft_path_header(std::ofstream& f)
{
    f << "step,julian_day,x_au,y_au,z_au\n";
}

/* ── write_mission_report ─────────────────────────────────── */
void write_mission_report(const MissionResult& r, const MissionConfig& cfg)
{
    /* ── CSV row ─── */
    const char* path_csv = "results/missions/mission_report.csv";
    bool exists = false;
    {
        std::ifstream chk(path_csv);
        exists = chk.good();
    }
    std::ofstream f(path_csv, std::ios::app);
    if (!f) { std::cerr << "  [ERROR] Cannot write mission_report.csv\n"; return; }
    if (!exists) {
        f << "launch_jd,target_body,mission_mode,"
          << "dv1_SI,dv2_SI,dv_total_SI,fuel_mass_kg,"
          << "transfer_time_days,landing_lat_deg,landing_lon_deg,"
          << "landing_time_jd,closest_approach_AU\n";
    }
    const char* tgt  = (cfg.target == TargetBody::MOON) ? "Moon" : "Mars";
    const char* mode = (cfg.mode   == MissionMode::HOHMANN) ? "Hohmann" : "GravityAssist";

    f << std::fixed << std::setprecision(4)
      << r.launch_jd     << "," << tgt << "," << mode     << ","
      << std::setprecision(2)
      << r.dv1_SI        << "," << r.dv2_SI << "," << r.dv_total_SI << ","
      << r.fuel_mass_kg  << ","
      << std::setprecision(3)
      << r.transfer_time_days << ","
      << r.landing_lat_deg    << "," << r.landing_lon_deg << ","
      << std::setprecision(4)
      << r.landing_time_jd    << ","
      << std::setprecision(8)
      << r.closest_approach_AU << "\n";

    /* ── Human-readable summary ─── */
    std::ofstream s("results/missions/mission_summary.txt", std::ios::app);
    if (s) {
        s << "========================================\n"
          << "  Mission to " << tgt << " (" << mode << ")\n"
          << "  Launch JD   : " << std::fixed << std::setprecision(1) << r.launch_jd << "\n"
          << "  Total dv    : " << std::setprecision(1) << r.dv_total_SI << " m/s ("
          << r.dv_total_SI/1000.0 << " km/s)\n"
          << "  Fuel needed : " << r.fuel_mass_kg << " kg\n"
          << "  Transfer    : " << r.transfer_time_days << " days\n";
        if (r.success) {
            s << "  Landing lat : " << std::setprecision(4) << r.landing_lat_deg << " deg\n"
              << "  Landing lon : " << r.landing_lon_deg << " deg\n";
        } else {
            s << "  Status      : Did not reach target within simulation window\n";
        }
        s << "========================================\n\n";
    }
    std::cout << "  [OK] Mission report written to results/missions/\n";
}

/* ── write_benchmark_header ───────────────────────────────── */
void write_benchmark_header()
{
    const char* path = "results/benchmarks/thread_sweep_results.csv";
    std::ifstream chk(path);
    if (chk.good()) return;   /* already exists from Fortran writes */
    std::ofstream f(path);
    if (f) {
        f << "kernel_name,n_threads,run1_sec,run2_sec,run3_sec,run4_sec,run5_sec,"
          << "mean_sec,stddev_sec,speedup,efficiency_pct\n";
        std::cout << "  [OK] Created results/benchmarks/thread_sweep_results.csv\n";
    }
}

/* ── write_mpi_scaling_row ────────────────────────────────── */
void write_mpi_scaling_row(int n_ranks, double wall_time,
                             double speedup, double efficiency)
{
    const char* path = "results/benchmarks/mpi_scaling.csv";
    bool exists = false;
    { std::ifstream chk(path); exists = chk.good(); }
    std::ofstream f(path, std::ios::app);
    if (!f) return;
    if (!exists)
        f << "n_ranks,wall_time_sec,speedup,efficiency_pct\n";
    f << std::fixed << std::setprecision(6)
      << n_ranks << "," << wall_time << ","
      << std::setprecision(4) << speedup << ","
      << std::setprecision(2) << efficiency << "\n";
}

/* ── get_git_hash ─────────────────────────────────────────── */
static std::string get_git_hash()
{
    FILE* pipe = popen("git rev-parse --short HEAD 2>/dev/null", "r");
    if (!pipe) return "unknown";
    char buf[16] = {0};
    if (fgets(buf, sizeof(buf), pipe)) {
        pclose(pipe);
        std::string s(buf);
        if (!s.empty() && s.back() == '\n') s.pop_back();
        return s.empty() ? "unknown" : s;
    }
    pclose(pipe);
    return "unknown";
}

/* ── write_run_metadata ───────────────────────────────────── */
void write_run_metadata(const MissionConfig& cfg,
                         int n_mpi_ranks, int max_omp_threads)
{
    char hostname[256] = {0};
    if (gethostname(hostname, sizeof(hostname)) != 0)
        strncpy(hostname, "unknown", 255);

    std::time_t now = std::time(nullptr);
    char ts[64];
    std::strftime(ts, sizeof(ts), "%Y-%m-%d %H:%M:%S UTC", std::gmtime(&now));

    std::ofstream f("results/run_metadata.txt");
    if (!f) return;

    const char* tgt  = (cfg.target == TargetBody::MOON) ? "Moon" : "Mars";
    const char* mode = (cfg.mode   == MissionMode::HOHMANN) ? "Hohmann" : "GravityAssist";

    f << "=== SolarHPC Run Metadata ===\n"
      << "Timestamp        : " << ts                   << "\n"
      << "Hostname         : " << hostname              << "\n"
      << "Git commit       : " << get_git_hash()        << "\n"
      << "MPI ranks        : " << n_mpi_ranks           << "\n"
      << "OMP threads (max): " << max_omp_threads       << "\n"
      << "Compiler flags   : -O3 -fopenmp -fdefault-real-8\n"
      << "Target body      : " << tgt                   << "\n"
      << "Mission mode     : " << mode                  << "\n"
      << "Launch site      : " << cfg.site.name         << "\n"
      << "Spacecraft mass  : " << cfg.spacecraft_mass_kg << " kg\n"
      << "Sim years        : " << cfg.n_simulation_years << "\n"
      << "Timestep         : " << cfg.timestep_days     << " day\n"
      << "Perturbations    : " << (cfg.enable_perturbations ? "ON" : "OFF") << "\n";

    std::cout << "  [OK] Run metadata: results/run_metadata.txt\n";
}
