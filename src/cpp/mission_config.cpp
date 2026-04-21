/* ============================================================
   mission_config.cpp  --  SolarHPC
   Parses command-line flags into MissionConfig.
   Provides the 5 built-in launch sites.
   ============================================================ */
#include "solar_types.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <stdexcept>
#include <cstring>

/* ── Built-in launch sites ───────────────────────────────── */
static const LaunchSite LAUNCH_SITES[] = {
    {"ISRO_SRIHARIKOTA",    13.7199,   80.2304,   0.0},
    {"NASA_KENNEDY",         28.5729,  -80.6490,   3.0},
    {"ROSCOSMOS_BAIKONUR",   45.9200,   63.3420, 100.0},
    {"ESA_KOUROU",            5.2320,  -52.7683,  14.0},
    {"JAXA_TANEGASHIMA",     30.4005,  130.9749,  30.0},
};
static const int N_SITES = 5;

/* ── get_launch_site ──────────────────────────────────────── */
LaunchSite get_launch_site(const std::string& name)
{
    for (int i = 0; i < N_SITES; ++i)
        if (name == LAUNCH_SITES[i].name) return LAUNCH_SITES[i];

    std::string msg = "Unknown launch site: " + name + "\nAvailable: ";
    for (int i = 0; i < N_SITES; ++i) {
        msg += LAUNCH_SITES[i].name;
        if (i < N_SITES-1) msg += ", ";
    }
    throw std::invalid_argument(msg);
}

/* ── parse_config ─────────────────────────────────────────── */
MissionConfig parse_config(int argc, char** argv)
{
    MissionConfig cfg = default_config();

    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);

        auto val_of = [&](const std::string& prefix) -> std::string {
            if (arg.rfind(prefix, 0) == 0) return arg.substr(prefix.size());
            return "";
        };

        std::string v;
        if (!(v = val_of("--site=")).empty()) {
            cfg.site = get_launch_site(v);
        } else if (!(v = val_of("--target=")).empty()) {
            if      (v == "MOON") cfg.target = TargetBody::MOON;
            else if (v == "MARS") cfg.target = TargetBody::MARS;
            else throw std::invalid_argument("--target must be MOON or MARS");
        } else if (!(v = val_of("--mode=")).empty()) {
            if      (v == "HOHMANN")        cfg.mode = MissionMode::HOHMANN;
            else if (v == "GRAVITY_ASSIST") cfg.mode = MissionMode::GRAVITY_ASSIST;
            else throw std::invalid_argument("--mode must be HOHMANN or GRAVITY_ASSIST");
        } else if (!(v = val_of("--mass=")).empty()) {
            cfg.spacecraft_mass_kg = std::stod(v);
        } else if (!(v = val_of("--payload=")).empty()) {
            cfg.payload_mass_kg = std::stod(v);
        } else if (!(v = val_of("--years=")).empty()) {
            cfg.n_simulation_years = std::stoi(v);
        } else if (!(v = val_of("--dt=")).empty()) {
            cfg.timestep_days = std::stoi(v);
        } else if (!(v = val_of("--decimate=")).empty()) {
            cfg.output_decimation = std::stoi(v);
        } else if (!(v = val_of("--search-start=")).empty()) {
            cfg.search_window_start_jd = std::stod(v);
        } else if (!(v = val_of("--search-end=")).empty()) {
            cfg.search_window_end_jd = std::stod(v);
        } else if (arg == "--perturb") {
            cfg.enable_perturbations = true;
        } else if (arg == "--no-perturb") {
            cfg.no_perturb = true;
        } else if (arg == "--benchmark-only") {
            cfg.benchmark_only = true;
        } else if (arg == "--dry-run") {
            cfg.dry_run = true;
        } else if (arg == "--help" || arg == "-h") {
            std::cout << "\nSolarHPC -- usage:\n"
                      << "  ./solarhpc [options]\n\n"
                      << "  --site=<NAME>          Launch site (default: ISRO_SRIHARIKOTA)\n"
                      << "  --target=<MOON|MARS>   Target body (default: MOON)\n"
                      << "  --mode=<HOHMANN|GRAVITY_ASSIST>\n"
                      << "  --mass=<kg>            Spacecraft mass (default: 5000)\n"
                      << "  --payload=<kg>         Payload mass (default: 1000)\n"
                      << "  --years=<N>            Simulation years (default: 2)\n"
                      << "  --dt=<days>            Timestep (default: 1)\n"
                      << "  --perturb              Enable J2/SRP/GR perturbations\n"
                      << "  --benchmark-only       Skip to thread sweep\n"
                      << "  --dry-run              Validate only, no simulation\n\n"
                      << "  Available sites: ISRO_SRIHARIKOTA, NASA_KENNEDY,\n"
                      << "                   ROSCOSMOS_BAIKONUR, ESA_KOUROU, JAXA_TANEGASHIMA\n\n";
            exit(0);
        } else {
            std::cerr << "  [WARN] Unknown flag: " << arg << " (ignored)\n";
        }
    }
    return cfg;
}

/* ── print_config ─────────────────────────────────────────── */
void print_config(const MissionConfig& cfg)
{
    const char* target_str = (cfg.target == TargetBody::MOON) ? "Moon" : "Mars";
    const char* mode_str   = (cfg.mode   == MissionMode::HOHMANN) ? "Hohmann transfer"
                                                                    : "Gravity assist";
    int n_steps = cfg.n_simulation_years * 365 / cfg.timestep_days;

    std::cout << "\n  ============================================\n"
              << "   SolarHPC Mission Configuration\n"
              << "  ============================================\n"
              << std::fixed << std::setprecision(4)
              << "   Launch site  : " << cfg.site.name
              << " (" << cfg.site.latitude_deg << "N, "
              << cfg.site.longitude_deg << "E)\n"
              << "   Target body  : " << target_str << "\n"
              << "   Mission mode : " << mode_str << "\n"
              << std::setprecision(1)
              << "   Spacecraft   : " << cfg.spacecraft_mass_kg << " kg"
              << "  (payload: " << cfg.payload_mass_kg << " kg)\n"
              << "   Sim duration : " << cfg.n_simulation_years << " years"
              << "  (" << n_steps << " steps, dt=" << cfg.timestep_days << " day)\n"
              << "   Perturbations: " << (cfg.enable_perturbations ? "enabled" : "disabled") << "\n"
              << "  ============================================\n\n";
}
