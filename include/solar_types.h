#ifndef SOLAR_TYPES_H
#define SOLAR_TYPES_H
/* ============================================================
   solar_types.h  --  SolarHPC
   All shared C++ types, enums, and structs.
   Included by every C++ source file.
   ============================================================ */

#include <cstring>
#include <cmath>

/* ── Physical constants (mirror solar_constants.f90 exactly) ── */
static const double G_SI          = 6.674e-11;
static const double AU_TO_M       = 1.495978707e11;
static const double DAY_TO_SEC    = 86400.0;
static const double G_AU_DAY      = 2.95912e-4;   /* AU^3 kg^-1 day^-2 */
static const double AU_PER_DAY_TO_M_PER_S = AU_TO_M / DAY_TO_SEC;
static const double J2000_JD      = 2451545.0;
static const int    N_BODIES_CPP  = 10;

/* ── Body indices (0-based, C++ convention) ── */
enum BodyIdx {
    SUN=0, MERCURY=1, VENUS=2, EARTH=3, MOON_B=4,
    MARS=5, JUPITER=6, SATURN=7, URANUS=8, NEPTUNE=9
};

static const char* BODY_NAMES_CPP[10] = {
    "Sun","Mercury","Venus","Earth","Moon",
    "Mars","Jupiter","Saturn","Uranus","Neptune"
};

/* ── Target body and mission mode ── */
enum class TargetBody  { MOON=5, MARS=6 };
enum class MissionMode { HOHMANN=1, GRAVITY_ASSIST=2 };

/* ── Per-body state vector ── */
struct BodyState {
    char   name[32];
    double mass_kg;
    double pos[3];       /* AU, barycentric */
    double vel[3];       /* AU/day, barycentric */
    int    body_index;   /* 0-based */
};

/* ── Launch site ── */
struct LaunchSite {
    char   name[64];
    double latitude_deg;
    double longitude_deg;
    double altitude_m;
};

/* ── Mission configuration ── */
struct MissionConfig {
    LaunchSite  site;
    TargetBody  target;
    MissionMode mode;
    double spacecraft_mass_kg;
    double payload_mass_kg;
    double search_window_start_jd;
    double search_window_end_jd;
    int    timestep_days;
    int    n_simulation_years;
    int    output_decimation;
    bool   enable_perturbations;
    bool   benchmark_only;
    bool   dry_run;
    bool   no_perturb;
};

/* ── Mission result ── */
struct MissionResult {
    double launch_jd;
    double dv1_SI;
    double dv2_SI;
    double dv_total_SI;
    double fuel_mass_kg;
    double transfer_time_days;
    double landing_lat_deg;
    double landing_lon_deg;
    double landing_time_jd;
    double closest_approach_AU;
    bool   success;
    int    mpi_rank;
};

/* ── Default MissionConfig initialiser ── */
inline MissionConfig default_config() {
    MissionConfig c{};
    strncpy(c.site.name, "ISRO_SRIHARIKOTA", 63);
    c.site.latitude_deg   = 13.7199;
    c.site.longitude_deg  = 80.2304;
    c.site.altitude_m     = 0.0;
    c.target              = TargetBody::MOON;
    c.mode                = MissionMode::HOHMANN;
    c.spacecraft_mass_kg  = 5000.0;
    c.payload_mass_kg     = 1000.0;
    c.search_window_start_jd = J2000_JD;
    c.search_window_end_jd   = J2000_JD + 730.0;
    c.timestep_days       = 1;
    c.n_simulation_years  = 2;
    c.output_decimation   = 10;
    c.enable_perturbations= false;
    c.benchmark_only      = false;
    c.dry_run             = false;
    c.no_perturb          = false;
    return c;
}

#endif /* SOLAR_TYPES_H */
