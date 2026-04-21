#ifndef FORTRAN_IFACE_H
#define FORTRAN_IFACE_H
/* ============================================================
   fortran_iface.h  --  SolarHPC
   gfortran module procedure name mangling:
     __modulename_MOD_subname
   (double underscore prefix, _MOD_ separator, no trailing _)
   Verified against nm output of compiled .o files.
   ============================================================ */

#ifdef __cplusplus
extern "C" {
#endif

/* ── nbody_force ─────────────────────────────────────────── */
void __nbody_force_MOD_compute_forces(int* n, double* pos,
                                       double* mass, double* acc);
void __nbody_force_MOD_verify_forces(int* passed);
void __nbody_force_MOD_run_force_sweep(int* n, double* pos, double* mass);

/* ── integrators ─────────────────────────────────────────── */
void __integrators_MOD_leapfrog_step(int* n, double* pos, double* vel,
                                      double* mass, double* dt);
void __integrators_MOD_rk4_step(int* n, double* pos, double* vel,
                                  double* mass, double* dt);
void __integrators_MOD_compute_energy(int* n, double* pos, double* vel,
                                       double* mass, double* E_total);
void __integrators_MOD_run_integrator_sweep(int* n, double* pos, double* vel,
                                             double* mass, double* dt,
                                             int* n_steps);

/* ── eclipse_geom ────────────────────────────────────────── */
void __eclipse_geom_MOD_scan_eclipses_omp(int* n_steps, double* pos_history,
                                           double* dt_days, double* t_start_jd,
                                           int* n_solar, int* n_lunar);
void __eclipse_geom_MOD_verify_eclipse_2024(char* eclipse_file, int* file_len,
                                              int* passed);
void __eclipse_geom_MOD_run_eclipse_sweep(int* n_steps, double* pos_history,
                                           double* dt_days, double* t_start_jd);

/* ── trajectory_prop ─────────────────────────────────────── */
void __trajectory_prop_MOD_hohmann_transfer(double* r1, double* r2,
                                             double* dv1_SI, double* dv2_SI,
                                             double* transfer_days);
void __trajectory_prop_MOD_tsiolkovsky(double* dv_SI, double* isp,
                                        double* m_payload, double* m_fuel);
void __trajectory_prop_MOD_find_optimal_launch_date(
    double* pos0, double* vel0, double* mass, int* target_idx,
    double* jd_start, double* jd_end, double* dt_sim,
    double* best_jd, double* best_dv_SI);
void __trajectory_prop_MOD_run_trajectory_sweep(int* n_pl, double* pos,
                                                 double* vel, double* mass);

/* ── validation_suite ────────────────────────────────────── */
void __validation_suite_MOD_run_all_tests(double* pos0, double* vel0,
    double* mass, int* all_passed, int* n_passed, int* n_total);

#ifdef __cplusplus
}
#endif

/* ── Convenience macros matching old names in main.cpp ────── */
#define nbody_force__compute_forces_        __nbody_force_MOD_compute_forces
#define nbody_force__verify_forces_         __nbody_force_MOD_verify_forces
#define nbody_force__run_force_sweep_       __nbody_force_MOD_run_force_sweep
#define integrators__leapfrog_step_         __integrators_MOD_leapfrog_step
#define integrators__rk4_step_              __integrators_MOD_rk4_step
#define integrators__compute_energy_        __integrators_MOD_compute_energy
#define integrators__run_integrator_sweep_  __integrators_MOD_run_integrator_sweep
#define eclipse_geom__scan_eclipses_omp_    __eclipse_geom_MOD_scan_eclipses_omp
#define eclipse_geom__verify_eclipse_2024_  __eclipse_geom_MOD_verify_eclipse_2024
#define eclipse_geom__run_eclipse_sweep_    __eclipse_geom_MOD_run_eclipse_sweep
#define trajectory_prop__hohmann_transfer_  __trajectory_prop_MOD_hohmann_transfer
#define trajectory_prop__tsiolkovsky_       __trajectory_prop_MOD_tsiolkovsky
#define trajectory_prop__find_optimal_launch_date_ \
        __trajectory_prop_MOD_find_optimal_launch_date
#define trajectory_prop__run_trajectory_sweep_ \
        __trajectory_prop_MOD_run_trajectory_sweep
#define validation_suite__run_all_tests_    __validation_suite_MOD_run_all_tests

#endif /* FORTRAN_IFACE_H */
