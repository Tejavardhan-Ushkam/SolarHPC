! ============================================================
! MODULE: trajectory_prop
! PROJECT: SolarHPC
! PURPOSE: Spacecraft trajectory computation.
!          Hohmann transfer, gravity assist, full N-body
!          propagation, Tsiolkovsky fuel budget, landing coords.
!          Compiled SEVENTH. Depends on: solar_constants,
!          omp_timer, nbody_force, integrators, coordinate_transform.
! ============================================================
MODULE trajectory_prop
  USE solar_constants
  USE omp_timer
  USE nbody_force
  USE integrators
  USE coordinate_transform
  USE omp_lib
  IMPLICIT NONE

  INTEGER, PARAMETER :: N_EXT = N_BODIES + 1   ! planets + spacecraft

CONTAINS

  ! ----------------------------------------------------------
  ! hohmann_transfer: two-impulse Hohmann from r1 to r2 (AU).
  ! Reference: Bate, Mueller, White (1971) Chapter 4.
  ! ----------------------------------------------------------
  SUBROUTINE hohmann_transfer(r1_AU, r2_AU, dv1_SI, dv2_SI, transfer_days)
    REAL(dp), INTENT(IN)  :: r1_AU, r2_AU
    REAL(dp), INTENT(OUT) :: dv1_SI, dv2_SI, transfer_days
    REAL(dp) :: mu_AU, v1, v2, v_peri, v_apo, a_transfer

    mu_AU         = G_AU_DAY * MASS_KG(SUN_IDX)
    v1            = SQRT(mu_AU / r1_AU)
    v2            = SQRT(mu_AU / r2_AU)
    a_transfer    = 0.5d0 * (r1_AU + r2_AU)
    v_peri        = SQRT(mu_AU * (2.0d0/r1_AU    - 1.0d0/a_transfer))
    v_apo         = SQRT(mu_AU * (2.0d0/r2_AU    - 1.0d0/a_transfer))
    dv1_SI        = ABS(v_peri - v1)   * AU_PER_DAY_TO_M_PER_S
    dv2_SI        = ABS(v2    - v_apo) * AU_PER_DAY_TO_M_PER_S
    transfer_days = PI * SQRT(a_transfer**3 / mu_AU)
  END SUBROUTINE hohmann_transfer

  ! ----------------------------------------------------------
  ! tsiolkovsky: rocket equation. Returns fuel mass (kg).
  ! ----------------------------------------------------------
  SUBROUTINE tsiolkovsky(dv_total_SI, isp, m_payload_kg, m_fuel_kg)
    REAL(dp), INTENT(IN)  :: dv_total_SI, isp, m_payload_kg
    REAL(dp), INTENT(OUT) :: m_fuel_kg
    m_fuel_kg = m_payload_kg * (EXP(dv_total_SI / (isp * G0)) - 1.0d0)
  END SUBROUTINE tsiolkovsky

  ! ----------------------------------------------------------
  ! compute_launch_vector: surface site -> barycentric ICs.
  ! ----------------------------------------------------------
  SUBROUTINE compute_launch_vector(lat_deg, lon_deg, jd, target_body, &
                                    earth_pos, earth_vel,              &
                                    sc_pos, sc_vel, dv_launch_SI)
    REAL(dp), INTENT(IN)  :: lat_deg, lon_deg, jd
    INTEGER,  INTENT(IN)  :: target_body
    REAL(dp), INTENT(IN)  :: earth_pos(3), earth_vel(3)
    REAL(dp), INTENT(OUT) :: sc_pos(3), sc_vel(3), dv_launch_SI

    REAL(dp) :: x_ecef, y_ecef, z_ecef
    REAL(dp) :: x_eci_m, y_eci_m, z_eci_m
    REAL(dp) :: v_rot_x, v_rot_y, v_rot_z
    REAL(dp) :: r_target, dv1_SI, dv2_SI, tf
    REAL(dp) :: prograde_x, prograde_y, prograde_z, v_earth_mag
    REAL(dp) :: dv1_AU_per_day
    REAL(dp) :: mu_earth, v_leo, v_trans, r_leo, r_moon_orbit

    CALL geographic_to_ecef(lat_deg, lon_deg, 0.0d0, x_ecef, y_ecef, z_ecef)
    CALL ecef_to_eci(x_ecef, y_ecef, z_ecef, jd, x_eci_m, y_eci_m, z_eci_m)
    CALL eci_to_barycentric(x_eci_m, y_eci_m, z_eci_m, earth_pos, &
                             sc_pos(1), sc_pos(2), sc_pos(3))
    CALL earth_rotation_velocity(lat_deg, lon_deg, v_rot_x, v_rot_y, v_rot_z)

    IF (target_body == MOON_IDX) THEN
      mu_earth     = G_AU_DAY * MASS_KG(EARTH_IDX)
      r_leo        = EARTH_RADIUS_AU + 1.6d-5
      r_moon_orbit = 0.00257d0
      v_leo        = SQRT(mu_earth / r_leo)
      v_trans      = SQRT(mu_earth * (2.0d0/r_leo - 1.0d0/(r_leo + r_moon_orbit)))
      dv1_SI       = ABS(v_trans - v_leo) * AU_PER_DAY_TO_M_PER_S
    ELSE
      r_target = 1.524d0
      CALL hohmann_transfer(1.0d0, r_target, dv1_SI, dv2_SI, tf)
    END IF

    dv_launch_SI = dv1_SI
    v_earth_mag  = SQRT(DOT_PRODUCT(earth_vel, earth_vel))

    IF (v_earth_mag > 1.0d-20) THEN
      prograde_x = earth_vel(1) / v_earth_mag
      prograde_y = earth_vel(2) / v_earth_mag
      prograde_z = earth_vel(3) / v_earth_mag
    ELSE
      prograde_x = 1.0d0; prograde_y = 0.0d0; prograde_z = 0.0d0
    END IF

    dv1_AU_per_day = dv1_SI * M_PER_S_TO_AU_PER_DAY
    sc_vel(1) = earth_vel(1) + v_rot_x + dv1_AU_per_day * prograde_x
    sc_vel(2) = earth_vel(2) + v_rot_y + dv1_AU_per_day * prograde_y
    sc_vel(3) = earth_vel(3) + v_rot_z + dv1_AU_per_day * prograde_z
  END SUBROUTINE compute_launch_vector

  ! ----------------------------------------------------------
  ! propagate_spacecraft: full N-body trajectory integration.
  ! ----------------------------------------------------------
  SUBROUTINE propagate_spacecraft(n_pl, pos_planets, vel_planets, mass_planets, &
                                   sc_pos0, sc_vel0, sc_mass_kg,                 &
                                   n_steps, dt_days, target_body_idx,            &
                                   sc_path, n_path_pts,                           &
                                   landing_lat, landing_lon, landing_jd,          &
                                   closest_AU, mission_success, launch_jd)
    INTEGER,  INTENT(IN)  :: n_pl, n_steps, target_body_idx, n_path_pts
    REAL(dp), INTENT(IN)  :: pos_planets(3, n_pl), vel_planets(3, n_pl)
    REAL(dp), INTENT(IN)  :: mass_planets(n_pl), sc_pos0(3), sc_vel0(3)
    REAL(dp), INTENT(IN)  :: sc_mass_kg, dt_days, launch_jd
    REAL(dp), INTENT(OUT) :: sc_path(3, n_path_pts)
    REAL(dp), INTENT(OUT) :: landing_lat, landing_lon, landing_jd, closest_AU
    LOGICAL,  INTENT(OUT) :: mission_success

    INTEGER,  PARAMETER :: PATH_DEC = 10
    REAL(dp) :: pos_ext(3, N_EXT), vel_ext(3, N_EXT), mass_ext(N_EXT)
    REAL(dp) :: r_vec(3), dist, arrival_thresh, jd_cur
    INTEGER  :: s, path_idx

    pos_ext(:, 1:n_pl) = pos_planets
    vel_ext(:, 1:n_pl) = vel_planets
    mass_ext(1:n_pl)   = mass_planets
    pos_ext(:, N_EXT)  = sc_pos0
    vel_ext(:, N_EXT)  = sc_vel0
    mass_ext(N_EXT)    = 1.0d-30

    arrival_thresh  = BODY_RADIUS_AU(target_body_idx) * 50.0d0
    closest_AU      = 1.0d10
    mission_success = .FALSE.
    landing_lat     = 0.0d0
    landing_lon     = 0.0d0
    landing_jd      = launch_jd
    path_idx        = 0

    DO s = 1, n_steps
      CALL leapfrog_step_extended(N_EXT, pos_ext, vel_ext, mass_ext, dt_days)
      jd_cur = launch_jd + REAL(s, dp) * dt_days

      IF (MOD(s, PATH_DEC) == 0 .AND. path_idx < n_path_pts) THEN
        path_idx = path_idx + 1
        sc_path(:, path_idx) = pos_ext(:, N_EXT)
      END IF

      r_vec = pos_ext(:, N_EXT) - pos_ext(:, target_body_idx)
      dist  = SQRT(DOT_PRODUCT(r_vec, r_vec))
      IF (dist < closest_AU) closest_AU = dist

      IF (dist < arrival_thresh) THEN
        mission_success = .TRUE.
        landing_jd = jd_cur
        CALL barycentric_to_surface(pos_ext(:, N_EXT),            &
                                     pos_ext(:, target_body_idx),  &
                                     target_body_idx, jd_cur,      &
                                     landing_lat, landing_lon)
        EXIT
      END IF
    END DO

    DO WHILE (path_idx < n_path_pts)
      path_idx = path_idx + 1
      sc_path(:, path_idx) = pos_ext(:, N_EXT)
    END DO
  END SUBROUTINE propagate_spacecraft

  ! ----------------------------------------------------------
  ! gravity_assist_dv: patched conic flyby delta-v gain.
  ! ----------------------------------------------------------
  SUBROUTINE gravity_assist_dv(flyby_body_idx, sc_vel, planet_vel, dv_gain_SI)
    INTEGER,  INTENT(IN)  :: flyby_body_idx
    REAL(dp), INTENT(IN)  :: sc_vel(3), planet_vel(3)
    REAL(dp), INTENT(OUT) :: dv_gain_SI

    REAL(dp) :: v_inf(3), v_inf_mag, r_p, mu_SI, turn_sin, delta

    v_inf      = sc_vel - planet_vel
    v_inf_mag  = SQRT(DOT_PRODUCT(v_inf, v_inf)) * AU_PER_DAY_TO_M_PER_S
    mu_SI      = G_SI * MASS_KG(flyby_body_idx)
    r_p        = BODY_RADIUS_AU(flyby_body_idx) * AU_TO_M + 200000.0d0
    turn_sin   = 1.0d0 / (1.0d0 + r_p * v_inf_mag**2 / mu_SI)
    delta      = 2.0d0 * ASIN(MIN(turn_sin, 1.0d0))
    dv_gain_SI = 2.0d0 * v_inf_mag * SIN(delta / 2.0d0)
  END SUBROUTINE gravity_assist_dv

  ! ----------------------------------------------------------
  ! find_optimal_launch_date: serial 1-day sweep for best Δv.
  ! ----------------------------------------------------------
  SUBROUTINE find_optimal_launch_date(pos0, vel0, mass, target_idx, &
                                       jd_start, jd_end, dt_sim,    &
                                       best_jd, best_dv_SI)
    REAL(dp), INTENT(IN)  :: pos0(3, N_BODIES), vel0(3, N_BODIES)
    REAL(dp), INTENT(IN)  :: mass(N_BODIES)
    INTEGER,  INTENT(IN)  :: target_idx
    REAL(dp), INTENT(IN)  :: jd_start, jd_end, dt_sim
    REAL(dp), INTENT(OUT) :: best_jd, best_dv_SI

    REAL(dp) :: pos_c(3, N_BODIES), vel_c(3, N_BODIES)
    REAL(dp) :: r_earth(3), r_target(3), r1, r2
    REAL(dp) :: dv1, dv2, tf, cand_dv, jd_cur
    INTEGER  :: n_days, i

    pos_c  = pos0
    vel_c  = vel0
    n_days = INT(jd_end - jd_start)

    best_dv_SI = 1.0d15
    best_jd    = jd_start

    DO i = 1, n_days
      jd_cur = jd_start + REAL(i-1, dp)
      CALL leapfrog_step(N_BODIES, pos_c, vel_c, mass, dt_sim)

      r_earth = pos_c(:, EARTH_IDX) - pos_c(:, SUN_IDX)
      r1      = SQRT(DOT_PRODUCT(r_earth, r_earth))

      IF (target_idx == MARS_IDX) THEN
        r_target = pos_c(:, MARS_IDX) - pos_c(:, SUN_IDX)
        r2       = SQRT(DOT_PRODUCT(r_target, r_target))
        CALL hohmann_transfer(r1, r2, dv1, dv2, tf)
        cand_dv = dv1 + dv2
      ELSE
        cand_dv = 3900.0d0 + 200.0d0 * SIN(TWO_PI * jd_cur / 27.32d0)
      END IF

      IF (cand_dv < best_dv_SI) THEN
        best_dv_SI = cand_dv
        best_jd    = jd_cur
      END IF
    END DO
  END SUBROUTINE find_optimal_launch_date

  ! ----------------------------------------------------------
  ! time_trajectory: timing wrapper for thread sweep.
  ! ----------------------------------------------------------
  SUBROUTINE time_trajectory(n_pl, pos, vel, mass, sc_pos, sc_vel, &
                               sc_mass, target_idx, n_steps, dt, elapsed_sec)
    INTEGER,  INTENT(IN)  :: n_pl, target_idx, n_steps
    REAL(dp), INTENT(IN)  :: pos(3,n_pl), vel(3,n_pl), mass(n_pl)
    REAL(dp), INTENT(IN)  :: sc_pos(3), sc_vel(3), sc_mass, dt
    REAL(dp), INTENT(OUT) :: elapsed_sec

    INTEGER,  PARAMETER :: N_PATH = 100
    REAL(dp) :: sc_path(3, N_PATH)
    REAL(dp) :: land_lat, land_lon, land_jd, closest
    REAL(dp) :: t0
    LOGICAL  :: success

    CALL timer_start(t0)
    CALL propagate_spacecraft(n_pl, pos, vel, mass, sc_pos, sc_vel,  &
                               sc_mass, n_steps, dt, target_idx,     &
                               sc_path, N_PATH, land_lat, land_lon,  &
                               land_jd, closest, success, J2000_JD)
    CALL timer_stop(t0, elapsed_sec)
  END SUBROUTINE time_trajectory

  ! ----------------------------------------------------------
  ! run_trajectory_sweep: thread sweep for trajectory kernel.
  ! ----------------------------------------------------------
  SUBROUTINE run_trajectory_sweep(n_pl, pos, vel, mass)
    INTEGER,  INTENT(IN) :: n_pl
    REAL(dp), INTENT(IN) :: pos(3,n_pl), vel(3,n_pl), mass(n_pl)

    INTEGER,  PARAMETER  :: N_TRAJ_STEPS = 60
    REAL(dp), PARAMETER  :: DT_TRAJ      = 1.0d0

    REAL(dp) :: sc_pos(3), sc_vel(3), dummy_dv
    REAL(dp) :: run_times(N_SWEEP_REPEATS)
    REAL(dp) :: means(N_THREAD_COUNTS), stddevs(N_THREAD_COUNTS)
    REAL(dp) :: speedups(N_THREAD_COUNTS), efficiencies(N_THREAD_COUNTS)
    REAL(dp) :: t_serial, elapsed
    INTEGER  :: tc, r

    CALL compute_launch_vector(28.57d0, -80.65d0, J2000_JD, MOON_IDX, &
                                pos(:, EARTH_IDX), vel(:, EARTH_IDX),  &
                                sc_pos, sc_vel, dummy_dv)

    WRITE(*,'(A)') '  [SWEEP] Starting trajectory_prop thread sweep...'

    DO tc = 1, N_THREAD_COUNTS
      CALL omp_set_num_threads(THREAD_COUNTS(tc))
      DO r = 1, N_SWEEP_REPEATS
        CALL time_trajectory(n_pl, pos, vel, mass, sc_pos, sc_vel,    &
                              5000.0d0, MOON_IDX, N_TRAJ_STEPS, DT_TRAJ, elapsed)
        run_times(r) = elapsed
      END DO
      CALL compute_stats(run_times, N_SWEEP_REPEATS, means(tc), stddevs(tc))
    END DO

    t_serial = means(1)
    DO tc = 1, N_THREAD_COUNTS
      CALL compute_speedup(t_serial, means(tc), speedups(tc))
      CALL compute_efficiency(speedups(tc), THREAD_COUNTS(tc), efficiencies(tc))
      CALL append_sweep_csv('trajectory_prop', THREAD_COUNTS(tc), run_times, &
                            means(tc), stddevs(tc), speedups(tc), efficiencies(tc))
    END DO

    CALL print_sweep_table('trajectory_prop', THREAD_COUNTS, means, stddevs, &
                            speedups, efficiencies, N_THREAD_COUNTS)
  END SUBROUTINE run_trajectory_sweep

END MODULE trajectory_prop
