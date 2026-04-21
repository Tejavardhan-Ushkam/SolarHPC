! ============================================================
! MODULE: eclipse_geom
! PROJECT: SolarHPC
! PURPOSE: Solar and lunar eclipse detection from trajectory
!          data. Shadow cone geometry. Saros cycle detection.
!          OpenMP-parallelised timestep scan.
!          Compiled SIXTH. Depends on: solar_constants,
!          omp_timer, coordinate_transform, omp_lib.
! ============================================================
MODULE eclipse_geom
  USE solar_constants
  USE omp_timer
  USE coordinate_transform
  USE omp_lib
  IMPLICIT NONE

  ! Eclipse event record
  TYPE :: EclipseEvent
    REAL(dp) :: jd
    INTEGER  :: eclipse_type    ! 0=solar, 1=lunar
    INTEGER  :: event_class     ! 1=partial, 2=total, 3=annular
    REAL(dp) :: umbra_fraction
    REAL(dp) :: duration_min
  END TYPE EclipseEvent

CONTAINS

  ! ----------------------------------------------------------
  ! check_solar_eclipse: shadow cone method.
  ! ----------------------------------------------------------
  SUBROUTINE check_solar_eclipse(pos, event_flag, umbra_fraction, duration_min)
    REAL(dp), INTENT(IN)  :: pos(3, N_BODIES)
    INTEGER,  INTENT(OUT) :: event_flag
    REAL(dp), INTENT(OUT) :: umbra_fraction, duration_min

    REAL(dp) :: r_sm(3), r_me(3)      ! Sun->Moon, Moon->Earth vectors
    REAL(dp) :: r_sm_mag, r_me_mag
    REAL(dp) :: alpha_umbra, alpha_pen
    REAL(dp) :: r_sm_unit(3)
    REAL(dp) :: cross_prod(3)
    REAL(dp) :: d_axis
    REAL(dp) :: r_shadow_umbra, r_shadow_pen
    REAL(dp) :: sun_ang_radius, moon_ang_radius

    event_flag     = 0
    umbra_fraction = 0.0d0
    duration_min   = 0.0d0

    r_sm = pos(:, SUN_IDX)  - pos(:, MOON_IDX)
    r_me = pos(:, EARTH_IDX) - pos(:, MOON_IDX)

    r_sm_mag = SQRT(DOT_PRODUCT(r_sm, r_sm))
    r_me_mag = SQRT(DOT_PRODUCT(r_me, r_me))

    IF (r_sm_mag < 1.0d-20 .OR. r_me_mag < 1.0d-20) RETURN

    r_sm_unit = r_sm / r_sm_mag

    ! Perpendicular distance from Earth to Moon-Sun axis
    cross_prod(1) = r_sm_unit(2)*r_me(3) - r_sm_unit(3)*r_me(2)
    cross_prod(2) = r_sm_unit(3)*r_me(1) - r_sm_unit(1)*r_me(3)
    cross_prod(3) = r_sm_unit(1)*r_me(2) - r_sm_unit(2)*r_me(1)
    d_axis = SQRT(DOT_PRODUCT(cross_prod, cross_prod))

    ! Shadow cone half-angles
    alpha_umbra = ASIN((SUN_RADIUS_AU - MOON_RADIUS_AU) / r_sm_mag)
    alpha_pen   = ASIN((SUN_RADIUS_AU + MOON_RADIUS_AU) / r_sm_mag)

    ! Shadow radii at Earth's distance
    r_shadow_umbra = MOON_RADIUS_AU - r_me_mag * TAN(MAX(alpha_umbra, 0.0d0))
    r_shadow_pen   = MOON_RADIUS_AU + r_me_mag * TAN(alpha_pen)

    ! Angular radii as seen from Earth
    sun_ang_radius  = ASIN(SUN_RADIUS_AU  / r_sm_mag)
    moon_ang_radius = ASIN(MOON_RADIUS_AU / r_me_mag)

    IF (d_axis < r_shadow_pen + EARTH_RADIUS_AU) THEN
      IF (d_axis < ABS(r_shadow_umbra)) THEN
        ! Umbra touches Earth
        IF (moon_ang_radius >= sun_ang_radius) THEN
          event_flag = 2   ! total
        ELSE
          event_flag = 3   ! annular
        END IF
        umbra_fraction = 1.0d0
      ELSE
        event_flag = 1   ! partial
        umbra_fraction = 1.0d0 - d_axis / (r_shadow_pen + EARTH_RADIUS_AU)
        umbra_fraction = MAX(0.0d0, MIN(1.0d0, umbra_fraction))
      END IF
      ! Rough duration estimate (minutes)
      duration_min = (2.0d0 * MOON_RADIUS_AU / 0.5d0) * DAY_TO_SEC / 60.0d0
    END IF
  END SUBROUTINE check_solar_eclipse

  ! ----------------------------------------------------------
  ! check_lunar_eclipse: Earth's shadow covers Moon.
  ! ----------------------------------------------------------
  SUBROUTINE check_lunar_eclipse(pos, event_flag, umbra_fraction, duration_min)
    REAL(dp), INTENT(IN)  :: pos(3, N_BODIES)
    INTEGER,  INTENT(OUT) :: event_flag
    REAL(dp), INTENT(OUT) :: umbra_fraction, duration_min

    REAL(dp) :: r_se(3), r_em(3)
    REAL(dp) :: r_se_mag, r_em_mag
    REAL(dp) :: r_se_unit(3), cross_prod(3)
    REAL(dp) :: d_axis, alpha_umbra, alpha_pen
    REAL(dp) :: r_umbra_at_moon, r_pen_at_moon

    event_flag     = 0
    umbra_fraction = 0.0d0
    duration_min   = 0.0d0

    r_se = pos(:, EARTH_IDX) - pos(:, SUN_IDX)
    r_em = pos(:, MOON_IDX)  - pos(:, EARTH_IDX)

    r_se_mag = SQRT(DOT_PRODUCT(r_se, r_se))
    r_em_mag = SQRT(DOT_PRODUCT(r_em, r_em))

    IF (r_se_mag < 1.0d-20 .OR. r_em_mag < 1.0d-20) RETURN

    r_se_unit = r_se / r_se_mag

    ! Distance from Moon to Sun-Earth axis
    cross_prod(1) = r_se_unit(2)*r_em(3) - r_se_unit(3)*r_em(2)
    cross_prod(2) = r_se_unit(3)*r_em(1) - r_se_unit(1)*r_em(3)
    cross_prod(3) = r_se_unit(1)*r_em(2) - r_se_unit(2)*r_em(1)
    d_axis = SQRT(DOT_PRODUCT(cross_prod, cross_prod))

    alpha_umbra = ASIN((SUN_RADIUS_AU - EARTH_RADIUS_AU) / r_se_mag)
    alpha_pen   = ASIN((SUN_RADIUS_AU + EARTH_RADIUS_AU) / r_se_mag)

    r_umbra_at_moon = EARTH_RADIUS_AU - r_em_mag * TAN(MAX(alpha_umbra, 0.0d0))
    r_pen_at_moon   = EARTH_RADIUS_AU + r_em_mag * TAN(alpha_pen)

    IF (d_axis < r_pen_at_moon + MOON_RADIUS_AU) THEN
      IF (d_axis < r_umbra_at_moon - MOON_RADIUS_AU) THEN
        event_flag     = 2
        umbra_fraction = 1.0d0
      ELSE IF (d_axis < r_umbra_at_moon + MOON_RADIUS_AU) THEN
        event_flag = 1
        umbra_fraction = (r_umbra_at_moon + MOON_RADIUS_AU - d_axis) / &
                         (2.0d0 * MOON_RADIUS_AU)
        umbra_fraction = MAX(0.0d0, MIN(1.0d0, umbra_fraction))
      ELSE
        event_flag = 1
        umbra_fraction = 0.1d0
      END IF
      duration_min = (2.0d0 * r_umbra_at_moon / 0.02d0) * DAY_TO_SEC / 60.0d0
      duration_min = MIN(duration_min, 200.0d0)
    END IF
  END SUBROUTINE check_lunar_eclipse

  ! ----------------------------------------------------------
  ! scan_eclipses_omp: OpenMP scan over stored trajectory.
  ! Writes results/eclipses/eclipse_predictions.csv.
  ! ----------------------------------------------------------
  SUBROUTINE scan_eclipses_omp(n_steps, pos_history, dt_days, t_start_jd, &
                                n_solar, n_lunar)
    INTEGER,  INTENT(IN)  :: n_steps
    REAL(dp), INTENT(IN)  :: pos_history(3, N_BODIES, n_steps)
    REAL(dp), INTENT(IN)  :: dt_days, t_start_jd
    INTEGER,  INTENT(OUT) :: n_solar, n_lunar

    INTEGER,  PARAMETER  :: MAX_EVENTS = 5000
    TYPE(EclipseEvent)   :: events(MAX_EVENTS)
    INTEGER              :: n_events

    INTEGER  :: i, flag_s, flag_l, unit_csv, ios
    REAL(dp) :: uf_s, uf_l, dur_s, dur_l, jd_i
    REAL(dp) :: last_solar_jd, last_lunar_jd
    INTEGER  :: yr, mo, dy, hr, mn
    REAL(dp) :: sc
    CHARACTER(LEN=20) :: date_str

    n_events      = 0
    n_solar       = 0
    n_lunar       = 0
    last_solar_jd = -1.0d10
    last_lunar_jd = -1.0d10

    ! Note: the loop is over timesteps. Eclipse detection is fast
    ! (pure arithmetic), so we parallelise but use CRITICAL for output.
    !$OMP PARALLEL DO PRIVATE(i, flag_s, flag_l, uf_s, uf_l, dur_s, dur_l, jd_i) &
    !$OMP&            SCHEDULE(STATIC) DEFAULT(SHARED)
    DO i = 1, n_steps
      jd_i = t_start_jd + REAL(i-1, dp) * dt_days

      CALL check_solar_eclipse(pos_history(:,:,i), flag_s, uf_s, dur_s)
      CALL check_lunar_eclipse(pos_history(:,:,i), flag_l, uf_l, dur_l)

      !$OMP CRITICAL(eclipse_write)
      IF (flag_s > 0 .AND. (jd_i - last_solar_jd) > MIN_ECLIPSE_SEPARATION_DAYS) THEN
        IF (n_events < MAX_EVENTS) THEN
          n_events = n_events + 1
          events(n_events)%jd            = jd_i
          events(n_events)%eclipse_type  = 0
          events(n_events)%event_class   = flag_s
          events(n_events)%umbra_fraction= uf_s
          events(n_events)%duration_min  = dur_s
          n_solar       = n_solar + 1
          last_solar_jd = jd_i
        END IF
      END IF
      IF (flag_l > 0 .AND. (jd_i - last_lunar_jd) > MIN_ECLIPSE_SEPARATION_DAYS) THEN
        IF (n_events < MAX_EVENTS) THEN
          n_events = n_events + 1
          events(n_events)%jd            = jd_i
          events(n_events)%eclipse_type  = 1
          events(n_events)%event_class   = flag_l
          events(n_events)%umbra_fraction= uf_l
          events(n_events)%duration_min  = dur_l
          n_lunar       = n_lunar + 1
          last_lunar_jd = jd_i
        END IF
      END IF
      !$OMP END CRITICAL(eclipse_write)
    END DO
    !$OMP END PARALLEL DO

    ! Write CSV
    unit_csv = 55
    OPEN(UNIT=unit_csv, FILE='results/eclipses/eclipse_predictions.csv', &
         STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
    IF (ios /= 0) THEN
      WRITE(*,*) '[ERROR] Cannot write eclipse_predictions.csv'
      RETURN
    END IF
    WRITE(unit_csv,'(A)') &
      'julian_day,gregorian_date,eclipse_type,event_class,umbra_fraction,duration_minutes'

    DO i = 1, n_events
      CALL jd_to_gregorian(events(i)%jd, yr, mo, dy, hr, mn, sc)
      WRITE(date_str,'(I4,A,I2.2,A,I2.2)') yr, '-', mo, '-', dy
      WRITE(unit_csv,'(F14.2,A,A,A,I1,A,I1,A,F8.4,A,F8.1)') &
        events(i)%jd, ',', TRIM(date_str), ',', &
        events(i)%eclipse_type, ',', events(i)%event_class, ',', &
        events(i)%umbra_fraction, ',', events(i)%duration_min
    END DO
    CLOSE(unit_csv)
  END SUBROUTINE scan_eclipses_omp

  ! ----------------------------------------------------------
  ! verify_eclipse_2024: check 2024-Apr-08 total solar eclipse
  ! ----------------------------------------------------------
  SUBROUTINE verify_eclipse_2024(eclipse_file, passed)
    CHARACTER(LEN=*), INTENT(IN)  :: eclipse_file
    LOGICAL,          INTENT(OUT) :: passed

    REAL(dp), PARAMETER :: TARGET_JD  = 2460409.0d0   ! 2024-Apr-08
    REAL(dp), PARAMETER :: TOLERANCE  = 2.0d0          ! days

    INTEGER  :: unit_in, ios, etype, eclass
    REAL(dp) :: jd_read, uf, dur
    CHARACTER(LEN=64) :: line, date_str

    passed   = .FALSE.
    unit_in  = 56

    OPEN(UNIT=unit_in, FILE=TRIM(eclipse_file), STATUS='OLD', &
         ACTION='READ', IOSTAT=ios)
    IF (ios /= 0) THEN
      WRITE(*,'(A)') '  [TEST] Eclipse 2024-04-08: SKIP (no eclipse file found)'
      RETURN
    END IF

    READ(unit_in, '(A)', IOSTAT=ios) line   ! skip header
    DO
      READ(unit_in, *, IOSTAT=ios) jd_read, date_str, etype, eclass, uf, dur
      IF (ios /= 0) EXIT
      IF (etype == 0 .AND. ABS(jd_read - TARGET_JD) < TOLERANCE) THEN
        passed = .TRUE.
        EXIT
      END IF
    END DO
    CLOSE(unit_in)

    IF (passed) THEN
      WRITE(*,'(A,F10.1,A,I1)') '  [TEST] Eclipse 2024-04-08: PASS  found at JD=', &
        jd_read, '  class=', eclass
    ELSE
      WRITE(*,'(A)') '  [TEST] Eclipse 2024-04-08: FAIL (not found in eclipse list)'
    END IF
  END SUBROUTINE verify_eclipse_2024

  ! ----------------------------------------------------------
  ! time_eclipse_scan: timing wrapper for thread sweep.
  ! ----------------------------------------------------------
  SUBROUTINE time_eclipse_scan(n_steps, pos_history, dt_days, t_start_jd, &
                                elapsed_sec)
    INTEGER,  INTENT(IN)  :: n_steps
    REAL(dp), INTENT(IN)  :: pos_history(3, N_BODIES, n_steps)
    REAL(dp), INTENT(IN)  :: dt_days, t_start_jd
    REAL(dp), INTENT(OUT) :: elapsed_sec

    INTEGER  :: n_s, n_l
    REAL(dp) :: t0

    CALL timer_start(t0)
    CALL scan_eclipses_omp(n_steps, pos_history, dt_days, t_start_jd, n_s, n_l)
    CALL timer_stop(t0, elapsed_sec)
  END SUBROUTINE time_eclipse_scan

  ! ----------------------------------------------------------
  ! run_eclipse_sweep: thread sweep for eclipse scan kernel.
  ! ----------------------------------------------------------
  SUBROUTINE run_eclipse_sweep(n_steps, pos_history, dt_days, t_start_jd)
    INTEGER,  INTENT(IN)  :: n_steps
    REAL(dp), INTENT(IN)  :: pos_history(3, N_BODIES, n_steps)
    REAL(dp), INTENT(IN)  :: dt_days, t_start_jd

    REAL(dp) :: run_times(N_SWEEP_REPEATS)
    REAL(dp) :: means(N_THREAD_COUNTS), stddevs(N_THREAD_COUNTS)
    REAL(dp) :: speedups(N_THREAD_COUNTS), efficiencies(N_THREAD_COUNTS)
    REAL(dp) :: t_serial, elapsed
    INTEGER  :: tc, r

    WRITE(*,'(A)') '  [SWEEP] Starting eclipse_scan thread sweep...'

    DO tc = 1, N_THREAD_COUNTS
      CALL omp_set_num_threads(THREAD_COUNTS(tc))
      DO r = 1, N_SWEEP_REPEATS
        CALL time_eclipse_scan(n_steps, pos_history, dt_days, t_start_jd, elapsed)
        run_times(r) = elapsed
      END DO
      CALL compute_stats(run_times, N_SWEEP_REPEATS, means(tc), stddevs(tc))
    END DO

    t_serial = means(1)
    DO tc = 1, N_THREAD_COUNTS
      CALL compute_speedup(t_serial, means(tc), speedups(tc))
      CALL compute_efficiency(speedups(tc), THREAD_COUNTS(tc), efficiencies(tc))
      CALL append_sweep_csv('eclipse_scan', THREAD_COUNTS(tc), run_times, &
                            means(tc), stddevs(tc), speedups(tc), efficiencies(tc))
    END DO

    CALL print_sweep_table('eclipse_scan', THREAD_COUNTS, means, stddevs, &
                            speedups, efficiencies, N_THREAD_COUNTS)
  END SUBROUTINE run_eclipse_sweep

END MODULE eclipse_geom
