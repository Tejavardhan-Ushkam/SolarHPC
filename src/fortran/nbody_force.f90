! ============================================================
! MODULE: nbody_force
! PROJECT: SolarHPC
! PURPOSE: Gravitational force computation for N bodies.
!          Most computationally expensive kernel — primary
!          OpenMP parallelisation target.
!          Compiled THIRD. Depends on: solar_constants, omp_timer.
! ============================================================
MODULE nbody_force
  USE solar_constants
  USE omp_timer
  USE omp_lib
  IMPLICIT NONE

CONTAINS

  ! ----------------------------------------------------------
  ! compute_forces: O(N^2) direct summation with softening.
  ! OpenMP parallelised on outer i-loop (DYNAMIC scheduling).
  ! Thread count set externally via omp_set_num_threads().
  ! ----------------------------------------------------------
  SUBROUTINE compute_forces(n, pos, mass, acc)
    INTEGER,  INTENT(IN)  :: n
    REAL(dp), INTENT(IN)  :: pos(3, n)
    REAL(dp), INTENT(IN)  :: mass(n)
    REAL(dp), INTENT(OUT) :: acc(3, n)

    INTEGER  :: i, j
    REAL(dp) :: dx, dy, dz, r2, r3, factor
    REAL(dp) :: eps2

    eps2 = SOFTENING_EPS * SOFTENING_EPS

    acc(:, :) = 0.0d0

    !$OMP PARALLEL DO PRIVATE(i, j, dx, dy, dz, r2, r3, factor) &
    !$OMP&            SCHEDULE(DYNAMIC) DEFAULT(NONE)             &
    !$OMP&            SHARED(n, pos, mass, acc, eps2)
    DO i = 1, n
      DO j = 1, n
        IF (j == i) CYCLE
        dx = pos(1, j) - pos(1, i)
        dy = pos(2, j) - pos(2, i)
        dz = pos(3, j) - pos(3, i)
        r2 = dx*dx + dy*dy + dz*dz + eps2
        r3 = r2 * SQRT(r2)
        factor = G_AU_DAY * mass(j) / r3
        acc(1, i) = acc(1, i) + factor * dx
        acc(2, i) = acc(2, i) + factor * dy
        acc(3, i) = acc(3, i) + factor * dz
      END DO
    END DO
    !$OMP END PARALLEL DO
  END SUBROUTINE compute_forces

  ! ----------------------------------------------------------
  ! compute_forces_extended: same as compute_forces but
  ! handles n_ext bodies (n solar system + spacecraft).
  ! Spacecraft treated as massless (mass ~ 0 for force ON
  ! other bodies, but receives full gravitational pull).
  ! ----------------------------------------------------------
  SUBROUTINE compute_forces_extended(n_ext, pos, mass, acc)
    INTEGER,  INTENT(IN)  :: n_ext
    REAL(dp), INTENT(IN)  :: pos(3, n_ext)
    REAL(dp), INTENT(IN)  :: mass(n_ext)
    REAL(dp), INTENT(OUT) :: acc(3, n_ext)

    INTEGER  :: i, j
    REAL(dp) :: dx, dy, dz, r2, r3, factor
    REAL(dp) :: eps2

    eps2 = SOFTENING_EPS * SOFTENING_EPS
    acc(:, :) = 0.0d0

    !$OMP PARALLEL DO PRIVATE(i, j, dx, dy, dz, r2, r3, factor) &
    !$OMP&            SCHEDULE(DYNAMIC) DEFAULT(NONE)             &
    !$OMP&            SHARED(n_ext, pos, mass, acc, eps2)
    DO i = 1, n_ext
      DO j = 1, n_ext
        IF (j == i) CYCLE
        dx = pos(1, j) - pos(1, i)
        dy = pos(2, j) - pos(2, i)
        dz = pos(3, j) - pos(3, i)
        r2 = dx*dx + dy*dy + dz*dz + eps2
        r3 = r2 * SQRT(r2)
        ! Spacecraft (index > N_BODIES) has negligible mass effect on planets
        IF (j <= N_BODIES) THEN
          factor = G_AU_DAY * mass(j) / r3
        ELSE
          factor = 0.0d0
        END IF
        acc(1, i) = acc(1, i) + factor * dx
        acc(2, i) = acc(2, i) + factor * dy
        acc(3, i) = acc(3, i) + factor * dz
      END DO
    END DO
    !$OMP END PARALLEL DO
  END SUBROUTINE compute_forces_extended

  ! ----------------------------------------------------------
  ! verify_forces: unit test.
  ! Sun at origin, Earth at 1 AU. Check force magnitude.
  ! Expected: |a| = G_AU_DAY * M_sun / 1.0^2
  ! ----------------------------------------------------------
  SUBROUTINE verify_forces(passed)
    LOGICAL, INTENT(OUT) :: passed

    INTEGER,  PARAMETER :: nt = 2
    REAL(dp) :: pos(3, nt), mass_t(nt), acc(3, nt)
    REAL(dp) :: expected_mag, computed_mag, rel_err

    pos(:, 1)  = [0.0d0, 0.0d0, 0.0d0]   ! Sun at origin
    pos(:, 2)  = [1.0d0, 0.0d0, 0.0d0]   ! Earth at 1 AU
    mass_t(1)  = MASS_KG(SUN_IDX)
    mass_t(2)  = MASS_KG(EARTH_IDX)

    CALL compute_forces(nt, pos, mass_t, acc)

    expected_mag = G_AU_DAY * MASS_KG(SUN_IDX) / 1.0d0**2
    computed_mag = SQRT(acc(1,2)**2 + acc(2,2)**2 + acc(3,2)**2)
    rel_err      = ABS(computed_mag - expected_mag) / expected_mag

    passed = (rel_err < 1.0d-10)
    IF (passed) THEN
      WRITE(*,'(A,ES12.4)') '  [TEST] nbody verify_forces: PASS  rel_err=', rel_err
    ELSE
      WRITE(*,'(A,ES12.4)') '  [TEST] nbody verify_forces: FAIL  rel_err=', rel_err
    END IF
  END SUBROUTINE verify_forces

  ! ----------------------------------------------------------
  ! time_force_kernel: timing wrapper for thread sweep.
  ! Calls compute_forces n_repeat times.
  ! ----------------------------------------------------------
  SUBROUTINE time_force_kernel(n, pos, mass, n_repeat, elapsed_sec)
    INTEGER,  INTENT(IN)  :: n, n_repeat
    REAL(dp), INTENT(IN)  :: pos(3, n), mass(n)
    REAL(dp), INTENT(OUT) :: elapsed_sec

    REAL(dp) :: acc(3, n)
    REAL(dp) :: t0
    INTEGER  :: r

    CALL timer_start(t0)
    DO r = 1, n_repeat
      CALL compute_forces(n, pos, mass, acc)
    END DO
    CALL timer_stop(t0, elapsed_sec)
  END SUBROUTINE time_force_kernel

  ! ----------------------------------------------------------
  ! run_force_sweep: full thread-count sweep for force kernel.
  ! Tests N_THREAD_COUNTS thread counts, N_SWEEP_REPEATS times each.
  ! Writes rows to results/benchmarks/thread_sweep_results.csv.
  ! ----------------------------------------------------------
  SUBROUTINE run_force_sweep(n, pos, mass)
    INTEGER,  INTENT(IN) :: n
    REAL(dp), INTENT(IN) :: pos(3, n), mass(n)

    INTEGER,  PARAMETER  :: N_REPEAT_KERNEL = 20   ! repeats per timing call
    REAL(dp) :: run_times(N_SWEEP_REPEATS)
    REAL(dp) :: means(N_THREAD_COUNTS), stddevs(N_THREAD_COUNTS)
    REAL(dp) :: speedups(N_THREAD_COUNTS), efficiencies(N_THREAD_COUNTS)
    REAL(dp) :: t_serial, elapsed
    INTEGER  :: tc, r

    WRITE(*,'(A)') '  [SWEEP] Starting force_computation thread sweep...'

    DO tc = 1, N_THREAD_COUNTS
      CALL omp_set_num_threads(THREAD_COUNTS(tc))
      DO r = 1, N_SWEEP_REPEATS
        CALL time_force_kernel(n, pos, mass, N_REPEAT_KERNEL, elapsed)
        run_times(r) = elapsed / REAL(N_REPEAT_KERNEL, dp)
      END DO
      CALL compute_stats(run_times, N_SWEEP_REPEATS, means(tc), stddevs(tc))
    END DO

    t_serial = means(1)
    DO tc = 1, N_THREAD_COUNTS
      CALL compute_speedup(t_serial, means(tc), speedups(tc))
      CALL compute_efficiency(speedups(tc), THREAD_COUNTS(tc), efficiencies(tc))
      CALL append_sweep_csv('force_computation', THREAD_COUNTS(tc), &
                            run_times, means(tc), stddevs(tc), &
                            speedups(tc), efficiencies(tc))
    END DO

    CALL print_sweep_table('force_computation', THREAD_COUNTS, means, &
                            stddevs, speedups, efficiencies, N_THREAD_COUNTS)
  END SUBROUTINE run_force_sweep

END MODULE nbody_force
