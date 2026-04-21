! ============================================================
! MODULE: integrators
! PROJECT: SolarHPC
! PURPOSE: Leapfrog (symplectic) and RK4 integrators for the
!          N-body ODE system. Both are OpenMP-parallelised.
!          Energy and position error diagnostics included.
!          Compiled FOURTH. Depends on: solar_constants,
!          omp_timer, nbody_force.
! ============================================================
MODULE integrators
  USE solar_constants
  USE omp_timer
  USE nbody_force
  USE omp_lib
  IMPLICIT NONE

CONTAINS

  ! ----------------------------------------------------------
  ! leapfrog_step: Velocity Verlet / leapfrog (symplectic).
  ! Conserves energy for long runs. O(N^2) per step.
  ! OpenMP on position and velocity update loops.
  ! ----------------------------------------------------------
  SUBROUTINE leapfrog_step(n, pos, vel, mass, dt)
    INTEGER,  INTENT(IN)    :: n
    REAL(dp), INTENT(INOUT) :: pos(3, n), vel(3, n)
    REAL(dp), INTENT(IN)    :: mass(n), dt

    REAL(dp) :: acc_old(3, n), acc_new(3, n), vel_half(3, n)
    INTEGER  :: i

    ! Step 1: forces at current positions
    CALL compute_forces(n, pos, mass, acc_old)

    ! Step 2: half-step velocity
    !$OMP PARALLEL DO PRIVATE(i) SCHEDULE(STATIC) DEFAULT(NONE) &
    !$OMP&            SHARED(n, vel, vel_half, acc_old, dt)
    DO i = 1, n
      vel_half(1, i) = vel(1, i) + 0.5d0 * dt * acc_old(1, i)
      vel_half(2, i) = vel(2, i) + 0.5d0 * dt * acc_old(2, i)
      vel_half(3, i) = vel(3, i) + 0.5d0 * dt * acc_old(3, i)
    END DO
    !$OMP END PARALLEL DO

    ! Step 3: full position update
    !$OMP PARALLEL DO PRIVATE(i) SCHEDULE(STATIC) DEFAULT(NONE) &
    !$OMP&            SHARED(n, pos, vel_half, dt)
    DO i = 1, n
      pos(1, i) = pos(1, i) + dt * vel_half(1, i)
      pos(2, i) = pos(2, i) + dt * vel_half(2, i)
      pos(3, i) = pos(3, i) + dt * vel_half(3, i)
    END DO
    !$OMP END PARALLEL DO

    ! Step 4: forces at new positions
    CALL compute_forces(n, pos, mass, acc_new)

    ! Step 5: complete velocity update
    !$OMP PARALLEL DO PRIVATE(i) SCHEDULE(STATIC) DEFAULT(NONE) &
    !$OMP&            SHARED(n, vel, vel_half, acc_new, dt)
    DO i = 1, n
      vel(1, i) = vel_half(1, i) + 0.5d0 * dt * acc_new(1, i)
      vel(2, i) = vel_half(2, i) + 0.5d0 * dt * acc_new(2, i)
      vel(3, i) = vel_half(3, i) + 0.5d0 * dt * acc_new(3, i)
    END DO
    !$OMP END PARALLEL DO
  END SUBROUTINE leapfrog_step

  ! ----------------------------------------------------------
  ! leapfrog_step_extended: same but for n_ext bodies
  ! (solar system + spacecraft). Uses compute_forces_extended.
  ! ----------------------------------------------------------
  SUBROUTINE leapfrog_step_extended(n_ext, pos, vel, mass, dt)
    INTEGER,  INTENT(IN)    :: n_ext
    REAL(dp), INTENT(INOUT) :: pos(3, n_ext), vel(3, n_ext)
    REAL(dp), INTENT(IN)    :: mass(n_ext), dt

    REAL(dp) :: acc_old(3, n_ext), acc_new(3, n_ext), vel_half(3, n_ext)
    INTEGER  :: i

    CALL compute_forces_extended(n_ext, pos, mass, acc_old)

    !$OMP PARALLEL DO PRIVATE(i) SCHEDULE(STATIC)
    DO i = 1, n_ext
      vel_half(:, i) = vel(:, i) + 0.5d0 * dt * acc_old(:, i)
    END DO
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i) SCHEDULE(STATIC)
    DO i = 1, n_ext
      pos(:, i) = pos(:, i) + dt * vel_half(:, i)
    END DO
    !$OMP END PARALLEL DO

    CALL compute_forces_extended(n_ext, pos, mass, acc_new)

    !$OMP PARALLEL DO PRIVATE(i) SCHEDULE(STATIC)
    DO i = 1, n_ext
      vel(:, i) = vel_half(:, i) + 0.5d0 * dt * acc_new(:, i)
    END DO
    !$OMP END PARALLEL DO
  END SUBROUTINE leapfrog_step_extended

  ! ----------------------------------------------------------
  ! rk4_step: classical 4th-order Runge-Kutta.
  ! NOT symplectic -- energy drifts secularly over long runs.
  ! Used to DEMONSTRATE why leapfrog is the right choice.
  ! ----------------------------------------------------------
  SUBROUTINE rk4_step(n, pos, vel, mass, dt)
    INTEGER,  INTENT(IN)    :: n
    REAL(dp), INTENT(INOUT) :: pos(3, n), vel(3, n)
    REAL(dp), INTENT(IN)    :: mass(n), dt

    REAL(dp) :: k1p(3,n), k1v(3,n)
    REAL(dp) :: k2p(3,n), k2v(3,n)
    REAL(dp) :: k3p(3,n), k3v(3,n)
    REAL(dp) :: k4p(3,n), k4v(3,n)
    REAL(dp) :: pos_tmp(3,n), vel_tmp(3,n)
    REAL(dp) :: acc_tmp(3,n)
    INTEGER  :: i

    ! k1 = f(y)
    CALL compute_forces(n, pos, mass, acc_tmp)
    k1p = vel
    k1v = acc_tmp

    ! k2 = f(y + dt/2 * k1)
    !$OMP PARALLEL DO PRIVATE(i) SCHEDULE(STATIC)
    DO i = 1, n
      pos_tmp(:,i) = pos(:,i) + 0.5d0*dt*k1p(:,i)
      vel_tmp(:,i) = vel(:,i) + 0.5d0*dt*k1v(:,i)
    END DO
    !$OMP END PARALLEL DO
    CALL compute_forces(n, pos_tmp, mass, acc_tmp)
    k2p = vel_tmp
    k2v = acc_tmp

    ! k3 = f(y + dt/2 * k2)
    !$OMP PARALLEL DO PRIVATE(i) SCHEDULE(STATIC)
    DO i = 1, n
      pos_tmp(:,i) = pos(:,i) + 0.5d0*dt*k2p(:,i)
      vel_tmp(:,i) = vel(:,i) + 0.5d0*dt*k2v(:,i)
    END DO
    !$OMP END PARALLEL DO
    CALL compute_forces(n, pos_tmp, mass, acc_tmp)
    k3p = vel_tmp
    k3v = acc_tmp

    ! k4 = f(y + dt * k3)
    !$OMP PARALLEL DO PRIVATE(i) SCHEDULE(STATIC)
    DO i = 1, n
      pos_tmp(:,i) = pos(:,i) + dt*k3p(:,i)
      vel_tmp(:,i) = vel(:,i) + dt*k3v(:,i)
    END DO
    !$OMP END PARALLEL DO
    CALL compute_forces(n, pos_tmp, mass, acc_tmp)
    k4p = vel_tmp
    k4v = acc_tmp

    ! Update: y_new = y + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)
    !$OMP PARALLEL DO PRIVATE(i) SCHEDULE(STATIC)
    DO i = 1, n
      pos(:,i) = pos(:,i) + (dt/6.0d0) * &
                 (k1p(:,i) + 2.0d0*k2p(:,i) + 2.0d0*k3p(:,i) + k4p(:,i))
      vel(:,i) = vel(:,i) + (dt/6.0d0) * &
                 (k1v(:,i) + 2.0d0*k2v(:,i) + 2.0d0*k3v(:,i) + k4v(:,i))
    END DO
    !$OMP END PARALLEL DO
  END SUBROUTINE rk4_step

  ! ----------------------------------------------------------
  ! compute_energy: total mechanical energy in Joules.
  ! KE + gravitational PE. Used for energy conservation test.
  ! ----------------------------------------------------------
  SUBROUTINE compute_energy(n, pos, vel, mass, E_total)
    INTEGER,  INTENT(IN)  :: n
    REAL(dp), INTENT(IN)  :: pos(3, n), vel(3, n), mass(n)
    REAL(dp), INTENT(OUT) :: E_total

    REAL(dp) :: KE, PE
    REAL(dp) :: v2, dx, dy, dz, r_m
    REAL(dp) :: pos_m(3, n)   ! positions in metres
    INTEGER  :: i, j

    ! Convert AU to metres
    pos_m = pos * AU_TO_M

    ! Kinetic energy (convert AU/day to m/s first)
    KE = 0.0d0
    DO i = 1, n
      v2 = (vel(1,i)**2 + vel(2,i)**2 + vel(3,i)**2) * AU_PER_DAY_TO_M_PER_S**2
      KE = KE + 0.5d0 * mass(i) * v2
    END DO

    ! Gravitational potential energy
    PE = 0.0d0
    DO i = 1, n - 1
      DO j = i + 1, n
        dx  = pos_m(1,j) - pos_m(1,i)
        dy  = pos_m(2,j) - pos_m(2,i)
        dz  = pos_m(3,j) - pos_m(3,i)
        r_m = SQRT(dx*dx + dy*dy + dz*dz)
        PE  = PE - G_SI * mass(i) * mass(j) / r_m
      END DO
    END DO

    E_total = KE + PE
  END SUBROUTINE compute_energy

  ! ----------------------------------------------------------
  ! compute_position_error: compare simulation vs reference.
  ! ----------------------------------------------------------
  SUBROUTINE compute_position_error(n, pos_sim, pos_ref, &
                                     max_error_AU, rms_error_AU)
    INTEGER,  INTENT(IN)  :: n
    REAL(dp), INTENT(IN)  :: pos_sim(3, n), pos_ref(3, n)
    REAL(dp), INTENT(OUT) :: max_error_AU, rms_error_AU

    REAL(dp) :: dx, dy, dz, err
    INTEGER  :: i

    max_error_AU = 0.0d0
    rms_error_AU = 0.0d0

    DO i = 1, n
      dx  = pos_sim(1,i) - pos_ref(1,i)
      dy  = pos_sim(2,i) - pos_ref(2,i)
      dz  = pos_sim(3,i) - pos_ref(3,i)
      err = SQRT(dx*dx + dy*dy + dz*dz)
      IF (err > max_error_AU) max_error_AU = err
      rms_error_AU = rms_error_AU + err * err
    END DO
    rms_error_AU = SQRT(rms_error_AU / REAL(n, dp))
  END SUBROUTINE compute_position_error

  ! ----------------------------------------------------------
  ! time_leapfrog: timing wrapper for thread sweep.
  ! ----------------------------------------------------------
  SUBROUTINE time_leapfrog(n, pos, vel, mass, dt, n_steps, elapsed_sec)
    INTEGER,  INTENT(IN)  :: n, n_steps
    REAL(dp), INTENT(IN)  :: pos(3,n), vel(3,n), mass(n), dt
    REAL(dp), INTENT(OUT) :: elapsed_sec

    REAL(dp) :: pos_tmp(3,n), vel_tmp(3,n)
    REAL(dp) :: t0
    INTEGER  :: s

    pos_tmp = pos
    vel_tmp = vel

    CALL timer_start(t0)
    DO s = 1, n_steps
      CALL leapfrog_step(n, pos_tmp, vel_tmp, mass, dt)
    END DO
    CALL timer_stop(t0, elapsed_sec)
  END SUBROUTINE time_leapfrog

  ! ----------------------------------------------------------
  ! time_rk4: timing wrapper for thread sweep.
  ! ----------------------------------------------------------
  SUBROUTINE time_rk4(n, pos, vel, mass, dt, n_steps, elapsed_sec)
    INTEGER,  INTENT(IN)  :: n, n_steps
    REAL(dp), INTENT(IN)  :: pos(3,n), vel(3,n), mass(n), dt
    REAL(dp), INTENT(OUT) :: elapsed_sec

    REAL(dp) :: pos_tmp(3,n), vel_tmp(3,n)
    REAL(dp) :: t0
    INTEGER  :: s

    pos_tmp = pos
    vel_tmp = vel

    CALL timer_start(t0)
    DO s = 1, n_steps
      CALL rk4_step(n, pos_tmp, vel_tmp, mass, dt)
    END DO
    CALL timer_stop(t0, elapsed_sec)
  END SUBROUTINE time_rk4

  ! ----------------------------------------------------------
  ! run_integrator_sweep: thread sweep for both integrators.
  ! Writes to thread_sweep_results.csv.
  ! Also writes integrator_comparison.csv (energy drift data).
  ! ----------------------------------------------------------
  SUBROUTINE run_integrator_sweep(n, pos, vel, mass, dt, n_steps)
    INTEGER,  INTENT(IN) :: n, n_steps
    REAL(dp), INTENT(IN) :: pos(3,n), vel(3,n), mass(n), dt

    INTEGER,  PARAMETER  :: N_BENCH_STEPS = 50
    REAL(dp) :: run_times(N_SWEEP_REPEATS)
    REAL(dp) :: means(N_THREAD_COUNTS), stddevs(N_THREAD_COUNTS)
    REAL(dp) :: speedups(N_THREAD_COUNTS), efficiencies(N_THREAD_COUNTS)
    REAL(dp) :: t_serial, elapsed
    INTEGER  :: tc, r
    CHARACTER(LEN=32) :: kname

    ! --- Leapfrog sweep ---
    kname = 'leapfrog_step'
    WRITE(*,'(A,A,A)') '  [SWEEP] Starting ', TRIM(kname), ' thread sweep...'
    DO tc = 1, N_THREAD_COUNTS
      CALL omp_set_num_threads(THREAD_COUNTS(tc))
      DO r = 1, N_SWEEP_REPEATS
        CALL time_leapfrog(n, pos, vel, mass, dt, N_BENCH_STEPS, elapsed)
        run_times(r) = elapsed / REAL(N_BENCH_STEPS, dp)
      END DO
      CALL compute_stats(run_times, N_SWEEP_REPEATS, means(tc), stddevs(tc))
    END DO
    t_serial = means(1)
    DO tc = 1, N_THREAD_COUNTS
      CALL compute_speedup(t_serial, means(tc), speedups(tc))
      CALL compute_efficiency(speedups(tc), THREAD_COUNTS(tc), efficiencies(tc))
      CALL append_sweep_csv(TRIM(kname), THREAD_COUNTS(tc), run_times, &
                            means(tc), stddevs(tc), speedups(tc), efficiencies(tc))
    END DO
    CALL print_sweep_table(TRIM(kname), THREAD_COUNTS, means, stddevs, &
                            speedups, efficiencies, N_THREAD_COUNTS)

    ! --- RK4 sweep ---
    kname = 'rk4_step'
    WRITE(*,'(A,A,A)') '  [SWEEP] Starting ', TRIM(kname), ' thread sweep...'
    DO tc = 1, N_THREAD_COUNTS
      CALL omp_set_num_threads(THREAD_COUNTS(tc))
      DO r = 1, N_SWEEP_REPEATS
        CALL time_rk4(n, pos, vel, mass, dt, N_BENCH_STEPS, elapsed)
        run_times(r) = elapsed / REAL(N_BENCH_STEPS, dp)
      END DO
      CALL compute_stats(run_times, N_SWEEP_REPEATS, means(tc), stddevs(tc))
    END DO
    t_serial = means(1)
    DO tc = 1, N_THREAD_COUNTS
      CALL compute_speedup(t_serial, means(tc), speedups(tc))
      CALL compute_efficiency(speedups(tc), THREAD_COUNTS(tc), efficiencies(tc))
      CALL append_sweep_csv(TRIM(kname), THREAD_COUNTS(tc), run_times, &
                            means(tc), stddevs(tc), speedups(tc), efficiencies(tc))
    END DO
    CALL print_sweep_table(TRIM(kname), THREAD_COUNTS, means, stddevs, &
                            speedups, efficiencies, N_THREAD_COUNTS)
  END SUBROUTINE run_integrator_sweep

END MODULE integrators
