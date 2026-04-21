! ============================================================
! MODULE: omp_timer
! PROJECT: SolarHPC
! PURPOSE: Timing, statistics, Amdahl fitting, and formatted
!          table printing for all benchmarking code.
!          Compiled SECOND. Depends on: solar_constants.
! ============================================================
MODULE omp_timer
  USE solar_constants
  USE omp_lib
  IMPLICIT NONE

CONTAINS

  ! ----------------------------------------------------------
  ! Record current wall-clock time
  ! ----------------------------------------------------------
  SUBROUTINE timer_start(t0)
    REAL(dp), INTENT(OUT) :: t0
    t0 = OMP_GET_WTIME()
  END SUBROUTINE timer_start

  ! ----------------------------------------------------------
  ! Compute elapsed time since t0
  ! ----------------------------------------------------------
  SUBROUTINE timer_stop(t0, elapsed)
    REAL(dp), INTENT(IN)  :: t0
    REAL(dp), INTENT(OUT) :: elapsed
    elapsed = OMP_GET_WTIME() - t0
  END SUBROUTINE timer_stop

  ! ----------------------------------------------------------
  ! Print a labelled timing result
  ! ----------------------------------------------------------
  SUBROUTINE timer_print(label, elapsed)
    CHARACTER(LEN=*), INTENT(IN) :: label
    REAL(dp),         INTENT(IN) :: elapsed
    WRITE(*,'(A,A30,A,F14.8,A)') '  [TIMER] ', TRIM(label), ' : ', elapsed, ' sec'
  END SUBROUTINE timer_print

  ! ----------------------------------------------------------
  ! Sample mean and standard deviation (n-1 denominator)
  ! ----------------------------------------------------------
  SUBROUTINE compute_stats(times, n, mean_t, stddev_t)
    INTEGER,  INTENT(IN)  :: n
    REAL(dp), INTENT(IN)  :: times(n)
    REAL(dp), INTENT(OUT) :: mean_t, stddev_t
    REAL(dp) :: sumv, sumsq, diff
    INTEGER  :: i

    sumv = 0.0d0
    DO i = 1, n
      sumv = sumv + times(i)
    END DO
    mean_t = sumv / REAL(n, dp)

    IF (n < 2) THEN
      stddev_t = 0.0d0
      RETURN
    END IF

    sumsq = 0.0d0
    DO i = 1, n
      diff  = times(i) - mean_t
      sumsq = sumsq + diff * diff
    END DO
    stddev_t = SQRT(sumsq / REAL(n - 1, dp))
  END SUBROUTINE compute_stats

  ! ----------------------------------------------------------
  ! Speedup: S = T_serial / T_parallel  (guarded against zero)
  ! ----------------------------------------------------------
  SUBROUTINE compute_speedup(t_serial, t_parallel, speedup)
    REAL(dp), INTENT(IN)  :: t_serial, t_parallel
    REAL(dp), INTENT(OUT) :: speedup
    speedup = t_serial / MAX(t_parallel, 1.0d-14)
  END SUBROUTINE compute_speedup

  ! ----------------------------------------------------------
  ! Parallel efficiency E = S/N * 100  (percent)
  ! ----------------------------------------------------------
  SUBROUTINE compute_efficiency(speedup, n_threads, efficiency_pct)
    REAL(dp), INTENT(IN)  :: speedup
    INTEGER,  INTENT(IN)  :: n_threads
    REAL(dp), INTENT(OUT) :: efficiency_pct
    efficiency_pct = speedup / REAL(n_threads, dp) * 100.0d0
  END SUBROUTINE compute_efficiency

  ! ----------------------------------------------------------
  ! Fit Amdahl's Law: S(N) = 1 / (f + (1-f)/N)
  ! Grid search over f in [0.001, 0.999] minimising SSR
  ! ----------------------------------------------------------
  SUBROUTINE fit_amdahl(thread_counts, speedups, n_points, serial_frac)
    INTEGER,  INTENT(IN)  :: n_points
    INTEGER,  INTENT(IN)  :: thread_counts(n_points)
    REAL(dp), INTENT(IN)  :: speedups(n_points)
    REAL(dp), INTENT(OUT) :: serial_frac

    REAL(dp) :: f, best_f, ssr, best_ssr, predicted, residual
    INTEGER  :: k, j
    INTEGER  :: nf
    REAL(dp) :: df

    nf       = 999
    df       = 0.001d0
    best_f   = 0.5d0
    best_ssr = 1.0d30

    DO k = 1, nf
      f   = REAL(k, dp) * df
      ssr = 0.0d0
      DO j = 1, n_points
        predicted = 1.0d0 / (f + (1.0d0 - f) / REAL(thread_counts(j), dp))
        residual  = speedups(j) - predicted
        ssr       = ssr + residual * residual
      END DO
      IF (ssr < best_ssr) THEN
        best_ssr = ssr
        best_f   = f
      END IF
    END DO

    serial_frac = best_f
  END SUBROUTINE fit_amdahl

  ! ----------------------------------------------------------
  ! Print a formatted ASCII sweep table to stdout
  ! ----------------------------------------------------------
  SUBROUTINE print_sweep_table(kernel_name, thread_counts, means, &
                                stddevs, speedups, efficiencies, n)
    CHARACTER(LEN=*), INTENT(IN) :: kernel_name
    INTEGER,          INTENT(IN) :: n
    INTEGER,          INTENT(IN) :: thread_counts(n)
    REAL(dp),         INTENT(IN) :: means(n), stddevs(n), speedups(n), efficiencies(n)
    INTEGER :: i

    WRITE(*,*)
    WRITE(*,'(A,A)') '  Kernel: ', TRIM(kernel_name)
    WRITE(*,'(A)') '  +----------+-----------+-----------+---------+------------+'
    WRITE(*,'(A)') '  | Threads  | Mean (s)  | Std (s)   | Speedup | Effic. (%) |'
    WRITE(*,'(A)') '  +----------+-----------+-----------+---------+------------+'
    DO i = 1, n
      WRITE(*,'(A,I5,A,F10.6,A,F10.6,A,F8.3,A,F9.1,A)') &
        '  |', thread_counts(i), '     |', means(i), ' |', &
        stddevs(i), ' |', speedups(i), ' |', efficiencies(i), '     |'
    END DO
    WRITE(*,'(A)') '  +----------+-----------+-----------+---------+------------+'
    WRITE(*,*)
  END SUBROUTINE print_sweep_table

  ! ----------------------------------------------------------
  ! Append one row to thread_sweep_results.csv
  ! FORMAT matches contract exactly:
  ! kernel_name,n_threads,run1..run5,mean,stddev,speedup,efficiency
  ! ----------------------------------------------------------
  SUBROUTINE append_sweep_csv(kernel_name, n_threads, run_times, &
                               mean_t, stddev_t, speedup, efficiency)
    CHARACTER(LEN=*), INTENT(IN) :: kernel_name
    INTEGER,          INTENT(IN) :: n_threads
    REAL(dp),         INTENT(IN) :: run_times(N_SWEEP_REPEATS)
    REAL(dp),         INTENT(IN) :: mean_t, stddev_t, speedup, efficiency

    INTEGER            :: ios, unit_csv
    LOGICAL            :: file_exists
    CHARACTER(LEN=256) :: filepath

    filepath = 'results/benchmarks/thread_sweep_results.csv'
    unit_csv = 42

    INQUIRE(FILE=TRIM(filepath), EXIST=file_exists)
    IF (.NOT. file_exists) THEN
      OPEN(UNIT=unit_csv, FILE=TRIM(filepath), STATUS='NEW', &
           ACTION='WRITE', IOSTAT=ios)
      IF (ios /= 0) THEN
        WRITE(*,*) '[ERROR] Cannot create ', TRIM(filepath)
        RETURN
      END IF
      WRITE(unit_csv,'(A)') &
        'kernel_name,n_threads,run1_sec,run2_sec,run3_sec,run4_sec,run5_sec,'// &
        'mean_sec,stddev_sec,speedup,efficiency_pct'
      CLOSE(unit_csv)
    END IF

    OPEN(UNIT=unit_csv, FILE=TRIM(filepath), STATUS='OLD', &
         ACTION='WRITE', POSITION='APPEND', IOSTAT=ios)
    IF (ios /= 0) THEN
      WRITE(*,*) '[ERROR] Cannot append to ', TRIM(filepath)
      RETURN
    END IF
    WRITE(unit_csv,'(A,A,I3,5(A,F14.8),2(A,F14.8),A,F10.4,A,F8.2)') &
      TRIM(kernel_name), ',', n_threads, &
      ',', run_times(1), ',', run_times(2), ',', run_times(3), &
      ',', run_times(4), ',', run_times(5), &
      ',', mean_t, ',', stddev_t, ',', speedup, ',', efficiency
    CLOSE(unit_csv)
  END SUBROUTINE append_sweep_csv

END MODULE omp_timer
