! ============================================================
! PROGRAM: test_energy
! PROJECT: SolarHPC
! PURPOSE: Energy conservation over 10 years.
!          Writes results/benchmarks/integrator_comparison.csv
! ============================================================
PROGRAM test_energy
  USE solar_constants
  USE nbody_force
  USE integrators
  IMPLICIT NONE

  INTEGER,  PARAMETER :: N    = N_BODIES
  REAL(dp), PARAMETER :: DT   = 1.0d0
  INTEGER,  PARAMETER :: NYRS = 10

  REAL(dp) :: pos0(3,N), vel0(3,N), mass(N)
  REAL(dp) :: pos_lf(3,N), vel_lf(3,N)
  REAL(dp) :: pos_rk(3,N), vel_rk(3,N)
  REAL(dp) :: E0, E_lf, E_rk, drift_lf, drift_rk
  INTEGER  :: s, yr, ios
  INTEGER  :: unit_ic, unit_csv
  CHARACTER(LEN=32) :: bname
  REAL(dp) :: mrd

  unit_ic = 10
  OPEN(UNIT=unit_ic, FILE='data/initial_conditions.dat', STATUS='OLD', &
       ACTION='READ', IOSTAT=ios)
  IF (ios /= 0) THEN
    WRITE(*,'(A)') '  [ERROR] data/initial_conditions.dat not found.'
    WRITE(*,'(A)') '  Run: python3 fetch_jpl.py'
    STOP 1
  END IF
  READ(unit_ic, '(A)') bname  ! header
  DO s = 1, N
    READ(unit_ic, *, IOSTAT=ios) bname, mrd, &
      pos0(1,s), pos0(2,s), pos0(3,s), vel0(1,s), vel0(2,s), vel0(3,s)
    IF (ios /= 0) THEN
      WRITE(*,'(A,I2)') '  [ERROR] Failed reading body ', s; STOP 1
    END IF
    mass(s) = mrd
  END DO
  CLOSE(unit_ic)

  pos_lf = pos0; vel_lf = vel0
  pos_rk = pos0; vel_rk = vel0

  CALL compute_energy(N, pos0, vel0, mass, E0)

  WRITE(*,'(A)') ''
  WRITE(*,'(A)') '  Energy conservation test (10 years, dt=1 day, N=10 bodies)'
  WRITE(*,'(A,ES14.6,A)') '  E0 = ', E0, ' J'
  WRITE(*,'(A)') ''
  WRITE(*,'(A)') '  Year | Leapfrog drift      | RK4 drift'
  WRITE(*,'(A)') '  -----+---------------------+---------------------'

  unit_csv = 20
  OPEN(UNIT=unit_csv, FILE='results/benchmarks/integrator_comparison.csv', &
       STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
  IF (ios == 0) WRITE(unit_csv,'(A)') &
    'sim_year,leapfrog_energy_drift,rk4_energy_drift,leapfrog_pos_error_AU,rk4_pos_error_AU'

  DO yr = 1, NYRS
    DO s = 1, 365
      CALL leapfrog_step(N, pos_lf, vel_lf, mass, DT)
      CALL rk4_step(N, pos_rk, vel_rk, mass, DT)
    END DO
    CALL compute_energy(N, pos_lf, vel_lf, mass, E_lf)
    CALL compute_energy(N, pos_rk, vel_rk, mass, E_rk)
    drift_lf = ABS(E_lf - E0) / ABS(E0)
    drift_rk = ABS(E_rk - E0) / ABS(E0)
    WRITE(*,'(A,I4,A,ES16.6,A,ES16.6)') '  ', yr, '   | ', drift_lf, '   | ', drift_rk
    IF (ios == 0) WRITE(unit_csv,'(I4,4(A,ES14.6))') yr,',',drift_lf,',',drift_rk,',0.0,0.0'
  END DO
  IF (ios == 0) CLOSE(unit_csv)

  WRITE(*,'(A)') ''
  IF (drift_lf < 1.0d-4) THEN
    WRITE(*,'(A,ES10.3)') '  [PASS] Leapfrog drift < 1e-4 after 10yr: ', drift_lf
  ELSE
    WRITE(*,'(A,ES10.3)') '  [WARN] Leapfrog drift: ', drift_lf
  END IF
  IF (drift_lf < drift_rk) THEN
    WRITE(*,'(A)') '  [PASS] Leapfrog conserves energy better than RK4'
  ELSE
    WRITE(*,'(A)') '  [WARN] Unexpected: RK4 drift <= leapfrog'
  END IF
  WRITE(*,'(A)') '  Written: results/benchmarks/integrator_comparison.csv'
  WRITE(*,'(A)') ''
END PROGRAM test_energy
