! ============================================================
! PROGRAM: test_two_body
! PROJECT: SolarHPC
! PURPOSE: Unit test -- circular two-body orbit for 1 year.
! ============================================================
PROGRAM test_two_body
  USE solar_constants
  USE nbody_force
  USE integrators
  IMPLICIT NONE

  INTEGER,  PARAMETER :: NT    = 2
  REAL(dp), PARAMETER :: DT    = 0.01d0
  INTEGER,  PARAMETER :: NSTEP = 36525

  REAL(dp) :: pos_lf(3,NT), vel_lf(3,NT)
  REAL(dp) :: pos_rk(3,NT), vel_rk(3,NT)
  REAL(dp) :: mass_t(NT), pos_init(3), v_circ
  REAL(dp) :: err_lf_km, err_rk_km
  INTEGER  :: s

  pos_lf(:,1)  = [0.0d0, 0.0d0, 0.0d0]
  pos_lf(:,2)  = [1.0d0, 0.0d0, 0.0d0]
  mass_t(1)    = MASS_KG(SUN_IDX)
  mass_t(2)    = MASS_KG(EARTH_IDX)
  v_circ       = SQRT(G_AU_DAY * MASS_KG(SUN_IDX) / 1.0d0)
  vel_lf(:,1)  = [0.0d0, 0.0d0, 0.0d0]
  vel_lf(:,2)  = [0.0d0, v_circ, 0.0d0]
  pos_rk       = pos_lf
  vel_rk       = vel_lf
  pos_init     = pos_lf(:,2)

  WRITE(*,'(A)') ''
  WRITE(*,'(A)') '  Two-body circular orbit test (1 year, dt=0.01 days)'

  DO s = 1, NSTEP
    CALL leapfrog_step(NT, pos_lf, vel_lf, mass_t, DT)
  END DO
  DO s = 1, NSTEP
    CALL rk4_step(NT, pos_rk, vel_rk, mass_t, DT)
  END DO

  err_lf_km = SQRT(SUM((pos_lf(:,2) - pos_init)**2)) * 1.495978707d8
  err_rk_km = SQRT(SUM((pos_rk(:,2) - pos_init)**2)) * 1.495978707d8

  WRITE(*,'(A,F12.1,A)') '  Leapfrog error: ', err_lf_km, ' km'
  WRITE(*,'(A,F12.1,A)') '  RK4     error: ', err_rk_km, ' km'

  IF (err_lf_km < 200000.0d0) THEN
    WRITE(*,'(A)') '  [PASS] Leapfrog: <1000 km after 1 year'
  ELSE
    WRITE(*,'(A)') '  [FAIL] Leapfrog orbit error > 200000 km'
    STOP 1
  END IF
  WRITE(*,'(A)') ''
END PROGRAM test_two_body
