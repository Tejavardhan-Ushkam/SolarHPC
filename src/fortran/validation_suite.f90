! ============================================================
! MODULE: validation_suite
! PROJECT: SolarHPC
! PURPOSE: Formal test suite. Every test has PASS/FAIL outcome
!          with numeric evidence. Writes test_report.txt.
!          Compiled EIGHTH. Depends on all preceding modules.
! ============================================================
MODULE validation_suite
  USE solar_constants
  USE nbody_force
  USE integrators
  USE coordinate_transform
  USE trajectory_prop
  IMPLICIT NONE

CONTAINS

  ! ----------------------------------------------------------
  ! test_two_body_orbit: circular orbit for 1 year.
  ! Earth should return within <1000 km of starting position.
  ! ----------------------------------------------------------
  SUBROUTINE test_two_body_orbit(passed, pos_error_km)
    LOGICAL,  INTENT(OUT) :: passed
    REAL(dp), INTENT(OUT) :: pos_error_km

    INTEGER,  PARAMETER :: NT = 2
    REAL(dp), PARAMETER :: DT = 0.01d0
    INTEGER,  PARAMETER :: NSTEPS = 36525   ! 365.25 / 0.01

    REAL(dp) :: pos(3, NT), vel(3, NT), mass_t(NT)
    REAL(dp) :: pos_init(3), v_circ, pos_error_AU
    INTEGER  :: s

    pos(:,1) = [0.0d0, 0.0d0, 0.0d0]
    pos(:,2) = [1.0d0, 0.0d0, 0.0d0]
    mass_t(1) = MASS_KG(SUN_IDX)
    mass_t(2) = MASS_KG(EARTH_IDX)

    v_circ   = SQRT(G_AU_DAY * MASS_KG(SUN_IDX) / 1.0d0)
    vel(:,1) = [0.0d0, 0.0d0, 0.0d0]
    vel(:,2) = [0.0d0, v_circ, 0.0d0]

    pos_init = pos(:,2)

    DO s = 1, NSTEPS
      CALL leapfrog_step(NT, pos, vel, mass_t, DT)
    END DO

    pos_error_AU = SQRT(SUM((pos(:,2) - pos_init)**2))
    pos_error_km = pos_error_AU * 1.495978707d8

    passed = (pos_error_km < 1000.0d0)
    IF (passed) THEN
      WRITE(*,'(A,F10.1,A)') &
        '  [TEST 1] Two-body orbit: PASS  error=', pos_error_km, ' km'
    ELSE
      WRITE(*,'(A,F10.1,A)') &
        '  [TEST 1] Two-body orbit: FAIL  error=', pos_error_km, ' km'
    END IF
  END SUBROUTINE test_two_body_orbit

  ! ----------------------------------------------------------
  ! test_hohmann_delta_v: Earth-Mars Hohmann.
  ! Reference: Bate, Mueller, White (1971): ~5591 m/s total.
  ! ----------------------------------------------------------
  SUBROUTINE test_hohmann_delta_v(passed, dv_error_pct)
    LOGICAL,  INTENT(OUT) :: passed
    REAL(dp), INTENT(OUT) :: dv_error_pct

    REAL(dp), PARAMETER :: R1_AU  = 1.0d0
    REAL(dp), PARAMETER :: R2_AU  = 1.524d0
    REAL(dp), PARAMETER :: DV_REF = 5591.0d0   ! m/s  (Bate et al.)
    REAL(dp), PARAMETER :: TF_REF = 258.9d0    ! days

    REAL(dp) :: dv1, dv2, tf, dv_total, tf_err_pct

    CALL hohmann_transfer(R1_AU, R2_AU, dv1, dv2, tf)
    dv_total    = dv1 + dv2
    dv_error_pct= ABS(dv_total - DV_REF) / DV_REF * 100.0d0
    tf_err_pct  = ABS(tf       - TF_REF) / TF_REF * 100.0d0

    passed = (dv_error_pct < 1.0d0 .AND. tf_err_pct < 1.0d0)
    IF (passed) THEN
      WRITE(*,'(A,F8.1,A,F6.3,A)') &
        '  [TEST 2] Hohmann Earth-Mars: PASS  dv=', dv_total, ' m/s  err=', dv_error_pct, '%'
    ELSE
      WRITE(*,'(A,F8.1,A,F6.3,A)') &
        '  [TEST 2] Hohmann Earth-Mars: FAIL  dv=', dv_total, ' m/s  err=', dv_error_pct, '%'
    END IF
  END SUBROUTINE test_hohmann_delta_v

  ! ----------------------------------------------------------
  ! test_tsiolkovsky: known analytical answer.
  ! dv=5590 m/s, Isp=450s, payload=1000kg -> fuel=2760.5 kg
  ! ----------------------------------------------------------
  SUBROUTINE test_tsiolkovsky_eqn(passed, error_pct)
    LOGICAL,  INTENT(OUT) :: passed
    REAL(dp), INTENT(OUT) :: error_pct

    REAL(dp), PARAMETER :: DV_IN   = 5590.0d0
    REAL(dp), PARAMETER :: ISP_IN  = 450.0d0
    REAL(dp), PARAMETER :: M_PAY   = 1000.0d0
    REAL(dp), PARAMETER :: M_REF   = 2760.5d0  ! kg  (analytical)

    REAL(dp) :: m_fuel

    CALL tsiolkovsky(DV_IN, ISP_IN, M_PAY, m_fuel)
    error_pct = ABS(m_fuel - M_REF) / M_REF * 100.0d0
    passed    = (error_pct < 0.01d0)

    IF (passed) THEN
      WRITE(*,'(A,F8.1,A,F8.4,A)') &
        '  [TEST 3] Tsiolkovsky: PASS  fuel=', m_fuel, ' kg  err=', error_pct, '%'
    ELSE
      WRITE(*,'(A,F8.1,A,F8.4,A)') &
        '  [TEST 3] Tsiolkovsky: FAIL  fuel=', m_fuel, ' kg  err=', error_pct, '%'
    END IF
  END SUBROUTINE test_tsiolkovsky_eqn

  ! ----------------------------------------------------------
  ! test_energy_conservation: 1-year leapfrog, 10-body system.
  ! Drift < 0.01 % is a pass.
  ! ----------------------------------------------------------
  SUBROUTINE test_energy_conservation(pos0, vel0, mass, passed, drift_pct)
    REAL(dp), INTENT(IN)  :: pos0(3, N_BODIES), vel0(3, N_BODIES), mass(N_BODIES)
    LOGICAL,  INTENT(OUT) :: passed
    REAL(dp), INTENT(OUT) :: drift_pct

    INTEGER,  PARAMETER :: NSTEPS = 365
    REAL(dp), PARAMETER :: DT     = 1.0d0

    REAL(dp) :: pos(3, N_BODIES), vel(3, N_BODIES)
    REAL(dp) :: E0, E1
    INTEGER  :: s

    pos = pos0
    vel = vel0

    CALL compute_energy(N_BODIES, pos, vel, mass, E0)
    DO s = 1, NSTEPS
      CALL leapfrog_step(N_BODIES, pos, vel, mass, DT)
    END DO
    CALL compute_energy(N_BODIES, pos, vel, mass, E1)

    drift_pct = ABS(E1 - E0) / ABS(E0) * 100.0d0
    passed    = (drift_pct < 0.01d0)

    IF (passed) THEN
      WRITE(*,'(A,ES10.3,A)') &
        '  [TEST 4] Energy conservation (1yr): PASS  drift=', drift_pct, '%'
    ELSE
      WRITE(*,'(A,ES10.3,A)') &
        '  [TEST 4] Energy conservation (1yr): FAIL  drift=', drift_pct, '%'
    END IF
  END SUBROUTINE test_energy_conservation

  ! ----------------------------------------------------------
  ! test_jd_conversion: round-trip JD <-> Gregorian.
  ! J2000.0 = 2451545.0 exactly.
  ! ----------------------------------------------------------
  SUBROUTINE test_jd_conversion(passed, jd_error)
    LOGICAL,  INTENT(OUT) :: passed
    REAL(dp), INTENT(OUT) :: jd_error

    REAL(dp) :: jd_out
    CALL gregorian_to_jd(2000, 1, 1, 12, 0, 0.0d0, jd_out)
    jd_error = ABS(jd_out - J2000_JD)
    passed   = (jd_error < 1.0d-6)

    IF (passed) THEN
      WRITE(*,'(A,ES10.3)') '  [TEST 5] JD conversion: PASS  error=', jd_error
    ELSE
      WRITE(*,'(A,ES10.3)') '  [TEST 5] JD conversion: FAIL  error=', jd_error
    END IF
  END SUBROUTINE test_jd_conversion

  ! ----------------------------------------------------------
  ! test_force_unit: Sun + Earth at 1 AU, check force magnitude.
  ! ----------------------------------------------------------
  SUBROUTINE test_force_unit(passed)
    LOGICAL, INTENT(OUT) :: passed
    CALL verify_forces(passed)
    IF (.NOT. passed) WRITE(*,*) '  [TEST 6] Force unit: FAIL'
  END SUBROUTINE test_force_unit

  ! ----------------------------------------------------------
  ! run_all_tests: master runner. Writes test_report.txt.
  ! ----------------------------------------------------------
  SUBROUTINE run_all_tests(pos0, vel0, mass, all_passed, n_passed, n_total)
    REAL(dp), INTENT(IN)  :: pos0(3, N_BODIES), vel0(3, N_BODIES), mass(N_BODIES)
    LOGICAL,  INTENT(OUT) :: all_passed
    INTEGER,  INTENT(OUT) :: n_passed, n_total

    LOGICAL  :: p1, p2, p3, p4, p5, p6
    REAL(dp) :: v1, v2, v3, v4, v5

    INTEGER  :: unit_rpt, ios
    CHARACTER(LEN=20) :: result_str

    n_total  = 6
    n_passed = 0

    WRITE(*,*)
    WRITE(*,'(A)') '  ========================================'
    WRITE(*,'(A)') '  SolarHPC Validation Test Suite'
    WRITE(*,'(A)') '  ========================================'

    CALL test_force_unit(p6)
    CALL test_two_body_orbit(p1, v1)
    CALL test_hohmann_delta_v(p2, v2)
    CALL test_tsiolkovsky_eqn(p3, v3)
    CALL test_energy_conservation(pos0, vel0, mass, p4, v4)
    CALL test_jd_conversion(p5, v5)

    IF (p1) n_passed = n_passed + 1
    IF (p2) n_passed = n_passed + 1
    IF (p3) n_passed = n_passed + 1
    IF (p4) n_passed = n_passed + 1
    IF (p5) n_passed = n_passed + 1
    IF (p6) n_passed = n_passed + 1

    all_passed = (n_passed == n_total)

    WRITE(*,'(A)') '  ========================================'
    WRITE(*,'(A,I1,A,I1,A)') '  Results: ', n_passed, '/', n_total, ' passed'
    IF (all_passed) THEN
      WRITE(*,'(A)') '  STATUS: ALL TESTS PASSED'
    ELSE
      WRITE(*,'(A)') '  STATUS: SOME TESTS FAILED -- check physics constants'
    END IF
    WRITE(*,'(A)') '  ========================================'
    WRITE(*,*)

    ! Write report file
    unit_rpt = 99
    OPEN(UNIT=unit_rpt, FILE='results/validation/test_report.txt', &
         STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
    IF (ios /= 0) RETURN

    WRITE(unit_rpt,'(A)') '=== SolarHPC Validation Test Report ==='
    WRITE(unit_rpt,'(A,I1,A,I1)') 'Tests passed: ', n_passed, ' / ', n_total
    WRITE(unit_rpt,*)
    WRITE(unit_rpt,'(A,A,F12.1,A)') 'TEST 1 Two-body orbit (1yr): ', &
      MERGE('PASS','FAIL',p1), v1, ' km position error'
    WRITE(unit_rpt,'(A,A,F8.3,A)') 'TEST 2 Hohmann Earth-Mars:  ', &
      MERGE('PASS','FAIL',p2), v2, ' % dv error'
    WRITE(unit_rpt,'(A,A,F10.4,A)') 'TEST 3 Tsiolkovsky:         ', &
      MERGE('PASS','FAIL',p3), v3, ' % fuel error'
    WRITE(unit_rpt,'(A,A,ES10.3,A)') 'TEST 4 Energy conservation: ', &
      MERGE('PASS','FAIL',p4), v4, ' % drift in 1yr'
    WRITE(unit_rpt,'(A,A,ES10.3,A)') 'TEST 5 JD conversion:       ', &
      MERGE('PASS','FAIL',p5), v5, ' days error'
    WRITE(unit_rpt,'(A,A)') 'TEST 6 Force unit:          ', &
      MERGE('PASS','FAIL',p6)
    CLOSE(unit_rpt)
    WRITE(*,'(A)') '  [OK] Test report written to results/validation/test_report.txt'
  END SUBROUTINE run_all_tests

END MODULE validation_suite
