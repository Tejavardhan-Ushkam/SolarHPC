! ============================================================
! PROGRAM: test_eclipse
! PROJECT: SolarHPC
! PURPOSE: Eclipse geometry unit tests + JD conversion.
! ============================================================
PROGRAM test_eclipse
  USE solar_constants
  USE coordinate_transform
  USE eclipse_geom
  IMPLICIT NONE

  REAL(dp) :: pos(3, N_BODIES)
  INTEGER  :: event_flag, i
  REAL(dp) :: umbra_frac, dur_min
  REAL(dp) :: jd_out, jd_err
  INTEGER  :: yr, mo, dy, hr, mn
  REAL(dp) :: sc

  DO i = 1, N_BODIES
    pos(:,i) = [REAL(i*3,dp), REAL(i*2,dp), 0.0d0]
  END DO

  pos(:, SUN_IDX)   = [0.0d0, 0.0d0, 0.0d0]
  pos(:, EARTH_IDX) = [1.0d0, 0.0d0, 0.0d0]
  pos(:, MOON_IDX)  = [1.0d0 - MOON_RADIUS_AU*80.0d0, 0.0d0, 0.0d0]

  WRITE(*,'(A)') ''
  WRITE(*,'(A)') '  Test 1: Aligned eclipse geometry'
  CALL check_solar_eclipse(pos, event_flag, umbra_frac, dur_min)
  WRITE(*,'(A,I1,A,F6.3)') '  event_flag=', event_flag, '  umbra=', umbra_frac
  IF (event_flag >= 1) THEN
    WRITE(*,'(A)') '  [PASS] Eclipse detected in aligned geometry'
  ELSE
    WRITE(*,'(A)') '  [FAIL] Eclipse not detected'
    STOP 1
  END IF

  WRITE(*,'(A)') '  Test 2: Off-axis (should be no eclipse)'
  pos(:, MOON_IDX) = [1.0d0, 0.1d0, 0.0d0]
  CALL check_solar_eclipse(pos, event_flag, umbra_frac, dur_min)
  IF (event_flag == 0) THEN
    WRITE(*,'(A)') '  [PASS] No false-positive eclipse'
  ELSE
    WRITE(*,'(A)') '  [FAIL] False positive eclipse'
    STOP 1
  END IF

  WRITE(*,'(A)') '  Test 3: JD round-trip (J2000.0 = 2451545.0)'
  CALL gregorian_to_jd(2000, 1, 1, 12, 0, 0.0d0, jd_out)
  jd_err = ABS(jd_out - J2000_JD)
  WRITE(*,'(A,F14.4,A,ES10.3)') '  JD = ', jd_out, '  err=', jd_err
  IF (jd_err < 1.0d-6) THEN
    WRITE(*,'(A)') '  [PASS] JD conversion'
  ELSE
    WRITE(*,'(A)') '  [FAIL] JD conversion error'
    STOP 1
  END IF
  WRITE(*,'(A)') ''
END PROGRAM test_eclipse
