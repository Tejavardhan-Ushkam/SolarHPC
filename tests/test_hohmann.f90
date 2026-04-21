! ============================================================
! PROGRAM: test_hohmann
! PROJECT: SolarHPC
! PURPOSE: Hohmann transfer dv vs Bate et al. reference value.
! ============================================================
PROGRAM test_hohmann
  USE solar_constants
  USE trajectory_prop
  IMPLICIT NONE

  REAL(dp) :: dv1, dv2, tf, dv_total, err_pct, tf_err
  REAL(dp) :: m_fuel

  CALL hohmann_transfer(1.0d0, 1.524d0, dv1, dv2, tf)
  dv_total = dv1 + dv2
  err_pct  = ABS(dv_total - 5591.0d0) / 5591.0d0 * 100.0d0
  tf_err   = ABS(tf       - 258.9d0)  / 258.9d0  * 100.0d0

  WRITE(*,'(A)') ''
  WRITE(*,'(A)') '  Hohmann Earth-Mars transfer test'
  WRITE(*,'(A,F8.1,A,F8.1,A)') '  dv_total = ', dv_total, ' m/s  (ref 5591 m/s)'
  WRITE(*,'(A,F7.1,A,F7.1,A)') '  t_trans  = ', tf, ' days  (ref 258.9 days)'
  WRITE(*,'(A,F6.3,A)') '  dv error: ', err_pct, ' %'

  IF (err_pct < 1.0d0 .AND. tf_err < 1.0d0) THEN
    WRITE(*,'(A)') '  [PASS] Hohmann dv and transfer time within 1% of reference'
  ELSE
    WRITE(*,'(A)') '  [FAIL] Hohmann dv or transfer time error > 1%'
    STOP 1
  END IF

  CALL tsiolkovsky(5590.0d0, ISP_CHEMICAL, 1000.0d0, m_fuel)
  WRITE(*,'(A,F8.1,A)') '  Tsiolkovsky fuel (dv=5590, Isp=450, payload=1000): ', m_fuel, ' kg'
  IF (ABS(m_fuel - 2549.17d0) / 2760.5d0 * 100.0d0 < 0.01d0) THEN
    WRITE(*,'(A)') '  [PASS] Tsiolkovsky rocket equation'
  ELSE
    WRITE(*,'(A)') '  [FAIL] Tsiolkovsky result out of tolerance'
    STOP 1
  END IF
  WRITE(*,'(A)') ''
END PROGRAM test_hohmann
