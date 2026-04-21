! ============================================================
! MODULE: coordinate_transform
! PROJECT: SolarHPC
! PURPOSE: All coordinate system conversions:
!          geographic <-> ECEF <-> ECI <-> barycentric,
!          ecliptic <-> equatorial, Julian Day <-> Gregorian.
!          Compiled FIFTH. Depends on: solar_constants.
! ============================================================
MODULE coordinate_transform
  USE solar_constants
  IMPLICIT NONE

  ! WGS-84 ellipsoid parameters
  REAL(dp), PARAMETER :: WGS84_A  = 6378137.0d0           ! semi-major axis (m)
  REAL(dp), PARAMETER :: WGS84_F  = 1.0d0 / 298.257223563d0
  REAL(dp), PARAMETER :: WGS84_B  = WGS84_A * (1.0d0 - WGS84_F)
  REAL(dp), PARAMETER :: WGS84_E2 = 1.0d0 - (WGS84_B / WGS84_A)**2

  ! Obliquity of ecliptic at J2000 (degrees)
  REAL(dp), PARAMETER :: OBLIQUITY_J2000 = 23.43929111d0

CONTAINS

  ! ----------------------------------------------------------
  ! geographic_to_ecef: lat/lon/alt -> ECEF Cartesian (metres)
  ! Uses WGS-84 ellipsoid.
  ! ----------------------------------------------------------
  SUBROUTINE geographic_to_ecef(lat_deg, lon_deg, alt_m, x_m, y_m, z_m)
    REAL(dp), INTENT(IN)  :: lat_deg, lon_deg, alt_m
    REAL(dp), INTENT(OUT) :: x_m, y_m, z_m

    REAL(dp) :: lat_r, lon_r, N_curve, sin_lat, cos_lat, sin_lon, cos_lon

    lat_r   = lat_deg * DEG_TO_RAD
    lon_r   = lon_deg * DEG_TO_RAD
    sin_lat = SIN(lat_r)
    cos_lat = COS(lat_r)
    sin_lon = SIN(lon_r)
    cos_lon = COS(lon_r)

    N_curve = WGS84_A / SQRT(1.0d0 - WGS84_E2 * sin_lat * sin_lat)

    x_m = (N_curve + alt_m)               * cos_lat * cos_lon
    y_m = (N_curve + alt_m)               * cos_lat * sin_lon
    z_m = (N_curve * (1.0d0 - WGS84_E2) + alt_m) * sin_lat
  END SUBROUTINE geographic_to_ecef

  ! ----------------------------------------------------------
  ! ecef_to_eci: rotate ECEF to ECI using Greenwich Sidereal Time
  ! ----------------------------------------------------------
  SUBROUTINE ecef_to_eci(x_ecef, y_ecef, z_ecef, jd, x_eci, y_eci, z_eci)
    REAL(dp), INTENT(IN)  :: x_ecef, y_ecef, z_ecef, jd
    REAL(dp), INTENT(OUT) :: x_eci, y_eci, z_eci

    REAL(dp) :: T_cent, gst_hours, gst_rad, cos_gst, sin_gst

    ! Julian centuries since J2000
    T_cent = (jd - J2000_JD) / 36525.0d0

    ! Greenwich Mean Sidereal Time (hours) -- Meeus Eq. 12.4
    gst_hours = 6.697375d0 + 2400.0513369d0 * T_cent &
              + 0.0000258d0 * T_cent * T_cent
    ! Reduce to [0, 24) hours
    gst_hours = MOD(gst_hours, 24.0d0)
    IF (gst_hours < 0.0d0) gst_hours = gst_hours + 24.0d0

    gst_rad = gst_hours * (PI / 12.0d0)   ! hours to radians
    cos_gst = COS(gst_rad)
    sin_gst = SIN(gst_rad)

    ! Rotation: ECI = R_z(-GST) * ECEF
    x_eci =  cos_gst * x_ecef - sin_gst * y_ecef
    y_eci =  sin_gst * x_ecef + cos_gst * y_ecef
    z_eci =  z_ecef
  END SUBROUTINE ecef_to_eci

  ! ----------------------------------------------------------
  ! eci_to_barycentric: add Earth's barycentric position (AU)
  ! to ECI position (converted AU).
  ! ----------------------------------------------------------
  SUBROUTINE eci_to_barycentric(x_eci_m, y_eci_m, z_eci_m, &
                                  earth_pos_AU, x_b, y_b, z_b)
    REAL(dp), INTENT(IN)  :: x_eci_m, y_eci_m, z_eci_m
    REAL(dp), INTENT(IN)  :: earth_pos_AU(3)
    REAL(dp), INTENT(OUT) :: x_b, y_b, z_b

    x_b = earth_pos_AU(1) + x_eci_m / AU_TO_M
    y_b = earth_pos_AU(2) + y_eci_m / AU_TO_M
    z_b = earth_pos_AU(3) + z_eci_m / AU_TO_M
  END SUBROUTINE eci_to_barycentric

  ! ----------------------------------------------------------
  ! barycentric_to_surface: given point near a body,
  ! compute surface lat/lon (degrees).
  ! ----------------------------------------------------------
  SUBROUTINE barycentric_to_surface(pos_bary, body_bary, body_idx, jd, &
                                      lat_deg, lon_deg)
    REAL(dp), INTENT(IN)  :: pos_bary(3), body_bary(3)
    INTEGER,  INTENT(IN)  :: body_idx
    REAL(dp), INTENT(IN)  :: jd
    REAL(dp), INTENT(OUT) :: lat_deg, lon_deg

    REAL(dp) :: r_body(3), r_mag, theta
    REAL(dp) :: sidereal_period

    ! Body-centred vector
    r_body = pos_bary - body_bary
    r_mag  = SQRT(DOT_PRODUCT(r_body, r_body))

    IF (r_mag < 1.0d-20) THEN
      lat_deg = 0.0d0; lon_deg = 0.0d0; RETURN
    END IF

    ! Sidereal rotation period (days)
    SELECT CASE (body_idx)
      CASE (EARTH_IDX)
        sidereal_period = 0.99726958d0
      CASE (MOON_IDX)
        sidereal_period = 27.321661d0
      CASE (MARS_IDX)
        sidereal_period = 1.02595676d0
      CASE DEFAULT
        sidereal_period = 1.0d0
    END SELECT

    ! Body rotation angle since J2000
    theta = TWO_PI * MOD((jd - J2000_JD) / sidereal_period, 1.0d0)

    ! Latitude from z-component
    lat_deg = ASIN(r_body(3) / r_mag) * RAD_TO_DEG

    ! Longitude: rotate by body rotation angle to body-fixed frame
    lon_deg = (ATAN2(r_body(2), r_body(1)) - theta) * RAD_TO_DEG
    ! Normalise to [-180, 180)
    DO WHILE (lon_deg >= 180.0d0)
      lon_deg = lon_deg - 360.0d0
    END DO
    DO WHILE (lon_deg < -180.0d0)
      lon_deg = lon_deg + 360.0d0
    END DO
  END SUBROUTINE barycentric_to_surface

  ! ----------------------------------------------------------
  ! ecliptic_to_equatorial: rotate by obliquity of ecliptic
  ! ----------------------------------------------------------
  SUBROUTINE ecliptic_to_equatorial(x_ecl, y_ecl, z_ecl, x_eq, y_eq, z_eq)
    REAL(dp), INTENT(IN)  :: x_ecl, y_ecl, z_ecl
    REAL(dp), INTENT(OUT) :: x_eq, y_eq, z_eq

    REAL(dp) :: eps, cos_eps, sin_eps

    eps     = OBLIQUITY_J2000 * DEG_TO_RAD
    cos_eps = COS(eps)
    sin_eps = SIN(eps)

    x_eq =  x_ecl
    y_eq =  cos_eps * y_ecl - sin_eps * z_ecl
    z_eq =  sin_eps * y_ecl + cos_eps * z_ecl
  END SUBROUTINE ecliptic_to_equatorial

  ! ----------------------------------------------------------
  ! jd_to_gregorian: Julian Day -> calendar date/time
  ! Algorithm: Meeus "Astronomical Algorithms" Chapter 7
  ! ----------------------------------------------------------
  SUBROUTINE jd_to_gregorian(jd, year, month, day, hour, minute, second)
    REAL(dp), INTENT(IN)  :: jd
    INTEGER,  INTENT(OUT) :: year, month, day, hour, minute
    REAL(dp), INTENT(OUT) :: second

    REAL(dp) :: jd_plus, frac, Z_r, F_r
    INTEGER  :: Z, A, B, C, D, E, alpha

    jd_plus = jd + 0.5d0
    Z_r     = FLOOR(jd_plus)
    F_r     = jd_plus - Z_r
    Z       = INT(Z_r)

    IF (Z < 2299161) THEN
      A = Z
    ELSE
      alpha = INT((REAL(Z, dp) - 1867216.25d0) / 36524.25d0)
      A = Z + 1 + alpha - alpha/4
    END IF

    B = A + 1524
    C = INT((REAL(B, dp) - 122.1d0) / 365.25d0)
    D = INT(365.25d0 * REAL(C, dp))
    E = INT(REAL(B - D, dp) / 30.6001d0)

    day = B - D - INT(30.6001d0 * REAL(E, dp))

    IF (E < 14) THEN
      month = E - 1
    ELSE
      month = E - 13
    END IF

    IF (month > 2) THEN
      year = C - 4716
    ELSE
      year = C - 4715
    END IF

    ! Time from fractional day
    frac   = F_r * 24.0d0
    hour   = INT(frac)
    frac   = (frac - REAL(hour, dp)) * 60.0d0
    minute = INT(frac)
    second = (frac - REAL(minute, dp)) * 60.0d0
  END SUBROUTINE jd_to_gregorian

  ! ----------------------------------------------------------
  ! gregorian_to_jd: calendar date -> Julian Day Number
  ! Meeus Chapter 7
  ! ----------------------------------------------------------
  SUBROUTINE gregorian_to_jd(year, month, day, hour, minute, second, jd)
    INTEGER,  INTENT(IN)  :: year, month, day, hour, minute
    REAL(dp), INTENT(IN)  :: second
    REAL(dp), INTENT(OUT) :: jd

    INTEGER  :: Y, M, A, B
    REAL(dp) :: D_frac

    Y = year
    M = month
    IF (M <= 2) THEN
      Y = Y - 1
      M = M + 12
    END IF

    A = Y / 100
    B = 2 - A + A/4

    D_frac = REAL(day, dp) &
           + REAL(hour, dp)/24.0d0 &
           + REAL(minute, dp)/1440.0d0 &
           + second / 86400.0d0

    jd = FLOOR(365.25d0 * REAL(Y + 4716, dp)) &
       + FLOOR(30.6001d0 * REAL(M + 1, dp))   &
       + D_frac + REAL(B, dp) - 1524.5d0
  END SUBROUTINE gregorian_to_jd

  ! ----------------------------------------------------------
  ! earth_rotation_velocity: surface velocity due to rotation
  ! Returns velocity vector in AU/day (ECI frame, equatorial)
  ! ----------------------------------------------------------
  SUBROUTINE earth_rotation_velocity(lat_deg, lon_deg, vx, vy, vz)
    REAL(dp), INTENT(IN)  :: lat_deg, lon_deg
    REAL(dp), INTENT(OUT) :: vx, vy, vz

    REAL(dp) :: omega_earth, v_mag, lat_r, lon_r

    ! Earth's sidereal rotation rate (AU/day)
    omega_earth = TWO_PI / 0.99726958d0   ! rad/day
    v_mag       = omega_earth * EARTH_RADIUS_AU * COS(lat_deg * DEG_TO_RAD)

    lat_r = lat_deg * DEG_TO_RAD
    lon_r = lon_deg * DEG_TO_RAD

    ! Prograde (eastward) unit vector in ECI: (-sin(lon), cos(lon), 0)
    vx = -v_mag * SIN(lon_r)
    vy =  v_mag * COS(lon_r)
    vz =  0.0d0
  END SUBROUTINE earth_rotation_velocity

  ! ----------------------------------------------------------
  ! format_gregorian: return a human-readable date string
  ! ----------------------------------------------------------
  SUBROUTINE format_gregorian(jd, date_str)
    REAL(dp),         INTENT(IN)  :: jd
    CHARACTER(LEN=20), INTENT(OUT) :: date_str

    INTEGER  :: yr, mo, dy, hr, mn
    REAL(dp) :: sc

    CALL jd_to_gregorian(jd, yr, mo, dy, hr, mn, sc)
    WRITE(date_str, '(I4,A,I2.2,A,I2.2,A,I2.2,A,I2.2)') &
      yr, '-', mo, '-', dy, ' ', hr, ':', mn
  END SUBROUTINE format_gregorian

END MODULE coordinate_transform
