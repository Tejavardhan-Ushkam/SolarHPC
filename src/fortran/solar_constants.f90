! ============================================================
! MODULE: solar_constants
! PROJECT: SolarHPC -- Parallel N-body Solar System Simulator
! PURPOSE: Single source of truth for all physical constants,
!          body indices, masses, and array dimensions.
!          Every other Fortran module USE-s this one.
!          Compiled FIRST -- no dependencies.
! ============================================================
MODULE solar_constants
  IMPLICIT NONE

  ! ---- Precision kind ----
  INTEGER, PARAMETER :: dp = KIND(1.0d0)   ! double precision

  ! ---- Fundamental physical constants ----
  REAL(dp), PARAMETER :: G_SI         = 6.674d-11      ! m^3 kg^-1 s^-2
  REAL(dp), PARAMETER :: AU_TO_M      = 1.495978707d11 ! metres per AU
  REAL(dp), PARAMETER :: DAY_TO_SEC   = 86400.0d0       ! seconds per day
  REAL(dp), PARAMETER :: SPEED_LIGHT  = 2.99792458d8    ! m/s

  ! G in simulation units: AU^3 kg^-1 day^-2
  REAL(dp), PARAMETER :: G_AU_DAY = G_SI * DAY_TO_SEC**2 / AU_TO_M**3

  ! Conversion factors
  REAL(dp), PARAMETER :: AU_PER_DAY_TO_M_PER_S = AU_TO_M / DAY_TO_SEC
  REAL(dp), PARAMETER :: M_PER_S_TO_AU_PER_DAY = DAY_TO_SEC / AU_TO_M

  ! ---- Mathematical constants ----
  REAL(dp), PARAMETER :: PI         = 3.14159265358979323846d0
  REAL(dp), PARAMETER :: TWO_PI     = 2.0d0 * PI
  REAL(dp), PARAMETER :: DEG_TO_RAD = PI / 180.0d0
  REAL(dp), PARAMETER :: RAD_TO_DEG = 180.0d0 / PI

  ! ---- Epoch ----
  REAL(dp), PARAMETER :: J2000_JD = 2451545.0d0   ! Julian Day of J2000.0

  ! ---- Simulation array dimensions ----
  INTEGER, PARAMETER :: N_BODIES  = 10    ! Sun + 9 bodies
  INTEGER, PARAMETER :: MAX_BODIES = 20   ! headroom for spacecraft

  ! ---- Body indices (1-based, Fortran convention) ----
  INTEGER, PARAMETER :: SUN_IDX     = 1
  INTEGER, PARAMETER :: MERCURY_IDX = 2
  INTEGER, PARAMETER :: VENUS_IDX   = 3
  INTEGER, PARAMETER :: EARTH_IDX   = 4
  INTEGER, PARAMETER :: MOON_IDX    = 5
  INTEGER, PARAMETER :: MARS_IDX    = 6
  INTEGER, PARAMETER :: JUPITER_IDX = 7
  INTEGER, PARAMETER :: SATURN_IDX  = 8
  INTEGER, PARAMETER :: URANUS_IDX  = 9
  INTEGER, PARAMETER :: NEPTUNE_IDX = 10

  ! ---- Body names (for output formatting) ----
  CHARACTER(LEN=10), PARAMETER :: BODY_NAMES(N_BODIES) = [ &
    'Sun       ', 'Mercury   ', 'Venus     ', 'Earth     ', 'Moon      ', &
    'Mars      ', 'Jupiter   ', 'Saturn    ', 'Uranus    ', 'Neptune   '  ]

  ! ---- Body masses (kg) ----
  REAL(dp), PARAMETER :: MASS_KG(N_BODIES) = [ &
    1.989000d+30, &   ! Sun
    3.285000d+23, &   ! Mercury
    4.867000d+24, &   ! Venus
    5.972000d+24, &   ! Earth
    7.342000d+22, &   ! Moon
    6.417000d+23, &   ! Mars
    1.898000d+27, &   ! Jupiter
    5.683000d+26, &   ! Saturn
    8.681000d+25, &   ! Uranus
    1.024000d+26  ]   ! Neptune

  ! ---- Body radii (AU) ----
  REAL(dp), PARAMETER :: BODY_RADIUS_AU(N_BODIES) = [ &
    4.65047d-3,  &   ! Sun
    1.63083d-5,  &   ! Mercury
    4.04540d-5,  &   ! Venus
    4.26352d-5,  &   ! Earth
    1.16139d-5,  &   ! Moon
    2.26874d-5,  &   ! Mars
    4.77895d-4,  &   ! Jupiter
    4.02866d-4,  &   ! Saturn
    1.70738d-4,  &   ! Uranus
    1.65507d-4   ]   ! Neptune

  ! Convenient named aliases for key radii
  REAL(dp), PARAMETER :: SUN_RADIUS_AU   = 4.65047d-3
  REAL(dp), PARAMETER :: EARTH_RADIUS_AU = 4.26352d-5
  REAL(dp), PARAMETER :: MOON_RADIUS_AU  = 1.16139d-5
  REAL(dp), PARAMETER :: MARS_RADIUS_AU  = 2.26874d-5

  ! ---- Numerical parameters ----
  REAL(dp), PARAMETER :: SOFTENING_EPS = 1.0d-10  ! AU, avoids singularity

  ! ---- Spacecraft / propulsion ----
  REAL(dp), PARAMETER :: ISP_CHEMICAL = 450.0d0   ! seconds, specific impulse
  REAL(dp), PARAMETER :: G0           = 9.80665d0 ! m/s^2 standard gravity

  ! ---- Eclipse detection ----
  REAL(dp), PARAMETER :: MIN_ECLIPSE_SEPARATION_DAYS = 20.0d0
  REAL(dp), PARAMETER :: SAROS_PERIOD_DAYS           = 6585.3211d0

  ! ---- Thread sweep parameters ----
  INTEGER, PARAMETER :: N_THREAD_COUNTS                = 5
  INTEGER, PARAMETER :: THREAD_COUNTS(N_THREAD_COUNTS) = [1, 2, 4, 8, 16]
  INTEGER, PARAMETER :: N_SWEEP_REPEATS                = 5

  ! ---- Output decimation (write every N steps) ----
  INTEGER, PARAMETER :: DEFAULT_DECIMATE = 10

END MODULE solar_constants
