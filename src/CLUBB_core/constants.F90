!-----------------------------------------------------------------------------
! $Id$
!=============================================================================
module constants

  ! Description:
  ! Contains frequently occuring model constants

  ! References:
  ! None
  !---------------------------------------------------------------------------

  use stats_precision, only:  & 
      time_precision ! Variable(s)

  implicit none
 
  public :: fstderr, fstdin, fstdout, var_length, &
            pi_dp, pi, sqrt_2pi, sqrt_2, three_halves, &
            Cp, Lv, Ls, Lf, Rd, Rv, ep, ep1, ep2, &
            kappa, grav, p0, vonk, rho_lw, &
            wtol, thltol, rttol, s_mellor_tol, thltol_mfl, rttol_mfl, & 
            wtol_sqd, rc_tol, Nc_tol, rr_tol, Nr_tol, emin, &
            eps, zero_threshold, max_mag_correlation, sec_per_day, &
            sec_per_hr, sec_per_min, g_per_kg, T_freeze_K, &
            Skw_max_mag, Skw_max_mag_sqd, stefan_boltzmann, &
            cm3_per_m3, pascal_per_mb,  &
            gamma_over_implicit_ts

  private ! Default scope

  ! Fortran file unit I/O constants
  integer, parameter ::  & 
    fstderr = 0, fstdin = 5, fstdout = 6, var_length = 30

  ! Mathematical Constants
  double precision, parameter ::  & 
    pi_dp = 3.14159265358979323846d0

  real, parameter ::  & 
    pi = 3.141592654 ! The ratio of radii to their circumference

  real, parameter :: &
    sqrt_2pi = 2.5066282746310005024, &  ! sqrt(2*pi)
    sqrt_2   = 1.4142135623730950488     ! sqrt(2)

  real, parameter :: &
    three_halves = 3.0/2.0   ! 3/2

  ! Physical constants
  real, parameter ::  & 
    Cp = 1004.67,  & ! Dry air specific heat at constant p [J/kg/K]
    Lv = 2.5e6,    & ! Latent heat of vaporization         [J/kg]
    Ls = 2.834e6,  & ! Latent heat of sublimation          [J/kg]
    Lf = 3.33e5,   & ! Latent heat of fusion               [J/kg]
    Rd = 287.04,   & ! Dry air gas constant                [J/kg/K]
    Rv = 461.5       ! Water vapor gas constant            [J/kg/K]


  real, parameter :: &
    stefan_boltzmann = 5.6704e-8 ! Stefan-Boltzmann constant [W/(m^2 K^4)]

  real, parameter :: &
    T_freeze_K = 273.15 ! Freezing point of water [K]

  ! Useful combinations of Rd and Rv
  real, parameter ::  & 
    ep  = Rd / Rv,    & ! ep  = 0.622  [-]
    ep1 = (1.0-ep)/ep,& ! ep1 = 0.61   [-]
    ep2 = 1.0/ep        ! ep2 = 1.61   [-]

  real, parameter :: & 
    kappa = Rd / Cp     ! kappa        [-]

  ! Changed g to grav to make it easier to find in the code 5/25/05
  ! real, parameter :: grav  = 9.80665 ! Gravitational acceleration [m/s^2]
  real, parameter :: & 
    grav = 9.81, & ! Gravitational acceleration     [m/s^2]
    p0   = 1.0e5   ! Reference pressure             [Pa]

  ! Von Karman's constant
  ! Constant of the logarithmic wind profile in the surface layer 
  real, parameter :: & 
    vonk   = 0.4,    & ! Accepted value is 0.40 (+/-) 0.01      [-]
    rho_lw = 1000.0    ! Density of liquid water                [kg/m^3]

  ! Tolerances below which we consider moments to be zero
  real, parameter ::  & 
    wtol    = 2.e-2,  & ! [m/s]
    thltol  = 1.e-2,  & ! [K]
    rttol   = 1.e-8,  & ! [kg/kg]
    s_mellor_tol   = 1.e-8     ! [kg/kg]

  ! Tolerances for use by the monatonic flux limiter.
  ! rttol_mfl is larger than rttol. rttol is extremely small
  ! (1e-8) to prevent spurious cloud formation aloft in LBA.
  ! rttol_mfl is larger (1e-4) to prevent the mfl from
  ! depositing moisture at the top of the domain.
  real, parameter :: &
    thltol_mfl = 1.e-2, & ! [K]
    rttol_mfl = 1.e-4     ! [kg/kg]
      
  ! The tolerance for w'^2 is the square of the tolerance for w.
  real, parameter :: &
    wtol_sqd = wtol**2 ! [m^2/s^2]

  real, parameter :: &
    Skw_max_mag = 4.5  ! Max magnitude of skewness     [-]

  real, parameter :: &
    Skw_max_mag_sqd = Skw_max_mag**2 ! Max mag. of Skw squared [-]

  ! Set tolerances for Khairoutdinov and Kogan rain microphysics to insure
  ! against numerical errors.  The tolerance values for Nc, rr, and Nr insure
  ! against underflow errors in computing the PDF for l_kk_rain.  Basically,
  ! they insure that those values squared won't be less then 10^-38, which is
  ! the lowest number that can be numerically represented.  However, the
  ! tolerance value for rc doubles as the lowest mixing ratio there can be to
  ! still officially have a cloud at that level.  This is figured to be about 
  ! 1.0 x 10^-7 kg/kg.  Brian; February 10, 2007.
  real, parameter :: & 
    rc_tol = 1.0E-7,  & ! [kg/kg]
    Nc_tol = 1.0E-10, & ! [#/kg]
    rr_tol = 1.0E-10, & ! [kg/kg]
    Nr_tol = 1.0E-10    ! [#/kg]

  ! Minimum value for em (turbulence kinetic energy)
  ! If anisotropic TKE is enabled, em = (1/2) * ( up2 + vp2 + wp2 );
  ! otherwise, em = (3/2) * wp2.  Since up2, vp2, and wp2 all have
  ! the same minimum threshold value of wtol_sqd, em cannot be less
  ! than (3/2) * wtol_sqd.  Thus, emin = (3/2) * wtol_sqd.
  real, parameter :: emin = 1.5 * wtol_sqd  ! [m^2/s^2]

  real, parameter ::  & 
    eps = 1.0e-10 ! Small value to prevent a divide by zero

  real, parameter ::  &
    zero_threshold = 0.0 ! Defining a threshold to be 0.

  ! The maximum absolute value (or magnitude) that a correlation is allowed to
  ! have.  Statistically, a correlation is not allowed to be less than -1 or 
  ! greater than 1, so the maximum magnitude would be 1.
  real, parameter :: &
    max_mag_correlation = 0.99

  ! "Over-implicit" weighted time step.
  !
  ! The weight of the implicit portion of a term is controlled by the factor
  ! gamma_over_implicit_ts (abbreviated "gamma" in the expression below).  A
  ! factor is added to the right-hand side of the equation in order to balance a
  ! weight that is not equal to 1, such that:
  !
  !      -y(t) * [ gamma * X(t+1) + ( 1 - gamma ) * X(t) ] + RHS;
  !
  ! where X is the variable that is being solved for in a predictive equation
  ! (such as w'^3, w'th_l', r_t'^2, etc), y(t) is the linearized portion of the
  ! term that gets treated implicitly, and RHS is the portion of the term that
  ! is always treated explicitly.  A weight of greater than 1 can be applied to
  ! make the term more numerically stable.
  !
  !    gamma_over_implicit_ts          Effect on term
  !
  !            0.0               Term becomes completely explicit
  !
  !            1.0               Standard implicit portion of the term;
  !                              as it was without the weighting factor.
  !
  !            1.5               Strongly weighted implicit portion of the term;
  !                              increased numerical stability.
  !
  !            2.0               More strongly weighted implicit portion of the
  !                              term; increased numerical stability.
  !
  ! Note:  The "over-implicit" weighted time step is only applied to terms that
  !        tend to significantly decrease the amount of numerical stability for
  !        variable X.
  !        The "over-implicit" weighted time step is applied to the turbulent
  !        advection term for the following variables:
  !           w'^3 (also applied to the turbulent production term), found in
  !           module advance_wp2_wp3_module;
  !           w'r_t', w'th_l', and w'sclr', found in
  !           module advance_xm_wpxp_module; and
  !           r_t'^2, th_l'^2, r_t'th_l', u'^2, v'^2, sclr'^2, sclr'r_t',
  !           and sclr'th_l', found in module advance_xp2_xpyp_module.
  real, parameter :: &
    gamma_over_implicit_ts = 1.50

  ! Useful conversion factors.
  real(kind=time_precision), parameter ::  & 
    sec_per_day = 86400.0, & ! Seconds in a day.
    sec_per_hr  = 3600.0,  & ! Seconds in an hour.
    sec_per_min = 60.0       ! Seconds in a minute.

  real, parameter :: & 
    g_per_kg = 1000.0     ! Grams in a kilogram.

  real, parameter :: &
    pascal_per_mb = 100.0 ! Pascals per Millibar

  real, parameter :: & 
    cm3_per_m3 = 1.e6 ! Cubic centimeters per cubic meter

!=============================================================================

end module constants
