!-----------------------------------------------------------------------
! $Id$

module constants

!       Description:
!       Contains frequently occuring model constants

!       References:
!       None
!-----------------------------------------------------------------------

use stats_precision, only:  & 
    time_precision ! Variable(s)

implicit none
 
public :: fstderr, fstdin, fstdout, pi_dp, pi, Cp, Lv, Ls, Lf, & 
          Rd, Rv, ep, ep1, ep2, kappa, grav, p0, vonk, rho_lw, & 
          sol_const, wtol, thltol, rttol, qttol, sstol, difftol, & 
          rc_tol, Nc_tol, rr_tol, Nr_tol, emin, eps,  & 
          sec_per_day, sec_per_hr, sec_per_min, g_per_kg,  & 
          Lscale_max, T_freeze_K
!    .            sclr_dim, hydromet_dim, sclrtol, 

private ! Default scope

! Fortran file unit I/O constants
integer, parameter ::  & 
fstderr = 0, fstdin = 5, fstdout = 6

! Mathematical Constants
double precision, parameter ::  & 
pi_dp = 3.14159265358979323846d0
real, parameter ::  & 
pi = 3.141592654 ! The ratio of radii to their circumference

! Physical constants
real, parameter ::  & 
Cp = 1004.67,  & ! Dry air specific heat at constant p [J/kg/K]
Lv = 2.5e6,    & ! Latent heat of vaporization         [J/kg]
Ls = 2.834e6,  & ! Latent heat of sublimation          [J/kg]
Lf = 3.33e5,   & ! Latent heat of fusion               [J/kg]
Rd = 287.04,   & ! Dry air gas constant                [J/kg/K]
Rv = 461.5    ! Water vapor gas constant            [J/kg/K]

real, parameter :: &
T_freeze_K = 273.15 ! Freezing point of water [K]


! Useful combinations of Rd and Rv
real, parameter ::  & 
ep  = Rd / Rv,     & ! ep  = 0.622  [-]
ep1 = (1.0-ep)/ep, & ! ep1 = 0.61   [-]
ep2 = 1.0/ep      ! ep2 = 1.61   [-]

real, parameter :: & 
kappa = Rd / Cp   ! kappa        [-]

! Changed g to grav to make it easier to find in the code 5/25/05
!       real, parameter :: grav  = 9.80665 ! Gravitational acceleration [m/s^2]
real, parameter ::  & 
grav = 9.81,    & ! Gravitational acceleration     [m/s^2]
p0   = 1.0e5   ! Reference pressure             [Pa]

! Von Karman's constant
! Constant of the logarithmic wind profile in the surface layer 
real, parameter ::  & 
vonk   = 0.4,     & ! Accepted value is 0.40 (+/-) 0.01      [-]
rho_lw = 1000.0  ! Density of liquid water                [kg/m^3]

! BUGSrad Constants
real, parameter ::  & 
sol_const = 1367.0 ! Solar constant     [W/m^2]

! Tolerances below which we consider moments to be zero
real, parameter ::  & 
wtol    = 2.e-2,  & ! [m/s]
thltol  = 1.e-2,  & ! [K]
rttol   = 1.e-5,  & ! [kg/kg]
qttol   = 1.e-4,  & ! [?]
sstol   = 1.e-8,  & ! [kg/kg]
difftol = 0.4    ! [?]


! Set tolerances for Khairoutdinov and Kogan rain microphysics
! to insure against numerical errors.
! The tolerance values for Nc, rr, and Nr insure against
! underflow errors in computing the PDF for l_kk_rain.  Basically,
! they insure that those values squared won't be less then 
! 10^-38, which is the lowest number that can be numerically
! represented.  However, the tolerance value for rc doubles as
! the lowest mixing ratio there can be to still officially have 
! a cloud at that level.  This is figured to be about 
! 1.0 x 10^-7 kg/kg.  Brian; February 10, 2007.
real, parameter ::  & 
rc_tol = 1.0E-7,   & ! [kg/kg]
Nc_tol = 1.0E-18,  & ! [#/m^3]
rr_tol = 1.0E-18,  & ! [kg/kg]
Nr_tol = 1.0E-18  ! [#/m^3]

! Minimum value for em (turbulence kinetic energy)
real, parameter :: emin = 1.0e-6  ! [m^2/s^2]

real, parameter ::  & 
eps = 1.0e-10 ! Small value to prevent a divide by zero

! Useful conversion factors.
real(kind=time_precision), parameter ::  & 
sec_per_day = 86400.0, & ! Seconds in a day.
sec_per_hr  = 3600.0,  & ! Seconds in an hour.
sec_per_min = 60.0    ! Seconds in a minute.

real, parameter :: & 
g_per_kg = 1000.0     ! Grams in a kilogram.

!       integer, parameter :: 
!    .  sclr_dim     = 2, ! Number of passive scalars. 
!    .  hydromet_dim = 5  ! Number of hydrometeor fields.

! Tolerance for new mixing scheme
! Currently it's elements are equal to rttol and thltol
!       real, dimension(sclr_dim), parameter :: 
!    .  sclrtol = (/thltol, rttol/)

! Maximum allowable value for Lscale [m].
! Value depends on whether the model is run by itself or as part
! of a host model.
real :: Lscale_max

end module constants
