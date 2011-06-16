! $Id$
module physconst
!
! Dummy module for importing variables into morrison-gettelman microphysics
!---------------------------------------------------------------------------------------------------

  use constants_clubb, only: &
    gravit => grav,      &   ! acceleration due to gravity                      [m s-2]
    rair => Rd,          &   ! dry air gas constant for air                     [J kg-1 K-1]
    tmelt => T_freeze_K, &   ! temperature of melting point for water           [K]
    cpair => Cp,         &   ! specific heat at constant pressure for dry air   [J kg-1 K-1]
    rh2o => Rv,          &   ! gas constant for water vapor                     [J kg-1 K-1]
    rhoh2o => rho_lw,    &   ! Density of liquid water                          [kg m-3]
    latvap => Lv,        &   ! latent heat of vaporization                      [J kg-1]
    latice => Lf,        &   ! latent heat of fusion                            [J kg-1]
    epsilo => ep             ! Ratio of h2o to dry air molecular weights        [-]

  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none

  private

  public :: gravit, rair, tmelt, cpair, rh2o, r_universal, mwh2o, rhoh2o, latvap, latice, epsilo
   
  real(r8) :: &
    r_universal = 8314.472_r8, &
    mwh2o = 18.016_r8

end module physconst
