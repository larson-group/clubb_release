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
    epsil => ep              ! Ratio of h2o to dry air molecular weights        [-]

  implicit none

  private

  public :: gravit, rair, tmelt, cpair, rh2o, r_universal, mwh2o, rhoh2o, latvap, latice, epsil
    
  ! These variables are not used anywhere in MG, they are just imported. Because of this we
  ! are setting them to dummy values
  integer :: &
    r_universal = 0, &
    mwh2o = 0

end module physconst
