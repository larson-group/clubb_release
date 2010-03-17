! $Id$
!----------------------------------------------------------------------
module cloud_rad_module

  implicit none

  public :: cloud_rad

  private ! Default Scope

  contains
!----------------------------------------------------------------------
  subroutine cloud_rad( rho, rcm, exner, Frad, radht,  & 
                        thlm_forcing )
! Description:
!   Subroutine to compute cloud IR radiation using a simple scheme
!   based on LWP

! References:
!   None

! Notes:
!  Based on GCSS ATEX intercomparison case
!----------------------------------------------------------------------

    use grid_class, only: gr ! Variable(s)

    use grid_class, only: ddzm ! Procedure(s)

    use constants, only: Cp ! Variable(s)

    use parameters_radiation, only: F0, kappa

    implicit none

    ! External
    intrinsic :: exp

    ! Input Variables
    real, dimension(gr%nnzp), intent(in) ::  & 
      rho,  & ! Density (thermo point)          [kg/m^3]
      rcm,  & ! Liquid water mixing ratio       [kg/kg]
      exner   ! Exner function                  [-]

    ! Output Variables
    real, dimension(gr%nnzp), intent(out) ::  & 
      Frad, & ! IR radiative flux               [W/m^2]
      radht   ! Radiative heating rate          [K/s]

    real, dimension(gr%nnzp), intent(inout) ::  & 
      thlm_forcing ! Radht + LS      [K/s]

    ! Local Variables
    real, dimension(1:gr%nnzp) :: LWP        ! Liquid Water Path

    integer :: k

    ! ---- Begin Code ----

    ! Compute liquid water path from top of the model
    ! We define liquid water path on momentum levels
    LWP(gr%nnzp) = 0.

    do k = gr%nnzp-1, 1, -1

      LWP(k) = LWP(k+1) + rho(k+1) * rcm(k+1) / gr%dzt(k+1)

    end do ! k=gr%nnzp..1

    ! Compute IR radiative flux

    do k = 1, gr%nnzp, 1

      Frad(k) = F0 * EXP( -kappa * 1.0 * LWP(k) )

    end do

    ! Compute IR heating rate

    radht(1:gr%nnzp) & 
    = ( -1.0/(Cp*rho(1:gr%nnzp) ) * ddzm( Frad(1:gr%nnzp) )  & 
      * 1.0 / exner(1:gr%nnzp) )

    radht(1)       = 0.
    radht(gr%nnzp) = 0.

    ! Note that for ATEX after 90 minutes, advect and clear air
    ! radiation must be added in from atex_tndcy
    thlm_forcing(1:gr%nnzp)  & 
    = thlm_forcing(1:gr%nnzp) + radht(1:gr%nnzp)

    return
  end subroutine cloud_rad

end module cloud_rad_module
