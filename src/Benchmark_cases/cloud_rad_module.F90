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

    use constants_clubb, only: Cp ! Variable(s)

    use parameters_radiation, only: F0, kappa

    use clubb_precision, only: core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: exp

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) ::  & 
      rho,  & ! Density (thermo point)          [kg/m^3]
      rcm,  & ! Liquid water mixing ratio       [kg/kg]
      exner   ! Exner function                  [-]

    ! Output Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(out) ::  & 
      Frad, & ! IR radiative flux               [W/m^2]
      radht   ! Radiative heating rate          [K/s]

    real( kind = core_rknd ), dimension(gr%nz), intent(inout) ::  & 
      thlm_forcing ! Radht + LS      [K/s]

    ! Local Variables
    real( kind = core_rknd ), dimension(1:gr%nz) :: LWP        ! Liquid Water Path

    integer :: k

    ! ---- Begin Code ----

    ! Compute liquid water path from top of the model
    ! We define liquid water path on momentum levels
    LWP(gr%nz) = 0._core_rknd

    do k = gr%nz-1, 1, -1

      LWP(k) = LWP(k+1) + rho(k+1) * rcm(k+1) / gr%invrs_dzt(k+1)

    end do ! k=gr%nz..1

    ! Compute IR radiative flux

    do k = 1, gr%nz, 1

      Frad(k) = F0 * EXP( -kappa * 1.0_core_rknd * LWP(k) )

    end do

    ! Compute IR heating rate

    radht(1:gr%nz) & 
    = ( -1.0_core_rknd/(Cp*rho(1:gr%nz) ) * ddzm( Frad(1:gr%nz) )  & 
      * 1.0_core_rknd / exner(1:gr%nz) )

    radht(1)       = 0._core_rknd
    radht(gr%nz) = 0._core_rknd

    ! Note that for ATEX after 90 minutes, advect and clear air
    ! radiation must be added in from atex_tndcy
    thlm_forcing(1:gr%nz)  & 
    = thlm_forcing(1:gr%nz) + radht(1:gr%nz)

    return
  end subroutine cloud_rad

end module cloud_rad_module
