!----------------------------------------------------------------------
! $Id$
module lba

  ! Description:
  !   Contains subroutines for the LBA case.
  ! References:
  !   http://www.mmm.ucar.edu/gcss-wg4/gcss/case4.html
  !----------------------------------------------------------------------

  implicit none

  private ! Default Scope

  public :: lba_tndcy, lba_sfclyr

  contains

  !----------------------------------------------------------------------
  subroutine lba_tndcy( ngrdcol, sclr_dim, edsclr_dim, sclr_idx, &
                        gr, thlm_forcing, rtm_forcing, & 
                        sclrm_forcing, edsclrm_forcing )
    !       Description:
    !       Subroutine to set theta and water tendencies for LBA case.

    !       References:
    !       http://www.mmm.ucar.edu/gcss-wg4/gcss/case4.html
    !----------------------------------------------------------------------

    use array_index, only: &
      sclr_idx_type

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    use grid_class, only: &
      grid

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
      ngrdcol, &
      sclr_dim, & 
      edsclr_dim

    type (sclr_idx_type), intent(in) :: &
      sclr_idx

    type(grid), intent(in) :: &
      gr

    !--------------------- Output Variables ---------------------
    real( kind = core_rknd ), intent(out), dimension(ngrdcol,gr%nzt) :: & 
      thlm_forcing, & ! Liquid water potential temperature tendency  [K/s]
      rtm_forcing     ! Total water mixing ratio tendency            [kg/kg/s]

    real( kind = core_rknd ), intent(out), dimension(ngrdcol,gr%nzt,sclr_dim) :: & 
      sclrm_forcing ! Passive scalar forcing [units vary]

    real( kind = core_rknd ), intent(out), dimension(ngrdcol,gr%nzt,edsclr_dim) :: & 
      edsclrm_forcing ! Passive eddy-scalar forcing [units vary]

    integer :: i, k

    !--------------------- Begin Code ---------------------

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, gr%nzt
      do i = 1, ngrdcol

        ! Large-scale temperature tendency
        thlm_forcing(i,k) = 0.0_core_rknd

        ! Large-scale advective moisture tendency
        rtm_forcing(i,k) = 0.0_core_rknd

      end do
    end do

    if ( sclr_dim > 0 ) then
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, gr%nzt
        do i = 1, ngrdcol
          ! Test scalars with thetal and rt if desired
          if ( sclr_idx%iisclr_thl > 0 ) sclrm_forcing(i,k,sclr_idx%iisclr_thl) = thlm_forcing(i,k)
          if ( sclr_idx%iisclr_rt  > 0 ) sclrm_forcing(i,k,sclr_idx%iisclr_rt)  = rtm_forcing(i,k)
        end do
      end do
    end if

    if ( edsclr_dim > 0 ) then
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, gr%nzt
        do i = 1, ngrdcol
          if ( sclr_idx%iiedsclr_thl > 0 ) edsclrm_forcing(i,k,sclr_idx%iiedsclr_thl) = thlm_forcing(i,k)
          if ( sclr_idx%iiedsclr_rt  > 0 ) edsclrm_forcing(i,k,sclr_idx%iiedsclr_rt)  = rtm_forcing(i,k)
        end do
      end do
    end if

    return
  end subroutine lba_tndcy

  !----------------------------------------------------------------------
  subroutine lba_sfclyr( ngrdcol, time_current, time_initial, & 
                         z, rho_sfc, thlm_sfc, ubar,  & 
                         wpthlp_sfc, wprtp_sfc, ustar )

    !       Description:
    !       This subroutine computes surface fluxes of horizontal momentum,
    !       heat and moisture according to GCSS BOMEX specifications

    !       References:
    !       Grabowski, et al. (2005)
    !       http://www.mmm.ucar.edu/gcss-wg4/gcss/case4.html
    !----------------------------------------------------------------------

    use constants_clubb, only: pi, grav, sec_per_hr ! Variable(s)

    use clubb_precision, only: time_precision, core_rknd ! Variable(s)

    use diag_ustar_module, only: diag_ustar ! Variable(s)

    use sfc_flux, only: convert_sens_ht_to_km_s, convert_latent_ht_to_m_s ! Procedure(s)

    implicit none

    intrinsic :: max, sqrt

    ! Constant Parameters
    real( kind = core_rknd ), parameter ::  & 
      z0    = 0.035_core_rknd  ! ARM mom. roughness height

    ! Input Variables
    integer :: &
      ngrdcol

    real(kind=time_precision), intent(in) ::  & 
      time_current, & ! Current time              [s]
      time_initial    ! Start time of model run   [s]

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) ::  & 
      z,         & ! Height at zt=2      [m] 
      rho_sfc,   & ! Density at zm=1     [kg/m^3] 
      thlm_sfc,  & ! thlm at (2)         [m/s]
      ubar

    ! Output variables
    real( kind = core_rknd ), dimension(ngrdcol), intent(out) ::  & 
      wpthlp_sfc,   & ! w'th_l' at (1)   [(m K)/s]  
      wprtp_sfc,    & ! w'r_t'(1) at (1) [(m kg)/(s kg)]
      ustar           ! surface friction velocity [m/s]

    ! Local variables
    real(kind=time_precision) ::  & 
      time      ! Elapsed time of model run    [s]

    real( kind = core_rknd ) :: ft, bflx

    integer :: i

    ! Compute heat and moisture fluxes
    ! From Table A.1.
    time = time_current - time_initial

    ft = real( max( 0._core_rknd,  & 
                 cos( 0.5_core_rknd * pi * ( (5.25_core_rknd - &
                 real( time,kind=core_rknd)/sec_per_hr) / 5.25_core_rknd ) ) & 
            ), kind = core_rknd ) ! Known magic number

    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol

      ! Known magic numbers
      wpthlp_sfc(i) =  convert_sens_ht_to_km_s( ( 270._core_rknd * ft**1.5_core_rknd ), rho_sfc(i) ) 
      wprtp_sfc(i)  =  convert_latent_ht_to_m_s( ( 554._core_rknd * ft**1.3_core_rknd ), rho_sfc(i) )

      bflx = grav / thlm_sfc(i) * wpthlp_sfc(i)

      ! Compute ustar
      ustar(i) = diag_ustar( z(i), bflx, ubar(i), z0 )

    end do

    return

  end subroutine lba_sfclyr


end module lba

