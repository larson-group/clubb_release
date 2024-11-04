!----------------------------------------------------------------------
! $Id$
module dycoms2_rf02

  ! Description:
  !   Contains subroutines for the DYCOMS II RF02 case.
  !----------------------------------------------------------------------

  implicit none

  public :: dycoms2_rf02_tndcy, dycoms2_rf02_sfclyr

  private ! Default Scope

  contains

  !----------------------------------------------------------------------
  subroutine dycoms2_rf02_tndcy( ngrdcol, sclr_dim, edsclr_dim, sclr_idx, &
                                 gr, wm_zt, wm_zm,    &
                                 thlm_forcing, rtm_forcing,  & 
                                 sclrm_forcing, edsclrm_forcing )
    ! Description:
    !   Compute thlm_ls, rtm_ls and adjust subsidence as needed.

    ! References:
    !  "A single-column model intercomparison of a heavily drizzling
    !  stratocumulus-topped boundary layer" Wyant, Matthew C., et al., (2007)
    !  J, Geophys. Res., 112, D24204,
    !  ftp://eos.atmos.washington.edu/pub/breth/papers/2007/GCSS-DYCOMS2-SCM.pdf
    !
    !  ``Dynamics and Chemistry of Marine Stratocumulus -- DYCOMS-II''
    !  Stevens, Bjorn, et al., (2003)
    !  Bull. Amer. Meteorol. Soc., 84, 579-593.
    !  http://www.atmos.ucla.edu/~bstevens/Documents/dycoms.pdf 
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

    !--------------------- InOut Variables ---------------------
    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), intent(inout) :: &
      wm_zt    ! W wind component at thermodynamic levels   [m/s]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzm), intent(inout) :: &
      wm_zm    ! W wind component at momentum levels        [m/s]

    !--------------------- Output Variables ---------------------
    real( kind = core_rknd ), intent(out), dimension(ngrdcol,gr%nzt) ::  & 
      thlm_forcing, & ! theta_l forcing                [K/s]
      rtm_forcing     ! r_t forcing                    [(kg/kg)/s] 

    real( kind = core_rknd ), intent(out), dimension(ngrdcol,gr%nzt,sclr_dim) :: & 
      sclrm_forcing    ! Passive scalar tendency        [units/s]

    real( kind = core_rknd ), intent(out), dimension(ngrdcol,gr%nzt,edsclr_dim) :: & 
      edsclrm_forcing  ! Eddy-passive scalar tendency   [units/s]

    integer :: i, k

    !--------------------- Begin Code ---------------------

    ! Enter the final thlm and rtm tendency

    ! Imposed large-scale subsidence at the uppermost level.
    ! CLUBB used a "one-sided" derivative method to compute mean advection at
    ! the uppermost thermodynamic level.  In order to avoid bringing in large
    ! amounts of various quantities from above the top of the domain, set wm_zt
    ! to 0 at level gr%nz.  To stay consistent, set wm_zm to 0 at level
    ! gr%nz.
    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      wm_zt(i,gr%nzt) = 0.0_core_rknd
      wm_zm(i,gr%nzm) = 0.0_core_rknd
    end do

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, gr%nzt
      do i = 1, ngrdcol

        thlm_forcing(i,k) = 0.0_core_rknd
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
  end subroutine dycoms2_rf02_tndcy


!----------------------------------------------------------------------

  subroutine dycoms2_rf02_sfclyr( ngrdcol, time, wpthlp_sfc, wprtp_sfc, ustar )
  ! Description:
  !   This subroutine computes surface fluxes of
  !   heat and moisture according to GCSS DYCOMS II RF 02 specifications

  ! References:
  !  ``Dynamics and Chemistry of Marine Stratocumulus -- DYCOMS-II''
  !  Stevens, Bjorn, et al., (2003)
  !  Bull. Amer. Meteorol. Soc., 84, 579-593.
  !  http://www.atmos.ucla.edu/~bstevens/Documents/dycoms.pdf 
  !----------------------------------------------------------------------

    use sfc_flux, only: convert_sens_ht_to_km_s, convert_latent_ht_to_m_s ! Procedure(s)

    use time_dependent_input, only: sens_ht_given, latent_ht_given, time_sfc_given,& ! Variable(s)
                                    time_select ! Procedure(s)

    use interpolation, only: linear_interp_factor ! Procedure(s)

    use clubb_precision, only: time_precision, core_rknd ! Variable(s)

    implicit none

    integer :: &
      ngrdcol

    real(time_precision), intent(in) :: &
      time ! The current time [s]

    ! Output
    real( kind = core_rknd ), dimension(ngrdcol), intent(out) ::  & 
      wpthlp_sfc,   & ! w'th_l' at (1)   [(m K)/s]  
      wprtp_sfc,    & ! w'r_t'(1) at (1) [(m kg)/(s kg)]
      ustar           ! surface friction velocity [m/s]

    ! External
    intrinsic :: sqrt

    ! Constant parameters
    real( kind = core_rknd ) ::  & 
      sens_ht, &   ! Sensible heat flux
      latent_ht, &   ! Latent heat fluxi
      time_frac ! The time fraction for interpolation

    integer :: &
      before_time, after_time, i ! The times used for interpolation

    real( kind = core_rknd ), parameter :: &
      rho_sfc = 1.21_core_rknd ! Air density at surface

    !------------------------BEGIN CODE-----------------------------------

    call time_select( time, size( time_sfc_given ), time_sfc_given, &
                      before_time, after_time, time_frac )

    sens_ht = linear_interp_factor( time_frac, sens_ht_given(after_time), &
                                    sens_ht_given(before_time) )
                                    
    latent_ht = linear_interp_factor( time_frac, latent_ht_given(after_time), &
                                      latent_ht_given(before_time) )

    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      ! Declare the value of ustar.
      ustar(i) = 0.25_core_rknd

      wpthlp_sfc(i) = convert_sens_ht_to_km_s( sens_ht, rho_sfc )
      wprtp_sfc(i)  = convert_latent_ht_to_m_s( latent_ht, rho_sfc )
    end do

    return
  end subroutine dycoms2_rf02_sfclyr

end module dycoms2_rf02
