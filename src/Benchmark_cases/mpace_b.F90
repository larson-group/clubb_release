!----------------------------------------------------------------------
! $Id$
module mpace_b

!       Description:
!       Contains subroutines for the mpace_b intercomparison.
!      
!       References:
!       http://science.arm.gov/wg/cpm/scm/scmic5/index.html
!----------------------------------------------------------------------

  implicit none

  public :: mpace_b_tndcy, mpace_b_sfclyr

  private ! Default Scope

  contains

  !----------------------------------------------------------------------
  subroutine mpace_b_tndcy( ngrdcol, sclr_dim, edsclr_dim, sclr_idx, &
                            gr, p_in_Pa, thvm, & 
                            wm_zt, wm_zm, thlm_forcing, rtm_forcing, & 
                            sclrm_forcing, edsclrm_forcing )

  !        Description:
  !          Subroutine to large-scale subsidence for mpace_b case (Michael
  !          Falk, 21 July 2006).  Added ice and radiation based on Adam
  !          Smith Nov 11 case, 27 July 2006.  Comments and documentation
  !          added 31 July 2006.
  !
  !        References:
  !          Liou, Wallace and Hobbs, Shettle and Weinman
  !          http://science.arm.gov/wg/cpm/scm/scmic5/index.html
  !-----------------------------------------------------------------------

    use constants_clubb, only: &
      Rd, & ! Variable(s)
      Cp, &
      grav, &
      sec_per_day, &
      g_per_kg 

    use grid_class, only: &
      grid ! Type

    use grid_class, only: &
      zt2zm_api ! Procedure(s)

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    use array_index, only: &
      sclr_idx_type

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
      ngrdcol, &
      sclr_dim, & 
      edsclr_dim

    type (sclr_idx_type), intent(in) :: &
      sclr_idx

    type (grid), intent(in) :: gr

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), intent(in) :: & 
      p_in_Pa, & ! Pressure                               [Pa]
      thvm       ! Virtual potential temperature          [K]

    !--------------------- Output Variables ---------------------
    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), intent(out) ::  & 
      wm_zt,        & ! Large-scale vertical motion on t grid   [m/s]
      thlm_forcing, & ! Large-scale thlm tendency               [K/s]
      rtm_forcing     ! Large-scale rtm tendency                [kg/kg/s]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzm), intent(out) ::  & 
      wm_zm           ! Large-scale vertical motion on m grid   [m/s]

    real( kind = core_rknd ), intent(out), dimension(ngrdcol,gr%nzt,sclr_dim) :: & 
      sclrm_forcing ! Passive scalar LS tendency            [units/s]

    real( kind = core_rknd ), intent(out), dimension(ngrdcol,gr%nzt,edsclr_dim) :: & 
      edsclrm_forcing ! Eddy-passive scalar LS tendency     [units/s]

    !--------------------- Local Variables ---------------------

    ! Local constants, subsidence
    real( kind = core_rknd ), parameter :: & 
      D     = 5.8e-6_core_rknd,  & ! 1/s
      p_sfc = 101000._core_rknd, & ! Pa
      pinv  = 85000._core_rknd     ! Pa; ditto

    real( kind = core_rknd ) :: & 
      velocity_omega

    real( kind = core_rknd ) :: &
      t_tendency

    integer :: k, i ! loop index

    !--------------------- Begin Code ---------------------

    !$acc enter data create( velocity_omega, t_tendency )

    ! Compute vertical motion
    !$acc parallel loop gang vector collapse(2) default(present)
    do k=1,gr%nzt
      do i = 1, ngrdcol
        velocity_omega = min( D*(p_sfc-p_in_Pa(i,k)), D*(p_sfc-pinv) )
        wm_zt(i,k) = -velocity_omega * Rd * thvm(i,k) / p_in_Pa(i,k) / grav
      end do
    end do

    ! Interpolate
    wm_zm = zt2zm_api( gr%nzm, gr%nzt, ngrdcol, gr, wm_zt )

    ! Boundary conditions
    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      wm_zm(i,1) = 0.0_core_rknd        ! At surface
      wm_zm(i,gr%nzm) = 0.0_core_rknd  ! Model top
    end do

    ! Compute large-scale tendencies
    !$acc parallel loop gang vector collapse(2) default(present)
    do k=1,gr%nzt
      do i = 1, ngrdcol
        t_tendency = min( -4._core_rknd,-15._core_rknd* &
          (1._core_rknd-((p_sfc-p_in_Pa(i,k))/21818._core_rknd)) ) ! K/day - known magic number
        thlm_forcing(i,k) = (t_tendency * ((p_sfc/p_in_Pa(i,k)) ** (Rd/Cp)))  & 
                          / sec_per_day ! K/s
        rtm_forcing(i,k)  = min( 0.164_core_rknd,-3._core_rknd* &
          (1._core_rknd-((p_sfc-p_in_Pa(i,k))/15171._core_rknd)) ) /  & 
                      g_per_kg / sec_per_day ! known magic number
                      ! g/kg/day -> kg/kg/s
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

    !$acc exit data delete( velocity_omega, t_tendency )

    return

  end subroutine mpace_b_tndcy

!----------------------------------------------------------------------
  subroutine mpace_b_sfclyr( ngrdcol, time, rho_sfc, & 
                             wpthlp_sfc, wprtp_sfc, ustar )

  !        Description:
  !          Surface forcing subroutine for mpace_b case.  Written July-
  !          November 2006 by Michael Falk.
  !
  !        References:
  !          http://science.arm.gov/wg/cpm/scm/scmic5/index.html
  !-----------------------------------------------------------------------

    use sfc_flux, only: convert_sens_ht_to_km_s, convert_latent_ht_to_m_s ! Functions(s)

    use time_dependent_input, only: sens_ht_given, latent_ht_given, time_sfc_given,& !Variable(s)
                                    time_select ! Procedure(s)

    use interpolation, only: linear_interp_factor ! Procedure(s)

    use clubb_precision, only: time_precision, core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: max, sqrt

    ! Input Variables
    integer, intent(in) :: &
      ngrdcol

    real(time_precision), intent(in) :: &
      time ! The current time [s]

    real( kind = core_rknd ), dimension(ngrdcol), intent(in)  :: & 
      rho_sfc    ! Air density at surface       [kg/m^3

    ! Output Variables
    real( kind = core_rknd ), dimension(ngrdcol), intent(out) ::  & 
      wpthlp_sfc,   & ! w'th_l' at (1)   [(m K)/s]  
      wprtp_sfc,    & ! w'r_t'(1) at (1) [(m kg)/(s kg)]
      ustar           ! surface friction velocity [m/s]

    real( kind = core_rknd ) :: &
    ! The values of these are from the mpace_b specification.
      sensible_heat_flx,  & ! Sensible Heat Flux     [W m^-2] 
      latent_heat_flx, &    ! Latent Heat Flux       [W m^-2] 
      time_frac ! The time fraction used for interpolation
    ! eMFc

    integer :: &
      before_time, after_time, i ! The time used for interpolation

    !-----------------------------------------------------------------------

    call time_select( time, size(time_sfc_given), time_sfc_given, &
                     before_time, after_time, time_frac )

    ! Get sens_ht and latent_ht from the input.
    sensible_heat_flx = linear_interp_factor( time_frac, sens_ht_given(after_time), &
                                       sens_ht_given(before_time) )
    latent_heat_flx = linear_interp_factor( time_frac, latent_ht_given(after_time), &
                                       latent_ht_given(before_time) )

    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      ! Declare the value of ustar.
      ustar(i) = 0.25_core_rknd

      ! Compute heat and moisture fluxes
      wpthlp_sfc(i) = convert_sens_ht_to_km_s( sensible_heat_flx, rho_sfc(i) )
      wprtp_sfc(i)  = convert_latent_ht_to_m_s( latent_heat_flx, rho_sfc(i) )
    end do

    return
  end subroutine mpace_b_sfclyr

end module mpace_b
