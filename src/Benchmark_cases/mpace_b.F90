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
  subroutine mpace_b_tndcy( p_in_Pa, thvm, & 
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

    use constants_clubb, only: Rd, Cp, Lv, p0, rc_tol, zero_threshold, &
                               grav, sec_per_day, g_per_kg ! Variable(s)

    use parameters_model, only: sclr_dim, edsclr_dim ! Variable(s)

    use grid_class, only: gr ! Variable(s)

    use grid_class, only: zt2zm ! Procedure(s)

    use stats_precision, only: time_precision ! Variable(s)

    use array_index, only: iiedsclr_rt, iiedsclr_thl, iisclr_rt, iisclr_thl ! Variable(s)

    implicit none

    ! Local constants, subsidence
    real, parameter :: & 
      D     = 5.8e-6,  & ! 1/s
      p_sfc  = 101000., & ! Pa
      pinv  = 85000.     ! Pa; ditto

    ! Input Variables
    real, dimension(gr%nnzp), intent(in) :: & 
      p_in_Pa, & ! Pressure                               [Pa]
      thvm       ! Virtual potential temperature          [K]

    ! Output Variables
    real, dimension(gr%nnzp), intent(out) ::  & 
      wm_zt,        & ! Large-scale vertical motion on t grid   [m/s]
      wm_zm,        & ! Large-scale vertical motion on m grid   [m/s]
      thlm_forcing, & ! Large-scale thlm tendency               [K/s]
      rtm_forcing     ! Large-scale rtm tendency                [kg/kg/s]

    real, intent(out), dimension(gr%nnzp,sclr_dim) :: & 
      sclrm_forcing ! Passive scalar LS tendency            [units/s]

    real, intent(out), dimension(gr%nnzp,edsclr_dim) :: & 
      edsclrm_forcing ! Eddy-passive scalar LS tendency     [units/s]

    ! Local Variables, subsidence scheme
    real :: & 
      velocity_omega


    ! Local Variables
    real :: t_tendency

    integer :: k ! loop index
    !-----------------------------------------------------------------------

    ! Compute vertical motion
    do k=2,gr%nnzp
      velocity_omega = min( D*(p_sfc-p_in_Pa(k)), D*(p_sfc-pinv) )
      wm_zt(k) = -velocity_omega * Rd * thvm(k) / p_in_Pa(k) / grav
    end do

    ! Boundary condition
    wm_zt(1) = 0.0        ! Below surface

    ! Interpolate
    wm_zm = zt2zm( wm_zt )

    ! Boundary conditions
    wm_zm(1) = 0.0        ! At surface
    wm_zm(gr%nnzp) = 0.0  ! Model top


    ! Compute large-scale tendencies
    do k=1,gr%nnzp
     t_tendency = min( -4.,-15.*(1.-((p_sfc-p_in_Pa(k))/21818.)) ) ! K/day - known magic number
     thlm_forcing(k) = (t_tendency * ((p_sfc/p_in_Pa(k)) ** (Rd/Cp)))  & 
                      / real(sec_per_day) ! K/s
     rtm_forcing(k)  = min( 0.164,-3.*(1.-((p_sfc-p_in_Pa(k))/15171.)) ) /  & 
                   g_per_kg / real(sec_per_day) ! g/kg/day -> kg/kg/s - known magic number
    end do

    ! Test scalars with thetal and rt if desired
    if ( iisclr_thl > 0 ) sclrm_forcing(:,iisclr_thl) = thlm_forcing
    if ( iisclr_rt  > 0 ) sclrm_forcing(:,iisclr_rt)  = rtm_forcing

    if ( iiedsclr_thl > 0 ) edsclrm_forcing(:,iiedsclr_thl) = thlm_forcing
    if ( iiedsclr_rt  > 0 ) edsclrm_forcing(:,iiedsclr_rt)  = rtm_forcing

    return
  end subroutine mpace_b_tndcy

!----------------------------------------------------------------------
  subroutine mpace_b_sfclyr( time, rho_sfc, & 
                           wpthlp_sfc, wprtp_sfc, ustar )

  !        Description:
  !          Surface forcing subroutine for mpace_b case.  Written July-
  !          November 2006 by Michael Falk.
  !
  !        References:
  !          http://science.arm.gov/wg/cpm/scm/scmic5/index.html
  !-----------------------------------------------------------------------

    use constants_clubb, only: Cp, Lv ! Variable(s)

    use surface_flux, only: convert_sens_ht_to_km_s, convert_latent_ht_to_m_s ! Functions(s)

    use time_dependent_input, only: sens_ht_given, latent_ht_given, time_sfc_given,& !Variable(s)
                                    time_select ! Procedure(s)

    use interpolation, only: factor_interp ! Procedure(s)

    use stats_precision, only: time_precision ! Variable(s)

    implicit none

    ! External
    intrinsic :: max, sqrt

    ! Input Variables
    real(time_precision), intent(in) :: &
      time ! The current time [s]

    real, intent(in)  :: & 
      rho_sfc    ! Air density at surface       [kg/m^3

    ! Output Variables
    real, intent(out) ::  & 
      wpthlp_sfc,   & ! w'th_l' at (1)   [(m K)/s]  
      wprtp_sfc,    & ! w'r_t'(1) at (1) [(m kg)/(s kg)]
      ustar           ! surface friction velocity [m/s]

    real :: &
    ! The values of these are from the mpace_b specification.
      sensible_heat_flx,  & ! Sensible Heat Flux     [W m^-2] 
      latent_heat_flx, &    ! Latent Heat Flux       [W m^-2] 
      time_frac ! The time fraction used for interpolation
    ! eMFc

    integer :: &
      before_time, after_time ! The time used for interpolation

    !-----------------------------------------------------------------------

    call time_select( time, size(time_sfc_given), time_sfc_given, &
                     before_time, after_time, time_frac )

    ! Get sens_ht and latent_ht from the input.
    sensible_heat_flx = factor_interp( time_frac, sens_ht_given(after_time), &
                                       sens_ht_given(before_time) )
    latent_heat_flx = factor_interp( time_frac, latent_ht_given(after_time), &
                                       latent_ht_given(before_time) )


    ! Declare the value of ustar.
    ustar = 0.25

    ! Compute heat and moisture fluxes
    wpthlp_sfc = convert_sens_ht_to_km_s( sensible_heat_flx, rho_sfc )
    wprtp_sfc  = convert_latent_ht_to_m_s( latent_heat_flx, rho_sfc )

    return
  end subroutine mpace_b_sfclyr

end module mpace_b
