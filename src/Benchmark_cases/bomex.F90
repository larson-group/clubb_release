!----------------------------------------------------------------------
! $Id$
module bomex

!       Description:
!       Contains subroutines for the GCSS BOMEX case.
!
!       References:
!       <http://www.knmi.nl/~siebesma/gcss/bomexcomp.init.html>
!----------------------------------------------------------------------

  implicit none

  public :: bomex_tndcy, bomex_sfclyr

  private ! Default Scope

  contains

!----------------------------------------------------------------------
  subroutine bomex_tndcy( rtm, & 
                          thlm_forcing, rtm_forcing, & 
                          sclrm_forcing, edsclrm_forcing )
!       Description:
!       Subroutine to set theta and water tendencies for BOMEX case

!       References:
!       <http://www.knmi.nl/~siebesma/gcss/bomexcomp.init.html>
!----------------------------------------------------------------------

    use grid_class, only: gr ! Variable(s)

    use grid_class, only: zt2zm ! Procedure(s)

    use spec_hum_to_mixing_ratio, only: &
        force_spec_hum_to_mixing_ratio ! Procedure(s)

    use parameters_model, only: sclr_dim, edsclr_dim ! Variable(s)

    use array_index, only: iisclr_rt, iisclr_thl, iiedsclr_rt, iiedsclr_thl ! Variable(s)

    implicit none

    ! Input Variable
    real, intent(in), dimension(gr%nzmax) :: &
      rtm    ! Total water mixing ratio (thermodynamic levels)        [kg/kg]

    ! Output Variables
    real, intent(out), dimension(gr%nzmax) :: & 
      thlm_forcing,  & ! Liquid water potential temperature tendency  [K/s]
      rtm_forcing      ! Total water mixing ratio tendency            [kg/kg/s]

    real, intent(out), dimension(gr%nzmax,sclr_dim) :: & 
      sclrm_forcing ! Passive scalar forcing        [units vary/s]

    real, intent(out), dimension(gr%nzmax,edsclr_dim) :: & 
      edsclrm_forcing ! Eddy-passive scalar forcing [units vary/s]

    ! Local Variables
    real, dimension(gr%nzmax) :: &
      qtm_forcing  ! Specified total water spec. humidity tendency    [kg/kg/s]

    integer :: k

    ! ---- Begin Code ----

    ! Radiative theta-l tendency
    thlm_forcing = 0.0

    ! Large scale advective moisture tendency
    ! The BOMEX specifications give large-scale advective moisture tendency in
    ! terms of total water specific humidity.
    do k = 2, gr%nzmax

      if ( gr%zt(k) >= 0. .and. gr%zt(k) < 300. ) then
        qtm_forcing(k) = -1.2e-8
      else if ( gr%zt(k) >= 300. .and. gr%zt(k) < 500. ) then
        qtm_forcing(k)  & 
          = - 1.2e-8  & 
              * ( 1. - ( gr%zt(k) - 300. )/( 500. - 300. ) ) !Known magic number
      else
        qtm_forcing(k) = 0.
      end if

      ! Convert forcings from terms of total water specific humidity to terms of
      ! total water mixing ratio.
      call force_spec_hum_to_mixing_ratio( rtm(k), qtm_forcing(k), rtm_forcing(k) )

    end do


    ! Boundary conditions
    thlm_forcing(1) = 0.0  ! Below surface
    rtm_forcing(1)  = 0.0  ! Below surface

    ! Test scalars with thetal and rt if desired
    if ( iisclr_thl > 0 ) sclrm_forcing(:,iisclr_thl) = thlm_forcing
    if ( iisclr_rt  > 0 ) sclrm_forcing(:,iisclr_rt)  = rtm_forcing

    if ( iiedsclr_thl > 0 ) edsclrm_forcing(:,iiedsclr_thl) = thlm_forcing
    if ( iiedsclr_rt  > 0 ) edsclrm_forcing(:,iiedsclr_rt)  = rtm_forcing

    return
  end subroutine bomex_tndcy

!----------------------------------------------------------------------
  subroutine bomex_sfclyr( time, rtm_sfc, & 
                           wpthlp_sfc, wprtp_sfc, ustar )

!       Description:
!       This subroutine computes surface fluxes of horizontal momentum,
!       heat and moisture according to GCSS BOMEX specifications

!       References:
!       <http://www.knmi.nl/~siebesma/gcss/bomexcomp.init.html>
!----------------------------------------------------------------------

    use spec_hum_to_mixing_ratio, only: &
        flux_spec_hum_to_mixing_ratio ! Procedure(s)

    use time_dependent_input, only: &
        time_select, &  ! Procedure(s)
        wpthlp_sfc_given, wpqtp_sfc_given, &
        time_sfc_given ! Variable(s)

    use interpolation, only: &
        linear_interp_factor ! Procedure(s)

    use clubb_precision, only: time_precision ! Variable(s)

    implicit none

    ! Input Variables
    real(time_precision), intent(in) ::  & 
      time    ! the current time [s]

    real, intent(in) :: &
      rtm_sfc    ! rtm(2) [kg/kg]

    ! Output variables
    real, intent(out) ::  & 
      wpthlp_sfc,   & ! w'th_l' at (1)   [(m K)/s]  
      wprtp_sfc,    & ! w'r_t' at (1)    [(m kg)/(s kg)]
      ustar           ! surface friction velocity [m/s]

    ! Local variables
    real :: wpqtp_sfc, &  ! w'q_t' at (1)         [(m kg)/(s kg)]   
            time_frac ! The time fraction used for interpolation.

    integer :: before_time, after_time ! The time bounds used for interpolation

    ! Declare the value of ustar.
    ustar = 0.28

    ! Compute heat and moisture fluxes

    call time_select(time, size(time_sfc_given), time_sfc_given, &
                before_time, after_time, time_frac)

    wpthlp_sfc = linear_interp_factor(time_frac, wpthlp_sfc_given(after_time), &
                               wpthlp_sfc_given(before_time))

    ! The BOMEX specifications give surface moisture flux in terms of total water
    ! specific humidity.
    wpqtp_sfc  = linear_interp_factor(time_frac, wpqtp_sfc_given(after_time), &
                               wpqtp_sfc_given(before_time))


    ! Convert flux from terms of total water specific humidity to terms of total
    ! water mixing ratio.
    call flux_spec_hum_to_mixing_ratio( rtm_sfc, wpqtp_sfc, wprtp_sfc )

    return
  end subroutine bomex_sfclyr

end module bomex
