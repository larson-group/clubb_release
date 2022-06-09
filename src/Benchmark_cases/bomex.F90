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
  subroutine bomex_tndcy( gr, rtm, & 
                          thlm_forcing, rtm_forcing, & 
                          sclrm_forcing, edsclrm_forcing )
!       Description:
!       Subroutine to set theta and water tendencies for BOMEX case

!       References:
!       <http://www.knmi.nl/~siebesma/gcss/bomexcomp.init.html>
!----------------------------------------------------------------------


    use grid_class, only: &
        zt2zm, & ! Procedure(s)
        grid     ! Type

    use spec_hum_to_mixing_ratio, only: &
        force_spec_hum_to_mixing_ratio ! Procedure(s)

    use parameters_model, only: sclr_dim, edsclr_dim ! Variable(s)

    use array_index, only: iisclr_rt, iisclr_thl, iiedsclr_rt, iiedsclr_thl ! Variable(s)

    use clubb_precision, only: core_rknd ! Variable(s)

    implicit none

    type (grid), target, intent(in) :: gr

    ! Input Variable
    real( kind = core_rknd ), intent(in), dimension(gr%nz) :: &
      rtm    ! Total water mixing ratio (thermodynamic levels)        [kg/kg]

    ! Output Variables
    real( kind = core_rknd ), intent(out), dimension(gr%nz) :: & 
      thlm_forcing,  & ! Liquid water potential temperature tendency  [K/s]
      rtm_forcing      ! Total water mixing ratio tendency            [kg/kg/s]

    real( kind = core_rknd ), intent(out), dimension(gr%nz,sclr_dim) :: & 
      sclrm_forcing ! Passive scalar forcing        [units vary/s]

    real( kind = core_rknd ), intent(out), dimension(gr%nz,edsclr_dim) :: & 
      edsclrm_forcing ! Eddy-passive scalar forcing [units vary/s]

    ! Local Variables
    real( kind = core_rknd ), dimension(gr%nz) :: &
      qtm_forcing  ! Specified total water spec. humidity tendency    [kg/kg/s]

    integer :: k

    ! ---- Begin Code ----

    ! Radiative theta-l tendency
    thlm_forcing = 0.0_core_rknd

    ! Large scale advective moisture tendency
    ! The BOMEX specifications give large-scale advective moisture tendency in
    ! terms of total water specific humidity.
    do k = 2, gr%nz

      if ( gr%zt(1,k) >= 0._core_rknd .and. gr%zt(1,k) < 300._core_rknd ) then
        qtm_forcing(k) = -1.2e-8_core_rknd
      else if ( gr%zt(1,k) >= 300._core_rknd .and. gr%zt(1,k) < 500._core_rknd ) then
        qtm_forcing(k)  & 
          = - 1.2e-8_core_rknd  & 
              * ( 1._core_rknd - ( gr%zt(1,k) - 300._core_rknd )/ &
              ( 500._core_rknd - 300._core_rknd ) ) !Known magic number
      else
        qtm_forcing(k) = 0._core_rknd
      end if

      ! Convert forcings from terms of total water specific humidity to terms of
      ! total water mixing ratio.
      call force_spec_hum_to_mixing_ratio( rtm(k), qtm_forcing(k), rtm_forcing(k) )

    end do


    ! Boundary conditions
    thlm_forcing(1) = 0.0_core_rknd  ! Below surface
    rtm_forcing(1)  = 0.0_core_rknd  ! Below surface

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

    use clubb_precision, only: time_precision, core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real(time_precision), intent(in) ::  & 
      time    ! the current time [s]

    real( kind = core_rknd ), intent(in) :: &
      rtm_sfc    ! rtm(2) [kg/kg]

    ! Output variables
    real( kind = core_rknd ), intent(out) ::  & 
      wpthlp_sfc,   & ! w'th_l' at (1)   [(m K)/s]  
      wprtp_sfc,    & ! w'r_t' at (1)    [(m kg)/(s kg)]
      ustar           ! surface friction velocity [m/s]

    ! Local variables
    real( kind = core_rknd ) :: wpqtp_sfc, &  ! w'q_t' at (1)         [(m kg)/(s kg)]   
            time_frac ! The time fraction used for interpolation.

    integer :: before_time, after_time ! The time bounds used for interpolation

    ! Declare the value of ustar.
    ustar = 0.28_core_rknd

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
