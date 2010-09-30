!$Id$
!-------------------------------------------------------------------------------
module gabls2

! Description:
!   Contains subroutines for the GABLS2 case.
!-------------------------------------------------------------------------------

  implicit none

  public :: gabls2_tndcy, gabls2_sfclyr

  private ! Default Scope

  contains

!-------------------------------------------------------------------------------
  subroutine gabls2_tndcy( time, time_initial, &
                           wm_zt, wm_zm, thlm_forcing, & 
                           rtm_forcing, & 
                           sclrm_forcing, edsclrm_forcing )

! Description:
!   Subroutine to apply case-specific forcings to GABLS2 case
!   (Michael Falk, 29 Dec 2006).
!
! References:
!   http://people.su.se/~gsven/gabls/
!-------------------------------------------------------------------------------

    use parameters_model, only: sclr_dim, edsclr_dim ! Variable(s)

    use grid_class, only: gr ! Variable(s)

    use grid_class, only: zt2zm ! Procedure(s)

    use stats_precision, only: time_precision ! Variable(s)

    use array_index, only: iisclr_rt, iisclr_thl, iiedsclr_rt, iiedsclr_thl

    implicit none

    ! Input Variables
    real(kind=time_precision), intent(in) :: & 
      time,        & ! Current length of timestep      [s]
      time_initial   ! Current length of timestep      [s]

    ! Output Variables
    real, dimension(gr%nnzp), intent(out) :: & 
      wm_zt,        & ! Large-scale vertical motion on t grid   [m/s]
      wm_zm,        & ! Large-scale vertical motion on m grid   [m/s]
      thlm_forcing, & ! Large-scale thlm tendency               [K/s]
      rtm_forcing     ! Large-scale rtm tendency                [kg/kg/s]

    real, intent(out), dimension(gr%nnzp,sclr_dim) :: & 
      sclrm_forcing ! Passive scalar LS tendency            [units/s]

    real, intent(out), dimension(gr%nnzp,edsclr_dim) :: & 
      edsclrm_forcing ! Eddy-passive scalar forcing         [units vary/s]

    ! Local Variables, general
    integer :: k ! Loop index

    ! Compute vertical motion
    if (time > (time_initial + 93600.)) then ! That is, after 26 hours of model time;
      ! per GABLS2 specification
      do k=1,gr%nnzp
        if (gr%zt(k) <= 1000) then
          wm_zt(k) = -0.005 * (gr%zt(k) / 1000)
        else
          wm_zt(k) = -0.005
        end if
      end do
    else
      do k=1,gr%nnzp
        wm_zt(k) = 0.
      end do
    end if

    wm_zm = zt2zm( wm_zt )

    ! Boundary conditions on vertical motion.
    wm_zt(1) = 0.0        ! Below surface
    wm_zm(1) = 0.0        ! At surface
    wm_zm(gr%nnzp) = 0.0  ! Model top


    ! Compute large-scale horizontal temperature advection
    do k=1,gr%nnzp
      thlm_forcing(k) = 0.
    end do


    ! Compute large-scale horizontal moisture advection [g/kg/s]
    do k=1,gr%nnzp
      rtm_forcing(k) = 0.
    end do


    ! Test scalars with thetal and rt if desired
    if ( iisclr_thl > 0 ) sclrm_forcing(:,iisclr_thl) = thlm_forcing
    if ( iisclr_rt  > 0 ) sclrm_forcing(:,iisclr_rt)  = rtm_forcing

    if ( iiedsclr_thl > 0 ) edsclrm_forcing(:,iiedsclr_thl) = thlm_forcing
    if ( iiedsclr_rt  > 0 ) edsclrm_forcing(:,iiedsclr_rt)  = rtm_forcing


    return
  end subroutine gabls2_tndcy
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
  subroutine gabls2_sfclyr( time, time_initial, &
                            lowest_level, p_sfc, & 
                            ubar, thlm, rtm, exner_sfc, &
                            wpthlp_sfc, wprtp_sfc, ustar, T_sfc )
! Description:
!   Surface forcing subroutine for GABLS2 case.  Written
!   29 December 2006 by Michael Falk.
!
! References:
!   <http://people.su.se/~gsven/gabls/>
!-------------------------------------------------------------------------------

    use constants_clubb, only: Cp, Rd, p0, grav ! Variable(s)

    use saturation, only: sat_mixrat_liq ! Procedure(s)

    use stats_precision, only: time_precision ! Variable(s)

    use diag_ustar_module, only: diag_ustar ! Variable(s)

    use surface_flux, only: compute_wprtp_sfc, compute_wpthlp_sfc ! Procedure(s)
    implicit none

    ! Local constants
    real, parameter ::     & 
      C_10    = 0.0013,    & ! Drag coefficient, defined by ATEX specification
      z0      = 0.03         ! Roughness length, defined by GABLS2 specification

    ! Input variables
    real(kind=time_precision), intent(in) :: & 
      time,                & ! Time elapsed since 0.0 s          [s]
      time_initial           ! Initial time of model integration [s]

    real, intent(in) ::    & 
      p_sfc,                & ! Surface pressure [Pa]
      lowest_level,        & ! Height of lowest above-ground gridpoint [m]
      ubar,                & ! Root (u^2 + v^2), per ATEX and RICO spec.
      thlm,                & ! theta-l at the lowest above-ground model level. 
                           ! (theta = theta-l because there's no liquid in this case)  [K]
      rtm,                 &   ! rt at the lowest above-ground model level.  [kg/kg]
      exner_sfc

    ! Output variables
    real, intent(out) :: & 
      wpthlp_sfc, & ! The turbulent upward flux of theta-l            [K m/s]
      wprtp_sfc,  & ! The turbulent upward flux of rtm (total water)  [kg/kg m/s]
      ustar,      & ! surface friction velocity                       [m/s]
      T_sfc          ! Sea surface temperature [K].
    ! Local variables
    real :: & 
      Cz,                  & ! C_10 scaled to the height of the lowest 
                           ! model level. (Per ATEX spec)
      time_in_hours,       & ! time in hours from 00 local on first day of experiment 
                           ! (experiment starts at 14)
      sstheta,             & ! Sea surface potential temperature [K].
      bflx                   ! Needed for diag_ustar; equal to wpthlp_sfc * (g/theta)

    Cz   = C_10 * ((log( 10/z0 ))/(log( lowest_level/z0 ))) * & 
           ((log( 10/z0 ))/(log( lowest_level/z0 ))) ! Modification in case
    ! lowest model level isn't at 10 m,
    ! from ATEX specification
    time_in_hours = real((time - time_initial) / 3600. + 14.) ! at initial time,
    ! time_in_hours = 14
    ! (14 local; 19 UTC)


    ! Compute sea surface temperature
    if (time_in_hours <= 17.4) then
      T_sfc = -10 - (25*cos(time_in_hours*0.22 + 0.2)) ! SST in celsius per GABLS2 spec
    else if (time_in_hours <= 30.0) then
      T_sfc = (-0.54 * time_in_hours) + 15.2
    else if (time_in_hours <= 41.9) then
      T_sfc = -7 - (25*cos(time_in_hours*0.21 + 1.8))
    else if (time_in_hours <= 53.3) then
      T_sfc = (-0.37 * time_in_hours) + 18.0
    else if (time_in_hours <= 65.6) then
      T_sfc = -4 - (25*cos(time_in_hours*0.22 + 2.5))
    else
      T_sfc = 4.4
    end if

    T_sfc     = T_sfc + 273.15
    sstheta = T_sfc * ((p0 / p_sfc)**(Rd/Cp))

    ! Compute heat and moisture fluxes
    wpthlp_sfc = compute_wpthlp_sfc( Cz, ubar, thlm, T_sfc, exner_sfc ) 
    wprtp_sfc = compute_wprtp_sfc( Cz, ubar, rtm, sat_mixrat_liq(p_sfc,T_sfc) ) * 0.025

    ! 2.5% factor from
    ! GABLS2 specification

    ! Compute momentum fluxes
    bflx  = wpthlp_sfc * grav / sstheta
    ustar = diag_ustar(lowest_level,bflx,ubar,z0)

    return
  end subroutine gabls2_sfclyr

end module gabls2
