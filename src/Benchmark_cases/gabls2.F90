!$Id$
!-------------------------------------------------------------------------------
module gabls2

! Description:
!   Contains subroutines for the GABLS2 case.
! References:
!   http://www.misu.su.se/~gunilla/gabls/
!-------------------------------------------------------------------------------

  implicit none

  public :: gabls2_tndcy, gabls2_sfclyr

  private ! Default Scope

  contains

!-------------------------------------------------------------------------------
  subroutine gabls2_tndcy( ngrdcol, sclr_dim, edsclr_dim, sclr_idx, &
                           gr, time, time_initial, &
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

    use grid_class, only: &
      grid ! Type

    use grid_class, only: &
      zt2zm ! Procedure(s)

    use clubb_precision, only: &
      time_precision, & ! Variable(s)
      core_rknd 

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

    type (grid), intent(in) :: &
      gr

    real(kind=time_precision), intent(in) :: & 
      time,        & ! Current length of timestep      [s]
      time_initial   ! Current length of timestep      [s]

    !--------------------- Output Variables ---------------------
    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), intent(out) :: & 
      wm_zt,        & ! Large-scale vertical motion on t grid   [m/s]
      thlm_forcing, & ! Large-scale thlm tendency               [K/s]
      rtm_forcing     ! Large-scale rtm tendency                [kg/kg/s]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzm), intent(out) :: & 
      wm_zm           ! Large-scale vertical motion on m grid   [m/s]

    real( kind = core_rknd ), intent(out), dimension(ngrdcol,gr%nzt,sclr_dim) :: & 
      sclrm_forcing ! Passive scalar LS tendency            [units/s]

    real( kind = core_rknd ), intent(out), dimension(ngrdcol,gr%nzt,edsclr_dim) :: & 
      edsclrm_forcing ! Eddy-passive scalar forcing         [units vary/s]

    !--------------------- Local Variables ---------------------
    integer :: i, k ! Loop index

    !--------------------- Begin Code ---------------------

    ! Compute vertical motion
    ! 93600 seconds = 26 hours of simulation time;

    if ( time > (time_initial + 93600._time_precision ) ) then 

      ! per GABLS2 specification
      !$acc parallel loop gang vector collapse(2) default(present)
      do k=1,gr%nzt
        do i = 1, ngrdcol

          if ( gr%zt(i,k) <= 1000._core_rknd ) then
            wm_zt(i,k) = -0.005_core_rknd * (gr%zt(i,k) / 1000._core_rknd ) ! Known magic number
          else
            wm_zt(i,k) = -0.005_core_rknd
          end if

        end do
      end do

    else

      !$acc parallel loop gang vector collapse(2) default(present)
      do k=1,gr%nzt
        do i = 1, ngrdcol
          wm_zt(i,k) = 0._core_rknd
        end do
      end do
      
    end if

    wm_zm = zt2zm( gr%nzm, gr%nzt, ngrdcol, gr, wm_zt )

    ! Boundary conditions on vertical motion.
    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      wm_zm(i,1) = 0.0_core_rknd        ! At surface
      wm_zm(i,gr%nzm) = 0.0_core_rknd  ! Model top
    end do

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, gr%nzt
      do i = 1, ngrdcol

        ! Compute large-scale horizontal temperature advection
        thlm_forcing(i,k) = 0._core_rknd

        ! Compute large-scale horizontal moisture advection [g/kg/s]
        rtm_forcing(i,k) = 0._core_rknd
        
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

  end subroutine gabls2_tndcy
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
  subroutine gabls2_sfclyr( ngrdcol, time, time_initial, &
                            lowest_level, p_sfc, & 
                            ubar, thlm, rtm, exner_sfc, &
                            saturation_formula, &
                            wpthlp_sfc, wprtp_sfc, ustar, T_sfc )
! Description:
!   Surface forcing subroutine for GABLS2 case.  Written
!   29 December 2006 by Michael Falk.
!
! References:
!   http://people.su.se/~gsven/gabls/
!-------------------------------------------------------------------------------

    use constants_clubb, only: Cp, Rd, p0, grav, sec_per_hr ! Variable(s)

    use saturation, only: sat_mixrat_liq ! Procedure(s)

    use clubb_precision, only: time_precision, core_rknd ! Variable(s)

    use diag_ustar_module, only: diag_ustar ! Variable(s)

    use sfc_flux, only: compute_wprtp_sfc, compute_wpthlp_sfc ! Procedure(s)
    implicit none

    ! Local constants
    real( kind = core_rknd ), parameter :: & 
      standard_flux_alt = 10._core_rknd,& ! Default height at which the surface flux is computed [m]
      C_10    = 0.0013_core_rknd,        & ! Drag coefficient, defined by ATEX specification
      z0      = 0.03_core_rknd             ! Roughness length, defined by GABLS2 specification

    ! Input variables
    integer, intent(in) :: &
      ngrdcol

    real(kind=time_precision), intent(in) :: & 
      time,                & ! Time elapsed since 0.0 s          [s]
      time_initial           ! Initial time of model integration [s]

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) ::    & 
      p_sfc,                & ! Surface pressure [Pa]
      lowest_level,        & ! Height of lowest above-ground gridpoint [m]
      ubar,                & ! Root (u^2 + v^2), per ATEX and RICO spec.
      thlm,                & ! theta-l at the lowest above-ground model level. 
                           ! (theta = theta-l because there's no liquid in this case)  [K]
      rtm,                 &   ! rt at the lowest above-ground model level.  [kg/kg]
      exner_sfc

    integer, intent(in) :: &
      saturation_formula ! Integer that stores the saturation formula to be used

    ! Output variables
    real( kind = core_rknd ), dimension(ngrdcol), intent(out) :: & 
      wpthlp_sfc, & ! The turbulent upward flux of theta-l            [K m/s]
      wprtp_sfc,  & ! The turbulent upward flux of rtm (total water)  [kg/kg m/s]
      ustar,      & ! surface friction velocity                       [m/s]
      T_sfc          ! Sea surface temperature [K].

    ! Local variables
    real( kind = core_rknd ) :: & 
      T_sfc_calc, &
      time_in_hours,       & ! time in hours from 00 local on first day of experiment 
                           ! (experiment starts at 14)
                           ! model level. (Per ATEX spec)
      sstheta,             & ! Sea surface potential temperature [K].
      bflx                   ! Needed for diag_ustar; equal to wpthlp_sfc * (g/theta)

    real( kind = core_rknd ), dimension(ngrdcol) :: & 
      Cz,                  & ! C_10 scaled to the height of the lowest 
      rsat

    real(kind=time_precision) :: & 
      time_in_hours_init ! Initial time in hours is 14

    integer :: i

    ! ---- Begin Code ----
    
    !$acc enter data create( rsat, Cz )

    ! lowest model level isn't at 10 m,
    ! from ATEX specification (Stevens, et al. 2000, eq 3)
    time_in_hours_init = 14._time_precision
    time_in_hours = real( (time - time_initial) / real(sec_per_hr,kind=time_precision) + &
      time_in_hours_init, kind = core_rknd ) 
    ! at initial time,
    ! time_in_hours = 14
    ! (14 local; 19 UTC)

    ! Compute sea surface temperature
    if ( time_in_hours <= 17.4_core_rknd ) then
      ! SST in celsius per GABLS2 spec
      ! Known magic number
      T_sfc_calc = -10._core_rknd &
                  - (25._core_rknd*cos(time_in_hours*0.22_core_rknd + 0.2_core_rknd))
    else if (time_in_hours <= 30.0_core_rknd) then
      ! Known magic number
      T_sfc_calc = (-0.54_core_rknd * time_in_hours) + 15.2_core_rknd
    else if (time_in_hours <= 41.9_core_rknd) then
      ! Known magic number
      T_sfc_calc = -7._core_rknd &
                  - (25._core_rknd*cos(time_in_hours*0.21_core_rknd + 1.8_core_rknd))
    else if (time_in_hours <= 53.3_core_rknd) then
      ! Known magic number
      T_sfc_calc = (-0.37_core_rknd * time_in_hours) + 18.0_core_rknd
    else if (time_in_hours <= 65.6_core_rknd) then
      ! Known magic number
      T_sfc_calc = -4._core_rknd &
                  - (25._core_rknd*cos(time_in_hours*0.22_core_rknd + 2.5_core_rknd))
    else
      T_sfc_calc = 4.4_core_rknd
    end if

    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol

      Cz(i) = C_10 * ((log( standard_flux_alt/z0 ))/(log( lowest_level(i)/z0 ))) * & 
            ((log( standard_flux_alt/z0 ))/(log( lowest_level(i)/z0 ))) ! Modification in case

      T_sfc(i) = T_sfc_calc + 273.15_core_rknd
      rsat(i) = sat_mixrat_liq( p_sfc(i), T_sfc(i), saturation_formula )
    end do

    ! Compute heat and moisture fluxes
    call compute_wpthlp_sfc( ngrdcol, Cz, ubar, thlm, T_sfc, exner_sfc, &
                             wpthlp_sfc ) 

    call compute_wprtp_sfc( ngrdcol, Cz, ubar, rtm, rsat, &
                            wprtp_sfc )

    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol

      ! The latent heat flux at the surface is 2.5% of its potential value
      wprtp_sfc(i) = wprtp_sfc(i) * 0.025_core_rknd ! Known magic number

      ! Compute momentum fluxes
      sstheta = T_sfc(i) * ((p0 / p_sfc(i))**(Rd/Cp))
      bflx  = wpthlp_sfc(i) * grav / sstheta
      ustar(i) = diag_ustar( lowest_level(i), bflx, ubar(i), z0 )
    end do

    !$acc exit data delete( rsat, Cz )

    return

  end subroutine gabls2_sfclyr

end module gabls2
