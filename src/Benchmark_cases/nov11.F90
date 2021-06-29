!----------------------------------------------------------------------
! $Id$
module nov11

! Description:
!   Contains subroutines for the Nov. 11 case.

! References:
!   See below.
!----------------------------------------------------------------------

  implicit none

  public :: nov11_altocu_rtm_adjust, nov11_altocu_read_t_dependent

  ! Note: bottom point not at surface, so there is no sfc
  ! subroutine

  ! Note: nov11_altocu_tndcy has been marked as private, as the subroutine
  !       is now obsolete.

  !private :: nov11_altocu_tndcy

  private ! Default Scope

  contains

!-----------------------------------------------------------------------
  subroutine nov11_altocu_rtm_adjust( gr, time, time_initial, dt, &
                                      rtm )

! Description:
!   This subroutine performs a one-time adjustment on the total water above
!   cloud at time = 3600 seconds after the start of the simulation.
!
!   This was moved from the nov11_altocu_tndcy subroutine below, as said
!   subroutine is now obsolete and no longer needs to calculate forcings.
!
! References:
!   Larson, V. E., A. J. Smith, M. J. Falk, K. E. Kotenberg
!   and J.-C. Golaz, 2006: What determines altocumulus lifetime?
!   J. Geophys. Res., 111, D19207, doi:10.1029/2005JD007002.
!   https://pantherfile.uwm.edu/vlarson/www/journal_articles/JGR_09_smith_clex_LES.pdf
!--------------------------------------------------------------------------

    use grid_class, only: grid
      
    use clubb_precision, only: &
      time_precision, & ! Variable(s)
      core_rknd

    implicit none

    type (grid), target, intent(in) :: gr

    ! Input variables
    real(kind=time_precision), intent(in) :: & 
      time,            & ! Current time          [s]
      time_initial       ! Initial time          [s]

    real(kind=core_rknd), intent(in) :: & 
      dt              ! Timestep              [s]

    ! Input/Output variables
    real( kind = core_rknd ), intent(inout), dimension(gr%nz) :: & 
      rtm     ! Total water mixing ratio      [kg/kg]

    ! Local variables
    integer :: &
      k       ! Used for iterating over the vertical domain

    ! ---- Begin Code ----

!-----------------------------------------------------------------------
! SPECIAL NOV.11 CONDITION FOR TOTAL WATER ABOVE CLOUD
! One hour after the initial time, the total water above cloud
! is adjusted to be 0.89 of what it previously was.
!
! The conditional statement here is set so that if the timestep
! is such that there is no timestep at exactly 3600.0 seconds,
! then the operation still happnens at the first timestep and
! only the first timestep after 3600.0 seconds.
!-----------------------------------------------------------------------
    if ( time >= time_initial + 3600.0_time_precision  .and. & 
         time <  time_initial + 3600.0_time_precision + real(dt,kind=time_precision) ) then

      do k = 1, gr%nz, 1
        if ( gr%zt(k) > ( 2900.0_core_rknd + gr%zm(1) ) ) then
          rtm(k) = 0.89_core_rknd * rtm(k) ! Known magic number
        end if
      end do

    end if

    return

  end subroutine nov11_altocu_rtm_adjust

!-----------------------------------------------------------------------
!  subroutine nov11_altocu_tndcy( time, time_initial, dt, & 
!                                 rtm, &
!                                 wm_zt, wm_zm, thlm_forcing, rtm_forcing, & 
!                                 sclrm_forcing, edsclrm_forcing )
!
! Description:
!   Compute subsidence, radiation, and large-scale tendencies.
!
! NOTE:
!   This subroutine is obsolete; the subsidence profile has been moved
!   to nov11_altocu_forcings.in; and the rtm adjustment has been moved to a
!   separate subroutine.
!   ~EIHoppe/20110104

! References:
!   Larson, V. E., A. J. Smith, M. J. Falk, K. E. Kotenberg
!   and J.-C. Golaz, 2006: What determines altocumulus lifetime?
!   J. Geophys. Res., 111, D19207, doi:10.1029/2005JD007002.
!   https://pantherfile.uwm.edu/vlarson/www/journal_articles/JGR_09_smith_clex_LES.pdf
!----------------------------------------------------------------------
!
!
!    use grid_class, only: zt2zm ! Procedure(s)
!
!    use constants_clubb, only: fstderr ! Variable(s)
!
!    use parameters_model, only: sclr_dim, edsclr_dim ! Variable(s)
!
!    use clubb_precision, only: time_precision, core_rknd ! Variable(s)
!
!    use interpolation, only: lin_interpolate_on_grid ! Procedure(s)
!
!    use error_code, only: clubb_at_least_debug_level  ! Procedure
!
!    use array_index, only:  & 
!        iisclr_thl, iisclr_rt, iiedsclr_thl, iiedsclr_rt ! Variable(s)
!
!    implicit none
!
!    ! Local constants
!    ! Toggles for activating/deactivating forcings
!    logical, parameter ::  & 
!      l_subs_on   = .true.
!
!    ! Input variables
!    real(kind=time_precision), intent(in) :: & 
!      time,            & ! Current time          [s]
!      time_initial       ! Initial time          [s]
!
!    real(kind=core_rknd), intent(in) :: & 
!      dt              ! Timestep              [s]
!
!    ! Input/Output variables
!    real( kind = core_rknd ), intent(inout), dimension(gr%nz) :: & 
!      rtm     ! Total water mixing ratio      [kg/kg]
!
!    ! Output variables
!    real( kind = core_rknd ), intent(out), dimension(gr%nz) :: & 
!      wm_zt,           & ! Mean vertical wind on the thermo. grid  [m/s]
!      wm_zm,           & ! Mean vertical wind on the moment. grid  [m/s]
!      thlm_forcing,    & ! Theta_l forcing                         [K/s]
!      rtm_forcing        ! Total water forcing                     [kg/kg/s]
!
!    real( kind = core_rknd ), intent(out), dimension(gr%nz,sclr_dim) :: & 
!      sclrm_forcing   ! Passive scalar forcing                  [units/s]
!
!    real( kind = core_rknd ), intent(out), dimension(gr%nz,edsclr_dim) :: & 
!      edsclrm_forcing   ! Passive scalar forcing                  [units/s]
!
!
!    ! Local variables
!
!    ! Working arrays for subsidence interpolation
!    real( kind = core_rknd ), dimension(7) ::  & 
!    zsubs, & ! Heights at which wm_zt data is supplied (used for subsidence interpolation) [m]
!    wt1      ! ONLY wt1 IS NEEDED FOR NOV.11 CASE
!
!    ! Subsidence constant and variables (for Nov.11 case only)
!    real( kind = core_rknd ) :: & 
!    wmax,  & ! Defines value of maximum subsidence in profile  [cm/s]
!    z_inversion,    & ! Defines approx. height of inversion within cloud 
!           ! (subsidence is equal to wmax at this height) [m]
!    daz_inversion,  & ! Defines height above inversion (above this height 
!           ! subsidence linearly tapers off to zero)     [m]
!    dbz_inversion,  & ! Defines height above inversion (below this height, 
!           ! subsidence linearly tapers off to zero)    [m]
!    dbc,   & ! Defines height below cloud (at / below this height, we have NO subsidence) [m]
!    dac      ! Defines height above cloud (at / above this height, we have NO subsidence) [m]
!
!    ! Variable used for working within vertical arrays
!
!    integer :: k
!
!    integer :: nparam ! input for lin_interpolate_on_grid subroutine
!
!    ! ---- Begin Code ----
!
!-----------------------------------------------------------------------
!
! Subsidence Parameters
!
! FOR NOV.11 CASE
! ---------------
! The Nov.11 case uses a constant subsidence profile, initiated after
! 1 hour of model runtime.  This initial hour is used to "spinup" the
! simulation and produce a realistic cloud.
!
! In jun25.F (in mjfalk's /coamps/mod/consolidated5 directory on
! condella), subsidence varies over time.  As a result, he uses a
! number of arrays defining subsidence profiles for different times.
! Since Nov.11's subsidence does not vary with time, we only need one
! of these arrays, and the rest have been removed.  The array listing
! different times for subsidence has been removed as well.
!
! Comment by Adam Smith on 26 June 2006
!
! NOTE:
! The constant subsidence profile has been moved to a forcings specification
! in the input directory.
! Please check input/case_setups/nov11_altocu_forcings.in for the subsidence
! profile. This code is now obsolete.
! ~EIHoppe/20110104
!-----------------------------------------------------------------------
!
!    !----------------------
!    ! Subsidence Parameters
!    !----------------------
!    wmax =  -0.03_core_rknd
!    z_inversion = 2500.00_core_rknd + gr%zm(1)
!    daz_inversion = 1500.0_core_rknd
!    dbz_inversion = 2000.0_core_rknd
!    dbc =  300.0_core_rknd
!    dac =  200.0_core_rknd
!
!    zsubs(1) = 0._core_rknd + gr%zm(1)
!    zsubs(2) = z_inversion-dbz_inversion-dbc
!    zsubs(3) = z_inversion-dbz_inversion
!    zsubs(4) = z_inversion
!    zsubs(5) = z_inversion+daz_inversion
!    zsubs(6) = z_inversion+daz_inversion+dac
!    zsubs(7) = 4500._core_rknd + gr%zm(1)
!
!    wt1(1) = 0._core_rknd
!    wt1(2) = 0._core_rknd
!    wt1(3) = wmax
!    wt1(4) = wmax
!    wt1(5) = wmax
!    wt1(6) = 0._core_rknd
!    wt1(7) = 0._core_rknd
!
!-----------------------------------------------------------------------
! SPECIAL NOV.11 CONDITION FOR TOTAL WATER ABOVE CLOUD
! One hour after the initial time, the total water above cloud
! is adjusted to be 0.89 of what it previously was.
!
! The conditional statement here is set so that if the timestep
! is such that there is no timestep at exactly 3600.0 seconds,
! then the operation still happnens at the first timestep and
! only the first timestep after 3600.0 seconds.
!-----------------------------------------------------------------------
!    if ( time >= time_initial + 3600.0_time_precision  .and. & 
!         time <  time_initial + 3600.0_time_precision + real(dt,kind=time_precision) ) then
!
!      do k = 1, gr%nz, 1
!        if ( gr%zt(k) > ( 2900.0_core_rknd + gr%zm(1) ) ) then
!          rtm(k) = 0.89_core_rknd * rtm(k) ! Known magic number
!        end if
!      end do
!
!    end if
!
!
!    ! Impose no large-scale tendency on thetal.
!    ! Radiation may be computing interactively or analytically elsewhere.
!    thlm_forcing = 0._core_rknd
!
!---------------------------------------------------------------------
!
! Using linear interpolation scheme to interpolate subsidence
!
! FOR NOV.11 CASE
! ---------------
! As mentioned above, we want to implement a constant subsidence
! profile throughout the entire simulation (except for the 1 hour
! initial "spinup" period).  Because we do not have variations in
! subsidence over time, all tsubs sections (used in Jun.25 case)
! have been removed below.  Only one loop remains below, which will
! implement the interpolation subroutine, then return the w_ls
! value.  Unlike jun25, we use this w_ls value directly instead of
! using interpolation to calculate a value between time steps.
!
! DIAGRAM OF NOV. 11 SUBSIDENCE PROFILE
! -------------------------------------
!       |      <- No subsidence in this region
!       |
!       ------------ Height = (z_inversion + daz_inversion + dac)
!        \
!         \    <- Subsidence tapers linearly in this region
!          \
!           -------- Height = (z_inversion + daz_inversion)
!           |
!           |  <- Subsidence equals wmax in this region
!           |
!           -------- Height = z_inversion
!           |
!           -------- Height = (z_inversion - dbz_inversion)
!          /
!         /    <- Subsidence tapers linearly in this region
!        /
!       ------------ Height = (z_inversion - dbz_inversion - dbc)
!       |
!       |      <- No subsidence in this region
!
!
! Comment by Adam Smith on 26 June 2006
!-----------------------------------------------------------------------
!
!    do k=2,gr%nz
!      if ( (time >= time_initial + 3600.0_time_precision ) .and. l_subs_on ) then
!        if ( gr%zt(k) <= zsubs(7) ) then
!          nparam = 7
!          call lin_interpolate_on_grid( nparam, zsubs, wt1, gr%zt(k), wm_zt(k) )
!        else
!          wm_zt(k) = 0.0_core_rknd
!          if ( clubb_at_least_debug_level( 1 ) ) then
!            write(fstderr,*) "Thermodynamic grid level", k, "at height",  &
!                             gr%zt(k), "m. is above the highest level ",  &
!                             "specified in the subsidence sounding, which ",  &
!                             "is at height", zsubs(7), "m."
!            write(fstderr,*) "The value of subsidence is being set to 0 at ",  &
!                             "this altitude."
!          endif
!        endif
!      else
!        ! If time is not yet one hour, we have no subsidence
!        wm_zt(k) = 0.0_core_rknd
!      end if
!
!      wm_zt(1) = wm_zt(2)
!    end do
!
!    wm_zm = zt2zm(wm_zt)
!
!    ! Enter the final rtm tendency
!    do k = 1, gr%nz, 1
!
!      rtm_forcing(k) = 0._core_rknd
!
!    end do
!
!    ! Test scalars with thetal and rt if desired
!    if ( iisclr_thl > 0 ) sclrm_forcing(:,iisclr_thl) = thlm_forcing
!    if ( iisclr_rt  > 0 ) sclrm_forcing(:,iisclr_rt)  = rtm_forcing
!
!    if ( iiedsclr_thl > 0 ) edsclrm_forcing(:,iiedsclr_thl) = thlm_forcing
!    if ( iiedsclr_rt  > 0 ) edsclrm_forcing(:,iiedsclr_rt)  = rtm_forcing
!
!   return
!  end subroutine nov11_altocu_tndcy

!-----------------------------------------------------------------------
  subroutine nov11_altocu_read_t_dependent( time, &
                                            sens_ht, latent_ht )

! Description:
!   This subroutine reads in the values from the _surface.in file for
!   this case.
!
! References:
!   None
!--------------------------------------------------------------------------

    use clubb_precision, only: &
      time_precision, core_rknd ! Variable(s)

    use time_dependent_input, only: time_sfc_given, &             ! Variable(s)
                                    sens_ht_given, latent_ht_given, &
                                    time_select                   ! Procedure(s)

    use interpolation, only: linear_interp_factor ! Procedure(s)

    implicit none

    ! Input variables
    real(kind=time_precision), intent(in) :: & 
      time             ! Current time          [s]

    ! Output variables
    real( kind = core_rknd ), intent(out) :: & 
      sens_ht,    &   ! sensible heat flux [W/m^2]
      latent_ht         ! latent heat flux [W/m^2]

    ! Local variables
    real( kind = core_rknd ) :: &
      time_frac ! time fraction used for interpolation

    integer :: &
      before_time, after_time  ! time indexes used for interpolation


    ! ---- Begin Code ----

    ! interpolate T_sfc from time_dependent_input

    call time_select( time, size(time_sfc_given), time_sfc_given, &
                      before_time, after_time, time_frac )

    sens_ht = linear_interp_factor( time_frac, sens_ht_given(after_time), &
                                   sens_ht_given(before_time) )
    
    latent_ht = linear_interp_factor( time_frac, latent_ht_given(after_time), &
                                   latent_ht_given(before_time) )

    return

  end subroutine nov11_altocu_read_t_dependent 


end module nov11
