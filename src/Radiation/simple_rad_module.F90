! $Id$
module simple_rad_module

  use clubb_precision, only: &
    core_rknd ! Variable(s)

  implicit none

  public :: simple_rad, simple_rad_bomex, &
    simple_rad_lba, simple_rad_lba_init

  private :: liq_water_path

  private

!-----------------------------------------------------------------------
! The simplified LW radiation code used by both various idealized cases is now
! consolidated into the simple_rad subroutine below.  The comments here are
! in regard to the differences between the COAMPS-LES simulations done in the
! past and how the idealized radiation is handled here in CLUBB.
! -dschanen 29 Sep 2010

! SPECIAL METHOD USED TO CALCULATE RADIATION
! Grid descriptions by Adam Smith, 27 June 2006
!
! In order to verify our CLUBB simulations are working properly, we
! have first developed a series of 3D simulations using the COAMPS-LES
! model.  This large-eddy simulation (LES) simulation uses specific
! methods to calculate radiation, subsidence, and other microphysical
! processes.  To make the two models simluate clouds as closely as
! possible, we use the same radiation scheme in both models.
!
! In COAMPS-LES, we use a separate subroutine, rad_lwsw, to implement
! all radiation code.  This allows the subroutine to be duplicated
! exactly in many different models.  However, the subroutine uses the
! COAMPS vertical grid.  Therefore, for CLUBB to implement this code
! correctly, we must modify some of our variable profiles before
! calling the radiation subroutine.
!
! The following diagram describes the differences in model grids
! (see "ADDITIONAL NOTE" for important grid information):
!
!       COAMPS-LES                                   CLUBB
!
!  k= 1   (m) ----------    <MODEL TOP>    ---------- k=kk+1      (m)
!
!  k= 1   (t) ----------                   ---------- k=kk        (t)
!
!  k= 2   (m) ----------                   ---------- k=kk        (m)
!
!  k= 2   (t) ----------                   ---------- k=kk-1      (t)
!
!            .                  .                     .
!            .                  .                     .
!            .                  .                     .
!
!  k=kk-1 (m) ----------  m = momentum     ---------- k=3         (m)
!                                  levels
!  k=kk-1 (t) ----------  t = thermo       ---------- k=2         (t)
!                                  levels
!  k=kk   (m) ----------                   ---------- k=2         (m)
!
!  k=kk   (t) ----------  kk = number of   ---------- k=1         (t)
!                              vertical
!  k=kk+1 (m) ----------       heights     ---------- k=1         (m)
!
! //////////////////////// MODEL SURFACE /////////////////////////////
!
! ADDITIONAL NOTE: In order to reconcile the COAMPS and CLUBB grids,
!                  the uppermost level of the COAMPS grid is ignored,
!                  eliminating the need to add level kk+1 to the CLUBB
!                  grid.  For the purposes of this code, the COAMPS
!                  grid levels are re-indexed to start from level 1
!                  at the uppermost useful COAMPS grid level (which
!                  was previously referred to as level 2 in the above
!                  diagram and is at the same altitude as CLUBB
!                  level kk).  Likewise, the lowermost COAMPS grid
!                  level is now indexed at kk (rather than kk+1 in the
!                  above diagram).  Brian Griffin; May 10, 2008.
!
! Also, the COAMPS grid indices are numbered from the top of the model
! downward, while the CLUBB grid indices are numbered from the bottom
! up.  Therefore, since we are using a COAMPS radiation scheme, we
! flip moisture and temperature profiles that are passed into the
! rad_lwsw subroutine.  The rad scheme will produce results in using
! the COAMPS grid scheme, so all radiation output will be flipped
! back to the CLUBB grid before being applied to the model.
!
! Finally, since the COAMPS scheme does not have a gridpoint below
! model surface, we add that point to all radiative output files once
! they are converted back to CLUBB setup.  This allows all averages and
! calculations to be done correctly.
!
!
! Computation of radiative fluxes on staggered grid
! Comments by Michael Falk, 16 February 2005.
!
! Frad (and its components Frad_LW and Frad_SW) should be computed on
! w points, not on mass points, which is apparent from its formulation
! and from its location in stats_sw instead of stats_sm.  The grid
! looks like this:
!
!
! -----Frad----------------------------------    k = 1  (w level)
!     /    \            |-dwm
! -LWP------radht----------------------------    k = 1  (mass level)
!     \    /            |
! -----Frad----------------------------------    k = 2  (w level)
!     /    \
! -LWP------radht----------------------------    k = 2  (mass level)
!     \    /
! -----Frad----------------------------------    k = 3  (w level)
!     /    \
! -LWP------radht----------------------------    k = 3  (mass level)
!
! If you consider Frad to take place on mass levels, then computing
! LWP is a forward difference and is only first-order accurate, while
! if Frad computed in between LWP levels, it is a centered difference
! which is second-order accurate.
!
! The coding implementation requires that Frad depend on LWP(k) and
! LWP(k-1) since the w level for a given k is at a higher altitude
! than the mass level.  radht, back on mass levels, depends on Frad(k)
! and Frad(k+1).
!
! ADDITIONAL NOTE: For clarification of terminology, a w level on the
!                  COAMPS grid is equivalent to a momentum (m) level
!                  on the CLUBB grid, and a mass level on the COAMPS
!                  grid is equivalent to a thermodynamic (t) level on
!                  the CLUBB grid.  Brian Griffin; May 10, 2008.
!
! Additionally, these computations assume that the distance between
! mass levels (dsigma) is constant, and that the w levels (spaced by
! dsigmw) always fall exactly halfway in between the mass levels.  If
! this is not the case, consider dwm to be the distance between a w
! level and the mass level below it, and dmw to be the distance
! between a mass level and the w level below it.  Then, the
! formulation for Frad_LW, for instance, would use a weighted average:
!
! (dwm/(dwm+dmw)) * lwp(k) + (dmw/(dwm+dmw)) * lwp(k-1)
! which, for dwm always == dmw, reduces to
! (1/2) * (lwp(k)) + (1/2) * (lwp(k-1))
! which is identical to the current formulation.
! ((lwp(k)+lwp(k-1))/2)
!
! ADDITIONAL NOTE: The CLUBB parameterization is now set up to be
!                  compatible with the use of a stretched
!                  (or unevenly-spaced) grid, as well as with the use
!                  of an evenly-spaced grid.  Interpolation functions
!                  are used to compute any weighted averages, rather
!                  than using general numbers such as (1/2), which is
!                  compatible only with an evenly-spaced grid.
!                  Brian Griffin; May 10, 2008.
!
!
!-----------------------------------------------------------------------

  ! These variables are for the LBA radiation
  ! Constant Parameters
  integer, parameter, private :: lba_ntimes = 36, lba_nzrad = 33
 
  real( kind = core_rknd ), dimension(lba_nzrad), private :: & 
    lba_zrad ! Altitudes        [m]
 
  real( kind = core_rknd ), dimension(lba_nzrad,lba_ntimes), private :: & 
    lba_krad ! Radiative tendencies     [K/s]

  contains

  !-------------------------------------------------------------------------------
  subroutine simple_rad( gr, ngrdcol, rho, rho_zm, rtm, rcm, exner,  &
                         stats, err_info, Frad_LW, radht_LW )
! Description:
!   A simplified radiation driver
! References:
!   None
!-------------------------------------------------------------------------------


    use grid_class, only: zt2zm_api, grid ! Procedure(s)

    use constants_clubb, only: fstderr, one, Cp, eps ! Variable(s)

    use error_code, only: &
        clubb_at_least_debug_level_api,  & ! Procedure
        clubb_fatal_error                 ! Constant

    use stats_netcdf, only: &
      stats_type, &
      stats_update

    use interpolation, only: lin_interpolate_two_points ! Procedure(s)

    use parameters_radiation, only: &
      F0,  & ! Variable(s)
      F1,  &
      l_rad_above_cloud, &
      kappa

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    use err_info_type_module, only: &
      err_info_type        ! Type

    implicit none

    ! External
    intrinsic :: exp

    ! Constant parameters

    real( kind = core_rknd ), parameter ::  &
      ls_div = 3.75e-6_core_rknd

    ! Input Variables
    type (grid), intent(in) :: gr

    integer, intent(in) :: ngrdcol

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,gr%nzt) :: &
      rho,    & ! Density on thermodynamic grid  [kg/m^3]
      rtm,    & ! Total water mixing ratio       [kg/kg]
      rcm,    & ! Cloud water mixing ratio       [kg/kg]
      exner     ! Exner function.                [-]

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,gr%nzm) :: &
      rho_zm    ! Density on momentum grid       [kg/m^3]

    type(stats_type), intent(inout) :: &
      stats

    type(err_info_type), intent(inout) :: &
      err_info        ! err_info struct containing err_code and err_header

    ! Output Variables
    real( kind = core_rknd ), intent(out), dimension(ngrdcol,gr%nzm) ::  &
      Frad_LW            ! Radiative flux                 [W/m^2]

    real( kind = core_rknd ), intent(out), dimension(ngrdcol,gr%nzt) ::  &
      radht_LW           ! Radiative heating rate         [K/s]

    ! Local Variables
    real( kind = core_rknd ), dimension(ngrdcol,gr%nzm) ::  &
      LWP,      & ! Liquid water path
      Heaviside

    real( kind = core_rknd ), dimension(ngrdcol) :: z_i

    integer :: i, k, k_iso, k_rtm

    ! ---- Begin Code ----

    !$acc enter data create( LWP, Heaviside, z_i )

    LWP(:,:) = liq_water_path( ngrdcol, gr%nzm, gr%nzt, rho(:,:), rcm(:,:), gr%invrs_dzt(:,:) )

    if ( F1 > eps ) then

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, gr%nzm
        do i = 1, ngrdcol
          Frad_LW(i,k) = F0 * exp( -kappa * LWP(i,k) ) &
                          + F1 * exp( -kappa * ( LWP(i,1) - LWP(i,k) ) )
        end do
      end do

    else ! Mathematically equivalent to the above, but computationally cheaper

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, gr%nzm
        do i = 1, ngrdcol
          Frad_LW(i,k) = F0 * exp( -kappa * LWP(i,k) )
        end do
      end do

    end if ! F1 /= 0

    if ( l_rad_above_cloud ) then

      ! Find the height of the isotherm rtm = 8.0 g/kg.

      !$acc parallel loop default(present)
      do i = 1, ngrdcol

        k_iso = 1

        do while ( k_iso <= gr%nzt .and. rtm(i,k_iso) > 8.0e-3_core_rknd )
          k_iso = k_iso + 1
        end do

        ! if ( clubb_at_least_debug_level_api( 0 ) ) then
        !   if ( k_iso == gr%nzt+1 .or. k_iso == 1 ) then
        !     write(fstderr,*) err_info%err_header_global
        !     write(fstderr,*) "Identification of 8.0 g/kg level failed"
        !     write(fstderr,*) "Subroutine: simple_rad. " &
        !       // "File: simple_rad_module.F90"
        !     write(fstderr,*) "k = ", k_iso
        !     k_rtm = max( 1, min( k_iso, gr%nzt ) )
        !     write(fstderr,*) "rtm(k) = ", rtm(i,k_rtm)
        !     ! General error -> set all entries to clubb_fatal_error
        !     if ( allocated(err_info%err_code) ) err_info%err_code(i) = clubb_fatal_error
        !     return
        !   end if
        ! end if

        z_i(i) = lin_interpolate_two_points( 8.0e-3_core_rknd, rtm(i,k_iso), rtm(i,k_iso-1), gr%zt(i,k_iso), &
                                             gr%zt(i,k_iso-1) )
      end do ! i=1..ngrdcol

      ! Compute the Heaviside step function for z - z_i.
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, gr%nzm, 1
        do i = 1, ngrdcol
          ! if gr%zm(i,k) > z_i
          if ( gr%zm(i,k)-z_i(i) < -eps ) then
            Heaviside(i,k) = 0.0_core_rknd
          ! if gr%zm(i,k) < z_i
          else if ( gr%zm(i,k)-z_i(i) > eps) then
            Heaviside(i,k) = 1.0_core_rknd
          else ! gr%zm(i,k) and z_i are equal within eps
            Heaviside(i,k) = 0.5_core_rknd
          end if
        end do
      end do

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, gr%nzm, 1
        do i = 1, ngrdcol
          if ( Heaviside(i,k) > 0.0_core_rknd ) then
            Frad_LW(i,k) = Frad_LW(i,k) &
                           + rho_zm(i,k) * Cp * ls_div * Heaviside(i,k) &
                             * ( 0.25_core_rknd * ((gr%zm(i,k)-z_i(i))**(4.0_core_rknd/3.0_core_rknd)) &
                                 + z_i(i) * ((gr%zm(i,k)-z_i(i))**(1.0_core_rknd/3.0_core_rknd)) )
          end if
        end do ! k=1..gr%nzm
      end do

      ! Update inversion-height statistics used by surface-radiation diagnostics.
      if ( stats%l_sample ) then
        !$acc update host( z_i )
        call stats_update( "z_inversion", z_i, stats )
      end if

    end if ! l_rad_above_cloud

    ! Compute the radiative heating rate.
    ! The radiative heating rate is defined on thermodynamic levels.

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, gr%nzt, 1
      do i = 1, ngrdcol
        radht_LW(i,k) = ( one / exner(i,k) ) * ( -one / (Cp*rho(i,k)) ) &
                        * ( Frad_LW(i,k+1) - Frad_LW(i,k) ) * gr%invrs_dzt(1,k)
      end do
    end do

    !$acc exit data delete( LWP, Heaviside, z_i )

    return

  end subroutine simple_rad

!-------------------------------------------------------------------------------
  subroutine simple_rad_bomex( gr, ngrdcol, radht )
! Description:
!   Compute radiation as per the GCSS BOMEX specification.
! References:
!   <http://www.knmi.nl/~siebesma/gcss/bomexcomp.init.html>
!-------------------------------------------------------------------------------

    use grid_class, only: grid

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    type (grid), intent(in) :: gr
    integer, intent(in) :: &
      ngrdcol

    ! Output Variables
    real( kind = core_rknd ), intent(out), dimension(ngrdcol,gr%nzt) :: &
      radht  ! Radiative heating rate [K/s]

    ! Local Variables
    integer :: k, i

    ! ---- Begin Code ----

    ! Radiative theta-l tendencysimple_rad_bomex
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, gr%nzt
      do i = 1, ngrdcol

        if ( gr%zt(i,k) >= 0._core_rknd .and. gr%zt(i,k) < 1500._core_rknd ) then
          radht(i,k) = -2.315e-5_core_rknd
        else if ( gr%zt(i,k) >= 1500._core_rknd .and. gr%zt(i,k) < 2500._core_rknd ) then
          ! From bomex specification, section 3.4
          radht(i,k) &
            = - 2.315e-5_core_rknd  &
              + 2.315e-5_core_rknd  &
                * ( gr%zt(i,k) - 1500._core_rknd ) / &
                ( 2500._core_rknd - 1500._core_rknd ) ! Known magic number
        else
          radht(i,k) = 0._core_rknd
        end if

      end do
    end do ! k=1..gr%nzt

    return
  end subroutine simple_rad_bomex

!-------------------------------------------------------------------------------
  subroutine simple_rad_lba( gr, ngrdcol, time_current, time_initial, radht )
! Description:
!   Compute radiation For the LBA TRMM case.  Uses a prescribed formula and
!   interpolates with respect to time.
! References:
!   None
!-------------------------------------------------------------------------------

    use interpolation, only: zlinterp_fnc ! Procedure(s)

    use grid_class, only: grid

    use interpolation, only: linear_interp_factor ! Procedure(s)

    use clubb_precision, only: time_precision, core_rknd ! Constant

    implicit none

    type (grid), intent(in) :: gr
    integer, intent(in) :: &
      ngrdcol

    ! Input Variables
    real(kind=time_precision), intent(in) :: &
      time_current, & ! Current time of model run   [s]
      time_initial    ! Start time of model run     [s]

    ! Output Variables
    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), intent(out) :: radht

    ! Local Variables
    real( kind = core_rknd ), dimension(lba_nzrad) :: radhtz
    real( kind = core_rknd ) :: a, time
    integer :: i1, i2, i

    time = real( time_current - time_initial, kind = core_rknd )

    ! Calculate radiative heating rate
    if ( time <=  600._core_rknd ) then

      radhtz = lba_krad(:,1)

    else if ( time >= real( lba_ntimes, kind = core_rknd ) * 600._core_rknd ) then

      radhtz = lba_krad(:,lba_ntimes)

    else

      i1 = 1
      do while ( i1 <= lba_ntimes-1 )

        i2 = i1 + 1

        if ( time >= 600._core_rknd * real( i1, kind = core_rknd ) .and. &
             time < 600._core_rknd * real( i2, kind = core_rknd )  ) then

          a  = ( time - 600._core_rknd * real( i1, kind = core_rknd ) ) & ! Known magic number
            /( 600._core_rknd * real( i2, kind = core_rknd ) - 600._core_rknd * &
            real( i1, kind = core_rknd )) ! Known magic number

          radhtz(:) = linear_interp_factor( a, lba_krad(:,i2), lba_krad(:,i1) )

          i1 = lba_ntimes

        end if

        i1 = i2

      end do
    end if ! time <= times(1)

    ! Radiative theta-l tendency
    do i = 1, ngrdcol
      radht(i,:) = zlinterp_fnc( gr%nzt, lba_nzrad, gr%zt(i,:), lba_zrad, radhtz )
    end do

    return
  end subroutine simple_rad_lba

  !----------------------------------------------------------------
  subroutine simple_rad_lba_init( iunit, file_path )
    !
    !       Description:
    !       This subroutine initializes the module by reading in forcing
    !       data used in the tndcy subroutine.
    !----------------------------------------------------------------

    use file_functions, only: file_read_1d, file_read_2d ! Procedure(s)

    implicit none

    ! Constant Parameters
    integer, parameter :: per_line = 5

    integer, intent(in) :: iunit ! File unit number

    character(len=*), intent(in) :: &
      file_path ! Path to the forcing files

    ! ---- Begin Code ----

    call file_read_1d( iunit, & 
      file_path//'lba_heights.dat', & 
      lba_nzrad, per_line, lba_zrad )

    call file_read_2d( iunit, & 
      file_path//'lba_rad.dat', & 
      lba_nzrad, lba_ntimes, per_line, lba_krad )

    return
  end subroutine simple_rad_lba_init

!-------------------------------------------------------------------------------
  pure function liq_water_path( ngrdcol, nzm, nzt, rho, rcm, invrs_dzt )

! Description:
!   Compute liquid water path
! References:
!   None
!-------------------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      ngrdcol, & ! Number of grid columns in the model
      nzm,     & ! Number of momentum vertical levels in the model
      nzt        ! Number of thermodynamic vertical levels in the model

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzt) :: &
      rho, &       ! Air Density                      [kg/m^3]
      rcm, &       ! Cloud water mixing ratio         [kg/kg]
      invrs_dzt    ! Inverse of distance per level    [1/m]

    ! Output Variables
    real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
      liq_water_path ! Liquid water path

    integer :: i, k

    ! ---- Begin Code ----

    ! Liquid water path is defined on the intermediate model levels between the
    ! rcm and rho levels (i.e. the momentum levels in CLUBB).
    !$acc parallel loop default(present)
    do i = 1, ngrdcol
      liq_water_path(i,nzm) = 0.0_core_rknd
      do k = nzm-1, 1, -1
        liq_water_path(i,k) = liq_water_path(i,k+1) + rcm(i,k) * rho(i,k) / invrs_dzt(i,k)
      end do
    end do ! k = nzm..1

    return
  end function liq_water_path

end module simple_rad_module
