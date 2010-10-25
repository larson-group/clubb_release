! $Id$
module simple_rad_module

  implicit none

  public :: simple_rad, simple_rad_bomex, sunray_sw_wrap, &
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
!  k= 1   (t) ----------                   ---------- k=kk+1      (t)
!
!  k= 2   (m) ----------                   ---------- k=kk        (m)
!
!  k= 2   (t) ----------                   ---------- k=kk        (t)
!
!            .                  .                     .
!            .                  .                     .
!            .                  .                     .
!
!  k=kk-1 (m) ----------  m = momentum     ---------- k=3         (m)
!                                  levels
!  k=kk-1 (t) ----------  t = thermo       ---------- k=3         (t)
!                                  levels
!  k=kk   (m) ----------                   ---------- k=2         (m)
!
!  k=kk   (t) ----------  kk = number of   ---------- k=2         (t)
!                              vertical
!  k=kk+1 (m) ----------       heights     ---------- k=1         (m)
!
! //////////////////////// MODEL SURFACE /////////////////////////////
!                                          ---------- k=1         (t)
!
!
! The major difference in the grids is that CLUBB uses an additional
! thermodynamic level below the model "surface".  This means that all
! CLUBB thermodynamic heights are shifted down one vertical level, and
! CLUBB also has one fewer momentum level than COAMPS.  Therefore, we
! use one additional vertical level in CLUBB, to make sure that the
! vertical domain matches in both models.
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
 
  real, dimension(lba_nzrad), private :: & 
    lba_zrad ! Altitudes        [m]
 
  real, dimension(lba_nzrad,lba_ntimes), private :: & 
    lba_krad ! Radiative tendencies     [K/s]

  contains

!-------------------------------------------------------------------------------
  subroutine simple_rad( rho, rho_zm, rtm, rcm, exner,  & 
                         err_code, Frad_LW, radht_LW )
! Description:
!   A simplified radiation driver
! References:
!   None
!-------------------------------------------------------------------------------

    use grid_class, only: gr ! Variable(s)

    use grid_class, only: zt2zm ! Procedure(s)

    use constants_clubb, only: fstderr, Cp, rc_tol ! Variable(s)

    use error_code, only: clubb_rtm_level_not_found ! Variable(s)

    use stats_type, only: stat_update_var_pt ! Procedure(s)

    use stats_variables, only:  & 
        iz_inversion, sfc, l_stats_samp ! Variable(s)

    use interpolation, only: lin_int ! Procedure(s)

    use parameters_radiation, only: &
      F0,  & ! Variable(s)
      F1,  &
      l_rad_above_cloud, &
      kappa

    implicit none

    ! External
    intrinsic :: exp

    ! Constant parameters

    real, parameter ::  & 
      ls_div = 3.75e-6

    ! Input Variables

    real, intent(in), dimension(gr%nnzp) :: & 
      rho,    & ! Density on thermodynamic grid  [kg/m^3] 
      rho_zm, & ! Density on momentum grid       [kg/m^3]
      rtm,    & ! Total water mixing ratio       [kg/kg]
      rcm,    & ! Cloud water mixing ratio       [kg/kg]
      exner     ! Exner function.                [-]

    integer, intent(inout) :: err_code

    ! Output Variables
    real, intent(out), dimension(gr%nnzp) ::  & 
      Frad_LW,         & ! Radiative flux                 [W/m^2]
      radht_LW           ! Radiative heating rate         [K/s]

    ! Local Variables
    real, dimension(gr%nnzp) ::  & 
      LWP,      & ! Liquid water path
      Heaviside

    real :: z_i

    integer :: k

    ! ---- Begin Code ----

    LWP(1:gr%nnzp) = liq_water_path( gr%nnzp, rho, rcm, gr%invrs_dzt )

    do k = 1, gr%nnzp, 1

      if ( F1 /= 0 ) then
        Frad_LW(k) = F0 * exp( -kappa * LWP(k) ) & 
                + F1 * exp( -kappa * (LWP(1) - LWP(k)) )

      else ! Mathematically equivalent to the above, but computationally cheaper
        Frad_LW(k) = F0 * exp( -kappa * LWP(k) )

      end if ! F1 /= 0

    end do 

    if ( l_rad_above_cloud ) then
      ! Find the height of the isotherm rtm = 8.0 g/kg.

      k = 2
      do while ( k <= gr%nnzp .and. rtm(k) > 8.0e-3 )
        k = k + 1
      end do
      if ( k == gr%nnzp+1 .or. k == 2 ) then
        write(fstderr,*) "Identification of 8.0 g/kg level failed"
        write(fstderr,*) "Subroutine: simple_rad. " & 
          // "File: simple_rad_module.F90"
        write(fstderr,*) "k = ", k
        write(fstderr,*) "rtm(k) = ", rtm(k)
        err_code = clubb_rtm_level_not_found
        return
      end if

      z_i = lin_int( 8.0e-3, rtm(k), rtm(k-1), gr%zt(k), gr%zt(k-1) )

      ! Compute the Heaviside step function for z - z_i.
      do k = 1, gr%nnzp, 1
        if ( gr%zm(k) - z_i  <  0.0 ) then
          Heaviside(k) = 0.0
        else if ( gr%zm(k) - z_i  ==  0.0 ) then
          Heaviside(k) = 0.5
        else if ( gr%zm(k) - z_i  >  0.0 ) then
          Heaviside(k) = 1.0
        end if
      end do

      do k = 1, gr%nnzp, 1
        if ( Heaviside(k) > 0.0 ) then
          Frad_LW(k) = Frad_LW(k) & 
                  + rho_zm(k) * Cp * ls_div * Heaviside(k) & 
                    * ( 0.25 * ((gr%zm(k)-z_i)**(4.0/3.0)) & 
                  + z_i * ((gr%zm(k)-z_i)**(1.0/3.0)) )
        end if
      end do ! k=1..gr%nnzp

      ! Update surface statistics
      if ( l_stats_samp ) then

        call stat_update_var_pt( iz_inversion, 1, z_i, sfc )

      end if

    end if ! l_rad_above_cloud

    ! Compute the radiative heating rate.
    ! The radiative heating rate is defined on thermodynamic levels.

    do k = 2, gr%nnzp, 1
      radht_LW(k) = ( 1.0 / exner(k) ) * ( -1.0/(Cp*rho(k)) ) & 
               * ( Frad_LW(k) - Frad_LW(k-1) ) * gr%invrs_dzt(k)
    end do
    radht_LW(1) = radht_LW(2)

    return
  end subroutine simple_rad

!-------------------------------------------------------------------------------
  subroutine simple_rad_bomex( radht )
! Description:
!   Compute radiation as per the GCSS BOMEX specification.
! References:
!   <http://www.knmi.nl/~siebesma/gcss/bomexcomp.init.html>
!-------------------------------------------------------------------------------
    use grid_class, only: gr ! Type

    implicit none

    ! Output Variables
    real, intent(out), dimension(gr%nnzp) :: & 
      radht  ! Radiative heating rate [K/s]

    ! Local Variables
    integer :: k

    ! ---- Begin Code ----

    ! Radiative theta-l tendency
    do k = 2, gr%nnzp

      if ( gr%zt(k) >= 0. .and. gr%zt(k) < 1500. ) then
        radht(k) = -2.315e-5
      else if ( gr%zt(k) >= 1500. .and. gr%zt(k) < 2500. ) then
        radht(k) & 
          = - 2.315e-5  & 
            + 2.315e-5  & 
              * ( gr%zt(k) - 1500. ) / ( 2500. - 1500. )
      else
        radht(k) = 0.
      end if

    end do ! k=2..gr%nnzp

    ! Boundary condition
    radht(1) = 0.0

    return
  end subroutine simple_rad_bomex

!-------------------------------------------------------------------------------
  subroutine simple_rad_lba( time, radht )
! Description:
!   Compute radiation For the LBA TRMM case.  Uses a prescribed formula and
!   interpolates with respect to time.
! References:
!   None
!-------------------------------------------------------------------------------
    use grid_class, only: gr ! Type

    use interpolation, only: zlinterp_fnc ! Procedure(s)

    use interpolation, only: factor_interp ! Procedure(s)

    use stats_precision, only: time_precision ! Constant

    implicit none

    ! Input Variables
    real(kind=time_precision), intent(in) :: &
      time ! Model time [s]

    ! Output Variables
    real, dimension(gr%nnzp), intent(out) :: radht

    ! Local Variables
    real, dimension(lba_nzrad) :: radhtz
    real :: a
    integer :: i1, i2

    ! Calculate radiative heating rate
    if ( time <=  600. ) then
      radhtz = lba_krad(:,1)

    else if ( time >= lba_ntimes * 600. ) then
      radhtz = lba_krad(:,lba_ntimes)

    else
      i1 = 1
      do while ( i1 <= lba_ntimes-1 )
        i2 = i1 + 1
        if ( time >= 600. * i1 .and. time < 600. * i2  ) then
          a  = real(( time - 600. * i1 )/( 600. * i2 - 600. * i1))
          radhtz(:) = factor_interp( a, lba_krad(:,i2), lba_krad(:,i1) )
          i1     = lba_ntimes
        end if
        i1 = i2
      end do
    end if ! time <= times(1)

    ! Radiative theta-l tendency
    radht = zlinterp_fnc( gr%nnzp, lba_nzrad, gr%zt, lba_zrad, radhtz )

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
  pure function liq_water_path( nnzp, rho, rcm, invrs_dzt )

! Description:
!   Compute liquid water path
! References:
!   None
!-------------------------------------------------------------------------------

    implicit none

    ! Input Variables
    integer, intent(in) :: nnzp

    real, intent(in), dimension(nnzp) :: &
      rho, &       ! Air Density                      [kg/m^3]
      rcm, &       ! Cloud water mixing ratio         [kg/kg]
      invrs_dzt    ! Inverse of distance per level    [1/m]

    ! Output Variables
    real, dimension(nnzp) :: &
      liq_water_path ! Liquid water path

    integer :: k

    ! ---- Begin Code ----

    liq_water_path(nnzp) = 0.0

    ! Liquid water path is defined on the intermediate model levels between the
    ! rcm and rho levels (i.e. the momentum levels in CLUBB).
    do k = nnzp-1, 1, -1
       liq_water_path(k) = liq_water_path(k+1) + rcm(k+1)*rho(k+1) / invrs_dzt(k+1)
    end do ! k = nnzp..1

    return
  end function liq_water_path

!-------------------------------------------------------------------------------
  subroutine sunray_sw_wrap( Fs0, amu0, rho, rcm, & 
                             Frad_SW, radht_SW )
! Description:
!   Wrapper for the Geert Lenderink code.

! References:
!  See subroutine sunray_sw
!-------------------------------------------------------------------------------

    use grid_class, only: gr ! Variable(s)

    use grid_class, only: ddzm ! Procedure(s)

    use constants_clubb, only: Cp ! Variable(s)

    use parameters_radiation, only: &
      eff_drop_radius, & ! Variable(s)
      omega, &
      alvdr, &
      gc

    use rad_lwsw_module, only: &
      sunray_sw ! Procedure

    implicit none

    ! Constant parameters
    ! Toggle for centered/forward differencing (in interpolations)
    ! To use centered differencing, set the toggle to .true.
    ! To use forward differencing, set the toggle to  .false.
    logical, parameter :: &
      l_center = .true. 

    ! Input Variables
    real, intent(in) :: &
      Fs0, & ! [W/m^2]
      amu0   ! Cosine of the solar zenith angle [-]

    real, intent(in), dimension(gr%nnzp) :: & 
      rho,    & ! Density on thermodynamic grid  [kg/m^3] 
      rcm       ! Cloud water mixing ratio       [kg/kg]

    ! Output Variables
    real, intent(out), dimension(gr%nnzp) ::  & 
      Frad_SW,  & ! SW Radiative flux                 [W/m^2]
      radht_SW    ! SW Radiative heating rate         [K/s]

    ! Local Variables
    real, dimension(gr%nnzp-1) ::  & 
      rcm_flipped, &
      rho_flipped, &
      dzt_flipped   

    real, dimension(gr%nnzp) ::  & 
      zt_flipped, &
      zm_flipped, &
      Frad_SW_flipped

    integer :: k, kflip

    ! ---- Begin Code ----

    ! Certain arrays in sunray_sw lack a ghost point and are therefore
    ! dimension nnzp-1
    do k = 1, gr%nnzp-1
      kflip = gr%nnzp+1-k
      rcm_flipped(k) = rcm(kflip)
      rho_flipped(k) = rho(kflip)
      dzt_flipped(k) = 1.0 / gr%invrs_dzt(kflip)
    end do

    ! The zt array does have a ghost point, but it looks like it's not
    ! referenced within the sunray_sw code.  We set it anyway just in case.
    do k = 1, gr%nnzp
      kflip = gr%nnzp+1-k
      zt_flipped(k) = gr%zt(kflip)
      zm_flipped(k) = gr%zm(kflip)
    end do

    ! Call the old sunray_sw code
    call sunray_sw( rcm_flipped, rho_flipped, amu0, dzt_flipped, gr%nnzp-1, &
                    zm_flipped, zt_flipped, &
                    eff_drop_radius, real( alvdr ), gc, Fs0, omega, l_center, &
                    Frad_SW_flipped )

    ! Return the radiation flux to the CLUBB grid
    do k = 1, gr%nnzp
      Frad_SW(k) = Frad_SW_flipped(gr%nnzp+1-k)
    end do

    ! Take the derivative of the flux to compute radht_SW (see comment above)
    radht_SW(:) = (-1.0 / (rho(:) * Cp) ) * ddzm( Frad_SW(:) )

    return
  end subroutine sunray_sw_wrap

end module simple_rad_module
