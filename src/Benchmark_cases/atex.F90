!----------------------------------------------------------------------
!$Id$
module atex

  !       Description:
  !       Contains subroutines for the GCSS ATEX case.
  !----------------------------------------------------------------------

  use clubb_precision, only: core_rknd

  use grid_class, only: &
    grid ! Type
  
  implicit none

  public :: atex_tndcy, atex_sfclyr

  private ! Default Scope

  contains

  subroutine calc_forcings( ngrdcol, gr, z_inversion, &
                            thlm_forcing, rtm_forcing )

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: ngrdcol
    type (grid), intent(in) :: gr
    real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: z_inversion

    !--------------------- Output Variables ---------------------
    real( kind = core_rknd ), intent(out), dimension(ngrdcol,gr%nzt) :: & 
      thlm_forcing, & ! Liquid water potential temperature tendency [K/s]
      rtm_forcing     ! Total water mixing ratio tendency           [kg/kg/s]

    !--------------------- Local Variables ---------------------
    integer :: i, k

    !--------------------- Begin Code ---------------------
    ! Theta-l tendency
    !$acc parallel loop gang vector collapse(2) default(present)
    do i = 1, ngrdcol
      do k = 1, gr%nzt

        if ( gr%zt(i,k) > 0._core_rknd .and. gr%zt(i,k) < z_inversion(i) ) then
          ! Known magic number
          thlm_forcing(i,k) = -1.1575e-5_core_rknd &
                              * ( 3._core_rknd - gr%zt(i,k)/z_inversion(i) ) 
        else if ( gr%zt(i,k) > z_inversion(i) .and. gr%zt(i,k) <= z_inversion(i)+300._core_rknd ) then
          ! Known magic number
          thlm_forcing(i,k) = -2.315e-5_core_rknd &
                              * ( 1._core_rknd - (gr%zt(i,k)-z_inversion(i))/300._core_rknd ) 
        else
          thlm_forcing(i,k) = 0.0_core_rknd
        end if

        ! Moisture tendency
        if ( gr%zt(i,k) > 0._core_rknd .and. gr%zt(i,k) < z_inversion(i) ) then
          ! Brian - known magic number
          rtm_forcing(i,k) = -1.58e-8_core_rknd * ( 1._core_rknd - gr%zt(i,k)/z_inversion(i) )  
        else
          ! Brian
          rtm_forcing(i,k) = 0.0_core_rknd       
        end if
      end do
    end do

  end subroutine calc_forcings

  !======================================================================
  subroutine atex_tndcy( ngrdcol, sclr_dim, edsclr_dim, sclr_idx, &
                         gr, time, time_initial, &
                         rtm, &
                         interp_from_dycore_grid_method, &
                         dycore_gr, rho_ds_zm_dycore, &
                         wm_zt, wm_zm, & 
                         thlm_forcing, rtm_forcing, & 
                         sclrm_forcing, edsclrm_forcing )
  ! Description:
  !   Subroutine to set theta-l and water tendencies for ATEX case

  ! References:
  !   B. Stevens et al., 2000: Simulations of trade-wind cumuli 
  !   under a strong inversion, J. Atmos. Sci, 58,1870-1891.
  !   http://www.atmos.washington.edu/~breth/GCSS/Stevens_etal_
  !   ATEX_JAS2001.pdf

  !----------------------------------------------------------------------

  use constants_clubb, only: &
    fstderr ! Constant(s)

  use grid_class, only: &
    grid ! Type

  use grid_class, only: &
    zt2zm ! Procedure(s)

  use clubb_precision, only: &
    time_precision, & ! Variable(s)
    core_rknd

  use error_code, only: &
    clubb_at_least_debug_level, &   ! Procedure
    clubb_fatal_error, &            ! Constant
    err_code                        ! Error indicator

  use array_index, only: &
    sclr_idx_type   
  
  use grid_adaptation_module, only: &
    lin_interpolate, check_consistent, &
    interpolate_forcings, check_mass_conservation, &
    check_conservation, check_remap_for_consistency

  implicit none

  !--------------------- Input Variables ---------------------
  integer, intent(in) :: &
    ngrdcol, &
    sclr_dim, & 
    edsclr_dim

  type (sclr_idx_type), intent(in) :: &
    sclr_idx

  type (grid), intent(in) :: gr

  real(kind=time_precision), intent(in) ::  & 
    time,         & ! Current time     [s]
    time_initial ! Initial time     [s]

  real( kind = core_rknd ), intent(in), dimension(ngrdcol,gr%nzt) :: & 
    rtm      ! Total water mixing ratio        [kg/kg]

  integer, intent(in) :: &
    interp_from_dycore_grid_method

  type( grid ), intent(in) :: &
    dycore_gr

  real( kind = core_rknd ), intent(in), dimension(ngrdcol,dycore_gr%nzm) :: & 
    rho_ds_zm_dycore ! Dry, static density on momentum levels on dycore grid [kg/m^3]
                     ! use this to assume the exact linear spline as the rho_ds profile

  !--------------------- Output Variables ---------------------
  real( kind = core_rknd ), intent(out), dimension(ngrdcol,gr%nzt) :: & 
    wm_zt,        & ! w wind on thermodynamic grid                [m/s]
    thlm_forcing, & ! Liquid water potential temperature tendency [K/s]
    rtm_forcing     ! Total water mixing ratio tendency           [kg/kg/s]

  real( kind = core_rknd ), intent(out), dimension(ngrdcol,gr%nzm) :: & 
    wm_zm           ! w wind on momentum grid                     [m/s]

  real( kind = core_rknd ), intent(out), dimension(ngrdcol,gr%nzt, sclr_dim) :: & 
    sclrm_forcing   ! Passive scalar tendency         [units/s]

  real( kind = core_rknd ), intent(out), dimension(ngrdcol,gr%nzt, edsclr_dim) :: & 
    edsclrm_forcing ! Eddy-passive scalar tendency    [units/s]

  !--------------------- Local Variables ---------------------
  integer :: i, k

  integer, dimension(ngrdcol) :: &
    z_lev

  integer, dimension(ngrdcol) :: &
    z_lev_dycore

  real( kind = core_rknd ), dimension(ngrdcol) :: &
    z_inversion, z_inversion_dycore

  real( kind = core_rknd ), dimension(ngrdcol,dycore_gr%nzt) :: &
    thlm_forcing_dycore, & ! Liquid water potential temperature tendency [K/s]
    rtm_forcing_dycore     ! Total water mixing ratio tendency           [kg/kg/s]

  real( kind = core_rknd ), dimension(ngrdcol,dycore_gr%nzt) :: & 
    rtm_dycore             ! Total water mixing ratio                      [kg/kg]

  !--------------------- Begin Code ---------------------

  !$acc enter data create( z_lev, z_inversion )

  if ( time >= time_initial + 5400.0_time_precision ) then

    ! Identify height of 6.5 g/kg moisture level

    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      z_lev(i) = 1
      do while ( z_lev(i) <= gr%nzt .and. rtm(i,z_lev(i)) > 6.5e-3_core_rknd )
        z_lev(i) = z_lev(i) + 1
      end do
    end do

    if ( clubb_at_least_debug_level(2) ) then

      !$acc update host( z_lev, rtm )

      do i = 1, ngrdcol
        if ( z_lev(i) == gr%nzt+1 .or. z_lev(i) == 1 ) then
          write(fstderr,*) "Identification of 6.5 g/kg level failed"
          write(fstderr,*) "Subroutine: atex_tndcy. File: atex.F"
          write(fstderr,*) "k = ", z_lev(i), " i = ", i
          write(fstderr,*) "rtm(k) = ",rtm(i,z_lev(i))
          err_code = clubb_fatal_error
          return
        end if
      end do
    end if

    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      z_inversion(i) = gr%zt(i,z_lev(i)-1)
    end do

    ! Large scale subsidence
    !$acc parallel loop gang vector collapse(2) default(present)
    do i = 1, ngrdcol
      do k = 1, gr%nzt

        if ( gr%zt(i,k) > 0._core_rknd .and. gr%zt(i,k) <= z_inversion(i) ) then
          wm_zt(i,k) = -0.0065_core_rknd * gr%zt(i,k) / z_inversion(i) ! Known magic number
        else if ( gr%zt(i,k) > z_inversion(i) .and. gr%zt(i,k) <= z_inversion(i)+300._core_rknd ) then
          wm_zt(i,k) = - 0.0065_core_rknd * ( 1._core_rknd - (gr%zt(i,k)-z_inversion(i)) &
                                                             / 300._core_rknd ) ! Known magic number
        else
          wm_zt(i,k) = 0._core_rknd
        end if

      end do
    end do

    wm_zm = zt2zm( gr%nzm, gr%nzt, ngrdcol, gr, wm_zt )

    ! Boundary conditions.
    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      wm_zm(i,1) = 0.0_core_rknd       ! At surface
      wm_zm(i,gr%nzm) = 0.0_core_rknd  ! Model top
    end do

    if ( interp_from_dycore_grid_method > 0 ) then

      do i = 1, ngrdcol
        rtm_dycore(i,:) = lin_interpolate( dycore_gr%nzt, gr%nzt, &
                                           dycore_gr%zt, gr%zt, &
                                           rtm(i,:) )
      end do

      ! calculate z_inversion for dycore grid (z_inversion_dycore) with interpolated rtm values on 
      ! dycore grid
      do i = 1, ngrdcol
        z_lev_dycore(i) = 1
        do while ( z_lev_dycore(i) <= dycore_gr%nzt .and. rtm_dycore(i,z_lev_dycore(i)) > 6.5e-3_core_rknd )
          z_lev_dycore(i) = z_lev_dycore(i) + 1
        end do
      end do

      if ( clubb_at_least_debug_level(2) ) then

        do i = 1, ngrdcol
          if ( z_lev_dycore(i) == dycore_gr%nzt+1 .or. z_lev_dycore(i) == 1 ) then
            write(fstderr,*) "Identification of 6.5 g/kg level failed on dycore grid"
            write(fstderr,*) "Subroutine: atex_tndcy. File: atex.F"
            write(fstderr,*) "k = ", z_lev_dycore(i), " i = ", i
            write(fstderr,*) "rtm_dycore(k) = ",rtm_dycore(i,z_lev_dycore(i))
            err_code = clubb_fatal_error
            return
          end if
        end do
      end if

      do i = 1, ngrdcol
        z_inversion_dycore(i) = dycore_gr%zt(i,z_lev_dycore(i)-1)
      end do

      call calc_forcings( ngrdcol, dycore_gr, z_inversion_dycore, &  ! intent(in)
                          thlm_forcing_dycore, rtm_forcing_dycore )  ! intent(out)

      call interpolate_forcings( ngrdcol, dycore_gr, gr, &                   ! intent(in)
                                 dycore_gr%nzm, rho_ds_zm_dycore, &          ! intent(in)
                                 dycore_gr%zm, &                             ! intent(in)
                                 thlm_forcing_dycore, rtm_forcing_dycore, &  ! intent(in)
                                 thlm_forcing, rtm_forcing )                 ! intent(out)

      if ( clubb_at_least_debug_level( 2 ) ) then
        ! TODO, if grid adjusts, also check over time steps

        do i = 1, ngrdcol
          ! checks if the mass over the physics and dycore grid is the same
          call check_mass_conservation( gr, dycore_gr, &                         ! intent(in)
                                        dycore_gr%nzm, rho_ds_zm_dycore(i,:), &  ! intent(in)
                                        dycore_gr%zm )                           ! intent(in)
          ! checks whether the vertical integral of thlm_forcing and rtm_forcing is the same on
          ! the physics and dycore grid
          call check_conservation( dycore_gr%nzm, rho_ds_zm_dycore(i,:), &       ! intent(in)
                                   dycore_gr%zm, &                               ! intent(in)
                                   gr%nzm, dycore_gr%nzm, &                      ! intent(in)
                                   gr%zm, dycore_gr%zm, &                        ! intent(in)
                                   thlm_forcing(i,:), thlm_forcing_dycore(i,:) ) ! intent(in)
          call check_conservation( dycore_gr%nzm, rho_ds_zm_dycore(i,:), &      ! intent(in)
                                   dycore_gr%zm, &                              ! intent(in)
                                   gr%nzm, dycore_gr%nzm, &                     ! intent(in)
                                   gr%zm, dycore_gr%zm, &                       ! intent(in)
                                   rtm_forcing(i,:), rtm_forcing_dycore(i,:) )  ! intent(in)
          call check_remap_for_consistency( dycore_gr, gr, &                        ! intent(in)
                                            dycore_gr%nzm, rho_ds_zm_dycore(i,:), & ! intent(in)
                                            dycore_gr%zm )                          ! intent(in)
        end do
      end if

    else

      call calc_forcings( ngrdcol, gr, z_inversion, &                ! intent(in)
                          thlm_forcing, rtm_forcing )                ! intent(out)
    
    end if

  else 

    ! Forcings are applied only after t = 5400 s
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, gr%nzt
      do i = 1, ngrdcol
        wm_zt(i,k)        = 0._core_rknd
        thlm_forcing(i,k) = 0._core_rknd
        rtm_forcing(i,k)  = 0._core_rknd
      end do
    end do

    if (interp_from_dycore_grid_method > 0 ) then
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, dycore_gr%nzt
        do i = 1, ngrdcol
          thlm_forcing_dycore(i,k) = 0._core_rknd
          rtm_forcing_dycore(i,k)  = 0._core_rknd
        end do
      end do
    end if

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, gr%nzm
      do i = 1, ngrdcol
        wm_zm(i,k) = 0._core_rknd
      end do
    end do

  end if ! time >= time_initial + 5400.0_core_rknd

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

  !$acc exit data delete( z_lev, z_inversion )

  return

  end subroutine atex_tndcy

  !======================================================================
  subroutine atex_sfclyr( ngrdcol, time, ubar, & 
                          thlm_sfc, rtm_sfc, exner_sfc, & 
                          wpthlp_sfc, wprtp_sfc, ustar, T_sfc )
  ! Description:
  !   This subroutine computes surface fluxes of
  !   heat and moisture according to GCSS ATEX specifications

  ! References:   
  !   B. Stevens et al., 2000: Simulations of trade-wind cumuli 
  !   under a strong inversion, J. Atmos. Sci, 58,1870-1891.
  !   http://www.atmos.washington.edu/~breth/GCSS/Stevens_etal_
  !   ATEX_JAS2001.pdf
  !----------------------------------------------------------------------

  use sfc_flux, only: compute_wpthlp_sfc, compute_wprtp_sfc

  use interpolation, only: linear_interp_factor ! Procedure(s)

  use time_dependent_input, only: time_sfc_given, T_sfc_given, & ! Variable(s)
                                  time_select                    ! Procedure(s)

  use clubb_precision, only: time_precision, core_rknd ! Variable(s)

  implicit none

  ! Input variables
  integer, intent(in) :: &
    ngrdcol

  real(time_precision), intent(in) :: &
    time       ! the current time [s]

  real( kind = core_rknd ), dimension(ngrdcol), intent(in) ::  &
    ubar,    & ! mean sfc wind speed                           [m/s]
    thlm_sfc,& ! theta_l at first model layer                  [K]
    rtm_sfc, & ! Total water mixing ratio at first model layer [kg/kg]
    exner_sfc  ! Exner function                                [-]

  ! Output variables
  real( kind = core_rknd ), dimension(ngrdcol), intent(out) ::  & 
    wpthlp_sfc,  & ! w'theta_l' surface flux   [(m K)/s]
    wprtp_sfc,   &    ! w'rt' surface flux        [(m kg)/(kg s)]
    ustar,       &
    T_sfc          ! Surface temperature                           [K]
    
  ! Local Variable
  real( kind = core_rknd ) :: & 
    time_frac, & ! Time fraction used for interpolation
    T_sfc_interp

  real( kind = core_rknd ), dimension(ngrdcol) :: & 
    C_10, &      ! Coefficient
    adjustment   ! The adjustment for compute_wprtp_sfc

  integer :: &
    before_time, after_time, i ! The 

  !-----------------BEGIN CODE-----------------------

  !$acc enter data create( C_10, adjustment )

  ! Interpolate T_sfc from time_dependent_input

  call time_select( time, size(time_sfc_given), time_sfc_given, &
                    before_time, after_time, time_frac )

  T_sfc_interp = linear_interp_factor( time_frac, T_sfc_given(after_time), &
                                       T_sfc_given(before_time) )
                
  !$acc parallel loop gang vector default(present)
  do i = 1, ngrdcol
    C_10(i)       = 0.0013_core_rknd
    adjustment(i) = 0.0198293_core_rknd
    ustar(i)      = 0.3_core_rknd
    T_sfc(i)      = T_sfc_interp
  end do

  ! Compute wpthlp_sfc and wprtp_sfc
  call compute_wpthlp_sfc( ngrdcol, C_10, ubar, thlm_sfc, T_sfc, exner_sfc, &
                           wpthlp_sfc ) 

  call compute_wprtp_sfc( ngrdcol, C_10, ubar, rtm_sfc, adjustment, &
                          wprtp_sfc )

  !$acc exit data delete( C_10, adjustment )

  return

  end subroutine atex_sfclyr

!----------------------------------------------------------------------
end module atex
