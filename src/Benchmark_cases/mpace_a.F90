!----------------------------------------------------------------------
! $Id$
module mpace_a

! Description:
!   Contains subroutines for the mpace_a intercomparison.
! References:
!   http://science.arm.gov/wg/cpm/scm/scmic5/index.html
!----------------------------------------------------------------------

  use clubb_precision, only: &
    core_rknd ! Variable(s)

  implicit none

  public :: mpace_a_tndcy, mpace_a_sfclyr, mpace_a_init

  private ! Default Scope


  ! These variables were moved here so that they could be
  ! accessible to all subroutines in mpace_a
  ! Joshua Fasching December 2007

  integer, parameter :: file_ntimes = 139
  integer, parameter :: file_nlevels = 38
  integer, parameter :: per_line = 5

  real( kind = core_rknd ), dimension(file_nlevels) :: file_pressure
  real( kind = core_rknd ), dimension(file_nlevels) :: file_heights
  real( kind = core_rknd ), dimension(file_ntimes) :: file_times


! Michael Falk is, on 28 September 2007, removing omega.  We are going to try
! to force the model without specifying it, so we can do the temperature and
! moisture forcings the way Steve Klein wants us to.

!       real, dimension(file_nlevels,file_ntimes) :: omega_forcing ! mb/s
  real( kind = core_rknd ), dimension(file_nlevels,file_ntimes) :: dTdt_forcing  ! K/hr
  real( kind = core_rknd ), dimension(file_nlevels,file_ntimes) :: dqdt_forcing  ! g/kg/hr
  real( kind = core_rknd ), dimension(file_nlevels,file_ntimes) :: vertT_forcing  ! K/hr
  real( kind = core_rknd ), dimension(file_nlevels,file_ntimes) :: vertq_forcing  ! g/kg/hr
  real( kind = core_rknd ), dimension(file_nlevels,file_ntimes) :: um_obs  ! m/s
  real( kind = core_rknd ), dimension(file_nlevels,file_ntimes) :: vm_obs  ! m/s
  real( kind = core_rknd ), dimension(file_ntimes) :: file_latent_ht
  real( kind = core_rknd ), dimension(file_ntimes) :: file_sens_ht

  contains

!----------------------------------------------------------------------
  subroutine mpace_a_tndcy( ngrdcol, sclr_dim, edsclr_dim, sclr_idx, &
                            gr, time, p_in_Pa, & 
                            wm_zt, wm_zm, thlm_forcing, rtm_forcing, & 
                            um_hoc_grid, vm_hoc_grid, & 
                            sclrm_forcing, edsclrm_forcing )

!        Description:
!
!        References:
!          Liou, Wallace and Hobbs, Shettle and Weinman
!          http://science.arm.gov/wg/cpm/scm/scmic5/index.html
!-----------------------------------------------------------------------

    use constants_clubb, only: &
      Cp, &       ! Variable(s)
      Rd, & 
      sec_per_hr, &
      g_per_kg, &
      fstderr

    use grid_class, only: &
      grid ! Type

    use grid_class, only: &
      zt2zm_api ! Procedure(s)

    use interpolation, only: &
      zlinterp_fnc, &   ! Procedure(s)
      linear_interp_factor 

    use clubb_precision, only: &
      time_precision, &  ! Variable(s)
      core_rknd

    use error_code, only: &
      clubb_at_least_debug_level_api  ! Procedure

    use array_index, only: &
      sclr_idx_type

    ! Note that this subroutine is from the time_dependent_input module, but
    ! mpace_a does not have time dependent input.
    use time_dependent_input, only: &
      time_select ! Procedure(s)

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

    real(kind=time_precision), intent(in) ::  & 
      time  ! Current time of simulation      [s]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), intent(in) :: & 
      p_in_Pa  ! Pressure                               [Pa]

    ! Output Variables
    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), intent(out) ::  & 
      wm_zt,        & ! Large-scale vertical motion on t grid   [m/s]
      thlm_forcing, & ! Large-scale thlm tendency               [K/s]
      rtm_forcing     ! Large-scale rtm tendency                [kg/kg/s]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzm), intent(out) ::  & 
      wm_zm           ! Large-scale vertical motion on m grid   [m/s]

    real( kind = core_rknd ), intent(out), dimension(ngrdcol,gr%nzt,sclr_dim) :: & 
      sclrm_forcing ! Passive scalar LS tendency            [units/s]

    real( kind = core_rknd ), intent(out), dimension(ngrdcol,gr%nzt,edsclr_dim) :: & 
      edsclrm_forcing ! Eddy-passive scalar forcing         [units/s]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), intent(out) ::  & 
      um_hoc_grid,       & ! Observed wind, for nudging         [m/s]
      vm_hoc_grid       ! Observed wind, for nudging         [m/s]

    !--------------------- Local Variables ---------------------

    ! Local constants, subsidence
    real( kind = core_rknd ), parameter :: & 
      p_sfc  = 101000._core_rknd   ! Pa

    integer :: i, k ! Loop indices

    ! Local Variables, subsidence scheme
    ! real :: velocity_omega [Pa/s]

    ! Open external files (21 Aug 2007, Michael Falk)

    integer before_time,after_time
    real( kind = core_rknd ) :: ratio

    ! real, dimension(file_nlevels) :: omega_column
    real( kind = core_rknd ), dimension(file_nlevels) :: dTdt_column
    real( kind = core_rknd ), dimension(file_nlevels) :: dqdt_column
    real( kind = core_rknd ), dimension(file_nlevels) :: vertT_column
    real( kind = core_rknd ), dimension(file_nlevels) :: vertq_column
    real( kind = core_rknd ), dimension(file_nlevels) :: um_column
    real( kind = core_rknd ), dimension(file_nlevels) :: vm_column

    ! real, dimension(gr%nz) :: omega_hoc_grid
    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt) :: &
      dqdt_hoc_grid, &
      vertT_hoc_grid, &
      vertq_hoc_grid, &
      dTdt_hoc_grid

    !--------------------- Begin Code ---------------------
   
    before_time = -1
    after_time = -1

    ! Use time_select to get the indexes before and after the specified time
    ! and to get the ratio necessary for interpolation.
    call time_select( time, file_ntimes, file_times, &
                      before_time, after_time, ratio )

    ! Sanity check to ensure that time_times is sorted.
    if( clubb_at_least_debug_level_api( 1 ) .and. (before_time == -1 .or. after_time == -1) ) then
      write(fstderr,*) "file_times not sorted in mpace_a_tndcy."
    endif 

    do k = 1, file_nlevels
!        omega_column(k) = ratio *			       ! Do linear interpolation in time
!     .                      (omega_forcing(k,after_time)
!     .                      -omega_forcing(k,before_time))
!     .                     + omega_forcing(k,before_time)

      dTdt_column(k)  = linear_interp_factor( ratio, dTdt_forcing(k, after_time), &
                                              dTdt_forcing(k, before_time) )

      dqdt_column(k)  = linear_interp_factor( ratio, dqdt_forcing(k, after_time), &
                                              dqdt_forcing(k, before_time) )

      vertT_column(k) = linear_interp_factor( ratio, vertT_forcing(k,after_time), &
                                              vertT_forcing(k,before_time) )

      vertq_column(k) = linear_interp_factor( ratio, vertq_forcing(k,after_time), &
                                              vertq_forcing(k,before_time) )

      um_column(k)    = linear_interp_factor( ratio, um_obs(k, after_time), um_obs(k, before_time) )

      vm_column(k)    = linear_interp_factor( ratio, vm_obs(k, after_time), vm_obs(k, before_time) )

    end do

!     Do linear interpolation in space
!     using zlinterp_fnc
    do i = 1, ngrdcol

      dTdt_hoc_grid(i,:)  = zlinterp_fnc( gr%nzt, file_nlevels, gr%zt(i,:), & 
                                          file_heights, dTdt_column )

      dqdt_hoc_grid(i,:)  = zlinterp_fnc( gr%nzt, file_nlevels, gr%zt(i,:), & 
                                          file_heights, dqdt_column )

      vertT_hoc_grid(i,:)  = zlinterp_fnc( gr%nzt, file_nlevels, gr%zt(i,:), & 
                                            file_heights, vertT_column )

      vertq_hoc_grid(i,:)  = zlinterp_fnc( gr%nzt, file_nlevels, gr%zt(i,:), & 
                                          file_heights, vertq_column )

      um_hoc_grid(i,:)  = zlinterp_fnc( gr%nzt, file_nlevels, gr%zt(i,:), & 
                                        file_heights, um_column )

      vm_hoc_grid(i,:)  = zlinterp_fnc( gr%nzt, file_nlevels, gr%zt(i,:), & 
                                        file_heights, vm_column )
    end do

! eMFc

!-----------------------------------------------------------------------


    ! Compute vertical motion
    ! do k = 1, gr%nzt
    !   do i = 1, ngrdcol
    !     !  velocity_omega = omega_hoc_grid(k) * 100 / 3600 ! convering mb/hr to Pa/s
    !     !  wm_zt(i,k) = -velocity_omega * Rd * thvm(i,k) / p_in_Pa(i,k) / grav
    !   end do
    ! end do
    ! End of Michael Falk's obliteration of omega.

    ! Interpolation
    ! no need to interpolate since wm_zt is set to 0 above
    !wm_zm = zt2zm_api( gr%nzm, gr%nzt, ngrdcol, gr, wm_zt )

    ! Boundary condition
    do k = 1, gr%nzm
      do i = 1, ngrdcol
        ! wm_zm(i,1) = 0.0_core_rknd        ! At surface
        ! wm_zm(i,gr%nzm) = 0.0_core_rknd  ! Model top
        wm_zm(i,k) = 0.0_core_rknd        ! At surface
      end do
    end do


    ! Compute large-scale tendencies
    do k = 1,gr%nzt
      do i = 1, ngrdcol

        wm_zt(i,k) = 0._core_rknd

        thlm_forcing(i,k) = ((dTdt_hoc_grid(i,k) + vertT_hoc_grid(i,k)) & 
                        * ((p_sfc/p_in_Pa(i,k)) ** (Rd/Cp))) & 
                        / sec_per_hr ! K/s
        rtm_forcing(i,k)  = (dqdt_hoc_grid(i,k)+vertq_hoc_grid(i,k)) & 
        / g_per_kg / sec_per_hr ! g/kg/hr -> kg/kg/s

        ! Test scalars with thetal and rt if desired
        if ( sclr_idx%iisclr_thl > 0 ) sclrm_forcing(i,k,sclr_idx%iisclr_thl) = thlm_forcing(i,k)
        if ( sclr_idx%iisclr_rt  > 0 ) sclrm_forcing(i,k,sclr_idx%iisclr_rt)  = rtm_forcing(i,k)

        if ( sclr_idx%iiedsclr_thl > 0 ) edsclrm_forcing(i,k,sclr_idx%iiedsclr_thl) = thlm_forcing(i,k)
        if ( sclr_idx%iiedsclr_rt  > 0 ) edsclrm_forcing(i,k,sclr_idx%iiedsclr_rt)  = rtm_forcing(i,k)

      end do
    end do  

    return
  end subroutine mpace_a_tndcy

!----------------------------------------------------------------------
  subroutine mpace_a_sfclyr( ngrdcol, time, rho_sfc, & 
                             wpthlp_sfc, wprtp_sfc, ustar )
!        Description:
!          Surface forcing subroutine for mpace_a case.  Written
!          October 2007 by Michael Falk.
!
!        References:
!          http://science.arm.gov/wg/cpm/scm/scmic5/index.html
!-----------------------------------------------------------------------

    use constants_clubb, only: Cp, Lv, fstderr ! Constant(s)

    use clubb_precision, only: time_precision, core_rknd ! Variable(s)

    use interpolation, only: linear_interp_factor ! Procedure(s)

    use error_code, only: clubb_at_least_debug_level_api  ! Procedure

    ! Note that this subroutine is from time_dependent_input, but 
    ! mpace_a does not have time_dependent input.
    use time_dependent_input, only: time_select ! Procedure(s)

    implicit none

    ! External
    intrinsic :: max, sqrt, present

    ! Input Variables
    integer, intent(in) :: &
      ngrdcol

    real(kind=time_precision), intent(in) :: & 
      time     ! current model time           [s]

    real( kind = core_rknd ), dimension(ngrdcol), intent(in)  :: & 
    rho_sfc     ! Air density at surface       [kg/m^3]

    ! Output Variables
    real( kind = core_rknd ), dimension(ngrdcol), intent(out) ::  & 
      wpthlp_sfc,   & ! w'th_l' at (1)   [(m K)/s]  
      wprtp_sfc,    & ! w'r_t' at (1)    [(m kg)/(s kg)]
      ustar           ! surface friction velocity [m/s]

    ! Local Variables
    real( kind = core_rknd ) :: & 
      latent_heat_flx, & 
      sensible_heat_flx

    integer :: & 
      before_time, after_time, i

    real( kind = core_rknd ) :: ratio
    !-----------------------------------------------------------------------

    before_time = -1
    after_time = -1

    ! choose which times to use
    call time_select( time, file_ntimes, file_times, &
                      before_time, after_time, ratio )

    ! Sanity check to make certain that the values read into
    ! file_times are sorted. Joshua Fasching June 2008
    if ( clubb_at_least_debug_level_api( 1 ) .and. (before_time == -1 .or. after_time == -1) ) then
      write(fstderr,*) "file_times not sorted in MPACE_A"
    endif

    latent_heat_flx = linear_interp_factor( ratio, file_latent_ht(after_time), &
                                            file_latent_ht(before_time) )

    sensible_heat_flx = linear_interp_factor( ratio, file_sens_ht(after_time), &
                                              file_sens_ht(before_time) )

    ! Compute heat and moisture fluxes
    do i = 1, ngrdcol
      wpthlp_sfc(i) = sensible_heat_flx / ( rho_sfc(i) * Cp )
      wprtp_sfc(i)  = latent_heat_flx / ( rho_sfc(i) * Lv )

      ! Declare the value of ustar.
      ustar(i) = 0.25_core_rknd
    end do

    return

  end subroutine mpace_a_sfclyr

!----------------------------------------------------------------
  subroutine mpace_a_init( iunit, file_path )
!
!       Description:
!       This subroutine initializes the module by reading in forcing
!       data used in the tndcy and sfclyr subroutines.
!----------------------------------------------------------------
    use file_functions, only: file_read_1d, file_read_2d ! Procedure(s)

    implicit none

    integer, intent(in) :: iunit ! File unit number

    character(len=*), intent(in) :: &
      file_path ! Path to the forcing files

    ! ---- Begin Code ----

    call file_read_1d( iunit, & 
      file_path//'mpace_a_press.dat', & 
      file_nlevels, per_line, file_pressure )

    call file_read_1d( iunit, & 
      file_path//'mpace_a_heights.dat', & 
      file_nlevels, per_line, file_heights )

    call file_read_1d( iunit, & 
      file_path//'mpace_a_times.dat', & 
      file_ntimes, per_line, file_times )

    !      call file_read_2d( iunit,
    !     . file_path//'mpace_a_omega.dat',
    !     . file_nlevels, file_ntimes, per_line, omega_forcing)

    call file_read_2d( iunit, & 
      file_path//'mpace_a_dTdt.dat', & 
      file_nlevels, file_ntimes, per_line, dTdt_forcing )

    call file_read_2d( iunit, & 
      file_path//'mpace_a_dqdt_horiz.dat', & 
      file_nlevels, file_ntimes, per_line, dqdt_forcing )

    call file_read_2d( iunit, & 
      file_path//'mpace_a_verts.dat', & 
      file_nlevels, file_ntimes, per_line, vertT_forcing )

    call file_read_2d( iunit, & 
      file_path//'mpace_a_vertq.dat', & 
      file_nlevels, file_ntimes, per_line, vertq_forcing )

    call file_read_2d( iunit, & 
      file_path//'mpace_a_um_obs.dat', & 
      file_nlevels, file_ntimes, per_line, um_obs )

    call file_read_2d( iunit, & 
      file_path//'mpace_a_vm_obs.dat', & 
      file_nlevels, file_ntimes, per_line, vm_obs )

    call file_read_1d( iunit, & 
      file_path//'mpace_a_lh.dat', & 
      file_ntimes, per_line, file_latent_ht )

    call file_read_1d( iunit, & 
      file_path//'mpace_a_sh.dat', & 
      file_ntimes, per_line, file_sens_ht )

    return
  end subroutine mpace_a_init

end module mpace_a
