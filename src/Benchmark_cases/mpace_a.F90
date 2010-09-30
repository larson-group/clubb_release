!----------------------------------------------------------------------
! $Id$
module mpace_a

! Description:
!   Contains subroutines for the mpace_a intercomparison.
!----------------------------------------------------------------------

  implicit none

  public :: mpace_a_tndcy, mpace_a_sfclyr, mpace_a_init

  private ! Default Scope


  ! These variables were moved here so that they could be
  ! accessible to all subroutines in mpace_a
  ! Joshua Fasching December 2007

  integer, parameter :: file_ntimes = 139
  integer, parameter :: file_nlevels = 38
  integer, parameter :: per_line = 5

  real, dimension(file_nlevels) :: file_pressure
  real, dimension(file_nlevels) :: file_heights
  real, dimension(file_ntimes) :: file_times


! Michael Falk is, on 28 September 2007, removing omega.  We are going to try
! to force the model without specifying it, so we can do the temperature and
! moisture forcings the way Steve Klein wants us to.

!       real, dimension(file_nlevels,file_ntimes) :: omega_forcing ! mb/s
  real, dimension(file_nlevels,file_ntimes) :: dTdt_forcing  ! K/hr
  real, dimension(file_nlevels,file_ntimes) :: dqdt_forcing  ! g/kg/hr
  real, dimension(file_nlevels,file_ntimes) :: vertT_forcing  ! K/hr
  real, dimension(file_nlevels,file_ntimes) :: vertq_forcing  ! g/kg/hr
  real, dimension(file_nlevels,file_ntimes) :: um_obs  ! m/s
  real, dimension(file_nlevels,file_ntimes) :: vm_obs  ! m/s
  real, dimension(file_ntimes) :: file_LH
  real, dimension(file_ntimes) :: file_SH

  contains

!----------------------------------------------------------------------
  subroutine mpace_a_tndcy( time, p_in_Pa, & 
                            wm_zt, wm_zm, thlm_forcing, rtm_forcing, & 
                            um_hoc_grid, vm_hoc_grid, & 
                            sclrm_forcing, edsclrm_forcing )

!        Description:
!
!        References:
!          Liou, Wallace and Hobbs, Shettle and Weinman
!-----------------------------------------------------------------------

    use constants_clubb, only: Cp, Rd, Lv, p0, rc_tol, & ! Variable(s)
                         zero_threshold, fstderr

    use parameters_model, only: sclr_dim, edsclr_dim ! Variable(s)

    use grid_class, only: gr ! Variable(s)

    use grid_class, only: zt2zm ! Procedure(s)

    use interpolation, only: zlinterp_fnc, factor_interp ! Procedure(s)

    use stats_precision, only: time_precision ! Variable(s)

    use array_index, only: iisclr_rt, iisclr_thl, iiedsclr_rt, iiedsclr_thl ! Variable(s)

    use error_code, only: clubb_debug, clubb_at_least_debug_level ! Procedure(s)

    implicit none

    ! Local constants, subsidence
    real, parameter :: & 
    p_sfc  = 101000.   ! Pa

    ! Input Variables
    real(kind=time_precision), intent(in) ::  & 
    time  ! Current time of simulation      [s]

    real, dimension(gr%nnzp), intent(in) :: & 
    p_in_Pa  ! Pressure                               [Pa]

    ! Output Variables
    real, dimension(gr%nnzp), intent(out) ::  & 
    wm_zt,        & ! Large-scale vertical motion on t grid   [m/s]
    wm_zm,        & ! Large-scale vertical motion on m grid   [m/s]
    thlm_forcing, & ! Large-scale thlm tendency               [K/s]
    rtm_forcing     ! Large-scale rtm tendency                [kg/kg/s]

    real, intent(out), dimension(gr%nnzp,sclr_dim) :: & 
    sclrm_forcing ! Passive scalar LS tendency            [units/s]

    real, intent(out), dimension(gr%nnzp,edsclr_dim) :: & 
    edsclrm_forcing ! Eddy-passive scalar forcing         [units/s]

    ! Local Variables, general
    integer :: i, k ! Loop indices

    ! Local Variables, subsidence scheme
!        real :: velocity_omega [Pa/s]

    ! Open external files (21 Aug 2007, Michael Falk)

    integer left_time,right_time
    real :: ratio

!      real, dimension(file_nlevels) :: omega_column
    real, dimension(file_nlevels) :: dTdt_column
    real, dimension(file_nlevels) :: dqdt_column
    real, dimension(file_nlevels) :: vertT_column
    real, dimension(file_nlevels) :: vertq_column
    real, dimension(file_nlevels) :: um_column
    real, dimension(file_nlevels) :: vm_column

!      real, dimension(gr%nnzp) :: omega_hoc_grid
    real, dimension(gr%nnzp) :: dTdt_hoc_grid
    real, dimension(gr%nnzp) :: dqdt_hoc_grid
    real, dimension(gr%nnzp) :: vertT_hoc_grid
    real, dimension(gr%nnzp) :: vertq_hoc_grid

    real, dimension(gr%nnzp), intent(out) ::  & 
    um_hoc_grid,       & ! Observed wind, for nudging         [m/s]
    vm_hoc_grid       ! Observed wind, for nudging         [m/s]

! This code block takes the model time, finds the time before it and the time after it on
! the list, and marks them left_time and right_time for interpolation.  If the time is
! before the first or after the last time in the file, it just uses the first or last
! time without interpolation.

    left_time = -1
    right_time = -1

    if (time <= file_times(1)) then
      if ( clubb_at_least_debug_level( 1 ) ) then
        write(fstderr,*) 'Time is at or before the first time in the list.'
      endif
      left_time = 1
      right_time = 1
    else if (time >= file_times(file_ntimes)) then
      if ( clubb_at_least_debug_level( 1 ) ) then
        write(fstderr,*) 'Time is at or after the last time in the list.'
      endif
      left_time = file_ntimes
      right_time = file_ntimes
    else
      do k=1,file_ntimes-1
        if ((time > file_times(k)) .AND. & 
            (time <=file_times(k+1))) then
          left_time = k
          right_time = k+1
        end if
      end do
    end if

    if( left_time == -1 .or. right_time == -1 ) then
      call clubb_debug(1, "file_times not sorted in mpace_a_tndcy.")
    endif

! This is the ratio "a" needed for linear interpolation in time.
    ratio = real((time - file_times(left_time)) /  &          ! at the first time a=0;
            (file_times(right_time) - file_times(left_time))) ! at the second time a=1.

    do k=1,file_nlevels
!        omega_column(k) = ratio *			       ! Do linear interpolation in time
!     .                      (omega_forcing(k,right_time)
!     .                      -omega_forcing(k,left_time))
!     .                     + omega_forcing(k,left_time)

      dTdt_column(k)  = factor_interp( ratio, dTdt_forcing(k, right_time), &
                                       dTdt_forcing(k, left_time) )
      dqdt_column(k)  = factor_interp( ratio, dqdt_forcing(k, right_time), &
                                       dqdt_forcing(k, left_time) )
      vertT_column(k) = factor_interp( ratio, vertT_forcing(k,right_time), &
                                       vertT_forcing(k,left_time) )
      vertq_column(k) = factor_interp( ratio, vertq_forcing(k,right_time), &
                                       vertq_forcing(k,left_time) )
      um_column(k)    = factor_interp( ratio, um_obs(k, right_time), um_obs(k, left_time) )
      vm_column(k)    = factor_interp( ratio, vm_obs(k, right_time), vm_obs(k, left_time) )
    end do

!     Do linear interpolation in space
!     using zlinterp_fnc
    dTdt_hoc_grid  = zlinterp_fnc(gr%nnzp, file_nlevels, gr%zt, & 
                             file_heights,dTdt_column)
    dqdt_hoc_grid  = zlinterp_fnc(gr%nnzp, file_nlevels, gr%zt, & 
                             file_heights,dqdt_column)
    vertT_hoc_grid  = zlinterp_fnc(gr%nnzp, file_nlevels, gr%zt, & 
                             file_heights,vertT_column)
    vertq_hoc_grid  = zlinterp_fnc(gr%nnzp, file_nlevels, gr%zt, & 
                             file_heights,vertq_column)
    um_hoc_grid  = zlinterp_fnc(gr%nnzp, file_nlevels, gr%zt, & 
                             file_heights,um_column)
    vm_hoc_grid  = zlinterp_fnc(gr%nnzp, file_nlevels, gr%zt, & 
                             file_heights,vm_column)

    um_hoc_grid (1) = um_hoc_grid(2)
    vm_hoc_grid (1) = vm_hoc_grid(2)

! eMFc

!-----------------------------------------------------------------------


    ! Compute vertical motion
    do i=2,gr%nnzp
!          velocity_omega = omega_hoc_grid(i) * 100 / 3600 ! convering mb/hr to Pa/s
!          wm_zt(i) = -velocity_omega * Rd * thvm(i) / p_in_Pa(i) / grav
      wm_zt(i) = 0.
! End of Michael Falk's obliteration of omega.
    end do

    ! Boundary condition
    wm_zt(1) = 0.0        ! Below surface

    ! Interpolation
    wm_zm = zt2zm( wm_zt )

    ! Boundary condition
    wm_zm(1) = 0.0        ! At surface
    wm_zm(gr%nnzp) = 0.0  ! Model top


    ! Compute large-scale tendencies
    do i=1,gr%nnzp
      thlm_forcing(i) = ((dTdt_hoc_grid(i) + vertT_hoc_grid(i)) & 
                       * ((p_sfc/p_in_Pa(i)) ** (Rd/Cp))) & 
                       / 3600. ! K/s
      rtm_forcing(i)  = (dqdt_hoc_grid(i)+vertq_hoc_grid(i)) & 
       / 1000. / 3600. ! g/kg/hr -> kg/kg/s
    end do

    ! Test scalars with thetal and rt if desired
    if ( iisclr_thl > 0 ) sclrm_forcing(:,iisclr_thl) = thlm_forcing
    if ( iisclr_rt  > 0 ) sclrm_forcing(:,iisclr_rt)  = rtm_forcing

    if ( iiedsclr_thl > 0 ) edsclrm_forcing(:,iiedsclr_thl) = thlm_forcing
    if ( iiedsclr_rt  > 0 ) edsclrm_forcing(:,iiedsclr_rt)  = rtm_forcing

    return
  end subroutine mpace_a_tndcy

!----------------------------------------------------------------------
  subroutine mpace_a_sfclyr( time, rho0, & 
                             wpthlp_sfc, wprtp_sfc, ustar )
!        Description:
!          Surface forcing subroutine for mpace_a case.  Written
!          October 2007 by Michael Falk.
!
!        References:
!          mpace_a specification from arm.gov
!-----------------------------------------------------------------------

    use constants_clubb, only: Cp, Lv, fstderr ! Variable(s)

    use stats_precision, only: time_precision ! Variable(s)

    use error_code, only: clubb_debug, clubb_at_least_debug_level ! Procedure(s)

    use interpolation, only: factor_interp ! Procedure(s)

    implicit none

    ! External
    intrinsic :: max, sqrt, present

    ! Input Variables
    real(kind=time_precision), intent(in) :: & 
    time     ! current model time           [s]

    real, intent(in)  :: & 
    rho0     ! Air density at surface       [kg/m^3]

    ! Output Variables
    real, intent(out) ::  & 
    wpthlp_sfc,   & ! w'th_l' at (1)   [(m K)/s]  
    wprtp_sfc,    & ! w'r_t' at (1)    [(m kg)/(s kg)]
    ustar           ! surface friction velocity [m/s]

    ! Local Variables
    real :: & 
      latent_heat_flx, & 
      sensible_heat_flx

    integer :: k

    integer :: & 
      left_time, right_time

    real :: ratio
    !-----------------------------------------------------------------------

    left_time = -1
    right_time = -1

    ! choose which times to use
    if (time <= file_times(1)) then
      if ( clubb_at_least_debug_level( 1 ) ) then
        write(fstderr,*) 'Time is at or before the first time in the list.'
      endif
      left_time = 1
      right_time = 1
    else if (time >= file_times(file_ntimes)) then
      if ( clubb_at_least_debug_level( 1 ) ) then
        write(fstderr,*) 'Time is at or after the last time in the list.'
      endif
      left_time = file_ntimes
      right_time = file_ntimes
    else
      do k=1,file_ntimes-1
        if ((time > file_times(k)) .AND. & 
            (time <= file_times(k+1))) then
          left_time = k
          right_time = k+1
        end if
      end do
    end if

    ! Sanity check to make certain that the values read into
    ! file_times are sorted. Joshua Fasching June 2008
    if ( left_time == -1 .or. right_time == -1 ) then
      call clubb_debug(1, "file_times not sorted in MPACE_A")
    endif

    ratio = real(((time-file_times(left_time)) /  & 
         (file_times(right_time)-file_times(left_time))))

    latent_heat_flx = factor_interp( ratio, file_LH(right_time), file_LH(left_time) )

    sensible_heat_flx = factor_interp( ratio, file_SH(right_time), file_SH(left_time) )

    ! Compute heat and moisture fluxes
    wpthlp_sfc = sensible_heat_flx/(rho0*Cp)
    wprtp_sfc  = latent_heat_flx/(rho0*Lv)

    ! Declare the value of ustar.
    ustar = 0.25

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
      file_ntimes, per_line, file_LH )

    call file_read_1d( iunit, & 
      file_path//'mpace_a_sh.dat', & 
      file_ntimes, per_line, file_SH )

    return
  end subroutine mpace_a_init

end module mpace_a
