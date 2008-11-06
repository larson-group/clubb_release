!-----------------------------------------------------------------------
!  $Id$
module stats_subs

  implicit none
      
  private ! Set Default Scope
      
  public :: stats_init, stats_begin_timestep, stats_end_timestep, & 
    stats_accumulate, stats_finalize
      
  private :: stats_zero, stats_avg
      
  contains
      
!-----------------------------------------------------------------------
  subroutine stats_init( iunit, fname_prefix, l_stats_in, stats_fmt_in, stats_tsamp_in, &
                         stats_tout_in, fnamelist, nnzp, gzt, gzm, & 
                         day, month, year, rlat, rlon, time_current, delt )


!     Description: Initializes the statistics saving functionality of
!     the HOC model.
!-----------------------------------------------------------------------
    use stats_variables, only: & 
      zt,      & ! Variables
      ztscr01, & 
      ztscr02, & 
      ztscr03, & 
      ztscr04, & 
      ztscr05, & 
      ztscr06, & 
      ztscr07, & 
      ztscr08, & 
      ztscr09, & 
      ztscr10, & 
      ztscr11, & 
      ztscr12, & 
      ztscr13, & 
      ztscr14, & 
      ztscr15, & 
      ztscr16, & 
      ztscr17, & 
      ztscr18, & 
      ztscr19, & 
      ztscr20, & 
      ztscr21, & 
      zm, & 
      zmscr01, & 
      zmscr02, & 
      zmscr03, & 
      zmscr04, & 
      zmscr05, & 
      zmscr06, & 
      zmscr07, & 
      zmscr08, & 
      zmscr09, & 
      zmscr10, & 
      zmscr11, & 
      zmscr12, & 
      zmscr13, & 
      zmscr14, & 
      zmscr15, &
      zmscr16, &
      zmscr17, &
      sfc, & 
      l_stats, & 
      stats_tsamp, & 
      stats_tout, & 
      l_stats_samp, & 
      l_stats_first, & 
      l_stats_last, & 
      fname_zt, & 
      fname_zm, & 
      fname_sfc, & 
      l_netcdf, & 
      l_grads
    use stats_precision, only: & 
      time_precision   ! Variable(s)
    use output_grads, only: & 
      open_grads  ! Procedure
#ifdef NETCDF
    use output_netcdf, only: & 
      open_netcdf     ! Procedure
#endif
    use stats_zm, only: & 
      stats_init_zm ! Procedure
    use stats_zt, only: & 
      stats_init_zt ! Procedure
    use stats_sfc, only: & 
      stats_init_sfc ! Procedure

    use error_code, only: &
      clubb_at_least_debug_level ! Function

    use constants, only: &
      fstdout, fstderr ! Constants

    implicit none

    ! Constant Parameters
 
    integer, parameter :: nvarmax = 250  ! Max variables

    ! Input Variables

    integer, intent(in) :: iunit  ! File unit for fnamelist

    character(len=*), intent(in) ::  & 
      fname_prefix    ! Start of the stats filenames

    logical, intent(in) :: l_stats_in ! Stats on? T/F

    character(len=*), intent(in) :: &
      stats_fmt_in    ! Format of the stats file output

    real(kind=time_precision), intent(in) ::  & 
      stats_tsamp_in,  & ! Sampling interval   [s]
      stats_tout_in      ! Output interval     [s]

    character(len=*), intent(in) :: &
      fnamelist          ! Filename holding the &statsnl

    integer, intent(in) :: nnzp ! Grid points in the vertical [count]

    real, intent(in), dimension(nnzp) ::  & 
      gzt, gzm  ! Thermodynamic and momentum levels           [m]

    integer, intent(in) :: day, month, year  ! Time of year

    real, intent(in) ::  & 
      rlat, rlon   ! Latitude and Longitude             [Degrees N/E]

    real(kind=time_precision), intent(in) ::  & 
      time_current ! Model time                         [s]

    real(kind=time_precision), intent(in) ::  & 
      delt         ! Timestep (dtmain in HOC)           [s]


    ! Local Variables

    ! Namelist Variables

    character(len=10) :: stats_fmt  ! File storage convention

    character(len=30), dimension(nvarmax) ::  & 
      vars_zt,   & ! Variables on the thermodynamic levels
      vars_zm,   & ! Variables on the momentum levels
      vars_sfc  ! Variables at the model surface

    namelist /statsnl/ & 
      vars_zt, & 
      vars_zm, & 
      vars_sfc

    ! Local Variables

    logical :: l_error

    character(len=200) :: fdir, fname

    integer :: i, ntot

    ! Initialize
    l_error = .false.

    ! Set stats_variables variables with inputs from calling subroutine
    l_stats = l_stats_in

    stats_tsamp = stats_tsamp_in
    stats_tsamp = stats_tsamp_in
    stats_tout  = stats_tout_in
    stats_fmt   = trim( stats_fmt_in )

    if ( .not. l_stats ) then
      l_stats_samp  = .false.
      l_stats_first = .false.
      l_stats_last  = .false.
      return
    end if

    ! Initialize namelist variables

    vars_zt  = ''
    vars_zm  = ''
    vars_sfc = ''

      ! Read namelist

    open(unit=iunit, file=fnamelist)
    read(unit=iunit, nml=statsnl, end=100)
    close(unit=iunit)

    if ( clubb_at_least_debug_level( 1 ) ) then
      write(fstdout,*) "--------------------------------------------------"

      write(fstdout,*) "Statistics"

      write(fstdout,*) "--------------------------------------------------"
      write(fstdout,*) "vars_zt = "
      i = 1
      do while ( vars_zt(i) /= '' )
        write(fstdout,*) vars_zt(i)
        i = i + 1
      end do
 
      write(fstdout,*) "vars_zm = "
      i = 1
      do while ( vars_zm(i) /= '' )
        write(fstdout,*) vars_zm(i)
        i = i + 1
      end do

      write(fstdout,*) "vars_sfc = "
      i = 1
      do while ( vars_sfc(i) /= '' )
        write(fstdout,*) vars_sfc(i)
        i = i + 1
      end do

      write(fstdout,*) "--------------------------------------------------"
    end if ! clubb_at_least_debug_level 1

    ! Determine file names for GrADS or NetCDF files
    fname_zt  = trim( fname_prefix )//"_zt"
    fname_zm  = trim( fname_prefix )//"_zm"
    fname_sfc = trim( fname_prefix )//"_sfc"

    ! Parse the file type for stats output.  Currently only GrADS and
    ! NetCDF v3 are supported by this code.

    select case( trim( stats_fmt ) ) 
    case( "GrADS", "grads", "gr" )
      l_netcdf = .false.
      l_grads  = .true.

    case ( "NetCDF", "netcdf", "nc" )
      l_netcdf = .true.
      l_grads  = .false.

    case default
      write(fstderr,*) "Invalid data format "//trim( stats_fmt )
      stop

    end select

    ! Check sampling and output frequencies

    if ( abs( stats_tsamp/delt - floor(stats_tsamp/delt) )  & 
           > 1.e-8 ) then
      l_error = .true.
      write(fstderr,*) 'Error: stats_tsamp should be a multiple of delt'
      write(fstderr,*) 'stats_tsamp = ',stats_tsamp
      write(fstderr,*) 'delt = ',delt
    end if

    if ( abs( stats_tout/stats_tsamp - floor(stats_tout/stats_tsamp) ) & 
           > 1.e-8 ) then
      l_error = .true.
      write(0,*)  'Error: stats_tout should be a multiple of stats_tsamp'
      write(0,*) 'stats_tout = ',stats_tout
      write(0,*) 'stats_tsamp = ',stats_tsamp
    end if

    ! Initialize zt (mass points)

    i = 1
    do while ( ichar(vars_zt(i)(1:1)) /= 0  & 
               .and. len_trim(vars_zt(i)) /= 0 & 
               .and. i <= nvarmax )
      i = i + 1
    end do

    ntot = i - 1
    if ( ntot == nvarmax ) then
      write(fstderr,*) 'WARNING: check nvarmax in statistics.f'
    end if
      
    zt%nn = ntot
    zt%kk = nnzp
!      write(*,*) 'Number of variables for zt ',zt%nn

    allocate( zt%z( zt%kk ) )
    zt%z = gzt

    allocate( zt%x( zt%kk, zt%nn ) )
    allocate( zt%n( zt%kk, zt%nn ) )
    allocate( zt%l_in_update( zt%kk, zt%nn ) )
    call stats_zero( zt%kk, zt%nn, zt%x, zt%n, zt%l_in_update )

    allocate( zt%f%var( zt%nn ) )
    allocate( zt%f%z( zt%kk ) )

      ! Allocate scratch space

    allocate( ztscr01(zt%kk) )
    allocate( ztscr02(zt%kk) )
    allocate( ztscr03(zt%kk) )
    allocate( ztscr04(zt%kk) )
    allocate( ztscr05(zt%kk) )
    allocate( ztscr06(zt%kk) )
    allocate( ztscr07(zt%kk) )
    allocate( ztscr08(zt%kk) )
    allocate( ztscr09(zt%kk) )
    allocate( ztscr10(zt%kk) )
    allocate( ztscr11(zt%kk) )
    allocate( ztscr12(zt%kk) )
    allocate( ztscr13(zt%kk) )
    allocate( ztscr14(zt%kk) )
    allocate( ztscr15(zt%kk) )
    allocate( ztscr16(zt%kk) )
    allocate( ztscr17(zt%kk) )
    allocate( ztscr18(zt%kk) )
    allocate( ztscr19(zt%kk) )
    allocate( ztscr20(zt%kk) )
    allocate( ztscr21(zt%kk) )

    ztscr01 = 0.0
    ztscr02 = 0.0
    ztscr03 = 0.0
    ztscr04 = 0.0
    ztscr05 = 0.0
    ztscr06 = 0.0
    ztscr07 = 0.0
    ztscr08 = 0.0
    ztscr09 = 0.0
    ztscr10 = 0.0
    ztscr11 = 0.0
    ztscr12 = 0.0
    ztscr13 = 0.0
    ztscr14 = 0.0
    ztscr15 = 0.0
    ztscr16 = 0.0
    ztscr17 = 0.0
    ztscr18 = 0.0
    ztscr19 = 0.0
    ztscr20 = 0.0
    ztscr21 = 0.0

    fdir = "./"
    fname = trim( fname_zt )

    if ( l_grads ) then

        ! Open GrADS file
      call open_grads( iunit, fdir, fname,  & 
                       1, zt%kk, zt%z, & 
                       day, month, year, rlat, rlon, & 
                       time_current+stats_tout, stats_tout, & 
                       zt%nn,zt%f )

    else ! Open NetCDF file
#ifdef NETCDF
      call open_netcdf( iunit, fdir, fname,  & 
                        1, zt%kk, zt%z, & 
                        day, month, year, rlat, rlon, & 
                        time_current+stats_tout, stats_tout, & 
                        zt%nn, zt%f )
#else
      stop "netCDF support was not compiled into this build."
#endif

    end if

    ! Default initialization for array indices for zt

    call stats_init_zt( vars_zt, l_error )

    ! Initialize zm (momentum points)

    i = 1
    do while ( ichar(vars_zm(i)(1:1)) /= 0  & 
               .and. len_trim(vars_zm(i)) /= 0 & 
               .and. i <= nvarmax )
      i = i + 1
    end do
    ntot = i - 1
    if ( ntot == nvarmax ) write(fstderr,*) 'WARNING: check nvarmax in statistics.f'

    zm%nn = ntot
    zm%kk = nnzp
!   write(*,*) 'Number of variables for zm ',zm%nn

    allocate( zm%z( zm%kk ) )
    zm%z = gzm

    allocate( zm%x( zm%kk, zm%nn ) )
    allocate( zm%n( zm%kk, zm%nn ) )
    allocate( zm%l_in_update( zm%kk, zm%nn ) )
      
    call stats_zero( zm%kk, zm%nn, zm%x, zm%n, zm%l_in_update )

    allocate( zm%f%var( zm%nn ) )
    allocate( zm%f%z( zm%kk ) )

      ! Allocate scratch space

    allocate( zmscr01(zm%kk) )
    allocate( zmscr02(zm%kk) )
    allocate( zmscr03(zm%kk) )
    allocate( zmscr04(zm%kk) )
    allocate( zmscr05(zm%kk) )
    allocate( zmscr06(zm%kk) )
    allocate( zmscr07(zm%kk) )
    allocate( zmscr08(zm%kk) )
    allocate( zmscr09(zm%kk) )
    allocate( zmscr10(zm%kk) )
    allocate( zmscr11(zm%kk) )
    allocate( zmscr12(zm%kk) )
    allocate( zmscr13(zm%kk) )
    allocate( zmscr14(zm%kk) )
    allocate( zmscr15(zm%kk) )
    allocate( zmscr16(zm%kk) )
    allocate( zmscr17(zm%kk) )

    zmscr01 = 0.0
    zmscr02 = 0.0
    zmscr03 = 0.0
    zmscr04 = 0.0
    zmscr05 = 0.0
    zmscr06 = 0.0
    zmscr07 = 0.0
    zmscr08 = 0.0
    zmscr09 = 0.0
    zmscr10 = 0.0
    zmscr11 = 0.0
    zmscr12 = 0.0
    zmscr13 = 0.0
    zmscr14 = 0.0
    zmscr15 = 0.0
    zmscr16 = 0.0
    zmscr17 = 0.0


    fdir = "./"
    fname = trim( fname_zm )
    if ( l_grads ) then

      ! Open GrADS files
      call open_grads( iunit, fdir, fname,  & 
                       1, zm%kk, zm%z, & 
                       day, month, year, rlat, rlon, & 
                       time_current+stats_tout, stats_tout, & 
                       zm%nn, zm%f )

    else ! Open NetCDF file
#ifdef NETCDF
      call open_netcdf( iunit, fdir, fname,  & 
                        1, zm%kk, zm%z, & 
                        day, month, year, rlat, rlon, & 
                        time_current+stats_tout, stats_tout, & 
                        zm%nn, zm%f )

#else
      stop "netCDF support was not compiled into this build."
#endif
    end if

      call stats_init_zm( vars_zm, l_error )

      ! Initialize sfc (surface point)

    i = 1
    do while ( ichar(vars_sfc(i)(1:1)) /= 0  & 
               .and. len_trim(vars_sfc(i)) /= 0 & 
               .and. i <= nvarmax )
      i = i + 1
    end do

    ntot = i - 1

    if ( ntot == nvarmax ) write(fstderr,*) 'WARNING: check nvarmax in statistics.f'

    sfc%nn = ntot
    sfc%kk = 1
!   write(*,*) 'Number of variables for sfc ',sfc%nn

    allocate( sfc%z( sfc%kk ) )
    sfc%z = gzm(1)

    allocate( sfc%x( sfc%kk, sfc%nn ) )
    allocate( sfc%n( sfc%kk, sfc%nn ) )
    allocate( sfc%l_in_update( sfc%kk, sfc%nn ) )
      
    call stats_zero( sfc%kk, sfc%nn, sfc%x, sfc%n, sfc%l_in_update )

    allocate( sfc%f%var( sfc%nn ) )
    allocate( sfc%f%z( sfc%kk ) )

    fdir = "./"
    fname = trim( fname_sfc )

    if ( l_grads ) then

        ! Open GrADS files
      call open_grads( iunit, fdir, fname,  & 
                       1, sfc%kk, sfc%z, & 
                       day, month, year, rlat, rlon, & 
                         time_current+stats_tout, stats_tout, & 
                         sfc%nn, sfc%f )

    else ! Open NetCDF files
#ifdef NETCDF
      call open_netcdf( iunit, fdir, fname,  & 
                        1, sfc%kk, sfc%z, & 
                        day, month, year, rlat, rlon, & 
                        time_current+stats_tout, stats_tout, & 
                        sfc%nn, sfc%f )

#else
      stop "netCDF support was not compiled into this build."
#endif
    end if

    call stats_init_sfc( vars_sfc, l_error )

    ! Check for errors

    if ( l_error ) then
      write(fstderr,*) 'stats_init: errors found'
      stop
    end if

    return

    ! If namelist was not found in input file, turn off statistics

100 continue
    write(fstderr,*) 'Error with statsnl, statistics is turned off'
    l_stats       = .false.
    l_stats_samp  = .false.
    l_stats_first = .false.
    l_stats_last  = .false.

    return
  end subroutine stats_init
!-----------------------------------------------------------------------
      subroutine stats_zero( kk, nn, x, n, l_in_update )

!     Description:
!     Initialize stats to zero
!-----------------------------------------------------------------------
      use stats_precision, only: & 
          stat_rknd,   & ! Variable(s)
          stat_nknd

      implicit none

      ! Input
      integer, intent(in) :: kk, nn

      ! Output
      real(kind=stat_rknd), dimension(kk,nn), intent(out)    :: x
      integer(kind=stat_nknd), dimension(kk,nn), intent(out) :: n
      logical, dimension(kk,nn), intent(out) :: l_in_update

      ! Zero out arrays

      if ( nn > 0 ) then
        x(:,:) = 0.0
        n(:,:) = 0
        l_in_update(:,:) = .false.
      end if
      
      return
      end subroutine stats_zero

!-----------------------------------------------------------------------
      subroutine stats_avg( kk, nn, x, n )

!     Description:
!     Compute the average of stats fields
!-----------------------------------------------------------------------
      use stats_precision, only: & 
          stat_rknd,   & ! Variable(s)
          stat_nknd

      implicit none

      ! Input
      integer, intent(in) :: nn, kk
      integer(kind=stat_nknd), dimension(kk,nn), intent(in) :: n

      ! Output
      real(kind=stat_rknd), dimension(kk,nn), intent(inout)  :: x

      ! Internal

      integer k,m

      ! Compute averages

      do m=1,nn
        do k=1,kk

          if ( n(k,m) > 0 ) then
            x(k,m) = x(k,m) / n(k,m)
          end if

        end do
      end do

      return
      end subroutine stats_avg

!-----------------------------------------------------------------------
      subroutine stats_begin_timestep( time_elapsed, delt )

!     Description:
!     Begin sampling for the current timestep.
!-----------------------------------------------------------------------

      use stats_variables, only: & 
          l_stats,  & ! Variable(s)
          l_stats_samp, & 
          l_stats_first, & 
          l_stats_last, & 
          stats_tsamp, & 
          stats_tout
      use stats_precision, only: & 
          time_precision ! Variable(s)

      implicit none

      ! Input

      real(kind=time_precision), intent(in) ::  & 
        time_elapsed ! Elapsed model time       [s]

      real(kind=time_precision), intent(in) ::  & 
        delt         ! Model time step          [s]

      if ( .not. l_stats ) return

      ! Set sample this time step flag
      if ( mod( time_elapsed, stats_tsamp ) < 1.e-8 ) then
        l_stats_samp = .true.
      else
        l_stats_samp = .false.
      end if

      ! Set first time step flag

      if ( mod( time_elapsed - delt, stats_tout ) < 1.e-8 ) then
        l_stats_first = .true.
      else
        l_stats_first = .false.
      end if

      ! Set last time step flag

      if ( mod( time_elapsed, stats_tout ) < 1.e-8 ) then
        l_stats_last = .true.
      else
        l_stats_last = .false.
      end if

      return

      end subroutine stats_begin_timestep

!-----------------------------------------------------------------------
      subroutine stats_end_timestep( )

!     Description:
!-----------------------------------------------------------------------

      use stats_variables, only: & 
          zt,  & ! Variable(s)
          zm, & 
          sfc, & 
          l_stats_last, & 
          stats_tsamp, & 
          stats_tout, & 
          l_grads

      use stats_precision, only: & 
          time_precision ! Variable(s)

      use output_grads, only: & 
          write_grads ! Procedure(s)

#ifdef NETCDF
      use output_netcdf, only: & 
          write_netcdf ! Procedure(s)
#endif

      implicit none

      ! Local Variables

      integer :: i, k

      ! Check if it is time to write to file

      if ( .not. l_stats_last ) return

      ! Check number of sampling points

      do i=1,zt%nn
       do k=1,zt%kk
         if ( zt%n(k,i) /= 0  & 
              .and. zt%n(k,i) /= floor(stats_tout/stats_tsamp) ) then
           write(0,*) 'Possible sampling error for variable ', & 
                           trim(zt%f%var(i)%name),' in zt ','at k =',k, & 
                           ' zt%n(',k,',',i,')=',zt%n(k,i)
!           write(0,*) 'Possible sampling error for zt ',i,k,zt%n(k,i)
           !pause
         end if
        end do
      end do
      
      do i=1,zm%nn
       do k=1,zm%kk
         if ( zm%n(k,i) /= 0  & 
              .and. zm%n(k,i) /= floor(stats_tout/stats_tsamp) ) then
           write(0,*) 'Possible sampling error for variable ', & 
                           trim(zm%f%var(i)%name),' in zm ','at k =',k, & 
                           ' zm%n(',k,',',i,')=',zm%n(k,i)
!           write(0,*) 'Possible sampling error for zm ',i,k,zm%n(k,i)
!           Made error message more descriptive
!           Joshua Fasching July 2008

           !pause
         end if
        end do
      end do
      
      do i=1,sfc%nn
       do k=1,sfc%kk
         if ( sfc%n(k,i) /= 0  & 
              .and. sfc%n(k,i) /= floor(stats_tout/stats_tsamp) ) then
           write(0,*) 'Possible sampling error for variable ', & 
                          trim(sfc%f%var(i)%name),' in sfc ','at k =',k, & 
                           ' sfc%n(',k,',',i,')=',sfc%n(k,i)
!           write(0,*) 'Possible sampling error for sfc ',i,k,sfc%n(k,i)
!           Made error message more descriptive
!           Joshua Fasching July 2008

           !pause
         end if
        end do
      end do
      
      ! Compute averages

      call stats_avg( zt%kk, zt%nn, zt%x, zt%n )
      call stats_avg( zm%kk, zm%nn, zm%x, zm%n )
      call stats_avg( sfc%kk, sfc%nn, sfc%x, sfc%n )

      ! Write to file
      if ( l_grads ) then
        call write_grads( zt%f  )
        call write_grads( zm%f  )
        call write_grads( sfc%f  )
      else ! l_netcdf
#ifdef NETCDF
        call write_netcdf( zt%f  )
        call write_netcdf( zm%f  )
        call write_netcdf( sfc%f  )
#else 
        stop "This program was not compiled with netCDF support"
#endif
      end if

      ! Reset sample fields
      call stats_zero( zt%kk, zt%nn, zt%x, zt%n, zt%l_in_update )
      call stats_zero( zm%kk, zm%nn, zm%x, zm%n, zm%l_in_update )
      call stats_zero( sfc%kk, sfc%nn, sfc%x, sfc%n, sfc%l_in_update )


      return
      end subroutine stats_end_timestep

!----------------------------------------------------------------------
subroutine stats_accumulate & 
                 ( um, vm, upwp, vpwp, up2, vp2, thlm, & 
                   rtm, wprtp, wpthlp, wp2, wp3, rtp2, thlp2, rtpthlp, & 
                   p_in_Pa, exner, rho, rho_zm, & 
                   wm_zt, sigma_sqd_w, tau_zm, rcm, cf, & 
                   sclrm, edsclrm, sclrm_forcing, wpsclrp )

! Description:
! Accumulate those stats variables that are preserved in HOC from timestep to 
! timestep, but not those stats that are not, (e.g. budget terms, longwave and 
! shortwave components, etc. )
!----------------------------------------------------------------------

use stats_variables, only: & 
    zt,      & ! Variables
    zm, & 
    sfc, & 
    l_stats_samp, & 
    ithlm, & 
    iT_in_K, & 
    ithvm, & 
    irtm, & 
    ircm, & 
    ium, & 
    ivm, & 
    iwm_zt, & 
    iug, & 
    ivg, & 
    icf, & 
    ip_in_Pa, & 
    iexner, & 
    iLscale, & 
    iwp3, & 
    iwpthlp2, & 
    iwp2thlp,  & 
    iwprtp2, & 
    iwp2rtp, & 
    iLscale_up, & 
    iLscale_down, & 
    itau_zt, & 
    iKh_zt, & 
    iwp2thvp, & 
    iwp2rcp, & 
    iwprtpthlp, & 
    isigma_sqd_w_zt,          & 
    irho, & 
    irsat, & 
    iAKm, & 
    iAKm_est, & 
    iradht, & 
    ia, & 
    iw1, & 
    iw2, & 
    isw1, & 
    isw2, & 
    ithl1, & 
    ithl2, & 
    isthl1, & 
    isthl2, & 
    irt1, & 
    irt2, & 
    isrt1, & 
    isrt2, & 
    irc1, & 
    irc2, & 
    irsl1, & 
    irsl2, & 
    iR1, & 
    iR2, & 
    is1, & 
    is2, & 
    iss1, & 
    iss2, & 
    irrtthl

use stats_variables, only: & 
    iwp2_zt, & 
    ithlp2_zt, & 
    iwpthlp_zt, & 
    iwprtp_zt, & 
    irtp2_zt, & 
    irtpthlp_zt, & 
    iwp2, & 
    irtp2, & 
    ithlp2, & 
    irtpthlp, & 
    iwprtp,  & 
    iwpthlp, & 
    iwp4,  & 
    iwpthvp, & 
    irtpthvp, & 
    ithlpthvp, & 
    itau_zm, & 
    iKh_zm, & 
    iwprcp, & 
    ithlprcp, & 
    irtprcp, & 
    ircp2, & 
    iupwp, & 
    ivpwp, & 
    iup2, & 
    ivp2, & 
    irho_zm, & 
    isigma_sqd_w, & 
    iem, & 
    ishear, & 
    iFrad, & 
    icc, & 
    izb, & 
    ilwp, &
    ithlm_vert_avg, &
    irtm_vert_avg, &
    ium_vert_avg, &
    ivm_vert_avg

use stats_variables, only: & 
    isclram, & 
    isclram_f, & 
    isclrbm, & 
    isclrbm_f, & 
    iedsclram, & 
    iedsclrbm, & 
    isclraprtp, & 
    isclrbprtp, & 
    isclrap2, & 
    isclrbp2, & 
    isclrapthvp, & 
    isclrbpthvp, & 
    isclrapthlp, & 
    isclrbpthlp, & 
    isclraprcp, & 
    isclrbprcp, & 
    iwpsclrap, & 
    iwpsclrbp, & 
    iwp2sclrap, & 
    iwp2sclrbp, & 
    iwpsclrap2, & 
    iwpsclrbp2, & 
    iwpsclraprtp, & 
    iwpsclrbprtp, & 
    iwpsclrapthlp, & 
    iwpsclrbpthlp, & 
    iwpedsclrap, & 
    iwpedsclrbp

use grid_class, only: & 
    gr ! Variable

use variables_diagnostic_module, only: & 
    pdf_parms,  & ! Variable(s)
    thvm, & 
    ug, & 
    vg, & 
    Lscale, & 
    wpthlp2, & 
    wp2thlp, & 
    wprtp2, & 
    wp2rtp, & 
    Lscale_up, & 
    Lscale_down, & 
    tau_zt, & 
    Kh_zt, & 
    wp2thvp, & 
    wp2rcp, & 
    wprtpthlp, & 
    sigma_sqd_w_zt, & 
    rsat, & 
    Akm, & 
    Akm_est, & 
    radht, & 
    wp2_zt, & 
    thlp2_zt, & 
    wpthlp_zt, & 
    wprtp_zt, & 
    rtp2_zt, & 
    rtpthlp_zt, & 
    wp4, & 
    wpthvp, & 
    rtpthvp, & 
    thlpthvp, & 
    Kh_zm, & 
    wprcp, & 
    thlprcp, & 
    rtprcp, & 
    rcp2, & 
    em, & 
    shear, & 
    Frad, & 
    sclrprtp, & 
    sclrp2, & 
    sclrpthvp, & 
    sclrpthlp, & 
    sclrprcp, & 
    wp2sclrp, & 
    wpsclrp2, & 
    wpsclrprtp, & 
    wpsclrpthlp, & 
    wpedsclrp   
    
use model_flags, only: & 
    l_LH_on ! Variable(s)

use T_in_K_mod, only: & 
    thlm2T_in_K ! Procedure

use constants, only: & 
    rc_tol

use parameters_tunable, only: & 
    sclr_dim  ! Variable(s)

use stats_type, only: & 
    stat_update_var,  & ! Procedure(s)
    stat_update_var_pt

use fill_holes, only: &
    vertical_avg

use interpolation, only: & 
    lin_int ! Procedure

implicit none

! Input Variable
real, intent(in), dimension(gr%nnzp) :: & 
  um,      & ! u wind                        [m/s]
  vm,      & ! v wind                        [m/s]
  upwp,    & ! vertical u momentum flux      [m^2/s^2]
  vpwp,    & ! vertical v momentum flux      [m^2/s^2]
  up2,     & ! u'^2                          [m^2/s^2]
  vp2,     & ! v'^2                          [m^2/s^2]
  thlm,    & ! liquid potential temperature  [K]
  rtm,     & ! total water mixing ratio      [kg/kg]
  wprtp,   & ! w'rt'                         [m kg/s kg]
  wpthlp,  & ! w'thl'                        [m K /s]
  wp2,     & ! w'^2                          [m^2/s^2]
  wp3,     & ! w'^3                          [m^3/s^3]
  rtp2,    & ! rt'^2                         [kg/kg]
  thlp2,   & ! thl'^2                        [K^2]
  rtpthlp    ! rt'thl'                       [kg/kg K]

real, intent(in), dimension(gr%nnzp) :: & 
  p_in_Pa,      & ! Pressure (Pa) on thermodynamic points    [Pa]
  exner,        & ! Exner function = ( p / p0 ) ** kappa     [-]
  rho,          & ! Density                                  [kg/m^3]
  rho_zm,       & ! Density                                  [kg/m^3]
  wm_zt,        & ! w on thermodynamic levels                [m/s]
  sigma_sqd_w,  & ! PDF width paramter                       [-]
  tau_zm          ! Dissipation time                         [s]

real, intent(in), dimension(gr%nnzp) :: & 
  rcm,   & ! Cloud water mixing ratio                [kg/kg]
!  Ncm,   & ! Cloud droplet number concentration      [num/kg]
!  Ncnm,  & ! Cloud nuclei number concentration       [num/m^3]
!  Nim,   & ! Ice nuclei number concentration         [num/m^3]
  cf       ! Cloud fraction                          [%]

real, intent(in), dimension(gr%nnzp,sclr_dim) :: & 
  sclrm,           & ! High-Order Passive scalar     [units vary]
  edsclrm,         & ! Eddy-diff Passive scalar      [units vary] 
  sclrm_forcing,   & ! Large-scale forcing of scalar [units/s]
  wpsclrp         ! w'sclr'                       [units m/s]

! Prognostic drizzle variable array
!real, intent(in), dimension(gr%nnzp,hydromet_dim) :: hydromet
! Contains:
! 1 rrainm   Rain water mixing ratio               [kg/kg]
! 2 Nrm      Rain droplet number concentration     [num/kg]
! 3 rsnow    Snow water mixing ratio               [kg/kg]
! 4 rice     Ice water mixing ratio                [kg/kg]
! 5 rgraupel Graupel water mixing ratio            [kg/kg]

! Local Variables

integer :: i, k

real :: xtmp

! Sample fields

if ( l_stats_samp ) then

   ! zt variables

   call stat_update_var( ithlm, thlm, zt )
   call stat_update_var( iT_in_K,  & 
                         thlm2T_in_K( thlm, exner, rcm), zt )
   call stat_update_var( ithvm, thvm, zt )
   call stat_update_var( irtm, rtm, zt )       
   call stat_update_var( ircm, rcm, zt )
   call stat_update_var( ium, um, zt )
   call stat_update_var( ivm, vm, zt )
   call stat_update_var( iwm_zt, wm_zt, zt )
   call stat_update_var( iug, ug, zt )
   call stat_update_var( ivg, vg, zt )
   call stat_update_var( icf, cf, zt )
   call stat_update_var( ip_in_Pa, p_in_Pa, zt )
   call stat_update_var( iexner, exner, zt )
   call stat_update_var( iLscale, Lscale, zt )
   call stat_update_var( iwp3, wp3, zt )
   call stat_update_var( iwpthlp2, wpthlp2, zt )
   call stat_update_var( iwp2thlp, wp2thlp, zt )
   call stat_update_var( iwprtp2, wprtp2, zt )
   call stat_update_var( iwp2rtp, wp2rtp, zt )
   call stat_update_var( iLscale_up, Lscale_up, zt )
   call stat_update_var( iLscale_down, Lscale_down, zt )
   call stat_update_var( itau_zt, tau_zt, zt )
   call stat_update_var( iKh_zt, Kh_zt, zt )
   call stat_update_var( iwp2thvp, wp2thvp, zt )
   call stat_update_var( iwp2rcp, wp2rcp, zt )
   call stat_update_var( iwprtpthlp, wprtpthlp, zt )
   call stat_update_var( isigma_sqd_w_zt, sigma_sqd_w_zt, zt )
   call stat_update_var( irho, rho, zt )
!   call stat_update_var( iNcm, Ncm, zt )
!   call stat_update_var( iNcnm, Ncnm, zt )
!   call stat_update_var( iNim, Nim, zt )
!   if ( l_cloud_sed ) then
!      call stat_update_var( ised_rcm, sed_rcm, zt )
!   endif
   call stat_update_var( irsat, rsat, zt )
!   call stat_update_var( irrainm, hydromet(:,1), zt )
!   call stat_update_var( iNrm, hydromet(:,2), zt )
!   call stat_update_var( irsnowm, hydromet(:,3), zt )
!   call stat_update_var( iricem, hydromet(:,4), zt )
!   call stat_update_var( irgraupelm, hydromet(:,5), zt )

   if ( l_LH_on ) then
      call stat_update_var( iAKm, AKm, zt )
      call stat_update_var( iAkm_est, AKm_est, zt)
   endif

   call stat_update_var( iradht, radht, zt )
   call stat_update_var( ia, pdf_parms(:,13), zt )
   call stat_update_var( iw1, pdf_parms(:,1), zt )
   call stat_update_var( iw2, pdf_parms(:,2), zt )
   call stat_update_var( isw1, pdf_parms(:,3), zt )
   call stat_update_var( isw2, pdf_parms(:,4), zt )
   call stat_update_var( ithl1, pdf_parms(:,9), zt )
   call stat_update_var( ithl2, pdf_parms(:,10), zt )
   call stat_update_var( isthl1, pdf_parms(:,11), zt )
   call stat_update_var( isthl2, pdf_parms(:,12), zt )
   call stat_update_var( irt1, pdf_parms(:,5), zt )
   call stat_update_var( irt2, pdf_parms(:,6), zt )
   call stat_update_var( isrt1, pdf_parms(:,7), zt )
   call stat_update_var( isrt2, pdf_parms(:,8), zt )
   call stat_update_var( irc1, pdf_parms(:,14), zt )
   call stat_update_var( irc2, pdf_parms(:,15), zt )
   call stat_update_var( irsl1, pdf_parms(:,16), zt )
   call stat_update_var( irsl2, pdf_parms(:,17), zt )
   call stat_update_var( iR1, pdf_parms(:,18), zt )
   call stat_update_var( iR2, pdf_parms(:,19), zt )
   call stat_update_var( is1, pdf_parms(:,20), zt )
   call stat_update_var( is2, pdf_parms(:,21), zt )
   call stat_update_var( iss1, pdf_parms(:,22), zt )
   call stat_update_var( iss2, pdf_parms(:,23), zt )
   call stat_update_var( irrtthl, pdf_parms(:,24), zt )
   call stat_update_var( iwp2_zt, wp2_zt, zt )
   call stat_update_var( ithlp2_zt, thlp2_zt, zt )
   call stat_update_var( iwpthlp_zt, wpthlp_zt, zt )
   call stat_update_var( iwprtp_zt, wprtp_zt, zt )
   call stat_update_var( irtp2_zt, rtp2_zt, zt )
   call stat_update_var( irtpthlp_zt, rtpthlp_zt, zt )

   if ( sclr_dim > 0 ) then
      call stat_update_var( isclram, sclrm(:,1), zt )
      call stat_update_var( isclram_f, sclrm_forcing(:,1),  zt )
      call stat_update_var( iedsclram, edsclrm(:,1), zt )
   endif

   if ( sclr_dim > 1 ) then
      call stat_update_var( isclrbm, sclrm(:,2), zt )
      call stat_update_var( isclrbm_f, sclrm_forcing(:,2), zt )
      call stat_update_var( iedsclrbm, edsclrm(:,2), zt )
   endif


   ! zm variables

   call stat_update_var( iwp2, wp2, zm )
   call stat_update_var( irtp2, rtp2, zm )
   call stat_update_var( ithlp2, thlp2, zm )
   call stat_update_var( irtpthlp, rtpthlp, zm )
   call stat_update_var( iwprtp, wprtp, zm )
   call stat_update_var( iwpthlp, wpthlp, zm )
   call stat_update_var( iwp4, wp4, zm )
   call stat_update_var( iwpthvp, wpthvp, zm )
   call stat_update_var( irtpthvp, rtpthvp, zm )
   call stat_update_var( ithlpthvp, thlpthvp, zm )
   call stat_update_var( itau_zm, tau_zm, zm )
   call stat_update_var( iKh_zm, Kh_zm, zm )
   call stat_update_var( iwprcp, wprcp, zm )
   call stat_update_var( ithlprcp, thlprcp, zm )
   call stat_update_var( irtprcp, rtprcp, zm )
   call stat_update_var( ircp2, rcp2, zm )
   call stat_update_var( iupwp, upwp, zm )
   call stat_update_var( ivpwp, vpwp, zm )
   call stat_update_var( ivp2, vp2, zm )
   call stat_update_var( iup2, up2, zm )
   call stat_update_var( irho_zm, rho_zm, zm )
   call stat_update_var( isigma_sqd_w, sigma_sqd_w, zm )
   call stat_update_var( iem, em, zm )
   call stat_update_var( ishear, shear, zm )
   call stat_update_var( iFrad, Frad, zm )
!   if ( l_cloud_sed ) then
!      call stat_update_var( iFcsed, Fcsed, zm )
!   endif

   if ( sclr_dim > 0 ) then
      call stat_update_var( isclraprtp, sclrprtp(:,1), zm )
      call stat_update_var( isclrap2, sclrp2(:,1), zm )
      call stat_update_var( isclrapthvp, sclrpthvp(:,1), zm )
      call stat_update_var( isclrapthlp, sclrpthlp(:,1), zm )
      call stat_update_var( isclraprcp, sclrprcp(:,1), zm ) 
      call stat_update_var( iwpsclrap, wpsclrp(:,1), zm )
      call stat_update_var( iwp2sclrap, wp2sclrp(:,1), zm )
      call stat_update_var( iwpsclrap2, wpsclrp2(:,1), zm )
      call stat_update_var( iwpsclraprtp, wpsclrprtp(:,1), zm )
      call stat_update_var( iwpsclrapthlp, wpsclrpthlp(:,1), zm )
      call stat_update_var( iwpedsclrap, wpedsclrp(:,1), zm )
   endif 

   if ( sclr_dim > 1 ) then
      call stat_update_var( isclrbprtp, sclrprtp(:,2), zm )
      call stat_update_var( isclrbp2, sclrp2(:,2), zm )
      call stat_update_var( isclrbpthvp, sclrpthvp(:,2), zm )
      call stat_update_var( isclrbpthlp, sclrpthlp(:,2), zm )
      call stat_update_var( isclrbprcp, sclrprcp(:,2), zm )
      call stat_update_var( iwpsclrbp, wpsclrp(:,2), zm )
      call stat_update_var( iwp2sclrbp, wp2sclrp(:,2), zm )
      call stat_update_var( iwpsclrbp2, wpsclrp2(:,2), zm )
      call stat_update_var( iwpsclrbprtp, wpsclrprtp(:,2), zm )
      call stat_update_var( iwpsclrbpthlp, wpsclrpthlp(:,2), zm )
      call stat_update_var( iwpedsclrbp, wpedsclrp(:,2), zm )
   endif 
        

   ! sfc variables

   ! Cloud cover
   call stat_update_var_pt( icc, 1, maxval( cf(1:gr%nnzp) ), sfc )

   ! Cloud base
   if ( izb > 0 ) then

      k = 1
      do while ( rcm(k) < rc_tol .and. k < gr%nnzp )
         k = k + 1
      enddo

      if ( k > 1 .AND. k < gr%nnzp) then

         ! Use linear interpolation to find the exact height of the 
         ! rc_tol kg/kg level.  Brian.
         call stat_update_var_pt( izb, 1, lin_int( rc_tol, rcm(k),  &
                                  rcm(k-1), gr%zt(k), gr%zt(k-1) ), sfc )

      else

         ! Mark the cloud base at -10 m. if it's clear.
         call stat_update_var_pt( izb, 1, -10.0 , sfc )

      endif

   endif

   ! Liquid Water Path
   if ( ilwp > 0 ) then

      xtmp = 0.
      do i = gr%nnzp-1, 1, -1
         xtmp = xtmp + rho(i+1) * rcm(i+1) / gr%dzt(i+1)
      enddo
          
      call stat_update_var_pt( ilwp, 1, xtmp, sfc )

   endif

   ! Vertical average of thermodynamic level variables.

   ! Find the vertical average of thermodynamic level variables, averaged from 
   ! level 2 (the first thermodynamic level above model surface) through 
   ! level gr%nnzp (the top of the model).  Use the vertical averaging function
   ! found in fill_holes.F90.

   ! Vertical average of thlm.
   call stat_update_var_pt( ithlm_vert_avg, 1,  &
        vertical_avg( 2, gr%nnzp, "zt", thlm(2:gr%nnzp) ), sfc )

   ! Vertical average of rtm.
   call stat_update_var_pt( irtm_vert_avg, 1,  &
        vertical_avg( 2, gr%nnzp, "zt", rtm(2:gr%nnzp) ), sfc )

   ! Vertical average of um.
   call stat_update_var_pt( ium_vert_avg, 1,  &
        vertical_avg( 2, gr%nnzp, "zt", um(2:gr%nnzp) ), sfc )

   ! Vertical average of vm.
   call stat_update_var_pt( ivm_vert_avg, 1,  &
        vertical_avg( 2, gr%nnzp, "zt", vm(2:gr%nnzp) ), sfc )

   ! Note:  currently, a vertical average cannot be taken properly for 
   !        hydrometeor (microphysical) variables due to the fact that function
   !        vertical_avg does not include the k=1 level, which is below the 
   !        model surface for thermodynamic level variables.  In general, this 
   !        is desired, as hole-filling should not consider the k=1 level.  
   !        However, hydrometeor variables do include level 1 in eddy-diffusion
   !        effects.  This should be included in a vertical average for 
   !        statistical output purposes.

   ! Note:  currently, a vertical average cannot be taken properly for momentum 
   !        level variables due to the fact that function vertical_avg does not 
   !        include the k=1 level, which is at the model surface for momentum 
   !        level variables.  In general, this is desired, as hole-filling 
   !        should not consider the k=1 level.  Momentum level variables have 
   !        set surface values that should not be changed due to hole-filling.
   !        Likewise, function vertical_avg does not include the k=gr%nnzp level
   !        for momentum level variables, as momentum level variables have set 
   !        values at the upper boundary.  However, values at the k=1 and 
   !        k=gr%nnzp levels need to be included in a vertical average for 
   !        statistical output purposes.

endif  ! l_stats_samp


return
end subroutine stats_accumulate

!-----------------------------------------------------------------------
      subroutine stats_finalize( )

!     Description:
!     Close NetCDF files and deallocate scratch space and 
!     stats file structures.
!-----------------------------------------------------------------------

      use stats_variables, only: & 
          zt,  & ! Variable(s)
          zm, & 
          sfc, & 
          ztscr01, & 
          ztscr02, & 
          ztscr03, & 
          ztscr04, & 
          ztscr05, & 
          ztscr06, & 
          ztscr07, & 
          ztscr08, & 
          ztscr09, & 
          ztscr10, & 
          ztscr11, & 
          ztscr12, & 
          ztscr13, & 
          ztscr14, & 
          ztscr15, & 
          ztscr16, & 
          ztscr17, & 
          ztscr18, & 
          ztscr19, & 
          ztscr20, & 
          ztscr21, & 
          zmscr01, & 
          zmscr02, & 
          zmscr03, & 
          zmscr04, & 
          zmscr05, & 
          zmscr06, & 
          zmscr07, & 
          zmscr08, & 
          zmscr09, & 
          zmscr10, & 
          zmscr11, & 
          zmscr12, & 
          zmscr13, & 
          zmscr14, & 
          zmscr15, & 
          zmscr16, & 
          zmscr17, & 
          l_netcdf, & 
          l_stats
#ifdef NETCDF
      use output_netcdf, only:  & 
          close_netcdf ! Procedure
#endif

      implicit none

      if ( l_stats .and. l_netcdf ) then
#ifdef NETCDF
        call close_netcdf( zt%f )
        call close_netcdf( zm%f )
        call close_netcdf( sfc%f )
#else
        stop "This program was not compiled with netCDF support"
#endif
      end if

      if ( l_stats ) then
        ! De-allocate all zt variables
        deallocate( zt%z )

        deallocate( zt%x )

        deallocate( zt%n )
        deallocate( zt%l_in_update )

        
        deallocate( zt%f%var )
        deallocate( zt%f%z )

        deallocate ( ztscr01 )
        deallocate ( ztscr02 )
        deallocate ( ztscr03 )
        deallocate ( ztscr04 )
        deallocate ( ztscr05 )
        deallocate ( ztscr06 )
        deallocate ( ztscr07 )
        deallocate ( ztscr08 )
        deallocate ( ztscr09 )
        deallocate ( ztscr10 )
        deallocate ( ztscr11 )
        deallocate ( ztscr12 )
        deallocate ( ztscr13 )
        deallocate ( ztscr14 )
        deallocate ( ztscr15 )
        deallocate ( ztscr16 )
        deallocate ( ztscr17 )
        deallocate ( ztscr18 )
        deallocate ( ztscr19 )
        deallocate ( ztscr20 )
        deallocate ( ztscr21 )

        ! De-allocate all zm variables
        deallocate( zm%z )

        deallocate( zm%x )
        deallocate( zm%n )

        deallocate( zm%f%var )
        deallocate( zm%f%z )
        deallocate( zm%l_in_update )
        
        deallocate ( zmscr01 )
        deallocate ( zmscr02 )
        deallocate ( zmscr03 )
        deallocate ( zmscr04 )
        deallocate ( zmscr05 )
        deallocate ( zmscr06 )
        deallocate ( zmscr07 )
        deallocate ( zmscr08 )
        deallocate ( zmscr09 )
        deallocate ( zmscr10 )
        deallocate ( zmscr11 )
        deallocate ( zmscr12 )
        deallocate ( zmscr13 )
        deallocate ( zmscr14 )
        deallocate ( zmscr15 )
        deallocate ( zmscr16 )
        deallocate ( zmscr17 )

        ! De-allocate all sfc variables
        deallocate( sfc%z )

        deallocate( sfc%x )
        deallocate( sfc%n )
        deallocate( sfc%l_in_update )

        deallocate( sfc%f%var )
        deallocate( sfc%f%z )
      end if ! l_stats

      return
      end subroutine stats_finalize

end module stats_subs
