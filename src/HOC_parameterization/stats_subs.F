!-----------------------------------------------------------------------
!  $Id: stats_subs.F,v 1.8 2008-04-05 19:54:32 griffinb Exp $
#ifdef STATS
      module stats_subs
 
      implicit none
      
      private
      
      public :: stats_init, stats_begin_timestep, stats_end_timestep,
     .          stats_accumulate, stats_finalize
      
      
      contains
      
!-----------------------------------------------------------------------
      subroutine stats_init( iunit, fnamelist, nnzp, gzt, gzm,
     .                       day, month, year, rlat, rlon, 
     .                       time_current, delt )

!     Description: Initializes the statistics saving functionality of
!     the HOC model.
!-----------------------------------------------------------------------
      use stats_hoc
      use stats_prec
      use outputgrads
      use outputnc
      use stats_zm
      use stats_zt
      use stats_sfc

      implicit none

      ! Constant Parameters
      ! Because the units are only open and closed on at a time, they can
      ! all use the same file unit with the scope of write_grads()
!     integer, parameter :: 
!    .  zt_funit  = 20,
!    .  zm_funit  = 20,
!    .  sfc_funit = 20
      ! dschanen changed this for threadsafe-ish I/O

      integer, parameter :: nvarmax = 250  ! Max variables

      ! Input Variables
      integer, intent(in) :: iunit  ! File unit for fnamelist

      character(len=*) :: fnamelist ! Filename holding the &statsnl

      integer, intent(in) :: nnzp ! Grid points in the vertical [count]

      real, intent(in), dimension(nnzp) :: 
     .  gzt, gzm  ! Thermodynamic and momentum levels           [m]

      integer, intent(in) :: day, month, year  ! Time of year

      real, intent(in) :: 
     .  rlat, rlon   ! Latitude and Longitude             [Degrees N/E]

      real(kind=time_prec), intent(in) :: 
     .  delt         ! Timestep (dtmain in HOC)           [s]

      real(kind=time_prec), intent(in) :: 
     .  time_current ! Model time                         [s]

      ! Local Variables

      ! Namelist Variables

      character(len=10) :: stats_fmt  ! File storage convention

      character(len=20), dimension(nvarmax) :: 
     .  vars_zt,  ! Variables on the thermodynamic levels
     .  vars_zm,  ! Variables on the momentum levels
     .  vars_sfc  ! Variables at the model surface

      namelist /statsnl/
     .  lstats, stats_fmt, stats_tsamp, stats_tout,
     .  fname_zt,  vars_zt,
     .  fname_zm,  vars_zm,
     .  fname_sfc, vars_sfc

      ! Local Variables

      logical :: lerror

      character(len=200) :: fdir, fname

      integer :: i, k, ntot

      ! Initialize
      lerror = .false.

      ! Read namelist

      open(unit=iunit, file=fnamelist)
      read(unit=iunit, nml=statsnl, end=100)
      close(unit=iunit)

      if ( .not. lstats ) then
        lstats_samp  = .false.
        lstats_first = .false.
        lstats_last  = .false.
        return
      end if

      ! Parse the file type for stats output.  Currently only GrADS and
      ! NetCDF v3 are supported by this code.

      select case( trim( stats_fmt ) ) 
      case( "GrADS", "grads", "gr" )
        lnetcdf = .false.
        lgrads  = .true.

      case ( "NetCDF", "netcdf", "nc" )
        lnetcdf = .true.
        lgrads  = .false.

      case default
        write(0,*) "Invalid data format "//trim( stats_fmt )
        stop

      end select

      ! Check sampling and output frequencies

      if ( abs( stats_tsamp/delt - int(stats_tsamp/delt) ) 
     .     > 1.e-8 ) then
         lerror = .true.
         write(0,*) 'Error: stats_tsamp should be a multiple of delt'
         write(0,*) 'stats_tsamp = ',stats_tsamp
         write(0,*) 'delt = ',delt
      end if

      if ( abs( stats_tout/stats_tsamp - int(stats_tout/stats_tsamp) ) 
     .     > 1.e-8 ) then
         lerror = .true.
         write(0,*) 
     .      'Error: stats_tout should be a multiple of stats_tsamp'
         write(0,*) 'stats_tout = ',stats_tout
         write(0,*) 'stats_tsamp = ',stats_tsamp
      end if

      ! Initialize zt (mass points)

      i = 1
      do while ( ichar(vars_zt(i)(1:1)) /= 0 
     .           .and. len_trim(vars_zt(i)) /= 0
     .           .and. i <= nvarmax )
        i = i + 1
      end do
      ntot = i - 1
      if ( ntot == nvarmax )
     .  write(0,*) 'WARNING: check nvarmax in statistics.f'

      zt%nn = ntot
      zt%kk = nnzp
!      write(*,*) 'Number of variables for zt ',zt%nn

      allocate( zt%z( zt%kk ) )
      zt%z = gzt

      allocate( zt%x( zt%kk, zt%nn ) )
      allocate( zt%n( zt%kk, zt%nn ) )
      call stats_zero( zt%kk, zt%nn, zt%x, zt%n )

      allocate( zt%f%var( zt%nn ) )
      allocate( zt%f%z( zt%kk ) )

      ! Allocate scratch space

      allocate ( ztscr01(zt%kk) )
      allocate ( ztscr02(zt%kk) )
      allocate ( ztscr03(zt%kk) )
      allocate ( ztscr04(zt%kk) )
      allocate ( ztscr05(zt%kk) )
      allocate ( ztscr06(zt%kk) )
      allocate ( ztscr07(zt%kk) )
      allocate ( ztscr08(zt%kk) )
      allocate ( ztscr09(zt%kk) )
      allocate ( ztscr10(zt%kk) )
      allocate ( ztscr11(zt%kk) )
      allocate ( ztscr12(zt%kk) )
      allocate ( ztscr13(zt%kk) )
      allocate ( ztscr14(zt%kk) )
      allocate ( ztscr15(zt%kk) )

      fdir = "./"
      fname = trim( fname_zt )

      if ( lgrads ) then

        ! Open GrADS file
        call open_grads( zt%f, iunit, fdir, fname, 
     .                   1, zt%kk, zt%z,
     .                   day, month, year, rlat, rlon,
     .                   time_current+stats_tout, stats_tout,
     .                   zt%nn )

      else ! Open NetCDF file
        call open_nc( zt%f, iunit, fdir, fname, 
     .                1, zt%kk, zt%z,
     .                day, month, year, rlat, rlon,
     .                time_current+stats_tout, stats_tout,
     .                zt%nn )

      end if

      ! Default initialization for array indices for zt

      call stats_init_zt( vars_zt, lerror )

      ! Initialize zm (momentum points)

      i = 1
      do while ( ichar(vars_zm(i)(1:1)) /= 0 
     .           .and. len_trim(vars_zm(i)) /= 0
     .           .and. i <= nvarmax )
        i = i + 1
      end do
      ntot = i - 1
      if ( ntot == nvarmax )
     .  write(0,*) 'WARNING: check nvarmax in statistics.f'

      zm%nn = ntot
      zm%kk = nnzp
!      write(*,*) 'Number of variables for zm ',zm%nn

      allocate( zm%z( zm%kk ) )
      zm%z = gzm

      allocate( zm%x( zm%kk, zm%nn ) )
      allocate( zm%n( zm%kk, zm%nn ) )

      call stats_zero( zm%kk, zm%nn, zm%x, zm%n )

      allocate( zm%f%var( zm%nn ) )
      allocate( zm%f%z( zm%kk ) )

      ! Allocate scratch space

      allocate ( zmscr01(zm%kk) )
      allocate ( zmscr02(zm%kk) )
      allocate ( zmscr03(zm%kk) )
      allocate ( zmscr04(zm%kk) )
      allocate ( zmscr05(zm%kk) )
      allocate ( zmscr06(zm%kk) )
      allocate ( zmscr07(zm%kk) )
      allocate ( zmscr08(zm%kk) )
      allocate ( zmscr09(zm%kk) )
      allocate ( zmscr10(zm%kk) )
      allocate ( zmscr11(zm%kk) )
      allocate ( zmscr12(zm%kk) )
      allocate ( zmscr13(zm%kk) )


      fdir = "./"
      fname = trim(fname_zm)
      if ( lgrads ) then

        ! Open GrADS files
        call open_grads( zm%f, iunit, fdir, fname, 
     .                   1, zm%kk, zm%z,
     .                   day, month, year, rlat, rlon,
     .                   time_current+stats_tout, stats_tout,
     .                   zm%nn )

      else ! Open NetCDF file
        call open_nc( zm%f, iunit, fdir, fname, 
     .                1, zm%kk, zm%z,
     .                day, month, year, rlat, rlon,
     .                time_current+stats_tout, stats_tout,
     .                zm%nn )
      end if

      call stats_init_zm( vars_zm, lerror )

      ! Initialize sfc (surface point)

      i = 1
      do while ( ichar(vars_sfc(i)(1:1)) /= 0 
     .           .and. len_trim(vars_sfc(i)) /= 0
     .           .and. i <= nvarmax )
        i = i + 1
      end do
      ntot = i - 1
      if ( ntot == nvarmax )
     .  write(0,*) 'WARNING: check nvarmax in statistics.f'

      sfc%nn = ntot
      sfc%kk = 1
!      write(*,*) 'Number of variables for sfc ',sfc%nn

      allocate( sfc%z( sfc%kk ) )
      sfc%z = gzm(1)

      allocate( sfc%x( sfc%kk, sfc%nn ) )
      allocate( sfc%n( sfc%kk, sfc%nn ) )
      call stats_zero( sfc%kk, sfc%nn, sfc%x, sfc%n )

      allocate( sfc%f%var( sfc%nn ) )
      allocate( sfc%f%z( sfc%kk ) )

      fdir = "./"
      fname = trim( fname_sfc )

      if ( lgrads ) then

        ! Open GrADS files
        call open_grads( sfc%f, iunit, fdir, fname, 
     .                   1, sfc%kk, sfc%z,
     .                   day, month, year, rlat, rlon,
     .                   time_current+stats_tout, stats_tout,
     .                   sfc%nn )

      else ! Open NetCDF files
        call open_nc( sfc%f, iunit, fdir, fname, 
     .                1, sfc%kk, sfc%z,
     .                day, month, year, rlat, rlon,
     .                time_current+stats_tout, stats_tout,
     .                sfc%nn )
      end if

      call stats_init_sfc( vars_sfc, lerror )

!     Check for errors

      if ( lerror ) then
        write(0,*) 'stats_init: errors found'
        stop
      end if

      return

!     If namelist was not found in input file, turn off statistics

100   continue
      write(0,*) 'Error with statsnl, statistics is turned off'
      lstats       = .false.
      lstats_samp  = .false.
      lstats_first = .false.
      lstats_last  = .false.

      return
      end subroutine stats_init
!-----------------------------------------------------------------------
      subroutine stats_zero( kk, nn, x, n )

!     Description:
!     Initialize stats to zero
!-----------------------------------------------------------------------
      use stats_prec

      implicit none

      ! Input
      integer, intent(in) :: kk, nn

      ! Output
      real(kind=stat_rknd), dimension(kk,nn), intent(out)    :: x
      integer(kind=stat_nknd), dimension(kk,nn), intent(out) :: n

      ! Zero out arrays

      if ( nn > 0 ) then
        x(:,:) = 0.0
        n(:,:) = 0
      end if

      return
      end subroutine stats_zero

!-----------------------------------------------------------------------
      subroutine stats_avg( kk, nn, x, n )

!     Description:
!     Compute the average of stats fields
!-----------------------------------------------------------------------
      use stats_prec

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

      use stats_hoc
      use stats_prec

      implicit none

      ! Input

      real(kind=time_prec), intent(in) :: 
     .  time_elapsed ! Elapsed model time       [s]

      real(kind=time_prec), intent(in) :: 
     .  delt         ! Model time step          [s]

      if ( .not. lstats ) return

      ! Set sample this time step flag
      if ( mod( time_elapsed, stats_tsamp ) < 1.e-8 ) then
        lstats_samp = .true.
      else
        lstats_samp = .false.
      end if

      ! Set first time step flag

      if ( mod( time_elapsed - delt, stats_tout ) < 1.e-8 ) then
        lstats_first = .true.
      else
        lstats_first = .false.
      end if

      ! Set last time step flag

      if ( mod( time_elapsed, stats_tout ) < 1.e-8 ) then
        lstats_last = .true.
      else
        lstats_last = .false.
      end if

      return

      end subroutine stats_begin_timestep

!-----------------------------------------------------------------------
      subroutine stats_end_timestep( time_elapsed, delt )

!     Description:
!-----------------------------------------------------------------------

      use stats_hoc
      use stats_prec
      use outputgrads
      use outputnc

      implicit none

      ! Input

      real(kind=time_prec), intent(in) :: 
     .  time_elapsed ! Time elapsed             [s]

      real(kind=time_prec), intent(in) :: 
     .  delt         ! timestep                 [s]

      ! Local Variables

      integer :: i, k

      ! Check if it is time to write to file

      if ( .not. lstats_last ) return

      ! Check number of sampling points

      do i=1,zt%nn
       do k=1,zt%kk
         if ( zt%n(k,i) /= 0 
     .        .and. zt%n(k,i) /= int(stats_tout/stats_tsamp) ) then
           write(0,*) 'Possible sampling error for zt ',i,k,zt%n(k,i)
           pause
         end if
        end do
      end do
      
      do i=1,zm%nn
       do k=1,zm%kk
         if ( zm%n(k,i) /= 0 
     .        .and. zm%n(k,i) /= int(stats_tout/stats_tsamp) ) then
           write(0,*) 'Possible sampling error for zm ',i,k,zm%n(k,i)
           pause
         end if
        end do
      end do
      
      ! Compute averages

      call stats_avg( zt%kk, zt%nn, zt%x, zt%n )
      call stats_avg( zm%kk, zm%nn, zm%x, zm%n )
      call stats_avg( sfc%kk, sfc%nn, sfc%x, sfc%n )

      ! Write to file
      if ( lgrads ) then
        call write_grads( zt%f  )
        call write_grads( zm%f  )
        call write_grads( sfc%f  )
      else ! lnetcdf
        call write_nc( zt%f  )
        call write_nc( zm%f  )
        call write_nc( sfc%f  )
      end if

      ! Reset sample fields
      call stats_zero( zt%kk, zt%nn, zt%x, zt%n )
      call stats_zero( zm%kk, zm%nn, zm%x, zm%n )
      call stats_zero( sfc%kk, sfc%nn, sfc%x, sfc%n )

      return
      end subroutine stats_end_timestep

!----------------------------------------------------------------------
      subroutine stats_accumulate
     .           ( um, vm, upwp, vpwp, up2, vp2, thlm,
     .             rtm, wprtp, wpthlp, wp2, wp3, rtp2, thlp2, rtpthlp,
     .             p, exner, rhot, rhom, thlm_forcing, rtm_forcing, wmm,
     .             wmt, Scm, taum, rcm, cf,
     .             sclrm, edsclrm, sclrm_forcing )

!     Description:
!     Accumulate those stats variables that are preserved in HOC from
!     timestep to timestep, but not those stats that are not, (e.g. 
!     budget terms, longwave and shortwave components, etc. )
!----------------------------------------------------------------------

      use stats_hoc
      use grid_class
      !use prognostic_variables
      use diagnostic_variables
      use model_flags
      use constants

      implicit none

      ! Input Variable
      real, intent(in), dimension(gr%nnzp) ::
     .  um,     ! u wind                        [m/s]
     .  vm,     ! v wind                        [m/s]
     .  upwp,   ! vertical u momentum flux      [m^2/s^2]
     .  vpwp,   ! vertical v momentum flux      [m^2/s^2]
     .  up2,    ! u'^2                          [m^2/s^2]
     .  vp2,    ! v'^2                          [m^2/s^2]
     .  thlm,   ! liquid potential temperature  [K]
     .  rtm,    ! total water mixing ratio      [kg/kg]
     .  wprtp,  ! w'rt'                         [m kg/s kg]
     .  wpthlp, ! w'thl'                        [m K /s]
     .  wp2,    ! w'^2                          [m^2/s^2]
     .  wp3,    ! w'^3                          [m^3/s^3]
     .  rtp2,   ! rt'^2                         [kg/kg]
     .  thlp2,  ! thl'^2                        [K^2]
     .  rtpthlp ! rt'thl'                       [kg/kg K]

      real, intent(in), dimension(gr%nnzp) ::
     .  p,           ! Pressure (Pa) on thermodynamic points    [Pa]
     .  exner,       ! Exner function = ( p / p0 ) ** kappa     [-]
     .  rhot,        ! Density                                  [kg/m^3]
     .  rhom,        ! Density                                  [kg/m^3]
     .  thlm_forcing,! thlm large-scale forcing                 [K/s]
     .  rtm_forcing, ! rtm large-scale forcing                  [kg/kg/s]
     .  wmm,         ! w on momentum levels                     [m/s]
     .  wmt,         ! w on thermodynamic levels                [m/s]
     .  Scm,         ! PDF width paramter                       [-]
     .  taum         ! Dissipation time                         [s]

      real, intent(in), dimension(gr%nnzp) ::
     .  rcm,  ! Cloud water mixing ratio                [kg/kg]
!    .  Ncm,  ! Cloud droplet number concentration      [num/kg]
!    .  Ncnm, ! Cloud nuclei number concentration       [num/m^3]
!    .  Nim,  ! Ice nuclei number concentration         [num/m^3]
     .  cf    ! Cloud fraction                          [%]

      real, intent(in), optional, dimension(gr%nnzp,sclr_dim) ::
     .  sclrm,          ! High-Order Passive scalar     [units vary]
     .  edsclrm,        ! Eddy-diff Passive scalar      [units vary] 
     .  sclrm_forcing   ! Large-scale forcing of scalar [units/s]

      ! Prognostic drizzle variable array
!     real, intent(in), dimension(gr%nnzp,hydromet_dim) :: hydromet
      ! Contains:
      ! 1 rrm      Rain water mixing ratio               [kg/kg]
      ! 2 Nrm      Rain droplet number concentration     [num/kg]
      ! 3 rsnow    Snow water mixing ratio               [kg/kg]
      ! 4 rice     Ice water mixing ratio                [kg/kg]
      ! 5 rgraupel Graupel water mixing ratio            [kg/kg]

      ! Local Variables

      integer :: i, k

      real :: xtmp

      ! Sample fields

      if ( lstats_samp ) then

        ! zt variables

        if ( ithlm > 0 ) then
          zt%x(:,ithlm) = zt%x(:,ithlm) + thlm
          zt%n(:,ithlm) = zt%n(:,ithlm) + 1
        end if
        if ( ithvm > 0 ) then
          zt%x(:,ithvm) = zt%x(:,ithvm) + thvm
          zt%n(:,ithvm) = zt%n(:,ithvm) + 1
        end if
        if ( irtm > 0 ) then
          zt%x(:,irtm) = zt%x(:,irtm) + rtm
          zt%n(:,irtm) = zt%n(:,irtm) + 1
        end if
        if ( ircm > 0 ) then
          zt%x(:,ircm) = zt%x(:,ircm) + rcm
          zt%n(:,ircm) = zt%n(:,ircm) + 1
        end if
        if ( ium > 0 ) then
          zt%x(:,ium) = zt%x(:,ium) + um
          zt%n(:,ium) = zt%n(:,ium) + 1
        end if
        if ( ivm > 0 ) then
          zt%x(:,ivm) = zt%x(:,ivm) + vm
          zt%n(:,ivm) = zt%n(:,ivm) + 1
        end if
        if ( iwmt > 0 ) then
          zt%x(:,iwmt) = zt%x(:,iwmt) + wmt
          zt%n(:,iwmt) = zt%n(:,iwmt) + 1
        end if
        if ( iug > 0 ) then
          zt%x(:,iug) = zt%x(:,iug) + ug
          zt%n(:,iug) = zt%n(:,iug) + 1
        end if
        if ( ivg > 0 ) then
          zt%x(:,ivg) = zt%x(:,ivg) + vg
          zt%n(:,ivg) = zt%n(:,ivg) + 1
        end if
        if ( icf > 0 ) then
          zt%x(:,icf) = zt%x(:,icf) + cf
          zt%n(:,icf) = zt%n(:,icf) + 1
        end if
        if ( ip > 0 ) then
          zt%x(:,ip) = zt%x(:,ip) + p
          zt%n(:,ip) = zt%n(:,ip) + 1
        end if
        if ( iexner > 0 ) then
          zt%x(:,iexner) = zt%x(:,iexner) + exner
          zt%n(:,iexner) = zt%n(:,iexner) + 1
        end if
        if ( iLscale > 0 ) then
          zt%x(:,iLscale) = zt%x(:,iLscale) + Lscale
          zt%n(:,iLscale) = zt%n(:,iLscale) + 1
        end if
        if ( iwp3 > 0 ) then
          zt%x(:,iwp3) = zt%x(:,iwp3) + wp3
          zt%n(:,iwp3) = zt%n(:,iwp3) + 1
        end if
        if ( iwpthlp2 > 0 ) then
          zt%x(:,iwpthlp2) = zt%x(:,iwpthlp2) + wpthlp2
          zt%n(:,iwpthlp2) = zt%n(:,iwpthlp2) + 1
        end if
        if ( iwp2thlp > 0 ) then
          zt%x(:,iwp2thlp) = zt%x(:,iwp2thlp) + wp2thlp
          zt%n(:,iwp2thlp) = zt%n(:,iwp2thlp) + 1
        end if
        if ( iwprtp2 > 0 ) then
          zt%x(:,iwprtp2) = zt%x(:,iwprtp2) + wprtp2
          zt%n(:,iwprtp2) = zt%n(:,iwprtp2) + 1
        end if
        if ( iwp2rtp > 0 ) then
          zt%x(:,iwp2rtp) = zt%x(:,iwp2rtp) + wp2rtp
          zt%n(:,iwp2rtp) = zt%n(:,iwp2rtp) + 1
        end if
        if ( iLup > 0 ) then
          zt%x(:,iLup) = zt%x(:,iLup) + Lup
          zt%n(:,iLup) = zt%n(:,iLup) + 1
        end if
        if ( iLdown > 0 ) then
          zt%x(:,iLdown) = zt%x(:,iLdown) + Ldown
          zt%n(:,iLdown) = zt%n(:,iLdown) + 1
        end if
        if ( itaut > 0 ) then
          zt%x(:,itaut) = zt%x(:,itaut) + taut
          zt%n(:,itaut) = zt%n(:,itaut) + 1
        end if
        if ( iKht > 0 ) then
          zt%x(:,iKht) = zt%x(:,iKht) + Kht
          zt%n(:,iKht) = zt%n(:,iKht) + 1
        end if
        if ( iwp2thvp > 0 ) then
          zt%x(:,iwp2thvp) = zt%x(:,iwp2thvp) + wp2thvp
          zt%n(:,iwp2thvp) = zt%n(:,iwp2thvp) + 1
        end if
        if ( iwp2rcp > 0 ) then
          zt%x(:,iwp2rcp) = zt%x(:,iwp2rcp) + wp2rcp
          zt%n(:,iwp2rcp) = zt%n(:,iwp2rcp) + 1
        end if
        if ( iwprtpthlp > 0 ) then
          zt%x(:,iwprtpthlp) = zt%x(:,iwprtpthlp) + wprtpthlp
          zt%n(:,iwprtpthlp) = zt%n(:,iwprtpthlp) + 1
        end if
        if ( isct > 0 ) then
          zt%x(:,isct) = zt%x(:,isct) + sct
          zt%n(:,isct) = zt%n(:,isct) + 1
        end if
        if ( irhot > 0 ) then
          zt%x(:,irhot) = zt%x(:,irhot) + rhot
          zt%n(:,irhot) = zt%n(:,irhot) + 1
        end if
!       if ( iNcm > 0 ) then           ! Brian
!         zt%x(:,iNcm) = zt%x(:,iNcm) + Ncm
!         zt%n(:,iNcm) = zt%n(:,iNcm) + 1
!       end if
!       if ( iNcnm > 0 ) then
!         zt%x(:,iNcnm) = zt%x(:,iNcnm) + Ncnm
!         zt%n(:,iNcnm) = zt%n(:,iNcnm) + 1
!       end if
!       if ( iNim > 0 ) then
!         zt%x(:,iNim) = zt%x(:,iNim) + Nim
!         zt%n(:,iNim) = zt%n(:,iNim) + 1
!       end if
!       if ( ised_rcm > 0 .and. cloud_sed ) then       ! Brian
!         zt%x(:,ised_rcm) = zt%x(:,ised_rcm) + sed_rcm
!         zt%n(:,ised_rcm) = zt%n(:,ised_rcm) + 1
!       end if
        if ( irsm > 0 ) then           ! Brian
          zt%x(:,irsm) = zt%x(:,irsm) + rsm
          zt%n(:,irsm) = zt%n(:,irsm) + 1
        end if

!       if ( irrm > 0 ) then     ! Brian
!         zt%x(:,irrm) = zt%x(:,irrm) + hydromet(:,1)
!         zt%n(:,irrm) = zt%n(:,irrm) + 1
!       end if

!       if ( iNrm > 0 ) then     ! Brian
!         zt%x(:,iNrm) = zt%x(:,iNrm) + hydromet(:,2)
!         zt%n(:,iNrm) = zt%n(:,iNrm) + 1
!       end if

!       if ( irsnowm > 0 ) then
!         zt%x(:,irsnowm) = zt%x(:,irsnowm) + hydromet(:,3)
!         zt%n(:,irsnowm) = zt%n(:,irsnowm) + 1
!       end if

!       if ( iricem > 0 ) then
!         zt%x(:,iricem) = zt%x(:,iricem) + hydromet(:,4)
!         zt%n(:,iricem) = zt%n(:,iricem) + 1
!       end if

!       if ( irgraupelm > 0 ) then
!         zt%x(:,irgraupelm) = zt%x(:,irgraupelm) + hydromet(:,5)
!         zt%n(:,irgraupelm) = zt%n(:,irgraupelm) + 1
!       end if

        if ( iAKm > 0 .and. LH_on ) then         ! Vince Larson 22 May 2005
          zt%x(:,iAKm) = zt%x(:,iAKm) + AKm
          zt%n(:,iAKm) = zt%n(:,iAKm) + 1
        end if

        if ( iAKm_est > 0 .and. LH_on ) then     ! Vince Larson 22 May 2005
          zt%x(:,iAKm_est) = zt%x(:,iAKm_est) + AKm_est
          zt%n(:,iAKm_est) = zt%n(:,iAKm_est) + 1
        end if

        if ( iradht > 0 ) then
          zt%x(:,iradht) = zt%x(:,iradht) + radht
          zt%n(:,iradht) = zt%n(:,iradht) + 1
        end if

        if ( ia > 0 ) then
          zt%x(:,ia) = zt%x(:,ia) + pdf_parms(:,13)
          zt%n(:,ia) = zt%n(:,ia) + 1
        end if

        if ( iw1 > 0 ) then
          zt%x(:,iw1) = zt%x(:,iw1) + pdf_parms(:,1)
          zt%n(:,iw1) = zt%n(:,iw1) + 1
        end if

        if ( iw2 > 0 ) then
          zt%x(:,iw2) = zt%x(:,iw2) + pdf_parms(:,2)
          zt%n(:,iw2) = zt%n(:,iw2) + 1
        end if

        if ( isw1 > 0 ) then
          zt%x(:,isw1) = zt%x(:,isw1) + pdf_parms(:,3)
          zt%n(:,isw1) = zt%n(:,isw1) + 1
        end if

        if ( isw2 > 0 ) then
          zt%x(:,isw2) = zt%x(:,isw2) + pdf_parms(:,4)
          zt%n(:,isw2) = zt%n(:,isw2) + 1
        end if

        if ( ithl1 > 0 ) then
          zt%x(:,ithl1) = zt%x(:,ithl1) + pdf_parms(:,9)
          zt%n(:,ithl1) = zt%n(:,ithl1) + 1
        end if

        if ( ithl2 > 0 ) then
          zt%x(:,ithl2) = zt%x(:,ithl2) + pdf_parms(:,10)
          zt%n(:,ithl2) = zt%n(:,ithl2) + 1
        end if

        if ( isthl1 > 0 ) then
          zt%x(:,isthl1) = zt%x(:,isthl1) + pdf_parms(:,11)
          zt%n(:,isthl1) = zt%n(:,isthl1) + 1
        end if

        if ( isthl2 > 0 ) then
          zt%x(:,isthl2) = zt%x(:,isthl2) + pdf_parms(:,12)
          zt%n(:,isthl2) = zt%n(:,isthl2) + 1
        end if

        if ( irt1 > 0 ) then
          zt%x(:,irt1) = zt%x(:,irt1) + pdf_parms(:,5)
          zt%n(:,irt1) = zt%n(:,irt1) + 1
        end if

        if ( irt2 > 0 ) then
          zt%x(:,irt2) = zt%x(:,irt2) + pdf_parms(:,6)
          zt%n(:,irt2) = zt%n(:,irt2) + 1
        end if

        if ( isrt1 > 0 ) then
          zt%x(:,isrt1) = zt%x(:,isrt1) + pdf_parms(:,7)
          zt%n(:,isrt1) = zt%n(:,isrt1) + 1
        end if

        if ( isrt2 > 0 ) then
          zt%x(:,isrt2) = zt%x(:,isrt2) + pdf_parms(:,8)
          zt%n(:,isrt2) = zt%n(:,isrt2) + 1
        end if

        if ( irc1 > 0 ) then
          zt%x(:,irc1) = zt%x(:,irc1) + pdf_parms(:,14)
          zt%n(:,irc1) = zt%n(:,irc1) + 1
        end if

        if ( irc2 > 0 ) then
          zt%x(:,irc2) = zt%x(:,irc2) + pdf_parms(:,15)
          zt%n(:,irc2) = zt%n(:,irc2) + 1
        end if

        if ( irsl1 > 0 ) then
          zt%x(:,irsl1) = zt%x(:,irsl1) + pdf_parms(:,16)
          zt%n(:,irsl1) = zt%n(:,irsl1) + 1
        end if

        if ( irsl2 > 0 ) then
          zt%x(:,irsl2) = zt%x(:,irsl2) + pdf_parms(:,17)
          zt%n(:,irsl2) = zt%n(:,irsl2) + 1
        end if

        if ( iR1 > 0 ) then
          zt%x(:,iR1) = zt%x(:,iR1) + pdf_parms(:,18)
          zt%n(:,iR1) = zt%n(:,iR1) + 1
        end if

        if ( iR2 > 0 ) then
          zt%x(:,iR2) = zt%x(:,iR2) + pdf_parms(:,19)
          zt%n(:,iR2) = zt%n(:,iR2) + 1
        end if

        if ( is1 > 0 ) then
          zt%x(:,is1) = zt%x(:,is1) + pdf_parms(:,20)
          zt%n(:,is1) = zt%n(:,is1) + 1
        end if

        if ( is2 > 0 ) then
          zt%x(:,is2) = zt%x(:,is2) + pdf_parms(:,21)
          zt%n(:,is2) = zt%n(:,is2) + 1
        end if

        if ( iss1 > 0 ) then
          zt%x(:,iss1) = zt%x(:,iss1) + pdf_parms(:,22)
          zt%n(:,iss1) = zt%n(:,iss1) + 1
        end if

        if ( iss2 > 0 ) then
          zt%x(:,iss2) = zt%x(:,iss2) + pdf_parms(:,23)
          zt%n(:,iss2) = zt%n(:,iss2) + 1
        end if

        if ( irrtthl > 0 ) then
          zt%x(:,irrtthl) = zt%x(:,irrtthl) + pdf_parms(:,24)
          zt%n(:,irrtthl) = zt%n(:,irrtthl) + 1
        end if

        if ( iwp2zt > 0 ) then
          zt%x(:,iwp2zt) = zt%x(:,iwp2zt) + wp2zt
          zt%n(:,iwp2zt) = zt%n(:,iwp2zt) + 1
        end if

        if ( ithlp2zt > 0 ) then
          zt%x(:,ithlp2zt) = zt%x(:,ithlp2zt) + thlp2zt
          zt%n(:,ithlp2zt) = zt%n(:,ithlp2zt) + 1
        end if

        if ( iwpthlpzt > 0 ) then
          zt%x(:,iwpthlpzt) = zt%x(:,iwpthlpzt) + wpthlpzt
          zt%n(:,iwpthlpzt) = zt%n(:,iwpthlpzt) + 1
        end if

        if ( iwprtpzt > 0 ) then
          zt%x(:,iwprtpzt) = zt%x(:,iwprtpzt) + wprtpzt
          zt%n(:,iwprtpzt) = zt%n(:,iwprtpzt) + 1
        end if

        if ( irtp2zt > 0 ) then
          zt%x(:,irtp2zt) = zt%x(:,irtp2zt) + rtp2zt
          zt%n(:,irtp2zt) = zt%n(:,irtp2zt) + 1
        end if

        if ( irtpthlpzt > 0 ) then
          zt%x(:,irtpthlpzt) = zt%x(:,irtpthlpzt) + rtpthlpzt
          zt%n(:,irtpthlpzt) = zt%n(:,irtpthlpzt) + 1
        end if
#ifdef SCALARS
        if ( isclram > 0 ) then
          zt%x(:,isclram) = zt%x(:,isclram) + sclrm(:,1)
          zt%n(:,isclram) = zt%n(:,isclram) + 1
        end if

        if ( isclram_f > 0 ) then
          zt%x(:,isclram_f) = zt%x(:,isclram_f) + sclrm_forcing(:,1)
          zt%n(:,isclram_f) = zt%n(:,isclram_f) + 1
        end if

        if ( isclrbm > 0 ) then
          zt%x(:,isclrbm) = zt%x(:,isclrbm) + sclrm(:,2)
          zt%n(:,isclrbm) = zt%n(:,isclrbm) + 1
        end if

        if ( isclrbm_f > 0 ) then
          zt%x(:,isclrbm_f) = zt%x(:,isclrbm_f) + sclrm_forcing(:,2)
          zt%n(:,isclrbm_f) = zt%n(:,isclrbm_f) + 1
        end if

        if ( iedsclram > 0 ) then
          zt%x(:,iedsclram) = zt%x(:,iedsclram) + edsclrm(:,1)
          zt%n(:,iedsclram) = zt%n(:,iedsclram) + 1
        end if

        if ( iedsclrbm > 0 ) then
          zt%x(:,iedsclrbm) = zt%x(:,iedsclrbm) + edsclrm(:,2)
          zt%n(:,iedsclrbm) = zt%n(:,iedsclrbm) + 1
        end if

#endif /*SCALARS*/

!       zm variables

        if ( iwp2 > 0 ) then
          zm%x(:,iwp2) = zm%x(:,iwp2) + wp2
          zm%n(:,iwp2) = zm%n(:,iwp2) + 1
        end if

        if ( irtp2 > 0 ) then
          zm%x(:,irtp2) = zm%x(:,irtp2) + rtp2
          zm%n(:,irtp2) = zm%n(:,irtp2) + 1
        end if

        if ( ithlp2 > 0 ) then
          zm%x(:,ithlp2) = zm%x(:,ithlp2) + thlp2
          zm%n(:,ithlp2) = zm%n(:,ithlp2) + 1
        end if

        if ( irtpthlp > 0 ) then
          zm%x(:,irtpthlp) = zm%x(:,irtpthlp) + rtpthlp
          zm%n(:,irtpthlp) = zm%n(:,irtpthlp) + 1
        end if

        if ( iwprtp > 0 ) then
          zm%x(:,iwprtp) = zm%x(:,iwprtp) + wprtp
          zm%n(:,iwprtp) = zm%n(:,iwprtp) + 1
        end if

        if ( iwpthlp > 0 ) then
          zm%x(:,iwpthlp) = zm%x(:,iwpthlp) + wpthlp
          zm%n(:,iwpthlp) = zm%n(:,iwpthlp) + 1
        end if

        if ( iwp4 > 0 ) then
          zm%x(:,iwp4) = zm%x(:,iwp4) + wp4
          zm%n(:,iwp4) = zm%n(:,iwp4) + 1
        end if

        if ( iwpthvp > 0 ) then
          zm%x(:,iwpthvp) = zm%x(:,iwpthvp) + wpthvp
          zm%n(:,iwpthvp) = zm%n(:,iwpthvp) + 1
        end if

        if ( irtpthvp > 0 ) then
          zm%x(:,irtpthvp) = zm%x(:,irtpthvp) + rtpthvp
          zm%n(:,irtpthvp) = zm%n(:,irtpthvp) + 1
        end if

        if ( ithlpthvp > 0 ) then
          zm%x(:,ithlpthvp) = zm%x(:,ithlpthvp) + thlpthvp
          zm%n(:,ithlpthvp) = zm%n(:,ithlpthvp) + 1
        end if

        if ( itaum > 0 ) then
          zm%x(:,itaum) = zm%x(:,itaum) + taum
          zm%n(:,itaum) = zm%n(:,itaum) + 1
        end if

        if ( iKhm > 0 ) then
          zm%x(:,iKhm) = zm%x(:,iKhm) + Khm
          zm%n(:,iKhm) = zm%n(:,iKhm) + 1
        end if

        if ( iwprcp > 0 ) then
          zm%x(:,iwprcp) = zm%x(:,iwprcp) + wprcp
          zm%n(:,iwprcp) = zm%n(:,iwprcp) + 1
        end if

        if ( ithlprcp > 0 ) then
          zm%x(:,ithlprcp) = zm%x(:,ithlprcp) + thlprcp
          zm%n(:,ithlprcp) = zm%n(:,ithlprcp) + 1
        end if

        if ( irtprcp > 0 ) then
          zm%x(:,irtprcp) = zm%x(:,irtprcp) + rtprcp
          zm%n(:,irtprcp) = zm%n(:,irtprcp) + 1
        end if

        if ( ircp2 > 0 ) then
          zm%x(:,ircp2) = zm%x(:,ircp2) + rcp2
          zm%n(:,ircp2) = zm%n(:,ircp2) + 1
        end if

        if ( iupwp > 0 ) then
          zm%x(:,iupwp) = zm%x(:,iupwp) + upwp
          zm%n(:,iupwp) = zm%n(:,iupwp) + 1
        end if

        if ( ivpwp > 0 ) then
          zm%x(:,ivpwp) = zm%x(:,ivpwp) + vpwp
          zm%n(:,ivpwp) = zm%n(:,ivpwp) + 1
        end if

        if ( ivp2 > 0 ) then
          zm%x(:,ivp2) = zm%x(:,ivp2) + vp2
          zm%n(:,ivp2) = zm%n(:,ivp2) + 1
        end if

        if ( iup2 > 0 ) then
          zm%x(:,iup2) = zm%x(:,iup2) + up2
          zm%n(:,iup2) = zm%n(:,iup2) + 1
        end if

        if ( irhom > 0 ) then
          zm%x(:,irhom) = zm%x(:,irhom) + rhom
          zm%n(:,irhom) = zm%n(:,irhom) + 1
        end if

        if ( iscm > 0 ) then
          zm%x(:,iscm) = zm%x(:,iscm) + scm
          zm%n(:,iscm) = zm%n(:,iscm) + 1
        end if

        if ( iem > 0 ) then
          zm%x(:,iem) = zm%x(:,iem) + em
          zm%n(:,iem) = zm%n(:,iem) + 1
        end if

        if ( ishear > 0 ) then    ! Brian
          zm%x(:,ishear) = zm%x(:,ishear) + shear
          zm%n(:,ishear) = zm%n(:,ishear) + 1
        end if

        if ( iFrad > 0 ) then
          zm%x(:,iFrad) = zm%x(:,iFrad) + Frad
          zm%n(:,iFrad) = zm%n(:,iFrad) + 1
        end if

!       if ( iFcsed > 0 .and. cloud_sed ) then    ! Brian
!         zm%x(:,iFcsed) = zm%x(:,iFcsed) + Fcsed
!         zm%n(:,iFcsed) = zm%n(:,iFcsed) + 1
!       end if

#ifdef SCALARS
        if ( isclraprtp > 0 ) then
          zm%x(:,isclraprtp) = zm%x(:,isclraprtp) + sclrprtp(:,1)
          zm%n(:,isclraprtp) = zm%n(:,isclraprtp) + 1
        end if
        if ( isclrbprtp > 0 ) then
          zm%x(:,isclrbprtp) = zm%x(:,isclrbprtp) + sclrprtp(:,2)
          zm%n(:,isclrbprtp) = zm%n(:,isclrbprtp) + 1
        end if
        if ( isclrap2 > 0 ) then
          zm%x(:,isclrap2) = zm%x(:,isclrap2) + sclrp2(:,1)
          zm%n(:,isclrap2) = zm%n(:,isclrap2) + 1
        end if
        if ( isclrbp2 > 0 ) then
          zm%x(:,isclrbp2) = zm%x(:,isclrbp2) + sclrp2(:,2)
          zm%n(:,isclrbp2) = zm%n(:,isclrbp2) + 1
        end if
        if ( isclrapthvp > 0 ) then
          zm%x(:,isclrapthvp) = zm%x(:,isclrapthvp) 
     .                            + sclrpthvp(:,1)
          zm%n(:,isclrapthvp) = zm%n(:,isclrapthvp) + 1
        end if
        if ( isclrbpthvp > 0 ) then
          zm%x(:,isclrbpthvp) = zm%x(:,isclrbpthvp) 
     .                            + sclrpthvp(:,2)
          zm%n(:,isclrbpthvp) = zm%n(:,isclrbpthvp) + 1
        end if
        if ( isclrapthlp > 0 ) then
          zm%x(:,isclrapthlp) = zm%x(:,isclrapthlp) 
     .                            + sclrpthlp(:,1)
          zm%n(:,isclrapthlp) = zm%n(:,isclrapthlp) + 1
        end if
        if ( isclrbpthlp > 0 ) then
          zm%x(:,isclrbpthlp) = zm%x(:,isclrbpthlp) 
     .                            + sclrpthlp(:,2)
          zm%n(:,isclrbpthlp) = zm%n(:,isclrbpthlp) + 1
        end if
        if ( isclraprcp > 0 ) then
          zm%x(:,isclraprcp) = zm%x(:,isclraprcp) + sclrprcp(:,1)
          zm%n(:,isclraprcp) = zm%n(:,isclraprcp) + 1
        end if
        if ( isclrbprcp > 0 ) then
          zm%x(:,isclrbprcp) = zm%x(:,isclrbprcp) + sclrprcp(:,2)
          zm%n(:,isclrbprcp) = zm%n(:,isclrbprcp) + 1
        end if
        if ( iwpsclrap > 0 ) then
          zm%x(:,iwpsclrap) = zm%x(:,iwpsclrap) + wpsclrp(:,1)
          zm%n(:,iwpsclrap) = zm%n(:,iwpsclrap) + 1
        end if
        if ( iwpsclrbp > 0 ) then
          zm%x(:,iwpsclrbp) = zm%x(:,iwpsclrbp) + wpsclrp(:,2)
          zm%n(:,iwpsclrbp) = zm%n(:,iwpsclrbp) + 1
        end if
        if ( iwp2sclrap > 0 ) then
          zm%x(:,iwp2sclrap) = zm%x(:,iwp2sclrap) + wp2sclrp(:,1)
          zm%n(:,iwp2sclrap) = zm%n(:,iwp2sclrap) + 1
        end if
        if ( iwp2sclrbp > 0 ) then
          zm%x(:,iwp2sclrbp) = zm%x(:,iwp2sclrbp) + wp2sclrp(:,2)
          zm%n(:,iwp2sclrbp) = zm%n(:,iwp2sclrbp) + 1
        end if
        if ( iwpsclrap2 > 0 ) then
          zm%x(:,iwpsclrap2) = zm%x(:,iwpsclrap2) + wpsclrp2(:,1)
          zm%n(:,iwpsclrap2) = zm%n(:,iwpsclrap2) + 1
        end if
        if ( iwpsclrbp2 > 0 ) then
          zm%x(:,iwpsclrbp2) = zm%x(:,iwpsclrbp2) + wpsclrp2(:,2)
          zm%n(:,iwpsclrbp2) = zm%n(:,iwpsclrbp2) + 1
        end if
        if ( iwpsclraprtp > 0 ) then
          zm%x(:,iwpsclraprtp) = zm%x(:,iwpsclraprtp) + wpsclrprtp(:,1)
          zm%n(:,iwpsclraprtp) = zm%n(:,iwpsclraprtp) + 1
        end if
        if ( iwpsclrbprtp > 0 ) then
          zm%x(:,iwpsclrbprtp) = zm%x(:,iwpsclrbprtp) + wpsclrprtp(:,2)
          zm%n(:,iwpsclrbprtp) = zm%n(:,iwpsclrbprtp) + 1
        end if
        if ( iwpsclrapthlp > 0 ) then
          zm%x(:,iwpsclrapthlp) = zm%x(:,iwpsclrapthlp) 
     .                            + wpsclrpthlp(:,1)
          zm%n(:,iwpsclrapthlp) = zm%n(:,iwpsclrapthlp) + 1
        end if
        if ( iwpsclrbpthlp > 0 ) then
          zm%x(:,iwpsclrbpthlp) = zm%x(:,iwpsclrbpthlp) 
     .                            + wpsclrpthlp(:,2)
          zm%n(:,iwpsclrbpthlp) = zm%n(:,iwpsclrbpthlp) + 1
        end if
        if ( iwpedsclrap > 0 ) then
          zm%x(:,iwpedsclrap) = zm%x(:,iwpedsclrap) 
     .                            + wpedsclrp(:,1)
          zm%n(:,iwpedsclrap) = zm%n(:,iwpedsclrap) + 1
        end if
        if ( iwpedsclrbp > 0 ) then
          zm%x(:,iwpedsclrbp) = zm%x(:,iwpedsclrbp) 
     .                            + wpedsclrp(:,2)
          zm%n(:,iwpedsclrbp) = zm%n(:,iwpedsclrbp) + 1
        end if
        
#endif /*SCALARS*/

! sfc variables

!       Cloud cover
        if ( icc > 0 ) then
          sfc%x(1,icc) = sfc%x(1,icc) + maxval( cf(1:gr%nnzp) )
          sfc%n(1,icc) = sfc%n(1,icc) + 1
        end if

!       Cloud base
        if ( izb > 0 ) then
          k = 1
          do while ( rcm(k) < rc_tol .and. k < gr%nnzp )
            k = k + 1
          end do

!          if ( k < gr%nnzp) then
!            sfc%x(1,izb) = sfc%x(1,izb) + gr%zt(k)
          if ( k > 1 .AND. k < gr%nnzp) then
            sfc%x(1,izb) = sfc%x(1,izb) + (                     ! Use linear interpolation
     .         ( ( rc_tol - rcm(k-1) )/( rcm(k) - rcm(k-1) ) )  ! to find the exact height
     .        * ( gr%zt(k) - gr%zt(k-1) ) + gr%zt(k-1)          ! of the rc_tol kg/kg level.
     .                                    )                     ! Brian.
          else
            sfc%x(1,izb) = sfc%x(1,izb) - 10
          end if
          sfc%n(1,izb) = sfc%n(1,izb) + 1
        end if

!       LWP
        if ( ilwp > 0 ) then
          xtmp = 0.
          do i = gr%nnzp-1, 1, -1
           xtmp = xtmp + rhot(i+1) * rcm(i+1) / gr%dzt(i+1)
          end do
          sfc%x(1,ilwp) = sfc%x(1,ilwp) + xtmp
          sfc%n(1,ilwp) = sfc%n(1,ilwp) + 1
        end if

      end if

      return
      end subroutine stats_accumulate

!-----------------------------------------------------------------------
      subroutine stats_finalize( )

!     Description:
!     Close NetCDF files and deallocate scratch space and 
!     stats file structures.
!-----------------------------------------------------------------------

      use stats_hoc
      use outputnc, only: close_nc

      implicit none

      integer :: i

      if ( lstats .and. lnetcdf ) then
        call close_nc( zt%f )
        call close_nc( zm%f )
        call close_nc( sfc%f )
      end if

      if ( lstats ) then
        ! De-allocate all zt variables
        deallocate( zt%z )

        deallocate( zt%x, stat=i )
        deallocate( zt%n )

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

        ! De-allocate all zm variables
        deallocate( zm%z )

        deallocate( zm%x )
        deallocate( zm%n )

        deallocate( zm%f%var )
        deallocate( zm%f%z )

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

        ! De-allocate all sfc variables
        deallocate( sfc%z )

        deallocate( sfc%x )
        deallocate( sfc%n )

        deallocate( sfc%f%var )
        deallocate( sfc%f%z )
      end if ! lstats

      return
      end subroutine stats_finalize


      end module stats_subs
          
#endif /*STATS*/
