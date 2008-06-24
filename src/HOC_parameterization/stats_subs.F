!-----------------------------------------------------------------------
!  $Id: stats_subs.F,v 1.22 2008-06-24 19:44:02 faschinj Exp $
      module stats_subs
#ifdef STATS
 
      implicit none
      
      private ! Set Default Scope
      
      public :: stats_init, stats_begin_timestep, stats_end_timestep,
     .          stats_accumulate, stats_finalize
      
      private :: stats_zero, stats_avg
      
      contains
      
!-----------------------------------------------------------------------
      subroutine stats_init
     .           ( iunit, fname_prefix,
     .             lstats_in, stats_fmt_in, 
     .             stats_tsamp_in, stats_tout_in,
     .             fnamelist, nnzp, gzt, gzm,
     .             day, month, year, rlat, rlon, 
     .             time_current, delt )


!     Description: Initializes the statistics saving functionality of
!     the HOC model.
!-----------------------------------------------------------------------
      use stats_variables, only:
     .    zt,     ! Variables
     .    ztscr01,
     .    ztscr02,
     .    ztscr03,
     .    ztscr04,
     .    ztscr05,
     .    ztscr06,
     .    ztscr07,
     .    ztscr08,
     .    ztscr09,
     .    ztscr10,
     .    ztscr11,
     .    ztscr12,
     .    ztscr13,
     .    ztscr14,
     .    ztscr15,
     .    ztscr16,
     .    zm,
     .    zmscr01,
     .    zmscr02,
     .    zmscr03,
     .    zmscr04,
     .    zmscr05,
     .    zmscr06,
     .    zmscr07,
     .    zmscr08,
     .    zmscr09,
     .    zmscr10,
     .    zmscr11,
     .    zmscr12,
     .    zmscr13,
     .    zmscr14,
     .    sfc,
     .    lstats,
     .    stats_tsamp,
     .    stats_tout,
     .    lstats_samp,
     .    lstats_first,
     .    lstats_last,
     .    fname_zt,
     .    fname_zm,
     .    fname_sfc,
     .    lnetcdf,
     .    lgrads
      use stats_precision, only:
     .    time_precision   ! Variable(s)
      use outputgrads, only:
     .    open_grads  ! Procedure
#ifdef NETCDF
      use outputnc, only:
     .    open_nc     ! Procedure
#endif
      use stats_zm, only:
     .    stats_init_zm ! Procedure
      use stats_zt, only:
     .    stats_init_zt ! Procedure
      use stats_sfc, only:
     .    stats_init_sfc ! Procedure
      implicit none

      ! Constant Parameters
 
      integer, parameter :: nvarmax = 250  ! Max variables

      ! Input Variables

      integer, intent(in) :: iunit  ! File unit for fnamelist

      character(len=*) :: 
     .  fnamelist,   ! Filename holding the &statsnl
     .  fname_prefix,! Start of the stats filenames
     .  stats_fmt_in ! Format of the stats file output

      logical, intent(in) :: lstats_in ! Stats on? T/F

      real(kind=time_precision) :: 
     .  stats_tsamp_in, ! Sampling interval   [s]
     .  stats_tout_in   ! Output interval     [s]

      integer, intent(in) :: nnzp ! Grid points in the vertical [count]

      real, intent(in), dimension(nnzp) :: 
     .  gzt, gzm  ! Thermodynamic and momentum levels           [m]

      integer, intent(in) :: day, month, year  ! Time of year

      real, intent(in) :: 
     .  rlat, rlon   ! Latitude and Longitude             [Degrees N/E]

      real(kind=time_precision), intent(in) :: 
     .  delt         ! Timestep (dtmain in HOC)           [s]

      real(kind=time_precision), intent(in) :: 
     .  time_current ! Model time                         [s]

      ! Local Variables

      ! Namelist Variables

      character(len=10) :: stats_fmt  ! File storage convention

      character(len=20), dimension(nvarmax) :: 
     .  vars_zt,  ! Variables on the thermodynamic levels
     .  vars_zm,  ! Variables on the momentum levels
     .  vars_sfc  ! Variables at the model surface

      namelist /statsnl/
     .  vars_zt,
     .  vars_zm,
     .  vars_sfc

      ! Local Variables

      logical :: lerror

      character(len=200) :: fdir, fname

      integer :: i, k, ntot

      ! Initialize
      lerror = .false.

      ! Set stats_variables variables with inputs from calling subroutine
      lstats = lstats_in

      stats_tsamp = stats_tsamp_in
      stats_tsamp = stats_tsamp_in
      stats_tout  = stats_tout_in
      stats_fmt   = trim( stats_fmt_in )

      if ( .not. lstats ) then
        lstats_samp  = .false.
        lstats_first = .false.
        lstats_last  = .false.
        return
      end if

      ! Read namelist

      open(unit=iunit, file=fnamelist)
      read(unit=iunit, nml=statsnl, end=100)
      close(unit=iunit)


      ! Determine file names for GrADS or NetCDF files
      fname_zt  = trim( fname_prefix )//"_zt"
      fname_zm  = trim( fname_prefix )//"_zm"
      fname_sfc = trim( fname_prefix )//"_sfc"

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
      allocate( zt%in_update( zt%kk, zt%nn ) )
      call stats_zero( zt%kk, zt%nn, zt%x, zt%n, zt%in_update )

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
      allocate ( ztscr16(zt%kk) )

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
#ifdef NETCDF
        call open_nc( zt%f, iunit, fdir, fname, 
     .                1, zt%kk, zt%z,
     .                day, month, year, rlat, rlon,
     .                time_current+stats_tout, stats_tout,
     .                zt%nn )
#else
        stop "netCDF support was not compiled into this build."
#endif

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
      allocate( zm%in_update( zm%kk, zm%nn ) )
      
      call stats_zero( zm%kk, zm%nn, zm%x, zm%n, zm%in_update )

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
      allocate ( zmscr14(zm%kk) )


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
#ifdef NETCDF
        call open_nc( zm%f, iunit, fdir, fname, 
     .                1, zm%kk, zm%z,
     .                day, month, year, rlat, rlon,
     .                time_current+stats_tout, stats_tout,
     .                zm%nn )

#else
        stop "netCDF support was not compiled into this build."
#endif
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
      allocate( sfc%in_update( sfc%kk, sfc%nn ) )
      
      call stats_zero( sfc%kk, sfc%nn, sfc%x, sfc%n, sfc%in_update )

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
#ifdef NETCDF
        call open_nc( sfc%f, iunit, fdir, fname, 
     .                1, sfc%kk, sfc%z,
     .                day, month, year, rlat, rlon,
     .                time_current+stats_tout, stats_tout,
     .                sfc%nn )

#else
        stop "netCDF support was not compiled into this build."
#endif
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
      subroutine stats_zero( kk, nn, x, n, in_update )

!     Description:
!     Initialize stats to zero
!-----------------------------------------------------------------------
      use stats_precision, only:
     .    stat_rknd,  ! Variable(s)
     .    stat_nknd

      implicit none

      ! Input
      integer, intent(in) :: kk, nn

      ! Output
      real(kind=stat_rknd), dimension(kk,nn), intent(out)    :: x
      integer(kind=stat_nknd), dimension(kk,nn), intent(out) :: n
      logical, dimension(kk,nn), intent(out)              :: in_update

      ! Zero out arrays

      if ( nn > 0 ) then
        x(:,:) = 0.0
        n(:,:) = 0
        in_update(:,:) = .false.
      end if
      
      return
      end subroutine stats_zero

!-----------------------------------------------------------------------
      subroutine stats_avg( kk, nn, x, n )

!     Description:
!     Compute the average of stats fields
!-----------------------------------------------------------------------
      use stats_precision, only:
     .    stat_rknd,  ! Variable(s)
     .    stat_nknd

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

      use stats_variables, only:
     .    lstats, ! Variable(s)
     .    lstats_samp,
     .    lstats_first,
     .    lstats_last,
     .    stats_tsamp,
     .    stats_tout
      use stats_precision, only:
     .    time_precision ! Variable(s)

      implicit none

      ! Input

      real(kind=time_precision), intent(in) :: 
     .  time_elapsed ! Elapsed model time       [s]

      real(kind=time_precision), intent(in) :: 
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

      use stats_variables, only:
     .    zt, ! Variable(s)
     .    zm,
     .    sfc,
     .    lstats_last,
     .    stats_tsamp,
     .    stats_tout,
     .    lgrads

      use stats_precision, only:
     .    time_precision ! Variable(s)

      use outputgrads, only:
     .    write_grads ! Procedure(s)

#ifdef NETCDF
      use outputnc, only:
     .    write_nc ! Procedure(s)
#endif

      implicit none

      ! Input

      real(kind=time_precision), intent(in) :: 
     .  time_elapsed ! Time elapsed             [s]

      real(kind=time_precision), intent(in) :: 
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
           !pause
         end if
        end do
      end do
      
      do i=1,zm%nn
       do k=1,zm%kk
         if ( zm%n(k,i) /= 0 
     .        .and. zm%n(k,i) /= int(stats_tout/stats_tsamp) ) then
           write(0,*) 'Possible sampling error for zm ',i,k,zm%n(k,i)
           !pause
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
#ifdef NETCDF
        call write_nc( zt%f  )
        call write_nc( zm%f  )
        call write_nc( sfc%f  )
#else 
        stop "This program was not compiled with netCDF support"
#endif
      end if

      ! Reset sample fields
      call stats_zero( zt%kk, zt%nn, zt%x, zt%n, zt%in_update )
      call stats_zero( zm%kk, zm%nn, zm%x, zm%n, zt%in_update )
      call stats_zero( sfc%kk, sfc%nn, sfc%x, sfc%n, sfc%in_update )

      return
      end subroutine stats_end_timestep

!----------------------------------------------------------------------
      subroutine stats_accumulate
     .           ( um, vm, upwp, vpwp, up2, vp2, thlm,
     .             rtm, wprtp, wpthlp, wp2, wp3, rtp2, thlp2, rtpthlp,
     .             p, exner, rhot, rhom, thlm_forcing, rtm_forcing, wmm,
     .             wmt, Scm, taum, rcm, cf,
     .             sclrm, edsclrm, sclrm_forcing, wpsclrp )

!     Description:
!     Accumulate those stats variables that are preserved in HOC from
!     timestep to timestep, but not those stats that are not, (e.g. 
!     budget terms, longwave and shortwave components, etc. )
!----------------------------------------------------------------------

      use stats_variables, only:
     .    zt,     ! Variables
     .    zm,
     .    sfc,
     .    lstats_samp,
     .    ithlm,
     .    itemp_in_K,
     .    ithvm,
     .    irtm,
     .    ircm,
     .    ium,
     .    ivm,
     .    iwmt,
     .    iug,
     .    ivg,
     .    icf,
     .    ip,
     .    iexner,
     .    iLscale,
     .    iwp3,
     .    iwpthlp2,
     .    iwp2thlp, 
     .    iwprtp2,
     .    iwp2rtp,
     .    iLup,
     .    iLdown,
     .    itaut,
     .    iKht,
     .    iwp2thvp,
     .    iwp2rcp,
     .    iwprtpthlp,
     .    isct,         
     .    irhot,
     .    irsm,
     .    iAKm,
     .    iAKm_est,
     .    iradht,
     .    ia,
     .    iw1,
     .    iw2,
     .    isw1,
     .    isw2,
     .    ithl1,
     .    ithl2,
     .    isthl1,
     .    isthl2,
     .    irt1,
     .    irt2,
     .    isrt1,
     .    isrt2,
     .    irc1,
     .    irc2,
     .    irsl1,
     .    irsl2,
     .    iR1,
     .    iR2,
     .    is1,
     .    is2,
     .    iss1,
     .    iss2,
     .    irrtthl

      use stats_variables, only:
     .    iwp2zt,
     .    ithlp2zt,
     .    iwpthlpzt,
     .    iwprtpzt,
     .    irtp2zt,
     .    irtpthlpzt,
     .    iwp2,
     .    irtp2,
     .    ithlp2,
     .    irtpthlp,
     .    iwprtp, 
     .    iwpthlp,
     .    iwp4, 
     .    iwpthvp,
     .    irtpthvp,
     .    ithlpthvp,
     .    itaum,
     .    iKhm,
     .    iwprcp,
     .    ithlprcp,
     .    irtprcp,
     .    ircp2,
     .    iupwp,
     .    ivpwp,
     .    iup2,
     .    ivp2,
     .    irhom,
     .    iscm,
     .    iem,
     .    ishear,
     .    iFrad,
     .    icc,
     .    izb,
     .    ilwp

#ifdef SCALARS
      use stats_variables, only:
     .    isclram,
     .    isclram_f,
     .    isclrbm,
     .    isclrbm_f,
     .    iedsclram,
     .    iedsclrbm,
     .    isclraprtp,
     .    isclrbprtp,
     .    isclrap2,
     .    isclrbp2,
     .    isclrapthvp,
     .    isclrbpthvp,
     .    isclrapthlp,
     .    isclrbpthlp,
     .    isclraprcp,
     .    isclrbprcp,
     .    iwpsclrap,
     .    iwpsclrbp,
     .    iwp2sclrap,
     .    iwp2sclrbp,
     .    iwpsclrap2,
     .    iwpsclrbp2,
     .    iwpsclraprtp,
     .    iwpsclrbprtp,
     .    iwpsclrapthlp,
     .    iwpsclrbpthlp,
     .    iwpedsclrap,
     .    iwpedsclrbp
#endif /*SCALARS*/

      use grid_class, only:
     .    gr ! Variable
      !use prognostic_variables
      use diagnostic_variables, only:
     .    pdf_parms, ! Variable(s)
     .    thvm,
     .    ug,
     .    vg,
     .    Lscale,
     .    wpthlp2,
     .    wp2thlp,
     .    wprtp2,
     .    wp2rtp,
     .    lup,
     .    ldown,
     .    taut,
     .    Kht,
     .    wp2thvp,
     .    wp2rcp,
     .    wprtpthlp,
     .    sct,
     .    rsm,
     .    Akm,
     .    Akm_est,
     .    radht,
     .    wp2zt,
     .    thlp2zt,
     .    wpthlpzt,
     .    wprtpzt,
     .    rtp2zt,
     .    rtpthlpzt,
     .    wp4,
     .    wpthvp,
     .    rtpthvp,
     .    thlpthvp,
     .    khm,
     .    wprcp,
     .    thlprcp,
     .    rtprcp,
     .    rcp2,
     .    em,
     .    shear,
     .    Frad
#ifdef SCALARS
     .   ,sclrprtp,
     .    sclrp2,
     .    sclrpthvp,
     .    sclrpthlp,
     .    sclrprcp,
     .    wp2sclrp,
     .    wpsclrp2,
     .    wpsclrprtp,
     .    wpsclrpthlp,
     .    wpedsclrp   
#endif /*SCALARS*/

      
      use model_flags, only:
     .    LH_on ! Variable(s)

      use temp_in_K_mod, only:
     .    thlm2temp_in_K ! Procedure

      use constants, only:
     .    sclr_dim, ! Variable(s)
     .    rc_tol

      use stats_type, only:
     .    stat_update_var, ! Procedure(s)
     ,    stat_update_var_pt

      use interpolation, only:
     .    linint ! Procedure

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
     .  sclrm_forcing,  ! Large-scale forcing of scalar [units/s]
     .  wpsclrp         ! w'sclr'                       [units m/s]

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
        call stat_update_var( ithlm, thlm, zt )
        call stat_update_var( itemp_in_K, 
     .         thlm2temp_in_K( thlm, exner, rcm), zt )
        call stat_update_var( ithvm, thvm, zt )
        call stat_update_var( irtm, rtm, zt )       
        call stat_update_var( ircm, rcm, zt )
        call stat_update_var( ium, um, zt )
        call stat_update_var( ivm, vm, zt )
        call stat_update_var( iwmt, wmt, zt )
        call stat_update_var( iug, ug, zt )
        call stat_update_var( ivg, vg, zt )
        call stat_update_var( icf, cf, zt )
        call stat_update_var( ip, p, zt )
        call stat_update_var( iexner, exner, zt )
        call stat_update_var( iLscale, Lscale, zt )
        call stat_update_var( iwp3, wp3, zt )
        call stat_update_var( iwpthlp2, wpthlp2, zt )
        call stat_update_var( iwp2thlp, wp2thlp, zt )
        call stat_update_var( iwprtp2, wprtp2, zt )
        call stat_update_var( iwp2rtp, wp2rtp, zt )
        call stat_update_var( iLup, Lup, zt )
        call stat_update_var( iLdown, Ldown, zt )
        call stat_update_var( itaut, taut, zt )
        call stat_update_var( iKht, Kht, zt )
        call stat_update_var( iwp2thvp, wp2thvp, zt )
        call stat_update_var( iwp2rcp, wp2rcp, zt )
        call stat_update_var( iwprtpthlp, wprtpthlp, zt )
        call stat_update_var( isct, sct, zt )
        call stat_update_var( irhot, rhot, zt )
!        call stat_update_var( iNcm, Ncm, zt )
!        call stat_update_var( iNcnm, Ncnm, zt )
!        call stat_update_var( iNim, Nim, zt )
!       if ( cloud_sed ) then
!        call stat_update_var( ised_rcm, sed_rcm, zt )
!       endif
        call stat_update_var( irsm, rsm, zt )
!        call stat_update_var( irrm, hydromet(:,1), zt )
!        call stat_update_var( iNrm, hydromet(:,2), zt )
!        call stat_update_var( irsnowm, hydromet(:,3), zt )
!        call stat_update_var( iricem, hydromet(:,4), zt )
!        call stat_update_var( irgraupelm, hydromet(:,5), zt )

        if ( LH_on ) then
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
        call stat_update_var( iwp2zt, wp2zt, zt )
        call stat_update_var( ithlp2zt, thlp2zt, zt )
        call stat_update_var( iwpthlpzt, wpthlpzt, zt )
        call stat_update_var( iwprtpzt, wprtpzt, zt )
        call stat_update_var( irtp2zt, rtp2zt, zt )
        call stat_update_var( irtpthlpzt, rtpthlpzt, zt )

#ifdef SCALARS
        call stat_update_var( isclram, sclrm(:,1), zt )
        call stat_update_var( isclram_f, sclrm_forcing(:,1),  zt )
        call stat_update_var( isclrbm, sclrm(:,2), zt )
        call stat_update_var( isclrbm_f, sclrm_forcing(:,2), zt )
        call stat_update_var( iedsclram, edsclrm(:,1), zt )
        call stat_update_var( iedsclrbm, edsclrm(:,2), zt )


#endif /*SCALARS*/

!       zm variables

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
        call stat_update_var( itaum, taum, zm )
        call stat_update_var( iKhm, Khm, zm )
        call stat_update_var( iwprcp, wprcp, zm )
        call stat_update_var( ithlprcp, thlprcp, zm )
        call stat_update_var( irtprcp, rtprcp, zm )
        call stat_update_var( ircp2, rcp2, zm )
        call stat_update_var( iupwp, upwp, zm )
        call stat_update_var( ivpwp, vpwp, zm )
        call stat_update_var( ivp2, vp2, zm )
        call stat_update_var( iup2, up2, zm )
        call stat_update_var( irhom, rhom, zm )
        call stat_update_var( iscm, scm, zm )
        call stat_update_var( iem, em, zm )
        call stat_update_var( ishear, shear, zm )
        call stat_update_var( iFrad, Frad, zm )
!        if ( cloud_sed ) then
!          call stat_update_var( iFcsed, Fcsed, zm )
!        endif

#ifdef SCALARS
        call stat_update_var( isclraprtp,sclrprtp(:,1), zm )
        call stat_update_var( isclrbprtp, sclrprtp(:,2), zm )
        call stat_update_var( isclrap2, sclrp2(:,1), zm )
        call stat_update_var( isclrbp2, sclrp2(:,2), zm )
        call stat_update_var( isclrapthvp, sclrpthvp(:,1), zm )
        call stat_update_var( isclrbpthvp, sclrpthvp(:,2), zm )
        call stat_update_var( isclrapthlp, sclrpthlp(:,1), zm )
        call stat_update_var( isclrbpthlp, sclrpthlp(:,2), zm )
        call stat_update_var( isclraprcp, sclrprcp(:,1), zm ) 
        call stat_update_var( isclrbprcp, sclrprcp(:,2), zm )
        call stat_update_var( iwpsclrap, wpsclrp(:,1), zm )
        call stat_update_var( iwpsclrbp, wpsclrp(:,2), zm )
        call stat_update_var( iwp2sclrap, wp2sclrp(:,1), zm )
        call stat_update_var( iwp2sclrbp, wp2sclrp(:,2), zm )
        call stat_update_var( iwpsclrap2, wpsclrp2(:,1), zm )
        call stat_update_var( iwpsclrbp2, wpsclrp2(:,2), zm )
        call stat_update_var( iwpsclraprtp, wpsclrprtp(:,1), zm )
        call stat_update_var( iwpsclrbprtp, wpsclrprtp(:,2), zm )
        call stat_update_var( iwpsclrapthlp, wpsclrpthlp(:,1), zm )
        call stat_update_var( iwpsclrbpthlp, wpsclrpthlp(:,2), zm )
        call stat_update_var( iwpedsclrap, wpedsclrp(:,1), zm )
        call stat_update_var( iwpedsclrbp, wpedsclrp(:,2), zm )
        
#endif /*SCALARS*/

! sfc variables

!       Cloud cover
        call stat_update_var_pt( icc, 1, maxval( cf(1:gr%nnzp) ), sfc )

!       Cloud base
        if ( izb > 0 ) then
          k = 1
          do while ( rcm(k) < rc_tol .and. k < gr%nnzp )
            k = k + 1
          end do

!          if ( k < gr%nnzp) then
!            sfc%x(1,izb) = sfc%x(1,izb) + gr%zt(k)

!        pure real function linint( zmid, ztop, zbot, vtop, vbot )
!        linint = ( (zmid-zbot)/(ztop-zbot) ) * (vtop-vbot) + vbot

          if ( k > 1 .AND. k < gr%nnzp) then
            call stat_update_var_pt( izb, 1, linint( rc_tol,rcm(k), ! Use linear interpolation
     .               rcm(k-1), gr%zt(k) , gr%zt(k-1) ), sfc )      ! to find the exact height
                                                                   ! of the rc_tol kg/kg level.
                                                                   ! Brian.
          else
            call stat_update_var_pt( izb, 1, -10.0 , sfc )

          end if

        end if

!       LWP
        if ( ilwp > 0 ) then
          xtmp = 0.
          do i = gr%nnzp-1, 1, -1
           xtmp = xtmp + rhot(i+1) * rcm(i+1) / gr%dzt(i+1)
          end do
          
          call stat_update_var_pt( ilwp, 1, xtmp, sfc )

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

      use stats_variables, only:
     .    zt, ! Variable(s)
     .    zm,
     .    sfc,
     .    ztscr01,
     .    ztscr02,
     .    ztscr03,
     .    ztscr04,
     .    ztscr05,
     .    ztscr06,
     .    ztscr07,
     .    ztscr08,
     .    ztscr09,
     .    ztscr10,
     .    ztscr11,
     .    ztscr12,
     .    ztscr13,
     .    ztscr14,
     .    ztscr15,
     .    ztscr16,
     .    zmscr01,
     .    zmscr02,
     .    zmscr03,
     .    zmscr04,
     .    zmscr05,
     .    zmscr06,
     .    zmscr07,
     .    zmscr08,
     .    zmscr09,
     .    zmscr10,
     .    zmscr11,
     .    zmscr12,
     .    zmscr13,
     .    zmscr14,
     .    lnetcdf,
     .    lstats
#ifdef NETCDF
      use outputnc, only: 
     .    close_nc ! Procedure
#endif

      implicit none

      integer :: i

      if ( lstats .and. lnetcdf ) then
#ifdef NETCDF
        call close_nc( zt%f )
        call close_nc( zm%f )
        call close_nc( sfc%f )
#else
        stop "This program was not compiled with netCDF support"
#endif
      end if

      if ( lstats ) then
        ! De-allocate all zt variables
        deallocate( zt%z )

        deallocate( zt%x, stat=i )
        deallocate( zt%n )
        deallocate( zt%in_update )

        
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

        ! De-allocate all zm variables
        deallocate( zm%z )

        deallocate( zm%x )
        deallocate( zm%n )

        deallocate( zm%f%var )
        deallocate( zm%f%z )
        deallocate( zm%in_update )
        
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

        ! De-allocate all sfc variables
        deallocate( sfc%z )

        deallocate( sfc%x )
        deallocate( sfc%n )
        deallocate( sfc%in_update )

        deallocate( sfc%f%var )
        deallocate( sfc%f%z )
      end if ! lstats

      return
      end subroutine stats_finalize

#endif /*STATS*/

      end module stats_subs
