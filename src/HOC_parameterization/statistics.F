!-----------------------------------------------------------------------
!  $Id: statistics.F,v 1.4 2006-05-12 16:55:12 dschanen Exp $

!  module stats_type
!  contains:
!    derived data type 'stats'
!-----------------------------------------------------------------------
#ifdef STATS
      module stats_type

        use outputfile_class

c       Derived data types to store GrADS statistics

        type stats 

c         Number of fields to sample

          integer nn

c         Vertical extend of variable

          integer kk

c         Vertical levels

          real,    pointer :: z (:)

c         Array to store sampled fields

          real,    pointer :: x (:,:)
          integer, pointer :: n (:,:)

c         Data for GrADS output

          type (outputfile) f

        end type stats

      end module stats_type

!-----------------------------------------------------------------------
!  module statistics

!  holds pointers to variables to be written to GrADS files
!-----------------------------------------------------------------------
      module statistics

        use stats_type

c       Sampling and output frequencies

        logical lstats         ! Main flag to turn statistics on/off
        logical lnetcdf        ! Output to netcdf format
        logical lgrads         ! Output to GrADS format
        real stats_tsamp       ! Sampling interval (s)
        real stats_tout        ! Output interval (s)

        logical lstats_samp    ! Sample flag for current time step
        logical lstats_first   ! First time step of output period
        logical lstats_last    ! Last time step of output period

        character(len=200) :: fname_zt
        character(len=200) :: fname_zm
        character(len=200) :: fname_sfc

c       Indices for statistics in zt file

        integer
     .     ithlm,
     .     ithvm,
     .     irtm,
     .     ircm,
     .     ium,
     .     ivm,
     .     iwmt,
     .     iug,
     .     ivg,
     .     icf,
     .     ip,
     .     iLscale,
     .     ithlm_forcing,
     .     irtm_forcing,
     .     iwp3,
     .     iwpthlp2,
     .     iwp2thlp,
     .     iwprtp2,
     .     iwp2rtp,
     .     ilup,
     .     ildown,
     .     itaut,
     .     ikht,
     .     iwp2thvp,
     .     iwp2rcp,
     .     iwprtpthlp,
     .     isct,
     .     irhot,
     .     iNcm,           ! Brian
     .     ised_rcm,       ! Brian
     .     irsm,           ! Brian
     .     irrm,           ! Brian
     .     iNrm,           ! Brian
     .     iVrr,           ! Brian
     .     imean_vol_rad,  ! Brian
     .     irain_rate,     ! Brian
     .     iAKm,           ! analytic Kessler.  Vince Larson 22 May 2005 
     .     iAKm_est,       ! LH Kessler.  Vince Larson  22 May 2005
     .     iradht,         ! Radiative heating.  Dave Schanen 7 Oct 2005
     .     iradht_LW,      !   "           "   Long-wave
     .     iradht_SW       !   "           "   Short-wave


        integer
     .     iwp3_bt,
     .     iwp3_ma,
     .     iwp3_ta,
     .     iwp3_tp,
     .     iwp3_ac,
     .     iwp3_bp,
     .     iwp3_pr1,
     .     iwp3_pr2,
     .     iwp3_dp1,
     .     iwp3_cl

        integer
     .     irrm_bt,
     .     irrm_ma,
     .     irrm_sd,
     .     irrm_dff,
     .     irrm_cond,
     .     irrm_auto,
     .     irrm_accr,
     .     iNrm_bt,
     .     iNrm_ma,
     .     iNrm_sd,
     .     iNrm_dff,
     .     iNrm_cond,
     .     iNrm_auto

        integer
     .     ia,
     .     iw1,
     .     iw2,
     .     isw1,
     .     isw2,
     .     ithl1,
     .     ithl2,
     .     isthl1,
     .     isthl2,
     .     irt1,
     .     irt2,
     .     isrt1,
     .     isrt2,
     .     irc1,
     .     irc2,
     .     irsl1,
     .     irsl2,
     .     iR1,
     .     iR2,
     .     is1,
     .     is2,
     .     iss1,
     .     iss2,
     .     irrtthl 

        integer
     .     iwp2zt,
     .     ithlp2zt,
     .     iwpthlpzt,
     .     iwprtpzt,
     .     irtp2zt,
     .     irtpthlpzt

#ifdef SCALARS
        integer
     .     isclram,          ! New mixing scheme ( mean rt )
     .     isclrbm,          ! New mixing scheme ( mean th_l )
     .     isclram_f,        ! New mixing scheme ( mean rt forcing)
     .     isclrbm_f         ! New mixing scheme ( mean th_l forcing)

        integer
     .     iedsclram,         ! Eddy-diff. term ( rt )
     .     iedsclrbm          ! Eddy-diff. term ( theta_l )
#endif /*SCALARS*/

c       Indices for statistics in zm file

        integer
     .     iwp2,
     .     irtp2,
     .     ithlp2,
     .     irtpthlp,
     .     iwprtp,
     .     iwpthlp,
     .     iwp4,
     .     iwpthvp,
     .     irtpthvp,
     .     ithlpthvp,
     .     itaum,
     .     ikhm,
     .     iwprcp,
     .     ithlprcp,
     .     irtprcp,
     .     iupwp,
     .     ivpwp,
     .     irhom,
     .     iscm,
     .     iem,
     .     ishear,    ! Brian
     .     iFrad,
     .     iFrad_LW,  ! Brian
     .     iFrad_SW,  ! Brian
     .     iFprec,    ! Brian
     .     iFcsed     ! Brian

        integer
     .     iwp2_bt,
     .     iwp2_ma,
     .     iwp2_ta,
     .     iwp2_ac,
     .     iwp2_bp,
     .     iwp2_pr1,
     .     iwp2_pr2,
     .     iwp2_pr3,
     .     iwp2_dp1,
     .     iwp2_dp2,
     .     iwp2_cl

!    Changed qt to rt to avoid confusion -dschanen 2 Aug 05
        integer
     .     iwprtp_bt,
     .     iwprtp_ma,
     .     iwprtp_ta,
     .     iwprtp_tp,
     .     iwprtp_ac,
     .     iwprtp_bp,
     .     iwprtp_pr1,
     .     iwprtp_pr2,
     .     iwprtp_pr3,
     .     iwprtp_dp1

        integer
     .     iwpthlp_bt,
     .     iwpthlp_ma,
     .     iwpthlp_ta,
     .     iwpthlp_tp,
     .     iwpthlp_ac,
     .     iwpthlp_bp,
     .     iwpthlp_pr1,
     .     iwpthlp_pr2,
     .     iwpthlp_pr3,
     .     iwpthlp_dp1

!    Dr. Golaz's new budget terms
!    qt was changed to rt to avoid confusion

        integer
     .     irtp2_bt,
     .     irtp2_ma,
     .     irtp2_ta,
     .     irtp2_tp,
     .     irtp2_dp1,
     .     irtp2_dp2,
     .     irtp2_cl

        integer
     .     ithlp2_bt,
     .     ithlp2_ma,
     .     ithlp2_ta,
     .     ithlp2_tp,
     .     ithlp2_dp1,
     .     ithlp2_dp2,
     .     ithlp2_cl

        integer
     .     irtpthlp_bt,
     .     irtpthlp_ma,
     .     irtpthlp_ta,
     .     irtpthlp_tp1,
     .     irtpthlp_tp2,
     .     irtpthlp_dp1,
     .     irtpthlp_dp2

#ifdef SCALARS
!       New mixing scheme.  Note that floating point roundoff may make 
!       mathematically equivalent variables different.
        integer
     .     isclraprtp,          ! sclr'(1)rt'     / rt'^2
     .     isclrbprtp,          ! sclr'(2)rt'     / thl'rt'
     .     isclrap2,            ! sclr'(1)^2      / rt'^2
     .     isclrbp2,            ! sclr'(2)^2      / thl'^2
     .     isclrapthvp,         ! sclr'(1)th_v'   / rt'th_v' 
     .     isclrbpthvp,         ! sclr'(2)th_v'   / th_l' th_v'
     .     isclrapthlp,         ! sclr'(1)th_l'   / rt'th_l' 
     .     isclrbpthlp,         ! sclr'(2)th_l'   / th_l'^2
     .     isclraprcp,          ! sclr'(1)rc'     / rt'rc'
     .     isclrbprcp,          ! sclr'(2)rc'     / th_l'rc'
     .     iwpsclrap,           ! w'slcr'(1)      / w'rt'
     .     iwpsclrbp,           ! w'sclr'(2)      / w'th_l'
     .     iwp2sclrap,          ! w'^2 sclr'(1)   / w'^2 rt'
     .     iwp2sclrbp,          ! w'^2 sclr'(2)   / w'^2 th_l'
     .     iwpsclrap2,          ! w'sclr'(1)^2    / w'rt'^2
     .     iwpsclrbp2,          ! w'sclr'(2)^2    / w'th_l'^2
     .     iwpsclraprtp,        ! w'sclr'(1)rt'   / w'rt'^2
     .     iwpsclrbprtp,        ! w'sclr'(2)rt'   / w'thl'rt'
     .     iwpsclrapthlp,       ! w'sclr'(1)th_l' / w'rt'th_l' 
     .     iwpsclrbpthlp        ! w'sclr'(2)th_l' / w'th_l'^2

        integer
     .     iwpedsclrap,  ! eddy sclr'(1)w'
     .     iwpedsclrbp   ! eddy sclr'(2)w'
#endif /*SCALARS*/

c       Indices for statistics in sfc file

        integer
     .     iustar,
     .     ilh,
     .     ish,
     .     icc,
     .     ilwp,
     .     izb,
     .     izi,
     .     irain,   ! Brian
     .     ipflux   ! Brian
  
c       Variables that contains all the statistics

        type (stats), target :: zt,  ! zt grid
     .                          zm,  ! zm grid
     .                          sfc  ! sfc


c       Scratch space

        real, allocatable :: ztscr01(:), ztscr02(:), ztscr03(:), 
     .                       ztscr04(:), ztscr05(:), ztscr06(:), 
     .                       ztscr07(:), ztscr08(:), ztscr09(:),
     .                       ztscr10(:), ztscr11(:), ztscr12(:),
     .                       ztscr13(:), ztscr14(:), ztscr15(:)

        real, allocatable :: zmscr01(:), zmscr02(:), zmscr03(:), 
     .                       zmscr04(:), zmscr05(:), zmscr06(:),
     .                       zmscr07(:), zmscr08(:), zmscr09(:),
     .                       zmscr10(:), zmscr11(:), zmscr12(:),
     .                       zmscr13(:)

      end module statistics

c     -----------------------------------------------------------------
      subroutine stats_init( iunit, fnamelist, nnzp, gzt, gzm,
     .                       day, month, year, rlat, rlon, time, delt )

      use statistics
      use outputgrads
      use outputnc

      implicit none

!     Input

      integer, intent(in) :: iunit
      character(len=*)    :: fnamelist
      integer, intent(in) :: nnzp
      real,    intent(in) :: gzt(nnzp), gzm(nnzp)
      integer, intent(in) :: day, month, year
      real, intent(in)    :: rlat, rlon
      real, intent(in)    :: time
      real, intent(in)    :: delt

!     Namelist

      integer, parameter :: nvarmax = 250
      character(len=10) stats_fmt
      character(len=20) vars_zt  (nvarmax)
      character(len=20) vars_zm  (nvarmax)
      character(len=20) vars_sfc (nvarmax)

      namelist /statsnl/
     .  lstats, stats_fmt, stats_tsamp, stats_tout,
     .  fname_zt,  vars_zt,
     .  fname_zm,  vars_zm,
     .  fname_sfc, vars_sfc

!     Internal

      logical lerror

      character(len=200) :: fdir, fname

      integer :: i, k, ntot

!     Read namelist

      lerror = .false.
      open(unit=iunit, file=fnamelist)
      read(iunit,statsnl,end=100)
      close(iunit)

      if ( .not. lstats ) then
        lstats_samp  = .false.
        lstats_first = .false.
        lstats_last  = .false.
        return
      endif

      if ( trim(stats_fmt) =='grads'.or. 
     .     trim(stats_fmt) =='gr'  ) then
        lnetcdf = .false.
        lgrads  = .true.
      else if ( trim( stats_fmt ) == 'netcdf'.or.
     .          trim( stats_fmt ) == 'nc' ) then
        lnetcdf = .true.
        lgrads  = .false.
      else
        write(*,*) "invalid data format "//trim( stats_fmt )
        stop
      endif

!     Check sampling and output frequencies

      if ( abs( stats_tsamp/delt - int(stats_tsamp/delt) ) 
     .     > 1.e-8 ) then
         lerror = .true.
         write(*,*) 'Error: stats_tsamp should be a multiple of delt'
         write(*,*) 'stats_tsamp = ',stats_tsamp
         write(*,*) 'delt = ',delt
      endif

      if ( abs( stats_tout/stats_tsamp - int(stats_tout/stats_tsamp) ) 
     .     > 1.e-8 ) then
         lerror = .true.
         write(*,*) 
     .      'Error: stats_tout should be a multiple of stats_tsamp'
         write(*,*) 'stats_tout = ',stats_tout
         write(*,*) 'stats_tsamp = ',stats_tsamp
      endif

!     Initialize zt (mass points)

      i = 1
      do while ( ichar(vars_zt(i)(1:1)) /= 0 
     .           .and. len_trim(vars_zt(i)) /= 0
     .           .and. i <= nvarmax )
        i = i + 1
      end do
      ntot = i - 1
      if ( ntot == nvarmax )
     .  write(*,*) 'WARNING: check nvarmax in statistics.f'

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

!     Allocate scratch space

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

!     Open GrADS files
      fdir = "./"
      fname = trim(fname_zt)
      if ( lgrads ) then
        call open_grads( zt%f, 20, fdir, fname, 
     .                   1, zt%kk, zt%z,
     .                   day, month, year, rlat, rlon,
     .                   time+stats_tout, stats_tout,
     .                   zt%nn )
      else ! lnetcdf
        call open_nc( zt%f, 20, fdir, fname, 
     .                1, zt%kk, zt%z,
     .                day, month, year, rlat, rlon,
     .                time+stats_tout, stats_tout,
     .                zt%nn )

      endif

!     Default initialization for array indices for zt

      ithlm         = 0
      ithvm         = 0
      irtm          = 0
      ircm          = 0
      ium           = 0
      ivm           = 0
      iwmt          = 0
      iug           = 0
      ivg           = 0
      icf           = 0
      ip            = 0
      iLscale       = 0
      ithlm_forcing = 0
      irtm_forcing  = 0
      iwp3          = 0
      iwpthlp2      = 0
      iwp2thlp      = 0
      iwprtp2       = 0
      iwp2rtp       = 0
      ilup          = 0
      ildown        = 0
      itaut         = 0
      ikht          = 0
      iwp2thvp      = 0
      iwp2rcp       = 0
      iwprtpthlp    = 0
      isct          = 0
      irhot         = 0
      iNcm          = 0  ! Brian
      ised_rcm      = 0  ! Brian
      irsm          = 0  ! Brian
      irrm          = 0  ! Brian
      iNrm          = 0  ! Brian
      iVrr          = 0  ! Brian
      imean_vol_rad = 0  ! Brian
      irain_rate    = 0  ! Brian
      iAKm          = 0  ! analytic Kessler.  Vince Larson 22 May 2005
      iAKm_est      = 0  ! LH Kessler.  Vince Larson 22 May 2005
      iradht        = 0
      iradht_LW     = 0
      iradht_SW     = 0

      iwp3_bt       = 0
      iwp3_ma       = 0
      iwp3_ta       = 0
      iwp3_tp       = 0
      iwp3_ac       = 0
      iwp3_bp       = 0
      iwp3_pr1      = 0
      iwp3_pr2      = 0
      iwp3_dp1      = 0
      iwp3_cl       = 0

      irrm_bt       = 0
      irrm_ma       = 0
      irrm_sd       = 0
      irrm_dff      = 0
      irrm_cond     = 0
      irrm_auto     = 0
      irrm_accr     = 0
      iNrm_bt       = 0
      iNrm_ma       = 0
      iNrm_sd       = 0
      iNrm_dff      = 0
      iNrm_cond     = 0
      iNrm_auto     = 0

      ia            = 0
      iw1           = 0
      iw2           = 0
      isw1          = 0
      isw2          = 0
      ithl1         = 0
      ithl2         = 0
      isthl1        = 0
      isthl2        = 0
      irt1          = 0
      irt2          = 0
      isrt1         = 0
      isrt2         = 0
      irc1          = 0
      irc2          = 0
      irsl1         = 0
      irsl2         = 0
      iR1           = 0
      iR2           = 0
      is1           = 0
      is2           = 0
      iss1          = 0
      iss2          = 0
      irrtthl       = 0

      iwp2zt     = 0
      ithlp2zt   = 0
      iwpthlpzt  = 0
      iwprtpzt   = 0
      irtp2zt    = 0
      irtpthlpzt = 0

#ifdef SCALARS
      isclram     = 0    ! new mixing scheme
      isclram_f   = 0    !  "             "
      isclrbm     = 0    ! new mixing scheme
      isclrbm_f   = 0    !  "             "

      iedsclram   = 0
      iedsclrbm   = 0
#endif /*SCALARS*/

!     Assign pointers for statistics variables zt

      k = 1
      do i=1,zt%nn

        select case ( trim(vars_zt(i)) )
        case ('thlm')
          ithlm = k
          zt%f%var(ithlm)%ptr => zt%x(:,k)
          zt%f%var(ithlm)%name = "thlm"
          zt%f%var(ithlm)%description = "thetal (K)"
          zt%f%var(ithlm)%units = "K"
          k = k + 1

        case ('thvm')
          ithvm = k
          zt%f%var(ithvm)%ptr => zt%x(:,k)
          zt%f%var(ithvm)%name = "thvm"
          zt%f%var(ithvm)%description 
     .    = "virtual potential temperature (K)"
          zt%f%var(ithvm)%units = "K"
          k = k + 1

        case ('rtm')
          irtm = k
          zt%f%var(irtm)%ptr => zt%x(:,k)
          zt%f%var(irtm)%name = "rtm"
          zt%f%var(irtm)%description 
     .    = "total water mixing ratio (kg/kg)"
          zt%f%var(irtm)%units = "kg/kg"
          k = k + 1

        case ('rcm')
          ircm = k
          zt%f%var(ircm)%ptr => zt%x(:,k)
          zt%f%var(ircm)%name = "rcm"
          zt%f%var(ircm)%description 
     .    = "liquid water mixing ratio (kg/kg)"
          zt%f%var(ircm)%units = "kg/kg"
          k = k + 1

        case ('um')
          ium = k
          zt%f%var(ium)%ptr => zt%x(:,k)
          zt%f%var(ium)%name = "um"
          zt%f%var(ium)%description = "u wind (m/s)"
          zt%f%var(ium)%units = "m/s"
          k = k + 1

        case ('vm')
          ivm = k
          zt%f%var(ivm)%ptr => zt%x(:,k)
          zt%f%var(ivm)%name = "vm"
          zt%f%var(ivm)%description = "v wind (m/s)"
          zt%f%var(ivm)%units = "m/s"
          k = k + 1

        case ('wmt')
          iwmt = k
          zt%f%var(iwmt)%ptr => zt%x(:,k)
          zt%f%var(iwmt)%name = "wm"
          zt%f%var(iwmt)%description = "w wind (m/s)"
          zt%f%var(iwmt)%units = "m/s"
          k = k + 1

        case ('ug')
          iug = k
          zt%f%var(iug)%ptr => zt%x(:,k)
          zt%f%var(iug)%name = "ug"
          zt%f%var(iug)%description = "u geostrophic wind (m/s)"
          zt%f%var(iug)%units = "m/s"
          k = k + 1

        case ('vg')
          ivg = k
          zt%f%var(ivg)%ptr => zt%x(:,k)
          zt%f%var(ivg)%name = "vg"
          zt%f%var(ivg)%description = "v geostrophic wind (m/s)"
          zt%f%var(ivg)%units = "m/s"
          k = k + 1

        case ('cf')
          icf = k
          zt%f%var(icf)%ptr => zt%x(:,k)
          zt%f%var(icf)%name = "cf"
          zt%f%var(icf)%description = "cloud fraction"
          zt%f%var(icf)%units = "count"
          k = k + 1

        case ('p')
          ip = k
          zt%f%var(ip)%ptr => zt%x(:,k)
          zt%f%var(ip)%name = "p"
          zt%f%var(ip)%description = "pressure (Pa)"
          zt%f%var(ip)%units = "Pa"
          k = k + 1

        case ('Lscale')
          iLscale = k
          zt%f%var(iLscale)%ptr => zt%x(:,k)
          zt%f%var(iLscale)%name = "Lscale"
          zt%f%var(iLscale)%description = "mixing length"
          zt%f%var(iLscale)%units = "m"
          k = k + 1

        case ('thlm_forcing')
          ithlm_forcing = k
          zt%f%var(ithlm_forcing)%ptr => zt%x(:,k)
          zt%f%var(ithlm_forcing)%name = "thlm_f"
          zt%f%var(ithlm_forcing)%description = "thetal forcing"
          zt%f%var(ithlm_forcing)%units = "K/s"
          k = k + 1

        case ('rtm_forcing')
          irtm_forcing = k
          zt%f%var(irtm_forcing)%ptr => zt%x(:,k)
          zt%f%var(irtm_forcing)%name = "rtm_f"
          zt%f%var(irtm_forcing)%description = "rt forcing"
          zt%f%var(irtm_forcing)%units = "kg/(kg s)"
          k = k + 1

        case ('wp3')
          iwp3 = k
          zt%f%var(iwp3)%ptr => zt%x(:,k)
          zt%f%var(iwp3)%name = "wp3"
          zt%f%var(iwp3)%description = "w third order moment"
          zt%f%var(iwp3)%units = "(m^3)/(s^3)"
          k = k + 1

        case ('wpthlp2')
          iwpthlp2 = k
          zt%f%var(iwpthlp2)%ptr => zt%x(:,k)
          zt%f%var(iwpthlp2)%name = "wpthlp2"
          zt%f%var(iwpthlp2)%description = "wpthlp2 covariance"
          zt%f%var(iwpthlp2)%units = "(m K^2)/s"
          k = k + 1

        case ('wp2thlp')
          iwp2thlp = k
          zt%f%var(iwp2thlp)%ptr => zt%x(:,k)
          zt%f%var(iwp2thlp)%name = "wp2thlp"
          zt%f%var(iwp2thlp)%description = "wp2thlp covariance"
          zt%f%var(iwp2thlp)%units = "(m^2 K)/(s^2)"
          k = k + 1

        case ('wprtp2')
          iwprtp2 = k
          zt%f%var(iwprtp2)%ptr => zt%x(:,k)
          zt%f%var(iwprtp2)%name = "wprtp2"
          zt%f%var(iwprtp2)%description = "wprtp2 covariance"
          zt%f%var(iwprtp2)%units = "(m kg^2)/(s kg^2)"
          k = k + 1

        case ('wp2rtp')
          iwp2rtp = k
          zt%f%var(iwp2rtp)%ptr => zt%x(:,k)
          zt%f%var(iwp2rtp)%name = "wp2rtp"
          zt%f%var(iwp2rtp)%description = "wp2rtp covariance"
          zt%f%var(iwp2rtp)%units = "(m2 kg)/(s2 kg)"
          k = k + 1

        case ('lup')
          ilup = k
          zt%f%var(ilup)%ptr => zt%x(:,k)
          zt%f%var(ilup)%name = "lup"
          zt%f%var(ilup)%description = "upward mixing length"
          zt%f%var(ilup)%units = "m"
          k = k + 1

        case ('ldown')
          ildown = k
          zt%f%var(ildown)%ptr => zt%x(:,k)
          zt%f%var(ildown)%name = "ldown"
          zt%f%var(ildown)%description = "downward mixing length"
          zt%f%var(ildown)%units = "m"
          k = k + 1

        case ('taut')
          itaut = k
          zt%f%var(itaut)%ptr => zt%x(:,k)
          zt%f%var(itaut)%name = "taut"
          zt%f%var(itaut)%description = "dissipation time"
          zt%f%var(itaut)%units = "s"
          k = k + 1

        case ('kht')
          ikht = k
          zt%f%var(ikht)%ptr => zt%x(:,k)
          zt%f%var(ikht)%name = "Kht"
          zt%f%var(ikht)%description = "eddy diffusivity"
          zt%f%var(ikht)%units = "m2/s "
          k = k + 1

        case ('wp2thvp')
          iwp2thvp = k
          zt%f%var(iwp2thvp)%ptr => zt%x(:,k)
          zt%f%var(iwp2thvp)%name = "wp2thvp"
          zt%f%var(iwp2thvp)%description = "wp2thvp"
          zt%f%var(iwp2thvp)%units = "(m2 K)/s2"
          k = k + 1

        case ('wp2rcp')
          iwp2rcp = k
          zt%f%var(iwp2rcp)%ptr => zt%x(:,k)
          zt%f%var(iwp2rcp)%name = "wp2rcp"
          zt%f%var(iwp2rcp)%description = "wp2rcp"
          zt%f%var(iwp2rcp)%units = "(m2 kg)/(s2 kg)"
          k = k + 1

        case ('wprtpthlp')
          iwprtpthlp = k
          zt%f%var(iwprtpthlp)%ptr => zt%x(:,k)
          zt%f%var(iwprtpthlp)%name = "wprtpthlp"
          zt%f%var(iwprtpthlp)%description = "wprtpthlp"
          zt%f%var(iwprtpthlp)%units = "(m kg K)/(s kg)"
          k = k + 1

        case ('sc')
          isct = k
          zt%f%var(isct)%ptr => zt%x(:,k)
          zt%f%var(isct)%name = "Sc"
          zt%f%var(isct)%description = "Sc"
          zt%f%var(isct)%units = "count"
          k = k + 1

        case ('rhot')
          irhot = k
          zt%f%var(irhot)%ptr => zt%x(:,k)
          zt%f%var(irhot)%name = "rhot"
          zt%f%var(irhot)%description = "density"
          zt%f%var(irhot)%units = "kg / m^3"
          k = k + 1

        case ('Ncm')           ! Brian
          iNcm = k
          zt%f%var(iNcm)%ptr => zt%x(:,k)
          zt%f%var(iNcm)%name = "Ncm"
          zt%f%var(iNcm)%description
     .    = "cloud droplet number concentration (num/m^3)"
          zt%f%var(iNcm)%units = "count/m^3"
          k = k + 1

        case ('sed_rcm')       ! Brian
          ised_rcm = k
          zt%f%var(ised_rcm)%ptr => zt%x(:,k)
          zt%f%var(ised_rcm)%name = "sed_rcm"
          zt%f%var(ised_rcm)%description
     .    = "d(rcm)/dt due to cloud sedimentation (kg / [m^2 s])"
          zt%f%var(ised_rcm)%units = "kg/(m^2 s)"
          k = k + 1

        case ('rsm')           ! Brian
          irsm = k
          zt%f%var(irsm)%ptr => zt%x(:,k)
          zt%f%var(irsm)%name = "rsm"
          zt%f%var(irsm)%description
     .    = "saturation mixing ratio (kg/kg)"
          zt%f%var(irsm)%units = "kg/kg"
          k = k + 1

        case ('rrm')           ! Brian
          irrm = k
          zt%f%var(irrm)%ptr => zt%x(:,k)
          zt%f%var(irrm)%name = "rrm"
          zt%f%var(irrm)%description
     .    = "rain water mixing ratio (kg/kg)"
          zt%f%var(irrm)%units = "kg/kg"
          k = k + 1

        case ('Nrm')           ! Brian
          iNrm = k
          zt%f%var(iNrm)%ptr => zt%x(:,k)
          zt%f%var(iNrm)%name = "Nrm"
          zt%f%var(iNrm)%description
     .    = "rain droplet number concentration (num/m^3)"
          zt%f%var(iNrm)%units = "num/(m^3)"
          k = k + 1

        case ('Vrr')           ! Brian
          iVrr = k
          zt%f%var(iVrr)%ptr => zt%x(:,k)
          zt%f%var(iVrr)%name = "Vrr"
          zt%f%var(iVrr)%description
     .    = "rrm sedimentation velocity (m/s)"
          zt%f%var(iVrr)%units = "m/s"
          k = k + 1

        case ('mean_vol_rad')  ! Brian
          imean_vol_rad = k
          zt%f%var(imean_vol_rad)%ptr => zt%x(:,k)
          zt%f%var(imean_vol_rad)%name = "mvr"
          zt%f%var(imean_vol_rad)%description
     .    = "rain drop mean volume radius (m)"
          zt%f%var(imean_vol_rad)%units = "m"
          k = k + 1

        case ('rain_rate')     ! Brian
          irain_rate = k
          zt%f%var(irain_rate)%ptr => zt%x(:,k)
          zt%f%var(irain_rate)%name = "rain_rate"
          zt%f%var(irain_rate)%description = "rain rate (mm/day)"
          zt%f%var(irain_rate)%units = "mm/day"
          k = k + 1
 
        case ('AKm')           ! Vince Larson 22 May 2005
          iAKm = k
          zt%f%var(iAKm)%ptr => zt%x(:,k)
          zt%f%var(iAKm)%name = "AKm"
          zt%f%var(iAKm)%description = "Analytic Kessler ac [kg/kg]"
          zt%f%var(iAKm)%units = "kg/kg"
          k = k + 1
 
        case ('AKm_est')       ! Vince Larson 22 May 2005
          iAKm_est = k
          zt%f%var(iAKm_est)%ptr => zt%x(:,k)
          zt%f%var(iAKm_est)%name = "AKm_est"
          zt%f%var(iAKm_est)%description = "LH Kessler estimate [kg/kg]"
          zt%f%var(iAKm_est)%units = "kg/kg"
          k = k + 1

        case ('radht')
          iradht = k
          zt%f%var(iradht)%ptr => zt%x(:,k)
          zt%f%var(iradht)%name = "radht"
          zt%f%var(iradht)%description = "Heating rate"
          zt%f%var(iradht)%units = "K/s"
          k = k + 1

        case ('radht_LW')
          iradht_LW = k
          zt%f%var(iradht_LW)%ptr => zt%x(:,k)
          zt%f%var(iradht_LW)%name = "radht_LW"
          zt%f%var(iradht_LW)%description = "Long-wave heating rate"
          zt%f%var(iradht_LW)%units = "K/s"
          k = k + 1

        case ('radht_SW')
          iradht_SW = k
          zt%f%var(iradht_SW)%ptr => zt%x(:,k)
          zt%f%var(iradht_SW)%name = "radht_SW"
          zt%f%var(iradht_SW)%description = "Short-wave heating rate"
          zt%f%var(iradht_SW)%units = "K/s"
          k = k + 1
 
        case ('wp3_bt')
          iwp3_bt = k
          zt%f%var(iwp3_bt)%ptr => zt%x(:,k)
          zt%f%var(iwp3_bt)%name = "wp3_bt"
          zt%f%var(iwp3_bt)%description = "wp3 budget"
          zt%f%var(iwp3_bt)%units = "(m^3)/(s^4)"
          k = k + 1
 
        case ('wp3_ma')
          iwp3_ma = k
          zt%f%var(iwp3_ma)%ptr => zt%x(:,k)
          zt%f%var(iwp3_ma)%name = "wp3_ma"
          zt%f%var(iwp3_ma)%description = "wp3 mean advection"
          zt%f%var(iwp3_ma)%units = "(m^3)/(s^4)"
          k = k + 1
 
        case ('wp3_ta')
          iwp3_ta = k
          zt%f%var(iwp3_ta)%ptr => zt%x(:,k)
          zt%f%var(iwp3_ta)%name = "wp3_ta"
          zt%f%var(iwp3_ta)%description = "wp3 turbulent advection"
          zt%f%var(iwp3_ta)%units = "(m^3)/(s^4)"
          k = k + 1
 
        case ('wp3_tp')
          iwp3_tp = k
          zt%f%var(iwp3_tp)%ptr => zt%x(:,k)
          zt%f%var(iwp3_tp)%name = "wp3_tp"
          zt%f%var(iwp3_tp)%description = "wp3 turbulent transport"
          zt%f%var(iwp3_tp)%units = "(m^3)/(s^4)"
          k = k + 1
 
        case ('wp3_ac')
          iwp3_ac = k
          zt%f%var(iwp3_ac)%ptr => zt%x(:,k)
          zt%f%var(iwp3_ac)%name = "wp3_ac"
          zt%f%var(iwp3_ac)%description = "wp3 accumulation term"
          zt%f%var(iwp3_ac)%units = "(m^3)/(s^4)"
          k = k + 1
 
        case ('wp3_bp')
          iwp3_bp = k
          zt%f%var(iwp3_bp)%ptr => zt%x(:,k)
          zt%f%var(iwp3_bp)%name = "wp3_bp"
          zt%f%var(iwp3_bp)%description = "wp3 buoyancy production"
          zt%f%var(iwp3_bp)%units = "(m^3)/(s^4)"
          k = k + 1
 
        case ('wp3_pr1')
          iwp3_pr1 = k
          zt%f%var(iwp3_pr1)%ptr => zt%x(:,k)
          zt%f%var(iwp3_pr1)%name = "wp3_pr1"
          zt%f%var(iwp3_pr1)%description = "wp3 pressure term 1"
          zt%f%var(iwp3_pr1)%units = "(m^3)/(s^4)"
          k = k + 1
 
        case ('wp3_pr2')
          iwp3_pr2 = k
          zt%f%var(iwp3_pr2)%ptr => zt%x(:,k)
          zt%f%var(iwp3_pr2)%name = "wp3_pr2"
          zt%f%var(iwp3_pr2)%description = "wp3 pressure term 2"
          zt%f%var(iwp3_pr2)%units = "(m^3)/(s^4)"
          k = k + 1
 
        case ('wp3_dp1')
          iwp3_dp1 = k
          zt%f%var(iwp3_dp1)%ptr => zt%x(:,k)
          zt%f%var(iwp3_dp1)%name = "wp3_dp1"
          zt%f%var(iwp3_dp1)%description = "wp3 dissipation term 1"
          zt%f%var(iwp3_dp1)%units = "(m^3)/(s^4)"
          k = k + 1
 
        case ('wp3_cl')
          iwp3_cl = k
          zt%f%var(iwp3_cl)%ptr => zt%x(:,k)
          zt%f%var(iwp3_cl)%name = "wp3_cl"
          zt%f%var(iwp3_cl)%description = "wp3 clipping term"
          zt%f%var(iwp3_cl)%units = "(m^3)/(s^4)"
          k = k + 1
 
        case ('rrm_bt')
          irrm_bt = k
          zt%f%var(irrm_bt)%ptr => zt%x(:,k)
          zt%f%var(irrm_bt)%name = "rrm_bt"
          zt%f%var(irrm_bt)%description = "rrm budget"
          zt%f%var(irrm_bt)%units = "(kg/kg)/(s)"
          k = k + 1
 
        case ('rrm_ma')
          irrm_ma = k
          zt%f%var(irrm_ma)%ptr => zt%x(:,k)
          zt%f%var(irrm_ma)%name = "rrm_ma"
          zt%f%var(irrm_ma)%description = "rrm mean advection"
          zt%f%var(irrm_ma)%units = "(kg/kg)/(s)"
          k = k + 1
 
        case ('rrm_sd')
          irrm_sd = k
          zt%f%var(irrm_sd)%ptr => zt%x(:,k)
          zt%f%var(irrm_sd)%name = "rrm_sd"
          zt%f%var(irrm_sd)%description = "rrm sedimentation"
          zt%f%var(irrm_sd)%units = "(kg/kg)/(s)"
          k = k + 1
 
        case ('rrm_dff')
          irrm_dff = k
          zt%f%var(irrm_dff)%ptr => zt%x(:,k)
          zt%f%var(irrm_dff)%name = "rrm_dff"
          zt%f%var(irrm_dff)%description = "rrm diffusion"
          zt%f%var(irrm_dff)%units = "(kg/kg)/(s)"
          k = k + 1
 
        case ('rrm_cond')
          irrm_cond = k
          zt%f%var(irrm_cond)%ptr => zt%x(:,k)
          zt%f%var(irrm_cond)%name = "rrm_cond"
          zt%f%var(irrm_cond)%description = "rrm cond/evap"
          zt%f%var(irrm_cond)%units = "(kg/kg)/(s)"
          k = k + 1
 
        case ('rrm_auto')
          irrm_auto = k
          zt%f%var(irrm_auto)%ptr => zt%x(:,k)
          zt%f%var(irrm_auto)%name = "rrm_auto"
          zt%f%var(irrm_auto)%description = "rrm autoconversion"
          zt%f%var(irrm_auto)%units = "(kg/kg)/(s)"
          k = k + 1
 
        case ('rrm_accr')
          irrm_accr = k
          zt%f%var(irrm_accr)%ptr => zt%x(:,k)
          zt%f%var(irrm_accr)%name = "rrm_accr"
          zt%f%var(irrm_accr)%description = "rrm accretion"
          zt%f%var(irrm_accr)%units = "(kg/kg)/(s)"
          k = k + 1
 
        case ('Nrm_bt')
          iNrm_bt = k
          zt%f%var(iNrm_bt)%ptr => zt%x(:,k)
          zt%f%var(iNrm_bt)%name = "Nrm_bt"
          zt%f%var(iNrm_bt)%description = "Nrm budget"
          zt%f%var(iNrm_bt)%units = "(num/m^3)/(s)"
          k = k + 1
 
        case ('Nrm_ma')
          iNrm_ma = k
          zt%f%var(iNrm_ma)%ptr => zt%x(:,k)
          zt%f%var(iNrm_ma)%name = "Nrm_ma"
          zt%f%var(iNrm_ma)%description = "Nrm mean advection"
          zt%f%var(iNrm_ma)%units = "(num/m^3)/(s)"
          k = k + 1
 
        case ('Nrm_sd')
          iNrm_sd = k
          zt%f%var(iNrm_sd)%ptr => zt%x(:,k)
          zt%f%var(iNrm_sd)%name = "Nrm_sd"
          zt%f%var(iNrm_sd)%description = "Nrm sedimentation"
          zt%f%var(iNrm_sd)%units = "(num/m^3)/(s)"
          k = k + 1
 
        case ('Nrm_dff')
          iNrm_dff = k
          zt%f%var(iNrm_dff)%ptr => zt%x(:,k)
          zt%f%var(iNrm_dff)%name = "Nrm_dff"
          zt%f%var(iNrm_dff)%description = "Nrm diffusion"
          zt%f%var(iNrm_dff)%units = "(num/m^3)/(s)"
          k = k + 1
 
        case ('Nrm_cond')
          iNrm_cond = k
          zt%f%var(iNrm_cond)%ptr => zt%x(:,k)
          zt%f%var(iNrm_cond)%name = "Nrm_cond"
          zt%f%var(iNrm_cond)%description = "Nrm cond/evap"
          zt%f%var(iNrm_cond)%units = "(num/m^3)/(s)"
          k = k + 1
 
        case ('Nrm_auto')
          iNrm_auto = k
          zt%f%var(iNrm_auto)%ptr => zt%x(:,k)
          zt%f%var(iNrm_auto)%name = "Nrm_auto"
          zt%f%var(iNrm_auto)%description = "Nrm autoconversion"
          zt%f%var(iNrm_auto)%units = "(num/m^3)/(s)"
          k = k + 1
 
        case ('a')
          ia = k
          zt%f%var(ia)%ptr => zt%x(:,k)
          zt%f%var(ia)%name = "a"
          zt%f%var(ia)%description = "pdf parameter a"
          zt%f%var(ia)%units = "count"
          k = k + 1
 
        case ('w1')
          iw1 = k
          zt%f%var(iw1)%ptr => zt%x(:,k)
          zt%f%var(iw1)%name = "w1"
          zt%f%var(iw1)%description = "pdf parameter w1"
          zt%f%var(iw1)%units = "m/s"
          k = k + 1
 
        case ('w2')
          iw2 = k
          zt%f%var(iw2)%ptr => zt%x(:,k)
          zt%f%var(iw2)%name = "w2"
          zt%f%var(iw2)%description = "pdf parameter w2"
          zt%f%var(iw2)%units = "m/s"
          k = k + 1
 
        case ('sw1')
          isw1 = k
          zt%f%var(isw1)%ptr => zt%x(:,k)
          zt%f%var(isw1)%name = "sw1"
          zt%f%var(isw1)%description = "pdf parameter sw1"
          zt%f%var(isw1)%units = "m2/s2"
          k = k + 1
 
        case ('sw2')
          isw2 = k
          zt%f%var(isw2)%ptr => zt%x(:,k)
          zt%f%var(isw2)%name = "sw2"
          zt%f%var(isw2)%description = "pdf parameter sw2"
          zt%f%var(isw2)%units = "m2/s2"
          k = k + 1
 
        case ('thl1')
          ithl1 = k
          zt%f%var(ithl1)%ptr => zt%x(:,k)
          zt%f%var(ithl1)%name = "thl1"
          zt%f%var(ithl1)%description = "pdf parameter thl1"
          zt%f%var(ithl1)%units = "K"
          k = k + 1
 
        case ('thl2')
          ithl2 = k
          zt%f%var(ithl2)%ptr => zt%x(:,k)
          zt%f%var(ithl2)%name = "thl2"
          zt%f%var(ithl2)%description = "pdf parameter thl2"
          zt%f%var(ithl2)%units = "K"
          k = k + 1
 
        case ('sthl1')
          isthl1 = k
          zt%f%var(isthl1)%ptr => zt%x(:,k)
          zt%f%var(isthl1)%name = "sthl1"
          zt%f%var(isthl1)%description = "pdf parameter sthl1"
          zt%f%var(isthl1)%units = "K^2"
          k = k + 1
 
        case ('sthl2')
          isthl2 = k
          zt%f%var(isthl2)%ptr => zt%x(:,k)
          zt%f%var(isthl2)%name = "sthl2"
          zt%f%var(isthl2)%description = "pdf parameter sthl2"
          zt%f%var(isthl2)%units = "K^2"
          k = k + 1
 
        case ('rt1')
          irt1 = k
          zt%f%var(irt1)%ptr => zt%x(:,k)
          zt%f%var(irt1)%name = "rt1"
          zt%f%var(irt1)%description = "pdf parameter rt1"
          zt%f%var(irt1)%units = "kg/kg"
          k = k + 1
 
        case ('rt2')
          irt2 = k
          zt%f%var(irt2)%ptr => zt%x(:,k)
          zt%f%var(irt2)%name = "rt2"
          zt%f%var(irt2)%description = "pdf parameter rt2"
          zt%f%var(irt2)%units = "kg/kg"
          k = k + 1
 
        case ('srt1')
          isrt1 = k
          zt%f%var(isrt1)%ptr => zt%x(:,k)
          zt%f%var(isrt1)%name = "srt1"
          zt%f%var(isrt1)%description = "pdf parameter srt1"
          zt%f%var(isrt1)%units = "(kg^2)/(kg^2)"
          k = k + 1
 
        case ('srt2')
          isrt2 = k
          zt%f%var(isrt2)%ptr => zt%x(:,k)
          zt%f%var(isrt2)%name = "srt2"
          zt%f%var(isrt2)%description = "pdf parameter srt2"
          zt%f%var(isrt2)%units = "(kg^2)/(kg^2)"
          k = k + 1
 
        case ('rc1')
          irc1 = k
          zt%f%var(irc1)%ptr => zt%x(:,k)
          zt%f%var(irc1)%name = "rc1"
          zt%f%var(irc1)%description = "pdf parameter rc1"
          zt%f%var(irc1)%units = "kg/kg"
          k = k + 1
 
        case ('rc2')
          irc2 = k
          zt%f%var(irc2)%ptr => zt%x(:,k)
          zt%f%var(irc2)%name = "rc2"
          zt%f%var(irc2)%description = "pdf parameter rc2"
          zt%f%var(irc2)%units = "kg/kg"
          k = k + 1
 
        case ('rsl1')
          irsl1 = k
          zt%f%var(irsl1)%ptr => zt%x(:,k)
          zt%f%var(irsl1)%name = "rsl1"
          zt%f%var(irsl1)%description = "pdf parameter rsl1"
          zt%f%var(irsl1)%units = "kg/kg"
          k = k + 1
 
        case ('rsl2')
          irsl2 = k
          zt%f%var(irsl2)%ptr => zt%x(:,k)
          zt%f%var(irsl2)%name = "rsl2"
          zt%f%var(irsl2)%description = "pdf parameter rsl2"
          zt%f%var(irsl2)%units = "kg/kg"
          k = k + 1
 
        case ('R1')
          iR1 = k
          zt%f%var(iR1)%ptr => zt%x(:,k)
          zt%f%var(iR1)%name = "R1"
          zt%f%var(iR1)%description = "pdf parameter R1"
          zt%f%var(iR1)%units = "count"
          k = k + 1
 
        case ('R2')
          iR2 = k
          zt%f%var(iR2)%ptr => zt%x(:,k)
          zt%f%var(iR2)%name = "R2"
          zt%f%var(iR2)%description = "pdf parameter R2"
          zt%f%var(iR2)%units = "count"
          k = k + 1
 
        case ('s1')
          is1 = k
          zt%f%var(is1)%ptr => zt%x(:,k)
          zt%f%var(is1)%name = "s1"
          zt%f%var(is1)%description = "pdf parameter s1"
          zt%f%var(is1)%units = "kg/kg"
          k = k + 1
 
        case ('s2')
          is2 = k
          zt%f%var(is2)%ptr => zt%x(:,k)
          zt%f%var(is2)%name = "s2"
          zt%f%var(is2)%description = "pdf parameter s2"
          zt%f%var(is2)%units = "kg/kg"
          k = k + 1
 
        case ('ss1')
          iss1 = k
          zt%f%var(iss1)%ptr => zt%x(:,k)
          zt%f%var(iss1)%name = "ss1"
          zt%f%var(iss1)%description = "pdf parameter ss1"
          zt%f%var(iss1)%units = "kg/kg"
          k = k + 1
 
        case ('ss2')
          iss2 = k
          zt%f%var(iss2)%ptr => zt%x(:,k)
          zt%f%var(iss2)%name = "ss2"
          zt%f%var(iss2)%description = "pdf parameter ss2"
          zt%f%var(iss2)%units = "kg/kg"
          k = k + 1
 
        case ('rrtthl')
          irrtthl = k
          zt%f%var(irrtthl)%ptr => zt%x(:,k)
          zt%f%var(irrtthl)%name = "rrtthl"
          zt%f%var(irrtthl)%description = "pdf parameter rrtthl"
          zt%f%var(irrtthl)%units = "count"
          k = k + 1

        case('wp2zt')
          iwp2zt = k
          zt%f%var(iwp2zt)%ptr => zt%x(:,k)
          zt%f%var(iwp2zt)%name = "wp2zt"
          zt%f%var(iwp2zt)%description ="wp2zt"
          zt%f%var(iwp2zt)%units = "m2/s2"
          k = k + 1

        case('thlp2zt')
          ithlp2zt = k
          zt%f%var(ithlp2zt)%ptr => zt%x(:,k)
          zt%f%var(ithlp2zt)%name = "thlp2zt"
          zt%f%var(ithlp2zt)%description ="thlp2zt"
          zt%f%var(ithlp2zt)%units = "K^2"
          k = k + 1

        case('wpthlpzt')
          iwpthlpzt = k
          zt%f%var(iwpthlpzt)%ptr => zt%x(:,k)
          zt%f%var(iwpthlpzt)%name = "wpthlpzt"
          zt%f%var(iwpthlpzt)%description ="wpthlpzt"
          zt%f%var(iwpthlpzt)%units = "(m K)/s"
          k = k + 1

        case('wprtpzt')   
          iwprtpzt = k
          zt%f%var(iwprtpzt)%ptr => zt%x(:,k)
          zt%f%var(iwprtpzt)%name = "wprtpzt"
          zt%f%var(iwprtpzt)%description = "wprtpzt"
          zt%f%var(iwprtpzt)%units = "(m kg)/(s kg)"
          k = k + 1

        case('rtp2zt')   
          irtp2zt = k
          zt%f%var(irtp2zt)%ptr => zt%x(:,k)
          zt%f%var(irtp2zt)%name = "rtp2zt"
          zt%f%var(irtp2zt)%description = "rtp2zt"
          zt%f%var(irtp2zt)%units = "kg/kg"
          k = k + 1

        case('rtpthlpzt')   
          irtpthlpzt = k
          zt%f%var(irtpthlpzt)%ptr => zt%x(:,k)
          zt%f%var(irtpthlpzt)%name = "rtpthlpzt"
          zt%f%var(irtpthlpzt)%description = "rtpthlpzt"
          zt%f%var(irtpthlpzt)%units = "(kg K)/kg"
          k = k + 1
#ifdef SCALARS 
        case ('sclram')
          isclram = k
          zt%f%var(isclram)%ptr => zt%x(:,k)
          zt%f%var(isclram)%name = "sclram"
          zt%f%var(isclram)%description 
     .    = "passive scalar (currently thlm)"
          zt%f%var(isclram)%units = "K"
          k = k + 1

        case ('sclram_f')
          isclram_f = k
          zt%f%var(isclram_f)%ptr => zt%x(:,k)
          zt%f%var(isclram_f)%name = "sclram_f"
          zt%f%var(isclram_f)%description
     .    = "passive scalar forcing (currently thlm_forcing)"
          zt%f%var(isclram_f)%units = "K/s"
          k = k + 1

        case ('sclrbm')
          isclrbm = k
          zt%f%var(isclrbm)%ptr => zt%x(:,k)
          zt%f%var(isclrbm)%name = "sclrbm"
          zt%f%var(isclrbm)%description 
     .    = "passive scalar (currently rtm)"
          zt%f%var(isclrbm)%units = "kg/kg"
          k = k + 1

        case ('sclrbm_f')
          isclrbm_f = k
          zt%f%var(isclrbm_f)%ptr => zt%x(:,k)
          zt%f%var(isclrbm_f)%name = "sclrbm_f"
          zt%f%var(isclrbm_f)%description 
     .    = "passive scalar forcing (currently rtm_forcing)"
          zt%f%var(isclrbm)%units = "kg/(kg s)"
          k = k + 1

        case ('edsclram')
          iedsclram = k
          zt%f%var(iedsclram)%ptr => zt%x(:,k)
          zt%f%var(iedsclram)%name = "edsclram"
          zt%f%var(iedsclram)%description 
     .    = "eddy-diff scalar (currently thlm)"
          zt%f%var(iedsclram)%units = "K"
          k = k + 1

        case ('edsclrbm')
          iedsclrbm = k
          zt%f%var(iedsclrbm)%ptr => zt%x(:,k)
          zt%f%var(iedsclrbm)%name = "edsclrbm"
          zt%f%var(iedsclrbm)%description 
     .    = "eddy-diff scalar (currently rt)"
          zt%f%var(iedsclrbm)%units = "kg/kg"
          k = k + 1
#endif /*SCALARS*/

        case default
          write(*,*) 'Error: unrecognized variable ',trim(vars_zt(i))
          lerror = .true.

        end select

      end do

!     Initialize zm (momentum points)

      i = 1
      do while ( ichar(vars_zm(i)(1:1)) /= 0 
     .           .and. len_trim(vars_zm(i)) /= 0
     .           .and. i <= nvarmax )
        i = i + 1
      end do
      ntot = i - 1
      if ( ntot == nvarmax )
     .  write(*,*) 'WARNING: check nvarmax in statistics.f'

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

!     Allocate scratch space

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

!     Open GrADS files

      fdir = "./"
      fname = trim(fname_zm)
      if ( lgrads ) then
        call open_grads( zm%f, 20, fdir, fname, 
     .                   1, zm%kk, zm%z,
     .                   day, month, year, rlat, rlon,
     .                   time+stats_tout, stats_tout,
     .                   zm%nn )
      else ! lnetcdf
        call open_nc( zm%f, 20, fdir, fname, 
     .                1, zm%kk, zm%z,
     .                day, month, year, rlat, rlon,
     .                time+stats_tout, stats_tout,
     .                zm%nn )
      end if

!     Default initialization for array indices for zm

      iwp2      = 0
      irtp2     = 0
      ithlp2    = 0
      irtpthlp  = 0
      iwprtp    = 0
      iwpthlp   = 0
      iwp4      = 0
      iwpthvp   = 0
      irtpthvp  = 0
      ithlpthvp = 0
      itaum     = 0
      ikhm      = 0
      iwprcp    = 0
      ithlprcp  = 0
      irtprcp   = 0
      iupwp     = 0
      ivpwp     = 0
      irhom     = 0
      iscm      = 0
      iem       = 0
      ishear    = 0    ! Brian
      iFrad     = 0
      iFrad_LW  = 0    ! Brian
      iFrad_SW  = 0    ! Brian
      iFprec    = 0    ! Brian
      iFcsed    = 0    ! Brian

      iwp2_bt   = 0
      iwp2_ma   = 0
      iwp2_ta   = 0
      iwp2_ac   = 0
      iwp2_bp   = 0
      iwp2_pr1  = 0
      iwp2_pr2  = 0
      iwp2_pr3  = 0
      iwp2_dp1  = 0
      iwp2_dp2  = 0
      iwp2_cl   = 0

      iwprtp_bt  = 0
      iwprtp_ma  = 0
      iwprtp_ta  = 0
      iwprtp_tp  = 0
      iwprtp_ac  = 0
      iwprtp_bp  = 0
      iwprtp_pr1 = 0
      iwprtp_pr2 = 0
      iwprtp_pr3 = 0
      iwprtp_dp1 = 0

      iwpthlp_bt  = 0
      iwpthlp_ma  = 0
      iwpthlp_ta  = 0
      iwpthlp_tp  = 0
      iwpthlp_ac  = 0
      iwpthlp_bp  = 0
      iwpthlp_pr1 = 0
      iwpthlp_pr2 = 0
      iwpthlp_pr3 = 0
      iwpthlp_dp1 = 0

!     Dr. Golaz's new budget terms

      irtp2_bt    = 0
      irtp2_ma    = 0
      irtp2_ta    = 0
      irtp2_tp    = 0
      irtp2_dp1   = 0
      irtp2_dp2   = 0
      irtp2_cl    = 0

      ithlp2_bt    = 0
      ithlp2_ma    = 0
      ithlp2_ta    = 0
      ithlp2_tp    = 0
      ithlp2_dp1   = 0
      ithlp2_dp2   = 0
      ithlp2_cl    = 0

      irtpthlp_bt  = 0
      irtpthlp_ma  = 0
      irtpthlp_ta  = 0
      irtpthlp_tp1 = 0
      irtpthlp_tp2 = 0
      irtpthlp_dp1 = 0
      irtpthlp_dp2 = 0

#ifdef SCALARS
      isclraprtp    = 0
      isclrbprtp    = 0
      isclrap2      = 0
      isclrbp2      = 0
      isclrapthvp   = 0
      isclrbpthvp   = 0
      isclrapthlp   = 0
      isclrbpthlp   = 0
      isclraprcp    = 0
      isclrbprcp    = 0
      iwpsclrap     = 0
      iwpsclrbp     = 0
      iwp2sclrap    = 0
      iwp2sclrbp    = 0
      iwpsclrap2    = 0
      iwpsclrbp2    = 0
      iwpsclraprtp  = 0
      iwpsclrbprtp  = 0
      iwpsclrapthlp = 0
      iwpsclrbpthlp = 0

      iwpedsclrap   = 0
      iwpedsclrbp   = 0
#endif /*SCALARS*/

!     Assign pointers for statistics variables zm

      k = 1
      do i=1,zm%nn

        select case ( trim(vars_zm(i)) )

        case ('wp2')
          iwp2 = k
          zm%f%var(iwp2)%ptr => zm%x(:,k)
          zm%f%var(iwp2)%name = "wp2"
          zm%f%var(iwp2)%description = "wp2"
          zm%f%var(iwp2)%units = "m2/s2"
          k = k + 1

        case ('rtp2')
          irtp2 = k
          zm%f%var(irtp2)%ptr => zm%x(:,k)
          zm%f%var(irtp2)%name = "rtp2"
          zm%f%var(irtp2)%description = "rtp2"
          zm%f%var(irtp2)%units = "(kg/kg)^2"
          k = k + 1

        case ('thlp2')
          ithlp2 = k
          zm%f%var(ithlp2)%ptr => zm%x(:,k)
          zm%f%var(ithlp2)%name = "thlp2"
          zm%f%var(ithlp2)%description = "thlp2"
          zm%f%var(ithlp2)%units = "K^2"
          k = k + 1

        case ('rtpthlp')
          irtpthlp = k
          zm%f%var(irtpthlp)%ptr => zm%x(:,k)
          zm%f%var(irtpthlp)%name = "rtpthlp"
          zm%f%var(irtpthlp)%description = "rtpthlp"
          zm%f%var(irtpthlp)%units = "(kg K)/kg"
          k = k + 1

        case ('wprtp')
          iwprtp = k
          zm%f%var(iwprtp)%ptr => zm%x(:,k)
          zm%f%var(iwprtp)%name = "wprtp"
          zm%f%var(iwprtp)%description = "wprtp"
          zm%f%var(iwprtp)%units = "(m kg)/(s kg)"
          k = k + 1

        case ('wpthlp')
          iwpthlp = k
          zm%f%var(iwpthlp)%ptr => zm%x(:,k)
          zm%f%var(iwpthlp)%name = "wpthlp"
          zm%f%var(iwpthlp)%description = "wpthlp"
          zm%f%var(iwpthlp)%units = "(m K)/s"
          k = k + 1

        case ('wp4')
          iwp4 = k
          zm%f%var(iwp4)%ptr => zm%x(:,k)
          zm%f%var(iwp4)%name = "wp4"
          zm%f%var(iwp4)%description = "wp4"
          zm%f%var(iwp4)%units = "(m^4)/(s^4)"
          k = k + 1

        case ('wpthvp')
          iwpthvp = k
          zm%f%var(iwpthvp)%ptr => zm%x(:,k)
          zm%f%var(iwpthvp)%name = "wpthvp"
          zm%f%var(iwpthvp)%description = "Buoyancy flux (K m/s)"
          zm%f%var(iwpthvp)%units = "(K m)/s"
          k = k + 1

        case ('rtpthvp')
          irtpthvp = k
          zm%f%var(irtpthvp)%ptr => zm%x(:,k)
          zm%f%var(irtpthvp)%name = "rtpthvp"
          zm%f%var(irtpthvp)%description = "rtpthvp"
          zm%f%var(irtpthvp)%units = "(kg K)/kg"
          k = k + 1

        case ('thlpthvp')
          ithlpthvp = k
          zm%f%var(ithlpthvp)%ptr => zm%x(:,k)
          zm%f%var(ithlpthvp)%name = "thlpthvp"
          zm%f%var(ithlpthvp)%description = "thlpthvp"
          zm%f%var(ithlpthvp)%units = "K^2"
          k = k + 1

        case ('taum')
          itaum = k
          zm%f%var(itaum)%ptr => zm%x(:,k)
          zm%f%var(itaum)%name = "taum"
          zm%f%var(itaum)%description = "dissipation time"
          zm%f%var(itaum)%units = "s"
          k = k + 1

        case ('khm')
          ikhm = k
          zm%f%var(ikhm)%ptr => zm%x(:,k)
          zm%f%var(ikhm)%name = "khm"
          zm%f%var(ikhm)%description = "eddy diffusivity"
          zm%f%var(ikhm)%units = "m2/s"
          k = k + 1

        case ('wprcp')
          iwprcp = k
          zm%f%var(iwprcp)%ptr => zm%x(:,k)
          zm%f%var(iwprcp)%name = "wprcp"
          zm%f%var(iwprcp)%description = "wprcp"
          zm%f%var(iwprcp)%units = "(m kg)/(s kg)"
          k = k + 1

        case ('thlprcp')
          ithlprcp = k
          zm%f%var(ithlprcp)%ptr => zm%x(:,k)
          zm%f%var(ithlprcp)%name = "thlprcp"
          zm%f%var(ithlprcp)%description = "thlprcp"
          zm%f%var(ithlprcp)%units = "(K kg)/(kg)"
          k = k + 1

        case ('rtprcp')
          irtprcp = k
          zm%f%var(irtprcp)%ptr => zm%x(:,k)
          zm%f%var(irtprcp)%name = "rtprcp"
          zm%f%var(irtprcp)%description = "rtprcp"
          zm%f%var(irtprcp)%units = "(kg^2)/(kg^2)"
          k = k + 1

        case ('upwp')
          iupwp = k
          zm%f%var(iupwp)%ptr => zm%x(:,k)
          zm%f%var(iupwp)%name = "upwp"
          zm%f%var(iupwp)%description = "upwp"
          zm%f%var(iupwp)%units = "m2/s2"
          k = k + 1

        case ('vpwp')
          ivpwp = k
          zm%f%var(ivpwp)%ptr => zm%x(:,k)
          zm%f%var(ivpwp)%name = "vpwp"
          zm%f%var(ivpwp)%description = "vpwp"
          zm%f%var(ivpwp)%units = "m2/s2"
          k = k + 1

        case ('rhom')
          irhom = k
          zm%f%var(irhom)%ptr => zm%x(:,k)
          zm%f%var(irhom)%name = "rhom"
          zm%f%var(irhom)%description = "density"
          zm%f%var(irhom)%units = "kg/(m^3)"
          k = k + 1

        case ('sc')
          iscm = k
          zm%f%var(iscm)%ptr => zm%x(:,k)
          zm%f%var(iscm)%name = "scm"
          zm%f%var(iscm)%description = "scm"
          zm%f%var(iscm)%units = "count"
          k = k + 1

        case ('em')
          iem = k
          zm%f%var(iem)%ptr => zm%x(:,k)
          zm%f%var(iem)%name = "em"
          zm%f%var(iem)%description = "em"
          zm%f%var(iem)%units = "m2/s2"
          k = k + 1

        case ('shear')      ! Brian
          ishear = k
          zm%f%var(ishear)%ptr => zm%x(:,k)
          zm%f%var(ishear)%name = "shear"
          zm%f%var(ishear)%description = "wind shear term (m^2/s^3)"
          zm%f%var(ishear)%units = "m2/s3"
          k = k + 1

        case ('Frad')
          iFrad = k
          zm%f%var(iFrad)%ptr => zm%x(:,k)
          zm%f%var(iFrad)%name = "Frad"
          zm%f%var(iFrad)%description = "radiative flux"
          zm%f%var(iFrad)%units = "W/(m^-2)"
          k = k + 1

        case ('Frad_LW')    ! Brian
          iFrad_LW = k
          zm%f%var(iFrad_LW)%ptr => zm%x(:,k)
          zm%f%var(iFrad_LW)%name = "Frad_LW"
          zm%f%var(iFrad_LW)%description = "Long-wave radiative flux"
          zm%f%var(iFrad_LW)%units = "W/(m^-2)"
          k = k + 1

        case ('Frad_SW')    ! Brian
          iFrad_SW = k
          zm%f%var(iFrad_SW)%ptr => zm%x(:,k)
          zm%f%var(iFrad_SW)%name = "Frad_SW"
          zm%f%var(iFrad_SW)%description = "Short-wave radiative flux"
          zm%f%var(iFrad_SW)%units = "W/(m^-2)"
          k = k + 1

        case ('Fprec')      ! Brian
          iFprec = k
          zm%f%var(iFprec)%ptr => zm%x(:,k)
          zm%f%var(iFprec)%name = "Fprec"
          zm%f%var(iFprec)%description = "precipitation flux (W/m^2)"
          zm%f%var(iFprec)%units = "W/(m^2)"
          k = k + 1
 
        case ('Fcsed')      ! Brian
          iFcsed = k
          zm%f%var(iFcsed)%ptr => zm%x(:,k)
          zm%f%var(iFcsed)%name = "Fcsed"
          zm%f%var(iFcsed)%description
     .    = "cloud water sedimentation flux (kg/(s*m^2))"
          zm%f%var(iFcsed)%units = "kg/(s m2)"
          k = k + 1
 
        case ('wp2_bt')
          iwp2_bt = k
          zm%f%var(iwp2_bt)%ptr => zm%x(:,k)
          zm%f%var(iwp2_bt)%name = "wp2_bt"
          zm%f%var(iwp2_bt)%description = "wp2 budget"
          zm%f%var(iwp2_bt)%units = "m2/s3"
          k = k + 1

        case ('wp2_ma')
          iwp2_ma = k
          zm%f%var(iwp2_ma)%ptr => zm%x(:,k)
          zm%f%var(iwp2_ma)%name = "wp2_ma"
          zm%f%var(iwp2_ma)%description = "wp2 mean advection"
          zm%f%var(iwp2_ma)%units = "m2/s3"
          k = k + 1

        case ('wp2_ta')
          iwp2_ta = k
          zm%f%var(iwp2_ta)%ptr => zm%x(:,k)
          zm%f%var(iwp2_ta)%name = "wp2_ta"
          zm%f%var(iwp2_ta)%description = "wp2 turbulent advection"
          zm%f%var(iwp2_ta)%units = "m2/s3"
          k = k + 1

        case ('wp2_ac')
          iwp2_ac = k
          zm%f%var(iwp2_ac)%ptr => zm%x(:,k)
          zm%f%var(iwp2_ac)%name = "wp2_ac"
          zm%f%var(iwp2_ac)%description = "wp2 accumulation term"
          zm%f%var(iwp2_ac)%units = "m2/s3"
          k = k + 1

        case ('wp2_bp')
          iwp2_bp = k
          zm%f%var(iwp2_bp)%ptr => zm%x(:,k)
          zm%f%var(iwp2_bp)%name = "wp2_bp"
          zm%f%var(iwp2_bp)%description = "wp2 buoyancy production"
          zm%f%var(iwp2_bp)%units = "m2/s3"
          k = k + 1

        case ('wp2_pr1')
          iwp2_pr1 = k
          zm%f%var(iwp2_pr1)%ptr => zm%x(:,k)
          zm%f%var(iwp2_pr1)%name = "wp2_pr1"
          zm%f%var(iwp2_pr1)%description = "wp2 pressure term 1"
          zm%f%var(iwp2_pr1)%units = "m2/s3"
          k = k + 1

        case ('wp2_pr2')
          iwp2_pr2 = k
          zm%f%var(iwp2_pr2)%ptr => zm%x(:,k)
          zm%f%var(iwp2_pr2)%name = "wp2_pr2"
          zm%f%var(iwp2_pr2)%description = "wp2 pressure term 2"
          zm%f%var(iwp2_pr2)%units = "m2/s3"
          k = k + 1

        case ('wp2_pr3')
          iwp2_pr3 = k
          zm%f%var(iwp2_pr3)%ptr => zm%x(:,k)
          zm%f%var(iwp2_pr3)%name = "wp2_pr3"
          zm%f%var(iwp2_pr3)%description = "wp2 pressure term 3"
          zm%f%var(iwp2_pr3)%units = "m2/s3"
          k = k + 1

        case ('wp2_dp1')
          iwp2_dp1 = k
          zm%f%var(iwp2_dp1)%ptr => zm%x(:,k)
          zm%f%var(iwp2_dp1)%name = "wp2_dp1"
          zm%f%var(iwp2_dp1)%description = "wp2 dissipation term 1"
          zm%f%var(iwp2_dp1)%units = "m2/s3"
          k = k + 1

        case ('wp2_dp2')
          iwp2_dp2 = k
          zm%f%var(iwp2_dp2)%ptr => zm%x(:,k)
          zm%f%var(iwp2_dp2)%name = "wp2_dp2"
          zm%f%var(iwp2_dp2)%description = "wp2 dissipation term 2"
          zm%f%var(iwp2_dp2)%units = "m2/s3"
          k = k + 1

        case ('wp2_cl')
          iwp2_cl = k
          zm%f%var(iwp2_cl)%ptr => zm%x(:,k)
          zm%f%var(iwp2_cl)%name = "wp2_cl"
          zm%f%var(iwp2_cl)%description = "wp2 clipping term"
          zm%f%var(iwp2_cl)%units = "m2/s3"
          k = k + 1

        case ('wprtp_bt')
          iwprtp_bt = k
          zm%f%var(iwprtp_bt)%ptr => zm%x(:,k)
          zm%f%var(iwprtp_bt)%name = "wprtp_bt"
          zm%f%var(iwprtp_bt)%description = "wprtp budget"
          zm%f%var(iwprtp_bt)%units = "(m kg)/(s2 kg)"
          k = k + 1

        case ('wprtp_ma')
          iwprtp_ma = k
          zm%f%var(iwprtp_ma)%ptr => zm%x(:,k)
          zm%f%var(iwprtp_ma)%name = "wprtp_ma"
          zm%f%var(iwprtp_ma)%description = "wprtp mean advection"
          zm%f%var(iwprtp_ma)%units = "(m kg)/(s2 kg)"
          k = k + 1

        case ('wprtp_ta')
          iwprtp_ta = k
          zm%f%var(iwprtp_ta)%ptr => zm%x(:,k)
          zm%f%var(iwprtp_ta)%name = "wprtp_ta"
          zm%f%var(iwprtp_ta)%description = "wprtp turbulent advection"
          zm%f%var(iwprtp_ta)%units = "(m kg)/(s2 kg)"
          k = k + 1

        case ('wprtp_tp')
          iwprtp_tp = k
          zm%f%var(iwprtp_tp)%ptr => zm%x(:,k)
          zm%f%var(iwprtp_tp)%name = "wprtp_tp"
          zm%f%var(iwprtp_tp)%description = "wprtp turbulent production"
          zm%f%var(iwprtp_tp)%units = "(m kg)/(s2 kg)"
          k = k + 1

        case ('wprtp_ac')
          iwprtp_ac = k
          zm%f%var(iwprtp_ac)%ptr => zm%x(:,k)
          zm%f%var(iwprtp_ac)%name = "wprtp_ac"
          zm%f%var(iwprtp_ac)%description = "wprtp accumulation term"
          zm%f%var(iwprtp_ac)%units = "(m kg)/(s2 kg)"
          k = k + 1

        case ('wprtp_bp')
          iwprtp_bp = k
          zm%f%var(iwprtp_bp)%ptr => zm%x(:,k)
          zm%f%var(iwprtp_bp)%name = "wprtp_bp"
          zm%f%var(iwprtp_bp)%description = "wprtp buoyancy production"
          zm%f%var(iwprtp_bp)%units = "(m kg)/(s2 kg)"
          k = k + 1

        case ('wprtp_pr1')
          iwprtp_pr1 = k
          zm%f%var(iwprtp_pr1)%ptr => zm%x(:,k)
          zm%f%var(iwprtp_pr1)%name = "wprtp_pr1"
          zm%f%var(iwprtp_pr1)%description = "wprtp pressure term 1"
          zm%f%var(iwprtp_pr1)%units = "(m kg)/(s2 kg)"
          k = k + 1

        case ('wprtp_pr2')
          iwprtp_pr2 = k
          zm%f%var(iwprtp_pr2)%ptr => zm%x(:,k)
          zm%f%var(iwprtp_pr2)%name = "wprtp_pr2"
          zm%f%var(iwprtp_pr2)%description = "wprtp pressure term 2"
          zm%f%var(iwprtp_pr2)%units = "(m kg)/(s2 kg)"
          k = k + 1

        case ('wprtp_pr3')
          iwprtp_pr3 = k
          zm%f%var(iwprtp_pr3)%ptr => zm%x(:,k)
          zm%f%var(iwprtp_pr3)%name = "wprtp_pr3"
          zm%f%var(iwprtp_pr3)%description = "wprtp pressure term 3"
          zm%f%var(iwprtp_pr3)%units = "(m kg)/(s2 kg)"
          k = k + 1

        case ('wprtp_dp1')
          iwprtp_dp1 = k
          zm%f%var(iwprtp_dp1)%ptr => zm%x(:,k)
          zm%f%var(iwprtp_dp1)%name = "wprtp_dp1"
          zm%f%var(iwprtp_dp1)%description = "wprtp dissipation term 1"
          zm%f%var(iwprtp_dp1)%units = "(m kg)/(s2 kg)"
          k = k + 1

        case ('wpthlp_bt')
          iwpthlp_bt = k
          zm%f%var(iwpthlp_bt)%ptr => zm%x(:,k)
          zm%f%var(iwpthlp_bt)%name = "wpthlp_bt"
          zm%f%var(iwpthlp_bt)%description = "wpthlp budget"
          zm%f%var(iwpthlp_bt)%units = "(m K)/s2"
          k = k + 1

        case ('wpthlp_ma')
          iwpthlp_ma = k
          zm%f%var(iwpthlp_ma)%ptr => zm%x(:,k)
          zm%f%var(iwpthlp_ma)%name = "wpthlp_ma"
          zm%f%var(iwpthlp_ma)%description = "wpthlp mean advection"
          zm%f%var(iwpthlp_ma)%units = "(m K)/s2"
          k = k + 1

        case ('wpthlp_ta')
          iwpthlp_ta = k
          zm%f%var(iwpthlp_ta)%ptr => zm%x(:,k)
          zm%f%var(iwpthlp_ta)%name = "wpthlp_ta"
          zm%f%var(iwpthlp_ta)%description 
     .    = "wpthlp turbulent advection"
          zm%f%var(iwpthlp_ta)%units = "(m K)/s2"
          k = k + 1

        case ('wpthlp_tp')
          iwpthlp_tp = k
          zm%f%var(iwpthlp_tp)%ptr => zm%x(:,k)
          zm%f%var(iwpthlp_tp)%name = "wpthlp_tp"
          zm%f%var(iwpthlp_tp)%description
     .    = "wpthlp turbulent production"
          zm%f%var(iwpthlp_tp)%units = "(m K)/s2"
          k = k + 1

        case ('wpthlp_ac')
          iwpthlp_ac = k
          zm%f%var(iwpthlp_ac)%ptr => zm%x(:,k)
          zm%f%var(iwpthlp_ac)%name = "wpthlp_ac"
          zm%f%var(iwpthlp_ac)%description
     .    = "wpthlp accumulation term"
          zm%f%var(iwpthlp_ac)%units = "(m K)/s2"
          k = k + 1

        case ('wpthlp_bp')
          iwpthlp_bp = k
          zm%f%var(iwpthlp_bp)%ptr => zm%x(:,k)
          zm%f%var(iwpthlp_bp)%name = "wpthlp_bp"
          zm%f%var(iwpthlp_bp)%description
     .    = "wpthlp buoyancy production"
          zm%f%var(iwpthlp_bp)%units = "(m K)/s2"
          k = k + 1

        case ('wpthlp_pr1')
          iwpthlp_pr1 = k
          zm%f%var(iwpthlp_pr1)%ptr => zm%x(:,k)
          zm%f%var(iwpthlp_pr1)%name = "wpthlp_pr1"
          zm%f%var(iwpthlp_pr1)%description
     .    = "wpthlp pressure term 1"
          zm%f%var(iwpthlp_pr1)%units = "(m K)/s2"
          k = k + 1

        case ('wpthlp_pr2')
          iwpthlp_pr2 = k
          zm%f%var(iwpthlp_pr2)%ptr => zm%x(:,k)
          zm%f%var(iwpthlp_pr2)%name = "wpthlp_pr2"
          zm%f%var(iwpthlp_pr2)%description
     .    = "wpthlp pressure term 2"
          zm%f%var(iwpthlp_pr2)%units = "(m K)/s2"
          k = k + 1

        case ('wpthlp_pr3')
          iwpthlp_pr3 = k
          zm%f%var(iwpthlp_pr3)%ptr => zm%x(:,k)
          zm%f%var(iwpthlp_pr3)%name = "wpthlp_pr3"
          zm%f%var(iwpthlp_pr3)%description
     .    = "wpthlp pressure term 3"
          zm%f%var(iwpthlp_pr3)%units = "(m K)/s2"
          k = k + 1

        case ('wpthlp_dp1')
          iwpthlp_dp1 = k
          zm%f%var(iwpthlp_dp1)%ptr => zm%x(:,k)
          zm%f%var(iwpthlp_dp1)%name = "wpthlp_dp1"
          zm%f%var(iwpthlp_dp1)%description
     .    = "wpthlp dissipation term 1"
          zm%f%var(iwpthlp_dp1)%units = "(m K)/s2"
          k = k + 1

!     Dr. Golaz's new budget terms
        case ('rtp2_bt')
          irtp2_bt = k
          zm%f%var(irtp2_bt)%ptr => zm%x(:,k)
          zm%f%var(irtp2_bt)%name = "rtp2_bt"
          zm%f%var(irtp2_bt)%description = "rtp2 budget"
          zm%f%var(irtp2_bt)%units = "kg/(kg s)"
          k = k + 1

        case ('rtp2_ma')
          irtp2_ma = k
          zm%f%var(irtp2_ma)%ptr => zm%x(:,k)
          zm%f%var(irtp2_ma)%name = "rtp2_ma"
          zm%f%var(irtp2_ma)%description = "rtp2 mean advection"
          zm%f%var(irtp2_ma)%units = "kg/(kg s)"
          k = k + 1

        case ('rtp2_ta')
          irtp2_ta = k
          zm%f%var(irtp2_ta)%ptr => zm%x(:,k)
          zm%f%var(irtp2_ta)%name = "rtp2_ta"
          zm%f%var(irtp2_ta)%description
     .    = "rtp2 turbulent advection"
          zm%f%var(irtp2_ta)%units = "kg/(kg s)"
          k = k + 1

        case ('rtp2_tp')
          irtp2_tp = k
          zm%f%var(irtp2_tp)%ptr => zm%x(:,k)
          zm%f%var(irtp2_tp)%name = "rtp2_tp"
          zm%f%var(irtp2_tp)%description
     .    = "rtp2 turbulent production"
          zm%f%var(irtp2_tp)%units = "kg/(kg s)"
          k = k + 1

        case ('rtp2_dp1')
          irtp2_dp1 = k
          zm%f%var(irtp2_dp1)%ptr => zm%x(:,k)
          zm%f%var(irtp2_dp1)%name = "rtp2_dp1"
          zm%f%var(irtp2_dp1)%description
     .    = "rtp2 dissipation term 1"
          zm%f%var(irtp2_dp1)%units = "kg/(kg s)"
          k = k + 1

        case ('rtp2_dp2')
          irtp2_dp2 = k
          zm%f%var(irtp2_dp2)%ptr => zm%x(:,k)
          zm%f%var(irtp2_dp2)%name = "rtp2_dp2"
          zm%f%var(irtp2_dp2)%description
     .    = "rtp2 dissipation term 2"
          zm%f%var(irtp2_dp2)%units = "kg/(kg s)"
          k = k + 1

        case ('rtp2_cl')
          irtp2_cl = k
          zm%f%var(irtp2_cl)%ptr => zm%x(:,k)
          zm%f%var(irtp2_cl)%name = "rtp2_cl"
          zm%f%var(irtp2_cl)%description
     .    = "rtp2 clipping term"
          zm%f%var(irtp2_cl)%units = "kg/(kg s)"
          k = k + 1

        case ('thlp2_bt')
          ithlp2_bt = k
          zm%f%var(ithlp2_bt)%ptr => zm%x(:,k)
          zm%f%var(ithlp2_bt)%name = "thlp2_bt"
          zm%f%var(ithlp2_bt)%description = "thlp2 budget"
          zm%f%var(ithlp2_bt)%units = "(K^2)/s"
          k = k + 1

        case ('thlp2_ma')
          ithlp2_ma = k
          zm%f%var(ithlp2_ma)%ptr => zm%x(:,k)
          zm%f%var(ithlp2_ma)%name = "thlp2_ma"
          zm%f%var(ithlp2_ma)%description = "thlp2 mean advection"
          zm%f%var(ithlp2_ma)%units = "(K^2)/s"
          k = k + 1

        case ('thlp2_ta')
          ithlp2_ta = k
          zm%f%var(ithlp2_ta)%ptr => zm%x(:,k)
          zm%f%var(ithlp2_ta)%name = "thlp2_ta"
          zm%f%var(ithlp2_ta)%description
     .    = "thlp2 turbulent advection"
          zm%f%var(ithlp2_ta)%units = "(K^2)/s"
          k = k + 1

        case ('thlp2_tp')
          ithlp2_tp = k
          zm%f%var(ithlp2_tp)%ptr => zm%x(:,k)
          zm%f%var(ithlp2_tp)%name = "thlp2_tp"
          zm%f%var(ithlp2_tp)%description
     .    = "thlp2 turbulent production"
          zm%f%var(ithlp2_tp)%units = "(K^2)/s"
          k = k + 1

        case ('thlp2_dp1')
          ithlp2_dp1 = k
          zm%f%var(ithlp2_dp1)%ptr => zm%x(:,k)
          zm%f%var(ithlp2_dp1)%name = "thlp2_dp1"
          zm%f%var(ithlp2_dp1)%description
     .    = "thlp2 dissipation term 1"
          zm%f%var(ithlp2_dp1)%units = "(K^2)/s"
          k = k + 1

        case ('thlp2_dp2')
          ithlp2_dp2 = k
          zm%f%var(ithlp2_dp2)%ptr => zm%x(:,k)
          zm%f%var(ithlp2_dp2)%name = "thlp2_dp2"
          zm%f%var(ithlp2_dp2)%description
     .    = "thlp2 dissipation term 2"
          zm%f%var(ithlp2_dp2)%units = "(K^2)/s"
          k = k + 1

        case ('thlp2_cl')
          ithlp2_cl = k
          zm%f%var(ithlp2_cl)%ptr => zm%x(:,k)
          zm%f%var(ithlp2_cl)%name = "thlp2_cl"
          zm%f%var(ithlp2_cl)%description
     .    = "thlp2 clipping term"
          zm%f%var(ithlp2_cl)%units = "(K^2)/s"
          k = k + 1

        case ('rtpthlp_bt')
          irtpthlp_bt = k
          zm%f%var(irtpthlp_bt)%ptr => zm%x(:,k)
          zm%f%var(irtpthlp_bt)%name = "rtpthlp_bt"
          zm%f%var(irtpthlp_bt)%description = "rtpthlp budget"
          zm%f%var(irtpthlp_bt)%units = "(kg K)/(kg s)"
          k = k + 1

        case ('rtpthlp_ma')
          irtpthlp_ma = k
          zm%f%var(irtpthlp_ma)%ptr => zm%x(:,k)
          zm%f%var(irtpthlp_ma)%name = "rtpthlp_ma"
          zm%f%var(irtpthlp_ma)%description = "rtpthlp mean advection"
          zm%f%var(irtpthlp_ma)%units = "(kg K)/(kg s)"
          k = k + 1

        case ('rtpthlp_ta')
          irtpthlp_ta = k
          zm%f%var(irtpthlp_ta)%ptr => zm%x(:,k)
          zm%f%var(irtpthlp_ta)%name = "rtpthlp_ta"
          zm%f%var(irtpthlp_ta)%description
     .    = "rtpthlp turbulent advection"
          zm%f%var(irtpthlp_ta)%units = "(kg K)/(kg s)"
          k = k + 1

        case ('rtpthlp_tp1')
          irtpthlp_tp1 = k
          zm%f%var(irtpthlp_tp1)%ptr => zm%x(:,k)
          zm%f%var(irtpthlp_tp1)%name = "rtpthlp_tp1"
          zm%f%var(irtpthlp_tp1)%description
     .    = "rtpthlp turbulent production 1"
          zm%f%var(irtpthlp_tp1)%units = "(kg K)/(kg s)"
          k = k + 1

        case ('rtpthlp_tp2')
          irtpthlp_tp2 = k
          zm%f%var(irtpthlp_tp2)%ptr => zm%x(:,k)
          zm%f%var(irtpthlp_tp2)%name = "rtpthlp_tp2"
          zm%f%var(irtpthlp_tp2)%description
     .    = "rtpthlp turbulent production 2"
          zm%f%var(irtpthlp_tp2)%units = "(kg K)/(kg s)"
          k = k + 1

        case ('rtpthlp_dp1')
          irtpthlp_dp1 = k
          zm%f%var(irtpthlp_dp1)%ptr => zm%x(:,k)
          zm%f%var(irtpthlp_dp1)%name = "rtpthlp_dp1"
          zm%f%var(irtpthlp_dp1)%description
     .    = "rtpthlp dissipation term 1"
          zm%f%var(irtpthlp_dp1)%units = "(kg K)/(kg s)"
          k = k + 1

        case ('rtpthlp_dp2')
          irtpthlp_dp2 = k
          zm%f%var(irtpthlp_dp2)%ptr => zm%x(:,k)
          zm%f%var(irtpthlp_dp2)%name = "rtpthlp_dp2"
          zm%f%var(irtpthlp_dp2)%description
     .    = "rtpthlp dissipation term 2"
          zm%f%var(irtpthlp_dp2)%units = "(kg K)/(kg s)"
          k = k + 1

#ifdef SCALARS
        case ('sclraprtp')
          isclraprtp = k
          zm%f%var(isclraprtp)%ptr => zm%x(:,k)
          zm%f%var(isclraprtp)%name = "sclraprtp"
          zm%f%var(isclraprtp)%description 
     .    = "scalar(2)'rt'"
          zm%f%var(isclraprtp)%units = "unknown"
          k = k + 1
        case ('sclrbprtp')
          isclrbprtp = k
          zm%f%var(isclrbprtp)%ptr => zm%x(:,k)
          zm%f%var(isclrbprtp)%name = "sclrbprtp"
          zm%f%var(isclrbprtp)%description
     .    = "scalar(2)'rt'"
          zm%f%var(isclrbprtp)%units = "unknown"
          k = k + 1
        case ('sclrap2')
          isclrap2 = k
          zm%f%var(isclrap2)%ptr => zm%x(:,k)
          zm%f%var(isclrap2)%name = "sclrap2"
          zm%f%var(isclrap2)%description 
     .    = "scalar(1)'^2 "
          zm%f%var(isclrap2)%units = "unknown"
          k = k + 1
        case ('sclrbp2')
          isclrbp2 = k
          zm%f%var(isclrbp2)%ptr => zm%x(:,k)
          zm%f%var(isclrbp2)%name = "sclrbp2"
          zm%f%var(isclrbp2)%description 
     .    = "scalar(2)'^2"
          zm%f%var(isclrbp2)%units = "unknown"
          k = k + 1
        case ('sclrapthvp')
          isclrapthvp = k
          zm%f%var(isclrapthvp)%ptr => zm%x(:,k)
          zm%f%var(isclrapthvp)%name = "sclrapthvp"
          zm%f%var(isclrapthvp)%description
     .    = "scalar(1)'th_v'"
          zm%f%var(isclrapthvp)%units = "unknown"
          k = k + 1
        case ('sclrbpthvp')
          isclrbpthvp = k
          zm%f%var(isclrbpthvp)%ptr => zm%x(:,k)
          zm%f%var(isclrbpthvp)%name = "sclrbpthvp"
          zm%f%var(isclrbpthvp)%description
     .    = "scalar(2)'th_v'"
          zm%f%var(isclrbpthvp)%units = "unknown"
          k = k + 1
        case ('sclrapthlp')
          isclrapthlp = k
          zm%f%var(isclrapthlp)%ptr => zm%x(:,k)
          zm%f%var(isclrapthlp)%name = "sclrapthlp"
          zm%f%var(isclrapthlp)%description
     .    = "scalar(1)'th_l'"
          zm%f%var(isclrapthlp)%units = "unknown"
          k = k + 1
        case ('sclrbpthlp')
          isclrbpthlp = k
          zm%f%var(isclrbpthlp)%ptr => zm%x(:,k)
          zm%f%var(isclrbpthlp)%name = "sclrbpthlp"
          zm%f%var(isclrbpthlp)%description
     .    = "scalar(2)'th_l'"
          zm%f%var(isclrbpthlp)%units = "unknown"
          k = k + 1
        case ('sclraprcp')
          isclraprcp = k
          zm%f%var(isclraprcp)%ptr => zm%x(:,k)
          zm%f%var(isclraprcp)%name = "sclraprcp"
          zm%f%var(isclraprcp)%description 
     .    = "scalar(1)'rc'"
          zm%f%var(isclraprcp)%units = "unknown"
          k = k + 1
        case ('sclrbprcp')
          isclrbprcp = k
          zm%f%var(isclrbprcp)%ptr => zm%x(:,k)
          zm%f%var(isclrbprcp)%name = "sclrbprcp"
          zm%f%var(isclrbprcp)%description
     .    = "scalar(2)'rc'"
          zm%f%var(isclrbprcp)%units = "unknown"
          k = k + 1
        case ('wpsclrap')
          iwpsclrap = k
          zm%f%var(iwpsclrap)%ptr => zm%x(:,k)
          zm%f%var(iwpsclrap)%name = "wpsclrap"
          zm%f%var(iwpsclrap)%description 
     .    = "w'scalar(1)'"
          zm%f%var(iwpsclrap)%units = "unknown"
          k = k + 1
        case ('wpsclrbp')
          iwpsclrbp = k
          zm%f%var(iwpsclrbp)%ptr => zm%x(:,k)
          zm%f%var(iwpsclrbp)%name = "wpsclrbp"
          zm%f%var(iwpsclrbp)%description 
     .    = "w'scalar(2)'"
          zm%f%var(iwpsclrbp)%units = "unknown"
          k = k + 1
        case ('wpsclrap2')
          iwpsclrap2 = k
          zm%f%var(iwpsclrap2)%ptr => zm%x(:,k)
          zm%f%var(iwpsclrap2)%name = "wpsclrap2"
          zm%f%var(iwpsclrap2)%description 
     .    = "w'scalar(1)'^2 "
          zm%f%var(iwpsclrap2)%units = "unknown"
          k = k + 1
        case ('wpsclrbp2')
          iwpsclrbp2 = k
          zm%f%var(iwpsclrbp2)%ptr => zm%x(:,k)
          zm%f%var(iwpsclrbp2)%name = "wpsclrbp2"
          zm%f%var(iwpsclrbp2)%description 
     .    = "w'scalar(2)'^2 "
          zm%f%var(iwpsclrbp2)%units = "unknown"
          k = k + 1
        case ('wp2sclrap')
          iwp2sclrap = k
          zm%f%var(iwp2sclrap)%ptr => zm%x(:,k)
          zm%f%var(iwp2sclrap)%name = "wp2sclrap"
          zm%f%var(iwp2sclrap)%description 
     .    = "w'^2 scalar(1)'"
          zm%f%var(iwp2sclrap)%units = "unknown"
          k = k + 1
        case ('wp2sclrbp')
          iwp2sclrbp = k
          zm%f%var(iwp2sclrbp)%ptr => zm%x(:,k)
          zm%f%var(iwp2sclrbp)%name = "wp2sclrbp"
          zm%f%var(iwp2sclrbp)%description 
     .    = "w'^2 scalar(2)'"
          zm%f%var(iwp2sclrbp)%units = "unknown"
          k = k + 1
        case ('wpsclraprtp')
          iwpsclraprtp = k
          zm%f%var(iwpsclraprtp)%ptr => zm%x(:,k)
          zm%f%var(iwpsclraprtp)%name = "wpsclraprtp"
          zm%f%var(iwpsclraprtp)%description 
     .    = "w'scalar(1)'rt'"
          zm%f%var(iwpsclraprtp)%units = "unknown"
          k = k + 1
        case ('wpsclrbprtp')
          iwpsclrbprtp = k
          zm%f%var(iwpsclrbprtp)%ptr => zm%x(:,k)
          zm%f%var(iwpsclrbprtp)%name = "wpsclrbprtp"
          zm%f%var(iwpsclrbprtp)%description
     .    = "w'scalar(2)'rt'"
          zm%f%var(iwpsclrbprtp)%units = "unknown"
          k = k + 1
        case ('wpsclrapthlp')
          iwpsclrapthlp = k
          zm%f%var(iwpsclrapthlp)%ptr => zm%x(:,k)
          zm%f%var(iwpsclrapthlp)%name = "wpsclrapthlp"
          zm%f%var(iwpsclrapthlp)%description 
     .    = "w'scalar(1)'th_l'"
          zm%f%var(iwpsclrapthlp)%units = "unknown"
          k = k + 1
        case ('wpsclrbpthlp')
          iwpsclrbpthlp = k
          zm%f%var(iwpsclrbpthlp)%ptr => zm%x(:,k)
          zm%f%var(iwpsclrbpthlp)%name = "wpsclrbpthlp"
          zm%f%var(iwpsclrbpthlp)%description
     .    = "w'scalar(2)'th_l'"
          zm%f%var(iwpsclrbpthlp)%units = "unknown"
          k = k + 1
        case ('wpedsclrap')
          iwpedsclrap = k
          zm%f%var(iwpedsclrap)%ptr => zm%x(:,k)
          zm%f%var(iwpedsclrap)%name = "wpedsclrap"
          zm%f%var(iwpedsclrap)%description
     .    = "eddy scalar(1)'w'"
          zm%f%var(iwpedsclrap)%units = "unknown"
          k = k + 1
        case ('wpedsclrbp')
          iwpedsclrbp = k
          zm%f%var(iwpedsclrbp)%ptr => zm%x(:,k)
          zm%f%var(iwpedsclrbp)%name = "wpedsclrbp"
          zm%f%var(iwpedsclrbp)%description
     .    = "eddy scalar(2)'w'"
          zm%f%var(iwpedsclrbp)%units = "unknown"
          k = k + 1

#endif /*SCALARS*/

        case default
          write(*,*) 'Error: unrecognized variable ', trim(vars_zm(i))
          lerror = .true.

        end select

      end do

!     Initialize sfc

      i = 1
      do while ( ichar(vars_sfc(i)(1:1)) /= 0 
     .           .and. len_trim(vars_sfc(i)) /= 0
     .           .and. i <= nvarmax )
        i = i + 1
      end do
      ntot = i - 1
      if ( ntot == nvarmax )
     .  write(*,*) 'WARNING: check nvarmax in statistics.f'

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

!     Open GrADS files

      fdir = "./"
      fname = trim(fname_sfc)
      if ( lgrads ) then
        call open_grads( sfc%f, 20, fdir, fname, 
     .                   1, sfc%kk, sfc%z,
     .                   day, month, year, rlat, rlon,
     .                   time+stats_tout, stats_tout,
     .                   sfc%nn )
      else ! netcdf
        call open_nc( sfc%f, 20, fdir, fname, 
     .                1, sfc%kk, sfc%z,
     .                day, month, year, rlat, rlon,
     .                time+stats_tout, stats_tout,
     .                sfc%nn )
      end if

!     Default initialization for array indices for sfc

      iustar = 0
      ilh    = 0
      ish    = 0
      icc    = 0
      ilwp   = 0
      izb    = 0
      izi    = 0
      irain  = 0   ! Brian
      ipflux = 0   ! Brian

!     Assign pointers for statistics variables sfc

      k = 1
      do i=1,sfc%nn

        select case ( trim(vars_sfc(i)) )

        case ('ustar')
          iustar = k
          sfc%f%var(iustar)%ptr => sfc%x(:,k)
          sfc%f%var(iustar)%name = "ustar"
          sfc%f%var(iustar)%description = "ustar"
          sfc%f%var(iustar)%units = "m/s"
          k = k + 1

        case ('lh')
          ilh = k
          sfc%f%var(ilh)%ptr => sfc%x(:,k)
          sfc%f%var(ilh)%name = "lh"
          sfc%f%var(ilh)%description = "lh"
          sfc%f%var(ilh)%units = "W/m2"
          k = k + 1

        case ('sh')
          ish = k
          sfc%f%var(ish)%ptr => sfc%x(:,k)
          sfc%f%var(ish)%name = "sh"
          sfc%f%var(ish)%description = "sh"
          sfc%f%var(ish)%units = "W/m2"
          k = k + 1

        case ('cc')
          icc = k
          sfc%f%var(icc)%ptr => sfc%x(:,k)
          sfc%f%var(icc)%name = "cc"
          sfc%f%var(icc)%description = "cc"
          sfc%f%var(icc)%units = "count"
          k = k + 1

        case ('lwp')
          ilwp = k
          sfc%f%var(ilwp)%ptr => sfc%x(:,k)
          sfc%f%var(ilwp)%name = "lwp"
          sfc%f%var(ilwp)%description = "lwp"
          sfc%f%var(ilwp)%units = "kg/m2"
          k = k + 1

        case ('zb')
          izb = k
          sfc%f%var(izb)%ptr => sfc%x(:,k)
          sfc%f%var(izb)%name = "zb"
          sfc%f%var(izb)%description = "zb"
          sfc%f%var(izb)%units = "m"
          k = k + 1

        case ('zi')
          izi = k
          sfc%f%var(izi)%ptr => sfc%x(:,k)
          sfc%f%var(izi)%name = "zi"
          sfc%f%var(izi)%description = "zi"
          sfc%f%var(izi)%units = "m"
          k = k + 1

        case ('rain')          ! Brian
          irain = k
          sfc%f%var(irain)%ptr => sfc%x(:,k)
          sfc%f%var(irain)%name = "rain_rate"
          sfc%f%var(irain)%description = "rainfall rate (mm/day)"
          sfc%f%var(irain)%units = "mm/day"
          k = k + 1

        case ('pflux')         ! Brian
          ipflux = k
          sfc%f%var(ipflux)%ptr => sfc%x(:,k)
          sfc%f%var(ipflux)%name = "prec_flux"
          sfc%f%var(ipflux)%description = "Precipitation flux (W/m^2)"
          sfc%f%var(ipflux)%units = "W/m2"
          k = k + 1

        case default
          write(*,*) 'Error: unrecognized variable ',trim(vars_sfc(i))
          lerror = .true.

        end select

      end do

!     Check for errors

      if ( lerror ) then
        write(*,*) 'stats_init: errors found'
        stop
      endif

      return

!     If namelist was not found in input file, turn off statistics

100   continue
      write(*,*) 'Error with statsnl, statistics is turned off'
      lstats       = .false.
      lstats_samp  = .false.
      lstats_first = .false.
      lstats_last  = .false.

      return
      end subroutine stats_init

c     -----------------------------------------------------------------
      subroutine stats_zero( kk, nn, x, n )

      implicit none

!     Input
      integer, intent(in) :: kk, nn

!     Output
      real, dimension(kk,nn), intent(out)    :: x
      integer, dimension(kk,nn), intent(out) :: n

!     Zero out arrays

      if ( nn > 0 ) then
        x = 0.0
        n = 0
      endif

      return
      end subroutine stats_zero

c     -----------------------------------------------------------------
      subroutine stats_avg( kk, nn, x, n )

c     use stats_type
      implicit none

!     Input/output

!     type (stats), intent(inout) :: s

!     Input
      integer, intent(in) :: nn, kk

!     Output
      real, dimension(kk,nn), intent(inout)  :: x
      integer, dimension(kk,nn), intent(in)  :: n

!     Internal

      integer k,m

!     Compute averages

      do m=1,nn
        do k=1,kk

          if ( n(k,m) > 0 ) then
            x(k,m) = x(k,m) / n(k,m)
          end if

        enddo
      enddo

      return
      end subroutine stats_avg

c     -----------------------------------------------------------------
      subroutine stats_begin_timestep( time, delt )

      use statistics

      implicit none

!     Input

      real, intent(in) :: time, delt

      if ( .not. lstats ) return

!     Set sample this time step flag
      if ( mod( time, stats_tsamp ) < 1.e-8 ) then
         lstats_samp = .true.
      else
         lstats_samp = .false.
      endif

!     Set first time step flag

      if ( mod( time - delt, stats_tout ) < 1.e-8 ) then
         lstats_first = .true.
      else
         lstats_first = .false.
      endif

!     Set last time step flag

      if ( mod( time, stats_tout ) < 1.e-8 ) then
         lstats_last = .true.
      else
         lstats_last = .false.
      endif

      return

      end subroutine stats_begin_timestep

c     -----------------------------------------------------------------
      subroutine stats_end_timestep( time, delt )

      use statistics
      use outputgrads
      use outputnc

      implicit none

!     Input

      real, intent(in) :: time, delt

!     Internal

      integer i, k

!     Check if it is time to write to file

      if ( .not. lstats_last ) return

!     Check number of sampling points

      do i=1,zt%nn
       do k=1,zt%kk
         if ( zt%n(k,i) /= 0 
     .        .and. zt%n(k,i) /= int(stats_tout/stats_tsamp) ) then
           write(*,*) 'Possible sampling error for zt ',i,k,zt%n(k,i)
         endif
        enddo
      enddo
      
      do i=1,zm%nn
       do k=1,zm%kk
         if ( zm%n(k,i) /= 0 
     .        .and. zm%n(k,i) /= int(stats_tout/stats_tsamp) ) then
           write(*,*) 'Possible sampling error for zm ',i,k,zm%n(k,i)
         endif
        enddo
      enddo
      
!     Compute averages

      call stats_avg( zt%kk, zt%nn, zt%x, zt%n )
      call stats_avg( zm%kk, zm%nn, zm%x, zm%n )
      call stats_avg( sfc%kk, sfc%nn, sfc%x, sfc%n )

!     Write to file
      if ( lgrads ) then
        call write_grads( zt%f  )
        call write_grads( zm%f  )
        call write_grads( sfc%f  )
      else
        call write_nc( zt%f  )
        call write_nc( zm%f  )
        call write_nc( sfc%f  )
      endif

!     Reset sample fields
      call stats_zero( zt%kk, zt%nn, zt%x, zt%n )
      call stats_zero( zm%kk, zm%nn, zm%x, zm%n )
      call stats_zero( sfc%kk, sfc%nn, sfc%x, sfc%n )

      return
      end subroutine stats_end_timestep

c     -----------------------------------------------------------------
!     subroutine stats_accumulate( time, delt )
      subroutine stats_accumulate( )

      use statistics
      use grid_class
      use prognostic_variables
      use diagnostic_variables
      use model_flags

      implicit none

!     Input

!     real, intent(in) :: time, delt

!     Internal

      integer i, k
      real xtmp

!     Sample fields

      if ( lstats_samp ) then

!       zt variables

        if ( ithlm > 0 ) then
          zt%x(:,ithlm) = zt%x(:,ithlm) + thlm
          zt%n(:,ithlm) = zt%n(:,ithlm) + 1
        endif
        if ( ithvm > 0 ) then
          zt%x(:,ithvm) = zt%x(:,ithvm) + thvm
          zt%n(:,ithvm) = zt%n(:,ithvm) + 1
        endif
        if ( irtm > 0 ) then
          zt%x(:,irtm) = zt%x(:,irtm) + rtm
          zt%n(:,irtm) = zt%n(:,irtm) + 1
        endif
        if ( ircm > 0 ) then
          zt%x(:,ircm) = zt%x(:,ircm) + rcm
          zt%n(:,ircm) = zt%n(:,ircm) + 1
        endif
        if ( ium > 0 ) then
          zt%x(:,ium) = zt%x(:,ium) + um
          zt%n(:,ium) = zt%n(:,ium) + 1
        endif
        if ( ivm > 0 ) then
          zt%x(:,ivm) = zt%x(:,ivm) + vm
          zt%n(:,ivm) = zt%n(:,ivm) + 1
        endif
        if ( iwmt > 0 ) then
          zt%x(:,iwmt) = zt%x(:,iwmt) + wmt
          zt%n(:,iwmt) = zt%n(:,iwmt) + 1
        endif
        if ( iug > 0 ) then
          zt%x(:,iug) = zt%x(:,iug) + ug
          zt%n(:,iug) = zt%n(:,iug) + 1
        endif
        if ( ivg > 0 ) then
          zt%x(:,ivg) = zt%x(:,ivg) + vg
          zt%n(:,ivg) = zt%n(:,ivg) + 1
        endif
        if ( icf > 0 ) then
          zt%x(:,icf) = zt%x(:,icf) + cf
          zt%n(:,icf) = zt%n(:,icf) + 1
        endif
        if ( ip > 0 ) then
          zt%x(:,ip) = zt%x(:,ip) + p
          zt%n(:,ip) = zt%n(:,ip) + 1
        endif
        if ( iLscale > 0 ) then
          zt%x(:,iLscale) = zt%x(:,iLscale) + Lscale
          zt%n(:,iLscale) = zt%n(:,iLscale) + 1
        endif
        if ( ithlm_forcing > 0 ) then
          zt%x(:,ithlm_forcing) = zt%x(:,ithlm_forcing) + thlm_forcing
          zt%n(:,ithlm_forcing) = zt%n(:,ithlm_forcing) + 1
        endif
        if ( irtm_forcing > 0 ) then
          zt%x(:,irtm_forcing) = zt%x(:,irtm_forcing) + rtm_forcing
          zt%n(:,irtm_forcing) = zt%n(:,irtm_forcing) + 1
        endif
        if ( iwp3 > 0 ) then
          zt%x(:,iwp3) = zt%x(:,iwp3) + wp3
          zt%n(:,iwp3) = zt%n(:,iwp3) + 1
        endif
        if ( iwpthlp2 > 0 ) then
          zt%x(:,iwpthlp2) = zt%x(:,iwpthlp2) + wpthlp2
          zt%n(:,iwpthlp2) = zt%n(:,iwpthlp2) + 1
        endif
        if ( iwp2thlp > 0 ) then
          zt%x(:,iwp2thlp) = zt%x(:,iwp2thlp) + wp2thlp
          zt%n(:,iwp2thlp) = zt%n(:,iwp2thlp) + 1
        endif
        if ( iwprtp2 > 0 ) then
          zt%x(:,iwprtp2) = zt%x(:,iwprtp2) + wprtp2
          zt%n(:,iwprtp2) = zt%n(:,iwprtp2) + 1
        endif
        if ( iwp2rtp > 0 ) then
          zt%x(:,iwp2rtp) = zt%x(:,iwp2rtp) + wp2rtp
          zt%n(:,iwp2rtp) = zt%n(:,iwp2rtp) + 1
        endif
        if ( ilup > 0 ) then
          zt%x(:,ilup) = zt%x(:,ilup) + lup
          zt%n(:,ilup) = zt%n(:,ilup) + 1
        endif
        if ( ildown > 0 ) then
          zt%x(:,ildown) = zt%x(:,ildown) + ldown
          zt%n(:,ildown) = zt%n(:,ildown) + 1
        endif
        if ( itaut > 0 ) then
          zt%x(:,itaut) = zt%x(:,itaut) + taut
          zt%n(:,itaut) = zt%n(:,itaut) + 1
        endif
        if ( ikht > 0 ) then
          zt%x(:,ikht) = zt%x(:,ikht) + kht
          zt%n(:,ikht) = zt%n(:,ikht) + 1
        endif
        if ( iwp2thvp > 0 ) then
          zt%x(:,iwp2thvp) = zt%x(:,iwp2thvp) + wp2thvp
          zt%n(:,iwp2thvp) = zt%n(:,iwp2thvp) + 1
        endif
        if ( iwp2rcp > 0 ) then
          zt%x(:,iwp2rcp) = zt%x(:,iwp2rcp) + wp2rcp
          zt%n(:,iwp2rcp) = zt%n(:,iwp2rcp) + 1
        endif
        if ( iwprtpthlp > 0 ) then
          zt%x(:,iwprtpthlp) = zt%x(:,iwprtpthlp) + wprtpthlp
          zt%n(:,iwprtpthlp) = zt%n(:,iwprtpthlp) + 1
        endif
        if ( isct > 0 ) then
          zt%x(:,isct) = zt%x(:,isct) + sct
          zt%n(:,isct) = zt%n(:,isct) + 1
        endif
        if ( irhot > 0 ) then
          zt%x(:,irhot) = zt%x(:,irhot) + rhot
          zt%n(:,irhot) = zt%n(:,irhot) + 1
        endif
        if ( iNcm > 0 ) then           ! Brian
          zt%x(:,iNcm) = zt%x(:,iNcm) + Ncm
          zt%n(:,iNcm) = zt%n(:,iNcm) + 1
        endif
        if ( ised_rcm > 0 .and. cloud_sed ) then       ! Brian
          zt%x(:,ised_rcm) = zt%x(:,ised_rcm) + sed_rcm
          zt%n(:,ised_rcm) = zt%n(:,ised_rcm) + 1
        endif
        if ( irsm > 0 ) then           ! Brian
          zt%x(:,irsm) = zt%x(:,irsm) + rsm
          zt%n(:,irsm) = zt%n(:,irsm) + 1
        endif
        if ( irrm > 0 .and. kk_rain ) then     ! Brian
          zt%x(:,irrm) = zt%x(:,irrm) + rrm
          zt%n(:,irrm) = zt%n(:,irrm) + 1
        endif
        if ( iNrm > 0 .and. kk_rain ) then     ! Brian
          zt%x(:,iNrm) = zt%x(:,iNrm) + Nrm
          zt%n(:,iNrm) = zt%n(:,iNrm) + 1
        endif
        if ( iVrr > 0 .and. kk_rain ) then      ! Brian
          zt%x(:,iVrr) = zt%x(:,iVrr) + Vrr
          zt%n(:,iVrr) = zt%n(:,iVrr) + 1
        endif
        if ( imean_vol_rad > 0 .and. kk_rain ) then  ! Brian
          zt%x(:,imean_vol_rad) = zt%x(:,imean_vol_rad) + mean_vol_rad
          zt%n(:,imean_vol_rad) = zt%n(:,imean_vol_rad) + 1
        endif
        if ( irain_rate > 0 .and. kk_rain ) then     ! Brian
          zt%x(:,irain_rate) = zt%x(:,irain_rate) + rain_rate
          zt%n(:,irain_rate) = zt%n(:,irain_rate) + 1
        endif

        if ( iAKm > 0 .and. LH_on ) then         ! Vince Larson 22 May 2005
          zt%x(:,iAKm) = zt%x(:,iAKm) + AKm
          zt%n(:,iAKm) = zt%n(:,iAKm) + 1
        endif

        if ( iAKm_est > 0 .and. LH_on ) then     ! Vince Larson 22 May 2005
          zt%x(:,iAKm_est) = zt%x(:,iAKm_est) + AKm_est
          zt%n(:,iAKm_est) = zt%n(:,iAKm_est) + 1
        endif

        if ( iradht > 0 ) then
          zt%x(:,iradht) = zt%x(:,iradht) + radht
          zt%n(:,iradht) = zt%n(:,iradht) + 1
        endif

        if ( iradht_LW > 0 ) then
          zt%x(:,iradht_LW) = zt%x(:,iradht_LW) + radht_LW
          zt%n(:,iradht_LW) = zt%n(:,iradht_LW) + 1
        endif

        if ( iradht_SW > 0 ) then
          zt%x(:,iradht_SW) = zt%x(:,iradht_SW) + radht_SW
          zt%n(:,iradht_SW) = zt%n(:,iradht_SW) + 1
        endif

        if ( ia > 0 ) then
          zt%x(:,ia) = zt%x(:,ia) + pdf_parms(:,13)
          zt%n(:,ia) = zt%n(:,ia) + 1
        endif

        if ( iw1 > 0 ) then
          zt%x(:,iw1) = zt%x(:,iw1) + pdf_parms(:,1)
          zt%n(:,iw1) = zt%n(:,iw1) + 1
        endif

        if ( iw2 > 0 ) then
          zt%x(:,iw2) = zt%x(:,iw2) + pdf_parms(:,2)
          zt%n(:,iw2) = zt%n(:,iw2) + 1
        endif

        if ( isw1 > 0 ) then
          zt%x(:,isw1) = zt%x(:,isw1) + pdf_parms(:,3)
          zt%n(:,isw1) = zt%n(:,isw1) + 1
        endif

        if ( isw2 > 0 ) then
          zt%x(:,isw2) = zt%x(:,isw2) + pdf_parms(:,4)
          zt%n(:,isw2) = zt%n(:,isw2) + 1
        endif

        if ( ithl1 > 0 ) then
          zt%x(:,ithl1) = zt%x(:,ithl1) + pdf_parms(:,9)
          zt%n(:,ithl1) = zt%n(:,ithl1) + 1
        endif

        if ( ithl2 > 0 ) then
          zt%x(:,ithl2) = zt%x(:,ithl2) + pdf_parms(:,10)
          zt%n(:,ithl2) = zt%n(:,ithl2) + 1
        endif

        if ( isthl1 > 0 ) then
          zt%x(:,isthl1) = zt%x(:,isthl1) + pdf_parms(:,11)
          zt%n(:,isthl1) = zt%n(:,isthl1) + 1
        endif

        if ( isthl2 > 0 ) then
          zt%x(:,isthl2) = zt%x(:,isthl2) + pdf_parms(:,12)
          zt%n(:,isthl2) = zt%n(:,isthl2) + 1
        endif

        if ( irt1 > 0 ) then
          zt%x(:,irt1) = zt%x(:,irt1) + pdf_parms(:,5)
          zt%n(:,irt1) = zt%n(:,irt1) + 1
        endif

        if ( irt2 > 0 ) then
          zt%x(:,irt2) = zt%x(:,irt2) + pdf_parms(:,6)
          zt%n(:,irt2) = zt%n(:,irt2) + 1
        endif

        if ( isrt1 > 0 ) then
          zt%x(:,isrt1) = zt%x(:,isrt1) + pdf_parms(:,7)
          zt%n(:,isrt1) = zt%n(:,isrt1) + 1
        endif

        if ( isrt2 > 0 ) then
          zt%x(:,isrt2) = zt%x(:,isrt2) + pdf_parms(:,8)
          zt%n(:,isrt2) = zt%n(:,isrt2) + 1
        endif

        if ( irc1 > 0 ) then
          zt%x(:,irc1) = zt%x(:,irc1) + pdf_parms(:,14)
          zt%n(:,irc1) = zt%n(:,irc1) + 1
        endif

        if ( irc2 > 0 ) then
          zt%x(:,irc2) = zt%x(:,irc2) + pdf_parms(:,15)
          zt%n(:,irc2) = zt%n(:,irc2) + 1
        endif

        if ( irsl1 > 0 ) then
          zt%x(:,irsl1) = zt%x(:,irsl1) + pdf_parms(:,16)
          zt%n(:,irsl1) = zt%n(:,irsl1) + 1
        endif

        if ( irsl2 > 0 ) then
          zt%x(:,irsl2) = zt%x(:,irsl2) + pdf_parms(:,17)
          zt%n(:,irsl2) = zt%n(:,irsl2) + 1
        endif

        if ( iR1 > 0 ) then
          zt%x(:,iR1) = zt%x(:,iR1) + pdf_parms(:,18)
          zt%n(:,iR1) = zt%n(:,iR1) + 1
        endif

        if ( iR2 > 0 ) then
          zt%x(:,iR2) = zt%x(:,iR2) + pdf_parms(:,19)
          zt%n(:,iR2) = zt%n(:,iR2) + 1
        endif

        if ( is1 > 0 ) then
          zt%x(:,is1) = zt%x(:,is1) + pdf_parms(:,20)
          zt%n(:,is1) = zt%n(:,is1) + 1
        endif

        if ( is2 > 0 ) then
          zt%x(:,is2) = zt%x(:,is2) + pdf_parms(:,21)
          zt%n(:,is2) = zt%n(:,is2) + 1
        endif

        if ( iss1 > 0 ) then
          zt%x(:,iss1) = zt%x(:,iss1) + pdf_parms(:,22)
          zt%n(:,iss1) = zt%n(:,iss1) + 1
        endif

        if ( iss2 > 0 ) then
          zt%x(:,iss2) = zt%x(:,iss2) + pdf_parms(:,23)
          zt%n(:,iss2) = zt%n(:,iss2) + 1
        endif

        if ( irrtthl > 0 ) then
          zt%x(:,irrtthl) = zt%x(:,irrtthl) + pdf_parms(:,24)
          zt%n(:,irrtthl) = zt%n(:,irrtthl) + 1
        endif

        if ( iwp2zt > 0 ) then
          zt%x(:,iwp2zt) = zt%x(:,iwp2zt) + wp2zt
          zt%n(:,iwp2zt) = zt%n(:,iwp2zt) + 1
        endif

        if ( ithlp2zt > 0 ) then
          zt%x(:,ithlp2zt) = zt%x(:,ithlp2zt) + thlp2zt
          zt%n(:,ithlp2zt) = zt%n(:,ithlp2zt) + 1
        endif

        if ( iwpthlpzt > 0 ) then
          zt%x(:,iwpthlpzt) = zt%x(:,iwpthlpzt) + wpthlpzt
          zt%n(:,iwpthlpzt) = zt%n(:,iwpthlpzt) + 1
        endif

        if ( iwprtpzt > 0 ) then
          zt%x(:,iwprtpzt) = zt%x(:,iwprtpzt) + wprtpzt
          zt%n(:,iwprtpzt) = zt%n(:,iwprtpzt) + 1
        endif

        if ( irtp2zt > 0 ) then
          zt%x(:,irtp2zt) = zt%x(:,irtp2zt) + rtp2zt
          zt%n(:,irtp2zt) = zt%n(:,irtp2zt) + 1
        endif

        if ( irtpthlpzt > 0 ) then
          zt%x(:,irtpthlpzt) = zt%x(:,irtpthlpzt) + rtpthlpzt
          zt%n(:,irtpthlpzt) = zt%n(:,irtpthlpzt) + 1
        endif
#ifdef SCALARS
        if ( isclram > 0 ) then
          zt%x(:,isclram) = zt%x(:,isclram) + sclrm(:,1)
          zt%n(:,isclram) = zt%n(:,isclram) + 1
        endif

        if ( isclram_f > 0 ) then
          zt%x(:,isclram_f) = zt%x(:,isclram_f) + sclrm_forcing(:,1)
          zt%n(:,isclram_f) = zt%n(:,isclram_f) + 1
        endif

        if ( isclrbm > 0 ) then
          zt%x(:,isclrbm) = zt%x(:,isclrbm) + sclrm(:,2)
          zt%n(:,isclrbm) = zt%n(:,isclrbm) + 1
        endif

        if ( isclrbm_f > 0 ) then
          zt%x(:,isclrbm_f) = zt%x(:,isclrbm_f) + sclrm_forcing(:,2)
          zt%n(:,isclrbm_f) = zt%n(:,isclrbm_f) + 1
        endif

        if ( iedsclram > 0 ) then
          zt%x(:,iedsclram) = zt%x(:,iedsclram) + edsclrm(:,1)
          zt%n(:,iedsclram) = zt%n(:,iedsclram) + 1
        endif

        if ( iedsclrbm > 0 ) then
          zt%x(:,iedsclrbm) = zt%x(:,iedsclrbm) + edsclrm(:,2)
          zt%n(:,iedsclrbm) = zt%n(:,iedsclrbm) + 1
        endif

#endif /*SCALARS*/

!       zm variables

        if ( iwp2 > 0 ) then
          zm%x(:,iwp2) = zm%x(:,iwp2) + wp2
          zm%n(:,iwp2) = zm%n(:,iwp2) + 1
        endif

        if ( irtp2 > 0 ) then
          zm%x(:,irtp2) = zm%x(:,irtp2) + rtp2
          zm%n(:,irtp2) = zm%n(:,irtp2) + 1
        endif

        if ( ithlp2 > 0 ) then
          zm%x(:,ithlp2) = zm%x(:,ithlp2) + thlp2
          zm%n(:,ithlp2) = zm%n(:,ithlp2) + 1
        endif

        if ( irtpthlp > 0 ) then
          zm%x(:,irtpthlp) = zm%x(:,irtpthlp) + rtpthlp
          zm%n(:,irtpthlp) = zm%n(:,irtpthlp) + 1
        endif

        if ( iwprtp > 0 ) then
          zm%x(:,iwprtp) = zm%x(:,iwprtp) + wprtp
          zm%n(:,iwprtp) = zm%n(:,iwprtp) + 1
        endif

        if ( iwpthlp > 0 ) then
          zm%x(:,iwpthlp) = zm%x(:,iwpthlp) + wpthlp
          zm%n(:,iwpthlp) = zm%n(:,iwpthlp) + 1
        endif

        if ( iwp4 > 0 ) then
          zm%x(:,iwp4) = zm%x(:,iwp4) + wp4
          zm%n(:,iwp4) = zm%n(:,iwp4) + 1
        endif

        if ( iwpthvp > 0 ) then
          zm%x(:,iwpthvp) = zm%x(:,iwpthvp) + wpthvp
          zm%n(:,iwpthvp) = zm%n(:,iwpthvp) + 1
        endif

        if ( irtpthvp > 0 ) then
          zm%x(:,irtpthvp) = zm%x(:,irtpthvp) + rtpthvp
          zm%n(:,irtpthvp) = zm%n(:,irtpthvp) + 1
        endif

        if ( ithlpthvp > 0 ) then
          zm%x(:,ithlpthvp) = zm%x(:,ithlpthvp) + thlpthvp
          zm%n(:,ithlpthvp) = zm%n(:,ithlpthvp) + 1
        endif

        if ( itaum > 0 ) then
          zm%x(:,itaum) = zm%x(:,itaum) + taum
          zm%n(:,itaum) = zm%n(:,itaum) + 1
        endif

        if ( ikhm > 0 ) then
          zm%x(:,ikhm) = zm%x(:,ikhm) + khm
          zm%n(:,ikhm) = zm%n(:,ikhm) + 1
        endif

        if ( iwprcp > 0 ) then
          zm%x(:,iwprcp) = zm%x(:,iwprcp) + wprcp
          zm%n(:,iwprcp) = zm%n(:,iwprcp) + 1
        endif

        if ( ithlprcp > 0 ) then
          zm%x(:,ithlprcp) = zm%x(:,ithlprcp) + thlprcp
          zm%n(:,ithlprcp) = zm%n(:,ithlprcp) + 1
        endif

        if ( irtprcp > 0 ) then
          zm%x(:,irtprcp) = zm%x(:,irtprcp) + rtprcp
          zm%n(:,irtprcp) = zm%n(:,irtprcp) + 1
        endif

        if ( iupwp > 0 ) then
          zm%x(:,iupwp) = zm%x(:,iupwp) + upwp
          zm%n(:,iupwp) = zm%n(:,iupwp) + 1
        endif

        if ( ivpwp > 0 ) then
          zm%x(:,ivpwp) = zm%x(:,ivpwp) + vpwp
          zm%n(:,ivpwp) = zm%n(:,ivpwp) + 1
        endif

        if ( irhom > 0 ) then
          zm%x(:,irhom) = zm%x(:,irhom) + rhom
          zm%n(:,irhom) = zm%n(:,irhom) + 1
        endif

        if ( iscm > 0 ) then
          zm%x(:,iscm) = zm%x(:,iscm) + scm
          zm%n(:,iscm) = zm%n(:,iscm) + 1
        endif

        if ( iem > 0 ) then
          zm%x(:,iem) = zm%x(:,iem) + em
          zm%n(:,iem) = zm%n(:,iem) + 1
        endif

        if ( ishear > 0 ) then    ! Brian
          zm%x(:,ishear) = zm%x(:,ishear) + shear
          zm%n(:,ishear) = zm%n(:,ishear) + 1
        endif

        if ( iFrad > 0 ) then
          zm%x(:,iFrad) = zm%x(:,iFrad) + Frad
          zm%n(:,iFrad) = zm%n(:,iFrad) + 1
        endif

        if ( iFrad_LW > 0 ) then  ! Brian
          zm%x(:,iFrad_LW) = zm%x(:,iFrad_LW) + Frad_LW
          zm%n(:,iFrad_LW) = zm%n(:,iFrad_LW) + 1
        endif

        if ( iFrad_SW > 0 ) then  ! Brian
          zm%x(:,iFrad_SW) = zm%x(:,iFrad_SW) + Frad_SW
          zm%n(:,iFrad_SW) = zm%n(:,iFrad_SW) + 1
        endif

        if ( iFprec > 0 .and. kk_rain ) then    ! Brian
          zm%x(:,iFprec) = zm%x(:,iFprec) + Fprec
          zm%n(:,iFprec) = zm%n(:,iFprec) + 1
        endif

        if ( iFcsed > 0 .and. cloud_sed ) then    ! Brian
          zm%x(:,iFcsed) = zm%x(:,iFcsed) + Fcsed
          zm%n(:,iFcsed) = zm%n(:,iFcsed) + 1
        endif

#ifdef SCALARS
        if ( isclraprtp > 0 ) then
          zm%x(:,isclraprtp) = zm%x(:,isclraprtp) + sclrprtp(:,1)
          zm%n(:,isclraprtp) = zm%n(:,isclraprtp) + 1
        endif
        if ( isclrbprtp > 0 ) then
          zm%x(:,isclrbprtp) = zm%x(:,isclrbprtp) + sclrprtp(:,2)
          zm%n(:,isclrbprtp) = zm%n(:,isclrbprtp) + 1
        endif
        if ( isclrap2 > 0 ) then
          zm%x(:,isclrap2) = zm%x(:,isclrap2) + sclrp2(:,1)
          zm%n(:,isclrap2) = zm%n(:,isclrap2) + 1
        endif
        if ( isclrbp2 > 0 ) then
          zm%x(:,isclrbp2) = zm%x(:,isclrbp2) + sclrp2(:,2)
          zm%n(:,isclrbp2) = zm%n(:,isclrbp2) + 1
        endif
        if ( isclrapthvp > 0 ) then
          zm%x(:,isclrapthvp) = zm%x(:,isclrapthvp) 
     .                            + sclrpthvp(:,1)
          zm%n(:,isclrapthvp) = zm%n(:,isclrapthvp) + 1
        endif
        if ( isclrbpthvp > 0 ) then
          zm%x(:,isclrbpthvp) = zm%x(:,isclrbpthvp) 
     .                            + sclrpthvp(:,2)
          zm%n(:,isclrbpthvp) = zm%n(:,isclrbpthvp) + 1
        endif
        if ( isclrapthlp > 0 ) then
          zm%x(:,isclrapthlp) = zm%x(:,isclrapthlp) 
     .                            + sclrpthlp(:,1)
          zm%n(:,isclrapthlp) = zm%n(:,isclrapthlp) + 1
        endif
        if ( isclrbpthlp > 0 ) then
          zm%x(:,isclrbpthlp) = zm%x(:,isclrbpthlp) 
     .                            + sclrpthlp(:,2)
          zm%n(:,isclrbpthlp) = zm%n(:,isclrbpthlp) + 1
        endif
        if ( isclraprcp > 0 ) then
          zm%x(:,isclraprcp) = zm%x(:,isclraprcp) + sclrprcp(:,1)
          zm%n(:,isclraprcp) = zm%n(:,isclraprcp) + 1
        endif
        if ( isclrbprcp > 0 ) then
          zm%x(:,isclrbprcp) = zm%x(:,isclrbprcp) + sclrprcp(:,2)
          zm%n(:,isclrbprcp) = zm%n(:,isclrbprcp) + 1
        endif
        if ( iwpsclrap > 0 ) then
          zm%x(:,iwpsclrap) = zm%x(:,iwpsclrap) + wpsclrp(:,1)
          zm%n(:,iwpsclrap) = zm%n(:,iwpsclrap) + 1
        endif
        if ( iwpsclrbp > 0 ) then
          zm%x(:,iwpsclrbp) = zm%x(:,iwpsclrbp) + wpsclrp(:,2)
          zm%n(:,iwpsclrbp) = zm%n(:,iwpsclrbp) + 1
        endif
        if ( iwp2sclrap > 0 ) then
          zm%x(:,iwp2sclrap) = zm%x(:,iwp2sclrap) + wp2sclrp(:,1)
          zm%n(:,iwp2sclrap) = zm%n(:,iwp2sclrap) + 1
        endif
        if ( iwp2sclrbp > 0 ) then
          zm%x(:,iwp2sclrbp) = zm%x(:,iwp2sclrbp) + wp2sclrp(:,2)
          zm%n(:,iwp2sclrbp) = zm%n(:,iwp2sclrbp) + 1
        endif
        if ( iwpsclrap2 > 0 ) then
          zm%x(:,iwpsclrap2) = zm%x(:,iwpsclrap2) + wpsclrp2(:,1)
          zm%n(:,iwpsclrap2) = zm%n(:,iwpsclrap2) + 1
        endif
        if ( iwpsclrbp2 > 0 ) then
          zm%x(:,iwpsclrbp2) = zm%x(:,iwpsclrbp2) + wpsclrp2(:,2)
          zm%n(:,iwpsclrbp2) = zm%n(:,iwpsclrbp2) + 1
        endif
        if ( iwpsclraprtp > 0 ) then
          zm%x(:,iwpsclraprtp) = zm%x(:,iwpsclraprtp) + wpsclrprtp(:,1)
          zm%n(:,iwpsclraprtp) = zm%n(:,iwpsclraprtp) + 1
        endif
        if ( iwpsclrbprtp > 0 ) then
          zm%x(:,iwpsclrbprtp) = zm%x(:,iwpsclrbprtp) + wpsclrprtp(:,2)
          zm%n(:,iwpsclrbprtp) = zm%n(:,iwpsclrbprtp) + 1
        endif
        if ( iwpsclrapthlp > 0 ) then
          zm%x(:,iwpsclrapthlp) = zm%x(:,iwpsclrapthlp) 
     .                            + wpsclrpthlp(:,1)
          zm%n(:,iwpsclrapthlp) = zm%n(:,iwpsclrapthlp) + 1
        endif
        if ( iwpsclrbpthlp > 0 ) then
          zm%x(:,iwpsclrbpthlp) = zm%x(:,iwpsclrbpthlp) 
     .                            + wpsclrpthlp(:,2)
          zm%n(:,iwpsclrbpthlp) = zm%n(:,iwpsclrbpthlp) + 1
        endif
        if ( iwpedsclrap > 0 ) then
          zm%x(:,iwpedsclrap) = zm%x(:,iwpedsclrap) 
     .                            + wpedsclrp(:,1)
          zm%n(:,iwpedsclrap) = zm%n(:,iwpedsclrap) + 1
        endif
        if ( iwpedsclrbp > 0 ) then
          zm%x(:,iwpedsclrbp) = zm%x(:,iwpedsclrbp) 
     .                            + wpedsclrp(:,2)
          zm%n(:,iwpedsclrbp) = zm%n(:,iwpedsclrbp) + 1
        endif
        
#endif /*SCALARS*/

! sfc variables

!       Cloud cover
        if ( icc > 0 ) then
          sfc%x(1,icc) = sfc%x(1,icc) + maxval( cf(1:gr%nnzp) )
          sfc%n(1,icc) = sfc%n(1,icc) + 1
        endif

!       Cloud base
        if ( izb > 0 ) then
          k = 1
          do while ( rcm(k) < 1.e-6 .and. k < gr%nnzp )
            k = k + 1
          end do

!          if ( k < gr%nnzp) then
!            sfc%x(1,izb) = sfc%x(1,izb) + gr%zt(k)
          if ( k > 1 .AND. k < gr%nnzp) then
            sfc%x(1,izb) = sfc%x(1,izb) + (                     ! Use linear interpolation
     .         ( ( 1.0E-6 - rcm(k-1) )/( rcm(k) - rcm(k-1) ) )  ! to find the exact height
     .        * ( gr%zt(k) - gr%zt(k-1) ) + gr%zt(k-1)          ! of the 1.0E-6 kg/kg level.
     .                                    )                     ! Brian.
          else
            sfc%x(1,izb) = sfc%x(1,izb) - 10
          endif
          sfc%n(1,izb) = sfc%n(1,izb) + 1
        endif

!       LWP
        if ( ilwp > 0 ) then
          xtmp = 0.
          do i = gr%nnzp-1, 1, -1
           xtmp = xtmp + rhot(i+1) * rcm(i+1) / gr%dzt(i+1)
          enddo
          sfc%x(1,ilwp) = sfc%x(1,ilwp) + xtmp
          sfc%n(1,ilwp) = sfc%n(1,ilwp) + 1
        end if

      endif

      return
      end subroutine stats_accumulate

c     -----------------------------------------------------------------
      subroutine stats_finalize

      use statistics
      use outputnc

      implicit none

      if ( lstats .and. lnetcdf ) then
        call close_nc( zt%f )
        call close_nc( zm%f )
        call close_nc( sfc%f )
      end if

      if ( lstats ) then
!       De-allocate all zt variables
        deallocate( zt%z )

        deallocate( zt%x )
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

!       De-allocate all zm variables
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

!       De-allocate all sfc variables
        deallocate( sfc%z )

        deallocate( sfc%x )
        deallocate( sfc%n )

        deallocate( sfc%f%var )
        deallocate( sfc%f%z )
      endif ! lstats

      return
      end subroutine stats_finalize

#endif /*STATS*/
