#ifdef STATS

      module stats_type

        use outputfile_class

c       Derived data types to store statistics

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

      module statistics

        use stats_type

c       Sampling and output frequencies

        logical lstats         ! Main flag to turn statistics on/off
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
     .     irsm,           ! Brian
     .     irrm,           ! Brian
     .     iNrm,           ! Brian
     .     icond_rrm,      ! Brian
     .     iauto_rrm,      ! Brian
     .     iaccr_rrm,      ! Brian
     .     icond_Nrm,      ! Brian
     .     iauto_Nrm,      ! Brian
     .     ised_rrm,       ! Brian
     .     idiff_rrm,      ! Brian
     .     isubs_rrm,      ! Brian
     .     iVr,            ! Brian
     .     imean_vol_rad,  ! Brian
     .     irain_rate,     ! Brian
     .     iF_precip,       ! Brian
     .     iAKm,            ! analytic Kessler.  Vince Larson 22 May 2005 
     .     iAKm_est        ! LH Kessler.  Vince Larson  22 May 2005


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
     .     ir1,
     .     ir2,
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
     .     irtp2zt

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
     .     ifrad,
     .     iscm,
     .     iem,
     .     ishear     ! Brian

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

        integer
     .     iwpqtp_bt,
     .     iwpqtp_ma,
     .     iwpqtp_ta,
     .     iwpqtp_tp,
     .     iwpqtp_ac,
     .     iwpqtp_bp,
     .     iwpqtp_pr1,
     .     iwpqtp_pr2,
     .     iwpqtp_pr3,
     .     iwpqtp_dp1

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
     .                       ztscr10(:), ztscr11(:)

        real, allocatable :: zmscr01(:), zmscr02(:), zmscr03(:), 
     .                       zmscr04(:), zmscr05(:), zmscr06(:),
     .                       zmscr07(:), zmscr08(:), zmscr09(:),
     .                       zmscr10(:), zmscr11(:), zmscr12(:),
     .                       zmscr13(:)

      end module statistics

c     -----------------------------------------------------------------
      subroutine stats_init(iunit,fnamelist,nnzp,gzt,gzm,
     .                      day,month,year,time,delt)

      use statistics

      implicit none

!     Input

      integer, intent(in) :: iunit
      character(len=*)    :: fnamelist
      integer, intent(in) :: nnzp
      real,    intent(in) :: gzt(nnzp), gzm(nnzp)
      integer, intent(in) :: day, month, year
      real, intent(in)    :: time
      real, intent(in)    :: delt

!     Namelist

      integer, parameter :: nvarmax = 250
      character(len=20) vars_zt  (nvarmax)
      character(len=20) vars_zm  (nvarmax)
      character(len=20) vars_sfc (nvarmax)

      namelist /statsnl/
     .  lstats, stats_tsamp, stats_tout,
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
        lstats_samp = .false.
        lstats_first = .false.
        lstats_last = .false.
        return
      endif

!     Check sampling and output frequencies

      if ( abs( stats_tsamp/delt - int(stats_tsamp/delt) ) 
     .     .gt. 1.e-8 ) then
         lerror = .true.
         write(*,*) 'Error: stats_tsamp should be a multiple of delt'
         write(*,*) 'stats_tsamp = ',stats_tsamp
         write(*,*) 'delt = ',delt
      endif

      if ( abs( stats_tout/stats_tsamp - int(stats_tout/stats_tsamp) ) 
     .     .gt. 1.e-8 ) then
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
      if ( ntot .eq. nvarmax )
     .  write(*,*) 'WARNING: check nvarmax in statistics.f'

      zt%nn = ntot
      zt%kk = nnzp
!      write(*,*) 'Number of variables for zt ',zt%nn

      allocate( zt%z( zt%kk ) )
      zt%z = gzt

      allocate( zt%x( zt%kk, zt%nn ) )
      allocate( zt%n( zt%kk, zt%nn ) )
      call stats_zero( zt )

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

!     Open GrADS files

      fdir = "./"
      fname = trim(fname_zt)
      call open_grads( zt%f, 20, fdir, fname, 
     .                 1, zt%kk, zt%z,
!     .                 day, month, year, time+stats_tout, stats_tout,
     .                 day, month, year, time, stats_tout,   ! Brian
     .                 zt%nn )

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
      irsm          = 0  ! Brian
      irrm          = 0  ! Brian
      iNrm          = 0  ! Brian
      icond_rrm     = 0  ! Brian
      iauto_rrm     = 0  ! Brian
      iaccr_rrm     = 0  ! Brian
      icond_Nrm     = 0  ! Brian
      iauto_Nrm     = 0  ! Brian
      ised_rrm      = 0  ! Brian
      idiff_rrm     = 0  ! Brian
      isubs_rrm     = 0  ! Brian
      iVr           = 0  ! Brian
      imean_vol_rad = 0  ! Brian
      irain_rate    = 0  ! Brian
      iF_precip     = 0  ! Brian
      iAKm          = 0  ! analytic Kessler.  Vince Larson 22 May 2005
      iAKm_est      = 0  ! LH Kessler.  Vince Larson 22 May 2005

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
      ir1           = 0
      ir2           = 0
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

!     Assign pointers for statistics variables zt

      k = 1
      do i=1,zt%nn

        select case ( trim(vars_zt(i)) )
        case ('thlm')
          ithlm = k
          zt%f%var(ithlm)%ptr => zt%x(:,k)
          zt%f%var(ithlm)%name = "thlm"
          zt%f%var(ithlm)%description = "thetal (K)"
          k = k + 1

        case ('thvm')
          ithvm = k
          zt%f%var(ithvm)%ptr => zt%x(:,k)
          zt%f%var(ithvm)%name = "thvm"
          zt%f%var(ithvm)%description 
     .    = "virtual potential temperature (K)"
          k = k + 1

        case ('rtm')
          irtm = k
          zt%f%var(irtm)%ptr => zt%x(:,k)
          zt%f%var(irtm)%name = "rtm"
          zt%f%var(irtm)%description 
     .    = "total water mixing ratio (kg/kg)"
          k = k + 1

        case ('rcm')
          ircm = k
          zt%f%var(ircm)%ptr => zt%x(:,k)
          zt%f%var(ircm)%name = "rcm"
          zt%f%var(ircm)%description 
     .    = "liquid water mixing ratio (kg/kg)"
          k = k + 1

        case ('um')
          ium = k
          zt%f%var(ium)%ptr => zt%x(:,k)
          zt%f%var(ium)%name = "um"
          zt%f%var(ium)%description = "u wind (m/s)"
          k = k + 1

        case ('vm')
          ivm = k
          zt%f%var(ivm)%ptr => zt%x(:,k)
          zt%f%var(ivm)%name = "vm"
          zt%f%var(ivm)%description = "v wind (m/s)"
          k = k + 1

        case ('wmt')
          iwmt = k
          zt%f%var(iwmt)%ptr => zt%x(:,k)
          zt%f%var(iwmt)%name = "wm"
          zt%f%var(iwmt)%description = "w wind (m/s)"
          k = k + 1

        case ('ug')
          iug = k
          zt%f%var(iug)%ptr => zt%x(:,k)
          zt%f%var(iug)%name = "ug"
          zt%f%var(iug)%description = "u geostrophic wind (m/s)"
          k = k + 1

        case ('vg')
          ivg = k
          zt%f%var(ivg)%ptr => zt%x(:,k)
          zt%f%var(ivg)%name = "vg"
          zt%f%var(ivg)%description = "v geostrophic wind (m/s)"
          k = k + 1

        case ('cf')
          icf = k
          zt%f%var(icf)%ptr => zt%x(:,k)
          zt%f%var(icf)%name = "cf"
          zt%f%var(icf)%description = "cloud fraction"
          k = k + 1

        case ('p')
          ip = k
          zt%f%var(ip)%ptr => zt%x(:,k)
          zt%f%var(ip)%name = "p"
          zt%f%var(ip)%description = "pressure (Pa)"
          k = k + 1

        case ('Lscale')
          iLscale = k
          zt%f%var(iLscale)%ptr => zt%x(:,k)
          zt%f%var(iLscale)%name = "Lscale"
          zt%f%var(iLscale)%description = "mixing length"
          k = k + 1

        case ('thlm_forcing')
          ithlm_forcing = k
          zt%f%var(ithlm_forcing)%ptr => zt%x(:,k)
          zt%f%var(ithlm_forcing)%name = "thlm_f"
          zt%f%var(ithlm_forcing)%description = "thetal forcing"
          k = k + 1

        case ('rtm_forcing')
          irtm_forcing = k
          zt%f%var(irtm_forcing)%ptr => zt%x(:,k)
          zt%f%var(irtm_forcing)%name = "rtm_f"
          zt%f%var(irtm_forcing)%description = "rt forcing"
          k = k + 1

        case ('wp3')
          iwp3 = k
          zt%f%var(iwp3)%ptr => zt%x(:,k)
          zt%f%var(iwp3)%name = "wp3"
          zt%f%var(iwp3)%description = "w third order moment"
          k = k + 1

        case ('wpthlp2')
          iwpthlp2 = k
          zt%f%var(iwpthlp2)%ptr => zt%x(:,k)
          zt%f%var(iwpthlp2)%name = "wpthlp2"
          zt%f%var(iwpthlp2)%description = "wpthlp2 covariance"
          k = k + 1

        case ('wp2thlp')
          iwp2thlp = k
          zt%f%var(iwp2thlp)%ptr => zt%x(:,k)
          zt%f%var(iwp2thlp)%name = "wp2thlp"
          zt%f%var(iwp2thlp)%description = "wp2thlp covariance"
          k = k + 1

        case ('wprtp2')
          iwprtp2 = k
          zt%f%var(iwprtp2)%ptr => zt%x(:,k)
          zt%f%var(iwprtp2)%name = "wprtp2"
          zt%f%var(iwprtp2)%description = "wprtp2 covariance"
          k = k + 1

        case ('wp2rtp')
          iwp2rtp = k
          zt%f%var(iwp2rtp)%ptr => zt%x(:,k)
          zt%f%var(iwp2rtp)%name = "wp2rtp"
          zt%f%var(iwp2rtp)%description = "wp2rtp covariance"
          k = k + 1

        case ('lup')
          ilup = k
          zt%f%var(ilup)%ptr => zt%x(:,k)
          zt%f%var(ilup)%name = "lup"
          zt%f%var(ilup)%description = "upward mixing length"
          k = k + 1

        case ('ldown')
          ildown = k
          zt%f%var(ildown)%ptr => zt%x(:,k)
          zt%f%var(ildown)%name = "ldown"
          zt%f%var(ildown)%description = "downward mixing length"
          k = k + 1

        case ('taut')
          itaut = k
          zt%f%var(itaut)%ptr => zt%x(:,k)
          zt%f%var(itaut)%name = "taut"
          zt%f%var(itaut)%description = "dissipation time"
          k = k + 1

        case ('kht')
          ikht = k
          zt%f%var(ikht)%ptr => zt%x(:,k)
          zt%f%var(ikht)%name = "Kht"
          zt%f%var(ikht)%description = "eddy diffusivity"
          k = k + 1

        case ('wp2thvp')
          iwp2thvp = k
          zt%f%var(iwp2thvp)%ptr => zt%x(:,k)
          zt%f%var(iwp2thvp)%name = "wp2thvp"
          zt%f%var(iwp2thvp)%description = "wp2thvp"
          k = k + 1

        case ('wp2rcp')
          iwp2rcp = k
          zt%f%var(iwp2rcp)%ptr => zt%x(:,k)
          zt%f%var(iwp2rcp)%name = "wp2rcp"
          zt%f%var(iwp2rcp)%description = "wp2rcp"
          k = k + 1

        case ('wprtpthlp')
          iwprtpthlp = k
          zt%f%var(iwprtpthlp)%ptr => zt%x(:,k)
          zt%f%var(iwprtpthlp)%name = "wprtpthlp"
          zt%f%var(iwprtpthlp)%description = "wprtpthlp"
          k = k + 1

        case ('sc')
          isct = k
          zt%f%var(isct)%ptr => zt%x(:,k)
          zt%f%var(isct)%name = "Sc"
          zt%f%var(isct)%description = "Sc"
          k = k + 1

        case ('rhot')
          irhot = k
          zt%f%var(irhot)%ptr => zt%x(:,k)
          zt%f%var(irhot)%name = "rhot"
          zt%f%var(irhot)%description = "density"
          k = k + 1

        case ('rsm')           ! Brian
          irsm = k
          zt%f%var(irsm)%ptr => zt%x(:,k)
          zt%f%var(irsm)%name = "rsm"
          zt%f%var(irsm)%description
     .    = "saturation mixing ratio (kg/kg)"
          k = k + 1

        case ('rrm')           ! Brian
          irrm = k
          zt%f%var(irrm)%ptr => zt%x(:,k)
          zt%f%var(irrm)%name = "rrm"
          zt%f%var(irrm)%description
     .    = "rain water mixing ratio (kg/kg)"
          k = k + 1

        case ('Nrm')           ! Brian
          iNrm = k
          zt%f%var(iNrm)%ptr => zt%x(:,k)
          zt%f%var(iNrm)%name = "Nrm"
          zt%f%var(iNrm)%description
     .    = "rain droplet number concentration (num/m^3)"
          k = k + 1

        case ('cond_rrm')      ! Brian
          icond_rrm = k
          zt%f%var(icond_rrm)%ptr => zt%x(:,k)
          zt%f%var(icond_rrm)%name = "cond_rrm"
          zt%f%var(icond_rrm)%description
     .    = "d(rrm)/dt due to condensation ((kg/kg)/s)"
          k = k + 1

        case ('auto_rrm')      ! Brian
          iauto_rrm = k
          zt%f%var(iauto_rrm)%ptr => zt%x(:,k)
          zt%f%var(iauto_rrm)%name = "auto_rrm"
          zt%f%var(iauto_rrm)%description
     .    = "d(rrm)/dt due to autoconversion ((kg/kg)/s)"
          k = k + 1

        case ('accr_rrm')      ! Brian
          iaccr_rrm = k
          zt%f%var(iaccr_rrm)%ptr => zt%x(:,k)
          zt%f%var(iaccr_rrm)%name = "accr_rrm"
          zt%f%var(iaccr_rrm)%description
     .    = "d(rrm)/dt due to accretion ((kg/kg)/s)"
          k = k + 1

        case ('cond_Nrm')      ! Brian
          icond_Nrm = k
          zt%f%var(icond_Nrm)%ptr => zt%x(:,k)
          zt%f%var(icond_Nrm)%name = "cond_Nrm"
          zt%f%var(icond_Nrm)%description
     .    = "d(Nrm)/dt due to condensation ((num/m^3)/s)"
          k = k + 1

        case ('auto_Nrm')      ! Brian
          iauto_Nrm = k
          zt%f%var(iauto_Nrm)%ptr => zt%x(:,k)
          zt%f%var(iauto_Nrm)%name = "auto_Nrm"
          zt%f%var(iauto_Nrm)%description
     .    = "d(Nrm)/dt due to autoconversion ((num/m^3)/s)"
          k = k + 1

        case ('sed_rrm')       ! Brian
          ised_rrm = k
          zt%f%var(ised_rrm)%ptr => zt%x(:,k)
          zt%f%var(ised_rrm)%name = "sed_rrm"
          zt%f%var(ised_rrm)%description
     .    = "d(rrm)/dt due to sedimentation ((kg/kg)/s)"
          k = k + 1

        case ('diff_rrm')      ! Brian
          idiff_rrm = k
          zt%f%var(idiff_rrm)%ptr => zt%x(:,k)
          zt%f%var(idiff_rrm)%name = "diff_rrm"
          zt%f%var(idiff_rrm)%description
     .    = "d(rrm)/dt due to diffusion ((kg/kg)/s)"
          k = k + 1

        case ('subs_rrm')      ! Brian
          isubs_rrm = k
          zt%f%var(isubs_rrm)%ptr => zt%x(:,k)
          zt%f%var(isubs_rrm)%name = "subs_rrm"
          zt%f%var(isubs_rrm)%description
     .    = "d(rrm)/dt due to subsidence ((kg/kg)/s)"
          k = k + 1

        case ('Vr')            ! Brian
          iVr = k
          zt%f%var(iVr)%ptr => zt%x(:,k)
          zt%f%var(iVr)%name = "Vr"
          zt%f%var(iVr)%description
     .    = "sedimentation velocity (m/s)"
          k = k + 1

        case ('mean_vol_rad')  ! Brian
          imean_vol_rad = k
          zt%f%var(imean_vol_rad)%ptr => zt%x(:,k)
          zt%f%var(imean_vol_rad)%name = "mvr"
          zt%f%var(imean_vol_rad)%description
     .    = "rain drop mean volume radius (m)"
          k = k + 1

        case ('rain_rate')     ! Brian
          irain_rate = k
          zt%f%var(irain_rate)%ptr => zt%x(:,k)
          zt%f%var(irain_rate)%name = "rain_rate"
          zt%f%var(irain_rate)%description = "rain rate (mm/day)"
          k = k + 1

        case ('F_precip')      ! Brian
          iF_precip = k
          zt%f%var(iF_precip)%ptr => zt%x(:,k)
          zt%f%var(iF_precip)%name = "F_precip"
          zt%f%var(iF_precip)%description = "precipitation flux (W/m^2)"
          k = k + 1
 
        case ('AKm')           ! Vince Larson 22 May 2005
          iAKm = k
          zt%f%var(iAKm)%ptr => zt%x(:,k)
          zt%f%var(iAKm)%name = "AKm"
          zt%f%var(iAKm)%description = "Analytic Kessler ac [kg/kg]"
          k = k + 1
 
        case ('AKm_est')       ! Vince Larson 22 May 2005
          iAKm_est = k
          zt%f%var(iAKm_est)%ptr => zt%x(:,k)
          zt%f%var(iAKm_est)%name = "AKm_est"
          zt%f%var(iAKm_est)%description = "LH Kessler estimate [kg/kg]"
          k = k + 1
 
        case ('wp3_bt')
          iwp3_bt = k
          zt%f%var(iwp3_bt)%ptr => zt%x(:,k)
          zt%f%var(iwp3_bt)%name = "wp3_bt"
          zt%f%var(iwp3_bt)%description = "wp3 budget"
          k = k + 1
 
        case ('wp3_ma')
          iwp3_ma = k
          zt%f%var(iwp3_ma)%ptr => zt%x(:,k)
          zt%f%var(iwp3_ma)%name = "wp3_ma"
          zt%f%var(iwp3_ma)%description = "wp3 mean advection"
          k = k + 1
 
        case ('wp3_ta')
          iwp3_ta = k
          zt%f%var(iwp3_ta)%ptr => zt%x(:,k)
          zt%f%var(iwp3_ta)%name = "wp3_ta"
          zt%f%var(iwp3_ta)%description = "wp3 turbulent advection"
          k = k + 1
 
        case ('wp3_tp')
          iwp3_tp = k
          zt%f%var(iwp3_tp)%ptr => zt%x(:,k)
          zt%f%var(iwp3_tp)%name = "wp3_tp"
          zt%f%var(iwp3_tp)%description = "wp3 turbulent transport"
          k = k + 1
 
        case ('wp3_ac')
          iwp3_ac = k
          zt%f%var(iwp3_ac)%ptr => zt%x(:,k)
          zt%f%var(iwp3_ac)%name = "wp3_ac"
          zt%f%var(iwp3_ac)%description = "wp3 accumulation term"
          k = k + 1
 
        case ('wp3_bp')
          iwp3_bp = k
          zt%f%var(iwp3_bp)%ptr => zt%x(:,k)
          zt%f%var(iwp3_bp)%name = "wp3_bp"
          zt%f%var(iwp3_bp)%description = "wp3 buoyancy production"
          k = k + 1
 
        case ('wp3_pr1')
          iwp3_pr1 = k
          zt%f%var(iwp3_pr1)%ptr => zt%x(:,k)
          zt%f%var(iwp3_pr1)%name = "wp3_pr1"
          zt%f%var(iwp3_pr1)%description = "wp3 pressure term 1"
          k = k + 1
 
        case ('wp3_pr2')
          iwp3_pr2 = k
          zt%f%var(iwp3_pr2)%ptr => zt%x(:,k)
          zt%f%var(iwp3_pr2)%name = "wp3_pr2"
          zt%f%var(iwp3_pr2)%description = "wp3 pressure term 2"
          k = k + 1
 
        case ('wp3_dp1')
          iwp3_dp1 = k
          zt%f%var(iwp3_dp1)%ptr => zt%x(:,k)
          zt%f%var(iwp3_dp1)%name = "wp3_dp1"
          zt%f%var(iwp3_dp1)%description = "wp3 dissipation term 1"
          k = k + 1
 
        case ('wp3_cl')
          iwp3_cl = k
          zt%f%var(iwp3_cl)%ptr => zt%x(:,k)
          zt%f%var(iwp3_cl)%name = "wp3_cl"
          zt%f%var(iwp3_cl)%description = "wp3 clipping term"
          k = k + 1
 
        case ('a')
          ia = k
          zt%f%var(ia)%ptr => zt%x(:,k)
          zt%f%var(ia)%name = "a"
          zt%f%var(ia)%description = "pdf parameter a"
          k = k + 1
 
        case ('w1')
          iw1 = k
          zt%f%var(iw1)%ptr => zt%x(:,k)
          zt%f%var(iw1)%name = "w1"
          zt%f%var(iw1)%description = "pdf parameter w1"
          k = k + 1
 
        case ('w2')
          iw2 = k
          zt%f%var(iw2)%ptr => zt%x(:,k)
          zt%f%var(iw2)%name = "w2"
          zt%f%var(iw2)%description = "pdf parameter w2"
          k = k + 1
 
        case ('sw1')
          isw1 = k
          zt%f%var(isw1)%ptr => zt%x(:,k)
          zt%f%var(isw1)%name = "sw1"
          zt%f%var(isw1)%description = "pdf parameter sw1"
          k = k + 1
 
        case ('sw2')
          isw2 = k
          zt%f%var(isw2)%ptr => zt%x(:,k)
          zt%f%var(isw2)%name = "sw2"
          zt%f%var(isw2)%description = "pdf parameter sw2"
          k = k + 1
 
        case ('thl1')
          ithl1 = k
          zt%f%var(ithl1)%ptr => zt%x(:,k)
          zt%f%var(ithl1)%name = "thl1"
          zt%f%var(ithl1)%description = "pdf parameter thl1"
          k = k + 1
 
        case ('thl2')
          ithl2 = k
          zt%f%var(ithl2)%ptr => zt%x(:,k)
          zt%f%var(ithl2)%name = "thl2"
          zt%f%var(ithl2)%description = "pdf parameter thl2"
          k = k + 1
 
        case ('sthl1')
          isthl1 = k
          zt%f%var(isthl1)%ptr => zt%x(:,k)
          zt%f%var(isthl1)%name = "sthl1"
          zt%f%var(isthl1)%description = "pdf parameter sthl1"
          k = k + 1
 
        case ('sthl2')
          isthl2 = k
          zt%f%var(isthl2)%ptr => zt%x(:,k)
          zt%f%var(isthl2)%name = "sthl2"
          zt%f%var(isthl2)%description = "pdf parameter sthl2"
          k = k + 1
 
        case ('rt1')
          irt1 = k
          zt%f%var(irt1)%ptr => zt%x(:,k)
          zt%f%var(irt1)%name = "rt1"
          zt%f%var(irt1)%description = "pdf parameter rt1"
          k = k + 1
 
        case ('rt2')
          irt2 = k
          zt%f%var(irt2)%ptr => zt%x(:,k)
          zt%f%var(irt2)%name = "rt2"
          zt%f%var(irt2)%description = "pdf parameter rt2"
          k = k + 1
 
        case ('srt1')
          isrt1 = k
          zt%f%var(isrt1)%ptr => zt%x(:,k)
          zt%f%var(isrt1)%name = "srt1"
          zt%f%var(isrt1)%description = "pdf parameter srt1"
          k = k + 1
 
        case ('srt2')
          isrt2 = k
          zt%f%var(isrt2)%ptr => zt%x(:,k)
          zt%f%var(isrt2)%name = "srt2"
          zt%f%var(isrt2)%description = "pdf parameter srt2"
          k = k + 1
 
        case ('rc1')
          irc1 = k
          zt%f%var(irc1)%ptr => zt%x(:,k)
          zt%f%var(irc1)%name = "rc1"
          zt%f%var(irc1)%description = "pdf parameter rc1"
          k = k + 1
 
        case ('rc2')
          irc2 = k
          zt%f%var(irc2)%ptr => zt%x(:,k)
          zt%f%var(irc2)%name = "rc2"
          zt%f%var(irc2)%description = "pdf parameter rc2"
          k = k + 1
 
        case ('rsl1')
          irsl1 = k
          zt%f%var(irsl1)%ptr => zt%x(:,k)
          zt%f%var(irsl1)%name = "rsl1"
          zt%f%var(irsl1)%description = "pdf parameter rsl1"
          k = k + 1
 
        case ('rsl2')
          irsl2 = k
          zt%f%var(irsl2)%ptr => zt%x(:,k)
          zt%f%var(irsl2)%name = "rsl2"
          zt%f%var(irsl2)%description = "pdf parameter rsl2"
          k = k + 1
 
        case ('r1')
          ir1 = k
          zt%f%var(ir1)%ptr => zt%x(:,k)
          zt%f%var(ir1)%name = "r1"
          zt%f%var(ir1)%description = "pdf parameter r1"
          k = k + 1
 
        case ('r2')
          ir2 = k
          zt%f%var(ir2)%ptr => zt%x(:,k)
          zt%f%var(ir2)%name = "r2"
          zt%f%var(ir2)%description = "pdf parameter r2"
          k = k + 1
 
        case ('s1')
          is1 = k
          zt%f%var(is1)%ptr => zt%x(:,k)
          zt%f%var(is1)%name = "s1"
          zt%f%var(is1)%description = "pdf parameter s1"
          k = k + 1
 
        case ('s2')
          is2 = k
          zt%f%var(is2)%ptr => zt%x(:,k)
          zt%f%var(is2)%name = "s2"
          zt%f%var(is2)%description = "pdf parameter s2"
          k = k + 1
 
        case ('ss1')
          iss1 = k
          zt%f%var(iss1)%ptr => zt%x(:,k)
          zt%f%var(iss1)%name = "ss1"
          zt%f%var(iss1)%description = "pdf parameter ss1"
          k = k + 1
 
        case ('ss2')
          iss2 = k
          zt%f%var(iss2)%ptr => zt%x(:,k)
          zt%f%var(iss2)%name = "ss2"
          zt%f%var(iss2)%description = "pdf parameter ss2"
          k = k + 1
 
        case ('rrtthl')
          irrtthl = k
          zt%f%var(irrtthl)%ptr => zt%x(:,k)
          zt%f%var(irrtthl)%name = "rrtthl"
          zt%f%var(irrtthl)%description = "pdf parameter rrtthl"
          k = k + 1

        case('wp2zt')
          iwp2zt = k
          zt%f%var(iwp2zt)%ptr => zt%x(:,k)
          zt%f%var(iwp2zt)%name = "wp2zt"
          zt%f%var(iwp2zt)%description ="wp2zt"
          k = k + 1

        case('thlp2zt')
          ithlp2zt = k
          zt%f%var(ithlp2zt)%ptr => zt%x(:,k)
          zt%f%var(ithlp2zt)%name = "thlp2zt"
          zt%f%var(ithlp2zt)%description ="thlp2zt"
          k = k + 1

        case('wpthlpzt')
          iwpthlpzt = k
          zt%f%var(iwpthlpzt)%ptr => zt%x(:,k)
          zt%f%var(iwpthlpzt)%name = "wpthlpzt"
          zt%f%var(iwpthlpzt)%description ="wpthlpzt"
          k = k + 1

        case('wprtpzt')   
          iwprtpzt = k
          zt%f%var(iwprtpzt)%ptr => zt%x(:,k)
          zt%f%var(iwprtpzt)%name = "wprtpzt"
          zt%f%var(iwprtpzt)%description = "wprtpzt"
          k = k + 1

        case('rtp2zt')   
          irtp2zt = k
          zt%f%var(irtp2zt)%ptr => zt%x(:,k)
          zt%f%var(irtp2zt)%name = "rtp2zt"
          zt%f%var(irtp2zt)%description = "rtp2zt"
          k = k + 1
 
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
      if ( ntot .eq. nvarmax )
     .  write(*,*) 'WARNING: check nvarmax in statistics.f'

      zm%nn = ntot
      zm%kk = nnzp
!      write(*,*) 'Number of variables for zm ',zm%nn

      allocate( zm%z( zm%kk ) )
      zm%z = gzm

      allocate( zm%x( zm%kk, zm%nn ) )
      allocate( zm%n( zm%kk, zm%nn ) )
      call stats_zero( zm )

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
      call open_grads( zm%f, 20, fdir, fname, 
     .                 1, zm%kk, zm%z,
!     .                 day, month, year, time+stats_tout, stats_tout,
     .                 day, month, year, time, stats_tout,   ! Brian
     .                 zm%nn )

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
      ifrad     = 0
      iscm      = 0
      iem       = 0
      ishear    = 0    ! Brian

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

      iwpqtp_bt  = 0
      iwpqtp_ma  = 0
      iwpqtp_ta  = 0
      iwpqtp_tp  = 0
      iwpqtp_ac  = 0
      iwpqtp_bp  = 0
      iwpqtp_pr1 = 0
      iwpqtp_pr2 = 0
      iwpqtp_pr3 = 0
      iwpqtp_dp1 = 0

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

!     Assign pointers for statistics variables zm

      k = 1
      do i=1,zm%nn

        select case ( trim(vars_zm(i)) )

        case ('wp2')
          iwp2 = k
          zm%f%var(iwp2)%ptr => zm%x(:,k)
          zm%f%var(iwp2)%name = "wp2"
          zm%f%var(iwp2)%description = "wp2"
          k = k + 1

        case ('rtp2')
          irtp2 = k
          zm%f%var(irtp2)%ptr => zm%x(:,k)
          zm%f%var(irtp2)%name = "rtp2"
          zm%f%var(irtp2)%description = "rtp2"
          k = k + 1

        case ('thlp2')
          ithlp2 = k
          zm%f%var(ithlp2)%ptr => zm%x(:,k)
          zm%f%var(ithlp2)%name = "thlp2"
          zm%f%var(ithlp2)%description = "thlp2"
          k = k + 1

        case ('rtpthlp')
          irtpthlp = k
          zm%f%var(irtpthlp)%ptr => zm%x(:,k)
          zm%f%var(irtpthlp)%name = "rtpthlp"
          zm%f%var(irtpthlp)%description = "rtpthlp"
          k = k + 1

        case ('wprtp')
          iwprtp = k
          zm%f%var(iwprtp)%ptr => zm%x(:,k)
          zm%f%var(iwprtp)%name = "wprtp"
          zm%f%var(iwprtp)%description = "wprtp"
          k = k + 1

        case ('wpthlp')
          iwpthlp = k
          zm%f%var(iwpthlp)%ptr => zm%x(:,k)
          zm%f%var(iwpthlp)%name = "wpthlp"
          zm%f%var(iwpthlp)%description = "wpthlp"
          k = k + 1

        case ('wp4')
          iwp4 = k
          zm%f%var(iwp4)%ptr => zm%x(:,k)
          zm%f%var(iwp4)%name = "wp4"
          zm%f%var(iwp4)%description = "wp4"
          k = k + 1

        case ('wpthvp')
          iwpthvp = k
          zm%f%var(iwpthvp)%ptr => zm%x(:,k)
          zm%f%var(iwpthvp)%name = "wpthvp"
          zm%f%var(iwpthvp)%description = "Buoyancy flux (K m/s)"
          k = k + 1

        case ('rtpthvp')
          irtpthvp = k
          zm%f%var(irtpthvp)%ptr => zm%x(:,k)
          zm%f%var(irtpthvp)%name = "rtpthvp"
          zm%f%var(irtpthvp)%description = "rtpthvp"
          k = k + 1

        case ('thlpthvp')
          ithlpthvp = k
          zm%f%var(ithlpthvp)%ptr => zm%x(:,k)
          zm%f%var(ithlpthvp)%name = "thlpthvp"
          zm%f%var(ithlpthvp)%description = "thlpthvp"
          k = k + 1

        case ('taum')
          itaum = k
          zm%f%var(itaum)%ptr => zm%x(:,k)
          zm%f%var(itaum)%name = "taum"
          zm%f%var(itaum)%description = "dissipation time"
          k = k + 1

        case ('khm')
          ikhm = k
          zm%f%var(ikhm)%ptr => zm%x(:,k)
          zm%f%var(ikhm)%name = "khm"
          zm%f%var(ikhm)%description = "eddy diffusivity"
          k = k + 1

        case ('wprcp')
          iwprcp = k
          zm%f%var(iwprcp)%ptr => zm%x(:,k)
          zm%f%var(iwprcp)%name = "wprcp"
          zm%f%var(iwprcp)%description = "wprcp"
          k = k + 1

        case ('thlprcp')
          ithlprcp = k
          zm%f%var(ithlprcp)%ptr => zm%x(:,k)
          zm%f%var(ithlprcp)%name = "thlprcp"
          zm%f%var(ithlprcp)%description = "thlprcp"
          k = k + 1

        case ('rtprcp')
          irtprcp = k
          zm%f%var(irtprcp)%ptr => zm%x(:,k)
          zm%f%var(irtprcp)%name = "rtprcp"
          zm%f%var(irtprcp)%description = "rtprcp"
          k = k + 1

        case ('upwp')
          iupwp = k
          zm%f%var(iupwp)%ptr => zm%x(:,k)
          zm%f%var(iupwp)%name = "upwp"
          zm%f%var(iupwp)%description = "upwp"
          k = k + 1

        case ('vpwp')
          ivpwp = k
          zm%f%var(ivpwp)%ptr => zm%x(:,k)
          zm%f%var(ivpwp)%name = "vpwp"
          zm%f%var(ivpwp)%description = "vpwp"
          k = k + 1

        case ('rhom')
          irhom = k
          zm%f%var(irhom)%ptr => zm%x(:,k)
          zm%f%var(irhom)%name = "rhom"
          zm%f%var(irhom)%description = "density"
          k = k + 1

        case ('frad')
          ifrad = k
          zm%f%var(ifrad)%ptr => zm%x(:,k)
          zm%f%var(ifrad)%name = "frad"
          zm%f%var(ifrad)%description = "radiative flux"
          k = k + 1

        case ('sc')
          iscm = k
          zm%f%var(iscm)%ptr => zm%x(:,k)
          zm%f%var(iscm)%name = "scm"
          zm%f%var(iscm)%description = "scm"
          k = k + 1

        case ('em')
          iem = k
          zm%f%var(iem)%ptr => zm%x(:,k)
          zm%f%var(iem)%name = "em"
          zm%f%var(iem)%description = "em"
          k = k + 1

        case ('shear')      ! Brian
          ishear = k
          zm%f%var(ishear)%ptr => zm%x(:,k)
          zm%f%var(ishear)%name = "shear"
          zm%f%var(ishear)%description = "wind shear term (m^2/s^3)"
          k = k + 1

        case ('wp2_bt')
          iwp2_bt = k
          zm%f%var(iwp2_bt)%ptr => zm%x(:,k)
          zm%f%var(iwp2_bt)%name = "wp2_bt"
          zm%f%var(iwp2_bt)%description = "wp2 budget"
          k = k + 1

        case ('wp2_ma')
          iwp2_ma = k
          zm%f%var(iwp2_ma)%ptr => zm%x(:,k)
          zm%f%var(iwp2_ma)%name = "wp2_ma"
          zm%f%var(iwp2_ma)%description = "wp2 mean advection"
          k = k + 1

        case ('wp2_ta')
          iwp2_ta = k
          zm%f%var(iwp2_ta)%ptr => zm%x(:,k)
          zm%f%var(iwp2_ta)%name = "wp2_ta"
          zm%f%var(iwp2_ta)%description = "wp2 turbulent advection"
          k = k + 1

        case ('wp2_ac')
          iwp2_ac = k
          zm%f%var(iwp2_ac)%ptr => zm%x(:,k)
          zm%f%var(iwp2_ac)%name = "wp2_ac"
          zm%f%var(iwp2_ac)%description = "wp2 accumulation term"
          k = k + 1

        case ('wp2_bp')
          iwp2_bp = k
          zm%f%var(iwp2_bp)%ptr => zm%x(:,k)
          zm%f%var(iwp2_bp)%name = "wp2_bp"
          zm%f%var(iwp2_bp)%description = "wp2 buoyancy production"
          k = k + 1

        case ('wp2_pr1')
          iwp2_pr1 = k
          zm%f%var(iwp2_pr1)%ptr => zm%x(:,k)
          zm%f%var(iwp2_pr1)%name = "wp2_pr1"
          zm%f%var(iwp2_pr1)%description = "wp2 pressure term 1"
          k = k + 1

        case ('wp2_pr2')
          iwp2_pr2 = k
          zm%f%var(iwp2_pr2)%ptr => zm%x(:,k)
          zm%f%var(iwp2_pr2)%name = "wp2_pr2"
          zm%f%var(iwp2_pr2)%description = "wp2 pressure term 2"
          k = k + 1

        case ('wp2_pr3')
          iwp2_pr3 = k
          zm%f%var(iwp2_pr3)%ptr => zm%x(:,k)
          zm%f%var(iwp2_pr3)%name = "wp2_pr3"
          zm%f%var(iwp2_pr3)%description = "wp2 pressure term 3"
          k = k + 1

        case ('wp2_dp1')
          iwp2_dp1 = k
          zm%f%var(iwp2_dp1)%ptr => zm%x(:,k)
          zm%f%var(iwp2_dp1)%name = "wp2_dp1"
          zm%f%var(iwp2_dp1)%description = "wp2 dissipation term 1"
          k = k + 1

        case ('wp2_dp2')
          iwp2_dp2 = k
          zm%f%var(iwp2_dp2)%ptr => zm%x(:,k)
          zm%f%var(iwp2_dp2)%name = "wp2_dp2"
          zm%f%var(iwp2_dp2)%description = "wp2 dissipation term 2"
          k = k + 1

        case ('wp2_cl')
          iwp2_cl = k
          zm%f%var(iwp2_cl)%ptr => zm%x(:,k)
          zm%f%var(iwp2_cl)%name = "wp2_cl"
          zm%f%var(iwp2_cl)%description = "wp2 clipping term"
          k = k + 1

        case ('wpqtp_bt')
          iwpqtp_bt = k
          zm%f%var(iwpqtp_bt)%ptr => zm%x(:,k)
          zm%f%var(iwpqtp_bt)%name = "wpqtp_bt"
          zm%f%var(iwpqtp_bt)%description = "wpqtp budget"
          k = k + 1

        case ('wpqtp_ma')
          iwpqtp_ma = k
          zm%f%var(iwpqtp_ma)%ptr => zm%x(:,k)
          zm%f%var(iwpqtp_ma)%name = "wpqtp_ma"
          zm%f%var(iwpqtp_ma)%description = "wpqtp mean advection"
          k = k + 1

        case ('wpqtp_ta')
          iwpqtp_ta = k
          zm%f%var(iwpqtp_ta)%ptr => zm%x(:,k)
          zm%f%var(iwpqtp_ta)%name = "wpqtp_ta"
          zm%f%var(iwpqtp_ta)%description = "wpqtp turbulent advection"
          k = k + 1

        case ('wpqtp_tp')
          iwpqtp_tp = k
          zm%f%var(iwpqtp_tp)%ptr => zm%x(:,k)
          zm%f%var(iwpqtp_tp)%name = "wpqtp_tp"
          zm%f%var(iwpqtp_tp)%description = "wpqtp turbulent production"
          k = k + 1

        case ('wpqtp_ac')
          iwpqtp_ac = k
          zm%f%var(iwpqtp_ac)%ptr => zm%x(:,k)
          zm%f%var(iwpqtp_ac)%name = "wpqtp_ac"
          zm%f%var(iwpqtp_ac)%description = "wpqtp accumulation term"
          k = k + 1

        case ('wpqtp_bp')
          iwpqtp_bp = k
          zm%f%var(iwpqtp_bp)%ptr => zm%x(:,k)
          zm%f%var(iwpqtp_bp)%name = "wpqtp_bp"
          zm%f%var(iwpqtp_bp)%description = "wpqtp buoyancy production"
          k = k + 1

        case ('wpqtp_pr1')
          iwpqtp_pr1 = k
          zm%f%var(iwpqtp_pr1)%ptr => zm%x(:,k)
          zm%f%var(iwpqtp_pr1)%name = "wpqtp_pr1"
          zm%f%var(iwpqtp_pr1)%description = "wpqtp pressure term 1"
          k = k + 1

        case ('wpqtp_pr2')
          iwpqtp_pr2 = k
          zm%f%var(iwpqtp_pr2)%ptr => zm%x(:,k)
          zm%f%var(iwpqtp_pr2)%name = "wpqtp_pr2"
          zm%f%var(iwpqtp_pr2)%description = "wpqtp pressure term 2"
          k = k + 1

        case ('wpqtp_pr3')
          iwpqtp_pr3 = k
          zm%f%var(iwpqtp_pr3)%ptr => zm%x(:,k)
          zm%f%var(iwpqtp_pr3)%name = "wpqtp_pr3"
          zm%f%var(iwpqtp_pr3)%description = "wpqtp pressure term 3"
          k = k + 1

        case ('wpqtp_dp1')
          iwpqtp_dp1 = k
          zm%f%var(iwpqtp_dp1)%ptr => zm%x(:,k)
          zm%f%var(iwpqtp_dp1)%name = "wpqtp_dp1"
          zm%f%var(iwpqtp_dp1)%description = "wpqtp dissipation term 1"
          k = k + 1

        case ('wpthlp_bt')
          iwpthlp_bt = k
          zm%f%var(iwpthlp_bt)%ptr => zm%x(:,k)
          zm%f%var(iwpthlp_bt)%name = "wpthlp_bt"
          zm%f%var(iwpthlp_bt)%description = "wpthlp budget"
          k = k + 1

        case ('wpthlp_ma')
          iwpthlp_ma = k
          zm%f%var(iwpthlp_ma)%ptr => zm%x(:,k)
          zm%f%var(iwpthlp_ma)%name = "wpthlp_ma"
          zm%f%var(iwpthlp_ma)%description = "wpthlp mean advection"
          k = k + 1

        case ('wpthlp_ta')
          iwpthlp_ta = k
          zm%f%var(iwpthlp_ta)%ptr => zm%x(:,k)
          zm%f%var(iwpthlp_ta)%name = "wpthlp_ta"
          zm%f%var(iwpthlp_ta)%description 
     .    = "wpthlp turbulent advection"
          k = k + 1

        case ('wpthlp_tp')
          iwpthlp_tp = k
          zm%f%var(iwpthlp_tp)%ptr => zm%x(:,k)
          zm%f%var(iwpthlp_tp)%name = "wpthlp_tp"
          zm%f%var(iwpthlp_tp)%description
     .    = "wpthlp turbulent production"
          k = k + 1

        case ('wpthlp_ac')
          iwpthlp_ac = k
          zm%f%var(iwpthlp_ac)%ptr => zm%x(:,k)
          zm%f%var(iwpthlp_ac)%name = "wpthlp_ac"
          zm%f%var(iwpthlp_ac)%description
     .    = "wpthlp accumulation term"
          k = k + 1

        case ('wpthlp_bp')
          iwpthlp_bp = k
          zm%f%var(iwpthlp_bp)%ptr => zm%x(:,k)
          zm%f%var(iwpthlp_bp)%name = "wpthlp_bp"
          zm%f%var(iwpthlp_bp)%description
     .    = "wpthlp buoyancy production"
          k = k + 1

        case ('wpthlp_pr1')
          iwpthlp_pr1 = k
          zm%f%var(iwpthlp_pr1)%ptr => zm%x(:,k)
          zm%f%var(iwpthlp_pr1)%name = "wpthlp_pr1"
          zm%f%var(iwpthlp_pr1)%description
     .    = "wpthlp pressure term 1"
          k = k + 1

        case ('wpthlp_pr2')
          iwpthlp_pr2 = k
          zm%f%var(iwpthlp_pr2)%ptr => zm%x(:,k)
          zm%f%var(iwpthlp_pr2)%name = "wpthlp_pr2"
          zm%f%var(iwpthlp_pr2)%description
     .    = "wpthlp pressure term 2"
          k = k + 1

        case ('wpthlp_pr3')
          iwpthlp_pr3 = k
          zm%f%var(iwpthlp_pr3)%ptr => zm%x(:,k)
          zm%f%var(iwpthlp_pr3)%name = "wpthlp_pr3"
          zm%f%var(iwpthlp_pr3)%description
     .    = "wpthlp pressure term 3"
          k = k + 1

        case ('wpthlp_dp1')
          iwpthlp_dp1 = k
          zm%f%var(iwpthlp_dp1)%ptr => zm%x(:,k)
          zm%f%var(iwpthlp_dp1)%name = "wpthlp_dp1"
          zm%f%var(iwpthlp_dp1)%description
     .    = "wpthlp dissipation term 1"
          k = k + 1

        case default
          write(*,*) 'Error: unrecognized variable ',trim(vars_zm(i))
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
      if ( ntot .eq. nvarmax )
     .  write(*,*) 'WARNING: check nvarmax in statistics.f'

      sfc%nn = ntot
      sfc%kk = 1
!      write(*,*) 'Number of variables for sfc ',sfc%nn

      allocate( sfc%z( sfc%kk ) )
      sfc%z = gzm(1)

      allocate( sfc%x( sfc%kk, sfc%nn ) )
      allocate( sfc%n( sfc%kk, sfc%nn ) )
      call stats_zero( sfc )

      allocate( sfc%f%var( sfc%nn ) )
      allocate( sfc%f%z( sfc%kk ) )

!     Open GrADS files

      fdir = "./"
      fname = trim(fname_sfc)
      call open_grads( sfc%f, 20, fdir, fname, 
     .                 1, sfc%kk, sfc%z,
!     .                 day, month, year, time+stats_tout, stats_tout,
     .                 day, month, year, time, stats_tout,   ! Brian
     .                 sfc%nn )

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
          k = k + 1

        case ('lh')
          ilh = k
          sfc%f%var(ilh)%ptr => sfc%x(:,k)
          sfc%f%var(ilh)%name = "lh"
          sfc%f%var(ilh)%description = "lh"
          k = k + 1

        case ('sh')
          ish = k
          sfc%f%var(ish)%ptr => sfc%x(:,k)
          sfc%f%var(ish)%name = "sh"
          sfc%f%var(ish)%description = "sh"
          k = k + 1

        case ('cc')
          icc = k
          sfc%f%var(icc)%ptr => sfc%x(:,k)
          sfc%f%var(icc)%name = "cc"
          sfc%f%var(icc)%description = "cc"
          k = k + 1

        case ('lwp')
          ilwp = k
          sfc%f%var(ilwp)%ptr => sfc%x(:,k)
          sfc%f%var(ilwp)%name = "lwp"
          sfc%f%var(ilwp)%description = "lwp"
          k = k + 1

        case ('zb')
          izb = k
          sfc%f%var(izb)%ptr => sfc%x(:,k)
          sfc%f%var(izb)%name = "zb"
          sfc%f%var(izb)%description = "zb"
          k = k + 1

        case ('zi')
          izi = k
          sfc%f%var(izi)%ptr => sfc%x(:,k)
          sfc%f%var(izi)%name = "zi"
          sfc%f%var(izi)%description = "zi"
          k = k + 1

        case ('rain')          ! Brian
          irain = k
          sfc%f%var(irain)%ptr => sfc%x(:,k)
          sfc%f%var(irain)%name = "rain_rate"
          sfc%f%var(irain)%description = "rainfall rate (mm/day)"
          k = k + 1

        case ('pflux')         ! Brian
          ipflux = k
          sfc%f%var(ipflux)%ptr => sfc%x(:,k)
          sfc%f%var(ipflux)%name = "prec_flux"
          sfc%f%var(ipflux)%description = "Precipitation flux (W/m^2)"
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
      subroutine stats_zero( s )

      use stats_type
      implicit none

!     Input/output

      type (stats), intent(inout) :: s

!     Zero out arrays

      if ( s%nn > 0 ) then
        s%x = 0.0
        s%n = 0
      endif

      return
      end subroutine stats_zero

c     -----------------------------------------------------------------
      subroutine stats_avg( s )

      use stats_type
      implicit none

!     Input/output

      type (stats), intent(inout) :: s

!     Internal

      integer k,n

!     Compute averages

      do n=1,s%nn
        do k=1,s%kk

          if ( s%n(k,n) .gt. 0 ) then
            s%x(k,n) = s%x(k,n) / s%n(k,n)
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

!     Set sample this time step flag
      if ( mod( time, stats_tsamp ) .lt. 1.e-8 ) then
         lstats_samp = .true.
      else
         lstats_samp = .false.
      endif

!     Set first time step flag

      if ( mod( time-delt, stats_tout ) .lt. 1.e-8 ) then
         lstats_first = .true.
      else
         lstats_first = .false.
      endif

!     Set last time step flag

      if ( mod( time, stats_tout ) .lt. 1.e-8 ) then
         lstats_last = .true.
      else
         lstats_last = .false.
      endif

      return

      end subroutine stats_begin_timestep

c     -----------------------------------------------------------------
!      subroutine stats_end_timestep( time, delt )
      subroutine stats_end_timestep( time, delt, deltclose )   ! Brian

      use statistics

      implicit none

!     Input

      real, intent(in) :: time, delt, deltclose   ! Brian

!     Internal

      integer i, k

!     Check if it is time to write to file

      if ( .not. lstats_last ) return

!     Check number of sampling points

      do i=1,zt%nn
       do k=1,zt%kk
!         if ( zt%n(k,i).ne.0 
!     .        .and. zt%n(k,i).ne.stats_tout/stats_tsamp ) then
!         if ( zt%n(k,i).ne.0 
!     .        .and. zt%n(k,i).ne.1 ) then
         if ( zt%n(k,i).ne.0
     .        .and. zt%n(k,i).ne.1
     .        .and. zt%n(k,i) .ne. int(delt/deltclose) ) then  ! Brian
           write(*,*) 'Possible sampling error for zt ',i,k,zt%n(k,i)
         endif
        enddo
      enddo
      
      do i=1,zm%nn
       do k=1,zm%kk
!         if ( zm%n(k,i).ne.0 
!     .        .and. zm%n(k,i).ne.stats_tout/stats_tsamp ) then
!         if ( zm%n(k,i).ne.0 
!     .        .and. zm%n(k,i).ne.1 ) then
         if ( zm%n(k,i).ne.0
     .        .and. zm%n(k,i).ne.1
     .        .and. zm%n(k,i) .ne. int(delt/deltclose) ) then  ! Brian
           write(*,*) 'Possible sampling error for zm ',i,k,zm%n(k,i)
         endif
        enddo
      enddo
      
!     Compute averages

      call stats_avg( zt )
      call stats_avg( zm )
      call stats_avg( sfc )

!     Write to file

      call write_grads( zt%f  )
      call write_grads( zm%f  )
      call write_grads( sfc%f  )

!     Reset sample fields

      call stats_zero( zt )
      call stats_zero( zm )
      call stats_zero( sfc )

      return
      end subroutine stats_end_timestep

c     -----------------------------------------------------------------
      subroutine stats_accumulate( time, delt )

      use statistics
      use grid_class
      use model_arrays

      implicit none

!     Input

      real, intent(in) :: time, delt

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
        if ( irsm > 0 ) then           ! Brian
          zt%x(:,irsm) = zt%x(:,irsm) + rsm
          zt%n(:,irsm) = zt%n(:,irsm) + 1
        endif
        if ( irrm > 0 ) then           ! Brian
          zt%x(:,irrm) = zt%x(:,irrm) + rrm
          zt%n(:,irrm) = zt%n(:,irrm) + 1
        endif
        if ( iNrm > 0 ) then           ! Brian
          zt%x(:,iNrm) = zt%x(:,iNrm) + Nrm
          zt%n(:,iNrm) = zt%n(:,iNrm) + 1
        endif
        if ( icond_rrm > 0 ) then      ! Brian
          zt%x(:,icond_rrm) = zt%x(:,icond_rrm) + cond_rrm
          zt%n(:,icond_rrm) = zt%n(:,icond_rrm) + 1
        endif
        if ( iauto_rrm > 0 ) then      ! Brian
          zt%x(:,iauto_rrm) = zt%x(:,iauto_rrm) + auto_rrm
          zt%n(:,iauto_rrm) = zt%n(:,iauto_rrm) + 1
        endif
        if ( iaccr_rrm > 0 ) then      ! Brian
          zt%x(:,iaccr_rrm) = zt%x(:,iaccr_rrm) + accr_rrm
          zt%n(:,iaccr_rrm) = zt%n(:,iaccr_rrm) + 1
        endif
        if ( icond_Nrm > 0 ) then      ! Brian
          zt%x(:,icond_Nrm) = zt%x(:,icond_Nrm) + cond_Nrm
          zt%n(:,icond_Nrm) = zt%n(:,icond_Nrm) + 1
        endif
        if ( iauto_Nrm > 0 ) then      ! Brian
          zt%x(:,iauto_Nrm) = zt%x(:,iauto_Nrm) + auto_Nrm
          zt%n(:,iauto_Nrm) = zt%n(:,iauto_Nrm) + 1
        endif
        if ( ised_rrm > 0 ) then       ! Brian
          zt%x(:,ised_rrm) = zt%x(:,ised_rrm) + sed_rrm
          zt%n(:,ised_rrm) = zt%n(:,ised_rrm) + 1
        endif
        if ( idiff_rrm > 0 ) then      ! Brian
          zt%x(:,idiff_rrm) = zt%x(:,idiff_rrm) + diff_rrm
          zt%n(:,idiff_rrm) = zt%n(:,idiff_rrm) + 1
        endif
        if ( isubs_rrm > 0 ) then      ! Brian
          zt%x(:,isubs_rrm) = zt%x(:,isubs_rrm) + subs_rrm
          zt%n(:,isubs_rrm) = zt%n(:,isubs_rrm) + 1
        endif
        if ( iVr > 0 ) then            ! Brian
          zt%x(:,iVr) = zt%x(:,iVr) + Vr
          zt%n(:,iVr) = zt%n(:,iVr) + 1
        endif
        if ( imean_vol_rad > 0 ) then  ! Brian
          zt%x(:,imean_vol_rad) = zt%x(:,imean_vol_rad) + mean_vol_rad
          zt%n(:,imean_vol_rad) = zt%n(:,imean_vol_rad) + 1
        endif
        if ( irain_rate > 0 ) then     ! Brian
          zt%x(:,irain_rate) = zt%x(:,irain_rate) + rain_rate
          zt%n(:,irain_rate) = zt%n(:,irain_rate) + 1
        endif
        if ( iF_precip > 0 ) then      ! Brian
          zt%x(:,iF_precip) = zt%x(:,iF_precip) + F_precip
          zt%n(:,iF_precip) = zt%n(:,iF_precip) + 1
        endif

        if ( iAKm > 0 ) then           ! Vince Larson 22 May 2005
          zt%x(:,iAKm) = zt%x(:,iAKm) + AKm
          zt%n(:,iAKm) = zt%n(:,iAKm) + 1
        endif

        if ( iAKm_est > 0 ) then           ! Vince Larson 22 May 2005
          zt%x(:,iAKm_est) = zt%x(:,iAKm_est) + AKm_est
          zt%n(:,iAKm_est) = zt%n(:,iAKm_est) + 1
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

        if ( ir1 > 0 ) then
          zt%x(:,ir1) = zt%x(:,ir1) + pdf_parms(:,18)
          zt%n(:,ir1) = zt%n(:,ir1) + 1
        endif

        if ( ir2 > 0 ) then
          zt%x(:,ir2) = zt%x(:,ir2) + pdf_parms(:,19)
          zt%n(:,ir2) = zt%n(:,ir2) + 1
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
        if ( ifrad > 0 ) then
          zm%x(:,ifrad) = zm%x(:,ifrad) + frad
          zm%n(:,ifrad) = zm%n(:,ifrad) + 1
        endif
        if ( iscm > 0 ) then
          zm%x(:,iscm) = zm%x(:,iscm) + scm
          zm%n(:,iscm) = zm%n(:,iscm) + 1
        endif
        if ( iem > 0 ) then
          zm%x(:,iem) = zm%x(:,iem) + em
          zm%n(:,iem) = zm%n(:,iem) + 1
        endif
        if ( ishear > 0 ) then   ! Brian
          zm%x(:,ishear) = zm%x(:,ishear) + shear
          zm%n(:,ishear) = zm%n(:,ishear) + 1
        endif

! sfc variables

!       Cloud cover
        if ( icc > 0 ) then
          sfc%x(1,icc) = sfc%x(1,icc) + maxval( cf(1:gr%nnzp) )
          sfc%n(1,icc) = sfc%n(1,icc) + 1
        endif

!       Cloud base
        if ( izb > 0 ) then
          k = 1
!         do while ( rcm(k) < 1.e-6 .and. k <= gr%nnzp ) ! array out of bounds
!							  dschanen 24 May 2005
          do while ( rcm(k) < 1.e-6 .and. k < gr%nnzp )
            k = k + 1
          end do
          if ( k < gr%nnzp) then
           sfc%x(1,izb) = sfc%x(1,izb) + gr%zt(k)
          else
           sfc%x(1,izb) = sfc%x(1,izb) - 10
          endif
          sfc%n(1,izb) = sfc%n(1,izb) + 1
        endif

!       LWP
        if ( ilwp > 0 ) then
          xtmp = 0.
          do i=gr%nnzp-1,2,-1
           xtmp = xtmp + rhot(i+1) * rcm(i+1) / gr%dzt(i+1)
          enddo
          sfc%x(1,ilwp) = sfc%x(1,ilwp) + xtmp
          sfc%n(1,ilwp) = sfc%n(1,ilwp) + 1
        end if

      endif

      return
      end subroutine stats_accumulate

#endif STATS
