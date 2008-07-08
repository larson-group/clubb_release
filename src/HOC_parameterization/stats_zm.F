!-----------------------------------------------------------------------
! $Id: stats_zm.F,v 1.15 2008-07-08 04:44:25 dschanen Exp $
      module stats_zm
#ifdef STATS
       
      implicit none
      
      private ! Default Scope
      
      public :: stats_init_zm
      
      contains

!-----------------------------------------------------------------------
      subroutine stats_init_zm( vars_zm, lerror )

!     Description: 
!     Initializes array indices for zm
!-----------------------------------------------------------------------

      use stats_variables, only:
     .    zm,
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
     .    irhom,
     .    iscm,
     .    iem,
     .    ishear,
     .    iFrad,
     .    iFrad_LW,
     .    iFrad_SW,
     .    iFprec,
     .    iFcsed,
     .    iup2_plus_vp2,
     .    iup2,
     .    ivp2,
     .    iup2_bt,
     .    iup2_ta,
     .    iup2_tp,
     .    iup2_ma,
     .    iup2_dp1,
     .    iup2_dp2,
     .    iup2_pr1,
     .    iup2_pr2,
     .    iup2_cl,
     .    iup2_pd,
     .    ivp2_bt,
     .    ivp2_ta,
     .    ivp2_tp,
     .    ivp2_ma,
     .    ivp2_dp1,
     .    ivp2_dp2,
     .    ivp2_pr1,
     .    ivp2_pr2,
     .    ivp2_cl,
     .    ivp2_pd,
     .    iVrr,
     .    iVNr,
     .    iVice,
     .    iVgraupel,
     .    iVsnow,
     .    iwp2_bt,
     .    iwp2_ma,
     .    iwp2_ta,
     .    iwp2_ac,
     .    iwp2_bp,
     .    iwp2_pr1,
     .    iwp2_pr2,
     .    iwp2_pr3,
     .    iwp2_dp1,
     .    iwp2_dp2,
     .    iwp2_cl,
     .    iwp2_pd,
     .    iwprtp_bt,
     .    iwprtp_ma,
     .    iwprtp_ta,
     .    iwprtp_tp,
     .    iwprtp_ac,
     .    iwprtp_bp,
     .    iwprtp_pr1,
     .    iwprtp_pr2,
     .    iwprtp_pr3,
     .    iwprtp_dp1,
     .    iwprtp_cl,
     .    iwprtp_sicl,
     .    iwprtp_pd,
     .    iwpthlp_bt,
     .    iwpthlp_ma,
     .    iwpthlp_ta,
     .    iwpthlp_tp,
     .    iwpthlp_ac,
     .    iwpthlp_bp,
     .    iwpthlp_pr1,
     .    iwpthlp_pr2,
     .    iwpthlp_pr3,
     .    iwpthlp_dp1,
     .    iwpthlp_cl,
     .    iwpthlp_sicl

      use stats_variables, only:
     .    irtp2_bt,
     .    irtp2_ma,
     .    irtp2_ta,
     .    irtp2_tp,
     .    irtp2_dp1,
     .    irtp2_dp2,
     .    irtp2_cl,
     .    irtp2_pd,
     .    ithlp2_bt,
     .    ithlp2_ma,
     .    ithlp2_ta,
     .    ithlp2_tp,
     .    ithlp2_dp1,
     .    ithlp2_dp2,
     .    ithlp2_cl,
     .    ithlp2_pd,
     .    irtpthlp_bt,
     .    irtpthlp_ma,
     .    irtpthlp_ta,
     .    irtpthlp_tp1,
     .    irtpthlp_tp2,
     .    irtpthlp_dp1,
     .    irtpthlp_dp2,
     .    irtpthlp_dp2,
     .    irtpthlp_cl,
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

      use stats_type, only:
     .    stat_assign ! Procedure

      implicit none

      integer, parameter :: nvarmax = 250  ! Max variables

      ! Input Variable
      ! zm variable names
      
      character(len= * ), dimension(nvarmax), intent(in) :: vars_zm

      ! Output Variable	
      logical, intent(inout) :: lerror
      
      ! Local Varables
      integer :: i, k
      
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
      iKhm      = 0
      iwprcp    = 0
      ithlprcp  = 0
      irtprcp   = 0
      ircp2     = 0
      iupwp     = 0
      ivpwp     = 0
      irhom     = 0
      iscm      = 0
      iem       = 0
      ishear    = 0  ! Brian
      iFrad     = 0
      iFrad_LW  = 0  ! Brian
      iFrad_SW  = 0  ! Brian
      iFprec    = 0  ! Brian
      iFcsed    = 0  ! Brian

      iup2_plus_vp2 = 0 ! Michael

      iup2 = 0
      ivp2 = 0

      iup2_bt  = 0
      iup2_ta  = 0
      iup2_tp  = 0
      iup2_ma  = 0
      iup2_dp1 = 0
      iup2_dp2 = 0
      iup2_pr1 = 0
      iup2_pr2 = 0
      iup2_cl  = 0

      ivp2_bt  = 0
      ivp2_ta  = 0
      ivp2_tp  = 0
      ivp2_ma  = 0
      ivp2_dp1 = 0
      ivp2_dp2 = 0
      ivp2_pr1 = 0
      ivp2_pr2 = 0
      ivp2_cl  = 0

      ! Sedimentation velocities
      iVrr      = 0  ! Brian
      iVNr      = 0  !  " "
      iVice     = 0  ! COAMPS
      iVgraupel = 0  !  " "
      iVsnow    = 0  !  " "

      ! Vertical velocity budgets
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
      iwp2_pd   = 0

      ! Flux budgets
      iwprtp_bt   = 0
      iwprtp_ma   = 0
      iwprtp_ta   = 0
      iwprtp_tp   = 0
      iwprtp_ac   = 0
      iwprtp_bp   = 0
      iwprtp_pr1  = 0
      iwprtp_pr2  = 0
      iwprtp_pr3  = 0
      iwprtp_dp1  = 0
      iwprtp_cl   = 0
      iwprtp_sicl = 0
      iwprtp_pd   = 0

      iwpthlp_bt   = 0
      iwpthlp_ma   = 0
      iwpthlp_ta   = 0
      iwpthlp_tp   = 0
      iwpthlp_ac   = 0
      iwpthlp_bp   = 0
      iwpthlp_pr1  = 0
      iwpthlp_pr2  = 0
      iwpthlp_pr3  = 0
      iwpthlp_dp1  = 0
      iwpthlp_cl   = 0
      iwpthlp_sicl = 0

      ! Variance budgets
      irtp2_bt    = 0
      irtp2_ma    = 0
      irtp2_ta    = 0
      irtp2_tp    = 0
      irtp2_dp1   = 0
      irtp2_dp2   = 0
      irtp2_cl    = 0
      irtp2_pd    = 0

      ithlp2_bt    = 0
      ithlp2_ma    = 0
      ithlp2_ta    = 0
      ithlp2_tp    = 0
      ithlp2_dp1   = 0
      ithlp2_dp2   = 0
      ithlp2_cl    = 0
      ithlp2_pd    = 0

      irtpthlp_bt  = 0
      irtpthlp_ma  = 0
      irtpthlp_ta  = 0
      irtpthlp_tp1 = 0
      irtpthlp_tp2 = 0
      irtpthlp_dp1 = 0
      irtpthlp_dp2 = 0
      irtpthlp_cl  = 0

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

!     Assign pointers for statistics variables zm

      k = 1
      do i=1,zm%nn

        select case ( trim(vars_zm(i)) )

        case ('wp2')
          iwp2 = k
          call stat_assign(iwp2,"wp2",
     .         "wp2","m^2/s^2",zm)
          k = k + 1

        case ('rtp2')
          irtp2 = k
          call stat_assign(irtp2,"rtp2",
     .         "rtp2","(kg/kg)^2",zm)
          k = k + 1

        case ('thlp2')
          ithlp2 = k
          call stat_assign(ithlp2,"thlp2",
     .         "thlp2","K^2",zm)
          k = k + 1

        case ('rtpthlp')
          irtpthlp = k
          call stat_assign(irtpthlp,"rtpthlp",
     .         "rtpthlp","(kg K)/kg",zm)
          k = k + 1

        case ('wprtp')
          iwprtp = k

          call stat_assign(iwprtp,"wprtp",
     .         "wprtp","(m kg)/(s kg)",zm)
          k = k + 1

        case ('wpthlp')
          iwpthlp = k

          call stat_assign(iwpthlp,"wpthlp",
     .         "wpthlp","(m K)/s",zm)
          k = k + 1

        case ('wp4')
          iwp4 = k
          call stat_assign(iwp4,"wp4",
     .         "wp4","(m^4)/(s^4)",zm)
          k = k + 1

        case ('wpthvp')
          iwpthvp = k
          call stat_assign(iwpthvp,"wpthvp",
     .         "Buoyancy flux (K m/s)","(K m)/s",zm)
          k = k + 1

        case ('rtpthvp')
          irtpthvp = k
          call stat_assign(irtpthvp,"rtpthvp",
     .         "rtpthvp","(kg K)/kg",zm)
          k = k + 1

        case ('thlpthvp')
          ithlpthvp = k
          call stat_assign(ithlpthvp,"thlpthvp",
     .         "thlpthvp","K^2",zm)
          k = k + 1

        case ('taum')
          itaum = k

          call stat_assign(itaum,"taum",
     .         "Dissipation time","s",zm)
          k = k + 1

        case ('khm')
          iKhm = k

          call stat_assign(iKhm,"Khm",
     .         "Eddy diffusivity","m^2/s",zm)
          k = k + 1

        case ('wprcp')
          iwprcp = k
          call stat_assign(iwprcp,"wprcp",
     .         "wprcp","(m kg)/(s kg)",zm)
          k = k + 1

        case ('thlprcp')
          ithlprcp = k
          call stat_assign(ithlprcp,"thlprcp",
     .         "thlprcp","(K kg)/(kg)",zm)
          k = k + 1

        case ('rtprcp')
          irtprcp = k

          call stat_assign(irtprcp,"rtprcp",
     .         "rtprcp","(kg^2)/(kg^2)",zm)
          k = k + 1

        case ('rcp2')
          ircp2 = k
          call stat_assign(ircp2,"rcp2",
     .         "rcp2","(kg^2)/(kg^2)",zm)
          k = k + 1
        case ('upwp')
          iupwp = k
          call stat_assign(iupwp,"upwp",
     .         "upwp","m^2/s^2",zm)
          k = k + 1
        case ('vpwp')
          ivpwp = k
          call stat_assign(ivpwp,"vpwp",
     .         "vpwp","m^2/s^2",zm)
          k = k + 1
        case ('rhom')
          irhom = k
          call stat_assign(irhom,"rhom",
     .         "density","kg/(m^3)",zm)
          k = k + 1
        case ('sc')
          iscm = k
          call stat_assign(iscm,"scm",
     .         "scm","count",zm)
          k = k + 1
        case ('em')
          iem = k
          call stat_assign(iem,"em",
     .         "em","m2/s2",zm)
          k = k + 1
        case ('shear')      ! Brian
          ishear = k
          call stat_assign(ishear,"shear",
     .         "wind shear term (m^2/s^3)","m^2/s^3",zm)
          k = k + 1
        case ('Frad')
          iFrad = k
          call stat_assign(iFrad,"Frad",
     .         "radiative flux","W/m^2",zm)
          k = k + 1
        case ('Frad_LW')    ! Brian
          iFrad_LW = k
          call stat_assign(iFrad_LW,"Frad_LW",
     .         "Long-wave radiative flux","W/m^2",zm)
          k = k + 1
        case ('Frad_SW')    ! Brian
          iFrad_SW = k

          call stat_assign(iFrad_SW,"Frad_SW",
     .         "Short-wave radiative flux","W/m^2",zm)
          k = k + 1

        case ('Fprec')      ! Brian
          iFprec = k

          call stat_assign(iFprec,"Fprec",
     .         "precipitation flux (W/m^2)","W/m^2",zm)
          k = k + 1
 
        case ('Fcsed')      ! Brian
          iFcsed = k

          call stat_assign(iFcsed,"Fcsed",
     .         "cloud water sedimentation flux (kg/(s*m^2))",
     .         "kg/(s m2)",zm)
          k = k + 1

        case ('Vrr')           ! Brian
          iVrr = k

          call stat_assign(iVrr,"Vrr",
     .         "rrm sedimentation velocity (m/s)","m/s",zm)
          k = k + 1

        case ('VNr')           ! Brian
          iVNr = k

          call stat_assign(iVNr,"VNr",
     .         "Nrm sedimentation velocity (m/s)","m/s",zm)
          k = k + 1

        case ('Vsnow')
          iVsnow = k

          call stat_assign(iVsnow,"Vsnow",
     .         "Snow sedimentation velocity (m/s)","m/s",zm)
          k = k + 1

        case ('Vgraupel')
          iVgraupel = k

          call stat_assign(iVgraupel,"Vgraupel",
     .         "Graupel sedimentation velocity (m/s)","m/s",zm)
          k = k + 1

        case ('Vice')
          iVice = k

          call stat_assign(iVice,"Vice",
     .         "Pristine ice sedimentation velocity (m/s)","m/s",zm)
          k = k + 1
 
        case ('wp2_bt')
          iwp2_bt = k

          call stat_assign(iwp2_bt,"wp2_bt",
     .         "wp2 budget","m2/s3",zm)
          k = k + 1

        case ('wp2_ma')
          iwp2_ma = k

          call stat_assign(iwp2_ma,"wp2_ma",
     .         "wp2 mean advection","m2/s3",zm)
          k = k + 1

        case ('wp2_ta')
          iwp2_ta = k

          call stat_assign(iwp2_ta,"wp2_ta",
     .         "wp2 turbulent advection","m2/s3",zm)
          k = k + 1

        case ('wp2_ac')
          iwp2_ac = k

          call stat_assign(iwp2_ac,"wp2_ac",
     .         "wp2 accumulation term","m2/s3",zm)
          k = k + 1

        case ('wp2_bp')
          iwp2_bp = k
     
          call stat_assign(iwp2_bp,"wp2_bp",
     .         "wp2 buoyancy production","m2/s3",zm)
          k = k + 1

        case ('wp2_pr1')
          iwp2_pr1 = k

          call stat_assign(iwp2_pr1,"wp2_pr1",
     .         "wp2 pressure term 1","m2/s3",zm)
          k = k + 1

        case ('wp2_pr2')
          iwp2_pr2 = k
          call stat_assign(iwp2_pr2,"wp2_pr2",
     .         "wp2 pressure term 2","m2/s3",zm)
          k = k + 1

        case ('wp2_pr3')
          iwp2_pr3 = k
          call stat_assign(iwp2_pr3,"wp2_pr3",
     .         "wp2 pressure term 3","m2/s3",zm)

          k = k + 1

        case ('wp2_dp1')
          iwp2_dp1 = k
          call stat_assign(iwp2_dp1,"wp2_dp1",
     .         "wp2 dissipation term 1","m2/s3",zm)
          k = k + 1

        case ('wp2_dp2')
          iwp2_dp2 = k
          call stat_assign(iwp2_dp2,"wp2_dp2",
     .         "wp2 dissipation term 2","m2/s3",zm)

          k = k + 1

        case ('wp2_cl')
          iwp2_cl = k

          call stat_assign(iwp2_cl,"wp2_cl",
     .         "wp2 clipping term","m2/s3",zm)

          k = k + 1

        case ('wp2_pd')
          iwp2_pd = k

          call stat_assign(iwp2_pd,"wp2_pd",
     .         "wp2 positive definite adjustment","m2/s3",zm)

          k = k + 1

        case ('wprtp_bt')
          iwprtp_bt = k
          call stat_assign(iwprtp_bt,"wprtp_bt",
     .         "wprtp budget","(m kg)/(s2 kg)",zm)
          k = k + 1

        case ('wprtp_ma')
          iwprtp_ma = k

          call stat_assign(iwprtp_ma,"wprtp_ma",
     .         "wprtp mean advection","(m kg)/(s2 kg)",zm)
          k = k + 1

        case ('wprtp_ta')
          iwprtp_ta = k

          call stat_assign(iwprtp_ta,"wprtp_ta",
     .         "wprtp turbulent advection","(m kg)/(s2 kg)",zm)
          k = k + 1

        case ('wprtp_tp')
          iwprtp_tp = k

          call stat_assign(iwprtp_tp,"wprtp_tp",
     .         "wprtp turbulent production","(m kg)/(s2 kg)",zm)
          k = k + 1

        case ('wprtp_ac')
          iwprtp_ac = k

          call stat_assign(iwprtp_ac,"wprtp_ac",
     .         "wprtp accumulation term","(m kg)/(s2 kg)",zm)
          k = k + 1

        case ('wprtp_bp')
          iwprtp_bp = k

          call stat_assign(iwprtp_bp,"wprtp_bp",
     .         "wprtp buoyancy production","(m kg)/(s2 kg)",zm)
          k = k + 1

        case ('wprtp_pr1')
          iwprtp_pr1 = k

          call stat_assign(iwprtp_pr1,"wprtp_pr1",
     .         "wprtp pressure term 1","(m kg)/(s2 kg)",zm)
          k = k + 1

        case ('wprtp_pr2')
          iwprtp_pr2 = k

          call stat_assign(iwprtp_pr2,"wprtp_pr2",
     .         "wprtp pressure term 2","(m kg)/(s2 kg)",zm)
          k = k + 1

        case ('wprtp_pr3')
          iwprtp_pr3 = k

          call stat_assign(iwprtp_pr3,"wprtp_pr3",
     .         "wprtp pressure term 3","(m kg)/(s2 kg)",zm)
          k = k + 1

        case ('wprtp_dp1')
          iwprtp_dp1 = k

          call stat_assign(iwprtp_dp1,"wprtp_dp1",
     .         "wprtp dissipation term 1","(m kg)/(s2 kg)",zm)
          k = k + 1

        case ('wprtp_cl')
          iwprtp_cl = k

          call stat_assign(iwprtp_cl,"wprtp_cl",
     .         "wprtp clipping term","(m kg)/(s2 kg)",zm)
          k = k + 1

        case ('wprtp_sicl')
          iwprtp_sicl = k

          call stat_assign(iwprtp_sicl,"wprtp_sicl",
     .         "wprtp semi-implicit clipping term","(m kg)/(s2 kg)",zm)
          k = k + 1

        case ('wprtp_pd')
          iwprtp_pd = k

          call stat_assign(iwprtp_pd,"wprtp_pd",
     .         "wprtp flux corrected trans. term","(m kg)/(s2 kg)",zm)
          k = k + 1

        case ('wpthlp_bt')
          iwpthlp_bt = k

          call stat_assign(iwpthlp_bt,"wpthlp_bt",
     .         "wpthlp budget","(m K)/s2",zm)
          k = k + 1

        case ('wpthlp_ma')
          iwpthlp_ma = k
          call stat_assign(iwpthlp_ma,"wpthlp_ma",
     .         "wpthlp mean advection","(m K)/s2",zm)

          k = k + 1

        case ('wpthlp_ta')
          iwpthlp_ta = k
          call stat_assign(iwpthlp_ta,"wpthlp_ta",
     .         "wpthlp turbulent advection","(m K)/s2",zm)

          k = k + 1

        case ('wpthlp_tp')
          iwpthlp_tp = k
          call stat_assign(iwpthlp_tp,"wpthlp_tp",
     .         "wpthlp turbulent production","(m K)/s2",zm)

          k = k + 1

        case ('wpthlp_ac')
          iwpthlp_ac = k
          call stat_assign(iwpthlp_ac,"wpthlp_ac",
     .         "wpthlp accumulation term","(m K)/s2",zm)

          k = k + 1

        case ('wpthlp_bp')
          iwpthlp_bp = k
          call stat_assign(iwpthlp_bp,"wpthlp_bp",
     .         "wpthlp buoyancy production","(m K)/s2",zm)
          k = k + 1

        case ('wpthlp_pr1')
          iwpthlp_pr1 = k

          call stat_assign(iwpthlp_pr1,"wpthlp_pr1",
     .         "wpthlp pressure term 1","(m K)/s2",zm)
          k = k + 1

        case ('wpthlp_pr2')
          iwpthlp_pr2 = k

          call stat_assign(iwpthlp_pr2,"wpthlp_pr2",
     .         "wpthlp pressure term 2","(m K)/s2",zm)
          k = k + 1

        case ('wpthlp_pr3')
          iwpthlp_pr3 = k
          call stat_assign(iwpthlp_pr3,"wpthlp_pr3",
     .         "wpthlp pressure term 3","(m K)/s2",zm)
          k = k + 1

        case ('wpthlp_dp1')
          iwpthlp_dp1 = k
          call stat_assign(iwpthlp_dp1,"wpthlp_dp1",
     .         "wpthlp dissipation term 1","(m K)/s2",zm)
          k = k + 1

        case ('wpthlp_cl')
          iwpthlp_cl = k
          call stat_assign(iwpthlp_cl,"wpthlp_cl",
     .         "wpthlp clipping term","(m K)/s2",zm)
          k = k + 1

        case ('wpthlp_sicl')
          iwpthlp_sicl = k
          call stat_assign(iwpthlp_sicl,"wpthlp_sicl",
     .         "wpthlp semi-implicit clipping term","(m K)/s2",zm)
          k = k + 1

        ! Variance budgets
        case ('rtp2_bt')
          irtp2_bt = k
          call stat_assign(irtp2_bt,"rtp2_bt",
     .         "rtp2 budget","kg/(kg s)",zm)
          k = k + 1
        case ('rtp2_ma')
          irtp2_ma = k
          call stat_assign(irtp2_ma,"rtp2_ma",
     .         "rtp2 mean advection","kg/(kg s)",zm)
          k = k + 1
        case ('rtp2_ta')
          irtp2_ta = k
          call stat_assign(irtp2_ta,"rtp2_ta",
     .         "rtp2 turbulent advection","kg/(kg s)",zm)
          k = k + 1
        case ('rtp2_tp')
          irtp2_tp = k
          call stat_assign(irtp2_tp,"rtp2_tp",
     .         "rtp2 turbulent production","kg/(kg s)",zm)
          k = k + 1
        case ('rtp2_dp1')
          irtp2_dp1 = k
          call stat_assign(irtp2_dp1,"rtp2_dp1",
     .         "rtp2 dissipation term 1","kg/(kg s)",zm)
          k = k + 1
        case ('rtp2_dp2')
          irtp2_dp2 = k
          call stat_assign(irtp2_dp2,"rtp2_dp2",
     .         "rtp2 dissipation term 2","kg/(kg s)",zm)
          k = k + 1
        case ('rtp2_cl')
          irtp2_cl = k
          call stat_assign(irtp2_cl,"rtp2_cl",
     .         "rtp2 clipping term","kg/(kg s)",zm)
          k = k + 1

        case ('rtp2_pd')
          irtp2_pd = k
          call stat_assign( irtp2_pd, "rtp2_pd",
     .         "rtp2 positive definite adjustment", "m^2/s^2", zm )
          k = k + 1

        case ('thlp2_bt')
          ithlp2_bt = k
          call stat_assign(ithlp2_bt,"thlp2_bt",
     .         "thlp2 budget","(K^2)/s",zm)
          k = k + 1
        case ('thlp2_ma')
          ithlp2_ma = k
          call stat_assign(ithlp2_ma,"thlp2_ma",
     .         "thlp2 mean advection","(K^2)/s",zm)
          k = k + 1
        case ('thlp2_ta')
          ithlp2_ta = k
          call stat_assign(ithlp2_ta,"thlp2_ta",
     .         "thlp2 turbulent advection","(K^2)/s",zm)
          k = k + 1
        case ('thlp2_tp')
          ithlp2_tp = k
          call stat_assign(ithlp2_tp,"thlp2_tp",
     .         "thlp2 turbulent production","(K^2)/s",zm)
          k = k + 1
        case ('thlp2_dp1')
          ithlp2_dp1 = k
          call stat_assign(ithlp2_dp1,"thlp2_dp1",
     .         "thlp2 dissipation term 1","(K^2)/s",zm)
          k = k + 1
        case ('thlp2_dp2')
          ithlp2_dp2 = k
          call stat_assign(ithlp2_dp2,"thlp2_dp2",
     .         "thlp2 dissipation term 2","(K^2)/s",zm)
          k = k + 1
        case ('thlp2_cl')
          ithlp2_cl = k
          call stat_assign(ithlp2_cl,"thlp2_cl",
     .         "thlp2 clipping term","(K^2)/s",zm)
          k = k + 1

        case ('thlp2_pd')
          ithlp2_pd = k
          call stat_assign( ithlp2_pd, "thlp2_pd",
     .         "thlp2 positive definite adjustment", "m^2/s^2", zm )
          k = k + 1

        case ('rtpthlp_bt')
          irtpthlp_bt = k
          call stat_assign(irtpthlp_bt,"rtpthlp_bt",
     .         "rtpthlp budget","(kg K)/(kg s)",zm)
          k = k + 1
        case ('rtpthlp_ma')
          irtpthlp_ma = k
          call stat_assign(irtpthlp_ma,"rtpthlp_ma",
     .         "rtpthlp mean advection","(kg K)/(kg s)",zm)
          k = k + 1
        case ('rtpthlp_ta')
          irtpthlp_ta = k
          call stat_assign(irtpthlp_ta,"rtpthlp_ta",
     .         "rtpthlp turbulent advection","(kg K)/(kg s)",zm)
          k = k + 1
        case ('rtpthlp_tp1')
          irtpthlp_tp1 = k
          call stat_assign(irtpthlp_tp1,"rtpthlp_tp1",
     .         "rtpthlp turbulent production 1","(kg K)/(kg s)",zm)
          k = k + 1
        case ('rtpthlp_tp2')
          irtpthlp_tp2 = k
          call stat_assign(irtpthlp_tp2,"rtpthlp_tp2",
     .         "rtpthlp turbulent production 2","(kg K)/(kg s)",zm)
          k = k + 1
        case ('rtpthlp_dp1')
          irtpthlp_dp1 = k
          call stat_assign(irtpthlp_dp1,"rtpthlp_dp1",
     .         "rtpthlp dissipation term 1","(kg K)/(kg s)",zm)
          k = k + 1
        case ('rtpthlp_dp2')
          irtpthlp_dp2 = k
          call stat_assign(irtpthlp_dp2,"rtpthlp_dp2",
     .         "rtpthlp dissipation term 2","(kg K)/(kg s)",zm)
          k = k + 1
        case ('rtpthlp_cl')
          irtpthlp_cl = k
          call stat_assign(irtpthlp_cl,"rtpthlp_cl",
     .         "rtpthlp clipping term","(kg K)/(kg s)",zm)
          k = k + 1

        case ('up2_plus_vp2')
          iup2_plus_vp2 = k
          call stat_assign(iup2_plus_vp2,"up2_plus_vp2",
     .         "up2_plus_vp2","m^2/s^2",zm)
          k = k + 1

        case ('up2')
          iup2 = k
          call stat_assign(iup2,"up2",
     .         "up2","m^2/s^2",zm)
          k = k + 1

        case ('vp2')
          ivp2 = k
          call stat_assign(ivp2,"vp2",
     .         "vp2","m^2/s^2",zm)
          k = k + 1

        case ('up2_bt')
          iup2_bt = k
          call stat_assign(iup2_bt,"up2_bt",
     .         "up2 time tendency","m^2/s^2",zm)
          k = k + 1

        case ('up2_ma')
          iup2_ma = k
          call stat_assign(iup2_ma,"up2_ma",
     .         "up2 mean advection","m^2/s^2",zm)
          k = k + 1

        case ('up2_ta')
          iup2_ta = k
          call stat_assign(iup2_ta,"up2_ta",
     .         "up2 turbulent advection","m^2/s^2",zm)
          k = k + 1

        case ('up2_tp')
          iup2_tp = k
          call stat_assign(iup2_tp,"up2_tp",
     .         "up2 turbulent production","m^2/s^2",zm)
          k = k + 1

        case ('up2_dp1')
          iup2_dp1 = k
          call stat_assign(iup2_dp1,"up2_dp1",
     .         "up2 dissipation term 1","m^2/s^2",zm)
          k = k + 1

        case ('up2_dp2')
          iup2_dp2 = k
          call stat_assign(iup2_dp2,"up2_dp2",
     .         "up2 dissipation term 2","m^2/s^2",zm)
          k = k + 1

        case ('up2_pr1')
          iup2_pr1 = k
          call stat_assign(iup2_pr1,"up2_pr1",
     .         "up2 pressure term 1","m^2/s^2",zm)
          k = k + 1

        case ('up2_pr2')
          iup2_pr2 = k
          call stat_assign(iup2_pr2,"up2_pr2",
     .         "up2 pressure term 2","m^2/s^2",zm)
          k = k + 1

        case ('up2_cl')
          iup2_cl = k
          call stat_assign(iup2_cl,"up2_cl",
     .         "up2 clipping","m^2/s^2",zm)
          k = k + 1

        case ('up2_pd')
          iup2_pd = k
          call stat_assign( iup2_pd, "up2_pd",
     .         "up2 positive definite adjustment", "m^2/s^2", zm )
          k = k + 1

        case ('vp2_bt')
          ivp2_bt = k
          call stat_assign(ivp2_bt,"vp2_bt",
     .         "vp2 time tendency","m^2/s^2",zm)
          k = k + 1

        case ('vp2_ma')
          ivp2_ma = k
          call stat_assign(ivp2_ma,"vp2_ma",
     .         "vp2 mean advection","m^2/s^2",zm)
          k = k + 1

        case ('vp2_ta')
          ivp2_ta = k
          call stat_assign(ivp2_ta,"vp2_ta",
     .         "vp2 turbulent advection","m^2/s^2",zm)
          k = k + 1

        case ('vp2_tp')
          ivp2_tp = k
          call stat_assign(ivp2_tp,"vp2_tp",
     .         "vp2 turbulent production","m^2/s^2",zm)
          k = k + 1

        case ('vp2_dp1')
          ivp2_dp1 = k
          call stat_assign(ivp2_dp1,"vp2_dp1",
     .         "vp2 dissipation term 1","m^2/s^2",zm)
          k = k + 1

        case ('vp2_dp2')
          ivp2_dp2 = k
          call stat_assign(ivp2_dp2,"vp2_dp2",
     .         "vp2 dissipation term 2","m^2/s^2",zm)
          k = k + 1

        case ('vp2_pr1')
          ivp2_pr1 = k
          call stat_assign(ivp2_pr1,"vp2_pr1",
     .         "vp2 pressure term 1","m^2/s^2",zm)
          k = k + 1

        case ('vp2_pr2')
          ivp2_pr2 = k
          call stat_assign(ivp2_pr2,"vp2_pr2",
     .         "vp2 pressure term 2","m^2/s^2",zm)
          k = k + 1

        case ('vp2_cl')
          ivp2_cl = k
          call stat_assign(ivp2_cl,"vp2_cl",
     .         "vp2 clipping","m^2/s^2",zm)
          k = k + 1

        case ('vp2_pd')
          ivp2_pd = k
          call stat_assign( ivp2_pd, "vp2_pd",
     .         "vp2 positive definite adjustment", "m^2/s^2", zm )
          k = k + 1

        case ('sclraprtp')
          isclraprtp = k

          call stat_assign(isclraprtp,"sclraprtp",
     .         "scalar(2)'rt'","unknown",zm)
          k = k + 1
        case ('sclrbprtp')
          isclrbprtp = k
          call stat_assign(isclrbprtp,"sclrbprtp",
     .         "scalar(2)'rt'","unknown",zm)

          k = k + 1
        case ('sclrap2')
          isclrap2 = k
          call stat_assign(isclrap2,"sclrap2",
     .         "scalar(1)'^2'","unknown",zm)
          k = k + 1
        case ('sclrbp2')
          isclrbp2 = k
          call stat_assign(isclrbp2,"sclrbp2",
     .         "scalar(2)'^2'","unknown",zm)
          k = k + 1
        case ('sclrapthvp')
          isclrapthvp = k

          call stat_assign(isclrapthvp,"sclrathvp",
     .         "scalar(1)'th_v'","unknown",zm)
          k = k + 1
        case ('sclrbpthvp')
          isclrbpthvp = k
          call stat_assign(isclrbpthvp,"sclrbthvp",
     .         "scalar(2)'th_v'","unknown",zm)
          k = k + 1
        case ('sclrapthlp')
          isclrapthlp = k

          call stat_assign(isclrapthlp,"sclrathlp",
     .         "scalar(1)'th_l'","unknown",zm)
          k = k + 1
        case ('sclrbpthlp')
          isclrbpthlp = k

          call stat_assign(isclrbpthlp,"sclrbthlp",
     .         "scalar(2)'th_l'","unknown",zm)
          k = k + 1
        case ('sclraprcp')
          isclraprcp = k

          call stat_assign(isclraprcp,"sclraprcp",
     .         "scalar(1)'rc'","unknown",zm)
          k = k + 1
        case ('sclrbprcp')
          isclrbprcp = k

          call stat_assign(isclrbprcp,"sclrbprcp",
     .         "scalar(2)'rc'","unknown",zm)
          k = k + 1
        case ('wpsclrap')
          iwpsclrap = k

          call stat_assign(iwpsclrap,"wpsclrap",
     .         "'w'scalar(1)","unknown",zm)
          k = k + 1
        case ('wpsclrbp')
          iwpsclrbp = k

          call stat_assign(iwpsclrbp,"wpsclrbp",
     .         "'w'scalar(2)","unknown",zm)
          k = k + 1
        case ('wpsclrap2')
          iwpsclrap2 = k

          call stat_assign(iwpsclrap2,"wpsclrap2",
     .         "'w'scalar(1)'^2'","unknown",zm)
          k = k + 1
        case ('wpsclrbp2')
          iwpsclrbp2 = k

          call stat_assign(iwpsclrbp2,"wpsclrbp2",
     .         "'w'scalar(2)'^2'","unknown",zm)
          k = k + 1
        case ('wp2sclrap')
          iwp2sclrap = k

          call stat_assign(iwp2sclrap,"wp2sclrap",
     .         "'w'^2 scalar(1)","unknown",zm)
          k = k + 1
        case ('wp2sclrbp')
          iwp2sclrbp = k

          call stat_assign(iwp2sclrbp,"wp2sclrbp",
     .         "'w'^2 scalar(2)","unknown",zm)
          k = k + 1
        case ('wpsclraprtp')
          iwpsclraprtp = k

          call stat_assign(iwpsclraprtp,"wpsclraprtp",
     .         "'w' scalar(1)'rt'","unknown",zm)
          k = k + 1
        case ('wpsclrbprtp')
          iwpsclrbprtp = k

          call stat_assign(iwpsclrbprtp,"wpsclrbprtp",
     .         "'w' scalar(2)'rt'","unknown",zm)
          k = k + 1
        case ('wpsclrapthlp')
          iwpsclrapthlp = k

          call stat_assign(iwpsclrapthlp,"wpsclrapthlp",
     .         "'w' scalar(1)'th_l'","unknown",zm)
          k = k + 1
        case ('wpsclrbpthlp')
          iwpsclrbpthlp = k

          call stat_assign(iwpsclrbpthlp,"wpsclrbpthlp",
     .         "'w' scalar(2)'th_l'","unknown",zm)
          k = k + 1
        case ('wpedsclrap')
          iwpedsclrap = k

          call stat_assign(iwpedsclrap,"wpedsclrap",
     .         "eddy scalar(1)'w'","unknown",zm)
          k = k + 1
        case ('wpedsclrbp')
          iwpedsclrbp = k

          call stat_assign(iwpedsclrbp,"wpedsclrbp",
     .         "eddy scalar(2)'w'","unknown",zm)
          k = k + 1

        case default
          write(0,*) 'Error: unrecognized variable in vars_zm: ',
     .       trim(vars_zm(i))
          lerror = .true.

        end select

      end do

      return
      
      end subroutine stats_init_zm

#endif /*STATS*/
      end module stats_zm
