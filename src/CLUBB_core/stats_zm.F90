!-----------------------------------------------------------------------
! $Id$
module stats_zm

  implicit none

  private ! Default Scope

  public :: stats_init_zm

  ! Constant parameters
  integer, parameter, public :: nvarmax_zm = 250  ! Maximum variables allowed

  contains

!-----------------------------------------------------------------------
  subroutine stats_init_zm( vars_zm, l_error )

! Description:
!   Initializes array indices for zm
!-----------------------------------------------------------------------

    use constants, only: &
        fstderr ! Constant(s)

    use stats_variables, only: & 
          zm, & 
          iwp2, & 
          irtp2, & 
          ithlp2, & 
          irtpthlp, & 
          iwprtp, & 
          iwpthlp, & 
          iwp4, & 
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
          irho_zm, & 
          isigma_sqd_w, & 
          iem, & 
          ishear, &
          imean_w_up, &
          imean_w_down, & 
          iFrad, & 
          iFrad_LW, & 
          iFrad_SW, & 
          iFrad_LW_up, & 
          iFrad_SW_up, & 
          iFrad_LW_down, & 
          iFrad_SW_down, & 
          iFprec, & 
          iFcsed

    use stats_variables, only: & 
          iup2, & 
          ivp2, & 
          iup2_bt, & 
          iup2_ta, & 
          iup2_tp, & 
          iup2_ma, & 
          iup2_dp1, & 
          iup2_dp2, & 
          iup2_pr1, & 
          iup2_pr2, & 
          iup2_cl, & 
          iup2_pd, & 
          ivp2_bt, & 
          ivp2_ta, & 
          ivp2_tp, & 
          ivp2_ma, & 
          ivp2_dp1, & 
          ivp2_dp2, & 
          ivp2_pr1, & 
          ivp2_pr2, & 
          ivp2_cl, & 
          ivp2_pd, & 
          iVrr, & 
          iVNr, & 
          iVice, & 
          iVgraupel, & 
          iVsnow

    use stats_variables, only: & 
          iwp2_bt, & 
          iwp2_ma, & 
          iwp2_ta, & 
          iwp2_ac, & 
          iwp2_bp, & 
          iwp2_pr1, & 
          iwp2_pr2, & 
          iwp2_pr3, & 
          iwp2_dp1, & 
          iwp2_dp2, &
          iwp2_4hd, & 
          iwp2_cl, & 
          iwp2_pd

    use stats_variables, only: & 
          iwprtp_bt, & 
          iwprtp_ma, & 
          iwprtp_ta, & 
          iwprtp_tp, & 
          iwprtp_ac, & 
          iwprtp_bp, & 
          iwprtp_pr1, & 
          iwprtp_pr2, & 
          iwprtp_pr3, & 
          iwprtp_dp1, & 
          iwprtp_mfl, & 
          iwprtp_cl, & 
          iwprtp_sicl, & 
          iwprtp_pd, & 
          iwpthlp_bt, & 
          iwpthlp_ma, & 
          iwpthlp_ta, & 
          iwpthlp_tp, & 
          iwpthlp_ac, & 
          iwpthlp_bp, & 
          iwpthlp_pr1, & 
          iwpthlp_pr2, & 
          iwpthlp_pr3, & 
          iwpthlp_dp1, & 
          iwpthlp_mfl, & 
          iwpthlp_cl, & 
          iwpthlp_sicl

    use stats_variables, only: & 
        irtp2_bt, & 
        irtp2_ma, & 
        irtp2_ta, & 
        irtp2_tp, & 
        irtp2_dp1, & 
        irtp2_dp2, & 
        irtp2_cl, & 
        irtp2_pd, & 
        ithlp2_bt, & 
        ithlp2_ma, & 
        ithlp2_ta, & 
        ithlp2_tp, & 
        ithlp2_dp1, & 
        ithlp2_dp2, & 
        ithlp2_cl, & 
        ithlp2_pd, & 
        irtpthlp_bt, & 
        irtpthlp_ma, & 
        irtpthlp_ta, & 
        irtpthlp_tp1, & 
        irtpthlp_tp2, & 
        irtpthlp_dp1, & 
        irtpthlp_dp2, & 
        irtpthlp_cl
    
    use stats_variables, only: & 
        iwpthlp_enter_mfl, &
        iwpthlp_exit_mfl, &
        iwpthlp_mfl_lower_lim, &
        iwpthlp_mfl_upper_lim, &
        iwprtp_enter_mfl, &
        iwprtp_exit_mfl, &
        iwprtp_mfl_lower_lim, &
        iwprtp_mfl_upper_lim

    use stats_variables, only: & 
        isclrprtp, & 
        isclrp2, & 
        isclrpthvp, & 
        isclrpthlp, & 
        isclrprcp, & 
        iwpsclrp, & 
        iwp2sclrp, & 
        iwpsclrp2, & 
        iwpsclrprtp, & 
        iwpsclrpthlp, & 
        iwpedsclrp

    use stats_type, only: & 
        stat_assign ! Procedure

    use parameters_model, only: &
        sclr_dim, &
        edsclr_dim

!   use error_code, only: &
!       clubb_at_least_debug_level ! Function

    implicit none

    ! Input Variable
    ! zm variable names

    character(len= * ), dimension(nvarmax_zm), intent(in) :: vars_zm

    ! Output Variable
    logical, intent(inout) :: l_error

    ! Local Varables
    integer :: i,j, k

    logical :: l_found

    character(len=50) :: sclr_idx

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
    itau_zm   = 0
    iKh_zm    = 0
    iwprcp    = 0
    ithlprcp  = 0
    irtprcp   = 0
    ircp2     = 0
    iupwp     = 0
    ivpwp     = 0
    irho_zm   = 0
    isigma_sqd_w  = 0
    iem           = 0
    ishear        = 0  ! Brian
    imean_w_up    = 0
    imean_w_down  = 0
    iFrad         = 0
    iFrad_LW      = 0  ! Brian
    iFrad_SW      = 0  ! Brian
    iFrad_LW_up   = 0  ! Brian
    iFrad_SW_up   = 0  ! Brian
    iFrad_LW_down = 0  ! Brian
    iFrad_SW_down = 0  ! Brian
    iFprec        = 0  ! Brian
    iFcsed        = 0  ! Brian


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
    iwp2_4hd  = 0
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
    iwprtp_mfl  = 0
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
    iwpthlp_mfl  = 0
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

    !Monatonic flux limiter diagnostic output
    iwpthlp_mfl_lower_lim = 0
    iwpthlp_mfl_upper_lim = 0
    iwpthlp_enter_mfl = 0
    iwpthlp_exit_mfl = 0
    iwprtp_mfl_lower_lim = 0
    iwprtp_mfl_upper_lim = 0
    iwprtp_enter_mfl = 0
    iwprtp_exit_mfl = 0

    allocate(isclrprtp(1:sclr_dim))
    allocate(isclrp2(1:sclr_dim))
    allocate(isclrpthvp(1:sclr_dim))
    allocate(isclrpthlp(1:sclr_dim))
    allocate(isclrprcp(1:sclr_dim))
    allocate(iwpsclrp(1:sclr_dim))
    allocate(iwp2sclrp(1:sclr_dim))
    allocate(iwpsclrp2(1:sclr_dim))
    allocate(iwpsclrprtp(1:sclr_dim))
    allocate(iwpsclrpthlp(1:sclr_dim))

    allocate(iwpedsclrp(1:edsclr_dim))

!     Assign pointers for statistics variables zm

    isclrprtp    = 0
    isclrp2      = 0
    isclrpthvp   = 0
    isclrpthlp   = 0
    isclrprcp    = 0
    iwpsclrp     = 0
    iwp2sclrp    = 0
    iwpsclrp2    = 0
    iwpsclrprtp  = 0
    iwpsclrpthlp = 0

    iwpedsclrp   = 0

!     Assign pointers for statistics variables zm

    k = 1
    do i=1,zm%nn

      select case ( trim(vars_zm(i)) )

      case ('wp2')
        iwp2 = k
        call stat_assign(iwp2,"wp2", & 
             "wp2","m^2/s^2",zm)
        k = k + 1

      case ('rtp2')
        irtp2 = k
        call stat_assign(irtp2,"rtp2", & 
             "rtp2","(kg/kg)^2",zm)
        k = k + 1

      case ('thlp2')
        ithlp2 = k
        call stat_assign(ithlp2,"thlp2", & 
             "thlp2","K^2",zm)
        k = k + 1

      case ('rtpthlp')
        irtpthlp = k
        call stat_assign(irtpthlp,"rtpthlp", & 
             "rtpthlp","(kg K)/kg",zm)
        k = k + 1

      case ('wprtp')
        iwprtp = k

        call stat_assign(iwprtp,"wprtp", & 
             "wprtp","(m kg)/(s kg)",zm)
        k = k + 1

      case ('wpthlp')
        iwpthlp = k

        call stat_assign(iwpthlp,"wpthlp", & 
             "wpthlp","(m K)/s",zm)
        k = k + 1

      case ('wp4')
        iwp4 = k
        call stat_assign(iwp4,"wp4", & 
             "wp4","(m^4)/(s^4)",zm)
        k = k + 1

      case ('wpthvp')
        iwpthvp = k
        call stat_assign(iwpthvp,"wpthvp", & 
             "Buoyancy flux (K m/s)","(K m)/s",zm)
        k = k + 1

      case ('rtpthvp')
        irtpthvp = k
        call stat_assign(irtpthvp,"rtpthvp", & 
             "rtpthvp","(kg K)/kg",zm)
        k = k + 1

      case ('thlpthvp')
        ithlpthvp = k
        call stat_assign(ithlpthvp,"thlpthvp", & 
             "thlpthvp","K^2",zm)
        k = k + 1

      case ('tau_zm')
        itau_zm = k

        call stat_assign(itau_zm,"tau_zm", & 
             "Dissipation time","s",zm)
        k = k + 1

      case ('Kh_zm')
        iKh_zm = k

        call stat_assign(iKh_zm,"Kh_zm", & 
             "Eddy diffusivity","m^2/s",zm)
        k = k + 1

      case ('wprcp')
        iwprcp = k
        call stat_assign(iwprcp,"wprcp", & 
             "wprcp","(m kg)/(s kg)",zm)
        k = k + 1

      case ('thlprcp')
        ithlprcp = k
        call stat_assign(ithlprcp,"thlprcp", & 
             "thlprcp","(K kg)/(kg)",zm)
        k = k + 1

      case ('rtprcp')
        irtprcp = k

        call stat_assign(irtprcp,"rtprcp", & 
             "rtprcp","(kg^2)/(kg^2)",zm)
        k = k + 1

      case ('rcp2')
        ircp2 = k
        call stat_assign(ircp2,"rcp2", & 
             "rcp2","(kg^2)/(kg^2)",zm)
        k = k + 1
      case ('upwp')
        iupwp = k
        call stat_assign(iupwp,"upwp", & 
             "upwp","m^2/s^2",zm)
        k = k + 1
      case ('vpwp')
        ivpwp = k
        call stat_assign(ivpwp,"vpwp", & 
             "vpwp","m^2/s^2",zm)
        k = k + 1
      case ('rho_zm')
        irho_zm = k
        call stat_assign(irho_zm,"rho_zm", & 
             "density","kg/(m^3)",zm)
        k = k + 1
      case ('sigma_sqd_w')
        isigma_sqd_w = k
        call stat_assign(isigma_sqd_w,"sigma_sqd_w", & 
             "sigma_sqd_w","count",zm)
        k = k + 1
      case ('em')
        iem = k
        call stat_assign(iem,"em", & 
             "em","m2/s2",zm)
        k = k + 1
      case ('shear')      ! Brian
        ishear = k
        call stat_assign(ishear,"shear", & 
             "wind shear term (m^2/s^3)","m^2/s^3",zm)
        k = k + 1
      case ('mean_w_up')
        imean_w_up = k
        call stat_assign(imean_w_up, "mean_w_up", & 
             "mean of all values of w >= 0", "m/s", zm)
        k = k + 1
      case ('mean_w_down')
        imean_w_down = k
        call stat_assign(imean_w_down, "mean_w_down", & 
             "mean of all values of w <= 0", "m/s", zm)
        k = k + 1
      case ('Frad')
        iFrad = k
        call stat_assign(iFrad,"Frad", & 
             "radiative flux","W/m^2",zm)
        k = k + 1
      case ('Frad_LW')    ! Brian
        iFrad_LW = k
        call stat_assign(iFrad_LW,"Frad_LW", & 
             "Net Long-wave radiative flux (W/m^2)","W/m^2",zm)
        k = k + 1
      case ('Frad_SW')    ! Brian
        iFrad_SW = k

        call stat_assign(iFrad_SW,"Frad_SW", & 
             "Net Short-wave radiative flux (W/m^2)","W/m^2",zm)
        k = k + 1

      case ('Frad_LW_up')    ! Brian
        iFrad_LW_up = k
        call stat_assign(iFrad_LW_up,"Frad_LW_up", & 
             "Long-wave upwelling radiative flux (W/m^2)","W/m^2",zm)
        k = k + 1
      case ('Frad_SW_up')    ! Brian
        iFrad_SW_up = k

        call stat_assign(iFrad_SW_up,"Frad_SW_up", & 
             "Short-wave upwelling radiative flux (W/m^2)","W/m^2",zm)
        k = k + 1

      case ('Frad_LW_down')    ! Brian
        iFrad_LW_down = k
        call stat_assign(iFrad_LW_down,"Frad_LW_down", & 
        "Long-wave downwelling radiative flux, defined positive down (W/m^2)", "W/m^2", zm )
        k = k + 1
      case ('Frad_SW_down')    ! Brian
        iFrad_SW_down = k

        call stat_assign(iFrad_SW_down,"Frad_SW_down", & 
        "Short-wave downwelling radiative flux, defined positive down (W/m^2)", "W/m^2", zm )
        k = k + 1


      case ('Fprec')      ! Brian
        iFprec = k

        call stat_assign(iFprec,"Fprec", & 
             "precipitation flux (W/m^2)","W/m^2",zm)
        k = k + 1

      case ('Fcsed')      ! Brian
        iFcsed = k

        call stat_assign(iFcsed,"Fcsed", & 
             "cloud water sedimentation flux (kg/(s*m^2))", & 
             "kg/(s m2)",zm)
        k = k + 1

      case ('Vrr')           ! Brian
        iVrr = k

        call stat_assign(iVrr,"Vrr", & 
             "rrainm sedimentation velocity (m/s)","m/s",zm)
        k = k + 1

      case ('VNr')           ! Brian
        iVNr = k

        call stat_assign(iVNr,"VNr", & 
             "Nrm sedimentation velocity (m/s)","m/s",zm)
        k = k + 1

      case ('Vsnow')
        iVsnow = k

        call stat_assign(iVsnow,"Vsnow", & 
             "Snow sedimentation velocity (m/s)","m/s",zm)
        k = k + 1

      case ('Vgraupel')
        iVgraupel = k

        call stat_assign(iVgraupel,"Vgraupel", & 
             "Graupel sedimentation velocity (m/s)","m/s",zm)
        k = k + 1

      case ('Vice')
        iVice = k

        call stat_assign(iVice,"Vice", & 
             "Pristine ice sedimentation velocity (m/s)","m/s",zm)
        k = k + 1

      case ('wp2_bt')
        iwp2_bt = k

        call stat_assign(iwp2_bt,"wp2_bt", & 
             "wp2 budget","m2/s3",zm)
        k = k + 1

      case ('wp2_ma')
        iwp2_ma = k

        call stat_assign(iwp2_ma,"wp2_ma", & 
             "wp2 mean advection","m2/s3",zm)
        k = k + 1

      case ('wp2_ta')
        iwp2_ta = k

        call stat_assign(iwp2_ta,"wp2_ta", & 
             "wp2 turbulent advection","m2/s3",zm)
        k = k + 1

      case ('wp2_ac')
        iwp2_ac = k

        call stat_assign(iwp2_ac,"wp2_ac", & 
             "wp2 accumulation term","m2/s3",zm)
        k = k + 1

      case ('wp2_bp')
        iwp2_bp = k

        call stat_assign(iwp2_bp,"wp2_bp", & 
             "wp2 buoyancy production","m2/s3",zm)
        k = k + 1

      case ('wp2_pr1')
        iwp2_pr1 = k

        call stat_assign(iwp2_pr1,"wp2_pr1", & 
             "wp2 pressure term 1","m2/s3",zm)
        k = k + 1

      case ('wp2_pr2')
        iwp2_pr2 = k
        call stat_assign(iwp2_pr2,"wp2_pr2", & 
             "wp2 pressure term 2","m2/s3",zm)
        k = k + 1

      case ('wp2_pr3')
        iwp2_pr3 = k
        call stat_assign(iwp2_pr3,"wp2_pr3", & 
             "wp2 pressure term 3","m2/s3",zm)

        k = k + 1

      case ('wp2_dp1')
        iwp2_dp1 = k
        call stat_assign(iwp2_dp1,"wp2_dp1", & 
             "wp2 dissipation term 1","m2/s3",zm)
        k = k + 1

      case ('wp2_dp2')
        iwp2_dp2 = k
        call stat_assign(iwp2_dp2,"wp2_dp2", & 
             "wp2 dissipation term 2","m2/s3",zm)

        k = k + 1

      case ('wp2_4hd')
        iwp2_4hd = k
        call stat_assign(iwp2_4hd,"wp2_4hd", & 
             "wp2 4th-order hyper-diffusion","m2/s3",zm)

        k = k + 1

      case ('wp2_cl')
        iwp2_cl = k

        call stat_assign(iwp2_cl,"wp2_cl", & 
             "wp2 clipping term","m2/s3",zm)

        k = k + 1

      case ('wp2_pd')
        iwp2_pd = k

        call stat_assign(iwp2_pd,"wp2_pd", & 
             "wp2 positive definite adjustment","m2/s3",zm)

        k = k + 1

      case ('wprtp_bt')
        iwprtp_bt = k
        call stat_assign(iwprtp_bt,"wprtp_bt", & 
             "wprtp budget","(m kg)/(s2 kg)",zm)
        k = k + 1

      case ('wprtp_ma')
        iwprtp_ma = k

        call stat_assign(iwprtp_ma,"wprtp_ma", & 
             "wprtp mean advection","(m kg)/(s2 kg)",zm)
        k = k + 1

      case ('wprtp_ta')
        iwprtp_ta = k

        call stat_assign(iwprtp_ta,"wprtp_ta", & 
             "wprtp turbulent advection","(m kg)/(s2 kg)",zm)
        k = k + 1

      case ('wprtp_tp')
        iwprtp_tp = k

        call stat_assign(iwprtp_tp,"wprtp_tp", & 
             "wprtp turbulent production","(m kg)/(s2 kg)",zm)
        k = k + 1

      case ('wprtp_ac')
        iwprtp_ac = k

        call stat_assign(iwprtp_ac,"wprtp_ac", & 
             "wprtp accumulation term","(m kg)/(s2 kg)",zm)
        k = k + 1

      case ('wprtp_bp')
        iwprtp_bp = k

        call stat_assign(iwprtp_bp,"wprtp_bp", & 
             "wprtp buoyancy production","(m kg)/(s2 kg)",zm)
        k = k + 1

      case ('wprtp_pr1')
        iwprtp_pr1 = k

        call stat_assign(iwprtp_pr1,"wprtp_pr1", & 
             "wprtp pressure term 1","(m kg)/(s2 kg)",zm)
        k = k + 1

      case ('wprtp_pr2')
        iwprtp_pr2 = k

        call stat_assign(iwprtp_pr2,"wprtp_pr2", & 
             "wprtp pressure term 2","(m kg)/(s2 kg)",zm)
        k = k + 1

      case ('wprtp_pr3')
        iwprtp_pr3 = k

        call stat_assign(iwprtp_pr3,"wprtp_pr3", & 
             "wprtp pressure term 3","(m kg)/(s2 kg)",zm)
        k = k + 1

      case ('wprtp_dp1')
        iwprtp_dp1 = k

        call stat_assign(iwprtp_dp1,"wprtp_dp1", & 
             "wprtp dissipation term 1","(m kg)/(s2 kg)",zm)
        k = k + 1

      case ('wprtp_mfl')
        iwprtp_mfl = k

        call stat_assign(iwprtp_mfl,"wprtp_mfl", & 
             "wprtp monotonic flux limiter","(m kg)/(s2 kg)",zm)
        k = k + 1

      case ('wprtp_cl')
        iwprtp_cl = k

        call stat_assign(iwprtp_cl,"wprtp_cl", & 
             "wprtp clipping term","(m kg)/(s2 kg)",zm)
        k = k + 1

      case ('wprtp_sicl')
        iwprtp_sicl = k

        call stat_assign(iwprtp_sicl,"wprtp_sicl", & 
             "wprtp semi-implicit clipping term","(m kg)/(s2 kg)",zm)
        k = k + 1

      case ('wprtp_pd')
        iwprtp_pd = k

        call stat_assign(iwprtp_pd,"wprtp_pd", & 
             "wprtp flux corrected trans. term","(m kg)/(s2 kg)",zm)
        k = k + 1

      case ('wpthlp_bt')
        iwpthlp_bt = k

        call stat_assign(iwpthlp_bt,"wpthlp_bt", & 
             "wpthlp budget","(m K)/s2",zm)
        k = k + 1

      case ('wpthlp_ma')
        iwpthlp_ma = k
        call stat_assign(iwpthlp_ma,"wpthlp_ma", & 
             "wpthlp mean advection","(m K)/s2",zm)

        k = k + 1

      case ('wpthlp_ta')
        iwpthlp_ta = k
        call stat_assign(iwpthlp_ta,"wpthlp_ta", & 
             "wpthlp turbulent advection","(m K)/s2",zm)

        k = k + 1

      case ('wpthlp_tp')
        iwpthlp_tp = k
        call stat_assign(iwpthlp_tp,"wpthlp_tp", & 
             "wpthlp turbulent production","(m K)/s2",zm)

        k = k + 1

      case ('wpthlp_ac')
        iwpthlp_ac = k
        call stat_assign(iwpthlp_ac,"wpthlp_ac", & 
             "wpthlp accumulation term","(m K)/s2",zm)

        k = k + 1

      case ('wpthlp_bp')
        iwpthlp_bp = k
        call stat_assign(iwpthlp_bp,"wpthlp_bp", & 
             "wpthlp buoyancy production","(m K)/s2",zm)
        k = k + 1

      case ('wpthlp_pr1')
        iwpthlp_pr1 = k

        call stat_assign(iwpthlp_pr1,"wpthlp_pr1", & 
             "wpthlp pressure term 1","(m K)/s2",zm)
        k = k + 1

      case ('wpthlp_pr2')
        iwpthlp_pr2 = k

        call stat_assign(iwpthlp_pr2,"wpthlp_pr2", & 
             "wpthlp pressure term 2","(m K)/s2",zm)
        k = k + 1

      case ('wpthlp_pr3')
        iwpthlp_pr3 = k
        call stat_assign(iwpthlp_pr3,"wpthlp_pr3", & 
             "wpthlp pressure term 3","(m K)/s2",zm)
        k = k + 1

      case ('wpthlp_dp1')
        iwpthlp_dp1 = k
        call stat_assign(iwpthlp_dp1,"wpthlp_dp1", & 
             "wpthlp dissipation term 1","(m K)/s2",zm)
        k = k + 1

      case ('wpthlp_mfl')
        iwpthlp_mfl = k
        call stat_assign(iwpthlp_mfl,"wpthlp_mfl", & 
             "wpthlp monotonic flux limiter","(m K)/s2",zm)
        k = k + 1

      case ('wpthlp_cl')
        iwpthlp_cl = k
        call stat_assign(iwpthlp_cl,"wpthlp_cl", & 
             "wpthlp clipping term","(m K)/s2",zm)
        k = k + 1

      case ('wpthlp_sicl')
        iwpthlp_sicl = k
        call stat_assign(iwpthlp_sicl,"wpthlp_sicl", & 
             "wpthlp semi-implicit clipping term","(m K)/s2",zm)
        k = k + 1

        ! Variance budgets
      case ('rtp2_bt')
        irtp2_bt = k
        call stat_assign(irtp2_bt,"rtp2_bt", & 
             "rtp2 budget","kg/(kg s)",zm)
        k = k + 1
      case ('rtp2_ma')
        irtp2_ma = k
        call stat_assign(irtp2_ma,"rtp2_ma", & 
             "rtp2 mean advection","kg/(kg s)",zm)
        k = k + 1
      case ('rtp2_ta')
        irtp2_ta = k
        call stat_assign(irtp2_ta,"rtp2_ta", & 
             "rtp2 turbulent advection","kg/(kg s)",zm)
        k = k + 1
      case ('rtp2_tp')
        irtp2_tp = k
        call stat_assign(irtp2_tp,"rtp2_tp", & 
             "rtp2 turbulent production","kg/(kg s)",zm)
        k = k + 1
      case ('rtp2_dp1')
        irtp2_dp1 = k
        call stat_assign(irtp2_dp1,"rtp2_dp1", & 
             "rtp2 dissipation term 1","kg/(kg s)",zm)
        k = k + 1
      case ('rtp2_dp2')
        irtp2_dp2 = k
        call stat_assign(irtp2_dp2,"rtp2_dp2", & 
             "rtp2 dissipation term 2","kg/(kg s)",zm)
        k = k + 1
      case ('rtp2_cl')
        irtp2_cl = k
        call stat_assign(irtp2_cl,"rtp2_cl", & 
             "rtp2 clipping term","kg/(kg s)",zm)
        k = k + 1

      case ('rtp2_pd')
        irtp2_pd = k
        call stat_assign( irtp2_pd, "rtp2_pd", & 
             "rtp2 positive definite adjustment", "m^2/s^2", zm )
        k = k + 1

      case ('thlp2_bt')
        ithlp2_bt = k
        call stat_assign(ithlp2_bt,"thlp2_bt", & 
             "thlp2 budget","(K^2)/s",zm)
        k = k + 1
      case ('thlp2_ma')
        ithlp2_ma = k
        call stat_assign(ithlp2_ma,"thlp2_ma", & 
             "thlp2 mean advection","(K^2)/s",zm)
        k = k + 1
      case ('thlp2_ta')
        ithlp2_ta = k
        call stat_assign(ithlp2_ta,"thlp2_ta", & 
             "thlp2 turbulent advection","(K^2)/s",zm)
        k = k + 1
      case ('thlp2_tp')
        ithlp2_tp = k
        call stat_assign(ithlp2_tp,"thlp2_tp", & 
             "thlp2 turbulent production","(K^2)/s",zm)
        k = k + 1
      case ('thlp2_dp1')
        ithlp2_dp1 = k
        call stat_assign(ithlp2_dp1,"thlp2_dp1", & 
             "thlp2 dissipation term 1","(K^2)/s",zm)
        k = k + 1
      case ('thlp2_dp2')
        ithlp2_dp2 = k
        call stat_assign(ithlp2_dp2,"thlp2_dp2", & 
             "thlp2 dissipation term 2","(K^2)/s",zm)
        k = k + 1
      case ('thlp2_cl')
        ithlp2_cl = k
        call stat_assign(ithlp2_cl,"thlp2_cl", & 
             "thlp2 clipping term","(K^2)/s",zm)
        k = k + 1

      case ('thlp2_pd')
        ithlp2_pd = k
        call stat_assign( ithlp2_pd, "thlp2_pd", & 
             "thlp2 positive definite adjustment", "m^2/s^2", zm )
        k = k + 1

      case ('rtpthlp_bt')
        irtpthlp_bt = k
        call stat_assign(irtpthlp_bt,"rtpthlp_bt", & 
             "rtpthlp budget","(kg K)/(kg s)",zm)
        k = k + 1
      case ('rtpthlp_ma')
        irtpthlp_ma = k
        call stat_assign(irtpthlp_ma,"rtpthlp_ma", & 
             "rtpthlp mean advection","(kg K)/(kg s)",zm)
        k = k + 1
      case ('rtpthlp_ta')
        irtpthlp_ta = k
        call stat_assign(irtpthlp_ta,"rtpthlp_ta", & 
             "rtpthlp turbulent advection","(kg K)/(kg s)",zm)
        k = k + 1
      case ('rtpthlp_tp1')
        irtpthlp_tp1 = k
        call stat_assign(irtpthlp_tp1,"rtpthlp_tp1", & 
             "rtpthlp turbulent production 1","(kg K)/(kg s)",zm)
        k = k + 1
      case ('rtpthlp_tp2')
        irtpthlp_tp2 = k
        call stat_assign(irtpthlp_tp2,"rtpthlp_tp2", & 
             "rtpthlp turbulent production 2","(kg K)/(kg s)",zm)
        k = k + 1
      case ('rtpthlp_dp1')
        irtpthlp_dp1 = k
        call stat_assign(irtpthlp_dp1,"rtpthlp_dp1", & 
             "rtpthlp dissipation term 1","(kg K)/(kg s)",zm)
        k = k + 1
      case ('rtpthlp_dp2')
        irtpthlp_dp2 = k
        call stat_assign(irtpthlp_dp2,"rtpthlp_dp2", & 
             "rtpthlp dissipation term 2","(kg K)/(kg s)",zm)
        k = k + 1
      case ('rtpthlp_cl')
        irtpthlp_cl = k
        call stat_assign(irtpthlp_cl,"rtpthlp_cl", & 
             "rtpthlp clipping term","(kg K)/(kg s)",zm)
        k = k + 1

      case ('up2')
        iup2 = k
        call stat_assign(iup2,"up2", & 
             "up2","m^2/s^2",zm)
        k = k + 1

      case ('vp2')
        ivp2 = k
        call stat_assign(ivp2,"vp2", & 
             "vp2","m^2/s^2",zm)
        k = k + 1

      case ('up2_bt')
        iup2_bt = k
        call stat_assign(iup2_bt,"up2_bt", & 
             "up2 time tendency","m^2/s^2",zm)
        k = k + 1

      case ('up2_ma')
        iup2_ma = k
        call stat_assign(iup2_ma,"up2_ma", & 
             "up2 mean advection","m^2/s^2",zm)
        k = k + 1

      case ('up2_ta')
        iup2_ta = k
        call stat_assign(iup2_ta,"up2_ta", & 
             "up2 turbulent advection","m^2/s^2",zm)
        k = k + 1

      case ('up2_tp')
        iup2_tp = k
        call stat_assign(iup2_tp,"up2_tp", & 
             "up2 turbulent production","m^2/s^2",zm)
        k = k + 1

      case ('up2_dp1')
        iup2_dp1 = k
        call stat_assign(iup2_dp1,"up2_dp1", & 
             "up2 dissipation term 1","m^2/s^2",zm)
        k = k + 1

      case ('up2_dp2')
        iup2_dp2 = k
        call stat_assign(iup2_dp2,"up2_dp2", & 
             "up2 dissipation term 2","m^2/s^2",zm)
        k = k + 1

      case ('up2_pr1')
        iup2_pr1 = k
        call stat_assign(iup2_pr1,"up2_pr1", & 
             "up2 pressure term 1","m^2/s^2",zm)
        k = k + 1

      case ('up2_pr2')
        iup2_pr2 = k
        call stat_assign(iup2_pr2,"up2_pr2", & 
             "up2 pressure term 2","m^2/s^2",zm)
        k = k + 1

      case ('up2_cl')
        iup2_cl = k
        call stat_assign(iup2_cl,"up2_cl", & 
             "up2 clipping","m^2/s^2",zm)
        k = k + 1

      case ('up2_pd')
        iup2_pd = k
        call stat_assign( iup2_pd, "up2_pd", & 
             "up2 positive definite adjustment", "m^2/s^2", zm )
        k = k + 1

      case ('vp2_bt')
        ivp2_bt = k
        call stat_assign(ivp2_bt,"vp2_bt", & 
             "vp2 time tendency","m^2/s^2",zm)
        k = k + 1

      case ('vp2_ma')
        ivp2_ma = k
        call stat_assign(ivp2_ma,"vp2_ma", & 
             "vp2 mean advection","m^2/s^2",zm)
        k = k + 1

      case ('vp2_ta')
        ivp2_ta = k
        call stat_assign(ivp2_ta,"vp2_ta", & 
             "vp2 turbulent advection","m^2/s^2",zm)
        k = k + 1

      case ('vp2_tp')
        ivp2_tp = k
        call stat_assign(ivp2_tp,"vp2_tp", & 
             "vp2 turbulent production","m^2/s^2",zm)
        k = k + 1

      case ('vp2_dp1')
        ivp2_dp1 = k
        call stat_assign(ivp2_dp1,"vp2_dp1", & 
             "vp2 dissipation term 1","m^2/s^2",zm)
        k = k + 1

      case ('vp2_dp2')
        ivp2_dp2 = k
        call stat_assign(ivp2_dp2,"vp2_dp2", & 
             "vp2 dissipation term 2","m^2/s^2",zm)
        k = k + 1

      case ('vp2_pr1')
        ivp2_pr1 = k
        call stat_assign(ivp2_pr1,"vp2_pr1", & 
             "vp2 pressure term 1","m^2/s^2",zm)
        k = k + 1

      case ('vp2_pr2')
        ivp2_pr2 = k
        call stat_assign(ivp2_pr2,"vp2_pr2", & 
             "vp2 pressure term 2","m^2/s^2",zm)
        k = k + 1

      case ('vp2_cl')
        ivp2_cl = k
        call stat_assign(ivp2_cl,"vp2_cl", & 
             "vp2 clipping","m^2/s^2",zm)
        k = k + 1

      case ('vp2_pd')
        ivp2_pd = k
        call stat_assign( ivp2_pd, "vp2_pd", & 
             "vp2 positive definite adjustment", "m^2/s^2", zm )
        k = k + 1

      case ('wpthlp_enter_mfl')
        iwpthlp_enter_mfl = k
        call stat_assign( iwpthlp_enter_mfl, "wpthlp_in_mfl", & 
             "Wpthlp entering flux limiter ((m K)/s)", "(m K)/s", zm )
        k = k + 1

      case ('wpthlp_exit_mfl')
        iwpthlp_exit_mfl = k
        call stat_assign( iwpthlp_exit_mfl, "wpthlp_out_mfl", & 
             "Wpthlp exiting flux limiter ((m K)/s)", "(m K)/s", zm )
        k = k + 1

      case ('wpthlp_mfl_lower_lim')
        iwpthlp_mfl_lower_lim = k
        call stat_assign( iwpthlp_mfl_lower_lim, "wpthlp_mfl_min", & 
             "Minimum allowable wpthlp ((m K)/s)", "(m K)/s", zm )
        k = k + 1

      case ('wpthlp_mfl_upper_lim')
        iwpthlp_mfl_upper_lim = k
        call stat_assign( iwpthlp_mfl_upper_lim, "wpthlp_mfl_max", & 
             "Maximum allowable wpthlp ((m K)/s)", "(m K)/s", zm )
        k = k + 1

      case ('wprtp_mfl_lower_lim')
        iwprtp_mfl_lower_lim = k
        call stat_assign( iwprtp_mfl_lower_lim, "wprtp_mfl_min", & 
             "Minimum allowable wprtp ((m kg)/(s kg))", "(m kg)/(s kg)", zm )
        k = k + 1

      case ('wprtp_mfl_upper_lim')
        iwprtp_mfl_upper_lim = k
        call stat_assign( iwprtp_mfl_upper_lim, "wprtp_mfl_max", & 
             "Maximum allowable wprtp ((m kg)/(s kg))", "(m kg)/(s kg)", zm )
        k = k + 1

      case ('wprtp_enter_mfl')
        iwprtp_enter_mfl = k
        call stat_assign( iwprtp_enter_mfl, "wprtp_in_mfl", & 
             "Wprtp entering flux limiter ((m kg)/(s kg))", "(m kg)/(s kg)", zm )
        k = k + 1

      case ('wprtp_exit_mfl')
        iwprtp_exit_mfl = k
        call stat_assign( iwprtp_exit_mfl, "wprtp_out_mfl", & 
             "Wprtp exiting flux limiter ((m kg)/(s kg))", "(m kg)/(s kg)", zm )
        k = k + 1        

      case default
        l_found = .false.

        j = 1

        do while( j <= sclr_dim .and. .not. l_found )
          write( sclr_idx, * ) j
          sclr_idx = adjustl(sclr_idx)

          if( trim(vars_zm(i)) == 'sclr'//trim(sclr_idx)//'prtp'.and. .not. l_found ) then
            isclrprtp(j) = k

            call stat_assign(isclrprtp(j),"sclr"//trim(sclr_idx)//"prtp", & 
               "scalar("//trim(sclr_idx)//")'rt'","unknown",zm)
            k = k + 1
            l_found = .true.
          end if
          if( trim(vars_zm(i)) == 'sclr'//trim(sclr_idx)//'p2'.and. .not. l_found ) then
            isclrp2(j) = k
            call stat_assign(isclrp2(j) ,"sclr"//trim(sclr_idx)//"p2", & 
               "scalar("//trim(sclr_idx)//")'^2'","unknown",zm)
            k = k + 1
            l_found = .true.
          end if
          if( trim(vars_zm(i)) == 'sclr'//trim(sclr_idx)//'pthvp'.and. .not. l_found ) then
            isclrpthvp(j) = k
            call stat_assign(isclrpthvp(j),"sclr"//trim(sclr_idx)//"pthvp", & 
               "scalar("//trim(sclr_idx)//")'th_v'","unknown",zm)
            k = k + 1
            l_found = .true.
          end if
          if( trim(vars_zm(i)) == 'sclr'//trim(sclr_idx)//'pthlp'.and. .not. l_found ) then
            isclrpthlp(j) = k

            call stat_assign(isclrpthlp(j),"sclr"//trim(sclr_idx)//"pthlp", & 
               "scalar("//trim(sclr_idx)//")'th_l'","unknown",zm)
            k = k + 1
            l_found = .true.
          end if
          if( trim(vars_zm(i)) == 'sclr'//trim(sclr_idx)//'prcp'.and. .not. l_found ) then

            isclrprcp(j) = k

            call stat_assign(isclrprcp(j),"sclr"//trim(sclr_idx)//"prcp", & 
               "scalar("//trim(sclr_idx)//")'rc'","unknown",zm)
            k = k + 1
            l_found = .true.
          end if
          if( trim(vars_zm(i)) == 'wpsclr'//trim(sclr_idx)//'p'.and. .not. l_found ) then
            iwpsclrp(j) = k

            call stat_assign(iwpsclrp(j),"wpsclr"//trim(sclr_idx)//"p", & 
               "'w'scalar("//trim(sclr_idx)//")","unknown",zm)
            k = k + 1
            l_found = .true.
          end if
          if( trim(vars_zm(i)) == 'wpsclr'//trim(sclr_idx)//'p2'.and. .not. l_found ) then

            iwpsclrp2(j) = k

            call stat_assign(iwpsclrp2(j),"wpsclr"//trim(sclr_idx)//"p2", & 
               "'w'scalar("//trim(sclr_idx)//")'^2'","unknown",zm)
            k = k + 1
            l_found = .true.
          end if
          if( trim(vars_zm(i)) == 'wp2sclr'//trim(sclr_idx)//'p'.and. .not. l_found ) then

            iwp2sclrp(j) = k

            call stat_assign(iwp2sclrp(j) ,"wp2sclr"//trim(sclr_idx)//"p", & 
                 "'w'^2 scalar("//trim(sclr_idx)//")","unknown",zm)
            k = k + 1
            l_found = .true.
          end if
          if( trim(vars_zm(i)) == 'wpsclr'//trim(sclr_idx)//'prtp'.and. .not. l_found ) then
            iwpsclrprtp(j) = k

            call stat_assign( iwpsclrprtp(j),"wpsclr"//trim(sclr_idx)//"prtp", & 
               "'w' scalar("//trim(sclr_idx)//")'rt'","unknown",zm )
            k = k + 1
            l_found = .true.
          end if
          if( trim(vars_zm(i)) == 'wpsclr'//trim(sclr_idx)//'pthlp'.and. .not. l_found ) then
            iwpsclrpthlp(j) = k

            call stat_assign(iwpsclrpthlp(j),"wpsclr"//trim(sclr_idx)//"pthlp", & 
               "'w' scalar("//trim(sclr_idx)//")'th_l'","unknown",zm)
            k = k + 1
            l_found = .true.
          end if
          j = j + 1
        end do

        j = 1

        do while( j <= edsclr_dim .and. .not. l_found )

          write( sclr_idx, * ) j
          sclr_idx = adjustl(sclr_idx)

          if( trim(vars_zm(i)) == 'wpedsclr'//trim(sclr_idx)//'p'.and. .not. l_found ) then
            iwpedsclrp(j) = k

            call stat_assign(iwpedsclrp(j),"wpedsclr"//trim(sclr_idx)//"p", & 
               "eddy scalar("//trim(sclr_idx)//")'w'","unknown",zm)
            k = k + 1
            l_found = .true.
          end if

          j = j + 1

        end do

        if( .not. l_found ) then
          write(fstderr,*) 'Error:  unrecognized variable in vars_zm:  ',  trim(vars_zm(i))
          l_error = .true.  ! This will stop the run.
        end if
      end select

    end do

!   Non-interative diagnostics (zm)
!   iwp4, ircp2

!   if ( .not. clubb_at_least_debug_level( 1 ) ) then
!     if ( iwp4 + ircp2 + ishear > 0 ) then
!       write(fstderr,'(a)') &
!         "Warning: at debug level 0.  Non-interactive diagnostics will not be computed, "
!       write(fstderr,'(a)') "but some appear in the stats_zm namelist variable."
!     end if
!   end if

    return
  end subroutine stats_init_zm

end module stats_zm
