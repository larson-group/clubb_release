!-----------------------------------------------------------------------
! $Id$
!===============================================================================
module pdf_parameter_module

  ! Description:
  ! This module defines the derived type pdf_parameter.

  ! References:
  !   None
  !-----------------------------------------------------------------------

  use clubb_precision, only: &
      core_rknd

  implicit none

  private ! Default scope

  public :: pdf_parameter,                 & ! Variable Type(s)
            implicit_coefs_terms,          &
            init_pdf_params,               & ! Procedure(s)
            zero_pdf_params_api,               &
            print_pdf_params,              &
            init_pdf_implicit_coefs_terms_api, &
            zero_pdf_implicit_coefs_terms_api

  ! CLUBB's PDF parameters.
  type pdf_parameter

    integer :: & 
      ngrdcol, &  ! Dimensions of variables
      nz

    real( kind = core_rknd ), dimension(:,:), allocatable :: &
      w_1,             & ! Mean of w (1st PDF component)                   [m/s]
      w_2,             & ! Mean of w (2nd PDF component)                   [m/s]
      varnce_w_1,      & ! Variance of w (1st PDF component)           [m^2/s^2]
      varnce_w_2,      & ! Variance of w (2nd PDF component)           [m^2/s^2]
      rt_1,            & ! Mean of r_t (1st PDF component)               [kg/kg]
      rt_2,            & ! Mean of r_t (2nd PDF component)               [kg/kg]
      varnce_rt_1,     & ! Variance of r_t (1st PDF component)       [kg^2/kg^2]
      varnce_rt_2,     & ! Variance of r_t (2nd PDF component)       [kg^2/kg^2]
      thl_1,           & ! Mean of th_l (1st PDF component)                  [K]
      thl_2,           & ! Mean of th_l (2nd PDF component)                  [K]
      varnce_thl_1,    & ! Variance of th_l (1st PDF component)            [K^2]
      varnce_thl_2,    & ! Variance of th_l (2nd PDF component)            [K^2]
      corr_w_rt_1,     & ! Correlation of w and r_t (1st PDF component)      [-]
      corr_w_rt_2,     & ! Correlation of w and r_t (2nd PDF component)      [-]
      corr_w_thl_1,    & ! Correlation of w and th_l (1st PDF component)     [-]
      corr_w_thl_2,    & ! Correlation of w and th_l (2nd PDF component)     [-]
      corr_rt_thl_1,   & ! Correlation of r_t and th_l (1st PDF component)   [-]
      corr_rt_thl_2,   & ! Correlation of r_t and th_l (2nd PDF component)   [-]
      alpha_thl,       & ! Factor relating to normalized variance for th_l   [-]
      alpha_rt,        & ! Factor relating to normalized variance for r_t    [-]
      crt_1,           & ! r_t coef. in chi/eta eqns. (1st PDF comp.)        [-]
      crt_2,           & ! r_t coef. in chi/eta eqns. (2nd PDF comp.)        [-]
      cthl_1,          & ! th_l coef.: chi/eta eqns. (1st PDF comp.) [(kg/kg)/K]
      cthl_2,          & ! th_l coef.: chi/eta eqns. (2nd PDF comp.) [(kg/kg)/K]
      chi_1,           & ! Mean of chi (old s) (1st PDF component)       [kg/kg]
      chi_2,           & ! Mean of chi (old s) (2nd PDF component)       [kg/kg]
      stdev_chi_1,     & ! Standard deviation of chi (1st PDF component) [kg/kg]
      stdev_chi_2,     & ! Standard deviation of chi (2nd PDF component) [kg/kg]
      stdev_eta_1,     & ! Standard dev. of eta (old t) (1st PDF comp.)  [kg/kg]
      stdev_eta_2,     & ! Standard dev. of eta (old t) (2nd PDF comp.)  [kg/kg]
      covar_chi_eta_1, & ! Covariance of chi and eta (1st PDF comp.) [kg^2/kg^2]
      covar_chi_eta_2, & ! Covariance of chi and eta (2nd PDF comp.) [kg^2/kg^2]
      corr_w_chi_1,    & ! Correlation of w and chi (1st PDF component)      [-]
      corr_w_chi_2,    & ! Correlation of w and chi (2nd PDF component)      [-]
      corr_w_eta_1,    & ! Correlation of w and eta (1st PDF component)      [-]
      corr_w_eta_2,    & ! Correlation of w and eta (2nd PDF component)      [-]
      corr_chi_eta_1,  & ! Correlation of chi and eta (1st PDF component)    [-]
      corr_chi_eta_2,  & ! Correlation of chi and eta (2nd PDF component)    [-]
      rsatl_1,         & ! Saturation mixing ratio r_sat(mu_Tl_1,p)      [kg/kg]
      rsatl_2,         & ! Saturation mixing ratio r_sat(mu_Tl_2,p)      [kg/kg]
      rc_1,            & ! Mean of r_c (1st PDF component)               [kg/kg]
      rc_2,            & ! Mean of r_c (2nd PDF component)               [kg/kg]
      cloud_frac_1,    & ! Cloud fraction (1st PDF component)                [-]
      cloud_frac_2,    & ! Cloud fraction (2nd PDF component)                [-]
      mixt_frac,       & ! Weight of 1st PDF component (Sk_w dependent)      [-]
      ice_supersat_frac_1, & ! Ice supersaturation fraction (1st PDF comp.)  [-]
      ice_supersat_frac_2    ! Ice supersaturation fraction (2nd PDF comp.)  [-]

  end type pdf_parameter

  ! The implicit coefficients, semi-implicit coefficients and terms, and
  ! explicit terms for turbulent advection of turbulent fields are calculated
  ! from the PDF and the resulting PDF parameters.
  type implicit_coefs_terms

    integer :: & 
      ngrdcol, &  ! Dimensions of variables
      nz, &
      sclr_dim

    real( kind = core_rknd ), dimension(:,:), allocatable :: &
      coef_wp4_implicit    ! <w'^4> = coef_wp4_implicit * <w'^2>^2       [-]

    ! <w'^2 rt'> = coef_wp2rtp_implicit * <w'rt'> + term_wp2rtp_explicit
    real( kind = core_rknd ), dimension(:,:), allocatable :: &
      coef_wp2rtp_implicit, & ! Coefficient that is multiplied by <w'rt'>  [m/s]
      term_wp2rtp_explicit    ! Term that is on the RHS          [m^2/s^2 kg/kg]

    ! <w'^2 thl'> = coef_wp2thlp_implicit * <w'thl'> + term_wp2thlp_explicit
    real( kind = core_rknd ), dimension(:,:), allocatable :: &
      coef_wp2thlp_implicit, & ! Coef. that is multiplied by <w'thl'>      [m/s]
      term_wp2thlp_explicit    ! Term that is on the RHS             [m^2/s^2 K]

    ! <w'^2 u'> = coef_wp2up_implicit * <u'w'> + term_wp2up_explicit
    real( kind = core_rknd ), dimension(:,:), allocatable :: &
      coef_wp2up_implicit, & ! Coefficient that is multiplied by <u'w'>    [m/s]
      term_wp2up_explicit    ! Term that is on the RHS                 [m^3/s^3]

    ! <w'^2 v'> = coef_wp2vp_implicit * <v'w'> + term_wp2vp_explicit
    real( kind = core_rknd ), dimension(:,:), allocatable :: &
      coef_wp2vp_implicit, & ! Coefficient that is multiplied by <v'w'>    [m/s]
      term_wp2vp_explicit    ! Term that is on the RHS                 [m^3/s^3]

    ! <w'rt'^2> = coef_wprtp2_implicit * <rt'^2> + term_wprtp2_explicit
    real( kind = core_rknd ), dimension(:,:), allocatable :: &
      coef_wprtp2_implicit, & ! Coefficient that is multiplied by <rt'^2>  [m/s]
      term_wprtp2_explicit    ! Term that is on the RHS          [m/s kg^2/kg^2]

    ! <w'thl'^2> = coef_wpthlp2_implicit * <thl'^2> + term_wpthlp2_explicit
    real( kind = core_rknd ), dimension(:,:), allocatable :: &
      coef_wpthlp2_implicit, & ! Coef. that is multiplied by <thl'^2>      [m/s]
      term_wpthlp2_explicit    ! Term that is on the RHS               [m/s K^2]

    ! <w'rt'thl'> = coef_wprtpthlp_implicit*<rt'thl'> + term_wprtpthlp_explicit
    real( kind = core_rknd ), dimension(:,:), allocatable :: &
      coef_wprtpthlp_implicit, & ! Coef. that is multiplied by <rt'thl'>   [m/s]
      term_wprtpthlp_explicit    ! Term that is on the RHS         [m/s(kg/kg)K]

    ! <w'u'^2> = coef_wpup2_implicit * <u'^2> + term_wpup2_explicit
    real( kind = core_rknd ), dimension(:,:), allocatable :: &
      coef_wpup2_implicit, & ! Coefficient that is multiplied by <u'^2>    [m/s]
      term_wpup2_explicit    ! Term that is on the RHS                 [m^3/s^3]

    ! <w'v'^2> = coef_wpvp2_implicit * <v'^2> + term_wpvp2_explicit
    real( kind = core_rknd ), dimension(:,:), allocatable :: &
      coef_wpvp2_implicit, & ! Coefficient that is multiplied by <v'^2>    [m/s]
      term_wpvp2_explicit    ! Term that is on the RHS                 [m^3/s^3]

    ! <w'^2 sclr'> = coef_wp2sclrp_implicit * <w'sclr'> + term_wp2sclrp_explicit
    real( kind = core_rknd ), dimension(:,:,:), allocatable :: &
      coef_wp2sclrp_implicit, & ! Coef. that is multiplied by <w'sclr'>    [m/s]
      term_wp2sclrp_explicit    ! Term that is on the RHS   [m^2/s^2 (un. vary)]

    ! <w'sclr'^2> = coef_wpsclrp2_implicit * <sclr'^2> + term_wpsclrp2_explicit
    real( kind = core_rknd ), dimension(:,:,:), allocatable :: &
      coef_wpsclrp2_implicit, & ! Coef. that is multiplied by <sclr'^2>    [m/s]
      term_wpsclrp2_explicit    ! Term that is on the RHS    [m/s(units vary)^2]

    ! <w'rt'sclr'> = coef_wprtpsclrp_implicit * <sclr'rt'>
    !                + term_wprtpsclrp_explicit
    real( kind = core_rknd ), dimension(:,:,:), allocatable :: &
      coef_wprtpsclrp_implicit, & ! Coef. that is multiplied by <sclr'rt'> [m/s]
      term_wprtpsclrp_explicit    ! Term that is on the RHS [m/s(kg/kg)(un. v.)]

    ! <w'thl'sclr'> = coef_wpthlpsclrp_implicit * <sclr'thl'>
    !                 + term_wpthlpsclrp_explicit
    real( kind = core_rknd ), dimension(:,:,:), allocatable :: &
      coef_wpthlpsclrp_implicit, & ! Coef. that is mult. by <sclr'thl'>    [m/s]
      term_wpthlpsclrp_explicit    ! Term that is on the RHS  [(m/s)K(un. vary)]

  end type implicit_coefs_terms

  contains

  !=============================================================================
  subroutine init_pdf_params( nz, ngrdcol, &
                              pdf_params )

    ! Description:
    ! Initializes all PDF parameters in the variable type pdf_parameter.

    ! References:
    !--------------------------------------------------------------------

    implicit none

    !--------------------------- Input Variables ---------------------------
    integer, intent(in) :: &
      nz,   & ! Number of vertical grid levels    [-]
      ngrdcol ! Number of grid columns            [-]

    !--------------------------- Output Variables ---------------------------
    type(pdf_parameter), intent(out) :: &
      pdf_params    ! PDF parameters            [units vary]

    allocate( pdf_params%w_1(ngrdcol,nz), &
              pdf_params%w_2(ngrdcol,nz), &
              pdf_params%varnce_w_1(ngrdcol,nz), &
              pdf_params%varnce_w_2(ngrdcol,nz), &
              pdf_params%rt_1(ngrdcol,nz), &
              pdf_params%rt_2(ngrdcol,nz), &
              pdf_params%varnce_rt_1(ngrdcol,nz), &
              pdf_params%varnce_rt_2(ngrdcol,nz), &
              pdf_params%thl_1(ngrdcol,nz), &
              pdf_params%thl_2(ngrdcol,nz), &
              pdf_params%varnce_thl_1(ngrdcol,nz), &
              pdf_params%varnce_thl_2(ngrdcol,nz), &
              pdf_params%corr_w_rt_1(ngrdcol,nz), &
              pdf_params%corr_w_rt_2(ngrdcol,nz), &
              pdf_params%corr_w_thl_1(ngrdcol,nz), &
              pdf_params%corr_w_thl_2(ngrdcol,nz), &
              pdf_params%corr_rt_thl_1(ngrdcol,nz), &
              pdf_params%corr_rt_thl_2(ngrdcol,nz), &
              pdf_params%alpha_thl(ngrdcol,nz), &
              pdf_params%alpha_rt(ngrdcol,nz), &
              pdf_params%crt_1(ngrdcol,nz), &
              pdf_params%crt_2(ngrdcol,nz), &
              pdf_params%cthl_1(ngrdcol,nz), &
              pdf_params%cthl_2(ngrdcol,nz), &
              pdf_params%chi_1(ngrdcol,nz), &
              pdf_params%chi_2(ngrdcol,nz), &
              pdf_params%stdev_chi_1(ngrdcol,nz), &
              pdf_params%stdev_chi_2(ngrdcol,nz), &
              pdf_params%stdev_eta_1(ngrdcol,nz), &
              pdf_params%stdev_eta_2(ngrdcol,nz), &
              pdf_params%covar_chi_eta_1(ngrdcol,nz), &
              pdf_params%covar_chi_eta_2(ngrdcol,nz), &
              pdf_params%corr_w_chi_1(ngrdcol,nz), & 
              pdf_params%corr_w_chi_2(ngrdcol,nz), &
              pdf_params%corr_w_eta_1(ngrdcol,nz), & 
              pdf_params%corr_w_eta_2(ngrdcol,nz), & 
              pdf_params%corr_chi_eta_1(ngrdcol,nz), & 
              pdf_params%corr_chi_eta_2(ngrdcol,nz), &
              pdf_params%rsatl_1(ngrdcol,nz), &
              pdf_params%rsatl_2(ngrdcol,nz), &
              pdf_params%rc_1(ngrdcol,nz), &
              pdf_params%rc_2(ngrdcol,nz), &
              pdf_params%cloud_frac_1(ngrdcol,nz), &
              pdf_params%cloud_frac_2(ngrdcol,nz), &
              pdf_params%mixt_frac(ngrdcol,nz), &
              pdf_params%ice_supersat_frac_1(ngrdcol,nz), &
              pdf_params%ice_supersat_frac_2(ngrdcol,nz) )

    ! Set metadata values 
    pdf_params%ngrdcol       = ngrdcol
    pdf_params%nz            = nz

    call zero_pdf_params_api( pdf_params )

    return

  end subroutine init_pdf_params

  !=============================================================================
  subroutine zero_pdf_params_api( pdf_params )

    ! Description:
    ! Zeros all PDF parameters in the variable type pdf_parameter.

    ! References:
    !--------------------------------------------------------------------

    use constants_clubb, only: &
      zero    ! Constant(s)

    implicit none

    !--------------------------- Input/Output Variables ---------------------------
    !----------------------------------- Inout Variable(s) -----------------------------------
    type(pdf_parameter), intent(inout) :: &
      pdf_params    ! PDF parameters            [units vary]

    !----------------------------------- Local Variables -----------------------------------
    integer :: i, k

    !----------------------------------- Begin Code -----------------------------------

    !$acc data copyin( pdf_params ) &
    !$acc      copy( pdf_params%w_1, pdf_params%w_2, pdf_params%varnce_w_1, &
    !$acc            pdf_params%varnce_w_2, pdf_params%mixt_frac, &
    !$acc            pdf_params%rt_1, pdf_params%rt_2, &
    !$acc            pdf_params%varnce_rt_1, pdf_params%varnce_rt_2, pdf_params%thl_1, &
    !$acc            pdf_params%thl_2, pdf_params%varnce_thl_1, pdf_params%varnce_thl_2, &
    !$acc            pdf_params%corr_w_rt_1, pdf_params%corr_w_rt_2, pdf_params%corr_w_thl_1, &
    !$acc            pdf_params%corr_w_thl_2, pdf_params%corr_rt_thl_1, pdf_params%corr_rt_thl_2, &
    !$acc            pdf_params%alpha_thl, pdf_params%alpha_rt, &
    !$acc            pdf_params%crt_1, pdf_params%crt_2, &
    !$acc            pdf_params%cthl_1, pdf_params%cthl_2, pdf_params%chi_1, pdf_params%chi_2, &
    !$acc            pdf_params%stdev_chi_1, pdf_params%stdev_chi_2, pdf_params%stdev_eta_1, &
    !$acc            pdf_params%stdev_eta_2, pdf_params%covar_chi_eta_1, &
    !$acc            pdf_params%covar_chi_eta_2, &
    !$acc            pdf_params%corr_w_chi_1, pdf_params%corr_w_chi_2, pdf_params%corr_w_eta_1, &
    !$acc            pdf_params%corr_w_eta_2, pdf_params%corr_chi_eta_1, &
    !$acc            pdf_params%corr_chi_eta_2, &
    !$acc            pdf_params%rsatl_1, pdf_params%rsatl_2, pdf_params%rc_1, pdf_params%rc_2, &
    !$acc            pdf_params%cloud_frac_1, pdf_params%cloud_frac_2, &
    !$acc            pdf_params%ice_supersat_frac_1, pdf_params%ice_supersat_frac_2 )

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, pdf_params%nz
      do i = 1, pdf_params%ngrdcol
        pdf_params%w_1(i,k)                 = zero
        pdf_params%w_2(i,k)                 = zero
        pdf_params%varnce_w_1(i,k)          = zero
        pdf_params%varnce_w_2(i,k)          = zero
        pdf_params%rt_1(i,k)                = zero
        pdf_params%rt_2(i,k)                = zero
        pdf_params%varnce_rt_1(i,k)         = zero
        pdf_params%varnce_rt_2(i,k)         = zero
        pdf_params%thl_1(i,k)               = zero
        pdf_params%thl_2(i,k)               = zero
        pdf_params%varnce_thl_1(i,k)        = zero
        pdf_params%varnce_thl_2(i,k)        = zero
        pdf_params%corr_w_rt_1(i,k)         = zero
        pdf_params%corr_w_rt_2(i,k)         = zero
        pdf_params%corr_w_thl_1(i,k)        = zero
        pdf_params%corr_w_thl_2(i,k)        = zero
        pdf_params%corr_rt_thl_1(i,k)       = zero
        pdf_params%corr_rt_thl_2(i,k)       = zero
        pdf_params%alpha_thl(i,k)           = zero
        pdf_params%alpha_rt(i,k)            = zero
        pdf_params%crt_1(i,k)               = zero
        pdf_params%crt_2(i,k)               = zero
        pdf_params%cthl_1(i,k)              = zero
        pdf_params%cthl_2(i,k)              = zero
        pdf_params%chi_1(i,k)               = zero
        pdf_params%chi_2(i,k)               = zero
        pdf_params%stdev_chi_1(i,k)         = zero
        pdf_params%stdev_chi_2(i,k)         = zero
        pdf_params%stdev_eta_1(i,k)         = zero
        pdf_params%stdev_eta_2(i,k)         = zero
        pdf_params%covar_chi_eta_1(i,k)     = zero
        pdf_params%covar_chi_eta_2(i,k)     = zero
        pdf_params%corr_w_chi_1(i,k)        = zero 
        pdf_params%corr_w_chi_2(i,k)        = zero 
        pdf_params%corr_w_eta_1(i,k)        = zero 
        pdf_params%corr_w_eta_2(i,k)        = zero 
        pdf_params%corr_chi_eta_1(i,k)      = zero 
        pdf_params%corr_chi_eta_2(i,k)      = zero 
        pdf_params%rsatl_1(i,k)             = zero
        pdf_params%rsatl_2(i,k)             = zero
        pdf_params%rc_1(i,k)                = zero
        pdf_params%rc_2(i,k)                = zero
        pdf_params%cloud_frac_1(i,k)        = zero
        pdf_params%cloud_frac_2(i,k)        = zero
        pdf_params%mixt_frac(i,k)           = zero
        pdf_params%ice_supersat_frac_1(i,k) = zero
        pdf_params%ice_supersat_frac_2(i,k) = zero
      end do
    end do

    !$acc end data

    return

  end subroutine zero_pdf_params_api

  !=============================================================================
  subroutine init_pdf_implicit_coefs_terms_api( nz, ngrdcol, sclr_dim, &
                                                pdf_implicit_coefs_terms )

    ! Description:
    ! Initializes all PDF implicit coefficients and explicit terms in the
    ! variable type implicit_coefs_terms.

    ! References:
    !--------------------------------------------------------------------

    implicit none

    !--------------------------- Input Variables ---------------------------
    integer, intent(in) :: &
      nz,       & ! Number of vertical grid levels    [-]
      ngrdcol,  & ! Number of grid columns            [-]
      sclr_dim    ! Number of scalar variables        [-]

    !--------------------------- Output Variables ---------------------------
    type(implicit_coefs_terms), intent(out) :: &
      pdf_implicit_coefs_terms    ! Implicit coefs / explicit terms [units vary]


    ! Allocate pdf_implicit_coefs_terms
    allocate( pdf_implicit_coefs_terms%coef_wp4_implicit(1:ngrdcol,1:nz), &
              pdf_implicit_coefs_terms%coef_wp2rtp_implicit(1:ngrdcol,1:nz), &
              pdf_implicit_coefs_terms%term_wp2rtp_explicit(1:ngrdcol,1:nz), &
              pdf_implicit_coefs_terms%coef_wp2thlp_implicit(1:ngrdcol,1:nz), &
              pdf_implicit_coefs_terms%term_wp2thlp_explicit(1:ngrdcol,1:nz), &
              pdf_implicit_coefs_terms%coef_wp2up_implicit(1:ngrdcol,1:nz), &
              pdf_implicit_coefs_terms%term_wp2up_explicit(1:ngrdcol,1:nz), &
              pdf_implicit_coefs_terms%coef_wp2vp_implicit(1:ngrdcol,1:nz), &
              pdf_implicit_coefs_terms%term_wp2vp_explicit(1:ngrdcol,1:nz), &
              pdf_implicit_coefs_terms%coef_wprtp2_implicit(1:ngrdcol,1:nz), &
              pdf_implicit_coefs_terms%term_wprtp2_explicit(1:ngrdcol,1:nz), &
              pdf_implicit_coefs_terms%coef_wpthlp2_implicit(1:ngrdcol,1:nz), &
              pdf_implicit_coefs_terms%term_wpthlp2_explicit(1:ngrdcol,1:nz), &
              pdf_implicit_coefs_terms%coef_wprtpthlp_implicit(1:ngrdcol,1:nz), &
              pdf_implicit_coefs_terms%term_wprtpthlp_explicit(1:ngrdcol,1:nz), &
              pdf_implicit_coefs_terms%coef_wpup2_implicit(1:ngrdcol,1:nz), &
              pdf_implicit_coefs_terms%term_wpup2_explicit(1:ngrdcol,1:nz), &
              pdf_implicit_coefs_terms%coef_wpvp2_implicit(1:ngrdcol,1:nz), &
              pdf_implicit_coefs_terms%term_wpvp2_explicit(1:ngrdcol,1:nz) )

    if ( sclr_dim > 0 ) then
       allocate( &
          pdf_implicit_coefs_terms%coef_wp2sclrp_implicit(1:ngrdcol,1:nz,1:sclr_dim), &
          pdf_implicit_coefs_terms%term_wp2sclrp_explicit(1:ngrdcol,1:nz,1:sclr_dim), &
          pdf_implicit_coefs_terms%coef_wpsclrp2_implicit(1:ngrdcol,1:nz,1:sclr_dim), &
          pdf_implicit_coefs_terms%term_wpsclrp2_explicit(1:ngrdcol,1:nz,1:sclr_dim), &
          pdf_implicit_coefs_terms%coef_wprtpsclrp_implicit(1:ngrdcol,1:nz,1:sclr_dim), &
          pdf_implicit_coefs_terms%term_wprtpsclrp_explicit(1:ngrdcol,1:nz,1:sclr_dim), &
          pdf_implicit_coefs_terms%coef_wpthlpsclrp_implicit(1:ngrdcol,1:nz,1:sclr_dim), &
          pdf_implicit_coefs_terms%term_wpthlpsclrp_explicit(1:ngrdcol,1:nz,1:sclr_dim) )
    endif ! sclr_dim > 0

    ! Set metadata values 
    pdf_implicit_coefs_terms%ngrdcol       = ngrdcol
    pdf_implicit_coefs_terms%nz            = nz
    pdf_implicit_coefs_terms%sclr_dim      = sclr_dim

    ! Zero pdf_implicit_coefs_terms
    call zero_pdf_implicit_coefs_terms_api( pdf_implicit_coefs_terms)

    return

  end subroutine init_pdf_implicit_coefs_terms_api

  subroutine zero_pdf_implicit_coefs_terms_api( pdf_implicit_coefs_terms )

    ! Description:
    ! Zeros all PDF implicit coefficients and explicit terms in the
    ! variable type implicit_coefs_terms.

    ! References:
    !--------------------------------------------------------------------

    use constants_clubb, only: &
        zero    ! Constant(s)

    implicit none

    !--------------------------- Input/Output Variables ---------------------------
    !----------------------------------- Inout Variable(s) -----------------------------------
    type(implicit_coefs_terms), intent(inout) :: &
      pdf_implicit_coefs_terms    ! Implicit coefs / explicit terms [units vary]

    !----------------------------------- Local Variables -----------------------------------
    integer :: i, k, sclr

    !----------------------------------- Begin Code -----------------------------------

    ! Initialize pdf_implicit_coefs_terms
    do k = 1, pdf_implicit_coefs_terms%nz
      do i = 1, pdf_implicit_coefs_terms%ngrdcol
        pdf_implicit_coefs_terms%coef_wp4_implicit(i,k) = zero
        pdf_implicit_coefs_terms%coef_wp2rtp_implicit(i,k) = zero
        pdf_implicit_coefs_terms%term_wp2rtp_explicit(i,k) = zero
        pdf_implicit_coefs_terms%coef_wp2thlp_implicit(i,k) = zero
        pdf_implicit_coefs_terms%term_wp2thlp_explicit(i,k) = zero
        pdf_implicit_coefs_terms%coef_wp2up_implicit(i,k) = zero
        pdf_implicit_coefs_terms%term_wp2up_explicit(i,k) = zero
        pdf_implicit_coefs_terms%coef_wp2vp_implicit(i,k) = zero
        pdf_implicit_coefs_terms%term_wp2vp_explicit(i,k) = zero
        pdf_implicit_coefs_terms%coef_wprtp2_implicit(i,k) = zero
        pdf_implicit_coefs_terms%term_wprtp2_explicit(i,k) = zero
        pdf_implicit_coefs_terms%coef_wpthlp2_implicit(i,k) = zero
        pdf_implicit_coefs_terms%term_wpthlp2_explicit(i,k) = zero
        pdf_implicit_coefs_terms%coef_wprtpthlp_implicit(i,k) = zero
        pdf_implicit_coefs_terms%term_wprtpthlp_explicit(i,k) = zero
        pdf_implicit_coefs_terms%coef_wpup2_implicit(i,k) = zero
        pdf_implicit_coefs_terms%term_wpup2_explicit(i,k) = zero
        pdf_implicit_coefs_terms%coef_wpvp2_implicit(i,k) = zero
        pdf_implicit_coefs_terms%term_wpvp2_explicit(i,k) = zero
      end do
    end do

    if ( pdf_implicit_coefs_terms%sclr_dim > 0 ) then
      do sclr = 1, pdf_implicit_coefs_terms%sclr_dim
        do k = 1, pdf_implicit_coefs_terms%nz
          do i = 1, pdf_implicit_coefs_terms%ngrdcol
            pdf_implicit_coefs_terms%coef_wp2sclrp_implicit(i,k,sclr) = zero
            pdf_implicit_coefs_terms%term_wp2sclrp_explicit(i,k,sclr) = zero
            pdf_implicit_coefs_terms%coef_wpsclrp2_implicit(i,k,sclr) = zero
            pdf_implicit_coefs_terms%term_wpsclrp2_explicit(i,k,sclr) = zero
            pdf_implicit_coefs_terms%coef_wprtpsclrp_implicit(i,k,sclr) = zero
            pdf_implicit_coefs_terms%term_wprtpsclrp_explicit(i,k,sclr) = zero
            pdf_implicit_coefs_terms%coef_wpthlpsclrp_implicit(i,k,sclr) = zero
            pdf_implicit_coefs_terms%term_wpthlpsclrp_explicit(i,k,sclr) = zero
          end do
        end do
      end do
    endif ! sclr_dim > 0

    return

  end subroutine zero_pdf_implicit_coefs_terms_api

  !=============================================================================
  subroutine print_pdf_params( pdf_params, ngrdcol )

    ! Description:
    ! Prints every value found in type pdf_parameter

    !-------------------------------------------------------------------

    use constants_clubb, only: &
        fstderr

    implicit none

    !--------------------------- Input Variables ---------------------------
    type(pdf_parameter), intent(in) :: &
      pdf_params    ! PDF parameters            [units vary]

    integer, intent(in) :: &
      ngrdcol    ! Number of horizontal grid columns

    ! Local Variable(s)
    integer :: i    ! Loop index


    do i = 1, ngrdcol

       write(fstderr,*) "PDF parameters for column ", i
       write(fstderr,*) "---------------------------------------------------"
       write(fstderr,*) "w_1(i,:) = ", pdf_params%w_1(i,:)
       write(fstderr,*) "w_2(i,:) = ", pdf_params%w_2(i,:)
       write(fstderr,*) "varnce_w_1(i,:) = ", pdf_params%varnce_w_1(i,:)
       write(fstderr,*) "varnce_w_2(i,:) = ", pdf_params%varnce_w_2(i,:)
       write(fstderr,*) "rt_1(i,:) = ", pdf_params%rt_1(i,:)
       write(fstderr,*) "rt_2(i,:) = ", pdf_params%rt_2(i,:)
       write(fstderr,*) "varnce_rt_1(i,:) = ", pdf_params%varnce_rt_1(i,:)
       write(fstderr,*) "varnce_rt_2(i,:) = ", pdf_params%varnce_rt_2(i,:)
       write(fstderr,*) "thl_1(i,:) = ", pdf_params%thl_1(i,:)
       write(fstderr,*) "thl_2(i,:) = ", pdf_params%thl_2(i,:)
       write(fstderr,*) "varnce_thl_1(i,:) = ", pdf_params%varnce_thl_1(i,:)
       write(fstderr,*) "varnce_thl_2(i,:) = ", pdf_params%varnce_thl_2(i,:)
       write(fstderr,*) "corr_w_rt_1(i,:) = ", pdf_params%corr_w_rt_1(i,:)
       write(fstderr,*) "corr_w_rt_2(i,:) = ", pdf_params%corr_w_rt_2(i,:)
       write(fstderr,*) "corr_w_thl_1(i,:) = ", pdf_params%corr_w_thl_1(i,:)
       write(fstderr,*) "corr_w_thl_2(i,:) = ", pdf_params%corr_w_thl_2(i,:)
       write(fstderr,*) "corr_rt_thl_1(i,:) = ", pdf_params%corr_rt_thl_1(i,:)
       write(fstderr,*) "corr_rt_thl_2(i,:) = ", pdf_params%corr_rt_thl_2(i,:)
       write(fstderr,*) "alpha_thl(i,:) = ", pdf_params%alpha_thl(i,:)
       write(fstderr,*) "alpha_rt(i,:) = ", pdf_params%alpha_rt(i,:)
       write(fstderr,*) "crt_1(i,:) = ", pdf_params%crt_1(i,:)
       write(fstderr,*) "crt_2(i,:) = ", pdf_params%crt_2(i,:)
       write(fstderr,*) "cthl_1(i,:) = ", pdf_params%cthl_1(i,:)
       write(fstderr,*) "cthl_2(i,:) = ", pdf_params%cthl_2(i,:)
       write(fstderr,*) "chi_1(i,:) = ", pdf_params%chi_1(i,:)
       write(fstderr,*) "chi_2(i,:) = ", pdf_params%chi_2(i,:)
       write(fstderr,*) "stdev_chi_1(i,:) = ", pdf_params%stdev_chi_1(i,:)
       write(fstderr,*) "stdev_chi_2(i,:) = ", pdf_params%stdev_chi_2(i,:)
       write(fstderr,*) "stdev_eta_1(i,:) = ", pdf_params%stdev_eta_1(i,:)
       write(fstderr,*) "stdev_eta_2(i,:) = ", pdf_params%stdev_eta_2(i,:)
       write(fstderr,*) "covar_chi_eta_1(i,:) = ", &
                        pdf_params%covar_chi_eta_1(i,:)
       write(fstderr,*) "covar_chi_eta_2(i,:) = ", &
                        pdf_params%covar_chi_eta_2(i,:)
       write(fstderr,*) "corr_w_chi_1(i,:) = ", pdf_params%corr_w_chi_1(i,:)
       write(fstderr,*) "corr_w_chi_2(i,:) = ", pdf_params%corr_w_chi_2(i,:)
       write(fstderr,*) "corr_w_eta_1(i,:) = ", pdf_params%corr_w_eta_1(i,:)
       write(fstderr,*) "corr_w_eta_2(i,:) = ", pdf_params%corr_w_eta_2(i,:)
       write(fstderr,*) "corr_chi_eta_1(i,:) = ", &
                        pdf_params%corr_chi_eta_1(i,:)
       write(fstderr,*) "corr_chi_eta_2(i,:) = ", &
                        pdf_params%corr_chi_eta_2(i,:)
       write(fstderr,*) "rsatl_1(i,:) = ", pdf_params%rsatl_1(i,:)
       write(fstderr,*) "rsatl_2(i,:) = ", pdf_params%rsatl_2(i,:)
       write(fstderr,*) "rc_1(i,:) = ", pdf_params%rc_1(i,:)
       write(fstderr,*) "rc_2(i,:) = ", pdf_params%rc_2(i,:)
       write(fstderr,*) "cloud_frac_1(i,:) = ", pdf_params%cloud_frac_1(i,:)
       write(fstderr,*) "cloud_frac_2(i,:) = ", pdf_params%cloud_frac_2(i,:)
       write(fstderr,*) "mixt_frac(i,:) = ", pdf_params%mixt_frac(i,:)
       write(fstderr,*) "ice_supersat_frac_1(i,:) = ", &
                        pdf_params%ice_supersat_frac_1(i,:)
       write(fstderr,*) "ice_supersat_frac_2(i,:) = ", &
                        pdf_params%ice_supersat_frac_2(i,:)

    enddo ! i = 1, ngrdcol


    return

  end subroutine print_pdf_params

end module pdf_parameter_module
