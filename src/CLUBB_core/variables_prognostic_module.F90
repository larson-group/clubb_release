!-----------------------------------------------------------------------
! $Id$
module variables_prognostic_module

!       This module contains definitions of all prognostic
!       arrays used in the single column model, as well as subroutines
!       to allocate, deallocate and initialize them.

!       Note that while these are all same dimension, there is a
!       thermodynamic grid and a momentum grid, and the grids have
!       different points.
!-----------------------------------------------------------------------

  implicit none

  private ! Set Default Scoping

  public ::  & 
    setup_prognostic_variables,  & 
    cleanup_prognostic_variables

  ! Prognostic variables
  real, target, allocatable, dimension(:), public :: & 
    um,      & ! u wind                        [m/s]
    vm,      & ! v wind                        [m/s]
    upwp,    & ! vertical u momentum flux      [m^2/s^2]
    vpwp,    & ! vertical v momentum flux      [m^2/s^2]
    up2,     & ! u'^2                          [m^2/s^2]
    vp2,     & ! v'^2                          [m^2/s^2]
    thlm,    & ! liquid potential temperature  [K]
    rtm,     & ! total water mixing ratio      [kg/kg]
    wprtp,   & ! w'rt'                         [(kg/kg) m/s]
    wpthlp,  & ! w'thl'                        [m K/s]
    wprcp,   & ! w'rc'                         [(kg/kg) m/s]
    wp2,     & ! w'^2                          [m^2/s^2]
    wp3,     & ! w'^3                          [m^3/s^3]
    rtp2,    & ! rt'^2                         [(kg/kg)^2]
    thlp2,   & ! thl'^2                        [K^2]
    rtpthlp    ! rt'thl'                       [kg/kg K]

!$omp   threadprivate(um, vm, upwp, vpwp, up2, vp2)
!$omp   threadprivate(thlm, rtm, wprtp, wpthlp, wprcp)
!$omp   threadprivate(wp2, wp3, rtp2, thlp2, rtpthlp)

  real, target, allocatable, dimension(:), public :: & 
    p_in_Pa,         & ! Pressure (Pa) (thermodynamic levels)          [Pa]
    exner,           & ! Exner function = ( p / p0 ) ** kappa          [-]
    rho,             & ! Density (thermodynamic levels)                [kg/m^3]
    rho_zm,          & ! Density on momentum levels                    [kg/m^3]
    rho_ds_zm,       & ! Dry, static density (momentum levels)         [kg/m^3]
    rho_ds_zt,       & ! Dry, static density (thermodynamic levels)    [kg/m^3]
    invrs_rho_ds_zm, & ! Inverse dry, static density (momentum levs.)  [m^3/kg]
    invrs_rho_ds_zt, & ! Inverse dry, static density (thermo. levs.)   [m^3/kg]
    thv_ds_zm,       & ! Dry, base-state theta_v (momentum levels)     [K]
    thv_ds_zt,       & ! Dry, base-state theta_v (thermodynamic levs.) [K]
    thlm_forcing,    & ! thlm large-scale forcing                      [K/s]
    rtm_forcing,     & ! rtm large-scale forcing                       [kg/kg/s]
    um_forcing,      & ! u wind forcing                                [m/s/s] 
    vm_forcing         ! v wind forcing                                [m/s/s]

!$omp   threadprivate(p_in_Pa, exner, rho, rho_zm, rho_ds_zm, &
!$omp     rho_ds_zt, invrs_rho_ds_zm, invrs_rho_ds_zt, thv_ds_zm, &
!$omp     thv_ds_zt, thlm_forcing, rtm_forcing, um_forcing, vm_forcing)

  ! Imposed large scale w
  real, target, allocatable, dimension(:), public :: & 
    wm_zm, & ! w on momentum levels              [m/s]
    wm_zt    ! w on thermodynamic levels         [m/s]

!$omp   threadprivate(wm_zm, wm_zt)

  ! Cloud water variables
  real, target, allocatable, dimension(:), public :: & 
    rcm,          & ! Cloud water mixing ratio                 [kg/kg]
    cloud_frac,   & ! Cloud fraction                           [-]
    rcm_in_layer, & ! Cloud water mixing ratio in cloud layer  [kg/kg]
    cloud_cover     ! Cloud cover                              [-]

!$omp   threadprivate(rcm, cloud_frac, rcm_in_layer, cloud_cover)

  ! Surface fluxes
  real, public ::  & 
    wpthlp_sfc,        & ! w'thl'      [m K/s]
    wprtp_sfc,         & ! w'rt'       [m kg/(kg s)]
    upwp_sfc, vpwp_sfc   ! u'w' & v'w' [m^2/s^2]

!$omp   threadprivate(wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc)

  ! Surface fluxes for passive scalars
  real, dimension(:), allocatable, public :: & 
    wpsclrp_sfc,     & ! w'sclr' at surface    [units m/s]
    wpedsclrp_sfc      ! w'edsclr' at surface  [units m/s]

!$omp   threadprivate(wpsclrp_sfc, wpedsclrp_sfc)

  ! More surface data
  real, public ::  & 
    Tsfc,  & ! surface temperature     [K]
    psfc,  & ! surface pressure        [Pa]
    SE,    & ! sensible heat flux      [K/s]
    LE       ! latent heat flux        [1/s]

!$omp   threadprivate(Tsfc, psfc, SE, LE)

  ! Passive scalars
  real, target, allocatable, dimension(:,:), public :: & 
    sclrm,           & ! Mean passive scalars           [units vary]
    sclrp2,          & ! sclr'^2                        [units^2]
    sclrprtp,        & ! sclr'rt'                       [units kg/kg]
    sclrpthlp,       & ! sclr'th_l'                     [units K]
    sclrm_forcing,   & ! Scalars' forcing               [units/s]
    edsclrm,         & ! Mean eddy-diffusivity scalars  [units vary]
    edsclrm_forcing, & ! Eddy-diff. scalars forcing     [units/s]
    wpsclrp            ! w'sclr'                        [units vary m/s]

!$omp   threadprivate(sclrm, sclrp2, sclrprtp, sclrpthlp, sclrm_forcing, &
!$omp     edsclrm, edsclrm_forcing, wpsclrp)

! PDF parameters
  public :: pdf_parameter

  type pdf_parameter
    real, pointer, dimension(:) ::  &
      w1,          & ! Mean of w for 1st normal distribution                 [m/s]
      w2,          & ! Mean of w for 2nd normal distribution                 [m/s]
      varnce_w1,   & ! Variance of w for 1st normal distribution         [m^2/s^2]
      varnce_w2,   & ! Variance of w for 2nd normal distribution         [m^2/s^2]
      rt1,         & ! Mean of r_t for 1st normal distribution             [kg/kg]
      rt2,         & ! Mean of r_t for 2nd normal distribution             [kg/kg]
      varnce_rt1,  & ! Variance of r_t for 1st normal distribution     [kg^2/kg^2]
      varnce_rt2,  & ! Variance of r_t for 2nd normal distribution     [kg^2/kg^2]
      crt1,        & ! Coefficient for s'                                      [-]
      crt2,        & ! Coefficient for s'                                      [-]
      cthl1,       & ! Coefficient for s'                                    [1/K]
      cthl2,       & ! Coefficient for s'                                    [1/K]
      thl1,        & ! Mean of th_l for 1st normal distribution                [K]
      thl2,        & ! Mean of th_l for 2nd normal distribution                [K]
      varnce_thl1, & ! Variance of th_l for 1st normal distribution          [K^2]
      varnce_thl2, & ! Variance of th_l for 2nd normal distribution          [K^2]
      mixt_frac,   & ! Weight of 1st normal distribution (Sk_w dependent)      [-]
      rc1,         & ! Mean of r_c for 1st normal distribution             [kg/kg]
      rc2,         & ! Mean of r_c for 2nd normal distribution             [kg/kg]
      rsl1,        & ! Mean of r_sl for 1st normal distribution            [kg/kg]
      rsl2,        & ! Mean of r_sl for 2nd normal distribution            [kg/kg]
      cloud_frac1, & ! Cloud fraction for 1st normal distribution              [-]
      cloud_frac2, & ! Cloud fraction for 2nd normal distribution              [-]
      s1,          & ! Mean of s for 1st normal distribution               [kg/kg]
      s2,          & ! Mean of s for 2nd normal distribution               [kg/kg]
      stdev_s1,    & ! Standard deviation of s for 1st normal distribution [kg/kg]
      stdev_s2,    & ! Standard deviation of s for 2nd normal distribution [kg/kg]
      rrtthl,      & ! Within-a-normal correlation of r_t and th_l             [-]
      alpha_thl,   & ! Factor relating to normalized variance for th_l         [-]
      alpha_rt       ! Factor relating to normalized variance for r_t          [-]
  end type

  type(pdf_parameter), target, public :: &
    pdf_params

!$omp threadprivate(pdf_params)

  contains
!-----------------------------------------------------------------------
  subroutine setup_prognostic_variables( nzmax )

!       Description:
!       Allocates and Initializes prognostic scalar and array variables
!       for the HOC model code

!       References:
!       None
!-----------------------------------------------------------------------
    use constants, only:  & 
        rttol, &
        thltol, &
        wtol_sqd

    use parameters_model, only: & 
        sclr_dim,  & ! Variable(s)
        edsclr_dim

    implicit none

    integer, intent(in) :: nzmax ! Number of grid levels [-]

    integer :: i

!   --- Allocation ---

! Prognostic variables

    allocate( um(1:nzmax) )        ! u wind
    allocate( vm(1:nzmax) )        ! v wind

    allocate( upwp(1:nzmax) )      ! vertical u momentum flux
    allocate( vpwp(1:nzmax) )      ! vertical v momentum flux

    allocate( up2(1:nzmax) )
    allocate( vp2(1:nzmax) )

    allocate( thlm(1:nzmax) )      ! liquid potential temperature
    allocate( rtm(1:nzmax) )       ! total water mixing ratio
    allocate( wprtp(1:nzmax) )     ! w'rt'
    allocate( wpthlp(1:nzmax) )    ! w'thl'
    allocate( wprcp(1:nzmax) )     ! w'rc'
    allocate( wp2(1:nzmax) )       ! w'^2
    allocate( wp3(1:nzmax) )       ! w'^3
    allocate( rtp2(1:nzmax) )      ! rt'^2
    allocate( thlp2(1:nzmax) )     ! thl'^2
    allocate( rtpthlp(1:nzmax) )   ! rt'thlp'

    allocate( p_in_Pa(1:nzmax) )         ! pressure (pascals)
    allocate( exner(1:nzmax) )           ! exner function
    allocate( rho(1:nzmax) )             ! density: t points
    allocate( rho_zm(1:nzmax) )          ! density: m points
    allocate( rho_ds_zm(1:nzmax) )       ! dry, static density: m-levs
    allocate( rho_ds_zt(1:nzmax) )       ! dry, static density: t-levs
    allocate( invrs_rho_ds_zm(1:nzmax) ) ! inv. dry, static density: m-levs
    allocate( invrs_rho_ds_zt(1:nzmax) ) ! inv. dry, static density: t-levs
    allocate( thv_ds_zm(1:nzmax) )       ! dry, base-state theta_v: m-levs
    allocate( thv_ds_zt(1:nzmax) )       ! dry, base-state theta_v: t-levs

    allocate( thlm_forcing(1:nzmax) )    ! thlm ls forcing
    allocate( rtm_forcing(1:nzmax) )     ! rtm ls forcing
    allocate( um_forcing(1:nzmax) )      ! u forcing
    allocate( vm_forcing(1:nzmax) )      ! v forcing

    ! Imposed large scale w

    allocate( wm_zm(1:nzmax) )       ! momentum levels
    allocate( wm_zt(1:nzmax) )       ! thermodynamic levels

    ! Cloud water variables

    allocate( rcm(1:nzmax) )
    allocate( cloud_frac(1:nzmax) )
    allocate( rcm_in_layer(1:nzmax) )
    allocate( cloud_cover(1:nzmax) )

    ! Passive scalar variables
    ! Note that sclr_dim can be 0
    allocate( wpsclrp_sfc(1:sclr_dim) )
    allocate( sclrm(1:nzmax, 1:sclr_dim) )
    allocate( sclrp2(1:nzmax, 1:sclr_dim) )
    allocate( sclrm_forcing(1:nzmax, 1:sclr_dim) )
    allocate( sclrprtp(1:nzmax, 1:sclr_dim) )
    allocate( sclrpthlp(1:nzmax, 1:sclr_dim) )

    allocate( wpedsclrp_sfc(1:edsclr_dim) )
    allocate( edsclrm_forcing(1:nzmax, 1:edsclr_dim) )

    allocate( edsclrm(1:nzmax, 1:edsclr_dim) )
    allocate( wpsclrp(1:nzmax, 1:sclr_dim) )

    ! Variables for pdf closure scheme
    allocate( pdf_params%w1(1:nzmax),          pdf_params%w2(1:nzmax),  &
              pdf_params%varnce_w1(1:nzmax),   pdf_params%varnce_w2(1:nzmax),  &
              pdf_params%rt1(1:nzmax),         pdf_params%rt2(1:nzmax),  &
              pdf_params%varnce_rt1(1:nzmax),  pdf_params%varnce_rt2(1:nzmax),  &
              pdf_params%thl1(1:nzmax),        pdf_params%thl2(1:nzmax),  &
              pdf_params%varnce_thl1(1:nzmax), pdf_params%varnce_thl2(1:nzmax),  &
              pdf_params%mixt_frac(1:nzmax),   pdf_params%rrtthl(1:nzmax),  &
              pdf_params%rc1(1:nzmax),         pdf_params%rc2(1:nzmax),  &
              pdf_params%rsl1(1:nzmax),        pdf_params%rsl2(1:nzmax),  &
              pdf_params%cloud_frac1(1:nzmax), pdf_params%cloud_frac2(1:nzmax),  &
              pdf_params%s1(1:nzmax),          pdf_params%s2(1:nzmax),  &
              pdf_params%stdev_s1(1:nzmax),    pdf_params%stdev_s2(1:nzmax),  &
              pdf_params%alpha_thl(1:nzmax),   pdf_params%alpha_rt(1:nzmax), &
              pdf_params%crt1(1:nzmax),        pdf_params%crt2(1:nzmax), &
              pdf_params%cthl1(1:nzmax),       pdf_params%cthl2(1:nzmax) )



!--------- Set initial values for array variables ---------

    ! Prognostic variables

    um(1:nzmax)      = 0.0     ! u wind
    vm (1:nzmax)     = 0.0     ! v wind

    upwp(1:nzmax)    = 0.0     ! vertical u momentum flux
    vpwp(1:nzmax)    = 0.0     ! vertical v momentum flux

    up2(1:nzmax)     = wtol_sqd ! u'^2
    vp2(1:nzmax)     = wtol_sqd ! v'^2
    wp2(1:nzmax)     = wtol_sqd ! w'^2

    thlm(1:nzmax)    = 0.0         ! liquid potential temperature
    rtm(1:nzmax)     = 0.0         ! total water mixing ratio
    wprtp(1:nzmax)   = 0.0         ! w'rt'
    wpthlp(1:nzmax)  = 0.0         ! w'thl'
    wprcp(1:nzmax)   = 0.0         ! w'rc'
    wp3(1:nzmax)     = 0.0         ! w'^3
    rtp2(1:nzmax)    = rttol**2    ! rt'^2
    thlp2(1:nzmax)   = thltol**2   ! thl'^2
    rtpthlp(1:nzmax) = 0.0         ! rt'thl'

    p_in_Pa(1:nzmax)= 0.0           ! pressure (Pa)
    exner(1:nzmax) = 0.0            ! exner
    rho(1:nzmax)  = 0.0             ! density on thermo. levels
    rho_zm(1:nzmax)  = 0.0          ! density on moment. levels
    rho_ds_zm(1:nzmax) = 0.0        ! dry, static density: m-levs
    rho_ds_zt(1:nzmax) = 0.0        ! dry, static density: t-levs
    invrs_rho_ds_zm(1:nzmax) = 0.0  ! inv. dry, static density: m-levs
    invrs_rho_ds_zt(1:nzmax) = 0.0  ! inv. dry, static density: t-levs
    thv_ds_zm(1:nzmax) = 0.0        ! dry, base-state theta_v: m-levs
    thv_ds_zt(1:nzmax) = 0.0        ! dry, base-state theta_v: t-levs

    thlm_forcing(1:nzmax) = 0.0     ! thlm large-scale forcing
    rtm_forcing(1:nzmax)  = 0.0     ! rtm large-scale forcing
    um_forcing(1:nzmax) = 0.0       ! u forcing
    vm_forcing(1:nzmax) = 0.0       ! v forcing

    ! Imposed large scale w

    wm_zm(1:nzmax) = 0.0      ! Momentum levels
    wm_zt(1:nzmax) = 0.0      ! Thermodynamic levels

    ! Cloud water variables

    rcm(1:nzmax)          = 0.0
    cloud_frac(1:nzmax)   = 0.0
    rcm_in_layer(1:nzmax) = 0.0
    cloud_cover(1:nzmax)  = 0.0

    ! Variables for PDF closure scheme
    pdf_params%w1          = 0.0
    pdf_params%w2          = 0.0
    pdf_params%varnce_w1   = 0.0
    pdf_params%varnce_w2   = 0.0
    pdf_params%rt1         = 0.0
    pdf_params%rt2         = 0.0
    pdf_params%varnce_rt1  = 0.0
    pdf_params%varnce_rt2  = 0.0
    pdf_params%thl1        = 0.0
    pdf_params%thl2        = 0.0
    pdf_params%varnce_thl1 = 0.0
    pdf_params%varnce_thl2 = 0.0
    pdf_params%mixt_frac   = 0.0
    pdf_params%rc1         = 0.0
    pdf_params%rc2         = 0.0
    pdf_params%rsl1        = 0.0
    pdf_params%rsl2        = 0.0
    pdf_params%cloud_frac1 = 0.0
    pdf_params%cloud_frac2 = 0.0
    pdf_params%s1          = 0.0
    pdf_params%s2          = 0.0
    pdf_params%stdev_s1    = 0.0
    pdf_params%stdev_s2    = 0.0
    pdf_params%rrtthl      = 0.0
    pdf_params%alpha_thl   = 0.0
    pdf_params%alpha_rt    = 0.0
    pdf_params%crt1        = 0.0
    pdf_params%crt2        = 0.0
    pdf_params%cthl1       = 0.0
    pdf_params%cthl2       = 0.0

    ! Surface fluxes
    wpthlp_sfc = 0.0
    wprtp_sfc  = 0.0
    upwp_sfc   = 0.0
    vpwp_sfc   = 0.0

    ! Passive scalars
    do i = 1, sclr_dim, 1
      wpsclrp_sfc(i)   = 0.0

      sclrm(1:nzmax,i)         = 0.0
      sclrp2(1:nzmax,i)        = 0.0
      sclrprtp(1:nzmax,i)      = 0.0
      sclrpthlp(1:nzmax,i)     = 0.0
      sclrm_forcing(1:nzmax,i) = 0.0
      wpsclrp(1:nzmax,i)         = 0.0
    end do

    do i = 1, edsclr_dim, 1
      wpedsclrp_sfc(i) = 0.0

      edsclrm(1:nzmax,i)         = 0.0
      edsclrm_forcing(1:nzmax,i) = 0.0
    end do

    return
  end subroutine setup_prognostic_variables
!-----------------------------------------------------------------------
  subroutine cleanup_prognostic_variables
    implicit none

    ! Prognostic variables

    deallocate( um )        ! u wind
    deallocate( vm )        ! v wind

    deallocate( upwp )      ! vertical u momentum flux
    deallocate( vpwp )      ! vertical v momentum flux

    deallocate( up2, vp2 )

    deallocate( thlm )      ! liquid potential temperature
    deallocate( rtm )       ! total water mixing ratio
    deallocate( wprtp )     ! w'rt'
    deallocate( wpthlp )    ! w'thl'
    deallocate( wprcp )     ! w'rc'
    deallocate( wp2 )       ! w'^2
    deallocate( wp3 )       ! w'^3
    deallocate( rtp2 )      ! rt'^2
    deallocate( thlp2 )     ! thl'^2
    deallocate( rtpthlp )   ! rt'thl'

    deallocate( p_in_Pa )         ! pressure
    deallocate( exner )           ! exner
    deallocate( rho )             ! density: t points
    deallocate( rho_zm )          ! density: m points
    deallocate( rho_ds_zm )       ! dry, static density: m-levs
    deallocate( rho_ds_zt )       ! dry, static density: t-levs
    deallocate( invrs_rho_ds_zm ) ! inv. dry, static density: m-levs
    deallocate( invrs_rho_ds_zt ) ! inv. dry, static density: t-levs
    deallocate( thv_ds_zm )       ! dry, base-state theta_v: m-levs
    deallocate( thv_ds_zt )       ! dry, base-state theta_v: t-levs

    deallocate( thlm_forcing )
    deallocate( rtm_forcing )
    deallocate( um_forcing )
    deallocate( vm_forcing )

    ! Imposed large scale w

    deallocate( wm_zm )     ! momentum levels
    deallocate( wm_zt )     ! thermodynamic levels

    ! Cloud water variables

    deallocate( rcm )
    deallocate( cloud_frac )
    deallocate( rcm_in_layer )
    deallocate( cloud_cover )

    ! Variable for pdf closure scheme
    deallocate( pdf_params%w1,          pdf_params%w2,  &
                pdf_params%varnce_w1,   pdf_params%varnce_w2,  &
                pdf_params%rt1,         pdf_params%rt2,  &
                pdf_params%varnce_rt1,  pdf_params%varnce_rt2,  &
                pdf_params%thl1,        pdf_params%thl2,  &
                pdf_params%varnce_thl1, pdf_params%varnce_thl2,  &
                pdf_params%mixt_frac,   pdf_params%rrtthl,  &
                pdf_params%rc1,         pdf_params%rc2,  &
                pdf_params%rsl1,        pdf_params%rsl2,  &
                pdf_params%cloud_frac1, pdf_params%cloud_frac2,  &
                pdf_params%s1,          pdf_params%s2,  &
                pdf_params%stdev_s1,    pdf_params%stdev_s2,  &
                pdf_params%alpha_thl,   pdf_params%alpha_rt, &
                pdf_params%crt1,        pdf_params%crt2, &
                pdf_params%cthl1,       pdf_params%cthl2 )

    ! Passive scalars
    deallocate( wpsclrp_sfc, wpedsclrp_sfc )
    deallocate( sclrm )
    deallocate( sclrp2 )
    deallocate( sclrprtp )
    deallocate( sclrpthlp )
    deallocate( sclrm_forcing )
    deallocate( wpsclrp )

    deallocate( edsclrm )
    deallocate( edsclrm_forcing )

    return
  end subroutine cleanup_prognostic_variables

end module variables_prognostic_module
