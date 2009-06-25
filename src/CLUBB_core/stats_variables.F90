!-----------------------------------------------------------------------
! $Id$
!-----------------------------------------------------------------------
!  module stats_variables

!  holds pointers to variables to be written to GrADS files
!-----------------------------------------------------------------------
module stats_variables


  use stats_type, only:  & 
      stats ! Type

  use stats_precision, only:  & 
      time_precision ! Variable


  implicit none

  private ! Set Default Scope

  ! Sampling and output frequencies
  real(kind=time_precision), public :: &
    stats_tsamp, & ! Sampling interval   [s]
    stats_tout     ! Output interval     [s]

!$omp   threadprivate(stats_tsamp, stats_tout)

  logical, public ::  & 
    l_stats,  & ! Main flag to turn statistics on/off
    l_netcdf, & ! Output to NetCDF format
    l_grads     ! Output to GrADS format

!$omp   threadprivate(l_stats, l_netcdf, l_grads)

  logical, public :: & 
    l_stats_samp,   & ! Sample flag for current time step
    l_stats_first,  & ! First time step of output period
    l_stats_last      ! Last time step of output period

!$omp   threadprivate(l_stats_samp, l_stats_first, l_stats_last)

  character(len=200), public ::  & 
  fname_zt,  & ! Name of the stats file for thermodynamic grid fields
  fname_zm,  & ! Name of the stats file for momentum grid fields
  fname_sfc    ! Name of the stats file for surface only fields

!$omp   threadprivate(fname_zt, fname_zm, fname_sfc)

!       Indices for statistics in zt file

  integer, public :: & 
     ithlm, & 
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
     iwp2thlp, & 
     iwprtp2, & 
     iwp2rtp, & 
     iLscale_up, & 
     iLscale_down, & 
     itau_zt, & 
     iKh_zt, & 
     iwp2thvp, & 
     iwp2rcp, & 
     iwprtpthlp, & 
     isigma_sqd_w_zt, & 
     irho

  integer, public :: & 
     iNcm,            & ! Brian
     iNcnm, & 
     isnowslope,      & ! Adam Smith, 22 April 2008
     ised_rcm,        & ! Brian
     irsat,           & ! Brian
     irrainm,         & ! Brian
     imean_vol_rad_rain,   & ! Brian
     imean_vol_rad_cloud,  & ! COAMPS only. dschanen 6 Dec 2006
     irain_rate,      & ! Brian
     iAKm,            & ! analytic Kessler.  Vince Larson 22 May 2005 
     iAKm_est,        & ! LH Kessler.  Vince Larson  22 May 2005
     iradht,          & ! Radiative heating.
     iradht_LW,       & !   "           "   Long-wave component
     iradht_SW          !   "           "   Short-wave component

  integer, public :: & 
     iAKstd,     & 
     iAKstd_cld, & 
     iAKm_rcm, & 
     iAKm_rcc

  integer, public :: & 
     iNrm,       & ! Rain droplet number concentration
     iNim,       & ! Ice number concentration
     iNsnowm,    & ! Snow number concentration
     iNgraupelm    ! Graupel number concentration

  integer, public :: & 
     iT_in_K      ! Absolute temperature


!$omp   threadprivate(ithlm, ithvm, irtm, ircm, ium, ivm, iwm_zt, iug)
!$omp   threadprivate(ivg, icf, ip_in_Pa, iexner, iLscale, iwp3, iwpthlp2)
!$omp   threadprivate(iwp2thlp, iwprtp2, iwp2rtp, iLscale_up, iLscale_down, itau_zt)
!$omp   threadprivate(iKh_zt, iwp2thvp, iwp2rcp, iwprtpthlp, isigma_sqd_w_zt, irho)
!$omp   threadprivate(iNcm, iNcnm, isnowslope)
!$omp   threadprivate(ised_rcm, irsat, irrainm)
!$omp   threadprivate(imean_vol_rad_rain, imean_vol_rad_cloud)
!$omp   threadprivate(irain_rate, iAKm, iAKm_est)
!$omp   threadprivate(iradht, iradht_LW, iradht_SW)
!$omp   threadprivate(iT_in_K)
!$omp   threadprivate(iNrm, iNim, iNsnowm, iNgraupelm)
!$omp   threadprivate(iAKstd, iAKstd_cld, iAKm_rcm, iAKm_rcc)

  integer, public :: &
    ieff_rad_cloud, &
    ieff_rad_ice, &
    ieff_rad_snow, &
    ieff_rad_rain, &
    ieff_rad_graupel


  integer, public :: & 
    irsnowm, & 
    irgraupelm, & 
    iricem, & 
    idiam,           & ! Diameter of ice crystal           [m]
    imass_ice_cryst, & ! Mass of a single ice crystal      [kg]
    ircm_icedfs,     & ! Change in liquid water due to ice [kg/kg/s]
    iu_T_cm         ! Fallspeed of ice crystal in cm/s  [cm s^{-1}]

!$omp   threadprivate(irsnowm, irgraupelm, iricem, idiam)
!$omp   threadprivate(imass_ice_cryst, ircm_icedfs, iu_T_cm)


  ! thlm/rtm budget terms
  integer, public :: & 
    irtm_bt,       & ! rtm total time tendency
    irtm_ma,       & ! rtm mean advect. term
    irtm_ta,       & ! rtm turb. advect. term
    irtm_forcing,  & ! rtm large scale forcing term
    irtm_mc,       & ! rtm change from microphysics
    irvm_mc,       & ! rvm change from microphysics
    ircm_mc,       & ! rcm change from microphysics
    irtm_mfl,      & ! rtm change due to monotonic flux limiter
    irtm_tacl,     & ! rtm correction from turbulent advection (wprtp) clipping
    irtm_cl,       & ! rtm clipping term
    irtm_pd,       & ! thlm postive definite adj term
    ithlm_bt,      & ! thlm total time tendency
    ithlm_ma,      & ! thlm mean advect. term
    ithlm_ta,      & ! thlm turb. advect. term
    ithlm_forcing, & ! thlm large scale forcing term
    ithlm_mc,      & ! thlm change from microphysics
    ithlm_mfl,     & ! thlm change due to monotonic flux limiter
    ithlm_tacl,    & ! thlm correction from turbulent advection (wpthlp) clipping
    ithlm_cl         ! thlm clipping term

!$omp   threadprivate(irtm_bt, irtm_ma, irtm_ta, irtm_forcing, &
!$omp     irtm_mc, irtm_mfl, irtm_tacl, irtm_cl, irtm_pd, &
!$omp     irvm_mc, ircm_mc, &
!$omp     ithlm_bt, ithlm_ma, ithlm_ta, ithlm_forcing, &
!$omp     ithlm_mc, ithlm_mfl, ithlm_tacl, ithlm_cl)


  integer, public :: & 
     iwp3_bt, & 
     iwp3_ma, & 
     iwp3_ta, & 
     iwp3_tp, & 
     iwp3_ac, & 
     iwp3_bp, & 
     iwp3_pr1, & 
     iwp3_pr2, & 
     iwp3_dp1, &
     iwp3_4hd, & 
     iwp3_cl

!$omp   threadprivate(iwp3_bt, iwp3_ma, iwp3_ta, iwp3_tp, iwp3_ac)
!$omp   threadprivate(iwp3_bp, iwp3_pr1, iwp3_pr2, iwp3_dp1, iwp3_4hd, iwp3_cl)

  ! Rain mixing ratio budgets
  integer, public :: & 
     irrainm_bt, &
     irrainm_ma, &
     irrainm_sd, &
     irrainm_dff, &
     irrainm_cond, &
     irrainm_auto, &
     irrainm_accr, &
     irrainm_cond_adj, &
     irrainm_src_adj, &
     irrainm_mc, &
     irrainm_cl

!$omp   threadprivate(irrainm_bt, irrainm_ma, irrainm_sd, irrainm_dff)
!$omp   threadprivate(irrainm_cond, irrainm_auto, irrainm_accr)
!$omp   threadprivate(irrainm_cond_adj, irrainm_src_adj, irrainm_mc, irrainm_cl)

  integer, public :: &
     iNrm_bt, &
     iNrm_ma, &
     iNrm_sd, &
     iNrm_dff, &
     iNrm_cond, &
     iNrm_auto, &
     iNrm_cond_adj, &
     iNrm_src_adj, &
     iNrm_mc, &
     iNrm_cl

!$omp   threadprivate(iNrm_bt, iNrm_ma, iNrm_sd, iNrm_dff, iNrm_cond)
!$omp   threadprivate(iNrm_auto, iNrm_cond_adj, iNrm_src_adj, iNrm_mc, iNrm_cl)


  ! Snow/Ice/Graupel mixing ratio budgets
  integer, public :: & 
     irsnowm_bt, & 
     irsnowm_ma, & 
     irsnowm_sd, & 
     irsnowm_dff, & 
     irsnowm_mc, & 
     irsnowm_cl

!$omp   threadprivate(irsnowm_bt, irsnowm_ma, irsnowm_sd, irsnowm_dff)
!$omp   threadprivate(irsnowm_mc, irsnowm_cl)

  integer, public :: & 
     irgraupelm_bt, & 
     irgraupelm_ma, & 
     irgraupelm_sd, & 
     irgraupelm_dff, & 
     irgraupelm_mc, & 
     irgraupelm_cl

!$omp   threadprivate(irgraupelm_bt, irgraupelm_ma, irgraupelm_sd)
!$omp   threadprivate(irgraupelm_dff, irgraupelm_mc, irgraupelm_cl)

  integer, public :: & 
     iricem_bt, & 
     iricem_ma, & 
     iricem_sd, & 
     iricem_dff, & 
     iricem_mc, & 
     iricem_cl

!$omp   threadprivate(iricem_bt, iricem_ma, iricem_sd, iricem_dff)
!$omp   threadprivate(iricem_mc, iricem_cl)

  integer, public :: &
    iNsnowm_bt,  &
    iNsnowm_ma,  &
    iNsnowm_sd,  &
    iNsnowm_dff, &
    iNsnowm_mc,  &
    iNsnowm_cl

!$omp threadprivate(iNsnowm_bt, iNsnowm_ma, iNsnowm_sd, iNsnowm_dff, &
!$omp   iNsnowm_mc, iNsnowm_cl)

  integer, public :: &
    iNgraupelm_bt,  &
    iNgraupelm_ma,  &
    iNgraupelm_sd,  &
    iNgraupelm_dff, &
    iNgraupelm_mc,  &
    iNgraupelm_cl

!$omp threadprivate(iNgraupelm_bt, iNgraupelm_ma, iNgraupelm_sd, &
!$omp   iNgraupelm_dff, iNgraupelm_mc, iNgraupelm_cl)

  integer, public :: &
    iNim_bt,  &
    iNim_ma,  &
    iNim_sd,  &
    iNim_dff, &
    iNim_mc,  &
    iNim_cl

!$omp threadprivate(iNim_bt, iNim_ma, iNim_sd, iNim_dff, &
!$omp   iNim_mc, iNim_cl)

  integer, public :: &
    iNcm_bt,  &
    iNcm_ma,  &
    iNcm_dff, &
    iNcm_mc,  &
    iNcm_cl

!$omp threadprivate(iNcm_bt, iNcm_ma, iNcm_dff, &
!$omp   iNcm_mc, iNcm_cl)


  ! Wind budgets
  integer, public :: & 
     ivm_bt, & 
     ivm_ma, & 
     ivm_ta, & 
     ivm_gf, & 
     ivm_cf, &
     ivm_f

!$omp   threadprivate(ivm_bt, ivm_ma, ivm_ta, ivm_gf, ivm_cf)

  integer, public :: & 
     ium_bt, & 
     ium_ma, & 
     ium_ta, & 
     ium_gf, & 
     ium_cf, & 
     ium_f

!$omp   threadprivate(ium_bt, ium_ma, ium_ta, ium_gf, ium_cf)


  ! PDF parameters
  integer, public :: & 
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

!$omp   threadprivate(ia, iw1, iw2, isw1, isw2, ithl1, ithl2, isthl1)
!$omp   threadprivate(isthl2, irt1, irt2, isrt1, isrt2, irc1, irc2)
!$omp   threadprivate(irsl1, irsl2, iR1, iR2, is1, is2, iss1, iss2)
!$omp   threadprivate(irrtthl)

  integer, public :: & 
     iwp2_zt, & 
     ithlp2_zt, & 
     iwpthlp_zt, & 
     iwprtp_zt, & 
     irtp2_zt, & 
     irtpthlp_zt

!$omp   threadprivate(iwp2_zt, ithlp2_zt, iwpthlp_zt, irtp2_zt, irtpthlp_zt)

  integer, target, allocatable, dimension(:), public :: & 
    isclrm,    & ! Passive scalar mean (1)
    isclrm_f     ! Passive scalar forcing (1)

!$omp   threadprivate(isclrm, isclrm_f)

  integer, target, allocatable, dimension(:), public :: & 
    iedsclrm,   & ! Eddy-diff. scalar term (1)
    iedsclrm_f    ! Eddy-diffusivity scalar forcing (1)

!$omp   threadprivate(iedsclrm, iedsclrm_f)

  integer, public :: &
    iLH_rrainm_mc_est, & ! Latin hypercube estimate of rrainm_mc
    iLH_Nrm_mc_est,    & ! Latin hypercube estimate of Nrm_mc
    iLH_rcm_mc_est,    & ! Latin hypercube estimate of rcm_mc
    iLH_rvm_mc_est       ! Latin hypercube estimate of rvm_mc
!$omp   threadprivate(iLH_rrainm_mc_est,  iLH_Nrm_mc_est, &
!$omp                 iLH_rcm_mc_est, iLH_rvm_mc_est)

!       Indices for statistics in zm file

  integer, public :: & 
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
     ishear,     & ! Brian
     imean_w_up, &
     imean_w_down, &
     iFrad, & 
     iFrad_LW,   & ! Brian
     iFrad_SW,   & ! Brian
     iFrad_LW_up,   & 
     iFrad_SW_up,   & 
     iFrad_LW_down,   & 
     iFrad_SW_down,   & 
     iFprec,          & ! Brian
     iFcsed             ! Brian

!$omp   threadprivate(iwp2, irtp2, ithlp2, irtpthlp, iwprtp, iwpthlp)
!$omp   threadprivate(iwp4, iwpthvp, irtpthvp, ithlpthvp, itau_zm, iKh_zm)
!$omp   threadprivate(iwprcp, ithlprcp, irtprcp, ircp2, iupwp, ivpwp)
!$omp   threadprivate(irho_zm, isigma_sqd_w, iem, ishear, iFrad, iFrad_LW)
!$omp   threadprivate(iFrad_SW, iFrad_SW_up, iFrad_SW_down)
!$omp   threadprivate(iFrad_LW_up,iFrad_LW_down,iFprec, iFcsed)

  ! Sedimentation velocities
  integer, public :: & 
    iVrr,      & ! Brian
    iVNr,      & !  " "
    iVsnow,    & ! COAMPS
    iVice,     & !  " "
    iVgraupel    !  " "

!$omp   threadprivate(iVrr, iVNr, iVsnow, iVice, iVgraupel)

  integer, public :: & 
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
     iwp2_pd, & 
     iwp2_cl

!$omp   threadprivate(iwp2_bt, iwp2_ma, iwp2_ta, iwp2_ac, iwp2_bp)
!$omp   threadprivate(iwp2_pr1, iwp2_pr2, iwp2_pr3)
!$omp   threadprivate(iwp2_dp1, iwp2_dp2, iwp2_4hd)
!$omp   threadprivate(iwp2_pd, iwp2_cl)

  integer, public :: & 
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
     iwprtp_pd

!$omp   threadprivate(iwprtp_bt, iwprtp_ma, iwprtp_ta, iwprtp_tp)
!$omp   threadprivate(iwprtp_ac, iwprtp_bp, iwprtp_pr1, iwprtp_pr2)
!$omp   threadprivate(iwprtp_pr3, iwprtp_dp1, iwprtp_mfl, iwprtp_cl)
!$omp   threadprivate(iwprtp_sicl, iwprtp_pd)

  integer, public :: & 
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

!$omp   threadprivate(iwpthlp_bt, iwpthlp_ma, iwpthlp_ta, iwpthlp_tp)
!$omp   threadprivate(iwpthlp_ac, iwpthlp_bp, iwpthlp_pr1, iwpthlp_pr2)
!$omp   threadprivate(iwpthlp_pr3, iwpthlp_dp1, iwpthlp_mfl, iwpthlp_cl)
!$omp   threadprivate(iwpthlp_sicl)

!    Dr. Golaz's new variance budget terms
!    qt was changed to rt to avoid confusion

  integer, public :: & 
     irtp2_bt, & 
     irtp2_ma, & 
     irtp2_ta, & 
     irtp2_tp, & 
     irtp2_dp1, & 
     irtp2_dp2, & 
     irtp2_pd, & 
     irtp2_cl
!$omp   threadprivate(irtp2_bt, irtp2_ma, irtp2_ta, irtp2_tp)
!$omp   threadprivate(irtp2_dp1, irtp2_dp2, irtp2_pd, irtp2_cl)

  integer, public :: & 
     ithlp2_bt, & 
     ithlp2_ma, & 
     ithlp2_ta, & 
     ithlp2_tp, & 
     ithlp2_dp1, & 
     ithlp2_dp2, & 
     ithlp2_pd, & 
     ithlp2_cl

!$omp   threadprivate(ithlp2_bt, ithlp2_ma, ithlp2_ta, ithlp2_tp)
!$omp   threadprivate(ithlp2_dp1, ithlp2_dp2, ithlp2_pd, ithlp2_cl)

  integer, public :: & 
    irtpthlp_bt, & 
    irtpthlp_ma, & 
    irtpthlp_ta, & 
    irtpthlp_tp1, & 
    irtpthlp_tp2, & 
    irtpthlp_dp1, & 
    irtpthlp_dp2, & 
    irtpthlp_cl

!$omp   threadprivate(irtpthlp_bt, irtpthlp_ma, irtpthlp_ta)
!$omp   threadprivate(irtpthlp_tp1, irtpthlp_tp2, irtpthlp_dp1)
!$omp   threadprivate(irtpthlp_dp2, irtpthlp_cl)

  integer, public :: & 
    iup2, & 
    ivp2

!$omp   threadprivate(iup2, ivp2)

  integer, public :: & 
    iup2_bt, & 
    iup2_ta, & 
    iup2_tp, & 
    iup2_ma, & 
    iup2_dp1, & 
    iup2_dp2, & 
    iup2_pr1, & 
    iup2_pr2, & 
    iup2_pd, & 
    iup2_cl, & 
    ivp2_bt, & 
    ivp2_ta, & 
    ivp2_tp, & 
    ivp2_ma, & 
    ivp2_dp1, & 
    ivp2_dp2, & 
    ivp2_pr1, & 
    ivp2_pr2, & 
    ivp2_pd, & 
    ivp2_cl

!$omp   threadprivate(iup2_bt, iup2_ta, iup2_tp, iup2_ma, iup2_dp1)
!$omp   threadprivate(iup2_dp2, iup2_pr1, iup2_pr2, iup2_cl)
!$omp   threadprivate(ivp2_bt, ivp2_ta, ivp2_tp, ivp2_ma, ivp2_dp1)
!$omp   threadprivate(ivp2_dp2, ivp2_pr1, ivp2_pr2, ivp2_cl)
!$omp   threadprivate(iup2_pd, ivp2_pd)

!       Passive scalars.  Note that floating point roundoff may make
!       mathematically equivalent variables different values.
  integer,target, allocatable, dimension(:), public :: & 
    isclrprtp,           & ! sclr'(1)rt'     / rt'^2
    isclrp2,             & ! sclr'(1)^2      / rt'^2
    isclrpthvp,          & ! sclr'(1)th_v'   / rt'th_v' 
    isclrpthlp,          & ! sclr'(1)th_l'   / rt'th_l' 
    isclrprcp,           & ! sclr'(1)rc'     / rt'rc'
    iwpsclrp,            & ! w'slcr'(1)      / w'rt'
    iwp2sclrp,           & ! w'^2 sclr'(1)   / w'^2 rt'
    iwpsclrp2,           & ! w'sclr'(1)^2    / w'rt'^2
    iwpsclrprtp,         & ! w'sclr'(1)rt'   / w'rt'^2
    iwpsclrpthlp           ! w'sclr'(1)th_l' / w'rt'th_l'

  integer, target, allocatable, dimension(:), public :: & 
     iwpedsclrp ! eddy sclr'(1)w'

  ! Indices for statistics in sfc file

  integer, public :: & 
    iustar, &
    isoil_heat_flux,&
    iveg_T_in_K,&
    isfc_soil_T_in_K, &
    ideep_soil_T_in_K,& 
    ilh, & 
    ish, & 
    icc, & 
    ilwp, &
    ivwp, &        ! nielsenb
    iiwp, &        ! nielsenb
    iswp, &        ! nielsenb
    izb, & 
    izi, & 
    irain,    &    ! Brian
    irain_flux,   &    ! Brian
    irrainm_sfc, & ! Brian
    iwpthlp_sfc, &
    iwprtp_sfc, &
    iupwp_sfc, &
    ivpwp_sfc, &
    ithlm_vert_avg, &
    irtm_vert_avg, &
    ium_vert_avg, &
    ivm_vert_avg, &
    iwp2_vert_avg, & ! nielsenb
    iup2_vert_avg, &
    ivp2_vert_avg, &
    irtp2_vert_avg, &
    ithlp2_vert_avg

  integer, public :: & 
    iwp23_matrix_condt_num, & 
    irtm_matrix_condt_num, & 
    ithlm_matrix_condt_num, & 
    irtp2_matrix_condt_num, & 
    ithlp2_matrix_condt_num, & 
    irtpthlp_matrix_condt_num, & 
    iup2_vp2_matrix_condt_num, & 
    iwindm_matrix_condt_num

  integer, public :: & 
    imorr_rain_rate, &
    imorr_snow_rate


!$omp threadprivate(iustar, isoil_heat_flux, iveg_T_in_K, isfc_soil_T_in_K, ideep_soil_T_in_K, &
!$omp   ilh, ish, icc, ilwp, izb, izi, &
!$omp   irain, irain_flux, irrainm_sfc, &
!$omp   iwpthlp_sfc, iwprtp_sfc, iupwp_sfc, ivpwp_sfc, &
!$omp   ithlm_vert_avg, irtm_vert_avg, ium_vert_avg, ivm_vert_avg, &
!$omp   iwp2_vert_avg, iup2_vert_avg, ivp2_vert_avg, irtp2_vert_avg, ithlp2_vert_avg, &
!$omp   iwp23_matrix_condt_num, irtm_matrix_condt_num, ithlm_matrix_condt_num, &
!$omp   irtp2_matrix_condt_num, ithlp2_matrix_condt_num, irtpthlp_matrix_condt_num, &
!$omp   iup2_vp2_matrix_condt_num, iwindm_matrix_condt_num, &
!$omp   imorr_rain_rate, imorr_snow_rate)

  ! Variables that contains all the statistics

  type (stats), target, public :: zt,   & ! zt grid
                                  zm,   & ! zm grid
                                  sfc     ! sfc

!$omp   threadprivate(zt, zm, sfc)

  ! Scratch space

  real, allocatable, public :: ztscr01(:), ztscr02(:), ztscr03(:), & 
                               ztscr04(:), ztscr05(:), ztscr06(:), & 
                               ztscr07(:), ztscr08(:), ztscr09(:), & 
                               ztscr10(:), ztscr11(:), ztscr12(:), & 
                               ztscr13(:), ztscr14(:), ztscr15(:), & 
                               ztscr16(:), ztscr17(:), ztscr18(:), &
                               ztscr19(:), ztscr20(:), ztscr21(:)

!$omp   threadprivate(ztscr01, ztscr02, ztscr03, ztscr04, ztscr05)
!$omp   threadprivate(ztscr06, ztscr07, ztscr08, ztscr09, ztscr10)
!$omp   threadprivate(ztscr11, ztscr12, ztscr13, ztscr14, ztscr15)
!$omp   threadprivate(ztscr16, ztscr17, ztscr18, ztscr19, ztscr20)
!$omp   threadprivate(ztscr21)

  real, allocatable, public :: zmscr01(:), zmscr02(:), zmscr03(:), &
                               zmscr04(:), zmscr05(:), zmscr06(:), & 
                               zmscr07(:), zmscr08(:), zmscr09(:), & 
                               zmscr10(:), zmscr11(:), zmscr12(:), & 
                               zmscr13(:), zmscr14(:), zmscr15(:), &
                               zmscr16(:), zmscr17(:)

!$omp   threadprivate(zmscr01, zmscr02, zmscr03, zmscr04, zmscr05)
!$omp   threadprivate(zmscr06, zmscr07, zmscr08, zmscr09, zmscr10)
!$omp   threadprivate(zmscr11, zmscr12, zmscr13, zmscr14, zmscr15)
!$omp   threadprivate(zmscr16, zmscr17)


end module stats_variables

