!-----------------------------------------------------------------------
! $Id: hoc_stats.F,v 1.1 2007-09-08 20:57:57 faschinj Exp $
#ifdef STATS
!-----------------------------------------------------------------------
!  module hoc_stats

!  holds pointers to variables to be written to GrADS files
!-----------------------------------------------------------------------
      module hoc_stats

        use stats_type

        ! Sampling and output frequencies
        real :: stats_tsamp       ! Sampling interval   [s]
        real :: stats_tout        ! Output interval     [s]

        logical :: 
     .  lstats,  ! Main flag to turn statistics on/off
     .  lnetcdf, ! Output to NetCDF format
     .  lgrads   ! Output to GrADS format

!$omp   threadprivate(lstats, lnetcdf, lgrads)
!$omp   threadprivate(stats_tsamp, stats_tout)

        logical ::
     .  lstats_samp,   ! Sample flag for current time step
     .  lstats_first,  ! First time step of output period
     .  lstats_last    ! Last time step of output period

!$omp   threadprivate(lstats_samp, lstats_first, lstats_last)

        character(len=200) :: fname_zt
        character(len=200) :: fname_zm
        character(len=200) :: fname_sfc

!$omp   threadprivate(fname_zt, fname_zm, fname_sfc)

!       Indices for statistics in zt file

        integer ::
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
     .     iexner,
     .     iLscale,
     .     iwp3,
     .     iwpthlp2,
     .     iwp2thlp,
     .     iwprtp2,
     .     iwp2rtp,
     .     ilup,
     .     ildown,
     .     itaut,
     .     iKht,
     .     iwp2thvp,
     .     iwp2rcp,
     .     iwprtpthlp,
     .     isct,
     .     irhot,
     .     iNcm,           ! Brian
     .     iNcnm,
     .     iNim,
     .     ised_rcm,       ! Brian
     .     irsm,           ! Brian
     .     irrm,           ! Brian
     .     iNrm,           ! Brian
     .     imean_vol_rad_rain,  ! Brian
     .     imean_vol_rad_cloud, ! COAMPS only. dschanen 6 Dec 2006
     .     irain_rate,     ! Brian
     .     iAKm,           ! analytic Kessler.  Vince Larson 22 May 2005 
     .     iAKm_est,       ! LH Kessler.  Vince Larson  22 May 2005
     .     iradht,         ! Radiative heating. dschanen 7 Oct 2005
     .     iradht_LW,      !   "           "   Long-wave
     .     iradht_SW       !   "           "   Short-wave

!$omp   threadprivate(ithlm, ithvm, irtm, ircm, ium, ivm, iwmt, iug)
!$omp   threadprivate(ivg, icf, ip, iexner, iLscale, iwp3, iwpthlp2)
!$omp   threadprivate(iwp2thlp, iwprtp2, iwp2rtp, ilup, ildown, itaut)
!$omp   threadprivate(iKht, iwp2thvp, iwp2rcp, iwprtpthlp, isct, irhot)
!$omp   threadprivate(iNcm, iNcnm, iNim, ised_rcm, irsm, irrm, iNrm)
!$omp   threadprivate(imean_vol_rad_rain, imean_vol_rad_cloud)
!$omp   threadprivate(irain_rate, iAKm, iAKm_est)
!$omp   threadprivate(iradht, iradht_LW, iradht_SW)

        integer ::
     .    irsnowm,
     .    irgraupelm,
     .    iricem,
     .    idiam,          ! Necessary for icedfs scheme.  Michael Falk
     .    im,             ! Necessary for icedfs scheme.  Michael Falk
     .    idqc_dt_icedfs, ! Necessary for icedfs scheme.  Michael Falk
     .    iu_T_cm         ! Necessary for icedfs scheme.  Michael Falk

!$omp   threadprivate(irsnowm, irgraupelm, iricem, idiam)
!$omp   threadprivate(im, idqc_dt_icedfs, iu_T_cm)


        ! thlm/rtm budget terms
        integer ::
     .    irtm_bt,      ! rtm total time tendency
     .    irtm_ma,      ! rtm mean advect. term
     .    irtm_ta,      ! rtm turb. advect. term
     .    irtm_forcing, ! rtm large scale forcing term
     .    irtm_mc,      ! rtm change from microphysics
     .    ithlm_bt,     ! thlm total time tendency
     .    ithlm_ma,     ! thlm mean advect. term
     .    ithlm_ta,     ! thlm turb. advect. term
     .    ithlm_forcing,! thlm large scale forcing term
     .    ithlm_mc      ! thlm change from microphysics

!$omp   threadprivate(irtm_bt, irtm_ma, irtm_ta, irtm_forcing, irtm_mc)
!$omp   threadprivate(ithlm_bt, ithlm_ma, ithlm_ta)
!$omp   threadprivate(ithlm_forcing, ithlm_mc)


        integer ::
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

!$omp   threadprivate(iwp3_bt, iwp3_ma, iwp3_ta, iwp3_tp, iwp3_ac)
!$omp   threadprivate(iwp3_bp, iwp3_pr1, iwp3_pr2, iwp3_dp1, iwp3_cl)

        ! Rain mixing ratio budgets
        integer ::
     .     irrm_bt,
     .     irrm_ma,
     .     irrm_sd,
     .     irrm_dff,
     .     irrm_cond,
     .     irrm_auto,
     .     irrm_accr,
     .     irrm_cond_adj,
     .     irrm_mc,
     .     irrm_cl

!$omp   threadprivate(irrm_bt, irrm_ma, irrm_sd, irrm_dff, irrm_cond)
!$omp   threadprivate(irrm_auto, irrm_accr, irrm_cond_adj, irrm_mc)
!$omp   threadprivate(irrm_cl)

        integer ::
     .     iNrm_bt,
     .     iNrm_ma,
     .     iNrm_sd,
     .     iNrm_dff,
     .     iNrm_cond,
     .     iNrm_auto,
     .     iNrm_cond_adj,
     .     iNrm_mc,
     .     iNrm_cl

!$omp   threadprivate(iNrm_bt, iNrm_ma, iNrm_sd, iNrm_dff, iNrm_cond)
!$omp   threadprivate(iNrm_auto, iNrm_cond_adj, iNrm_mc, iNrm_cl)


        ! Snow/Ice/Graupel mixing ratio budgets
        integer ::
     .     irsnowm_bt,
     .     irsnowm_ma,
     .     irsnowm_sd,
     .     irsnowm_dff,
     .     irsnowm_mc,
     .     irsnowm_cl

!$omp   threadprivate(irsnowm_bt, irsnowm_ma, irsnowm_sd, irsnowm_dff)
!$omp   threadprivate(irsnowm_mc, irsnowm_cl)

        integer ::
     .     irgraupelm_bt,
     .     irgraupelm_ma,
     .     irgraupelm_sd,
     .     irgraupelm_dff,
     .     irgraupelm_mc,
     .     irgraupelm_cl

!$omp   threadprivate(irgraupelm_bt, irgraupelm_ma, irgraupelm_sd)
!$omp   threadprivate(irgraupelm_dff, irgraupelm_mc, irgraupelm_cl)

        integer ::
     .     iricem_bt,
     .     iricem_ma,
     .     iricem_sd,
     .     iricem_dff,
     .     iricem_mc,
     .     iricem_cl

!$omp   threadprivate(iricem_bt, iricem_ma, iricem_sd, iricem_dff)
!$omp   threadprivate(iricem_mc, iricem_cl)

        ! Wind budgets
        integer ::
     .     ivm_bt,
     .     ivm_ma,
     .     ivm_ta,
     .     ivm_gf,
     .     ivm_cf

!$omp   threadprivate(ivm_bt, ivm_ma, ivm_ta, ivm_gf, ivm_cf)

        integer ::
     .     ium_bt,
     .     ium_ma,
     .     ium_ta,
     .     ium_gf,
     .     ium_cf

!$omp   threadprivate(ium_bt, ium_ma, ium_ta, ium_gf, ium_cf)


        ! PDF parameters
        integer ::
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

!$omp   threadprivate(ia, iw1, iw2, isw1, isw2, ithl1, ithl2, isthl1)
!$omp   threadprivate(isthl2, irt1, irt2, isrt1, isrt2, irc1, irc2)
!$omp   threadprivate(irsl1, irsl2, iR1, iR2, is1, is2, iss1, iss2)
!$omp   threadprivate(irrtthl)

        integer ::
     .     iwp2zt,
     .     ithlp2zt,
     .     iwpthlpzt,
     .     iwprtpzt,
     .     irtp2zt,
     .     irtpthlpzt

!$omp   threadprivate(iwp2zt, ithlp2zt, iwpthlpzt, irtp2zt, irtpthlpzt)

#ifdef SCALARS
        integer ::
     .    isclram,   ! Passive scalar mean (currently mean thl)
     .    isclrbm,   ! Passive scalar mean (currently mean rt)
     .    isclram_f, ! Passive scalar tendency (currently mean thl forcing)
     .    isclrbm_f  ! Passive scalar tendency (currently mean rt forcing)

        integer ::
     .    iedsclram, ! Eddy-diff. scalar term (currently rt)
     .    iedsclrbm  ! Eddy-diff. scalar term (currently thl)

!$omp   threadprivate(isclram, isclrbm, isclram_f, isclrbm_f)
!$omp   threadprivate(iedsclram, iedsclrbm)

#endif /*SCALARS*/

!       Indices for statistics in zm file

        integer ::
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
     .     iKhm,
     .     iwprcp,
     .     ithlprcp,
     .     irtprcp,
     .     ircp2,
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
     .     iFcsed,    ! Brian
     .     iup2_plus_vp2   ! Michael Falk, for GABLS2 case

!$omp   threadprivate(iwp2, irtp2, ithlp2, irtpthlp, iwprtp, iwpthlp)
!$omp   threadprivate(iwp4, iwpthvp, irtpthvp, ithlpthvp, itaum, iKhm)
!$omp   threadprivate(iwprcp, ithlprcp, irtprcp, ircp2, iupwp, ivpwp)
!$omp   threadprivate(irhom, iscm, iem, ishear, iFrad, iFrad_LW)
!$omp   threadprivate(iFrad_SW, iFprec, iFcsed, iup2_plus_vp2)

        ! Sedimentation velocities
        integer ::
     .    iVrr,     ! Brian
     .    iVNr,     !  " "
     .    iVsnow,   ! COAMPS
     .    iVice,    !  " "
     .    iVgraupel !  " " 

!$omp   threadprivate(iVrr, iVNr, iVsnow, iVice, iVgraupel)

        integer ::
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

!$omp   threadprivate(iwp2_bt, iwp2_ma, iwp2_ta, iwp2_ac, iwp2_bp)
!$omp   threadprivate(iwp2_pr1, iwp2_pr2, iwp2_pr3, iwp2_dp1, iwp2_dp2)
!$omp   threadprivate(iwp2_cl)

        integer ::
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

!$omp   threadprivate(iwprtp_bt, iwprtp_ma, iwprtp_ta, iwprtp_tp)
!$omp   threadprivate(iwprtp_ac, iwprtp_bp, iwprtp_pr1, iwprtp_pr2)
!$omp   threadprivate(iwprtp_pr3, iwprtp_dp1)

        integer ::
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

!$omp   threadprivate(iwpthlp_bt, iwpthlp_ma, iwpthlp_ta, iwpthlp_tp)
!$omp   threadprivate(iwpthlp_ac, iwpthlp_bp, iwpthlp_pr1, iwpthlp_pr2)
!$omp   threadprivate(iwpthlp_pr3, iwpthlp_dp1)

!    Dr. Golaz's new variance budget terms
!    qt was changed to rt to avoid confusion

        integer ::
     .     irtp2_bt,
     .     irtp2_ma,
     .     irtp2_ta,
     .     irtp2_tp,
     .     irtp2_dp1,
     .     irtp2_dp2,
     .     irtp2_cl
!$omp   threadprivate(irtp2_bt, irtp2_ma, irtp2_ta, irtp2_tp)
!$omp   threadprivate(irtp2_dp1, irtp2_dp2, irtp2_cl)

        integer ::
     .     ithlp2_bt,
     .     ithlp2_ma,
     .     ithlp2_ta,
     .     ithlp2_tp,
     .     ithlp2_dp1,
     .     ithlp2_dp2,
     .     ithlp2_cl

!$omp   threadprivate(ithlp2_bt, ithlp2_ma, ithlp2_ta, ithlp2_tp)
!$omp   threadprivate(ithlp2_dp1, ithlp2_dp2, ithlp2_cl)

        integer ::
     .    irtpthlp_bt,
     .    irtpthlp_ma,
     .    irtpthlp_ta,
     .    irtpthlp_tp1,
     .    irtpthlp_tp2,
     .    irtpthlp_dp1,
     .    irtpthlp_dp2

!$omp   threadprivate(irtpthlp_bt, irtpthlp_ma, irtpthlp_ta)
!$omp   threadprivate(irtpthlp_tp1, irtpthlp_tp2, irtpthlp_dp1)
!$omp   threadprivate(irtpthlp_dp2, irtpthlp_cl)

        integer ::
     .    iup2,
     .    ivp2

!$omp   threadprivate(iup2, ivp2)

        integer ::
     .    iup2_bt,
     .    iup2_ta,
     .    iup2_tp,
     .    iup2_ma,
     .    iup2_dp1,
     .    iup2_dp2,
     .    iup2_pr1,
     .    iup2_pr2,
     .    ivp2_bt,
     .    ivp2_ta,
     .    ivp2_tp,
     .    ivp2_ma,
     .    ivp2_dp1,
     .    ivp2_dp2,
     .    ivp2_pr1,
     .    ivp2_pr2

!$omp   threadprivate(iup2_bt, iup2_ta, iup2_tp, iup2_ma, iup2_dp1)
!$omp   threadprivate(iup2_dp2, iup2_pr1, iup2_pr2)
!$omp   threadprivate(ivp2_bt, ivp2_ta, ivp2_tp, ivp2_ma, ivp2_dp1)
!$omp   threadprivate(ivp2_dp2, ivp2_pr1, ivp2_pr2)

#ifdef SCALARS
!       Passive scalars.  Note that floating point roundoff may make 
!       mathematically equivalent variables different values.
        integer ::
     .    isclraprtp,          ! sclr'(1)rt'     / rt'^2
     .    isclrbprtp,          ! sclr'(2)rt'     / thl'rt'
     .    isclrap2,            ! sclr'(1)^2      / rt'^2
     .    isclrbp2,            ! sclr'(2)^2      / thl'^2
     .    isclrapthvp,         ! sclr'(1)th_v'   / rt'th_v' 
     .    isclrbpthvp,         ! sclr'(2)th_v'   / th_l' th_v'
     .    isclrapthlp,         ! sclr'(1)th_l'   / rt'th_l' 
     .    isclrbpthlp,         ! sclr'(2)th_l'   / th_l'^2
     .    isclraprcp,          ! sclr'(1)rc'     / rt'rc'
     .    isclrbprcp,          ! sclr'(2)rc'     / th_l'rc'
     .    iwpsclrap,           ! w'slcr'(1)      / w'rt'
     .    iwpsclrbp,           ! w'sclr'(2)      / w'th_l'
     .    iwp2sclrap,          ! w'^2 sclr'(1)   / w'^2 rt'
     .    iwp2sclrbp,          ! w'^2 sclr'(2)   / w'^2 th_l'
     .    iwpsclrap2,          ! w'sclr'(1)^2    / w'rt'^2
     .    iwpsclrbp2,          ! w'sclr'(2)^2    / w'th_l'^2
     .    iwpsclraprtp,        ! w'sclr'(1)rt'   / w'rt'^2
     .    iwpsclrbprtp,        ! w'sclr'(2)rt'   / w'thl'rt'
     .    iwpsclrapthlp,       ! w'sclr'(1)th_l' / w'rt'th_l' 
     .    iwpsclrbpthlp        ! w'sclr'(2)th_l' / w'th_l'^2

        integer
     .     iwpedsclrap,  ! eddy sclr'(1)w'
     .     iwpedsclrbp   ! eddy sclr'(2)w'
#endif /*SCALARS*/

        ! Indices for statistics in sfc file

        integer ::
     .    iustar,
     .    ilh,
     .    ish,
     .    icc,
     .    ilwp,
     .    izb,
     .    izi,
     .    irain,   ! Brian
     .    ipflux,  ! Brian
     .    irrm_sfc ! Brian

        integer ::
     .    iwp23_cn,
     .    irtm_cn,
     .    ithlm_cn,
     .    irtp2_cn,
     .    ithlp2_cn,
     .    irtpthlp_cn,
     .    iup2_cn,
     .    ivp2_cn

!$omp   threadprivate(iustar, ilh, ish, icc, ilwp, izb, izi)
!$omp   threadprivate(irain, ipflux, irrm_sfc)
!$omp   threadprivate(iwp23_cn, irtm_cn, ithlm_cn)
!$omp   threadprivate(irtp2_cn, ithlp2_cn, irtpthlp_cn)
!$omp   threadprivate(iup2_cn, ivp2_cn)

        ! Variables that contains all the statistics

        type (stats), target :: zt,  ! zt grid
     .                          zm,  ! zm grid
     .                          sfc  ! sfc

!$omp   threadprivate(zt, zm, sfc)

        ! Scratch space

        real, allocatable :: ztscr01(:), ztscr02(:), ztscr03(:), 
     .                       ztscr04(:), ztscr05(:), ztscr06(:), 
     .                       ztscr07(:), ztscr08(:), ztscr09(:),
     .                       ztscr10(:), ztscr11(:), ztscr12(:),
     .                       ztscr13(:), ztscr14(:), ztscr15(:)

!$omp   threadprivate(ztscr01, ztscr02, ztscr03, ztscr04, ztscr05)
!$omp   threadprivate(ztscr06, ztscr07, ztscr08, ztscr09, ztscr10)
!$omp   threadprivate(ztscr11, ztscr12, ztscr13, ztscr14, ztscr15)

        real, allocatable :: zmscr01(:), zmscr02(:), zmscr03(:), 
     .                       zmscr04(:), zmscr05(:), zmscr06(:),
     .                       zmscr07(:), zmscr08(:), zmscr09(:),
     .                       zmscr10(:), zmscr11(:), zmscr12(:),
     .                       zmscr13(:)

!$omp   threadprivate(zmscr01, zmscr02, zmscr03, zmscr04, zmscr05) 
!$omp   threadprivate(zmscr06, zmscr07, zmscr08, zmscr09, zmscr10)
!$omp   threadprivate(zmscr11, zmscr12, zmscr13)

      end module hoc_stats

#endif /*STATS*/
