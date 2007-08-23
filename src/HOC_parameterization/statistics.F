!-----------------------------------------------------------------------
!  $Id: statistics.F,v 1.29 2007-08-23 16:06:38 faschinj Exp $

#ifdef STATS
      module stats_type

!     Description:
!     Contains derived data type 'stats'
!-----------------------------------------------------------------------

        use outputfile_class

        ! Derived data types to store GrADS statistics

        type stats 

          ! Number of fields to sample
          integer :: nn

          ! Vertical extent of variable
          integer :: kk

          ! Vertical levels
          real, pointer, dimension(:) :: z

!         Array to store sampled fields

          real, pointer, dimension(:,:) :: x

          integer, pointer, dimension(:,:) :: n

!         Data for GrADS output

          type (outputfile) f

        end type stats
	contains
!-----------------------------------------------------------------------
        subroutine stat_assign( location, name, 
     .                     description, units, grid )

!     Description: Assigns pointers for statistics variables in grid
!    
!-----------------------------------------------------------------------
           !input
             integer,intent(in) :: location !
             character(len = *), intent(in) :: name        ! Variable name
             character(len = *), intent(in) :: description ! Variable description
             character(len = *), intent(in) :: units       ! Variable units

           !output
             type(stats), intent(out) :: grid

             grid%f%var(location)%ptr => grid%x(:,location)
             grid%f%var(location)%name = name
             grid%f%var(location)%description = description
             grid%f%var(location)%units = units
	
             !Example of the old format 
             !changed by Joshua Fasching 23 August 2007
             
             !zt%f%var(ithlm)%ptr => zt%x(:,k)
             !zt%f%var(ithlm)%name = "thlm"
             !zt%f%var(ithlm)%description = "thetal (K)"
             !zt%f%var(ithlm)%units = "K"
      
             return
           end subroutine stat_assign
      end module stats_type

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

!-----------------------------------------------------------------------
      subroutine stats_init( iunit, fnamelist, nnzp, gzt, gzm,
     .                       day, month, year, rlat, rlon, time, delt )

!     Description: Initializes the statistics saving functionality of
!     the HOC model.
!-----------------------------------------------------------------------
      use hoc_stats
      use outputgrads
      use outputnc

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
     .  rlat, rlon,! Latitude and Longitude             [Degrees N/E]
     .  time,      ! Model time                         [s]
     .  delt       ! Timestep (dtmain in HOC)           [s]

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
     .                   time+stats_tout, stats_tout,
     .                   zt%nn )

      else ! Open NetCDF file
        call open_nc( zt%f, iunit, fdir, fname, 
     .                1, zt%kk, zt%z,
     .                day, month, year, rlat, rlon,
     .                time+stats_tout, stats_tout,
     .                zt%nn )

      end if

      ! Default initialization for array indices for zt

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
      iexner        = 0
      iLscale       = 0
      iwp3          = 0
      iwpthlp2      = 0
      iwp2thlp      = 0
      iwprtp2       = 0
      iwp2rtp       = 0
      ilup          = 0
      ildown        = 0
      itaut         = 0
      iKht          = 0
      iwp2thvp      = 0
      iwp2rcp       = 0
      iwprtpthlp    = 0
      isct          = 0
      irhot         = 0
      iNcm          = 0  ! Brian
      iNcnm         = 0
      iNim          = 0
      ised_rcm      = 0  ! Brian
      irsm          = 0  ! Brian
      irrm          = 0  ! Brian
      iNrm          = 0  ! Brian
      irain_rate    = 0  ! Brian
      iAKm          = 0  ! analytic Kessler.  Vince Larson 22 May 2005
      iAKm_est      = 0  ! LH Kessler.  Vince Larson 22 May 2005
      iradht        = 0
      iradht_LW     = 0
      iradht_SW     = 0
      idiam         = 0
      im            = 0
      idqc_dt_icedfs= 0
      iu_T_cm       = 0

      imean_vol_rad_rain  = 0  ! Brian
      imean_vol_rad_cloud = 0

      irsnowm       = 0
      irgraupelm    = 0
      iricem        = 0

      irtm_bt       = 0
      irtm_ma       = 0
      irtm_ta       = 0
      irtm_forcing  = 0
      irtm_mc       = 0
      ithlm_bt      = 0
      ithlm_ma      = 0
      ithlm_ta      = 0
      ithlm_forcing = 0
      ithlm_mc      = 0

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
      irrm_cond_adj = 0
      irrm_mc       = 0
      irrm_cl       = 0

      iNrm_bt       = 0
      iNrm_ma       = 0
      iNrm_sd       = 0
      iNrm_dff      = 0
      iNrm_cond     = 0
      iNrm_auto     = 0
      iNrm_cond_adj = 0
      iNrm_mc       = 0
      iNrm_cl       = 0

      irsnowm_bt    = 0
      irsnowm_ma    = 0
      irsnowm_sd    = 0
      irsnowm_dff   = 0
      irsnowm_mc    = 0
      irsnowm_cl    = 0

      irgraupelm_bt = 0
      irgraupelm_ma = 0
      irgraupelm_sd = 0
      irgraupelm_dff= 0
      irgraupelm_mc = 0
      irgraupelm_cl = 0

      iricem_bt     = 0
      iricem_ma     = 0
      iricem_sd     = 0
      iricem_dff    = 0
      iricem_mc     = 0
      iricem_cl     = 0

      ivm_bt = 0
      ivm_ma = 0
      ivm_gf = 0
      ivm_cf = 0
      ivm_ta = 0

      ium_bt = 0
      ium_ma = 0
      ium_gf = 0
      ium_cf = 0
      ium_ta = 0

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
      isclram     = 0
      isclram_f   = 0
      isclrbm     = 0
      isclrbm_f   = 0

      iedsclram   = 0
      iedsclrbm   = 0
#endif /*SCALARS*/

!     Assign pointers for statistics variables zt

      k = 1
      do i=1,zt%nn

        select case ( trim(vars_zt(i)) )
        case ('thlm')
          ithlm = k
          call stat_assign( ithlm, "thlm", 
     .          "thetal (K)", "K", zt)
          k = k + 1

        case ('thvm')
          ithvm = k
          call stat_assign( ithvm, "thvm",
     .          "virtual potential temperature (K)","K",zt)
          k = k + 1

        case ('rtm')
          irtm = k

          call stat_assign( irtm, "rtm",
     .          "total water mixing ratio (kg/kg)","kg/kg",zt)
          
	  !zt%f%var(irtm)%ptr => zt%x(:,k)
          !zt%f%var(irtm)%name = "rtm"
          !zt%f%var(irtm)%description 
     .    != "total water mixing ratio (kg/kg)"
          !zt%f%var(irtm)%units = "kg/kg"
          
	  k = k + 1

        case ('rcm')
          ircm = k
          call stat_assign( ircm, "rcm",
     .          "liquid water mixing ratio (kg/kg)","kg/kg",zt)
          k = k + 1
        case ('um')
          ium = k
          call stat_assign( ium, "um",
     .          "u wind (m/s)","m/s",zt)
          k = k + 1
        case ('vm')
          ivm = k
          call stat_assign( ivm,  "vm",
     .          "v wind (m/s)","m/s",zt)
          k = k + 1
        case ('wmt')
          iwmt = k
          call stat_assign( iwmt, "wm",
     .          "w wind (m/s)","m/s",zt)
          k = k + 1
        case ('ug')
          iug = k
          call stat_assign( iug,"ug",
     .         "u geostrophic wind (m/s)", "m/s", zt)
          k = k + 1
        case ('vg')
          ivg = k
          call stat_assign(ivg,"vg",
     .         "v geostrophic wind (m/s)", "m/s",zt)
          k = k + 1
        case ('cf')
          icf = k
          call stat_assign(icf,"cf",
     .         "cloud fraction", "count",zt)
          k = k + 1
        case ('p')
          ip = k
          call stat_assign(ip,"p",
     .         "pressure (Pa)","Pa",zt)
          k = k + 1
        case ('exner')
          iexner = k
          call stat_assign(iexner,"exner",
     .         "Exner","count",zt)
          k = k + 1
        case ('Lscale')
          iLscale = k
          call stat_assign(iLscale,"Lscale",
     .         "Mixing length","m",zt)
          k = k + 1
        case ('thlm_forcing')
          ithlm_forcing = k
          call stat_assign(ithlm_forcing,"thlm_f",
     .         "thetal forcing", "K/s",zt)
          k = k + 1
        case ('thlm_mc')
          ithlm_mc = k
          call stat_assign(ithlm_mc,"thlm_mc",
     .         "thetal micro", "K/s",zt)
          k = k + 1
        case ('rtm_forcing')
          irtm_forcing = k
	  call stat_assign(irtm_forcing,"rtm_f",
     .         "rt forcing", "kg/(kg s)",zt)
          k = k + 1
        case ('rtm_mc')
          irtm_mc = k
          call stat_assign(irtm_mc,"rtm_mc",
     .         "rt micro", "kg/(kg s)",zt)
          k = k + 1

        case ('wp3')
          iwp3 = k
          call stat_assign(iwp3,"wp3",
     .         "w third order moment", "(m^3)/(s^3)",zt)
          k = k + 1

        case ('wpthlp2')
          iwpthlp2 = k
          call stat_assign(iwpthlp2,"wpthlp2",
     .         "wpthlp2 covariance", "(m K^2)/s",zt)
          k = k + 1

        case ('wp2thlp')
          iwp2thlp = k
          call stat_assign(iwp2thlp,"wp2thlp",
     .         "wp2thlp covariance", "(m^2 K)/(s^2)",zt)
          k = k + 1

        case ('wprtp2')
          iwprtp2 = k
          call stat_assign(iwprtp2,"wprtp2",
     .         "wprtp2 covariance", "(m kg^2)/(s kg^2)",zt)
          k = k + 1

        case ('wp2rtp')
          iwp2rtp = k
          call stat_assign(iwp2rtp,"wp2rtp",
     .         "wp2rtp covariance", "(m2 kg)/(s2 kg)",zt)
          k = k + 1

        case ('lup')
          ilup = k
          call stat_assign(ilup,"lup",
     .         "Upward mixing length","m",zt)
          k = k + 1

        case ('ldown')
          ildown = k

          call stat_assign(ildown,"ldown",
     .         "Downward mixing length","m",zt)
          k = k + 1

        case ('taut')
          itaut = k
          call stat_assign(itaut,"taut",
     .         "Dissipation time","s",zt)
          k = k + 1

        case ('kht')
          iKht = k
          call stat_assign(iKht,"Kht",
     .         "Eddy diffusivity","m^2/s",zt)
          k = k + 1

        case ('wp2thvp')
          iwp2thvp = k
          call stat_assign(iwp2thvp,"wp2thvp",
     .         "wp2thvp","(m2 K)/s2",zt)
          k = k + 1

        case ('wp2rcp')
          iwp2rcp = k
          call stat_assign(iwp2rcp,"wp2rcp",
     .         "wp2rcp","(m2 kg)/(s2 kg)",zt)
          k = k + 1

        case ('wprtpthlp')
          iwprtpthlp = k
          call stat_assign(iwprtpthlp,"wprtpthlp",
     .         "wprtpthlp","(m kg K)/(s kg)",zt)
          k = k + 1

        case ('sc')
          isct = k
          call stat_assign(isct,"Sc",
     .         "Sc","count",zt)
          k = k + 1

        case ('rhot')
          irhot = k
          call stat_assign(irhot,"rhot",
     .         "density","kg/m^3",zt)
          k = k + 1

        case ('Ncm')           ! Brian
          iNcm = k
          call stat_assign(iNcm,"Ncm",
     .         "Cloud droplet number concentration (num/m^3)",
     .         "count/m^3",zt)
          k = k + 1

        case ('Ncnm')
          iNcnm = k
          call stat_assign(iNcnm,"Ncnm",
     .         "Cloud nuclei number concentration (num/m^3)",
     .         "count/m^3",zt)
          k = k + 1

        case ('Nim')           ! Brian
          iNim = k
          call stat_assign(iNim,"Nim",
     .         "Ice crystal number concentration (num/m^3)",
     .         "count/m^3",zt)
          k = k + 1

        case ('sed_rcm')       ! Brian
          ised_rcm = k
          call stat_assign(ised_rcm,"sed_rcm",
     .         "d(rcm)/dt due to cloud sedimentation (kg / [m^2 s])",
     .         "kg/(m^2 s)",zt)
          k = k + 1

        case ('rsm')           ! Brian
          irsm = k
          call stat_assign(irsm,"rsm",
     .         "Saturation mixing ratio (kg/kg)","kg/kg",zt)
          k = k + 1

        case ('rrm')           ! Brian
          irrm = k
          call stat_assign(irrm,"rrm",
     .         "Rain water mixing ratio (kg/kg)","kg/kg",zt)
          k = k + 1

        case ('rsnowm')
          irsnowm = k
          call stat_assign(irsnowm,"rsnowm",
     .         "Snow water mixing ratio (kg/kg)","kg/kg",zt)
          k = k + 1

        case ('ricem')
          iricem = k
          call stat_assign(iricem,"ricem",
     .         "Pristine ice water mixing ratio (kg/kg)","kg/kg",zt)
          k = k + 1

        case ('rgraupelm')
          irgraupelm = k
          call stat_assign(irgraupelm,"rgraupelm",
     .         "Graupel water mixing ratio (kg/kg)","kg/kg",zt)
          k = k + 1

        case ('Nrm')           ! Brian
          iNrm = k
          call stat_assign(iNrm,"Nrm",
     .         "Rain droplet number concentration (num/m^3)",
     .         "count/(m^3)",zt)
          k = k + 1

        case ('mean_vol_rad_rain')  ! Brian
          imean_vol_rad_rain = k
          call stat_assign(imean_vol_rad_rain,"mvrr",
     .         "Rain drop mean volume radius (m)","m",zt)
          k = k + 1

        case ('mean_vol_rad_cloud')
          imean_vol_rad_cloud = k

          call stat_assign(imean_vol_rad_cloud,"mvrc",
     .         "Cloud drop mean volume radius (m)","m",zt)
          k = k + 1

        case ('rain_rate')     ! Brian
          irain_rate = k

          call stat_assign(irain_rate,"rain_rate",
     .         "Rain rate (mm/day)","mm/day",zt)
          k = k + 1
 
        case ('AKm')           ! Vince Larson 22 May 2005
          iAKm = k
          call stat_assign(iAKm,"AKm",
     .         "Analytic Kessler ac [kg/kg]","kg/kg",zt)
          k = k + 1
 
        case ('AKm_est')       ! Vince Larson 22 May 2005
          iAKm_est = k

          call stat_assign(iAKm_est,"AKm_est",
     .         "LH Kessler estimate [kg/kg]","kg/kg",zt)
          k = k + 1

        case ('radht')
          iradht = k

          call stat_assign(iradht,"radht",
     .         "Heating rate","K/s",zt)
          k = k + 1

        case ('radht_LW')
          iradht_LW = k

          call stat_assign(iradht_LW,"radht_LW",
     .         "Long-wave heating rate","K/s",zt)

          k = k + 1

        case ('radht_SW')
          iradht_SW = k
          call stat_assign(iradht_SW,"radht_SW",
     .         "Short-wave heating rate","K/s",zt)
          k = k + 1
        case ('diam')
          idiam = k

          call stat_assign(idiam,"diam",
     .         "Ice crystal diameter","m",zt)
          k = k + 1

        case ('m')
          im = k
          call stat_assign(im,"m",
     .         "Ice crystal mass","kg",zt)
          k = k + 1

        case ('dqc_dt_icedfs')

          idqc_dt_icedfs = k
          call stat_assign(idqc_dt_icedfs,"dqc_dt_icedfs",
     .         "Change in liquid due to ice","kg/kg/s",zt)
          k = k + 1

        case ('u_T_cm')
          iu_T_cm = k
          call stat_assign(iu_T_cm,"u_T_cm",
     .         "Ice crystal fallspeed","cm/s",zt)
          k = k + 1

        case ('rtm_bt')
          irtm_bt = k

          call stat_assign(irtm_bt,"rtm_bt",
     .         "rtm budget","kg/kg/s",zt)
          k = k + 1

        case ('rtm_ma')
          irtm_ma = k

          call stat_assign(irtm_ma,"rtm_ma",
     .         "rtm ma","kg/kg/s",zt)
          k = k + 1

        case ('rtm_ta')
          irtm_ta = k

          call stat_assign(irtm_ta,"rtm_ta",
     .         "rtm ta","kg/kg/s",zt)
          k = k + 1

        case ('thlm_bt')
          ithlm_bt = k

          call stat_assign(ithlm_bt,"thlm_bt",
     .         "thlm bt","kg/kg/s",zt)
          k = k + 1

        case ('thlm_ma')
          ithlm_ma = k

          call stat_assign(ithlm_ma,"thlm_ma",
     .         "thlm ma","kg/kg/s",zt)
          k = k + 1

        case ('thlm_ta')
          ithlm_ta = k

          call stat_assign(ithlm_ta,"thlm_ta",
     .         "thlm ta","kg/kg/s",zt)
          k = k + 1
 
        case ('wp3_bt')
          iwp3_bt = k

          call stat_assign(iwp3_bt,"wp3_bt",
     .         "wp3 budget","(m^3)/(s^4)",zt)
          k = k + 1
 
        case ('wp3_ma')
          iwp3_ma = k

          call stat_assign(iwp3_ma,"wp3_ma",
     .         "wp3 mean advection","(m^3)/(s^4)",zt)
          k = k + 1
 
        case ('wp3_ta')
          iwp3_ta = k

          call stat_assign(iwp3_ta,"wp3_ta",
     .         "wp3 turbulent advection","(m^3)/(s^4)",zt)

          k = k + 1
 
        case ('wp3_tp')
          iwp3_tp = k
          call stat_assign(iwp3_tp,"wp3_tp",
     .         "wp3 turbulent transport","(m^3)/(s^4)",zt)
          k = k + 1
 
        case ('wp3_ac')
          iwp3_ac = k
          call stat_assign(iwp3_ac,"wp3_ac",
     .         "wp3 accumulation term","(m^3)/(s^4)",zt)
          k = k + 1
 
        case ('wp3_bp')
          iwp3_bp = k
          call stat_assign(iwp3_bp,"wp3_bp",
     .         "wp3 buoyancy production","(m^3)/(s^4)",zt)
          k = k + 1
 
        case ('wp3_pr1')
          iwp3_pr1 = k
          call stat_assign(iwp3_pr1,"wp3_pr1",
     .         "wp3 pressure term 1","(m^3)/(s^4)",zt)
          k = k + 1
 
        case ('wp3_pr2')
          iwp3_pr2 = k
          call stat_assign(iwp3_pr2,"wp3_pr2",
     .         "wp3 pressure term 2","(m^3)/(s^4)",zt)

          k = k + 1

        case ('wp3_dp1')
          iwp3_dp1 = k
          call stat_assign(iwp3_dp1,"wp3_dp1",
     .         "wp3 dissipation term 1","(m^3)/(s^4)",zt)
          k = k + 1
 
        case ('wp3_cl')
          iwp3_cl = k
          call stat_assign(iwp3_cl,"wp3_cl",
     .         "wp3 clipping term","(m^3)/(s^4)",zt)
          k = k + 1
 
        case ('rrm_bt')
          irrm_bt = k
          call stat_assign(irrm_bt,"rrm_bt",
     .         "rrm budget","(kg/kg)/(s)",zt)
          k = k + 1
 
        case ('rrm_ma')
          irrm_ma = k

          call stat_assign(irrm_ma,"rrm_ma",
     .         "rrm mean advection","(kg/kg)/(s)",zt)
          k = k + 1
 
        case ('rrm_sd')
          irrm_sd = k

          call stat_assign(irrm_sd,"rrm_sd",
     .         "rrm sedimentation","(kg/kg)/(s)",zt)
          k = k + 1
 
        case ('rrm_dff')
          irrm_dff = k

          call stat_assign(irrm_dff,"rrm_dff",
     .         "rrm diffusion","(kg/kg)/(s)",zt)
          k = k + 1
 
        case ('rrm_cond')
          irrm_cond = k

          call stat_assign(irrm_cond,"rrm_cond",
     .         "rrm cond/evap","(kg/kg)/(s)",zt)
          k = k + 1
 
        case ('rrm_auto')
          irrm_auto = k

          call stat_assign(irrm_auto,"rrm_auto",
     .         "rrm autoconversion","(kg/kg)/(s)",zt)
          k = k + 1
 
        case ('rrm_accr')
          irrm_accr = k
          call stat_assign(irrm_accr,"rrm_accr",
     .         "rrm accretion","(kg/kg)/(s)",zt)
          k = k + 1

        case ('rrm_cond_adj')
          irrm_cond_adj = k

          call stat_assign(irrm_cond_adj,"rrm_cond_adj",
     .         "rrm cond/evap adjustment due to over-evaporation",
     .         "(kg/kg)/(s)",zt)
          k = k + 1
 
        case ('rrm_mc')
          irrm_mc = k

          call stat_assign(irrm_mc,"rrm_mc",
     .         "rrm total microphysical tendency","(kg/kg)/(s)",zt)

          k = k + 1

        case ('rrm_cl')
          irrm_cl = k
          call stat_assign(irrm_cl,"rrm_cl",
     .         "rrm clipping term","(kg/kg)/(s)",zt)

          k = k + 1
 
        case ('Nrm_bt')
          iNrm_bt = k
          call stat_assign(iNrm_bt,"Nrm_bt",
     .         "Nrm budget","(count/m^3)/s",zt)

          k = k + 1
 
        case ('Nrm_ma')
          iNrm_ma = k

          call stat_assign(iNrm_ma,"Nrm_ma",
     .         "Nrm mean advection","(count/m^3)/s",zt)
          k = k + 1
 
        case ('Nrm_sd')
          iNrm_sd = k

          call stat_assign(iNrm_sd,"Nrm_sd",
     .         "Nrm sedimentation","(count/m^3)/s",zt)

          k = k + 1
 
        case ('Nrm_dff')
          iNrm_dff = k
          call stat_assign(iNrm_dff,"Nrm_dff",
     .         "Nrm diffusion","(count/m^3)/s",zt)

          k = k + 1
 
        case ('Nrm_cond')
          iNrm_cond = k

          call stat_assign(iNrm_cond,"Nrm_cond",
     .         "Nrm cond/evap","(count/m^3)/s",zt)
          k = k + 1
 
        case ('Nrm_auto')
          iNrm_auto = k

          call stat_assign(iNrm_auto,"Nrm_auto",
     .         "Nrm autoconversion","(count/m^3)/s",zt)

          k = k + 1

        case ('Nrm_cond_adj')
          iNrm_cond_adj = k

          call stat_assign(iNrm_cond_adj,"Nrm_cond_adj",
     .         "Nrm cond/evap adjustment due to over-evaporation",
     .         "(count/m^3)/s",zt)
          k = k + 1
 
        case ('Nrm_mc')
          iNrm_mc = k
          call stat_assign(iNrm_mc,"Nrm_mc",
     .         "Nrm micro","(count/m^3)/s",zt)

          k = k + 1

        case ('Nrm_cl')
          iNrm_cl = k
          call stat_assign(iNrm_cl,"Nrm_cl",
     .         "Nrm clipping term","(count/m^3)/s",zt)
          k = k + 1

        case ('rsnowm_bt')
          irsnowm_bt = k
          call stat_assign(irsnowm_bt,"rsnowm_bt",
     .         "rsnowm budget","(kg/kg)/s",zt)

          k = k + 1
 
        case ('rsnowm_ma')
          irsnowm_ma = k

          call stat_assign(irsnowm_ma,"rsnowm_ma",
     .         "rsnowm mean advection","(kg/kg)/s",zt)
          k = k + 1
 
        case ('rsnowm_sd')
          irsnowm_sd = k
          call stat_assign(irsnowm_sd,"rsnowm_sd",
     .         "rsnowm sedimentation","(kg/kg)/s",zt)
          k = k + 1
 
        case ('rsnowm_dff')
          irsnowm_dff = k

          call stat_assign(irsnowm_dff,"rsnowm_dff",
     .         "rsnowm diffusion","(kg/kg)/s",zt)
          k = k + 1

        case ('rsnowm_mc')
          irsnowm_mc = k

          call stat_assign(irsnowm_mc,"rsnowm_mc",
     .         "rsnowm micro","(kg/kg)/s",zt)
          k = k + 1

        case ('rsnowm_cl')
          irsnowm_cl = k

          call stat_assign(irsnowm_cl,"rsnowm_cl",
     .         "rsnowm clipping term","(kg/kg)/s",zt)
          k = k + 1

       case ('ricem_bt')
          iricem_bt = k

          call stat_assign(iricem_bt,"ricem_bt",
     .         "ricem budget","(kg/kg)/s",zt)
 
          k = k + 1
 
        case ('ricem_ma')
          iricem_ma = k

          call stat_assign(iricem_ma,"ricem_ma",
     .         "ricem mean advection","(kg/kg)/s",zt)
          k = k + 1
 
        case ('ricem_sd')
          iricem_sd = k

          call stat_assign(iricem_sd,"ricem_sd",
     .         "ricem sedimentation","(kg/kg)/s",zt)
          k = k + 1
 
        case ('ricem_dff')
          iricem_dff = k

          call stat_assign(iricem_dff,"ricem_dff",
     .         "ricem diffusion","(kg/kg)/s",zt)
          k = k + 1

        case ('ricem_mc')
          iricem_mc = k

          call stat_assign(iricem_mc,"ricem_mc",
     .         "ricem micro","(kg/kg)/s",zt)
          k = k + 1

        case ('ricem_cl')
          iricem_cl = k

          call stat_assign(iricem_cl,"ricem_cl",
     .         "ricem clipping term","(kg/kg)/s",zt)
          k = k + 1

        case ('rgraupelm_bt')
          irgraupelm_bt = k

          call stat_assign(irgraupelm_bt,"rgraupelm_bt",
     .         "rgraupelm budget","(kg/kg)/s",zt)
          k = k + 1
 
        case ('rgraupelm_ma')
          irgraupelm_ma = k

          call stat_assign(irgraupelm_ma,"rgraupelm_ma",
     .         "rgraupelm mean advection","(kg/kg)/s",zt)
          k = k + 1
 
        case ('rgraupelm_sd')
          irgraupelm_sd = k

          call stat_assign(irgraupelm_sd,"rgraupelm_sd",
     .         "rgraupelm sedimentation","(kg/kg)/s",zt)
          k = k + 1
 
        case ('rgraupelm_dff')
          irgraupelm_dff = k

          call stat_assign(irgraupelm_dff,"rgraupelm_dff",
     .         "rgraupelm diffusion","(kg/kg)/s",zt)
          k = k + 1

        case ('rgraupelm_mc')
          irgraupelm_mc = k

          call stat_assign(irgraupelm_mc,"rgraupelm_mc",
     .         "rgraupelm micro","(kg/kg)/s",zt)
          k = k + 1

        case ('rgraupelm_cl')
          irgraupelm_cl = k

          call stat_assign(irgraupelm_cl,"rgraupelm_cl",
     .         "rgraupelm clipping term","(kg/kg)/s",zt)
          k = k + 1

        case ('vm_bt')
          ivm_bt = k

          call stat_assign(ivm_bt,"vm_bt",
     .         "vm time tendency(note: does not balance at zt=2)",
     .         "m/s",zt)
          k = k + 1

        case ('vm_ma')
          ivm_ma = k
          call stat_assign(ivm_ma,"vm_ma",
     .         "vm mean advection","m/s",zt)
          k = k + 1

        case ('vm_gf')
          ivm_gf = k

          call stat_assign(ivm_gf,"vm_gf",
     .         "vm geostrophic forcing","m/s",zt)
          k = k + 1

        case ('vm_cf')
          ivm_cf = k

          call stat_assign(ivm_cf,"vm_cf",
     .         "vm coriolis forcing","m/s",zt)
          k = k + 1

        case ('vm_ta')
          ivm_ta = k

          call stat_assign(ivm_ta,"vm_ta",
     .         "vm turbulent transport","m/s",zt)
          k = k + 1

        case ('um_bt')
          ium_bt = k

          call stat_assign(ium_bt,"um_bt",
     .         "um time tendency (note: does not balance at zt=2)",
     .         "m/s",zt)
          k = k + 1

        case ('um_ma')
          ium_ma = k

          call stat_assign(ium_ma,"um_ma",
     .         "um mean advection","m/s",zt)
          k = k + 1

        case ('um_gf')
          ium_gf = k
          call stat_assign(ium_gf,"um_gf",
     .         "um geostrophic forcing","m/s",zt)
          k = k + 1

        case ('um_cf')
          ium_cf = k
          call stat_assign(ium_cf,"um_cf",
     .         "um coriolis forcing","m/s",zt)
          k = k + 1

        case ('um_ta')
          ium_ta = k
          call stat_assign(ium_ta,"um_ta",
     .         "um turbulent transport","m/s",zt)
          k = k + 1
 
        case ('a')
          ia = k
          call stat_assign(ia,"a",
     .         "pdf parameter a","count",zt)
          k = k + 1
 
        case ('w1')
          iw1 = k
          call stat_assign(iw1,"w1",
     .         "pdf parameter w1","m/s",zt)

          k = k + 1
 
        case ('w2')
          iw2 = k

          call stat_assign(iw2,"w2",
     .         "pdf parameter w2","m/s",zt)
          k = k + 1
 
        case ('sw1')
          isw1 = k
          call stat_assign(isw1,"sw1",
     .         "pdf parameter sw1","m^2/s^2",zt)

          k = k + 1
 
        case ('sw2')
          isw2 = k

          call stat_assign(isw2,"sw2",
     .         "pdf parameter sw2","m^2/s^2",zt)
          k = k + 1
 
        case ('thl1')
          ithl1 = k

          call stat_assign(ithl1,"thl1",
     .         "pdf parameter thl1","K",zt)

          k = k + 1
 
        case ('thl2')
          ithl2 = k

          call stat_assign(ithl2,"thl2",
     .         "pdf parameter thl2","K",zt)
          k = k + 1
 
        case ('sthl1')
          isthl1 = k

          call stat_assign(isthl1,"sthl1",
     .         "pdf parameter sthl1","K^2",zt)

          k = k + 1
 
        case ('sthl2')
          isthl2 = k
          call stat_assign(isthl2,"sthl2",
     .         "pdf parameter sthl2","K^2",zt)

          k = k + 1
 
        case ('rt1')
          irt1 = k
          call stat_assign(irt1,"rt1",
     .         "pdf parameter rt1","kg/kg",zt)

          k = k + 1
 
        case ('rt2')
          irt2 = k

          call stat_assign(irt2,"rt2",
     .         "pdf parameter rt2","kg/kg",zt)
          k = k + 1
 
        case ('srt1')
          isrt1 = k
          call stat_assign(isrt1,"srt1",
     .         "pdf parameter srt1","(kg^2)/(kg^2)",zt)
          k = k + 1
 
        case ('srt2')
          isrt2 = k

          call stat_assign(isrt2,"srt2",
     .         "pdf parameter srt2","(kg^2)/(kg^2)",zt)
          k = k + 1
 
        case ('rc1')
          irc1 = k

          call stat_assign(irc1,"rc1",
     .         "pdf parameter rc1","kg/kg",zt)
          k = k + 1
 
        case ('rc2')
          irc2 = k

          call stat_assign(irc2,"rc2",
     .         "pdf parameter rc2","kg/kg",zt)
          k = k + 1
 
        case ('rsl1')
          irsl1 = k

          call stat_assign(irsl1,"rsl1",
     .         "pdf parameter rsl1","kg/kg",zt)
          k = k + 1
 
        case ('rsl2')
          irsl2 = k

          call stat_assign(irsl2,"rsl2",
     .         "pdf parameter rsl2","kg/kg",zt)
          k = k + 1
 
        case ('R1')
          iR1 = k
          call stat_assign(iR1,"R1",
     .         "pdf parameter R1","count",zt)
          k = k + 1
 
        case ('R2')
          iR2 = k

          call stat_assign(iR2,"R2",
     .         "pdf parameter R2","count",zt)
          k = k + 1
 
        case ('s1')
          is1 = k

          call stat_assign(is1,"s1",
     .         "pdf parameter s1","kg/kg",zt)
          k = k + 1
 
        case ('s2')
          is2 = k

          call stat_assign(is2,"s2",
     .         "pdf parameter s2","kg/kg",zt)
          k = k + 1
 
        case ('ss1')
          iss1 = k

          call stat_assign(iss1,"ss1",
     .         "pdf parameter ss1","kg/kg",zt)
          k = k + 1
 
        case ('ss2')
          iss2 = k

          call stat_assign(iss2,"ss2",
     .         "pdf parameter ss2","kg/kg",zt)
          k = k + 1
 
        case ('rrtthl')
          irrtthl = k

          call stat_assign(irrtthl,"rrtthl",
     .         "pdf parameter rrtthl","count",zt)
          k = k + 1

        case('wp2zt')
          iwp2zt = k

          call stat_assign(iwp2zt,"wp2zt",
     .         "wp2zt","m^2/s^2",zt)
          k = k + 1

        case('thlp2zt')
          ithlp2zt = k

          call stat_assign(ithlp2zt,"thlp2zt",
     .         "thlp2zt","K^2",zt)
          k = k + 1

        case('wpthlpzt')
          iwpthlpzt = k

          call stat_assign(iwpthlpzt,"wpthlpzt",
     .         "wpthlpzt","(m K)/s",zt)
          k = k + 1

        case('wprtpzt')   
          iwprtpzt = k

          call stat_assign(iwprtpzt,"wprtpzt",
     .         "wprtpzt","(m kg)/(s kg)",zt)
          k = k + 1

        case('rtp2zt')   
          irtp2zt = k

          call stat_assign(irtp2zt,"rtp2zt",
     .         "rtp2zt","kg/kg",zt)
          k = k + 1

        case('rtpthlpzt')   
          irtpthlpzt = k

          call stat_assign(irtpthlpzt,"rtpthlpzt",
     .         "rtpthlpzt","(kg K)/kg",zt)
          k = k + 1
#ifdef SCALARS 
        case ('sclram')
          isclram = k

          call stat_assign(isclram,"sclram",
     .         "passive scalar (currently thlm)","K",zt)
          k = k + 1

        case ('sclram_f')
          isclram_f = k

          call stat_assign(isclram_f,"sclram_f",
     .         "passive scalar (currently thlm_forcing)","K/s",zt)
          k = k + 1

        case ('sclrbm')
          isclrbm = k

          call stat_assign(isclrbm,"sclrbm",
     .         "passive scalar (currently rtm)","kg/kg",zt)
          k = k + 1

        case ('sclrbm_f')
          isclrbm_f = k
          call stat_assign(isclrbm_f,"sclrbm_f",
     .         "passive scalar (currently rtm_forcing)","kg/kg",zt)
          k = k + 1

        case ('edsclram')
          iedsclram = k

          call stat_assign(iedsclram,"edsclram",
     .         "eddy-diff scalar (currently thlm)","K",zt)
          k = k + 1

        case ('edsclrbm')
          iedsclrbm = k

          call stat_assign(iedsclrbm,"edsclrbm",
     .         "eddy-diff scalar (currently rt)","kg/kg",zt)
          k = k + 1
#endif /*SCALARS*/

        case default
          write(0,*) 'Error: unrecognized variable ',trim(vars_zt(i))
          lerror = .true.

        end select

      end do

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
     .                   time+stats_tout, stats_tout,
     .                   zm%nn )

      else ! Open NetCDF file
        call open_nc( zm%f, iunit, fdir, fname, 
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

      ivp2_bt  = 0
      ivp2_ta  = 0
      ivp2_tp  = 0
      ivp2_ma  = 0
      ivp2_dp1 = 0
      ivp2_dp2 = 0
      ivp2_pr1 = 0
      ivp2_pr2 = 0

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

      ! Flux budgets
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

      ! Variance budgets
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
     .         "up2 dissipation term","m^2/s^2",zm)
          k = k + 1

        case ('up2_dp2')
          iup2_dp2 = k
          call stat_assign(iup2_dp2,"up2_dp2",
     .         "up2 dissipation term","m^2/s^2",zm)
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
     .         "vp2 dissipating term","m^2/s^2",zm)
          k = k + 1

        case ('vp2_dp2')
          ivp2_dp2 = k

          call stat_assign(ivp2_dp2,"vp2_dp2",
     .         "vp2 dissipating term","m^2/s^2",zm)
          k = k + 1

        case ('vp2_pr1')
          ivp2_pr1 = k

          call stat_assign(ivp2_pr1,"vp2_pr1",
     .         "vp2 pressure term 1","m^2/s^2",zm)
          k = k + 1

        case ('vp2_pr2')
          ivp2_pr2 = k

          call stat_assign(ivp2_pr1,"vp2_pr1",
     .         "vp2 pressure term 1","m^2/s^2",zm)
          k = k + 1

#ifdef SCALARS
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

          call stat_assign(isclrathvp,"sclrathvp",
     .         "scalar(1)'th_v'","unknown",zm)
          k = k + 1
        case ('sclrbpthvp')
          isclrbpthvp = k
          call stat_assign(isclrbthvp,"sclrbthvp",
     .         "scalar(2)'th_v'","unknown",zm)
          k = k + 1
        case ('sclrapthlp')
          isclrapthlp = k

          call stat_assign(isclrathlp,"sclrathlp",
     .         "scalar(1)'th_l'","unknown",zm)
          k = k + 1
        case ('sclrbpthlp')
          isclrbpthlp = k

          call stat_assign(isclrbthlp,"sclrbthlp",
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

#endif /*SCALARS*/

        case default
          write(0,*) 'Error: unrecognized variable ', trim(vars_zm(i))
          lerror = .true.

        end select

      end do

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
     .                   time+stats_tout, stats_tout,
     .                   sfc%nn )

      else ! Open NetCDF files
        call open_nc( sfc%f, iunit, fdir, fname, 
     .                1, sfc%kk, sfc%z,
     .                day, month, year, rlat, rlon,
     .                time+stats_tout, stats_tout,
     .                sfc%nn )
      end if

      ! Default initialization for array indices for sfc

      iustar   = 0
      ilh      = 0
      ish      = 0
      icc      = 0
      ilwp     = 0
      izb      = 0
      izi      = 0
      irain    = 0   ! Brian
      ipflux   = 0   ! Brian
      irrm_sfc = 0   ! Brian

      ! These are estimates of the condition number on each implicit
      ! matrices, and not located at the surface of the domain.
      iwp23_cn    = 0
      irtm_cn     = 0
      ithlm_cn    = 0
      irtp2_cn    = 0
      ithlp2_cn   = 0
      irtpthlp_cn = 0
      iup2_cn     = 0
      ivp2_cn     = 0

      ! Assign pointers for statistics variables sfc

      k = 1
      do i=1,sfc%nn

        select case ( trim(vars_sfc(i)) )

        case ('ustar')
          iustar = k

          call stat_assign(iustar,"ustar",
     .         "ustar","m/s",sfc)
          k = k + 1

        case ('lh')
          ilh = k
          call stat_assign(ilh,"lh",
     .         "lh","W/m2",sfc)
          k = k + 1

        case ('sh')
          ish = k
          call stat_assign(ish,"sh",
     .         "sh","W/m2",sfc)
          k = k + 1

        case ('cc')
          icc = k
          call stat_assign(icc,"cc",
     .         "cc","count",sfc)
          k = k + 1

        case ('lwp')
          ilwp = k
          call stat_assign(ilwp,"lwp",
     .         "lwp","kg/m2",sfc)
          k = k + 1

        case ('zb')
          izb = k
          call stat_assign(izb,"zb",
     .         "zb","m",sfc)
          k = k + 1

        case ('zi')
          izi = k
          call stat_assign(izi,"zi",
     .         "zi","m",sfc)
          k = k + 1

        case ('rain')          ! Brian
          irain = k
          call stat_assign(irain,"rain_rate",
     .         "rainfall rate (mm/day)","mm/day",sfc)
          k = k + 1

        case ('pflux')         ! Brian
          ipflux = k

          call stat_assign(ipflux,"prec_flux",
     .         "Precipitation flux (W/m^2)","W/m^2",sfc)
          k = k + 1

        case ('rrm_sfc')       ! Brian
          irrm_sfc = k

          call stat_assign(irrm_sfc,"rrm_sfc",
     .         "Surface rain water mixing ratio","kg/kg",sfc)
          k = k + 1

        case ('wp23_cn')
          iwp23_cn = k
          call stat_assign(iwp23_cn,"wp23_cn",
     .         "Estimate of the condition number for wp2/3","count",sfc)
          k = k + 1

        case ('thlm_cn')
          ithlm_cn = k
	  call stat_assign(ithlm_cn,"thlm_cn",
     .         "Estimate of the condition number for thlm/wpthlp",
     .         "count",sfc)
          k = k + 1

        case ('rtm_cn')
          irtm_cn = k

          call stat_assign(irtm_cn,"rtm_cn",
     .         "Estimate of the condition number for rtm/wprtp",
     .         "count",sfc)
          k = k + 1

        case ('thlp2_cn')
          ithlp2_cn = k

          call stat_assign(ithlp2_cn,"thlp2_cn",
     .         "Estimate of the condition number for thlp2",
     .         "count",sfc)
          k = k + 1

        case ('rtp2_cn')
          irtp2_cn = k
  	  call stat_assign(irtp2_cn,"rtp2_cn",
     .         "Estimate of the condition number for rtp2",
     .         "count",sfc)
          k = k + 1

        case ('rtpthlp_cn')
          irtpthlp_cn = k
  	  call stat_assign(irtpthlp_cn,"rtpthlp_cn",
     .         "Estimate of the condition number for rtpthlp",
     .         "count",sfc)
          k = k + 1

        case ('up2_cn')
          iup2_cn = k
          call stat_assign(iup2_cn,"up2_cn",
     .         "Estimate of the condition number for up2","count",sfc)
          k = k + 1

        case ('vp2_cn')
          ivp2_cn = k  
	  call stat_assign(ivp2_cn,"vp2_cn",
     .         "Estimate of the condition number for vp2","count",sfc)

          k = k + 1

        case default
          write(0,*) 'Error: unrecognized variable ',trim(vars_sfc(i))
          lerror = .true.

        end select

      end do

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

      implicit none

      ! Input
      integer, intent(in) :: kk, nn

      ! Output
      real, dimension(kk,nn), intent(out)    :: x
      integer, dimension(kk,nn), intent(out) :: n

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

      implicit none

      ! Input
      integer, intent(in) :: nn, kk

      ! Output
      real, dimension(kk,nn), intent(inout)  :: x
      integer, dimension(kk,nn), intent(in)  :: n

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

      use hoc_stats

      implicit none

      ! Input

      real, intent(in) :: 
     .  time_elapsed,! Elapsed model time       [s]
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
      subroutine stats_end_timestep( time, delt )

!     Description:
!-----------------------------------------------------------------------

      use hoc_stats
      use outputgrads
      use outputnc

      implicit none

      ! Input

      real, intent(in) :: time, delt

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
         end if
        end do
      end do
      
      do i=1,zm%nn
       do k=1,zm%kk
         if ( zm%n(k,i) /= 0 
     .        .and. zm%n(k,i) /= int(stats_tout/stats_tsamp) ) then
           write(0,*) 'Possible sampling error for zm ',i,k,zm%n(k,i)
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
     .             wmt, Scm, taum, rcm, Ncm, Ncnm, Nim, cf, hydromet,
     .             sclrm, edsclrm, sclrm_forcing )

!     Description:
!     Accumulate those stats variables that are preserved in HOC from
!     timestep to timestep, but not those stats that are not, (e.g. 
!     budget terms, longwave and shortwave components, etc. )
!----------------------------------------------------------------------

      use hoc_stats
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
     .  Ncm,  ! Cloud droplet number concentration      [num/m^3]
     .  Ncnm, ! Cloud nuclei number concentration       [num/m^3]
     .  Nim,  ! Ice nuclei number concentration         [num/m^3]
     .  cf    ! Cloud fraction                          [%]

      real, intent(in), optional, dimension(gr%nnzp,sclr_dim) ::
     .  sclrm,          ! High-Order Passive scalar     [units vary]
     .  edsclrm,        ! Eddy-diff Passive scalar      [units vary] 
     .  sclrm_forcing   ! Large-scale forcing of scalar [units/s]

      ! Prognostic drizzle variable array
      real, intent(in), dimension(gr%nnzp,hydromet_dim) :: hydromet
      ! Contains:
      ! 1 rrm      Rain water mixing ratio               [kg/kg]
      ! 2 Nrm      Rain droplet number concentration     [num/m^3]
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
        if ( ilup > 0 ) then
          zt%x(:,ilup) = zt%x(:,ilup) + lup
          zt%n(:,ilup) = zt%n(:,ilup) + 1
        end if
        if ( ildown > 0 ) then
          zt%x(:,ildown) = zt%x(:,ildown) + ldown
          zt%n(:,ildown) = zt%n(:,ildown) + 1
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
        if ( iNcm > 0 ) then           ! Brian
          zt%x(:,iNcm) = zt%x(:,iNcm) + Ncm
          zt%n(:,iNcm) = zt%n(:,iNcm) + 1
        end if
        if ( iNcnm > 0 ) then
          zt%x(:,iNcnm) = zt%x(:,iNcnm) + Ncnm
          zt%n(:,iNcnm) = zt%n(:,iNcnm) + 1
        end if
        if ( iNim > 0 ) then
          zt%x(:,iNim) = zt%x(:,iNim) + Nim
          zt%n(:,iNim) = zt%n(:,iNim) + 1
        end if
        if ( ised_rcm > 0 .and. cloud_sed ) then       ! Brian
          zt%x(:,ised_rcm) = zt%x(:,ised_rcm) + sed_rcm
          zt%n(:,ised_rcm) = zt%n(:,ised_rcm) + 1
        end if
        if ( irsm > 0 ) then           ! Brian
          zt%x(:,irsm) = zt%x(:,irsm) + rsm
          zt%n(:,irsm) = zt%n(:,irsm) + 1
        end if

        if ( irrm > 0 ) then     ! Brian
          zt%x(:,irrm) = zt%x(:,irrm) + hydromet(:,1)
          zt%n(:,irrm) = zt%n(:,irrm) + 1
        end if

        if ( iNrm > 0 ) then     ! Brian
          zt%x(:,iNrm) = zt%x(:,iNrm) + hydromet(:,2)
          zt%n(:,iNrm) = zt%n(:,iNrm) + 1
        end if

        if ( irsnowm > 0 ) then
          zt%x(:,irsnowm) = zt%x(:,irsnowm) + hydromet(:,3)
          zt%n(:,irsnowm) = zt%n(:,irsnowm) + 1
        end if

        if ( iricem > 0 ) then
          zt%x(:,iricem) = zt%x(:,iricem) + hydromet(:,4)
          zt%n(:,iricem) = zt%n(:,iricem) + 1
        end if

        if ( irgraupelm > 0 ) then
          zt%x(:,irgraupelm) = zt%x(:,irgraupelm) + hydromet(:,5)
          zt%n(:,irgraupelm) = zt%n(:,irgraupelm) + 1
        end if

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

        if ( iFcsed > 0 .and. cloud_sed ) then    ! Brian
          zm%x(:,iFcsed) = zm%x(:,iFcsed) + Fcsed
          zm%n(:,iFcsed) = zm%n(:,iFcsed) + 1
        end if

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

      use hoc_stats
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

#endif /*STATS*/
