! CVS:  $Id: bugs_rad.F,v 1.2 2005-10-27 20:11:37 dschanen Exp $
! CVS:  $Name: not supported by cvs2svn $

!------------------------------------------------------------------------

#ifndef radoffline
      subroutine bugs_rad
     +               (   nlen ,      len ,     ibeg ,     nlm
     +,                    mm ,   qpintp ,  qpterp1 ,      pl2
     +,                  plk2 ,    phil2 ,       pl ,      plk
     +,                   dpl ,       tl ,       ql ,     qcwl
     +,                  qcil ,     qril ,    dlcup ,    dicup
     +,                   o3l ,    o3top ,    cldpc ,       ts
     +,                  tbmn ,     amu0 ,      slr ,    alvdf
     +,                 alndf ,    alvdr ,    alndr ,       s0
     +,                  ptop ,   tbgmin ,   tbgmax ,     grav
     +,                    cp ,      day ,    pcos1 ,    psin1
     +,                 pcos2 ,    psin2 ,      asl ,      atl
     +,                 LWtoa ,  LWdnsfc , LWnetsfc , SWnettoa
     +,               SWdnsfc , SWnetsfc ,    SWinc ,   radvbc
     +,                radvdc ,   radnbc ,   radndc ,   cldt2d
     +               )
#else
      subroutine bugs_rad
     +               (   nlen ,      len ,      nlm ,      pl2
     +,                    pl ,      dpl ,       tl ,       ql
     +,                  qcwl ,     qcil ,     qril ,      o3l
     +,                    ts ,     amu0 ,      slr ,    alvdf
     +,                 alndf ,    alvdr ,    alndr ,       s0
     +,                  grav ,       cp ,      asl ,      atl
     +,                  fdsw ,     fusw ,     fdlw ,     fulw
     +,                  acld )
#endif
      
      use kinds
      use variables_radiation_module, only: fdswcl, fuswcl 
     +,                                     fdlwcl, fulwcl

#ifdef USE_BUGSrad_ocast_random
      use bugs_cloud_correlate
#endif /* USE_BUGSrad_ocast_random */
      
      implicit none

!------------------------------------------------------------------------
! bugs_rad is the driver for the parameterization of the cloud fraction, 
! cloud optical properties, and long and short wave radiative fluxes.

! REFERENCES:
! Laura D. Fowler, and Phil Partain /slikrock (01-11-98).

! MODIFICATIONS:
! * added radoffline compiler directive; Phil Partain (04-04-00).

! * changed declarations to adapt the code from BUGS4 to BUGS5.
!   Laura D. Fowler/slikrock (02-01-00).

! SUBROUTINES CALLED:
!     EAUcf     :Calculation of the cloud fraction.
!     EAUcfQPs  :Calculation of diagnostics related to EAUcf. 
!     BUGSlwr   :Calculation of LW radiative heating rates and fluxes.
!     BUGSswr   :Calculation of SW radiative heating rates and fluxes.
!     BUGSradQPs:Calculation of diagnostics related to BUGSlwr and
!                BUGSswr.
! note: EAUcf, EAUcfQPs, and BUGSradQPs are not called in the offline 
!       version of this code.

! FUNCTIONS CALLED:
!     none.

! INCLUDED COMMON BLOCKS:
!     none.

! ARGUMENT LIST VARIABLES:
!     INPUT ARGUMENTS:
!     ----------------
      integer (kind=int_kind), intent(in)::
     &  nlen       !Length of total domain.
     &, len        !Length of sub domain.
     &, nlm        !Number of layers.

      real (kind=dbl_kind), intent(in)::
     &  grav       !Gravitational constant
     &, cp         !Constant pressure spec. heat of dry air  (J/(K-kg)).
     &, s0         !Solar constant                             (W/m^-2).

      real (kind=dbl_kind), intent(in), dimension(nlen)::
     &  ts         !Surface temperature                             (K).

      real (kind=dbl_kind), intent(in)::
     &  amu0       !Cosine of solar zenith angle                    (-).
     &, slr        !Fraction of daylight                            (-).
     &, alvdr      !Visible direct surface albedo                   (-).
     &, alndr      !Near-IR direct surface albedo                   (-).
     &, alvdf      !Visible diffuse surface albedo                  (-).
     &, alndf      !Near-IR diffuse surface albedo                  (-).

      real (kind=dbl_kind), intent(in), dimension(nlen,nlm)::
     &  pl         !Layer pressure                                (hPa).
     &, dpl        !Layer thickness                               (hPa).
     &, tl         !Temperature                                     (K).
     &, ql         !Specific humidity                           (kg/kg).
     &, qcwl       !Cloud water mixing ratio                    (kg/kg).
     &, qcil       !Cloud ice mixing ratio                      (kg/kg).
     &, qril       !Snow mixing ratio                           (kg/kg).
     &, o3l        !Ozone mixing ratio                          (kg/kg).
      real (kind=dbl_kind), intent(in), dimension(nlen,nlm+1)::
     &  pl2        !Level pressure                                (hPa).

#ifndef radoffline
      logical (kind=log_kind), intent(in)::
     &  qpintp     ! If true, starts computing QP diagnostics.

      integer (kind=int_kind), intent(in)::
     &  mm         !Number of interpolated layers.
     &, ibeg       !Beginning of the array. 

      real (kind=dbl_kind), intent(in)::
     &  day        !
     &, tbgmin     !Temperature below which only cloud ice forms    (K).
     &, tbgmax     !Temperature above which only cloud water forms  (K).
     &, pcos1      !Cosine of first diurnal harmonic                (-).
     &, psin1      !Sine of second diurnal harmonic                 (-).
     &, pcos2      !Cosine of first diurnal harmonic                (-).
     &, psin2      !Sine of second diurnal harmonic                 (-).
     &, ptop       !Top-of-model pressure level                   (hPa).

      real (kind=dbl_kind), intent(in), dimension(nlen)::
     &  tbmn       !PBL top temperature                             (K).
     &, o3top      !Ozone mixing ratio above the model top      (kg/kg).
     &, cldpc      !Fraction of PBL clouds                          (-).

      real (kind=dbl_kind), intent(in), dimension(nlen,nlm)::
     &  phil2      !Level geopotential height                       (m).
     &, plk        !(pl/p0)**kapa                                   (-).
     &, dlcup      !Cumulus detrainment of cloud water        (kg/kg/s).
     &, dicup      !Cumulus detrainment of cloud ice          (kg/kg/s).

      real (kind=dbl_kind), intent(in), dimension(nlen,nlm+1)::
     &  plk2       !(pl2/p0)**kapa                                  (-).

      real (kind=dbl_kind), intent(in), dimension(nlen,mm)::
     &  qpterp1
#else
      real (kind=dbl_kind), intent(in), dimension(nlen,nlm)::
     &  acld       !Radiative cloud fraction                        (-).
#endif


!     OUTPUT ARGUMENTS:
!     ------------------
      real (kind=dbl_kind), dimension(nlen)::
     &  radvbc     !SFC visible direct net SW radiation        (W/m^-2).
     &, radvdc     !SFC visible diffuse net SW radiation       (W/m^-2).
     &, radnbc     !SFC near-ir direct net SW radiation        (W/m^-2).
     &, radndc     !SFC near-ir diffuse net SW radiation       (W/m^-2).
     &, radvbcc    !SFC visible direct clear net SW radiation  (W/m^-2).
     &, radvdcc    !SFC visible diffuse clear net SW radiation (W/m^-2).
     &, radnbcc    !SFC near-ir direct clear net SW radiation  (W/m^-2).
     &, radndcc    !SFC near-ir diffuse clear net SW radiation (W/m^-2).

#ifndef radoffline
      real (kind=dbl_kind), dimension(len,nlm+2):: 
#else
      real (kind=dbl_kind), dimension(len,nlm+1):: 
#endif
     &  fdsw       !Downward SW flux                           (W/m^-2).
     &, fusw       !Upward SW flux                             (W/m^-2).
     &, fdlw       !Downward LW flux                           (W/m^-2).
     &, fulw       !Upward LW flux                             (W/m^-2).
      real (kind=dbl_kind), dimension(nlen,nlm)::
     &  atl        !All-sky LW radiative heating rate             (K/s).
     &, asl        !All-sky SW radiative heating rate             (K/s).      
     &, atlcl      !Clear-sky LW radiative heating rate           (K/s).
     &, aslcl      !Clear-sky SW radiative heating rate           (K/s).
#ifndef radoffline
      real (kind=dbl_kind), intent(out), dimension(nlen)::
     &  cldt2d     !Vertically-integrated cloud fraction            (-).
     &, SWinc      !TOA incident solar radiation               (W/m^-2).
     &, SWnettoa   !TOA net SW radiation                       (W/m^-2).
     &, SWnetsfc   !SFC net SW radiation                       (W/m^-2).
     &, SWdnsfc    !SFC downward SW radiation                  (W/m^-2).
     &, LWtoa      !TOA outgoing LW radiation                  (W/m^-2).
     &, LWdnsfc    !SFC downward LW radiation                  (W/m^-2).
     &, LWnetsfc   !SFC net LW radiation                       (W/m^-2).
     &  SWtoa      !TOA all-sky net SW radiation               (W/m^-2).
     &, SWnettoacl !TOA clear-sky net SW radiation             (W/m^-2).
     &, SWdnsfccl  !SFC clear-sky net SW radiation             (W/m^-2).
     &, SWnetsfccl !SFC clear-sky net SW radiation             (W/m^-2).
     &, SWtom      !TOM all-sky net SW radiation               (W/m^-2).
     &, SWtomcl    !TOM clear-sky net SW radiation             (W/m^-2).
     &, LWtoacl    !TOA clear-sky outgoing LW radiation        (W/m^-2).
     &, LWupsfc    !SFC upward LW radiation                    (W/m^-2).
     &, LWdnsfccl  !SFC clear-sky downward LW radiation        (W/m^-2).
     &, LWnetsfccl !SFC clear-sky net LW radiation             (W/m^-2).
     &, LWtom      !TOM all-sky net LW radiation               (W/m^-2).
     &, LWtomcl    !TOM clear-sky net LW radiation             (W/m^-2).
#endif

! LOCAL VARIABLES:

      character(len=80), save::
     &  cvs_version1
     &, cvs_version2

      logical (kind=log_kind)::
     &  sel_rules_lw
     &, sel_rules_sw

      logical (kind=log_kind), dimension(len)::
     &  bitx      !Scans daytime grid-points                       (-).
     
      integer (kind=int_kind)::
     &  i, l, ll
     &, nnp        !Number of layers plus one (TOA "sponge" layer).
     &, nday       !Number of daytime grid-points .

      integer (kind=int_kind), dimension(len)::
     &  iday       !Location of daytime grid-points.
     
      real (kind=dbl_kind)::
     &  heat_fac   !Mutiplying factor         .
     &, umco2      !Concentration of CO2 (ppm).
     &, umch4      !Concentration of CH4      .
     &, umn2o      !Concentration of N2O      .

! Kurt Kotenberg says on 05 APR 05 that to turn off gases, the data
! below must be set to 0.0

      data umco2,umch4,umn2o /330.,1.6,0.28/
!      data umco2,umch4,umn2o /0.0,0.0,0.0/

! End comments

      real (kind=dbl_kind), dimension(len)::
     &  ts_loc     !Local surface temperature                       (K).
     &, amu0_loc   !Local cosine of solar zenith angle              (-).
     &, slr_loc    !Local fraction of daylight                      (-).
     &, alvdr_loc  !Local visible direct surface albedo             (-).
     &, alndr_loc  !Local near-IR direct surface albedo             (-).
     &, alvdf_loc  !Local visible diffuse surface albedo            (-).
     &, alndf_loc  !Local near-IR diffuse surface albedo            (-).
     &, acldmx     !Local Maximum cloud fraction                    (-).
     &, den        !Multiplying factor                                 .
     &, delf       !Net flux
     &, delfcl     !Net flux (clear)
     &, radvbc_loc !Local SFC visible direct net SW radiation  (W/m^-2).
     &, radvdc_loc !Local SFC visible direct net SW radiation  (W/m^-2).
     &, radnbc_loc !Local SFC near-ir direct net SW radiation  (W/m^-2).
     &, radndc_loc !Local SFC near-ir direct net SW radiation  (W/m^-2).
     &, radvbcc_loc!As radvbc_loc, but clear-sky               (W/m^-2).
     &, radvdcc_loc!As radvbc_loc, but clear-sky               (W/m^-2).
     &, radnbcc_loc!As radncb_loc, but clear-sky               (W/m^-2).
     &, radndcc_loc!As radncb_loc, but clear-sky               (W/m^-2).

#ifndef radoffline
      real (kind=dbl_kind), dimension(len,nlm+1)::
#else
      real (kind=dbl_kind), dimension(len,nlm)::
#endif
     &  pl_loc     !Local layer pressure                          (hPa).
     &, dpl_loc    !Local layer thickness                         (hPa).
     &, tl_loc     !Local temperature                               (K).
     &, ql_loc     !Local specific humidity                     (kg/kg).
     &, qcwl_loc   !Local cloud water mixing ratio              (kg/kg).
     &, qcil_loc   !Local cloud ice mixing ratio                (kg/kg).
     &, o3l_loc    !Local ozone mixing ratio                    (kg/kg).
     &, acld_loc   !Local radiative cloud fraction                  (-).
     &, rmix       !Water vapor mixing ratio                    (kg/kg).
     &, cwrho      !Density of cloud water                     (g/m^-3).
     &, cirho      !Density of cloud ice                       (g/m^-3).
     &, o3mix      !Ozone mixing ratio                          (kg/kg).
     &, b1         !Cloud overlap parameter                         (-).
     &, b2         !Cloud overlap parameter                         (-).
     &, b3         !Cloud overlap parameter                         (-).
     &, b4         !Cloud overlap parameter                         (-).

#ifndef radoffline
      real (kind=dbl_kind), dimension(len,nlm+2)::
#else
      real (kind=dbl_kind), dimension(len,nlm+1)::
#endif
     &  pl2_loc    !Local level pressure                          (hPa).

      real (kind=dbl_kind), dimension(:), allocatable::
     &  ts_day     !As ts_loc, but for daytime grid-points          (K).
     &, amu0_day   !As amu0_loc,but for daytime grid-points         (-).
     &, slr_day    !As slr_loc,but for daytime grid-points          (-).
     &, alvdf_day  !As alvdf_loc,but for daytime grid-points        (-).
     &, alndf_day  !As alndf_loc,but for daytime grid-points        (-).
     &, alvdr_day  !As alvdr_loc,but for daytime grid-points        (-).
     &, alndr_day  !As alndr_loc,but for daytime grid-points        (-).
     &, acldmx_day !As acldmx_loc,but for daytime grid-points       (-).

      real (kind=dbl_kind), dimension(:,:), allocatable::
     &  pl_day     !As pl,but for daytime grid-points             (hPa).
     &, dpl_day    !As dpl,but for daytime grid-points            (hPa).
     &, tl_day     !As tl,but for daytime grid-points               (K).
     &, rmix_day   !As rmix,but for daytime grid-points         (kg/kg).
     &, cwrho_day  !As cwrho,but for daytime grid-points       (g/m^-3).
     &, cirho_day  !As cirho,but for daytime grid-points       (g/m^-3).
     &, o3mix_day  !As o3mix,but for daytime grid-points        (kg/kg).
     &, acld_day   !As acld,but for daytime grid-points             (-).

      real (kind=dbl_kind), dimension(:,:), allocatable::
     &  pl2_day    !As pl2_loc, but for daytime grid-points       (hPa).
     &, fdsw_day   !As fdsw, but for daytime grid points       (W/m^-2).
     &, fdswcl_day !As fdswcl, but for daytime grid points     (W/m^-2).
     &, fusw_day   !As fusw, but for daytime grid points       (W/m^-2).
     &, fuswcl_day !As fuswcl, but for daytime grid points     (W/m^-2).
     &, b1_day     !As b1, but for daytime grid points              (-).
     &, b2_day     !As b2, but for daytime grid points              (-).
     &, b3_day     !As b3, but for daytime grid points              (-).
     &, b4_day     !As b4, but for daytime grid points              (-).

      real (kind=dbl_kind), dimension(:), allocatable::      
     &  radvbc_day !As radvbc_loc, but for daytime grid-points (W/m^-2).
     &, radvdc_day !As radvbc_loc, but for daytime grid-points (W/m^-2).
     &, radnbc_day !As radnbc_loc, but for daytime grid-points (W/m^-2).
     &, radndc_day !As radnbc_loc, but for daytime grid-points (W/m^-2).
     &, radvbcc_day!As radvbcc, but for daytime grid-points    (W/m^-2).
     &, radvdcc_day!As radvbcc, but for daytime grid-points    (W/m^-2).
     &, radnbcc_day!As radnbc_loc, but for daytime grid-points (W/m^-2).
     &, radndcc_day!As radnbc_loc, but for daytime grid-points (W/m^-2).

#ifndef radoffline
      real (kind=dbl_kind), dimension(len)::
     &  o3top_loc  !Local ozone mixing ratio above TOM          (kg/kg).
      real (kind=dbl_kind), dimension(nlen,nlm)::
     &  clu        !Anvil cloud fraction                            (-).
     &, cld        !Large-scale saturation cloud fraction           (-).
     &, cldw       !Fraction of water clouds                        (-).
     &, cldi       !Fraction of ice clouds                          (-).
     &, cldm       !Fraction of mixed-phase clouds                  (-).
     &, cldt       !Total cloud fraction                            (-).      
#endif

#ifdef USE_BUGSrad_ocast_random
      real (kind=dbl_kind), dimension(len) ::
     & corr_length    !Corr. length for total cloud calculation     (-).
     &,ctot_corr      !Correlated total cloud amount                (-).
     &,c_maximal      !                                             (-).
#ifndef radoffline
      real (kind=dbl_kind), dimension(len,nlm+1) ::
#else
      real (kind=dbl_kind), dimension(len,nlm)   ::
#endif
     & cf_max         !                                             (-).
     &,cf_rndm        !                                             (-).

       real (kind=dbl_kind), dimension(:), allocatable ::
     & c_maximal_day  !

       real (kind=dbl_kind), dimension(:,:), allocatable ::
     & cf_max_day     !
     &,cf_rndm_day    
#endif

!-----------------------------------------------------------------------
!     print*,'---- enter subroutine bugs_rad:'

!---- 1.0 initialization of output variables and arrays:

      acld_loc(:,:)      = 0.
      fdsw(:,:)          = 0.
      fdswcl(:,:)        = 0.
      fusw(:,:)          = 0.
      fuswcl(:,:)        = 0.
      fdlw(:,:)          = 0.
      fdlwcl(:,:)        = 0.
      fulw(:,:)          = 0.
      fulwcl(:,:)        = 0.
      atlcl(1:len,:)     = 0.
      aslcl(1:len,:)     = 0.
      atl(1:len,:)       = 0.
      asl(1:len,:)       = 0.
      radvbc(1:len)      = 0.
      radvdc(1:len)      = 0.
      radnbc(1:len)      = 0.
      radndc(1:len)      = 0.
      radvbcc(1:len)     = 0.
      radvdcc(1:len)     = 0.
      radnbcc(1:len)     = 0.
      radndcc(1:len)     = 0.
      radvbc_loc(1:len)  = 0.
      radvdc_loc(1:len)  = 0.
      radnbc_loc(1:len)  = 0.
      radndc_loc(1:len)  = 0.
      radvbcc_loc(1:len) = 0.
      radvdcc_loc(1:len) = 0.
      radnbcc_loc(1:len) = 0.
      radndcc_loc(1:len) = 0.

#ifndef radoffline
      SWtoa(1:len)       = 0.
      SWnettoacl(1:len)  = 0.
      SWdnsfccl(1:len)   = 0.
      SWnetsfccl(1:len)  = 0.
      SWtom(1:len)       = 0.
      SWtomcl(1:len)     = 0.
      LWtoacl(1:len)     = 0.
      LWupsfc(1:len)     = 0.
      LWdnsfccl(1:len)   = 0.
      LWnetsfccl(1:len)  = 0.
      LWtom(1:len)       = 0.
      LWtomcl(1:len)     = 0.
      cldt2d(1:len)      = 0.
      SWinc(1:len)       = 0.
      SWnettoa(1:len)    = 0.
      SWdnsfc(1:len)     = 0.
      SWnetsfc(1:len)    = 0.
      LWtoa(1:len)       = 0.
      LWdnsfc(1:len)     = 0.
      LWnetsfc(1:len)    = 0.
      do i = 1, len
         if(amu0_loc(i) .ge. .01) SWinc(i)= s0*slr_loc(i)*amu0_loc(i)
      enddo   
#endif
                 
!---- 1.1 initialization of local scalars

      cvs_version1 =
     & "Version: 
     & $Id: bugs_rad.F,v 1.2 2005-10-27 20:11:37 dschanen Exp $"
      cvs_version2 =
     & "Version:$Name: not supported by cvs2svn $"

      sel_rules_sw = .false.
      sel_rules_lw = .false.

      heat_fac = grav*0.01/cp

!---- 1.2 initialization of local arrays

      ts_loc(:)          = ts(1:len)
      amu0_loc(:)        = amu0
      slr_loc(:)         = slr
      alvdr_loc(:)       = alvdr
      alndr_loc(:)       = alndr
      alvdf_loc(:)       = alvdf
      alndf_loc(:)       = alndf

#ifndef radoffline
      nnp = nlm+1
      pl2_loc(:,1)       = 0.0001
      pl2_loc(:,2:nnp+1) = pl2(1:len,1:nlm+1)
      pl_loc(:,1)        = 
     +   pl2_loc(1:len,2)-(pl2_loc(1:len,2)-pl2_loc(1:len,1))/2.
      pl_loc(:,2:nnp)    = pl(1:len,1:nlm)
      dpl_loc(:,1)       = pl2(1:len,1)
      dpl_loc(:,2:nnp)   = dpl(1:len,1:nlm)
      tl_loc(:,1)        = tl(1:len,1)
      tl_loc(:,2:nnp)    = tl(1:len,1:nlm)
      ql_loc(:,1)        = 3.e-06
      ql_loc(:,2:nnp)    = ql(1:len,1:nlm)
      qcwl_loc(:,1)      = 0.
      qcwl_loc(:,2:nnp)  = qcwl(1:len,1:nlm)
      qcil_loc(:,1)      = 0.
      qcil_loc(:,2:nnp)  = qcil(1:len,1:nlm)
      o3l_loc(:,1)       = o3top(1:len)
      o3l_loc(:,2:nnp)   = o3l(1:len,1:nlm)

      !Initialize acld_loc:  cloud fraction and diagnostics obtained from GCM data
      call EAUcf(nlen,len,nlm,tbgmin,tbgmax,cldpc,tl,qcwl,qcil,qril
     +,          dlcup,dicup,clu,cld,cldt,cldw,cldi,cldm,cldt2d)
      
      if(qpintp)
     +   call EAUcfQPs(nlen,len,ibeg,pcos1,psin1,pcos2,psin2,ptop
     +,            qpterp1,pl2,plk2,plk,tl,cld,clu,cldt,cldw,cldi
     +,            cldm,cldpc)

      acld_loc(:,1)     = 0.
      acld_loc(:,2:nnp) = cld(:,1:nlm) + clu(:,1:nlm)
      do l = 1, nnp
         do i = 1, len
            if(acld_loc(i,l) .ge. 0.5) acld_loc(i,l) = 1.
            if(acld_loc(i,l) .lt. 0.5) acld_loc(i,l) = 0.
         enddo
      enddo

#else
      nnp = nlm
      pl2_loc(:,:) = pl2(:,:)
      pl_loc(:,:)  = pl(:,:)
      dpl_loc(:,:) = dpl(:,:)
      tl_loc(:,:)  = tl(:,:)
      ql_loc(:,:)  = ql(:,:)
      qcwl_loc(:,:)= qcwl(:,:)
      qcil_loc(:,:)= qcil(:,:)
      o3l_loc(:,:) = o3l(:,:)
      acld_loc(:,:)= acld(1:len,1:nlm)
#endif

      acldmx(1:len) = 0.0 ! Un-initialized memory -dschanen 15 Dec 2009
      do l = 1, nlm
         acldmx(1:len) = max(acldmx(1:len),acld_loc(1:len,l))
      enddo



!--- 2. computation of cloud overlap parameters
      do i = 1, len
        b1(i,1) = 1.0 - acld_loc(i,1)
        b3(i,1) = 1.0
        do l = 2,nlm
          if (acld_loc(i,l-1).eq.1.0) then
            b1(i,l) = 1.0
          else
            b1(i,l) = (1.0 - max(acld_loc(i,l),acld_loc(i,l-1))) /
     *              (1.0 - acld_loc(i,l-1))
          endif
          if (acld_loc(i,l-1).eq.0.0) then
            b3(i,l) = 1.0
          else
            b3(i,l) = min(acld_loc(i,l),acld_loc(i,l-1)) /
     &                acld_loc(i,l-1)
          endif
        enddo
 
        b2(i,nlm) = 1.0 - acld_loc(i,nlm)
        b4(i,nlm) = 1.0
        do l = 1,nlm-1
          if (acld_loc(i,l+1).eq.1.0) then
            b2(i,l) = 1.0
          else
            b2(i,l) = (1.0 - max(acld_loc(i,l),acld_loc(i,l+1))) /
     *              (1.0 - acld_loc(i,l+1))
          endif
          if (acld_loc(i,l+1).eq.0.0) then
            b4(i,l) = 1.0
          else
            b4(i,l) = min(acld_loc(i,l),acld_loc(i,l+1)) /
     &                acld_loc(i,l+1)
          endif
        enddo
      enddo
!      do l = 1,nlm
!        print *,l,acld_loc(1,l),b1(1,l),b2(1,l),b3(1,l),b4(1,l)
!      enddo


      do l = 1, nnp
         den(1:len)     = pl_loc(1:len,l)*100/(287.*tl_loc(1:len,l))
         rmix(1:len,l)  = ql_loc(1:len,l)/(1.-ql_loc(1:len,l))
#ifdef nooverlap
         cwrho(1:len,l) = den*1000.*qcwl_loc(1:len,l)*acld_loc(1:len,l)
         cirho(1:len,l) = den*1000.*qcil_loc(1:len,l)*acld_loc(1:len,l)
#else
         cwrho(1:len,l) = den*1000.*qcwl_loc(1:len,l)
         cirho(1:len,l) = den*1000.*qcil_loc(1:len,l)
#endif
         o3mix(1:len,l) = o3l_loc(1:len,l)

! Kurt Kotenberg added the following print statement in efforts to create the needed info
! for bugradDEN.dat, which is used in BUGSrad_profile.m 

!       print*, l, den, acld(1,l)
!       print *, l, tl(1,l), ql(1,l)
! End Kurt Kotenberg's additions


      enddo                  

#ifdef USE_BUGSrad_ocast_random
!Initialize the cloud amounts for overcast random
       corr_length = 3000. !m
       call bugs_ctot(nlen,len,nnp,pl2,tl,acld,corr_length,ctot_corr)
       call bugs_cloudfit(len, nnp, ctot_corr, acld, c_maximal,
     &                    cf_max, cf_rndm)
#endif


!--   3. call to the bugs_rad longwave radiation code:
#ifdef USE_BUGSrad_ocast_random
      call bugs_lwr
     +          (        len ,     nnp , pl2_loc , pl_loc
     +,              dpl_loc ,  tl_loc ,    rmix ,  cwrho
     +,                cirho ,   o3mix ,  ts_loc ,acld_loc
     *,               acldmx ,      b1 ,      b2 ,     b3
     +,                   b4 ,   umco2 ,   umch4 ,  umn2o
     +,                 fdlw ,  fdlwcl ,    fulw , fulwcl
     +,         sel_rules_lw ,         c_maximal   
     +,               cf_max , cf_rndm )
#else
!     print*,'---- enter subroutine bugs_lwr:'
      call bugs_lwr
     +          (         len ,    nnp , pl2_loc , pl_loc
     +,               dpl_loc , tl_loc ,    rmix ,  cwrho
     +,                 cirho ,  o3mix ,  ts_loc ,acld_loc
     +,                acldmx ,     b1 ,      b2 ,     b3
     +,                    b4 ,  umco2 ,   umch4 ,  umn2o
     +,                  fdlw , fdlwcl ,    fulw , fulwcl
     +,          sel_rules_lw
     +           )
#endif /* USE_BUGSrad_ocast_random */

!     print*,'---- exit subroutine bugs_lwr:'
         
!--   4. call to the bugs_rad shortwave radiation code:

!--   note: slr needs to be modified to accomodate the difference in
!           s0 between the original radiation code and that in BUGS2.
      slr_loc(1:len) = slr_loc(1:len) * s0/1339.945

!--   note: computation of the shortwave radiative heating rates and
!     fluxes are made for daytime grid-points only: 
      bitx(1:len) = amu0_loc(1:len) .ge. 0.01
      nday = 0
      do i = 1, len
         if(bitx(i)) then 
            nday       = nday + 1
            iday(nday) = i 
         endif
      enddo
!      print *, "nday", nday
!      print *, "amu0_loc", amu0_loc
!      print *, "amu0", amu0
      if (nday .eq. 0) goto 1000

      allocate(ts_day(nday))
      allocate(amu0_day(nday))
      allocate(slr_day(nday))
      allocate(alvdr_day(nday))      
      allocate(alndr_day(nday))      
      allocate(alvdf_day(nday))
      allocate(alndf_day(nday))      
      allocate(acldmx_day(nday))
      allocate(pl_day(nday,nnp))
      allocate(dpl_day(nday,nnp))
      allocate(tl_day(nday,nnp))
      allocate(rmix_day(nday,nnp))
      allocate(cwrho_day(nday,nnp))
      allocate(cirho_day(nday,nnp))
      allocate(o3mix_day(nday,nnp))
      allocate(acld_day(nday,nnp))
      allocate(b1_day(nday,nnp))
      allocate(b2_day(nday,nnp))
      allocate(b3_day(nday,nnp))
      allocate(b4_day(nday,nnp))

      allocate(pl2_day(nday,nnp+1))

      allocate(fdsw_day(nday,nnp+1))
      allocate(fdswcl_day(nday,nnp+1))      
      allocate(fusw_day(nday,nnp+1))
      allocate(fuswcl_day(nday,nnp+1))           
      allocate(radvbc_day(nday))
      allocate(radvdc_day(nday))
      allocate(radnbc_day(nday))
      allocate(radndc_day(nday))
      allocate(radvbcc_day(nday))
      allocate(radvdcc_day(nday))
      allocate(radnbcc_day(nday))
      allocate(radndcc_day(nday))

      ts_day(1:nday)          = ts_loc(iday(1:nday))
      amu0_day(1:nday)        = amu0_loc(iday(1:nday))
      slr_day(1:nday)         = slr_loc(iday(1:nday))
      alvdf_day(1:nday)       = alvdf_loc(iday(1:nday))
      alndf_day(1:nday)       = alndf_loc(iday(1:nday))
      alvdr_day(1:nday)       = alvdr_loc(iday(1:nday))
      alndr_day(1:nday)       = alndr_loc(iday(1:nday))
      acldmx_day(1:nday)      = acldmx(iday(1:nday))
      pl_day(1:nday,1:nnp)    = pl_loc(iday(1:nday),1:nnp)
      dpl_day(1:nday,1:nnp)   = dpl_loc(iday(1:nday),1:nnp)
      tl_day(1:nday,1:nnp)    = tl_loc(iday(1:nday),1:nnp)
      rmix_day(1:nday,1:nnp)  = rmix(iday(1:nday),1:nnp)
      cwrho_day(1:nday,1:nnp) = cwrho(iday(1:nday),1:nnp)
      cirho_day(1:nday,1:nnp) = cirho(iday(1:nday),1:nnp)
      o3mix_day(1:nday,1:nnp) = o3mix(iday(1:nday),1:nnp)
      acld_day(1:nday,1:nnp)  = acld_loc(iday(1:nday),1:nnp)
      b1_day(1:nday,1:nnp)    = b1(iday(1:nday),1:nnp)
      b2_day(1:nday,1:nnp)    = b2(iday(1:nday),1:nnp)
      b3_day(1:nday,1:nnp)    = b3(iday(1:nday),1:nnp)
      b4_day(1:nday,1:nnp)    = b4(iday(1:nday),1:nnp)
      pl2_day(1:nday,1:nnp+1) = pl2_loc(iday(1:nday),1:nnp+1)

#ifdef USE_BUGSrad_ocast_random
      allocate(c_maximal_day(nday))
      allocate(cf_max_day(nday,nlm))
      allocate(cf_rndm_day(nday,nlm))

      c_maximal_day(1:nday) = c_maximal(iday(1:nday))
      cf_max_day(1:nday,1:nlm) = cf_max(iday(1:nday),1:nlm)
      cf_rndm_day(1:nday,1:nlm) = cf_rndm(iday(1:nday),1:nlm)
#endif

!     print*,'---- enter subroutine bugs_swr:'
#ifdef USE_BUGSrad_ocast_random

!      print *, "Calling bugs_swr"

      call bugs_swr
     +          (      nday ,         nnp ,      pl2_day ,      pl_day
     +,             dpl_day ,      tl_day ,     rmix_day ,   cwrho_day
     +,           cirho_day ,   o3mix_day ,       ts_day ,    amu0_day
     +,             slr_day ,   alvdf_day ,    alndf_day ,   alvdr_day
     +,           alndr_day ,    acld_day ,   acldmx_day ,       umco2
     +,               umch4 ,       umn2o ,       b1_day ,      b2_day
     +,              b3_day ,      b4_day ,     fdsw_day ,  fdswcl_day
     +,            fusw_day ,  fuswcl_day ,   radvbc_day , radvbcc_day
     +,          radvdc_day , radvdcc_day ,   radnbc_day , radnbcc_day
     +,          radndc_day , radndcc_day , sel_rules_sw
     +,       c_maximal_day ,  cf_max_day ,  cf_rndm_day
     +          )
#else
!      print *, "Calling bugs_swr"

      call bugs_swr
     +          (      nday ,         nnp ,      pl2_day ,      pl_day
     +,             dpl_day ,      tl_day ,     rmix_day ,   cwrho_day
     +,           cirho_day ,   o3mix_day ,       ts_day ,    amu0_day
     +,             slr_day ,   alvdf_day ,    alndf_day ,   alvdr_day
     +,           alndr_day ,    acld_day ,   acldmx_day ,       umco2
     +,               umch4 ,       umn2o ,       b1_day ,      b2_day
     +,              b3_day ,      b4_day ,     fdsw_day ,  fdswcl_day
     +,            fusw_day ,  fuswcl_day ,   radvbc_day , radvbcc_day
     +,          radvdc_day , radvdcc_day ,   radnbc_day , radnbcc_day
     +,          radndc_day , radndcc_day , sel_rules_sw
     +          )
#endif  /* USE_BUGSrad_ocast_random */

!     print*,'---- end subroutine bugs_swr:'

      radvbc_loc(iday(1:nday))  = radvbc_day(1:nday)
      radvdc_loc(iday(1:nday))  = radvdc_day(1:nday)
      radnbc_loc(iday(1:nday))  = radnbc_day(1:nday)
      radndc_loc(iday(1:nday))  = radndc_day(1:nday)
      radvbcc_loc(iday(1:nday)) = radvbcc_day(1:nday)
      radvdcc_loc(iday(1:nday)) = radvdcc_day(1:nday)
      radnbcc_loc(iday(1:nday)) = radnbcc_day(1:nday)
      radndcc_loc(iday(1:nday)) = radndcc_day(1:nday)
      fdsw(iday(1:nday),:)      = fdsw_day(1:nday,:)
      fdswcl(iday(1:nday),:)    = fdswcl_day(1:nday,:)
      fusw(iday(1:nday),:)      = fusw_day(1:nday,:)
      fuswcl(iday(1:nday),:)    = fuswcl_day(1:nday,:)

      deallocate(ts_day)
      deallocate(amu0_day)
      deallocate(slr_day)
      deallocate(alvdr_day)      
      deallocate(alndr_day)      
      deallocate(alvdf_day)
      deallocate(alndf_day)      
      deallocate(acldmx_day)
      deallocate(pl_day)
      deallocate(dpl_day)
      deallocate(tl_day)
      deallocate(rmix_day)
      deallocate(cwrho_day)
      deallocate(cirho_day)
      deallocate(o3mix_day)
      deallocate(acld_day)
      deallocate(pl2_day)
      deallocate(b1_day)
      deallocate(b2_day)
      deallocate(b3_day)
      deallocate(b4_day)
      deallocate(radvbc_day)
      deallocate(radvdc_day)
      deallocate(radnbc_day)
      deallocate(radndc_day)
      deallocate(radvbcc_day)
      deallocate(radvdcc_day)
      deallocate(radnbcc_day)
      deallocate(radndcc_day)
      deallocate(fdsw_day)
      deallocate(fdswcl_day)      
      deallocate(fusw_day)
      deallocate(fuswcl_day)

#ifdef USE_BUGSrad_ocast_random
      deallocate(c_maximal_day, cf_max_day, cf_rndm_day)
#endif /* USE_BUGSrad_ocast_random */

 1000 continue

!---- 5. computation of long and short wave radiative heating rates:
      do l = 1, nlm
#ifndef radoffline
         ll = l + 1
#else
         ll = l
#endif
         delf(1:len)    = fulw(1:len,ll)-fdlw(1:len,ll)
     +                  - fulw(1:len,ll+1)+fdlw(1:len,ll+1)
         delfcl(1:len)  = fulwcl(1:len,ll)-fdlwcl(1:len,ll)
     +                  - fulwcl(1:len,ll+1)+ fdlwcl(1:len,ll+1)
        
!         print *, l, dpl(1,l)

         atl(1:len,l)   = -heat_fac*delf/dpl(1:len,l)
         atlcl(1:len,l) = -heat_fac*delfcl/dpl(1:len,l)

         delf(1:len)    = fusw(1:len,ll)-fdsw(1:len,ll)
     +                  - fusw(1:len,ll+1)+fdsw(1:len,ll+1)
         delfcl(1:len)  = fuswcl(1:len,ll)-fdswcl(1:len,ll)
     +                  - fuswcl(1:len,ll+1) + fdswcl(1:len,ll+1)
         asl(1:len,l)   = -heat_fac*delf/dpl(1:len,l)
         aslcl(1:len,l) = -heat_fac*delfcl/dpl(1:len,l)             

      enddo 

#ifndef radoffline
!---- top-of-the atmosphere and surface lw and sw radiative fluxes:
      LWtoa(1:len)      = fulw(1:len,1)
      LWtoacl(1:len)    = fulwcl(1:len,1)        
      LWupsfc(1:len)    = fulw(1:len,nnp+1)
      LWdnsfc(1:len)    = fdlw(1:len,nnp+1)
      LWdnsfccl(1:len)  = fdlwcl(1:len,nnp+1)
      LWnetsfc(1:len)   = fulw(1:len,nnp+1)-fdlw(1:len,nnp+1)
      LWnetsfccl(1:len) = fulwcl(1:len,nnp+1)-fdlwcl(1:len,nnp+1)

      SWtoa(1:len)      = fdsw(1:len,1) 
      SWnettoa(1:len)   = fdsw(1:len,1)-fusw(1:len,1)
      SWnettoacl(1:len) = fdswcl(1:len,1)-fuswcl(1:len,1)
      SWdnsfc(1:len)    = fdsw(1:len,nnp+1) 
      SWdnsfccl(1:len)  = fdswcl(1:len,nnp+1)
      SWnetsfc(1:len)   = fdsw(1:len,nnp+1)-fusw(1:len,nnp+1) 
      SWnetsfccl(1:len) = fdswcl(1:len,nnp+1)-fuswcl(1:len,nnp+1)  

!---- top-of-model lw and sw radiative fluxes:
      LWtom(1:len)      = fulw(1:len,2)-fdlw(1:len,2)
      LWtomcl(1:len)    = fulwcl(1:len,2)-fdlwcl(1:len,2)
      SWtom(1:len)      = fdsw(1:len,2)-fusw(1:len,2)
      SWtomcl(1:len)    = fdswcl(1:len,2)-fuswcl(1:len,2)
#endif

!---- back to full arrays:
      radvbc(1:len)  = radvbc_loc(1:len)
      radvdc(1:len)  = radvdc_loc(1:len)
      radnbc(1:len)  = radnbc_loc(1:len)
      radndc(1:len)  = radndc_loc(1:len)
      radvbcc(1:len) = radvbcc_loc(1:len)
      radvdcc(1:len) = radvdcc_loc(1:len)
      radnbcc(1:len) = radnbcc_loc(1:len)
      radndcc(1:len) = radndcc_loc(1:len)
      
#ifndef radoffline
!---- begin calculations of diagnostics:
!     print*,'---- enter subroutine BUGSradQPs:'
      if(qpintp) 
     +   call BUGSradQPs(nlen,len,ibeg,day,pcos1,psin1,pcos2,psin2
     +,        qpterp1,LWtoa,LWtoacl,LWupsfc
     +,        LWdnsfc,LWdnsfccl,LWnetsfc,LWnetsfccl,LWtom,LWtomcl
     +,        SWnettoa,SWnettoacl,SWdnsfc,SWdnsfccl,SWnetsfc,SWnetsfccl
     +,        SWinc,SWtom,SWtomcl,radvbc,radvdc,radnbc,radndc,radvbcc
     +,        radvdcc,radnbcc,radndcc,alvdf,alvdr,alndf,alndr,atl,atlcl
     +,        asl,aslcl)
!     print*,'---- exit subroutine BUGSradQPs:'
#endif

!     print*,'---- exit subroutine bugs_rad:'
!      print *, umco2,umch4,umn2o

      return
      end subroutine bugs_rad

!----------------------------------------------------------------------
