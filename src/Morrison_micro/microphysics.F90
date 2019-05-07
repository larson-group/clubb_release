module microphysics

! main interface to Morrison microphysics.
! original implementation by Peter Blossey, UW

use params, only: lcond, lsub, fac_cond, fac_sub, ggr

use grid, only: nx,ny,nzm,nz, &  !grid dimensions; nzm = nz-1 # of scalar lvls
     dimx1_s,dimx2_s,dimy1_s,dimy2_s, & ! actual scalar-array dimensions in x,y
     dz, adz, dostatis, masterproc, &
     doSAMconditionals, dosatupdnconditionals

use vars, only: pres, rho, dtn, w, t, tlatqi, condavg_mask, &
     ncondavg, condavgname, condavglongname
use params, only: doprecip, docloud, doclubb

use module_mp_GRAUPEL, only: GRAUPEL_INIT, M2005MICRO_GRAUPEL, &
      doicemicro, &         ! use ice species (snow/cloud ice/graupel)
      dograupel, &          ! use graupel
      dohail, &             ! use graupel
      dosb_warm_rain, &     ! use Seifert & Beheng (2001) warm rain parameterization
      dopredictNc, &        ! prediction of cloud droplet number
      aerosol_mode, &   ! specify two modes of (sulfate) aerosol
      dosubgridw, &         ! input estimate of subgrid w to microphysics
      doarcticicenucl,&     ! use arctic parameter values for ice nucleation
      docloudedgeactivation,&! activate droplets at cloud edges as well as base
      Nc0,            &     ! initial/specified cloud droplet number conc (#/cm3)
      ccnconst, ccnexpnt, & ! parameters for aerosol_mode=1 (powerlaw CCN)
      aer_rm1, aer_rm2, &   ! two modes of aerosol for aerosol_mode=2
      aer_n1, aer_n2, &     ! rm=geometric mean radius (um), n=aerosol conc. (#/cm3)
      aer_sig1, aer_sig2, & ! sig=geom standard deviation of aerosol size distn.
      dofix_pgam, pgam_fixed ! option to specify pgam (exponent of cloud water's gamma distn)

implicit none

logical :: isallocatedMICRO = .false.

integer :: nmicro_fields ! total number of prognostic water vars

#ifdef SILHS
real, allocatable, target, dimension(:,:,:,:) :: micro_field  ! holds mphys quantities
#else
real, allocatable, dimension(:,:,:,:) :: micro_field  ! holds mphys quantities
#endif

! indices of water quantities in micro_field, e.g. qv = micro_field(:,:,:,iqv)
integer :: iqv, iqcl, iqci, iqr, iqs, iqg, incl, inci, inr, ins, ing
integer :: index_water_vapor ! separate water vapor index used by SAM

real, allocatable, dimension(:) :: lfac
integer, allocatable, dimension(:) :: flag_wmass, flag_precip, flag_number
integer, allocatable, dimension(:) :: flag_micro3Dout

integer, parameter :: index_cloud_ice = -1 ! historical variable (don't change)

real, allocatable, dimension(:,:,:) :: fluxbmk, fluxtmk !surface/top fluxes
real, allocatable, dimension(:,:,:) :: reffc, reffi
#ifdef SILHS
real, allocatable, target, dimension(:,:,:) :: cloudliq
#else
real, allocatable, dimension(:,:,:) :: cloudliq
#endif

real, allocatable, dimension(:,:) :: & ! statistical arrays
     mkwle, & ! resolved vertical flux
     mkwsb, & ! SGS vertical flux
     mksed, & ! sedimentation vertical flux
     mkadv, & ! tendency due to vertical advection
     mkdiff, &! tendency due to vertical diffusion
     mklsadv, & ! tendency due to large-scale vertical advection
     mfrac, & ! fraction of domain with microphysical quantity > 1.e-6
     stend, & ! tendency due to sedimentation
     mtend, & ! tendency due to microphysical processes (other than sedimentation)
     mstor, & ! storage terms of microphysical variables 
     trtau    ! optical depths of various species

real, allocatable, dimension(:) :: tmtend

real :: sfcpcp, sfcicepcp

! arrays with names/units for microphysical outputs in statistics.
character*3, allocatable, dimension(:) :: mkname
character*80, allocatable, dimension(:) :: mklongname
character*10, allocatable, dimension(:) :: mkunits
real, allocatable, dimension(:) :: mkoutputscale
logical douse_reffc, douse_reffi

! You can also have some additional, diagnostic, arrays, for example, total
! nonprecipitating cloud water, etc:

!bloss: array which holds temperature tendency due to microphysics
real, allocatable, dimension(:,:,:), SAVE :: tmtend3d

real, dimension(:,:,:), pointer :: conc, qn

#ifdef UWM_STATS
integer :: &
  idx_cor_w_ncl, idx_cor_w_qr, idx_cor_w_nr, idx_cor_w_qci, idx_cor_w_nci, idx_cor_w_qs, idx_cor_w_ns, &
  idx_cor_w_qg, idx_cor_w_ng, idx_cor_ncl_qr, idx_cor_ncl_nr, idx_cor_ncl_qci, idx_cor_ncl_nci, &
  idx_cor_ncl_qs, idx_cor_ncl_ns, idx_cor_ncl_qg, idx_cor_ncl_ng, idx_cor_qr_nr, idx_cor_qr_qci, &
  idx_cor_qr_nci, idx_cor_qr_qs, idx_cor_qr_ns, idx_cor_qr_qg, idx_cor_qr_ng, idx_cor_nr_qci, &
  idx_cor_nr_nci, idx_cor_nr_qs, idx_cor_nr_ns, idx_cor_nr_qg, idx_cor_nr_ng, idx_cor_qci_nci, &
  idx_cor_qci_qs, idx_cor_qci_ns, idx_cor_qci_qg, idx_cor_qci_ng, idx_cor_nci_qs, idx_cor_nci_ns, &
  idx_cor_nci_qg, idx_cor_nci_ng, idx_cor_qs_ns, idx_cor_qs_qg, idx_cor_qs_ng, idx_cor_ns_qg, &
  idx_cor_ns_ng, idx_cor_qg_ng,&
  
  !Fractions
  idx_cld_frac, idx_rain_frac, idx_cldic_frac, idx_snw_frac, idx_grpl_frac 
#endif /*UWM_STATS*/

! Microphysical process rates  +++mhwang
real, allocatable, dimension(:) :: mPRC, mPRA, mPSMLT, mEVPMS, mPRACS, mEVPMG, mPRACG, mPRE, mPGMLT, &
                        mMNUCCC, mPSACWS, mPSACWI, mQMULTS, mQMULTG, mPSACWG, mPGSACW, &
                        mPRD, mPRCI, mPRAI, mQMULTR, mQMULTRG, mMNUCCD, mPRACI, mPRACIS, mEPRD, &
                        mMNUCCR, mPIACR, mPIACRS, mPGRACS, &
                        mPRDS, mEPRDS, mPSACR, &
                        mPRDG, mEPRDG, &
                        mNPRC1, mNRAGG, mNPRACG, mNSUBR, mNSMLTR, &
                        mNGMLTR, mNPRACS, mNNUCCR, mNIACR, mNIACRS, mNGRACS, &
                        mNSMLTS, mNSAGG, mNPRCI, mNSCNG, mNSUBS, mPCC, &
                        mNNUCCC, mNPSACWS, mNPRA, mNPRC, mNPSACWI, &
                        mNPSACWG, mNPRAI, mNMULTS, mNMULTG, mNMULTR, &
                        mNMULTRG, mNNUCCD, mNSUBI, mNGMLTG, mNSUBG, mNACT, &
                        mSIZEFIX_NR, mSIZEFIX_NC, mSIZEFIX_NI, mSIZEFIX_NS, mSIZEFIX_NG, &
                        mNEGFIX_NR, mNEGFIX_NC, mNEGFIX_NI, mNEGFIX_NS, mNEGFIX_NG, &
                        mNIM_MORR_CL, mQC_INST, mQR_INST, mQI_INST, mQS_INST, mQG_INST, &
                        mNC_INST, mNR_INST, mNI_INST, mNS_INST, mNG_INST  




CONTAINS

!----------------------------------------------------------------------
!!! Read microphysical options from prm file and allocate variables
!
subroutine micro_setparm()
  use vars
#ifdef CLUBB
  use module_mp_graupel, only: NNUCCD_REDUCE_COEF, NNUCCC_REDUCE_COEF
#endif
  implicit none

  integer ierr, ios, ios_missing_namelist, place_holder
  
   NAMELIST /MICRO_M2005/ &
#ifdef CLUBB
      NNUCCD_REDUCE_COEF, NNUCCC_REDUCE_COEF, &
#endif
      doicemicro, &         ! use ice species (snow/cloud ice/graupel)
      dograupel, &          ! use graupel
      dohail, &             ! graupel species has qualities of hail
      dosb_warm_rain, &     ! use Seifert & Beheng (2001) warm rain parameterization in place of KK(2000)
      dopredictNc, &        ! prediction of cloud droplet number
      aerosol_mode, &   ! specify two modes of (sulfate) aerosol
      dosubgridw, &         ! input estimate of subgrid w to microphysics
      doarcticicenucl,&     ! use arctic parameter values for ice nucleation
      docloudedgeactivation,&! activate droplets at cloud edges as well as base
      Nc0,            &     ! initial/specified cloud droplet number conc (#/cm3)
      ccnconst, ccnexpnt, & ! parameters for aerosol_mode=1 (powerlaw CCN)
      aer_rm1, aer_rm2, &   ! two modes of aerosol for aerosol_mode=2
      aer_n1, aer_n2, &     ! rm=geometric mean radius (um), n=aerosol conc. (#/cm3)
      aer_sig1, aer_sig2, & ! sig=geom standard deviation of aerosol size distn.
      dofix_pgam, pgam_fixed, & ! option to specify pgam (exponent of cloud water's gamma distn)
      douse_reffc, &        ! use computed effective radius in radiation computation
      douse_reffi           ! use computed effective ice size in radiation computation

   !bloss: Create dummy namelist, so that we can figure out error code
   !       for a mising namelist.  This lets us differentiate between
   !       missing namelists and those with an error within the namelist.
   NAMELIST /BNCUIODSBJCB/ place_holder
#ifdef SILHS
   integer :: k
#endif

   ! define default values for namelist variables
   doicemicro = .true.        ! use ice
   dograupel = .true.         ! use graupel
   dohail = .false.           ! graupel species has properties of graupel
   dosb_warm_rain = .false.   ! use KK (2000) warm rain scheme by default
   dopredictNc = .false.       ! prognostic cloud droplet number (default = .false. !bloss Apr 09)
#ifdef CLUBB
   dosubgridw = .true.        ! Use clubb's w'^2 for sgs w
   aerosol_mode = 2           ! use lognormal CCN relationship
   doarcticicenucl = .false.  ! use mid-latitude parameters
#else
   aerosol_mode = 1           ! use powerlaw CCN relationship
   dosubgridw=.false.         ! don't bother with estimating w_sgs for now
   doarcticicenucl = .false.  ! use mid-latitude parameters
#endif
   docloudedgeactivation  = .false. ! activate droplets at cloud base, not edges
   douse_reffc = .false.  ! use computed effective radius in rad computations?
   douse_reffi = .false.  ! use computed effective radius in rad computations?

   Nc0 = 100. ! default droplet number concentration
   
   ccnconst = 120.            ! maritime value (/cm3), adapted from Rasmussen 
   ccnexpnt = 0.4             !   et al (2002) by Hugh Morrison et al.  Values
                              !   of 1000. and 0.5 suggested for continental
#ifdef CLUBB
   ! Aerosol for RF02 from Mikhail
   aer_rm1  = 0.011
   aer_sig1 = 1.2
   aer_n1   = 125.
   aer_rm2  = 0.06
   aer_sig2 = 1.7
   aer_n2   = 65.

#else
   aer_rm1 = 0.052           ! two aerosol mode defaults from MPACE (from Hugh)
   aer_sig1 = 2.04
   aer_n1 = 72.2
   aer_rm2 = 1.3
   aer_sig2 = 2.5
   aer_n2 = 1.8
#endif /* CLUBB */
   
   dofix_pgam = .false.
   pgam_fixed = 5. ! middle range value -- corresponds to radius dispersion ~ 0.4

  !----------------------------------
  !  Read namelist for microphysics options from prm file:
  !------------
  open(55,file='./'//trim(case)//'/prm', status='old',form='formatted') 
  
  !bloss: get error code for missing namelist (by giving the name for
  !       a namelist that doesn't exist in the prm file).
  read (UNIT=55,NML=BNCUIODSBJCB,IOSTAT=ios_missing_namelist)
  rewind(55) !note that one must rewind before searching for new namelists

  !bloss: read in MICRO_M2005 namelist
  read (55,MICRO_M2005,IOSTAT=ios)

  if (ios.ne.0) then
     !namelist error checking
     if(ios.ne.ios_missing_namelist) then
        write(*,*) '****** ERROR: bad specification in MICRO_M2005 namelist'
        call task_abort()
     elseif(masterproc) then
        write(*,*) '****************************************************'
        write(*,*) '****** No MICRO_M2005 namelist in prm file *********'
        write(*,*) '****************************************************'
     end if
  end if
  close(55)

   if(.not.doicemicro) dograupel=.false.

   if(dohail.and..NOT.dograupel) then
      if(masterproc) write(*,*) 'dograupel must be .true. for dohail to be used.'
      call task_abort()
   end if

#ifndef CLUBB /* Disable this for UWM simulations */
   ! write namelist values out to file for documentation
   if(masterproc) then
      open(unit=55,file='./'//trim(case)//'/'//trim(case)//'_'//trim(caseid)//'.namelists', form='formatted', position='append')
      write (unit=55,nml=MICRO_M2005,IOSTAT=ios)
      write(55,*) ' '
      close(unit=55)
   end if
#endif

   ! scale values of parameters for m2005micro
   aer_rm1 = 1.e-6*aer_rm1 ! convert from um to m
   aer_rm2 = 1.e-6*aer_rm2 
   aer_n1 = 1.e6*aer_n1 ! convert from #/cm3 to #/m3
   aer_n2 = 1.e6*aer_n2
#ifdef MICRO_RESTART /* Save all fields to allow for restarting with micro enabled */
  nmicro_fields = 10
  iqv  = 1   ! total water (vapor + cloud liq) mass mixing ratio [kg H2O / kg dry air]
  incl = 2  ! cloud water number mixing ratio [#/kg dry air]
  iqr  = 3   ! rain mass mixing ratio [kg H2O / kg dry air]
  inr  = 4   ! rain number mixing ratio [#/kg dry air]
  iqci = 5  ! cloud ice mass mixing ratio [kg H2O / kg dry air]
  inci = 6  ! cloud ice number mixing ratio [#/kg dry air]
  iqs  = 7   ! snow mass mixing ratio [kg H2O / kg dry air]
  ins  = 8   ! snow number mixing ratio [#/kg dry air]
  iqg  = 9  ! graupel mass mixing ratio [kg H2O / kg dry air]
  ing  = 10  ! graupel number mixing ratio [#/kg dry air]
#else
  nmicro_fields = 1 ! start with water vapor and cloud water mass mixing ratio
#ifdef CLUBB
  if(docloud.or.doclubb) then
#else
  if(docloud) then
#endif
!bloss/qt     nmicro_fields = nmicro_fields + 1 ! add cloud water mixing ratio
     if(dopredictNc) nmicro_fields = nmicro_fields + 1 ! add cloud water number concentration (if desired)
  end if
  if(doprecip)    nmicro_fields = nmicro_fields + 2 ! add rain mass and number (if desired)
  if(doicemicro)  nmicro_fields = nmicro_fields + 4 ! add snow and cloud ice number and mass (if desired)
  if(dograupel)   nmicro_fields = nmicro_fields + 2 ! add graupel mass and number (if desired)

  ! specify index of various quantities in micro_field array
  !  *** note that not all of these may be used if(.not.doicemicro) ***
  iqv = 1   ! total water (vapor + cloud liq) mass mixing ratio [kg H2O / kg dry air]
!bloss/qt  iqcl = 2  ! cloud water mass mixing ratio [kg H2O / kg dry air]
  
!bloss/qt: cloud liquid water no longer prognosed
  if(dopredictNc) then
     incl = 2  ! cloud water number mixing ratio [#/kg dry air]
     iqr = 3   ! rain mass mixing ratio [kg H2O / kg dry air]
     inr = 4   ! rain number mixing ratio [#/kg dry air]
     iqci = 5  ! cloud ice mass mixing ratio [kg H2O / kg dry air]
     inci = 6  ! cloud ice number mixing ratio [#/kg dry air]
     iqs = 7   ! snow mass mixing ratio [kg H2O / kg dry air]
     ins = 8   ! snow number mixing ratio [#/kg dry air]
     iqg = 9  ! graupel mass mixing ratio [kg H2O / kg dry air]
     ing = 10  ! graupel number mixing ratio [#/kg dry air]
  else
     iqr = 2   ! rain mass mixing ratio [kg H2O / kg dry air]
     inr = 3   ! rain number mixing ratio [#/kg dry air]
     iqci = 4  ! cloud ice mass mixing ratio [kg H2O / kg dry air]
     inci = 5  ! cloud ice number mixing ratio [#/kg dry air]
     iqs = 6   ! snow mass mixing ratio [kg H2O / kg dry air]
     ins = 7   ! snow number mixing ratio [#/kg dry air]
     iqg = 8   ! graupel mass mixing ratio [kg H2O / kg dry air]
     ing = 9  ! graupel number mixing ratio [#/kg dry air]
#ifdef UWM_STATS
     ! Even if dopredictNc is false, we want to monitor the output of tmpncl
     ! for certainty. 
     incl = 10  ! cloud water number mixing ratio [#/kg dry air]
     nmicro_fields = nmicro_fields + 1
#endif
  end if
#endif /* MICRO_RESTART */
  ! stop if icemicro is specified without precip -- we don't support this right now.
  if((doicemicro).and.(.not.doprecip)) then
     if(masterproc) write(*,*) 'Morrison 2005 Microphysics does not support both doice and .not.doprecip'
     call task_abort()
  end if
  index_water_vapor = iqv ! set SAM water vapor flag

  if(.not.isallocatedMICRO) then
     ! allocate microphysical variables
     allocate(micro_field(dimx1_s:dimx2_s,dimy1_s:dimy2_s,nzm,nmicro_fields), &
          fluxbmk(nx,ny,nmicro_fields), fluxtmk(nx,ny,nmicro_fields), &
          reffc(nx,ny,nzm), reffi(nx,ny,nzm), &
          mkwle(nz,nmicro_fields), mkwsb(nz,nmicro_fields), &
          mkadv(nz,nmicro_fields), mkdiff(nz,nmicro_fields), &
          mklsadv(nz,nmicro_fields), &
          stend(nzm,nmicro_fields), mtend(nzm,nmicro_fields), &
          mfrac(nzm,nmicro_fields), trtau(nzm,nmicro_fields), &
          mksed(nzm,nmicro_fields), tmtend(nzm), &
          mstor(nzm,nmicro_fields),  &
          cloudliq(nx,ny,nzm), &
          tmtend3d(nx,ny,nzm), flag_micro3Dout(nmicro_fields), &
          flag_wmass(nmicro_fields), flag_precip(nmicro_fields), &
          flag_number(nmicro_fields), lfac(nmicro_fields), &
          mkname(nmicro_fields), mklongname(nmicro_fields), &
          mkunits(nmicro_fields), mkoutputscale(nmicro_fields), STAT=ierr)

 
! Microphysical process rates +++mhwang 
      allocate(mPRC(nzm), mPRA(nzm), mPSMLT(nzm), mEVPMS(nzm), mPRACS(nzm), mEVPMG(nzm), mPRACG(nzm), mPRE(nzm), mPGMLT(nzm),&
                        mMNUCCC(nzm), mPSACWS(nzm), mPSACWI(nzm), mQMULTS(nzm), mQMULTG(nzm), mPSACWG(nzm), mPGSACW(nzm), & 
                        mPRD(nzm), mPRCI(nzm), mPRAI(nzm), mQMULTR(nzm), mQMULTRG(nzm), mMNUCCD(nzm), mPRACI(nzm), &
                        mPRACIS(nzm), mEPRD(nzm), mMNUCCR(nzm), mPIACR(nzm), mPIACRS(nzm), mPGRACS(nzm), &
                        mPRDS(nzm), mEPRDS(nzm), mPSACR(nzm), mPRDG(nzm), mEPRDG(nzm), mNPRC1(nzm), mNRAGG(nzm), &
                        mNPRACG(nzm), mNSUBR(nzm), mNSMLTR(nzm), mNGMLTR(nzm), mNPRACS(nzm), &
                        mNNUCCR(nzm), mNIACR(nzm), mNIACRS(nzm), mNGRACS(nzm), mNSMLTS(nzm), &
                        mNSAGG(nzm), mNPRCI(nzm), mNSCNG(nzm), mNSUBS(nzm), mPCC(nzm), &
                        mNNUCCC(nzm), mNPSACWS(nzm), mNPRA(nzm), mNPRC(nzm), mNPSACWI(nzm), &
                        mNPSACWG(nzm), mNPRAI(nzm), mNMULTS(nzm), mNMULTG(nzm), mNMULTR(nzm), &
                        mNMULTRG(nzm), mNNUCCD(nzm), mNSUBI(nzm), mNGMLTG(nzm), mNSUBG(nzm), mNACT(nzm), & 
                        mSIZEFIX_NR(nzm), mSIZEFIX_NC(nzm), mSIZEFIX_NI(nzm), mSIZEFIX_NS(nzm), &
                        mSIZEFIX_NG(nzm), mNEGFIX_NR(nzm), mNEGFIX_NC(nzm), mNEGFIX_NI(nzm), &
                        mNEGFIX_NS(nzm), mNEGFIX_NG(nzm), mNIM_MORR_CL(nzm), &
                        mQC_INST(nzm), mQR_INST(nzm), mQI_INST(nzm), mQS_INST(nzm), mQG_INST(nzm), &
                        mNC_INST(nzm), mNR_INST(nzm), mNI_INST(nzm), mNS_INST(nzm), mNG_INST(nzm), &
                        STAT=ierr)

     if(ierr.ne.0) then
        write(*,*) 'Failed to allocate microphysical arrays on proc ', rank
        call task_abort()
     else
        isallocatedMICRO = .true.
     end if

     ! zero out statistics variables associated with cloud ice sedimentation
     !   in Marat's default SAM microphysics
     tlatqi = 0.

     ! initialize these arrays
     micro_field = 0.
     cloudliq = 0. !bloss/qt: auxially cloud liquid water variable, analogous to qn in MICRO_SAM1MOM
     fluxbmk = 0.
     fluxtmk = 0.
     mkwle = 0.
     mkwsb = 0.
     mkadv = 0.
     mkdiff = 0.
     mklsadv = 0.
  mstor =0.

    ! initialize flag arrays to all mass, no number, no precip
     flag_wmass = 1
     flag_number = 0
     flag_precip = 0
     flag_micro3Dout = 0

  end if

  compute_reffc = douse_reffc
  compute_reffi = douse_reffi
#ifdef SILHS
  if ( dopredictnc ) then
    conc => micro_field(1:nx,1:ny,1:nzm,incl)
  else
    allocate( conc(nx,ny,nzm) )
  end if
  qn => cloudliq
#endif 
end subroutine micro_setparm

!----------------------------------------------------------------------
!!! Initialize microphysics:
!
! this one is guaranteed to be called by SAM at the 
!   beginning of each run, initial or restart:
subroutine micro_init()

  use vars
  
  implicit none
  
  real, dimension(nzm) :: qc0, qi0

! Commented out by dschanen UWM 23 Nov 2009 to avoid a linking error
! real, external :: satadj_water 
  integer :: k
  integer :: i, j, n

#ifndef MICRO_RESTART
  ! initialize flag arrays
  if(dopredictNc) then
     ! Cloud droplet number concentration is a prognostic variable
     if(doicemicro) then
        if(dograupel) then
          !bloss/qt: qt, Nc, qr, Nr, qi, Ni, qs, Ns, qg, Ng
           flag_wmass  = (/1,0,1,0,1,0,1,0,1,0/)
           flag_precip = (/0,0,1,1,0,0,1,1,1,1/)
           flag_number = (/0,1,0,1,0,1,0,1,0,1/)
        else
          !bloss/qt: qt, Nc, qr, Nr, qi, Ni, qs, Ns
           flag_wmass  = (/1,0,1,0,1,0,1,0/)
           flag_precip = (/0,0,1,1,0,0,1,1/)
           flag_number = (/0,1,0,1,0,1,0,1/)
        end if
     else
        if(doprecip) then
          !bloss/qt: qt, Nc, qr, Nr
           flag_wmass  = (/1,0,1,0/)
           flag_precip = (/0,0,1,1/)
           flag_number = (/0,1,0,1/)
        else
          !bloss/qt: qt, Nc
           flag_wmass  = (/1,0/)
           flag_precip = (/0,0/)
           flag_number = (/0,1/)
        end if
     end if
  else
     ! Cloud droplet number concentration is NOT a prognostic variable
     if(doicemicro) then
        if(dograupel) then
          !bloss/qt: qt, qr, Nr, qi, Ni, qs, Ns, qg, Ng
           flag_wmass  = (/1,1,0,1,0,1,0,1,0/)
           flag_precip = (/0,1,1,0,0,1,1,1,1/)
           flag_number = (/0,0,1,0,1,0,1,0,1/)
        else
          !bloss/qt: qt, qr, Nr, qi, Ni, qs, Ns
           flag_wmass  = (/1,1,0,1,0,1,0/)
           flag_precip = (/0,1,1,0,0,1,1/)
           flag_number = (/0,0,1,0,1,0,1/)
        end if
     else
        if(doprecip) then
          !bloss/qt: qt, qr, Nr
           flag_wmass  = (/1,1,0/)
           flag_precip = (/0,1,1/)
           flag_number = (/0,0,1/)
        else
          !bloss/qt: only total water variable is needed for no-precip, 
          !            fixed droplet number, warm cloud and no cloud simulations.
           flag_wmass  = (/1/)
           flag_precip = (/0/)
           flag_number = (/0/)
        end if
     end if
  end if
#else
flag_wmass  = (/1,0,1,0,1,0,1,0,1,0/)
flag_precip = (/0,0,1,1,0,0,1,1,1,1/)
flag_number = (/0,1,0,1,0,1,0,1,0,1/)
#endif /* MICRO_RESTART */
  ! output all microphysical fields to 3D output files if using more than
  !   just docloud.  Otherwise, rely on basic SAM outputs
#ifdef CLUBB
  if((docloud.OR.doclubb).AND.(doprecip.OR.dopredictNc)) then
#else
  if(docloud.AND.(doprecip.OR.dopredictNc)) then
#endif
     flag_micro3Dout = 1
  end if

  ! initialize factor for latent heat
  lfac(:) = 1. ! use one as default for number species
  lfac(iqv) = lcond
!bloss/qt  if(docloud) lfac(iqcl) = lcond
  if(doprecip) lfac(iqr) = lcond
  if(doicemicro) then
     lfac(iqci) = lsub
     lfac(iqs) = lsub
     if(dograupel) lfac(iqg) = lsub
  end if

  call graupel_init() ! call initialization routine within mphys module

  if(nrestart.eq.0) then

 ! compute initial profiles of liquid water - M.K.
      call satadj_liquid(nzm,tabs0,q0,qc0,pres*100.)

     ! initialize microphysical quantities
     q0 = q0 + qc0
     do k = 1,nzm
        micro_field(:,:,k,iqv) = q0(k)
        cloudliq(:,:,k) = qc0(k)
        tabs(:,:,k) = tabs0(k)
     end do
     if(dopredictNc) then ! initialize concentration somehow...
       do k = 1,nzm
         if(q0(k).gt.0.) then
            micro_field(:,:,k,incl) = 0.5*ccnconst*1.e6
         end if
       end do
     end if
#ifdef CLUBB
     if(docloud.or.doclubb)  call micro_diagnose()   ! leave this line here
#else
     if(docloud) call micro_diagnose()   ! leave this here
#endif


  end if

  ! set up names, units and scales for these microphysical quantities
  mkname(iqv) = 'QTO'
  mklongname(iqv) = 'TOTAL WATER (VAPOR + CLOUD LIQUID)'
  mkunits(iqv) = 'g/kg'
  mkoutputscale(iqv) = 1.e3

#ifdef MICRO_RESTART
  mkname(incl) = 'NC'
  mklongname(incl) = 'CLOUD WATER NUMBER CONCENTRATION'
  mkunits(incl) = '#/cm3'
  mkoutputscale(incl) = 1.e-6

  mkname(iqr) = 'QR'
  mklongname(iqr) = 'RAIN'
  mkunits(iqr) = 'g/kg'
  mkoutputscale(iqr) = 1.e3

  mkname(inr) = 'NR'
  mklongname(inr) = 'RAIN NUMBER CONCENTRATION'
  mkunits(inr) = '#/cm3'
  mkoutputscale(inr) = 1.e-6

  mkname(iqci) = 'QI'
  mklongname(iqci) = 'CLOUD ICE'
  mkunits(iqci) = 'g/kg'
  mkoutputscale(iqci) = 1.e3

  mkname(inci) = 'NI'
  mklongname(inci) = 'CLOUD ICE NUMBER CONCENTRATION'
  mkunits(inci) = '#/cm3'
  mkoutputscale(inci) = 1.e-6

  mkname(iqs) = 'QS'
  mklongname(iqs) = 'SNOW'
  mkunits(iqs) = 'g/kg'
  mkoutputscale(iqs) = 1.e3

  mkname(ins) = 'NS'
  mklongname(ins) = 'SNOW NUMBER CONCENTRATION'
  mkunits(ins) = '#/cm3'
  mkoutputscale(ins) = 1.e-6

  mkname(iqg) = 'QG'
  mklongname(iqg) = 'GRAUPEL'
  mkunits(iqg) = 'g/kg'
  mkoutputscale(iqg) = 1.e3

  mkname(ing) = 'NG'
  mklongname(ing) = 'GRAUPEL NUMBER CONCENTRATION'
  mkunits(ing) = '#/cm3'
  mkoutputscale(ing) = 1.e-6
#else
#ifdef CLUBB
  if(docloud.or.doclubb) then
#else
  if(docloud) then
#endif
!bloss/qt     mkname(iqcl) = 'QC'
!bloss/qt     mklongname(iqcl) = 'CLOUD WATER'
!bloss/qt     mkunits(iqcl) = 'g/kg'
!bloss/qt     mkoutputscale(iqcl) = 1.e3

! We want to output NC even if dopredictNc is false for certainty. 
#ifndef UWM_STATS
     if(dopredictNc) then
#endif
        mkname(incl) = 'NC'
        mklongname(incl) = 'CLOUD WATER NUMBER CONCENTRATION'
        mkunits(incl) = '#/cm3'
        mkoutputscale(incl) = 1.e-6
#ifndef UWM_STATS
     end if
#endif
  end if

  if(doprecip) then
     mkname(iqr) = 'QR'
     mklongname(iqr) = 'RAIN'
     mkunits(iqr) = 'g/kg'
     mkoutputscale(iqr) = 1.e3

     mkname(inr) = 'NR'
     mklongname(inr) = 'RAIN NUMBER CONCENTRATION'
     mkunits(inr) = '#/cm3'
     mkoutputscale(inr) = 1.e-6
  end if

  if(doicemicro) then
     mkname(iqci) = 'QI'
     mklongname(iqci) = 'CLOUD ICE'
     mkunits(iqci) = 'g/kg'
     mkoutputscale(iqci) = 1.e3

     mkname(inci) = 'NI'
     mklongname(inci) = 'CLOUD ICE NUMBER CONCENTRATION'
     mkunits(inci) = '#/cm3'
     mkoutputscale(inci) = 1.e-6

     mkname(iqs) = 'QS'
     mklongname(iqs) = 'SNOW'
     mkunits(iqs) = 'g/kg'
     mkoutputscale(iqs) = 1.e3

     mkname(ins) = 'NS'
     mklongname(ins) = 'SNOW NUMBER CONCENTRATION'
     mkunits(ins) = '#/cm3'
     mkoutputscale(ins) = 1.e-6

     if(dograupel) then
        mkname(iqg) = 'QG'
        mklongname(iqg) = 'GRAUPEL'
        mkunits(iqg) = 'g/kg'
        mkoutputscale(iqg) = 1.e3

        mkname(ing) = 'NG'
        mklongname(ing) = 'GRAUPEL NUMBER CONCENTRATION'
        mkunits(ing) = '#/cm3'
        mkoutputscale(ing) = 1.e-6
     end if
end if
#endif /* MICRO_RESTART */

! set mstor to be the inital microphysical mixing ratios    
     do n=1, nmicro_fields
     do k=1, nzm
         mstor(k, n) = SUM(micro_field(1:nx,1:ny,k,n))
     end do
     end do

end subroutine micro_init

!----------------------------------------------------------------------
!!! fill-in surface and top boundary fluxes:
!
! Obviously, for liquid/ice water variables those fluxes are zero. They are not zero
! only for water vapor variable and, possibly, for CCN and IN if you have those.

subroutine micro_flux()

use vars, only: fluxbq, fluxtq
#ifdef CLUBB
use params, only: doclubb, doclubb_sfc_fluxes
#endif

fluxbmk(:,:,:) = 0. ! initialize all fluxes at surface to zero
fluxtmk(:,:,:) = 0. ! initialize all fluxes at top of domain to zero
#ifdef CLUBB
if ( doclubb .and. doclubb_sfc_fluxes ) then
  fluxbmk(:,:,index_water_vapor) = 0.0 ! surface qv (latent heat) flux
else
  fluxbmk(:,:,index_water_vapor) = fluxbq(:,:) ! surface qv (latent heat) flux
end if
#else
fluxbmk(:,:,index_water_vapor) = fluxbq(:,:) ! surface qv (latent heat) flux
#endif
fluxtmk(:,:,index_water_vapor) = fluxtq(:,:) ! top of domain qv flux

end subroutine micro_flux

!----------------------------------------------------------------------
!!! compute local microphysics processes (beyond advection and SGS diffusion):
!
!  This is the place where the condensation/sublimation, accretion, coagulation, freezing,
!  melting, etc., that is  all the microphysics processes except for the spatial transport happen.

! IMPORTANT: You need to use the thermodynamic constants like specific heat, or
! specific heat of condensation, gas constant, etc, the same as in file params.f90
! Also, you should assume that the conservative thermodynamic variable during these
! proceses is the liquid/ice water static energy: t = tabs + gz - Lc (qc+qr) - Ls (qi+qs+qg) 
! It should not be changed during all of your point microphysical processes!

subroutine micro_proc()

use params, only: fac_cond, fac_sub, rgas
use grid, only: z, zi
use vars, only: t,  gamaz, precsfc, precflux, qpfall, tlat, prec_xy, &
     nstep, nstatis, icycle, total_water_prec

#ifdef UWM_STATS
!weberjk(UWM), to compute budget statistics on q2, t2, tw, qw, include the
!effect of precipitation (microphysics)    
use vars, only: t2leprec, q2leprec, qwleprec, twleprec
use module_mp_GRAUPEL, only: Nc0
#endif /*UWM_STATS*/

#ifdef CLUBB
use params, only: doclubb, docloud, dosmoke
use grid, only: nz
use error_code, only: clubb_at_least_debug_level
use fill_holes, only: fill_holes_vertical
use clubbvars, only: wp2, cloud_frac, rho_ds_zt, rho_ds_zm ! are used, but not modified here
use vars, only: qcl ! Used here and updated in micro_diagnose
use vars, only: prespot ! exner^-1
use module_mp_GRAUPEL, only: &
  cloud_frac_thresh ! Threshold for using sgs cloud fraction to weight 
                    ! microphysical quantities [%]
use clubb_precision, only: core_rknd

#endif
#ifdef SILHS
use parameters_microphys, only: &
  LH_microphys_type, & ! Variables
  LH_microphys_interactive, &
  LH_microphys_non_interactive, &
  LH_microphys_disabled
#endif

real, dimension(nzm) :: &
     tmpqcl, tmpqci, tmpqr, tmpqs, tmpqg, tmpqv, &
     tmpncl, tmpnci, tmpnr, tmpns, tmpng,  &
     tmpw, tmpwsub, tmppres, tmpdz, tmptabs, &
     tmtend1d, &
     mtendqcl, mtendqci, mtendqr, mtendqs, mtendqg, mtendqv, &
     mtendncl, mtendnci, mtendnr, mtendns, mtendng,  &
     stendqcl, stendqci, stendqr, stendqs, stendqg, stendqv, &
     stendncl, stendnci, stendnr, stendns, stendng,  &
     effg1d, effr1d, effs1d, effc1d, effi1d

real df(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
real f0(nzm),df0(nzm)

#ifdef CLUBB
real(kind=core_rknd), dimension(nz) :: &
     qv_clip, qcl_clip

!real, dimension(nzm) :: & These variables were renamed PRC, PRA, and PRE respectively as per ticket #552 and are declared below
!  qr_auto, qr_accr, qr_evap
real, dimension(nzm) :: &
  qr_auto, qr_accr, qr_evap, rain_vel
real, dimension(nzm) :: cloud_frac_in
#endif /*CLUBB*/

! Microphysical rate +++mwang  
real, dimension(nzm) :: PRC, PRA, PSMLT, EVPMS, PRACS, EVPMG, PRACG, PRE, PGMLT,  & 
                        MNUCCC, PSACWS, PSACWI, QMULTS, QMULTG, PSACWG, PGSACW, & 
                        PRD, PRCI, PRAI, QMULTR, QMULTRG, MNUCCD, PRACI, PRACIS, EPRD, & 
                        MNUCCR, PIACR, PIACRS, PGRACS, & 
                        PRDS, EPRDS, PSACR, & 
                        PRDG, EPRDG, &
                        NPRC1, NRAGG, NPRACG, NSUBR, NSMLTR, NGMLTR, &
                        NPRACS, NNUCCR, NIACR, NIACRS, NGRACS, &
                        NSMLTS, NSAGG, NPRCI, NSCNG, NSUBS, PCC, NNUCCC, &
                        NPSACWS, NPRA, NPRC, NPSACWI, NPSACWG, NPRAI, &
                        NMULTS, NMULTG, NMULTR, NMULTRG, NNUCCD, NSUBI, NGMLTG, NSUBG, NACT, &
                        SIZEFIX_NR, SIZEFIX_NC, SIZEFIX_NI, SIZEFIX_NS, SIZEFIX_NG, &
                        NEGFIX_NR, NEGFIX_NC, NEGFIX_NI, NEGFIX_NS, NEGFIX_NG, &
                        NIM_MORR_CL, QC_INST, QR_INST, QI_INST, QS_INST, QG_INST, &
                        NC_INST, NR_INST, NI_INST, NS_INST, NG_INST 



real, dimension(nzm,nmicro_fields) :: stend1d, mtend1d
real :: tmpc, tmpr, tmpi, tmps, tmpg
integer :: i1, i2, j1, j2, i, j, k, m, n

real(8) :: tmp_total, tmptot

#ifdef UWM_STATS
!weberjk(UWM), declare arrays to compute statistics. Make them less ambiguous
!than 'f' and 'df'
real t_before(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) 
real qt_before(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) 
real t_avg(nzm), t_before_avg(nzm) 
real qt_avg(nzm), qt_avg_before(nzm)
#endif /*UWM_STATS*/

call t_startf ('micro_proc')

if(mod(nstep-1,nstatis).eq.0.and.icycle.eq.1) then
   do j=1,ny
      do i=1,nx
         precsfc(i,j)=0.
      end do
   end do
   do k=1,nzm
      precflux(k) = 0.
   end do
end if

if(dostatis) then ! initialize arrays for statistics
   mfrac(:,:) = 0.
   mtend(:,:) = 0.
   trtau(:,:) = 0.
   qpfall(:)=0.
   tlat(:) = 0.
   tmtend3d(:,:,:) = 0.


! Microphysic process rates +++mhwang 
   mPRC=0.0  
   mPRA=0.0  
   mPSMLT=0.0  
   mEVPMS=0.0  
   mPRACS=0.0  
   mEVPMG=0.0  
   mPRACG=0.0  
   mPRE=0.0  
   mPGMLT=0.0   
   mMNUCCC=0.0  
   mPSACWS=0.0  
   mPSACWI=0.0  
   mQMULTS=0.0  
   mQMULTG=0.0  
   mPSACWG=0.0  
   mPGSACW=0.0  
   mPRD=0.0  
   mPRCI=0.0  
   mPRAI=0.0  
   mQMULTR=0.0  
   mQMULTRG=0.0  
   mMNUCCD=0.0  
   mPRACI=0.0  
   mPRACIS=0.0  
   mEPRD=0.0 
   mMNUCCR=0.0  
   mPIACR=0.0 
   mPIACRS=0.0  
   mPGRACS=0.0 
   mPRDS=0.0  
   mEPRDS=0.0  
   mPSACR=0.0 
   mPRDG=0.0  
   mEPRDG=0.0

   mNPRC1=0.0
   mNRAGG=0.0
   mNPRACG=0.0 
   mNSUBR=0.0 
   mNSMLTR=0.0 
   mNGMLTR=0.0 
   mNPRACS=0.0 
   mNNUCCR=0.0 
   mNIACR=0.0 
   mNIACRS=0.0 
   mNGRACS=0.0
   mNSMLTS=0.0 
   mNSAGG=0.0 
   mNPRCI=0.0 
   mNSCNG=0.0 
   mNSUBS=0.0
   
   mPCC=0.0
   mNNUCCC=0.0 
   mNPSACWS=0.0 
   mNPRA=0.0 
   mNPRC=0.0 
   mNPSACWI=0.0 
   mNPSACWG=0.0 
   mNPRAI=0.0
   mNMULTS=0.0 
   mNMULTG=0.0 
   mNMULTR=0.0 
   mNMULTRG=0.0 
   mNNUCCD=0.0
   mNSUBI=0.0 
   mNGMLTG=0.0 
   mNSUBG=0.0 
   mNACT=0.0 

   mSIZEFIX_NR=0.0 
   mSIZEFIX_NC=0.0
   mSIZEFIX_NI=0.0
   mSIZEFIX_NS=0.0
   mSIZEFIX_NG=0.0
   mNEGFIX_NR=0.0
   mNEGFIX_NC=0.0
   mNEGFIX_NI=0.0
   mNEGFIX_NS=0.0
   mNEGFIX_NG=0.0
   mNIM_MORR_CL = 0.0

   mQC_INST=0.0 
   mQR_INST=0.0 
   mQI_INST=0.0 
   mQS_INST=0.0 
   mQG_INST=0.0
   mNC_INST=0.0 
   mNR_INST=0.0 
   mNI_INST=0.0 
   mNS_INST=0.0 
   mNG_INST=0.0


#ifdef UWM_STATS
   !weberjk(UWM) Set 'before' values for budget statistics        
   do k=1,nzm 
     do j=dimy1_s,dimy2_s
       do i=dimx1_s,dimx2_s
         t_before(i,j,k) = t(i,j,k)
         qt_before(i,j,k) = micro_field(i,j,k,iqv)
       end do 
     end do 
   end do 
#endif /*UWM_STATS*/
end if
stend(:,:) = 0.
mksed(:,:) = 0.

!!$if(doprecip) total_water_prec = total_water_prec + total_water()
 
do j = 1,ny
   do i = 1,nx

      ! zero out mixing ratios of microphysical species
      tmpqv(:) = 0.
      tmpqcl(:) = 0.
      tmpncl(:) = 0.
      tmpqr(:) = 0.
      tmpnr(:) = 0.
      tmpqci(:) = 0.
      tmpnci(:) = 0.
      tmpqs(:) = 0.
      tmpns(:) = 0.
      tmpqg(:) = 0.
      tmpng(:) = 0.

      ! get microphysical quantities in this grid column
      tmpqv(:) = micro_field(i,j,:,iqv) !bloss/qt: This is total water (qv+qcl)

!bloss/qt: compute below from saturation adjustment.
!bloss/qt      tmpqcl(:) = micro_field(i,j,:,iqcl)
#ifndef UWM_STATS
      if(dopredictNc) tmpncl(:) = micro_field(i,j,:,incl)
#else
      if(dopredictNc) then
        tmpncl(:) = micro_field(i,j,:,incl)
      else
        tmpncl(:) = Nc0
      endif
#endif
      if(doprecip) then
         tmpqr(:) = micro_field(i,j,:,iqr)
         tmpnr(:) = micro_field(i,j,:,inr)
      end if

      if(doicemicro) then
         tmpqci(:) = micro_field(i,j,:,iqci)
         tmpnci(:) = micro_field(i,j,:,inci)
         tmpqs(:) = micro_field(i,j,:,iqs)
         tmpns(:) = micro_field(i,j,:,ins)
         if(dograupel) then
            tmpqg(:) = micro_field(i,j,:,iqg)
            tmpng(:) = micro_field(i,j,:,ing)
         end if
      end if

      ! get absolute temperature in this column
      !bloss/qt: before saturation adjustment for liquid,
      !          this is Tcl = T - (L/Cp)*qcl (the cloud liquid water temperature)
      tmptabs(:) = t(i,j,:)  &           ! liquid water-ice static energy over Cp
           - gamaz(:) &                                   ! potential energy
           + fac_cond * (tmpqr(:)) &    ! bloss/qt: liquid latent energy due to rain only
           + fac_sub  * (tmpqci(:) + tmpqs(:) + tmpqg(:)) ! ice latent energy

      tmpdz = adz(:)*dz
!      tmpw = 0.5*(w(i,j,1:nzm) + w(i,j,2:nz))  ! MK: changed for stretched grids 
      tmpw = ((zi(2:nz)-z(1:nzm))*w(i,j,1:nzm)+ &
             (z(1:nzm)-zi(1:nzm))*w(i,j,2:nz))/(zi(2:nz)-zi(1:nzm))
#ifdef CLUBB
      ! Added by dschanen on 4 Nov 2008 to account for w_sgs 
      if ( doclubb .and. dosubgridw ) then
        ! Compute w_sgs.  Formula is consistent with that used with 
        ! TKE from MYJ pbl scheme in WRF (see module_mp_graupel.f90).
        tmpwsub = sqrt( LIN_INT( real( wp2(i,j,2:nz) ), real( wp2(i,j,1:nzm) ), &
                                  zi(2:nz), zi(1:nzm), z(1:nzm) ) )
      else
        tmpwsub = 0.
      end if
#else /* Old code */
      tmpwsub = 0.
#endif
#ifdef CLUBB
      if ( doclubb ) then
        cloud_frac_in(1:nzm) = cloud_frac(i,j,2:nz)
      else
        cloud_frac_in(1:nzm) = 0.0
      end if
#endif

      tmppres(:) = 100.*pres(1:nzm)

      !bloss/qt: saturation adjustment to compute cloud liquid water content.
      !          Note: tmpqv holds qv+qcl on input, qv on output.
      !                tmptabs hold T-(L/Cp)*qcl on input, T on output.
      !                tmpqcl hold qcl on output.
      !                tmppres is unchanged on output, should be in Pa.
#ifdef CLUBB
      ! In the CLUBB case, we want to call the microphysics on sub-saturated grid
      ! boxes and weight by cloud fraction, therefore we use the CLUBB value of 
      ! liquid water. -dschanen 23 Nov 2009
      if ( .not. ( docloud .or. dosmoke ) ) then
        tmpqcl  = cloudliq(i,j,:) ! Liquid updated by CLUBB just prior to this
        tmpqv   = tmpqv - tmpqcl ! Vapor
        tmptabs = tmptabs + fac_cond * tmpqcl ! Update temperature
      else
        call satadj_liquid(nzm,tmptabs,tmpqv,tmpqcl,tmppres)
      end if
#else
      call satadj_liquid(nzm,tmptabs,tmpqv,tmpqcl,tmppres)
#endif

      i1 = 1 ! dummy variables used by WRF convention in subroutine call
      i2 = 1
      j1 = 1
      j2 = 1

      mtendqv = 0.
      mtendqcl = 0.
      mtendqr = 0.
      mtendqci = 0.
      mtendqs = 0.
      mtendqg = 0.
      mtendncl = 0.
      mtendnr = 0.
      mtendnci = 0.
      mtendns = 0.
      mtendng = 0.

      tmtend1d = 0.

      sfcpcp = 0.
      sfcicepcp = 0.

      effc1d(:) = 10. ! default liquid and ice effective radii
      effi1d(:) = 75.

! Zero out microphysical process rates  +++mhwang 
      PRC=0.0 
      PRA=0.0 
      PSMLT=0.0 
      EVPMS=0.0 
      PRACS=0.0 
      EVPMG=0.0 
      PRACG=0.0 
      PRE=0.0 
      PGMLT=0.0 
      MNUCCC=0.0 
      PSACWS=0.0 
      PSACWI=0.0 
      QMULTS=0.0 
      QMULTG=0.0 
      PSACWG=0.0 
      PGSACW=0.0 
      PRD=0.0 
      PRCI=0.0 
      PRAI=0.0 
      QMULTR=0.0 
      QMULTRG=0.0 
      MNUCCD=0.0 
      PRACI=0.0 
      PRACIS=0.0 
      EPRD=0.0 
      MNUCCR=0.0 
      PIACR=0.0 
      PIACRS=0.0 
      PGRACS=0.0 
      PRDS=0.0 
      EPRDS=0.0 
      PSACR=0.0 
      PRDG=0.0 
      EPRDG=0.0 

      NPRC1=0.0 
      NRAGG=0.0 
      NPRACG=0.0 
      NSUBR=0.0 
      NSMLTR=0.0 
      NGMLTR=0.0 
      NPRACS=0.0 
      NNUCCR=0.0 
      NIACR=0.0 
      NIACRS=0.0 
      NGRACS=0.0
      NSMLTS=0.0
      NSAGG=0.0 
      NPRCI=0.0 
      NSCNG=0.0 
      NSUBS=0.0

      PCC=0.0
      NNUCCC=0.0
      NPSACWS=0.0
      NPRA=0.0
      NPRC=0.0
      NPSACWI=0.0
      NPSACWG=0.0
      NPRAI=0.0
      NMULTS=0.0
      NMULTG=0.0
      NMULTR=0.0
      NMULTRG=0.0
      NNUCCD=0.0
      NSUBI=0.0
      NGMLTG=0.0
      NSUBG=0.0
      NACT=0.0

      SIZEFIX_NR=0.0
      SIZEFIX_NC=0.0
      SIZEFIX_NI=0.0
      SIZEFIX_NS=0.0
      SIZEFIX_NG=0.0
      NEGFIX_NR=0.0
      NEGFIX_NC=0.0
      NEGFIX_NI=0.0
      NEGFIX_NS=0.0
      NEGFIX_NG=0.0
      NIM_MORR_CL = 0.0

      QC_INST=0.0
      QR_INST=0.0
      QI_INST=0.0
      QS_INST=0.0
      QG_INST=0.0
      NC_INST=0.0
      NR_INST=0.0
      NI_INST=0.0
      NS_INST=0.0
      NG_INST=0.0


#ifdef CLUBB
      if ( any( tmpqv < 0. ) ) then
        qv_clip(2:nz) = tmpqv(1:nzm)
        qv_clip(1) = 0.0_core_rknd
        if ( clubb_at_least_debug_level( 1 ) ) then
          write(0,*) "M2005 has received a negative water vapor"
        end if
        call fill_holes_vertical( 2, 0._core_rknd, "zt", rho_ds_zt, rho_ds_zm, qv_clip )
        tmpqv = qv_clip(2:nz)
      end if
      if ( any( tmpqcl < 0. ) ) then
        qcl_clip(2:nz) = tmpqcl(1:nzm)
        qcl_clip(1) = 0.0_core_rknd
        if ( clubb_at_least_debug_level( 1 ) ) then
          write(0,*) "M2005 has received a negative liquid water"
        end if
        call fill_holes_vertical( 2, 0._core_rknd, "zt", rho_ds_zt, rho_ds_zm, qcl_clip )
        tmpqcl = qcl_clip(2:nz)
      end if

      ! Set autoconversion and accretion rates to 0;  these are diagnostics and
      ! don't feed back into the calculation.
      PRC = 0.
      PRA = 0.
      PRE = 0.
#endif /*CLUBB*/
      
      ! explanation of variable names:
      !   mtend1d: array of 1d profiles of microphysical tendencies (w/o sed.)
      !   stend1d: array of 1d profiles of sedimentation tendencies for q*
      !   tmp**: on input, current value of **.  On output, new value of **.
      !   eff*1d: one-dim. profile of effective raduis for *
      call m2005micro_graupel(&
           mtendqcl,mtendqci,mtendqs,mtendqr, &
           mtendncl,mtendnci,mtendns,mtendnr, &
           tmpqcl,tmpqci,tmpqs,tmpqr, &
           tmpncl,tmpnci,tmpns,tmpnr, &
           tmtend1d,mtendqv, &
           tmptabs,tmpqv,tmppres,rho,tmpdz,tmpw,tmpwsub, &
#ifdef CLUBB           
           rain_vel,&
#endif /*CLUBB*/
           sfcpcp, sfcicepcp, &
           effc1d,effi1d,effs1d,effr1d, &
           dtn, &
           i1,i2, j1,j2, 1,nzm, i1,i2, j1,j2, 1,nzm, &
           mtendqg,mtendng,tmpqg,tmpng,effg1d,stendqg, &
           stendqr,stendqci,stendqs,stendqcl, &
           stendng, stendnr, stendnci, stendns, stendncl, &
#ifdef CLUBB
           cloud_frac_in, & ! cloud_frac added by dschanen UWM
#endif /*CLUBB*/
           ! These variables were renamed during the merging for
           ! clubb:ticket:552
           ! qr_auto = PRC, qr_accr = PRA, qr_evap = PRE
           PRC,  PRA,  &
           PSMLT, EVPMS, PRACS, EVPMG, PRACG, PRE, PGMLT,  &
           MNUCCC, PSACWS, PSACWI, QMULTS, QMULTG, PSACWG, PGSACW, &
           PRD, PRCI, PRAI, QMULTR, QMULTRG, MNUCCD, PRACI, PRACIS, EPRD, &
           MNUCCR, PIACR, PIACRS, PGRACS, PRDS, EPRDS, PSACR, &
           PRDG, EPRDG, NPRC1, NRAGG, NPRACG, NSUBR, NSMLTR, NGMLTR, &
           NPRACS, NNUCCR, NIACR, NIACRS, NGRACS, NSMLTS, NSAGG, NPRCI, NSCNG, NSUBS, &
           PCC, NNUCCC, NPSACWS, NPRA, NPRC, NPSACWI, NPSACWG, NPRAI, &
           NMULTS, NMULTG, NMULTR, NMULTRG, NNUCCD, NSUBI, NGMLTG, NSUBG, NACT, &
           SIZEFIX_NR, SIZEFIX_NC, SIZEFIX_NI, SIZEFIX_NS, SIZEFIX_NG, &
           NEGFIX_NR, NEGFIX_NC, NEGFIX_NI, NEGFIX_NS, NEGFIX_NG, &
           NIM_MORR_CL, QC_INST, QR_INST, QI_INST, QS_INST, QG_INST, &
           NC_INST, NR_INST, NI_INST, NS_INST, NG_INST   )
#ifdef CLUBB
      if ( any( tmpqv < 0. ) ) then
        qv_clip(2:nz) = tmpqv(1:nzm)
        qv_clip(1) = 0.0_core_rknd
        if ( clubb_at_least_debug_level( 1 ) ) then
          write(0,*) "M2005 has produced a negative water vapor"
        end if
        call fill_holes_vertical( 2, 0._core_rknd, "zt", rho_ds_zt, rho_ds_zm, qv_clip )
        tmpqv = qv_clip(2:nz)
      end if
      if ( any( tmpqcl < 0. ) ) then
        qcl_clip(2:nz) = tmpqcl(1:nzm)
        qcl_clip(1) = 0.0_core_rknd
        if ( clubb_at_least_debug_level( 1 ) ) then
          write(0,*) "M2005 has produced a negative liquid water"
        end if
        call fill_holes_vertical( 2, 0._core_rknd, "zt", rho_ds_zt, rho_ds_zm, qcl_clip )
        tmpqcl = qcl_clip(2:nz)
      end if
#endif /*CLUBB*/

#ifdef UWM_STATS
! We want the tendencies to be explicity: (before - after) / dt
! Therefore, we overwrite the output from Morrison
if(dostatis) then

  mtendqcl(:) = ( tmpqcl(:) - cloudliq(i,j,:) ) / dtn
  mtendqv(:) = ( tmpqv(:) - ( micro_field(i,j,:,iqv) - cloudliq(i,j,:) ) ) / dtn

if(doprecip .EQV. .FALSE.) then
  tmpqcl(:) = tmpqcl(:) + tmpqr(:) ! add rain mass back to cloud water
  tmpncl(:) = tmpncl(:) + tmpnr(:) ! add rain number back to cloud water

  ! zero out rain 
  tmpqr(:) = 0.
  tmpnr(:) = 0.
end if !doprecip is false

  if(dopredictNc) then
    mtendncl(:) = (tmpncl(:) - micro_field(i,j,:,incl) ) / dtn
  end if

  if(doprecip) then
    mtendnr(:) = (tmpnr(:) - micro_field(i,j,:,inr) )/ dtn
    mtendqr(:) = (tmpqr(:) - micro_field(i,j,:,iqr) )/ dtn
  end if

  if(doicemicro) then
    mtendnci(:) = (tmpnci(:) - micro_field(i,j,:,inci) ) / dtn
    mtendqci(:) = (tmpqci(:) - micro_field(i,j,:,iqci) ) / dtn
    
    mtendns(:) = (tmpns(:) - micro_field(i,j,:,ins) ) / dtn
    mtendqs(:) = (tmpqs(:) - micro_field(i,j,:,iqs) ) / dtn
    
    if(dograupel) then
      mtendng(:) = (tmpng(:) - micro_field(i,j,:,ing) ) / dtn
      mtendqg(:) = (tmpqg(:) - micro_field(i,j,:,iqg) ) / dtn
    endif

  endif
endif
#endif /*UWM_STATS*/
     ! update microphysical quantities in this grid column
      if(doprecip) then
         total_water_prec = total_water_prec + sfcpcp

         ! take care of surface precipitation
         precsfc(i,j) = precsfc(i,j) + sfcpcp/dz
         prec_xy(i,j) = prec_xy(i,j) + sfcpcp/dz

         ! update rain
         micro_field(i,j,:,iqr) = tmpqr(:)
         micro_field(i,j,:,inr) = tmpnr(:)
      else
         ! add rain to cloud
         tmpqcl(:) = tmpqcl(:) + tmpqr(:) ! add rain mass back to cloud water
         tmpncl(:) = tmpncl(:) + tmpnr(:) ! add rain number back to cloud water

         ! zero out rain 
         tmpqr(:) = 0.
         tmpnr(:) = 0.

         ! add rain tendencies to cloud
         stendqcl(:) = stendqcl(:) + stendqr(:)
         stendncl(:) = stendncl(:) + stendnr(:)
         mtendqcl(:) = mtendqcl(:) + mtendqr(:)
         mtendncl(:) = mtendncl(:) + mtendnr(:)

         ! zero out rain tendencies
         stendqr(:) = 0.
         stendnr(:) = 0.
         mtendqr(:) = 0.
         mtendnr(:) = 0.
      end if

      !bloss/qt: update total water and cloud liquid.
      !          Note: update of total water moved to after if(doprecip),
      !                  since no precip moves rain --> cloud liq.
      micro_field(i,j,:,iqv) = tmpqv(:) + tmpqcl(:) !bloss/qt: total water
      cloudliq(i,j,:) = tmpqcl(:) !bloss/qt: auxilliary cloud liquid water variable
#ifndef UWM_STATS
      if(dopredictNc) micro_field(i,j,:,incl) = tmpncl(:)
#else
      ! This will store the value of tmpncl after the call to Morrison
      ! If Nc0 is prescribed, it should be very close to the prescribed 
      ! Nc0.
      micro_field(i,j,:,incl) = tmpncl(:)
#endif
      reffc(i,j,:) = effc1d(:)

      if(doicemicro) then
         micro_field(i,j,:,iqci) = tmpqci(:)
         micro_field(i,j,:,inci) = tmpnci(:)
         micro_field(i,j,:,iqs) = tmpqs(:)
         micro_field(i,j,:,ins) = tmpns(:)
         if(dograupel) then
            micro_field(i,j,:,iqg) = tmpqg(:)
            micro_field(i,j,:,ing) = tmpng(:)
         end if
         reffi(i,j,:) = effi1d(:)  
      end if

      !=====================================================
      ! update liquid-ice static energy due to precipitation
      t(i,j,:) = t(i,j,:) &
           - dtn*fac_cond*(stendqcl+stendqr) &
           - dtn*fac_sub*(stendqci+stendqs+stendqg)
      !=====================================================
      if(dostatis) then
!bloss/qt: total water microphysical tendency includes qv and qcl
         mtend(:,iqv) = mtend(:,iqv) + mtendqv + mtendqcl

!bloss/qt         mtend(:,iqcl) = mtend(:,iqcl) + mtendqcl
         if(dopredictNc) mtend(:,incl) = mtend(:,incl) + mtendncl
         if(doprecip) then
            mtend(:,iqr) = mtend(:,iqr) + mtendqr
            mtend(:,inr) = mtend(:,inr) + mtendnr
         end if

         if(doicemicro) then
            mtend(:,iqci) = mtend(:,iqci) + mtendqci
            mtend(:,inci) = mtend(:,inci) + mtendnci
            !bloss            stend(:,inci) = stend(:,inci) + stendnci

            mtend(:,iqs) = mtend(:,iqs) + mtendqs
            mtend(:,ins) = mtend(:,ins) + mtendns
            !bloss            stend(:,ins) = stend(:,ins) + stendns

            if(dograupel) then
               mtend(:,iqg) = mtend(:,iqg) + mtendqg
               mtend(:,ing) = mtend(:,ing) + mtendng
               !bloss            stend(:,ing) = stend(:,ing) + stendng
            end if
         end if

         do n = 1,nmicro_fields
            do k = 1,nzm
               if(micro_field(i,j,k,n).ge.1.e-6) mfrac(k,n) = mfrac(k,n)+1.
            end do
         end do

         ! approximate optical depth = 0.0018*lwp/effrad
         !  integrated up to level at which output
         tmpc = 0.
         tmpr = 0.
         tmpi = 0.
         tmps = 0.
         tmpg = 0.

         do k = 1,nzm
            tmpc = tmpc + 0.0018*rho(k)*dz*adz(k)*tmpqcl(k)/(1.e-20+1.e-6*effc1d(k))
            tmpr = tmpr + 0.0018*rho(k)*dz*adz(k)*tmpqr(k)/(1.e-20+1.e-6*effr1d(k))
            !bloss/qt: put cloud liquid optical depth in trtau(:,iqv)
            trtau(k,iqv) = trtau(k,iqv) + tmpc
            if(doprecip) trtau(k,iqr) = trtau(k,iqr) + tmpr

            if(doicemicro) then
               tmpi = tmpi + 0.0018*rho(k)*dz*adz(k)*tmpqci(k)/(1.e-20+1.e-6*effi1d(k))
               tmps = tmps + 0.0018*rho(k)*dz*adz(k)*tmpqs(k)/(1.e-20+1.e-6*effs1d(k))
               tmpg = tmpg + 0.0018*rho(k)*dz*adz(k)*tmpqg(k)/(1.e-20+1.e-6*effg1d(k))

               trtau(k,iqci) = trtau(k,iqci) + tmpi
               trtau(k,iqs) = trtau(k,iqs) + tmps
#ifdef CLUBB /* Bug fix -dschanen 9 Mar 2012 */
               if ( dograupel ) then
                 trtau(k,iqg) = trtau(k,iqg) + tmpg
               end if
#else
               trtau(k,iqg) = trtau(k,iqg) + tmpg
#endif /* CLUBB */
            end if
         end do

         tlat(1:nzm) = tlat(1:nzm) &
              - dtn*fac_cond*(stendqcl+stendqr) &
              - dtn*fac_sub*(stendqci+stendqs+stendqg)
         qpfall(1:nzm) = qpfall(1:nzm) + dtn*(stendqr+stendqs+stendqg)

         !bloss: temperature tendency (sensible heating) due to phase changes
         tmtend3d(i,j,1:nzm) = tmtend1d(1:nzm)

! Microphysical process rates +++mhwang 
         mPRC=mPRC + PRC 
         mPRA=mPRA + PRA 
         mPSMLT=mPSMLT + PSMLT 
         mEVPMS=mEVPMS + EVPMS 
         mPRACS=mPRACS + PRACS 
         mEVPMG=mEVPMG + EVPMG 
         mPRACG=mPRACG + PRACG 
         mPRE=mPRE + PRE 
         mPGMLT=mPGMLT + PGMLT 
         mMNUCCC=mMNUCCC + MNUCCC 
         mPSACWS=mPSACWS + PSACWS 
         mPSACWI=mPSACWI + PSACWI 
         mQMULTS=mQMULTS + QMULTS 
         mQMULTG=mQMULTG + QMULTG 
         mPSACWG=mPSACWG + PSACWG 
         mPGSACW=mPGSACW + PGSACW 
         mPRD=mPRD + PRD 
         mPRCI=mPRCI + PRCI 
         mPRAI=mPRAI + PRAI 
         mQMULTR=mQMULTR + QMULTR 
         mQMULTRG=mQMULTRG + QMULTRG 
         mMNUCCD=mMNUCCD + MNUCCD 
         mPRACI=mPRACI + PRACI 
         mPRACIS=mPRACIS + PRACIS 
         mEPRD=mEPRD + EPRD  
         mMNUCCR=mMNUCCR + MNUCCR 
         mPIACR=mPIACR + PIACR 
         mPIACRS=mPIACRS + PIACRS 
         mPGRACS=mPGRACS + PGRACS 
         mPRDS=mPRDS + PRDS 
         mEPRDS=mEPRDS + EPRDS 
         mPSACR=mPSACR + PSACR 
         mPRDG=mPRDG + PRDG 
         mEPRDG=mEPRDG + EPRDG 

         mNPRC1=mNPRC1+NPRC1
         mNRAGG=mNRAGG+NRAGG
         mNPRACG=mNPRACG+NPRACG
         mNSUBR=mNSUBR+NSUBR
         mNSMLTR=mNSMLTR+NSMLTR
         mNGMLTR=mNGMLTR+NGMLTR
           
         mNPRACS=mPRACS+NPRACS
         mNNUCCR=mNNUCCR+NNUCCR
         mNIACR=mNIACR+NIACR
         mNIACRS=mNIACRS+NIACRS
         mNGRACS=mNGRACS+NGRACS
         mNSMLTS=mNSMLTS+NSMLTS
         mNSAGG=mNSAGG+NSAGG
         mNPRCI=mNPRCI+NPRCI
         mNSCNG=mNSCNG+NSCNG
         mNSUBS=mNSUBS+NSUBS
         mPCC=mPCC+PCC
         mNNUCCC=mNNUCCC+NNUCCC
         mNPSACWS=mNPSACWS+NPSACWS
         mNPRA=mNPRA+NPRA
         mNPRC=mNPRC+NPRC
         mNPSACWI=mNPSACWI+NPSACWI
         mNPSACWG=mNPSACWG+NPSACWG
         mNPRAI=mNPRAI+NPRAI
         mNMULTS=mNMULTS+NMULTS
         mNMULTG=mNMULTG+NMULTG
         mNMULTR=mNMULTR+NMULTR
         mNMULTRG=mNMULTRG+NMULTRG
         mNNUCCD=mNNUCCD+NNUCCD
         mNSUBI=mNSUBI+NSUBI
         mNGMLTG=mNGMLTG+NGMLTG
         mNSUBG=mNSUBG+NSUBG
         mNACT=mNACT+NACT
         mSIZEFIX_NR=mSIZEFIX_NR+SIZEFIX_NR
         mSIZEFIX_NC=mSIZEFIX_NC+SIZEFIX_NC
         mSIZEFIX_NI=mSIZEFIX_NI+SIZEFIX_NI
         mSIZEFIX_NS=mSIZEFIX_NS+SIZEFIX_NS
         mSIZEFIX_NG=mSIZEFIX_NG+SIZEFIX_NG
         mNEGFIX_NR=mNEGFIX_NR+NEGFIX_NR
         mNEGFIX_NC=mNEGFIX_NC+NEGFIX_NC
         mNEGFIX_NI=mNEGFIX_NI+NEGFIX_NI
         mNEGFIX_NS=mNEGFIX_NS+NEGFIX_NS
         mNEGFIX_NG=mNEGFIX_NG+NEGFIX_NG
         mNIM_MORR_CL=mNIM_MORR_CL+NIM_MORR_CL    

         mQC_INST=mQC_INST+QC_INST
         mQR_INST=mQR_INST+QR_INST
         mQI_INST=mQI_INST+QI_INST
         mQS_INST=mQS_INST+QS_INST
         mQG_INST=mQG_INST+QG_INST
         mNC_INST=mNC_INST+NC_INST
         mNR_INST=mNR_INST+NR_INST
         mNI_INST=mNI_INST+NI_INST
         mNS_INST=mNS_INST+NS_INST
         mNG_INST=mNG_INST+NG_INST
         
    
endif !do statis

      stend(:,iqv) = stend(:,iqv) + stendqcl !bloss/qt: iqcl --> iqv
      if(dopredictNc) stend(:,incl) = stend(:,incl) + stendncl
      if(doprecip) then
         stend(:,iqr) = stend(:,iqr) + stendqr
         stend(:,inr) = stend(:,inr) + stendnr
      end if

      if(doicemicro) then
         stend(:,iqci) = stend(:,iqci) + stendqci
         stend(:,iqs) = stend(:,iqs) + stendqs
         stend(:,inci) = stend(:,inci) + stendnci
         stend(:,ins) = stend(:,ins) + stendns
         if(dograupel) stend(:,iqg) = stend(:,iqg) + stendqg
         if(dograupel) stend(:,ing) = stend(:,ing) + stendng
      end if

   end do ! i = 1,nx
end do ! j = 1,ny

#ifdef UWM_STATS      
if(dostatis) then
        !------------------------------------------------------------------------------------------
        !weberjk(UWM), compute the microphysical effects on t2, q2, tw, and qw
        !budgets
        call stat_varscalar(t,t_before,t_avg,t_before_avg,t2leprec)
        call stat_varscalar(micro_field(:,:,:,iqv),qt_before,qt_avg,qt_avg_before,q2leprec)
         
        call setvalue(twleprec,nzm,0.)
        call stat_sw2(t,t_before,twleprec)
         
        call setvalue(qwleprec,nzm,0.)
        call stat_sw2(micro_field(:,:,:,iqv),qt_before,qwleprec)

endif !do statis
#endif /*UWM_STATS*/

! back sedimentation flux out from sedimentation tendencies
tmpc = 0.
do k = 1,nzm
   m = nz-k
   tmpc = tmpc + stend(m,iqv)*rho(m)*dz*adz(m)  !bloss/qt: iqcl --> iqv
   mksed(m,iqv) = tmpc
end do
precflux(1:nzm) = precflux(1:nzm) - mksed(:,iqv)*dtn/dz

if(doprecip) then
   tmpr = 0.
   do k = 1,nzm
      m = nz-k
      tmpr = tmpr + stend(m,iqr)*rho(m)*dz*adz(m)
      mksed(m,iqr) = tmpr
   end do
   precflux(1:nzm) = precflux(1:nzm) - mksed(:,iqr)*dtn/dz
end if

if(doicemicro) then
   tmpi = 0.
   tmps = 0.
   tmpg = 0.
   do k = 1,nzm
      m = nz-k
      tmpi = tmpi + stend(m,iqci)*rho(m)*dz*adz(m)
      tmps = tmps + stend(m,iqs)*rho(m)*dz*adz(m)
#ifdef CLUBB /* Bug fix -dschanen 9 Mar 2012 */
      if ( dograupel ) then
        tmpg = tmpg + stend(m,iqg)*rho(m)*dz*adz(m)
      else
        tmpg = 0.
      end if
#else
      tmpg = tmpg + stend(m,iqg)*rho(m)*dz*adz(m)
#endif
      mksed(m,iqci) = tmpi
      mksed(m,iqs) = tmps
#ifdef CLUBB /* Bug fix -dschanen 9 Mar 2012 */
      if ( dograupel ) then
        mksed(m,iqg) = tmpg
      end if
#else
      mksed(m,iqg) = tmpg
#endif
   end do
#ifdef CLUBB /* Bug fix -dschanen 9 Mar 2012 */
   if ( dograupel ) then
     precflux(1:nzm) = precflux(1:nzm) &
          - (mksed(:,iqci) + mksed(:,iqs) + mksed(:,iqg))*dtn/dz
   else
     precflux(1:nzm) = precflux(1:nzm) &
          - (mksed(:,iqci) + mksed(:,iqs))*dtn/dz
   end if
#else
   precflux(1:nzm) = precflux(1:nzm) &
        - (mksed(:,iqci) + mksed(:,iqs) + mksed(:,iqg))*dtn/dz
#endif
end if

!!$if(doprecip) total_water_prec = total_water_prec - total_water()

#ifdef CLUBB
if (docloud.or.doclubb)  call micro_diagnose()   ! leave this line here
#else
if (docloud)  call micro_diagnose()   ! leave this line here
#endif

call t_stopf ('micro_proc')

end subroutine micro_proc

!----------------------------------------------------------------------
!!! Diagnose arrays nessesary for dynamical core and radiation:
!
!  This is the pace where the microphysics field that SAM actually cares about
!  are diagnosed.

subroutine micro_diagnose()

use vars
#ifdef CLUBB
use error_code, only: clubb_at_least_debug_level ! Procedure
use constants_clubb, only: fstderr, zero_threshold
implicit none
#endif

real omn, omp
integer i,j,k

! water vapor = total water - cloud liquid
qv(1:nx,1:ny,1:nzm) = micro_field(1:nx,1:ny,1:nzm,iqv) &
     - cloudliq(1:nx,1:ny,1:nzm)

#ifdef CLUBB
do i = 1, nx
  do j = 1, ny
    do k = 1, nzm
      ! Apply local hole-filling to vapor by converting liquid to vapor. Moist
      ! static energy should be conserved, so updating temperature is not
      ! needed here. -dschanen 31 August 2011
      if ( qv(i,j,k) < zero_threshold ) then
        cloudliq(i,j,k) = cloudliq(i,j,k) + qv(i,j,k)
        qv(i,j,k) = zero_threshold
        if ( cloudliq(i,j,k) < zero_threshold ) then
          if ( clubb_at_least_debug_level( 1 ) ) then
            write(fstderr,*) "Total water at", "i =", i, "j =", j, "k =", k, "is negative.", &
              "Applying non-conservative hard clipping."
          end if
          cloudliq(i,j,k) = zero_threshold
        end if ! cloud_liq < 0
      end if ! qv < 0
    end do ! 1.. nzm
  end do ! 1.. ny
end do ! 1.. nx
#endif /* CLUBB */
! cloud liquid water
qcl(1:nx,1:ny,1:nzm) = cloudliq(1:nx,1:ny,1:nzm)

! rain water
if(doprecip) qpl(1:nx,1:ny,1:nzm) = micro_field(1:nx,1:ny,1:nzm,iqr)

! cloud ice 
if(doicemicro) then
   qci(1:nx,1:ny,1:nzm) = micro_field(1:nx,1:ny,1:nzm,iqci)

   if(dograupel) then
      qpi(1:nx,1:ny,1:nzm) = micro_field(1:nx,1:ny,1:nzm,iqs) &
           + micro_field(1:nx,1:ny,1:nzm,iqg)
   else
      qpi(1:nx,1:ny,1:nzm) = micro_field(1:nx,1:ny,1:nzm,iqs)
   end if
end if

end subroutine micro_diagnose

#ifdef CLUBB
!---------------------------------------------------------------------
subroutine micro_update()

! Description:
! This subroutine essentially does what micro_proc does but does not
! call any microphysics subroutines.  We need to do this for the 
! single-moment bulk microphysics (SAM1MOM) so that CLUBB gets a
! properly updated value of ice fed in. 
!
! -dschanen UWM
!---------------------------------------------------------------------

  ! Update the dynamical core variables (e.g. qv, qcl) with the value in
  ! micro_field.  Diffusion, advection, and other processes are applied to
  ! micro_field but not the variables in vars.f90
  call micro_diagnose()

  return
end subroutine micro_update

!---------------------------------------------------------------------
subroutine micro_adjust( new_qv, new_qc )
! Description:
!   Adjust total water in SAM based on values from CLUBB.
! References:
!   None
!---------------------------------------------------------------------

  use vars, only: qci

  implicit none

  real, dimension(nx,ny,nzm), intent(in) :: &
    new_qv, & ! Water vapor mixing ratio that has been adjusted by CLUBB [kg/kg]
    new_qc    ! Cloud water mixing ratio that has been adjusted by CLUBB [kg/kg]

  ! Total water mixing ratio
  micro_field(1:nx,1:ny,1:nzm,iqv) = new_qv(1:nx,1:ny,1:nzm) &
                                   + new_qc(1:nx,1:ny,1:nzm)

  ! Cloud water mixing ratio
  cloudliq(1:nx,1:ny,1:nzm) = new_qc(1:nx,1:ny,1:nzm) 

  return
end subroutine micro_adjust

#endif /*CLUBB*/
!----------------------------------------------------------------------
!!! functions to compute terminal velocity for precipitating variables:
!
! you need supply functions to compute terminal velocity for all of your 
! precipitating prognostic variables. Note that all functions should
! compute vertical velocity given two microphysics parameters var1, var2, 
! and temperature, and water vapor (single values, not arrays). Var1 and var2 
! are some microphysics variables like water content and concentration.
! Don't change the number of arguments or their meaning!

!!$real function term_vel_qr(qr,nr,tabs,rho)
!!$! .......  
!!$end function term_vel_qr
!!$
!!$real function term_vel_Nr(qr,nr,tabs,rho)
!!$! .......  
!!$end function term_vel_Nr
!!$
!!$real function term_vel_qs(qs,ns,tabs,rho)
!!$! .......  
!!$end function term_vel_qs

! etc.

!----------------------------------------------------------------------
!!! compute sedimentation 
!
!  The perpose of this subroutine is to prepare variables needed to call
! the precip_all() for each of the falling hydrometeor varibles
subroutine micro_precip_fall()

! before calling precip_fall() for each of falling prognostic variables,
! you need to set hydro_type and omega(:,:,:) variables.
! hydro_type can have four values:
! 0 - variable is liquid water mixing ratio
! 1 - hydrometeor is ice mixing ratio
! 2 - hydrometeor is mixture-of-liquid-and-ice mixing ratio. (As in original SAM microphysics).
! 3 - variable is not mixing ratio, but, for example, rain drop concentration
! OMEGA(:,:,:) is used only for hydro_type=2, and is the fraction of liquid phase (0-1).
! for hour hypothetical case, there is no mixed hydrometeor, so omega is not actually used.

integer hydro_type
real omega(nx,ny,nzm) 

integer i,j,k

return ! do not need this routine -- sedimentation done in m2005micro.

!!$! Initialize arrays that accumulate surface precipitation flux
!!$
!!$ if(mod(nstep-1,nstatis).eq.0.and.icycle.eq.1) then
!!$   do j=1,ny
!!$    do i=1,nx
!!$     precsfc(i,j)=0.
!!$    end do
!!$   end do
!!$   do k=1,nzm
!!$    precflux(k) = 0.
!!$   end do
!!$ end if
!!$
!!$ do k = 1,nzm ! Initialize arrays which hold precipitation fluxes for stats.
!!$    qpfall(k)=0.
!!$    tlat(k) = 0.
!!$ end do
!!$   
!!$! Compute sedimentation of falling variables:
!!$
!!$ hydro_type=0
!!$ call precip_fall(qr, term_vel_qr, hydro_type, omega)
!!$ hydro_type=3
!!$ call precip_fall(Nr, term_vel_Nr, hydro_type, omega)
!!$ hydro_type=1
!!$ call precip_fall(qs, term_vel_qs, hydro_type, omega)
!!$ hydro_type=3
!!$ call precip_fall(Ns, term_vel_Ns, hydro_type, omega)
!!$ hydro_type=1
!!$ call precip_fall(qg, term_vel_qg, hydro_type, omega)
!!$ hydro_type=3
!!$ call precip_fall(Ng, term_vel_Ng, hydro_type, omega)
!!$


end subroutine micro_precip_fall

!----------------------------------------------------------------------
! called when stepout() called

subroutine micro_print()
  implicit none
  integer :: k

  ! print out min/max values of all microphysical variables
  do k=1,nmicro_fields
     call fminmax_print(trim(mkname(k))//':', &
          micro_field(:,:,:,k),dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm)
  end do

end subroutine micro_print

!----------------------------------------------------------------------
!!! Initialize the list of microphysics statistics that will be outputted
!!  to *.stat statistics file

subroutine micro_hbuf_init(namelist,deflist,unitlist,status,average_type,count,microcount)


#ifdef UWM_STATS
use compute_correlation_module, only: &
    corr_avg
#endif /*UWM_STATS*/
character(*) namelist(*), deflist(*), unitlist(*)
integer status(*),average_type(*),count,microcount, n, ii, jj, ncond

#ifdef UWM_STATS
character*20 name
#else
character*10 name
#endif /*UWM_STATS*/
character*80 longname
character*10 units

microcount = 0

name = 'QTFLUX'
longname = 'Total (resolved + subgrid) total water (vapor+cloud) flux'
units = 'W/m2'
call add_to_namelist(count,microcount,name,longname,units,0)

do n = 1,nmicro_fields
!bloss/qt   if(n.ne.iqv) then
  ! add mean value of microphysical field to statistics
  !   EXCEPT for water vapor (added in statistics.f90)
  name = trim(mkname(n))
  longname = trim(mklongname(n))
  units = trim(mkunits(n))
  call add_to_namelist(count,microcount,name,longname,units,0)
  if(n.eq.iqv) then
      ! add variance of ONLY total water (vapor + cloud liq) field to statistics
      !   cloud water variance and cloud ice variance
      !   already output in statistics.f90
      name = trim(mkname(n))//'2'
      longname = 'Variance of '//trim(mklongname(n))
      units = '('//trim(mkunits(n))//')^2'
      call add_to_namelist(count,microcount,name,longname,units,0)
   end if

   ! add vertical advective tendency
   name = trim(mkname(n))//'ADV'
   longname = 'Tendency of '//trim(mklongname(n))// &
        ' due to resolved vertical advection'
   units = trim(mkunits(n))//'/day'
   call add_to_namelist(count,microcount,name,longname,units,0)

   ! add vertical diffusive tendency
   name = trim(mkname(n))//'DIFF'
   longname = 'Tendency of '//trim(mklongname(n))// &
        ' due to vertical SGS transport'
   units = trim(mkunits(n))//'/day'
   call add_to_namelist(count,microcount,name,longname,units,0)

   ! add tendency due to large-scale vertical advection
   name = trim(mkname(n))//'LSADV'
   longname = 'Tendency of '//trim(mklongname(n))// &
        ' due to large-scale vertical advection'
   units = trim(mkunits(n))//'/day'
   call add_to_namelist(count,microcount,name,longname,units,0)

   ! add tendency due to microphysical processes
   name = trim(mkname(n))//'MPHY'
   longname = 'Tendency of '//trim(mklongname(n))// &
        ' due to microphysical processes'
   units = trim(mkunits(n))//'/day'
   call add_to_namelist(count,microcount,name,longname,units,0)

   ! add vertical diffusive tendency
   name = trim(mkname(n))//'SED'
   longname = 'Tendency of '//trim(mklongname(n))//' due to sedimentation'
   units = trim(mkunits(n))//'/day'
   call add_to_namelist(count,microcount,name,longname,units,0)

   ! add storage terms of microphysical variables 
   name = trim(mkname(n))//'STO'
   longname = 'Storage term of '//trim(mklongname(n))
   units = trim(mkunits(n))//'/day'
   call add_to_namelist(count,microcount,name,longname,units,0)

   if(flag_wmass(n).gt.0) then
      ! fluxes output in W/m2 for mass mixing ratios
      units = 'W/m2'
   else
      ! fluxes output in #/m2/s for number concentrations
      units = '#/m2/s'
   end if

   ! add flux of microphysical fields to scalar
   name = trim(mkname(n))//'FLX'
   longname = 'Total flux (Resolved + Subgrid) of '//trim(mklongname(n))
   call add_to_namelist(count,microcount,name,longname,units,0)

   ! add subgrid flux of microphysical fields to scalar
   name = trim(mkname(n))//'FLXS'
   longname = 'Subgrid flux of '//trim(mklongname(n))
   call add_to_namelist(count,microcount,name,longname,units,0)

   ! add sedimentation flux of microphysical fields to scalar
   name = trim(mkname(n))//'SDFLX'
   longname = 'Sedimentation flux of '//trim(mklongname(n))
   call add_to_namelist(count,microcount,name,longname,units,0)

   if((flag_wmass(n).gt.0).and.(n.ne.iqv)) then
      ! add area fraction of microphysical field to statistics
      name = trim(mkname(n))//'FRAC'
      longname = trim(mklongname(n))//' FRACTION'
      units = '1'
      call add_to_namelist(count,microcount,name,longname,units,0)

      ! add approximate optical depth of hydrometeor fields
      name = 'TAU'//trim(mkname(n))
      longname = 'Approx optical depth of '//trim(mklongname(n))
      units = '1'
      call add_to_namelist(count,microcount,name,longname,units,0)

!bloss (Apr 09): Eliminate this output.  It is unreliable when
!            hydrometeor fractions are variable across processors
!            or in time.  You can still compute this from 
!            TAU* and Q* values in the statistics file.
!bloss      ! add approximate optical depth of hydrometeor fields
!bloss      name = 'EFFR'//trim(mkname(n))
!bloss      longname = 'Effective radius of '//trim(mklongname(n))
!bloss      units = 'microns'
!bloss      call add_to_namelist(count,microcount,name,longname,units,0)

      ! add field which can be used to recover mean effective radius.
      name = trim(mkname(n))//'OEFFR'
      longname = 'Mixing ratio of '//trim(mklongname(n)) &
           //' over effective radius, EFFR = ' &
           //trim(mkname(n))//'/'//trim(mkname(n))//'OEFFR'
      units = 'g/kg/microns'
      call add_to_namelist(count,microcount,name,longname,units,0)
   end if

end do

!bloss/qt: add output for cloud liquid water (not included explicitly in 
!  total water formulation).
call add_to_namelist(count,microcount,'QC', &
     'Cloud liquid water mass mixing ratio', 'g/kg',0)

! add approximate optical depth of cloud liquid water
name = 'TAUQC'
longname = 'Approx optical depth of cloud liquid water'
units = '1'
call add_to_namelist(count,microcount,name,longname,units,0)

! add field which can be used to recover mean effective radius.
name = 'QCOEFFR'
longname = 'Mixing ratio of QC'// &
     ' over effective radius, EFFR = QC/QCOEFFR'
units = 'g/kg/microns'
call add_to_namelist(count,microcount,name,longname,units,0)

!bloss/qt: Can't be computed reliably in total water formulation.
! add temperature tendency (sensible energy) tendency due to mphys
!bloss call add_to_namelist(count,microcount,'QLAT', &
!bloss      'Sensible energy tendency due to phase changes', 'K/day',0)
 
do ncond = 1,ncondavg
   ! add conditional averages of hydrometeor fields
!bloss/qt: Can't be computed reliably in total water formulation.
!bloss   call add_to_namelist(count,microcount,'QLAT' // TRIM(condavgname(ncond)), &
!bloss        'Sensible energy tendency due to phase changes in ' // TRIM(condavglongname(ncond)), &
!bloss        'K/day',ncond)
   !bloss/qt: add conditional averages for water vapor and cloud liquid water
   call add_to_namelist(count,microcount,'QV' // TRIM(condavgname(ncond)), &
        'Water vapor mixing ratio in ' // TRIM(condavglongname(ncond)),'kg/kg',ncond)
   call add_to_namelist(count,microcount,'QC' // TRIM(condavgname(ncond)), &
        'Cloud liquid water mixing ratio in ' // TRIM(condavglongname(ncond)),'kg/kg',ncond)
   do n = 1,nmicro_fields
      call add_to_namelist(count,microcount,trim(mkname(n)) // TRIM(condavgname(ncond)), &
           trim(mklongname(n)) // ' in ' // TRIM(condavglongname(ncond)), &
           trim(mkunits(n)),ncond)
   end do
end do

! add microphysical process rates into model statstics +++mhwang
name = 'PRC'
longname = 'AUTOCONVERSION DROPLETS'
units = 'kg/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'PRA'
longname = 'ACCRETION DROPLETS BY RAIN'
units = 'kg/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'PSMLT'
longname = 'CHANGE Q MELTING SNOW TO RAIN'
units = 'kg/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'EVPMS'
longname = 'CHNAGE Q MELTING SNOW EVAPORATING'
units = 'kg/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'PRACS'
longname = 'CHANGE Q RAIN-SNOW COLLECTION'
units = 'kg/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'EVPMG'
longname = 'CHANGE Q MELTING OF GRAUPEL AND EVAPORATION'
units = 'kg/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'PRACG'
longname = 'CHANGE IN Q COLLECTION RAIN BY GRAUPEL'
units = 'kg/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'PRE'
longname = 'EAPORATION OF RAIN'
units = 'kg/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'PGMLT'
longname = 'CHANGE Q MELTING OF GRAUPEL'
units = 'kg/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'MNUCCC'
longname = 'CHANGE Q DUE TO CONTACT FREEZ DROPLETS'
units = 'kg/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'PSACWS'
longname = 'CHANGE Q DROPLET ACCRETION BY SNOW'
units = 'kg/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'PSACWI'
longname = 'CHANGE Q DROPLET ACCRETION BY CLOUD ICE'
units = 'kg/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'QMULTS'
longname = 'CHANGE Q DUE TO ICE MULT DROPLETS/SNOW'
units = 'kg/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'QMULTG'
longname = 'CHANGE Q DUE TO ICE MULT DROPLETS/GRAUPEL'
units = 'kg/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'PSACWG'
longname = 'CHANGE IN Q COLLECTION DROPLETS BY GRAUPEL'
units = 'kg/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'PGSACW'
longname = 'CONVERSION Q TO GRAUPEL DUE TO COLLECTION DROPLETS BY SNOW'
units = 'kg/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'PRD'
longname = 'DEP CLOUD ICE'
units = 'kg/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'PRCI'
longname = 'CHANGE Q AUTOCONVERSION CLOUD ICE BY SNOW'
units = 'kg/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'PRAI'
longname = 'CHANGE Q ACCRETION CLOUD ICE'
units = 'kg/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'QMULTR'
longname = 'CHANGE Q DUE TO ICE RAIN/SNOW'
units = 'kg/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'QMULTRG'
longname = 'CHANGE Q DUE TO ICE MULT RAIN/GRAUPEL'
units = 'kg/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'MNUCCD'
longname = 'CHANGE Q FREEZING AEROSOL (PRIM ICE NUCLEATION)'
units = 'kg/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'PRACI'
longname = 'CHANGE QI, ICE-RAIN COLLECTION'
units = 'kg/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'PRACIS'
longname = 'CHANGE QI, ICE RAIN COLLISION, ADDED TO SNOW'
units = 'kg/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'EPRD'
longname = 'SUBLIMATION CLOUD ICE'
units = 'kg/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'MNUCCR'
longname = 'CHANGE Q DUE TO CONTACT FREEZ RAIN'
units = 'kg/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'PIACR'
longname = 'CHANGE QR, ICE-RAIN COLLECTION'
units = 'kg/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'PIACRS'
longname = 'CHANGE QR, ICE RAIN COLLISION, ADDED TO SNOW'
units = 'kg/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'PGRACS'
longname = 'CONVERSION Q TO GRAUPEL DUE TO COLLECTION RAIN BY SNOW'
units = 'kg/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'PRDS'
longname = 'DEP OF SNOW'
units = 'kg/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'EPRDS'
longname = 'SUBLIMATION OF SNOW'
units = 'kg/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'PSACR'
longname = 'CONVERSION DUE TO COLL OF SNOW BY RAIN'
units = 'kg/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'PRDG'
longname = 'DEP OF GRAUPEL'
units = 'kg/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'EPRDG'
longname = 'SUBLIMATION OF GRAUPEL'
units = 'kg/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'NPRC1'
longname = 'CHANGE NR AUTOCONVERSION DROPLETS'
units = '#/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'NRAGG'
longname = 'SELF-COLLECTION OF RAIN'
units = '#/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'NPRACG'
longname = 'CHANGE N COLLECTION RAIN BY GRAUPEL'
units = '#/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'NSUBR'
longname = 'LOSS OF NR DURING EVAP'
units = '#/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'NSMLTR'
longname = 'CHANGE N MELTING SNOW TO RAIN'
units = '#/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'NGMLTR'
longname = 'CHANGE N MELTING GRAUPEL TO RAIN'
units = '#/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'NPRACS'
longname = 'CHANGE N RAIN-SNOW COLLECTION'
units = '#/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'NNUCCR'
longname = 'CHANGE N DUE TO CONTACT FREEZ RAIN'
units = '#/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'NIACR'
longname = 'CHANGE N, ICE-RAIN COLLECTION'
units = '#/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'NIACRS'
longname = 'CHANGE N, ICE RAIN COLLISION, ADDED TO SNOW'
units = '#/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'NGRACS'
longname = 'CHANGE N CONVERSION TO GRAUPEL DUE TO COLLECTION RAIN BY SNOW'
units = '#/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'NSMLTS'
longname = 'CHANGE N MELTING SNOW' 
units = '#/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'NSAGG'
longname = 'SELF-COLLECTION OF SNOW'
units = '#/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'NPRCI'
longname = 'CHANGE N AUTOCONVERSION CLOUD ICE BY SNOW'
units = '#/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'NSCNG'
longname = 'CHANGE N CONVERSION TO GRAUPEL DUE TO COLLECTION DROPLETS BY SNOW'
units = '#/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'NSUBS'
longname = 'LOSS OF NS DURING SUB.'
units = '#/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'PCC'
longname = 'SATURATION ADJUSTMENT'
units = 'kg/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'NNUCCC'
longname = 'CONTACT FREEZING OF CLOUD DROPS'
units = '#/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'NPSACWS'
longname = 'DROPLET ACCRETION BY SNOW'
units = '#/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'NPRA'
longname = 'DROPLET ACCRETION BY SNOW'
units = '#/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'NPRC'
longname = 'DROPLET ACCRETION BY RAIN'
units = '#/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'NPSACWI'
longname = 'DROPLET ACCRETION BY CLOUD ICE'
units = '#/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'NPSACWG'
longname = 'COLLECTION OF CLOUD DROPS BY GRAUPEL'
units = '#/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'NPRAI'
longname = 'ACCRETION OF CLOUD ICE BY SNOW'
units = '#/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'NMULTS'
longname = 'ICE MULT., RIMING OF CLOUD DROPS BY SNOW'
units = '#/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'NMULTG'
longname = 'ICE MULT., ACCRETION OF CLOUD DROPS BY GRAUPEL'
units = '#/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'NMULTR'
longname = 'ICE MULT., RIMING OF RAIN BY SNOW'
units = '#/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'NMULTRG'
longname = 'ICE MULT., ACCRETION OF RAIN BY GRAUPEL'
units = '#/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'NNUCCD'
longname = 'PRIMARY ICE NUC., FREEZING OF AEROSOLS'
units = '#/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'NSUBI'
longname = 'LOSS OF ICE DUE TO SUBLIMATION'
units = '#/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'NGMLTG'
longname = 'LOSS OF GRUPEL DUE TO MELTING'
units = '#/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'NSUBG'
longname = 'LOSS OF GRAUPEL DUE TO SUBLIMATION'
units = '#/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'NACT'
longname = 'CLOUD DROP FORMATION BY AEROSOL ACTIVATION'
units = '#/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'SIZEFIX_NR'
longname = 'ADJUST RAIN # FOR LARGE/SMALL DROPS'
units = '#/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'SIZEFIX_NC'
longname = 'ADJUST CLOUD # FOR LARGE/SMALL DROPS'
units = '#/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'SIZEFIX_NI'
longname = 'ADJUST ICE # FOR LARGE/SMALL DROPS'
units = '#/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'SIZEFIX_NS'
longname = 'ADJUST SNOW # FOR LARGE/SMALL DROPS'
units = '#/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'SIZEFIX_NG'
longname = 'ADJUST GRAUPEL # FOR LARGE/SMALL DROPS'
units = '#/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'NEGFIX_NR'
longname = 'ADJUST RAIN # TO REMOVE NEGATIVE VALUES'
units = '#/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'NEGFIX_NC'
longname = 'ADJUST CLOUD # TO REMOVE NEGATIVE VALUES'
units = '#/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'NEGFIX_NI'
longname = 'ADJUST ICE # TO REMOVE NEGATIVE VALUES'
units = '#/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'NEGFIX_NS'
longname = 'ADJUST SNOW # TO REMOVE NEGATIVE VALUES'
units = '#/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'NEGFIX_NG'
longname = 'ADJUST GRAUPEL # TO REMOVE NEGATIVE VALUES'
units = '#/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'NIM_MORR_CL'
longname = 'CLIP LARGE ICE NUMBER CONCENTRATIONS'
units = '#/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'QC_INST'
longname = 'INSTANTANEOUS CHANGES TO QC'
units = 'kg/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'QR_INST'
longname = 'INSTANTANEOUS CHANGES TO QR'
units = 'kg/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'QI_INST'
longname = 'INSTANTANEOUS CHANGES TO QI'
units = 'kg/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'QS_INST'
longname = 'INSTANTANEOUS CHANGES TO QNI'
units = 'kg/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'QG_INST'
longname = 'INSTANTANEOUS CHANGES TO QG'
units = 'kg/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'NC_INST'
longname = 'INSTANTANEOUS CHANGES TO NC'
units = '#/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'NR_INST'
longname = 'INSTANTANEOUS CHANGES TO NR'
units = '#/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'NI_INST'
longname = 'INSTANTANEOUS CHANGES TO NI'
units = '#/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'NS_INST'
longname = 'INSTANTANEOUS CHANGES TO NS'
units = '#/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'NG_INST'
longname = 'INSTANTANEOUS CHANGES TO NG'
units = '#/kg/s'
call add_to_namelist(count,microcount,name,longname,units,0)



#ifndef UWM_STATS
if(masterproc) then
   write(*,*) 'Added ', microcount, ' arrays to statistics for M2005 microphysics'
end if
#endif /*UWM_STATS*/

#ifdef UWM_STATS
    

!----------------------------
! Correlations 
!----------------------------

if(dopredictNc) then
  name = 'corr_w_ncl'
  longname = 'Correlation of vertical velocity and cloud droplet concentration'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_cor_w_ncl = count
end if ! do predictNc
    

  if (doprecip) then
    
    name = 'corr_w_qr'
    longname = 'Correlation of vertical velocity and rain water mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_w_qr = count
    
    name = 'corr_w_nr'
    longname = 'Correlation of vertical velocity and rain drop concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_w_nr = count
    
  end if ! doprecip

  if (doicemicro) then
    
    name = 'corr_w_qci'
    longname = 'Correlation of vertical velocity and cloud ice mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_w_qci = count
    
    name = 'corr_w_nci'
    longname = 'Correlation of vertical velocity and cloud ice concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_w_nci = count
    
    name = 'corr_w_qs'
    longname = 'Correlation of vertical velocity and snow mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_w_qs = count
    
    name = 'corr_w_ns'
    longname = 'Correlation of vertical velocity and snowflake concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_w_ns = count
    

    if (dograupel) then
    
      name = 'corr_w_qg'
      longname = 'Correlation of vertical velocity and graupel mixing ratio'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_cor_w_qg = count
    
      name = 'corr_w_ng'
      longname = 'Correlation of vertical velocity and graupel concentration'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_cor_w_ng = count
    
    end if ! dograupel
  end if ! doicemicro

if(dopredictNc) then
  if (doprecip) then
        
    name = 'corr_ncl_qr'
    longname = 'Correlation of cloud droplet conc. and rain water mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_ncl_qr = count
        
    name = 'corr_ncl_nr'
    longname = 'Correlation of cloud droplet conc. and rain drop conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_ncl_nr = count
        
  end if ! doprecip
    
  if (doicemicro) then
        
    name = 'corr_ncl_qci'
    longname = 'Correlation of cloud droplet conc. and cloud ice mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_ncl_qci = count
        
    name = 'corr_ncl_nci'
    longname = 'Correlation of cloud droplet conc. and cloud ice conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_ncl_nci = count
       
    name = 'corr_ncl_qs'
    longname = 'Correlation of cloud droplet conc. and snow mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_ncl_qs = count
      
    name = 'corr_ncl_ns'
    longname = 'Correlation of cloud droplet conc. and snowflake conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_ncl_ns = count
        
    if (dograupel) then
        
      name = 'corr_ncl_qg'
      longname = 'Correlation of cloud droplet conc. and graupel mixing ratio'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_cor_ncl_qg = count
        
      name = 'corr_ncl_ng'
      longname = 'Correlation of cloud droplet conc. and graupel conc.'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_cor_ncl_ng = count
        
    end if ! dograupel
  end if ! doicemicro 
end if ! dopredictNc

if (doprecip) then

  if (dopredictNc) then
    name = 'corr_qr_nr'
    longname = 'Correlation of rain water mixing ratio and rain drop concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_qr_nr = count
  end if
  
  if(doicemicro) then

    name = 'corr_qr_qci'
    longname = 'Correlation of rain water mixing ratio and cloud ice mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_qr_qci = count
        
    name = 'corr_qr_nci'
    longname = 'Correlation of rain water mixing ratio and cloud ice conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_qr_nci = count
        
    name = 'corr_qr_qs'
    longname = 'Correlation of rain water mixing ratio and snow mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_qr_qs = count
        
    name = 'corr_qr_ns'
    longname = 'Correlation of rain water mixing ratio and snowflake conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_qr_ns = count
        

    if(dograupel) then

      name = 'corr_qr_qg'
      longname = 'Correlation of rain water mixing ratio and graupel mixing ratio'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_cor_qr_qg = count
        
      name = 'corr_qr_ng'
      longname = 'Correlation of rain water mixing ratio and graupel conc.'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_cor_qr_ng = count
        
    end if !dograupel
        
    name = 'corr_nr_qci'
    longname = 'Correlation of rain drop conc. and cloud ice mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_nr_qci = count
        
    name = 'corr_nr_nci'
    longname = 'Correlation of rain drop conc. and cloud ice conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_nr_nci = count
        
    name = 'corr_nr_qs'
    longname = 'Correlation of rain drop conc. and snow mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_nr_qs = count
        
    name = 'corr_nr_ns'
    longname = 'Correlation of rain drop conc. and snowflake conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_nr_ns = count
        
      if(dograupel) then

        name = 'corr_nr_qg'
        longname = 'Correlation of rain drop conc. and graupel mixing ratio'
        units = '-'
        call add_to_namelist(count,microcount,name,longname,units,corr_avg)
        ! Note:  it is important to place these after the call to add_to_namelist.
        idx_cor_nr_qg = count
        
        name = 'corr_nr_ng'
        longname = 'Correlation of rain drop conc. and graupel conc.'
        units = '-'
        call add_to_namelist(count,microcount,name,longname,units,corr_avg)
        ! Note:  it is important to place these after the call to add_to_namelist.
        idx_cor_nr_ng = count
        
      end if !dograupel
    end if !doicemicro
end if !doprecip

if (doicemicro) then
    
  name = 'corr_qci_nci'
  longname = 'Correlation of cloud ice mixing ratio and cloud ice conc.'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_cor_qci_nci = count
    
  name = 'corr_qci_qs'
  longname = 'Correlation of cloud ice mixing ratio and snow mixing ratio'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_cor_qci_qs = count
    
  name = 'corr_qci_ns'
  longname = 'Correlation of cloud ice mixing ratio and snowflake conc.'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_cor_qci_ns = count
    
  if (dograupel) then

    name = 'corr_qci_qg'
    longname = 'Correlation of cloud ice mixing ratio and graupel mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_qci_qg = count
        
    name = 'corr_qci_ng'
    longname = 'Correlation of cloud ice mixing ratio and graupel conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_qci_ng = count
        
  end if !dograupel
    
  name = 'corr_nci_qs'
  longname = 'Correlation of cloud ice concentration and snow mixing ratio'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_cor_nci_qs = count
    
  name = 'corr_nci_ns'
  longname = 'Correlation of cloud ice concentration and snowflake conc.'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_cor_nci_ns = count
    
  if (dograupel) then

    name = 'corr_nci_qg'
    longname = 'Correlation of cloud ice conc. and graupel mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_nci_qg = count
        
    name = 'corr_nci_ng'
    longname = 'Correlation of cloud ice conc. and graupel conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_nci_ng = count
        
  end if !dograupel
  
  name = 'corr_qs_ns'
  longname = 'Correlation of snow mixing ratio and snowflake concentration'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_cor_qs_ns = count
    
  if (dograupel) then

    name = 'corr_qs_qg'
    longname = 'Correlation of snow mixing ratio and graupel mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_qs_qg = count
        
    name = 'corr_qs_ng'
    longname = 'Correlation of snow mixing ratio and graupel concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_qs_ng = count
        
    name = 'corr_ns_qg'
    longname = 'Correlation of snowflake concentration and graupel mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_ns_qg = count
        
    name = 'corr_ns_ng'
    longname = 'Correlation of snowflake concentration and graupel conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_ns_ng = count
        
    name = 'corr_qg_ng'
    longname = 'Correlation of graupel mixing ratio and graupel concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_qg_ng = count
        
  end if !dograupel

end if !doicemicro

!----------------------------
! Fractions
! 
! The _em# specifies the threshold (kg/kg) for each species.
!----------------------------


name = 'cloudliq_frac_em8'
longname = 'cloud liquid fraction'
units = '-'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'cloudliq_frac_em7'
longname = 'cloud liquid fraction'
units = '-'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'cloudliq_frac_em6'
longname = 'cloud liquid fraction'
units = '-'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'cloudliq_frac_em5'
longname = 'cloud liquid fraction'
units = '-'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'cloudliq_frac_em4'
longname = 'cloud liquid fraction'
units = '-'
call add_to_namelist(count,microcount,name,longname,units,0)


if (doprecip) then
  
    name = 'rain_frac_em8'
    longname = 'rain fraction'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,0)

    name = 'rain_frac_em7'
    longname = 'rain fraction'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,0)

    name = 'rain_frac_em6'
    longname = 'rain fraction'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,0)

    name = 'rain_frac_em5'
    longname = 'rain fraction'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,0)

    name = 'rain_frac_em4'
    longname = 'rain fraction'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,0)

    if (doicemicro) then
      
      name = 'cloudice_frac_em8'
      longname = 'cloud ice fraction'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,0)
    
      name = 'cloudice_frac_em7'
      longname = 'cloud ice fraction'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,0)
    
      name = 'cloudice_frac_em6'
      longname = 'cloud ice fraction'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,0)
    
      name = 'cloudice_frac_em5'
      longname = 'cloud ice fraction'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,0)
    
      name = 'cloudice_frac_em4'
      longname = 'cloud ice fraction'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,0)
    
      
      name = 'snow_frac_em8'
      longname = 'snow fraction'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,0)

      name = 'snow_frac_em7'
      longname = 'snow fraction'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,0)

      name = 'snow_frac_em6'
      longname = 'snow fraction'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,0)

      name = 'snow_frac_em5'
      longname = 'snow fraction'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,0)

      name = 'snow_frac_em4'
      longname = 'snow fraction'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,0)

      if (dograupel) then
    
        name = 'graupel_frac_em8'
        longname = 'graupel fraction'
        units = '-'
        call add_to_namelist(count,microcount,name,longname,units,0)

        name = 'graupel_frac_em7'
        longname = 'graupel fraction'
        units = '-'
        call add_to_namelist(count,microcount,name,longname,units,0)

        name = 'graupel_frac_em6'
        longname = 'graupel fraction'
        units = '-'
        call add_to_namelist(count,microcount,name,longname,units,0)

        name = 'graupel_frac_em5'
        longname = 'graupel fraction'
        units = '-'
        call add_to_namelist(count,microcount,name,longname,units,0)

        name = 'graupel_frac_em4'
        longname = 'graupel fraction'
        units = '-'
        call add_to_namelist(count,microcount,name,longname,units,0)

      end if !dograupel
    end if !doicemicro
end if !doprecip

!----------------------------
! Domain-wide variances
!----------------------------

name = 'ncloudliqp2'
longname = 'Domain-wide variance of cloud liquid number concentration'
units = '[(#/kg)^2]'
call add_to_namelist(count,microcount,name,longname,units,0)

if (doprecip) then
  
    name = 'nrainp2'
    longname = 'Domain-wide variance of rain number concentration'
    units = '[(#/kg)^2]'
    call add_to_namelist(count,microcount,name,longname,units,0)

    name = 'qrainp2'
    longname = 'Domain-wide variance of rain mixing ratio'
    units = '[(kg/kg)^2]'
    call add_to_namelist(count,microcount,name,longname,units,0)

    if (doicemicro) then
      
      name = 'ncloudicep2'
      longname = 'Domain-wide variance of cloud ice number concentration'
      units = '[(#/kg)^2]'
      call add_to_namelist(count,microcount,name,longname,units,0)
    
      name = 'qcloudicep2'
      longname = 'Domain-wide variance of cloud ice mixing ratio'
      units = '[(kg/kg)^2]'
      call add_to_namelist(count,microcount,name,longname,units,0)
      
      name = 'nsnowp2'
      longname = 'Domain-wide variance of snow number concentration'
      units = '[(#/kg)^2]'
      call add_to_namelist(count,microcount,name,longname,units,0)

      name = 'qsnowp2'
      longname = 'Domain-wide variance of snow mixing ratio'
      units = '[(kg/kg)^2]'
      call add_to_namelist(count,microcount,name,longname,units,0)

      if (dograupel) then
    
        name = 'ngraupelp2'
        longname = 'Domain-wide variance of graupel number concentration'
        units = '[(#/kg)^2]'
        call add_to_namelist(count,microcount,name,longname,units,0)

        name = 'qgraupelp2'
        longname = 'Domain-wide variance of graupel mixing ratio'
        units = '[(kg/kg)^2]'
        call add_to_namelist(count,microcount,name,longname,units,0)

      end if !dograupel
    end if !doicemicro
  end if !doprecip

!----------------------------
! In-precip means and variances
!----------------------------

name = 'ncloudliqp2_ip'
longname = 'Within cloud variance of cloud liquid number concentration'
units = '[(#/kg)^2]'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'ncloudliqm_ip'
longname = 'Within cloud mean of cloud liquid number concentration'
units = '[(#/kg)]'
call add_to_namelist(count,microcount,name,longname,units,0)

if (doprecip) then
  
    name = 'nrainp2_ip'
    longname = 'Within rain variance of rain number concentration'
    units = '[(#/kg)^2]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
    name = 'nrainm_ip'
    longname = 'Within rain mean of rain number concentration'
    units = '[(#/kg)]'
    call add_to_namelist(count,microcount,name,longname,units,0)

    name = 'qrainp2_ip'
    longname = 'Within rain variance of rain mixing ratio'
    units = '[(kg/kg)^2]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
    name = 'qrainm_ip'
    longname = 'Within rain mean of rain mixing ratio'
    units = '[(kg/kg)]'
    call add_to_namelist(count,microcount,name,longname,units,0)

    if (doicemicro) then
      
      name = 'ncloudicep2_ip'
      longname = 'Within cloudice variance of cloud ice number concentration'
      units = '[(#/kg)^2]'
      call add_to_namelist(count,microcount,name,longname,units,0)
      
      name = 'ncloudicem_ip'
      longname = 'Within cloudice mean of cloud ice number concentration'
      units = '[(#/kg)]'
      call add_to_namelist(count,microcount,name,longname,units,0)
    
      name = 'qcloudicep2_ip'
      longname = 'Within cloudice variance of cloud ice mixing ratio'
      units = '[(kg/kg)^2]'
      call add_to_namelist(count,microcount,name,longname,units,0)
      
      name = 'qcloudicem_ip'
      longname = 'Within cloudice mean of cloud ice mixing ratio'
      units = '[(kg/kg)]'
      call add_to_namelist(count,microcount,name,longname,units,0)
      
      name = 'nsnowp2_ip'
      longname = 'Within snow variance of snow number concentration'
      units = '[(#/kg)^2]'
      call add_to_namelist(count,microcount,name,longname,units,0)
      
      name = 'nsnowm_ip'
      longname = 'Within snow mean of snow number concentration'
      units = '[(#/kg)]'
      call add_to_namelist(count,microcount,name,longname,units,0)

      name = 'qsnowp2_ip'
      longname = 'Within snow variance of snow mixing ratio'
      units = '[(kg/kg)^2]'
      call add_to_namelist(count,microcount,name,longname,units,0)
      
      name = 'qsnowm_ip'
      longname = 'Within snow mean of snow mixing ratio'
      units = '[(kg/kg)]'
      call add_to_namelist(count,microcount,name,longname,units,0)

      if (dograupel) then
    
        name = 'ngraupelp2_ip'
        longname = 'Within graupel variance of graupel number concentration'
        units = '[(#/kg)^2]'
        call add_to_namelist(count,microcount,name,longname,units,0)
        
        name = 'ngraupelm_ip'
        longname = 'Within graupel mean of graupel number concentration'
        units = '[(#/kg)]'
        call add_to_namelist(count,microcount,name,longname,units,0)

        name = 'qgraupelp2_ip'
        longname = 'Within graupel variance of graupel mixing ratio'
        units = '[(kg/kg)^2]'
        call add_to_namelist(count,microcount,name,longname,units,0)
        
        name = 'qgraupelm_ip'
        longname = 'Within graupel mean of graupel mixing ratio'
        units = '[(kg/kg)]'
        call add_to_namelist(count,microcount,name,longname,units,0)

      end if !dograupel
    end if !doicemicro
  end if !doprecip

!----------------------------
! Covariances
!---------------------------

if(dopredictNc) then
  
  name = 'covarnce_w_ncl'
  longname = 'Covariance of vertical velocity and cloud droplet concentration'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,0)

end if

if (doprecip) then
    
  name = 'covarnce_w_qr'
  longname = 'Covariance of vertical velocity and rain water mixing ratio'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'covarnce_w_nr'
  longname = 'Covariance of vertical velocity and rain drop concentration'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,0)

end if ! doprecip

if (doicemicro) then
    
  name = 'covarnce_w_qci'
  longname = 'Covariance of vertical velocity and cloud ice mixing ratio'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'covarnce_w_nci'
  longname = 'Covariance of vertical velocity and cloud ice concentration'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'covarnce_w_qs'
  longname = 'Covariance of vertical velocity and snow mixing ratio'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'covarnce_w_ns'
  longname = 'Covariance of vertical velocity and snowflake concentration'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,0)

end if ! doicemicro

if (dograupel) then
    
  name = 'covarnce_w_qg'
  longname = 'Covariance of vertical velocity and graupel mixing ratio'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'covarnce_w_ng'
  longname = 'Covariance of vertical velocity and graupel concentration'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,0)

end if ! dograupel

if(dopredictNc) then
    if (doprecip) then
        
      name = 'covarnce_ncl_qr'
      longname = 'Covariance of cloud droplet conc. and rain water mixing ratio'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,0)
        
      name = 'covarnce_ncl_nr'
      longname = 'Covariance of cloud droplet conc. and rain drop conc.'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,0)
    
    end if ! doprecip
    
    if (doicemicro) then
        
      name = 'covarnce_ncl_qci'
      longname = 'Covariance of cloud droplet conc. and cloud ice mixing ratio'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,0)
        
      name = 'covarnce_ncl_nci'
      longname = 'Covariance of cloud droplet conc. and cloud ice conc.'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,0)
        
      name = 'covarnce_ncl_qs'
      longname = 'Covariance of cloud droplet conc. and snow mixing ratio'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,0)
    
      name = 'covarnce_ncl_ns'
      longname = 'Covariance of cloud droplet conc. and snowflake conc.'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,0)
    
    end if ! doicemicro
    
    if (dograupel) then
        
      name = 'covarnce_ncl_qg'
      longname = 'Covariance of cloud droplet conc. and graupel mixing ratio'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,0)
        
      name = 'covarnce_ncl_ng'
      longname = 'Covariance of cloud droplet conc. and graupel conc.'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,0)

    end if ! dograupel
end if ! do predictNc

if (doprecip) then

  if (dopredictNc) then
    name = 'covarnce_qr_nr'
    longname = 'Covariance of rain water mixing ratio and rain drop concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,0)
  end if 

  if(doicemicro) then

    name = 'covarnce_qr_qci'
    longname = 'Covariance of rain water mixing ratio and cloud ice mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,0)
        
    name = 'covarnce_qr_nci'
    longname = 'Covariance of rain water mixing ratio and cloud ice conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
    name = 'covarnce_qr_qs'
    longname = 'Covariance of rain water mixing ratio and snow mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
    name = 'covarnce_qr_ns'
    longname = 'Covariance of rain water mixing ratio and snowflake conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,0)
  
  end if !doicemicro

  if(dograupel) then

    name = 'covarnce_qr_qg'
    longname = 'Covariance of rain water mixing ratio and graupel mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,0)
        
    name = 'covarnce_qr_ng'
    longname = 'Covariance of rain water mixing ratio and graupel conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
  end if !dograupel

  if(doicemicro) then
        
    name = 'covarnce_nr_qci'
    longname = 'Covariance of rain drop conc. and cloud ice mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,0)
        
    name = 'covarnce_nr_nci'
    longname = 'Covariance of rain drop conc. and cloud ice conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
    name = 'covarnce_nr_qs'
    longname = 'Covariance of rain drop conc. and snow mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
    name = 'covarnce_nr_ns'
    longname = 'Covariance of rain drop conc. and snowflake conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,0)
  
  end if !doicemicro

  if(dograupel) then

    name = 'covarnce_nr_qg'
    longname = 'Covariance of rain drop conc. and graupel mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,0)
        
    name = 'covarnce_nr_ng'
    longname = 'Covariance of rain drop conc. and graupel conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
  end if !dograupel
end if !doprecip

if (doicemicro) then
    
  name = 'covarnce_qci_nci'
  longname = 'Covariance of cloud ice mixing ratio and cloud ice conc.'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'covarnce_qci_qs'
  longname = 'Covariance of cloud ice mixing ratio and snow mixing ratio'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'covarnce_qci_ns'
  longname = 'Covariance of cloud ice mixing ratio and snowflake conc.'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,0)

  if (dograupel) then

    name = 'covarnce_qci_qg'
    longname = 'Covariance of cloud ice mixing ratio and graupel mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,0)

    name = 'covarnce_qci_ng'
    longname = 'Covariance of cloud ice mixing ratio and graupel conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
  end if !dograupel
    
  name = 'covarnce_nci_qs'
  longname = 'Covariance of cloud ice concentration and snow mixing ratio'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'covarnce_nci_ns'
  longname = 'Covariance of cloud ice concentration and snowflake conc.'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,0)

  if (dograupel) then

    name = 'covarnce_nci_qg'
    longname = 'Covariance of cloud ice conc. and graupel mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,0)

    name = 'covarnce_nci_ng'
    longname = 'Covariance of cloud ice conc. and graupel conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,0)
  end if !dograupel
    
  name = 'covarnce_qs_ns'
  longname = 'Covariance of snow mixing ratio and snowflake concentration'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,0)

  if (dograupel) then

    name = 'covarnce_qs_qg'
    longname = 'Covariance of snow mixing ratio and graupel mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,0)

    name = 'covarnce_qs_ng'
    longname = 'Covariance of snow mixing ratio and graupel concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
    name = 'covarnce_ns_qg'
    longname = 'Covariance of snowflake concentration and graupel mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,0)

    name = 'covarnce_ns_ng'
    longname = 'Covariance of snowflake concentration and graupel conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
    name = 'covarnce_qg_ng'
    longname = 'Covariance of graupel mixing ratio and graupel concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
  end if !dograupel
end if !doicemicro

if(masterproc) then
   write(*,*) 'Added ', microcount, ' arrays to statistics for M2005 microphysics'
end if
#endif /*UWM_STATS*/

end subroutine micro_hbuf_init

!----------------------------------------------------------------------
!!!! Collect microphysics history statistics (vertical profiles)
!! Note that only the fields declared in micro_hbuf_init() are allowed to
! be collected

subroutine micro_statistics()

use vars
#ifdef UWM_STATS
use hbuffer, only: hbuf_put, hbuf_put_level
#else
use hbuffer, only: hbuf_put
#endif /*UWM_STATS*/
use params, only : lcond

#ifdef UWM_STATS
use compute_correlation_module, only: corravg_count, sum_value_xy_domain,& 
                                      mean_xy_domain, variance_xy_domain,&
                                      mean_ip_xy_domain, variance_ip_xy_domain,&
                                      covariance_xy_domain, compute_correlation, undef_corr  

!Microphysics correlations and covariances matricies
real, dimension(11,11,nzm) :: micro_correlations 
real, dimension(11,11,nzm) :: micro_covarnce 

!Vertical velocity, interpolated, domain means, and variances
real, dimension(nx,ny,nzm) :: w_zt !vertical velocity interpolated on the scalar grid
real, dimension(nzm) :: domain_mean_w_zt !domain average of w_zt
real, dimension(nzm) :: domain_varnce_w_zt !domain variance of w_zt

!Microphysical domain-wide means and variances 
real, dimension(nmicro_fields,nzm) :: domain_mean_micro !domain averages of micro_fields
real, dimension(nmicro_fields,nzm) :: domain_varnce_micro !domain variance of micro_fields

!Microphysical within-'precip' means and variances 
real, dimension(nmicro_fields,nzm) :: mean_micro_ip ! within-'precip' averages of micro_fields
real, dimension(nmicro_fields,nzm) :: varnce_micro_ip ! within-'precip' of micro_fields

integer :: idx_w,& !index for vertical velocity
         idx_microstart !index for first microfield variable (2) 
!--------------------------------------------------------------------------------------------------
!Ensemble fractions 
!
!  These are used to see the effect of tolerance values on the in-'precip'
!  variances and fractions. The magic number 5 refers to the number of tolerance
!  values we test. Starting from SAM's defaut 1e-6 [kg/kg] for all micro.
!  species, we modulate the threshold +/- 2 orders of magnitude. 

!  Any variables that are reliant on these thresholds are defined in this
!  section

!Binary mask of in(1) and out(0) of micro. species.
real, dimension(nx,ny,nzm,5) :: cloudliq_mask, rain_mask, cloudice_mask, snow_mask, graupel_mask

!Number of gridpoints within micro. species
integer, dimension(nzm,5) :: cloudliq_sum, rain_sum, cloudice_sum, snow_sum, graupel_sum

!Micro. fractions
real, dimension(nzm,5) :: cloudliq_frac, rain_frac, cloudice_frac, snow_frac, graupel_frac 

!Indicies for vars computed with 1e-8 [kg/kg, 1e-7, etc. thresholds
integer :: em8, em7, em6, em5, em4, thresh_index

!To make code more readable.
real ::  frac_threshold, frac_threshold_init
integer :: order_of_magnitude, ensemble_frac_iter_start, ensemble_frac_iter_end

!A switch on which threshold to compute in-'precip' variances
integer :: thresh_out
!--------------------------------------------------------------------------------------------------
#endif /*UWM_STATS*/
real, dimension(nzm) :: tr0, tr2

real tmp(2), factor_xy
integer i,j,k,m, n, ii, jj, nn, ncond

call t_startf ('micro_statistics')

factor_xy = 1./float(nx*ny)

do n = 1,nmicro_fields
   do k = 1,nzm
      tmp(1) = dz
      tmp(2) = dz/dtn
      tr0(k) = SUM(micro_field(1:nx,1:ny,k,n))
      tr2(k) = SUM(micro_field(1:nx,1:ny,k,n)*micro_field(1:nx,1:ny,k,n))
      mkwle(k,n) = mkwle(k,n)*tmp(2)*lfac(n) ! resolved flux
      mkwsb(k,n) = mkwsb(k,n)*tmp(1)*lfac(n) ! subgrid flux
      mksed(k,n) = mksed(k,n)*lfac(n) ! sedimentation flux

      mstor(k, n) = SUM(micro_field(1:nx,1:ny,k,n))-mstor(k,n)
   end do

   if(flag_wmass(n).lt.1) then
      ! remove factor of rho from number concentrations
      tr0(:) = tr0(:)*rho(:)
      tr2(:) = tr2(:)*rho(:)**2
      mkadv(1:nzm,n) = mkadv(1:nzm,n)*rho(:)
      mkdiff(1:nzm,n) = mkdiff(1:nzm,n)*rho(:)
      mtend(1:nzm,n) = mtend(1:nzm,n)*rho(:)
      stend(1:nzm,n) = stend(1:nzm,n)*rho(:)
      mklsadv(1:nzm,n) = mklsadv(1:nzm,n)*rho(:)

      mstor(1:nzm,n) = mstor(1:nzm,n)*rho(:)
   end if

!bloss/qt: output all microphysical fields
!   if(n.ne.iqv) then
   ! mean microphysical field
   call hbuf_put(trim(mkname(n)),tr0,mkoutputscale(n)*factor_xy)
  if(n.eq.iqv) then
      ! variance of microphysical field,  only for QTO (qv+qcl)
      call hbuf_put(trim(mkname(n))//'2',tr2,mkoutputscale(n)**2*factor_xy)
   end if

   ! do not rescale fluxes
   call hbuf_put(trim(mkname(n))//'FLX',mkwle(1,n),factor_xy)
   call hbuf_put(trim(mkname(n))//'FLXS',mkwsb(1,n),factor_xy)
   call hbuf_put(trim(mkname(n))//'SDFLX',mksed(1,n),factor_xy)

   ! tendencies
   call hbuf_put(trim(mkname(n))//'ADV', &
        mkadv(:,n),mkoutputscale(n)*factor_xy*86400./dtn)
   call hbuf_put(trim(mkname(n))//'DIFF', &
        mkdiff(:,n),mkoutputscale(n)*factor_xy*86400./dtn)
   call hbuf_put(trim(mkname(n))//'LSADV', &
        mklsadv(:,n),mkoutputscale(n)*factor_xy*86400.)
   call hbuf_put(trim(mkname(n))//'MPHY', &
        mtend(:,n),mkoutputscale(n)*factor_xy*86400.)
   call hbuf_put(trim(mkname(n))//'SED', &
        stend(:,n),mkoutputscale(n)*factor_xy*86400.)

   ! Storage terms
   call hbuf_put(trim(mkname(n))//'STO', &
        mstor(:,n),mkoutputscale(n)*factor_xy*86400./dtn)

   if((flag_wmass(n).gt.0).and.(n.ne.iqv)) then
      ! fractional area of microphysical field > 1.e-6
      call hbuf_put(trim(mkname(n))//'FRAC',mfrac(1,n),factor_xy)

      ! approx optical depth
      call hbuf_put('TAU'//trim(mkname(n)),trtau(:,n),factor_xy)

      !bloss (Apr 09):  This measure of effective radius is unreliable if the 
      !          microphysical fraction is not roughly uniform across
      !          the processors in a MPI run.  As a result, I am
      !          removing it from the outputs.  It is reliable if computed from
      !          the quantities TAU* and Q* in the output file.
!bloss      ! effective radius
!bloss      tr2(:) = 0.
!bloss      if(trtau(1,n).gt.0.) then
!bloss         tr2(1) = 1.e6*0.0018*rho(1)*dz*adz(1)*tr0(1)/trtau(1,n)
!bloss      end if

!bloss      do k = 2,nzm
!bloss         if(trtau(k,n).gt.trtau(k-1,n)) then
!bloss            tr2(k) = 1.e6*0.0018*rho(k)*dz*adz(k)*tr0(k)/(trtau(k,n)-trtau(k-1,n))
!bloss         end if
!bloss      end do
!bloss      call hbuf_put('EFFR'//trim(mkname(n)),tr2,1.)

      !bloss (Apr 09): Make an alternate statistic that can be used
      ! to easily compute the mean effective radius in a consistent
      ! way from optical depth.  This quantity Q*OEFFR is essentially
      ! the layer optical depth scaled in such a way that
      !
      !    EFFR = <Q*> / <Q*OEFFR>
      !
      ! where <.> is a time- and horizontal average.
      tr2(:) = 0.
      tr2(1) = trtau(1,n) / (1.e6*0.0018*rho(1)*dz*adz(1)*1.e-3)
      do k = 2,nzm
            tr2(k) = (trtau(k,n)-trtau(k-1,n)) / (1.e6*0.0018*rho(k)*dz*adz(k)*1.e-3) 
      end do
      call hbuf_put(trim(mkname(n))//'OEFFR',tr2,factor_xy)
      
   end if

   do ncond = 1,ncondavg
      do k = 1,nzm
         tr0(k) = SUM(micro_field(1:nx,1:ny,k,n)*condavg_mask(1:nx,1:ny,k,ncond))
      end do
      if(flag_number(n).eq.1) tr0(:) = tr0(:)*rho(:) ! remove factor of rho from number concentrations
      call hbuf_put(TRIM(mkname(n)) // TRIM(condavgname(ncond)), &
           tr0,mkoutputscale(n))
   end do

end do

!bloss/qt: in total water formulation, fluxes of qv and qcl computed together.
tr0(:) = mkwle(1:nzm,iqv) + mkwsb(1:nzm,iqv) ! qv + qcl tendencies
if(doicemicro) then
   tr0(:) = tr0(:) + mkwle(1:nzm,iqci) + mkwsb(1:nzm,iqci)
end if
call hbuf_put('QTFLUX',tr0,factor_xy)

!bloss/qt: Can't be computed reliably in total water formulation.
!bloss do k = 1,nzm
!bloss    tr0(k) = SUM(tmtend3d(1:nx,1:ny,k))
!bloss end do
!bloss call hbuf_put('QLAT',tr0,factor_xy*86400.)

!bloss do ncond = 1,ncondavg
!bloss    do k = 1,nzm
!bloss       tr0(k) = SUM(tmtend3d(1:nx,1:ny,k)*condavg_mask(1:nx,1:ny,k,ncond))
!bloss    end do
!bloss    call hbuf_put('QLAT' // TRIM(condavgname(ncond)),tr0,86400.)
!bloss end do

!bloss/qt: add separate output for cloud liquid water
!           and approx cloud liquid optical depth.
do k = 1,nzm
  tr0(k) = SUM(cloudliq(1:nx,1:ny,k))
end do
call hbuf_put('QC',tr0,factor_xy*mkoutputscale(iqv))

call hbuf_put('TAUQC',trtau(:,iqv),factor_xy)

!bloss (Apr 09): Make an alternate statistic that can be used
! to easily compute the mean effective radius in a consistent
! way from optical depth.  This quantity Q*OEFFR is essentially
! the layer optical depth scaled in such a way that
!
!    EFFR = <Q*> / <Q*OEFFR>
!
! where <.> is a time- and horizontal average.
tr2(:) = 0.
tr2(1) = trtau(1,iqv) / (1.e6*0.0018*rho(1)*dz*adz(1)*1.e-3)
do k = 2,nzm
  tr2(k) = (trtau(k,iqv)-trtau(k-1,iqv)) / (1.e6*0.0018*rho(k)*dz*adz(k)*1.e-3) 
end do
call hbuf_put('QCOEFFR',tr2,factor_xy)

!bloss/qt: add separate conditional averages for cloud liquid water and vapor.
do ncond = 1,ncondavg
   do k = 1,nzm
      tr0(k) = SUM(cloudliq(1:nx,1:ny,k)*condavg_mask(1:nx,1:ny,k,ncond))
   end do
   call hbuf_put('QC' // TRIM(condavgname(ncond)),tr0,mkoutputscale(iqv))
   do k = 1,nzm
      tr0(k) = SUM((micro_field(1:nx,1:ny,k,iqv)-cloudliq(1:nx,1:ny,k))*condavg_mask(1:nx,1:ny,k,ncond))
   end do
   call hbuf_put('QV' // TRIM(condavgname(ncond)),tr0,mkoutputscale(iqv))
end do

if(dopredictNc) then
  nn = 0.
  tmp(1)=0.
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      if(qcl(i,j,k).gt.0.) then
         tmp(1) = tmp(1) + micro_field(i,j,k,incl)*1.e-6
         nn = nn + 1
       end if
    end do
   end do      
  end do
  if (nn.gt.0) ncmn = ncmn + tmp(1)/dble(nn)
else
  ncmn = Nc0
end if
if(doprecip) then
  nn = 0.
  tmp(1)=0.
  do k=1,nzm
   do j=1,ny
    do i=1,nx 
      if(micro_field(i,j,k,iqr).gt.0.) then 
         tmp(1) = tmp(1) + micro_field(i,j,k,inr)*1.e-6
         nn = nn + 1
       end if
    end do
   end do
  end do
  if (nn.gt.0) then
      nrainy = nrainy + 1
      nrmn = nrmn + tmp(1)/dble(nn)
  end if
else
  nrmn = 0.
end if


#ifdef UWM_STATS
!Ensemble fractions. Indicies for thresholds of 1e-8 [kg/kg], 1e-7, etc.
 em8 = 1 
 em7 = 2
 em6 = 3
 em5 = 4
 em4 = 5
 
 !Compute fractions and in-'precip' variances using the 1e-6 [kg/kg] threshold
 thresh_out = em6
 
 frac_threshold_init = 1e-8
 ensemble_frac_iter_start= 1
 ensemble_frac_iter_end  = 5

 order_of_magnitude = 10

!Initial
frac_threshold = frac_threshold_init

do thresh_index=ensemble_frac_iter_start,ensemble_frac_iter_end

  do k=1,nzm
    do j=1,ny
      do i=1,nx
        ! We only want cloud water, not cloud water + rain water.
        if ( cloudliq(i,j,k)  > frac_threshold) then
          cloudliq_mask(i,j,k,thresh_index) = 1.
        else
          cloudliq_mask(i,j,k,thresh_index) = 0.
        end if
        
        if (doprecip) then
          if (micro_field(i,j,k,iqr) > frac_threshold) then
            rain_mask(i,j,k,thresh_index) = 1.
          else
            rain_mask(i,j,k,thresh_index) = 0.
          end if
          
          if (doicemicro) then
            if (micro_field(i,j,k,iqci) > frac_threshold) then
              cloudice_mask(i,j,k,thresh_index) = 1.
            else
              cloudice_mask(i,j,k,thresh_index) = 0.
            end if

            if (micro_field(i,j,k,iqs) > frac_threshold) then
              snow_mask(i,j,k,thresh_index) = 1.
            else
              snow_mask(i,j,k,thresh_index) = 0.
            end if
          
            if (dograupel) then
              if (micro_field(i,j,k,iqg) > frac_threshold) then
                graupel_mask(i,j,k,thresh_index) = 1.
              else
                graupel_mask(i,j,k,thresh_index) = 0.
              end if
            
            end if!do graupel  
          end if!do icemicro
        end if!doprecip
      
      end do !i
    end do !j
  end do !k
  
  cloudliq_sum(:,thresh_index) = sum_value_xy_domain(nzm,cloudliq_mask(1:nx,1:ny,1:nzm,thresh_index))
  rain_sum(:,thresh_index)     = sum_value_xy_domain(nzm,rain_mask(1:nx,1:ny,1:nzm,thresh_index))
  cloudice_sum(:,thresh_index) = sum_value_xy_domain(nzm,cloudice_mask(1:nx,1:ny,1:nzm,thresh_index))
  snow_sum(:,thresh_index)     = sum_value_xy_domain(nzm,snow_mask(1:nx,1:ny,1:nzm,thresh_index))
  graupel_sum(:,thresh_index)  = sum_value_xy_domain(nzm,graupel_mask(1:nx,1:ny,1:nzm,thresh_index))

  cloudliq_frac(:,thresh_index) =mean_xy_domain(nzm,cloudliq_mask(1:nx,1:ny,1:nzm,thresh_index))
  rain_frac(:,thresh_index) = mean_xy_domain(nzm,rain_mask(1:nx,1:ny,1:nzm,thresh_index))
  cloudice_frac(:,thresh_index) = mean_xy_domain(nzm,cloudice_mask(1:nx,1:ny,1:nzm,thresh_index))
  snow_frac(:,thresh_index) = mean_xy_domain(nzm,snow_mask(1:nx,1:ny,1:nzm,thresh_index))
  graupel_frac(:,thresh_index) = mean_xy_domain(nzm,graupel_mask(1:nx,1:ny,1:nzm,thresh_index))

! Increase the threshold by a factor of 10
frac_threshold = frac_threshold * order_of_magnitude

end do !thresh_index

!Below is used to compute covariances and correlations of microphysical quantities
!=================================================================================
idx_w =1 
idx_microstart = 2

  !Interpolate w to the scalar grid
  do k=1,nz-1
    do i=1,nx
      do j = 1,ny
        w_zt(i,j,k) =  LIN_INT( w(i,j,k+1), w(i,j,k), zi(k+1), zi(k), z(k) ) 
      end do
    end do
  end do

  !Find the domain wide means and variances
  domain_mean_w_zt = mean_xy_domain(nzm, w_zt)
  domain_varnce_w_zt = variance_xy_domain(nzm, w_zt, domain_mean_w_zt)

  do n=idx_microstart,nmicro_fields
  
    domain_mean_micro(n,:) = mean_xy_domain(nzm, micro_field(1:nx,1:ny,1:nzm,n))
    domain_varnce_micro(n,:) = variance_xy_domain(nzm, micro_field(1:nx,1:ny,1:nzm,n),&
                                                domain_mean_micro(n,:))
  end do
 
  !----------------
  ! Within precip means
  !----------------
  mean_micro_ip(incl,:) = mean_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,incl),&
                                          cloudliq_mask(1:nx,1:ny,1:nzm,thresh_out),cloudliq_sum(1:nzm,thresh_out))
  if(doprecip) then
    ! Rain 
    mean_micro_ip(iqr,:) = mean_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,iqr),&
                                          rain_mask(1:nx,1:ny,1:nzm,thresh_out),rain_sum(1:nzm,thresh_out))
  
    mean_micro_ip(inr,:) = mean_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,inr),&
                                          rain_mask(1:nx,1:ny,1:nzm,thresh_out),rain_sum(1:nzm,thresh_out))
    if(doicemicro) then  
      ! Cloud ice
      mean_micro_ip(iqci,:) = mean_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,iqci),&
                                          cloudice_mask(1:nx,1:ny,1:nzm,thresh_out),cloudice_sum(1:nzm,thresh_out))
  
      mean_micro_ip(inci,:) = mean_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,inci),&
                                          cloudice_mask(1:nx,1:ny,1:nzm,thresh_out),cloudice_sum(1:nzm,thresh_out))
      ! Snow 
      mean_micro_ip(iqs,:) = mean_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,iqs),&
                                          snow_mask(1:nx,1:ny,1:nzm,thresh_out),snow_sum(1:nzm,thresh_out))
  
      mean_micro_ip(ins,:) = mean_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,ins),&
                                          snow_mask(1:nx,1:ny,1:nzm,thresh_out),snow_sum(1:nzm,thresh_out))
      if(dograupel) then
      ! Graupel 
        mean_micro_ip(iqg,:) = mean_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,iqg),&
                                          graupel_mask(1:nx,1:ny,1:nzm,thresh_out),graupel_sum(1:nzm,thresh_out))
  
        mean_micro_ip(ing,:) = mean_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,ing),&
                                          graupel_mask(1:nx,1:ny,1:nzm,thresh_out),graupel_sum(1:nzm,thresh_out))
      endif !dograupel
    endif !doice
  endif !doprecip
  
  !----------------
  ! Within precip variances
  !----------------
  varnce_micro_ip(incl,:) = variance_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,incl),&
                                          cloudliq_mask(1:nx,1:ny,1:nzm,thresh_out),mean_micro_ip(incl,:),&
                                          cloudliq_sum(1:nzm,thresh_out))
  if(doprecip) then
    ! Rain 
    varnce_micro_ip(iqr,:) = variance_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,iqr),&
                                          rain_mask(1:nx,1:ny,1:nzm,thresh_out),mean_micro_ip(iqr,:),&
                                          rain_sum(1:nzm,thresh_out))
  
    varnce_micro_ip(inr,:) = variance_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,inr),&
                                          rain_mask(1:nx,1:ny,1:nzm,thresh_out),mean_micro_ip(inr,:),&
                                          rain_sum(1:nzm,thresh_out))
    if(doicemicro) then  
      ! Cloud ice
      varnce_micro_ip(iqci,:) = variance_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,iqci),&
                                          cloudice_mask(1:nx,1:ny,1:nzm,thresh_out),mean_micro_ip(iqci,:),&
                                          cloudice_sum(1:nzm,thresh_out))
  
      varnce_micro_ip(inci,:) = variance_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,inci),&
                                          cloudice_mask(1:nx,1:ny,1:nzm,thresh_out),mean_micro_ip(inci,:),&
                                          cloudice_sum(1:nzm,thresh_out))
      ! Snow 
      varnce_micro_ip(iqs,:) = variance_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,iqs),&
                                          snow_mask(1:nx,1:ny,1:nzm,thresh_out),mean_micro_ip(iqs,:),&
                                          snow_sum(1:nzm,thresh_out))
  
      varnce_micro_ip(ins,:) = variance_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,ins),&
                                          snow_mask(1:nx,1:ny,1:nzm,thresh_out),mean_micro_ip(ins,:),&
                                          snow_sum(1:nzm,thresh_out))
      if(dograupel) then
      ! Graupel 
        varnce_micro_ip(iqg,:) = variance_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,iqg),&
                                          graupel_mask(1:nx,1:ny,1:nzm,thresh_out),mean_micro_ip(iqg,:),&
                                          graupel_sum(1:nzm,thresh_out))
  
        varnce_micro_ip(ing,:) = variance_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,ing),&
                                          graupel_mask(1:nx,1:ny,1:nzm,thresh_out),mean_micro_ip(ing,:),&
                                          graupel_sum(1:nzm,thresh_out))
      endif !dograupel
    endif !doice
  endif !doprecip
 
  !Compute covariances and correlations, first for w(m=1)
  do n = idx_microstart,nmicro_fields 
      micro_covarnce(idx_w,n,:) = covariance_xy_domain(nzm, w_zt, micro_field(1:nx,1:ny,1:nzm,n),&
                              domain_mean_w_zt, domain_mean_micro(n,:) )

      micro_correlations(idx_w,n,:) = compute_correlation(nzm, domain_varnce_w_zt,&
                                  domain_varnce_micro(n,:), micro_covarnce(1,n,:) )

  end do

  !Continue to fill the covariance and correlation matricies, starting at m=2, micro index of incl)  
  do m = idx_microstart, nmicro_fields 
    do n = m, nmicro_fields
  
      micro_covarnce(m,n,:) = covariance_xy_domain(nzm, micro_field(1:nx,1:ny,1:nzm,m), &
          micro_field(1:nx,1:ny,1:nzm,n),domain_mean_micro(m,:), domain_mean_micro(n,:) )

      micro_correlations(m,n,:) = compute_correlation(nzm, domain_varnce_micro(m,:),&
                                  domain_varnce_micro(n,:), micro_covarnce(m,n,:) )
    end do
  end do

#endif UWM_STATS

! Microphysical process rates +++mhwnag
   call hbuf_put('PRC', mPRC, factor_xy)
   call hbuf_put('PRA', mPRA, factor_xy)
   call hbuf_put('PSMLT', mPSMLT, factor_xy)
   call hbuf_put('EVPMS', mEVPMS, factor_xy)
   call hbuf_put('PRACS', mPRACS, factor_xy)
   call hbuf_put('EVPMG', mEVPMG, factor_xy)
   call hbuf_put('PRACG', mPRACG, factor_xy)
   call hbuf_put('PRE', mPRE, factor_xy)
   call hbuf_put('PGMLT', mPGMLT, factor_xy)
   call hbuf_put('MNUCCC', mMNUCCC, factor_xy)
   call hbuf_put('PSACWS', mPSACWS, factor_xy)
   call hbuf_put('PSACWI', mPSACWI, factor_xy)
   call hbuf_put('QMULTS', mQMULTS, factor_xy)
   call hbuf_put('QMULTG', mQMULTG, factor_xy)
   call hbuf_put('PSACWG', mPSACWG, factor_xy)
   call hbuf_put('PGSACW', mPGSACW, factor_xy)

   call hbuf_put('PRD', mPRD, factor_xy)
   call hbuf_put('PRCI', mPRCI, factor_xy)
   call hbuf_put('PRAI', mPRAI, factor_xy)
   call hbuf_put('QMULTR', mQMULTR, factor_xy)
   call hbuf_put('QMULTRG', mQMULTRG, factor_xy)
   call hbuf_put('MNUCCD', mMNUCCD, factor_xy)
   call hbuf_put('PRACI', mPRACI, factor_xy)

   call hbuf_put('PRACIS', mPRACIS, factor_xy)
   call hbuf_put('EPRD', mEPRD, factor_xy)
   call hbuf_put('MNUCCR', mMNUCCR, factor_xy)
   call hbuf_put('PIACR', mPIACR, factor_xy)
   call hbuf_put('PIACRS', mPIACRS, factor_xy)
   call hbuf_put('PGRACS', mPGRACS, factor_xy)
   call hbuf_put('PRDS', mPRDS, factor_xy)
   call hbuf_put('EPRDS', mEPRDS, factor_xy)
   call hbuf_put('PSACR', mPSACR, factor_xy)
   call hbuf_put('PRDG', mPRDG, factor_xy)
   call hbuf_put('EPRDG', mEPRDG, factor_xy)
   call hbuf_put('NPRC1', mNPRC1, factor_xy)
   call hbuf_put('NRAGG', mNRAGG, factor_xy)
   call hbuf_put('NPRACG', mNPRACG, factor_xy)
   call hbuf_put('NSUBR', mNSUBR, factor_xy)
   call hbuf_put('NSMLTR', mNSMLTR, factor_xy)
   call hbuf_put('NGMLTR', mNGMLTR, factor_xy)
   call hbuf_put('NPRACS', mNPRACS, factor_xy)
   call hbuf_put('NNUCCR', mNNUCCR, factor_xy)
   call hbuf_put('NIACR', mNIACR, factor_xy)
   call hbuf_put('NIACRS', mNIACRS, factor_xy)
   call hbuf_put('NGRACS', mNGRACS, factor_xy)
   call hbuf_put('NSMLTS', mNSMLTS, factor_xy)
   call hbuf_put('NSAGG', mNSAGG, factor_xy)
   call hbuf_put('NPRCI', mNPRCI, factor_xy)
   call hbuf_put('NSCNG', mNSCNG, factor_xy)
   call hbuf_put('NSUBS', mNSUBS, factor_xy)

   call hbuf_put('PCC', mPCC, factor_xy)
   call hbuf_put('NNUCCC', mNNUCCC, factor_xy)
   call hbuf_put('NPSACWS', mNPSACWS, factor_xy)
   call hbuf_put('NPRA', mNPRA, factor_xy)
   call hbuf_put('NPRC', mNPRC, factor_xy)
   call hbuf_put('NPSACWI', mNPSACWI, factor_xy)
   call hbuf_put('NPSACWG', mNPSACWG, factor_xy)
   call hbuf_put('NPRAI', mNPRAI, factor_xy)
   call hbuf_put('NMULTS', mNMULTS, factor_xy)
   call hbuf_put('NMULTG', mNMULTG, factor_xy)
   call hbuf_put('NMULTR', mNMULTR, factor_xy)
   call hbuf_put('NMULTRG', mNMULTRG, factor_xy)
   call hbuf_put('NNUCCD', mNNUCCD, factor_xy)
   call hbuf_put('NSUBI', mNSUBI, factor_xy)
   call hbuf_put('NGMLTG', mNGMLTG, factor_xy)
   call hbuf_put('NSUBG', mNSUBG, factor_xy)
   call hbuf_put('NACT', mNACT, factor_xy)

   call hbuf_put('SIZEFIX_NR', mSIZEFIX_NR, factor_xy)
   call hbuf_put('SIZEFIX_NC', mSIZEFIX_NC, factor_xy)
   call hbuf_put('SIZEFIX_NI', mSIZEFIX_NI, factor_xy)
   call hbuf_put('SIZEFIX_NS', mSIZEFIX_NS, factor_xy)
   call hbuf_put('SIZEFIX_NG', mSIZEFIX_NG, factor_xy)
   call hbuf_put('NEGFIX_NR', mNEGFIX_NR, factor_xy)
   call hbuf_put('NEGFIX_NC', mNEGFIX_NC, factor_xy)
   call hbuf_put('NEGFIX_NI', mNEGFIX_NI, factor_xy)
   call hbuf_put('NEGFIX_NS', mNEGFIX_NS, factor_xy)
   call hbuf_put('NEGFIX_NG', mNEGFIX_NG, factor_xy)
   call hbuf_put('NIM_MORR_CL', mNIM_MORR_CL, factor_xy)

   call hbuf_put('QC_INST', mQC_INST, factor_xy)
   call hbuf_put('QR_INST', mQR_INST, factor_xy)
   call hbuf_put('QI_INST', mQI_INST, factor_xy)
   call hbuf_put('QS_INST', mQS_INST, factor_xy)
   call hbuf_put('QG_INST', mQG_INST, factor_xy)
   call hbuf_put('NC_INST', mNC_INST, factor_xy)
   call hbuf_put('NR_INST', mNR_INST, factor_xy)
   call hbuf_put('NI_INST', mNI_INST, factor_xy)
   call hbuf_put('NS_INST', mNS_INST, factor_xy)
   call hbuf_put('NG_INST', mNG_INST, factor_xy)



#ifdef UWM_STATS
  call hbuf_put('cloudliq_frac_em8', cloudliq_frac(:,em8), 1.)
  call hbuf_put('cloudliq_frac_em7', cloudliq_frac(:,em7), 1.)
  call hbuf_put('cloudliq_frac_em6', cloudliq_frac(:,em6), 1.)
  call hbuf_put('cloudliq_frac_em5', cloudliq_frac(:,em5), 1.)
  call hbuf_put('cloudliq_frac_em4', cloudliq_frac(:,em4), 1.)
  
  if (doprecip) then
    call hbuf_put('rain_frac_em8', rain_frac(:,em8), 1.)
    call hbuf_put('rain_frac_em7', rain_frac(:,em7), 1.)
    call hbuf_put('rain_frac_em6', rain_frac(:,em6), 1.)
    call hbuf_put('rain_frac_em5', rain_frac(:,em5), 1.)
    call hbuf_put('rain_frac_em4', rain_frac(:,em4), 1.)
    
    if (doicemicro) then
      call hbuf_put('cloudice_frac_em8', cloudice_frac(:,em8), 1.)
      call hbuf_put('cloudice_frac_em7', cloudice_frac(:,em7), 1.)
      call hbuf_put('cloudice_frac_em6', cloudice_frac(:,em6), 1.)
      call hbuf_put('cloudice_frac_em5', cloudice_frac(:,em5), 1.)
      call hbuf_put('cloudice_frac_em4', cloudice_frac(:,em4), 1.)
      
      call hbuf_put('snow_frac_em8', snow_frac(:,em8), 1.)
      call hbuf_put('snow_frac_em7', snow_frac(:,em7), 1.)
      call hbuf_put('snow_frac_em6', snow_frac(:,em6), 1.)
      call hbuf_put('snow_frac_em5', snow_frac(:,em5), 1.)
      call hbuf_put('snow_frac_em4', snow_frac(:,em4), 1.)
      
      if (dograupel) then  
        call hbuf_put('graupel_frac_em8', graupel_frac(:,em8), 1.)
        call hbuf_put('graupel_frac_em7', graupel_frac(:,em7), 1.)
        call hbuf_put('graupel_frac_em6', graupel_frac(:,em6), 1.)
        call hbuf_put('graupel_frac_em5', graupel_frac(:,em5), 1.)
        call hbuf_put('graupel_frac_em4', graupel_frac(:,em4), 1.)
      end if
    end if 
  end if 

!----------------------
! Domain-wide variance
!----------------------
  call hbuf_put('ncloudliqp2', domain_varnce_micro(incl,:), 1.)
  
  if (doprecip) then
    call hbuf_put('nrainp2', domain_varnce_micro(inr,:), 1.)
    call hbuf_put('qrainp2', domain_varnce_micro(iqr,:), 1.)
    
    if (doicemicro) then
      call hbuf_put('ncloudicep2', domain_varnce_micro(inci,:), 1.)
      call hbuf_put('qcloudicep2', domain_varnce_micro(iqci,:), 1.)
      
      call hbuf_put('nsnowp2', domain_varnce_micro(ins,:), 1.)
      call hbuf_put('qsnowp2', domain_varnce_micro(iqs,:), 1.)
      
      if (dograupel) then  
        call hbuf_put('ngraupelp2', domain_varnce_micro(ing,:), 1.)
        call hbuf_put('qgraupelp2', domain_varnce_micro(iqg,:), 1.)
      end if
    end if 
  endif

!----------------------
! Within-precip mean
!----------------------
  call hbuf_put('ncloudliqm_ip', mean_micro_ip(incl,:), 1.)
  
  if (doprecip) then
    call hbuf_put('nrainm_ip', mean_micro_ip(inr,:), 1.)
    call hbuf_put('qrainm_ip', mean_micro_ip(iqr,:), 1.)
    
    if (doicemicro) then
      call hbuf_put('ncloudicem_ip', mean_micro_ip(inci,:), 1.)
      call hbuf_put('qcloudicem_ip', mean_micro_ip(iqci,:), 1.)
      
      call hbuf_put('nsnowm_ip', mean_micro_ip(ins,:), 1.)
      call hbuf_put('qsnowm_ip', mean_micro_ip(iqs,:), 1.)
      
      if (dograupel) then  
        call hbuf_put('ngraupelm_ip', mean_micro_ip(ing,:), 1.)
        call hbuf_put('qgraupelm_ip', mean_micro_ip(iqg,:), 1.)
      end if
    end if 
  endif


!----------------------
! Within-precip variance
!----------------------
  call hbuf_put('ncloudliqp2_ip', varnce_micro_ip(incl,:), 1.)
  
  if (doprecip) then
    call hbuf_put('nrainp2_ip', varnce_micro_ip(inr,:), 1.)
    call hbuf_put('qrainp2_ip', varnce_micro_ip(iqr,:), 1.)
    
    if (doicemicro) then
      call hbuf_put('ncloudicep2_ip', varnce_micro_ip(inci,:), 1.)
      call hbuf_put('qcloudicep2_ip', varnce_micro_ip(iqci,:), 1.)
      
      call hbuf_put('nsnowp2_ip', varnce_micro_ip(ins,:), 1.)
      call hbuf_put('qsnowp2_ip', varnce_micro_ip(iqs,:), 1.)
      
      if (dograupel) then  
        call hbuf_put('ngraupelp2_ip', varnce_micro_ip(ing,:), 1.)
        call hbuf_put('qgraupelp2_ip', varnce_micro_ip(iqg,:), 1.)
      end if
    end if 
  endif

! weberjk(UWM), Microphysical correlations

   do k = 1, nzm, 1
      if(dopredictNc) then
        if ( micro_correlations(idx_w,incl,k) /= undef_corr ) then
           call hbuf_put_level('corr_w_ncl',micro_correlations(idx_w,incl,k),1.,k)
           corravg_count(idx_cor_w_ncl,k) = corravg_count(idx_cor_w_ncl,k) + 1
        else ! micro_correlations(1,incl,k) == undef_corr
           call hbuf_put_level('corr_w_ncl',0.0,1.,k)
        endif ! micro_correlations(1,incl,k) /= undef_corr
      end if

      if (doprecip) then
      if ( micro_correlations(idx_w,iqr,k) /= undef_corr ) then
           call hbuf_put_level('corr_w_qr',micro_correlations(idx_w,iqr,k),1.,k)
           corravg_count(idx_cor_w_qr,k) = corravg_count(idx_cor_w_qr,k) + 1
        else ! micro_correlations(1,iqr,k) == undef_corr
           call hbuf_put_level('corr_w_qr',0.0,1.,k)
        endif ! micro_correlations(1,iqr,k) /= undef_corr
        if ( micro_correlations(idx_w,inr,k) /= undef_corr ) then
           call hbuf_put_level('corr_w_nr',micro_correlations(idx_w,inr,k),1.,k)
           corravg_count(idx_cor_w_nr,k) = corravg_count(idx_cor_w_nr,k) + 1
        else ! micro_correlations(1,inr,k) == undef_corr
           call hbuf_put_level('corr_w_nr',0.0,1.,k)
        endif ! ! micro_correlations(1,inr,k) /= undef_corr
      end if ! doprecip

      if (doicemicro) then
      if ( micro_correlations(idx_w,iqci,k) /= undef_corr ) then
           call hbuf_put_level('corr_w_qci',micro_correlations(idx_w,iqci,k),1.,k)
           corravg_count(idx_cor_w_qci,k) = corravg_count(idx_cor_w_qci,k) + 1
        else ! micro_correlations(1,iqci,k) == undef_corr
           call hbuf_put_level('corr_w_qci',0.0,1.,k)
        endif ! micro_correlations(1,iqci,k) /= undef_corr
        if ( micro_correlations(idx_w,inci,k) /= undef_corr ) then
           call hbuf_put_level('corr_w_nci',micro_correlations(idx_w,inci,k),1.,k)
           corravg_count(idx_cor_w_nci,k) = corravg_count(idx_cor_w_nci,k) + 1
        else ! micro_correlations(1,inci,k) == undef_corr
           call hbuf_put_level('corr_w_nci',0.0,1.,k)
        endif ! micro_correlations(1,inci,k) /= undef_corr
        if ( micro_correlations(idx_w,iqs,k) /= undef_corr ) then
           call hbuf_put_level('corr_w_qs',micro_correlations(idx_w,iqs,k),1.,k)
           corravg_count(idx_cor_w_qs,k) = corravg_count(idx_cor_w_qs,k) + 1
        else ! micro_correlations(1,iqs,k) == undef_corr
           call hbuf_put_level('corr_w_qs',0.0,1.,k)
        endif ! micro_correlations(1,iqs,k) /= undef_corr
        if ( micro_correlations(idx_w,ins,k) /= undef_corr ) then
           call hbuf_put_level('corr_w_ns',micro_correlations(idx_w,ins,k),1.,k)
           corravg_count(idx_cor_w_ns,k) = corravg_count(idx_cor_w_ns,k) + 1
        else ! micro_correlations(1,ins,k) == undef_corr
           call hbuf_put_level('corr_w_ns',0.0,1.,k)
        endif ! micro_correlations(1,ins,k) /= undef_corr
      end if ! doicemicro

      if (dograupel) then
        if ( micro_correlations(idx_w,iqg,k) /= undef_corr ) then
           call hbuf_put_level('corr_w_qg',micro_correlations(idx_w,iqg,k),1.,k)
           corravg_count(idx_cor_w_qg,k) = corravg_count(idx_cor_w_qg,k) + 1
        else ! micro_correlations(1,iqg,k) == undef_corr
           call hbuf_put_level('corr_w_qg',0.0,1.,k)
        endif ! micro_correlations(1,iqg,k) /= undef_corr
        if ( micro_correlations(idx_w,ing,k) /= undef_corr ) then
           call hbuf_put_level('corr_w_ng',micro_correlations(idx_w,ing,k),1.,k)
           corravg_count(idx_cor_w_ng,k) = corravg_count(idx_cor_w_ng,k) + 1
        else ! micro_correlations(1,ing,k) == undef_corr
           call hbuf_put_level('corr_w_ng',0.0,1.,k)
        endif ! micro_correlations(1,ing,k) /= undef_corr
      end if ! dograupel

      if (dopredictNc) then
     
        if (doprecip) then
           if ( micro_correlations(incl,iqr,k) /= undef_corr ) then
              call hbuf_put_level('corr_ncl_qr',micro_correlations(incl,iqr,k),1.,k)
              corravg_count(idx_cor_ncl_qr,k) = corravg_count(idx_cor_ncl_qr,k) + 1
           else ! micro_correlations(incl,iqr,k) == undef_corr
              call hbuf_put_level('corr_ncl_qr',0.0,1.,k)
           endif ! micro_correlations(incl,iqr,k) /= undef_corr
           if ( micro_correlations(incl,inr,k) /= undef_corr ) then
              call hbuf_put_level('corr_ncl_nr',micro_correlations(incl,inr,k),1.,k)
              corravg_count(idx_cor_ncl_nr,k) = corravg_count(idx_cor_ncl_nr,k) + 1
           else ! micro_correlations(incl,inr,k) == undef_corr
              call hbuf_put_level('corr_ncl_nr',0.0,1.,k)
           endif ! micro_correlations(incl,inr,k) /= undef_corr
        end if
     
        if (doicemicro) then
           if ( micro_correlations(incl,iqci,k) /= undef_corr ) then
              call hbuf_put_level('corr_ncl_qci',micro_correlations(incl,iqci,k),1.,k) 
              corravg_count(idx_cor_ncl_qci,k) = corravg_count(idx_cor_ncl_qci,k) + 1
           else ! micro_correlations(incl,iqci,k) == undef_corr
              call hbuf_put_level('corr_ncl_qci',0.0,1.,k)
           endif ! micro_correlations(incl,iqci,k) /= undef_corr
           if ( micro_correlations(incl,inci,k) /= undef_corr ) then
              call hbuf_put_level('corr_ncl_nci',micro_correlations(incl,inci,k),1.,k)
              corravg_count(idx_cor_ncl_nci,k) = corravg_count(idx_cor_ncl_nci,k) + 1
           else ! micro_correlations(incl,inci,k) == undef_corr
              call hbuf_put_level('corr_ncl_nci',0.0,1.,k)
           endif ! micro_correlations(incl,inci,k) /= undef_corr
           if ( micro_correlations(incl,iqs,k) /= undef_corr ) then
              call hbuf_put_level('corr_ncl_qs',micro_correlations(incl,iqs,k),1.,k)
              corravg_count(idx_cor_ncl_qs,k) = corravg_count(idx_cor_ncl_qs,k) + 1
           else ! micro_correlations(incl,iqs,k) == undef_corr
              call hbuf_put_level('corr_ncl_qs',0.0,1.,k)
           endif ! micro_correlations(incl,iqs,k) /= undef_corr
           if ( micro_correlations(incl,ins,k) /= undef_corr ) then
              call hbuf_put_level('corr_ncl_ns',micro_correlations(incl,ins,k),1.,k)
              corravg_count(idx_cor_ncl_ns,k) = corravg_count(idx_cor_ncl_ns,k) + 1
           else ! micro_correlations(incl,ins,k) == undef_corr
              call hbuf_put_level('corr_ncl_ns',0.0,1.,k)
           endif ! micro_correlations(incl,ins,k) /= undef_corr
        end if    

        if (dograupel) then
           if ( micro_correlations(incl,iqg,k) /= undef_corr ) then
              call hbuf_put_level('corr_ncl_qg',micro_correlations(incl,iqg,k),1.,k)
              corravg_count(idx_cor_ncl_qg,k) = corravg_count(idx_cor_ncl_qg,k) + 1
           else ! micro_correlations(incl,iqg,k) == undef_corr
              call hbuf_put_level('corr_ncl_qg',0.0,1.,k)
           endif ! micro_correlations(incl,iqg,k) /= undef_corr
           if ( micro_correlations(incl,ing,k) /= undef_corr ) then
              call hbuf_put_level('corr_ncl_ng',micro_correlations(incl,ing,k),1.,k)
              corravg_count(idx_cor_ncl_ng,k) = corravg_count(idx_cor_ncl_ng,k) + 1
           else ! micro_correlations(incl,ing,k) == undef_corr
              call hbuf_put_level('corr_ncl_ng',0.0,1.,k)
           endif ! micro_correlations(incl,ing,k) /= undef_corr
        end if   

      end if !dopredictNc
  
      if (doprecip) then 
      
      if (dopredictNc) then
        call hbuf_put_level('covarnce_qr_nr',micro_covarnce(iqr,inr,k),1.,k)
      
          if ( micro_correlations(iqr,inr,k) /= undef_corr ) then
             call hbuf_put_level('corr_qr_nr',micro_correlations(iqr,inr,k),1.,k)
             corravg_count(idx_cor_qr_nr,k) = corravg_count(idx_cor_qr_nr,k) + 1
          else ! micro_correlations(iqr,inr,k) == undef_corr
             call hbuf_put_level('corr_qr_nr',0.0,1.,k)
          endif ! micro_correlations(iqr,inr,k) /= undef_corr
      end if

        if (doicemicro) then
           if ( micro_correlations(iqr,iqci,k) /= undef_corr ) then
              call hbuf_put_level('corr_qr_qci',micro_correlations(iqr,iqci,k),1.,k)
              corravg_count(idx_cor_qr_qci,k) = corravg_count(idx_cor_qr_qci,k) + 1
           else ! micro_correlations(iqr,iqci,k) == undef_corr
              call hbuf_put_level('corr_qr_qci',0.0,1.,k)
           endif ! micro_correlations(iqr,iqci,k) /= undef_corr
           if ( micro_correlations(iqr,inci,k) /= undef_corr ) then
              call hbuf_put_level('corr_qr_nci',micro_correlations(iqr,inci,k),1.,k)
              corravg_count(idx_cor_qr_nci,k) = corravg_count(idx_cor_qr_nci,k) + 1
           else ! micro_correlations(iqr,inci,k) == undef_corr
              call hbuf_put_level('corr_qr_nci',0.0,1.,k)
           endif ! micro_correlations(iqr,inci,k) /= undef_corr
           if ( micro_correlations(iqr,iqs,k) /= undef_corr ) then
              call hbuf_put_level('corr_qr_qs',micro_correlations(iqr,iqs,k),1.,k)
              corravg_count(idx_cor_qr_qs,k) = corravg_count(idx_cor_qr_qs,k) + 1
           else ! micro_correlations(iqr,iqs,k) == undef_corr
              call hbuf_put_level('corr_qr_qs',0.0,1.,k)
           endif ! micro_correlations(iqr,iqs,k) /= undef_corr
           if ( micro_correlations(iqr,ins,k) /= undef_corr ) then
              call hbuf_put_level('corr_qr_ns',micro_correlations(iqr,ins,k),1.,k)
              corravg_count(idx_cor_qr_ns,k) = corravg_count(idx_cor_qr_ns,k) + 1
           else ! micro_correlations(iqr,ins,k) == undef_corr
              call hbuf_put_level('corr_qr_ns',0.0,1.,k)
           endif ! micro_correlations(iqr,ins,k) /= undef_corr
        end if    

        if (dograupel) then
           if ( micro_correlations(iqr,iqg,k) /= undef_corr ) then
              call hbuf_put_level('corr_qr_qg',micro_correlations(iqr,iqg,k),1.,k)
              corravg_count(idx_cor_qr_qg,k) = corravg_count(idx_cor_qr_qg,k) + 1
           else ! micro_correlations(iqr,iqg,k) == undef_corr
              call hbuf_put_level('corr_qr_qg',0.0,1.,k)
           endif ! micro_correlations(iqr,iqg,k) /= undef_corr
           if ( micro_correlations(iqr,ing,k) /= undef_corr ) then
              call hbuf_put_level('corr_qr_ng',micro_correlations(iqr,ing,k),1.,k)
              corravg_count(idx_cor_qr_ng,k) = corravg_count(idx_cor_qr_ng,k) + 1
           else ! micro_correlations(iqr,ing,k) == undef_corr
              call hbuf_put_level('corr_qr_ng',0.0,1.,k)
           endif ! micro_correlations(iqr,ing,k) /= undef_corr
        end if

        if (doicemicro) then
           if ( micro_correlations(inr,iqci,k) /= undef_corr ) then
              call hbuf_put_level('corr_nr_qci',micro_correlations(inr,iqci,k),1.,k)
              corravg_count(idx_cor_nr_qci,k) = corravg_count(idx_cor_nr_qci,k) + 1
           else ! micro_correlations(inr,iqci,k) == undef_corr
              call hbuf_put_level('corr_nr_qci',0.0,1.,k)
           endif ! micro_correlations(inr,iqci,k) /= undef_corr
           if ( micro_correlations(inr,inci,k) /= undef_corr ) then
              call hbuf_put_level('corr_nr_nci',micro_correlations(inr,inci,k),1.,k)
              corravg_count(idx_cor_nr_nci,k) = corravg_count(idx_cor_nr_nci,k) + 1
           else ! micro_correlations(inr,inci,k) == undef_corr
              call hbuf_put_level('corr_nr_nci',0.0,1.,k)
           endif ! micro_correlations(inr,inci,k) /= undef_corr
           if ( micro_correlations(inr,iqs,k) /= undef_corr ) then
              call hbuf_put_level('corr_nr_qs',micro_correlations(inr,iqs,k),1.,k)
              corravg_count(idx_cor_nr_qs,k) = corravg_count(idx_cor_nr_qs,k) + 1
           else ! micro_correlations(inr,iqs,k) == undef_corr
              call hbuf_put_level('corr_nr_qs',0.0,1.,k)
           endif ! micro_correlations(inr,iqs,k) /= undef_corr
           if ( micro_correlations(inr,ins,k) /= undef_corr ) then
              call hbuf_put_level('corr_nr_ns',micro_correlations(inr,ins,k),1.,k)
              corravg_count(idx_cor_nr_ns,k) = corravg_count(idx_cor_nr_ns,k) + 1
           else ! micro_correlations(inr,ins,k) == undef_corr
              call hbuf_put_level('corr_nr_ns',0.0,1.,k)
           endif ! micro_correlations(inr,ins,k) /= undef_corr
        end if

        if (dograupel) then
        if ( micro_correlations(inr,iqg,k) /= undef_corr ) then
              call hbuf_put_level('corr_nr_qg',micro_correlations(inr,iqg,k),1.,k)
              corravg_count(idx_cor_nr_qg,k) = corravg_count(idx_cor_nr_qg,k) + 1
           else ! micro_correlations(inr,iqg,k) == undef_corr
              call hbuf_put_level('corr_nr_qg',0.0,1.,k)
           endif ! micro_correlations(inr,iqg,k) /= undef_corr
           if ( micro_correlations(inr,ing,k) /= undef_corr ) then
              call hbuf_put_level('corr_nr_ng',micro_correlations(inr,ing,k),1.,k)
              corravg_count(idx_cor_nr_ng,k) = corravg_count(idx_cor_nr_ng,k) + 1
           else ! micro_correlations(inr,ing,k) == undef_corr
              call hbuf_put_level('corr_nr_ng',0.0,1.,k)
           endif ! micro_correlations(inr,ing,k) /= undef_corr
        end if

      end if !doprecip

      if (doicemicro) then   
         if ( micro_correlations(iqci,inci,k) /= undef_corr ) then
            call hbuf_put_level('corr_qci_nci',micro_correlations(iqci,inci,k),1.,k)
            corravg_count(idx_cor_qci_nci,k) = corravg_count(idx_cor_qci_nci,k) + 1
         else ! micro_correlations(iqci,inci,k) == undef_corr
            call hbuf_put_level('corr_qci_nci',0.0,1.,k)
         endif ! micro_correlations(iqci,inci,k) /= undef_corr
         if ( micro_correlations(iqci,iqs,k) /= undef_corr ) then
            call hbuf_put_level('corr_qci_qs',micro_correlations(iqci,iqs,k),1.,k)
            corravg_count(idx_cor_qci_qs,k) = corravg_count(idx_cor_qci_qs,k) + 1
         else ! micro_correlations(iqci,iqs,k) == undef_corr
            call hbuf_put_level('corr_qci_qs',0.0,1.,k)
         endif ! micro_correlations(iqci,iqs,k) /= undef_corr
         if ( micro_correlations(iqci,ins,k) /= undef_corr ) then
            call hbuf_put_level('corr_qci_ns',micro_correlations(iqci,ins,k),1.,k)
            corravg_count(idx_cor_qci_ns,k) = corravg_count(idx_cor_qci_ns,k) + 1
         else ! micro_correlations(iqci,ins,k) == undef_corr
            call hbuf_put_level('corr_qci_ns',0.0,1.,k)
         endif ! micro_correlations(iqci,ins,k) /= undef_corr
      end if   

      if (dograupel) then
         if ( micro_correlations(iqci,iqg,k) /= undef_corr ) then
            call hbuf_put_level('corr_qci_qg',micro_correlations(iqci,iqg,k),1.,k)
            corravg_count(idx_cor_qci_qg,k) = corravg_count(idx_cor_qci_qg,k) + 1
         else ! micro_correlations(iqci,iqg,k) == undef_corr
            call hbuf_put_level('corr_qci_qg',0.0,1.,k)
         endif ! micro_correlations(iqci,iqg,k) /= undef_corr
         if ( micro_correlations(iqci,ing,k) /= undef_corr ) then
            call hbuf_put_level('corr_qci_ng',micro_correlations(iqci,ing,k),1.,k)
            corravg_count(idx_cor_qci_ng,k) = corravg_count(idx_cor_qci_ng,k) + 1
         else ! micro_correlations(iqci,ing,k) == undef_corr
            call hbuf_put_level('corr_qci_ng',0.0,1.,k)
         endif ! micro_correlations(iqci,ing,k) /= undef_corr
      end if
 
      if (doicemicro) then   
         if ( micro_correlations(inci,iqs,k) /= undef_corr ) then
            call hbuf_put_level('corr_nci_qs',micro_correlations(inci,iqs,k),1.,k)
            corravg_count(idx_cor_nci_qs,k) = corravg_count(idx_cor_nci_qs,k) + 1
         else ! micro_correlations(inci,iqs,k) == undef_corr
            call hbuf_put_level('corr_nci_qs',0.0,1.,k)
         endif ! micro_correlations(inci,iqs,k) /= undef_corr
         if ( micro_correlations(inci,ins,k) /= undef_corr ) then
            call hbuf_put_level('corr_nci_ns',micro_correlations(inci,ins,k),1.,k)
            corravg_count(idx_cor_nci_ns,k) = corravg_count(idx_cor_nci_ns,k) + 1
         else ! micro_correlations(inci,ins,k) == undef_corr
            call hbuf_put_level('corr_nci_ns',0.0,1.,k)
         endif ! micro_correlations(inci,ins,k) /= undef_corr
      end if

      if (dograupel) then     
         if ( micro_correlations(inci,iqg,k) /= undef_corr ) then
            call hbuf_put_level('corr_nci_qg',micro_correlations(inci,iqg,k),1.,k)
            corravg_count(idx_cor_nci_qg,k) = corravg_count(idx_cor_nci_qg,k) + 1
         else ! micro_correlations(inci,iqg,k) == undef_corr
            call hbuf_put_level('corr_nci_qg',0.0,1.,k)
         endif ! micro_correlations(inci,iqg,k) /= undef_corr
         if ( micro_correlations(inci,ing,k) /= undef_corr ) then
            call hbuf_put_level('corr_nci_ng',micro_correlations(inci,ing,k),1.,k)
            corravg_count(idx_cor_nci_ng,k) = corravg_count(idx_cor_nci_ng,k) + 1
         else ! micro_correlations(inci,ing,k) == undef_corr
            call hbuf_put_level('corr_nci_ng',0.0,1.,k)
         endif ! micro_correlations(inci,ing,k) /= undef_corr
      end if 

      if (doicemicro) then
         if ( micro_correlations(iqs,ins,k) /= undef_corr ) then
            call hbuf_put_level('corr_qs_ns',micro_correlations(iqs,ins,k),1.,k)
            corravg_count(idx_cor_qs_ns,k) = corravg_count(idx_cor_qs_ns,k) + 1
         else ! micro_correlations(iqs,ins,k) == undef_corr
            call hbuf_put_level('corr_qs_ns',0.0,1.,k)
         endif ! micro_correlations(iqs,ins,k) /= undef_corr
       
          if (dograupel) then   
             if ( micro_correlations(iqs,iqg,k) /= undef_corr ) then
                call hbuf_put_level('corr_qs_qg',micro_correlations(iqs,iqg,k),1.,k)
                corravg_count(idx_cor_qs_qg,k) = corravg_count(idx_cor_qs_qg,k) + 1
             else ! micro_correlations(iqs,iqg,k) == undef_corr
                call hbuf_put_level('corr_qs_qg',0.0,1.,k)
             endif ! micro_correlations(iqs,iqg,k) /= undef_corr
             if ( micro_correlations(iqs,ing,k) /= undef_corr ) then
                call hbuf_put_level('corr_qs_ng',micro_correlations(iqs,ing,k),1.,k)
                corravg_count(idx_cor_qs_ng,k) = corravg_count(idx_cor_qs_ng,k) + 1
             else ! micro_correlations(iqs,ing,k) == undef_corr
                call hbuf_put_level('corr_qs_ng',0.0,1.,k)
             endif ! micro_correlations(iqs,ing,k) /= undef_corr
             if ( micro_correlations(ins,iqg,k) /= undef_corr ) then
                call hbuf_put_level('corr_ns_qg',micro_correlations(ins,iqg,k),1.,k)
                corravg_count(idx_cor_ns_qg,k) = corravg_count(idx_cor_ns_qg,k) + 1
             else ! micro_correlations(ins,iqg,k) == undef_corr
                call hbuf_put_level('corr_ns_qg',0.0,1.,k)
             endif ! micro_correlations(ins,iqg,k) /= undef_corr
             if ( micro_correlations(ins,ing,k) /= undef_corr ) then
                call hbuf_put_level('corr_ns_ng',micro_correlations(ins,ing,k),1.,k)
                corravg_count(idx_cor_ns_ng,k) = corravg_count(idx_cor_ns_ng,k) + 1
             else ! micro_correlations(ins,ing,k) == undef_corr
                call hbuf_put_level('corr_ns_ng',0.0,1.,k)
             endif ! micro_correlations(ins,ing,k) /= undef_corr
          end if
      end if

      if (dograupel) then      
        if ( micro_correlations(iqg,ing,k) /= undef_corr ) then
            call hbuf_put_level('corr_qg_ng',micro_correlations(iqg,ing,k),1.,k)
            corravg_count(idx_cor_qg_ng,k) = corravg_count(idx_cor_qg_ng,k) + 1
         else ! micro_correlations(iqg,ing,k) == undef_corr
            call hbuf_put_level('corr_qg_ng',0.0,1.,k)
         endif ! micro_correlations(iqg,ing,k) /= undef_corr
      endif

   enddo ! k = 1, nzm, 1
      
!---------------------------------
! Covariances
!--------------------------------

if(dopredictNc) then
  call hbuf_put('covarnce_w_ncl',micro_covarnce(idx_w,incl,:),1.)
endif !dopredictNc

if (doprecip) then
  call hbuf_put('covarnce_w_qr',micro_covarnce(idx_w,iqr,:),1.)
  call hbuf_put('covarnce_w_nr',micro_covarnce(idx_w,inr,:),1.)
end if ! doprecip

if (doicemicro) then
  call hbuf_put('covarnce_w_qci',micro_covarnce(idx_w,iqci,:),1.)
  call hbuf_put('covarnce_w_nci',micro_covarnce(idx_w,inci,:),1.)
  call hbuf_put('covarnce_w_qs',micro_covarnce(idx_w,iqs,:),1.)
  call hbuf_put('covarnce_w_ns',micro_covarnce(idx_w,ins,:),1.)
end if ! doicemicro

if (dograupel) then
  call hbuf_put('covarnce_w_qg',micro_covarnce(idx_w,iqg,:),1.)
  call hbuf_put('covarnce_w_ng',micro_covarnce(idx_w,ing,:),1.)
end if ! dograupel

if(dopredictNc) then  
  if (doprecip) then
    call hbuf_put('covarnce_ncl_qr',micro_covarnce(incl,iqr,:),1.)
    call hbuf_put('covarnce_ncl_nr',micro_covarnce(incl,inr,:),1.)
  end if ! doprecip
     
  if (doicemicro) then
    call hbuf_put('covarnce_ncl_qci',micro_covarnce(incl,iqci,:),1.) 
    call hbuf_put('covarnce_ncl_nci',micro_covarnce(incl,inci,:),1.)
    call hbuf_put('covarnce_ncl_qs',micro_covarnce(incl,iqs,:),1.)
    call hbuf_put('covarnce_ncl_ns',micro_covarnce(incl,ins,:),1.)
  end if ! doicemicro    

  if (dograupel) then
    call hbuf_put('covarnce_ncl_qg',micro_covarnce(incl,iqg,:),1.)
    call hbuf_put('covarnce_ncl_ng',micro_covarnce(incl,ing,:),1.)
  end if ! dograupel

end if !dopredictNc
  
if (doprecip) then 
  if (dopredictNc) then      
  call hbuf_put('covarnce_qr_nr',micro_covarnce(iqr,inr,:),1.)
  endif

  if (doicemicro) then
    call hbuf_put('covarnce_qr_qci',micro_covarnce(iqr,iqci,:),1.)
    call hbuf_put('covarnce_qr_nci',micro_covarnce(iqr,inci,:),1.)
    call hbuf_put('covarnce_qr_qs',micro_covarnce(iqr,iqs,:),1.)
    call hbuf_put('covarnce_qr_ns',micro_covarnce(iqr,ins,:),1.)

    if (dograupel) then
      call hbuf_put('covarnce_qr_qg',micro_covarnce(iqr,iqg,:),1.)
      call hbuf_put('covarnce_qr_ng',micro_covarnce(iqr,ing,:),1.)
    end if !dograupel

    call hbuf_put('covarnce_nr_qci',micro_covarnce(inr,iqci,:),1.)
    call hbuf_put('covarnce_nr_nci',micro_covarnce(inr,inci,:),1.)
    call hbuf_put('covarnce_nr_qs',micro_covarnce(inr,iqs,:),1.)
    call hbuf_put('covarnce_nr_ns',micro_covarnce(inr,ins,:),1.)

    if (dograupel) then
      call hbuf_put('covarnce_nr_qg',micro_covarnce(inr,iqg,:),1.)
      call hbuf_put('covarnce_nr_ng',micro_covarnce(inr,ing,:),1.)
    end if !dograupel
  
  end if !doicemicro

end if !doprecip

if (doicemicro) then   
  call hbuf_put('covarnce_qci_nci',micro_covarnce(iqci,inci,:),1.)
  call hbuf_put('covarnce_qci_qs',micro_covarnce(iqci,iqs,:),1.)
  call hbuf_put('covarnce_qci_ns',micro_covarnce(iqci,ins,:),1.)

  if (dograupel) then
    call hbuf_put('covarnce_qci_qg',micro_covarnce(iqci,iqg,:),1.)
    call hbuf_put('covarnce_qci_ng',micro_covarnce(iqci,ing,:),1.)
  end if !dograupel
 
  call hbuf_put('covarnce_nci_qs',micro_covarnce(inci,iqs,:),1.)
  call hbuf_put('covarnce_nci_ns',micro_covarnce(inci,ins,:),1.)
      
  if (dograupel) then     
    call hbuf_put('covarnce_nci_qg',micro_covarnce(inci,iqg,:),1.)
    call hbuf_put('covarnce_nci_ng',micro_covarnce(inci,ing,:),1.)
  end if !dograupel

  call hbuf_put('covarnce_qs_ns',micro_covarnce(iqs,ins,:),1.)
       
  if (dograupel) then   
    call hbuf_put('covarnce_qs_qg',micro_covarnce(iqs,iqg,:),1.)
    call hbuf_put('covarnce_qs_ng',micro_covarnce(iqs,ing,:),1.)
    call hbuf_put('covarnce_ns_qg',micro_covarnce(ins,iqg,:),1.)
    call hbuf_put('covarnce_ns_ng',micro_covarnce(ins,ing,:),1.)
    call hbuf_put('covarnce_qg_ng',micro_covarnce(iqg,ing,:),1.)
  endif !dograupel
endif !doicemicro

#endif /*UWM_STATS*/

call t_stopf ('micro_statistics')

end subroutine micro_statistics

!-----------------------------------------
subroutine satadj_liquid(nzm,tabs,qt,qc,pres)
  !bloss/qt: Utility routine based on cloud.f90 in 
  !  MICRO_SAM1MOM that was written by Marat Khairoutdinov.
  !  This routine performs a saturation adjustment for
  !  cloud liquid water only using a Newton method.
  !  While 20 iterations are allowed, most often this
  !  routine should exit in five iterations or less.
  !  Only a single calculation of the saturation vapor
  !  pressure is required in subsaturated air.

  use module_mp_GRAUPEL, only: polysvp
  use params, only: cp, lcond, rv, fac_cond
  implicit none

  integer, intent(in) :: nzm
  real, intent(inout), dimension(nzm) :: tabs ! absolute temperature, K
  real, intent(inout), dimension(nzm) :: qt  ! on input: qt; on output: qv
  real, intent(out), dimension(nzm) :: qc ! cloud liquid water, kg/kg
  real, intent(in), dimension(nzm) :: pres ! pressure, Pa

  real tabs1, dtabs, thresh, esat1, qsat1, fff, dfff
  integer k, niter

  integer, parameter :: maxiter = 20

  !bloss/qt: quick saturation adjustment to compute cloud liquid water content.
  do k = 1,nzm
    tabs1 = tabs(k) 
    esat1 = polysvp(tabs1,0)
    qsat1 = 0.622*esat1/ (pres(k) - esat1)
    qc(k) = 0. ! no cloud unless qt > qsat
    
    if (qt(k).gt.qsat1) then

      ! if unsaturated, nothing to do (i.e., qv=qt, T=Tl) --> just exit.
      ! if saturated, do saturation adjustment 
      !    (modeled after Marat's cloud.f90).

      ! generate initial guess based on above calculation of qsat
      dtabs = + fac_cond*MAX(0.,qt(k) - qsat1) &
           / ( 1. + lcond**2*qsat1/(cp*rv*tabs1**2) )
      tabs1 = tabs1 + dtabs
      niter = 1

      ! convergence threshold: min of 0.01K and latent heating due to
      !    condensation of 1% of saturation mixing ratio.
      thresh = MIN(0.01, 0.01*fac_cond*qsat1)

      ! iterate while temperature increment > thresh and niter < maxiter
      do while((ABS(dtabs).GT.thresh) .AND. (niter.lt.maxiter))

        esat1 = polysvp(tabs1,0)
        qsat1 = 0.622*esat1/ (pres(k) - esat1) ! saturation mixing ratio

        fff = tabs(k) - tabs1 + fac_cond*MAX(0.,qt(k) - qsat1)
        dfff = 1. + lcond**2*qsat1/(cp*rv*tabs1**2)
        dtabs = fff/dfff
        tabs1 = tabs1 + dtabs

        niter = niter + 1

      end do

      qc(k) = MAX( 0.,tabs1 - tabs(k) )/fac_cond ! cloud liquid mass mixing ratio
      qt(k) = qt(k) - qc(k) ! This now holds the water vapor mass mixing ratio.
      tabs(k) = tabs1 ! update temperature.
      
      if(niter.gt.maxiter-1) write(*,*) 'Reached iteration limit in satadj_liquid'

    end if ! qt_in > qsat

  end do ! k = 1,nzm

end subroutine satadj_liquid

!-----------------------------------------------------------------------
! Supply function that computes total water in a domain:
!
real(8) function total_water()

  use vars, only : nstep,nprint,adz,dz,rho
  real(8) tmp
  integer i,j,k,m

  total_water = 0.
  do m=1,nmicro_fields
   if(flag_wmass(m).eq.1) then
    do k=1,nzm
      tmp = 0.
      do j=1,ny
        do i=1,nx
          tmp = tmp + micro_field(i,j,k,m)
        end do
      end do
      total_water = total_water + tmp*adz(k)*dz*rho(k)
    end do
   end if
  end do

end function total_water

function Get_reffc() ! liquid water
  real, dimension(nx,ny,nzm) :: Get_reffc
  Get_reffc = reffc
end function Get_reffc

function Get_reffi() ! ice
  real, dimension(nx,ny,nzm) :: Get_reffi
  Get_reffi = reffi
end function Get_reffi

#if defined(CLUBB) || defined(UWM_STATS)

!-------------------------------------------------------------------------------
ELEMENTAL FUNCTION LIN_INT( var_high, var_low, height_high, height_low, height_int )

! This function computes a linear interpolation of the value of variable.
! Given two known values of a variable at two height values, the value
! of that variable at a height between those two height levels (rather 
! than a height outside of those two height levels) is computed.
!
! Here is a diagram:
!
!  ################################ Height high, know variable value
!
!
!
!  -------------------------------- Height to be interpolated to; linear interpolation
!
!
!
!
!
!  ################################ Height low, know variable value
!
!
! FORMULA:
!
! variable(@ Height interpolation) =
!
! [ (variable(@ Height high) - variable(@ Height low)) / (Height high - Height low) ]
! * (Height interpolation - Height low)  +  variable(@ Height low)

! Author: Brian Griffin, UW-Milwaukee
! Modifications: Dave Schanen added the elemental attribute 4 Nov 2008
! References: None

IMPLICIT NONE

! Input Variables
REAL, INTENT(IN):: var_high
REAL, INTENT(IN):: var_low
REAL, INTENT(IN):: height_high
REAL, INTENT(IN):: height_low
REAL, INTENT(IN):: height_int

! Output Variable
REAL:: LIN_INT

LIN_INT = ( var_high - var_low ) / ( height_high - height_low ) &
         * ( height_int - height_low ) + var_low


END FUNCTION LIN_INT
#endif

end module microphysics



