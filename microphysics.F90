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
use params, only: doprecip, docloud
#ifdef CLUBB
use sgs_params, only: doclubb
#endif

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

#ifdef CLUBB
real, allocatable, target, dimension(:,:,:,:) :: micro_field  ! holds mphys quantities
#else
real, allocatable, dimension(:,:,:,:) :: micro_field  ! holds mphys quantities
#endif

! indices of water quantities in micro_field, e.g. qv = micro_field(:,:,:,iqv)
integer :: iqv, iqcl, iqci, iqr, iqs, iqg, incl, inci, inr, ins, ing
integer :: index_water_vapor ! separate water vapor index used by SAM

#ifdef PNNL_STATS
integer :: indx_qr, indx_qi, indx_qs, indx_qg ! For advection theta-l/qtog stats
#endif /*PNNL_STATS*/

real, allocatable, dimension(:) :: lfac
integer, allocatable, dimension(:) :: flag_wmass, flag_precip, flag_number
integer, allocatable, dimension(:) :: flag_micro3Dout

integer, parameter :: index_cloud_ice = -1 ! historical variable (don't change)

real, allocatable, dimension(:,:,:) :: fluxbmk, fluxtmk !surface/top fluxes
real, allocatable, dimension(:,:,:) :: reffc, reffi
#ifdef CLUBB
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
#ifdef CLUBB 
real, dimension(:,:,:), pointer :: conc, qn
#endif

#ifdef UWM_STATS
integer ::                                                                                         &
  idx_cor_chi_w, idx_cor_chi_Nc, idx_cor_chi_rr, idx_cor_chi_Nr, idx_cor_chi_ri, idx_cor_chi_Ni,   &
  idx_cor_chi_rs, idx_cor_chi_Ns, idx_cor_chi_rg, idx_cor_chi_Ng,                                  &
  idx_cor_w_Nc, idx_cor_w_rr, idx_cor_w_Nr, idx_cor_w_ri, idx_cor_w_Ni, idx_cor_w_rs, idx_cor_w_Ns,&
  idx_cor_w_rg, idx_cor_w_Ng, idx_cor_Nc_rr, idx_cor_Nc_Nr, idx_cor_Nc_ri, idx_cor_Nc_Ni,          &
  idx_cor_Nc_rs, idx_cor_Nc_Ns, idx_cor_Nc_rg, idx_cor_Nc_Ng, idx_cor_rr_Nr, idx_cor_rr_ri,        &
  idx_cor_rr_Ni, idx_cor_rr_rs, idx_cor_rr_Ns, idx_cor_rr_rg, idx_cor_rr_Ng, idx_cor_Nr_ri,        &
  idx_cor_Nr_Ni, idx_cor_Nr_rs, idx_cor_Nr_Ns, idx_cor_Nr_rg, idx_cor_Nr_Ng, idx_cor_ri_Ni,        &
  idx_cor_ri_rs, idx_cor_ri_Ns, idx_cor_ri_rg, idx_cor_ri_Ng, idx_cor_Ni_rs, idx_cor_Ni_Ns,        &
  idx_cor_Ni_rg, idx_cor_Ni_Ng, idx_cor_rs_Ns, idx_cor_rs_rg, idx_cor_rs_Ng, idx_cor_Ns_rg,        &
  idx_cor_Ns_Ng, idx_cor_rg_Ng,                                                                    &
  
  idx_ic_cor_chi_w, idx_ic_cor_chi_Nc, idx_ic_cor_chi_rr, idx_ic_cor_chi_Nr, idx_ic_cor_chi_ri,    &
  idx_ic_cor_chi_Ni, idx_ic_cor_chi_rs, idx_ic_cor_chi_Ns, idx_ic_cor_chi_rg, idx_ic_cor_chi_Ng,   &
  idx_ic_cor_w_Nc, idx_ic_cor_w_rr, idx_ic_cor_w_Nr, idx_ic_cor_w_ri, idx_ic_cor_w_Ni,             &
  idx_ic_cor_w_rs, idx_ic_cor_w_Ns, idx_ic_cor_w_rg, idx_ic_cor_w_Ng, idx_ic_cor_Nc_rr,            &
  idx_ic_cor_Nc_Nr, idx_ic_cor_Nc_ri, idx_ic_cor_Nc_Ni, idx_ic_cor_Nc_rs, idx_ic_cor_Nc_Ns,        &
  idx_ic_cor_Nc_rg, idx_ic_cor_Nc_Ng, idx_ic_cor_rr_Nr, idx_ic_cor_rr_ri, idx_ic_cor_rr_Ni,        &
  idx_ic_cor_rr_rs, idx_ic_cor_rr_Ns, idx_ic_cor_rr_rg, idx_ic_cor_rr_Ng, idx_ic_cor_Nr_ri,        &
  idx_ic_cor_Nr_Ni, idx_ic_cor_Nr_rs, idx_ic_cor_Nr_Ns, idx_ic_cor_Nr_rg, idx_ic_cor_Nr_Ng,        &
  idx_ic_cor_ri_Ni, idx_ic_cor_ri_rs, idx_ic_cor_ri_Ns, idx_ic_cor_ri_rg, idx_ic_cor_ri_Ng,        &
  idx_ic_cor_Ni_rs, idx_ic_cor_Ni_Ns, idx_ic_cor_Ni_rg, idx_ic_cor_Ni_Ng, idx_ic_cor_rs_Ns,        &
  idx_ic_cor_rs_rg, idx_ic_cor_rs_Ng, idx_ic_cor_Ns_rg, idx_ic_cor_Ns_Ng, idx_ic_cor_rg_Ng,        &

  idx_oc_cor_chi_w, idx_oc_cor_chi_Nc, idx_oc_cor_chi_rr, idx_oc_cor_chi_Nr, idx_oc_cor_chi_ri,    &
  idx_oc_cor_chi_Ni, idx_oc_cor_chi_rs, idx_oc_cor_chi_Ns, idx_oc_cor_chi_rg, idx_oc_cor_chi_Ng,   &
  idx_oc_cor_w_Nc, idx_oc_cor_w_rr, idx_oc_cor_w_Nr, idx_oc_cor_w_ri, idx_oc_cor_w_Ni,             &
  idx_oc_cor_w_rs, idx_oc_cor_w_Ns, idx_oc_cor_w_rg, idx_oc_cor_w_Ng, idx_oc_cor_Nc_rr,            &
  idx_oc_cor_Nc_Nr, idx_oc_cor_Nc_ri, idx_oc_cor_Nc_Ni, idx_oc_cor_Nc_rs, idx_oc_cor_Nc_Ns,        &
  idx_oc_cor_Nc_rg, idx_oc_cor_Nc_Ng, idx_oc_cor_rr_Nr, idx_oc_cor_rr_ri, idx_oc_cor_rr_Ni,        &
  idx_oc_cor_rr_rs, idx_oc_cor_rr_Ns, idx_oc_cor_rr_rg, idx_oc_cor_rr_Ng, idx_oc_cor_Nr_ri,        &
  idx_oc_cor_Nr_Ni, idx_oc_cor_Nr_rs, idx_oc_cor_Nr_Ns, idx_oc_cor_Nr_rg, idx_oc_cor_Nr_Ng,        &
  idx_oc_cor_ri_Ni, idx_oc_cor_ri_rs, idx_oc_cor_ri_Ns, idx_oc_cor_ri_rg, idx_oc_cor_ri_Ng,        &
  idx_oc_cor_Ni_rs, idx_oc_cor_Ni_Ns, idx_oc_cor_Ni_rg, idx_oc_cor_Ni_Ng, idx_oc_cor_rs_Ns,        &
  idx_oc_cor_rs_rg, idx_oc_cor_rs_Ng, idx_oc_cor_Ns_rg, idx_oc_cor_Ns_Ng, idx_oc_cor_rg_Ng,        &
  
  !Fractions
  idx_cld_frac, idx_rain_frac, idx_cldic_frac, idx_snw_frac, idx_grpl_frac 

! Use SAM's default method for updating liquid/ice static energy
logical :: l_update_mse_using_state_vars = .false.

! Default values for the ensemble fractions. 
logical :: doensemble_fractions = .false.  
real :: frac_threshold_init = 1e-4
integer :: nfractions = 1
integer ::  nfrac_fields = 1

character*10, allocatable, dimension(:) :: frac_name
character*80, allocatable, dimension(:) :: frac_longname
character*10, allocatable, dimension(:) :: frac_in_char
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

#ifdef UWM_STATS
! weberjk(UWM), 3D QR budget
real, allocatable, dimension(:,:,:) :: mPRE_3D, mPRA_3D, mPRC_3D, mPRACS_3D, mMNUCCR_3D, &
                        mQMULTR_3D, mQMULTRG_3D, mPIACR_3D, mPIACRS_3D, mPRACG_3D, mPGRACS_3D, &
                        mPSMLT_3D, mPGMLT_3D, mQR_INST_3D

! weberjk(UWM), 3D NR budget
real, allocatable, dimension(:,:,:) :: mSIZEFIX_NR_3D, mNEGFIX_NR_3D, mNSUBR_3D, mNSMLTR_3D, &
                        mNGMLTR_3D, mNPRC1_3D, mNPRACS_3D, mNNUCCR_3D, mNRAGG_3D, mNIACR_3D, &
                        mNIACRS_3D, mNPRACG_3D, mNGRACS_3D, mNR_INST_3D

! weberjk(UWM), Ancilliary micro variables
real, allocatable, dimension(:,:,:) :: rain_vel_3D, EFFR_3D

! weberjk(UWM), liquid water potential temperature, chi (s_mellor),
! and eta (t_mellor)
real, allocatable, dimension(:,:,:) :: theta_l, chi, eta

! Store the residual between updating 't' using sedimentation or state
! variables
real, allocatable, dimension(:,:,:) :: t_res 

#endif /*UWM_STATS*/


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
#ifdef UWM_STATS
      doensemble_fractions,&  ! Rather than a fixed threshold to compute
                              ! in-'precip' means, variances, etc, use an ensemble of thresholds
      frac_threshold_init, &  ! Initial threshold value
      nfractions         , &  ! Number of fractions to compute
      l_update_mse_using_state_vars, & ! Update liquid/ice static energy using
                                       ! state variables (like CLUBB). By default, SAM uses sedimentation instead.
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
  end if
#endif /* MICRO_RESTART */
#ifdef PNNL_STATS
if(doprecip) then
  indx_qr = iqr  ! Index for rain water mixing ratio (theta-l/qtog adv. budgets)
endif

if(doicemicro) then
  indx_qi = iqci ! Index for ice mixing ratio (theta-l/qtog adv. budgets)
  indx_qs = iqs  ! Index for snow mixing ratio (theta-l/qtog adv. budgets)
endif

if(dograupel) then
  indx_qg = iqg  ! Index for graupel mixing ratio (theta-l/qtog adv. budgets)
endif
#endif /* PNNL_STATS */
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

#ifdef UWM_STATS
! weberjk(UWM), 3D QR budget
      allocate(mPRE_3D(nx,ny,nzm), mPRA_3D(nx,ny,nzm), mPRC_3D(nx,ny,nzm), mPRACS_3D(nx,ny,nzm), mMNUCCR_3D(nx,ny,nzm), &
                        mQMULTR_3D(nx,ny,nzm), mQMULTRG_3D(nx,ny,nzm), mPIACR_3D(nx,ny,nzm), mPIACRS_3D(nx,ny,nzm), &
                        mPRACG_3D(nx,ny,nzm), mPGRACS_3D(nx,ny,nzm), mPSMLT_3D(nx,ny,nzm), mPGMLT_3D(nx,ny,nzm), &
                        mQR_INST_3D(nx,ny,nzm), STAT=ierr)

! weberjk(UWM), 3D NR budget
      allocate(mSIZEFIX_NR_3D(nx,ny,nzm), mNEGFIX_NR_3D(nx,ny,nzm), mNSUBR_3D(nx,ny,nzm), mNSMLTR_3D(nx,ny,nzm), &
                        mNGMLTR_3D(nx,ny,nzm), mNPRC1_3D(nx,ny,nzm), mNPRACS_3D(nx,ny,nzm), mNNUCCR_3D(nx,ny,nzm),&
                        mNRAGG_3D(nx,ny,nzm), mNIACR_3D(nx,ny,nzm), mNIACRS_3D(nx,ny,nzm), mNPRACG_3D(nx,ny,nzm),&
                        mNGRACS_3D(nx,ny,nzm), mNR_INST_3D(nx,ny,nzm), STAT=ierr)

! weberjk(UWM), 3D NR budget
      allocate(rain_vel_3D(nx,ny,nzm), EFFR_3D(nx,ny,nzm), STAT=ierr)

! weberjk(UWM), liquid water potential temperature, chi (s_mellor),
! and eta (t_mellor).
      allocate(theta_l(nx,ny,nzm), chi(nx,ny,nzm), eta(nx,ny,nzm), STAT=ierr)

! Store the residual between calculating liquid/ice static energy using
! sedimentation versus state variables     
      allocate(t_res(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm), STAT=ierr)

! weberjk(UWM), store fraction names and units to create them within SAM. Create
! 5 elements for the (possibly) 5 prognostic variables (cloud, rain, ice, snow,
! graupel)
      allocate(frac_name(5), frac_longname(5), frac_in_char(5), STAT=ierr)
#endif /*UWM_STATS*/
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
#ifdef CLUBB
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
  
  real, dimension(nzm) :: qc0 !, qi0  !MWSWong:qi0 not used, comment out

! Commented out by dschanen UWM 23 Nov 2009 to avoid a linking error
! real, external :: satadj_water 
  integer :: k
  integer :: i, j, n

#ifndef MICRO_RESTART
  write(6,*) 'microphysics: '
  write(6,*) 'microphysics: dopredictNc', dopredictNc
  write(6,*) 'microphysics: doicemicro', doicemicro
  write(6,*) 'microphysics: doprecip', doprecip
  write(6,*) 'microphysics: dograupel', dograupel
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
     if(dopredictNc) then
        mkname(incl) = 'NC'
        mklongname(incl) = 'CLOUD WATER NUMBER CONCENTRATION'
        mkunits(incl) = '#/cm3'
        mkoutputscale(incl) = 1.e-6
     end if
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
         !mstor(k, n) = SUM(micro_field(1:nx,1:ny,k,n))
         mstor(k, n) = SUM(dble(micro_field(1:nx,1:ny,k,n))) !MWSWong: cast to double
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
use sgs_params, only: doclubb, doclubb_sfc_fluxes
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
use vars, only: t2leprec, q2leprec, qwleprec, twleprec, prespot, qsatw
use grid, only: nsave3D, nsave3dstart, nsave3dend
use module_mp_GRAUPEL, only: Nc0
use compute_chi_module, only: compute_chi_eta
use calc_vars_util, only: t2thetal
#endif /*UWM_STATS*/

#ifdef PNNL_STATS
!MWSWong: budget for THL, QTOG, and QTHL
use vars, only: thel2leprec, thelwleprec, thel, thellat, &
                qtog2leprec, qtogwleprec, qtog, &
                qthelleprec
#endif /*PNNL_STATS*/

#ifdef CLUBB
use params, only: docloud, dosmoke
use sgs_params, only: doclubb
use grid, only: nz
use clubb_api_module, only: &
  clubb_at_least_debug_level_api, &
  fill_holes_vertical_api, &
  core_rknd
use clubbvars, only: wp2, cloud_frac, rho_ds_zt, rho_ds_zm ! are used, but not modified here
use vars, only: qcl ! Used here and updated in micro_diagnose
use vars, only: prespot ! exner^-1
use module_mp_GRAUPEL, only: &
  cloud_frac_thresh ! Threshold for using sgs cloud fraction to weight 
                    ! microphysical quantities [%]

#endif
#ifdef SILHS
use clubb_silhs_vars, only: &
  lh_microphys_type, & ! Variables
  lh_microphys_interactive, &
  lh_microphys_non_interactive, &
  lh_microphys_disabled
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
  qr_auto, qr_accr, qr_evap
real, dimension(nzm) :: cloud_frac_in
#endif /*CLUBB*/

#if defined(CLUBB) || defined(UWM_STATS)
real, dimension(nzm) :: rain_vel
#endif


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

! Compare the difference between updating 't' using sedimentation or state
! variables
real t_state(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) 
real t_sed(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) 

#endif /*UWM_STATS*/

#ifdef PNNL_STATS
!MWSWong:budgets for THL and QTOG
real thel_before(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) 
real thel_avg(nzm), thel_before_avg(nzm)
real qtog_before(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm) 
real qtog_avg(nzm), qtog_before_avg(nzm)
#endif /*PNNL_STATS*/

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
#ifdef PNNL_STATS
   thellat(:) = 0.
#endif /*PNNL_STATS*/
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

   mPRE_3D=0.0
   mPRA_3D=0.0
   mPRC_3D=0.0
   mPRACS_3D=0.0
   mMNUCCR_3D=0.0
   mQMULTR_3D=0.0
   mQMULTRG_3D=0.0
   mPIACR_3D=0.0
   mPIACRS_3D=0.0
   mPRACG_3D=0.0
   mPGRACS_3D=0.0
   mPSMLT_3D=0.0
   mPGMLT_3D=0.0
   mQR_INST_3D=0.0

   mSIZEFIX_NR_3D=0.0
   mNEGFIX_NR_3D=0.0
   mNSUBR_3D=0.0
   mNSMLTR_3D=0.0
   mNGMLTR_3D=0.0
   mNPRC1_3D=0.0
   mNPRACS_3D=0.0
   mNNUCCR_3D=0.0
   mNRAGG_3D=0.0
   mNIACR_3D=0.0
   mNIACRS_3D=0.0
   mNPRACG_3D=0.0
   mNGRACS_3D=0.0
   mNR_INST_3D=0.0
   
   rain_vel_3D=0.0
   EFFR_3D=0.0

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
      if(dopredictNc) tmpncl(:) = micro_field(i,j,:,incl)
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

#ifdef PNNL_STATS
      do k=1,nzm 
!MWSWong: add budget for THL and QTOG
         thel_before(i,j,k) = t2thetal( t(i,j,k), gamaz(k), tmpqr(k), &
                                        tmpqci(k), tmpqs(k)+tmpqg(k), &
                                        prespot(k) )
         qtog_before(i,j,k) = micro_field(i,j,k,iqv)
         do n=2,nmicro_fields ! prevents accessing unreferenced memory 
            qtog_before(i,j,k) = qtog_before(i,j,k) + flag_wmass(n)*micro_field(i,j,k,n)
         end do
      end do 
#endif /*PNNL_STATS*/ 
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
        if ( clubb_at_least_debug_level_api( 1 ) ) then
          write(0,*) "M2005 has received a negative water vapor"
        end if
        call fill_holes_vertical_api( 2, 0._core_rknd, "zt", rho_ds_zt, rho_ds_zm, qv_clip )
        tmpqv = qv_clip(2:nz)
      end if
      if ( any( tmpqcl < 0. ) ) then
        qcl_clip(2:nz) = tmpqcl(1:nzm)
        qcl_clip(1) = 0.0_core_rknd
        if ( clubb_at_least_debug_level_api( 1 ) ) then
          write(0,*) "M2005 has received a negative liquid water"
        end if
        call fill_holes_vertical_api( 2, 0._core_rknd, "zt", rho_ds_zt, rho_ds_zm, qcl_clip )
        tmpqcl = qcl_clip(2:nz)
      end if

      ! Set autoconversion and accretion rates to 0;  these are diagnostics and
      ! don't feed back into the calculation.
      PRC = 0.
      PRA = 0.
      PRE = 0.
#endif /*CLUBB*/
      
      ! explanation of variable names:
      !   mtend1d: array of 1d profiles of microphysical tendencies (w/ sed.)
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
#if defined(CLUBB) || defined(UWM_STATS)           
           rain_vel,&
#endif /*CLUBB or UWM_STATS*/
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
        if ( clubb_at_least_debug_level_api( 1 ) ) then
          write(0,*) "M2005 has produced a negative water vapor"
        end if
        call fill_holes_vertical_api( 2, 0._core_rknd, "zt", rho_ds_zt, rho_ds_zm, qv_clip )
        tmpqv = qv_clip(2:nz)
      end if
      if ( any( tmpqcl < 0. ) ) then
        qcl_clip(2:nz) = tmpqcl(1:nzm)
        qcl_clip(1) = 0.0_core_rknd
        if ( clubb_at_least_debug_level_api( 1 ) ) then
          write(0,*) "M2005 has produced a negative liquid water"
        end if
        call fill_holes_vertical_api( 2, 0._core_rknd, "zt", rho_ds_zt, rho_ds_zm, qcl_clip )
        tmpqcl = qcl_clip(2:nz)
      end if
#endif /*CLUBB*/

#ifdef UWM_STATS
! We want the tendencies to be explicity: (before - after) / dt
! Therefore, we overwrite the output from Morrison before the fields are updated
! below. Note: mtend does not include sedimentation here. 
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
#endif /*UWM_STATS*/
endif

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
      if(dopredictNc) micro_field(i,j,:,incl) = tmpncl(:)
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

#ifdef UWM_STATS
      ! Compare the difference between updating 't' using sedimentation or state
      ! variables

      ! update t using state variables
      t_state(i,j,:) = tmptabs(:) + gamaz(:) &
                      - fac_cond*( tmpqcl(:) + tmpqr(:) ) &
                      - fac_sub*( tmpqci(:) + tmpqs(:) + tmpqg(:) )

      ! update liquid-ice static energy due to precipitation
      t_sed(i,j,:) = t(i,j,:) &
           - dtn*fac_cond*(stendqcl+stendqr) &
           - dtn*fac_sub*(stendqci+stendqs+stendqg)

      ! Compute the residual between the two
      t_res(i,j,:) = t_state(i,j,:) - t_before(i,j,:) &
                    - dtn*fac_cond*(stendqcl+stendqr) &
                    - dtn*fac_sub*(stendqci+stendqs+stendqg)

      ! Update t
      if(l_update_mse_using_state_vars) then
          t(i,j,:) = t_state(i,j,:)
      else
          t(i,j,:) = t_sed(i,j,:)
      endif
#else
      !=====================================================
      ! update liquid-ice static energy due to precipitation
      t(i,j,:) = t(i,j,:) &
           - dtn*fac_cond*(stendqcl+stendqr) &
           - dtn*fac_sub*(stendqci+stendqs+stendqg)
      !=====================================================
#endif /*UWM_STATS*/
      if(dostatis) then

#ifdef PNNL_STATS
      !MWSWong: THL budget
      !thel(i,j,1:nzm) = prespot(1:nzm)*( t(i,j,1:nzm) - gamaz(1:nzm)  &
      !                 + fac_cond * micro_field(i,j,:,iqr) &
      !                 + fac_sub  * (micro_field(i,j,:,iqci)+micro_field(i,j,:,iqs)+micro_field(i,j,:,iqg)))
      thel(i,j,1:nzm) = t2thetal( t(i,j,1:nzm), gamaz(1:nzm), tmpqr(1:nzm), &
                                  tmpqci(1:nzm), tmpqs(1:nzm)+tmpqg(1:nzm), &
                                  prespot(1:nzm) )

      thellat(1:nzm) = thellat(1:nzm) + (thel(i,j,1:nzm)-thel_before(i,j,1:nzm))
#endif /*PNNL_STATS*/

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
#if defined( CLUBB ) || defined( UWM_MISC ) /* Bug fix -dschanen 9 Mar 2012 */
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
           
         mNPRACS=mNPRACS+NPRACS
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

#ifdef UWM_STATS
   !3d Budgets
   mPRE_3D(i,j,:)=mPRE_3D(i,j,:) + PRE
   mPRA_3D(i,j,:)=mPRA_3D(i,j,:) + PRA
   mPRC_3D(i,j,:)=mPRC_3D(i,j,:) + PRC
   mPRACS_3D(i,j,:)=mPRACS_3D(i,j,:) + PRACS
   mMNUCCR_3D(i,j,:)=mMNUCCR_3D(i,j,:) + MNUCCR
   mQMULTR_3D(i,j,:)=mQMULTR_3D(i,j,:) + QMULTR
   mQMULTRG_3D(i,j,:)=mQMULTRG_3D(i,j,:) + QMULTRG
   mPIACR_3D(i,j,:)=mPIACR_3D(i,j,:) + PIACR
   mPIACRS_3D(i,j,:)=mPIACRS_3D(i,j,:) + PIACRS
   mPRACG_3D(i,j,:)=mPRACG_3D(i,j,:) + PRACG
   mPGRACS_3D(i,j,:)=mPGRACS_3D(i,j,:) + PGRACS
   mPSMLT_3D(i,j,:)=mPSMLT_3D(i,j,:) + PSMLT
   mPGMLT_3D(i,j,:)=mPGMLT_3D(i,j,:) + PGMLT
   mQR_INST_3D(i,j,:)= mQR_INST_3D(i,j,:) + QR_INST

   mSIZEFIX_NR_3D(i,j,:)=mSIZEFIX_NR_3D(i,j,:) + SIZEFIX_NR
   mNEGFIX_NR_3D(i,j,:)=mNEGFIX_NR_3D(i,j,:) + NEGFIX_NR
   mNSUBR_3D(i,j,:)=mNSUBR_3D(i,j,:) + NSUBR
   mNSMLTR_3D(i,j,:)=mNSMLTR_3D(i,j,:) + NSMLTR
   mNGMLTR_3D(i,j,:)=mNGMLTR_3D(i,j,:) + NGMLTR
   mNPRC1_3D(i,j,:)=mNPRC1_3D(i,j,:) + NPRC1 
   mNPRACS_3D(i,j,:)=mNPRACS_3D(i,j,:) + NPRACS
   mNNUCCR_3D(i,j,:)=mNNUCCR_3D(i,j,:) + NNUCCR
   mNRAGG_3D(i,j,:)=mNRAGG_3D(i,j,:) + NRAGG
   mNIACR_3D(i,j,:)=mNIACR_3D(i,j,:) + NIACR
   mNIACRS_3D(i,j,:)=mNIACRS_3D(i,j,:) + NIACRS
   mNPRACG_3D(i,j,:)=mNPRACG_3D(i,j,:) + NPRACG
   mNGRACS_3D(i,j,:)=mNGRACS_3D(i,j,:) + NGRACS
   mNR_INST_3D(i,j,:)=mNR_INST_3D(i,j,:) + NR_INST

   ! Additional micro. variables
   rain_vel_3D(i,j,:)=rain_vel_3D(i,j,:) + rain_vel
   EFFR_3D(i,j,:)=EFFR_3D(i,j,:) + effr1d
#endif /*UWM_STATS*/
    
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
do j = 1,ny
   do i = 1,nx
      if (doprecip) then
         tmpqr(:) = micro_field(i,j,:,iqr)
      else
         tmpqr(:) = 0.0
      endif
      ! Compute liquid water potential temperature
      theta_l(i,j,:) = t2thetal( t(i,j,:), gamaz(:), tmpqr(:), &
      ! pass in zero for the ice species ==> tmpqci(:), tmpqs(:) + tmpqg(:), prespot(:) )
                                 0.0,        0.0,                  prespot(:) )
   enddo ! i = 1, nx
enddo ! j = 1, ny

! Compute extended liquid water mixing ratio
call compute_chi_eta( theta_l, micro_field(1:nx,1:ny,1:nzm,iqv), pres, prespot,&
                      chi, eta )

#endif /*UWM_STATS*/
if(mod(nstep,nsave3D).eq.0.and.nstep.ge.nsave3Dstart.and.nstep.le.nsave3Dend) then
  call write_3d_micro_fields()
endif

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

#ifdef PNNL_STATS
        !MWSWong: add budget for THL
        do k=1,nzm 
          do j=dimy1_s,dimy2_s
            do i=dimx1_s,dimx2_s
              qtog(i,j,k) = micro_field(i,j,k,index_water_vapor)
              do n=2,nmicro_fields ! prevents accessing unreferenced memory 
                qtog(i,j,k) = qtog(i,j,k) + flag_wmass(n)*micro_field(i,j,k,n)
              end do
            end do 
          end do 
        end do 
        call stat_varscalar(thel,thel_before,thel_avg,thel_before_avg,thel2leprec)
        call stat_varscalar(qtog,qtog_before,qtog_avg,qtog_before_avg,qtog2leprec)

        call setvalue(thelwleprec,nzm,0.)
        call stat_sw2(thel,thel_before,thelwleprec)

        call setvalue(qtogwleprec,nzm,0.)
        call stat_sw2(qtog,qtog_before,qtogwleprec)

        !MWSWong: add budget for covariance r'thl' 
        do k=1,nzm
           qthelleprec(k)=0.
           do j=1,ny
              do i=1,nx
                 qthelleprec(k)=qthelleprec(k) + (thel(i,j,k)-thel_avg(k))*(micro_field(i,j,k,iqv)-qt_avg(k)) &
                              - (thel_before(i,j,k)-thel_before_avg(k))*(qt_before(i,j,k)-qt_avg_before(k))
              end do
           end do
           qthelleprec(k)=qthelleprec(k)*(1./(dtn*nx*ny))
        end do

#endif /*PNNL_STATS*/

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
#if defined( CLUBB ) || defined( UWM_MISC ) /* Bug fix -dschanen 9 Mar 2012 */
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
#if defined( CLUBB ) || defined( UWM_MISC ) /* Bug fix -dschanen 9 Mar 2012 */
      if ( dograupel ) then
        mksed(m,iqg) = tmpg
      end if
#else
      mksed(m,iqg) = tmpg
#endif
   end do
#if defined( CLUBB ) || defined( UWM_MISC ) /* Bug fix -dschanen 9 Mar 2012 */
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
use clubb_api_module, only: &
  clubb_at_least_debug_level_api, & ! Procedure
  fstderr, zero_threshold
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
          if ( clubb_at_least_debug_level_api( 1 ) ) then
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
integer :: idx_frac_fields ! Index for each hydrometeor fraction
#endif /*UWM_STATS*/
character(*) namelist(*), deflist(*), unitlist(*)
integer status(*),average_type(*),count,microcount, n, ii, jj, ncond

#ifdef UWM_STATS
character*20 name
integer m
character*10 frac_in_char_temp
real curr_frac
#else
character*8 name
#endif /*UWM_STATS*/
character*80 longname
character*10 units

microcount = 0

#ifdef UWM_STATS
name = 'hl_on_Cp_res'
longname = 'Change in hl over Cp when using state variables instead of precipitation'
units = 'K'
call add_to_namelist(count,microcount,name,longname,units,0)
#endif /*UWM_STATS*/
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
! In-cloud  Correlations 
!----------------------------

!----------------------------
! chi(S_mellor) 
!----------------------------
  name = 'corr_chi_w'
  longname = '  Correlation of chi and vertical velocity'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_cor_chi_w = count

if(dopredictNc) then
  name = 'corr_chi_Nc'
  longname = '  Correlation of chi and cloud droplet concentration'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_cor_chi_Nc = count
end if ! do predictNc
    

  if (doprecip) then
    
    name = 'corr_chi_rr'
    longname = '  Correlation of chi and rain water mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_chi_rr = count
    
    name = 'corr_chi_Nr'
    longname = '  Correlation of chi and rain drop concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_chi_Nr = count
    
  end if ! doprecip

  if (doicemicro) then
    
    name = 'corr_chi_ri'
    longname = '  Correlation of chi and cloud ice mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_chi_ri = count
    
    name = 'corr_chi_Ni'
    longname = '  Correlation of chi and cloud ice concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_chi_Ni = count
    
    name = 'corr_chi_rs'
    longname = '  Correlation of chi and snow mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_chi_rs = count
    
    name = 'corr_chi_Ns'
    longname = '  Correlation of chi and snowflake concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_chi_Ns = count
    

    if (dograupel) then
    
      name = 'corr_chi_rg'
      longname = '  Correlation of chi and graupel mixing ratio'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_cor_chi_rg = count
    
      name = 'corr_chi_Ng'
      longname = '  Correlation of chi and graupel concentration'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_cor_chi_Ng = count
    
    end if ! dograupel
  endif ! doicemicro
!----------------------------
! Vertical Velocity 
!----------------------------

if(dopredictNc) then
  name = 'corr_w_Nc'
  longname = '  Correlation of vertical velocity and cloud droplet concentration'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_cor_w_Nc = count
end if ! do predictNc
    

  if (doprecip) then
    
    name = 'corr_w_rr'
    longname = '  Correlation of vertical velocity and rain water mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_w_rr = count
    
    name = 'corr_w_Nr'
    longname = '  Correlation of vertical velocity and rain drop concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_w_Nr = count
    
  end if ! doprecip

  if (doicemicro) then
    
    name = 'corr_w_ri'
    longname = '  Correlation of vertical velocity and cloud ice mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_w_ri = count
    
    name = 'corr_w_Ni'
    longname = '  Correlation of vertical velocity and cloud ice concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_w_Ni = count
    
    name = 'corr_w_rs'
    longname = '  Correlation of vertical velocity and snow mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_w_rs = count
    
    name = 'corr_w_Ns'
    longname = '  Correlation of vertical velocity and snowflake concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_w_Ns = count
    

    if (dograupel) then
    
      name = 'corr_w_rg'
      longname = '  Correlation of vertical velocity and graupel mixing ratio'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_cor_w_rg = count
    
      name = 'corr_w_Ng'
      longname = '  Correlation of vertical velocity and graupel concentration'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_cor_w_Ng = count
    
    end if ! dograupel
  end if ! doicemicro

!----------------------------
! Cloud droplet number concentration 
!----------------------------
if(dopredictNc) then
  if (doprecip) then
        
    name = 'corr_Nc_rr'
    longname = '  Correlation of cloud droplet conc. and rain water mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_Nc_rr = count
        
    name = 'corr_Nc_Nr'
    longname = '  Correlation of cloud droplet conc. and rain drop conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_Nc_Nr = count
        
  end if ! doprecip
    
  if (doicemicro) then
        
    name = 'corr_Nc_ri'
    longname = '  Correlation of cloud droplet conc. and cloud ice mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_Nc_ri = count
        
    name = 'corr_Nc_Ni'
    longname = '  Correlation of cloud droplet conc. and cloud ice conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_Nc_Ni = count
       
    name = 'corr_Nc_rs'
    longname = '  Correlation of cloud droplet conc. and snow mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_Nc_rs = count
      
    name = 'corr_Nc_Ns'
    longname = '  Correlation of cloud droplet conc. and snowflake conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_Nc_Ns = count
        
    if (dograupel) then
        
      name = 'corr_Nc_rg'
      longname = '  Correlation of cloud droplet conc. and graupel mixing ratio'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_cor_Nc_rg = count
        
      name = 'corr_Nc_Ng'
      longname = '  Correlation of cloud droplet conc. and graupel conc.'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_cor_Nc_Ng = count
        
    end if ! dograupel
  end if ! doicemicro 
end if ! dopredictNc

!----------------------------
! Rainwater mixing ratio 
!----------------------------
if (doprecip) then

    name = 'corr_rr_Nr'
    longname = '  Correlation of rain water mixing ratio and rain drop concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_rr_Nr = count
  
  if(doicemicro) then

    name = 'corr_rr_ri'
    longname = '  Correlation of rain water mixing ratio and cloud ice mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_rr_ri = count
        
    name = 'corr_rr_Ni'
    longname = '  Correlation of rain water mixing ratio and cloud ice conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_rr_Ni = count
        
    name = 'corr_rr_rs'
    longname = '  Correlation of rain water mixing ratio and snow mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_rr_rs = count
        
    name = 'corr_rr_Ns'
    longname = '  Correlation of rain water mixing ratio and snowflake conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_rr_Ns = count
        

    if(dograupel) then

      name = 'corr_rr_rg'
      longname = '  Correlation of rain water mixing ratio and graupel mixing ratio'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_cor_rr_rg = count
        
      name = 'corr_rr_Ng'
      longname = '  Correlation of rain water mixing ratio and graupel conc.'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_cor_rr_Ng = count
        
    end if !dograupel

!----------------------------
! Rainwater number concentration 
!----------------------------
        
    name = 'corr_Nr_ri'
    longname = '  Correlation of rain drop conc. and cloud ice mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_Nr_ri = count
        
    name = 'corr_Nr_Ni'
    longname = '  Correlation of rain drop conc. and cloud ice conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_Nr_Ni = count
        
    name = 'corr_Nr_rs'
    longname = '  Correlation of rain drop conc. and snow mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_Nr_rs = count
        
    name = 'corr_Nr_Ns'
    longname = '  Correlation of rain drop conc. and snowflake conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_Nr_Ns = count
        
      if(dograupel) then

        name = 'corr_Nr_rg'
        longname = '  Correlation of rain drop conc. and graupel mixing ratio'
        units = '-'
        call add_to_namelist(count,microcount,name,longname,units,corr_avg)
        ! Note:  it is important to place these after the call to add_to_namelist.
        idx_cor_Nr_rg = count
        
        name = 'corr_Nr_Ng'
        longname = '  Correlation of rain drop conc. and graupel conc.'
        units = '-'
        call add_to_namelist(count,microcount,name,longname,units,corr_avg)
        ! Note:  it is important to place these after the call to add_to_namelist.
        idx_cor_Nr_Ng = count
        
      end if !dograupel
    end if !doicemicro
end if !doprecip

!----------------------------
! Cloud-ice mixing ratio 
!----------------------------
if (doicemicro) then
    
  name = 'corr_ri_Ni'
  longname = '  Correlation of cloud ice mixing ratio and cloud ice conc.'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_cor_ri_Ni = count
    
  name = 'corr_ri_rs'
  longname = '  Correlation of cloud ice mixing ratio and snow mixing ratio'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_cor_ri_rs = count
    
  name = 'corr_ri_Ns'
  longname = '  Correlation of cloud ice mixing ratio and snowflake conc.'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_cor_ri_Ns = count
    
  if (dograupel) then

    name = 'corr_ri_rg'
    longname = '  Correlation of cloud ice mixing ratio and graupel mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_ri_rg = count
        
    name = 'corr_ri_Ng'
    longname = '  Correlation of cloud ice mixing ratio and graupel conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_ri_Ng = count
        
  end if !dograupel

!----------------------------
! Cloud-ice number concentration 
!----------------------------
    
  name = 'corr_Ni_rs'
  longname = '  Correlation of cloud ice concentration and snow mixing ratio'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_cor_Ni_rs = count
    
  name = 'corr_Ni_Ns'
  longname = '  Correlation of cloud ice concentration and snowflake conc.'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_cor_Ni_Ns = count
    
  if (dograupel) then

    name = 'corr_Ni_rg'
    longname = '  Correlation of cloud ice conc. and graupel mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_Ni_rg = count
        
    name = 'corr_Ni_Ng'
    longname = '  Correlation of cloud ice conc. and graupel conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_Ni_Ng = count
        
  end if !dograupel

!----------------------------
! Snow mixing ratio 
!----------------------------
  
  name = 'corr_rs_Ns'
  longname = '  Correlation of snow mixing ratio and snowflake concentration'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_cor_rs_Ns = count
    
  if (dograupel) then

    name = 'corr_rs_rg'
    longname = '  Correlation of snow mixing ratio and graupel mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_rs_rg = count
        
    name = 'corr_rs_Ng'
    longname = '  Correlation of snow mixing ratio and graupel concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_rs_Ng = count
        
!----------------------------
! Snow number concentration 
!----------------------------
    name = 'corr_Ns_rg'
    longname = '  Correlation of snowflake concentration and graupel mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_Ns_rg = count
        
    name = 'corr_Ns_Ng'
    longname = '  Correlation of snowflake concentration and graupel conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_Ns_Ng = count

!----------------------------
! Graupel mixing ratio / number concentration 
!----------------------------
        
    name = 'corr_rg_Ng'
    longname = '  Correlation of graupel mixing ratio and graupel concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_cor_rg_Ng = count
        
  end if !dograupel
end if !doicemicro



!---------------------------
! In-cloud ic_correlations 
!----------------------------

!----------------------------
! chi(S_mellor) 
!----------------------------
  name = 'ic_corr_chi_w'
  longname = 'In-Cloud Correlation of chi and vertical velocity'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_ic_cor_chi_w = count

if(dopredictNc) then
  name = 'ic_corr_chi_Nc'
  longname = 'In-Cloud Correlation of chi and cloud droplet concentration'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_ic_cor_chi_Nc = count
end if ! do predictNc
    

  if (doprecip) then
    
    name = 'ic_corr_chi_rr'
    longname = 'In-Cloud Correlation of chi and rain water mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_chi_rr = count
    
    name = 'ic_corr_chi_Nr'
    longname = 'In-Cloud Correlation of chi and rain drop concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_chi_Nr = count
    
  end if ! doprecip

  if (doicemicro) then
    
    name = 'ic_corr_chi_ri'
    longname = 'In-Cloud Correlation of chi and cloud ice mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_chi_ri = count
    
    name = 'ic_corr_chi_Ni'
    longname = 'In-Cloud Correlation of chi and cloud ice concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_chi_Ni = count
    
    name = 'ic_corr_chi_rs'
    longname = 'In-Cloud Correlation of chi and snow mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_chi_rs = count
    
    name = 'ic_corr_chi_Ns'
    longname = 'In-Cloud Correlation of chi and snowflake concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_chi_Ns = count
    

    if (dograupel) then
    
      name = 'ic_corr_chi_rg'
      longname = 'In-Cloud Correlation of chi and graupel mixing ratio'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_ic_cor_chi_rg = count
    
      name = 'ic_corr_chi_Ng'
      longname = 'In-Cloud Correlation of chi and graupel concentration'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_ic_cor_chi_Ng = count
    
    end if ! dograupel
  endif ! doicemicro
!----------------------------
! Vertical Velocity 
!----------------------------

if(dopredictNc) then
  name = 'ic_corr_w_Nc'
  longname = 'In-Cloud Correlation of vertical velocity and cloud droplet concentration'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_ic_cor_w_Nc = count
end if ! do predictNc
    

  if (doprecip) then
    
    name = 'ic_corr_w_rr'
    longname = 'In-Cloud Correlation of vertical velocity and rain water mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_w_rr = count
    
    name = 'ic_corr_w_Nr'
    longname = 'In-Cloud Correlation of vertical velocity and rain drop concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_w_Nr = count
    
  end if ! doprecip

  if (doicemicro) then
    
    name = 'ic_corr_w_ri'
    longname = 'In-Cloud Correlation of vertical velocity and cloud ice mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_w_ri = count
    
    name = 'ic_corr_w_Ni'
    longname = 'In-Cloud Correlation of vertical velocity and cloud ice concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_w_Ni = count
    
    name = 'ic_corr_w_rs'
    longname = 'In-Cloud Correlation of vertical velocity and snow mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_w_rs = count
    
    name = 'ic_corr_w_Ns'
    longname = 'In-Cloud Correlation of vertical velocity and snowflake concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_w_Ns = count
    

    if (dograupel) then
    
      name = 'ic_corr_w_rg'
      longname = 'In-Cloud Correlation of vertical velocity and graupel mixing ratio'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_ic_cor_w_rg = count
    
      name = 'ic_corr_w_Ng'
      longname = 'In-Cloud Correlation of vertical velocity and graupel concentration'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_ic_cor_w_Ng = count
    
    end if ! dograupel
  end if ! doicemicro

!----------------------------
! Cloud droplet number concentration 
!----------------------------
if(dopredictNc) then
  if (doprecip) then
        
    name = 'ic_corr_Nc_rr'
    longname = 'In-Cloud Correlation of cloud droplet conc. and rain water mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_Nc_rr = count
        
    name = 'ic_corr_Nc_Nr'
    longname = 'In-Cloud Correlation of cloud droplet conc. and rain drop conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_Nc_Nr = count
        
  end if ! doprecip
    
  if (doicemicro) then
        
    name = 'ic_corr_Nc_ri'
    longname = 'In-Cloud Correlation of cloud droplet conc. and cloud ice mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_Nc_ri = count
        
    name = 'ic_corr_Nc_Ni'
    longname = 'In-Cloud Correlation of cloud droplet conc. and cloud ice conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_Nc_Ni = count
       
    name = 'ic_corr_Nc_rs'
    longname = 'In-Cloud Correlation of cloud droplet conc. and snow mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_Nc_rs = count
      
    name = 'ic_corr_Nc_Ns'
    longname = 'In-Cloud Correlation of cloud droplet conc. and snowflake conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_Nc_Ns = count
        
    if (dograupel) then
        
      name = 'ic_corr_Nc_rg'
      longname = 'In-Cloud Correlation of cloud droplet conc. and graupel mixing ratio'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_ic_cor_Nc_rg = count
        
      name = 'ic_corr_Nc_Ng'
      longname = 'In-Cloud Correlation of cloud droplet conc. and graupel conc.'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_ic_cor_Nc_Ng = count
        
    end if ! dograupel
  end if ! doicemicro 
end if ! dopredictNc

!----------------------------
! Rainwater mixing ratio 
!----------------------------
if (doprecip) then

    name = 'ic_corr_rr_Nr'
    longname = 'In-Cloud Correlation of rain water mixing ratio and rain drop concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_rr_Nr = count
  
  if(doicemicro) then

    name = 'ic_corr_rr_ri'
    longname = 'In-Cloud Correlation of rain water mixing ratio and cloud ice mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_rr_ri = count
        
    name = 'ic_corr_rr_Ni'
    longname = 'In-Cloud Correlation of rain water mixing ratio and cloud ice conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_rr_Ni = count
        
    name = 'ic_corr_rr_rs'
    longname = 'In-Cloud Correlation of rain water mixing ratio and snow mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_rr_rs = count
        
    name = 'ic_corr_rr_Ns'
    longname = 'In-Cloud Correlation of rain water mixing ratio and snowflake conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_rr_Ns = count
        

    if(dograupel) then

      name = 'ic_corr_rr_rg'
      longname = 'In-Cloud Correlation of rain water mixing ratio and graupel mixing ratio'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_ic_cor_rr_rg = count
        
      name = 'ic_corr_rr_Ng'
      longname = 'In-Cloud Correlation of rain water mixing ratio and graupel conc.'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_ic_cor_rr_Ng = count
        
    end if !dograupel

!----------------------------
! Rainwater number concentration 
!----------------------------
        
    name = 'ic_corr_Nr_ri'
    longname = 'In-Cloud Correlation of rain drop conc. and cloud ice mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_Nr_ri = count
        
    name = 'ic_corr_Nr_Ni'
    longname = 'In-Cloud Correlation of rain drop conc. and cloud ice conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_Nr_Ni = count
        
    name = 'ic_corr_Nr_rs'
    longname = 'In-Cloud Correlation of rain drop conc. and snow mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_Nr_rs = count
        
    name = 'ic_corr_Nr_Ns'
    longname = 'In-Cloud Correlation of rain drop conc. and snowflake conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_Nr_Ns = count
        
      if(dograupel) then

        name = 'ic_corr_Nr_rg'
        longname = 'In-Cloud Correlation of rain drop conc. and graupel mixing ratio'
        units = '-'
        call add_to_namelist(count,microcount,name,longname,units,corr_avg)
        ! Note:  it is important to place these after the call to add_to_namelist.
        idx_ic_cor_Nr_rg = count
        
        name = 'ic_corr_Nr_Ng'
        longname = 'In-Cloud Correlation of rain drop conc. and graupel conc.'
        units = '-'
        call add_to_namelist(count,microcount,name,longname,units,corr_avg)
        ! Note:  it is important to place these after the call to add_to_namelist.
        idx_ic_cor_Nr_Ng = count
        
      end if !dograupel
    end if !doicemicro
end if !doprecip

!----------------------------
! Cloud-ice mixing ratio 
!----------------------------
if (doicemicro) then
    
  name = 'ic_corr_ri_Ni'
  longname = 'In-Cloud Correlation of cloud ice mixing ratio and cloud ice conc.'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_ic_cor_ri_Ni = count
    
  name = 'ic_corr_ri_rs'
  longname = 'In-Cloud Correlation of cloud ice mixing ratio and snow mixing ratio'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_ic_cor_ri_rs = count
    
  name = 'ic_corr_ri_Ns'
  longname = 'In-Cloud Correlation of cloud ice mixing ratio and snowflake conc.'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_ic_cor_ri_Ns = count
    
  if (dograupel) then

    name = 'ic_corr_ri_rg'
    longname = 'In-Cloud Correlation of cloud ice mixing ratio and graupel mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_ri_rg = count
        
    name = 'ic_corr_ri_Ng'
    longname = 'In-Cloud Correlation of cloud ice mixing ratio and graupel conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_ri_Ng = count
        
  end if !dograupel

!----------------------------
! Cloud-ice number concentration 
!----------------------------
    
  name = 'ic_corr_Ni_rs'
  longname = 'In-Cloud Correlation of cloud ice concentration and snow mixing ratio'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_ic_cor_Ni_rs = count
    
  name = 'ic_corr_Ni_Ns'
  longname = 'In-Cloud Correlation of cloud ice concentration and snowflake conc.'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_ic_cor_Ni_Ns = count
    
  if (dograupel) then

    name = 'ic_corr_Ni_rg'
    longname = 'In-Cloud Correlation of cloud ice conc. and graupel mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_Ni_rg = count
        
    name = 'ic_corr_Ni_Ng'
    longname = 'In-Cloud Correlation of cloud ice conc. and graupel conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_Ni_Ng = count
        
  end if !dograupel

!----------------------------
! Snow mixing ratio 
!----------------------------
  
  name = 'ic_corr_rs_Ns'
  longname = 'In-Cloud Correlation of snow mixing ratio and snowflake concentration'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_ic_cor_rs_Ns = count
    
  if (dograupel) then

    name = 'ic_corr_rs_rg'
    longname = 'In-Cloud Correlation of snow mixing ratio and graupel mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_rs_rg = count
        
    name = 'ic_corr_rs_Ng'
    longname = 'In-Cloud Correlation of snow mixing ratio and graupel concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_rs_Ng = count
        
!----------------------------
! Snow number concentration 
!----------------------------
    name = 'ic_corr_Ns_rg'
    longname = 'In-Cloud Correlation of snowflake concentration and graupel mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_Ns_rg = count
        
    name = 'ic_corr_Ns_Ng'
    longname = 'In-Cloud Correlation of snowflake concentration and graupel conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_Ns_Ng = count

!----------------------------
! Graupel mixing ratio / number concentration 
!----------------------------
        
    name = 'ic_corr_rg_Ng'
    longname = 'In-Cloud Correlation of graupel mixing ratio and graupel concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_ic_cor_rg_Ng = count
        
  end if !dograupel
end if !doicemicro

!---------------------------
! Out of cloud oc_correlations 
!----------------------------

!----------------------------
! chi(S_mellor) 
!----------------------------
  name = 'oc_corr_chi_w'
  longname = 'In-Cloud Correlation of chi and vertical velocity'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_oc_cor_chi_w = count

if(dopredictNc) then
  name = 'oc_corr_chi_Nc'
  longname = 'In-Cloud Correlation of chi and cloud droplet concentration'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_oc_cor_chi_Nc = count
end if ! do predictNc
    

  if (doprecip) then
    
    name = 'oc_corr_chi_rr'
    longname = 'In-Cloud Correlation of chi and rain water mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_chi_rr = count
    
    name = 'oc_corr_chi_Nr'
    longname = 'In-Cloud Correlation of chi and rain drop concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_chi_Nr = count
    
  end if ! doprecip

  if (doicemicro) then
    
    name = 'oc_corr_chi_ri'
    longname = 'In-Cloud Correlation of chi and cloud ice mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_chi_ri = count
    
    name = 'oc_corr_chi_Ni'
    longname = 'In-Cloud Correlation of chi and cloud ice concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_chi_Ni = count
    
    name = 'oc_corr_chi_rs'
    longname = 'In-Cloud Correlation of chi and snow mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_chi_rs = count
    
    name = 'oc_corr_chi_Ns'
    longname = 'In-Cloud Correlation of chi and snowflake concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_chi_Ns = count
    

    if (dograupel) then
    
      name = 'oc_corr_chi_rg'
      longname = 'In-Cloud Correlation of chi and graupel mixing ratio'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_oc_cor_chi_rg = count
    
      name = 'oc_corr_chi_Ng'
      longname = 'In-Cloud Correlation of chi and graupel concentration'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_oc_cor_chi_Ng = count
    
    end if ! dograupel
  endif ! doicemicro
!----------------------------
! Vertical Velocity 
!----------------------------

if(dopredictNc) then
  name = 'oc_corr_w_Nc'
  longname = 'In-Cloud Correlation of vertical velocity and cloud droplet concentration'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_oc_cor_w_Nc = count
end if ! do predictNc
    

  if (doprecip) then
    
    name = 'oc_corr_w_rr'
    longname = 'In-Cloud Correlation of vertical velocity and rain water mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_w_rr = count
    
    name = 'oc_corr_w_Nr'
    longname = 'In-Cloud Correlation of vertical velocity and rain drop concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_w_Nr = count
    
  end if ! doprecip

  if (doicemicro) then
    
    name = 'oc_corr_w_ri'
    longname = 'In-Cloud Correlation of vertical velocity and cloud ice mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_w_ri = count
    
    name = 'oc_corr_w_Ni'
    longname = 'In-Cloud Correlation of vertical velocity and cloud ice concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_w_Ni = count
    
    name = 'oc_corr_w_rs'
    longname = 'In-Cloud Correlation of vertical velocity and snow mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_w_rs = count
    
    name = 'oc_corr_w_Ns'
    longname = 'In-Cloud Correlation of vertical velocity and snowflake concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_w_Ns = count
    

    if (dograupel) then
    
      name = 'oc_corr_w_rg'
      longname = 'In-Cloud Correlation of vertical velocity and graupel mixing ratio'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_oc_cor_w_rg = count
    
      name = 'oc_corr_w_Ng'
      longname = 'In-Cloud Correlation of vertical velocity and graupel concentration'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_oc_cor_w_Ng = count
    
    end if ! dograupel
  end if ! doicemicro

!----------------------------
! Cloud droplet number concentration 
!----------------------------
if(dopredictNc) then
  if (doprecip) then
        
    name = 'oc_corr_Nc_rr'
    longname = 'In-Cloud Correlation of cloud droplet conc. and rain water mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_Nc_rr = count
        
    name = 'oc_corr_Nc_Nr'
    longname = 'In-Cloud Correlation of cloud droplet conc. and rain drop conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_Nc_Nr = count
        
  end if ! doprecip
    
  if (doicemicro) then
        
    name = 'oc_corr_Nc_ri'
    longname = 'In-Cloud Correlation of cloud droplet conc. and cloud ice mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_Nc_ri = count
        
    name = 'oc_corr_Nc_Ni'
    longname = 'In-Cloud Correlation of cloud droplet conc. and cloud ice conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_Nc_Ni = count
       
    name = 'oc_corr_Nc_rs'
    longname = 'In-Cloud Correlation of cloud droplet conc. and snow mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_Nc_rs = count
      
    name = 'oc_corr_Nc_Ns'
    longname = 'In-Cloud Correlation of cloud droplet conc. and snowflake conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_Nc_Ns = count
        
    if (dograupel) then
        
      name = 'oc_corr_Nc_rg'
      longname = 'In-Cloud Correlation of cloud droplet conc. and graupel mixing ratio'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_oc_cor_Nc_rg = count
        
      name = 'oc_corr_Nc_Ng'
      longname = 'In-Cloud Correlation of cloud droplet conc. and graupel conc.'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_oc_cor_Nc_Ng = count
        
    end if ! dograupel
  end if ! doicemicro 
end if ! dopredictNc

!----------------------------
! Rainwater mixing ratio 
!----------------------------
if (doprecip) then

    name = 'oc_corr_rr_Nr'
    longname = 'In-Cloud Correlation of rain water mixing ratio and rain drop concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_rr_Nr = count
  
  if(doicemicro) then

    name = 'oc_corr_rr_ri'
    longname = 'In-Cloud Correlation of rain water mixing ratio and cloud ice mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_rr_ri = count
        
    name = 'oc_corr_rr_Ni'
    longname = 'In-Cloud Correlation of rain water mixing ratio and cloud ice conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_rr_Ni = count
        
    name = 'oc_corr_rr_rs'
    longname = 'In-Cloud Correlation of rain water mixing ratio and snow mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_rr_rs = count
        
    name = 'oc_corr_rr_Ns'
    longname = 'In-Cloud Correlation of rain water mixing ratio and snowflake conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_rr_Ns = count
        

    if(dograupel) then

      name = 'oc_corr_rr_rg'
      longname = 'In-Cloud Correlation of rain water mixing ratio and graupel mixing ratio'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_oc_cor_rr_rg = count
        
      name = 'oc_corr_rr_Ng'
      longname = 'In-Cloud Correlation of rain water mixing ratio and graupel conc.'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,corr_avg)
      ! Note:  it is important to place these after the call to add_to_namelist.
      idx_oc_cor_rr_Ng = count
        
    end if !dograupel

!----------------------------
! Rainwater number concentration 
!----------------------------
        
    name = 'oc_corr_Nr_ri'
    longname = 'In-Cloud Correlation of rain drop conc. and cloud ice mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_Nr_ri = count
        
    name = 'oc_corr_Nr_Ni'
    longname = 'In-Cloud Correlation of rain drop conc. and cloud ice conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_Nr_Ni = count
        
    name = 'oc_corr_Nr_rs'
    longname = 'In-Cloud Correlation of rain drop conc. and snow mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_Nr_rs = count
        
    name = 'oc_corr_Nr_Ns'
    longname = 'In-Cloud Correlation of rain drop conc. and snowflake conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_Nr_Ns = count
        
      if(dograupel) then

        name = 'oc_corr_Nr_rg'
        longname = 'In-Cloud Correlation of rain drop conc. and graupel mixing ratio'
        units = '-'
        call add_to_namelist(count,microcount,name,longname,units,corr_avg)
        ! Note:  it is important to place these after the call to add_to_namelist.
        idx_oc_cor_Nr_rg = count
        
        name = 'oc_corr_Nr_Ng'
        longname = 'In-Cloud Correlation of rain drop conc. and graupel conc.'
        units = '-'
        call add_to_namelist(count,microcount,name,longname,units,corr_avg)
        ! Note:  it is important to place these after the call to add_to_namelist.
        idx_oc_cor_Nr_Ng = count
        
      end if !dograupel
    end if !doicemicro
end if !doprecip

!----------------------------
! Cloud-ice mixing ratio 
!----------------------------
if (doicemicro) then
    
  name = 'oc_corr_ri_Ni'
  longname = 'In-Cloud Correlation of cloud ice mixing ratio and cloud ice conc.'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_oc_cor_ri_Ni = count
    
  name = 'oc_corr_ri_rs'
  longname = 'In-Cloud Correlation of cloud ice mixing ratio and snow mixing ratio'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_oc_cor_ri_rs = count
    
  name = 'oc_corr_ri_Ns'
  longname = 'In-Cloud Correlation of cloud ice mixing ratio and snowflake conc.'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_oc_cor_ri_Ns = count
    
  if (dograupel) then

    name = 'oc_corr_ri_rg'
    longname = 'In-Cloud Correlation of cloud ice mixing ratio and graupel mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_ri_rg = count
        
    name = 'oc_corr_ri_Ng'
    longname = 'In-Cloud Correlation of cloud ice mixing ratio and graupel conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_ri_Ng = count
        
  end if !dograupel

!----------------------------
! Cloud-ice number concentration 
!----------------------------
    
  name = 'oc_corr_Ni_rs'
  longname = 'In-Cloud Correlation of cloud ice concentration and snow mixing ratio'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_oc_cor_Ni_rs = count
    
  name = 'oc_corr_Ni_Ns'
  longname = 'In-Cloud Correlation of cloud ice concentration and snowflake conc.'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_oc_cor_Ni_Ns = count
    
  if (dograupel) then

    name = 'oc_corr_Ni_rg'
    longname = 'In-Cloud Correlation of cloud ice conc. and graupel mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_Ni_rg = count
        
    name = 'oc_corr_Ni_Ng'
    longname = 'In-Cloud Correlation of cloud ice conc. and graupel conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_Ni_Ng = count
        
  end if !dograupel

!----------------------------
! Snow mixing ratio 
!----------------------------
  
  name = 'oc_corr_rs_Ns'
  longname = 'In-Cloud Correlation of snow mixing ratio and snowflake concentration'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,corr_avg)
  ! Note:  it is important to place these after the call to add_to_namelist.
  idx_oc_cor_rs_Ns = count
    
  if (dograupel) then

    name = 'oc_corr_rs_rg'
    longname = 'In-Cloud Correlation of snow mixing ratio and graupel mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_rs_rg = count
        
    name = 'oc_corr_rs_Ng'
    longname = 'In-Cloud Correlation of snow mixing ratio and graupel concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_rs_Ng = count
        
!----------------------------
! Snow number concentration 
!----------------------------
    name = 'oc_corr_Ns_rg'
    longname = 'In-Cloud Correlation of snowflake concentration and graupel mixing ratio'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_Ns_rg = count
        
    name = 'oc_corr_Ns_Ng'
    longname = 'In-Cloud Correlation of snowflake concentration and graupel conc.'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_Ns_Ng = count

!----------------------------
! Graupel mixing ratio / number concentration 
!----------------------------
        
    name = 'oc_corr_rg_Ng'
    longname = 'In-Cloud Correlation of graupel mixing ratio and graupel concentration'
    units = '-'
    call add_to_namelist(count,microcount,name,longname,units,corr_avg)
    ! Note:  it is important to place these after the call to add_to_namelist.
    idx_oc_cor_rg_Ng = count
        
  end if !dograupel
end if !doicemicro

!----------------------------
! Fractions
!----------------------------


! This section creates the names and the number of hydrometeor fraction fields
! we have
frac_name(1) = 'rc_frac'
frac_longname(1) = 'Cloud liquid water fraction'
nfrac_fields = 1 

if (doprecip) then 

  frac_name(2) = 'rr_frac'
  frac_longname(2) = 'Rain water fraction'
  nfrac_fields = nfrac_fields + 1

  if (doicemicro) then

    frac_name(3) = 'ri_frac'
    frac_longname(3) = 'Cloud ice fraction'
    nfrac_fields = nfrac_fields + 1

    frac_name(4) = 'rs_frac'
    frac_longname(4) = 'Snow fraction'
    nfrac_fields = nfrac_fields + 1

    if (dograupel) then

      frac_name(5) = 'rg_frac'
      frac_longname(5) = 'Graupel fraction'
      nfrac_fields = nfrac_fields + 1

    endif ! dograupel
  endif ! doicemicro
endif ! doprecip

do n = 1,nfrac_fields ! Loop through number of hydrometeor fraction fields
  do m = 1,nfractions ! Loop through the number of fractions the user requests,
                      ! which is set in the namelist
    
    ! The user specifies an initial threshold to determine fractional coverage.
    ! Increment each subsequet fraction an order of magnitude
    curr_frac = (10**(m-1))*frac_threshold_init

    ! Convert the number to a string
    WRITE(frac_in_char_temp, '(ES8.1)' ) curr_frac 

    ! Kludgy, but replace the decimal (.) and minus (-) with 'p' and 'm'
    ! respectively.
    frac_in_char_temp(3:3) = 'p'
    frac_in_char_temp(6:6) = 'm'

    ! Stor in array
    frac_in_char(m) = frac_in_char_temp  

    ! Create hydrometeor fraction name, append the threshold value used. 
    name = trim(frac_name(n))//'_'//trim(adjustl(frac_in_char(m)))
    longname = trim(frac_longname(n))
    write(*,*) name
    units = '[fraction]'
    call add_to_namelist(count,microcount,name,longname,units,0)

  end do ! nfrac_fields
end do ! nfractions

!----------------------------
! Domain-wide means
!----------------------------
name = 'Ncm'
longname = 'Domain-wide mean of Nc'
units = '[(kg/kg)]'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'chim'
longname = 'Domain-wide mean of chi'
units = '[(kg/kg)]'
call add_to_namelist(count,microcount,name,longname,units,0)

!----------------------------
! Domain-wide variances
!----------------------------

name = 'chip2'
longname = 'Domain-wide variance of cloud liquid number concentration'
units = '[(kg/kg)^2]'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'Ncp2'
longname = 'Domain-wide variance of cloud liquid number concentration'
units = '[(#/kg)^2]'
call add_to_namelist(count,microcount,name,longname,units,0)

if (doprecip) then
  
    name = 'Nrp2'
    longname = 'Domain-wide variance of rain number concentration'
    units = '[(#/kg)^2]'
    call add_to_namelist(count,microcount,name,longname,units,0)

    name = 'rrp2'
    longname = 'Domain-wide variance of rain mixing ratio'
    units = '[(kg/kg)^2]'
    call add_to_namelist(count,microcount,name,longname,units,0)

    if (doicemicro) then
      
      name = 'Nip2'
      longname = 'Domain-wide variance of cloud ice number concentration'
      units = '[(#/kg)^2]'
      call add_to_namelist(count,microcount,name,longname,units,0)
    
      name = 'rip2'
      longname = 'Domain-wide variance of cloud ice mixing ratio'
      units = '[(kg/kg)^2]'
      call add_to_namelist(count,microcount,name,longname,units,0)
      
      name = 'Nsp2'
      longname = 'Domain-wide variance of snow number concentration'
      units = '[(#/kg)^2]'
      call add_to_namelist(count,microcount,name,longname,units,0)

      name = 'rsp2'
      longname = 'Domain-wide variance of snow mixing ratio'
      units = '[(kg/kg)^2]'
      call add_to_namelist(count,microcount,name,longname,units,0)

      if (dograupel) then
    
        name = 'Ngp2'
        longname = 'Domain-wide variance of graupel number concentration'
        units = '[(#/kg)^2]'
        call add_to_namelist(count,microcount,name,longname,units,0)

        name = 'rgp2'
        longname = 'Domain-wide variance of graupel mixing ratio'
        units = '[(kg/kg)^2]'
        call add_to_namelist(count,microcount,name,longname,units,0)

      end if !dograupel
    end if !doicemicro
  end if !doprecip


!----------------------------
! In-cloud means and variances
!----------------------------
name = 'chim_ic'
longname = 'In-cloud mean of chi'
units = '[(kg/kg)]'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'chip2_ic'
longname = 'In-cloud variance of chi'
units = '[(kg/kg)^2]'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'wm_zt_ic'
longname = 'In-cloud mean of w (on zt grid)'
units = '[(m/s)]'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'wp2_zt_ic'
longname = 'In-cloud variance of w (on zt grid)'
units = '[(kg/kg)^2]'
call add_to_namelist(count,microcount,name,longname,units,0)



!----------------------------
! In-precip means and variances
!----------------------------

name = 'Ncp2_ip'
longname = 'Within cloud variance of cloud liquid number concentration'
units = '[(#/kg)^2]'
call add_to_namelist(count,microcount,name,longname,units,0)

name = 'Ncm_ip'
longname = 'Within cloud mean of cloud liquid number concentration'
units = '[(#/kg)]'
call add_to_namelist(count,microcount,name,longname,units,0)

if (doprecip) then
  
    name = 'Nrp2_ip'
    longname = 'Within rain variance of rain number concentration'
    units = '[(#/kg)^2]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
    name = 'Nrm_ip'
    longname = 'Within rain mean of rain number concentration'
    units = '[(#/kg)]'
    call add_to_namelist(count,microcount,name,longname,units,0)

    name = 'rrp2_ip'
    longname = 'Within rain variance of rain mixing ratio'
    units = '[(kg/kg)^2]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
    name = 'rrm_ip'
    longname = 'Within rain mean of rain mixing ratio'
    units = '[(kg/kg)]'
    call add_to_namelist(count,microcount,name,longname,units,0)

    if (doicemicro) then
      
      name = 'Nip2_ip'
      longname = 'Within cloudice variance of cloud ice number concentration'
      units = '[(#/kg)^2]'
      call add_to_namelist(count,microcount,name,longname,units,0)
      
      name = 'Nim_ip'
      longname = 'Within cloudice mean of cloud ice number concentration'
      units = '[(#/kg)]'
      call add_to_namelist(count,microcount,name,longname,units,0)
    
      name = 'rip2_ip'
      longname = 'Within cloudice variance of cloud ice mixing ratio'
      units = '[(kg/kg)^2]'
      call add_to_namelist(count,microcount,name,longname,units,0)
      
      name = 'rim_ip'
      longname = 'Within cloudice mean of cloud ice mixing ratio'
      units = '[(kg/kg)]'
      call add_to_namelist(count,microcount,name,longname,units,0)
      
      name = 'Nsp2_ip'
      longname = 'Within snow variance of snow number concentration'
      units = '[(#/kg)^2]'
      call add_to_namelist(count,microcount,name,longname,units,0)
      
      name = 'Nsm_ip'
      longname = 'Within snow mean of snow number concentration'
      units = '[(#/kg)]'
      call add_to_namelist(count,microcount,name,longname,units,0)

      name = 'rsp2_ip'
      longname = 'Within snow variance of snow mixing ratio'
      units = '[(kg/kg)^2]'
      call add_to_namelist(count,microcount,name,longname,units,0)
      
      name = 'rsm_ip'
      longname = 'Within snow mean of snow mixing ratio'
      units = '[(kg/kg)]'
      call add_to_namelist(count,microcount,name,longname,units,0)

      if (dograupel) then
    
        name = 'Ngp2_ip'
        longname = 'Within graupel variance of graupel number concentration'
        units = '[(#/kg)^2]'
        call add_to_namelist(count,microcount,name,longname,units,0)
        
        name = 'Ngm_ip'
        longname = 'Within graupel mean of graupel number concentration'
        units = '[(#/kg)]'
        call add_to_namelist(count,microcount,name,longname,units,0)

        name = 'rgp2_ip'
        longname = 'Within graupel variance of graupel mixing ratio'
        units = '[(kg/kg)^2]'
        call add_to_namelist(count,microcount,name,longname,units,0)
        
        name = 'rgm_ip'
        longname = 'Within graupel mean of graupel mixing ratio'
        units = '[(kg/kg)]'
        call add_to_namelist(count,microcount,name,longname,units,0)

      end if !dograupel
    end if !doicemicro
  end if !doprecip

!----------------------------
! Covariances
!---------------------------

!----------------------------
! chi(s_mellor)
!---------------------------
name = 'covarnce_chi_w'
longname = 'Covariance of chi and cloud droplet concentration'
units = '[(kg kg^-1)(m s^-1) ]'
call add_to_namelist(count,microcount,name,longname,units,0)

if(dopredictNc) then
  
  name = 'covarnce_chi_Nc'
  longname = 'Covariance of chi and cloud droplet concentration'
  units = '[(kg kg^-1)(#/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

end if

if (doprecip) then
    
  name = 'covarnce_chi_rr'
  longname = 'Covariance of chi and rain water mixing ratio'
  units = '[(kg kg^-1)(kg/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'covarnce_chi_Nr'
  longname = 'Covariance of chi and rain drop concentration'
  units = '[(kg kg^-1)(#/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

end if ! doprecip

if (doicemicro) then
    
  name = 'covarnce_chi_ri'
  longname = 'Covariance of chi and cloud ice mixing ratio'
  units = '[(kg kg^-1)(kg/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'covarnce_chi_Ni'
  longname = 'Covariance of chi and cloud ice concentration'
  units = '[(kg kg^-1)(#/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'covarnce_chi_rs'
  longname = 'Covariance of chi and snow mixing ratio'
  units = '[(kg kg^-1)(kg/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'covarnce_chi_Ns'
  longname = 'Covariance of chi and snowflake concentration'
  units = '[(kg kg^-1)(#/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

end if ! doicemicro

if (dograupel) then
    
  name = 'covarnce_chi_rg'
  longname = 'Covariance of chi and graupel mixing ratio'
  units = '[(kg kg^-1)(kg/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'covarnce_chi_Ng'
  longname = 'Covariance of chi and graupel concentration'
  units = '[(kg kg^-1)(#/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

endif !doicemicro

!----------------------------
! Vertical velocity
!---------------------------
if(dopredictNc) then
  
  name = 'covarnce_w_Nc'
  longname = 'Covariance of vertical velocity and cloud droplet concentration'
  units = '[(m s^-1)(# kg^-1) ]'
  call add_to_namelist(count,microcount,name,longname,units,0)

end if

if (doprecip) then
    
  name = 'covarnce_w_rr'
  longname = 'Covariance of vertical velocity and rain water mixing ratio'
  units = '[(m s^-1)(kg kg^-1) ]'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'covarnce_w_Nr'
  longname = 'Covariance of vertical velocity and rain drop concentration'
  units = '[(m s^-1)(# kg^-1) ]'
  call add_to_namelist(count,microcount,name,longname,units,0)

end if ! doprecip

if (doicemicro) then
    
  name = 'covarnce_w_ri'
  longname = 'Covariance of vertical velocity and cloud ice mixing ratio'
  units = '[(m s^-1)(kg kg^-1) ]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'covarnce_w_Ni'
  longname = 'Covariance of vertical velocity and cloud ice concentration'
  units = '[(m s^-1)(# kg^-1) ]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'covarnce_w_rs'
  longname = 'Covariance of vertical velocity and snow mixing ratio'
  units = '[(m s^-1)(kg kg^-1) ]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'covarnce_w_Ns'
  longname = 'Covariance of vertical velocity and snowflake concentration'
  units = '[(m s^-1)(# kg^-1) ]'
  call add_to_namelist(count,microcount,name,longname,units,0)

end if ! doicemicro

if (dograupel) then
    
  name = 'covarnce_w_rg'
  longname = 'Covariance of vertical velocity and graupel mixing ratio'
  units = '[(m s^-1)(kg kg^-1) ]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'covarnce_w_Ng'
  longname = 'Covariance of vertical velocity and graupel concentration'
  units = '[(m s^-1)(# kg^-1) ]'
  call add_to_namelist(count,microcount,name,longname,units,0)

end if ! dograupel

!----------------------------
! Cloud droplet number concentration
!---------------------------

if(dopredictNc) then
    if (doprecip) then
        
      name = 'covarnce_Nc_rr'
      longname = 'Covariance of cloud droplet conc. and rain water mixing ratio'
      units = '[(# kg^-1)(kg kg^-1)]'
      call add_to_namelist(count,microcount,name,longname,units,0)
        
      name = 'covarnce_Nc_Nr'
      longname = 'Covariance of cloud droplet conc. and rain drop conc.'
      units = '[(# kg^-1)(# kg^-1)]'
      call add_to_namelist(count,microcount,name,longname,units,0)
    
    end if ! doprecip
    
    if (doicemicro) then
        
      name = 'covarnce_Nc_ri'
      longname = 'Covariance of cloud droplet conc. and cloud ice mixing ratio'
      units = '[(# kg^-1)(kg kg^-1)]'
      call add_to_namelist(count,microcount,name,longname,units,0)
        
      name = 'covarnce_Nc_Ni'
      longname = 'Covariance of cloud droplet conc. and cloud ice conc.'
      units = '[(# kg^-1)(# kg^-1)]'
      call add_to_namelist(count,microcount,name,longname,units,0)
        
      name = 'covarnce_Nc_rs'
      longname = 'Covariance of cloud droplet conc. and snow mixing ratio'
      units = '[(# kg^-1)(kg kg^-1)]'
      call add_to_namelist(count,microcount,name,longname,units,0)
    
      name = 'covarnce_Nc_Ns'
      longname = 'Covariance of cloud droplet conc. and snowflake conc.'
      units = '[(# kg^-1)(# kg^-1)]'
      call add_to_namelist(count,microcount,name,longname,units,0)
    
    end if ! doicemicro
    
    if (dograupel) then
        
      name = 'covarnce_Nc_rg'
      longname = 'Covariance of cloud droplet conc. and graupel mixing ratio'
      units = '[(# kg^-1)(kg kg^-1)]'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,0)
        
      name = 'covarnce_Nc_Ng'
      longname = 'Covariance of cloud droplet conc. and graupel conc.'
      units = '[(# kg^-1)(# kg^-1)]'
      call add_to_namelist(count,microcount,name,longname,units,0)

    end if ! dograupel
end if ! do predictNc

!----------------------------
! Rainwater mixing ratio
!---------------------------
if (doprecip) then

    name = 'covarnce_rr_Nr'
    longname = 'Covariance of rain water mixing ratio and rain drop concentration'
    units = '[(kg kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)

  if(doicemicro) then

    name = 'covarnce_rr_ri'
    longname = 'Covariance of rain water mixing ratio and cloud ice mixing ratio'
    units = '[(kg kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
        
    name = 'covarnce_rr_Ni'
    longname = 'Covariance of rain water mixing ratio and cloud ice conc.'
    units = '[(kg kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
    name = 'covarnce_rr_rs'
    longname = 'Covariance of rain water mixing ratio and snow mixing ratio'
    units = '[(kg kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
    name = 'covarnce_rr_Ns'
    longname = 'Covariance of rain water mixing ratio and snowflake conc.'
    units = '[(kg kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
  
  end if !doicemicro

  if(dograupel) then

    name = 'covarnce_rr_rg'
    longname = 'Covariance of rain water mixing ratio and graupel mixing ratio'
    units = '[(kg kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
        
    name = 'covarnce_rr_Ng'
    longname = 'Covariance of rain water mixing ratio and graupel conc.'
    units = '[(kg kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
  end if !dograupel

!----------------------------
! Rainwater number concentration
!---------------------------
  if(doicemicro) then
        
    name = 'covarnce_Nr_ri'
    longname = 'Covariance of rain drop conc. and cloud ice mixing ratio'
    units = '[(# kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
        
    name = 'covarnce_Nr_Ni'
    longname = 'Covariance of rain drop conc. and cloud ice conc.'
    units = '[(# kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
    name = 'covarnce_Nr_rs'
    longname = 'Covariance of rain drop conc. and snow mixing ratio'
    units = '[(# kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
    name = 'covarnce_Nr_Ns'
    longname = 'Covariance of rain drop conc. and snowflake conc.'
    units = '[(# kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
  
  end if !doicemicro

  if(dograupel) then

    name = 'covarnce_Nr_rg'
    longname = 'Covariance of rain drop conc. and graupel mixing ratio'
    units = '[(# kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
        
    name = 'covarnce_Nr_Ng'
    longname = 'Covariance of rain drop conc. and graupel conc.'
    units = '[(# kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
  end if !dograupel
end if !doprecip
!----------------------------
! Cloudice mixing ratio
!---------------------------

if (doicemicro) then
    
  name = 'covarnce_ri_Ni'
  longname = 'Covariance of cloud ice mixing ratio and cloud ice conc.'
  units = '[(kg kg^-1)(# kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'covarnce_ri_rs'
  longname = 'Covariance of cloud ice mixing ratio and snow mixing ratio'
  units = '[(kg kg^-1)(kg kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'covarnce_ri_Ns'
  longname = 'Covariance of cloud ice mixing ratio and snowflake conc.'
  units = '[(kg kg^-1)(# kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  if (dograupel) then

    name = 'covarnce_ri_rg'
    longname = 'Covariance of cloud ice mixing ratio and graupel mixing ratio'
    units = '[(kg kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)

    name = 'covarnce_ri_Ng'
    longname = 'Covariance of cloud ice mixing ratio and graupel conc.'
    units = '[(kg kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
  end if !dograupel

!----------------------------
! Cloudice number concentration
!---------------------------
    
  name = 'covarnce_Ni_rs'
  longname = 'Covariance of cloud ice concentration and snow mixing ratio'
  units = '[(# kg^-1)(kg kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'covarnce_Ni_Ns'
  longname = 'Covariance of cloud ice concentration and snowflake conc.'
  units = '[(# kg^-1)(# kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  if (dograupel) then

    name = 'covarnce_Ni_rg'
    longname = 'Covariance of cloud ice conc. and graupel mixing ratio'
    units = '[(# kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)

    name = 'covarnce_Ni_Ng'
    longname = 'Covariance of cloud ice conc. and graupel conc.'
    units = '[(# kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
  end if !dograupel

!----------------------------
! Snow mixing ratio
!---------------------------
    
  name = 'covarnce_rs_Ns'
  longname = 'Covariance of snow mixing ratio and snowflake concentration'
  units = '[(kg kg^-1)(# kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  if (dograupel) then

    name = 'covarnce_rs_rg'
    longname = 'Covariance of snow mixing ratio and graupel mixing ratio'
    units = '[(kg kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)

    name = 'covarnce_rs_Ng'
    longname = 'Covariance of snow mixing ratio and graupel concentration'
    units = '[(kg kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)

!----------------------------
! Snow number concentration
!---------------------------
    
    name = 'covarnce_Ns_rg'
    longname = 'Covariance of snowflake concentration and graupel mixing ratio'
    units = '[(# kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)

    name = 'covarnce_Ns_Ng'
    longname = 'Covariance of snowflake concentration and graupel conc.'
    units = '[(# kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)

!----------------------------
! Graupel mixing ratio / number concentration
!---------------------------
    
    name = 'covarnce_rg_Ng'
    longname = 'Covariance of graupel mixing ratio and graupel concentration'
    units = '[(kg kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
  end if !dograupel
end if !doicemicro

!----------------------------
! In-cloud covariances
!---------------------------

!----------------------------
! chi(s_mellor)
!---------------------------
name = 'ic_covarnce_chi_w'
longname = 'In-cloud Covariance of chi and cloud droplet concentration'
units = '[(kg kg^-1)(m s^-1) ]'
call add_to_namelist(count,microcount,name,longname,units,0)

if(dopredictNc) then
  
  name = 'ic_covarnce_chi_Nc'
  longname = 'In-cloud Covariance of chi and cloud droplet concentration'
  units = '[(kg kg^-1)(#/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

end if

if (doprecip) then
    
  name = 'ic_covarnce_chi_rr'
  longname = 'In-cloud Covariance of chi and rain water mixing ratio'
  units = '[(kg kg^-1)(kg/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'ic_covarnce_chi_Nr'
  longname = 'In-cloud Covariance of chi and rain drop concentration'
  units = '[(kg kg^-1)(#/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

end if ! doprecip

if (doicemicro) then
    
  name = 'ic_covarnce_chi_ri'
  longname = 'In-cloud Covariance of chi and cloud ice mixing ratio'
  units = '[(kg kg^-1)(kg/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'ic_covarnce_chi_Ni'
  longname = 'In-cloud Covariance of chi and cloud ice concentration'
  units = '[(kg kg^-1)(#/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'ic_covarnce_chi_rs'
  longname = 'In-cloud Covariance of chi and snow mixing ratio'
  units = '[(kg kg^-1)(kg/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'ic_covarnce_chi_Ns'
  longname = 'In-cloud Covariance of chi and snowflake concentration'
  units = '[(kg kg^-1)(#/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

end if ! doicemicro

if (dograupel) then
    
  name = 'ic_covarnce_chi_rg'
  longname = 'In-cloud Covariance of chi and graupel mixing ratio'
  units = '[(kg kg^-1)(kg/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'ic_covarnce_chi_Ng'
  longname = 'In-cloud Covariance of chi and graupel concentration'
  units = '[(kg kg^-1)(#/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

endif !doicemicro

!----------------------------
! Vertical velocity
!---------------------------
if(dopredictNc) then
  
  name = 'ic_covarnce_w_Nc'
  longname = 'In-cloud Covariance of vertical velocity and cloud droplet concentration'
  units = '[(m s^-1)(# kg^-1) ]'
  call add_to_namelist(count,microcount,name,longname,units,0)

end if

if (doprecip) then
    
  name = 'ic_covarnce_w_rr'
  longname = 'In-cloud Covariance of vertical velocity and rain water mixing ratio'
  units = '[(m s^-1)(kg kg^-1) ]'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'ic_covarnce_w_Nr'
  longname = 'In-cloud Covariance of vertical velocity and rain drop concentration'
  units = '[(m s^-1)(# kg^-1) ]'
  call add_to_namelist(count,microcount,name,longname,units,0)

end if ! doprecip

if (doicemicro) then
    
  name = 'ic_covarnce_w_ri'
  longname = 'In-cloud Covariance of vertical velocity and cloud ice mixing ratio'
  units = '[(m s^-1)(kg kg^-1) ]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'ic_covarnce_w_Ni'
  longname = 'In-cloud Covariance of vertical velocity and cloud ice concentration'
  units = '[(m s^-1)(# kg^-1) ]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'ic_covarnce_w_rs'
  longname = 'In-cloud Covariance of vertical velocity and snow mixing ratio'
  units = '[(m s^-1)(kg kg^-1) ]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'ic_covarnce_w_Ns'
  longname = 'In-cloud Covariance of vertical velocity and snowflake concentration'
  units = '[(m s^-1)(# kg^-1) ]'
  call add_to_namelist(count,microcount,name,longname,units,0)

end if ! doicemicro

if (dograupel) then
    
  name = 'ic_covarnce_w_rg'
  longname = 'In-cloud Covariance of vertical velocity and graupel mixing ratio'
  units = '[(m s^-1)(kg kg^-1) ]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'ic_covarnce_w_Ng'
  longname = 'In-cloud Covariance of vertical velocity and graupel concentration'
  units = '[(m s^-1)(# kg^-1) ]'
  call add_to_namelist(count,microcount,name,longname,units,0)

end if ! dograupel

!----------------------------
! Cloud droplet number concentration
!---------------------------

if(dopredictNc) then
    if (doprecip) then
        
      name = 'ic_covarnce_Nc_rr'
      longname = 'In-cloud Covariance of cloud droplet conc. and rain water mixing ratio'
      units = '[(# kg^-1)(kg kg^-1)]'
      call add_to_namelist(count,microcount,name,longname,units,0)
        
      name = 'ic_covarnce_Nc_Nr'
      longname = 'In-cloud Covariance of cloud droplet conc. and rain drop conc.'
      units = '[(# kg^-1)(# kg^-1)]'
      call add_to_namelist(count,microcount,name,longname,units,0)
    
    end if ! doprecip
    
    if (doicemicro) then
        
      name = 'ic_covarnce_Nc_ri'
      longname = 'In-cloud Covariance of cloud droplet conc. and cloud ice mixing ratio'
      units = '[(# kg^-1)(kg kg^-1)]'
      call add_to_namelist(count,microcount,name,longname,units,0)
        
      name = 'ic_covarnce_Nc_Ni'
      longname = 'In-cloud Covariance of cloud droplet conc. and cloud ice conc.'
      units = '[(# kg^-1)(# kg^-1)]'
      call add_to_namelist(count,microcount,name,longname,units,0)
        
      name = 'ic_covarnce_Nc_rs'
      longname = 'In-cloud Covariance of cloud droplet conc. and snow mixing ratio'
      units = '[(# kg^-1)(kg kg^-1)]'
      call add_to_namelist(count,microcount,name,longname,units,0)
    
      name = 'ic_covarnce_Nc_Ns'
      longname = 'In-cloud Covariance of cloud droplet conc. and snowflake conc.'
      units = '[(# kg^-1)(# kg^-1)]'
      call add_to_namelist(count,microcount,name,longname,units,0)
    
    end if ! doicemicro
    
    if (dograupel) then
        
      name = 'ic_covarnce_Nc_rg'
      longname = 'In-cloud Covariance of cloud droplet conc. and graupel mixing ratio'
      units = '[(# kg^-1)(kg kg^-1)]'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,0)
        
      name = 'ic_covarnce_Nc_Ng'
      longname = 'In-cloud Covariance of cloud droplet conc. and graupel conc.'
      units = '[(# kg^-1)(# kg^-1)]'
      call add_to_namelist(count,microcount,name,longname,units,0)

    end if ! dograupel
end if ! do predictNc

!----------------------------
! Rainwater mixing ratio
!---------------------------
if (doprecip) then

    name = 'ic_covarnce_rr_Nr'
    longname = 'In-cloud Covariance of rain water mixing ratio and rain drop concentration'
    units = '[(kg kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)

  if(doicemicro) then

    name = 'ic_covarnce_rr_ri'
    longname = 'In-cloud Covariance of rain water mixing ratio and cloud ice mixing ratio'
    units = '[(kg kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
        
    name = 'ic_covarnce_rr_Ni'
    longname = 'In-cloud Covariance of rain water mixing ratio and cloud ice conc.'
    units = '[(kg kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
    name = 'ic_covarnce_rr_rs'
    longname = 'In-cloud Covariance of rain water mixing ratio and snow mixing ratio'
    units = '[(kg kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
    name = 'ic_covarnce_rr_Ns'
    longname = 'In-cloud Covariance of rain water mixing ratio and snowflake conc.'
    units = '[(kg kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
  
  end if !doicemicro

  if(dograupel) then

    name = 'ic_covarnce_rr_rg'
    longname = 'In-cloud Covariance of rain water mixing ratio and graupel mixing ratio'
    units = '[(kg kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
        
    name = 'ic_covarnce_rr_Ng'
    longname = 'In-cloud Covariance of rain water mixing ratio and graupel conc.'
    units = '[(kg kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
  end if !dograupel

!----------------------------
! Rainwater number concentration
!---------------------------
  if(doicemicro) then
        
    name = 'ic_covarnce_Nr_ri'
    longname = 'In-cloud Covariance of rain drop conc. and cloud ice mixing ratio'
    units = '[(# kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
        
    name = 'ic_covarnce_Nr_Ni'
    longname = 'In-cloud Covariance of rain drop conc. and cloud ice conc.'
    units = '[(# kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
    name = 'ic_covarnce_Nr_rs'
    longname = 'In-cloud Covariance of rain drop conc. and snow mixing ratio'
    units = '[(# kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
    name = 'ic_covarnce_Nr_Ns'
    longname = 'In-cloud Covariance of rain drop conc. and snowflake conc.'
    units = '[(# kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
  
  end if !doicemicro

  if(dograupel) then

    name = 'ic_covarnce_Nr_rg'
    longname = 'In-cloud Covariance of rain drop conc. and graupel mixing ratio'
    units = '[(# kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
        
    name = 'ic_covarnce_Nr_Ng'
    longname = 'In-cloud Covariance of rain drop conc. and graupel conc.'
    units = '[(# kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
  end if !dograupel
end if !doprecip
!----------------------------
! Cloudice mixing ratio
!---------------------------

if (doicemicro) then
    
  name = 'ic_covarnce_ri_Ni'
  longname = 'In-cloud Covariance of cloud ice mixing ratio and cloud ice conc.'
  units = '[(kg kg^-1)(# kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'ic_covarnce_ri_rs'
  longname = 'In-cloud Covariance of cloud ice mixing ratio and snow mixing ratio'
  units = '[(kg kg^-1)(kg kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'ic_covarnce_ri_Ns'
  longname = 'In-cloud Covariance of cloud ice mixing ratio and snowflake conc.'
  units = '[(kg kg^-1)(# kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  if (dograupel) then

    name = 'ic_covarnce_ri_rg'
    longname = 'In-cloud Covariance of cloud ice mixing ratio and graupel mixing ratio'
    units = '[(kg kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)

    name = 'ic_covarnce_ri_Ng'
    longname = 'In-cloud Covariance of cloud ice mixing ratio and graupel conc.'
    units = '[(kg kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
  end if !dograupel

!----------------------------
! Cloudice number concentration
!---------------------------
    
  name = 'ic_covarnce_Ni_rs'
  longname = 'In-cloud Covariance of cloud ice concentration and snow mixing ratio'
  units = '[(# kg^-1)(kg kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'ic_covarnce_Ni_Ns'
  longname = 'In-cloud Covariance of cloud ice concentration and snowflake conc.'
  units = '[(# kg^-1)(# kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  if (dograupel) then

    name = 'ic_covarnce_Ni_rg'
    longname = 'In-cloud Covariance of cloud ice conc. and graupel mixing ratio'
    units = '[(# kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)

    name = 'ic_covarnce_Ni_Ng'
    longname = 'In-cloud Covariance of cloud ice conc. and graupel conc.'
    units = '[(# kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
  end if !dograupel

!----------------------------
! Snow mixing ratio
!---------------------------
    
  name = 'ic_covarnce_rs_Ns'
  longname = 'In-cloud Covariance of snow mixing ratio and snowflake concentration'
  units = '[(kg kg^-1)(# kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  if (dograupel) then

    name = 'ic_covarnce_rs_rg'
    longname = 'In-cloud Covariance of snow mixing ratio and graupel mixing ratio'
    units = '[(kg kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)

    name = 'ic_covarnce_rs_Ng'
    longname = 'In-cloud Covariance of snow mixing ratio and graupel concentration'
    units = '[(kg kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)

!----------------------------
! Snow number concentration
!---------------------------
    
    name = 'ic_covarnce_Ns_rg'
    longname = 'In-cloud Covariance of snowflake concentration and graupel mixing ratio'
    units = '[(# kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)

    name = 'ic_covarnce_Ns_Ng'
    longname = 'In-cloud Covariance of snowflake concentration and graupel conc.'
    units = '[(# kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)

!----------------------------
! Graupel mixing ratio / number concentration
!---------------------------
    
    name = 'ic_covarnce_rg_Ng'
    longname = 'In-cloud Covariance of graupel mixing ratio and graupel concentration'
    units = '[(kg kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
  end if !dograupel
end if !doicemicro

!----------------------------
!  Out of cloud covariances
!---------------------------

!----------------------------
! chi(s_mellor)
!---------------------------
name = 'oc_covarnce_chi_w'
longname = 'In-cloud Covariance of chi and cloud droplet concentration'
units = '[(kg kg^-1)(m s^-1) ]'
call add_to_namelist(count,microcount,name,longname,units,0)

if(dopredictNc) then
  
  name = 'oc_covarnce_chi_Nc'
  longname = 'In-cloud Covariance of chi and cloud droplet concentration'
  units = '[(kg kg^-1)(#/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

end if

if (doprecip) then
    
  name = 'oc_covarnce_chi_rr'
  longname = 'In-cloud Covariance of chi and rain water mixing ratio'
  units = '[(kg kg^-1)(kg/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'oc_covarnce_chi_Nr'
  longname = 'In-cloud Covariance of chi and rain drop concentration'
  units = '[(kg kg^-1)(#/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

end if ! doprecip

if (doicemicro) then
    
  name = 'oc_covarnce_chi_ri'
  longname = 'In-cloud Covariance of chi and cloud ice mixing ratio'
  units = '[(kg kg^-1)(kg/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'oc_covarnce_chi_Ni'
  longname = 'In-cloud Covariance of chi and cloud ice concentration'
  units = '[(kg kg^-1)(#/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'oc_covarnce_chi_rs'
  longname = 'In-cloud Covariance of chi and snow mixing ratio'
  units = '[(kg kg^-1)(kg/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'oc_covarnce_chi_Ns'
  longname = 'In-cloud Covariance of chi and snowflake concentration'
  units = '[(kg kg^-1)(#/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

end if ! doicemicro

if (dograupel) then
    
  name = 'oc_covarnce_chi_rg'
  longname = 'In-cloud Covariance of chi and graupel mixing ratio'
  units = '[(kg kg^-1)(kg/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'oc_covarnce_chi_Ng'
  longname = 'In-cloud Covariance of chi and graupel concentration'
  units = '[(kg kg^-1)(#/kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

endif !doicemicro

!----------------------------
! Vertical velocity
!---------------------------
if(dopredictNc) then
  
  name = 'oc_covarnce_w_Nc'
  longname = 'In-cloud Covariance of vertical velocity and cloud droplet concentration'
  units = '[(m s^-1)(# kg^-1) ]'
  call add_to_namelist(count,microcount,name,longname,units,0)

end if

if (doprecip) then
    
  name = 'oc_covarnce_w_rr'
  longname = 'In-cloud Covariance of vertical velocity and rain water mixing ratio'
  units = '[(m s^-1)(kg kg^-1) ]'
  units = '-'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'oc_covarnce_w_Nr'
  longname = 'In-cloud Covariance of vertical velocity and rain drop concentration'
  units = '[(m s^-1)(# kg^-1) ]'
  call add_to_namelist(count,microcount,name,longname,units,0)

end if ! doprecip

if (doicemicro) then
    
  name = 'oc_covarnce_w_ri'
  longname = 'In-cloud Covariance of vertical velocity and cloud ice mixing ratio'
  units = '[(m s^-1)(kg kg^-1) ]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'oc_covarnce_w_Ni'
  longname = 'In-cloud Covariance of vertical velocity and cloud ice concentration'
  units = '[(m s^-1)(# kg^-1) ]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'oc_covarnce_w_rs'
  longname = 'In-cloud Covariance of vertical velocity and snow mixing ratio'
  units = '[(m s^-1)(kg kg^-1) ]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'oc_covarnce_w_Ns'
  longname = 'In-cloud Covariance of vertical velocity and snowflake concentration'
  units = '[(m s^-1)(# kg^-1) ]'
  call add_to_namelist(count,microcount,name,longname,units,0)

end if ! doicemicro

if (dograupel) then
    
  name = 'oc_covarnce_w_rg'
  longname = 'In-cloud Covariance of vertical velocity and graupel mixing ratio'
  units = '[(m s^-1)(kg kg^-1) ]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'oc_covarnce_w_Ng'
  longname = 'In-cloud Covariance of vertical velocity and graupel concentration'
  units = '[(m s^-1)(# kg^-1) ]'
  call add_to_namelist(count,microcount,name,longname,units,0)

end if ! dograupel

!----------------------------
! Cloud droplet number concentration
!---------------------------

if(dopredictNc) then
    if (doprecip) then
        
      name = 'oc_covarnce_Nc_rr'
      longname = 'In-cloud Covariance of cloud droplet conc. and rain water mixing ratio'
      units = '[(# kg^-1)(kg kg^-1)]'
      call add_to_namelist(count,microcount,name,longname,units,0)
        
      name = 'oc_covarnce_Nc_Nr'
      longname = 'In-cloud Covariance of cloud droplet conc. and rain drop conc.'
      units = '[(# kg^-1)(# kg^-1)]'
      call add_to_namelist(count,microcount,name,longname,units,0)
    
    end if ! doprecip
    
    if (doicemicro) then
        
      name = 'oc_covarnce_Nc_ri'
      longname = 'In-cloud Covariance of cloud droplet conc. and cloud ice mixing ratio'
      units = '[(# kg^-1)(kg kg^-1)]'
      call add_to_namelist(count,microcount,name,longname,units,0)
        
      name = 'oc_covarnce_Nc_Ni'
      longname = 'In-cloud Covariance of cloud droplet conc. and cloud ice conc.'
      units = '[(# kg^-1)(# kg^-1)]'
      call add_to_namelist(count,microcount,name,longname,units,0)
        
      name = 'oc_covarnce_Nc_rs'
      longname = 'In-cloud Covariance of cloud droplet conc. and snow mixing ratio'
      units = '[(# kg^-1)(kg kg^-1)]'
      call add_to_namelist(count,microcount,name,longname,units,0)
    
      name = 'oc_covarnce_Nc_Ns'
      longname = 'In-cloud Covariance of cloud droplet conc. and snowflake conc.'
      units = '[(# kg^-1)(# kg^-1)]'
      call add_to_namelist(count,microcount,name,longname,units,0)
    
    end if ! doicemicro
    
    if (dograupel) then
        
      name = 'oc_covarnce_Nc_rg'
      longname = 'In-cloud Covariance of cloud droplet conc. and graupel mixing ratio'
      units = '[(# kg^-1)(kg kg^-1)]'
      units = '-'
      call add_to_namelist(count,microcount,name,longname,units,0)
        
      name = 'oc_covarnce_Nc_Ng'
      longname = 'In-cloud Covariance of cloud droplet conc. and graupel conc.'
      units = '[(# kg^-1)(# kg^-1)]'
      call add_to_namelist(count,microcount,name,longname,units,0)

    end if ! dograupel
end if ! do predictNc

!----------------------------
! Rainwater mixing ratio
!---------------------------
if (doprecip) then

    name = 'oc_covarnce_rr_Nr'
    longname = 'In-cloud Covariance of rain water mixing ratio and rain drop concentration'
    units = '[(kg kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)

  if(doicemicro) then

    name = 'oc_covarnce_rr_ri'
    longname = 'In-cloud Covariance of rain water mixing ratio and cloud ice mixing ratio'
    units = '[(kg kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
        
    name = 'oc_covarnce_rr_Ni'
    longname = 'In-cloud Covariance of rain water mixing ratio and cloud ice conc.'
    units = '[(kg kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
    name = 'oc_covarnce_rr_rs'
    longname = 'In-cloud Covariance of rain water mixing ratio and snow mixing ratio'
    units = '[(kg kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
    name = 'oc_covarnce_rr_Ns'
    longname = 'In-cloud Covariance of rain water mixing ratio and snowflake conc.'
    units = '[(kg kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
  
  end if !doicemicro

  if(dograupel) then

    name = 'oc_covarnce_rr_rg'
    longname = 'In-cloud Covariance of rain water mixing ratio and graupel mixing ratio'
    units = '[(kg kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
        
    name = 'oc_covarnce_rr_Ng'
    longname = 'In-cloud Covariance of rain water mixing ratio and graupel conc.'
    units = '[(kg kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
  end if !dograupel

!----------------------------
! Rainwater number concentration
!---------------------------
  if(doicemicro) then
        
    name = 'oc_covarnce_Nr_ri'
    longname = 'In-cloud Covariance of rain drop conc. and cloud ice mixing ratio'
    units = '[(# kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
        
    name = 'oc_covarnce_Nr_Ni'
    longname = 'In-cloud Covariance of rain drop conc. and cloud ice conc.'
    units = '[(# kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
    name = 'oc_covarnce_Nr_rs'
    longname = 'In-cloud Covariance of rain drop conc. and snow mixing ratio'
    units = '[(# kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
    name = 'oc_covarnce_Nr_Ns'
    longname = 'In-cloud Covariance of rain drop conc. and snowflake conc.'
    units = '[(# kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
  
  end if !doicemicro

  if(dograupel) then

    name = 'oc_covarnce_Nr_rg'
    longname = 'In-cloud Covariance of rain drop conc. and graupel mixing ratio'
    units = '[(# kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
        
    name = 'oc_covarnce_Nr_Ng'
    longname = 'In-cloud Covariance of rain drop conc. and graupel conc.'
    units = '[(# kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
  end if !dograupel
end if !doprecip
!----------------------------
! Cloudice mixing ratio
!---------------------------

if (doicemicro) then
    
  name = 'oc_covarnce_ri_Ni'
  longname = 'In-cloud Covariance of cloud ice mixing ratio and cloud ice conc.'
  units = '[(kg kg^-1)(# kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'oc_covarnce_ri_rs'
  longname = 'In-cloud Covariance of cloud ice mixing ratio and snow mixing ratio'
  units = '[(kg kg^-1)(kg kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'oc_covarnce_ri_Ns'
  longname = 'In-cloud Covariance of cloud ice mixing ratio and snowflake conc.'
  units = '[(kg kg^-1)(# kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  if (dograupel) then

    name = 'oc_covarnce_ri_rg'
    longname = 'In-cloud Covariance of cloud ice mixing ratio and graupel mixing ratio'
    units = '[(kg kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)

    name = 'oc_covarnce_ri_Ng'
    longname = 'In-cloud Covariance of cloud ice mixing ratio and graupel conc.'
    units = '[(kg kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
    
  end if !dograupel

!----------------------------
! Cloudice number concentration
!---------------------------
    
  name = 'oc_covarnce_Ni_rs'
  longname = 'In-cloud Covariance of cloud ice concentration and snow mixing ratio'
  units = '[(# kg^-1)(kg kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  name = 'oc_covarnce_Ni_Ns'
  longname = 'In-cloud Covariance of cloud ice concentration and snowflake conc.'
  units = '[(# kg^-1)(# kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  if (dograupel) then

    name = 'oc_covarnce_Ni_rg'
    longname = 'In-cloud Covariance of cloud ice conc. and graupel mixing ratio'
    units = '[(# kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)

    name = 'oc_covarnce_Ni_Ng'
    longname = 'In-cloud Covariance of cloud ice conc. and graupel conc.'
    units = '[(# kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)
  end if !dograupel

!----------------------------
! Snow mixing ratio
!---------------------------
    
  name = 'oc_covarnce_rs_Ns'
  longname = 'In-cloud Covariance of snow mixing ratio and snowflake concentration'
  units = '[(kg kg^-1)(# kg^-1)]'
  call add_to_namelist(count,microcount,name,longname,units,0)

  if (dograupel) then

    name = 'oc_covarnce_rs_rg'
    longname = 'In-cloud Covariance of snow mixing ratio and graupel mixing ratio'
    units = '[(kg kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)

    name = 'oc_covarnce_rs_Ng'
    longname = 'In-cloud Covariance of snow mixing ratio and graupel concentration'
    units = '[(kg kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)

!----------------------------
! Snow number concentration
!---------------------------
    
    name = 'oc_covarnce_Ns_rg'
    longname = 'In-cloud Covariance of snowflake concentration and graupel mixing ratio'
    units = '[(# kg^-1)(kg kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)

    name = 'oc_covarnce_Ns_Ng'
    longname = 'In-cloud Covariance of snowflake concentration and graupel conc.'
    units = '[(# kg^-1)(# kg^-1)]'
    call add_to_namelist(count,microcount,name,longname,units,0)

!----------------------------
! Graupel mixing ratio / number concentration
!---------------------------
    
    name = 'oc_covarnce_rg_Ng'
    longname = 'In-cloud Covariance of graupel mixing ratio and graupel concentration'
    units = '[(kg kg^-1)(# kg^-1)]'
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
use hbuffer, only: hbuf_put, hbuf_put_level, hbuf_avg_put
#else
use hbuffer, only: hbuf_put
#endif /*UWM_STATS*/
use params, only : lcond

#ifdef UWM_STATS
use compute_correlation_module, only: corravg_count, sum_value_xy_domain,& 
                                      mean_xy_domain, variance_xy_domain,&
                                      mean_ip_xy_domain, variance_ip_xy_domain,&
                                      covariance_xy_domain, covariance_ip_xy_domain,&
                                      compute_correlation, undef_corr  

implicit none

!Microphysics correlations and covariances matricies
real, dimension(11,11,nzm) :: micro_correlations ! Magic numbers are the number
                                                 ! of variables in the correlation
!                                                  and covariance matricies
real, dimension(11,11,nzm) :: micro_covarnce 

!Microphysics correlations and covariances matricies in-cloud
real, dimension(11,11,nzm) :: ic_micro_correlations 
real, dimension(11,11,nzm) :: ic_micro_covarnce 

!Microphysics correlations and covariances matricies out of cloud
real, dimension(11,11,nzm) :: oc_micro_correlations 
real, dimension(11,11,nzm) :: oc_micro_covarnce 

! Mean and variance of chi(s_mellor)
real, dimension(nzm) :: domain_mean_chi 
real, dimension(nzm) :: domain_varnce_chi 

! Mean and variance of chi(s_mellor) in-cloud
real, dimension(nzm) :: ic_mean_chi 
real, dimension(nzm) :: ic_varnce_chi 

! Mean and variance of chi(s_mellor) out of cloud
real, dimension(nzm) :: oc_mean_chi 
real, dimension(nzm) :: oc_varnce_chi 

! Vertical velocity, interpolated, domain means, and variances
real, dimension(nx,ny,nzm) :: w_zt !vertical velocity interpolated on the scalar grid
real, dimension(nzm) :: domain_mean_w_zt !domain average of w_zt
real, dimension(nzm) :: domain_varnce_w_zt !domain variance of w_zt

real, dimension(nzm) :: ic_mean_w_zt !domain average of w_zt in-cloud
real, dimension(nzm) :: ic_varnce_w_zt !domain variance of w_zt in-cloud

real, dimension(nzm) :: oc_mean_w_zt !domain average of w_zt out of cloud
real, dimension(nzm) :: oc_varnce_w_zt !domain variance of w_zt out of cloud

!Microphysical domain-wide means and variances 
real, dimension(nmicro_fields,nzm) :: domain_mean_micro !domain averages of micro_fields
real, dimension(nmicro_fields,nzm) :: domain_varnce_micro !domain variance of micro_fields

!Microphysical within-'precip' means and variances 
real, dimension(nmicro_fields,nzm) :: ip_mean_micro ! within-'precip' averages of micro_fields
real, dimension(nmicro_fields,nzm) :: ip_varnce_micro ! within-'precip' of micro_fields

!Microphysical within-cloud means and variances 
real, dimension(nmicro_fields,nzm) :: ic_mean_micro ! within-cloud averages of micro_fields
real, dimension(nmicro_fields,nzm) :: ic_varnce_micro ! within-cloud of micro_fields

!Microphysical out of cloud means and variances 
real, dimension(nmicro_fields,nzm) :: oc_mean_micro ! within-cloud averages of micro_fields
real, dimension(nmicro_fields,nzm) :: oc_varnce_micro ! within-cloud of micro_fields

integer :: idx_s, idx_w, idx_Nc, idx_rr, idx_Nr, idx_ri, idx_Ni,&
           idx_rs, idx_Ns, idx_rg, idx_Ng, micro_indx_start, &
           offset ! Morrison has indices for Nc, rr, etc. This is the offset
                  ! between Morrison's and the correlation/covariance matrix indicies
   
!--------------------------------------------------------------------------------------------------
!Ensemble fractions 
!
!  These are used to see the effect of tolerance values on the in-'precip'
!  variances and fractions. The magic number 5 refers to the number of tolerance
!  values we test. Starting from SAM's defaut 1e-6 [kg/kg] for all micro.
!  species, we modulate the threshold +/- 2 orders of magnitude. 

!  Any variables that are reliant on these thresholds are defined in this
!  section

! Binary mask of in(1) and out(0) of micro. species. Add '1' for 'out of cloud'
! field
real, dimension(nx,ny,nzm,nfrac_fields,nfractions) :: micro_mask
real, dimension(nx,ny,nzm) :: out_cloud_mask

! Number of gridpoints within micro. species. Add '1' for 'out of cloud'
! field
integer, dimension(nzm,nfrac_fields,nfractions) :: micro_sum
integer, dimension(nzm) :: out_cloud_sum

! Micro. fractions. Add '1' for 'out of cloud'
! field
real, dimension(nzm, nfrac_fields, nfractions) :: micro_frac 

! Cloud, Rain, Ice, Snow, Graupel points as defined by the threshhold,
! (thresh_index [kg kg^-1])
integer :: rc_pts, rr_pts, ri_pts, rs_pts, rg_pts, thresh_index

real ::  curr_thresh ! Current threshold value 
integer :: order_of_magnitude ! 10

! Threshold value used to compute in-'precip' means/variance. Set to nfractions
integer :: thresh_out

! Name for each fraction written to disk. 
character*20 name
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

      !mstor(k, n) = SUM(micro_field(1:nx,1:ny,k,n))-mstor(k,n)
      mstor(k, n) = SUM(dble(micro_field(1:nx,1:ny,k,n)))-mstor(k,n) !MWSWong: cast to double
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
#ifdef UWM_STATS
!mtend already includes sedimentation
        (mtend(:,n)-stend(:,n)),mkoutputscale(n)*factor_xy*86400.)
#else
        mtend(:,n),mkoutputscale(n)*factor_xy*86400.)
#endif
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
! In micro arrays, indicies for hydrometeor
 rc_pts = 1 
 rr_pts = 2
 ri_pts = 3
 rs_pts = 4
 rg_pts = 5
 
 ! Compute in-precip means/variance with the final threshold
 thresh_out = nfractions
 order_of_magnitude = 10

! Current threshold value
curr_thresh = frac_threshold_init

do thresh_index=1,nfractions

  do k=1,nzm
    do j=1,ny
      do i=1,nx
        ! We only want cloud water, not cloud water + rain water.
        if ( cloudliq(i,j,k)  > curr_thresh) then
          micro_mask(i,j,k,rc_pts,thresh_index) = 1
          out_cloud_mask(i,j,k) = 0
        else
          micro_mask(i,j,k,rc_pts,thresh_index) = 0
          out_cloud_mask(i,j,k) = 1
        end if 
        
        if (doprecip) then
          if (micro_field(i,j,k,iqr) > curr_thresh) then
            micro_mask(i,j,k,rr_pts,thresh_index) = 1
          else
            micro_mask(i,j,k,rr_pts,thresh_index) = 0
          end if
          
          if (doicemicro) then
            if (micro_field(i,j,k,iqci) > curr_thresh) then
              micro_mask(i,j,k,ri_pts,thresh_index) = 1
            else
              micro_mask(i,j,k,ri_pts,thresh_index) = 0
            end if

            if (micro_field(i,j,k,iqs) > curr_thresh) then
              micro_mask(i,j,k,rs_pts,thresh_index) = 1
            else
              micro_mask(i,j,k,rs_pts,thresh_index) = 0
            end if
          
            if (dograupel) then
              if (micro_field(i,j,k,iqg) > curr_thresh) then
                micro_mask(i,j,k,rg_pts,thresh_index) = 1
              else
                micro_mask(i,j,k,rg_pts,thresh_index) = 0
              end if
            
            end if!do graupel  
          end if!do icemicro
        end if!doprecip
      
      end do !i
    end do !j
  end do !k
 
! Increase the threshold by a factor of 10
curr_thresh = curr_thresh * order_of_magnitude
end do !thresh_index

! Compute the fractions by summing and finding the mean of the binary mask above.
do n=1,nfrac_fields
  do m=1,nfractions
    micro_sum(:,n,m) = sum_value_xy_domain(nzm,micro_mask(1:nx,1:ny,1:nzm,n,m))
    micro_frac(:,n,m) = mean_xy_domain(nzm,micro_mask(1:nx,1:ny,1:nzm,n,m))
  end do
end do
! Take care of 'out of' cloud points
out_cloud_sum(:) = sum_value_xy_domain( nzm,out_cloud_mask(1:nx,1:ny,1:nzm) )

!Below is used to compute covariances and correlations of microphysical quantities
!=================================================================================

! These are the indicies in the correlation/covariance matrix for each of the
! respective variables. 

! If all the hydrometeors are present, the correlation matrix would look like:
!     |  chi  |  w  |  Nc  |  rr  |  Nr  |  ri  |  Ni  |  rs  |  Ns  |  rg  |  Ng  |
! chi |   1      #     #      #      #      #      #      #      #      #      #
! w   |          1     #      #      #      #      #      #      #      #      #
! Nc  |                1      #      #      #      #      #      #      #      #
! rr  |                       1      #      #      #      #      #      #      #
! Nr  |                              1      #      #      #      #      #      #
! ri  |                                     1      #      #      #      #      #
! Ni  |                                            1      #      #      #      #
! rs  |                                                   1      #      #      #
! Ns  |                                                          1      #      #
! rg  |                                                                 1      #
! Ng  |                                                                        1

! SAM does not calculate the corr(chi,eta) since it would be difficult to compare
! CLUBB's 'two-component' correlation of chi,eta with SAM's domain-wide
! correlation.


idx_s  = 1  ! chi
idx_w  = 2  ! Vertical velocity
if(dopredictNc) then
  idx_Nc = 3  
  idx_rr = 4  
  idx_Nr = 5  
  idx_ri = 6 
  idx_Ni = 7 
  idx_rs = 8  
  idx_Ns = 9  
  idx_rg = 10 
  idx_Ng = 11 
  ! When computing means, variances, covariances, and correlations, we loop
  ! through Morrison's indicies (iqr) and write to our matricies (idx_rr)
  ! These two are related by the 'offset', weberjk(UWM)
  offset = idx_Nc - incl
  micro_indx_start = incl
else
  idx_rr = 3  
  idx_Nr = 4  
  idx_ri = 5 
  idx_Ni = 6 
  idx_rs = 7  
  idx_Ns = 8  
  idx_rg = 9 
  idx_Ng = 10 
  offset = idx_rr - iqr
  micro_indx_start = iqr
endif
  !Interpolate w to the scalar grid
  do k=1,nz-1
    do i=1,nx
      do j = 1,ny
        w_zt(i,j,k) =  LIN_INT( w(i,j,k+1), w(i,j,k), zi(k+1), zi(k), z(k) ) 
      end do
    end do
  end do

  !Find the domain wide means and variances of chi and w_zt
  domain_mean_chi(:) = mean_xy_domain(nzm, chi)
  domain_varnce_chi(:) = variance_xy_domain(nzm, chi, domain_mean_chi)
  
  domain_mean_w_zt(:) = mean_xy_domain(nzm, w_zt)
  domain_varnce_w_zt(:) = variance_xy_domain(nzm, w_zt, domain_mean_w_zt)

  !Find the in-cloud means and variances of chi and w_zt
  ic_mean_chi(:) = mean_ip_xy_domain( nzm, chi, micro_mask(1:nx,1:ny,1:nzm,rc_pts,thresh_out),&
                                   micro_sum(1:nzm,rc_pts,thresh_out) )
  ic_varnce_chi(:) = variance_ip_xy_domain( nzm, chi, micro_mask(1:nx,1:ny,1:nzm,rc_pts,thresh_out),&
                                            ic_mean_chi(:), micro_sum(1:nzm,rc_pts,thresh_out) )

  ic_mean_w_zt(:) = mean_ip_xy_domain( nzm, w_zt, micro_mask(1:nx,1:ny,1:nzm,rc_pts,thresh_out),&
                                   micro_sum(1:nzm,rc_pts,thresh_out) )
  ic_varnce_w_zt(:) = variance_ip_xy_domain( nzm, w_zt, micro_mask(1:nx,1:ny,1:nzm,rc_pts,thresh_out),&
                                            ic_mean_w_zt(:), micro_sum(1:nzm,rc_pts,thresh_out) )

  !Find the out of cloud means and variances of chi and w_zt
  oc_mean_chi(:) = mean_ip_xy_domain( nzm, chi, out_cloud_mask(1:nx,1:ny,1:nzm),&
                                   out_cloud_sum(1:nzm) )
  oc_varnce_chi(:) = variance_ip_xy_domain( nzm, chi, out_cloud_mask(1:nx,1:ny,1:nzm),&
                                            oc_mean_chi(:), out_cloud_sum(1:nzm) )

  oc_mean_w_zt(:) = mean_ip_xy_domain( nzm, w_zt, out_cloud_mask(1:nx,1:ny,1:nzm),&
                                   out_cloud_sum(1:nzm) )
  oc_varnce_w_zt(:) = variance_ip_xy_domain( nzm, w_zt, out_cloud_mask(1:nx,1:ny,1:nzm),&
                                            oc_mean_w_zt(:), out_cloud_sum(1:nzm) )
  
  do n=micro_indx_start,nmicro_fields
    ! Domain mean and variances stored according to Morrison indicies
    domain_mean_micro(n,:) = mean_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,n))
    domain_varnce_micro(n,:) = variance_xy_domain(nzm, micro_field(1:nx,1:ny,1:nzm,n),&
                                                domain_mean_micro(n,:))

    ! in-cloud means and variances
    ic_mean_micro(n,:) = mean_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,n),&
                                          micro_mask(1:nx,1:ny,1:nzm,rc_pts,thresh_out),micro_sum(1:nzm,rc_pts,thresh_out))

    ic_varnce_micro(n,:) = variance_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,n),&
                                          micro_mask(1:nx,1:ny,1:nzm,rc_pts,thresh_out),ic_mean_micro(n,:),&
                                          micro_sum(1:nzm,rc_pts,thresh_out))

    ! out of cloud means and variances
    oc_mean_micro(n,:) = mean_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,n),&
                                          out_cloud_mask(1:nx,1:ny,1:nzm),out_cloud_sum(1:nzm) )

    oc_varnce_micro(n,:) = variance_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,n),&
                                          out_cloud_mask(1:nx,1:ny,1:nzm),oc_mean_micro(n,:),&
                                          out_cloud_sum(1:nzm) )
  end do
  !----------------
  ! Within precip means
  !----------------
if(dopredictNc) then
  ip_mean_micro(incl,:) = mean_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,incl),&
                                          micro_mask(1:nx,1:ny,1:nzm,rc_pts,thresh_out),micro_sum(1:nzm,rc_pts,thresh_out))
endif

  if(doprecip) then
     
    ip_mean_micro(iqr,:) = mean_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,iqr),&
                                          micro_mask(1:nx,1:ny,1:nzm,rr_pts,thresh_out),micro_sum(1:nzm,rr_pts,thresh_out))
  
    ip_mean_micro(inr,:) = mean_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,inr),&
                                          micro_mask(1:nx,1:ny,1:nzm,rr_pts,thresh_out),micro_sum(1:nzm,rr_pts,thresh_out))
    if(doicemicro) then  
      
      ip_mean_micro(iqci,:) = mean_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,iqci),&
                                          micro_mask(1:nx,1:ny,1:nzm,ri_pts,thresh_out),micro_sum(1:nzm,ri_pts,thresh_out))
  
      ip_mean_micro(inci,:) = mean_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,inci),&
                                          micro_mask(1:nx,1:ny,1:nzm,ri_pts,thresh_out),micro_sum(1:nzm,ri_pts,thresh_out))
      
      ip_mean_micro(iqs,:) = mean_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,iqs),&
                                          micro_mask(1:nx,1:ny,1:nzm,rs_pts,thresh_out),micro_sum(1:nzm,rs_pts,thresh_out))
  
      ip_mean_micro(ins,:) = mean_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,ins),&
                                          micro_mask(1:nx,1:ny,1:nzm,rs_pts,thresh_out),micro_sum(1:nzm,rs_pts,thresh_out))
      if(dograupel) then
        
        ip_mean_micro(iqg,:) = mean_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,iqg),&
                                          micro_mask(1:nx,1:ny,1:nzm,rg_pts,thresh_out),micro_sum(1:nzm,rg_pts,thresh_out))
  
        ip_mean_micro(ing,:) = mean_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,ing),&
                                          micro_mask(1:nx,1:ny,1:nzm,rg_pts,thresh_out),micro_sum(1:nzm,rg_pts,thresh_out))
      endif !dograupel
    endif !doice
  endif !doprecip
  
  !----------------
  ! Within precip variances
  !----------------
if(dopredictNc) then
  ip_varnce_micro(incl,:) = variance_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,incl),&
                                          micro_mask(1:nx,1:ny,1:nzm,1,thresh_out),ip_mean_micro(incl,:),&
                                          micro_sum(1:nzm,rc_pts,thresh_out))
endif

  if(doprecip) then
    ! Rain 
    ip_varnce_micro(iqr,:) = variance_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,iqr),&
                                          micro_mask(1:nx,1:ny,1:nzm,2,thresh_out),ip_mean_micro(iqr,:),&
                                          micro_sum(1:nzm,rr_pts,thresh_out))
  
    ip_varnce_micro(inr,:) = variance_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,inr),&
                                          micro_mask(1:nx,1:ny,1:nzm,2,thresh_out),ip_mean_micro(inr,:),&
                                          micro_sum(1:nzm,rr_pts,thresh_out))
    if(doicemicro) then  
      ! Cloud ice
      ip_varnce_micro(iqci,:) = variance_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,iqci),&
                                          micro_mask(1:nx,1:ny,1:nzm,3,thresh_out),ip_mean_micro(iqci,:),&
                                          micro_sum(1:nzm,ri_pts,thresh_out))
  
      ip_varnce_micro(inci,:) = variance_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,inci),&
                                          micro_mask(1:nx,1:ny,1:nzm,3,thresh_out),ip_mean_micro(inci,:),&
                                          micro_sum(1:nzm,ri_pts,thresh_out))
      ! Snow 
      ip_varnce_micro(iqs,:) = variance_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,iqs),&
                                          micro_mask(1:nx,1:ny,1:nzm,rs_pts,thresh_out),ip_mean_micro(iqs,:),&
                                          micro_sum(1:nzm,rs_pts,thresh_out))
  
      ip_varnce_micro(ins,:) = variance_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,ins),&
                                          micro_mask(1:nx,1:ny,1:nzm,rs_pts,thresh_out),ip_mean_micro(ins,:),&
                                          micro_sum(1:nzm,rs_pts,thresh_out))
      if(dograupel) then
      ! Graupel 
        ip_varnce_micro(iqg,:) = variance_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,iqg),&
                                          micro_mask(1:nx,1:ny,1:nzm,rg_pts,thresh_out),ip_mean_micro(iqg,:),&
                                          micro_sum(1:nzm,rg_pts,thresh_out))
  
        ip_varnce_micro(ing,:) = variance_ip_xy_domain(nzm,micro_field(1:nx,1:ny,1:nzm,ing),&
                                          micro_mask(1:nx,1:ny,1:nzm,rg_pts,thresh_out),ip_mean_micro(ing,:),&
                                          micro_sum(1:nzm,rg_pts,thresh_out))
      endif !dograupel
    endif !doice
  endif !doprecip
 
!-----------------------------------------
! Domain-wide covariances and correlations
!-----------------------------------------

  ! Compute covariances and correlations, first for chi(s_mellor). Manually for
  ! vertical velocity. 
  micro_covarnce(idx_s,idx_w,:) = covariance_xy_domain(nzm,&
                                    chi(1:nx,1:ny,1:nzm), w_zt(1:nx,1:ny,1:nzm),&
                                    domain_mean_chi(:), domain_mean_w_zt(:) )

  micro_correlations(idx_s,idx_w,:) = compute_correlation(nzm, domain_varnce_chi(:),&
                                        domain_varnce_w_zt(:), micro_covarnce(idx_s,idx_w,:) )
  

  ! Compute covariances and correlations, first for chi(s_mellor). Now, loop through
  ! all hydrometeor fields using Morrison indicies, but write in our
  ! covarnce/corr matrix by offset.
  do n = micro_indx_start,nmicro_fields 
      micro_covarnce(idx_s,n+offset,:) = covariance_xy_domain(nzm, chi(1:nx,1:ny,1:nzm), &
                                    micro_field(1:nx,1:ny,1:nzm,n), domain_mean_chi(:), &
                                    domain_mean_micro(n,:) )

      micro_correlations(idx_s,n+offset,:) = compute_correlation(nzm, domain_varnce_chi(:),&
                                  domain_varnce_micro(n,:), micro_covarnce(idx_s,n+offset,:) )
!                                                                                ^ in covarnce matrix
!                                                                                Morrison index is offset
  end do
  
  ! Compute covariances and correlations, now for vertical velocity and all
  ! hydrometeor variables.
  do n = micro_indx_start,nmicro_fields 
      micro_covarnce(idx_w,n+offset,:) = covariance_xy_domain(nzm, w_zt, micro_field(1:nx,1:ny,1:nzm,n),&
                              domain_mean_w_zt, domain_mean_micro(n,:) )

      micro_correlations(idx_w,n+offset,:) = compute_correlation(nzm, domain_varnce_w_zt,&
                                  domain_varnce_micro(n,:), micro_covarnce(idx_w,n+offset,:) )
  end do

  ! Continue to fill the covariance and correlation matricies  
  do m = micro_indx_start, nmicro_fields 
    do n = m, nmicro_fields
  
      micro_covarnce(m+offset,n+offset,:) = covariance_xy_domain(nzm, micro_field(1:nx,1:ny,1:nzm,m), &
          micro_field(1:nx,1:ny,1:nzm,n),domain_mean_micro(m,:), domain_mean_micro(n,:) )

      micro_correlations(m+offset,n+offset,:) = compute_correlation(nzm, domain_varnce_micro(m,:),&
                                  domain_varnce_micro(n,:), micro_covarnce(m+offset,n+offset,:) )
    end do
  end do

!-----------------------------------------
! In-cloud covariances and correlations
!-----------------------------------------

  ic_micro_covarnce(idx_s,idx_w,:) = covariance_ip_xy_domain(nzm,&
                                    chi(1:nx,1:ny,1:nzm), w_zt(1:nx,1:ny,1:nzm),&
                                    micro_mask(1:nx,1:ny,1:nzm,rc_pts,thresh_out),&
                                    ic_mean_chi(:), ic_mean_w_zt(:),&
                                    micro_sum(:,rc_pts,thresh_out) )

  ic_micro_correlations(idx_s,idx_w,:) = compute_correlation(nzm, ic_varnce_chi(:),&
                                        ic_varnce_w_zt(:), ic_micro_covarnce(idx_s,idx_w,:) )
  

  do n = micro_indx_start,nmicro_fields 
      ic_micro_covarnce(idx_s,n+offset,:) = covariance_ip_xy_domain(nzm,&
                                             chi(1:nx,1:ny,1:nzm), micro_field(1:nx,1:ny,1:nzm,n),&
                                             micro_mask(1:nx,1:ny,1:nzm,rc_pts,thresh_out),&
                                             ic_mean_chi(:), ic_mean_micro(n,:),&
                                             micro_sum(:,rc_pts,thresh_out) )

      ic_micro_correlations(idx_s,n+offset,:) = compute_correlation(nzm, ic_varnce_chi(:),&
                                        ic_varnce_micro(n,:), ic_micro_covarnce(idx_s,n+offset,:) )
  end do
  
  ! Compute covariances and correlations, now for vertical velocity and all
  ! hydrometeor variables.
  do n = micro_indx_start,nmicro_fields 
      ic_micro_covarnce(idx_w,n+offset,:) = covariance_ip_xy_domain(nzm,&
                                            w_zt, micro_field(1:nx,1:ny,1:nzm,n),&
                                            micro_mask(1:nx,1:ny,1:nzm,rc_pts,thresh_out),&
                                            ic_mean_w_zt, ic_mean_micro(n,:), &
                                            micro_sum(:,rc_pts,thresh_out) )

      ic_micro_correlations(idx_w,n+offset,:) = compute_correlation(nzm, ic_varnce_w_zt,&
                                  ic_varnce_micro(n,:), ic_micro_covarnce(idx_w,n+offset,:) )
  end do

  ! Continue to fill the covariance and correlation matricies  
  do m = micro_indx_start, nmicro_fields 
    do n = m, nmicro_fields
  
      ic_micro_covarnce(m+offset,n+offset,:) = covariance_ip_xy_domain(nzm,&
                                 micro_field(1:nx,1:ny,1:nzm,m), micro_field(1:nx,1:ny,1:nzm,n),&
                                 micro_mask(1:nx,1:ny,1:nzm,rc_pts,thresh_out),&
                                 ic_mean_micro(m,:), ic_mean_micro(n,:),&
                                 micro_sum(:,rc_pts,thresh_out) )

      ic_micro_correlations(m+offset,n+offset,:) = compute_correlation(nzm, ic_varnce_micro(m,:),&
                                  ic_varnce_micro(n,:), ic_micro_covarnce(m+offset,n+offset,:) )
    end do
  end do

!-----------------------------------------
! Out of cloud covariances and correlations
!-----------------------------------------

  oc_micro_covarnce(idx_s,idx_w,:) = covariance_ip_xy_domain(nzm,&
                                    chi(1:nx,1:ny,1:nzm), w_zt(1:nx,1:ny,1:nzm),&
                                    out_cloud_mask(1:nx,1:ny,1:nzm),&
                                    oc_mean_chi(:), oc_mean_w_zt(:),&
                                    out_cloud_sum(:) )

  oc_micro_correlations(idx_s,idx_w,:) = compute_correlation(nzm, oc_varnce_chi(:),&
                                        oc_varnce_w_zt(:), oc_micro_covarnce(idx_s,idx_w,:) )
  

  do n = micro_indx_start,nmicro_fields 
      oc_micro_covarnce(idx_s,n+offset,:) = covariance_ip_xy_domain(nzm,&
                                             chi(1:nx,1:ny,1:nzm), micro_field(1:nx,1:ny,1:nzm,n),&
                                             out_cloud_mask(1:nx,1:ny,1:nzm),&
                                             oc_mean_chi(:), oc_mean_micro(n,:),&
                                             out_cloud_sum(:) )

      oc_micro_correlations(idx_s,n+offset,:) = compute_correlation(nzm, oc_varnce_chi(:),&
                                        oc_varnce_micro(n,:), oc_micro_covarnce(idx_s,n+offset,:) )
  end do
  
  ! Compute covariances and correlations, now for vertical velocity and all
  ! hydrometeor variables.
  do n = micro_indx_start,nmicro_fields 
      oc_micro_covarnce(idx_w,n+offset,:) = covariance_ip_xy_domain(nzm,&
                                            w_zt, micro_field(1:nx,1:ny,1:nzm,n),&
                                            out_cloud_mask(1:nx,1:ny,1:nzm),&
                                            oc_mean_w_zt, oc_mean_micro(n,:), &
                                            out_cloud_sum(:) )

      oc_micro_correlations(idx_w,n+offset,:) = compute_correlation(nzm, oc_varnce_w_zt,&
                                  oc_varnce_micro(n,:), oc_micro_covarnce(idx_w,n+offset,:) )
  end do

  ! Continue to fill the covariance and correlation matricies  
  do m = micro_indx_start, nmicro_fields 
    do n = m, nmicro_fields
  
      oc_micro_covarnce(m+offset,n+offset,:) = covariance_ip_xy_domain(nzm,&
                                 micro_field(1:nx,1:ny,1:nzm,m), micro_field(1:nx,1:ny,1:nzm,n),&
                                 out_cloud_mask(1:nx,1:ny,1:nzm),&
                                 oc_mean_micro(m,:), oc_mean_micro(n,:),&
                                 out_cloud_sum(:) )

      oc_micro_correlations(m+offset,n+offset,:) = compute_correlation(nzm, oc_varnce_micro(m,:),&
                                  oc_varnce_micro(n,:), oc_micro_covarnce(m+offset,n+offset,:) )
    end do
  end do
  
  
#endif /* UWM_STATS */

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
call hbuf_avg_put('hl_on_Cp_res',t_res,dimx1_s,dimx2_s,dimy1_s,dimy2_s,nzm,1.)

! Write out each hydrometeor fraction
do n = 1,nfrac_fields
  do m = 1,nfractions
    name = trim(frac_name(n))//'_'//trim(adjustl(frac_in_char(m)))

    call hbuf_put(name,micro_frac(:,n,m),1.)

  end do
end do

!----------------------
! Domain-wide mean
!----------------------
  call hbuf_put('chim', domain_mean_chi(:), 1.)
if(dopredictNc) then
  call hbuf_put('Ncm', domain_mean_micro(incl,:), 1.)
endif
!----------------------
! Domain-wide variance
!----------------------
  call hbuf_put('chip2', domain_varnce_chi(:), 1.)

if(dopredictNc) then
  call hbuf_put('Ncp2', domain_varnce_micro(incl,:), 1.)
endif

  if (doprecip) then
    call hbuf_put('Nrp2', domain_varnce_micro(inr,:), 1.)
    call hbuf_put('rrp2', domain_varnce_micro(iqr,:), 1.)
    
    if (doicemicro) then
      call hbuf_put('Nip2', domain_varnce_micro(inci,:), 1.)
      call hbuf_put('rip2', domain_varnce_micro(iqci,:), 1.)
      
      call hbuf_put('Nsp2', domain_varnce_micro(ins,:), 1.)
      call hbuf_put('rsp2', domain_varnce_micro(iqs,:), 1.)
      
      if (dograupel) then  
        call hbuf_put('Ngp2', domain_varnce_micro(ing,:), 1.)
        call hbuf_put('rgp2', domain_varnce_micro(iqg,:), 1.)
      end if
    end if 
  endif

!----------------------
! Within-precip mean
!----------------------
if(dopredictNc) then
  call hbuf_put('Ncm_ip', ip_mean_micro(incl,:), 1.)
endif

  if (doprecip) then
    call hbuf_put('Nrm_ip', ip_mean_micro(inr,:), 1.)
    call hbuf_put('rrm_ip', ip_mean_micro(iqr,:), 1.)
    
    if (doicemicro) then
      call hbuf_put('Nim_ip', ip_mean_micro(inci,:), 1.)
      call hbuf_put('rim_ip', ip_mean_micro(iqci,:), 1.)
      
      call hbuf_put('Nsm_ip', ip_mean_micro(ins,:), 1.)
      call hbuf_put('rsm_ip', ip_mean_micro(iqs,:), 1.)
      
      if (dograupel) then  
        call hbuf_put('Ngm_ip', ip_mean_micro(ing,:), 1.)
        call hbuf_put('rgm_ip', ip_mean_micro(iqg,:), 1.)
      end if
    end if 
  endif


!----------------------
! Within-precip variance
!----------------------
if(dopredictNc) then
  call hbuf_put('Ncp2_ip', ip_varnce_micro(incl,:), 1.)
endif

  if (doprecip) then
    call hbuf_put('Nrp2_ip', ip_varnce_micro(inr,:), 1.)
    call hbuf_put('rrp2_ip', ip_varnce_micro(iqr,:), 1.)
    
    if (doicemicro) then
      call hbuf_put('Nip2_ip', ip_varnce_micro(inci,:), 1.)
      call hbuf_put('rip2_ip', ip_varnce_micro(iqci,:), 1.)
      
      call hbuf_put('Nsp2_ip', ip_varnce_micro(ins,:), 1.)
      call hbuf_put('rsp2_ip', ip_varnce_micro(iqs,:), 1.)
      
      if (dograupel) then  
        call hbuf_put('Ngp2_ip', ip_varnce_micro(ing,:), 1.)
        call hbuf_put('rgp2_ip', ip_varnce_micro(iqg,:), 1.)
      end if
    end if 
  endif

! weberjk(UWM), Microphysical correlations

   do k = 1, nzm, 1
!-------------
! chi
!-------------
      if ( micro_correlations(idx_s,idx_w,k) /= undef_corr ) then
         call hbuf_put_level('corr_chi_w',micro_correlations(idx_s,idx_w,k),1.,k)
         corravg_count(idx_cor_chi_w,k) = corravg_count(idx_cor_chi_w,k) + 1
      else ! micro_correlations(1,incl,k) == undef_corr
         call hbuf_put_level('corr_chi_w',0.0,1.,k)
      endif ! micro_correlations(1,incl,k) /= undef_corr
      
      if(dopredictNc) then
        if ( micro_correlations(idx_s,idx_Nc,k) /= undef_corr ) then
           call hbuf_put_level('corr_chi_Nc',micro_correlations(idx_s,idx_Nc,k),1.,k)
           corravg_count(idx_cor_chi_Nc,k) = corravg_count(idx_cor_chi_Nc,k) + 1
        else ! micro_correlations(1,incl,k) == undef_corr
           call hbuf_put_level('corr_chi_Nc',0.0,1.,k)
        endif ! micro_correlations(1,incl,k) /= undef_corr
      end if

      if (doprecip) then
      if ( micro_correlations(idx_s,idx_rr,k) /= undef_corr ) then
           call hbuf_put_level('corr_chi_rr',micro_correlations(idx_s,idx_rr,k),1.,k)
           corravg_count(idx_cor_chi_rr,k) = corravg_count(idx_cor_chi_rr,k) + 1
        else ! micro_correlations(1,iqr,k) == undef_corr
           call hbuf_put_level('corr_chi_rr',0.0,1.,k)
        endif ! micro_correlations(1,iqr,k) /= undef_corr
        if ( micro_correlations(idx_s,idx_Nr,k) /= undef_corr ) then
           call hbuf_put_level('corr_chi_Nr',micro_correlations(idx_s,idx_Nr,k),1.,k)
           corravg_count(idx_cor_chi_Nr,k) = corravg_count(idx_cor_chi_Nr,k) + 1
        else ! micro_correlations(1,inr,k) == undef_corr
           call hbuf_put_level('corr_chi_Nr',0.0,1.,k)
        endif ! ! micro_correlations(1,inr,k) /= undef_corr
      end if ! doprecip

      if (doicemicro) then
      if ( micro_correlations(idx_s,idx_ri,k) /= undef_corr ) then
           call hbuf_put_level('corr_chi_ri',micro_correlations(idx_s,idx_ri,k),1.,k)
           corravg_count(idx_cor_chi_ri,k) = corravg_count(idx_cor_chi_ri,k) + 1
        else ! micro_correlations(1,iqci,k) == undef_corr
           call hbuf_put_level('corr_chi_ri',0.0,1.,k)
        endif ! micro_correlations(1,iqci,k) /= undef_corr
        if ( micro_correlations(idx_s,idx_Ni,k) /= undef_corr ) then
           call hbuf_put_level('corr_chi_Ni',micro_correlations(idx_s,idx_Ni,k),1.,k)
           corravg_count(idx_cor_chi_Ni,k) = corravg_count(idx_cor_chi_Ni,k) + 1
        else ! micro_correlations(1,inci,k) == undef_corr
           call hbuf_put_level('corr_chi_Ni',0.0,1.,k)
        endif ! micro_correlations(1,inci,k) /= undef_corr
        if ( micro_correlations(idx_s,idx_rs,k) /= undef_corr ) then
           call hbuf_put_level('corr_chi_rs',micro_correlations(idx_s,idx_rs,k),1.,k)
           corravg_count(idx_cor_chi_rs,k) = corravg_count(idx_cor_chi_rs,k) + 1
        else ! micro_correlations(1,iqs,k) == undef_corr
           call hbuf_put_level('corr_chi_rs',0.0,1.,k)
        endif ! micro_correlations(1,iqs,k) /= undef_corr
        if ( micro_correlations(idx_s,idx_Ns,k) /= undef_corr ) then
           call hbuf_put_level('corr_chi_Ns',micro_correlations(idx_s,idx_Ns,k),1.,k)
           corravg_count(idx_cor_chi_Ns,k) = corravg_count(idx_cor_chi_Ns,k) + 1
        else ! micro_correlations(1,ins,k) == undef_corr
           call hbuf_put_level('corr_chi_Ns',0.0,1.,k)
        endif ! micro_correlations(1,ins,k) /= undef_corr
      end if ! doicemicro

      if (dograupel) then
        if ( micro_correlations(idx_s,idx_rg,k) /= undef_corr ) then
           call hbuf_put_level('corr_chi_rg',micro_correlations(idx_s,idx_rg,k),1.,k)
           corravg_count(idx_cor_chi_rg,k) = corravg_count(idx_cor_chi_rg,k) + 1
        else ! micro_correlations(1,iqg,k) == undef_corr
           call hbuf_put_level('corr_chi_rg',0.0,1.,k)
        endif ! micro_correlations(1,iqg,k) /= undef_corr
        if ( micro_correlations(idx_s,idx_Ng,k) /= undef_corr ) then
           call hbuf_put_level('corr_chi_Ng',micro_correlations(idx_s,idx_Ng,k),1.,k)
           corravg_count(idx_cor_chi_Ng,k) = corravg_count(idx_cor_chi_Ng,k) + 1
        else ! micro_correlations(1,ing,k) == undef_corr
           call hbuf_put_level('corr_chi_Ng',0.0,1.,k)
        endif ! micro_correlations(1,ing,k) /= undef_corr
      end if ! dograupel

!-------------
! Vertical Velocity
!-------------
      if(dopredictNc) then
        if ( micro_correlations(idx_w,idx_Nc,k) /= undef_corr ) then
           call hbuf_put_level('corr_w_Nc',micro_correlations(idx_w,idx_Nc,k),1.,k)
           corravg_count(idx_cor_w_Nc,k) = corravg_count(idx_cor_w_Nc,k) + 1
        else ! micro_correlations(1,incl,k) == undef_corr
           call hbuf_put_level('corr_w_Nc',0.0,1.,k)
        endif ! micro_correlations(1,incl,k) /= undef_corr
      end if

      if (doprecip) then
      if ( micro_correlations(idx_w,idx_rr,k) /= undef_corr ) then
           call hbuf_put_level('corr_w_rr',micro_correlations(idx_w,idx_rr,k),1.,k)
           corravg_count(idx_cor_w_rr,k) = corravg_count(idx_cor_w_rr,k) + 1
        else ! micro_correlations(1,iqr,k) == undef_corr
           call hbuf_put_level('corr_w_rr',0.0,1.,k)
        endif ! micro_correlations(1,iqr,k) /= undef_corr
        if ( micro_correlations(idx_w,idx_Nr,k) /= undef_corr ) then
           call hbuf_put_level('corr_w_Nr',micro_correlations(idx_w,idx_Nr,k),1.,k)
           corravg_count(idx_cor_w_nr,k) = corravg_count(idx_cor_w_nr,k) + 1
        else ! micro_correlations(1,inr,k) == undef_corr
           call hbuf_put_level('corr_w_Nr',0.0,1.,k)
        endif ! ! micro_correlations(1,inr,k) /= undef_corr
      end if ! doprecip

      if (doicemicro) then
      if ( micro_correlations(idx_w,idx_ri,k) /= undef_corr ) then
           call hbuf_put_level('corr_w_ri',micro_correlations(idx_w,idx_ri,k),1.,k)
           corravg_count(idx_cor_w_ri,k) = corravg_count(idx_cor_w_ri,k) + 1
        else ! micro_correlations(1,iqci,k) == undef_corr
           call hbuf_put_level('corr_w_ri',0.0,1.,k)
        endif ! micro_correlations(1,iqci,k) /= undef_corr
        if ( micro_correlations(idx_w,idx_Ni,k) /= undef_corr ) then
           call hbuf_put_level('corr_w_Ni',micro_correlations(idx_w,idx_Ni,k),1.,k)
           corravg_count(idx_cor_w_Ni,k) = corravg_count(idx_cor_w_Ni,k) + 1
        else ! micro_correlations(1,inci,k) == undef_corr
           call hbuf_put_level('corr_w_Ni',0.0,1.,k)
        endif ! micro_correlations(1,inci,k) /= undef_corr
        if ( micro_correlations(idx_w,idx_rs,k) /= undef_corr ) then
           call hbuf_put_level('corr_w_rs',micro_correlations(idx_w,idx_rs,k),1.,k)
           corravg_count(idx_cor_w_rs,k) = corravg_count(idx_cor_w_rs,k) + 1
        else ! micro_correlations(1,iqs,k) == undef_corr
           call hbuf_put_level('corr_w_rs',0.0,1.,k)
        endif ! micro_correlations(1,iqs,k) /= undef_corr
        if ( micro_correlations(idx_w,idx_Ns,k) /= undef_corr ) then
           call hbuf_put_level('corr_w_Ns',micro_correlations(idx_w,idx_Ns,k),1.,k)
           corravg_count(idx_cor_w_Ns,k) = corravg_count(idx_cor_w_Ns,k) + 1
        else ! micro_correlations(1,ins,k) == undef_corr
           call hbuf_put_level('corr_w_Ns',0.0,1.,k)
        endif ! micro_correlations(1,ins,k) /= undef_corr
      end if ! doicemicro

      if (dograupel) then
        if ( micro_correlations(idx_w,idx_rg,k) /= undef_corr ) then
           call hbuf_put_level('corr_w_rg',micro_correlations(idx_w,idx_rg,k),1.,k)
           corravg_count(idx_cor_w_rg,k) = corravg_count(idx_cor_w_rg,k) + 1
        else ! micro_correlations(1,iqg,k) == undef_corr
           call hbuf_put_level('corr_w_rg',0.0,1.,k)
        endif ! micro_correlations(1,iqg,k) /= undef_corr
        if ( micro_correlations(idx_w,idx_Ng,k) /= undef_corr ) then
           call hbuf_put_level('corr_w_Ng',micro_correlations(idx_w,idx_Ng,k),1.,k)
           corravg_count(idx_cor_w_Ng,k) = corravg_count(idx_cor_w_Ng,k) + 1
        else ! micro_correlations(1,ing,k) == undef_corr
           call hbuf_put_level('corr_w_Ng',0.0,1.,k)
        endif ! micro_correlations(1,ing,k) /= undef_corr
      end if ! dograupel

!-------------
! Cloud droplet number concentration
!-------------
      if (dopredictNc) then
     
        if (doprecip) then
           if ( micro_correlations(idx_Nc,idx_rr,k) /= undef_corr ) then
              call hbuf_put_level('corr_Nc_rr',micro_correlations(idx_Nc,idx_rr,k),1.,k)
              corravg_count(idx_cor_Nc_rr,k) = corravg_count(idx_cor_Nc_rr,k) + 1
           else ! micro_correlations(incl,iqr,k) == undef_corr
              call hbuf_put_level('corr_Nc_rr',0.0,1.,k)
           endif ! micro_correlations(incl,iqr,k) /= undef_corr
           if ( micro_correlations(idx_Nc,idx_Nr,k) /= undef_corr ) then
              call hbuf_put_level('corr_Nc_Nr',micro_correlations(idx_Nc,idx_Nr,k),1.,k)
              corravg_count(idx_cor_Nc_nr,k) = corravg_count(idx_cor_Nc_nr,k) + 1
           else ! micro_correlations(incl,inr,k) == undef_corr
              call hbuf_put_level('corr_Nc_Nr',0.0,1.,k)
           endif ! micro_correlations(incl,inr,k) /= undef_corr
        end if
     
        if (doicemicro) then
           if ( micro_correlations(idx_Nc,idx_ri,k) /= undef_corr ) then
              call hbuf_put_level('corr_Nc_ri',micro_correlations(idx_Nc,idx_ri,k),1.,k) 
              corravg_count(idx_cor_Nc_ri,k) = corravg_count(idx_cor_Nc_ri,k) + 1
           else ! micro_correlations(incl,iqci,k) == undef_corr
              call hbuf_put_level('corr_Nc_ri',0.0,1.,k)
           endif ! micro_correlations(incl,iqci,k) /= undef_corr
           if ( micro_correlations(idx_Nc,idx_Ni,k) /= undef_corr ) then
              call hbuf_put_level('corr_Nc_Ni',micro_correlations(idx_Nc,idx_Ni,k),1.,k)
              corravg_count(idx_cor_Nc_Ni,k) = corravg_count(idx_cor_Nc_Ni,k) + 1
           else ! micro_correlations(incl,inci,k) == undef_corr
              call hbuf_put_level('corr_Nc_Ni',0.0,1.,k)
           endif ! micro_correlations(incl,inci,k) /= undef_corr
           if ( micro_correlations(idx_Nc,idx_rs,k) /= undef_corr ) then
              call hbuf_put_level('corr_Nc_rs',micro_correlations(idx_Nc,idx_rs,k),1.,k)
              corravg_count(idx_cor_Nc_rs,k) = corravg_count(idx_cor_Nc_rs,k) + 1
           else ! micro_correlations(incl,iqs,k) == undef_corr
              call hbuf_put_level('corr_Nc_rs',0.0,1.,k)
           endif ! micro_correlations(incl,iqs,k) /= undef_corr
           if ( micro_correlations(idx_Nc,idx_Ns,k) /= undef_corr ) then
              call hbuf_put_level('corr_Nc_Ns',micro_correlations(idx_Nc,idx_Ns,k),1.,k)
              corravg_count(idx_cor_Nc_Ns,k) = corravg_count(idx_cor_Nc_Ns,k) + 1
           else ! micro_correlations(incl,ins,k) == undef_corr
              call hbuf_put_level('corr_Nc_Ns',0.0,1.,k)
           endif ! micro_correlations(incl,ins,k) /= undef_corr
        end if    

        if (dograupel) then
           if ( micro_correlations(idx_Nc,idx_rg,k) /= undef_corr ) then
              call hbuf_put_level('corr_Nc_rg',micro_correlations(idx_Nc,idx_rg,k),1.,k)
              corravg_count(idx_cor_Nc_rg,k) = corravg_count(idx_cor_Nc_rg,k) + 1
           else ! micro_correlations(incl,iqg,k) == undef_corr
              call hbuf_put_level('corr_Nc_rg',0.0,1.,k)
           endif ! micro_correlations(incl,iqg,k) /= undef_corr
           if ( micro_correlations(idx_Nc,idx_Ng,k) /= undef_corr ) then
              call hbuf_put_level('corr_Nc_Ng',micro_correlations(idx_Nc,idx_Ng,k),1.,k)
              corravg_count(idx_cor_Nc_Ng,k) = corravg_count(idx_cor_Nc_Ng,k) + 1
           else ! micro_correlations(incl,ing,k) == undef_corr
              call hbuf_put_level('corr_Nc_Ng',0.0,1.,k)
           endif ! micro_correlations(incl,ing,k) /= undef_corr
        end if   

      end if !dopredictNc
  
!-------------
! Rainwater mixing ratio
!-------------
      if (doprecip) then 
      
         if ( micro_correlations(idx_rr,idx_Nr,k) /= undef_corr ) then
             call hbuf_put_level('corr_rr_Nr',micro_correlations(idx_rr,idx_Nr,k),1.,k)
             corravg_count(idx_cor_rr_nr,k) = corravg_count(idx_cor_rr_nr,k) + 1
          else ! micro_correlations(iqr,inr,k) == undef_corr
             call hbuf_put_level('corr_rr_Nr',0.0,1.,k)
          endif ! micro_correlations(iqr,inr,k) /= undef_corr

        if (doicemicro) then
           if ( micro_correlations(idx_rr,idx_ri,k) /= undef_corr ) then
              call hbuf_put_level('corr_rr_ri',micro_correlations(idx_rr,idx_ri,k),1.,k)
              corravg_count(idx_cor_rr_ri,k) = corravg_count(idx_cor_rr_ri,k) + 1
           else ! micro_correlations(iqr,iqci,k) == undef_corr
              call hbuf_put_level('corr_rr_ri',0.0,1.,k)
           endif ! micro_correlations(iqr,iqci,k) /= undef_corr
           if ( micro_correlations(idx_rr,idx_Ni,k) /= undef_corr ) then
              call hbuf_put_level('corr_rr_Ni',micro_correlations(idx_rr,idx_Ni,k),1.,k)
              corravg_count(idx_cor_rr_Ni,k) = corravg_count(idx_cor_rr_Ni,k) + 1
           else ! micro_correlations(iqr,inci,k) == undef_corr
              call hbuf_put_level('corr_rr_Ni',0.0,1.,k)
           endif ! micro_correlations(iqr,inci,k) /= undef_corr
           if ( micro_correlations(idx_rr,idx_rs,k) /= undef_corr ) then
              call hbuf_put_level('corr_rr_rs',micro_correlations(idx_rr,idx_rs,k),1.,k)
              corravg_count(idx_cor_rr_rs,k) = corravg_count(idx_cor_rr_rs,k) + 1
           else ! micro_correlations(iqr,iqs,k) == undef_corr
              call hbuf_put_level('corr_rr_rs',0.0,1.,k)
           endif ! micro_correlations(iqr,iqs,k) /= undef_corr
           if ( micro_correlations(idx_rr,idx_Ns,k) /= undef_corr ) then
              call hbuf_put_level('corr_rr_Ns',micro_correlations(idx_rr,idx_Ns,k),1.,k)
              corravg_count(idx_cor_rr_Ns,k) = corravg_count(idx_cor_rr_Ns,k) + 1
           else ! micro_correlations(iqr,ins,k) == undef_corr
              call hbuf_put_level('corr_rr_Ns',0.0,1.,k)
           endif ! micro_correlations(iqr,ins,k) /= undef_corr
        end if    

        if (dograupel) then
           if ( micro_correlations(idx_rr,idx_rg,k) /= undef_corr ) then
              call hbuf_put_level('corr_rr_rg',micro_correlations(idx_rr,idx_rg,k),1.,k)
              corravg_count(idx_cor_rr_rg,k) = corravg_count(idx_cor_rr_rg,k) + 1
           else ! micro_correlations(iqr,iqg,k) == undef_corr
              call hbuf_put_level('corr_rr_rg',0.0,1.,k)
           endif ! micro_correlations(iqr,iqg,k) /= undef_corr
           if ( micro_correlations(idx_rr,idx_Ng,k) /= undef_corr ) then
              call hbuf_put_level('corr_rr_Ng',micro_correlations(idx_rr,idx_Ng,k),1.,k)
              corravg_count(idx_cor_rr_Ng,k) = corravg_count(idx_cor_rr_Ng,k) + 1
           else ! micro_correlations(iqr,ing,k) == undef_corr
              call hbuf_put_level('corr_rr_Ng',0.0,1.,k)
           endif ! micro_correlations(iqr,ing,k) /= undef_corr
        end if
!-------------
! Rainwater number concentration
!-------------

        if (doicemicro) then
           if ( micro_correlations(idx_Nr,idx_ri,k) /= undef_corr ) then
              call hbuf_put_level('corr_Nr_ri',micro_correlations(idx_Nr,idx_ri,k),1.,k)
              corravg_count(idx_cor_Nr_ri,k) = corravg_count(idx_cor_Nr_ri,k) + 1
           else ! micro_correlations(inr,iqci,k) == undef_corr
              call hbuf_put_level('corr_Nr_ri',0.0,1.,k)
           endif ! micro_correlations(inr,iqci,k) /= undef_corr
           if ( micro_correlations(idx_Nr,idx_Ni,k) /= undef_corr ) then
              call hbuf_put_level('corr_Nr_Ni',micro_correlations(idx_Nr,idx_Ni,k),1.,k)
              corravg_count(idx_cor_Nr_Ni,k) = corravg_count(idx_cor_Nr_Ni,k) + 1
           else ! micro_correlations(inr,inci,k) == undef_corr
              call hbuf_put_level('corr_Nr_Ni',0.0,1.,k)
           endif ! micro_correlations(inr,inci,k) /= undef_corr
           if ( micro_correlations(idx_Nr,idx_rs,k) /= undef_corr ) then
              call hbuf_put_level('corr_Nr_rs',micro_correlations(idx_Nr,idx_rs,k),1.,k)
              corravg_count(idx_cor_Nr_rs,k) = corravg_count(idx_cor_Nr_rs,k) + 1
           else ! micro_correlations(inr,iqs,k) == undef_corr
              call hbuf_put_level('corr_Nr_rs',0.0,1.,k)
           endif ! micro_correlations(inr,iqs,k) /= undef_corr
           if ( micro_correlations(idx_Nr,idx_Ns,k) /= undef_corr ) then
              call hbuf_put_level('corr_Nr_Ns',micro_correlations(idx_Nr,idx_Ns,k),1.,k)
              corravg_count(idx_cor_Nr_Ns,k) = corravg_count(idx_cor_Nr_Ns,k) + 1
           else ! micro_correlations(inr,ins,k) == undef_corr
              call hbuf_put_level('corr_Nr_Ns',0.0,1.,k)
           endif ! micro_correlations(inr,ins,k) /= undef_corr
        end if

        if (dograupel) then
        if ( micro_correlations(idx_Nr,idx_rg,k) /= undef_corr ) then
              call hbuf_put_level('corr_Nr_rg',micro_correlations(idx_Nr,idx_rg,k),1.,k)
              corravg_count(idx_cor_Nr_rg,k) = corravg_count(idx_cor_Nr_rg,k) + 1
           else ! micro_correlations(inr,iqg,k) == undef_corr
              call hbuf_put_level('corr_Nr_rg',0.0,1.,k)
           endif ! micro_correlations(inr,iqg,k) /= undef_corr
           if ( micro_correlations(idx_Nr,idx_Ng,k) /= undef_corr ) then
              call hbuf_put_level('corr_Nr_Ng',micro_correlations(idx_Nr,idx_Ng,k),1.,k)
              corravg_count(idx_cor_Nr_Ng,k) = corravg_count(idx_cor_Nr_Ng,k) + 1
           else ! micro_correlations(inr,ing,k) == undef_corr
              call hbuf_put_level('corr_Nr_Ng',0.0,1.,k)
           endif ! micro_correlations(inr,ing,k) /= undef_corr
        end if

      end if !doprecip
!-------------
! Cloud-ice mixing ratio
!-------------

      if (doicemicro) then   
         if ( micro_correlations(idx_ri,idx_Ni,k) /= undef_corr ) then
            call hbuf_put_level('corr_ri_Ni',micro_correlations(idx_ri,idx_Ni,k),1.,k)
            corravg_count(idx_cor_ri_Ni,k) = corravg_count(idx_cor_ri_Ni,k) + 1
         else ! micro_correlations(iqci,inci,k) == undef_corr
            call hbuf_put_level('corr_ri_Ni',0.0,1.,k)
         endif ! micro_correlations(iqci,inci,k) /= undef_corr
         if ( micro_correlations(idx_ri,idx_rs,k) /= undef_corr ) then
            call hbuf_put_level('corr_ri_rs',micro_correlations(idx_ri,idx_rs,k),1.,k)
            corravg_count(idx_cor_ri_rs,k) = corravg_count(idx_cor_ri_rs,k) + 1
         else ! micro_correlations(iqci,iqs,k) == undef_corr
            call hbuf_put_level('corr_ri_rs',0.0,1.,k)
         endif ! micro_correlations(iqci,iqs,k) /= undef_corr
         if ( micro_correlations(idx_ri,idx_Ns,k) /= undef_corr ) then
            call hbuf_put_level('corr_ri_Ns',micro_correlations(idx_ri,idx_Ns,k),1.,k)
            corravg_count(idx_cor_ri_Ns,k) = corravg_count(idx_cor_ri_Ns,k) + 1
         else ! micro_correlations(iqci,ins,k) == undef_corr
            call hbuf_put_level('corr_ri_Ns',0.0,1.,k)
         endif ! micro_correlations(iqci,ins,k) /= undef_corr
      end if   

      if (dograupel) then
         if ( micro_correlations(idx_ri,idx_rg,k) /= undef_corr ) then
            call hbuf_put_level('corr_ri_rg',micro_correlations(idx_ri,idx_rg,k),1.,k)
            corravg_count(idx_cor_ri_rg,k) = corravg_count(idx_cor_ri_rg,k) + 1
         else ! micro_correlations(iqci,iqg,k) == undef_corr
            call hbuf_put_level('corr_ri_rg',0.0,1.,k)
         endif ! micro_correlations(iqci,iqg,k) /= undef_corr
         if ( micro_correlations(idx_ri,idx_Ng,k) /= undef_corr ) then
            call hbuf_put_level('corr_ri_Ng',micro_correlations(idx_ri,idx_Ng,k),1.,k)
            corravg_count(idx_cor_ri_Ng,k) = corravg_count(idx_cor_ri_Ng,k) + 1
         else ! micro_correlations(iqci,ing,k) == undef_corr
            call hbuf_put_level('corr_ri_Ng',0.0,1.,k)
         endif ! micro_correlations(iqci,ing,k) /= undef_corr
      end if
 
!-------------
! Cloud-ice number concentration
!-------------
      if (doicemicro) then   
         if ( micro_correlations(idx_Ni,idx_rs,k) /= undef_corr ) then
            call hbuf_put_level('corr_Ni_rs',micro_correlations(idx_Ni,idx_rs,k),1.,k)
            corravg_count(idx_cor_Ni_rs,k) = corravg_count(idx_cor_Ni_rs,k) + 1
         else ! micro_correlations(inci,iqs,k) == undef_corr
            call hbuf_put_level('corr_Ni_rs',0.0,1.,k)
         endif ! micro_correlations(inci,iqs,k) /= undef_corr
         if ( micro_correlations(idx_Ni,idx_Ns,k) /= undef_corr ) then
            call hbuf_put_level('corr_Ni_Ns',micro_correlations(idx_Ni,idx_Ns,k),1.,k)
            corravg_count(idx_cor_Ni_Ns,k) = corravg_count(idx_cor_Ni_Ns,k) + 1
         else ! micro_correlations(inci,ins,k) == undef_corr
            call hbuf_put_level('corr_Ni_Ns',0.0,1.,k)
         endif ! micro_correlations(inci,ins,k) /= undef_corr
      end if

      if (dograupel) then     
         if ( micro_correlations(idx_Ni,idx_rg,k) /= undef_corr ) then
            call hbuf_put_level('corr_Ni_rg',micro_correlations(idx_Ni,idx_rg,k),1.,k)
            corravg_count(idx_cor_Ni_rg,k) = corravg_count(idx_cor_Ni_rg,k) + 1
         else ! micro_correlations(inci,iqg,k) == undef_corr
            call hbuf_put_level('corr_Ni_rg',0.0,1.,k)
         endif ! micro_correlations(inci,iqg,k) /= undef_corr
         if ( micro_correlations(idx_Ni,idx_Ng,k) /= undef_corr ) then
            call hbuf_put_level('corr_Ni_Ng',micro_correlations(idx_Ni,idx_Ng,k),1.,k)
            corravg_count(idx_cor_Ni_Ng,k) = corravg_count(idx_cor_Ni_Ng,k) + 1
         else ! micro_correlations(inci,ing,k) == undef_corr
            call hbuf_put_level('corr_Ni_Ng',0.0,1.,k)
         endif ! micro_correlations(inci,ing,k) /= undef_corr
      end if 
!-------------
! Snow mixing ratio
!-------------

      if (doicemicro) then
         if ( micro_correlations(idx_rs,idx_Ns,k) /= undef_corr ) then
            call hbuf_put_level('corr_rs_Ns',micro_correlations(idx_rs,idx_Ns,k),1.,k)
            corravg_count(idx_cor_rs_Ns,k) = corravg_count(idx_cor_rs_Ns,k) + 1
         else ! micro_correlations(iqs,ins,k) == undef_corr
            call hbuf_put_level('corr_rs_Ns',0.0,1.,k)
         endif ! micro_correlations(iqs,ins,k) /= undef_corr
       
          if (dograupel) then   
             if ( micro_correlations(idx_rs,idx_rg,k) /= undef_corr ) then
                call hbuf_put_level('corr_rs_rg',micro_correlations(idx_rs,idx_rg,k),1.,k)
                corravg_count(idx_cor_rs_rg,k) = corravg_count(idx_cor_rs_rg,k) + 1
             else ! micro_correlations(iqs,iqg,k) == undef_corr
                call hbuf_put_level('corr_rs_rg',0.0,1.,k)
             endif ! micro_correlations(iqs,iqg,k) /= undef_corr
             if ( micro_correlations(idx_rs,idx_Ng,k) /= undef_corr ) then
                call hbuf_put_level('corr_rs_Ng',micro_correlations(idx_rs,idx_Ng,k),1.,k)
                corravg_count(idx_cor_rs_Ng,k) = corravg_count(idx_cor_rs_Ng,k) + 1
             else ! micro_correlations(iqs,ing,k) == undef_corr
                call hbuf_put_level('corr_rs_Ng',0.0,1.,k)
             endif ! micro_correlations(iqs,ing,k) /= undef_corr
!-------------
! Snow number concentration
!-------------
             if ( micro_correlations(idx_Ns,idx_rg,k) /= undef_corr ) then
                call hbuf_put_level('corr_Ns_rg',micro_correlations(idx_Ns,idx_rg,k),1.,k)
                corravg_count(idx_cor_Ns_rg,k) = corravg_count(idx_cor_Ns_rg,k) + 1
             else ! micro_correlations(ins,iqg,k) == undef_corr
                call hbuf_put_level('corr_Ns_rg',0.0,1.,k)
             endif ! micro_correlations(ins,iqg,k) /= undef_corr
             if ( micro_correlations(idx_Ns,idx_Ng,k) /= undef_corr ) then
                call hbuf_put_level('corr_Ns_Ng',micro_correlations(idx_Ns,idx_Ng,k),1.,k)
                corravg_count(idx_cor_Ns_Ng,k) = corravg_count(idx_cor_Ns_Ng,k) + 1
             else ! micro_correlations(ins,ing,k) == undef_corr
                call hbuf_put_level('corr_Ns_Ng',0.0,1.,k)
             endif ! micro_correlations(ins,ing,k) /= undef_corr
          end if
      end if

!-------------
! Graupel mixing ratio / number concentration
!-------------
      if (dograupel) then      
        if ( micro_correlations(idx_rg,idx_Ng,k) /= undef_corr ) then
            call hbuf_put_level('corr_rg_Ng',micro_correlations(idx_rg,idx_Ng,k),1.,k)
            corravg_count(idx_cor_rg_Ng,k) = corravg_count(idx_cor_rg_Ng,k) + 1
         else ! micro_correlations(iqg,ing,k) == undef_corr
            call hbuf_put_level('corr_rg_Ng',0.0,1.,k)
         endif ! micro_correlations(iqg,ing,k) /= undef_corr
      endif

   enddo ! k = 1, nzm, 1

!------------
! In-cloud correlations
!------------
   do k = 1, nzm, 1
!-------------
! chi
!-------------
      if ( ic_micro_correlations(idx_s,idx_w,k) /= undef_corr ) then
         call hbuf_put_level('ic_corr_chi_w',ic_micro_correlations(idx_s,idx_w,k),1.,k)
         corravg_count(idx_ic_cor_chi_w,k) = corravg_count(idx_ic_cor_chi_w,k) + 1
      else ! ic_micro_correlations(1,incl,k) == undef_corr
         call hbuf_put_level('ic_corr_chi_w',0.0,1.,k)
      endif ! ic_micro_correlations(1,incl,k) /= undef_corr
      
      if(dopredictNc) then
        if ( ic_micro_correlations(idx_s,idx_Nc,k) /= undef_corr ) then
           call hbuf_put_level('ic_corr_chi_Nc',ic_micro_correlations(idx_s,idx_Nc,k),1.,k)
           corravg_count(idx_ic_cor_chi_Nc,k) = corravg_count(idx_ic_cor_chi_Nc,k) + 1
        else ! ic_micro_correlations(1,incl,k) == undef_corr
           call hbuf_put_level('ic_corr_chi_Nc',0.0,1.,k)
        endif ! ic_micro_correlations(1,incl,k) /= undef_corr
      end if

      if (doprecip) then
      if ( ic_micro_correlations(idx_s,idx_rr,k) /= undef_corr ) then
           call hbuf_put_level('ic_corr_chi_rr',ic_micro_correlations(idx_s,idx_rr,k),1.,k)
           corravg_count(idx_ic_cor_chi_rr,k) = corravg_count(idx_ic_cor_chi_rr,k) + 1
        else ! ic_micro_correlations(1,iqr,k) == undef_corr
           call hbuf_put_level('ic_corr_chi_rr',0.0,1.,k)
        endif ! ic_micro_correlations(1,iqr,k) /= undef_corr
        if ( ic_micro_correlations(idx_s,idx_Nr,k) /= undef_corr ) then
           call hbuf_put_level('ic_corr_chi_Nr',ic_micro_correlations(idx_s,idx_Nr,k),1.,k)
           corravg_count(idx_ic_cor_chi_Nr,k) = corravg_count(idx_ic_cor_chi_Nr,k) + 1
        else ! ic_micro_correlations(1,inr,k) == undef_corr
           call hbuf_put_level('ic_corr_chi_Nr',0.0,1.,k)
        endif ! ! ic_micro_correlations(1,inr,k) /= undef_corr
      end if ! doprecip

      if (doicemicro) then
      if ( ic_micro_correlations(idx_s,idx_ri,k) /= undef_corr ) then
           call hbuf_put_level('ic_corr_chi_ri',ic_micro_correlations(idx_s,idx_ri,k),1.,k)
           corravg_count(idx_ic_cor_chi_ri,k) = corravg_count(idx_ic_cor_chi_ri,k) + 1
        else ! ic_micro_correlations(1,iqci,k) == undef_corr
           call hbuf_put_level('ic_corr_chi_ri',0.0,1.,k)
        endif ! ic_micro_correlations(1,iqci,k) /= undef_corr
        if ( ic_micro_correlations(idx_s,idx_Ni,k) /= undef_corr ) then
           call hbuf_put_level('ic_corr_chi_Ni',ic_micro_correlations(idx_s,idx_Ni,k),1.,k)
           corravg_count(idx_ic_cor_chi_Ni,k) = corravg_count(idx_ic_cor_chi_Ni,k) + 1
        else ! ic_micro_correlations(1,inci,k) == undef_corr
           call hbuf_put_level('ic_corr_chi_Ni',0.0,1.,k)
        endif ! ic_micro_correlations(1,inci,k) /= undef_corr
        if ( ic_micro_correlations(idx_s,idx_rs,k) /= undef_corr ) then
           call hbuf_put_level('ic_corr_chi_rs',ic_micro_correlations(idx_s,idx_rs,k),1.,k)
           corravg_count(idx_ic_cor_chi_rs,k) = corravg_count(idx_ic_cor_chi_rs,k) + 1
        else ! ic_micro_correlations(1,iqs,k) == undef_corr
           call hbuf_put_level('ic_corr_chi_rs',0.0,1.,k)
        endif ! ic_micro_correlations(1,iqs,k) /= undef_corr
        if ( ic_micro_correlations(idx_s,idx_Ns,k) /= undef_corr ) then
           call hbuf_put_level('ic_corr_chi_Ns',ic_micro_correlations(idx_s,idx_Ns,k),1.,k)
           corravg_count(idx_ic_cor_chi_Ns,k) = corravg_count(idx_ic_cor_chi_Ns,k) + 1
        else ! ic_micro_correlations(1,ins,k) == undef_corr
           call hbuf_put_level('ic_corr_chi_Ns',0.0,1.,k)
        endif ! ic_micro_correlations(1,ins,k) /= undef_corr
      end if ! doicemicro

      if (dograupel) then
        if ( ic_micro_correlations(idx_s,idx_rg,k) /= undef_corr ) then
           call hbuf_put_level('ic_corr_chi_rg',ic_micro_correlations(idx_s,idx_rg,k),1.,k)
           corravg_count(idx_ic_cor_chi_rg,k) = corravg_count(idx_ic_cor_chi_rg,k) + 1
        else ! ic_micro_correlations(1,iqg,k) == undef_corr
           call hbuf_put_level('ic_corr_chi_rg',0.0,1.,k)
        endif ! ic_micro_correlations(1,iqg,k) /= undef_corr
        if ( ic_micro_correlations(idx_s,idx_Ng,k) /= undef_corr ) then
           call hbuf_put_level('ic_corr_chi_Ng',ic_micro_correlations(idx_s,idx_Ng,k),1.,k)
           corravg_count(idx_ic_cor_chi_Ng,k) = corravg_count(idx_ic_cor_chi_Ng,k) + 1
        else ! ic_micro_correlations(1,ing,k) == undef_corr
           call hbuf_put_level('ic_corr_chi_Ng',0.0,1.,k)
        endif ! ic_micro_correlations(1,ing,k) /= undef_corr
      end if ! dograupel

!-------------
! Vertical Velocity
!-------------
      if(dopredictNc) then
        if ( ic_micro_correlations(idx_w,idx_Nc,k) /= undef_corr ) then
           call hbuf_put_level('ic_corr_w_Nc',ic_micro_correlations(idx_w,idx_Nc,k),1.,k)
           corravg_count(idx_ic_cor_w_Nc,k) = corravg_count(idx_ic_cor_w_Nc,k) + 1
        else ! ic_micro_correlations(1,incl,k) == undef_corr
           call hbuf_put_level('ic_corr_w_Nc',0.0,1.,k)
        endif ! ic_micro_correlations(1,incl,k) /= undef_corr
      end if

      if (doprecip) then
      if ( ic_micro_correlations(idx_w,idx_rr,k) /= undef_corr ) then
           call hbuf_put_level('ic_corr_w_rr',ic_micro_correlations(idx_w,idx_rr,k),1.,k)
           corravg_count(idx_ic_cor_w_rr,k) = corravg_count(idx_ic_cor_w_rr,k) + 1
        else ! ic_micro_correlations(1,iqr,k) == undef_corr
           call hbuf_put_level('ic_corr_w_rr',0.0,1.,k)
        endif ! ic_micro_correlations(1,iqr,k) /= undef_corr
        if ( ic_micro_correlations(idx_w,idx_Nr,k) /= undef_corr ) then
           call hbuf_put_level('ic_corr_w_Nr',ic_micro_correlations(idx_w,idx_Nr,k),1.,k)
           corravg_count(idx_ic_cor_w_nr,k) = corravg_count(idx_ic_cor_w_nr,k) + 1
        else ! ic_micro_correlations(1,inr,k) == undef_corr
           call hbuf_put_level('ic_corr_w_Nr',0.0,1.,k)
        endif ! ! ic_micro_correlations(1,inr,k) /= undef_corr
      end if ! doprecip

      if (doicemicro) then
      if ( ic_micro_correlations(idx_w,idx_ri,k) /= undef_corr ) then
           call hbuf_put_level('ic_corr_w_ri',ic_micro_correlations(idx_w,idx_ri,k),1.,k)
           corravg_count(idx_ic_cor_w_ri,k) = corravg_count(idx_ic_cor_w_ri,k) + 1
        else ! ic_micro_correlations(1,iqci,k) == undef_corr
           call hbuf_put_level('ic_corr_w_ri',0.0,1.,k)
        endif ! ic_micro_correlations(1,iqci,k) /= undef_corr
        if ( ic_micro_correlations(idx_w,idx_Ni,k) /= undef_corr ) then
           call hbuf_put_level('ic_corr_w_Ni',ic_micro_correlations(idx_w,idx_Ni,k),1.,k)
           corravg_count(idx_ic_cor_w_Ni,k) = corravg_count(idx_ic_cor_w_Ni,k) + 1
        else ! ic_micro_correlations(1,inci,k) == undef_corr
           call hbuf_put_level('ic_corr_w_Ni',0.0,1.,k)
        endif ! ic_micro_correlations(1,inci,k) /= undef_corr
        if ( ic_micro_correlations(idx_w,idx_rs,k) /= undef_corr ) then
           call hbuf_put_level('ic_corr_w_rs',ic_micro_correlations(idx_w,idx_rs,k),1.,k)
           corravg_count(idx_ic_cor_w_rs,k) = corravg_count(idx_ic_cor_w_rs,k) + 1
        else ! ic_micro_correlations(1,iqs,k) == undef_corr
           call hbuf_put_level('ic_corr_w_rs',0.0,1.,k)
        endif ! ic_micro_correlations(1,iqs,k) /= undef_corr
        if ( ic_micro_correlations(idx_w,idx_Ns,k) /= undef_corr ) then
           call hbuf_put_level('ic_corr_w_Ns',ic_micro_correlations(idx_w,idx_Ns,k),1.,k)
           corravg_count(idx_ic_cor_w_Ns,k) = corravg_count(idx_ic_cor_w_Ns,k) + 1
        else ! ic_micro_correlations(1,ins,k) == undef_corr
           call hbuf_put_level('ic_corr_w_Ns',0.0,1.,k)
        endif ! ic_micro_correlations(1,ins,k) /= undef_corr
      end if ! doicemicro

      if (dograupel) then
        if ( ic_micro_correlations(idx_w,idx_rg,k) /= undef_corr ) then
           call hbuf_put_level('ic_corr_w_rg',ic_micro_correlations(idx_w,idx_rg,k),1.,k)
           corravg_count(idx_ic_cor_w_rg,k) = corravg_count(idx_ic_cor_w_rg,k) + 1
        else ! ic_micro_correlations(1,iqg,k) == undef_corr
           call hbuf_put_level('ic_corr_w_rg',0.0,1.,k)
        endif ! ic_micro_correlations(1,iqg,k) /= undef_corr
        if ( ic_micro_correlations(idx_w,idx_Ng,k) /= undef_corr ) then
           call hbuf_put_level('ic_corr_w_Ng',ic_micro_correlations(idx_w,idx_Ng,k),1.,k)
           corravg_count(idx_ic_cor_w_Ng,k) = corravg_count(idx_ic_cor_w_Ng,k) + 1
        else ! ic_micro_correlations(1,ing,k) == undef_corr
           call hbuf_put_level('ic_corr_w_Ng',0.0,1.,k)
        endif ! ic_micro_correlations(1,ing,k) /= undef_corr
      end if ! dograupel

!-------------
! Cloud droplet number concentration
!-------------
      if (dopredictNc) then
     
        if (doprecip) then
           if ( ic_micro_correlations(idx_Nc,idx_rr,k) /= undef_corr ) then
              call hbuf_put_level('ic_corr_Nc_rr',ic_micro_correlations(idx_Nc,idx_rr,k),1.,k)
              corravg_count(idx_ic_cor_Nc_rr,k) = corravg_count(idx_ic_cor_Nc_rr,k) + 1
           else ! ic_micro_correlations(incl,iqr,k) == undef_corr
              call hbuf_put_level('ic_corr_Nc_rr',0.0,1.,k)
           endif ! ic_micro_correlations(incl,iqr,k) /= undef_corr
           if ( ic_micro_correlations(idx_Nc,idx_Nr,k) /= undef_corr ) then
              call hbuf_put_level('ic_corr_Nc_Nr',ic_micro_correlations(idx_Nc,idx_Nr,k),1.,k)
              corravg_count(idx_ic_cor_Nc_nr,k) = corravg_count(idx_ic_cor_Nc_nr,k) + 1
           else ! ic_micro_correlations(incl,inr,k) == undef_corr
              call hbuf_put_level('ic_corr_Nc_Nr',0.0,1.,k)
           endif ! ic_micro_correlations(incl,inr,k) /= undef_corr
        end if
     
        if (doicemicro) then
           if ( ic_micro_correlations(idx_Nc,idx_ri,k) /= undef_corr ) then
              call hbuf_put_level('ic_corr_Nc_ri',ic_micro_correlations(idx_Nc,idx_ri,k),1.,k) 
              corravg_count(idx_ic_cor_Nc_ri,k) = corravg_count(idx_ic_cor_Nc_ri,k) + 1
           else ! ic_micro_correlations(incl,iqci,k) == undef_corr
              call hbuf_put_level('ic_corr_Nc_ri',0.0,1.,k)
           endif ! ic_micro_correlations(incl,iqci,k) /= undef_corr
           if ( ic_micro_correlations(idx_Nc,idx_Ni,k) /= undef_corr ) then
              call hbuf_put_level('ic_corr_Nc_Ni',ic_micro_correlations(idx_Nc,idx_Ni,k),1.,k)
              corravg_count(idx_ic_cor_Nc_Ni,k) = corravg_count(idx_ic_cor_Nc_Ni,k) + 1
           else ! ic_micro_correlations(incl,inci,k) == undef_corr
              call hbuf_put_level('ic_corr_Nc_Ni',0.0,1.,k)
           endif ! ic_micro_correlations(incl,inci,k) /= undef_corr
           if ( ic_micro_correlations(idx_Nc,idx_rs,k) /= undef_corr ) then
              call hbuf_put_level('ic_corr_Nc_rs',ic_micro_correlations(idx_Nc,idx_rs,k),1.,k)
              corravg_count(idx_ic_cor_Nc_rs,k) = corravg_count(idx_ic_cor_Nc_rs,k) + 1
           else ! ic_micro_correlations(incl,iqs,k) == undef_corr
              call hbuf_put_level('ic_corr_Nc_rs',0.0,1.,k)
           endif ! ic_micro_correlations(incl,iqs,k) /= undef_corr
           if ( ic_micro_correlations(idx_Nc,idx_Ns,k) /= undef_corr ) then
              call hbuf_put_level('ic_corr_Nc_Ns',ic_micro_correlations(idx_Nc,idx_Ns,k),1.,k)
              corravg_count(idx_ic_cor_Nc_Ns,k) = corravg_count(idx_ic_cor_Nc_Ns,k) + 1
           else ! ic_micro_correlations(incl,ins,k) == undef_corr
              call hbuf_put_level('ic_corr_Nc_Ns',0.0,1.,k)
           endif ! ic_micro_correlations(incl,ins,k) /= undef_corr
        end if    

        if (dograupel) then
           if ( ic_micro_correlations(idx_Nc,idx_rg,k) /= undef_corr ) then
              call hbuf_put_level('ic_corr_Nc_rg',ic_micro_correlations(idx_Nc,idx_rg,k),1.,k)
              corravg_count(idx_ic_cor_Nc_rg,k) = corravg_count(idx_ic_cor_Nc_rg,k) + 1
           else ! ic_micro_correlations(incl,iqg,k) == undef_corr
              call hbuf_put_level('ic_corr_Nc_rg',0.0,1.,k)
           endif ! ic_micro_correlations(incl,iqg,k) /= undef_corr
           if ( ic_micro_correlations(idx_Nc,idx_Ng,k) /= undef_corr ) then
              call hbuf_put_level('ic_corr_Nc_Ng',ic_micro_correlations(idx_Nc,idx_Ng,k),1.,k)
              corravg_count(idx_ic_cor_Nc_Ng,k) = corravg_count(idx_ic_cor_Nc_Ng,k) + 1
           else ! ic_micro_correlations(incl,ing,k) == undef_corr
              call hbuf_put_level('ic_corr_Nc_Ng',0.0,1.,k)
           endif ! ic_micro_correlations(incl,ing,k) /= undef_corr
        end if   

      end if !dopredictNc
  
!-------------
! Rainwater mixing ratio
!-------------
      if (doprecip) then 
      
         if ( ic_micro_correlations(idx_rr,idx_Nr,k) /= undef_corr ) then
             call hbuf_put_level('ic_corr_rr_Nr',ic_micro_correlations(idx_rr,idx_Nr,k),1.,k)
             corravg_count(idx_ic_cor_rr_nr,k) = corravg_count(idx_ic_cor_rr_nr,k) + 1
          else ! ic_micro_correlations(iqr,inr,k) == undef_corr
             call hbuf_put_level('ic_corr_rr_Nr',0.0,1.,k)
          endif ! ic_micro_correlations(iqr,inr,k) /= undef_corr

        if (doicemicro) then
           if ( ic_micro_correlations(idx_rr,idx_ri,k) /= undef_corr ) then
              call hbuf_put_level('ic_corr_rr_ri',ic_micro_correlations(idx_rr,idx_ri,k),1.,k)
              corravg_count(idx_ic_cor_rr_ri,k) = corravg_count(idx_ic_cor_rr_ri,k) + 1
           else ! ic_micro_correlations(iqr,iqci,k) == undef_corr
              call hbuf_put_level('ic_corr_rr_ri',0.0,1.,k)
           endif ! ic_micro_correlations(iqr,iqci,k) /= undef_corr
           if ( ic_micro_correlations(idx_rr,idx_Ni,k) /= undef_corr ) then
              call hbuf_put_level('ic_corr_rr_Ni',ic_micro_correlations(idx_rr,idx_Ni,k),1.,k)
              corravg_count(idx_ic_cor_rr_Ni,k) = corravg_count(idx_ic_cor_rr_Ni,k) + 1
           else ! ic_micro_correlations(iqr,inci,k) == undef_corr
              call hbuf_put_level('ic_corr_rr_Ni',0.0,1.,k)
           endif ! ic_micro_correlations(iqr,inci,k) /= undef_corr
           if ( ic_micro_correlations(idx_rr,idx_rs,k) /= undef_corr ) then
              call hbuf_put_level('ic_corr_rr_rs',ic_micro_correlations(idx_rr,idx_rs,k),1.,k)
              corravg_count(idx_ic_cor_rr_rs,k) = corravg_count(idx_ic_cor_rr_rs,k) + 1
           else ! ic_micro_correlations(iqr,iqs,k) == undef_corr
              call hbuf_put_level('ic_corr_rr_rs',0.0,1.,k)
           endif ! ic_micro_correlations(iqr,iqs,k) /= undef_corr
           if ( ic_micro_correlations(idx_rr,idx_Ns,k) /= undef_corr ) then
              call hbuf_put_level('ic_corr_rr_Ns',ic_micro_correlations(idx_rr,idx_Ns,k),1.,k)
              corravg_count(idx_ic_cor_rr_Ns,k) = corravg_count(idx_ic_cor_rr_Ns,k) + 1
           else ! ic_micro_correlations(iqr,ins,k) == undef_corr
              call hbuf_put_level('ic_corr_rr_Ns',0.0,1.,k)
           endif ! ic_micro_correlations(iqr,ins,k) /= undef_corr
        end if    

        if (dograupel) then
           if ( ic_micro_correlations(idx_rr,idx_rg,k) /= undef_corr ) then
              call hbuf_put_level('ic_corr_rr_rg',ic_micro_correlations(idx_rr,idx_rg,k),1.,k)
              corravg_count(idx_ic_cor_rr_rg,k) = corravg_count(idx_ic_cor_rr_rg,k) + 1
           else ! ic_micro_correlations(iqr,iqg,k) == undef_corr
              call hbuf_put_level('ic_corr_rr_rg',0.0,1.,k)
           endif ! ic_micro_correlations(iqr,iqg,k) /= undef_corr
           if ( ic_micro_correlations(idx_rr,idx_Ng,k) /= undef_corr ) then
              call hbuf_put_level('ic_corr_rr_Ng',ic_micro_correlations(idx_rr,idx_Ng,k),1.,k)
              corravg_count(idx_ic_cor_rr_Ng,k) = corravg_count(idx_ic_cor_rr_Ng,k) + 1
           else ! ic_micro_correlations(iqr,ing,k) == undef_corr
              call hbuf_put_level('ic_corr_rr_Ng',0.0,1.,k)
           endif ! ic_micro_correlations(iqr,ing,k) /= undef_corr
        end if
!-------------
! Rainwater number concentration
!-------------

        if (doicemicro) then
           if ( ic_micro_correlations(idx_Nr,idx_ri,k) /= undef_corr ) then
              call hbuf_put_level('ic_corr_Nr_ri',ic_micro_correlations(idx_Nr,idx_ri,k),1.,k)
              corravg_count(idx_ic_cor_Nr_ri,k) = corravg_count(idx_ic_cor_Nr_ri,k) + 1
           else ! ic_micro_correlations(inr,iqci,k) == undef_corr
              call hbuf_put_level('ic_corr_Nr_ri',0.0,1.,k)
           endif ! ic_micro_correlations(inr,iqci,k) /= undef_corr
           if ( ic_micro_correlations(idx_Nr,idx_Ni,k) /= undef_corr ) then
              call hbuf_put_level('ic_corr_Nr_Ni',ic_micro_correlations(idx_Nr,idx_Ni,k),1.,k)
              corravg_count(idx_ic_cor_Nr_Ni,k) = corravg_count(idx_ic_cor_Nr_Ni,k) + 1
           else ! ic_micro_correlations(inr,inci,k) == undef_corr
              call hbuf_put_level('ic_corr_Nr_Ni',0.0,1.,k)
           endif ! ic_micro_correlations(inr,inci,k) /= undef_corr
           if ( ic_micro_correlations(idx_Nr,idx_rs,k) /= undef_corr ) then
              call hbuf_put_level('ic_corr_Nr_rs',ic_micro_correlations(idx_Nr,idx_rs,k),1.,k)
              corravg_count(idx_ic_cor_Nr_rs,k) = corravg_count(idx_ic_cor_Nr_rs,k) + 1
           else ! ic_micro_correlations(inr,iqs,k) == undef_corr
              call hbuf_put_level('ic_corr_Nr_rs',0.0,1.,k)
           endif ! ic_micro_correlations(inr,iqs,k) /= undef_corr
           if ( ic_micro_correlations(idx_Nr,idx_Ns,k) /= undef_corr ) then
              call hbuf_put_level('ic_corr_Nr_Ns',ic_micro_correlations(idx_Nr,idx_Ns,k),1.,k)
              corravg_count(idx_ic_cor_Nr_Ns,k) = corravg_count(idx_ic_cor_Nr_Ns,k) + 1
           else ! ic_micro_correlations(inr,ins,k) == undef_corr
              call hbuf_put_level('ic_corr_Nr_Ns',0.0,1.,k)
           endif ! ic_micro_correlations(inr,ins,k) /= undef_corr
        end if

        if (dograupel) then
        if ( ic_micro_correlations(idx_Nr,idx_rg,k) /= undef_corr ) then
              call hbuf_put_level('ic_corr_Nr_rg',ic_micro_correlations(idx_Nr,idx_rg,k),1.,k)
              corravg_count(idx_ic_cor_Nr_rg,k) = corravg_count(idx_ic_cor_Nr_rg,k) + 1
           else ! ic_micro_correlations(inr,iqg,k) == undef_corr
              call hbuf_put_level('ic_corr_Nr_rg',0.0,1.,k)
           endif ! ic_micro_correlations(inr,iqg,k) /= undef_corr
           if ( ic_micro_correlations(idx_Nr,idx_Ng,k) /= undef_corr ) then
              call hbuf_put_level('ic_corr_Nr_Ng',ic_micro_correlations(idx_Nr,idx_Ng,k),1.,k)
              corravg_count(idx_ic_cor_Nr_Ng,k) = corravg_count(idx_ic_cor_Nr_Ng,k) + 1
           else ! ic_micro_correlations(inr,ing,k) == undef_corr
              call hbuf_put_level('ic_corr_Nr_Ng',0.0,1.,k)
           endif ! ic_micro_correlations(inr,ing,k) /= undef_corr
        end if

      end if !doprecip
!-------------
! Cloud-ice mixing ratio
!-------------

      if (doicemicro) then   
         if ( ic_micro_correlations(idx_ri,idx_Ni,k) /= undef_corr ) then
            call hbuf_put_level('ic_corr_ri_Ni',ic_micro_correlations(idx_ri,idx_Ni,k),1.,k)
            corravg_count(idx_ic_cor_ri_Ni,k) = corravg_count(idx_ic_cor_ri_Ni,k) + 1
         else ! ic_micro_correlations(iqci,inci,k) == undef_corr
            call hbuf_put_level('ic_corr_ri_Ni',0.0,1.,k)
         endif ! ic_micro_correlations(iqci,inci,k) /= undef_corr
         if ( ic_micro_correlations(idx_ri,idx_rs,k) /= undef_corr ) then
            call hbuf_put_level('ic_corr_ri_rs',ic_micro_correlations(idx_ri,idx_rs,k),1.,k)
            corravg_count(idx_ic_cor_ri_rs,k) = corravg_count(idx_ic_cor_ri_rs,k) + 1
         else ! ic_micro_correlations(iqci,iqs,k) == undef_corr
            call hbuf_put_level('ic_corr_ri_rs',0.0,1.,k)
         endif ! ic_micro_correlations(iqci,iqs,k) /= undef_corr
         if ( ic_micro_correlations(idx_ri,idx_Ns,k) /= undef_corr ) then
            call hbuf_put_level('ic_corr_ri_Ns',ic_micro_correlations(idx_ri,idx_Ns,k),1.,k)
            corravg_count(idx_ic_cor_ri_Ns,k) = corravg_count(idx_ic_cor_ri_Ns,k) + 1
         else ! ic_micro_correlations(iqci,ins,k) == undef_corr
            call hbuf_put_level('ic_corr_ri_Ns',0.0,1.,k)
         endif ! ic_micro_correlations(iqci,ins,k) /= undef_corr
      end if   

      if (dograupel) then
         if ( ic_micro_correlations(idx_ri,idx_rg,k) /= undef_corr ) then
            call hbuf_put_level('ic_corr_ri_rg',ic_micro_correlations(idx_ri,idx_rg,k),1.,k)
            corravg_count(idx_ic_cor_ri_rg,k) = corravg_count(idx_ic_cor_ri_rg,k) + 1
         else ! ic_micro_correlations(iqci,iqg,k) == undef_corr
            call hbuf_put_level('ic_corr_ri_rg',0.0,1.,k)
         endif ! ic_micro_correlations(iqci,iqg,k) /= undef_corr
         if ( ic_micro_correlations(idx_ri,idx_Ng,k) /= undef_corr ) then
            call hbuf_put_level('ic_corr_ri_Ng',ic_micro_correlations(idx_ri,idx_Ng,k),1.,k)
            corravg_count(idx_ic_cor_ri_Ng,k) = corravg_count(idx_ic_cor_ri_Ng,k) + 1
         else ! ic_micro_correlations(iqci,ing,k) == undef_corr
            call hbuf_put_level('ic_corr_ri_Ng',0.0,1.,k)
         endif ! ic_micro_correlations(iqci,ing,k) /= undef_corr
      end if
 
!-------------
! Cloud-ice number concentration
!-------------
      if (doicemicro) then   
         if ( ic_micro_correlations(idx_Ni,idx_rs,k) /= undef_corr ) then
            call hbuf_put_level('ic_corr_Ni_rs',ic_micro_correlations(idx_Ni,idx_rs,k),1.,k)
            corravg_count(idx_ic_cor_Ni_rs,k) = corravg_count(idx_ic_cor_Ni_rs,k) + 1
         else ! ic_micro_correlations(inci,iqs,k) == undef_corr
            call hbuf_put_level('ic_corr_Ni_rs',0.0,1.,k)
         endif ! ic_micro_correlations(inci,iqs,k) /= undef_corr
         if ( ic_micro_correlations(idx_Ni,idx_Ns,k) /= undef_corr ) then
            call hbuf_put_level('ic_corr_Ni_Ns',ic_micro_correlations(idx_Ni,idx_Ns,k),1.,k)
            corravg_count(idx_ic_cor_Ni_Ns,k) = corravg_count(idx_ic_cor_Ni_Ns,k) + 1
         else ! ic_micro_correlations(inci,ins,k) == undef_corr
            call hbuf_put_level('ic_corr_Ni_Ns',0.0,1.,k)
         endif ! ic_micro_correlations(inci,ins,k) /= undef_corr
      end if

      if (dograupel) then     
         if ( ic_micro_correlations(idx_Ni,idx_rg,k) /= undef_corr ) then
            call hbuf_put_level('ic_corr_Ni_rg',ic_micro_correlations(idx_Ni,idx_rg,k),1.,k)
            corravg_count(idx_ic_cor_Ni_rg,k) = corravg_count(idx_ic_cor_Ni_rg,k) + 1
         else ! ic_micro_correlations(inci,iqg,k) == undef_corr
            call hbuf_put_level('ic_corr_Ni_rg',0.0,1.,k)
         endif ! ic_micro_correlations(inci,iqg,k) /= undef_corr
         if ( ic_micro_correlations(idx_Ni,idx_Ng,k) /= undef_corr ) then
            call hbuf_put_level('ic_corr_Ni_Ng',ic_micro_correlations(idx_Ni,idx_Ng,k),1.,k)
            corravg_count(idx_ic_cor_Ni_Ng,k) = corravg_count(idx_ic_cor_Ni_Ng,k) + 1
         else ! ic_micro_correlations(inci,ing,k) == undef_corr
            call hbuf_put_level('ic_corr_Ni_Ng',0.0,1.,k)
         endif ! ic_micro_correlations(inci,ing,k) /= undef_corr
      end if 
!-------------
! Snow mixing ratio
!-------------

      if (doicemicro) then
         if ( ic_micro_correlations(idx_rs,idx_Ns,k) /= undef_corr ) then
            call hbuf_put_level('ic_corr_rs_Ns',ic_micro_correlations(idx_rs,idx_Ns,k),1.,k)
            corravg_count(idx_ic_cor_rs_Ns,k) = corravg_count(idx_ic_cor_rs_Ns,k) + 1
         else ! ic_micro_correlations(iqs,ins,k) == undef_corr
            call hbuf_put_level('ic_corr_rs_Ns',0.0,1.,k)
         endif ! ic_micro_correlations(iqs,ins,k) /= undef_corr
       
          if (dograupel) then   
             if ( ic_micro_correlations(idx_rs,idx_rg,k) /= undef_corr ) then
                call hbuf_put_level('ic_corr_rs_rg',ic_micro_correlations(idx_rs,idx_rg,k),1.,k)
                corravg_count(idx_ic_cor_rs_rg,k) = corravg_count(idx_ic_cor_rs_rg,k) + 1
             else ! ic_micro_correlations(iqs,iqg,k) == undef_corr
                call hbuf_put_level('ic_corr_rs_rg',0.0,1.,k)
             endif ! ic_micro_correlations(iqs,iqg,k) /= undef_corr
             if ( ic_micro_correlations(idx_rs,idx_Ng,k) /= undef_corr ) then
                call hbuf_put_level('ic_corr_rs_Ng',ic_micro_correlations(idx_rs,idx_Ng,k),1.,k)
                corravg_count(idx_ic_cor_rs_Ng,k) = corravg_count(idx_ic_cor_rs_Ng,k) + 1
             else ! ic_micro_correlations(iqs,ing,k) == undef_corr
                call hbuf_put_level('ic_corr_rs_Ng',0.0,1.,k)
             endif ! ic_micro_correlations(iqs,ing,k) /= undef_corr
!-------------
! Snow number concentration
!-------------
             if ( ic_micro_correlations(idx_Ns,idx_rg,k) /= undef_corr ) then
                call hbuf_put_level('ic_corr_Ns_rg',ic_micro_correlations(idx_Ns,idx_rg,k),1.,k)
                corravg_count(idx_ic_cor_Ns_rg,k) = corravg_count(idx_ic_cor_Ns_rg,k) + 1
             else ! ic_micro_correlations(ins,iqg,k) == undef_corr
                call hbuf_put_level('ic_corr_Ns_rg',0.0,1.,k)
             endif ! ic_micro_correlations(ins,iqg,k) /= undef_corr
             if ( ic_micro_correlations(idx_Ns,idx_Ng,k) /= undef_corr ) then
                call hbuf_put_level('ic_corr_Ns_Ng',ic_micro_correlations(idx_Ns,idx_Ng,k),1.,k)
                corravg_count(idx_ic_cor_Ns_Ng,k) = corravg_count(idx_ic_cor_Ns_Ng,k) + 1
             else ! ic_micro_correlations(ins,ing,k) == undef_corr
                call hbuf_put_level('ic_corr_Ns_Ng',0.0,1.,k)
             endif ! ic_micro_correlations(ins,ing,k) /= undef_corr
          end if
      end if

!-------------
! Graupel mixing ratio / number concentration
!-------------
      if (dograupel) then      
        if ( ic_micro_correlations(idx_rg,idx_Ng,k) /= undef_corr ) then
            call hbuf_put_level('ic_corr_rg_Ng',ic_micro_correlations(idx_rg,idx_Ng,k),1.,k)
            corravg_count(idx_ic_cor_rg_Ng,k) = corravg_count(idx_ic_cor_rg_Ng,k) + 1
         else ! ic_micro_correlations(iqg,ing,k) == undef_corr
            call hbuf_put_level('ic_corr_rg_Ng',0.0,1.,k)
         endif ! ic_micro_correlations(iqg,ing,k) /= undef_corr
      endif

    end do ! do k = 1, nzm, 1

!------------
! Out of cloud correlations
!------------
   do k = 1, nzm, 1
!-------------
! chi
!-------------
      if ( oc_micro_correlations(idx_s,idx_w,k) /= undef_corr ) then
         call hbuf_put_level('oc_corr_chi_w',oc_micro_correlations(idx_s,idx_w,k),1.,k)
         corravg_count(idx_oc_cor_chi_w,k) = corravg_count(idx_oc_cor_chi_w,k) + 1
      else ! oc_micro_correlations(1,incl,k) == undef_corr
         call hbuf_put_level('oc_corr_chi_w',0.0,1.,k)
      endif ! oc_micro_correlations(1,incl,k) /= undef_corr
      
      if(dopredictNc) then
        if ( oc_micro_correlations(idx_s,idx_Nc,k) /= undef_corr ) then
           call hbuf_put_level('oc_corr_chi_Nc',oc_micro_correlations(idx_s,idx_Nc,k),1.,k)
           corravg_count(idx_oc_cor_chi_Nc,k) = corravg_count(idx_oc_cor_chi_Nc,k) + 1
        else ! oc_micro_correlations(1,incl,k) == undef_corr
           call hbuf_put_level('oc_corr_chi_Nc',0.0,1.,k)
        endif ! oc_micro_correlations(1,incl,k) /= undef_corr
      end if

      if (doprecip) then
      if ( oc_micro_correlations(idx_s,idx_rr,k) /= undef_corr ) then
           call hbuf_put_level('oc_corr_chi_rr',oc_micro_correlations(idx_s,idx_rr,k),1.,k)
           corravg_count(idx_oc_cor_chi_rr,k) = corravg_count(idx_oc_cor_chi_rr,k) + 1
        else ! oc_micro_correlations(1,iqr,k) == undef_corr
           call hbuf_put_level('oc_corr_chi_rr',0.0,1.,k)
        endif ! oc_micro_correlations(1,iqr,k) /= undef_corr
        if ( oc_micro_correlations(idx_s,idx_Nr,k) /= undef_corr ) then
           call hbuf_put_level('oc_corr_chi_Nr',oc_micro_correlations(idx_s,idx_Nr,k),1.,k)
           corravg_count(idx_oc_cor_chi_Nr,k) = corravg_count(idx_oc_cor_chi_Nr,k) + 1
        else ! oc_micro_correlations(1,inr,k) == undef_corr
           call hbuf_put_level('oc_corr_chi_Nr',0.0,1.,k)
        endif ! ! oc_micro_correlations(1,inr,k) /= undef_corr
      end if ! doprecip

      if (doicemicro) then
      if ( oc_micro_correlations(idx_s,idx_ri,k) /= undef_corr ) then
           call hbuf_put_level('oc_corr_chi_ri',oc_micro_correlations(idx_s,idx_ri,k),1.,k)
           corravg_count(idx_oc_cor_chi_ri,k) = corravg_count(idx_oc_cor_chi_ri,k) + 1
        else ! oc_micro_correlations(1,iqci,k) == undef_corr
           call hbuf_put_level('oc_corr_chi_ri',0.0,1.,k)
        endif ! oc_micro_correlations(1,iqci,k) /= undef_corr
        if ( oc_micro_correlations(idx_s,idx_Ni,k) /= undef_corr ) then
           call hbuf_put_level('oc_corr_chi_Ni',oc_micro_correlations(idx_s,idx_Ni,k),1.,k)
           corravg_count(idx_oc_cor_chi_Ni,k) = corravg_count(idx_oc_cor_chi_Ni,k) + 1
        else ! oc_micro_correlations(1,inci,k) == undef_corr
           call hbuf_put_level('oc_corr_chi_Ni',0.0,1.,k)
        endif ! oc_micro_correlations(1,inci,k) /= undef_corr
        if ( oc_micro_correlations(idx_s,idx_rs,k) /= undef_corr ) then
           call hbuf_put_level('oc_corr_chi_rs',oc_micro_correlations(idx_s,idx_rs,k),1.,k)
           corravg_count(idx_oc_cor_chi_rs,k) = corravg_count(idx_oc_cor_chi_rs,k) + 1
        else ! oc_micro_correlations(1,iqs,k) == undef_corr
           call hbuf_put_level('oc_corr_chi_rs',0.0,1.,k)
        endif ! oc_micro_correlations(1,iqs,k) /= undef_corr
        if ( oc_micro_correlations(idx_s,idx_Ns,k) /= undef_corr ) then
           call hbuf_put_level('oc_corr_chi_Ns',oc_micro_correlations(idx_s,idx_Ns,k),1.,k)
           corravg_count(idx_oc_cor_chi_Ns,k) = corravg_count(idx_oc_cor_chi_Ns,k) + 1
        else ! oc_micro_correlations(1,ins,k) == undef_corr
           call hbuf_put_level('oc_corr_chi_Ns',0.0,1.,k)
        endif ! oc_micro_correlations(1,ins,k) /= undef_corr
      end if ! doicemicro

      if (dograupel) then
        if ( oc_micro_correlations(idx_s,idx_rg,k) /= undef_corr ) then
           call hbuf_put_level('oc_corr_chi_rg',oc_micro_correlations(idx_s,idx_rg,k),1.,k)
           corravg_count(idx_oc_cor_chi_rg,k) = corravg_count(idx_oc_cor_chi_rg,k) + 1
        else ! oc_micro_correlations(1,iqg,k) == undef_corr
           call hbuf_put_level('oc_corr_chi_rg',0.0,1.,k)
        endif ! oc_micro_correlations(1,iqg,k) /= undef_corr
        if ( oc_micro_correlations(idx_s,idx_Ng,k) /= undef_corr ) then
           call hbuf_put_level('oc_corr_chi_Ng',oc_micro_correlations(idx_s,idx_Ng,k),1.,k)
           corravg_count(idx_oc_cor_chi_Ng,k) = corravg_count(idx_oc_cor_chi_Ng,k) + 1
        else ! oc_micro_correlations(1,ing,k) == undef_corr
           call hbuf_put_level('oc_corr_chi_Ng',0.0,1.,k)
        endif ! oc_micro_correlations(1,ing,k) /= undef_corr
      end if ! dograupel

!-------------
! Vertical Velocity
!-------------
      if(dopredictNc) then
        if ( oc_micro_correlations(idx_w,idx_Nc,k) /= undef_corr ) then
           call hbuf_put_level('oc_corr_w_Nc',oc_micro_correlations(idx_w,idx_Nc,k),1.,k)
           corravg_count(idx_oc_cor_w_Nc,k) = corravg_count(idx_oc_cor_w_Nc,k) + 1
        else ! oc_micro_correlations(1,incl,k) == undef_corr
           call hbuf_put_level('oc_corr_w_Nc',0.0,1.,k)
        endif ! oc_micro_correlations(1,incl,k) /= undef_corr
      end if

      if (doprecip) then
      if ( oc_micro_correlations(idx_w,idx_rr,k) /= undef_corr ) then
           call hbuf_put_level('oc_corr_w_rr',oc_micro_correlations(idx_w,idx_rr,k),1.,k)
           corravg_count(idx_oc_cor_w_rr,k) = corravg_count(idx_oc_cor_w_rr,k) + 1
        else ! oc_micro_correlations(1,iqr,k) == undef_corr
           call hbuf_put_level('oc_corr_w_rr',0.0,1.,k)
        endif ! oc_micro_correlations(1,iqr,k) /= undef_corr
        if ( oc_micro_correlations(idx_w,idx_Nr,k) /= undef_corr ) then
           call hbuf_put_level('oc_corr_w_Nr',oc_micro_correlations(idx_w,idx_Nr,k),1.,k)
           corravg_count(idx_oc_cor_w_nr,k) = corravg_count(idx_oc_cor_w_nr,k) + 1
        else ! oc_micro_correlations(1,inr,k) == undef_corr
           call hbuf_put_level('oc_corr_w_Nr',0.0,1.,k)
        endif ! ! oc_micro_correlations(1,inr,k) /= undef_corr
      end if ! doprecip

      if (doicemicro) then
      if ( oc_micro_correlations(idx_w,idx_ri,k) /= undef_corr ) then
           call hbuf_put_level('oc_corr_w_ri',oc_micro_correlations(idx_w,idx_ri,k),1.,k)
           corravg_count(idx_oc_cor_w_ri,k) = corravg_count(idx_oc_cor_w_ri,k) + 1
        else ! oc_micro_correlations(1,iqci,k) == undef_corr
           call hbuf_put_level('oc_corr_w_ri',0.0,1.,k)
        endif ! oc_micro_correlations(1,iqci,k) /= undef_corr
        if ( oc_micro_correlations(idx_w,idx_Ni,k) /= undef_corr ) then
           call hbuf_put_level('oc_corr_w_Ni',oc_micro_correlations(idx_w,idx_Ni,k),1.,k)
           corravg_count(idx_oc_cor_w_Ni,k) = corravg_count(idx_oc_cor_w_Ni,k) + 1
        else ! oc_micro_correlations(1,inci,k) == undef_corr
           call hbuf_put_level('oc_corr_w_Ni',0.0,1.,k)
        endif ! oc_micro_correlations(1,inci,k) /= undef_corr
        if ( oc_micro_correlations(idx_w,idx_rs,k) /= undef_corr ) then
           call hbuf_put_level('oc_corr_w_rs',oc_micro_correlations(idx_w,idx_rs,k),1.,k)
           corravg_count(idx_oc_cor_w_rs,k) = corravg_count(idx_oc_cor_w_rs,k) + 1
        else ! oc_micro_correlations(1,iqs,k) == undef_corr
           call hbuf_put_level('oc_corr_w_rs',0.0,1.,k)
        endif ! oc_micro_correlations(1,iqs,k) /= undef_corr
        if ( oc_micro_correlations(idx_w,idx_Ns,k) /= undef_corr ) then
           call hbuf_put_level('oc_corr_w_Ns',oc_micro_correlations(idx_w,idx_Ns,k),1.,k)
           corravg_count(idx_oc_cor_w_Ns,k) = corravg_count(idx_oc_cor_w_Ns,k) + 1
        else ! oc_micro_correlations(1,ins,k) == undef_corr
           call hbuf_put_level('oc_corr_w_Ns',0.0,1.,k)
        endif ! oc_micro_correlations(1,ins,k) /= undef_corr
      end if ! doicemicro

      if (dograupel) then
        if ( oc_micro_correlations(idx_w,idx_rg,k) /= undef_corr ) then
           call hbuf_put_level('oc_corr_w_rg',oc_micro_correlations(idx_w,idx_rg,k),1.,k)
           corravg_count(idx_oc_cor_w_rg,k) = corravg_count(idx_oc_cor_w_rg,k) + 1
        else ! oc_micro_correlations(1,iqg,k) == undef_corr
           call hbuf_put_level('oc_corr_w_rg',0.0,1.,k)
        endif ! oc_micro_correlations(1,iqg,k) /= undef_corr
        if ( oc_micro_correlations(idx_w,idx_Ng,k) /= undef_corr ) then
           call hbuf_put_level('oc_corr_w_Ng',oc_micro_correlations(idx_w,idx_Ng,k),1.,k)
           corravg_count(idx_oc_cor_w_Ng,k) = corravg_count(idx_oc_cor_w_Ng,k) + 1
        else ! oc_micro_correlations(1,ing,k) == undef_corr
           call hbuf_put_level('oc_corr_w_Ng',0.0,1.,k)
        endif ! oc_micro_correlations(1,ing,k) /= undef_corr
      end if ! dograupel

!-------------
! Cloud droplet number concentration
!-------------
      if (dopredictNc) then
     
        if (doprecip) then
           if ( oc_micro_correlations(idx_Nc,idx_rr,k) /= undef_corr ) then
              call hbuf_put_level('oc_corr_Nc_rr',oc_micro_correlations(idx_Nc,idx_rr,k),1.,k)
              corravg_count(idx_oc_cor_Nc_rr,k) = corravg_count(idx_oc_cor_Nc_rr,k) + 1
           else ! oc_micro_correlations(incl,iqr,k) == undef_corr
              call hbuf_put_level('oc_corr_Nc_rr',0.0,1.,k)
           endif ! oc_micro_correlations(incl,iqr,k) /= undef_corr
           if ( oc_micro_correlations(idx_Nc,idx_Nr,k) /= undef_corr ) then
              call hbuf_put_level('oc_corr_Nc_Nr',oc_micro_correlations(idx_Nc,idx_Nr,k),1.,k)
              corravg_count(idx_oc_cor_Nc_nr,k) = corravg_count(idx_oc_cor_Nc_nr,k) + 1
           else ! oc_micro_correlations(incl,inr,k) == undef_corr
              call hbuf_put_level('oc_corr_Nc_Nr',0.0,1.,k)
           endif ! oc_micro_correlations(incl,inr,k) /= undef_corr
        end if
     
        if (doicemicro) then
           if ( oc_micro_correlations(idx_Nc,idx_ri,k) /= undef_corr ) then
              call hbuf_put_level('oc_corr_Nc_ri',oc_micro_correlations(idx_Nc,idx_ri,k),1.,k) 
              corravg_count(idx_oc_cor_Nc_ri,k) = corravg_count(idx_oc_cor_Nc_ri,k) + 1
           else ! oc_micro_correlations(incl,iqci,k) == undef_corr
              call hbuf_put_level('oc_corr_Nc_ri',0.0,1.,k)
           endif ! oc_micro_correlations(incl,iqci,k) /= undef_corr
           if ( oc_micro_correlations(idx_Nc,idx_Ni,k) /= undef_corr ) then
              call hbuf_put_level('oc_corr_Nc_Ni',oc_micro_correlations(idx_Nc,idx_Ni,k),1.,k)
              corravg_count(idx_oc_cor_Nc_Ni,k) = corravg_count(idx_oc_cor_Nc_Ni,k) + 1
           else ! oc_micro_correlations(incl,inci,k) == undef_corr
              call hbuf_put_level('oc_corr_Nc_Ni',0.0,1.,k)
           endif ! oc_micro_correlations(incl,inci,k) /= undef_corr
           if ( oc_micro_correlations(idx_Nc,idx_rs,k) /= undef_corr ) then
              call hbuf_put_level('oc_corr_Nc_rs',oc_micro_correlations(idx_Nc,idx_rs,k),1.,k)
              corravg_count(idx_oc_cor_Nc_rs,k) = corravg_count(idx_oc_cor_Nc_rs,k) + 1
           else ! oc_micro_correlations(incl,iqs,k) == undef_corr
              call hbuf_put_level('oc_corr_Nc_rs',0.0,1.,k)
           endif ! oc_micro_correlations(incl,iqs,k) /= undef_corr
           if ( oc_micro_correlations(idx_Nc,idx_Ns,k) /= undef_corr ) then
              call hbuf_put_level('oc_corr_Nc_Ns',oc_micro_correlations(idx_Nc,idx_Ns,k),1.,k)
              corravg_count(idx_oc_cor_Nc_Ns,k) = corravg_count(idx_oc_cor_Nc_Ns,k) + 1
           else ! oc_micro_correlations(incl,ins,k) == undef_corr
              call hbuf_put_level('oc_corr_Nc_Ns',0.0,1.,k)
           endif ! oc_micro_correlations(incl,ins,k) /= undef_corr
        end if    

        if (dograupel) then
           if ( oc_micro_correlations(idx_Nc,idx_rg,k) /= undef_corr ) then
              call hbuf_put_level('oc_corr_Nc_rg',oc_micro_correlations(idx_Nc,idx_rg,k),1.,k)
              corravg_count(idx_oc_cor_Nc_rg,k) = corravg_count(idx_oc_cor_Nc_rg,k) + 1
           else ! oc_micro_correlations(incl,iqg,k) == undef_corr
              call hbuf_put_level('oc_corr_Nc_rg',0.0,1.,k)
           endif ! oc_micro_correlations(incl,iqg,k) /= undef_corr
           if ( oc_micro_correlations(idx_Nc,idx_Ng,k) /= undef_corr ) then
              call hbuf_put_level('oc_corr_Nc_Ng',oc_micro_correlations(idx_Nc,idx_Ng,k),1.,k)
              corravg_count(idx_oc_cor_Nc_Ng,k) = corravg_count(idx_oc_cor_Nc_Ng,k) + 1
           else ! oc_micro_correlations(incl,ing,k) == undef_corr
              call hbuf_put_level('oc_corr_Nc_Ng',0.0,1.,k)
           endif ! oc_micro_correlations(incl,ing,k) /= undef_corr
        end if   

      end if !dopredictNc
  
!-------------
! Rainwater mixing ratio
!-------------
      if (doprecip) then 
      
         if ( oc_micro_correlations(idx_rr,idx_Nr,k) /= undef_corr ) then
             call hbuf_put_level('oc_corr_rr_Nr',oc_micro_correlations(idx_rr,idx_Nr,k),1.,k)
             corravg_count(idx_oc_cor_rr_nr,k) = corravg_count(idx_oc_cor_rr_nr,k) + 1
          else ! oc_micro_correlations(iqr,inr,k) == undef_corr
             call hbuf_put_level('oc_corr_rr_Nr',0.0,1.,k)
          endif ! oc_micro_correlations(iqr,inr,k) /= undef_corr

        if (doicemicro) then
           if ( oc_micro_correlations(idx_rr,idx_ri,k) /= undef_corr ) then
              call hbuf_put_level('oc_corr_rr_ri',oc_micro_correlations(idx_rr,idx_ri,k),1.,k)
              corravg_count(idx_oc_cor_rr_ri,k) = corravg_count(idx_oc_cor_rr_ri,k) + 1
           else ! oc_micro_correlations(iqr,iqci,k) == undef_corr
              call hbuf_put_level('oc_corr_rr_ri',0.0,1.,k)
           endif ! oc_micro_correlations(iqr,iqci,k) /= undef_corr
           if ( oc_micro_correlations(idx_rr,idx_Ni,k) /= undef_corr ) then
              call hbuf_put_level('oc_corr_rr_Ni',oc_micro_correlations(idx_rr,idx_Ni,k),1.,k)
              corravg_count(idx_oc_cor_rr_Ni,k) = corravg_count(idx_oc_cor_rr_Ni,k) + 1
           else ! oc_micro_correlations(iqr,inci,k) == undef_corr
              call hbuf_put_level('oc_corr_rr_Ni',0.0,1.,k)
           endif ! oc_micro_correlations(iqr,inci,k) /= undef_corr
           if ( oc_micro_correlations(idx_rr,idx_rs,k) /= undef_corr ) then
              call hbuf_put_level('oc_corr_rr_rs',oc_micro_correlations(idx_rr,idx_rs,k),1.,k)
              corravg_count(idx_oc_cor_rr_rs,k) = corravg_count(idx_oc_cor_rr_rs,k) + 1
           else ! oc_micro_correlations(iqr,iqs,k) == undef_corr
              call hbuf_put_level('oc_corr_rr_rs',0.0,1.,k)
           endif ! oc_micro_correlations(iqr,iqs,k) /= undef_corr
           if ( oc_micro_correlations(idx_rr,idx_Ns,k) /= undef_corr ) then
              call hbuf_put_level('oc_corr_rr_Ns',oc_micro_correlations(idx_rr,idx_Ns,k),1.,k)
              corravg_count(idx_oc_cor_rr_Ns,k) = corravg_count(idx_oc_cor_rr_Ns,k) + 1
           else ! oc_micro_correlations(iqr,ins,k) == undef_corr
              call hbuf_put_level('oc_corr_rr_Ns',0.0,1.,k)
           endif ! oc_micro_correlations(iqr,ins,k) /= undef_corr
        end if    

        if (dograupel) then
           if ( oc_micro_correlations(idx_rr,idx_rg,k) /= undef_corr ) then
              call hbuf_put_level('oc_corr_rr_rg',oc_micro_correlations(idx_rr,idx_rg,k),1.,k)
              corravg_count(idx_oc_cor_rr_rg,k) = corravg_count(idx_oc_cor_rr_rg,k) + 1
           else ! oc_micro_correlations(iqr,iqg,k) == undef_corr
              call hbuf_put_level('oc_corr_rr_rg',0.0,1.,k)
           endif ! oc_micro_correlations(iqr,iqg,k) /= undef_corr
           if ( oc_micro_correlations(idx_rr,idx_Ng,k) /= undef_corr ) then
              call hbuf_put_level('oc_corr_rr_Ng',oc_micro_correlations(idx_rr,idx_Ng,k),1.,k)
              corravg_count(idx_oc_cor_rr_Ng,k) = corravg_count(idx_oc_cor_rr_Ng,k) + 1
           else ! oc_micro_correlations(iqr,ing,k) == undef_corr
              call hbuf_put_level('oc_corr_rr_Ng',0.0,1.,k)
           endif ! oc_micro_correlations(iqr,ing,k) /= undef_corr
        end if
!-------------
! Rainwater number concentration
!-------------

        if (doicemicro) then
           if ( oc_micro_correlations(idx_Nr,idx_ri,k) /= undef_corr ) then
              call hbuf_put_level('oc_corr_Nr_ri',oc_micro_correlations(idx_Nr,idx_ri,k),1.,k)
              corravg_count(idx_oc_cor_Nr_ri,k) = corravg_count(idx_oc_cor_Nr_ri,k) + 1
           else ! oc_micro_correlations(inr,iqci,k) == undef_corr
              call hbuf_put_level('oc_corr_Nr_ri',0.0,1.,k)
           endif ! oc_micro_correlations(inr,iqci,k) /= undef_corr
           if ( oc_micro_correlations(idx_Nr,idx_Ni,k) /= undef_corr ) then
              call hbuf_put_level('oc_corr_Nr_Ni',oc_micro_correlations(idx_Nr,idx_Ni,k),1.,k)
              corravg_count(idx_oc_cor_Nr_Ni,k) = corravg_count(idx_oc_cor_Nr_Ni,k) + 1
           else ! oc_micro_correlations(inr,inci,k) == undef_corr
              call hbuf_put_level('oc_corr_Nr_Ni',0.0,1.,k)
           endif ! oc_micro_correlations(inr,inci,k) /= undef_corr
           if ( oc_micro_correlations(idx_Nr,idx_rs,k) /= undef_corr ) then
              call hbuf_put_level('oc_corr_Nr_rs',oc_micro_correlations(idx_Nr,idx_rs,k),1.,k)
              corravg_count(idx_oc_cor_Nr_rs,k) = corravg_count(idx_oc_cor_Nr_rs,k) + 1
           else ! oc_micro_correlations(inr,iqs,k) == undef_corr
              call hbuf_put_level('oc_corr_Nr_rs',0.0,1.,k)
           endif ! oc_micro_correlations(inr,iqs,k) /= undef_corr
           if ( oc_micro_correlations(idx_Nr,idx_Ns,k) /= undef_corr ) then
              call hbuf_put_level('oc_corr_Nr_Ns',oc_micro_correlations(idx_Nr,idx_Ns,k),1.,k)
              corravg_count(idx_oc_cor_Nr_Ns,k) = corravg_count(idx_oc_cor_Nr_Ns,k) + 1
           else ! oc_micro_correlations(inr,ins,k) == undef_corr
              call hbuf_put_level('oc_corr_Nr_Ns',0.0,1.,k)
           endif ! oc_micro_correlations(inr,ins,k) /= undef_corr
        end if

        if (dograupel) then
        if ( oc_micro_correlations(idx_Nr,idx_rg,k) /= undef_corr ) then
              call hbuf_put_level('oc_corr_Nr_rg',oc_micro_correlations(idx_Nr,idx_rg,k),1.,k)
              corravg_count(idx_oc_cor_Nr_rg,k) = corravg_count(idx_oc_cor_Nr_rg,k) + 1
           else ! oc_micro_correlations(inr,iqg,k) == undef_corr
              call hbuf_put_level('oc_corr_Nr_rg',0.0,1.,k)
           endif ! oc_micro_correlations(inr,iqg,k) /= undef_corr
           if ( oc_micro_correlations(idx_Nr,idx_Ng,k) /= undef_corr ) then
              call hbuf_put_level('oc_corr_Nr_Ng',oc_micro_correlations(idx_Nr,idx_Ng,k),1.,k)
              corravg_count(idx_oc_cor_Nr_Ng,k) = corravg_count(idx_oc_cor_Nr_Ng,k) + 1
           else ! oc_micro_correlations(inr,ing,k) == undef_corr
              call hbuf_put_level('oc_corr_Nr_Ng',0.0,1.,k)
           endif ! oc_micro_correlations(inr,ing,k) /= undef_corr
        end if

      end if !doprecip
!-------------
! Cloud-ice mixing ratio
!-------------

      if (doicemicro) then   
         if ( oc_micro_correlations(idx_ri,idx_Ni,k) /= undef_corr ) then
            call hbuf_put_level('oc_corr_ri_Ni',oc_micro_correlations(idx_ri,idx_Ni,k),1.,k)
            corravg_count(idx_oc_cor_ri_Ni,k) = corravg_count(idx_oc_cor_ri_Ni,k) + 1
         else ! oc_micro_correlations(iqci,inci,k) == undef_corr
            call hbuf_put_level('oc_corr_ri_Ni',0.0,1.,k)
         endif ! oc_micro_correlations(iqci,inci,k) /= undef_corr
         if ( oc_micro_correlations(idx_ri,idx_rs,k) /= undef_corr ) then
            call hbuf_put_level('oc_corr_ri_rs',oc_micro_correlations(idx_ri,idx_rs,k),1.,k)
            corravg_count(idx_oc_cor_ri_rs,k) = corravg_count(idx_oc_cor_ri_rs,k) + 1
         else ! oc_micro_correlations(iqci,iqs,k) == undef_corr
            call hbuf_put_level('oc_corr_ri_rs',0.0,1.,k)
         endif ! oc_micro_correlations(iqci,iqs,k) /= undef_corr
         if ( oc_micro_correlations(idx_ri,idx_Ns,k) /= undef_corr ) then
            call hbuf_put_level('oc_corr_ri_Ns',oc_micro_correlations(idx_ri,idx_Ns,k),1.,k)
            corravg_count(idx_oc_cor_ri_Ns,k) = corravg_count(idx_oc_cor_ri_Ns,k) + 1
         else ! oc_micro_correlations(iqci,ins,k) == undef_corr
            call hbuf_put_level('oc_corr_ri_Ns',0.0,1.,k)
         endif ! oc_micro_correlations(iqci,ins,k) /= undef_corr
      end if   

      if (dograupel) then
         if ( oc_micro_correlations(idx_ri,idx_rg,k) /= undef_corr ) then
            call hbuf_put_level('oc_corr_ri_rg',oc_micro_correlations(idx_ri,idx_rg,k),1.,k)
            corravg_count(idx_oc_cor_ri_rg,k) = corravg_count(idx_oc_cor_ri_rg,k) + 1
         else ! oc_micro_correlations(iqci,iqg,k) == undef_corr
            call hbuf_put_level('oc_corr_ri_rg',0.0,1.,k)
         endif ! oc_micro_correlations(iqci,iqg,k) /= undef_corr
         if ( oc_micro_correlations(idx_ri,idx_Ng,k) /= undef_corr ) then
            call hbuf_put_level('oc_corr_ri_Ng',oc_micro_correlations(idx_ri,idx_Ng,k),1.,k)
            corravg_count(idx_oc_cor_ri_Ng,k) = corravg_count(idx_oc_cor_ri_Ng,k) + 1
         else ! oc_micro_correlations(iqci,ing,k) == undef_corr
            call hbuf_put_level('oc_corr_ri_Ng',0.0,1.,k)
         endif ! oc_micro_correlations(iqci,ing,k) /= undef_corr
      end if
 
!-------------
! Cloud-ice number concentration
!-------------
      if (doicemicro) then   
         if ( oc_micro_correlations(idx_Ni,idx_rs,k) /= undef_corr ) then
            call hbuf_put_level('oc_corr_Ni_rs',oc_micro_correlations(idx_Ni,idx_rs,k),1.,k)
            corravg_count(idx_oc_cor_Ni_rs,k) = corravg_count(idx_oc_cor_Ni_rs,k) + 1
         else ! oc_micro_correlations(inci,iqs,k) == undef_corr
            call hbuf_put_level('oc_corr_Ni_rs',0.0,1.,k)
         endif ! oc_micro_correlations(inci,iqs,k) /= undef_corr
         if ( oc_micro_correlations(idx_Ni,idx_Ns,k) /= undef_corr ) then
            call hbuf_put_level('oc_corr_Ni_Ns',oc_micro_correlations(idx_Ni,idx_Ns,k),1.,k)
            corravg_count(idx_oc_cor_Ni_Ns,k) = corravg_count(idx_oc_cor_Ni_Ns,k) + 1
         else ! oc_micro_correlations(inci,ins,k) == undef_corr
            call hbuf_put_level('oc_corr_Ni_Ns',0.0,1.,k)
         endif ! oc_micro_correlations(inci,ins,k) /= undef_corr
      end if

      if (dograupel) then     
         if ( oc_micro_correlations(idx_Ni,idx_rg,k) /= undef_corr ) then
            call hbuf_put_level('oc_corr_Ni_rg',oc_micro_correlations(idx_Ni,idx_rg,k),1.,k)
            corravg_count(idx_oc_cor_Ni_rg,k) = corravg_count(idx_oc_cor_Ni_rg,k) + 1
         else ! oc_micro_correlations(inci,iqg,k) == undef_corr
            call hbuf_put_level('oc_corr_Ni_rg',0.0,1.,k)
         endif ! oc_micro_correlations(inci,iqg,k) /= undef_corr
         if ( oc_micro_correlations(idx_Ni,idx_Ng,k) /= undef_corr ) then
            call hbuf_put_level('oc_corr_Ni_Ng',oc_micro_correlations(idx_Ni,idx_Ng,k),1.,k)
            corravg_count(idx_oc_cor_Ni_Ng,k) = corravg_count(idx_oc_cor_Ni_Ng,k) + 1
         else ! oc_micro_correlations(inci,ing,k) == undef_corr
            call hbuf_put_level('oc_corr_Ni_Ng',0.0,1.,k)
         endif ! oc_micro_correlations(inci,ing,k) /= undef_corr
      end if 
!-------------
! Snow mixing ratio
!-------------

      if (doicemicro) then
         if ( oc_micro_correlations(idx_rs,idx_Ns,k) /= undef_corr ) then
            call hbuf_put_level('oc_corr_rs_Ns',oc_micro_correlations(idx_rs,idx_Ns,k),1.,k)
            corravg_count(idx_oc_cor_rs_Ns,k) = corravg_count(idx_oc_cor_rs_Ns,k) + 1
         else ! oc_micro_correlations(iqs,ins,k) == undef_corr
            call hbuf_put_level('oc_corr_rs_Ns',0.0,1.,k)
         endif ! oc_micro_correlations(iqs,ins,k) /= undef_corr
       
          if (dograupel) then   
             if ( oc_micro_correlations(idx_rs,idx_rg,k) /= undef_corr ) then
                call hbuf_put_level('oc_corr_rs_rg',oc_micro_correlations(idx_rs,idx_rg,k),1.,k)
                corravg_count(idx_oc_cor_rs_rg,k) = corravg_count(idx_oc_cor_rs_rg,k) + 1
             else ! oc_micro_correlations(iqs,iqg,k) == undef_corr
                call hbuf_put_level('oc_corr_rs_rg',0.0,1.,k)
             endif ! oc_micro_correlations(iqs,iqg,k) /= undef_corr
             if ( oc_micro_correlations(idx_rs,idx_Ng,k) /= undef_corr ) then
                call hbuf_put_level('oc_corr_rs_Ng',oc_micro_correlations(idx_rs,idx_Ng,k),1.,k)
                corravg_count(idx_oc_cor_rs_Ng,k) = corravg_count(idx_oc_cor_rs_Ng,k) + 1
             else ! oc_micro_correlations(iqs,ing,k) == undef_corr
                call hbuf_put_level('oc_corr_rs_Ng',0.0,1.,k)
             endif ! oc_micro_correlations(iqs,ing,k) /= undef_corr
!-------------
! Snow number concentration
!-------------
             if ( oc_micro_correlations(idx_Ns,idx_rg,k) /= undef_corr ) then
                call hbuf_put_level('oc_corr_Ns_rg',oc_micro_correlations(idx_Ns,idx_rg,k),1.,k)
                corravg_count(idx_oc_cor_Ns_rg,k) = corravg_count(idx_oc_cor_Ns_rg,k) + 1
             else ! oc_micro_correlations(ins,iqg,k) == undef_corr
                call hbuf_put_level('oc_corr_Ns_rg',0.0,1.,k)
             endif ! oc_micro_correlations(ins,iqg,k) /= undef_corr
             if ( oc_micro_correlations(idx_Ns,idx_Ng,k) /= undef_corr ) then
                call hbuf_put_level('oc_corr_Ns_Ng',oc_micro_correlations(idx_Ns,idx_Ng,k),1.,k)
                corravg_count(idx_oc_cor_Ns_Ng,k) = corravg_count(idx_oc_cor_Ns_Ng,k) + 1
             else ! oc_micro_correlations(ins,ing,k) == undef_corr
                call hbuf_put_level('oc_corr_Ns_Ng',0.0,1.,k)
             endif ! oc_micro_correlations(ins,ing,k) /= undef_corr
          end if
      end if

!-------------
! Graupel mixing ratio / number concentration
!-------------
      if (dograupel) then      
        if ( oc_micro_correlations(idx_rg,idx_Ng,k) /= undef_corr ) then
            call hbuf_put_level('oc_corr_rg_Ng',oc_micro_correlations(idx_rg,idx_Ng,k),1.,k)
            corravg_count(idx_oc_cor_rg_Ng,k) = corravg_count(idx_oc_cor_rg_Ng,k) + 1
         else ! oc_micro_correlations(iqg,ing,k) == undef_corr
            call hbuf_put_level('oc_corr_rg_Ng',0.0,1.,k)
         endif ! oc_micro_correlations(iqg,ing,k) /= undef_corr
      endif

    end do ! do k = 1, nzm, 1

!---------------------------------
! Covariances
!--------------------------------

!---------------------------------
! chi
!--------------------------------

call hbuf_put('covarnce_chi_w',micro_covarnce(idx_s,idx_w,:),1.)

if(dopredictNc) then
  call hbuf_put('covarnce_chi_Nc',micro_covarnce(idx_s,idx_Nc,:),1.)
endif !dopredictNc

if (doprecip) then
  call hbuf_put('covarnce_chi_rr',micro_covarnce(idx_s,idx_rr,:),1.)
  call hbuf_put('covarnce_chi_Nr',micro_covarnce(idx_s,idx_Nr,:),1.)
end if ! doprecip

if (doicemicro) then
  call hbuf_put('covarnce_chi_ri',micro_covarnce(idx_s,idx_ri,:),1.)
  call hbuf_put('covarnce_chi_Ni',micro_covarnce(idx_s,idx_Ni,:),1.)
  call hbuf_put('covarnce_chi_rs',micro_covarnce(idx_s,idx_rs,:),1.)
  call hbuf_put('covarnce_chi_Ns',micro_covarnce(idx_s,idx_Ns,:),1.)
end if ! doicemicro

if (dograupel) then
  call hbuf_put('covarnce_chi_rg',micro_covarnce(idx_s,idx_rg,:),1.)
  call hbuf_put('covarnce_chi_Ng',micro_covarnce(idx_s,idx_Ng,:),1.)
end if ! dograupel


!---------------------------------
! Vertical velocity
!--------------------------------
if(dopredictNc) then
  call hbuf_put('covarnce_w_Nc',micro_covarnce(idx_w,idx_Nc,:),1.)
endif !dopredictNc

if (doprecip) then
  call hbuf_put('covarnce_w_rr',micro_covarnce(idx_w,idx_rr,:),1.)
  call hbuf_put('covarnce_w_Nr',micro_covarnce(idx_w,idx_Nr,:),1.)
end if ! doprecip

if (doicemicro) then
  call hbuf_put('covarnce_w_ri',micro_covarnce(idx_w,idx_ri,:),1.)
  call hbuf_put('covarnce_w_Ni',micro_covarnce(idx_w,idx_Ni,:),1.)
  call hbuf_put('covarnce_w_rs',micro_covarnce(idx_w,idx_rs,:),1.)
  call hbuf_put('covarnce_w_Ns',micro_covarnce(idx_w,idx_Ns,:),1.)
end if ! doicemicro

if (dograupel) then
  call hbuf_put('covarnce_w_rg',micro_covarnce(idx_w,idx_rg,:),1.)
  call hbuf_put('covarnce_w_Ng',micro_covarnce(idx_w,idx_Ng,:),1.)
end if ! dograupel

!---------------------------------
! Cloud liquid number concentration
!--------------------------------

if(dopredictNc) then  
  if (doprecip) then
    call hbuf_put('covarnce_Nc_rr',micro_covarnce(idx_Nc,idx_rr,:),1.)
    call hbuf_put('covarnce_Nc_Nr',micro_covarnce(idx_Nc,idx_Nr,:),1.)
  end if ! doprecip
     
  if (doicemicro) then
    call hbuf_put('covarnce_Nc_ri',micro_covarnce(idx_Nc,idx_ri,:),1.) 
    call hbuf_put('covarnce_Nc_Ni',micro_covarnce(idx_Nc,idx_Ni,:),1.)
    call hbuf_put('covarnce_Nc_rs',micro_covarnce(idx_Nc,idx_rs,:),1.)
    call hbuf_put('covarnce_Nc_Ns',micro_covarnce(idx_Nc,idx_Ns,:),1.)
  end if ! doicemicro    

  if (dograupel) then
    call hbuf_put('covarnce_Nc_rg',micro_covarnce(idx_Nc,idx_rg,:),1.)
    call hbuf_put('covarnce_Nc_Ng',micro_covarnce(idx_Nc,idx_Ng,:),1.)
  end if ! dograupel

end if !dopredictNc

!---------------------------------
! Rainwater mixing ratio
!--------------------------------
  
if (doprecip) then 
  call hbuf_put('covarnce_rr_Nr',micro_covarnce(idx_rr,idx_Nr,:),1.)

  if (doicemicro) then
    call hbuf_put('covarnce_rr_ri',micro_covarnce(idx_rr,idx_ri,:),1.)
    call hbuf_put('covarnce_rr_Ni',micro_covarnce(idx_rr,idx_Ni,:),1.)
    call hbuf_put('covarnce_rr_rs',micro_covarnce(idx_rr,idx_rs,:),1.)
    call hbuf_put('covarnce_rr_Ns',micro_covarnce(idx_rr,idx_Ns,:),1.)

    if (dograupel) then
      call hbuf_put('covarnce_rr_rg',micro_covarnce(idx_rr,idx_rg,:),1.)
      call hbuf_put('covarnce_rr_Ng',micro_covarnce(idx_rr,idx_Ng,:),1.)
    end if !dograupel

!---------------------------------
! Rainwater number concentration
!--------------------------------

    call hbuf_put('covarnce_Nr_ri',micro_covarnce(idx_Nr,idx_ri,:),1.)
    call hbuf_put('covarnce_Nr_Ni',micro_covarnce(idx_Nr,idx_Ni,:),1.)
    call hbuf_put('covarnce_Nr_rs',micro_covarnce(idx_Nr,idx_rs,:),1.)
    call hbuf_put('covarnce_Nr_Ns',micro_covarnce(idx_Nr,idx_Ns,:),1.)

    if (dograupel) then
      call hbuf_put('covarnce_Nr_rg',micro_covarnce(idx_Nr,idx_rg,:),1.)
      call hbuf_put('covarnce_Nr_Ng',micro_covarnce(idx_Nr,idx_Ng,:),1.)
    end if !dograupel
  
  end if !doicemicro

end if !doprecip

!---------------------------------
! Cloud-ice mixing ratio
!--------------------------------

if (doicemicro) then   
  call hbuf_put('covarnce_ri_Ni',micro_covarnce(idx_ri,idx_Ni,:),1.)
  call hbuf_put('covarnce_ri_rs',micro_covarnce(idx_ri,idx_rs,:),1.)
  call hbuf_put('covarnce_ri_Ns',micro_covarnce(idx_ri,idx_Ns,:),1.)

  if (dograupel) then
    call hbuf_put('covarnce_ri_rg',micro_covarnce(idx_ri,idx_rg,:),1.)
    call hbuf_put('covarnce_ri_Ng',micro_covarnce(idx_ri,idx_Ng,:),1.)
  end if !dograupel

!---------------------------------
! Cloud-ice number concentration
!--------------------------------
 
  call hbuf_put('covarnce_Ni_rs',micro_covarnce(idx_Ni,idx_rs,:),1.)
  call hbuf_put('covarnce_Ni_Ns',micro_covarnce(idx_Ni,idx_Ns,:),1.)
      
  if (dograupel) then     
    call hbuf_put('covarnce_Ni_rg',micro_covarnce(idx_Ni,idx_rg,:),1.)
    call hbuf_put('covarnce_Ni_Ng',micro_covarnce(idx_Ni,idx_Ng,:),1.)
  end if !dograupel

!---------------------------------
! Snow mixing ratio
!--------------------------------

  call hbuf_put('covarnce_rs_Ns',micro_covarnce(idx_rs,idx_Ns,:),1.)
       
  if (dograupel) then   
    call hbuf_put('covarnce_rs_rg',micro_covarnce(idx_rs,idx_rg,:),1.)
    call hbuf_put('covarnce_rs_Ng',micro_covarnce(idx_rs,idx_Ng,:),1.)

!---------------------------------
! Snow number concentration
!--------------------------------
    call hbuf_put('covarnce_Ns_rg',micro_covarnce(idx_Ns,idx_rg,:),1.)
    call hbuf_put('covarnce_Ns_Ng',micro_covarnce(idx_Ns,idx_Ng,:),1.)

!---------------------------------
! Graupel mixing ratio / number concentration
!--------------------------------
    call hbuf_put('covarnce_rg_Ng',micro_covarnce(idx_rg,idx_Ng,:),1.)
  endif !dograupel
endif !doicemicro


!---------------------------------
! In-cloud covariances
!--------------------------------

!---------------------------------
! chi
!--------------------------------

call hbuf_put('ic_covarnce_chi_w',ic_micro_covarnce(idx_s,idx_w,:),1.)

if(dopredictNc) then
  call hbuf_put('ic_covarnce_chi_Nc',ic_micro_covarnce(idx_s,idx_Nc,:),1.)
endif !dopredictNc

if (doprecip) then
  call hbuf_put('ic_covarnce_chi_rr',ic_micro_covarnce(idx_s,idx_rr,:),1.)
  call hbuf_put('ic_covarnce_chi_Nr',ic_micro_covarnce(idx_s,idx_Nr,:),1.)
end if ! doprecip

if (doicemicro) then
  call hbuf_put('ic_covarnce_chi_ri',ic_micro_covarnce(idx_s,idx_ri,:),1.)
  call hbuf_put('ic_covarnce_chi_Ni',ic_micro_covarnce(idx_s,idx_Ni,:),1.)
  call hbuf_put('ic_covarnce_chi_rs',ic_micro_covarnce(idx_s,idx_rs,:),1.)
  call hbuf_put('ic_covarnce_chi_Ns',ic_micro_covarnce(idx_s,idx_Ns,:),1.)
end if ! doicemicro

if (dograupel) then
  call hbuf_put('ic_covarnce_chi_rg',ic_micro_covarnce(idx_s,idx_rg,:),1.)
  call hbuf_put('ic_covarnce_chi_Ng',ic_micro_covarnce(idx_s,idx_Ng,:),1.)
end if ! dograupel


!---------------------------------
! Vertical velocity
!--------------------------------
if(dopredictNc) then
  call hbuf_put('ic_covarnce_w_Nc',ic_micro_covarnce(idx_w,idx_Nc,:),1.)
endif !dopredictNc

if (doprecip) then
  call hbuf_put('ic_covarnce_w_rr',ic_micro_covarnce(idx_w,idx_rr,:),1.)
  call hbuf_put('ic_covarnce_w_Nr',ic_micro_covarnce(idx_w,idx_Nr,:),1.)
end if ! doprecip

if (doicemicro) then
  call hbuf_put('ic_covarnce_w_ri',ic_micro_covarnce(idx_w,idx_ri,:),1.)
  call hbuf_put('ic_covarnce_w_Ni',ic_micro_covarnce(idx_w,idx_Ni,:),1.)
  call hbuf_put('ic_covarnce_w_rs',ic_micro_covarnce(idx_w,idx_rs,:),1.)
  call hbuf_put('ic_covarnce_w_Ns',ic_micro_covarnce(idx_w,idx_Ns,:),1.)
end if ! doicemicro

if (dograupel) then
  call hbuf_put('ic_covarnce_w_rg',ic_micro_covarnce(idx_w,idx_rg,:),1.)
  call hbuf_put('ic_covarnce_w_Ng',ic_micro_covarnce(idx_w,idx_Ng,:),1.)
end if ! dograupel

!---------------------------------
! Cloud liquid number concentration
!--------------------------------

if(dopredictNc) then  
  if (doprecip) then
    call hbuf_put('ic_covarnce_Nc_rr',ic_micro_covarnce(idx_Nc,idx_rr,:),1.)
    call hbuf_put('ic_covarnce_Nc_Nr',ic_micro_covarnce(idx_Nc,idx_Nr,:),1.)
  end if ! doprecip
     
  if (doicemicro) then
    call hbuf_put('ic_covarnce_Nc_ri',ic_micro_covarnce(idx_Nc,idx_ri,:),1.) 
    call hbuf_put('ic_covarnce_Nc_Ni',ic_micro_covarnce(idx_Nc,idx_Ni,:),1.)
    call hbuf_put('ic_covarnce_Nc_rs',ic_micro_covarnce(idx_Nc,idx_rs,:),1.)
    call hbuf_put('ic_covarnce_Nc_Ns',ic_micro_covarnce(idx_Nc,idx_Ns,:),1.)
  end if ! doicemicro    

  if (dograupel) then
    call hbuf_put('ic_covarnce_Nc_rg',ic_micro_covarnce(idx_Nc,idx_rg,:),1.)
    call hbuf_put('ic_covarnce_Nc_Ng',ic_micro_covarnce(idx_Nc,idx_Ng,:),1.)
  end if ! dograupel

end if !dopredictNc

!---------------------------------
! Rainwater mixing ratio
!--------------------------------
  
if (doprecip) then 
  call hbuf_put('ic_covarnce_rr_Nr',ic_micro_covarnce(idx_rr,idx_Nr,:),1.)

  if (doicemicro) then
    call hbuf_put('ic_covarnce_rr_ri',ic_micro_covarnce(idx_rr,idx_ri,:),1.)
    call hbuf_put('ic_covarnce_rr_Ni',ic_micro_covarnce(idx_rr,idx_Ni,:),1.)
    call hbuf_put('ic_covarnce_rr_rs',ic_micro_covarnce(idx_rr,idx_rs,:),1.)
    call hbuf_put('ic_covarnce_rr_Ns',ic_micro_covarnce(idx_rr,idx_Ns,:),1.)

    if (dograupel) then
      call hbuf_put('ic_covarnce_rr_rg',ic_micro_covarnce(idx_rr,idx_rg,:),1.)
      call hbuf_put('ic_covarnce_rr_Ng',ic_micro_covarnce(idx_rr,idx_Ng,:),1.)
    end if !dograupel

!---------------------------------
! Rainwater number concentration
!--------------------------------

    call hbuf_put('ic_covarnce_Nr_ri',ic_micro_covarnce(idx_Nr,idx_ri,:),1.)
    call hbuf_put('ic_covarnce_Nr_Ni',ic_micro_covarnce(idx_Nr,idx_Ni,:),1.)
    call hbuf_put('ic_covarnce_Nr_rs',ic_micro_covarnce(idx_Nr,idx_rs,:),1.)
    call hbuf_put('ic_covarnce_Nr_Ns',ic_micro_covarnce(idx_Nr,idx_Ns,:),1.)

    if (dograupel) then
      call hbuf_put('ic_covarnce_Nr_rg',ic_micro_covarnce(idx_Nr,idx_rg,:),1.)
      call hbuf_put('ic_covarnce_Nr_Ng',ic_micro_covarnce(idx_Nr,idx_Ng,:),1.)
    end if !dograupel
  
  end if !doicemicro

end if !doprecip

!---------------------------------
! Cloud-ice mixing ratio
!--------------------------------

if (doicemicro) then   
  call hbuf_put('ic_covarnce_ri_Ni',ic_micro_covarnce(idx_ri,idx_Ni,:),1.)
  call hbuf_put('ic_covarnce_ri_rs',ic_micro_covarnce(idx_ri,idx_rs,:),1.)
  call hbuf_put('ic_covarnce_ri_Ns',ic_micro_covarnce(idx_ri,idx_Ns,:),1.)

  if (dograupel) then
    call hbuf_put('ic_covarnce_ri_rg',ic_micro_covarnce(idx_ri,idx_rg,:),1.)
    call hbuf_put('ic_covarnce_ri_Ng',ic_micro_covarnce(idx_ri,idx_Ng,:),1.)
  end if !dograupel

!---------------------------------
! Cloud-ice number concentration
!--------------------------------
 
  call hbuf_put('ic_covarnce_Ni_rs',ic_micro_covarnce(idx_Ni,idx_rs,:),1.)
  call hbuf_put('ic_covarnce_Ni_Ns',ic_micro_covarnce(idx_Ni,idx_Ns,:),1.)
      
  if (dograupel) then     
    call hbuf_put('ic_covarnce_Ni_rg',ic_micro_covarnce(idx_Ni,idx_rg,:),1.)
    call hbuf_put('ic_covarnce_Ni_Ng',ic_micro_covarnce(idx_Ni,idx_Ng,:),1.)
  end if !dograupel

!---------------------------------
! Snow mixing ratio
!--------------------------------

  call hbuf_put('ic_covarnce_rs_Ns',ic_micro_covarnce(idx_rs,idx_Ns,:),1.)
       
  if (dograupel) then   
    call hbuf_put('ic_covarnce_rs_rg',ic_micro_covarnce(idx_rs,idx_rg,:),1.)
    call hbuf_put('ic_covarnce_rs_Ng',ic_micro_covarnce(idx_rs,idx_Ng,:),1.)

!---------------------------------
! Snow number concentration
!--------------------------------
    call hbuf_put('ic_covarnce_Ns_rg',ic_micro_covarnce(idx_Ns,idx_rg,:),1.)
    call hbuf_put('ic_covarnce_Ns_Ng',ic_micro_covarnce(idx_Ns,idx_Ng,:),1.)

!---------------------------------
! Graupel mixing ratio / number concentration
!--------------------------------
    call hbuf_put('ic_covarnce_rg_Ng',ic_micro_covarnce(idx_rg,idx_Ng,:),1.)
  endif !dograupel
endif !doicemicro

!---------------------------------
! Out of cloud covariances
!--------------------------------

!---------------------------------
! chi
!--------------------------------

call hbuf_put('oc_covarnce_chi_w',oc_micro_covarnce(idx_s,idx_w,:),1.)

if(dopredictNc) then
  call hbuf_put('oc_covarnce_chi_Nc',oc_micro_covarnce(idx_s,idx_Nc,:),1.)
endif !dopredictNc

if (doprecip) then
  call hbuf_put('oc_covarnce_chi_rr',oc_micro_covarnce(idx_s,idx_rr,:),1.)
  call hbuf_put('oc_covarnce_chi_Nr',oc_micro_covarnce(idx_s,idx_Nr,:),1.)
end if ! doprecip

if (doicemicro) then
  call hbuf_put('oc_covarnce_chi_ri',oc_micro_covarnce(idx_s,idx_ri,:),1.)
  call hbuf_put('oc_covarnce_chi_Ni',oc_micro_covarnce(idx_s,idx_Ni,:),1.)
  call hbuf_put('oc_covarnce_chi_rs',oc_micro_covarnce(idx_s,idx_rs,:),1.)
  call hbuf_put('oc_covarnce_chi_Ns',oc_micro_covarnce(idx_s,idx_Ns,:),1.)
end if ! doicemicro

if (dograupel) then
  call hbuf_put('oc_covarnce_chi_rg',oc_micro_covarnce(idx_s,idx_rg,:),1.)
  call hbuf_put('oc_covarnce_chi_Ng',oc_micro_covarnce(idx_s,idx_Ng,:),1.)
end if ! dograupel


!---------------------------------
! Vertical velocity
!--------------------------------
if(dopredictNc) then
  call hbuf_put('oc_covarnce_w_Nc',oc_micro_covarnce(idx_w,idx_Nc,:),1.)
endif !dopredictNc

if (doprecip) then
  call hbuf_put('oc_covarnce_w_rr',oc_micro_covarnce(idx_w,idx_rr,:),1.)
  call hbuf_put('oc_covarnce_w_Nr',oc_micro_covarnce(idx_w,idx_Nr,:),1.)
end if ! doprecip

if (doicemicro) then
  call hbuf_put('oc_covarnce_w_ri',oc_micro_covarnce(idx_w,idx_ri,:),1.)
  call hbuf_put('oc_covarnce_w_Ni',oc_micro_covarnce(idx_w,idx_Ni,:),1.)
  call hbuf_put('oc_covarnce_w_rs',oc_micro_covarnce(idx_w,idx_rs,:),1.)
  call hbuf_put('oc_covarnce_w_Ns',oc_micro_covarnce(idx_w,idx_Ns,:),1.)
end if ! doicemicro

if (dograupel) then
  call hbuf_put('oc_covarnce_w_rg',oc_micro_covarnce(idx_w,idx_rg,:),1.)
  call hbuf_put('oc_covarnce_w_Ng',oc_micro_covarnce(idx_w,idx_Ng,:),1.)
end if ! dograupel

!---------------------------------
! Cloud liquid number concentration
!--------------------------------

if(dopredictNc) then  
  if (doprecip) then
    call hbuf_put('oc_covarnce_Nc_rr',oc_micro_covarnce(idx_Nc,idx_rr,:),1.)
    call hbuf_put('oc_covarnce_Nc_Nr',oc_micro_covarnce(idx_Nc,idx_Nr,:),1.)
  end if ! doprecip
     
  if (doicemicro) then
    call hbuf_put('oc_covarnce_Nc_ri',oc_micro_covarnce(idx_Nc,idx_ri,:),1.) 
    call hbuf_put('oc_covarnce_Nc_Ni',oc_micro_covarnce(idx_Nc,idx_Ni,:),1.)
    call hbuf_put('oc_covarnce_Nc_rs',oc_micro_covarnce(idx_Nc,idx_rs,:),1.)
    call hbuf_put('oc_covarnce_Nc_Ns',oc_micro_covarnce(idx_Nc,idx_Ns,:),1.)
  end if ! doicemicro    

  if (dograupel) then
    call hbuf_put('oc_covarnce_Nc_rg',oc_micro_covarnce(idx_Nc,idx_rg,:),1.)
    call hbuf_put('oc_covarnce_Nc_Ng',oc_micro_covarnce(idx_Nc,idx_Ng,:),1.)
  end if ! dograupel

end if !dopredictNc

!---------------------------------
! Rainwater mixing ratio
!--------------------------------
  
if (doprecip) then 
  call hbuf_put('oc_covarnce_rr_Nr',oc_micro_covarnce(idx_rr,idx_Nr,:),1.)

  if (doicemicro) then
    call hbuf_put('oc_covarnce_rr_ri',oc_micro_covarnce(idx_rr,idx_ri,:),1.)
    call hbuf_put('oc_covarnce_rr_Ni',oc_micro_covarnce(idx_rr,idx_Ni,:),1.)
    call hbuf_put('oc_covarnce_rr_rs',oc_micro_covarnce(idx_rr,idx_rs,:),1.)
    call hbuf_put('oc_covarnce_rr_Ns',oc_micro_covarnce(idx_rr,idx_Ns,:),1.)

    if (dograupel) then
      call hbuf_put('oc_covarnce_rr_rg',oc_micro_covarnce(idx_rr,idx_rg,:),1.)
      call hbuf_put('oc_covarnce_rr_Ng',oc_micro_covarnce(idx_rr,idx_Ng,:),1.)
    end if !dograupel

!---------------------------------
! Rainwater number concentration
!--------------------------------

    call hbuf_put('oc_covarnce_Nr_ri',oc_micro_covarnce(idx_Nr,idx_ri,:),1.)
    call hbuf_put('oc_covarnce_Nr_Ni',oc_micro_covarnce(idx_Nr,idx_Ni,:),1.)
    call hbuf_put('oc_covarnce_Nr_rs',oc_micro_covarnce(idx_Nr,idx_rs,:),1.)
    call hbuf_put('oc_covarnce_Nr_Ns',oc_micro_covarnce(idx_Nr,idx_Ns,:),1.)

    if (dograupel) then
      call hbuf_put('oc_covarnce_Nr_rg',oc_micro_covarnce(idx_Nr,idx_rg,:),1.)
      call hbuf_put('oc_covarnce_Nr_Ng',oc_micro_covarnce(idx_Nr,idx_Ng,:),1.)
    end if !dograupel
  
  end if !doicemicro

end if !doprecip

!---------------------------------
! Cloud-ice mixing ratio
!--------------------------------

if (doicemicro) then   
  call hbuf_put('oc_covarnce_ri_Ni',oc_micro_covarnce(idx_ri,idx_Ni,:),1.)
  call hbuf_put('oc_covarnce_ri_rs',oc_micro_covarnce(idx_ri,idx_rs,:),1.)
  call hbuf_put('oc_covarnce_ri_Ns',oc_micro_covarnce(idx_ri,idx_Ns,:),1.)

  if (dograupel) then
    call hbuf_put('oc_covarnce_ri_rg',oc_micro_covarnce(idx_ri,idx_rg,:),1.)
    call hbuf_put('oc_covarnce_ri_Ng',oc_micro_covarnce(idx_ri,idx_Ng,:),1.)
  end if !dograupel

!---------------------------------
! Cloud-ice number concentration
!--------------------------------
 
  call hbuf_put('oc_covarnce_Ni_rs',oc_micro_covarnce(idx_Ni,idx_rs,:),1.)
  call hbuf_put('oc_covarnce_Ni_Ns',oc_micro_covarnce(idx_Ni,idx_Ns,:),1.)
      
  if (dograupel) then     
    call hbuf_put('oc_covarnce_Ni_rg',oc_micro_covarnce(idx_Ni,idx_rg,:),1.)
    call hbuf_put('oc_covarnce_Ni_Ng',oc_micro_covarnce(idx_Ni,idx_Ng,:),1.)
  end if !dograupel

!---------------------------------
! Snow mixing ratio
!--------------------------------

  call hbuf_put('oc_covarnce_rs_Ns',oc_micro_covarnce(idx_rs,idx_Ns,:),1.)
       
  if (dograupel) then   
    call hbuf_put('oc_covarnce_rs_rg',oc_micro_covarnce(idx_rs,idx_rg,:),1.)
    call hbuf_put('oc_covarnce_rs_Ng',oc_micro_covarnce(idx_rs,idx_Ng,:),1.)

!---------------------------------
! Snow number concentration
!--------------------------------
    call hbuf_put('oc_covarnce_Ns_rg',oc_micro_covarnce(idx_Ns,idx_rg,:),1.)
    call hbuf_put('oc_covarnce_Ns_Ng',oc_micro_covarnce(idx_Ns,idx_Ng,:),1.)

!---------------------------------
! Graupel mixing ratio / number concentration
!--------------------------------
    call hbuf_put('oc_covarnce_rg_Ng',oc_micro_covarnce(idx_rg,idx_Ng,:),1.)
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

#ifdef UWM_STATS
!--------------------------------------------------------------------------------------------------
subroutine write_3d_micro_fields()

! Description:
!   This subroutine was created in response to sam_clubb:ticket:82.
!   We would like to look at 3D fields of various microphysical budgets, terms,
!   etc. However, the variables are local to Morrison. Therefore, we write a
!   seperate 3D field for the micro variables.
!
!   This is largely 'write_3D_fields.F90', but reused here. 
!--------------------------------------------------------------------------------------------------

use vars, only: &
    lenstr,  & ! Variable(s)
    t,       &
    qv,      &
    qcl,     &
    w,       &
    qpl,     &
    gamaz,   &
    prespot

use domain, only: &
    nsubdomains_x, & ! Variable(s)
    nsubdomains_y

use grid, only:  &
    masterproc,  & ! Variable(s)
    output_sep,  &
    rank,        &
    nsubdomains, &
    nstep,       &
    RUN3D,       &
    save3Dbin,   &
    case,        &
    caseid,      &
    save3Dsep,   &
    nrestart,    &
    notopened3D, &
    nx,          &
    ny,          &
    nzm,         &
    z,           &
    pres,        &
    dx,          &
    dy,          &
    nstep,       &
    dt,          &
    day0,        &
    dompi,       &
    dogzip3D

use calc_vars_util, only: &
    t2thetal  ! Procedure(s)

use compute_chi_module, only: &
    compute_chi_eta  ! Procedure(s)

implicit none
character *120 filename
character *80 long_name
character *8 name
character *10 timechar
character *4 rankchar
character *5 sepchar
character *12 filetype
character *10 units
character *12 c_z(nzm),c_p(nzm),c_dx, c_dy, c_time
integer i,j,k,n,nfields,nfields1
real tmp(nx,ny,nzm)

real, dimension(nx,ny,nzm) :: &
  thl, & ! Liquid water potential temperature                 [K]
  rt,  & ! Total water mixing ratio                           [kg/kg]
  chi, & ! Extended liquid water mixing ratio                 [kg/kg]
  eta    ! Coordinate orthogonal to chi in PDF transformation [kg/kg]

call t_startf('3D_out')

nfields=37 ! number of 3D fields to save
nfields1=0 ! assertion check

if(masterproc.or.output_sep) then
  
  if(output_sep) then
    write(rankchar,'(i4)') rank
    sepchar="_"//rankchar(5-lenstr(rankchar):4)
  else
     sepchar=""
  end if
  
  write(rankchar,'(i4)') nsubdomains
  write(timechar,'(i10)') nstep
  
  do k=1,11-lenstr(timechar)-1
    timechar(k:k)='0'
  end do

  if(RUN3D) then
    if(save3Dbin) then
      filetype = '_micro.bin3D'
    else
      filetype = '_micro.com3D'
    end if
    
    filename='./OUT_3D/'//trim(case)//'_'//trim(caseid)//'_'// &
        rankchar(5-lenstr(rankchar):4)//'_'//timechar(1:10)//filetype//sepchar
    
    open(46,file=filename,status='unknown',form='unformatted')

  else
    
    if(save3Dbin) then
      if(save3Dsep) then
        filetype = '_micro.bin3D'
      else
        filetype = '_micro.bin2D'
      end if
    else
      if(save3Dsep) then
        filetype = '_micro.com3D'
      else
        filetype = '_micro.com2D'
      end if
    end if
  
    if(save3Dsep) then
      filename='./OUT_3D/'//trim(case)//'_'//trim(caseid)//'_'// &
      rankchar(5-lenstr(rankchar):4)//'_'//timechar(1:10)//filetype//sepchar
      open(46,file=filename,status='unknown',form='unformatted')        
    else
      filename='./OUT_3D/'//trim(case)//'_'//trim(caseid)//'_'// &
      rankchar(5-lenstr(rankchar):4)//filetype//sepchar
      if(nrestart.eq.0.and.notopened3D) then
        open(46,file=filename,status='unknown',form='unformatted')     
      else
        open(46,file=filename,status='unknown', &
                              form='unformatted', position='append')
      end if
    notopened3D=.false.
    end if!save3Dsep
  end if!RUN3D

  if(masterproc) then

    if(save3Dbin) then

      write(46) nx,ny,nzm,nsubdomains,nsubdomains_x,nsubdomains_y,nfields
      do k=1,nzm
        write(46) z(k)
      end do
      do k=1,nzm
        write(46) pres(k)
      end do
      write(46) dx
      write(46) dy
      write(46) nstep*dt/(3600.*24.)+day0

    else
      
      write(long_name,'(8i4)') nx,ny,nzm,nsubdomains, &
                                   nsubdomains_x,nsubdomains_y,nfields
      do k=1,nzm
        write(c_z(k),'(f12.3)') z(k)
      end do
      do k=1,nzm
        write(c_p(k),'(f12.3)') pres(k)
      end do
      write(c_dx,'(f12.0)') dx
      write(c_dy,'(f12.0)') dy
      write(c_time,'(f12.5)') nstep*dt/(3600.*24.)+day0
        
      write(46) long_name(1:32)
      write(46) c_time,c_dx,c_dy, (c_z(k),k=1,nzm),(c_p(k),k=1,nzm)

    end if ! save3Dbin
  end if ! masterproc
end if ! masterproc.or.output_sep


!--------------------------------------
! Micro fields
!--------------------------------------

do i = 1, nx, 1
   do j = 1, ny, 1
      do k = 1, nzm, 1

         ! Calculate rt
         rt(i,j,k) = qv(i,j,k) + qcl(i,j,k)

         ! Calculate thetal
         thl(i,j,k) = t2thetal( t(i,j,k), gamaz(k), qpl(i,j,k), &
                                0.0, 0.0, prespot(k) )

      enddo ! k = 1, nzm, 1
   enddo ! j = 1, ny, 1
enddo ! i = 1, nx, 1

! Calculate the values of chi and eta.
call compute_chi_eta( thl, rt, pres, prespot, &
                      chi, eta )

nfields1=nfields1+1
do k = 1, nzm
   do j = 1, ny
      do i = 1, nx
         tmp(i,j,k) = 0.5 * ( w(i,j,k) + w(i,j,k+1) )
      enddo
   enddo
enddo
name = 'W'
long_name = 'Vertical Velocity'
units = 'm/s'
call compress3D( tmp, nx, ny, nzm, name, long_name, units, &
                 save3Dbin, dompi, rank, nsubdomains )

nfields1=nfields1+1
do k = 1, nzm
   do j = 1, ny
      do i = 1, nx
         tmp(i,j,k) = rt(i,j,k)
      enddo
   enddo
enddo
name = 'RT'
long_name = 'Total water mixing ratio (vapor+cloud)'
units = 'kg/kg'
call compress3D( tmp, nx, ny, nzm, name, long_name, units, &
                 save3Dbin, dompi, rank, nsubdomains )

nfields1=nfields1+1
do k = 1, nzm
   do j = 1, ny
      do i = 1, nx
         tmp(i,j,k) = thl(i,j,k)
      enddo
   enddo
enddo
name = 'THL'
long_name = 'Liquid water potential temperature'
units = 'K'
call compress3D( tmp, nx, ny, nzm, name, long_name,units, &
                 save3Dbin, dompi, rank, nsubdomains )

nfields1=nfields1+1
do k = 1, nzm
   do j = 1, ny
      do i = 1, nx
         tmp(i,j,k) = chi(i,j,k)
      enddo
   enddo
enddo
name = 'CHI'
long_name = 'Extended liquid water mixing ratio, chi'
units = 'kg/kg'
call compress3D( tmp, nx, ny, nzm, name, long_name, units, &
                 save3Dbin, dompi, rank, nsubdomains )

nfields1=nfields1+1
do k = 1, nzm
   do j = 1, ny
      do i = 1, nx
         tmp(i,j,k) = eta(i,j,k)
      enddo
   enddo
enddo
name = 'ETA'
long_name = 'Eta (orthogonal to chi in PDF trans.)'
units = 'kg/kg'
call compress3D( tmp, nx, ny, nzm, name, long_name, units, &
                 save3Dbin, dompi, rank, nsubdomains )

nfields1=nfields1+1
do k = 1, nzm
   do j = 1, ny
      do i = 1, nx
         tmp(i,j,k) = qcl(i,j,k)
      enddo
   enddo
enddo
name = 'RC'
long_name = 'Cloud water mixing ratio'
units = 'kg/kg'
call compress3D( tmp, nx, ny, nzm, name, long_name, units, &
                 save3Dbin, dompi, rank, nsubdomains )



if (doprecip) then

  nfields1=nfields1+1
  do k = 1, nzm
     do j = 1, ny
        do i = 1, nx
           tmp(i,j,k) = qpl(i,j,k)
        enddo
     enddo
  enddo
  name = 'RR'
  long_name = 'Rain water mixing ratio'
  units = 'kg/kg'
  call compress3D( tmp, nx, ny, nzm, name, long_name, units, &
                 save3Dbin, dompi, rank, nsubdomains )

  nfields1=nfields1+1
  do k = 1, nzm
     do j = 1, ny
        do i = 1, nx
           tmp(i,j,k) = micro_field(i,j,k,inr)
        enddo
     enddo
  enddo
  name = 'NR'
  long_name = 'Rain drop concentration'
  units = 'num/kg'
  call compress3D( tmp, nx, ny, nzm, name, long_name, units, &
                 save3Dbin, dompi, rank, nsubdomains )
endif ! doprecip

  !--------------------------------------
  ! RR Budget Terms
  !--------------------------------------
nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mPRE_3D(i,j,k)
    end do
   end do
  end do
  name='PRE'
  long_name='Evaporation of rain'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mPRA_3D(i,j,k)
    end do
   end do
  end do
  name='PRA'
  long_name='Accretion of cloud droplets by rain'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mPRC_3D(i,j,k)
    end do
   end do
  end do
  name='PRC'
  long_name='Autoconversion of cloud droplets to rain'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mPRACS_3D(i,j,k)
    end do
   end do
  end do
  name='PRACS'
  long_name='Collection of rain by snow to form snow'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mMNUCCR_3D(i,j,k)
    end do
   end do
  end do
  name='MNUCCR'
  long_name='Contact freezing of rain drops'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mQMULTR_3D(i,j,k)
    end do
   end do
  end do
  name='QMULTR'
  long_name='Splintering  from rain drops accreted onto snow'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mQMULTRG_3D(i,j,k)
    end do
   end do
  end do
  name='QMULTRG'
  long_name='Splintering from rain drops accreted onto graupel'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mPIACR_3D(i,j,k)
    end do
   end do
  end do
  name='PIACR'
  long_name='Collection of cloud ice by rain to form graupel'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mPIACRS_3D(i,j,k)
    end do
   end do
  end do
  name='PIACRS'
  long_name='Collection of cloud ice by rain to form snow'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mPRACG_3D(i,j,k)
    end do
   end do
  end do
  name='PRACG'
  long_name='Collection of rain by graupel'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mPGRACS_3D(i,j,k)
    end do
   end do
  end do
  name='PGRACS'
  long_name='Collection of rain by snow to form snow'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mPSMLT_3D(i,j,k)
    end do
   end do
  end do
  name='PSMLT'
  long_name='Freezing of rain to form snow'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mPGMLT_3D(i,j,k)
    end do
   end do
  end do
  name='PGMLT'
  long_name='Freezing of rain to form graupel'
  units='kg/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)


!--------------------------------------
! NR Budget Terms
!--------------------------------------
nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mSIZEFIX_NR_3D(i,j,k)
    end do
   end do
  end do
  name='SIZEFIX_NR'
  long_name='Adjust rain number when size is too large/small'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mNEGFIX_NR_3D(i,j,k)
    end do
   end do
  end do
  name='NEGFIX_NR'
  long_name='Removal of negative rain drop number concentration'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mNSUBR_3D(i,j,k)
    end do
   end do
  end do
  name='NSUBR'
  long_name='Evaporation of rain'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mNSMLTR_3D(i,j,k)
    end do
   end do
  end do
  name='NSMLTR'
  long_name='Melting of snow to form rain'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mNGMLTR_3D(i,j,k)
    end do
   end do
  end do
  name='NGMLTR'
  long_name='Melting of graupel to form rain'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mNPRC1_3D(i,j,k)
    end do
   end do
  end do
  name='NPRC1'
  long_name='Change in rain due to autoconversion'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mNPRACS_3D(i,j,k)
    end do
   end do
  end do
  name='NPRACS'
  long_name='Collection of rain by snow'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mNNUCCR_3D(i,j,k)
    end do
   end do
  end do
  name='NNUCCR'
  long_name='Contact freezing of rain'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mNRAGG_3D(i,j,k)
    end do
   end do
  end do
  name='NRAGG'
  long_name='Self collection of rain drops'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mNIACR_3D(i,j,k)
    end do
   end do
  end do
  name='NIACR'
  long_name='Collection of ice by rain to form graupel'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mNIACRS_3D(i,j,k)
    end do
   end do
  end do
  name='NIACRS'
  long_name='Collection of ice by rain to form snow'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mNPRACG_3D(i,j,k)
    end do
   end do
  end do
  name='NPRACG'
  long_name='Collection of rain by graupel'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mNGRACS_3D(i,j,k)
    end do
   end do
  end do
  name='NGRACS'
  long_name='Collection of rain by snow'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=mNR_INST_3D(i,j,k)
    end do
   end do
  end do
  name='NR_INST'
  long_name='Instantaneous processes'
  units='#/kg/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)
nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=rain_vel_3D(i,j,k)
    end do
   end do
  end do
  name='rain_vel'
  long_name='Rain drop velocity'
  units='m/s'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

nfields1=nfields1+1
  do k=1,nzm
   do j=1,ny
    do i=1,nx
      tmp(i,j,k)=EFFR_3D(i,j,k)
    end do
   end do
  end do
  name='EFFR'
  long_name='Effective radius of rain'
  units='micron'
  call compress3D(tmp,nx,ny,nzm,name,long_name,units, &
                                 save3Dbin,dompi,rank,nsubdomains)

call task_barrier()

if(nfields.ne.nfields1) then
  if( masterproc ) print*,'write_fields3D error: nfields=',nfields,'nfields1=',nfields1
  call task_abort()
endif
if ( masterproc ) then
   close (46)
   if(RUN3D .or. save3Dsep) then
     if(dogzip3D) call systemf('gzip -f '//filename)
     print*, 'Writting 3D data. file:'//filename
   else
     print*, 'Appending 3D data. file:'//filename
  end if
endif


call t_stopf('3D_OUT')

end subroutine write_3d_micro_fields
#endif /*UWM_STATS*/


end module microphysics



