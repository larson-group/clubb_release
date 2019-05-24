MODULE CLUBB_driver_mod

!=======================================================================
!   one-dimensional boundary layer cloud pameterization
!=======================================================================
! For a detailed description of the model code see:

!``A PDF-Based Model for Boundary Layer Clouds. Part I:
!  Method and Model Description'' Golaz, et al. (2002)
 ! JAS, Vol. 59, pp. 3540--3551. 
!---------------------------------------------------------------------

! <DATASET NAME="CLUBB.res">
!   native format of the restart file
! </DATASET>
! <DATASET NAME="CLUBB.res.nc">
!   netcdf format of the restart file
! </DATASET>

use            mpp_mod, only: mpp_pe, mpp_root_pe, stdlog, mpp_chksum,  &
                              mpp_clock_id, mpp_clock_begin, mpp_clock_end,  &
                              CLOCK_MODULE_DRIVER
use   Diag_Manager_Mod, ONLY: register_diag_field, send_data
use   Time_Manager_Mod, ONLY: time_type, get_time, set_time, get_date, set_date, operator(+), operator(-)
use            fms_mod, only: write_version_number, open_file,        &
                                open_namelist_file, check_nml_error,  &
                                FILE_EXIST, ERROR_MESG,     &
                                CLOSE_FILE, FATAL,  NOTE,   &
                                read_data, write_data,      &
                                mpp_error, mpp_chksum

use   field_manager_mod,only: MODEL_ATMOS
use  tracer_manager_mod,only: get_number_tracers,           &
                                get_tracer_index,           &
                                get_tracer_names

use  rad_utilities_mod, only: aerosol_type
use  aer_ccn_act_mod,   only: aer_ccn_act_init
use  aer_ccn_act_k_mod, only: aer_ccn_act_k

#ifdef SCM
USE    scm_forc_mod,    only: experiment
#endif

use parameter_indices,  only: nparams              ! Variable(s)

use parameters_tunable, only: params_list          ! Variable(s) 
use parameters_tunable, only: taumax, taumin, c_K  ! Variable(s)
use parameters_tunable, only: read_parameters

use parameters_model,   only: T0, sclr_dim,  cloud_frac_min      

use variables_diagnostic_module, only: setup_diagnostic_variables ! Procedure
use variables_prognostic_module, only: setup_prognostic_variables ! Procedure
      
USE stats_variables, ONLY: l_stats, l_stats_samp        ! Main flag to turn statistics on/off

!cjg   begin
use clubb_precision, only: time_precision
use stats_subs,      only: stats_init, stats_begin_timestep, stats_end_timestep, stats_finalize
!cjg   end

use grid_class, only:  gr ! Variable(s)
use grid_class, only:  setup_grid_heights, zt2zm, zm2zt ! Procedure(s)


use variables_diagnostic_module, only: ug, vg, em,  & ! Variable(s)
      tau_zm, tau_zt, thvm, Lscale, Kh_zm, Kh_zt,   &
      um_ref, vm_ref, Ncnm, wp2_zt,                 &
      hydromet, thlm_ref, rtm_ref, wpthvp

use variables_prognostic_module, only:   & 
      thlm, rtm,     & ! Variable(s)
      um, vm, wp2, rcm, wm_zt, wm_zm, exner,          &
      p_in_Pa, rho_zm, upwp, vpwp, wpthlp,            &
      wprcp, rho, wprtp, wpthlp_sfc, wprtp_sfc,       &
      upwp_sfc, vpwp_sfc, rho_ds_zm, rho_ds_zt,       &
      invrs_rho_ds_zm, invrs_rho_ds_zt, thv_ds_zm,    &
      thv_ds_zt, thlm_forcing, rtm_forcing,           &
      um_forcing, vm_forcing, wprtp_forcing,          &
      wpthlp_forcing, rtp2_forcing, thlp2_forcing,    &
      rtpthlp_forcing, up2, vp2, wp3, rtp2,           &
      thlp2, rtpthlp, pdf_params, cloud_frac,         &
      temp_clubb, rcm_in_layer, cloud_cover

use variables_prognostic_module, only:                   & 
      sclrm, sclrp2, sclrprtp, sclrpthlp, sclrm_forcing, & ! Variables
      wpsclrp, wpsclrp_sfc,                              & 
      edsclrm, edsclrm_forcing, wpedsclrp_sfc

use variables_prognostic_module, only: RH_crit  ! critical relative humidity for droplet and ice nucleation

USE  CLUBB_3D_var,  only:    &
      setup_CLUBB_3D_var,    &
      cleanup_CLUBB_3D_var

USE  CLUBB_3D_var, only:     &
      upwp_3D,               &
      vpwp_3D,               &
      up2_3D,                &
      vp2_3D,                &
      wp2_3D,                &
      wprtp_3D,              &
      wpthlp_3D,             &
      rtp2_3D,               &
      thlp2_3D,              &
      rtpthlp_3D,            &
      wp3_3D,                &          
      RH_crit_clubb_3D

use T_in_K_module, only: thlm2T_in_K ! Procedure

use clubb_core, only:        &
      advance_clubb_core,    &  ! Procedure
      cleanup_clubb_core,    &
      setup_clubb_core

!use numerical_check, only: invalid_model_arrays ! Procedure(s) 
use error_code, only: fatal_error,  set_clubb_debug_level ! Procedure(s)
use constants_mod, only: RAD_TO_DEG

! use constants_mod, only:  grav
use constants_clubb, only:  &
      Cp,                   &
      Lv,                   &
      ep2,                  &
      ep1,                  &
      grav,                 &
      em_min,               &
      w_tol_sqd,            &           
      kappa,                & 
      p0,                   &         
      Rd,                   &
      zero_threshold

use model_flags, only: l_tke_aniso

use alt_cloud_mod, only: alt_cloud

implicit none 
public ::        & 
    CLUBB_SETUP, & 
    CLUBB_INIT,  & 
    CLUBB,       & 
    CLUBB_END
    
public ::             & 
    host2CLUBB_full,  &
    host2CLUBB_half,  &
    CLUBB2host_full,  &
    CLUBB2host_half,  &   
    tndy_CLUBB

public :: CLUBB_3D_2_1D, CLUBB_1D_2_3D

logical              :: module_is_initialized = .false.
character(len=32)    :: tracer_units, tracer_name
character(len=128)   :: diaglname

character(len=9)     :: mod_name = 'CLUBB_BLC'
logical              :: do_netcdf_restart = .true.

real                 :: missing_value = -999.
logical              :: used

!for the calculation of surface momentum fluxes from u_star
real, parameter :: gust_const =  1.0

real, dimension(nparams) :: params  ! Array of the model constants, initialize in "CLUBB_SETUP"

! the vertical levels only are limited within boundary layer!
real ::  CLUBB_height  = 3000.  !  vertical domain height for CLUBB [m]
! the maximum vertical grid spacing for CLUBB
real ::  dz_CLUBB_max = 40.0    !  vertical grid spacing for CLUBB [m]     

integer, save :: nsublevel      ! # of interpolated vertical sub-levels.
                                ! nsublevel=0: no vertical sub-lelvels
integer       :: iz_shift       ! # of shifted vertical levels

real ::  CLUBB_dt  = 60.0       ! time-step used in CLUBB [s]

integer, parameter        :: sclr_max = 10     ! Maximum number of passive scalars (arbitrary)
real, dimension(sclr_max) :: sclr_tol = 1.e-2  ! Thresholds on the passive scalars [units vary]

! minmum cloud fraction
real ::  cloud_frac_min_in = 0.008

logical   ::                       &
   use_sclr_HOC       = .false.,   &
   use_online_aerosol = .true.,    &
   use_sub_seasalt    = .true.,    &
   use_cloud_cover    = .false.,   &
   do_dust_berg       = .false.,   &
   do_drp_evap_cf     = .true.,    &
   do_liq_num         = .true.,    &
   spread_host_tdcy   = .true.

real   :: sea_salt_scale =  0.1
real   :: om_to_oc       =  1.67
real   :: w_min   = 0.001        !  min. W for aerosol activation

logical   ::                        &
   do_BL_gauss         =  .false.,  &
   do_diffK_gauss      =  .false.,  &
   do_quadrature_gauss =  .true.

logical   ::                        &
   l_uv_nudge_in      = .false.

!
logical   ::                                &
   l_host_applies_sfc_fluxes_in =  .true.,  &
   I_sat_sphum_in               =  .true.

real :: host_dx = 4.0e5
real :: host_dy = 4.0e5

integer :: icheck_temp = 0

real :: ts_nudge_in          = 86400.0
real :: T0_in                = 300.0

logical :: do_aeromass_clubb_const = .false.
real    :: aeromass_clubb_const = 2.25e-12

real ::  Var_w = 0.49  ! variance of vertical velocity (w) m2/s2

! for the ice nucleation
logical         ::  do_ice_nucl_wpdf       = .true.
logical         ::  do_ice_nucl_ss_wpdf_in = .false.
real            ::  d_sulf_in =  0.04e-6
real            ::  d_bc_in =  0.0236e-6
logical         ::  do_het_in = .true.
logical         ::  use_dust_instead_of_bc_in = .false.
logical         ::  limit_immersion_frz_in    = .false.
logical         ::  limit_rhil_in             = .false.

real            ::  dust_frac_in = 1.0   
real            ::  rh_crit_het_in = 1.2

real            ::  dust_frac_min_in = 0.0
real            ::  dust_frac_max_in = 1.e15
real            ::  dust_surf_in = 0.5
integer         ::  dust_opt_in = 1
real            ::  rh_dust_max_in = 150.
real            ::  cf_thresh_nucl = 0.98
integer         ::  rh_act_opt = 2   

character(len=6)::  saturation_formula_in = "GFDL  " ! "bolton" approx., "flatau" approx, or "GFDL" approx.

integer ::   debug_level = 2     ! Amount of debugging information

integer, dimension(6)   ::  init_date = (/ 0, 0, 0, 0, 0, 0 /)
integer, save           ::  current_days0,  current_sec0

integer                 ::  do_conv_flag_clubb = 1
real                    ::  avg_deltaz = 50.0

integer                 ::  do_alt_cloud = 0

!--------------------- version number ----------------------------------
!
character(len=128) :: version = '$Id: CLUBB_driver_SCM.F90,v 16.0 2009/02/05 22:09:50 fms Exp $'
character(len=128) :: tagname = '$Name: perth $'
integer, dimension(1) :: restart_versions = (/ 1 /)

! Definition of namelists
namelist /CLUBB_setting_nml/                                &
         use_sclr_HOC,                                      &
         use_online_aerosol,                                &
         use_sub_seasalt,                                   &
         use_cloud_cover,                                   &
         do_dust_berg,                                      &
         do_liq_num,                                        &
         do_drp_evap_cf,                                    &
         cloud_frac_min_in,                                 &
         l_uv_nudge_in,   ts_nudge_in,                      &
         saturation_formula_in,                             &
         w_min,                                             &
         do_quadrature_gauss,  do_diffK_gauss, do_BL_gauss, &
         Var_w,                                             &
         CLUBB_height,                                      &
         CLUBB_dt,                                          &
         dz_clubb_max,                                      &
         debug_level,                                       &
         do_aeromass_clubb_const, aeromass_clubb_const,     &
         do_ice_nucl_wpdf,                                  &
         l_host_applies_sfc_fluxes_in,                      &
         host_dx, host_dy, icheck_temp,                     &
         init_date,                                         &
         do_conv_flag_clubb,                                &
         spread_host_tdcy,                                  &
         avg_deltaz,                                        &
         do_alt_cloud

!cjg begins
! For CLUBB stats
type(time_type) :: Time_stats_init
integer :: unit_stats

logical :: do_stats = .false. ! Whether statistics are computed and output to disk

real, dimension(10) :: lon_stats = -9999.0,  &  ! lon and lat of grid point selected for statistical output
                       lat_stats = -9999.0      ! in degrees.

integer :: istats = -99999,  &     ! index for stats column
           jstats = -99999
integer :: fstats = -99999         ! index for file name

character(len=100), dimension(10) :: & 
      fname_prefix = "clubb"       ! Prefix of stats filenames, to be followed by, for example "_zt"

real(kind=time_precision) :: & 
      stats_tsamp = 0.0,   & ! Stats sampling interval [s]
      stats_tout = 0.0       ! Stats output interval   [s]

character(len=10) :: & 
      stats_fmt = "grads"  ! File format for stats; typically GrADS.

namelist /CLUBB_stats_setting_nml/ & 
      do_stats, lon_stats, lat_stats, fname_prefix, stats_tsamp, stats_tout, stats_fmt
!cjg ends

! local arrays
REAL, dimension(:), allocatable :: &
       momentum_heights, thermodynamic_heights

Integer :: nz_CLUBB,      iz_CLUBB         ! vertical levels for CLUBB
Integer :: nz_host_model, iz_host_model    ! vertical levels for host model

integer :: nsphum, nql, nqi, nqa, nqn, nqni  ! tracer indices for stratiform clouds

! save tendencies due to CLUBB (high order closure)
integer   :: nt                       ! total no. of tracers
integer   :: ntp                      ! total no. of prognostic tracers
integer   :: it                       ! total no. of tracers

! ---> h1g
! dump different aerosol mass concentration
integer ::  id_sulfate,         &
            id_seasalt_sub,     &
            id_seasalt_sup,     &
            id_om
! <--- h1g

! ---> h1g check  ccolumn masss  and energy conservation 
logical, private ::  do_CLUBB_conservation_checks  =  .true.
integer  ::  id_enth_CLUBB_col,    id_wat_CLUBB_col
integer, private , dimension(:), allocatable :: id_tracer_CLUBB_col

real, private, dimension(:, :, :), allocatable :: pmass_3d
real, private, dimension(:, : ),   allocatable :: tmp2D_check
! <---  h1g

integer :: id_udt_CLUBB, id_vdt_CLUBB, id_tdt_CLUBB 
integer, dimension(:), allocatable :: id_tracer_CLUBB

integer :: id_drp_evap_CLUBB
integer :: id_aer_ccn_act_CLUBB
integer :: id_Ndrop_act_CLUBB

integer :: id_drp_flux_CLUBB
integer :: id_qndt_CLUBB_trsport_only

integer :: id_ice_evap_CLUBB
integer :: id_aer_ice_act_CLUBB
integer :: id_Icedrop_act_CLUBB

integer :: id_icedrop_flux_CLUBB
integer :: id_qnidt_CLUBB_trsport_only

integer :: id_wm_CLUBB,   id_omega_CLUBB

!output for higher order moments, fluxes
!at full levels
integer :: id_wp3_CLUBB

!at half levels
integer :: id_wp2_CLUBB,      &
           id_upwp_CLUBB,     &
           id_vpwp_CLUBB,     &
           id_up2_CLUBB,      &
           id_vp2_CLUBB,      &
           id_wprtp_CLUBB,    &
           id_wpthlp_CLUBB,   &
           id_rtp2_CLUBB,     &
           id_thlp2_CLUBB,    &
           id_rtpthlp_CLUBB

integer :: clubb_core_clock

!---------------------------------------------------------------------

 contains


 subroutine clubb(is, ie, js, je, lon, lat,                  &
                  Time_next,                                 &
                  dtmain,                                    &
                  phalf, pfull, zhalf, zfull, omega_avg,     &
                  t, q, r, u, v,                             &
                  u_star, b_star, q_star,                    &
                  tdt, qdt, rdt, udt, vdt,                   &
                  dcond_ls_liquid, dcond_ls_ice,             &
                  Ndrop_act_CLUBB, Icedrop_act_CLUBB,        &
                  diff_t_clubb,                              &
                  Aerosol, mask,                             &
                  mc_full,                                   &
                  conv_frac_clubb,                           &
                  convective_humidity_ratio_clubb)

!-----------------------------------------------------------------------
!
!      is,ie,js,je    starting/ending subdomain i,j indices of data in 
!                     the physics_window being integrated

!      Time_next      next time, used for diagnostics   (time_type)!
!
!         dtmain    time step (from t(n-1) to t(n+1) if leapfrog)
!                       in seconds   [real]
!
!         phalf      pressure at half levels in pascals
!                      [real, dimension(nlon,nlat,nlev+1)]
!
!         pfull      pressure at full levels in pascals
!                      [real, dimension(nlon,nlat,nlev)]
!
!         omega_avg     grid average omega (vertical velocity) at full levels
!                       in pascals per second
!                       [real, dimension(nlon,nlat,nlev)]
!
!         t, q       temperature (t) [deg k] and specific humidity
!                    of water vapor (q) [kg/kg] at full model levels,
!                    at the current time step if leapfrog scheme
!                      [real, dimension(nlon,nlat,nlev)]
!
!         r          tracer fields at full model levels, unit varies,
!                    at the current time step if leapfrog 
!                      [real, dimension(nlon,nlat,nlev,ntrace)]
!
!         u, v,      zonal and meridional wind [m/s] at full model levels,
!                    at the current time step if leapfrog scheme
!                      [real, dimension(nlon,nlat,nlev)]
!  
!        u_star,     friction velocity [m/s]
!        b_star,     buoyancy scale    [m/s^2/K]
!        q_star,     moisture scale    dimensionless

! inout:  tdt, qdt   temperature (tdt) [deg k/sec] and specific
!                    humidity of water vapor (qdt) tendency [1/sec]
!                      [real, dimension(nlon,nlat,nlev)]
!
!         rdt        tracer tendencies , unit varies, 
!                      [real, dimension(nlon,nlat,nlev,ntrace)]
!
!         udt, vdt   zonal and meridional wind tendencies [m/s/s]

!       ---------------
  !       optional INPUT:
  !       ---------------
  !
  !         variable              definition                  unit
  !       ------------   -----------------------------   ---------------
  !
  !       mc_full (kg/m2/s)            (total) net convective mass flux 
  !  adjust environmental vertical velocity due to convection
  !           omega = omega_avg + mc_full*grav

  !       conv_frac_clubb    total convective cloud fraction
  
  !       convective_humidity_ratio_clubb: ratio of 
  !          the grid average specific humidity 
  !          to environmental specific humidity 
  !
  !       Aerosol     different aerosol species mass concentration
  !
  !       mask           real array indicating the 
  !                      point is above the surface
  !                      if equal to 1.0 and 
  !                      indicating the point is below
  !                      the surface if equal to 0.
!--------------------------------------------------------------------------------------------------

integer,            intent(in)          ::  is, ie, js, je
real, intent(in),  dimension(:,:)       ::  lon, lat
type(time_type),    intent(in)          ::  Time_next
type(time_type)                             Time

real, intent(in)                        ::  dtmain
real, intent(in) , dimension(:,:,:)     ::  phalf, pfull, zhalf, zfull, omega_avg
real, intent(in) , dimension(:,:,:)     ::  t, q, u, v
real, intent(in) , dimension(:,:,:,:)   ::  r
real,dimension(:,:),    intent(in)      ::  u_star, b_star, q_star

real, intent(inOUT) , dimension(:,:,:)  :: tdt, qdt, udt, vdt
real, intent(inOUT) , dimension(:,:,:,:):: rdt

! get online aerosol mass
real, dimension(size(t,1),size(t,2),size(t,3))     ::  airdens, concen_dust_sub
real, dimension(size(t,1),size(t,2),size(t,3), 4)  ::  totalmass1    

! large-scale liquid condensation rate
real, intent(OUT) , dimension(:,:,:)            :: dcond_ls_liquid
! large-scale ice condensation rate
real, intent(OUT) , dimension(:,:,:)            :: dcond_ls_ice

!   **********  in current version   **********
! in-cloud (NOT DOMAIN)  averaged activated  droplet number concentration  (#/kg)
real, intent(OUT) , dimension(:,:,:)           :: Ndrop_act_CLUBB
! in-cloud (NOT DOMAIN)  averaged activated  ice-crystal number concentration  (#/kg)
real, intent(OUT) , dimension(:,:,:)           :: Icedrop_act_CLUBB

! eddy diffusivity mixing coefficients from clubb
real, intent(OUT), dimension(:,:,:)            :: diff_t_clubb
! 1d temporary array used to compute time averaged eddy diffusivity
! coefficients. Defined on clubb's grid.
real, dimension( gr%nz )                :: diff_t_1d

real, dimension( is:ie, js:je, size(t,3) )     :: part_liquid
real, dimension( is:ie, js:je, size(t,3) )     :: part_ice
 
type(aerosol_type), intent (in), optional      :: Aerosol
real, intent (in), optional, dimension(:,:,:)  :: mask


real, intent (in), optional, dimension(:,:,:)  :: mc_full
real, intent (in), optional, dimension(:,:,:)  :: conv_frac_clubb
real, intent (in), optional, dimension(:,:,:)  :: convective_humidity_ratio_clubb

! local variables
real, dimension( size(omega_avg,1), size(omega_avg,2), size(omega_avg,3) ) :: omega
real, dimension( size(omega_avg,1), size(omega_avg,2), size(omega_avg,3) ) :: env_qv_scale
real, dimension( size(omega_avg,1), size(omega_avg,2), size(omega_avg,3) ) :: env_condensate_scale
real, dimension( size(rdt,1), size(rdt,2), size(rdt,3), size(rdt,4) )      :: rdt_orig

real, dimension( size(t,3) )                   :: tmp_host
real, dimension( gr%nz )                     :: tmp_CLUBB

! Passive scalar concentration due to pure transport [{units vary}/s]
real, dimension( gr%nz , sclr_dim )          :: sclrm_trsport_only

! updated tracers concentration after CLUBB
real, dimension( gr%nz )                     :: trs_CLUBB

! Pressure on momentum levels                    [Pa]
real, dimension(gr%nz) ::  P_in_Pa_zm
real, dimension(gr%nz) ::  exner_zm
real, dimension(gr%nz) ::  thvm_zm

!aerosol mass concentration
real,  dimension( gr%nz, 4 )               :: aeromass_clubb

real,  dimension( gr%nz )                  :: Ndrop_max
real,  dimension( gr%nz )                  :: rcm_before_clubb
real,  dimension( gr%nz )                  :: cloud_frac_before_clubb
real,  dimension( gr%nz )                  :: drp_before_clubb
real,  dimension( is:ie, js:je, size(t,3) )  :: aer_ccn_act_CLUBB, drp_evap_CLUBB

real,  dimension( is:ie, js:je, size(t,3)+1 ):: drp_flux_CLUBB

! 1D ice crystal number concentration after nucleation
real,  dimension( gr%nz )                  :: Ncrystal_max
real,  dimension( gr%nz )                  :: ice_before_clubb
real,  dimension( is:ie, js:je, size(t,3) )  :: aer_ice_act_CLUBB, ice_evap_CLUBB

real,  dimension( is:ie, js:je, size(t,3)+1 ):: icedrop_flux_CLUBB
 
real    ::  sum_sclrm_trsport
integer ::  isum_sclrm

Integer ::         &
      isub_CLUBB,  &
      nsub_CLUBB,  &
      i_dtmain,    &
      i_CLUBB_dt

integer :: err_code    ! valid run?

Integer ::  ix_host_model, iy_host_model
Integer ::  ix_CLUBB,      iy_CLUBB

integer :: k

! tendencies
real, dimension(is:ie, js:je, size(t,3)) :: &
            uten_CLUBB,           &     ! u-component tendency
            vten_CLUBB,           &     ! u-component tendency
            tten_CLUBB                  ! air temperature tendency

!tracer tendencies
real, dimension(is:ie, js:je, size(r,3), size(r,4)) :: rdt_CLUBB    ! tracer transport tendency

!tracer transport-only tendencies
real, dimension(is:ie, js:je, size(r,3) ) :: qndt_CLUBB_trsport_only      !drop number
real, dimension(is:ie, js:je, size(r,3) ) :: qnidt_CLUBB_trsport_only     !ice number


! higher order terms and fluxes from CLUBB
!at full levels
real, dimension(is:ie, js:je, size(t,3)) :: wp3_CLUBB
real, dimension(is:ie, js:je, size(t,3)) :: wm_CLUBB

!at half levels
real, dimension(is:ie, js:je, size(t,3)+1) :: &
                 wp2_CLUBB,      &
                 upwp_CLUBB,     &
                 vpwp_CLUBB,     &
                 up2_CLUBB,      &
                 vp2_CLUBB,      &
                 wprtp_CLUBB,    &
                 wpthlp_CLUBB,   &
                 rtp2_CLUBB,     &
                 thlp2_CLUBB,    &
                 rtpthlp_CLUBB

! temperature fix to force exact  enthalpy conservation
real    :: tten_CLUBB_fix

!cjg begins
type(time_type)           :: Time_CLUBB
integer                   :: itime_elapsed
real(kind=time_precision) :: time_elapsed
!cjg ends


!-----------------------------------------------------------------------------------------------------------------
! --- Begin Code --- !
! initialization
tmp_host         = 0.0
trs_CLUBB        = 0.0

uten_CLUBB       = 0.0
vten_CLUBB       = 0.0
tten_CLUBB       = 0.0
rdt_CLUBB        = 0.0

wp3_CLUBB        = 0.0
wp2_CLUBB        = 0.0
upwp_CLUBB       = 0.0
vpwp_CLUBB       = 0.0
up2_CLUBB        = 0.0
vp2_CLUBB        = 0.0
wprtp_CLUBB      = 0.0
wpthlp_CLUBB     = 0.0
rtp2_CLUBB       = 0.0
thlp2_CLUBB      = 0.0
rtpthlp_CLUBB    = 0.0

! setup vertical levels for CLUBB
nz_CLUBB = gr%nz
!find # of vertical levels for hostmodel
nz_host_model = size( zhalf, 3 )

i_dtmain    = dtmain
i_CLUBB_dt  = CLUBB_dt

Time = Time_next - set_time( i_dtmain )

if( mod(i_dtmain, i_CLUBB_dt) /= 0 ) then
    print*,  "ERROR: the time step in host atmosphere must be an intergal multiple of the time step in CLUBB"
    STOP
endif

nsub_CLUBB = i_dtmain/i_CLUBB_dt
if(  nsub_CLUBB < 1 ) then
    print*,  "ERROR: the time step in host atmosphere must be larger than or equal to the time step in CLUBB"
    STOP
endif

if( nqn > 0 ) then 
       aeromass_clubb           = 0.0
       Ndrop_max                = 0.0
       rcm_before_clubb         = 0.0
       cloud_frac_before_clubb  = 0.0
       drp_before_clubb         = 0.0
       aer_ccn_act_CLUBB        = 0.0
       drp_evap_CLUBB           = 0.0
       Ndrop_act_CLUBB          = 0.0
       qndt_clubb_trsport_only  = 0.0
       sclrm_trsport_only(:, 1) = 0.0
endif

if( nqni > 0 ) then
      if(do_ice_nucl_wpdf == .false.)  &
        call error_mesg ('CLUBB_driver_mod', 'nqni > 0, but do_ice_nucl_wpdf is false', FATAL)
       Icedrop_act_CLUBB          = 0.0
       Ncrystal_max               = 0.0
       ice_before_clubb          = 0.0
       ice_evap_CLUBB            = 0.0
       qnidt_clubb_trsport_only   = 0.0
       sclrm_trsport_only (:, 2)  = 0.0
endif
 
! initialize large-scale condensation (liquid + ice), which is the difference of
dcond_ls_liquid  = 0.0
dcond_ls_ice     = 0.0

concen_dust_sub = 0.0
totalmass1 = 0.0

airdens = pfull/(Rd*t*(1.0-r(:,:,:,nql)-r(:,:,:,nqi)))
! get aerosol mass
if( nqn> 0 )  call get_aer_mass_host (is, js,  Time_next, phalf, airdens, t, &
                       concen_dust_sub, totalmass1, Aerosol, mask)

omega                = omega_avg
env_qv_scale         = 1.0
env_condensate_scale = 1.0
rdt_orig             = rdt

if ( do_conv_flag_clubb > 0 ) then
   if ( present(mc_full) ) then
      omega = omega + mc_full*grav
   endif
   if ( present( conv_frac_clubb ) ) then
      omega = omega / ( 1.0 -  conv_frac_clubb )
   endif


   if( do_conv_flag_clubb >1 ) then
      if ( present( convective_humidity_ratio_clubb ) ) then
          env_qv_scale =  convective_humidity_ratio_clubb
          where ( convective_humidity_ratio_clubb .lt. 0. )
            env_qv_scale = 1.0
          end where
      endif

      if ( present( conv_frac_clubb ) ) then
          env_condensate_scale = 1.0 -  conv_frac_clubb
      endif
   endif
endif


do iy_host_model =  js, je
 do ix_host_model = is, ie
 
!re-adjust model height (from host model)
   !zhalf: 0, dz, 3dz, 3dz, ...
     call host2CLUBB_half(zhalf(ix_host_model, iy_host_model, :),    &    !intent (in)
                                momentum_heights)                         !intent (out)
 
 !zfull: 0.5dz, 1.5dz, 2.5dz, 3.5dz, ...
     call host2CLUBB_full (zfull( ix_host_model, iy_host_model, :), &     !intent (in)
                           thermodynamic_heights)                         !intent (out)   

     call setup_grid_heights(.true., 2, avg_deltaz, momentum_heights(1), momentum_heights, thermodynamic_heights)
! end of re-adjust model height (from host model)

 ! load higher order terms and high-res results
     ix_CLUBB=ix_host_model - is +1
     iy_CLUBB=iy_host_model - js +1

     call CLUBB_3D_2_1D( ix_CLUBB, iy_CLUBB)
! get pressure on thermodynamic points[Pa]
     call host2CLUBB_full(pfull( ix_host_model, iy_host_model, :), &    !intent (in)
                          p_in_Pa)                                      !intent (out)
     p_in_Pa( 1 ) = phalf(ix_host_model, iy_host_model,  nz_host_model) ! surface pressure  [Pa]
! get exner  on thermodynamic points
     exner = ( p_in_Pa/p0 )**kappa
     exner(1) = exner(2)
 
! get pressure on momentum points[Pa]
     call host2CLUBB_half(phalf( ix_host_model, iy_host_model, :),     &    !intent (in)
                          p_in_Pa_zm)                                       !intent (out)
! get exner  on momentum points
     exner_zm = ( p_in_Pa_zm/p0 )**kappa

! calculate tempearture on thermodynamic levels
     call host2CLUBB_full(t( ix_host_model, iy_host_model, :), &    !intent (in)
                          temp_clubb(:) )                           !intent (out) 
     temp_clubb(1) = temp_clubb(2)
          
! calculate liquid potential tempearture on thermodynamic levels
     tmp_host(:)= t(ix_host_model, iy_host_model, :) - Lv/Cp &
                  *(r(ix_host_model, iy_host_model, :, nql) + r(ix_host_model, iy_host_model, :, nqi)) &
                   /env_condensate_scale(ix_host_model, iy_host_model, :)
     call host2CLUBB_full (tmp_host, &     !intent (in)
                           thlm )          !intent (out)
     thlm      = thlm/exner
     thlm( 1 ) = thlm( 2 )
 
! calculate total water content on thermodynamic levels [kg/kg]
     tmp_host(:) = r(ix_host_model, iy_host_model, : , nsphum)/env_qv_scale(ix_host_model, iy_host_model, :)       &
                  +r(ix_host_model, iy_host_model, : , nql)/env_condensate_scale(ix_host_model, iy_host_model, :)  &
                  +r(ix_host_model, iy_host_model, : , nqi)/env_condensate_scale(ix_host_model, iy_host_model, :)

     call host2CLUBB_full(tmp_host, &         ! intent (in)
                          rtm )               ! intent (out)
     rtm( 1 ) = rtm( 2 )

! calculate liquid water content on thermodynamic levels [kg/kg], 
! add liquid and ice together as liquid as inputs for CLUBB
     tmp_host(:) =  r(ix_host_model, iy_host_model, :, nql)/env_condensate_scale(ix_host_model, iy_host_model, :) &
                   +r(ix_host_model, iy_host_model, :, nqi)/env_condensate_scale(ix_host_model, iy_host_model, :)
     call host2CLUBB_full(tmp_host, &         ! intent (in)
                          rcm )               ! intent (out)        
     rcm( 1 ) = rcm( 2 )
 
 ! calculate cloud fraction thermodynamic levels
     tmp_host(:)= r(ix_host_model, iy_host_model, :, nqa)
     call host2CLUBB_full(tmp_host,        &         ! intent (in)
                          cloud_frac )               ! intent (out)
     cloud_frac( 1 ) = cloud_frac( 2 )
 
     thvm = thlm + ep1 * T0 * rtm + ( Lv/(Cp*exner) - ep2*T0 ) * rcm
     do  iz_clubb = 1, nz_clubb-1
        thvm_zm( iz_clubb ) = zt2zm(thlm, iz_clubb)  + ep1 * T0 *  zt2zm(rtm, iz_clubb)  &
                             +( Lv/(Cp*exner_zm(iz_clubb) ) - ep2*T0) * zt2zm(rcm, iz_clubb)
     enddo
     thvm_zm( nz_clubb ) = 2.0*thvm_zm(nz_clubb-1) - thvm_zm(nz_clubb-2)

! ---> h1g, 2010-10-01
     thv_ds_zt = temp_clubb/exner * ( 1.0+ep2*(rtm-rcm) )**kappa
     do iz_clubb = 1, nz_clubb-1
       thv_ds_zm(iz_clubb) = zt2zm(temp_clubb/exner, iz_clubb) &
                       * ( 1.0 + ep2 * max( zt2zm( rtm - rcm, iz_clubb), &
                                            zero_threshold ) )**kappa
     enddo
     iz_clubb = nz_clubb
     thv_ds_zm(iz_clubb) = 2.0*thv_ds_zm(iz_clubb-1)-thv_ds_zm(iz_clubb-2)
! <--- h1g, 2010-10-01

! calculate droplet number concentration
      if( nqn > 0 ) then
! domain average droplet number concentration
        tmp_host(:) = r(ix_host_model, iy_host_model, :, nqn) /env_condensate_scale(ix_host_model, iy_host_model, :)
        if(use_sclr_HOC) then
            call host2CLUBB_full(tmp_host,           &       ! intent (in)
                                 sclrm(:, 1))                ! intent (out)
            do  iz_clubb = 1, nz_clubb
                if( cloud_frac(iz_clubb) <= cloud_frac_min ) sclrm( iz_clubb, 1) = 0.0
            enddo
          Else
            call host2CLUBB_full(tmp_host,           &       ! intent (in)
                                 edsclrm(:, 1))              ! intent (out)
            do  iz_clubb = 1, nz_clubb
                if( cloud_frac(iz_clubb) <= cloud_frac_min ) edsclrm( iz_clubb, 1) = 0.0
            enddo
        endif
      endif


 ! calculate ice number concentration
      if( nqni > 0 ) then
! domain average ice number concentration
         tmp_host(:)= r( ix_host_model, iy_host_model, :, nqni) /env_condensate_scale(ix_host_model, iy_host_model, :) 
         if( use_sclr_HOC ) then
              call host2CLUBB_full(tmp_host,        &       ! intent (in)
                                   sclrm( :, 2))            ! intent (out)
              do  iz_clubb = 1, nz_clubb
                 if( cloud_frac(iz_clubb) <= cloud_frac_min) sclrm( iz_clubb, 2) = 0.0
              enddo
         Else
              call host2CLUBB_full(tmp_host,         &      ! intent (in)
                                   edsclrm( : , 2))         ! intent (out)      
              do  iz_clubb = 1, nz_clubb
                 if( cloud_frac(iz_clubb) <= cloud_frac_min) edsclrm( iz_clubb, 2) = 0.0
              enddo
         endif !  use_sclr_HOC
      endif  ! nqni > 0

! ---> h1g, 2010-06-11
! in order to have better conservation properties, air density is calculated 
!    using hydro-static approximation, rather than gas state equation
! at thermo-dynamic levels
     do iz_clubb = 2, nz_clubb
         rho(iz_clubb) = - (p_in_Pa_zm (iz_clubb) - p_in_Pa_zm (iz_clubb-1)) &
                          * gr%invrs_dzt(iz_clubb)/grav
     enddo
     rho(1) =  rho(2)! surface air density  [kg/m3]

! at momentum levels
     do iz_clubb = 2, nz_clubb-1
         rho_zm( iz_clubb ) = - (p_in_Pa(iz_clubb+1) - p_in_Pa(iz_clubb)) &
                              * gr%invrs_dzm(iz_clubb)/grav
     enddo
     rho_zm(1)        = rho(1)
     rho_zm(nz_clubb) = rho(nz_clubb)
! <--- h1g, 2010-06-11

     rho_ds_zt  =  rho
     rho_ds_zm  =  rho_zm

     invrs_rho_ds_zt  = 1.0 / rho_ds_zt
     invrs_rho_ds_zm  = 1.0 / rho_ds_zm

! horizontal velocity on thermodynamic points
     call host2CLUBB_full( u(ix_host_model, iy_host_model, :) , &   ! intent (in)
                           um)                                      ! intent (out)
     um(1) = um(2)

     call host2CLUBB_full( v(ix_host_model, iy_host_model, :) , &   ! intent (in)
                           vm )                                     ! intent (out)
     vm(1) = vm(2)

! vertical velocity on thermodynamic points
     call host2CLUBB_full( omega(ix_host_model, iy_host_model, :), &   !intent (in)
                           wm_zt )                                     !intent (out)
! Boundary conditions on subsidence (thermodynamic grid)
     wm_zt    = -wm_zt/grav/rho
     wm_zt(1) = 0.0        ! Below surface

! vertical velocity on momentum points
     wm_zm = zt2zm(wm_zt)

! Boundary conditions on subsidence (mom. grid)
     wm_zm(1)        = 0.0               ! At surface
     wm_zm(nz_CLUBB) = 0.0               ! Model top

   ! Apply host tendencies at once
   if (.not.spread_host_tdcy) then
     call add_host_tdcy(                                                    &
            ix_host_model, iy_host_model, dtmain,                           & ! in
            udt, vdt, tdt, rdt, env_qv_scale, env_condensate_scale,         & ! in
            u_star, b_star, q_star,                                         & ! in
            um, vm, thlm, thvm, rtm, cloud_frac, rcm, edsclrm, sclrm,       & ! inout
            upwp_sfc, vpwp_sfc, wpthlp_sfc, wprtp_sfc )                       ! out
   end if

   err_code = 0
   do isub_CLUBB=1, nsub_CLUBB  

     Time_CLUBB = Time + set_time( int( isub_CLUBB*CLUBB_dt ) )

     ! Apply host tendencies by spreading them over the time step
     if (spread_host_tdcy) then
       call add_host_tdcy(                                                    &
              ix_host_model, iy_host_model, CLUBB_dt,                         & ! in
              udt, vdt, tdt, rdt, env_qv_scale, env_condensate_scale,         & ! in
              u_star, b_star, q_star,                                         & ! in
              um, vm, thlm, thvm, rtm, cloud_frac, rcm, edsclrm, sclrm,       & ! inout
              upwp_sfc, vpwp_sfc, wpthlp_sfc, wprtp_sfc )                       ! out
     end if

     if ( nqn > 0 ) then
 ! droplet evaporation
        rcm_before_clubb        = rcm
        cloud_frac_before_clubb = cloud_frac
        if( use_sclr_HOC ) then
            drp_before_clubb( : ) = max( 0.0, sclrm( : , 1) )
        else
            drp_before_clubb( : ) = max( 0.0, edsclrm( : , 1) )
        endif
     endif !nqn > 0

     if ( nqni > 0 ) then
 ! ice particle evaporation
        if( use_sclr_HOC ) then
            ice_before_clubb( : ) = max( 0.0, sclrm( : , 2) )
        else
            ice_before_clubb( : ) = max( 0.0, edsclrm( : , 2) )
        endif
     endif !nqni > 0

     ! Activate CLUBB internal stats
     if ( do_stats .and. ix_host_model.eq.istats .and. iy_host_model.eq.jstats ) then
       call get_time( Time_CLUBB - Time_stats_init, itime_elapsed )
       time_elapsed = itime_elapsed
       l_stats = .true.
       call stats_begin_timestep( time_elapsed )
     else
       l_stats = .false.
       l_stats_samp = .false.
     end if

     thlm_forcing    = 0.0
     rtm_forcing     = 0.0
     um_forcing      = 0.0
     vm_forcing      = 0.0
     wprtp_forcing   = 0.0
     wpthlp_forcing  = 0.0
     rtp2_forcing    = 0.0
     thlp2_forcing   = 0.0
     rtpthlp_forcing = 0.0

     wpsclrp_sfc(1: sclr_dim)         = 0.0
     wpedsclrp_sfc(1: sclr_dim)       = 0.0
     sclrm_forcing(: , 1: sclr_dim)   = 0.0
     edsclrm_forcing(: , 1: sclr_dim) = 0.0

     call mpp_clock_begin ( clubb_core_clock )

     call advance_clubb_core( .true. , CLUBB_dt, 0.0,  momentum_heights(1), &          ! Intent(in)
               thlm_forcing, rtm_forcing, um_forcing, vm_forcing, &    ! Intent(in)
               sclrm_forcing, edsclrm_forcing, wprtp_forcing, &        ! Intent(in)
               wpthlp_forcing, rtp2_forcing, thlp2_forcing, &          ! Intent(in)
               rtpthlp_forcing, wm_zm, wm_zt, &                        ! Intent(in)
               wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, &            ! Intent(in)
               wpsclrp_sfc, wpedsclrp_sfc, &                           ! Intent(in)
               p_in_Pa, rho_zm, rho, exner, &                          ! Intent(in)
               rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &                ! Intent(in)
               invrs_rho_ds_zt, thv_ds_zm, thv_ds_zt, &                ! Intent(in)
               um, vm, upwp, vpwp, up2, vp2, &                         ! Intent(inout)
               thlm, rtm, wprtp, wpthlp,  &                            ! Intent(inout)
               wp2, wp3, rtp2, thlp2, rtpthlp, &                       ! Intent(inout)
               sclrm,   &                                              ! Intent(inout)
#ifdef GFDL
               sclrm_trsport_only,  &  ! h1g, 2010-06-16               ! Intent(inout)
#endif
               sclrp2, sclrprtp, sclrpthlp, &                          ! Intent(inout)
               wpsclrp, edsclrm, err_code, &                           ! Intent(inout)
#ifdef GFDL
               RH_crit, &  ! h1g, 2010-06-16                           ! Intent(out)
#endif
               rcm, wprcp, cloud_frac, &                               ! Intent(out)
               rcm_in_layer, cloud_cover, &                            ! Intent(out)
               pdf_params )                                            ! Intent(out)

     call mpp_clock_end ( clubb_core_clock )

!    Option to use cloud_cover and rcm_in_layer instead of cloud_frac and rcm
!    suggested by Vince
     if (use_cloud_cover) then
       cloud_frac = cloud_cover
       rcm = rcm_in_layer
     endif

     thlm( 1 ) =  thlm( 2 )

     if ( fatal_error( err_code ) ) exit

     ! End of timestep for CLUBB internal stats
     if ( do_stats .and. ix_host_model.eq.istats .and. iy_host_model.eq.jstats ) then
       call stats_end_timestep()
     end if

!convert liquid potential temperature to air temperature, and calculate temperature tendency
     temp_clubb = thlm2T_in_K( thlm, exner, rcm)

! for droplet number evaporation due to macro-physics
     if(nqn>0)  then 
      iz_shift = nsublevel/2

      if(use_sclr_HOC) then
       if(do_drp_evap_cf) then 
          do iz_CLUBB=2+iz_shift, nz_CLUBB, nsublevel+1
             iz_host_model = nz_host_model - (iz_CLUBB -1+iz_shift)/(nsublevel+1)  
             if( ( cloud_frac_before_clubb(iz_CLUBB) - cloud_frac(iz_CLUBB) ) > 0.0 .and. cloud_frac(iz_CLUBB)>0.0) then
! droplet evaporation
                 drp_evap_CLUBB( ix_host_model, iy_host_model,  iz_host_model) = &
                    drp_evap_CLUBB( ix_host_model, iy_host_model,  iz_host_model)   &
                   - min( sclrm( iz_CLUBB, 1 ), drp_before_clubb( iz_CLUBB )    &
                      * ( cloud_frac_before_clubb(iz_CLUBB) -cloud_frac(iz_CLUBB) ) &
                      /cloud_frac_before_clubb( iz_CLUBB ) )

                 sclrm( iz_CLUBB, 1 ) = sclrm( iz_CLUBB, 1 ) &
                                    - min( sclrm( iz_CLUBB, 1 ), drp_before_clubb( iz_CLUBB )    &
                                       * ( cloud_frac_before_clubb(iz_CLUBB) -cloud_frac(iz_CLUBB) ) &
                                      /cloud_frac_before_clubb( iz_CLUBB ) )

             endif

             if( rcm(iz_CLUBB) <= 0.0 .or. cloud_frac(iz_CLUBB) <=cloud_frac_min ) then
! in clear areas, droplet concentration should be 0
                 drp_evap_CLUBB( ix_host_model, iy_host_model,  iz_host_model) = &
                   drp_evap_CLUBB( ix_host_model, iy_host_model,  iz_host_model)   &
                  - sclrm( iz_CLUBB, 1 )
 
                 sclrm( iz_CLUBB, 1 ) = 0.0

             endif
          enddo !enddo iz_CLUBB
       else !NOT do_drp_evap_cf
          do iz_CLUBB=2+iz_shift, nz_CLUBB, nsublevel+1
             iz_host_model = nz_host_model - (iz_CLUBB -1+iz_shift)/(nsublevel+1)  
             if( ( rcm_before_clubb(iz_CLUBB) - rcm(iz_CLUBB) ) > 0.0 .and. rcm(iz_CLUBB)>0.0) then
! droplet evaporation
               drp_evap_CLUBB( ix_host_model, iy_host_model,  iz_host_model) = &
                 drp_evap_CLUBB( ix_host_model, iy_host_model,  iz_host_model)   &
                 - min( sclrm( iz_CLUBB, 1 ), drp_before_clubb( iz_CLUBB )    &
                   * ( rcm_before_clubb(iz_CLUBB) -rcm(iz_CLUBB) ) &
                   /rcm_before_clubb( iz_CLUBB ) )

               sclrm( iz_CLUBB, 1 ) = sclrm( iz_CLUBB, 1 ) &
                                    - min( sclrm( iz_CLUBB, 1 ), drp_before_clubb( iz_CLUBB )    &
                                      * ( rcm_before_clubb(iz_CLUBB) -rcm(iz_CLUBB) ) &
                                      /rcm_before_clubb( iz_CLUBB ) )

             endif

             if( rcm(iz_CLUBB) <= 0.0 .or. cloud_frac(iz_CLUBB) <=cloud_frac_min ) then
! in clear areas, droplet concentration should be 0
                drp_evap_CLUBB( ix_host_model, iy_host_model,  iz_host_model) = &
                  drp_evap_CLUBB( ix_host_model, iy_host_model,  iz_host_model)   &
                 - sclrm( iz_CLUBB, 1 )
 
                sclrm( iz_CLUBB, 1 ) = 0.0

             endif
          enddo !enddo iz_CLUBB
       endif !endif do_drp_evap_cf

     else ! not use scalar HOC

       if(do_drp_evap_cf) then 
          do iz_CLUBB=2+iz_shift, nz_CLUBB, nsublevel+1
            iz_host_model = nz_host_model - (iz_CLUBB -1+iz_shift)/(nsublevel+1)            
            if( ( cloud_frac_before_clubb(iz_CLUBB) - cloud_frac(iz_CLUBB) ) > 0.0 .and. cloud_frac(iz_CLUBB)>0.0) then
! droplet evaporation
                 drp_evap_CLUBB( ix_host_model, iy_host_model,  iz_host_model) = &
                    drp_evap_CLUBB( ix_host_model, iy_host_model,  iz_host_model)   &
                   - min( edsclrm( iz_CLUBB, 1 ), drp_before_clubb( iz_CLUBB )    &
                      * ( cloud_frac_before_clubb(iz_CLUBB) -cloud_frac(iz_CLUBB) ) &
                         /cloud_frac_before_clubb( iz_CLUBB ) )

                  edsclrm( iz_CLUBB, 1 ) = edsclrm( iz_CLUBB, 1 ) &
                                            - min( edsclrm( iz_CLUBB, 1 ), drp_before_clubb( iz_CLUBB )    &
                                              * ( cloud_frac_before_clubb(iz_CLUBB) -cloud_frac(iz_CLUBB) ) &
                                             /cloud_frac_before_clubb( iz_CLUBB ) )

            endif

            if( rcm(iz_CLUBB) <= 0.0 .or. cloud_frac(iz_CLUBB) <=cloud_frac_min ) then
! in clear areas, droplet concentration should be 0
               drp_evap_CLUBB( ix_host_model, iy_host_model,  iz_host_model) = &
                   drp_evap_CLUBB( ix_host_model, iy_host_model,  iz_host_model)   &
                  - edsclrm( iz_CLUBB, 1 )
 
               edsclrm( iz_CLUBB, 1 ) = 0.0

             endif
          enddo !enddo iz_CLUBB
        else !NOT do_drp_evap_cf
          do iz_CLUBB=2+iz_shift, nz_CLUBB, nsublevel+1
             iz_host_model = nz_host_model - (iz_CLUBB -1+iz_shift)/(nsublevel+1)
             if( ( rcm_before_clubb(iz_CLUBB) - rcm(iz_CLUBB) ) > 0.0 .and. rcm(iz_CLUBB)>0.0) then
! droplet evaporation
                  drp_evap_CLUBB( ix_host_model, iy_host_model,  iz_host_model) = &
                    drp_evap_CLUBB( ix_host_model, iy_host_model,  iz_host_model)   &
                   - min( edsclrm( iz_CLUBB, 1 ), drp_before_clubb( iz_CLUBB )    &
                      * ( rcm_before_clubb(iz_CLUBB) -rcm(iz_CLUBB) ) &
                      /rcm_before_clubb( iz_CLUBB ) )

                   edsclrm( iz_CLUBB, 1 ) = edsclrm( iz_CLUBB, 1 ) &
                                                 - min( edsclrm( iz_CLUBB, 1 ), drp_before_clubb( iz_CLUBB )    &
                                                 * ( rcm_before_clubb(iz_CLUBB) -rcm(iz_CLUBB) ) &
                                                 /rcm_before_clubb( iz_CLUBB ) )

             endif

             if( rcm(iz_CLUBB) <= 0.0 .or. cloud_frac(iz_CLUBB) <=cloud_frac_min ) then
! in clear areas, droplet concentration should be 0
                 drp_evap_CLUBB( ix_host_model, iy_host_model,  iz_host_model) = &
                   drp_evap_CLUBB( ix_host_model, iy_host_model,  iz_host_model)   &
                  - edsclrm( iz_CLUBB, 1 )
 
                  edsclrm( iz_CLUBB, 1 ) = 0.0

             endif
          enddo
        endif !endif do_drp_evap_cf
      endif !use_sclr_HOC
     endif !nqn > 0


! for ice crystal evaporation due to macro-physics
     if(nqni>0)  then 
      iz_shift = nsublevel/2

      if(use_sclr_HOC) then
       if(do_drp_evap_cf) then 
          do iz_CLUBB=2+iz_shift, nz_CLUBB, nsublevel+1
             iz_host_model = nz_host_model - (iz_CLUBB -1+iz_shift)/(nsublevel+1)  
             if( ( cloud_frac_before_clubb(iz_CLUBB) - cloud_frac(iz_CLUBB) ) > 0.0 .and. cloud_frac(iz_CLUBB)>0.0) then
! ice crystal evaporation
                 ice_evap_CLUBB( ix_host_model, iy_host_model,  iz_host_model) = &
                    ice_evap_CLUBB( ix_host_model, iy_host_model,  iz_host_model)   &
                   - min( sclrm( iz_CLUBB, 2 ), ice_before_clubb( iz_CLUBB )    &
                      * ( cloud_frac_before_clubb(iz_CLUBB) -cloud_frac(iz_CLUBB) ) &
                      /cloud_frac_before_clubb( iz_CLUBB ) )

                 sclrm( iz_CLUBB, 2 ) = sclrm( iz_CLUBB, 2 ) &
                                    - min( sclrm( iz_CLUBB, 2 ), ice_before_clubb( iz_CLUBB )    &
                                       * ( cloud_frac_before_clubb(iz_CLUBB) -cloud_frac(iz_CLUBB) ) &
                                      /cloud_frac_before_clubb( iz_CLUBB ) )

             endif

             if( rcm(iz_CLUBB) <= 0.0 .or. cloud_frac(iz_CLUBB) <=cloud_frac_min ) then
! in clear areas, ice crystal concentration should be 0
                 ice_evap_CLUBB( ix_host_model, iy_host_model,  iz_host_model) = &
                   ice_evap_CLUBB( ix_host_model, iy_host_model,  iz_host_model)   &
                  - sclrm( iz_CLUBB, 2 )

                 sclrm( iz_CLUBB, 2 ) = 0.0

             endif
          enddo !enddo iz_CLUBB
       else !NOT do_drp_evap_cf
          do iz_CLUBB=2+iz_shift, nz_CLUBB, nsublevel+1
             iz_host_model = nz_host_model - (iz_CLUBB -1+iz_shift)/(nsublevel+1)  
             if( ( rcm_before_clubb(iz_CLUBB) - rcm(iz_CLUBB) ) > 0.0 .and. rcm(iz_CLUBB)>0.0) then
! ice crystal evaporation
                 ice_evap_CLUBB( ix_host_model, iy_host_model,  iz_host_model) = &
                    ice_evap_CLUBB( ix_host_model, iy_host_model,  iz_host_model)   &
                 - min( sclrm( iz_CLUBB, 2 ), ice_before_clubb( iz_CLUBB )    &
                  * ( rcm_before_clubb(iz_CLUBB) -rcm(iz_CLUBB) ) &
                   /rcm_before_clubb( iz_CLUBB ) )

                 sclrm( iz_CLUBB, 2 ) = sclrm( iz_CLUBB, 2 ) &
                                    - min( sclrm( iz_CLUBB, 2 ), ice_before_clubb( iz_CLUBB )    &
                                       * ( rcm_before_clubb(iz_CLUBB) -rcm(iz_CLUBB) ) &
                                        /rcm_before_clubb( iz_CLUBB ) )
             endif

             if( rcm(iz_CLUBB) <= 0.0 .or. cloud_frac(iz_CLUBB) <=cloud_frac_min ) then
! in clear areas, ice crystal concentration should be 0
                 ice_evap_CLUBB( ix_host_model, iy_host_model,  iz_host_model) = &
                   ice_evap_CLUBB( ix_host_model, iy_host_model,  iz_host_model)   &
                  - sclrm( iz_CLUBB, 2 )

                 sclrm( iz_CLUBB, 2 ) = 0.0
             endif
          enddo !enddo iz_CLUBB
       endif !endif do_drp_evap_cf

     else ! not use scalar HOC

       if(do_drp_evap_cf) then 
          do iz_CLUBB=2+iz_shift, nz_CLUBB, nsublevel+1
            iz_host_model = nz_host_model - (iz_CLUBB -1+iz_shift)/(nsublevel+1)            
            if( ( cloud_frac_before_clubb(iz_CLUBB) - cloud_frac(iz_CLUBB) ) > 0.0 .and. cloud_frac(iz_CLUBB)>0.0) then
! ice crystal evaporation
                 ice_evap_CLUBB( ix_host_model, iy_host_model,  iz_host_model) = &
                    ice_evap_CLUBB( ix_host_model, iy_host_model,  iz_host_model)   &
                   - min( edsclrm( iz_CLUBB, 2 ), ice_before_clubb( iz_CLUBB )    &
                      * ( cloud_frac_before_clubb(iz_CLUBB) -cloud_frac(iz_CLUBB) ) &
                      /cloud_frac_before_clubb( iz_CLUBB ) )

                 edsclrm( iz_CLUBB, 2 ) = edsclrm( iz_CLUBB, 2 ) &
                                    - min( edsclrm( iz_CLUBB, 2 ), ice_before_clubb( iz_CLUBB )    &
                                       * ( cloud_frac_before_clubb(iz_CLUBB) -cloud_frac(iz_CLUBB) ) &
                                      /cloud_frac_before_clubb( iz_CLUBB ) )
            endif

            if( rcm(iz_CLUBB) <= 0.0 .or. cloud_frac(iz_CLUBB) <=cloud_frac_min ) then
! in clear areas, ice crystal concentration should be 0
               ice_evap_CLUBB( ix_host_model, iy_host_model,  iz_host_model) = &
                   ice_evap_CLUBB( ix_host_model, iy_host_model,  iz_host_model)   &
                  - edsclrm( iz_CLUBB, 2 )

               edsclrm( iz_CLUBB, 2 ) = 0.0
             endif
          enddo !enddo iz_CLUBB
        else !NOT do_drp_evap_cf
          do iz_CLUBB=2+iz_shift, nz_CLUBB, nsublevel+1
             iz_host_model = nz_host_model - (iz_CLUBB -1+iz_shift)/(nsublevel+1)
             if( ( rcm_before_clubb(iz_CLUBB) - rcm(iz_CLUBB) ) > 0.0 .and. rcm(iz_CLUBB)>0.0) then
! ice crystal evaporation
                 ice_evap_CLUBB( ix_host_model, iy_host_model,  iz_host_model) = &
                    ice_evap_CLUBB( ix_host_model, iy_host_model,  iz_host_model)   &
                                   - min( edsclrm( iz_CLUBB, 2 ), ice_before_clubb( iz_CLUBB )    &
                                      * ( rcm_before_clubb(iz_CLUBB) -rcm(iz_CLUBB) ) &
                                      /rcm_before_clubb( iz_CLUBB ) )

                 edsclrm( iz_CLUBB, 2 ) = edsclrm( iz_CLUBB, 2 ) &
                                   - min( edsclrm( iz_CLUBB, 2 ), ice_before_clubb( iz_CLUBB )    &
                                      * ( rcm_before_clubb(iz_CLUBB) -rcm(iz_CLUBB) ) &
                                      /rcm_before_clubb( iz_CLUBB ) )

             endif

             if( rcm(iz_CLUBB) <= 0.0 .or. cloud_frac(iz_CLUBB) <=cloud_frac_min ) then
! in clear areas, ice crystal concentration should be 0
                 ice_evap_CLUBB( ix_host_model, iy_host_model,  iz_host_model) = &
                   ice_evap_CLUBB( ix_host_model, iy_host_model,  iz_host_model)   &
                  - edsclrm( iz_CLUBB, 2 )

                 edsclrm( iz_CLUBB, 2 ) = 0.0
             endif
          enddo
        endif !endif do_drp_evap_cf
      endif !use_sclr_HOC
     endif !nqni > 0


   end do ! isub_CLUBB=1, ..., nsub_CLUBB, end of CLUBB calculation

! **********************************************************
!calculate u component tendencies
   call tndy_CLUBB (dtmain, um, u( ix_host_model, iy_host_model, : ),     & ! Intent(in) 
                    udt( ix_host_model, iy_host_model, : ),               & ! Intent(inout) 
                    uten_CLUBB ( ix_host_model, iy_host_model, :) )         ! Intent(out) 

!calculate v component tendencies
   call tndy_CLUBB (dtmain, vm, v( ix_host_model, iy_host_model, : ),  & ! Intent(in) 
                    vdt( ix_host_model, iy_host_model, : ),            & ! Intent(inout) 
                    vten_CLUBB ( ix_host_model, iy_host_model, :) )      ! Intent(out) 

!calculate temperature tendencies
   tten_CLUBB ( ix_host_model, iy_host_model, :) = 0.0       
   call tndy_CLUBB (dtmain, temp_clubb, t( ix_host_model, iy_host_model, : ),  & ! Intent(in)
                    tdt( ix_host_model, iy_host_model, : ),                    & ! Intent(inout)
                    tten_CLUBB ( ix_host_model, iy_host_model, :) )              ! Intent(out)

   if ( nqn > 0 )   then
! tendency for droplet evaporation
        drp_evap_CLUBB   = drp_evap_CLUBB/dtmain
        sclrm_trsport_only(:,1) = sclrm_trsport_only(:,1)/dtmain

! calculate droplet number tendency due to pure transport in CLUBB
!convert from clubb coordinate to host SCM coordinate
        do iz_CLUBB=2+iz_shift, nz_CLUBB, nsublevel+1
           iz_host_model = nz_host_model - (iz_CLUBB -1+iz_shift)/(nsublevel+1)
           qndt_clubb_trsport_only( ix_host_model, iy_host_model,  iz_host_model) = sclrm_trsport_only( iz_CLUBB, 1)
       enddo

   endif !nqn > 0

   if ( nqni > 0 )   then
! tendency for droplet evaporation
        ice_evap_CLUBB   = ice_evap_CLUBB/dtmain
        sclrm_trsport_only(:,2) = sclrm_trsport_only(:,2)/dtmain
! calculate ice number tendency due to pure transport in CLUBB
!convert from clubb coordinate to host SCM coordinate
     do iz_CLUBB=2+iz_shift, nz_CLUBB, nsublevel+1
        iz_host_model = nz_host_model - (iz_CLUBB -1+iz_shift)/(nsublevel+1)
        qnidt_clubb_trsport_only( ix_host_model, iy_host_model,  iz_host_model) = sclrm_trsport_only( iz_CLUBB, 2)
     enddo
   endif !nqni > 0

   ! initiallize eddy diffusivity coefficients to zero
   diff_t_1d = 0.0

! calculate tracer tendencies
   do it = 1, ntp
       if( it == nsphum )  trs_CLUBB =  rtm - rcm
       if( it == nql )     trs_CLUBB =  rcm
       if( it == nqa )     trs_CLUBB =  cloud_frac

!passive tracer evaporation during transport 
       if ( nqn > 0 ) then
         if ( it == nqn )   then
           if( use_sclr_HOC ) then
             trs_CLUBB( : ) = sclrm( : , 1 )
           else
             trs_CLUBB( : ) = edsclrm( : , 1 )
           endif !use_sclr_HOC
           rdt( ix_host_model , iy_host_model , : , it) &
          = rdt( ix_host_model , iy_host_model , : , it) + drp_evap_CLUBB( ix_host_model , iy_host_model , : )
         endif ! it == nqn
       endif ! nqn > 0

       if ( nqni > 0 ) then
         if ( it == nqni )   then
            if( use_sclr_HOC ) then
               trs_CLUBB( : ) = sclrm( : , 2 )
            else
               trs_CLUBB( : ) = edsclrm( : , 2 )
            endif !use_sclr_HOC
            rdt( ix_host_model , iy_host_model , : , it) &
          = rdt( ix_host_model , iy_host_model , : , it) + ice_evap_CLUBB( ix_host_model , iy_host_model , : )

         endif ! it == nqni
       endif ! nqni > 0

!  ---> h1g 12/29/2009
!  water vapor (nsphum), cloud fraction (nqa), droplet number concentration (nqn), ice number concentration (nqni)
!  tendendy in CLUBB will be handled here
       rdt_CLUBB ( ix_host_model, iy_host_model, : , it) = 0.0
       if ( it==nqa ) then
         call tndy_CLUBB (dtmain, trs_CLUBB, r( ix_host_model, iy_host_model, : , it),    & ! Intent(in)
                          rdt( ix_host_model, iy_host_model, : , it),                     & ! Intent(inout)
                          rdt_CLUBB ( ix_host_model, iy_host_model, : , it) )               ! Intent(out)
       endif

       if ( it==nsphum .or. it==nqn .or. it==nqni ) then
         call tndy_CLUBB (dtmain, trs_CLUBB, r( ix_host_model, iy_host_model, : , it),    & ! Intent(in)
                          rdt( ix_host_model, iy_host_model, : , it),                     & ! Intent(inout)
                          rdt_CLUBB ( ix_host_model, iy_host_model, : , it) )               ! Intent(out)
! ---> h1g, scaled by enviromental fraction (1.0 - convective_fraction)
         rdt_CLUBB ( ix_host_model, iy_host_model, : , it) = rdt_CLUBB ( ix_host_model, iy_host_model, : , it)     &
                                                            *env_condensate_scale(ix_host_model, iy_host_model, :)

         rdt( ix_host_model, iy_host_model, : , it)        = rdt_orig( ix_host_model, iy_host_model, : , it)   &
                                                            +rdt_CLUBB ( ix_host_model, iy_host_model, : , it)
! <--- h1g, scaled by enviromental fraction (1.0 - convective_fraction)
       endif

       if ( it==nql ) then
          tmp_host = rdt( ix_host_model, iy_host_model, : , nql) + rdt( ix_host_model, iy_host_model, : , nqi)
          call tndy_CLUBB (dtmain, trs_CLUBB, & 
               r( ix_host_model, iy_host_model, : , nql) + r( ix_host_model, iy_host_model, : , nqi),   & ! Intent(in)
                                      tmp_host,                                                         & ! Intent(inout)
                                      rdt_CLUBB ( ix_host_model, iy_host_model, : , nql) )                ! Intent(out)
               rdt( ix_host_model, iy_host_model, : , nql) = tmp_host - rdt( ix_host_model, iy_host_model, : , nqi)
               rdt_CLUBB ( ix_host_model, iy_host_model, : , nql) = rdt_CLUBB ( ix_host_model, iy_host_model, : , nql)    &
                                                                   *env_condensate_scale(ix_host_model, iy_host_model, :)
               rdt( ix_host_model, iy_host_model, : , nql) = rdt_orig( ix_host_model, iy_host_model, : , nql)             &
                                                            +rdt_CLUBB ( ix_host_model, iy_host_model, : , nql)
       endif 

       ! store eddy diffusivity coefficient
       diff_t_1d = diff_t_1d + Kh_zt

   enddo !do it = 1, ntp

! compute average eddy diffusivity cefficients and copy to host grid
   diff_t_1d = diff_t_1d / nt
   call CLUBB2host_full ( diff_t_1d,                                     & ! Intent(in)
                          diff_t_clubb(ix_host_model, iy_host_model, :))   ! Intent(out)

! dump higher-order terms
   call CLUBB2host_full ( wp3,                                        & ! Intent(in)
                          wp3_CLUBB(ix_host_model, iy_host_model, :))   ! Intent(out)

   call CLUBB2host_half ( wp2,                                        & ! Intent(in)
                          wp2_CLUBB(ix_host_model, iy_host_model, :))   ! Intent(out)

   call CLUBB2host_half ( upwp,                                        & ! Intent(in)
                          upwp_CLUBB(ix_host_model, iy_host_model, :))   ! Intent(out)

   call CLUBB2host_half ( vpwp,                                        & ! Intent(in)
                          vpwp_CLUBB(ix_host_model, iy_host_model, :))   ! Intent(out)

   call CLUBB2host_half ( up2,                                         & ! Intent(in)
                          up2_CLUBB(ix_host_model, iy_host_model, :))    ! Intent(out)
  
   call CLUBB2host_half ( vp2,                                         & ! Intent(in)
                          vp2_CLUBB(ix_host_model, iy_host_model, :))    ! Intent(out)

   call CLUBB2host_half ( wprtp,                                       & ! Intent(in)
                          wprtp_CLUBB(ix_host_model, iy_host_model, :))  ! Intent(out)
 
   call CLUBB2host_half ( wpthlp,                                      & ! Intent(in)
                          wpthlp_CLUBB(ix_host_model, iy_host_model,:))  ! Intent(out)
 
   call CLUBB2host_half ( rtp2,                                        & ! Intent(in)
                          rtp2_CLUBB(ix_host_model, iy_host_model, :))   ! Intent(out)
       
   call CLUBB2host_half ( thlp2,                                       & ! Intent(in)
                          thlp2_CLUBB(ix_host_model, iy_host_model, :))  ! Intent(out)
 
   call CLUBB2host_half ( rtpthlp,                                      & ! Intent(in)
                          rtpthlp_CLUBB(ix_host_model, iy_host_model,:))  ! Intent(out)

   call CLUBB2host_full ( wm_zt,                                      & ! Intent(in)
                          wm_CLUBB(ix_host_model, iy_host_model,:))     ! Intent(out)

! if do droplet
   if ( nqn > 0 )   then
! convert aerosol mass from host from clubb levels
     call host2CLUBB_full( totalmass1(ix_host_model, iy_host_model,:, 1), &    !intent (in)
                           aeromass_clubb (:, 1) )                             !intent (out)
     call host2CLUBB_full( totalmass1(ix_host_model, iy_host_model,:, 2), &    !intent (in)
                           aeromass_clubb (:, 2) )                             !intent (out)
     call host2CLUBB_full( totalmass1(ix_host_model, iy_host_model,:, 3), &    !intent (in)
                           aeromass_clubb (:, 3) )                             !intent (out)
     call host2CLUBB_full( totalmass1(ix_host_model, iy_host_model,:, 4), &    !intent (in)
                           aeromass_clubb (:, 4) )                             !intent (out)

     if (do_aeromass_clubb_const) then
       aeromass_clubb (:, :) = aeromass_clubb_const
     end if
 
   !---> h1g, 2010-07-15
   ! +++++++++++++++++++++++++++++++++ 
   ! ---> temporary turn off the  " do_ice_nucl_wpdf " 
   ! in order to do sensitivity tests of differene aerosol activation scheme, e.g., Diffusion-K, Boucher-Lohmann  

! Option to check for unphysical temperature values
     if (icheck_temp .gt. 0 ) then
       do k=1,nz_host_model
         if (temp_clubb(k).lt.-150.0+273.15 .or. temp_clubb(k).gt.90+273.15) then
           write(*,'(a,f6.2,3i4,2f16.8)') 'CLUBB_driver: bad temperature ',   &
            temp_clubb(k),ix_host_model,iy_host_model,k,  &
            lon(ix_host_model,iy_host_model),lat(ix_host_model,iy_host_model)
         end if
       end do
     end if
 
     IF ( do_ice_nucl_wpdf ) THEN
! have both warm and ice nucleation 
            call  warm_ice_nucl_mns_clubb( aeromass_clubb,     &    !  Intent(in)
                                           temp_clubb,         &    !  Intent(in)
                                           Ndrop_max,          &    !  Intent(out)
                                           Ncrystal_max,       &    !  Intent(out)
                                           RH_crit  )               !  Intent(out)
     else
           Ncrystal_max = 0.0
           RH_crit = 1.0

           if( do_quadrature_gauss )   &
           call aer_act_clubb_quadrature_Gauss(Time_next, aeromass_clubb, temp_clubb, &   ! Intent(in)
                                               Ndrop_max)                                 ! Intent(out)

           if( do_diffK_gauss )   &
           call aer_act_clubb_diffK_Gauss( aeromass_clubb, temp_clubb,                &   ! Intent(in)
                                           Ndrop_max )                                    ! Intent(out)

           if( do_BL_gauss )   &
           call  aer_act_clubb_BL_Gauss( aeromass_clubb,        &   ! Intent(in)
                                         Ndrop_max )                ! Intent(out)

      endif !do_ice_nucl_wpdf
      !  +++++++++++++++++++++++++++++++++ 
      ! <--- h1g, 2010-07-15

! convert the unit from #/cm3 to #/kg air 
      Ndrop_max = Ndrop_max * 1.0e6 / rho
   
! update tracer concentrations: droplet # concentration
      do iz_CLUBB =2, nz_CLUBB

! ---> h1g, 2011-04-20,   no liquid drop nucleation if T < -40 C
         if( temp_clubb(iz_CLUBB) <= 233.15 )  Ndrop_max( iz_CLUBB ) = 0.0  ! if T<-40C, no liquid drop nucleation
! <--- h1g, 2011-04-20

         if(rcm(iz_CLUBB) > 0.0 .and. cloud_frac(iz_CLUBB) > cloud_frac_min) then
            if( use_sclr_HOC ) then
                sclrm( iz_CLUBB, 1 ) = max(Ndrop_max( iz_CLUBB ), sclrm( iz_CLUBB, 1 ) )
            else
                edsclrm( iz_CLUBB, 1 ) = max(Ndrop_max( iz_CLUBB ), edsclrm( iz_CLUBB, 1 ) )
            endif
         endif
      enddo

      do iz_CLUBB =2+iz_shift, nz_CLUBB, nsublevel+1
         iz_host_model = nz_host_model - (iz_CLUBB -1+iz_shift)/(nsublevel+1)
         if(rcm(iz_CLUBB) > 0.0 .and. cloud_frac(iz_CLUBB) > cloud_frac_min) then
!  Ndrop_act_CLUBB:    in-cloud  averaged activated droplet number concentration (in the unit of #/kg)
           Ndrop_act_CLUBB(ix_host_model, iy_host_model, iz_host_model) = Ndrop_max(iz_CLUBB)/cloud_frac(iz_CLUBB)
         endif
      enddo

      if( use_sclr_HOC ) then
         trs_CLUBB(:) =  sclrm(:,1)
      else
         trs_CLUBB(:) =  edsclrm(:,1)
      endif

      call tndy_CLUBB (dtmain, trs_CLUBB, r(ix_host_model, iy_host_model,:, nqn),  & ! Intent(in)
                       rdt( ix_host_model, iy_host_model, : , nqn),                & ! Intent(inout)
                       aer_ccn_act_CLUBB( ix_host_model, iy_host_model,:))           ! Intent(out)

!calculate drop concentration net flux
      drp_flux_CLUBB(  ix_host_model, iy_host_model, 1 ) = 0.0
      do  iz_host_model = 2, nz_host_model
          iz_CLUBB = (nz_host_model - iz_host_model)*(nsublevel+1) - iz_shift + 1 + 1

          if(  iz_CLUBB  .gt.  (nz_CLUBB - nsublevel -1) ) then
               drp_flux_CLUBB(  ix_host_model, iy_host_model,  iz_host_model ) = 0.0
          else
               sum_sclrm_trsport = 0.0
               do isum_sclrm = iz_CLUBB+1,  iz_CLUBB+nsublevel+1
                   sum_sclrm_trsport = sum_sclrm_trsport + sclrm_trsport_only( isum_sclrm , 1)
               enddo ! isum_sclrm
 
               drp_flux_CLUBB(  ix_host_model, iy_host_model, iz_host_model ) =           &
                            drp_flux_CLUBB(  ix_host_model, iy_host_model, iz_host_model-1 )   &
                             + sum_sclrm_trsport/( nsublevel+1 )                                              &
                               *( zhalf( ix_host_model, iy_host_model, iz_host_model-1) - zhalf( ix_host_model, iy_host_model, iz_host_model) )
          endif !   iz_CLUBB  .gt.  (nz_CLUBB - nsublevel -1)
      enddo !  iz_host_model

      IF ( do_ice_nucl_wpdf ) THEN
! convert the unit from #/cm3 to #/kg air for ice crystal number concentration
         Ncrystal_max =  Ncrystal_max / rho
  
         do iz_CLUBB =2, nz_CLUBB
            if(rcm(iz_CLUBB) > 0.0 .and. cloud_frac(iz_CLUBB) > cloud_frac_min) then
                if( use_sclr_HOC ) then
                    sclrm( iz_CLUBB, 2 ) = max(Ncrystal_max( iz_CLUBB ), sclrm( iz_CLUBB, 2 ) )
                else
                    edsclrm( iz_CLUBB, 2 ) = max(Ncrystal_max( iz_CLUBB ), edsclrm( iz_CLUBB, 2 ) )
               endif
            endif
         enddo

!   **********  in current version   ********** 
         do iz_CLUBB =2+iz_shift, nz_CLUBB, nsublevel+1
            iz_host_model = nz_host_model - (iz_CLUBB -1+iz_shift)/(nsublevel+1)  
            if(rcm(iz_CLUBB) > 0.0 .and. cloud_frac(iz_CLUBB) > cloud_frac_min) then
! Icedrop_act_CLUBB:    in-cloud  averaged nucleated ice-crystal number concentration (in the unit of #/kg)  
                Icedrop_act_CLUBB( ix_host_model, iy_host_model, iz_host_model ) = Ncrystal_max( iz_CLUBB )/cloud_frac(iz_CLUBB) 
            endif
         enddo

         if( use_sclr_HOC ) then
             trs_CLUBB( : ) =  sclrm( : , 2 )
         else
             trs_CLUBB( : ) =  edsclrm( : , 2 )
         endif

         aer_ice_act_CLUBB( ix_host_model, iy_host_model, : )  = 0.0
          
         !call tndy_CLUBB (dtmain, trs_CLUBB, r( ix_host_model, iy_host_model, : , nqni),    & ! Intent(in)
         !                 rdt( ix_host_model, iy_host_model, : , nqni),                     & ! Intent(inout)
         !                 aer_ice_act_CLUBB( ix_host_model, iy_host_model, : ) )              ! Intent(out)

!calculate ice crystal concentration flux
         Icedrop_flux_CLUBB(  ix_host_model, iy_host_model, 1 ) = 0.0
         do  iz_host_model = 2, nz_host_model
             iz_CLUBB = (nz_host_model - iz_host_model)*(nsublevel+1) - iz_shift + 1 + 1

             if( iz_CLUBB  .gt.  (nz_CLUBB - nsublevel -1) ) then
                 Icedrop_flux_CLUBB(  ix_host_model, iy_host_model,  iz_host_model ) = 0.0
             else
                 sum_sclrm_trsport = 0.0
                 do isum_sclrm = iz_CLUBB+1,  iz_CLUBB+nsublevel+1
                    sum_sclrm_trsport = sum_sclrm_trsport + sclrm_trsport_only( isum_sclrm , 2)
                 enddo ! isum_sclrm
 
                Icedrop_flux_CLUBB(  ix_host_model, iy_host_model, iz_host_model ) =           &
                            Icedrop_flux_CLUBB(  ix_host_model, iy_host_model, iz_host_model-1 )   &
                           + sum_sclrm_trsport/( nsublevel+1 )                                              &
                           *( zhalf( ix_host_model, iy_host_model, iz_host_model-1) - zhalf( ix_host_model, iy_host_model, iz_host_model) )
             endif !   iz_CLUBB  .gt.  (nz_CLUBB - nsublevel -1)
         enddo !  iz_host_model
      endif !do_ice_nucl_wpdf

   endif  ! nqn > 0
!end passive tracers

   call CLUBB_1D_2_3D( ix_CLUBB, iy_CLUBB )

  enddo ! enddo of  ix_host_model
enddo  ! enddo of iy_host_model

!-->cjg: test option to override CLUBB cloud fraction for ice clouds
if (do_alt_cloud > 0) then

  call alt_cloud( do_alt_cloud,dtmain,pfull,  &
                  t,q,r(:,:,:,nqa),r(:,:,:,nql),r(:,:,:,nqi), &
                  tdt,qdt,rdt(:,:,:,nqa),rdt(:,:,:,nql),rdt(:,:,:,nqi) )

endif
!<--cjg

if (id_udt_CLUBB > 0) then
    used = send_data ( id_udt_CLUBB,   uten_CLUBB, &
                       Time_next, is, js, 1 )
endif

if (id_vdt_CLUBB > 0) then
    used = send_data ( id_vdt_CLUBB,   vten_CLUBB, &
                       Time_next, is, js, 1 )
endif

do it = 1, ntp
   if ( id_tracer_CLUBB(it)>0 ) then
        used = send_data ( id_tracer_CLUBB(it),   rdt_CLUBB ( : , : , : , it),   &
                           Time_next, is, js, 1 )
   endif 
enddo

if( do_CLUBB_conservation_checks )  then
    do iz_host_model = 1,  nz_host_model - 1
       do iy_host_model =  js, je
         do ix_host_model = is, ie
            iy_CLUBB =   iy_host_model -  js + 1
            ix_CLUBB =   ix_host_model -  is + 1
            pmass_3d( ix_CLUBB, iy_CLUBB,  iz_host_model) = (phalf( ix_host_model, iy_host_model,  iz_host_model+1) &
                                                             - phalf( ix_host_model, iy_host_model, iz_host_model) )/grav
         enddo
       enddo
    enddo

    if ( id_wat_CLUBB_col > 0 ) then
         tmp2d_check = 0.0
         do iz_host_model =  1,  nz_host_model - 1
           do  iy_host_model =  js, je
             do  ix_host_model =  is, ie
                 iy_CLUBB =   iy_host_model -  js + 1
                 ix_CLUBB =   ix_host_model -  is + 1
                 tmp2d_check(  ix_CLUBB, iy_CLUBB ) =   tmp2d_check(  ix_CLUBB, iy_CLUBB ) &
                      +   pmass_3d(ix_CLUBB, iy_CLUBB,  iz_host_model) &
                         *(   rdt_CLUBB(  ix_host_model, iy_host_model,  iz_host_model , nsphum) &
                            + rdt_CLUBB(  ix_host_model, iy_host_model,  iz_host_model , nql) )
             enddo
           enddo
         enddo
         if (id_wat_CLUBB_col > 0 ) &
             used = send_data ( id_wat_CLUBB_col, tmp2d_check,   &
                                Time_next, is, js )
         endif

         if (id_enth_CLUBB_col > 0 ) then
             tmp2d_check = 0.0
             do iz_clubb= 2+iz_shift, nz_CLUBB, nsublevel+1
                iz_host_model = nz_host_model - (iz_clubb -1+iz_shift)/(nsublevel+1)
                do  iy_host_model =  js, je 
                  do  ix_host_model =  is, ie
                    iy_CLUBB =   iy_host_model -  js + 1
                    ix_CLUBB =   ix_host_model -  is + 1
                    tmp2d_check(  ix_CLUBB, iy_CLUBB ) =   tmp2d_check(  ix_CLUBB, iy_CLUBB ) &
                      +   pmass_3d(ix_CLUBB, iy_CLUBB,  iz_host_model) &
                         *(   tten_CLUBB(  ix_host_model, iy_host_model,  iz_host_model )*Cp &
                             - rdt_CLUBB(  ix_host_model, iy_host_model,  iz_host_model , nql)*Lv )
                  enddo
                enddo
              enddo

! ---> h1g    Very slightly adjust tendencies to force exact   ***
!                     enthalpy  conservation, 2010-06-15
              do  iy_host_model =  js, je 
                do  ix_host_model =  is, ie
                    iy_CLUBB =   iy_host_model -  js + 1
                    ix_CLUBB =   ix_host_model -  is + 1

! temperature fix is only applied for the CLUBB vertical domain
                     tten_CLUBB_fix =  tmp2d_check( ix_CLUBB, iy_CLUBB ) & 
                                       / ( phalf(ix_host_model, iy_host_model, nz_host_model) - phalf(ix_host_model, iy_host_model, iz_host_model)) &
                                       * grav / Cp

                     tmp2d_check(ix_CLUBB, iy_CLUBB) = 0.0

! recalculate column enthalpy tendency
                     do  iz_clubb= 2+iz_shift, nz_CLUBB, nsublevel+1
                       iz_host_model = nz_host_model - (iz_CLUBB -1+iz_shift)/(nsublevel+1)   
                       tten_CLUBB(ix_host_model, iy_host_model, iz_host_model) = tten_CLUBB(ix_host_model, iy_host_model, iz_host_model) - tten_CLUBB_fix

! add temperature fix to the temperature tendency  "tdt"
                       tdt(ix_host_model, iy_host_model, iz_host_model) = tdt(ix_host_model, iy_host_model, iz_host_model) - tten_CLUBB_fix

                       tmp2d_check(ix_CLUBB, iy_CLUBB) = tmp2d_check(ix_CLUBB, iy_CLUBB) &
                          +   pmass_3d(ix_CLUBB, iy_CLUBB,  iz_host_model) &
                             *(   tten_CLUBB( ix_host_model, iy_host_model, iz_host_model)*Cp &
                                 - rdt_CLUBB( ix_host_model, iy_host_model, iz_host_model, nql)*Lv)
                    enddo
                enddo
              enddo
! <--- h1g,  2010-06-15 

              used = send_data(id_enth_CLUBB_col, tmp2d_check, Time_next, is, js)
             endif

             if (id_tdt_CLUBB > 0) then
                used = send_data ( id_tdt_CLUBB,   tten_CLUBB, Time_next, is, js, 1 )
             endif

            do it = 1, ntp
               tmp2d_check = 0.0
              do iz_host_model =  1,  nz_host_model - 1
                do  iy_host_model =  js, je 
                  do  ix_host_model =  is, ie
                    iy_CLUBB = iy_host_model -  js + 1
                    ix_CLUBB = ix_host_model -  is + 1
                    tmp2d_check(ix_CLUBB, iy_CLUBB) =   tmp2d_check(  ix_CLUBB, iy_CLUBB ) &
                      +   pmass_3d(ix_CLUBB, iy_CLUBB,  iz_host_model) &
                         * rdt_CLUBB(ix_host_model, iy_host_model, iz_host_model, it)
                  enddo
                enddo
              enddo

              if ( id_tracer_CLUBB_col(it)>0 ) &
              used = send_data (id_tracer_CLUBB_col(it), tmp2d_check, Time_next, is, js)
              enddo
 
         endif

!        Adjust cloud fraction for ice clouds


         if( nqn > 0 )   then
           if ( id_drp_evap_CLUBB>0) then
                used = send_data(id_drp_evap_CLUBB, drp_evap_CLUBB, Time_next, is, js, 1)
           endif

           if ( id_aer_ccn_act_CLUBB > 0) then
                used = send_data(id_aer_ccn_act_CLUBB, aer_ccn_act_CLUBB, Time_next, is, js, 1)
           endif

           if ( id_Ndrop_act_CLUBB > 0) then
                used = send_data(id_Ndrop_act_CLUBB, Ndrop_act_CLUBB, Time_next, is, js, 1 )
           endif

!calculate droplet number turbulence flux (m/cm3/s)
           if ( id_drp_flux_CLUBB > 0 ) then
                used = send_data(id_drp_flux_CLUBB, drp_flux_CLUBB, Time_next, is, js, 1 )
           endif !  id_drp_flux_CLUBB > 0

           if ( id_qndt_CLUBB_trsport_only > 0 ) then
                used = send_data(id_qndt_CLUBB_trsport_only, qndt_CLUBB_trsport_only, Time_next, is, js, 1)
           endif
         endif   !  nqn > 0
  

         if( nqni>0 )   then
           if ( id_ice_evap_CLUBB > 0) then
                used = send_data(id_ice_evap_CLUBB, ice_evap_CLUBB, Time_next, is, js, 1)
           endif

           if ( id_aer_ice_act_CLUBB>0 ) then
                 used = send_data(id_aer_ice_act_CLUBB, aer_ice_act_CLUBB, Time_next, is, js, 1)
           endif

           if ( id_Icedrop_act_CLUBB>0 ) then
                used = send_data(id_Icedrop_act_CLUBB, Icedrop_act_CLUBB, Time_next, is, js, 1)
           endif

!calculate droplet number turbulence flux (m/cm3/s)
           if ( id_icedrop_flux_CLUBB > 0 ) then
                used = send_data(id_icedrop_flux_CLUBB, icedrop_flux_CLUBB, Time_next, is, js, 1)
           endif !  id_drp_flux_CLUBB > 0

           if ( id_qnidt_CLUBB_trsport_only > 0 ) then
                used = send_data ( id_qnidt_CLUBB_trsport_only, qnidt_CLUBB_trsport_only, Time_next, is, js, 1)
           endif
         endif   !  nqni > 0

!output higher order terms and fluxes
!at full levels
         if ( id_wp3_CLUBB > 0 ) then
          used = send_data(id_wp3_CLUBB, wp3_CLUBB, Time_next, is, js, 1)
         endif

         if ( id_wm_CLUBB > 0 ) then
          used = send_data(id_wm_CLUBB, wm_CLUBB, Time_next, is, js, 1)
         endif

         if ( id_omega_CLUBB > 0 ) then
          used = send_data(id_omega_CLUBB, omega, Time_next, is, js, 1)
         endif

!at half levels
         if ( id_wp2_CLUBB > 0 ) then
          used = send_data ( id_wp2_CLUBB, wp2_CLUBB, Time_next, is, js, 1)
         endif

         if ( id_upwp_CLUBB > 0 ) then
           used = send_data ( id_upwp_CLUBB, upwp_CLUBB, Time_next, is, js, 1)
         endif

         if ( id_vpwp_CLUBB > 0 ) then
           used = send_data ( id_vpwp_CLUBB, vpwp_CLUBB, Time_next, is, js, 1)
         endif

         if ( id_up2_CLUBB > 0 ) then
           used = send_data ( id_up2_CLUBB, up2_CLUBB, Time_next, is, js, 1)
         endif

         if ( id_vp2_CLUBB > 0 ) then
           used = send_data ( id_vp2_CLUBB, vp2_CLUBB, Time_next, is, js, 1)
         endif

         if ( id_wprtp_CLUBB > 0 ) then
           used = send_data ( id_wprtp_CLUBB, wprtp_CLUBB, Time_next, is, js, 1)
         endif

         if ( id_wpthlp_CLUBB > 0 ) then
           used = send_data ( id_wpthlp_CLUBB, wpthlp_CLUBB, Time_next, is, js, 1)
         endif

         if ( id_rtp2_CLUBB > 0 ) then
           used = send_data ( id_rtp2_CLUBB, rtp2_CLUBB, Time_next, is, js, 1)
         endif

         if ( id_thlp2_CLUBB > 0 ) then
           used = send_data ( id_thlp2_CLUBB, thlp2_CLUBB, Time_next, is, js, 1)
         endif

         if ( id_rtpthlp_CLUBB > 0 ) then
           used = send_data ( id_rtpthlp_CLUBB, rtpthlp_CLUBB, Time_next, is, js, 1)
         endif

    return
    end subroutine clubb
!#########################################################



 
SUBROUTINE CLUBB_INIT(id, jd,  kd,  lon, lat,  &
                      axes, Time,   phalf   )
!=======================================================================
! ***** INITIALIZE and set-upCLUBB
!=======================================================================
!-----------------------------------------------------------------------
!         is:     index for "longitude"
!         js:     index for "latitude"
!         kd:          vertical dimension
!         lon:       model longitudes at cell center [radians]
!         lat:       model latitudes cell center [radians]

!         axes       axis indices, (/x,y,pf,ph/)
!                (returned from diag axis manager)
!         Time       current time (time_type)
!
!         phalf       pressure at half levels in pascals
!                       [real, dimension(nlon,nlat,nlev+1)]
!-----------------------------------------------------------------------     
IMPLICIT NONE
integer,              intent(in)     :: id, jd, kd
real,dimension(:,:),  intent(in)     :: lon, lat

integer,dimension(4),   intent(in)   :: axes
type(time_type),        intent(in)   :: Time
    
real, intent(in), dimension(:,:,:)   :: phalf
 

integer, dimension(3) :: half = (/1,2,4/)
INTEGER :: ierr, io, unit            !open namelist error

integer :: i, j, n

Integer :: isclr_dim
character(len=64)            :: fname='INPUT/CLUBB.res.nc'
character(len=64)            :: sclr_prefix,  sclr_txt
integer :: year, month, day, hour, minute, second
real(kind=time_precision) :: time_current
type(time_type)           :: Time_init

!-----------------------------------------------------------------------------------------------------------------
! --- Begin Code --- !
if(module_is_initialized) then
   return
else
   module_is_initialized = .true.
endif

clubb_core_clock = mpp_clock_id( '   clubb_driver: core',    &
                                grain=CLOCK_MODULE_DRIVER )

!---------------------------------------------------------------------
! --- Read namelist
!---------------------------------------------------------------------
if( FILE_EXIST( 'input.nml' ) ) then
! -------------------------------------
    unit = OPEN_NAMELIST_FILE ( )
    ierr = 1        
    do while( ierr /= 0 )
       READ ( unit,  nml = CLUBB_setting_nml, iostat = io, end = 10 ) 
       ierr = check_nml_error (io, 'CLUBB_setting_nml')
    end do
10  CALL CLOSE_FILE( unit )
! -------------------------------------

!cjg begins
    unit = OPEN_NAMELIST_FILE ( )
    ierr = 1
    do while( ierr /= 0 )
       READ ( unit,  nml = CLUBB_stats_setting_nml, iostat = io, end = 20 ) 
       ierr = check_nml_error (io, 'CLUBB_stats_setting_nml')
    end do
20  CALL CLOSE_FILE( unit )
! -------------------------------------
!cjg ends
end if

! ---> h1g, 2010-08-24, get initial time
Time_init = set_date ( init_date(1),  init_date(2),  init_date(3),  init_date(4),  init_date(5),   init_date(6) )
call get_time( Time_init, current_sec0, current_days0)
! <--- h1g, 2010-08-24

! get tracer indices
nsphum = get_tracer_index ( MODEL_ATMOS, 'sphum' )
nql    = get_tracer_index ( MODEL_ATMOS, 'liq_wat' )
nqi    = get_tracer_index ( MODEL_ATMOS, 'ice_wat' )
nqa    = get_tracer_index ( MODEL_ATMOS, 'cld_amt' )
nqn    = get_tracer_index ( MODEL_ATMOS, 'liq_drp' )
nqni   = get_tracer_index ( MODEL_ATMOS, 'ice_num' )

if( do_liq_num .and. nqn<= 0 )  &
    call error_mesg ('CLUBB_driver_mod','do_liq_num is true, but nqn<= 0', FATAL)
 
!--------- write version number and namelist --------
call write_version_number ( version, tagname )
if(mpp_pe() == mpp_root_pe() ) write(stdlog(),nml=CLUBB_setting_nml)

! set-up CLUBB
call CLUBB_SETUP(id, jd, phalf)

if (file_exist('INPUT/CLUBB.res.nc') ) then
         if(mpp_pe() == mpp_root_pe() ) call mpp_error ('CLUBB_mod', &
         'Reading netCDF formatted restart file: INPUT/CLUBB.res.nc', NOTE)
           call read_data (fname, 'wp2_3D',           wp2_3D)
           call read_data (fname, 'wp3_3D',           wp3_3D)
           call read_data (fname, 'upwp_3D',          upwp_3D)
           call read_data (fname, 'vpwp_3D',          vpwp_3D)
           call read_data (fname, 'wprtp_3D',         wprtp_3D)
           call read_data (fname, 'wpthlp_3D',        wpthlp_3D)
           call read_data (fname, 'rtp2_3D',          rtp2_3D)
           call read_data (fname, 'thlp2_3D',         thlp2_3D)
           call read_data (fname, 'rtpthlp_3D',       rtpthlp_3D)
           call read_data (fname, 'up2_3D',           up2_3D)
           call read_data (fname, 'vp2_3D',           vp2_3D)

            if( sclr_dim > 0 ) then
              call read_data (fname, 'RH_crit_3D1',   RH_crit_clubb_3D( : , : , : ,1) )
              call read_data (fname, 'RH_crit_3D2',   RH_crit_clubb_3D( : , : , : ,2) )
            endif ! sclr_dim > 0

endif ! do_netcdf_restart

call aer_ccn_act_init

!set-up diagosis 
call get_number_tracers(MODEL_ATMOS, num_tracers=nt, num_prog=ntp)

if( do_CLUBB_conservation_checks  ) then
    if( .not. allocated(  id_tracer_CLUBB ) ) then
         allocate ( id_tracer_CLUBB(nt) ) 
    else
         deallocate(  id_tracer_CLUBB )
         allocate ( id_tracer_CLUBB(nt) )
    endif


    if( .not. allocated( id_tracer_CLUBB_col ) ) then
         ALLOCATE ( id_tracer_CLUBB_col(nt) )
    else
         deallocate( id_tracer_CLUBB_col )
         ALLOCATE( id_tracer_CLUBB_col(nt) )
    endif


    if( .not. allocated( pmass_3d ) ) then
         ALLOCATE( pmass_3d( size(phalf,1),   size(phalf,2),  size(phalf,3)-1 ) )
    else
         deallocate( pmass_3d )
         ALLOCATE( pmass_3d( size(phalf,1),   size(phalf,2),  size(phalf,3)-1 ) )
    endif

    if( .not. allocated( tmp2d_check ) ) then
         allocate( tmp2d_check( size(phalf,1),   size(phalf,2) ) )
    else
         deallocate( tmp2d_check )
         allocate( tmp2d_check( size(phalf,1),   size(phalf,2) ) )
    endif     
endif

id_udt_CLUBB = register_diag_field ( mod_name,    &
   'udt_CLUBB', axes(1:3), Time,                  &
   'u-component tendency from CLUBB (boundary layer cloud)', &
   'm/s/s', missing_value=missing_value)

id_vdt_CLUBB = register_diag_field ( mod_name,    &
  'vdt_CLUBB', axes(1:3), Time,                   &
  'v-component tendency from CLUBB (boundary layer cloud)', &
  'm/s/s', missing_value=missing_value)

id_tdt_CLUBB = register_diag_field ( mod_name,    &
  'tdt_CLUBB', axes(1:3), Time,                   &
  'temperature tendency from CLUBB (boundary layer cloud)', &
   'K/s', missing_value=missing_value)
 
do it = 1, ntp
   call get_tracer_names (MODEL_ATMOS, it, name = tracer_name,  &
                          units = tracer_units)
        
   diaglname = trim(tracer_name)//  &
              ' tendency from CLUBB (boundary layer cloud)'
   id_tracer_CLUBB(it) =    &
        register_diag_field ( mod_name, &
        TRIM(tracer_name)//'dt_CLUBB',  &
        axes(1:3), Time, trim(diaglname), &
        TRIM(tracer_units)//'/s',  &
        missing_value=missing_value)
enddo

do it = 1, ntp
   call get_tracer_names (MODEL_ATMOS, it, name = tracer_name,  &
                          units = tracer_units)

   diaglname = trim(tracer_name)//' column tendency from CLUBB'  

   id_tracer_CLUBB_col(it) =    &
                         register_diag_field ( mod_name, &
                         TRIM(tracer_name)//'_CLUBB_col',  &
                         axes(1:2), Time, trim(diaglname), &
                         TRIM(tracer_units)//'/s',  &
                         missing_value=missing_value) 
enddo

         id_enth_CLUBB_col = register_diag_field ( mod_name, &
        'enth_clubb_col', axes(1:2), Time, &
        'Column enthalpy tendency from CLUBB','W/m2' )
 
         id_wat_CLUBB_col = register_diag_field ( mod_name, &
       'wat_clubb_col', axes(1:2), Time, &
       'Column total water tendency from CLUBB','kg/m2/s' )

         id_drp_evap_CLUBB = register_diag_field ( mod_name,    &
         'drp_evapdt_CLUBB', axes(1:3), Time,                  &
         'droplet evaporation rate from CLUBB (boundary layer cloud)', &
         '#/kg/s', missing_value=missing_value)

         id_aer_ccn_act_CLUBB = register_diag_field ( mod_name,    &
         'aer_ccn_actdt_CLUBB', axes(1:3), Time,                  &
         'droplet activation rate from CLUBB (boundary layer cloud)', &
         '#/kg/s', missing_value=missing_value)

         id_Ndrop_act_CLUBB = register_diag_field ( mod_name,    &
         'Ndrop_act_CLUBB', axes(1:3), Time,                  &
         'Maximum in-cloud droplet activation number concentration in CLUBB (boundary layer cloud)', &
         '#/kg', missing_value=missing_value)

         id_drp_flux_CLUBB = register_diag_field ( mod_name,    &
         'drp_flux_CLUBB', axes(half), Time,                  &
         'droplet number turbulent flux from CLUBB (boundary layer cloud)', &
         'm/kg/s', missing_value=missing_value)

         id_qndt_CLUBB_trsport_only = register_diag_field ( mod_name,    &
         'qndt_CLUBB_trsport_only', axes(1:3), Time,                  &
         'tendency of droplet number from pure transport in CLUBB (boundary layer cloud)', &
         'm/kg/s', missing_value=missing_value)

         id_ice_evap_CLUBB = register_diag_field ( mod_name,    &
         'ice_evapdt_CLUBB', axes(1:3), Time,                  &
         'ice crystal evaporation rate from CLUBB (boundary layer cloud)', &
         '#/kg/s', missing_value=missing_value)

         id_aer_ice_act_CLUBB = register_diag_field ( mod_name,    &
         'aer_ice_actdt_CLUBB', axes(1:3), Time,                  &
         'ice nucleation rate from CLUBB (boundary layer cloud)', &
         '#/kg/s', missing_value=missing_value)

         id_Icedrop_act_CLUBB = register_diag_field ( mod_name,    &
         'Icedrop_act_CLUBB', axes(1:3), Time,                  &
         'Maximum ice nucleation number concentration in CLUBB (boundary layer cloud)', &
         '#/kg', missing_value=missing_value)

         id_icedrop_flux_CLUBB = register_diag_field ( mod_name,    &
         'icedrop_flux_CLUBB', axes(half), Time,                  &
         'ice number turbulent flux from CLUBB (boundary layer cloud)', &
         'm/kg/s', missing_value=missing_value)

         id_qnidt_CLUBB_trsport_only = register_diag_field ( mod_name,    &
         'qnidt_CLUBB_trsport_only', axes(1:3), Time,                  &
         'tendency of ice number from pure transport in CLUBB (boundary layer cloud)', &
         'm/kg/s', missing_value=missing_value)

         id_sulfate = register_diag_field ( mod_name, 'sulfate', &
         axes(1:3), Time, 'sulfate mass conentration',     &
         'ug so4/m3', missing_value=missing_value )
 
         id_seasalt_sub = register_diag_field ( mod_name, 'seasalt_sub', &
         axes(1:3), Time, 'sub-micron sea salt mass conentration',     &
         'ug/m3', missing_value=missing_value )
 
         id_seasalt_sup = register_diag_field ( mod_name, 'seasalt_sup', &
         axes(1:3), Time, 'super-micron sea salt mass conentration',     &
         'ug/m3', missing_value=missing_value )
 
         id_om = register_diag_field ( mod_name, 'OM', &
         axes(1:3), Time, 'OM mass conentration',     &
         'ug/m3', missing_value=missing_value )


!register higher order terms and fluxes
!at full levels
         id_wp3_CLUBB = register_diag_field ( mod_name,    &
         'wp3_CLUBB', axes(1:3), Time,                  &
         'w third order moment from CLUBB (boundary layer cloud)', &
         'm^3/s^3', missing_value=missing_value)

         id_wm_CLUBB = register_diag_field ( mod_name,    &
         'wm_CLUBB', axes(1:3), Time,                  &
         'vertical velocity W adjusted by (donner) convective mass flux', &
         'm/s', missing_value=missing_value)

         id_omega_CLUBB = register_diag_field ( mod_name,    &
         'omega_CLUBB', axes(1:3), Time,                  &
         'vertical velocity Omega adjusted by (donner) convective mass flux', &
         'Pa/s', missing_value=missing_value)

!at half levels
         id_wp2_CLUBB = register_diag_field ( mod_name,    &
         'wp2_CLUBB', axes(half), Time,                  &
         'w second order moment from CLUBB (boundary layer cloud) at half levels', &
         'm^2/s^2', missing_value=missing_value)

         id_upwp_CLUBB = register_diag_field ( mod_name,    &
         'upwp_CLUBB', axes(half), Time,                  &
         'u and w co-variance from CLUBB (boundary layer cloud) at half levels', &
         'm^2/s^2', missing_value=missing_value)

         id_vpwp_CLUBB = register_diag_field ( mod_name,    &
         'vpwp_CLUBB', axes(half), Time,                  &
         'v and w co-variance from CLUBB (boundary layer cloud) at half levels', &
         'm^2/s^2', missing_value=missing_value)

         id_up2_CLUBB = register_diag_field ( mod_name,    &
         'up2_CLUBB', axes(half), Time,                  &
         'u second order moment from CLUBB (boundary layer cloud) at half levels', &
         'm^2/s^2', missing_value=missing_value)
 
         id_vp2_CLUBB = register_diag_field ( mod_name,    &
         'vp2_CLUBB', axes(half), Time,                  &
         'v second order moment from CLUBB (boundary layer cloud) at half levels', &
         'm^2/s^2', missing_value=missing_value)

         id_wprtp_CLUBB = register_diag_field ( mod_name,    &
         'wprtp_CLUBB', axes(half), Time,                  &
         'w and rtm co-variance from CLUBB (boundary layer cloud) at half levels', &
         'm/s', missing_value=missing_value)

         id_wpthlp_CLUBB = register_diag_field ( mod_name,    &
         'wpthlp_CLUBB', axes(half), Time,                  &
         'w and thlm co-variance from CLUBB (boundary layer cloud) at half levels', &
         'm/s K', missing_value=missing_value)

         id_rtp2_CLUBB = register_diag_field ( mod_name,    &
         'rtp2_CLUBB', axes(half), Time,                  &
         'rtm  second order moment from CLUBB (boundary layer cloud) at half levels', &
         ' ', missing_value=missing_value)

         id_thlp2_CLUBB = register_diag_field ( mod_name,    &
         'thlp2_CLUBB', axes(half), Time,                  &
         'thlm second order moment from CLUBB (boundary layer cloud) at half levels', &
         'K^2', missing_value=missing_value)

         id_rtpthlp_CLUBB = register_diag_field ( mod_name,    &
         'rtpthlp_CLUBB', axes(half), Time,                  &
         'rtm and thlm co-variance from CLUBB (boundary layer cloud) at half levels', &
         'K^2', missing_value=missing_value)

! Initialize CLUBB internal stats
if (do_stats) then

  ! For GCM identify column of interest, for SCM set to 1,1
#ifdef SCM
  istats = 1
  jstats = 1
  fstats = 1
#else
  do i=1,id
    do j=1,jd
      do n=1,size(lon_stats)
        if ( abs(lon(i,j)*RAD_TO_DEG-lon_stats(n)).lt.0.5 .and.  &
             abs(lat(i,j)*RAD_TO_DEG-lat_stats(n)).lt.0.5 ) then
          istats = i
          jstats = j
          fstats = n
        end if
      end do
    end do
  end do
#endif

  ! If stats is active on current PE
  if ( istats.ge.0 .and. jstats.ge.0 ) then

    write(*,'(a,2i4,2f12.4,a)') 'CLUBB_driver: stats active for ',   &
    istats, jstats,                                                  &
    RAD_TO_DEG*lon(istats,jstats), RAD_TO_DEG*lat(istats,jstats),    &
    trim(fname_prefix(fstats))

    unit_stats = open_namelist_file()
    close(unit_stats)

    Time_stats_init = Time
    call get_date( Time_stats_init, year, month, day, hour, minute, second )
    time_current = 3600.0*hour + 60.0*minute + second

    call stats_init( unit_stats, fname_prefix(fstats), "",                       &
                     do_stats, stats_fmt, stats_tsamp, stats_tout, 'input.nml',  &
                     gr%nz, gr%zt, gr%zm,                                      &
                     gr%nz, gr%zt, gr%nz, gr%zm,                             &
                     day, month, year,                                           & 
                     RAD_TO_DEG*lat(istats:istats,jstats),                       &
                     RAD_TO_DEG*lon(istats:istats,jstats),                       &
                     time_current, CLUBB_dt)

  end if

end if

return
end subroutine CLUBB_INIT
!#######################################################################



!#######################################################################
SUBROUTINE CLUBB_SETUP(id, jd, phalf)
!=======================================================================
! ***** set-up CLUBB
!=======================================================================
!-----------------------------------------------------------------------
!         is:     index for "longitude"
!         js:     index for "latitude"
!
!         phalf       pressure at half levels in pascals
!                       [real, dimension(nlon,nlat,nlev+1)]
!-----------------------------------------------------------------------     
integer,         intent(in)          :: id, jd
real, intent(in), dimension(:,:,:)   ::  phalf

real, parameter          :: Hscale = 7500.0
! total # of number concentrations
integer :: ntot_num_trs
integer :: err_code    ! valid run?
     
!------------------------------------------------------------------------------------------
! --- Begin Code --- !

! Define model constant parameters
ntot_num_trs = 0
if(nqn > 0)    ntot_num_trs = ntot_num_trs + 1
if(nqni > 0)   ntot_num_trs = ntot_num_trs + 1

! setup vertical levels for CLUBB
nz_host_model = size( phalf, 3 )
nz_CLUBB  = nz_host_model
nsublevel = 0

!allocate momentum and thermodynamic heights
if ( .not. allocated( momentum_heights ) ) then
      allocate( momentum_heights(1:nz_CLUBB) )
else
      deallocate( momentum_heights )
      allocate( momentum_heights(1:nz_CLUBB) )
endif
    
if ( .not. allocated( thermodynamic_heights ) ) then
      allocate( thermodynamic_heights(1:nz_CLUBB) )
else
      deallocate( thermodynamic_heights )
      allocate( thermodynamic_heights(1:nz_CLUBB) )
endif

!using scale height (Hscale=7500) to get a generic momentum height
do iz_clubb =1, nz_CLUBB-1
   momentum_heights(iz_clubb)= Hscale*log(phalf(id, jd,nz_CLUBB)/phalf(id, jd,nz_CLUBB+1-iz_clubb))
enddo
iz_clubb=nz_CLUBB
momentum_heights(iz_clubb)=2.0*momentum_heights(iz_clubb-1)-momentum_heights(iz_clubb-2)

!using linear interpolation to get a  generic thermodynamic height
do iz_clubb =2, nz_CLUBB
   thermodynamic_heights(iz_clubb)=0.5*( momentum_heights(iz_clubb) &
                                        +momentum_heights(iz_clubb-1) )
enddo
iz_clubb=1
thermodynamic_heights(iz_clubb)=  2.0*thermodynamic_heights(iz_clubb+1) &
                                 -thermodynamic_heights(iz_clubb+2)

! Read in model parameter values
call read_parameters(1,"input.nml", params )   ! intent (out)

call set_clubb_debug_level( debug_level ) ! Intent(in)

call  setup_clubb_core & 
       ( nz_CLUBB, T0_in, ts_nudge_in, & ! In
         0,     ntot_num_trs,          & ! In
         sclr_tol,  ntot_num_trs,      &
         params,                       &  ! In
         l_host_applies_sfc_fluxes_in,                           &        ! intent(in)
         l_uv_nudge_in, saturation_formula_in,                   &        ! intent(in)
         I_sat_sphum_in,                                         &        ! intent(in)
         .true.,  2, avg_deltaz,  momentum_heights(1),   momentum_heights(nz_CLUBB),  &  ! intent(in) 
         momentum_heights,  thermodynamic_heights,  & ! In
         host_dx, host_dy,  momentum_heights(1),    & ! In
         cloud_frac_min_in ,                        & ! intent(in)
         err_code ) ! Out

! setup 3D variables to save higher-order moments in CLUBB 
call setup_CLUBB_3D_var(id, jd, nz_CLUBB )

return
end subroutine CLUBB_SETUP
!#######################################################################





!#######################################################################
SUBROUTINE host2CLUBB_full (var_host, var_CLUBB)
!=======================================================================
real, intent(in) ,  dimension(:)     ::  var_host
real, intent(out) , dimension(:)     ::  var_CLUBB

real, dimension( size(var_host)+2 )  :: var_host_tmp
integer ::      iz_host_modelp1
real    ::      weight_host
       
! nz_host_model = size(zhalf,3): total # of half levels
var_host_tmp(1)               = 2.0* var_host(  nz_host_model-1)- var_host( nz_host_model-2 )
var_host_tmp(nz_host_model+1) = 2.0* var_host(2) - var_host(1)

do iz_host_model = 2, nz_host_model
   var_host_tmp(iz_host_model) = var_host(nz_host_model +1 - iz_host_model)
enddo

iz_shift = nsublevel/2
do  iz_CLUBB=1, nz_clubb
    iz_host_model   = (iz_clubb - 1+ iz_shift)/(nsublevel+1)+1
    iz_host_modelp1 =  iz_host_model + 1
    weight_host     = real( mod(iz_clubb-1+ iz_shift, nsublevel+1) )/(nsublevel+1)
    var_CLUBB( iz_CLUBB) = var_host_tmp( iz_host_model ) * (1.0-weight_host)  &
                          +var_host_tmp( iz_host_modelp1 ) * weight_host
enddo
   
return
end subroutine host2CLUBB_full
!#######################################################################



!#######################################################################
SUBROUTINE host2CLUBB_half (var_host, var_CLUBB)
!=======================================================================
real, intent(in),  dimension(:)     ::  var_host
real, intent(out), dimension(:)     ::  var_CLUBB
integer ::      iz_host_modelm1
real::  weight_host, weight_hostm1

do iz_CLUBB=1, nz_CLUBB - 1
   iz_host_model   = nz_host_model - (iz_clubb-1)/(nsublevel+1)
   iz_host_modelm1 = iz_host_model - 1

   weight_hostm1 = real(mod(iz_clubb-1, nsublevel+1))/(nsublevel+1)
   weight_host   = 1.0 - weight_hostm1
   var_CLUBB( iz_CLUBB) = var_host( iz_host_model ) * weight_host &
                         +var_host( iz_host_modelm1 ) * weight_hostm1
enddo

iz_clubb             = nz_clubb
iz_host_model        = nz_host_model - (iz_clubb-1)/(nsublevel+1)
var_CLUBB( iz_CLUBB) = var_host( iz_host_model )

return
end subroutine host2CLUBB_half
!#######################################################################


!#######################################################################
SUBROUTINE CLUBB2host_full (var_CLUBB, var_host)
!=======================================================================
real, intent(in),  dimension(:)     ::  var_CLUBB
real, intent(out), dimension(:)     ::  var_host

iz_shift = nsublevel/2
do iz_CLUBB=2+iz_shift, nz_CLUBB, nsublevel+1
   iz_host_model = nz_host_model - (iz_CLUBB -1+iz_shift)/(nsublevel+1)  
   var_host( iz_host_model ) = var_CLUBB( iz_CLUBB)      
enddo

return
end subroutine CLUBB2host_full 
!#######################################################################


!#######################################################################
SUBROUTINE CLUBB2host_half (var_CLUBB, var_host)
!=======================================================================
real, intent(in),  dimension(:)     ::  var_CLUBB
real, intent(out), dimension(:)     ::  var_host
    
do iz_CLUBB=1, nz_CLUBB,  nsublevel+1    
   iz_host_model = nz_host_model  - (iz_CLUBB-1)/(nsublevel+1) 
   var_host( iz_host_model ) = var_CLUBB( iz_CLUBB)
enddo

return
end subroutine CLUBB2host_half
!#######################################################################




!#######################################################################
SUBROUTINE CLUBB_3D_2_1D( ix_CLUBB, iy_CLUBB )
implicit none
integer,  intent(in) ::  ix_CLUBB, iy_CLUBB
!=======================================================================
wp2( : )          = wp2_3D( ix_CLUBB, iy_CLUBB, : )
wp3( : )          = wp3_3D( ix_CLUBB, iy_CLUBB, : )
upwp( : )         = upwp_3D( ix_CLUBB, iy_CLUBB, : )
vpwp( : )         = vpwp_3D( ix_CLUBB, iy_CLUBB, : )
wprtp( : )        = wprtp_3D( ix_CLUBB, iy_CLUBB, : )
wpthlp( : )       = wpthlp_3D( ix_CLUBB, iy_CLUBB, : )
rtp2( : )         = rtp2_3D( ix_CLUBB, iy_CLUBB, : )
thlp2( : )        = thlp2_3D( ix_CLUBB, iy_CLUBB, : )
rtpthlp( : )      = rtpthlp_3D( ix_CLUBB, iy_CLUBB, : )
up2( : )          = up2_3D( ix_CLUBB, iy_CLUBB, : )
vp2( : )          = vp2_3D( ix_CLUBB, iy_CLUBB, : )

if ( sclr_dim > 0) then
     RH_crit( : , 1, 1:2 ) = RH_crit_clubb_3D( ix_CLUBB, iy_CLUBB, : , 1:2 )
endif

return
end subroutine CLUBB_3D_2_1D
!#######################################################################    


!#######################################################################
SUBROUTINE CLUBB_1D_2_3D( ix_CLUBB, iy_CLUBB )
implicit none
integer,  intent(in) ::  ix_CLUBB, iy_CLUBB
!=======================================================================
wp2_3D( ix_CLUBB, iy_CLUBB, : )         = wp2( : )
wp3_3D( ix_CLUBB, iy_CLUBB, : )         = wp3( : )
upwp_3D( ix_CLUBB, iy_CLUBB, : )        = upwp( : )
vpwp_3D( ix_CLUBB, iy_CLUBB, : )        = vpwp( : )
wprtp_3D( ix_CLUBB, iy_CLUBB, : )       = wprtp( : )
wpthlp_3D( ix_CLUBB, iy_CLUBB, : )      = wpthlp( : )
rtp2_3D( ix_CLUBB, iy_CLUBB, : )        = rtp2( : )
thlp2_3D( ix_CLUBB, iy_CLUBB, : )       = thlp2( : )
rtpthlp_3D( ix_CLUBB, iy_CLUBB, : )     = rtpthlp( : )
up2_3D( ix_CLUBB, iy_CLUBB, : )         = up2( : )
vp2_3D( ix_CLUBB, iy_CLUBB, : )         = vp2( : )
         
if ( sclr_dim > 0) then
     RH_crit_clubb_3D( ix_CLUBB, iy_CLUBB, : , 1:2 ) = RH_crit( : , 1, 1:2) 
endif

return
end subroutine CLUBB_1D_2_3D
!#######################################################################    



!#######################################################################
SUBROUTINE tndy_CLUBB ( dtmain, var_CLUBB,     &
                        var_host, vardt_host,  &
                        vardt_CLUBB )
!=======================================================================
real :: dtmain
real, intent(in),    dimension(:)    ::   var_CLUBB
real, intent(in),    dimension(:)    ::   var_host
real, intent(inout), dimension(:)    ::   vardt_host
real, intent(out),   dimension(:)    ::   vardt_CLUBB
 
iz_shift = nsublevel/2
do iz_CLUBB=2+iz_shift, nz_CLUBB, nsublevel+1
   iz_host_model = nz_host_model - (iz_CLUBB -1+iz_shift)/(nsublevel+1)  
   vardt_CLUBB( iz_host_model ) =                                             &
            ( var_CLUBB( iz_CLUBB ) - var_host( iz_host_model ) )/dtmain   &
           - vardt_host( iz_host_model)
 
   vardt_host( iz_host_model) = ( var_CLUBB( iz_CLUBB ) - var_host( iz_host_model ) )/dtmain
enddo
  
return
end subroutine tndy_CLUBB
!#######################################################################    



recursive function erff(x) RESULT(y)
! Error function from Numerical Recipes.
! erf(x) = 1 - erfc(x)

real dumerfc, x
real t, z, y

z = abs(x)
t = 1.0 / ( 1.0 + 0.5 * z )

dumerfc =     t * exp(-z * z - 1.26551223 + t *      &
          ( 1.00002368 + t * ( 0.37409196 + t *    &
          ( 0.09678418 + t * (-0.18628806 + t *    &
          ( 0.27886807 + t * (-1.13520398 + t *    &
          ( 1.48851587 + t * (-0.82215223 + t * 0.17087277 )))))))))

if ( x.lt.0.0 ) dumerfc = 2.0 - dumerfc
 
y = 1.0 - dumerfc

end function erff
!#######################################################################    



!#######################################################################    
 SUBROUTINE aer_act_clubb_BL_Gauss( aeromass_clubb, & ! Intent(in)
                                    Ndrop_max )       ! Intent(out)
! using Boucher an Lohmann empirical relationship between Nd and sulface mass concentration
! Nd = 10^(2.21+0.41 log(mso4))
!=======================================================================
use  aer_ccn_act_k_mod,   only: aer_ccn_act_k, aer_ccn_act_wpdf_k 

implicit none
real, intent(inout),  dimension(:, :)    ::    aeromass_clubb
real, intent(out),    dimension(:)       ::    Ndrop_max
real ::   drop
  
!=======================================================================
!  if      aeromass_clubb = 1 ug/m3,  drop = 162 /cm3
!  if      aeromass_clubb = 1 ug/m3,  drop = 314 /cm3
do iz_clubb = 2, nz_clubb
   drop = 2.21 + 0.41 * log10(  aeromass_clubb ( iz_clubb, 1 )*1.e12  )
   drop = exp( log(10.0) * drop )
   Ndrop_max ( iz_clubb ) = drop  &
          * ( pdf_params( iz_clubb)%mixt_frac * pdf_params( iz_clubb)%cloud_frac1  &
             +( 1.0-pdf_params( iz_clubb)%mixt_frac ) * pdf_params( iz_clubb)%cloud_frac2 )     
enddo  
return
end subroutine  aer_act_clubb_BL_Gauss
!#######################################################################    


SUBROUTINE aer_act_clubb_diffK_Gauss( aeromass_clubb, temp_clubb_act, &   ! Intent(in)
                                                      Ndrop_max )         ! Intent(out)
use aer_ccn_act_k_mod,   only: aer_ccn_act_k, aer_ccn_act_wpdf_k 

implicit none
real,  intent(in),      dimension(:)     ::    temp_clubb_act
real, intent(inout),    dimension(:, :)  ::    aeromass_clubb
real, intent(out),    dimension(:)       ::    Ndrop_max

real                ::  drop
integer             ::  Tym, ier
character(len=256)  ::  ermesg
!=======================================================================
Tym = size(aeromass_clubb, 2)

do iz_clubb = 2, nz_clubb
   call aer_ccn_act_wpdf_k(temp_clubb_act( iz_clubb), p_in_Pa( iz_clubb), &! intent (in)
                           wm_zt( iz_clubb),          Var_w,              &! intent (in)
                           aeromass_clubb( iz_clubb, : ), Tym,            &! intent (in)
                           drop,   ier,   ermesg )                         ! intent (out)

    Ndrop_max ( iz_clubb ) =   drop &
            * ( pdf_params(iz_clubb)%mixt_frac * pdf_params(iz_clubb)%cloud_frac1  &
               +( 1.0-pdf_params(iz_clubb)%mixt_frac)*pdf_params(iz_clubb)%cloud_frac2 )
enddo  
return
end subroutine  aer_act_clubb_diffK_Gauss
!###################################################################    






SUBROUTINE aer_act_clubb_quadrature_Gauss(Time_next, aeromass_clubb, temp_clubb_act, & ! Intent(in)
                                          Ndrop_max   )                                ! Intent(out)
!=======================================================================
use aer_ccn_act_k_mod,   only: aer_ccn_act_k, aer_ccn_act_wpdf_k 

implicit none
! ---> h1g, 2010-08-24, dumping Nact
type(time_type),         intent(in)      ::    Time_next
! <--- h1g, 2010-08-24

real,  intent(in),    dimension(:)       ::    temp_clubb_act
real, intent(inout),  dimension(:, :)    ::    aeromass_clubb
real, intent(out),    dimension(:)       ::    Ndrop_max

real                ::  drop
integer             ::  Tym, ier
character(len=256)  ::  ermesg

real               :: P1_updraft,   P2_updraft  ! probability of updraft
real, parameter    :: P_updraft_eps = 1.e-16    ! updraft probability threshold
real, parameter    :: wp2_eps = 0.0001          ! w variance threshold
 
!=======================================================================
Tym = size(aeromass_clubb, 2)

do iz_clubb = 2, nz_clubb

   if( pdf_params( iz_clubb)%varnce_w1 > wp2_eps) then
       P1_updraft = 0.5+0.5*erff( pdf_params( iz_clubb)%w1/sqrt(2.0*pdf_params( iz_clubb)%varnce_w1) )
       P1_updraft = P1_updraft * pdf_params( iz_clubb)%mixt_frac * pdf_params( iz_clubb)%cloud_frac1
   else
      if( pdf_params( iz_clubb)%w1 > 0.0) &
          P1_updraft = pdf_params( iz_clubb)%mixt_frac * pdf_params( iz_clubb)%cloud_frac1
   endif


   if( pdf_params( iz_clubb)%varnce_w2 > wp2_eps) then
       P2_updraft = 0.5+0.5*erff( pdf_params( iz_clubb)%w2/sqrt(2.0*pdf_params( iz_clubb)%varnce_w2) )
       P2_updraft = P2_updraft * ( 1.0-pdf_params( iz_clubb )%mixt_frac ) * pdf_params( iz_clubb)%cloud_frac2
   else
       if( pdf_params( iz_clubb)%w2 > 0.0) &
           P2_updraft = ( 1.0-pdf_params( iz_clubb )%mixt_frac ) * pdf_params( iz_clubb)%cloud_frac2
   endif

   if( P1_updraft + P2_updraft   > P_updraft_eps  ) then
       P1_updraft =  P1_updraft / (  P1_updraft + P2_updraft  )
       P2_updraft =  P2_updraft / (  P1_updraft + P2_updraft  )
   else
       P1_updraft = 0.0
       P2_updraft = 0.0
   endif
 
   call aer_ccn_act_wpdf_k( temp_clubb_act( iz_clubb),  p_in_Pa( iz_clubb),              &! intent (in)
                            pdf_params( iz_clubb)%w1,   pdf_params( iz_clubb)%varnce_w1, &! intent (in)
                            aeromass_clubb( iz_clubb, : ), Tym,                          &! intent (in)
                            drop,   ier,   ermesg )           ! intent (out)
    
    Ndrop_max ( iz_clubb ) = drop * P1_updraft

    call aer_ccn_act_wpdf_k( temp_clubb_act( iz_clubb),  p_in_Pa( iz_clubb),               &! intent (in)
                             pdf_params( iz_clubb)%w2,   pdf_params( iz_clubb)%varnce_w2,  &! intent (in)
                             aeromass_clubb( iz_clubb, : ), Tym,                           &! intent (in)
                             drop,   ier,   ermesg)           ! intent (out)
! in-cloud activated droplet concentration
    Ndrop_max ( iz_clubb ) = Ndrop_max( iz_clubb ) + drop * P2_updraft

! get the layer-averaged activated droplet concentration (/cm3)
    Ndrop_max ( iz_clubb ) = Ndrop_max ( iz_clubb ) *  &
                 (      pdf_params( iz_clubb)%mixt_frac  * pdf_params( iz_clubb)%cloud_frac1 + &
                   (1.- pdf_params( iz_clubb)%mixt_frac) * pdf_params( iz_clubb)%cloud_frac2 )

enddo
return
end subroutine aer_act_clubb_quadrature_Gauss
!#######################################################################    



subroutine warm_ice_nucl_mns_clubb( aeromass_clubb,      &    !  Intent(in)  
                                    temp_clubb_act,      &    !  Intent(in)  
                                    Ndrop_max,           &    !  Intent(out)  
                                    Ncrystal_max,        &    !  Intent(out)
                                    RH_crit_clubb )           !  Intent(out)

use  aer_ccn_act_k_mod,   only: aer_ccn_act_k,  aer_ccn_act_wpdf_k 
use  sat_vapor_pres_mod,  only: compute_qs
USE  polysvp_mod,         ONLY: polysvp_l,  polysvp_i
USE  ice_nucl_mod,        ONLY: ice_nucl_wpdf


implicit none          
real,  intent(in),    dimension(:)       ::    temp_clubb_act
real,  intent(inOUT), dimension(:, :)    ::    aeromass_clubb
     
real, intent(out),    dimension(:)       ::    Ndrop_max     ! domain warm cloud droplet number concentration
real, intent(out),    dimension(:)       ::    Ncrystal_max  ! domain ice crystal number concentration
real, intent(INout),  dimension( gr%nz, 1:min(1,sclr_dim), 2 )  :: RH_crit_clubb       ! critical relative humidity for nucleation
! -------------------------------------------------------------------------------------------------------------
!local variables for droplet nucleation
real               :: drop
integer            :: Tym, ier
character(len=256) :: ermesg

real               :: P1_updraft,   P2_updraft ! probability of updraft
real, parameter    :: P_updraft_eps = 1.e-16   ! updraft probability threshold
real, parameter    :: wp2_eps = 0.0001         ! w variance threshold

! -------------------------------------------------------------------------------------------------------------
!local variables for ice nucleation
logical         ::      do_ice_nucl_ss_wpdf
real            ::      d_sulf
real            ::      d_bc
logical         ::      do_het
logical         ::      use_dust_instead_of_bc
logical         ::      limit_immersion_frz
logical         ::      limit_rhil

real            ::      dust_frac
real            ::      rh_crit_het
real            ::      dust_frac_min
real            ::      dust_frac_max
Real            ::      dust_surf
integer         ::      dust_opt
real            ::      rh_dust_max

real            ::      rh_crit_1d
real            ::      rh_crit_min_1d

real            ::      ni_sulf, ni_dust, ni_bc
      
integer, parameter  :: n_totmass_in  = 4
integer, parameter  :: n_imass_in     = 12

integer :: n_totmass, n_imass

REAL, DIMENSION(  n_totmass_in  ) :: totalmass !r1 , mass_ratio
REAL, DIMENSION( n_imass_in ) :: imass
real, parameter :: d378 = 0.378

real  ::   qs
real  ::   eslt,  qs_d, qvsl,  esit,  qvsi,  qvt,  u_i,  u_l
real  ::    Ncrystal1, Ncrystal2
   
REAL, DIMENSION( nz_clubb ) ::  hom
!===========================================================

n_totmass  = n_totmass_in
n_imass     = n_imass_in

do_ice_nucl_ss_wpdf    = do_ice_nucl_ss_wpdf_in
d_sulf                 = d_sulf_in
d_bc                   = d_bc_in
do_het                 = do_het_in
use_dust_instead_of_bc = use_dust_instead_of_bc_in
limit_immersion_frz    = limit_immersion_frz_in
limit_rhil             = limit_rhil_in

dust_frac              = dust_frac_in
rh_crit_het            = rh_crit_het_in

dust_frac_min          = dust_frac_min_in
dust_frac_max          = dust_frac_max_in
dust_surf              = dust_surf_in
dust_opt               = dust_opt_in
rh_dust_max            = rh_dust_max_in


Tym = size(aeromass_clubb, 2)

do iz_clubb = 2, nz_clubb

  P1_updraft = 0.0
  P2_updraft = 0.0
  if( pdf_params( iz_clubb)%varnce_w1 > wp2_eps) then
    P1_updraft = 0.5+0.5*erff( pdf_params( iz_clubb)%w1/sqrt(2.0*pdf_params( iz_clubb)%varnce_w1) )
    P1_updraft = P1_updraft * pdf_params( iz_clubb)%mixt_frac * pdf_params( iz_clubb)%cloud_frac1
  else
    if( pdf_params( iz_clubb)%w1 > 0.0) &
                 P1_updraft = pdf_params( iz_clubb)%mixt_frac * pdf_params( iz_clubb)%cloud_frac1
  endif


  if( pdf_params( iz_clubb)%varnce_w2 > wp2_eps) then
    P2_updraft = 0.5+0.5*erff( pdf_params( iz_clubb)%w2/sqrt(2.0*pdf_params( iz_clubb)%varnce_w2) )
    P2_updraft = P2_updraft * ( 1.0-pdf_params( iz_clubb )%mixt_frac ) * pdf_params( iz_clubb)%cloud_frac2
  else
    if( pdf_params( iz_clubb)%w2 > 0.0) &
                 P2_updraft = ( 1.0-pdf_params( iz_clubb )%mixt_frac ) * pdf_params( iz_clubb)%cloud_frac2
  endif

  if( P1_updraft + P2_updraft > P_updraft_eps ) then
    P1_updraft =  P1_updraft / (  P1_updraft + P2_updraft  )
    P2_updraft =  P2_updraft / (  P1_updraft + P2_updraft  )
  endif

! for the first Gaussian distribution portion of the warm cloud nucleation 
  call aer_ccn_act_wpdf_k ( temp_clubb_act( iz_clubb), p_in_Pa( iz_clubb), &  !  intent(in)
                            pdf_params( iz_clubb)%w1,                      &  !  intent(in)
                            pdf_params( iz_clubb)%varnce_w1,               &  !  intent(in)
                            aeromass_clubb( iz_clubb, : ), Tym,            &  !  intent(in)
                            drop,   ier,   ermesg )                           !  intent(out)

  Ndrop_max ( iz_clubb ) = drop * P1_updraft


! for the first Gaussian distribution portion of the ice cloud nucleation 
  imass(:)     = aeromass_clubb(iz_clubb,1)
  totalmass(:) = aeromass_clubb(iz_clubb,:)

  eslt         =  polysvp_l( temp_clubb(iz_clubb) )  !satuarated vapor pressure with respect to liquid
  qs_d         =  p_in_Pa(iz_clubb) - d378*eslt
  qs_d         =  max(qs_d,eslt)
  qvsl         = 0.622 *eslt / qs_d

  esit         = polysvp_i( temp_clubb(iz_clubb) )  !satuarated vapor pressure with respect to ice
  qs_d         = p_in_Pa(iz_clubb) - d378*esit
  qs_d         = max(qs_d,esit)
  qvsi         = 0.622 *esit / qs_d

!cms 4/21/2009 changed nothing
  IF ( rh_act_opt .EQ. 1 ) THEN
     qvt = rtm(iz_clubb) - rcm(iz_clubb)
  ELSE
!cms 4/23/2009  changed 
!environmental qv
     IF ( cloud_frac(iz_clubb) .LT. cf_thresh_nucl ) THEN
          call compute_qs( temp_clubb(iz_clubb),  p_in_Pa(iz_clubb), qs)
          qvt =  (  rtm(iz_clubb) - rcm(iz_clubb) - cloud_frac(iz_clubb) * qs ) / (1. - cloud_frac(iz_clubb) )
     ELSE
          qvt =  rtm(iz_clubb) - rcm(iz_clubb)
     ENDIF
  endif

  u_i =  qvt/qvsi
  u_l =  qvt/qvsl

  rh_crit_1d     = 1.
  rh_crit_min_1d = 1.

  hom     = 0.0
  ni_sulf = 0.
  ni_dust = 0.
  ni_bc = 0.

  call ice_nucl_wpdf( 1, 1, iz_clubb, &! intent (in) 
                      do_ice_nucl_ss_wpdf,   temp_clubb_act(iz_clubb),  u_i    , u_l,  &! intent (in) 
                      pdf_params( iz_clubb )%w1,         &! intent (in) 
                      pdf_params( iz_clubb )%varnce_w1,  &! intent (in)
                      thermodynamic_heights( iz_clubb ), &! intent (in)
                      totalmass, imass,                  &! intent (in)
                      n_totmass, n_imass,                &! intent (in)
                      Ncrystal1,                         &! intent (inout)
                      drop,                              &! intent (inout)
                      hom(iz_clubb),                     &! intent (inout)
                      d_sulf, d_bc, do_het,              &! intent (in)
                      use_dust_instead_of_bc, limit_immersion_frz, limit_rhil, &! intent (in)
                      dust_frac, rh_crit_het,            &! intent (in)
                      rh_crit_1d, rh_crit_min_1d,        &! intent (inout)
                      dust_frac_min, dust_frac_max,      &! intent (in)
                      dust_surf, dust_opt, rh_dust_max,  &! intent (in)
                      ni_sulf, ni_dust, ni_bc )           ! intent (out)

  if( sclr_dim > 0) then
    if( temp_clubb( iz_clubb ) .LT. 250. ) THEN
        RH_crit_clubb(  iz_clubb, 1, 1 ) = MAX(rh_crit_1d, 1.)
    ELSE
        RH_crit_clubb(  iz_clubb, 1, 1 ) = 1.
    END IF
  endif



! for the second Gaussian distribution portion of the warm cloud nucleation 
  call aer_ccn_act_wpdf_k ( temp_clubb_act( iz_clubb), p_in_Pa( iz_clubb), &   ! intent (in)
                            pdf_params( iz_clubb)%w2,                      &   ! intent (in) 
                            pdf_params( iz_clubb)%varnce_w2,               &   ! intent (in)
                            aeromass_clubb( iz_clubb, : ), Tym,            &   ! intent (in)
                            drop, ier, ermesg )                                ! intent (out)

! in-cloud activated droplet concentration
  Ndrop_max ( iz_clubb ) = Ndrop_max ( iz_clubb ) + drop *  P2_updraft

! get the layer-averaged activated droplet concentration (/cm3)
  Ndrop_max ( iz_clubb ) = Ndrop_max ( iz_clubb ) *  &
                   (    pdf_params( iz_clubb)%mixt_frac * pdf_params( iz_clubb)%cloud_frac1 + &
                   (1.- pdf_params( iz_clubb)%mixt_frac)* pdf_params( iz_clubb)%cloud_frac2 )


! for the second Gaussian distribution portion of the ice cloud nucleation 
  rh_crit_1d     = 1.
  rh_crit_min_1d = 1.

  call ice_nucl_wpdf( 1, 1, iz_clubb, &! intent (in) 
                      do_ice_nucl_ss_wpdf,   temp_clubb_act(iz_clubb),  u_i    , u_l,  &! intent (in) 
                      pdf_params( iz_clubb)%w2,  &! intent (in) 
                      pdf_params( iz_clubb)%varnce_w2, &! intent (in)
                      thermodynamic_heights( iz_clubb), &! intent (in)
                      totalmass, imass,  &! intent (in)
                      n_totmass, n_imass,&! intent (in)
                      Ncrystal2,&! intent (inout)
                      drop,&! intent (inout)
                      hom(iz_clubb), &! intent (inout)
                      d_sulf, d_bc, do_het,  &! intent (in)
                      use_dust_instead_of_bc, limit_immersion_frz, limit_rhil, &! intent (in)
                      dust_frac, rh_crit_het, &! intent (in)
                      rh_crit_1d, rh_crit_min_1d,&! intent (inout)
                      dust_frac_min, dust_frac_max,&! intent (in)
                      dust_surf, dust_opt, rh_dust_max, &! intent (in)
                      ni_sulf, ni_dust, ni_bc )  ! intent (out)

  if( sclr_dim > 0) then
    if( temp_clubb( iz_clubb ) .LT. 250. ) THEN
        RH_crit_clubb(  iz_clubb, 1, 2 ) = MAX(rh_crit_1d, 1.)
    ELSE
        RH_crit_clubb(  iz_clubb, 1, 2 ) = 1.
    END IF
  endif

  Ncrystal_max(iz_clubb) = Ncrystal1 * pdf_params( iz_clubb)%mixt_frac * pdf_params( iz_clubb) %cloud_frac1&
                                                          + Ncrystal2 * (1.0 - pdf_params( iz_clubb)%mixt_frac) * pdf_params( iz_clubb)%cloud_frac2
enddo

return
end subroutine  warm_ice_nucl_mns_clubb
!#######################################################################    




!#####################################################################
subroutine get_aer_mass_host ( is, js, Time, phalf, airdens, T, &
                            concen_dust_sub, totalmass1, Aerosol, mask)

integer, intent (in)                           :: is,  js
type(time_type), intent (in)                   :: Time
real, dimension(:,:,:), intent(in )            :: phalf,  airdens,  T 
real, dimension(:,:,:), intent(out)            :: concen_dust_sub
real, dimension(:,:,:,:), intent(out)          :: totalmass1
type(aerosol_type), intent (in), optional      :: Aerosol  
real, intent (in), optional, dimension(:,:,:)  :: mask

real, dimension(size(T,1),size(T,2),size(T,3)) :: pthickness
real, dimension(size(T,1),size(T,2),size(T,3)) :: concen,                      &
                                                  concen_all_sub,              &
                                                  concen_ss_sub, concen_ss_sup,&
                                                  concen_om, concen_na, concen_an
integer  :: i,  j,  k,  na , s
integer  :: idim,   jdim,   kdim
logical  :: used

idim = size(T,1)
jdim = size(T,2)
kdim = size(T,3)

concen_dust_sub(:,:,:) = 0.
totalmass1(:,:,:,:)    = 0.

if (id_sulfate > 0) then
    concen_an(:,:,:) = 0.
    concen_na(:,:,:) = 0.
    concen(:,:,:) = 0.
endif

if (use_online_aerosol) then
    concen_ss_sub(:,:,:) = 0.
    concen_ss_sup(:,:,:) = 0.
    concen_all_sub(:,:,:) = 0.
endif

do k = 1, kdim
 do j = 1, jdim
  do i = 1, idim
   if (phalf(i,j,k) < 1.0) then
       pthickness(i,j,k) = (phalf(i,j,k+1) - phalf(i,j,k))/&
                           grav/airdens(i,j,k)
   else
       pthickness(i,j,k) = log(phalf(i,j,k+1)/ &
                            phalf(i,j,k))*8.314*T(i,j,k)/(9.8*0.02888)
   end if
  end do
 end do
end do
 
if (present (Aerosol)) then
    if (do_liq_num) then
        if (use_online_aerosol) then
            do na = 1,size(Aerosol%aerosol,4)
               if (trim(Aerosol%aerosol_names(na)) == 'so4' .or. &
                   trim(Aerosol%aerosol_names(na)) == 'so4_anthro' .or.&
                   trim(Aerosol%aerosol_names(na)) == 'so4_natural')  &
                                                                 then
                do k = 1,kdim
                 do j = 1,jdim
                   do i = 1,idim
                     totalmass1(i,j,k,1) = totalmass1(i,j,k,1) + &
                                           Aerosol%aerosol(i,j,k,na)
                   end do
                 end do
                end do
               else if(trim(Aerosol%aerosol_names(na)) == 'omphilic' .or.&
                     trim(Aerosol%aerosol_names(na)) == 'omphobic') &
                                                                 then
                do k = 1,kdim
                 do j = 1,jdim
                   do i = 1,idim
                     totalmass1(i,j,k,4) = totalmass1(i,j,k,4) +  &
                                           Aerosol%aerosol(i,j,k,na)
                   end do
                 end do
                end do
               else if(trim(Aerosol%aerosol_names(na)) == 'seasalt1' .or.&
                     trim(Aerosol%aerosol_names(na)) == 'seasalt2') &
                                                                   then
                do k = 1,kdim
                 do j = 1,jdim
                   do i = 1,idim
                     concen_ss_sub(i,j,k) = concen_ss_sub(i,j,k) +  &
                                            Aerosol%aerosol(i,j,k,na)
                   end do
                 end do
                end do
               else if(trim(Aerosol%aerosol_names(na)) == 'seasalt3' .or.&
                     trim(Aerosol%aerosol_names(na)) == 'seasalt4' .or.&
                     trim(Aerosol%aerosol_names(na)) == 'seasalt5')  &
                                                                  then
                do k = 1,kdim
                 do j = 1,jdim
                   do i = 1,idim
                     concen_ss_sup(i,j,k) = concen_ss_sup(i,j,k) +  &
                                            Aerosol%aerosol(i,j,k,na)
                   end do
                 end do
                end do
               else if(trim(Aerosol%aerosol_names(na)) == 'bcphilic' .or.&
                     trim(Aerosol%aerosol_names(na)) == 'bcphobic' .or.&
                     trim(Aerosol%aerosol_names(na)) == 'dust1' .or.&
                     trim(Aerosol%aerosol_names(na)) == 'dust2' .or.&
                     trim(Aerosol%aerosol_names(na)) == 'dust3')  &
                                                                  then
                do k = 1,kdim
                 do j = 1,jdim
                   do i = 1,idim
                     concen_all_sub(i,j,k) = concen_all_sub(i,j,k) +  &
                                             Aerosol%aerosol(i,j,k,na)
                   end do
                 end do
                end do
               endif
               if (do_dust_berg) then
                if (trim(Aerosol%aerosol_names(na)) == 'dust1' .or. &
                    trim(Aerosol%aerosol_names(na)) == 'dust2' .or. &
                    trim( Aerosol%aerosol_names(na)) == 'dust3') then
                do k = 1,kdim
                  do j = 1,jdim
                    do i = 1,idim
                       concen_dust_sub(i,j,k) =    &
                                           concen_dust_sub(i,j,k) +   &
                                              Aerosol%aerosol(i,j,k,na)
                    end do
                  end do
                end do
               endif
             endif
           end do
!        endif
!      endif
          
!       if (do_liq_num) then
!         if (use_online_aerosol) then
           do k = 1,kdim
             do j = 1,jdim
               do i = 1,idim
                 totalmass1(i,j,k,3) = concen_ss_sub(i,j,k)
                 totalmass1(i,j,k,2) = concen_all_sub(i,j,k) + &
                                       totalmass1(i,j,k,4) + &
                                       concen_ss_sub(i,j,k)
               end do
             end do
           end do
           if (use_sub_seasalt) then
           else
             do k = 1,kdim
               do j = 1,jdim
                 do i = 1,idim
                   totalmass1(i,j,k,3) = concen_ss_sub(i,j,k) +  &
                                                  concen_ss_sup(i,j,k)
                 end do
               end do
             end do
           endif

           if (id_sulfate > 0) then
             do k = 1,kdim
               do j = 1,jdim
                 do i = 1,idim
                   concen(i,j,k) = 0.7273*totalmass1(i,j,k,1)/  &
                                       pthickness(i,j,k)*1.0e9
                 end do
               end do
             end do
           endif
         
           do k = 1,kdim
             do j = 1,jdim
               do i = 1,idim
                 concen_ss_sub(i,j,k) = concen_ss_sub(i,j,k)/  &
                                              pthickness(i,j,k)*1.0e9
                 concen_ss_sup(i,j,k) = concen_ss_sup(i,j,k)/  &
                                              pthickness(i,j,k)*1.0e9
               end do
             end do
           end do

         else  ! (use_online_aerosol)
           if (do_dust_berg) then
!     YMice submicron dust (NO. 14 to NO. 18)
             do s = 14,18
               do k = 1,kdim
                 do j = 1,jdim
                   do i = 1,idim
                     concen_dust_sub(i,j,k) = concen_dust_sub(i,j,k)+ &
                                              Aerosol%aerosol(i,j,k,s)
                   end do
                 end do
               end do
             end do
           endif

           if (id_sulfate > 0) then
             do k = 1,kdim
               do j = 1,jdim
                 do i = 1,idim
!     anthro. and natural sulfate concentration (ug so4/m3)
                   concen_an(i,j,k) = 0.7273*Aerosol%aerosol(i,j,k,1)/&
                                                pthickness(i,j,k)*1.0e9
                   concen_na(i,j,k) = 0.7273*Aerosol%aerosol(i,j,k,2)/&
                                                pthickness(i,j,k)*1.0e9
                   concen(i,j,k) = concen_an(i,j,k) + concen_na(i,j,k)
                 end do
               end do
             end do
           endif

           do k = 1,kdim
             do j = 1,jdim
               do i = 1,idim
!offline
! NO. 1 natural Sulfate; NO. 2 anthro. sulfate; NO. 3 Sea Salt; NO. 4 Organics
                 totalmass1(i,j,k,1) = Aerosol%aerosol(i,j,k,2)
                 totalmass1(i,j,k,2) = Aerosol%aerosol(i,j,k,1)
                 totalmass1(i,j,k,3) = sea_salt_scale*  &
                                       Aerosol%aerosol(i,j,k,5)
                 totalmass1(i,j,k,4) = om_to_oc*  &
                                       Aerosol%aerosol(i,j,k,3)
               end do
             end do
           end do
         endif ! (use_online_aerosol)

         do na = 1, 4
           do k = 1,kdim
             do j = 1,jdim
               do i = 1,idim
                 totalmass1(i,j,k,na) = totalmass1(i,j,k,na)/  &
                                        pthickness(i,j,k)*1.0e9*1.0e-12
               end do
             end do
           end do
         end do
         if (do_dust_berg) then
           do k = 1,kdim
             do j = 1,jdim
               do i = 1,idim
! submicron dust concentration (ug/m3) (NO. 2 to NO. 4)
                 concen_dust_sub(i,j,k) = concen_dust_sub(i,j,k)/ &
                                              pthickness(i,j,k)*1.0e9 
               end do
             end do
           end do
         endif

         if (id_sulfate > 0) then
           used = send_data ( id_sulfate, concen, Time, is, js, 1,&
                              rmask=mask )
         end if

         if (use_online_aerosol) then
           if (id_seasalt_sub > 0) then
             used = send_data (id_seasalt_sub, concen_ss_sub, Time, &
                               is, js, 1, rmask=mask )
           endif
   
           if (id_seasalt_sup > 0) then
             used = send_data (id_seasalt_sup, concen_ss_sup, Time, &
                               is, js, 1, rmask=mask )
           endif
         endif

         if (id_om > 0) then
           do k = 1,kdim
             do j = 1,jdim
               do i = 1,idim
                 concen_om(i,j,k) = totalmass1(i,j,k,2)*1.0e12
               end do
             end do
           end do
           used = send_data (id_om, concen_om, Time, is, js, 1,&
                             rmask=mask )
         endif
       endif  ! (do_liq_num)
     endif ! (Present(Aerosol))

!----------------------------------------------------------------------
end subroutine   get_aer_mass_host









SUBROUTINE get_aer_mass(  ix_host_model, iy_host_model,  Aerosol  , aeromass_clubb)

integer, intent (in)                         :: ix_host_model, iy_host_model
real, dimension( : , : ), intent(out)        :: aeromass_clubb
type(aerosol_type), intent (in), optional    :: Aerosol  

!local variables
integer     :: na, s
real, dimension( size(aeromass_clubb, 1))    :: concen_ss_sub_clubb
real, dimension( size(aeromass_clubb, 1))    :: concen_ss_sup_clubb
real, dimension( size(aeromass_clubb, 1))    :: concen_all_sub_clubb
real, dimension( size(aeromass_clubb, 1))    :: concen_dust_sub_clubb
real, dimension( size(aeromass_clubb, 1))    :: concen_an_clubb
real, dimension( size(aeromass_clubb, 1))    :: concen_na_clubb  
real, dimension( size(aeromass_clubb, 1))    :: concen_clubb
real, dimension( size(aeromass_clubb, 1))    :: concen_om_clubb

integer :: id_sulfate, id_om
!=======================================================================
aeromass_clubb        = 0.0
concen_dust_sub_clubb = 0.0

if (id_sulfate > 0) then
    concen_an_clubb   = 0.0
    concen_na_clubb   = 0.0
    concen_clubb      = 0.0  
endif

if (use_online_aerosol) then
    concen_ss_sub_clubb   = 0.0
    concen_ss_sup_clubb   = 0.0
    concen_all_sub_clubb  = 0.0
endif

if (present (Aerosol)) then
    if( nqn > 0 ) then
        if (use_online_aerosol) then
           do na = 1,size(Aerosol%aerosol,4)               
             if (trim(Aerosol%aerosol_names(na)) == 'so4' .or. &
                 trim(Aerosol%aerosol_names(na)) == 'so4_anthro' .or.&
                 trim(Aerosol%aerosol_names(na)) == 'so4_natural')  &
                                                                 then
               do iz_clubb = 2,  nz_clubb
                   iz_host_model = nz_host_model +1 - iz_CLUBB
                   aeromass_clubb( iz_clubb, 1 ) = aeromass_clubb( iz_clubb, 1 ) + &
                   Aerosol%aerosol( ix_host_model, iy_host_model, iz_host_model, na)
               end do

             else if(trim(Aerosol%aerosol_names(na)) == 'omphilic' .or. &
                     trim(Aerosol%aerosol_names(na)) == 'omphobic') &
                                                                 then
               do iz_clubb = 2,  nz_clubb
                   iz_host_model = nz_host_model +1 - iz_CLUBB
                   aeromass_clubb( iz_clubb, 4 ) = aeromass_clubb( iz_clubb, 4 ) + &
                   Aerosol%aerosol( ix_host_model, iy_host_model, iz_host_model, na)
               end do

             else if(trim(Aerosol%aerosol_names(na)) == 'seasalt1' .or.&
                     trim(Aerosol%aerosol_names(na)) == 'seasalt2') &
                                                                   then
 
               do iz_clubb = 2,  nz_clubb
                   iz_host_model = nz_host_model +1 - iz_CLUBB
                   concen_ss_sub_clubb( iz_clubb) = concen_ss_sub_clubb( iz_clubb) + &
                   Aerosol%aerosol( ix_host_model, iy_host_model, iz_host_model, na)
               end do

             else if(trim(Aerosol%aerosol_names(na)) == 'seasalt3' .or.&
                     trim(Aerosol%aerosol_names(na)) == 'seasalt4' .or.&
                     trim(Aerosol%aerosol_names(na)) == 'seasalt5')  &
                                                                  then
 
               do iz_clubb = 2,  nz_clubb
                   iz_host_model = nz_host_model +1 - iz_CLUBB
                   concen_ss_sup_clubb( iz_clubb) = concen_ss_sup_clubb( iz_clubb) + &
                   Aerosol%aerosol( ix_host_model, iy_host_model, iz_host_model, na)
               end do

            else if(trim(Aerosol%aerosol_names(na)) == 'bcphilic' .or.&
                     trim(Aerosol%aerosol_names(na)) == 'bcphobic' .or.&
                     trim(Aerosol%aerosol_names(na)) == 'dust1' .or.&
                     trim(Aerosol%aerosol_names(na)) == 'dust2' .or.&
                     trim(Aerosol%aerosol_names(na)) == 'dust3')  &
                                                                  then

               do iz_clubb = 2,  nz_clubb
                   iz_host_model = nz_host_model +1 - iz_CLUBB
                   concen_all_sub_clubb( iz_clubb) = concen_all_sub_clubb( iz_clubb) + &
                   Aerosol%aerosol( ix_host_model, iy_host_model, iz_host_model, na)
               end do
             endif  ! trim (so4, dust3)

             if (do_dust_berg) then
               if (trim(Aerosol%aerosol_names(na)) == 'dust1' .or. &
                   trim(Aerosol%aerosol_names(na)) == 'dust2' .or. &
                   trim( Aerosol%aerosol_names(na)) == 'dust3') then

                 do iz_clubb = 2,  nz_clubb
                   iz_host_model = nz_host_model +1 - iz_CLUBB
                   concen_dust_sub_clubb( iz_clubb) = concen_dust_sub_clubb( iz_clubb) + &
                   Aerosol%aerosol( ix_host_model, iy_host_model, iz_host_model, na)
                 end do
               endif !trim (bcphilic, ..., dust3)
             endif   ! do_dust_berg
            enddo   ! na

           aeromass_clubb( : , 3 ) = concen_ss_sub_clubb( : )
           aeromass_clubb( : , 2 ) = concen_all_sub_clubb( : ) + &
                                       aeromass_clubb( : , 4) + &
                                       concen_ss_sub_clubb( : )

           if (use_sub_seasalt) then
           else
            aeromass_clubb( : , 3 ) = concen_ss_sub_clubb( : ) + &
                                                    concen_ss_sup_clubb( : )
            endif
 
            if (id_sulfate > 0) then
              concen_clubb( : ) = 0.7273 * aeromass_clubb( : , 1)/&
                                       gr%invrs_dzt( iz_clubb )*1.0e9
            endif
           
	   do iz_clubb = 2, nz_clubb
             concen_ss_sub_clubb ( iz_clubb ) = concen_ss_sub_clubb ( iz_clubb)/&
                                       gr%invrs_dzt( iz_clubb )*1.0e9
             concen_ss_sup_clubb ( iz_clubb ) = concen_ss_sup_clubb ( iz_clubb)/&
                                       gr%invrs_dzt( iz_clubb )*1.0e9
           enddo
                  !(the above use_online_aerosol)

         else  ! (the following use_offline_aerosol)
           if (do_dust_berg) then
!     YMice submicron dust (NO. 14 to NO. 18)
             do s = 14,18
               do iz_clubb = 2,  nz_clubb
                     iz_host_model = nz_host_model +1 - iz_CLUBB
                     concen_dust_sub_clubb (  iz_clubb ) = concen_dust_sub_clubb (  iz_clubb ) + &
                                              Aerosol%aerosol( ix_host_model, iy_host_model, iz_host_model, s)
                end do
             end do
           endif


           if (id_sulfate > 0) then
             do iz_clubb = 2,  nz_clubb
                   iz_host_model = nz_host_model +1 - iz_CLUBB
 !    anthro. and natural sulfate concentration (ug so4/m3)
                   concen_an_clubb(  iz_clubb ) = 0.7273* &
                   Aerosol%aerosol( ix_host_model, iy_host_model, iz_host_model, 1)/&
                                                gr%invrs_dzt( iz_clubb )*1.0e9
                   concen_na_clubb ( iz_clubb ) = 0.7273* &
                   Aerosol%aerosol( ix_host_model, iy_host_model, iz_host_model, 2)/&
                                                gr%invrs_dzt( iz_clubb )*1.0e9
                   concen_clubb(  iz_clubb ) =  concen_an_clubb(  iz_clubb ) + concen_na_clubb ( iz_clubb )
             end do
          endif


            do iz_clubb = 2,  nz_clubb
!off- line
! NO. 1 natural Sulfate; NO. 2 anthro. sulfate; NO. 3 Sea Salt; NO. 4 Or organics
               iz_host_model = nz_host_model +1 - iz_CLUBB
               aeromass_clubb( iz_clubb , 1) = Aerosol%aerosol( ix_host_model, iy_host_model, iz_host_model,2)
               aeromass_clubb( iz_clubb , 2) = Aerosol%aerosol( ix_host_model, iy_host_model, iz_host_model,1)
               aeromass_clubb( iz_clubb , 3) = sea_salt_scale*  &
                                                        Aerosol%aerosol( ix_host_model, iy_host_model, iz_host_model,5)
               aeromass_clubb( iz_clubb , 4) = om_to_oc*  &
                                                        Aerosol%aerosol( ix_host_model, iy_host_model, iz_host_model,3)
             end do

           endif ! (use_online_aerosol)

           do iz_clubb = 2,  nz_clubb
              aeromass_clubb( iz_clubb , : ) = aeromass_clubb( iz_clubb , : )/&
                                        gr%invrs_dzt( iz_clubb )*1.0e9*1.0e-12
           end do

           if (do_dust_berg) then
             do iz_clubb = 2,  nz_clubb
               concen_dust_sub_clubb(  iz_clubb ) = concen_dust_sub_clubb(  iz_clubb )/ &
                                              gr%invrs_dzt( iz_clubb )*1.0e9 
             end do
           endif

           if (id_om > 0) then
             concen_om_clubb( : ) = aeromass_clubb( : , 2)*1.0e12
           endif

  
        endif  ! (nqn > 0)
endif ! (Present(Aerosol))

return
end subroutine get_aer_mass
 !#######################################################################    
    





!#######################################################################
SUBROUTINE CLUBB_END

use clubb_core, only: cleanup_clubb_core  ! Procedure

implicit none
character(len=64)     :: fname='RESTART/CLUBB.res.nc'
character(len=64)     :: sclr_prefix,  sclr_txt
integer               ::  isclr_dim
   
!---------------------------------------------------------------------
!    verify that the module is initialized.
!---------------------------------------------------------------------
if ( .not. module_is_initialized) then
     call error_mesg('CLUBB_driver_mod','module has not been initialized', FATAL)
     return
else
     module_is_initialized = .false.
endif

! If stats is active on current PE
if ( istats.ge.0 .and. jstats.ge.0 ) then
     l_stats = .true.
     call stats_finalize()
endif

if( do_netcdf_restart) then
    if (mpp_pe() == mpp_root_pe()) then
       call mpp_error ('CLUBB_mod', 'Writing netCDF formatted restart file: RESTART/CLUBB.res.nc', NOTE)
    endif
    call write_data (fname, 'vers', restart_versions(size(restart_versions(:))) , no_domain=.true.)	   
    call write_data (fname, 'wp2_3D',        wp2_3D)
    call write_data (fname, 'wp3_3D',        wp3_3D)
    call write_data (fname, 'upwp_3D',       upwp_3D)
    call write_data (fname, 'vpwp_3D',       vpwp_3D)
    call write_data (fname, 'wprtp_3D',      wprtp_3D)
    call write_data (fname, 'wpthlp_3D',     wpthlp_3D)
    call write_data (fname, 'rtp2_3D',       rtp2_3D)
    call write_data (fname, 'thlp2_3D',      thlp2_3D)
    call write_data (fname, 'rtpthlp_3D',    rtpthlp_3D)
    call write_data (fname, 'up2_3D',        up2_3D)
    call write_data (fname, 'vp2_3D',        vp2_3D)

    if( sclr_dim > 0 ) then
      call write_data (fname, 'RH_crit_3D1', RH_crit_clubb_3D( : , : , : ,1))
      call write_data (fname, 'RH_crit_3D2', RH_crit_clubb_3D( : , : , : ,2))
    endif ! sclr_dim > 0
endif ! do_netcdf_restart

call cleanup_CLUBB_3D_var( )

call cleanup_clubb_core( .true. )
! De-allocate the array for the passive scalar tolerances

deallocate( momentum_heights )
deallocate( thermodynamic_heights )       
if( do_CLUBB_conservation_checks  ) then
    if( allocated( id_tracer_CLUBB ) )      deallocate (id_tracer_CLUBB)
    if( allocated( id_tracer_CLUBB_col ) )  deallocate (id_tracer_CLUBB_col)
    if( allocated( pmass_3d ) )             deallocate (pmass_3d)
    if( allocated( tmp2D_check ) )          deallocate (tmp2D_check)
endif
    
!---------------------------------------------------------------------
!    mark the module as uninitialized.
!---------------------------------------------------------------------
module_is_initialized = .false.

return
end subroutine CLUBB_END
!#######################################################################



  ! ----------------------------------------------------------------------------
  ! This subroutine adds model tendencies to the CLUBB prognostic variables
  
  subroutine add_host_tdcy(                                                    &
               ix_host_model, iy_host_model, CLUBB_dt,                         &
               udt, vdt, tdt, rdt, env_qv_scale, env_condensate_scale,         &
               u_star, b_star, q_star,                                         &
               um, vm, thlm, thvm, rtm, cloud_frac, rcm, edsclrm, sclrm,       &
               upwp_sfc, vpwp_sfc, wpthlp_sfc, wprtp_sfc )

  implicit none

  ! Calling arguments
  integer, intent(in)                  :: ix_host_model, iy_host_model
  real, intent(in)                     :: CLUBB_dt   ! time-step to apply tendencies
  real, intent(in), dimension(:,:,:)   :: udt, vdt, tdt
  real, intent(in), dimension(:,:,:,:) :: rdt
  real, intent(in), dimension(:,:,:)   :: env_qv_scale, env_condensate_scale
  real, intent(in), dimension(:,:)     :: u_star, b_star, q_star

  real, intent(inout), dimension(:)    :: um, vm, thlm, thvm, rtm, cloud_frac, rcm
  real, intent(inout), dimension(:,:)  :: edsclrm, sclrm
  real, intent(out)                    :: upwp_sfc, vpwp_sfc, wpthlp_sfc, wprtp_sfc

  ! Local variables
  real, dimension(size(thlm,1)) :: tmp_CLUBB

  ! Cloud fraction
  tmp_CLUBB = 0.0
  call host2CLUBB_full(rdt(ix_host_model, iy_host_model, :, nqa),  &  ! intent (in)
                       tmp_CLUBB)                                     ! intent (out)
  cloud_frac    = cloud_frac + CLUBB_dt*tmp_CLUBB
  cloud_frac(1) = cloud_frac(2)
  where ( cloud_frac(:) <=cloud_frac_min )
    cloud_frac(:) = 0.0
  end where

  ! Water vapor
  tmp_CLUBB = 0.0
  call host2CLUBB_full ( rdt(ix_host_model, iy_host_model, : , nsphum)    &
                         /env_qv_scale(ix_host_model, iy_host_model, :),  &  !intent (in)
                         tmp_CLUBB )                                         !intent (out)
  rtm = rtm + CLUBB_dt * tmp_CLUBB

  ! Cloud water
  tmp_CLUBB = 0.0
  call host2CLUBB_full ( rdt( ix_host_model, iy_host_model, : , nql)              &
                         /env_condensate_scale(ix_host_model, iy_host_model, :),  &  !intent (in)
                         tmp_CLUBB )                                                 !intent (out) 
  rcm = rcm + CLUBB_dt * tmp_CLUBB
  rtm = rtm + CLUBB_dt * tmp_CLUBB

  ! Cloud ice
  tmp_CLUBB = 0.0
  call host2CLUBB_full ( rdt(ix_host_model, iy_host_model, : , nqi)               &
                         /env_condensate_scale(ix_host_model, iy_host_model, :),  &  !intent (in)
                         tmp_CLUBB )                                                 !intent (out)
  ! regard ice water as liquid water
  rcm = rcm + CLUBB_dt * tmp_CLUBB
  rtm = rtm + CLUBB_dt * tmp_CLUBB
  rcm(1) = rcm(2)
  rtm(1) = rtm(2)

  ! Droplet number concentration
  if (nqn > 0) then
    tmp_CLUBB = 0.0
    call host2CLUBB_full( rdt(ix_host_model, iy_host_model, : , nqn)               &
                          /env_condensate_scale(ix_host_model, iy_host_model, :),  & !intent (in)
                          tmp_CLUBB )                                                !intent (out)
    if (use_sclr_HOC) then
      sclrm( :, 1 ) = sclrm( :, 1 ) + CLUBB_dt * tmp_CLUBB(:)
      sclrm( 1, 1 ) = sclrm( 2, 1 )
      do  iz_clubb = 1, nz_clubb
        if( cloud_frac(iz_clubb) <= cloud_frac_min ) sclrm( iz_clubb , 1) = 0.0
      enddo
    else
      edsclrm( :, 1) = edsclrm( : , 1 ) + CLUBB_dt * tmp_CLUBB(:)
      edsclrm( 1, 1) = edsclrm( 2, 1 )
      do  iz_clubb = 1, nz_clubb
        if( cloud_frac(iz_clubb) <= cloud_frac_min ) edsclrm( iz_clubb , 1) = 0.0
      enddo              
    endif
  endif

  ! Ice number concentration
  if (nqni > 0) then
    tmp_CLUBB = 0.0
    call host2CLUBB_full( rdt( ix_host_model, iy_host_model, : , nqni)             &
                          /env_condensate_scale(ix_host_model, iy_host_model, :),  & !intent (in)
                          tmp_CLUBB)                                                 !intent (out)
    if (use_sclr_HOC) then
      sclrm( : , 2 ) =  sclrm( : , 2 ) + CLUBB_dt * tmp_CLUBB( : )
      sclrm( 1 , 2 ) = sclrm( 2 , 2 )
      do  iz_clubb = 1, nz_clubb
        if ( cloud_frac(iz_clubb) <= cloud_frac_min ) sclrm( iz_clubb , 2 ) = 0.0
      enddo
    else
      edsclrm( : , 2) = edsclrm( : , 2) + CLUBB_dt * tmp_CLUBB( :)
      edsclrm( 1 , 2) = edsclrm( 2 , 2)
      do  iz_clubb = 1, nz_clubb
        if( cloud_frac(iz_clubb) <= cloud_frac_min ) edsclrm( iz_clubb , 2) = 0.0
      enddo
    endif
  endif

  ! Temperature
  call host2CLUBB_full( tdt(ix_host_model, iy_host_model, :),  & !intent (in)
                        tmp_CLUBB )                              !intent (out)
  temp_clubb = temp_clubb + CLUBB_dt * tmp_CLUBB

  ! Compute thlm, thvm
  thlm = (temp_clubb - Lv*rcm/Cp)/exner
  thlm(1) = thlm(2)
  thvm = thlm + ep1*T0*rtm + (Lv/(Cp*exner) - ep2*T0) * rcm

  ! Calculate surface u (upwp_sfc) and v (vpwp_sfc) momentum fluxes from u_star (input)
  upwp_sfc = -um(2) * u_star(ix_host_model, iy_host_model)**2  &
             / (max(gust_const, sqrt( um(2)*um(2) + vm(2)*vm(2)) ))
  vpwp_sfc = -vm(2) * u_star(ix_host_model, iy_host_model)**2  &
             / (max(gust_const, sqrt( um(2)*um(2) + vm(2)*vm(2)) ))

  ! Calculate surface sensible (wpthlp_sfc) and latend (wprtp_sfc) fluxes from u_star, b_star, and q_star (inputs)
  wpthlp_sfc = b_star( ix_host_model, iy_host_model ) * u_star( ix_host_model, iy_host_model ) *292./grav
  wprtp_sfc  = q_star( ix_host_model, iy_host_model ) * u_star( ix_host_model, iy_host_model )

  ! u wind [m/s]
  tmp_CLUBB = 0.0
  call host2CLUBB_full( udt( ix_host_model, iy_host_model, : ),  & !intent (in)
                        tmp_CLUBB)                                 !intent (out)
  um = um + CLUBB_dt * tmp_CLUBB
  um( 1 ) = um( 2 )

  ! v wind [m/s]
  tmp_CLUBB = 0.0
  call host2CLUBB_full( vdt( ix_host_model, iy_host_model, : ),  & !intent (in)
                        tmp_CLUBB )                                !intent (out)
  vm = vm + CLUBB_dt * tmp_CLUBB
  vm( 1 ) = vm( 2 )

  return
  end subroutine add_host_tdcy

end module  CLUBB_driver_mod
