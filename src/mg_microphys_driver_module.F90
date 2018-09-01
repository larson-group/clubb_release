! $Id$
module mg_microphys_driver_module

  implicit none

  public :: mg_microphys_driver

  private

  contains
!-------------------------------------------------------------------------------
  subroutine mg_microphys_driver &
             ( dt, nz, l_stats_samp, dzt, thlm, p_in_Pa, exner, &
               rho, cloud_frac, rcm, Ncm, rvm, Nccnm, pdf_params, hydromet, &
               hydromet_mc, hydromet_vel, rcm_mc, rvm_mc, thlm_mc )
! Description:
!   Wrapper for the Morrison-Gettelman microphysics
! References:
!   None
!-------------------------------------------------------------------------------

    use parameters_model, only: hydromet_dim

    use T_in_K_module, only: &
      T_in_K2thlm, & ! Procedure(s)
      thlm2T_in_K

    use array_index, only:  &
      iiri, iiNi

    use constants_clubb, only: &
      T_freeze_K, &
!      Lv, &
      Cp
      
    use stats_variables, only: & 
      stats_zt,            &  ! Variables
      stats_zm,            &
      stats_sfc,           &
      irsm,       &
      irrm,       &
      iswp,          &
      iprecip_rate_sfc, &
      ieff_rad_cloud, &
      ieff_rad_ice,  &
      ieff_rad_rain, &
      ieff_rad_snow, &
      irrm_auto,  &
      irrm_accr,  &
      irwp,          &
      ircm_in_cloud, &
      iNsm,       &
      iNrm,          &
      iVNr,          &
      iVrr,          &
      iVNc,          &
      iVrc,          &
      iVNs,       &
      iVrs,       &
      iVNi,        &
      iVri,        &
      ircm_sd_mg_morr,&
      irim_sd_mg_morr

    use stats_type_utilities, only:  & 
      stat_update_var, &
      stat_update_var_pt
      
    use cldwat2m_micro, only: &
      mmicro_pcond ! Procedure
      
    use microp_aero, only: &
      microp_aero_ts ! Procedure

    use variables_diagnostic_module, only: &
      Kh_zm, &
      em

    use wv_saturation, only: &
      gestbl ! Procedure
      
    use physconst, only: &
      cpair,  &
      rh2o,   &
      tmelt,  &
      latice, &
      latvap, &
      epsilo
      
    use ppgrid, only: &
      init_ppgrid ! Procedure

    use ppgrid, only: &
      pcols ! Constant(s)
      
    use phys_buffer, only: &
      pbuf        ! Variable
      
    use grid_class, only: &
      flip ! Procedure
      
    use constants_clubb, only: &
      rc_tol, & ! Variable
      sec_per_day, &
      mm_per_m

    use shr_kind_mod, only: r8 => shr_kind_r8 

    use clubb_precision, only: &
      core_rknd, & ! Variable(s)
      dp

    use fill_holes, only: &
      vertical_integral ! Procedure(s)

    use pdf_parameter_module, only: &
        pdf_parameter  ! Variable(s)

    implicit none

    ! External
    intrinsic :: max

    ! Constant Parameters
    integer, parameter :: lchnk = 0
    logical, parameter :: sub_column = .false.

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: dt ! Model timestep        [s]

    integer, intent(in) :: nz ! Points in the Vertical        [-]

    logical, intent(in) :: &
      l_stats_samp  ! Whether to accumulate statistics [T/F]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      dzt,        & ! Inverse of the grid spacing        [1/m]
      thlm,       & ! Liquid potential temperature       [K]
      p_in_Pa,    & ! Pressure                           [Pa]
      exner,      & ! Exner function                     [-]
      rho,        & ! Density on thermo. grid            [kg/m^3]
      cloud_frac    ! Cloud fraction                     [-]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      rcm,   & ! Liquid water mixing ratio                [kg/kg]
      Ncm,   & ! Cloud droplet number concentration       [count/kg]
      rvm,   & ! Vapor water mixing ratio                 [kg/kg]
      Nccnm    ! Cloud condensation nuclei concentration  [count/kg]

    type(pdf_parameter), dimension(nz), intent(in) :: &
      pdf_params    ! PDF parameters                   [units vary]

    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(in) :: &
      hydromet      ! Hydrometeor species    [units vary]

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(out) :: &
      hydromet_mc, &  ! Hydrometeor time tendency          [(units vary)/s]
      hydromet_vel    ! Hydrometeor sedimentation velocity [m/s]

    ! Output Variables
    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      rcm_mc, & ! Time tendency of liquid water mixing ratio    [kg/kg/s]
      rvm_mc, & ! Time tendency of vapor water mixing ratio     [kg/kg/s]
      thlm_mc   ! Time tendency of liquid potential temperature [K/s]

    ! Local Variables

    logical :: l_microp_aero_ts ! Use microp_aero_ts to determine aerosols

    real(r8) :: cldfsnow ! Radiative cloud fraction at sfc level for calcuating swp   [-]

    real( kind = core_rknd ), dimension(nz) :: & 
      T_in_K,       & ! Temperature                                                   [K]
      T_in_K_new,   & ! Temperature after microphysics                                [K]
!      T_in_K_mc,    & ! Temperature tendency due to microphysics                      [K/s]
      rcm_new,      & ! Cloud water mixing ratio after microphysics                   [kg/kg]
      rsm,       & ! Snow mixing ratio (not in hydromet because it is diagnostic)  [kg/kg]
      rrm,       & ! Rain mixing ratio (not in hydromet because it is diagnostic)  [kg/kg]
      effc,         & ! Droplet effective radius                                      [μ]
      effi,         & ! cloud ice effective radius                                    [μ]
      reff_rain,    & ! rain effective radius                                         [μ]
      reff_snow,    & ! snow effective radius                                         [μ]
      rim_sten,   & ! Ice sedimentation tendency                                    [kg/kg/s]
      rcm_sten        ! Cloud dropet sedimentation tendency                           [kg/kg/s]

    real(r8), dimension(pcols,nz) :: & 
      turbtype_flip,& ! Turbulence type at each interface                             [-]
      smaw_flip,    & ! Normalized instability function of momentum                   [???]
      wsub_flip,    & ! Diagnosed sub-grid vertical velocity st. dev.                 [m/s]
      wsubi_flip      ! Diagnosed sub-grid vertical velocity ice                      [m/s]
      
    real(r8), dimension(1, nz-1, 1) :: &
      aer_mmr_flip

    ! Note that the MG grid is flipped with respect to CLUBB.
    ! For this reason, all MG variables are marked as "flip"

    ! MG Input Variables on zm grid
    real(r8), dimension(pcols,nz) :: &
      Kh_zm_flip,   & ! Eddy diffusivity coefficient on momentum levels      [m^2/s]
      em_flip         ! Turbulent Kinetic Energy (TKE)                       [m^2/s^2]

    ! MG Input Variables on zt grid
    real(r8), dimension(pcols,nz-1) :: &
      T_in_K_flip,  & ! Air temperature                                      [K]
      rvm_flip,     & ! Vapor water mixing ratio                             [kg/kg]
      rcm_flip,     & ! Cloud water mixing ratio                             [kg/kg]
      p_in_Pa_flip, & ! Pressure                                             [Pa]
      pdel_flip,    & ! difference in air pressure between vertical levels   [Pa]
      cldn_flip,    & ! Relative Humidity Cloud fraction                     [-]
      liqcldf_flip, & ! Liquid cloud fraction                                [-]
      icecldf_flip, & ! Ice cloud fraction                                   [-]
      naai_flip,    & ! number of activated ice nuclei                       [1/kg]
      rho_flip,     & ! Density on thermo. grid                              [kg/m^3]
      npccn_flip,   & ! Number of cloud nuclei                               [count/kg]
      unused_in       ! Represents MG variables that are not used in the current code
      
    real(r8), dimension(pcols,nz-1,4) :: &
      rndst_flip,   & ! radius of 4 dust bins for contact freezing [units unknown, our guess is mm]
      nacon_flip      ! number of 4 dust bins for contact nucleation [units unknown]
      
    ! MG Output Variables
    real(r8), dimension(pcols) :: &
      prect,        & ! Surface precipitation rate                                    [m/s]
      preci           ! Unused

    real(r8), dimension(pcols,nz-1) :: &
      tlat_flip,         & ! Latent heating rate                       [W/kg]
      rcm_mc_flip,       & ! Time tend. of liquid water mixing ratio   [kg/kg/s]
      rvm_mc_flip,       & ! Time tendency of vapor water mixing ratio [kg/kg/s]
      cldo_flip,         & ! Old cloud fraction.                       [-]
      qsout_flip,        & ! Snow mixing ratio                         [kg/kg]
      effc_flip,         & ! Droplet effective radius                  [μ]
      effi_flip,         & ! cloud ice effective radius                [μ]
      reff_rain_flip,    & ! rain effective radius                     [μ]
      reff_snow_flip,    & ! snow effective radius                     [μ]
      qrout_flip,        & ! rain mixing ratio                         [kg/kg]
      T_in_K_flip_inout, & ! Absolute temperature                      [K]
      T_in_K_flip_out      ! Absolute temperature (before reset)       [K]
      
    ! MG zt variables that are not used in CLUBB.
    real(r8), dimension(pcols,nz-1) :: &
      rate1ord_cw2pr_st_flip, effc_fn_flip,                                    &
      nevapr_flip, evapsnow_flip, prain_flip, prodsnow_flip, cmeout_flip,      &
      deffi_flip, pgamrad_flip, lamcrad_flip, dsout_flip, qcsevap_flip,        &
      qisevap_flip, qvres_flip, cmeiout_flip, vtrmc_flip, vtrmi_flip,          &
      qcsedten_flip, qisedten_flip, prao_flip, prco_flip, mnuccco_flip,        &
      mnuccto_flip, msacwio_flip, psacwso_flip, bergso_flip, bergo_flip,       &
      melto_flip, homoo_flip, qcreso_flip, qireso_flip,&
      mnuccro_flip, pracso_flip, meltsdt_flip, frzrdt_flip, mnuccdo_flip,      &
      naai_hom_flip

    ! MG zt variables that are diagnostics.
    real(r8), dimension(pcols,nz-1) :: &
      prcio_flip, praio_flip, &
      qcic_flip, &  ! In-cloud cloud liquid mixing ratio
      nsic_flip, &  ! In-precip snow number concentration
      nric_flip, &  ! In-precip rain number concentration
      cldmax        ! Precip frac assuming max overlap

    ! The above variables flipped to the CLUBB zt grid
    real( kind = core_rknd ), dimension(nz) :: &
      rrm_auto, &  ! Autoconversion rate                 [kg/kg/s]
      rrm_accr, &  ! Accretion rate                      [kg/kg/s]
      rcm_in_cloud, & ! Liquid water in cloud               [kg/kg]
      Nsm,       & ! Snow particle number concentration  [num/kg]
      Nrm             ! Rain droplet number concentration   [num/kg]

    ! MG zm variables that are diagnostics.
    real(r8), dimension(nz-1) :: &
      uns_flip,  &  ! Number-weighted snow fallspeed
      ums_flip,  &  ! Mass-weighted snow fallspeed
      unr_flip,  &  ! Number-weighted rain fallspeed
      umr_flip      ! Mass-weighted rain fallspeed
      
    real(r8) ::  &
      uni,  &  ! Number-weighted cloud ice fallspeed
      umi,  &  ! Mass-weighted cloud ice fallspeed
      unc,  &  ! Number-weighted cloud droplet fallspeed
      umc      ! Mass-weighted cloud droplet fallspeed

    ! The above variables flipped to the CLUBB zm grid
    real( kind = core_rknd ), dimension(nz) :: &
      VNs,       & ! Number-weighted snow sedimentation velocity          [m/s]
      Vrs,       & ! Mass-weighted snow sedimentation velocity            [m/s]
      VNr,          & ! Number-weighted rain sedimentation velocity          [m/s]
      Vrr,          & ! Mass-weighted rain sedimentation velocity            [m/s]
      VNi,        & ! Number-weighted cloud ice sedimentation velocity     [m/s]
      Vri,        & ! Mass-weighted cloud ice sedimentation velocity       [m/s]
      VNc,          & ! Number-weighted cloud droplet sedimentation velocity [m/s]
      Vrc             ! Mass-weighted cloud droplet sedimentation velocity   [m/s]


    ! MG zm variables that are not used in CLUBB
    real(r8), dimension(pcols,nz) :: &
      rflx_flip, sflx_flip
      
    real(r8), dimension(pcols,nz-1) :: &
      Ncm_flip, & ! Flipped array for cloud droplet number concentration   [#/kg]
      Ncm_mc_flip

    real(r8), dimension(pcols,nz-1,hydromet_dim) :: &
       hydromet_flip,    &  ! Hydrometeor species                            [units vary]
       hydromet_mc_flip     ! Hydrometeor time tendency                      [units vary]

    type(pdf_parameter), dimension(nz-1) :: &
      pdf_params_flip    ! PDF parameters flipped for MG grid       [units vary]

    real( kind = core_rknd ) :: xtmp

    integer :: i, k, ncols, icol

    ! ---- Begin Code ----

    ncols = pcols ! This should be 1

    ! icol is the current subcolumn;  since sub column sampling is disabled icol is always 1
    icol = 1 

    l_microp_aero_ts = .true.
    
    unused_in(1:pcols,:) = -9999.9_r8
    
    ! Initialize output arrays to zero
    naai_flip = 0.0_r8
    npccn_flip = 0.0_r8
    rndst_flip = 0.0_r8
    nacon_flip = 0.0_r8
    effc(:) = 0.0_core_rknd
    effi(:) = 0.0_core_rknd
    reff_rain(:) = 0.0_core_rknd
    reff_snow(:) = 0.0_core_rknd
    rsm(:) = 0.0_core_rknd
    rrm(:) = 0.0_core_rknd
    rcm_new(:) = 0.0_core_rknd
    T_in_K_new(:) = 0.0_core_rknd
!    T_in_K_mc(:) = 0.0_core_rknd
    rcm_mc(:) = 0.0_core_rknd
    rvm_mc(:) = 0.0_core_rknd
    hydromet_mc(:,:) = 0.0_core_rknd
    rcm_sten(:) = 0.0_core_rknd
    rim_sten(:) = 0.0_core_rknd
    hydromet_mc_flip(:,:,:) = 0.0_r8

    ! Determine temperature
    T_in_K = thlm2T_in_K( thlm, exner, rcm )
    
    ! MG's grid is flipped with respect to CLUBB.
    ! Flip CLUBB variables before inputting them into MG.
    Kh_zm_flip(icol,1:nz) = real( flip( real(Kh_zm(1:nz), kind=dp ), nz ), kind=r8 )
    em_flip(icol,1:nz) = real( flip( real(em(1:nz), kind=dp ), nz ), kind=r8 )

    ! MG does not include zt(1) as CLUBB does, so it is removed during the flip
    T_in_K_flip(icol,1:nz-1) = real( flip( real(T_in_K(2:nz),kind=dp ), nz-1 ), kind=r8 )
    rvm_flip(icol,1:nz-1) = real( flip( real(rvm(2:nz),kind=dp ), nz-1 ), kind=r8 )
    rcm_flip(icol,1:nz-1) = real( flip( real(rcm(2:nz),kind=dp ), nz-1 ), kind=r8 )
    Ncm_flip(icol,1:nz-1) = real( flip( real(Ncm(2:nz),kind=dp ), nz-1 ), kind=r8 )
    p_in_Pa_flip(icol,1:nz-1) = real( flip( real(p_in_Pa(2:nz),kind=dp ), nz-1 ), kind=r8 )
    liqcldf_flip(icol,1:nz-1) = real( flip( real(cloud_frac(2:nz),kind=dp ), nz-1 ), kind=r8 )

    ! Hydromet is 2 dimensional, so flip function doesn't work
    do i = 1, hydromet_dim, 1
      hydromet_flip(icol,1:nz-1, i) = real( hydromet(nz:2:-1, i), kind=r8 )
    end do

    ! Flip the CLUBB PDF parameters.
    pdf_params_flip(1:nz-1) = pdf_params(nz:2:-1)
    
    ! Initialize MG input variables. The top level is skipped because MG's grid is flipped with
    ! respect to CLUBB's, and MG doesn't include CLUBB's lowest level.
    do k = 1, nz-1, 1
    
      ! Find difference in air pressure between vertical levels
      if ( k /= nz-1 ) then
        pdel_flip(icol,k) = p_in_Pa_flip(icol,k+1) - p_in_Pa_flip(icol,k)
      else
        pdel_flip(icol,k) = pdel_flip(icol,k-1)
      end if
      
      ! Cloud fraction. In MG there is no difference between ice cloud fraction and
      ! liquid cloud fraction. However, just setting icecldf equal to CLUBBs cloud_frac
      ! won't work for all ice, no liquid clouds.
      ! In CAM, cldn is equal to max(liqcldf,icecldf).
      if ( rcm_flip(icol,k) > 0._r8 ) then
        cldn_flip(icol,k) = liqcldf_flip(icol,k)
      else if ( rcm_flip(icol,k) < real( rc_tol, kind=r8 ) &
                .and. hydromet_flip(icol,k,iiri) > real( rc_tol, kind=r8 ) ) then
        cldn_flip(icol,k) = 1._r8
      else
        cldn_flip(icol,k) = 0._r8
      end if

      liqcldf_flip(icol,k) = cldn_flip(icol,k)
      icecldf_flip(icol,k) = cldn_flip(icol,k)

    end do ! k = 1 .. nz-1
    
    ! Initialize grid variables. These are imported in the MG code.
    call init_ppgrid( nz )

    call gestbl( 173.16_r8, 375.16_r8, 20.00_r8, .true., real( epsilo, kind=r8), &
                 real( latvap, kind=r8), real( latice, kind=r8), real( rh2o, kind=r8), &
                 real( cpair, kind=r8), real( tmelt, kind=r8) ) ! Known magic flag

    ! Ensure no hydromet arrays are input as 0, because this makes MG crash.
    ! dschanen found this was no longer necessary, due a prior bug fix? 5 Jan 2011
!   do i=1, hydromet_dim, 1
!     hydromet_flip(:,i) = max( 1e-8_r8, hydromet_flip(:,i) ) ! Known magic number
!   end do

    ! These variables are unused in CLUBB, so they are initialized to 1 before input into
    ! microp_aero_ts. 
    aer_mmr_flip(icol,:,:) = 1._r8
    turbtype_flip(icol,:) = 1._r8
    smaw_flip(icol,:) = 1._r8
   
    if( l_microp_aero_ts ) then
      ! Need to prevent a floating-point exception below;  the result doesn't
      ! appear to feedback into the model
      aer_mmr_flip(icol,:,:) = 0._r8

      ! Calculate aerosol activiation, dust size, and number for contact nucleation
      call microp_aero_ts &
         ( lchnk, ncols, real( dt, kind=r8), T_in_K_flip, unused_in, &              ! in
         rvm_flip, rcm_flip, hydromet_flip(:,:,iiri), &                          ! in
         Ncm_flip, hydromet_flip(:,:,iiNi), p_in_Pa_flip, pdel_flip, cldn_flip, &  ! in
         liqcldf_flip, icecldf_flip, &                                              ! in
         cldo_flip, unused_in, unused_in, unused_in, unused_in, &                   ! in
         aer_mmr_flip, &
         pbuf, &
         Kh_zm_flip, em_flip, turbtype_flip, smaw_flip, wsub_flip, wsubi_flip, &
         naai_flip, naai_hom_flip, npccn_flip, rndst_flip, nacon_flip )             ! out

    else

      rho_flip(icol,1:nz-1) = real( flip( real(rho(2:nz),kind=dp ), nz-1 ), kind=r8 )
      npccn_flip(icol,1:nz-1) &
         = real( flip( real( Nccnm(2:nz), kind=dp ), nz-1 ), kind=r8 )

      ! Determine ice nulceation number using Meyers formula found in the Morrison microphysics
      do k = 1, nz-1
        naai_flip(icol,k) = exp( -2.80_r8 + 0.262_r8 * (  real( T_freeze_K,kind=dp )  &
                                                        - real( T_in_K_flip(icol,k),kind=dp ) ) &
                               ) &
                       * 1000.0_r8 / rho_flip(icol,k)
      end do ! k = 1 ... nz-1

      ! Set values for radius of 4 dust bins for contact freezing. Currently, we assume
      ! nacon = 0 therby ommitting contact nucleation, but this is set regardless.
      rndst_flip(icol,:,1) = 0.001_r8
      rndst_flip(icol,:,2) = 0.01_r8
      rndst_flip(icol,:,3) = 0.1_r8
      rndst_flip(icol,:,4) = 1.0_r8
      nacon_flip(icol,:,1:4) = 0.0_r8

    end if 

    qcic_flip(icol,:) = 0.0_r8 ! This is needed in case qc is 0 everywhere, due to a goto statement
    nsic_flip(icol,:) = 0.0_r8
    nric_flip(icol,:) = 0.0_r8
    uni = 0.0_r8
    umi = 0.0_r8
    uns_flip(:) = 0.0_r8
    ums_flip(:) = 0.0_r8
    unr_flip(:) = 0.0_r8
    umr_flip(:) = 0.0_r8
    unc = 0.0_r8
    umc = 0.0_r8

    T_in_K_flip_inout = T_in_K_flip ! Set the output variable to the current values
    T_in_K_flip_out   = T_in_K_flip ! Set the output variable to the current values

    ! Call the Morrison-Gettelman microphysics
    call mmicro_pcond                                                                        &! in
         ( sub_column,                                                                       &! in
           lchnk, ncols, real( dt, kind=r8), T_in_K_flip,                                    &! in
           rvm_flip, rcm_flip, hydromet_flip(:,:,iiri),                                    &! in
           Ncm_flip, hydromet_flip(:,:,iiNi), p_in_Pa_flip, pdel_flip, cldn_flip,            &
           liqcldf_flip, icecldf_flip,                                                        &! in
           cldo_flip,                                                                         &! in
           pdf_params_flip,                                                                   & ! in
           rate1ord_cw2pr_st_flip,                                                            &! out
           naai_flip, npccn_flip, rndst_flip, nacon_flip,                                     &! in
           tlat_flip, rvm_mc_flip,                                                            &! out
           rcm_mc_flip, hydromet_mc_flip(:,:,iiri), Ncm_mc_flip,                           &! out
           hydromet_mc_flip(:,:,iiNi), effc_flip,                                            &! out
           effc_fn_flip, effi_flip, prect, preci,                                             &! out
           nevapr_flip, evapsnow_flip,                                                        &! out
           prain_flip, prodsnow_flip, cmeout_flip, deffi_flip, pgamrad_flip,                  &! out
           lamcrad_flip, qsout_flip, dsout_flip,                                              &! out
           rflx_flip, sflx_flip, qrout_flip, reff_rain_flip, reff_snow_flip,                  &! out
           qcsevap_flip, qisevap_flip, qvres_flip, cmeiout_flip,                              &! out
           vtrmc_flip, vtrmi_flip, qcsedten_flip, qisedten_flip,                              &! out
           prao_flip, prco_flip, mnuccco_flip, mnuccto_flip, msacwio_flip, psacwso_flip,      &! out
           bergso_flip, bergo_flip, melto_flip, homoo_flip, qcreso_flip, prcio_flip, praio_flip, &
           qireso_flip,                                                                       &! out
           mnuccro_flip, pracso_flip, meltsdt_flip, frzrdt_flip, mnuccdo_flip,                &! out
           T_in_K_flip_out,                                                                   &! out
           qcic_flip, T_in_K_flip_inout, nsic_flip, nric_flip, uni, umi,                      &! out
           uns_flip, ums_flip, unr_flip, umr_flip, unc, umc, cldmax )                          ! out

    ! Flip MG variables into CLUBB grid
    rcm_mc(2:nz) = real( flip( real(rcm_mc_flip(icol,1:nz-1),kind=dp ), nz-1 ), kind = core_rknd )
    ! Since we fix Ncm the change in Ncm isn't needed
!   Ncm_mc(2:nz) = real( flip( real(Ncm_mc_flip(icol,1:nz-1),kind=dp ), nz-1 ), kind = core_rknd )
    rvm_mc(2:nz) = real( flip( real(rvm_mc_flip(icol,1:nz-1),kind=dp ), nz-1 ), kind = core_rknd )
    rcm_new(2:nz) = real( flip( real(rcm_flip(icol,1:nz-1),kind=dp ), nz-1 ), kind = core_rknd )
    effc(2:nz) = real( flip( real(effc_flip(icol,1:nz-1),kind=dp ), nz-1 ), kind = core_rknd )
    effi(2:nz) = real( flip( real(effi_flip(icol,1:nz-1),kind=dp ), nz-1 ), kind = core_rknd )
    reff_rain(2:nz) = real( flip( real(reff_rain_flip(icol,1:nz-1),kind=dp ), nz-1 ), &
                            kind = core_rknd )
    reff_snow(2:nz) = real( flip( real(reff_snow_flip(icol,1:nz-1),kind=dp ), nz-1 ), &
                            kind = core_rknd )
    rsm(2:nz) = real( flip( real(qsout_flip(icol,1:nz-1),kind=dp ), nz-1 ), kind = core_rknd )
    rrm(2:nz) = real( flip( real(qrout_flip(icol,1:nz-1),kind=dp ), nz-1 ), kind = core_rknd )
    rcm_sten(2:nz) = real( flip( real(qcsedten_flip(icol, 1:nz-1),kind=dp ), nz-1 ), &
                            kind = core_rknd )
    rim_sten(2:nz) = real( flip( real(qisedten_flip(icol, 1:nz-1),kind=dp ), nz-1 ), &
                               kind = core_rknd )
    
    do i = 1, hydromet_dim, 1      
      hydromet_mc(2:nz, i) = real( hydromet_mc_flip(icol,nz-1:1:-1, i), kind = core_rknd )
    end do

    ! Compute the snow water path
    cldfsnow = cldn_flip(icol,nz-1) ! Only need sfc level
    if ( ( cldfsnow > 1.e-4_r8 ) .and. ( rcm_flip(icol,nz-1) < 1e-10_r8 ) ) then
      cldfsnow = 0._r8
    else if ( ( cldfsnow < 1.e-4_r8 ) .and. ( qsout_flip(icol,nz-1) > 1.e-6_r8 ) ) then
      cldfsnow = 0.25_r8
    end if
    
    ! Update thetal based on absolute temperature. We use this rather than tlat
    ! because using the latter seemed to cause spurious heating effects.
    T_in_K_new(2:nz) = real( flip( real(T_in_K_flip_out(icol,1:nz-1),kind=dp ), nz-1 ), &
                             kind = core_rknd )
    T_in_K_new(1) = T_in_K(1)
   
    ! Compute total change in thlm using ( thlm_new - thlm_old ) / dt
    thlm_mc = ( T_in_K2thlm( T_in_K_new, exner, rcm_new ) - thlm ) / dt

    ! Rate of change of thlm due to microphysics
    !
    ! The rate of change of mean temperature with respect to time due to
    ! microphysics, T_in_K_mc, is based on tlat, where dT/dt = tlat / Cp.
    ! The rate of change of thlm with respect to time due to microphysics,
    ! thlm_mc, is given by:
    !
    ! dthlm/dt|_mc = ( 1 / exner ) * dT/dt - ( Lv / ( Cp *exner ) ) * drcm/dt;
    !
    ! where drcm/dt is the rate of change of mean cloud water mixing ratio
    ! with respect to time, rcm_mc.  The equation for dthlm/dt|_mc can be
    ! written in terms of tlat:
    !
    ! dthlm/dt|_mc = ( 1 / ( Cp * exner ) ) * ( tlat - Lv * drcm/dt )

!    T_in_K_mc(2:nz) = real( flip( real( tlat_flip(icol,1:nz-1),kind=dp ), nz-1 ), &
!                            kind = core_rknd ) / Cp
!    T_in_K_mc(1) = 0.0_core_rknd
!
!    thlm_mc = ( T_in_K_mc / exner ) - ( Lv / ( Cp * exner ) ) * rcm_mc
    
    if ( l_stats_samp ) then

      ! Update rain and snow water mixing ratio (computed diagnostically)
      call stat_update_var( irsm, rsm(:), stats_zt)
      call stat_update_var( irrm, rrm(:), stats_zt)
      
      ! Effective radii of hydrometeor species
      call stat_update_var( ieff_rad_cloud, effc(:), stats_zt )
      call stat_update_var( ieff_rad_ice, effi(:), stats_zt )
      call stat_update_var( ieff_rad_rain, reff_rain(:), stats_zt)
      call stat_update_var( ieff_rad_snow, reff_snow(:), stats_zt)

      ! Rain rates at the bottom of the domain, in mm/day
      call stat_update_var_pt( iprecip_rate_sfc, 1, &
                               real( prect(icol), kind = core_rknd ) * mm_per_m * &
                               sec_per_day, stats_sfc )

     ! Snow water path is updated in stats_subs.F90
     ! call stat_update_var_pt( iswp, 1, real( hydromet_mc(3,iirs) / max( 0.0001, cldfsnow ) * &
     !                          pdel_flip(nz-1) / gravit ), stats_sfc )

      ! Compute autoconversion
      rrm_auto(2:nz) = real( flip( real( prco_flip(icol,1:nz-1),kind=dp ), nz-1 ), &
                                kind = core_rknd )
      rrm_auto(1) = 0.0_core_rknd
      call stat_update_var( irrm_auto, rrm_auto, stats_zt )

      ! Compute accretion
      rrm_accr(2:nz) = real( flip( real( prao_flip(icol,1:nz-1),kind=dp ), nz-1 ), &
                                kind = core_rknd )
      rrm_accr(1) = 0.0_core_rknd
      call stat_update_var( irrm_accr, rrm_accr, stats_zt )

      ! Compute cloud water mixing ratio in cloud (based on a threshold in the
      ! MG microphysics)
      rcm_in_cloud(2:nz) = real( flip( real( qcic_flip(icol,1:nz-1),kind=dp ), nz-1 ), &
                                 kind = core_rknd )
      rcm_in_cloud(1) = 0.0_core_rknd
      call stat_update_var( ircm_in_cloud, rcm_in_cloud, stats_zt )

      ! Compute in-precipitation snow number concentration
      Nsm(2:nz) = real( flip( real( nsic_flip(icol, 1:nz-1),kind=dp ), nz-1 ), &
                           kind = core_rknd ) * real( cldmax(icol, 1:nz-1), kind = core_rknd )
      Nsm(1) = 0.0_core_rknd
      call stat_update_var( iNsm, Nsm, stats_zt )

      ! Compute in precipitation rain number concentration
      Nrm(2:nz) = real( flip( real( nric_flip(icol, 1:nz-1),kind=dp ), nz-1 ), kind = core_rknd ) &
                  * real( cldmax(icol, 1:nz-1), kind = core_rknd )
      Nrm(1) = 0.0_core_rknd
      call stat_update_var( iNrm, Nrm, stats_zt )

      ! Compute snow sedimentation velocity
      ! Using number concentration
      VNs(2:nz) = real( flip( real( uns_flip(1:nz-1),kind=dp ), nz-1 ), kind = core_rknd )
      VNs(1) = 0.0_core_rknd
      call stat_update_var( iVNs, VNs, stats_zm )
      ! Using mixing ratio
      Vrs(2:nz) = real( flip( real( ums_flip(1:nz-1),kind=dp ), nz-1 ), kind = core_rknd )
      Vrs(1) = 0.0_core_rknd
      call stat_update_var( iVrs, Vrs, stats_zm )

      ! Compute rain sedimentation velocity
      ! Using number concentration
      VNr(2:nz) = real( flip( real( unr_flip(1:nz-1),kind=dp ), nz-1 ), kind = core_rknd )
      VNr(1) = 0.0_core_rknd
      call stat_update_var( iVNr, VNr, stats_zm )
      ! Using mixing ratio
      Vrr(2:nz) = real( flip( real( umr_flip(1:nz-1),kind=dp ), nz-1 ), kind = core_rknd )
      Vrr(1) = 0.0_core_rknd
      call stat_update_var( iVrr, Vrr, stats_zm )

      ! Compute ice sedimentation velocity
      ! Using number concentration
      VNi(:) = real( uni, kind = core_rknd )
      call stat_update_var( iVNi, VNi, stats_zm )
      ! Using mixing ratio
      Vri(:) = real( umi, kind = core_rknd )
      call stat_update_var( iVri, Vri, stats_zm )

      ! Compute cloud droplet sedimentation
      ! Using number concentration
      VNc(:) = real( unc, kind = core_rknd )
      call stat_update_var( iVNc, VNc, stats_zm )
      ! Using mixing ratio
      Vrc(:) = real( umc, kind = core_rknd )
      call stat_update_var( iVrc, Vrc, stats_zm )

      ! Output sedimentation tendencies
      call stat_update_var( ircm_sd_mg_morr, rcm_sten, stats_zt )
      call stat_update_var( irim_sd_mg_morr, rim_sten, stats_zt ) 

      ! Sedimentation is handled within the MG microphysics
      hydromet_vel(:,:) = 0.0_core_rknd

      ! Compute Rain Water Path
      if ( irwp > 0 ) then

        xtmp = vertical_integral &
               ( (nz - 2 + 1), rho(2:nz), &
                 rrm(2:nz), dzt(2:nz) )

        call stat_update_var_pt( irwp, 1, xtmp, stats_sfc )

      end if

      ! Compute Snow Water Path
      if ( iswp > 0 ) then

        xtmp = vertical_integral &
               ( (nz - 2 + 1), rho(2:nz), &
                 rsm(2:nz), dzt(2:nz) )

        call stat_update_var_pt( iswp, 1, xtmp, stats_sfc )

      end if

    end if ! l_stats_samp

    return
  end subroutine mg_microphys_driver

end module mg_microphys_driver_module
