! $Id$
module mg_micro_driver_module

  implicit none

  public :: mg_microphys_driver

  private

  contains
!-------------------------------------------------------------------------------
  subroutine mg_microphys_driver &
             ( dt, nz, l_stats_samp, l_local_kk, l_latin_hypercube, &
               thlm, p_in_Pa, exner, rho, pdf_params, &
               rcm, rvm, Ncnm, hydromet, hydromet_mc, &
               hydromet_vel, rcm_mc, rvm_mc, thlm_mc, wp2_zt)
! Description:
!   Wrapper for the Morrison-Gettelman microphysics
! References:
!   None
!-------------------------------------------------------------------------------

    use pdf_parameter_module, only: pdf_parameter

    use parameters_model, only: hydromet_dim

    use T_in_K_module, only: &
      T_in_K2thlm, & ! Procedure(s)
      thlm2T_in_K

    use array_index, only:  &
      iiricem, iiNim, iiNcm

    use constants_clubb, only: &
      zero_threshold, &
      T_freeze_K, &
      Cp
      
    use stats_variables, only: & 
      zt, &  ! Variables
      irsnowm, &
      ieff_rad_cloud, &
      ieff_rad_ice

    use stats_type, only:  & 
      stat_update_var
      
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
      
    use phys_buffer, only: & ! Used for placing wp2_zt in morrison-gettelman microphysics
      pbuf_add,            &
      pbuf_allocate

    use phys_buffer, only: &
      pbuf, &      ! Variable
      pbuf_setval  ! Procedure
      
    use grid_class, only: &
      flip ! Procedure
      
    use constants_clubb, only: &
      rc_tol, & ! Variable
      rt_tol

    use shr_kind_mod, only: r8 => shr_kind_r8 

    implicit none

    ! External
    intrinsic :: max

    ! Input Variables
    real, intent(in) :: dt ! Model timestep        [s]

    integer, intent(in) :: nz ! Points in the Vertical        [-]

    logical, intent(in) :: &
      l_stats_samp,     & ! Whether to accumulate statistics [T/F]
      l_local_kk,       & ! Whether we're using the local formulas
      l_latin_hypercube   ! Whether we're using latin hypercube sampling

    real, dimension(nz), intent(in) :: &
      thlm,       & ! Liquid potential temperature       [K]
      p_in_Pa,    & ! Pressure                           [Pa]
      exner,      & ! Exner function                     [-]
      rho           ! Density on thermo. grid            [kg/m^3]

    type(pdf_parameter), dimension(nz), intent(in) :: &
      pdf_params ! PDF parameters

    real, dimension(nz), intent(in) :: &
      rcm,      & ! Liquid water mixing ratio          [kg/kg]
      rvm,      & ! Vapor water mixing ratio           [kg/kg]
      Ncnm,     & ! Cloud nuclei number concentration  [count/m^3]
      wp2_zt      ! w'^2 on the thermo. grid           [m^2/s^2]

    real, dimension(nz,hydromet_dim), intent(in) :: &
      hydromet ! Hydrometeor species    [units vary]

    ! Input/Output Variables
    real, dimension(nz,hydromet_dim), intent(inout) :: &
      hydromet_mc,   & ! Hydrometeor time tendency          [(units vary)/s]
      hydromet_vel     ! Hydrometeor sedimentation velocity [m/s]

    ! Output Variables
    real, dimension(nz), intent(out) :: &
      rcm_mc, & ! Time tendency of liquid water mixing ratio    [kg/kg/s]
      rvm_mc, & ! Time tendency of vapor water mixing ratio     [kg/kg/s]
      thlm_mc   ! Time tendency of liquid potential temperature [K/s]

    ! Local Variables
    real, dimension(nz) :: & 
      T_in_K,       & ! Temperature                                                   [K]
      T_in_K_new,   & ! Temperature after microphysics                                [K]
      tlat,         & ! Latent heating rate                                           [W/kg]
      cloud_frac,   & ! Liquid Cloud fraction                                         [-]
      rcm_new,      & ! Cloud water mixing ratio after microphysics                   [kg/kg]
      rsnowm,       & ! Snow mixing ratio (not in hydromet, output straight to stats) [kg/kg]
      effc,         & ! Droplet effective radius                                      [μ]
      effi            ! cloud ice effective radius                                    [μ]

    real(r8), dimension(nz) :: & 
      turbtype_flip,& ! Turbulence type at each interface                             [-]
      smaw_flip,    & ! Normalized instability function of momentum                   [???]
      wsub_flip,    & ! Diagnosed sub-grid vertical velocity st. dev.                 [m/s]
      wsubi_flip      ! Diagnosed sub-grid vertical velocity ice                      [m/s]
      
    real(r8), dimension(1, nz-1, 1) :: &
      aer_mmr_flip

    ! Note that the MG grid is flipped with respect to CLUBB.
    ! For this reason, all MG variables are marked as "flip"

    ! MG Input Variables on zm grid
    real(r8), dimension(nz) :: &
      Kh_zm_flip,   & ! Eddy diffusivity coefficient on momentum levels      [m^2/s]
      em_flip         ! Turbulent Kinetic Energy (TKE)                       [m^2/s^2]

    ! MG Input Variables on zt grid
    real(r8), dimension(nz-1) :: &
      T_in_K_flip,  & ! Air temperature                                      [K]
      rvm_flip,     & ! Vapor water mixing ratio                             [kg/kg]
      rcm_flip,     & ! Cloud water mixing ratio                             [kg/kg]
      p_in_Pa_flip, & ! Pressure                                             [Pa]
      pdel_flip,    & ! difference in air pressure between vertical levels   [Pa]
      cldn_flip,    & ! Cloud fraction                                       [-]
      liqcldf_flip, & ! Liquid cloud fraction                                [-]
      icecldf_flip, & ! Ice cloud fraction                                   [-]
      naai_flip,    & ! number of activated ice nuclei                       [1/kg]
      !rho_flip,     & ! Density on thermo. grid                              [kg/m^3]
      npccn_flip,   & ! Number of cloud nuclei                               [count/m^3]
      unused_in       ! Represents MG variables that are not used in the current code
      
    real(r8), dimension(nz-1,4) :: &
      rndst_flip,   & ! radius of 4 dust bins for contact freezing [units unknown, our guess is mm]
      nacon_flip      ! number of 4 dust bins for contact nucleation [units unknown]
      
    ! MG Output Variables
    real(r8), dimension(nz-1) :: &
      tlat_flip,    & ! Latent heating rate                                  [W/kg]
      rcm_mc_flip,  & ! Time tendency of liquid water mixing ratio           [kg/kg/s]
      rvm_mc_flip,  & ! Time tendency of vapor water mixing ratio            [kg/kg/s]
      cldo_flip,    & ! Old cloud fraction.                                  [-]
      rsnowm_flip,  & ! snow mixing ratio                                    [kg/kg]
      effc_flip,    & ! Droplet effective radius                             [μ]
      effi_flip       ! cloud ice effective radius                           [μ]
      
    ! MG zt variables that are not used in CLUBB.
    real(r8), dimension(nz-1) :: &
      rate1ord_cw2pr_st_flip, effc_fn_flip, prect_flip, preci_flip,            &
      nevapr_flip, evapsnow_flip, prain_flip, prodsnow_flip, cmeout_flip,      &
      deffi_flip, pgamrad_flip, lamcrad_flip, dsout_flip, qcsevap_flip,        &
      qisevap_flip, qvres_flip, cmeiout_flip, vtrmc_flip, vtrmi_flip,          &
      qcsedten_flip, qisedten_flip, prao_flip, prco_flip, mnuccco_flip,        &
      mnuccto_flip, msacwio_flip, psacwso_flip, bergso_flip, bergo_flip,       &
      melto_flip, homoo_flip, qcreso_flip, prcio_flip, praio_flip, qireso_flip,&
      mnuccro_flip, pracso_flip, meltsdt_flip, frzrdt_flip, mnuccdo_flip,      &
      qrout_flip, reff_rain_flip, reff_snow_flip, naai_hom_flip, dum

    ! MG zm variables that are not used in CLUBB
    real(r8), dimension(nz) :: &
      rflx_flip, sflx_flip
      
    real(r8), dimension(nz-1,hydromet_dim) :: &
       hydromet_flip, &  ! Hydrometeor species                               [units vary]
       hydromet_mc_flip  ! Hydrometeor time tendency                         [units vary]

    integer :: i

    ! ---- Begin Code ----
    ! Some dummy assignments to make compiler warnings go away...
    if ( .false. ) then
      if ( l_local_kk ) stop
    end if
    
    unused_in(:) = -9999.9_r8
    
    ! Initialize output arrays to zero
    naai_flip = 0.0_r8
    npccn_flip = 0.0_r8
    rndst_flip = 0.0_r8
    nacon_flip = 0.0_r8
    effc(:) = 0.0
    effi(:) = 0.0
    rsnowm(:) = 0.0
    tlat(:) = 0.0
    rcm_new(:) = 0.0
    T_in_K_new(:) = 0.0
    rcm_mc(1:nz) = 0.0
    rvm_mc(1:nz) = 0.0
    hydromet_mc(1:nz,:) = 0.0
    hydromet_mc_flip(1:nz-1,:) = 0.0_r8
    hydromet_vel(:,:) = 0.0

    ! Determine temperature
    T_in_K = thlm2T_in_K( thlm, exner, rcm )
    
    if ( l_latin_hypercube ) then
      ! Don't use sgs cloud fraction to weight the tendencies
      cloud_frac(1:nz) = 0.0

    else 
      ! Use sgs cloud fraction to weight tendencies
      cloud_frac(1:nz) = max( zero_threshold, &
                        pdf_params%mixt_frac * pdf_params%cloud_frac1 &
                        + (1.-pdf_params%mixt_frac) * pdf_params%cloud_frac2 )
    end if
    
    ! MG's grid is flipped with respect to CLUBB.
    ! Flip CLUBB variables before inputting them into MG.
    Kh_zm_flip(1:nz) = real( flip( dble(Kh_zm(1:nz) ), nz ), kind=r8 )
    em_flip(1:nz) = real( flip( dble(em(1:nz) ), nz ), kind=r8 )

    ! MG does not include zt(1) as CLUBB does, so it is removed during the flip
    T_in_K_flip(1:nz-1) = real( flip( dble(T_in_K(2:nz) ), nz-1 ), kind=r8 )
    rvm_flip(1:nz-1) = real( flip( dble(rvm(2:nz) ), nz-1 ), kind=r8 )
    rcm_flip(1:nz-1) = real( flip( dble(rcm(2:nz) ), nz-1 ), kind=r8 )
    p_in_Pa_flip(1:nz-1) = real( flip( dble(p_in_Pa(2:nz) ), nz-1 ), kind=r8 )
    liqcldf_flip(1:nz-1) = real( flip( dble(cloud_frac(2:nz) ), nz-1 ), kind=r8 )
    
    ! Hydromet is 2 dimensional, so flip function doesn't work
    do i = 1, hydromet_dim, 1
      hydromet_flip(1:nz-1, i) = real( hydromet(nz:2:-1, i), kind=r8 )
    end do
    
    ! Initialize MG input variables. The top level is skipped because MG's grid is flipped with
    ! respect to CLUBB's, and MG doesn't include CLUBB's lowest level.
    do i = 1, nz-1, 1
    
      ! Find difference in air pressure between vertical levels
      if ( i /= nz-1 ) then
        pdel_flip(i) = p_in_Pa_flip(i) - p_in_Pa_flip(i+1)
      else
        pdel_flip(i) = p_in_Pa_flip(i)
      end if
      
      ! Cloud fraction. In MG there is no difference between ice cloud fraction and
      ! liquid cloud fraction. However, just setting icecldf equal to CLUBBs cloud_frac
      ! won't work for all ice, no liquid clouds.
      if ( rcm_flip(i) > 0._r8 ) then
        cldn_flip(i) = liqcldf_flip(i)
      else if ( rcm_flip(i) < real( rc_tol, kind=r8 ) &
                .and. hydromet_flip(i,iiricem) > real( rc_tol, kind=r8 ) ) then
        cldn_flip(i) = 1._r8
      else
        cldn_flip(i) = 0._r8
      end if
      liqcldf_flip(i) = cldn_flip(i)
      icecldf_flip(i) = cldn_flip(i)

    end do
    
    ! Initialize grid variables. These are imported in the MG code.
    call init_ppgrid( nz )

    ! Place wp2 into the dummy phys_buffer module to import it into microp_aero_ts.
    call pbuf_add( 'WP2', 1, nz, 1 )
    call pbuf_allocate()
    call pbuf_setval( 'WP2', real( wp2_zt, kind=r8 ) )
    
    call gestbl(173.16_r8, 375.16_r8, 20.00_r8, .true., real( epsilo, kind=r8), &
                 real( latvap, kind=r8), real( latice, kind=r8), real( rh2o, kind=r8), &
                 real( cpair, kind=r8), real( tmelt, kind=r8) ) ! Known magic flag

    ! Ensure no hydromet arrays are input as 0, because this makes MG crash.
    do i=1, hydromet_dim, 1
      hydromet_flip(:,i) = max( 1e-8_r8, hydromet_flip(:,i) ) ! Known magic number
    end do
    
    ! Prescribe droplet concentration
    hydromet_flip(:,iiNcm) = 1e8_r8
    
    !---------------------------------------------------------------------------------------
    ! This code block initializes outputs of the microp_aero_ts subroutine in case we
    ! decide not to use it. If microp_aero_ts is called, this code block can be commented
    ! out.
    !
    !rho_flip(1:nz-1) = real( flip( dble(rho(2:nz) ), nz-1 ) )
    !npccn_flip(1:nz-1) = real( flip( dble(Ncnm(2:nz) ), nz-1 ) )
    !
    ! Determine ice nulceation number using Meyers formula found in the Morrison microphysics
    !naai_flip(i) = exp( -2.80 + 0.262 * ( T_freeze_K - T_in_K_flip(i) ) ) * 1000 / rho_flip(i)
    !
    ! Set values for radius of 4 dust bins for contact freezing. Currently, we assume
    ! nacon = 0 therby ommitting contact nucleation, but this is set regardless.
    !rndst_flip(:,1) = 0.001
    !rndst_flip(:,2) = 0.01
    !rndst_flip(:,3) = 0.1
    !rndst_flip(:,4) = 1
    !nacon_flip(:,1:4) = 0
    !----------------------------------------------------------------------------------------

    ! These variables are unused in CLUBB, so they are initialized to 1 before input into
    ! microp_aero_ts. 
    aer_mmr_flip(:,:,:) = 1._r8
    turbtype_flip(:) = 1._r8
    smaw_flip(:) = 1._r8
    
    ! Calculate aerosol activiation, dust size, and number for contact nucleation
    call microp_aero_ts &
         ( 1, 1, real( dt, kind=r8), T_in_K_flip, unused_in, &                                ! in
         rvm_flip, rcm_flip, hydromet_flip(:,iiricem), &                                      ! in
         hydromet_flip(:,iiNcm), hydromet_flip(:,iiNim), p_in_Pa_flip, pdel_flip, cldn_flip, &! in
         liqcldf_flip, icecldf_flip, &                                                        ! in
         cldo_flip, unused_in, unused_in, unused_in, unused_in, &                             ! in
         aer_mmr_flip, &
         pbuf, &
         Kh_zm_flip, em_flip, turbtype_flip, smaw_flip, wsub_flip, wsubi_flip, &
         naai_flip, naai_hom_flip, npccn_flip, rndst_flip, nacon_flip )                       ! out

    ! Call the Morrison-Gettelman microphysics
    call mmicro_pcond &
         (.false., &                                                                          ! in
         0, 1, real( dt, kind=r8), T_in_K_flip, &                                             ! in
         rvm_flip, rcm_flip, hydromet_flip(:,iiricem), &                                      ! in
         hydromet_flip(:,iiNcm), hydromet_flip(:,iiNim), p_in_Pa_flip, pdel_flip, cldn_flip, &! in
         liqcldf_flip, icecldf_flip, &                                                        ! in
         cldo_flip, &                                                                         ! in
         rate1ord_cw2pr_st_flip, &                                                            ! out
         naai_flip, npccn_flip, rndst_flip, nacon_flip, &                                     ! in
         tlat_flip, rvm_mc_flip, &                                                            ! out
         rcm_mc_flip, hydromet_mc_flip(:,iiricem), hydromet_mc_flip(:,iiNcm), &               ! out
         hydromet_mc_flip(:,iiNim), effc_flip, &                                              ! out
         effc_fn_flip, effi_flip, prect_flip, preci_flip,             &                       ! out
         nevapr_flip, evapsnow_flip,      &                                                   ! out
         prain_flip, prodsnow_flip, cmeout_flip, deffi_flip, pgamrad_flip, &                  ! out
         lamcrad_flip, rsnowm_flip, dsout_flip, &                                             ! out
         rflx_flip, sflx_flip, qrout_flip, reff_rain_flip, reff_snow_flip, &                  ! out
         qcsevap_flip, qisevap_flip, qvres_flip, cmeiout_flip, &                              ! out
         vtrmc_flip, vtrmi_flip, qcsedten_flip, qisedten_flip, &                              ! out
         prao_flip, prco_flip, mnuccco_flip, mnuccto_flip, msacwio_flip, psacwso_flip,&       ! out
         bergso_flip, bergo_flip, melto_flip, homoo_flip, qcreso_flip, prcio_flip, praio_flip, &
         qireso_flip,&                                                                        ! out
         mnuccro_flip, pracso_flip, meltsdt_flip, frzrdt_flip, mnuccdo_flip)                  ! out

    ! Flip MG variables into CLUBB grid
    rcm_mc(2:nz) = real( flip( dble(rcm_mc_flip(1:nz-1) ), nz-1 ) )
    rvm_mc(2:nz) = real( flip( dble(rvm_mc_flip(1:nz-1) ), nz-1 ) )
    rcm_new(2:nz) = real( flip( dble(rcm_flip(1:nz-1) ), nz-1 ) )
    effc(2:nz) = real( flip( dble(effc_flip(1:nz-1) ), nz-1 ) )
    effi(2:nz) = real( flip( dble(effi_flip(1:nz-1) ), nz-1 ) )
    tlat(2:nz) = real( flip( dble(tlat_flip(1:nz-1) ), nz-1 ) )
    rsnowm(2:nz) = real( flip( dble(rsnowm_flip(1:nz-1) ), nz-1 ) )
      
    do i = 1, hydromet_dim, 1      
      hydromet_mc(2:nz, i) = real( hydromet_mc_flip(nz-1:1:-1, i) )
    end do
    
    ! Update thetal based on absolute temperature
    T_in_K_new = T_in_K + (tlat/Cp) * real( dt )
    
    ! TODO: To remove compile warnings:
    dum = real( T_in_K_new(1:nz-1), kind=r8 )
    ! TODO: T_in_K is not changed within MG, so the change in temperature will need to be
    ! calculated another way. However, using the eqution above does not seem to work correctly.
    thlm_mc = ( T_in_K2thlm( T_in_K, exner, rcm_new ) - thlm ) / real( dt )
    
    if ( l_stats_samp ) then
      call stat_update_var( irsnowm, rsnowm, zt )
      
      ! Effective radii of hydrometeor species
      call stat_update_var( ieff_rad_cloud, effc(:), zt )
      call stat_update_var( ieff_rad_ice, effi(:), zt )
      
    end if ! l_stats_samp

    return
  end subroutine mg_microphys_driver

end module mg_micro_driver_module
