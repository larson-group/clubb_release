! $Id$
module mg_micro_driver_module

  implicit none

  public :: mg_microphys_driver

  private

  contains
!-------------------------------------------------------------------------------
  subroutine mg_microphys_driver &
             ( dt, nnzp, l_stats_samp, l_local_kk, l_latin_hypercube, &
               thlm, p_in_Pa, exner, rho, pdf_params, &
               rcm, rvm, Ncnm, hydromet, hydromet_mc, &
               hydromet_vel, rcm_mc, rvm_mc, thlm_mc)
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
      
    use grid_class, only: &
      flip ! Procedure
      
    use constants_clubb, only: &
      rc_tol ! Variable

    use shr_kind_mod, only: r8 => shr_kind_r8 

    implicit none

    ! External
    intrinsic :: max

    ! Input Variables
    real, intent(in) :: dt ! Model timestep        [s]

    integer, intent(in) :: nnzp ! Points in the Vertical        [-]

    logical, intent(in) :: &
      l_stats_samp,     & ! Whether to accumulate statistics [T/F]
      l_local_kk,       & ! Whether we're using the local formulas
      l_latin_hypercube   ! Whether we're using latin hypercube sampling

    real, dimension(nnzp), intent(in) :: &
      thlm,       & ! Liquid potential temperature       [K]
      p_in_Pa,    & ! Pressure                           [Pa]
      exner,      & ! Exner function                     [-]
      rho           ! Density on thermo. grid            [kg/m^3]

    type(pdf_parameter), dimension(nnzp), intent(in) :: &
      pdf_params ! PDF parameters

    real, dimension(nnzp), intent(in) :: &
      rcm,      & ! Liquid water mixing ratio          [kg/kg]
      rvm,      & ! Vapor water mixing ratio           [kg/kg]
      Ncnm        ! Cloud nuclei number concentration  [count/m^3]

    real, dimension(nnzp,hydromet_dim), intent(in) :: &
      hydromet ! Hydrometeor species    [units vary]

    ! Input/Output Variables
    real, dimension(nnzp,hydromet_dim), intent(inout) :: &
      hydromet_mc,   & ! Hydrometeor time tendency          [(units vary)/s]
      hydromet_vel     ! Hydrometeor sedimentation velocity [m/s]

    ! Output Variables
    real, dimension(nnzp), intent(out) :: &
      rcm_mc, & ! Time tendency of liquid water mixing ratio    [kg/kg/s]
      rvm_mc, & ! Time tendency of vapor water mixing ratio     [kg/kg/s]
      thlm_mc   ! Time tendency of liquid potential temperature [K/s]

    ! Local Variables
    real, dimension(nnzp) :: & 
      T_in_K,       & ! Temperature                                                   [K]
      T_in_K_new,   & ! Temperature after microphysics                                [K]
      tlat,         & ! Latent heating rate                                           [W/kg]
      cloud_frac,   & ! Liquid Cloud fraction                                         [-]
      rcm_new,      & ! Cloud water mixing ratio after microphysics                   [kg/kg]
      rsnowm,       & ! Snow mixing ratio (not in hydromet, output straight to stats) [kg/kg]
      effc,         & ! Droplet effective radius                                      [μ]
      effi            ! cloud ice effective radius                                    [μ]

    real(r8), dimension(nnzp) :: & 
      turbtype_flip,& ! Turbulence type at each interface                             [-]
      smaw_flip,    & ! Normalized instability function of momentum                   [???]
      wsub_flip,    & ! Diagnosed sub-grid vertical velocity st. dev.                 [m/s]
      wsubi_flip      ! Diagnosed sub-grid vertical velocity ice                      [m/s]
      
    real(r8), dimension(1, nnzp-1, 0) :: &
      aer_mmr_flip

    ! MG Input Variables
    ! Note that the MG grid is flipped with respect to CLUBB.
    real(r8), dimension(nnzp-1) :: &
      Kh_zm_flip,   & ! Eddy diffusivity coefficient on momentum levels      [m^2/s]
      em_flip,      & ! Turbulent Kinetic Energy (TKE)                       [m^2/s^2]
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
      
    real(r8), dimension(nnzp-1,4) :: &
      rndst_flip,   & ! radius of 4 dust bins for contact freezing [units unknown, our guess is mm]
      nacon_flip      ! number of 4 dust bins for contact nucleation [units unknown]
      
    ! MG Output Variables
    ! Note that the MG grid is flipped with respect to CLUBB
    real(r8), dimension(nnzp-1) :: &
      tlat_flip,    & ! Latent heating rate                                  [W/kg]
      rcm_mc_flip,  & ! Time tendency of liquid water mixing ratio           [kg/kg/s]
      rvm_mc_flip,  & ! Time tendency of vapor water mixing ratio            [kg/kg/s]
      cldo_flip,    & ! Old cloud fraction.                                  [-]
      rsnowm_flip,  & ! snow mixing ratio                                    [kg/kg]
      effc_flip,    & ! Droplet effective radius                             [μ]
      effi_flip       ! cloud ice effective radius                           [μ]
      
    ! MG variables that are not used in CLUBB.
    ! A separate variables needs to be passed in for every intent_out variable in MG,
    ! otherwise Fortran treats all the variables are if they point to the same memory
    ! location.
    real(r8), dimension(nnzp-1) :: &
      unused_out01, unused_out02, unused_out03, &
      unused_out04, unused_out05, unused_out06, &
      unused_out07, unused_out08, unused_out09, &
      unused_out10, unused_out11, unused_out12, &
      unused_out13, unused_out14, unused_out15, &
      unused_out16, unused_out17, unused_out18, &
      unused_out19, unused_out20, unused_out21, &
      unused_out22, unused_out23, unused_out24, &
      unused_out25, unused_out26, unused_out27, &
      unused_out28, unused_out29, unused_out30, &
      unused_out31, unused_out32, unused_out33, &
      unused_out34, unused_out35, unused_out36, &
      unused_out37, unused_out38, unused_out39
      
    real(r8), dimension(nnzp-1,hydromet_dim) :: &
       hydromet_flip, &  ! Hydrometeor species                               [units vary]
       hydromet_mc_flip  ! Hydrometeor time tendency                         [units vary]

    integer :: i

    ! ---- Begin Code ----
    ! Some dummy assignments to make compiler warnings go away...
    if ( .false. ) then
      if ( l_local_kk ) stop
    end if
    
    unused_in(:) = -9999.9
    
    ! Initialize output arrays to zero
    naai_flip = 0.0_r8
    npccn_flip = 0.0_r8
    rndst_flip = 0.0_r8
    nacon_flip = 0.0_r8
    effc(:) = 0.0_r8
    effi(:) = 0.0_r8
    rsnowm(:) = 0.0_r8
    tlat(:) = 0.0_r8
    rcm_new(:) = 0.0_r8
    T_in_K_new(:) = 0.0_r8
    rcm_mc(1:nnzp) = 0.0_r8
    rvm_mc(1:nnzp) = 0.0_r8
    hydromet_mc(1:nnzp,:) = 0.0_r8
    hydromet_mc_flip(1:nnzp-1,:) = 0.0_r8
    hydromet_vel(:,:) = 0.0_r8

    ! Determine temperature
    T_in_K = thlm2T_in_K( thlm, exner, rcm )
    
    if ( l_latin_hypercube ) then
      ! Don't use sgs cloud fraction to weight the tendencies
      cloud_frac(1:nnzp) = 0.0

    else 
      ! Use sgs cloud fraction to weight tendencies
      cloud_frac(1:nnzp) = max( zero_threshold, &
                        pdf_params%mixt_frac * pdf_params%cloud_frac1 &
                        + (1.-pdf_params%mixt_frac) * pdf_params%cloud_frac2 )
    end if
    
    ! MG's grid is flipped with respect to CLUBB.
    ! Flip CLUBB variables before inputting them into MG.
    ! In addition, MG does not include CLUBB's bottom level
    Kh_zm_flip(1:nnzp-1) = real( flip( dble(Kh_zm(2:nnzp) ), nnzp-1 ) )
    em_flip(1:nnzp-1) = real( flip( dble(em(2:nnzp) ), nnzp-1 ) )
    T_in_K_flip(1:nnzp-1) = real( flip( dble(T_in_K(2:nnzp) ), nnzp-1 ) )
    rvm_flip(1:nnzp-1) = real( flip( dble(rvm(2:nnzp) ), nnzp-1 ) )
    rcm_flip(1:nnzp-1) = real( flip( dble(rcm(2:nnzp) ), nnzp-1 ) )
    p_in_Pa_flip(1:nnzp-1) = real( flip( dble(p_in_Pa(2:nnzp) ), nnzp-1 ) )
    liqcldf_flip(1:nnzp-1) = real( flip( dble(cloud_frac(2:nnzp) ), nnzp-1 ) )
    
    ! Hydromet is 2 dimensional, so flip function doesn't work
    do i = 1, hydromet_dim, 1
      hydromet_flip(1:nnzp-1, i) = hydromet(nnzp:2:-1, i)
    end do
    
    ! Initialize MG input variables. The top level is skipped because MG's grid is flipped with
    ! respect to CLUBB's, and MG doesn't include CLUBB's lowest level.
    do i = 1, nnzp-1, 1
    
      ! Find difference in air pressure between vertical levels
      if ( i /= nnzp-1 ) then
        pdel_flip(i) = p_in_Pa_flip(i) - p_in_Pa_flip(i+1)
      else
        pdel_flip(i) = p_in_Pa_flip(i)
      end if
      
      ! Cloud fraction. In MG there is no difference between ice cloud fraction and
      ! liquid cloud fraction. However, just setting icecldf equal to CLUBBs cloud_frac
      ! won't work for all ice, no liquid clouds.
      if ( rcm_flip(i) > 0 ) then
        cldn_flip(i) = liqcldf_flip(i)
      else if ( rcm_flip(i) < rc_tol .and. hydromet_flip(i,iiricem) > rc_tol ) then
        cldn_flip(i) = 1
      else
        cldn_flip(i) = 0
      end if
      liqcldf_flip(i) = cldn_flip(i)
      icecldf_flip(i) = cldn_flip(i)

    end do
    
    ! Initialize grid variables. These are imported in the MG code.
    call init_ppgrid()
    
    call gestbl(173.16_r8, 375.16_r8, 20.00_r8, .true., real( epsilo, kind=r8), &
                 real( latvap, kind=r8), real( latice, kind=r8), real( rh2o, kind=r8), &
                 real( cpair, kind=r8), real( tmelt, kind=r8) )
                 
    ! Ensure no hydromet arrays are input as 0, because this makes MG crash.
    do i=1, hydromet_dim, 1
      hydromet_flip(:,i) = max(1e-8, hydromet_flip(:,i))
    end do
    
    ! Prescribe droplet concentration
    hydromet_flip(:,iiNcm) = 1e8
    
    !---------------------------------------------------------------------------------------
    ! This code block initializes outputs of the microp_aero_ts subroutine in case we
    ! decide not to use it. If microp_aero_ts is called, this code block can be commented
    ! out.
    !
    !rho_flip(1:nnzp-1) = real( flip( dble(rho(2:nnzp) ), nnzp-1 ) )
    !npccn_flip(1:nnzp-1) = real( flip( dble(Ncnm(2:nnzp) ), nnzp-1 ) )
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
    aer_mmr_flip(:,:,:) = 1
    turbtype_flip(:) = 1
    smaw_flip(:) = 1
    
    ! Calculate aerosol activiation, dust size, and number for contact nucleation
    call microp_aero_ts &
         ( 0, 1, real( dt, kind=r8), T_in_K_flip, unused_in, &                                ! in
         rvm_flip, rcm_flip, hydromet_flip(:,iiricem), &                                      ! in
         hydromet_flip(:,iiNcm), hydromet_flip(:,iiNim), p_in_Pa_flip, pdel_flip, cldn_flip, &! in
         liqcldf_flip, icecldf_flip, &                                                        ! in
         cldo_flip, unused_in, unused_in, unused_in, unused_in, &                             ! in
         aer_mmr_flip, &
         Kh_zm_flip, em_flip, turbtype_flip, smaw_flip, wsub_flip, wsubi_flip, &
         naai_flip, npccn_flip, rndst_flip, nacon_flip )                                     ! out

    ! Call the Morrison-Gettelman microphysics
    call mmicro_pcond &
         ( 0, 1, real( dt, kind=r8), T_in_K_flip, unused_in, &                                ! in
         rvm_flip, unused_in, unused_in, rcm_flip, hydromet_flip(:,iiricem), &                ! in
         hydromet_flip(:,iiNcm), hydromet_flip(:,iiNim), p_in_Pa_flip, pdel_flip, cldn_flip, &! in
         liqcldf_flip, icecldf_flip, &                                                        ! in
         cldo_flip, unused_in, unused_in, unused_in, unused_in, &                             ! in
         unused_out01, &
         naai_flip, npccn_flip, rndst_flip, nacon_flip, &                                     ! in
         unused_in, unused_in, unused_in, tlat_flip, rvm_mc_flip, &
         rcm_mc_flip, hydromet_mc_flip(:,iiricem), hydromet_mc_flip(:,iiNcm), &               ! out
         hydromet_mc_flip(:,iiNim), effc_flip, &                                              ! out
         unused_out02, effi_flip, unused_out03, unused_out04,             &
         unused_out05, unused_out06,      &
         unused_out07, unused_out08, unused_out09, unused_out10, unused_out11, &
         unused_out12, rsnowm_flip, unused_out13, & !out
         unused_out14,unused_out15,unused_out16,unused_out17, &
         unused_out18, unused_out19, unused_out20, unused_out21, &
         unused_out22,unused_out23,unused_out24,unused_out25,unused_out26,unused_out27,& 
         unused_out28,unused_out29,unused_out30,unused_out31,unused_out32,unused_out33, &
         unused_out34,unused_out35,&
         unused_out36,unused_out37,unused_out38,unused_out39 )

    ! Flip MG variables into CLUBB grid
    rcm_mc(2:nnzp) = real( flip( dble(rcm_mc_flip(1:nnzp-1) ), nnzp-1 ) )
    rvm_mc(2:nnzp) = real( flip( dble(rvm_mc_flip(1:nnzp-1) ), nnzp-1 ) )
    rcm_new(2:nnzp) = real( flip( dble(rcm_flip(1:nnzp-1) ), nnzp-1 ) )
    effc(2:nnzp) = real( flip( dble(effc_flip(1:nnzp-1) ), nnzp-1 ) )
    effi(2:nnzp) = real( flip( dble(effi_flip(1:nnzp-1) ), nnzp-1 ) )
    tlat(2:nnzp) = real( flip( dble(tlat_flip(1:nnzp-1) ), nnzp-1 ) )
    rsnowm(2:nnzp) = real( flip( dble(rsnowm_flip(1:nnzp-1) ), nnzp-1 ) )
      
    do i = 1, hydromet_dim, 1      
      hydromet_mc(2:nnzp, i) = real( hydromet_mc_flip(nnzp-1:1:-1, i) )
    end do
    
    ! Update thetal based on absolute temperature
    T_in_K_new = T_in_K + (tlat/Cp) * real( dt )
    
    ! TODO: To remove compile warnings:
    unused_out01 = T_in_K_new(1:nnzp-1)
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
