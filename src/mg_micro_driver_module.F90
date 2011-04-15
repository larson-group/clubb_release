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
      iiricem, iiNim, iiNcm, iirsnowm

    use constants_clubb, only: &
      zero_threshold, &
      T_freeze_K
      
    use cldwat2m_micro, only: &
      mmicro_pcond ! Procedure
      
    use wv_saturation, only: &
      gestbl ! Procedure
      
    use physconst, only: &
      cpair,  &
      rh2o,   &
      tmelt,  &
      latice, &
      latvap, &
      epsil
      
    use ppgrid, only: &
      init_ppgrid ! Procedure
      
    use grid_class, only: &
      flip ! Procedure
      
    use constants_clubb, only: &
      rc_tol ! Variable

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
      effc, effi ! Effective droplet radii [μ]

    real, dimension(nnzp) :: & 
      T_in_K,       & ! Temperature                             [K]
      cloud_frac,   & ! Liquid Cloud fraction                   [-]
      rcm_tmp         ! Temporary variable                      [kg/kg]
      
    real, dimension(nnzp,hydromet_dim) :: & 
      hydromet_tmp     ! Temporary variable
      
    ! MG Input Variables
    ! Note that the MG grid is flipped with respect to CLUBB.
    real, dimension(nnzp-1) :: &
      T_in_K_flip,  & ! Air temperature                                      [K]
      rvm_flip,     & ! Vapor water mixing ratio                             [kg/kg]
      rcm_flip,     & ! Cloud water mixing ratio                             [kg/kg]
      p_in_Pa_flip, & ! Pressure                                             [Pa]
      pdel_flip,    & ! difference in air pressure between vertical levels   [Pa]
      cldn_flip,    & ! Cloud fraction                                       [-]
      liqcldf_flip, & ! Liquid cloud fraction                                [-]
      icecldf_flip, & ! Ice cloud fraction                                   [-]
      naai_flip,    & ! number of activated ice nuclei                       [1/kg]
      rho_flip,     & ! Density on thermo. grid                              [kg/m^3]
      Ncnm_flip,    & ! Number of cloud nuclei                               [count/m^3]
      unused          ! Represents MG variables that are not used in the current code
      
    ! MG Output Variables
    ! Note that the MG grid is flipped with respect to CLUBB
    real, dimension(nnzp-1) :: &
      rcm_mc_flip     ! Time tendency of liquid water mixing ratio           [kg/kg/s]
      
    real, dimension(nnzp-1,hydromet_dim) :: &
       hydromet_flip, &  ! Hydrometeor species                               [units vary]
       hydromet_mc_flip  ! Hydrometeor time tendency                         [units vary]
      
    real, dimension(nnzp-1,4) :: &
      rndst_flip,   & ! radius of 4 dust bins for contact freezing [units unknown, our guess is mm]
      nacon_flip      ! number of 4 dust bins for contact nucleation [units unknown]

    integer :: i
    
    ! TODO temp output variables
    real, dimension(nnzp-1) :: &
      cldo, & ! Old cloud fraction
      rate1ord_cw2pr_st, & ! first order rate for direct cw to precip conversion
      tlat, &           !latent heating rate from microphysics					W kg-1
      qvlat, &  !qv tendency from microphysics						kg/kg/s
      effc_fn, & !	droplet effective radius, assuming nc = 1.e8 kg-1				micron
      prect, &  !surface precip rate							m s-1
      preci, &  !cloud ice/snow surface precip rate					m s-1
      nevapr, & !	evaporation rate of rain+snow						s-1
      evapsnow, &!	sublimation rate of snow						s-1
      prain, &  !production of rain+snow						s-1
      prodsnow, &!	production rate of snow						s-1
      cmeout, & !	dep/sublimation of cloud ice						s-1
      deffi, &  !effective ice diameter for Mitchell optics				micron
      pgamrad, & ! 	droplet size distribution shape parameter
      lamcrad, & !slope of droplet size distribution					m-1
      dsout, &  !snow effective diameter for Mitchell optics (after bug fix…)		micron
      qcsevap, & !	evaporation rate of droplets during sedimentation			s-1
      qisevap, & !		sublimation rate of ice during sedimentation				s-1
      qvres, &  !cond rate due to residual removal of excess supersaturation		s-1
      cmeiout, & !	cloud ice sub/dep rate	(same as cmeout)				s-1
      qcsedten, &!	tendency of cloud water due to sedimentation				s-1
      qisedten, & !	tendency of cloud ice due to sedimentation				s-1
      prao, & ! accretion of cloud by rain 
      prco, & ! autoconversion of cloud to rain
      mnuccco, & ! mixing rat tend due to immersion freezing
      mnuccto, & ! mixing ratio tend due to contact freezing
      msacwio, & ! mixing ratio tend due to H-M splintering
      psacwso, & ! collection of cloud water by snow
      melto, & ! melting of cloud ice
      homoo, & ! homogeneos freezign cloud water
      qcreso, & ! residual cloud condensation due to removal of excess supersat
      prcio, & ! autoconversion of cloud ice to snow
      praio, & ! accretion of cloud ice by snow
      qireso, & ! residual ice deposition due to removal of excess supersat
      mnuccro, & ! mixing ratio tendency due to heterogeneous freezing of rain to snow (1/s)
      pracso , & ! mixing ratio tendency due to accretion of rain by snow (1/s)
      meltsdt, & ! latent heating rate due to melting of snow  (W/kg)
      frzrdt  ! latent heating rate due to homogeneous freezing of rain (W/kg


    ! ---- Begin Code ----
    ! Some dummy assignments to make compiler warnings go away...
    if ( .false. ) then
      if ( l_local_kk ) stop
    end if

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

    ! Initialize tendencies to zero
    rcm_mc(1:nnzp) = 0.0
    rvm_mc(1:nnzp) = 0.0
    hydromet_mc(1:nnzp,:) = 0.0
    
    ! MG's grid is flipped with respect to CLUBB.
    ! Flip CLUBB variables before inputting them into MG.
    ! In addition, MG does not include CLUBB's bottom level
    T_in_K_flip(1:nnzp-1) = real( flip( dble(T_in_K(2:nnzp) ), nnzp-1 ) )
    rvm_flip(1:nnzp-1) = real( flip( dble(rvm(2:nnzp) ), nnzp-1 ) )
    rcm_flip(1:nnzp-1) = real( flip( dble(rcm(2:nnzp) ), nnzp-1 ) )
    p_in_Pa_flip(1:nnzp-1) = real( flip( dble(p_in_Pa(2:nnzp) ), nnzp-1 ) )
    liqcldf_flip(1:nnzp-1) = real( flip( dble(cloud_frac(2:nnzp) ), nnzp-1 ) )
    rho_flip(1:nnzp-1) = real( flip( dble(rho(2:nnzp) ), nnzp-1 ) )
    Ncnm_flip(1:nnzp-1) = real( flip( dble(Ncnm(2:nnzp) ), nnzp-1 ) )
    
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

      ! Determine ice nulceation number using Meyers formula found in the Morrison microphysics
      naai_flip(i) = exp( -2.80 + 0.262 * ( T_freeze_K - T_in_K_flip(i) ) ) * 1000 / rho_flip(i)
      
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
      
      ! Set values for radius of 4 dust bins for contact freezing. Currently, we assume
      ! nacon = 0 therby ommitting contact nucleation, but this is set regardless.
      rndst_flip(i,1) = 0.001
      rndst_flip(i,2) = 0.01
      rndst_flip(i,3) = 0.1
      rndst_flip(i,4) = 1
      
      nacon_flip(i,1:4) = 0
      
      ! Initialize unused MG variables
      unused(i) = -9999.9
      
    end do
    
    ! Initialize grid variables. These are imported in the MG code.
    call init_ppgrid()
    
    call gestbl(173.16, 375.16, 20.00, .true., epsil, &
                 latvap, latice, rh2o, cpair, tmelt)

    ! Call the Morrison-Gettelman microphysics
    call mmicro_pcond &
         ( 0, 1, dt, T_in_K_flip, unused, & ! Input
         rvm_flip, unused, unused, rcm_flip, & !in
         hydromet_flip(:,iiricem), hydromet_flip(:,iiNcm), hydromet_flip(:,iiNim), & !in/out
         p_in_Pa_flip, pdel_flip, cldn_flip, & !in
         liqcldf_flip, icecldf_flip, &
         cldo, unused, unused, unused, unused, & !output
         rate1ord_cw2pr_st, & !out
         naai_flip, Ncnm_flip, rndst_flip, nacon_flip, & !in
         unused, unused, unused, &
         tlat, qvlat, & ! Output
         rcm_mc_flip, hydromet_mc_flip(:,iiricem), hydromet_mc_flip(:,iiNcm), &
         hydromet_mc_flip(:,iiNim), effc, &
         effc_fn, effi, prect, preci,             &  
         nevapr, evapsnow,      &
         prain, prodsnow, cmeout, deffi, pgamrad, &
         lamcrad, hydromet_flip(:,iirsnowm), dsout, &
         qcsevap,qisevap,qvres,cmeiout, &
         unused,unused,qcsedten,qisedten, &
         prao,prco,mnuccco,mnuccto,msacwio,psacwso,&
         unused,unused,melto,homoo,qcreso,prcio,praio,qireso,&
         mnuccro,pracso,meltsdt,frzrdt )

    rcm_mc(2:nnzp) = real( flip( dble(rcm_mc_flip(1:nnzp-1) ), nnzp-1 ) )
    T_in_K(2:nnzp) = real( flip( dble(T_in_K_flip(1:nnzp-1) ), nnzp-1 ) )
    rcm_tmp(2:nnzp) = real( flip( dble(rcm_flip(1:nnzp-1) ), nnzp-1 ) )
      
    do i = 1, hydromet_dim, 1
      ! MG doesn't take CLUBB's lowest level into account
      hydromet_tmp(1, i) = 0
      hydromet_mc(1, i) = 0
      
      hydromet_tmp(2:nnzp, i) = hydromet_flip(nnzp-1:1:-1, i)
      hydromet_mc(2:nnzp, i) = hydromet_mc_flip(nnzp-1:1:-1, i)
    end do
    

    ! Update thetal based on absolute temperature
    thlm_mc = ( T_in_K2thlm( T_in_K, exner, rcm_tmp ) - thlm ) / real( dt )
    
    ! Sedimentation is handled within the MG microphysics
    hydromet_vel(:,:) = 0.0
    
    if ( l_stats_samp ) then
    end if ! l_stats_samp

    return
  end subroutine mg_microphys_driver

end module mg_micro_driver_module
