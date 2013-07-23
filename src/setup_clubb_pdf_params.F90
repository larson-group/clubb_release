! $Id$
!===============================================================================
module setup_clubb_pdf_params

  implicit none

  private

  public :: setup_pdf_parameters, &
            unpack_pdf_params,    &
            normalize_pdf_params, &
            compute_mean_stdev,   &
            compute_corr

  private :: component_means_rain,   &
             precip_fraction,        &
             pdf_param_hm_stats,     &
             pdf_param_log_hm_stats, &
             pack_pdf_params

  contains

  !=============================================================================
  subroutine setup_pdf_parameters( nz, rrainm, Nrm, Ncnm, rho, rcm, &       ! Intent(in)
                                   cloud_frac, w_std_dev, wphydrometp, &    ! Intent(in)
                                   corr_array_cloud, corr_array_below, &    ! Intent(in)
                                   pdf_params, l_stats_samp, d_variables, & ! Intent(in)
                                   corr_array_1, corr_array_2, &            ! Intent(inout)
                                   mu_x_1, mu_x_2, sigma_x_1, sigma_x_2, &  ! Intent(out)
                                   hydromet_pdf_params )                    ! Intent(out)

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        zm2zt, & ! Procedure(s)
        gr

    use constants_clubb, only: &
        one,     & ! Constant(s)
        zero,    &
        rc_tol,  &
        rr_tol,  &
        Nr_tol,  &
        Ncn_tol, &
        eps

    use pdf_parameter_module, only: &
        pdf_parameter  ! Variable(s)

    use hydromet_pdf_parameter_module, only: &
        hydromet_pdf_parameter  ! Variable(s)

    use parameters_model, only: &
        hydromet_dim  ! Variable(s)

    use array_index, only: &
        iirrainm, & ! Constant(s)
        iiNrm

    use model_flags, only: &
        l_use_precip_frac, & ! Flag(s)
        l_calc_w_corr, &
        l_use_modified_corr

    use advance_windm_edsclrm_module, only: &
        xpwp_fnc

    use variables_diagnostic_module, only: &
        Kh_zm

    use parameters_tunable, only: &
        c_Krrainm

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    use parameters_microphys, only: &
        rrp2_on_rrm2_cloud,     & ! Variable(s)
        rrp2_on_rrm2_below,     &
        Nrp2_on_Nrm2_cloud,     &
        Nrp2_on_Nrm2_below,     &
        Ncnp2_on_Ncnm2_cloud

    use stats_type, only: &
        stat_update_var ! Procedure(s)

    use stats_variables, only: &
        irr1,             & ! Variable(s)
        irr2,             &
        iNr1,             &
        iNr2,             &
        iprecip_frac,     &
        iprecip_frac_1,   &
        iprecip_frac_2,   &
        zt

    use model_flags, only: &
        l_diagnose_correlations ! Variable(s)

    use diagnose_correlations_module, only: &
        diagnose_correlations, & ! Procedure(s)
        setup_corr_cholesky_mtx, &
        cholesky_to_corr_mtx_approx, &
        rearrange_corr_array, &
        corr_array_assertion_checks

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz,          & ! Number of model vertical grid levels
      d_variables    ! Number of variables in the correlation array

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      rrainm,     & ! Mean rain water mixing ratio, < r_r >           [kg/kg]
      Nrm,        & ! Mean rain drop concentration, < N_r >           [num/kg]
      Ncnm,       & ! Mean cloud nuclei concentration, < N_cn >       [num/kg]
      rho,        & ! Density                                         [kg/m^3]
      rcm,        & ! Mean cloud water mixing ratio, < r_c >          [kg/kg]
      cloud_frac, & ! Cloud fraction                                  [-]
      w_std_dev     ! Standard deviation of vertical velocity, w      [m/s]

    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(in) :: &
      wphydrometp    ! Covariance < w'h_m' > (momentum levels)     [(m/s)units]

    real( kind = core_rknd ), dimension(d_variables,d_variables), intent(in) :: &
      corr_array_cloud, & ! Prescribed correlation array in cloud      [-]
      corr_array_below    ! Prescribed correlation array below cloud   [-]

    type(pdf_parameter), dimension(nz), intent(in) :: &
      pdf_params    ! PDF parameters                               [units vary]

    logical, intent(in) :: &
      l_stats_samp    ! Flag to sample statistics

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(d_variables,d_variables,nz), &
    intent(inout) :: &
      corr_array_1, & ! Correlation array for the 1st PDF component   [-]
      corr_array_2    ! Correlation array for the 2nd PDF component   [-]

    ! Output Variables
    real( kind = core_rknd ), dimension(d_variables, nz), intent(out) :: &
      mu_x_1,    & ! Mean array for the 1st PDF component                 [units vary]
      mu_x_2,    & ! Mean array for the 2nd PDF component                 [units vary]
      sigma_x_1, & ! Standard deviation array for the 1st PDF component   [units vary]
      sigma_x_2    ! Standard deviation array for the 2nd PDF component   [units vary]

    type(hydromet_pdf_parameter), dimension(nz), intent(out) :: &
      hydromet_pdf_params    ! Hydrometeor PDF parameters        [units vary]

    ! Local Variables
    real( kind = core_rknd ), dimension(d_variables,d_variables,nz) :: &
      corr_matrix_approx     ! correlation matrix                     [-]

    real( kind = core_rknd ), dimension(d_variables,d_variables) :: &
      corr_cholesky_mtx_t,     & ! transposed correlation cholesky matrix [-]
      corr_matrix_approx_swap, & ! Swapped correlation matrix             [-]
      corr_array_swap            ! Swapped correlation matrix             [-]

    real( kind = core_rknd ), dimension(nz) ::  &
      wprrp, & ! Covariance of w and r_r, < w'r_r' > (m-levs.)  [(m/s)(kg/kg)]
      wpNrp    ! Covariance of w and N_r, < w'N_r' > (m-levs.)  [(m/s)(num/kg)]

    real( kind = core_rknd ), dimension(nz) :: &
      rc1,         & ! Mean of r_c (1st PDF component)              [kg/kg]
      rc2,         & ! Mean of r_c (2nd PDF component)              [kg/kg]
      cloud_frac1, & ! Cloud fraction (1st PDF component)           [-]
      cloud_frac2, & ! Cloud fraction (2nd PDF component)           [-]
      mixt_frac      ! Mixture fraction                             [-]

    real( kind = core_rknd ) :: &
      mu_w_1,        & ! Mean of w (1st PDF component)                     [m/s]
      mu_w_2,        & ! Mean of w (2nd PDF component)                     [m/s]
      mu_s_1,        & ! Mean of s (1st PDF component)                   [kg/kg]
      mu_s_2,        & ! Mean of s (2nd PDF component)                   [kg/kg]
      mu_t_1,        & ! Mean of t (1st PDF component)                   [kg/kg]
      mu_t_2,        & ! Mean of t (2nd PDF component)                   [kg/kg]
      mu_rr_1,       & ! Mean of rr (1st PDF component) in-precip (ip)   [kg/kg]
      mu_rr_2,       & ! Mean of rr (2nd PDF component) ip               [kg/kg]
      mu_Nr_1,       & ! Mean of Nr (1st PDF component) ip              [num/kg]
      mu_Nr_2,       & ! Mean of Nr (2nd PDF component) ip              [num/kg]
      mu_Ncn_1,      & ! Mean of Ncn (1st PDF component)                [num/kg]
      mu_Ncn_2,      & ! Mean of Ncn (2nd PDF component)                [num/kg]
      mu_rr_1_n,     & ! Mean of ln rr (1st PDF component) ip        [ln(kg/kg)]
      mu_rr_2_n,     & ! Mean of ln rr (2nd PDF component) ip        [ln(kg/kg)]
      mu_Nr_1_n,     & ! Mean of ln Nr (1st PDF component) ip       [ln(num/kg)]
      mu_Nr_2_n,     & ! Mean of ln Nr (2nd PDF component) ip       [ln(num/kg)]
      mu_Ncn_1_n,    & ! Mean of ln Ncn (1st PDF component)         [ln(num/kg)]
      mu_Ncn_2_n,    & ! Mean of ln Ncn (2nd PDF component)         [ln(num/kg)]
      sigma_w_1,     & ! Standard deviation of w (1st PDF component)       [m/s]
      sigma_w_2,     & ! Standard deviation of w (2nd PDF component)       [m/s]
      sigma_s_1,     & ! Standard deviation of s (1st PDF component)     [kg/kg]
      sigma_s_2,     & ! Standard deviation of s (2nd PDF component)     [kg/kg]
      sigma_t_1,     & ! Standard deviation of t (1st PDF component)     [kg/kg]
      sigma_t_2,     & ! Standard deviation of t (2nd PDF component)     [kg/kg]
      sigma_rr_1,    & ! Standard deviation of rr (1st PDF component) ip [kg/kg]
      sigma_rr_2,    & ! Standard deviation of rr (2nd PDF component) ip [kg/kg]
      sigma_Nr_1,    & ! Standard deviation of Nr (1st PDF comp.) ip    [num/kg]
      sigma_Nr_2,    & ! Standard deviation of Nr (2nd PDF comp.) ip    [num/kg]
      sigma_Ncn_1,   & ! Standard deviation of Ncn (1st PDF component)  [num/kg]
      sigma_Ncn_2,   & ! Standard deviation of Ncn (2nd PDF component)  [num/kg]
      sigma_rr_1_n,  & ! Standard dev. of ln rr (1st PDF comp.) ip   [ln(kg/kg)]
      sigma_rr_2_n,  & ! Standard dev. of ln rr (2nd PDF comp.) ip   [ln(kg/kg)]
      sigma_Nr_1_n,  & ! Standard dev. of ln Nr (1st PDF comp.) ip  [ln(num/kg)]
      sigma_Nr_2_n,  & ! Standard dev. of ln Nr (2nd PDF comp.) ip  [ln(num/kg)]
      sigma_Ncn_1_n, & ! Standard dev. of ln Ncn (1st PDF comp.)    [ln(num/kg)]
      sigma_Ncn_2_n    ! Standard dev. of ln Ncn (2nd PDF comp.)    [ln(num/kg)]

    real( kind = core_rknd ) :: &
      corr_ws_1,     & ! Correlation between w and s (1st PDF component)     [-]
      corr_ws_2,     & ! Correlation between w and s (2nd PDF component)     [-]
      corr_wrr_1,    & ! Correlation between w and rr (1st PDF component) ip [-]
      corr_wrr_2,    & ! Correlation between w and rr (2nd PDF component) ip [-]
      corr_wNr_1,    & ! Correlation between w and Nr (1st PDF component) ip [-]
      corr_wNr_2,    & ! Correlation between w and Nr (2nd PDF component) ip [-]
      corr_wNcn_1,   & ! Correlation between w and Ncn (1st PDF component)   [-]
      corr_wNcn_2,   & ! Correlation between w and Ncn (2nd PDF component)   [-]
      corr_st_1,     & ! Correlation between s and t (1st PDF component)     [-]
      corr_st_2,     & ! Correlation between s and t (2nd PDF component)     [-]
      corr_srr_1,    & ! Correlation between s and rr (1st PDF component) ip [-]
      corr_srr_2,    & ! Correlation between s and rr (2nd PDF component) ip [-]
      corr_sNr_1,    & ! Correlation between s and Nr (1st PDF component) ip [-]
      corr_sNr_2,    & ! Correlation between s and Nr (2nd PDF component) ip [-]
      corr_sNcn_1,   & ! Correlation between s and Ncn (1st PDF component)   [-]
      corr_sNcn_2,   & ! Correlation between s and Ncn (2nd PDF component)   [-]
      corr_trr_1,    & ! Correlation between t and rr (1st PDF component) ip [-]
      corr_trr_2,    & ! Correlation between t and rr (2nd PDF component) ip [-]
      corr_tNr_1,    & ! Correlation between t and Nr (1st PDF component) ip [-]
      corr_tNr_2,    & ! Correlation between t and Nr (2nd PDF component) ip [-]
      corr_tNcn_1,   & ! Correlation between t and Ncn (1st PDF component)   [-]
      corr_tNcn_2,   & ! Correlation between t and Ncn (2nd PDF component)   [-]
      corr_rrNr_1,   & ! Correlation between rr & Nr (1st PDF component) ip  [-]
      corr_rrNr_2      ! Correlation between rr & Nr (2nd PDF component) ip  [-]

    real( kind = core_rknd ) :: &
      corr_wrr_1_n,  & ! Correlation between w and ln rr (1st PDF comp.) ip  [-]
      corr_wrr_2_n,  & ! Correlation between w and ln rr (2nd PDF comp.) ip  [-]
      corr_wNr_1_n,  & ! Correlation between w and ln Nr (1st PDF comp.) ip  [-]
      corr_wNr_2_n,  & ! Correlation between w and ln Nr (2nd PDF comp.) ip  [-]
      corr_wNcn_1_n, & ! Correlation between w and ln Ncn (1st PDF comp.)    [-]
      corr_wNcn_2_n, & ! Correlation between w and ln Ncn (2nd PDF comp.)    [-]
      corr_srr_1_n,  & ! Correlation between s and ln rr (1st PDF comp.) ip  [-]
      corr_srr_2_n,  & ! Correlation between s and ln rr (2nd PDF comp.) ip  [-]
      corr_sNr_1_n,  & ! Correlation between s and ln Nr (1st PDF comp.) ip  [-]
      corr_sNr_2_n,  & ! Correlation between s and ln Nr (2nd PDF comp.) ip  [-]
      corr_sNcn_1_n, & ! Correlation between s and ln Ncn (1st PDF comp.)    [-]
      corr_sNcn_2_n, & ! Correlation between s and ln Ncn (2nd PDF comp.)    [-]
      corr_trr_1_n,  & ! Correlation between t and ln rr (1st PDF comp.) ip  [-]
      corr_trr_2_n,  & ! Correlation between t and ln rr (2nd PDF comp.) ip  [-]
      corr_tNr_1_n,  & ! Correlation between t and ln Nr (1st PDF comp.) ip  [-]
      corr_tNr_2_n,  & ! Correlation between t and ln Nr (2nd PDF comp.) ip  [-]
      corr_tNcn_1_n, & ! Correlation between t and ln Ncn (1st PDF comp.)    [-]
      corr_tNcn_2_n, & ! Correlation between t and ln Ncn (2nd PDF comp.)    [-]
      corr_rrNr_1_n, & ! Correlation btwn. ln rr & ln Nr (1st PDF comp.) ip  [-]
      corr_rrNr_2_n    ! Correlation btwn. ln rr & ln Nr (2nd PDF comp.) ip  [-]

    real( kind = core_rknd ), dimension(nz) ::  &
      wpsp_zm,     & ! Covariance of s and w (momentum levels)   [(m/s)(kg/kg)]
      wpNcnp_zm,   & ! Covariance of N_cn and w (momentum levs.) [(m/s)(num/kg)]
      wpsp_zt,     & ! Covariance of s and w on t-levs           [(m/s)(kg/kg)]
      wprrp_ip_zt, & ! Covar. of r_r and w (in-precip) on t-levs [(m/s)(kg/kg)]
      wpNrp_ip_zt, & ! Covar. of N_r and w (in-precip) on t-levs [(m/s)(num/kg)]
      wpNcnp_zt      ! Covariance of N_cn and w on t-levs        [(m/s)(num/kg)]

     real( kind = core_rknd ), dimension(nz) :: &
      rr1, & ! Mean rain water mixing ratio (1st PDF component)      [kg/kg]
      rr2, & ! Mean rain water mixing ratio (2nd PDF component)      [kg/kg]
      Nr1, & ! Mean rain drop concentration (1st PDF component)      [num/kg]
      Nr2    ! Mean rain drop concentration (2nd PDF component)      [num/kg]

    real( kind = core_rknd ), dimension(nz) :: &
      precip_frac,   & ! Precipitation fraction (overall)           [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component) [-]
      precip_frac_2    ! Precipitation fraction (2nd PDF component) [-]

    integer :: k  ! Loop index

    ! ---- Begin Code ----

    ! Covariance of vertical velocity and a hydrometeor
    ! (< w'r_r' > and < w'N_r' >).
    wprrp = wphydrometp(:,iirrainm)
    wpNrp = wphydrometp(:,iiNrm)

    ! Setup some of the PDF parameters
    rc1         = pdf_params%rc1
    rc2         = pdf_params%rc2
    cloud_frac1 = pdf_params%cloud_frac1
    cloud_frac2 = pdf_params%cloud_frac2
    mixt_frac   = pdf_params%mixt_frac

    ! Component mean values for r_r and N_r, and precipitation fraction.
    if ( l_use_precip_frac ) then

       call component_means_rain( nz, rrainm, Nrm, rho, rc1, rc2, & ! Intent(in)
                                  mixt_frac, l_stats_samp, &        ! Intent(in)
                                  rr1, rr2, Nr1, Nr2 )              ! Intent(out)

       call precip_fraction( nz, rrainm, rr1, rr2, Nrm, Nr1, Nr2, &     ! Intent(in)
                             cloud_frac, cloud_frac1, mixt_frac, &      ! Intent(in)
                             precip_frac, precip_frac_1, precip_frac_2 )! Intent(out)

    else

       rr1 = rrainm
       rr2 = rrainm
       Nr1 = Nrm
       Nr2 = Nrm

       precip_frac   = one
       precip_frac_1 = one
       precip_frac_2 = one

    endif

    ! Statistics
    if ( l_stats_samp ) then

       if ( irr1 > 0 ) then
          ! Mean rain water mixing ratio in PDF component 1.
          call stat_update_var( irr1, rr1, zt )
       endif

       if ( irr2 > 0 ) then
          ! Mean rain water mixing ratio in PDF component 2.
          call stat_update_var( irr2, rr2, zt )
       endif

       if ( iNr1 > 0 ) then
          ! Mean rain drop concentration in PDF component 1.
          call stat_update_var( iNr1, Nr1, zt )
       endif

       if ( iNr2 > 0 ) then
          ! Mean rain drop concentration in PDF component 2.
          call stat_update_var( iNr2, Nr2, zt )
       endif

       if ( iprecip_frac > 0 ) then
          ! Overall precipitation fraction.
          call stat_update_var( iprecip_frac, precip_frac, zt )
       endif

       if ( iprecip_frac_1 > 0 ) then
          ! Precipitation fraction in PDF component 1.
          call stat_update_var( iprecip_frac_1, precip_frac_1, zt )
       endif

       if ( iprecip_frac_2 > 0 ) then
          ! Precipitation fraction in PDF component 2.
          call stat_update_var( iprecip_frac_2, precip_frac_2, zt )
       endif

    endif

    ! Calculate correlations involving w by first calculating total covariances
    ! involving w (<w'r_r'>, etc.) using the down-gradient approximation.
    if ( l_calc_w_corr ) then

       ! Calculate the covariances of w with the hydrometeors
       do k = 1, nz
          wpsp_zm(k) = pdf_params(k)%mixt_frac &
                       * ( one - pdf_params(k)%mixt_frac ) &
                       * ( pdf_params(k)%s1 - pdf_params(k)%s2 ) &
                       * ( pdf_params(k)%w1 - pdf_params(k)%w2 )
       enddo

       wpNcnp_zm(1:nz-1) = xpwp_fnc( -c_Krrainm * Kh_zm(1:nz-1), Ncnm(1:nz-1), &
                                     Ncnm(2:nz), gr%invrs_dzm(1:nz-1) )

       ! Boundary conditions; We are assuming zero flux at the top.
       wpNcnp_zm(nz) = zero

       ! Interpolate the covariances to thermodynamic grid levels.
       wpsp_zt     = zm2zt(wpsp_zm)
       wprrp_ip_zt = zm2zt(wprrp) / max( precip_frac, eps )
       wpNrp_ip_zt = zm2zt(wpNrp) / max( precip_frac, eps )
       wpNcnp_zt   = zm2zt(wpNcnp_zm)

       ! When the mean value of a hydrometeor is below tolerance value, it is
       ! considered to have a value of 0, and does not vary over the grid level.
       ! Any covariance involving that hydrometeor also has a value of 0 at that
       ! grid level.
       do k = 1, nz, 1
          if ( rrainm(k) <= rr_tol ) then
             wprrp_ip_zt(k) = zero
          endif
          if ( Nrm(k) <= Nr_tol ) then
             wpNrp_ip_zt(k) = zero
          endif
          if ( Ncnm(k) <= Ncn_tol ) then
             wpNcnp_zt(k) = zero
          endif
       enddo

    endif


    !!! Setup PDF parameters loop.
    ! Loop over all model thermodynamic level above the model lower boundary.
    do k = 2, nz, 1

       !!! Calculate the means, standard deviations, and necessary correlations
       !!! involving w, s, t, r_r (in-precip), N_r (in-precip), and N_cn for
       !!! each PDF component.
       call compute_mean_stdev( rcm(k), rrainm(k), Nrm(k), Ncnm(k), &         ! Intent(in)
                                rr1(k), rr2(k), Nr1(k), Nr2(k), rc1(k), &     ! Intent(in)
                                rc2(k), cloud_frac1(k), cloud_frac2(k), &     ! Intent(in)
                                precip_frac_1(k), precip_frac_2(k), &         ! Intent(in)
                                wpsp_zt(k), wprrp_ip_zt(k), wpNrp_ip_zt(k), & ! Intent(in)
                                wpNcnp_zt(k), w_std_dev(k), mixt_frac(k), &   ! Intent(in)
                                pdf_params(k), &                              ! Intent(in)
                                mu_w_1, mu_w_2, mu_s_1, mu_s_2, mu_t_1, &     ! Intent(out)
                                mu_t_2, mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, & ! Intent(out)
                                mu_Ncn_1, mu_Ncn_2, sigma_w_1, sigma_w_2, &   ! Intent(out)
                                sigma_s_1, sigma_s_2, sigma_t_1, sigma_t_2, & ! Intent(out)
                                sigma_rr_1, sigma_rr_2, sigma_Nr_1, &         ! Intent(out)
                                sigma_Nr_2, sigma_Ncn_1, sigma_Ncn_2 )        ! Intent(out)

       if ( l_diagnose_correlations ) then

!          if ( rrainm > rr_tol ) then
!             rrp2_on_rrm2 = (sigma_rr_1/mu_rr_1)**2
!          else
!             ! The ratio is undefined; set it equal to 0.
!             rrp2_on_rrm2 = zero
!          endif
!
!          if ( Nrm > Nr_tol ) then
!             Nrp2_on_Nrm2 = (sigma_Nr_1/mu_Nr_1)**2
!          else
!             ! The ratio is undefined; set it equal to 0.
!             Nrp2_on_Nrm2 = zero
!          endif
!
!          if ( Ncnm > Ncn_tol ) then
!             Ncnp2_on_Ncnm2 = (sigma_Ncn_1/mu_Ncn_1)**2
!          else
!             ! The ratio is undefined; set it equal to 0.
!             Ncnp2_on_Ncnm2 = zero
!          endif

          if ( rcm(k) > rc_tol ) then

             call diagnose_correlations( d_variables, corr_array_cloud, & ! Intent(in)
                                         corr_ws_1, corr_wrr_1, &         ! Intent(out)
                                         corr_wNr_1, corr_wNcn_1, &       ! Intent(out)
                                         corr_st_1, corr_srr_1, &         ! Intent(out)
                                         corr_sNr_1, corr_sNcn_1, &       ! Intent(out)
                                         corr_trr_1, corr_tNr_1, &        ! Intent(out)
                                         corr_tNcn_1, corr_rrNr_1 )       ! Intent(out)

             call diagnose_correlations( d_variables, corr_array_cloud, & ! Intent(in)
                                         corr_ws_2, corr_wrr_2, &         ! Intent(out)
                                         corr_wNr_2, corr_wNcn_2, &       ! Intent(out)
                                         corr_st_2, corr_srr_2, &         ! Intent(out)
                                         corr_sNr_2, corr_sNcn_2, &       ! Intent(out)
                                         corr_trr_2, corr_tNr_2, &        ! Intent(out)
                                         corr_tNcn_2, corr_rrNr_2 )       ! Intent(out)

          else

             call diagnose_correlations( d_variables, corr_array_below, & ! Intent(in)
                                         corr_ws_1, corr_wrr_1, &         ! Intent(out)
                                         corr_wNr_1, corr_wNcn_1, &       ! Intent(out)
                                         corr_st_1, corr_srr_1, &         ! Intent(out)
                                         corr_sNr_1, corr_sNcn_1, &       ! Intent(out)
                                         corr_trr_1, corr_tNr_1, &        ! Intent(out)
                                         corr_tNcn_1, corr_rrNr_1 )       ! Intent(out)

             call diagnose_correlations( d_variables, corr_array_below, & ! Intent(in)
                                         corr_ws_2, corr_wrr_2, &         ! Intent(out)
                                         corr_wNr_2, corr_wNcn_2, &       ! Intent(out)
                                         corr_st_2, corr_srr_2, &         ! Intent(out)
                                         corr_sNr_2, corr_sNcn_2, &       ! Intent(out)
                                         corr_trr_2, corr_tNr_2, &        ! Intent(out)
                                         corr_tNcn_2, corr_rrNr_2 )       ! Intent(out)

          endif

       else ! if .not. l_diagnose_correlations

          call compute_corr( rcm(k), rrainm(k), Nrm(k), Ncnm(k), &         ! Intent(in)
                             rr1(k), rr2(k), Nr1(k), Nr2(k), rc1(k), &     ! Intent(in)
                             rc2(k), cloud_frac1(k), cloud_frac2(k), &     ! Intent(in)
                             precip_frac_1(k), precip_frac_2(k), &         ! Intent(in)
                             wpsp_zt(k), wprrp_ip_zt(k), wpNrp_ip_zt(k), & ! Intent(in)
                             wpNcnp_zt(k), w_std_dev(k), mixt_frac(k), &   ! Intent(in)
                             sigma_rr_1, sigma_Nr_1, sigma_Ncn_1, &        ! Intent(in)
                             pdf_params(k), &                              ! Intent(in)
                             corr_ws_1, corr_ws_2, corr_wrr_1, &           ! Intent(out)
                             corr_wrr_2, corr_wNr_1, corr_wNr_2, &         ! Intent(out)
                             corr_wNcn_1, corr_wNcn_2, corr_st_1, &        ! Intent(out)
                             corr_st_2, corr_srr_1, corr_srr_2, &          ! Intent(out)
                             corr_sNr_1, corr_sNr_2, corr_sNcn_1, &        ! Intent(out)
                             corr_sNcn_2, corr_trr_1, corr_trr_2, &        ! Intent(out)
                             corr_tNr_1, corr_tNr_2, corr_tNcn_1, &        ! Intent(out)
                             corr_tNcn_2, corr_rrNr_1, corr_rrNr_2 )       ! Intent(out)

       endif ! l_diagnose_correlations

       !!! Calculate the mean, standard deviations, and correlations involving
       !!! ln r_r, ln N_r, and ln N_cn for each PDF component.
       call normalize_pdf_params( rr1(k), rr2(k), Nr1(k), Nr2(k), Ncnm(k), &     ! Intent(in)
                                   mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, &         ! Intent(in)
                                   mu_Ncn_1, mu_Ncn_2, sigma_rr_1, sigma_rr_2, & ! Intent(in)
                                   sigma_Nr_1, sigma_Nr_2, sigma_Ncn_1, &        ! Intent(in)
                                   sigma_Ncn_2, corr_wrr_1, corr_wrr_2, &        ! Intent(in)
                                   corr_wNr_1, corr_wNr_2, corr_wNcn_1, &        ! Intent(in)
                                   corr_wNcn_2, corr_srr_1, corr_srr_2, &        ! Intent(in)
                                   corr_sNr_1, corr_sNr_2, corr_sNcn_1, &        ! Intent(in)
                                   corr_sNcn_2, corr_trr_1, corr_trr_2, &        ! Intent(in)
                                   corr_tNr_1, corr_tNr_2, corr_tNcn_1, &        ! Intent(in)
                                   corr_tNcn_2, corr_rrNr_1, corr_rrNr_2, &      ! Intent(in)
                                   mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, &            ! Intent(out)
                                   mu_Nr_2_n, mu_Ncn_1_n, mu_Ncn_2_n, &          ! Intent(out)
                                   sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, &   ! Intent(out)
                                   sigma_Nr_2_n, sigma_Ncn_1_n, sigma_Ncn_2_n, & ! Intent(out)
                                   corr_wrr_1_n, corr_wrr_2_n, corr_wNr_1_n, &   ! Intent(out)
                                   corr_wNr_2_n, corr_wNcn_1_n, corr_wNcn_2_n, & ! Intent(out)
                                   corr_srr_1_n, corr_srr_2_n, corr_sNr_1_n, &   ! Intent(out)
                                   corr_sNr_2_n, corr_sNcn_1_n, corr_sNcn_2_n, & ! Intent(out)
                                   corr_trr_1_n, corr_trr_2_n, corr_tNr_1_n, &   ! Intent(out)
                                   corr_tNr_2_n, corr_tNcn_1_n, corr_tNcn_2_n, & ! Intent(out)
                                   corr_rrNr_1_n, corr_rrNr_2_n )                ! Intent(out)

       !!! Statistics
       call pdf_param_hm_stats(  mu_rr_1, mu_rr_2, mu_Nr_1, &
                                 mu_Nr_2, mu_Ncn_1, mu_Ncn_2, &
                                 sigma_rr_1, sigma_rr_2, sigma_Nr_1, &
                                 sigma_Nr_2, sigma_Ncn_1, sigma_Ncn_2, &
                                 corr_wrr_1, corr_wrr_2, corr_wNr_1, &
                                 corr_wNr_2, corr_wNcn_1, corr_wNcn_2, &
                                 corr_srr_1, corr_srr_2, corr_sNr_1, &
                                 corr_sNr_2, corr_sNcn_1, corr_sNcn_2, &
                                 corr_trr_1, corr_trr_2, corr_tNr_1, &
                                 corr_tNr_2, corr_tNcn_1, corr_tNcn_2, &
                                 corr_rrNr_1, corr_rrNr_2, &
                                 k, l_stats_samp )

       call pdf_param_log_hm_stats( mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, &
                                    mu_Nr_2_n, mu_Ncn_1_n, mu_Ncn_2_n, &
                                    sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, &
                                    sigma_Nr_2_n, sigma_Ncn_1_n, sigma_Ncn_2_n, &
                                    corr_wrr_1_n, corr_wrr_2_n, corr_wNr_1_n, &
                                    corr_wNr_2_n, corr_wNcn_1_n, corr_wNcn_2_n, &
                                    corr_srr_1_n, corr_srr_2_n, corr_sNr_1_n, &
                                    corr_sNr_2_n, corr_sNcn_1_n, corr_sNcn_2_n, &
                                    corr_trr_1_n, corr_trr_2_n, corr_tNr_1_n, &
                                    corr_tNr_2_n, corr_tNcn_1_n, corr_tNcn_2_n, &
                                    corr_rrNr_1_n, corr_rrNr_2_n, k, &
                                    l_stats_samp )

       !!! Pack the PDF parameters
       call pack_pdf_params( mu_w_1, mu_w_2, mu_s_1, mu_s_2, &                     ! Intent(in)
                             mu_t_1, mu_t_2, mu_rr_1, mu_rr_2, &                   ! Intent(in)
                             mu_Nr_1, mu_Nr_2, mu_Ncn_1, mu_Ncn_2, &               ! Intent(in)
                             mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, &                    ! Intent(in)
                             mu_Nr_2_n, mu_Ncn_1_n, mu_Ncn_2_n, &                  ! Intent(in)
                             sigma_w_1, sigma_w_2, sigma_s_1, &                    ! Intent(in)
                             sigma_s_2, sigma_t_1, sigma_t_2, &                    ! Intent(in)
                             sigma_rr_1, sigma_rr_2, sigma_Nr_1, &                 ! Intent(in)
                             sigma_Nr_2, sigma_Ncn_1, sigma_Ncn_2, &               ! Intent(in)
                             sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, &           ! Intent(in)
                             sigma_Nr_2_n, sigma_Ncn_1_n, sigma_Ncn_2_n, &         ! Intent(in)
                             corr_ws_1, corr_ws_2, corr_st_1, corr_st_2, &         ! Intent(in)
                             corr_wrr_1_n, corr_wrr_2_n, corr_wNr_1_n, &           ! Intent(in)
                             corr_wNr_2_n, corr_wNcn_1_n, corr_wNcn_2_n, &         ! Intent(in)
                             corr_srr_1_n, corr_srr_2_n, corr_sNr_1_n, &           ! Intent(in)
                             corr_sNr_2_n, corr_sNcn_1_n, corr_sNcn_2_n, &         ! Intent(in)
                             corr_trr_1_n, corr_trr_2_n, corr_tNr_1_n, &           ! Intent(in)
                             corr_tNr_2_n, corr_tNcn_1_n, corr_tNcn_2_n, &         ! Intent(in)
                             corr_rrNr_1_n, corr_rrNr_2_n, &                       ! Intent(in)
                             rr1(k), rr2(k), Nr1(k), Nr2(k), &                     ! Intent(in)
                             precip_frac(k), precip_frac_1(k), precip_frac_2(k), & ! Intent(in)
                             d_variables, &                                        ! Intent(in)
                             corr_array_1(:,:,k), corr_array_2(:,:,k), &           ! Intent(inout)
                             mu_x_1(:,k), mu_x_2(:,k), &                           ! Intent(out)
                             sigma_x_1(:,k), sigma_x_2(:,k), &                     ! Intent(out)
                             hydromet_pdf_params(k) )                              ! Intent(out)

       if ( l_use_modified_corr ) then

          call rearrange_corr_array( d_variables, corr_array_1(:,:,k), & ! Intent(in)
                                     corr_array_swap)                    ! Intent(inout)

          call setup_corr_cholesky_mtx( d_variables, corr_array_swap, & ! intent(in)
                                        corr_cholesky_mtx_t )           ! intent(out)

          call cholesky_to_corr_mtx_approx( d_variables, corr_cholesky_mtx_t, & ! intent(in)
                                            corr_matrix_approx_swap )           ! intent(out)

          call rearrange_corr_array( d_variables, corr_matrix_approx_swap, & ! Intent(in)
                                     corr_matrix_approx(:,:,k))              ! Intent(inout)

          call corr_array_assertion_checks( d_variables, corr_matrix_approx(:, :, k) )

       endif

    enddo  ! Setup PDF parameters loop: k = 2, nz, 1

    return

  end subroutine setup_pdf_parameters

  !=============================================================================
  subroutine component_means_rain( nz, rrainm, Nrm, rho, rc1, rc2, &
                                   mixt_frac, l_stats_samp, &
                                   rr1, rr2, Nr1, Nr2 )

    ! Description:
    ! The values of grid-level mean rain water mixing ratio, <r_r>, and
    ! grid-level mean rain drop concentration, <N_r>, are solved as part of the
    ! CLUBB model's predictive equation set.  However, CLUBB has a two component
    ! PDF.  The grid-level means of r_r and N_r must be subdivided into
    ! component means for each PDF component.  The equations relating the
    ! overall means to the component means are:
    !
    ! <r_r> = a * rr1 + (1-a) * rr2, and
    ! <N_r> = a * Nr1 + (1-a) * Nr2;
    !
    ! where "a" is the mixture fraction (weight of the 1st PDF component), rr1
    ! is the mean rain water mixing ratio in PDF component 1, rr2 is the mean
    ! rain water mixing ratio in PDF component 2, Nr1 is the mean rain drop
    ! concentration in PDF component 1, and Nr2 is the mean rain drop
    ! concentration in PDF component 2.  These equations can be rewritten as:
    !
    ! <r_r> = rr1 * ( a + (1-a) * rr2/rr1 ), and
    ! <N_r> = Nr1 * ( a + (1-a) * Nr2/Nr1 ).
    !
    ! One way to solve for a component mean is to relate the ratios rr2/rr1 and
    ! Nr2/Nr1 to other factors.  For now, these ratios based on other factors
    ! will be called rr2_rr1_ratio (for rr2/rr1) and Nr2_Nr1_ratio
    ! (for Nr2/Nr1).  These ratios are entered into the above equations,
    ! allowing the equations to be solved for rr1 and Nr1:
    !
    ! rr1 = <r_r> / ( a + (1-a) * rr2_rr1_ratio ), and
    ! Nr1 = <N_r> / ( a + (1-a) * Nr2_Nr1_ratio ).
    !
    ! Once that rr1 and Nr1 have been solved, rr2 and Nr2 can be solved by:
    !
    ! rr2 = ( <r_r> - a * rr1 ) / (1-a); and
    ! Nr2 = ( <N_r> - a * Nr1 ) / (1-a).
    !
    ! At a grid level that is at least mostly cloudy, the simplest way to handle
    ! the ratios rr2/rr1 and Nr2/Nr1 is to set them equal to the ratio rc2/rc1,
    ! where rc1 is the mean cloud water mixing ratio in PDF component 1 and rc2
    ! is the mean cloud water mixing ratio in PDF component 2.  However, rain
    ! sediments, falling from higher altitudes downwards.  The values of cloud
    ! water mixing ratio at a given grid level are not necessarily indicative
    ! of the amount of cloud water at higher levels, which has already produced
    ! rain which has fallen downwards to the given grid level.  Additionally,
    ! using grid-level cloud water mixing ratio especially does not work for
    ! rain below cloud base (near the ground).
    !
    ! However, an alternative to component cloud water mixing ratio is component
    ! liquid water path.  Liquid water path accounts for the cloud water mixing
    ! ratio at the given grid level and at all grid levels higher in altitude.
    !
    ! In a stratocumulus case, the cloud water is spread out over all or almost
    ! all of the horizontal domain over a group of vertical levels.  At a given
    ! vertical level, the component mean cloud water mixing ratios should be
    ! almost equal, although usually slightly larger in the component with the
    ! larger component mean extended liquid water mixing ratio, s.  Likewise,
    ! the component liquid water paths should be nearly equal, with one
    ! component having a slightly larger liquid water path than the other
    ! component.
    ! 
    ! In a case of cumulus rising into stratocumulus, the upper portion of the
    ! cloudy domain will be very similar to the stratocumulus case described
    ! above, with similar cloud water mixing ratio and liquid water path
    ! results.  However, below the base of the stratocumulus clouds, where the
    ! cumulus clouds are found, the horizontal domain at each vertical level is
    ! only partially cloudy.  At these levels, rain produced in the
    ! stratocumulus clouds above is evaporating in the clear-air portions, while
    ! rain is not evaporating in the cloudy portions.  Additionally, more rain
    ! is being produced in the cloudy portions.  The rain in the cloudy portions
    ! becomes significantly larger than the rain in the clear portions.  The
    ! partiallly cloudy levels usually have a PDF where one component is
    ! significantly more saturated than the other component.  By the time the
    ! cloud base of the cumulus clouds is reached, the liquid water path for one
    ! PDF component should be significantly greater than the liquid water path
    ! for the other PDF component.
    ! 
    ! In a cumulus case, the horizontal domain at each level is usually partly
    ! cloudy.  Throughout the entire vertical domain, at every vertical level,
    ! one component usually is much more saturated than the other component.
    ! The liquid water path for one component is much greater than the liquid
    ! water path in the other component.  Likewise, rain that is formed in cloud
    ! and falls preferentially through cloud will have large values in a portion
    ! of the horizontal domain and very small or 0 values over the rest of the 
    ! horizontal domain.
    !
    ! In order to estimate the amount of rain in each PDF component, the ratios
    ! rr2/rr1 and Nr2/Nr1 are going to be set equal to the ratio LWP2/LWP1,
    ! where LWP1 is the liquid water path in PDF component 1 and LWP2 is the
    ! liquid water path in PDF component 2.  LWP1 will be computed by taking the
    ! vertical integral of cloud water (see equation below) through the 1st PDF
    ! component from the given vertical level all the way to the top of the
    ! model.  LWP2 will be computed in the same manner.   It should be noted
    ! that this method makes the poor assumption that PDF component 1 always
    ! overlaps PDF component 1 between vertical levels, and likewise for PDF
    ! component 2.
    !
    ! Total liquid water path, LWP, is given by the following equation:
    !
    ! LWP(z) = INT(z:z_top) rho_a <r_c> dz';
    !
    ! where z is the altitude of the vertical level for which LWP is desired,
    ! z_top is the altitude at the top of the model domain, and z' is the
    ! dummy variable of integration.  Mean cloud water mixing ratio can be
    ! written as:
    !
    ! <r_c> = a * rc1 + (1-a) * rc2.
    !
    ! The equation for liquid water path is rewritten as:
    !
    ! LWP(z) = INT(z:z_top) rho_a ( a rc1 + (1-a) rc2 ) dz'; or
    !
    ! LWP(z) = INT(z:z_top) a rho_a rc1 dz'
    !          + INT(z:z_top) (1-a) rho_a rc2 dz'.
    !
    ! This can be rewritten as:
    !
    ! LWP(z) = LWP1(z) + LWP2(z);
    !
    ! where:
    !
    ! LWP1(z) = INT(z:z_top) a rho_a rc1 dz'; and
    ! LWP2(z) = INT(z:z_top) (1-a) rho_a rc2 dz'.
    !
    ! The trapezoidal rule will be used to numerically integrate for LWP1
    ! and LWP2.
    

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        gr    ! Variable(s)

    use constants_clubb, only: &
        one,      & ! Constant(s)
        one_half, &
        zero,     &
        rr_tol,   &
        Nr_tol

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    use stats_type, only: &
        stat_update_var  ! Procedure(s)

    use stats_variables, only : &
        iLWP1, & ! Variable(s)
        iLWP2, &
        zt

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz           ! Number of model vertical grid levels

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      rrainm,    & ! Overall mean rain water mixing ratio               [kg/kg]
      Nrm,       & ! Overall mean rain drop concentration               [num/kg]
      rho,       & ! Air density                                        [kg/m^3]
      rc1,       & ! Mean cloud water mixing ratio (1st PDF component)  [kg/kg]
      rc2,       & ! Mean cloud water mixing ratio (2nd PDF component)  [kg/kg]
      mixt_frac    ! Mixture fraction                                   [-]

    logical, intent(in) :: &
      l_stats_samp     ! Flag to record statistical output.

    ! Output Variables
    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      rr1, & ! Mean rain water mixing ratio (1st PDF component)      [kg/kg]
      rr2, & ! Mean rain water mixing ratio (2nd PDF component)      [kg/kg]
      Nr1, & ! Mean rain drop concentration (1st PDF component)      [num/kg]
      Nr2    ! Mean rain drop concentration (2nd PDF component)      [num/kg]

    ! Local Variable
    real( kind = core_rknd ), dimension(nz) :: &
      LWP1, & ! Liquid water path (1st PDF component) on thermo. levs.  [kg/m^2]
      LWP2    ! Liquid water path (2nd PDF component) on thermo. levs.  [kg/m^2]

    integer :: k  ! Array index

    real( kind = core_rknd ), parameter :: &
      LWP_tol = 5.0e-7_core_rknd  ! Tolerance value for component LWP


    !!! Compute component liquid water paths using trapezoidal rule for
    !!! numerical integration.

    ! At the uppermost thermodynamic level (k = nz), use the trapezoidal rule:
    !
    ! 0.5 * (integrand_a + integrand_b) * delta_z,
    !
    ! where integrand_a is the integrand at thermodynamic level k = nz,
    ! integrand_b is the integrand at momentum level k = nz (model upper
    ! boundary), and delta_z = zm(nz) - zt(nz).  At the upper boundary, r_c is
    ! set to 0, and the form of the trapezoidal rule is simply:
    !
    ! 0.5 * integrand_a * delta_z.

    ! Liquid water path in PDF component 1.
    LWP1(nz) &
    = one_half * mixt_frac(nz) * rho(nz) * rc1(nz) * ( gr%zm(nz) - gr%zt(nz) )

    ! Liquid water path in PDF component 2.
    LWP2(nz) &
    = one_half * ( one - mixt_frac(nz) ) * rho(nz) * rc2(nz) &
      * ( gr%zm(nz) - gr%zt(nz) )

    ! At all other thermodynamic levels, compute liquid water path using the
    ! trapezoidal rule:
    !
    ! 0.5 * (integrand_a + integrand_b) * delta_z,
    !
    ! where integrand_a is the integrand at thermodynamic level k, integrand_b
    ! is the integrand at thermodynamic level k+1, and
    ! delta_z = zt(k+1) - zt(k), or 1/invrs_dzm(k).  The total for the segment
    ! is added to the sum total of all higher vertical segments to compute the
    ! total vertical integral.
    do k = nz-1, 1, -1

       ! Liquid water path in PDF component 1.
       LWP1(k) &
       = LWP1(k+1) &
         + one_half * ( mixt_frac(k+1) * rho(k+1) * rc1(k+1) &
                        + mixt_frac(k) * rho(k) * rc1(k) ) / gr%invrs_dzm(k)

       ! Liquid water path in PDF component 2.
       LWP2(k) &
       = LWP2(k+1) &
         + one_half * ( ( one - mixt_frac(k+1) ) * rho(k+1) * rc2(k+1) &
                        + ( one - mixt_frac(k) ) * rho(k) * rc2(k) ) &
           / gr%invrs_dzm(k)

    enddo ! k = nz-1, 1, -1


    !!! Find rr1, rr2, Nr1, and Nr2 based on the ratio of LWP2/LWP1, such that:
    !!! rr2/rr1 = Nr2/Nr1 = LWP2/LWP1.
    do k = 1, nz, 1

       !!! Calculate the component means for rain water mixing ratio.
       if ( rrainm(k) > rr_tol ) then

          if ( LWP1(k) <= LWP_tol .and. LWP2(k) <= LWP_tol ) then

             ! Both LWP1 and LWP2 are 0 (or an insignificant amount).
             !
             ! There is rain at this level, yet no cloud at or above the
             ! current level.  This is usually due to a numerical artifact.
             ! For example, rain is diffused above cloud top.  Simply set
             ! each component mean equal to the overall mean.
             rr1(k) = rrainm(k)
             rr2(k) = rrainm(k)

          elseif ( LWP1(k) > LWP_tol .and. LWP2(k) <= LWP_tol ) then

             ! LWP1 is (significantly) greater than 0, while LWP2 is 0 (or an
             ! insignificant amount).
             !
             ! There is rain at this level, and all cloud water at or above
             ! this level is found in the 1st PDF component.  All rain water
             ! mixing ratio is found in the 1st PDF component.
             rr1(k) = rrainm(k) / mixt_frac(k)
             rr2(k) = zero

          elseif ( LWP2(k) > LWP_tol .and. LWP1(k) <= LWP_tol ) then

             ! LWP2 is (significantly) greater than 0, while LWP1 is 0 (or an
             ! insignificant amount).
             !
             ! There is rain at this level, and all cloud water at or above
             ! this level is found in the 2nd PDF component.  All rain water
             ! mixing ratio is found in the 2nd PDF component.
             rr1(k) = zero
             rr2(k) = rrainm(k) / ( one - mixt_frac(k) )

          else ! LWP1(k) > LWP_tol and LWP2(k) > LWP_tol

             ! Both LWP1 and LWP2 are (significantly) greater than 0.
             !
             ! There is rain at this level, and there is sufficient cloud water
             ! at or above this level in both PDF components to find rain in
             ! both PDF components.  Delegate rain water mixing ratio between
             ! the 1st and 2nd PDF components according to the above equations.
             rr1(k) &
             = rrainm(k) &
               / ( mixt_frac(k) + ( one - mixt_frac(k) ) * LWP2(k)/LWP1(k) )

             rr2(k) &
             = ( rrainm(k) - mixt_frac(k) * rr1(k) ) / ( one - mixt_frac(k) )

          endif


       else ! rrainm(k) <= rr_tol

          ! Overall mean rain water mixing ratio is either 0 or below tolerance
          ! value (any postive value is considered to be a numerical artifact).
          ! Simply set each component mean equal to the overall mean.
           rr1(k) = rrainm(k)
           rr2(k) = rrainm(k)

       endif


       !!! Calculate the component means for rain drop concentration.
       if ( Nrm(k) > Nr_tol ) then

          if ( LWP1(k) <= LWP_tol .and. LWP2(k) <= LWP_tol ) then

             ! Both LWP1 and LWP2 are 0 (or an insignificant amount).
             !
             ! There is rain at this level, yet no cloud at or above the
             ! current level.  This is usually due to a numerical artifact.
             ! For example, rain is diffused above cloud top.  Simply set
             ! each component mean equal to the overall mean.
             Nr1(k) = Nrm(k)
             Nr2(k) = Nrm(k)

          elseif ( LWP1(k) > LWP_tol .and. LWP2(k) <= LWP_tol ) then

             ! LWP1 is (significantly) greater than 0, while LWP2 is 0 (or an
             ! insignificant amount).
             !
             ! There is rain at this level, and all cloud water at or above
             ! this level is found in the 1st PDF component.  All rain drop
             ! concentration is found in the 1st PDF component.
             Nr1(k) = Nrm(k) / mixt_frac(k)
             Nr2(k) = zero

          elseif ( LWP2(k) > LWP_tol .and. LWP1(k) <= LWP_tol ) then

             ! LWP2 is (significantly) greater than 0, while LWP1 is 0 (or an
             ! insignificant amount).
             !
             ! There is rain at this level, and all cloud water at or above
             ! this level is found in the 2nd PDF component.  All rain drop
             ! concentration is found in the 2nd PDF component.
             Nr1(k) = zero
             Nr2(k) = Nrm(k) / ( one - mixt_frac(k) )

          else ! LWP1(k) > LWP_tol and LWP2(k) > LWP_tol

             ! Both LWP1 and LWP2 are (significantly) greater than 0.
             !
             ! There is rain at this level, and there is sufficient cloud water
             ! at or above this level in both PDF components to find rain in
             ! both PDF components.  Delegate rain drop concentration between
             ! the 1st and 2nd PDF components according to the above equations.
             Nr1(k) &
             = Nrm(k) &
               / ( mixt_frac(k) + ( one - mixt_frac(k) ) * LWP2(k)/LWP1(k) )

             Nr2(k) &
             = ( Nrm(k) - mixt_frac(k) * Nr1(k) ) / ( one - mixt_frac(k) )

          endif


       else ! Nrm(k) <= Nr_tol

          ! Overall mean rain drop concentration is either 0 or below tolerance
          ! value (any postive value is considered to be a numerical artifact).
          ! Simply set each component mean equal to the overall mean.
           Nr1(k) = Nrm(k)
           Nr2(k) = Nrm(k)

       endif


    enddo ! k = 1, nz, 1

    ! Statistics
    if ( l_stats_samp ) then

       if ( iLWP1 > 0 ) then
          ! Liquid water path in PDF component 1.
          call stat_update_var( iLWP1, LWP1, zt )
       endif

       if ( iLWP2 > 0 ) then
          ! Liquid water path in PDF component 2.
          call stat_update_var( iLWP2, LWP2, zt )
       endif
       
    endif


    return

  end subroutine component_means_rain

  !=============================================================================
  subroutine precip_fraction( nz, rrainm, rr1, rr2, Nrm, Nr1, Nr2, &
                              cloud_frac, cloud_frac1, mixt_frac, &
                              precip_frac, precip_frac_1, precip_frac_2 )

    ! Description:
    ! Determines (overall) precipitation fraction over the horizontal domain, as
    ! well as the precipitation fraction within each PDF component, at every
    ! vertical grid level.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one,            & ! Constant(s)
        zero,           &
        rr_tol,         &
        Nr_tol,         &
        cloud_frac_min

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz          ! Number of model vertical grid levels

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      rrainm,      & ! Mean rain water mixing ratio (overall)           [kg/kg]
      rr1,         & ! Mean rain water mixing ratio (1st PDF component) [kg/kg]
      rr2,         & ! Mean rain water mixing ratio (2nd PDF component) [kg/kg]
      Nrm,         & ! Mean rain drop concentration (overall)           [num/kg]
      Nr1,         & ! Mean rain drop concentration (1st PDF component) [num/kg]
      Nr2,         & ! Mean rain drop concentration (2nd PDF component) [num/kg]
      cloud_frac,  & ! Cloud fraction (overall)                         [-] 
      cloud_frac1, & ! Cloud fraction (1st PDF component)               [-]
      mixt_frac      ! Mixture fraction                                 [-]

    ! Output Variables
    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      precip_frac,   & ! Precipitation fraction (overall)               [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component)     [-]
      precip_frac_2    ! Precipitation fraction (2nd PDF component)     [-]

    ! Local Variables
    real( kind = core_rknd ), dimension(nz) :: &
      weighted_pfrac1    ! Product of mixt_frac and precip_frac_1       [-]

    real( kind = core_rknd ), parameter :: &
      precip_frac_tol = cloud_frac_min  ! Minimum precip. frac.         [-]
    
    integer :: &
      k   ! Loop index


    !!! Find overall precipitation fraction.
    do k = nz, 1, -1

       ! The precipitation fraction is the greatest cloud fraction at or above a
       ! vertical level.
       if ( k < nz ) then
          precip_frac(k) = max( precip_frac(k+1), cloud_frac(k) )
       else  ! k = nz
          precip_frac(k) = cloud_frac(k)
       endif

       if ( ( rrainm(k) > rr_tol .or. Nrm(k) > Nr_tol ) &
            .and. precip_frac(k) < precip_frac_tol ) then

          ! In a scenario where we find rain at this grid level, but no cloud at
          ! or above this grid level, set precipitation fraction to a minimum
          ! threshold value.
          precip_frac(k) = precip_frac_tol

       elseif ( ( rrainm(k) < rr_tol .and. Nrm(k) < Nr_tol ) &
                .and. precip_frac(k) < precip_frac_tol ) then

          ! Mean (overall) rain water mixing ratio and mean (overall) rain drop
          ! concentration are both less than their respective tolerance amounts.
          ! They are both considered to have values of 0.  There is not any rain
          ! at this grid level.  There is also no cloud at or above this grid
          ! level, so set precipitation fraction to 0.
          precip_frac(k) = zero

       endif

    enddo ! Overall precipitation fraction loop: k = nz, 1, -1.


    !!! Find precipitation fraction within each PDF component.
    !
    ! The overall precipitation fraction, f_p, is given by the equation:
    !
    ! f_p = a * f_p(1) + ( 1 - a ) * f_p(2);
    !
    ! where a is the mixture fraction (weight of PDF component 1), f_p(1) is
    ! the precipitation fraction within PDF component 1, and f_p(2) is the
    ! precipitation fraction within PDF component 2.  Overall precipitation
    ! fraction is found according the method above, and mixture fraction is
    ! already determined, leaving f_p(1) and f_p(2) to be solved for.  The
    ! values for f_p(1) and f_p(2) must satisfy the above equation.
    if ( .false. ) then

       !!! Find precipitation fraction within PDF component 1.
       ! The method used to find overall precipitation fraction will also be to
       ! find precipitation fraction within PDF component 1.  In order to do so,
       ! it is assumed (poorly) that PDF component 1 overlaps PDF component 1 at
       ! every vertical level in the vertical profile.
       do k = nz, 1, -1

          ! The weighted precipitation fraction (PDF component 1) is the
          ! greatest value of the product of mixture fraction and cloud fraction
          ! (PDF component 1) at or above a vertical level.
          if ( k < nz ) then
             weighted_pfrac1(k) = max( weighted_pfrac1(k+1), &
                                       mixt_frac(k) * cloud_frac1(k) )
          else  ! k = nz
             weighted_pfrac1(k) = mixt_frac(k) * cloud_frac1(k)
          endif

          precip_frac_1(k) = weighted_pfrac1(k) / mixt_frac(k)

          ! Special cases for precip_frac_1.
          if ( precip_frac_1(k) > one ) then

             ! Using the above method, it is possible for precip_frac_1 to be
             ! greater than 1.  For example, the mixture fraction at level k+1
             ! is 0.10 and the cloud_frac1 at level k+1 is 1, resulting in a
             ! weighted_pfrac1 of 0.10.  This product is greater than the
             ! product of mixt_frac and cloud_frac1 at level k.  The mixture
             ! fraction at level k is 0.05, resulting in a precip_frac_1 of 2.
             ! The value of precip_frac_1 is limited at 1.  The leftover
             ! precipitation fraction (a result of the decreasing weight of PDF
             ! component 1 between the levels) is applied to PDF component 2.
             precip_frac_1(k) = one

          elseif ( ( rr1(k) > rr_tol .or. Nr1(k) > Nr_tol ) &
                   .and. precip_frac_1(k) <= precip_frac_tol ) then

             ! In a scenario where we find rain in the 1st PDF component at this
             ! grid level, but no cloud in the 1st PDF component at or above
             ! this grid level, set precipitation fraction (in the 1st PDF
             ! component) to a minimum threshold value.
             precip_frac_1(k) = precip_frac_tol

          elseif ( ( rr1(k) <= rr_tol .and. Nr1(k) <= Nr_tol ) &
                   .and. precip_frac_1(k) <= precip_frac_tol ) then

             ! Mean rain water mixing ratio and mean rain drop concentration in
             ! the 1st PDF component are both less than their respective
             ! tolerance amounts.  They are both considered to have values of 0.
             ! There is not any rain in the 1st PDF component at this grid
             ! level.  There is also no cloud at or above this grid level, so
             ! set precipitation fraction (in the 1st PDF component) to 0.
             precip_frac_1(k) = zero

          endif

       enddo ! Precipitation fraction (1st PDF component) loop: k = nz, 1, -1.


       !!! Find precipitation fraction within PDF component 2.
       ! The equation for precipitation fraction within PDF component 2 is:
       !
       ! f_p(2) = ( f_p - a * f_p(1) ) / ( 1 - a );
       !
       ! given the overall precipitation fraction, f_p (calculated above), the
       ! precipitation fraction within PDF component 1, f_p(1) (calculated
       ! above), and mixture fraction, a.  Any leftover precipitation fraction
       ! from precip_frac_1 will be included in this calculation of
       ! precip_frac_2.
       do k = 1, nz, 1

          precip_frac_2(k) &
          = ( precip_frac(k) - mixt_frac(k) * precip_frac_1(k) ) &
            / ( one - mixt_frac(k) )

          ! Special cases for precip_frac_2.
          if ( precip_frac_2(k) > one ) then

             ! Again, it is possible for precip_frac_2 to be greater than 1.
             ! For example, the mixture fraction at level k+1 is 0.10 and the
             ! cloud_frac1 at level k+1 is 1, resulting in a weighted_pfrac1 of
             ! 0.10.  This product is greater than the product of mixt_frac and
             ! cloud_frac1 at level k.  Additionally, precip_frac (overall) is 1
             ! for level k.  The mixture fraction at level k is 0.5, resulting
             ! in a precip_frac_1 of 0.2.  Using the above equation,
             ! precip_frac_2 is calculated to be 1.8.  The value of
             ! precip_frac_2 is limited at 1.  The leftover precipitation
             ! fraction (as a result of the increasing weight of component 1
             ! between the levels) is applied to PDF component 1.
             precip_frac_2(k) = one

             ! Recalculate the precipitation fraction in PDF component 1.
             precip_frac_1(k) &
             = ( precip_frac(k) - ( one - mixt_frac(k) ) * precip_frac_2(k) ) &
               / mixt_frac(k)

             ! Double check for errors in PDF component 1.
             if ( precip_frac_1(k) > one ) then
                precip_frac_1(k) = one
             elseif ( ( rr1(k) > rr_tol .or. Nr1(k) > Nr_tol ) &
                      .and. precip_frac_1(k) <= precip_frac_tol ) then
                precip_frac_1(k) = precip_frac_tol
             elseif ( ( rr1(k) <= rr_tol .and. Nr1(k) <= Nr_tol ) &
                      .and. precip_frac_1(k) <= precip_frac_tol ) then
                precip_frac_1(k) = zero
             endif

          elseif ( ( rr2(k) > rr_tol .or. Nr2(k) > Nr_tol ) &
                   .and. precip_frac_2(k) <= precip_frac_tol ) then

             ! In a scenario where we find rain in the 2nd PDF component at this
             ! grid level, but no cloud in the 2nd PDF component at or above
             ! this grid level, set precipitation fraction (in the 2nd PDF
             ! component) to a minimum threshold value.
             precip_frac_2(k) = precip_frac_tol

          elseif ( ( rr2(k) <= rr_tol .and. Nr2(k) <= Nr_tol ) &
                   .and. precip_frac_2(k) <= precip_frac_tol ) then

             ! Mean rain water mixing ratio and mean rain drop concentration in
             ! the 2nd PDF component are both less than their respective
             ! tolerance amounts.  They are both considered to have values of 0.
             ! There is not any rain in the 2nd PDF component at this grid
             ! level.  There is also no cloud at or above this grid level, so
             ! set precipitation fraction (in the 2nd PDF component) to 0.
             precip_frac_2(k) = zero

          endif

       enddo ! Precipitation fraction (2nd PDF component) loop: k = 1, nz, 1.


    else  ! .true.

       ! Precipitation fraction in each PDF component is based on mean rain
       ! water mixing ratio in each PDF component.  The ratio it is based on is:
       !
       ! rr1/f_p(1) = rr2/f_p(2);
       !
       ! which can be rewritten as:
       !
       ! f_p(2)/f_p(1) = rr2/rr1.
       !
       ! Since overall precipitation fraction is given by the equation:
       !
       ! f_p = a f_p(1) + (1-a) f_p(2);
       !
       ! it can be rewritten as:
       !
       ! f_p = f_p(1) ( a + (1-a) f_p(2)/f_p(1) ).
       !
       ! Substituting the ratio rr2/rr1 for the ratio f_p(2)/f_p(1), the above
       ! equation can be solved for f_p(1):
       !
       ! f_p(1) = f_p / ( a + (1-a) rr2/rr1 ).
       !
       ! Then, f_p(2) can be solved for according to the equation:
       !
       ! f_p(2) = ( f_p - a f_p(1) ) / (1-a).
       do k = 1, nz, 1

          if ( ( rr1(k) <= rr_tol .and. Nr1(k) <= Nr_tol ) &
               .and. ( rr2(k) <= rr_tol .and. Nr2(k) <= Nr_tol ) ) then

             ! There is no rain in each PDF component.  Precipitation fraction
             ! within each component is set to 0.
             precip_frac_1(k) = zero
             precip_frac_2(k) = zero

          elseif ( ( rr1(k) > rr_tol .or. Nr1(k) > Nr_tol ) &
                   .and. ( rr2(k) <= rr_tol .and. Nr2(k) <= Nr_tol ) ) then

             ! All the rain is within the 1st PDF component.
             precip_frac_1(k) = precip_frac(k) / mixt_frac(k)
             precip_frac_2(k) = zero

          elseif ( ( rr2(k) > rr_tol .or. Nr2(k) > Nr_tol ) &
                   .and. ( rr1(k) <= rr_tol .and. Nr1(k) <= Nr_tol ) ) then

             ! All the rain is within the 2nd PDF component.
             precip_frac_1(k) = zero
             precip_frac_2(k) = precip_frac(k) / ( one - mixt_frac(k) )

          else ! rr1(k) > rr_tol or Nr1(k) > Nr_tol
               ! AND rr2(k) > rr_tol or Nr2(k) > Nr_tol

             ! Rain within both PDF components.

             !!! Find precipitation fraction within PDF component 1.
             if ( rr1(k) > rr_tol ) then
                precip_frac_1(k) &
                = precip_frac(k) &
                  / ( mixt_frac(k) + ( one - mixt_frac(k) ) * rr2(k)/rr1(k) )
             else ! Nr1 > Nr_tol
                precip_frac_1(k) &
                = precip_frac(k) &
                  / ( mixt_frac(k) + ( one - mixt_frac(k) ) * Nr2(k)/Nr1(k) )
             endif

             ! Using the above method, it is possible for precip_frac_1 to be
             ! greater than 1.  The value of precip_frac_1 is limited at 1.
             if ( precip_frac_1(k) > one ) then
                precip_frac_1(k) = one
             endif

             !!! Find precipitation fraction within PDF component 2.
             precip_frac_2(k) &
             = ( precip_frac(k) - mixt_frac(k) *  precip_frac_1(k) ) &
               / ( one - mixt_frac(k) )

             ! Using the above method, it is possible for precip_frac_2 to be
             ! greater than 1.  The value of precip_frac_2 is limited at 1.
             if ( precip_frac_2(k) > one ) then

                precip_frac_2(k) = one

                ! Recalculate the precipitation fraction in PDF component 1.
                precip_frac_1(k) &
                = ( precip_frac(k) &
                    - ( one - mixt_frac(k) ) * precip_frac_2(k) ) &
                  / mixt_frac(k)

             endif

          endif


          ! Special cases for PDF component 1.
          if ( ( rr1(k) > rr_tol .or. Nr1(k) > Nr_tol ) &
               .and. precip_frac_1(k) <= precip_frac_tol ) then

             ! In a scenario where we find rain in the 1st PDF component at this
             ! grid level, but no cloud in the 1st PDF component at or above
             ! this grid level, set precipitation fraction (in the 1st PDF
             ! component) to a minimum threshold value.
             precip_frac_1(k) = precip_frac_tol

          elseif ( ( rr1(k) <= rr_tol .and. Nr1(k) <= Nr_tol ) &
                   .and. precip_frac_1(k) <= precip_frac_tol ) then

             ! Mean rain water mixing ratio and mean rain drop concentration in
             ! the 1st PDF component are both less than their respective
             ! tolerance amounts.  They are both considered to have values of 0.
             ! There is not any rain in the 1st PDF component at this grid
             ! level.  There is also no cloud at or above this grid level, so
             ! set precipitation fraction (in the 1st PDF component) to 0.
             precip_frac_1(k) = zero

          endif


          ! Special cases for PDF component 2.
          if ( ( rr2(k) > rr_tol .or. Nr2(k) > Nr_tol ) &
               .and. precip_frac_2(k) <= precip_frac_tol ) then

             ! In a scenario where we find rain in the 2nd PDF component at this
             ! grid level, but no cloud in the 2nd PDF component at or above
             ! this grid level, set precipitation fraction (in the 2nd PDF
             ! component) to a minimum threshold value.
             precip_frac_2(k) = precip_frac_tol

          elseif ( ( rr2(k) <= rr_tol .and. Nr2(k) <= Nr_tol ) &
                   .and. precip_frac_2(k) <= precip_frac_tol ) then

             ! Mean rain water mixing ratio and mean rain drop concentration in
             ! the 2nd PDF component are both less than their respective
             ! tolerance amounts.  They are both considered to have values of 0.
             ! There is not any rain in the 2nd PDF component at this grid
             ! level.  There is also no cloud at or above this grid level, so
             ! set precipitation fraction (in the 2nd PDF component) to 0.
             precip_frac_2(k) = zero

          endif


       enddo ! Component precipitation fraction loop: k = 1, nz, 1.


    endif ! Select component precipitation fraction method.


    return

  end subroutine precip_fraction

  !=============================================================================
  subroutine compute_mean_stdev( rcm, rrainm, Nrm, Ncnm, &                     ! Intent(in)
                                 rr1, rr2, Nr1, Nr2, rc1, &                    ! Intent(in)
                                 rc2, cloud_frac1, cloud_frac2, &              ! Intent(in)
                                 precip_frac_1, precip_frac_2, &               ! Intent(in)
                                 wpsp, wprrp_ip, wpNrp_ip, &                   ! Intent(in)
                                 wpNcnp, stdev_w, mixt_frac, &                 ! Intent(in)
                                 pdf_params, &                                 ! Intent(in)
                                 mu_w_1, mu_w_2, mu_s_1, mu_s_2, mu_t_1, &     ! Intent(out)
                                 mu_t_2, mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, & ! Intent(out)
                                 mu_Ncn_1, mu_Ncn_2, sigma_w_1, sigma_w_2, &   ! Intent(out)
                                 sigma_s_1, sigma_s_2, sigma_t_1, sigma_t_2, & ! Intent(out)
                                 sigma_rr_1, sigma_rr_2, sigma_Nr_1, &         ! Intent(out)
                                 sigma_Nr_2, sigma_Ncn_1, sigma_Ncn_2 )        ! Intent(out)
       
    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        rc_tol,       & ! Constant(s)
        rr_tol,       &
        Nr_tol,       &
        Ncn_tol,      &
        w_tol,        & ! [m/s]
        s_mellor_tol, & ! [kg/kg]
        one,          &
        zero

    use parameters_microphys, only: &
        rrp2_on_rrm2_cloud,     & ! Variable(s)
        rrp2_on_rrm2_below,     &
        Nrp2_on_Nrm2_cloud,     &
        Nrp2_on_Nrm2_below,     &
        Ncnp2_on_Ncnm2_cloud,   &
        l_fix_s_t_correlations

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    use pdf_parameter_module, only: &
        pdf_parameter  ! Variable(s) type

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      rcm,           & ! Mean cloud water mixing ratio (overall)         [kg/kg]
      rrainm,        & ! Mean rain water mixing ratio (overall)          [kg/kg]
      Nrm,           & ! Mean rain drop concentration (overall)         [num/kg]
      Ncnm,          & ! Mean cloud nuclei concentration                [num/kg]
      rr1,           & ! Mean rain water mixing ratio (1st PDF comp.)    [kg/kg]
      rr2,           & ! Mean rain water mixing ratio (2nd PDF comp.)    [kg/kg]
      Nr1,           & ! Mean rain drop concentration (1st PDF comp.)   [num/kg]
      Nr2,           & ! Mean rain drop concentration (2nd PDF comp.)   [num/kg]
      rc1,           & ! Mean of r_c (1st PDF component)                 [kg/kg]
      rc2,           & ! Mean of r_c (2nd PDF component)                 [kg/kg]
      cloud_frac1,   & ! Cloud fraction (1st PDF component)                  [-]
      cloud_frac2,   & ! Cloud fraction (2nd PDF component)                  [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component)          [-]
      precip_frac_2, & ! Precipitation fraction (2nd PDF component)          [-]
      wpsp,          & ! Covariance of w and s                      [(m/s)kg/kg]
      wprrp_ip,      & ! Covariance of w and r_r (overall) ip       [(m/s)kg/kg]
      wpNrp_ip,      & ! Covariance of w and N_r (overall) ip      [(m/s)num/kg]
      wpNcnp,        & ! Covariance of w and N_cn                  [(m/s)num/kg]
      stdev_w,       & ! Standard deviation of w                           [m/s]
      mixt_frac        ! Mixture fraction                                    [-]

    type(pdf_parameter), intent(in) :: &
      pdf_params    ! PDF parameters                                [units vary]

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
      mu_w_1,      & ! Mean of w (1st PDF component)                       [m/s]
      mu_w_2,      & ! Mean of w (2nd PDF component)                       [m/s]
      mu_s_1,      & ! Mean of s (1st PDF component)                     [kg/kg]
      mu_s_2,      & ! Mean of s (2nd PDF component)                     [kg/kg]
      mu_t_1,      & ! Mean of t (1st PDF component)                     [kg/kg]
      mu_t_2,      & ! Mean of t (2nd PDF component)                     [kg/kg]
      mu_rr_1,     & ! Mean of rr (1st PDF component) in-precip (ip)     [kg/kg]
      mu_rr_2,     & ! Mean of rr (2nd PDF component) ip                 [kg/kg]
      mu_Nr_1,     & ! Mean of Nr (1st PDF component) ip                [num/kg]
      mu_Nr_2,     & ! Mean of Nr (2nd PDF component) ip                [num/kg]
      mu_Ncn_1,    & ! Mean of Ncn (1st PDF component)                  [num/kg]
      mu_Ncn_2,    & ! Mean of Ncn (2nd PDF component)                  [num/kg]
      sigma_w_1,   & ! Standard deviation of w (1st PDF component)         [m/s]
      sigma_w_2,   & ! Standard deviation of w (2nd PDF component)         [m/s]
      sigma_s_1,   & ! Standard deviation of s (1st PDF component)       [kg/kg]
      sigma_s_2,   & ! Standard deviation of s (2nd PDF component)       [kg/kg]
      sigma_t_1,   & ! Standard deviation of t (1st PDF component)       [kg/kg]
      sigma_t_2,   & ! Standard deviation of t (2nd PDF component)       [kg/kg]
      sigma_rr_1,  & ! Standard deviation of rr (1st PDF component) ip   [kg/kg]
      sigma_rr_2,  & ! Standard deviation of rr (2nd PDF component) ip   [kg/kg]
      sigma_Nr_1,  & ! Standard deviation of Nr (1st PDF component) ip  [num/kg]
      sigma_Nr_2,  & ! Standard deviation of Nr (2nd PDF component) ip  [num/kg]
      sigma_Ncn_1, & ! Standard deviation of Ncn (1st PDF component)    [num/kg]
      sigma_Ncn_2    ! Standard deviation of Ncn (2nd PDF component)    [num/kg]

    ! Local Variables

    ! The component correlations of w and r_t and the component correlations of
    ! w and theta_l are both set to be 0 within the CLUBB model code.  In other
    ! words, w and r_t (theta_l) have overall covariance w'r_t' (w'theta_l'),
    ! but the single component covariance and correlation are defined to be 0.
    ! Likewise, the single component correlation and covariance of w and s, as
    ! well as w and t, are defined to be 0.
    logical, parameter :: &
      l_follow_CLUBB_PDF_standards = .true.

    ! Prescribed parameters are set to in-cloud or outside-cloud (below-cloud)
    ! values based on whether or not cloud water mixing ratio has a value of at
    ! least rc_tol.  However, this does not take into account the amount of
    ! cloudiness in a component, just whether or not there is any cloud in the
    ! component.  The option l_interp_prescribed_params allows for an
    ! interpolated value between the in-cloud and below-cloud parameter value
    ! based on the component cloud fraction.
    logical, parameter :: &
      l_interp_prescribed_params = .false.

    real( kind = core_rknd ) :: &
      rrp2_on_rrm2,    & ! Ratio of < r_r'^2 > to < r_r >^2              [-]
      Nrp2_on_Nrm2,    & ! Ratio of < N_r'^2 > to < N_r >^2              [-]
      Ncnp2_on_Ncnm2,  & ! Ratio of < N_cn'^2 > to < N_cn >^2            [-]
      s_mellor_m,      & ! Mean of s_mellor                              [kg/kg]
      stdev_s_mellor,  & ! Standard deviation of s_mellor                [kg/kg]
      corr_ws,         & ! Correlation between w and s                   [-]
      corr_wrr,        & ! Correlation between w and rr ip               [-]
      corr_wNr,        & ! Correlation between w and Nr ip               [-]
      corr_wNcn          ! Correlation between w and Ncn                 [-]


    !!! Enter the PDF parameters.

    !!! Means.

    ! Mean of vertical velocity, w, in PDF component 1.
    mu_w_1 = pdf_params%w1

    ! Mean of vertical velocity, w, in PDF component 2.
    mu_w_2 = pdf_params%w2

    ! Mean of extended liquid water mixing ratio, s, in PDF component 1.
    mu_s_1 = pdf_params%s1

    ! Mean of extended liquid water mixing ratio, s, in PDF component 2.
    mu_s_2 = pdf_params%s2

    ! Mean of t in PDF component 1.
    ! Set the component mean values of t to 0.
    ! The component mean values of t are not important.  They can be set to
    ! anything.  They cancel out in the model code.  However, the best thing to
    ! do is to set them to 0 and avoid any kind of numerical error.
    mu_t_1 = zero

    ! Mean of t in PDF component 2.
    ! Set the component mean values of t to 0.
    ! The component mean values of t are not important.  They can be set to
    ! anything.  They cancel out in the model code.  However, the best thing to
    ! do is to set them to 0 and avoid any kind of numerical error.
    mu_t_2 = zero

    ! Mean of in-precip rain water mixing ratio in PDF component 1.
    if ( rr1 > rr_tol ) then
       mu_rr_1 = rr1 / precip_frac_1
    else
       ! Mean in-precip rain water mixing ratio in PDF component 1 is less than
       ! the tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 1st PDF component at this grid level.
       mu_rr_1 = zero
    endif

    ! Mean of in-precip rain water mixing ratio in PDF component 2.
    if ( rr2 > rr_tol ) then
       mu_rr_2 = rr2 / precip_frac_2
    else
       ! Mean in-precip rain water mixing ratio in PDF component 2 is less than
       ! the tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 2nd PDF component at this grid level.
       mu_rr_2 = zero
    endif

    ! Mean of in-precip rain drop concentration in PDF component 1.
    if ( Nr1 > Nr_tol ) then
       mu_Nr_1 = Nr1 / precip_frac_1
    else
       ! Mean in-precip rain drop concentration in PDF component 1 is less than
       ! the tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 1st PDF component at this grid level.
       mu_Nr_1 = zero
    endif

    ! Mean of in-precip rain drop concentration in PDF component 2.
    if ( Nr2 > Nr_tol ) then
       mu_Nr_2 = Nr2 / precip_frac_2
    else
       ! Mean in-precip rain drop concentration in PDF component 2 is less than
       ! the tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 2nd PDF component at this grid level.
       mu_Nr_2 = zero
    endif

    ! Mean of cloud nuclei concentration in PDF component 1.
    if ( Ncnm > Ncn_tol ) then
       mu_Ncn_1 = Ncnm
    else
       ! Mean cloud nuclei concentration is less than the tolerance amount.  It
       ! is considered to have a value of 0.  There are not any cloud nuclei or
       ! cloud at this grid level.
       mu_Ncn_1 = zero
    endif

    ! Mean of cloud nuclei concentration in PDF component 2.
    if ( Ncnm > Ncn_tol ) then
       mu_Ncn_2 = Ncnm
    else
       ! Mean cloud nuclei concentration is less than the tolerance amount.  It
       ! is considered to have a value of 0.  There are not any cloud nuclei or
       ! cloud at this grid level.
       mu_Ncn_2 = zero
    endif


    !!! Standard deviations.

    ! Standard deviation of vertical velocity, w, in PDF component 1.
    sigma_w_1 = sqrt( pdf_params%varnce_w1 )

    ! Standard deviation of vertical velocity, w, in PDF component 2.
    sigma_w_2 = sqrt( pdf_params%varnce_w2 )

    ! Standard deviation of extended liquid water mixing ratio, s,
    ! in PDF component 1.
    sigma_s_1 = pdf_params%stdev_s1

    ! Standard deviation of extended liquid water mixing ratio, s,
    ! in PDF component 2.
    sigma_s_2 = pdf_params%stdev_s2

    ! Standard deviation of t in PDF component 1.
    sigma_t_1 = pdf_params%stdev_t1

    ! Standard deviation of t in PDF component 2.
    sigma_t_2 = pdf_params%stdev_t2

    ! Set up the values of the statistical correlations and variances.  Since we
    ! currently do not have enough variables to compute the correlations and
    ! variances directly, we have obtained these values by analyzing LES runs of
    ! certain cases.  We have divided those results into an inside-cloud average
    ! and an outside-cloud (or below-cloud) average.  This coding leaves the
    ! software architecture in place in case we ever have the variables in place
    ! to compute these values directly.  It also allows us to use separate
    ! inside-cloud and outside-cloud parameter values.
    ! Brian Griffin; February 3, 2007.

    ! Standard deviation of in-precip rain water mixing ratio
    ! in PDF component 1.
    if ( rr1 > rr_tol ) then
       if ( l_interp_prescribed_params ) then
          sigma_rr_1 = sqrt( cloud_frac1 * rrp2_on_rrm2_cloud &
                             + ( one - cloud_frac1 ) * rrp2_on_rrm2_below ) &
                       * mu_rr_1
       else
          if ( rc1 > rc_tol ) then
             sigma_rr_1 = sqrt( rrp2_on_rrm2_cloud ) * mu_rr_1
          else
             sigma_rr_1 = sqrt( rrp2_on_rrm2_below ) * mu_rr_1
          endif
       endif
    else
       ! Mean in-precip rain water mixing ratio in PDF component 1 is less than
       ! the tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 1st PDF component at this grid level.  The standard
       ! deviation is simply 0 since rain water mixing ratio does not vary in
       ! this component at this grid level.
       sigma_rr_1 = zero
    endif

    ! Standard deviation of in-precip rain water mixing ratio
    ! in PDF component 2.
    if ( rr2 > rr_tol ) then
       if ( l_interp_prescribed_params ) then
          sigma_rr_2 = sqrt( cloud_frac2 * rrp2_on_rrm2_cloud &
                             + ( one - cloud_frac2 ) * rrp2_on_rrm2_below ) &
                       * mu_rr_2
       else
          if ( rc2 > rc_tol ) then
             sigma_rr_2 = sqrt( rrp2_on_rrm2_cloud ) * mu_rr_2
          else
             sigma_rr_2 = sqrt( rrp2_on_rrm2_below ) * mu_rr_2
          endif
       endif
    else
       ! Mean in-precip rain water mixing ratio in PDF component 2 is less than
       ! the tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 2nd PDF component at this grid level.  The standard
       ! deviation is simply 0 since rain water mixing ratio does not vary in
       ! this component at this grid level.
       sigma_rr_2 = zero
    endif

    ! Standard deviation of in-precip rain drop concentration
    ! in PDF component 1.
    if ( Nr1 > Nr_tol ) then
       if ( l_interp_prescribed_params ) then
          sigma_Nr_1 = sqrt( cloud_frac1 * Nrp2_on_Nrm2_cloud &
                             + ( one - cloud_frac1 ) * Nrp2_on_Nrm2_below ) &
                       * mu_Nr_1
       else
          if ( rc1 > rc_tol ) then
             sigma_Nr_1 = sqrt( Nrp2_on_Nrm2_cloud ) * mu_Nr_1
          else
             sigma_Nr_1 = sqrt( Nrp2_on_Nrm2_below ) * mu_Nr_1
          endif
       endif
    else
       ! Mean in-precip rain drop concentration in PDF component 1 is less than
       ! the tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 1st PDF component at this grid level.  The standard
       ! deviation is simply 0 since rain drop concentration does not vary in
       ! this component at this grid level.
       sigma_Nr_1 = zero
    endif

    ! Standard deviation of in-precip rain drop concentration
    ! in PDF component 2.
    if ( Nr2 > Nr_tol ) then
       if ( l_interp_prescribed_params ) then
          sigma_Nr_2 = sqrt( cloud_frac2 * Nrp2_on_Nrm2_cloud &
                             + ( one - cloud_frac2 ) * Nrp2_on_Nrm2_below ) &
                       * mu_Nr_2
       else
          if ( rc2 > rc_tol ) then
             sigma_Nr_2 = sqrt( Nrp2_on_Nrm2_cloud ) * mu_Nr_2
          else
             sigma_Nr_2 = sqrt( Nrp2_on_Nrm2_below ) * mu_Nr_2
          endif
       endif
    else
       ! Mean in-precip rain drop concentration in PDF component 2 is less than
       ! the tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 2nd PDF component at this grid level.  The standard
       ! deviation is simply 0 since rain drop concentration does not vary in
       ! this component at this grid level.
       sigma_Nr_2 = zero
    endif

    ! Standard deviation of cloud nuclei concentration in PDF component 1.
    if ( Ncnm > Ncn_tol ) then
       sigma_Ncn_1 = sqrt( Ncnp2_on_Ncnm2_cloud ) * mu_Ncn_1
    else
       ! Mean cloud nuclei concentration is less than the tolerance amount.  It
       ! is considered to have a value of 0.  There are not any cloud nuclei or
       ! cloud at this grid level.  The standard deviation is simply 0 since
       ! cloud nuclei concentration does not vary at this grid level.
       sigma_Ncn_1 = zero
    endif

    ! Standard deviation of cloud nuclei concentration in PDF component 2.
    if ( Ncnm > Ncn_tol ) then
       sigma_Ncn_2 = sqrt( Ncnp2_on_Ncnm2_cloud ) * mu_Ncn_2
    else
       ! Mean cloud nuclei concentration is less than the tolerance amount.  It
       ! is considered to have a value of 0.  There are not any cloud nuclei or
       ! cloud at this grid level.  The standard deviation is simply 0 since
       ! cloud nuclei concentration does not vary at this grid level.
       sigma_Ncn_2 = zero
    endif

    return

  end subroutine compute_mean_stdev

!=============================================================================
  subroutine compute_corr( rcm, rrainm, Nrm, Ncnm, &               ! Intent(in)
                           rr1, rr2, Nr1, Nr2, rc1, &              ! Intent(in)
                           rc2, cloud_frac1, cloud_frac2, &        ! Intent(in)
                           precip_frac_1, precip_frac_2, &         ! Intent(in)
                           wpsp, wprrp_ip, wpNrp_ip, &             ! Intent(in)
                           wpNcnp, stdev_w, mixt_frac, &           ! Intent(in)
                           sigma_rr_1, sigma_Nr_1, sigma_Ncn_1, &  ! Intent(in)
                           pdf_params, &                           ! Intent(in)
                           corr_ws_1, corr_ws_2, corr_wrr_1, &     ! Intent(out)
                           corr_wrr_2, corr_wNr_1, corr_wNr_2, &   ! Intent(out)
                           corr_wNcn_1, corr_wNcn_2, corr_st_1, &  ! Intent(out)
                           corr_st_2, corr_srr_1, corr_srr_2, &    ! Intent(out)
                           corr_sNr_1, corr_sNr_2, corr_sNcn_1, &  ! Intent(out)
                           corr_sNcn_2, corr_trr_1, corr_trr_2, &  ! Intent(out)
                           corr_tNr_1, corr_tNr_2, corr_tNcn_1, &  ! Intent(out)
                           corr_tNcn_2, corr_rrNr_1, corr_rrNr_2 ) ! Intent(out)
       
    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        rc_tol,       & ! Constant(s)
        rr_tol,       &
        Nr_tol,       &
        Ncn_tol,      &
        w_tol,        & ! [m/s]
        s_mellor_tol, & ! [kg/kg]
        one,          &
        zero

    use parameters_microphys, only: &
        l_fix_s_t_correlations ! Variable(s)

    use KK_fixed_correlations, only: &
        corr_sw_NN_cloud,   & ! Variable(s)
        corr_wrr_NL_cloud,  &
        corr_wNr_NL_cloud,  &
        corr_wNcn_NL_cloud, &
        corr_st_NN_cloud,   &
        corr_srr_NL_cloud,  &
        corr_sNr_NL_cloud,  &
        corr_sNcn_NL_cloud, &
        corr_trr_NL_cloud,  &
        corr_tNr_NL_cloud,  &
        corr_tNcn_NL_cloud, &
        corr_rrNr_LL_cloud

    use KK_fixed_correlations, only: &
        corr_sw_NN_below,   & ! Variable(s)
        corr_wrr_NL_below,  &
        corr_wNr_NL_below,  &
        corr_wNcn_NL_below, &
        corr_st_NN_below,   &
        corr_srr_NL_below,  &
        corr_sNr_NL_below,  &
        corr_sNcn_NL_below, &
        corr_trr_NL_below,  &
        corr_tNr_NL_below,  &
        corr_tNcn_NL_below, &
        corr_rrNr_LL_below

    use model_flags, only: &
        l_calc_w_corr

    use diagnose_correlations_module, only: &
        calc_mean,        & ! Procedure(s)
        calc_w_corr

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    use pdf_parameter_module, only: &
        pdf_parameter  ! Variable(s) type

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      rcm,           & ! Mean cloud water mixing ratio (overall)         [kg/kg]
      rrainm,        & ! Mean rain water mixing ratio (overall)          [kg/kg]
      Nrm,           & ! Mean rain drop concentration (overall)         [num/kg]
      Ncnm,          & ! Mean cloud nuclei concentration                [num/kg]
      rr1,           & ! Mean rain water mixing ratio (1st PDF comp.)    [kg/kg]
      rr2,           & ! Mean rain water mixing ratio (2nd PDF comp.)    [kg/kg]
      Nr1,           & ! Mean rain drop concentration (1st PDF comp.)   [num/kg]
      Nr2,           & ! Mean rain drop concentration (2nd PDF comp.)   [num/kg]
      rc1,           & ! Mean of r_c (1st PDF component)                 [kg/kg]
      rc2,           & ! Mean of r_c (2nd PDF component)                 [kg/kg]
      cloud_frac1,   & ! Cloud fraction (1st PDF component)                  [-]
      cloud_frac2,   & ! Cloud fraction (2nd PDF component)                  [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component)          [-]
      precip_frac_2, & ! Precipitation fraction (2nd PDF component)          [-]
      wpsp,          & ! Covariance of w and s                      [(m/s)kg/kg]
      wprrp_ip,      & ! Covariance of w and r_r (overall) ip       [(m/s)kg/kg]
      wpNrp_ip,      & ! Covariance of w and N_r (overall) ip      [(m/s)num/kg]
      wpNcnp,        & ! Covariance of w and N_cn                  [(m/s)num/kg]
      stdev_w,       & ! Standard deviation of w                           [m/s]
      mixt_frac        ! Mixture fraction                                    [-]

    real( kind = core_rknd ), intent(in) :: &
      sigma_rr_1, &
      sigma_Nr_1, &
      sigma_Ncn_1

    type(pdf_parameter), intent(in) :: &
      pdf_params    ! PDF parameters                                [units vary]

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
      corr_ws_1,   & ! Correlation between w and s (1st PDF component)       [-]
      corr_ws_2,   & ! Correlation between w and s (2nd PDF component)       [-]
      corr_wrr_1,  & ! Correlation between w and rr (1st PDF component) ip   [-]
      corr_wrr_2,  & ! Correlation between w and rr (2nd PDF component) ip   [-]
      corr_wNr_1,  & ! Correlation between w and Nr (1st PDF component) ip   [-]
      corr_wNr_2,  & ! Correlation between w and Nr (2nd PDF component) ip   [-]
      corr_wNcn_1, & ! Correlation between w and Ncn (1st PDF component)     [-]
      corr_wNcn_2, & ! Correlation between w and Ncn (2nd PDF component)     [-]
      corr_st_1,   & ! Correlation between s and t (1st PDF component)       [-]
      corr_st_2,   & ! Correlation between s and t (2nd PDF component)       [-]
      corr_srr_1,  & ! Correlation between s and rr (1st PDF component) ip   [-]
      corr_srr_2,  & ! Correlation between s and rr (2nd PDF component) ip   [-]
      corr_sNr_1,  & ! Correlation between s and Nr (1st PDF component) ip   [-]
      corr_sNr_2,  & ! Correlation between s and Nr (2nd PDF component) ip   [-]
      corr_sNcn_1, & ! Correlation between s and Ncn (1st PDF component)     [-]
      corr_sNcn_2, & ! Correlation between s and Ncn (2nd PDF component)     [-]
      corr_trr_1,  & ! Correlation between t and rr (1st PDF component) ip   [-]
      corr_trr_2,  & ! Correlation between t and rr (2nd PDF component) ip   [-]
      corr_tNr_1,  & ! Correlation between t and Nr (1st PDF component) ip   [-]
      corr_tNr_2,  & ! Correlation between t and Nr (2nd PDF component) ip   [-]
      corr_tNcn_1, & ! Correlation between t and Ncn (1st PDF component)     [-]
      corr_tNcn_2, & ! Correlation between t and Ncn (2nd PDF component)     [-]
      corr_rrNr_1, & ! Correlation between rr and Nr (1st PDF component) ip  [-]
      corr_rrNr_2    ! Correlation between rr and Nr (2nd PDF component) ip  [-]

    ! Local Variables

    ! The component correlations of w and r_t and the component correlations of
    ! w and theta_l are both set to be 0 within the CLUBB model code.  In other
    ! words, w and r_t (theta_l) have overall covariance w'r_t' (w'theta_l'),
    ! but the single component covariance and correlation are defined to be 0.
    ! Likewise, the single component correlation and covariance of w and s, as
    ! well as w and t, are defined to be 0.
    logical, parameter :: &
      l_follow_CLUBB_PDF_standards = .true.

    ! Prescribed parameters are set to in-cloud or outside-cloud (below-cloud)
    ! values based on whether or not cloud water mixing ratio has a value of at
    ! least rc_tol.  However, this does not take into account the amount of
    ! cloudiness in a component, just whether or not there is any cloud in the
    ! component.  The option l_interp_prescribed_params allows for an
    ! interpolated value between the in-cloud and below-cloud parameter value
    ! based on the component cloud fraction.
    logical, parameter :: &
      l_interp_prescribed_params = .false.

    real( kind = core_rknd ) :: &
      rrp2_on_rrm2,    & ! Ratio of < r_r'^2 > to < r_r >^2              [-]
      Nrp2_on_Nrm2,    & ! Ratio of < N_r'^2 > to < N_r >^2              [-]
      Ncnp2_on_Ncnm2,  & ! Ratio of < N_cn'^2 > to < N_cn >^2            [-]
      s_mellor_m,      & ! Mean of s_mellor                              [kg/kg]
      stdev_s_mellor,  & ! Standard deviation of s_mellor                [kg/kg]
      corr_ws,         & ! Correlation between w and s                   [-]
      corr_wrr,        & ! Correlation between w and rr ip               [-]
      corr_wNr,        & ! Correlation between w and Nr ip               [-]
      corr_wNcn          ! Correlation between w and Ncn                 [-]


    !!! Enter the PDF parameters.

    !!! Correlations

    ! Calculate correlations involving w by first calculating total covariances
    ! involving w (<w'r_r'>, etc.) using the down-gradient approximation.
    if ( l_calc_w_corr ) then

       s_mellor_m &
       = calc_mean( pdf_params%mixt_frac, pdf_params%s1, pdf_params%s2 )

       stdev_s_mellor &
       = sqrt( pdf_params%mixt_frac &
               * ( ( pdf_params%s1 - s_mellor_m )**2 &
                   + pdf_params%stdev_s1**2 ) &
             + ( 1 - pdf_params%mixt_frac ) &
               * ( ( pdf_params%s2 - s_mellor_m )**2 &
                   + pdf_params%stdev_s2**2 ) &
             )

       corr_ws &
       = calc_w_corr( wpsp, stdev_w, stdev_s_mellor, w_tol, s_mellor_tol )

       corr_wrr = calc_w_corr( wprrp_ip, stdev_w, sigma_rr_1, w_tol, rr_tol )

       corr_wNr = calc_w_corr( wpNrp_ip, stdev_w, sigma_Nr_1, w_tol, Nr_tol )

       corr_wNcn = calc_w_corr( wpNcnp, stdev_w, sigma_Ncn_1, w_tol, Ncn_tol )

    endif

    ! Correlation between w and s in PDF component 1.
    ! The component correlations of w and r_t and the component correlations of
    ! w and theta_l are both set to be 0 within the CLUBB model code.  In other
    ! words, w and r_t (theta_l) have overall covariance w'r_t' (w'theta_l'),
    ! but the single component covariance and correlation are defined to be 0.
    ! Likewise, the single component correlation and covariance of w and s, as
    ! well as w and t, are defined to be 0.
    if ( l_follow_CLUBB_PDF_standards ) then
       corr_ws_1 = zero
    else ! not following CLUBB PDF standards
       ! WARNING:  the standards used in the generation of the two-component
       !           CLUBB PDF are not being obeyed.  The use of this code is
       !           inconsistent with the rest of CLUBB's PDF.
       if ( l_calc_w_corr ) then
          corr_ws_1 = corr_ws
       else ! use prescribed parameter values
          if ( l_interp_prescribed_params ) then
             corr_ws_1 = cloud_frac1 * corr_sw_NN_cloud &
                         + ( one - cloud_frac1 ) * corr_sw_NN_below
          else
             if ( rc1 > rc_tol ) then
                corr_ws_1 = corr_sw_NN_cloud
             else
                corr_ws_1 = corr_sw_NN_below
             endif
          endif ! l_interp_prescribed_params
       endif ! l_calc_w_corr
    endif ! l_follow_CLUBB_PDF_standards

    ! Correlation between w and s in PDF component 2.
    ! The component correlations of w and r_t and the component correlations of
    ! w and theta_l are both set to be 0 within the CLUBB model code.  In other
    ! words, w and r_t (theta_l) have overall covariance w'r_t' (w'theta_l'),
    ! but the single component covariance and correlation are defined to be 0.
    ! Likewise, the single component correlation and covariance of w and s, as
    ! well as w and t, are defined to be 0.
    if ( l_follow_CLUBB_PDF_standards ) then
       corr_ws_2 = zero
    else ! not following CLUBB PDF standards
       ! WARNING:  the standards used in the generation of the two-component
       !           CLUBB PDF are not being obeyed.  The use of this code is
       !           inconsistent with the rest of CLUBB's PDF.
       if ( l_calc_w_corr ) then
          corr_ws_2 = corr_ws
       else ! use prescribed parameter values
          if ( l_interp_prescribed_params ) then
             corr_ws_2 = cloud_frac2 * corr_sw_NN_cloud &
                         + ( one - cloud_frac2 ) * corr_sw_NN_below
          else
             if ( rc2 > rc_tol ) then
                corr_ws_2 = corr_sw_NN_cloud
             else
                corr_ws_2 = corr_sw_NN_below
             endif
          endif ! l_interp_prescribed_params
       endif ! l_calc_w_corr
    endif ! l_follow_CLUBB_PDF_standards

    ! Correlation (in-precip) between w and r_r in PDF component 1.
    if ( rr1 > rr_tol ) then
       if ( l_calc_w_corr ) then
          corr_wrr_1 = corr_wrr
       else ! use prescribed parameter values
          if ( l_interp_prescribed_params ) then
             corr_wrr_1 = cloud_frac1 * corr_wrr_NL_cloud &
                          + ( one - cloud_frac1 ) * corr_wrr_NL_below
          else
             if ( rc1 > rc_tol ) then
                corr_wrr_1 = corr_wrr_NL_cloud
             else
                corr_wrr_1 = corr_wrr_NL_below
             endif
          endif ! l_interp_prescribed_params
       endif ! ! l_calc_w_corr
    else
       ! Mean in-precip rain water mixing ratio in PDF component 1 is less than
       ! the tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 1st PDF component at this grid level.  The
       ! correlations involving rain water mixing ratio in the 1st PDF component
       ! are 0 since rain water mixing ratio does not vary in this component at
       ! this grid level.
       corr_wrr_1 = zero
    endif

    ! Correlation (in-precip) between w and r_r in PDF component 2.
    if ( rr2 > rr_tol ) then
       if ( l_calc_w_corr ) then
          corr_wrr_2 = corr_wrr
       else ! use prescribed parameter values
          if ( l_interp_prescribed_params ) then
             corr_wrr_2 = cloud_frac2 * corr_wrr_NL_cloud &
                          + ( one - cloud_frac2 ) * corr_wrr_NL_below
          else
             if ( rc2 > rc_tol ) then
                corr_wrr_2 = corr_wrr_NL_cloud
             else
                corr_wrr_2 = corr_wrr_NL_below
             endif
          endif ! l_interp_prescribed_params
       endif ! l_calc_w_corr
    else
       ! Mean in-precip rain water mixing ratio in PDF component 2 is less than
       ! the tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 2nd PDF component at this grid level.  The
       ! correlations involving rain water mixing ratio in the 2nd PDF component
       ! are 0 since rain water mixing ratio does not vary in this component at
       ! this grid level.
       corr_wrr_2 = zero
    endif

    ! Correlation (in-precip) between w and N_r in PDF component 1.
    if ( Nr1 > Nr_tol ) then
       if ( l_calc_w_corr ) then
          corr_wNr_1 = corr_wNr
       else ! use prescribed parameter values
          if ( l_interp_prescribed_params ) then
             corr_wNr_1 = cloud_frac1 * corr_wNr_NL_cloud &
                          + ( one - cloud_frac1 ) * corr_wNr_NL_below
          else
             if ( rc1 > rc_tol ) then
                corr_wNr_1 = corr_wNr_NL_cloud
             else
                corr_wNr_1 = corr_wNr_NL_below
             endif
          endif ! l_interp_prescribed_params
       endif ! l_calc_w_corr
    else
       ! Mean in-precip rain drop concentration in PDF component 1 is less than
       ! the tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 1st PDF component at this grid level.  The
       ! correlations involving rain drop concentration in the 1st PDF component
       ! are 0 since rain drop concentration does not vary in this component at
       ! this grid level.
       corr_wNr_1 = zero
    endif

    ! Correlation (in-precip) between w and N_r in PDF component 2.
    if ( Nr2 > Nr_tol ) then
       if ( l_calc_w_corr ) then
          corr_wNr_2 = corr_wNr
       else ! use prescribed parameter values
          if ( l_interp_prescribed_params ) then
             corr_wNr_2 = cloud_frac2 * corr_wNr_NL_cloud &
                          + ( one - cloud_frac2 ) * corr_wNr_NL_below
          else
             if ( rc2 > rc_tol ) then
                corr_wNr_2 = corr_wNr_NL_cloud
             else
                corr_wNr_2 = corr_wNr_NL_below
             endif
          endif ! l_interp_prescribed_params
       endif ! l_calc_w_corr
    else
       ! Mean in-precip rain drop concentration in PDF component 2 is less than
       ! the tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 2nd PDF component at this grid level.  The
       ! correlations involving rain drop concentration in the 2nd PDF component
       ! are 0 since rain drop concentration does not vary in this component at
       ! this grid level.
       corr_wNr_2 = zero
    endif

    ! Correlation between w and N_cn in PDF component 1.
    if ( Ncnm > Ncn_tol ) then
       if ( l_calc_w_corr ) then
          corr_wNcn_1 = corr_wNcn
       else ! use prescribed parameter values
          corr_wNcn_1 = corr_wNcn_NL_cloud
       endif ! l_calc_w_corr
    else
       ! Mean cloud nuclei concentration is less than the tolerance amount.  It
       ! is considered to have a value of 0.  There are not any cloud nuclei or
       ! cloud at this grid level.  The correlations involving cloud nuclei
       ! concentration are 0 since cloud nuclei concentration does not vary at
       ! this grid level.
       corr_wNcn_1 = zero
    endif

   ! Correlation between w and N_cn in PDF component 2.
    if ( Ncnm > Ncn_tol ) then
       if ( l_calc_w_corr ) then
          corr_wNcn_2 = corr_wNcn
       else ! use prescribed parameter values
          corr_wNcn_2 = corr_wNcn_NL_cloud
       endif ! l_calc_w_corr
    else
       ! Mean cloud nuclei concentration is less than the tolerance amount.  It
       ! is considered to have a value of 0.  There are not any cloud nuclei or
       ! cloud at this grid level.  The correlations involving cloud nuclei
       ! concentration are 0 since cloud nuclei concentration does not vary at
       ! this grid level.
       corr_wNcn_2 = zero
    endif

    ! Correlation between s and t in PDF component 1.
    ! The PDF variables s and t result from a transformation of the PDF
    ! involving r_t and theta_l.  The correlation between s and t depends on the
    ! correlation between r_t and theta_l, as well as the variances of r_t and
    ! theta_l, and other factors.  The correlation between s and t is subject to
    ! change at every vertical level and model time step, and is calculated as
    ! part of the CLUBB PDF parameters.
    if ( .not. l_fix_s_t_correlations ) then
       corr_st_1 = pdf_params%corr_st_1
    else ! fix the correlation between s and t.
       ! WARNING:  this code is inconsistent with the rest of CLUBB's PDF.  This
       !           code is necessary because SIHLS is lazy and wussy, and only
       !           wants to declare correlation arrays at the start of the model
       !           run, rather than updating them throughout the model run.
       if ( l_interp_prescribed_params ) then
          corr_st_1 = cloud_frac1 * corr_st_NN_cloud &
                      + ( one - cloud_frac1 ) * corr_st_NN_below
       else
          if ( rc1 > rc_tol ) then
             corr_st_1 = corr_st_NN_cloud
          else
             corr_st_1 = corr_st_NN_below
          endif
       endif
    endif

    ! Correlation between s and t in PDF component 2.
    ! The PDF variables s and t result from a transformation of the PDF
    ! involving r_t and theta_l.  The correlation between s and t depends on the
    ! correlation between r_t and theta_l, as well as the variances of r_t and
    ! theta_l, and other factors.  The correlation between s and t is subject to
    ! change at every vertical level and model time step, and is calculated as
    ! part of the CLUBB PDF parameters.
    if ( .not. l_fix_s_t_correlations ) then
       corr_st_2 = pdf_params%corr_st_2
    else ! fix the correlation between s and t.
       ! WARNING:  this code is inconsistent with the rest of CLUBB's PDF.  This
       !           code is necessary because SIHLS is lazy and wussy, and only
       !           wants to declare correlation arrays at the start of the model
       !           run, rather than updating them throughout the model run.
       if ( l_interp_prescribed_params ) then
          corr_st_2 = cloud_frac2 * corr_st_NN_cloud &
                      + ( one - cloud_frac2 ) * corr_st_NN_below
       else
          if ( rc2 > rc_tol ) then
             corr_st_2 = corr_st_NN_cloud
          else
             corr_st_2 = corr_st_NN_below
          endif
       endif
    endif

    ! Correlation (in-precip) between s and r_r in PDF component 1.
    if ( rr1 > rr_tol ) then
       if ( l_interp_prescribed_params ) then
          corr_srr_1 = cloud_frac1 * corr_srr_NL_cloud &
                       + ( one - cloud_frac1 ) * corr_srr_NL_below
       else
          if ( rc1 > rc_tol ) then
             corr_srr_1 = corr_srr_NL_cloud
          else
             corr_srr_1 = corr_srr_NL_below
          endif
       endif
    else
       ! Mean in-precip rain water mixing ratio in PDF component 1 is less than
       ! the tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 1st PDF component at this grid level.  The
       ! correlations involving rain water mixing ratio in the 1st PDF component
       ! are 0 since rain water mixing ratio does not vary in this component at
       ! this grid level.
       corr_srr_1 = zero
    endif

    ! Correlation (in-precip) between s and r_r in PDF component 2.
    if ( rr2 > rr_tol ) then
       if ( l_interp_prescribed_params ) then
          corr_srr_2 = cloud_frac2 * corr_srr_NL_cloud &
                       + ( one - cloud_frac2 ) * corr_srr_NL_below
       else
          if ( rc2 > rc_tol ) then
             corr_srr_2 = corr_srr_NL_cloud
          else
             corr_srr_2 = corr_srr_NL_below
          endif
       endif
    else
       ! Mean in-precip rain water mixing ratio in PDF component 2 is less than
       ! the tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 2nd PDF component at this grid level.  The
       ! correlations involving rain water mixing ratio in the 2nd PDF component
       ! are 0 since rain water mixing ratio does not vary in this component at
       ! this grid level.
       corr_srr_2 = zero
    endif

    ! Correlation (in-precip) between s and N_r in PDF component 1.
    if ( Nr1 > Nr_tol ) then
       if ( l_interp_prescribed_params ) then
          corr_sNr_1 = cloud_frac1 * corr_sNr_NL_cloud &
                       + ( one - cloud_frac1 ) * corr_sNr_NL_below
       else
          if ( rc1 > rc_tol ) then
             corr_sNr_1 = corr_sNr_NL_cloud
          else
             corr_sNr_1 = corr_sNr_NL_below
          endif
       endif
    else
       ! Mean in-precip rain drop concentration in PDF component 1 is less than
       ! the tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 1st PDF component at this grid level.  The
       ! correlations involving rain drop concentration in the 1st PDF component
       ! are 0 since rain drop concentration does not vary in this component at
       ! this grid level.
       corr_sNr_1 = zero
    endif

    ! Correlation (in-precip) between s and N_r in PDF component 2.
    if ( Nr2 > Nr_tol ) then
       if ( l_interp_prescribed_params ) then
          corr_sNr_2 = cloud_frac2 * corr_sNr_NL_cloud &
                       + ( one - cloud_frac2 ) * corr_sNr_NL_below
       else
          if ( rc2 > rc_tol ) then
             corr_sNr_2 = corr_sNr_NL_cloud
          else
             corr_sNr_2 = corr_sNr_NL_below
          endif
       endif
    else
       ! Mean in-precip rain drop concentration in PDF component 2 is less than
       ! the tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 2nd PDF component at this grid level.  The
       ! correlations involving rain drop concentration in the 2nd PDF component
       ! are 0 since rain drop concentration does not vary in this component at
       ! this grid level.
       corr_sNr_2 = zero
    endif

    ! Correlation between s and N_cn in PDF component 1.
    if ( Ncnm > Ncn_tol ) then
       corr_sNcn_1 = corr_sNcn_NL_cloud
    else
       ! Mean cloud nuclei concentration is less than the tolerance amount.  It
       ! is considered to have a value of 0.  There are not any cloud nuclei or
       ! cloud at this grid level.  The correlations involving cloud nuclei
       ! concentration are 0 since cloud nuclei concentration does not vary at
       ! this grid level.
       corr_sNcn_1 = zero
    endif

    ! Correlation between s and N_cn in PDF component 2.
    if ( Ncnm > Ncn_tol ) then
       corr_sNcn_2 = corr_sNcn_NL_cloud
    else
       ! Mean cloud nuclei concentration is less than the tolerance amount.  It
       ! is considered to have a value of 0.  There are not any cloud nuclei or
       ! cloud at this grid level.  The correlations involving cloud nuclei
       ! concentration are 0 since cloud nuclei concentration does not vary at
       ! this grid level.
       corr_sNcn_2 = zero
    endif

    ! Correlation (in-precip) between t and r_r in PDF component 1.
    if ( rr1 > rr_tol ) then
       if ( l_interp_prescribed_params ) then
          corr_trr_1 = cloud_frac1 * corr_trr_NL_cloud &
                       + ( one - cloud_frac1 ) * corr_trr_NL_below
       else
          if ( rc1 > rc_tol ) then
             corr_trr_1 = corr_trr_NL_cloud
          else
             corr_trr_1 = corr_trr_NL_below
          endif
       endif
    else
       ! Mean in-precip rain water mixing ratio in PDF component 1 is less than
       ! the tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 1st PDF component at this grid level.  The
       ! correlations involving rain water mixing ratio in the 1st PDF component
       ! are 0 since rain water mixing ratio does not vary in this component at
       ! this grid level.
       corr_trr_1 = zero
    endif

    ! Correlation (in-precip) between t and r_r in PDF component 2.
    if ( rr2 > rr_tol ) then
       if ( l_interp_prescribed_params ) then
          corr_trr_2 = cloud_frac2 * corr_trr_NL_cloud &
                       + ( one - cloud_frac2 ) * corr_trr_NL_below
       else
          if ( rc2 > rc_tol ) then
             corr_trr_2 = corr_trr_NL_cloud
          else
             corr_trr_2 = corr_trr_NL_below
          endif
       endif
    else
       ! Mean in-precip rain water mixing ratio in PDF component 2 is less than
       ! the tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 2nd PDF component at this grid level.  The
       ! correlations involving rain water mixing ratio in the 2nd PDF component
       ! are 0 since rain water mixing ratio does not vary in this component at
       ! this grid level.
       corr_trr_2 = zero
    endif

    ! Correlation (in-precip) between t and N_r in PDF component 1.
    if ( Nr1 > Nr_tol ) then
       if ( l_interp_prescribed_params ) then
          corr_tNr_1 = cloud_frac1 * corr_tNr_NL_cloud &
                       + ( one - cloud_frac1 ) * corr_tNr_NL_below
       else
          if ( rc1 > rc_tol ) then
             corr_tNr_1 = corr_tNr_NL_cloud
          else
             corr_tNr_1 = corr_tNr_NL_below
          endif
       endif
    else
       ! Mean in-precip rain drop concentration in PDF component 1 is less than
       ! the tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 1st PDF component at this grid level.  The
       ! correlations involving rain drop concentration in the 1st PDF component
       ! are 0 since rain drop concentration does not vary in this component at
       ! this grid level.
       corr_tNr_1 = zero
    endif

    ! Correlation (in-precip) between t and N_r in PDF component 2.
    if ( Nr2 > Nr_tol ) then
       if ( l_interp_prescribed_params ) then
          corr_tNr_2 = cloud_frac2 * corr_tNr_NL_cloud &
                       + ( one - cloud_frac2 ) * corr_tNr_NL_below
       else
          if ( rc2 > rc_tol ) then
             corr_tNr_2 = corr_tNr_NL_cloud
          else
             corr_tNr_2 = corr_tNr_NL_below
          endif
       endif
    else
       ! Mean in-precip rain drop concentration in PDF component 2 is less than
       ! the tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 2nd PDF component at this grid level.  The
       ! correlations involving rain drop concentration in the 2nd PDF component
       ! are 0 since rain drop concentration does not vary in this component at
       ! this grid level.
       corr_tNr_2 = zero
    endif

    ! Correlation between t and N_cn in PDF component 1.
    if ( Ncnm > Ncn_tol ) then
       corr_tNcn_1 = corr_tNcn_NL_cloud
    else
       ! Mean cloud nuclei concentration is less than the tolerance amount.  It
       ! is considered to have a value of 0.  There are not any cloud nuclei or
       ! cloud at this grid level.  The correlations involving cloud nuclei
       ! concentration are 0 since cloud nuclei concentration does not vary at
       ! this grid level.
       corr_tNcn_1 = zero
    endif

    ! Correlation between t and N_cn in PDF component 2.
    if ( Ncnm > Ncn_tol ) then
       corr_tNcn_2 = corr_tNcn_NL_cloud
    else
       ! Mean cloud nuclei concentration is less than the tolerance amount.  It
       ! is considered to have a value of 0.  There are not any cloud nuclei or
       ! cloud at this grid level.  The correlations involving cloud nuclei
       ! concentration are 0 since cloud nuclei concentration does not vary at
       ! this grid level.
       corr_tNcn_2 = zero
    endif

    ! Correlation (in-precip) between r_r and N_r in PDF component 1.
    if ( rr1 > rr_tol .and. Nr1 > Nr_tol ) then
       if ( l_interp_prescribed_params ) then
          corr_rrNr_1 = cloud_frac1 * corr_rrNr_LL_cloud &
                        + ( one - cloud_frac1 ) * corr_rrNr_LL_below
       else
          if ( rc1 > rc_tol ) then
             corr_rrNr_1 = corr_rrNr_LL_cloud
          else
             corr_rrNr_1 = corr_rrNr_LL_below
          endif
       endif
    else
       ! Mean in-precip rain water mixing ratio in PDF component 1 and (or) mean
       ! in-precip rain drop concentration in PDF component 1 are (is) less than
       ! their (its) respective tolerance amount(s), and are (is) considered to
       ! have a value of 0.  There is not any rain in the 1st PDF component at
       ! this grid level.  The correlation is 0 since rain does not vary in this
       ! component at this grid level.
       corr_rrNr_1 = zero
    endif

    ! Correlation (in-precip) between r_r and N_r in PDF component 2.
    if ( rr2 > rr_tol .and. Nr2 > Nr_tol ) then
       if ( l_interp_prescribed_params ) then
          corr_rrNr_2 = cloud_frac2 * corr_rrNr_LL_cloud &
                        + ( one - cloud_frac2 ) * corr_rrNr_LL_below
       else
          if ( rc2 > rc_tol ) then
             corr_rrNr_2 = corr_rrNr_LL_cloud
          else
             corr_rrNr_2 = corr_rrNr_LL_below
          endif
       endif
    else
       ! Mean in-precip rain water mixing ratio in PDF component 2 and (or) mean
       ! in-precip rain drop concentration in PDF component 2 are (is) less than
       ! their (its) respective tolerance amount(s), and are (is) considered to
       ! have a value of 0.  There is not any rain in the 2nd PDF component at
       ! this grid level.  The correlation is 0 since rain does not vary in this
       ! component at this grid level.
       corr_rrNr_2 = zero
    endif

    return

  end subroutine compute_corr

  !=============================================================================
  subroutine normalize_pdf_params( rr1, rr2, Nr1, Nr2, Ncnm, &                   ! Intent(in)
                                   mu_rr_1, mu_rr_2, mu_Nr_1, mu_Nr_2, &         ! Intent(in)
                                   mu_Ncn_1, mu_Ncn_2, sigma_rr_1, sigma_rr_2, & ! Intent(in)
                                   sigma_Nr_1, sigma_Nr_2, sigma_Ncn_1, &        ! Intent(in)
                                   sigma_Ncn_2, corr_wrr_1, corr_wrr_2, &        ! Intent(in)
                                   corr_wNr_1, corr_wNr_2, corr_wNcn_1, &        ! Intent(in)
                                   corr_wNcn_2, corr_srr_1, corr_srr_2, &        ! Intent(in)
                                   corr_sNr_1, corr_sNr_2, corr_sNcn_1, &        ! Intent(in)
                                   corr_sNcn_2, corr_trr_1, corr_trr_2, &        ! Intent(in)
                                   corr_tNr_1, corr_tNr_2, corr_tNcn_1, &        ! Intent(in)
                                   corr_tNcn_2, corr_rrNr_1, corr_rrNr_2, &      ! Intent(in)
                                   mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, &            ! Intent(out)
                                   mu_Nr_2_n, mu_Ncn_1_n, mu_Ncn_2_n, &          ! Intent(out)
                                   sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, &   ! Intent(out)
                                   sigma_Nr_2_n, sigma_Ncn_1_n, sigma_Ncn_2_n, & ! Intent(out)
                                   corr_wrr_1_n, corr_wrr_2_n, corr_wNr_1_n, &   ! Intent(out)
                                   corr_wNr_2_n, corr_wNcn_1_n, corr_wNcn_2_n, & ! Intent(out)
                                   corr_srr_1_n, corr_srr_2_n, corr_sNr_1_n, &   ! Intent(out)
                                   corr_sNr_2_n, corr_sNcn_1_n, corr_sNcn_2_n, & ! Intent(out)
                                   corr_trr_1_n, corr_trr_2_n, corr_tNr_1_n, &   ! Intent(out)
                                   corr_tNr_2_n, corr_tNcn_1_n, corr_tNcn_2_n, & ! Intent(out)
                                   corr_rrNr_1_n, corr_rrNr_2_n )                ! Intent(out)

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        rr_tol, & ! Constant(s)
        Nr_tol, &
        Ncn_tol, &
        zero

    use KK_utilities, only: &
        mean_L2N,   & ! Procedure(s)
        stdev_L2N,  &
        corr_NL2NN, &
        corr_LL2NN

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      rr1,         & ! Mean rain water mixing ratio (1st PDF component)  [kg/kg]
      rr2,         & ! Mean rain water mixing ratio (2nd PDF component)  [kg/kg]
      Nr1,         & ! Mean rain drop concentration (1st PDF component) [num/kg]
      Nr2,         & ! Mean rain drop concentration (2nd PDF component) [num/kg]
      Ncnm,        & ! Mean cloud nuclei concentration (overall)        [num/kg]
      mu_rr_1,     & ! Mean of rr (1st PDF component) in-precip (ip)     [kg/kg]
      mu_rr_2,     & ! Mean of rr (2nd PDF component) ip                 [kg/kg]
      mu_Nr_1,     & ! Mean of Nr (1st PDF component) ip                [num/kg]
      mu_Nr_2,     & ! Mean of Nr (2nd PDF component) ip                [num/kg]
      mu_Ncn_1,    & ! Mean of Ncn (1st PDF component)                  [num/kg]
      mu_Ncn_2,    & ! Mean of Ncn (2nd PDF component)                  [num/kg]
      sigma_rr_1,  & ! Standard deviation of rr (1st PDF component) ip   [kg/kg]
      sigma_rr_2,  & ! Standard deviation of rr (2nd PDF component) ip   [kg/kg]
      sigma_Nr_1,  & ! Standard deviation of Nr (1st PDF component) ip  [num/kg]
      sigma_Nr_2,  & ! Standard deviation of Nr (2nd PDF component) ip  [num/kg]
      sigma_Ncn_1, & ! Standard deviation of Ncn (1st PDF component)    [num/kg]
      sigma_Ncn_2, & ! Standard deviation of Ncn (2nd PDF component)    [num/kg]
      corr_wrr_1,  & ! Correlation between w and rr (1st PDF component) ip   [-]
      corr_wrr_2,  & ! Correlation between w and rr (2nd PDF component) ip   [-]
      corr_wNr_1,  & ! Correlation between w and Nr (1st PDF component) ip   [-]
      corr_wNr_2,  & ! Correlation between w and Nr (2nd PDF component) ip   [-]
      corr_wNcn_1, & ! Correlation between w and Ncn (1st PDF component)     [-]
      corr_wNcn_2, & ! Correlation between w and Ncn (2nd PDF component)     [-]
      corr_srr_1,  & ! Correlation between s and rr (1st PDF component) ip   [-]
      corr_srr_2,  & ! Correlation between s and rr (2nd PDF component) ip   [-]
      corr_sNr_1,  & ! Correlation between s and Nr (1st PDF component) ip   [-]
      corr_sNr_2,  & ! Correlation between s and Nr (2nd PDF component) ip   [-]
      corr_sNcn_1, & ! Correlation between s and Ncn (1st PDF component)     [-]
      corr_sNcn_2, & ! Correlation between s and Ncn (2nd PDF component)     [-]
      corr_trr_1,  & ! Correlation between t and rr (1st PDF component) ip   [-]
      corr_trr_2,  & ! Correlation between t and rr (2nd PDF component) ip   [-]
      corr_tNr_1,  & ! Correlation between t and Nr (1st PDF component) ip   [-]
      corr_tNr_2,  & ! Correlation between t and Nr (2nd PDF component) ip   [-]
      corr_tNcn_1, & ! Correlation between t and Ncn (1st PDF component)     [-]
      corr_tNcn_2, & ! Correlation between t and Ncn (2nd PDF component)     [-]
      corr_rrNr_1, & ! Correlation between rr & Nr (1st PDF component) ip    [-]
      corr_rrNr_2    ! Correlation between rr & Nr (2nd PDF component) ip    [-]

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
      mu_rr_1_n,     & ! Mean of ln rr (1st PDF component) ip        [ln(kg/kg)]
      mu_rr_2_n,     & ! Mean of ln rr (2nd PDF component) ip        [ln(kg/kg)]
      mu_Nr_1_n,     & ! Mean of ln Nr (1st PDF component) ip       [ln(num/kg)]
      mu_Nr_2_n,     & ! Mean of ln Nr (2nd PDF component) ip       [ln(num/kg)]
      mu_Ncn_1_n,    & ! Mean of ln Ncn (1st PDF component)         [ln(num/kg)]
      mu_Ncn_2_n,    & ! Mean of ln Ncn (2nd PDF component)         [ln(num/kg)]
      sigma_rr_1_n,  & ! Standard dev. of ln rr (1st PDF comp.) ip   [ln(kg/kg)]
      sigma_rr_2_n,  & ! Standard dev. of ln rr (2nd PDF comp.) ip   [ln(kg/kg)]
      sigma_Nr_1_n,  & ! Standard dev. of ln Nr (1st PDF comp.) ip  [ln(num/kg)]
      sigma_Nr_2_n,  & ! Standard dev. of ln Nr (2nd PDF comp.) ip  [ln(num/kg)]
      sigma_Ncn_1_n, & ! Standard dev. of ln Ncn (1st PDF comp.)    [ln(num/kg)]
      sigma_Ncn_2_n, & ! Standard dev. of ln Ncn (2nd PDF comp.)    [ln(num/kg)]
      corr_wrr_1_n,  & ! Correlation between w and ln rr (1st PDF comp.) ip  [-]
      corr_wrr_2_n,  & ! Correlation between w and ln rr (2nd PDF comp.) ip  [-]
      corr_wNr_1_n,  & ! Correlation between w and ln Nr (1st PDF comp.) ip  [-]
      corr_wNr_2_n,  & ! Correlation between w and ln Nr (2nd PDF comp.) ip  [-]
      corr_wNcn_1_n, & ! Correlation between w and ln Ncn (1st PDF comp.)    [-]
      corr_wNcn_2_n, & ! Correlation between w and ln Ncn (2nd PDF comp.)    [-]
      corr_srr_1_n,  & ! Correlation between s and ln rr (1st PDF comp.) ip  [-]
      corr_srr_2_n,  & ! Correlation between s and ln rr (2nd PDF comp.) ip  [-]
      corr_sNr_1_n,  & ! Correlation between s and ln Nr (1st PDF comp.) ip  [-]
      corr_sNr_2_n,  & ! Correlation between s and ln Nr (2nd PDF comp.) ip  [-]
      corr_sNcn_1_n, & ! Correlation between s and ln Ncn (1st PDF comp.)    [-]
      corr_sNcn_2_n, & ! Correlation between s and ln Ncn (2nd PDF comp.)    [-]
      corr_trr_1_n,  & ! Correlation between t and ln rr (1st PDF comp.) ip  [-]
      corr_trr_2_n,  & ! Correlation between t and ln rr (2nd PDF comp.) ip  [-]
      corr_tNr_1_n,  & ! Correlation between t and ln Nr (1st PDF comp.) ip  [-]
      corr_tNr_2_n,  & ! Correlation between t and ln Nr (2nd PDF comp.) ip  [-]
      corr_tNcn_1_n, & ! Correlation between t and ln Ncn (1st PDF comp.)    [-]
      corr_tNcn_2_n, & ! Correlation between t and ln Ncn (2nd PDF comp.)    [-]
      corr_rrNr_1_n, & ! Correlation btwn. ln rr & ln Nr (1st PDF comp.) ip  [-]
      corr_rrNr_2_n    ! Correlation btwn. ln rr & ln Nr (2nd PDF comp.) ip  [-]


    ! --- Begin Code ---

    !!! Calculate the normalized mean of variables that have an assumed
    !!! lognormal distribution, given the mean and variance of those variables.

    ! Normalized mean of in-precip rain water mixing ratio in PDF component 1.
    if ( rr1 > rr_tol ) then
       mu_rr_1_n = mean_L2N( mu_rr_1, sigma_rr_1**2 )
    else
       ! Mean rain water mixing ratio in PDF component 1 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 1st PDF component at this grid level.  The mean
       ! in-precip rain water mixing ratio (1st PDF component) is also 0.  The
       ! value of mu_rr_1_n should be -inf.  It will be set to -huge for
       ! purposes of assigning it a value.
       mu_rr_1_n = -huge( mu_rr_1_n )
    endif

    ! Normalized mean of in-precip rain water mixing ratio in PDF component 2.
    if ( rr2 > rr_tol ) then
       mu_rr_2_n = mean_L2N( mu_rr_2, sigma_rr_2**2 )
    else
       ! Mean rain water mixing ratio in PDF component 2 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 2nd PDF component at this grid level.  The mean
       ! in-precip rain water mixing ratio (2nd PDF component) is also 0.  The
       ! value of mu_rr_2_n should be -inf.  It will be set to -huge for
       ! purposes of assigning it a value.
       mu_rr_2_n = -huge( mu_rr_2_n )
    endif

    ! Normalized mean of in-precip rain drop concentration in PDF component 1.
    if ( Nr1 > Nr_tol ) then
       mu_Nr_1_n = mean_L2N( mu_Nr_1, sigma_Nr_1**2 )
    else
       ! Mean rain drop concentration in PDF component 1 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 1st PDF component at this grid level.  The mean
       ! in-precip rain drop concentration (1st PDF component) is also 0.  The
       ! value of mu_Nr_1_n should be -inf.  It will be set to -huge for
       ! purposes of assigning it a value.
       mu_Nr_1_n = -huge( mu_Nr_1_n )
    endif

    ! Normalized mean of in-precip rain drop concentration in PDF component 2.
    if ( Nr2 > Nr_tol ) then
       mu_Nr_2_n = mean_L2N( mu_Nr_2, sigma_Nr_2**2 )
    else
       ! Mean rain drop concentration in PDF component 2 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 2nd PDF component at this grid level.  The mean
       ! in-precip rain drop concentration (2nd PDF component) is also 0.  The
       ! value of mu_Nr_2_n should be -inf.  It will be set to -huge for
       ! purposes of assigning it a value.
       mu_Nr_2_n = -huge( mu_Nr_2_n )
    endif

    ! Normalized mean of cloud nuclei concentration in PDF component 1.
    if ( Ncnm > Ncn_tol ) then
       mu_Ncn_1_n = mean_L2N( mu_Ncn_1, sigma_Ncn_1**2 )
    else
       ! Mean cloud nuclei concentration is less than the tolerance amount.  It
       ! is considered to have a value of 0.  There are not any cloud nuclei or
       ! cloud at this grid level.  The value of mu_Ncn_1_n should be -inf.  It
       ! will be set to -huge for purposes of assigning it a value.
       mu_Ncn_1_n = -huge( mu_Ncn_1_n )
    endif

    ! Normalized mean of cloud nuclei concentration in PDF component 2.
    if ( Ncnm > Ncn_tol ) then
       mu_Ncn_2_n = mean_L2N( mu_Ncn_2, sigma_Ncn_2**2 )
    else
       ! Mean cloud nuclei concentration is less than the tolerance amount.  It
       ! is considered to have a value of 0.  There are not any cloud nuclei or
       ! cloud at this grid level.  The value of mu_Ncn_2_n should be -inf.  It
       ! will be set to -huge for purposes of assigning it a value.
       mu_Ncn_2_n = -huge( mu_Ncn_2_n )
    endif

    !!! Calculate the normalized standard deviation of variables that have
    !!! an assumed lognormal distribution, given the mean and variance of
    !!! those variables.

    ! Normalized standard deviation of in-precip rain water mixing ratio
    ! in PDF component 1.
    if ( rr1 > rr_tol ) then
       sigma_rr_1_n = stdev_L2N( mu_rr_1, sigma_rr_1**2 )
    else
       ! Mean rain water mixing ratio in PDF component 1 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 1st PDF component at this grid level.  The mean
       ! in-precip rain water mixing ratio (1st PDF component) is also 0.  The
       ! standard deviation is simply 0 since rain water mixing ratio does not
       ! vary in this component at this grid level.
       sigma_rr_1_n = zero
    endif

    ! Normalized standard deviation of in-precip rain water mixing ratio
    ! in PDF component 2.
    if ( rr2 > rr_tol ) then
       sigma_rr_2_n = stdev_L2N( mu_rr_2, sigma_rr_2**2 )
    else
       ! Mean rain water mixing ratio in PDF component 2 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 2nd PDF component at this grid level.  The mean
       ! in-precip rain water mixing ratio (2nd PDF component) is also 0.  The
       ! standard deviation is simply 0 since rain water mixing ratio does not
       ! vary in this component at this grid level.
       sigma_rr_2_n = zero
    endif

    ! Normalized standard deviation of in-precip rain drop concentration
    ! in PDF component 1.
    if ( Nr1 > Nr_tol ) then
       sigma_Nr_1_n = stdev_L2N( mu_Nr_1, sigma_Nr_1**2 )
    else
       ! Mean rain drop concentration in PDF component 1 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 1st PDF component at this grid level.  The mean
       ! in-precip rain drop concentration (1st PDF component) is also 0.  The
       ! standard deviation is simply 0 since rain water mixing ratio does not
       ! vary in this component at this grid level.
       sigma_Nr_1_n = zero
    endif

    ! Normalized standard deviation of in-precip rain drop concentration
    ! in PDF component 2.
    if ( Nr2 > Nr_tol ) then
       sigma_Nr_2_n = stdev_L2N( mu_Nr_2, sigma_Nr_2**2 )
    else
       ! Mean rain drop concentration in PDF component 2 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 2nd PDF component at this grid level.  The mean
       ! in-precip rain drop concentration (2nd PDF component) is also 0.  The
       ! standard deviation is simply 0 since rain water mixing ratio does not
       ! vary in this component at this grid level.
       sigma_Nr_2_n = zero
    endif

    ! Normalized standard deviation of cloud nuclei concentration
    ! in PDF component 1.
    if ( Ncnm > Ncn_tol ) then
       sigma_Ncn_1_n = stdev_L2N( mu_Ncn_1, sigma_Ncn_1**2 )
    else
       ! Mean cloud nuclei concentration is less than the tolerance amount.  It
       ! is considered to have a value of 0.  There are not any cloud nuclei or
       ! cloud at this grid level.  The standard deviation is simply 0 since
       ! cloud nuclei concentration does not vary at this grid level.
       sigma_Ncn_1_n = zero
    endif

    ! Normalized standard deviation of cloud nuclei concentration
    ! in PDF component 2.
    if ( Ncnm > Ncn_tol ) then
       sigma_Ncn_2_n = stdev_L2N( mu_Ncn_2, sigma_Ncn_2**2 )
    else
       ! Mean cloud nuclei concentration is less than the tolerance amount.  It
       ! is considered to have a value of 0.  There are not any cloud nuclei or
       ! cloud at this grid level.  The standard deviation is simply 0 since
       ! cloud nuclei concentration does not vary at this grid level.
       sigma_Ncn_2_n = zero
    endif

    !!! Calculate the normalized correlation between variables that have
    !!! an assumed normal distribution and variables that have an assumed
    !!! lognormal distribution for the ith PDF component, given their
    !!! correlation and the normalized standard deviation of the variable with
    !!! the assumed lognormal distribution.

    ! Normalize the correlation (in-precip) between w and r_r
    ! in PDF component 1.
    if ( rr1 > rr_tol ) then
       corr_wrr_1_n = corr_NL2NN( corr_wrr_1, sigma_rr_1_n )
    else
       ! Mean rain water mixing ratio in PDF component 1 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 1st PDF component at this grid level.  The mean
       ! in-precip rain water mixing ratio (1st PDF component) is also 0.  The
       ! correlations involving in-precip rain water mixing ratio (1st PDF
       ! component) are 0 since in-precip rain water mixing ratio does not vary
       ! in this component at this grid level.
       corr_wrr_1_n = zero
    endif

    ! Normalize the correlation (in-precip) between w and r_r
    ! in PDF component 2.
    if ( rr2 > rr_tol ) then
       corr_wrr_2_n = corr_NL2NN( corr_wrr_2, sigma_rr_2_n )
    else
       ! Mean rain water mixing ratio in PDF component 2 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 2nd PDF component at this grid level.  The mean
       ! in-precip rain water mixing ratio (2nd PDF component) is also 0.  The
       ! correlations involving in-precip rain water mixing ratio (2nd PDF
       ! component) are 0 since in-precip rain water mixing ratio does not vary
       ! in this component at this grid level.
       corr_wrr_2_n = zero
    endif

    ! Normalize the correlation (in-precip) between w and N_r
    ! in PDF component 1.
    if ( Nr1 > Nr_tol ) then
       corr_wNr_1_n = corr_NL2NN( corr_wNr_1, sigma_Nr_1_n )
    else
       ! Mean rain drop concentration in PDF component 1 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 1st PDF component at this grid level.  The mean
       ! in-precip rain drop concentration (1st PDF component) is also 0.  The
       ! correlations involving in-precip rain drop concentration (1st PDF
       ! component) are 0 since in-precip rain drop concentration does not vary
       ! in this component at this grid level.
       corr_wNr_1_n = zero
    endif

    ! Normalize the correlation (in-precip) between w and N_r
    ! in PDF component 2.
    if ( Nr2 > Nr_tol ) then
       corr_wNr_2_n = corr_NL2NN( corr_wNr_2, sigma_Nr_2_n )
    else
       ! Mean rain drop concentration in PDF component 2 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 2nd PDF component at this grid level.  The mean
       ! in-precip rain drop concentration (2nd PDF component) is also 0.  The
       ! correlations involving in-precip rain drop concentration (2nd PDF
       ! component) are 0 since in-precip rain drop concentration does not vary
       ! in this component at this grid level.
       corr_wNr_2_n = zero
    endif

    ! Normalize the correlation between w and N_cn in PDF component 1.
    if ( Ncnm > Ncn_tol ) then
       corr_wNcn_1_n = corr_NL2NN( corr_wNcn_1, sigma_Ncn_1_n )
    else
       ! Mean cloud nuclei concentration is less than the tolerance amount.  It
       ! is considered to have a value of 0.  There are not any cloud nuclei or
       ! cloud at this grid level.  The correlations involving cloud nuclei
       ! concentration are 0 since cloud nuclei concentration does not vary at
       ! this grid level.
       corr_wNcn_1_n = zero
    endif

    ! Normalize the correlation between w and N_cn in PDF component 2.
    if ( Ncnm > Ncn_tol ) then
       corr_wNcn_2_n = corr_NL2NN( corr_wNcn_2, sigma_Ncn_2_n )
    else
       ! Mean cloud nuclei concentration is less than the tolerance amount.  It
       ! is considered to have a value of 0.  There are not any cloud nuclei or
       ! cloud at this grid level.  The correlations involving cloud nuclei
       ! concentration are 0 since cloud nuclei concentration does not vary at
       ! this grid level.
       corr_wNcn_2_n = zero
    endif

    ! Normalize the correlation (in-precip) between s and r_r
    ! in PDF component 1.
    if ( rr1 > rr_tol ) then
       corr_srr_1_n = corr_NL2NN( corr_srr_1, sigma_rr_1_n )
    else
       ! Mean rain water mixing ratio in PDF component 1 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 1st PDF component at this grid level.  The mean
       ! in-precip rain water mixing ratio (1st PDF component) is also 0.  The
       ! correlations involving in-precip rain water mixing ratio (1st PDF
       ! component) are 0 since in-precip rain water mixing ratio does not vary
       ! in this component at this grid level.
       corr_srr_1_n = zero
    endif

    ! Normalize the correlation (in-precip) between s and r_r
    ! in PDF component 2.
    if ( rr2 > rr_tol ) then
       corr_srr_2_n = corr_NL2NN( corr_srr_2, sigma_rr_2_n )
    else
       ! Mean rain water mixing ratio in PDF component 2 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 2nd PDF component at this grid level.  The mean
       ! in-precip rain water mixing ratio (2nd PDF component) is also 0.  The
       ! correlations involving in-precip rain water mixing ratio (2nd PDF
       ! component) are 0 since in-precip rain water mixing ratio does not vary
       ! in this component at this grid level.
       corr_srr_2_n = zero
    endif

    ! Normalize the correlation (in-precip) between s and N_r
    ! in PDF component 1.
    if ( Nr1 > Nr_tol ) then
       corr_sNr_1_n = corr_NL2NN( corr_sNr_1, sigma_Nr_1_n )
    else
       ! Mean rain drop concentration in PDF component 1 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 1st PDF component at this grid level.  The mean
       ! in-precip rain drop concentration (1st PDF component) is also 0.  The
       ! correlations involving in-precip rain drop concentration (1st PDF
       ! component) are 0 since in-precip rain drop concentration does not vary
       ! in this component at this grid level.
       corr_sNr_1_n = zero
    endif

    ! Normalize the correlation (in-precip) between s and N_r
    ! in PDF component 2.
    if ( Nr2 > Nr_tol ) then
       corr_sNr_2_n = corr_NL2NN( corr_sNr_2, sigma_Nr_2_n )
    else
       ! Mean rain drop concentration in PDF component 2 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 2nd PDF component at this grid level.  The mean
       ! in-precip rain drop concentration (2nd PDF component) is also 0.  The
       ! correlations involving in-precip rain drop concentration (2nd PDF
       ! component) are 0 since in-precip rain drop concentration does not vary
       ! in this component at this grid level.
       corr_sNr_2_n = zero
    endif

    ! Normalize the correlation between s and N_cn in PDF component 1.
    if ( Ncnm > Ncn_tol ) then
       corr_sNcn_1_n = corr_NL2NN( corr_sNcn_1, sigma_Ncn_1_n )
    else
       ! Mean cloud nuclei concentration is less than the tolerance amount.  It
       ! is considered to have a value of 0.  There are not any cloud nuclei or
       ! cloud at this grid level.  The correlations involving cloud nuclei
       ! concentration are 0 since cloud nuclei concentration does not vary at
       ! this grid level.
       corr_sNcn_1_n = zero
    endif

    ! Normalize the correlation between s and N_cn in PDF component 2.
    if ( Ncnm > Ncn_tol ) then
       corr_sNcn_2_n = corr_NL2NN( corr_sNcn_2, sigma_Ncn_2_n )
    else
       ! Mean cloud nuclei concentration is less than the tolerance amount.  It
       ! is considered to have a value of 0.  There are not any cloud nuclei or
       ! cloud at this grid level.  The correlations involving cloud nuclei
       ! concentration are 0 since cloud nuclei concentration does not vary at
       ! this grid level.
       corr_sNcn_2_n = zero
    endif

    ! Normalize the correlation (in-precip) between t and r_r
    ! in PDF component 1.
    if ( rr1 > rr_tol ) then
       corr_trr_1_n = corr_NL2NN( corr_trr_1, sigma_rr_1_n )
    else
       ! Mean rain water mixing ratio in PDF component 1 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 1st PDF component at this grid level.  The mean
       ! in-precip rain water mixing ratio (1st PDF component) is also 0.  The
       ! correlations involving in-precip rain water mixing ratio (1st PDF
       ! component) are 0 since in-precip rain water mixing ratio does not vary
       ! in this component at this grid level.
       corr_trr_1_n = zero
    endif

    ! Normalize the correlation (in-precip) between t and r_r
    ! in PDF component 2.
    if ( rr2 > rr_tol ) then
       corr_trr_2_n = corr_NL2NN( corr_trr_2, sigma_rr_2_n )
    else
       ! Mean rain water mixing ratio in PDF component 2 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 2nd PDF component at this grid level.  The mean
       ! in-precip rain water mixing ratio (2nd PDF component) is also 0.  The
       ! correlations involving in-precip rain water mixing ratio (2nd PDF
       ! component) are 0 since in-precip rain water mixing ratio does not vary
       ! in this component at this grid level.
       corr_trr_2_n = zero
    endif

    ! Normalize the correlation (in-precip) between t and N_r
    ! in PDF component 1.
    if ( Nr1 > Nr_tol ) then
       corr_tNr_1_n = corr_NL2NN( corr_tNr_1, sigma_Nr_1_n )
    else
       ! Mean rain drop concentration in PDF component 1 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 1st PDF component at this grid level.  The mean
       ! in-precip rain drop concentration (1st PDF component) is also 0.  The
       ! correlations involving in-precip rain drop concentration (1st PDF
       ! component) are 0 since in-precip rain drop concentration does not vary
       ! in this component at this grid level.
       corr_tNr_1_n = zero
    endif

    ! Normalize the correlation (in-precip) between t and N_r
    ! in PDF component 2.
    if ( Nr2 > Nr_tol ) then
       corr_tNr_2_n = corr_NL2NN( corr_tNr_2, sigma_Nr_2_n )
    else
       ! Mean rain drop concentration in PDF component 2 is less than the
       ! tolerance amount.  It is considered to have a value of 0.  There is
       ! not any rain in the 2nd PDF component at this grid level.  The mean
       ! in-precip rain drop concentration (2nd PDF component) is also 0.  The
       ! correlations involving in-precip rain drop concentration (2nd PDF
       ! component) are 0 since in-precip rain drop concentration does not vary
       ! in this component at this grid level.
       corr_tNr_2_n = zero
    endif

    ! Normalize the correlation between t and N_cn in PDF component 1.
    if ( Ncnm > Ncn_tol ) then
       corr_tNcn_1_n = corr_NL2NN( corr_tNcn_1, sigma_Ncn_1_n )
    else
       ! Mean cloud nuclei concentration is less than the tolerance amount.  It
       ! is considered to have a value of 0.  There are not any cloud nuclei or
       ! cloud at this grid level.  The correlations involving cloud nuclei
       ! concentration are 0 since cloud nuclei concentration does not vary at
       ! this grid level.
       corr_tNcn_1_n = zero
    endif

    ! Normalize the correlation between t and N_cn in PDF component 2.
    if ( Ncnm > Ncn_tol ) then
       corr_tNcn_2_n = corr_NL2NN( corr_tNcn_2, sigma_Ncn_2_n )
    else
       ! Mean cloud nuclei concentration is less than the tolerance amount.  It
       ! is considered to have a value of 0.  There are not any cloud nuclei or
       ! cloud at this grid level.  The correlations involving cloud nuclei
       ! concentration are 0 since cloud nuclei concentration does not vary at
       ! this grid level.
       corr_tNcn_2_n = zero
    endif

    !!! Calculate the normalized correlation between two variables that both
    !!! have an assumed lognormal distribution, given their correlation and both
    !!! of their normalized standard deviations.

    ! Normalize the correlation (in-precip) between r_r and N_r
    ! in PDF component 1.
    if ( rr1 > rr_tol .and. Nr1 > Nr_tol ) then
       corr_rrNr_1_n = corr_LL2NN( corr_rrNr_1, sigma_rr_1_n, sigma_Nr_1_n )
    else
       ! Mean rain water mixing ratio in PDF component 1 and (or) mean rain drop
       ! concentration in PDF component 1 are (is) less than their (its)
       ! respective tolerance amount(s), and are (is) considered to have a value
       ! of 0.  There is not any rain at this grid level.  The mean in-precip
       ! rain water mixing ratio (1st PDF component) and (or) mean in-precip
       ! rain drop concentration (1st PDF component) are also considered to have
       ! a value of 0.  The correlation is 0 since rain does not vary in this
       ! component at this grid level.
       corr_rrNr_1_n = zero
    endif

    ! Normalize the correlation (in-precip) between r_r and N_r
    ! in PDF component 2.
    if ( rr2 > rr_tol .and. Nr2 > Nr_tol ) then
       corr_rrNr_2_n = corr_LL2NN( corr_rrNr_2, sigma_rr_2_n, sigma_Nr_2_n )
    else
       ! Mean rain water mixing ratio in PDF component 2 and (or) mean rain drop
       ! concentration in PDF component 2 are (is) less than their (its)
       ! respective tolerance amount(s), and are (is) considered to have a value
       ! of 0.  There is not any rain at this grid level.  The mean in-precip
       ! rain water mixing ratio (2nd PDF component) and (or) mean in-precip
       ! rain drop concentration (2nd PDF component) are also considered to have
       ! a value of 0.  The correlation is 0 since rain does not vary in this
       ! component at this grid level.
       corr_rrNr_2_n = zero
    endif


    return

  end subroutine normalize_pdf_params

  !=============================================================================
  subroutine pdf_param_hm_stats( mu_rr_1, mu_rr_2, mu_Nr_1, &
                                 mu_Nr_2, mu_Ncn_1, mu_Ncn_2, &
                                 sigma_rr_1, sigma_rr_2, sigma_Nr_1, &
                                 sigma_Nr_2, sigma_Ncn_1, sigma_Ncn_2, &
                                 corr_wrr_1, corr_wrr_2, corr_wNr_1, &
                                 corr_wNr_2, corr_wNcn_1, corr_wNcn_2, &
                                 corr_srr_1, corr_srr_2, corr_sNr_1, &
                                 corr_sNr_2, corr_sNcn_1, corr_sNcn_2, &
                                 corr_trr_1, corr_trr_2, corr_tNr_1, &
                                 corr_tNr_2, corr_tNcn_1, corr_tNcn_2, &
                                 corr_rrNr_1, corr_rrNr_2, &
                                 level, l_stats_samp )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd   ! Variable(s)

    use stats_type, only: &
        stat_update_var_pt  ! Procedure(s)

    use stats_variables, only : &
        imu_rr_1,       & ! Variable(s)
        imu_rr_2,       &
        imu_Nr_1,       &
        imu_Nr_2,       &
        imu_Ncn_1,      &
        imu_Ncn_2,      &
        isigma_rr_1,    &
        isigma_rr_2,    &
        isigma_Nr_1,    &
        isigma_Nr_2,    &
        isigma_Ncn_1,   &
        isigma_Ncn_2


    use stats_variables, only : &
        icorr_wrr_1,    & ! Variable(s)
        icorr_wrr_2,    &
        icorr_wNr_1,    &
        icorr_wNr_2,    &
        icorr_wNcn_1,   &
        icorr_wNcn_2,   &
        icorr_srr_1,    &
        icorr_srr_2,    &
        icorr_sNr_1,    &
        icorr_sNr_2,    &
        icorr_sNcn_1,   &
        icorr_sNcn_2,   &
        icorr_trr_1,    &
        icorr_trr_2,    &
        icorr_tNr_1,    &
        icorr_tNr_2,    &
        icorr_tNcn_1,   &
        icorr_tNcn_2,   &
        icorr_rrNr_1,   &
        icorr_rrNr_2, &
        zt

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_rr_1,       & ! Mean of rr (1st PDF component) in-precip (ip)   [kg/kg]
      mu_rr_2,       & ! Mean of rr (2nd PDF component) ip               [kg/kg]
      mu_Nr_1,       & ! Mean of Nr (1st PDF component) ip              [num/kg]
      mu_Nr_2,       & ! Mean of Nr (2nd PDF component) ip              [num/kg]
      mu_Ncn_1,      & ! Mean of Ncn (1st PDF component)                [num/kg]
      mu_Ncn_2,      & ! Mean of Ncn (2nd PDF component)                [num/kg]
      sigma_rr_1,    & ! Standard deviation of rr (1st PDF component) ip [kg/kg]
      sigma_rr_2,    & ! Standard deviation of rr (2nd PDF component) ip [kg/kg]
      sigma_Nr_1,    & ! Standard deviation of Nr (1st PDF comp.) ip    [num/kg]
      sigma_Nr_2,    & ! Standard deviation of Nr (2nd PDF comp.) ip    [num/kg]
      sigma_Ncn_1,   & ! Standard deviation of Ncn (1st PDF component)  [num/kg]
      sigma_Ncn_2      ! Standard deviation of Ncn (2nd PDF component)  [num/kg]


    real( kind = core_rknd ), intent(in) :: &
      corr_wrr_1,    & ! Correlation between w and rr (1st PDF component) ip [-]
      corr_wrr_2,    & ! Correlation between w and rr (2nd PDF component) ip [-]
      corr_wNr_1,    & ! Correlation between w and Nr (1st PDF component) ip [-]
      corr_wNr_2,    & ! Correlation between w and Nr (2nd PDF component) ip [-]
      corr_wNcn_1,   & ! Correlation between w and Ncn (1st PDF component)   [-]
      corr_wNcn_2,   & ! Correlation between w and Ncn (2nd PDF component)   [-]
      corr_srr_1,    & ! Correlation between s and rr (1st PDF component) ip [-]
      corr_srr_2,    & ! Correlation between s and rr (2nd PDF component) ip [-]
      corr_sNr_1,    & ! Correlation between s and Nr (1st PDF component) ip [-]
      corr_sNr_2,    & ! Correlation between s and Nr (2nd PDF component) ip [-]
      corr_sNcn_1,   & ! Correlation between s and Ncn (1st PDF component)   [-]
      corr_sNcn_2,   & ! Correlation between s and Ncn (2nd PDF component)   [-]
      corr_trr_1,    & ! Correlation between t and rr (1st PDF component) ip [-]
      corr_trr_2,    & ! Correlation between t and rr (2nd PDF component) ip [-]
      corr_tNr_1,    & ! Correlation between t and Nr (1st PDF component) ip [-]
      corr_tNr_2,    & ! Correlation between t and Nr (2nd PDF component) ip [-]
      corr_tNcn_1,   & ! Correlation between t and Ncn (1st PDF component)   [-]
      corr_tNcn_2,   & ! Correlation between t and Ncn (2nd PDF component)   [-]
      corr_rrNr_1,   & ! Correlation between rr & Nr (1st PDF component) ip  [-]
      corr_rrNr_2      ! Correlation between rr & Nr (2nd PDF component) ip  [-]

    integer, intent(in) :: &
      level   ! Vertical level index 

    logical, intent(in) :: &
      l_stats_samp     ! Flag to record statistical output.


    !!! Output the statistics for upscaled KK.

    ! Statistics
    if ( l_stats_samp ) then

       ! Mean of in-precip rain drop concentration in PDF component 1.
       if ( imu_rr_1 > 0 ) then
          call stat_update_var_pt( imu_rr_1, level, mu_rr_1, zt )
       endif

       ! Mean of in-precip rain drop concentration in PDF component 2.
       if ( imu_rr_2 > 0 ) then
          call stat_update_var_pt( imu_rr_2, level, mu_rr_2, zt )
       endif

       ! Mean of in-precip rain drop concentration in PDF component 1.
       if ( imu_Nr_1 > 0 ) then
          call stat_update_var_pt( imu_Nr_1, level, mu_Nr_1, zt )
       endif

       ! Mean of in-precip rain drop concentration in PDF component 2.
       if ( imu_Nr_2 > 0 ) then
          call stat_update_var_pt( imu_Nr_2, level, mu_Nr_2, zt )
       endif

       ! Mean of cloud nuclei concentration in PDF component 1.
       if ( imu_Ncn_1 > 0 ) then
          call stat_update_var_pt( imu_Ncn_1, level, mu_Ncn_1, zt )
       endif

       ! Mean of cloud nuclei concentration in PDF component 2.
       if ( imu_Ncn_2 > 0 ) then
          call stat_update_var_pt( imu_Ncn_2, level, mu_Ncn_2, zt )
       endif

       ! Standard deviation of in-precip rain water mixing ratio
       ! in PDF component 1.
       if ( isigma_rr_1 > 0 ) then
          call stat_update_var_pt( isigma_rr_1, level, sigma_rr_1, zt )
       endif

       ! Standard deviation of in-precip rain water mixing ratio
       ! in PDF component 2.
       if ( isigma_rr_2 > 0 ) then
          call stat_update_var_pt( isigma_rr_2, level, sigma_rr_2, zt )
       endif

       ! Standard deviation of in-precip rain drop concentration
       ! in PDF component 1.
       if ( isigma_Nr_1 > 0 ) then
          call stat_update_var_pt( isigma_Nr_1, level, sigma_Nr_1, zt )
       endif

       ! Standard deviation of in-precip rain drop concentration
       ! in PDF component 2.
       if ( isigma_Nr_2 > 0 ) then
          call stat_update_var_pt( isigma_Nr_2, level, sigma_Nr_2, zt )
       endif

       ! Standard deviation of cloud nuclei concentration in PDF component 1.
       if ( isigma_Ncn_1 > 0 ) then
          call stat_update_var_pt( isigma_Ncn_1, level, sigma_Ncn_1, zt )
       endif

       ! Standard deviation of cloud nuclei concentration in PDF component 2.
       if ( isigma_Ncn_2 > 0 ) then
          call stat_update_var_pt( isigma_Ncn_2, level, sigma_Ncn_2, zt )
       endif

       ! Correlation (in-precip) between w and r_r in PDF component 1.
       if ( icorr_wrr_1 > 0 ) then
          call stat_update_var_pt( icorr_wrr_1, level, corr_wrr_1, zt )
       endif

       ! Correlation (in-precip) between w and r_r in PDF component 2.
       if ( icorr_wrr_2 > 0 ) then
          call stat_update_var_pt( icorr_wrr_2, level, corr_wrr_2, zt )
       endif

       ! Correlation (in-precip) between w and N_r in PDF component 1.
       if ( icorr_wNr_1 > 0 ) then
          call stat_update_var_pt( icorr_wNr_1, level, corr_wNr_1, zt )
       endif

       ! Correlation (in-precip) between w and N_r in PDF component 2.
       if ( icorr_wNr_2 > 0 ) then
          call stat_update_var_pt( icorr_wNr_2, level, corr_wNr_2, zt )
       endif

       ! Correlation between w and N_cn in PDF component 1.
       if ( icorr_wNcn_1 > 0 ) then
          call stat_update_var_pt( icorr_wNcn_1, level, corr_wNcn_1, zt )
       endif

       ! Correlation between w and N_cn in PDF component 2.
       if ( icorr_wNcn_2 > 0 ) then
          call stat_update_var_pt( icorr_wNcn_2, level, corr_wNcn_2, zt )
       endif

       ! Correlation (in-precip) between s and r_r in PDF component 1.
       if ( icorr_srr_1 > 0 ) then
          call stat_update_var_pt( icorr_srr_1, level, corr_srr_1, zt )
       endif

       ! Correlation (in-precip) between s and r_r in PDF component 2.
       if ( icorr_srr_2 > 0 ) then
          call stat_update_var_pt( icorr_srr_2, level, corr_srr_2, zt )
       endif

       ! Correlation (in-precip) between s and N_r in PDF component 1.
       if ( icorr_sNr_1 > 0 ) then
          call stat_update_var_pt( icorr_sNr_1, level, corr_sNr_1, zt )
       endif

       ! Correlation (in-precip) between s and N_r in PDF component 2.
       if ( icorr_sNr_2 > 0 ) then
          call stat_update_var_pt( icorr_sNr_2, level, corr_sNr_2, zt )
       endif

       ! Correlation between s and N_cn in PDF component 1.
       if ( icorr_sNcn_1 > 0 ) then
          call stat_update_var_pt( icorr_sNcn_1, level, corr_sNcn_1, zt )
       endif

       ! Correlation between s and N_cn in PDF component 2.
       if ( icorr_sNcn_2 > 0 ) then
          call stat_update_var_pt( icorr_sNcn_2, level, corr_sNcn_2, zt )
       endif

       ! Correlation (in-precip) between t and r_r in PDF component 1.
       if ( icorr_trr_1 > 0 ) then
          call stat_update_var_pt( icorr_trr_1, level, corr_trr_1, zt )
       endif

       ! Correlation (in-precip) between t and r_r in PDF component 2.
       if ( icorr_trr_2 > 0 ) then
          call stat_update_var_pt( icorr_trr_2, level, corr_trr_2, zt )
       endif

       ! Correlation (in-precip) between t and N_r in PDF component 1.
       if ( icorr_tNr_1 > 0 ) then
          call stat_update_var_pt( icorr_tNr_1, level, corr_tNr_1, zt )
       endif

       ! Correlation (in-precip) between t and N_r in PDF component 2.
       if ( icorr_tNr_2 > 0 ) then
          call stat_update_var_pt( icorr_tNr_2, level, corr_tNr_2, zt )
       endif

       ! Correlation between t and N_cn in PDF component 1.
       if ( icorr_tNcn_1 > 0 ) then
          call stat_update_var_pt( icorr_tNcn_1, level, corr_tNcn_1, zt )
       endif

       ! Correlation between t and N_cn in PDF component 2.
       if ( icorr_tNcn_2 > 0 ) then
          call stat_update_var_pt( icorr_tNcn_2, level, corr_tNcn_2, zt )
       endif

       ! Correlation (in-precip) between r_r and N_r in PDF component 1.
       if ( icorr_rrNr_1 > 0 ) then
          call stat_update_var_pt( icorr_rrNr_1, level, corr_rrNr_1, zt )
       endif

       ! Correlation (in-precip) between r_r and N_r in PDF component 2.
       if ( icorr_rrNr_2 > 0 ) then
          call stat_update_var_pt( icorr_rrNr_2, level, corr_rrNr_2, zt )
       endif

    endif ! l_stats_samp

    return


  end subroutine pdf_param_hm_stats

  !=============================================================================
  subroutine pdf_param_log_hm_stats( mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, &
                                     mu_Nr_2_n, mu_Ncn_1_n, mu_Ncn_2_n, &
                                     sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, &
                                     sigma_Nr_2_n, sigma_Ncn_1_n, sigma_Ncn_2_n, &
                                     corr_wrr_1_n, corr_wrr_2_n, corr_wNr_1_n, &
                                     corr_wNr_2_n, corr_wNcn_1_n, corr_wNcn_2_n, &
                                     corr_srr_1_n, corr_srr_2_n, corr_sNr_1_n, &
                                     corr_sNr_2_n, corr_sNcn_1_n, corr_sNcn_2_n, &
                                     corr_trr_1_n, corr_trr_2_n, corr_tNr_1_n, &
                                     corr_tNr_2_n, corr_tNcn_1_n, corr_tNcn_2_n, &
                                     corr_rrNr_1_n, corr_rrNr_2_n, level, &
                                     l_stats_samp )

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd   ! Variable(s)

    use stats_type, only: &
        stat_update_var_pt  ! Procedure(s)

    use stats_variables, only : &
        imu_rr_1_n,     & ! Variable(s)
        imu_rr_2_n,     &
        imu_Nr_1_n,     &
        imu_Nr_2_n,     &
        imu_Ncn_1_n,    &
        imu_Ncn_2_n,    &
        isigma_rr_1_n,  &
        isigma_rr_2_n,  &
        isigma_Nr_1_n,  &
        isigma_Nr_2_n,  &
        isigma_Ncn_1_n, &
        isigma_Ncn_2_n

    use stats_variables, only : &
        icorr_wrr_1_n,  & ! Variables
        icorr_wrr_2_n,  &
        icorr_wNr_1_n,  &
        icorr_wNr_2_n,  &
        icorr_wNcn_1_n, &
        icorr_wNcn_2_n, &
        icorr_srr_1_n,  &
        icorr_srr_2_n,  &
        icorr_sNr_1_n,  &
        icorr_sNr_2_n,  &
        icorr_sNcn_1_n, &
        icorr_sNcn_2_n, &
        icorr_trr_1_n,  &
        icorr_trr_2_n,  &
        icorr_tNr_1_n,  &
        icorr_tNr_2_n,  &
        icorr_tNcn_1_n, &
        icorr_tNcn_2_n, &
        icorr_rrNr_1_n, &
        icorr_rrNr_2_n, &
        zt

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_rr_1_n,     & ! Mean of ln rr (1st PDF component) ip        [ln(kg/kg)]
      mu_rr_2_n,     & ! Mean of ln rr (2nd PDF component) ip        [ln(kg/kg)]
      mu_Nr_1_n,     & ! Mean of ln Nr (1st PDF component) ip       [ln(num/kg)]
      mu_Nr_2_n,     & ! Mean of ln Nr (2nd PDF component) ip       [ln(num/kg)]
      mu_Ncn_1_n,    & ! Mean of ln Ncn (1st PDF component)         [ln(num/kg)]
      mu_Ncn_2_n,    & ! Mean of ln Ncn (2nd PDF component)         [ln(num/kg)]
      sigma_rr_1_n,  & ! Standard dev. of ln rr (1st PDF comp.) ip   [ln(kg/kg)]
      sigma_rr_2_n,  & ! Standard dev. of ln rr (2nd PDF comp.) ip   [ln(kg/kg)]
      sigma_Nr_1_n,  & ! Standard dev. of ln Nr (1st PDF comp.) ip  [ln(num/kg)]
      sigma_Nr_2_n,  & ! Standard dev. of ln Nr (2nd PDF comp.) ip  [ln(num/kg)]
      sigma_Ncn_1_n, & ! Standard dev. of ln Ncn (1st PDF comp.)    [ln(num/kg)]
      sigma_Ncn_2_n    ! Standard dev. of ln Ncn (2nd PDF comp.)    [ln(num/kg)]

    real( kind = core_rknd ), intent(in) :: &
      corr_wrr_1_n,  & ! Correlation between w and ln rr (1st PDF comp.) ip  [-]
      corr_wrr_2_n,  & ! Correlation between w and ln rr (2nd PDF comp.) ip  [-]
      corr_wNr_1_n,  & ! Correlation between w and ln Nr (1st PDF comp.) ip  [-]
      corr_wNr_2_n,  & ! Correlation between w and ln Nr (2nd PDF comp.) ip  [-]
      corr_wNcn_1_n, & ! Correlation between w and ln Ncn (1st PDF comp.)    [-]
      corr_wNcn_2_n, & ! Correlation between w and ln Ncn (2nd PDF comp.)    [-]
      corr_srr_1_n,  & ! Correlation between s and ln rr (1st PDF comp.) ip  [-]
      corr_srr_2_n,  & ! Correlation between s and ln rr (2nd PDF comp.) ip  [-]
      corr_sNr_1_n,  & ! Correlation between s and ln Nr (1st PDF comp.) ip  [-]
      corr_sNr_2_n,  & ! Correlation between s and ln Nr (2nd PDF comp.) ip  [-]
      corr_sNcn_1_n, & ! Correlation between s and ln Ncn (1st PDF comp.)    [-]
      corr_sNcn_2_n, & ! Correlation between s and ln Ncn (2nd PDF comp.)    [-]
      corr_trr_1_n,  & ! Correlation between t and ln rr (1st PDF comp.) ip  [-]
      corr_trr_2_n,  & ! Correlation between t and ln rr (2nd PDF comp.) ip  [-]
      corr_tNr_1_n,  & ! Correlation between t and ln Nr (1st PDF comp.) ip  [-]
      corr_tNr_2_n,  & ! Correlation between t and ln Nr (2nd PDF comp.) ip  [-]
      corr_tNcn_1_n, & ! Correlation between t and ln Ncn (1st PDF comp.)    [-]
      corr_tNcn_2_n, & ! Correlation between t and ln Ncn (2nd PDF comp.)    [-]
      corr_rrNr_1_n, & ! Correlation btwn. ln rr & ln Nr (1st PDF comp.) ip  [-]
      corr_rrNr_2_n    ! Correlation btwn. ln rr & ln Nr (2nd PDF comp.) ip  [-]

    integer, intent(in) :: &
      level   ! Vertical level index

    logical, intent(in) :: &
      l_stats_samp     ! Flag to record statistical output.


    !!! Output the statistics for upscaled KK.

    ! Statistics
    if ( l_stats_samp ) then

       ! Mean (in-precip) of ln r_r in PDF component 1.
       if ( imu_rr_1_n > 0 ) then
          if ( mu_rr_1_n > -huge( 0.0 ) ) then
             call stat_update_var_pt( imu_rr_1_n, level, mu_rr_1_n, zt )
          else
             ! When rr1 is 0 (or below tolerance value), mu_rr_1_n is -inf, and
             ! is set to -huge for the default CLUBB kind.  Some compilers have
             ! issues outputting to stats files (in single precision) when the
             ! default CLUBB kind is in double precision.
             ! Set to -huge for single precision.
             call stat_update_var_pt( imu_rr_1_n, level, &
                                      real( -huge( 0.0 ), kind = core_rknd ), &
                                      zt )
          endif
       endif

       ! Mean (in-precip) of ln r_r in PDF component 2.
       if ( imu_rr_2_n > 0 ) then
          if ( mu_rr_2_n > -huge( 0.0 ) ) then
             call stat_update_var_pt( imu_rr_2_n, level, mu_rr_2_n, zt )
          else
             ! When rr2 is 0 (or below tolerance value), mu_rr_2_n is -inf, and
             ! is set to -huge for the default CLUBB kind.  Some compilers have
             ! issues outputting to stats files (in single precision) when the
             ! default CLUBB kind is in double precision.
             ! Set to -huge for single precision.
             call stat_update_var_pt( imu_rr_2_n, level, &
                                      real( -huge( 0.0 ), kind = core_rknd ), &
                                      zt )
          endif
       endif

       ! Mean (in-precip) of ln N_r in PDF component 1.
       if ( imu_Nr_1_n > 0 ) then
          if ( mu_Nr_1_n > -huge( 0.0 ) ) then
             call stat_update_var_pt( imu_Nr_1_n, level, mu_Nr_1_n, zt )
          else
             ! When Nr1 is 0 (or below tolerance value), mu_Nr_1_n is -inf, and
             ! is set to -huge for the default CLUBB kind.  Some compilers have
             ! issues outputting to stats files (in single precision) when the
             ! default CLUBB kind is in double precision.
             ! Set to -huge for single precision.
             call stat_update_var_pt( imu_Nr_1_n, level, &
                                      real( -huge( 0.0 ), kind = core_rknd ), &
                                      zt )
          endif
       endif

       ! Mean (in-precip) of ln N_r in PDF component 2.
       if ( imu_Nr_2_n > 0 ) then
          if ( mu_Nr_2_n > -huge( 0.0 ) ) then
             call stat_update_var_pt( imu_Nr_2_n, level, mu_Nr_2_n, zt )
          else
             ! When Nr2 is 0 (or below tolerance value), mu_Nr_2_n is -inf, and
             ! is set to -huge for the default CLUBB kind.  Some compilers have
             ! issues outputting to stats files (in single precision) when the
             ! default CLUBB kind is in double precision.
             ! Set to -huge for single precision.
             call stat_update_var_pt( imu_Nr_2_n, level, &
                                      real( -huge( 0.0 ), kind = core_rknd ), &
                                      zt )
          endif
       endif

       ! Mean of ln N_cn in PDF component 1.
       if ( imu_Ncn_1_n > 0 ) then
          if ( mu_Ncn_1_n > -huge( 0.0 ) ) then
             call stat_update_var_pt( imu_Ncn_1_n, level, mu_Ncn_1_n, zt )
          else
             ! When Ncnm is 0 (or below tolerance value), mu_Ncn_1_n is -inf,
             ! and is set to -huge for the default CLUBB kind.  Some compilers
             ! have issues outputting to stats files (in single precision) when
             ! the default CLUBB kind is in double precision.
             ! Set to -huge for single precision.
             call stat_update_var_pt( imu_Ncn_1_n, level, &
                                      real( -huge( 0.0 ), kind = core_rknd ), &
                                      zt )
          endif
       endif

       ! Mean of ln N_cn in PDF component 2.
       if ( imu_Ncn_2_n > 0 ) then
          if ( mu_Ncn_2_n > -huge( 0.0 ) ) then
             call stat_update_var_pt( imu_Ncn_2_n, level, mu_Ncn_2_n, zt )
          else
             ! When Ncnm is 0 (or below tolerance value), mu_Ncn_2_n is -inf,
             ! and is set to -huge for the default CLUBB kind.  Some compilers
             ! have issues outputting to stats files (in single precision) when
             ! the default CLUBB kind is in double precision.
             ! Set to -huge for single precision.
             call stat_update_var_pt( imu_Ncn_2_n, level, &
                                      real( -huge( 0.0 ), kind = core_rknd ), &
                                      zt )
          endif
       endif

       ! Standard deviation (in-precip) of ln r_r in PDF component 1.
       if ( isigma_rr_1_n > 0 ) then
          call stat_update_var_pt( isigma_rr_1_n, level, sigma_rr_1_n, zt )
       endif

       ! Standard deviation (in-precip) of ln r_r in PDF component 2.
       if ( isigma_rr_2_n > 0 ) then
          call stat_update_var_pt( isigma_rr_2_n, level, sigma_rr_2_n, zt )
       endif

       ! Standard deviation (in-precip) of ln N_r in PDF component 1.
       if ( isigma_Nr_1_n > 0 ) then
          call stat_update_var_pt( isigma_Nr_1_n, level, sigma_Nr_1_n, zt )
       endif

       ! Standard deviation (in-precip) of ln N_r in PDF component 2.
       if ( isigma_Nr_2_n > 0 ) then
          call stat_update_var_pt( isigma_Nr_2_n, level, sigma_Nr_2_n, zt )
       endif

       ! Standard deviation of ln N_cn in PDF component 1.
       if ( isigma_Ncn_1_n > 0 ) then
          call stat_update_var_pt( isigma_Ncn_1_n, level, sigma_Ncn_1_n, zt )
       endif

       ! Standard deviation of ln N_cn in PDF component 2.
       if ( isigma_Ncn_2_n > 0 ) then
          call stat_update_var_pt( isigma_Ncn_2_n, level, sigma_Ncn_2_n, zt )
       endif

       ! Correlation (in-precip) between w and ln r_r in PDF component 1.
       if ( icorr_wrr_1_n > 0 ) then
          call stat_update_var_pt( icorr_wrr_1_n, level, corr_wrr_1_n, zt )
       endif

       ! Correlation (in-precip) between w and ln r_r in PDF component 2.
       if ( icorr_wrr_2_n > 0 ) then
          call stat_update_var_pt( icorr_wrr_2_n, level, corr_wrr_2_n, zt )
       endif

       ! Correlation (in-precip) between w and ln N_r in PDF component 1.
       if ( icorr_wNr_1_n > 0 ) then
          call stat_update_var_pt( icorr_wNr_1_n, level, corr_wNr_1_n, zt )
       endif

       ! Correlation (in-precip) between w and ln N_r in PDF component 2.
       if ( icorr_wNr_2_n > 0 ) then
          call stat_update_var_pt( icorr_wNr_2_n, level, corr_wNr_2_n, zt )
       endif

       ! Correlation between w and ln N_cn in PDF component 1.
       if ( icorr_wNcn_1_n > 0 ) then
          call stat_update_var_pt( icorr_wNcn_1_n, level, corr_wNcn_1_n, zt )
       endif

       ! Correlation between w and ln N_cn in PDF component 2.
       if ( icorr_wNcn_2_n > 0 ) then
          call stat_update_var_pt( icorr_wNcn_2_n, level, corr_wNcn_2_n, zt )
       endif

       ! Correlation (in-precip) between s and ln r_r in PDF component 1.
       if ( icorr_srr_1_n > 0 ) then
          call stat_update_var_pt( icorr_srr_1_n, level, corr_srr_1_n, zt )
       endif

       ! Correlation (in-precip) between s and ln r_r in PDF component 2.
       if ( icorr_srr_2_n > 0 ) then
          call stat_update_var_pt( icorr_srr_2_n, level, corr_srr_2_n, zt )
       endif

       ! Correlation (in-precip) between s and ln N_r in PDF component 1.
       if ( icorr_sNr_1_n > 0 ) then
          call stat_update_var_pt( icorr_sNr_1_n, level, corr_sNr_1_n, zt )
       endif

       ! Correlation (in-precip) between s and ln N_r in PDF component 2.
       if ( icorr_sNr_2_n > 0 ) then
          call stat_update_var_pt( icorr_sNr_2_n, level, corr_sNr_2_n, zt )
       endif

       ! Correlation between s and ln N_cn in PDF component 1.
       if ( icorr_sNcn_1_n > 0 ) then
          call stat_update_var_pt( icorr_sNcn_1_n, level, corr_sNcn_1_n, zt )
       endif

       ! Correlation between s and ln N_cn in PDF component 2.
       if ( icorr_sNcn_2_n > 0 ) then
          call stat_update_var_pt( icorr_sNcn_2_n, level, corr_sNcn_2_n, zt )
       endif

       ! Correlation (in-precip) between t and ln r_r in PDF component 1.
       if ( icorr_trr_1_n > 0 ) then
          call stat_update_var_pt( icorr_trr_1_n, level, corr_trr_1_n, zt )
       endif

       ! Correlation (in-precip) between t and ln r_r in PDF component 2.
       if ( icorr_trr_2_n > 0 ) then
          call stat_update_var_pt( icorr_trr_2_n, level, corr_trr_2_n, zt )
       endif

       ! Correlation (in-precip) between t and ln N_r in PDF component 1.
       if ( icorr_tNr_1_n > 0 ) then
          call stat_update_var_pt( icorr_tNr_1_n, level, corr_tNr_1_n, zt )
       endif

       ! Correlation (in-precip) between t and ln N_r in PDF component 2.
       if ( icorr_tNr_2_n > 0 ) then
          call stat_update_var_pt( icorr_tNr_2_n, level, corr_tNr_2_n, zt )
       endif

       ! Correlation between t and ln N_cn in PDF component 1.
       if ( icorr_tNcn_1_n > 0 ) then
          call stat_update_var_pt( icorr_tNcn_1_n, level, corr_tNcn_1_n, zt )
       endif

       ! Correlation between t and ln N_cn in PDF component 2.
       if ( icorr_tNcn_2_n > 0 ) then
          call stat_update_var_pt( icorr_tNcn_2_n, level, corr_tNcn_2_n, zt )
       endif

       ! Correlation (in-precip) between ln r_r and ln N_r in PDF component 1.
       if ( icorr_rrNr_1_n > 0 ) then
          call stat_update_var_pt( icorr_rrNr_1_n, level, corr_rrNr_1_n, zt )
       endif

       ! Correlation (in-precip) between ln r_r and ln N_r in PDF component 2.
       if ( icorr_rrNr_2_n > 0 ) then
          call stat_update_var_pt( icorr_rrNr_2_n, level, corr_rrNr_2_n, zt )
       endif

    endif ! l_stats_samp

    return


  end subroutine pdf_param_log_hm_stats

  !=============================================================================
  subroutine pack_pdf_params( mu_w_1, mu_w_2, mu_s_1, mu_s_2, &             ! Intent(in)
                              mu_t_1, mu_t_2, mu_rr_1, mu_rr_2, &           ! Intent(in)
                              mu_Nr_1, mu_Nr_2, mu_Ncn_1, mu_Ncn_2, &       ! Intent(in)
                              mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, &            ! Intent(in)
                              mu_Nr_2_n, mu_Ncn_1_n, mu_Ncn_2_n, &          ! Intent(in)
                              sigma_w_1, sigma_w_2, sigma_s_1, &            ! Intent(in)
                              sigma_s_2, sigma_t_1, sigma_t_2, &            ! Intent(in)
                              sigma_rr_1, sigma_rr_2, sigma_Nr_1, &         ! Intent(in)
                              sigma_Nr_2, sigma_Ncn_1, sigma_Ncn_2, &       ! Intent(in)
                              sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, &   ! Intent(in)
                              sigma_Nr_2_n, sigma_Ncn_1_n, sigma_Ncn_2_n, & ! Intent(in)
                              corr_ws_1, corr_ws_2, corr_st_1, corr_st_2, & ! Intent(in)
                              corr_wrr_1_n, corr_wrr_2_n, corr_wNr_1_n, &   ! Intent(in)
                              corr_wNr_2_n, corr_wNcn_1_n, corr_wNcn_2_n, & ! Intent(in)
                              corr_srr_1_n, corr_srr_2_n, corr_sNr_1_n, &   ! Intent(in)
                              corr_sNr_2_n, corr_sNcn_1_n, corr_sNcn_2_n, & ! Intent(in)
                              corr_trr_1_n, corr_trr_2_n, corr_tNr_1_n, &   ! Intent(in)
                              corr_tNr_2_n, corr_tNcn_1_n, corr_tNcn_2_n, & ! Intent(in)
                              corr_rrNr_1_n, corr_rrNr_2_n, &               ! Intent(in)
                              rr1, rr2, Nr1, Nr2, &                         ! Intent(in)
                              precip_frac, precip_frac_1, precip_frac_2, &  ! Intent(in)
                              d_variables, &                                ! Intent(in)
                              corr_array_1, corr_array_2, &                 ! Intent(inout)
                              mu_x_1, mu_x_2, sigma_x_1, sigma_x_2, &       ! Intent(out)
                              hydromet_pdf_params )                         ! Intent(out)

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use hydromet_pdf_parameter_module, only: &
        hydromet_pdf_parameter  ! Variable(s)

    use corr_matrix_module, only: &
        iiLH_w,        & ! Variable(s)
        iiLH_s_mellor, &
        iiLH_t_mellor, &
        iiLH_rrain,    &
        iiLH_Nr,       &
        iiLH_Ncn

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_w_1,        & ! Mean of w (1st PDF component)                     [m/s]
      mu_w_2,        & ! Mean of w (2nd PDF component)                     [m/s]
      mu_s_1,        & ! Mean of s (1st PDF component)                   [kg/kg]
      mu_s_2,        & ! Mean of s (2nd PDF component)                   [kg/kg]
      mu_t_1,        & ! Mean of t (1st PDF component)                   [kg/kg]
      mu_t_2,        & ! Mean of t (2nd PDF component)                   [kg/kg]
      mu_rr_1,       & ! Mean of rr (1st PDF component) in-precip (ip)   [kg/kg]
      mu_rr_2,       & ! Mean of rr (2nd PDF component) ip               [kg/kg]
      mu_Nr_1,       & ! Mean of Nr (1st PDF component) ip              [num/kg]
      mu_Nr_2,       & ! Mean of Nr (2nd PDF component) ip              [num/kg]
      mu_Ncn_1,      & ! Mean of Ncn (1st PDF component)                [num/kg]
      mu_Ncn_2,      & ! Mean of Ncn (2nd PDF component)                [num/kg]
      mu_rr_1_n,     & ! Mean of ln rr (1st PDF component) ip        [ln(kg/kg)]
      mu_rr_2_n,     & ! Mean of ln rr (2nd PDF component) ip        [ln(kg/kg)]
      mu_Nr_1_n,     & ! Mean of ln Nr (1st PDF component) ip       [ln(num/kg)]
      mu_Nr_2_n,     & ! Mean of ln Nr (2nd PDF component) ip       [ln(num/kg)]
      mu_Ncn_1_n,    & ! Mean of ln Ncn (1st PDF component)         [ln(num/kg)]
      mu_Ncn_2_n,    & ! Mean of ln Ncn (2nd PDF component)         [ln(num/kg)]
      sigma_w_1,     & ! Standard deviation of w (1st PDF component)       [m/s]
      sigma_w_2,     & ! Standard deviation of w (2nd PDF component)       [m/s]
      sigma_s_1,     & ! Standard deviation of s (1st PDF component)     [kg/kg]
      sigma_s_2,     & ! Standard deviation of s (2nd PDF component)     [kg/kg]
      sigma_t_1,     & ! Standard deviation of t (1st PDF component)     [kg/kg]
      sigma_t_2,     & ! Standard deviation of t (2nd PDF component)     [kg/kg]
      sigma_rr_1,    & ! Standard deviation of rr (1st PDF component) ip [kg/kg]
      sigma_rr_2,    & ! Standard deviation of rr (2nd PDF component) ip [kg/kg]
      sigma_Nr_1,    & ! Standard deviation of Nr (1st PDF comp.) ip    [num/kg]
      sigma_Nr_2,    & ! Standard deviation of Nr (2nd PDF comp.) ip    [num/kg]
      sigma_Ncn_1,   & ! Standard deviation of Ncn (1st PDF component)  [num/kg]
      sigma_Ncn_2,   & ! Standard deviation of Ncn (2nd PDF component)  [num/kg]
      sigma_rr_1_n,  & ! Standard dev. of ln rr (1st PDF comp.) ip   [ln(kg/kg)]
      sigma_rr_2_n,  & ! Standard dev. of ln rr (2nd PDF comp.) ip   [ln(kg/kg)]
      sigma_Nr_1_n,  & ! Standard dev. of ln Nr (1st PDF comp.) ip  [ln(num/kg)]
      sigma_Nr_2_n,  & ! Standard dev. of ln Nr (2nd PDF comp.) ip  [ln(num/kg)]
      sigma_Ncn_1_n, & ! Standard dev. of ln Ncn (1st PDF comp.)    [ln(num/kg)]
      sigma_Ncn_2_n    ! Standard dev. of ln Ncn (2nd PDF comp.)    [ln(num/kg)]

    real( kind = core_rknd ), intent(in) :: &
      corr_ws_1,     & ! Correlation between w and s (1st PDF component)     [-]
      corr_ws_2,     & ! Correlation between w and s (2nd PDF component)     [-]
      corr_st_1,     & ! Correlation between s and t (1st PDF component)     [-]
      corr_st_2        ! Correlation between s and t (2nd PDF component)     [-]

    real( kind = core_rknd ), intent(in) :: &
      corr_wrr_1_n,  & ! Correlation between w and ln rr (1st PDF comp.) ip  [-]
      corr_wrr_2_n,  & ! Correlation between w and ln rr (2nd PDF comp.) ip  [-]
      corr_wNr_1_n,  & ! Correlation between w and ln Nr (1st PDF comp.) ip  [-]
      corr_wNr_2_n,  & ! Correlation between w and ln Nr (2nd PDF comp.) ip  [-]
      corr_wNcn_1_n, & ! Correlation between w and ln Ncn (1st PDF comp.)    [-]
      corr_wNcn_2_n, & ! Correlation between w and ln Ncn (2nd PDF comp.)    [-]
      corr_srr_1_n,  & ! Correlation between s and ln rr (1st PDF comp.) ip  [-]
      corr_srr_2_n,  & ! Correlation between s and ln rr (2nd PDF comp.) ip  [-]
      corr_sNr_1_n,  & ! Correlation between s and ln Nr (1st PDF comp.) ip  [-]
      corr_sNr_2_n,  & ! Correlation between s and ln Nr (2nd PDF comp.) ip  [-]
      corr_sNcn_1_n, & ! Correlation between s and ln Ncn (1st PDF comp.)    [-]
      corr_sNcn_2_n, & ! Correlation between s and ln Ncn (2nd PDF comp.)    [-]
      corr_trr_1_n,  & ! Correlation between t and ln rr (1st PDF comp.) ip  [-]
      corr_trr_2_n,  & ! Correlation between t and ln rr (2nd PDF comp.) ip  [-]
      corr_tNr_1_n,  & ! Correlation between t and ln Nr (1st PDF comp.) ip  [-]
      corr_tNr_2_n,  & ! Correlation between t and ln Nr (2nd PDF comp.) ip  [-]
      corr_tNcn_1_n, & ! Correlation between t and ln Ncn (1st PDF comp.)    [-]
      corr_tNcn_2_n, & ! Correlation between t and ln Ncn (2nd PDF comp.)    [-]
      corr_rrNr_1_n, & ! Correlation btwn. ln rr & ln Nr (1st PDF comp.) ip  [-]
      corr_rrNr_2_n    ! Correlation btwn. ln rr & ln Nr (2nd PDF comp.) ip  [-]

    real( kind = core_rknd ), intent(in) :: &
      rr1, & ! Mean rain water mixing ratio (1st PDF component)      [kg/kg]
      rr2, & ! Mean rain water mixing ratio (2nd PDF component)      [kg/kg]
      Nr1, & ! Mean rain drop concentration (1st PDF component)      [num/kg]
      Nr2    ! Mean rain drop concentration (2nd PDF component)      [num/kg]

    real( kind = core_rknd ), intent(in) :: &
      precip_frac,   & ! Precipitation fraction (overall)           [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component) [-]
      precip_frac_2    ! Precipitation fraction (2nd PDF component) [-]


    integer, intent(in) :: &
      d_variables    ! Number of variables in the correlation array.

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(d_variables,d_variables), &
    intent(inout) :: &
      corr_array_1, & ! Correlation array for the 1st PDF component   [-]
      corr_array_2    ! Correlation array for the 2nd PDF component   [-]

    ! Output Variables
    real( kind = core_rknd ), dimension(d_variables), intent(out) :: &
      mu_x_1,    & ! Mean array for the 1st PDF component                 [units vary]
      mu_x_2,    & ! Mean array for the 2nd PDF component                 [units vary]
      sigma_x_1, & ! Standard deviation array for the 1st PDF component   [units vary]
      sigma_x_2    ! Standard deviation array for the 2nd PDF component   [units vary]

    type(hydromet_pdf_parameter), intent(out) :: &
      hydromet_pdf_params    ! Hydrometeor PDF parameters        [units vary]

    ! ---- Begin Code ----

    ! Pack Means and Standard Deviations into arrays
    mu_x_1(iiLh_w)        = mu_w_1
    mu_x_2(iiLh_w)        = mu_w_2
    mu_x_1(iiLh_s_mellor) = mu_s_1
    mu_x_2(iiLh_s_mellor) = mu_s_2
    mu_x_1(iiLh_t_mellor) = mu_t_1
    mu_x_2(iiLh_t_mellor) = mu_t_2
    mu_x_1(iiLh_rrain)    = mu_rr_1_n
    mu_x_2(iiLh_rrain)    = mu_rr_2_n
    mu_x_1(iiLh_Nr)       = mu_Nr_1_n
    mu_x_2(iiLh_Nr)       = mu_Nr_2_n
    mu_x_1(iiLh_Ncn)      = mu_Ncn_1_n
    mu_x_2(iiLh_Ncn)      = mu_Ncn_2_n
    sigma_x_1(iiLh_w)        = sigma_w_1
    sigma_x_2(iiLh_w)        = sigma_w_2
    sigma_x_1(iiLh_s_mellor) = sigma_s_1
    sigma_x_2(iiLh_s_mellor) = sigma_s_2
    sigma_x_1(iiLh_t_mellor) = sigma_t_1
    sigma_x_2(iiLh_t_mellor) = sigma_t_2
    sigma_x_1(iiLh_rrain)    = sigma_rr_1_n
    sigma_x_2(iiLh_rrain)    = sigma_rr_2_n
    sigma_x_1(iiLh_Nr)       = sigma_Nr_1_n
    sigma_x_2(iiLh_Nr)       = sigma_Nr_2_n
    sigma_x_1(iiLh_Ncn)      = sigma_Ncn_1_n
    sigma_x_2(iiLh_Ncn)      = sigma_Ncn_2_n

    ! Pack remaining variables into hydromet_pdf_params
    hydromet_pdf_params%mu_rr_1     = mu_rr_1
    hydromet_pdf_params%mu_rr_2     = mu_rr_2
    hydromet_pdf_params%mu_Nr_1     = mu_Nr_1
    hydromet_pdf_params%mu_Nr_2     = mu_Nr_2
    hydromet_pdf_params%mu_Ncn_1    = mu_Ncn_1
    hydromet_pdf_params%mu_Ncn_2    = mu_Ncn_2
    hydromet_pdf_params%sigma_rr_1  = sigma_rr_1
    hydromet_pdf_params%sigma_rr_2  = sigma_rr_2
    hydromet_pdf_params%sigma_Nr_1  = sigma_Nr_1
    hydromet_pdf_params%sigma_Nr_2  = sigma_Nr_2
    hydromet_pdf_params%sigma_Ncn_1 = sigma_Ncn_1
    hydromet_pdf_params%sigma_Ncn_2 = sigma_Ncn_2

    hydromet_pdf_params%rr1           = rr1
    hydromet_pdf_params%rr2           = rr2
    hydromet_pdf_params%Nr1           = Nr1
    hydromet_pdf_params%Nr2           = Nr2
    hydromet_pdf_params%precip_frac   = precip_frac
    hydromet_pdf_params%precip_frac_1 = precip_frac_1
    hydromet_pdf_params%precip_frac_2 = precip_frac_2

    ! Pack correlations (1st PDF component) into corr_array_1.
    corr_array_1(iiLH_t_mellor, iiLH_s_mellor) = corr_st_1
    corr_array_1(iiLH_s_mellor,iiLH_w)        = corr_ws_1
    corr_array_1(iiLH_rrain, iiLH_s_mellor)    = corr_srr_1_n
    corr_array_1(iiLH_Nr, iiLH_s_mellor)       = corr_sNr_1_n
    corr_array_1(iiLH_Ncn, iiLH_s_mellor)      = corr_sNcn_1_n
    corr_array_1(iiLH_rrain, iiLH_t_mellor)    = corr_trr_1_n
    corr_array_1(iiLH_Nr, iiLH_t_mellor)       = corr_tNr_1_n
    corr_array_1(iiLH_Ncn, iiLH_t_mellor)      = corr_tNcn_1_n
    corr_array_1(iiLH_rrain, iiLH_w)           = corr_wrr_1_n
    corr_array_1(iiLH_Nr, iiLH_w)              = corr_wNr_1_n
    corr_array_1(iiLH_Ncn, iiLH_w)             = corr_wNcn_1_n
    corr_array_1(iiLH_rrain,iiLH_Nr)          = corr_rrNr_1_n

    ! Pack correlations (2nd PDF component) into corr_array_2.
    corr_array_2(iiLH_t_mellor, iiLH_s_mellor) = corr_st_2
    corr_array_2(iiLH_s_mellor,iiLH_w)        = corr_ws_2
    corr_array_2(iiLH_rrain, iiLH_s_mellor)    = corr_srr_2_n
    corr_array_2(iiLH_Nr, iiLH_s_mellor)       = corr_sNr_2_n
    corr_array_2(iiLH_Ncn, iiLH_s_mellor)      = corr_sNcn_2_n
    corr_array_2(iiLH_rrain, iiLH_t_mellor)    = corr_trr_2_n
    corr_array_2(iiLH_Nr, iiLH_t_mellor)       = corr_tNr_2_n
    corr_array_2(iiLH_Ncn, iiLH_t_mellor)      = corr_tNcn_2_n
    corr_array_2(iiLH_rrain, iiLH_w)           = corr_wrr_2_n
    corr_array_2(iiLH_Nr, iiLH_w)              = corr_wNr_2_n
    corr_array_2(iiLH_Ncn, iiLH_w)             = corr_wNcn_2_n
    corr_array_2(iiLH_rrain,iiLH_Nr)          = corr_rrNr_2_n

    return

  end subroutine pack_pdf_params

  !=============================================================================
  subroutine unpack_pdf_params( d_variables, corr_array_1, corr_array_2, &    ! Intent(in)
                                mu_x_1, mu_x_2, sigma_x_1, sigma_x_2, &       ! Intent(in)
                                hydromet_pdf_params, &                        ! Intent(in)
                                mu_w_1, mu_w_2, mu_s_1, mu_s_2, &             ! Intent(out)
                                mu_t_1, mu_t_2, mu_rr_1, mu_rr_2, &           ! Intent(out)
                                mu_Nr_1, mu_Nr_2, mu_Ncn_1, mu_Ncn_2, &       ! Intent(out)
                                mu_rr_1_n, mu_rr_2_n, mu_Nr_1_n, &            ! Intent(out)
                                mu_Nr_2_n, mu_Ncn_1_n, mu_Ncn_2_n, &          ! Intent(out)
                                sigma_w_1, sigma_w_2, sigma_s_1, &            ! Intent(out)
                                sigma_s_2, sigma_t_1, sigma_t_2, &            ! Intent(out)
                                sigma_rr_1, sigma_rr_2, sigma_Nr_1, &         ! Intent(out)
                                sigma_Nr_2, sigma_Ncn_1, sigma_Ncn_2, &       ! Intent(out)
                                sigma_rr_1_n, sigma_rr_2_n, sigma_Nr_1_n, &   ! Intent(out)
                                sigma_Nr_2_n, sigma_Ncn_1_n, sigma_Ncn_2_n, & ! Intent(out)
                                corr_ws_1, corr_ws_2, corr_st_1, corr_st_2, & ! Intent(out)
                                corr_wrr_1_n, corr_wrr_2_n, corr_wNr_1_n, &   ! Intent(out)
                                corr_wNr_2_n, corr_wNcn_1_n, corr_wNcn_2_n, & ! Intent(out)
                                corr_srr_1_n, corr_srr_2_n, corr_sNr_1_n, &   ! Intent(out)
                                corr_sNr_2_n, corr_sNcn_1_n, corr_sNcn_2_n, & ! Intent(out)
                                corr_trr_1_n, corr_trr_2_n, corr_tNr_1_n, &   ! Intent(out)
                                corr_tNr_2_n, corr_tNcn_1_n, corr_tNcn_2_n, & ! Intent(out)
                                corr_rrNr_1_n, corr_rrNr_2_n , &              ! Intent(out)
                                rr1, rr2, Nr1, Nr2, &                         ! Intent(out)
                                precip_frac, precip_frac_1, precip_frac_2 )   ! Intent(out)

    ! Description:

    ! References:
    !-----------------------------------------------------------------------

    use hydromet_pdf_parameter_module, only: &
        hydromet_pdf_parameter  ! Variable(s)

    use corr_matrix_module, only: &
        iiLH_w,        & ! Variable(s)
        iiLH_s_mellor, &
        iiLH_t_mellor, &
        iiLH_rrain,    &
        iiLH_Nr,       &
        iiLH_Ncn

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      d_variables    ! Number of variables in the correlation array.

    real( kind = core_rknd ), dimension(d_variables,d_variables), &
    intent(in) :: &
      corr_array_1, & ! Correlation array for the 1st PDF component   [-]
      corr_array_2    ! Correlation array for the 2nd PDF component   [-]

    real( kind = core_rknd ), dimension(d_variables), intent(in) :: &
      mu_x_1,    & ! Mean array for the 1st PDF component                 [units vary]
      mu_x_2,    & ! Mean array for the 2nd PDF component                 [units vary]
      sigma_x_1, & ! Standard deviation array for the 1st PDF component   [units vary]
      sigma_x_2    ! Standard deviation array for the 2nd PDF component   [units vary]

    type(hydromet_pdf_parameter), intent(in) :: &
      hydromet_pdf_params    ! Hydrometeor PDF parameters        [units vary]

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
      mu_w_1,        & ! Mean of w (1st PDF component)                     [m/s]
      mu_w_2,        & ! Mean of w (2nd PDF component)                     [m/s]
      mu_s_1,        & ! Mean of s (1st PDF component)                   [kg/kg]
      mu_s_2,        & ! Mean of s (2nd PDF component)                   [kg/kg]
      mu_t_1,        & ! Mean of t (1st PDF component)                   [kg/kg]
      mu_t_2,        & ! Mean of t (2nd PDF component)                   [kg/kg]
      mu_rr_1,       & ! Mean of rr (1st PDF component) in-precip (ip)   [kg/kg]
      mu_rr_2,       & ! Mean of rr (2nd PDF component) ip               [kg/kg]
      mu_Nr_1,       & ! Mean of Nr (1st PDF component) ip              [num/kg]
      mu_Nr_2,       & ! Mean of Nr (2nd PDF component) ip              [num/kg]
      mu_Ncn_1,      & ! Mean of Ncn (1st PDF component)                [num/kg]
      mu_Ncn_2,      & ! Mean of Ncn (2nd PDF component)                [num/kg]
      mu_rr_1_n,     & ! Mean of ln rr (1st PDF component) ip        [ln(kg/kg)]
      mu_rr_2_n,     & ! Mean of ln rr (2nd PDF component) ip        [ln(kg/kg)]
      mu_Nr_1_n,     & ! Mean of ln Nr (1st PDF component) ip       [ln(num/kg)]
      mu_Nr_2_n,     & ! Mean of ln Nr (2nd PDF component) ip       [ln(num/kg)]
      mu_Ncn_1_n,    & ! Mean of ln Ncn (1st PDF component)         [ln(num/kg)]
      mu_Ncn_2_n,    & ! Mean of ln Ncn (2nd PDF component)         [ln(num/kg)]
      sigma_w_1,     & ! Standard deviation of w (1st PDF component)       [m/s]
      sigma_w_2,     & ! Standard deviation of w (2nd PDF component)       [m/s]
      sigma_s_1,     & ! Standard deviation of s (1st PDF component)     [kg/kg]
      sigma_s_2,     & ! Standard deviation of s (2nd PDF component)     [kg/kg]
      sigma_t_1,     & ! Standard deviation of t (1st PDF component)     [kg/kg]
      sigma_t_2,     & ! Standard deviation of t (2nd PDF component)     [kg/kg]
      sigma_rr_1,    & ! Standard deviation of rr (1st PDF component) ip [kg/kg]
      sigma_rr_2,    & ! Standard deviation of rr (2nd PDF component) ip [kg/kg]
      sigma_Nr_1,    & ! Standard deviation of Nr (1st PDF comp.) ip    [num/kg]
      sigma_Nr_2,    & ! Standard deviation of Nr (2nd PDF comp.) ip    [num/kg]
      sigma_Ncn_1,   & ! Standard deviation of Ncn (1st PDF component)  [num/kg]
      sigma_Ncn_2,   & ! Standard deviation of Ncn (2nd PDF component)  [num/kg]
      sigma_rr_1_n,  & ! Standard dev. of ln rr (1st PDF comp.) ip   [ln(kg/kg)]
      sigma_rr_2_n,  & ! Standard dev. of ln rr (2nd PDF comp.) ip   [ln(kg/kg)]
      sigma_Nr_1_n,  & ! Standard dev. of ln Nr (1st PDF comp.) ip  [ln(num/kg)]
      sigma_Nr_2_n,  & ! Standard dev. of ln Nr (2nd PDF comp.) ip  [ln(num/kg)]
      sigma_Ncn_1_n, & ! Standard dev. of ln Ncn (1st PDF comp.)    [ln(num/kg)]
      sigma_Ncn_2_n    ! Standard dev. of ln Ncn (2nd PDF comp.)    [ln(num/kg)]

    real( kind = core_rknd ), intent(out) :: &
      corr_ws_1,     & ! Correlation between w and s (1st PDF component)     [-]
      corr_ws_2,     & ! Correlation between w and s (2nd PDF component)     [-]
      corr_st_1,     & ! Correlation between s and t (1st PDF component)     [-]
      corr_st_2        ! Correlation between s and t (2nd PDF component)     [-]

   real( kind = core_rknd ), intent(out) :: &
      corr_wrr_1_n,  & ! Correlation between w and ln rr (1st PDF comp.) ip  [-]
      corr_wrr_2_n,  & ! Correlation between w and ln rr (2nd PDF comp.) ip  [-]
      corr_wNr_1_n,  & ! Correlation between w and ln Nr (1st PDF comp.) ip  [-]
      corr_wNr_2_n,  & ! Correlation between w and ln Nr (2nd PDF comp.) ip  [-]
      corr_wNcn_1_n, & ! Correlation between w and ln Ncn (1st PDF comp.)    [-]
      corr_wNcn_2_n, & ! Correlation between w and ln Ncn (2nd PDF comp.)    [-]
      corr_srr_1_n,  & ! Correlation between s and ln rr (1st PDF comp.) ip  [-]
      corr_srr_2_n,  & ! Correlation between s and ln rr (2nd PDF comp.) ip  [-]
      corr_sNr_1_n,  & ! Correlation between s and ln Nr (1st PDF comp.) ip  [-]
      corr_sNr_2_n,  & ! Correlation between s and ln Nr (2nd PDF comp.) ip  [-]
      corr_sNcn_1_n, & ! Correlation between s and ln Ncn (1st PDF comp.)    [-]
      corr_sNcn_2_n, & ! Correlation between s and ln Ncn (2nd PDF comp.)    [-]
      corr_trr_1_n,  & ! Correlation between t and ln rr (1st PDF comp.) ip  [-]
      corr_trr_2_n,  & ! Correlation between t and ln rr (2nd PDF comp.) ip  [-]
      corr_tNr_1_n,  & ! Correlation between t and ln Nr (1st PDF comp.) ip  [-]
      corr_tNr_2_n,  & ! Correlation between t and ln Nr (2nd PDF comp.) ip  [-]
      corr_tNcn_1_n, & ! Correlation between t and ln Ncn (1st PDF comp.)    [-]
      corr_tNcn_2_n, & ! Correlation between t and ln Ncn (2nd PDF comp.)    [-]
      corr_rrNr_1_n, & ! Correlation btwn. ln rr & ln Nr (1st PDF comp.) ip  [-]
      corr_rrNr_2_n    ! Correlation btwn. ln rr & ln Nr (2nd PDF comp.) ip  [-]

    real( kind = core_rknd ), intent(out) :: &
      rr1, & ! Mean rain water mixing ratio (1st PDF component)      [kg/kg]
      rr2, & ! Mean rain water mixing ratio (2nd PDF component)      [kg/kg]
      Nr1, & ! Mean rain drop concentration (1st PDF component)      [num/kg]
      Nr2    ! Mean rain drop concentration (2nd PDF component)      [num/kg]

    real( kind = core_rknd ), intent(out) :: &
      precip_frac,   & ! Precipitation fraction (overall)           [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component) [-]
      precip_frac_2    ! Precipitation fraction (2nd PDF component) [-]


    ! Unpack mu_x_i and sigma_x_i into Means and Standard Deviations.
    mu_w_1        = mu_x_1(iiLh_w)
    mu_w_2        = mu_x_2(iiLh_w)
    mu_s_1        = mu_x_1(iiLh_s_mellor)
    mu_s_2        = mu_x_2(iiLh_s_mellor)
    mu_t_1        = mu_x_1(iiLh_t_mellor)
    mu_t_2        = mu_x_2(iiLh_t_mellor)
    mu_rr_1_n     = mu_x_1(iiLh_rrain)
    mu_rr_2_n     = mu_x_2(iiLh_rrain)
    mu_Nr_1_n     = mu_x_1(iiLh_Nr)
    mu_Nr_2_n     = mu_x_2(iiLh_Nr)
    mu_Ncn_1_n    = mu_x_1(iiLh_Ncn)
    mu_Ncn_2_n    = mu_x_2(iiLh_Ncn)
    sigma_w_1     = sigma_x_1(iiLh_w)
    sigma_w_2     = sigma_x_2(iiLh_w)
    sigma_s_1     = sigma_x_1(iiLh_s_mellor)
    sigma_s_2     = sigma_x_2(iiLh_s_mellor)
    sigma_t_1     = sigma_x_1(iiLh_t_mellor)
    sigma_t_2     = sigma_x_2(iiLh_t_mellor)
    sigma_rr_1_n  = sigma_x_1(iiLh_rrain)
    sigma_rr_2_n  = sigma_x_2(iiLh_rrain)
    sigma_Nr_1_n  = sigma_x_1(iiLh_Nr)
    sigma_Nr_2_n  = sigma_x_2(iiLh_Nr)
    sigma_Ncn_1_n = sigma_x_1(iiLh_Ncn)
    sigma_Ncn_2_n = sigma_x_2(iiLh_Ncn)

    ! Unpack variables from hydromet_pdf_params
    mu_rr_1     = hydromet_pdf_params%mu_rr_1
    mu_rr_2     = hydromet_pdf_params%mu_rr_2
    mu_Nr_1     = hydromet_pdf_params%mu_Nr_1
    mu_Nr_2     = hydromet_pdf_params%mu_Nr_2
    mu_Ncn_1    = hydromet_pdf_params%mu_Ncn_1
    mu_Ncn_2    = hydromet_pdf_params%mu_Ncn_2
    sigma_rr_1  = hydromet_pdf_params%sigma_rr_1
    sigma_rr_2  = hydromet_pdf_params%sigma_rr_2
    sigma_Nr_1  = hydromet_pdf_params%sigma_Nr_1
    sigma_Nr_2  = hydromet_pdf_params%sigma_Nr_2
    sigma_Ncn_1 = hydromet_pdf_params%sigma_Ncn_1
    sigma_Ncn_2 = hydromet_pdf_params%sigma_Ncn_2

    rr1           = hydromet_pdf_params%rr1
    rr2           = hydromet_pdf_params%rr2
    Nr1           = hydromet_pdf_params%Nr1
    Nr2           = hydromet_pdf_params%Nr2
    precip_frac   = hydromet_pdf_params%precip_frac
    precip_frac_1 = hydromet_pdf_params%precip_frac_1
    precip_frac_2 = hydromet_pdf_params%precip_frac_2

    ! Unpack corr_array_1 into correlations (1st PDF component).
    corr_st_1     = corr_array_1(iiLH_t_mellor, iiLH_s_mellor)
    corr_ws_1     = corr_array_1(iiLH_s_mellor,iiLH_w)
    corr_srr_1_n  = corr_array_1(iiLH_rrain, iiLH_s_mellor)
    corr_sNr_1_n  = corr_array_1(iiLH_Nr, iiLH_s_mellor)
    corr_sNcn_1_n = corr_array_1(iiLH_Ncn, iiLH_s_mellor)
    corr_trr_1_n  = corr_array_1(iiLH_rrain, iiLH_t_mellor)
    corr_tNr_1_n  = corr_array_1(iiLH_Nr, iiLH_t_mellor)
    corr_tNcn_1_n = corr_array_1(iiLH_Ncn, iiLH_t_mellor)
    corr_wrr_1_n  = corr_array_1(iiLH_rrain, iiLH_w)
    corr_wNr_1_n  = corr_array_1(iiLH_Nr, iiLH_w)
    corr_wNcn_1_n = corr_array_1(iiLH_Ncn, iiLH_w)
    corr_rrNr_1_n = corr_array_1(iiLH_Nr, iiLH_rrain)

    ! Unpack corr_array_2 into correlations (2nd PDF component).
    corr_st_2     = corr_array_2(iiLH_t_mellor, iiLH_s_mellor)
    corr_ws_2     = corr_array_2(iiLH_s_mellor,iiLH_w)
    corr_srr_2_n  = corr_array_2(iiLH_rrain, iiLH_s_mellor)
    corr_sNr_2_n  = corr_array_2(iiLH_Nr, iiLH_s_mellor)
    corr_sNcn_2_n = corr_array_2(iiLH_Ncn, iiLH_s_mellor)
    corr_trr_2_n  = corr_array_2(iiLH_rrain, iiLH_t_mellor)
    corr_tNr_2_n  = corr_array_2(iiLH_Nr, iiLH_t_mellor)
    corr_tNcn_2_n = corr_array_2(iiLH_Ncn, iiLH_t_mellor)
    corr_wrr_2_n  = corr_array_2(iiLH_rrain, iiLH_w)
    corr_wNr_2_n  = corr_array_2(iiLH_Nr, iiLH_w)
    corr_wNcn_2_n = corr_array_2(iiLH_Ncn, iiLH_w)
    corr_rrNr_2_n = corr_array_2(iiLH_Nr, iiLH_rrain)


    return

  end subroutine unpack_pdf_params

!===============================================================================

end module setup_clubb_pdf_params
