! $Id$
!===============================================================================
module mixed_moment_PDF_integrals

  implicit none

  private

  public :: hydrometeor_mixed_moments

  private :: xphmp_integral_covar,         &
             hmxphmyp_integral_covar,      &
             xp_a_hmpb_integrals_all_MM,   &
             bivar_NL_x_hm_all_MM_comp_eq, &
             bivar_NL_int_PDF_comp_all_MM, &
             univar_N_int_PDF_comp_all_MM, &
             univar_L_int_PDF_comp_all_MM

  contains

  !=============================================================================
  subroutine hydrometeor_mixed_moments( nz, d_variables, hydromet, &
                                        mu_x_1_n, mu_x_2_n, &
                                        sigma_x_1_n, sigma_x_2_n, &
                                        corr_array_1_n, corr_array_2_n, &
                                        pdf_params, hydromet_pdf_params, &
                                        rtphmp_zt, thlphmp_zt, wp2hmp ) 

    ! Description:
    ! Calculates <rt'hm'>, <thl'hm'>, and <w'^2 hm'>, for all hydrometeors, hm.
    ! These terms are used in the liquid/ice water loading term as part of the
    ! buoyancy term in some of the CLUBB predictive equations.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        zt2zm    ! Procedure(s)

    use constants_clubb, only: &
        one,     & ! Constant(s)
        zero,    &
        w_tol,   &
        rt_tol,  &
        thl_tol

    use pdf_parameter_module, only: &
        pdf_parameter   ! Variable(s)

    use hydromet_pdf_parameter_module, only: &
        hydromet_pdf_parameter  ! Variable(s)

    use index_mapping, only: &
        hydromet2pdf_idx   ! Procedure(s)

    use pdf_utilities, only: &
        compute_mean_binormal, & ! Procedure(s)
        calc_corr_rt_x,        &
        calc_corr_thl_x

    use array_index, only: &
        hydromet_tol   ! Variable(s)

    use parameters_model, only: &
        hydromet_dim   ! Variable(s)

    use corr_varnce_module, only: &
        iiPDF_chi, & ! Variable(s)
        iiPDF_eta, &
        iiPDF_w

    use stats_type_utilities, only: &
        stat_update_var    ! Procedure(s)

    use stats_variables, only: &
        iwp2hmp,      & ! Variable(s)
        irtphmp,      &
        ithlphmp,     &
        ihmxphmyp,    &
        stats_zt,     &
        stats_zm,     &
        l_stats_samp

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz,          & ! Number of model vertical grid levels
      d_variables    ! Number of variables in the correlation array

    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(in) :: &
      hydromet    ! Mean of hydrometeor, hm (overall) (t-levels)    [units vary]

    real( kind = core_rknd ), dimension(d_variables, nz), intent(in) :: &
      mu_x_1_n,    & ! Mean array (normal space): PDF vars. (comp. 1) [un. vary]
      mu_x_2_n,    & ! Mean array (normal space): PDF vars. (comp. 2) [un. vary]
      sigma_x_1_n, & ! Std. dev. array (normal space): PDF vars (comp. 1) [u.v.]
      sigma_x_2_n    ! Std. dev. array (normal space): PDF vars (comp. 2) [u.v.]

    real( kind = core_rknd ), dimension(d_variables,d_variables,nz), &
    intent(in) :: &
      corr_array_1_n, & ! Corr. array (normal space) of PDF vars. (comp. 1)  [-]
      corr_array_2_n    ! Corr. array (normal space) of PDF vars. (comp. 2)  [-]

    type(pdf_parameter), dimension(nz), intent(in) :: &
      pdf_params    ! PDF parameters                                [units vary]

    type(hydromet_pdf_parameter), dimension(nz), intent(in) :: &
      hydromet_pdf_params    ! Hydrometeor PDF parameters           [units vary]

    ! Output Variables
    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(out) :: &
      wp2hmp,     & ! Higher-order mixed moment:  < w'^2 hm' > [(m/s)^2<hm un.>]
      rtphmp_zt,  & ! Covariance of rt and hm (on t-levs.)     [(kg/kg)<hm un.>]
      thlphmp_zt    ! Covariance of thl and hm (on t-levs.)    [K<hm units>]

    ! Local Variables
    ! Unpacked parameters.
    real( kind = core_rknd ) :: &
      mu_w_1,       & ! Mean of w (1st PDF component)                      [m/s]
      mu_w_2,       & ! Mean of w (2nd PDF component)                      [m/s]
      mu_rt_1,      & ! Mean of rt (1st PDF component)                   [kg/kg]
      mu_rt_2,      & ! Mean of rt (2nd PDF component)                   [kg/kg]
      mu_thl_1,     & ! Mean of thl (1st PDF component)                      [K]
      mu_thl_2,     & ! Mean of thl (2nd PDF component)                      [K]
      mu_hm_1,      & ! Mean of hm (1st PDF component) in-precip (ip) [hm units]
      mu_hm_2,      & ! Mean of hm (2nd PDF component) ip             [hm units]
      mu_hmy_1,     & ! Mean of hmy (1st PDF component) ip           [hmy units]
      mu_hmy_2,     & ! Mean of hmy (2nd PDF component) ip           [hmy units]
      mu_hm_1_n,    & ! Mean of ln hm (1st PDF component) ip      [ln(hm units)]
      mu_hm_2_n,    & ! Mean of ln hm (2nd PDF component) ip      [ln(hm units)]
      sigma_w_1,    & ! Standard deviation of w (1st PDF component)        [m/s]
      sigma_w_2,    & ! Standard deviation of w (2nd PDF component)        [m/s]
      sigma_rt_1,   & ! Standard deviation of rt (1st PDF component)     [kg/kg]
      sigma_rt_2,   & ! Standard deviation of rt (2nd PDF component)     [kg/kg]
      sigma_thl_1,  & ! Standard deviation of thl (1st PDF component)        [K]
      sigma_thl_2,  & ! Standard deviation of thl (2nd PDF component)        [K]
      sigma_chi_1,  & ! Standard deviation of chi (1st PDF component)    [kg/kg]
      sigma_chi_2,  & ! Standard deviation of chi (2nd PDF component)    [kg/kg]
      sigma_eta_1,  & ! Standard deviation of eta (1st PDF component)    [kg/kg]
      sigma_eta_2,  & ! Standard deviation of eta (2nd PDF component)    [kg/kg]
      sigma_hm_1,   & ! Standard deviation of hm (1st PDF comp.) ip   [hm units]
      sigma_hm_2,   & ! Standard deviation of hm (2nd PDF comp.) ip   [hm units]
      sigma_hmy_1,  & ! Standard deviation of hmy (1st PDF comp.) ip [hmy units]
      sigma_hmy_2,  & ! Standard deviation of hmy (2nd PDF comp.) ip [hmy units]
      sigma_hm_1_n, & ! Standard deviation of ln hm (1st PDF component) ip   [-]
      sigma_hm_2_n    ! Standard deviation of ln hm (2nd PDF component) ip   [-]

    real( kind = core_rknd ) :: &
      corr_chi_hm_1, & ! Correlation of chi and hm (1st PDF component) ip    [-]
      corr_chi_hm_2, & ! Correlation of chi and hm (2nd PDF component) ip    [-]
      corr_eta_hm_1, & ! Correlation of eta and hm (1st PDF component) ip    [-]
      corr_eta_hm_2, & ! Correlation of eta and hm (2nd PDF component) ip    [-]
      corr_hm_hmy_1, & ! Correlation of hm and hmy (1st PDF component) ip    [-]
      corr_hm_hmy_2, & ! Correlation of hm and hmy (2nd PDF component) ip    [-]
      corr_w_hm_1_n, & ! Correlation of w and ln hm (1st PDF component) ip   [-]
      corr_w_hm_2_n, & ! Correlation of w and ln hm (2nd PDF component) ip   [-]
      mixt_frac,     & ! Mixture fraction                                    [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component)          [-]
      precip_frac_2, & ! Precipitation fraction (2nd PDF component)          [-]
      crt_1,         & ! Coef. of rt in chi/eta eqns. (1st PDF component)    [-]
      crt_2,         & ! Coef. of rt in chi/eta eqns. (2nd PDF component)    [-]
      cthl_1,        & ! Coef. of thl: chi/eta eqns. (1st PDF comp.) [(kg/kg)/K]
      cthl_2,        & ! Coef. of thl: chi/eta eqns. (2nd PDF comp.) [(kg/kg)/K]
      hm_mean,       & ! Mean of hydrometeor, hm (overall)            [hm units]
      hmy_mean,      & ! Mean of second hydrometeor, hmy (overall)   [hmy units]
      hm_tol,        & ! Tolerance value of hydrometeor, hm           [hm units]
      hmy_tol          ! Tolerance value of second hydrometeor, hmy  [hmy units]

    ! Calculated or recalculated values.
    real( kind = core_rknd ) :: &
      corr_rt_hm_1,  & ! Correlation of rt and hm (1st PDF component)        [-]
      corr_rt_hm_2,  & ! Correlation of rt and hm (2nd PDF component)        [-]
      corr_thl_hm_1, & ! Correlation of thl and hm (1st PDF component)       [-]
      corr_thl_hm_2, & ! Correlation of thl and hm (2nd PDF component)       [-]
      wm,            & ! Mean of w (overall)                               [m/s]
      rtm,           & ! Mean of rt (overall)                            [kg/kg]
      thlm             ! Mean of thl (overall)                               [K]

    real( kind = core_rknd ), dimension(nz,hydromet_dim,hydromet_dim) :: &
      hmxphmyp_zt    ! Covariance (overall) of two hydrometeors  [hmx*hmy units]

    integer :: &
      hm_idx,  & ! Index of a hydrometeor in the hydrometeor set of indices
      hmy_idx, & ! Index of second hydrometeor in the hydrometeor set of indices
      pdf_idx, & ! Index of a hydrometeor in the PDF set of indices
      a_exp,   & ! Exponent on w' in < w'^a hm'^b >
      b_exp,   & ! Exponent on hm' in < w'^a hm'^b >
      k          ! Index of a vertical level


    ! Loop over all thermodynamic levels between the model lower and upper
    ! boundaries (thermodynamic levels 2 to gr%nz).
    do k = 2, nz, 1

       ! Unpack the means of w, rt, and thl in each PDF component.
       mu_w_1   = mu_x_1_n(iiPDF_w,k)
       mu_w_2   = mu_x_2_n(iiPDF_w,k)
       mu_rt_1  = pdf_params(k)%rt_1
       mu_rt_2  = pdf_params(k)%rt_2
       mu_thl_1 = pdf_params(k)%thl_1
       mu_thl_2 = pdf_params(k)%thl_2

       ! Unpack the standard deviations of w, rt, and thl in each PDF component.
       sigma_w_1   = sigma_x_1_n(iiPDF_w,k)
       sigma_w_2   = sigma_x_2_n(iiPDF_w,k)
       sigma_rt_1  = sqrt( pdf_params(k)%varnce_rt_1 )
       sigma_rt_2  = sqrt( pdf_params(k)%varnce_rt_2 )
       sigma_thl_1 = sqrt( pdf_params(k)%varnce_thl_1 )
       sigma_thl_2 = sqrt( pdf_params(k)%varnce_thl_2 )

       ! Unpack the standard deviations of chi and eta in each PDF component.
       sigma_chi_1 = sigma_x_1_n(iiPDF_chi,k)
       sigma_chi_2 = sigma_x_2_n(iiPDF_chi,k)
       sigma_eta_1 = sigma_x_1_n(iiPDF_eta,k)
       sigma_eta_2 = sigma_x_2_n(iiPDF_eta,k)

       ! Unpack the mixture fraction.
       mixt_frac = pdf_params(k)%mixt_frac

       ! Unpack the precipitation fraction in each PDF component.
       precip_frac_1 = hydromet_pdf_params(k)%precip_frac_1
       precip_frac_2 = hydromet_pdf_params(k)%precip_frac_2

       ! Unpack the coefficients of rt and thl in the chi/eta PDF transformation
       ! equations for each PDF component.
       crt_1  = pdf_params(k)%crt_1
       crt_2  = pdf_params(k)%crt_2
       cthl_1 = pdf_params(k)%cthl_1
       cthl_2 = pdf_params(k)%cthl_2

       ! Re-calculate rtm, thlm, and wm from PDF parameters.
       ! This needs to be done because rtm and thlm have been advanced since
       ! the PDF parameters have been calculated.  It is necessary to use values
       ! of the mean fields (rtm, thlm, and wm) that are consistent with the
       ! PDF parameters.  This does not need to be done for hydromet because
       ! hydrometeors have not been advanced since the hydrometeor PDF
       ! parameters were set up.
       wm   = compute_mean_binormal( mu_w_1, mu_w_2, mixt_frac )
       rtm  = compute_mean_binormal( mu_rt_1, mu_rt_2, mixt_frac )
       thlm = compute_mean_binormal( mu_thl_1, mu_thl_2, mixt_frac )


       ! Calculate <rt'hm'>, <thl'hm'>, and <w'^2 hm'> for each hydrometeor
       ! species.
       do hm_idx = 1, hydromet_dim, 1

          ! Unpack the mean (in-precip) of hm in each PDF component.
          mu_hm_1 = hydromet_pdf_params(k)%mu_hm_1(hm_idx)
          mu_hm_2 = hydromet_pdf_params(k)%mu_hm_2(hm_idx)

          ! Unpack the standard deviation (in-precip) of hm in each PDF
          ! component.
          sigma_hm_1 = hydromet_pdf_params(k)%sigma_hm_1(hm_idx)
          sigma_hm_2 = hydromet_pdf_params(k)%sigma_hm_2(hm_idx)

          ! Calculate the correlation (in-precip) of rt/thl and hm for each PDF
          ! component.  Since CLUBB uses a PDF transformation from rt and
          ! theta-l coordinates to chi and eta coordinates for each PDF
          ! component, the correlation arrays are written in terms of chi and
          ! eta correlations.  This makes a calculation necessary for these
          ! correlations.
          corr_chi_hm_1 = hydromet_pdf_params(k)%corr_chi_hm_1(hm_idx)
          corr_chi_hm_2 = hydromet_pdf_params(k)%corr_chi_hm_2(hm_idx)
          corr_eta_hm_1 = hydromet_pdf_params(k)%corr_eta_hm_1(hm_idx)
          corr_eta_hm_2 = hydromet_pdf_params(k)%corr_eta_hm_2(hm_idx)

          corr_rt_hm_1 &
          = calc_corr_rt_x( crt_1, sigma_rt_1, sigma_chi_1, &
                            sigma_eta_1, corr_chi_hm_1, corr_eta_hm_1 )

          corr_rt_hm_2 &
          = calc_corr_rt_x( crt_2, sigma_rt_2, sigma_chi_2, &
                            sigma_eta_2, corr_chi_hm_2, corr_eta_hm_2 )

          corr_thl_hm_1 &
          = calc_corr_thl_x( cthl_1, sigma_thl_1, sigma_chi_1, &
                             sigma_eta_1, corr_chi_hm_1, corr_eta_hm_1 )

          corr_thl_hm_2 &
          = calc_corr_thl_x( cthl_2, sigma_thl_2, sigma_chi_2, &
                             sigma_eta_2, corr_chi_hm_2, corr_eta_hm_2 )

          ! Unpack the tolerance value for the hydrometeor, hm.
          hm_tol = hydromet_tol(hm_idx)

          ! Calculate <rt'hm'>.
          rtphmp_zt(k,hm_idx) &
          = xphmp_integral_covar( mu_rt_1, mu_rt_2, mu_hm_1, mu_hm_2, &
                                  sigma_rt_1, sigma_rt_2, sigma_hm_1, &
                                  sigma_hm_2, corr_rt_hm_1, corr_rt_hm_2, &
                                  mixt_frac, precip_frac_1, &
                                  precip_frac_2, rtm, rt_tol, hm_tol )

          ! Calculate <thl'hm'>.
          thlphmp_zt(k,hm_idx) &
          = xphmp_integral_covar( mu_thl_1, mu_thl_2, mu_hm_1, mu_hm_2, &
                                  sigma_thl_1, sigma_thl_2, sigma_hm_1, &
                                  sigma_hm_2, corr_thl_hm_1, corr_thl_hm_2, &
                                  mixt_frac, precip_frac_1, &
                                  precip_frac_2, thlm, thl_tol, hm_tol )

          ! Find the index of hydrometeor in the PDF indices.
          pdf_idx = hydromet2pdf_idx(hm_idx)

          ! Unpack the mean (in-precip) of ln hm in each PDF component.
          mu_hm_1_n = mu_x_1_n(pdf_idx,k)
          mu_hm_2_n = mu_x_2_n(pdf_idx,k)

          ! Unpack the standard deviation (in-precip) of ln hm in each PDF
          ! component.
          sigma_hm_1_n = sigma_x_1_n(pdf_idx,k)
          sigma_hm_2_n = sigma_x_2_n(pdf_idx,k)

          ! Unpack the correlation (in-precip) of w and ln hm in each PDF
          ! component.
          corr_w_hm_1_n = corr_array_1_n(pdf_idx,iiPDF_w,k)
          corr_w_hm_2_n = corr_array_2_n(pdf_idx,iiPDF_w,k)

          ! Unpack the mean (overall) value of the hydrometeor.
          hm_mean = hydromet(k,hm_idx)

          ! The general form of the mixed moment equation is <w'^a hm'^b>.
          ! For <w'^2 hm'>, a = 2 and b = 1.
          a_exp = 2
          b_exp = 1

          ! Calculate <w'^2 hm'>.
          wp2hmp(k,hm_idx) &
          = xp_a_hmpb_integrals_all_MM( mu_w_1, mu_w_2, mu_hm_1, mu_hm_2, &
                                        mu_hm_1_n, mu_hm_2_n, sigma_w_1, &
                                        sigma_w_2, sigma_hm_1, sigma_hm_2, &
                                        sigma_hm_1_n, sigma_hm_2_n, &
                                        corr_w_hm_1_n, corr_w_hm_2_n, &
                                        mixt_frac, precip_frac_1, &
                                        precip_frac_2, wm, hm_mean, &
                                        w_tol, hm_tol, a_exp, b_exp )

          ! Calculate the covariance (overall) of two hydrometeors, <hmx'hmy'>,
          ! for each unique set of two different hydrometeors.
          do hmy_idx = hm_idx+1, hydromet_dim, 1

             ! Unpack the mean (in-precip) of the second hydrometeor, hmy, in
             ! each PDF component.
             mu_hmy_1 = hydromet_pdf_params(k)%mu_hm_1(hmy_idx)
             mu_hmy_2 = hydromet_pdf_params(k)%mu_hm_2(hmy_idx)

             ! Unpack the standard deviation (in-precip) of hmy in each PDF
             ! component.
             sigma_hmy_1 = hydromet_pdf_params(k)%sigma_hm_1(hmy_idx)
             sigma_hmy_2 = hydromet_pdf_params(k)%sigma_hm_2(hmy_idx)

             ! Unpack the correlation (in-precip) of hm and hmy in each PDF
             ! component.
             corr_hm_hmy_1 &
             = hydromet_pdf_params(k)%corr_hmx_hmy_1(hm_idx,hmy_idx)
             corr_hm_hmy_2 &
             = hydromet_pdf_params(k)%corr_hmx_hmy_2(hm_idx,hmy_idx)

             ! Unpack the mean (overall) value of hmy.
             hmy_mean = hydromet(k,hmy_idx)

             ! Unpack the tolerance value for the second hydrometeor, hmy.
             hmy_tol = hydromet_tol(hmy_idx)

             ! Calculate the covariance <hmx'hmy'>.
             hmxphmyp_zt(k,hmy_idx,hm_idx) &
             = hmxphmyp_integral_covar( mu_hm_1, mu_hm_2, mu_hmy_1, mu_hmy_2, &
                                        sigma_hm_1, sigma_hm_2, sigma_hmy_1, &
                                        sigma_hmy_2, corr_hm_hmy_1, &
                                        corr_hm_hmy_2, mixt_frac, &
                                        precip_frac_1, precip_frac_2, hm_mean, &
                                        hmy_mean, hm_tol, hmy_tol )

          enddo ! hmy_idx = hm_idx+1, hydromet_dim, 1

       enddo ! hm_idx = 1, hydromet_dim, 1

    enddo ! k = 2, nz, 1

    ! Lower boundary
    ! The k = 1 thermodynamic is below the model lower boundary.
    ! Set all moments to 0.
    rtphmp_zt(1,:)     = zero
    thlphmp_zt(1,:)    = zero
    wp2hmp(1,:)        = zero
    hmxphmyp_zt(1,:,:) = zero


    ! Statistics
    if ( l_stats_samp ) then

       do hm_idx = 1, hydromet_dim, 1

          if ( iwp2hmp(hm_idx) > 0 ) then
             call stat_update_var( iwp2hmp(hm_idx), wp2hmp(:,hm_idx), stats_zt )
          endif ! iwp2hmp(hm_idx) > 0

          if ( irtphmp(hm_idx) > 0 ) then
             call stat_update_var( irtphmp(hm_idx), &
                                   zt2zm( rtphmp_zt(:,hm_idx) ), stats_zm )
          endif ! irtphmp(hm_idx) > 0

          if ( ithlphmp(hm_idx) > 0 ) then
             call stat_update_var( ithlphmp(hm_idx), &
                                   zt2zm( thlphmp_zt(:,hm_idx) ), stats_zm )
          endif ! ithlphmp(hm_idx) > 0

          do hmy_idx = hm_idx+1, hydromet_dim, 1
             if ( ihmxphmyp(hmy_idx,hm_idx) > 0 ) then
                call stat_update_var( ihmxphmyp(hmy_idx,hm_idx), &
                                      zt2zm( hmxphmyp_zt(:,hmy_idx,hm_idx) ), &
                                      stats_zm )
             endif ! ihmxphmyp(hmy_idx,hm_idx) > 0
          enddo ! hmy_idx = hm_idx+1, hydromet_dim, 1

       enddo ! hm_idx = 1, hydromet_dim, 1

    endif ! l_stats_samp


    return

  end subroutine hydrometeor_mixed_moments

  !=============================================================================
  function xphmp_integral_covar( mu_x_1, mu_x_2, mu_hm_1, mu_hm_2, &
                                 sigma_x_1, sigma_x_2, sigma_hm_1, &
                                 sigma_hm_2, corr_x_hm_1, corr_x_hm_2, &
                                 mixt_frac, precip_frac_1, &
                                 precip_frac_2, x_mean, x_tol, hm_tol )

    ! Description:
    ! Solves for the covariance < x'hm' >.  The variable "x" is a variable that
    ! is distributed binormally (meaning that it has a normally-distributed
    ! individual marginal in each PDF component).  This applies to such
    ! variables as w, rt, thl, sclr, etc.  The variable "hm" stands for any
    ! precipitating hydrometeor.  This applies to such variables as rr, Nr, ri,
    ! Ni, etc.
    !
    ! The covariance < x'hm' > can also be found by passing a = 1 and b = 1 into
    ! function xp_a_hmpb_integrals_all_MM, which is found below.  However,
    ! the code found here has been streamlined for the special case of
    ! covariances.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        one  ! Constant(s)

    use clubb_precision, only: &
        core_rknd  ! Variable(s) 

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_x_1,        & ! Mean of x (1st PDF component)                 [x units]
      mu_x_2,        & ! Mean of x (2nd PDF component)                 [x units]
      mu_hm_1,       & ! Mean of hm (1st PDF component)               [hm units]
      mu_hm_2,       & ! Mean of hm (2nd PDF component)               [hm units]
      sigma_x_1,     & ! Standard deviation of x (1st PDF component)   [x units]
      sigma_x_2,     & ! Standard deviation of x (2nd PDF component)   [x units]
      sigma_hm_1,    & ! Standard deviation of hm (1st PDF component) [hm units]
      sigma_hm_2,    & ! Standard deviation of hm (2nd PDF component) [hm units]
      corr_x_hm_1,   & ! Correlation of x and hm (1st PDF component)         [-]
      corr_x_hm_2,   & ! Correlation of x and hm (2nd PDF component)         [-]
      mixt_frac,     & ! Mixture fraction                                    [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component)          [-]
      precip_frac_2, & ! Precipitation fraction (2nd PDF component)          [-]
      x_mean,        & ! Mean of x (overall)                           [x units]
      x_tol,         & ! Tolerance value of x                          [x units]
      hm_tol           ! Tolerance value for hm                       [hm units]

    ! Return Variable
    real( kind = core_rknd ) :: &
      xphmp_integral_covar


    if ( ( sigma_x_1 <= x_tol .or. sigma_hm_1 <= hm_tol ) .and. &
         ( sigma_x_2 <= x_tol .or. sigma_hm_2 <= hm_tol ) ) then

       xphmp_integral_covar &
       = mixt_frac * precip_frac_1 * ( mu_x_1 - x_mean ) * mu_hm_1 &
         + ( one - mixt_frac ) * precip_frac_2 * ( mu_x_2 - x_mean ) * mu_hm_2


    elseif ( sigma_x_1 <= x_tol .or. sigma_hm_1 <= hm_tol ) then

       xphmp_integral_covar &
       = mixt_frac * precip_frac_1 &
         * ( mu_x_1 - x_mean ) * mu_hm_1 &
         + ( one - mixt_frac ) * precip_frac_2 &
           * ( ( mu_x_2 - x_mean ) * mu_hm_2 &
               + corr_x_hm_2 * sigma_x_2 * sigma_hm_2 )


    elseif ( sigma_x_2 <= x_tol .or. sigma_hm_2 <= hm_tol ) then

       xphmp_integral_covar &
       = mixt_frac * precip_frac_1 &
         * ( ( mu_x_1 - x_mean ) * mu_hm_1 &
             + corr_x_hm_1 * sigma_x_1 * sigma_hm_1 ) &
         + ( one - mixt_frac ) * precip_frac_2 &
           * ( mu_x_2 - x_mean ) * mu_hm_2


    else ! sigma_x_1 > x_tol and sigma_hm_1 > hm_tol
         ! and sigma_x_2 > x_tol and sigma_hm_2 > hm_tol

       xphmp_integral_covar &
       = mixt_frac * precip_frac_1 &
         * ( ( mu_x_1 - x_mean ) * mu_hm_1 &
             + corr_x_hm_1 * sigma_x_1 * sigma_hm_1 ) &
         + ( one - mixt_frac ) * precip_frac_2 &
           * ( ( mu_x_2 - x_mean ) * mu_hm_2 &
               + corr_x_hm_2 * sigma_x_2 * sigma_hm_2 )


    endif


    return

  end function xphmp_integral_covar

  !=============================================================================
  function hmxphmyp_integral_covar( mu_hmx_1, mu_hmx_2, mu_hmy_1, mu_hmy_2, &
                                    sigma_hmx_1, sigma_hmx_2, sigma_hmy_1, &
                                    sigma_hmy_2, corr_hmx_hmy_1, &
                                    corr_hmx_hmy_2, mixt_frac, &
                                    precip_frac_1, precip_frac_2, hmx_mean, &
                                    hmy_mean, hmx_tol, hmy_tol )

    ! Description:
    ! Solves for the covariance of two precipitating hydrometoers < hmx'hmy' >.
    ! The variable "hmx" stands for any precipitating hydrometeor.  This applies
    ! to such variables as rr, Nr, ri, Ni, etc.  The variable "hmy" stands for
    ! any different precipitating hydrometeor.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        one  ! Constant(s)

    use clubb_precision, only: &
        core_rknd  ! Variable(s) 

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_hmx_1,       & ! Mean of hmx (1st PDF component)            [hmx units]
      mu_hmx_2,       & ! Mean of hmx (2nd PDF component)            [hmx units]
      mu_hmy_1,       & ! Mean of hmy (1st PDF component)            [hmy units]
      mu_hmy_2,       & ! Mean of hmy (2nd PDF component)            [hmy units]
      sigma_hmx_1,    & ! Standard deviation of hmx (1st PDF comp.)  [hmx units]
      sigma_hmx_2,    & ! Standard deviation of hmx (2nd PDF comp.)  [hmx units]
      sigma_hmy_1,    & ! Standard deviation of hmy (1st PDF comp.)  [hmy units]
      sigma_hmy_2,    & ! Standard deviation of hmy (2nd PDF comp.)  [hmy units]
      corr_hmx_hmy_1, & ! Correlation of hmx and hmy (1st PDF component)     [-]
      corr_hmx_hmy_2, & ! Correlation of hmx and hmy (2nd PDF component)     [-]
      mixt_frac,      & ! Mixture fraction                                   [-]
      precip_frac_1,  & ! Precipitation fraction (1st PDF component)         [-]
      precip_frac_2,  & ! Precipitation fraction (2nd PDF component)         [-]
      hmx_mean,       & ! Mean of hmx (overall)                      [hmx units]
      hmy_mean,       & ! Mean of hmy (overall)                      [hmy units]
      hmx_tol,        & ! Tolerance value for hmx                    [hmx units]
      hmy_tol           ! Tolerance value for hmy                    [hmy units]

    ! Return Variable
    real( kind = core_rknd ) :: &
      hmxphmyp_integral_covar    ! Covariance of hmx and hmy     [hmx*hmy units]


    if ( ( sigma_hmx_1 <= hmx_tol .or. sigma_hmy_1 <= hmy_tol ) .and. &
         ( sigma_hmx_2 <= hmx_tol .or. sigma_hmy_2 <= hmy_tol ) ) then

       hmxphmyp_integral_covar &
       = mixt_frac * precip_frac_1 * mu_hmx_1 * mu_hmy_1 &
         + ( one - mixt_frac ) * precip_frac_2 * mu_hmx_2 * mu_hmy_2 &
         - hmx_mean * hmy_mean


    elseif ( sigma_hmx_1 <= hmx_tol .or. sigma_hmy_1 <= hmy_tol ) then

       hmxphmyp_integral_covar &
       = mixt_frac * precip_frac_1 * mu_hmx_1 * mu_hmy_1 &
         + ( one - mixt_frac ) * precip_frac_2 &
           * ( mu_hmx_2 * mu_hmy_2 &
               + corr_hmx_hmy_2 * sigma_hmx_2 * sigma_hmy_2 ) &
         - hmx_mean * hmy_mean


    elseif ( sigma_hmx_2 <= hmx_tol .or. sigma_hmy_2 <= hmy_tol ) then

       hmxphmyp_integral_covar &
       = mixt_frac * precip_frac_1 &
         * ( mu_hmx_1 * mu_hmy_1 &
             + corr_hmx_hmy_1 * sigma_hmx_1 * sigma_hmy_1 ) &
         + ( one - mixt_frac ) * precip_frac_2 * mu_hmx_2 * mu_hmy_2 &
         - hmx_mean * hmy_mean


    else ! sigma_hmx_1 > hmx_tol and sigma_hmy_1 > hmy_tol
         ! and sigma_hmx_2 > hmx_tol and sigma_hmy_2 > hmy_tol

       hmxphmyp_integral_covar &
       = mixt_frac * precip_frac_1 &
         * ( mu_hmx_1 * mu_hmy_1 &
             + corr_hmx_hmy_1 * sigma_hmx_1 * sigma_hmy_1 ) &
         + ( one - mixt_frac ) * precip_frac_2 &
           * ( mu_hmx_2 * mu_hmy_2 &
               + corr_hmx_hmy_2 * sigma_hmx_2 * sigma_hmy_2 ) &
         - hmx_mean * hmy_mean


    endif


    return

  end function hmxphmyp_integral_covar

  !=============================================================================
  function xp_a_hmpb_integrals_all_MM( mu_x_1, mu_x_2, mu_hm_1, mu_hm_2, &
                                       mu_hm_1_n, mu_hm_2_n, sigma_x_1, &
                                       sigma_x_2, sigma_hm_1, sigma_hm_2, &
                                       sigma_hm_1_n, sigma_hm_2_n, &
                                       corr_x_hm_1_n, corr_x_hm_2_n, &
                                       mixt_frac, precip_frac_1, &
                                       precip_frac_2, x_mean, hm_mean, &
                                       x_tol, hm_tol, a_exp, b_exp )

    ! Description:
    ! Solves for any mixed moment < x'^a hm'^b >, where a and b are both
    ! integers with values greater than or equal to 0.  The variable "x" is a
    ! variable that is distributed binormally (meaning that it has a
    ! normally-distributed individual marginal in each PDF component).  This
    ! applies to such variables as w, rt, thl, sclr, etc.  The variable "hm"
    ! stands for any precipitating hydrometeor.  This applies to such variables
    ! as rr, Nr, ri, Ni, etc.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        one  ! Constant(s)

    use clubb_precision, only: &
        core_rknd  ! Variable(s) 

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_x_1,        & ! Mean of x (1st PDF component)                 [x units]
      mu_x_2,        & ! Mean of x (2nd PDF component)                 [x units]
      mu_hm_1,       & ! Mean of hm (1st PDF component)               [hm units]
      mu_hm_2,       & ! Mean of hm (2nd PDF component)               [hm units]
      mu_hm_1_n,     & ! Mean of ln hm (1st PDF component)        [ln(hm units)]
      mu_hm_2_n,     & ! Mean of ln hm (2nd PDF component)        [ln(hm units)]
      sigma_x_1,     & ! Standard deviation of x (1st PDF component)   [x units]
      sigma_x_2,     & ! Standard deviation of x (2nd PDF component)   [x units]
      sigma_hm_1,    & ! Standard deviation of hm (1st PDF component) [hm units]
      sigma_hm_2,    & ! Standard deviation of hm (2nd PDF component) [hm units]
      sigma_hm_1_n,  & ! Standard deviation of ln hm (1st PDF component)     [-]
      sigma_hm_2_n,  & ! Standard deviation of ln hm (2nd PDF component)     [-]
      corr_x_hm_1_n, & ! Correlation of x and ln hm (1st PDF component)      [-]
      corr_x_hm_2_n, & ! Correlation of x and ln hm (2nd PDF component)      [-]
      mixt_frac,     & ! Mixture fraction                                    [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component)          [-]
      precip_frac_2, & ! Precipitation fraction (2nd PDF component)          [-]
      x_mean,        & ! Mean of x (overall)                           [x units]
      hm_mean,       & ! Mean of hm (overall)                         [hm units]
      x_tol,         & ! Tolerance value of x                          [x units]
      hm_tol           ! Tolerance value for hm                       [hm units]

    integer, intent(in) :: &
      a_exp, & ! Order prime of x1 - < x1 >                                  [-]
      b_exp    ! Order prime of x2 - < x2 >                                  [-]

    ! Return Variable
    real( kind = core_rknd ) ::  &
      xp_a_hmpb_integrals_all_MM


    xp_a_hmpb_integrals_all_MM &
    = mixt_frac &
      * bivar_NL_x_hm_all_MM_comp_eq( mu_x_1, mu_hm_1, mu_hm_1_n, &
                                      sigma_x_1, sigma_hm_1, sigma_hm_1_n, &
                                      corr_x_hm_1_n, precip_frac_1, x_mean, &
                                      hm_mean, x_tol, hm_tol, a_exp, b_exp ) &
      + ( one - mixt_frac ) &
        * bivar_NL_x_hm_all_MM_comp_eq( mu_x_2, mu_hm_2, mu_hm_2_n, &
                                        sigma_x_2, sigma_hm_2, sigma_hm_2_n, &
                                        corr_x_hm_2_n, precip_frac_2, x_mean, &
                                        hm_mean, x_tol, hm_tol, a_exp, b_exp )


    return

  end function xp_a_hmpb_integrals_all_MM

  !=============================================================================
  function bivar_NL_x_hm_all_MM_comp_eq( mu_x_i, mu_hm_i, mu_hm_i_n, &
                                         sigma_x_i, sigma_hm_i, sigma_hm_i_n, &
                                         corr_x_hm_i_n, precip_frac_i, x_mean, &
                                         hm_mean, x_tol, hm_tol, a_exp, b_exp )

    ! Description:
    ! Solves the portion of the integral for < x'^a hm'^b > that relates to the
    ! ith PDF component.  This takes into account the portions of the component
    ! that are within-precipitation and outside-precipitation.  The equation
    ! and function that are used depends on whether or not x and/or hm vary
    ! inside of a PDF component.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        one  ! Constant(s)

    use clubb_precision, only: &
        core_rknd  ! Variable(s) 

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_x_i,        & ! Mean of x (ith PDF component)                 [x units]
      mu_hm_i,       & ! Mean of hm (ith PDF component)               [hm units]
      mu_hm_i_n,     & ! Mean of ln hm (ith PDF component)        [ln(hm units)]
      sigma_x_i,     & ! Standard deviation of x (ith PDF component)   [x units]
      sigma_hm_i,    & ! Standard deviation of hm (ith PDF component) [hm units]
      sigma_hm_i_n,  & ! Standard deviation of ln hm (ith PDF component)     [-]
      corr_x_hm_i_n, & ! Correlation of x and ln hm (ith PDF component)      [-]
      precip_frac_i, & ! Precipitation fraction (ith PDF component)          [-]
      x_mean,        & ! Mean of x (overall)                           [x units]
      hm_mean,       & ! Mean of hm (overall)                         [hm units]
      x_tol,         & ! Tolerance value of x                          [x units]
      hm_tol           ! Tolerance value for hm                       [hm units]

    integer, intent(in) :: &
      a_exp, & ! Order prime of x1 - < x1 >                                  [-]
      b_exp    ! Order prime of x2 - < x2 >                                  [-]

    ! Return Variable
    real( kind = core_rknd ) ::  &
      bivar_NL_x_hm_all_MM_comp_eq


    if ( sigma_x_i <= x_tol .and. sigma_hm_i <= hm_tol ) then

       ! Both x and hm do not vary in the ith PDF component.

       bivar_NL_x_hm_all_MM_comp_eq &
       = ( mu_x_i - x_mean )**a_exp &
         * ( precip_frac_i * ( mu_hm_i - hm_mean )**b_exp &
             + ( one - precip_frac_i ) * ( -hm_mean )**b_exp )


    elseif ( sigma_x_i <= x_tol ) then

       ! Only x does not vary in the ith PDF component.

       bivar_NL_x_hm_all_MM_comp_eq &
       = ( mu_x_i - x_mean )**a_exp &
         * ( precip_frac_i &
             * univar_L_int_PDF_comp_all_MM( mu_hm_i_n, sigma_hm_i_n, hm_mean, &
                                             b_exp ) &
             + ( one - precip_frac_i ) * ( -hm_mean )**b_exp )     


    elseif ( sigma_hm_i <= hm_tol ) then

       ! Only hm does not vary in the ith PDF component.

       bivar_NL_x_hm_all_MM_comp_eq &
       = ( precip_frac_i * ( mu_hm_i - hm_mean )**b_exp &
           + ( one - precip_frac_i ) * ( -hm_mean )**b_exp ) &
         * univar_N_int_PDF_comp_all_MM( mu_x_i, sigma_x_i, x_mean, a_exp )


    else ! sigma_x_i > x_tol and sigma_hm_i > hm_tol

       ! Both x and hm vary in the ith PDF component.

       bivar_NL_x_hm_all_MM_comp_eq &
       = precip_frac_i &
         * bivar_NL_int_PDF_comp_all_MM( mu_x_i, mu_hm_i_n, sigma_x_i, &
                                         sigma_hm_i_n, corr_x_hm_i_n, &
                                         x_mean, hm_mean, a_exp, b_exp ) &
         + ( one - precip_frac_i ) &
           * ( -hm_mean )**b_exp &
           * univar_N_int_PDF_comp_all_MM( mu_x_i, sigma_x_i, x_mean, a_exp )

 
    endif


    return

  end function bivar_NL_x_hm_all_MM_comp_eq

  !=============================================================================
  function bivar_NL_int_PDF_comp_all_MM( mu_x1_i, mu_x2_i_n, sigma_x1_i, &
                                         sigma_x2_i_n, corr_x1_x2_i_n, &
                                         x1_mean, x2_mean, a_exp, b_exp )

    ! Description:
    ! This function is the evaluated form of the following integral:
    !
    ! INT(-inf:inf) INT(0:inf) ( x1 - <x1> )^a ( x2 - <x2> )^b
    !                          * P_NL_i( x1, x2 ) dx2 dx1;
    !
    ! where P_NL_i( x1, x2 ) is the functional form of a bivariate
    ! normal-lognormal PDF (in the PDF component), x1 is a variable that has an
    ! individual marginal that is distributed normally (in the PDF component),
    ! and x2 is a variable that has an individual marginal that is distributed
    ! lognormally (in the PDF component).  Additionally, <x1> is the overall
    ! mean of x1, <x2> is the overall mean of x2, "a" is the integer (>= 0)
    ! order of the mixed moment with respect to x1, and b is the integer (>= 0)
    ! order of the mixed moment with respect to x2.
    !
    ! When the integral is evaluated, the equation is:
    !
    ! INT(-inf:inf) INT(0:inf) ( x1 - <x1> )^a ( x2 - <x2> )^b
    !                          * P_NL_i( x1, x2 ) dx2 dx1
    ! = SUM( p = 0:floor(a/2) ) SUM ( q = 0:b )
    !   ( a! / ( ( a - 2p )! p! ) ) * ( b! / ( ( b - q )! q! ) )
    !   * [ (1/2) * sigma_x1_i**2 ]**p
    !   * ( mu_x1_i - <x1>
    !       + corr_x1_x2_i_n * sigma_x1_i * sigma_x2_i_n * q )**(a-2p)
    !   * ( -<x2> )**(b-q)
    !   * exp{ mu_x2_i_n * q + (1/2) * sigma_x2_i_n**2 * q**2 }.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        one_half, & ! Constant(s)
        zero

    use KK_utilities, only:  &
        factorial  ! Procedure(s)

    use clubb_precision, only: &
        core_rknd  ! Variable(s) 

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_x1_i,        & ! Mean of x1 (ith PDF component)              [x1 units]
      mu_x2_i_n,      & ! Mean of ln x2 (ith PDF component)       [ln(x2 units)]
      sigma_x1_i,     & ! Standard deviation of x1 (ith PDF comp.)    [x1 units]
      sigma_x2_i_n,   & ! Standard deviation of ln x2 (ith PDF component)    [-]
      corr_x1_x2_i_n, & ! Correlation of x1 and ln x2 (ith PDF component)    [-]
      x1_mean,        & ! Mean of x1 (overall)                        [x1 units]
      x2_mean           ! Mean of x2 (overall)                        [x2 units]
    
    integer, intent(in) :: &
      a_exp, & ! Order prime of x1 - < x1 >                                  [-]
      b_exp    ! Order prime of x2 - < x2 >                                  [-]

    ! Return Variable
    real( kind = core_rknd ) ::  &
      bivar_NL_int_PDF_comp_all_MM

    ! Local Variables
    integer :: &
      p, & ! Index relating to exponent a (a_exp)
      q    ! Index relating to exponent b (b_exp)

    real( kind = core_rknd ) ::  &
      sigma_sum,      & ! Running sum total value
      q_fp,           & ! Floating point version of index q
      factorial_term    ! Floating point version of factorial term


    ! Initialize sigma_sum
    sigma_sum = zero

    do p = 0, a_exp/2, 1
       do q = 0, b_exp, 1

          q_fp = real( q, kind = core_rknd )

          factorial_term &
          = real( ( factorial( a_exp ) &
                    / ( factorial( a_exp - 2*p ) * factorial( p ) ) ) &
                  * ( factorial( b_exp ) &
                      / ( factorial( b_exp - q ) * factorial( q ) ) ), &
                  kind = core_rknd )
 
          sigma_sum &
          = sigma_sum &
            + factorial_term &
              * ( one_half * sigma_x1_i**2 )**p &
              * ( mu_x1_i - x1_mean &
                  + corr_x1_x2_i_n * sigma_x1_i &
                                   * sigma_x2_i_n * q_fp )**(a_exp-2*p) &
              * ( -x2_mean )**(b_exp-q) &
              * exp( mu_x2_i_n * q_fp + one_half * sigma_x2_i_n**2 * q_fp**2 )

       enddo
    enddo 


    bivar_NL_int_PDF_comp_all_MM = sigma_sum


    return

  end function bivar_NL_int_PDF_comp_all_MM

  !=============================================================================
  function univar_N_int_PDF_comp_all_MM( mu_x_i, sigma_x_i, x_mean, a_exp )

    ! Description:
    ! This function is the evaluated form of the following integral:
    !
    ! INT(-inf:inf) ( x - <x> )^a * P_N_i( x ) dx;
    !
    ! where P_N_i( x ) is the functional form of a single-variable normal PDF
    ! (in the PDF component), and x is a variable that has an individual
    ! marginal that is distributed normally (in the PDF component).
    ! Additionally, <x> is the overall mean of x, and "a" is the integer (>= 0)
    ! order of the mixed moment with respect to x.
    !
    ! When the integral is evaluated, the equation is:
    !
    ! INT(-inf:inf) ( x - <x> )^a * P_N_i( x ) dx
    ! = SUM( p = 0:floor(a/2) ) ( a! / ( ( a - 2p )! p! ) )
    !   * [ (1/2) * sigma_x_i**2 ]**p * ( mu_x_i - <x> )**(a-2p).

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        one_half, & ! Constant(s)
        zero

    use KK_utilities, only:  &
        factorial  ! Procedure(s)

    use clubb_precision, only: &
        core_rknd  ! Variable(s) 

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_x_i,    & ! Mean of x (ith PDF component)                  [x units]
      sigma_x_i, & ! Standard deviation of x (ith PDF component)    [x units]
      x_mean       ! Mean of x (overall)                            [x units]
    
    integer, intent(in) :: &
      a_exp    ! Order prime of x - < x >                           [-]

    ! Return Variable
    real( kind = core_rknd ) ::  &
      univar_N_int_PDF_comp_all_MM

    ! Local Variables
    integer :: &
      p    ! Index relating to exponent a (a_exp)

    real( kind = core_rknd ) ::  &
      sigma_sum,      & ! Running sum total value
      factorial_term    ! Floating point version of factorial term


    ! Initialize sigma_sum
    sigma_sum = zero

    do p = 0, a_exp/2, 1

       factorial_term &
       = real( factorial( a_exp ) &
               / ( factorial( a_exp - 2*p ) * factorial( p ) ), &
               kind = core_rknd )
 
       sigma_sum &
       = sigma_sum &
         + factorial_term &
           * ( one_half * sigma_x_i**2 )**p  &
           * ( mu_x_i - x_mean )**(a_exp-2*p)

    enddo 


    ! Total integral within a component.
    univar_N_int_PDF_comp_all_MM = sigma_sum


    return

  end function univar_N_int_PDF_comp_all_MM

  !=============================================================================
  function univar_L_int_PDF_comp_all_MM( mu_x_i_n, sigma_x_i_n, x_mean, b_exp )

    ! Description:
    ! This function is the evaluated form of the following integral:
    !
    ! INT(0:inf) ( x - <x> )^b * P_L_i( x ) dx;
    !
    ! where P_L_i( x ) is the functional form of a single-variable lognormal PDF
    ! (in the PDF component), and x is a variable that has an individual
    ! marginal that is distributed lognormally (in the PDF component).
    ! Additionally, <x> is the overall mean of x, and b is the integer (>= 0)
    ! order of the mixed moment with respect to x.
    !
    ! When the integral is evaluated, the equation is:
    !
    ! INT(0:inf) ( x - <x> )^b * P_L_i( x ) dx;
    ! = SUM ( q = 0:b ) ( b! / ( ( b - q )! q! ) )
    !   * ( -<x> )**(b-q)
    !   * exp{ mu_x_i_n * q + (1/2) * sigma_x_i_n**2 * q**2 }.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        one_half, & ! Constant(s)
        zero

    use KK_utilities, only:  &
        factorial  ! Procedure(s)

    use clubb_precision, only: &
        core_rknd  ! Variable(s) 

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_x_i_n,    & ! Mean of ln x (ith PDF component)            [ln(x units)]
      sigma_x_i_n, & ! Standard deviation of ln x (ith PDF component)        [-]
      x_mean         ! Mean of x (overall)                             [x units]
    
    integer, intent(in) :: &
      b_exp    ! Order prime of x - < x >                                    [-]

    ! Return Variable
    real( kind = core_rknd ) ::  &
      univar_L_int_PDF_comp_all_MM

    ! Local Variables
    integer :: &
      q    ! Index relating to exponent b (b_exp)

    real( kind = core_rknd ) ::  &
      sigma_sum,      & ! Running sum total value
      q_fp,           & ! Floating point version of index q
      factorial_term    ! Floating point version of factorial term


    ! Initialize sigma_sum
    sigma_sum = zero

    do q = 0, b_exp, 1

       q_fp = real( q, kind = core_rknd )

       factorial_term &
       = real( factorial( b_exp ) &
               / ( factorial( b_exp - q ) * factorial( q ) ), &
               kind = core_rknd )
 
       sigma_sum &
       = sigma_sum &
         + factorial_term &
           * ( -x_mean )**(b_exp-q) &
           * exp( mu_x_i_n * q_fp + one_half * sigma_x_i_n**2 * q_fp**2 )

    enddo


    univar_L_int_PDF_comp_all_MM = sigma_sum


    return

  end function univar_L_int_PDF_comp_all_MM

!===============================================================================

end module mixed_moment_PDF_integrals
