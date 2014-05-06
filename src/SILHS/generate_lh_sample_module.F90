!-------------------------------------------------------------------------------
!$Id$
!===============================================================================
module generate_lh_sample_module

  use clubb_precision, only: &
    dp ! double precision

  implicit none

  public :: generate_lh_sample, generate_uniform_sample, &
     ltqnorm, multiply_Cholesky, generate_lh_sample_mod

  private :: sample_points, gaus_mixt_points, & 
    truncate_gaus_mixt, &
    st_2_rtthl, log_sqd_normalized, choose_permuted_random, &
    set_min_varnce_and_mean, construct_gaus_LN_element, &
    construct_LN_LN_element, corr_LN_to_covar_gaus, &
    corr_gaus_LN_to_covar_gaus, mu_LN_to_mu_gaus, &
    sigma_LN_to_sigma_gaus

  private ! Default scope

  contains

!-------------------------------------------------------------------------------
  subroutine generate_lh_sample &
             ( d_variables, hydromet_dim, & 
               wm, rcm, Ncnm, rvm, thlm, & 
               mixt_frac, &
               w1_in, w2_in, rc1_in, rc2_in, &
               varnce_w1_in, varnce_w2_in, &
               thl1_in, thl2_in, varnce_thl1_in, varnce_thl2_in, &
               rt1_in, rt2_in, varnce_rt1_in, varnce_rt2_in, &
               s1_in, s2_in, &
               stdev_s1_in, stdev_s2_in, &
               stdev_t1_in, stdev_t2_in, &
               covar_st_1, covar_st_2, &
               crt1, crt2, cthl1, cthl2, &
               hydromet, sigma2_on_mu2_ip_array_cloud, sigma2_on_mu2_ip_array_below, &
               corr_array_cloud, corr_array_below, &
               X_u_one_lev, X_mixt_comp_one_lev, &
               LH_rt, LH_thl, X_nl_one_lev ) ! Out
! Description:
!   This subroutine generates a Latin Hypercube sample.

! Assumptions:
!   The l_fix_s_t_correlations = false code does not set the correlation 
!   between Ncn and the other variates (i.e. it assumes they are all zero).
!   We do this is because while we have data for the 
!   correlation of e.g. s & Ncn and s & rr, we do not know the correlation of
!   Ncn and rr.
!   It would not be possible to decompose a covariance matrix with zero
!   correlation between rr and Ncn when the correlation between s and Ncn is
!   non-zero, and the code would have to halt.
!
!   One implication of this is that if l_fix_s_t_correlations = false 
!   then the correlation of s and Ncn must be set to 
!   zero in the correlation file to check the convergence of a non-interactive
!   SILHS solution against the analytic K&K solution.
!
!   The l_fix_s_t_correlations = true code does not have the above limitation
!   but will use a value for the covariance of the s and t that is not necessarily 
!   equal to the one computed by the PDF, so setting the correlation of 
!   s and Ncn to zero is not needed.
!   It will also fix the value of the correlation between s and t in the 
!   analytic K&K code, which should allow for convergence between the two solutions.
!   If it does not, then there is probably a new bug in the code.

! References:
!   ``Supplying Local Microphysical Parameterizations with Information about
!     Subgrid Variability: Latin Hypercube Sampling'', JAS Vol. 62,
!     p. 4010--4026, Larson, et al. 2005
!-------------------------------------------------------------------------------

    use constants_clubb, only:  &
      max_mag_correlation, &  ! Constant
      s_mellor_tol,  &  ! s tolerance in kg/kg
      t_mellor_tol,  &  ! t tolerance in kg/kg
      rt_tol, &         ! rt tolerance in kg/kg
      thl_tol, &        ! thetal tolerance in K
      w_tol_sqd, &      ! w^2 tolerance in m^2/s^2
      rc_tol, &         ! rc tolerance in kg/kg
      rr_tol, &         ! rr tolerance in kg/kg
      Nr_tol            ! Nr tolerance in #/kg

    use array_index, only: &
      iiNim,    & ! Variables
      iiNsnowm, &
      iiNrm,    &
      iirrainm, &
      iiricem, &
      iirsnowm, &
      iiNgraupelm, &
      iirgraupelm

    use corr_matrix_module, only: &
      iiPDF_rrain, & ! Variables
      iiPDF_rsnow, &
      iiPDF_rice, &
      iiPDF_rgraupel, &
      iiPDF_Nr, &
      iiPDF_Ncn, &
      iiPDF_Ni, &
      iiPDF_Nsnow, &
      iiPDF_Ngraupel, &
      iiPDF_s_mellor, &
      iiPDF_t_mellor, &
      iiPDF_w

    use latin_hypercube_arrays, only: &
      l_fixed_corr_initialized, & ! Variable(s)
      corr_stw_cloud_Cholesky, & 
      corr_stw_below_Cholesky, &
      corr_stw_cloud_scaling, & 
      corr_stw_below_scaling, &
      l_corr_stw_cloud_scaling, & 
      l_corr_stw_below_scaling

    use matrix_operations, only: &
      set_lower_triangular_matrix_dp, & ! Procedures
      get_lower_triangular_matrix, &
      row_mult_lower_tri_matrix, &
      print_lower_triangular_matrix

    use matrix_operations, only: Cholesky_factor ! Procedure(s)

    use matrix_operations, only: &
      symm_covar_matrix_2_corr_matrix ! Procedure(s)

    use error_code, only:  &
      clubb_at_least_debug_level  ! Procedure(s)

    use constants_clubb, only:  &
      fstderr, &  ! Constant(s)
      max_mag_correlation 

    use parameters_microphys, only: &
      l_fix_s_t_correlations ! Varible(s)

    use clubb_precision, only: &
      dp, & ! double precision
      core_rknd

    implicit none

    ! External
    intrinsic :: min, max, sqrt, null

    ! Input Variables
    integer, intent(in) :: &
      d_variables,   & ! `d' Number of variates (normally 3 + microphysics specific variables)
      hydromet_dim     ! Number of hydrometeor species

    real( kind = core_rknd ), dimension(hydromet_dim), intent(in) :: &
      hydromet ! Hydrometeor species [units vary]

    real( kind = core_rknd ), intent(in) :: &
      wm,   & ! Vertical velocity                                      [m/s]
      rcm,  & ! Mean liquid water mixing ratio                         [kg/kg]
      Ncnm, & ! Cloud nuclei concentration (simplified); Nc=Ncn*H(s)   [#/kg]
      rvm,  & ! Mean vapor water mixing ratio                          [kg/kg]
      thlm    ! Mean liquid potential temperature                      [K]

    real( kind = core_rknd ), intent(in) :: &
      mixt_frac,      & ! Mixture fraction                                        [-]
      rc1_in,         & ! Mean of rc for 1st normal distribution              [kg/kg]
      rc2_in,         & ! Mean of rc for 2nd normal distribution              [kg/kg]
      w1_in,          & ! Mean of w for 1st normal distribution                 [m/s]
      w2_in,          & ! Mean of w for 2nd normal distribution                 [m/s]
      varnce_w1_in,   & ! Variance of w for 1st normal distribution         [m^2/s^2]
      varnce_w2_in,   & ! Variance of w for 2nd normal distribution         [m^2/s^2]
      thl1_in,        & ! Mean of th_l for 1st normal distribution                [K]
      thl2_in,        & ! Mean of th_l for 2nd normal distribution                [K]
      varnce_thl1_in, & ! Variance of th_l for 1st normal distribution          [K^2]
      varnce_thl2_in, & ! Variance of th_l for 2nd normal distribution          [K^2]
      rt1_in,         & ! Mean of r_t for 1st normal distribution             [kg/kg]
      rt2_in,         & ! Mean of r_t for 2nd normal distribution             [kg/kg]
      varnce_rt1_in,  & ! Variance of r_t for 1st normal distribution     [kg^2/kg^2]
      varnce_rt2_in,  & ! Variance of r_t for 2nd normal distribution     [kg^2/kg^2]
      s1_in,          & ! Mean of s for 1st normal distribution               [kg/kg]
      s2_in,          & ! Mean of s for 2nd normal distribution               [kg/kg]
      stdev_s1_in,    & ! Standard deviation of s for 1st normal distribution [kg/kg]
      stdev_s2_in,    & ! Standard deviation of s for 2nd normal distribution [kg/kg]
      stdev_t1_in,    & ! Standard deviation of t for 1st normal distribution [kg/kg]
      stdev_t2_in,    & ! Standard deviation of t for 2nd normal distribution [kg/kg]
      covar_st_1,     & ! Covariance of s and t for 1st normal distribution   [kg/kg]
      covar_st_2,     & ! Covariance of s and t for 2nd normal distribution   [kg/kg]
      crt1,           & ! Coefficient for s'                                      [-]
      crt2,           & ! Coefficient for s'                                      [-]
      cthl1,          & ! Coefficient for s'                                    [1/K]
      cthl2             ! Coefficient for s'                                    [1/K]

    ! From the KK_microphys_module
    real( kind = core_rknd ), dimension(d_variables), target, intent(in) :: &
      sigma2_on_mu2_ip_array_cloud, & ! Ratio:  sigma_x^2/mu_x^2 ip; cloudy levs. [-]
      sigma2_on_mu2_ip_array_below    ! Ratio:  sigma_x^2/mu_x^2 ip; clear levs.  [-]

    real( kind = core_rknd ), dimension(d_variables,d_variables), target, intent(in) :: &
      corr_array_cloud, & ! Correlations in cloud for sampled variables  [-]
      corr_array_below    ! Correlations out of cloud                    [-]

    real( kind = dp ), intent(in), dimension(d_variables+1) :: & 
      X_u_one_lev ! Sample drawn from uniform distribution from a particular grid level

    integer, intent(in) :: &
      X_mixt_comp_one_lev ! Whether we're in the 1st or 2nd mixture component

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
      LH_rt, & ! Total water mixing ratio          [kg/kg]
      LH_thl   ! Liquid potential temperature      [K]

    real( kind = dp ), intent(out), dimension(d_variables) :: &
      X_nl_one_lev ! Sample that is transformed ultimately to normal-lognormal

    ! Local Variables

    logical, dimension(d_variables) :: &
      l_d_variable_lognormal ! Whether a given variable in X_nl has a lognormal dist.

    real( kind = core_rknd ) :: &
      rtm,         & ! Mean total water mixing ratio                           [kg/kg]
      s_mellor,    & ! Mean s_mellor (for when stdev_s1 < s_mellor_tol)        [kg/kg]
      w1,          & ! Mean of w for 1st normal distribution                     [m/s]
      w2,          & ! Mean of w for 2nd normal distribution                     [m/s]
      varnce_w1,   & ! Variance of w for 1st normal distribution             [m^2/s^2]
      varnce_w2,   & ! Variance of w for 2nd normal distribution             [m^2/s^2]
      thl1,        & ! Mean of th_l for 1st normal distribution                    [K]
      thl2,        & ! Mean of th_l for 2nd normal distribution                    [K]
      varnce_thl1, & ! Variance of th_l for 1st normal distribution              [K^2]
      varnce_thl2, & ! Variance of th_l for 2nd normal distribution              [K^2]
      rt1,         & ! Mean of r_t for 1st normal distribution                 [kg/kg]
      rt2,         & ! Mean of r_t for 2nd normal distribution                 [kg/kg]
      varnce_rt1,  & ! Variance of r_t for 1st normal distribution         [kg^2/kg^2]
      varnce_rt2,  & ! Variance of r_t for 2nd normal distribution         [kg^2/kg^2]
      s1,          & ! Mean of s for 1st normal distribution                   [kg/kg]
      s2,          & ! Mean of s for 2nd normal distribution                   [kg/kg]
      t1,          & ! Mean of t for 1st normal distribution                   [kg/kg]
      t2,          & ! Mean of t for 2nd normal distribution                   [kg/kg]
      stdev_s1,    & ! Standard deviation of s for 1st normal distribution     [kg/kg]
      stdev_s2,    & ! Standard deviation of s for 2nd normal distribution     [kg/kg]
      stdev_t1,    & ! Standard deviation of t for the 1st normal distribution [kg/kg]
      stdev_t2       ! Standard deviation of t for the 1st normal distribution [kg/kg]

    ! Means of s, t, w, & hydrometeors for plumes 1 and 2
    real( kind = core_rknd ), dimension(d_variables) :: &
      mu1, mu2

    ! Columns of Sigma_stw, X_nl_one_lev:  1   2   3   4 ... d_variables
    !                                      s   t   w   hydrometeors
    real( kind = dp ), dimension(d_variables,d_variables) :: &
      Sigma_stw_1, & ! Covariance of s,t, w + hydrometeors for plume 1
      Sigma_stw_2    ! Covariance of s,t, w + hydrometeors for plume 2

    real( kind = dp ) :: &
      var_Ncn1, & ! PDF param for width of plume 1.              [(#/kg)^2]
      var_Ncn2, & ! PDF param for width of plume 2.              [(#/kg^2]
      Nrm,      & ! Rain droplet number concentration.           [#/kg]
      var_Nr1,  & ! PDF param for width of plume 1.              [(#/kg)^2]
      var_Nr2     ! PDF param for width of plume 2.              [(#/kg^2]

    real( kind = core_rknd ) :: corr_rrNr, covar_rrNr1, covar_rrNr2, corr_srr, corr_sNr, &
            covar_sNr1, covar_sNr2, covar_srr1, covar_srr2

    real( kind = dp ) :: covar_trr1, covar_trr2, covar_tNr2, covar_tNr1

!   real :: &
!     stdev_Ncn, & ! Standard deviation of Ncn   [#/kg]
!     corr_tNcn, & ! Correlation between t and Ncn [-]
!     corr_sNcn, & ! Correlation between s and Ncn [-]
!     covar_tNcn1,    & ! Covariance of t and Ncn1      []
!     covar_tNcn2,    & ! Covariance of t and Ncn2      []
!     covar_sNcn1,    & ! Covariance of s and Ncn1      [# kg/kg^2]
!     covar_sNcn2       ! Covariance of s and Ncn2      [# kg/kg^2]

!   real( kind = dp ), dimension(2,2) :: corr_st_mellor_1, corr_st_mellor_2

    ! rr = specific rain content. [rr] = kg rain / kg air
    real( kind = dp ) :: &
      rrainm, &  ! rain water mixing ratio         [kg/kg]
      var_rr1, & ! PDF param for width of plume 1     [(kg/kg)^2]
      var_rr2    ! PDF param for width of plume 2.   [(kg/kg)^2]

    real( kind = dp ) :: &
      tp2_mellor_2, sp2_mellor_2,  & ! Variance of s,t         [(kg/kg)^2]
      sptp_mellor_2,               & ! Covariance of s and t   [kg/kg]
      tp2_mellor_1, sp2_mellor_1,  & ! Variance of s,t         [(kg/kg)^2]
      sptp_mellor_1                  ! Covariance of s and t   [kg/kg]


    real( kind = dp ), dimension(d_variables,d_variables) :: &
      Sigma1_Cholesky, Sigma2_Cholesky ! Cholesky factorization of Sigma1,2

    real( kind = dp ), dimension(d_variables) :: &
       Sigma1_scaling, & ! Scaling factors for Sigma1 for accuracy [units vary]
       Sigma2_scaling    ! Scaling factors for Sigma2 for accuracy [units vary]

    logical :: &
      l_Sigma1_scaling, l_Sigma2_scaling ! Whether we're scaling Sigma1 or Sigma2

    real( kind = dp ), dimension(d_variables,d_variables) :: &
      Corr_stw_1, Corr_stw_2 ! Correlation matrix for Sigma_stw_1,2

    real( kind = dp ), dimension(:,:), allocatable :: &
      corr_stw_matrix ! Correlation matrix      [-]

    real( kind = dp ), dimension(3) :: &
      temp_3_elements

    real( kind = core_rknd ), pointer, dimension(:) :: &
      sigma2_on_mu2_ip_array => null() ! Pointer for the sigma_x^2/mu_x^2 array

    real( kind = core_rknd ), pointer, dimension(:,:) :: &
      corr_array => null()  ! Correlation array pointer

    real( kind = dp ), pointer, dimension(:,:) :: &
      corr_stw_matrix_Cholesky => null() ! Pointer to the correct Cholesky factorization

    logical :: l_in_cloud

    integer :: i, index1, index2, ivar1, ivar2

    ! ---- Begin Code ----

    ! Determine which variables are a lognormal distribution
    i = max( iiPDF_s_mellor, iiPDF_t_mellor, iiPDF_w )
    l_d_variable_lognormal(1:i) = .false. ! The 1st 3 variates
    l_d_variable_lognormal(i+1:d_variables) = .true.  ! Hydrometeors

    ! Input pdf parameters.

    if ( l_fix_s_t_correlations ) then

      ! For fixed correlations, these don't appear in the correlation matrix, so
      ! we don't need them to be over some threshold.  In deep convective cases
      ! we don't want e.g. the variance of rt aloft to be rt_tol^2 necessarily,
      ! since this can lead to negative values of total water, so using the fixed s,t
      ! code may work better for those cases. -dschanen 15 Oct 2012

      ! Set means
      w1 = w1_in
      w2 = w2_in
      rt1 = rt1_in
      rt2 = rt2_in
      thl1 = thl1_in
      thl2 = thl2_in
      s1 = s1_in
      s2 = s2_in

      ! Set variances
      varnce_w1 = varnce_w1_in
      varnce_w2 = varnce_w2_in
      varnce_rt1 = varnce_rt1_in
      varnce_rt2 = varnce_rt2_in
      varnce_thl1 = varnce_thl1_in
      varnce_thl2 = varnce_thl2_in

      ! Set standard deviation of s1/s2
      stdev_s1 = stdev_s1_in
      stdev_s2 = stdev_s2_in

      stdev_t1 = stdev_t1_in
      stdev_t2 = stdev_t2_in

    else ! don't fixed the correlation of s and t.

      ! In this case we may need to set the variance to minimum value or the
      ! variance matrix cannot be decomposed

      call set_min_varnce_and_mean &
          ( wm, w_tol_sqd, w1_in, varnce_w1_in, & ! In
            varnce_w1, w1 ) ! Out

      call set_min_varnce_and_mean &
          ( wm, w_tol_sqd, w2_in, varnce_w2_in, & ! In
            varnce_w2, w2 ) ! Out

      rtm = rvm + rcm

      call set_min_varnce_and_mean &
          ( rtm, rt_tol**2, rt1_in, varnce_rt1_in, & ! In
            varnce_rt1, rt1 ) ! Out

      call set_min_varnce_and_mean &
          ( rtm, rt_tol**2, rt2_in, varnce_rt2_in, & ! In
            varnce_rt2, rt2 ) ! Out

      call set_min_varnce_and_mean &
          ( thlm, thl_tol**2, thl1_in, varnce_thl1_in, & ! In
            varnce_thl1, thl1 ) ! Out

      call set_min_varnce_and_mean &
          ( thlm, thl_tol**2, thl2_in, varnce_thl2_in, & ! In
            varnce_thl2, thl2 ) ! Out

      ! Compute the mean of s1 and s2
      s_mellor = s1_in * mixt_frac + (1.0_core_rknd-mixt_frac) * s2_in

      ! Here the subroutine name is a little misleading since we're imposing the
      ! threshold on a standard deviation rather than a variance.
      call set_min_varnce_and_mean &
          ( s_mellor, s_mellor_tol, s1_in, stdev_s1_in, & ! In
            stdev_s1, s1 ) ! Out

      ! See comment above.
      call set_min_varnce_and_mean &
          ( s_mellor, s_mellor_tol, s2_in, stdev_s2_in, & ! In
            stdev_s2, s2 ) ! Out

      ! The mean of t is zero;  we set the standard deviation to allow the 
      ! matrix to be decomposed for the t element
      stdev_t1 = max( stdev_t1_in, t_mellor_tol )
      stdev_t2 = max( stdev_t2_in, t_mellor_tol )

    end if ! l_fix_s_t_correlations

    !---------------------------------------------------------------------------
    ! Generate a set of sample points for a microphysics/radiation scheme
    !---------------------------------------------------------------------------

    ! We prognose rt-thl-w,
    !    but we set means, covariance of hydrometeors (e.g. rrain, Ncn) to constants.


    ! Standard sample for testing purposes when n=2
    ! X_u_one_lev(1,1:(d+1)) = ( / 0.0001_dp, 0.46711825945881_dp, &
    !             0.58015016959859_dp, 0.61894015386778_dp, 0.1_dp, 0.1_dp  / )
    ! X_u_one_lev(2,1:(d+1)) = ( / 0.999_dp, 0.63222458307464_dp, &
    !             0.43642762850981_dp, 0.32291562498749_dp, 0.1_dp, 0.1_dp  / )

    ! Select the in-cloud or out of cloud values if the correlations and
    ! x'^2 /  xm^2 terms.

    ! We define in cloud to be those points where mean liquid water is greater
    ! than rc_tol.  This should be done consistently with analytic K&K code.
    if ( X_mixt_comp_one_lev == 1 .and. rc1_in > rc_tol ) then
      sigma2_on_mu2_ip_array => sigma2_on_mu2_ip_array_cloud
      corr_array => corr_array_cloud
      l_in_cloud = .true.
    else if ( X_mixt_comp_one_lev == 2 .and. rc2_in > rc_tol ) then
      sigma2_on_mu2_ip_array => sigma2_on_mu2_ip_array_cloud
      corr_array => corr_array_cloud
      l_in_cloud = .true.
    else
      sigma2_on_mu2_ip_array => sigma2_on_mu2_ip_array_below
      corr_array => corr_array_below
      l_in_cloud = .false.

    end if

    ! Compute PDF parameters for Ncn, rr.
    ! Assume that Ncn, rr obey single-lognormal distributions

    ! Ncn  = cloud nuclei conc. (simplfied) / extended cloud droplet conc.;
    !        Nc = Ncn * H(s).  [Ncn] = number / kg air
    ! Ncnm  = mean of Ncn
    ! Ncnp2_on_Ncnm2 = variance of Ncn divided by Ncnm^2
    !  We must have a Ncnp2_on_Ncnm2 >= machine epsilon for the matrix
    ! Ncn1  = PDF parameter for mean of plume 1. [Ncn1] = (#/kg)
    ! Ncn2  = PDF parameter for mean of plume 2. [Ncn2] = (#/kg)

    if ( iiPDF_Ncn > 0 ) then
      call add_mu_element_LN &
           ( d_variables, iiPDF_Ncn, real( Ncnm, kind = dp ), & ! In
             sigma2_on_mu2_ip_array, & ! In
             mu1, mu2 ) ! In/out
    end if

    ! rr = specific rain content. [rr] = kg rain / kg air
    ! rrainm  = mean of rr; rrp2 = variance of rr, must have rrp2>0.
    ! rr1  = PDF parameter for mean of plume 1. [rr1] = (kg/kg)
    ! rr2  = PDF parameter for mean of plume 2. [rr2] = (kg/kg)

    if ( iiPDF_rrain > 0 ) then
      call add_mu_element_LN &
           ( d_variables, iiPDF_rrain, real(hydromet(iirrainm), kind = dp), & ! In
             sigma2_on_mu2_ip_array, &! In
             mu1, mu2 ) ! In/out
    end if

    if ( iiPDF_Nr > 0 ) then
      call add_mu_element_LN &
           ( d_variables, iiPDF_Nr, real(hydromet(iiNrm), kind = dp), & ! In
             sigma2_on_mu2_ip_array, & ! In
             mu1, mu2 ) ! In/out
    end if

    if ( iiPDF_rsnow > 0 ) then
      call add_mu_element_LN &
           ( d_variables, iiPDF_rsnow, real(hydromet(iirsnowm),kind = dp), & ! In
             sigma2_on_mu2_ip_array, & ! In
             mu1, mu2 ) ! In/out
    end if

    if ( iiPDF_Nsnow > 0 ) then
      call add_mu_element_LN &
           ( d_variables, iiPDF_Nsnow, real(hydromet(iiNsnowm),kind = dp), & ! In
             sigma2_on_mu2_ip_array, & ! In
             mu1, mu2 ) ! In/out
    end if

    if ( iiPDF_rice > 0 ) then
      call add_mu_element_LN &
           ( d_variables, iiPDF_rice, real(hydromet(iiricem),kind = dp), & ! In
             sigma2_on_mu2_ip_array, & ! In
             mu1, mu2 ) ! In/out
    end if

    if ( iiPDF_Ni > 0 ) then
      call add_mu_element_LN &
           ( d_variables, iiPDF_Ni, real(hydromet(iiNim),kind = dp), & ! In
             sigma2_on_mu2_ip_array, & ! In
             mu1, mu2 ) ! In/out
    end if

    if ( iiPDF_rgraupel > 0 ) then
      call add_mu_element_LN &
           ( d_variables, iiPDF_rgraupel, real(hydromet(iirgraupelm),kind = dp), & ! In
             sigma2_on_mu2_ip_array, & ! In
             mu1, mu2 ) ! In/out
    end if

    if ( iiPDF_Ngraupel > 0 ) then
      call add_mu_element_LN &
           ( d_variables, iiPDF_Ngraupel, real(hydromet(iiNgraupelm),kind = dp), & ! In
             sigma2_on_mu2_ip_array, & ! In
             mu1, mu2 ) ! In/out
    end if


    ! Means of s, t, w, Ncn, Nr, rr for Gaussians 1 and 2

    ! The mean of t is always 0.
    t1 = 0._core_rknd
    t2 = 0._core_rknd

    mu1((/iiPDF_s_mellor,iiPDF_t_mellor,iiPDF_w/)) &
      = (/ s1, t1, w1 /)
    mu2((/iiPDF_s_mellor,iiPDF_t_mellor,iiPDF_w/)) &
      = (/ s2, t2, w2 /)

    ! Define the variance of s and t
    tp2_mellor_1 = real(stdev_t1, kind = dp)**2
    tp2_mellor_2 = real(stdev_t2, kind = dp)**2

    sp2_mellor_1 = real(stdev_s1, kind = dp)**2
    sp2_mellor_2 = real(stdev_s2, kind = dp)**2

    ! An old subroutine, gaus_rotate, couldn't handle large correlations;
    !   I assume the replacement, gaus_condt, has equal trouble.
    !   Therefore we input smaller correlations
    !   The current code uses a Cholesky decomposition, which also cannot handle
    !   a correlation of exactly 1 without using the modified method -dschanen 11 Oct 2012
    ! max_mag_correlation = 0.99_core_rknd in constants.F90

    sptp_mellor_1 = real(min( max( -max_mag_correlation * stdev_t1 * stdev_s1, covar_st_1 ), &
      max_mag_correlation * stdev_t1 * stdev_s1 ), kind = dp)
    sptp_mellor_2 = real(min( max( -max_mag_correlation * stdev_t2 * stdev_s2, covar_st_2 ), &
      max_mag_correlation * stdev_t2 * stdev_s2 ), kind = dp)

    if ( .not. l_fix_s_t_correlations ) then

      ! Covariance (not correlation) matrices of rt-thl-w
      !    for Gaussians 1 and 2
      ! For now, assume no within-plume correlation of w with
      !    any other variables when the s and t correlations are not fixed.

      ! If l_Sigma_scaling = .true., we are dealing with a correlation Cholesky
      ! matrix here. But we need the covariance Cholesky matrix in sample_points.
      ! Therefore the results are rescaled in sample_points.

      ! Sigma_stw_1,2
      Sigma_stw_1 = 0._dp ! Start with no covariance, and add matrix elements
      Sigma_stw_2 = 0._dp

      ! Convert each Gaussian from rt-thl-w variables to s-t-w vars.

      ! Setup the Sigma matrices for s,t
      Sigma_stw_1(iiPDF_s_mellor,iiPDF_s_mellor) = sp2_mellor_1
      Sigma_stw_1(iiPDF_t_mellor,iiPDF_t_mellor) = tp2_mellor_1
      call set_lower_triangular_matrix_dp( 2, iiPDF_s_mellor, iiPDF_t_mellor, sptp_mellor_1, &
                                           Sigma_stw_1(1:2,1:2) )

      Sigma_stw_2(iiPDF_s_mellor,iiPDF_s_mellor) = sp2_mellor_2
      Sigma_stw_2(iiPDF_t_mellor,iiPDF_t_mellor) = tp2_mellor_2
      call set_lower_triangular_matrix_dp( 2, iiPDF_s_mellor, iiPDF_t_mellor, sptp_mellor_2, &
                                           Sigma_stw_2(1:2,1:2) )
      ! Add the w element
      Sigma_stw_1(iiPDF_w,iiPDF_w) = real(varnce_w1, kind = dp)

      Sigma_stw_2(iiPDF_w,iiPDF_w) = real(varnce_w2, kind = dp)

      if ( iiPDF_Ncn > 0 ) then
        ! var_Ncn1,2 = PDF param for width of plume 1,2. [var_Ncn1,2] = (#/kg)**2
        var_Ncn1 = log( 1._dp+ real( sigma2_on_mu2_ip_array(iiPDF_Ncn), kind=dp ) )
        var_Ncn2 = var_Ncn1
        Sigma_stw_1(iiPDF_Ncn,iiPDF_Ncn) = var_Ncn1
        Sigma_stw_2(iiPDF_Ncn,iiPDF_Ncn) = var_Ncn2
      end if

      if ( iiPDF_Nr > 0 ) then
        ! var_Nr1,2 = PDF param for width of plume 1,2. [var_Nr1,2] = (#/kg)**2
        var_Nr1 = log( 1._dp+ real( sigma2_on_mu2_ip_array(iiPDF_Nr), kind=dp ) )
        var_Nr2 = var_Nr1
        Sigma_stw_1(iiPDF_Nr,iiPDF_Nr) = var_Nr1
        Sigma_stw_2(iiPDF_Nr,iiPDF_Nr) = var_Nr2
      end if

      if ( iiPDF_rrain > 0 ) then
        ! var_rr1,2 = PDF param for width of plume 1,2. [var_rr1,2] = (kg/kg)**2
        var_rr1 = log( 1._dp+ real( sigma2_on_mu2_ip_array(iiPDF_rrain), kind=dp ) )
        var_rr2 = var_rr1
        Sigma_stw_1(iiPDF_rrain,iiPDF_rrain) = var_rr1
        Sigma_stw_2(iiPDF_rrain,iiPDF_rrain) = var_rr2
      end if

      if ( iiPDF_rrain > 0 .and. iiPDF_Nr > 0 ) then

        rrainm = real(hydromet(iirrainm), kind = dp)
        Nrm = real(hydromet(iiNrm), kind = dp)

        index1 = iiPDF_rrain
        index2 = iiPDF_Nr

        ! Covariance between rain water mixing ratio rain number concentration
        if ( rrainm > real(rr_tol, kind = dp) .and. Nrm > real(Nr_tol, kind = dp) ) then

          call get_lower_triangular_matrix &
               ( d_variables, index1, index2, corr_array, & ! In
                 corr_rrNr ) ! Out

          call construct_LN_LN_element &
               ( corr_rrNr, sigma2_on_mu2_ip_array(index1), & ! In
                 sigma2_on_mu2_ip_array(index2), & ! In
                 covar_rrNr1 ) ! Out

          ! rr1 = rr2 and Nr1 = Nr2, so we can just set covar_rrNr2 here
          covar_rrNr2 = covar_rrNr1

          call set_lower_triangular_matrix_dp &
               ( d_variables, index1, index2, real(covar_rrNr1, kind=dp), & ! In
                 Sigma_stw_1 ) ! In/out
          call set_lower_triangular_matrix_dp &
               ( d_variables, index1, index2, real(covar_rrNr2, kind=dp), & ! In
                 Sigma_stw_2 ) ! In/out
        end if

        index1 = iiPDF_s_mellor
        index2 = iiPDF_Nr
        ! Covariances involving s and Nr & rr
        if ( stdev_s1 > s_mellor_tol .and. Nrm > real(Nr_tol, kind = dp) ) then
          call get_lower_triangular_matrix &
               ( d_variables, index1, index2, corr_array, & ! In
                 corr_sNr ) ! Out

          ! Covariance between s and rain number conc.
          call construct_gaus_LN_element &
               ( corr_sNr, stdev_s1, sigma2_on_mu2_ip_array(index2), & ! In
                 covar_sNr1 ) ! Out

          call set_lower_triangular_matrix_dp &
               ( d_variables, index1, index2, real(covar_sNr1, kind=dp), & ! In
                 Sigma_stw_1 ) ! In/out

          ! Approximate the covariance of t and Nr
          ! This formula relies on the fact that iiPDF_s_mellor < iiPDF_t_mellor
          covar_tNr1 = ( Sigma_stw_1(iiPDF_t_mellor,iiPDF_s_mellor) &
            * real(covar_sNr1, kind = dp) ) / real(stdev_s1, kind = dp)**2

          call set_lower_triangular_matrix_dp &
               ( d_variables, iiPDF_t_mellor, iiPDF_Nr, real(covar_tNr1, kind=dp) , & ! In
                 Sigma_stw_1 ) ! In/out
        end if

        if ( stdev_s2 > s_mellor_tol .and. Nrm > real(Nr_tol, kind = dp) ) then

          call get_lower_triangular_matrix &
               ( d_variables, index1, index2, corr_array, & ! In
                 corr_sNr ) ! Out

          call construct_gaus_LN_element &
               ( corr_sNr, stdev_s2, sigma2_on_mu2_ip_array(index2), & ! In
                 covar_sNr2 ) ! Out

          call set_lower_triangular_matrix_dp &
               ( d_variables, index1, index2, real(covar_sNr2, kind=dp), & ! In
                 Sigma_stw_2 ) ! In/out

          ! Approximate the covariance of t and Nr
          ! This formula relies on the fact that iiPDF_s_mellor < iiPDF_t_mellor
          covar_tNr2 = ( Sigma_stw_2(iiPDF_t_mellor,iiPDF_s_mellor) &
            * real(covar_sNr2, kind = dp) ) / real(stdev_s2, kind = dp)**2

          call set_lower_triangular_matrix_dp &
               ( d_variables, iiPDF_t_mellor, iiPDF_Nr, real(covar_tNr2, kind = dp), & ! In
                 Sigma_stw_2 ) ! In/out
        end if

        index1 = iiPDF_s_mellor
        index2 = iiPDF_rrain
        ! Covariances involving s and Nr & rr
        if ( stdev_s1 > s_mellor_tol .and. rrainm > real(rr_tol, kind = dp) ) then

          call get_lower_triangular_matrix &
               ( d_variables, index1, index2, corr_array, & ! In
                 corr_srr ) ! Out

          ! Covariance between s and rain water mixing ratio
          call construct_gaus_LN_element &
               ( corr_srr, stdev_s1, sigma2_on_mu2_ip_array(index2), & ! In
                 covar_srr1 ) ! Out

          call set_lower_triangular_matrix_dp &
               ( d_variables, iiPDF_s_mellor, iiPDF_rrain, real(covar_srr1, kind=dp), & ! In
                 Sigma_stw_1 ) ! In/out

          ! Approximate the covariance of t and rr
          ! This formula relies on the fact that iiPDF_s_mellor < iiPDF_t_mellor
          covar_trr1 = ( Sigma_stw_1(iiPDF_t_mellor,iiPDF_s_mellor) &
            * real(covar_srr1, kind = dp) ) / real(stdev_s1, kind = dp)**2

          call set_lower_triangular_matrix_dp &
               ( d_variables, iiPDF_t_mellor, iiPDF_rrain, real(covar_trr1, kind = dp), & ! In
                 Sigma_stw_1 ) ! In/out
        end if

        if ( stdev_s2 > s_mellor_tol .and. rrainm > real( rr_tol, kind = dp ) ) then

          call get_lower_triangular_matrix &
               ( d_variables, index1, index2, corr_array, & ! In
                 corr_srr ) ! Out

          call construct_gaus_LN_element &
               ( corr_srr, stdev_s2, sigma2_on_mu2_ip_array(index2), & ! In
                 covar_srr2 ) ! Out

          call set_lower_triangular_matrix_dp &
               ( d_variables, index1, index2, real(covar_srr2, kind=dp), & ! In
                 Sigma_stw_2 ) ! In/out

          ! Approximate the covariance of t and rr
          ! This formula relies on the fact that iiPDF_s_mellor < iiPDF_t_mellor
          covar_trr2 = ( Sigma_stw_2(iiPDF_t_mellor,iiPDF_s_mellor) &
            * real(covar_srr2, kind = dp) ) / real(stdev_s2, kind = dp)**2

          call set_lower_triangular_matrix_dp &
               ( d_variables, iiPDF_t_mellor, iiPDF_rrain, real(covar_trr2, kind = dp), & ! In
                 Sigma_stw_2 ) ! In/out
        end if

      end if ! if iiPDF_rrain > 0 .and. iiPDF_Nr > 0

!     if ( iiPDF_Ncn > 0 ) then

      ! Covariances involving s and Ncn (currently disabled)
!       corr_sNcn = corr_array(iiPDF_s_mellor,iiPDF_Ncn)
!       stdev_Ncn = real( Ncnm, kind = core_rknd ) * sqrt( sigma2_on_mu2_ip_array(iiPDF_Ncn) )

!       if ( stdev_s1 > s_mellor_tol .and. Ncnm > real(Ncn_tol, kind = dp) ) then
!         ! The variable s is already Gaussian
!         stdev_sNcn1 = corr_gaus_LN_to_covar_gaus &
!                 ( corr_sNcn, &
!                   stdev_s1, &
!                   sigma_LN_to_sigma_gaus( sigma2_on_mu2_ip_array(iiPDF_Ncn) ) )

!         Sigma_stw_1(iiPDF_s_mellor,iiPDF_Ncn) = real(stdev_sNcn1, kind = dp)
!         Sigma_stw_1(iiPDF_Ncn,iiPDF_s_mellor) = real(stdev_sNcn1, kind = dp)

!         ! Approximate the covariance of t and Ncn
!         covar_tNcn1 = ( Sigma_stw_1(iiPDF_t_mellor,iiPDF_s_mellor) * covar_sNcn1 ) / stdev_s1**2

!         Sigma_stw_1(iiPDF_t_mellor,iiPDF_Ncn) = real(covar_tNcn1, kind = dp)
!         Sigma_stw_2(iiPDF_Ncn,iiPDF_t_mellor) = real(covar_tNcn2, kind = dp)

!       end if

!       if ( stdev_s2 > s_mellor_tol .and. Ncnm > real(Ncn_tol, kind = dp) ) then
!         stdev_sNcn2 = corr_gaus_LN_to_covar_gaus &
!                 ( corr_sNcn, &
!                   stdev_s2, &
!                   sigma_LN_to_sigma_gaus( sigma2_on_mu2_ip_array(iiPDF_Ncn) ) )

!         Sigma_stw_2(iiPDF_s_mellor,iiPDF_Ncn) = real(stdev_sNcn2, kind = dp)
!         Sigma_stw_2(iiPDF_Ncn,iiPDF_s_mellor) = real(stdev_sNcn2, kind = dp)

!         ! Approximate the covariance of t and Ncn
!         covar_tNcn2 = ( Sigma_stw_2(iiPDF_t_mellor,iiPDF_s_mellor) * covar_sNcn2 ) / stdev_s2**2

!         Sigma_stw_2(iiPDF_t_mellor,iiPDF_Ncn) = real(stNcn2, kind = dp)
!         Sigma_stw_2(iiPDF_Ncn,iiPDF_t_mellor) = real(stNcn2, kind = dp)

!       end if

!     end if ! iiPDF_Ncn > 0

      if ( clubb_at_least_debug_level( 2 ) ) then

        call symm_covar_matrix_2_corr_matrix( d_variables, Sigma_stw_1, Corr_stw_1 )
        call symm_covar_matrix_2_corr_matrix( d_variables, Sigma_stw_2, Corr_stw_2 )

        if ( any( Corr_stw_1 > 1.0_dp ) .or. any( Corr_stw_1 < -1.0_dp ) ) then
          write(fstderr,*) "Sigma_stw_1 has a correlation > 1 or < -1"
          call print_lower_triangular_matrix( fstderr, d_variables, &
            real( Corr_stw_1, kind = core_rknd ) )
        end if
        if ( any( Corr_stw_2 > 1.0_dp ) .or. any( Corr_stw_2 < -1.0_dp ) ) then
          write(fstderr,*) "Sigma_stw_2 has a correlation > 1 or < -1"
          call print_lower_triangular_matrix( fstderr, d_variables, &
            real( Corr_stw_2, kind = core_rknd ) )
        end if

      end if ! clubb_at_least_debug_level( 2 )

      ! Compute cholesky factorization Sigma_stw_1 / Sigma_stw_2
      if ( X_mixt_comp_one_lev == 1 ) then
        call Cholesky_factor( d_variables, Sigma_stw_1, & ! In
                              Sigma1_scaling, Sigma1_Cholesky, l_Sigma1_scaling ) ! Out
      end if

      if ( X_mixt_comp_one_lev == 2 ) then
        call Cholesky_factor( d_variables, real(Sigma_stw_2, kind = dp), & ! In
                              Sigma2_scaling, Sigma2_Cholesky, l_Sigma2_scaling ) ! Out
      end if

    else ! Using fixed correlations

      ! Here we are dealing with the correlation Cholesky matrix. Thus there is no scaling
      ! involved here, since the condition number of a correlation matrix is always one.
      ! Hence l_Sigma_scaling should always be .false. here. But since we need the covariance
      ! Cholesky matrix for the sampling, we have to convert the correlation Cholesky matrix
      ! before we feed it to sample points.
      !
      ! Attention: corr_stw_matrix_Cholesky is not a correlation matrix. It is rather a mixture
      ! of a correlation and covariance matrix (see description of costruct_corr_stw_matrix).

      ! Compute the Cholesky factorization of the correlations if it's not
      ! already computed.
      if ( .not. l_fixed_corr_initialized ) then

        allocate( corr_stw_matrix(d_variables,d_variables), &
                  corr_stw_cloud_Cholesky(d_variables,d_variables), &
                  corr_stw_below_Cholesky(d_variables,d_variables), &
                  corr_stw_cloud_scaling(d_variables), &
                  corr_stw_below_scaling(d_variables) )

        call construct_corr_stw_matrix &
             ( d_variables, corr_array_cloud, & ! In
               sigma2_on_mu2_ip_array_cloud, & ! In
               corr_stw_matrix ) ! Out

        ! Compute choleksy factorization for the correlation matrix (in cloud)
        call Cholesky_factor( d_variables, real(corr_stw_matrix, kind = dp), & ! In
                              corr_stw_cloud_scaling, corr_stw_cloud_Cholesky, & ! Out
                              l_corr_stw_cloud_scaling ) ! Out

        ! This subroutine constructs a matrix where the first 3x3 elements are correlations
        ! and the other elements are LN covariances
        call construct_corr_stw_matrix &
             ( d_variables, corr_array_below, & ! In
               sigma2_on_mu2_ip_array_below, & ! In
               corr_stw_matrix ) ! Out

        ! Compute choleksy factorization for the correlation matrix (out of cloud)
        call Cholesky_factor( d_variables, real(corr_stw_matrix, kind = dp), & ! In
                              corr_stw_below_scaling, corr_stw_below_Cholesky, &  ! Out
                              l_corr_stw_below_scaling ) ! Out


        deallocate( corr_stw_matrix )

        l_fixed_corr_initialized = .true.

      end if

      ! Determine if the point is in or out of cloud for the purposes of picking
      ! the values of the correlations to use
      if ( l_in_cloud ) then
        l_Sigma1_scaling = l_corr_stw_cloud_scaling
        l_Sigma2_scaling = l_corr_stw_cloud_scaling
        corr_stw_matrix_Cholesky => corr_stw_cloud_Cholesky
        Sigma1_scaling = corr_stw_cloud_scaling
        Sigma2_scaling = corr_stw_cloud_scaling
      else
        l_Sigma1_scaling = l_corr_stw_below_scaling
        l_Sigma2_scaling = l_corr_stw_below_scaling
        corr_stw_matrix_Cholesky => corr_stw_below_Cholesky
        Sigma1_scaling = corr_stw_below_scaling
        Sigma2_scaling = corr_stw_below_scaling
      end if

      if ( X_mixt_comp_one_lev == 1 ) then

        Sigma1_Cholesky = 0._dp ! Initialize the variance to zero

        temp_3_elements = (/ real(stdev_s1, kind = dp), real(stdev_t1, kind = dp),&
                             sqrt( real(varnce_w1, kind = dp) ) /)

        ! Multiply the first three elements of the variance matrix by the
        ! values of the standard deviation of s1, t1, and w1
        call row_mult_lower_tri_matrix &
             ( 3, temp_3_elements, corr_stw_matrix_Cholesky(1:3,1:3), & ! In
               Sigma1_Cholesky(1:3,1:3) ) ! Out

        ! Set the remaining elements (the lognormal variates) to the value
        ! contained in the matrix, since they don't vary in space in time
        do ivar1 = 4, d_variables
          do ivar2 = 4, ivar1
            Sigma1_Cholesky(ivar1,ivar2) = corr_stw_matrix_Cholesky(ivar1,ivar2)
          end do
        end do
      end if ! X_mixt_comp_one_lev == 1

      if ( X_mixt_comp_one_lev == 2 ) then
        Sigma2_Cholesky = 0._dp

        temp_3_elements = (/ real(stdev_s2, kind = dp), real(stdev_t2, kind = dp),&
                             sqrt( real(varnce_w2, kind = dp) ) /)

        ! Multiply the first three elements of the variance matrix by the
        ! values of the standard deviation of s2, t2, and w2
        call row_mult_lower_tri_matrix &
             ( 3, temp_3_elements, corr_stw_matrix_Cholesky(1:3,1:3), & ! In
               Sigma2_Cholesky(1:3,1:3) ) ! Out

        ! Set the remaining elements (the lognormal variates) to the value
        ! contained in the matrix, since they don't vary in space in time
        do ivar1 = 4, d_variables
          do ivar2 = 4, ivar1
            Sigma2_Cholesky(ivar1,ivar2) = corr_stw_matrix_Cholesky(ivar1,ivar2)
          end do
        end do
      end if ! X_mixt_comp_one_lev == 2

    end if ! l_fix_s_t_correlations

    ! Compute the new set of sample points using the update variance matrices
    ! for this level
    call sample_points( d_variables, 1, &  ! intent(in
                        real(rt1, kind = dp), real(thl1, kind = dp), &  ! intent(in)
                        real(rt2, kind = dp), real(thl2, kind = dp), &  ! intent(in)
                        real(crt1, kind = dp), real(cthl1, kind = dp), &  ! intent(in)
                        real(crt2, kind = dp), real(cthl2, kind = dp), &  ! intent(in)
                        mu1, mu2, &  ! intent(in)
                        l_d_variable_lognormal, & ! intent(in)
                        X_u_one_lev, & ! intent(in)
                        X_mixt_comp_one_lev, & ! intent(in)
                        Sigma1_Cholesky, Sigma2_Cholesky, & ! intent(in)
                        Sigma1_scaling, Sigma2_scaling, & ! intent(in)
                        l_Sigma1_scaling, l_Sigma2_scaling, & ! intent(in)
                        LH_rt, LH_thl, X_nl_one_lev ) ! intent(out)

    return
  end subroutine generate_lh_sample

!-------------------------------------------------------------------------------
  subroutine generate_lh_sample_mod &
             ( d_variables, d_uniform_extra, & ! In
               thl1, thl2, rt1, rt2, & ! In
               crt1, crt2, cthl1, cthl2, & ! In
               mu1, mu2, sigma1, sigma2, & ! In
               corr_stw_matrix_Cholesky_1, & ! In
               corr_stw_matrix_Cholesky_2, & ! In
               X_u_one_lev, X_mixt_comp_one_lev, & ! In
               l_in_precip_one_lev, & ! In
               LH_rt, LH_thl, X_nl_one_lev ) ! Out
! Description:
!   This subroutine generates a Latin Hypercube sample.

! Assumptions:
!   The l_fix_s_t_correlations = false code does not set the correlation
!   between Ncn and the other variates (i.e. it assumes they are all zero).
!   We do this is because while we have data for the
!   correlation of e.g. s & Ncn and s & rr, we do not know the correlation of
!   Ncn and rr.
!   It would not be possible to decompose a covariance matrix with zero
!   correlation between rr and Ncn when the correlation between s and Ncn is
!   non-zero, and the code would have to halt.
!
!   One implication of this is that if l_fix_s_t_correlations = false
!   then the correlation of s and Ncn must be set to
!   zero in the correlation file to check the convergence of a non-interactive
!   SILHS solution against the analytic K&K solution.
!
!   The l_fix_s_t_correlations = true code does not have the above limitation
!   but will use a value for the covariance of the s and t that is not necessarily
!   equal to the one computed by the PDF, so setting the correlation of
!   s and Ncn to zero is not needed.
!   It will also fix the value of the correlation between s and t in the
!   analytic K&K code, which should allow for convergence between the two solutions.
!   If it does not, then there is probably a new bug in the code.

! References:
!   ``Supplying Local Microphysical Parameterizations with Information about
!     Subgrid Variability: Latin Hypercube Sampling'', JAS Vol. 62,
!     p. 4010--4026, Larson, et al. 2005
!-------------------------------------------------------------------------------

    use corr_matrix_module, only: &
      iiPDF_s_mellor, & ! Variable(s)
      iiPDF_t_mellor, &
      iiPDF_w

    use matrix_operations, only: &
      row_mult_lower_tri_matrix ! Procedures


    use constants_clubb, only:  &
      one

    use clubb_precision, only: &
      dp, & ! double precision
      core_rknd

    implicit none

    ! External
    intrinsic :: max

    ! Input Variables
    integer, intent(in) :: &
      d_variables, &  ! `d' Number of variates (normally 3 + microphysics specific variables)
      d_uniform_extra ! Number of variates included in uniform sample only (often 2)

    real( kind = core_rknd ), intent(in) :: &
      rt1,         & ! Mean of r_t for 1st normal distribution   [kg/kg]
      rt2,         & ! Mean of r_t for 2nd normal distribution   [kg/kg]
      thl1,        & ! Mean of th_l for 1st normal distribution  [K]
      thl2,        & ! Mean of th_l for 2nd normal distribution  [K]
      crt1,        & ! Coefficient for s'                        [-]
      crt2,        & ! Coefficient for s'                        [-]
      cthl1,       & ! Coefficient for s'                        [1/K]
      cthl2          ! Coefficient for s'                        [1/K]

    real( kind = dp ), dimension(d_variables,d_variables), intent(in) :: &
      corr_stw_matrix_Cholesky_1, & ! Correlations Cholesky matrix (1st comp.)  [-]
      corr_stw_matrix_Cholesky_2    ! Correlations Cholesky matrix (2nd comp.)  [-]

    real( kind = core_rknd ), dimension(d_variables), intent(in) :: &
      mu1,    & ! Means of the hydrometeors, 1st comp. (s, t, w, <hydrometeors>)  [units vary]
      mu2,    & ! Means of the hydrometeors, 2nd comp. (s, t, w, <hydrometeors>)  [units vary]
      sigma1, & ! Stdevs of the hydrometeors, 1st comp. (s, t, w, <hydrometeors>) [units vary]
      sigma2    ! Stdevs of the hydrometeors, 2nd comp. (s, t, w, <hydrometeors>) [units vary]

    real( kind = dp ), intent(in), dimension(d_variables+d_uniform_extra) :: &
      X_u_one_lev ! Sample drawn from uniform distribution from a particular grid level

    integer, intent(in) :: &
      X_mixt_comp_one_lev ! Whether we're in the 1st or 2nd mixture component

    logical, intent(in) :: &
      l_in_precip_one_lev ! Whether we are in precipitation (T/F)

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
      LH_rt, & ! Total water mixing ratio          [kg/kg]
      LH_thl   ! Liquid potential temperature      [K]

    real( kind = dp ), intent(out), dimension(d_variables) :: &
      X_nl_one_lev ! Sample that is transformed ultimately to normal-lognormal

    ! Local Variables

    logical, dimension(d_variables) :: &
      l_d_variable_lognormal ! Whether a given variable in X_nl has a lognormal dist.

    real( kind = dp ), dimension(d_variables,d_variables) :: &
      Sigma1_Cholesky, Sigma2_Cholesky ! Cholesky factorization of Sigma1,2

    real( kind = dp ), dimension(d_variables) :: &
       Sigma1_scaling, & ! Scaling factors for Sigma1 for accuracy [units vary]
       Sigma2_scaling    ! Scaling factors for Sigma2 for accuracy [units vary]

    logical :: &
      l_Sigma1_scaling, l_Sigma2_scaling ! Whether we're scaling Sigma1 or Sigma2

!    real( kind = dp ), dimension(3) :: &
!      temp_3_elements

    integer :: i !, ivar1, ivar2

    ! ---- Begin Code ----

    ! Determine which variables are a lognormal distribution
    i = max( iiPDF_s_mellor, iiPDF_t_mellor, iiPDF_w )
    l_d_variable_lognormal(1:i) = .false. ! The 1st 3 variates
    l_d_variable_lognormal(i+1:d_variables) = .true.  ! Hydrometeors

    !---------------------------------------------------------------------------
    ! Generate a set of sample points for a microphysics/radiation scheme
    !---------------------------------------------------------------------------

    ! We prognose rt-thl-w,
    !    but we set means, covariance of hydrometeors (e.g. rrain, Ncn) to constants.


    ! Standard sample for testing purposes when n=2
    ! X_u_one_lev(1,1:(d+1)) = ( / 0.0001_dp, 0.46711825945881_dp, &
    !             0.58015016959859_dp, 0.61894015386778_dp, 0.1_dp, 0.1_dp  / )
    ! X_u_one_lev(2,1:(d+1)) = ( / 0.999_dp, 0.63222458307464_dp, &
    !             0.43642762850981_dp, 0.32291562498749_dp, 0.1_dp, 0.1_dp  / )

    ! The old code was dealing with Covariance matrices. These were rescaled by
    ! Sigma'(i,j) = 1/sqrt(Sigma(i,i)) * Sigma(i,j) * 1/Sigma(j,j)
    ! to reduce the condition number of the matrices. Sigma' is the correlation
    ! matrix. This code deals directly with the correlation matrix. Hence we don't
    ! need any rescaling here.
    l_Sigma1_scaling = .false.
    l_Sigma2_scaling = .false.
    Sigma1_scaling = one
    Sigma2_scaling = one

    if ( X_mixt_comp_one_lev == 1 ) then

      Sigma1_Cholesky = 0._dp ! Initialize the variance to zero

      ! Multiply the first three elements of the variance matrix by the
      ! values of the standard deviation of s1, t1, and w1
      call row_mult_lower_tri_matrix &
           ( d_variables, real( sigma1, kind = dp ), corr_stw_matrix_Cholesky_1, & ! In
             Sigma1_Cholesky ) ! Out

    elseif ( X_mixt_comp_one_lev == 2 ) then
      Sigma2_Cholesky = 0._dp

      ! Multiply the first three elements of the variance matrix by the
      ! values of the standard deviation of s2, t2, and w2
      call row_mult_lower_tri_matrix &
           ( d_variables, real( sigma2, kind = dp ), corr_stw_matrix_Cholesky_2, & ! In
             Sigma2_Cholesky ) ! Out

    end if ! X_mixt_comp_one_lev == 1

    ! Compute the new set of sample points using the update variance matrices
    ! for this level
    call sample_points( d_variables, d_uniform_extra, &  ! intent(in)
                        real(rt1, kind = dp), real(thl1, kind = dp), &  ! intent(in)
                        real(rt2, kind = dp), real(thl2, kind = dp), &  ! intent(in)
                        real(crt1, kind = dp), real(cthl1, kind = dp), &  ! intent(in)
                        real(crt2, kind = dp), real(cthl2, kind = dp), &  ! intent(in)
                        mu1, mu2, &  ! intent(in)
                        l_d_variable_lognormal, & ! intent(in)
                        X_u_one_lev, & ! intent(in)
                        X_mixt_comp_one_lev, & ! intent(in)
                        Sigma1_Cholesky, Sigma2_Cholesky, & ! intent(in)
                        Sigma1_scaling, Sigma2_scaling, & ! intent(in)
                        l_Sigma1_scaling, l_Sigma2_scaling, & ! intent(in)
                        LH_rt, LH_thl, X_nl_one_lev ) ! intent(out)

    ! Zero precipitation hydrometeors if not in precipitation
    if ( .not. l_in_precip_one_lev ) then

      call zero_precip_hydromets( d_variables, & ! Intent(in)
                                X_nl_one_lev ) ! Intent(inout)

    end if


    return
  end subroutine generate_lh_sample_mod

!---------------------------------------------------------------------------------------------------
  subroutine zero_precip_hydromets( d_variables, X_nl_one_lev )

  ! Description:
  !   Sets the sample values of the precipitating hydrometeors to zero

  ! References:
  !   None
  !-----------------------------------------------------------------------------

    use constants_clubb, only: &
      zero_dp

    use corr_matrix_module, only: &
      iiPDF_Ncn ! Variable(s)

    implicit none

    ! Input Variables

    integer, intent(in) :: &
      d_variables  ! Number of hydrometeors                                        [count]

    ! Input/Output Variables

    real( kind = dp ), intent(inout), dimension(d_variables) :: &
      X_nl_one_lev      ! Sample of hydrometeors (normal-lognormal space)          [units vary]

    ! Local Variables
    integer :: &
      ivar              ! Loop counter                                             [count]
  !-----------------------------------------------------------------------------

    !----- Begin Code -----

    do ivar = iiPDF_Ncn+1, d_variables
      X_nl_one_lev(ivar) = zero_dp
    end do

  end subroutine zero_precip_hydromets

!---------------------------------------------------------------------------------------------------
  subroutine sample_points( d_variables, d_uniform_extra, &
                            rt1, thl1, rt2, thl2, &
                            crt1, cthl1, crt2, cthl2, &
                            mu1, mu2,  &
                            l_d_variable_lognormal, &
                            X_u_one_lev, &
                            X_mixt_comp_one_lev, &
                            Sigma1_Cholesky, Sigma2_Cholesky, &
                            Sigma1_scaling, Sigma2_scaling, &
                            l_Sigma1_scaling, l_Sigma2_scaling, &
                            LH_rt, LH_thl, X_nl_one_lev )

! Description:
!   Generates n random samples from a d-dim Gaussian-mixture PDF.
!   Uses Latin hypercube method.

!   Original formulation takes samples only from the cloudy part of the grid box.
!   Revised formulation samples in and out of cloud.

!   We use MKS units on all variates.

! References:
!   None
!----------------------------------------------------------------------

    use corr_matrix_module, only: &
      iiPDF_s_mellor, & ! Variables
      iiPDF_t_mellor

    use clubb_precision, only: &
      core_rknd, & ! Variable(s)
      dp

    implicit none

    ! Input variables
    integer, intent(in) :: &
      d_variables, &    ! Number of variates
      d_uniform_extra   ! Variates included in uniform sample only  

    !rt1, thl1 = mean of rt, thl for Gaus comp 1
    !rt2, thl2 = mean of rt, thl for Gaus comp 2
    real( kind = dp ), intent(in) :: rt1, thl1, rt2, thl2

    ! Thermodynamic constants for plumes 1 and 2, units of kg/kg
    real( kind = dp ), intent(in) :: &
      crt1,  & ! coefficient relating rt, s and t for Gaus comp 1
      cthl1, & ! coeff relating thl, s and t for component 1
      crt2,  & ! coefficient relating rt, s and t for component 2
      cthl2    ! coefficient relating thl, s and t for comp. 2

    ! Latin hypercube variables, i.e. s, t, w, etc.
    real( kind = core_rknd ), intent(in), dimension(d_variables) :: &
      mu1, mu2 ! d-dimensional column vector of means of 1st, 2nd components

    logical, intent(in), dimension(d_variables) :: &
      l_d_variable_lognormal ! Whether a given element of X_nl is lognormal

    real( kind = dp ), intent(in), dimension(d_variables+d_uniform_extra) :: &
      X_u_one_lev ! Sample drawn from uniform distribution from particular grid level [-]

    integer, intent(in) :: &
      X_mixt_comp_one_lev ! Whether we're in the 1st or 2nd mixture component

    ! Columns of Sigma_Cholesky, X_nl_one_lev:  1   2   3   4 ... d_variables
    !                                           s   t   w   hydrometeors
    real( kind = dp ), intent(in), dimension(d_variables,d_variables) :: &
      Sigma1_Cholesky, & ! [units vary]
      Sigma2_Cholesky

    real( kind = dp ), intent(in), dimension(d_variables) :: &
      Sigma1_scaling, Sigma2_scaling ! Scaling factors on Sigma1,2 [units vary]

    logical, intent(in) :: &
      l_Sigma1_scaling, l_Sigma2_scaling ! Whether we're scaling Sigma1 or Sigma2

    ! Output Variables
    ! Total water, theta_l: mean plus perturbations
    real( kind = core_rknd ), intent(out) :: &
      LH_rt,  & ! Total water   [kg/kg]
      LH_thl    ! Liquid potential temperature  [K]

    real( kind = dp ), intent(out), dimension(d_variables) :: &
      X_nl_one_lev ! Sample that is transformed ultimately to normal-lognormal

    ! ---- Begin Code ----

    ! Generate n samples of a d-variate Gaussian mixture
    ! by transforming Latin hypercube points, X_u_one_lev.
    call gaus_mixt_points( d_variables, d_uniform_extra, mu1, mu2, &  ! intent(in)
                           Sigma1_Cholesky, Sigma2_Cholesky, & ! intent(in)
                           Sigma1_scaling, Sigma2_scaling, & ! intent(in)
                           l_Sigma1_scaling, l_Sigma2_scaling, & ! intent(in)
                           X_u_one_lev, X_mixt_comp_one_lev, & ! intent(in)
                           X_nl_one_lev ) ! intent(out)

! Transform s (column 1) and t (column 2) back to rt and thl
!   This is only needed if you need total water mixing ratio and liquid potential
!   samples points.
!     call sptp_2_rtpthlp &
!          ( 1, d_variables, mixt_frac, crt1, cthl1, crt2, cthl2, &
!            cloud_frac1, cloud_frac2, X_nl_one_lev(1), &
!            X_nl_one_lev(2), &
!            X_u_one_lev, rtp, thlp )
    call st_2_rtthl( rt1, thl1, rt2, thl2, & ! intent(in)
                     crt1, cthl1, crt2, cthl2, & ! intent(in)
                     real(mu1(iiPDF_s_mellor), kind = dp), & ! intent(in)
                     real(mu2(iiPDF_s_mellor), kind = dp), & ! intent(in)
                     X_nl_one_lev(iiPDF_s_mellor), & ! intent(in)
                     X_nl_one_lev(iiPDF_t_mellor), & ! intent(in)
                     X_mixt_comp_one_lev, & ! intent(in)
                     LH_rt, LH_thl ) ! intent(out)

    ! Convert lognormal variates (e.g. Ncn and rr) to lognormal
    where ( l_d_variable_lognormal )
      X_nl_one_lev(:) = exp( X_nl_one_lev(:) )
    end where

    return
  end subroutine sample_points

!-------------------------------------------------------------------------------

!-----------------------------------------------------------------------
  subroutine rtpthlp_2_sptp( stdev_s_mellor, varnce_rt, varnce_thl, &
                             rrtthl_covar, crt, cthl, &
                             tp2, sp2, sptp ) ! Out

! Description:
!   Transform covariance matrix from rt', theta_l' coordinates
!   to s', t' coordinates.
!   Use linear approximation for s', t'.
!
! References:
!   ``Supplying Local Microphysics Parameterizations with Information about
!     Subgrid Variability: Latin Hypercube Sampling'', V.E. Larson et al.,
!     JAS 62 pp. 4015
!
! Notes:
!   We no longer use this subroutine since the code in pdf_closure will compute
!   the value of the variance of t, s, and the covariance of the two directly.
!-----------------------------------------------------------------------

    use constants_clubb, only: &
      max_mag_correlation ! Constant

    implicit none

    ! External
    intrinsic :: min, max, sqrt

    ! Input Variables

    real( kind = dp ), intent(in) :: &
      stdev_s_mellor, & ! Standard deviation of s_mellor [(kg/kg)^2]
      varnce_rt,      & ! Variance of rt1/rt2
      varnce_thl,     & ! Variance of thl1/thl2
      rrtthl_covar,   & ! Covariance of rt, thl
      crt, cthl         ! Coefficients that define s', t'

    ! Output Variables

    real( kind = dp ), intent(out) :: &
      tp2, sp2,  &    ! Variance of s,t         [(kg/kg)^2]
      sptp            ! Covariance of s and t   [kg/kg]

    ! Local Variables

    real( kind = dp ) :: crt_sqd, cthl_sqd

    real( kind = dp ) :: &
      sqrt_sp2_tp2, & ! sqrt of the product of the variances of s and t [kg/kg]
      max_mag_corr_dp

    ! ---- Begin Code ----

    ! Simplified formula. Here we compute the variance t_mellor and covariance
    ! of s_mellor and t_mellor using formula's derived from the matrix
    ! multiplication on Larson, et al. See figure 14.
    crt_sqd = crt**2
    cthl_sqd = cthl**2
    sptp = crt_sqd * varnce_rt - cthl_sqd * varnce_thl
    tp2 = crt_sqd * varnce_rt + 2._dp * crt * cthl * rrtthl_covar &
        + cthl_sqd * varnce_thl
    sp2 = stdev_s_mellor**2

    ! Reduce the correlation of s and t Mellor if it's greater than 0.99
    sqrt_sp2_tp2 = sqrt( sp2 * tp2 )
    max_mag_corr_dp = real( max_mag_correlation, kind=dp )
    sptp = min( max( -max_mag_corr_dp * sqrt_sp2_tp2, sptp ), &
                max_mag_corr_dp * sqrt_sp2_tp2 )

!   sptp = 0.3 * sqrt( sp2 ) * sqrt( tp2 )

    return
  end subroutine rtpthlp_2_sptp
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
  subroutine generate_uniform_sample( num_samples, nt_repeat, n_vars, p_matrix, X_u_one_lev )

! Description:
!   Generates a matrix X that contains a Latin Hypercube sample.
!   The sample is uniformly distributed.
! References:
!   See Art B. Owen (2003), ``Quasi-Monte Carlo Sampling,"
!      a chapter from SIGGRAPH 2003
!-------------------------------------------------------------------------------

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      num_samples,  & ! `n'  Number of samples generated
      nt_repeat,    & ! `n_t' Num. random samples before sequence repeats
      n_vars           ! Number of uniform variables to generate

    integer, intent(in), dimension(num_samples,n_vars) :: &
      p_matrix    !num_samples x n_vars array of permuted integers

    ! Output Variables

    real(kind=dp), intent(out), dimension(num_samples,n_vars) :: &
      X_u_one_lev ! num_samples by n_vars matrix, X
                  ! each row of which is a n_vars-dimensional sample

    ! Local Variables

    integer :: j, k

    ! ---- Begin Code ----

!  Compute random permutation row by row
!       do j=1,dp1
!       ! Generate a column vector of integers from 0 to n-1,
!       !    whose order is random.
!         call rand_permute( n, p_matrix(1:n,j) )
!       end do

    ! Choose values of sample using permuted vector and random number generator
    do j = 1,num_samples
      do k = 1,n_vars
        X_u_one_lev(j,k) = choose_permuted_random( nt_repeat, p_matrix(j,k) )
      end do
    end do

    return
  end subroutine generate_uniform_sample

!----------------------------------------------------------------------
  function choose_permuted_random( nt_repeat, p_matrix_element )

! Description:
!   Chooses a permuted random using the Mersenne Twister algorithm.
!
! References:
!   None
!----------------------------------------------------------------------

    use mt95, only: genrand_real3 ! Procedure(s)

    use mt95, only: genrand_real ! Constants

    implicit none

    ! Input Variables
    integer, intent(in) :: & 
      nt_repeat,        & ! Number of samples before the sequence repeats
      p_matrix_element    ! Permuted integer

    ! Output Variable
    real(kind=genrand_real) :: choose_permuted_random

    ! Local Variable
    real(kind=genrand_real) :: & 
      rand ! Random float with a range of (0,1)

    ! ---- Begin Code ----

    call genrand_real3( rand ) ! genrand_real3's range is (0,1)

    choose_permuted_random = (1.0_genrand_real/real( nt_repeat, kind=genrand_real) ) &
       *( real( p_matrix_element, kind=genrand_real ) + rand )

    return
  end function choose_permuted_random

!----------------------------------------------------------------------
  subroutine gaus_mixt_points( d_variables, d_uniform_extra, mu1, mu2, & ! Intent(in)
                               Sigma1_Cholesky, Sigma2_Cholesky, & ! Intent(in)
                               Sigma1_scaling, Sigma2_scaling, & ! Intent(in)
                               l_Sigma1_scaling, l_Sigma2_scaling, & ! Intent(in)
                               X_u_one_lev, X_mixt_comp_one_lev, & ! Intent(in)
                               X_nl_one_lev ) ! Intent(out)
! Description:
!   Generates n random samples from a d-dimensional Gaussian-mixture PDF.
!   Uses Latin hypercube method.
!
! References:
!   None
!----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! double precision

    implicit none

    ! Input Variables

    integer, intent(in) :: &
      d_variables, &    ! Number of variates
      d_uniform_extra   ! Variates included in uniform sample only

    real( kind = core_rknd ), intent(in), dimension(d_variables) :: &
      mu1, mu2 ! d-dimensional column vector of means of 1st, 2nd Gaussians

    ! Latin hypercube sample from uniform distribution from a particular grid level
    real( kind = dp ), intent(in), dimension(d_variables+d_uniform_extra) :: &
      X_u_one_lev

    real( kind = dp ), dimension(d_variables,d_variables), intent(in) :: &
      Sigma1_Cholesky, Sigma2_Cholesky ! Cholesky factorization of Sigma1,2

    real( kind = dp ), dimension(d_variables), intent(in) :: &
      Sigma1_scaling, & ! Scaling factors for Sigma1 for accuracy [units vary]
      Sigma2_scaling    ! Scaling factors for Sigma2 for accuracy [units vary]

    logical, intent(in) :: &
      l_Sigma1_scaling, l_Sigma2_scaling ! Whether we're scaling Sigma1 or Sigma2

    integer, intent(in) :: &
      X_mixt_comp_one_lev ! Which mixture component we're in

    ! Output Variables

    real( kind = dp ), intent(out), dimension(d_variables) :: &
      X_nl_one_lev ! [n by d] matrix, each row of which is a d-dimensional sample

    ! Local Variables

    real( kind = dp ), dimension(d_variables) :: &
      std_normal  ! Standard normal multiplied by the factorized Sigma    [-]

    integer :: ivar ! Loop iterators

    ! ---- Begin Code ----

    ! From Latin hypercube sample, generate standard normal sample
    do ivar = 1, d_variables
      std_normal(ivar) = ltqnorm( X_u_one_lev(ivar) )
    end do


      ! Determine which mixture fraction we are in.
    if ( X_mixt_comp_one_lev == 1 ) then

      call multiply_Cholesky &
          ( d_variables, std_normal, & ! In
            mu1, Sigma1_Cholesky, &  ! In
            Sigma1_scaling, l_Sigma1_scaling, & ! In
            X_nl_one_lev(1:d_variables) ) ! Out

    else if ( X_mixt_comp_one_lev == 2 ) then

      call multiply_Cholesky &
           ( d_variables, std_normal, & ! In
             mu2, Sigma2_Cholesky, &  ! In
             Sigma2_scaling, l_Sigma2_scaling, & ! In
             X_nl_one_lev(1:d_variables) ) ! Out

    else
      stop "Error determining mixture component in gaus_mixt_points"

    end if ! X_mixt_comp_one_lev

    return
  end subroutine gaus_mixt_points

!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
  subroutine truncate_gaus_mixt( d_variables, col, mixt_frac, mu1, mu2, & 
                  Sigma1, Sigma2, cloud_frac1, cloud_frac2, X_u_one_lev, &
                  X_mixt_comp_one_lev, truncated_column )
! Description:
!   Converts sample points drawn from a uniform distribution
!    to truncated Gaussian points.
!
! References:
!   None
!-------------------------------------------------------------------------------

    use constants_clubb, only:  &
      fstderr  ! Constant(s)

    use error_code, only:  &
      clubb_at_least_debug_level  ! Procedure(s)

    use clubb_precision, only: &
      core_rknd ! double precision

    implicit none

    ! Input Variables

    integer, intent(in) :: &
      d_variables,   &  ! Number of variates (normally=5)
      col               ! Scalar indicated which column of X_nl_one_lev to truncate

    real( kind = dp ), intent(in) :: &
      mixt_frac,    & ! Mixture fraction of Gaussians
      cloud_frac1, cloud_frac2  ! Cloud fraction associated w/ 1st, 2nd mixture component

    real( kind = core_rknd ), intent(in), dimension(d_variables) :: &
      mu1, mu2 ! d-dimensional column vector of means of 1st, 2nd Gaussians

    real( kind = dp ), intent(in), dimension(d_variables,d_variables) :: &
      Sigma1, Sigma2 ! dxd dimensional covariance matrices

    ! Latin hypercube sample from uniform distribution from a particular grid level
    real( kind = dp ), intent(in), dimension(d_variables+1) :: &
      X_u_one_lev

    integer, intent(in) :: &
      X_mixt_comp_one_lev ! Whether we're in the first or second mixture component

    ! Output Variables

    ! A column vector of length n that is transformed from a Gaussian PDF to truncated Gaussian PDF.
    real( kind = dp ), intent(out) :: truncated_column

    ! Local Variables

    real( kind = dp ) :: s_std

    ! ---- Begin Code ----

    ! Handle some possible errors re: proper ranges of mixt_frac,
    ! cloud_frac1, cloud_frac2.
    if ( (mixt_frac > 1.0_dp) .or. (mixt_frac < 0.0_dp) ) then
      write(fstderr,*) 'Error in truncate_gaus_mixt:  ',  &
                       'mixture fraction, mixt_frac, does not lie in [0,1].'
      stop
    end if
    if ( (cloud_frac1 > 1.0_dp) .or. (cloud_frac1 < 0.0_dp) ) then
      write(fstderr,*) 'Error in truncate_gaus_mixt:  ',  &
                       'cloud fraction 1, cloud_frac1, does not lie in [0,1].'
      stop
    end if
    if ( (cloud_frac2 > 1.0_dp) .or. (cloud_frac2 < 0.0_dp) ) then
      write(fstderr,*) 'Error in truncate_gaus_mixt:  ',  &
                       'cloud fraction 2, cloud_frac2, does not lie in [0,1].'
      stop
    end if

    ! Make sure there is some cloud.
    if (mixt_frac*cloud_frac1 < 0.001_dp .and. (1._dp-mixt_frac) * cloud_frac2 < 0.001_dp) then
      if ( clubb_at_least_debug_level( 1 ) ) then
        write(fstderr,*) 'Error in truncate_gaus_mixt:  ',  &
                         'there is no cloud or almost no cloud!'
      end if
    end if

    ! Make s PDF (1st column) a truncated Gaussian.
    ! This allows us to sample solely from the cloud points.

    ! Choose which mixture fraction we are in.
    ! Account for cloud fraction.
    ! Follow M. E. Johnson (1987), p. 56.
    if  ( X_mixt_comp_one_lev == 1 ) then
      ! Replace first dimension (s) with
      !  sample from cloud (i.e. truncated standard Gaussian)
      s_std = ltqnorm( X_u_one_lev(col) * cloud_frac1 &
                 + (1._dp - cloud_frac1 ) )
      ! Convert to nonstandard normal with mean mu1 and variance Sigma1
      truncated_column =  & 
                 s_std * sqrt( Sigma1(col,col) ) + real(mu1(col), kind = dp)
    else if ( X_mixt_comp_one_lev == 2 ) then

      ! Replace first dimension (s) with
      !   sample from cloud (i.e. truncated Gaussian)
      s_std = ltqnorm( (X_u_one_lev(col) * cloud_frac2) &
            + (1._dp - cloud_frac2) )

        ! Convert to nonstandard normal with mean mu2 and variance Sigma2
      truncated_column =  & 
                    s_std * sqrt( Sigma2(col,col) ) + real(mu2(col), kind = dp)
    else
      stop "Error in truncate_gaus_mixt"
    end if

    return
  end subroutine truncate_gaus_mixt

!-----------------------------------------------------------------------
  function ltqnorm( p )
! Description:
!   This function is ported to Fortran from the same function written in Matlab,
!    see the following description of this function.  Hongli Jiang, 2/17/2004
!   Converted to double precision by Vince Larson 2/22/2004;
!    this improves results for input values of p near 1.

! LTQNORM Lower tail quantile for standard normal distribution.
!
!   Z = LTQNORM(P) returns the lower tail quantile for the standard normal
!   distribution function.  I.e., it returns the Z satisfying Pr{X < Z} = P,
!   where X has a standard normal distribution.
!
!   LTQNORM(P) is the same as SQRT(2) * ERFINV(2*P-1), but the former returns a
!   more accurate value when P is close to zero.

!   The algorithm uses a minimax approximation by rational functions and the
!   result has a relative error less than 1.15e-9. A last refinement by
!   Halley's rational method is applied to achieve full machine precision.

!   Author:      Peter J. Acklam
!   Time-stamp:  2003-04-23 08:26:51 +0200
!   E-mail:      pjacklam@online.no
!   URL:         http://home.online.no/~pjacklam
!-----------------------------------------------------------------------

    ! This is commented out because it is used by old versions of this code that
    ! have been commented out, but not by the new code.

    ! use constants_clubb, only: Pi_DP ! Variable(s)

    use clubb_precision, only: &
      dp ! double precision

    implicit none

    ! External

    intrinsic :: log, sqrt

    ! Input Variable(s)

    real( kind = dp ), intent(in) :: p

    ! Return Variable

    real( kind = dp ) :: ltqnorm


    ! Local Variable(s)

    real( kind = dp ) a1, a2, a3, a4, a5, a6, b1, b2, b3, b4, b5, & 
                     c1, c2, c3, c4, c5, c6, d1, d2, d3, d4

    real( kind = dp ) q, r, z, z1, plow, phigh

!       double preciseion e, erf_dp, u

!       Occurs in constants.F now.  Isn't actually used currently.
!        double precision, parameter :: pi=3.1415926_dp

! Coefficients in rational approximations.
! equivalent: a(1)=a1, a(2)=a2, and etc, when a(1) is in Matlab.
! Similarly for b, c, and d's
    parameter (a1 = -3.969683028665376E+01_dp,  &
               a2 = 2.209460984245205E+02_dp, &
               a3 = -2.759285104469687E+02_dp,  &
               a4 = 1.383577518672690E+02_dp, &
               a5 = -3.066479806614716E+01_dp,  &
               a6 = 2.506628277459239E+00_dp)
    parameter (b1 = -5.447609879822406E+01_dp,  &
               b2 = 1.615858368580409E+02_dp, &
               b3 = -1.556989798598866E+02_dp,  &
               b4 = 6.680131188771972E+01_dp, &
               b5 = -1.328068155288572E+01_dp)
    parameter (c1 = -7.784894002430293E-03_dp,  &
               c2 = -3.223964580411365E-01_dp, &
               c3 = -2.400758277161838E+00_dp,  &
               c4 = -2.549732539343734E+00_dp, &
               c5 =  4.374664141464968E+00_dp,  &
               c6 =  2.938163982698783E+00_dp)
    parameter (d1 =  7.784695709041462E-03_dp,  &
               d2 =  3.224671290700398E-01_dp, &
               d3 =  2.445134137142996E+00_dp,  &
               d4 =  3.754408661907416E+00_dp)

    ! Default initialization
    z = 0.0_dp

!  Define break-points.
    plow  = 0.02425_dp
    phigh = 1._dp - plow

!  Initialize output array. Don't need this in Fortran
!   z = zeros(size(p));

!  Rational approximation for lower region:
    if (p > 0._dp .and. p < plow) then
      q = sqrt( -2._dp * log( p ) )
      z = (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/ & 
                ((((d1*q+d2)*q+d3)*q+d4)*q+1._dp)
!  Rational approximation for central region:
    else if (p >= plow .and. p <= phigh) then
      q = p - 0.5_dp
      r = q * q
      z = (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q & 
                 /(((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1._dp)
! Rational approximation for upper region:
    else if (p > phigh .and. p < 1._dp) then
      q  = sqrt( -2._dp * log(1._dp - p) )
      z  = -(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) & 
                  /((((d1*q+d2)*q+d3)*q+d4)*q+1._dp)
    end if

!  Case when P = 0: z = -inf, to create inf z =-1.0.
!     to create NaN's inf*inf.
    z1 = 0._dp
    if (p == 0._dp) then
      z = (-1._dp)/z1
    end if

! Case when P = 1:, z=inf
    if(p == 1._dp)then
      z = 1._dp/z1
    end if

!  Cases when output will be NaN:
!   k = p < 0 | p > 1 | isnan(p);
! usually inf*inf --> NaN's.
    if (p < 0._dp .or. p > 1._dp) then
      z = (1._dp/z1)**2
    end if

!  The relative error of the approximation has absolute value less
!  than 1.15e-9. One iteration of Halley's rational method (third
!  order) gives full machine precision.
! V. Larson 20Feb04: Don't use the following if-end if loop.
!   The value of e is very different than what MATLAB produces,
!   possibly because of
!   poor values of erf from Numerical Recipes.
!   The value is close to MATLAB's
!   if I omit the following if-end if loop.
! End V. Larson comment
!!   k = 0 < p & p < 1;
!       if (p.gt.0 .and. p.lt.1)then
!         e = 0.5_core_rknd*(1.0_core_rknd - erf_dp(-z/sqrt(2._core_rknd)) - p          ! error
!         u = e * sqrt(2*pi_dp) * exp(z**2/2)       ! f(z)/df(z)
!         z = z - u/( 1 + z*u/2 )               ! Halley's method
!       end if

! return z as double precision:
    ltqnorm = z

    return
  end function ltqnorm

!-------------------------------------------------------------------------------
  subroutine multiply_Cholesky( d_variables, std_normal, mu, Sigma_Cholesky, & 
                                  Sigma_scaling, l_scaled, &
                                  nonstd_normal )
! Description:
!   Computes the nonstd_normal from the Cholesky factorization of Sigma,
!   std_normal, and mu.
!   nonstd_normal = Sigma_Cholesky * std_normal + mu.

! References:
!   M. E. Johnson (1987), ``Multivariate Normal and Related Distributions'' p50-55
!-------------------------------------------------------------------------------

    use clubb_precision, only: &
      dp, & ! double precision
      core_rknd

    implicit none

    external :: dtrmv ! BLAS upper/lower triangular multiplication subroutine

    ! Parameters
    integer, parameter :: &
      incx = 1 ! Increment for x in dtrmv

    ! Input Variables
    integer, intent(in) :: d_variables ! Number of variates (normally=5)

    real( kind = dp ), intent(in), dimension(d_variables) :: &
      std_normal ! vector of d-variate standard normal distribution [-]

    real( kind = core_rknd ), intent(in), dimension(d_variables) :: &
      mu ! d-dimensional column vector of means of Gaussian     [units vary]

    real( kind = dp ), intent(in), dimension(d_variables,d_variables) :: &
      Sigma_Cholesky ! Cholesky factorization of the Sigma matrix [units vary]

    real( kind = dp ), intent(in), dimension(d_variables) :: &
      Sigma_scaling ! Scaling for Sigma / mu    [units vary]

    logical, intent(in) :: l_scaled ! Whether any scaling was done to Sigma

    ! Output Variables

    ! nxd matrix of n samples from d-variate normal distribution
    !   with mean mu and covariance structure Sigma
    real( kind = dp ), intent(out) :: &
      nonstd_normal(d_variables)

    ! Local Variables
    real( kind = dp ), dimension(d_variables) :: &
      Sigma_times_std_normal ! Sigma * std_normal [units vary]

    ! --- Begin Code ---

    Sigma_times_std_normal = std_normal ! Copy std_normal into 'x'

    ! Call the level 2 BLAS subroutine to multiply std_normal by Sigma_Cholesky
    call dtrmv( 'Lower', 'N', 'N', d_variables, Sigma_Cholesky, d_variables, & ! In
                Sigma_times_std_normal, & ! In/out
                incx ) ! In

    if ( l_scaled ) then
      ! Add mu to Sigma * std_normal (scaled)
      nonstd_normal = Sigma_times_std_normal + real(mu, kind = dp) * Sigma_scaling
      ! Determine 'y' vector by removing the scaling factors
      nonstd_normal = nonstd_normal / Sigma_scaling
    else
      ! Add mu to Sigma * std_normal
      nonstd_normal = Sigma_times_std_normal + real(mu, kind = dp)
    end if

    return
  end subroutine multiply_Cholesky
!-----------------------------------------------------------------------
  subroutine st_2_rtthl( rt1, thl1, rt2, thl2, & 
                         crt1, cthl1, crt2, cthl2, & 
                         mu_s1, mu_s2, &
                         s_mellor, t_mellor, X_mixt_comp_one_lev, &
                         LH_rt, LH_thl )
! Description:
!   Converts from s, t variables to rt, thl.  Also sets a limit on the value
!   of cthl1 and cthl2 to prevent extreme values of temperature.
!
! References:
!   None
!-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! double precision

    implicit none

    ! External

    intrinsic :: max, real

    ! Constant Parameters

    ! Reduce the below value if model seems to crashing due excessive
    ! lh_thlp2_zt and it will limit the extremes of the samples.
!   real(kind = dp), parameter :: &
!     cthl_thresh = 1e-5_dp ! Threshold on cthl1 and cthl2 [kg/kg/K]

    real(kind = dp), parameter :: &
      thl_dev_lim = 5.0_dp ! Max deviation from mean thetal [K]

    ! Input Variables

    real( kind = dp ), intent(in) :: &
      rt1, rt2,    & ! n dimensional column vector of rt         [kg/kg]
      thl1, thl2,  & ! n dimensional column vector of thetal     [K]
      crt1, crt2,  & ! Constants from plumes 1 & 2 of rt
      cthl1, cthl2   ! Constants from plumes 1 & 2 of thetal

    real( kind = dp ), intent(in) :: &
      mu_s1, mu_s2 ! Mean for s1 and s2         [kg/kg]

    ! n-dimensional column vector of Mellor's s and t, including mean and perturbation
    real( kind = dp ), intent(in) :: &
      s_mellor, &  ! [kg/kg]
      t_mellor     ! [-]

    integer, intent(in) :: &
      X_mixt_comp_one_lev ! Whether we're in the first or second mixture component

    ! Output variables

    real( kind = core_rknd ), intent(out) :: &
      LH_rt, LH_thl ! n-dimensional column vectors of rt and thl, including mean and perturbation

    ! Local Variables

!   real( kind= dp ) :: cthl1_clip, cthl2_clip, & ! Clipped values of cthl1,2 [kg/kg/K]
    real( kind= dp ) :: LH_dev_thl_lim ! Limited value of the deviation on thetal [K]

    ! ---- Begin Code ----

    ! Clip the value of cthl1,2.  This prevents large values of theta_l when the
    ! saturation point is low and limits the chance of instability.
    ! See ticket #527 on the CLUBB TRAC
!   cthl1_clip = max( cthl1, cthl_thresh )
!   cthl2_clip = max( cthl2, cthl_thresh )

    ! Choose which mixture fraction we are in.
    ! Account for cloud fraction.
    ! Follow M. E. Johnson (1987), p. 56.
!     fraction_1     = mixt_frac*cloud_frac1 / &
!                      (mixt_frac*cloud_frac1+(1-mixt_frac)*cloud_frac2)

    if ( X_mixt_comp_one_lev == 1 ) then
      LH_rt  = real( rt1 + (0.5_dp/crt1)*(s_mellor-mu_s1) +  & 
                             (0.5_dp/crt1)*t_mellor, kind=core_rknd )

      ! Limit the quantity that temperature can vary by (in K)
      LH_dev_thl_lim = (-0.5_dp/cthl1)*(s_mellor-mu_s1) & 
                     + (0.5_dp/cthl1)*t_mellor

      LH_dev_thl_lim = max( min( LH_dev_thl_lim, thl_dev_lim ), -thl_dev_lim )

      LH_thl = real( thl1 + LH_dev_thl_lim, kind=core_rknd )

        ! Old code
!       LH_thl = real( thl1 + (-0.5_dp/cthl1_clip)*(s_mellor-mu_s1) +  & 
!                              (0.5_dp/cthl1_clip)*t_mellor, kind=core_rknd )

    else if ( X_mixt_comp_one_lev == 2 ) then
        ! mixture fraction 2
      LH_rt = real( rt2 + (0.5_dp/crt2)*(s_mellor-mu_s2) +  & 
                             (0.5_dp/crt2)*t_mellor, kind=core_rknd )

      ! Limit the quantity that temperature can vary by (in K)
      LH_dev_thl_lim = (-0.5_dp/cthl2)*(s_mellor-mu_s2) & 
                     + (0.5_dp/cthl2)*t_mellor

      LH_dev_thl_lim = max( min( LH_dev_thl_lim, thl_dev_lim ), -thl_dev_lim )

      LH_thl = real( thl2 + LH_dev_thl_lim, kind=core_rknd )

      ! Old code
!     LH_thl = real( thl2 + (-0.5_dp/cthl2_clip)*(s_mellor-mu_s2) +  & 
!                           (0.5_dp/cthl2_clip)*t_mellor, kind=core_rknd )
    else
      stop "Error determining mixture fraction in st_2_rtthl"

    end if

    return
  end subroutine st_2_rtthl
!-------------------------------------------------------------------------------
  subroutine log_sqd_normalized( Xm, sigma2_on_mu2_ip, &
                                 X1, X2 )
! Description:
!
! References:
!   None
!-------------------------------------------------------------------------------

    use clubb_precision, only: &
        dp ! double precision

    implicit none

    ! External
    intrinsic :: log, epsilon, max

    ! Input Variables
    real( kind = dp ), intent(in) :: &
      Xm,         & ! Mean X          [units vary]
      sigma2_on_mu2_ip    ! X'^2 / X^2      [-]

    ! Output Variables
    real( kind = dp ), intent(out) :: &
      X1, X2  ! PDF parameters for mean of plume 1, 2   [units vary]

    ! ---- Begin Code ----

    ! Here the variables X1 & X2 have ambiguous units.  After the exp function
    ! is applied to the result at the very end of this code the units will be
    ! correct because of the 0.5 coefficient. I.e. sqrt( Xm^2 ) = Xm.
    ! Here we use epsilon to impose a limit on the numerator to prevent
    ! taking the log of 0 while still imposing an upper bound.
    X1 = 0.5_dp * log( max( Xm, epsilon( Xm ) )**2 / ( 1._dp + sigma2_on_mu2_ip ) )
    X2 = X1

    return
  end subroutine log_sqd_normalized

!-------------------------------------------------------------------------------
  subroutine construct_gaus_LN_element( corr_sy, stdev_s, y_sigma2_on_mu2, &
                                        covar_sy )

! Description:
!   Compute the covariance of s_mellor and a lognormal variate,
!   converting from lognormal to gaussian space as required.
!
! References:
!   None
!-------------------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    real( kind = core_rknd ), intent(in) :: &
      corr_sy,         & ! Correlation between x and y                       [-]
      stdev_s,         & ! Standard deviation of s                       [kg/kg]
      y_sigma2_on_mu2    ! Ratio:  sigma_y^2 over mu_y^2 (ith PDF comp.) ip  [-]

    real( kind = core_rknd ), intent(out) :: covar_sy

    real( kind = core_rknd ) :: y_sigma2_on_mu2_gaus

    ! ---- Begin Code ----

    y_sigma2_on_mu2_gaus = sigma_LN_to_sigma_gaus( y_sigma2_on_mu2 )

    covar_sy = corr_gaus_LN_to_covar_gaus( corr_sy, stdev_s, y_sigma2_on_mu2_gaus )

    return
  end subroutine construct_gaus_LN_element
!-------------------------------------------------------------------------------
  subroutine construct_LN_LN_element( corr_xy, x_sigma2_on_mu2_ip, &
                                      y_sigma2_on_mu2_ip, &
                                      covar_xy )

! Description:
!   Compute the covariance of 2 variables, converting from lognormal
!   to gaussian space as required.
! References:
!   None
!-------------------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    real( kind = core_rknd ), intent(in) :: &
      corr_xy,            & ! Correlation btwn. x and y (ith PDF comp.) ip  [-]
      x_sigma2_on_mu2_ip, & ! Ratio:  sigma_x^2 / mu_x^2 (ith PDF comp.) ip [-]
      y_sigma2_on_mu2_ip    ! Ratio:  sigma_y^2 / mu_y^2 (ith PDF comp.) ip [-]

    real( kind = core_rknd ), intent(out) :: covar_xy

    real( kind = core_rknd ) :: sigma_x_gaus, sigma_y_gaus

    ! ---- Begin Code ----

    sigma_x_gaus = sigma_LN_to_sigma_gaus( x_sigma2_on_mu2_ip )
    sigma_y_gaus = sigma_LN_to_sigma_gaus( y_sigma2_on_mu2_ip )

    covar_xy = corr_LN_to_covar_gaus( corr_xy, sigma_x_gaus, sigma_y_gaus )

    return
  end subroutine construct_LN_LN_element

!-------------------------------------------------------------------------------
  subroutine set_min_varnce_and_mean( Xmean, varnce_X_tol, Xn, varnce_Xn, &
                                      varnce_Xn_out, Xn_out )
! Description:
!   Here we impose a threshold on the variance of a term (usually from the PDF)
!   so as to avoid numerical instability.  If the variance is too small we then
!   set the variable associated with the nth mixture component to the mean
!   value.  E.g., if varnce_rt1 is small we set it to rt_tol^2, then set rt1 = rtm.
! References:
!   None
!-------------------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    real( kind = core_rknd ), intent(in) :: &
      Xmean,            & ! Mean value of variable X             [units vary]
      varnce_X_tol,     & ! Min tolerance on the variance of X   [units vary]
      Xn,               & ! Value of X for the 1st norm. dist.   [units vary]
      varnce_Xn           ! Variance of X for the 2nd norm. dist.[units vary]

    real( kind = core_rknd ), intent(out) :: &
      varnce_Xn_out, & ! [units vary]
      Xn_out           ! [units vary]

    ! --- Begin Code ----

    if ( varnce_Xn > varnce_X_tol ) then
      ! The variance is large enough, so we use it.
      varnce_Xn_out = varnce_Xn
      ! Keep X1
      Xn_out = Xn
    else
      ! Set the variance to a small number to prevent a singular matrix
      varnce_Xn_out = varnce_X_tol
      ! Set X to the mean value
      Xn_out = Xmean
    end if

    return
  end subroutine set_min_varnce_and_mean

!-------------------------------------------------------------------------------
  subroutine construct_corr_stw_matrix &
             ( d_variables, corr_array, &
               sigma2_on_mu2_ip_array, &
               corr_stw_matrix )
! Description:
!   Construct a correlation matrix containing s,t,w and the lognormal variates.
!   This code is only called when l_fix_s_t_correlations is true.  It does not
!   assume zero correlation between w and the other variates.
!
!   The matrix is i.e. a mixture, where the first 3 rows contain the correlations for s,t,w and
!   the remaining elements are normalized covariances. So the matrix looks like
!
!   --                     --
!   |   1                   |
!   |   c   1               |
!   |   c   c   1           |
!   |   k   k   k   v       |
!   |   k   k   k   k   v   |
!   |      [  . . .  ]      |
!   --                     --
!
!   where c = correlation, k = covariance (normalized) and v = variance (normalized).
!
! References:
!   None.
!-------------------------------------------------------------------------------
    use corr_matrix_module, only: &
      iiPDF_s_mellor, & ! Variable(s)
      iiPDF_t_mellor, &
      iiPDF_w

    use matrix_operations, only: &
      set_lower_triangular_matrix_dp, & ! Procedures
      get_lower_triangular_matrix

    use clubb_precision, only: &
      dp, &! double precision
      core_rknd

    implicit none

    ! External
    intrinsic :: log

    ! Input Variables
    integer, intent(in) :: d_variables ! Number of variates

    real( kind = core_rknd ), dimension(d_variables,d_variables), intent(in) :: &
      corr_array ! Correlations between variates

    real( kind = core_rknd ), dimension(d_variables), intent(in) :: &
      sigma2_on_mu2_ip_array  ! sigma_^2 / mu^2 (ith PDF comp.) in-precip.  [-]

    ! Output variables
    real( kind = dp ), dimension(d_variables,d_variables), intent(out) :: &
      corr_stw_matrix ! Correlations between variates with some terms in lognormal space

    ! Local Variables
    real( kind = core_rknd ) :: &
      corr_st, & ! Correlation for s,t  [-]
      corr_sw, & ! Correlation for w,s  [-]
      corr_tw    ! Correlation for w,t  [-]

    integer :: i, index1, index2, LN_index

    ! ---- Begin Code ----

    ! LN_index is the location of the first variate which is lognormally
    ! distributed (e.g. rain water mixing ratio).
    LN_index = max( iiPDF_s_mellor, iiPDF_t_mellor, iiPDF_w )+1

    corr_stw_matrix = 0.0_dp ! Initialize to 0

    do i = 1, LN_index-1, 1
      ! Set main diagonal to 1
      corr_stw_matrix(i,i) = 1.0_dp
    end do

    ! Set the correlation of s and t. For this part of the code we assume a
    ! fixed correlation in order to only compute the Cholesky factorization
    ! once per simulation.
    index1 = iiPDF_s_mellor
    index2 = iiPDF_t_mellor

    call get_lower_triangular_matrix &
         ( d_variables, index1, index2, corr_array, & ! In
            corr_st ) ! Out

    call set_lower_triangular_matrix_dp &
         ( d_variables, index1, index2, real(corr_st, kind = dp), & ! In
           corr_stw_matrix ) ! In/out

    ! The correlation between w and s,t is typically not fixed either, but for
    ! the reasons listed above we compute it using a fixed value.
    index1 = iiPDF_s_mellor
    index2 = iiPDF_w

    call get_lower_triangular_matrix &
         ( d_variables, index1, index2, corr_array, & ! In
            corr_sw ) ! Out
    call set_lower_triangular_matrix_dp &
         ( d_variables, index1, index2, real(corr_sw, kind = dp), & ! In
           corr_stw_matrix ) ! In/out

    ! Obtain the fixed value for the correlation between t and w.
    index1 = iiPDF_t_mellor
    index2 = iiPDF_w
    call get_lower_triangular_matrix &
         ( d_variables, index1, index2, corr_array, & ! In
            corr_tw ) ! Out

    ! Add the correlation to the matrix
    call set_lower_triangular_matrix_dp &
         ( d_variables, index1, index2, real(corr_tw, kind = dp), & ! In
           corr_stw_matrix ) ! In/out

    ! Compute the main diagonal for each lognormal variate
    forall ( i = LN_index:d_variables )
      corr_stw_matrix(i,i)&
      = real(log( 1._core_rknd + sigma2_on_mu2_ip_array(i) ), kind = dp)
    end forall

    do index1 = LN_index, d_variables
      do index2 = LN_index, index1
        ! Add all lognormal covariances
        call add_corr_to_matrix_LN_LN &
             ( d_variables, index1, index2, & ! In
               sigma2_on_mu2_ip_array, corr_array, & ! In
               corr_stw_matrix ) ! In/Out
      end do
    end do

    ! Correlations involving s, t, w and the lognormal variates
    do index1 = LN_index, d_variables
      call add_corr_to_matrix_gaus_LN &
           ( d_variables, iiPDF_s_mellor, & ! In
             iiPDF_t_mellor, iiPDF_w, index1, & ! In
             sigma2_on_mu2_ip_array, corr_array, & ! In
             corr_stw_matrix ) ! In/Out
    end do

    return
  end subroutine construct_corr_stw_matrix

!-------------------------------------------------------------------------------
  subroutine add_corr_to_matrix_LN_LN( d_variables, index1, index2, &
                                       sigma2_on_mu2_ip_array, corr_array, &
                                       corr_stw_matrix )
! Description:
!   Added a correlation between two lognormally distributed variates to a
!   correlation matrix.
! References:
!   None
!-------------------------------------------------------------------------------
    use matrix_operations, only: &
      get_lower_triangular_matrix, & ! Procedure(s)
      set_lower_triangular_matrix_dp

    use clubb_precision, only: &
      dp, & ! double precision
      core_rknd

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      d_variables, & ! Total variates
      index1, index2 ! Index of the 2 variates

    real( kind = core_rknd ), dimension(d_variables), intent(in) :: &
      sigma2_on_mu2_ip_array  ! sigma_^2 / mu^2 (ith PDF comp.) in-precip.  [-]

    real( kind = core_rknd ), dimension(d_variables,d_variables), intent(in) :: &
      corr_array ! Array of correlations        [-]

    ! Input/Output Variables
    real( kind = dp ), dimension(d_variables,d_variables), intent(inout) :: &
      corr_stw_matrix ! Correlation matrix      [-]

    ! Local Variables
    real( kind = core_rknd ) :: &
      corr_xy, & ! Correlation between two variates
      covar_xy   ! Lognormal covariance (not a dimensional covariance)

    ! ---- Begin Code ----

    ! Lognormal covariance between two lognormal variates (e.g. rrain and Nr )

    call get_lower_triangular_matrix &
         ( d_variables, index1, index2, corr_array, & ! In
           corr_xy ) ! Out

    if ( corr_xy /= 0._core_rknd ) then
      call construct_LN_LN_element &
           ( corr_xy, sigma2_on_mu2_ip_array(index1), sigma2_on_mu2_ip_array(index2), & ! In
             covar_xy ) ! Out
    else
      covar_xy = 0._core_rknd
    end if

    call set_lower_triangular_matrix_dp &
         ( d_variables, index1, index2, real(covar_xy, kind = dp), & ! In
           corr_stw_matrix ) ! In/out

    return
  end subroutine add_corr_to_matrix_LN_LN

!-------------------------------------------------------------------------------
  subroutine add_corr_to_matrix_gaus_LN( d_variables, iiPDF_s_mellor, &
                                         iiPDF_t_mellor, iiPDF_w, index1, &
                                         sigma2_on_mu2_ip_array, corr_array, &
                                         corr_stw_matrix )
! Description:
!   Add a correlation between s,t Mellor, w and a lognormal variate to a
!   correlation matrix.
! References:
!   None
!-------------------------------------------------------------------------------

    use matrix_operations, only: &
      get_lower_triangular_matrix, & ! Procedure(s)
      set_lower_triangular_matrix_dp

    use clubb_precision, only: &
      dp, & ! double precision
      core_rknd

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      d_variables, & ! Total variates
      iiPDF_s_mellor, & ! Index of s_mellor
      iiPDF_t_mellor, & ! Index of t_mellor
      iiPDF_w, &        ! Index of w (vertical velocity)
      index1           ! Index of the lognormal variate

    real( kind = core_rknd ), dimension(d_variables), intent(in) :: &
      sigma2_on_mu2_ip_array  ! sigma_^2 / mu^2 (ith PDF comp.) in-precip.  [-]

    real( kind = core_rknd ), dimension(d_variables,d_variables), intent(in) :: &
      corr_array ! Array of correlations        [-]

    ! Input/Output Variables
    real( kind = dp ), dimension(d_variables,d_variables), intent(inout) :: &
      corr_stw_matrix ! Correlation matrix      [-]

    ! Local Variables
    real( kind = core_rknd ) :: &
      corr_sx, &  ! Correlation between s and a lognormal variate
      corr_wx, &  ! Correlation between w and a lognormal variate
      covar_sx, & ! Lognormal covariance of s_mellor and x
      covar_wx    ! Lognormal covariance of w and x

    real( kind = dp ) :: &
      covar_tx    ! Lognormal covariance of t_mellor and x

    ! ---- Begin Code ----

    ! Correlations involving s and lognormal variate x

    call get_lower_triangular_matrix &
         ( d_variables, iiPDF_s_mellor, index1, corr_array, & ! In
           corr_sx ) ! Out

    if ( corr_sx /= 0._core_rknd ) then
      ! Covariance between s and lognormal variate x
      ! The variable x could be rrain, Nr, Ncn, et cetera.
      call construct_gaus_LN_element &
           ( corr_sx, 1.0_core_rknd, sigma2_on_mu2_ip_array(index1), & ! In
             covar_sx ) ! Out
    else
      covar_sx = 0._core_rknd
    end if

    call set_lower_triangular_matrix_dp &
         ( d_variables, iiPDF_s_mellor, index1, real(covar_sx, kind = dp), & ! In
           corr_stw_matrix ) ! In/out

    if ( corr_sx /= 0._core_rknd ) then
      ! Approximate the covariance of t and x
      ! This formula relies on the fact that iiPDF_s_mellor < iiPDF_t_mellor
      covar_tx = corr_stw_matrix(iiPDF_t_mellor,iiPDF_s_mellor) * real(covar_sx, kind = dp)
    else
      covar_tx = 0._dp
    end if

    call set_lower_triangular_matrix_dp &
         ( d_variables, iiPDF_t_mellor, index1, covar_tx, & ! In
           corr_stw_matrix ) ! In/out

    ! Correlations involving w and lognormal variate x

    call get_lower_triangular_matrix &
         ( d_variables, iiPDF_w, index1, corr_array, & ! In
           corr_wx ) ! Out

    if ( corr_wx /= 0._core_rknd ) then
      ! Covariance between w and lognormal variate x
      ! The variable x could be rrain, Nr, Ncn, et cetera.
      call construct_gaus_LN_element &
           ( corr_wx, 1.0_core_rknd, sigma2_on_mu2_ip_array(index1), & ! In
             covar_wx ) ! Out
    else
      covar_wx = 0._core_rknd
    end if

    call set_lower_triangular_matrix_dp &
         ( d_variables, iiPDF_w, index1, real(covar_wx, kind = dp), & ! In
           corr_stw_matrix ) ! In/out

    return
  end subroutine add_corr_to_matrix_gaus_LN

!-------------------------------------------------------------------------------
  subroutine add_mu_element_LN( d_variables, index1, xm, sigma2_on_mu2_ip, mu1, mu2 )

! Description:
!   Compute an element of mu1 and mu2 for a lognormal variate.

! References:
!   None
!-------------------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd, & ! variable(s)
      dp

    implicit none

    ! External

    intrinsic :: real

    ! Input Variables

    integer, intent(in) :: &
      d_variables, & ! Number of variates
      index1         ! Index of x in mu1 and mu1

    real( kind = dp ), intent(in) :: &
      Xm ! Mean X  [kg/kg or #/kg]

    real( kind = core_rknd ), dimension(d_variables), intent(in) :: & 
      sigma2_on_mu2_ip ! sigma_^2 / mu^2 (ith PDF comp.) in-precip.  [-]

    ! Input /Output Variables

    real( kind = core_rknd ), dimension(d_variables), intent(inout) :: &
      mu1, mu2 ! Mu 1 and 2     [-]

    ! Local Variables

    real( kind = dp ) :: &
      sigma2_on_mu2_ip_element, & ! sigma_^2 / mu^2 (ith PDF comp.) ip  [-]
      X1, X2  ! PDF parameter for mean of plume 1 and 2.

    ! ---- Begin Code ----

    sigma2_on_mu2_ip_element = real( sigma2_on_mu2_ip(index1), kind = dp )

    call log_sqd_normalized( Xm, sigma2_on_mu2_ip_element, & ! In
                             X1, X2 ) ! Out

    mu1(index1) = real( X1, kind = core_rknd )
    mu2(index1) = real( X2, kind = core_rknd )

    return
  end subroutine add_mu_element_LN

  !-----------------------------------------------------------------------------
  pure function corr_LN_to_covar_gaus( corr_xy, sigma_x_gaus, sigma_y_gaus ) &
    result( covar_xy_gaus )

  ! Description:

  ! References:
  !-----------------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: sqrt, exp, log

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      corr_xy,      & ! Correlation of x and y    [-]
      sigma_x_gaus, & ! Normalized std dev of first term 'x' [-]
      sigma_y_gaus    ! Normalized std dev second term 'y'   [-]

    real( kind = core_rknd ) :: covar_xy_gaus ! Covariance for a gaussian dist. [-]

    ! ---- Begin Code ----

    covar_xy_gaus = log( 1.0_core_rknd + corr_xy * sqrt( exp( sigma_x_gaus**2 ) - 1.0_core_rknd ) &
                                     * sqrt( exp( sigma_y_gaus**2 ) - 1.0_core_rknd ) &
                     )

    return
  end function corr_LN_to_covar_gaus

  !-----------------------------------------------------------------------------
  pure function corr_gaus_LN_to_covar_gaus( corr_sy, sigma_s, sigma_y_gaus ) &
    result( covar_sy_gaus )
  ! Description:

  ! References:

  !-----------------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: sqrt, exp

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      corr_sy,     & ! Correlation of s and y    [-]
      sigma_s,     & ! Normalized std dev of first term (usually Gaussian 's') [units vary]
      sigma_y_gaus   ! Normalized std dev second term 'y'   [-]

    real( kind = core_rknd ) :: covar_sy_gaus ! Covariance for a gaussian dist. [units vary]

    ! ---- Begin Code ----

    covar_sy_gaus = corr_sy * sigma_s * sqrt( exp( sigma_y_gaus**2 ) - 1.0_core_rknd )

    return
  end function corr_gaus_LN_to_covar_gaus

  !-----------------------------------------------------------------------------
  pure function mu_LN_to_mu_gaus( mu, sigma2_on_mu2 ) &
    result( mu_gaus )

  ! Description:
  !
  ! References:
  ! 
  !-----------------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: sqrt, log

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu, &         ! Mean term 'x'           [-]
      sigma2_on_mu2 ! sigma_x^2 over mu_x^2   [-]

    real( kind = core_rknd ) :: mu_gaus ! Mean field converted to gaussian  [-]

    ! ---- Begin Code ----

    mu_gaus = log( mu / sqrt( 1.0_core_rknd + ( sigma2_on_mu2 ) ) )

    return
  end function mu_LN_to_mu_gaus

  !-----------------------------------------------------------------------------
  pure function sigma_LN_to_sigma_gaus( sigma2_on_mu2 ) result( sigma_gaus )

  ! Description:
  !
  ! References:
  ! 
  !-----------------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: sqrt, log

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      sigma2_on_mu2 ! sigma_x^2 over mu_x^2   [-]

    real( kind = core_rknd ) :: sigma_gaus ! Sigma converted to gaussian dist. [-]

    ! ---- Begin Code ----

    sigma_gaus = sqrt( log( 1.0_core_rknd + ( sigma2_on_mu2 ) ) )

    return
  end function sigma_LN_to_sigma_gaus

!-----------------------------------------------------------------------------

end module generate_lh_sample_module
