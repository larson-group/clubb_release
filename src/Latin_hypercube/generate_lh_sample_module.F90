!$Id$
module generate_lh_sample_module

  implicit none

  public :: generate_lh_sample, generate_uniform_sample

  private :: sample_points, gaus_mixt_points, & 
             truncate_gaus_mixt, ltqnorm, gaus_condt, & 
             st_2_rtthl, log_sqd_normalized, choose_permuted_random


  private ! Default scope

  contains

!-------------------------------------------------------------------------------
  subroutine generate_lh_sample &
             ( n_micro_calls, d_variables, hydromet_dim, & 
               cloud_frac, wm, rtm, thlm, pdf_params, level, & 
               hydromet, correlation_array, X_u_one_lev, &
               X_mixt_comp_one_lev, &
               LH_rt, LH_thl, X_nl_one_lev )
! Description:
!   This subroutine generates a Latin Hypercube sample.

! References:
!   ``Supplying Local Microphysical Parameterizations with Information about
!     Subgrid Variability: Latin Hypercube Sampling'', JAS Vol. 62,
!     p. 4010--4026, Larson, et al. 2005.
!-------------------------------------------------------------------------------

    use KK_microphys_module, only: &
      corr_LN_to_cov_gaus, & ! Procedure(s)
      corr_gaus_LN_to_cov_gaus, &
      sigma_LN_to_sigma_gaus

    use constants, only:  &
      max_mag_correlation, &  ! Constant
!     s_mellor_tol,  &  ! s tolerance in kg/kg
      rttol, &      ! rt tolerance in kg/kg
      thltol, &     ! thetal tolerance in K
      wtol_sqd, &   ! w^2 tolerance in m^2/s^2
      rr_tol, &     ! rr tolerance in kg/kg
      Nr_tol, &     ! Nr tolerance in #/kg
      Nc_tol        ! Nc tolerance in #/kg

    use variables_prognostic_module, only:  &
      pdf_parameter  ! type

    use array_index, only: &
      iiNcm,    & ! Variables
      iiNrm,    &
      iirrainm, &
      iiLH_rrain, &
      iiLH_Nr, &
      iiLH_Nc, &
      iiLH_rt, &
      iiLH_thl, &
      iiLH_w

    use mt95, only: genrand_real ! Constants

    implicit none

    ! External
    intrinsic :: dble, min, max, sqrt

    ! Constant Parameters
    logical, parameter :: &
      l_sample_out_of_cloud = .true.

    ! Tolerance on the standard deviation of s for latin hypercube only
    real, parameter :: &
      LH_stdev_s_tol = 1e-6  ! [kg/kg]

    ! Input Variables
    integer, intent(in) :: &
      n_micro_calls, & ! `n' Number of calls to microphysics (normally=2)
      d_variables,   & ! `d' Number of variates (normally 3 + microphysics specific variables)
      hydromet_dim     ! Number of hydrometeor species

    real, dimension(hydromet_dim), intent(in) :: &
      hydromet ! Hydrometeor species [units vary]

    real, intent(in) :: &
      cloud_frac, & ! Cloud fraction, 0 <= cloud_frac <= 1
      wm,         & ! Vertical velocity                   [m/s]
      rtm,        & ! Mean total water mixing ratio       [kg/kg]
      thlm          ! Mean liquid potential temperature   [K]

    type(pdf_parameter), intent(in) :: &
      pdf_params ! PDF parameters output by closure_new [units vary]

    integer, intent(in) :: level  ! Level info. for PDF parameters.

    ! From the KK_microphys_module
    real, dimension(d_variables,d_variables), intent(in) :: &
      correlation_array ! Correlations for sampled variables    [-]

    real(kind=genrand_real), intent(in), dimension(n_micro_calls,d_variables+1) :: &
      X_u_one_lev ! Sample drawn from uniform distribution from a particular grid level

    integer, intent(in), dimension(n_micro_calls) :: &
      X_mixt_comp_one_lev ! Whether we're in the 1st or 2nd mixture component

    ! Output Variables
    real, intent(out), dimension(n_micro_calls) :: &
      LH_rt, & ! Total water mixing ratio          [kg/kg]
      LH_thl   ! Liquid potential temperature      [K]

    double precision, intent(out), dimension(n_micro_calls,d_variables) :: &
      X_nl_one_lev ! Sample that is transformed ultimately to normal-lognormal

    ! A true/false flag that determines whether
    ! the PDF allows us to construct a sample
    logical :: l_sample_flag

    ! Local Variables

    logical, dimension(d_variables) :: &
      l_d_variable_lognormal ! Whether a given variable in X_nl has a lognormal dist.

    real :: &
      mixt_frac,   & ! Weight of 1st normal distribution (Sk_w dependent)      [-]
      w1,          & ! Mean of w for 1st normal distribution                 [m/s]
      w2,          & ! Mean of w for 2nd normal distribution                 [m/s]
      varnce_w1,   & ! Variance of w for 1st normal distribution         [m^2/s^2]
      varnce_w2,   & ! Variance of w for 2nd normal distribution         [m^2/s^2]
      thl1,        & ! Mean of th_l for 1st normal distribution                [K]
      thl2,        & ! Mean of th_l for 2nd normal distribution                [K]
      varnce_thl1, & ! Variance of th_l for 1st normal distribution          [K^2]
      varnce_thl2, & ! Variance of th_l for 2nd normal distribution          [K^2]
      rt1,         & ! Mean of r_t for 1st normal distribution             [kg/kg]
      rt2,         & ! Mean of r_t for 2nd normal distribution             [kg/kg]
      varnce_rt1,  & ! Variance of r_t for 1st normal distribution     [kg^2/kg^2]
      varnce_rt2,  & ! Variance of r_t for 2nd normal distribution     [kg^2/kg^2]
      s1,          & ! Mean of s for 1st normal distribution               [kg/kg]
      s2,          & ! Mean of s for 2nd normal distribution               [kg/kg]
      stdev_s1,    & ! Standard deviation of s for 1st normal distribution [kg/kg]
      stdev_s2,    & ! Standard deviation of s for 2nd normal distribution [kg/kg]
      crt1,        & ! Coefficient for s'                                      [-]
      crt2,        & ! Coefficient for s'                                      [-]
      cthl1,       & ! Coefficient for s'                                    [1/K]
      cthl2          ! Coefficient for s'                                    [1/K]

    real :: &
      cloud_frac1, & ! Cloud fraction for 1st normal distribution              [-]
      cloud_frac2    ! Cloud fraction for 2nd normal distribution              [-]

    ! sub-plume correlation coefficient between rt, thl
    ! varies between -1 < rrtthl < 1
    real :: rrtthl

    ! Clip the magnitude of the correlation between rt and thl
    real :: rrtthl_reduced

    double precision :: rrtthl_reduced1, rrtthl_reduced2

    ! Means of s, t, w, & hydrometeors for plumes 1 and 2
    real, dimension(d_variables) :: &
      mu1, mu2

    ! Covariance (not correlation) matrix of rt, thl, w for plumes 1 and 2
    double precision, dimension(3,3) :: &
      Sigma_rtthlw_1,  & 
      Sigma_rtthlw_2

    ! Columns of Sigma_stw, X_nl_one_lev:  1   2   3   4 ... d_variables
    !                                      s   t   w   hydrometeors
    double precision, dimension(d_variables,d_variables) :: &
      Sigma_stw_1, & ! Covariance of s,t, w + hydrometeors for plume 1
      Sigma_stw_2    ! Covariance of s,t, w + hydrometeors for plume 2

    double precision :: &
      Ncm,     & ! Cloud droplet number concentration.[number / kg air]
      Nc1,     & ! PDF parameter for mean of plume 1. [#/kg]
      Nc2,     & ! PDF parameter for mean of plume 2. [#/kg]
      var_Nc1, & ! PDF param for width of plume 1.    [(#/kg)^2]
      var_Nc2, & ! PDF param for width of plume 2.    [(#/kg^2]
      Nrm,     & ! Rain droplet number concentration. [number / kg air]
      Nr1,     & ! PDF parameter for mean of plume 1. [#/kg]
      Nr2,     & ! PDF parameter for mean of plume 2. [#/kg]
      var_Nr1, & ! PDF param for width of plume 1.    [(#/kg)^2]
      var_Nr2    ! PDF param for width of plume 2.    [(#/kg^2]

    real :: corr_rrNr, covar_rrNr1, covar_rrNr2, corr_srr, corr_sNr, &
            covar_sNr1, covar_sNr2, covar_srr1, covar_srr2

    double precision :: covar_trr1, covar_trr2, covar_tNr2, covar_tNr1

!   real :: & 
!     stdev_Nc, & ! Standard deviation of Nc   [#/kg]
!     corr_tNc, & ! Correlation between t and Nc [-]
!     corr_sNc, & ! Correlation between s and Nc [-]
!     covar_tNc1,    & ! Covariance of t and Nc1      []
!     covar_tNc2,    & ! Covariance of t and Nc2      []
!     covar_sNc1,    & ! Covariance of s and Nc1      [kg^2/kg^2]
!     covar_sNc2       ! Covariance of s and Nc2      [kg^2/kg^2]

!   double precision, dimension(2,2) :: corr_st_mellor_1, corr_st_mellor_2

    real :: stdev_rr, stdev_Nr

    ! rr = specific rain content. [rr] = kg rain / kg air
    double precision :: &
      rrainm, &  ! rain water mixing ratio         [kg/kg]
      rr1, &  ! PDF parameter for mean of plume 1. [kg/kg]
      rr2, &  ! PDF parameter for mean of plume 2. [kg/kg]
      var_rr1, & ! PDF param for width of plume 1     [(kg/kg)^2]
      var_rr2    ! PDF param for width of plume 2.    [(kg/kg)^2]

    double precision :: &
      Ncp2_on_Ncm2, & ! = Ncp2 divided by Ncm^2    [-]
      Nrp2_on_Nrm2, & ! = Nrp2 divided by Nrm^2    [-]
      rrp2_on_rrainm2 ! = rrp2 divided by rrainm^2 [-]

    integer :: i, iiLH_s_mellor, iiLH_t_mellor

    ! ---- Begin Code ----
    ! Determine which variables are a lognormal distribution
    i = max( iiLH_rt, iiLH_thl, iiLH_w )
    l_d_variable_lognormal(1:i) = .false. ! The 1st 3 variates
    l_d_variable_lognormal(i+1:d_variables) = .true.  ! Hydrometeors

    ! Input pdf parameters.

    if ( pdf_params%varnce_w1(level) > wtol_sqd ) then
      varnce_w1 = pdf_params%varnce_w1(level)
      w1 = pdf_params%w1(level)
    else
      varnce_w1 = wtol_sqd
      w1 = wm
    end if
    if ( pdf_params%varnce_w2(level) > wtol_sqd ) then
      varnce_w2 = pdf_params%varnce_w2(level)
      w2 = pdf_params%w2(level)
    else
      varnce_w2 = wtol_sqd
      w2 = wm
    end if
    if ( pdf_params%varnce_rt1(level) > rttol**2 ) then
      varnce_rt1 = pdf_params%varnce_rt1(level)
      rt1 = pdf_params%rt1(level)
    else
      varnce_rt1 = rttol**2
      rt1 = rtm
    end if
    if ( pdf_params%varnce_rt2(level) > rttol**2 ) then
      varnce_rt2 = pdf_params%varnce_rt2(level)
      rt2 = pdf_params%rt2(level)
    else
      varnce_rt2 = rttol**2
      rt2 = rtm
    end if
    if ( pdf_params%varnce_thl1(level) > thltol**2 ) then 
      varnce_thl1 = pdf_params%varnce_thl1(level)
      thl1  = pdf_params%thl1(level)
    else
      varnce_thl1 = thltol**2
      thl1  = thlm
    end if
    if ( pdf_params%varnce_thl2(level) > thltol**2 ) then
      varnce_thl2 = pdf_params%varnce_thl2(level)
      thl2  = pdf_params%thl2(level)
    else
      varnce_thl2 = thltol**2
      thl2  = thlm
    end if
    if ( pdf_params%stdev_s1(level) > LH_stdev_s_tol ) then
      stdev_s1 = pdf_params%stdev_s1(level)
      s1       = pdf_params%s1(level)
    else
      ! Use a larger value than s_mellor_tol, for reasons of numerical stability
      stdev_s1 = LH_stdev_s_tol
      s1       = pdf_params%s1(level) * pdf_params%mixt_frac(level) &
               + (1.0-pdf_params%mixt_frac(level)) * pdf_params%s2(level)
    end if
    if ( pdf_params%stdev_s2(level) > LH_stdev_s_tol ) then
      stdev_s2 = pdf_params%stdev_s2(level)
      s2       = pdf_params%s2(level)
    else
      ! Use a larger value than s_mellor_tol, as above.
      stdev_s2 = LH_stdev_s_tol
      s2       = pdf_params%s1(level) * pdf_params%mixt_frac(level) &
               + (1.0-pdf_params%mixt_frac(level)) * pdf_params%s2(level)
    end if

    crt1 = pdf_params%crt1(level)
    crt2 = pdf_params%crt2(level)
    cthl1 = pdf_params%cthl1(level)
    cthl2 = pdf_params%cthl2(level)

    mixt_frac   = pdf_params%mixt_frac(level)
    cloud_frac1 = pdf_params%cloud_frac1(level)
    cloud_frac2 = pdf_params%cloud_frac2(level)
    rrtthl      = pdf_params%rrtthl(level)

    !---------------------------------------------------------------------------
    ! Generate a set of sample points for a microphysics scheme
    !---------------------------------------------------------------------------

    ! We prognose rt-thl-w,
    !    but we set means, covariance of hydrometeors (e.g. rrain, Nc) to constants.

    l_sample_flag = .true.
    if ( l_sample_out_of_cloud ) then
      ! Sample non-cloudy grid boxes as well -dschanen 3 June 2009
      cloud_frac1 = 1.0
      cloud_frac2 = 1.0

    else if ( cloud_frac < 0.001 ) then
      ! In this case there are essentially no cloudy points to sample;
      ! Set sample points to zero.
      X_nl_one_lev(:,:)   = -999.0
      l_sample_flag = .false.

    end if ! l_sample_out_of_cloud

    ! Standard sample for testing purposes when n=2
    ! X_u_one_lev(1,1:(d+1)) = ( / 0.0001d0, 0.46711825945881d0, &
    !             0.58015016959859d0, 0.61894015386778d0, 0.1d0, 0.1d0  / )
    ! X_u_one_lev(2,1:(d+1)) = ( / 0.999d0, 0.63222458307464d0, &
    !             0.43642762850981d0, 0.32291562498749d0, 0.1d0, 0.1d0  / )

    if ( l_sample_flag ) then


      ! Compute PDF parameters for Nc, rr.
      ! Assume that Nc, rr obey single-lognormal distributions

      ! Nc  = droplet number concentration.  [Nc] = number / kg air
      ! Ncm  = mean of Nc
      ! Ncp2_on_Ncm2 = variance of Nc divided by Ncm^2
      !  We must have a Ncp2_on_Ncm2 >= machine epsilon
      ! Nc1  = PDF parameter for mean of plume 1. [Nc1] = (#/kg)
      ! Nc2  = PDF parameter for mean of plume 2. [Nc2] = (#/kg)
      ! var_Nc1,2 = PDF param for width of plume 1,2. [var_Nc1,2] = (#/kg)**2

      if ( iiLH_Nc > 0 ) then 
        Ncm = dble( hydromet(iiNcm) )
        Ncp2_on_Ncm2 = dble( correlation_array(iiLH_Nc,iiLH_Nc) )

        call log_sqd_normalized( Ncm, Ncp2_on_Ncm2, dble( Nc_tol ), & ! In
                                 Nc1, Nc2, var_Nc1, var_Nc2 ) ! Out
      end if

      ! rr = specific rain content. [rr] = kg rain / kg air
      ! rrainm  = mean of rr; rrp2 = variance of rr, must have rrp2>0.
      ! rr1  = PDF parameter for mean of plume 1. [rr1] = (kg/kg)
      ! rr2  = PDF parameter for mean of plume 2. [rr2] = (kg/kg)
      ! var_rr1,2 = PDF param for width of plume 1,2. [var_rr1,2] = (kg/kg)**2

      if ( iiLH_rrain > 0 ) then 
        rrainm = dble( hydromet(iirrainm) )
        rrp2_on_rrainm2 = dble( correlation_array(iiLH_rrain,iiLH_rrain) )
        call log_sqd_normalized( rrainm, rrp2_on_rrainm2, dble( rr_tol ), & ! In
                                 rr1, rr2, var_rr1, var_rr2 ) ! Out
      end if

      if ( iiLH_Nr > 0 ) then 
        Nrm = dble( hydromet(iiNrm) )
        Nrp2_on_Nrm2 = dble( correlation_array(iiLH_Nr,iiLH_Nr) )

        call log_sqd_normalized( Nrm, Nrp2_on_Nrm2, dble( Nr_tol ), & ! In
                                 Nr1, Nr2, var_Nr1, var_Nr2 ) ! Out
      end if

      ! Means of s, t, w, Nc, Nr, rr for Gaussians 1 and 2

      mu1((/iiLH_rt,iiLH_thl,iiLH_w/)) &
        = (/ s1, 0., w1 /)
      mu2((/iiLH_rt,iiLH_thl,iiLH_w/)) &
        = (/ s2, 0., w2 /)

      if ( iiLH_rrain > 0 ) then
        mu1(iiLH_rrain) = rr1
        mu2(iiLH_rrain) = rr2
      end if

      if ( iiLH_Nc > 0 ) then
        mu1(iiLH_Nc) = Nc1
        mu2(iiLH_Nc) = Nc2
      end if

      if ( iiLH_Nr > 0 ) then
        mu1(iiLH_Nr) = Nr1
        mu2(iiLH_Nr) = Nr2
      end if

      ! An old subroutine, gaus_rotate, couldn't handle large correlations;
      !   I assume the replacement, gaus_condt, has equal trouble.
      !   Therefore we input smaller correlations
      ! max_mag_correlation = 0.99 in constants.F90
      rrtthl_reduced = min( max_mag_correlation, max( rrtthl, -max_mag_correlation ) )

      ! Within-plume rt-thl correlation terms with rt in kg/kg
      rrtthl_reduced1 = dble( rrtthl_reduced*sqrt( varnce_rt1*varnce_thl1 ) )
      rrtthl_reduced2 = dble( rrtthl_reduced*sqrt( varnce_rt2*varnce_thl2 ) )

      ! Covariance (not correlation) matrices of rt-thl-w
      !    for Gaussians 1 and 2
      ! For now, assume no within-plume correlation of w with
      !    any other variables.

      ! Sigma_rtthlw_1,2
      Sigma_rtthlw_1 = 0.d0 ! Start with no covariance, and add matrix elements
      Sigma_rtthlw_2 = 0.d0

      ! Sigma_stw_1,2
      Sigma_stw_1 = 0.d0 ! Start with no covariance, and add matrix elements
      Sigma_stw_2 = 0.d0

      Sigma_rtthlw_1(iiLH_rt,(/iiLH_rt,iiLH_thl/))  = (/ dble( varnce_rt1 ), rrtthl_reduced1 /)

      Sigma_rtthlw_2(iiLH_rt,(/iiLH_rt,iiLH_thl/))  = (/ dble( varnce_rt2 ), rrtthl_reduced2 /)

      Sigma_rtthlw_1(iiLH_thl,(/iiLH_rt,iiLH_thl/)) = (/ rrtthl_reduced1, dble( varnce_thl1 ) /)

      Sigma_rtthlw_2(iiLH_thl,(/iiLH_rt,iiLH_thl/)) = (/ rrtthl_reduced2, dble( varnce_thl2 ) /)

      Sigma_rtthlw_1(iiLH_w,iiLH_w) = dble( varnce_w1 )

      Sigma_rtthlw_2(iiLH_w,iiLH_w) = dble( varnce_w2 )

      ! Convert each Gaussian from rt-thl-w variables to s-t-w vars.
      call rtpthlp_2_sptp( 3, Sigma_rtthlw_1(1:3,1:3), dble( crt1 ), dble( cthl1 ), & ! In
                           Sigma_stw_1(1:3,1:3) ) ! Out
      call rtpthlp_2_sptp( 3, Sigma_rtthlw_2(1:3,1:3), dble( crt2 ), dble( cthl2 ), & ! In
                           Sigma_stw_2(1:3,1:3) ) ! Out

      ! The s and t elements in Sigma_stw correspond to the rt thl element
      !   in Sigma_rtthlw
      iiLH_s_mellor = iiLH_rt 
      iiLH_t_mellor = iiLH_thl

      ! Determine the correlation of s and t for the purposes of approximating
      ! the correlation of t and the other samples
!     call covar_matrix_2_corr_matrix( 2, Sigma_stw_1(1:2,1:2), corr_st_mellor_1 )
!     call covar_matrix_2_corr_matrix( 2, Sigma_stw_2(1:2,1:2), corr_st_mellor_2 )

!     stdev_t1 = sqrt( real( Sigma_stw_1(iiLH_t_mellor,iiLH_t_mellor) ) ) ! Std dev of t (1st plume)
!     stdev_t2 = sqrt( real( Sigma_stw_2(iiLH_t_mellor,iiLH_t_mellor) ) ) ! Std dev of t (2nd plume)

      if ( iiLH_Nc > 0 ) then
        Sigma_stw_1(iiLH_Nc,iiLH_Nc) = var_Nc1
        Sigma_stw_2(iiLH_Nc,iiLH_Nc) = var_Nc2
      end if

      if ( iiLH_Nr > 0 ) then
        Sigma_stw_1(iiLH_Nr,iiLH_Nr) = var_Nr1
        Sigma_stw_2(iiLH_Nr,iiLH_Nr) = var_Nr2
      end if

      if ( iiLH_rrain > 0 ) then
        Sigma_stw_1(iiLH_rrain,iiLH_rrain) = var_rr1
        Sigma_stw_2(iiLH_rrain,iiLH_rrain) = var_rr2
      end if

      if ( iiLH_rrain > 0 .and. iiLH_Nr > 0 ) then 
        ! Compute standard deviation of Nr & rrain
        stdev_rr = real( rrainm ) * sqrt( correlation_array(iiLH_rrain,iiLH_rrain) )
        stdev_Nr = real( Nrm ) * sqrt( correlation_array(iiLH_Nr,iiLH_Nr) )

        if ( rrainm > dble( rr_tol ) .and. Nrm > dble( Nr_tol ) ) then

          corr_rrNr = correlation_array(iiLH_rrain,iiLH_Nr)

          ! Covariance between rain water mixing ratio rain number concentration
          covar_rrNr1 = corr_LN_to_cov_gaus &
                   ( corr_rrNr, &
                     sigma_LN_to_sigma_gaus( real( stdev_rr ), real( rrainm ) ), &
                     sigma_LN_to_sigma_gaus( real( stdev_Nr ), real( Nrm ) ) )
          covar_rrNr2 = covar_rrNr1

          Sigma_stw_1(iiLH_rrain,iiLH_Nr) = dble( covar_rrNr1 )
          Sigma_stw_1(iiLH_Nr,iiLH_rrain) = dble( covar_rrNr1 )
          Sigma_stw_2(iiLH_rrain,iiLH_Nr) = dble( covar_rrNr2 )
          Sigma_stw_2(iiLH_Nr,iiLH_rrain) = dble( covar_rrNr2 )
        end if

        ! Covariances involving s and Nr & rr
        if ( stdev_s1 > LH_stdev_s_tol .and. Nrm > dble( Nr_tol ) ) then
          corr_sNr = correlation_array(iiLH_s_mellor,iiLH_Nr)

          ! Covariance between s and rain number conc.
          covar_sNr1 = corr_gaus_LN_to_cov_gaus &
                   ( corr_sNr, &
                     stdev_s1, &
                     sigma_LN_to_sigma_gaus( real( stdev_Nr ), real( Nrm ) ) )

          Sigma_stw_1(iiLH_s_mellor,iiLH_Nr) = dble( covar_sNr1 )
          Sigma_stw_1(iiLH_Nr,iiLH_s_mellor) = dble( covar_sNr1 )

          ! Approximate the covariance of t and Nr
          covar_tNr1 = ( Sigma_stw_1(iiLH_t_mellor,iiLH_s_mellor) &
            * dble( covar_sNr1 ) ) / dble( stdev_s1 )**2

          Sigma_stw_1(iiLH_t_mellor,iiLH_Nr) = dble( covar_tNr1 )
          Sigma_stw_1(iiLH_Nr,iiLH_t_mellor) = dble( covar_tNr1 )
        end if

        if ( stdev_s2 > LH_stdev_s_tol .and. Nrm > dble( Nr_tol ) ) then

          corr_sNr = correlation_array(iiLH_s_mellor,iiLH_Nr)

          covar_sNr2 = corr_gaus_LN_to_cov_gaus &
                   ( corr_sNr, &
                     stdev_s2, &
                     sigma_LN_to_sigma_gaus( real( stdev_Nr ), real( Nrm ) ) )

          Sigma_stw_2(iiLH_s_mellor,iiLH_Nr) = dble( covar_sNr2 )
          Sigma_stw_2(iiLH_Nr,iiLH_s_mellor) = dble( covar_sNr2 )

          ! Approximate the covariance of t and Nr
          covar_tNr2 = ( Sigma_stw_2(iiLH_t_mellor,iiLH_s_mellor) &
            * dble( covar_sNr2 ) ) / dble( stdev_s2 )**2

          Sigma_stw_2(iiLH_t_mellor,iiLH_Nr) = dble( covar_tNr2 )
          Sigma_stw_2(iiLH_Nr,iiLH_t_mellor) = dble( covar_tNr2 )
        end if

        if ( stdev_s1 > LH_stdev_s_tol .and. rrainm > dble( rr_tol ) ) then

          corr_srr = correlation_array(iiLH_s_mellor,iiLH_rrain)

          ! Covariance between s and rain water mixing ratio
          covar_srr1 = corr_gaus_LN_to_cov_gaus &
                   ( corr_srr, &
                     stdev_s1, &
                     sigma_LN_to_sigma_gaus( real( stdev_rr ), real( rrainm ) ) )

          Sigma_stw_1(iiLH_s_mellor,iiLH_rrain) = dble( covar_srr1 )
          Sigma_stw_1(iiLH_rrain,iiLH_s_mellor) = dble( covar_srr1 )

          ! Approximate the covariance of t and rr
          covar_trr1 = ( Sigma_stw_1(iiLH_t_mellor,iiLH_s_mellor) &
            * dble( covar_srr1 ) ) / dble( stdev_s1 )**2

          Sigma_stw_1(iiLH_t_mellor,iiLH_rrain) = dble( covar_trr1 )
          Sigma_stw_1(iiLH_rrain,iiLH_t_mellor) = dble( covar_trr1 )
        end if

        if ( stdev_s2 > LH_stdev_s_tol .and. rrainm > dble( rr_tol ) ) then

          corr_srr = correlation_array(iiLH_s_mellor,iiLH_rrain)

          covar_srr2 = corr_gaus_LN_to_cov_gaus &
                   ( corr_srr, &
                     stdev_s2, &
                     sigma_LN_to_sigma_gaus( real( stdev_rr ), real( rrainm ) ) )

          Sigma_stw_2(iiLH_s_mellor,iiLH_rrain) = dble( covar_srr2 )
          Sigma_stw_2(iiLH_rrain,iiLH_s_mellor) = dble( covar_srr2 )

          ! Approximate the covariance of t and rr
          covar_trr2 = ( Sigma_stw_2(iiLH_t_mellor,iiLH_s_mellor) &
            * dble( covar_srr2 ) ) / dble( stdev_s2 )**2

          Sigma_stw_2(iiLH_t_mellor,iiLH_rrain) = dble( covar_trr2 )
          Sigma_stw_2(iiLH_rrain,iiLH_t_mellor) = dble( covar_trr2 )
        end if
      end if ! if iiLH_rrain > 0 .and. iiLH_Nr > 0 

!     if ( iiLH_Nc > 0 ) then

        ! Covariances involving s and Nc (currently disabled)
!       corr_sNc = correlation_array(iiLH_s_mellor,iiLH_Nc)
!       stdev_Nc = real( Ncm ) * sqrt( correlation_array(iiLH_Nc,iiLH_Nc) )

!       if ( stdev_s1 > LH_stdev_s_tol .and. Ncm > dble( Nc_tol ) ) then
!         ! The variable s is already Gaussian
!         stdev_sNc1 = corr_gaus_LN_to_cov_gaus &
!                 ( corr_sNc, &
!                   stdev_s1, &
!                   sigma_LN_to_sigma_gaus( real( stdev_Nc ), real( Ncm ) ) )

!         Sigma_stw_1(iiLH_s_mellor,iiLH_Nc) = dble( stdev_sNc1 )
!         Sigma_stw_1(iiLH_Nc,iiLH_s_mellor) = dble( stdev_sNc1 )

!         ! Approximate the covariance of t and Nc
!         covar_tNc1 = ( Sigma_stw_1(iiLH_t_mellor,iiLH_s_mellor) * covar_sNc1 ) / stdev_s1**2

!         Sigma_stw_1(iiLH_t_mellor,iiLH_Nc) = dble( covar_tNc1 )
!         Sigma_stw_2(iiLH_Nc,iiLH_t_mellor) = dble( covar_tNc2 )

!       end if

!       if ( stdev_s2 > LH_stdev_s_tol .and. Ncm > dble( Nc_tol ) ) then
!         stdev_sNc2 = corr_gaus_LN_to_cov_gaus &
!                 ( corr_sNc, &
!                   stdev_s2, &
!                   sigma_LN_to_sigma_gaus( real( stdev_Nc ), real( Ncm ) ) )

!         Sigma_stw_2(iiLH_s_mellor,iiLH_Nc) = dble( stdev_sNc2 )
!         Sigma_stw_2(iiLH_Nc,iiLH_s_mellor) = dble( stdev_sNc2 )

!         ! Approximate the covariance of t and Nc
!         covar_tNc2 = ( Sigma_stw_2(iiLH_t_mellor,iiLH_s_mellor) * covar_sNc2 ) / stdev_s2**2

!         Sigma_stw_2(iiLH_t_mellor,iiLH_Nc) = dble( stNc2 )
!         Sigma_stw_2(iiLH_Nc,iiLH_t_mellor) = dble( stNc2 )

!       end if

!     end if ! iiLH_Nc > 0 

      call sample_points( n_micro_calls, d_variables, dble( mixt_frac ), &  ! In
                          dble( rt1 ), dble( thl1 ), &  ! In
                          dble( rt2 ), dble( thl2 ), &  ! In
                          dble( crt1 ), dble( cthl1 ), &  ! In
                          dble( crt2 ), dble( cthl2 ), &  ! In
                          mu1, mu2, &  ! In
                          dble( cloud_frac1 ), dble( cloud_frac2 ), & ! In 
                          l_d_variable_lognormal, & ! In
                          X_u_one_lev, & ! In
                          X_mixt_comp_one_lev, & ! In
                          Sigma_stw_1, Sigma_stw_2, & ! In/Out
                          LH_rt, LH_thl, X_nl_one_lev ) ! Out

      ! End of overall if-then statement for Latin hypercube code
    end if

    return
  end subroutine generate_lh_sample

!---------------------------------------------------------------------------------------------------
  subroutine sample_points( n_micro_calls, d_variables, mixt_frac, & 
                            rt1, thl1, rt2, thl2, & 
                            crt1, cthl1, crt2, cthl2, & 
                            mu1, mu2,  & 
                            cloud_frac1, cloud_frac2, & 
                            l_d_variable_lognormal, &
                            X_u_one_lev, &
                            X_mixt_comp_one_lev, &
                            Sigma_stw_1, Sigma_stw_2, & 
                            LH_rt, LH_thl, X_nl_one_lev )

! Description:
!   Generates n random samples from a d-dim Gaussian-mixture PDF.
!   Uses Latin hypercube method.
!   To be called from pdf_closure of CLUBB.

!   We take samples only from the cloudy part of the grid box.
!   We use units of kg/kg.

! References:
!   None
!----------------------------------------------------------------------

    use constants, only:  &
        fstderr  ! Constant(s)

    use array_index, only: &
      iiLH_rt, & ! Variable(s)
      iiLH_thl

    use matrix_operations, only: &
      covar_matrix_2_corr_matrix ! Procedure(s)

    use error_code, only:  &
      clubb_at_least_debug_level  ! Procedure(s)

    implicit none

    ! Input variables
    integer, intent(in) :: &
      n_micro_calls, & ! `n'   Number of calls to microphysics (normally=2)
      d_variables      ! Number of variates (normally=5)

    ! Weight of 1st Gaussian, 0 <= mixt_frac <= 1
    double precision, intent(in) :: mixt_frac

    !rt1, thl1 = mean of rt, thl for Gaus comp 1
    !rt2, thl2 = mean of rt, thl for Gaus comp 2
    double precision, intent(in) :: rt1, thl1, rt2, thl2

    ! Thermodynamic constants for plumes 1 and 2, units of kg/kg
    double precision, intent(in) :: &
      crt1,  & ! coefficient relating rt, s and t for Gaus comp 1
      cthl1, & ! coeff relating thl, s and t for component 1
      crt2,  & ! coefficient relating rt, s and t for component 2
      cthl2    ! coefficient relating thl, s and t for comp. 2

    ! Latin hypercube variables, i.e. s, t, w, etc.
    real, intent(in), dimension(d_variables) :: &
      mu1, mu2 ! d-dimensional column vector of means of 1st, 2nd components

    ! Cloud fractions for components 1 and 2
    double precision, intent(in) :: &
      cloud_frac1, cloud_frac2 ! cloud fraction associated w/ 1st, 2nd mixture component

    logical, intent(in), dimension(d_variables) :: &
      l_d_variable_lognormal ! Whether a given element of X_nl is lognormal

    double precision, intent(in), dimension(n_micro_calls,d_variables+1) :: &
      X_u_one_lev ! Sample drawn from uniform distribution from particular grid level

    integer, intent(in), dimension(n_micro_calls) :: &
      X_mixt_comp_one_lev ! Whether we're in the 1st or 2nd mixture component

    ! Columns of Sigma_stw, X_nl_one_lev:  1   2   3   4 ... d_variables
    !                                      s   t   w   hydrometeors
    double precision, intent(inout), dimension(d_variables,d_variables) :: &
      Sigma_stw_1, &
      Sigma_stw_2

    ! Output Variables
    ! Total water, theta_l: mean plus perturbations
    real, intent(out), dimension(n_micro_calls) :: LH_rt, LH_thl


    double precision, intent(out), dimension(n_micro_calls,d_variables) :: &
      X_nl_one_lev ! Sample that is transformed ultimately to normal-lognormal

    ! Local Variables
    integer :: col, sample

    double precision, dimension(d_variables,d_variables) :: &
      Sigma_stw_1_corr, Sigma_stw_2_corr

    ! Sample of s points that is drawn only from normal distribution
    double precision, dimension(n_micro_calls) :: s_pts

    integer :: iiLH_s_mellor, iiLH_t_mellor
!   integer :: i, j

    ! ---- Begin Code ----

    iiLH_s_mellor = iiLH_rt  ! Mellor's s is at the same index as rt in the Sigma_stw arrays
    iiLH_t_mellor = iiLH_thl ! Mellor's t is at the same index as thl "  "

    if ( clubb_at_least_debug_level( 2 ) ) then

      call covar_matrix_2_corr_matrix( d_variables, Sigma_stw_1, Sigma_stw_1_corr )
      call covar_matrix_2_corr_matrix( d_variables, Sigma_stw_2, Sigma_stw_2_corr )

      if ( any( Sigma_stw_1_corr > 1.0 ) .or. any( Sigma_stw_2_corr < -1.0 ) ) then
        write(fstderr,*) "Sigma_stw_1 has a correlation > 1 or < -1"
      end if
      if ( any( Sigma_stw_1_corr > 1.0 ) .or. any( Sigma_stw_1_corr < -1.0 ) ) then
        write(fstderr,*) "Sigma_stw_2 has a correlation > 1 or < -1"
      end if
      ! Kluge to reduce the covariance when the correlation is too large.
      ! Doesn't help very much. -dschanen 20 Aug 2009
!     do i = 1, d_variables 
!       do j = 1, d_variables 
!         if ( Sigma_stw_1_corr(i,j) > 1.0 .or. Sigma_stw_1_corr(i,j) < -1.0 ) then
!           Sigma_stw_1(i,j) = Sigma_stw_1(i,j) / ( abs( Sigma_stw_1_corr(i,j) ) )
!         end if
!         if ( Sigma_stw_2_corr(i,j) > 1.0 .or. Sigma_stw_2_corr(i,j) < -1.0 ) then
!           Sigma_stw_2(i,j) = Sigma_stw_2(i,j) / ( abs( Sigma_stw_2_corr(i,j) ) )
!         end if
!       end do
!     end do

    end if ! clubb_at_least_debug_level( 2 )

    ! Let s PDF (1st column) be a truncated Gaussian.
    ! Take sample solely from cloud points.
    col = iiLH_s_mellor
    call truncate_gaus_mixt( n_micro_calls, d_variables, col, mixt_frac, mu1, mu2, &  ! In
                             Sigma_stw_1, Sigma_stw_2, cloud_frac1, cloud_frac2, X_u_one_lev, & ! In
                             X_mixt_comp_one_lev, & ! In
                             s_pts ) ! Out

    ! Generate n samples of a d-variate Gaussian mixture
    ! by transforming Latin hypercube points, X_u_one_lev.
    call gaus_mixt_points( n_micro_calls, d_variables, mixt_frac, mu1, mu2, &  ! In
                           Sigma_stw_1, Sigma_stw_2, & ! In
                           cloud_frac1, cloud_frac2, X_u_one_lev, s_pts, & ! In
                           X_mixt_comp_one_lev, & ! In
                           X_nl_one_lev ) ! Out

! Transform s (column 1) and t (column 2) back to rt and thl
! This is only needed if you need rt, thl in your microphysics.
!     call sptp_2_rtpthlp &
!          ( n_micro_calls, d_variables, mixt_frac, crt1, cthl1, crt2, cthl2, &
!            cloud_frac1, cloud_frac2, X_nl_one_lev(1:n_micro_calls,1), &
!            X_nl_one_lev(1:n_micro_calls,2), &
!            X_u_one_lev, rtp, thlp )
      call st_2_rtthl( n_micro_calls, mixt_frac, rt1, thl1, rt2, thl2, & ! In
                       crt1, cthl1, crt2, cthl2, & ! In
                       cloud_frac1, cloud_frac2, mu1(iiLH_s_mellor), mu2(iiLH_s_mellor), & ! In
                       X_nl_one_lev(1:n_micro_calls,iiLH_s_mellor), & ! In
                       X_nl_one_lev(1:n_micro_calls,iiLH_t_mellor), & ! In
                       X_mixt_comp_one_lev, & ! In
                       LH_rt, LH_thl ) ! Out

! Compute some diagnostics
!       print*, 'C=', mixt_frac*cloud_frac1 + (1-mixt_frac)*cloud_frac2
!       print*, 'rtm_anl=', mixt_frac*rt1+(1-mixt_frac)*rt2, 'rtm_est=', mean(rt(1:n),n)
!       print*, 'thl_anl=',mixt_frac*thl1+(1-mixt_frac)*thl2, 'thlm_est=',mean(thl(1:n),n)
!       print*, 'rtpthlp_coef_est=', corrcoef(rt,thl,n)

    ! Convert lognormal variates (e.g. Nc and rr) to lognormal
    forall ( sample = 1:n_micro_calls )
      where ( l_d_variable_lognormal )
        X_nl_one_lev(sample,:) = exp( X_nl_one_lev(sample,:) )
      end where
    end forall

! Test diagnostics
!       print*, 'mean X_nl_one_lev(:,2)=', mean(X_nl_one_lev(1:n,2),n)
!       print*, 'mean X_nl_one_lev(:,3)=', mean(X_nl_one_lev(1:n,3),n)
!       print*, 'mean X_nl_one_lev(:,4)=', mean(X_nl_one_lev(1:n,4),n)
!       print*, 'mean X_nl_one_lev(:,5)=', mean(X_nl_one_lev(1:n,5),n)

!       print*, 'std X_nl_one_lev(:,2)=', std(X_nl_one_lev(1:n,2),n)
!       print*, 'std X_nl_one_lev(:,3)=', std(X_nl_one_lev(1:n,3),n)
!       print*, 'std X_nl_one_lev(:,4)=', std(X_nl_one_lev(1:n,4),n)
!       print*, 'std X_nl_one_lev(:,5)=', std(X_nl_one_lev(1:n,5),n)

!       print*, 'corrcoef X_nl_one_lev(:,2:3)=',
!     .            corrcoef( X_nl_one_lev(1:n,2), X_nl_one_lev(1:n,3), n )


    return
  end subroutine sample_points

!-------------------------------------------------------------------------------

!-----------------------------------------------------------------------
  subroutine rtpthlp_2_sptp( d_variables, Sigma_rtthlw, crt, cthl, Sigma_stw )

! Description:
!   Transform covariance matrix from rt', theta_l' coordinates
!   to s', t' coordinates.
!   Use linear approximation for s', t'.
! References:
!   None
!-----------------------------------------------------------------------

    use constants, only: fstderr ! Constant(s)

    implicit none

    ! External

    intrinsic :: matmul

    ! Input Variables

    integer, intent(in) :: d_variables ! Number of variates

    double precision, intent(in), dimension(d_variables,d_variables) :: &
      Sigma_rtthlw ! D x D dimensional covariance matrix of rt, thl, w, ...

    double precision, intent(in) :: &
      crt, cthl ! Coefficients that define s', t'

    ! Output Variables

    double precision, intent(out), dimension(d_variables,d_variables) :: &
      Sigma_stw ! Covariance matrix in terms of s', t' ordering of Sigma_stw is s, t, w, ...

    ! Local Variables

    integer j, k
    double precision, dimension(d_variables,d_variables) :: &
      T_transpose, T, M_int

    ! ---- Begin Code ----

    ! Check that matrix is large enough (at least 2x2)
    if ( d_variables < 3 ) then
      write(fstderr,*) 'Error: Input matrix too small in rtpthlp_2_stpthlp.'
      stop
    end if

    ! Transform rt-thl-w matrix to s-t-w
    !    according to Sigma_stw = T*Sigma_rtthlw*T_transpose.

    ! Set up transformation matrix, T_transpose
    T_transpose(1,1) = crt
    T_transpose(2,1) = -1.d0*cthl
    T_transpose(1,2) = crt
    T_transpose(2,2) = 1.d0*cthl

    ! Put zeros in the two off-diagonal blocks
    if (d_variables > 2) then
      do j=1,2
        do k=3,d_variables
          T_transpose(j,k) = 0.d0
          T_transpose(k,j) = 0.d0
        end do
      end do
    end if

    ! Put identity matrix in the lower diagonal block.
    if (d_variables > 2) then
      do j = 3, d_variables
        do k = 3, d_variables
          if (j == k) then
            T_transpose(j,k) = 1.d0
          else
            T_transpose(j,k) = 0.d0
          end if
        end do
      end do
    end if

    ! Multiply to obtain M_int = Sigma_rtthlw * T_transpose
    M_int = matmul( Sigma_rtthlw, T_transpose )

    ! Set up other transformation matrix, T
    T      = T_transpose
    T(1,2) = -1.d0*cthl
    T(2,1) = crt

    ! Perform final matrix multiplication: Sigma_stw = T * M_int
    Sigma_stw = matmul( T, M_int )

    return
  end subroutine rtpthlp_2_sptp
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
  subroutine generate_uniform_sample( n_micro_calls, nt_repeat, dp1, p_matrix, X_u_one_lev )

! Description:
!   Generates a matrix X that contains a Latin Hypercube sample.
!   The sample is uniformly distributed.
! References:
!   See Art B. Owen (2003), ``Quasi-Monte Carlo Sampling,"
!      a chapter from SIGGRAPH 2003
!-------------------------------------------------------------------------------

    use mt95, only: genrand_real ! Constants

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      n_micro_calls, & ! `n'   Number of calls to microphysics (normally=2)
      nt_repeat,     & ! `n_t' Num. random samples before sequence repeats (normally=10)
      dp1              !  d+1  Number of variates plus 1 (normally=6)

    integer, intent(in), dimension(n_micro_calls,dp1) :: &
      p_matrix    !n x dp1 array of permuted integers

    ! Output Variables

    real(kind=genrand_real), intent(out), dimension(n_micro_calls,dp1) :: &
      X_u_one_lev ! n by dp1 matrix, X, each row of which is a dp1-dimensional sample

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
    do j = 1,n_micro_calls
      do k = 1,dp1
        X_u_one_lev(j,k) = choose_permuted_random( nt_repeat, p_matrix(j,k) )
      end do
    end do

!        print*, 'p_matrix(:,1)= ', p_matrix(:,1)
!        print*, 'p_matrix(:,dp1)= ', p_matrix(:,dp1)
!        print*, 'X(:,1)= ', X(:,:)

    return
  end subroutine generate_uniform_sample

!----------------------------------------------------------------------
  function choose_permuted_random( nt_repeat, p_matrix_element )

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

    choose_permuted_random = (1.0_genrand_real/nt_repeat)*(p_matrix_element + rand )

    return
  end function choose_permuted_random

!----------------------------------------------------------------------
  subroutine gaus_mixt_points( n_micro_calls, d_variables, mixt_frac, mu1, mu2, Sigma1, Sigma2, &
                               cloud_frac1, cloud_frac2, X_u_one_lev, s_pts, &
                               X_mixt_comp_one_lev, &
                               X_gm )
! Description:
!   Generates n random samples from a d-dimensional Gaussian-mixture PDF.
!   Uses Latin hypercube method.
! References:
!   None
!----------------------------------------------------------------------

    use constants, only:  &
      fstderr  ! Constant(s)

    use error_code, only:  &
      clubb_at_least_debug_level  ! Procedure(s)

    implicit none

    ! Input Variables

    integer, intent(in) :: &
      n_micro_calls, &  ! Number of calls to microphysics (normally=2) 
      d_variables       ! Number of variates (normally=5)

    double precision, intent(in) :: &
      mixt_frac,     & ! Mixture fraction of Gaussians
      cloud_frac1, cloud_frac2   ! Cloud fraction associated w/ 1st, 2nd mixture component

    real, intent(in), dimension(d_variables) :: &
      mu1, mu2 ! d-dimensional column vector of means of 1st, 2nd Gaussians

    double precision, intent(in), dimension(d_variables,d_variables) :: &
      Sigma1, Sigma2 ! dxd dimensional covariance matrices

    ! Latin hypercube sample from uniform distribution from a particular grid level
    double precision, intent(in), dimension(n_micro_calls,d_variables+1) :: &
      X_u_one_lev 

    double precision, intent(in), dimension(n_micro_calls) :: &
      s_pts ! n-dimensional vector giving values of s

    integer, intent(in), dimension(n_micro_calls) :: &
      X_mixt_comp_one_lev ! Which mixture component we're in

    ! Output Variables

    double precision, intent(out), dimension(n_micro_calls,d_variables) :: &
      X_gm ! [n by d] matrix, X_gm, each row of which is a d-dimensional sample

    ! Local Variables

    integer :: j, sample
!   double precision, dimension(n_micro_calls) :: std_normal
    double precision, dimension(d_variables) :: std_normal
!   double precision :: fraction_1

    ! ---- Begin Code ----

    ! Handle some possible errors re: proper ranges of mixt_frac, 
    ! cloud_frac1, cloud_frac2.
    if (mixt_frac > 1.0d0 .or. mixt_frac < 0.0d0) then
      write(fstderr,*) 'Error in gaus_mixt_points:  ',  &
                       'mixture fraction, mixt_frac, does not lie in [0,1].'
      stop
    end if
    if (cloud_frac1 > 1.0d0 .or. cloud_frac1 < 0.0d0) then
      write(fstderr,*) 'Error in gaus_mixt_points:  ',  &
                       'cloud fraction 1, cloud_frac1, does not lie in [0,1].'
      stop
    end if
    if (cloud_frac2 > 1.0d0 .or. cloud_frac2 < 0.0d0) then
      write(fstderr,*) 'Error in gaus_mixt_points:  ',  &
                       'cloud fraction 2, cloud_frac2, does not lie in [0,1].'
      stop
    end if

    ! Make sure there is some cloud.
    if (mixt_frac*cloud_frac1 < 0.001d0 .and. (1-mixt_frac)*cloud_frac2 < 0.001d0) then
      if ( clubb_at_least_debug_level( 1 ) ) then
        write(fstderr,*) 'Error in gaus_mixt_points:  ',  &
                         'there is no cloud or almost no cloud!'
      end if
    end if

    do sample = 1, n_micro_calls

      ! From Latin hypercube sample, generate standard normal sample
      do j = 1, d_variables
        std_normal(j) = ltqnorm( X_u_one_lev(sample,j) )
      end do

      ! Determine which mixture fraction we are in.
      if ( X_mixt_comp_one_lev(sample) == 1 ) then
        call gaus_condt( d_variables, &
                         std_normal, mu1, Sigma1, s_pts(sample), &  ! In
                         X_gm(sample, 1:d_variables) ) ! Out

      else if ( X_mixt_comp_one_lev(sample) == 2 ) then
        call gaus_condt( d_variables, &
                         std_normal, mu2, Sigma2, s_pts(sample), &  ! In
                         X_gm(sample, 1:d_variables) ) ! Out

      else
        stop "Error determining mixture component in gaus_mixt_points"

      end if

      ! Loop to get new sample
    end do

    return
  end subroutine gaus_mixt_points

!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
  subroutine truncate_gaus_mixt( n_micro_calls, d_variables, col, mixt_frac, mu1, mu2, & 
                  Sigma1, Sigma2, cloud_frac1, cloud_frac2, X_u_one_lev, &
                  X_mixt_comp_one_lev, truncated_column )
! Description:
!   Converts sample points drawn from a uniform distribution
!    to truncated Gaussian points.
! References:
!   None
!-------------------------------------------------------------------------------

    use constants, only:  &
      fstderr  ! Constant(s)

    use error_code, only:  &
      clubb_at_least_debug_level  ! Procedure(s)

    implicit none

    ! Input Variables

    integer, intent(in) :: &
      n_micro_calls, &  ! Number of calls to microphysics (normally=2) 
      d_variables,   &  ! Number of variates (normally=5)
      col               ! Scalar indicated which column of X_nl_one_lev to truncate

    double precision, intent(in) :: &
      mixt_frac,    & ! Mixture fraction of Gaussians
      cloud_frac1, cloud_frac2  ! Cloud fraction associated w/ 1st, 2nd mixture component

    real, intent(in), dimension(d_variables) :: &
      mu1, mu2 ! d-dimensional column vector of means of 1st, 2nd Gaussians

    double precision, intent(in), dimension(d_variables,d_variables) :: &
      Sigma1, Sigma2 ! dxd dimensional covariance matrices

    ! Latin hypercube sample from uniform distribution from a particular grid level
    double precision, intent(in), dimension(n_micro_calls,d_variables+1) :: &
      X_u_one_lev 

    integer, intent(in), dimension(n_micro_calls) :: &
      X_mixt_comp_one_lev ! Whether we're in the first or second mixture component

    ! Output Variables

    ! A column vector of length n that is transformed from a Gaussian PDF to truncated Gaussian PDF.
    double precision, intent(out), dimension(n_micro_calls) :: truncated_column

    ! Local Variables

    integer :: sample
    double precision :: s_std
!   double precision :: fraction_1

    ! ---- Begin Code ----

    ! Handle some possible errors re: proper ranges of mixt_frac, 
    ! cloud_frac1, cloud_frac2.
    if ( (mixt_frac > 1.0d0) .or. (mixt_frac < 0.0d0) ) then
      write(fstderr,*) 'Error in truncate_gaus_mixt:  ',  &
                       'mixture fraction, mixt_frac, does not lie in [0,1].'
      stop
    end if
    if ( (cloud_frac1 > 1.0d0) .or. (cloud_frac1 < 0.0d0) ) then
      write(fstderr,*) 'Error in truncate_gaus_mixt:  ',  &
                       'cloud fraction 1, cloud_frac1, does not lie in [0,1].'
      stop
    end if
    if ( (cloud_frac2 > 1.0d0) .or. (cloud_frac2 < 0.0d0) ) then
      write(fstderr,*) 'Error in truncate_gaus_mixt:  ',  &
                       'cloud fraction 2, cloud_frac2, does not lie in [0,1].'
      stop
    end if

    ! Make sure there is some cloud.
    if (mixt_frac*cloud_frac1 < 0.001d0 .and. (1-mixt_frac) * cloud_frac2 < 0.001d0) then
      if ( clubb_at_least_debug_level( 1 ) ) then
        write(fstderr,*) 'Error in truncate_gaus_mixt:  ',  &
                         'there is no cloud or almost no cloud!'
      end if
    end if

    ! Make s PDF (1st column) a truncated Gaussian.
    ! This allows us to sample solely from the cloud points.
    do sample = 1, n_micro_calls

      ! Choose which mixture fraction we are in.
      ! Account for cloud fraction.
      ! Follow M. E. Johnson (1987), p. 56.
!     fraction_1 = mixt_frac * cloud_frac1 / &
!                  ( mixt_frac * cloud_frac1 + (1.d0 - mixt_frac) *cloud_frac2 )
!     if ( X_u_one_lev( sample, d_variables+1 ) < fraction_1 ) then
      if  ( X_mixt_comp_one_lev(sample) == 1 ) then
        ! Replace first dimension (s) with
        !  sample from cloud (i.e. truncated standard Gaussian)
        s_std = ltqnorm( X_u_one_lev( sample, col ) * cloud_frac1 + (1.d0 - &
          cloud_frac1) )
        ! Convert to nonstandard normal with mean mu1 and variance Sigma1
        truncated_column(sample) =  & 
                   s_std * sqrt( Sigma1(col,col) ) + mu1(col)
      else if ( X_mixt_comp_one_lev(sample) == 2 ) then

        ! Replace first dimension (s) with
        !   sample from cloud (i.e. truncated Gaussian)
        s_std = ltqnorm( X_u_one_lev( sample, col ) * cloud_frac2 + (1.d0 - &
          cloud_frac2) )

        ! Convert to nonstandard normal with mean mu2 and variance Sigma2
        truncated_column(sample) =  & 
                      s_std * sqrt( Sigma2(col,col) ) + mu2(col)
      else
        stop "Error in truncate_gaus_mixt"
      end if

      ! Loop to get new sample
    end do

    return
  end subroutine truncate_gaus_mixt

!-----------------------------------------------------------------------
  double precision function ltqnorm( p )
! Description:
!   This function is ported to Fortran from the same function written in Matlab, see the following
!   description of this function.  Hongli Jiang, 2/17/2004
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
!   result has a relative error less than 1.15e-9.  A last refinement by
!   Halley's rational method is applied to achieve full machine precision.

!   Author:      Peter J. Acklam
!   Time-stamp:  2003-04-23 08:26:51 +0200
!   E-mail:      pjacklam@online.no
!   URL:         http://home.online.no/~pjacklam
!-----------------------------------------------------------------------

    use constants, only: Pi_DP ! Variable(s)

    implicit none

    ! External

    intrinsic :: log, sqrt

    ! Input Variable(s)

    double precision, intent(in) :: p


    ! Local Variable(s)

    double precision a1, a2, a3, a4, a5, a6, b1, b2, b3, b4, b5, & 
                     c1, c2, c3, c4, c5, c6, d1, d2, d3, d4

    double precision q, r, z, z1, plow, phigh

!       double preciseion e, erf_dp, u

!       Occurs in constants.F now.  Isn't actually used currently.
!        double precision, parameter :: pi=3.1415926d0

! Coefficients in rational approximations.
! equivalent: a(1)=a1, a(2)=a2, and etc, when a(1) is in Matlab.
! Similarly for b, c, and d's
    parameter (a1 = -3.969683028665376d+01,  & 
               a2 = 2.209460984245205d+02, & 
               a3 = -2.759285104469687d+02,  & 
               a4 = 1.383577518672690d+02, & 
               a5 = -3.066479806614716d+01,  & 
               a6 = 2.506628277459239d+00)
    parameter (b1 = -5.447609879822406d+01,  & 
               b2 = 1.615858368580409d+02, & 
               b3 = -1.556989798598866d+02,  & 
               b4 = 6.680131188771972d+01, & 
               b5 = -1.328068155288572d+01)
    parameter (c1 = -7.784894002430293d-03,  & 
               c2 = -3.223964580411365d-01, & 
               c3 = -2.400758277161838d+00,  & 
               c4 = -2.549732539343734d+00, & 
               c5 =  4.374664141464968d+00,  & 
               c6 =  2.938163982698783d+00)
    parameter (d1 =  7.784695709041462d-03,  & 
               d2 =  3.224671290700398d-01, & 
               d3 =  2.445134137142996d+00,  & 
               d4 =  3.754408661907416d+00)

    ! Default initialization
    z = 0.0

!  Define break-points.
    plow  = 0.02425d0
    phigh = 1.d0 - plow

!  Initialize output array. Don't need this in Fortran
!   z = zeros(size(p));

!  Rational approximation for lower region:
    if (p > 0.d0 .and. p < plow) then
      q = sqrt( -2 * log(p) )
      z = (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/ & 
                ((((d1*q+d2)*q+d3)*q+d4)*q+1.d0)
!  Rational approximation for central region:
    else if (p >= plow .and. p <= phigh) then
      q = p - 0.5d0
      r = q * q
      z = (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q & 
                 /(((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1.d0)
! Rational approximation for upper region:
    else if (p > phigh .and. p < 1.d0) then
      q  = sqrt( -2.d0 * log(1.d0 - p) )
      z  = -(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) & 
                  /((((d1*q+d2)*q+d3)*q+d4)*q+1.d0)
    end if

!  Case when P = 0: z = -inf, to create inf z =-1./0.,
!     to create NaN's inf*inf.
    z1 = 0.d0
    if (p == 0.d0) then
      z = (-1.d0)/z1
    end if

! Case when P = 1:, z=inf
    if(p == 1.d0)then
      z = 1.d0/z1
    end if

!  Cases when output will be NaN:
!   k = p < 0 | p > 1 | isnan(p);
! usually inf*inf --> NaN's.
    if (p < 0.d0 .or. p > 1d0) then
      z = (1.d0/z1)**2
    end if

!  The relative error of the approximation has absolute value less
!  than 1.15e-9.  One iteration of Halley's rational method (third
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
!         e = 0.5*(1.0 - erf_dp(-z/sqrt(2.))) - p          ! error
!         u = e * sqrt(2*pi_dp) * exp(z**2/2)       ! f(z)/df(z)
!         z = z - u/( 1 + z*u/2 )               ! Halley's method
!       end if

! return z as double precision:
    ltqnorm = z

    return
  end function ltqnorm

!----------------------------------------------------------------------
  subroutine gaus_condt( d_variables, std_normal, mu, Sigma, s_pt, & 
                         nonstd_normal )
! Description:
!   Using Gaussian conditional distributions given s,
!   convert a standard, uncorrelated Gaussian to one with
!   mean mu and covariance structure Sigma.

! References:
!   Follow M. E. Johnson, ``Multivariate Statistical Simulation," p50.
!----------------------------------------------------------------------

    use matrix_operations, only: linear_symm_upper_eqn_solve ! Procedure(s)

!   use matrix_operations, only: linear_eqn_solve ! Procedure(s)

    use error_code, only: clubb_at_least_debug_level

    implicit none

    ! External
    intrinsic :: sqrt, matmul, reshape

    ! Input Variables

    integer, intent(in) :: d_variables ! Number of variates (normally=5)

    double precision, intent(in), dimension(d_variables) :: &
      std_normal ! nxd matrix of n independent samples from d-variate standard normal distribution

    real, intent(in), dimension(d_variables,1) :: &
      mu ! d-dimensional column vector of means of Gaussian

    double precision, intent(in), dimension(d_variables,d_variables) :: &
      Sigma ! dxd dimensional covariance matrix

    double precision, intent(in) :: s_pt ! Value of Mellor's s

    ! Output Variables

    ! nxd matrix of n samples from d-variate normal distribution
    !   with mean mu and covariance structure Sigma
    double precision, intent(out) :: &
      nonstd_normal(d_variables)

    ! Local Variables

    integer :: v, j ! Loop iterators

    double precision :: mu_two, mu_condt_const

    double precision, dimension(d_variables,1) :: mu_one

    double precision, dimension(d_variables,d_variables) :: &
      Sigma_oneone

    double precision, dimension(d_variables,1) :: &
      Sigma_onetwo

    double precision, dimension(1,d_variables) :: &
      Sigma_twoone

    double precision :: Sigma_twotwo, Sigma_condt

    ! Local intermediate quantities
    double precision, dimension(1,d_variables) :: Sigma_int
    double precision, dimension(1,1) :: dum

!   do i = 1, d_variables
!     do j = 1, d_variables
!       write(6,'(g12.4)',advance='no') Sigma(i,j)
!     end do
!     write(6,*) ""
!   end do
!   write(6,*) ""
!   pause

    ! First store value of s
    nonstd_normal(1) = s_pt
    ! Loop over variables t, w, . . .
    ! [Conditional distribution of
    !      X_two given X_one] = x_one = std_normal(1:n,v)
    ! is N [ mu_two + Sigma_twoone*Inverse(Sigma_oneone)*(x_one-mu_one),
    !    Sigma_twotwo - Sigma_twoone*Inverse(Sigma_oneone)*Sigma_onetwo]
    ! Here 'one' refers to given variables, 'two' refers to values to find.
    ! Loop over all variables other than s_mellor:
    do v=2, d_variables

      ! Means
      mu_one(1:(v-1),1) = mu(1:(v-1),1)
      mu_two            = mu(v,1)

      ! Matrices used to compute variances
      Sigma_twoone(1,1:(v-1)) = Sigma(v,1:(v-1))
      Sigma_onetwo(1:(v-1),1) = Sigma(1:(v-1),v)
      Sigma_oneone( 1:(v-1) , 1:(v-1) )  & 
           = Sigma( 1:(v-1) , 1:(v-1) )
      Sigma_twotwo = Sigma(v,v)

      ! Check whether input matrix is symmetric.
      if ( clubb_at_least_debug_level( 2 ) ) then
        do j=1,(v-1)
          if ( abs( Sigma_twoone(1,j) - Sigma_onetwo(j,1) ) > epsilon( dum ) ) then
            write(0,*) 'Error: matrix Sigma not symmetric in gaus_condt.'
            write(0,*) 'Sigma in gaus_condt=', Sigma
          end if
        end do
      end if

      ! Compute an intermediate matrix, Sigma_int(1,1:(v-1)),
      !    that is needed several times below.
      ! Solve A * X = B for X, where X here is sigma_int.
      call linear_symm_upper_eqn_solve &
           ( n = v-1, a = Sigma_oneone(1:v-1,1:v-1), b = Sigma_twoone(1,1:v-1), &
             x = Sigma_int(1,1:v-1) )

      ! Constant scalar needed to compute mean of non-standard normal.
      ! mu_condt_const = mu_two - Sigma_twoone*Sigma_oneone_inv*mu_one

      dum = matmul( Sigma_int(:,1:v-1), mu_one(1:v-1,:) )

      mu_condt_const = mu_two - dum(1,1)

      ! Variance of non-standard normal for variable v.
! Sigma_condt is a scalar.
! Sigma_condt = Sigma_twotwo -
!                    Sigma_twoone*Sigma_oneone_inv*Sigma_onetwo
      dum = matmul( Sigma_int(:,1:v-1), Sigma_onetwo(1:v-1,:) )
      Sigma_condt = Sigma_twotwo - dum(1,1)

      ! Compute next element of nonstd_normal from prior elements
! nonstd_normal(v) = sqrt(Sigma_condt)*std_normal(v)
!      + mu_condt_const ...
!      + Sigma_twoone*Sigma_oneone_inv*nonstd_normal(1:(v-1))
      dum = matmul( Sigma_int(:,1:v-1), &
                    reshape( nonstd_normal(1:v-1), (/v-1, 1/) ) )
      nonstd_normal(v) = sqrt( Sigma_condt )*std_normal(v)  & 
         + mu_condt_const  &
         + dum(1,1)

      ! Loop to obtain new variable
    end do ! 2..d_variables

    return
  end subroutine gaus_condt

!-----------------------------------------------------------------------
  subroutine st_2_rtthl( n_micro_calls, mixt_frac, rt1, thl1, rt2, thl2, & 
                         crt1, cthl1, crt2, cthl2, & 
                         cloud_frac1, cloud_frac2, mu_s1, mu_s2, &
                         s_mellor, t_mellor, X_mixt_comp_one_lev, &
                         LH_rt, LH_thl )
! Description:
!   Converts from s, t variables to rt, thl
! References:
!   None
!-----------------------------------------------------------------------

    use constants, only:  &
        fstderr  ! Constant(s)

    use error_code, only:  &
        clubb_at_least_debug_level  ! Procedure(s)

    implicit none

    ! Input

    integer, intent(in) :: &
      n_micro_calls   ! Number of calls to microphysics (normally=2)

    double precision, intent(in) :: &
      mixt_frac,   & ! Mixture fraction of Gaussians 'mixt_frac' [-]
      rt1, rt2,    & ! n dimensional column vector of rt         [kg/kg]
      thl1, thl2,  & ! n dimensional column vector of thetal     [K]
      crt1, crt2,  & ! Constants from plumes 1 & 2 of rt
      cthl1, cthl2   ! Constants from plumes 1 & 2 of thetal

    double precision, intent(in) :: &
      cloud_frac1, cloud_frac2 ! Cloud fraction associated with 1st / 2nd mixture component

    real, intent(in) :: &
      mu_s1, mu_s2

    ! n-dimensional column vector of Mellor's s and t, including mean and perturbation
    double precision, intent(in), dimension(n_micro_calls) :: &
      s_mellor, & 
      t_mellor 

    integer, dimension(n_micro_calls), intent(in) :: &
      X_mixt_comp_one_lev ! Whether we're in the first or second mixture component

    ! Output variables

    real, dimension(n_micro_calls), intent(out) :: &
      LH_rt, LH_thl ! n-dimensional column vectors of rt and thl, including mean and perturbation

    ! Local

    integer :: sample

    ! ---- Begin Code ----

    ! Handle some possible errors re: proper ranges of mixt_frac,
    ! cloud_frac1, cloud_frac2.

    if ( mixt_frac > 1.0d0 .or. mixt_frac < 0.0d0 ) then
      write(fstderr,*) 'Error in st_2_rtthl:  ',  &
                       'mixture fraction, mixt_frac, does not lie in [0,1].'
      stop
    end if
    if ( cloud_frac1 > 1.0d0 .or. cloud_frac1 < 0.0d0 ) then
      write(fstderr,*) 'Error in st_2_rtthl:  ',  &
                       'cloud fraction 1, cloud_frac1, does not lie in [0,1].'
      stop
    end if
    if ( cloud_frac2 > 1.0d0 .or. cloud_frac2 < 0.0d0 ) then
      write(fstderr,*) 'Error in st_2_rtthl:  ',  &
                       'cloud fraction 2, cloud_frac2, does not lie in [0,1].'
      stop
    end if

    ! Make sure there is some cloud.
    if ( mixt_frac*cloud_frac1 < 0.001d0 .and. (1-mixt_frac)*cloud_frac2 < 0.001d0 ) then
      if ( clubb_at_least_debug_level( 1 ) ) then
        write(fstderr,*) 'Error in st_2_rtthl:  ',  &
                         'there is no cloud or almost no cloud!'
      end if
    end if

    do sample = 1, n_micro_calls

      ! Choose which mixture fraction we are in.
      ! Account for cloud fraction.
      ! Follow M. E. Johnson (1987), p. 56.
!     fraction_1     = mixt_frac*cloud_frac1 / &
!                      (mixt_frac*cloud_frac1+(1-mixt_frac)*cloud_frac2)

!     if ( in_mixt_frac_1( X_u_one_lev(sample,d_variables+1), fraction_1 ) ) then
      if ( X_mixt_comp_one_lev(sample) == 1 ) then
        LH_rt(sample)  = rt1 + (0.5d0/crt1)*(s_mellor(sample)-mu_s1) +  & 
                           (0.5d0/crt1)*t_mellor(sample)
        LH_thl(sample) = thl1 + (-0.5d0/cthl1)*(s_mellor(sample)-mu_s1) +  & 
                           (0.5d0/cthl1)*t_mellor(sample)

      else if ( X_mixt_comp_one_lev(sample) == 2 ) then
      ! mixture fraction 2
        LH_rt(sample)  = rt2 + (0.5d0/crt2)*(s_mellor(sample)-mu_s2) +  & 
                           (0.5d0/crt2)*t_mellor(sample)
        LH_thl(sample) = thl2 + (-0.5d0/cthl2)*(s_mellor(sample)-mu_s2) +  & 
                           (0.5d0/cthl2)*t_mellor(sample)

      else 
        stop "Error determining mixture fraction in st_2_rtthl"

      end if

      ! Loop to get new sample
    end do

    return
  end subroutine st_2_rtthl
!-------------------------------------------------------------------------------
  subroutine log_sqd_normalized( Xm, Xp2_on_Xm2, X_tol, &
                                 X1, X2, sX1, sX2 )
! Description:
!
! References:
!   None
!-------------------------------------------------------------------------------
    implicit none

    ! External
    intrinsic :: log

    ! Input Variables
    double precision, intent(in) :: &
      Xm,           & ! Mean X          [units vary]
      Xp2_on_Xm2,   & ! X'^2 / X^2      [-]
      X_tol           ! tolerance on X  [units vary]

    ! Output Variables
    double precision, intent(out) :: &
      X1, X2,  & ! PDF parameters for mean of plume 1, 2   [units vary]
      sX1, sX2   ! PDF parameters for width of plume 1, 2  [units vary]

    ! ---- Begin Code ----

    if ( Xm >= X_tol ) then
      sX1 = log( 1. + Xp2_on_Xm2 )
      sX2 = sX1
      X1 = 0.5 * log( Xm**2 / ( 1. + Xp2_on_Xm2 ) )
      X2 = X1
    else
      sX1 = log( 1. + epsilon( Xp2_on_Xm2 ) )
      sX2 = sX1
      X1 = 0.5 * log( ( X_tol * 0.999 )**2 / ( 1. + epsilon( Xp2_on_Xm2 ) ) )
      X2 = X1
    end if

    return
  end subroutine log_sqd_normalized
  
end module generate_lh_sample_module

