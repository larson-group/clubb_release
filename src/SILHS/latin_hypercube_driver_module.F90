!-------------------------------------------------------------------------------
! $Id$
!===============================================================================
module latin_hypercube_driver_module

  implicit none

  ! Constant Parameters
  logical, parameter, private :: &
    l_diagnostic_iter_check      = .true., &  ! Check for a problem in iteration
    l_output_2D_lognormal_dist   = .false., & ! Output a 2D netCDF file of the lognormal variates
    l_output_2D_uniform_dist     = .false.    ! Output a 2D netCDF file of the uniform distribution

  integer, private :: &
    prior_iter ! Prior iteration number (for diagnostic purposes)
!$omp threadprivate( prior_iter )

  integer, public :: k_lh_start ! For an assertion check
!$omp threadprivate( k_lh_start )


  private ! Default scope

#ifdef SILHS
  public :: lh_subcolumn_generator, lh_microphys_driver, latin_hypercube_2D_output, &
    latin_hypercube_2D_close, stats_accumulate_lh, lh_subcolumn_generator_mod

  private :: stats_accumulate_uniform_lh

  contains

!-------------------------------------------------------------------------------
  subroutine lh_subcolumn_generator &
             ( iter, d_variables, num_samples, sequence_length, nz, &
               thlm, pdf_params, wm_zt, delta_zm, rcm, Ncnm, rvm, &
               hydromet, sigma2_on_mu2_ip_array_cloud, sigma2_on_mu2_ip_array_below, &
               corr_array_cloud, corr_array_below, Lscale_vert_avg, &
               X_nl_all_levs, X_mixt_comp_all_levs, lh_rt, lh_thl, &
               lh_sample_point_weights )

! Description:
!   Generate sample points of moisture, temperature, et cetera for the purpose
!   of computing tendencies with a microphysics or radiation scheme.
!
! References:
!   None
!-------------------------------------------------------------------------------

    use corr_matrix_module, only: &
      iiPDF_chi    ! Variables

    use latin_hypercube_arrays, only: &
      height_time_matrix ! Variables

    use parameters_model, only: hydromet_dim ! Variable

    use permute_height_time_module, only: & 
      permute_height_time ! Procedure

    use generate_lh_sample_module, only: & 
      generate_lh_sample, & ! Procedure
      generate_uniform_sample

    use output_2D_samples_module, only: &
      output_2D_lognormal_dist_file, & ! Procedure(s)
      output_2D_uniform_dist_file

    use pdf_parameter_module, only: &
      pdf_parameter  ! Type

    use constants_clubb, only: &
      fstderr, & ! Constant
      zero_threshold

    use parameters_microphys, only: &
      l_lh_vert_overlap, &  ! Variables
      l_lh_cloud_weighted_sampling

    use error_code, only: &
      clubb_at_least_debug_level ! Procedure

    use mt95, only: genrand_real ! Constants

    use mt95, only: genrand_init ! Procedure

    use clubb_precision, only: & 
      dp, & ! double precision
      core_rknd, &
      stat_rknd

    implicit none

    ! External
    intrinsic :: allocated, mod, maxloc, epsilon

    ! Parameter Constants
    real( kind = core_rknd ), parameter :: &
      cloud_frac_thresh = 0.001_core_rknd ! Threshold for sampling preferentially within cloud

    ! Find in and out of cloud points using the rejection method rather than scaling
    logical, parameter :: &
      l_use_rejection_method = .false.

    ! Input Variables
    integer, intent(in) :: &
      iter,            & ! Model iteration number
      d_variables,     & ! Number of variables to sample
      num_samples,     & ! Number of samples per variable
      sequence_length, & ! nt_repeat/num_samples; number of timesteps before sequence repeats.
      nz                 ! Number of vertical model levels

    type(pdf_parameter), dimension(nz), intent(in) :: & 
      pdf_params ! PDF parameters       [units vary]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      thlm,      & ! Liquid potential temperature       [K]
      wm_zt,     & ! Mean w                             [m/s]
      delta_zm     ! Difference in moment. altitudes    [m]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      rcm,  & ! Liquid water mixing ratio                              [kg/kg]
      Ncnm, & ! Cloud nuclei concentration (simplified); Nc=Ncn*H(s)   [#/kg]
      rvm     ! Vapor water mixing ratio                               [kg/kg]

    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(in) :: &
      hydromet ! Hydrometeor species    [units vary]

    real( kind = core_rknd ), dimension(d_variables), intent(in) :: &
      sigma2_on_mu2_ip_array_cloud, & ! sigma_x^2/mu_x^2 ip; cloudy levs. [-]
      sigma2_on_mu2_ip_array_below    ! sigma_x^2/mu_x^2 ip; clear levs.  [-]

    real( kind = core_rknd ), dimension(d_variables,d_variables), intent(in) :: &
      corr_array_cloud, & ! Correlation for hydrometeor species [-]
      corr_array_below

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      Lscale_vert_avg ! 3pt vertical average of Lscale  [m]

    ! Output Variables
    real( kind = dp ), intent(out), dimension(nz,num_samples,d_variables) :: &
      X_nl_all_levs ! Sample that is transformed ultimately to normal-lognormal

    integer, intent(out), dimension(nz,num_samples) :: &
      X_mixt_comp_all_levs ! Which mixture component we're in

    real( kind = core_rknd ), intent(out), dimension(nz,num_samples) :: &
      lh_rt, lh_thl ! Sample of total water and liquid potential temperature [kg/kg],[K]

    real( kind = core_rknd ), intent(out), dimension(num_samples) :: &
      lh_sample_point_weights

    ! Local variables

    real( kind = dp ), dimension(nz,num_samples,(d_variables+1)) :: &
      X_u_all_levs ! Sample drawn from uniform distribution

    integer :: p_matrix(num_samples,d_variables+1)

    real(kind=core_rknd) :: lh_start_cloud_frac ! Cloud fraction at k_lh_start [-]

    ! Try to obtain 12 digit accuracy for a diagnostic mean
    real(kind=dp) :: mean_weight

    real(kind=dp), dimension(num_samples) :: &
      X_u_dp1_k_lh_start, X_u_chi_k_lh_start

    real(kind=dp), dimension(nz) :: &
      X_vert_corr ! Vertical correlation of a variate   [-]

    real(kind=dp) :: X_u_temp

    ! Number of random samples before sequence of repeats (normally=10)
    integer :: nt_repeat

    integer :: &
      i_rmd, &   ! Remainder of ( iter-1 / sequence_length )
      k, sample  ! Loop iterators

    integer :: ivar ! Loop iterator

    real( kind = core_rknd ), parameter :: &
      cloud_frac_max_weighted_smpl = 0.5_core_rknd ! Use cloud weighted sampling only if cloud
                                                       ! fraction is less than
                                                       ! cloud_frac_max_weighted_smpl

    logical :: l_cloudy_sample ! Whether a sample point is cloudy or clear air

    logical :: &
      l_small_nonzero_cloud_frac ! True if cloud fraction is greater than cloud_frac_thresh
                                 ! and less than cloud_frac_max_weighted_smpl

    integer, dimension(1) :: tmp_loc

    ! ---- Begin Code ----

    ! Get rid of compiler warning
    l_small_nonzero_cloud_frac = .false.

    nt_repeat = num_samples * sequence_length

    if ( .not. allocated( height_time_matrix ) ) then
      ! If this is first time latin_hypercube_driver is called, then allocate
      ! the height_time_matrix and set the prior iteration number for debugging
      ! purposes.
      allocate( height_time_matrix(nz, nt_repeat, d_variables+1) )

      prior_iter = iter

      ! Check for a bug where the iteration number isn't incrementing correctly,
      ! which will lead to improper sampling.
    else if ( l_diagnostic_iter_check ) then

      if ( prior_iter /= iter-1 ) then
        write(fstderr,*) "The iteration number in latin_hypercube_driver is"// &
        " not incrementing properly."

      else
        prior_iter = iter

      end if

    end if ! First call to the driver

    ! Sanity checks for l_lh_cloud_weighted_sampling
    if ( l_lh_cloud_weighted_sampling .and. mod( num_samples, 2 ) /= 0 ) then
      write(fstderr,*) "Cloud weighted sampling requires num_samples to be divisible by 2."
      stop "Fatal error."
    end if
    if ( l_lh_cloud_weighted_sampling .and. sequence_length /= 1 ) then
      write(fstderr,*) "Cloud weighted sampling requires sequence length be equal to 1."
      stop "Fatal error."
    end if

    ! Initialize the sample point weights to 1.0
    lh_sample_point_weights(1:num_samples)  = 1.0_core_rknd

    ! Latin hypercube sample generation
    ! Generate height_time_matrix, an nz x nt_repeat x d_variables array of random integers
    i_rmd = mod( iter-1, sequence_length )

    if ( i_rmd == 0 ) then
      call permute_height_time( nz, nt_repeat, d_variables+1, & ! intent(in)
                                height_time_matrix )            ! intent(out)
    end if
    ! End Latin hypercube sample generation

    !--------------------------------------------------------------
    ! Latin hypercube sampling
    !--------------------------------------------------------------


    ! For a 100 level fixed grid, this looks to be about the middle of the cloud for RICO
!   k_lh_start = 50
    tmp_loc    = maxloc( rcm )
    k_lh_start = tmp_loc(1) ! Attempt using the maximal value of rcm for now

    ! If there's no cloud k_lh_start appears to end up being 1. Check if
    ! k_lh_start is 1 or nz and set it to the middle of the domain in that
    ! case.
    if ( k_lh_start == nz .or. k_lh_start == 1 ) then
      k_lh_start = nz / 2
    end if

    if ( l_lh_cloud_weighted_sampling ) then

      ! Determine cloud fraction at k_lh_start
      lh_start_cloud_frac = &
      pdf_params(k_lh_start)%mixt_frac * pdf_params(k_lh_start)%cloud_frac1 &
        + (1.0_core_rknd-pdf_params(k_lh_start)%mixt_frac) * pdf_params(k_lh_start)%cloud_frac2

      ! Determine p_matrix at k_lh_start
      p_matrix(1:num_samples,1:(d_variables+1)) = &
        height_time_matrix(k_lh_start, &
        (num_samples*i_rmd+1):(num_samples*i_rmd+num_samples), &
        1:(d_variables+1))
    else
      lh_start_cloud_frac = -9999.0_core_rknd  
         ! This assignment eliminates a g95 compiler warning for an uninitialized variable. 
         !-meyern

    end if ! l_lh_cloud_weighted_sampling

    if ( l_lh_vert_overlap ) then

      ! Choose which rows of LH sample to feed into closure at the k_lh_start level
      p_matrix(1:num_samples,1:(d_variables+1)) = &
        height_time_matrix(k_lh_start, &
        (num_samples*i_rmd+1):(num_samples*i_rmd+num_samples), &
        1:(d_variables+1))

      ! Generate the uniform distribution using the Mersenne twister at the k_lh_start level
      !  X_u has one extra dimension for the mixture component.
      call generate_uniform_sample( num_samples, nt_repeat, d_variables+1, p_matrix, & ! In
                                    X_u_all_levs(k_lh_start,:,:) ) ! Out

      do sample = 1, num_samples

        ! Get the uniform sample of d+1 (telling us if we're in component 1 or
        ! component 2), and the uniform sample for chi
        X_u_dp1_k_lh_start(sample)      = X_u_all_levs(k_lh_start,sample,d_variables+1)
        X_u_chi_k_lh_start(sample) = X_u_all_levs(k_lh_start,sample,iiPDF_chi)

        if ( l_lh_cloud_weighted_sampling ) then

          ! Save the cloud fraction as a weight for averaging preferentially
          ! within cloud
          l_small_nonzero_cloud_frac = (lh_start_cloud_frac < cloud_frac_max_weighted_smpl &
                                  .and. lh_start_cloud_frac > cloud_frac_thresh)
          if ( .not. l_small_nonzero_cloud_frac ) then
            lh_sample_point_weights(sample)  = 1.0_core_rknd
            ! There's no cloud or cloud fraction is >= cloud_frac_max_weighted_smpl, 
            ! so we do nothing

          else ! If we're in a partly cloud gridbox then we continue to the code below

            ! Detect which half of the sample points are in clear air and which half are in
            ! the cloudy air
            if ( p_matrix(sample,iiPDF_chi) < ( num_samples / 2 ) ) then

              l_cloudy_sample = .false.
              lh_sample_point_weights(sample) = 2._core_rknd* &
                ( 1.0_core_rknd - lh_start_cloud_frac )
            else
              l_cloudy_sample = .true.
              lh_sample_point_weights(sample) = 2._core_rknd* lh_start_cloud_frac
            end if

            if ( l_use_rejection_method ) then
              ! Use the rejection method to select points that are in or out of cloud
              call choose_X_u_reject &
                   ( l_cloudy_sample, pdf_params(k_lh_start)%cloud_frac1, & ! In
                     pdf_params(k_lh_start)%cloud_frac2, pdf_params(k_lh_start)%mixt_frac, & !In
                     cloud_frac_thresh, & ! In
                     X_u_dp1_k_lh_start(sample), X_u_chi_k_lh_start(sample) ) ! In/out

            else ! Transpose and scale the points to be in or out of cloud
              call choose_X_u_scaled &
                   ( l_cloudy_sample, & ! In
                     p_matrix(sample,iiPDF_chi), num_samples, & ! In 
                     pdf_params(k_lh_start)%cloud_frac1, pdf_params(k_lh_start)%cloud_frac2, & ! In
                     pdf_params(k_lh_start)%mixt_frac, & !In
                     X_u_dp1_k_lh_start(sample), X_u_chi_k_lh_start(sample) ) ! In/out

            end if

          end if ! Cloud fraction is between cloud_frac_thresh and 50%

        end if ! l_lh_cloud_weighted_sampling

      end do ! 1..num_samples

      ! Use a fixed number for the vertical correlation.
!     X_vert_corr(1:nz) = 0.95_dp

      ! Compute vertical correlation using a formula based on Lscale, the
      ! the difference in height levels, and an empirical constant
      X_vert_corr(1:nz) = &
        real(compute_vert_corr( nz, delta_zm, Lscale_vert_avg ), kind = dp)

      ! Assertion check for the vertical correlation
      if ( clubb_at_least_debug_level( 1 ) ) then
        if ( any( X_vert_corr > 1.0_dp ) .or. any( X_vert_corr < 0.0_dp ) ) then
          write(fstderr,*) "The vertical correlation in latin_hypercube_driver"// &
            "is not in the correct range"
          do k = 1, nz
            write(fstderr,*) "k = ", k,  "Vert. correlation = ", X_vert_corr(k)
          end do
        end if ! Some correlation isn't between [0,1]
      end if ! clubb_at_least_debug_level 1

      do sample = 1, num_samples
        ! Correlate chi vertically
        call compute_arb_overlap &
             ( nz, k_lh_start, &  ! In
               X_u_chi_k_lh_start(sample), X_vert_corr, & ! In
               X_u_all_levs(:,sample,iiPDF_chi) ) ! Out
        ! Correlate the d+1 variate vertically (used to compute the mixture
        ! component later)
        call compute_arb_overlap &
             ( nz, k_lh_start, &  ! In
               X_u_dp1_k_lh_start(sample), X_vert_corr, & ! In
               X_u_all_levs(:,sample,d_variables+1) ) ! Out

        ! Use these lines to make all variates vertically correlated, using the
        ! same correlation we used above for chi and the d+1 variate
        do ivar = 1, d_variables
          if ( ivar /= iiPDF_chi ) then
            X_u_temp = X_u_all_levs(k_lh_start,sample,ivar)
            call compute_arb_overlap &
                 ( nz, k_lh_start, &  ! In
                   X_u_temp, X_vert_corr, & ! In
                   X_u_all_levs(:,sample,ivar) ) ! Out
          end if
        end do ! 1..d_variables
      end do ! 1..num_samples
      ! %% Debug %%
      ! Testing what happens when we clip uniformally distributed variates to
      ! avoid extreme values.
!     where ( X_u_all_levs (:,:,2:d_variables) > 0.99_core_rknd ) &
!       X_u_all_levs(:,:,2:d_variables) = 0.99_core_rknd
      ! %% End Debug %%

    else ! Random overlap

      do k = 1, nz
        ! Choose which rows of LH sample to feed into closure.
        p_matrix(1:num_samples,1:(d_variables+1)) = &
          height_time_matrix(k, num_samples*i_rmd+1:num_samples*i_rmd+num_samples, &
                             1:d_variables+1)

        ! Generate the uniform distribution using the Mersenne twister
        !  X_u has one extra dimension for the mixture component.
        call generate_uniform_sample( num_samples, nt_repeat, d_variables+1, p_matrix, & ! In
                                      X_u_all_levs(k,:,:) ) ! Out
      end do ! 1..nz

    end if ! l_lh_vert_overlap

    ! Determine mixture component for all levels
    do k = 1, nz

      where ( in_mixt_comp_1( X_u_all_levs(k,:,d_variables+1), &
           real(pdf_params(k)%mixt_frac, kind = dp) ) )
        X_mixt_comp_all_levs(k,:) = 1
      else where
        X_mixt_comp_all_levs(k,:) = 2
      end where

    end do ! k = 1 .. nz

    ! Assertion check for whether half of sample points are cloudy.
    ! This is for the uniform sample only.  Another assertion check is in the
    ! estimate_kessler_microphys_module for X_nl_all_levs.
    if ( l_lh_cloud_weighted_sampling ) then
      if ( clubb_at_least_debug_level( 2 ) .and. l_small_nonzero_cloud_frac ) then
        call assert_check_half_cloudy &
             ( num_samples, pdf_params(k_lh_start)%cloud_frac1, &
               pdf_params(k_lh_start)%cloud_frac2, X_mixt_comp_all_levs(k_lh_start,:), &
               X_u_all_levs(k_lh_start,:,iiPDF_chi) )

      end if ! Maximal overlap, debug_level 2, and cloud-weighted averaging
    end if ! l_lh_cloud_weighted_sampling

    ! Assertion check to ensure that the sample point weights sum to approximately 1
    if ( l_lh_cloud_weighted_sampling .and. clubb_at_least_debug_level( 2 ) ) then
      mean_weight = 0._dp
      do sample = 1, num_samples
        mean_weight = mean_weight + real( lh_sample_point_weights(sample), kind=dp )
      end do
      mean_weight = mean_weight / real( num_samples, kind=dp )

      ! Using more precision for mean_weight should make this work out.
      ! The formula below could probably be redefined to estimate maximal ulps
      ! given the precision of lh_sample_point_weights and the number of
      ! num_samples, but the formula below seems to be an ok approximation
      ! when we're using 4 or 8 byte precision floats.
      ! -dschanen 19 Nov 2010
      if ( abs( mean_weight - 1.0_dp ) > &
          real(real( num_samples, kind=core_rknd ) &
          * epsilon( lh_sample_point_weights ), kind = dp) ) then
        write(fstderr,*) "Error in cloud weighted sampling code ", "mean_weight = ", mean_weight
        stop
      end if

    end if ! l_lh_cloud_weighted_sampling .and. clubb_at_least_debug_level( 2 )

    ! Upwards loop
    do k = k_lh_start, nz, 1
      ! Generate LH sample, represented by X_u and X_nl, for level k
      do sample = 1, num_samples, 1
        call generate_lh_sample &
             ( d_variables, hydromet_dim, &  ! In
               wm_zt(k), rcm(k), Ncnm(k), rvm(k), thlm(k), & ! In
               pdf_params(k)%mixt_frac, & ! In
               pdf_params(k)%w1, pdf_params(k)%w2, pdf_params(k)%rc1, pdf_params(k)%rc2, & ! In
               pdf_params(k)%varnce_w1, pdf_params(k)%varnce_w2, & ! In
               pdf_params(k)%thl1, pdf_params(k)%thl2, & ! In
               pdf_params(k)%varnce_thl1, pdf_params(k)%varnce_thl2, & ! In
               pdf_params(k)%rt1, pdf_params(k)%rt2, & ! In
               pdf_params(k)%varnce_rt1, pdf_params(k)%varnce_rt2, & ! In
               pdf_params(k)%chi_1, pdf_params(k)%chi_2, & ! In
               pdf_params(k)%stdev_chi_1, pdf_params(k)%stdev_chi_2, & ! In
               pdf_params(k)%stdev_eta_1, pdf_params(k)%stdev_eta_2, & ! In
               pdf_params(k)%covar_chi_eta_1, pdf_params(k)%covar_chi_eta_2, & ! In
               pdf_params(k)%crt1, pdf_params(k)%crt2, & ! In
               pdf_params(k)%cthl1, pdf_params(k)%cthl2, & ! In
               hydromet(k,:), sigma2_on_mu2_ip_array_cloud, sigma2_on_mu2_ip_array_below, & ! In
               corr_array_cloud, corr_array_below, & ! In
               X_u_all_levs(k,sample,:), X_mixt_comp_all_levs(k,sample), & ! In
               lh_rt(k,sample), lh_thl(k,sample), X_nl_all_levs(k,sample,:) ) ! Out
      end do ! sample = 1, num_samples, 1
    end do ! k = k_lh_start..nz

      ! Downwards loop
    do k = k_lh_start-1, 1, -1
      do sample = 1, num_samples, 1
        call generate_lh_sample &
             ( d_variables, hydromet_dim, &  ! In
               wm_zt(k), rcm(k), Ncnm(k), rvm(k), thlm(k), & ! In
               pdf_params(k)%mixt_frac, & ! In
               pdf_params(k)%w1, pdf_params(k)%w2, pdf_params(k)%rc1, pdf_params(k)%rc2, & ! In
               pdf_params(k)%varnce_w1, pdf_params(k)%varnce_w2, & ! In
               pdf_params(k)%thl1, pdf_params(k)%thl2, & ! In
               pdf_params(k)%varnce_thl1, pdf_params(k)%varnce_thl2, & ! In
               pdf_params(k)%rt1, pdf_params(k)%rt2, & ! In
               pdf_params(k)%varnce_rt1, pdf_params(k)%varnce_rt2, & ! In
               pdf_params(k)%chi_1, pdf_params(k)%chi_2, & ! In
               pdf_params(k)%stdev_chi_1, pdf_params(k)%stdev_chi_2, & ! In
               pdf_params(k)%stdev_eta_1, pdf_params(k)%stdev_eta_2, & ! In
               pdf_params(k)%covar_chi_eta_1, pdf_params(k)%covar_chi_eta_2, & ! In
               pdf_params(k)%crt1, pdf_params(k)%crt2, & ! In
               pdf_params(k)%cthl1, pdf_params(k)%cthl2, & ! In
               hydromet(k,:), sigma2_on_mu2_ip_array_cloud, sigma2_on_mu2_ip_array_below, & ! In
               corr_array_cloud, corr_array_below, & ! In
               X_u_all_levs(k,sample,:), X_mixt_comp_all_levs(k,sample), & ! In
               lh_rt(k,sample), lh_thl(k,sample), X_nl_all_levs(k,sample,:) ) ! Out
      end do ! sample = 1, num_samples, 1
    end do ! k_lh_start-1..1

    if ( l_output_2D_lognormal_dist ) then
      call output_2D_lognormal_dist_file( nz, num_samples, d_variables, &
                                          real(X_nl_all_levs, kind = stat_rknd), &
                                          lh_rt, lh_thl )
    end if
    if ( l_output_2D_uniform_dist ) then
      call output_2D_uniform_dist_file( nz, num_samples, d_variables+1, &
                                        real(X_u_all_levs, kind = genrand_real), &
                                        X_mixt_comp_all_levs, p_matrix )
    end if

    ! Verify total water isn't negative
    if ( any( lh_rt < 0._core_rknd) ) then
      if ( clubb_at_least_debug_level( 1 ) ) then
        write(fstderr,*) "Total water negative in LH sample point"
        write(fstderr,*) "Applying non-conservative hard clipping to rv sample point."
      end if ! clubb_at_least_debug_level( 1 )
      where ( lh_rt < 0._core_rknd)
        lh_rt = zero_threshold
      end where
    end if ! Some rv_all_points(:,sample) < 0

    return
  end subroutine lh_subcolumn_generator


!-------------------------------------------------------------------------------
  subroutine lh_subcolumn_generator_mod &
             ( iter, d_variables, num_samples, sequence_length, nz, & ! In
               pdf_params, delta_zm, rcm, Lscale_vert_avg, & ! In
               mu1, mu2, sigma1, sigma2, & ! In
               corr_cholesky_mtx_1, corr_cholesky_mtx_2, & ! In
               hydromet_pdf_params, & ! In
               X_nl_all_levs, X_mixt_comp_all_levs, lh_rt, lh_thl, & ! Out
               lh_sample_point_weights ) ! Out

! Description:
!   Generate sample points of moisture, temperature, et cetera for the purpose
!   of computing tendencies with a microphysics or radiation scheme.
!
! References:
!   None
!-------------------------------------------------------------------------------

    use corr_matrix_module, only: &
      iiPDF_chi    ! Variables

    use latin_hypercube_arrays, only: &
      height_time_matrix ! Variables

    use permute_height_time_module, only: &
      permute_height_time ! Procedure

    use generate_lh_sample_module, only: &
      generate_lh_sample_mod, & ! Procedure
      generate_uniform_sample

    use output_2D_samples_module, only: &
      output_2D_lognormal_dist_file, & ! Procedure(s)
      output_2D_uniform_dist_file

    use pdf_parameter_module, only: &
      pdf_parameter  ! Type

    use hydromet_pdf_parameter_module, only: &
      hydromet_pdf_parameter ! Type

    use constants_clubb, only: &
      fstderr, & ! Constant
      zero_threshold

    use parameters_microphys, only: &
      l_lh_vert_overlap, &  ! Variables
      l_lh_cloud_weighted_sampling

    use error_code, only: &
      clubb_at_least_debug_level ! Procedure

    use mt95, only: genrand_real ! Constants

    use mt95, only: genrand_init ! Procedure

    use clubb_precision, only: &
      dp, & ! double precision
      core_rknd, &
      stat_rknd

    implicit none

    ! External
    intrinsic :: allocated, mod, maxloc, epsilon, transpose

    ! Parameter Constants
    real( kind = core_rknd ), parameter :: &
      cloud_frac_thresh = 0.001_core_rknd ! Threshold for sampling preferentially within cloud

    ! Find in and out of cloud points using the rejection method rather than scaling
    logical, parameter :: &
      l_use_rejection_method = .false.

    integer, parameter :: &
      d_uniform_extra = 2   ! Number of variables that are included in the uniform sample but not in
                            ! the lognormal sample. Currently:
                            !
                            ! d_variables+1: Mixture component, for choosing PDF component
                            ! d_variables+2: Precipitation fraction, for determining precipitation

    ! Input Variables
    integer, intent(in) :: &
      iter,            & ! Model iteration number
      d_variables,     & ! Number of variables to sample
      num_samples,     & ! Number of samples per variable
      sequence_length, & ! nt_repeat/num_samples; number of timesteps before sequence repeats.
      nz                 ! Number of vertical model levels

    type(pdf_parameter), dimension(nz), intent(in) :: &
      pdf_params ! PDF parameters       [units vary]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      delta_zm, &  ! Difference in moment. altitudes    [m]
      rcm          ! Liquid water mixing ratio          [kg/kg]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      Lscale_vert_avg ! 3pt vertical average of Lscale  [m]

    ! Output Variables
    real( kind = dp ), intent(out), dimension(nz,num_samples,d_variables) :: &
      X_nl_all_levs ! Sample that is transformed ultimately to normal-lognormal

    integer, intent(out), dimension(nz,num_samples) :: &
      X_mixt_comp_all_levs ! Which mixture component we're in

    real( kind = core_rknd ), intent(out), dimension(nz,num_samples) :: &
      lh_rt, lh_thl ! Sample of total water and liquid potential temperature [kg/kg],[K]

    real( kind = core_rknd ), intent(out), dimension(num_samples) :: &
      lh_sample_point_weights

    ! More Input Variables!
    real( kind = dp ), dimension(d_variables,d_variables,nz), intent(in) :: &
      corr_cholesky_mtx_1, & ! Correlations Cholesky matrix (1st comp.)  [-]
      corr_cholesky_mtx_2    ! Correlations Cholesky matrix (2nd comp.)  [-]

    real( kind = core_rknd ), dimension(d_variables,nz), intent(in) :: &
      mu1,    & ! Means of the hydrometeors, 1st comp. (s, t, w, <hydrometeors>)  [units vary]
      mu2,    & ! Means of the hydrometeors, 2nd comp. (s, t, w, <hydrometeors>)  [units vary]
      sigma1, & ! Stdevs of the hydrometeors, 1st comp. (s, t, w, <hydrometeors>) [units vary]
      sigma2    ! Stdevs of the hydrometeors, 2nd comp. (s, t, w, <hydrometeors>) [units vary]

    type(hydromet_pdf_parameter), dimension(nz), intent(in) :: &
      hydromet_pdf_params

    ! Local variables

    real( kind = dp ), dimension(nz,num_samples,(d_variables+d_uniform_extra)) :: &
      X_u_all_levs ! Sample drawn from uniform distribution

    integer :: p_matrix(num_samples,d_variables+d_uniform_extra)

    real(kind=core_rknd) :: lh_start_cloud_frac ! Cloud fraction at k_lh_start [-]

    ! Try to obtain 12 digit accuracy for a diagnostic mean
    real(kind=dp) :: mean_weight

    real(kind=dp), dimension(num_samples) :: &
      X_u_dp1_k_lh_start, X_u_chi_k_lh_start

    real(kind=dp), dimension(nz) :: &
      X_vert_corr ! Vertical correlation of a variate   [-]

    real(kind=dp) :: X_u_temp

    ! Number of random samples before sequence of repeats (normally=10)
    integer :: nt_repeat

    integer :: &
      i_rmd, &   ! Remainder of ( iter-1 / sequence_length )
      k, sample  ! Loop iterators

    integer :: ivar ! Loop iterator

    real( kind = core_rknd ), parameter :: &
      cloud_frac_max_weighted_smpl = 0.5_core_rknd ! Use cloud weighted sampling only if cloud
                                                       ! fraction is less than
                                                       ! cloud_frac_max_weighted_smpl

    logical :: l_cloudy_sample ! Whether a sample point is cloudy or clear air

    logical :: &
      l_small_nonzero_cloud_frac ! True if cloud fraction is greater than cloud_frac_thresh
                                 ! and less than cloud_frac_max_weighted_smpl

    integer, dimension(1) :: tmp_loc

    logical, dimension(nz,num_samples) :: &
      l_in_precip   ! Whether sample is in precipitation

    ! Precipitation fraction in a component of the PDF, for each sample
    real( kind = dp ), dimension(num_samples) :: precip_frac_i

    ! ---- Begin Code ----

    ! Get rid of compiler warning
    l_small_nonzero_cloud_frac = .false.

    nt_repeat = num_samples * sequence_length

    if ( .not. allocated( height_time_matrix ) ) then
      ! If this is first time latin_hypercube_driver is called, then allocate
      ! the height_time_matrix and set the prior iteration number for debugging
      ! purposes.
      allocate( height_time_matrix(nz, nt_repeat, d_variables+d_uniform_extra) )

      prior_iter = iter

      ! Check for a bug where the iteration number isn't incrementing correctly,
      ! which will lead to improper sampling.
    else if ( l_diagnostic_iter_check ) then

      if ( prior_iter /= iter-1 ) then
        write(fstderr,*) "The iteration number in latin_hypercube_driver is"// &
        " not incrementing properly."

      else
        prior_iter = iter

      end if

    end if ! First call to the driver

    ! Sanity checks for l_lh_cloud_weighted_sampling
    if ( l_lh_cloud_weighted_sampling .and. mod( num_samples, 2 ) /= 0 ) then
      write(fstderr,*) "Cloud weighted sampling requires num_samples to be divisible by 2."
      stop "Fatal error."
    end if
    if ( l_lh_cloud_weighted_sampling .and. sequence_length /= 1 ) then
      write(fstderr,*) "Cloud weighted sampling requires sequence length be equal to 1."
      stop "Fatal error."
    end if

    ! Initialize the sample point weights to 1.0
    lh_sample_point_weights(1:num_samples)  = 1.0_core_rknd

    ! Latin hypercube sample generation
    ! Generate height_time_matrix, an nz x nt_repeat x d_variables array of random integers
    i_rmd = mod( iter-1, sequence_length )

    if ( i_rmd == 0 ) then
      call permute_height_time( nz, nt_repeat, d_variables+d_uniform_extra, & ! intent(in)
                                height_time_matrix )                          ! intent(out)
    end if
    ! End Latin hypercube sample generation

    !--------------------------------------------------------------
    ! Latin hypercube sampling
    !--------------------------------------------------------------


    ! For a 100 level fixed grid, this looks to be about the middle of the cloud for RICO
!   k_lh_start = 50
    tmp_loc    = maxloc( rcm )
    k_lh_start = tmp_loc(1) ! Attempt using the maximal value of rcm for now

    ! If there's no cloud k_lh_start appears to end up being 1. Check if
    ! k_lh_start is 1 or nz and set it to the middle of the domain in that
    ! case.
    if ( k_lh_start == nz .or. k_lh_start == 1 ) then
      k_lh_start = nz / 2
    end if

    if ( l_lh_cloud_weighted_sampling ) then

      ! Determine cloud fraction at k_lh_start
      lh_start_cloud_frac = &
      pdf_params(k_lh_start)%mixt_frac * pdf_params(k_lh_start)%cloud_frac1 &
        + (1.0_core_rknd-pdf_params(k_lh_start)%mixt_frac) * pdf_params(k_lh_start)%cloud_frac2

      ! Determine p_matrix at k_lh_start
      p_matrix(1:num_samples,1:(d_variables+d_uniform_extra)) = &
        height_time_matrix(k_lh_start, &
        (num_samples*i_rmd+1):(num_samples*i_rmd+num_samples), &
        1:(d_variables+d_uniform_extra))
    else
      lh_start_cloud_frac = -9999.0_core_rknd
         ! This assignment eliminates a g95 compiler warning for an uninitialized variable.
         !-meyern

    end if ! l_lh_cloud_weighted_sampling

    if ( l_lh_vert_overlap ) then

      ! Choose which rows of LH sample to feed into closure at the k_lh_start level
      p_matrix(1:num_samples,1:(d_variables+d_uniform_extra)) = &
        height_time_matrix(k_lh_start, &
        (num_samples*i_rmd+1):(num_samples*i_rmd+num_samples), &
        1:(d_variables+d_uniform_extra))

      ! Generate the uniform distribution using the Mersenne twister at the k_lh_start level
      call generate_uniform_sample( num_samples, nt_repeat, d_variables+d_uniform_extra, & ! In
                                    p_matrix, & ! In
                                    X_u_all_levs(k_lh_start,:,:) ) ! Out

      do sample = 1, num_samples

        ! Get the uniform sample of d+1 (telling us if we're in component 1 or
        ! component 2), and the uniform sample for chi
        X_u_dp1_k_lh_start(sample)      = X_u_all_levs(k_lh_start,sample,d_variables+1)
        X_u_chi_k_lh_start(sample) = X_u_all_levs(k_lh_start,sample,iiPDF_chi)

        if ( l_lh_cloud_weighted_sampling ) then

          ! Save the cloud fraction as a weight for averaging preferentially
          ! within cloud
          l_small_nonzero_cloud_frac = (lh_start_cloud_frac < cloud_frac_max_weighted_smpl &
                                  .and. lh_start_cloud_frac > cloud_frac_thresh)
          if ( .not. l_small_nonzero_cloud_frac ) then
            lh_sample_point_weights(sample)  = 1.0_core_rknd
            ! There's no cloud or cloud fraction is >= cloud_frac_max_weighted_smpl,
            ! so we do nothing

          else ! If we're in a partly cloud gridbox then we continue to the code below

            ! Detect which half of the sample points are in clear air and which half are in
            ! the cloudy air
            if ( p_matrix(sample,iiPDF_chi) < ( num_samples / 2 ) ) then

              l_cloudy_sample = .false.
              lh_sample_point_weights(sample) = 2._core_rknd* &
                ( 1.0_core_rknd - lh_start_cloud_frac )
            else
              l_cloudy_sample = .true.
              lh_sample_point_weights(sample) = 2._core_rknd* lh_start_cloud_frac
            end if

            if ( l_use_rejection_method ) then
              ! Use the rejection method to select points that are in or out of cloud
              call choose_X_u_reject &
                   ( l_cloudy_sample, pdf_params(k_lh_start)%cloud_frac1, & ! In
                     pdf_params(k_lh_start)%cloud_frac2, pdf_params(k_lh_start)%mixt_frac, & !In
                     cloud_frac_thresh, & ! In
                     X_u_dp1_k_lh_start(sample), X_u_chi_k_lh_start(sample) ) ! Out

            else ! Transpose and scale the points to be in or out of cloud
              call choose_X_u_scaled &
                   ( l_cloudy_sample, & ! In
                     p_matrix(sample,iiPDF_chi), num_samples, & ! In
                     pdf_params(k_lh_start)%cloud_frac1, pdf_params(k_lh_start)%cloud_frac2, & ! In
                     pdf_params(k_lh_start)%mixt_frac, & !In
                     X_u_dp1_k_lh_start(sample), X_u_chi_k_lh_start(sample) ) ! Out

            end if

          end if ! Cloud fraction is between cloud_frac_thresh and 50%

        end if ! l_lh_cloud_weighted_sampling

      end do ! 1..num_samples

      ! Use a fixed number for the vertical correlation.
!     X_vert_corr(1:nz) = 0.95_dp

      ! Compute vertical correlation using a formula based on Lscale, the
      ! the difference in height levels, and an empirical constant
      X_vert_corr(1:nz) = &
        real(compute_vert_corr( nz, delta_zm, Lscale_vert_avg ), kind = dp)

      ! Assertion check for the vertical correlation
      if ( clubb_at_least_debug_level( 1 ) ) then
        if ( any( X_vert_corr > 1.0_dp ) .or. any( X_vert_corr < 0.0_dp ) ) then
          write(fstderr,*) "The vertical correlation in latin_hypercube_driver"// &
            "is not in the correct range"
          do k = 1, nz
            write(fstderr,*) "k = ", k,  "Vert. correlation = ", X_vert_corr(k)
          end do
        end if ! Some correlation isn't between [0,1]
      end if ! clubb_at_least_debug_level 1

      do sample = 1, num_samples
        ! Correlate chi vertically
        call compute_arb_overlap &
             ( nz, k_lh_start, &  ! In
               X_u_chi_k_lh_start(sample), X_vert_corr, & ! In
               X_u_all_levs(:,sample,iiPDF_chi) ) ! Out
        ! Correlate the d+1 variate vertically (used to compute the mixture
        ! component later)
        call compute_arb_overlap &
             ( nz, k_lh_start, &  ! In
               X_u_dp1_k_lh_start(sample), X_vert_corr, & ! In
               X_u_all_levs(:,sample,d_variables+1) ) ! Out

        ! Correlate the d+2 variate vertically (used to determine precipitation
        ! later)
        X_u_temp = X_u_all_levs(k_lh_start,sample,d_variables+2)
        call compute_arb_overlap &
             ( nz, k_lh_start, &  ! In
               X_u_temp, X_vert_corr, & ! In
               X_u_all_levs(:,sample,d_variables+2) )

        ! Use these lines to make all variates vertically correlated, using the
        ! same correlation we used above for chi and the d+1 variate
        do ivar = 1, d_variables
          if ( ivar /= iiPDF_chi ) then
            X_u_temp = X_u_all_levs(k_lh_start,sample,ivar)
            call compute_arb_overlap &
                 ( nz, k_lh_start, &  ! In
                   X_u_temp, X_vert_corr, & ! In
                   X_u_all_levs(:,sample,ivar) ) ! Out
          end if
        end do ! 1..d_variables
      end do ! 1..num_samples
      ! %% Debug %%
      ! Testing what happens when we clip uniformally distributed variates to
      ! avoid extreme values.
!     where ( X_u_all_levs (:,:,2:d_variables) > 0.99_core_rknd ) &
!       X_u_all_levs(:,:,2:d_variables) = 0.99_core_rknd
      ! %% End Debug %%

    else ! Random overlap

      do k = 1, nz
        ! Choose which rows of LH sample to feed into closure.
        p_matrix(1:num_samples,1:(d_variables+d_uniform_extra)) = &
          height_time_matrix(k, num_samples*i_rmd+1:num_samples*i_rmd+num_samples, &
                             1:d_variables+d_uniform_extra)

        ! Generate the uniform distribution using the Mersenne twister
        !  X_u has one extra dimension for the mixture component.
        call generate_uniform_sample( num_samples, nt_repeat, d_variables+d_uniform_extra, &
                                      p_matrix, & ! In
                                      X_u_all_levs(k,:,:) ) ! Out
      end do ! 1..nz

    end if ! l_lh_vert_overlap

    do k = 1, nz
      ! Determine mixture component for all levels
      where ( in_mixt_comp_1( X_u_all_levs(k,:,d_variables+1), &
           real(pdf_params(k)%mixt_frac, kind = dp) ) )
        X_mixt_comp_all_levs(k,:) = 1
      else where
        X_mixt_comp_all_levs(k,:) = 2
      end where

      ! Determine precipitation fraction
      where ( X_mixt_comp_all_levs(k,:) == 1 )
        precip_frac_i(:) = real( hydromet_pdf_params(k)%precip_frac_1, kind=dp )
      else where
        precip_frac_i(:) = real( hydromet_pdf_params(k)%precip_frac_2, kind=dp )
      end where

      ! Determine precipitation for all levels
      where ( in_precipitation( X_u_all_levs(k,:,d_variables+2), &
                  precip_frac_i(:) ) )
        l_in_precip(k,:) = .true.
      else where
        l_in_precip(k,:) = .false.
      end where

    end do ! k = 1 .. nz

    ! Assertion check for whether half of sample points are cloudy.
    ! This is for the uniform sample only.  Another assertion check is in the
    ! estimate_kessler_microphys_module for X_nl_all_levs.
    if ( l_lh_cloud_weighted_sampling ) then
      if ( clubb_at_least_debug_level( 2 ) .and. l_small_nonzero_cloud_frac ) then
        call assert_check_half_cloudy &
             ( num_samples, pdf_params(k_lh_start)%cloud_frac1, &
               pdf_params(k_lh_start)%cloud_frac2, X_mixt_comp_all_levs(k_lh_start,:), &
               X_u_all_levs(k_lh_start,:,iiPDF_chi) )

      end if ! Maximal overlap, debug_level 2, and cloud-weighted averaging
    end if ! l_lh_cloud_weighted_sampling

    ! Assertion check to ensure that the sample point weights sum to approximately 1
    if ( l_lh_cloud_weighted_sampling .and. clubb_at_least_debug_level( 2 ) ) then
      mean_weight = 0._dp
      do sample = 1, num_samples
        mean_weight = mean_weight + real( lh_sample_point_weights(sample), kind=dp )
      end do
      mean_weight = mean_weight / real( num_samples, kind=dp )

      ! Using more precision for mean_weight should make this work out.
      ! The formula below could probably be redefined to estimate maximal ulps
      ! given the precision of lh_sample_point_weights and the number of
      ! num_samples, but the formula below seems to be an ok approximation
      ! when we're using 4 or 8 byte precision floats.
      ! -dschanen 19 Nov 2010
      if ( abs( mean_weight - 1.0_dp ) > &
          real( num_samples, kind=dp ) * max( epsilon( mean_weight ), &
            real( epsilon( lh_sample_point_weights ), kind=dp ) ) ) then
        write(fstderr,*) "Error in cloud weighted sampling code ", "mean_weight = ", mean_weight
        stop
      end if

    end if ! l_lh_cloud_weighted_sampling .and. clubb_at_least_debug_level( 2 )

    call stats_accumulate_uniform_lh( nz, num_samples, l_in_precip, X_mixt_comp_all_levs, &
                                      lh_sample_point_weights)

    ! Sample loop
    do k = 1, nz
      ! Generate LH sample, represented by X_u and X_nl, for level k
      do sample = 1, num_samples, 1
        call generate_lh_sample_mod &
             ( d_variables, d_uniform_extra, & ! In
               pdf_params(k)%thl1, pdf_params(k)%thl2, & ! In
               pdf_params(k)%rt1, pdf_params(k)%rt2, & ! In
               pdf_params(k)%crt1, pdf_params(k)%crt2, & ! In
               pdf_params(k)%cthl1, pdf_params(k)%cthl2, & ! In
               mu1(:,k), mu2(:,k), sigma1(:,k), sigma2(:,k), & ! In
               corr_cholesky_mtx_1(:,:,k), & ! In
               corr_cholesky_mtx_2(:,:,k), & ! In
               X_u_all_levs(k,sample,:), X_mixt_comp_all_levs(k,sample), & ! In
               l_in_precip(k,sample), & ! In
               lh_rt(k,sample), lh_thl(k,sample), X_nl_all_levs(k,sample,:) ) ! Out
      end do ! sample = 1, num_samples, 1
    end do ! k = 1, nz

    if ( l_output_2D_lognormal_dist ) then
      call output_2D_lognormal_dist_file( nz, num_samples, d_variables, &
                                          real(X_nl_all_levs, kind = stat_rknd), &
                                          lh_rt, lh_thl )
    end if
    if ( l_output_2D_uniform_dist ) then
      call output_2D_uniform_dist_file( nz, num_samples, d_variables+1, &
                                        real(X_u_all_levs(:,:,1:d_variables+1), &
                                         kind = genrand_real), &
                                        X_mixt_comp_all_levs, p_matrix )
    end if

    ! Verify total water isn't negative
    if ( any( lh_rt < 0._core_rknd) ) then
      if ( clubb_at_least_debug_level( 1 ) ) then
        write(fstderr,*) "Total water negative in LH sample point"
        write(fstderr,*) "Applying non-conservative hard clipping to rv sample point."
      end if ! clubb_at_least_debug_level( 1 )
      where ( lh_rt < 0._core_rknd)
        lh_rt = zero_threshold
      end where
    end if ! Some rv_all_points(:,sample) < 0

    return
  end subroutine lh_subcolumn_generator_mod


  !=============================================================================
  subroutine lh_microphys_driver &
             ( dt, nz, num_samples, d_variables, &
               X_nl_all_levs, lh_rt, lh_thl, lh_sample_point_weights, &
               pdf_params, p_in_Pa, exner, rho, &
               rcm, w_std_dev, delta_zt, cloud_frac, &
               hydromet, X_mixt_comp_all_levs, Nc_in_cloud, &
               lh_hydromet_mc, lh_hydromet_vel, lh_Ncm_mc, &
               lh_rcm_mc, lh_rvm_mc, lh_thlm_mc, &
               lh_rtp2_mc, lh_thlp2_mc, lh_wprtp_mc, &
               lh_wpthlp_mc, lh_rtpthlp_mc, &
               microphys_sub )

    ! Description:
    !   Computes an estimate of the change due to microphysics given a set of
    !   subcolumns of thlm, rtm, et cetera from the subcolumn generator.
    !
    ! References:
    !   None
    !---------------------------------------------------------------------------

    use variables_diagnostic_module, only: & 
      lh_AKm,  & 
      AKm, & 
      AKstd, & 
      AKstd_cld, & 
      AKm_rcm, & 
      AKm_rcc, & 
      lh_rcm_avg

    use parameters_model, only: hydromet_dim ! Variable

    use pdf_parameter_module, only: &
      pdf_parameter  ! Type

    use est_kessler_microphys_module, only: &
      est_kessler_microphys

    use clubb_precision, only: &
      dp, & ! double precision
      core_rknd, &
      time_precision

    use error_code, only: &
      clubb_at_least_debug_level ! Procedure

    use estimate_scm_microphys_module, only: &
      est_single_column_tndcy

    implicit none

    ! Interface block
#include "microphys_interface.inc"

    ! Input Variables
    real( kind = time_precision ), intent(in) :: &
      dt ! Model timestep       [s]

    integer, intent(in) :: &
      d_variables,     & ! Number of variables to sample
      num_samples,   & ! Number of calls to microphysics per timestep (normally=2)
      nz               ! Number of vertical model levels

    ! Input Variables
    real( kind = dp ), intent(in), dimension(nz,num_samples,d_variables) :: &
      X_nl_all_levs ! Sample that is transformed ultimately to normal-lognormal

    integer, intent(in), dimension(nz,num_samples) :: &
      X_mixt_comp_all_levs ! Which mixture component we're in

    real( kind = core_rknd ), intent(in), dimension(nz,num_samples) :: &
      lh_rt, lh_thl ! Sample of total water and liquid potential temperature [kg/kg],[K]

    real( kind = core_rknd ), intent(in), dimension(num_samples) :: &
      lh_sample_point_weights ! Weight given the individual sample points

    type(pdf_parameter), dimension(nz), intent(in) :: & 
      pdf_params ! PDF parameters       [units vary]

    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(in) :: &
      hydromet ! Hydrometeor species    [units vary]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      cloud_frac,  & ! Cloud fraction               [-]
      w_std_dev,   & ! Standard deviation of w      [m/s]
      delta_zt,    & ! Change in meters with height [m]
      rcm,         & ! Liquid water mixing ratio    [kg/kg]
      p_in_Pa,     & ! Pressure                     [Pa]
      exner,       & ! Exner function               [-]
      rho            ! Density on thermo. grid      [kg/m^3]

    real( kind = core_rknd), dimension(nz), intent(in) :: &
      ! Constant value of N_c within cloud, to be used with l_const_Nc_in_cloud
      Nc_in_cloud

    ! Output Variables
    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(out) :: &
      lh_hydromet_mc, & ! LH estimate of hydrometeor time tendency          [(units vary)/s]
      lh_hydromet_vel   ! LH estimate of hydrometeor sedimentation velocity [m/s]

    ! Output Variables
    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      lh_Ncm_mc,     & ! LH estimate of time tndcy. of cloud droplet conc.     [num/kg/s]
      lh_rcm_mc,     & ! LH estimate of time tndcy. of liq. water mixing ratio [kg/kg/s]
      lh_rvm_mc,     & ! LH estimate of time tndcy. of vapor water mix. ratio  [kg/kg/s]
      lh_thlm_mc,    & ! LH estimate of time tndcy. of liquid potential temp.  [K/s]
      lh_rtp2_mc,    & ! LH microphysics tendency for <rt'^2>                  [(kg/kg)^2/s]
      lh_thlp2_mc,   & ! LH microphysics tendency for <thl'^2>                 [K^2/s]
      lh_wprtp_mc,   & ! LH microphysics tendency for <w'rt'>                  [m*(kg/kg)/s^2]
      lh_wpthlp_mc,  & ! LH microphysics tendency for <w'thl'>                 [m*K/s^2]
      lh_rtpthlp_mc    ! LH microphysics tendency for <rt'thl'>                [K*(kg/kg)/s]

    ! ---- Begin Code ----

    ! Perform LH and analytic microphysical calculations
    ! As a test of SILHS, compute an estimate of Kessler microphysics
    if ( clubb_at_least_debug_level( 2 ) ) then
       call est_kessler_microphys &
            ( nz, num_samples, d_variables, &                  ! Intent(in)
              X_nl_all_levs, pdf_params, rcm, cloud_frac, &      ! Intent(in)
              X_mixt_comp_all_levs, lh_sample_point_weights, &   ! Intent(in)
              lh_AKm, AKm, AKstd, AKstd_cld, &                   ! Intent(out)
              AKm_rcm, AKm_rcc, lh_rcm_avg )                     ! Intent(out)
    end if

    ! Call the latin hypercube microphysics driver for microphys_sub
    call est_single_column_tndcy &
         ( dt, nz, num_samples, d_variables, &                     ! Intent(in)
           k_lh_start, lh_rt, lh_thl, &                              ! Intent(in)
           X_nl_all_levs, lh_sample_point_weights, &                 ! Intent(in) 
           p_in_Pa, exner, rho, cloud_frac, w_std_dev, &             ! Intent(in)
           delta_zt, hydromet, rcm, Nc_in_cloud, &                   ! Intent(in)
           lh_hydromet_mc, lh_hydromet_vel, lh_Ncm_mc, &             ! Intent(out)
           lh_rvm_mc, lh_rcm_mc, lh_thlm_mc, &                       ! Intent(out)
           lh_rtp2_mc, lh_thlp2_mc, lh_wprtp_mc, &                   ! Intent(out)
           lh_wpthlp_mc, lh_rtpthlp_mc, &                            ! Intent(out)
           microphys_sub )                                           ! Intent(Procedure)

    return
  end subroutine lh_microphys_driver
!-------------------------------------------------------------------------------
  subroutine latin_hypercube_2D_output &
             ( fname_prefix, fdir, stats_tout, nz, &
               zt, time_initial )
!-------------------------------------------------------------------------------

    use corr_matrix_module, only: &
      iiPDF_chi, & ! Variables
      iiPDF_eta, &
      iiPDF_w, &
      iiPDF_rrain, & 
      iiPDF_rice, &
      iiPDF_rsnow, &
      iiPDF_rgraupel, &
      iiPDF_Nr, &
      iiPDF_Ni, &
      iiPDF_Nsnow, &
      iiPDF_Ngraupel, &
      iiPDF_Ncn

    use parameters_microphys, only: &
      lh_microphys_calls ! Variable

    use clubb_precision, only: &
      time_precision, & ! Constant
      core_rknd

    use output_2D_samples_module, only: &
      open_2D_samples_file ! Procedure

    use output_2D_samples_module, only: &
      lognormal_sample_file, & ! Instance of a type
      uniform_sample_file

    use corr_matrix_module, only: &
      d_variables ! Variable


    implicit none

    ! Input Variables
    character(len=*), intent(in) :: &
      fname_prefix, & ! Prefix for file name
      fdir            ! Directory for output

    real(kind=time_precision), intent(in) :: &
      stats_tout, & ! Frequency to write to disk        [s]
      time_initial  ! Initial time                      [s]

    integer, intent(in) :: &
      nz ! Number of vertical levels

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      zt ! Altitudes [m]

    ! Local Variables
    character(len=100), allocatable, dimension(:) :: &
      variable_names, variable_descriptions, variable_units

    integer :: i

    ! ---- Begin Code ----

    if ( l_output_2D_lognormal_dist ) then

      allocate( variable_names(d_variables+2), variable_descriptions(d_variables+2), &
                variable_units(d_variables+2) )

      variable_names(iiPDF_chi)        = "chi"
      variable_descriptions(iiPDF_chi) = "The variable 's' from Mellor 1977"
      variable_units(iiPDF_chi)        = "kg/kg"

      variable_names(iiPDF_eta)        = "eta"
      variable_descriptions(iiPDF_eta) = "The variable 't' from Mellor 1977"
      variable_units(iiPDF_eta)        = "kg/kg"

      variable_names(iiPDF_w)        = "w"
      variable_descriptions(iiPDF_w) = "Vertical velocity"
      variable_units(iiPDF_w)        = "m/s"

      if ( iiPDF_rrain > 0 ) then
        variable_names(iiPDF_rrain)        = "rrain"
        variable_descriptions(iiPDF_rrain) = "Rain water mixing ratio"
        variable_units(iiPDF_rrain)        = "kg/kg"
      end if
      if ( iiPDF_rice > 0 ) then
        variable_names(iiPDF_rice)        = "rice"
        variable_descriptions(iiPDF_rice) = "Ice water mixing ratio"
        variable_units(iiPDF_rice)        = "kg/kg"
      end if
      if ( iiPDF_rsnow > 0 ) then
        variable_names(iiPDF_rsnow)        = "rsnow"
        variable_descriptions(iiPDF_rsnow) = "Snow water mixing ratio"
        variable_units(iiPDF_rsnow)        = "kg/kg"
      end if
      if ( iiPDF_rgraupel > 0 ) then
        variable_names(iiPDF_rgraupel)        = "rgraupel"
        variable_descriptions(iiPDF_rgraupel) = "Graupel water mixing ratio"
        variable_units(iiPDF_rgraupel)        = "kg/kg"
      end if

      if ( iiPDF_Nr > 0 ) then
        variable_names(iiPDF_Nr)        = "Nr"
        variable_descriptions(iiPDF_Nr) = "Rain droplet number concentration"
        variable_units(iiPDF_Nr)        = "count/kg"
      end if
      if ( iiPDF_Ncn > 0 ) then
        variable_names(iiPDF_Ncn)        = "Ncn"
        variable_descriptions(iiPDF_Ncn) = "Cloud nuclei concentration (simplified)"
        variable_units(iiPDF_Ncn)        = "count/kg"
      end if
      if ( iiPDF_Ni > 0 ) then
        variable_names(iiPDF_Ni)        = "Ni"
        variable_descriptions(iiPDF_Ni) = "Ice number concentration"
        variable_units(iiPDF_Ni)        = "count/kg"
      end if
      if ( iiPDF_Nsnow > 0 ) then
        variable_names(iiPDF_Nsnow)        = "Nsnow"
        variable_descriptions(iiPDF_Nsnow) = "Snow number concentration"
        variable_units(iiPDF_Nsnow)        = "count/kg"
      end if
      if ( iiPDF_Ngraupel > 0 ) then
        variable_names(iiPDF_Ngraupel)        = "Ngraupel"
        variable_descriptions(iiPDF_Ngraupel) = "Graupel number concentration"
        variable_units(iiPDF_Ngraupel)        = "count/kg"
      end if

      i = d_variables + 1
      variable_names(i)        = "rt"
      variable_descriptions(i) = "Total water mixing ratio"
      variable_units(i)        = "kg/kg"

      i = d_variables + 2
      variable_names(i)        = "thl"
      variable_descriptions(i) = "Liquid potential temperature"
      variable_units(i)        = "K"

      call open_2D_samples_file( nz, lh_microphys_calls, d_variables+2, & ! In
                                 trim( fname_prefix )//"_nl", fdir, & ! In
                                 time_initial, stats_tout, zt, variable_names, & ! In
                                 variable_descriptions, variable_units, & ! In
                                 lognormal_sample_file ) ! In/Out

      deallocate( variable_names, variable_descriptions, variable_units )

    end if

    if ( l_output_2D_uniform_dist ) then

      allocate( variable_names(d_variables+3), variable_descriptions(d_variables+3), &
                variable_units(d_variables+3) )

      ! The uniform distribution corresponds to all the same variables as X_nl,
      ! except the d+1 component is the mixture component.

      variable_names(iiPDF_chi)        = "chi"
      variable_descriptions(iiPDF_chi) = "Uniform dist of the variable 's' from Mellor 1977"

      variable_names(iiPDF_eta)        = "eta"
      variable_descriptions(iiPDF_eta) = "Uniform dist of the variable 't' from Mellor 1977"

      variable_names(iiPDF_w)        = "w"
      variable_descriptions(iiPDF_w) = "Uniform dist of the vertical velocity"


      if ( iiPDF_rrain > 0 ) then
        variable_names(iiPDF_rrain)        = "rrain"
        variable_descriptions(iiPDF_rrain) = "Rain water mixing ratio"
        variable_units(iiPDF_rrain)        = "kg/kg"
      end if
      if ( iiPDF_rice > 0 ) then
        variable_names(iiPDF_rice)        = "rice"
        variable_descriptions(iiPDF_rice) = "Ice water mixing ratio"
        variable_units(iiPDF_rice)        = "kg/kg"
      end if
      if ( iiPDF_rsnow > 0 ) then
        variable_names(iiPDF_rsnow)        = "rsnow"
        variable_descriptions(iiPDF_rsnow) = "Snow water mixing ratio"
        variable_units(iiPDF_rsnow)        = "kg/kg"
      end if
      if ( iiPDF_rgraupel > 0 ) then
        variable_names(iiPDF_rgraupel)        = "rgraupel"
        variable_descriptions(iiPDF_rgraupel) = "Graupel water mixing ratio"
        variable_units(iiPDF_rgraupel)        = "kg/kg"
      end if

      if ( iiPDF_Nr > 0 ) then
        variable_names(iiPDF_Nr)        = "Nr"
        variable_descriptions(iiPDF_Nr) = "Rain droplet number concentration"
        variable_units(iiPDF_Nr)        = "count/kg"
      end if
      if ( iiPDF_Ncn > 0 ) then
        variable_names(iiPDF_Ncn)        = "Ncn"
        variable_descriptions(iiPDF_Ncn) = "Cloud nuclei concentration (simplified)"
        variable_units(iiPDF_Ncn)        = "count/kg"
      end if
      if ( iiPDF_Ni > 0 ) then
        variable_names(iiPDF_Ni)        = "Ni"
        variable_descriptions(iiPDF_Ni) = "Ice number concentration"
        variable_units(iiPDF_Ni)        = "count/kg"
      end if
      if ( iiPDF_Nsnow > 0 ) then
        variable_names(iiPDF_Nsnow)        = "Nsnow"
        variable_descriptions(iiPDF_Nsnow) = "Snow number concentration"
        variable_units(iiPDF_Nsnow)        = "count/kg"
      end if
      if ( iiPDF_Ngraupel > 0 ) then
        variable_names(iiPDF_Ngraupel)        = "Ngraupel"
        variable_descriptions(iiPDF_Ngraupel) = "Graupel number concentration"
        variable_units(iiPDF_Ngraupel)        = "count/kg"
      end if

      i = d_variables + 1
      variable_names(i) = "dp1"
      variable_descriptions(i) = "Uniform distribution for the mixture component"

      i = d_variables + 2
      variable_names(i) = "X_mixt_comp"
      variable_descriptions(i) = "Mixture component (should be 1 or 2)"

      i = d_variables + 3
      variable_names(i) = "p_matrix"
      variable_descriptions(i) = "P matrix elements at k_lh_start"

      ! Set all the units
      variable_units(:) = "count" ! Unidata units format for a dimensionless quantity

      call open_2D_samples_file( nz, lh_microphys_calls, i, & ! In
                                 trim( fname_prefix )//"_u", fdir, & ! In
                                 time_initial, stats_tout, zt, & ! In
                                 variable_names(1:i), variable_descriptions(1:i), & ! In
                                 variable_units(1:i), & ! In
                                 uniform_sample_file ) ! In/Out

      deallocate( variable_names, variable_descriptions, variable_units )

    end if

    return
  end subroutine latin_hypercube_2D_output

!-------------------------------------------------------------------------------
  subroutine latin_hypercube_2D_close
! Description:
!   Close a 2D sample file

! References:
!   None
!-------------------------------------------------------------------------------
    use output_2D_samples_module, only: &
      close_2D_samples_file ! Procedure

    use output_2D_samples_module, only: &
      lognormal_sample_file, & ! Variable(s)
      uniform_sample_file

    implicit none

    ! ---- Begin Code ----

    if ( l_output_2D_lognormal_dist ) then
      call close_2D_samples_file( lognormal_sample_file )
    end if
    if ( l_output_2D_uniform_dist ) then
      call close_2D_samples_file( uniform_sample_file )
    end if

    return
  end subroutine latin_hypercube_2D_close

!-------------------------------------------------------------------------------
  subroutine choose_X_u_reject &
             ( l_cloudy_sample, cloud_frac1, &
              cloud_frac2, mixt_frac, cloud_frac_thresh, &
              X_u_dp1_element, X_u_chi_element )

! Description:
!   Find a clear or cloudy point for sampling using the rejection method.
!
! References:
!   None
!-------------------------------------------------------------------------------
    use mt95, only: genrand_real3 ! Procedure

    use mt95, only: genrand_real ! Constant

    use constants_clubb, only: &
      fstderr ! Constant

    use clubb_precision, only: &
      core_rknd, & ! Variable(s)
      dp

    implicit none

    ! External
    intrinsic :: ceiling

    ! Constant parameters
    real( kind = core_rknd), parameter :: &
      max_iter_numerator = 100._core_rknd

    ! Input Variables
    logical, intent(in) :: &
      l_cloudy_sample ! Whether his is a cloudy or clear air sample point

    real( kind = core_rknd ), intent(in) :: &
      cloud_frac1, &    ! Cloud fraction associated with mixture component 1     [-]
      cloud_frac2, &    ! Cloud fraction associated with mixture component 2     [-]
      mixt_frac, &      ! Mixture fraction                                       [-]
      cloud_frac_thresh ! Minimum threshold for cloud fraction                   [-]

    ! Output Variables
    real(kind=dp), intent(out) :: &
      X_u_dp1_element, X_u_chi_element ! Elements from X_u (uniform dist.)

    ! Local Variables
    real(kind=core_rknd) :: cloud_frac_i

    real(kind=genrand_real) :: rand ! Random number

    ! Maximum iterations searching for the cloudy/clear part of the gridbox
    integer :: itermax

    integer :: i

!   integer :: X_mixt_comp_one_lev ! Whether we're in the first or second mixture component

    ! ---- Begin code ----

    ! Maximum iterations searching for the cloudy/clear part of the gridbox
    ! This should't appear in a parameter statement because it's set based on
    ! a floating-point calculation, and apparently that's not ISO Fortran
    itermax = ceiling( max_iter_numerator / cloud_frac_thresh )

    ! Find some new random numbers between (0,1)
    call genrand_real3( rand )
    X_u_dp1_element      = real(rand, kind = dp)
    call genrand_real3( rand )
    X_u_chi_element = real(rand, kind = dp)
    ! Here we use the rejection method to find a value in either the
    ! clear or cloudy part of the grid box
    do i = 1, itermax

      if ( in_mixt_comp_1( X_u_dp1_element, real(mixt_frac, kind = dp) ) ) then
        ! Component 1
        cloud_frac_i = cloud_frac1
!       X_mixt_comp_one_lev = 1
      else
        ! Component 2
        cloud_frac_i = cloud_frac2
!       X_mixt_comp_one_lev = 2
      end if

      if ( X_u_chi_element >= 1._dp-real(cloud_frac_i, kind = dp) .and. l_cloudy_sample ) then
        ! If we're looking for the cloudy part of the grid box, then exit this loop
        exit
      else if ( X_u_chi_element < ( 1._dp-real(cloud_frac_i, kind = dp) ) &
                .and. .not. l_cloudy_sample ) then
        ! If we're looking for the clear part of the grid box, then exit this loop
        exit
      else
        ! To prevent infinite loops we have this check here.
        ! Theoretically some seed might result in never picking the
        ! point we want after many iterations, but it's highly unlikely
        ! given that our current itermax is 100 / cloud_frac_thresh.
        ! -dschanen 19 March 2010
        if ( i == itermax ) then
          write(fstderr,*) "Maximum iteration reached in latin_hypercube driver."
          stop "Fatal error"
        else
          ! Find some new test values within the interval (0,1)
          call genrand_real3( rand )
          X_u_dp1_element      = real(rand, kind = dp)
          call genrand_real3( rand )
          X_u_chi_element = real(rand, kind = dp)
        end if
      end if ! Looking for a clear or cloudy point

    end do ! Loop until we either find what we want or reach itermax

    return
  end subroutine choose_X_u_reject

!-------------------------------------------------------------------------------
  subroutine choose_X_u_scaled &
             ( l_cloudy_sample, &
               p_matrix_element, num_samples, &
               cloud_frac1, cloud_frac2, &
               mixt_frac, &
               X_u_dp1_element, X_u_chi_element )

! Description:
!   Find a clear or cloudy point for sampling.
!
! References:
!   None
!-------------------------------------------------------------------------------
    use mt95, only: genrand_real3 ! Procedure

    use mt95, only: genrand_real ! Constant
    use mt95, only: r8 => genrand_real ! Constant

    use clubb_precision, only: &
      core_rknd, & ! Variable(s)
      dp

    implicit none

    ! Parameter Constants
    logical, parameter :: &
      l_use_p_matrix = .true.

    ! Input Variables
    logical, intent(in) :: &
      l_cloudy_sample ! Whether his is a cloudy or clear air sample point

    integer, intent(in) :: &
      p_matrix_element, & ! Integer from 0..num_samples for this sample
      num_samples       ! Total number of calls to the microphysics

    real( kind = core_rknd ), intent(in) :: &
      cloud_frac1, &    ! Cloud fraction associated with mixture component 1     [-]
      cloud_frac2, &    ! Cloud fraction associated with mixture component 2     [-]
      mixt_frac         ! Mixture fraction                                       [-]

    ! Output Variables
    real(kind=dp), intent(out) :: &
      X_u_dp1_element, X_u_chi_element ! Elements from X_u (uniform dist.)

    ! Local Variables
    real(kind=dp) :: cloud_frac_i, cloud_weighted_mixt_frac, clear_weighted_mixt_frac

    real(kind=genrand_real) :: rand, rand1, rand2 ! Random numbers

!   integer :: X_mixt_comp_one_lev

    ! ---- Begin code ----

    ! Pick a new mixture component value between (0,1)
    call genrand_real3( rand1 )

    call genrand_real3( rand2 ) ! Determine a 2nd rand for the if ... then

    if ( l_cloudy_sample ) then
      cloud_weighted_mixt_frac = real(mixt_frac*cloud_frac1, kind = dp) / &
                   real(mixt_frac*cloud_frac1 + (1._core_rknd-mixt_frac)*cloud_frac2, kind = dp)

      if ( in_mixt_comp_1( real( rand1, kind=dp ), cloud_weighted_mixt_frac ) ) then
        ! Component 1
        cloud_frac_i = real( cloud_frac1, kind=dp )
!       X_mixt_comp_one_lev = 1
        X_u_dp1_element = real( mixt_frac, kind=dp ) * rand2
      else
        ! Component 2
        cloud_frac_i = real( cloud_frac2, kind=dp )
!       X_mixt_comp_one_lev = 2
        X_u_dp1_element = real( mixt_frac, kind=dp ) &
                        + real(1._core_rknd-mixt_frac, kind=dp) * real(rand2, kind = dp)
      end if

      call genrand_real3( rand ) ! Rand between (0,1)

      ! Scale and translate sample point to reside in cloud
      if ( l_use_p_matrix ) then
        ! New formula based on p_matrix
        X_u_chi_element = 1._dp + 2._dp &
          * (real(p_matrix_element, kind = dp)/real(num_samples, kind = dp) &
          - 1._dp) * cloud_frac_i &
          + real(rand, kind = dp) * ( 2._dp/real( num_samples, kind=dp ) ) * cloud_frac_i
      else
        X_u_chi_element = cloud_frac_i * real(rand, kind = dp) &
          + (1._dp-cloud_frac_i)
      end if

    else ! Clear air sample
      clear_weighted_mixt_frac = ( ( 1._dp - real(cloud_frac1, kind = dp) ) &
        * real(mixt_frac, kind = dp) ) / ( ( 1._dp-real(cloud_frac1, kind = dp) ) &
        * real(mixt_frac, kind = dp) + ( 1._dp-real(cloud_frac2, kind = dp) )&
        *( 1._dp-real(mixt_frac, kind = dp) ) )

      if ( in_mixt_comp_1( real( rand1, kind=dp ), clear_weighted_mixt_frac ) ) then
        ! Component 1
        cloud_frac_i = real(cloud_frac1, kind = dp)
!       X_mixt_comp_one_lev = 1
        X_u_dp1_element = real(mixt_frac, kind = dp) * real(rand2, kind = dp)
      else
        ! Component 2
        cloud_frac_i = real(cloud_frac2, kind = dp)
!       X_mixt_comp_one_lev = 2
        X_u_dp1_element = real(mixt_frac, kind = dp) &
        + (1._dp-real(mixt_frac, kind = dp)) * real(rand2, kind = dp)
      end if

      call genrand_real3( rand ) ! Rand between (0,1)

      ! Scale and translate sample point to reside in clear air (no cloud)
      if ( l_use_p_matrix ) then
        ! New formula based on p_matrix
        X_u_chi_element = real( p_matrix_element, kind=dp ) &
          * (2._dp/real(num_samples, kind = dp) ) * (1._dp-cloud_frac_i) &
          + (2._dp/real(num_samples, kind = dp) * (1._dp-cloud_frac_i) &
          * real(rand, kind = dp))
      else
        X_u_chi_element = (1._dp-cloud_frac_i) &
          * real(rand, kind = dp)
      end if

    end if

    return
  end subroutine choose_X_u_scaled

!----------------------------------------------------------------------
  elemental function in_mixt_comp_1( X_u_dp1_element, frac )

! Description:
!   Determine if we're in mixture component 1

! References:
!   None
!----------------------------------------------------------------------

    use clubb_precision, only: &
      dp ! Variable(s)

    implicit none

    real(kind=dp), intent(in) :: &
      X_u_dp1_element, & ! Element of X_u telling us which mixture component we're in
      frac               ! The mixture fraction

    logical :: in_mixt_comp_1

    ! ---- Begin Code ----

    if ( X_u_dp1_element < frac ) then
      in_mixt_comp_1 = .true.
    else
      in_mixt_comp_1 = .false.
    end if

    return
  end function in_mixt_comp_1

!-------------------------------------------------------------------------------
  elemental function in_precipitation( rnd, precip_frac ) result( l_in_precip )

  ! Description:
  !   Determines if a sample is in precipitation

  ! References:
  !   None
  !-----------------------------------------------------------------------------

    use clubb_precision, only: dp

    implicit none

    ! Input Variables
    real( kind=dp ), intent(in) :: &
      rnd, &         ! Random number between 0 and 1
      precip_frac    ! Precipitation fraction

    ! Output Variable
    logical :: &
      l_in_precip    ! Whether the sample is in precipitation

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    if ( rnd < precip_frac ) then
      l_in_precip = .true.
    else
      l_in_precip = .false.
    end if

    return
  end function in_precipitation
!-------------------------------------------------------------------------------
  subroutine compute_arb_overlap( nz, k_lh_start, &
                                  X_u_one_var_k_lh_start, vert_corr, &
                                  X_u_one_var_all_levs )
! Description:
!   Re-computes X_u (uniform sample) for a single variate (e.g. chi) using
!   an arbitrary correlation specified by the input vert_corr variable (which
!   can vary with height).
!   This is an improved algorithm that doesn't require us to convert from a
!   unifrom distribution to a Gaussian distribution and back again.

! References:
!   None
!-------------------------------------------------------------------------------
    use mt95, only: &
      genrand_real ! Constant

    use mt95, only: &
      genrand_real3 ! Procedure

    use clubb_precision, only: &
      dp ! Variable(s)

    implicit none

    ! Parameter Constants

    ! Input Variables
    integer, intent(in) :: &
      nz,      & ! Number of vertical levels [-]
      k_lh_start   ! Starting k level          [-]

    real(kind=dp), intent(in) :: &
      X_u_one_var_k_lh_start  ! Uniform distribution of 1 variate (e.g. chi) at k_lh_start [-]

    real(kind=dp), dimension(nz), intent(in) :: &
      vert_corr ! Vertical correlation between k points in range [0,1]   [-]

    ! Output Variables
    real(kind=dp), dimension(nz), intent(out) :: &
      X_u_one_var_all_levs ! Uniform distribution of 1 variate at all levels [-]

    ! Local Variables
    real(kind=genrand_real) :: rand ! random number in the range (0,1)

    real(kind=dp) :: min_val, half_width, offset, unbounded_point

    integer :: k, kp1, km1 ! Loop iterators

    ! ---- Begin Code ----

    ! Set the value at k_lh_start using a special value we've selected
    X_u_one_var_all_levs(k_lh_start) = X_u_one_var_k_lh_start

    ! Upwards loop
    do k = k_lh_start, nz-1

      kp1 = k+1 ! This is the level we're computing

      if ( vert_corr(kp1) < 0._dp .or. vert_corr(kp1) > 1._dp ) then
        stop "vert_corr(kp1) not between 0 and 1"
      end if

      half_width = 1.0_dp - vert_corr(kp1)
      min_val = X_u_one_var_all_levs(k) - half_width

      call genrand_real3( rand ) ! (0,1)
      offset = 2.0_dp * half_width * real(rand, kind = dp)

      unbounded_point = min_val + offset

      ! If unbounded_point lies outside the range [0,1],
      ! fold it back so that it is between [0,1]
      if ( unbounded_point > 1.0_dp ) then
        X_u_one_var_all_levs(kp1) = 2.0_dp - unbounded_point
      else if ( unbounded_point < 0.0_dp ) then
        X_u_one_var_all_levs(kp1) = - unbounded_point
      else
        X_u_one_var_all_levs(kp1) = unbounded_point
      end if

    end do ! k_lh_start..nz-1

    ! Downwards loop
    do k = k_lh_start, 2, -1

      km1 = k-1 ! This is the level we're computing

      if ( vert_corr(km1) < 0._dp .or. vert_corr(km1) > 1._dp ) then
        stop "vert_corr(km1) not between 0 and 1"
      end if

      half_width = 1.0_dp - vert_corr(km1)
      min_val = X_u_one_var_all_levs(k) - half_width

      call genrand_real3( rand ) ! (0,1)
      offset = 2.0_dp * half_width * real(rand, kind = dp)

      unbounded_point = min_val + offset

      ! If unbounded_point lies outside the range [0,1],
      ! fold it back so that it is between [0,1]
      if ( unbounded_point > 1.0_dp ) then
        X_u_one_var_all_levs(km1) = 2.0_dp - unbounded_point
      else if ( unbounded_point < 0.0_dp ) then
        X_u_one_var_all_levs(km1) = - unbounded_point
      else
        X_u_one_var_all_levs(km1) = unbounded_point
      end if

    end do ! k_lh_start..2 decrementing

    return
  end subroutine compute_arb_overlap

!-------------------------------------------------------------------------------
  subroutine assert_check_half_cloudy &
             ( num_samples, cloud_frac1, &
               cloud_frac2, X_mixt_comp_k_lh_start, &
               X_u_chi_k_lh_start )
! Description:
!   Verify that half the points are in cloud if cloud weighted sampling is
!   enabled and other conditions are met. We stop rather than exiting gracefully
!   if this assertion check fails.

! References:
!   None.
!-------------------------------------------------------------------------------

    use constants_clubb, only: &
      fstderr ! Constant

    use clubb_precision, only: &
      core_rknd, & ! Variable(s)
      dp

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      num_samples ! Total calls to the microphysics

    real( kind = core_rknd ), intent(in) :: &
      cloud_frac1, cloud_frac2 ! Cloud fraction associatated with component 1 & 2 [-]

    integer, dimension(num_samples), intent(in) :: &
      X_mixt_comp_k_lh_start ! Mixture components at k_lh_start

    real(kind=dp), dimension(num_samples), intent(in) :: &
      X_u_chi_k_lh_start   ! Uniform distribution for chi at k_lh_start [-]

    ! Local Variables
    real(kind = core_rknd) :: cloud_frac_i

    integer :: number_cloudy_samples

    integer :: sample  ! Loop iterator

    ! ---- Begin Code ----

    number_cloudy_samples = 0

    do sample = 1, num_samples
      if ( X_mixt_comp_k_lh_start(sample) == 1 ) then
        cloud_frac_i = cloud_frac1
      else
        cloud_frac_i = cloud_frac2
      end if
      if ( X_u_chi_k_lh_start(sample) >= real(1._core_rknd-cloud_frac_i, kind = dp) ) then
        number_cloudy_samples = number_cloudy_samples + 1
      else
        ! Do nothing, the air is clear
      end if
    end do
    if ( number_cloudy_samples /= ( num_samples / 2 ) ) then
      write(fstderr,*) "Error, half of all samples aren't in cloud"
      write(fstderr,*) "X_u chi random = ", &
        X_u_chi_k_lh_start(:), "X_mixt_comp = ", X_mixt_comp_k_lh_start, &
        "cloudy samples =", number_cloudy_samples
      write(fstderr,*) "cloud_frac1 = ", cloud_frac1
      write(fstderr,*) "cloud_frac2 = ", cloud_frac2
      !stop "Fatal Error"
    end if

    return
  end subroutine assert_check_half_cloudy

!-------------------------------------------------------------------------------
  function compute_vert_corr( nz, delta_zm, Lscale_vert_avg ) result( vert_corr )
! Description:
!   This function computes the vertical correlation for arbitrary overlap, using
!   density weighted 3pt averaged Lscale and the difference in height levels
!   (delta_zm).
! References:
!   None
!-------------------------------------------------------------------------------

    use clubb_precision, only: &
      dp, & ! Variable(s)
      core_rknd

    implicit none

    ! External
    intrinsic :: exp

    ! Parameter Constants
    real( kind = dp ), parameter :: &
      vert_decorr_coef = 0.03_dp ! Empirically defined de-correlation constant [-]

    ! Input Variables
    integer, intent(in) :: &
      nz ! Number of vertical levels  [-]

    real( kind = core_rknd ), intent(in), dimension(nz) :: &
      delta_zm, &     ! Difference between altitudes    [m]
      Lscale_vert_avg ! Vertically averaged Lscale      [m]

    ! Output Variable
    real( kind = dp ), dimension(nz) :: &
      vert_corr ! The vertical correlation      [-]

    ! ---- Begin Code ----
    vert_corr(1:nz) = exp( -vert_decorr_coef &
                            * real( delta_zm(1:nz) / Lscale_vert_avg(1:nz), kind=dp ) )

    return
  end function compute_vert_corr

!-------------------------------------------------------------------------------
  subroutine stats_accumulate_lh &
             ( nz, num_samples, d_variables, rho_ds_zt, &
               lh_sample_point_weights, X_nl_all_levs, &
               lh_thl, lh_rt, Nc_in_cloud )

! Description:
!   Clip subcolumns from latin hypercube and create stats for diagnostic
!   purposes.

! References:
!   None
!-------------------------------------------------------------------------------

    use parameters_model, only: hydromet_dim ! Variable

    use grid_class, only: gr

    use parameters_microphys, only: &
        l_const_Nc_in_cloud ! Variable(s)

    use stats_variables, only: &
      l_stats_samp, & ! Variable(s)
      ilh_rrainm, &
      ilh_Nrm, &
      ilh_ricem, &
      ilh_Nim, &
      ilh_rsnowm, &
      ilh_Nsnowm, &
      ilh_rgraupelm, &
      ilh_Ngraupelm, &
      ilh_thlm, &
      ilh_rcm, &
      ilh_Ncm, &
      ilh_Ncnm, &
      ilh_rvm, &
      ilh_wm, &
      ilh_cloud_frac, &
      ilh_chi, &
      ilh_chip2, &
      ilh_eta

    use stats_variables, only: &
      ilh_wp2_zt, &  ! Variable(s)
      ilh_Nrp2_zt, &
      ilh_Ncp2_zt, &
      ilh_Ncnp2_zt, &
      ilh_rcp2_zt, &
      ilh_rtp2_zt, &
      ilh_thlp2_zt, &
      ilh_rrainp2_zt, &
      ilh_vwp, &
      ilh_lwp, &
      lh_zt, &
      lh_sfc

    use math_utilities, only: &
      compute_sample_mean, & ! Procedure(s)
      compute_sample_variance

    use stats_type_utilities, only: &
      stat_update_var, & ! Procedure(s)
      stat_update_var_pt

    use array_index, only: &
      iirrainm, & ! Variables
      iirsnowm, & 
      iiricem, & 
      iirgraupelm, & 
      iiNrm, &
      iiNsnowm, &
      iiNim, &
      iiNgraupelm

    use corr_matrix_module, only: &
      iiPDF_chi, & ! Variable(s)
      iiPDF_eta, &
      iiPDF_w, &
      iiPDF_Ncn

    use estimate_scm_microphys_module, only: &
      copy_X_nl_into_hydromet_all_pts ! Procedure(s)

    use constants_clubb, only: &
      zero_threshold    ! Constant(s)

    use clubb_precision, only: & 
      core_rknd, & ! Variable(s)
      dp

   use fill_holes, only: &
     vertical_integral ! Procedure(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      d_variables,     & ! Number of variables to sample
      num_samples,   & ! Number of calls to microphysics per timestep (normally=2)
      nz                 ! Number of vertical model levels

    real( kind = core_rknd ), intent(in), dimension(nz) :: &
      rho_ds_zt  ! Dry, static density (thermo. levs.) [kg/m^3]

    real( kind = core_rknd ), intent(in), dimension(num_samples) :: &
      lh_sample_point_weights

    real( kind = dp ), intent(in), dimension(nz,num_samples,d_variables) :: &
      X_nl_all_levs ! Sample that is transformed ultimately to normal-lognormal

    real( kind = core_rknd ), intent(in), dimension(nz,num_samples) :: &
      lh_thl, & ! Sample of liquid potential temperature [K]
      lh_rt     ! Sample of total water mixing ratio     [kg/kg]

    real( kind = core_rknd), dimension(nz), intent(in) :: &
      ! Constant value of N_c within cloud, to be used with l_const_Nc_in_cloud
      Nc_in_cloud

    ! Local variables
    real( kind = core_rknd ), dimension(nz,num_samples) :: &
      rc_all_points,  & ! Cloud water mixing ratio for all levels        [kg/kg]
      Nc_all_points,  & ! Cloud droplet conc. for all levels             [#/kg]
      Ncn_all_points, & ! Cloud nuclei conc. for all levels; Nc=Ncn*H(s) [#/kg]
      rv_all_points     ! Vapor mixing ratio for all levels              [kg/kg]

    real( kind = core_rknd ), dimension(nz,num_samples,hydromet_dim) :: &
      hydromet_all_points ! Hydrometeor species    [units vary]

    real( kind = core_rknd ), dimension(nz,hydromet_dim) :: &
      lh_hydromet ! Average value of the latin hypercube est. of all hydrometeors [units vary]

    real( kind = core_rknd ), dimension(nz) :: &
      lh_thlm,       & ! Average value of the latin hypercube est. of theta_l           [K]
      lh_rcm,        & ! Average value of the latin hypercube est. of rc                [kg/kg]
      lh_Ncm,        & ! Average value of the latin hypercube est. of Nc                [num/kg]
      lh_Ncnm,       & ! Average value of the latin hypercube est. of Ncn               [num/kg]
      lh_rvm,        & ! Average value of the latin hypercube est. of rv                [kg/kg]
      lh_wm,         & ! Average value of the latin hypercube est. of vertical velocity [m/s]
      lh_wp2_zt,     & ! Average value of the variance of the LH est. of vert. vel.     [m^2/s^2]
      lh_rrainp2_zt, & ! Average value of the variance of the LH est. of rrain.         [(kg/kg)^2]
      lh_rcp2_zt,    & ! Average value of the variance of the LH est. of rc.            [(kg/kg)^2]
      lh_rtp2_zt,    & ! Average value of the variance of the LH est. of rt             [kg^2/kg^2]
      lh_thlp2_zt,   & ! Average value of the variance of the LH est. of thetal         [K^2]
      lh_Nrp2_zt,    & ! Average value of the variance of the LH est. of Nr.            [#^2/kg^2]
      lh_Ncp2_zt,    & ! Average value of the variance of the LH est. of Nc.            [#^2/kg^2]
      lh_Ncnp2_zt,   & ! Average value of the variance of the LH est. of Ncn.           [#^2/kg^2]
      lh_cloud_frac, & ! Average value of the latin hypercube est. of cloud fraction    [-]
      lh_chi,   & ! Average value of the latin hypercube est. of Mellor's s        [kg/kg]
      lh_eta,   & ! Average value of the latin hypercube est. of Mellor's t        [kg/kg]
      lh_chip2           ! Average value of the variance of the LH est. of chi       [kg/kg]


    real(kind=core_rknd) :: xtmp

    integer :: sample, ivar, k

    ! ---- Begin Code ----

    ! Clip 's' from Mellor to obtain cloud-water mixing ratio
    rc_all_points = max( zero_threshold, real( X_nl_all_levs(:,:,iiPDF_chi), kind=core_rknd ) )

    if ( l_stats_samp ) then

      ! For all cases where l_lh_cloud_weighted_sampling is false, the weights
      ! will be 1 (all points equally weighted)

      if ( ilh_rcm + ilh_rcp2_zt + ilh_lwp > 0 ) then
        lh_rcm = compute_sample_mean( nz, num_samples, lh_sample_point_weights, &
                                      rc_all_points )
        call stat_update_var( ilh_rcm, lh_rcm, lh_zt )

        if ( ilh_lwp > 0 ) then
          xtmp &
          = vertical_integral &
               ( (gr%nz - 2 + 1), rho_ds_zt(2:gr%nz), &
                 lh_rcm(2:gr%nz), gr%invrs_dzt(2:gr%nz) )

          call stat_update_var_pt( ilh_lwp, 1, xtmp, lh_sfc )
        end if
      end if

      if ( ilh_thlm + ilh_thlp2_zt > 0 ) then
        lh_thlm = compute_sample_mean( nz, num_samples, lh_sample_point_weights, &
                                       real( lh_thl, kind = core_rknd ) )
        call stat_update_var( ilh_thlm, lh_thlm, lh_zt )
      end if

      if ( ilh_rvm + ilh_rtp2_zt > 0 ) then
        rv_all_points = lh_rt - rc_all_points
        lh_rvm = compute_sample_mean( nz, num_samples, lh_sample_point_weights, &
                                      rv_all_points )
        call stat_update_var( ilh_rvm, lh_rvm, lh_zt )
        if ( ilh_vwp > 0 ) then
          xtmp &
          = vertical_integral &
               ( (gr%nz - 2 + 1), rho_ds_zt(2:gr%nz), &
                 lh_rvm(2:gr%nz), gr%invrs_dzt(2:gr%nz) )

          call stat_update_var_pt( ilh_vwp, 1, xtmp, lh_sfc )
        end if
      end if

      if ( ilh_wm + ilh_wp2_zt > 0 ) then
        lh_wm  = compute_sample_mean( nz, num_samples, lh_sample_point_weights, &
                                      real( X_nl_all_levs(:,:,iiPDF_w), kind = core_rknd) )
        call stat_update_var( ilh_wm, lh_wm, lh_zt )
      end if

      if ( ilh_rrainm + ilh_Nrm + ilh_ricem + ilh_Nim + ilh_rsnowm + ilh_Nsnowm + &
           ilh_rgraupelm + ilh_Ngraupelm + ilh_Ncnm + ilh_Ncm > 0 ) then

        lh_hydromet = 0._core_rknd
        call copy_X_nl_into_hydromet_all_pts( nz, d_variables, num_samples, & ! In
                                      X_nl_all_levs, &  ! In
                                      lh_hydromet, & ! In
                                      hydromet_all_points, &  ! Out
                                      Ncn_all_points ) ! Out

        if ( l_const_Nc_in_cloud ) then
           ! For l_const_Nc_in_cloud, we want to use the same value of Nc for
           ! all sample points.
           do k = 1, nz, 1
              where ( rc_all_points(k,:) > 0.0_core_rknd )
                 Nc_all_points(k,:) = Nc_in_cloud(k)
              else where
                 Nc_all_points(k,:) = 0.0_core_rknd
              end where
           enddo ! k = 1, nz, 1
        else ! Nc varies in-cloud
           where ( rc_all_points > 0.0_core_rknd )
              Nc_all_points = Ncn_all_points
           else where
              Nc_all_points = 0.0_core_rknd
           end where
        endif

        ! Get rid of an annoying compiler warning.
        ivar = 1
        ivar = ivar

        forall ( ivar = 1:hydromet_dim )
          lh_hydromet(:,ivar) = compute_sample_mean( nz, num_samples, lh_sample_point_weights,&
                                                     hydromet_all_points(:,:,ivar) )
        end forall ! 1..hydromet_dim

      end if

      if ( ilh_Ncnm > 0 ) then
        lh_Ncnm = compute_sample_mean( nz, num_samples, lh_sample_point_weights, &
                                       Ncn_all_points(:,:) )
        call stat_update_var( ilh_Ncnm, lh_Ncnm, lh_zt )
      end if

      if ( ilh_Ncm > 0 ) then
        lh_Ncm = compute_sample_mean( nz, num_samples, lh_sample_point_weights, &
                                      Nc_all_points(:,:) )
        call stat_update_var( ilh_Ncm, lh_Ncm, lh_zt )
      end if

      ! Latin hypercube estimate of cloud fraction
      if ( ilh_cloud_frac > 0 ) then
        lh_cloud_frac(:) = 0._core_rknd
        do sample = 1, num_samples
          where ( X_nl_all_levs(:,sample,iiPDF_chi) > 0._dp )
            lh_cloud_frac(:) = lh_cloud_frac(:) + 1.0_core_rknd * lh_sample_point_weights(sample)
          end where
        end do
        lh_cloud_frac(:) = lh_cloud_frac(:) / real( num_samples, kind = core_rknd )

        call stat_update_var( ilh_cloud_frac, lh_cloud_frac, lh_zt )
      end if

      ! Latin hypercube estimate of chi
      if ( ilh_chi > 0 ) then
        lh_chi(1:nz) &
        = compute_sample_mean( nz, num_samples, lh_sample_point_weights, &
                               real( X_nl_all_levs(1:nz, 1:num_samples, iiPDF_chi), &
                                     kind = core_rknd ) )
        call stat_update_var( ilh_chi, lh_chi, lh_zt )
      end if

      ! Latin hypercube estimate of variance of chi
      if ( ilh_chip2 > 0 ) then
        lh_chip2(1:nz) &
        = compute_sample_variance( nz, num_samples, &
                                   real( X_nl_all_levs(:,:,iiPDF_chi), kind = core_rknd ), &
                                   lh_sample_point_weights, lh_chi(1:nz) )
        call stat_update_var( ilh_chip2, lh_chip2, lh_zt )
      end if

      ! Latin hypercube estimate of eta
      if ( ilh_eta > 0 ) then
        lh_eta(1:nz) &
        = compute_sample_mean( nz, num_samples, lh_sample_point_weights, &
                               real( X_nl_all_levs(1:nz, 1:num_samples, iiPDF_eta), &
                                     kind = core_rknd ) )
        call stat_update_var( ilh_eta, lh_eta, lh_zt )
      end if

      if ( ilh_wp2_zt > 0 ) then
        ! Compute the variance of vertical velocity
        lh_wp2_zt = compute_sample_variance( nz, num_samples, &
                                             real( X_nl_all_levs(:,:,iiPDF_w), kind = core_rknd ), &
                                             lh_sample_point_weights, lh_wm )
        call stat_update_var( ilh_wp2_zt, lh_wp2_zt, lh_zt )
      end if

      if ( ilh_rcp2_zt  > 0 ) then
        ! Compute the variance of cloud water mixing ratio
        lh_rcp2_zt = compute_sample_variance &
                     ( nz, num_samples, rc_all_points, &
                       lh_sample_point_weights, lh_rcm )
        call stat_update_var( ilh_rcp2_zt, lh_rcp2_zt, lh_zt )
      end if

      if ( ilh_rtp2_zt > 0 ) then
        ! Compute the variance of total water
        lh_rtp2_zt = compute_sample_variance &
                     ( nz, num_samples, &
                       real( lh_rt, kind = core_rknd ), lh_sample_point_weights, lh_rvm+lh_rcm )
        call stat_update_var( ilh_rtp2_zt, lh_rtp2_zt, lh_zt )
      end if

      if ( ilh_thlp2_zt > 0 ) then
        ! Compute the variance of liquid potential temperature
        lh_thlp2_zt = compute_sample_variance( nz, num_samples, &
                        real( lh_thl, kind = core_rknd ), lh_sample_point_weights, lh_thlm )
        call stat_update_var( ilh_thlp2_zt, lh_thlp2_zt, lh_zt )
      end if

      ! Compute the variance of rain water mixing ratio
      if ( iirrainm > 0 .and. ilh_rrainp2_zt > 0 ) then
        lh_rrainp2_zt = compute_sample_variance &
                        ( nz, num_samples, hydromet_all_points(:,:,iirrainm), &
                          lh_sample_point_weights, lh_hydromet(:,iirrainm) )
        call stat_update_var( ilh_rrainp2_zt, lh_rrainp2_zt, lh_zt )
      end if

      ! Compute the variance of cloud nuclei concentration (simplifed)
      if ( iiPDF_Ncn > 0 .and. ilh_Ncnp2_zt > 0 ) then
        lh_Ncnp2_zt = compute_sample_variance &
                      ( nz, num_samples, Ncn_all_points(:,:), &
                        lh_sample_point_weights, lh_Ncnm(:) )
        call stat_update_var( ilh_Ncnp2_zt, lh_Ncnp2_zt, lh_zt )
      end if

      ! Compute the variance of cloud droplet concentration
      if ( ilh_Ncp2_zt > 0 ) then
        lh_Ncp2_zt = compute_sample_variance &
                     ( nz, num_samples, Nc_all_points(:,:), &
                       lh_sample_point_weights, lh_Ncm(:) )
        call stat_update_var( ilh_Ncp2_zt, lh_Ncp2_zt, lh_zt )
      end if

      ! Compute the variance of rain droplet number concentration
      if ( iiNrm > 0 .and. ilh_Nrp2_zt > 0 ) then
        lh_Nrp2_zt = compute_sample_variance( nz, num_samples, hydromet_all_points(:,:,iiNrm),&
                                              lh_sample_point_weights, lh_hydromet(:,iiNrm) )
        call stat_update_var( ilh_Nrp2_zt, lh_Nrp2_zt, lh_zt )
      end if

      ! Averages of points being fed into the microphysics
      ! These are for diagnostic purposes, and are not needed for anything
      if ( iirrainm > 0 ) then
        call stat_update_var( ilh_rrainm, lh_hydromet(:,iirrainm), lh_zt )
      end if
      if ( iiNrm > 0 ) then
        call stat_update_var( ilh_Nrm, lh_hydromet(:,iiNrm), lh_zt )
      end if
      if ( iiricem > 0 ) then
        call stat_update_var( ilh_ricem, lh_hydromet(:,iiricem), lh_zt )
      end if
      if ( iiNim > 0 ) then
        call stat_update_var( ilh_Nim, lh_hydromet(:,iiNim), lh_zt )
      end if
      if ( iirsnowm > 0 ) then
        call stat_update_var( ilh_rsnowm, lh_hydromet(:,iirsnowm), lh_zt )
      end if
      if ( iiNsnowm > 0 ) then
        call stat_update_var( ilh_Nsnowm, lh_hydromet(:,iiNsnowm), lh_zt )
      end if
      if ( iirgraupelm > 0 ) then
        call stat_update_var( ilh_rgraupelm, lh_hydromet(:,iirgraupelm), lh_zt )
      end if
      if ( iiNgraupelm > 0 ) then
        call stat_update_var( ilh_Ngraupelm, lh_hydromet(:,iiNgraupelm), lh_zt )
      end if

    end if ! l_stats_samp

    return
  end subroutine stats_accumulate_lh

  !-----------------------------------------------------------------------
  subroutine stats_accumulate_uniform_lh( nz, num_samples, l_in_precip_all_levs, &
                                          X_mixt_comp_all_levs, &
                                          lh_sample_point_weights )

  ! Description:
  !   Samples statistics that cannot be deduced from the normal-lognormal
  !   SILHS sample (X_nl_all_levs)

  ! References:
  !   None
  !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd

    use stats_type_utilities, only: &
      stat_update_var

    use stats_variables, only: &
      l_stats_samp, &
      ilh_precip_frac, &
      ilh_mixt_frac, &
      lh_zt

    use math_utilities, only: &
      compute_sample_mean
      

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz, &         ! Number of vertical levels
      num_samples ! Number of SILHS sample points

    logical, dimension(nz,num_samples), intent(in) :: &
      l_in_precip_all_levs ! Boolean variables indicating whether a sample is in
                           ! precipitation at a given height level

    integer, dimension(nz,num_samples), intent(in) :: &
      X_mixt_comp_all_levs ! Integers indicating which mixture component a
                           ! sample is in at a given height level
    
    real( kind = core_rknd ), dimension(num_samples), intent(in) :: &
      lh_sample_point_weights ! The weight of each sample

    ! Local Variables
    real( kind = core_rknd ), dimension(nz) :: &
      lh_precip_frac, &
      lh_mixt_frac

    real( kind = core_rknd ), dimension(nz,num_samples) :: &
      int_in_precip, & ! '1' for samples in precipitation, '0' otherwise
      int_mixt_comp    ! '1' for samples in the first PDF component, '0' otherwise

  !-----------------------------------------------------------------------

    !----- Begin Code -----
    if ( l_stats_samp ) then
      ! Estimate of lh_precip_frac
      where ( l_in_precip_all_levs )
        int_in_precip = 1.0_core_rknd
      else where
        int_in_precip = 0.0_core_rknd
      end where
      lh_precip_frac(:) = compute_sample_mean( nz, num_samples, lh_sample_point_weights, &
                                               int_in_precip )
      call stat_update_var( ilh_precip_frac, lh_precip_frac, lh_zt )

      ! Estimate of lh_mixt_frac
      where ( X_mixt_comp_all_levs == 1 )
        int_mixt_comp = 1.0_core_rknd
      else where
        int_mixt_comp = 0.0_core_rknd
      end where
      lh_mixt_frac(:) = compute_sample_mean( nz, num_samples, lh_sample_point_weights, &
                                             int_mixt_comp )
      call stat_update_var( ilh_mixt_frac, lh_mixt_frac, lh_zt )

    end if ! l_stats_samp

    return
  end subroutine stats_accumulate_uniform_lh

#endif /* SILHS */

end module latin_hypercube_driver_module
