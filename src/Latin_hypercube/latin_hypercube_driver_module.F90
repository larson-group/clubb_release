! $Id$
!-------------------------------------------------------------------------------
module latin_hypercube_driver_module

  implicit none

  ! Constant Parameters
  logical, parameter, private :: &
    l_diagnostic_iter_check      = .true., &  ! Check for a problem in iteration
    l_output_2D_lognormal_dist   = .false., & ! Output a 2D netCDF file of the lognormal variates
    l_output_2D_uniform_dist     = .false.    ! Output a 2D netCDF file of the uniform distribution

  integer, allocatable, dimension(:,:,:), private :: & 
    height_time_matrix ! matrix of rand ints

  integer, private :: &
    prior_iter ! Prior iteration number (for diagnostic purposes)

  private ! Default scope

#ifdef LATIN_HYPERCUBE
  public :: LH_subcolumn_generator, LH_microphys_driver, latin_hypercube_2D_output, &
    latin_hypercube_2D_close

  contains

!-------------------------------------------------------------------------------
  subroutine LH_subcolumn_generator &
             ( iter, d_variables, n_micro_calls, sequence_length, nnzp, &
               thlm, pdf_params, wm_zt, delta_zm, rcm, rvm, &
               hydromet, xp2_on_xm2_array_cloud, xp2_on_xm2_array_below, &
               corr_array_cloud, corr_array_below, Lscale_vert_avg, &
               X_nl_all_levs, X_mixt_comp_all_levs, LH_rt, LH_thl, &
               LH_sample_point_weights )

! Description:
!   Call a microphysics scheme or generate an estimate of Kessler autoconversion
!   using latin hypercube sampling.
! References:
!   None
!-------------------------------------------------------------------------------

    use latin_hypercube_arrays, only: &
      iiLH_s_mellor   ! Variables

    use parameters_model, only: hydromet_dim ! Variable

    use permute_height_time_module, only: & 
      permute_height_time ! Procedure

    use generate_lh_sample_module, only: & 
      generate_lh_sample, & ! Procedure
      generate_uniform_sample

    use estimate_lh_micro_module, only: & 
      k_lh_start ! Variable

    use output_2D_samples_module, only: &
      output_2D_lognormal_dist_file, & ! Procedure(s)
      output_2D_uniform_dist_file

    use pdf_parameter_module, only: &
      pdf_parameter  ! Type

    use constants_clubb, only: &
      fstderr, & ! Constant
      cm3_per_m3

    use parameters_microphys, only: &
      l_lh_vert_overlap, &  ! Variables
      l_lh_cloud_weighted_sampling

    use error_code, only: &
      clubb_at_least_debug_level ! Procedure

    use mt95, only: genrand_real, genrand_intg ! Constants

    use mt95, only: genrand_init ! Procedure

    implicit none

    ! External
    intrinsic :: allocated, mod, maxloc, dble, epsilon

    ! Parameter Constants
    real, parameter :: &
      cloud_frac_thresh = 0.001 ! Threshold for sampling preferentially within cloud

    ! Find in and out of cloud points using the rejection method rather than scaling
    logical, parameter :: &
      l_use_rejection_method = .false., &
      l_re_seed = .true. ! Re-seed the Mersenne twister algorithm

    ! Input Variables
    integer, intent(in) :: &
      iter,            & ! Model iteration number
      d_variables,     & ! Number of variables to sample
      n_micro_calls,   & ! Number of calls to microphysics per timestep (normally=2)
      sequence_length, & ! nt_repeat/n_micro_call; number of timesteps before sequence repeats.
      nnzp               ! Number of vertical model levels

    type(pdf_parameter), dimension(nnzp), intent(in) :: & 
      pdf_params ! PDF parameters       [units vary]

    real, dimension(nnzp), intent(in) :: &
      thlm,      & ! Liquid potential temperature         [K]
      wm_zt,     & ! Mean w                             [m/s]
      delta_zm     ! Difference in moment. altitudes    [m]

    real, dimension(nnzp), intent(in) :: &
      rcm, & ! Liquid water mixing ratio        [kg/kg]
      rvm    ! Vapor water mixing ratio         [kg/kg]

    real, dimension(nnzp,hydromet_dim), intent(in) :: &
      hydromet ! Hydrometeor species    [units vary]

    real, dimension(d_variables), intent(in) :: &
      xp2_on_xm2_array_cloud, &! Variances over mean values squared [-]
      xp2_on_xm2_array_below

    real, dimension(d_variables,d_variables), intent(in) :: &
      corr_array_cloud, & ! Correlation for hydrometeor species [-]
      corr_array_below

    real, dimension(nnzp), intent(in) :: &
      Lscale_vert_avg ! 3pt vertical average of Lscale  [m]

    ! Output Variables
    double precision, intent(out), dimension(nnzp,n_micro_calls,d_variables) :: &
      X_nl_all_levs ! Sample that is transformed ultimately to normal-lognormal

    integer, intent(out), dimension(nnzp,n_micro_calls) :: &
      X_mixt_comp_all_levs ! Which mixture component we're in

    real, intent(out), dimension(nnzp,n_micro_calls) :: &
      LH_rt, LH_thl ! Sample of total water and liquid potential temperature [kg/kg],[K]

    real, intent(out), dimension(n_micro_calls) :: &
      LH_sample_point_weights

    ! Local variables

    real(kind=genrand_real), dimension(nnzp,n_micro_calls,(d_variables+1)) :: &
      X_u_all_levs ! Sample drawn from uniform distribution

    integer :: p_matrix(n_micro_calls,d_variables+1)

    real :: lh_start_cloud_frac ! Cloud fraction at k_lh_start [-]

    ! Try to obtain 12 digit accuracy for a diagnostic mean
    real(kind=selected_real_kind( p=12 ) ) :: mean_weight

    double precision :: mixt_frac_dp

    real(kind=genrand_real), dimension(n_micro_calls) :: &
      X_u_dp1_k_lh_start, X_u_s_mellor_k_lh_start

    real(kind=genrand_real), dimension(nnzp) :: &
      X_vert_corr ! Vertical correlation of a variate   [-]

    real(kind=genrand_real) :: X_u_temp

    ! Number of random samples before sequence of repeats (normally=10)
    integer :: nt_repeat

    integer :: &
      i_rmd, &   ! Remainder of ( iter-1 / sequence_length )
      k, sample  ! Loop iterators

    integer :: ivar ! Loop iterator

    logical :: l_cloudy_sample ! Whether a sample point is cloudy or clear air

    integer, dimension(1) :: tmp_loc

    ! ---- Begin Code ----

    nt_repeat = n_micro_calls * sequence_length

    if ( .not. allocated( height_time_matrix ) ) then
      ! If this is first time latin_hypercube_driver is called, then allocate
      ! the height_time_matrix and set the prior iteration number for debugging
      ! purposes.
      allocate( height_time_matrix(nnzp, nt_repeat, d_variables+1) )

      prior_iter = iter

      ! Re-seed
      if ( l_re_seed ) then
        ! This is the default seed in mt95.  Change this number to try
        ! non-default values.
        call genrand_init( put=5489_genrand_intg )
      end if

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
    if ( l_lh_cloud_weighted_sampling .and. mod( n_micro_calls, 2 ) /= 0 ) then
      write(fstderr,*) "Cloud weighted sampling requires micro calls to be divisible by 2."
      stop "Fatal error."
    end if
    if ( l_lh_cloud_weighted_sampling .and. sequence_length /= 1 ) then
      write(fstderr,*) "Cloud weighted sampling requires sequence length be equal to 1."
      stop "Fatal error."
    end if

    ! Initialize the sample point weights to 1.0
    LH_sample_point_weights(1:n_micro_calls)  = 1.0

    ! Latin hypercube sample generation
    ! Generate height_time_matrix, an nnzp x nt_repeat x d_variables array of random integers
    i_rmd = mod( iter-1, sequence_length )

    if ( i_rmd == 0 ) then
      call permute_height_time( nnzp, nt_repeat, d_variables+1, & ! intent(in)
                                height_time_matrix )              ! intent(out)
    end if
    ! End Latin hypercube sample generation

    ! print*, 'latin_hypercube_driver: i_rmd=', i_rmd

    !--------------------------------------------------------------
    ! Latin hypercube sampling
    !--------------------------------------------------------------


    ! For a 100 level fixed grid, this looks to be about the middle of the cloud for RICO
!   k_lh_start = 50
    tmp_loc    = maxloc( rcm )
    k_lh_start = tmp_loc(1) ! Attempt using the maximal value of rcm for now

    ! If there's no cloud k_lh_start appears to end up being 1.  Check if
    ! k_lh_start is 1 or nnzp and set it to the middle of the domain in that
    ! case.
    if ( k_lh_start == nnzp .or. k_lh_start == 1 ) then
      k_lh_start = nnzp / 2
    end if

    if ( l_lh_cloud_weighted_sampling ) then

      ! Determine cloud fraction at k_lh_start
      lh_start_cloud_frac = &
      pdf_params(k_lh_start)%mixt_frac * pdf_params(k_lh_start)%cloud_frac1 &
        + (1.0-pdf_params(k_lh_start)%mixt_frac) * pdf_params(k_lh_start)%cloud_frac2

      ! Determine p_matrix at k_lh_start
      p_matrix(1:n_micro_calls,1:(d_variables+1)) = &
        height_time_matrix(k_lh_start, n_micro_calls*i_rmd+1:n_micro_calls*i_rmd+n_micro_calls, &
                           1:d_variables+1)
    else
      lh_start_cloud_frac = -9999.0  ! This assignment eliminates a g95 compiler warning for an
      ! uninitialized variable. -meyern

    end if ! l_lh_cloud_weighted_sampling

    if ( l_lh_vert_overlap ) then

      ! Choose which rows of LH sample to feed into closure at the k_lh_start level
      p_matrix(1:n_micro_calls,1:(d_variables+1)) = &
        height_time_matrix(k_lh_start, n_micro_calls*i_rmd+1:n_micro_calls*i_rmd+n_micro_calls, &
                           1:d_variables+1)

      ! Generate the uniform distribution using the Mersenne twister at the k_lh_start level
      !  X_u has one extra dimension for the mixture component.
      call generate_uniform_sample( n_micro_calls, nt_repeat, d_variables+1, p_matrix, & ! In
                                    X_u_all_levs(k_lh_start,:,:) ) ! Out

      do sample = 1, n_micro_calls

        ! Get the uniform sample of d+1 (telling us if we're in component 1 or
        ! component 2), and the uniform sample for s_mellor
        X_u_dp1_k_lh_start(sample)      = X_u_all_levs(k_lh_start,sample,d_variables+1)
        X_u_s_mellor_k_lh_start(sample) = X_u_all_levs(k_lh_start,sample,iiLH_s_mellor)

        if ( l_lh_cloud_weighted_sampling ) then

          ! Save the cloud fraction as a weight for averaging preferentially
          ! within cloud
          if ( lh_start_cloud_frac >= 0.5 .or. lh_start_cloud_frac <= cloud_frac_thresh ) then
            LH_sample_point_weights(sample)  = 1.0
            ! There's no cloud or cloud fraction is >= 50%, so we do nothing

          else ! If we're in a partly cloud gridbox then we continue to the code below

            ! Detect which half of the sample points are in clear air and which half are in
            ! the cloudy air
!           if ( sample-1 < ( n_micro_calls / 2 ) ) then
!           if ( mod( sample, 2 ) == 0 ) then
            if ( p_matrix(sample,iiLH_s_mellor) < ( n_micro_calls / 2 ) ) then

              l_cloudy_sample = .false.
              LH_sample_point_weights(sample) = 2. * ( 1.0 - lh_start_cloud_frac )
            else
              l_cloudy_sample = .true.
              LH_sample_point_weights(sample) = 2. * lh_start_cloud_frac
            end if

            if ( l_use_rejection_method ) then
              ! Use the rejection method to select points that are in or out of cloud
              call choose_X_u_reject &
                   ( l_cloudy_sample, pdf_params(k_lh_start)%cloud_frac1, & ! In
                     pdf_params(k_lh_start)%cloud_frac2, pdf_params(k_lh_start)%mixt_frac, & !In
                     cloud_frac_thresh, & ! In
                     X_u_dp1_k_lh_start(sample), X_u_s_mellor_k_lh_start(sample) ) ! In/out

            else ! Transpose and scale the points to be in or out of cloud
              call choose_X_u_scaled &
                   ( l_cloudy_sample, & ! In
                     p_matrix(sample,iiLH_s_mellor), n_micro_calls, & ! In 
                     pdf_params(k_lh_start)%cloud_frac1, pdf_params(k_lh_start)%cloud_frac2, & ! In
                     pdf_params(k_lh_start)%mixt_frac, & !In
                     X_u_dp1_k_lh_start(sample), X_u_s_mellor_k_lh_start(sample) ) ! In/out

            end if

          end if ! Cloud fraction is between cloud_frac_thresh and 50%

        end if ! l_lh_cloud_weighted_sampling

      end do ! 1..n_micro_calls

      ! Use a fixed number for the vertical correlation.
!     X_vert_corr(1:nnzp) = 0.95_genrand_real

      ! Compute vertical correlation using a formula based on Lscale, the
      ! the difference in height levels, and an empirical constant
      X_vert_corr(1:nnzp) = &
        real( compute_vert_corr( nnzp, delta_zm, Lscale_vert_avg ), kind=genrand_real )

      ! Assertion check for the vertical correlation
      if ( clubb_at_least_debug_level( 1 ) ) then
        if ( any( X_vert_corr > 1.0 ) .or. any( X_vert_corr < 0.0 ) ) then
          write(fstderr,*) "The vertical correlation in latin_hypercube_driver"// &
            "is not in the correct range"
          do k = 1, nnzp
            write(fstderr,*) "k = ", k,  "Vert. correlation = ", X_vert_corr(k)
          end do
        end if ! Some correlation isn't between [0,1]
      end if ! clubb_at_least_debug_level 1

      do sample = 1, n_micro_calls
        ! Correlate s_mellor vertically
        call compute_arb_overlap &
             ( nnzp, k_lh_start, &  ! In
               X_u_s_mellor_k_lh_start(sample), X_vert_corr, & ! In
               X_u_all_levs(:,sample,iiLH_s_mellor) ) ! Out
        ! Correlate the d+1 variate vertically (used to compute the mixture
        ! component later)
        call compute_arb_overlap &
             ( nnzp, k_lh_start, &  ! In
               X_u_dp1_k_lh_start(sample), X_vert_corr, & ! In
               X_u_all_levs(:,sample,d_variables+1) ) ! Out

        ! Use these lines to make all variates vertically correlated, using the
        ! same correlation we used above for s_mellor and the d+1 variate
        do ivar = 1, d_variables
          if ( ivar /= iiLH_s_mellor ) then
            X_u_temp = X_u_all_levs(k_lh_start,sample,ivar)
            call compute_arb_overlap &
                 ( nnzp, k_lh_start, &  ! In
                   X_u_temp, X_vert_corr, & ! In
                   X_u_all_levs(:,sample,ivar) ) ! Out
          end if
        end do ! 1..d_variables
      end do ! 1..n_micro_calls
      ! %% Debug %%
      ! Testing what happens when we clip uniformally distributed variates to
      ! avoid extreme values.
!     where ( X_u_all_levs (:,:,2:d_variables) > 0.99 ) X_u_all_levs(:,:,2:d_variables) = 0.99
      ! %% Debug %%
    else ! Random overlap

      do k = 1, nnzp
        ! Choose which rows of LH sample to feed into closure.
        p_matrix(1:n_micro_calls,1:(d_variables+1)) = &
          height_time_matrix(k, n_micro_calls*i_rmd+1:n_micro_calls*i_rmd+n_micro_calls, &
                             1:d_variables+1)

        ! Generate the uniform distribution using the Mersenne twister
        !  X_u has one extra dimension for the mixture component.
        call generate_uniform_sample( n_micro_calls, nt_repeat, d_variables+1, p_matrix, & ! In
                                      X_u_all_levs(k,:,:) ) ! Out
      end do ! 1..nnzp

    end if ! l_lh_vert_overlap

    ! Determine mixture component for all levels
    do k = 1, nnzp

      mixt_frac_dp = dble( pdf_params(k)%mixt_frac )

      where ( in_mixt_comp_1(X_u_all_levs(k,:,d_variables+1), mixt_frac_dp ) )
        X_mixt_comp_all_levs(k,:) = 1
      else where
        X_mixt_comp_all_levs(k,:) = 2
      end where

    end do ! k = 1 .. nnzp

    ! Assertion check for whether half of sample points are cloudy.
    ! This is for the uniform sample only.  Another assertion check is in the
    ! estimate_lh_micro_module for X_nl_all_levs.
    if ( l_lh_cloud_weighted_sampling ) then
      if(clubb_at_least_debug_level( 2 ) .and. &
         lh_start_cloud_frac < 0.5 .and. lh_start_cloud_frac > cloud_frac_thresh ) then

        call assert_check_half_cloudy &
             ( n_micro_calls, pdf_params(k_lh_start)%cloud_frac1, &
               pdf_params(k_lh_start)%cloud_frac2, X_mixt_comp_all_levs(k_lh_start,:), &
               X_u_all_levs(k_lh_start,:,iiLH_s_mellor) )

      end if ! Maximal overlap, debug_level 2, and cloud-weighted averaging
    end if ! l_lh_cloud_weighted_sampling

    ! Assertion check to ensure that the sample point weights sum to approximately 1
    if ( l_lh_cloud_weighted_sampling .and. clubb_at_least_debug_level( 2 ) ) then
      mean_weight = 0.
      do sample = 1, n_micro_calls
        mean_weight = mean_weight + LH_sample_point_weights(sample)
      end do
      mean_weight = mean_weight / real( n_micro_calls )

      ! Using more precision for mean_weight should make this work out.
      ! The formula below could probably be redefined to estimate maximal ulps
      ! given the precision of LH_sample_point_weights and the number of
      ! n_micro_calls, but the formula below seems to be an ok approximation
      ! when we're using 4 or 8 byte precision floats.
      ! -dschanen 19 Nov 2010
      if ( abs( mean_weight - 1.0 ) > &
           real( n_micro_calls ) * epsilon( LH_sample_point_weights ) ) then
        write(fstderr,*) "Error in cloud weighted sampling code ", "mean_weight = ", mean_weight
        stop
      end if

    end if ! l_lh_cloud_weighted_sampling .and. clubb_at_least_debug_level( 2 )

    ! Upwards loop
    do k = k_lh_start, nnzp, 1
      ! Generate LH sample, represented by X_u and X_nl, for level k
      call generate_lh_sample &
           ( n_micro_calls, d_variables, hydromet_dim, &  ! In
             wm_zt(k), rcm(k), rvm(k), thlm(k), & ! In
             pdf_params(k)%mixt_frac, pdf_params(k)%rrtthl, & ! In
             pdf_params(k)%w1, pdf_params(k)%w2, & ! In
             pdf_params(k)%varnce_w1, pdf_params(k)%varnce_w2, & ! In
             pdf_params(k)%thl1, pdf_params(k)%thl2, & ! In
             pdf_params(k)%varnce_thl1, pdf_params(k)%varnce_thl2, & ! In
             pdf_params(k)%rt1, pdf_params(k)%rt2, & ! In
             pdf_params(k)%varnce_rt1, pdf_params(k)%varnce_rt2, & ! In
             pdf_params(k)%s1, pdf_params(k)%s2, & ! In
             pdf_params(k)%stdev_s1, pdf_params(k)%stdev_s2, & ! In
             pdf_params(k)%crt1, pdf_params(k)%crt2, & ! In
             pdf_params(k)%cthl1, pdf_params(k)%cthl2, & ! In
             hydromet(k,:), xp2_on_xm2_array_cloud, xp2_on_xm2_array_below, & ! In
             corr_array_cloud, corr_array_below, & ! In
             X_u_all_levs(k,:,:), X_mixt_comp_all_levs(k,:), & ! In
             LH_rt(k,:), LH_thl(k,:), X_nl_all_levs(k,:,:) ) ! Out
    end do ! k = k_lh_start..nnzp

    ! Downwards loop
    do k = k_lh_start-1, 1, -1
      call generate_lh_sample &
           ( n_micro_calls, d_variables, hydromet_dim, &  ! In
             wm_zt(k), rcm(k), rvm(k), thlm(k), & ! In
             pdf_params(k)%mixt_frac, pdf_params(k)%rrtthl, & ! In
             pdf_params(k)%w1, pdf_params(k)%w2, & ! In
             pdf_params(k)%varnce_w1, pdf_params(k)%varnce_w2, & ! In
             pdf_params(k)%thl1, pdf_params(k)%thl2, & ! In
             pdf_params(k)%varnce_thl1, pdf_params(k)%varnce_thl2, & ! In
             pdf_params(k)%rt1, pdf_params(k)%rt2, & ! In
             pdf_params(k)%varnce_rt1, pdf_params(k)%varnce_rt2, & ! In
             pdf_params(k)%s1, pdf_params(k)%s2, & ! In
             pdf_params(k)%stdev_s1, pdf_params(k)%stdev_s2, & ! In
             pdf_params(k)%crt1, pdf_params(k)%crt2, & ! In
             pdf_params(k)%cthl1, pdf_params(k)%cthl2, & ! In
             hydromet(k,:), xp2_on_xm2_array_cloud, xp2_on_xm2_array_below, & ! In
             corr_array_cloud, corr_array_below, & ! In
             X_u_all_levs(k,:,:), X_mixt_comp_all_levs(k,:), & ! In
             LH_rt(k,:), LH_thl(k,:), X_nl_all_levs(k,:,:) ) ! Out
    end do ! k_lh_start-1..1

    if ( l_output_2D_lognormal_dist ) then
      call output_2D_lognormal_dist_file( nnzp, n_micro_calls, d_variables, &
                                          X_nl_all_levs, LH_rt, LH_thl )
    end if
    if ( l_output_2D_uniform_dist ) then
      call output_2D_uniform_dist_file( nnzp, n_micro_calls, d_variables+1, &
                                        X_u_all_levs, X_mixt_comp_all_levs, &
                                        p_matrix )
    end if

    call LH_clipping_and_stats( nnzp, n_micro_calls, d_variables, & ! In
                                LH_sample_point_weights,  X_nl_all_levs, LH_thl, & ! In
                                LH_rt ) ! In/out

    return
  end subroutine LH_subcolumn_generator

  !=============================================================================
  subroutine LH_microphys_driver &
             ( dt, nnzp, n_micro_calls, d_variables, &
               X_nl_all_levs, LH_rt, LH_thl, LH_sample_point_weights, &
               pdf_params, p_in_Pa, exner, rho, &
               rcm, w_std_dev, delta_zt, cloud_frac, &
               hydromet, X_mixt_comp_all_levs, &
               LH_hydromet_mc, LH_hydromet_vel, &
               LH_rcm_mc, LH_rvm_mc, LH_thlm_mc, &
               microphys_sub )

    ! Description:
    !   Computes an estimate of the change due to microphysics given a set of
    !   subcolumns of thlm, rtm, et cetera from the subcolumn generator
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

    use estimate_lh_micro_module, only: & 
      estimate_lh_micro ! Procedure(s)

    implicit none

    ! Interface block
#include "../microphys_interface.inc"

    ! Input Variables
    real, intent(in) :: &
      dt ! Model timestep       [s]

    integer, intent(in) :: &
      d_variables,     & ! Number of variables to sample
      n_micro_calls,   & ! Number of calls to microphysics per timestep (normally=2)
      nnzp               ! Number of vertical model levels

    ! Input Variables
    double precision, intent(in), dimension(nnzp,n_micro_calls,d_variables) :: &
      X_nl_all_levs ! Sample that is transformed ultimately to normal-lognormal

    integer, intent(in), dimension(nnzp,n_micro_calls) :: &
      X_mixt_comp_all_levs ! Which mixture component we're in

    real, intent(in), dimension(nnzp,n_micro_calls) :: &
      LH_rt, LH_thl ! Sample of total water and liquid potential temperature [kg/kg],[K]

    real, intent(in), dimension(n_micro_calls) :: &
      LH_sample_point_weights ! Weight given the individual sample points

    type(pdf_parameter), dimension(nnzp), intent(in) :: & 
      pdf_params ! PDF parameters       [units vary]

    real, dimension(nnzp,hydromet_dim), intent(in) :: &
      hydromet ! Hydrometeor species    [units vary]

    real, dimension(nnzp), intent(in) :: &
      cloud_frac,  & ! Cloud fraction               [-]
      w_std_dev,   & ! Standard deviation of w      [m/s]
      delta_zt,    & ! Change in meters with height [m]
      rcm,         & ! Liquid water mixing ratio    [kg/kg]
      p_in_Pa,     & ! Pressure                     [Pa]
      exner,       & ! Exner function               [-]
      rho            ! Density on thermo. grid      [kg/m^3]

    ! Input/Output Variables
    real, dimension(nnzp,hydromet_dim), intent(inout) :: &
      LH_hydromet_mc, & ! LH estimate of hydrometeor time tendency          [(units vary)/s]
      LH_hydromet_vel   ! LH estimate of hydrometeor sedimentation velocity [m/s]

    ! Output Variables
    real, dimension(nnzp), intent(out) :: &
      LH_rcm_mc, & ! LH estimate of time tendency of liquid water mixing ratio    [kg/kg/s]
      LH_rvm_mc, & ! LH estimate of time tendency of vapor water mixing ratio     [kg/kg/s]
      LH_thlm_mc   ! LH estimate of time tendency of liquid potential temperature [K/s]

    ! ---- Begin Code ----

    ! Perform LH and analytic microphysical calculations
    call estimate_lh_micro &
         ( dt, nnzp, n_micro_calls, d_variables, &  ! intent(in)
           X_nl_all_levs, &                         ! intent(in)
           LH_rt, LH_thl, pdf_params, &             ! intent(in)
           p_in_Pa, exner, rho, &                   ! intent(in)
           rcm, w_std_dev, delta_zt, &              ! intent(in)
           cloud_frac, hydromet, &                  ! intent(in)
           X_mixt_comp_all_levs, LH_sample_point_weights, & ! intent(in)
           LH_hydromet_mc, LH_hydromet_vel, &       ! intent(inout)
           LH_rcm_mc, LH_rvm_mc, LH_thlm_mc, &      ! intent(out)
           lh_AKm, AKm, AKstd, AKstd_cld, &         ! intent(out)
           AKm_rcm, AKm_rcc, LH_rcm_avg, &          ! intent(out)
           microphys_sub )  ! Procedure

    return
  end subroutine LH_microphys_driver
!-------------------------------------------------------------------------------
  subroutine latin_hypercube_2D_output &
             ( fname_prefix, fdir, stats_tout, nnzp, &
               zt, time_initial )
!-------------------------------------------------------------------------------

    use latin_hypercube_arrays, only: &
      iiLH_s_mellor, & ! Variables
      iiLH_t_mellor, &
      iiLH_w, &
      iiLH_rrain, & 
      iiLH_rice, &
      iiLH_rsnow, &
      iiLH_rgraupel, &
      iiLH_Nr, &
      iiLH_Ni, &
      iiLH_Nsnow, &
      iiLH_Ngraupel, &
      iiLH_Nc

    use parameters_microphys, only: &
      LH_microphys_calls ! Variable

    use stats_precision, only: &
      time_precision ! Constant

    use output_2D_samples_module, only: &
      open_2D_samples_file ! Procedure

    use output_2D_samples_module, only: &
      lognormal_sample_file, & ! Instance of a type
      uniform_sample_file

    use latin_hypercube_arrays, only: &
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
      nnzp ! Number of vertical levels

    real, dimension(nnzp), intent(in) :: &
      zt ! Altitudes [m]

    ! Local Variables
    character(len=100), allocatable, dimension(:) :: &
      variable_names, variable_descriptions, variable_units

    integer :: i

    ! ---- Begin Code ----

    if ( l_output_2D_lognormal_dist ) then

      allocate( variable_names(d_variables+2), variable_descriptions(d_variables+2), &
                variable_units(d_variables+2) )

      variable_names(iiLH_s_mellor)        = "s_mellor"
      variable_descriptions(iiLH_s_mellor) = "The variable 's' from Mellor 1977"
      variable_units(iiLH_s_mellor)        = "kg/kg"

      variable_names(iiLH_t_mellor)        = "t_mellor"
      variable_descriptions(iiLH_t_mellor) = "The variable 't' from Mellor 1977"
      variable_units(iiLH_t_mellor)        = "kg/kg"

      variable_names(iiLH_w)        = "w"
      variable_descriptions(iiLH_w) = "Vertical velocity"
      variable_units(iiLH_w)        = "m/s"

      if ( iiLH_rrain > 0 ) then
        variable_names(iiLH_rrain)        = "rrain"
        variable_descriptions(iiLH_rrain) = "Rain water mixing ratio"
        variable_units(iiLH_rrain)        = "kg/kg"
      end if
      if ( iiLH_rice > 0 ) then
        variable_names(iiLH_rice)        = "rice"
        variable_descriptions(iiLH_rice) = "Ice water mixing ratio"
        variable_units(iiLH_rice)        = "kg/kg"
      end if
      if ( iiLH_rsnow > 0 ) then
        variable_names(iiLH_rsnow)        = "rsnow"
        variable_descriptions(iiLH_rsnow) = "Snow water mixing ratio"
        variable_units(iiLH_rsnow)        = "kg/kg"
      end if
      if ( iiLH_rgraupel > 0 ) then
        variable_names(iiLH_rgraupel)        = "rgraupel"
        variable_descriptions(iiLH_rgraupel) = "Graupel water mixing ratio"
        variable_units(iiLH_rgraupel)        = "kg/kg"
      end if

      if ( iiLH_Nr > 0 ) then
        variable_names(iiLH_Nr)        = "Nr"
        variable_descriptions(iiLH_Nr) = "Rain droplet number concentration"
        variable_units(iiLH_Nr)        = "count/kg"
      end if
      if ( iiLH_Nc > 0 ) then
        variable_names(iiLH_Nc)        = "Nc"
        variable_descriptions(iiLH_Nc) = "Cloud droplet number concentration"
        variable_units(iiLH_Nc)        = "count/kg"
      end if
      if ( iiLH_Ni > 0 ) then
        variable_names(iiLH_Ni)        = "Ni"
        variable_descriptions(iiLH_Ni) = "Ice number concentration"
        variable_units(iiLH_Ni)        = "count/kg"
      end if
      if ( iiLH_Nsnow > 0 ) then
        variable_names(iiLH_Nsnow)        = "Nsnow"
        variable_descriptions(iiLH_Nsnow) = "Snow number concentration"
        variable_units(iiLH_Nsnow)        = "count/kg"
      end if
      if ( iiLH_Ngraupel > 0 ) then
        variable_names(iiLH_Ngraupel)        = "Ngraupel"
        variable_descriptions(iiLH_Ngraupel) = "Graupel number concentration"
        variable_units(iiLH_Ngraupel)        = "count/kg"
      end if

      i = d_variables + 1
      variable_names(i)        = "rt"
      variable_descriptions(i) = "Total water mixing ratio"
      variable_units(i)        = "kg/kg"

      i = d_variables + 2
      variable_names(i)        = "thl"
      variable_descriptions(i) = "Liquid potential temperature"
      variable_units(i)        = "K"

      call open_2D_samples_file( nnzp, LH_microphys_calls, d_variables+2, & ! In
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

      variable_names(iiLH_s_mellor)        = "s_mellor"
      variable_descriptions(iiLH_s_mellor) = "Uniform dist of the variable 's' from Mellor 1977"

      variable_names(iiLH_t_mellor)        = "t_mellor"
      variable_descriptions(iiLH_t_mellor) = "Uniform dist of the variable 't' from Mellor 1977"

      variable_names(iiLH_w)        = "w"
      variable_descriptions(iiLH_w) = "Uniform dist of the vertical velocity"


      if ( iiLH_rrain > 0 ) then
        variable_names(iiLH_rrain)        = "rrain"
        variable_descriptions(iiLH_rrain) = "Rain water mixing ratio"
        variable_units(iiLH_rrain)        = "kg/kg"
      end if
      if ( iiLH_rice > 0 ) then
        variable_names(iiLH_rice)        = "rice"
        variable_descriptions(iiLH_rice) = "Ice water mixing ratio"
        variable_units(iiLH_rice)        = "kg/kg"
      end if
      if ( iiLH_rsnow > 0 ) then
        variable_names(iiLH_rsnow)        = "rsnow"
        variable_descriptions(iiLH_rsnow) = "Snow water mixing ratio"
        variable_units(iiLH_rsnow)        = "kg/kg"
      end if
      if ( iiLH_rgraupel > 0 ) then
        variable_names(iiLH_rgraupel)        = "rgraupel"
        variable_descriptions(iiLH_rgraupel) = "Graupel water mixing ratio"
        variable_units(iiLH_rgraupel)        = "kg/kg"
      end if

      if ( iiLH_Nr > 0 ) then
        variable_names(iiLH_Nr)        = "Nr"
        variable_descriptions(iiLH_Nr) = "Rain droplet number concentration"
        variable_units(iiLH_Nr)        = "count/kg"
      end if
      if ( iiLH_Nc > 0 ) then
        variable_names(iiLH_Nc)        = "Nc"
        variable_descriptions(iiLH_Nc) = "Cloud droplet number concentration"
        variable_units(iiLH_Nc)        = "count/kg"
      end if
      if ( iiLH_Ni > 0 ) then
        variable_names(iiLH_Ni)        = "Ni"
        variable_descriptions(iiLH_Ni) = "Ice number concentration"
        variable_units(iiLH_Ni)        = "count/kg"
      end if
      if ( iiLH_Nsnow > 0 ) then
        variable_names(iiLH_Nsnow)        = "Nsnow"
        variable_descriptions(iiLH_Nsnow) = "Snow number concentration"
        variable_units(iiLH_Nsnow)        = "count/kg"
      end if
      if ( iiLH_Ngraupel > 0 ) then
        variable_names(iiLH_Ngraupel)        = "Ngraupel"
        variable_descriptions(iiLH_Ngraupel) = "Graupel number concentration"
        variable_units(iiLH_Ngraupel)        = "count/kg"
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

      call open_2D_samples_file( nnzp, LH_microphys_calls, i, & ! In
                                 trim( fname_prefix )//"_u", fdir, & ! In
                                 time_initial, stats_tout, zt, &! In
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
              X_u_dp1_element, X_u_s_mellor_element )

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

    implicit none

    ! External
    intrinsic :: ceiling, dble

    ! Input Variables
    logical, intent(in) :: &
      l_cloudy_sample ! Whether his is a cloudy or clear air sample point

    real, intent(in) :: &
      cloud_frac1, &    ! Cloud fraction associated with mixture component 1     [-]
      cloud_frac2, &    ! Cloud fraction associated with mixture component 2     [-]
      mixt_frac, &      ! Mixture fraction                                       [-]
      cloud_frac_thresh ! Minimum threshold for cloud fraction                   [-]

    ! Input/Output Variables
    real(kind=genrand_real), intent(inout) :: &
      X_u_dp1_element, X_u_s_mellor_element ! Elements from X_u (uniform dist.)

    ! Local Variables
    real :: cloud_frac_n

    real(kind=genrand_real) :: rand ! Random number

    ! Maximum iterations searching for the cloudy/clear part of the gridbox
    integer :: itermax

    integer :: i

!   integer :: X_mixt_comp_one_lev ! Whether we're in the first or second mixture component

    ! ---- Begin code ----

    ! Maximum iterations searching for the cloudy/clear part of the gridbox
    ! This should't appear in a parameter statement because it's set based on
    ! a floating-point calculation, and apparently that's not ISO Fortran
    itermax = ceiling( 100. / cloud_frac_thresh )

    ! Find some new random numbers between (0,1)
    call genrand_real3( rand )
    X_u_dp1_element      = rand
    call genrand_real3( rand )
    X_u_s_mellor_element = rand
    ! Here we use the rejection method to find a value in either the
    ! clear or cloudy part of the grid box
    do i = 1, itermax

      if ( in_mixt_comp_1( X_u_dp1_element, real( mixt_frac, kind=genrand_real ) ) ) then
        ! Component 1
        cloud_frac_n = cloud_frac1
!       X_mixt_comp_one_lev = 1
      else
        ! Component 2
        cloud_frac_n = cloud_frac2
!       X_mixt_comp_one_lev = 2
      end if

      if ( X_u_s_mellor_element >= (1.-cloud_frac_n) .and. l_cloudy_sample ) then
        ! If we're looking for the cloudy part of the grid box, then exit this loop
        exit
      else if ( X_u_s_mellor_element < (1.-cloud_frac_n) & 
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
          X_u_dp1_element      = rand
          call genrand_real3( rand )
          X_u_s_mellor_element = rand
        end if
      end if ! Looking for a clear or cloudy point

    end do ! Loop until we either find what we want or reach itermax

    return
  end subroutine choose_X_u_reject

!-------------------------------------------------------------------------------
  subroutine choose_X_u_scaled &
             ( l_cloudy_sample, &
               p_matrix_element, n_micro_calls, &
               cloud_frac1, cloud_frac2, &
               mixt_frac, &
               X_u_dp1_element, X_u_s_mellor_element )

! Description:
!   Find a clear or cloudy point for sampling.
!
! References:
!   None
!-------------------------------------------------------------------------------
    use mt95, only: genrand_real3 ! Procedure

    use mt95, only: genrand_real ! Constant
    use mt95, only: r8 => genrand_real ! Constant

    use constants_clubb, only: &
      fstderr ! Constant

    implicit none

    ! External
    intrinsic :: dble

    ! Parameter Constants
    logical, parameter :: &
      l_use_p_matrix = .true.

    ! Input Variables
    logical, intent(in) :: &
      l_cloudy_sample ! Whether his is a cloudy or clear air sample point

    integer, intent(in) :: &
      p_matrix_element, & ! Integer from 0..n_micro_calls for this sample
      n_micro_calls       ! Total number of calls to the microphysics

    real, intent(in) :: &
      cloud_frac1, &    ! Cloud fraction associated with mixture component 1     [-]
      cloud_frac2, &    ! Cloud fraction associated with mixture component 2     [-]
      mixt_frac         ! Mixture fraction                                       [-]

    ! Input/Output Variables
    real(kind=genrand_real), intent(inout) :: &
      X_u_dp1_element, X_u_s_mellor_element ! Elements from X_u (uniform dist.)

    ! Local Variables
    real(kind=genrand_real) :: cloud_frac_n, cloud_weighted_mixt_frac, clear_weighted_mixt_frac

    real(kind=genrand_real) :: rand, rand1, rand2 ! Random numbers

!   integer :: X_mixt_comp_one_lev

    ! ---- Begin code ----

    ! Pick a new mixture component value between (0,1)
    call genrand_real3( rand1 )

    call genrand_real3( rand2 ) ! Determine a 2nd rand for the if ... then

    if ( l_cloudy_sample ) then
      cloud_weighted_mixt_frac = ( mixt_frac*cloud_frac1 ) / &
                   ( mixt_frac*cloud_frac1 + (1._r8-mixt_frac)*cloud_frac2 )

      if ( in_mixt_comp_1( rand1, real( cloud_weighted_mixt_frac, kind=genrand_real ) ) ) then
        ! Component 1
        cloud_frac_n = cloud_frac1
!       X_mixt_comp_one_lev = 1
        X_u_dp1_element = mixt_frac * rand2
      else
        ! Component 2
        cloud_frac_n = cloud_frac2
!       X_mixt_comp_one_lev = 2
        X_u_dp1_element = mixt_frac + (1._r8-mixt_frac) * rand2
      end if

      call genrand_real3( rand ) ! Rand between (0,1)

      ! Scale and translate sample point to reside in cloud
      if ( l_use_p_matrix ) then
        ! New formula based on p_matrix
        X_u_s_mellor_element = 1._r8 &
          + 2._r8 * (real( p_matrix_element )/real( n_micro_calls ) - 1._r8) * cloud_frac_n &
          + rand * ( 2._r8/real( n_micro_calls ) ) * cloud_frac_n
      else
        X_u_s_mellor_element = cloud_frac_n * rand + (1._r8-cloud_frac_n)
      end if

    else ! Clear air sample
      clear_weighted_mixt_frac = ( ( 1._r8 - cloud_frac1 ) * mixt_frac ) &
        / ( ( 1._r8-cloud_frac1 ) * mixt_frac + ( 1._r8-cloud_frac2 )*( 1._r8-mixt_frac ) )

      if ( in_mixt_comp_1( rand1, real( clear_weighted_mixt_frac, kind=genrand_real ) ) ) then
        ! Component 1
        cloud_frac_n = cloud_frac1
!       X_mixt_comp_one_lev = 1
        X_u_dp1_element = mixt_frac * rand2
      else
        ! Component 2
        cloud_frac_n = cloud_frac2
!       X_mixt_comp_one_lev = 2
        X_u_dp1_element = mixt_frac + (1._r8-mixt_frac) * rand2
      end if

      call genrand_real3( rand ) ! Rand between (0,1)

      ! Scale and translate sample point to reside in clear air (no cloud)
      if ( l_use_p_matrix ) then
        ! New formula based on p_matrix
        X_u_s_mellor_element = real( p_matrix_element ) &
            * (2._r8/real(n_micro_calls) ) * (1._r8-cloud_frac_n) &
          + (2._r8/real(n_micro_calls)) * (1._r8-cloud_frac_n) * rand
      else
        X_u_s_mellor_element = (1._r8-cloud_frac_n) * rand
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

    use mt95, only: genrand_real ! Constant

    implicit none

    real(kind=genrand_real), intent(in) :: &
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
  subroutine compute_arb_overlap( nnzp, k_lh_start, &
                                  X_u_one_var_k_lh_start, vert_corr, &
                                  X_u_one_var_all_levs )
! Description:
!   Re-computes X_u (uniform sample) for a single variate (e.g. s_mellor) using
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

    implicit none

    ! Parameter Constants

    ! Input Variables
    integer, intent(in) :: &
      nnzp,      & ! Number of vertical levels [-]
      k_lh_start   ! Starting k level          [-]

    real(kind=genrand_real), intent(in) :: &
      X_u_one_var_k_lh_start  ! Uniform distribution of 1 variate (e.g. s_mellor) at k_lh_start [-]

    real(kind=genrand_real), dimension(nnzp), intent(in) :: &
      vert_corr ! Vertical correlation between k points in range [0,1]   [-]

    ! Output Variables
    real(kind=genrand_real), dimension(nnzp), intent(out) :: &
      X_u_one_var_all_levs ! Uniform distribution of 1 variate at all levels [-]

    ! Local Variables
    real(kind=genrand_real) :: rand ! random number in the range (0,1)

    real(kind=genrand_real) :: min_val, half_width, offset, unbounded_point

    integer :: k, kp1, km1 ! Loop iterators

    ! ---- Begin Code ----

    ! Set the value at k_lh_start using a special value we've selected
    X_u_one_var_all_levs(k_lh_start) = X_u_one_var_k_lh_start

    ! Upwards loop
    do k = k_lh_start, nnzp-1

      kp1 = k+1 ! This is the level we're computing

      if ( vert_corr(kp1) < 0. .or. vert_corr(kp1) > 1. ) then
        stop "vert_corr(kp1) not between 0 and 1"
      end if

      half_width = 1.0 - vert_corr(kp1)
      min_val = X_u_one_var_all_levs(k) - half_width

      call genrand_real3( rand ) ! (0,1)
      offset = 2.0 * half_width * rand

      unbounded_point = min_val + offset

      ! If unbounded_point lies outside the range [0,1],
      ! fold it back so that it is between [0,1]
      if ( unbounded_point > 1.0 ) then
        X_u_one_var_all_levs(kp1) = 2.0 - unbounded_point
      else if ( unbounded_point < 0.0 ) then
        X_u_one_var_all_levs(kp1) = - unbounded_point
      else
        X_u_one_var_all_levs(kp1) = unbounded_point
      end if

!     print *, k, X_u_one_var_all_levs(k), kp1, X_u_one_var_all_levs(kp1)

    end do ! k_lh_start..nnzp-1

    ! Downwards loop
    do k = k_lh_start, 2, -1

      km1 = k-1 ! This is the level we're computing

      if ( vert_corr(km1) < 0. .or. vert_corr(km1) > 1. ) then
        stop "vert_corr(km1) not between 0 and 1"
      end if

      half_width = 1.0 - vert_corr(km1)
      min_val = X_u_one_var_all_levs(k) - half_width

      call genrand_real3( rand ) ! (0,1)
      offset = 2.0 * half_width * rand

      unbounded_point = min_val + offset

      ! If unbounded_point lies outside the range [0,1],
      ! fold it back so that it is between [0,1]
      if ( unbounded_point > 1.0 ) then
        X_u_one_var_all_levs(km1) = 2.0 - unbounded_point
      else if ( unbounded_point < 0.0 ) then
        X_u_one_var_all_levs(km1) = - unbounded_point
      else
        X_u_one_var_all_levs(km1) = unbounded_point
      end if

!     print *, k, X_u_one_var_all_levs(k), km1, X_u_one_var_all_levs(km1)

    end do ! k_lh_start..2 decrementing

    return
  end subroutine compute_arb_overlap

!-------------------------------------------------------------------------------
  subroutine assert_check_half_cloudy &
             ( n_micro_calls, cloud_frac1, &
               cloud_frac2, X_mixt_comp_k_lh_start, &
               X_u_s_mellor_k_lh_start )
! Description:
!   Verify that half the points are in cloud if cloud weighted sampling is
!   enabled and other conditions are met. We stop rather than exiting gracefully
!   if this assertion check fails.

! References:
!   None.
!-------------------------------------------------------------------------------

    use mt95, only: &
      genrand_real ! Constant

    use constants_clubb, only: &
      fstderr ! Constant

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      n_micro_calls ! Total calls to the microphysics

    real, intent(in) :: &
      cloud_frac1, cloud_frac2 ! Cloud fraction associatated with component 1 & 2 [-]

    integer, dimension(n_micro_calls), intent(in) :: &
      X_mixt_comp_k_lh_start ! Mixture components at k_lh_start

    real(kind=genrand_real), dimension(n_micro_calls), intent(in) :: &
      X_u_s_mellor_k_lh_start   ! Uniform distribution for s_mellor at k_lh_start [-]

    ! Local Variables
    real :: cloud_frac_n

    integer :: number_cloudy_samples

    integer :: sample  ! Loop iterator

    ! ---- Begin Code ----

    number_cloudy_samples = 0

    do sample = 1, n_micro_calls
      if ( X_mixt_comp_k_lh_start(sample) == 1 ) then
        cloud_frac_n = cloud_frac1
      else
        cloud_frac_n = cloud_frac2
      end if
      if ( X_u_s_mellor_k_lh_start(sample) >= 1.-cloud_frac_n ) then
        number_cloudy_samples = number_cloudy_samples + 1
      else
        ! Do nothing, the air is clear
      end if
    end do
    if ( number_cloudy_samples /= ( n_micro_calls / 2 ) ) then
      write(fstderr,*) "Error, half of all samples aren't in cloud"
      write(fstderr,*) "X_u s_mellor random = ", &
        X_u_s_mellor_k_lh_start(:), "X_mixt_comp = ", X_mixt_comp_k_lh_start, &
        "cloudy samples =", number_cloudy_samples
      write(fstderr,*) "cloud_frac1 = ", cloud_frac1
      write(fstderr,*) "cloud_frac2 = ", cloud_frac2
      stop "Fatal Error"
    end if

    return
  end subroutine assert_check_half_cloudy

!-------------------------------------------------------------------------------
  function compute_vert_corr( nnzp, delta_zm, Lscale_vert_avg ) result( vert_corr )
! Description:
!   This function computes the vertical correlation for arbitrary overlap, using
!   density weighted 3pt averaged Lscale and the difference in height levels
!   (delta_zm).
! References:
!   None
!-------------------------------------------------------------------------------

    implicit none

    ! External
    intrinsic :: exp

    ! Parameter Constants
    real, parameter :: &
      vert_corr_coef = 0.03 ! Empirically defined correlation constant [-]

    ! Input Variables
    integer, intent(in) :: &
      nnzp ! Number of vertical levels  [-]

    real, intent(in), dimension(nnzp) :: &
      delta_zm, &     ! Difference between altitudes    [m]
      Lscale_vert_avg ! Vertically averaged Lscale      [m]

    ! Output Variable
    real, dimension(nnzp) :: &
      vert_corr ! The vertical correlation      [-]

    ! ---- Begin Code ----
    vert_corr(1:nnzp) = exp( -vert_corr_coef * ( delta_zm(1:nnzp) / Lscale_vert_avg(1:nnzp) ) )

    return
  end function compute_vert_corr

!-------------------------------------------------------------------------------
  subroutine LH_clipping_and_stats( nnzp, n_micro_calls, d_variables, &
                                    LH_sample_point_weights, X_nl_all_levs, LH_thl, &
                                    LH_rt )
! Description:
!   Clip subcolumns from latin hypercube and create stats for diagnostic
!   purposes.

! References:
!   None
!-------------------------------------------------------------------------------

    use parameters_model, only: hydromet_dim ! Variable

    use stats_variables, only: &
      l_stats_samp, & ! Variable(s)
      iLH_rrainm, &
      iLH_Nrm, &
      iLH_ricem, &
      iLH_Nim, &
      iLH_rsnowm, &
      iLH_Nsnowm, &
      iLH_rgraupelm, &
      iLH_Ngraupelm, &
      iLH_thlm, &
      iLH_rcm, &
      iLH_Ncm, &
      iLH_rvm, &
      iLH_wm, &
      iLH_cloud_frac

    use stats_variables, only: &
      iLH_wp2_zt, &  ! Variable(s)
      iLH_Nrp2_zt, &
      iLH_Ncp2_zt, &
      iLH_rcp2_zt, &
      iLH_rtp2_zt, &
      iLH_thlp2_zt, &
      iLH_rrainp2_zt, &
      zt

    use math_utilities, only: &
      compute_sample_mean, & ! Procedure(s)
      compute_sample_variance

    use stats_type, only: &
      stat_update_var ! Procedure(s)

    use array_index, only: &
      iirrainm, & ! Variables
      iirsnowm, & 
      iiricem, & 
      iirgraupelm, & 
      iiNrm, &
      iiNsnowm, &
      iiNim, &
      iiNgraupelm, &
      iiNcm

    use latin_hypercube_arrays, only: &
      iiLH_rrain, & ! Variable(s)
      iiLH_rsnow, &
      iiLH_rice, &
      iiLH_rgraupel, &
      iiLH_Nr, &
      iiLH_Nsnow, &
      iiLH_Ngraupel, &
      iiLH_Nc, &
      iiLH_Ni, &
      iiLH_s_mellor, &
      iiLH_w

    use error_code, only: &
      clubb_at_least_debug_level ! Procedure(s)

    use estimate_scm_microphys_module, only: &
      copy_X_nl_into_hydromet ! Procedure(s)

    use constants_clubb, only: &
      zero_threshold, & ! Constant(s)
      fstderr

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      d_variables,     & ! Number of variables to sample
      n_micro_calls,   & ! Number of calls to microphysics per timestep (normally=2)
      nnzp               ! Number of vertical model levels

    real, intent(in), dimension(n_micro_calls) :: &
      LH_sample_point_weights

    double precision, intent(in), dimension(nnzp,n_micro_calls,d_variables) :: &
       X_nl_all_levs ! Sample that is transformed ultimately to normal-lognormal

    real, intent(in), dimension(nnzp,n_micro_calls) :: &
      LH_thl ! Sample of liquid potential temperature [K]

    ! Input/Output Variables
    real, intent(inout), dimension(nnzp,n_micro_calls) :: &
      LH_rt ! Sample of total water mixing ratio    [kg/kg]

    ! Local variables
    real, dimension(nnzp,n_micro_calls) :: &
      rc_all_points, & ! Cloud water mixing ratio for all levels   [kg/kg]
      rv_all_points    ! Vapor mixing ratio for all levels   [kg/kg]

    real, dimension(nnzp,n_micro_calls,hydromet_dim) :: &
      hydromet_all_points ! Hydrometeor species    [units vary]

    real, dimension(nnzp,hydromet_dim) :: &
      LH_hydromet ! Average value of the latin hypercube est. of all hydrometeors [units vary]

    real, dimension(nnzp) :: &
      LH_thlm,       & ! Average value of the latin hypercube est. of theta_l           [K]
      LH_rcm,        & ! Average value of the latin hypercube est. of rc                [kg/kg]
      LH_rvm,        & ! Average value of the latin hypercube est. of rv                [kg/kg]
      LH_wm,         & ! Average value of the latin hypercube est. of vertical velocity [m/s]
      LH_wp2_zt,     & ! Average value of the variance of the LH est. of vert. vel.     [m^2/s^2]
      LH_rrainp2_zt, & ! Average value of the variance of the LH est. of rrain.         [(kg/kg)^2]
      LH_rcp2_zt,    & ! Average value of the variance of the LH est. of rc.            [(kg/kg)^2]
      LH_rtp2_zt,    & ! Average value of the variance of the LH est. of rt             [kg^2/kg^2]
      LH_thlp2_zt,   & ! Average value of the variance of the LH est. of thetal         [K^2]
      LH_Nrp2_zt,    & ! Average value of the variance of the LH est. of Nr.            [#^2/kg^2]
      LH_Ncp2_zt,    & ! Average value of the variance of the LH est. of Nc.            [#^2/kg^2]
      LH_cloud_frac    ! Average value of the latin hypercube est. of cloud fraction    [-]

    integer :: k, sample, ivar

    ! ---- Begin Code ----

    ! Clip 's' from Mellor to obtain cloud-water mixing ratio
    rc_all_points = max( zero_threshold, X_nl_all_levs(:,:,iiLH_s_mellor) )

    ! Verify total water isn't negative
    if ( any( LH_rt < 0. ) ) then
      if ( clubb_at_least_debug_level( 1 ) ) then
        write(fstderr,*) "Total water negative in LH sample"
        write(fstderr,*) "Applying non-conservative hard clipping to rv sample."
      end if ! clubb_at_least_debug_level( 1 )
      where ( LH_rt < 0. )
        LH_rt = zero_threshold
      end where
    end if ! Some rv_all_points(:,sample) < 0

    if ( l_stats_samp ) then

      ! For all cases where l_lh_cloud_weighted_sampling is false, the weights
      ! will be 1 (all points equally weighted)

      if ( iLH_rcm + iLH_rcp2_zt > 0 ) then
        LH_rcm = compute_sample_mean( nnzp, n_micro_calls, LH_sample_point_weights, &
                                      rc_all_points )
      end if

      if ( iLH_thlm + iLH_thlp2_zt > 0 ) then
        LH_thlm = compute_sample_mean( nnzp, n_micro_calls, LH_sample_point_weights, &
                                       real( LH_thl ) )
      end if

      if ( iLH_rvm + iLH_rtp2_zt > 0 ) then
        rv_all_points = LH_rt - rc_all_points
        LH_rvm = compute_sample_mean( nnzp, n_micro_calls, LH_sample_point_weights, &
                                      rv_all_points )
      end if

      if ( iLH_wm + iLH_wp2_zt > 0 ) then
        LH_wm  = compute_sample_mean( nnzp, n_micro_calls, LH_sample_point_weights, &
                                      real( X_nl_all_levs(:,:,iiLH_w) ) )
      end if

      if ( iLH_rrainm + iLH_Nrm + iLH_ricem + iLH_Nim + iLH_rsnowm + iLH_Nsnowm + &
           iLH_rgraupelm + iLH_Ngraupelm + iLH_Ncm > 0 ) then

        LH_hydromet = 0.
        call copy_X_nl_into_hydromet( nnzp, d_variables, n_micro_calls, & ! In
                                      X_nl_all_levs, &  ! In
                                      LH_hydromet, & ! In
                                      hydromet_all_points ) ! Out

        forall ( ivar = 1:hydromet_dim )
          LH_hydromet(:,ivar) = compute_sample_mean( nnzp, n_micro_calls, LH_sample_point_weights, &
                                                     hydromet_all_points(:,:,ivar) )
        end forall ! 1..hydromet_dim
      end if

      ! Latin hypercube estimate of cloud fraction
      if ( iLH_cloud_frac > 0 ) then
        LH_cloud_frac = 0.
        do sample = 1, n_micro_calls
          where ( X_nl_all_levs(:,sample,iiLH_s_mellor) > 0. )
            LH_cloud_frac(:) = LH_cloud_frac(:) + 1.0 * LH_sample_point_weights(sample)
          end where
        end do
        LH_cloud_frac(:) = LH_cloud_frac(:) / real( n_micro_calls )
      end if

      if ( iLH_wp2_zt > 0 ) then
        ! Compute the variance of vertical velocity
        LH_wp2_zt = compute_sample_variance( nnzp, n_micro_calls, &
                                             real( X_nl_all_levs(:,:,iiLH_w) ), &
                                             LH_sample_point_weights, LH_wm )
      end if

      if ( iLH_rcp2_zt  > 0 ) then
        ! Compute the variance of cloud water mixing ratio
        LH_rcp2_zt = compute_sample_variance &
                     ( nnzp, n_micro_calls, rc_all_points, &
                       LH_sample_point_weights, LH_rcm )
      end if

      if ( iLH_rtp2_zt > 0 ) then
        ! Compute the variance of total water
        LH_rtp2_zt = compute_sample_variance &
                     ( nnzp, n_micro_calls, &
                       real( LH_rt ), LH_sample_point_weights, LH_rvm+LH_rcm )
      end if

      if ( iLH_thlp2_zt > 0 ) then
        ! Compute the variance of liquid potential temperature
        LH_thlp2_zt = compute_sample_variance( nnzp, n_micro_calls, &
                                               real( LH_thl ), LH_sample_point_weights, LH_thlm )
      end if

      ! Compute the variance of rain water mixing ratio
      if ( iirrainm > 0 .and. iLH_rrainp2_zt > 0 ) then
        LH_rrainp2_zt = compute_sample_variance &
                        ( nnzp, n_micro_calls, hydromet_all_points(:,:,iirrainm), &
                          LH_sample_point_weights, LH_hydromet(:,iirrainm) )
      end if

      ! Compute the variance of cloud droplet number concentration
      if ( iiNcm > 0 .and. iLH_Ncp2_zt > 0 ) then
        LH_Ncp2_zt = compute_sample_variance &
                     ( nnzp, n_micro_calls, hydromet_all_points(:,:,iiNcm), &
                       LH_sample_point_weights, LH_hydromet(:,iiNcm) )
      end if

      ! Compute the variance of rain droplet number concentration
      if ( iiNrm > 0 .and. iLH_Nrp2_zt > 0 ) then
        LH_Nrp2_zt = compute_sample_variance( nnzp, n_micro_calls, hydromet_all_points(:,:,iiNrm), &
                                              LH_sample_point_weights, LH_hydromet(:,iiNrm) )
      end if

      ! Averages of points being fed into the microphysics
      ! These are for diagnostic purposes, and are not needed for anything
      if ( iirrainm > 0 ) then
        call stat_update_var( iLH_rrainm, LH_hydromet(:,iirrainm), zt )
      end if
      if ( iiNrm > 0 ) then
        call stat_update_var( iLH_Nrm, LH_hydromet(:,iiNrm), zt )
      end if
      if ( iiricem > 0 ) then
        call stat_update_var( iLH_ricem, LH_hydromet(:,iiricem), zt )
      end if
      if ( iiNim > 0 ) then
        call stat_update_var( iLH_Nim, LH_hydromet(:,iiNim), zt )
      end if
      if ( iirsnowm > 0 ) then
        call stat_update_var( iLH_rsnowm, LH_hydromet(:,iirsnowm), zt )
      end if
      if ( iiNsnowm > 0 ) then
        call stat_update_var( iLH_Nsnowm, LH_hydromet(:,iiNsnowm), zt )
      end if
      if ( iirgraupelm > 0 ) then
        call stat_update_var( iLH_rgraupelm, LH_hydromet(:,iirgraupelm), zt )
      end if
      if ( iiNgraupelm > 0 ) then
        call stat_update_var( iLH_Ngraupelm, LH_hydromet(:,iiNgraupelm), zt )
      end if
      if ( iiNcm > 0 ) then
        call stat_update_var( iLH_Ncm, LH_hydromet(:,iiNcm), zt )
      end if

      call stat_update_var( iLH_rcm, LH_rcm, zt )
      call stat_update_var( iLH_thlm, LH_thlm, zt )
      call stat_update_var( iLH_rvm, LH_rvm, zt )
      call stat_update_var( iLH_wm, LH_wm, zt )
      call stat_update_var( iLH_wp2_zt, LH_wp2_zt, zt )
      call stat_update_var( iLH_Ncp2_zt, LH_Ncp2_zt, zt )
      call stat_update_var( iLH_Nrp2_zt, LH_Nrp2_zt, zt )
      call stat_update_var( iLH_rcp2_zt, LH_rcp2_zt, zt )
      call stat_update_var( iLH_rtp2_zt, LH_rtp2_zt, zt )
      call stat_update_var( iLH_thlp2_zt, LH_thlp2_zt, zt )
      call stat_update_var( iLH_rrainp2_zt, LH_rrainp2_zt, zt )
      call stat_update_var( iLH_cloud_frac, LH_cloud_frac, zt )

    end if ! l_stats_samp

    return
  end subroutine LH_clipping_and_stats

#endif /*LATIN_HYPERCUBE*/

end module latin_hypercube_driver_module
