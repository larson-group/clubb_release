! $Id$
!-------------------------------------------------------------------------------
module latin_hypercube_driver_module

  implicit none

  public :: latin_hypercube_driver, latin_hypercube_2D_output, &
    latin_hypercube_2D_close

  private ! Default scope

  ! Constant Parameters
  logical, parameter, private :: &
    l_diagnostic_iter_check      = .true., & ! Check for a problem in iteration
    l_output_2D_lognormal_dist   = .false., & ! Output a 2D netCDF file of the lognormal variates
    l_output_2D_uniform_dist     = .false.    ! Output a 2D netCDF file of the uniform distribution

  integer, allocatable, dimension(:,:,:), private :: & 
    height_time_matrix ! matrix of rand ints

  integer, private :: &
    prior_iter ! Prior iteration number (for diagnostic purposes)

  contains

!-------------------------------------------------------------------------------
  subroutine latin_hypercube_driver &
             ( dt, iter, d_variables, n_micro_calls, sequence_length, nnzp, &
               cloud_frac, thlm, p_in_Pa, exner, &
               rho, pdf_params, wm, w_std_dev, dzq, rcm, rvm, &
               hydromet, correlation_array, LH_hydromet_mc, LH_hydromet_vel, LH_rcm_mc, &
               LH_rvm_mc, LH_thlm_mc, microphys_sub )

! Description:
!   Call a microphysics scheme or generate a estimate of Kessler autoconversion
!   using latin hypercube sampling.
! References:
!   None
!-------------------------------------------------------------------------------
    use array_index, only: & 
      iirrainm,    & ! Variables  
      iirsnowm,    &
      iirgraupelm, &
      iiricem,     &
      iiNcm,       & 
      iiNrm,       & 
      iiNim,       & 
      iiNsnowm,    & 
      iiNgraupelm

    use array_index, only: & 
      iiLH_rt   ! Variables

    use parameters_model, only: hydromet_dim ! Variable

#ifdef UNRELEASED_CODE
    use permute_height_time_module, only: & 
      permute_height_time ! Procedure

    use generate_lh_sample_module, only: & 
      generate_lh_sample, & ! Procedure
      generate_uniform_sample

    use estimate_lh_micro_module, only: & 
      estimate_lh_micro, & ! Procedure
      k_lh_start ! Variable

    use output_2D_samples_module, only: &
      output_2D_lognormal_dist_file, & ! Procedure(s)
      output_2D_uniform_dist_file 
#endif

    use variables_prognostic_module, only: &
      pdf_parameter  ! Type

    use constants, only: &
      fstderr, & ! Constant
      cm3_per_m3

    use variables_diagnostic_module, only: & 
      lh_AKm,  & 
      AKm, & 
      AKstd, & 
      AKstd_cld, & 
      AKm_rcm, & 
      AKm_rcc, & 
      lh_rcm_avg

    use stats_variables, only: &
      l_stats_samp, & ! Variables
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
      iLH_wp2_zt, &
      iLH_Nrp2_zt, &
      iLH_Ncp2_zt, &
      iLH_rrainp2_zt, &
      iLH_rcp2_zt, &
      iLH_rtp2_zt, &
      iLH_thlp2_zt, &
      iLH_cloud_frac, &
      zt

    use stats_type, only: &
      stat_update_var ! Procedure(s)

    use parameters_microphys, only: &
      l_lh_vert_overlap, &  ! Variables
      LH_sample_point_weights, &
      l_lh_cloud_weighted_sampling

    use error_code, only: &
      clubb_at_least_debug_level ! Procedure

    use mt95, only: genrand_real ! Constants

    implicit none

    ! External
    intrinsic :: allocated, mod, maxloc, dble, epsilon

    ! Interface block
#include "./microphys_interface.inc"

    ! Parameter Constants
    real, parameter :: &
      cloud_frac_thresh = 0.01 ! Threshold for sampling preferentially within cloud

    ! Find in and out of cloud points using the rejection method rather than scaling
    logical, parameter :: &
      l_use_rejection_method = .true.

    ! Input Variables
    real, intent(in) :: &
      dt ! Model timestep       [s]

    integer, intent(in) :: &
      iter,            & ! Model iteration number
      d_variables,     & ! Number of variables to sample
      n_micro_calls,   & ! Number of calls to microphysics per timestep (normally=2)
      sequence_length, & ! nt_repeat/n_micro_call; number of timesteps before sequence repeats.
      nnzp               ! Number of vertical model levels

    real, dimension(nnzp), intent(in) :: &
      cloud_frac, & ! Cloud fraction               [-]
      thlm,       & ! Liquid potential temperature [K]
      p_in_Pa,    & ! Pressure                     [Pa]
      exner,      & ! Exner function               [-]
      rho           ! Density on thermo. grid      [kg/m^3]

    real, dimension(nnzp), intent(in) :: &
      wm, &        ! Mean w                     [m/s]
      w_std_dev, & ! Standard deviation of w    [m/s]
      dzq          ! Difference in altitudes    [m]

    real, dimension(nnzp), intent(in) :: &
      rcm, & ! Liquid water mixing ratio        [kg/kg]
      rvm    ! Vapor water mixing ratio         [kg/kg]

    real, dimension(nnzp,hydromet_dim), intent(in) :: &
      hydromet ! Hydrometeor species    [units vary]

    real, dimension(nnzp,d_variables,d_variables), intent(in) :: &
      correlation_array ! Correlation for hydrometeor species [-]

    ! Input/Output Variables
    real, dimension(nnzp,hydromet_dim), intent(inout) :: &
      LH_hydromet_mc, & ! LH estimate of hydrometeor time tendency          [(units vary)/s]
      LH_hydromet_vel   ! LH estimate of hydrometeor sedimentation velocity [m/s]

    ! Output Variables
    real, dimension(nnzp), intent(out) :: &
      LH_rcm_mc, & ! LH estimate of time tendency of liquid water mixing ratio    [kg/kg/s]
      LH_rvm_mc, & ! LH estimate of time tendency of vapor water mixing ratio     [kg/kg/s]
      LH_thlm_mc   ! LH estimate of time tendency of liquid potential temperature [K/s]

    type(pdf_parameter), intent(in) :: pdf_params

    ! Local variables
    integer :: p_matrix(n_micro_calls,d_variables+1)

    real(kind=genrand_real), dimension(nnzp,n_micro_calls,(d_variables+1)) :: &
      X_u_all_levs ! Sample drawn from uniform distribution

    double precision, dimension(nnzp,n_micro_calls,d_variables) :: &
      X_nl_all_levs ! Sample that is transformed ultimately to normal-lognormal

    integer, dimension(nnzp,n_micro_calls) :: &
      X_mixt_comp_all_levs ! Which mixture component we're in

    real, dimension(nnzp,n_micro_calls) :: &
      LH_rt, LH_thl ! Sample of total water and liquid potential temperature [kg/kg],[K]

    real, dimension(nnzp,hydromet_dim) :: &
      lh_hydromet ! Average value of the latin hypercube est. of all hydrometeors [units vary]

    real, dimension(nnzp) :: &
      lh_thlm,       & ! Average value of the latin hypercube est. of theta_l           [K]
      lh_rcm,        & ! Average value of the latin hypercube est. of rc                [kg/kg]
      lh_rvm,        & ! Average value of the latin hypercube est. of rv                [kg/kg]
      lh_wm,         & ! Average value of the latin hypercube est. of vertical velocity [m/s]
      lh_wp2_zt,     & ! Average value of the variance of the LH est. of vert. vel.     [m^2/s^2]
      lh_rrainp2_zt, & ! Average value of the variance of the LH est. of rrain.         [(kg/kg)^2]
      lh_rcp2_zt,    & ! Average value of the variance of the LH est. of rc.            [(kg/kg)^2]
      lh_rtp2_zt,    & ! Average value of the variance of the LH est. of rt             [kg^2/kg^2]
      lh_thlp2_zt,   & ! Average value of the variance of the LH est. of thetal         [K^2]
      lh_Nrp2_zt,    & ! Average value of the variance of the LH est. of Nr.            [#^2/kg^2]
      lh_Ncp2_zt,    & ! Average value of the variance of the LH est. of Nc.            [#^2/kg^2]
      lh_cloud_frac    ! Average value of the latin hypercube est. of cloud fraction    [-]

    real :: lh_start_cloud_frac ! Cloud fraction at k_lh_start [-]

    ! Try to obtain 12 digit accuracy for a diagnostic mean
    real(kind=selected_real_kind( p=12 ) ) :: mean_weight

    double precision :: mixt_frac_dp

    real(kind=genrand_real) :: X_u_dp1_element, X_u_s_mellor_element

    ! Number of random samples before sequence of repeats (normally=10)
    integer :: nt_repeat

    real :: cloud_frac_n ! Temporary variable for cloud_frac1/cloud_frac2

    integer :: &
      i_rmd, &   ! Remainder of ( iter-1 / sequence_length )
      k, sample  ! Loop iterators

    integer :: &
      iiLH_s_mellor, &      ! Index of Mellor's s (extended rcm)
      number_cloudy_samples ! Diagnostic sum

    logical :: l_cloudy_sample ! Whether a sample point is cloudy or clear air

    integer, dimension(1) :: tmp_loc

#ifdef UNRELEASED_CODE
    ! ---- Begin Code ----

    nt_repeat = n_micro_calls * sequence_length

    if ( .not. allocated( height_time_matrix ) ) then
      ! If this is first time latin_hypercube_driver is called, then allocate
      ! the height_time_matrix and set the prior iteration number for debugging
      ! purposes.
      allocate( height_time_matrix(nnzp, nt_repeat, d_variables+1) )

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
    if ( l_lh_cloud_weighted_sampling .and. mod( n_micro_calls, 2 ) /= 0 ) then
      write(fstderr,*) "Cloud weighted sampling requires micro calls to be divisible by 2."
      stop "Fatal error."
    end if
    if ( l_lh_cloud_weighted_sampling .and. sequence_length /= 1 ) then
      write(fstderr,*) "Cloud weighted sampling requires sequence length be equal to 1."
      stop "Fatal error."
    end if

    ! Determine the s_mellor element of the uniform/lognormal distribution arrays
    iiLH_s_mellor = iiLH_rt

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

    do k = 1, nnzp
      ! Choose which rows of LH sample to feed into closure.
      p_matrix(1:n_micro_calls,1:(d_variables+1)) = &
        height_time_matrix(k, n_micro_calls*i_rmd+1:n_micro_calls*i_rmd+n_micro_calls, &
                           1:d_variables+1)

      ! print*, 'latin_hypercube_sampling: got past p_matrix'

      ! Generate the uniform distribution using the Mersenne twister (not
      ! currently configured to re-seed).
      !  X_u has one extra dimension for the mixture component.
      call generate_uniform_sample( n_micro_calls, nt_repeat, d_variables+1, p_matrix, & ! In
                                    X_u_all_levs(k,:,:) ) ! Out

    end do ! 1..nnzp

    ! For a 100 level fixed grid, this looks to be about the middle of the cloud for RICO
    k_lh_start = 50
!   tmp_loc    = maxloc( rcm )
!   k_lh_start = tmp_loc(1) ! Attempt using the maximal value of rcm for now

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

    do k = 1, nnzp
      ! Choose which rows of LH sample to feed into closure.
      p_matrix(1:n_micro_calls,1:(d_variables+1)) = &
        height_time_matrix(k, n_micro_calls*i_rmd+1:n_micro_calls*i_rmd+n_micro_calls, &
                           1:d_variables+1)

      ! print*, 'latin_hypercube_sampling: got past p_matrix'

      ! Generate the uniform distribution using the Mersenne twister (not
      ! currently configured to re-seed).
      !  X_u has one extra dimension for the mixture component.
      call generate_uniform_sample( n_micro_calls, nt_repeat, d_variables+1, p_matrix, & ! In
                                    X_u_all_levs(k,:,:) ) ! Out

    end do ! 1..nnzp


    ! For a 100 level fixed grid, this looks to be about the middle of the cloud for RICO
!   k_lh_start = 50
    tmp_loc    = maxloc( rcm )
    k_lh_start = tmp_loc(1) ! Attempt using the maximal value of rcm for now

    if ( l_lh_cloud_weighted_sampling ) then

      ! Determine cloud fraction at k_lh_start
      lh_start_cloud_frac = &
        pdf_params%mixt_frac(k_lh_start) * pdf_params%cloud_frac1(k_lh_start) &
        + (1.0-pdf_params%mixt_frac(k_lh_start)) * pdf_params%cloud_frac2(k_lh_start)

      ! Determine p_matrix at k_lh_start
      p_matrix(1:n_micro_calls,1:(d_variables+1)) = &
        height_time_matrix(k_lh_start, n_micro_calls*i_rmd+1:n_micro_calls*i_rmd+n_micro_calls, &
                           1:d_variables+1)

    end if ! l_lh_cloud_weighted_sampling

    if ( l_lh_vert_overlap ) then

      do sample = 1, n_micro_calls

        ! Get the uniform sample of d+1 (telling us if we're in component 1 or
        ! component 2), and the uniform sample for s_mellor
        X_u_dp1_element      = X_u_all_levs(k_lh_start,sample,d_variables+1)
        X_u_s_mellor_element = X_u_all_levs(k_lh_start,sample,iiLH_s_mellor)

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
                   ( l_cloudy_sample, pdf_params%cloud_frac1(k_lh_start), & ! In
                     pdf_params%cloud_frac2(k_lh_start), pdf_params%mixt_frac(k_lh_start), & !In
                     cloud_frac_thresh, & ! In
                     X_u_dp1_element, X_u_s_mellor_element ) ! In/out

            else ! Transpose and scale the points to be in or out of cloud
              call choose_X_u_scaled &
                   ( l_cloudy_sample, pdf_params%cloud_frac1(k_lh_start), & ! In
                     pdf_params%cloud_frac2(k_lh_start), pdf_params%mixt_frac(k_lh_start), & !In
                     X_u_dp1_element, X_u_s_mellor_element ) ! In/out

            end if

          end if ! Cloud fraction is between cloud_frac_thresh and 50%

        end if ! l_lh_cloud_weighted_sampling

        do k = 1, nnzp

          ! Overwrite 2 elements for maximal overlap and cloud
          ! weighted sampling (if it's enabled).
          X_u_all_levs(k,sample,iiLH_s_mellor) =  X_u_s_mellor_element ! s_mellor
          X_u_all_levs(k,sample,d_variables+1) = X_u_dp1_element ! Mixture comp.

          ! Use this line to have all variates maximally correlated.  This
          ! probably isn't a good idea for dealing with all cloud types.
!         X_u_all_levs(k,sample,:) = X_u_all_levs(k_lh_start,sample,:)


        end do ! 1..nnzp

      end do ! 1..n_micro_calls

    end if ! l_lh_vert_overlap

    ! Determine mixture component for all levels
    do k = 1, nnzp

      mixt_frac_dp = dble( pdf_params%mixt_frac(k) )

      where ( in_mixt_frac_1(X_u_all_levs(k,:,d_variables+1), mixt_frac_dp ) )
        X_mixt_comp_all_levs(k,:) = 1
      else where
        X_mixt_comp_all_levs(k,:) = 2
      end where

    end do ! k = 1 .. nnzp

    ! Assertion check for whether half of sample points are cloudy.
    ! This is for the uniform sample only.  Another assertion check is in the
    ! estimate_lh_micro_module for X_nl_all_levs.
    ! TODO: put this in its own subroutine
    if ( l_lh_cloud_weighted_sampling .and. clubb_at_least_debug_level( 2 ) .and. &
         lh_start_cloud_frac < 0.5 .and. lh_start_cloud_frac > cloud_frac_thresh ) then

      number_cloudy_samples = 0

      do sample = 1, n_micro_calls
        if ( X_mixt_comp_all_levs(k_lh_start,sample) == 1 ) then
          cloud_frac_n = pdf_params%cloud_frac1(k_lh_start)
        else
          cloud_frac_n = pdf_params%cloud_frac2(k_lh_start)
        end if
        if ( X_u_all_levs(k_lh_start,sample,iiLH_s_mellor) >= 1.-cloud_frac_n ) then
          number_cloudy_samples = number_cloudy_samples + 1
        else
          ! Do nothing, the air is clear
        end if
      end do
      if ( number_cloudy_samples /= ( n_micro_calls / 2 ) ) then
        write(fstderr,*) "Error, half of all samples aren't in cloud"
        write(fstderr,*) "X_u s_smellor random = ", &
          X_u_all_levs(k_lh_start,:,iiLH_s_mellor), "cloudy samples =", number_cloudy_samples
        write(fstderr,*) "cloud_frac1 = ", pdf_params%cloud_frac1(k_lh_start)
        write(fstderr,*) "cloud_frac2 = ", pdf_params%cloud_frac2(k_lh_start)
        write(fstderr,*) "X_u d+1 element = ", X_u_all_levs(k_lh_start,:,d_variables+1)
        write(fstderr,*) "mixture fraction = ", mixt_frac_dp
        stop "Fatal Error"
      end if

    end if ! Maximal overlap, debug_level 2, and cloud-weighted averaging

    ! Assertion check to ensure that the sample point weights sum to approximately 1
    if ( l_lh_cloud_weighted_sampling .and. clubb_at_least_debug_level( 2 ) ) then
      mean_weight = 0.
      do sample = 1, n_micro_calls
        mean_weight = mean_weight + LH_sample_point_weights(sample)
      end do
      mean_weight = mean_weight / real( n_micro_calls )

      ! Using more precision for mean_weight should make this work out
      if ( abs( mean_weight - 1.0 ) > epsilon( LH_sample_point_weights ) ) then
        write(fstderr,*) "Error in cloud weighted sampling code ", "mean_weight = ", mean_weight
        stop
      end if

    end if ! l_lh_cloud_weighted_sampling .and. clubb_at_least_debug_level( 2 )

    ! Upwards loop
    do k = k_lh_start, nnzp, 1
      ! Generate LH sample, represented by X_u and X_nl, for level k
      call generate_lh_sample &
           ( n_micro_calls, d_variables, hydromet_dim, &        ! In
             cloud_frac(k), wm(k), rcm(k)+rvm(k), thlm(k), pdf_params, k, & ! In
             hydromet(k,:), correlation_array(k,:,:), X_u_all_levs(k,:,:), & ! In
             X_mixt_comp_all_levs(k,:), & ! In
             LH_rt(k,:), LH_thl(k,:), X_nl_all_levs(k,:,:) ) ! Out
    end do ! k = k_lh_start..nnzp

    ! Downwards loop
    do k = k_lh_start-1, 1, -1
      call generate_lh_sample &
           ( n_micro_calls, d_variables, hydromet_dim, &        ! In
             cloud_frac(k), wm(k), rcm(k)+rvm(k), thlm(k), pdf_params, k, & ! In
             hydromet(k,:), correlation_array(k,:,:), X_u_all_levs(k,:,:), &  !  In
             X_mixt_comp_all_levs(k,:), & ! In
             LH_rt(k,:), LH_thl(k,:), X_nl_all_levs(k,:,:) ) ! Out
    end do ! k_lh_start-1..1

    if ( l_output_2D_lognormal_dist ) then
      call output_2D_lognormal_dist_file( nnzp, n_micro_calls, d_variables, &
                                          X_nl_all_levs, LH_rt, LH_thl )
    end if
    if ( l_output_2D_uniform_dist ) then
      call output_2D_uniform_dist_file( nnzp, n_micro_calls, d_variables+1, &
                                        X_u_all_levs )
    end if
    ! Perform LH and analytic microphysical calculations
    call estimate_lh_micro &
         ( dt, nnzp, n_micro_calls, d_variables, &  ! intent(in)
           X_u_all_levs, X_nl_all_levs, &           ! intent(in)
           LH_rt, LH_thl, pdf_params, &             ! intent(in)
           p_in_Pa, exner, rho, &                   ! intent(in)
           rcm, w_std_dev, dzq, &                   ! intent(in)
           cloud_frac, hydromet, &                  ! intent(in)
           X_mixt_comp_all_levs, &                  ! intent(in)
           LH_hydromet_mc, LH_hydromet_vel, &       ! intent(in)
           LH_rcm_mc, LH_rvm_mc, LH_thlm_mc, &      ! intent(out)
           lh_AKm, AKm, AKstd, AKstd_cld, &         ! intent(out)
           AKm_rcm, AKm_rcc, lh_rcm_avg, &          ! intent(out)
           lh_hydromet, lh_thlm, lh_rcm, lh_rvm, &  ! intent(out)
           lh_wm, lh_Ncp2_zt, lh_Nrp2_zt, lh_rrainp2_zt, lh_rcp2_zt, &  ! intent(out)
           lh_wp2_zt, lh_rtp2_zt, lh_thlp2_zt, & ! intent(out)
           lh_cloud_frac, & ! intent(out)
           microphys_sub )  ! Procedure

    ! print*, 'latin_hypercube_driver: AKm=', AKm
    ! print*, 'latin_hypercube_driver: lh_AKm=', lh_AKm

    if ( l_stats_samp ) then

      ! Averages of points being fed into the microphysics
      ! These are for diagnostic purposes, and are not needed for anything
      if ( iirrainm > 0 ) then
        call stat_update_var( iLH_rrainm, lh_hydromet(:,iirrainm), zt )
      end if
      if ( iiNrm > 0 ) then
        call stat_update_var( iLH_Nrm, lh_hydromet(:,iiNrm), zt )
      end if
      if ( iiricem > 0 ) then
        call stat_update_var( iLH_ricem, lh_hydromet(:,iiricem), zt )
      end if
      if ( iiNim > 0 ) then
        call stat_update_var( iLH_Nim, lh_hydromet(:,iiNim), zt )
      end if
      if ( iirsnowm > 0 ) then
        call stat_update_var( iLH_rsnowm, lh_hydromet(:,iirsnowm), zt )
      end if
      if ( iiNsnowm > 0 ) then
        call stat_update_var( iLH_Nsnowm, lh_hydromet(:,iiNsnowm), zt )
      end if
      if ( iirgraupelm > 0 ) then
        call stat_update_var( iLH_rgraupelm, lh_hydromet(:,iirgraupelm), zt )
      end if
      if ( iiNgraupelm > 0 ) then
        call stat_update_var( iLH_Ngraupelm, lh_hydromet(:,iiNgraupelm), zt )
      end if

      call stat_update_var( iLH_rcm, lh_rcm, zt )
      if ( iiNcm > 0 ) then
        call stat_update_var( iLH_Ncm, lh_hydromet(:,iiNcm), zt )
      end if
      call stat_update_var( iLH_thlm, lh_thlm, zt )
      call stat_update_var( iLH_rvm, lh_rvm, zt )
      call stat_update_var( iLH_wm, lh_wm, zt )
      call stat_update_var( iLH_wp2_zt, lh_wp2_zt, zt )
      call stat_update_var( iLH_Ncp2_zt, lh_Ncp2_zt, zt )
      call stat_update_var( iLH_Nrp2_zt, lh_Nrp2_zt, zt )
      call stat_update_var( iLH_rcp2_zt, lh_rcp2_zt, zt )
      call stat_update_var( iLH_rtp2_zt, lh_rtp2_zt, zt )
      call stat_update_var( iLH_thlp2_zt, lh_thlp2_zt, zt )
      call stat_update_var( iLH_rrainp2_zt, lh_rrainp2_zt, zt )
      call stat_update_var( iLH_cloud_frac, lh_cloud_frac, zt )

    end if ! l_stats_samp

    return

#else
    stop "This code was not compiled with support for Latin Hypercube sampling"

    ! This is simply to avoid a compiler warning
    LH_rcm_mc  = -999.999
    LH_rvm_mc  = -999.999
    LH_thlm_mc = -999.999

#endif /* UNRELEASED_CODE */

  end subroutine latin_hypercube_driver

!-------------------------------------------------------------------------------
  subroutine latin_hypercube_2D_output &
             ( fname_prefix, fdir, stats_tout, nnzp, &
               zt, time_initial )
!-------------------------------------------------------------------------------

    use array_index, only: &
      iiLH_rrain, & ! Variables
      iiLH_Nr, &
      iiLH_Nc

    use parameters_model, only: &
      hydromet_dim ! Variable

    use parameters_microphys, only: &
      LH_microphys_calls ! Variable

    use stats_precision, only: &
      time_precision ! Constant

#ifdef UNRELEASED_CODE
    use output_2D_samples_module, only: &
      open_2D_samples_file ! Procedure

    use output_2D_samples_module, only: &
      lognormal_sample_file, & ! Instance of a type
      uniform_sample_file 

#endif /*UNRELEASED_CODE*/

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

    if ( l_output_2D_lognormal_dist .or. l_output_2D_uniform_dist ) then

      allocate( variable_names(hydromet_dim+5), variable_descriptions(hydromet_dim+5), &
                variable_units(hydromet_dim+5) )

      variable_names(1)        = "s_mellor"
      variable_descriptions(1) = "The variable 's' from Mellor 1977"
      variable_units(1)        = "kg/kg"

      variable_names(2)        = "t_mellor"
      variable_descriptions(2) = "The variable 't' from Mellor 1977"
      variable_units(2)        = "kg/kg"

      variable_names(3)        = "w"
      variable_descriptions(3) = "Vertical velocity"
      variable_units(3)        = "m/s"

      i = 3 ! Use i to determine the position of rt, thl in the output

      if ( iiLH_Nr > 0 ) then
        i = i + 1
        variable_names(iiLH_Nr)        = "Nr"
        variable_descriptions(iiLH_Nr) = "Rain droplet number concentration"
        variable_units(iiLH_Nr)        = "count/kg"
      end if
      if ( iiLH_Nc > 0 ) then
        i = i + 1
        variable_names(iiLH_Nc)        = "Nc"
        variable_descriptions(iiLH_Nc) = "Cloud droplet number concentration"
        variable_units(iiLH_Nc)        = "count/kg"
      end if
      if ( iiLH_rrain > 0 ) then
        i = i + 1
        variable_names(iiLH_rrain)        = "rrain"
        variable_descriptions(iiLH_rrain) = "Rain water mixing ratio"
        variable_units(iiLH_rrain)        = "kg/kg"
      end if

      i = i + 1
      variable_names(i)        = "rt"
      variable_descriptions(i) = "Total water mixing ratio"
      variable_units(i)        = "kg/kg"

      i = i + 1
      variable_names(i)        = "thl"
      variable_descriptions(i) = "Liquid potential temperature"
      variable_units(i)        = "K"
    end if
#ifdef UNRELEASED_CODE
    if ( l_output_2D_lognormal_dist ) then
      call open_2D_samples_file( nnzp, LH_microphys_calls, hydromet_dim+5, & ! In
                                 trim( fname_prefix )//"_nl", fdir, & ! In
                                 time_initial, stats_tout, zt, variable_names, & ! In
                                 variable_descriptions, variable_units, & ! In
                                 lognormal_sample_file ) ! In/Out
    end if

    if ( l_output_2D_uniform_dist ) then

      ! The uniform distribution corresponds to all the same variables as X_nl,
      ! except the d+1 component is the mixture component.
      ! We also derive LH_rt or LH_thl from s and t. (hence hydromet_dim+4)

      ! Overwrite the units
      variable_units(:) = "count" ! Unidata units format for a dimensionless quantity

      ! Overwrite for the mixture component
      variable_names(hydromet_dim+4) = "dp1"
      variable_descriptions(hydromet_dim+4) = "Uniform distribution for the mixture component"
      call open_2D_samples_file( nnzp, LH_microphys_calls, hydromet_dim+4, & ! In
                                 trim( fname_prefix )//"_u", fdir, & ! In
                                 time_initial, stats_tout, zt, &! In
                                 variable_names(1:hydromet_dim+4), & ! In
                                 variable_descriptions(1:hydromet_dim+4), & ! In
                                 variable_units(1:hydromet_dim+4), & ! In
                                 uniform_sample_file ) ! In/Out
    end if

    ! These should deallocate when we leave the scope of this subroutine, this
    ! is just in case.
    deallocate( variable_names, variable_descriptions, variable_units )

    return
#else
    stop "This code was not compiled with support for Latin Hypercube sampling"
#endif

  end subroutine latin_hypercube_2D_output

!-------------------------------------------------------------------------------
  subroutine latin_hypercube_2D_close
! Description:
!   Close a 2D sample file

! References:
!   None
!-------------------------------------------------------------------------------
#ifdef UNRELEASED_CODE
    use output_2D_samples_module, only: &
      close_2D_samples_file ! Procedure

    use output_2D_samples_module, only: &
      lognormal_sample_file, & ! Variable(s)
      uniform_sample_file 
#endif

    implicit none

    ! ---- Begin Code ----

#ifdef UNRELEASED_CODE
    if ( l_output_2D_lognormal_dist ) then
      call close_2D_samples_file( lognormal_sample_file )
    end if
    if ( l_output_2D_uniform_dist ) then
      call close_2D_samples_file( uniform_sample_file )
    end if
#else
    stop "This code was not compiled with support for Latin Hypercube sampling"
#endif

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

    use constants, only: &
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

#ifdef UNRELEASED_CODE
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

      if ( in_mixt_frac_1( X_u_dp1_element, real( mixt_frac, kind=genrand_real ) ) ) then
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

#endif

    return
  end subroutine choose_X_u_reject

!-------------------------------------------------------------------------------
  subroutine choose_X_u_scaled &
             ( l_cloudy_sample, cloud_frac1, &
              cloud_frac2, mixt_frac, &
              X_u_dp1_element, X_u_s_mellor_element )

! Description:
!   Find a clear or cloudy point for sampling.
!
! References:
!   None
!-------------------------------------------------------------------------------
    use mt95, only: genrand_real3 ! Procedure

    use mt95, only: genrand_real ! Constant

    use constants, only: &
      fstderr ! Constant

    implicit none

    ! External
    intrinsic :: dble

    ! Input Variables
    logical, intent(in) :: &
      l_cloudy_sample ! Whether his is a cloudy or clear air sample point

    real, intent(in) :: &
      cloud_frac1, &    ! Cloud fraction associated with mixture component 1     [-]
      cloud_frac2, &    ! Cloud fraction associated with mixture component 2     [-]
      mixt_frac         ! Mixture fraction                                       [-]

    ! Input/Output Variables
    real(kind=genrand_real), intent(inout) :: &
      X_u_dp1_element, X_u_s_mellor_element ! Elements from X_u (uniform dist.)

    ! Local Variables
    real :: cloud_frac_n, cloud_weighted_mixt_frac, clear_weighted_mixt_frac

    real(kind=genrand_real) :: rand, rand1, rand2 ! Random numbers

!   integer :: X_mixt_comp_one_lev


#ifdef UNRELEASED_CODE
    ! ---- Begin code ----

    ! Pick a new mixture component value between (0,1)
    call genrand_real3( rand1 )

    call genrand_real3( rand2 ) ! Determine a 2nd rand the if ... then

    if ( l_cloudy_sample ) then
      cloud_weighted_mixt_frac = ( mixt_frac*cloud_frac1 ) / &
                   ( mixt_frac*cloud_frac1 + (1.-mixt_frac)*cloud_frac2 )

      if ( in_mixt_frac_1( rand1, real( cloud_weighted_mixt_frac, kind=genrand_real ) ) ) then
        ! Component 1
        cloud_frac_n = cloud_frac1
!       X_mixt_comp_one_lev = 1
        X_u_dp1_element = mixt_frac * rand2
      else
        ! Component 2
        cloud_frac_n = cloud_frac2
!       X_mixt_comp_one_lev = 2
        X_u_dp1_element = mixt_frac + (1.-mixt_frac) * rand2
      end if
      call genrand_real3( rand ) ! Rand between (0,1)
      ! Scale and translate sample point to reside in cloud
      X_u_s_mellor_element = dble( cloud_frac_n * rand + (1.-cloud_frac_n) )

    else ! Clear air sample
      clear_weighted_mixt_frac = ( ( 1. - cloud_frac1 ) * mixt_frac ) &
        / ( ( 1.-cloud_frac1 ) * mixt_frac + ( 1.-cloud_frac2 )*( 1.-mixt_frac ) )

      if ( in_mixt_frac_1( rand1, dble( clear_weighted_mixt_frac ) ) ) then
        ! Component 1
        cloud_frac_n = cloud_frac1
!       X_mixt_comp_one_lev = 1
        X_u_dp1_element = mixt_frac * rand2
      else
        ! Component 2
        cloud_frac_n = cloud_frac2
!       X_mixt_comp_one_lev = 2
        X_u_dp1_element = mixt_frac + (1.-mixt_frac) * rand2
      end if
      call genrand_real3( rand ) ! Rand between (0,1)
      ! Scale and translate sample point to reside in clear air (no cloud)
      X_u_s_mellor_element = dble( (1.-cloud_frac_n) * rand )

    end if

#endif

    return
  end subroutine choose_X_u_scaled

!----------------------------------------------------------------------
  elemental function in_mixt_frac_1( X_u_dp1_element, frac )

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

    logical :: in_mixt_frac_1

    ! ---- Begin Code ----

    if ( X_u_dp1_element < frac ) then
      in_mixt_frac_1 = .true.
    else
      in_mixt_frac_1 = .false.
    end if

    return
  end function in_mixt_frac_1

end module latin_hypercube_driver_module
