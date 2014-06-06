!-------------------------------------------------------------------------------
! $Id$
!===============================================================================
module estimate_scm_microphys_module

  implicit none

  public :: est_single_column_tndcy, copy_X_nl_into_hydromet_all_pts, &
            copy_X_nl_into_rc_all_pts

  private ! Default scope

  contains

!-------------------------------------------------------------------------------
  subroutine est_single_column_tndcy &
             ( dt, nz, num_samples, d_variables, &
               k_lh_start, lh_rt, lh_thl, &
               X_nl_all_levs, lh_sample_point_weights, &
               p_in_Pa, exner, rho, cloud_frac, w_std_dev, &
               dzq, hydromet, rcm, Nc_in_cloud,  &
               lh_hydromet_mc, lh_hydromet_vel, lh_Ncm_mc, &
               lh_rvm_mc, lh_rcm_mc, lh_thlm_mc, &
               lh_rtp2_mc, lh_thlp2_mc, lh_wprtp_mc, &
               lh_wpthlp_mc, lh_rtpthlp_mc, &
               microphys_sub )
! Description:
!   Estimate the tendency of a microphysics scheme via latin hypercube sampling
!
! References:
!   None
!-------------------------------------------------------------------------------

    use constants_clubb, only:  &
      fstderr, &  ! Constant(s)
      zero_threshold

    use parameters_model, only: &
      hydromet_dim ! Variable

    use parameters_microphys, only: &
      l_lh_cloud_weighted_sampling, & ! Variable(s)
      l_const_Nc_in_cloud, &
      l_var_covar_src


    use corr_matrix_module, only: &
      iiPDF_s_mellor, &
      iiPDF_w

    use error_code, only: &
      clubb_at_least_debug_level ! Procedure

    use clubb_precision, only: &
      dp, & ! double precision
      core_rknd, &
      time_precision

    use stats_variables, only: &
      zt,           &
      sfc,          &
      l_stats_samp ! Variable(s)

    use microphys_stats_vars_module, only: &
      microphys_stats_vars_type, &
      microphys_stats_accumulate, &
      microphys_stats_cleanup

    use parameters_microphys, only: &
      l_silhs_KK_convergence_adj_mean ! Variable(s)

    use math_utilities, only: &
      compute_sample_mean             ! Procedure

#ifdef SILHS_KK_CONVERGENCE_TEST
    use array_index, only: &
      iiNrm, & ! Variable(s)
      iirrainm
#endif

    implicit none

    ! External
#include "microphys_interface.inc"

    intrinsic :: real, all, any

    ! Constant parameters
    logical, parameter :: &
      l_latin_hypercube = .true. ! We are the Latin hypercube!

    logical, parameter :: &
      l_check_lh_cloud_weighting = .true. ! Verify every other sample point is out of cloud

    ! Input Variables
    real( kind = time_precision ), intent(in) :: &
      dt ! Model timestep       [s]

    integer, intent(in) :: &
      nz,            & ! Number of vertical levels
      num_samples, & ! Number of calls to microphysics (normally=2)
      d_variables,   & ! Number of variates (normally=5) 
      k_lh_start       ! Starting level for computing arbitrary overlap

    real( kind = core_rknd ), dimension(nz,num_samples), intent(in) :: &
      lh_rt, & ! num_samples values of total water mixing ratio     [kg/kg]
      lh_thl   ! num_samples values of liquid potential temperature [K]

    real( kind = dp ), target, dimension(nz,num_samples,d_variables), intent(in) :: &
      X_nl_all_levs ! Sample that is transformed ultimately to normal-lognormal

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      p_in_Pa,    & ! Pressure                 [Pa]
      exner,      & ! Exner function           [-]
      rho,        & ! Density on thermo. grid  [kg/m^3]
      cloud_frac, & ! Cloud fraction           [-]
      w_std_dev,  & ! Standard deviation of w    [m/s]
      dzq,        & ! Difference in height per gridbox   [m]
      rcm           ! Mean liquid water mixing ratio     [kg/kg]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      ! Constant value of N_c within cloud, to be used with l_const_Nc_in_cloud
      Nc_in_cloud 

    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(in) :: &
      hydromet ! Hydrometeor species    [units vary]

    real( kind = core_rknd ), dimension(num_samples), intent(in) :: &
       lh_sample_point_weights ! Weight for cloud weighted sampling

    ! Output Variables

    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(out) :: &
      lh_hydromet_mc, & ! LH estimate of hydrometeor time tendency          [(units vary)/s]
      lh_hydromet_vel   ! LH estimate of hydrometeor sedimentation velocity [m/s]

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


    ! Local Variables
    real( kind = core_rknd ), dimension(nz,hydromet_dim,num_samples) :: &
      lh_hydromet_mc_all, & ! LH est of hydrometeor time tendency          [(units vary)/s]
      lh_hydromet_vel_all   ! LH est of hydrometeor sedimentation velocity [m/s]

    real( kind = core_rknd ), dimension(nz,num_samples) :: &
      lh_Ncm_mc_all,       & ! LH est of time tendency of cloud droplet concentration  [#/kg/s]
      lh_rcm_mc_all,       & ! LH est of time tendency of liquid water mixing ratio    [kg/kg/s]
      lh_rvm_mc_all,       & ! LH est of time tendency of vapor water mixing ratio     [kg/kg/s]
      lh_thlm_mc_all         ! LH est of time tendency of liquid potential temperature     [K/s]

    real( kind = core_rknd ), dimension(nz,hydromet_dim) :: &
      hydromet_all_points ! Hydrometeor species                    [units vary]

    real( kind = core_rknd ), dimension(nz) :: &
      Ncn_all_points    ! Cloud Nuclei conc. (simplified); Nc=Ncn*H(s)   [#/kg]

    real( kind = core_rknd ), dimension(nz) :: &
      s_mellor_column    ! 's' (Mellor 1977)                    [kg/kg]

    real( kind = core_rknd ), dimension(nz) :: &
      rv_column,  & ! Vapor water                               [kg/kg]
      thl_column, & ! Liquid potential temperature              [K]
      rc_column,  & ! Liquid water                              [kg/kg]
      w_column,   & ! Vertical velocity                         [m/s]
      Nc            ! Cloud droplet concentration               [#/kg]

    real( kind = core_rknd ), dimension(nz) :: &
      lh_rtp2_before_microphys,    & ! <rt'^2> before microphys_sub    [(kg/kg)^2]
      lh_rtp2_after_microphys,     & ! <rt'^2> after microphys_sub     [(kg/kg)^2]
      lh_thlp2_before_microphys,   & ! <thl'^2> before microphys_sub   [K^2]
      lh_thlp2_after_microphys,    & ! <thl'^2> after microphys_sub    [K^2]
      lh_wprtp_before_microphys,   & ! <w'rt'> before microphys_sub    [m*(kg/kg)/s]
      lh_wprtp_after_microphys,    & ! <w'rt'> after microphys_sub     [m*(kg/kg)/s]
      lh_wpthlp_before_microphys,  & ! <w'thl'> before microphys_sub   [m*K/s]
      lh_wpthlp_after_microphys,   & ! <w'thl'> after microphys_sub    [m*K/s]
      lh_rtpthlp_before_microphys, & ! <rt'thl'> before microphys_sub  [K*(kg/kg)/s]
      lh_rtpthlp_after_microphys     ! <rt'thl'> after microphys_sub   [K*(kg/kg)/s]

    real( kind = core_rknd ), dimension(nz,num_samples) :: &
      rt_all_samples,  & ! Columns used to calculate covariances [kg/kg]
      thl_all_samples, & ! Columns used to calculate covariances [K]
      w_all_samples      ! Columns used to calculate covariances [m/s] 

    real( kind = dp ), pointer, dimension(:,:) :: &
      s_mellor_all_points,  & ! num_samples values of 's' (Mellor 1977)      [kg/kg]
      w_all_points            ! num_samples values of vertical velocity      [m/s]

    type(microphys_stats_vars_type), dimension(num_samples) :: &
      microphys_stats_zt_all, &   ! Statistics variables output from microphysics on zt grid
      microphys_stats_sfc_all     ! Statistics variables on sfc grid

    type(microphys_stats_vars_type) :: &
      microphys_stats_zt_avg, &   ! Average of zt statistics over all sample points
      microphys_stats_sfc_avg     ! Average of sfc statistics over all sample points

    integer :: ivar, k, sample

    integer :: &
      in_cloud_points, &
      out_of_cloud_points

    ! ---- Begin Code ----

    ! Mellor's 's' is hardwired elsewhere to be the first column
    s_mellor_all_points => X_nl_all_levs(:,:,iiPDF_s_mellor)
    w_all_points        => X_nl_all_levs(:,:,iiPDF_w)

    ! Assertion check
    if ( l_check_lh_cloud_weighting .and. l_lh_cloud_weighted_sampling .and. &
         all( lh_sample_point_weights(:) /= 1.0_core_rknd ) ) then 
        ! The 1.0 indicates cloud_frac is > 0.5

      ! Verify every other sample point is out of cloud if we're doing
      ! cloud weighted sampling
      in_cloud_points     = 0
      out_of_cloud_points = 0
      do sample = 1, num_samples, 1
        if ( s_mellor_all_points(k_lh_start,sample) > 0._dp ) then
          in_cloud_points = in_cloud_points + 1
        else if ( s_mellor_all_points(k_lh_start,sample) <= 0._dp ) then
          out_of_cloud_points = out_of_cloud_points + 1
        end if
      end do ! 1..num_samples
      if ( in_cloud_points /= out_of_cloud_points ) then
        if ( clubb_at_least_debug_level( 2 ) ) then
          write(fstderr,*) "In est_single_column_tndcy:"
          write(fstderr,*) "The cloudy sample points do not equal the out of cloud points"
          write(fstderr,*) "in_cloud_points =", in_cloud_points
          write(fstderr,*) "out_of_cloud_points =", out_of_cloud_points
          write(fstderr,*)  "cloud fraction = ", cloud_frac(k_lh_start)
          write(fstderr,*) "k_lh_start = ", k_lh_start, "nz = ", nz
          write(fstderr,'(4X,A,A)')  "s_mellor_all_points  ", "weight   "
          do sample = 1, num_samples, 1
            write(fstderr,'(I4,2G20.4)') &
              sample, s_mellor_all_points(k_lh_start,sample), lh_sample_point_weights(sample)
          end do
        end if ! clubb_at_least_debug_level 2
        stop "Fatal Error in est_single_column_tndcy "
     end if  ! in_cloud_points /= out_of_cloud_points
    end if ! l_check_lh_cloud_weighting .and. l_lh_cloud_weighted_sampling

    lh_hydromet_vel(:,:) = 0._core_rknd

    ! Initialize microphysical tendencies for each mixture component
    lh_hydromet_mc_all(:,:,:) = 0._core_rknd

    lh_hydromet_vel_all(:,:,:) = 0._core_rknd

    lh_Ncm_mc_all(:,:) = 0._core_rknd

    lh_rcm_mc_all(:,:) = 0._core_rknd

    lh_rvm_mc_all(:,:) = 0._core_rknd

    lh_thlm_mc_all(:,:) = 0._core_rknd

    lh_rtp2_mc(:) = 0.0_core_rknd
    lh_thlp2_mc(:) = 0.0_core_rknd
    lh_wprtp_mc(:) = 0.0_core_rknd
    lh_wpthlp_mc(:) = 0.0_core_rknd
    lh_rtpthlp_mc(:) = 0.0_core_rknd

    lh_rtp2_before_microphys(:) = 0.0_core_rknd
    lh_thlp2_before_microphys(:) = 0.0_core_rknd
    lh_wprtp_before_microphys(:) = 0.0_core_rknd
    lh_wpthlp_before_microphys(:) = 0.0_core_rknd
    lh_rtpthlp_before_microphys(:) = 0.0_core_rknd
    lh_rtp2_after_microphys(:) = 0.0_core_rknd
    lh_thlp2_after_microphys(:) = 0.0_core_rknd
    lh_wprtp_after_microphys(:) = 0.0_core_rknd
    lh_wpthlp_after_microphys(:) = 0.0_core_rknd
    lh_rtpthlp_after_microphys(:) = 0.0_core_rknd


    if ( l_var_covar_src ) then
      rt_all_samples = lh_rt
      thl_all_samples = lh_thl
      w_all_samples = w_all_points

      call lh_moments ( num_samples, lh_sample_point_weights, nz, &           ! Intent (in)
                       rt_all_samples, thl_all_samples, w_all_samples, &        ! Intent (in)
                       lh_rtp2_before_microphys, lh_thlp2_before_microphys, &   ! Intent (out)
                       lh_wprtp_before_microphys, lh_wpthlp_before_microphys, & ! Intent (out)
                       lh_rtpthlp_before_microphys )                            ! Intent (out)
    end if

    do sample = 1, num_samples

      s_mellor_column = real( s_mellor_all_points(:,sample), kind = core_rknd )

      where ( s_mellor_all_points(:,sample) > 0.0_dp )
        rc_column = real( s_mellor_all_points(:,sample), kind = core_rknd )
      else where
        rc_column = 0.0_core_rknd
      end where

      w_column   = real( w_all_points(:,sample), kind = core_rknd )
      rv_column  = real( lh_rt(:,sample), kind = core_rknd ) - rc_column
      ! Verify total water isn't negative
      if ( any( rv_column < 0._core_rknd) ) then
        if ( clubb_at_least_debug_level( 1 ) ) then
          write(fstderr,*) "rv negative, LH sample number = ", sample
          write(fstderr,'(a3,3a20)') "k", "rt", "rv", "rc"
          do k = 1, nz
            if ( rv_column(k) < 0._core_rknd) then
              write(6,'(i3,3g20.7)')  k, lh_rt(k,sample), rv_column(k), &
                rc_column(k)
            end if
          end do
        end if ! clubb_at_least_debug_level( 1 )
        write(fstderr,*) "Applying non-conservative hard clipping to rv sample."
        where ( rv_column < 0._core_rknd) rv_column = zero_threshold
      end if ! Some rv_column element < 0

      thl_column = real( lh_thl(:,sample), kind = core_rknd )

      call copy_X_nl_into_hydromet_all_pts( nz, d_variables, 1, & ! In
                                    X_nl_all_levs(:,sample,:), & ! In
                                    hydromet, & ! In
                                    hydromet_all_points, &  ! Out
                                    Ncn_all_points ) ! Out

      ! Set Nc = Ncn * H(s).  For l_const_Nc_in_cloud, Ncn should have a
      ! constant value of Nc_in_cloud.
      if ( l_const_Nc_in_cloud ) then
         ! For l_const_Nc_in_cloud, we want to use the same value of Nc for all
         ! sample points.
         where ( s_mellor_column > 0.0_core_rknd )
            Nc = Nc_in_cloud
         else where
            Nc = 0.0_core_rknd
         end where
      else ! Nc varies in-cloud
         where ( s_mellor_column > 0.0_core_rknd )
            Nc = Ncn_all_points
         else where
            Nc = 0.0_core_rknd
         end where
      endif

      ! Call the microphysics scheme to obtain a sample point
      call microphys_sub &
           ( dt, nz, & ! In
             l_latin_hypercube, thl_column, w_column, p_in_Pa, & ! In
             exner, rho, cloud_frac, w_std_dev, & ! In
             dzq, rc_column, Nc, s_mellor_column, rv_column, & ! In
             hydromet_all_points, & ! In
             lh_hydromet_mc_all(:,:,sample), lh_hydromet_vel_all(:,:,sample), & ! Out
             lh_Ncm_mc_all(:,sample), & ! Out
             lh_rcm_mc_all(:,sample), lh_rvm_mc_all(:,sample), lh_thlm_mc_all(:,sample), & ! Out
             microphys_stats_zt_all(sample), microphys_stats_sfc_all(sample) ) ! Out

      rt_all_samples(:,sample) = rc_column + rv_column + &
        real( dt, kind=core_rknd ) * ( lh_rcm_mc_all(:,sample) + lh_rvm_mc_all(:,sample) )

      thl_all_samples(:,sample) = thl_column + real( dt, kind=core_rknd ) * lh_thlm_mc_all(:,sample)

      ! Loop to get new sample
    end do ! sample = 1, num_samples

    if ( l_var_covar_src ) then
      call lh_moments ( num_samples, lh_sample_point_weights, nz, &           ! Intent (in)
                         rt_all_samples, thl_all_samples, w_all_samples, &      ! Intent (in)
                         lh_rtp2_after_microphys, lh_thlp2_after_microphys, &   ! Intent (out)
                         lh_wprtp_after_microphys, lh_wpthlp_after_microphys, & ! Intent (out)
                         lh_rtpthlp_after_microphys )                           ! Intent (out)

      lh_wpthlp_mc = ( lh_wpthlp_after_microphys - lh_wpthlp_before_microphys ) / dt
      lh_wprtp_mc = ( lh_wprtp_after_microphys - lh_wprtp_before_microphys ) / dt
      lh_rtp2_mc = ( lh_rtp2_after_microphys - lh_rtp2_before_microphys ) / dt 
      lh_thlp2_mc = ( lh_thlp2_after_microphys - lh_thlp2_before_microphys) / dt
      lh_rtpthlp_mc = ( lh_rtpthlp_after_microphys - lh_rtpthlp_before_microphys) / dt

    end if

    ! Grid box average.
    forall( ivar = 1:hydromet_dim )

      lh_hydromet_vel(:,ivar) &
        = compute_sample_mean( nz, num_samples, lh_sample_point_weights, &
                               lh_hydromet_vel_all(:,ivar,:) )

      lh_hydromet_mc(:,ivar) &
        = compute_sample_mean( nz, num_samples, lh_sample_point_weights, &
                               lh_hydromet_mc_all(:,ivar,:) )

    end forall

    lh_Ncm_mc = compute_sample_mean( nz, num_samples, lh_sample_point_weights, lh_Ncm_mc_all(:,:) )
    lh_rcm_mc = compute_sample_mean( nz, num_samples, lh_sample_point_weights, lh_rcm_mc_all(:,:) )
    lh_rvm_mc = compute_sample_mean( nz, num_samples, lh_sample_point_weights, lh_rvm_mc_all(:,:) )
    lh_thlm_mc= compute_sample_mean( nz, num_samples, lh_sample_point_weights, lh_thlm_mc_all(:,:))

    ! Sample variables from microphys_stats_vars objects for statistics

    microphys_stats_zt_avg =  silhs_microphys_stats_avg( num_samples, microphys_stats_zt_all, &
                                                         lh_sample_point_weights )
    microphys_stats_sfc_avg = silhs_microphys_stats_avg( num_samples, microphys_stats_sfc_all, &
                                                         lh_sample_point_weights )

    call microphys_stats_accumulate( microphys_stats_zt_avg, l_stats_samp, zt )
    call microphys_stats_accumulate( microphys_stats_sfc_avg, l_stats_samp, sfc )

#ifdef SILHS_KK_CONVERGENCE_TEST
    ! Adjust the mean if l_silhs_KK_convergence_adj_mean is true
    if ( l_silhs_KK_convergence_adj_mean ) then
      call adjust_KK_src_means( dt, nz, exner, rcm, hydromet(:,iirrainm),           & ! intent(in)
                                hydromet(:,iiNrm),                                  & ! intent(in)
                                microphys_stats_zt_avg, l_stats_samp,               & ! intent(in)
                                lh_hydromet_mc(:,iirrainm), lh_hydromet_mc(:,iiNrm),& ! intent(out)
                                lh_rvm_mc, lh_rcm_mc, lh_thlm_mc )                    ! intent(out)
    end if
#else
    ! Eliminate the resulting compiler warning
    if (.false.) then
      Nc(1) = rcm(1)
    end if
#endif

    ! Cleanup microphys_stats_vars objects
    do ivar=1, num_samples
      call microphys_stats_cleanup( microphys_stats_zt_all(ivar) )
      call microphys_stats_cleanup( microphys_stats_sfc_all(ivar) )
    end do
    call microphys_stats_cleanup( microphys_stats_zt_avg )
    call microphys_stats_cleanup( microphys_stats_sfc_avg )

    return
  end subroutine est_single_column_tndcy

  !-----------------------------------------------------------------------
  function silhs_microphys_stats_avg &
           ( num_samples, microphys_stats_all, lh_sample_point_weights ) &
           result( microphys_stats_avg )

  ! Description:
  !   Computes a weighted average of all statistical variables 

  ! References:
  !   None
  !-----------------------------------------------------------------------

    ! Included Modules
    use clubb_precision, only: &
      core_rknd            ! Compile-time constant

    use microphys_stats_vars_module, only: &
      microphys_stats_vars_type, & ! Type
      microphys_put_var,         & ! Procedure(s)
      microphys_stats_alloc

    use math_utilities, only: &
      compute_sample_mean

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      num_samples           ! Number of SILHS sample points

    type(microphys_stats_vars_type), dimension(num_samples), intent(in) :: &
      microphys_stats_all   ! Statistics variables

    real( kind = core_rknd ), dimension(num_samples), intent(in) :: &
      lh_sample_point_weights ! Weight of each SILHS sample point

    ! Output Variables
    type(microphys_stats_vars_type) :: &
      microphys_stats_avg   ! Average of statistical variables

    ! Local Variables
    integer :: ivar, svar, nz, num_vars, stats_index

    real( kind = core_rknd ), dimension(microphys_stats_all(1)%nz,num_samples) :: &
      stats_var_all

    real( kind = core_rknd ), dimension(microphys_stats_all(1)%nz) :: &
      stats_var_avg

  !-----------------------------------------------------------------------
    !----- Begin Code-----

    nz = microphys_stats_all(1)%nz
    num_vars = microphys_stats_all(1)%num_vars

    call microphys_stats_alloc( nz, num_vars, microphys_stats_avg )

    do ivar=1, num_vars

      stats_index = microphys_stats_all(1)%stats_indices(ivar)

      do svar=1, num_samples
        stats_var_all(:,svar) = microphys_stats_all(svar)%output_values(:,ivar)
      end do

      stats_var_avg = compute_sample_mean( nz, num_samples, lh_sample_point_weights, stats_var_all )

      call microphys_put_var( stats_index, stats_var_avg, microphys_stats_avg )

    end do

    return
  end function silhs_microphys_stats_avg
  !-----------------------------------------------------------------------

#ifdef SILHS_KK_CONVERGENCE_TEST
  !-----------------------------------------------------------------------------
  subroutine adjust_KK_src_means( dt, nz, exner, rcm, rrainm, Nrm,         &
                                  microphys_stats_zt, l_stats_samp,        &
                                  rrainm_mc, Nrm_mc,                       &
                                  rvm_mc, rcm_mc, thlm_mc )

  ! Description:
  !   Adjusts the means of microphysics terms for KK microphysics by calling the
  !   KK microphysics adjustment subroutine

  ! References:
  !   clubb:ticket:558
  !-----------------------------------------------------------------------------
    use KK_Nrm_tendencies, only: &
      KK_Nrm_auto_mean, & ! Procedure(s)
      KK_Nrm_evap_local_mean

    use KK_microphys_module, only: &
      KK_microphys_adjust, &      ! Procedure
      KK_microphys_adj_terms_type ! Type

    use clubb_precision, only: &
      time_precision, &
      core_rknd

    use constants_clubb, only: &
      rr_tol, & ! Constant(s)
      Nr_tol, &
      zero

    use stats_type, only: &
      stat_update_var ! Procedure

    use stats_variables, only: &
      zt, &
      irrainm_auto, &
      irrainm_accr, &
      irrainm_cond, &
      iNrm_auto, &
      iNrm_cond, &
      irrainm_src_adj, &
      iNrm_src_adj, &
      irrainm_cond_adj, &
      iNrm_cond_adj

    use microphys_stats_vars_module, only: &
      microphys_stats_vars_type, &     ! Type
      microphys_get_var                ! Procedure

    implicit none

    ! Local Constants
    logical, parameter :: &
      ! Whether to adjust rrainm_source to not over-deplete cloud water
      l_src_adj_enabled = .true.,  &
      ! Whether to adjust rrainm_evap to not over-evaporate rain
      l_evap_adj_enabled = .true., &
      ! This subroutine is called from Latin hypercube.
      l_latin_hypercube = .true.

    ! Input variables
    real( kind = time_precision ), intent(in) :: &
      dt   ! Model timestep

    integer, intent(in) :: &
      nz   ! Number of vertical grid levels

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      exner,       & ! Exner function                            [-]
      rcm,         & ! Mean liquid water mixing ratio            [kg/kg]
      rrainm,      & ! Rain water mixing ration                  [kg/kg]
      Nrm            ! Rain drop concentration                   [num/kg]

    type(microphys_stats_vars_type), intent(in) :: &
      microphys_stats_zt ! Statistics variables                  [units vary]

    logical, intent(in) :: &
      l_stats_samp   ! Whether to sample this timestep

    ! Output variables
    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      rrainm_mc, & ! Mean change in rain due to microphysics [(kg/kg)/s] 
      Nrm_mc,    & ! Mean change in Nrm due to microphysics  [(kg/kg)/s]
      rvm_mc,    & ! Time tendency of rvm                    [(kg/kg)/s]
      rcm_mc,    & ! Time tendency of rcm                    [(kg/kg)/s]
      thlm_mc      ! Time tendency of thlm                   [(kg/kg)/s]

    ! Local Variables
    real( kind = core_rknd ), dimension(nz) :: &
      rrainm_evap, & ! Mean change in rain due to evap           [(kg/kg)/s]
      rrainm_auto, & ! Mean change in rain due to autoconversion [(kg/kg)/s]
      rrainm_accr, & ! Mean change in rain due to accretion      [(kg/kg)/s]
      Nrm_auto,    & ! Mean change in Nrm due to autoconversion  [(num/kg)/s]
      Nrm_evap       ! Mean change in Nrm due to evaporation     [(num/kg)/s]

    type(KK_microphys_adj_terms_type), dimension(nz) :: &
      adj_terms    ! Adjustment terms returned from the adjustment routine

    integer :: k

    !----- Begin code -----

    ! Initialize output
    rrainm_mc = zero
    Nrm_mc = zero
    rvm_mc = zero
    rcm_mc = zero
    thlm_mc = zero

    rrainm_auto = microphys_get_var( irrainm_auto, microphys_stats_zt )
    rrainm_accr = microphys_get_var( irrainm_accr, microphys_stats_zt )
    rrainm_evap = microphys_get_var( irrainm_cond, microphys_stats_zt )

    Nrm_auto    = microphys_get_var( iNrm_auto,    microphys_stats_zt )
    Nrm_evap    = microphys_get_var( iNrm_cond,    microphys_stats_zt )

    ! Loop over each vertical level above the lower boundary
    do k = 2, nz, 1

      ! We call KK_microphys_adjust to adjust the means of the mc terms
      call KK_microphys_adjust( dt, exner(k), rcm(k), rrainm(k), Nrm(k),   & !intent(in)
                                rrainm_evap(k), rrainm_auto(k),            & !intent(in)
                                rrainm_accr(k), Nrm_evap(k),               & !intent(in)
                                Nrm_auto(k), l_src_adj_enabled,            & !intent(in)
                                l_evap_adj_enabled,                        & !intent(in)
                                l_latin_hypercube, k,                      & !intent(in)
                                rrainm_mc(k), Nrm_mc(k),                   & !intent(out)
                                rvm_mc(k), rcm_mc(k), thlm_mc(k),          & !intent(out)
                                adj_terms(k) )                               !intent(out)
    end do ! k = 2, nz, 1

    ! Set boundary conditions
    rrainm_mc(1) = zero
    rrainm_mc(nz) = zero

    Nrm_mc(1) = zero
    Nrm_mc(nz) = zero

    rvm_mc(1) = zero
    rvm_mc(nz) = zero

    rcm_mc(1) = zero
    rcm_mc(nz) = zero

    thlm_mc(1) = zero
    thlm_mc(nz) = zero

    ! Statistical sampling
    if ( l_stats_samp ) then

      if ( irrainm_src_adj > 0 ) then
        call stat_update_var( irrainm_src_adj, adj_terms%rrainm_src_adj, zt )
      end if

      if ( iNrm_src_adj > 0 ) then
        call stat_update_var( iNrm_src_adj, adj_terms%Nrm_src_adj, zt )
      end if

      if ( irrainm_cond_adj > 0 ) then
        call stat_update_var( irrainm_cond_adj, adj_terms%rrainm_cond_adj, zt )
      end if

      if ( iNrm_cond_adj > 0 ) then
        call stat_update_var( iNrm_cond_adj, adj_terms%Nrm_cond_adj, zt )
      end if

    end if ! l_stats_samp

  end subroutine adjust_KK_src_means
#endif /*SILHS_KK_CONVERGENCE_TEST*/
  !-----------------------------------------------------------------------------
  subroutine copy_X_nl_into_hydromet_all_pts( nz, d_variables, num_samples, &
                                      X_nl_all_levs, &
                                      hydromet, &
                                      hydromet_all_points, &
                                      Ncn_all_points )

  ! Description:
  !   Copy the points from the latin hypercube sample to an array with just the
  !   hydrometeors
  ! References:
  !   None
  !-----------------------------------------------------------------------------
    use parameters_model, only: &
      hydromet_dim ! Variable

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
      iiPDF_rrain, &
      iiPDF_rsnow, &
      iiPDF_rice, &
      iiPDF_rgraupel, &
      iiPDF_Nr, &
      iiPDF_Nsnow, &
      iiPDF_Ngraupel, &
      iiPDF_Ncn, &
      iiPDF_Ni

    use clubb_precision, only: &
      dp, & ! double precision
      core_rknd

    implicit none

    integer, intent(in) :: &
      nz,            & ! Number of vertical levels
      d_variables,   & ! Number of variates
      num_samples    ! Number of calls to microphysics

    real( kind = dp ), dimension(nz,num_samples,d_variables), intent(in) :: &
      X_nl_all_levs ! Sample that is transformed ultimately to normal-lognormal

    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(in) :: &
      hydromet ! Hydrometeor species    [units vary]

    real( kind = core_rknd ), dimension(nz,num_samples,hydromet_dim), intent(out) :: &
      hydromet_all_points ! Hydrometeor species    [units vary]

    real( kind = core_rknd ), dimension(nz,num_samples), intent(out) :: &
      Ncn_all_points    ! Cloud nuclei conc. (simplified); Nc=Ncn*H(s)   [#/kg]

    integer :: sample, ivar

    do sample = 1, num_samples
      ! Copy the sample points into the temporary arrays
      do ivar = 1, hydromet_dim, 1
        if ( ivar == iirrainm .and. iiPDF_rrain > 0 ) then
          ! Use a sampled value of rain water mixing ratio
          hydromet_all_points(:,sample,ivar) = &
            real( X_nl_all_levs(:,sample,iiPDF_rrain), kind = core_rknd )

        else if ( ivar == iirsnowm .and. iiPDF_rsnow > 0 ) then
          ! Use a sampled value of rain water mixing ratio
          hydromet_all_points(:,sample,ivar) = &
            real( X_nl_all_levs(:,sample,iiPDF_rsnow), kind = core_rknd )

        else if ( ivar == iiricem .and. iiPDF_rice > 0 ) then
          ! Use a sampled value of rain water mixing ratio
          hydromet_all_points(:,sample,ivar) = &
            real( X_nl_all_levs(:,sample,iiPDF_rice), kind = core_rknd )

        else if ( ivar == iirgraupelm .and. iiPDF_rgraupel > 0 ) then
          ! Use a sampled value of rain water mixing ratio
          hydromet_all_points(:,sample,ivar) = &
            real( X_nl_all_levs(:,sample,iiPDF_rgraupel), kind = core_rknd )

        else if ( ivar == iiNrm .and. iiPDF_Nr > 0 ) then
          ! Use a sampled value of rain droplet number concentration
          hydromet_all_points(:,sample,ivar) = &
            real( X_nl_all_levs(:,sample,iiPDF_Nr), kind = core_rknd )

        else if ( ivar == iiNsnowm .and. iiPDF_Nsnow > 0 ) then
          ! Use a sampled value of rain droplet number concentration
          hydromet_all_points(:,sample,ivar) = &
            real( X_nl_all_levs(:,sample,iiPDF_Nsnow), kind = core_rknd )

        else if ( ivar == iiNgraupelm .and. iiPDF_Ngraupel > 0 ) then
          ! Use a sampled value of rain droplet number concentration
          hydromet_all_points(:,sample,ivar) = &
            real( X_nl_all_levs(:,sample,iiPDF_Ngraupel), kind = core_rknd )

        else if ( ivar == iiNim .and. iiPDF_Ni > 0 ) then
          ! Use a sampled value of rain droplet number concentration
          hydromet_all_points(:,sample,ivar) = &
            real( X_nl_all_levs(:,sample,iiPDF_Ni), kind = core_rknd )

        else ! Use the mean field, rather than a sample point
          ! This is the case for hail and graupel in the Morrison microphysics
          ! currently -dschanen 23 March 2010
          hydromet_all_points(:,sample,ivar) = hydromet(:,ivar)

        end if
      end do ! 1..hydromet_dim
      ! Copy Ncn into Ncn all points
      if ( iiPDF_Ncn > 0 ) then
        Ncn_all_points(:,sample) = &
          real( X_nl_all_levs(:,sample,iiPDF_Ncn), kind=core_rknd )
      end if
    end do ! 1..num_samples

    return
  end subroutine copy_X_nl_into_hydromet_all_pts
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  subroutine copy_X_nl_into_rc_all_pts &
             ( nz, d_variables, num_samples, X_nl_all_levs, &
               lh_rc )
  ! Description:
  !   Extracts a sample of rc from X_nl_all_levs, where rc = s * H(s), for each
  !   subcolumn.

  ! References:
  !   none
  !-----------------------------------------------------------------------------
    use clubb_precision, only: &
      dp, &
      core_rknd

    use corr_matrix_module, only: &
      iiPDF_s_mellor        ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz, &            ! Number of vertical levels
      d_variables, &   ! Number of lognormal variates
      num_samples      ! Number of SILHS samples per variate

    real( kind = dp ), dimension(nz,num_samples,d_variables), intent(in) :: &
      X_nl_all_levs    ! Normal-lognormal SILHS sample   [units vary]

    ! Output variables
    real( kind = core_rknd ), dimension(nz,num_samples), intent(out) :: &
      lh_rc            ! SILHS samples of rc       [kg/kg]

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    where ( X_nl_all_levs(:,:,iiPDF_s_mellor) >= 0.0_dp )
      lh_rc = real( X_nl_all_levs(:,:,iiPDF_s_mellor), kind=core_rknd )
    elsewhere
      lh_rc = 0.0_core_rknd
    end where

    return
  end subroutine copy_X_nl_into_rc_all_pts
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  subroutine lh_moments ( n_samples, lh_weights, nz, & 
                 rt_all_samples, thl_all_samples, w_all_samples, &
                 lh_rtp2, lh_thlp2, &
                 lh_wprtp, lh_wpthlp, &
                 lh_rtpthlp )

  ! Description:
  !   Calculates variances and covariances using LH sample columns
  !
  ! References:
  !   None
  !
  ! TODO: This code assumes nz == gr%nnzp since it references zt2zm;  this is
  ! not necessarily the case.
  !-----------------------------------------------------------------------------

    use grid_class, only: &
      zt2zm    ! Procedures

    use math_utilities, only: &
      compute_sample_mean, & ! functions
      compute_sample_variance, &
      compute_sample_covariance

    use clubb_precision, only: &
      core_rknd ! Constants

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz,            & ! Number of vertical levels
      n_samples        ! Number of sample columns from latin hypercube

    real( kind = core_rknd ), dimension(n_samples), intent(in) :: &
      lh_weights   ! Sample  point weights                   [-]

    real( kind = core_rknd ), dimension(nz,n_samples), intent(in) :: &
      rt_all_samples, &  ! rt columns from latin hypercube   [kg/kg]
      thl_all_samples, & ! thl columns from latin hypercube  [K]
      w_all_samples      ! w columns from latin hypercube    [m/s]

    ! Output Variables
    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      lh_rtp2, &    ! Latin hypercube estimate of <rt'^2>     [(kg/kg)^2]
      lh_thlp2, &   ! Latin hypercube estimate of <thl'^2>    [K^2]
      lh_wprtp, &   ! Latin hypercube estimate of <w'rt'>     [m*(kg/kg)/s]
      lh_wpthlp, &  ! Latin hypercube estimate of <w'thl'>    [m*K/s]
      lh_rtpthlp    ! Latin hypercube estimate of <rt'thl'>   [K*(kg/kg)]

    ! Local variables
    real( kind = core_rknd ), dimension(nz) :: &  
      rt_mean, &    ! Latin hypercube estimate of rtm         [kg/kg]
      thl_mean, &   ! Latin hypercube estimate of thlm        [K]
      w_mean, &     ! Latin hypercube estimate of wm          [m/s]
      rtp2_zt, &    ! Estimate of <rt'^2> on the zt grid      [(kg/kg)^2]
      thlp2_zt, &   ! Estimate of <thl'^2> on the zt grid     [K^2]
      wprtp_zt, &   ! Estimate of <w'rt'> on the zt grid      [m*(kg/kg)/s]
      wpthlp_zt, &  ! Estimate of <w'thl'> on the zt grid     [m*K/s]
      rtpthlp_zt    ! Estimate of <rt'thl'> on the zt grid    [K*(kg/kg)]


    ! ---- Begin code ----

    rt_mean = compute_sample_mean( nz, n_samples, lh_weights, rt_all_samples )
    thl_mean = compute_sample_mean( nz, n_samples, lh_weights, thl_all_samples )
    w_mean = compute_sample_mean( nz, n_samples, lh_weights, w_all_samples )

    rtp2_zt = compute_sample_variance( nz, n_samples, rt_all_samples, lh_weights, rt_mean )
    thlp2_zt = compute_sample_variance( nz, n_samples, thl_all_samples, lh_weights, thl_mean )
  
    wprtp_zt = compute_sample_covariance( nz, n_samples, lh_weights, &
                   w_all_samples, w_mean, rt_all_samples, rt_mean ) 
    wpthlp_zt = compute_sample_covariance( nz, n_samples, lh_weights, &
                   w_all_samples, w_mean, thl_all_samples, thl_mean )
    rtpthlp_zt = compute_sample_covariance( nz, n_samples, lh_weights, &
                   rt_all_samples, rt_mean, thl_all_samples, thl_mean ) 

    lh_rtp2 = zt2zm( rtp2_zt )
    lh_thlp2 = zt2zm( thlp2_zt )
    lh_wprtp = zt2zm( wprtp_zt )
    lh_wpthlp = zt2zm( wpthlp_zt )
    lh_rtpthlp = zt2zm( rtpthlp_zt )

    return
  end subroutine lh_moments
  !-----------------------------------------------------------------------------
end module estimate_scm_microphys_module
