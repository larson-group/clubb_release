! $Id$
module estimate_scm_microphys_module

  implicit none

  public :: est_single_column_tndcy, copy_X_nl_into_hydromet

  private ! Default scope

  contains

!-------------------------------------------------------------------------------
  subroutine est_single_column_tndcy &
             ( dt, nz, n_micro_calls, d_variables, &
               k_lh_start, LH_rt, LH_thl, &
               X_nl_all_levs, LH_sample_point_weights, &
               p_in_Pa, exner, rho, cloud_frac, w_std_dev, &
               dzq, pdf_params, hydromet,  &
               lh_rvm_mc, lh_rcm_mc, lh_hydromet_mc, &
               lh_hydromet_vel, lh_thlm_mc, &
               microphys_sub )
! Description:
!   Estimate the tendency of a microphysics scheme via latin hypercube sampling
!
! References:
!   None
!-------------------------------------------------------------------------------

    use constants_clubb, only:  &
      fstderr, &  ! Constant(s)
      zero_threshold, &
      rc_tol, &
      cm3_per_m3

    use parameters_model, only: &
      hydromet_dim ! Variable

!   use parameters_microphys, only: &
!     Ncm_initial

    use parameters_microphys, only: &
      l_lh_cloud_weighted_sampling

    use corr_matrix_module, only: &
      iiLH_s_mellor, &
      iiLH_w

    use pdf_parameter_module, only: &
      pdf_parameter ! Type

    use error_code, only: &
      clubb_at_least_debug_level ! Procedure

    use clubb_precision, only: &
      dp, & ! double precision
      core_rknd, &
      time_precision

    implicit none

    ! External
#include "../microphys_interface.inc"

    intrinsic :: real, dble

    ! Constant parameters
    logical, parameter :: &
      l_stats_samp                 = .false., &
      l_local_kk                   = .true., &
      l_latin_hypercube            = .true.

    logical, parameter :: &
      l_check_lh_cloud_weighting = .true. ! Verify every other sample point is out of cloud

    ! Input Variables
    real( kind = time_precision ), intent(in) :: &
      dt ! Model timestep       [s]

    integer, intent(in) :: &
      nz,          & ! Number of vertical levels
      n_micro_calls, & ! Number of calls to microphysics (normally=2)
      d_variables,   & ! Number of variates (normally=5) 
      k_lh_start       ! Starting level for computing arbitrary overlap

    real( kind = core_rknd ), dimension(nz,n_micro_calls), intent(in) :: &
      LH_rt, & ! n_micro_calls values of total water mixing ratio     [kg/kg]
      LH_thl   ! n_micro_calls values of liquid potential temperature [K]

    real( kind = dp ), target, dimension(nz,n_micro_calls,d_variables), intent(in) :: &
      X_nl_all_levs ! Sample that is transformed ultimately to normal-lognormal

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      p_in_Pa,    & ! Pressure                 [Pa]
      exner,      & ! Exner function           [-]
      rho,        & ! Density on thermo. grid  [kg/m^3]
      cloud_frac, & ! Cloud fraction           [-]
      w_std_dev,  & ! Standard deviation of w    [m/s]
      dzq           ! Difference in height per gridbox   [m]

    type(pdf_parameter), dimension(nz), intent(in) :: pdf_params

    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(in) :: &
      hydromet ! Hydrometeor species    [units vary]

    real( kind = core_rknd ), dimension(n_micro_calls), intent(in) :: &
       LH_sample_point_weights ! Weight for cloud weighted sampling

    ! Output Variables

    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(inout) :: &
      lh_hydromet_mc, & ! LH estimate of hydrometeor time tendency          [(units vary)/s]
      lh_hydromet_vel   ! LH estimate of hydrometeor sedimentation velocity [m/s]

    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      lh_rcm_mc, & ! LH estimate of time tendency of liquid water mixing ratio    [kg/kg/s]
      lh_rvm_mc, & ! LH estimate of time tendency of vapor water mixing ratio     [kg/kg/s]
      lh_thlm_mc   ! LH estimate of time tendency of liquid potential temperature [K/s]

    ! Local Variables
    real( kind = dp ), dimension(nz,hydromet_dim) :: &
      lh_hydromet_mc_sum, & ! LH est of hydrometeor time tendency          [(units vary)/s]
      lh_hydromet_vel_sum   ! LH est of hydrometeor sedimentation velocity [m/s]

    real( kind = dp ), dimension(nz) :: &
      lh_rcm_mc_sum,  & ! LH est of time tendency of liquid water mixing ratio    [kg/kg/s]
      lh_rvm_mc_sum,  & ! LH est of time tendency of vapor water mixing ratio     [kg/kg/s]
      lh_thlm_mc_sum    ! LH est of time tendency of liquid potential temperature [K/s]

    real( kind = core_rknd ), dimension(nz,hydromet_dim) :: &
      hydromet_columns ! Hydrometeor species    [units vary]

    real( kind = core_rknd ), dimension(nz) :: &
      s_mellor_column    ! 's' (Mellor 1977)            [kg/kg]

    real( kind = core_rknd ), dimension(nz) :: &
      rv_column,  & ! Vapor water                  [kg/kg]
      thl_column, & ! Liquid potential temperature [K]
      rc_column,  & ! Liquid water                [kg/kg]
      w_column      ! Vertical velocity           [m/s]

    real( kind = dp ), pointer, dimension(:,:) :: &
      s_mellor_all_points,  & ! n_micro_calls values of 's' (Mellor 1977)      [kg/kg]
      w_all_points            ! n_micro_calls values of vertical velocity      [m/s]

    real( kind = core_rknd ), dimension(nz) :: &
      lh_wprtp_mc_tndcy,   & ! LH micro. tendency for <w'rt'>   [m*(kg/kg)/s^2]
      lh_wpthlp_mc_tndcy,  & ! LH micro. tendency for <w'thl'>  [m*K/s^2]
      lh_rtp2_mc_tndcy,    & ! LH micro. tendency for <rt'^2>   [(kg/kg)^2/s]
      lh_thlp2_mc_tndcy,   & ! LH micro. tendency for <thl'^2>  [K^2/s]
      lh_rtpthlp_mc_tndcy    ! LH micro. tendency for <rt'thl'> [K*(kg/kg)/s]

    integer :: ivar, k, sample

    integer :: &
      in_cloud_points, &
      out_of_cloud_points

    ! ---- Begin Code ----

    ! Mellor's 's' is hardwired elsewhere to be the first column
    s_mellor_all_points => X_nl_all_levs(:,:,iiLH_s_mellor)
    w_all_points        => X_nl_all_levs(:,:,iiLH_w)

    ! Assertion check
    if ( clubb_at_least_debug_level( 2 ) ) then
      if ( l_check_lh_cloud_weighting .and. l_lh_cloud_weighted_sampling .and. &
           all( LH_sample_point_weights(:) /= 1.0_core_rknd ) ) then 
              ! The 1.0 indicates cloud_frac is > 0.5
        ! Verify every other sample point is out of cloud if we're doing
        ! cloud weighted sampling
        in_cloud_points     = 0
        out_of_cloud_points = 0
        do sample = 1, n_micro_calls, 1
          if ( s_mellor_all_points(k_lh_start,sample) > 0._core_rknd) then
            in_cloud_points = in_cloud_points + 1
          else if ( s_mellor_all_points(k_lh_start,sample) <= 0._core_rknd) then
            out_of_cloud_points = out_of_cloud_points + 1
          end if
        end do ! 1..n_micro_calls
        if ( in_cloud_points /= out_of_cloud_points ) then
          write(fstderr,*) "In est_single_column_tndcy:"
          write(fstderr,*) "The cloudy sample points do not equal the out of cloud points"
          write(fstderr,*) "in_cloud_points =", in_cloud_points
          write(fstderr,*) "out_of_cloud_points =", out_of_cloud_points
        end if
      end if ! l_check_lh_cloud_weighting .and. l_lh_cloud_weighted_sampling
    end if ! clubb_at_least_debug_level 2

    lh_hydromet_vel(:,:) = 0._core_rknd

    ! Initialize microphysical tendencies for each mixture component
    lh_hydromet_mc_sum(:,:) = 0._dp

    lh_hydromet_vel_sum(:,:) = 0._dp

    lh_rcm_mc_sum(:) = 0._dp

    lh_rvm_mc_sum(:) = 0._dp

    lh_thlm_mc_sum(:) = 0._dp

    do sample = 1, n_micro_calls

      s_mellor_column = real( s_mellor_all_points(:,sample), kind = core_rknd )

      where( s_mellor_all_points(:,sample) > 0.0_core_rknd )
        rc_column = real( s_mellor_all_points(:,sample), kind = core_rknd )
      else where
        rc_column = 0.0_core_rknd
      end where

      w_column   = real( w_all_points(:,sample), kind = core_rknd )
      rv_column  = real( LH_rt(:,sample), kind = core_rknd ) - rc_column
      ! Verify total water isn't negative
      if ( any( rv_column < 0._core_rknd) ) then
        if ( clubb_at_least_debug_level( 1 ) ) then
          write(fstderr,*) "rv negative, LH sample number = ", sample
          write(fstderr,'(a3,3a20)') "k", "rt", "rv", "rc"
          do k = 1, nz
            if ( rv_column(k) < 0._core_rknd) then
              write(6,'(i3,3g20.7)')  k, LH_rt(k,sample), rv_column(k), &
                rc_column(k)
            end if
          end do
        end if ! clubb_at_least_debug_level( 1 )
        write(fstderr,*) "Applying non-conservative hard clipping to rv sample."
        where ( rv_column < 0._core_rknd) rv_column = zero_threshold
      end if ! Some rv_column element < 0

      thl_column = real( LH_thl(:,sample), kind = core_rknd )

      call copy_X_nl_into_hydromet( nz, d_variables, 1, & ! In
                                    X_nl_all_levs(:,sample,:), & ! In
                                    hydromet, & ! In
                                    hydromet_columns )

      ! Call the microphysics scheme to obtain a sample point
      call microphys_sub &
           ( dt, nz, l_stats_samp, l_local_kk, & ! In
             l_latin_hypercube, thl_column, w_column, p_in_Pa, & ! In
             exner, rho, cloud_frac, pdf_params, w_std_dev, & ! In
             dzq, rc_column, s_mellor_column, rv_column, hydromet_columns, & ! In
             lh_hydromet_mc, lh_hydromet_vel, & ! Out
             lh_rcm_mc, lh_rvm_mc, lh_thlm_mc, & ! Out
             lh_wprtp_mc_tndcy, lh_wpthlp_mc_tndcy, & ! Out
             lh_rtp2_mc_tndcy, lh_thlp2_mc_tndcy, lh_rtpthlp_mc_tndcy ) ! Out

      if ( l_lh_cloud_weighted_sampling ) then
        ! Weight the output results depending on whether we're calling the
        ! microphysics on clear or cloudy air
        lh_hydromet_vel(:,:) = lh_hydromet_vel(:,:) * LH_sample_point_weights(sample)
        lh_hydromet_mc(:,:) = lh_hydromet_mc(:,:) * LH_sample_point_weights(sample)
        lh_rcm_mc(:) = lh_rcm_mc(:) * LH_sample_point_weights(sample)
        lh_rvm_mc(:) = lh_rvm_mc(:) * LH_sample_point_weights(sample)
        lh_thlm_mc(:) = lh_thlm_mc(:) * LH_sample_point_weights(sample)
      end if

      do ivar = 1, hydromet_dim
        lh_hydromet_vel_sum(:,ivar) = lh_hydromet_vel_sum(:,ivar) + lh_hydromet_vel(:,ivar)
        lh_hydromet_mc_sum(:,ivar) = lh_hydromet_mc_sum(:,ivar) + lh_hydromet_mc(:,ivar)
      end do

      lh_rcm_mc_sum(:) = lh_rcm_mc_sum(:) + lh_rcm_mc(:)
      lh_rvm_mc_sum(:) = lh_rvm_mc_sum(:) + lh_rvm_mc(:)
      lh_thlm_mc_sum(:) = lh_thlm_mc_sum(:) + lh_thlm_mc(:)

      ! Loop to get new sample
    end do ! sample = 1, n_micro_calls


    ! Grid box average.
    forall( ivar = 1:hydromet_dim )
      lh_hydromet_vel(:,ivar) = lh_hydromet_vel_sum(:,ivar) / real( n_micro_calls, kind=core_rknd )
      lh_hydromet_mc(:,ivar) = lh_hydromet_mc_sum(:,ivar) / real( n_micro_calls, kind=core_rknd )
    end forall

    lh_rcm_mc = lh_rcm_mc_sum / real( n_micro_calls, kind=core_rknd )
    lh_rvm_mc = lh_rvm_mc_sum / real( n_micro_calls, kind=core_rknd )
    lh_thlm_mc = lh_thlm_mc_sum / real( n_micro_calls, kind=core_rknd )

    return
  end subroutine est_single_column_tndcy

  !-----------------------------------------------------------------------------
  subroutine copy_X_nl_into_hydromet( nz, d_variables, n_micro_calls, &
                                      X_nl_all_levs, &
                                      hydromet, &
                                      hydromet_all_points )

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
      iiNgraupelm, &
      iiNcm

    use corr_matrix_module, only: &
      iiLH_rrain, &
      iiLH_rsnow, &
      iiLH_rice, &
      iiLH_rgraupel, &
      iiLH_Nr, &
      iiLH_Nsnow, &
      iiLH_Ngraupel, &
      iiLH_Nc, &
      iiLH_Ni

    use clubb_precision, only: &
      dp, & ! double precision
      core_rknd

    implicit none

    integer, intent(in) :: &
      nz,          & ! Number of vertical levels
      d_variables,   & ! Number of variates (normally=5) 
      n_micro_calls    ! Number of calls to microphysics (normally=2)

    real( kind = dp ), target, dimension(nz,n_micro_calls,d_variables), intent(in) :: &
      X_nl_all_levs ! Sample that is transformed ultimately to normal-lognormal

    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(in) :: &
      hydromet ! Hydrometeor species    [units vary]

    real( kind = core_rknd ), dimension(nz,n_micro_calls,hydromet_dim), intent(out) :: &
      hydromet_all_points ! Hydrometeor species    [units vary]

    integer :: sample, ivar

    do sample = 1, n_micro_calls
      ! Copy the sample points into the temporary arrays
      do ivar = 1, hydromet_dim, 1
        if ( ivar == iirrainm .and. iiLH_rrain > 0 ) then
          ! Use a sampled value of rain water mixing ratio
          hydromet_all_points(:,sample,ivar) = &
            real( X_nl_all_levs(:,sample,iiLH_rrain), kind = core_rknd )

        else if ( ivar == iirsnowm .and. iiLH_rsnow > 0 ) then
          ! Use a sampled value of rain water mixing ratio
          hydromet_all_points(:,sample,ivar) = &
            real( X_nl_all_levs(:,sample,iiLH_rsnow), kind = core_rknd )

        else if ( ivar == iiricem .and. iiLH_rice > 0 ) then
          ! Use a sampled value of rain water mixing ratio
          hydromet_all_points(:,sample,ivar) = &
            real( X_nl_all_levs(:,sample,iiLH_rice), kind = core_rknd )

        else if ( ivar == iirgraupelm .and. iiLH_rgraupel > 0 ) then
          ! Use a sampled value of rain water mixing ratio
          hydromet_all_points(:,sample,ivar) = &
            real( X_nl_all_levs(:,sample,iiLH_rgraupel), kind = core_rknd )

        else if ( ivar == iiNcm .and. iiLH_Nc > 0 ) then
          ! Kluge for when we don't have correlations between Nc, other variables
!         hydromet_all_points(:,iiNcm) = Ncm_initial * cm3_per_m3 / rho
          ! Use a sampled value of cloud droplet number concentration
          hydromet_all_points(:,sample,ivar) = &
            real( X_nl_all_levs(:,sample,iiLH_Nc), kind = core_rknd )

        else if ( ivar == iiNrm .and. iiLH_Nr > 0 ) then
          ! Use a sampled value of rain droplet number concentration
          hydromet_all_points(:,sample,ivar) = &
            real( X_nl_all_levs(:,sample,iiLH_Nr), kind = core_rknd )

        else if ( ivar == iiNsnowm .and. iiLH_Nsnow > 0 ) then
          ! Use a sampled value of rain droplet number concentration
          hydromet_all_points(:,sample,ivar) = &
            real( X_nl_all_levs(:,sample,iiLH_Nsnow), kind = core_rknd )

        else if ( ivar == iiNgraupelm .and. iiLH_Ngraupel > 0 ) then
          ! Use a sampled value of rain droplet number concentration
          hydromet_all_points(:,sample,ivar) = &
            real( X_nl_all_levs(:,sample,iiLH_Ngraupel), kind = core_rknd )

        else if ( ivar == iiNim .and. iiLH_Ni > 0 ) then
          ! Use a sampled value of rain droplet number concentration
          hydromet_all_points(:,sample,ivar) = &
            real( X_nl_all_levs(:,sample,iiLH_Ni), kind = core_rknd )

        else ! Use the mean field, rather than a sample point
          ! This is the case for ice phase fields in the Morrison microphysics
          ! currently -dschanen 23 March 2010
          hydromet_all_points(:,sample,ivar) = hydromet(:,ivar)

        end if
      end do ! 1..hydromet_dim
    end do ! 1..n_micro_calls

    return
  end subroutine copy_X_nl_into_hydromet
  !-----------------------------------------------------------------------------
end module estimate_scm_microphys_module
