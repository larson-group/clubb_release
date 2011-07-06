! $Id$
module estimate_scm_microphys_module

  implicit none

  public :: est_single_column_tndcy, copy_X_nl_into_hydromet

  private ! Default scope

  contains

!-------------------------------------------------------------------------------
  subroutine est_single_column_tndcy &
             ( dt, nnzp, n_micro_calls, d_variables, &
               k_lh_start, LH_rt, LH_thl, &
               X_nl_all_levs, LH_sample_point_weights, &
               p_in_Pa, exner, rho, w_std_dev, &
               dzq, pdf_params, hydromet,  &
               X_mixt_comp_all_levs, &
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
      iiLH_rrain, &
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

    use math_utilities, only: &
      compute_sample_variance, & ! Procedure
      compute_sample_mean

    use pdf_parameter_module, only: &
      pdf_parameter ! Type

    use error_code, only: &
      clubb_at_least_debug_level ! Procedure

    implicit none

    ! External
#include "../microphys_interface.inc"

    intrinsic :: real, dble

    ! Constant parameters
    logical, parameter :: &
      l_compute_diagnostic_average = .true., &
      l_stats_samp                 = .false., &
      l_local_kk                   = .true., &
      l_latin_hypercube            = .true.

    logical, parameter :: &
      l_check_lh_cloud_weighting = .true. ! Verify every other sample point is out of cloud

    ! Input Variables
    real, intent(in) :: &
      dt ! Model timestep       [s]

    integer, intent(in) :: &
      nnzp,          & ! Number of vertical levels
      n_micro_calls, & ! Number of calls to microphysics (normally=2)
      d_variables,   & ! Number of variates (normally=5) 
      k_lh_start       ! Starting level for computing arbitrary overlap

    real, dimension(nnzp,n_micro_calls), intent(in) :: &
      LH_rt, & ! n_micro_calls values of total water mixing ratio     [kg/kg]
      LH_thl   ! n_micro_calls values of liquid potential temperature [K]

    double precision, target, dimension(nnzp,n_micro_calls,d_variables), intent(in) :: &
      X_nl_all_levs ! Sample that is transformed ultimately to normal-lognormal

    real, dimension(nnzp), intent(in) :: &
      p_in_Pa,    & ! Pressure                 [Pa]
      exner,      & ! Exner function           [-]
      rho           ! Density on thermo. grid  [kg/m^3]

    real, dimension(nnzp), intent(in) :: &
      w_std_dev, & ! Standard deviation of w    [m/s]
      dzq          ! Difference in height per gridbox   [m]

    type(pdf_parameter), dimension(nnzp), intent(in) :: pdf_params

    real, dimension(nnzp,hydromet_dim), intent(in) :: &
      hydromet ! Hydrometeor species    [units vary]

    integer, dimension(nnzp,n_micro_calls), intent(in) :: &
      X_mixt_comp_all_levs ! Whether we're in the first or second mixture component

    real, dimension(n_micro_calls), intent(in) :: &
       LH_sample_point_weights ! Weight for cloud weighted sampling

    ! Output Variables

    real, dimension(nnzp,hydromet_dim), intent(inout) :: &
      lh_hydromet_mc, & ! LH estimate of hydrometeor time tendency          [(units vary)/s]
      lh_hydromet_vel   ! LH estimate of hydrometeor sedimentation velocity [m/s]

    real, dimension(nnzp), intent(out) :: &
      lh_rcm_mc, & ! LH estimate of time tendency of liquid water mixing ratio    [kg/kg/s]
      lh_rvm_mc, & ! LH estimate of time tendency of vapor water mixing ratio     [kg/kg/s]
      lh_thlm_mc   ! LH estimate of time tendency of liquid potential temperature [K/s]

    ! Local Variables
    double precision, dimension(nnzp,hydromet_dim) :: &
      lh_hydromet_mc_sum, & ! LH est of hydrometeor time tendency          [(units vary)/s]
      lh_hydromet_vel_sum   ! LH est of hydrometeor sedimentation velocity [m/s]

    double precision, dimension(nnzp) :: &
      lh_rcm_mc_sum,  & ! LH est of time tendency of liquid water mixing ratio    [kg/kg/s]
      lh_rvm_mc_sum,  & ! LH est of time tendency of vapor water mixing ratio     [kg/kg/s]
      lh_thlm_mc_sum    ! LH est of time tendency of liquid potential temperature [K/s]

    real, dimension(nnzp,hydromet_dim) :: &
      hydromet_columns ! Hydrometeor species    [units vary]

    real, dimension(nnzp) :: &
      s_mellor_column    ! 's' (Mellor 1977)            [kg/kg]

    real, dimension(nnzp) :: &
      rv_column,  & ! Vapor water                  [kg/kg]
      thl_column, & ! Liquid potential temperature [K]
      rc_column,  & ! Liquid water                [kg/kg]
      w_column      ! Vertical velocity           [m/s]

    integer, dimension(nnzp) :: n1, n2, zero

    double precision, pointer, dimension(:,:) :: &
      s_mellor_all_points,  & ! n_micro_calls values of 's' (Mellor 1977)      [kg/kg]
      w_all_points            ! n_micro_calls values of vertical velocity      [m/s]

    integer :: ivar, k, sample

    integer :: &
      in_cloud_points, &
      out_of_cloud_points

    logical :: l_error

    ! ---- Begin Code ----

    ! Mellor's 's' is hardwired elsewhere to be the first column
    s_mellor_all_points => X_nl_all_levs(:,:,iiLH_s_mellor)
    w_all_points        => X_nl_all_levs(:,:,iiLH_w)

    ! Assertion check
    if ( clubb_at_least_debug_level( 2 ) ) then
      if ( l_check_lh_cloud_weighting .and. l_lh_cloud_weighted_sampling .and. &
           all( LH_sample_point_weights(:) /= 1.0 ) ) then ! The 1.0 indicates cloud_frac is > 0.5
        ! Verify every other sample point is out of cloud if we're doing
        ! cloud weighted sampling
        in_cloud_points     = 0
        out_of_cloud_points = 0
        do sample = 1, n_micro_calls, 1
          if ( s_mellor_all_points(k_lh_start,sample) > 0. ) then
            in_cloud_points = in_cloud_points + 1
          else if ( s_mellor_all_points(k_lh_start,sample) <= 0. ) then
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

    lh_hydromet_vel(:,:) = 0.

    ! Initialize microphysical tendencies for each mixture component
    lh_hydromet_mc_sum(:,:) = 0.d0

    lh_hydromet_vel_sum(:,:) = 0.d0

    lh_rcm_mc_sum(:) = 0.d0

    lh_rvm_mc_sum(:) = 0.d0

    lh_thlm_mc_sum(:) = 0.d0

    do sample = 1, n_micro_calls

      s_mellor_column = real( s_mellor_all_points(:,sample) )

      where( s_mellor_all_points(:,sample) > 0.0 )
        rc_column = real( s_mellor_all_points(:,sample) )
      else where
        rc_column = 0.0
      end where

      w_column   = real( w_all_points(:,sample) )
      rv_column  = real( LH_rt(:,sample) ) - rc_column
      ! Verify total water isn't negative
      if ( any( rv_column < 0. ) ) then
        if ( clubb_at_least_debug_level( 1 ) ) then
          write(fstderr,*) "rv negative, LH sample number = ", sample
          write(fstderr,'(a3,3a20)') "k", "rt", "rv", "rc"
          do k = 1, nnzp
            if ( rv_column(k) < 0. ) then
              write(6,'(i3,3g20.7)')  k, LH_rt(k,sample), rv_column(k), &
                rc_column(k)
            end if
          end do
        end if ! clubb_at_least_debug_level( 1 )
        write(fstderr,*) "Applying non-conservative hard clipping to rv sample."
        where ( rv_column < 0. ) rv_column = zero_threshold
      end if ! Some rv_column element < 0

      thl_column = real( LH_thl(:,sample) )

      call copy_X_nl_into_hydromet( nnzp, d_variables, 1, & ! In
                                    X_nl_all_levs(:,sample,:), & ! In
                                    hydromet, & ! In
                                    hydromet_columns )

      ! Call the microphysics scheme to obtain a sample point
      call microphys_sub &
           ( dt, nnzp, l_stats_samp, l_local_kk, l_latin_hypercube, & ! In
             thl_column, p_in_Pa, exner, rho, pdf_params, & ! In
             w_column, w_std_dev, dzq, & ! In
             rc_column, s_mellor_column, & ! In
             rv_column, hydromet_columns,  & ! In
             lh_hydromet_mc, lh_hydromet_vel, lh_rcm_mc, lh_rvm_mc, lh_thlm_mc ) ! Out

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
      lh_hydromet_vel(:,ivar) = real( lh_hydromet_vel_sum(:,ivar) ) / real( n_micro_calls )
      lh_hydromet_mc(:,ivar) = real( lh_hydromet_mc_sum(:,ivar) ) / real( n_micro_calls )
    end forall

    lh_rcm_mc = real( lh_rcm_mc_sum  ) / real( n_micro_calls )
    lh_rvm_mc = real( lh_rvm_mc_sum ) / real( n_micro_calls )
    lh_thlm_mc = real( lh_thlm_mc_sum ) / real( n_micro_calls )

    return
  end subroutine est_single_column_tndcy

  !-----------------------------------------------------------------------------
  subroutine copy_X_nl_into_hydromet( nnzp, d_variables, n_micro_calls, &
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

    use latin_hypercube_arrays, only: &
      iiLH_rrain, &
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

    implicit none

    integer, intent(in) :: &
      nnzp,          & ! Number of vertical levels
      d_variables,   & ! Number of variates (normally=5) 
      n_micro_calls    ! Number of calls to microphysics (normally=2)

    double precision, target, dimension(nnzp,n_micro_calls,d_variables), intent(in) :: &
      X_nl_all_levs ! Sample that is transformed ultimately to normal-lognormal

    real, dimension(nnzp,hydromet_dim), intent(in) :: &
      hydromet ! Hydrometeor species    [units vary]

    real, dimension(nnzp,n_micro_calls,hydromet_dim), intent(out) :: &
      hydromet_all_points ! Hydrometeor species    [units vary]

    integer :: sample, ivar

    do sample = 1, n_micro_calls
      ! Copy the sample points into the temporary arrays
      do ivar = 1, hydromet_dim, 1
        if ( ivar == iirrainm .and. iiLH_rrain > 0 ) then
          ! Use a sampled value of rain water mixing ratio
          hydromet_all_points(:,sample,ivar) = real( X_nl_all_levs(:,sample,iiLH_rrain) )

        else if ( ivar == iirsnowm .and. iiLH_rsnow > 0 ) then
          ! Use a sampled value of rain water mixing ratio
          hydromet_all_points(:,sample,ivar) = real( X_nl_all_levs(:,sample,iiLH_rsnow) )

        else if ( ivar == iiricem .and. iiLH_rice > 0 ) then
          ! Use a sampled value of rain water mixing ratio
          hydromet_all_points(:,sample,ivar) = real( X_nl_all_levs(:,sample,iiLH_rice) )

        else if ( ivar == iirgraupelm .and. iiLH_rgraupel > 0 ) then
          ! Use a sampled value of rain water mixing ratio
          hydromet_all_points(:,sample,ivar) = real( X_nl_all_levs(:,sample,iiLH_rgraupel) )

        else if ( ivar == iiNcm .and. iiLH_Nc > 0 ) then
          ! Kluge for when we don't have correlations between Nc, other variables
!         hydromet_all_points(:,iiNcm) = Ncm_initial * cm3_per_m3 / rho
          ! Use a sampled value of cloud droplet number concentration
          hydromet_all_points(:,sample,ivar) = real( X_nl_all_levs(:,sample,iiLH_Nc) )

        else if ( ivar == iiNrm .and. iiLH_Nr > 0 ) then
          ! Use a sampled value of rain droplet number concentration
          hydromet_all_points(:,sample,ivar) = real( X_nl_all_levs(:,sample,iiLH_Nr) )

        else if ( ivar == iiNsnowm .and. iiLH_Nsnow > 0 ) then
          ! Use a sampled value of rain droplet number concentration
          hydromet_all_points(:,sample,ivar) = real( X_nl_all_levs(:,sample,iiLH_Nsnow) )

        else if ( ivar == iiNgraupelm .and. iiLH_Ngraupel > 0 ) then
          ! Use a sampled value of rain droplet number concentration
          hydromet_all_points(:,sample,ivar) = real( X_nl_all_levs(:,sample,iiLH_Ngraupel) )

        else if ( ivar == iiNim .and. iiLH_Ni > 0 ) then
          ! Use a sampled value of rain droplet number concentration
          hydromet_all_points(:,sample,ivar) = real( X_nl_all_levs(:,sample,iiLH_Ni) )

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
