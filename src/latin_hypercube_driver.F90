! $Id$
!-------------------------------------------------------------------------------
module latin_hypercube_mod

  implicit none
   
  public :: latin_hypercube_driver

  private ! Default scope

  logical, parameter, private :: &
    l_diagnostic_iter_check = .true.

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
               hydromet, correlation_array, hydromet_mc_est, hydromet_vel_est, rcm_mc_est, &
               rvm_mc_est, thlm_mc_est, microphys_sub )

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

    use parameters_model, only: hydromet_dim ! Variable

#ifdef UNRELEASED_CODE
    use permute_height_time_mod, only: & 
      permute_height_time ! Procedure

    use generate_lh_sample_mod, only: & 
      generate_lh_sample ! Procedure

    use estimate_lh_micro_mod, only: & 
      estimate_lh_micro ! Procedure
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
      iLH_Ncp2_zt, &
      iLH_cloud_frac, &
      zt

    use stats_type, only: &
      stat_update_var ! Procedure(s)

    implicit none

    ! External
    intrinsic :: allocated, mod

    ! Interface block
#include "./microphys_interface.inc"

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
      hydromet_mc_est, & ! LH estimate of hydrometeor time tendency          [(units vary)/s]
      hydromet_vel_est   ! LH estimate of hydrometeor sedimentation velocity [m/s]

    ! Output Variables
    real, dimension(nnzp), intent(out) :: &
      rcm_mc_est, & ! LH estimate of time tendency of liquid water mixing ratio    [kg/kg/s]
      rvm_mc_est, & ! LH estimate of time tendency of vapor water mixing ratio     [kg/kg/s]
      thlm_mc_est   ! LH estimate of time tendency of liquid potential temperature [K/s]

    type(pdf_parameter), intent(in) :: pdf_params

    ! Local variables
    integer :: p_matrix(n_micro_calls,d_variables+1)

    double precision, dimension(nnzp,n_micro_calls,(d_variables+1)) :: &
      X_u_all_levs ! Sample drawn from uniform distribution

    double precision, dimension(nnzp,n_micro_calls,d_variables) :: &
      X_nl_all_levs ! Sample that is transformed ultimately to normal-lognormal

    double precision, dimension(nnzp,n_micro_calls) :: &
      rt, thl ! Sample of total water and liquid potential temperature [g/kg],[K]

    real, dimension(nnzp,hydromet_dim) :: &
      lh_hydromet ! Average value of the latin hypercube est. of all hydrometeors [units vary]

    real, dimension(nnzp) :: &
      lh_thlm,    & ! Average value of the latin hypercube est. of theta_l           [K]
      lh_rcm,     & ! Average value of the latin hypercube est. of rc                [kg/kg]
      lh_rvm,     & ! Average value of the latin hypercube est. of rv                [kg/kg]
      lh_wm,      & ! Average value of the latin hypercube est. of vertical velocity [m/s]
      lh_wp2_zt,  & ! Average value of the variance of the LH est. of vert. vel.     [m^2/s^2]
      lh_Ncp2_zt, & ! Average value of the variance of the LH est. of Nc.            [#/kg]
      lh_cloud_frac ! Average value of the latin hypercube est. of cloud fraction    [-]

    ! A true/false flag that determines whether the PDF allows us to construct a sample
    logical, dimension(nnzp) :: l_sample_flag 

    ! Number of random samples before sequence of repeats (normally=10)
    integer :: nt_repeat

    integer :: i_rmd, k

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

      ! Generate LH sample, represented by X_u and X_nl, for level k
      call generate_lh_sample &
           ( n_micro_calls, nt_repeat, d_variables, hydromet_dim, &        ! intent(in)
             p_matrix, cloud_frac(k), pdf_params, k, &                     ! intent(in)
             hydromet(k,:), correlation_array(k,:,:), &                    ! intent(in)
             rt(k,:), thl(k,:), &                                          ! intent(out)
             X_u_all_levs(k,:,:), X_nl_all_levs(k,:,:), l_sample_flag(k) ) ! intent(out)

      ! print *, 'latin_hypercube_sampling: got past lh_sampler'
    end do ! 1..nnzp

    ! Perform LH and analytic microphysical calculations
    call estimate_lh_micro &
         ( dt, nnzp, n_micro_calls, d_variables, &  ! intent(in)
           X_u_all_levs, X_nl_all_levs, &           ! intent(in)
           rt, thl, l_sample_flag, pdf_params, &    ! intent(in)
           thlm, p_in_Pa, exner, rho, &             ! intent(in)
           wm, w_std_dev, dzq, rcm, rvm, &          ! intent(in)
           cloud_frac, hydromet, &                  ! intent(in)
           hydromet_mc_est, hydromet_vel_est, &     ! intent(in)
           rcm_mc_est, rvm_mc_est, thlm_mc_est, &   ! intent(out)
           lh_AKm, AKm, AKstd, AKstd_cld, &        ! intent(out)
           AKm_rcm, AKm_rcc, lh_rcm_avg, &          ! intent(out)
           lh_hydromet, lh_thlm, lh_rcm, lh_rvm, &  ! intent(out)
           lh_wm, lh_Ncp2_zt, lh_wp2_zt, lh_cloud_frac, &  ! intent(out)
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
      call stat_update_var( iLH_cloud_frac, lh_cloud_frac, zt )

    end if ! l_stats_samp

    return

#else
    stop "This code was not compiled with support for Latin Hypercube sampling"

    ! This is simply to avoid a compiler warning
    rcm_mc_est  = -999.999
    rvm_mc_est  = -999.999
    thlm_mc_est = -999.999

#endif /* UNRELEASED_CODE */

  end subroutine latin_hypercube_driver

end module latin_hypercube_mod
