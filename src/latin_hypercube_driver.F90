! $Id$
!-------------------------------------------------------------------------------
module latin_hypercube_mod

  implicit none
   
  public :: latin_hypercube_driver

  private ! Default scope

  integer, parameter :: &
    d_variables     = 5, & ! Number of variables to sample
    n_micro_calls   = 2, & ! Number of calls to microphysics per timestep (normally=2)
    sequence_length = 1     ! nt_repeat/n_micro_call; number of timesteps before sequence repeats.

  ! Number of random samples before sequence of repeats (normally=10)
  integer, parameter :: &
    nt_repeat = n_micro_calls * sequence_length 


  integer, allocatable, dimension(:,:,:) :: & 
    height_time_matrix ! matrix of rand ints

  contains

!-------------------------------------------------------------------------------
  subroutine latin_hypercube_driver &
             ( dt, iter, nnzp, cf, T_in_K, p_in_Pa, exner, &
               rho, pdf_params, wm, w_std_dev, altitudes, rcm, rvm, &
               hydromet, hydromet_mc_est, hydromet_vel_est, rcm_mc_est, &
               rvm_mc_est, thlm_mc_est, microphys_sub )

! Description:
!   Call a microphysics scheme or generate a estimate of Kessler autoconversion
!   using latin hypercube sampling.
! References:
!   None
!-------------------------------------------------------------------------------
    use array_index, only: & 
      iirrainm,    & ! Variables  
      iiNcm,       & 
      iirsnowm,    &
      iirgraupelm, &
      iiricem

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
      AKm_est,  & 
      AKm, & 
      AKstd, & 
      AKstd_cld, & 
      AKm_rcm, & 
      AKm_rcc, & 
      rcm_est

    use stats_variables, only: &
      l_stats_samp, & ! Variables
      iLH_rrainm, &
      iLH_ricem, &
      iLH_rsnowm, &
      iLH_rgraupelm, &
      iLH_thlm, &
      iLH_rcm, &
      iLH_rvm, &
      iLH_wm, &
      iLH_cf, &
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
      iter, & ! Model iteration number
      nnzp    ! Domain dimension

    real, dimension(nnzp), intent(in) :: &
      cf,         & ! Cloud fraction           [%]
      T_in_K,     & ! Temperature              [K]
      p_in_Pa,    & ! Pressure                 [Pa]
      exner,      & ! Exner function           [-]
      rho           ! Density on thermo. grid  [kg/m^3]

    real, dimension(nnzp), intent(in) :: &
      wm, &        ! Mean w                     [m/s]
      w_std_dev, & ! Standard deviation of w    [m/s]
      altitudes    ! Altitudes                  [m]

    real, dimension(nnzp), intent(in) :: &
      rcm, & ! Liquid water mixing ratio        [kg/kg]
      rvm    ! Vapor water mixing ratio         [kg/kg]

    real, dimension(nnzp,hydromet_dim), intent(in) :: &
      hydromet ! Hydrometeor species    [units vary]

    real, dimension(nnzp,hydromet_dim), intent(inout) :: &
      hydromet_mc_est, & ! LH estimate of hydrometeor time tendency          [(units vary)/s]
      hydromet_vel_est   ! LH estimate of hydrometeor sedimentation velocity [m/s]

    real, dimension(nnzp), intent(out) :: &
      rcm_mc_est, & ! LH estimate of time tendency of liquid water mixing ratio    [kg/kg/s]
      rvm_mc_est, & ! LH estimate of time tendency of vapor water mixing ratio     [kg/kg/s]
      thlm_mc_est   ! LH estimate of time tendency of liquid potential temperature [K/s]

    type(pdf_parameter), intent(in) :: pdf_params

    ! Local variables
    integer :: p_matrix(n_micro_calls,d_variables+1)

    ! Sample drawn from uniform distribution
    double precision, dimension(nnzp,n_micro_calls,(d_variables+1)) :: X_u

    ! Sample that is transformed ultimately to normal-lognormal
    double precision, dimension(nnzp,n_micro_calls,d_variables) :: X_nl

    double precision, dimension(nnzp,n_micro_calls) :: &
      rt, thl ! Sample of total water and liquid potential temperature [g/kg],[K]

    real, dimension(nnzp,hydromet_dim) :: &
      lh_hydromet ! Average value of the latin hypercube est. of all hydrometeors [units vary]

    real, dimension(nnzp) :: &
      lh_thlm, & ! Average value of the latin hypercube est. of theta_l           [K]
      lh_rcm,  & ! Average value of the latin hypercube est. of rc                [kg/kg]
      lh_rvm,  & ! Average value of the latin hypercube est. of rv                [kg/kg]
      lh_wm,   & ! Average value of the latin hypercube est. of vertical velocity [m/s]
      lh_cf      ! Average value of the latin hypercube est. of cloud fraction    [%]

    ! A true/false flag that determines whether the PDF allows us to construct a sample
    logical, dimension(nnzp) :: l_sample_flag 

    integer :: i_rmd, k

#ifdef UNRELEASED_CODE
    ! ---- Begin Code ----

    if ( .not. allocated( height_time_matrix ) ) then
      allocate( height_time_matrix(nnzp, nt_repeat, d_variables+1) )
    end if

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
    if ( iirrainm < 1 ) then
      write(fstderr,*) "Latin hypercube sampling is enabled, but there is"// &
      " no rain water mixing ratio being predicted."
      stop
    end if

    do k = 1, nnzp
      ! Choose which rows of LH sample to feed into closure.
      p_matrix(1:n_micro_calls,1:(d_variables+1)) = &
        height_time_matrix(k, n_micro_calls*i_rmd+1:n_micro_calls*i_rmd+n_micro_calls, &
                           1:d_variables+1)

      ! print*, 'latin_hypercube_sampling: got past p_matrix'
      ! Convert from number/kg of air to number/cc
      ! Ncm = dble( hydromet(k,iiNcm) * rho(k) / cm3_per_m3 )

      ! Generate LH sample, represented by X_u and X_nl, for level k
      call generate_lh_sample &
           ( n_micro_calls, nt_repeat, d_variables, p_matrix, & ! intent(in)
             cf(k), pdf_params, k, &                           ! intent(in)
             max( hydromet(k,iiNcm), 1.0 ), hydromet(k,iirrainm), &     ! intent(in)
             rt(k,:), thl(k,:), &       ! intent(out)
             X_u(k,:,:), X_nl(k,:,:), l_sample_flag(k) ) ! intent(out)

      ! print *, 'latin_hypercube_sampling: got past lh_sampler'
    end do ! 1..nnzp

    ! Perform LH and analytic microphysical calculations
    call estimate_lh_micro &
         ( dt, nnzp, n_micro_calls, d_variables, X_u, X_nl, & ! intent(in)
           rt, thl, l_sample_flag, pdf_params, &             ! intent(in)
           T_in_K, p_in_Pa, exner, rho, &                ! intent(in)
           wm, w_std_dev, altitudes, rcm, rvm, &         ! intent(in)
           cf, hydromet, &                               ! intent(in)
           hydromet_mc_est, hydromet_vel_est, &          ! intent(in)
           rcm_mc_est, rvm_mc_est, thlm_mc_est, &        ! intent(out)
           AKm_est, AKm, AKstd, AKstd_cld, &             ! intent(out)
           AKm_rcm, AKm_rcc, rcm_est, &                  ! intent(out)
           lh_hydromet, lh_thlm, lh_rcm, lh_rvm, &       ! intent(out)
           lh_wm, lh_cf, &                               ! intent(out)
           microphys_sub )  ! Procedure

    ! print*, 'latin_hypercube_driver: AKm=', AKm
    ! print*, 'latin_hypercube_driver: AKm_est=', AKm_est

    if ( l_stats_samp ) then
      if ( iirrainm > 0 ) then
        call stat_update_var( iLH_rrainm, lh_hydromet(:,iirrainm), zt )
      end if
      if ( iiricem > 0 ) then
        call stat_update_var( iLH_ricem, lh_hydromet(:,iiricem), zt )
      end if
      if ( iirsnowm > 0 ) then
        call stat_update_var( iLH_rsnowm, lh_hydromet(:,iirsnowm), zt )
      end if
      if ( iirgraupelm > 0 ) then
        call stat_update_var( iLH_rgraupelm, lh_hydromet(:,iirgraupelm), zt )
      end if

      call stat_update_var( iLH_thlm, lh_thlm, zt )
      call stat_update_var( iLH_rcm, lh_rcm, zt )
      call stat_update_var( iLH_rvm, lh_rvm, zt )
      call stat_update_var( iLH_wm, lh_wm, zt )
      call stat_update_var( iLH_cf, lh_cf, zt )
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
