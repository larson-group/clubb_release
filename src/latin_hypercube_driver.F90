! $Id$
!-------------------------------------------------------------------------------
module latin_hypercube_mod

  implicit none
   
  public :: latin_hypercube_driver, latin_hypercube_2D_output

  private ! Default scope

  logical, parameter, private :: &
    l_diagnostic_iter_check = .true., &
    l_output_2D_samples = .true.

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

    use parameters_model, only: hydromet_dim ! Variable

#ifdef UNRELEASED_CODE
    use permute_height_time_mod, only: & 
      permute_height_time ! Procedure

    use generate_lh_sample_mod, only: & 
      generate_lh_sample ! Procedure

    use estimate_lh_micro_mod, only: & 
      estimate_lh_micro ! Procedure

    use output_2D_samples_mod, only: &
      output_2D_samples_file ! Procedure(s)
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

    double precision, dimension(nnzp,n_micro_calls,(d_variables+1)) :: &
      X_u_all_levs ! Sample drawn from uniform distribution

    double precision, dimension(nnzp,n_micro_calls,d_variables) :: &
      X_nl_all_levs ! Sample that is transformed ultimately to normal-lognormal

    double precision, dimension(nnzp,n_micro_calls) :: &
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
             p_matrix, cloud_frac(k), wm(k), rcm(k)+rvm(k), thlm(k), pdf_params, k, & ! intent(in)
             hydromet(k,:), correlation_array(k,:,:), &                    ! intent(in)
             LH_rt(k,:), LH_thl(k,:), &                                    ! intent(out)
             X_u_all_levs(k,:,:), X_nl_all_levs(k,:,:) )                   ! intent(out)

      ! print *, 'latin_hypercube_sampling: got past lh_sampler'
    end do ! 1..nnzp

    if ( l_output_2D_samples ) then
      call output_2D_samples_file( nnzp, n_micro_calls, d_variables, &
                                   X_nl_all_levs )
    end if

    ! Perform LH and analytic microphysical calculations
    call estimate_lh_micro &
         ( dt, nnzp, n_micro_calls, d_variables, &  ! intent(in)
           X_u_all_levs, X_nl_all_levs, &           ! intent(in)
           LH_rt, LH_thl, pdf_params, &             ! intent(in)
           p_in_Pa, exner, rho, &                   ! intent(in)
           rcm, w_std_dev, dzq, &                   ! intent(in)
           cloud_frac, hydromet, &                  ! intent(in)
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
      hydromet_dim

    use parameters_microphys, only: &
      LH_microphys_calls

    use stats_precision, only: &
      time_precision ! Constant

#ifdef UNRELEASED_CODE
    use output_2D_samples_mod, only: &
      open_2D_samples_file ! Procedure
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

    if ( .not. l_output_2D_samples ) return

    allocate( variable_names(hydromet_dim+3), variable_descriptions(hydromet_dim+3), &
              variable_units(hydromet_dim+3) )

    variable_names(1)        = "s_mellor"
    variable_descriptions(1) = "The variable 's' from Mellor 1977"
    variable_units(1)        = "kg/kg"

    variable_names(2)        = "t_mellor"
    variable_descriptions(2) = "The variable 't' from Mellor 1977"
    variable_units(2)        = "kg/kg"

    variable_names(3)        = "w"
    variable_descriptions(3) = "Vertical velocity"
    variable_units(3)        = "m/s"

    i = 3

    if ( iiLH_Nr > 0 ) then
      i = i + 1
      variable_names(i)        = "Nr"
      variable_descriptions(i) = "Rain droplet number concentration"
      variable_units(i)        = "count/kg"
    end if
    if ( iiLH_Nc > 0 ) then
      i = i + 1
      variable_names(i)        = "Nc"
      variable_descriptions(i) = "Cloud droplet number concentration"
      variable_units(i)        = "count/kg"
    end if
    if ( iiLH_rrain > 0 ) then
      i = i + 1
      variable_names(i)        = "rrain"
      variable_descriptions(i) = "Rain water mixing ratio"
      variable_units(i)        = "kg/kg"
    end if

#ifdef UNRELEASED_CODE
    call open_2D_samples_file( nnzp, LH_microphys_calls, hydromet_dim+3, &
                               fname_prefix, fdir, &
                               time_initial, stats_tout, zt, variable_names, &
                               variable_descriptions, variable_units )
    return
#else
    stop "This code was not compiled with support for Latin Hypercube sampling"
#endif 

  end subroutine latin_hypercube_2D_output

end module latin_hypercube_mod
