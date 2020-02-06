!-------------------------------------------------------------------------------
! $Id$ 
!===============================================================================
module latin_hypercube_driver_module

  use clubb_precision, only: &
    core_rknd           ! Constant

  implicit none

  ! Constant Parameters
  logical, parameter, private :: &
    l_output_2D_lognormal_dist   = .false., & ! Output a 2D netCDF file of the lognormal variates
    l_output_2D_uniform_dist     = .false.    ! Output a 2D netCDF file of the uniform distribution

  private ! Default scope

#ifdef SILHS
  public :: latin_hypercube_2D_output, &
    latin_hypercube_2D_close, stats_accumulate_lh, generate_silhs_sample, &
    copy_X_nl_into_hydromet_all_pts, clip_transform_silhs_output

  private :: stats_accumulate_uniform_lh

  contains

!-------------------------------------------------------------------------------
  subroutine generate_silhs_sample &
             ( iter, pdf_dim, num_samples, sequence_length, nz, & ! intent(in)
               l_calc_weights_all_levs_itime, &                   ! intent(in)
               pdf_params, delta_zm, rcm, Lscale, &               ! intent(in)
               rho_ds_zt, mu1, mu2, sigma1, sigma2, &             ! intent(in)
               corr_cholesky_mtx_1, corr_cholesky_mtx_2, &        ! intent(in)
               hydromet_pdf_params, silhs_config_flags, &         ! intent(in)
               l_uv_nudge, &                                      ! intent(in)
               l_tke_aniso, &                                     ! intent(in)
               l_standard_term_ta, &                              ! intent(in)
               l_single_C2_Skw, &                                 ! intent(in)
               X_nl_all_levs, X_mixt_comp_all_levs, &             ! intent(out)
               lh_sample_point_weights )                          ! intent(out)

! Description:
!   Generate sample points of moisture, temperature, et cetera for the purpose
!   of computing tendencies with a microphysics or radiation scheme.
!
! References:
! https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:overview_silhs
!-------------------------------------------------------------------------------

    use array_index, only: &
      iiPDF_chi    ! Variables

    use transform_to_pdf_module, only: &
      transform_uniform_samples_to_pdf      ! Procedure

    use output_2D_samples_module, only: &
      output_2D_lognormal_dist_file, & ! Procedure(s)
      output_2D_uniform_dist_file

    use pdf_parameter_module, only: &
      pdf_parameter  ! Type

    use hydromet_pdf_parameter_module, only: &
      hydromet_pdf_parameter ! Type

    use constants_clubb, only: &
      fstderr, & ! Constant(s)
      zero, &
      one, &
      rc_tol

    use clubb_precision, only: &
      core_rknd, &
      stat_rknd

    use parameters_silhs, only: &
      vert_decorr_coef, & ! Variable(s)
      silhs_config_flags_type ! Type

    use error_code, only: &
      clubb_at_least_debug_level  ! Procedure
      
    use fill_holes, only: &
      vertical_avg  ! Procedure
      
    use grid_class, only: &
      gr
      
    use stats_variables, only: &
      l_stats_samp      ! Variable(s)

    implicit none

    ! External
    intrinsic :: allocated, mod, maxloc, epsilon, transpose

    ! Parameter Constants

    integer, parameter :: &
      d_uniform_extra = 2   ! Number of variables that are included in the uniform sample but not in
                            ! the lognormal sample. Currently:
                            !
                            ! pdf_dim+1: Mixture component, for choosing PDF component
                            ! pdf_dim+2: Precipitation fraction, for determining precipitation

    ! ---------------- Input Variables ----------------
    integer, intent(in) :: &
      iter,            & ! Model iteration (time step) number
      pdf_dim,         & ! Number of variables to sample
      num_samples,     & ! Number of samples per variable
      sequence_length, & ! nt_repeat/num_samples; number of timesteps before sequence repeats
      nz                 ! Number of vertical model levels

    type(pdf_parameter), intent(in) :: &
      pdf_params ! PDF parameters       [units vary]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      delta_zm, &  ! Difference in momentum altitudes    [m]
      rcm          ! Liquid water mixing ratio          [kg/kg]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      Lscale       ! Turbulent mixing length            [m]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      rho_ds_zt    ! Dry, static density on thermo. levels    [kg/m^3]

    logical, intent(in) :: &
      l_calc_weights_all_levs_itime ! determines if vertically correlated sample points are needed
      
    ! More Input Variables!
    real( kind = core_rknd ), dimension(pdf_dim,pdf_dim,nz), intent(in) :: &
      corr_cholesky_mtx_1, & ! Correlations Cholesky matrix (1st comp.)  [-]
      corr_cholesky_mtx_2    ! Correlations Cholesky matrix (2nd comp.)  [-]

    real( kind = core_rknd ), dimension(pdf_dim,nz), intent(in) :: &
      mu1,    & ! Means of the hydrometeors, 1st comp. (chi, eta, w, <hydrometeors>)  [units vary]
      mu2,    & ! Means of the hydrometeors, 2nd comp. (chi, eta, w, <hydrometeors>)  [units vary]
      sigma1, & ! Stdevs of the hydrometeors, 1st comp. (chi, eta, w, <hydrometeors>) [units vary]
      sigma2    ! Stdevs of the hydrometeors, 2nd comp. (chi, eta, w, <hydrometeors>) [units vary]

    type(hydromet_pdf_parameter), dimension(nz), intent(in) :: &
      hydromet_pdf_params ! Hydrometeor PDF parameters  [units vary]

    type(silhs_config_flags_type), intent(in) :: &
      silhs_config_flags ! Flags for the SILHS sampling code [-]

    logical, intent(in) :: &
      l_uv_nudge,         & ! For wind speed nudging.
      l_tke_aniso,        & ! For anisotropic turbulent kinetic energy, i.e.
                            ! TKE = 1/2 (u'^2 + v'^2 + w'^2)
      l_standard_term_ta, & ! Use the standard discretization for the turbulent advection terms.
                            ! Setting to .false. means that a_1 and a_3 are pulled outside of the
                            ! derivative in advance_wp2_wp3_module.F90 and in
                            ! advance_xp2_xpyp_module.F90.
      l_single_C2_Skw       ! Use a single Skewness dependent C2 for rtp2, thlp2, and rtpthlp
    
    ! ---------------- Output Variables ----------------
    real( kind = core_rknd ), intent(out), dimension(nz,num_samples,pdf_dim) :: &
      X_nl_all_levs ! Sample that is transformed ultimately to normal-lognormal

    integer, intent(out), dimension(nz,num_samples) :: &
      X_mixt_comp_all_levs ! Which mixture component we're in

    real( kind = core_rknd ), intent(out), dimension(nz,num_samples) :: &
      lh_sample_point_weights ! Weight of each sample point

    ! ---------------- Local variables ----------------
    real( kind = core_rknd ), dimension(nz,num_samples,(pdf_dim+d_uniform_extra)) :: &
      X_u_all_levs ! Sample drawn from uniform distribution

    integer :: &
      k_lh_start, & ! Height for preferentially sampling within cloud
      k, sample, i, j, km1, kp1  ! Loop iterators

    logical, dimension(nz,num_samples) :: &
      l_in_precip   ! Whether sample is in precipitation

    logical :: l_error, l_error_in_sub
    
    real( kind = core_rknd ), dimension(pdf_dim,pdf_dim,nz) :: &
      Sigma_Cholesky1, &  ! Correlations Cholesky matrix 1 [-]
      Sigma_Cholesky2     ! Correlations Cholesky matrix 2 [-]
      
    real( kind = core_rknd ), dimension(nz,num_samples) :: &
      cloud_frac ! Cloud fraction for grid level and sample
      
    real(kind = core_rknd), dimension(nz) :: &
      precip_frac_1, &  ! Array used to store hydromet_pdf_params(:)%precip_frac_1
      precip_frac_2     ! Array used to store hydromet_pdf_params(:)%precip_frac_2
      
    real( kind = core_rknd ), dimension(nz) :: &
      Lscale_vert_avg, &  ! 3pt vertical average of Lscale                    [m]
      X_vert_corr         ! Vertical correlations between height levels       [-]
      
    real(kind = core_rknd), dimension(nz,num_samples,pdf_dim+d_uniform_extra) :: &
      rand_pool ! Array of randomly generated numbers

    ! ---------------- Begin Code ----------------

    l_error = .false.

    ! Sanity checks for l_lh_importance_sampling
    if ( silhs_config_flags%l_lh_importance_sampling .and. sequence_length /= 1 ) then
      write(fstderr,*) "Cloud weighted sampling requires sequence length be equal to 1."
      stop "Fatal error."
    end if

    !--------------------------------------------------------------
    ! Latin hypercube sampling
    !--------------------------------------------------------------
        
#ifdef _OPENACC
    if ( .not. silhs_config_flags%l_lh_straight_mc ) then
      stop "CLUBB ERROR: Running SILHS with OpenACC requires lh_straight_mc=true"
    end if
#endif

    ! Copy type arrays to contiguous arrays, so they can be copied to the GPU
    precip_frac_1 = hydromet_pdf_params(:)%precip_frac_1
    precip_frac_2 = hydromet_pdf_params(:)%precip_frac_2

    ! Compute k_lh_start, the starting vertical grid level 
    !   for SILHS sampling
    k_lh_start = compute_k_lh_start( nz, rcm, pdf_params, &
                                     silhs_config_flags%l_rcm_in_cloud_k_lh_start, &
                                     silhs_config_flags%l_random_k_lh_start )
                                     
    ! Calculate possible Sigma_Cholesky values
    ! Row-wise multiply of the elements of a lower triangular matrix.
    do k = 1, nz
      do i = 1, pdf_dim
        do j = 1, i
          ! Calculate possible Sigma_Cholesky values
          Sigma_Cholesky1(i,j,k) = corr_cholesky_mtx_1(i,j,k) * sigma1(i,k)
          Sigma_Cholesky2(i,j,k) = corr_cholesky_mtx_2(i,j,k) * sigma2(i,k)
        end do
      end do
    end do
    
    ! Determine 3pt vertically averaged Lscale
    if ( silhs_config_flags%l_Lscale_vert_avg ) then
      do k = 1, nz, 1
        kp1 = min( k+1, nz )
        km1 = max( k-1, 1 )
        Lscale_vert_avg(k) = vertical_avg &
                             ( (kp1-km1+1), rho_ds_zt(km1:kp1), &
                               Lscale(km1:kp1), gr%dzt(km1:kp1) )
      end do
    else
        Lscale_vert_avg = Lscale 
    end if
    
    ! Compute the vertical correlation for arbitrary overlap, using
    !   density weighted 3pt averaged Lscale and the difference in height levels (delta_zm)
    X_vert_corr(1:nz) = exp( -vert_decorr_coef * ( delta_zm(1:nz) / Lscale_vert_avg(1:nz) ) )

    if ( silhs_config_flags%l_max_overlap_in_cloud ) then
      where ( rcm > rc_tol )
       X_vert_corr = one
      end where
    end if
    
    ! Assertion check for the vertical correlation
    if ( clubb_at_least_debug_level( 1 ) ) then
      if ( any( X_vert_corr > one ) .or. any( X_vert_corr < zero ) ) then
        write(fstderr,*) "The vertical correlation in latin_hypercube_driver"// &
          "is not in the correct range"
        do k = 1, nz
          write(fstderr,*) "k = ", k,  "Vert. correlation = ", X_vert_corr(k)
        end do
        stop "Fatal error in vertical_overlap_driver"
      end if ! Some correlation isn't between [0,1]
    end if ! clubb_at_least_debug_level 1
    
    !$acc data create( rand_pool, X_u_all_levs ) &
    !$acc&     copyin( X_vert_corr ) &
    !$acc& async(1)

    ! Generate pool of random numbers
    call generate_random_pool( nz, pdf_dim, num_samples, d_uniform_extra, & ! Intent(in)
                               rand_pool )                                  ! Intent(out)
                               
    ! Generate all uniform samples, based on the rand pool
    call generate_all_uniform_samples( &
           iter, pdf_dim, d_uniform_extra, num_samples, sequence_length, & ! Intent(in)
           nz, k_lh_start, X_vert_corr(:), rand_pool(:,:,:),             & ! Intent(in)
           pdf_params%cloud_frac_1(:),                                   & ! Intent(in)
           pdf_params%cloud_frac_2(:),                                   & ! Intent(in)
           pdf_params%mixt_frac(:), hydromet_pdf_params(:),              & ! Intent(in)
           silhs_config_flags%cluster_allocation_strategy,               & ! Intent(in)
           silhs_config_flags%l_lh_importance_sampling,                  & ! Intent(in)
           silhs_config_flags%l_lh_straight_mc,                          & ! Intent(in)
           silhs_config_flags%l_lh_clustered_sampling,                   & ! Intent(in)
           silhs_config_flags%l_lh_limit_weights,                        & ! Intent(in)
           silhs_config_flags%l_lh_var_frac,                             & ! Intent(in)
           silhs_config_flags%l_lh_normalize_weights,                    & ! Intent(in)
           l_calc_weights_all_levs_itime,                                & ! Intent(in)
           X_u_all_levs(:,:,:), lh_sample_point_weights(:,:) )             ! Intent(out)
    
    !$acc data create( cloud_frac, l_in_precip ) &
    !$acc&     copyin( pdf_params,pdf_params%mixt_frac,pdf_params%cloud_frac_1, &
    !$acc&             pdf_params%cloud_frac_2, precip_frac_1,precip_frac_2 ) &
    !$acc& async(3)
    
    !$acc parallel loop collapse(2) default(present) async(1) wait(3)
    do sample = 1, num_samples 
      do k = 1, nz
            
        ! Determine mixture component for all levels
        if ( X_u_all_levs(k,sample,pdf_dim+1) < pdf_params%mixt_frac(k) ) then
          
          ! Set pdf component indicator to 1 for this sample and vertical level
          X_mixt_comp_all_levs(k,sample) = 1
          
          ! Copy 1st component values
          cloud_frac(k,sample) = pdf_params%cloud_frac_1(k)
          
          ! Determine precipitation
          if ( X_u_all_levs(k,sample,pdf_dim+2) < precip_frac_1(k) ) then
            l_in_precip(k,sample) = .true.
          else
            l_in_precip(k,sample) = .false.
          end if
          
        else
          
          ! Set pdf component indicator to 2 for this sample and vertical level
          X_mixt_comp_all_levs(k,sample) = 2
          
          ! Copy 2nd component values
          cloud_frac(k,sample) = pdf_params%cloud_frac_2(k)
          
          ! Determine precipitation
          if ( X_u_all_levs(k,sample,pdf_dim+2) < precip_frac_2(k) ) then
            l_in_precip(k,sample) = .true.
          else
            l_in_precip(k,sample) = .false.
          end if
          
        end if

      end do 
    end do

    ! Generate LH sample, represented by X_u and X_nl, for level k
    ! Transform the uniformly distributed samples to
    !   ones distributed according to CLUBB's PDF.
    call transform_uniform_samples_to_pdf &
         ( nz, num_samples, pdf_dim, d_uniform_extra, & ! In
           Sigma_Cholesky1(:,:,:), Sigma_Cholesky2(:,:,:), &
           mu1(:,:), mu2(:,:), X_mixt_comp_all_levs(:,:), &
           X_u_all_levs(:,:,:), cloud_frac(:,:), & ! In
           l_in_precip(:,:), & ! In
           X_nl_all_levs(:,:,:) ) ! Out
           
    
    if ( l_stats_samp ) then
      !$acc update host(X_u_all_levs,l_in_precip,lh_sample_point_weights) wait
      call stats_accumulate_uniform_lh( nz, num_samples, l_in_precip, X_mixt_comp_all_levs, &
                                        X_u_all_levs(:,:,iiPDF_chi), pdf_params, &
                                        lh_sample_point_weights, k_lh_start )
    end if

    if ( l_output_2D_lognormal_dist ) then
      !$acc update host(X_nl_all_levs) wait
      
      ! Eric Raut removed lh_rt and lh_thl from call to output_2D_lognormal_dist_file
      ! because they are no longer generated in generate_silhs_sample.
      call output_2D_lognormal_dist_file( nz, num_samples, pdf_dim, &
                                          real(X_nl_all_levs, kind = stat_rknd), &
                                          l_uv_nudge, &
                                          l_tke_aniso, &
                                          l_standard_term_ta, &
                                          l_single_C2_Skw )
    end if
    
    if ( l_output_2D_uniform_dist ) then
      !$acc update host(X_u_all_levs,X_mixt_comp_all_levs,lh_sample_point_weights) wait
      call output_2D_uniform_dist_file( nz, num_samples, pdf_dim+2, &
                                        X_u_all_levs, &
                                        X_mixt_comp_all_levs, &
                                        lh_sample_point_weights, &
                                        l_uv_nudge, &
                                        l_tke_aniso, &
                                        l_standard_term_ta, &
                                        l_single_C2_Skw )
    end if

    ! Various nefarious assertion checks
    if ( clubb_at_least_debug_level( 2 ) ) then
      
      !$acc update host(X_u_all_levs,X_mixt_comp_all_levs,X_nl_all_levs) wait

      ! Simple assertion check to ensure uniform variates are in the appropriate
      ! range
      if ( any( X_u_all_levs <= zero .or. X_u_all_levs >= one ) ) then
        write(fstderr,*) "A uniform variate was not in the correct range."
        l_error = .true.
      end if

      do k=2, nz

        call assert_consistent_cloud_frac( pdf_params%chi_1(k), pdf_params%chi_2(k), & 
                                   pdf_params%cloud_frac_1(k), pdf_params%cloud_frac_2(k), &
                                   pdf_params%stdev_chi_1(k), pdf_params%stdev_chi_2(k), & 
                                   l_error_in_sub )
        l_error = l_error .or. l_error_in_sub

        ! Check for correct transformation in normal space
        call assert_correct_cloud_normal( num_samples, X_u_all_levs(k,:,iiPDF_chi), & ! In
                                          X_nl_all_levs(k,:,iiPDF_chi), & ! In
                                          X_mixt_comp_all_levs(k,:), & ! In
                                          pdf_params%cloud_frac_1(k), & ! In
                                          pdf_params%cloud_frac_2(k), & ! In
                                          l_error_in_sub ) ! Out
        l_error = l_error .or. l_error_in_sub

      end do ! k=2, nz

    end if ! clubb_at_least_debug_level( 2 )
    
    !$acc end data
    !$acc end data

    ! Stop the run if an error occurred
    if ( l_error ) then
      stop "Fatal error in generate_silhs_sample"
    end if

    return
  end subroutine generate_silhs_sample
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
  subroutine generate_random_pool( nz, pdf_dim, num_samples, d_uniform_extra, &
                                   rand_pool )
  ! Description:
  !     This subroutine populates rand_pool with random numbers. There are
  !     different requirements for generating random numbers on a CPU vs a
  !     GPU, so the operation of this procedure depends on where or not 
  !     openacc has been specified at compile time.
  !
  ! References:
  !   clubb ticket 869
  !
  ! Author: Gunther Huebler
  !----------------------------------------------------------------------------
    
    use clubb_precision, only: &
      core_rknd      ! Constant     
    
    use generate_uniform_sample_module, only: &
      rand_uniform_real ! Procedure
      
#ifdef _OPENACC
    use curand, only: &
      curandSetPseudoRandomGeneratorSeed, & ! Procedures
      curandCreateGenerator,              &
      curandGenerate,                     &
      curandGenerator,                    & ! Type
      CURAND_RNG_PSEUDO_DEFAULT             ! Parameter
#endif
    
    implicit none 
    
    ! ------------------- Input Variables -------------------
    
    integer, intent(in) :: &
      nz,               & ! Number of vertical levels
      pdf_dim,          & ! Variates
      num_samples,      & ! Number of samples
      d_uniform_extra     ! Uniform variates included in uniform sample but not
                          !  in normal/lognormal sample
      
    ! ------------------- Output Variables -------------------

    real(kind=core_rknd), dimension(nz,num_samples,pdf_dim+d_uniform_extra), intent(out) :: &
      rand_pool   ! Pool of random reals, generated from integers
      
    ! ------------------- Local Variables -------------------
    
#ifdef _OPENACC
    real(kind = core_rknd), dimension(nz,num_samples,pdf_dim+d_uniform_extra) :: &
      rand_pool_int ! Random intgers, stored as reals to avoid data conversion
      
    integer :: &
      r_status  ! Integer use to call curand functions
    
    logical, save :: &
      l_first_iter = .true.  ! First iteration indicator
                           ! The curand generator needs to be inialized, but only once

    type(curandGenerator) :: &
      cu_gen  ! curand generator variable
    
    ! Parameters stole from mt90 module, used to convert random integers
    ! to random reals between [0,1]
    real(kind=core_rknd), parameter :: p231       = 2147483648.0_core_rknd
    real(kind=core_rknd), parameter :: p232       = 4294967296.0_core_rknd
    real(kind=core_rknd), parameter :: pi232      = 1.0_core_rknd / p232
    real(kind=core_rknd), parameter :: p231_5d232 = ( p231 + 0.5_core_rknd ) / p232
#endif
    
    integer :: k, i, sample ! Loop variables
      
    ! ---------------- Begin Code ----------------
    
#ifdef _OPENACC

    ! Generate randoms on GPU
    
    ! If first iteration, intialize generator
    if ( l_first_iter ) then
      l_first_iter = .false.
      r_status = curandCreateGenerator( cu_gen, CURAND_RNG_PSEUDO_DEFAULT )
      r_status = curandSetPseudoRandomGeneratorSeed( cu_gen, 252435 )
    end if
    
    !$acc data create( rand_pool_int ) async(1)
    
    !$acc host_data use_device(rand_pool_int) 
    r_status = curandGenerate( cu_gen, rand_pool_int, nz*num_samples*(pdf_dim+d_uniform_extra) )
    !$acc end host_data

    ! Populate rand_pool with random reals, using random integers
    !$acc parallel loop collapse(3) default(present) async(1)
    do i=1, pdf_dim+d_uniform_extra
      do sample=1, num_samples
        do k = 1, nz
          rand_pool(k,sample,i) = rand_pool_int(k,sample,i) * pi232 + p231_5d232
        end do
      end do
    end do
    
    !$acc end data
    
#else

    ! Generate randoms on CPU

    ! Populate rand_pool with a generator designed for a CPU
    do i=1, pdf_dim+d_uniform_extra
      do sample=1, num_samples
        do k = 1, nz
          rand_pool(k,sample,i) = rand_uniform_real()
        end do
      end do
    end do

#endif
    
  end subroutine generate_random_pool

!-------------------------------------------------------------------------------
  subroutine generate_all_uniform_samples( &
               iter, pdf_dim, d_uniform_extra, num_samples, sequence_length, & ! Intent(in)
               nz, k_lh_start, X_vert_corr, rand_pool,                       & ! Intent(in)
               cloud_frac_1,                                                 & ! Intent(in)
               cloud_frac_2,                                                 & ! Intent(in)
               mixt_frac, hydromet_pdf_params,                               & ! Intent(in)
               cluster_allocation_strategy,                                  & ! Intent(in)
               l_lh_importance_sampling,                                     & ! Intent(in)
               l_lh_straight_mc,                                             & ! Intent(in)
               l_lh_clustered_sampling,                                      & ! Intent(in)
               l_lh_limit_weights,                                           & ! Intent(in)
               l_lh_var_frac,                                                & ! Intent(in)
               l_lh_normalize_weights,                                       & ! Intent(in)
               l_calc_weights_all_levs_itime,                                & ! Intent(in)
               X_u_all_levs, lh_sample_point_weights )                         ! Intent(out)
  ! Description:
  !   Generates uniform samples for all vertical levels, samples, and variates.
  !   Applys Latin Hypercude and importance sampling where conigured
  !
  ! References:
  !   V. E. Larson and D. P. Schanen, 2013. The Subgrid Importance Latin
  !   Hypercube Sampler (SILHS): a multivariate subcolumn generator
  !----------------------------------------------------------------------

    ! Included Modules
    use clubb_precision, only: &
      core_rknd                   ! Precision
      
    use parameters_silhs, only: &
      single_prec_thresh       ! Constant

    use constants_clubb, only: &
      one, fstderr                ! Constant(s)

    use hydromet_pdf_parameter_module, only: &
      hydromet_pdf_parameter      ! Type

    use generate_uniform_sample_module, only: &
      rand_uniform_real, &        ! Procedure(s)
      generate_uniform_lh_sample

    use silhs_importance_sample_module, only: &
      importance_sampling_driver, & ! Procedure(s)
      cloud_weighted_sampling_driver

    use latin_hypercube_arrays, only: &
      one_height_time_matrix      ! Variable

    use array_index, only: &
      iiPDF_chi                   ! Variable

    implicit none

    ! Local Constants
    logical, parameter :: &
      l_lh_old_cloud_weighted  = .false. ! Use the old method of importance sampling that
                                         ! places one point in cloud and one point out of
                                         ! cloud

    ! ------------------ Input Variables ------------------
    integer, intent(in) :: &
      iter,              &        ! Model iteration number
      pdf_dim,           &        ! Number of variates in CLUBB's PDF
      d_uniform_extra,   &        ! Uniform variates included in uniform sample but not
                                  !  in normal/lognormal sample
      num_samples,       &        ! Number of SILHS sample points
      sequence_length,   &        ! Number of timesteps before new sample points are picked
      k_lh_start,        &
      nz

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      cloud_frac_1, cloud_frac_2, &     ! The PDF parameters at k_lh_start
      mixt_frac, &
      X_vert_corr

    type(hydromet_pdf_parameter), dimension(nz), intent(in) :: &
      hydromet_pdf_params

    real(kind = core_rknd), dimension(nz,num_samples,pdf_dim+d_uniform_extra), intent(in) :: &
      rand_pool ! Array of randomly generated numbers

    integer, intent(in) :: &
      cluster_allocation_strategy ! Strategy for distributing sample points

    logical, intent(in) :: &
      l_lh_importance_sampling, & ! Do importance sampling (SILHS)
      l_lh_straight_mc, &         ! Do not apply LH or importance sampling at all (SILHS)
      l_lh_clustered_sampling, &  ! Use prescribed probability sampling with clusters (SILHS)
      l_lh_limit_weights, &       ! Ensure weights stay under a given value
      l_lh_var_frac, &            ! Prescribe variance fractions
      l_lh_normalize_weights, &   ! Normalize weights to sum to num_samples
      l_calc_weights_all_levs_itime

    ! ------------------ Output Variables ------------------
    real( kind = core_rknd ), dimension(nz,num_samples,pdf_dim+d_uniform_extra), intent(out) :: &
      X_u_all_levs              ! Uniform sample at k_lh_start

    real( kind = core_rknd ), dimension(nz,num_samples), intent(out) :: &
      lh_sample_point_weights     ! Weight of each sample point (all equal to one if importance
                                  ! sampling is not used)

    ! ------------------ Local Variables ------------------
    integer :: &
      k, i, sample

    !------------------ Begin Code ------------------

    ! Sanity check
    if ( l_lh_old_cloud_weighted .and. mod( num_samples, 2 ) /= 0 ) then
      write(fstderr,*) "Old cloud weighted sampling requires num_samples to be divisible by 2."
      stop "Fatal error."
    end if
    
    if ( .not. l_calc_weights_all_levs_itime ) then

      ! Generate random samples at k_lh_start, then use compute_arb_overlap to 
      ! populate the rest of the uniform samples

      if ( l_lh_straight_mc ) then

        ! Do a straight Monte Carlo sample without LH or importance sampling.
        
        !$acc parallel loop collapse(2) default(present) async(1)
        do i=1, pdf_dim+d_uniform_extra
          do sample=1, num_samples
            X_u_all_levs(k_lh_start,sample,i) = max( single_prec_thresh, &
                                                     min( one - single_prec_thresh, &
                                                          rand_pool(k_lh_start,sample,i) ) )
          end do
        end do

        !$acc parallel loop collapse(2) default(present) async(1)
        do sample = 1, num_samples
          do k = 1, nz
            ! Importance sampling is not performed, so all sample points have the same weight!!
            lh_sample_point_weights(k,sample) = one
          end do
        end do

      else ! .not. l_lh_straight_mc
      
        ! Generate a uniformly distributed Latin hypercube sample
        call generate_uniform_lh_sample( iter, num_samples, sequence_length, & ! Intent(in)
                                         pdf_dim+d_uniform_extra,            & ! Intent(in)
                                         X_u_all_levs(k_lh_start,:,:) )        ! Intent(out)
                              
        if ( l_lh_importance_sampling ) then

          if ( l_lh_old_cloud_weighted ) then

            call cloud_weighted_sampling_driver &
                 ( num_samples, one_height_time_matrix(:,iiPDF_chi), & ! In
                   one_height_time_matrix(:,pdf_dim+1), & ! In
                   cloud_frac_1(k_lh_start), cloud_frac_2(k_lh_start), & ! In
                   mixt_frac(k_lh_start), & ! In
                   X_u_all_levs(k_lh_start,:,iiPDF_chi), & ! In/Out
                   X_u_all_levs(k_lh_start,:,pdf_dim+1), & ! In/Out
                   lh_sample_point_weights(k_lh_start,:) ) ! Out

          else ! .not. l_lh_old_cloud_weighted

            call importance_sampling_driver &
                 ( num_samples,                                               & ! In
                   cloud_frac_1(k_lh_start), cloud_frac_2(k_lh_start),        & ! In
                   mixt_frac(k_lh_start), hydromet_pdf_params(k_lh_start),    & ! In
                   cluster_allocation_strategy, l_lh_clustered_sampling,      & ! In
                   l_lh_limit_weights, l_lh_var_frac, l_lh_normalize_weights, & ! In
                   X_u_all_levs(k_lh_start,:,iiPDF_chi),                      & ! In/Out
                   X_u_all_levs(k_lh_start,:,pdf_dim+1),                      & ! In/Out
                   X_u_all_levs(k_lh_start,:,pdf_dim+2),                      & ! In/Out
                   lh_sample_point_weights(k_lh_start,:) )                      ! Out

          end if ! l_lh_old_cloud_weighted
          
          ! Clip uniform sample points to expected range                                 
          X_u_all_levs(k_lh_start,:,:) = max( single_prec_thresh, &
                                min( one - single_prec_thresh, X_u_all_levs(k_lh_start,:,:) ) )

          
          do k = 1, nz
            lh_sample_point_weights(k,:) = lh_sample_point_weights(k_lh_start,:)
          end do

        else

          ! No importance sampling is performed, so all sample points have the same weight.
          lh_sample_point_weights(:,:) = one

        end if ! l_lh_importance_sampling


      end if ! l_lh_straight_mc
      
      ! Generate uniform sample at other grid levels by vertically correlating them
      ! https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:vert_corr
      call compute_arb_overlap( nz, num_samples, pdf_dim, d_uniform_extra, & ! In
                                k_lh_start, X_vert_corr, rand_pool,        & ! In 
                                X_u_all_levs )                               ! Out
                                
    else 
      
      ! Generate random samples for all vertical levels, samples, and variates
      
      if ( l_lh_straight_mc ) then

        ! Do a straight Monte Carlo sample without LH or importance sampling.
        
        !$acc parallel loop collapse(3) default(present) async(1)
        do i = 1, pdf_dim+d_uniform_extra
          do sample = 1, num_samples
            do k = 1, nz
              X_u_all_levs(k,sample,i) = max( single_prec_thresh, &
                                              min( one - single_prec_thresh, &
                                                   rand_pool(k,sample,i) ) )
            end do
          end do
        end do

        ! Importance sampling is not performed, so all sample points have the same weight!!
        lh_sample_point_weights(:,:)  =  one

      else ! .not. l_lh_straight_mc
        
        do k = 1, nz

          ! Generate a uniformly distributed Latin hypercube sample
          call generate_uniform_lh_sample( iter, num_samples, sequence_length, & ! Intent(in)
                                           pdf_dim+d_uniform_extra,            & ! Intent(in)
                                           X_u_all_levs(k,:,:) )                 ! Intent(out)

          if ( l_lh_importance_sampling ) then

            if ( l_lh_old_cloud_weighted ) then

              call cloud_weighted_sampling_driver &
                   ( num_samples, one_height_time_matrix(:,iiPDF_chi), & ! In
                     one_height_time_matrix(:,pdf_dim+1),              & ! In
                     cloud_frac_1(k), cloud_frac_2(k),                 & ! In
                     mixt_frac(k),                                     & ! In
                     X_u_all_levs(k,:,iiPDF_chi),                      & ! In/Out
                     X_u_all_levs(k,:,pdf_dim+1),                      & ! In/Out
                     lh_sample_point_weights(k,:) )                      ! Out

            else ! .not. l_lh_old_cloud_weighted

              call importance_sampling_driver &
                   ( num_samples,                                               & ! In
                     cloud_frac_1(k), cloud_frac_2(k),                          & ! In
                     mixt_frac(k), hydromet_pdf_params(k),                      & ! In
                     cluster_allocation_strategy, l_lh_clustered_sampling,      & ! In
                     l_lh_limit_weights, l_lh_var_frac, l_lh_normalize_weights, & ! In
                     X_u_all_levs(k,:,iiPDF_chi),                               & ! In/Out
                     X_u_all_levs(k,:,pdf_dim+1),                               & ! In/Out
                     X_u_all_levs(k,:,pdf_dim+2),                               & ! In/Out
                     lh_sample_point_weights(k,:) )                               ! Out

            end if ! l_lh_old_cloud_weighted

          else

            ! No importance sampling is performed, so all sample points have the same weight.
            lh_sample_point_weights(k,:) = one

          end if ! l_lh_importance_sampling
          
          ! Clip uniform sample points to expected range                                 
          X_u_all_levs(k,:,:) = max( single_prec_thresh, &
                                min( one - single_prec_thresh, X_u_all_levs(k,:,:) ) )
        
        end do

      end if ! l_lh_straight_mc
      
    end if ! .not. l_calc_weights_all_levs_itime

    return
  end subroutine generate_all_uniform_samples
!-------------------------------------------------------------------------------

!-----------------------------------------------------------------------
  function compute_k_lh_start( nz, rcm, pdf_params, &
                               l_rcm_in_cloud_k_lh_start, &
                               l_random_k_lh_start ) result( k_lh_start )

  ! Description:
  !   Determines the starting SILHS sample level

  ! References:
  !   None
  !-----------------------------------------------------------------------

    ! Included Modules
    use clubb_precision, only: &
      core_rknd      ! Constant

    use constants_clubb, only: &
      cloud_frac_min ! Constant

    use pdf_parameter_module, only: &
      pdf_parameter  ! Type

    use pdf_utilities, only: &
      compute_mean_binormal  ! Procedure

    use math_utilities, only: &
      rand_integer_in_range  ! Procedure

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz          ! Number of vertical levels

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      rcm         ! Liquid water mixing ratio               [kg/kg]

    type(pdf_parameter), intent(in) :: &
      pdf_params  ! PDF parameters       [units vary]

    logical, intent(in) :: &
      l_rcm_in_cloud_k_lh_start, & ! Determine k_lh_start based on maximum within-cloud rcm
      l_random_k_lh_start          ! k_lh_start found randomly between max rcm and rcm_in_cloud

    ! Output Variable
    integer :: &
      k_lh_start  ! Starting SILHS sample level

    ! Local Variables
    integer :: &
      k_lh_start_rcm_in_cloud, &
      k_lh_start_rcm

    real( kind = core_rknd ), dimension(nz) :: &
      rcm_pdf, cloud_frac_pdf

  !-----------------------------------------------------------------------

    !----- Begin Code -----
    if ( l_rcm_in_cloud_k_lh_start .or. l_random_k_lh_start ) then
      rcm_pdf = compute_mean_binormal( pdf_params%rc_1, pdf_params%rc_2, pdf_params%mixt_frac )
      cloud_frac_pdf = compute_mean_binormal( pdf_params%cloud_frac_1, pdf_params%cloud_frac_2, &
                                              pdf_params%mixt_frac )
      k_lh_start_rcm_in_cloud = maxloc( rcm_pdf / max( cloud_frac_pdf, cloud_frac_min ), 1 )
    end if

    if ( .not. l_rcm_in_cloud_k_lh_start .or. l_random_k_lh_start ) then
      k_lh_start_rcm    = maxloc( rcm, 1 )
    end if

    if ( l_random_k_lh_start ) then
      if ( k_lh_start_rcm_in_cloud == k_lh_start_rcm ) then
        k_lh_start = k_lh_start_rcm
      else
        ! Pick a random height level between k_lh_start_rcm and
        ! k_lh_start_rcm_in_cloud
        if ( k_lh_start_rcm_in_cloud > k_lh_start_rcm ) then
          k_lh_start = rand_integer_in_range( k_lh_start_rcm, k_lh_start_rcm_in_cloud )
        else if ( k_lh_start_rcm > k_lh_start_rcm_in_cloud ) then
          k_lh_start = rand_integer_in_range( k_lh_start_rcm_in_cloud, k_lh_start_rcm )
        end if
      end if
    else if ( l_rcm_in_cloud_k_lh_start ) then
      k_lh_start = k_lh_start_rcm_in_cloud
    else ! .not. l_random_k_lh_start .and. .not. l_rcm_in_cloud_k_lh_start
      k_lh_start = k_lh_start_rcm
    end if

    ! If there's no cloud k_lh_start appears to end up being 1. Check if
    ! k_lh_start is 1 or nz and set it to the middle of the domain in that
    ! case.
    if ( k_lh_start == nz .or. k_lh_start == 1 ) then
      k_lh_start = nz / 2
    end if

    return
  end function compute_k_lh_start
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  subroutine clip_transform_silhs_output( nz, num_samples,                & ! In
                                          pdf_dim, hydromet_dim,          & ! In
                                          X_mixt_comp_all_levs,           & ! In
                                          X_nl_all_levs,                  & ! Inout
                                          pdf_params, l_use_Ncn_to_Nc,    & ! In
                                          lh_rt_clipped, lh_thl_clipped,  & ! Out
                                          lh_rc_clipped, lh_rv_clipped,   & ! Out
                                          lh_Nc_clipped                   ) ! Out

  ! Description:
  !   Derives from the SILHS sampling structure X_nl_all_levs the variables
  !   rt, thl, rc, rv, and Nc, for all sample points and height levels.

  ! References:
  !   ticket:751
  !-----------------------------------------------------------------------

    ! Included Modules
    use clubb_precision, only: &
        core_rknd

    use constants_clubb, only: &
        zero, &   ! Constant(s)
        rt_tol

    use pdf_parameter_module, only: &
        pdf_parameter ! Type

    use index_mapping, only: &
        hydromet2pdf_idx    ! Procedure(s)

    use fill_holes, only: &
        clip_hydromet_conc_mvr    ! Procedure(s)

    use array_index, only: &
        iiPDF_chi, &  ! Variable(s)
        iiPDF_eta, &
        iiPDF_Ncn

    use transform_to_pdf_module, only: &
        chi_eta_2_rtthl ! Awesome procedure

    implicit none

    ! ------------------- Input Variables -------------------
    logical, intent(in) :: &
      l_use_Ncn_to_Nc  ! Whether to call Ncn_to_Nc (.true.) or not (.false.);
                       ! Ncn_to_Nc might cause problems with the MG microphysics
                       ! since the changes made here (Nc-tendency) are not fed
                       ! into the microphysics

    integer, intent(in) :: &
      nz,           & ! Number of vertical levels
      num_samples,  & ! Number of SILHS sample points
      pdf_dim,      & ! Number of variates in X_nl
      hydromet_dim    ! Number of hydrometeor species

    integer, dimension(nz,num_samples), intent(in) :: &
      X_mixt_comp_all_levs   ! Which component this sample is in (1 or 2)

    type(pdf_parameter), intent(in) :: &
      pdf_params             ! **The** PDF parameters!

    ! ------------------- Input/Output Variable -------------------
    
    real( kind = core_rknd ), dimension(nz,num_samples,pdf_dim), intent(inout) :: &
      X_nl_all_levs    ! SILHS sample points    [units vary]

    ! ------------------- Output Variables -------------------
    real( kind = core_rknd ), dimension(nz,num_samples), intent(out) :: &
      lh_rt_clipped,  & ! rt generated from silhs sample points
      lh_thl_clipped, & ! thl generated from silhs sample points
      lh_rc_clipped,  & ! rc generated from silhs sample points
      lh_rv_clipped,  & ! rv generated from silhs sample points
      lh_Nc_clipped     ! Nc generated from silhs sample points

    ! ------------------- Local Variables -------------------

    real( kind = core_rknd ), dimension(nz,hydromet_dim) :: &
      hydromet_pts,         & ! Sample point column of hydrometeors    [un vary]
      hydromet_pts_clipped    ! Clipped sample point column of hydromet   [un v]

    integer :: &
      sample, k, hm_idx

    ! Flag to clip sample points of hydrometeor concentrations.
    logical, parameter :: &
      l_clip_hydromet_samples = .false.

  !-----------------------------------------------------------------------

    ! Calculate (and clip) the SILHS sample point values of rt, thl, rc, rv,
    ! and Nc.
        
    ! Compute lh_rt and lh_thl
    call chi_eta_2_rtthl( nz, num_samples, &
                          pdf_params%rt_1(:), pdf_params%thl_1(:),    & ! Intent(in)
                          pdf_params%rt_2(:), pdf_params%thl_2(:),    & ! Intent(in)
                          pdf_params%crt_1(:), pdf_params%cthl_1(:),  & ! Intent(in)
                          pdf_params%crt_2(:), pdf_params%cthl_2(:),  & ! Intent(in)
                          pdf_params%chi_1(:), pdf_params%chi_2(:),   & ! Intent(in)
                          X_nl_all_levs(:,:,iiPDF_chi),               & ! Intent(in) 
                          X_nl_all_levs(:,:,iiPDF_eta),               & ! Intent(in)
                          X_mixt_comp_all_levs(:,:),                  & ! Intent(in)
                          lh_rt_clipped, lh_thl_clipped               ) ! Intent(out)
    
    !$acc parallel loop default(present) async(1)
    do sample = 1, num_samples
      ! These parameters are not computed at the model lower level.
      lh_rt_clipped(1,sample)  = zero
      lh_thl_clipped(1,sample) = zero
      lh_rc_clipped(1,sample)  = zero
      lh_rv_clipped(1,sample)  = zero
      lh_Nc_clipped(1,sample)  = zero
    end do
    
    !$acc parallel loop collapse(2) default(present) async(1)
    do sample = 1, num_samples
      do k = 2, nz
    
        ! If necessary, clip rt      
        lh_rt_clipped(k,sample) = max( lh_rt_clipped(k,sample), rt_tol )
        
        ! Compute lh_rc, rc = chi * H(chi), where H(x) is the Heaviside step function
        lh_rc_clipped(k,sample) = max( X_nl_all_levs(k,sample,iiPDF_chi), zero )
        
        ! Clip lh_rc.
        lh_rc_clipped(k,sample) = min( lh_rc_clipped(k,sample), &
                                            lh_rt_clipped(k,sample) - rt_tol )
        
        ! Compute lh_rv
        lh_rv_clipped(k,sample) = lh_rt_clipped(k,sample) - lh_rc_clipped(k,sample)
        
        if ( l_use_Ncn_to_Nc ) then
           ! Compute lh_Nc, Nc = Ncn * H(chi), where H(x) is the Heaviside step function
           if ( X_nl_all_levs(k,sample,iiPDF_chi) > zero ) then
             lh_Nc_clipped(k,sample) = X_nl_all_levs(k,sample,iiPDF_Ncn)
           else
             lh_Nc_clipped(k,sample) = zero
           end if
        else
           lh_Nc_clipped(k,sample) = X_nl_all_levs(k,sample,iiPDF_Ncn)
        endif ! l_use_Ncn_to_Nc
        
      end do ! sample = 1, num_samples
    end do ! k = 2, nz


    ! Clip the SILHS sample point values of hydrometeor concentrations.
    if ( l_clip_hydromet_samples ) then
      
#ifdef _OPENACC
       stop "CLUBB ERROR: Running SILHS with OpenACC requires l_clip_hydromet_samples=false"
#endif
       ! Loop over all sample columns.
       do sample = 1, num_samples, 1

          ! Pack the SILHS hydrometeor sample points, which are stored in arrays
          ! with the size pdf_dim, into arrays with the size hydromet_dim.
          do hm_idx = 1, hydromet_dim, 1
             hydromet_pts(:,hm_idx) &
             = X_nl_all_levs(:,sample,hydromet2pdf_idx(hm_idx))
          enddo ! hm_idx = 1, hydromet_dim, 1

          ! Clip the hydrometeor concentration sample points based on
          ! maintaining the value of the hydrometeor mixing ratio sample points
          ! and satisfying the maximum allowable mean volume radius for that
          ! hydrometeor species.
          call clip_hydromet_conc_mvr( hydromet_dim, hydromet_pts, & ! In
                                       hydromet_pts_clipped )        ! Out

          ! Unpack the clipped SILHS hydrometeor sample points, which are stored
          ! in arrays with the size hydromet_dim, back into arrays with the size
          ! pdf_dim.
          do hm_idx = 1, hydromet_dim, 1
             X_nl_all_levs(:,sample,hydromet2pdf_idx(hm_idx)) &
             = hydromet_pts_clipped(:,hm_idx)
          enddo ! hm_idx = 1, hydromet_dim, 1

       enddo ! isample = 1, num_samples, 1

    endif ! l_clip_hydromet_samples


    return

  end subroutine clip_transform_silhs_output

!-----------------------------------------------------------------------

!-------------------------------------------------------------------------------
  subroutine assert_consistent_cloud_frac( chi_1, chi_2, & 
                                           cloud_frac_1, cloud_frac_2, &
                                           stdev_chi_1, stdev_chi_2, &
                                           l_error )

  ! Description:
  !   Performs an assertion check that cloud_frac_i is consistent with chi_i and
  !   stdev_chi_i in pdf_params for each PDF component.
  ! 
  ! References:
  !   Eric Raut
  !-----------------------------------------------------------------------

    use constants_clubb, only: &
      fstderr          ! Constant

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      chi_1, chi_2, &                   ! Discription of the distribution of chi 
      cloud_frac_1, cloud_frac_2, &     ! and cloud fraction [units vary]
      stdev_chi_1, stdev_chi_2
                       

    ! Output Variables
    logical, intent(out) :: &
      l_error          ! True if the assertion check fails

    ! Local Variables
    logical :: &
      l_error_in_sub

  !-----------------------------------------------------------------------

    !----- Begin Code -----
    l_error = .false.

    ! Perform assertion check for PDF component 1
    call assert_consistent_cf_component &
         ( chi_1, stdev_chi_1, cloud_frac_1, & ! Intent(in)
           l_error_in_sub )                                                    ! Intent(out)

    l_error = l_error .or. l_error_in_sub
    if ( l_error_in_sub ) then
      write(fstderr,*) "Cloud fraction is inconsistent in PDF component 1"
    end if

    ! Perform assertion check for PDF component 2
    call assert_consistent_cf_component &
         ( chi_2, stdev_chi_2, cloud_frac_2, & ! Intent(in)
           l_error_in_sub )                                                    ! Intent(out)

    l_error = l_error .or. l_error_in_sub
    if ( l_error_in_sub ) then
      write(fstderr,*) "Cloud fraction is inconsistent in PDF component 2"
    end if

    return
  end subroutine assert_consistent_cloud_frac
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
  subroutine assert_consistent_cf_component( mu_chi_i, sigma_chi_i, cloud_frac_i, &
                                             l_error )

  ! Description:
  !   Performs an assertion check that cloud_frac_i is consistent with chi_i and
  !   stdev_chi_i for a PDF component.
  !
  !   The SILHS sample generation process relies on precisely a cloud_frac
  !   amount of mass in the cloudy portion of the PDF of chi, that is, where
  !   chi > 0. In other words, the probability that chi > 0 should be exactly
  !   cloud_frac.
  !
  !   Stated even more mathematically, CDF_chi(0) = 1 - cloud_frac, where
  !   CDF_chi is the cumulative distribution function of chi. This can be
  !   expressed as invCDF_chi(1 - cloud_frac) = zero.
  !
  !   This subroutine uses ltqnorm, which is apparently a fancy name for the
  !   inverse cumulative distribution function of the standard normal
  !   distribution.

  ! References:
  !   Eric Raut
  !-----------------------------------------------------------------------

    ! Included Modules
    use clubb_precision, only: &
        core_rknd        ! Constant(s)

    use constants_clubb, only: &
        fstderr, & ! Constant(s)
        zero,    &
        one,     &
        chi_tol, &
        eps

    use transform_to_pdf_module, only: &
        ltqnorm          ! Procedure

    implicit none

    ! Local Constants
    real( kind = core_rknd ), parameter :: &
      ! Values below ltqnorm_min_arg (or above 1-ltqnorm_min_arg) will not be
      ! supplied as arguments to the ltqnorm function.
      ltqnorm_min_arg = 1.0e-5_core_rknd, &
      ! It will be verified that all values of the chi uniform variate will not
      ! trigger the in/out of cloud assertion error, except for those values
      ! within a box with the following half-width, centered around the cloud
      ! fraction in this component.
      box_half_width = 5.0e-6_core_rknd

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mu_chi_i,      & ! Mean of chi in a PDF component
      sigma_chi_i,   & ! Standard deviation of chi in a PDF component
      cloud_frac_i     ! Cloud fraction in a PDF component

    ! Output Variables
    logical, intent(out) :: &
      l_error          ! True if the assertion check fails

    ! Local Variables
    real( kind = core_rknd ) :: chi, chi_std_normal, one_minus_ltqnorm_arg

  !-----------------------------------------------------------------------

    !----- Begin Code -----
    l_error = .false.

    ! The calculation of PDF component cloud fraction (cloud_frac_i) in
    ! pdf_closure now sets cloud_frac_i to 0 when the special condition is met
    ! that both | mu_chi_i | <= eps and sigma_chi_i <= chi_tol.  This check
    ! should be omitted under such conditions.
    ! 
    ! Note:  The conditions on | mu_chi_i | and sigma_chi_i need to match the
    !        conditions given in subroutine calc_cloud_frac_component in
    !        pdf_closure_module.F90 (a chi_at_sat of 0 can be used).  If those
    !        conditions are changed in that subroutine, they need to change here
    !        too.  Someday, a better solution might be to pass the value of
    !        pdf_description_i out of that part of the code and into this check.
    if ( abs( mu_chi_i ) <= eps .and. sigma_chi_i <= chi_tol ) then
       ! Return without performing this check.
       return
    endif

    ! Check left end of box.
    one_minus_ltqnorm_arg = cloud_frac_i - box_half_width
    ! Do not bother to check this end of the box if it dips below the minimum
    ! ltqnorm argument.
    if ( one_minus_ltqnorm_arg >= ltqnorm_min_arg ) then

      if ( one_minus_ltqnorm_arg > (one - ltqnorm_min_arg) ) then
        one_minus_ltqnorm_arg = one - ltqnorm_min_arg
      end if

      chi_std_normal = ltqnorm( one - one_minus_ltqnorm_arg )
      chi = chi_std_normal * sigma_chi_i + mu_chi_i
      if ( chi <= zero ) then
        l_error = .true.
        write(fstderr,*) 'chi (left side of box) = ', chi
      end if
    end if ! one_minus_ltqnorm_arg >= ltqnorm_min_arg

    ! Check right end of box.
    one_minus_ltqnorm_arg = cloud_frac_i + box_half_width
    ! Do not bother to check this end of the box if it exceeds the maximum
    ! ltqnorm argument
    if ( one_minus_ltqnorm_arg <= one-ltqnorm_min_arg ) then

      if ( one_minus_ltqnorm_arg < ltqnorm_min_arg ) then
        one_minus_ltqnorm_arg = ltqnorm_min_arg
      end if

      chi_std_normal = ltqnorm( one - one_minus_ltqnorm_arg )
      chi = chi_std_normal * sigma_chi_i + mu_chi_i
      if ( chi > zero ) then
        l_error = .true.
        write(fstderr,*) 'chi (right side of box) = ', chi
      end if
    end if ! one_minus_ltqnorm_arg >= ltqnorm_min_arg

    if ( l_error ) then
      write(fstderr,*) "In assert_consistent_cf_component, cloud_frac_i is inconsistent with &
                       &mu_chi_i and stdev_chi_i."
      write(fstderr,*) "mu_chi_i = ", mu_chi_i
      write(fstderr,*) "sigma_chi_i = ", sigma_chi_i
      write(fstderr,*) "cloud_frac_i = ", cloud_frac_i
    end if

    return
  end subroutine assert_consistent_cf_component
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
  subroutine assert_correct_cloud_normal( num_samples, X_u_chi, X_nl_chi, X_mixt_comp, &
                                          cloud_frac_1, cloud_frac_2, &
                                          l_error )

  ! Description:
  !   Asserts that all SILHS sample points that are in cloud in uniform space
  !   are in cloud in normal space, and that all SILHS sample points that are
  !   in clear air in uniform space are in clear air in normal space.
  
  ! References:
  !   None
  !-----------------------------------------------------------------------
  
    ! Included Modules
    use clubb_precision, only: &
      core_rknd

    use constants_clubb, only: &
      one, &      ! Constant(s)
      fstderr
      
    use parameters_silhs, only: &
      single_prec_thresh   ! Constant

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      num_samples            ! Number of SILHS sample points

    real( kind = core_rknd ), dimension(num_samples), intent(in) :: &
      X_u_chi,  &            ! Samples of chi in uniform space
      X_nl_chi               ! Samples of chi in normal space

    integer, dimension(num_samples), intent(in) :: &
      X_mixt_comp            ! PDF component of each sample

    real( kind = core_rknd ), intent(in) :: &
      cloud_frac_1,   &      ! Cloud fraction in PDF component 1
      cloud_frac_2           ! Cloud fraction in PDF component 2

    ! Output Variables
    logical, intent(out) :: &
      l_error                ! True if the assertion check fails

    ! Local Variables
    real( kind = core_rknd ) :: &
      cloud_frac_i

    integer :: sample

  !-----------------------------------------------------------------------

    !----- Begin Code -----
    l_error = .false.

    do sample = 1, num_samples, 1

      ! Determine the appropriate cloud fraction
      if ( X_mixt_comp(sample) == 1 ) then
        cloud_frac_i = cloud_frac_1
      else if ( X_mixt_comp(sample) == 2 ) then
        cloud_frac_i = cloud_frac_2
      end if

      if ( X_u_chi(sample) < (one - cloud_frac_i) ) then

        ! The uniform sample is in clear air
        if ( X_nl_chi(sample) > 1000._core_rknd * single_prec_thresh ) then
          l_error = .true.
          write(fstderr,*) "X_nl_chi(", sample, ") > ", single_prec_thresh
        end if

      else if ( X_u_chi(sample) >= (one - cloud_frac_i) .and. &
                X_u_chi(sample) < one ) then

        ! The uniform sample is in cloud
        if ( X_nl_chi(sample) <= - 1000._core_rknd * single_prec_thresh ) then
          l_error = .true.
          write(fstderr,*) "X_nl_chi(", sample, ") <= ", -single_prec_thresh
        end if

      else
        stop "X_u_chi not in correct range in assert_correct_cloud_normal"
      end if

    end do ! 1..num_samples

    if ( l_error ) then
      write(fstderr,*) "In assert_correct_cloud_normal:"
      write(fstderr,*) "The 'cloudiness' of points in uniform and normal space is not consistent"
      write(fstderr,'(4X,A,A)')  "X_u_chi         ", "X_nl_chi "
      do sample = 1, num_samples, 1
        write(fstderr,'(I4,2G20.4)') &
          sample, X_u_chi(sample), X_nl_chi(sample)
      end do
      ! This will hopefully stop the run at some unknown point in the future
      l_error = .true.
    end if  ! in_cloud_points /= out_of_cloud_points

    return
  end subroutine assert_correct_cloud_normal
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
  subroutine latin_hypercube_2D_output &
             ( fname_prefix, fdir, stats_tout, nz, &
               stats_zt, time_initial, num_samples )
!-------------------------------------------------------------------------------

    use array_index, only: &
      iiPDF_chi, & ! Variables
      iiPDF_eta, &
      iiPDF_w, &
      iiPDF_rr, & 
      iiPDF_ri, &
      iiPDF_rs, &
      iiPDF_rg, &
      iiPDF_Nr, &
      iiPDF_Ni, &
      iiPDF_Ns, &
      iiPDF_Ng, &
      iiPDF_Ncn

    use clubb_precision, only: &
      time_precision, & ! Constant
      core_rknd

    use output_2D_samples_module, only: &
      open_2D_samples_file ! Procedure

    use output_2D_samples_module, only: &
      lognormal_sample_file, & ! Instance of a type
      uniform_sample_file

    use corr_varnce_module, only: &
      pdf_dim! Variable


    implicit none

    ! Input Variables
    character(len=*), intent(in) :: &
      fname_prefix, & ! Prefix for file name
      fdir            ! Directory for output

    real(kind=core_rknd), intent(in) :: &
      stats_tout    ! Frequency to write to disk        [s]

    real(kind=time_precision), intent(in) :: &
      time_initial  ! Initial time                      [s]

    integer, intent(in) :: &
      nz ! Number of vertical levels

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      stats_zt ! Altitudes [m]

    integer, intent(in) :: num_samples

    ! Local Variables
    character(len=100), allocatable, dimension(:) :: &
      variable_names, variable_descriptions, variable_units

    integer :: i

    ! ---- Begin Code ----

    if ( l_output_2D_lognormal_dist ) then

      allocate( variable_names(pdf_dim), variable_descriptions(pdf_dim), &
                variable_units(pdf_dim) )

      variable_names(iiPDF_chi)        = "chi"
      variable_descriptions(iiPDF_chi) = "The variable 's' from Mellor 1977"
      variable_units(iiPDF_chi)        = "kg/kg"

      variable_names(iiPDF_eta)        = "eta"
      variable_descriptions(iiPDF_eta) = "The variable 't' from Mellor 1977"
      variable_units(iiPDF_eta)        = "kg/kg"

      variable_names(iiPDF_w)        = "w"
      variable_descriptions(iiPDF_w) = "Vertical velocity"
      variable_units(iiPDF_w)        = "m/s"

      if ( iiPDF_rr > 0 ) then
        variable_names(iiPDF_rr)        = "rr"
        variable_descriptions(iiPDF_rr) = "Rain water mixing ratio"
        variable_units(iiPDF_rr)        = "kg/kg"
      end if
      if ( iiPDF_ri > 0 ) then
        variable_names(iiPDF_ri)        = "ri"
        variable_descriptions(iiPDF_ri) = "Ice water mixing ratio"
        variable_units(iiPDF_ri)        = "kg/kg"
      end if
      if ( iiPDF_rs > 0 ) then
        variable_names(iiPDF_rs)        = "rs"
        variable_descriptions(iiPDF_rs) = "Snow water mixing ratio"
        variable_units(iiPDF_rs)        = "kg/kg"
      end if
      if ( iiPDF_rg > 0 ) then
        variable_names(iiPDF_rg)        = "rg"
        variable_descriptions(iiPDF_rg) = "Graupel water mixing ratio"
        variable_units(iiPDF_rg)        = "kg/kg"
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
      if ( iiPDF_Ns > 0 ) then
        variable_names(iiPDF_Ns)        = "Ns"
        variable_descriptions(iiPDF_Ns) = "Snow number concentration"
        variable_units(iiPDF_Ns)        = "count/kg"
      end if
      if ( iiPDF_Ng > 0 ) then
        variable_names(iiPDF_Ng)        = "Ng"
        variable_descriptions(iiPDF_Ng) = "Graupel number concentration"
        variable_units(iiPDF_Ng)        = "count/kg"
      end if

      call open_2D_samples_file( nz, num_samples, pdf_dim, & ! In
                                 trim( fname_prefix )//"_nl", fdir, & ! In
                                 time_initial, stats_tout, stats_zt, variable_names, & ! In
                                 variable_descriptions, variable_units, & ! In
                                 lognormal_sample_file ) ! In/Out

      deallocate( variable_names, variable_descriptions, variable_units )

    end if

    if ( l_output_2D_uniform_dist ) then

      allocate( variable_names(pdf_dim+4), variable_descriptions(pdf_dim+4), &
                variable_units(pdf_dim+4) )

      ! The uniform distribution corresponds to all the same variables as X_nl,
      ! except the d+1 component is the mixture component.

      variable_names(iiPDF_chi)        = "chi"
      variable_descriptions(iiPDF_chi) = "Uniform dist of the variable 's' from Mellor 1977"

      variable_names(iiPDF_eta)        = "eta"
      variable_descriptions(iiPDF_eta) = "Uniform dist of the variable 't' from Mellor 1977"

      variable_names(iiPDF_w)        = "w"
      variable_descriptions(iiPDF_w) = "Uniform dist of the vertical velocity"


      if ( iiPDF_rr > 0 ) then
        variable_names(iiPDF_rr)        = "rr"
        variable_descriptions(iiPDF_rr) = "Rain water mixing ratio"
        variable_units(iiPDF_rr)        = "kg/kg"
      end if
      if ( iiPDF_ri > 0 ) then
        variable_names(iiPDF_ri)        = "ri"
        variable_descriptions(iiPDF_ri) = "Ice water mixing ratio"
        variable_units(iiPDF_ri)        = "kg/kg"
      end if
      if ( iiPDF_rs > 0 ) then
        variable_names(iiPDF_rs)        = "rs"
        variable_descriptions(iiPDF_rs) = "Snow water mixing ratio"
        variable_units(iiPDF_rs)        = "kg/kg"
      end if
      if ( iiPDF_rg > 0 ) then
        variable_names(iiPDF_rg)        = "rg"
        variable_descriptions(iiPDF_rg) = "Graupel water mixing ratio"
        variable_units(iiPDF_rg)        = "kg/kg"
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
      if ( iiPDF_Ns > 0 ) then
        variable_names(iiPDF_Ns)        = "Ns"
        variable_descriptions(iiPDF_Ns) = "Snow number concentration"
        variable_units(iiPDF_Ns)        = "count/kg"
      end if
      if ( iiPDF_Ng > 0 ) then
        variable_names(iiPDF_Ng)        = "Ng"
        variable_descriptions(iiPDF_Ng) = "Graupel number concentration"
        variable_units(iiPDF_Ng)        = "count/kg"
      end if

      i = pdf_dim+ 1
      variable_names(i) = "dp1"
      variable_descriptions(i) = "Uniform distribution for the mixture component"

      i = pdf_dim+ 2
      variable_names(i) = "dp2"
      variable_descriptions(i) = "Uniform variate used to determine precipitation!"

      i = pdf_dim+ 3
      variable_names(i) = "X_mixt_comp"
      variable_descriptions(i) = "Mixture component (should be 1 or 2)"

      i = pdf_dim+ 4
      variable_names(i) = "lh_sample_point_weights"
      variable_descriptions(i) = "Weight of each sample point"

      ! Set all the units
      variable_units(:) = "count" ! Unidata units format for a dimensionless quantity

      call open_2D_samples_file( nz, num_samples, i, & ! In
                                 trim( fname_prefix )//"_u", fdir, & ! In
                                 time_initial, stats_tout, stats_zt, & ! In
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
  subroutine compute_arb_overlap( nz, num_samples, pdf_dim, d_uniform_extra, &
                                  k_lh_start, vert_corr, rand_pool, &
                                  X_u_all_levs )
! Description:
!   Re-computes X_u (uniform sample) using an arbitrary correlation specified 
!   by X_vert_corr (which can vary with height).
!   This is an improved algorithm that doesn't require us to convert from a
!   unifrom distribution to a Gaussian distribution and back again.
!
! References:
!   None
!-------------------------------------------------------------------------------

    use generate_uniform_sample_module, only: &
      rand_uniform_real ! Procedure

    use clubb_precision, only: &
      core_rknd ! Precision

    use constants_clubb, only: &
      zero, &   ! Constants
      one, &
      two, &
      fstderr
      
    use parameters_silhs, only: &
      single_prec_thresh

    implicit none

    ! ---------------- Input Variables ----------------
    integer, intent(in) :: &
      nz,               & ! Vertical levels
      pdf_dim,          & ! Number of variates in CLUBB's PDF
      d_uniform_extra,  & ! Uniform variates included in uniform sample but not
                                  !  in normal/lognormal sample
      num_samples,       & ! Number of SILHS sample points
      k_lh_start

    real(kind=core_rknd), dimension(nz), intent(in) :: &
      vert_corr ! Vertical correlation between k points in range [0,1]   [-]
      
    real(kind = core_rknd), dimension(nz,num_samples,pdf_dim+d_uniform_extra) :: &
      rand_pool ! Array of randomly generated numbers

    ! ---------------- Output Variables ----------------
    real(kind=core_rknd), dimension(nz,num_samples,pdf_dim+d_uniform_extra), intent(inout) :: &
      X_u_all_levs ! Uniform distribution of 1 variate at all levels [-]
                   ! The value of this variate at k_lh_start should already be populated
                   ! in this array and will be used to fill in the other levels.

    ! ---------------- Local Variables ----------------
    real(kind=core_rknd) :: min_val, half_width, offset, unbounded_point

    integer :: k, sample, i ! Loop iterators

    ! ---------------- Begin Code ----------------
    
    !$acc wait(1) async(2)

    ! Recompute from k_lh_start to nz-1 for all samples and variates, upward loop
    !$acc parallel loop collapse(2) default(present) async(1)
    do i = 1, pdf_dim + d_uniform_extra                         
      do sample = 1, num_samples
        
        unbounded_point = X_u_all_levs(k_lh_start,sample,i)
        
         do k = k_lh_start, nz-1
           
           half_width = one - vert_corr(k+1)
           min_val = unbounded_point - half_width

           offset = two * half_width * rand_pool(k+1,sample,i)

           unbounded_point = min_val + offset
           
           ! If unbounded_point lies outside [single_prec_thresh,1-single_prec_thresh],
           ! fold it back so that it is between the valid range
           if ( unbounded_point > one - single_prec_thresh ) then
             unbounded_point = two - unbounded_point - two * single_prec_thresh
           else if ( unbounded_point < single_prec_thresh ) then
             unbounded_point = - unbounded_point + two * single_prec_thresh
           end if
           
           X_u_all_levs(k+1,sample,i) = unbounded_point
           
         end do ! k_lh_start..nz-1
      end do ! 1..num_samples
    end do ! 1..pdf_dim


    ! Recompute from k_lh_start down to 2 for all samples and variates, downward loop 
    !$acc parallel loop collapse(2) default(present) async(2)
    do i = 1, pdf_dim + d_uniform_extra                        
      do sample = 1, num_samples
        
        unbounded_point = X_u_all_levs(k_lh_start,sample,i)
        
         do k = k_lh_start, 2, -1

           half_width = one - vert_corr(k-1)
           min_val = unbounded_point - half_width

           offset = two * half_width * rand_pool(k-1,sample,i)

           unbounded_point = min_val + offset
           
           ! If unbounded_point lies outside [single_prec_thresh,1-single_prec_thresh],
           ! fold it back so that it is between the valid range 
           if ( unbounded_point > one - single_prec_thresh ) then
             unbounded_point = two - unbounded_point - two * single_prec_thresh
           else if ( unbounded_point < single_prec_thresh ) then
             unbounded_point = - unbounded_point + two * single_prec_thresh
           end if
           
           X_u_all_levs(k-1,sample,i) = unbounded_point

         end do ! k_lh_start..2 decrementing
      end do ! 1..num_samples
    end do ! 1..pdf_dim
    
    !$acc wait(2) async(1)

    return
  end subroutine compute_arb_overlap

!-------------------------------------------------------------------------------
  subroutine stats_accumulate_lh &
             ( nz, num_samples, pdf_dim, rho_ds_zt, &
               lh_sample_point_weights, X_nl_all_levs, &
               lh_rt_clipped, lh_thl_clipped, & 
               lh_rc_clipped, lh_rv_clipped, & 
               lh_Nc_clipped )

! Description:
!   Clip subcolumns from latin hypercube and create stats for diagnostic
!   purposes.

! References:
!   None
!-------------------------------------------------------------------------------

    use parameters_model, only: hydromet_dim ! Variable

    use grid_class, only: gr

    use stats_variables, only: &
      l_stats_samp, & ! Variable(s)
      ilh_rrm, &
      ilh_Nrm, &
      ilh_rim, &
      ilh_Nim, &
      ilh_rsm, &
      ilh_Nsm, &
      ilh_rgm, &
      ilh_Ngm, &
      ilh_thlm, &
      ilh_rcm, &
      ilh_Ncm, &
      ilh_Ncnm, &
      ilh_rvm, &
      ilh_wm, &
      ilh_cloud_frac, &
      ilh_cloud_frac_unweighted, &
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
      ilh_rrp2_zt, &
      ilh_vwp, &
      ilh_lwp, &
      ilh_sample_weights_sum, &
      ilh_sample_weights_avg, &
      stats_lh_zt, &
      stats_lh_sfc

    use math_utilities, only: &
      compute_sample_mean, & ! Procedure(s)
      compute_sample_variance

    use stats_type_utilities, only: &
      stat_update_var, & ! Procedure(s)
      stat_update_var_pt

    use array_index, only: &
      iirr, & ! Variables(s)
      iirs, & 
      iiri, & 
      iirg, & 
      iiNr, &
      iiNs, &
      iiNi, &
      iiNg, &
      iiPDF_chi, &
      iiPDF_eta, &
      iiPDF_w,   &
      iiPDF_Ncn

    use constants_clubb, only: & 
      zero, &            ! Constant(s)
      one

    use clubb_precision, only: & 
      core_rknd    ! Constant

   use fill_holes, only: &
     vertical_integral ! Procedure(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      pdf_dim,     & ! Number of variables to sample
      num_samples,   & ! Number of calls to microphysics per timestep (normally=2)
      nz                 ! Number of vertical model levels

    real( kind = core_rknd ), intent(in), dimension(nz) :: &
      rho_ds_zt  ! Dry, static density (thermo. levs.) [kg/m^3]

    real( kind = core_rknd ), intent(in), dimension(nz,num_samples) :: &
      lh_sample_point_weights

    real( kind = core_rknd ), intent(in), dimension(nz,num_samples,pdf_dim) :: &
      X_nl_all_levs ! Sample that is transformed ultimately to normal-lognormal
      
    real( kind = core_rknd ), dimension(nz,num_samples), intent(in) :: &
      lh_rt_clipped,  & ! rt generated from silhs sample points
      lh_thl_clipped, & ! thl generated from silhs sample points
      lh_rc_clipped,  & ! rc generated from silhs sample points
      lh_rv_clipped,  & ! rv generated from silhs sample points
      lh_Nc_clipped     ! Nc generated from silhs sample points

    ! Local variables
    real( kind = core_rknd ), dimension(nz,num_samples) :: &
      Ncn_all_points ! Cloud nuclei conc. for all levs.; Nc=Ncn*H(chi) [#/kg]

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
      lh_rrp2_zt, & ! Average value of the variance of the LH est. of rr.         [(kg/kg)^2]
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

    integer :: sample, ivar

    ! ---- Begin Code ----

    if ( l_stats_samp ) then

      ! For all cases where l_lh_importance_sampling is false, the weights
      ! will be 1 (all points equally weighted)

      if ( ilh_rcm + ilh_rcp2_zt + ilh_lwp > 0 ) then
        lh_rcm = compute_sample_mean( nz, num_samples, lh_sample_point_weights, &
                                      lh_rc_clipped )
        call stat_update_var( ilh_rcm, lh_rcm, stats_lh_zt )

        if ( ilh_lwp > 0 ) then
          xtmp &
          = vertical_integral &
               ( (nz - 2 + 1), rho_ds_zt(2:nz), &
                 lh_rcm(2:nz), gr%dzt(2:nz) )

          call stat_update_var_pt( ilh_lwp, 1, xtmp, stats_lh_sfc )
        end if
      end if

      if ( ilh_sample_weights_sum > 0 ) then
          xtmp = sum(lh_sample_point_weights(:,:))
          call stat_update_var_pt( ilh_sample_weights_sum, 1, xtmp, stats_lh_sfc )
      end if
      
      if ( ilh_sample_weights_avg > 0 ) then
          xtmp = sum(lh_sample_point_weights(:,:)) / real( num_samples*nz, kind = core_rknd )
          call stat_update_var_pt( ilh_sample_weights_avg, 1, xtmp, stats_lh_sfc )
      end if
        
      if ( ilh_thlm + ilh_thlp2_zt > 0 ) then
        lh_thlm = compute_sample_mean( nz, num_samples, lh_sample_point_weights, &
                                       real( lh_thl_clipped, kind = core_rknd ) )
        call stat_update_var( ilh_thlm, lh_thlm, stats_lh_zt )
      end if

      if ( ilh_rvm + ilh_rtp2_zt > 0 ) then
        lh_rvm = compute_sample_mean( nz, num_samples, lh_sample_point_weights, &
                                      lh_rv_clipped )
        call stat_update_var( ilh_rvm, lh_rvm, stats_lh_zt )
        if ( ilh_vwp > 0 ) then
          xtmp &
          = vertical_integral &
               ( (nz - 2 + 1), rho_ds_zt(2:nz), &
                 lh_rvm(2:nz), gr%dzt(2:nz) )

          call stat_update_var_pt( ilh_vwp, 1, xtmp, stats_lh_sfc )
        end if
      end if

      if ( ilh_wm + ilh_wp2_zt > 0 ) then
        lh_wm  = compute_sample_mean( nz, num_samples, lh_sample_point_weights, &
                                      real( X_nl_all_levs(:,:,iiPDF_w), kind = core_rknd) )
        call stat_update_var( ilh_wm, lh_wm, stats_lh_zt )
      end if

      if ( ilh_rrm + ilh_Nrm + ilh_rim + ilh_Nim + ilh_rsm + ilh_Nsm + &
           ilh_rgm + ilh_Ngm + ilh_Ncnm + ilh_Ncm > 0 ) then

        lh_hydromet = 0._core_rknd
        call copy_X_nl_into_hydromet_all_pts( nz, pdf_dim, num_samples, & ! In
                                      X_nl_all_levs, &  ! In
                                      lh_hydromet, & ! In
                                      hydromet_all_points, &  ! Out
                                      Ncn_all_points ) ! Out

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
        call stat_update_var( ilh_Ncnm, lh_Ncnm, stats_lh_zt )
      end if

      if ( ilh_Ncm > 0 ) then
        lh_Ncm = compute_sample_mean( nz, num_samples, lh_sample_point_weights, &
                                      lh_Nc_clipped(:,:) )
        call stat_update_var( ilh_Ncm, lh_Ncm, stats_lh_zt )
      end if

      ! Latin hypercube estimate of cloud fraction
      if ( ilh_cloud_frac > 0 ) then
        lh_cloud_frac(:) = zero
        do sample = 1, num_samples
          where ( X_nl_all_levs(:,sample,iiPDF_chi) > zero )
            lh_cloud_frac(:) = lh_cloud_frac(:) + one * lh_sample_point_weights(:,sample)
          end where
        end do
        lh_cloud_frac(:) = lh_cloud_frac(:) / real( num_samples, kind = core_rknd )

        call stat_update_var( ilh_cloud_frac, lh_cloud_frac, stats_lh_zt )
      end if

      ! Sample of lh_cloud_frac that is not weighted
      if ( ilh_cloud_frac_unweighted > 0 ) then
        lh_cloud_frac(:) = zero
        do sample = 1, num_samples
          where ( X_nl_all_levs(:,sample,iiPDF_chi) > zero )
            lh_cloud_frac(:) = lh_cloud_frac(:) + one
          end where
        end do
        lh_cloud_frac(:) = lh_cloud_frac(:) / real( num_samples, kind = core_rknd )

        call stat_update_var( ilh_cloud_frac_unweighted, lh_cloud_frac, stats_lh_zt )
      end if

      ! Latin hypercube estimate of chi
      if ( ilh_chi > 0 ) then
        lh_chi(1:nz) &
        = compute_sample_mean( nz, num_samples, lh_sample_point_weights, &
                               X_nl_all_levs(1:nz, 1:num_samples, iiPDF_chi) )
        call stat_update_var( ilh_chi, lh_chi, stats_lh_zt )
      end if

      ! Latin hypercube estimate of variance of chi
      if ( ilh_chip2 > 0 ) then
        lh_chip2(1:nz) &
        = compute_sample_variance( nz, num_samples, &
                                   X_nl_all_levs(:,:,iiPDF_chi), &
                                   lh_sample_point_weights, lh_chi(1:nz) )
        call stat_update_var( ilh_chip2, lh_chip2, stats_lh_zt )
      end if

      ! Latin hypercube estimate of eta
      if ( ilh_eta > 0 ) then
        lh_eta(1:nz) &
        = compute_sample_mean( nz, num_samples, lh_sample_point_weights, &
                               X_nl_all_levs(1:nz, 1:num_samples, iiPDF_eta) )

        call stat_update_var( ilh_eta, lh_eta, stats_lh_zt )
      end if

      if ( ilh_wp2_zt > 0 ) then
        ! Compute the variance of vertical velocity
        lh_wp2_zt = compute_sample_variance( nz, num_samples, &
                                             X_nl_all_levs(:,:,iiPDF_w), &
                                             lh_sample_point_weights, lh_wm )
        call stat_update_var( ilh_wp2_zt, lh_wp2_zt, stats_lh_zt )
      end if

      if ( ilh_rcp2_zt  > 0 ) then
        ! Compute the variance of cloud water mixing ratio
        lh_rcp2_zt = compute_sample_variance &
                     ( nz, num_samples, lh_rc_clipped, &
                       lh_sample_point_weights, lh_rcm )
        call stat_update_var( ilh_rcp2_zt, lh_rcp2_zt, stats_lh_zt )
      end if

      if ( ilh_rtp2_zt > 0 ) then
        ! Compute the variance of total water
        lh_rtp2_zt = compute_sample_variance &
                     ( nz, num_samples, &
                       lh_rt_clipped, lh_sample_point_weights, &
                       lh_rvm+lh_rcm )
        call stat_update_var( ilh_rtp2_zt, lh_rtp2_zt, stats_lh_zt )
      end if

      if ( ilh_thlp2_zt > 0 ) then
        ! Compute the variance of liquid potential temperature
        lh_thlp2_zt = compute_sample_variance( nz, num_samples, &
                        lh_thl_clipped, lh_sample_point_weights, &
                        lh_thlm )
        call stat_update_var( ilh_thlp2_zt, lh_thlp2_zt, stats_lh_zt )
      end if

      ! Compute the variance of rain water mixing ratio
      if ( iirr > 0 .and. ilh_rrp2_zt > 0 ) then
        lh_rrp2_zt = compute_sample_variance &
                        ( nz, num_samples, hydromet_all_points(:,:,iirr), &
                          lh_sample_point_weights, lh_hydromet(:,iirr) )
        call stat_update_var( ilh_rrp2_zt, lh_rrp2_zt, stats_lh_zt )
      end if

      ! Compute the variance of cloud nuclei concentration (simplifed)
      if ( iiPDF_Ncn > 0 .and. ilh_Ncnp2_zt > 0 ) then
        lh_Ncnp2_zt = compute_sample_variance &
                      ( nz, num_samples, Ncn_all_points(:,:), &
                        lh_sample_point_weights, lh_Ncnm(:) )
        call stat_update_var( ilh_Ncnp2_zt, lh_Ncnp2_zt, stats_lh_zt )
      end if

      ! Compute the variance of cloud droplet concentration
      if ( ilh_Ncp2_zt > 0 ) then
        lh_Ncp2_zt = compute_sample_variance &
                     ( nz, num_samples, lh_Nc_clipped(:,:), &
                       lh_sample_point_weights, lh_Ncm(:) )
        call stat_update_var( ilh_Ncp2_zt, lh_Ncp2_zt, stats_lh_zt )
      end if

      ! Compute the variance of rain droplet number concentration
      if ( iiNr > 0 .and. ilh_Nrp2_zt > 0 ) then
        lh_Nrp2_zt = compute_sample_variance( nz, num_samples, hydromet_all_points(:,:,iiNr),&
                                              lh_sample_point_weights, lh_hydromet(:,iiNr) )
        call stat_update_var( ilh_Nrp2_zt, lh_Nrp2_zt, stats_lh_zt )
      end if

      ! Averages of points being fed into the microphysics
      ! These are for diagnostic purposes, and are not needed for anything
      if ( iirr > 0 ) then
        call stat_update_var( ilh_rrm, lh_hydromet(:,iirr), stats_lh_zt )
      end if
      if ( iiNr > 0 ) then
        call stat_update_var( ilh_Nrm, lh_hydromet(:,iiNr), stats_lh_zt )
      end if
      if ( iiri > 0 ) then
        call stat_update_var( ilh_rim, lh_hydromet(:,iiri), stats_lh_zt )
      end if
      if ( iiNi > 0 ) then
        call stat_update_var( ilh_Nim, lh_hydromet(:,iiNi), stats_lh_zt )
      end if
      if ( iirs > 0 ) then
        call stat_update_var( ilh_rsm, lh_hydromet(:,iirs), stats_lh_zt )
      end if
      if ( iiNs > 0 ) then
        call stat_update_var( ilh_Nsm, lh_hydromet(:,iiNs), stats_lh_zt )
      end if
      if ( iirg > 0 ) then
        call stat_update_var( ilh_rgm, lh_hydromet(:,iirg), stats_lh_zt )
      end if
      if ( iiNg > 0 ) then
        call stat_update_var( ilh_Ngm, lh_hydromet(:,iiNg), stats_lh_zt )
      end if

    end if ! l_stats_samp

    return
  end subroutine stats_accumulate_lh

  !-----------------------------------------------------------------------
  subroutine stats_accumulate_uniform_lh( nz, num_samples, l_in_precip_all_levs, &
                                          X_mixt_comp_all_levs, X_u_chi_all_levs, pdf_params, &
                                          lh_sample_point_weights, k_lh_start )

  ! Description:
  !   Samples statistics that cannot be deduced from the normal-lognormal
  !   SILHS sample (X_nl_all_levs)

  ! References:
  !   None
  !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd            ! Constant

    use stats_type_utilities, only: &
      stat_update_var,   & ! Procedure(s)
      stat_update_var_pt

    use stats_variables, only: &
      ilh_precip_frac, &  ! Variable(s)
      ilh_mixt_frac, &
      ilh_precip_frac_unweighted, &
      ilh_mixt_frac_unweighted, &
      ik_lh_start, &
      ilh_samp_frac_category, &
      stats_lh_zt, &
      stats_lh_sfc

    use math_utilities, only: &
      compute_sample_mean ! Procedure

    use constants_clubb, only: &
      one, &              ! Constant(s)
      zero

    use pdf_parameter_module, only: &
      pdf_parameter       ! Type

    use silhs_importance_sample_module, only: &
      importance_category_type, &  ! Type
      num_importance_categories, & ! Constant
      define_importance_categories

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

    real( kind = core_rknd ), dimension(nz,num_samples), intent(in) :: &
      X_u_chi_all_levs     ! Uniform value of chi

    type(pdf_parameter), intent(in) :: &
      pdf_params           ! The official PDF parameters!

    real( kind = core_rknd ), dimension(nz,num_samples), intent(in) :: &
      lh_sample_point_weights ! The weight of each sample

    integer, intent(in) :: &
      k_lh_start           ! Vertical level for sampling preferentially within       [-]
                           ! cloud

    ! Local Variables
    type(importance_category_type), dimension(num_importance_categories) :: &
      importance_categories

    real( kind = core_rknd ), dimension(nz) :: &
      lh_precip_frac, &
      lh_mixt_frac

    real( kind = core_rknd ), dimension(nz,num_samples) :: &
      int_in_precip, & ! '1' for samples in precipitation, '0' otherwise
      int_mixt_comp    ! '1' for samples in the first PDF component, '0' otherwise

    real( kind = core_rknd ), dimension(num_samples) :: &
      one_weights

    real( kind = core_rknd ), dimension(nz,num_importance_categories) :: &
      lh_samp_frac

    real( kind = core_rknd ) :: &
      cloud_frac_i

    logical :: &
      l_in_cloud, &
      l_in_comp_1

    integer, dimension(num_importance_categories) :: &
      category_counts  ! Count of number of samples in each importance category

    integer :: k, isample, icategory
  !-----------------------------------------------------------------------

    !----- Begin Code -----
    
    ! Estimate of lh_precip_frac
    if ( ilh_precip_frac > 0 ) then
      where ( l_in_precip_all_levs )
        int_in_precip = 1.0_core_rknd
      else where
        int_in_precip = 0.0_core_rknd
      end where
      lh_precip_frac(:) = compute_sample_mean( nz, num_samples, lh_sample_point_weights, &
                                               int_in_precip )
      call stat_update_var( ilh_precip_frac, lh_precip_frac, stats_lh_zt )
    end if

    ! Unweighted estimate of lh_precip_frac
    if ( ilh_precip_frac_unweighted > 0 ) then
      where ( l_in_precip_all_levs )
        int_in_precip = 1.0_core_rknd
      else where
        int_in_precip = 0.0_core_rknd
      end where
      one_weights = one
      lh_precip_frac(:) = compute_sample_mean( nz, num_samples, one_weights, &
                                               int_in_precip )
      call stat_update_var( ilh_precip_frac_unweighted, lh_precip_frac, stats_lh_zt )
    end if

    ! Estimate of lh_mixt_frac
    if ( ilh_mixt_frac > 0 ) then
      where ( X_mixt_comp_all_levs == 1 )
        int_mixt_comp = 1.0_core_rknd
      else where
        int_mixt_comp = 0.0_core_rknd
      end where
      lh_mixt_frac(:) = compute_sample_mean( nz, num_samples, lh_sample_point_weights, &
                                             int_mixt_comp )
      call stat_update_var( ilh_mixt_frac, lh_mixt_frac, stats_lh_zt )
    end if

    ! Unweighted estimate of lh_mixt_frac
    if ( ilh_mixt_frac_unweighted > 0 ) then
      where ( X_mixt_comp_all_levs == 1 )
        int_mixt_comp = 1.0_core_rknd
      else where
        int_mixt_comp = 0.0_core_rknd
      end where
      one_weights = one
      lh_mixt_frac(:) = compute_sample_mean( nz, num_samples, one_weights, &
                                             int_mixt_comp )
      call stat_update_var( ilh_mixt_frac_unweighted, lh_mixt_frac, stats_lh_zt )
    end if

    ! k_lh_start is an integer, so it would be more appropriate to sample it
    ! as an integer, but as far as I can tell our current sampling
    ! infrastructure mainly supports sampling real numbers.
    call stat_update_var_pt( ik_lh_start, 1, real( k_lh_start, kind=core_rknd ), stats_lh_sfc )

    if ( allocated( ilh_samp_frac_category ) ) then
      if ( ilh_samp_frac_category(1) > 0 ) then

        importance_categories = define_importance_categories( )

        do k=1, nz
          category_counts(:) = 0

          do isample=1, num_samples

            if ( X_mixt_comp_all_levs(k,isample) == 1 ) then
              l_in_comp_1 = .true.
              cloud_frac_i = pdf_params%cloud_frac_1(k)
            else
              l_in_comp_1 = .false.
              cloud_frac_i = pdf_params%cloud_frac_2(k)
            end if

            l_in_cloud = X_u_chi_all_levs(k,isample) > (one - cloud_frac_i)

            do icategory=1, num_importance_categories
              if ( (l_in_cloud .eqv. importance_categories(icategory)%l_in_cloud) .and. &
                   (l_in_precip_all_levs(k,isample) .eqv. importance_categories(icategory)%&
                                                         l_in_precip) .and. &
                   (l_in_comp_1 .eqv. importance_categories(icategory)%l_in_component_1) ) then

                category_counts(icategory) = category_counts(icategory) + 1
                exit

              end if
            end do

          end do ! isample=1, num_samples

          do icategory=1, num_importance_categories
            lh_samp_frac(k,icategory) = real( category_counts(icategory), kind=core_rknd ) / &
                                        real( num_samples, kind=core_rknd )
          end do

        end do ! k=2, nz

        ! Microphysics is not run at lower level
        lh_samp_frac(1,:) = zero

        do icategory=1, num_importance_categories
          call stat_update_var( ilh_samp_frac_category(icategory), lh_samp_frac(:,icategory), &
                                stats_lh_zt )
        end do ! icategory=1, num_importance_categories

      end if ! ilh_samp_frac_category(1) > 0
    end if ! allocated( ilh_samp_frac_category )

    return
  end subroutine stats_accumulate_uniform_lh

  !-----------------------------------------------------------------------------
  subroutine copy_X_nl_into_hydromet_all_pts( nz, pdf_dim, num_samples, &
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
      iirr, & ! Variables
      iirs, & 
      iiri, & 
      iirg, & 
      iiNr, &
      iiNs, &
      iiNi, &
      iiNg, &
      iiPDF_rr, &
      iiPDF_rs, &
      iiPDF_ri, &
      iiPDF_rg, &
      iiPDF_Nr, &
      iiPDF_Ns, &
      iiPDF_Ng, &
      iiPDF_Ncn, &
      iiPDF_Ni

    use clubb_precision, only: &
      core_rknd

    implicit none

    integer, intent(in) :: &
      nz,            & ! Number of vertical levels
      pdf_dim,   & ! Number of variates
      num_samples    ! Number of calls to microphysics

    real( kind = core_rknd ), dimension(nz,num_samples,pdf_dim), intent(in) :: &
      X_nl_all_levs ! Sample that is transformed ultimately to normal-lognormal

    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(in) :: &
      hydromet ! Hydrometeor species    [units vary]

    real( kind = core_rknd ), dimension(nz,num_samples,hydromet_dim), intent(out) :: &
      hydromet_all_points ! Hydrometeor species    [units vary]

    real( kind = core_rknd ), dimension(nz,num_samples), intent(out) :: &
      Ncn_all_points    ! Cloud nuclei conc. (simplified); Nc=Ncn*H(chi)  [#/kg]

    integer :: sample, ivar

    do sample = 1, num_samples
      ! Copy the sample points into the temporary arrays
      do ivar = 1, hydromet_dim, 1
        if ( ivar == iirr .and. iiPDF_rr > 0 ) then
          ! Use a sampled value of rain water mixing ratio
          hydromet_all_points(:,sample,ivar) = &
            real( X_nl_all_levs(:,sample,iiPDF_rr), kind = core_rknd )

        else if ( ivar == iirs .and. iiPDF_rs > 0 ) then
          ! Use a sampled value of rain water mixing ratio
          hydromet_all_points(:,sample,ivar) = &
            real( X_nl_all_levs(:,sample,iiPDF_rs), kind = core_rknd )

        else if ( ivar == iiri .and. iiPDF_ri > 0 ) then
          ! Use a sampled value of rain water mixing ratio
          hydromet_all_points(:,sample,ivar) = &
            real( X_nl_all_levs(:,sample,iiPDF_ri), kind = core_rknd )

        else if ( ivar == iirg .and. iiPDF_rg > 0 ) then
          ! Use a sampled value of rain water mixing ratio
          hydromet_all_points(:,sample,ivar) = &
            real( X_nl_all_levs(:,sample,iiPDF_rg), kind = core_rknd )

        else if ( ivar == iiNr .and. iiPDF_Nr > 0 ) then
          ! Use a sampled value of rain droplet number concentration
          hydromet_all_points(:,sample,ivar) = &
            real( X_nl_all_levs(:,sample,iiPDF_Nr), kind = core_rknd )

        else if ( ivar == iiNs .and. iiPDF_Ns > 0 ) then
          ! Use a sampled value of rain droplet number concentration
          hydromet_all_points(:,sample,ivar) = &
            real( X_nl_all_levs(:,sample,iiPDF_Ns), kind = core_rknd )

        else if ( ivar == iiNg .and. iiPDF_Ng > 0 ) then
          ! Use a sampled value of rain droplet number concentration
          hydromet_all_points(:,sample,ivar) = &
            real( X_nl_all_levs(:,sample,iiPDF_Ng), kind = core_rknd )

        else if ( ivar == iiNi .and. iiPDF_Ni > 0 ) then
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

#endif /* SILHS */

end module latin_hypercube_driver_module
