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
  public :: latin_hypercube_2D_output_api, &
    latin_hypercube_2D_close_api, stats_accumulate_lh_api, generate_silhs_sample, &
    copy_X_nl_into_hydromet_all_pts, clip_transform_silhs_output

  private :: stats_accumulate_uniform_lh

  contains

!-------------------------------------------------------------------------------
  subroutine generate_silhs_sample( &
               iter, pdf_dim, num_samples, sequence_length, nzt, ngrdcol, & ! intent(in)
               l_calc_weights_all_levs_itime, &                            ! intent(in)
               gr, pdf_params, delta_zm, Lscale, &                         ! intent(in)
               lh_seed, hm_metadata, &                                     ! intent(in)
               !rho_ds_zt, &
               mu1, mu2, sigma1, sigma2, &                                 ! intent(in)
               corr_cholesky_mtx_1, corr_cholesky_mtx_2, &                 ! intent(in)
               precip_fracs, silhs_config_flags, &                         ! intent(in)
               vert_decorr_coef, &                                         ! intent(in)
               stats_metadata, &                                           ! intent(in)
               stats_lh_zt, stats_lh_sfc, err_info, &                      ! intent(inout)
               X_nl_all_levs, X_mixt_comp_all_levs, &                      ! intent(out)
               lh_sample_point_weights )                                   ! intent(out)

! Description:
!   Generate sample points of moisture, temperature, et cetera for the purpose
!   of computing tendencies with a microphysics or radiation scheme.
!
! References:
! https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:overview_silhs
!-------------------------------------------------------------------------------

    use grid_class, only: &
        grid    ! Type(s)

    use transform_to_pdf_module, only: &
      transform_uniform_samples_to_pdf      ! Procedure

    use output_2D_samples_module, only: &
      output_2D_lognormal_dist_file, & ! Procedure(s)
      output_2D_uniform_dist_file

    use pdf_parameter_module, only: &
      pdf_parameter  ! Type

    use pdf_utilities, only: &
      compute_mean_binormal  ! Procedure

    use hydromet_pdf_parameter_module, only: &
      precipitation_fractions ! Type

    use constants_clubb, only: &
      fstderr, & ! Constant(s)
      zero, &
      one, &
      rc_tol

    use clubb_precision, only: &
      core_rknd, &
      stat_rknd

    use parameters_silhs, only: &
      silhs_config_flags_type ! Type(s)

    use error_code, only: &
      clubb_at_least_debug_level  ! Procedure
      
    use advance_helper_module, only: &
      vertical_avg  ! Procedure
      
    use stats_variables, only: &
      stats_metadata_type

    use mt95, only: &
      genrand_intg  ! Type

    use stats_type, only: &
      stats ! Type 

    use corr_varnce_module, only: &
      hm_metadata_type

    use err_info_type_module, only: &
      err_info_type         ! Type

    use error_code, only: &
        clubb_fatal_error   ! Constant

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
      nzt,             & ! Number of thermodynamic vertical model levels
      ngrdcol            ! Number of grid columns

    type(grid), intent(in) :: &
      gr    ! Grid variable type

    type(pdf_parameter), intent(in) :: &
      pdf_params ! PDF parameters       [units vary]

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      delta_zm   ! Difference in momentum altitudes    [m]

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      Lscale       ! Turbulent mixing length            [m]
      
    integer( kind = genrand_intg ), intent(in) :: &
      lh_seed      ! Random number generator seed

    type (hm_metadata_type), intent(in) :: &
      hm_metadata

!   real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
!     rho_ds_zt    ! Dry, static density on thermo. levels    [kg/m^3]

    logical, intent(in) :: &
      l_calc_weights_all_levs_itime ! determines if vertically correlated sample points are needed
      
    ! More Input Variables!
    real( kind = core_rknd ), dimension(ngrdcol,nzt,pdf_dim,pdf_dim), intent(in) :: &
      corr_cholesky_mtx_1, & ! Correlations Cholesky matrix (1st comp.)  [-]
      corr_cholesky_mtx_2    ! Correlations Cholesky matrix (2nd comp.)  [-]

    real( kind = core_rknd ), dimension(ngrdcol,nzt,pdf_dim), intent(in) :: &
      mu1,    & ! Means of the hydrometeors, 1st comp. (chi, eta, w, <hydrometeors>)  [units vary]
      mu2,    & ! Means of the hydrometeors, 2nd comp. (chi, eta, w, <hydrometeors>)  [units vary]
      sigma1, & ! Stdevs of the hydrometeors, 1st comp. (chi, eta, w, <hydrometeors>) [units vary]
      sigma2    ! Stdevs of the hydrometeors, 2nd comp. (chi, eta, w, <hydrometeors>) [units vary]

    type(precipitation_fractions), intent(in) :: &
      precip_fracs           ! Precipitation fractions      [-]

    type(silhs_config_flags_type), intent(in) :: &
      silhs_config_flags ! Flags for the SILHS sampling code [-]

    real( kind = core_rknd ), intent(in) :: &
      vert_decorr_coef    ! Empirically defined de-correlation constant [-]

    type (stats_metadata_type), intent(in) :: &
      stats_metadata

    ! ---------------- InOut Variables ----------------
    type(stats), dimension(ngrdcol), intent(inout) :: &
      stats_lh_zt, &
      stats_lh_sfc

    type(err_info_type), intent(inout) :: &
      err_info          ! err_info struct containing err_code and err_header

    ! ---------------- Output Variables ----------------
    real( kind = core_rknd ), intent(out), dimension(ngrdcol,num_samples,nzt,pdf_dim) :: &
      X_nl_all_levs ! Sample that is transformed ultimately to normal-lognormal

    integer, intent(out), dimension(ngrdcol,num_samples,nzt) :: &
      X_mixt_comp_all_levs ! Which mixture component we're in

    real( kind = core_rknd ), intent(out), dimension(ngrdcol,num_samples,nzt) :: &
      lh_sample_point_weights ! Weight of each sample point

    ! ---------------- Local variables ----------------
    real( kind = core_rknd ), dimension(ngrdcol,num_samples,nzt,(pdf_dim+d_uniform_extra)) :: &
      X_u_all_levs ! Sample drawn from uniform distribution

    integer, dimension(ngrdcol) :: &
      k_lh_start ! Height for preferentially sampling within cloud
      
    integer :: &
      k, sample, p, i, j  ! Loop iterators

    logical, dimension(ngrdcol,num_samples,nzt) :: &
      l_in_precip   ! Whether sample is in precipitation

    logical :: l_error, l_error_in_sub
    
    real( kind = core_rknd ), dimension(pdf_dim,ngrdcol,nzt,pdf_dim) :: &
      Sigma_Cholesky1, &  ! Correlations Cholesky matrix 1 [-]
      Sigma_Cholesky2     ! Correlations Cholesky matrix 2 [-]
      
    real( kind = core_rknd ), dimension(ngrdcol,num_samples,nzt) :: &
      cloud_frac ! Cloud fraction for grid level and sample

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      rcm_pdf    ! Liquid water mixing ratio          [kg/kg]
      
    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      X_vert_corr         ! Vertical correlations between height levels       [-]
      
    real(kind = core_rknd), dimension(ngrdcol,num_samples,nzt,pdf_dim+d_uniform_extra) :: &
      rand_pool ! Array of randomly generated numbers

    ! ---------------- Begin Code ----------------

    !$acc enter data create( X_u_all_levs, k_lh_start, l_in_precip, Sigma_Cholesky1, &
    !$acc                    Sigma_Cholesky2, cloud_frac, rcm_pdf, X_vert_corr, rand_pool )
        
#ifdef CLUBB_GPU
    if ( .not. silhs_config_flags%l_lh_straight_mc ) then
      error stop "CLUBB ERROR: Running SILHS on GPUs requires lh_straight_mc=true"
    end if
#endif

    ! Sanity checks for l_lh_importance_sampling
    if ( silhs_config_flags%l_lh_importance_sampling .and. sequence_length /= 1 ) then
      write(fstderr,*) "Cloud weighted sampling requires sequence length be equal to 1."
      error stop "Fatal error."
    end if

    if ( silhs_config_flags%l_Lscale_vert_avg ) then
      !do k = 1, nzt, 1
      !  kp1 = min( k+1, nzt )
      !  km1 = max( k-1, 1 )
      !  Lscale_vert_avg(k) = vertical_avg &
      !                       ( (kp1-km1+1), rho_ds_zt(km1:kp1), &
      !                         Lscale(km1:kp1), gr%dzt(km1:kp1) )
      !end do
      error stop "CLUBB ERROR: l_Lscale_vert_avg has been depricated"
    end if

    l_error = .false.

    ! Compute rcm_pdf for use within SILHS
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzt
      do i = 1, ngrdcol
        rcm_pdf(i,k) = compute_mean_binormal( pdf_params%rc_1(i,k), pdf_params%rc_2(i,k), &
                                              pdf_params%mixt_frac(i,k) )
      end do
    end do

    !--------------------------------------------------------------
    ! Latin hypercube sampling
    !--------------------------------------------------------------

    ! Compute k_lh_start, the starting vertical grid level for SILHS sampling
    call compute_k_lh_start( gr, nzt, ngrdcol, rcm_pdf, pdf_params, &
                             silhs_config_flags%l_rcm_in_cloud_k_lh_start, &
                             silhs_config_flags%l_random_k_lh_start, &
                             k_lh_start )
                                     
    ! Calculate possible Sigma_Cholesky values
    ! Row-wise multiply of the elements of a lower triangular matrix.
    !$acc parallel loop gang vector collapse(4) default(present)
    do p = 1, pdf_dim
      do j = 1, pdf_dim
        do k = 1, nzt
          do i = 1, ngrdcol
            ! Calculate possible Sigma_Cholesky values
            Sigma_Cholesky1(j,i,k,p) = corr_cholesky_mtx_1(i,k,p,j) * sigma1(i,k,p)
            Sigma_Cholesky2(j,i,k,p) = corr_cholesky_mtx_2(i,k,p,j) * sigma2(i,k,p)
          end do
        end do
      end do
    end do
    
    ! Compute the vertical correlation for arbitrary overlap, using
    !   density weighted 3pt averaged Lscale and the difference in height levels (delta_zm)
    if ( silhs_config_flags%l_max_overlap_in_cloud ) then

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzt
        do i = 1, ngrdcol
          
          if ( rcm_pdf(i,k) > rc_tol ) then
            X_vert_corr(i,k) = one
          else
            X_vert_corr(i,k) &
            = exp( -vert_decorr_coef * ( gr%grid_dir * delta_zm(i,k) / Lscale(i,k) ) )
          end if
          
        end do
      end do

    else

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzt
        do i = 1, ngrdcol
          X_vert_corr(i,k) &
          = exp( -vert_decorr_coef * ( gr%grid_dir * delta_zm(i,k) / Lscale(i,k) ) )
        end do
      end do

    end if
    
    ! Assertion check for the vertical correlation
    if ( clubb_at_least_debug_level( 1 ) ) then

      !$acc update host( X_vert_corr )

      if ( any( X_vert_corr > one ) .or. any( X_vert_corr < zero ) ) then
        write(fstderr,*) "The vertical correlation in latin_hypercube_driver"// &
          "is not in the correct range"
        do k = 1, nzt
          do i = 1, ngrdcol
            write(fstderr,*) "i = ", i, "k = ", k, "Vert. correlation = ", X_vert_corr(i,k)
          end do
        end do
        error stop "Fatal error in vertical_overlap_driver"
      end if ! Some correlation isn't between [0,1]

    end if ! clubb_at_least_debug_level 1

    ! Generate pool of random numbers
    call generate_random_pool( nzt, ngrdcol, pdf_dim, num_samples, d_uniform_extra, & ! Intent(in)
                               lh_seed, gr,                                         & ! Intent(in)
                               rand_pool )                                            ! Intent(out)

    ! Generate all uniform samples, based on the rand pool
    call generate_all_uniform_samples( &
           iter, pdf_dim, d_uniform_extra, num_samples, sequence_length, & ! Intent(in)
           nzt, ngrdcol, k_lh_start, X_vert_corr, rand_pool,              & ! Intent(in)
           hm_metadata%iiPDF_chi,                                     & ! Intent(in)
           pdf_params%cloud_frac_1,                                      & ! Intent(in)
           pdf_params%cloud_frac_2,                                      & ! Intent(in)
           pdf_params%mixt_frac, precip_fracs,                           & ! Intent(in)
           silhs_config_flags%cluster_allocation_strategy,               & ! Intent(in)
           silhs_config_flags%l_lh_importance_sampling,                  & ! Intent(in)
           silhs_config_flags%l_lh_straight_mc,                          & ! Intent(in)
           silhs_config_flags%l_lh_clustered_sampling,                   & ! Intent(in)
           silhs_config_flags%l_lh_limit_weights,                        & ! Intent(in)
           silhs_config_flags%l_lh_var_frac,                             & ! Intent(in)
           silhs_config_flags%l_lh_normalize_weights,                    & ! Intent(in)
           l_calc_weights_all_levs_itime,                                & ! Intent(in)
           X_u_all_levs, lh_sample_point_weights )                         ! Intent(out)     
           
    !$acc parallel loop gang vector collapse(3) default(present)
    do k = 1, nzt
      do sample = 1, num_samples 
        do i = 1, ngrdcol
            
          ! Determine mixture component for all levels
          if ( X_u_all_levs(i,sample,k,pdf_dim+1) < pdf_params%mixt_frac(i,k) ) then
            
            ! Set pdf component indicator to 1 for this sample and vertical level
            X_mixt_comp_all_levs(i,sample,k) = 1
            
            ! Copy 1st component values
            cloud_frac(i,sample,k) = pdf_params%cloud_frac_1(i,k)
            
            ! Determine precipitation
            if ( X_u_all_levs(i,sample,k,pdf_dim+2) < precip_fracs%precip_frac_1(i,k) ) then
              l_in_precip(i,sample,k) = .true.
            else
              l_in_precip(i,sample,k) = .false.
            end if
            
          else
            
            ! Set pdf component indicator to 2 for this sample and vertical level
            X_mixt_comp_all_levs(i,sample,k) = 2
            
            ! Copy 2nd component values
            cloud_frac(i,sample,k) = pdf_params%cloud_frac_2(i,k)
            
            ! Determine precipitation
            if ( X_u_all_levs(i,sample,k,pdf_dim+2) < precip_fracs%precip_frac_2(i,k) ) then
              l_in_precip(i,sample,k) = .true.
            else
              l_in_precip(i,sample,k) = .false.
            end if
            
          end if

        end do
      end do 
    end do

    ! Generate LH sample, represented by X_u and X_nl, for level k
    ! Transform the uniformly distributed samples to
    !   ones distributed according to CLUBB's PDF.
    call transform_uniform_samples_to_pdf( &
           nzt, ngrdcol, num_samples, pdf_dim, d_uniform_extra, & ! In
           hm_metadata,                                     & ! In
           Sigma_Cholesky1, Sigma_Cholesky2,                   & ! In
           mu1, mu2, X_mixt_comp_all_levs,                     & ! In
           X_u_all_levs, cloud_frac,                           & ! In
           l_in_precip,                                        & ! In
           X_nl_all_levs )                                       ! Out
     
    
    if ( stats_metadata%l_stats_samp ) then
      !$acc update host(X_u_all_levs,l_in_precip,lh_sample_point_weights,X_mixt_comp_all_levs)
      call stats_accumulate_uniform_lh( nzt, num_samples, ngrdcol, l_in_precip(:,:,:), &
                                        X_mixt_comp_all_levs(:,:,:), &
                                        X_u_all_levs(:,:,:,hm_metadata%iiPDF_chi), pdf_params, &
                                        lh_sample_point_weights(:,:,:), k_lh_start(:), &
                                        stats_metadata, &
                                        stats_lh_zt, stats_lh_sfc )
    end if

    if ( l_output_2D_lognormal_dist ) then
      !$acc update host(X_nl_all_levs)
      
      ! Eric Raut removed lh_rt and lh_thl from call to output_2D_lognormal_dist_file
      ! because they are no longer generated in generate_silhs_sample.
      do i = 1, ngrdcol
        call output_2D_lognormal_dist_file( nzt, num_samples, pdf_dim, &
                                            real(X_nl_all_levs(i,:,:,:), kind = stat_rknd), &
                                            stats_metadata, err_info )
      end do

      if ( any(err_info%err_code == clubb_fatal_error) ) then
        write(fstderr, *) err_info%err_header_global
        write(fstderr, *) "Fatal error writing to the 2D LOGNORMAL sample netcdf file", &
                          " in CLUBB SILHS procedure generate_silhs_sample"
        return
      end if

    end if

    if ( l_output_2D_uniform_dist ) then
      !$acc update host(X_u_all_levs,X_mixt_comp_all_levs,lh_sample_point_weights)
      do i = 1, ngrdcol
        call output_2D_uniform_dist_file( nzt, num_samples, pdf_dim+2, &
                                          X_u_all_levs(i,:,:,:), &
                                          X_mixt_comp_all_levs(i,:,:), &
                                          lh_sample_point_weights(i,:,:), &
                                          stats_metadata, err_info )
      end do

      if ( any(err_info%err_code == clubb_fatal_error) ) then
        write(fstderr, *) err_info%err_header_global
        write(fstderr, *) "Fatal error writing to the 2D UNIFORM sample netcdf file", &
                          " in CLUBB SILHS procedure generate_silhs_sample"
        return
      end if

    end if

    ! Various nefarious assertion checks
    if ( clubb_at_least_debug_level( 2 ) ) then
      !$acc update host(X_u_all_levs,X_mixt_comp_all_levs,X_nl_all_levs)

      ! Simple assertion check to ensure uniform variates are in the appropriate
      ! range
      if ( any( X_u_all_levs <= zero .or. X_u_all_levs >= one ) ) then
        write(fstderr,*) "A uniform variate was not in the correct range."
        l_error = .true.
      end if

      do i = 1, ngrdcol
        do k = 1, nzt

          call assert_consistent_cloud_frac( pdf_params%chi_1(i,k), pdf_params%chi_2(i,k), & 
                                    pdf_params%cloud_frac_1(i,k), pdf_params%cloud_frac_2(i,k), &
                                    pdf_params%stdev_chi_1(i,k), pdf_params%stdev_chi_2(i,k), & 
                                    l_error_in_sub )
          l_error = l_error .or. l_error_in_sub

          ! Check for correct transformation in normal space
          call assert_correct_cloud_normal( num_samples, & ! In
                                            X_u_all_levs(i,:,k,hm_metadata%iiPDF_chi), & ! In
                                            X_nl_all_levs(i,:,k,hm_metadata%iiPDF_chi), & ! In
                                            X_mixt_comp_all_levs(i,:,k), & ! In
                                            pdf_params%cloud_frac_1(i,k), & ! In
                                            pdf_params%cloud_frac_2(i,k), & ! In
                                            l_error_in_sub ) ! Out
          l_error = l_error .or. l_error_in_sub

        end do ! k=1, nzt
      end do ! i=1, ngrdcol

    end if ! clubb_at_least_debug_level( 2 )

    ! Stop the run if an error occurred
    if ( l_error ) then
      error stop "Fatal error in generate_silhs_sample"
    end if

    !$acc exit data delete( X_u_all_levs, k_lh_start, l_in_precip, Sigma_Cholesky1, &
    !$acc                    Sigma_Cholesky2, cloud_frac, rcm_pdf, X_vert_corr, rand_pool )

    return

  end subroutine generate_silhs_sample
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
  subroutine generate_random_pool( nzt, ngrdcol, pdf_dim, num_samples, d_uniform_extra, &
                                   lh_seed, gr, &
                                   rand_pool )
  ! Description:
  !     This subroutine populates rand_pool with random numbers. There are
  !     different requirements for generating random numbers on a CPU vs a
  !     GPU, so the operation of this procedure depends on where or not 
  !     we compile with -DCLUBB_GPU has been specified at compile time.
  !     We use a CUDA routine to generate on a GPU when available, but if
  !     -DCUDA is not specified then we generate randoms on the GPU and
  !     copy the results to the GPU.
  !
  ! References:
  !   clubb ticket 869
  !
  ! Author: Gunther Huebler
  !----------------------------------------------------------------------------

    use grid_class, only: &
        grid    ! Type(s)

    use clubb_precision, only: &
      core_rknd      ! Constant     
    
    use generate_uniform_sample_module, only: &
      rand_uniform_real ! Procedure
      
    use mt95, only: &
      genrand_intg, &
      genrand_init_api  ! Procedure
      
#if defined(CUDA) && defined(CLUBB_GPU)
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
      nzt,               & ! Number of vertical levels
      ngrdcol,               & ! Number of grid columns
      pdf_dim,          & ! Variates
      num_samples,      & ! Number of samples
      d_uniform_extra     ! Uniform variates included in uniform sample but not
                          !  in normal/lognormal sample
                          
    integer( kind = genrand_intg ), intent(in) :: &
      lh_seed

    type(grid), intent(in) :: &
      gr    ! Grid variable type
      
    ! ------------------- Output Variables -------------------

    real(kind=core_rknd), dimension(ngrdcol,num_samples,nzt,pdf_dim+d_uniform_extra), intent(out)::&
      rand_pool   ! Pool of random reals, generated from integers
      
    ! ------------------- Local Variables -------------------
    
#if defined(CUDA) && defined(CLUBB_GPU)
      
    integer :: &
      r_status  ! Integer use to call curand functions

    type(curandGenerator) :: &
      cu_gen  ! curand generator variable

    logical, save :: &
      l_first_iter = .true.  ! First iteration indicator
      
#endif
    
    integer :: k, p, i, sample ! Loop variables
      
    ! ---------------- Begin Code ----------------
    
#if defined(CUDA) && defined(CLUBB_GPU)

    ! Generate randoms on GPU using a CUDA function
    
    ! If first iteration, intialize generator
    if ( l_first_iter ) then
      l_first_iter = .false.
      r_status = curandCreateGenerator( cu_gen, CURAND_RNG_PSEUDO_DEFAULT )
      r_status = curandSetPseudoRandomGeneratorSeed( cu_gen, lh_seed )
    end if
    
    !$acc host_data use_device(rand_pool) 
    r_status = curandGenerate(cu_gen, rand_pool, ngrdcol*num_samples*nzt*(pdf_dim+d_uniform_extra))
    !$acc end host_data
    
#else

    ! Generate randoms on CPU using an Mersenne Twister
    
    ! Intialize generator, this is required every timestep to enable restart runs to 
    ! produce bit-for-bit results
    call genrand_init_api( lh_seed )

    ! Populate rand_pool with a generator designed for a CPU
    do p=1, pdf_dim+d_uniform_extra
      do k = gr%k_lb_zt, gr%k_ub_zt, gr%grid_dir_indx
        do sample=1, num_samples
          do i = 1, ngrdcol
            rand_pool(i,sample,k,p) = rand_uniform_real()
          end do
        end do
      end do
    end do

#ifdef CLUBB_GPU
    ! If we are using the GPU, but generating random numbers on the CPU, then
    ! we need to copy the CPU generated randoms to the GPU
    !$acc update device( rand_pool )
#endif

#endif

#ifdef SILHS_MULTI_COL_RAND_DUPLICATE
    !$acc parallel loop gang vector collapse(4) default(present)
    do i = 2, ngrdcol
      do p=1, pdf_dim+d_uniform_extra
        do k = 1, nzt
          do sample=1, num_samples
              rand_pool(i,sample,k,p) = rand_pool(1,sample,k,p)
            end do
          end do
        end do
      end do
#endif
    
  end subroutine generate_random_pool

!-------------------------------------------------------------------------------
  subroutine generate_all_uniform_samples( &
               iter, pdf_dim, d_uniform_extra, num_samples, sequence_length, & ! Intent(in)
               nzt, ngrdcol, k_lh_start, X_vert_corr, rand_pool,              & ! Intent(in)
               iiPDF_chi,                                                    & ! Intent(in)
               cloud_frac_1,                                                 & ! Intent(in)
               cloud_frac_2,                                                 & ! Intent(in)
               mixt_frac, precip_fracs,                                      & ! Intent(in)
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
      precipitation_fractions      ! Type

    use generate_uniform_sample_module, only: &
      rand_uniform_real, &        ! Procedure(s)
      generate_uniform_lh_sample

    use silhs_importance_sample_module, only: &
      importance_sampling_driver, & ! Procedure(s)
      cloud_weighted_sampling_driver

    use latin_hypercube_arrays, only: &
      one_height_time_matrix      ! Variable

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
      nzt, ngrdcol
      
    integer, dimension(ngrdcol), intent(in) :: &
      k_lh_start

    integer, intent(in) :: &
      iiPDF_chi

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      cloud_frac_1, cloud_frac_2, &     ! The PDF parameters at k_lh_start
      mixt_frac, &
      X_vert_corr

    type(precipitation_fractions), intent(in) :: &
      precip_fracs           ! Precipitation fractions      [-]

    real(kind = core_rknd),dimension(ngrdcol,num_samples,nzt,pdf_dim+d_uniform_extra),intent(in)::&
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
    real(kind = core_rknd),dimension(ngrdcol,num_samples,nzt,pdf_dim+d_uniform_extra),intent(out)::&
      X_u_all_levs              ! Uniform sample at k_lh_start

    real( kind = core_rknd ), dimension(ngrdcol,num_samples,nzt), intent(out) :: &
      lh_sample_point_weights     ! Weight of each sample point (all equal to one if importance
                                  ! sampling is not used)

    ! ------------------ Local Variables ------------------
    integer :: &
      k, p, i, sample

    !------------------ Begin Code ------------------

    ! Sanity check
    if ( l_lh_old_cloud_weighted .and. mod( num_samples, 2 ) /= 0 ) then
      write(fstderr,*) "Old cloud weighted sampling requires num_samples to be divisible by 2."
      error stop "Fatal error."
    end if
    
    if ( .not. l_calc_weights_all_levs_itime ) then

      ! Generate random samples at k_lh_start, then use compute_arb_overlap to 
      ! populate the rest of the uniform samples

      if ( l_lh_straight_mc ) then

        ! Do a straight Monte Carlo sample without LH or importance sampling.
        
        !$acc parallel loop gang vector collapse(4) default(present) 
        do p = 1, pdf_dim+d_uniform_extra
          do k = 1, nzt
            do sample = 1, num_samples
              do i = 1, ngrdcol
                X_u_all_levs(i,sample,k,p) = max( single_prec_thresh, &
                                                min( one - single_prec_thresh, &
                                                     rand_pool(i,sample,k,p) ) )
              end do
            end do
          end do
        end do

        !$acc parallel loop gang vector collapse(3) default(present) 
        do k = 1, nzt
          do sample = 1, num_samples
            do i = 1, ngrdcol
              ! Importance sampling is not performed, so all sample points have the same weight!!
              lh_sample_point_weights(i,sample,k) = one
            end do
          end do
        end do

      else ! .not. l_lh_straight_mc
      
        do i = 1, ngrdcol
          ! Generate a uniformly distributed Latin hypercube sample
          call generate_uniform_lh_sample( iter, num_samples, sequence_length, & ! Intent(in)
                                           pdf_dim+d_uniform_extra,            & ! Intent(in)
                                           X_u_all_levs(i,:,k_lh_start(i),:) )   ! Intent(out)
                                
          if ( l_lh_importance_sampling ) then

            if ( l_lh_old_cloud_weighted ) then

              call cloud_weighted_sampling_driver &
                   ( num_samples, one_height_time_matrix(:,iiPDF_chi), & ! In
                     one_height_time_matrix(:,pdf_dim+1), & ! In
                     cloud_frac_1(i,k_lh_start(i)), cloud_frac_2(i,k_lh_start(i)), & ! In
                     mixt_frac(i,k_lh_start(i)), & ! In
                     X_u_all_levs(i,:,k_lh_start(i),iiPDF_chi), & ! In/Out
                     X_u_all_levs(i,:,k_lh_start(i),pdf_dim+1), & ! In/Out
                     lh_sample_point_weights(i,:,k_lh_start(i)) ) ! Out

            else ! .not. l_lh_old_cloud_weighted

              call importance_sampling_driver &
                   ( num_samples,                                                      & ! In
                     cloud_frac_1(i,k_lh_start(i)), cloud_frac_2(i,k_lh_start(i)),     & ! In
                     mixt_frac(i,k_lh_start(i)),                                       & ! In
                     precip_fracs%precip_frac_1(i,k_lh_start(i)),                      & ! In
                     precip_fracs%precip_frac_2(i,k_lh_start(i)),                      & ! In
                     cluster_allocation_strategy, l_lh_clustered_sampling,             & ! In
                     l_lh_limit_weights, l_lh_var_frac, l_lh_normalize_weights,        & ! In
                     X_u_all_levs(i,:,k_lh_start(i),iiPDF_chi),                        & ! In/Out
                     X_u_all_levs(i,:,k_lh_start(i),pdf_dim+1),                        & ! In/Out
                     X_u_all_levs(i,:,k_lh_start(i),pdf_dim+2),                        & ! In/Out
                     lh_sample_point_weights(i,:,k_lh_start(i)) )                        ! Out

            end if ! l_lh_old_cloud_weighted
            
            ! Clip uniform sample points to expected range                                 
            X_u_all_levs(i,:,k_lh_start(i),:) = max( single_prec_thresh, &
                            min( one - single_prec_thresh, X_u_all_levs(i,:,k_lh_start(i),:) ) )

            
            do k = 1, nzt
              lh_sample_point_weights(i,:,k) = lh_sample_point_weights(i,:,k_lh_start(i))
            end do

          else

            ! No importance sampling is performed, so all sample points have the same weight.
            lh_sample_point_weights(i,:,:) = one

          end if ! l_lh_importance_sampling
          
        end do


      end if ! l_lh_straight_mc
      
      ! Generate uniform sample at other grid levels by vertically correlating them
      ! https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:vert_corr
      call compute_arb_overlap( nzt, ngrdcol, num_samples, pdf_dim, d_uniform_extra, & ! In
                                k_lh_start, X_vert_corr, rand_pool,        & ! In 
                                X_u_all_levs )
                                
    else 
      
      ! Generate random samples for all vertical levels, samples, and variates
      
      if ( l_lh_straight_mc ) then

        ! Do a straight Monte Carlo sample without LH or importance sampling.
        
        !$acc parallel loop gang vector collapse(4) default(present) 
        do p = 1, pdf_dim+d_uniform_extra
          do k = 1, nzt
            do sample = 1, num_samples
              do i = 1, ngrdcol
                X_u_all_levs(i,sample,k,p) = max( single_prec_thresh, &
                                                min( one - single_prec_thresh, &
                                                     rand_pool(i,sample,k,p) ) )
              end do
            end do
          end do
        end do

        !$acc parallel loop gang vector collapse(3) default(present) 
        do k = 1, nzt
          do sample = 1, num_samples
            do i = 1, ngrdcol
              ! Importance sampling is not performed, so all sample points have the same weight!!
              lh_sample_point_weights(i,sample,k) = one
            end do
          end do
        end do

      else ! .not. l_lh_straight_mc

        do i = 1, ngrdcol
          do k = 1, nzt

            ! Generate a uniformly distributed Latin hypercube sample
            call generate_uniform_lh_sample( iter, num_samples, sequence_length, & ! Intent(in)
                                             pdf_dim+d_uniform_extra,            & ! Intent(in)
                                             X_u_all_levs(i,:,k,:) )               ! Intent(out)

            if ( l_lh_importance_sampling ) then

              if ( l_lh_old_cloud_weighted ) then

                call cloud_weighted_sampling_driver &
                     ( num_samples, one_height_time_matrix(:,iiPDF_chi), & ! In
                       one_height_time_matrix(:,pdf_dim+1),              & ! In
                       cloud_frac_1(i,k), cloud_frac_2(i,k),             & ! In
                       mixt_frac(i,k),                                   & ! In
                       X_u_all_levs(i,:,k,iiPDF_chi),                    & ! In/Out
                       X_u_all_levs(i,:,k,pdf_dim+1),                    & ! In/Out
                       lh_sample_point_weights(i,:,k) )                    ! Out

              else ! .not. l_lh_old_cloud_weighted

                call importance_sampling_driver &
                     ( num_samples,                                               & ! In
                       cloud_frac_1(i,k), cloud_frac_2(i,k),                      & ! In
                       mixt_frac(i,k),                                            & ! In
                       precip_fracs%precip_frac_1(i,k),                           & ! In
                       precip_fracs%precip_frac_2(i,k),                           & ! In
                       cluster_allocation_strategy, l_lh_clustered_sampling,      & ! In
                       l_lh_limit_weights, l_lh_var_frac, l_lh_normalize_weights, & ! In
                       X_u_all_levs(i,:,k,iiPDF_chi),                             & ! In/Out
                       X_u_all_levs(i,:,k,pdf_dim+1),                             & ! In/Out
                       X_u_all_levs(i,:,k,pdf_dim+2),                             & ! In/Out
                       lh_sample_point_weights(i,:,k) )                             ! Out

              end if ! l_lh_old_cloud_weighted

            else

              ! No importance sampling is performed, so all sample points have the same weight.
              lh_sample_point_weights(i,:,k) = one

            end if ! l_lh_importance_sampling
            
            ! Clip uniform sample points to expected range                                 
            X_u_all_levs(i,:,k,:) = max( single_prec_thresh, &
                                    min( one - single_prec_thresh, X_u_all_levs(i,:,k,:) ) )
          
          end do
        end do

      end if ! l_lh_straight_mc
      
    end if ! .not. l_calc_weights_all_levs_itime

    return
  end subroutine generate_all_uniform_samples
!-------------------------------------------------------------------------------

!-----------------------------------------------------------------------
  subroutine compute_k_lh_start( gr, nzt, ngrdcol, rcm_pdf, pdf_params, &
                                 l_rcm_in_cloud_k_lh_start, &
                                 l_random_k_lh_start, &
                                 k_lh_start )

  ! Description:
  !   Determines the starting SILHS sample level

  ! References:
  !   None
  !-----------------------------------------------------------------------

    ! Included Modules
    use grid_class, only: &
        grid    ! Type(s)

    use clubb_precision, only: &
      core_rknd      ! Constant

    use constants_clubb, only: &
      cloud_frac_min, & ! Constants
      zero

    use pdf_parameter_module, only: &
      pdf_parameter  ! Type

    use pdf_utilities, only: &
      compute_mean_binormal  ! Procedure

    use math_utilities, only: &
      rand_integer_in_range  ! Procedure

    implicit none

    ! Input Variables
    type(grid), intent(in) :: &
      gr    ! Grid variable type

    integer, intent(in) :: &
      nzt, &   ! Number of vertical levels
      ngrdcol ! Number of grid columns

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      rcm_pdf      ! Liquid water mixing ratio               [kg/kg]

    type(pdf_parameter), intent(in) :: &
      pdf_params  ! PDF parameters       [units vary]

    logical, intent(in) :: &
      l_rcm_in_cloud_k_lh_start, & ! Determine k_lh_start based on maximum within-cloud rcm
      l_random_k_lh_start          ! k_lh_start found randomly between max rcm and rcm_in_cloud

    ! Output Variable
    integer, dimension(ngrdcol), intent(out) :: &
      k_lh_start  ! Starting SILHS sample level

    ! Local Variables
    integer, dimension(ngrdcol) :: &
      k_lh_start_rcm_in_cloud, &
      k_lh_start_rcm

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      cloud_frac_pdf

    real( kind = core_rknd ) :: &
      rcm_in_cloud_max, &
      rcm_in_cloud, &
      rcm_max
      
    integer :: i, k  ! Loop iterator

  !------------------------------ Begin Code -----------------------------------

    !$acc enter data create( k_lh_start_rcm_in_cloud, k_lh_start_rcm, cloud_frac_pdf )

    if ( l_rcm_in_cloud_k_lh_start .or. l_random_k_lh_start ) then

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, nzt
        do i = 1, ngrdcol
          cloud_frac_pdf(i,k) = compute_mean_binormal( pdf_params%cloud_frac_1(i,k), &
                                                       pdf_params%cloud_frac_2(i,k), &
                                                       pdf_params%mixt_frac(i,k) )
        end do
      end do

      !$acc parallel loop gang vector default(present)
      do i = 1, ngrdcol

        ! Initialize k_lh_start_rcm_in_cloud to nzt/2
        k_lh_start_rcm_in_cloud(i) = nzt / 2
        ! Use the same level with a descending grid as would be used
        ! with an ascending grid.
        if ( gr%grid_dir_indx < 0 ) then
           k_lh_start_rcm_in_cloud(i) = nzt - k_lh_start_rcm_in_cloud(i) + 1
        endif
        rcm_in_cloud_max = zero

        ! Loop over vertical levels, if any rcm_pdf > 0 then set k_lh_start_rcm_in_cloud
        ! to the vertical level that results in the maximum value for rcm_in_cloud
        do k = 1, nzt

          rcm_in_cloud = rcm_pdf(i,k) / max( cloud_frac_pdf(i,k), cloud_frac_min  )

          if ( rcm_in_cloud > rcm_in_cloud_max ) then
            rcm_in_cloud_max = rcm_in_cloud
            k_lh_start_rcm_in_cloud(i) = k
          end if

        end do

      end do

    end if

    if ( .not. l_rcm_in_cloud_k_lh_start .or. l_random_k_lh_start ) then

      !$acc parallel loop gang vector default(present)
      do i = 1, ngrdcol

        ! Initialize k_lh_start_rcm_in_cloud to nzt/2
        k_lh_start_rcm(i) = nzt / 2
        ! Use the same level with a descending grid as would be used
        ! with an ascending grid.
        if ( gr%grid_dir_indx < 0 ) then
           k_lh_start_rcm(i) = nzt - k_lh_start_rcm(i) + 1
        endif
        rcm_max = zero

        do k = 1, nzt
          if ( rcm_pdf(i,k) > rcm_max ) then
            rcm_max = rcm_pdf(i,k)
            k_lh_start_rcm(i) = k
          endif
        end do

      end do

    end if

    if ( l_random_k_lh_start ) then

      !$acc update host( k_lh_start_rcm_in_cloud, k_lh_start_rcm )

      do i = 1, ngrdcol

        ! Pick a random height level between k_lh_start_rcm and k_lh_start_rcm_in_cloud
        if ( k_lh_start_rcm_in_cloud(i) > k_lh_start_rcm(i) ) then
          k_lh_start(i) = rand_integer_in_range( k_lh_start_rcm(i), k_lh_start_rcm_in_cloud(i) )
        else if ( k_lh_start_rcm(i) > k_lh_start_rcm_in_cloud(i) ) then
          k_lh_start(i) = rand_integer_in_range( k_lh_start_rcm_in_cloud(i), k_lh_start_rcm(i) )
        else
          ! k_lh_start_rcm_in_cloud = k_lh_start_rcm(i)
          k_lh_start(i) = k_lh_start_rcm(i)
        end if

        ! Keep the randomized index when the grid direction is descending
        ! to be the same level when the grid direction is ascending for
        ! purposes of comparison.
        if ( gr%grid_dir_indx < 0 ) then
          k_lh_start(i) &
          = k_lh_start_rcm_in_cloud(i) + k_lh_start_rcm(i) - k_lh_start(i)
        endif
          
      end do

      !$acc update device( k_lh_start )

    else if ( l_rcm_in_cloud_k_lh_start ) then

      !$acc parallel loop gang vector default(present)
      do i = 1, ngrdcol
        k_lh_start(i) = k_lh_start_rcm_in_cloud(i)
      end do

    else ! .not. l_random_k_lh_start .and. .not. l_rcm_in_cloud_k_lh_start

      !$acc parallel loop gang vector default(present)
      do i = 1, ngrdcol
        k_lh_start(i) = k_lh_start_rcm(i)
      end do
      
    end if

    !$acc exit data delete( k_lh_start_rcm_in_cloud, k_lh_start_rcm, cloud_frac_pdf )

    return

  end subroutine compute_k_lh_start
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  subroutine clip_transform_silhs_output( nzt, ngrdcol, num_samples,           & ! In
                                          pdf_dim, hydromet_dim, hm_metadata,  & ! In
                                          X_mixt_comp_all_levs,                   & ! In
                                          X_nl_all_levs,                          & ! Inout
                                          pdf_params, l_use_Ncn_to_Nc,            & ! In
                                          lh_rt_clipped, lh_thl_clipped,          & ! Out
                                          lh_rc_clipped, lh_rv_clipped,           & ! Out
                                          lh_Nc_clipped                           ) ! Out

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

    use corr_varnce_module, only: &
        hm_metadata_type

    use transform_to_pdf_module, only: &
        chi_eta_2_rtthl ! Awesome procedure

    use grid_class, only: grid

    implicit none

    ! ------------------- Input Variables -------------------
    logical, intent(in) :: &
      l_use_Ncn_to_Nc  ! Whether to call Ncn_to_Nc (.true.) or not (.false.);
                       ! Ncn_to_Nc might cause problems with the MG microphysics
                       ! since the changes made here (Nc-tendency) are not fed
                       ! into the microphysics

    integer, intent(in) :: &
      nzt,           & ! Number of vertical levels
      ngrdcol,           & ! Number of grid columns
      num_samples,  & ! Number of SILHS sample points
      pdf_dim,      & ! Number of variates in X_nl
      hydromet_dim    ! Number of hydrometeor species

    type (hm_metadata_type), intent(in) :: &
      hm_metadata

    integer, dimension(ngrdcol,num_samples,nzt), intent(in) :: &
      X_mixt_comp_all_levs   ! Which component this sample is in (1 or 2)

    type(pdf_parameter), intent(in) :: &
      pdf_params             ! **The** PDF parameters!

    ! ------------------- Input/Output Variable -------------------
    
    real( kind = core_rknd ), dimension(ngrdcol,num_samples,nzt,pdf_dim), intent(inout) :: &
      X_nl_all_levs    ! SILHS sample points    [units vary]

    ! ------------------- Output Variables -------------------
    real( kind = core_rknd ), dimension(ngrdcol,num_samples,nzt), intent(out) :: &
      lh_rt_clipped,  & ! rt generated from silhs sample points
      lh_thl_clipped, & ! thl generated from silhs sample points
      lh_rc_clipped,  & ! rc generated from silhs sample points
      lh_rv_clipped,  & ! rv generated from silhs sample points
      lh_Nc_clipped     ! Nc generated from silhs sample points

    ! ------------------- Local Variables -------------------

    real( kind = core_rknd ), dimension(nzt,hydromet_dim) :: &
      hydromet_pts,         & ! Sample point column of hydrometeors    [un vary]
      hydromet_pts_clipped    ! Clipped sample point column of hydromet   [un v]

    integer :: &
      sample, k, hm_idx, i

    ! Flag to clip sample points of hydrometeor concentrations.
    logical, parameter :: &
      l_clip_hydromet_samples = .false.
      
  !-----------------------------------------------------------------------

    ! Calculate (and clip) the SILHS sample point values of rt, thl, rc, rv,
    ! and Nc.
        
    ! Compute lh_rt and lh_thl
    call chi_eta_2_rtthl( nzt, ngrdcol, num_samples,                     & ! Intent(in)
                          pdf_params%rt_1, pdf_params%thl_1,            & ! Intent(in)
                          pdf_params%rt_2, pdf_params%thl_2,            & ! Intent(in)
                          pdf_params%crt_1, pdf_params%cthl_1,          & ! Intent(in)
                          pdf_params%crt_2, pdf_params%cthl_2,          & ! Intent(in)
                          pdf_params%chi_1, pdf_params%chi_2,           & ! Intent(in)
                          X_nl_all_levs(:,:,:,hm_metadata%iiPDF_chi),& ! Intent(in) 
                          X_nl_all_levs(:,:,:,hm_metadata%iiPDF_eta),& ! Intent(in)
                          X_mixt_comp_all_levs(:,:,:),                  & ! Intent(in)
                          lh_rt_clipped(:,:,:), lh_thl_clipped(:,:,:)   ) ! Intent(out)
    
    !$acc parallel loop gang vector collapse(3) default(present) 
    do k = 1, nzt
      do sample = 1, num_samples
        do i = 1, ngrdcol
    
          ! If necessary, clip rt      
          lh_rt_clipped(i,sample,k) = max( lh_rt_clipped(i,sample,k), rt_tol )
          
          ! Compute lh_rc, rc = chi * H(chi), where H(x) is the Heaviside step function
          lh_rc_clipped(i,sample,k) = max( X_nl_all_levs(i,sample,k,hm_metadata%iiPDF_chi), zero )
          
          ! Clip lh_rc.
          lh_rc_clipped(i,sample,k) = min( lh_rc_clipped(i,sample,k), &
                                              lh_rt_clipped(i,sample,k) - rt_tol )
          
          ! Compute lh_rv
          lh_rv_clipped(i,sample,k) = lh_rt_clipped(i,sample,k) - lh_rc_clipped(i,sample,k)
          
          if ( l_use_Ncn_to_Nc ) then
             ! Compute lh_Nc, Nc = Ncn * H(chi), where H(x) is the Heaviside step function
             if ( X_nl_all_levs(i,sample,k,hm_metadata%iiPDF_chi) > zero ) then
               lh_Nc_clipped(i,sample,k) = X_nl_all_levs(i,sample,k,hm_metadata%iiPDF_Ncn)
             else
               lh_Nc_clipped(i,sample,k) = zero
             end if
          else
             lh_Nc_clipped(i,sample,k) = X_nl_all_levs(i,sample,k,hm_metadata%iiPDF_Ncn)
          endif ! l_use_Ncn_to_Nc
        
        end do
      end do ! sample = 1, num_samples
    end do ! k = 1, nzt


    ! Clip the SILHS sample point values of hydrometeor concentrations.
    if ( l_clip_hydromet_samples ) then
      
#ifdef CLUBB_GPU
       error stop "CLUBB ERROR: Running SILHS on the GPUs requires l_clip_hydromet_samples=false"
#endif
       ! Loop over all sample columns.
       do sample = 1, num_samples, 1
         do i = 1, ngrdcol

           ! Pack the SILHS hydrometeor sample points, which are stored in arrays
           ! with the size pdf_dim, into arrays with the size hydromet_dim.
           do hm_idx = 1, hydromet_dim, 1
              hydromet_pts(:,hm_idx) &
              = X_nl_all_levs(i,sample,:,hydromet2pdf_idx(hm_idx,hm_metadata))
           enddo ! hm_idx = 1, hydromet_dim, 1

           ! Clip the hydrometeor concentration sample points based on
           ! maintaining the value of the hydrometeor mixing ratio sample points
           ! and satisfying the maximum allowable mean volume radius for that
           ! hydrometeor species.
           call clip_hydromet_conc_mvr( nzt, hydromet_dim, hm_metadata, & ! In
                                        hydromet_pts,                     & ! In
                                        hydromet_pts_clipped )              ! Out

           ! Unpack the clipped SILHS hydrometeor sample points, which are stored
           ! in arrays with the size hydromet_dim, back into arrays with the size
           ! pdf_dim.
           do hm_idx = 1, hydromet_dim, 1
             X_nl_all_levs(i,sample,:,hydromet2pdf_idx(hm_idx,hm_metadata)) &
              = hydromet_pts_clipped(:,hm_idx)
           enddo ! hm_idx = 1, hydromet_dim, 1

         end do ! i = 1, ngrdcol
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
           l_error_in_sub )                    ! Intent(out)

    l_error = l_error .or. l_error_in_sub
    if ( l_error_in_sub ) then
      write(fstderr,*) "Cloud fraction is inconsistent in PDF component 1"
    end if

    ! Perform assertion check for PDF component 2
    call assert_consistent_cf_component &
         ( chi_2, stdev_chi_2, cloud_frac_2, & ! Intent(in)
           l_error_in_sub )                    ! Intent(out)

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
  subroutine assert_correct_cloud_normal( num_samples, &
                                          X_u_chi, &
                                          X_nl_chi, &
                                          X_mixt_comp, &
                                          cloud_frac_1, &
                                          cloud_frac_2, &
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
        error stop "X_u_chi not in correct range in assert_correct_cloud_normal"
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
  subroutine latin_hypercube_2D_output_api( &
               fname_prefix, fdir, stats_tout, &
               nzt, pdf_dim, &
               stats_zt, time_initial, num_samples, &
               nlon, nlat, lon_vals, lat_vals, &
               hm_metadata, &
               clubb_params, sclr_dim, sclr_tol, &
               l_uv_nudge, &
               l_tke_aniso, &
               l_standard_term_ta, &
               err_info )
!-------------------------------------------------------------------------------

    use corr_varnce_module, only: &
      hm_metadata_type

    use clubb_precision, only: &
      time_precision, & ! Constant
      core_rknd

    use constants_clubb, only: &
      fstderr       ! Constant(s)

    use output_2D_samples_module, only: &
      open_2D_samples_file ! Procedure

    use output_2D_samples_module, only: &
      lognormal_sample_file, & ! Instance of a type
      uniform_sample_file

    use parameter_indices, only: &
      nparams    ! Variable(s)

    use err_info_type_module, only: &
      err_info_type         ! Type

    use error_code, only: &
        clubb_fatal_error   ! Constant

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
      nzt, & ! Number of vertical levels
      pdf_dim

    real( kind = core_rknd ), dimension(nzt), intent(in) :: &
      stats_zt ! Altitudes [m]

    integer, intent(in) :: num_samples

    integer, intent(in) :: &
      nlon, & ! Number of points in the X direction [-]
      nlat    ! Number of points in the Y direction [-]

    real( kind = core_rknd ), dimension(nlon), intent(in) ::  &
      lon_vals  ! Longitude values [Degrees E]

    real( kind = core_rknd ), dimension(nlat), intent(in) ::  &
      lat_vals  ! Latitude values  [Degrees N]

    type (hm_metadata_type), intent(in) :: &
      hm_metadata

    real( kind = core_rknd ), dimension(nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

    integer, intent(in) :: &
      sclr_dim 

    real( kind = core_rknd ), dimension(sclr_dim), intent(in) :: &
      sclr_tol

    logical, intent(in) :: &
      l_uv_nudge,         & ! For wind speed nudging
      l_tke_aniso,        & ! For anisotropic turbulent kinetic energy, i.e. TKE = 1/2
                            ! (u'^2 + v'^2 + w'^2)
      l_standard_term_ta    ! Use the standard discretization for the turbulent advection terms.
                            ! Setting to .false. means that a_1 and a_3 are pulled outside of the
                            ! derivative in advance_wp2_wp3_module.F90 and in
                            ! advance_xp2_xpyp_module.F90.

    ! Input/Output Variables
    type(err_info_type), intent(inout) :: &
      err_info              ! err_info struct containing err_code and err_header

    ! Local Variables
    character(len=100), allocatable, dimension(:) :: &
      variable_names, variable_descriptions, variable_units

    integer :: p

    integer :: &
      iiPDF_chi, &
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

    ! ---- Begin Code ----

    iiPDF_chi = hm_metadata%iiPDF_chi
    iiPDF_eta = hm_metadata%iiPDF_eta
    iiPDF_w   = hm_metadata%iiPDF_w
    iiPDF_rr  = hm_metadata%iiPDF_rr
    iiPDF_ri  = hm_metadata%iiPDF_ri
    iiPDF_rs  = hm_metadata%iiPDF_rs
    iiPDF_rg  = hm_metadata%iiPDF_rg
    iiPDF_Nr  = hm_metadata%iiPDF_Nr
    iiPDF_Ni  = hm_metadata%iiPDF_Ni
    iiPDF_Ns  = hm_metadata%iiPDF_Ns
    iiPDF_Ng  = hm_metadata%iiPDF_Ng
    iiPDF_Ncn = hm_metadata%iiPDF_Ncn

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

      call open_2D_samples_file( nzt, num_samples, pdf_dim, & ! In
                                 trim( fname_prefix )//"_nl", fdir, & ! In
                                 time_initial, stats_tout, stats_zt, variable_names, & ! In
                                 variable_descriptions, variable_units, & ! In
                                 nlon, nlat, lon_vals, lat_vals, & ! In
                                 clubb_params, sclr_dim, sclr_tol, &
                                 l_uv_nudge, &
                                 l_tke_aniso, &
                                 l_standard_term_ta, &
                                 lognormal_sample_file, err_info ) ! In/Out

      if ( any(err_info%err_code == clubb_fatal_error) ) then
        write(fstderr, *) err_info%err_header_global
        write(fstderr, *) "Fatal error calling open_2D_samples_file for LOGNORMAL samples", &
                          " in CLUBB SILHS procedure latin_hypercube_2D_output_api"
        return
      end if

      deallocate( variable_names, variable_descriptions, variable_units )

    end if ! l_output_2D_lognormal_dist

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

      p = pdf_dim+ 1
      variable_names(p) = "dp1"
      variable_descriptions(p) = "Uniform distribution for the mixture component"

      p = pdf_dim+ 2
      variable_names(p) = "dp2"
      variable_descriptions(p) = "Uniform variate used to determine precipitation!"

      p = pdf_dim+ 3
      variable_names(p) = "X_mixt_comp"
      variable_descriptions(p) = "Mixture component (should be 1 or 2)"

      p = pdf_dim+ 4
      variable_names(p) = "lh_sample_point_weights"
      variable_descriptions(p) = "Weight of each sample point"

      ! Set all the units
      variable_units(:) = "count" ! Unidata units format for a dimensionless quantity

      call open_2D_samples_file( nzt, num_samples, p, & ! In
                                 trim( fname_prefix )//"_u", fdir, & ! In
                                 time_initial, stats_tout, stats_zt, & ! In
                                 variable_names(1:p), variable_descriptions(1:p), & ! In
                                 variable_units(1:p), & ! In
                                 nlon, nlat, lon_vals, lat_vals, & ! In
                                 clubb_params, sclr_dim, sclr_tol, &
                                 l_uv_nudge, &
                                 l_tke_aniso, &
                                 l_standard_term_ta, &
                                 uniform_sample_file, err_info ) ! In/Out

      if ( any(err_info%err_code == clubb_fatal_error) ) then
        write(fstderr, *) err_info%err_header_global
        write(fstderr, *) "Fatal error calling open_2D_samples_file for UNIFORM samples", &
                          " in CLUBB SILHS procedure latin_hypercube_2D_output_api"
        return
      end if

      deallocate( variable_names, variable_descriptions, variable_units )

    end if ! l_output_2D_uniform_dist

    return
  end subroutine latin_hypercube_2D_output_api

!-------------------------------------------------------------------------------
  subroutine latin_hypercube_2D_close_api
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
  end subroutine latin_hypercube_2D_close_api

!-------------------------------------------------------------------------------
  subroutine compute_arb_overlap( nzt, ngrdcol, num_samples, pdf_dim, d_uniform_extra, &
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
      one, &
      two
      
    use parameters_silhs, only: &
      single_prec_thresh

    implicit none

    ! ---------------- Input Variables ----------------
    integer, intent(in) :: &
      nzt,               & ! Vertical levels
      ngrdcol,               & ! Columns
      pdf_dim,          & ! Number of variates in CLUBB's PDF
      d_uniform_extra,  & ! Uniform variates included in uniform sample but not
                                  !  in normal/lognormal sample
      num_samples         ! Number of SILHS sample points
      
    integer, dimension(ngrdcol), intent(in) :: &
      k_lh_start

    real(kind=core_rknd), dimension(ngrdcol,nzt), intent(in) :: &
      vert_corr ! Vertical correlation between k points in range [0,1]   [-]
      
    real(kind = core_rknd), dimension(ngrdcol,num_samples,nzt,pdf_dim+d_uniform_extra) :: &
      rand_pool ! Array of randomly generated numbers

    ! ---------------- Output Variables ----------------
    real(kind=core_rknd),dimension(ngrdcol,num_samples,nzt,pdf_dim+d_uniform_extra),intent(inout)::&
      X_u_all_levs ! Uniform distribution of 1 variate at all levels [-]
                   ! The value of this variate at k_lh_start should already be populated
                   ! in this array and will be used to fill in the other levels.

    ! ---------------- Local Variables ----------------
    real(kind=core_rknd) :: min_val, half_width, offset, unbounded_point

    integer :: k, sample, p, i, k_lh_start_i ! Loop iterators

    ! ---------------- Begin Code ----------------
    
    ! Recompute from k_lh_start to nzt-1 for all samples and variates, upward loop
    !$acc parallel loop gang vector collapse(3) default(present) 
    do p = 1, pdf_dim + d_uniform_extra                         
      do sample = 1, num_samples
        do i = 1, ngrdcol
        
          k_lh_start_i = k_lh_start(i)
          unbounded_point = X_u_all_levs(i,sample,k_lh_start_i,p)
        
          do k = k_lh_start_i, nzt-1
           
            half_width = one - vert_corr(i,k+1)
            min_val = unbounded_point - half_width

            offset = two * half_width * rand_pool(i,sample,k+1,p)

            unbounded_point = min_val + offset
           
            ! If unbounded_point lies outside [single_prec_thresh,1-single_prec_thresh],
            ! fold it back so that it is between the valid range
            if ( unbounded_point > one - single_prec_thresh ) then
              unbounded_point = two - unbounded_point - two * single_prec_thresh
            else if ( unbounded_point < single_prec_thresh ) then
              unbounded_point = - unbounded_point + two * single_prec_thresh
            end if
           
            X_u_all_levs(i,sample,k+1,p) = unbounded_point
           
          end do ! k_lh_start..nzt-1
        end do ! 1..ngrdcol
      end do ! 1..num_samples
    end do ! 1..pdf_dim


    ! Recompute from k_lh_start down to 2 for all samples and variates, downward loop 
    !$acc parallel loop gang vector collapse(3) default(present) 
    do p = 1, pdf_dim + d_uniform_extra                        
      do sample = 1, num_samples
        do i = 1, ngrdcol
        
          k_lh_start_i = k_lh_start(i)
          unbounded_point = X_u_all_levs(i,sample,k_lh_start_i,p)

          do k = k_lh_start_i, 2, -1

            half_width = one - vert_corr(i,k-1)
            min_val = unbounded_point - half_width

            offset = two * half_width * rand_pool(i,sample,k-1,p)

            unbounded_point = min_val + offset
           
            ! If unbounded_point lies outside [single_prec_thresh,1-single_prec_thresh],
            ! fold it back so that it is between the valid range 
            if ( unbounded_point > one - single_prec_thresh ) then
              unbounded_point = two - unbounded_point - two * single_prec_thresh
            else if ( unbounded_point < single_prec_thresh ) then
              unbounded_point = - unbounded_point + two * single_prec_thresh
            end if
           
            X_u_all_levs(i,sample,k-1,p) = unbounded_point

          end do ! k_lh_start..2 decrementing
        end do
      end do ! 1..num_samples
    end do ! 1..pdf_dim

    return

  end subroutine compute_arb_overlap

!-------------------------------------------------------------------------------
  subroutine stats_accumulate_lh_api( &
               gr, nzt, num_samples, pdf_dim, rho_ds_zt, &
               hydromet_dim, hm_metadata,&
               lh_sample_point_weights, X_nl_all_levs, &
               lh_rt_clipped, lh_thl_clipped, & 
               lh_rc_clipped, lh_rv_clipped, & 
               lh_Nc_clipped, &
               stats_metadata, &
               stats_lh_zt, stats_lh_sfc )

! Description:
!   Clip subcolumns from latin hypercube and create stats for diagnostic
!   purposes.

! References:
!   None
!-------------------------------------------------------------------------------

    use grid_class, only: grid

    use math_utilities, only: & ! Procedure(s)
      compute_sample_mean, & ! Procedure(s)
      compute_sample_variance

    use stats_type_utilities, only: &
!      stat_update_var, & ! Procedure(s)
      stat_update_var_pt

    use corr_varnce_module, only: &
      hm_metadata_type

    use constants_clubb, only: & 
      zero, &            ! Constant(s)
      one

    use clubb_precision, only: & 
      core_rknd    ! Constant

    use advance_helper_module, only: &
      vertical_integral ! Procedure(s)

    use stats_type, only: &
      stats ! Type

    use stats_variables, only: &
      stats_metadata_type

    implicit none

    !-------------------------- Input Variables --------------------------
    type (grid), intent(in) :: gr

    integer, intent(in) :: &
      pdf_dim,        & ! Number of variables to sample
      num_samples,    & ! Number of calls to microphysics per timestep (normally=2)
      nzt,             & ! Number of vertical model levels
      hydromet_dim      ! Number of hydrometeor species

    type (hm_metadata_type), intent(in) :: &
      hm_metadata

    real( kind = core_rknd ), intent(in), dimension(nzt) :: &
      rho_ds_zt  ! Dry, static density (thermo. levs.) [kg/m^3]

    real( kind = core_rknd ), intent(in), dimension(num_samples,nzt) :: &
      lh_sample_point_weights

    real( kind = core_rknd ), intent(in), dimension(num_samples,nzt,pdf_dim) :: &
      X_nl_all_levs ! Sample that is transformed ultimately to normal-lognormal

    real( kind = core_rknd ), dimension(num_samples,nzt), intent(in) :: &
      lh_rt_clipped,  & ! rt generated from silhs sample points
      lh_thl_clipped, & ! thl generated from silhs sample points
      lh_rc_clipped,  & ! rc generated from silhs sample points
      lh_rv_clipped,  & ! rv generated from silhs sample points
      lh_Nc_clipped     ! Nc generated from silhs sample points

    type (stats_metadata_type), intent(in) :: &
      stats_metadata

    !-------------------------- InOut Variables --------------------------
    type(stats), intent(inout) :: &
      stats_lh_zt, &
      stats_lh_sfc

    !-------------------------- Local variables --------------------------
    real( kind = core_rknd ), dimension(num_samples,nzt) :: &
      Ncn_all_points ! Cloud nuclei conc. for all levs.; Nc=Ncn*H(chi) [#/kg]

    real( kind = core_rknd ), dimension(num_samples,nzt,hydromet_dim) :: &
      hydromet_all_points ! Hydrometeor species    [units vary]

    real( kind = core_rknd ), dimension(nzt,hydromet_dim) :: &
      lh_hydromet ! Average value of the latin hypercube est. of all hydrometeors [units vary]

    real( kind = core_rknd ), dimension(nzt) :: &
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

    integer :: sample, ivar, k

    integer :: &
      iirr, & 
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

    !-------------------------- Begin Code --------------------------

    iirr      = hm_metadata%iirr
    iirs      = hm_metadata%iirs
    iiri      = hm_metadata%iiri
    iirg      = hm_metadata%iirg
    iiNr      = hm_metadata%iiNr
    iiNs      = hm_metadata%iiNs
    iiNi      = hm_metadata%iiNi
    iiNg      = hm_metadata%iiNg
    iiPDF_chi = hm_metadata%iiPDF_chi
    iiPDF_eta = hm_metadata%iiPDF_eta
    iiPDF_w   = hm_metadata%iiPDF_w
    iiPDF_Ncn = hm_metadata%iiPDF_Ncn

    if ( stats_metadata%l_stats_samp ) then

      ! For all cases where l_lh_importance_sampling is false, the weights
      ! will be 1 (all points equally weighted)

      if ( stats_metadata%ilh_rcm + stats_metadata%ilh_rcp2_zt + stats_metadata%ilh_lwp > 0 ) then
        lh_rcm = compute_sample_mean( nzt, num_samples, lh_sample_point_weights, &
                                      lh_rc_clipped )
        ! Switch back to using stat_update_var once the code is generalized
        ! to pass in the number of vertical levels.
!        call stat_update_var( stats_metadata%ilh_rcm, lh_rcm, stats_lh_zt )
        do k = 1, nzt
           call stat_update_var_pt( stats_metadata%ilh_rcm, k, lh_rcm(k), stats_lh_zt )
        enddo ! k = 1, nzt

        if ( stats_metadata%ilh_lwp > 0 ) then
          xtmp &
          = vertical_integral &
               ( nzt, rho_ds_zt(1:nzt), &
                 lh_rcm(1:nzt), gr%dzt(1,1:nzt) )

          call stat_update_var_pt( stats_metadata%ilh_lwp, 1, xtmp, stats_lh_sfc )
        end if
      end if

      if ( stats_metadata%ilh_sample_weights_sum > 0 ) then
          xtmp = sum(lh_sample_point_weights(:,:))
          call stat_update_var_pt( stats_metadata%ilh_sample_weights_sum, 1, xtmp, stats_lh_sfc )
      end if
      
      if ( stats_metadata%ilh_sample_weights_avg > 0 ) then
          xtmp = sum(lh_sample_point_weights(:,:)) / real( num_samples*nzt, kind = core_rknd )
          call stat_update_var_pt( stats_metadata%ilh_sample_weights_avg, 1, xtmp, stats_lh_sfc )
      end if
        
      if ( stats_metadata%ilh_thlm + stats_metadata%ilh_thlp2_zt > 0 ) then
        lh_thlm = compute_sample_mean( nzt, num_samples, lh_sample_point_weights, &
                                       real( lh_thl_clipped, kind = core_rknd ) )
        ! Switch back to using stat_update_var once the code is generalized
        ! to pass in the number of vertical levels.
!        call stat_update_var( stats_metadata%ilh_thlm, lh_thlm, stats_lh_zt )
        do k = 1, nzt
           call stat_update_var_pt( stats_metadata%ilh_thlm, k, lh_thlm(k), stats_lh_zt )
        enddo ! k = 1, nzt
      end if

      if ( stats_metadata%ilh_rvm + stats_metadata%ilh_rtp2_zt > 0 ) then
        lh_rvm = compute_sample_mean( nzt, num_samples, lh_sample_point_weights, &
                                      lh_rv_clipped )
        ! Switch back to using stat_update_var once the code is generalized
        ! to pass in the number of vertical levels.
!        call stat_update_var( stats_metadata%ilh_rvm, lh_rvm, stats_lh_zt )
        do k = 1, nzt
           call stat_update_var_pt( stats_metadata%ilh_rvm, k, lh_rvm(k), stats_lh_zt )
        enddo ! k = 1, nzt
        if ( stats_metadata%ilh_vwp > 0 ) then
          xtmp &
          = vertical_integral &
               ( nzt, rho_ds_zt(1:nzt), &
                 lh_rvm(1:nzt), gr%dzt(1,1:nzt) )

          call stat_update_var_pt( stats_metadata%ilh_vwp, 1, xtmp, stats_lh_sfc )
        end if
      end if

      if ( stats_metadata%ilh_wm + stats_metadata%ilh_wp2_zt > 0 ) then
        lh_wm  = compute_sample_mean( nzt, num_samples, lh_sample_point_weights, &
                                      real( X_nl_all_levs(:,:,iiPDF_w), kind = core_rknd) )
        ! Switch back to using stat_update_var once the code is generalized
        ! to pass in the number of vertical levels.
!        call stat_update_var( stats_metadata%ilh_wm, lh_wm, stats_lh_zt )
        do k = 1, nzt
           call stat_update_var_pt( stats_metadata%ilh_wm, k, lh_wm(k), stats_lh_zt )
        enddo ! k = 1, nzt
      end if

      if (   stats_metadata%ilh_rrm  + stats_metadata%ilh_Nrm &
           + stats_metadata%ilh_rim  + stats_metadata%ilh_Nim &
           + stats_metadata%ilh_rsm  + stats_metadata%ilh_Nsm &
           + stats_metadata%ilh_rgm  + stats_metadata%ilh_Ngm &
           + stats_metadata%ilh_Ncnm + stats_metadata%ilh_Ncm > 0 ) then

        lh_hydromet = 0._core_rknd
        call copy_X_nl_into_hydromet_all_pts( nzt, pdf_dim, num_samples,     & ! In
                                              X_nl_all_levs,                & ! In
                                              hydromet_dim, hm_metadata, & ! In
                                              lh_hydromet,                  & ! In
                                              hydromet_all_points,          & ! Out
                                              Ncn_all_points                ) ! Out

        ! Get rid of an annoying compiler warning.
        ivar = 1
        ivar = ivar

        forall ( ivar = 1:hydromet_dim )
          lh_hydromet(:,ivar) = compute_sample_mean( nzt, num_samples, lh_sample_point_weights,&
                                                     hydromet_all_points(:,:,ivar) )
        end forall ! 1..hydromet_dim

      end if

      ! Switch back to using stat_update_var once the code is generalized
      ! to pass in the number of vertical levels.
      if ( stats_metadata%ilh_Ncnm > 0 ) then
        lh_Ncnm = compute_sample_mean( nzt, num_samples, lh_sample_point_weights, &
                                       Ncn_all_points(:,:) )
!        call stat_update_var( stats_metadata%ilh_Ncnm, lh_Ncnm, stats_lh_zt )
        do k = 1, nzt
           call stat_update_var_pt( stats_metadata%ilh_Ncnm, k, lh_Ncnm(k), stats_lh_zt )
        enddo ! k = 1, nzt
      end if

      if ( stats_metadata%ilh_Ncm > 0 ) then
        lh_Ncm = compute_sample_mean( nzt, num_samples, lh_sample_point_weights, &
                                      lh_Nc_clipped(:,:) )
!        call stat_update_var( stats_metadata%ilh_Ncm, lh_Ncm, stats_lh_zt )
        do k = 1, nzt
           call stat_update_var_pt( stats_metadata%ilh_Ncm, k, lh_Ncm(k), stats_lh_zt )
        enddo ! k = 1, nzt
      end if

      ! Latin hypercube estimate of cloud fraction
      if ( stats_metadata%ilh_cloud_frac > 0 ) then
        lh_cloud_frac(:) = zero
        do sample = 1, num_samples
          where ( X_nl_all_levs(sample,:,iiPDF_chi) > zero )
            lh_cloud_frac(:) = lh_cloud_frac(:) + one * lh_sample_point_weights(sample,:)
          end where
        end do
        lh_cloud_frac(:) = lh_cloud_frac(:) / real( num_samples, kind = core_rknd )

!        call stat_update_var( stats_metadata%ilh_cloud_frac, lh_cloud_frac, stats_lh_zt )
        do k = 1, nzt
           call stat_update_var_pt( stats_metadata%ilh_cloud_frac, k, lh_cloud_frac(k), &
                                    stats_lh_zt )
        enddo ! k = 1, nzt
      end if

      ! Sample of lh_cloud_frac that is not weighted
      if ( stats_metadata%ilh_cloud_frac_unweighted > 0 ) then
        lh_cloud_frac(:) = zero
        do sample = 1, num_samples
          where ( X_nl_all_levs(sample,:,iiPDF_chi) > zero )
            lh_cloud_frac(:) = lh_cloud_frac(:) + one
          end where
        end do
        lh_cloud_frac(:) = lh_cloud_frac(:) / real( num_samples, kind = core_rknd )

!        call stat_update_var( stats_metadata%ilh_cloud_frac_unweighted, lh_cloud_frac, stats_lh_zt)
        do k = 1, nzt
           call stat_update_var_pt( stats_metadata%ilh_cloud_frac_unweighted, k, &
                                    lh_cloud_frac(k), stats_lh_zt )
        enddo ! k = 1, nzt
      end if

      ! Latin hypercube estimate of chi
      if ( stats_metadata%ilh_chi > 0 ) then
        lh_chi(1:nzt) &
        = compute_sample_mean( nzt, num_samples, lh_sample_point_weights, &
                               X_nl_all_levs(1:num_samples, 1:nzt, iiPDF_chi) )
!        call stat_update_var( stats_metadata%ilh_chi, lh_chi, stats_lh_zt )
        do k = 1, nzt
           call stat_update_var_pt( stats_metadata%ilh_chi, k, lh_chi(k), stats_lh_zt )
        enddo ! k = 1, nzt
      end if

      ! Latin hypercube estimate of variance of chi
      if ( stats_metadata%ilh_chip2 > 0 ) then
        lh_chip2(1:nzt) &
        = compute_sample_variance( nzt, num_samples, &
                                   X_nl_all_levs(:,:,iiPDF_chi), &
                                   lh_sample_point_weights, lh_chi(1:nzt) )
!        call stat_update_var( stats_metadata%ilh_chip2, lh_chip2, stats_lh_zt )
        do k = 1, nzt
           call stat_update_var_pt( stats_metadata%ilh_chip2, k, lh_chip2(k), stats_lh_zt )
        enddo ! k = 1, nzt
      end if

      ! Latin hypercube estimate of eta
      if ( stats_metadata%ilh_eta > 0 ) then
        lh_eta(1:nzt) &
        = compute_sample_mean( nzt, num_samples, lh_sample_point_weights, &
                               X_nl_all_levs(1:num_samples, 1:nzt, iiPDF_eta) )

!        call stat_update_var( stats_metadata%ilh_eta, lh_eta, stats_lh_zt )
        do k = 1, nzt
           call stat_update_var_pt( stats_metadata%ilh_eta, k, lh_eta(k), stats_lh_zt )
        enddo ! k = 1, nzt
      end if

      if ( stats_metadata%ilh_wp2_zt > 0 ) then
        ! Compute the variance of vertical velocity
        lh_wp2_zt = compute_sample_variance( nzt, num_samples, &
                                             X_nl_all_levs(:,:,iiPDF_w), &
                                             lh_sample_point_weights, lh_wm )
!        call stat_update_var( stats_metadata%ilh_wp2_zt, lh_wp2_zt, stats_lh_zt )
        do k = 1, nzt
           call stat_update_var_pt( stats_metadata%ilh_wp2_zt, k, lh_wp2_zt(k), stats_lh_zt )
        enddo ! k = 1, nzt
      end if

      if ( stats_metadata%ilh_rcp2_zt  > 0 ) then
        ! Compute the variance of cloud water mixing ratio
        lh_rcp2_zt = compute_sample_variance &
                     ( nzt, num_samples, lh_rc_clipped, &
                       lh_sample_point_weights, lh_rcm )
!        call stat_update_var( stats_metadata%ilh_rcp2_zt, lh_rcp2_zt, stats_lh_zt )
        do k = 1, nzt
           call stat_update_var_pt( stats_metadata%ilh_rcp2_zt, k, lh_rcp2_zt(k), &
                                    stats_lh_zt )
        enddo ! k = 1, nzt
      end if

      if ( stats_metadata%ilh_rtp2_zt > 0 ) then
        ! Compute the variance of total water
        lh_rtp2_zt = compute_sample_variance &
                     ( nzt, num_samples, &
                       lh_rt_clipped, lh_sample_point_weights, &
                       lh_rvm+lh_rcm )
!        call stat_update_var( stats_metadata%ilh_rtp2_zt, lh_rtp2_zt, stats_lh_zt )
        do k = 1, nzt
           call stat_update_var_pt( stats_metadata%ilh_rtp2_zt, k, lh_rtp2_zt(k), &
                                    stats_lh_zt )
        enddo ! k = 1, nzt
      end if

      if ( stats_metadata%ilh_thlp2_zt > 0 ) then
        ! Compute the variance of liquid potential temperature
        lh_thlp2_zt = compute_sample_variance( nzt, num_samples, &
                        lh_thl_clipped, lh_sample_point_weights, &
                        lh_thlm )
!        call stat_update_var( stats_metadata%ilh_thlp2_zt, lh_thlp2_zt, stats_lh_zt )
        do k = 1, nzt
           call stat_update_var_pt( stats_metadata%ilh_thlp2_zt, k, lh_thlp2_zt(k), &
                                    stats_lh_zt )
        enddo ! k = 1, nzt
      end if

      ! Compute the variance of rain water mixing ratio
      if ( iirr > 0 .and. stats_metadata%ilh_rrp2_zt > 0 ) then
        lh_rrp2_zt = compute_sample_variance &
                        ( nzt, num_samples, hydromet_all_points(:,:,iirr), &
                          lh_sample_point_weights, lh_hydromet(:,iirr) )
!        call stat_update_var( stats_metadata%ilh_rrp2_zt, lh_rrp2_zt, stats_lh_zt )
        do k = 1, nzt
           call stat_update_var_pt( stats_metadata%ilh_rrp2_zt, k, lh_rrp2_zt(k), &
                                    stats_lh_zt )
        enddo ! k = 1, nzt
      end if

      ! Compute the variance of cloud nuclei concentration (simplifed)
      if ( iiPDF_Ncn > 0 .and. stats_metadata%ilh_Ncnp2_zt > 0 ) then
        lh_Ncnp2_zt = compute_sample_variance &
                      ( nzt, num_samples, Ncn_all_points(:,:), &
                        lh_sample_point_weights, lh_Ncnm(:) )
!        call stat_update_var( stats_metadata%ilh_Ncnp2_zt, lh_Ncnp2_zt, stats_lh_zt )
        do k = 1, nzt
           call stat_update_var_pt( stats_metadata%ilh_Ncnp2_zt, k, lh_Ncnp2_zt(k), &
                                    stats_lh_zt )
        enddo ! k = 1, nzt
      end if

      ! Compute the variance of cloud droplet concentration
      if ( stats_metadata%ilh_Ncp2_zt > 0 ) then
        lh_Ncp2_zt = compute_sample_variance &
                     ( nzt, num_samples, lh_Nc_clipped(:,:), &
                       lh_sample_point_weights, lh_Ncm(:) )
!        call stat_update_var( stats_metadata%ilh_Ncp2_zt, lh_Ncp2_zt, stats_lh_zt )
        do k = 1, nzt
           call stat_update_var_pt( stats_metadata%ilh_Ncp2_zt, k, lh_Ncp2_zt(k), &
                                    stats_lh_zt )
        enddo ! k = 1, nzt
      end if

      ! Compute the variance of rain droplet number concentration
      if ( iiNr > 0 .and. stats_metadata%ilh_Nrp2_zt > 0 ) then
        lh_Nrp2_zt = compute_sample_variance( nzt, num_samples, hydromet_all_points(:,:,iiNr),&
                                              lh_sample_point_weights, lh_hydromet(:,iiNr) )
!        call stat_update_var( stats_metadata%ilh_Nrp2_zt, lh_Nrp2_zt, stats_lh_zt )
        do k = 1, nzt
           call stat_update_var_pt( stats_metadata%ilh_Nrp2_zt, k, lh_Nrp2_zt(k), &
                                    stats_lh_zt )
        enddo ! k = 1, nzt
      end if

      ! Averages of points being fed into the microphysics
      ! These are for diagnostic purposes, and are not needed for anything
      if ( iirr > 0 ) then
!        call stat_update_var( stats_metadata%ilh_rrm, lh_hydromet(:,iirr), stats_lh_zt )
        do k = 1, nzt
           call stat_update_var_pt( stats_metadata%ilh_rrm, k, lh_hydromet(k,iirr), &
                                    stats_lh_zt )
        enddo ! k = 1, nzt
      end if
      if ( iiNr > 0 ) then
!        call stat_update_var( stats_metadata%ilh_Nrm, lh_hydromet(:,iiNr), stats_lh_zt )
        do k = 1, nzt
           call stat_update_var_pt( stats_metadata%ilh_Nrm, k, lh_hydromet(k,iiNr), &
                                    stats_lh_zt )
        enddo ! k = 1, nzt
      end if
      if ( iiri > 0 ) then
!        call stat_update_var( stats_metadata%ilh_rim, lh_hydromet(:,iiri), stats_lh_zt )
        do k = 1, nzt
           call stat_update_var_pt( stats_metadata%ilh_rim, k, lh_hydromet(k,iiri), &
                                    stats_lh_zt )
        enddo ! k = 1, nzt
      end if
      if ( iiNi > 0 ) then
!        call stat_update_var( stats_metadata%ilh_Nim, lh_hydromet(:,iiNi), stats_lh_zt )
        do k = 1, nzt
           call stat_update_var_pt( stats_metadata%ilh_Nim, k, lh_hydromet(k,iiNi), &
                                    stats_lh_zt )
        enddo ! k = 1, nzt
      end if
      if ( iirs > 0 ) then
!        call stat_update_var( stats_metadata%ilh_rsm, lh_hydromet(:,iirs), stats_lh_zt )
        do k = 1, nzt
           call stat_update_var_pt( stats_metadata%ilh_rsm, k, lh_hydromet(k,iirs), &
                                    stats_lh_zt )
        enddo ! k = 1, nzt
      end if
      if ( iiNs > 0 ) then
!        call stat_update_var( stats_metadata%ilh_Nsm, lh_hydromet(:,iiNs), stats_lh_zt )
        do k = 1, nzt
           call stat_update_var_pt( stats_metadata%ilh_Nsm, k, lh_hydromet(k,iiNs), &
                                    stats_lh_zt )
        enddo ! k = 1, nzt
      end if
      if ( iirg > 0 ) then
!        call stat_update_var( stats_metadata%ilh_rgm, lh_hydromet(:,iirg), stats_lh_zt )
        do k = 1, nzt
           call stat_update_var_pt( stats_metadata%ilh_rgm, k, lh_hydromet(k,iirg), &
                                    stats_lh_zt )
        enddo ! k = 1, nzt
      end if
      if ( iiNg > 0 ) then
!        call stat_update_var( stats_metadata%ilh_Ngm, lh_hydromet(:,iiNg), stats_lh_zt )
        do k = 1, nzt
           call stat_update_var_pt( stats_metadata%ilh_Ngm, k, lh_hydromet(k,iiNg), &
                                    stats_lh_zt )
        enddo ! k = 1, nzt
      end if

    end if ! stats_metadata%l_stats_samp

    return
  end subroutine stats_accumulate_lh_api

  !-----------------------------------------------------------------------
  subroutine stats_accumulate_uniform_lh( nzt, num_samples, ngrdcol, l_in_precip_all_levs, &
                                          X_mixt_comp_all_levs, X_u_chi_all_levs, pdf_params, &
                                          lh_sample_point_weights, k_lh_start, &
                                          stats_metadata, &
                                          stats_lh_zt, stats_lh_sfc )

  ! Description:
  !   Samples statistics that cannot be deduced from the normal-lognormal
  !   SILHS sample (X_nl_all_levs)

  ! References:
  !   None
  !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd            ! Constant

    use stats_type_utilities, only: &
      !stat_update_var,   & ! Procedure(s)
      stat_update_var_pt

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

    use stats_type, only: &
      stats ! Type

    use stats_variables, only: &
      stats_metadata_type

    implicit none

    !------------------------------- Input Variables -------------------------------
    integer, intent(in) :: &
      nzt,          & ! Number of vertical levels
      num_samples, & ! Number of SILHS sample points
      ngrdcol        ! Number of grid columns

    logical, dimension(ngrdcol,num_samples,nzt), intent(in) :: &
      l_in_precip_all_levs ! Boolean variables indicating whether a sample is in
                           ! precipitation at a given height level

    integer, dimension(ngrdcol,num_samples,nzt), intent(in) :: &
      X_mixt_comp_all_levs ! Integers indicating which mixture component a
                           ! sample is in at a given height level

    real( kind = core_rknd ), dimension(ngrdcol,num_samples,nzt), intent(in) :: &
      X_u_chi_all_levs     ! Uniform value of chi

    type(pdf_parameter), intent(in) :: &
      pdf_params           ! The official PDF parameters!

    real( kind = core_rknd ), dimension(ngrdcol,num_samples,nzt), intent(in) :: &
      lh_sample_point_weights ! The weight of each sample

    integer, dimension(ngrdcol), intent(in) :: &
      k_lh_start           ! Vertical level for sampling preferentially within       [-]
                           ! cloud
    
    type (stats_metadata_type), intent(in) :: &
      stats_metadata

    !------------------------------- InOut Variables -------------------------------
    type(stats), dimension(ngrdcol), intent(inout) :: &
      stats_lh_zt, &
      stats_lh_sfc

    !------------------------------- Local Variables -------------------------------
    type(importance_category_type), dimension(num_importance_categories) :: &
      importance_categories

    real( kind = core_rknd ), dimension(nzt) :: &
      lh_precip_frac, &
      lh_mixt_frac

    real( kind = core_rknd ), dimension(num_samples,nzt) :: &
      int_in_precip, & ! '1' for samples in precipitation, '0' otherwise
      int_mixt_comp    ! '1' for samples in the first PDF component, '0' otherwise

    real( kind = core_rknd ), dimension(num_samples,nzt) :: &
      one_weights

    real( kind = core_rknd ), dimension(nzt,num_importance_categories) :: &
      lh_samp_frac

    real( kind = core_rknd ) :: &
      cloud_frac_i

    logical :: &
      l_in_cloud, &
      l_in_comp_1

    integer, dimension(num_importance_categories) :: &
      category_counts  ! Count of number of samples in each importance category

    integer :: k, isample, icategory, i

    !------------------------------- Begin Code -------------------------------
    
    do i = 1, ngrdcol

      ! Switch back to using stat_update_var once the code is generalized
      ! to pass in the number of vertical levels.

      ! Estimate of lh_precip_frac
      if ( stats_metadata%ilh_precip_frac > 0 ) then
        where ( l_in_precip_all_levs(i,:,:) )
          int_in_precip = 1.0_core_rknd
        else where
          int_in_precip = 0.0_core_rknd
        end where
        lh_precip_frac(:) = compute_sample_mean( nzt, num_samples, lh_sample_point_weights(i,:,:), &
                                                 int_in_precip )
  !      call stat_update_var( stats_metadata%ilh_precip_frac, lh_precip_frac, stats_lh_zt )
        do k = 1, nzt
           call stat_update_var_pt( stats_metadata%ilh_precip_frac, k, lh_precip_frac(k), &
                                    stats_lh_zt(i) )
        enddo ! k = 1, nzt
      end if

      ! Unweighted estimate of lh_precip_frac
      if ( stats_metadata%ilh_precip_frac_unweighted > 0 ) then
        where ( l_in_precip_all_levs(i,:,:) )
          int_in_precip = 1.0_core_rknd
        else where
          int_in_precip = 0.0_core_rknd
        end where
        one_weights = one
        lh_precip_frac(:) = compute_sample_mean( nzt, num_samples, one_weights, &
                                                 int_in_precip )
  !      call stat_update_var(stats_metadata%ilh_precip_frac_unweighted,lh_precip_frac,stats_lh_zt)
        do k = 1, nzt
           call stat_update_var_pt( stats_metadata%ilh_precip_frac_unweighted, k, &
                                    lh_precip_frac(k), stats_lh_zt(i) )
        enddo ! k = 1, nzt
      end if

      ! Estimate of lh_mixt_frac
      if ( stats_metadata%ilh_mixt_frac > 0 ) then
        where ( X_mixt_comp_all_levs(i,:,:) == 1 )
          int_mixt_comp = 1.0_core_rknd
        else where
          int_mixt_comp = 0.0_core_rknd
        end where
        lh_mixt_frac(:) = compute_sample_mean( nzt, num_samples, lh_sample_point_weights(i,:,:), &
                                               int_mixt_comp )
  !      call stat_update_var( stats_metadata%ilh_mixt_frac, lh_mixt_frac, stats_lh_zt )
        do k = 1, nzt
           call stat_update_var_pt( stats_metadata%ilh_mixt_frac, k, lh_mixt_frac(k), &
                                    stats_lh_zt(i) )
        enddo ! k = 1, nzt
      end if

      ! Unweighted estimate of lh_mixt_frac
      if ( stats_metadata%ilh_mixt_frac_unweighted > 0 ) then
        where ( X_mixt_comp_all_levs(i,:,:) == 1 )
          int_mixt_comp = 1.0_core_rknd
        else where
          int_mixt_comp = 0.0_core_rknd
        end where
        one_weights = one
        lh_mixt_frac(:) = compute_sample_mean( nzt, num_samples, one_weights, &
                                               int_mixt_comp )
  !      call stat_update_var( stats_metadata%ilh_mixt_frac_unweighted, lh_mixt_frac, stats_lh_zt )
        do k = 1, nzt
           call stat_update_var_pt( stats_metadata%ilh_mixt_frac_unweighted, k, &
                                    lh_mixt_frac(k), stats_lh_zt(i) )
        enddo ! k = 1, nzt
      end if

      ! k_lh_start is an integer, so it would be more appropriate to sample it
      ! as an integer, but as far as I can tell our current sampling
      ! infrastructure mainly supports sampling real numbers.
      call stat_update_var_pt( stats_metadata%ik_lh_start, 1, &
                               real( k_lh_start(i), kind=core_rknd ), stats_lh_sfc(i) )

      if ( allocated( stats_metadata%ilh_samp_frac_category ) ) then
        if ( stats_metadata%ilh_samp_frac_category(1) > 0 ) then

          importance_categories = define_importance_categories( )

          do k=1, nzt
            category_counts(:) = 0

            do isample=1, num_samples

              if ( X_mixt_comp_all_levs(i,isample,k) == 1 ) then
                l_in_comp_1 = .true.
                cloud_frac_i = pdf_params%cloud_frac_1(i,k)
              else
                l_in_comp_1 = .false.
                cloud_frac_i = pdf_params%cloud_frac_2(i,k)
              end if

              l_in_cloud = X_u_chi_all_levs(i,isample,k) > (one - cloud_frac_i)

              do icategory=1, num_importance_categories
                if ( (l_in_cloud .eqv. importance_categories(icategory)%l_in_cloud) .and. &
                     (l_in_precip_all_levs(i,isample,k) .eqv. importance_categories(icategory)%&
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

          end do ! k=1, nzt

          ! Microphysics is not run at lower level
          lh_samp_frac(1,:) = zero

          do icategory=1, num_importance_categories
  !          call stat_update_var( stats_metadata%ilh_samp_frac_category(icategory), &
  !                                lh_samp_frac(:,icategory), stats_lh_zt )
            do k = 1, nzt
               call stat_update_var_pt( stats_metadata%ilh_samp_frac_category(icategory), k, &
                                        lh_samp_frac(k,icategory), stats_lh_zt(i) )
            enddo ! k = 1, nzt
          end do ! icategory=1, num_importance_categories

        end if ! stats_metadata%ilh_samp_frac_category(1) > 0
      end if ! allocated( stats_metadata%ilh_samp_frac_category )
      
    end do

    return
  end subroutine stats_accumulate_uniform_lh

  !-----------------------------------------------------------------------------
  subroutine copy_X_nl_into_hydromet_all_pts( nzt, pdf_dim, num_samples, &
                                              X_nl_all_levs, &
                                              hydromet_dim, hm_metadata, &
                                              hydromet, &
                                              hydromet_all_points, &
                                              Ncn_all_points )

  ! Description:
  !   Copy the points from the latin hypercube sample to an array with just the
  !   hydrometeors
  ! References:
  !   None
  !-----------------------------------------------------------------------------
    use corr_varnce_module, only: &
      hm_metadata_type

    use clubb_precision, only: &
      core_rknd

    implicit none

    integer, intent(in) :: &
      nzt,             & ! Number of vertical levels
      pdf_dim,        & ! Number of variates
      num_samples,    & ! Number of calls to microphysics
      hydromet_dim      ! Number of hydrometeor species

    type (hm_metadata_type), intent(in) :: &
      hm_metadata

    real( kind = core_rknd ), dimension(num_samples,nzt,pdf_dim), intent(in) :: &
      X_nl_all_levs ! Sample that is transformed ultimately to normal-lognormal

    real( kind = core_rknd ), dimension(nzt,hydromet_dim), intent(in) :: &
      hydromet ! Hydrometeor species    [units vary]

    real( kind = core_rknd ), dimension(num_samples,nzt,hydromet_dim), intent(out) :: &
      hydromet_all_points ! Hydrometeor species    [units vary]

    real( kind = core_rknd ), dimension(num_samples,nzt), intent(out) :: &
      Ncn_all_points    ! Cloud nuclei conc. (simplified); Nc=Ncn*H(chi)  [#/kg]

    integer :: sample, ivar

    do sample = 1, num_samples
      ! Copy the sample points into the temporary arrays
      do ivar = 1, hydromet_dim, 1
        if ( ivar == hm_metadata%iirr .and. hm_metadata%iiPDF_rr > 0 ) then
          ! Use a sampled value of rain water mixing ratio
          hydromet_all_points(sample,:,ivar) = &
            real( X_nl_all_levs(sample,:,hm_metadata%iiPDF_rr), kind = core_rknd )

        else if ( ivar == hm_metadata%iirs .and. hm_metadata%iiPDF_rs > 0 ) then
          ! Use a sampled value of rain water mixing ratio
          hydromet_all_points(sample,:,ivar) = &
            real( X_nl_all_levs(sample,:,hm_metadata%iiPDF_rs), kind = core_rknd )

        else if ( ivar == hm_metadata%iiri .and. hm_metadata%iiPDF_ri > 0 ) then
          ! Use a sampled value of rain water mixing ratio
          hydromet_all_points(sample,:,ivar) = &
            real( X_nl_all_levs(sample,:,hm_metadata%iiPDF_ri), kind = core_rknd )

        else if ( ivar == hm_metadata%iirg .and. hm_metadata%iiPDF_rg > 0 ) then
          ! Use a sampled value of rain water mixing ratio
          hydromet_all_points(sample,:,ivar) = &
            real( X_nl_all_levs(sample,:,hm_metadata%iiPDF_rg), kind = core_rknd )

        else if ( ivar == hm_metadata%iiNr .and. hm_metadata%iiPDF_Nr > 0 ) then
          ! Use a sampled value of rain droplet number concentration
          hydromet_all_points(sample,:,ivar) = &
            real( X_nl_all_levs(sample,:,hm_metadata%iiPDF_Nr), kind = core_rknd )

        else if ( ivar == hm_metadata%iiNs .and. hm_metadata%iiPDF_Ns > 0 ) then
          ! Use a sampled value of rain droplet number concentration
          hydromet_all_points(sample,:,ivar) = &
            real( X_nl_all_levs(sample,:,hm_metadata%iiPDF_Ns), kind = core_rknd )

        else if ( ivar == hm_metadata%iiNg .and. hm_metadata%iiPDF_Ng > 0 ) then
          ! Use a sampled value of rain droplet number concentration
          hydromet_all_points(sample,:,ivar) = &
            real( X_nl_all_levs(sample,:,hm_metadata%iiPDF_Ng), kind = core_rknd )

        else if ( ivar == hm_metadata%iiNi .and. hm_metadata%iiPDF_Ni > 0 ) then
          ! Use a sampled value of rain droplet number concentration
          hydromet_all_points(sample,:,ivar) = &
            real( X_nl_all_levs(sample,:,hm_metadata%iiPDF_Ni), kind = core_rknd )

        else ! Use the mean field, rather than a sample point
          ! This is the case for hail and graupel in the Morrison microphysics
          ! currently -dschanen 23 March 2010
          hydromet_all_points(sample,:,ivar) = hydromet(:,ivar)

        end if
      end do ! 1..hydromet_dim
      ! Copy Ncn into Ncn all points
      if ( hm_metadata%iiPDF_Ncn > 0 ) then
        Ncn_all_points(sample,:) = &
          real( X_nl_all_levs(sample,:,hm_metadata%iiPDF_Ncn), kind=core_rknd )
      end if
    end do ! 1..num_samples

    return
  end subroutine copy_X_nl_into_hydromet_all_pts
  !-----------------------------------------------------------------------------

#endif /* SILHS */

end module latin_hypercube_driver_module
