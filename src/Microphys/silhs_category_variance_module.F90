!-------------------------------------------------------------------------------
! $Id$
!===============================================================================
module silhs_category_variance_module

  implicit none

  public :: silhs_category_variance_driver

  private ! Set Default Scope

contains

  !-----------------------------------------------------------------------------
  subroutine silhs_category_variance_driver( &
               ngrdcol, nzt, num_samples, pdf_dim, hydromet_dim, hm_metadata, &
               X_nl_all_levs, X_mixt_comp_all_levs, lh_hydromet_mc_all, &
               lh_sample_point_weights, pdf_params, precip_fracs, stats )

    ! Description:
    !   Compute the variance of a microphysics variable in each importance
    !   category for all model columns.
    !---------------------------------------------------------------------------

    use clubb_precision, only: core_rknd
    use corr_varnce_module, only: hm_metadata_type
    use pdf_parameter_module, only: pdf_parameter
    use hydromet_pdf_parameter_module, only: precipitation_fractions
    use silhs_importance_sample_module, only: num_importance_categories
    use stats_netcdf, only: stats_type, stats_update

    implicit none

    !--------------------------- Input Variables ---------------------------
    integer, intent(in) :: &
      ngrdcol, &      ! Number of model columns
      nzt,             &      ! Number of height levels
      num_samples,     &      ! Number of SILHS sample points
      pdf_dim,         &      ! Number of variates in X_nl
      hydromet_dim            ! Number of elements of hydromet array
    type(hm_metadata_type), intent(in) :: hm_metadata

    real( kind = core_rknd ), dimension(ngrdcol,num_samples,nzt,pdf_dim), intent(in) :: &
      X_nl_all_levs           ! SILHS samples at all height levels
    integer, dimension(ngrdcol,num_samples,nzt), intent(in) :: &
      X_mixt_comp_all_levs    ! Mixture component (1 or 2) of each sample point
    real( kind = core_rknd ), &
      dimension(ngrdcol,num_samples,nzt,hydromet_dim), intent(in) :: &
      lh_hydromet_mc_all   ! Tendencies of hydrometeors at all sample points
    real( kind = core_rknd ), dimension(ngrdcol,num_samples,nzt), intent(in) :: &
      lh_sample_point_weights ! Weight of SILHS sample points

    type(pdf_parameter), intent(in) :: pdf_params ! The PDF parameters!
    type(precipitation_fractions), intent(in) :: precip_fracs ! Precipitation fractions [-]
    !--------------------------- InOut Variables ---------------------------
    type(stats_type), intent(inout) :: stats

    !--------------------------- Local Variables ---------------------------
    real( kind = core_rknd ), dimension(ngrdcol,num_samples,nzt) :: samples_all
    real( kind = core_rknd ), dimension(ngrdcol,nzt,num_importance_categories) :: &
      root_weight_mean_sq_cat

    integer :: icat
    character(len=32) :: cat_name
    character(len=2) :: cat_num

    !--------------------------- Begin Code ---------------------------

    ! Sample rrm_mc from the hydrometeor tendencies generated for each SILHS point.
    samples_all = lh_hydromet_mc_all(:,:,:,hm_metadata%iirr)

    call silhs_sample_category_variance( &
           ngrdcol, nzt, num_samples, pdf_dim, X_nl_all_levs, & ! In
           X_mixt_comp_all_levs, samples_all, lh_sample_point_weights, & ! In
           pdf_params, precip_fracs, hm_metadata, & ! In
           root_weight_mean_sq_cat ) ! Out

    if ( stats%l_sample ) then
      do icat = 1, num_importance_categories
        write(cat_num,'(I1)') icat
        cat_name = "silhs_var_cat_" // trim(cat_num)
        call stats_update( cat_name, root_weight_mean_sq_cat(:,:,icat), stats )
      enddo
    endif

  end subroutine silhs_category_variance_driver

  !-----------------------------------------------------------------------------
  subroutine silhs_sample_category_variance( &
               ngrdcol, nzt, num_samples, pdf_dim, X_nl_all_levs, & ! In
               X_mixt_comp_all_levs, samples_all, & ! In
               lh_sample_point_weights, pdf_params, precip_fracs, & ! In
               hm_metadata, root_weight_mean_sq_cat ) ! Out

    ! Description:
    !   Compute category-conditioned root-mean-square values in every column.
    !---------------------------------------------------------------------------

    use clubb_precision, only: core_rknd
    use constants_clubb, only: one, zero
    use silhs_importance_sample_module, only: &
      num_importance_categories, &
      define_importance_categories, &
      importance_category_type
    use pdf_parameter_module, only: pdf_parameter
    use hydromet_pdf_parameter_module, only: precipitation_fractions
    use corr_varnce_module, only: hm_metadata_type

    implicit none

    !---------------------- Input Variables ----------------------
    integer, intent(in) :: &
      ngrdcol, &      ! Number of model columns
      nzt,             & ! Number of height levels
      num_samples,     & ! Number of SILHS sample points
      pdf_dim            ! Number of variates in X_nl
    real( kind = core_rknd ), dimension(ngrdcol,num_samples,nzt,pdf_dim), intent(in) :: &
      X_nl_all_levs        ! SILHS samples at all height levels
    integer, dimension(ngrdcol,num_samples,nzt), intent(in) :: &
      X_mixt_comp_all_levs  ! Mixture component (1 or 2) of each sample point
    real( kind = core_rknd ), dimension(ngrdcol,num_samples,nzt), intent(in) :: &
      samples_all, &         ! Sample points of variable to compute variance of
      lh_sample_point_weights ! Weight of SILHS sample points
    type(pdf_parameter), intent(in) :: pdf_params ! The PDF parameters!
    type(precipitation_fractions), intent(in) :: precip_fracs ! Precipitation fractions [-]
    type(hm_metadata_type), intent(in) :: hm_metadata

    real( kind = core_rknd ), &
      dimension(ngrdcol,nzt,num_importance_categories), intent(out) :: &
      root_weight_mean_sq_cat

    type(importance_category_type), dimension(num_importance_categories) :: &
      importance_categories
    integer, dimension(ngrdcol,num_samples) :: &
      int_sample_category     ! Category of each sample point
    real( kind = core_rknd ), dimension(ngrdcol,num_importance_categories) :: &
      category_real_probs     ! PDF probability of each category

    logical, dimension(ngrdcol) :: &
      l_in_cloud, &
      l_in_precip, &
      l_in_component_1

    integer, dimension(ngrdcol) :: found_category_index

    integer :: i, isample, icat, k

    !-------------------------- Begin Code --------------------------

    if ( hm_metadata%iiPDF_rr == -1 ) then
      error stop "iiPDF_rr must be greater than zero for the category sampler to work."
    endif

    ! Keep the sample classifications consistent with the categories returned
    ! by define_importance_categories.
    root_weight_mean_sq_cat = zero
    importance_categories = define_importance_categories()

    do k = 1, nzt
      do isample = 1, num_samples
        do i = 1, ngrdcol
          l_in_cloud(i) = X_nl_all_levs(i,isample,k,hm_metadata%iiPDF_chi) >= zero
          l_in_precip(i) = X_nl_all_levs(i,isample,k,hm_metadata%iiPDF_rr) > zero
          l_in_component_1(i) = X_mixt_comp_all_levs(i,isample,k) == 1
        enddo

        found_category_index = -1
        do icat = 1, num_importance_categories
          do i = 1, ngrdcol
            if ( found_category_index(i) == -1 .and. &
                 ( importance_categories(icat)%l_in_cloud .eqv. l_in_cloud(i) ) .and. &
                 ( importance_categories(icat)%l_in_precip .eqv. l_in_precip(i) ) .and. &
                 ( importance_categories(icat)%l_in_component_1 .eqv. &
                   l_in_component_1(i) ) ) then
              found_category_index(i) = icat
            endif
          enddo
        enddo

        do i = 1, ngrdcol
          if ( found_category_index(i) == -1 ) then
            error stop "Fatal error determining category in determine_sample_categories"
          endif

          int_sample_category(i,isample) = found_category_index(i)
        enddo
      enddo

      do icat = 1, num_importance_categories
        if ( importance_categories(icat)%l_in_component_1 ) then
          do i = 1, ngrdcol
            if ( importance_categories(icat)%l_in_cloud ) then
              if ( importance_categories(icat)%l_in_precip ) then
                category_real_probs(i,icat) = pdf_params%mixt_frac(i,k) &
                  * pdf_params%cloud_frac_1(i,k) * precip_fracs%precip_frac_1(i,k)
              else
                category_real_probs(i,icat) = pdf_params%mixt_frac(i,k) &
                  * pdf_params%cloud_frac_1(i,k) * ( one - precip_fracs%precip_frac_1(i,k) )
              endif
            else
              if ( importance_categories(icat)%l_in_precip ) then
                category_real_probs(i,icat) = pdf_params%mixt_frac(i,k) &
                  * ( one - pdf_params%cloud_frac_1(i,k) ) * precip_fracs%precip_frac_1(i,k)
              else
                category_real_probs(i,icat) = pdf_params%mixt_frac(i,k) &
                  * ( one - pdf_params%cloud_frac_1(i,k) ) &
                  * ( one - precip_fracs%precip_frac_1(i,k) )
              endif
            endif
          enddo
        else
          do i = 1, ngrdcol
            if ( importance_categories(icat)%l_in_cloud ) then
              if ( importance_categories(icat)%l_in_precip ) then
                category_real_probs(i,icat) = ( one - pdf_params%mixt_frac(i,k) ) &
                  * pdf_params%cloud_frac_2(i,k) * precip_fracs%precip_frac_2(i,k)
              else
                category_real_probs(i,icat) = ( one - pdf_params%mixt_frac(i,k) ) &
                  * pdf_params%cloud_frac_2(i,k) * ( one - precip_fracs%precip_frac_2(i,k) )
              endif
            else
              if ( importance_categories(icat)%l_in_precip ) then
                category_real_probs(i,icat) = ( one - pdf_params%mixt_frac(i,k) ) &
                  * ( one - pdf_params%cloud_frac_2(i,k) ) * precip_fracs%precip_frac_2(i,k)
              else
                category_real_probs(i,icat) = ( one - pdf_params%mixt_frac(i,k) ) &
                  * ( one - pdf_params%cloud_frac_2(i,k) ) &
                  * ( one - precip_fracs%precip_frac_2(i,k) )
              endif
            endif
          enddo
        endif
      enddo

      do isample = 1, num_samples
        do i = 1, ngrdcol
          icat = int_sample_category(i,isample)
          root_weight_mean_sq_cat(i,k,icat) = root_weight_mean_sq_cat(i,k,icat) + &
            lh_sample_point_weights(i,isample,k) * samples_all(i,isample,k)**2
        enddo
      enddo

      do icat = 1, num_importance_categories
        do i = 1, ngrdcol
          root_weight_mean_sq_cat(i,k,icat) = root_weight_mean_sq_cat(i,k,icat) &
            / real( num_samples, kind=core_rknd )
          if ( category_real_probs(i,icat) > zero ) then
            root_weight_mean_sq_cat(i,k,icat) = sqrt( &
              root_weight_mean_sq_cat(i,k,icat) / category_real_probs(i,icat) )
          else
            root_weight_mean_sq_cat(i,k,icat) = -999._core_rknd
          endif
        enddo
      enddo
    enddo

  end subroutine silhs_sample_category_variance

end module silhs_category_variance_module
