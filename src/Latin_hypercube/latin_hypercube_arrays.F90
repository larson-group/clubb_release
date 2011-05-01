! $Id$
module latin_hypercube_arrays

  implicit none

  public :: setup_corr_varnce_array, cleanup_latin_hypercube_arrays

  private

  integer, public :: d_variables
!omp threadprivate(d_variables)

  real, public, dimension(:), allocatable :: &
    xp2_on_xm2_array_cloud, &
    xp2_on_xm2_array_below

  real, public, dimension(:,:), allocatable :: &
    corr_array_cloud, &
    corr_array_below
!$omp threadprivate(xp2_on_xm2_array_cloud, xp2_on_xm2_array_below, &
!$omp   corr_array_cloud, corr_array_below)


  contains
!===============================================================================
  subroutine setup_corr_varnce_array( d_variables_in )
! Description:
!   Setup an array with the x'^2/xm^2 variables on the diagonal and the other
!   elements to be correlations between various variables.
! References:
!   None.
!-------------------------------------------------------------------------------

    use parameters_microphys, only: &
      rrp2_on_rrainm2_cloud, Nrp2_on_Nrm2_cloud, Ncp2_on_Ncm2_cloud, & ! Variables
      corr_rrNr_LL_cloud, corr_srr_NL_cloud, corr_sNr_NL_cloud, &
      corr_sNc_NL_cloud, &
      rrp2_on_rrainm2_below, Nrp2_on_Nrm2_below, Ncp2_on_Ncm2_below, &
      corr_rrNr_LL_below, corr_srr_NL_below, corr_sNr_NL_below, &
      corr_sNc_NL_below

    use parameters_microphys, only: &
      rsnowp2_on_rsnowm2_cloud, & ! Variables
      Nsnowp2_on_Nsnowm2_cloud, & 
      ricep2_on_ricem2_cloud, & 
      Nicep2_on_Nicem2_cloud, &
      rsnowp2_on_rsnowm2_below, & 
      Nsnowp2_on_Nsnowm2_below, & 
      ricep2_on_ricem2_below, & 
      Nicep2_on_Nicem2_below

    use parameters_microphys, only: &
      corr_srsnow_NL_cloud, &
      corr_sNsnow_NL_cloud, &
      corr_rsnowNsnow_LL_cloud, &
      corr_srice_NL_cloud, &
      corr_sNi_NL_cloud, &
      corr_riceNi_LL_cloud, &
      corr_wrice_NL_cloud, &
      corr_wNi_NL_cloud, &
      corr_wrsnow_NL_cloud, &
      corr_wNsnow_NL_cloud, &
      corr_sw_NN_cloud

    use parameters_microphys, only: &
      corr_srsnow_NL_below, &
      corr_sNsnow_NL_below, &
      corr_rsnowNsnow_LL_below, &
      corr_srice_NL_below, &
      corr_sNi_NL_below, &
      corr_riceNi_LL_below, &
      corr_wrice_NL_below, &
      corr_wNi_NL_below, &
      corr_wrsnow_NL_below, &
      corr_wNsnow_NL_below

    use array_index, only: &
      iiLH_s_mellor, & ! Variables
      iiLH_t_mellor, & 
      iiLH_w, & 
      iiLH_rrain, &
      iiLH_rsnow, &
      iiLH_rice, &
      iiLH_rgraupel, &
      iiLH_Nr, &
      iiLH_Nsnow, &
      iiLH_Ni, &
      iiLH_Ngraupel, &
      iiLH_Nc

    use matrix_operations, only: set_lower_triangular_matrix_sp ! Procedure(s)

    implicit none

    ! External
    intrinsic :: max, epsilon

    ! Constant Parameters
    real, parameter :: corr_s_t = 0.3

    ! Input Variables
    integer, intent(in) :: &
      d_variables_in ! Number of variates in the array

    integer :: i

    ! ---- Begin Code ----
    d_variables = d_variables_in

    allocate( corr_array_cloud(d_variables,d_variables) )
    allocate( corr_array_below(d_variables,d_variables) )

    allocate( xp2_on_xm2_array_cloud(d_variables) )
    allocate( xp2_on_xm2_array_below(d_variables) )

    ! Initializing to zero means that correlations we don't have
    ! (e.g. Nc and any variable other than s_mellor ) are assumed to be 0.
    corr_array_cloud(:,:) = 0.0 ! Initialize to 0
    corr_array_below(:,:) = 0.0 ! Initialize to 0

    xp2_on_xm2_array_cloud(:) = 0.0
    xp2_on_xm2_array_below(:) = 0.0

    ! Set main diagonal to 1
    do i = 1, d_variables
      corr_array_cloud(i,i) = 1.0
      corr_array_below(i,i) = 1.0
    end do

    ! Use a fixed value for the correlation between s and t.
    call set_lower_triangular_matrix_sp &
         ( d_variables, iiLH_s_mellor, iiLH_t_mellor, corr_s_t, &
           corr_array_cloud )
    call set_lower_triangular_matrix_sp &
         ( d_variables, iiLH_s_mellor, iiLH_t_mellor, corr_s_t, &
           corr_array_below )

    ! Set the value of the correlation between w and s.
    call set_lower_triangular_matrix_sp &
         ( d_variables, iiLH_s_mellor, iiLH_w, corr_sw_NN_cloud, &
           corr_array_cloud )

    call set_lower_triangular_matrix_sp &
         ( d_variables, iiLH_s_mellor, iiLH_w, corr_sw_NN_cloud, &
           corr_array_below )

    if ( iiLH_Nc > 0 ) then
      xp2_on_xm2_array_cloud(iiLH_Nc) = Ncp2_on_Ncm2_cloud
      call set_lower_triangular_matrix_sp &
           ( d_variables, iiLH_s_mellor, iiLH_Nc, corr_sNc_NL_cloud, &
             corr_array_cloud )
    end if

    if ( iiLH_rrain > 0 ) then
      xp2_on_xm2_array_cloud(iiLH_rrain) = rrp2_on_rrainm2_cloud
      call set_lower_triangular_matrix_sp &
           ( d_variables, iiLH_s_mellor, iiLH_rrain, corr_srr_NL_cloud, &
             corr_array_cloud )
      if ( iiLH_Nr > 0 ) then
        xp2_on_xm2_array_cloud(iiLH_Nr) = Nrp2_on_Nrm2_cloud
        call set_lower_triangular_matrix_sp &
             ( d_variables, iiLH_rrain, iiLH_Nr, corr_rrNr_LL_cloud, &
               corr_array_cloud )
        call set_lower_triangular_matrix_sp &
             ( d_variables, iiLH_s_mellor, iiLH_Nr, corr_sNr_NL_cloud, &
               corr_array_cloud )
      end if ! iiLH_Nr > 0
    end if ! iiLH_rrain > 0

    if ( iiLH_rsnow > 0 ) then
      xp2_on_xm2_array_cloud(iiLH_rsnow) = rsnowp2_on_rsnowm2_cloud

      ! Correlation with s
      call set_lower_triangular_matrix_sp &
           ( d_variables, iiLH_rsnow, iiLH_s_mellor, corr_srsnow_NL_cloud, &
             corr_array_cloud )

      ! Correlation with w
      call set_lower_triangular_matrix_sp &
           ( d_variables, iiLH_rsnow, iiLH_w, corr_wrsnow_NL_cloud, &
             corr_array_cloud )

      if ( iiLH_Nsnow > 0 ) then
        xp2_on_xm2_array_cloud(iiLH_Nsnow) = Nsnowp2_on_Nsnowm2_cloud

        call set_lower_triangular_matrix_sp &
             ( d_variables, iiLH_rsnow, iiLH_Nsnow, corr_rsnowNsnow_LL_cloud, &
               corr_array_cloud )

        call set_lower_triangular_matrix_sp &
             ( d_variables, iiLH_Nsnow, iiLH_s_mellor, corr_sNsnow_NL_cloud, &
               corr_array_cloud )

        call set_lower_triangular_matrix_sp &
             ( d_variables, iiLH_Nsnow, iiLH_w, corr_wNsnow_NL_cloud, &
               corr_array_cloud )

      end if ! iiLH_Nsnow > 0
    end if ! iiLH_rsnow > 0

    if ( iiLH_rice > 0 ) then
      xp2_on_xm2_array_cloud(iiLH_rice) = ricep2_on_ricem2_cloud

      call set_lower_triangular_matrix_sp &
           ( d_variables, iiLH_rice, iiLH_s_mellor, corr_srice_NL_cloud, &
             corr_array_cloud )

      call set_lower_triangular_matrix_sp &
           ( d_variables, iiLH_rice, iiLH_w, corr_wrice_NL_cloud, &
             corr_array_cloud )

      if ( iiLH_Ni > 0 ) then
        xp2_on_xm2_array_cloud(iiLH_Ni) = Nicep2_on_Nicem2_cloud

        call set_lower_triangular_matrix_sp &
             ( d_variables, iiLH_rice, iiLH_Ni, corr_riceNi_LL_cloud, &
               corr_array_cloud )

        call set_lower_triangular_matrix_sp &
             ( d_variables, iiLH_Ni, iiLH_s_mellor, corr_sNi_NL_cloud, &
               corr_array_cloud )

        call set_lower_triangular_matrix_sp &
             ( d_variables, iiLH_Ni, iiLH_w, corr_wNi_NL_cloud, &
               corr_array_cloud )
      end if ! iiLH_Ni > 0
    end if ! iiLH_rice > 0

    ! Sampling for graupel (disabled)
    if ( iiLH_rgraupel > 0 ) then
      xp2_on_xm2_array_cloud(iiLH_rgraupel) = -999.

      call set_lower_triangular_matrix_sp &
           ( d_variables, iiLH_rgraupel, iiLH_s_mellor, -999., &
             corr_array_cloud )

      call set_lower_triangular_matrix_sp &
           ( d_variables, iiLH_rgraupel, iiLH_w, -999., &
             corr_array_cloud )

      if ( iiLH_Ngraupel > 0 ) then
        xp2_on_xm2_array_cloud(iiLH_Ngraupel) = -999.

        call set_lower_triangular_matrix_sp &
             ( d_variables, iiLH_rgraupel, iiLH_Ngraupel, -999., &
               corr_array_cloud )

        call set_lower_triangular_matrix_sp &
             ( d_variables, iiLH_Ngraupel, iiLH_s_mellor, -999., &
               corr_array_cloud )

        call set_lower_triangular_matrix_sp &
             ( d_variables, iiLH_Ngraupel, iiLH_w, -999., &
               corr_array_cloud )

      end if ! iiLH_Ngraupel > 0
    end if ! iiLH_rgraupel > 0

    if ( iiLH_Nc > 0 ) then
      ! The epsilon is a kluge to prevent a singular matrix in generate_lh_sample
      xp2_on_xm2_array_below(iiLH_Nc) = &
        max( Ncp2_on_Ncm2_below, epsilon( Ncp2_on_Ncm2_below ) )
      call set_lower_triangular_matrix_sp &
           ( d_variables, iiLH_Nc, iiLH_s_mellor, corr_sNc_NL_below, &
             corr_array_below )
    end if

    if ( iiLH_rrain > 0 ) then
      xp2_on_xm2_array_below(iiLH_rrain) = rrp2_on_rrainm2_below

      call set_lower_triangular_matrix_sp &
           ( d_variables, iiLH_rrain, iiLH_s_mellor, corr_srr_NL_below, &
             corr_array_below )

      if ( iiLH_Nr > 0 ) then
        xp2_on_xm2_array_below(iiLH_Nr) = Nrp2_on_Nrm2_below
        call set_lower_triangular_matrix_sp &
             ( d_variables, iiLH_rrain, iiLH_Nr, corr_rrNr_LL_below, &
               corr_array_below )

        call set_lower_triangular_matrix_sp &
             ( d_variables, iiLH_Nr, iiLH_s_mellor, corr_sNr_NL_below, &
               corr_array_below )

      end if ! iiLH_Nr > 0
    end if ! iiLH_rrain > 0

    if ( iiLH_rsnow > 0 ) then
      xp2_on_xm2_array_below(iiLH_rsnow) = rsnowp2_on_rsnowm2_below

      call set_lower_triangular_matrix_sp &
           ( d_variables, iiLH_rsnow, iiLH_s_mellor, corr_srsnow_NL_below, &
             corr_array_below )

      call set_lower_triangular_matrix_sp &
           ( d_variables, iiLH_rsnow, iiLH_w, corr_wrsnow_NL_below, &
             corr_array_below )

      if ( iiLH_Nsnow > 0 ) then
        xp2_on_xm2_array_below(iiLH_Nsnow) = Nsnowp2_on_Nsnowm2_below

        call set_lower_triangular_matrix_sp &
            ( d_variables, iiLH_rsnow, iiLH_Nsnow, corr_rsnowNsnow_LL_below, &
              corr_array_below )

        call set_lower_triangular_matrix_sp &
            ( d_variables, iiLH_Nsnow, iiLH_s_mellor, corr_sNsnow_NL_below, &
              corr_array_below )

        call set_lower_triangular_matrix_sp &
            ( d_variables, iiLH_Nsnow, iiLH_w, corr_wNsnow_NL_below, &
              corr_array_below )

      end if ! iiLH_Nsnow > 0
    end if ! iiLH_rsnow > 0

    if ( iiLH_rice > 0 ) then
      xp2_on_xm2_array_below(iiLH_rice) = ricep2_on_ricem2_below

      call set_lower_triangular_matrix_sp &
           ( d_variables, iiLH_rice, iiLH_s_mellor, corr_srice_NL_below, &
             corr_array_below )

      call set_lower_triangular_matrix_sp &
           ( d_variables, iiLH_rice, iiLH_w, corr_wrice_NL_below, &
             corr_array_below )

      if ( iiLH_Ni > 0 ) then
        xp2_on_xm2_array_below(iiLH_Ni) =  Nicep2_on_Nicem2_below

        call set_lower_triangular_matrix_sp &
            ( d_variables, iiLH_rice, iiLH_Ni, corr_riceNi_LL_below, &
              corr_array_below )

        call set_lower_triangular_matrix_sp &
             ( d_variables, iiLH_Ni, iiLH_s_mellor, corr_sNi_NL_below, &
               corr_array_below )

        call set_lower_triangular_matrix_sp &
             ( d_variables, iiLH_Ni, iiLH_w, corr_wNi_NL_below, &
               corr_array_below )

      end if ! iiLH_Ni > 0

    end if ! iiLH_rice > 0

    if ( iiLH_rgraupel > 0 ) then
      xp2_on_xm2_array_below(iiLH_rgraupel) = -999.

      call set_lower_triangular_matrix_sp &
           ( d_variables, iiLH_rgraupel, iiLH_s_mellor, -999., &
             corr_array_below )

!       call set_lower_triangular_matrix_sp &
!            ( d_variables, iiLH_rgraupel, iiLH_w, -999., &
!              corr_array_below )

      if ( iiLH_Ngraupel > 0 ) then
        xp2_on_xm2_array_below(iiLH_Ngraupel) = -999.

        call set_lower_triangular_matrix_sp &
            ( d_variables, iiLH_rgraupel, iiLH_Ngraupel, -999., &
              corr_array_below )

        call set_lower_triangular_matrix_sp &
             ( d_variables, iiLH_rgraupel, iiLH_s_mellor, -999., &
               corr_array_below )

!         call set_lower_triangular_matrix_sp &
!              ( d_variables, iiLH_rgraupel, iiLH_w, -999., &
!                corr_array_below )
      end if ! iiLH_Ngraupel > 0
    end if ! iiLH_rgraupel > 0

    ! Here we've hardwired the location and values for the ice variates.
    ! TODO: parameterize these properly.

    ! MPACE-A
!     corr_array_cloud(iiLH_Nsnow,iiLH_rice) = 0.93
!     corr_array_below(iiLH_Nsnow,iiLH_rice) = 0.93
!     corr_array_cloud(iiLH_Nsnow,iiLH_Ni) = 0.02
!     corr_array_below(iiLH_Nsnow,iiLH_Ni) = 0.02
!     corr_array_cloud(iiLH_rsnow,iiLH_Ni) = -0.10
!     corr_array_below(iiLH_rsnow,iiLH_Ni) = -0.10
!     corr_array_cloud(iiLH_rsnow,iiLH_rice) = 0.74
!     corr_array_below(iiLH_rsnow,iiLH_rice) = 0.74

!     corr_array_below(iiLH_Nc,iiLH_w) = 0.38
!     corr_array_below(iiLH_Nc,iiLH_s_mellor) = 0.67
!     corr_array_below(iiLH_rsnow,iiLH_Nc) = 0.30
!     corr_array_below(iiLH_Nsnow,iiLH_Nc) = 0.41
!     corr_array_below(iiLH_rice,iiLH_Nc) = 0.49
!     corr_array_below(iiLH_Ni,iiLH_Nc) = 0.08

!     corr_array_cloud(iiLH_Nc,iiLH_w) = 0.38
!     corr_array_cloud(iiLH_Nc,iiLH_s_mellor) = 0.67
!     corr_array_cloud(iiLH_rsnow,iiLH_Nc) = 0.30
!     corr_array_cloud(iiLH_Nsnow,iiLH_Nc) = 0.41
!     corr_array_cloud(iiLH_rice,iiLH_Nc) = 0.49
!     corr_array_cloud(iiLH_Ni,iiLH_Nc) = 0.08

    ! MPACE-B
!     corr_array_cloud(iiLH_Nsnow,iiLH_rice) = 0.87
!     corr_array_below(iiLH_Nsnow,iiLH_rice) = 0.87
!     corr_array_cloud(iiLH_Nsnow,iiLH_Ni) = 0.66
!     corr_array_below(iiLH_Nsnow,iiLH_Ni) = 0.66
!     corr_array_cloud(iiLH_rsnow,iiLH_Ni) = 0.50
!     corr_array_below(iiLH_rsnow,iiLH_Ni) = 0.50
!     corr_array_cloud(iiLH_rsnow,iiLH_rice) = 0.84
!     corr_array_below(iiLH_rsnow,iiLH_rice) = 0.84

!     corr_array_below(iiLH_Nc,iiLH_w) = 0.24
!     corr_array_below(iiLH_Nc,iiLH_s_mellor) = 0.94
!     corr_array_below(iiLH_rsnow,iiLH_Nc) = 0.44
!     corr_array_below(iiLH_Nsnow,iiLH_Nc) = 0.58
!     corr_array_below(iiLH_rice,iiLH_Nc) = 0.61
!     corr_array_below(iiLH_Ni,iiLH_Nc) = 0.89

!     corr_array_cloud(iiLH_Nc,iiLH_w) = 0.24
!     corr_array_cloud(iiLH_Nc,iiLH_s_mellor) = 0.94
!     corr_array_cloud(iiLH_rsnow,iiLH_Nc) = 0.44
!     corr_array_cloud(iiLH_Nsnow,iiLH_Nc) = 0.58
!     corr_array_cloud(iiLH_rice,iiLH_Nc) = 0.61
!     corr_array_cloud(iiLH_Ni,iiLH_Nc) = 0.89

    ! ISDAC Values
    if ( iiLH_Nsnow > 0 .and. iiLH_Ni > 0 .and. iiLH_Nc > 0 ) then
      corr_array_cloud(iiLH_Nsnow,iiLH_rice) = 0.49
      corr_array_below(iiLH_Nsnow,iiLH_rice) = 0.49
      corr_array_cloud(iiLH_Nsnow,iiLH_Ni) = 0.60
      corr_array_below(iiLH_Nsnow,iiLH_Ni) = 0.60
      corr_array_cloud(iiLH_rsnow,iiLH_Ni) = 0.43
      corr_array_below(iiLH_rsnow,iiLH_Ni) = 0.43
      corr_array_cloud(iiLH_rsnow,iiLH_rice) = 0.29
      corr_array_below(iiLH_rsnow,iiLH_rice) = 0.29

      corr_array_below(iiLH_Nc,iiLH_w) = 0.34
      corr_array_below(iiLH_Nc,iiLH_s_mellor) = 0.09
      corr_array_below(iiLH_rsnow,iiLH_Nc) = 0.14
      corr_array_below(iiLH_Nsnow,iiLH_Nc) = 0.21
      corr_array_below(iiLH_rice,iiLH_Nc) = 0.39
      corr_array_below(iiLH_Ni,iiLH_Nc) = 0.29

      corr_array_cloud(iiLH_Nc,iiLH_w) = 0.34
      corr_array_cloud(iiLH_Nc,iiLH_s_mellor) = 0.09
      corr_array_cloud(iiLH_rsnow,iiLH_Nc) = 0.14
      corr_array_cloud(iiLH_Nsnow,iiLH_Nc) = 0.21
      corr_array_cloud(iiLH_rice,iiLH_Nc) = 0.39
      corr_array_cloud(iiLH_Ni,iiLH_Nc) = 0.29
    end if

    ! Assume rr and Nr are uncorrelated with w for now.
    if ( iiLH_rrain > 0 .and. iiLH_Nr > 0 ) then
      corr_array_cloud(iiLH_rrain,iiLH_w) = 0.0
      corr_array_cloud(iiLH_Nr,iiLH_w) = 0.0
    end if

    ! Fill in values for t_mellor using s_mellor correlations
    do i = 3, d_variables
      corr_array_cloud(i,iiLH_t_mellor) = corr_array_cloud(i,iiLH_s_mellor) &
        * corr_array_cloud(iiLH_t_mellor,iiLH_s_mellor)
      corr_array_below(i,iiLH_t_mellor) = corr_array_below(i,iiLH_s_mellor) &
        * corr_array_below(iiLH_t_mellor,iiLH_s_mellor)
    end do

    return
  end subroutine setup_corr_varnce_array

  !-----------------------------------------------------------------------------
  subroutine cleanup_latin_hypercube_arrays( )

  ! Description:
  !   De-allocate latin hypercube arrays
  ! References:
  !   None
  !-----------------------------------------------------------------------------
    implicit none

    ! External
    intrinsic :: allocated

    ! ---- Begin Code ----

    if ( allocated( corr_array_cloud ) ) then
      deallocate( corr_array_cloud )
    end if

    if ( allocated( corr_array_below ) ) then
      deallocate( corr_array_below )
    end if

    if ( allocated( xp2_on_xm2_array_cloud ) ) then
      deallocate( xp2_on_xm2_array_cloud )
    end if

    if ( allocated( xp2_on_xm2_array_below ) ) then
      deallocate( xp2_on_xm2_array_below )
    end if

    return
  end subroutine cleanup_latin_hypercube_arrays

end module latin_hypercube_arrays
