! $Id$
module latin_hypercube_arrays

  implicit none

  public :: setup_corr_varnce_array

  private

  real, public, dimension(:), allocatable :: &
    xp2_on_xm2_array_cloud, &
    xp2_on_xm2_array_below

  real, public, dimension(:,:), allocatable :: &
    corr_array_cloud, &
    corr_array_below

  contains
!===============================================================================
  subroutine setup_corr_varnce_array( d_variables )
! Description:
!   Setup an array with the x'^2/x variables on the diagonal and the other
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

    use array_index, only: &
      iiLH_s_mellor, & ! Variables
      iiLH_t_mellor, & 
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
      d_variables ! Number of variates in the array

    integer :: i

    ! ---- Begin Code ----

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
             ( d_variables, iiLH_s_mellor, iiLH_Nr, corr_srr_NL_cloud, &
               corr_array_cloud )
      end if ! iiLH_Nr > 0
    end if ! iiLH_rrain > 0

      ! Placeholder until we have actual numbers for the correlations/variances
      ! of ice phase hydrometeors.
      if ( iiLH_rsnow > 0 ) then
        xp2_on_xm2_array_cloud(iiLH_rsnow) = rrp2_on_rrainm2_cloud
        call set_lower_triangular_matrix_sp &
             ( d_variables, iiLH_rsnow, iiLH_s_mellor, corr_srr_NL_cloud, &
               corr_array_cloud )
        if ( iiLH_Nsnow > 0 ) then
          xp2_on_xm2_array_cloud(iiLH_Nsnow) = Nrp2_on_Nrm2_cloud
          call set_lower_triangular_matrix_sp &
               ( d_variables, iiLH_rsnow, iiLH_Nsnow, corr_rrNr_LL_cloud, &
                 corr_array_cloud )
          call set_lower_triangular_matrix_sp &
               ( d_variables, iiLH_Nsnow, iiLH_s_mellor, corr_sNr_NL_cloud, &
                 corr_array_cloud )
        end if ! iiLH_Nsnow > 0
      end if ! iiLH_rsnow > 0
      if ( iiLH_rice > 0 ) then
        xp2_on_xm2_array_cloud(iiLH_rice) = 1.0 ! Dimensionless made up value
        call set_lower_triangular_matrix_sp &
             ( d_variables, iiLH_rice, iiLH_s_mellor, corr_sNc_NL_cloud, &
               corr_array_cloud )
        if ( iiLH_Ni > 0 ) then
          xp2_on_xm2_array_cloud(iiLH_Ni) = Ncp2_on_Ncm2_cloud
          call set_lower_triangular_matrix_sp &
               ( d_variables, iiLH_Ni, iiLH_s_mellor, corr_sNc_NL_cloud, &
                 corr_array_cloud )
        end if ! iiLH_Ni > 0
      end if ! iiLH_rice > 0
      if ( iiLH_rgraupel > 0 ) then
        xp2_on_xm2_array_cloud(iiLH_rgraupel) = rrp2_on_rrainm2_cloud
        call set_lower_triangular_matrix_sp &
             ( d_variables, iiLH_rgraupel, iiLH_s_mellor, corr_srr_NL_cloud, &
               corr_array_cloud )
        if ( iiLH_Ngraupel > 0 ) then
          xp2_on_xm2_array_cloud(iiLH_Ngraupel) = Nrp2_on_Nrm2_cloud
          call set_lower_triangular_matrix_sp &
               ( d_variables, iiLH_rgraupel, iiLH_Ngraupel, corr_rrNr_LL_cloud, &
                 corr_array_cloud )
          call set_lower_triangular_matrix_sp &
               ( d_variables, iiLH_Ngraupel, iiLH_s_mellor, corr_sNr_NL_cloud, &
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

      ! Placeholder until we have actual numbers for the correlations/variances
      ! of ice phase hydrometeors.
      if ( iiLH_rsnow > 0 ) then
        xp2_on_xm2_array_below(iiLH_rsnow) = rrp2_on_rrainm2_below
        call set_lower_triangular_matrix_sp &
             ( d_variables, iiLH_rsnow, iiLH_s_mellor, corr_srr_NL_below, &
               corr_array_below )
        if ( iiLH_Nsnow > 0 ) then
          xp2_on_xm2_array_below(iiLH_Nsnow)    = Nrp2_on_Nrm2_below
          call set_lower_triangular_matrix_sp &
              ( d_variables, iiLH_rsnow, iiLH_Nsnow, corr_rrNr_LL_below, &
                corr_array_below )
          call set_lower_triangular_matrix_sp &
              ( d_variables, iiLH_Nsnow, iiLH_s_mellor, corr_sNr_NL_below, &
                corr_array_below )
        end if ! iiLH_Nsnow > 0
      end if ! iiLH_rsnow > 0
      if ( iiLH_rice > 0 ) then
        xp2_on_xm2_array_below(iiLH_rice) = 2.0 ! Dimensionless made up value
        call set_lower_triangular_matrix_sp &
             ( d_variables, iiLH_rice, iiLH_s_mellor, corr_sNc_NL_below, &
               corr_array_below )
        if ( iiLH_Ni > 0 ) then
          xp2_on_xm2_array_below(iiLH_Ni)     = Ncp2_on_Ncm2_below
          call set_lower_triangular_matrix_sp &
               ( d_variables, iiLH_Ni, iiLH_s_mellor, corr_sNc_NL_below, &
                 corr_array_below )
        end if ! iiLH_Ni > 0
      end if ! iiLH_rice > 0
      if ( iiLH_rgraupel > 0 ) then
        xp2_on_xm2_array_below(iiLH_rgraupel) = rrp2_on_rrainm2_below
        call set_lower_triangular_matrix_sp &
             ( d_variables, iiLH_rgraupel, iiLH_s_mellor, corr_srr_NL_below, &
               corr_array_below )
        if ( iiLH_Ngraupel > 0 ) then
          xp2_on_xm2_array_below(iiLH_Ngraupel) = Nrp2_on_Nrm2_below
          call set_lower_triangular_matrix_sp &
              ( d_variables, iiLH_rgraupel, iiLH_Ngraupel, corr_rrNr_LL_below, &
                corr_array_below )
          call set_lower_triangular_matrix_sp &
               ( d_variables, iiLH_rgraupel, iiLH_s_mellor, corr_sNr_NL_below, &
                 corr_array_below )
        end if ! iiLH_Ngraupel > 0
      end if ! iiLH_rgraupel > 0

      ! Approximate the correlation between t and other variates.
      ! Since the indices are always s < t < all other variates we can use
      ! this iterated multiplication rather than using set/get_lower_triangular_matrix
!     do i = 4, d_variables
!       corr_array_below(i,iiLH_t_mellor) = &
!         corr_array_below(iiLH_t_mellor,iiLH_s_mellor) * corr_array_below(i,iiLH_s_mellor)
!       corr_array_cloud(i,iiLH_t_mellor) = &
!         corr_array_cloud(iiLH_t_mellor,iiLH_s_mellor) * corr_array_cloud(i,iiLH_s_mellor)
!     end do

      return
    end subroutine setup_corr_varnce_array

  end module latin_hypercube_arrays
