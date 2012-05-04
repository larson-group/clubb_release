! $Id$
module latin_hypercube_arrays

  use clubb_precision, only: &
    dp, & ! double precision
    core_rknd

  implicit none

  public :: setup_corr_varnce_array, cleanup_latin_hypercube_arrays

  private

  integer, public :: d_variables
!omp threadprivate(d_variables)

  real( kind = core_rknd ), public, dimension(:), allocatable :: &
    xp2_on_xm2_array_cloud, &
    xp2_on_xm2_array_below

  real( kind = core_rknd ), public, dimension(:,:), allocatable :: &
    corr_array_cloud, &
    corr_array_below
!$omp threadprivate(xp2_on_xm2_array_cloud, xp2_on_xm2_array_below, &
!$omp   corr_array_cloud, corr_array_below)

  logical, public :: &
    l_fixed_corr_initialized = .false.

!$omp threadprivate(l_fixed_corr_initialized)

  real( kind = dp ), allocatable, dimension(:,:), target, public :: &
    corr_stw_cloud_Cholesky, & ! Cholesky factorization of the correlation matrix
    corr_stw_below_Cholesky    ! Cholesky factorization of the correlation matrix

!$omp threadprivate(corr_stw_cloud_Cholesky, corr_stw_below_Cholesky)

  real( kind = dp ), allocatable, dimension(:), public :: &
    corr_stw_cloud_scaling, & ! Scaling factors for the correlation matrix [-]
    corr_stw_below_scaling    ! Scaling factors for the correlation matrix [-]

!$omp threadprivate(corr_stw_cloud_scaling, corr_stw_below_scaling)

  logical, public :: &
    l_corr_stw_cloud_scaling, & ! Whether we're scaling the correlation matrix
    l_corr_stw_below_scaling

!$omp threadprivate(l_corr_stw_cloud_scaling, l_corr_stw_below_scaling)

  integer, allocatable, dimension(:,:,:), public :: & 
    height_time_matrix ! matrix of rand ints

!$omp threadprivate(height_time_matrix)

  contains
!===============================================================================
  subroutine setup_corr_varnce_array( iiNcm, iirrainm, iiNrm, iiricem, iiNim, iirsnowm, iiNsnowm, &
                                      l_ice_micro, &
                                      input_file_cloud, input_file_below, iunit )

! Description:
!   Setup an array with the x'^2/xm^2 variables on the diagonal and the other
!   elements to be correlations between various variables.

! References:
!   None.
!-------------------------------------------------------------------------------

    use parameters_microphys, only: &
      rrp2_on_rrainm2_cloud, Nrp2_on_Nrm2_cloud, Ncp2_on_Ncm2_cloud, & ! Variables
      rrp2_on_rrainm2_below, Nrp2_on_Nrm2_below, Ncp2_on_Ncm2_below

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
      l_fix_s_t_correlations ! Variable(s)

!   use matrix_operations, only: print_lower_triangular_matrix ! Procedure(s)

    use constants_clubb, only: &
      fstdout, & ! Constant(s)
      fstderr, &
      zero

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    use corr_matrix_module, only: &
      iiLH_s_mellor, & ! Variable(s)
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

    use corr_matrix_module, only: &
      read_correlation_matrix ! Procedure(s)

    implicit none

    ! External
    intrinsic :: max, epsilon, trim

    ! Input Variables
    integer, intent(in) :: &
      iiNcm,    & ! Index of cloud droplet number conc.
      iirrainm, & ! Index of rain water mixing ratio
      iiNrm,    & ! Index of rain droplet number conc.
      iiricem,  & ! Index of ice water mixing ratio
      iiNim,    & ! Index of ice crystal number conc.
      iirsnowm, & ! Index snow mixing ratio
      iiNsnowm    ! Index of snow number conc.

    logical, intent(in) :: &
      l_ice_micro  ! Whether the microphysics scheme will do ice

    integer, intent(in) :: &
      iunit ! The file unit

    character(len=*), intent(in) :: &
      input_file_cloud, & ! Path to the in cloud correlation file
      input_file_below    ! Path to the out of cloud correlation file

    ! Local variables
    character(len=1) :: response
    logical :: l_warning
    integer :: i

    ! ---- Begin Code ----

    iiLH_s_mellor = 1 ! Extended rcm
    iiLH_t_mellor = 2 ! 't' orthogonal to 's'
    iiLH_w        = 3 ! vertical velocity
    iiLH_Nc       = 4 ! Cloud droplet number concentration

    i = iiLH_Nc

    call return_LH_index( iirrainm, i, iiLH_rrain )
    call return_LH_index( iiNrm, i, iiLH_Nr )
    if ( l_ice_micro ) then
      call return_LH_index( iiricem, i, iiLH_rice )
      call return_LH_index( iiNim, i, iiLH_Ni )
      call return_LH_index( iirsnowm, i, iiLH_rsnow )
      call return_LH_index( iiNsnowm, i, iiLH_Nsnow )
    else
      iiLH_rice = -1
      iiLH_Ni = -1
      iiLH_rsnow = -1
      iiLH_Nsnow = -1
    end if
    ! Disabled until we have values for the correlations of graupel and
    ! other variates in the latin hypercube sampling.
    iiLH_rgraupel = -1
    iiLH_Ngraupel = -1

    d_variables = i

    allocate( corr_array_cloud(d_variables,d_variables) )
    allocate( corr_array_below(d_variables,d_variables) )

    allocate( xp2_on_xm2_array_cloud(d_variables) )
    allocate( xp2_on_xm2_array_below(d_variables) )

    xp2_on_xm2_array_cloud(:) = zero
    xp2_on_xm2_array_below(:) = zero

    call read_correlation_matrix( iunit, trim( input_file_cloud ), d_variables, & ! In
                                  corr_array_cloud ) ! Out

    call read_correlation_matrix( iunit, trim( input_file_below ), d_variables, & ! In
                                  corr_array_below ) ! Out

    ! Sanity check to avoid confusing non-convergence results.
    if ( .not. l_fix_s_t_correlations .and. iiLH_Nc > 0 ) then
      l_warning = .false.
      do i = 1, d_variables
        if ( ( corr_array_cloud(i,iiLH_Nc) /= zero .or.  &
               corr_array_below(i,iiLH_Nc) /= zero ) .and. &
             i /= iiLH_Nc ) then
          l_warning = .true.
          print *, i, corr_array_cloud(i,iiLH_Nc) 
        end if
      end do ! 1..d_variables
      if ( l_warning ) then
        write(fstderr,*) "Warning: the specified correlations for s and Nc are non-zero."
        write(fstderr,*) "The latin hypercube code will not converge to the analytic solution "// &
          "using these settings."
        write(fstderr,'(A)',advance='no') "Continue? "
        read(*,*) response
        if ( response(1:1) /= 'y' .and. response(1:1) /= 'Y' ) then
           stop "Exiting..."
        end if
      end if
    end if ! l_fix_s_t_correlations

    if ( iiLH_Nc > 0 ) then
      xp2_on_xm2_array_cloud(iiLH_Nc) = Ncp2_on_Ncm2_cloud
    end if

    if ( iiLH_rrain > 0 ) then
      xp2_on_xm2_array_cloud(iiLH_rrain) = rrp2_on_rrainm2_cloud
      if ( iiLH_Nr > 0 ) then
        xp2_on_xm2_array_cloud(iiLH_Nr) = Nrp2_on_Nrm2_cloud
      end if ! iiLH_Nr > 0
    end if ! iiLH_rrain > 0

    if ( iiLH_rsnow > 0 ) then
      xp2_on_xm2_array_cloud(iiLH_rsnow) = rsnowp2_on_rsnowm2_cloud


      if ( iiLH_Nsnow > 0 ) then
        xp2_on_xm2_array_cloud(iiLH_Nsnow) = Nsnowp2_on_Nsnowm2_cloud


      end if ! iiLH_Nsnow > 0
    end if ! iiLH_rsnow > 0

    if ( iiLH_rice > 0 ) then
      xp2_on_xm2_array_cloud(iiLH_rice) = ricep2_on_ricem2_cloud


      if ( iiLH_Ni > 0 ) then
        xp2_on_xm2_array_cloud(iiLH_Ni) = Nicep2_on_Nicem2_cloud

      end if ! iiLH_Ni > 0
    end if ! iiLH_rice > 0

    ! Sampling for graupel (disabled)
    if ( iiLH_rgraupel > 0 ) then
      xp2_on_xm2_array_cloud(iiLH_rgraupel) = -999._core_rknd


      if ( iiLH_Ngraupel > 0 ) then
        xp2_on_xm2_array_cloud(iiLH_Ngraupel) = -999._core_rknd


      end if ! iiLH_Ngraupel > 0
    end if ! iiLH_rgraupel > 0

    if ( iiLH_Nc > 0 ) then
      ! The epsilon is a kluge to prevent a singular matrix in generate_lh_sample
      xp2_on_xm2_array_below(iiLH_Nc) = &
        max( Ncp2_on_Ncm2_below, epsilon( Ncp2_on_Ncm2_below ) )

    end if

    if ( iiLH_rrain > 0 ) then
      xp2_on_xm2_array_below(iiLH_rrain) = rrp2_on_rrainm2_below



      if ( iiLH_Nr > 0 ) then
        xp2_on_xm2_array_below(iiLH_Nr) = Nrp2_on_Nrm2_below


      end if ! iiLH_Nr > 0
    end if ! iiLH_rrain > 0

    if ( iiLH_rsnow > 0 ) then
      xp2_on_xm2_array_below(iiLH_rsnow) = rsnowp2_on_rsnowm2_below


      if ( iiLH_Nsnow > 0 ) then
        xp2_on_xm2_array_below(iiLH_Nsnow) = Nsnowp2_on_Nsnowm2_below

      end if ! iiLH_Nsnow > 0
    end if ! iiLH_rsnow > 0

    if ( iiLH_rice > 0 ) then
      xp2_on_xm2_array_below(iiLH_rice) = ricep2_on_ricem2_below


      if ( iiLH_Ni > 0 ) then
        xp2_on_xm2_array_below(iiLH_Ni) =  Nicep2_on_Nicem2_below
      end if ! iiLH_Ni > 0

    end if ! iiLH_rice > 0

    if ( iiLH_rgraupel > 0 ) then
      xp2_on_xm2_array_below(iiLH_rgraupel) = -999._core_rknd


      if ( iiLH_Ngraupel > 0 ) then
        xp2_on_xm2_array_below(iiLH_Ngraupel) = -999._core_rknd


      end if ! iiLH_Ngraupel > 0
    end if ! iiLH_rgraupel > 0

    return
  end subroutine setup_corr_varnce_array

  !-----------------------------------------------------------------------------
  subroutine cleanup_latin_hypercube_arrays( )

    ! Description:
    !   De-allocate latin hypercube arrays
    ! References:
    !   None
    !---------------------------------------------------------------------------
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

    if ( allocated( corr_stw_cloud_Cholesky ) ) then
      deallocate( corr_stw_cloud_Cholesky )
    end if

    if ( allocated( corr_stw_below_Cholesky ) ) then
      deallocate( corr_stw_below_Cholesky )
    end if

    if ( allocated( corr_stw_cloud_scaling ) ) then
      deallocate( corr_stw_cloud_scaling )
    end if

    if ( allocated( corr_stw_below_scaling ) ) then
      deallocate( corr_stw_below_scaling )
    end if

    if ( allocated( height_time_matrix ) ) then
      deallocate( height_time_matrix )
    end if

    return
  end subroutine cleanup_latin_hypercube_arrays

!-------------------------------------------------------------------------------
  subroutine return_LH_index( hydromet_index, LH_count, LH_index )

    ! Description:
    !   Set the Latin hypercube variable index if the hydrometeor exists
    ! References:
    !   None
    !-------------------------------------------------------------------------

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      hydromet_index

    ! Input/Output Variables
    integer, intent(inout) :: &
      LH_count

    ! Output Variables
    integer, intent(out) :: &
      LH_index

    ! ---- Begin Code ----

    if ( hydromet_index > 0 ) then
      LH_count = LH_count + 1
      LH_index = LH_count
    else
      LH_index = -1
    end if

    return
  end subroutine return_LH_index

end module latin_hypercube_arrays
