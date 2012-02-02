!$Id$
!---------------------------------------------------------------------------------------------------
module KK_fixed_correlations

  use clubb_precision, only: &
    core_rknd

  implicit none

  public :: setup_KK_corr


  private ! Default scope

  ! Parameters for in-cloud (from SAM RF02 DO).
  real( kind = core_rknd ), public :: &       ! RF02 value
    corr_rrNr_LL_cloud,    & ! 0.786
    corr_srr_NL_cloud,     & ! 0.242
    corr_sNr_NL_cloud,     & ! 0.285
    corr_sNc_NL_cloud,     & ! 0.433
    corr_trr_NL_cloud,     & ! 0.260
    corr_tNr_NL_cloud,     & ! 0.204
    corr_tNc_NL_cloud        ! 0.165

!$omp threadprivate( corr_rrNr_LL_cloud, corr_srr_NL_cloud,  corr_sNr_NL_cloud, corr_sNc_NL_cloud, &
!$omp   corr_trr_NL_cloud, corr_tNr_NL_cloud, corr_tNc_NL_cloud )

  ! Parameters for below-cloud (from SAM RF02 DO).
  real( kind = core_rknd ), public :: &       ! RF02 value
    corr_rrNr_LL_below,    & ! 0.886
    corr_srr_NL_below,     & ! 0.056
    corr_sNr_NL_below,     & ! 0.015
    corr_sNc_NL_below,     & ! 0.00 ! Not applicable below cloud.
    corr_trr_NL_below,     & ! -0.066
    corr_tNr_NL_below,     & ! -0.094
    corr_tNc_NL_below        ! 0.00 ! Not applicable below cloud.

!$omp threadprivate( corr_rrNr_LL_below, corr_srr_NL_below, corr_sNr_NL_below, corr_sNc_NL_below, &
!$omp   corr_trr_NL_below, corr_tNr_NL_below, corr_tNc_NL_below )
  contains

!---------------------------------------------------------------------------------------------------
  subroutine setup_KK_corr( iunit, input_file_cloud, input_file_below, &
                            l_write_to_file, case_info_file )

! Description:
!   Set the correlations for variables needed in Khairoutdinov Kogan based on a
!   matrix.
!
! References:
!   None
!---------------------------------------------------------------------------------------------------
    use corr_matrix_module, only: &
      iiLH_s_mellor, & ! Variables
      iiLH_t_mellor, &
      iiLH_rrain, &
      iiLH_Nc, &
      iiLH_Nr

    use corr_matrix_module, only: &
      read_correlation_matrix ! Procedure(s)

    use error_code, only: &
      clubb_at_least_debug_level ! Procedure(s)

    use text_writer, only: &
      write_text ! Procedure(s)

    implicit none

    ! Constant parameters
    integer, parameter :: &
      d_variables = 5 ! Total variables in the correlation matrix for K&K

    ! Input Variables
    integer, intent(in) :: &
      iunit ! Fortran unit number for file I/O

    character(len=*), intent(in) :: &
      input_file_cloud, & ! Filename for correlations in cloud
      input_file_below    ! Filename for correlations out of cloud

    logical, intent(in) :: &
      l_write_to_file ! Write the case setup to a text file

    character(len=*), intent(in) :: &
      case_info_file ! File where we write to when l_write_to_file is true

    ! Local variables

    real(kind = core_rknd), dimension(d_variables,d_variables) :: &
      corr_matrix ! Matrix of all correlations
    
    ! ---- Begin Code ----

    ! Set the indices of the correlation array for the purpose of setting the
    ! variables (above)
    iiLH_s_mellor = 1
    iiLH_t_mellor = 2
    iiLH_rrain    = 3
    iiLH_Nc       = 4
    iiLH_Nr       = 5

    call read_correlation_matrix( iunit, input_file_cloud, d_variables, & ! In
                                  corr_matrix ) ! In/Out

    ! Note: Since corr_matrix contains only the lower elements of the matrix,
    ! the 2nd index must be the smaller of the 2 indices (set above)
    ! e.g. iiLH_s_mellor < iiLH_rrain.
    corr_srr_NL_cloud  = corr_matrix(iiLH_rrain,iiLH_s_mellor)
    corr_sNr_NL_cloud  = corr_matrix(iiLH_Nr,iiLH_s_mellor)
    corr_sNc_NL_cloud  = corr_matrix(iiLH_Nc,iiLH_s_mellor)

    corr_trr_NL_cloud  = corr_matrix(iiLH_rrain,iiLH_t_mellor)
    corr_tNr_NL_cloud  = corr_matrix(iiLH_Nr,iiLH_t_mellor)
    corr_tNc_NL_cloud  = corr_matrix(iiLH_Nc,iiLH_t_mellor)

    corr_rrNr_LL_cloud = corr_matrix(iiLH_Nr,iiLH_rrain)

    call read_correlation_matrix( iunit, input_file_below, d_variables, &
                                  corr_matrix )

    corr_srr_NL_below  = corr_matrix(iiLH_rrain,iiLH_s_mellor)
    corr_sNr_NL_below  = corr_matrix(iiLH_Nr,iiLH_s_mellor)
    corr_sNc_NL_below  = corr_matrix(iiLH_Nc,iiLH_s_mellor)

    corr_trr_NL_below  = corr_matrix(iiLH_rrain,iiLH_t_mellor)
    corr_tNr_NL_below  = corr_matrix(iiLH_Nr,iiLH_t_mellor)
    corr_tNc_NL_below  = corr_matrix(iiLH_Nc,iiLH_t_mellor)

    corr_rrNr_LL_below = corr_matrix(iiLH_Nr,iiLH_rrain)

    ! Set all indices back to -1, to avoid the introduction of bugs
    iiLH_s_mellor = -1; iiLH_t_mellor = -1; iiLH_rrain = -1; iiLH_Nc = -1; iiLH_Nr = -1

    ! Printing correlation values for debugging
    if ( clubb_at_least_debug_level( 1 ) ) then

      ! This will open the cases setup.txt file and append it to include the
      ! correlation values. This file was created and written to from clubb_driver previously.
      if ( l_write_to_file ) then 
        open(unit=iunit, file=case_info_file, status='old', action='write', position='append')
      end if

      call write_text( "---- Khairoutdinov Kogan correlations ----", l_write_to_file, iunit )

      ! Out of cloud
      call write_text( "corr_rrNr_LL_below = ", corr_rrNr_LL_below, &
                       l_write_to_file, iunit )
      call write_text( "corr_srr_NL_below = ", corr_srr_NL_below, &
                       l_write_to_file, iunit )
      call write_text( "corr_sNr_NL_below = ", corr_sNr_NL_below, &
                       l_write_to_file, iunit )
      call write_text( "corr_sNc_NL_below = ", corr_sNc_NL_below, &
                       l_write_to_file, iunit )
      call write_text( "corr_trr_NL_below = ", corr_trr_NL_below, &
                       l_write_to_file, iunit )
      call write_text( "corr_tNr_NL_below = ", corr_tNr_NL_below, &
                       l_write_to_file, iunit )
      call write_text( "corr_tNc_NL_below = ", corr_tNc_NL_below, &
                       l_write_to_file, iunit )

      ! In cloud
      call write_text( "corr_rrNr_LL_cloud = ", corr_rrNr_LL_cloud, &
                       l_write_to_file, iunit )
      call write_text( "corr_srr_NL_cloud = ", corr_srr_NL_cloud, &
                       l_write_to_file, iunit )
      call write_text( "corr_sNr_NL_cloud = ", corr_sNr_NL_cloud, &
                       l_write_to_file, iunit )
      call write_text( "corr_sNc_NL_cloud = ", corr_sNc_NL_cloud, &
                       l_write_to_file, iunit )
      call write_text( "corr_trr_NL_cloud = ", corr_trr_NL_cloud, &
                       l_write_to_file, iunit )
      call write_text( "corr_tNr_NL_cloud = ", corr_tNr_NL_cloud, &
                       l_write_to_file, iunit )
      call write_text( "corr_tNc_NL_cloud = ", corr_tNc_NL_cloud, &
                       l_write_to_file, iunit )

      ! Close the prior file
      if ( l_write_to_file ) close(unit=iunit)

    end if ! clubb_at_least_debug_level 1

    return
  end subroutine setup_KK_corr

end module KK_fixed_correlations
