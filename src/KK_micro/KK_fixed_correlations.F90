!$Id$
!===============================================================================
module KK_fixed_correlations

  use clubb_precision, only: &
      core_rknd

  implicit none

  private ! Default scope

  public :: setup_KK_corr

  ! Parameters for in-cloud (from SAM RF02 DO).
  real( kind = core_rknd ), public :: &  ! RF02 value
    corr_wrr_NL_cloud,  & ! Prescribed in-cloud correlation of w and rr     [-]
    corr_wNr_NL_cloud,  & ! Prescribed in-cloud correlation of w and Nr     [-]
    corr_wNcn_NL_cloud, & ! Prescribed in-cloud correlation of w and Ncn    [-]
    corr_sw_NN_cloud,   & ! Prescribed in-cloud correlation of s and w      [-]
    corr_srr_NL_cloud,  & ! Prescribed in-cloud correlation of s and rr     [-]
    corr_sNr_NL_cloud,  & ! Prescribed in-cloud correlation of s and Nr     [-]
    corr_sNcn_NL_cloud, & ! Prescribed in-cloud correlation of s and Ncn    [-]
    corr_trr_NL_cloud,  & ! Prescribed in-cloud correlation of t and rr     [-]
    corr_tNr_NL_cloud,  & ! Prescribed in-cloud correlation of t and Nr     [-]
    corr_tNcn_NL_cloud, & ! Prescribed in-cloud correlation of t and Ncn    [-]
    corr_rrNr_LL_cloud    ! Prescribed in-cloud correlation of rr and Nr    [-]

!$omp threadprivate( corr_wrr_NL_cloud, corr_wNr_NL_cloud, corr_wNcn_NL_cloud, &
!$omp                corr_srr_NL_cloud, corr_sNr_NL_cloud, corr_sNcn_NL_cloud, &
!$omp                corr_trr_NL_cloud, corr_tNr_NL_cloud, corr_tNcn_NL_cloud, &
!$omp                corr_rrNr_LL_cloud, corr_sw_NN_cloud )

  ! Parameters for below-cloud (from SAM RF02 DO).
  real( kind = core_rknd ), public :: &  ! RF02 value
    corr_wrr_NL_below,  & ! Prescribed below-cloud correlation of w and rr  [-]
    corr_wNr_NL_below,  & ! Prescribed below-cloud correlation of w and Nr  [-]
    corr_wNcn_NL_below, & ! Prescribed below-cloud correlation of w and Ncn [-]
    corr_sw_NN_below,   & ! Prescribed below-cloud correlation of s and w   [-]
    corr_srr_NL_below,  & ! Prescribed below-cloud correlation of s and rr  [-]
    corr_sNr_NL_below,  & ! Prescribed below-cloud correlation of s and Nr  [-]
    corr_sNcn_NL_below, & ! Prescribed below-cloud correlation of s and Ncn [-]
    corr_trr_NL_below,  & ! Prescribed below-cloud correlation of t and rr  [-]
    corr_tNr_NL_below,  & ! Prescribed below-cloud correlation of t and Nr  [-]
    corr_tNcn_NL_below, & ! Prescribed below-cloud correlation of t and Ncn [-]
    corr_rrNr_LL_below    ! Prescribed below-cloud correlation of rr and Nr [-]

!$omp threadprivate( corr_wrr_NL_below, corr_wNr_NL_below, corr_wNcn_NL_below, &
!$omp                corr_srr_NL_below, corr_sNr_NL_below, corr_sNcn_NL_below, &
!$omp                corr_trr_NL_below, corr_tNr_NL_below, corr_tNcn_NL_below, &
!$omp                corr_rrNr_LL_below, corr_sw_NN_below )

  ! Only needed when l_fix_s_t_correlations is true
  real( kind = core_rknd ), public :: &
    corr_st_NN_cloud, & ! Prescribed in-cloud correlation of s and t        [-]
    corr_st_NN_below    ! Prescribed below-cloud correlation of s and t     [-]

!$omp threadprivate( corr_st_NN_cloud, corr_st_NN_below )

  contains

  !=============================================================================
  subroutine setup_KK_corr( iunit, input_file_cloud, input_file_below, &
                            l_write_to_file, case_info_file )

    ! Description:
    ! Set the correlations for variables needed in Khairoutdinov Kogan based on
    ! a matrix.

    ! References:
    !   None
    !-----------------------------------------------------------------------

    use corr_matrix_module, only: &
        iiPDF_w,        & ! Variables
        iiPDF_s_mellor, &
        iiPDF_t_mellor, &
        iiPDF_rrain, &
        iiPDF_Nr, &
        iiPDF_Ncn

    use parameters_microphys, only: &
        l_fix_s_t_correlations ! Variable(s)

    use corr_matrix_module, only: &
        read_correlation_matrix ! Procedure(s)

    use error_code, only: &
        clubb_at_least_debug_level ! Procedure(s)

    use text_writer, only: &
        write_text ! Procedure(s)

    implicit none

    ! Constant Parameters
    integer, parameter :: &
      d_variables = 6   ! Total variables in the correlation matrix for K&K

    ! Input Variables
    integer, intent(in) :: &
      iunit   ! Fortran unit number for file I/O

    character(len=*), intent(in) :: &
      input_file_cloud, & ! Filename for correlations in cloud
      input_file_below    ! Filename for correlations out of cloud

    logical, intent(in) :: &
      l_write_to_file   ! Write the case setup to a text file

    character(len=*), intent(in) :: &
      case_info_file   ! File where we write to when l_write_to_file is true

    ! Local variables

    real(kind = core_rknd), dimension(d_variables,d_variables) :: &
      corr_matrix   ! Matrix of all correlations


    ! Set the indices of the correlation array for the purpose of setting the
    ! variables (above)
    iiPDF_s_mellor = 1
    iiPDF_t_mellor = 2
    iiPDF_w        = 3
    iiPDF_rrain    = 4
    iiPDF_Ncn      = 5
    iiPDF_Nr       = 6

    call read_correlation_matrix( iunit, input_file_cloud, d_variables, & ! In
                                  corr_matrix ) ! In/Out

    ! Note: Since corr_matrix contains only the lower elements of the matrix,
    ! the 2nd index must be the smaller of the 2 indices (set above)
    ! e.g. iiPDF_s_mellor < iiPDF_rrain.
    corr_wrr_NL_cloud  = corr_matrix(iiPDF_rrain,iiPDF_w)
    corr_wNr_NL_cloud  = corr_matrix(iiPDF_Nr,iiPDF_w)
    corr_wNcn_NL_cloud = corr_matrix(iiPDF_Ncn,iiPDF_w)

    corr_sw_NN_cloud   = corr_matrix(iiPDF_w,iiPDF_s_mellor)

    corr_srr_NL_cloud  = corr_matrix(iiPDF_rrain,iiPDF_s_mellor)
    corr_sNr_NL_cloud  = corr_matrix(iiPDF_Nr,iiPDF_s_mellor)
    corr_sNcn_NL_cloud = corr_matrix(iiPDF_Ncn,iiPDF_s_mellor)

    corr_trr_NL_cloud  = corr_matrix(iiPDF_rrain,iiPDF_t_mellor)
    corr_tNr_NL_cloud  = corr_matrix(iiPDF_Nr,iiPDF_t_mellor)
    corr_tNcn_NL_cloud = corr_matrix(iiPDF_Ncn,iiPDF_t_mellor)

    corr_rrNr_LL_cloud = corr_matrix(iiPDF_Nr,iiPDF_rrain)

    if ( l_fix_s_t_correlations ) then
       corr_st_NN_cloud = corr_matrix(iiPDF_t_mellor,iiPDF_s_mellor)
    else
        ! Don't fix the value of the correlation
       corr_st_NN_cloud = -999._core_rknd
    endif

    call read_correlation_matrix( iunit, input_file_below, d_variables, & ! In
                                  corr_matrix ) ! In/Out

    corr_wrr_NL_below  = corr_matrix(iiPDF_rrain,iiPDF_w)
    corr_wNr_NL_below  = corr_matrix(iiPDF_Nr,iiPDF_w)
    corr_wNcn_NL_below = corr_matrix(iiPDF_Ncn,iiPDF_w)

    corr_sw_NN_below   = corr_matrix(iiPDF_w,iiPDF_s_mellor)

    corr_srr_NL_below  = corr_matrix(iiPDF_rrain,iiPDF_s_mellor)
    corr_sNr_NL_below  = corr_matrix(iiPDF_Nr,iiPDF_s_mellor)
    corr_sNcn_NL_below = corr_matrix(iiPDF_Ncn,iiPDF_s_mellor)

    corr_trr_NL_below  = corr_matrix(iiPDF_rrain,iiPDF_t_mellor)
    corr_tNr_NL_below  = corr_matrix(iiPDF_Nr,iiPDF_t_mellor)
    corr_tNcn_NL_below = corr_matrix(iiPDF_Ncn,iiPDF_t_mellor)

    corr_rrNr_LL_below = corr_matrix(iiPDF_Nr,iiPDF_rrain)

    if ( l_fix_s_t_correlations ) then
       corr_st_NN_below = corr_matrix(iiPDF_t_mellor,iiPDF_s_mellor)
    else
       ! As above, we let this vary in space and time
       corr_st_NN_below = -999._core_rknd
    endif

    ! Set all indices back to -1, to avoid the introduction of bugs
    iiPDF_s_mellor = -1
    iiPDF_t_mellor = -1
    iiPDF_w = -1
    iiPDF_rrain = -1
    iiPDF_Ncn = -1
    iiPDF_Nr = -1

    ! Printing correlation values for debugging
    if ( clubb_at_least_debug_level( 1 ) ) then

      ! This will open the cases setup.txt file and append it to include the
      ! correlation values. This file was created and written to from
      ! clubb_driver previously.
      if ( l_write_to_file ) then 
         open( unit=iunit, file=case_info_file, status='old', action='write', &
               position='append' )
      endif

      call write_text( "---- Khairoutdinov Kogan correlations ----", &
                       l_write_to_file, iunit )

      ! Out of cloud
      call write_text( "corr_wrr_NL_below = ", corr_wrr_NL_below, &
                       l_write_to_file, iunit )
      call write_text( "corr_wNr_NL_below = ", corr_wNr_NL_below, &
                       l_write_to_file, iunit )
      call write_text( "corr_wNcn_NL_below = ", corr_wNcn_NL_below, &
                       l_write_to_file, iunit )
      call write_text( "corr_sw_NN_below = ", corr_sw_NN_below, &
                       l_write_to_file, iunit )
      call write_text( "corr_srr_NL_below = ", corr_srr_NL_below, &
                       l_write_to_file, iunit )
      call write_text( "corr_sNr_NL_below = ", corr_sNr_NL_below, &
                       l_write_to_file, iunit )
      call write_text( "corr_sNcn_NL_below = ", corr_sNcn_NL_below, &
                       l_write_to_file, iunit )
      call write_text( "corr_trr_NL_below = ", corr_trr_NL_below, &
                       l_write_to_file, iunit )
      call write_text( "corr_tNr_NL_below = ", corr_tNr_NL_below, &
                       l_write_to_file, iunit )
      call write_text( "corr_tNcn_NL_below = ", corr_tNcn_NL_below, &
                       l_write_to_file, iunit )
      call write_text( "corr_rrNr_LL_below = ", corr_rrNr_LL_below, &
                       l_write_to_file, iunit )
      if ( l_fix_s_t_correlations ) then
         call write_text( "corr_st_NN_below = ", corr_st_NN_below, &
                          l_write_to_file, iunit )
      endif

      ! In cloud
      call write_text( "corr_wrr_NL_cloud = ", corr_wrr_NL_cloud, &
                       l_write_to_file, iunit )
      call write_text( "corr_wNr_NL_cloud = ", corr_wNr_NL_cloud, &
                       l_write_to_file, iunit )
      call write_text( "corr_wNcn_NL_cloud = ", corr_wNcn_NL_cloud, &
                       l_write_to_file, iunit )
      call write_text( "corr_sw_NN_cloud = ", corr_sw_NN_below, &
                       l_write_to_file, iunit )
      call write_text( "corr_srr_NL_cloud = ", corr_srr_NL_cloud, &
                       l_write_to_file, iunit )
      call write_text( "corr_sNr_NL_cloud = ", corr_sNr_NL_cloud, &
                       l_write_to_file, iunit )
      call write_text( "corr_sNcn_NL_cloud = ", corr_sNcn_NL_cloud, &
                       l_write_to_file, iunit )
      call write_text( "corr_trr_NL_cloud = ", corr_trr_NL_cloud, &
                       l_write_to_file, iunit )
      call write_text( "corr_tNr_NL_cloud = ", corr_tNr_NL_cloud, &
                       l_write_to_file, iunit )
      call write_text( "corr_tNcn_NL_cloud = ", corr_tNcn_NL_cloud, &
                       l_write_to_file, iunit )
      call write_text( "corr_rrNr_LL_cloud = ", corr_rrNr_LL_cloud, &
                       l_write_to_file, iunit )
      if ( l_fix_s_t_correlations ) then
         call write_text( "corr_st_NN_cloud = ", corr_st_NN_cloud, &
                          l_write_to_file, iunit )
      endif

      ! Close the prior file
      if ( l_write_to_file ) close(unit=iunit)

    endif ! clubb_at_least_debug_level 1


    return

  end subroutine setup_KK_corr

!===============================================================================

end module KK_fixed_correlations
