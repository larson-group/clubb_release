!$Id$
!===============================================================================
module KK_fixed_correlations

  implicit none

  private ! default scope

  public :: setup_KK_corr

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

    use fixed_correlations, only: &
        corr_wrr_NL_cloud,  & ! Variable(s)
        corr_wNr_NL_cloud,  &
        corr_wNcn_NL_cloud, &
        corr_sw_NN_cloud,   &
        corr_srr_NL_cloud,  &
        corr_sNr_NL_cloud,  &
        corr_sNcn_NL_cloud, &
        corr_trr_NL_cloud,  &
        corr_tNr_NL_cloud,  &
        corr_tNcn_NL_cloud, &
        corr_rrNr_LL_cloud, &
        corr_wrr_NL_below,  &
        corr_wNr_NL_below,  &
        corr_wNcn_NL_below, &
        corr_sw_NN_below,   &
        corr_srr_NL_below,  &
        corr_sNr_NL_below,  &
        corr_sNcn_NL_below, &
        corr_trr_NL_below,  &
        corr_tNr_NL_below,  &
        corr_tNcn_NL_below, &
        corr_rrNr_LL_below, &
        corr_st_NN_cloud,   &
        corr_st_NN_below

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

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    use error_code, only: &
        clubb_at_least_debug_level ! Procedure(s)

    use text_writer, only: &
        write_text ! Procedure(s)

    use matrix_operations, only: &
        get_lower_triangular_matrix ! Procedure(s)

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
    ! The subroutine get_lower_triangular_matrix handles this automatically.

    call get_lower_triangular_matrix( d_variables, iiPDF_rrain, iiPDF_w, corr_matrix, & ! In
                                 corr_wrr_NL_cloud ) ! Out

    call get_lower_triangular_matrix( d_variables, iiPDF_Nr, iiPDF_w, corr_matrix, & ! In
                                 corr_wNr_NL_cloud ) ! Out

    call get_lower_triangular_matrix( d_variables, iiPDF_Ncn, iiPDF_w, corr_matrix, & ! In
                                 corr_wNcn_NL_cloud ) ! Out

    call get_lower_triangular_matrix( d_variables, iiPDF_w, iiPDF_s_mellor, corr_matrix, & ! In
                                 corr_sw_NN_cloud ) ! Out

    call get_lower_triangular_matrix( d_variables, iiPDF_rrain, iiPDF_s_mellor, corr_matrix, & ! In
                                 corr_srr_NL_cloud ) ! Out

    call get_lower_triangular_matrix( d_variables, iiPDF_Nr, iiPDF_s_mellor, corr_matrix, & ! In
                                 corr_sNr_NL_cloud ) ! Out

    call get_lower_triangular_matrix( d_variables, iiPDF_Ncn, iiPDF_s_mellor, corr_matrix, & ! In
                                 corr_sNcn_NL_cloud ) ! Out

    call get_lower_triangular_matrix( d_variables, iiPDF_rrain, iiPDF_t_mellor, corr_matrix, & ! In
                                 corr_trr_NL_cloud ) ! Out

    call get_lower_triangular_matrix( d_variables, iiPDF_Nr, iiPDF_t_mellor, corr_matrix, & ! In
                                 corr_tNr_NL_cloud ) ! Out

    call get_lower_triangular_matrix( d_variables, iiPDF_Ncn, iiPDF_t_mellor, corr_matrix, & ! In
                                 corr_tNcn_NL_cloud ) ! Out


    call get_lower_triangular_matrix( d_variables, iiPDF_Nr, iiPDF_rrain, corr_matrix, & ! In
                                 corr_rrNr_LL_cloud ) ! Out

    if ( l_fix_s_t_correlations ) then
       call get_lower_triangular_matrix( d_variables, iiPDF_t_mellor, iiPDF_s_mellor, corr_matrix, &
                                    corr_st_NN_cloud ) ! Out
    else
        ! Don't fix the value of the correlation
       corr_st_NN_cloud = -999._core_rknd

    end if ! l_fix_s_t_correlations

    ! Below cloud
    call read_correlation_matrix( iunit, input_file_below, d_variables, & ! In
                                  corr_matrix ) ! In/Out


    call get_lower_triangular_matrix( d_variables, iiPDF_rrain, iiPDF_w, corr_matrix, & ! In
                                 corr_wrr_NL_below ) ! Out

    call get_lower_triangular_matrix( d_variables, iiPDF_Nr, iiPDF_w, corr_matrix, & ! In
                                 corr_wNr_NL_below ) ! Out

    call get_lower_triangular_matrix( d_variables, iiPDF_Ncn, iiPDF_w, corr_matrix, & ! In
                                 corr_wNcn_NL_below ) ! Out

    call get_lower_triangular_matrix( d_variables, iiPDF_w, iiPDF_s_mellor, corr_matrix, & ! In
                                 corr_sw_NN_below ) ! Out


    call get_lower_triangular_matrix( d_variables, iiPDF_rrain, iiPDF_s_mellor, corr_matrix, & ! In
                                 corr_srr_NL_below ) ! Out

    call get_lower_triangular_matrix( d_variables, iiPDF_Nr, iiPDF_s_mellor, corr_matrix, & ! In
                                 corr_sNr_NL_below ) ! Out

    call get_lower_triangular_matrix( d_variables, iiPDF_Ncn, iiPDF_s_mellor, corr_matrix, & ! In
                                 corr_sNcn_NL_below ) ! Out

    call get_lower_triangular_matrix( d_variables, iiPDF_rrain, iiPDF_t_mellor, corr_matrix, & ! In
                                 corr_trr_NL_below ) ! Out

    call get_lower_triangular_matrix( d_variables, iiPDF_Nr, iiPDF_t_mellor, corr_matrix, & ! In
                                 corr_tNr_NL_below ) ! Out

    call get_lower_triangular_matrix( d_variables, iiPDF_Ncn, iiPDF_t_mellor, corr_matrix, & ! In
                                 corr_tNcn_NL_below ) ! Out


    call get_lower_triangular_matrix( d_variables, iiPDF_Nr, iiPDF_rrain, corr_matrix, & ! In
                                 corr_rrNr_LL_below ) ! Out

    if ( l_fix_s_t_correlations ) then
       call get_lower_triangular_matrix( d_variables, iiPDF_t_mellor, iiPDF_s_mellor, corr_matrix, &
                                    corr_st_NN_below ) ! Out
    else
        ! Don't fix the value of the correlation
       corr_st_NN_below = -999._core_rknd

    end if ! l_fix_s_t_correlations

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
