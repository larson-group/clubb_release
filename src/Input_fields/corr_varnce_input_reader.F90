!-----------------------------------------------------------------------
! $Id$
!-------------------------------------------------------------------------------
module corr_varnce_input_reader

  use clubb_precision, only: &
      core_rknd

  implicit none

  private

  public :: read_correlation_matrix

contains

  !-----------------------------------------------------------------------------
  subroutine read_correlation_matrix( iunit, input_file, &
                                      pdf_dim, hm_metadata, &
                                      corr_array_n )

    ! Description:
    !   Reads a correlation variance array from a file and stores it in an array.
    !-----------------------------------------------------------------------------

    use input_reader, only: &
        one_dim_read_var, & ! Variable(s)
        read_one_dim_file, deallocate_one_dim_vars, count_columns ! Procedure(s)

    use matrix_operations, only: set_lower_triangular_matrix ! Procedure(s)

    use constants_clubb, only: fstderr ! Variable(s)

    use corr_varnce_module, only: &
        hm_metadata_type, &
        get_corr_var_index

    implicit none

    ! Input Variable(s)
    integer, intent(in) :: &
      iunit, &   ! File I/O unit
      pdf_dim    ! number of variables in the array

    type (hm_metadata_type), intent(in) :: &
      hm_metadata

    character(len=*), intent(in) :: input_file ! Path to the file

    ! Input/Output Variable(s)
    real( kind = core_rknd ), dimension(pdf_dim,pdf_dim), intent(inout) :: &
      corr_array_n ! Normal space correlation array

    ! Local Variable(s)
    type(one_dim_read_var), allocatable, dimension(:) :: &
      retVars ! stores the variables read in from the corr_varnce.in file

    integer :: &
      var_index1, & ! variable index
      var_index2, & ! variable index
      nCols, &      ! the number of columns in the file
      i, j          ! Loop index

    !--------------------------- BEGIN CODE -------------------------

    nCols = count_columns( iunit, input_file )

    ! Allocate all arrays based on pdf_dim
    allocate( retVars(1:nCols) )

    ! Initializing to zero means that correlations we don't have are assumed to be 0.
    corr_array_n(:,:) = 0.0_core_rknd

    ! Set main diagonal to 1
    do i=1, pdf_dim
      corr_array_n(i,i) = 1.0_core_rknd
    end do

    ! Read the values from the specified file
    call read_one_dim_file( iunit, nCols, input_file, & ! intent(in)
                            retVars ) ! intent(out)

    if( size( retVars(1)%values ) /= nCols ) then
      write(fstderr, *) "Correlation matrix must have an equal number of rows and cols in file ", &
            input_file
      error stop "Bad data in correlation file."
    end if

    ! Start at 2 because the first index is always just 1.0 in the first row
    ! and the rest of the rows are ignored
    do i=2, nCols
      var_index1 = get_corr_var_index( retVars(i)%name, hm_metadata )
      if( var_index1 > -1 ) then
        do j=1, (i-1)
          var_index2 = get_corr_var_index( retVars(j)%name, hm_metadata )
          if( var_index2 > -1 ) then
            call set_lower_triangular_matrix &
                 ( pdf_dim, var_index1, var_index2, retVars(i)%values(j), & ! intent(in)
                   corr_array_n ) ! intent(inout)
          end if
        end do
      end if
    end do

    call deallocate_one_dim_vars( nCols, & ! intent(in)
                                  retVars ) ! intent(inout)

    return
  end subroutine read_correlation_matrix

end module corr_varnce_input_reader
