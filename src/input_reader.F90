!$Id$
module input_reader
!
!  This module is respondsible for the procedures and structures necessary to
!  read in "SAM-Like" case specific files. Currently only the
!  <casename>_sounding.in file is formatted to be used by this module.
!
!---------------------------------------------------------------------------------------------------


  implicit none

  private

  public :: one_dim_read_var, &
            read_one_dim_file, &
            two_dim_read_var, &
            read_two_dim_file, &
            fill_blanks_one_dim_vars

  ! Derived type for representing a rank 1 variable that has been read in by one
  ! of the procedures.
  type one_dim_read_var
    character(len=30) :: name
    character(len=30) :: dim_name
    real, dimension(:), pointer :: values
  end type one_dim_read_var

  ! Derived type for representing a rank 2 variable that has been read in by one
  ! of the procedures.
  type two_dim_read_var
    character(len=30) :: name
    character(len=30) :: dim1_name
    character(len=30) :: dim2_name
    real, dimension(:,:), pointer :: values
  end type two_dim_read_var

  contains

  !-------------------------------------------------------------------------------------------------
  subroutine read_two_dim_file( nCol, filename, read_vars, other_dim )
    !
    ! Description: This subroutine reads from a file containing data that varies
    !              in two dimensions. These are dimensions are typically height
    !              and time.
    !
    !-----------------------------------------------------------------------------------------------
    implicit none

    ! Input Variable(s)
    integer, intent(in) :: nCol ! Number of columns expected in the data file

    character(len=*), intent(in) :: filename ! Name of the file being read from

    ! Output Variable(s)
    type (two_dim_read_var), dimension(nCol),intent(out) :: read_vars ! Structured information
    !                                                                   from the file

    type (one_dim_read_var), intent(out) :: other_dim  ! Structured information
    !                                                    on the dimesion not stored in read_vars

    ! Local Variables
    character(len=256),dimension(nCol) :: names ! Names of variables

    integer nRowI ! Inner row

    integer nRowO ! Outer row

    integer :: k, j, i

    real, dimension(ncol) :: tmp

    ! First run through, take names and determine how large the data file is.
    open(unit=10, file=trim(filename), status = 'old' )

    read(10, fmt=*) names

    nRowO = 0
    do while(.true.)
      read(10, *, end=77, err=77) tmp(1), nRowI
      do k =1, nRowI
        read(10, *) tmp
      end do
      nRowO = nRowO + 1
    end do

    77 rewind(10)

    read(10, *) ! Skip the line

    ! Store the names into the structure and allocate accordingly
    do k =1, nCol
      read_vars(k)%name = names(k)
      read_vars(k)%dim1_name = "Time(s)"
      read_vars(k)%dim2_name = names(1)

      allocate( read_vars(k)%values(nRowI, nRowO) )
    end do
    other_dim%name = "Time(s)"
    other_dim%dim_name = "Time(s)"

    allocate( other_dim%values(nRowO) )

    ! Read in the data again to the newly allocated arrays
    do k=1, nRowO
      read(10,*) other_dim%values(k)
      do j=1, nRowI
        read(10,*) ( read_vars(i)%values(j,k), i=1, nCol)
      end do
    end do

    close(10)

    ! Avoiding compiler warning
    if(.false.) print *, tmp

  end subroutine read_two_dim_file

  !-------------------------------------------------------------------------------------------------
  subroutine read_one_dim_file( iunit, nCol, filename, read_vars )
    !
    ! Description: This subroutine reads from a file containing data that varies
    !              in one dimensions. This is typically time.
    !
    !-----------------------------------------------------------------------------------------------
    implicit none

    ! Input Variable(s)

    integer, intent(in) :: iunit ! I/O unit

    integer, intent(in) :: nCol ! Number of columns expected in the data file

    character(len=*), intent(in)  :: filename ! Name of the file being read from

    ! Output Variable(s)
    type (one_dim_read_var), dimension(nCol),intent(out) :: read_vars

    character(len=256),dimension(nCol) :: names

    integer nRow

    integer :: k,j

    real, dimension(ncol) :: tmp

    ! First run through, take names and determine how large the data file is.
    open(unit=iunit, file=trim(filename), status = 'old' )

    read(iunit, fmt=*) names

    nRow = 0
    do while(.true.)
      read(iunit, *, end=77) tmp
      nRow = nRow+1
    end do
    77 continue

    rewind(iunit)
    read(iunit, *) ! Skip the line

    ! Store the names into the structure and allocate accordingly
    do k =1, nCol
      read_vars(k)%name = names(k)
      read_vars(k)%dim_name = names(1)
      allocate( read_vars(k)%values(nRow) )
    end do

    ! Read in the data again to the newly allocated arrays
    do k=1, nRow
      read(iunit,*) ( read_vars(j)%values(k), j=1, nCol)
    end do

    close(iunit)

    ! Avoiding compiler warning
    if(.false.) print *, tmp

  end subroutine read_one_dim_file

  !-------------------------------------------------------------------------------------------------
  subroutine fill_blanks_one_dim_vars( num_vars, one_dim_vars )
    !
    !  Description: This subroutine fills in the blank spots (signified by -999.0)
    !  with values linearly interpolated using the first element of the array as a
    !  guide.
    !
    !-----------------------------------------------------------------------------------------------

    implicit none

    integer, intent(in) :: num_vars ! Number of elements in one_dim_vars

    ! Input/Output Variable(s)
    type(one_dim_read_var), dimension(num_vars), intent(inout):: one_dim_vars ! Read data
    !                                                                           that may have gaps.

    ! Local variables
    integer i

    do i=1, num_vars
      one_dim_vars(i)%values = linear_fill_blanks( size(one_dim_vars(i)%values), &
                                                    one_dim_vars(1)%values, one_dim_vars(i)%values )
    end do

  end subroutine fill_blanks_one_dim_vars

  !-------------------------------------------------------------------------------------------------
  function linear_fill_blanks( dim_grid, grid, var ) &
  !
  !  Description: This function fills blanks in array var using the grid
  !               as a guide. Blank values in var are signified by being
  !               less than or equal to -999.0
  !
  !-----------------------------------------------------------------------------------------------
  result( var_out )

    use interpolation, only: zlinterp_fnc

    implicit none

    ! Input Variable(s)
    integer, intent(in) :: dim_grid ! Size of grid

    real, dimension(dim_grid), intent(in) :: grid ! Array that var is being
    !                                               interpolated to.

    real, dimension(dim_grid), intent(in) :: var ! Array that may contain gaps.

    ! Output Variable(s)
    real, dimension(dim_grid) :: var_out ! Return variable

    ! Local Variables
    real, dimension(dim_grid) :: temp_grid
    real, dimension(dim_grid) :: temp_var

    integer :: i
    integer :: amt


    ! Essentially this code leverages the previously written zlinterp function.
    ! A smaller temporary grid and var variable are being created to pass to
    ! zlinterp. zlinterp then performs the work of taking the temporary var
    ! array and interpolating it to the actual grid array.

    amt = 0
    do i=1, dim_grid
      if(var(i) > -999.0) then
        amt = amt + 1
        temp_var(amt) = var(i)
        temp_grid(amt) = grid(i)
      end if
    end do
    if (amt < dim_grid) then
      var_out = zlinterp_fnc(dim_grid, amt, grid, temp_grid(1:amt), temp_var(1:amt))
    else
      var_out = var
    end if
    return
  end function linear_fill_blanks

  !subroutine deallocate_one_dim_vars( num_vars, one_dim_vars )
  !  integer, intent(in) :: num_vars ! Number of elements in one_dim_vars

  ! Input/Output Variable(s)
  !  type(one_dim_read_var), dimension(num_vars), intent(in):: one_dim_vars ! Read data
  !                                                                           that may have gaps.

  !  integer i

  ! do i=1, num_vars
  !   if(allocated(one_dim_vars(i)%values) ) then
  !     deallocate(one_dim_vars(i)%values)
  !   end if
  ! end do

  ! end subroutine deallocate_one_dim_vars

end module input_reader
