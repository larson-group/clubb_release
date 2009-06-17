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
            fill_blanks_one_dim_vars, &
            fill_blanks_two_dim_vars, &
            deallocate_one_dim_vars, &
            deallocate_two_dim_vars, &
            read_x_table

  ! Derived type for representing a rank 1 variable that has been read in by one
  ! of the procedures.
  type one_dim_read_var

    character(len=30) :: name               ! Name of the variable

    character(len=30) :: dim_name           ! Name of the dimension that the
    !                                         variable varies along

    real, dimension(:), pointer :: values   ! Values of that variable

  end type one_dim_read_var

  ! Derived type for representing a rank 2 variable that has been read in by one
  ! of the procedures.
  type two_dim_read_var

    character(len=30) :: name                ! Name of the variable

    character(len=30) :: dim1_name           ! Name of one of the dimensions
    !                                          that the variable varies along

    character(len=30) :: dim2_name           ! Name of the other variable that
    !                                          the variable varies along

    real, dimension(:,:), pointer :: values  ! Values of that variable

  end type two_dim_read_var

  contains

  !-------------------------------------------------------------------------------------------------
  subroutine read_two_dim_file( iunit, nCol, filename, read_vars, other_dim )
    !
    ! Description: This subroutine reads from a file containing data that varies
    !              in two dimensions. These are dimensions are typically height
    !              and time.
    !
    !-----------------------------------------------------------------------------------------------
    implicit none

    ! Input Variable(s)

    integer, intent(in) :: iunit ! File I/O unit

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

    logical :: isComment

    character(len=200) :: tmpline

    real, dimension(ncol) :: tmp

    ! Begin Code

    ! First run through, take names and determine how large the data file is.
    open(unit=iunit, file=trim(filename), status = 'old' )

    isComment = .true.

    ! Skip all the comments at the top of the file
    do while(isComment)
      read(iunit,fmt='(A)') tmpline
      k = index(tmpline, "!")
      isComment = .false.
      if(k > 0) then
        isComment = .true.
      end if
    end do

    ! Go back to the line that wasn't a comment.
    backspace(iunit)

    read(iunit, fmt=*) names

    print *, names

    nRowO = 0
    do while(.true.)
      read(iunit, *, end=77, err=77) tmp(1), nRowI

      if( nRowI < 1 ) then
        stop "Number of elements must be an integer and greater than zero in two-dim  input file."
      end if 

      do k =1, nRowI
        read(iunit, *) tmp
      end do
      nRowO = nRowO + 1
    end do

    77 continue

    do i=1, nRowO
    
      backspace(iunit)
  
      do j=1, nRowI
    
        backspace(iunit)
  
      end do

    end do

    backspace(iunit)

    ! Store the names into the structure and allocate accordingly
    do k =1, nCol
      read_vars(k)%name = names(k)
      read_vars(k)%dim1_name = "Time[s]"
      read_vars(k)%dim2_name = names(1)

      allocate( read_vars(k)%values(nRowI, nRowO) )
    end do

    other_dim%name = "Time[s]"
    other_dim%dim_name = "Time[s]"

    allocate( other_dim%values(nRowO) )

    ! Read in the data again to the newly allocated arrays
    do k=1, nRowO
      read(iunit,*) other_dim%values(k)
      do j=1, nRowI
        read(iunit,*) ( read_vars(i)%values(j,k), i=1, nCol)
      end do
    end do

    close(iunit)


    ! Avoiding compiler warning
    if(.false.) print *, tmp

    return
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

    type (one_dim_read_var), dimension(nCol),intent(out) :: read_vars ! Structured information 
    !                                                                   from the file

    ! Local Variable(s) 
    character(len=256),dimension(nCol) :: names  

    character(len=200) :: tmpline

    integer nRow

    integer :: k, j 

    real, dimension(ncol) :: tmp

    logical :: isComment

    ! Begin Code

    isComment = .true.

    ! First run through, take names and determine how large the data file is.
    open(unit=iunit, file=trim(filename), status = 'old' )

    ! Skip all the comments at the top of the file
    do while(isComment)
      read(iunit,fmt='(A)') tmpline
      k = index(tmpline, "!")
      isComment = .false.
      if(k > 0) then
        isComment = .true.
      end if
    end do

    ! Go back to the line that wasn't a comment.
    backspace(iunit)

    read(iunit, fmt=*) names

    ! Count up the number of rows
    nRow = 0
    do while(.true.)
      read(iunit, *, end=77) tmp
      nRow = nRow+1
    end do
    77 continue

    ! Rewind that many rows
    do k=0, nRow
      backspace(iunit)
    end do

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

    ! Input Variable(s)
    integer, intent(in) :: num_vars ! Number of elements in one_dim_vars

    ! Input/Output Variable(s)
    type(one_dim_read_var), dimension(num_vars), intent(inout):: one_dim_vars ! Read data
    !                                                                           that may have gaps.

    ! Local variable(s)
    integer i

    ! Begin Code

    do i=1, num_vars
      one_dim_vars(i)%values = linear_fill_blanks( size(one_dim_vars(i)%values), &
                                                   one_dim_vars(1)%values, one_dim_vars(i)%values, &
                                                   0.0 )
    end do

  end subroutine fill_blanks_one_dim_vars

  !-------------------------------------------------------------------------------------------------
  subroutine fill_blanks_two_dim_vars( num_vars, other_dim, two_dim_vars )
    !
    !  Description: This subroutine fills in the blank spots (signified by -999.0)
    !  with values linearly interpolated using the first element of the array
    !  and the values in the other_dim argument as a guide.
    !
    !  This is a two step process. First we assume that the other_dim values
    !  have no holes, but there are blanks for that variable across that
    !  dimension. Then we fill holes across the dimension whose values are first
    !  in the array of two_dim_vars.
    !
    !  Ex. Time is the 'other_dim' and Height in meters is the first element in
    !  two_dim_vars.
    !
    !
    !-----------------------------------------------------------------------------------------------

    implicit none

    ! Input Variable(s)

    integer, intent(in) :: num_vars ! Number of elements in one_dim_vars

    ! Input/Output Variable(s)
    type(one_dim_read_var), intent(in):: other_dim ! Read data


    type(two_dim_read_var), dimension(num_vars), intent(inout):: two_dim_vars ! Read data
    !                                                                           that may have gaps.

    ! Local variables
    integer i,j

    integer dim_size
    integer other_dim_size

    ! Begin Code

    dim_size = size(two_dim_vars(1)%values, 1)

    other_dim_size = size(other_dim%values )

    do i=2, num_vars
      ! Interpolate along main dim
      do j=1, other_dim_size
        two_dim_vars(i)%values(:,j) = linear_fill_blanks( dim_size, &
                                               two_dim_vars(1)%values(:,j), &
                                               two_dim_vars(i)%values(:,j), -999.9 )
      end do
      ! Interpopate along other dim
      do j=1, dim_size
        two_dim_vars(i)%values(j,:) = linear_fill_blanks( other_dim_size, &
                                               other_dim%values, two_dim_vars(i)%values(j,:), 0.0 )
      end do


    end do

  end subroutine fill_blanks_two_dim_vars


  !-------------------------------------------------------------------------------------------------
  function linear_fill_blanks( dim_grid, grid, var, default_value ) &
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

    real, intent(in) :: default_value ! Default value if entire profile is -999.9

    ! Output Variable(s)
    real, dimension(dim_grid) :: var_out ! Return variable

    ! Local Variables
    real, dimension(dim_grid) :: temp_grid
    real, dimension(dim_grid) :: temp_var

    integer :: i
    integer :: amt

    ! Begin Code

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


    if( amt == 0 ) then
      var_out = default_value
    else if (amt < dim_grid) then
      var_out = zlinterp_fnc(dim_grid, amt, grid, temp_grid(1:amt), temp_var(1:amt))
    else
      var_out = var
    end if
    return
  end function linear_fill_blanks
  !----------------------------------------------------------------------------
  subroutine deallocate_one_dim_vars( num_vars, one_dim_vars )
    !
    !  Description: This subroutine deallocates the pointer stored in
    !  one_dim_vars%value for the whole array
    !
    !------------------------------------------------------------------------------
    implicit none

    ! Input Variable(s)
    integer, intent(in) :: num_vars ! Number of elements in one_dim_vars

    type(one_dim_read_var), dimension(num_vars), intent(inout):: one_dim_vars ! Read data
    !                                                                        that may have gaps.

    ! External functions
    intrinsic :: associated

    ! Local Variable(s)
    integer i

    ! Begin Code

    do i=1, num_vars

      if(associated(one_dim_vars(i)%values) ) then

        deallocate(one_dim_vars(i)%values)

      end if

    end do

  end subroutine deallocate_one_dim_vars

  !-------------------------------------------------------------------------------------------------
  subroutine deallocate_two_dim_vars( num_vars, two_dim_vars, other_dim )
    !
    !  Description: This subroutine deallocates the pointer stored in
    !  two_dim_vars%value for the whole array
    !
    !-----------------------------------------------------------------------------------------------
    implicit none

    ! Input Variable(s)
    integer, intent(in) :: num_vars ! Number of elements in one_dim_vars

    type(one_dim_read_var), intent(in) :: other_dim
    type(two_dim_read_var), dimension(num_vars), intent(inout):: two_dim_vars ! Read data
    !                                                                        that may have gaps.

    ! External Functions
    intrinsic :: associated

    ! Local Variable(s)
    integer i

    ! Begin Code

    do i=1, num_vars

      if(associated(two_dim_vars(i)%values) ) then

        deallocate(two_dim_vars(i)%values)

      end if

    end do

    if(associated(other_dim%values) ) then

      deallocate(other_dim%values)

    end if



  end subroutine deallocate_two_dim_vars
  !-------------------------------------------------------------------------------------------------
  function read_x_table( nvar, xdim, ydim, target_name, retVars ) result(x)
    !
    !  Description: Searches for the variable specified by target_name in the
    !  collection of retVars. If the function finds the variable then it returns
    !  it. If it does not the program using this function will exit gracefully
    !  with a warning message.
    !
    !-----------------------------------------------------------------------------------------------
    implicit none


    ! Input Variable(s)
    integer, intent(in) :: nvar ! Number of variables in retVars

    integer, intent(in) :: xdim, ydim
    character(len=*), intent(in) :: target_name ! Variable that is being
    !                                             searched for

    type(two_dim_read_var), dimension(nvar), intent(in) :: retVars ! Collection
    !                                                                being searched through

    ! Output Variable(s)
    real, dimension( xdim, ydim ) :: x

    ! Local Variables
    integer i

    logical l_found

    ! Begin Code

    l_found = .false.

    i = 1

    do while( i <= nvar .and. .not. l_found)

      if( retVars(i)%name == target_name ) then
  
        l_found = .true.

        x = retVars(i)%values

      end if

      i=i+1

    end do

    if( .not. l_found ) then

      print *,target_name//' could not be found. '

      stop

    end if

  end function read_x_table
!------------------------------------------------------------------------------
end module input_reader
