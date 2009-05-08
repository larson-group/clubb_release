program sclr_sounding_convert
!
!       Simply this program makes default *_sclr_sounding.in and
!       *_edsclr_sounding.in files from a *_sounding.in file. Where * is the
!       runtype.
!
!
!
!--------------------------------------------------------------------------------

  ! Arbitrary max number of sounding points
  integer, parameter :: nmaxsnd = 600

  integer, parameter :: sclr_max = 1000

  ! Derived type for representing a rank 1 variable that has been read in by one
  ! of the procedures.
  type one_dim_read_var
    character(len=30) :: name
    character(len=30) :: dim_name
    real, dimension(:), pointer :: values
  end type one_dim_read_var
  ! Name of the case being converted.
  character(len=200) :: casename

  ! Dummy variable
  character(len=200) :: theta_type

  integer :: i, j
  integer :: sclr_dim = 2
  integer :: nlevels
  integer :: edsclr_dim = 2

  type(one_dim_read_var), dimension(7) :: retVars ! Collection being

  ! Sounding Information
  real, dimension(nmaxsnd, sclr_max) :: &
        sclr, edsclr


  integer, parameter :: iunit = 43
  integer, parameter :: ounit = 44

  open(unit = ounit, file="toconvert", status="old")

  do while( .true. ) ! While there are names in toconvert list
    read(ounit, end=55, fmt=*) casename ! Get the casename to convert

    call read_one_dim_file( iunit, 7, trim(casename)//'_sounding.in', retVars )
    nlevels = size(retVars(1)%values)

    sclr(:,1) = read_x_profile(7,'rt[kg\kg]', retVars)
    call read_theta_profile(7, retVars, theta_type, sclr(:,2))

    ! Write out namelist in new format
    open(unit = iunit, file = trim(casename)//'_sclr_sounding.in', status = 'new')
    write(iunit, * ) "rt[kg\kg]       ", "        thl[K]"
    do i=1, nlevels
      write(iunit, *) (sclr(i,j),j=1, sclr_dim)
    end do

    close(iunit)
    ! Write out namelist in new format
    open(unit = iunit, file = trim(casename)//'_edsclr_sounding.in', status = 'new')
    write(iunit, * ) "rt[kg\kg]       ", "        thl[K]"
    do i=1, nlevels
      write(iunit, *) (sclr(i,j),j=1, sclr_dim)
    end do

    close(iunit)


  end do
  55 close(ounit)

  do i=1, 7
    deallocate(retVars(i)%values)
  end do
  contains

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
  function read_x_profile( nvar, target_name, retVars ) result(x)
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

    character(len=*), intent(in) :: target_name ! Variable that is being
    !                                             searched for

    type(one_dim_read_var), dimension(nvar), intent(in) :: retVars ! Collection
    !                                                                being searched through

    ! Output Variable(s)
    real, dimension(nmaxsnd) :: x

    ! Local Variables
    integer i

    logical l_found

    print *, target_name
    l_found = .false.

    i = 1
    do while( i <= nvar .and. .not. l_found)
      if( retVars(i)%name == target_name ) then
        l_found = .true.
        x(1:size(retVars(i)%values)) = retVars(i)%values
      end if
      i=i+1
    end do

    if( .not. l_found ) then
      stop ' Profile could not be found. Check your sounding.in file.'
    end if

  end function read_x_profile
!-------------------------------------------------------------------------------------------------
  subroutine read_theta_profile(nvar, retVars, theta_type, theta)
    !
    !  Description: Searches for the variable specified by either 'thetal[K]' or 'theta[K]' in the
    !  collection of retVars. If the function finds the variable then it returns
    !  it. If it does not the program using this function will exit gracefully
    !  with a warning message.
    !
    !-----------------------------------------------------------------------------------------------

    implicit none

    character(len=*), parameter :: theta_name = 'theta[K]'

    character(len=*), parameter :: thetal_name = 'thetal[K]'
    ! Input Variable(s)
    integer, intent(in) :: nvar ! Number of elements in retVars

    type(one_dim_read_var), dimension(nvar), intent(in) :: retVars ! Collection being
    !                                                                searched through

    character(len=*), intent(out) :: theta_type

    ! Output Variable(s)
    real, dimension(nmaxsnd), intent(out) :: theta

    if( count( (/ any(retVars%name == theta_name), any(retVars%name == thetal_name) /)) <= 1) then
      if( any(retVars%name == theta_name))then
        theta_type = theta_name
      elseif( any(retVars%name == thetal_name))then
        theta_type = thetal_name
      else
        stop "Could not read theta compatable variable"
      endif
      theta = read_x_profile(nvar, theta_type, retVars)

    end if
  end subroutine read_theta_profile

end program sclr_sounding_convert
