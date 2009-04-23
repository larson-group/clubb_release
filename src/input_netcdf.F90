!-----------------------------------------------------------------------
! $Id$

module input_netcdf

  use constants, only:  & 
    fstdout,  & ! Variable(s) 
    fstderr

  implicit none

  private ! Default Scope

  public :: get_var, open_netcdf_read,  & 
    close_netcdf_read

  contains

!-----------------------------------------------------------------------
  subroutine open_netcdf_read( path, ncid )

! Description:
!   Open a netCDF data set in read-only mode
!-----------------------------------------------------------------------
    use netcdf, only: nf90_open, nf90_strerror ! Function(s)

    use netcdf, only: NF90_NOWRITE, NF90_NOERR ! Constant(s)

    implicit none

    ! Input Variables
    character(len=*), intent(in) ::  & 
      path ! The file name

    ! Output Variables
    integer, intent(out) :: ncid ! The netCDF id of the file

    ! Local Variables
    integer :: ierr

!-----------------------------------------------------------------------

    ! ---- Begin Code ----

    ierr = nf90_open( path=trim( path ), mode=NF90_NOWRITE, ncid=ncid )

    if ( ierr /= NF90_NOERR ) then
      write(fstderr,*) "Error opening netCDF file: "//trim( path )
      write(fstderr,*) nf90_strerror( ierr )
      stop "Fatal error."
    end if

    return
  end subroutine open_netcdf_read

!----------------------------------------------------------------------
  subroutine get_var( ncid, varname, itime, x, l_error )

! Description:
!   Read in values from a netCDF file.
! Assumptions:
!   We assume the variable in question is 4 bytes long
!----------------------------------------------------------------------

    use netcdf, only: &
      nf90_inq_varid, & ! Function(s)
      nf90_get_var, &
      nf90_inquire_variable, &
      nf90_strerror, &
      nf90_inquire_dimension

    use netcdf, only: &
      NF90_NOERR, & ! Constant(s)
      NF90_MAX_VAR_DIMS, &
      NF90_FLOAT, &
      NF90_MAX_NAME

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      ncid ! NetCDF file id.

    character(len=*), intent(in) ::  & 
      varname ! The variable name as it occurs in the control file

    integer, intent(in) :: & 
      itime ! Obtain variable varname at time itime [m]

    ! Output
    real, dimension(:), intent(out) ::  & 
      x ! Result variable
     
    logical, intent(out) :: l_error

    ! Local Variables
    integer, dimension(NF90_MAX_VAR_DIMS) :: dimIds

    integer :: varid, ierr, tdim, xdim, ydim, zdim, i, itmp, xtype, ndims

    character(len=NF90_MAX_NAME) :: dim_name

    real(kind=4), dimension(:,:,:,:), allocatable :: & 
      x4 ! The variable from the file


    ! ---- Begin Code ----

    ! Initialize l_error to false
    l_error = .false.

    ierr = nf90_inq_varid( ncid=ncid, name=varname, varid=varid )

    if ( ierr /= NF90_NOERR ) then
      write(fstderr,*) nf90_strerror( ierr )
      l_error = .true.
      return
    end if

    ierr = nf90_inquire_variable( ncid=ncid, varid=varid, xtype=xtype, ndims=ndims, dimIds=dimIds )

    if ( ierr /= NF90_NOERR .or. ndims /= 4 .or. xtype /= NF90_FLOAT ) then
            write(fstderr,*) "input_netcdf.get_var : The netCDF data doesn't "// &
            "conform to expected precision, shape, or dimensions"
      l_error = .true.
      return

    else
      xdim = -1
      ydim = -1
      zdim = -1
      tdim = -1
      do i = 1, ndims
        ierr = nf90_inquire_dimension( ncid=ncid, dimId=dimIds(i), name=dim_name, len=itmp )
        select case ( trim( dim_name ) )
        case ( "Z", "z", "altitude", "height" )
          zdim = itmp
        case ( "X", "x", "longitude" )
          xdim = itmp
        case ( "Y", "y", "latitude" )
          ydim = itmp
        case ( "T", "t", "time" )
          tdim = itmp
        case default
          l_error = .true.
          return
        end select

      end do

      if ( tdim > itmp ) then
        write(fstderr,*) "input_netcdf.get_var: The time specified exceeds the netCDF data"
        l_error = .true.
        return
      end if

      allocate( x4(xdim,ydim,zdim,1) )

    end if

    ! Obtain a vertical column at time itime
    ierr = nf90_get_var( ncid=ncid, varid=varid, values=x4, &
      start=(/1,1,1,itime/), count=(/1,1,zdim,1/) )

    if ( ierr /= NF90_NOERR ) then
      write(fstderr,*) nf90_strerror( ierr )
      l_error = .true.
      return
    else
      x = real( x4(1,:,1,1) )
    end if

    return
  end subroutine get_var


!-----------------------------------------------------------------------
  subroutine close_netcdf_read( ncid )

! Description:
!   Close a previously opened netCDF file
!-----------------------------------------------------------------------
    use netcdf, only: nf90_close, nf90_strerror ! Function(s)

    use netcdf, only: NF90_NOERR ! Constant(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      ncid
    integer :: ierr
!-----------------------------------------------------------------------
    ! Close file
    ierr = nf90_close( ncid=ncid )

    if ( ierr /= NF90_NOERR ) then
      write(fstderr,*) "Error closing a netCDF file"
      write(fstderr,*) nf90_strerror( ierr )
      stop "Fatal error."
    end if

    return
  end subroutine close_netcdf_read

end module input_netcdf
