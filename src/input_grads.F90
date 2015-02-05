!-------------------------------------------------------------------------------
! $Id$

module input_grads

! Description:
!   This module contains structure and subroutine definitions to
!   open a GrADS data file a read it.
!
! References:
!   None
!
! Original Author:
!   Chris Golaz, 9/12/2000
!
! Modifications:
!   * Uses functions rather than subroutines to get endian type.
!   * Other cosmetic changes.
!   * Added preprocesing for RECL
! Other modifications have been tracked via subversion.
!-------------------------------------------------------------------------------
#include "CLUBB_core/recl.inc"
  use endian, only: & 
    little_endian, & ! Variable(s)
    big_endian, & 
    byte_order_swap ! Procedure

  use constants_clubb, only:  & 
    fstdout,  & ! Variable(s) 
    fstderr,  &
    var_length

  use clubb_precision, only:  & 
    time_precision ! Variable(s)

  implicit none

  private ! Default Scope

  public :: get_grads_var, open_grads_read,  & 
            close_grads_read


  ! Overloaded interface for get_grads_var.  All GrADS files are assumed
  ! to store variable as 4 byte IEEE floats, but the model may be
  ! using double or extended precision.

  contains

!-------------------------------------------------------------------------------
  subroutine open_grads_read( unit_number, fname, grads_file, l_error )

! Description:
!   Open a GrADS data set in read-only mode
! References:
!   None
!-------------------------------------------------------------------------------

    use model_flags, only: l_byteswap_io

    use stat_file_module, only: stat_file

    use text_writer, only: write_text

    use constants_clubb, only: &
      sec_per_day, & ! Constant(s)
      sec_per_hr, &
      sec_per_min

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Constant Parameters
    logical, parameter :: &
      l_extra_debugging = .false.

    ! Input Variables
    integer, intent(in) :: &
      unit_number ! Fortran I/O unit

    character(len=*), intent(in) ::  & 
      fname ! The file name

    ! Input / Output
    type (stat_file), intent(inout) ::  & 
      grads_file ! Derived data type containing info on the GrADS file

    ! Output Variable(s)
    logical, intent(out) :: l_error

    ! Local Variables
    logical :: l_done, l_file_exist
    integer :: ierr

    character(len=256) ::  & 
      line, tmp, date, dt

    integer ::  & 
      i, nx, ny, nzmax, & 
      ihour, imin

!-------------------------------------------------------------------------------

    ! ---- Begin Code ----

    !  Initialize status booleans
    grads_file%l_byte_swapped = .false.
    l_error          = .false.
    l_done           = .false.
    l_file_exist     = .false.

    ! Check if the file we're trying to open exists, preventing runtime error
    inquire( file=fname, exist=l_file_exist )

    if( .not.( l_file_exist ) ) then
      ! if the file doesn't exist, return to the calling function to handle the error
      write(fstderr,*) "File "//trim( fname )// " not found."
      l_error = .true.
      return
    end if

    ! Open control file
    open( unit=unit_number, file=trim( fname ), status = 'old' )

    ! Read and process it
    read(unit_number,iostat=ierr,fmt='(a256)') line
    if ( ierr < 0 ) l_done = .true.

    do while ( .not. l_done )

      if ( index(line,'DSET') > 0 ) then

        read(unit=line,fmt=*) tmp, grads_file%fname
        if ( .false. ) dt = tmp ! Dummy code to eliminate a compiler warning
        if ( grads_file%fname(1:1) == '^' ) then
          ! Get the name of the associated .dat file
          grads_file%fname = grads_file%fname(2:len_trim(grads_file%fname))
          ! Figure out the file path
          i = index( fname, '/', back = .true. )
          if ( i > 0 ) then
            ! Construct a path for the .date file
            grads_file%fname = fname(1:i) // grads_file%fname
          end if
        end if

      else if ( index(line,'BYTESWAPPED') > 0 ) then

        grads_file%l_byte_swapped = .true.

      else if ( index(line,'BIG_ENDIAN') > 0 ) then

        ! Swap bytes if local machine is little_endian and file
        ! big_endian

        if ( little_endian .or. ( l_byteswap_io .and. big_endian ) ) then
          grads_file%l_byte_swapped = .true.
        end if

      else if ( index(line,'LITTLE_ENDIAN') > 0 ) then

        ! Swap bytes if local machine is big_endian and file
        ! little_endian

        if ( big_endian .or. ( l_byteswap_io .and. little_endian ) ) then
          grads_file%l_byte_swapped = .true.
        end if

      else if ( index(line,'XDEF') > 0 ) then

        read(unit=line,fmt=*) tmp, nx
        if ( nx /= 1 ) then
          write(unit=fstderr,fmt=*) 'Error: XDEF can only be 1'
          l_error = .true.
        end if

      else if ( index(line,'YDEF') > 0 ) then

        read(unit=line,fmt=*) tmp, ny
        if ( ny /= 1 ) then
          write(unit=fstderr,fmt=*) 'Error: YDEF can only be 1'
          l_error = .true.
        end if

      else if ( index(line,'ZDEF') > 0 ) then

        read(unit=line,fmt=*) tmp, grads_file%iz
        grads_file%ia = 1
        allocate( grads_file%z(grads_file%ia:grads_file%iz) )
        ! Implied Do Loop with the purpose of reading in
        ! altitudes
        if ( grads_file%iz == 1 ) then
          grads_file%z(1) = 1._core_rknd
        else
          read(unit=unit_number,fmt=*) (grads_file%z(i),i=grads_file%ia,grads_file%iz)
        end if
      else if ( index(line,'TDEF') > 0 ) then

        read(unit=line,fmt=*) tmp, grads_file%ntimes, tmp, date, dt
        read(unit=date(1:2),fmt=*) ihour
        read(unit=date(4:5),fmt=*) imin

        grads_file%time = real( ihour, kind=time_precision ) * &
                          real( sec_per_hr,kind=time_precision) &
                        + real( imin, kind=time_precision ) *  &
                          real(sec_per_min, kind=time_precision)

        read(unit=date(7:8),fmt=*) grads_file%day
        read(unit=date(12:15),fmt=*) grads_file%year

        select case( date(9:11) )
        case( 'JAN' )
          grads_file%month = 1
        case( 'FEB' )
          grads_file%month = 2
        case( 'MAR' )
          grads_file%month = 3
        case( 'APR' )
          grads_file%month = 4
        case( 'MAY' )
          grads_file%month = 5
        case( 'JUN' )
          grads_file%month = 6
        case( 'JUL' )
          grads_file%month = 7
        case( 'AUG' )
          grads_file%month = 8
        case( 'SEP' )
          grads_file%month = 9
        case( 'OCT' )
          grads_file%month = 10
        case( 'NOV' )
          grads_file%month = 11
        case( 'DEC' )
          grads_file%month = 12
        case default
          write(unit=fstderr,fmt=*) "Unknown month: "//date(9:11)
          l_error = .true.
        end select

        i = len_trim( dt )

        ! Read time variable from the string
        read(dt(1:i-2),*) grads_file%dtwrite

        ! Determine units on time
        select case ( dt(i-1:i) )
        case ( 'mn' )
          grads_file%dtwrite = grads_file%dtwrite * sec_per_min
        case ( 'hr' )
          grads_file%dtwrite = grads_file%dtwrite * sec_per_hr
        case ( 'dy' )
          grads_file%dtwrite = grads_file%dtwrite * sec_per_day
        case default
          write(unit=fstderr,fmt=*) "Unknown time increment: "//dt(i-1:i)
          l_error = .true.
        end select

      else if ( index(line,'ENDVARS') > 0 ) then

        l_done = .true.

      else if ( index(line,'VARS') > 0 ) then

        read(unit=line,fmt=*) tmp, grads_file%nvar
        allocate( grads_file%var(grads_file%nvar) )

        do i = 1, grads_file%nvar, 1
          
          read(unit=unit_number,iostat=ierr,fmt='(a256)') line
          read(unit=line,fmt=*) grads_file%var(i)%name, nzmax
          
          if ( nzmax /= grads_file%iz ) then
            write(unit=fstderr,fmt=*) "Error reading ",  & 
              trim( grads_file%var(i)%name )
            l_error = .true.
          end if

          grads_file%var(i)%indx = i

        end do ! 1..grads_file%nvar


      end if

      read(unit=unit_number,iostat=ierr,fmt='(a256)') line
      if ( ierr < 0 ) l_done = .true.

    end do

    close( unit_number )

!--------- Debugging information -----------------------------------------------
    if ( l_extra_debugging ) then
      write(*,*) 'grads_file%fname = ',trim(grads_file%fname)
      write(*,*) 'grads_file%l_byte_swapped = ',grads_file%l_byte_swapped
      write(*,*) 'grads_file%ia = ',grads_file%ia
      write(*,*) 'grads_file%iz = ',grads_file%iz
      write(*,'(8f8.1)') (grads_file%z(i),i=grads_file%ia,grads_file%iz)
      write(*,*) 'grads_file%ntimes = ',grads_file%ntimes
      write(*,*) 'grads_file%day = ',grads_file%day
      write(*,*) 'grads_file%month = ',grads_file%month
      write(*,*) 'grads_file%year = ',grads_file%year
      write(*,*) 'grads_file%time = ',grads_file%time
      write(*,*) 'grads_file%dtwrite = ',grads_file%dtwrite
      write(*,*) 'grads_file%nvar = ',grads_file%nvar
      do i=1,grads_file%nvar
         write(*,*) trim(grads_file%var(i)%name)
      end do
    end if ! l_extra_debugging
!--------- End debugging information -------------------------------------------

    if ( l_error ) then
      write(unit=fstderr,fmt=*)  & 
        'Fatal error encountered while reading control file in open_grads_read'
      write(unit=fstderr,fmt=*) 'Cannot do miracles...'
      return
    end if

    ! Open binary file for direct access

    grads_file%iounit = unit_number
    inquire( file = grads_file%fname, exist = l_file_exist )
    if ( .not. l_file_exist ) then
      write(fstderr,*) 'binary GrADS file does not exist'
      l_error = .true.
      return
    end if

    open( unit = grads_file%iounit, & 
          file = trim( grads_file%fname ), & 
          form = 'unformatted', access = 'direct',  & 
          recl = F_RECL, status='old', iostat=ierr )

    if ( ierr /= 0 ) then
      write(unit=fstderr,fmt=*)  & 
        "input_grads: error opening binary file"
      write(unit=fstderr,fmt=*) "iostat = ", ierr, "filename = ", grads_file%fname
      l_error = .true.
      return
    end if

    return
  end subroutine open_grads_read

!-------------------------------------------------------------------------------
  subroutine get_grads_var( grads_file, varname, itime, variable, l_error )

! Description:
!   Read binary data from a 1D GrADS file and return the result 'variable'
!   in CLUBB's default precision.
!-------------------------------------------------------------------------------
    use constants_clubb, only: fstderr, sec_per_min ! Constant(s)

    use stat_file_module, only: stat_file ! Variable

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: floor, trim, real, nint, selected_real_kind

    ! Constant Parameters
    integer, parameter :: r4 = selected_real_kind( p=5 )

    ! Input Variables
    type (stat_file), intent(in) :: & 
      grads_file ! The file to read data from

    character(len=*), intent(in) ::  & 
      varname ! The variable name as it occurs in the control file

    integer, intent(in) :: & 
      itime ! Obtain variable varname at time itime [min]

    ! Output Variables
    real( kind = core_rknd ), dimension(:), intent(out) ::  & 
      variable ! Result variable

    logical, intent(out) :: l_error

    ! Local Variables
    real(kind=r4), dimension(size( variable )) :: grads_variable

    logical :: l_done
    integer :: i, k, nrec, ivar

    ! ---- Begin Code ----

    ! Initialize l_error to false
    l_error = .false.

    ! Check time index
    ! Now assumes itime is in minutes
    if ( itime < 1 .or. (itime/floor( grads_file%dtwrite/sec_per_min )) > grads_file%ntimes ) then
      l_error = .true.
      write(unit=fstderr,fmt=*)  & 
        "get_4byte_var: itime < 1 .or. itime > grads_file%ntimes"
      write(unit=fstderr,fmt=*) "itime = ", itime
      write(unit=fstderr,fmt=*) "grads_file%ntimes = ", grads_file%ntimes

      return
    end if

    ! Look up variable in list
    l_done = .false.
    i    = 1
    ivar = -1 ! Initialization to avoid a compiler warning
    do while ( .not. l_done )
      if ( trim( varname ) == trim( grads_file%var(i)%name ) ) then
        ivar = i
        l_done = .true.
      else
        i = i + 1
        if ( i > grads_file%nvar ) l_done = .true.
      end if
    end do ! .not. l_done

    if ( i > grads_file%nvar ) then
      l_error = .true.
!     write(*,*) 'get_4byte_var: i > grads_file%nvar'
!     write(*,*) 'i = ',i
!     write(*,*) 'grads_file%nvar = ',grads_file%nvar
      write(fstderr,*) "input_grads get_var: "//trim( varname ), " variable not found."
      return
    end if

    ! Read variable from file

    ! dschanen changed this to take into account varying dtwrite
    ! numbers 22 March 2007
!         nrec = (itime-1)*grads_file%nvar*(grads_file%iz-grads_file%ia+1)
!    .         + (ivar-1)*(grads_file%iz-grads_file%ia+1) + 1

!          print *, "Division check", nint(itime/(grads_file%dtwrite/60._core_rknd))-1
!          print *, "grads_file%nvar", grads_file%nvar
!          print *, "varindex", ivar-1
!          print *, "nlevels", (grads_file%iz-grads_file%ia+1)

    ! Probably not the most elegant to round but it does allow
    ! cases like arm_3year with their default _stats.in and
    ! _model.in to restart
    ! Joshua Fasching March 2008

    nrec = ( max( nint( real( itime, kind=core_rknd) &
                        /(grads_file%dtwrite/sec_per_min) &
                      ), & ! nint &
              1 ) & ! max
            -1 ) &
      *grads_file%nvar*(grads_file%iz-grads_file%ia+1)  & 
         + (ivar-1)*(grads_file%iz-grads_file%ia+1)
    nrec = nrec + 1

    ! Debug
!          print *, varname
!          print *, "ivar = ", ivar
!          print *, "nrec = ", nrec

    do k=grads_file%ia,grads_file%iz
      read(unit=grads_file%iounit,rec=nrec) grads_variable(k)
      if ( grads_file%l_byte_swapped ) call byte_order_swap( grads_variable(k) )
      nrec = nrec + 1
    end do

    variable = real( grads_variable, kind = core_rknd ) ! Convert to default precision

    return
  end subroutine get_grads_var

!-------------------------------------------------------------------------------
  subroutine close_grads_read( grads_file )

! Description:
!   Close a previously opened GrADS file
!-------------------------------------------------------------------------------

    use stat_file_module, only: stat_file

    implicit none

    ! Input/Output Variables
    type (stat_file), intent(inout) :: &
      grads_file ! Derived data type with information on the GrADS file

!-------------------------------------------------------------------------------

    ! ---- Begin Code ----

    ! Close file
    close( unit=grads_file%iounit )

    ! Deallocate
    deallocate( grads_file%var )
    deallocate( grads_file%z )

    return
  end subroutine close_grads_read

end module input_grads
