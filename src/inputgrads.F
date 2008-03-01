!-----------------------------------------------------------------------
! $Id: inputgrads.F,v 1.15 2008-03-01 18:58:23 faschinj Exp $

        module inputfile_class

!       Description:
!       This module contains structure and subroutine definitions to
!       open a GrADS data file a read it.
!
!       The structure (inputgrads) contains all necessary information
!       about the file to read a specific variable at  specific time.
!
!       Chris Golaz, 9/12/2000

!       Modifications:
!       * Uses functions rather than subroutines to get endian type. 
!       * Other cosmetic changes.
!       * Overloaded subroutine get_var to allow for 8 byte real output.
!       * Added preprocesing for RECL 
!-----------------------------------------------------------------------
#include "recl.F"
        use endian
        use constants, only: fstdout, fstderr

        implicit none
        ! Overloaded interface for get_var.  All GrADS files are assumed
        ! to store variable as 4 byte IEEE floats, but the model may be
        ! using double or extended precision.
        interface get_var
          module procedure get_4byte_var, get_8byte_var
        end interface

        ! Structure to hold description of a variable
        type variable
          integer :: index
          character(len=15) :: name ! variable name
        end type variable

        ! Structure to hold description of a GrADS input file

        type inputgrads

          ! File information
          character(len=128) :: fname ! Binary data file name
          integer :: iounit           ! Fortran unit #
          logical :: lbyteswapped     ! Need to swap bytes? (T/F)

          ! Grid information
          integer ia, iz              ! Vertical extent
          real, pointer :: z(:)       ! Height of levels        [m]

          ! Time information
          integer day, month, year    ! date of starting time

          real :: 
     .    time,    ! Start time                 [s]
     .    dtwrite  ! Interval between output    [s]

          integer :: ntimes ! Number of time records in file

          ! Local Variables
          integer :: nvar

          type (variable), pointer :: var(:)

        end type inputgrads

        contains

!-----------------------------------------------------------------------
          subroutine open_grads_read( unit, fname, f )

!         Description:
!         Open a GrADS data set in read-only mode
!-----------------------------------------------------------------------
          implicit none

          ! Input Variables
          integer, intent(in) :: unit ! Fortran I/O unit

          character(len=*), intent(in) :: 
     .    fname ! The file name

          ! Input / Output
          type (inputgrads), intent(inout) :: 
     .    f ! The GrADS file

          ! Local Variables
          logical :: done, error
          integer :: ierr

          character(len=256) :: 
     .    line, tmp, date, dt

          integer :: 
     .    i, nx, ny, nz,
     .    ihour, imin

!-----------------------------------------------------------------------
          !  Initialize status booleans
          f%lbyteswapped = .false.
          error          = .false.
          done           = .false.

          ! Open control file
          open( unit=unit, file=trim( fname ), status = 'old' )

          ! Read and process it
          read(unit,iostat=ierr,fmt='(a256)') line
          if ( ierr < 0 ) done = .true.

          do while ( .not. done )

             if ( index(line,'DSET') > 0 ) then
        
               read(unit=line,fmt=*) tmp, f%fname
               if ( f%fname(1:1) == '^' ) then
                  ! Get the name of the associated .dat file
                  f%fname = f%fname(2:len_trim(f%fname)) 
                  ! Figure out the file path
                  i = index( fname, '/', back = .true. )
                  if ( i > 0 ) then
                     ! Construct a path for the .date file
                     f%fname = fname(1:i) // f%fname
                  end if
               end if

             else if ( index(line,'BYTESWAPPED') > 0 ) then

               f%lbyteswapped = .true.

             else if ( index(line,'BIG_ENDIAN') > 0 ) then

               ! Swap bytes if local machine is little_endian and file
               ! big_endian

               if ( little_endian( ) ) f%lbyteswapped = .true.

             else if ( index(line,'LITTLE_ENDIAN') > 0 ) then

               ! Swap bytes if local machine is big_endian and file
               ! little_endian

               if ( big_endian( ) ) f%lbyteswapped = .true.

             else if ( index(line,'XDEF') > 0 ) then

               read(unit=line,fmt=*) tmp, nx
               if ( nx /= 1 ) then
                 write(unit=fstderr,fmt=*) 'Error: XDEF can only be 1'
                 error = .true.
               end if

             else if ( index(line,'YDEF') > 0 ) then

               read(unit=line,fmt=*) tmp, ny
               if ( ny /= 1 ) then
                 write(unit=fstderr,fmt=*) 'Error: YDEF can only be 1'
                 error = .true.
               end if

             else if ( index(line,'ZDEF') > 0 ) then

               read(unit=line,fmt=*) tmp, f%iz
               f%ia = 1
               allocate( f%z(f%ia:f%iz) )
               ! Implied Do Loop with the purpose of reading in
               ! altitudes
               read(unit=unit,fmt=*) (f%z(i),i=f%ia,f%iz)

             else if ( index(line,'TDEF') > 0 ) then

               read(unit=line,fmt=*) tmp, f%ntimes, tmp, date, dt
               read(unit=date(1:2),fmt=*) ihour
               read(unit=date(4:5),fmt=*) imin

               f%time = ihour * 3600. + imin * 60.

               read(unit=date(7:8),fmt=*) f%day
               read(unit=date(12:15),fmt=*) f%year

               select case( date(9:11) )
               case( 'JAN' )
                 f%month = 1
               case( 'FEB' )
                 f%month = 2
               case( 'MAR' )
                 f%month = 3
               case( 'APR' )
                 f%month = 4
               case( 'MAY' )
                 f%month = 5
               case( 'JUN' )
                 f%month = 6
               case( 'JUL' )
                 f%month = 7
               case( 'AUG' )
                 f%month = 8
               case( 'SEP' )
                 f%month = 9
               case( 'OCT' )
                 f%month = 10
               case( 'NOV' )
                 f%month = 11
               case( 'DEC' )
                 f%month = 12
               case default
                 write(unit=fstderr,fmt=*) "Unknown month: "//date(9:11)
                 error = .true.
               end select
               ! Assumes minutes
               read(dt(1:len_trim(dt)-2),*) f%dtwrite
               f%dtwrite = f%dtwrite * 60.

             else if ( index(line,'ENDVARS') > 0 ) then

               done = .true.

             else if ( index(line,'VARS') > 0 ) then

               read(unit=line,fmt=*) tmp, f%nvar
               allocate( f%var(f%nvar) )

               do i = 1, f%nvar, 1

                  read(unit=unit,iostat=ierr,fmt='(a256)') line
                  read(unit=line,fmt=*) f%var(i)%name, nz

                  if ( nz /= f%iz ) then
                     write(unit=fstderr,fmt=*) "Error reading ", 
     .                 trim( f%var(i)%name )
                     error = .true.
                  end if

                  f%var(i)%index = i

               end do ! 1..f%nvar

             end if

             read(unit=unit,iostat=ierr,fmt='(a256)') line
             if ( ierr < 0 ) done = .true.

          end do
          
          close( unit )

!--------- Debug -------------------------------------------------------
!         write(*,*) 'f%fname = ',trim(f%fname)
!         write(*,*) 'f%lbyteswapped = ',f%lbyteswapped
!         write(*,*) 'f%ia = ',f%ia
!         write(*,*) 'f%iz = ',f%iz
!         write(*,'(8f8.1)') (f%z(i),i=f%ia,f%iz)
!         write(*,*) 'f%ntimes = ',f%ntimes
!         write(*,*) 'f%day = ',f%day
!         write(*,*) 'f%month = ',f%month
!         write(*,*) 'f%year = ',f%year
!         write(*,*) 'f%time = ',f%time
!         write(*,*) 'f%dtwrite = ',f%dtwrite
!         write(*,*) 'f%nvar = ',f%nvar
!         do i=1,f%nvar
!            write(*,*) trim(f%var(i)%name)
!         end do
!--------- Debug -------------------------------------------------------
          if ( error ) then
             write(unit=fstderr,fmt=*) 
     .         'Fatal error encountered while reading control file'
             write(unit=fstderr,fmt=*) 'Cannot do miracles...'
             stop
          end if

        ! Open binary file for direct access

          f%iounit = unit
          open( unit = f%iounit,
     .          file = trim( f%fname ),
     .          form = 'unformatted', access = 'direct', 
     .          recl = F_RECL, status='old', iostat=ierr )

          if ( ierr /= 0 ) then
            write(unit=fstderr,fmt=*) 
     .        "input_grads: error opening binary file"
            write(unit=fstderr,fmt=*) "iostat = ", ierr
            stop
          end if

          return
          end subroutine open_grads_read

!----------------------------------------------------------------------
          subroutine get_4byte_var( f, varname, itime, x, error )

!         Description:
!         Read binary data from file units and return the result as 
!         as 4 byte float 'x'
!----------------------------------------------------------------------

          implicit none

          ! Input Variables
          type (inputgrads), intent(in) ::
     .    f ! The file to read data from

          character(len=*) :: 
     .    varname ! The variable name as it occurs in the control file
          
          integer, intent(in) ::
     .    itime ! Obtain variable varname at time itime [m]

          ! Output
          real(kind=4), dimension(:), intent(out) :: 
     .    x ! Result variable

          logical, intent(out) :: error

          ! Internal
          logical :: flag, done
          integer :: i, k, nrec, ivar

          ! Initialize error to false
          error = .false.

          ! Check time index
          ! Now assumes itime is in minutes
          if ( itime < 1 .or. (itime/(f%dtwrite/60.)) > f%ntimes ) then
            error = .true.
            write(unit=fstderr,fmt=*) 
     .        "get_var: itime < 1 .or. itime > f%ntimes"
            write(unit=fstderr,fmt=*) "itime = ", itime
            write(unit=fstderr,fmt=*) "f%ntimes = ", f%ntimes

            return
          end if
                    
          ! Look up variable in list
          done = .false.
          i    = 1
          do while ( .not. done )
             if ( trim( varname ) == trim( f%var(i)%name ) ) then
               ivar = i
               done = .true.
             else
               i = i + 1
               if ( i > f%nvar ) done = .true.
             end if
          end do ! .not. done

          if ( i > f%nvar ) then
            error = .true.
!            write(*,*) 'get_var: i > f%nvar'
!            write(*,*) 'i = ',i
!            write(*,*) 'f%nvar = ',f%nvar
            return
          end if

          ! Read variable from file

          ! dschanen changed this to take into account varying dtwrite
          ! numbers 22 March 2007
!         nrec = (itime-1)*f%nvar*(f%iz-f%ia+1) 
!    .         + (ivar-1)*(f%iz-f%ia+1) + 1

!          print *, "Division check", itime/(f%dtwrite/60.)-1
!          print *, "f%nvar", f%nvar
!          print *, "varindex", ivar-1
!          print *, "nlevels", (f%iz-f%ia+1)
          nrec = (itime/(f%dtwrite/60.)-1)*f%nvar*(f%iz-f%ia+1) 
     .         + (ivar-1)*(f%iz-f%ia+1)
          nrec = nrec + 1
          ! Debug
!          print *, varname
!          print *, "ivar = ", ivar
!          print *, "nrec = ", nrec
          
          do k=f%ia,f%iz
            read(unit=f%iounit,rec=nrec) x(k)
            if ( f%lbyteswapped ) call byte_order_swap( x(k) )
            nrec = nrec + 1
          end do

          return
          end subroutine get_4byte_var

!----------------------------------------------------------------------
          subroutine get_8byte_var( f, varname, itime, x, error )

!         Description:
!         Takes the result from get_4byte_var and returns it as a double
!         precision type, allowing for compile time promotion.
!----------------------------------------------------------------------

          implicit none

          ! External
          intrinsic :: dble, size

          ! Input Variables
          type (inputgrads), intent(in) :: 
     .    f ! The GrADS file

          character(len=*) :: 
     .    varname ! The variable name as it occurs in the control file

          integer, intent(in) ::
     .    itime   ! Obtain variable varname at time itime 

          ! Output Variables
          real(kind=8), intent(out) :: x(:)

          logical, intent(out) :: error

          ! Internal
          real(kind=4), dimension(size( x )) :: tmp

          call get_4byte_var( f, varname, itime, tmp, error )

          x = dble( tmp )

          return
          end subroutine get_8byte_var

!-----------------------------------------------------------------------
          subroutine close_grads_read( f )

!         Description:
!         Close a previously opened GrADS file
!-----------------------------------------------------------------------

          implicit none

          ! Input Variables
          type (inputgrads), intent(INOUT) :: f

!-----------------------------------------------------------------------
          ! Close file 
          close( unit=f%iounit )

          ! Deallocate
          deallocate( f%var )
          deallocate( f%z )

          return
          end subroutine close_grads_read

        end module inputfile_class
