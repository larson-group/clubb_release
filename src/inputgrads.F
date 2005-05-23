!-----------------------------------------------------------------------
! This module contains structure and subroutine definitions to
! open a GrADS data file a read it.
!
! The structure (inputgrads) contains all necessary information
! about the file to read a specific variable at  specific time.
!
! Chris Golaz, 9/12/2000

! Minor modifications, dschanen 4/19/2005

        module inputfile_class

        use endian

        implicit none

! Structure to hold description of a variable

          type variable
            integer index
            character (len = 10) name             ! variable name
          end type variable

! Structure to hold description of a GrADS input file

          type inputgrads

!           File information

            character (len = 128) fname           ! binary data file name
            integer iounit                        ! unit #
            logical lbyteswapped                  ! need to swap bytes

!           Grid information

            integer ia, iz                        ! extent
            real, pointer :: z(:)                 ! height of levels

!           Time information

            integer day, month, year    ! date of starting time
            real time                   ! start time in seconds
            real dtwrite                ! interval between outputs in seconds
            integer ntimes              ! number of time records in file

!           Variables

            integer nvar
            type (variable), pointer :: var(:)

          end type inputgrads

        contains

!-----------------------------------------------------------------------
          subroutine open_grads_read( unit, fname, f )

          implicit none

!         Input / Output

          integer, intent(IN)              :: unit
          character(len=*), intent(IN)     :: fname
          type (inputgrads), intent(INOUT) :: f

!         Internal

          logical               :: done, error
          integer               :: ierr
          character (len = 256) :: line, tmp, date, dt

          integer               :: i, nx, ny, nz
          integer               :: ihour, imin

! Open control file
          open( unit, file = trim( fname ), status = 'old' )

!  Initialize status booleans
          f%lbyteswapped = .false.
          error          = .false.
          done           = .false.

! Read and process it
          read(unit,iostat=ierr,fmt='(a256)') line
          if ( ierr < 0 ) done = .true.

          do while ( .not. done )

             if ( index(line,'DSET') > 0 ) then
        
               read(line,*) tmp, f%fname
               if ( f%fname(1:1) == '^' ) then
                  f%fname = f%fname(2:len_trim(f%fname)) 
                  i = index( fname, '/', back = .true. )
                  if ( i > 0 ) then
                     f%fname = fname(1:i) // f%fname
                  end if
               endif

             elseif ( index(line,'BYTESWAPPED') > 0 ) then

               f%lbyteswapped = .true.

             elseif ( index(line,'BIG_ENDIAN') > 0 ) then

!              Swap bytes if local machine is little_endian and file
!              big_endian

               if ( little_endian( ) ) f%lbyteswapped = .true.

             elseif ( index(line,'LITTLE_ENDIAN') > 0 ) then

!              Swap bytes if local machine is big_endian and file
!              little_endian

               if ( big_endian( ) ) f%lbyteswapped = .true.

             elseif ( index(line,'XDEF') > 0 ) then

               read(line,*) tmp, nx
               if ( nx /= 1 ) then
                  write(*,*) 'Error: XDEF can only be 1'
                  error = .true.
               end if

             elseif ( index(line,'YDEF') > 0 ) then

               read(line,*) tmp, ny
               if ( ny /= 1 ) then
                  write(*,*) 'Error: YDEF can only be 1'
                  error = .true.
               end if

             elseif ( index(line,'ZDEF') > 0 ) then

               read(line,*) tmp, f%iz
               f%ia = 1
               allocate( f%z(f%ia:f%iz) )
               read(unit,*) (f%z(i),i=f%ia,f%iz)

             elseif ( index(line,'TDEF') > 0 ) then

               read(line,*) tmp, f%ntimes, tmp, date, dt
               read(date(1:2),*) ihour
               read(date(4:5),*) imin
               f%time = ihour * 3600. + imin * 60.
               read(date(7:8),*) f%day
               read(date(12:15),*) f%year
               if ( date(9:11) == 'JAN' ) then
                  f%month = 1
               else if ( date(9:11) == 'FEB' ) then
                  f%month = 2
               else if ( date(9:11) == 'MAR' ) then
                  f%month = 3
               else if ( date(9:11) == 'APR' ) then
                  f%month = 4
               else if ( date(9:11) == 'MAY' ) then
                  f%month = 5
               else if ( date(9:11) == 'JUN' ) then
                  f%month = 6
               else if ( date(9:11) == 'JUL' ) then
                  f%month = 7
               else if ( date(9:11) == 'AUG' ) then
                  f%month = 8
               else if ( date(9:11) == 'SEP' ) then
                  f%month = 9
               else if ( date(9:11) == 'OCT' ) then
                  f%month = 10
               else if ( date(9:11) == 'NOV' ) then
                  f%month = 11
               else if ( date(9:11) == 'DEC' ) then
                  f%month = 12
               end if

               read(dt(1:len_trim(dt)-2),*) f%dtwrite
               f%dtwrite = f%dtwrite * 60.

             elseif ( index(line,'ENDVARS') > 0 ) then

               done = .true.

             elseif ( index(line,'VARS') > 0 ) then

               read(line,*) tmp, f%nvar
               allocate( f%var(f%nvar) )
               do i=1,f%nvar
                  read(unit,iostat=ierr,fmt='(a256)') line
                  read(line,*) f%var(i)%name, nz
                  if ( nz /= f%iz ) then
                     write(*,*) 'Error reading ',trim(f%var(i)%name)
                     error = .true.
                  end if
                  f%var(i)%index = i
               end do
             endif

             read(unit,iostat=ierr,fmt='(a256)') line
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
             write(*,*) 
     .         'Fatal error encountered while reading control file'
             write(*,*) 'Cannot do miracles...'
             stop
          end if

! Open binary file for direct access

          f%iounit = unit
          open( f%iounit,
     .          file = trim(f%fname),
     .          form = 'unformatted', access = 'direct',
     .          recl = 4, status = 'old' )

          return

          end subroutine open_grads_read

!----------------------------------------------------------------------
          subroutine get_var( f, varname, itime, x, error )

          implicit none

!         Input
          type (inputgrads), intent(IN) :: f
          character(len=*)              :: varname
          integer, intent(IN)           :: itime

!         Output
          real, intent(OUT)             :: x(:)
          logical, intent(OUT)          :: error

!         Internal
          logical                       :: flag, done
          integer                       :: i, k, nrec, ivar

!         Initialize error to false
          error = .false.

!         Check time index

          if ( itime < 1 .or. itime > f%ntimes ) then
             error = .true.
!             write(*,*) 'get_var: itime < 1 .or. itime > f%ntimes'
!             write(*,*) 'itime = ',itime
!             write(*,*) 'f%ntimes = ',f%ntimes
             return
          end if

!         Look up variable in list
          done = .false.
          i    = 1
          do while ( .not. done )
             if ( trim(varname) == trim(f%var(i)%name) ) then
                ivar = i
                done = .true.
             else
                i = i + 1
                if ( i > f%nvar ) done = .true.
             end if
          end do

          if ( i > f%nvar ) then
             error = .true.
!             write(*,*) 'get_var: i > f%nvar'
!             write(*,*) 'i = ',i
!             write(*,*) 'f%nvar = ',f%nvar
             return
          end if

!         Read variable from file

          nrec = (itime-1)*f%nvar*(f%iz-f%ia+1) 
     .           + (ivar-1)*(f%iz-f%ia+1) + 1

          do k=f%ia,f%iz
             read(f%iounit,rec=nrec) x(k)
             if (f%lbyteswapped) call native_4byte_real( x(k) )
             nrec = nrec + 1
          end do

          return

          end subroutine get_var

          subroutine close_grads_read( f )

          implicit none

!         Input

          type (inputgrads), intent(IN) :: f

!         Close file 
          close( f%iounit )

!         Deallocate
          deallocate( f%var )
          deallocate( f%z )

          return

          end subroutine close_grads_read
!----------------------------------------------------------------------

        end module inputfile_class
