c------------------------------------------------------------------------
c This module contains structure and subroutine definitions to
c create GrADS output data files for one dimensional arrays.
c
c The structure type (outputfile) contains all necessay information
c to generate a GrADS file and a list of variables to be output
c in the data file.
c
c Two subroutines are needed to create a GrADS file
c
c        subroutine open_grads( f )     Open initialize structure
c                                       If GrADS files already exist,
c                                       open_grads will attempt to
c                                       append data to them
c
c        subroutine write_grads( f )    Write data to file and update
c                                       control file. Can be callled as 
c                                       many times as necessary
c
c Chris Golaz, updated 2/18/2003
c
        module outputfile_class

c         Undefined value

          real, parameter :: undef = -9.99e33

c         Structure to hold description of a variable

          type variable
            real, pointer :: ptr(:)               ! pointer to the array
!            character (len = 10) name             ! variable name
            character (len = 15) name             ! variable name
            character (len = 50) description      ! variable decsription
          end type variable

c         Structure to hole description of a GrADS output file

          type outputfile

c           File information

            character (len = 200) fname   ! file name without suffix
            character (len = 200) fdir    ! path where fname resides
            integer iounit                ! logical unit
            integer nrecord               ! number of records written
            integer ntimes                ! number of times written

c           Grid information

            integer ia, iz                ! extent
            real, pointer :: z(:)         ! height of vertical levels

c           Time information

            integer day, month, year      ! date of starting time
            real time                     ! start time in seconds
            real dtwrite                  ! interval between outputs in seconds

c           Variables

            integer nvar                  ! number of variables for this file
            type (variable), pointer :: var(:)
                                          ! list and variable description

          end type outputfile

        contains

          !---------------------------------------------------------
          subroutine open_grads( f, unit, fdir, fname, 
     .                           ia, iz, z,
     .                           day, month, year, time, dtwrite,
     .                           nvar )
          implicit none

c         Input

          type (outputfile), intent(INOUT)   :: f
          integer, intent(IN)               :: unit
          character*(*), intent(IN) :: fdir, fname
          integer, intent(IN) :: ia, iz
          real, intent(IN)    :: z(:)
          integer, intent(IN) :: day, month, year
          real, intent(IN)    :: time, dtwrite
          integer, intent(IN) :: nvar
        
c         Internal

          integer k
          logical lctl, ldat, lerror

c         Define parameters
          
          f%iounit = unit
          f%fdir   = fdir
          f%fname  = fname
          f%ia     = ia
          f%iz     = iz

          if ( ia .le. iz ) then
            do k=1,iz-ia+1
              f%z(k) = z(ia+k-1)
            enddo
          else
            do k=1,ia-iz+1
              f%z(k) = z(ia-k+1)
            enddo
          endif

          f%day = day
          f%month = month
          f%year = year

          f%dtwrite = dtwrite

          f%nvar = nvar

c         Check whether GrADS files already exists

!         We don't need this feature for the single-column model
!         since there is no history restart feature
!
!          inquire( file=trim(fdir)//trim(fname)//'.ctl', exist=lctl )
!          inquire( file=trim(fdir)//trim(fname)//'.dat', exist=ldat )
          lctl = .false.
          ldat = .false.

c         If none of the files exists, set ntimes and nrecord and
c         to initial values and return

          if ( .not.lctl .and. .not.ldat ) then

            f%time = time
            f%ntimes = 0
            f%nrecord = 1
            return

c         If both files exists, attempt to append to existing files

          elseif ( lctl .and. ldat ) then

c           check existing ctl file

            call check_grads( unit, fdir, fname, 
     .                        ia, iz, z,
     .                        day, month, year, time, dtwrite,
     .                        nvar, 
     .                        lerror, f%ntimes, f%nrecord, f%time )

            if ( lerror ) then
              write(*,*) 'Error in open_grads:'
              write(*,*) 'Attempt to append to existing files failed'
!              call stopcode('open_grads')
              stop 'open_grads'
            endif

            return

c         If one file exists, but not the other, give up

          else

              write(*,*) 'Error in open_grads:'
              write(*,*) 'Attempt to append to existing files failed'
              write(*,*) 
     .        'because only one of the two GrADS files was found'
!              call stopcode('open_grads')
              stop 'open_grads'

          endif

          return
          end subroutine open_grads

          !---------------------------------------------------------
          ! For existing GrADS file, this subroutine will attempt
          ! to determine whether data can be safely appended to
          ! existing file.

          subroutine check_grads( unit, fdir, fname, 
     .                            ia, iz, z,
     .                            day, month, year, time, dtwrite,
     .                            nvar, 
     .                            error, ntimes, nrecord, time_in )
          implicit none

c         Input

          integer, intent(IN)               :: unit
          character*(*), intent(IN) :: fdir, fname
          integer, intent(IN) :: ia, iz
          real, intent(IN)    :: z(:)
          integer, intent(IN) :: day, month, year
          real, intent(IN)    :: time, dtwrite
          integer, intent(IN) :: nvar

c         Output

          logical, intent(OUT) :: error
          integer, intent(OUT) :: ntimes, nrecord
          real, intent(OUT)    :: time_in
        
c         Internal

          logical done
          integer ierr
          character (len = 256) :: line, tmp, date, dt

          integer i, nx, ny, nz
          integer ihour, imin

          integer ia_in, iz_in, ntimes_in, nvar_in
          integer day_in, month_in, year_in
          real dtwrite_in
          real, allocatable :: z_in(:)
          type (variable), pointer :: var_in(:)

          integer ios

c Open control file

          error = .false.
          done = .false.
          open( unit,
     .          file = trim(fdir)//trim(fname)//'.ctl',
     .          status = 'old', iostat = ierr )
          if ( ierr < 0 ) done = .true.

c Read and process it

          read(unit,iostat=ierr,fmt='(a256)') line
          if ( ierr < 0 ) done = .true.

          do while ( .not. done )

             if ( index(line,'XDEF') > 0 ) then

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

               read(line,*) tmp, iz_in

               if ( index(line,'LEVELS') > 0 ) then
                 ia_in = 1
                 allocate( z_in(ia_in:iz_in) )
                 read(unit,*) (z_in(i),i=ia_in,iz_in)
               endif

             elseif ( index(line,'TDEF') > 0 ) then

               read(line,*) tmp, ntimes_in, tmp, date, dt
               read(date(1:2),*) ihour
               read(date(4:5),*) imin
               time_in = ihour * 3600. + imin * 60.
               read(date(7:8),*) day_in
               read(date(12:15),*) year_in
               if ( date(9:11) == 'JAN' ) then
                  month_in = 1
               else if ( date(9:11) == 'FEB' ) then
                  month_in = 2
               else if ( date(9:11) == 'MAR' ) then
                  month_in = 3
               else if ( date(9:11) == 'APR' ) then
                  month_in = 4
               else if ( date(9:11) == 'MAY' ) then
                  month_in = 5
               else if ( date(9:11) == 'JUN' ) then
                  month_in = 6
               else if ( date(9:11) == 'JUL' ) then
                  month_in = 7
               else if ( date(9:11) == 'AUG' ) then
                  month_in = 8
               else if ( date(9:11) == 'SEP' ) then
                  month_in = 9
               else if ( date(9:11) == 'OCT' ) then
                  month_in = 10
               else if ( date(9:11) == 'NOV' ) then
                  month_in = 11
               else if ( date(9:11) == 'DEC' ) then
                  month_in = 12
               end if

               read(dt(1:len_trim(dt)-2),*) dtwrite_in
               dtwrite_in = dtwrite_in * 60.

             elseif ( index(line,'ENDVARS') > 0 ) then

               done = .true.

             elseif ( index(line,'VARS') > 0 ) then

               read(line,*) tmp, nvar_in
               allocate( var_in(nvar_in) )
               do i=1,nvar_in
                  read(unit,iostat=ierr,fmt='(a256)') line
                  read(line,*) var_in(i)%name, nz
                  if ( nz /= iz_in ) then
                     write(*,*) 'Error reading ',trim(var_in(i)%name)
                     error = .true.
                  end if
               end do
             endif

             read(unit,iostat=ierr,fmt='(a256)') line
             if ( ierr < 0 ) done = .true.

          end do
          
          close( unit )

c         Perform some error check

          if ( abs(ia_in - iz_in) /= abs(ia - iz) ) then
            write(*,*) 'check_grads: size mismatch'
            error = .true.
          endif

          if ( day_in /= day ) then
            write(*,*) 'check_grads: day mismatch'
            error = .true.
          endif

          if ( month_in /= month ) then
            write(*,*) 'check_grads: month mismatch'
            error = .true.
          endif

          if ( year_in /= year ) then
            write(*,*) 'check_grads: year mismatch'
            error = .true.
          endif

          if ( int( time_in + ntimes_in*dtwrite_in ) 
     .         /= int( time ) ) then
            write(*,*) 'check_grads: time mismatch'
            error = .true.
          endif

          if ( int( dtwrite_in ) /= int( dtwrite) ) then
            write(*,*) 'check_grads: dtwrite mismatch'
            error = .true.
          endif

          if ( nvar_in /= nvar ) then
            write(*,*) 'check_grads: nvar mismatch'
            error = .true.
          endif

          if ( error ) then
            write(*,*) 'check_grads diagnostic'
            write(*,*) 'ia      = ',ia_in,ia
            write(*,*) 'iz      = ',iz_in,iz
            write(*,*) 'day     = ',day_in,day
            write(*,*) 'month   = ',month_in,month
            write(*,*) 'year    = ',year_in,year
            write(*,*) 'time    = ',time_in,time
            write(*,*) 'dtwrite = ',dtwrite_in,dtwrite
            write(*,*) 'nvar    = ',nvar_in,nvar
          endif

c         Set ntimes and nrecord to append to existing files

          ntimes = ntimes_in
          nrecord = ntimes_in * nvar_in * iz_in + 1

c         Searching for a memory leak.  This seems to fixed only 
c         part of the problem. -dschanen 4/4/2005

          deallocate ( z_in ) 

          return
          end subroutine check_grads

          !---------------------------------------------------------
          subroutine write_grads( f )
        
          use endian ! -dschanen 4/4/05

          implicit none

c         Input

          type (outputfile), intent(INOUT) :: f

c         Internal

          integer i,k
          integer ios
!          real*4 xtmp    
          real(kind=4) xtmp    
          character(len=15) date

c         Check number of variables and write nothing if less than 1

          if ( f%nvar < 1 ) return

c         Output data to file

          open( f%iounit,
     .          file = trim(f%fdir)//trim(f%fname)//'.dat',
     .          form = 'unformatted', access = 'direct',
     .          recl = 4, status = 'unknown', iostat = ios )
          if ( ios .gt. 0 ) then
             write(*,*) 'write_grads: error opening binary file'
             write(*,*) 'iostat = ',ios
             stop
          endif
          do i=1,f%nvar
             if ( f%ia .le. f%iz ) then
               do k=f%ia,f%iz
                  xtmp = f%var(i)%ptr(k)
                  write(f%iounit,rec=f%nrecord) xtmp
                  f%nrecord = f%nrecord + 1
               end do
             else
               do k=f%ia,f%iz,-1
                  xtmp = f%var(i)%ptr(k)
                  write(f%iounit,rec=f%nrecord) xtmp
                  f%nrecord = f%nrecord + 1
               end do
             end if
          end do
          close( f%iounit, iostat = ios )
          if ( ios .gt. 0 ) then
             write(*,*) 'write_grads: error closing binary file'
             write(*,*) 'iostat = ',ios
             stop
          endif

          f%ntimes = f%ntimes + 1

c         Write control file

          open( f%iounit,
     .          file = trim(f%fdir)//trim(f%fname)//'.ctl',
     .          status = 'unknown', iostat = ios )
          if ( ios .gt. 0 ) then
             write(*,*) 'write_grads: error opening control file'
             write(*,*) 'iostat = ',ios
             stop
          endif

c         File header

! Added Endian options to file header to avoid confusing gradsc
! 4/4/2005 - dschanen
          if ( big_endian( ) ) then
            write(f%iounit,'(a)') 'OPTIONS BIG_ENDIAN'
          else
            write(f%iounit,'(a)') 'OPTIONS LITTLE_ENDIAN'
          endif 
! end change 4/4/2005
          write(f%iounit,'(a)') 'DSET ^'//trim(f%fname)//'.dat'
          write(f%iounit,'(a,e11.5)') 'UNDEF ',undef
          write(f%iounit,'(a)') 'XDEF    1 LINEAR 0. 1.'
          write(f%iounit,'(a)') 'YDEF    1 LINEAR 0. 1.'
          if ( f%ia .eq. f%iz ) then
            write(f%iounit,'(a,i5,a)') 'ZDEF    1 LINEAR 0. 1.'
          elseif ( f%ia .lt. f%iz ) then
            write(f%iounit,'(a,i5,a)') 
     .        'ZDEF',abs(f%iz-f%ia)+1,' LEVELS '
            write(f%iounit,'(8f8.1)') (f%z(i-f%ia+1),i=f%ia,f%iz)
          else
            write(f%iounit,'(a,i5,a)') 
     .        'ZDEF',abs(f%iz-f%ia)+1,' LEVELS '
            write(f%iounit,'(8f8.1)') (f%z(f%ia-i+1),i=f%ia,f%iz,-1)
          endif

          call formatdate(f%day,f%month,f%year,f%time,date)

          write(f%iounit,'(a,i6,a,a,i3,a)') 'TDEF    '
     .     ,f%ntimes,' LINEAR ',date,max(1,int(f%dtwrite/60.)),'mn'

c         Variables description

          write(f%iounit,'(a,i5)') 'VARS',f%nvar

          do i=1,f%nvar
             write(f%iounit,'(a,i5,a,a)')
     .         f%var(i)%name(1:len_trim(f%var(i)%name)),
     .         abs(f%iz-f%ia)+1,' 99 ',
     .         f%var(i)%description(1:len_trim(f%var(i)%description))
          end do

          write(f%iounit,'(a)') 'ENDVARS'

          close( f%iounit, iostat = ios )
          if ( ios .gt. 0 ) then
             write(*,*) 'write_grads: error closing control file'
             write(*,*) 'iostat = ',ios
             stop
          endif

          return
          end subroutine write_grads

          !---------------------------------------------------------
          subroutine formatdate(iday,imonth,iyear,st_time,date)
          implicit none

          integer iday,imonth,iyear
          real st_time
          character(len=15) date

          character(len=3) month(12)
          data month /'JAN','FEB','MAR','APR','MAY','JUN',
     .                'JUL','AUG','SEP','OCT','NOV','DEC'/

          date = 'hh:mmZddmmmyyyy'
          write(date(7:8),'(i2.2)') iday
          write(date(9:11),'(a3)') month(imonth)
          write(date(12:15),'(i4.4)') iyear
          write(date(1:2),'(i2.2)') int(st_time/3600.)
          write(date(4:5),'(i2.2)') int( mod(int(st_time),3600) / 60. )

          return
          end subroutine formatdate

        end module outputfile_class

c------------------------------------------------------------------------
