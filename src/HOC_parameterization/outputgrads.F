!-----------------------------------------------------------------------
! $Id: outputgrads.F,v 1.8 2008-04-04 22:02:49 dschanen Exp $
#ifdef STATS 
        module outputgrads

!       Description:
!       This module contains structure and subroutine definitions to
!       create GrADS output data files for one dimensional arrays.
!
!       The structure type (outputfile) contains all necessay information
!       to generate a GrADS file and a list of variables to be output
!       in the data file.
!
!       Two subroutines are needed to create a GrADS file
!
!        subroutine open_grads( f )     Open initialize structure
!                                       If GrADS files already exist,
!                                       open_grads will attempt to
!                                       append data to them
!
!        subroutine write_grads( f )    Write data to file and update
!                                       control file. Can be callled as 
!                                       many times as necessary
!       Author:
!        Chris Golaz, updated 2/18/2003
!-----------------------------------------------------------------------

          use constants, only: fstdout, fstderr

          public  :: open_grads, write_grads
          private :: formatdate, check_grads

          ! Undefined value
          real, private, parameter :: undef = -9.99e33

          contains

!-----------------------------------------------------------------------
          subroutine open_grads( f, unit, fdir, fname, 
     .                           ia, iz, z,
     .                           day, month, year, rlat, rlon,
     .                           time, dtwrite,
     .                           nvar )

!-----------------------------------------------------------------------
          use constants, only: fstderr, fstdout
          use outputfile_class

          implicit none

!         ! Input Variables

          integer, intent(IN) :: unit
          character*(*), intent(IN) :: fdir, fname

          integer, intent(IN) :: ia, iz
          real, dimension(:), intent(IN) :: z

          integer, intent(IN) :: day, month, year

          real, intent(IN) :: rlat, rlon, dtwrite

          real(kind=time_prec), intent(IN) :: time 

          integer, intent(IN) :: nvar

          ! Input/Output Variables
          type (outputfile), intent(INOUT) :: f
        
          ! Local Variables

          integer :: k
          logical :: lctl, ldat, lerror

          ! Define parameters
          
          f%iounit = unit
          f%fdir   = fdir
          f%fname  = fname
          f%ia     = ia
          f%iz     = iz

          if ( ia <= iz ) then
            do k=1,iz-ia+1
              f%z(k) = z(ia+k-1)
            end do
          else
            do k=1,ia-iz+1
              f%z(k) = z(ia-k+1)
            end do
          end if

          f%day   = day
          f%month = month
          f%year  = year

          f%rlat  = rlat
          f%rlon  = rlon

          f%dtwrite = dtwrite

          f%nvar = nvar

!         Check whether GrADS files already exists

!         We don't need this feature for the single-column model
!         since there is no history restart feature
!
!          inquire( file=trim(fdir)//trim(fname)//'.ctl', exist=lctl )
!          inquire( file=trim(fdir)//trim(fname)//'.dat', exist=ldat )
          lctl = .false.
          ldat = .false.

          ! If none of the files exist, set ntimes and nrecord and
          ! to initial values and return

          if ( .not.lctl .and. .not.ldat ) then

            f%time = time
            f%ntimes = 0
            f%nrecord = 1
            return

          ! If both files exists, attempt to append to existing files

          else if ( lctl .and. ldat ) then

            !  Check existing ctl file

            call check_grads( unit, fdir, fname, 
     .                        ia, iz, z,
     .                        day, month, year, time, dtwrite,
     .                        nvar, 
     .                        lerror, f%ntimes, f%nrecord, f%time )

            if ( lerror ) then
              write(unit=fstderr,fmt=*) "Error in open_grads:"
              write(unit=fstderr,fmt=*) 
     .        "Attempt to append to existing files failed"
!              call stopcode('open_grads')
              stop 'open_grads'
            end if

            return

!         If one file exists, but not the other, give up

          else
            write(unit=fstderr,fmt=*) 'Error in open_grads:'
            write(unit=fstderr,fmt=*) 
     .        "Attempt to append to existing files failed,"// 
     .        " because only one of the two GrADS files was found."
            stop "open_grads"

          end if

          return
          end subroutine open_grads

!-----------------------------------------------------------------------
          subroutine check_grads( unit, fdir, fname, 
     .                            ia, iz, z,
     .                            day, month, year, time, dtwrite,
     .                            nvar, 
     .                            error, ntimes, nrecord, time_in )
!         Description:
!         For existing GrADS file, this subroutine will attempt
!         to determine whether data can be safely appended to
!         existing file.
!-----------------------------------------------------------------------
          use outputfile_class
          use stats_prec
          use constants, only: fstderr, fstdout

          implicit none

          ! Input Variables

          integer, intent(in) :: 
     .    unit,             ! Fortran file unit
     .    ia, iz,           ! First and last level
     .    day, month, year, ! Day, month and year numbers
     .    nvar              ! Number of variables in the file

          character*(*), intent(in) ::
     .    fdir, fname ! File directory and name

          real, dimension(:), intent(in) ::
     .    z(:) ! Altitudes              [m]
          real(kind=time_prec), intent(in) ::
     .    time    ! Current time        [s]
          real, intent(in) ::
     .    dtwrite ! Time interval between writes to the file    [s]

          ! Output Variables
          logical, intent(out) :: 
     .    error

          integer, intent(out) :: 
     .    ntimes, nrecord

          real(kind=time_prec), intent(out) :: time_in
        
          ! Local Variables
          logical :: done
          integer :: ierr
          character(len = 256) :: line, tmp, date, dt

          integer :: 
     .    i, nx, ny, nz,
     .    ihour, imin,
     .    ia_in, iz_in, ntimes_in, nvar_in,
     .    day_in, month_in, year_in,
     .    ios

          real :: dtwrite_in

          real, dimension(:), allocatable :: z_in

          type (variable), dimension(:), pointer :: var_in

!-----------------------------------------------------------------------

          ! Initialize logicals
          error = .false.
          done  = .false.

          ! Open control file
          open( unit,
     .          file = trim( fdir )//trim( fname )//'.ctl',
     .          status = 'old', iostat = ierr )
          if ( ierr < 0 ) done = .true.

          ! Read and process it

          read(unit=unit,iostat=ierr,fmt='(a256)') line
          if ( ierr < 0 ) done = .true.

          do while ( .not. done )

             if ( index(line,'XDEF') > 0 ) then

               read(unit=line,fmt=*) tmp, nx
               if ( nx /= 1 ) then
                 write(unit=fstderr,fmt=*) 'Error: XDEF can only be 1'
                 error = .true.
               end if

             else if ( index(line,'YDEF') > 0 ) then

               read(unit=line,fmt=*) tmp, ny
               if ( ny /= 1 ) then
                 write(unit=fstderr,fmt=*) "Error: YDEF can only be 1"
                 error = .true.
               end if

             else if ( index(line,'ZDEF') > 0 ) then

               read(unit=line,fmt=*) tmp, iz_in

               if ( index(line,'LEVELS') > 0 ) then
                 ia_in = 1
                 allocate( z_in(ia_in:iz_in) )
                 read(unit=unit,fmt=*) (z_in(i),i=ia_in,iz_in)
               end if

             else if ( index(line,'TDEF') > 0 ) then

               read(unit=line,fmt=*) tmp, ntimes_in, tmp, date, dt
               read(unit=date(1:2),fmt=*) ihour
               read(unit=date(4:5),fmt=*) imin
               time_in = ihour * 3600.0 + imin * 60.0
               read(unit=date(7:8),fmt=*) day_in
               read(unit=date(12:15),fmt=*) year_in

               select case( date(9:11) )
               case( 'JAN' )
                 month_in = 1
               case( 'FEB' )
                 month_in = 2
               case( 'MAR' )
                 month_in = 3
               case( 'APR' )
                 month_in = 4
               case( 'MAY' )
                 month_in = 5
               case( 'JUN' )
                 month_in = 6
               case( 'JUL' )
                 month_in = 7
               case( 'AUG' )
                 month_in = 8
               case( 'SEP' )
                 month_in = 9
               case( 'OCT' )
                 month_in = 10
               case( 'NOV' )
                 month_in = 11
               case( 'DEC' )
                 month_in = 12
               case default
                 write(unit=fstderr,fmt=*) "Unknown month: "//date(9:11)
                 error = .true.
               end select

               read(unit=dt(1:len_trim(dt)-2),fmt=*) dtwrite_in
               dtwrite_in = dtwrite_in * 60.0

             else if ( index(line,'ENDVARS') > 0 ) then

               done = .true.

             else if ( index(line,'VARS') > 0 ) then

               read(line,*) tmp, nvar_in
               allocate( var_in(nvar_in) )
               do i=1, nvar_in
                  read(unit=unit,iostat=ierr,fmt='(a256)') line
                  read(unit=line,fmt=*) var_in(i)%name, nz
                  if ( nz /= iz_in ) then
                    write(unit=fstderr,fmt=*) 
     .                "Error reading ", trim( var_in(i)%name )
                     error = .true.
                  end if ! nz /= iz_in
               end do ! 1..nvar_in
             end if

             read(unit,iostat=ierr,fmt='(a256)') line
             if ( ierr < 0 ) done = .true.

          end do ! while ( .not. done )
          
          close( unit=unit )

          ! Perform some error check

          if ( abs(ia_in - iz_in) /= abs(ia - iz) ) then
            write(unit=fstderr,fmt=*) "check_grads: size mismatch"
            error = .true.
          end if

          if ( day_in /= day ) then
            write(unit=fstderr,fmt=*) "check_grads: day mismatch"
            error = .true.
          end if

          if ( month_in /= month ) then
            write(unit=fstderr,fmt=*) "check_grads: month mismatch"
            error = .true.
          end if

          if ( year_in /= year ) then
            write(unit=fstderr,fmt=*) "check_grads: year mismatch"
            error = .true.
          end if

          if ( int( time_in + ntimes_in*dtwrite_in ) 
     .         /= int( time ) ) then
            write(unit=fstderr,fmt=*) "check_grads: time mismatch"
            error = .true.
          end if

          if ( int( dtwrite_in ) /= int( dtwrite) ) then
            write(unit=fstderr,fmt=*) 'check_grads: dtwrite mismatch'
            error = .true.
          end if

          if ( nvar_in /= nvar ) then
            write(unit=fstderr,fmt=*) 'check_grads: nvar mismatch'
            error = .true.
          end if

          if ( error ) then
            write(unit=fstderr,fmt=*) "check_grads diagnostic"
            write(unit=fstderr,fmt=*) "ia      = ", ia_in, ia
            write(unit=fstderr,fmt=*) "iz      = ", iz_in, iz
            write(unit=fstderr,fmt=*) "day     = ", day_in, day
            write(unit=fstderr,fmt=*) "month   = ", month_in, month
            write(unit=fstderr,fmt=*) "year    = ", year_in, year
            write(unit=fstderr,fmt=*) "time    = ", time_in, time
            write(unit=fstderr,fmt=*) "dtwrite = ", dtwrite_in, dtwrite
            write(unit=fstderr,fmt=*) "nvar    = ", nvar_in, nvar
          end if

          ! Set ntimes and nrecord to append to existing files

          ntimes  = ntimes_in
          nrecord = ntimes_in * nvar_in * iz_in + 1

          deallocate ( z_in ) 

          return
          end subroutine check_grads

!---------------------------------------------------------
          subroutine write_grads( f )

!         Description:
!         Write part of a GrADS file to disk
!---------------------------------------------------------
        
          use endian
          use outputfile_class

          implicit none

          ! Input Variables
          type (outputfile), intent(inout) :: f

          ! Local Variables
          integer :: 
     .    i, k, ! Loop indices
     .    ios   ! I/O status

          character(len=15) :: date

          ! Check number of variables and write nothing if less than 1

          if ( f%nvar < 1 ) return

#include "../recl.F"

          ! Output data to file
          open( unit=f%iounit,
     .          file=trim( f%fdir )//trim( f%fname )//'.dat',
     .          form='unformatted', access='direct',
     .          recl=F_RECL*abs( f%iz-f%ia+1 ),
     .          status='unknown', iostat=ios )
!         open( f%iounit,
!    .          file = trim(f%fdir)//trim(f%fname)//'.dat',
!    .          form = 'unformatted', access = 'direct',
!    .          recl = F_RECL, status = 'unknown', iostat = ios )
          if ( ios /= 0 ) then
            write(unit=fstderr,fmt=*) 
     .        "write_grads: error opening binary file"
            write(unit=fstderr,fmt=*) "iostat = ", ios
            stop
          end if

          if ( f%ia <= f%iz ) then
            do i=1,f%nvar
              write(f%iounit,rec=f%nrecord) 
     .          real( f%var(i)%ptr(f%ia:f%iz), kind=4)
              f%nrecord = f%nrecord + 1
            end do

          else
            do i=1, f%nvar
              write(f%iounit,rec=f%nrecord)
     .          real( f%var(i)%ptr(f%ia:f%iz:-1), kind=4)
              f%nrecord = f%nrecord + 1 
            end do

          end if ! f%ia <= f%iz

          close( f%iounit, iostat = ios )
          if ( ios /= 0 ) then
            write(unit=fstderr,fmt=*) 
     .        "write_grads: error closing binary file"
            write(unit=fstderr,fmt=*) "iostat = ", ios
            stop
          end if

          f%ntimes = f%ntimes + 1

          ! Write control file

          open(unit=f%iounit,
     .         file=trim( f%fdir )//trim( f%fname )//'.ctl',
     .         status='unknown', iostat=ios)
          if ( ios > 0 ) then
            write(unit=fstderr,fmt=*) 
     .        "write_grads: error opening control file"
            write(unit=fstderr,fmt=*) "iostat = ", ios
            stop
          end if

          ! Write file header
          if ( big_endian( ) ) then
            write(unit=f%iounit,fmt='(a)') 'OPTIONS BIG_ENDIAN'
          else
            write(unit=f%iounit,fmt='(a)') 'OPTIONS LITTLE_ENDIAN'
          end if 
          write(unit=f%iounit,fmt='(a)') 'DSET ^'//trim(f%fname)//'.dat'
          write(unit=f%iounit,fmt='(a,e11.5)') 'UNDEF ',undef
          write(unit=f%iounit,fmt='(a)') 'XDEF    1 LINEAR 0. 1.'
          write(unit=f%iounit,fmt='(a)') 'YDEF    1 LINEAR 0. 1.'
          if ( f%ia == f%iz ) then
            write(unit=f%iounit,fmt='(a,i5,a)') 'ZDEF    1 LINEAR 0. 1.'
          else if ( f%ia < f%iz ) then
            write(unit=f%iounit,fmt='(a,i5,a)') 
     .        'ZDEF', abs(f%iz-f%ia)+1,' LEVELS '
            write(unit=f%iounit,fmt='(8f8.1)') 
     .        (f%z(i-f%ia+1),i=f%ia,f%iz)
          else
            write(unit=f%iounit,fmt='(a,i5,a)') 
     .        'ZDEF',abs(f%iz-f%ia)+1,' LEVELS '
            write(f%iounit,'(8f8.1)') (f%z(f%ia-i+1),i=f%ia,f%iz,-1)
          end if

          call formatdate(f%day,f%month,f%year,f%time,date)

          write(unit=f%iounit,fmt='(a,i6,a,a,i5,a)') 'TDEF    ',
     .      f%ntimes, ' LINEAR ', date, max(1,int(f%dtwrite/60.)),'mn'

          ! Variables description
          write(unit=f%iounit,fmt='(a,i5)') 'VARS', f%nvar

          do i=1, f%nvar, 1
             write(unit=f%iounit,fmt='(a,i5,a,a)')
     .         f%var(i)%name(1:len_trim(f%var(i)%name)),
     .         abs(f%iz-f%ia)+1,' 99 ',
     .         f%var(i)%description(1:len_trim(f%var(i)%description))
          end do

          write(unit=f%iounit,fmt='(a)') 'ENDVARS'

          close( unit=f%iounit, iostat=ios )
          if ( ios > 0 ) then
            write(unit=fstderr,fmt=*) 
     .        "write_grads: error closing control file"
            write(unit=fstderr,fmt=*) "iostat = ",ios
            stop
          end if

          return
          end subroutine write_grads

          !---------------------------------------------------------
          subroutine formatdate( day_in, month_in, year_in, 
     .                           time_in, date )
          use stats_prec, only: time_prec
          implicit none

          ! Constant Parameters
          character(len=3), dimension(12), parameter ::
     .    month = (/'JAN','FEB','MAR','APR','MAY','JUN',
     .              'JUL','AUG','SEP','OCT','NOV','DEC'/)

          integer, dimension(24), parameter ::
     .    ndays = (/31, 28, 31, 30, 31, 30,
     .              31, 31, 30, 31, 30, 31,
     .              31, 29, 31, 30, 31, 30,
     .              31, 31, 30, 31, 30, 31 /)

          ! Input Variables
          integer, intent(in) :: 
     .    day_in, month_in, year_in

          real(kind=time_prec), intent(in) :: time_in

          ! Output Variables
          character(len=15), intent(out) :: date

          ! Local Variables
          integer :: iday, imonth, iyear, ileap
          real :: time


!         Copy input arguments into internal variables

          iday   = day_in
          imonth = month_in
          iyear  = year_in
          time   = time_in

!         If starting time is greater than 24 hours, update
!         date accordingly

          do while ( int(time) >= 86400 )

            time = time - 86400.0
            iday = iday + 1

!           determine if current year is leap year
!           no:  ileap = 0
!           yes: ileap = 12

!  Picky compilers will crap out on the code below
!  Something to do with logical .true./.false. not always being 1/0?
!           ileap = abs( (mod(iyear,4).eq.0)
!    .                    .and. (.not.(  mod(iyear,100).eq.0
!    .                                   .and. mod(iyear,400).ne.0 )) 
!    .                 )*12

!   Alternate code:
            if ( ( mod( iyear, 4 ) == 0)
     .            .and. (.not.( mod( iyear, 100 ) == 0 
     .                        .and. mod( iyear, 400 ) /= 0 ) )
     .         ) then 
              ileap = 12
            else
              ileap = 0
            end if

            if ( iday > ndays(imonth+ileap) ) then
              iday = 1
              imonth = imonth + 1
            end if

            if ( imonth > 12 ) then
              imonth = 1
              iyear = iyear + 1
            end if

          end do

          date = 'hh:mmZddmmmyyyy'
          write(unit=date(7:8),fmt='(i2.2)') iday
          write(unit=date(9:11),fmt='(a3)') month(imonth)
          write(unit=date(12:15),fmt='(i4.4)') iyear
          write(unit=date(1:2),fmt='(i2.2)') int( time/3600. )
          write(unit=date(4:5),fmt='(i2.2)') 
     .      int( mod( int( time ),3600 ) / 60. )

          return
          end subroutine formatdate

        end module outputgrads
#endif /*STATS*/
!------------------------------------------------------------------------
