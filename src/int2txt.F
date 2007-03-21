!-----------------------------------------------------------------------
! $Id: int2txt.F,v 1.3 2007-03-21 22:28:43 dschanen Exp $
        program int2txt

!       Description:
!       Takes Input of 34 bytes of binary data, and outputs ASCII text.
!       Notes:  Not very useful on machines without a /dev/random, 
!       though any 34 bytes of pseudo-random seed data would suffice.
!-----------------------------------------------------------------------
#include "recl.F"

        implicit none

        ! Seems to work on most Fortran compilers, but if it doesn't,
        ! try using the F2003 subrountine is get_command_argument() instead.
        external :: getarg

        ! As above, except the F2003 function is command_argument_count()
        integer, external :: iargc

        ! Parameter Constants
        integer, parameter :: 
     .  fstderr = 0, fstdout = 6, seed_dim = 34

        ! Local Variables
        integer(kind=4), dimension(seed_dim) :: seed ! Our seed data

        character(len=50) :: rand_source

        ! Test to check there is an argument given
        if ( iargc( ) < 1 ) then
          write(fstderr,*) "Usage: int2txt <filename>"
          stop
        end if

        ! Get the source of the random data (usually /dev/random)
        call getarg(1, rand_source)

        open(unit=10, file=trim( rand_source ), action='read', 
     .       recl=F_RECL*seed_dim, form='unformatted', access='direct')

        read(unit=10, rec=1) seed(1:seed_dim)

        close(unit=10)

        write(unit=fstdout,fmt='(34I12)') seed(1:seed_dim)

        end program int2txt
