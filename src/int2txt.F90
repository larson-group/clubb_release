!-----------------------------------------------------------------------
! $Id$
program int2txt

! Description:
! Takes 34 bytes of binary data as a command-line argument,
!   and outputs ASCII text.
! Used by the tuner; invoked by tune/generate_seed.bash.
! Notes:  Not very useful on machines without a /dev/random, 
!   though any 34 bytes of pseudo-random seed data would suffice.
!-----------------------------------------------------------------------
#include "CLUBB_core/recl.inc"

  use mt95, only: &
    genrand_intg ! Integer precision

#ifdef AbsoftUNIXFortran
  use unix_library, only: iargc, getarg
#endif
  implicit none

#ifdef __GFORTRAN__
#define getarg get_command_argument
#define iargc command_argument_count

#elif AbsoftUNIXFortran

#else
  ! Seems to work on most Fortran compilers, but if it doesn't,
  ! try using the F2003 subrountine is get_command_argument() instead.
  external :: getarg

  ! As above, except the F2003 function is command_argument_count()
  integer, external :: iargc
#endif

  intrinsic :: trim

  ! Parameter Constants
  integer, parameter ::  & 
    fstderr = 0, fstdout = 6, seed_dim = 34

  ! Local Variables
  integer(kind=genrand_intg), dimension(seed_dim) :: seed ! Our seed data

  character(len=50) :: rand_source

!-----------------------------------------------------------------------

  ! Test to check there is an argument given
  if ( iargc( ) < 1 ) then
    write(fstderr,*) "Usage: int2txt <filename>"
    error stop
  end if

  ! Get the source of the random data (usually /dev/random)
  call getarg(1, rand_source)

  open(unit=10, file=trim( rand_source ), action='read',  & 
    recl=F_RECL*seed_dim, form='unformatted', access='direct')

  read(unit=10, rec=1) seed(1:seed_dim)

  close(unit=10)

  write(unit=fstdout,fmt='(34I12)') seed(1:seed_dim)

end program int2txt
!-----------------------------------------------------------------------
