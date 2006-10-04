!-----------------------------------------------------------------------
! $Id: int2txt.F,v 1.1 2006-10-04 16:48:04 dschanen Exp $
        program int2txt

!       Description:
!       Take a binary file of 34 4 byte integers and generate an ASCII 
!       text file from it.
!       Takes Input of  34 bytes of binary data;  "rand_seed.bin".
!       Outputs ASCII text; "rand_seed.dat".
!       Notes:  Not very useful on machines without a /dev/random, 
!       though any 34 bytes of pseudo-random seed data will suffice.
!-----------------------------------------------------------------------

! The same rules as inputgrads
#if defined (__i386__) && __INTEL_COMPILER >= 900 /* Newer ifort on i386 */
#define F_RECL 1
#elif defined (__x86_64__) && __INTEL_COMPILER >= 900 /*Newer ifort on x86_64*/
#define F_RECL 1
#elif __alpha /* Assume that alpha users are using f95 */
#define F_RECL 1
#else
#define F_RECL 4
#endif

        implicit none

        ! Parameter Constants
        integer, parameter :: 
     .  file1 = 10, file2 = 20, seed_dim = 34

        ! Local Variables
        integer(kind=4), dimension(seed_dim) :: seed ! Our seed data

        integer :: i ! Loop index

        open(unit=file1, file='rand_seed.bin', action='read', 
     .       recl=F_RECL, form='unformatted', access='direct')

        open(unit=file2, file='rand_seed.dat')

        do i = 1, seed_dim, 1
          read(unit=file1, rec=i) seed(i)
        end do

        do i = 1, seed_dim, 1
          write(unit=file2, fmt=* ) seed(i)
        end do

        close(unit=file1)
        close(unit=file2)

        end program int2txt
