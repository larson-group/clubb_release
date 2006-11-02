!-----------------------------------------------------------------------
! $Id: int2txt.F,v 1.2 2006-11-02 01:36:10 dschanen Exp $
        program int2txt

!       Description:
!       Take a binary file of 34 4 byte integers and generate an ASCII 
!       text file from it.
!       Takes Input of  34 bytes of binary data;  "rand_seed.bin".
!       Outputs ASCII text; "rand_seed.dat".
!       Notes:  Not very useful on machines without a /dev/random, 
!       though any 34 bytes of pseudo-random seed data will suffice.
!-----------------------------------------------------------------------
#include "recl.F"

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
