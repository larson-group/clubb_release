!----------------------------------------------------------------------
! Take a binary file of 34 32bit integers and generate text file from it

! Input:  34 bytes of binary data;  "rand_seed.bin"

! Output: ASCII text; "rand_seed.dat"

! Notes:  Basically useless on machines without a /dev/random.
!         Free format fortran 90
!----------------------------------------------------------------------
program int2txt

implicit none

integer i
integer, dimension(34) :: x

open(unit=10, file='rand_seed.bin', action="read", recl=4, form='unformatted', access='direct')
open(unit=15, file='rand_seed.dat')

do i=1, 34
  read(unit=10, rec=i) x(i)
  write(15, *) x(i)
enddo

close(10)
close(15)

end program int2txt
