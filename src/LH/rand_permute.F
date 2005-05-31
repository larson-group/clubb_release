!----------------------------------------------------------------------
! subroutine rand_permute( )

! Generates a vector of length n
!    containing the integers 0, ... , n-1 in random order.
! We do not use a new seed.
! Follow `Quasi-Monte Carlo sampling' by Art Owen, Section 1.3
! He follows, in turn, Luc Devroye 'Non-uniform random ...' (1986)

! Input: n = number of elements to permute

! Output: pvect = array of n numbers in random order
!----------------------------------------------------------------------
        subroutine rand_permute( n, pvect )

        implicit none

! Input
  
        integer, intent(in) :: n

! Output

        integer, intent(out) :: pvect(1:n)

! Local

        integer j, k, temp, seed
        double precision ran2

! Continue the old string of random numbers by choosing seed>0
        seed = 1;

! Start with an ordered vector, pvect
        do j=1,n
          pvect(j) = j
        enddo

! Now re-arrange the elements
        do j=n,2,-1
          temp = pvect(j)
          ! choose an element randomly between 1 and j
          k = int(dble(j)*ran2(seed)+1.0);
          ! swap elements j and k 
          pvect(j) = pvect(k)
          pvect(k) = temp
        enddo

! Convert range of array from 1:n to 0:n-1
        do j=1,n
          pvect(j) = pvect(j) - 1
        enddo

        return
        end subroutine rand_permute
!------------------------------------------------------------------------
