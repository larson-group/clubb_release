! $Id$
module random

  implicit none

  public :: rand_permute

  private ! Default Scope

  contains
!----------------------------------------------------------------------

  subroutine rand_permute( n, pvect )
! Generates a vector of length n
!    containing the integers 0, ... , n-1 in random order.
! We do not use a new seed.
! Follow `Quasi-Monte Carlo sampling' by Art Owen, Section 1.3
! He follows, in turn, Luc Devroye 'Non-uniform random ...' (1986)
!----------------------------------------------------------------------

    use mt95, only: genrand_real3 ! Procedures

    use mt95, only: genrand_real ! Constants

    implicit none

    ! External

    intrinsic :: int, dble

    ! Input Variables

    integer, intent(in) :: n ! Number of elements to permute

    ! Output Variables

    integer, dimension(n), intent(out) :: &
      pvect ! Array of n numbers in random order

    ! Local Variables

    integer j, k, temp

    real(kind=genrand_real) :: rand ! Random float on interval (0,1)

    ! Start with an ordered vector, pvect
    do j=1,n
      pvect(j) = j
    end do

    ! Now re-arrange the elements
    do j=n,2,-1
      temp = pvect(j)
      call genrand_real3( rand ) ! real3 excludes 0 and 1.
      ! choose an element randomly between 1 and j
      k = int( real( j, kind=genrand_real )*rand+1.0_genrand_real )
      ! swap elements j and k
      pvect(j) = pvect(k)
      pvect(k) = temp
    end do

    ! Convert range of array from 1:n to 0:n-1
    do j=1,n
      pvect(j) = pvect(j) - 1
    end do

    return
  end subroutine rand_permute
!------------------------------------------------------------------------
end module random
