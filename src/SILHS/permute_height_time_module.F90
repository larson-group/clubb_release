!-----------------------------------------------------------------------
!$Id$
!===============================================================================
module permute_height_time_module

  implicit none

  public :: permute_height_time, rand_permute

  private ! Default Scope

  contains
!-----------------------------------------------------------------------

  subroutine permute_height_time( nt_repeat, n_vars, height_time_matrix )

! Description:
!   Generates a matrix height_time_matrix, which is a nt_repeat x n_vars
!   matrix whose 1st dimension is random permutations of the integer sequence 
!   (0,...,nt_repeat-1).

! References:
!   None
!-----------------------------------------------------------------------

    implicit none

    ! Input Variables

    integer, intent(in) :: &
      nt_repeat,  & ! Total number of sample points before sequence repeats.
      n_vars        ! The number of variates in the uniform sample

    ! Output Variables

    integer, dimension(nt_repeat,n_vars), intent(out) :: &
      height_time_matrix ! nt_repeat x n_vars matrix of integers

    ! Local Variables

    integer :: i

    ! Choose elements of height_time_matrix, with a random integer LH sample
    ! for each variate
    do i = 1, n_vars
      call rand_permute( nt_repeat, height_time_matrix(1:nt_repeat,i) )
    end do

    return
  end subroutine permute_height_time
!------------------------------------------------------------------------

!------------------------------------------------------------------------
  subroutine rand_permute( n, pvect )
! Description:
!   Generates a vector of length n
!      containing the integers 0, ... , n-1 in random order.
!   We do not use a new seed.

! References:
!   Follow `Quasi-Monte Carlo sampling' by Art Owen, Section 1.3
!   He follows, in turn, Luc Devroye 'Non-uniform random ...' (1986)
!----------------------------------------------------------------------

    use mt95, only: genrand_real3 ! Procedures

    use mt95, only: genrand_real ! Constants

    implicit none

    ! External

    intrinsic :: int

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

end module permute_height_time_module
