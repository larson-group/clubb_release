! $Id$
!-------------------------------------------------------------------------------
! Description:
!   Recursive Fortran 95 quicksort routine
!   sorts real numbers into ascending numerical order
!
! Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
!
! References:
!   Based on algorithm from Cormen et al., Introduction to Algorithms,
!   1997 printing
!
! Made F conformant by Walt Brainerd
! 
! Modified by David Schanen to sort a logical table corresponding to the real
! values.
!-------------------------------------------------------------------------------

module quicksort

  implicit none

  public :: Qsort_flags
  private :: Partition

  private

  contains

  recursive subroutine Qsort_flags( table, cost  )

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input/Output Variables
    logical, intent(inout), dimension(:,:) :: table
    real( kind = core_rknd ), intent(inout), dimension(:) :: cost

    ! Local Variables
    integer :: iq

    if ( size( cost ) > 1 ) then
      call Partition( table, cost, iq )
      call Qsort_flags( table(:iq-1,:), cost(:iq-1) )
      call Qsort_flags( table(iq:,:), cost(iq:) )
    end if

    return
  end subroutine Qsort_flags

  subroutine Partition( table, cost, marker )

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Intput/Output Variables
    logical, intent(inout), dimension(:,:) :: table
    real( kind = core_rknd ), intent(inout), dimension(:) :: cost

    ! Output Variables
    integer, intent(out) :: marker

    ! Local Variables
    integer :: i, j
    real( kind = core_rknd ) :: temp
    logical, dimension(size( table(1,:) )) :: logic_temp
    real( kind = core_rknd ) :: x  ! pivot point

    ! ---- Begin Code ----

    x = cost(1)
    i = 0
    j = size( cost ) + 1

    marker = -999 ! Initialize to avoid warnings

    do
      j = j-1
      do
        if (cost(j) <= x) exit
        j = j-1
      end do
      i = i+1
      do
        if (cost(i) >= x) exit
        i = i+1
      end do
      if ( i < j ) then
        ! exchange cost(i) and cost(j)
        temp = cost(i)
        logic_temp = table(i,:)
        cost(i) = cost(j)
        table(i,:) = table(j,:)
        cost(j) = temp
        table(j,:) = logic_temp
      else if (i == j) then
        marker = i+1
        return
      else
        marker = i
        return
      end if
    end do

  end subroutine Partition

end module quicksort


