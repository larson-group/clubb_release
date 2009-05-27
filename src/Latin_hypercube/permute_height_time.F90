!$Id$
module permute_height_time_mod

  implicit none

  public :: permute_height_time

  private :: generate_k_order

  private ! Default Scope

  contains
!-----------------------------------------------------------------------

  subroutine permute_height_time( nnzp, nt_repeat, dp1,  & 
                                  height_time_matrix )

! Description:
!   Generates a matrix height_time_matrix, which is a nnzp x nt matrix whose rows
!   are random permutations of the integer sequence (0,...,nt-1).
!   from 1 to sequence_length. k_order gives vertical ordering
!   of sample points; generate a new k_order every nt/n time steps.
!   First timestep must have i == 1

! References:
!   None
!-----------------------------------------------------------------------


    use random, only: rand_permute ! Procedure(s)

    implicit none

    ! Input Variables

    integer, intent(in) :: &
      nnzp,       & ! Total number of vertical levels in the model timestep. 
      nt_repeat,  & ! Total number of sample points before sequence repeats.
      dp1           ! The number of variates + 1

    ! Output Variables

    integer, dimension(nnzp,nt_repeat,dp1), intent(out) :: &
      height_time_matrix ! nnzp x nt x dp1 matrix of integers

    ! Local Variables

    integer :: i, k

    ! Choose elements of height_time_matrix, with a random integer LH sample
    ! for each altitude and for each variate
    do k = 1, nnzp
      do i = 1, dp1
        call rand_permute( nt_repeat, height_time_matrix(k,1:nt_repeat,i) )
      end do
    end do

    ! Make elements of height_time_matrix in the range [1,nt] inclusive
    !height_time_matrix = height_time_matrix + 1

    !print*, 'height_time_matrix in permute_height_time=', height_time_matrix

    return
  end subroutine permute_height_time
!------------------------------------------------------------------------

!----------------------------------------------------------------------
  subroutine generate_k_order( i, sequence_length, k_order )

! Description:
!   Generates a vector k_order, which is a random vector of integers
!   from 1 to sequence_length. k_order gives vertical ordering
!   of sample points; generate a new k_order every nt/n time steps.
! References:
!   None
!----------------------------------------------------------------------

    use random, only: rand_permute ! Procedure(s)

    implicit none

    ! External

    intrinsic :: mod

    ! Input Variables

    ! Timestep number (increases by one with each new timestep; first timestep must have i==1)
    integer, intent(in) :: i

    integer, intent(in) :: &
      sequence_length ! Number of timesteps before sequence repeats.

    ! Output Variables

    integer, dimension(sequence_length), intent(out) :: &
      k_order ! Vector of length sequence_length.

    ! i==1 must be the first timestep
    if ( mod( i-1, sequence_length ) == 0 ) then
      call rand_permute( sequence_length, k_order )
      k_order = k_order + 1
    end if

    return
  end subroutine generate_k_order
!------------------------------------------------------------------------


end module permute_height_time_mod
