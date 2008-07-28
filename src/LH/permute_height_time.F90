!$Id: permute_height_time.F90,v 1.2 2008-07-28 19:20:06 faschinj Exp $
module permute_height_time_mod

implicit none

public :: permute_height_time

private :: generate_k_order

private ! Default Scope

contains
!-----------------------------------------------------------------------
! subroutine permute_height_time( )

! Generates a matrix p_height_time, which is a nnzp x nt matrix whose rows
! are random permutations of the integer sequence (0,...,nt-1).
! from 1 to sequence_length. k_order gives vertical ordering 
! of sample points; generate a new k_order every nt/n time steps.  

! Input: nnzp = total number of vertical levels in the model 
!            timestep; first timestep must have i==1)
!        nt = total number of sample points before sequence 
!                          repeats. 
!        dp1 = number of variates + 1

! Output: p_height_time = nnzp x nt x dp1 matrix of integers 
!-----------------------------------------------------------------------


  subroutine permute_height_time( nnzp, nt, dp1,  & 
                                  p_height_time )
  
  use random, only: rand_permute ! Procedure(s)
  
  implicit none

! Input
  
  integer, intent(in) ::  nnzp, nt, dp1

! Output

  integer, intent(out) :: p_height_time(1:nnzp,1:nt,1:dp1)

! Local

  integer i, k

! Choose elements of p_height_time, with a random integer LH sample
! for each altitude and for each variate 
do k = 1, nnzp
  do i = 1, dp1
    call rand_permute(nt, p_height_time(k,1:nt,i))
  enddo
enddo

! Make elements of p_height_time in the range [1,nt] inclusive
!       p_height_time = p_height_time + 1

!       print*, 'p_height_time in permute_height_time=', p_height_time

return
end subroutine permute_height_time
!------------------------------------------------------------------------

!----------------------------------------------------------------------
! subroutine generate_k_order( )

! Generates a vector k_order, which is a random vector of integers
! from 1 to sequence_length. k_order gives vertical ordering 
! of sample points; generate a new k_order every nt/n time steps.  

! Input: i = timestep number (increases by one with each new 
!            timestep; first timestep must have i==1)
!        sequence_length = number of timesteps before sequence 
!                          repeats. 

! Output: k_order = vector of length sequence_length. 
!----------------------------------------------------------------------
  subroutine generate_k_order( i, sequence_length, k_order )

  use random, only: rand_permute ! Procedure(s)
  
  implicit none

! Input
  
  integer, intent(in) :: i, sequence_length

! Output

  integer, intent(out) :: k_order( 1:sequence_length )

! i==1 must be the first timestep 
if (mod(i-1,sequence_length) == 0) then 
  call rand_permute(sequence_length, k_order)
  k_order = k_order + 1
endif

return
end subroutine generate_k_order
!------------------------------------------------------------------------


end module permute_height_time_mod
