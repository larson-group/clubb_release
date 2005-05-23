
! Generates a vector k_order, which is a random vector of integers
! from 1 to sequence_length. k_order gives vertical ordering 
! of sample points; generate a new k_order every nt/n time steps.  
! Input: i = timestep number (increases by one with each new 
!            timestep; first timestep must have i==1)
!        sequence_length = number of timesteps before sequence 
!                          repeats. 
! Output: k_order = vector of length sequence_length. 


        subroutine generate_k_order(i,sequence_length,k_order)

        implicit none

! Input
	  
	  integer, intent(in) :: i, sequence_length

! Output

        integer, intent(out) :: k_order(1:sequence_length)


! i==1 must be the first timestep 
	if (mod(i-1,sequence_length) .eq. 0) then 
	   call rand_permute(sequence_length, k_order)
	   k_order = k_order + 1	 
	endif


	return
      end

!------------------------------------------------------------------------
