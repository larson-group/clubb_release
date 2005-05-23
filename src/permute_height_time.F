
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


        subroutine permute_height_time(nnzp,nt,dp1,p_height_time)

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
!	p_height_time = p_height_time + 1

!	print*, 'p_height_time in permute_height_time=', p_height_time

	return
      end

!------------------------------------------------------------------------
