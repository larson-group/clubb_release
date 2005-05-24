
! Generates a matrix X that contains a Latin Hypercube sample.
! The sample is uniformly distributed.  
! See Art B. Owen (2003), ``Quasi-Monte Carlo Sampling," 
!    a chapter from SIGGRAPH 2003
! Input: n = number of calls to microphysics (normally=2)
!        dp1 = number of variates plus 1 (normally=6) 
!        p_matrix = n x dp1 array of permuted integers
! Output: n by dp1 matrix, X, 
!    each row of which is a dp1-dimensional sample


        subroutine latin_hyper_sample(n,nt,dp1,p_matrix,X)

        implicit none

! Input
	  
	  integer, intent(in) :: n, nt, dp1, p_matrix(1:n,1:dp1)


! Output

        double precision, intent(out) :: X(1:n,1:dp1)

! Local

	  integer j, k, seed
!	  integer p_matrix(1:n,1:dp1)
	  double precision ran2

! Continue the old string of random numbers by choosing seed>0
	  seed = 1
!	  seed = 2

!! Compute random permutation row by row
!        do j=1,dp1
!          ! Generate a column vector of integers from 0 to n-1, 
!	    !    whose order is random.
!	     call rand_permute( n, p_matrix(1:n,j) )
!        enddo

! Choose values of sample 
!     using permuted vector and random number generator.
	  do j = 1,n
	     do k = 1,dp1

              X(j,k) = (1.0d0/nt)*(p_matrix(j,k) + ran2(seed));

           enddo
	  enddo

!	  print*, 'p_matrix(:,1)= ', p_matrix(:,1)
!        print*, 'p_matrix(:,dp1)= ', p_matrix(:,dp1)
!	  print*, 'X(:,1)= ', X(:,:)

	  return
        end

!------------------------------------------------------------------------
