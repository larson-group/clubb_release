
! Converts sample points drawn from a uniform distribution 
!    to truncated Gaussian points.
! Input: n = number of calls to microphysics (normally=2)
!        d = number of variates (normally=5) 
!        col = scalar indicating which column of X_nl to truncate
!        a = mixture fraction of Gaussians
!        mu1, mu2 = d-dimensional column vector 
!                         of means of 1st, 2nd Gaussians
!        Sigma1, Sigma2 = dxd dimensional covariance matrix
!        C1, C2 = cloud fraction associated w/ 1st, 2nd mixture 
!                                                      component
!        X_u = nxd Latin hypercube sample 
!                            from uniform distribution 
! Output: truncated_column = A column vector of length n 
! that is transformed from a Gaussian PDF to truncated Gaussian PDF.


      subroutine truncate_gaus_mixt(n,d,col,a,mu1,mu2,
     .                Sigma1,Sigma2,C1,C2,X_u,truncated_column)
      
      implicit none

! Input
	  
	integer, intent(in) :: n, d, col

	double precision, intent(in) :: a, C1, C2
	double precision, intent(in) :: mu1(1:d), mu2(1:d)
	double precision, intent(in) :: Sigma1(1:d,1:d), Sigma2(1:d,1:d)
	double precision, intent(in) :: X_u(1:n,1:(d+1))

! Output

      double precision, intent(out) :: truncated_column(1:n) 

! Local

	integer sample
	double precision s_std
	double precision fraction_1
	double precision ltqnorm


! Handle some possible errors re: proper ranges of a, C1, C2.
	if ( (a .gt. 1.0d0) .or. (a .lt. 0.0d0) ) then
         print*, 'Error in truncate_gaus_mixt: 
     .             mixture fraction, a, does not lie in [0,1].'
         stop
	endif
	if ( (C1 .gt. 1.0d0) .or. (C1 .lt. 0.0d0) ) then 
         print*, 'Error in truncate_gaus_mixt: 
     .           cloud fraction 1, C1, does not lie in [0,1].'
         stop
	endif
	if ( (C2 .gt. 1.0d0) .or. (C2 .lt. 0.0d0) ) then 
         print*, 'Error in truncate_gaus_mixt: 
     .           cloud fraction 2, C2, does not lie in [0,1].'
         stop
	endif

! Make sure there is some cloud.
	if (a*C1 .lt. 0.001d0 .and. (1-a)*C2 .lt. 0.001d0) then 
         print*, 'Error in truncate_gaus_mixt: 
     .              there is none or almost no cloud!'
	endif

! Make s PDF (1st column) a truncated Gaussian.
! This allows us to sample solely from the cloud points.
	do sample = 1,n

! Choose which mixture fraction we are in.  
! Account for cloud fraction.
! Follow M. E. Johnson (1987), p. 56.
         fraction_1 = a*C1/(a*C1+(1.d0-a)*C2)
         if ( X_u(sample,d+1) .lt. fraction_1 ) then
! Replace first dimension (s) with 
!   sample from cloud (i.e. truncated standard Gaussian)
           s_std = ltqnorm( X_u(sample,col)*C1 + (1.d0-C1) ) 
! Convert to nonstandard normal with mean mu1 and variance Sigma1
           truncated_column(sample) = 
     .               s_std * sqrt( Sigma1(col,col) ) + mu1(col)    
         else
! Replace first dimension (s) with 
!   sample from cloud (i.e. truncated Gaussian)
           s_std = ltqnorm( X_u(sample,col)*C2 + (1.d0-C2) ) 
! Convert to nonstandard normal with mean mu2 and variance Sigma2
           truncated_column(sample) = 
     .               s_std * sqrt( Sigma2(col,col) ) + mu2(col) 
         endif

! Loop to get new sample
	enddo

	return
      end

!------------------------------------------------------------------------

