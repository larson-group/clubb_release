!----------------------------------------------------------------------
! subroutine gaus_mixt_points( )

! Generates n random samples 
!     from a d-dimensional Gaussian-mixture PDF.
! Uses Latin hypercube method.

! Input: n = number of calls to microphysics (normally=2)
!        d = number of variates (normally=5) 
!        a = mixture fraction of Gaussians
!        mu1, mu2 = d-dimensional column vector of means 
!                                        of 1st, 2nd Gaussians
!        Sigma1, Sigma2 = dxd dimensional covariance matrix
!        C1, C2 = cloud fraction associated w/ 1st, 2nd 
!                                            mixture component
!        X_u = nxd Latin hypercube sample 
!                                  from uniform distribution
!        s_pts = n-dimensional vector giving values of s 

! Output: [n by d] matrix, X_gm, each row of which is 
!                                     a d-dimensional sample
!----------------------------------------------------------------------
      subroutine gaus_mixt_points( n, d, a, mu1, mu2, Sigma1, Sigma2,
     .                             C1, C2, X_u, s_pts, X_gm )
      
      implicit none

! Input
  
      integer, intent(in)          :: n, d

      double precision, intent(in) :: a, C1, C2
      double precision, intent(in) :: mu1(1:d), mu2(1:d)
      double precision, intent(in) :: Sigma1(1:d,1:d), Sigma2(1:d,1:d)
      double precision, intent(in) :: X_u(1:n,1:(d+1))
      double precision, intent(in) :: s_pts(1:n)

! Output

      double precision, intent(out) :: X_gm(1:n,1:d) 

! Local

      integer j, k, sample
      double precision std_normal(1:d)
      double precision fraction_1
      double precision ltqnorm

! Handle some possible errors re: proper ranges of a, C1, C2.
      if (a > 1.0d0 .or. a < 0.0d0) then
        print *, 'Error in gaus_mixt_points: 
     .             mixture fraction, a, does not lie in [0,1].'
         stop
      endif
      if (C1 > 1.0d0 .or. C1 < 0.0d0) then 
        print *, 'Error in gaus_mixt_points: 
     .           cloud fraction 1, C1, does not lie in [0,1].'
         stop
      endif
      if (C2 > 1.0d0 .or. C2 < 0.0d0) then 
        print *, 'Error in gaus_mixt_points: 
     .           cloud fraction 2, C2, does not lie in [0,1].'
        stop
      endif

! Make sure there is some cloud.
      if (a*C1 < 0.001d0 .and. (1-a)*C2 < 0.001d0) then 
        print *, 'Error in gaus_mixt_points: 
     .              there is none or almost no cloud!'
      endif

      do sample = 1, n

! From Latin hypercube sample, generate standard normal sample
        do j = 1, d
          std_normal(j) = ltqnorm( X_u(sample,j) )
        enddo
    
! Choose which mixture fraction we are in.  
! Account for cloud fraction.
! Follow M. E. Johnson (1987), p. 56.
        fraction_1 = ( a*C1 ) / ( a*C1 + (1-a)*C2 )
        if ( X_u(sample, d+1) < fraction_1 ) then
          call gaus_condt( n, d, std_normal, mu1, Sigma1, s_pts(sample),
     .                     X_gm(sample, 1:d) )    
        else
          call gaus_condt( n, d, std_normal, mu2, Sigma2, s_pts(sample),
     .                     X_gm(sample, 1:d) )   
        endif

! Loop to get new sample
      enddo

      return
      end subroutine gaus_mixt_points
!------------------------------------------------------------------------
