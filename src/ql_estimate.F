
! Compute an Monte Carlo estimate of grid box avg liquid water.
! Input: n = number of calls to microphysics (normally=2)
!        d = number of variates (normally=5) 
!        a = mixture fraction of Gaussians
!        C1, C2 = cloud fraction associated 
!                            w/ 1st, 2nd mixture component
!        ql = n in-cloud values of spec liq water content (g/kg).
!        w  = n in-cloud values of vertical velocity (m/s)
!        N_pts  = n in-cloud values of droplet number (#/mg air)
!        rr = n in-cloud values of specific rain content (g/kg) 
!        X_u = nxd Latin hypercube sample from uniform dist 
! Output: ql_m = a scalar representing 
!                      grid box avg specific liquid water;
!                      divide by total cloud fraction to obtain 
!                      within-cloud liquid water

      subroutine ql_estimate(n,d,a,C1,C2,ql,w,N_pts,rr,X_u,ql_m)
      
      implicit none

! Input
	  
	integer, intent(in) :: n, d

	double precision, intent(in) :: a, C1, C2
	double precision, intent(in) :: ql(1:n), w(1:n)
	double precision, intent(in) :: N_pts(1:n), rr(1:n)
	double precision, intent(in) :: X_u(1:n,1:(d+1))

! Output

      double precision, intent(out) :: ql_m 

! Local

	integer sample
	integer n1, n2
	double precision ql_m1, ql_m2
	double precision coeff, expn
	double precision fraction_1


! Handle some possible errors re: proper ranges of a, C1, C2.
	if (a .gt. 1.0d0 .or. a .lt. 0.0d0) then
         print*, 'Error in autoconv_driver: 
     .             mixture fraction, a, does not lie in [0,1].'
         stop
	endif
	if (C1 .gt. 1.0d0 .or. C1 .lt. 0.0d0) then 
         print*, 'Error in autoconv_driver: 
     .           cloud fraction 1, C1, does not lie in [0,1].'
         stop
	endif
	if (C2 .gt. 1.0d0 .or. C2 .lt. 0.0d0) then 
         print*, 'Error in autoconv_driver: 
     .           cloud fraction 2, C2, does not lie in [0,1].'
         stop
	endif

! Make sure there is some cloud.
	if (a*C1 .lt. 0.001d0 .and. (1-a)*C2 .lt. 0.001d0) then 
         print*, 'Error in autoconv_driver: 
     .              there is none or almost no cloud!'
	endif

! To compute liquid water, need to set coeff=expn=1.
	coeff = 1.d0
	expn = 1.d0

! Initialize liquid in each mixture component
	ql_m1 = 0.d0
	ql_m2 = 0.d0

! Initialize numbers of sample points corresponding 
!    to each mixture component
	n1 = 0
	n2 = 0

	do sample = 1,n

! Choose which mixture fraction we are in.  
! Account for cloud fraction.
! Follow M. E. Johnson (1987), p. 56.
         fraction_1 = a*C1/(a*C1+(1-a)*C2)
         if ( X_u(sample,d+1) .lt. fraction_1 ) then
! Use an idealized formula to compute liquid 
!      in mixture comp. 1
             ql_m1 = ql_m1 + coeff*(ql(sample))**expn
             n1 = n1 + 1    
         else
! Use an idealized formula to compute liquid 
!      in mixture comp. 2
             ql_m2 = ql_m2 + coeff*(ql(sample))**expn
             n2 = n2 + 1    
         endif

! Loop to get new sample
	enddo

!! Convert sums to averages.
!! Old code that underestimates if n1 or n2 = 0.
!	if (n1 .eq. 0) then
!	   ql_m1 = 0.d0
!	else
!	   ql_m1 = ql_m1/n1
!	endif

!	if (n2 .eq. 0) then
!         ql_m2 = 0.d0
!	else
!	   ql_m2 = ql_m2/n2
!	endif


! Convert sums to averages.
! If we have no sample points for a certain plume,
!    then we estimate the plume liquid water by the
!    other plume's value.
      if (n1 .eq. 0 .and. n2 .eq. 0) then
	   print*, 'Error: no sample points in ql_estimate'
	   stop
      endif

	if ( .not. (n1 .eq. 0) ) then
	   ql_m1 = ql_m1/n1
	endif

	if ( .not. (n2 .eq. 0) ) then
	   ql_m2 = ql_m2/n2
	endif

	if (n1 .eq. 0) then
	   ql_m1 = ql_m2
	endif	

	if (n2 .eq. 0) then
	   ql_m2 = ql_m1
	endif	

! Grid box average.
	ql_m = a*C1*ql_m1 + (1-a)*C2*ql_m2

	return
      end

!---------------------------------------------------------------
