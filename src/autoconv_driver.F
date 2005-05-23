
! Compute Kessler grid box avg autoconversion (g/kg)/s.
! Input: n = number of calls to microphysics (normally=2)
!        d = number of variates (normally=5) 
!        a = mixture fraction of Gaussians
!        R1, R2 = cloud fraction associated 
!                            w/ 1st, 2nd mixture component
!        ql = n in-cloud values of spec liq water content (g/kg).
!        w  = n in-cloud values of vertical velocity (m/s)
!        N_pts  = n in-cloud values of droplet number (#/mg air)
!        rr = n in-cloud values of specific rain content (g/kg) 
!        X_u = nxd Latin hypercube sample from uniform dist 
! Output: ac_m = a scalar representing 
!                      grid box average autoconversion;
!                      has same units as ql/s;
!                      divide by total cloud fraction to obtain 
!                      within-cloud autoconversion

      subroutine autoconv_driver(n,d,a,R1,R2,ql,w,N_pts,rr,X_u,ac_m)
      
      implicit none

! Input
	  
	integer, intent(in) :: n, d

	double precision, intent(in) :: a, R1, R2
	double precision, intent(in) :: ql(1:n), w(1:n)
	double precision, intent(in) :: N_pts(1:n), rr(1:n)
	double precision, intent(in) :: X_u(1:n,1:(d+1))

! Output

      double precision, intent(out) :: ac_m 

! Local

	integer sample
	integer n1, n2
	double precision ac_m1, ac_m2
	double precision coeff, expn, q_crit
	double precision fraction_1


! Handle some possible errors re: proper ranges of a, R1, R2.
	if (a .gt. 1.0d0 .or. a .lt. 0.0d0) then
         print*, 'Error in autoconv_driver: 
     .             mixture fraction, a, does not lie in [0,1].'
	   print*, 'a=', a
         stop
	endif
	if (R1 .gt. 1.0d0 .or. R1 .lt. 0.0d0) then 
         print*, 'Error in autoconv_driver: 
     .           cloud fraction 1, R1, does not lie in [0,1].'
	   print*, 'R1=', R1
         stop
	endif
	if (R2 .gt. 1.0d0 .or. R2 .lt. 0.0d0) then 
         print*, 'Error in autoconv_driver: 
     .           cloud fraction 2, R2, does not lie in [0,1].'
         print*, 'R2=', R2
         stop
	endif

! Make sure there is some cloud.
	if (a*R1 .lt. 0.001d0 .and. (1-a)*R2 .lt. 0.001d0) then 
         print*, 'Error in autoconv_driver: 
     .              there is none or almost no cloud!'
	endif

! Autoconversion formula prefactor and exponent.
! These are for Kessler autoconversion in (g/kg)/s.
	coeff = 1.d-3
	expn = 1.d0
!        q_crit = 0.3d0
!	q_crit = 0.7d0
        q_crit = 0.2d0

! Initialize autoconversion in each mixture component
	ac_m1 = 0.d0
	ac_m2 = 0.d0

! Initialize numbers of sample points corresponding 
!    to each mixture component
	n1 = 0
	n2 = 0

	do sample = 1,n

! Choose which mixture fraction we are in.  
! Account for cloud fraction.
! Follow M. E. Johnson (1987), p. 56.
         fraction_1 = a*R1/(a*R1+(1-a)*R2)
!	   print*, 'fraction_1= ', fraction_1 

! V. Larson change to try to fix sampling
!         if ( X_u(sample,d+1) .lt. fraction_1 ) then
!	   print*, '-1+2*int((sample+1)/2)= ', -1+2*int((sample+1)/2)
!	   print*, '-1+2*int((sample+1)/2)= ', int(sample)
         if ( X_u(1,d+1) .lt. fraction_1 ) then
! End of V. Larson fix

! Use an idealized formula to compute autoconversion 
!      in mixture comp. 1
! A_K = (1e-3/s)*(ql-0.5g/kg)*H(ql-0.5g/kg) 
! This is the first of two lines where 
!      a user must add a new microphysics scheme.
             ac_m1 = ac_m1 + coeff*max(0.d0,ql(sample)-q_crit)
             n1 = n1 + 1    
         else
! Use an idealized formula to compute autoconversion 
!      in mixture comp. 2
! A_K = (1e-3/s)*(ql-0.5g/kg)*H(ql-0.5g/kg) 
! This is the second and last line where 
!      a user must add a new microphysics scheme.
             ac_m2 = ac_m2 + coeff*max(0.d0,ql(sample)-q_crit)
             n2 = n2 + 1    
         endif

! Loop to get new sample
	enddo

!! Convert sums to averages.
!! Old code that underestimates if a plume has no sample points
!	if (n1 .eq. 0) then
!	   ac_m1 = 0.d0
!	else
!	   ac_m1 = ac_m1/n1
!	endif

!	if (n2 .eq. 0) then
!         ac_m2 = 0.d0
!	else
!	   ac_m2 = ac_m2/n2
!	endif

! Convert sums to averages.
! If we have no sample points for a certain plume,
!    then we estimate the plume liquid water by the
!    other plume's value.
      if (n1 .eq. 0 .and. n2 .eq. 0) then
	   print*, 'Error: no sample points in autoconv_driver'
	   stop
      endif

	if ( .not. (n1 .eq. 0) ) then
	   ac_m1 = ac_m1/n1
	endif

	if ( .not. (n2 .eq. 0) ) then
	   ac_m2 = ac_m2/n2
	endif

	if (n1 .eq. 0) then
	   ac_m1 = ac_m2
	endif	

	if (n2 .eq. 0) then
	   ac_m2 = ac_m1
	endif	

! Grid box average.
	ac_m = a*R1*ac_m1 + (1-a)*R2*ac_m2

	return
      end

!---------------------------------------------------------------
