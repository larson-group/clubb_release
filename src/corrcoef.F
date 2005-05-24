
! Various mathematical utilities

	function corrcoef(vect1,vect2,n)
! Correlation coefficient of two vectors

	implicit none

	integer n
	double precision vect1(1:n), vect2(1:n)

	double precision corrcoef

	double precision cov

	corrcoef = cov(vect1,vect2,n) /
     .         sqrt( cov(vect1,vect1,n) * cov(vect2,vect2,n) )

	return
	end	

! -----------------------------------------------------------


	function std(vector,n)

	implicit none

	integer n
	double precision std, vector(1:n)

	double precision cov

	std = sqrt( cov(vector,vector,n) * (n/(n-1)) )

	return
	end	


! ---------------------------------------------------------	


	function cov(vect1,vect2,n)
! Covariance of two vectors

	implicit none

	integer n
	double precision cov, vect1(1:n), vect2(1:n)


	double precision sum, avg1, avg2, mean
	integer j

	avg1 = mean(vect1,n)
	avg2 = mean(vect2,n)
	
	sum = 0.d0
	do j=1,n
	   sum = sum + (vect1(j) - avg1) * (vect2(j) - avg2)
	enddo

	cov = sum/n

	return
	end	

! ---------------------------------------------------------	

	function mean(vector,n)

	implicit none

	integer n
	double precision mean, vector(1:n)

	double precision sum
	integer j
	
	sum = 0.d0
	do j=1,n
	   sum = sum + vector(j)
	enddo

	mean = sum/n

	return
	end
