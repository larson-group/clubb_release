        module math_utilities         
!-----------------------------------------------------------------------
! Various mathematical utilities
!-----------------------------------------------------------------------
        implicit none
!-----------------------------------------------------------------------
!  function corrcoef( )

! Correlation coefficient of two vectors
!-----------------------------------------------------------------------

        public :: corrcoef, std, cov, mean

        private

        contains
        
        function corrcoef( vect1, vect2, n )

        implicit none

        ! External
        double precision, external :: cov

        ! Input
        integer, intent(in)          :: n
        double precision, intent(in) :: vect1(1:n), vect2(1:n)

        ! Return type
        double precision corrcoef

        corrcoef = cov( vect1, vect2, n ) /
     .         sqrt( cov( vect1, vect1, n ) * cov( vect2, vect2, n ) )

        return
        end function corrcoef

!-----------------------------------------------------------------------
!  function std( )
!-----------------------------------------------------------------------
        function std( vector, n ) 

        implicit none

        ! External
        double precision, external :: cov

        ! Input
        integer, intent(in)          :: n
        double precision, intent(in) :: vector(1:n)

        ! Return type
        double precision std

        std = sqrt( cov( vector, vector, n )*( n/(n-1) ) )

        return
        end function std

!-----------------------------------------------------------------------
!  function cov( )

! Covariance of two vectors
!-----------------------------------------------------------------------
        function cov( vect1, vect2, n )

        implicit none

!       External
        double precision, external   :: mean

!       Input
        integer, intent(in)          :: n
        double precision, intent(in) :: vect1(1:n), vect2(1:n)

!       Return type
        double precision cov

!       Internal
        double precision sum, avg1, avg2 
        integer j

        avg1 = mean( vect1, n )
        avg2 = mean( vect2, n )
        
        sum = 0.d0
        do j = 1, n
           sum = sum + (vect1(j) - avg1) * (vect2(j) - avg2)
        enddo

        cov = sum / n

        return
        end function cov

!-----------------------------------------------------------------------
!  function mean( )
!-----------------------------------------------------------------------
        function mean( vector, n )
        
        implicit none

!       Input 
        integer, intent(in)          :: n
        double precision, intent(in) :: vector(1:n)

!       Internal
        double precision sum
        integer j

!       Return type
        double precision mean
        
        sum = 0.d0
        do j=1, n
           sum = sum + vector(j)
        enddo

        mean = sum/n

        return
        end function mean
!-----------------------------------------------------------------------
        end module math_utilities 
