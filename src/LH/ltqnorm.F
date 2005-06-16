!-----------------------------------------------------------------------
! This function is ported to Fortran from the same function written in Matlab, see the following
! description of this function.  Hongli Jiang, 2/17/2004
! Converted to double precision by Vince Larson 2/22/2004;
!    this improves results for input values of p near 1.

! LTQNORM Lower tail quantile for standard normal distribution.
!
!   Z = LTQNORM(P) returns the lower tail quantile for the standard normal
!   distribution function.  I.e., it returns the Z satisfying Pr{X < Z} = P,
!   where X has a standard normal distribution.
!
!   LTQNORM(P) is the same as SQRT(2) * ERFINV(2*P-1), but the former returns a
!   more accurate value when P is close to zero.

!   The algorithm uses a minimax approximation by rational functions and the
!   result has a relative error less than 1.15e-9.  A last refinement by
!   Halley's rational method is applied to achieve full machine precision.

!   Author:      Peter J. Acklam
!   Time-stamp:  2003-04-23 08:26:51 +0200
!   E-mail:      pjacklam@online.no
!   URL:         http://home.online.no/~pjacklam
!-----------------------------------------------------------------------
        function ltqnorm( p )

        use constants, only: Pi_DP

        implicit none

        double precision a1, a2, a3, a4, a5, a6, b1, b2, b3, b4, b5,
     .                   c1, c2, c3, c4, c5, c6, d1, d2, d3, d4
        double precision p, q, r, z, z1, plow, phigh, e, u, erf_dp
        double precision ltqnorm
!       Occurs in constants.F now.  Isn't actually used currently.
!        double precision, parameter :: pi=3.1415926d0

! Coefficients in rational approximations.
! equivalent: a(1)=a1, a(2)=a2, and etc, when a(1) is in Matlab. 
! Similarly for b, c, and d's
        parameter (a1 = -3.969683028665376d+01, 
     .             a2 = 2.209460984245205d+02,
     .             a3 = -2.759285104469687d+02, 
     .             a4 = 1.383577518672690d+02,
     .             a5 = -3.066479806614716d+01, 
     .             a6 = 2.506628277459239d+00)
        parameter (b1 = -5.447609879822406d+01, 
     .             b2 = 1.615858368580409d+02,
     .             b3 = -1.556989798598866d+02, 
     .             b4 = 6.680131188771972d+01,
     .             b5 = -1.328068155288572d+01)
        parameter (c1 = -7.784894002430293d-03, 
     .             c2 = -3.223964580411365d-01,
     .             c3 = -2.400758277161838d+00, 
     .             c4 = -2.549732539343734d+00,
     .             c5 =  4.374664141464968d+00, 
     .             c6 =  2.938163982698783d+00)
        parameter (d1 =  7.784695709041462d-03, 
     .             d2 =  3.224671290700398d-01,
     .             d3 =  2.445134137142996d+00, 
     .             d4 =  3.754408661907416d+00)


!  Define break-points.
        plow  = 0.02425d0
        phigh = 1.d0 - plow

!  Initialize output array. Don't need this in Fortran
!   z = zeros(size(p));

!  Rational approximation for lower region:
        if (p > 0.d0 .and. p < plow) then 
           q = sqrt( -2 * log(p) )
           z = (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/
     .               ((((d1*q+d2)*q+d3)*q+d4)*q+1.d0)
!  Rational approximation for central region:
      elseif (p >= plow .and. p <= phigh) then 
         q = p - 0.5d0
         r = q * q
         z = (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q
     .              /(((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1.d0)
! Rational approximation for upper region:
      elseif (p > phigh .and. p < 1.d0) then
         q  = sqrt( -2.d0 * log(1.d0 - p) )
         z  = -(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)
     .               /((((d1*q+d2)*q+d3)*q+d4)*q+1.d0)
      endif

!  Case when P = 0: z = -inf, to create inf z =-1./0., 
!     to create NaN's inf*inf.
        z1 = 0.d0
        if (p == 0.d0) then
           z = (-1.d0)/z1
        endif

! Case when P = 1:, z=inf
        if(p == 1.d0)then
         z = 1.d0/z1
        endif

!  Cases when output will be NaN:
!   k = p < 0 | p > 1 | isnan(p);
! usually inf*inf --> NaN's. 
        if (p < 0.d0 .or. p > 1d0) then
         z = (1.d0/z1)**2
        endif

!  The relative error of the approximation has absolute value less
!  than 1.15e-9.  One iteration of Halley's rational method (third
!  order) gives full machine precision.
! V. Larson 20Feb04: Don't use the following if-endif loop.
!   The value of e is very different than what MATLAB produces, 
!   possibly because of
!   poor values of erf from Numerical Recipes.
!   The value is close to MATLAB's 
!   if I omit the following if-endif loop.
! End V. Larson comment
!!   k = 0 < p & p < 1;
!       if (p.gt.0 .and. p.lt.1)then
!         e = 0.5*(1.0 - erf_dp(-z/sqrt(2.))) - p          ! error
!         u = e * sqrt(2*pi_dp) * exp(z**2/2)       ! f(z)/df(z)
!         z = z - u/( 1 + z*u/2 )               ! Halley's method
!       endif

! return z as double precision:
        ltqnorm = z

        return
        end function ltqnorm
