!$Id$

module micro_calcs_mod

implicit none

public :: micro_calcs

private :: autoconv_driver, ql_estimate

private ! Default Scope

contains
!------------------------------------------------------------------------
! subroutine micro_calcs( )

! This subroutine computes microphysical grid box averages,
!    given a Latin Hypercube sample.  

! Input: n = number of calls to microphysics (normally=2)
!        d = number of variates (normally=5) 
!        X_u = nxd Latin hypercube sample from uniform dist 
!        X_nl = Sample from normal-lognormal distribution 
!        l_sample_flag = logical that tells whether PDF has non-zero micro
!        pdf_parms = pdf parameters output by closure_new

! Output: AKm_est = LH estimate of grid box avg Kessler autoconv [kg/kg]
!         AKm = exact grid box avg Kessler autoconversion [kg/kg]
!         AKstd = exact standard deviation of gba Kessler [kg/kg]
!         AKstd_cld = exact w/in cloud std of gba Kessler [kg/kg]
!         AKm_rcm = exact local gba Kessler auto based on rcm [kg/kg]
!         AKm_rcc = exact local gba Kessler based on w/in cloud rc [kg/kg]
!         rcm_est = LH estimate of grid box avg liquid water [kg/kg]
!------------------------------------------------------------------------

subroutine micro_calcs( n, d, X_u, X_nl, l_sample_flag, & 
                        pdf_parms, & 
                        AKm_est_k, AKm_k, AKstd_k, AKstd_cld_k, & 
                        AKm_rcm_k, AKm_rcc_k, rcm_est_k )

use constants, only: pi, zero_threshold ! Variable(s)

use anl_erf, only: erf ! Procedure(s)

implicit none

! Input

integer, intent(in) :: n, d

! Sample drawn from uniform distribution
double precision, intent(in) :: X_u(1:n,1:(d+1)) 

! Sample that is transformed ultimately to normal-lognormal 
double precision, intent(in) :: X_nl(1:n,1:d)

! Flag that determines whether we have a special case (false)
logical, intent(in) :: l_sample_flag

! PDF parameter array
real, intent(in)    :: pdf_parms(26)

! Output

! Monte Carlo estimate of Kessler autoconversion for kth vertical level. 
real, intent(out) :: AKm_est_k

! Exact Kessler autoconversion, AKm, for kth vertical level.
real, intent(out) :: AKm_k

! Exact std of Kessler autoconversion, AKstd, for kth level,
!      and exact within cloud standard deviation, AKstd_cld
real, intent(out) :: AKstd_k, AKstd_cld_k

! kth Kessler autoconversion, using grid box avg liquid, rcm, as input
real, intent(out) :: AKm_rcm_k

! kth Kessler ac, using within cloud liquid, rcm/cf, as input
real, intent(out) :: AKm_rcc_k

! For comparison, estimate kth liquid water using Monte Carlo
real, intent(out) :: rcm_est_k

! Internal 


! PDF parameters
real :: a
!      real :: w1, w2 
!      real :: sw1, sw2
!      real :: thl1, thl2, sthl1, sthl2
!      real :: rt1,rt2
!      real :: srt1, srt2
real :: ss1, ss2, s1, s2
real :: R1, R2
real :: rc1, rc2

! Cloud fraction 0<cf<1, mean liquid water mix ratio [kg/kg]
real :: cf, rcm

! Double precision version of Monte Carlo Kessler ac
double precision :: AKm_est_dp
! Double precision version of Monte Carlo avg liquid water
double precision :: rcm_est_dp

! Variables needed for exact Kessler autoconversion, AKm
real q_crit, K_one
real sn1_crit, R1_crit, sn2_crit, R2_crit
real AK1, AK2 

! Variables needed for exact std of Kessler autoconversion, AKstd
!      and within cloud standard deviation, AKstd_cld
real AK1var, AK2var 

! For comparison, compute within-cloud vertical velocity analytically.
!      real C_w_cld1, C_w_cld2, w_cld_avg

! Code begins -------------------------------------------


!       save pdf parameters

!      w1    = pdf_parms(1)
!      w2    = pdf_parms(2)
!      sw1   = pdf_parms(3)
!      sw2   = pdf_parms(4)
!      rt1   = pdf_parms(5)
!      rt2   = pdf_parms(6)
!      srt1  = pdf_parms(7)
!      srt2  = pdf_parms(8)
!      thl1  = pdf_parms(9)
!      thl2  = pdf_parms(10)
!      sthl1 = pdf_parms(11)
!      sthl2 = pdf_parms(12)
a     = pdf_parms(13)
rc1   = pdf_parms(14) 
rc2   = pdf_parms(15) 
R1    = pdf_parms(18)
R2    = pdf_parms(19)
s1    = pdf_parms(20)
s2    = pdf_parms(21)
ss1   = pdf_parms(22)
ss2   = pdf_parms(23)

!       compute mean cloud fraction and cloud water

cf    = a * R1 + (1-a) * R2
rcm   = a * rc1 + (1-a) * rc2

!------------------------------------------------------------------------
! Call Kessler autoconversion microphysics using Latin hypercube sample 
! This acts as an interface between the boundary layer scheme
!   and the microphysics.  To add a call to a microphysics scheme,
!   alter two lines in autoconversion_driver.f.
! Then we compute Kessler ac analytically.
!------------------------------------------------------------------------

! Use units of [g/kg] to ameliorate numerical roundoff errors.
! We prognose rt-thl-w,
!    but we set means, covariance of N, qr to constants.

if ( .not. l_sample_flag ) then  

! In this case, sample points could not be constructed.
! Set autoconversion to zero.
    AKm_est_k   = 0.0
    AKm_k       = 0.0
    AKm_rcm_k   = 0.0
    AKm_rcc_k   = 0.0
    rcm_est_k   = 0.0
    AKstd_k     = 0.0
    AKstd_cld_k = 0.0

else

! Call microphysics, i.e. Kessler autoconversion.
! A_K = (1e-3/s)*(ql-0.5g/kg)*H(ql-0.5g/kg) 
  call autoconv_driver(n,d,dble(a),dble(R1),dble(R2),X_nl(1:n,1), & 
         !X_nl(1:n,3),X_nl(1:n,4),X_nl(1:n,5),
         X_u,AKm_est_dp)
! Convert to real number
  AKm_est_k = real(AKm_est_dp)

! Compute Monte Carlo estimate of liquid for test purposes.
  call ql_estimate(n,d,dble(a),dble(R1),dble(R2),X_nl(1:n,1), & 
         !X_nl(1:n,3),X_nl(1:n,4),X_nl(1:n,5),
         X_u,rcm_est_dp)
! Convert to real number
  rcm_est_k = real(rcm_est_dp)

! Convert rcm_est back to (kg/kg) and AKm_est back to (kg/kg)/s.
  rcm_est_k = 1.0e-3*rcm_est_k
  AKm_est_k = 1.0e-3*AKm_est_k

! A test of the scheme:
! Compare exact (rcm) and Monte Carlo estimates (rcm_est) of 
!    specific liq water content, rcm.  
!      print*, 'rcm=', rcm
!      print*, 'rcm_est=', rcm_est

! Exact Kessler autoconversion in units of (kg/kg)/s
!        q_crit = 0.3e-3
!        q_crit = 0.7e-3
  q_crit   = 0.2e-3
  K_one    = 1.e-3
  sn1_crit = (s1-q_crit)/ss1
  R1_crit  = 0.5*(1+erf(sn1_crit/sqrt(2.0)))
  AK1      = K_one * ( (s1-q_crit)*R1_crit  & 
            + ss1*exp(-0.5*sn1_crit**2)/(sqrt(2*pi)) )
  sn2_crit = (s2-q_crit)/ss2
  R2_crit  = 0.5*(1+erf(sn2_crit/sqrt(2.0)))
  AK2      = K_one * ( (s2-q_crit)*R2_crit  & 
            + ss2*exp(-0.5*sn2_crit**2)/(sqrt(2*pi)) )
  AKm_k    = a * AK1 + (1-a) * AK2

! Exact Kessler standard deviation in units of (kg/kg)/s
! For some reason, sometimes AK1var, AK2var are negative
  AK1var   = max( zero_threshold, K_one * (s1-q_crit) * AK1  & 
           + K_one * K_one * (ss1**2) * R1_crit  & 
           - AK1**2  )
  AK2var   = max( zero_threshold, K_one * (s2-q_crit) * AK2  & 
           + K_one * K_one * (ss2**2) * R2_crit  & 
           - AK2**2  )  
! This formula is for a grid box average:
  AKstd_k  = sqrt(    a  * ( (AK1-AKm_k)**2 + AK1var ) & 
              + (1-a) * ( (AK2-AKm_k)**2 + AK2var ) & 
              )
! This formula is for a within-cloud average:
  AKstd_cld_k = sqrt( max( zero_threshold,   & 
            (1./cf) * ( a  * ( AK1**2 + AK1var ) & 
                      + (1-a) * ( AK2**2 + AK2var )  & 
                      ) & 
           - (AKm_k/cf)**2  ) & 
                  )

! Kessler autoconversion, using grid box avg liquid, rcm, as input
  AKm_rcm_k = K_one * max( zero_threshold, rcm-q_crit )

! Kessler ac, using within cloud liquid, rcm/cf, as input
  AKm_rcc_k = cf * K_one * max( zero_threshold, rcm/cf-q_crit ) 

!       print*, 'a=', a
!       print*, 's1=', s1
!       print*, 's2=', s2
!       print*, 'ss1=', ss1
!       print*, 'ss2=', ss2
!       print*, 'AKm_k =', AKm_k
!       print*, 'AKm_est_k (estimate) =', AKm_est_k
!       print*, 'AK1=', AK1
!       print*, 'AK2=', AK2
!       print*, 'AK1var=', AK1var
!       print*, 'AK2var=', AK2var  
!       print*, 'AKstd =', AKstd
!       print*, 'AKstd_cld =', AKstd_cld
!       print*, 'AKm_rcc =', AKm_rcc
!       print*, 'AKm_rcm =', AKm_rcm

!------------------------------------------------------------------------
! Compute within-cloud vertical velocity, avg over full layer
! This call is a kludge: I feed w values into ql variable
! in autoconv_driver.
! Only works if coeff=expn=1 in autoconversion_driver.
!------------------------------------------------------------------------
!     call autoconv_driver( n, d, a, R1, R2, X_nl( 1:n, 3 ),
!     .                     X_nl( 1:n, 3 ), X_nl( 1:n, 4 ), 
!     .                     X_nl(1:n, 5), X_u, AKm2 )

! Another test:
! Compute within-cloud vertical velocity, avgd over full domain.
!        C_w_cld1 =  R1*w1
!        C_w_cld2 =  R2*w2
!        w_cld_avg = a * C_w_cld1 + (1-a) * C_w_cld2

! The following two values should match
!       print*, 'w_cld_avg=', w_cld_avg
!       print*, 'ac_m2=', ac_m2

! End of overall if-then statement for Latin hypercube code
  endif

return
end subroutine micro_calcs
!-----------------------------------------------------------------------
! subroutine autoconv_driver( )

! Compute Kessler grid box avg autoconversion (g/kg)/s.

! Input: n = number of calls to microphysics (normally=2)
!        d = number of variates (normally=5) 
!        a = mixture fraction of Gaussians
!        R1, R2 = cloud fraction associated 
!                            w/ 1st, 2nd mixture component
!        ql = n in-cloud values of spec liq water content (g/kg).
!        w  = n in-cloud values of vertical velocity (m/s)
!        Nc  = n in-cloud values of droplet number (#/mg air)
!        rr = n in-cloud values of specific rain content (g/kg) 
!        X_u = nxd Latin hypercube sample from uniform dist 

! Output: ac_m = a scalar representing 
!                      grid box average autoconversion;
!                      has same units as ql/s;
!                      divide by total cloud fraction to obtain 
!                      within-cloud autoconversion
!-----------------------------------------------------------------------
subroutine autoconv_driver( n, d, a, R1, R2, ql, & !w,  & 
!                            Nc, rr, 
                            X_u, ac_m )

implicit none

! Input
  
  integer, intent(in) :: n, d

  double precision, intent(in) :: a, R1, R2
  double precision, intent(in) :: ql(1:n)
!  double precision, intent(in) :: w(1:n)
!  double precision, intent(in) :: Nc(1:n), rr(1:n)
  double precision, intent(in) :: X_u(1:n,1:(d+1))

! Output

double precision, intent(out) :: ac_m 

! Local

  integer sample
  integer n1, n2
  double precision ac_m1, ac_m2
  double precision coeff, q_crit
!        double precision expn
  double precision fraction_1


! Handle some possible errors re: proper ranges of a, R1, R2.
  if (a .gt. 1.0d0 .or. a .lt. 0.0d0) then
     print*, 'Error in autoconv_driver: ', & 
             'mixture fraction, a, does not lie in [0,1].'
     print*, 'a=', a
     stop
  endif
  if (R1 .gt. 1.0d0 .or. R1 .lt. 0.0d0) then 
     print*, 'Error in autoconv_driver: ', & 
           'cloud fraction 1, R1, does not lie in [0,1].'
     print*, 'R1=', R1
     stop
  endif
  if (R2 .gt. 1.0d0 .or. R2 .lt. 0.0d0) then 
     print*, 'Error in autoconv_driver: ', & 
           'cloud fraction 2, R2, does not lie in [0,1].'
     print*, 'R2=', R2
     stop
  endif

! Make sure there is some cloud.
  if (a*R1 .lt. 0.001d0 .and. (1-a)*R2 .lt. 0.001d0) then 
     print*, 'Error in autoconv_driver: ', & 
              'there is none or almost no cloud!'
  endif

! Autoconversion formula prefactor and exponent.
! These are for Kessler autoconversion in (g/kg)/s.
  coeff  = 1.d-3
!        expn   = 1.d0
!       q_crit = 0.3d0
!       q_crit = 0.7d0
  q_crit = 0.2d0

! Initialize autoconversion in each mixture component
  ac_m1 = 0.d0
  ac_m2 = 0.d0

! Initialize numbers of sample points corresponding 
!    to each mixture component
  n1 = 0
  n2 = 0

  do sample = 1, n

! Choose which mixture fraction we are in.  
! Account for cloud fraction.
! Follow M. E. Johnson (1987), p. 56.
     fraction_1 = a*R1/(a*R1+(1-a)*R2)
!          print*, 'fraction_1= ', fraction_1 

! V. Larson change to try to fix sampling
!          if ( X_u(sample,d+1) .lt. fraction_1 ) then
!          print*, '-1+2*int((sample+1)/2)= ', -1+2*int((sample+1)/2)
!          print*, '-1+2*int((sample+1)/2)= ', int(sample)
     if ( X_u(1,d+1) < fraction_1 ) then
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
!       if (n1 .eq. 0) then
!          ac_m1 = 0.d0
!       else
!          ac_m1 = ac_m1/n1
!       endif

!       if (n2 .eq. 0) then
!         ac_m2 = 0.d0
!       else
!          ac_m2 = ac_m2/n2
!       endif

! Convert sums to averages.
! If we have no sample points for a certain plume,
!    then we estimate the plume liquid water by the
!    other plume's value.
if (n1 == 0 .and. n2 == 0) then
  print*, 'Error: no sample points in autoconv_driver'
  stop
endif

 if ( .not. (n1 == 0) ) then
    ac_m1 = ac_m1/n1
 endif

 if ( .not. (n2 == 0) ) then
    ac_m2 = ac_m2/n2
 endif

 if (n1 == 0) then
    ac_m1 = ac_m2
 endif

 if (n2 == 0) then
    ac_m2 = ac_m1
 endif

! Grid box average.
 ac_m = a*R1*ac_m1 + (1-a)*R2*ac_m2

!        print*, 'autoconv_driver.F: acm=', ac_m

return
end subroutine autoconv_driver
!---------------------------------------------------------------
!----------------------------------------------------------------------
! subroutine ql_estimate( )

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
!----------------------------------------------------------------------
subroutine ql_estimate( n, d, a, C1, C2, ql, & ! w,   & 
                         !N_pts, rr, 
                         X_u, ql_m )

implicit none

! Input
  
integer, intent(in) :: n, d

double precision, intent(in) :: a, C1, C2
double precision, intent(in) :: ql(1:n) 
!double precision, intent(in) :: w(1:n)
!double precision, intent(in) :: N_pts(1:n), rr(1:n)
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
if (a > 1.0d0 .or. a < 0.0d0) then
   print*, 'Error in autoconv_driver: ', & 
             'mixture fraction, a, does not lie in [0,1].'
   stop
endif
if (C1 > 1.0d0 .or. C1 < 0.0d0) then 
   print*, 'Error in autoconv_driver: ', & 
           'cloud fraction 1, C1, does not lie in [0,1].'
   stop
endif
if (C2 > 1.0d0 .or. C2 < 0.0d0) then 
   print*, 'Error in autoconv_driver: ', & 
           'cloud fraction 2, C2, does not lie in [0,1].'
   stop
endif

! Make sure there is some cloud.
if (a*C1 < 0.001d0 .and. (1-a)*C2 < 0.001d0) then 
   print*, 'Error in autoconv_driver: ', & 
              'there is none or almost no cloud!'
endif

! To compute liquid water, need to set coeff=expn=1.
coeff = 1.d0
expn  = 1.d0

! Initialize liquid in each mixture component
ql_m1 = 0.d0
ql_m2 = 0.d0

! Initialize numbers of sample points corresponding 
!    to each mixture component
n1    = 0
n2    = 0

do sample = 1, n

! Choose which mixture fraction we are in.  
! Account for cloud fraction.
! Follow M. E. Johnson (1987), p. 56.
   fraction_1 = a*C1/(a*C1+(1-a)*C2)
   if ( X_u(sample,d+1) < fraction_1 ) then
! Use an idealized formula to compute liquid 
!      in mixture comp. 1
       ql_m1 = ql_m1 + coeff*(ql(sample))**expn
       n1    = n1 + 1    
   else
! Use an idealized formula to compute liquid 
!      in mixture comp. 2
       ql_m2 = ql_m2 + coeff*(ql(sample))**expn
       n2    = n2 + 1    
   endif

! Loop to get new sample
enddo

!! Convert sums to averages.
!! Old code that underestimates if n1 or n2 = 0.
!       if (n1 .eq. 0) then
!          ql_m1 = 0.d0
!       else
!          ql_m1 = ql_m1/n1
!       endif

!       if (n2 .eq. 0) then
!         ql_m2 = 0.d0
!       else
!          ql_m2 = ql_m2/n2
!       endif


! Convert sums to averages.
! If we have no sample points for a certain plume,
!    then we estimate the plume liquid water by the
!    other plume's value.
if (n1 == 0 .and. n2 == 0) then
   print*, 'Error: no sample points in ql_estimate'
   stop
endif

if ( .not. (n1 == 0) ) then
   ql_m1 = ql_m1/n1
endif

if ( .not. (n2 == 0) ) then
   ql_m2 = ql_m2/n2
endif

if (n1 == 0) then
   ql_m1 = ql_m2
endif

if (n2 == 0) then
   ql_m2 = ql_m1
endif

! Grid box average.
ql_m = a*C1*ql_m1 + (1-a)*C2*ql_m2

return
end subroutine ql_estimate
!---------------------------------------------------------------

end module micro_calcs_mod
