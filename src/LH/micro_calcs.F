!------------------------------------------------------------------------
! subroutine micro_calcs( )

! This subroutine computes microphysical grid box averages,
!    given a Latin Hypercube sample.  

! Input: n = number of calls to microphysics (normally=2)
!        d = number of variates (normally=5) 
!        X_u = nxd Latin hypercube sample from uniform dist 
!        X_nl = Sample from normal-lognormal distribution 
!        sample_flag = logical that tells whether PDF has non-zero micro
!        pdf_parms = pdf parameters output by closure_new

! Output: AKm_est = LH estimate of grid box avg Kessler autoconv [kg/kg]
!         AKm = exact grid box avg Kessler autoconversion [kg/kg]
!         AKstd = exact standard deviation of gba Kessler [kg/kg]
!         AKstd_cld = exact w/in cloud std of gba Kessler [kg/kg]
!         AKm_rcm = exact local gba Kessler auto based on rcm [kg/kg]
!         AKm_rcc = exact local gba Kessler based on w/in cloud rc [kg/kg]
!         rcm_est = LH estimate of grid box avg liquid water [kg/kg]
!------------------------------------------------------------------------

      subroutine micro_calcs( n, d, X_u, X_nl, sample_flag,
     .                        pdf_parms,
     .                        AKm_est_k, AKm_k, AKstd_k, AKstd_cld_k,
     .                        AKm_rcm_k, AKm_rcc_k, rcm_est_k )

      use constants
      use anl_erf, only: erf

      implicit none

! Input

      integer, intent(in) :: n, d

! Sample drawn from uniform distribution
      double precision, intent(in) :: X_u(1:n,1:(d+1)) 

! Sample that is transformed ultimately to normal-lognormal 
      double precision, intent(in) :: X_nl(1:n,1:d)

! Flag that determines whether we have a special case (false)
      logical, intent(in) :: sample_flag

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
      real :: w1, w2, sw1, sw2
      real :: thl1, thl2, sthl1, sthl2
      real :: rt1, rt2, srt1, srt2
      real :: crt1, crt2, cthl1, cthl2
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
      real C_w_cld1, C_w_cld2, w_cld_avg

! Code begins -------------------------------------------


!       save pdf parameters

      w1    = pdf_parms(1)
      w2    = pdf_parms(2)
      sw1   = pdf_parms(3)
      sw2   = pdf_parms(4)
      rt1   = pdf_parms(5)
      rt2   = pdf_parms(6)
      srt1  = pdf_parms(7)
      srt2  = pdf_parms(8)
      thl1  = pdf_parms(9)
      thl2  = pdf_parms(10)
      sthl1 = pdf_parms(11)
      sthl2 = pdf_parms(12)
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

      if ( .not. sample_flag ) then  

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
        call autoconv_driver(n,d,dble(a),dble(R1),dble(R2),X_nl(1:n,1),
     .         X_nl(1:n,3),X_nl(1:n,4),X_nl(1:n,5),X_u,AKm_est_dp)
! Convert to real number
        AKm_est_k = real(AKm_est_dp)

! Compute Monte Carlo estimate of liquid for test purposes.
        call ql_estimate(n,d,dble(a),dble(R1),dble(R2),X_nl(1:n,1),
     .         X_nl(1:n,3),X_nl(1:n,4),X_nl(1:n,5),X_u,rcm_est_dp)
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
        AK1      = K_one * ( (s1-q_crit)*R1_crit 
     .            + ss1*exp(-0.5*sn1_crit**2)/(sqrt(2*pi)) )
        sn2_crit = (s2-q_crit)/ss2
        R2_crit  = 0.5*(1+erf(sn2_crit/sqrt(2.0)))
        AK2      = K_one * ( (s2-q_crit)*R2_crit 
     .            + ss2*exp(-0.5*sn2_crit**2)/(sqrt(2*pi)) )
        AKm_k    = a * AK1 + (1-a) * AK2

! Exact Kessler standard deviation in units of (kg/kg)/s
! For some reason, sometimes AK1var, AK2var are negative
        AK1var   = max( 0.0, K_one * (s1-q_crit) * AK1 
     .           + K_one * K_one * (ss1**2) * R1_crit 
     .           - AK1**2  )
        AK2var   = max( 0.0, K_one * (s2-q_crit) * AK2 
     .           + K_one * K_one * (ss2**2) * R2_crit 
     .           - AK2**2  )  
! This formula is for a grid box average:
        AKstd_k  = sqrt(    a  * ( (AK1-AKm_k)**2 + AK1var )
     .              + (1-a) * ( (AK2-AKm_k)**2 + AK2var )
     .              )
! This formula is for a within-cloud average:
        AKstd_cld_k = sqrt( max( 0.0,  
     .            (1./cf) * ( a  * ( AK1**2 + AK1var )
     .                      + (1-a) * ( AK2**2 + AK2var ) 
     .                      )
     .           - (AKm_k/cf)**2  )
     .                  )

! Kessler autoconversion, using grid box avg liquid, rcm, as input
        AKm_rcm_k = K_one * max( 0.0, rcm-q_crit )

! Kessler ac, using within cloud liquid, rcm/cf, as input
        AKm_rcc_k = cf * K_one * max( 0.0, rcm/cf-q_crit ) 

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
