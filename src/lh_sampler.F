!------------------------------------------------------------------------
!
! This subroutine generates a Latin Hypercube sample.  
!
! Input:  n = number of calls to microphysics (normally=2)
!         nt = num. random samples before sequence repeats (normally=10)
!         d = number of variates (normally=5) 
!         p_matrix = nxd random matrix of integers that's fed to closure
!         cf = cloud fraction, 0 <= cf <= 1
!         pdf_parms = pdf parameters output by closure_new
!         crt1, crt2, cthl1, cthl2 = coefficients for s
!         rrm  = mean of rr; must have rrm>0. 
! Output: X_u = nxd Latin hypercube sample from uniform dist 
!         X_nl = Sample from normal-lognormal distribution 
!         sample_flag = logical that tells whether PDF has non-zero micro

!
      subroutine lh_sampler( n, nt, d, p_matrix,
     .                       cf,
     .                       pdf_parms,
     .                       crt1, crt2, cthl1, cthl2,
     .                       rrm,
     .                       X_u, X_nl, sample_flag)

!      use constants

      implicit none

! Input

	integer, intent(in) :: n, nt, d
! rrm  = mean of rr; must have rrm>0. 
      real, intent(in) :: rrm

! cloud fraction
      real, intent(in) :: cf

! nxd matrix of random integers
	integer, intent(in) :: p_matrix( 1:n , 1:(d+1) )
      real, intent(in)    :: pdf_parms(26)

! Output

! Sample drawn from uniform distribution
      double precision, intent(out) :: X_u(1:n,1:(d+1)) 

! Sample that is transformed ultimately to normal-lognormal 
	double precision, intent(out) :: X_nl(1:n,1:d)

! A true/false flag that determines whether
!     the PDF allows us to construct a sample
	logical, intent(out) :: sample_flag

! Internal 

      real :: a
      real :: w1, w2, sw1, sw2
      real :: thl1, thl2, sthl1, sthl2
      real :: rt1, rt2, srt1, srt2
!        sub-plume correlation coefficient between rt, thl
!        varies between -1 < rrtthl < 1
      real :: rrtthl

!        quantities needed to predict higher order moments

      real :: crt1, crt2, cthl1, cthl2
      real :: s1, s2
      real :: R1, R2

! Clip the magnitude of the correlation between rt and thl
	real :: rrtthl_reduced
	double precision :: rrtthl_reduced1, rrtthl_reduced2

! Means of s, t, w, N, rr for plumes 1 and 2
	double precision :: mu1(1:d), mu2(1:d)
! Covariance (not correlation) matrix of rt, thl, w, N, rr
!     for plumes 1 and 2
	double precision :: Sigma_rtthlw_1(1:d,1:d), 
     .                    Sigma_rtthlw_2(1:d,1:d)

! N   = droplet number concentration.  [N] = number / mg air
! N1  = PDF parameter for mean of plume 1. [N1] = (#/mg)
! N2  = PDF parameter for mean of plume 2. [N2] = (#/mg)
! sN1,2 = PDF param for width of plume 1,2. [sN1,2] = (#/mg)**2
! Ncm = 65 per cc, Ncp2_Ncm2 = 0.07 is for DYCOMS2 RF02 (in cloud)
      double precision, parameter :: Ncm = 0.065
      double precision, parameter :: Ncp2_Ncm2 = 0.07 ! someday from constants.F
      double precision :: N1, N2, sN1, sN2

! rr = specific rain content. [rr] = g rain / kg air
! rr1  = PDF parameter for mean of plume 1. [rr1] = (g/kg)
! rr2  = PDF parameter for mean of plume 2. [rr2] = (g/kg)
! srr1,2 = PDF param for width of plume 1,2. [srr1,2] = (g/kg)**2
! rrp2_rrm2 = rrp2 divided by rrm^2 []
! rrp2_rrm2 = 0.4 is for DYCOMS2 RF02 in cloud
      double precision, parameter :: rrp2_rrm2 = 0.4 
      double precision :: rr1, rr2, srr1, srr2


! Code begins -------------------------------------------


!       Input pdf parameters.

        w1 = pdf_parms(1)
        w2 = pdf_parms(2)
        sw1 = pdf_parms(3)
        sw2 = pdf_parms(4)
        rt1 = pdf_parms(5)
        rt2 = pdf_parms(6)
        srt1 = pdf_parms(7)
        srt2 = pdf_parms(8)
        thl1 = pdf_parms(9)
        thl2 = pdf_parms(10)
        sthl1 = pdf_parms(11)
        sthl2 = pdf_parms(12)
        a = pdf_parms(13) 
        R1 = pdf_parms(18)
        R2 = pdf_parms(19)
        s1 = pdf_parms(20)
        s2 = pdf_parms(21)
        rrtthl = pdf_parms(24)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Call Latin Hypercube sampler to compute microphysics. 
!    V. Larson Mar 2004.
! This acts as an interface between the boundary layer scheme
!   and the microphysics.  To add a call to a microphysics scheme,
!   alter two lines in autoconversion_driver.f.

! Use units of [g/kg] to ameliorate numerical roundoff errors.
! We prognose rt-thl-w,
!    but we set means, covariance of N, qr to constants.

      sample_flag = .true.
	if ( cf .lt. 0.001 ) then  
! In this case there are essentially no cloudy points to sample;
! Set sample points to zero.

	  X_u(:,:) = 0.0
          X_nl(:,:) = 0.0
	  sample_flag = .false.

	elseif ( srt1  .eq. 0. .or. srt2  .eq. 0. .or.
     .         sthl1 .eq. 0. .or. sthl2 .eq. 0. .or.
     .         sw1 .eq. 0. .or. sw2 .eq. 0. ) then

! In this case, Sigma_rtthlw matrix is ill-conditioned; 
!     then matrix operations will fail.

!	  print*,'srt1=', srt1
!	  print*,'srt2=', srt2
!	  print*,'sthl1=', sthl1
!	  print*,'sthl2=', sthl2
!	  print*,'sw1=', sw1
!	  print*,'sw2=', sw2

!	  print*, 'Covariance matrix of r-thl-w is ill-conditioned'

	  X_u(:,:) = 0.0
          X_nl(:,:) = 0.0
	  sample_flag = .false.

	else

! Compute PDF parameters for N, rr.
! Assume that N, rr obey single-lognormal distributions

! N   = droplet number concentration.  [N] = number / mg air
! Ncm  = mean of N; must have Ncm>0
! Ncp2_Ncm2 = variance of N divided by Ncm^2; must have Np2>0.
! N1  = PDF parameter for mean of plume 1. [N1] = (#/mg)
! N2  = PDF parameter for mean of plume 2. [N2] = (#/mg)
! sN1,2 = PDF param for width of plume 1,2. [sN1,2] = (#/mg)**2
          N1  = 0.5*log( (Ncm**2) / (1. + Ncp2_Ncm2) )
	  N2  = N1
	  sN1 = log( 1. + Ncp2_Ncm2 )
          sN2 = sN1

! rr = specific rain content. [rr] = g rain / kg air
! rrm  = mean of rr; rrp2 = variance of rr, must have rrp2>0.
! rr1  = PDF parameter for mean of plume 1. [rr1] = (g/kg)
! rr2  = PDF parameter for mean of plume 2. [rr2] = (g/kg)
! srr1,2 = PDF param for width of plume 1,2. [srr1,2] = (g/kg)**2
          rr1 = 0.5*log( (dble(rrm)**2) / (1. + rrp2_rrm2) )
          rr2 = rr1
 	  srr1 = log( 1. + rrp2_rrm2 )
          srr2 = srr1

! Means of s, t, w, N, rr for Gaussians 1 and 2
          mu1 = ( /  dble(1.e3*s1), 0.d0, dble(w1), N1, rr1  /)
          mu2 = ( /  dble(1.e3*s2), 0.d0, dble(w2), N2, rr2  /)

! An old subroutine, gaus_rotate, couldn't handle large correlations;
!   I assume the replacement, gaus_condt, has equal trouble.
!   Therefore we input smaller correlations
          rrtthl_reduced = min( 0.99, max( rrtthl, -0.99 ) )

! Within-plume rt-thl correlation terms with rt in g/kg
          rrtthl_reduced1 = dble(rrtthl_reduced*1.d3*sqrt(srt1*sthl1))
          rrtthl_reduced2 = dble(rrtthl_reduced*1.d3*sqrt(srt2*sthl2))

! Covariance (not correlation) matrices of rt-thl-w-N-qr
!    for Gaussians 1 and 2
! For now, assume no within-plume correlation of w,N,qr with
!    any other variables.
	  Sigma_rtthlw_1(1,1:d) = ( / 
     .             dble(1.e6*srt1),
     .             rrtthl_reduced1, 
     .             0.d0, 
     .             0.d0, 
     .             0.d0 
     .                        / )
	  Sigma_rtthlw_1(2,1:d) = ( /
     .             rrtthl_reduced1, 
     .             dble(sthl1),  
     .             0.d0, 
     .             0.d0, 
     .             0.d0 
     .                        / )     
	  Sigma_rtthlw_1(3,1:d) = ( /
     .             0.d0,   
     .             0.d0,   
     .             dble(sw1),
     .             0.d0, 
     .             0.d0        
     .                        / ) 
	  Sigma_rtthlw_1(4,1:d) = ( /
     .             0.d0,   
     .             0.d0,   
     .             0.d0,
     .             sN1, 
     .             0.d0        
     .                        / ) 
	  Sigma_rtthlw_1(5,1:d) = ( /
     .             0.d0,   
     .             0.d0,   
     .             0.d0,
     .             0.d0, 
     .             srr1        
     .                        / ) 

	  Sigma_rtthlw_2(1,1:d) = ( / 
     .             dble(1.e6*srt2),
     .             rrtthl_reduced2, 
     .             0.d0, 
     .             0.d0, 
     .             0.d0 
     .                        / )
	  Sigma_rtthlw_2(2,1:d) = ( /
     .             rrtthl_reduced2, 
     .             dble(sthl2),  
     .             0.d0, 
     .             0.d0, 
     .             0.d0 
     .                        / )     
	  Sigma_rtthlw_2(3,1:d) = ( /
     .             0.d0,   
     .             0.d0,   
     .             dble(sw2),
     .             0.d0, 
     .             0.d0        
     .                        / ) 
	  Sigma_rtthlw_2(4,1:d) = ( /
     .             0.d0,   
     .             0.d0,   
     .             0.d0,
     .             sN2, 
     .             0.d0        
     .                        / ) 
	  Sigma_rtthlw_2(5,1:d) = ( /
     .             0.d0,   
     .             0.d0,   
     .             0.d0,
     .             0.d0, 
     .             srr2        
     .                        / ) 

	  print*, 'New call to sample_points -----------------'

! Use units of [g/kg] to ameliorate numerical roundoff
        call sample_points(n, nt, d, p_matrix, dble(a),
     .                 dble(1.e3*rt1), dble(thl1), 
     .                 dble(1.e3*rt2), dble(thl2),
     .                 dble(crt1), dble(1.e3*cthl1), 
     .                 dble(crt2), dble(1.e3*cthl2),
     .                 dble(mu1), dble(mu2), 
     .                 Sigma_rtthlw_1, Sigma_rtthlw_2,
     .                 dble(R1), dble(R2),
     .                 X_u, X_nl)	

! End of overall if-then statement for Latin hypercube code
	endif

      return
      end
