!$Id: lh_sampler.F90,v 1.2 2008-07-24 14:11:21 faschinj Exp $
      module lh_sampler_mod

      implicit none

      public :: lh_sampler

      private :: sample_points, latin_hyper_sample, gaus_mixt_points, & 
                 truncate_gaus_mixt, ltqnorm, gaus_condt, & 
                 st_2_rtthl


      private ! Default scope

      contains
!------------------------------------------------------------------------
! subroutine lh_sampler( )

! This subroutine generates a Latin Hypercube sample.  

! Input:  n = number of calls to microphysics (normally=2)
!         nt = num. random samples before sequence repeats (normally=10)
!         d = number of variates (normally=5) 
!         p_matrix = nxd random matrix of integers that's fed to closure
!         cf = cloud fraction, 0 <= cf <= 1
!         pdf_parms = pdf parameters output by closure_new
!         crt1, crt2, cthl1, cthl2 = coefficients for s
!         rrainm  = mean of rr; must have rrainm>0. 

! Output: X_u = nxd Latin hypercube sample from uniform dist 
!         X_nl = Sample from normal-lognormal distribution 
!         sample_flag = logical that tells whether PDF has non-zero micro
!------------------------------------------------------------------------
      subroutine lh_sampler( n, nt, d, p_matrix, & 
                             cf, & 
                             pdf_parms, & 
                             crt1, crt2, cthl1, cthl2, & 
                             rrainm, & 
                             X_u, X_nl, sample_flag)

!      use constants

      implicit none

! Input

      integer, intent(in) :: n, nt, d
! rrainm  = mean of rr; must have rrainm>0. 
      real, intent(in) :: rrainm

! cloud fraction
      real, intent(in) :: cf

! nxd matrix of random integers
      integer, intent(in) :: p_matrix( 1:n , 1:(d+1) )
      real, intent(in)    :: pdf_parms(26)

!     quantities needed to predict higher order moments
      real, intent(in) :: crt1, crt2, cthl1, cthl2

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

      real :: s1, s2
      real :: R1, R2

! Clip the magnitude of the correlation between rt and thl
      real :: rrtthl_reduced
      double precision :: rrtthl_reduced1, rrtthl_reduced2

! Means of s, t, w, N, rr for plumes 1 and 2
      double precision :: mu1(1:d), mu2(1:d)
! Covariance (not correlation) matrix of rt, thl, w, N, rr
!     for plumes 1 and 2
      double precision :: Sigma_rtthlw_1(1:d,1:d),  & 
                          Sigma_rtthlw_2(1:d,1:d)

! N   = droplet number concentration.  [N] = number / mg air
! N1  = PDF parameter for mean of plume 1. [N1] = (#/mg)
! N2  = PDF parameter for mean of plume 2. [N2] = (#/mg)
! sN1,2 = PDF param for width of plume 1,2. [sN1,2] = (#/mg)**2
! Ncm = 65 per cc, Ncp2_Ncm2 = 0.07 is for DYCOMS2 RF02 (in cloud)
      double precision, parameter :: Ncm       = 0.065
      double precision, parameter :: Ncp2_Ncm2 = 0.07 ! someday from constants.F
      double precision :: N1, N2, sN1, sN2

! rr = specific rain content. [rr] = g rain / kg air
! rr1  = PDF parameter for mean of plume 1. [rr1] = (g/kg)
! rr2  = PDF parameter for mean of plume 2. [rr2] = (g/kg)
! srr1,2 = PDF param for width of plume 1,2. [srr1,2] = (g/kg)**2
! rrp2_rrainm2 = rrp2 divided by rrainm^2 []
! rrp2_rrainm2 = 0.4 is for DYCOMS2 RF02 in cloud
      double precision, parameter :: rrp2_rrainm2 = 0.4 
      double precision :: rr1, rr2, srr1, srr2


! Code begins -------------------------------------------


!       Input pdf parameters.

        w1     = pdf_parms(1)
        w2     = pdf_parms(2)
        sw1    = pdf_parms(3)
        sw2    = pdf_parms(4)
        rt1    = pdf_parms(5)
        rt2    = pdf_parms(6)
        srt1   = pdf_parms(7)
        srt2   = pdf_parms(8)
        thl1   = pdf_parms(9)
        thl2   = pdf_parms(10)
        sthl1  = pdf_parms(11)
        sthl2  = pdf_parms(12)
        a      = pdf_parms(13) 
        R1     = pdf_parms(18)
        R2     = pdf_parms(19)
        s1     = pdf_parms(20)
        s2     = pdf_parms(21)
        rrtthl = pdf_parms(24)

!-----------------------------------------------------------------------
!
! Call Latin Hypercube sampler to compute microphysics. 
!    V. Larson Mar 2004.
! This acts as an interface between the boundary layer scheme
!   and the microphysics.  To add a call to a microphysics scheme,
!   alter two lines in autoconversion_driver.f.
!-----------------------------------------------------------------------

! Use units of [g/kg] to ameliorate numerical roundoff errors.
! We prognose rt-thl-w,
!    but we set means, covariance of N, qr to constants.

      sample_flag = .true.
      if ( cf < 0.001 ) then  
! In this case there are essentially no cloudy points to sample;
! Set sample points to zero.

        X_u(:,:)    = 0.0
        X_nl(:,:)   = 0.0
        sample_flag = .false.

      elseif ( srt1  == 0. .or. srt2  == 0. .or. & 
               sthl1 == 0. .or. sthl2 == 0. .or. & 
               sw1   == 0. .or. sw2   == 0. ) then

! In this case, Sigma_rtthlw matrix is ill-conditioned; 
!     then matrix operations will fail.

!         print*,'srt1=', srt1
!         print*,'srt2=', srt2
!         print*,'sthl1=', sthl1
!         print*,'sthl2=', sthl2
!         print*,'sw1=', sw1
!         print*,'sw2=', sw2

!         print*, 'Covariance matrix of r-thl-w is ill-conditioned'

        X_u(:,:)    = 0.0
        X_nl(:,:)   = 0.0
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
! rrainm  = mean of rr; rrp2 = variance of rr, must have rrp2>0.
! rr1  = PDF parameter for mean of plume 1. [rr1] = (g/kg)
! rr2  = PDF parameter for mean of plume 2. [rr2] = (g/kg)
! srr1,2 = PDF param for width of plume 1,2. [srr1,2] = (g/kg)**2
        rr1  = 0.5*log( (dble(rrainm)**2) / (1. + rrp2_rrainm2) )
        rr2  = rr1
        srr1 = log( 1. + rrp2_rrainm2 )
        srr2 = srr1

! Means of s, t, w, N, rr for Gaussians 1 and 2
        mu1 = (/  dble(1.e3*s1), 0.d0, dble(w1), N1, rr1  /)
        mu2 = (/  dble(1.e3*s2), 0.d0, dble(w2), N2, rr2  /)

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
        Sigma_rtthlw_1(1,1:d)   = (/  & 
                 dble(1.e6*srt1), & 
                 rrtthl_reduced1, & 
                 0.d0, & 
                 0.d0, & 
                 0.d0 & 
                                   /)
          Sigma_rtthlw_1(2,1:d) = (/ & 
                   rrtthl_reduced1,  & 
                   dble(sthl1),   & 
                   0.d0, & 
                   0.d0, & 
                   0.d0 & 
                                   /)     
          Sigma_rtthlw_1(3,1:d) = (/ & 
                   0.d0, & 
                   0.d0, & 
                   dble(sw1), & 
                   0.d0, & 
                   0.d0 & 
                                   /) 
          Sigma_rtthlw_1(4,1:d) = (/ & 
                   0.d0, & 
                   0.d0, & 
                   0.d0, & 
                   sN1,  & 
                   0.d0 & 
                                   /) 
          Sigma_rtthlw_1(5,1:d) = (/ & 
                   0.d0, & 
                   0.d0, & 
                   0.d0, & 
                   0.d0, & 
                   srr1 & 
                                   /) 
          Sigma_rtthlw_2(1,1:d) = (/  & 
                   dble(1.e6*srt2), & 
                   rrtthl_reduced2,  & 
                   0.d0, & 
                   0.d0, & 
                   0.d0  & 
                                   /)
          Sigma_rtthlw_2(2,1:d) = (/ & 
                   rrtthl_reduced2,  & 
                   dble(sthl2),   & 
                   0.d0,  & 
                   0.d0,  & 
                   0.d0  & 
                                   /)     
          Sigma_rtthlw_2(3,1:d) = (/ & 
                   0.d0,    & 
                   0.d0,    & 
                   dble(sw2), & 
                   0.d0,  & 
                   0.d0         & 
                                   /) 
          Sigma_rtthlw_2(4,1:d) = (/ & 
                   0.d0,    & 
                   0.d0,    & 
                   0.d0, & 
                   sN2,  & 
                   0.d0         & 
                                   /) 
          Sigma_rtthlw_2(5,1:d) = (/ & 
                   0.d0, & 
                   0.d0, & 
                   0.d0, & 
                   0.d0, & 
                   srr2 & 
                                   /)


!         print*, 'New call to sample_points -----------------'

! Use units of [g/kg] to ameliorate numerical roundoff
        call sample_points( n, nt, d, p_matrix, dble(a), & 
                            dble(1.e3*rt1), dble(thl1),  & 
                            dble(1.e3*rt2), dble(thl2), & 
                            dble(crt1), dble(1.e3*cthl1),  & 
                            dble(crt2), dble(1.e3*cthl2), & 
                            dble(mu1), dble(mu2),  & 
                            Sigma_rtthlw_1, Sigma_rtthlw_2, & 
                            dble(R1), dble(R2), & 
                            X_u, X_nl)

! End of overall if-then statement for Latin hypercube code
       endif

      return
      end subroutine lh_sampler

!----------------------------------------------------------------------
! Generates n random samples from a d-dim Gaussian-mixture PDF.
! Uses Latin hypercube method.
! To be called from closure of hoc.

! Input:  n = number of calls to microphysics (normally=2)
!         d = number of variates (normally=5) 
!         p_matrix = n x d+1 matrix of random integers.
!         a = mixture fraction of Gaussians
!         rt1, thl1 = mean of rt, thl for Gaus comp 1
!         rt2, thl2 = mean of rt, thl for Gaus comp 2
!         crt1 = coefficient relating rt, s and t for Gaus comp 1
!         cthl1 = coeff relating thl, s and t for component 1
!         crt2 = coefficient relating rt, s and t for component 2
!         cthl2 = coefficient relating thl, s and t for comp. 2
!         mu1, mu2 = d-dimensional column vector of means of 
!                           1st, 2nd components
!         Sigma_rtthlw_1, Sigma_rtthlw_2 = 
!          dxd dimensional covariance matrix for components 1 & 2
!         C1, C2 = cloud fraction associated w/ 
!                            1st, 2nd mixture component

! Output: X_u = Sample drawn from uniform distribution
!         X_nl = Sample that is transformed ultimately to normal-lognormal 

! Columns of Sigma_rtthlw:     1   2   3   4   5
!                              rt  thl w   N   rr
! 
! Columns of Sigma_stw, X_nl:  1   2   3   4   5
!                              s   t   w   N   rr  

! We take samples only from the cloudy part 
!       of the grid box.
! We use units of g/kg.
!----------------------------------------------------------------------

      subroutine sample_points( n, nt, d, p_matrix, a, & 
                                rt1, thl1, rt2, thl2, & 
                                crt1, cthl1, crt2, cthl2, & 
                                mu1, mu2,  & 
                                Sigma_rtthlw_1, Sigma_rtthlw_2, & 
                                C1, C2, & 
                                X_u, X_nl )

!      use math_utilities, only:
!     .    corrcoef,     ! Procedures
!     .    mean,
!     .    std

      implicit none

! Input variables  ---------------------------------------

      integer, intent(in) :: n, nt, d, p_matrix(1:n,1:(d+1))

! Weight of 1st Gaussian, 0 <= a <= 1
      double precision, intent(in) :: a

! Thermodynamic constants for plumes 1 and 2, units of g/kg
      double precision, intent(in) :: rt1, thl1, rt2, thl2
      double precision, intent(in) :: crt1, cthl1, crt2, cthl2

! Latin hypercube variables, i.e. s, t, w, etc.
      double precision, intent(in) :: mu1(1:d), mu2(1:d)

! Covariance matrices of rt, thl, w for each Gaussian
!      Ordering of matrix is rt, thl, w
      double precision, intent(in) :: Sigma_rtthlw_1(1:d,1:d)
      double precision, intent(in) :: Sigma_rtthlw_2(1:d,1:d)

! Cloud fractions for components 1 and 2
      double precision, intent(in) :: C1, C2

! Output variables ---------------------------------------

! Sample drawn from uniform distribution
      double precision, intent(out) :: X_u(1:n,1:(d+1)) 

! Sample that is transformed ultimately to normal-lognormal 
      double precision, intent(out) :: X_nl(1:n,1:d)


! Local variables  ---------------------------------------

      integer col

! Covariance matrices for variables s, t, w for comps 1 & 2
      double precision Sigma_stw_1(1:d,1:d), Sigma_stw_2(1:d,1:d)

! Sample of s points that is drawn only from normal distribution
      double precision s_pts(1:n)

! Total water, theta_l: mean plus perturbations
!      double precision rt(1:n), thl(1:n)

! Beginning of code -----------------------------------------

! Convert each Gaussian from rt-thl-w variables to s-t-w vars.
      call rtpthlp_2_sptp( d, Sigma_rtthlw_1, crt1, cthl1, Sigma_stw_1 )
      call rtpthlp_2_sptp( d, Sigma_rtthlw_2, crt2, cthl2, Sigma_stw_2 )

! Generate Latin hypercube sample, with one extra dimension 
!    for mixture component.
      call latin_hyper_sample( n, nt, d+1, p_matrix, X_u )

! Standard sample for testing purposes when n=2
!     X_u(1,1:(d+1)) = ( / 0.0001d0, 0.46711825945881d0, 
!     .       0.58015016959859d0, 0.61894015386778d0, 0.1d0, 0.1d0  / )
!     X_u(2,1:(d+1)) = ( / 0.999d0, 0.63222458307464d0,     
!     .       0.43642762850981d0, 0.32291562498749d0, 0.1d0, 0.1d0  / )     


! Let s PDF (1st column) be a truncated Gaussian.
! Take sample solely from cloud points.
      col = 1
      call truncate_gaus_mixt( n, d, col, a, mu1, mu2, & 
                               Sigma_stw_1, Sigma_stw_2, C1, C2, X_u, & 
                               s_pts )

! Generate n samples of a d-variate Gaussian mixture
! by transforming Latin hypercube points, X_u.
      call gaus_mixt_points( n, d, a, mu1, mu2,  & 
                             Sigma_stw_1, Sigma_stw_2, & 
                             C1, C2, X_u, s_pts, X_nl )

! Transform s (column 1) and t (column 2) back to rt and thl
! This is only needed if you need rt, thl in your microphysics.
!     call sptp_2_rtpthlp(n,d,a,crt1,cthl1,crt2,cthl2,
!     .                C1,C2,X_nl(1:n,1),X_nl(1:n,2),X_u,rtp,thlp)
!     call st_2_rtthl(n,d,a,rt1,thl1,rt2,thl2,
!     .                crt1,cthl1,crt2,cthl2,
!     .                C1,C2,mu1(1),mu2(1),X_nl(1:n,1),X_nl(1:n,2),X_u,
!     .                rt,thl)

! Compute some diagnostics
!       print*, 'C=', a*C1 + (1-a)*C2
!       print*, 'rtm_anl=', a*rt1+(1-a)*rt2, 'rtm_est=', mean(rt(1:n),n)
!       print*, 'thl_anl=',a*thl1+(1-a)*thl2, 'thlm_est=',mean(thl(1:n),n)
!       print*, 'rtpthlp_coef_est=', corrcoef(rt,thl,n)

! Convert last two variates (columns) (usu N and rr) to lognormal
      if (d > 3) then
         X_nl(1:n,(d-1):d) = exp( X_nl(1:n,(d-1):d) )
      else
         print*, 'd<4 in sampling_driver: ', & 
                 'will not convert last two variates to lognormal'
      endif

! Test diagnostics
!       print*, 'mean X_nl(:,2)=', mean(X_nl(1:n,2),n)
!       print*, 'mean X_nl(:,3)=', mean(X_nl(1:n,3),n)
!       print*, 'mean X_nl(:,4)=', mean(X_nl(1:n,4),n)
!       print*, 'mean X_nl(:,5)=', mean(X_nl(1:n,5),n)

!       print*, 'std X_nl(:,2)=', std(X_nl(1:n,2),n)
!       print*, 'std X_nl(:,3)=', std(X_nl(1:n,3),n)
!       print*, 'std X_nl(:,4)=', std(X_nl(1:n,4),n)
!       print*, 'std X_nl(:,5)=', std(X_nl(1:n,5),n)

!       print*, 'corrcoef X_nl(:,2:3)=', 
!     .            corrcoef( X_nl(1:n,2), X_nl(1:n,3), n )


      return
      end subroutine sample_points

!------------------------------------------------------------------------

!-----------------------------------------------------------------------
! Transform covariance matrix from rt', theta_l' coordinates 
!         to s', t' coordinates.  
! Use linear approximation for s', t'.

! Input:  d = number of variates (normally=5) 
!         crt, cthl = coefficients that define s', t' 
!         Sigma_rtthlw = dxd dimensional covariance matrix of 
!                 rt, thl, w, ...
!           ordering of Sigma is rt, thl, w, ...

! Output: Sigma_stw = covariance matrix in terms of s', t' 
!           ordering of Sigma_stw is s, t, w, ...
!-----------------------------------------------------------------------
      subroutine rtpthlp_2_sptp( d, Sigma_rtthlw, crt, cthl, Sigma_stw )

      use matrix_operations, only: matmult ! Procedure(s)      

      implicit none

! Input

      integer, intent(in) :: d

      double precision, intent(in) :: Sigma_rtthlw(1:d,1:d)
      double precision, intent(in) :: crt, cthl

! Output

      double precision, intent(out) :: Sigma_stw(1:d,1:d) 

! Local

      integer j, k
      double precision T_transpose(d,d), T(d,d), M_int(d,d)

! Check that matrix is large enough (at least 2x2)
      if (d < 3) then
         print*, 'Error: Input matrix too small ', & 
                              'in rtpthlp_2_stpthlp.'
         stop
      endif

! Transform rt-thl-w matrix to s-t-w
!    according to Sigma_stw = T*Sigma_rtthlw*T_transpose.

! Set up transformation matrix, T_transpose
      T_transpose(1,1) = crt
      T_transpose(2,1) = -1.d0*cthl
      T_transpose(1,2) = crt
      T_transpose(2,2) = 1.d0*cthl

! Put zeros in the two off-diagonal blocks
      if (d > 2) then
         do j=1,2
            do k=3,d
               T_transpose(j,k) = 0.d0
               T_transpose(k,j) = 0.d0
            enddo
         enddo
      endif

! Put identity matrix in the lower diagonal block.
      if (d > 2) then
        do j = 3, d
          do k = 3, d
            if (j == k) then
              T_transpose(j,k) = 1.d0
            else
              T_transpose(j,k) = 0.d0
            endif
          enddo
        enddo
      endif

! Multiply to obtain M_int = Sigma_rtthlw * T_transpose  
      call matmult( Sigma_rtthlw, d, d, d, d, T_transpose,  & 
                    d, d, d, d, M_int )

! Set up other transformation matrix, T
      T      = T_transpose
      T(1,2) = -1.d0*cthl
      T(2,1) = crt

! Perform final matrix multiplication: Sigma_stw = T * M_int
      call matmult( T, d, d, d, d, M_int, d, d, d, d, Sigma_stw )
      
      return
      end subroutine rtpthlp_2_sptp

!------------------------------------------------------------------------
!------------------------------------------------------------------------
! subroutine latin_hyper_sample( )

! Generates a matrix X that contains a Latin Hypercube sample.
! The sample is uniformly distributed.  
! See Art B. Owen (2003), ``Quasi-Monte Carlo Sampling," 
!    a chapter from SIGGRAPH 2003

! Input: n = number of calls to microphysics (normally=2)
!        dp1 = number of variates plus 1 (normally=6) 
!        p_matrix = n x dp1 array of permuted integers

! Output: n by dp1 matrix, X, 
!    each row of which is a dp1-dimensional sample
!------------------------------------------------------------------------


        subroutine latin_hyper_sample( n, nt, dp1, p_matrix, X)

        use random, only: ran2 ! Procedure(s)

        implicit none
! Input
  
        integer, intent(in) :: n, nt, dp1, p_matrix(1:n,1:dp1)

! Output

        double precision, intent(out) :: X(1:n,1:dp1)

! Local

        integer j, k, seed
!       integer p_matrix(1:n,1:dp1)

! Continue the old string of random numbers by choosing seed>0
        seed = 1
!       seed = 2

!  Compute random permutation row by row
!       do j=1,dp1
!       ! Generate a column vector of integers from 0 to n-1, 
!       !    whose order is random.
!         call rand_permute( n, p_matrix(1:n,j) )
!       enddo

! Choose values of sample 
!     using permuted vector and random number generator.
        do j = 1,n
          do k = 1,dp1
              X(j,k) = (1.0d0/nt)*(p_matrix(j,k) + ran2(seed));
          enddo
        enddo

!        print*, 'p_matrix(:,1)= ', p_matrix(:,1)
!        print*, 'p_matrix(:,dp1)= ', p_matrix(:,dp1)
!        print*, 'X(:,1)= ', X(:,:)

        return
        end subroutine latin_hyper_sample

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
      subroutine gaus_mixt_points( n, d, a, mu1, mu2, Sigma1, Sigma2, & 
                                   C1, C2, X_u, s_pts, X_gm )
      
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

      integer j, sample
      double precision std_normal(1:d)
      double precision fraction_1


! Handle some possible errors re: proper ranges of a, C1, C2.
      if (a > 1.0d0 .or. a < 0.0d0) then
        print *, 'Error in gaus_mixt_points:  ',& 
                   'mixture fraction, a, does not lie in [0,1].'
         stop
      endif
      if (C1 > 1.0d0 .or. C1 < 0.0d0) then 
        print *, 'Error in gaus_mixt_points: ', & 
                 'cloud fraction 1, C1, does not lie in [0,1].'
         stop
      endif
      if (C2 > 1.0d0 .or. C2 < 0.0d0) then 
        print *, 'Error in gaus_mixt_points: ', & 
                 'cloud fraction 2, C2, does not lie in [0,1].'
        stop
      endif

! Make sure there is some cloud.
      if (a*C1 < 0.001d0 .and. (1-a)*C2 < 0.001d0) then 
        print *, 'Error in gaus_mixt_points: ', & 
                    'there is none or almost no cloud!'
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
          call gaus_condt( n, d, std_normal, mu1, Sigma1, s_pts(sample), & 
                           X_gm(sample, 1:d) )    
        else
          call gaus_condt( n, d, std_normal, mu2, Sigma2, s_pts(sample), & 
                           X_gm(sample, 1:d) )   
        endif

! Loop to get new sample
      enddo

      return
      end subroutine gaus_mixt_points

!------------------------------------------------------------------------
!-----------------------------------------------------------------------
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
!-----------------------------------------------------------------------
      subroutine truncate_gaus_mixt( n, d, col, a, mu1, mu2, & 
                      Sigma1, Sigma2, C1, C2, X_u, truncated_column )
      
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
      


! Handle some possible errors re: proper ranges of a, C1, C2.
      if ( (a > 1.0d0) .or. (a < 0.0d0) ) then
        print*, 'Error in truncate_gaus_mixt: ', & 
                   'mixture fraction, a, does not lie in [0,1].'
        stop
      endif
      if ( (C1 > 1.0d0) .or. (C1 < 0.0d0) ) then 
        print*, 'Error in truncate_gaus_mixt: ', & 
                 'cloud fraction 1, C1, does not lie in [0,1].'
        stop
      endif
      if ( (C2 > 1.0d0) .or. (C2 < 0.0d0) ) then 
        print*, 'Error in truncate_gaus_mixt: ', & 
                 'cloud fraction 2, C2, does not lie in [0,1].'
        stop
      endif

! Make sure there is some cloud.
      if (a*C1 < 0.001d0 .and. (1-a) * C2 < 0.001d0) then 
        print*, 'Error in truncate_gaus_mixt: ', & 
                    'there is none or almost no cloud!'
      endif

! Make s PDF (1st column) a truncated Gaussian.
! This allows us to sample solely from the cloud points.
      do sample = 1, n

! Choose which mixture fraction we are in.  
! Account for cloud fraction.
! Follow M. E. Johnson (1987), p. 56.
        fraction_1 = a * C1 / ( a * C1 + (1.d0 - a) *C2 )
        if ( X_u( sample, d+1 ) < fraction_1 ) then
! Replace first dimension (s) with 
!   sample from cloud (i.e. truncated standard Gaussian)
          s_std = ltqnorm( X_u( sample, col ) * C1 + (1.d0 - C1) ) 
! Convert to nonstandard normal with mean mu1 and variance Sigma1
          truncated_column(sample) =  & 
                     s_std * sqrt( Sigma1(col,col) ) + mu1(col)    
        else
! Replace first dimension (s) with 
!   sample from cloud (i.e. truncated Gaussian)
          s_std = ltqnorm( X_u( sample, col ) * C2 + (1.d0 - C2) ) 
! Convert to nonstandard normal with mean mu2 and variance Sigma2
          truncated_column(sample) =  & 
                        s_std * sqrt( Sigma2(col,col) ) + mu2(col) 
        endif

! Loop to get new sample
      enddo
        
      return
      end subroutine truncate_gaus_mixt

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
        double precision function ltqnorm( p )

        use constants, only: Pi_DP ! Variable(s)

        implicit none


        ! Input Variable(s)

        double precision, intent(in) :: p


        ! Local Variable(s)

        double precision a1, a2, a3, a4, a5, a6, b1, b2, b3, b4, b5, & 
                         c1, c2, c3, c4, c5, c6, d1, d2, d3, d4

        double precision q, r, z, z1, plow, phigh 
       
!       double preciseion e, erf_dp, u

!       Occurs in constants.F now.  Isn't actually used currently.
!        double precision, parameter :: pi=3.1415926d0

! Coefficients in rational approximations.
! equivalent: a(1)=a1, a(2)=a2, and etc, when a(1) is in Matlab. 
! Similarly for b, c, and d's
        parameter (a1 = -3.969683028665376d+01,  & 
                   a2 = 2.209460984245205d+02, & 
                   a3 = -2.759285104469687d+02,  & 
                   a4 = 1.383577518672690d+02, & 
                   a5 = -3.066479806614716d+01,  & 
                   a6 = 2.506628277459239d+00)
        parameter (b1 = -5.447609879822406d+01,  & 
                   b2 = 1.615858368580409d+02, & 
                   b3 = -1.556989798598866d+02,  & 
                   b4 = 6.680131188771972d+01, & 
                   b5 = -1.328068155288572d+01)
        parameter (c1 = -7.784894002430293d-03,  & 
                   c2 = -3.223964580411365d-01, & 
                   c3 = -2.400758277161838d+00,  & 
                   c4 = -2.549732539343734d+00, & 
                   c5 =  4.374664141464968d+00,  & 
                   c6 =  2.938163982698783d+00)
        parameter (d1 =  7.784695709041462d-03,  & 
                   d2 =  3.224671290700398d-01, & 
                   d3 =  2.445134137142996d+00,  & 
                   d4 =  3.754408661907416d+00)

        ! Default initialization
        z = 0.0

!  Define break-points.
        plow  = 0.02425d0
        phigh = 1.d0 - plow

!  Initialize output array. Don't need this in Fortran
!   z = zeros(size(p));

!  Rational approximation for lower region:
        if (p > 0.d0 .and. p < plow) then 
           q = sqrt( -2 * log(p) )
           z = (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/ & 
                     ((((d1*q+d2)*q+d3)*q+d4)*q+1.d0)
!  Rational approximation for central region:
      elseif (p >= plow .and. p <= phigh) then 
         q = p - 0.5d0
         r = q * q
         z = (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q & 
                    /(((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1.d0)
! Rational approximation for upper region:
      elseif (p > phigh .and. p < 1.d0) then
         q  = sqrt( -2.d0 * log(1.d0 - p) )
         z  = -(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) & 
                     /((((d1*q+d2)*q+d3)*q+d4)*q+1.d0)
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

!----------------------------------------------------------------------
! Using Gaussian conditional distributions given s, 
!     convert a standard, uncorrelated Gaussian to one with
!     mean mu and covariance structure Sigma.  
! Follow M. E. Johnson, ``Multivariate Statistical Simulation," p50.

! Input: d = number of variates (normally=5) 
!        std_normal = nxd matrix of n independent samples 
!                     from d-variate standard normal distribution
!        mu = d-dimensional column vector of means of Gaussian
!        Sigma = dxd dimensional covariance matrix
!        s_pt = value of Mellor's s

! Output: nonstd_normal = nxd matrix of n samples 
!                         from d-variate normal distribution 
!                         with mean mu and covariance structure Sigma
!----------------------------------------------------------------------
      subroutine gaus_condt( n, d, std_normal, mu, Sigma, s_pt, & 
                             nonstd_normal )

      use matrix_operations, only: gaussj, matmult ! Procedure(s)

      implicit none

! Input
  
        integer, intent(in) :: n, d

        double precision, intent(in) :: std_normal(1:d)
        double precision, intent(in) :: mu(1:d,1)
        double precision, intent(in) :: Sigma(1:d,1:d)
        double precision, intent(in) :: s_pt

! Output

      double precision, intent(out) :: nonstd_normal(1:d) 

! Local

        integer v
!       integer j

        double precision mu_one(1:d,1), mu_two, mu_condt_const
        double precision Sigma_oneone(1:d,1:d),  & 
                         Sigma_oneone_inv(1:d,1:d)
        double precision Sigma_onetwo(1:d,1), Sigma_twoone(1,1:d)
        double precision Sigma_twotwo, Sigma_condt
! Local intermediate quantities
        double precision Sigma_int(1,1:d)
        double precision dum(1,1)

! First store value of s
      nonstd_normal(1) = s_pt

! Loop over variables t, w, . . . 
! [Conditional distribution of 
!      X_two given X_one] = x_one = std_normal(1:n,v)
! is N [ mu_two + Sigma_twoone*Inverse(Sigma_oneone)*(x_one-mu_one), 
!    Sigma_twotwo - Sigma_twoone*Inverse(Sigma_oneone)*Sigma_onetwo]
! Here 'one' refers to given variables, 'two' refers to values to find. 
! Loop over variables t, w, N, rr:
        do v=2, d

         ! Means     
         mu_one(1:(v-1),1) = mu(1:(v-1),1)
         mu_two            = mu(v,1)

         ! Matrices used to compute variances  
         Sigma_twoone(1,1:(v-1)) = Sigma(v,1:(v-1))
         Sigma_onetwo(1:(v-1),1) = Sigma(1:(v-1),v)
         Sigma_oneone( 1:(v-1) , 1:(v-1) )  & 
              = Sigma( 1:(v-1) , 1:(v-1) )
         Sigma_twotwo = Sigma(v,v)

         ! Check whether input matrix is symmetric.
!         do j=1,(v-1)
!            if ( Sigma_twoone(1,j) .ne. Sigma_onetwo(j,1) ) then
!              print*,'Error: matrix Sigma not symmetric in gaus_condt.'
!               print*,'Sigma in gaus_condt=', Sigma
!            endif
!         enddo

! Compute a needed matrix inverse
         call gaussj(Sigma_oneone,v-1,d,Sigma_oneone_inv)

! Compute an intermediate matrix, Sigma_int(1,1:(v-1)), 
!    that is needed several times below.
         call matmult(Sigma_twoone,      1, v-1, 1,  d, & 
                      Sigma_oneone_inv,v-1, v-1, d,  d, Sigma_int)

! Constant scalar needed to compute mean of non-standard normal.
! mu_condt_const = mu_two - Sigma_twoone*Sigma_oneone_inv*mu_one
         call matmult(Sigma_int,         1, v-1, 1,  d, & 
                         mu_one,       v-1,   1, d,  1, dum)
         mu_condt_const = mu_two - dum(1,1)

! Variance of non-standard normal for variable v.
! Sigma_condt is a scalar.
! Sigma_condt = Sigma_twotwo - 
!                    Sigma_twoone*Sigma_oneone_inv*Sigma_onetwo
         call matmult(Sigma_int,         1, v-1, 1,  d, & 
                      Sigma_onetwo,    v-1,   1, d,  1, dum)
         Sigma_condt = Sigma_twotwo - dum(1,1)

! Compute next element of nonstd_normal from prior elements
! nonstd_normal(v) = sqrt(Sigma_condt)*std_normal(v) 
!      + mu_condt_const ...
!      + Sigma_twoone*Sigma_oneone_inv*nonstd_normal(1:(v-1))
         call matmult(Sigma_int,         1, v-1, 1,  d, & 
                      nonstd_normal,   v-1,   1, d,  1, dum)
         nonstd_normal(v) = sqrt(Sigma_condt)*std_normal(v)  & 
            + mu_condt_const  & 
            + dum(1,1)

      ! Loop to obtain new variable
        enddo

        return
        end subroutine gaus_condt

!-----------------------------------------------------------------------
! subroutine st_2_rtthl( )

! Converts from s, t variables to rt, thl 

! Input: n = number of calls to microphysics (normally=2)
!        d = number of variates (normally=5) 
!        a = mixture fraction of Gaussians
!        crt1, cthl1, crt2, cthl2 = constants from plumes 1 & 2
!        C1, C2 = cloud fraction associated w/ 1st, 2nd 
!                                            mixture component
!        s, t = n-dimensional column vector of Mellor's s and t,
!           including mean and perturbation
!        X_u = nxd Latin hypercube sample 
!                                  from uniform distribution 

! Output: rt, thl: n-dimensional column vectors of rt and thl,
!           including mean and perturbation
!-----------------------------------------------------------------------
      subroutine st_2_rtthl( n, d, a, rt1, thl1, rt2, thl2, & 
                             crt1, cthl1, crt2, cthl2, & 
                             C1, C2, s1, s2, s, t, X_u, rt, thl )
      
      implicit none

! Input

      integer, intent(in) :: n, d

      double precision, intent(in) :: a
      double precision, intent(in) :: rt1, thl1, rt2, thl2
      double precision, intent(in) :: crt1, cthl1, crt2, cthl2
      double precision, intent(in) :: C1, C2
      double precision, intent(in) :: s1, s2
      double precision, intent(in) :: s(1:n), t(1:n)
      double precision, intent(in) :: X_u(1:n,1:(d+1))

! Output

      double precision, intent(out) :: rt(1:n), thl(1:n)

! Local

      integer sample
      double precision fraction_1

! Handle some possible errors re: proper ranges of a, C1, C2.
      if (a > 1.0d0 .or. a < 0.0d0) then
         print*, 'Error in sptp_2_rtpthlp: ', & 
                   'mixture fraction, a, does not lie in [0,1].'
         stop
      endif
      if (C1 > 1.0d0 .or. C1 < 0.0d0) then 
         print*, 'Error in sptp_2_rtpthp: ', & 
                 'cloud fraction 1, C1, does not lie in [0,1].'
         stop
      endif
      if (C2 > 1.0d0 .or. C2 < 0.0d0) then 
         print*, 'Error in sptp_2_rtpthp: ', & 
                 'cloud fraction 2, C2, does not lie in [0,1].'
         stop
      endif

! Make sure there is some cloud.
      if (a*C1 < 0.001d0 .and. (1-a)*C2 < 0.001d0) then 
         print*, 'Error in sptp_2_rtpthp:  ', & 
                    'there is none or almost no cloud!'
      endif

      do sample = 1, n

! Choose which mixture fraction we are in.  
! Account for cloud fraction.
! Follow M. E. Johnson (1987), p. 56.
         fraction_1     = a*C1/(a*C1+(1-a)*C2)
         if ( X_u(sample,d+1) < fraction_1 ) then
            rt(sample)  = rt1 + (0.5d0/crt1)*(s(sample)-s1) +  & 
                               (0.5d0/crt1)*t(sample)
            thl(sample) = thl1 + (-0.5d0/cthl1)*(s(sample)-s1) +  & 
                               (0.5d0/cthl1)*t(sample) 
         else
            rt(sample)  = rt2 + (0.5d0/crt2)*(s(sample)-s2) +  & 
                               (0.5d0/crt2)*t(sample)
            thl(sample) = thl2 + (-0.5d0/cthl2)*(s(sample)-s2) +  & 
                               (0.5d0/cthl2)*t(sample) 
         endif

! Loop to get new sample
      enddo

      return
      end subroutine st_2_rtthl
!------------------------------------------------------------------------


      end module lh_sampler_mod

