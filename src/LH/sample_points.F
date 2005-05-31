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

      subroutine sample_points( n, nt, d, p_matrix, a,
     .                          rt1, thl1, rt2, thl2,
     .                          crt1, cthl1, crt2, cthl2,
     .                          mu1, mu2, 
     .                          Sigma_rtthlw_1, Sigma_rtthlw_2,
     .                          C1, C2,
     .                          X_u, X_nl )

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
      double precision rt(1:n), thl(1:n)

      double precision mean, std, corrcoef


! Beginning of code -----------------------------------------

! Convert each Gaussian from rt-thl-w variables to s-t-w vars.
      call rtpthlp_2_sptp( d, Sigma_rtthlw_1, crt1, cthl1, Sigma_stw_1 )
      call rtpthlp_2_sptp( d, Sigma_rtthlw_2, crt2, cthl2, Sigma_stw_2 )

! Generate Latin hypercube sample, with one extra dimension 
!    for mixture component.
      call latin_hyper_sample( n, nt, d+1, p_matrix, X_u )

! Standard sample for testing purposes when n=2
!	X_u(1,1:(d+1)) = ( / 0.0001d0, 0.46711825945881d0, 
!     .       0.58015016959859d0, 0.61894015386778d0, 0.1d0, 0.1d0  / )
!	X_u(2,1:(d+1)) = ( / 0.999d0, 0.63222458307464d0,     
!     .       0.43642762850981d0, 0.32291562498749d0, 0.1d0, 0.1d0  / )     


! Let s PDF (1st column) be a truncated Gaussian.
! Take sample solely from cloud points.
      col = 1
      call truncate_gaus_mixt( n, d, col, a, mu1, mu2,
     .                         Sigma_stw_1, Sigma_stw_2, C1, C2, X_u,
     .                         s_pts )

! Generate n samples of a d-variate Gaussian mixture
! by transforming Latin hypercube points, X_u.
      call gaus_mixt_points( n, d, a, mu1, mu2, 
     .                       Sigma_stw_1, Sigma_stw_2,
     .                       C1, C2, X_u, s_pts, X_nl )

! Transform s (column 1) and t (column 2) back to rt and thl
! This is only needed if you need rt, thl in your microphysics.
!	call sptp_2_rtpthlp(n,d,a,crt1,cthl1,crt2,cthl2,
!     .                C1,C2,X_nl(1:n,1),X_nl(1:n,2),X_u,rtp,thlp)
!	call st_2_rtthl(n,d,a,rt1,thl1,rt2,thl2,
!     .                crt1,cthl1,crt2,cthl2,
!     .                C1,C2,mu1(1),mu2(1),X_nl(1:n,1),X_nl(1:n,2),X_u,
!     .                rt,thl)

! Compute some diagnostics
!      print*, 'C=', a*C1 + (1-a)*C2
!	print*, 'rtm_anl=', a*rt1+(1-a)*rt2, 'rtm_est=', mean(rt(1:n),n)
!	print*, 'thl_anl=',a*thl1+(1-a)*thl2, 'thlm_est=',mean(thl(1:n),n)
!	print*, 'rtpthlp_coef_est=', corrcoef(rt,thl,n)

! Convert last two variates (columns) (usu N and rr) to lognormal
      if (d > 3) then
         X_nl(1:n,(d-1):d) = exp( X_nl(1:n,(d-1):d) )
      else
         print*, 'd<4 in sampling_driver: 
     .           will not convert last two variates to lognormal'
      endif

! Test diagnostics
!	print*, 'mean X_nl(:,2)=', mean(X_nl(1:n,2),n)
!	print*, 'mean X_nl(:,3)=', mean(X_nl(1:n,3),n)
!	print*, 'mean X_nl(:,4)=', mean(X_nl(1:n,4),n)
!	print*, 'mean X_nl(:,5)=', mean(X_nl(1:n,5),n)

!	print*, 'std X_nl(:,2)=', std(X_nl(1:n,2),n)
!	print*, 'std X_nl(:,3)=', std(X_nl(1:n,3),n)
!	print*, 'std X_nl(:,4)=', std(X_nl(1:n,4),n)
!	print*, 'std X_nl(:,5)=', std(X_nl(1:n,5),n)

!	print*, 'corrcoef X_nl(:,2:3)=', 
!     .            corrcoef( X_nl(1:n,2), X_nl(1:n,3), n )


      return
      end subroutine sample_points
!------------------------------------------------------------------------
