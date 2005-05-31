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
      subroutine gaus_condt( n, d, std_normal, mu, Sigma, s_pt,
     .                       nonstd_normal )
      
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

        integer v, j
        double precision mu_one(1:d,1), mu_two, mu_condt_const
        double precision Sigma_oneone(1:d,1:d), 
     .                   Sigma_oneone_inv(1:d,1:d)
        double precision Sigma_onetwo(1:d,1), Sigma_twoone(1,1:d)
        double precision Sigma_twotwo, Sigma_condt
! Local intermediate quantities
        double precision Sigma_int(1,1:d)
        double precision dum

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
         Sigma_oneone( 1:(v-1) , 1:(v-1) ) 
     .        = Sigma( 1:(v-1) , 1:(v-1) )
         Sigma_twotwo = Sigma(v,v)

         ! Check whether input matrix is symmetric.
!         do j=1,(v-1)
!            if ( Sigma_twoone(1,j) .ne. Sigma_onetwo(j,1) ) then
!              print*,'Error: matrix Sigma not symmetric in gaus_condt.'
!	        print*,'Sigma in gaus_condt=', Sigma
!            endif
!         enddo

! Compute a needed matrix inverse
         call gaussj(Sigma_oneone,v-1,d,Sigma_oneone_inv)

! Compute an intermediate matrix, Sigma_int(1,1:(v-1)), 
!    that is needed several times below.
         call matmult(Sigma_twoone,      1, v-1, 1,  d,
     .                Sigma_oneone_inv,v-1, v-1, d,  d, Sigma_int)

! Constant scalar needed to compute mean of non-standard normal.
! mu_condt_const = mu_two - Sigma_twoone*Sigma_oneone_inv*mu_one
         call matmult(Sigma_int,         1, v-1, 1,  d,
     .                   mu_one,       v-1,   1, d,  1, dum)
         mu_condt_const = mu_two - dum

! Variance of non-standard normal for variable v.
! Sigma_condt is a scalar.
! Sigma_condt = Sigma_twotwo - 
!                    Sigma_twoone*Sigma_oneone_inv*Sigma_onetwo
         call matmult(Sigma_int,         1, v-1, 1,  d,
     .                Sigma_onetwo,    v-1,   1, d,  1, dum)
         Sigma_condt = Sigma_twotwo - dum

! Compute next element of nonstd_normal from prior elements
! nonstd_normal(v) = sqrt(Sigma_condt)*std_normal(v) 
!      + mu_condt_const ...
!      + Sigma_twoone*Sigma_oneone_inv*nonstd_normal(1:(v-1))
         call matmult(Sigma_int,         1, v-1, 1,  d,
     .                nonstd_normal,   v-1,   1, d,  1, dum)
         nonstd_normal(v) = sqrt(Sigma_condt)*std_normal(v) 
     .      + mu_condt_const 
     .      + dum

      ! Loop to obtain new variable
        enddo

        return
        end subroutine gaus_condt
