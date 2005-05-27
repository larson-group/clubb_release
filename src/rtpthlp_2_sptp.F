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
         print*, 'Error: Input matrix too small 
     .                        in rtpthlp_2_stpthlp.'
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
      call matmult( Sigma_rtthlw, d, d, d, d, T_transpose, 
     .              d, d, d, d, M_int )

! Set up other transformation matrix, T
      T      = T_transpose
      T(1,2) = -1.d0*cthl
      T(2,1) = crt

! Perform final matrix multiplication: Sigma_stw = T * M_int
      call matmult( T, d, d, d, d, M_int, d, d, d, d, Sigma_stw )
      
      return
      end subroutine rtpthlp_2_sptp
!------------------------------------------------------------------------

