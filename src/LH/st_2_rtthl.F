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
      subroutine st_2_rtthl( n, d, a, rt1, thl1, rt2, thl2,
     .                       crt1, cthl1, crt2, cthl2,
     .                       C1, C2, s1, s2, s, t, X_u, rt, thl )
      
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
      double precision mean

! Handle some possible errors re: proper ranges of a, C1, C2.
      if (a > 1.0d0 .or. a < 0.0d0) then
         print*, 'Error in sptp_2_rtpthlp: 
     .             mixture fraction, a, does not lie in [0,1].'
         stop
      endif
      if (C1 > 1.0d0 .or. C1 < 0.0d0) then 
         print*, 'Error in sptp_2_rtpthp: 
     .           cloud fraction 1, C1, does not lie in [0,1].'
         stop
      endif
      if (C2 > 1.0d0 .or. C2 < 0.0d0) then 
         print*, 'Error in sptp_2_rtpthp: 
     .           cloud fraction 2, C2, does not lie in [0,1].'
         stop
      endif

! Make sure there is some cloud.
      if (a*C1 < 0.001d0 .and. (1-a)*C2 < 0.001d0) then 
         print*, 'Error in sptp_2_rtpthp: 
     .              there is none or almost no cloud!'
      endif

      do sample = 1, n

! Choose which mixture fraction we are in.  
! Account for cloud fraction.
! Follow M. E. Johnson (1987), p. 56.
         fraction_1     = a*C1/(a*C1+(1-a)*C2)
         if ( X_u(sample,d+1) < fraction_1 ) then
            rt(sample)  = rt1 + (0.5d0/crt1)*(s(sample)-s1) + 
     .                         (0.5d0/crt1)*t(sample)
            thl(sample) = thl1 + (-0.5d0/cthl1)*(s(sample)-s1) + 
     .                         (0.5d0/cthl1)*t(sample) 
         else
            rt(sample)  = rt2 + (0.5d0/crt2)*(s(sample)-s2) + 
     .                         (0.5d0/crt2)*t(sample)
            thl(sample) = thl2 + (-0.5d0/cthl2)*(s(sample)-s2) + 
     .                         (0.5d0/cthl2)*t(sample) 
         endif

! Loop to get new sample
      enddo

      return
      end subroutine st_2_rtthl
!------------------------------------------------------------------------
