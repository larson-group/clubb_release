! CVS:  $Id: two_rt_lw_bs.F,v 1.1 2005-10-27 20:06:50 dschanen Exp $
! CVS:  $Name: not supported by cvs2svn $

!-----------------------------------------------------------------------
 
      subroutine two_rt_lw_bs
     +              (
     +                    ncol ,    nlm ,  mbs ,   mbir
     +,                     ib , cldamt ,   wc ,  wcclr
     +,                   asym , asyclr ,  tau , tauclr
     +,                     es ,     bf ,   fu ,     fd
     +,              sel_rules ,     b1 ,   b2 ,     b3
     +,                     b4
     +              )

      use kinds
      use bandsolve
#ifdef usenewexp
      use newexp, only: exp
#endif
      implicit none

!-----------------------------------------------------------------------
! REFERENCES:
! two_rt_lw replaces two_rt and add written by G. Stephens. two_rt_lw
! computes the spectral fluxes using a two-stream approximation method.
! Philip Partain, Philip Gabriel, and Laura D. Fowler/graben (09-08-99).
 
! MODIFICATIONS:
! * changed declarations to adapt the code from BUGS4 to BUGS5.
!   Laura D. Fowler/slikrock (02-01-00).

! SUBROUTINES CALLED:
!     none.

! FUNCTIONS CALLED:
!     none.

! INCLUDED COMMONS:
!     none.
 
! ARGUMENT LIST VARIABLES:
! All arrays indexed as nlm correspond to variables defined in the
! middle of layers. All arrays indexed as nlm+1 correspond to variables
! defined at levels at the top and bottom of layers.
 
!     INPUT ARGUMENTS:
!     ----------------
      logical (kind=log_kind), intent(in)::
     &  sel_rules

      integer (kind=int_kind), intent(in)::
     &  ncol   !Length of sub-domain.
     &, nlm    !Number of layers.
     &, mbs    !Number of SW spectral intervals.
     &, mbir   !Number of IR spectral intervals.
     &, ib     !Index of spectral interval.
 
      real (kind=dbl_kind), intent(in), dimension(ncol,mbir)::
     &  es    !Spectral surface emissivity                         (-).

      real (kind=dbl_kind), intent(in), dimension(ncol,nlm)::
     &  cldamt !Cloud fraction                                     (-).
     &, wc     !All sky single scattering albedo                   (-).
     &, wcclr  !Clear sky single scattering albedo                 (-).
     &, asym   !All sky asymmetry factor                           (-).
     &, asyclr !Clear sky asymmetry factor                         (-).
     &, tau    !All sky optical depth                              (-).
     &, tauclr !Clear sky optical depth                            (-).
     &, b1     !Cloud overlap parameter                            (-).
     &, b2     !Cloud overlap parameter                            (-).
     &, b3     !Cloud overlap parameter                            (-).
     &, b4     !Cloud overlap parameter                            (-).

      real (kind=dbl_kind), intent(in), dimension(ncol,nlm+1)::
     &  bf    !Planck function                                 (W/m^2).
 
!     OUTPUT ARGUMENTS:
!     -----------------
      real (kind=dbl_kind), intent(out), dimension(ncol,nlm+1)::
     &  fd     !Spectral downward flux                         (W/m^2).
     &, fu     !Spectral upward flux                           (W/m^2).
 
! LOCAL VARIABLES:
      integer (kind=int_kind)
     &  i      !Horizontal index.
     &, l      !Vertical index.
     &, ibms   !Index of spectral interval.

      real (kind=dbl_kind), dimension(nlm)::
     &  rrcld    !All sky global reflection                        (-).
     &, rrclr    !Clear sky global reflection                      (-).
     &, trcld    !All sky global transmission                      (-).
     &, trclr    !Clear sky global transmission                    (-).
     &, sigucld  !All sky upwelling source                         (-).
     &, siguclr  !Clear sky upwelling source                       (-).
     &, sigdcld  !All sky downwelling source                       (-).
     &, sigdclr  !Clear sky downwelling source                     (-).
     &, exptau   !

      real (kind=dbl_kind), dimension(nlm*4+2,11)::
     &  a     !
     &, al
      real (kind=dbl_kind), dimension(nlm*4+2)::
     &  b    !

      integer (kind=int_kind), dimension(nlm*4+2)::
     &  indx


      real (kind=dbl_kind)
     &      aa ,    bb , beta0 ,     cc
     &, diffac , denom ,  fact , eggtau
     &,  ggtau
     &,  kappa ,   oms ,  prop ,r, rinf
     &,      t ,  taus ,d
      data diffac /2./
 
! SELECTION RULE VARIABLES

      logical (kind=log_kind)::
     &  fail

      real (kind=dbl_kind)::
     &  tausthresh
     &, wcthresh
     &, tauscat


!#ifdef usenewexp
!      real(kind=dbl_kind), external :: exp
!#endif

      data tausthresh / 0.001 /
      data wcthresh   / 0.975 /
 
!----------------------------------------------------------------------
!        print *,'two'
!        do l = 1,nlm
!          print *,l,cldamt(1,l)
!        enddo

      fd(:,1) = 0.
 
      ibms = ib - mbs
 
      do 1000 i = 1, ncol
 
!        if (sel_rules) then
!          fail = .false.
!          tauscat = 0.0
!          do l = nlm, 1, -1
!            if (wc(i,l).gt.wcthresh) fail = .true.
!            tauscat = tauscat + wc(i,l)*tau(i,l)
!          enddo
!          if (fail.and.tauscat.ge.tausthresh) goto 2000
! 
!!>> BEGIN SELECTION RULES <<
!!          print *,'selection rules'
!          do l=1,nlm
!            exptau(l) = exp(-2*tau(i,l))
!            if(tau(i,l) .lt. .8e-2) then
!              sigu(l) = (bf(i,l)+bf(i,l+1))*tau(i,l)
!              sigd(l) = sigu(l)
!            else
!              prop = (1.-exptau(l))/tau(i,l)
!              aa = 2.-prop
!              bb = -2.*exptau(l)+prop
!              cc = 0.5
!              sigu(l) = (aa*bf(i,l)+bb*bf(i,l+1))*cc
!              sigd(l) = (bb*bf(i,l)+aa*bf(i,l+1))*cc
!            endif
!            fd(i,l+1) = sigd(l) + exptau(l) * fd(i,l)
!          enddo
! 
!          fu(i,nlm+1) = bf(i,nlm+1)*es(i,ibms)
!!     &                  + fd(i,nlm+1)*(1.0-es(i,ibms))
! 
!          do l=nlm,1,-1
!            fu(i,l) = sigu(l) + exptau(l) * fu(i,l+1)
!          enddo
! 
!          cycle
!        endif
!!>> END SELECTION RULES <<
 
!>> BEGIN FULL CALCULATION <<
2000    continue
!        print *,'full up calculation'

!---- 2. DO LONGWAVE:
        do l = 1, nlm
          !ALL SKY
          fact = asym(i,l)*asym(i,l)
          oms  = ((1.-fact)*wc(i,l))/(1.-fact*wc(i,l))
          taus   = (1.-fact*wc(i,l))*tau(i,l)
 
          beta0 = (4.+asym(i,l))/(8.*(1.+asym(i,l)))
          t = diffac*(1.-oms*(1.-beta0))     !-0.25
          r = diffac*oms*beta0               !-0.25
          kappa = sqrt(t**2-r**2)
          rinf   = r/(kappa+t)
          ggtau  = kappa*taus
          eggtau = exp(-ggtau)
          denom  = (1.-rinf**2*eggtau**2)
          trcld(l) = (1.-rinf**2)*eggtau/denom
          rrcld(l) = rinf*(1.-eggtau**2)/denom

          if(taus .lt. .8e-2) then
            sigucld(l) = cldamt(i,l)*0.5*diffac*(bf(i,l)+
     *                    bf(i,l+1))*taus
            sigdcld(l) = cldamt(i,l)*sigucld(l)
          else
            aa =  (t+r)*(1.-rrcld(l))-(1.+rrcld(l)-trcld(l))/taus
            bb = -(t+r)*trcld(l)+(1.+rrcld(l)-trcld(l))/taus
            cc = diffac*(1.-oms)/kappa**2
            sigucld(l) = cldamt(i,l)*cc*(aa*bf(i,l)+bb*bf(i,l+1))
            sigdcld(l) = cldamt(i,l)*cc*(bb*bf(i,l)+aa*bf(i,l+1))
          endif

          !CLEAR SKY
          fact = asyclr(i,l)*asyclr(i,l)
          oms  = ((1.-fact)*wcclr(i,l))/(1.-fact*wcclr(i,l))
          taus   = (1.-fact*wcclr(i,l))*tauclr(i,l)
 
          beta0 = (4.+asyclr(i,l))/(8.*(1.+asyclr(i,l)))
          t = diffac*(1.-oms*(1.-beta0))     !-0.25
          r = diffac*oms*beta0               !-0.25
          kappa = sqrt(t**2-r**2)
          rinf   = r/(kappa+t)
          ggtau  = kappa*taus
          eggtau = exp(-ggtau)
          denom  = (1.-rinf**2*eggtau**2)
          trclr(l) = (1.-rinf**2)*eggtau/denom
          rrclr(l) = rinf*(1.-eggtau**2)/denom

          if(taus .lt. .8e-2) then
            siguclr(l) = (1.0-cldamt(i,l))*0.5*diffac*(bf(i,l)+
     *                    bf(i,l+1))*taus
            sigdclr(l) = (1.0-cldamt(i,l))*siguclr(l)
          else
            aa =  (t+r)*(1.-rrclr(l))-(1.+rrclr(l)-trclr(l))/taus
            bb = -(t+r)*trclr(l)+(1.+rrclr(l)-trclr(l))/taus
            cc = diffac*(1.-oms)/kappa**2
            siguclr(l) = (1.0-cldamt(i,l))*cc*(aa*bf(i,l)+bb*bf(i,l+1))
            sigdclr(l) = (1.0-cldamt(i,l))*cc*(bb*bf(i,l)+aa*bf(i,l+1))
          endif

        enddo
 
!---- 1. LOAD A MATRIX
        a(:,:) = 0.0

        a(1,6)  = -1.0
        a(1,10) = trcld(1) * b4(i,1)
        a(1,11) = trcld(1) * (1.-b2(i,1))

        a(2,6)  = -1.0
        a(2,9)  = trclr(1) * (1.-b4(i,1))
        a(2,10) = trclr(1) * b2(i,1)

        a(3,6)  = -1.0
        a(3,8)  = rrcld(1) * b4(i,1)
        a(3,9)  = rrcld(1) * (1.-b2(i,1))

        a(4,6)  = -1.0
        a(4,7)  = rrclr(1) * (1.-b4(i,1))
        a(4,8)  = rrclr(1) * b2(i,1)

        do l = 2,nlm
          a(l*4-3,4)  = rrcld(l) * b3(i,l)
          a(l*4-3,5)  = rrcld(l) * (1.-b1(i,l))
          a(l*4-3,6)  = -1.0
          a(l*4-3,10) = trcld(l) * b4(i,l)
          a(l*4-3,11) = trcld(l) * (1.-b2(i,l))

          a(l*4-2,3)  = rrclr(l) * (1.-b3(i,l))
          a(l*4-2,4)  = rrclr(l) * b1(i,l)
          a(l*4-2,6)  = -1.0
          a(l*4-2,9)  = trclr(l) * (1.-b4(i,l))
          a(l*4-2,10) = trclr(l) * b2(i,l)

          a(l*4-1,2)  = trcld(l) * b3(i,l)
          a(l*4-1,3)  = trcld(l) * (1.-b1(i,l))
          a(l*4-1,6)  = -1.0
          a(l*4-1,8)  = rrcld(l) * b4(i,l)
          a(l*4-1,9)  = rrcld(l) * (1.-b2(i,l))

          a(l*4,1)    = trclr(l) * (1.-b3(i,l))
          a(l*4,2)    = trclr(l) * b1(i,l)
          a(l*4,6)    = -1.0
          a(l*4,7)    = rrclr(l) * (1.-b4(i,l))
          a(l*4,8)    = rrclr(l) * b2(i,l)
        enddo

        a(nlm*4+1,4) = 0.0 !1.-esfc(i,ibms)
        a(nlm*4+1,6) = -1.0
        a(nlm*4+2,4) = 0.0 !1.-esfc(i,ibms)
        a(nlm*4+2,6) = -1.0

        do l = 1,nlm
!          b(l*4-3) = -1.*(b3(i,l)*sigucld(l)+(1.-b1(i,l))*siguclr(l))
!          b(l*4-2) = -1.*((1.-b3(i,l)) * sigucld(l)+b1(i,l)*siguclr(l))
!          b(l*4-1) = -1.*(b3(i,l)*sigdcld(l)+(1.-b1(i,l))*sigdclr(l))
!          b(l*4)   = -1.*((1.-b3(i,l)) * sigdcld(l)+b1(i,l)*sigdclr(l))
          b(l*4-3) = -sigucld(l)
          b(l*4-2) = -siguclr(l)
          b(l*4-1) = -sigdcld(l)
          b(l*4)   = -sigdclr(l)
        enddo
        b(nlm*4+1) = -1.*cldamt(i,nlm)*es(i,ibms)*bf(i,nlm+1)
        b(nlm*4+2) = -1.*(1.-cldamt(i,nlm))*es(i,ibms)*bf(i,nlm+1)
        
!        do l = 1,nlm*4+2
!          print *,(a(l,k),k=1,11)
!        enddo

!        do l = 1,nlm*4+2
!          print *,l,b(l)
!        enddo
!         stop

        call bandec(a,nlm*4+2,5,5,nlm*4+2,11,al,11,indx,d)
        call banbks(a,nlm*4+2,5,5,nlm*4+2,11,al,11,indx,b)

!---- 3. SUM CLEAR AND CLOUDY FLUXES
        do l = 1,nlm+1
          fu(i,l) = b(l*4-3)+b(l*4-2)
        enddo
        do l = 1,nlm
          fd(i,l+1) = b(l*4-1)+b(l*4)
        enddo
!        do l = 1,nlm*4+2
!          print *,l,b(l)
!        enddo
!        do l = 1,nlm+1
!          print *,fu(i,l),fd(i,l)
!        enddo
!        stop
 
 1000 continue
 
      return
      end subroutine two_rt_lw_bs
