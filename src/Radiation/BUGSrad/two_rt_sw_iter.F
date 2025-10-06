! CVS:  $Id: two_rt_sw_iter.F,v 1.1 2005-10-27 20:06:50 dschanen Exp $
! CVS:  $Name: not supported by cvs2svn $

!-----------------------------------------------------------------------
 
      subroutine two_rt_sw_iter
     +               ( ncol ,      nlm ,   mbs ,     ib
     +,                 slr ,     amu0 ,    wc ,  wcclr
     +,                asym ,   asyclr ,   tau , tauclr
     +,               asdir ,    asdif , fudif ,  fddir
     +,               fddif ,sel_rules ,    b1 ,     b2
     +,                  b3 ,       b4
     +               )
 

      use kinds
#ifdef usenewexp
      use newexp, only: exp
#endif
      implicit none

!-----------------------------------------------------------------------
! REFERENCES:
! two_rt_sw replaces two_rt and add written by G. Stephens. two_rt_sw
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
     &  ncol  !Length of sub-domain.
     &, nlm   !Number of layers.
     &, mbs   !Number of SW spectral intervals.
     &, ib    !Index of spectral interval.
 
      real (kind=dbl_kind), intent(in), dimension(ncol)::
     &  slr   !Fraction of daylight                                (-).
     &, amu0  !Cosine of the solar zenith angle                    (-).

      real (kind=dbl_kind), intent(in), dimension(ncol,mbs)::
     &  asdir !Spectral direct surface albedo                      (-).
     &, asdif !Spectral diffuse surface albedo                     (-).

      real (kind=dbl_kind), intent(in), dimension(ncol,nlm)::
     &  wc    !Single scattering albedo                            (-).
     &, wcclr !Single scattering albedo                            (-).
     &, asym  !Asymmetry factor                                    (-).
     &, asyclr!Asymmetry factor                                    (-).
     &, tau   !Optical depth                                       (-).
     &, tauclr!Optical depth                                       (-).
     &, b1    !Cloud overlap parameter                             (-).
     &, b2    !Cloud overlap parameter                             (-).
     &, b3    !Cloud overlap parameter                             (-).
     &, b4    !Cloud overlap parameter                             (-).
 
!     OUTPUT ARGUMENTS:
!     -----------------
      real (kind=dbl_kind), intent(out), dimension(ncol,nlm+1)::
     &  fddir !Spectral direct downward flux                   (W/m^2).
     &, fddif !Spectral diffuse downward flux                  (W/m^2).
     &, fudif !Spectral diffuse upward flux                    (W/m^2).
 
! LOCAL VARIABLES:

      integer (kind=int_kind)
     &  i     ! Horizontal index.
     &, l     ! Vertical index.
     &, nsmx
     &, nsr
     &, idc(16*nlm-6)
     &, n
     &, nir(4*nlm+2)       !
     &, j,ii,jj,kk,ir,iter

 
      real (kind=dbl_kind), dimension(nlm)::
     &  rrcld(nlm)       !
     &, rrclr(nlm)       !
     &, trcld(nlm)       !
     &, trclr(nlm)       !
     &, sigucld(nlm)     !
     &, siguclr(nlm)     !
     &, sigdcld(nlm)     !
     &, sigdclr(nlm)     !

      real (kind=dbl_kind), dimension(4*nlm+2)::
     &  b                !
     &, fvc              !
     &, error            !
!     &, old_fvc         !

      real (kind=dbl_kind), dimension(16*nlm-6)::
     &  smx


      real (kind=dbl_kind)
     &  exptaucld
     &, exptauclr
     &, directcld(nlm+1)
     &, directclr(nlm+1)
     &, omega

      real (kind=dbl_kind)
     &      aa ,    bb ,   cc , denom
     &, eggtau ,   eps ,   g3 ,    g4
     &,  ggtau , kappa ,    r ,  rinf
     &,      t ,   oms , taus ,  fact
     &,    asy
      data eps /1.e-02/

      real (kind=dbl_kind), dimension(nlm+1)::
     &  re
     &, vd
 
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

!     data tausthresh / 0.001 /
!     data wcthresh   / 0.975 /
      data tausthresh / 0.01 /
      data wcthresh   / 0.98 /
 
!----------------------------------------------------------------------
      nsr = 4*nlm + 2
      nsmx = 16*nlm - 6
!      wc(:,:) = 0.5
!      wcclr(:,:) = 0.5
 
      do 1000 i = 1, ncol

        fddif(i,:) = 0.0
 
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
!          fddir(i,1) = amu0(i)*slr(i)
!          do l=1,nlm
!            fddir(i,l+1) = fddir(i,l) * exp(-1.*tau(i,l)/amu0(i))
!          enddo
! 
!          fudif(i,nlm+1) = fddir(i,nlm+1) * asdir(i,ib)
! 
!          do l=nlm,1,-1
!            fudif(i,l) = fudif(i,l+1) * exp(-2*tau(i,l))
!          enddo
! 
!          cycle
!        endif
!!>> END SELECTION RULES <<

!> BEGIN FULL CALCULATION <<
2000    directcld(1) = 0.
        directclr(1) = 1.
        re(1) = 0.
        vd(1) = 0.
!        print *,'full up calculation'
 
!---- 1. DO SHORTWAVE:
        do l = 1, nlm
          !CLOUDY SKY --------------------------------------------------
          fact = asym(i,l)*asym(i,l)
          oms  = ((1.-fact)*wc(i,l))/(1.-fact*wc(i,l))
          taus   = (1.-fact*wc(i,l))*tau(i,l)
          asy = asym(i,l)/(1.+asym(i,l))
 
          exptaucld = exp(-taus/amu0(i))
 
!--- local coefficients:  delta-eddington
          t     = 0.25 * (7. - oms*(4.+3.*asy))
          r     = -0.25 * (1. - oms*(4.-3.*asy))
          kappa  = sqrt(t**2-r**2)
          rinf   = r/(kappa+t)
          ggtau  = kappa*taus
          eggtau = exp(-ggtau)
          denom  = (1.-rinf**2*eggtau**2)
          trcld(l) = (1.-rinf**2)*eggtau/denom
          rrcld(l) = rinf*(1.-eggtau**2)/denom
 
          if(abs(kappa**2-1./amu0(i)**2) .lt. eps) then
            fact = 1./eps
          else
            fact = 1./(kappa**2-1./amu0(i)**2)
          endif
          g3 = 0.5-0.75*asy*amu0(i)
          g4 = 1.-g3

          aa = g3*(t-1./amu0(i))+g4*r
          bb = g4*(t+1./amu0(i))+g3*r
          !cc = oms*slr(i)*fact
          cc = oms*fact

          sigucld(l) = cc*((aa-rrcld(l)*bb)-aa*trcld(l)*exptaucld) * 
     &                 (b3(i,l)*directcld(l)+(1.-b1(i,l))*directclr(l))
          sigdcld(l) = cc*(-bb*trcld(l)+(bb-rrcld(l)*aa)*exptaucld) * 
     &                 (b3(i,l)*directcld(l)+(1.-b1(i,l))*directclr(l))

          !CLEAR SKY ---------------------------------------------------
          fact = asyclr(i,l)*asyclr(i,l)
          oms  = ((1.-fact)*wcclr(i,l))/(1.-fact*wcclr(i,l))
          taus   = (1.-fact*wcclr(i,l))*tauclr(i,l)
          asy = asyclr(i,l)/(1.+asyclr(i,l))
 
          exptauclr = exp(-taus/amu0(i))
 
!--- local coefficients:  delta-eddington
          t     = 0.25 * (7. - oms*(4.+3.*asy))
          r     = -0.25 * (1. - oms*(4.-3.*asy))
          kappa  = sqrt(t**2-r**2)
          rinf   = r/(kappa+t)
          ggtau  = kappa*taus
          eggtau = exp(-ggtau)
          denom  = (1.-rinf**2*eggtau**2)
          trclr(l) = (1.-rinf**2)*eggtau/denom
          rrclr(l) = rinf*(1.-eggtau**2)/denom
 
          if(abs(kappa**2-1./amu0(i)**2) .lt. eps) then
            fact = 1./eps
          else
            fact = 1./(kappa**2-1./amu0(i)**2)
          endif
          !cc = oms*slr(i)*fact
          cc = oms*fact
          g3 = 0.5-0.75*asy*amu0(i)
          g4 = 1.-g3
          aa = g3*(t-1./amu0(i))+g4*r
          bb = g4*(t+1./amu0(i))+g3*r
          siguclr(l) = cc*((aa-rrclr(l)*bb)-aa*trclr(l)*exptauclr) * 
     &                 (b1(i,l)*directclr(l)+(1.-b3(i,l))*directcld(l))
          sigdclr(l) = cc*(-bb*trclr(l)+(bb-rrclr(l)*aa)*exptauclr) * 
     &                 (b1(i,l)*directclr(l)+(1.-b3(i,l))*directcld(l))

          directclr(l+1) = exptauclr *
     &        ((1.-b3(i,l))*directcld(l) + b1(i,l)*directclr(l))
          directcld(l+1) = exptaucld *
     &        (b3(i,l)*directcld(l) + (1.-b1(i,l))*directclr(l))

        enddo


!---- 1. LOAD SMX VECTOR
        nir(:) = 4
 
        idc(1) = 5
        idc(2) = 6
        smx(1) = -trcld(1) * b4(i,1)
        smx(2) = -trcld(1) * (1.-b2(i,1))
        nir(1) = 2
 
        idc(3) = 5
        idc(4) = 6
        smx(3) = -trclr(1) * (1.-b4(i,1))
        smx(4) = -trclr(1) * b2(i,1)
        nir(2) = 2
 
        idc(5) = 5
        idc(6) = 6
        smx(5) = -rrcld(1) * b4(i,1)
        smx(6) = -rrcld(1) * (1.-b2(i,1))
        nir(3) = 2
 
        idc(7) = 5
        idc(8) = 6
        smx(7) = -rrclr(1) * (1.-b4(i,1))
        smx(8) = -rrclr(1) * b2(i,1)
        nir(4) = 2
 
        do l = 1,nlm-1
          n = (l-1)*16 + 9
          ir = 4*l

          idc(n)   = ir-1
          idc(n+1) = ir
          idc(n+2) = ir+5
          idc(n+3) = ir+6
          smx(n)   = -rrcld(l+1) * b3(i,l+1)
          smx(n+1) = -rrcld(l+1) * (1.-b1(i,l+1))
          smx(n+2) = -trcld(l+1) * b4(i,l+1)
          smx(n+3) = -trcld(l+1) * (1.-b2(i,l+1))
 
          idc(n+4) = ir-1
          idc(n+5) = ir
          idc(n+6) = ir+5
          idc(n+7) = ir+6
          smx(n+4) = -rrclr(l+1) * (1.-b3(i,l+1))
          smx(n+5) = -rrclr(l+1) * b1(i,l+1)
          smx(n+6) = -trclr(l+1) * (1.-b4(i,l+1))
          smx(n+7) = -trclr(l+1) * b2(i,l+1)
 
          idc(n+8) = ir-1
          idc(n+9) = ir
          idc(n+10) = ir+5
          idc(n+11) = ir+6
          smx(n+8) = -trcld(l+1) * b3(i,l+1)
          smx(n+9) = -trcld(l+1) * (1.-b1(i,l+1))
          smx(n+10) = -rrcld(l+1) * b4(i,l+1)
          smx(n+11) = -rrcld(l+1) * (1.-b2(i,l+1))
 
          idc(n+12) = ir-1
          idc(n+13) = ir
          idc(n+14) = ir+5
          idc(n+15) = ir+6
          smx(n+12) = -trclr(l+1) * (1.-b3(i,l+1))
          smx(n+13) = -trclr(l+1) * b1(i,l+1)
          smx(n+14) = -rrclr(l+1) * (1.-b4(i,l+1))
          smx(n+15) = -rrclr(l+1) * b2(i,l+1)
        enddo

        ir = 4*nlm
 
        idc(16*nlm-7) = 4*nlm-1
        idc(16*nlm-6) = 4*nlm
        smx(16*nlm-7) = -asdif(i,ib)
        smx(16*nlm-6) = -asdif(i,ib)
        nir(ir+1) = 1
        nir(ir+2) = 1
 
        b(:) = 0.0
        do l = 1,nlm
          b(l*4-3) = sigucld(l)
          b(l*4-2) = siguclr(l)
          b(l*4-1) = sigdcld(l)
          b(l*4)   = sigdclr(l)
        enddo
        b(nlm*4+1) = asdir(i,ib)*directcld(nlm+1)
        b(nlm*4+2) = asdir(i,ib)*directclr(nlm+1)
 
 
 
 
!-------------- GAUSS SEIDEL ------------------------------------
!        fvc(:) = b(:)
!
!        do iter=1,200
!           old_fvc(:) = fvc(:)
!           kk = 1
!           do ii = 1,nsr
!             fvc(ii) = b(ii)
!             do j=1,nir(ii)
!               jj = idc(kk)
!               fvc(ii) = fvc(ii) - smx(kk)*fvc(jj)
!               kk = kk + 1
!             enddo
!           enddo
!        enddo
 
!-------------- GAUSS SEIDEL W/ OVERRELAXATION ------------------
        omega = 1. !1.08

        fvc(:) = b(:)
 
        do iter=1,200
           kk = 1
           do ii=1,nsr
             t = 0.0
             do j=1,nir(ii)
               jj = idc(kk)
               t = t + smx(kk) * fvc(jj)
               kk = kk + 1
             enddo
             t = t + fvc(ii)
             fvc(ii) = fvc(ii) + omega * (b(ii)-t)
             error(ii) = b(ii) - t
           enddo
 
           if (maxval(abs(error)) .le. 0.0001) then
!           if (maxval(abs(error)) .le. 0.000001) then
!             print *,omega,iter,' iterations'
             exit
           endif
        enddo

!-------------- JACOBI ------------------------------------------
!        fvc(:) = b(:)
!
!        do iter=1,200
!           kk = 1
!           do ii = 1,nsr
!             old_fvc(ii) = fvc(ii)
!             fvc(ii) = b(ii)
!             do j=1,nir(ii)
!               jj = idc(kk)
!              fvc(ii) = fvc(ii) - smx(kk)*old_fvc(jj)
!               kk = kk + 1
!             enddo
!           enddo
!
!           kk = 1
!           do ii=1,nsr
!             t = 0.0
!             do j=1,nir(ii)
!               jj = idc(kk)
!               t = t + smx(kk) * fvc(jj)
!               kk = kk + 1
!             enddo
!             t = t + fvc(ii)
!           enddo
!        enddo
!----------------------------------------------------------------

!---- 3. SUM CLEAR AND CLOUDY FLUXES
        do l = 1,nlm+1
          fudif(i,l) = slr(i)*(fvc(l*4-3)+fvc(l*4-2))
        enddo
        do l = 1,nlm
          fddif(i,l+1) = slr(i)*(fvc(l*4-1)+fvc(l*4))
        enddo

!!---- 3. FINALLY ADD DIRECT BEAM TO THE DOWNWARD FLUX:
        fddir(i,1) = amu0(i)*slr(i)
        do l = 1, nlm
          fddir(i,l+1) = amu0(i)*slr(i)*
     &                   (directcld(l+1)+directclr(l+1))
        enddo
!        print *,ib,wc(i,1),wcclr(i,1)
!        do l = 1, nlm+1
!          !print *,ib,l,directcld(l),directclr(l),fddir(l),fddif(l),fudif(l)
!          print *,ib,l,fddir(l),fddif(l),fudif(l)
!        enddo
!!>> END FULL CALCULATION <<
!
 1000 continue
 
      return
      end
 
!------------------------------------------------------------------------
