! CVS:  $Id: two_rt_sw_bs.F,v 1.3 2006-09-05 22:38:51 dschanen Exp $
! CVS:  $Name: not supported by cvs2svn $

!-----------------------------------------------------------------------
 
      subroutine two_rt_sw_bs
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
! DGBSV: LAPACK simple driver routine for solving band diagonal
!        systems w/ LU decomposition. -dschanen 4 Sept. 2006

! FUNCTIONS CALLED:
!     none.

! INCLUDED COMMONS:
!     none.
 
! ARGUMENT LIST VARIABLES:
! All arrays indexed as nlm correspond to variables defined in the
! middle of layers. All arrays indexed as nlm+1 correspond to variables
! defined at levels at the top and bottom of layers.

!     EXTERNAL CALLS:
!     ----------------
      external ::
     &  dgbsv
 
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
! ----------------
      integer (kind=int_kind)
     &  i     ! Horizontal index.
     &, l     ! Vertical index.

      integer (kind=int_kind), dimension(nlm*4+2)::
     &  indx  ! Pivot index used by LU decomposition

      integer (kind=int_kind)
     &  info  ! Diagnostic for LAPACK band solver
 
      real (kind=dbl_kind), dimension(nlm)::
     &  rrcld   !All sky global reflection
     &, rrclr   !Clear sky global reflection
     &, trcld   !All sky global transmission
     &, trclr   !Clear sky global transmission
     &, sigucld !All sky upwelling source
     &, siguclr !Clear sky upwelling source
     &, sigdcld !All sky downwelling source
     &, sigdclr !Clear sky downwelling source

      ! Modified by dschanen for LAPACK
      real (kind=dbl_kind), dimension(5+5+1,nlm*4+2)::
     &  a       !Diagonal matrix for LAPACK LU decomp.

      real (kind=dbl_kind), dimension(2*5+5+1,nlm*4+2)::
     &  al      !LU Decomposition of a

      real (kind=dbl_kind), dimension(nlm*4+2)::
     &  b       !Vector of sources (for A*x = b)

      real (kind=dbl_kind)
     &  exptaucld
     &, exptauclr
     &, directcld(nlm+1)
     &, directclr(nlm+1)

      real (kind=dbl_kind)
     &      aa ,    bb ,   cc , denom
     &, eggtau ,   eps ,   g3 ,    g4
     &,  ggtau , kappa ,    r ,  rinf
     &,      t ,   oms , taus ,  fact
     &,    asy
      data eps /1.e-02/

      real (kind=dbl_kind), dimension(nlm+1)::     
     &  re , vd
 
! SELECTION RULE VARIABLES:
! -------------------------
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
 
      do 1000 i = 1, ncol
 
!         if (sel_rules) then
!            fail = .false.
!            tauscat = 0.0
!            do l = nlm, 1, -1
!               if (wc(i,l).gt.wcthresh) fail = .true.
!               tauscat = tauscat + wc(i,l)*tau(i,l)
!            enddo
!            if (fail.and.tauscat.ge.tausthresh) goto 2000
! 
!!>> BEGIN SELECTION RULES <<
!!            print *,'selection rules'
!            fddir(i,1) = amu0(i)*slr(i)
!            fddif(i,:) = 0.0
!            do l=1,nlm
!               fddir(i,l+1) = fddir(i,l) * exp(-1.*tau(i,l)/amu0(i))
!            enddo
!
!            fudif(i,nlm+1) = fddir(i,nlm+1) * asdir(i,ib)
! 
!            do l=nlm,1,-1
!               fudif(i,l) = fudif(i,l+1) * exp(-2*tau(i,l))
!            enddo
! 
!            cycle
!         endif
!!>> END SELECTION RULES <<

!>> BEGIN FULL CALCULATION <<
2000     directcld(1) = 0.
         directclr(1) = 1.
         re(1) = 0.
         vd(1) = 0.
!         print *,'full up calculation'
 
!---- 1. DO SHORTWAVE:
         !ALL SKY
         do l = 1, nlm
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
            !cc = oms*slr(i)*fact
            cc = oms*fact
            g3 = 0.5-0.75*asy*amu0(i)
            g4 = 1.-g3
            aa = g3*(t-1./amu0(i))+g4*r
            bb = g4*(t+1./amu0(i))+g3*r
 
            sigucld(l) = cc*((aa-rrcld(l)*bb)-aa*trcld(l)*exptaucld) *
     &                 (b3(i,l)*directcld(l)+(1.-b1(i,l))*directclr(l))
            sigdcld(l) = cc*(-bb*trcld(l)+(bb-rrcld(l)*aa)*exptaucld) *
     &                 (b3(i,l)*directcld(l)+(1.-b1(i,l))*directclr(l))

            !CLEAR SKY
            fact = asyclr(i,l)*asyclr(i,l)
            oms  = ((1.-fact)*wcclr(i,l))/(1.-fact*wcclr(i,l))
            taus   = (1.-fact*wcclr(i,l))*tauclr(i,l)
            asy = asyclr(i,l)/(1.+asyclr(i,l))
 
            exptauclr = exp(-taus/amu0(i))

 
!--- local coefficients:  delta-eddington
            t      = 0.25 * (7. - oms*(4.+3.*asy))
            r      = -0.25 * (1. - oms*(4.-3.*asy))
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
     &          ((1.-b3(i,l))*directcld(l) + b1(i,l)*directclr(l))
            directcld(l+1) = exptaucld *
     &          (b3(i,l)*directcld(l) + (1.-b1(i,l))*directclr(l))
         enddo
 
!---- 2. LOAD A MATRIX, B MATRIX
         a(:,:) = 0.0

        ! Modified by dschanen for LAPACK 
         a(6,1)  = -1.0
         a(2,1) = trcld(1) * b4(i,1)
         a(1,1) = trcld(1) * (1.-b2(i,1))
 
         a(6,2)  = -1.0
         a(3,2)  = trclr(1) * (1.-b4(i,1))
         a(2,2) = trclr(1) * b2(i,1)
 
         a(6,3)  = -1.0
         a(4,3)  = rrcld(1) * b4(i,1)
         a(3,3)  = rrcld(1) * (1.-b2(i,1))
 
         a(6,4)  = -1.0
         a(5,4)  = rrclr(1) * (1.-b4(i,1))
         a(4,4)  = rrclr(1) * b2(i,1)
 
         do l = 2,nlm
            a(8,l*4-3)  = rrcld(l) * b3(i,l)
            a(7,l*4-3)  = rrcld(l) * (1.-b1(i,l))
            a(6,l*4-3)  = -1.0
            a(2,l*4-3)  = trcld(l) * b4(i,l)
            a(1,l*4-3)  = trcld(l) * (1.-b2(i,l))
 
            a(9,l*4-2)  = rrclr(l) * (1.-b3(i,l))
            a(8,l*4-2)  = rrclr(l) * b1(i,l)
            a(6,l*4-2)  = -1.0
            a(3,l*4-2)  = trclr(l) * (1.-b4(i,l))
            a(2,l*4-2)  = trclr(l) * b2(i,l)
 
            a(10,l*4-1) = trcld(l) * b3(i,l)
            a(9,l*4-1)  = trcld(l) * (1.-b1(i,l))
            a(6,l*4-1)  = -1.0
            a(4,l*4-1)  = rrcld(l) * b4(i,l)
            a(3,l*4-1)  = rrcld(l) * (1.-b2(i,l))
 
            a(11,l*4)   = trclr(l) * (1.-b3(i,l))
            a(10,l*4)   = trclr(l) * b1(i,l)
            a(6,l*4)    = -1.0
            a(5,l*4)    = rrclr(l) * (1.-b4(i,l))
            a(4,l*4)    = rrclr(l) * b2(i,l)
         enddo
 
         a(8,nlm*4+1) = asdif(i,ib)
         a(6,nlm*4+1) = -1.0
         a(8,nlm*4+2) = asdif(i,ib)
         a(6,nlm*4+2) = -1.0
         ! End LAPACK modification

         b(:) = 0.0
         do l = 1,nlm
!            b(l*4-3) = -1.*(b3(i,l)*sigucld(l) + (1.-b1(i,l)) *siguclr(l))
!            b(l*4-2) = -1.*((1.-b3(i,l))*sigucld(l) + b1(i,l) *siguclr(l))
!            b(l*4-1) = -1.*(b3(i,l)*sigdcld(l) + (1.-b1(i,l)) *sigdclr(l))
!            b(l*4)   = -1.*((1.-b3(i,l))*sigdcld(l) + b1(i,l) *sigdclr(l))
            b(l*4-3) = -sigucld(l)
            b(l*4-2) = -siguclr(l)
            b(l*4-1) = -sigdcld(l)
            b(l*4)   = -sigdclr(l)
         enddo
         b(nlm*4+1) = -asdir(i,ib)*directcld(nlm+1)*amu0(i)
         b(nlm*4+2) = -asdir(i,ib)*directclr(nlm+1)*amu0(i)
 
!         do l = 1,nlm*4+2
!            print '(15E15.7)',(a(l,k),k=1,11)
!         enddo
 
!         do l = 1,nlm*4+2
!            print *,l,b(l)
!         enddo

         ! Modified to use LAPACK -dschanen 4 Sept. 2006
         do l = 1, 5, 1
          a(5+1+l,1:nlm*4+2)
     .    = eoshift( a(5+1+l,1:nlm*4+2), l )
         end do 

         do l = 1, 5, 1
           a(5+1-l,1:nlm*4+2)
     .     = eoshift( a(5+1-l,1:nlm*4+2), -l )
         end do

         ! Copy LHS into decomposition scratch space
         al(5+1:2*5+5+1, 1:nlm*4+2) = a(1:5+5+1, 1:nlm*4+2)

         ! If dbl_kind isn't the same double precision as LAPACK, then
         ! this will fail
         call dgbsv( nlm*4+2, 5, 5, 1, al, 5*2+5+1, 
     &               indx, b, nlm*4+2, info )

         if ( info /= 0 ) then
           ! Matrix is singular, singular to working precision, or had
           ! incorrect arguments.
           stop "BUGSrad: two_rt_sw_bs LU decomp. failed"
         end if
         ! End LAPACK modification
 
!---- 3. SUM CLEAR AND CLOUDY FLUXES
         do l = 1,nlm+1
            fudif(i,l) = slr(i)*(b(l*4-3)+b(l*4-2))
         enddo
         do l = 1,nlm
            fddif(i,l+1) = slr(i)*(b(l*4-1)+b(l*4))
         enddo
 
         fddir(i,1) = amu0(i)*slr(i)
         do l = 1, nlm
            fddir(i,l+1) = amu0(i)*slr(i)*
     &                     (directcld(l+1)+directclr(l+1))
         enddo
!         do l = 1, nlm+1
!            print *,ib,l,fddir(l),fddif(l),fudif(l)
!         enddo
!>> END FULL CALCULATION <<

 1000 continue
 
      return
      end subroutine two_rt_sw_bs
