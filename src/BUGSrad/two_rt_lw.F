! CVS:  $Id: two_rt_lw.F,v 1.1 2005-10-27 20:06:50 dschanen Exp $
! CVS:  $Name: not supported by cvs2svn $

!-----------------------------------------------------------------------
 
      subroutine two_rt_lw
     +              (     ncol , nlm,  mbs , mbir
     +,                     ib ,  wc, asym ,  tau
     +,                     es ,  bf,   fu ,   fd
     +,              sel_rules
     +              )
 
      use kinds
#ifdef usenewexp
      use newexp, only: exp
#endif

      use constants_clubb, only:  fstdout  ! Constant(s)

      use error_code, only: clubb_at_least_debug_level_api  ! Procedure

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
     &  ncol  !Length of sub-domain.
     &, nlm   !Number of layers.
     &, mbs   !Number of SW spectral intervals.
     &, mbir  !Number of IR spectral intervals.
     &, ib    !Index of spectral interval.
 
      real (kind=dbl_kind), intent(in), dimension(ncol,mbir)::
     &  es    !Spectral surface emissivity                         (-).

      real (kind=dbl_kind), intent(in), dimension(ncol,nlm)::
     &  wc    !Single scattering albedo                            (-).
     &, asym  !Asymmetry factor                                    (-).
     &, tau   !Optical depth                                       (-).

      real (kind=dbl_kind), intent(in), dimension(ncol,nlm+1)::
     &  bf    !Planck function                                 (W/m^2).
 
!     OUTPUT ARGUMENTS:
!     -----------------
      real (kind=dbl_kind), intent(out), dimension(ncol,nlm+1)::
     &  fd    !Spectral downward flux                          (W/m^2).
     &, fu    !Spectral upward flux                            (W/m^2).
 
! LOCAL VARIABLES:

      integer (kind=int_kind)
     &  i     !Horizontal index.
     &, l     !Vertical index.
     &, ibms  !Index of spectral interval.
 
      real (kind=dbl_kind), dimension(nlm)::
     &  rr    !
     &, tr    !
     &, sigu  !
     &, sigd  !

      real (kind=dbl_kind)
     &      aa ,    bb , beta0 ,     cc
     &, diffac , denom ,  fact , eggtau
     &,  ggtau
     &,  kappa ,   oms ,  prop ,r, rinf
     &,      t ,  taus
      data diffac /2./

      real (kind=dbl_kind), dimension(nlm)::
     &  td , vu , exptau

      real (kind=dbl_kind), dimension(nlm+1)::
     &  re, vd
 
! SELECTION RULE VARIABLES

      logical (kind=log_kind)::
     &  fail

      real (kind=dbl_kind)::
     &  tausthresh
     &, wcthresh
     &, tauscat

      data tausthresh / 0.001 /
      data wcthresh   / 0.975 /
 
!----------------------------------------------------------------------
      fd(:,1) = 0.
 
      ibms = ib - mbs
       
      do 1000 i = 1, ncol
      
         if(sel_rules) then

            fail = .false.
            tauscat = 0.0
            do l = nlm, 1, -1
               if (wc(i,l).gt.wcthresh) fail = .true.
               tauscat = tauscat + wc(i,l)*tau(i,l)
            enddo
            if (fail.and.tauscat.ge.tausthresh) goto 2000
 
!>> BEGIN SELECTION RULES <<
           if ( clubb_at_least_debug_level_api( 1 ) ) then
              write(fstdout,*) 'selection rules'
           endif
            do l = 1, nlm
               exptau(l) = exp(-2*tau(i,l))
               if(tau(i,l) .lt. .8e-2) then
                  sigu(l) = (bf(i,l)+bf(i,l+1))*tau(i,l)
                  sigd(l) = sigu(l)
               else
                  prop = (1.-exptau(l))/tau(i,l)
                  aa = 2.-prop
                  bb = -2.*exptau(l)+prop
                  cc = 0.5
                  sigu(l) = (aa*bf(i,l)+bb*bf(i,l+1))*cc
                  sigd(l) = (bb*bf(i,l)+aa*bf(i,l+1))*cc
               endif
               fd(i,l+1) = sigd(l) + exptau(l) * fd(i,l)
            enddo
 
            fu(i,nlm+1) = bf(i,nlm+1)*es(i,ibms)
!    &                  + fd(i,nlm+1)*(1.0-es(i,ibms))
 
            do l = nlm , 1, -1
               fu(i,l) = sigu(l) + exptau(l) * fu(i,l+1)
            enddo
 
            cycle

         endif

!>> END SELECTION RULES <<
 
!>> BEGIN FULL CALCULATION <<
2000  re(1) = 0.
      vd(1) = 0.
!      print *,'full up calculation'
 
      do l = 1, nlm
         fact = asym(i,l)*asym(i,l)

! Kurt Kotenberg commented out the above statement and set 
! fact =0, 11 Mar 05 conversation with N. Wood
! (purpose of setting fact =0 is to turn off delta-eddington)
!         fact = 0.0
!         print *, 'fact', fact

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
         tr(l) = (1.-rinf**2)*eggtau/denom
         rr(l) = rinf*(1.-eggtau**2)/denom


! Kurt Kotenberg changed the value for the below threshold to include ALL levels for
! the thin cloud 11 NOV 99 case (profile.dat_coamps_lwc_0p1).  See notes for 
! further detail.  16 JUN 05
         
!         if(taus .lt. .8e-2) then 	!Original value
         if(taus .lt.  .1481e-2) then


            sigu(l) = 0.5*diffac*(bf(i,l)+bf(i,l+1))*taus
            sigd(l) = sigu(l)

!           print *, 'test', l
         else
            aa =  (t+r)*(1.-rr(l))-(1.+rr(l)-tr(l))/taus
            bb = -(t+r)*tr(l)+(1.+rr(l)-tr(l))/taus
            cc = diffac*(1.-oms)/kappa**2
            sigu(l) = cc*(aa*bf(i,l)+bb*bf(i,l+1))
            sigd(l) = cc*(bb*bf(i,l)+aa*bf(i,l+1))

!            print *, 'test1', l
         endif
      enddo
!---- 1. do adding, going from top down:

        do l = 1, nlm
           prop = 1. / (1. - re(l)*rr(l))
           re(l+1) = rr(l) + tr(l)**2*re(l)*prop
           vd(l+1) = sigd(l) + (tr(l)*vd(l)
     +             + tr(l)*re(l)*sigu(l))*prop
           vu(l)   = (rr(l)*vd(l) + sigu(l))*prop
           td(l)   = prop
        enddo
 
!---- 2. calculate fluxes going up through the layers:

        fu(i,nlm+1) = es(i,ibms)*bf(i,nlm+1)
 
        do l = nlm+1, 2, -1
           fd(i,l)   = re(l)*fu(i,l) + vd(l)
           fu(i,l-1) = tr(l-1)*fu(i,l)*td(l-1) + vu(l-1)
        enddo

!>> END FULL CALCULATION <<
 
 1000 continue

!       print*, 'ncol=', ncol, ' ib=', ib, 'fu=', fu
!       print *,'wc_rt_lw=', wc
!       print *,'asym_rt', asym
      return
      end subroutine two_rt_lw
 
!------------------------------------------------------------------------
