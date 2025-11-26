! CVS:  $Id: rayle.F,v 1.2 2006-01-05 20:09:59 dschanen Exp $
! CVS:  $Name: not supported by cvs2svn $

!-----------------------------------------------------------------------

      subroutine rayle
     +              (ncol ,  nlm , mbs ,  ib
     +,              amu0 ,   ri ,  pp , ppl
     +,              tray , wray
     +              )

      use kinds
        
      implicit none
             
!-----------------------------------------------------------------------
! REFERENCES:
! Multitasked version of rayle.f from G. Stephens. Rayle computes the
! optical depth and single scattering albedo due to Rayleigh scattering.
! Laura D. Fowler (slikrock. 08-20-97)

! send comments to laura@slikrock.atmos.colostate.edu and
! partain@atmos.colostate.edu

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
!     INPUT ARGUMENTS:
!     ----------------
      integer (kind=int_kind), intent(in)::
     &  ncol !Length of sub-domain.
     &, nlm  !Number of layers.
     &, mbs  !Number of SW spectral intervals.
     &, ib   !Index of spectral interval.

      real (kind=dbl_kind), intent(in), dimension(mbs)::
     &  ri   !Coefs. for computation of Rayleigh scattering.
    
      real (kind=dbl_kind), intent(in), dimension(ncol)::
     &  amu0 !Cosine of zenith angle                               (-).

      real (kind=dbl_kind), intent(in), dimension(ncol,nlm)::
     &  ppl  !Layer pressure                                     (hPa).

      real (kind=dbl_kind), intent(in), dimension(ncol,nlm+1)::
     &  pp   !Level pressure                                     (hPa).
     
!     OUTPUT ARGUMENTS:
!     -----------------
      real (kind=dbl_kind), intent(out), dimension(ncol,nlm)::
     &  wray !Rayleigh single scattering albedo                    (-).
     &, tray !Rayleigh optical depth                               (-).
     
! LOCAL VARIABLES:
      integer (kind=int_kind)::
     &  i    !Horizontal index.
     &, l    !Vertical index.

      real (kind=dbl_kind)::
     &  u0,x,trp

!-----------------------------------------------------------------------

!---- computes rayleigh absorption:

      do i = 1, ncol
         u0=amu0(i)
            if(ib.eq.1) then
               x = -3.902860e-6*u0*u0+6.120070e-6*u0+4.177440e-6
            else
               x = ri(ib)
            endif   
         do l = 1, nlm
               !*** NAG error
               !trp = 29.267*ppl(i,l)*alog(pp(i,l+1)/pp(i,l) ) 
               trp = 29.267*ppl(i,l)*log(pp(i,l+1)/pp(i,l) )
               tray(i,l) = x*trp
               wray(i,l) = 1.0
         enddo
      enddo

      return
      end subroutine rayle

c-----------------------------------------------------------------------
