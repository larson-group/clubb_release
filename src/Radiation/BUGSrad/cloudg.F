! CVS: $Id: cloudg.F,v 1.2 2006-01-05 20:09:59 dschanen Exp $
! CVS: $Name: not supported by cvs2svn $

!-----------------------------------------------------------------------

      subroutine cloudg
     +              (  ncol ,  nlm ,    mb ,   ib
     +,                  pp ,   tt , wcont ,   re
     +,               pdist , cnrw ,  cniw , cnri
     +,                cnii , xlam ,  tcld , wcld
     +,              asycld , flag
     +              )

      use kinds

      use constants_clubb, only:  fstderr  ! Constant(s)

      implicit none

!-----------------------------------------------------------------------
! REFERENCES:
! Cleaned up version of cloud.f from G. Stephens. Computes the cloud
! optical properties.

! This routine has been modified to include the the corrections to
! ADT for spherical particles based upon the work of David Mitchell
! DRI.  All the derivations have been carried out for the modified
! gamma distribution assuming that m=0.5 (em in program), a 
! parameter in eqn (5) of Mitchell (1994).
 
! tcld, wcld, asycld are the optical depth, single scattering albedo,
! and asymmetry parameter of cloud particles based on the use of
! ADT theory as used by Stephens et al (1990). Effective radius re  
! is input (in microns) and the water content is in g/m3.  The logical
! variable flag is .false. for water and .true. for ice.

! send comments to laura@slikrock.atmos.colostate.edu and
! partain@atmos.colostate.edu.

! MODIFICATIONS:
! * changed declarations to adapt the code from BUGS4 to BUGS5.
!   Laura D. Fowler/slikrock (02-01-00).

! SUBROUTINES CALLED:
!     none.

! FUNCTIONS CALLED:
!     RATGAM:

! INCLUDED COMMONS:
!     none.

! ARGUMENT LIST VARIABLES:
!     INPUT ARGUMENTS:
!     ----------------
      logical (kind=log_kind), intent(in)::
     &  flag   !If true, computes optical properties of ice clouds, of
!              of water clouds otherwise.

      integer (kind=int_kind), intent(in)::
     &  ncol   !Length of sub-domain.
     &, nlm    !Number of layers.
     &, mb     !Total number of spectral intervals.
     &, ib     !Index of spectral interval.

      real (kind=dbl_kind), intent(in), dimension(mb)::
     &  cnrw   !Real part of refractive index (Water clouds).
     &, cniw   !Imaginary part of refractive index (Water clouds).
     &, cnri   !Real part of refractive index (Ice clouds).
     &, cnii   !Imaginary part of refractive index (Ice clouds).
     &, xlam   !Center of spectral band.

      real (kind=dbl_kind), intent(in)::
     &  pdist  !

      real (kind=dbl_kind), intent(in), dimension(ncol,nlm)::
     &  tt     !Temperature                                        (K).
     &, wcont  !Cloud water/ice content                       (g/m^-3).
     &, re     !Cloud effective radius                            (mu).      

      real (kind=dbl_kind), intent(in), dimension(ncol,nlm+1)::
     &  pp     !Level pressure                                   (hPa).
     
!     OUTPUT ARGUMENTS:
!     -----------------
      real (kind=dbl_kind), intent(out), dimension(ncol,nlm)::
     &  tcld   !Cloud optical depth                                (-).
     &, wcld   !Cloud single scattering albedo                     (-).
     &, asycld !Cloud asymmetry factor                             (-).

! LOCAL VARIABLES:
      complex (kind=dbl_kind)::
     &  cm,um
          
      integer (kind=int_kind)::
     &  i, l

      real (kind=dbl_kind)::
     &    abs ,  area ,    c0 ,  c1
     &,   cnr ,   cni ,    dz , eps
     &,   ext ,    f2 ,    f3 ,  no
     &,    p0 ,    p1 ,    p2 ,  pi
     &,    rm ,    xm ,    vm
     
!-----------------------------------------------------------------------
!---- initialize local and output arrays:

      tcld(:,:)   = 0.
      wcld(:,:)   = 0.
      asycld(:,:) = 1.

!--   initialize miscellaneous constants and indices of refraction:

      pi  = acos(-1.)
      eps = 1.e-5
      if(flag) then
         cnr = cnri(ib)
         cni = cnii(ib)
      else
         cnr = cnrw(ib)
         cni = cniw(ib)
      endif

!--   constants depending upon the characteristic width of the distribu
!     tion.(these may be made to vary with hydrometeor species and thus 
!     pdist could be made to depend upon level and column numbers).

!     p0    = 0.
      p0    = pdist
      p1    = p0 + 1.
      p2    = p0 + 2.
      f2    = p1 * p0
      f3    = p2 * f2

!---- calculate cloud optical properties:

      do l = 1, nlm
         do i = 1, ncol 
            if(wcont(i,l) .gt. eps) then
               !*** NAG error ***
               !dz=29.286*alog(pp(i,l+1)/pp(i,l)) * tt(i,l) ! in m
               dz=29.286*log(pp(i,l+1)/pp(i,l)) * tt(i,l)
               rm = re(i,l)/p2
               no = wcont(i,l) / ( (4.*pi/3.)*f3*1.e-6*rm**3 )
               area = 1.e-6*pi*f2*no*rm**2
             c0 = 2.*area
             c1 = c0/f2
!               print *,'area',area,ib,re(i,l)
!      For band 1 (0.5 um) only compute the extinction.  Compute both
!      extinction and absorption coefficients for all other bands.
             xm = 2.*pi*rm/xlam(ib)
             cm = cmplx(cnr,-cni)
               if(ib .eq. 1) then
               um   = 2.*xm*(cnr-1.)*cmplx(0.d0,1.d0)
                 ext  = c0 + 2.*c1*real(p0/(um*(um+1.)**p1)
     +                 + 1./(um**2*(um+1.)**p0)-1./um**2)
              tcld(i,l) = ext*dz

! Kurt Kotenberg set wcld equal to 0 in order to test the new anl. solution.
               wcld(i,l) = 0.999999
!               wcld(i,l) = 0.

               asycld(i,l) = 0.85
               !print *,'wcont,ext,abs',wcont(i,l),ext,abs,cnr,cni
             else
               um    = 2.*xm*(cm-1.)*cmplx(0.d0,1.d0)
                 ext   = c0 + 2.*c1*real( p0/(um*(um+1.)**p1)
     +                 + 1./(um**2*(um+1.)**p0)-1./um**2)
                 vm = 4.*xm*cni
!                  abs = area + c1*( p0/(vm*(vm+1.)**p1)
!     +                     + 1./(vm**2*(vm+1.)**p0) - 1./vm**2 )
               abs = area + c1*sngl(
     +            p0/(dble(vm)*(dble(vm)+1.)**dble(p1))
     +            + 1./(dble(vm)**2*(dble(vm)+1.)**dble(p0))
     +            - 1./dble(vm)**2)
                     tcld(i,l) = ext*dz
!PHIL
               if (ext.lt.abs) ext = abs

! Kurt Kotenberg set wcld = 0 again

               wcld(i,l) = (ext-abs)/ext
!                wcld(i,l) = 0.


              ! print *, 'ib=', ib, 'ext', ext
              ! print *,'wcont,ext,abs',wcont(i,l),ext,abs,cnr,cni
              ! print *,'**',ib,wcld(i,l),area,c1,p0,p1,vm
               if(wcld(i,l) .lt. 0.) then
                  write(fstderr,*) wcld(i,l), ext, abs, wcont(i,l)
                  write(fstderr,*) pp(i,l), pp(i,l+1)
                  write(fstderr,*) tt(i,l)
                  write(fstderr,*) re(i,l)
                  stop
               endif

               asycld(i,l)=0.85
             endif
!               print *,ib,wcld(i,l)
           endif
        enddo
      enddo
!               print *, 'ib=', ib, 'wcld=', wcld
      return
      end subroutine cloudg
       
!-----------------------------------------------------------------------
