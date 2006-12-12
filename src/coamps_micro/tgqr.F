      subroutine tgqr(fallv,falln,qr,sr,nr,rbm,pterm4,pcut,
     1                rvr,tvr1,tvr2,tvr3,tvr4,ldrizzle,len)

#include "tvqr.prol"
c
      implicit none
c
c***********************************************************************
c           parameters:
c***********************************************************************
c
      integer len
c
      real fallv  (len)
      real falln  (len) ! Added for HOC on 29 Nov 2006
      real rvr    (len)
      real pcut
      real qr     (len)
      real sr     (len)
      real nr     (len)
      real rbm    (len)
      real pterm4 (len)
c
      logical ldrizzle
c
c***********************************************************************
c          local variables and dynamic storage:
c***********************************************************************
c
      integer i
      integer k
c
      real ainvsf
      real falspd  
      real nrplus
      real qrplus
      real rvcof
c     real rvr
      real temp
      real tvr1
      real tvr2
      real tvr3
      real tvr4
      real x5
c
c***********************************************************************
c          local constants
c***********************************************************************
c
      temp=2.387324e-4

c
c***********************************************************************
c          compute terminal velocity (rh 1983 eqn a3)
c          ainvsf holds inverse slope factor (rh 1983 eqn 3a)
c***********************************************************************
c
ccmic$  do all autoscope
ccmic$1 private(i,x5,ainvsf,falspd)
ccmic$1 shared(kk,npts,p00,prbm,exbm,qr,pcut,sloper,rbm)
ccmic$1 private(qrplus,nrplus,rvcof,rvr,temp,fallv)
ccmic$2 shared(nr,rbm,tvr1,tvr2,tvr3,tvr4,ldrizzle,len)
c
c
c***********************************************************************
c          loop over every grid point at each level
c***********************************************************************
c
      if(.not.ldrizzle)then
        do i=1,len
c
c***********************************************************************
c          If the mixing ratio of raindrops does not exeed cutting level
c          then set fall velocity to zero
c          Note: this velocity formula for rain requires sr to 
c                have units 1/cm, sfraincm passed in)
c***********************************************************************
c
          if (qr(i).le.pcut.or.sr(i).lt.0.) then
            fallv(i)=0.0
          else
            ainvsf=1.0/sr(i)
            falspd=tvr1+ainvsf*tvr2+(ainvsf**2)*tvr3+(ainvsf**3)*tvr4
            x5=-falspd*pterm4(i)
            fallv(i)=(x5-abs(x5))*0.5
          endif
        enddo
c
      else
c
c***********************************************************************
c          Khairoutdinov and Kogan (2000) drizzle parameterization (37)
c          note: mean volume radius (cm)
c***********************************************************************
c
         do i = 1,len
          if (qr(i).le.pcut.or.sr(i).lt.0.) then
            fallv(i)=0.0
          else
            qrplus = max (0.0, qr(i))
            nrplus = max (1.E-7, nr(i))
            rvcof=rbm(i)*temp
c           rvr=(rvcof*qrplus/nrplus)**0.333333
            ! Preserve rvr for stats -dschanen 5 Dec 2006
            rvr(i)=(rvcof*qrplus/nrplus)**0.333333 
c---> mjfalk noticed an error in order of operations
c           fallv(i) = min (0.0, -120.*(rvr(i)-0.2))
            fallv(i) = min (0.0, -120.*rvr(i)+0.2)
c<--- End change for order of operations
c----> Khairoutdinov and Kogan says rain droplet number concentration
c      fall should be calculated differently than rain mixing
c      ratio fall speeds, so we add this:
            falln(i) = min (0.0, -70.*rvr(i)+0.1)
C<---- End nr3 fall speed addition
          end if
         enddo

      endif

c***********************************************************************
c
      return
      end
