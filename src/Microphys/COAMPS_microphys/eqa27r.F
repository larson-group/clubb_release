      subroutine eqa27r(piacw,qc,qi,tt,vi,conci,rho,pi,eic
     1                 ,pcut,delt,rholiq,tice,len)
c
#include "eqa22.prol"
      implicit none
c
c***********************************************************************
c           parameters:
c***********************************************************************
c
      integer len
c
      real conci  (len)
      real delt
      real eic
      real pcut
      real pi
      real piacw  (len)
      real qc     (len)
      real qi     (len)
      real rho    (len)
      real rholiq
      real tice
      real tt     (len)
      real vi     (len)
c
c***********************************************************************
c          local variables and dynamic storage:
c***********************************************************************
c
      integer i
c
      real di         
      real dc
      real phi
      real rate
      real rhoi
      real thomo
      real tp1
      real tp2
c
c***********************************************************************
c
      tp1=1./delt
      tp2=pi/4
      rhoi=500.
      thomo=233.15
c
c***********************************************************************
c used fixed droplet size
c***********************************************************************
c
      dc = 10.e-6
c
       do i=1,len
        if (qi(i).gt.pcut.and.qc(i).gt.pcut.and.tt(i).gt.thomo.and.
     1         tt(i).lt.tice.and.conci(i).gt.0.)then 
           di=min(.001,(6.0*rho(i)*qi(i)/(rhoi*pi*conci(i)))**.3333)
           phi=dc*rholiq*vi(i)/(di*2.24e-4)
           if(phi.lt..25)then
            eic=0.0
           else
            eic=min(1.0,.572*log(phi-.25)+.967)
           endif
           rate = qc(i)*tp1
           piacw(i)=eic*conci(i)*tp2*di*di*abs(vi(i))*qc(i)
c
c***********************************************************************
c limite rate to max availiable cloud water
c***********************************************************************
c
           if(piacw(i).gt.rate)piacw(i)=rate
           if(delt*piacw(i).lt.pcut)piacw(i)=0.
c
         endif
       enddo
c
c***********************************************************************
c
      return
      end
