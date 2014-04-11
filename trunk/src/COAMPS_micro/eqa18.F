      subroutine eqa18(pdepi,qv,qi,qsati,qadj,tt,rho,conci
     1                ,appbpp,tice,delt,pcut,j,icon,len)
c
#include "eqa18.prol"
c SCCS IDENTIFICATION:  %W% %G%

      implicit none
c
c***********************************************************************
c           parameters:
c***********************************************************************
c
      integer icon
      integer j
      integer len
c
      real appbpp (len)
      real conci  (len)
      real delt
      real pcut
      real pdepi  (len)
      real qadj   (len)
      real qi     (len)
      real qsati  (len)
      real qv     (len)
      real rho    (len)
      real tice
      real tt     (len)
c
c***********************************************************************
c          local variables and dynamic storage:
c***********************************************************************
c
      integer i
c
      real avgmas
      real cni    (len)
      real dbari     
      real qimax
      real rate
      real rhi
      real tcon
      real tcut
      real tdep
      real tp1
      real tp2
c
c************************************************************************
c          local constants
c************************************************************************
c
      qimax=9.4e-10
      tdep=268.16
      tcon=271.16
      tp1=16.3
      tp2=1.0/delt
c
      if(icon.eq.4)then
       tcut=tdep
      elseif(icon.eq.2)then
       tcut=tcon
      else
       tcut=tice
      endif
c
c************************************************************************
c diagnose concentration from qi and maxx ice mass if
c conci(i)=0. when crystals present (may result from
c advective tendencies)
c************************************************************************
c
      do i=1,len
         cni(i)=0.
         if(qi(i).gt.pcut)then
          if(conci(i).le.0.)then
           cni(i)=rho(i)*qi(i)/qimax
          else
           cni(i)=conci(i)
          endif
         endif
      enddo
c
      do i=1,len
        if (qi(i).gt.pcut.and.tt(i).le.tcut.and.cni(i).gt.0.) then
c
         avgmas=qi(i)*rho(i)/cni(i)
         dbari=min(.01,tp1*sqrt(avgmas))
         rhi=qv(i)/qsati(i)-1.0
         pdepi(i)=4.0*dbari*rhi*cni(i)*appbpp(i)/rho(i)
         rate=(qv(i)-qadj(i))*tp2
c
c************************************************************************
c          limit rates  
c************************************************************************
c
         if (pdepi(i).ge.0.0.and.pdepi(i).ge.rate) pdepi(i)=rate
         if (pdepi(i).lt.0.0) pdepi(i)=MAX(pdepi(i),-qi(i)*tp2,rate)
c
         if(abs(delt*pdepi(i)).lt.pcut)pdepi(i)=0.
c
        endif
      enddo
c
c***********************************************************************
c
      return
      end
