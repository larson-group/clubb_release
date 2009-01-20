      subroutine eqa10g(pwacs,qs,qc,rho,tt,ss,ssi,ccn,icomp,pterm4
     1       ,tzero,snzero,rhosno,pi,esc,delt,pcut,aprpr,bsnow,len)
c
c#include <eqa10g.prol>
      implicit none
c
c***********************************************************************
c           parameters:
c***********************************************************************
c
      integer len
      integer icomp(len)
c
      real aprpr
      real bsnow
      real ccn    (len)
      real delt
      real esc
      real garg
      real pcut
      real pi
      real pterm4 (len)
      real pwacs  (len)
      real qc     (len)
      real qs     (len)
      real rho    (len)
      real rhosno
      real snzero
      real ss     (len)
      real ssi    (len)
      real tt     (len)
      real tzero
c
c***********************************************************************
c          local variables and dynamic storage:
c***********************************************************************
c
      integer i
      integer npts
c
      real gm6pb
      real qlim
      real tp1
      real tp2
c
c************************************************************************
c          local constants
c************************************************************************
c
      tp1=esc*pi*pi*rhosno*snzero*aprpr/24.0
      tp2=1.0/delt
      garg=6.0+bsnow
      call gamma(garg,gm6pb)
c
c************************************************************************
c          compress active points
c************************************************************************
c
      npts=0
      do i=1,len
        if (qs(i).gt.pcut.and.qc(i).gt.pcut.and.tt(i).le.tzero
     1       .and.ss(i).gt.0.) then
         npts=npts+1
         icomp(npts)=i
        endif
      enddo
c
      if(npts.eq.0)return
c
      do i=1,npts
       ssi(i)=1.0/(ss(icomp(i)))**(6.+bsnow)
      enddo
c
c************************************************************************
c          compute rate pwacs
c************************************************************************
c
      do i=1,npts
c
       pwacs(icomp(i))=tp1*ccn(icomp(i))*pterm4(icomp(i))*gm6pb*ssi(i)
     1                   /rho(icomp(i))
c
       qlim=qs(icomp(i))*tp2
       pwacs(icomp(i))=min(pwacs(icomp(i)),qlim)

       if(pwacs(icomp(i))*delt.lt.pcut)pwacs(icomp(i))=0.
c
      enddo
c
c***********************************************************************
c
      return
      end
