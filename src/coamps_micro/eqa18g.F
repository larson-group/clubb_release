      subroutine eqa18g(pgmlt,qg,tt,icomp,termg,rho,tzero
     2       ,delt,pcut,hkaolf,len)
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
      real delt
      real hkaolf
      real pcut
      real pgmlt  (len)
      real qg     (len)
      real rho    (len)
      real termg  (len)
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
      real qlim
      real tp1
c
c************************************************************************
c          compress active points
c************************************************************************
c
      tp1=1.0/delt
      
      npts=0
      do i=1,len
        if (qg(i).gt.pcut.and.tt(i).gt.tzero)then
         npts=npts+1
         icomp(npts)=i
        endif
      enddo
c
      if(npts.eq.0)return
c
c************************************************************************
c          compute rate pgmlt - RH 1984
c************************************************************************
c
      do i=1,npts
c
       pgmlt(icomp(i))=-hkaolf*termg(icomp(i))*(tt(icomp(i))-tzero)
     1                /rho(icomp(i))
c
c************************************************************************
c          limit rates
c************************************************************************
c
       qlim=-qg(icomp(i))*tp1
c
       pgmlt(icomp(i))=max(pgmlt(icomp(i)),qlim)
c
       if(-1.0*pgmlt(icomp(i))*delt.lt.pcut)pgmlt(icomp(i))=0.
c
      enddo
c
c***********************************************************************
c
      return
      end
