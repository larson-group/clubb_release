      subroutine eqa19g(pmltge,qg,tt,qv,qsat,sat,icomp,termg
     1       ,apbp,qadj,rho,tzero,delt,pcut,len)
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
      real apbp   (len)
      real delt
      real pcut
      real pmltge (len)
      real qadj   (len)
      real qg     (len)
      real qv     (len)
      real qsat   (len)
      real rho    (len)
      real sat    (len)
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
      real tp2
      real rate
c
      tp2=1.0/delt
c
c***********************************************************************
c          compute supersaturation
c***********************************************************************
c
      do i=1,len
       sat(i)=qv(i)/qsat(i)-1.0
      enddo
c
c************************************************************************
c          compress active points
c************************************************************************
c
      npts=0
      do i=1,len
        if (qg(i).gt.pcut.and.tt(i).gt.tzero.and.sat(i).lt.0.)then
         npts=npts+1
         icomp(npts)=i
        endif
      enddo
c
      if(npts.eq.0)return
c
c************************************************************************
c          compute rate pmltge - eqa 17a RH 1984
c************************************************************************
c

      do i=1,npts
c
       pmltge(icomp(i))=apbp(icomp(i))*sat(icomp(i))*termg(icomp(i))
     1                /rho(icomp(i))
c
       rate=(qv(icomp(i))-qadj(icomp(i)))*tp2
c
c************************************************************************
c          limit rates
c************************************************************************
c
       if (pmltge(icomp(i)).gt.0.0.and.pmltge(icomp(i)).ge.rate) 
     1         pmltge(icomp(i))=rate
       if (pmltge(icomp(i)).lt.0.0)pmltge(icomp(i))
     1         =MAX(pmltge(icomp(i)),-qg(i)*tp2,rate)

       if(abs( pmltge(icomp(i))*delt ).lt.pcut)pmltge(icomp(i))=0.
c
      enddo
c
c***********************************************************************
c
      return
      end
