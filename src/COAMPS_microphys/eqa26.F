      subroutine eqa26(psdep,qv,qs,qsati,qadj,terma,appbpp,tt
     1                ,rho,snzero,tzero,pcut,delt,len)
c
#include "eqa26.prol"
c SCCS IDENTIFICATION:  %W% %G%

      implicit none
c
c***********************************************************************
c           parameters:
c***********************************************************************
c
      integer len
c
      real appbpp (len)
      real delt
      real pcut
      real psdep  (len)
      real qadj   (len)
      real qs     (len)
      real qsati  (len)
      real qv     (len)
      real rho    (len)
      real snzero
      real terma  (len)
      real tt     (len)
      real tzero
c
c***********************************************************************
c          local variables and dynamic storage:
c***********************************************************************
c
      
      integer i
c
      real rate
      real tp1
      real tp2
c
c************************************************************************
c          local constants
c************************************************************************
c
      tp1=4.0*snzero
      tp2=1.0/delt
c
      do i=1,len
        if (qs(i).gt.pcut.and.tt(i).lt.tzero)then
          psdep(i)=(qv(i)/qsati(i)-1.0)
          rate=(qv(i)-qadj(i))*tp2
          psdep(i)=psdep(i)*appbpp(i)*terma(i)*tp1/rho(i)
c
c************************************************************************
c          limit rates
c************************************************************************
c
          if (psdep(i).gt.0.0.and.psdep(i).ge.rate) psdep(i)=rate
          if (psdep(i).lt.0.0)psdep(i)=MAX(psdep(i),-qs(i)*tp2,rate)
          if(abs(delt*psdep(i)).lt.pcut)psdep(i)=0.
        endif
      enddo
c
      return
      end
