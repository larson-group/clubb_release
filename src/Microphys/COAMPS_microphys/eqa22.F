      subroutine eqa22(psacw,qc,qs,tt,termb,esc,pcut,delt,tzero,len)
c
#include "eqa22.prol"
c SCCS IDENTIFICATION:  %W% %G%

      implicit none
c
c***********************************************************************
c           parameters:
c***********************************************************************
c
      integer len
      integer iflag
c
      real esc
      real delt
      real pcut
      real psacw  (len)
      real qc     (len)
      real qs     (len)
      real termb  (len)
      real tt     (len)
      real tzero
c
c***********************************************************************
c          local variables and dynamic storage:
c***********************************************************************
c
      integer i
c
      real thomo
      real tp1
c
      tp1=1./delt

      thomo=233.15
      do i=1,len
        if ((qs(i).gt.pcut.and.qc(i).gt.pcut).and.tt(i).gt.thomo)then 
          psacw(i)=qc(i)*esc*termb(i)
c
c***********************************************************************
c limit rates
c***********************************************************************
c
          psacw(i)=min(psacw(i),qc(i)*tp1)
          if(delt*psacw(i).lt.pcut)psacw(i)=0.
        endif
      enddo
c
c***********************************************************************
c
      return
      end
