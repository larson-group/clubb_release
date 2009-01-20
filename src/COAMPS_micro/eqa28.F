      subroutine eqa28(psmlti,qi,tt,pcut,delt,tzero,len)
c
#include "eqa28.prol"
c SCCS IDENTIFICATION:  %W% %G%

      implicit none
c
c***********************************************************************
c           parameters:
c***********************************************************************
c
      integer len
c
      real delt
      real pcut
      real psmlti (len)
      real qi     (len)
      real tt     (len)
      real tzero
c
c***********************************************************************
c          local variables and dynamic storage:
c***********************************************************************
c
      integer i
      real tp1
c
c************************************************************************
c          local constants
c************************************************************************
c
      tp1=1.0/delt
c
      do i=1,len
        if (qi(i).gt.pcut.and.tt(i).ge.tzero) psmlti(i)=qi(i)*tp1
        if (delt*psmlti(i).lt.pcut)psmlti(i)=0.
      enddo
c
c***********************************************************************
c
      return
      end
