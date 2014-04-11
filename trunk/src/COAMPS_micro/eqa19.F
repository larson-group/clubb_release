      subroutine eqa19(pconv,qi,rho,conci,tt,pcut,cimass,delt,len)
c
#include "eqa19.prol"
c SCCS IDENTIFICATION:  %W% %G%

      implicit none
c
c***********************************************************************
c           parameters:
c***********************************************************************
c
      integer len
c
      real cimass
      real conci  (len)
      real delt
      real pconv  (len)
      real pcut
      real qi     (len)
      real rho    (len)
      real tt     (len)
c
c***********************************************************************
c          local variables and dynamic storage:
c***********************************************************************
c
      integer i
c
      real qcimax
      real tp1
c
c***********************************************************************
c          local constants
c***********************************************************************
c
      tp1=1.0/delt
c
      do i=1,len
        if (conci(i).gt.0.0.and.qi(i).gt.pcut) then
          qcimax=cimass*conci(i)/rho(i)
          if(qi(i).gt.qcimax)then
           pconv(i)=(qi(i)-qcimax)*tp1
           if (delt*pconv(i).lt.pcut) pconv(i)=0.0
          endif
        endif
      enddo
c
c***********************************************************************
c
      return
      end
