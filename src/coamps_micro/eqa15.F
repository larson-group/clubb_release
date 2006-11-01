      subroutine eqa15(pint,tt,qv,qsati,qadj,rho,conci,tice
     1                ,pcut,delt,icon,len)
c
#include "eqa15.prol"
c SCCS IDENTIFICATION:  %W% %G%

      implicit none
c
c***********************************************************************
c           parameters:
c***********************************************************************
c
      integer icon
      integer len
c
      real conci  (len)
      real delt
      real pcut
      real pint   (len)
      real qadj   (len)
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
      real check
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
      tp1=1.0/delt
      tp2=1.0e-12/delt
c
      tdep=268.16
      tcon=271.16
c
      if(icon.eq.4)then
       tcut=tdep
      elseif(icon.eq.2)then
       tcut=tcon
      else
       tcut=tice
      endif

      do i=1,len
        if (qv(i).gt.qsati(i).and.tt(i).le.tcut.and.conci(i).gt.0.) then
          check=(qv(i)-qadj(i))*tp1
          pint(i)=tp2*conci(i)/rho(i)
          if (pint(i).ge.check) pint(i)=check
          if(delt*pint(i).lt.pcut)pint(i)=0.
        endif
      enddo
c
c***********************************************************************
c
      return
      end
