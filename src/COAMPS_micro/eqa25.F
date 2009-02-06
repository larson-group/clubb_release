      subroutine eqa25(psmlt,qs,tt,terma,rho,pi,snzero,pcut
     1                 ,therco,tzero,heatlf,delt,len)
c
#include "eqa25.prol"
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
      real heatlf
      real pcut
      real pi
      real psmlt  (len)
      real qs     (len)
      real rho    (len)
      real snzero
      real terma  (len)
      real therco
      real tt     (len)
      real tzero
c
c***********************************************************************
c          local variables and dynamic storage:
c***********************************************************************
c
      integer i
c
      real tp1
      real tp2
c
c************************************************************************
c          local constants
c************************************************************************
c
      tp1=2.0*pi*snzero*therco/heatlf
      tp2=1.0/delt
c
      do i=1,len
        if (tt(i).ge.tzero.and.qs(i).gt.pcut) then
          psmlt(i)=-tp1*(tt(i)-tzero)*terma(i)/rho(i)
          psmlt(i)=MAX(psmlt(i),-qs(i)*tp2)
          if(abs(delt*psmlt(i)).lt.pcut)psmlt(i)=0.
        endif
      enddo
c
c**********************************************************************
c
      return
      end
