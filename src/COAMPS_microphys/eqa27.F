      subroutine eqa27(psmlte,qs,qv,qsatv,qadj,tt,terma,apbp,rho
     1                ,snzero,tzero,pcut,delt,j,len)
c
#include "eqa27.prol"
c SCCS IDENTIFICATION:  %W% %G%

      implicit none
c
c***********************************************************************
c           parameters:
c***********************************************************************
c
      integer len
c
      integer j
c
      real apbp   (len)
      real delt  
      real pcut
      real psmlte (len)
      real qadj   (len)
      real qs     (len)
      real qsatv  (len)
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
c***********************************************************************
c          local constants
c***********************************************************************
c
      tp1=4.0*snzero
      tp2=1.0/delt
c
      do i=1,len
          if (qs(i).gt.pcut.and.tt(i).ge.tzero) then
            psmlte(i)=tp1*(qv(i)/qsatv(i)-1.0)*apbp(i)*terma(i)/rho(i)
            rate=(qv(i)-qadj(i))*tp2
c
c************************************************************************
c          limit rates
c************************************************************************
c
          if (psmlte(i).gt.0.0.and.psmlte(i).ge.rate) psmlte(i)=rate
          if (psmlte(i).lt.0.0)
     1         psmlte(i)=MAX(psmlte(i),-qs(i)*tp2,rate)
          if(abs(delt*psmlte(i)).lt.pcut)psmlte(i)=0.

          endif
      enddo
c
c***********************************************************************
c
      return
      end
