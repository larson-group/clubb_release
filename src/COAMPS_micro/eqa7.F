      subroutine eqa7(praut,qc,rho,ccn,pr1,pr2,pi
     1                ,rholiq,icond,delt,pcut,len)
c
#include "eqa7.prol"
c SCCS IDENTIFICATION:  %W% %G%

      implicit none
c
c***********************************************************************
c           parameters:
c***********************************************************************
c
      integer len
      integer icond
c
      real ccn    (len)
      real delt
      real pcut
      real pi
      real pr1
      real pr2
      real praut  (len)
      real qc     (len)
      real rho    (len)
      real rholiq
      real rv
c
c***********************************************************************
c          local variables and dynamic storage:
c***********************************************************************
c
      integer i
c
      real beta
      real delti
      real e1
      real e2
      real e3
      real e4
      real e5
c
c************************************************************************
c
      e1=7./3.
      e2=1./3.
      e3=2.47
      e4=-1.79
      e5=4.*pi*rholiq/3.
      beta=2.2
      delti=1/delt
c
      if(icond.eq.1)then
c
c***********************************************************************
c Rutledge and Hobbs
c***********************************************************************
c
       do i=1,len
         praut(i)=pr1*(qc(i)-pr2)
         praut(i)=(praut(i)+abs(praut(i)))*0.5
         if (delt*praut(i).lt.pcut) praut(i)=0.0
       enddo
c
      elseif(icond.eq.2)then
c
c***********************************************************************
c Manton and Cotton (1977)
c***********************************************************************
c
       do i=1,len
        if(qc(i).gt.pcut.and.ccn(i).gt.0.)then
         rv=(e5/rho(i))**(-e2)*qc(i)**e2*ccn(i)**(-e2)
         praut(i)=beta*qc(i)**e1*ccn(i)**(-e2)*(rv-pr2)
         praut(i)=(praut(i)+abs(praut(i)))*0.5
         if (delt*praut(i).lt.pcut) praut(i)=0.0
        endif
       enddo
c
c-->Addition from Jerry Schmidt
      elseif(icond.eq.3)then
c
c***********************************************************************
c Khairoutdinov and Kogan (2000) (eq. 29)
c***********************************************************************
c
       do i=1,len
        if(qc(i).gt.pcut.and.ccn(i).gt.1.)then
         praut(i)=1350.*qc(i)**2.47*ccn(i)**(-1.79)
         praut(i)=min(praut(i),qc(i)*delti)
c        print*,'praut ',praut(i)*delt,qc(i),ccn(i)
         if (delt*praut(i).lt.pcut) praut(i)=0.0
        endif
       enddo
c
c<--end of addition
c***********************************************************************
c prevent autoconversion if desired
c***********************************************************************
c
      elseif(icond.eq.4)then
       do i=1,len
         praut(i)=0.
       enddo
c
c-->Addition from Jerry Schmidt
      elseif(icond.eq.5)then
c
c***********************************************************************
c Khairoutdinov and Kogan (2000) (eq. 29) with diagnosed droplet
c        concentration
c***********************************************************************
c
       do i=1,len
        if(qc(i).gt.pcut.and.ccn(i).gt.0.)then
         praut(i)=1350.*qc(i)**2.47*ccn(i)**(-1.79)
         praut(i)=min(praut(i),qc(i)*delti)
         if (delt*praut(i).lt.pcut) praut(i)=0.0
        endif
       enddo
c
c<--end of addition
      endif
c     
c***********************************************************************
c
      return
      end
