      subroutine eqa7g(piacr,qi,qr,rho,tt,sr,icomp,conci
     1       ,pterm4,gm6,gm7,gm8,gm9,rnzero,rholiq,pi,eri
     1       ,a0,a1,a2,a3,delt,pcut,tzero,len)

c
c#include <eqa7g.prol>
      implicit none
c
c***********************************************************************
c           parameters:
c***********************************************************************
c
      integer len
      integer icomp(len)
c
     
      real a0
      real a1
      real a2
      real a3
      real conci  (len)
      real delt
      real eri
      real gm6
      real gm7
      real gm8
      real gm9
      real rnzero
      real pcut
      real pi
      real pterm4 (len)
      real piacr  (len)
      real qi     (len)
      real qr     (len)
      real rho    (len)
      real rholiq
      real sr     (len)
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
      real a0gm6
      real a1gm7
      real a2gm8
      real a3gm9
      real qlimr
      real sri(len)
      real tp1
      real tp2
      real z1
      real z2
      real z3
      real z4
c
c************************************************************************
c          local constants
c************************************************************************
c
      tp1=eri*pi*pi*rholiq*rnzero/24.0
      tp2=1.0/delt
      a0gm6=a0*gm6
      a1gm7=a1*gm7
      a2gm8=a2*gm8
      a3gm9=a3*gm9
c
c************************************************************************
c          compress active points
c************************************************************************
c
      npts=0
      do i=1,len
        if (qi(i).gt.pcut.and.qr(i).gt.pcut.and.tt(i).le.tzero
     1       .and.sr(i).gt.0.) then
         npts=npts+1
         icomp(npts)=i
        endif
      enddo
c
      if(npts.eq.0)return
c
      do i=1,npts
       sri(i)=1.0/(sr(icomp(i)))**6
      enddo
c
c************************************************************************
c          compute rate piacr
c************************************************************************
c
      do i=1,npts
c
       z1 =  a0gm6*sri(i)
       z2 =  a1gm7*sri(i)/ sr(icomp(i))
       z3 =  a2gm8*sri(i)/(sr(icomp(i))*sr(icomp(i)))
       z4 =  a3gm9*sri(i)/( sr(icomp(i))*sr(icomp(i))*sr(icomp(i)) )
       piacr(icomp(i))=z1+z2+z3+z4
c
       piacr(icomp(i))=tp1*conci(icomp(i))*piacr(icomp(i))
     1       *pterm4(icomp(i))/rho(icomp(i))
c
c************************************************************************
c          limit rates
c************************************************************************
c
       qlimr=qr(icomp(i))*tp2
       piacr(icomp(i))=max(0.,min(piacr(icomp(i)),qlimr))
c
       if(piacr(icomp(i))*delt.lt.pcut)piacr(icomp(i))=0.
c
      enddo
c
c***********************************************************************
c
      return
      end
