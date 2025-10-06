      subroutine eqa5g(praci,qi,qr,tt,icomp,gm3,gm4,gm5,gm6,pterm4
     1       ,rnzero,pi,eri,a0,a1,a2,a3,sr,delt,pcut,tzero,len)
c
c#include <eqa5g.prol>
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
      real delt
      real eri
      real gm3
      real gm4
      real gm5
      real gm6
      real pcut
      real pi
      real praci  (len)
      real pterm4 (len)
      real qi     (len)
      real qr     (len)
      real rnzero
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
      real a0gm3
      real a1gm4
      real a2gm5
      real a3gm6
      real qlim
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
      tp1=eri*.25*pi*rnzero
      tp2=1.0/delt
      a0gm3=a0*gm3
      a1gm4=a1*gm4
      a2gm5=a2*gm5
      a3gm6=a3*gm6
c
c************************************************************************
c          compress active points
c************************************************************************
c
      npts=0
      do i=1,len
        if (qi(i).gt.pcut.and.qr(i).gt.pcut
     1       .and.sr(i).gt.0..and.tt(i).lt.tzero) then
         npts=npts+1
         icomp(npts)=i
        endif
      enddo
c
      if(npts.eq.0)return
c
      do i=1,npts
       sri(i)=1.0/(sr(icomp(i)))**3
      enddo
c
c************************************************************************
c          compute rate praci eqa5 - RH 1984
c************************************************************************
c
      do i=1,npts
c
       z1 =  a0gm3*sri(i)
       z2 =  a1gm4*sri(i)/ sr(icomp(i))
       z3 =  a2gm5*sri(i)/(sr(icomp(i))*sr(icomp(i)))
       z4 =  a3gm6*sri(i)/( sr(icomp(i))*sr(icomp(i))*sr(icomp(i)) )
       praci(icomp(i))=z1+z2+z3+z4
c
       praci(icomp(i))=tp1*qi(icomp(i))*praci(icomp(i))
     1                 *pterm4(icomp(i))
c
c************************************************************************
c limit rates
c************************************************************************
c
       qlim=qi(icomp(i))*tp2
       praci(icomp(i))=min(qlim,praci(icomp(i)))
       if(praci(icomp(i))*delt.lt.pcut)praci(icomp(i))=0.
c
      enddo
c
c***********************************************************************
c
      return
      end
