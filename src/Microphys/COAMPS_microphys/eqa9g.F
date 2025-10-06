      subroutine eqa9g(pracs,qr,qs,rho,fallr,falls
     1       ,icomp,ssi,sri,ss,sr
     1       ,pterm4,rnzero,snzero,rhosno,pi,esr
     2       ,delt,pcut,len)
c
c#include <eqa8g.prol>
      implicit none
c
c***********************************************************************
c           parameters:
c***********************************************************************
c
      integer len
      integer icomp(len)
c
     
      real delt
      real esr
      real fallr  (len)
      real falls  (len)
      real pcut
      real pi
      real pterm4 (len)
      real pracs  (len)
      real qr     (len)
      real qs     (len)
      real rho    (len)
      real rhosno
      real rnzero
      real snzero
      real sr     (len)
      real sri    (len)
      real ss     (len)
      real ssi    (len)
c
c***********************************************************************
c          local variables and dynamic storage:
c***********************************************************************
c
      integer i
      integer npts
c
      real qlim
      real tp1
      real tp2
      real z1
      real z2
      real z3
c
c************************************************************************
c          local constants
c************************************************************************
c
      tp1=esr*pi*pi*rhosno*rnzero*snzero
      tp2=1./delt
c
c************************************************************************
c          compress active points
c************************************************************************
c
      npts=0
      do i=1,len
        if (qr(i).gt.pcut.and.qs(i).gt.pcut
     1       .and.sr(i).gt.0..and.ss(i).gt.0.) then
         npts=npts+1
         icomp(npts)=i
        endif
      enddo
c
      if(npts.eq.0)return
c
      do i=1,npts
       sri(i)=1.0/sr(icomp(i)) 
       ssi(i)=1.0/(ss(icomp(i)) )**4
      enddo
c
c************************************************************************
c          compute rate pracs - eqa8 RH 1984
c************************************************************************
c
      do i=1,npts
c
       z1 =  5./( ss(icomp(i))*ss(icomp(i)) ) 
       z2 =  2.*sri(i)/ss(icomp(i))
       z3 =  .5*sri(i)*sri(i)
       pracs(icomp(i))=z1+z2+z3
c
       pracs(icomp(i))=tp1*abs(fallr(icomp(i))-falls(icomp(i)))
     1      *pracs(icomp(i))*sri(i)*ssi(i)
     1      *pterm4(icomp(i))/rho(icomp(i))
c
c************************************************************************
c         limit rates
c************************************************************************
c
       qlim=qs(icomp(i))*tp2
       pracs(icomp(i))=min(pracs(icomp(i)),qlim)
c
       if(pracs(icomp(i))*delt.lt.pcut)pracs(icomp(i))=0.
c
      enddo
c
c***********************************************************************
c
      return
      end
