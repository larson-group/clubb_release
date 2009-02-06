      subroutine eqa8g(psacr,qr,qs,rho,tt,fallr,falls,sri,ssi,icomp,sr
     1       ,ss,pterm4,rnzero,snzero,rholiq,pi,esr,delt,pcut,tzero,len)
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
      real psacr  (len)
      real pterm4 (len)
      real qr     (len)
      real qs     (len)
      real rho    (len)
      real rholiq
      real rnzero
      real snzero
      real sr     (len)
      real sri    (len)
      real ss     (len)
      real ssi    (len)
      real tzero
      real tt     (len)
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
      tp1=esr*pi*pi*rholiq*rnzero*snzero
      tp2=1.0/delt
c
c************************************************************************
c          compress active points
c************************************************************************
c
      npts=0
      do i=1,len
        if (qr(i).gt.pcut.and.qs(i).gt.pcut.and.tt(i).le.tzero
     1       .and.sr(i).gt.0..and.ss(i).gt.0.) then
         npts=npts+1
         icomp(npts)=i
        endif
      enddo
c
      if(npts.eq.0)return
c
      do i=1,npts
       sri(i)=1.0/(sr(icomp(i)) )**4.
       ssi(i)=1.0/ss(icomp(i)) 
      enddo
c
c************************************************************************
c          compute rate psacr - eqa8 RH 1984
c************************************************************************
c
      do i=1,npts
c
       z1 =  5./( sr(icomp(i))*sr(icomp(i)) ) 
       z2 =  2.*ssi(i)/sr(icomp(i))
       z3 =  .5*ssi(i)*ssi(i)
       psacr(icomp(i))=z1+z2+z3
c
       psacr(icomp(i))=tp1*abs(falls(icomp(i))-fallr(icomp(i)))
     1      *psacr(icomp(i))*sri(i)*ssi(i)
     1      *pterm4(icomp(i))/rho(icomp(i))
c
c************************************************************************
c          limit rates
c************************************************************************
c
       qlim=qr(icomp(i))*tp2
       psacr(icomp(i))=min(qlim,psacr(icomp(i)))

       if(psacr(icomp(i))*delt.lt.pcut)psacr(icomp(i))=0.
c
      enddo
c
c***********************************************************************
c
      return
      end
