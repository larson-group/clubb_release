      subroutine eqa13g(pgacr,qr,qg,tt,rho,fallr,fallg
     1       ,sr,sg,icomp,sri,sggi
     1       ,pterm4,rnzero,gnzero,rholiq,pi,egr
     2       ,delt,pcut,tzero,len)
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
      real egr
      real fallr  (len)
      real fallg  (len)
      real gnzero
      real pcut
      real pi
      real pterm4 (len)
      real pgacr  (len)
      real qr     (len)
      real qg     (len)
      real rho    (len)
      real rholiq
      real rnzero
      real sr     (len)
      real sri    (len)
      real sg     (len)
      real sggi    (len)
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
      tp1=egr*pi*pi*rholiq*rnzero*gnzero
      tp2=1.0/delt
c
c************************************************************************
c          compress active points
c************************************************************************
c
      npts=0
      do i=1,len
        if (qr(i).gt.pcut.and.qg(i).gt.pcut
     1       .and.sr(i).gt.0..and.sg(i).gt.0.) then
         npts=npts+1
         icomp(npts)=i
        endif
      enddo
c
      if(npts.eq.0)return
c
      do i=1,npts
       sggi(i)=1.0/sg(icomp(i)) 
       sri(i)=1.0/(sr(icomp(i)) )**4
      enddo
c
c************************************************************************
c          compute rate pgacr - eqa8 RH 1984
c************************************************************************
c
      do i=1,npts
c
       z1 =  5./( sr(icomp(i))*sr(icomp(i)) ) 
       z2 =  2.*sggi(i)/sr(icomp(i))
       z3 =  .5*sggi(i)*sggi(i)
       pgacr(icomp(i))=z1+z2+z3
c
       pgacr(icomp(i))=tp1*abs(fallg(icomp(i))-fallr(icomp(i)))
     1      *pgacr(icomp(i))*sri(i)*sggi(i)
     1      *pterm4(icomp(i))/rho(icomp(i))
c
c************************************************************************
c          limit rates
c************************************************************************
c
       qlim=qr(icomp(i))*tp2
       pgacr(icomp(i))=min(pgacr(icomp(i)),qlim)

       if(pgacr(icomp(i))*delt.lt.pcut)pgacr(icomp(i))=0.
c
      enddo
c
c***********************************************************************
c

      return
      end
