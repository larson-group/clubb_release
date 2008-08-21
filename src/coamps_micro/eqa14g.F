      subroutine eqa14g(pgacs,qs,qg,tt,rho,falls,fallg,ss,sg,icomp
     1       ,ssi,sggi,pterm4,snzero,gnzero,rhosno,pi,egs,delt,pcut
     2       ,tzero,len)
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
      real egs
      real fallg  (len)
      real falls  (len)
      real gnzero
      real pcut
      real pgacs  (len)
      real pi
      real pterm4 (len)
      real qg     (len)
      real qs     (len)
      real rho    (len)
      real rhosno
      real sg     (len)
      real sggi    (len)
      real snzero
      real ss     (len)
      real ssi    (len)
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
      tp1=egs*pi*pi*rhosno*snzero*gnzero
      tp2=1./delt
c
c************************************************************************
c          compress active points
c************************************************************************
c
      npts=0
      do i=1,len
        if (qs(i).gt.pcut.and.qg(i).gt.pcut
     1       .and.ss(i).gt.0..and.sg(i).gt.0..and.tt(i).lt.tzero) then
         npts=npts+1
         icomp(npts)=i
        endif
      enddo
c
      if(npts.eq.0)return
c
      do i=1,npts
       sggi(i)=1.0/sg(icomp(i)) 
       ssi(i)=1.0/(ss(icomp(i)) )**4
      enddo
c
c************************************************************************
c          compute rate pgacs - eqa14 RH84
c************************************************************************
c
      do i=1,npts
c
       z1 =  5./( ss(icomp(i))*ss(icomp(i)) ) 
       z2 =  2.*sggi(i)/ss(icomp(i))
       z3 =  .5*sggi(i)*sggi(i)
       pgacs(icomp(i))=z1+z2+z3
c
       pgacs(icomp(i))=tp1*abs(fallg(icomp(i))-falls(icomp(i)))
     1      *pgacs(icomp(i))*ssi(i)*sggi(i)
     2      *pterm4(icomp(i))/rho(icomp(i))
c
c************************************************************************
c         limit rates
c************************************************************************
c
       qlim=qs(icomp(i))*tp2
       pgacs(icomp(i))=min(pgacs(icomp(i)),qlim)
c
       if(pgacs(icomp(i))*delt.lt.pcut)pgacs(icomp(i))=0.
c
      enddo
c
c***********************************************************************
c
      return
      end
