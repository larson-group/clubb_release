      subroutine eqa12g(pgaci,qi,qg,tt,sg,sggi,icomp,pterm4,tzero
     1       ,gnzero,pi,egi,delt,pcut,abar,bgrp,len)
c
c#include <eqa10g.prol>
      implicit none
c
c***********************************************************************
c           parameters:
c***********************************************************************
c 
      integer len
      integer icomp(len)
c
      real abar
      real bgrp
      real delt
      real egi
      real garg
      real gnzero
      real pcut
      real pi
      real pgaci  (len)
      real pterm4 (len)
      real qg     (len)
      real qi     (len)
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
      real gm3pb
      real qlim
      real tp1
      real tp2
c
c************************************************************************
c          local constants
c************************************************************************
c
      tp1=pi*abar*egi*gnzero/4.0
      tp2=1.0/delt
      garg=3.0+bgrp
      call gamma(garg,gm3pb)
c
c************************************************************************
c          compress active points
c************************************************************************
c
      npts=0
      do i=1,len
        if (qg(i).gt.pcut.and.qi(i).gt.pcut.and.tt(i).le.tzero
     1       .and.sg(i).gt.0.) then
         npts=npts+1
         icomp(npts)=i
        endif
      enddo
c
      if(npts.eq.0)return
c
      do i=1,npts
       sggi(i)=1.0/(sg(icomp(i)))**(3.+bgrp)
      enddo
c
c************************************************************************
c          compute rate pgaci - RH 1984
c************************************************************************
c

      do i=1,npts
c
       pgaci(icomp(i))=tp1*pterm4(icomp(i))*qi(icomp(i))*gm3pb*sggi(i)
c
c************************************************************************
c         limit rates
c************************************************************************
c
       qlim=qi(icomp(i))*tp2
       pgaci(icomp(i))=min(pgaci(icomp(i)),qlim)
       if(pgaci(icomp(i))*delt.lt.pcut)pgaci(icomp(i))=0.
c
      enddo
c
c***********************************************************************
c
      return
      end
