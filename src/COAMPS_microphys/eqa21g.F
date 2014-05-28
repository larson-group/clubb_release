      subroutine eqa21g(pgacrm,pgacr,tt,rho,qg,tzero
     2       ,delt,pcut,heatlf,cw,len)
c
c#include <eqa10g.prol>
      implicit none
c
c***********************************************************************
c           parameters:
c***********************************************************************
c
      integer len
c
      real cw
      real delt
      real heatlf
      real pcut
      real pgacrm  (len)
      real pgacr   (len)
      real qg      (len)
      real rho    (len)
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
      real cwolf
      real qlim
      real tp1
c
      cwolf=-cw/heatlf
c
c************************************************************************
c          compute rate pgacrm - eqa 21a RH 1984
c************************************************************************
c
      tp1=1.0/delt

      do i=1,len
c
c************************************************************************
c limit rate
c************************************************************************
c
       qlim=-qg(i)*tp1
c
       if(tt(i).gt.tzero.and.qg(i).gt.pcut)then
        pgacrm(i)=max(qlim,cwolf*(tt(i)-tzero)*pgacr(i))
       endif
c
       if(-1.0*pgacrm(i)*delt.lt.pcut)pgacrm(i)=0.
c
      enddo
c
c***********************************************************************
c
      return
      end
