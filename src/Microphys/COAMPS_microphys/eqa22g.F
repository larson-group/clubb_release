      subroutine eqa22g(pgacwm,pgacw,tt,rho,qg,tzero
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
      real pgacwm (len)
      real pgacw  (len)
      real qg     (len)
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
      tp1=1.0/delt
c
c************************************************************************
c          compute rate pgacwm - eqa 22a RH 1984
c************************************************************************
c

      do i=1,len
c
c limit rates
c
       qlim=-qg(i)*tp1
c
       if(tt(i).gt.tzero.and.qg(i).gt.pcut)then
        pgacwm(i)=max(qlim,cwolf*(tt(i)-tzero)*pgacw(i))
       endif
c
       if(-1.0*pgacwm(i)*delt.lt.pcut)pgacwm(i)=0.
c
      enddo
c
c***********************************************************************
c
      return
      end
