      subroutine eqa20g(pgshr,pgacr,pgacw,tt,rho,tzero
     2       ,delt,pcut,len)
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
      real delt
      real pcut
      real pgacr   (len)
      real pgacw   (len)
      real pgshr   (len)
      real rho     (len)
      real tt     (len)
      real tzero
c
c***********************************************************************
c          local variables and dynamic storage:
c***********************************************************************
c
      integer i
c
c************************************************************************
c          compute rate pgshr - eqa 20a RH84
c************************************************************************
c
      do i=1,len
c
       if(tt(i).gt.tzero)then
        pgshr(i)=pgacw(i)/rho(i)+pgacr(i)
       endif
c
       if(pgshr(i)*delt.lt.pcut)pgshr(i)=0.
c
      enddo
c
c***********************************************************************
c
      return
      end
