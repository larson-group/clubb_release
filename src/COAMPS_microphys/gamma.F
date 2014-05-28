      subroutine gamma(arg,val)

#include "gamma.prol"
c SCCS IDENTIFICATION:  %W% %G%
      implicit none
c
c***********************************************************************
c           parameters:
c***********************************************************************
c
      real arg
      real val
c
c***********************************************************************
c          local variables and dynamic storage:
c***********************************************************************
c
      integer k
      real c      (15)
      real fac
      real wamma
      real zarg
c
c***********************************************************************
c          gamma function of z
c***********************************************************************
c
      data c/1.000000,0.577210,-0.655870,-0.042000,0.166530
     1     ,-0.042190,-0.009620, 0.007210,-0.001160,-0.000210
     2     , 0.000120,-0.000020,-0.000001, 0.000001, 2.e-7   /
c
c
      zarg=arg
      fac=1.0
c
      if (zarg.le.2.0) go to 105
  100 continue
c
      zarg=zarg-1.0
      fac=fac*zarg
c
      if (zarg.gt.2.0) go to 100
  105 continue
c
      wamma=0.0
c
      do k=1,15
        !print *, k, wamma, c(k), zarg**k
        wamma=wamma+c(k)*zarg**k
      enddo
c
      val=fac/wamma

c***********************************************************************
c
      return
      end subroutine gamma
