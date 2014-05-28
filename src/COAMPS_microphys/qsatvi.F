      subroutine qsatvi(t,p,qvs,qvi,lice,len)

#include "qsatvi.prol"
c SCCS IDENTIFICATION:  %W% %G%
c
      implicit none
c
c***********************************************************************
c           parameters:
c***********************************************************************
c
      integer len
c
      logical lice
c
      real p      (len)
      real qvi    (len)
      real qvs    (len)
      real t      (len)
c
c***********************************************************************
c          local variables and dynamic storage:
c***********************************************************************
c
      integer i
c
c***********************************************************************
c          compute saturation mixing ratio over water
c***********************************************************************
c
      call esatv(t,qvs,len)
c
c***********************************************************************
c          saturation mixing ratio over ice
c***********************************************************************
c
      if (lice) then
        call esati(t,qvi,len)
        do i=1,len
          qvs(i)=0.62197*qvs(i)/(p(i)-qvs(i))
          qvi(i)=0.62197*qvi(i)/(p(i)-qvi(i))
        enddo
      else
        do i=1,len
          qvs(i)=0.62197*qvs(i)/(p(i)-qvs(i))
          qvi(i)=0.0
        enddo
      endif
c
c***********************************************************************
c
      return
      end
