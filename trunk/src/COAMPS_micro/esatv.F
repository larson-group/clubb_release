      subroutine esatv(t,w1,m)
c
#include "esatv.prol"
c SCCS IDENTIFICATION:  %W% %G%
      implicit none
c
c***********************************************************************
c           parameters:
c***********************************************************************
c
      integer m
c
      real t     (m)
      real w1    (m)
c
c***********************************************************************
c          local variables and dynamic storage:
c***********************************************************************
c
      integer i
      integer i1     (m)
c
      real epi
      real epv
      real esi
      real esv
      real w2     (m)
      real w3     (m)
c
      logical l_new
      logical l_ice
c
      data l_new /.true./
c
c***********************************************************************
c
      common /vap/esv(375),epv(375)
      common /ice/esi(375),epi(375)
c
c
      if (l_new) then
c
c       compute sat vapor pressure (100*mb or mb?) over water using 
c        Smithsonian table, based on GOFF-GRATCH formulation, 1946
c
        l_ice=.false.
        call esat_new(t,w1,m,l_ice)
c
      else
c
c       compute sat vapor pressure over water using Teton's
c        approximate formula 
c
      do i=1,m
        i1(i)=t(i)
      enddo
c
      do i=1,m
        w3(i)=t(i)-i1(i)
      enddo
c
      do i=1,m
        w1(i)=esv(i1(i))
        w2(i)=esv(i1(i)+1)
      enddo
c
      do i=1,m
        w2(i)=(w2(i)-w1(i))*w3(i)
      enddo
c
      do i=1,m
        w1(i)=w1(i)+w2(i)
      enddo
c
      endif
c
c***********************************************************************
c
      return
      end
