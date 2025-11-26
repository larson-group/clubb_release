      subroutine eqa21(psaci,qi,qs,ss,termb,esi,pcut,twodt,len)
c
#include "eqa21.prol"
c SCCS IDENTIFICATION:  %W% %G%

      implicit none
c
c***********************************************************************
c           parameters:
c***********************************************************************
c
      integer len
c
      real esi
      real pcut
      real psaci  (len)
      real qi     (len)
      real qs     (len)
      real ss     (len)
      real termb  (len)
      real twodt
c
c***********************************************************************
c          local variables and dynamic storage:
c***********************************************************************
c
      integer i
c
c***********************************************************************
c     compute psaci: eqa21 RH83
c***********************************************************************
c
      do i=1,len
        if (ss(i).ge.0.0) then
          if (qs(i).ge.pcut.and.qi(i).ge.pcut)then
             psaci(i)=qi(i)*esi*termb(i)
             if(psaci(i)*twodt.lt.pcut)psaci(i)=0.0
          endif
        endif
      enddo
c
c***********************************************************************
c
      return
      end
