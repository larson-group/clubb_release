      subroutine eqa6(pcond,tt,qadj,qv,qc,qsatv,qsati,pp,w1,w2
     1               ,iwx1,heatlv,cp,delt,pcut,j,len)
c
#include "eqa6.prol"
c SCCS IDENTIFICATION:  %W% %G%

c***********************************************************************
c
      implicit none
c
c***********************************************************************
c           parameters:
c***********************************************************************
c
      integer j
      integer len
      integer iwx1   (len)
c
      real cp
      real delt
      real heatlv
      real pcond  (len)
      real pcut
      real pp     (len)
      real qadj   (len)
      real qc     (len)
      real qsati  (len)
      real qsatv  (len)
      real qv     (len)
      real tt     (len)
      real w1     (len)
      real w2     (len)
c
c***********************************************************************
c          local variables and dynamic storage:
c***********************************************************************
c
      integer i
      integer iter
      integer nevap
      integer nsat
c
      real tp2
c
c***********************************************************************
c          local constants
c***********************************************************************
c
      tp2=1.0/delt
c
c***********************************************************************
c***********************************************************************
c        evaporation of cloud droplets
c        count number of points that are below supersaturation
c        and have clouds
c***********************************************************************
c***********************************************************************
c
      nevap=0
      do i=1,len
        if (qv(i).lt.qsatv(i).and.qc(i).gt.pcut) then
          nevap=nevap+1
          iwx1(nevap)=i
        endif
      enddo
c
c***********************************************************************
c     Compute cloud evaporation. 
c     Don't allow evaporation rate to exceed available cloud water.
c     Limit vapor to  (qv-qadj)/twodt per timestep
c***********************************************************************
c
c***********************************************************************
c
      if (nevap.gt.0) then
c
        do i=1,nevap
          w1(i)=qc(iwx1(i))*tp2
          w2(i)=(qv(iwx1(i))-qadj(iwx1(i)))*tp2
        enddo

        do i=1,nevap
          pcond(iwx1(i))=max(-w1(i),w2(i))
          if(-1.0*delt*pcond(iwx1(i)).lt.pcut)then
           pcond(iwx1(i))=0.
          endif
        enddo
c
      endif
c
c***********************************************************************
c***********************************************************************
c        condensational heating
c***********************************************************************
c***********************************************************************
c
c***********************************************************************
c        count number of supersaturated points
c***********************************************************************
c
      nsat=0
      do i=1,len
        if (qv(i).gt.qsatv(i))then
          nsat=nsat+1
          iwx1(nsat)=i
        endif
      enddo
c
c***********************************************************************
c        compute heating rate and expand points to full grid
c        Limit vapor to  (qv-qadj)/twodt per timestep
c***********************************************************************
c
      if (nsat.gt.0) then
        do i=1,nsat
          pcond(iwx1(i))=(qv(iwx1(i))-qadj(iwx1(i)))*tp2
          if(delt*pcond(iwx1(i)).lt.pcut)pcond(iwx1(i))=0.
        enddo
      endif
c
c***********************************************************************
c
      return
      end
