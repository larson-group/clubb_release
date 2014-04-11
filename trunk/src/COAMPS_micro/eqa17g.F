      subroutine eqa17g(pgdep,qg,tt,qv,qsati,sat,icomp,termg
     1       ,appbpp,qadj,rho,tzero,delt,pcut,len)
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
      real appbpp (len)
      real delt
      real pcut
      real pgdep  (len)
      real qadj   (len)
      real qg     (len)
      real qsati  (len)
      real qv     (len)
      real rho    (len)
      real sat    (len)
      real termg  (len)
      real tt     (len)
      real tzero
c
c***********************************************************************
c          local variables and dynamic storage:
c***********************************************************************
c
      integer i
      integer nedep
      integer nevap
c
      real dqmax
      real rate
      real tp2
c
      tp2=1.0/delt
c
c***********************************************************************
c          compute ice supersaturation
c***********************************************************************
c
      do i=1,len
       sat(i)=qv(i)/qsati(i)-1.0
      enddo
c
c************************************************************************
c          compress active points for evaporating graupel
c************************************************************************
c
      nevap=0
      do i=1,len
        if (qg(i).gt.pcut.and.tt(i).le.tzero.and.sat(i).lt.0.)then
         nevap=nevap+1
         icomp(nevap)=i
        endif
      enddo
c
c************************************************************************
c          compute sublimation rate pgdep - eqa 17a RH 1984
c************************************************************************
c
      do i=1,nevap
       dqmax=-qg(icomp(i))*tp2
c
       pgdep(icomp(i))=appbpp(icomp(i))*sat(icomp(i))*termg(icomp(i))
     1                /rho(icomp(i))
c
       rate=(qv(icomp(i))-qadj(icomp(i)))*tp2
c
c************************************************************************
c          limit rates
c************************************************************************
c
       pgdep(icomp(i))=MAX(rate,dqmax,pgdep(icomp(i)) )

       if(-1.0*pgdep(icomp(i))*delt.lt.pcut)pgdep(icomp(i))=0.
c
      enddo
c
c************************************************************************
c          compress active points for vapor depostion
c************************************************************************
c
      nedep=0
      do i=1,len
        if (qg(i).gt.pcut.and.tt(i).le.tzero.and.sat(i).ge.0.)then
         nedep=nedep+1
         icomp(nedep)=i
        endif
      enddo
c
c***********************************************************************
c
      do i=1,nedep
       pgdep(icomp(i))=appbpp(icomp(i))*sat(icomp(i))*termg(icomp(i))
     1                /rho(icomp(i))
c
       rate=(qv(icomp(i))-qadj(icomp(i)))*tp2
c
c************************************************************************
c          limit rates
c************************************************************************
c
       pgdep(icomp(i))=min(rate,pgdep(icomp(i)) )
c
       if(pgdep(icomp(i))*delt.lt.pcut)pgdep(icomp(i))=0.
      enddo
c
c************************************************************************
c
      return
      end
