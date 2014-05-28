      subroutine tgqi(fallv,qi,rbm,pr,conci,pcut,p00,cimass,dmax,len)

#include "tvqs.prol"
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
      real cimass
      real conci  (len)
      real dmax
      real fallv  (len)
      real p00
      real pcut
      real pr     (len)
      real qi     (len)
      real rbm    (len)
c
c***********************************************************************
c          local variables and dynamic storage:
c***********************************************************************
c
      integer i
c
      real avgmas
      real ex1
      real ex2
      real falspd
      real k1
      real k2
      real pfact
      real x1
      real x6
c
      real tp1
      real dbari
c
c***********************************************************************
c          local constants
c***********************************************************************
c
      ex1    = 0.5
      ex2    = 0.0889
      tp1    = 16.3
c
      do i=1,len
c
c***********************************************************************
c          if the mixing ratio of pristine ice does not exeed cutting level
c           then set fall velocity to zero
c
c           special check: if pristine ice exist in a region where 
c                          the diagnosed value is zero then
c                          diagnose crystal concentration
c                          from avgerage ice mass assuming max crystal size
c                          of 1 cm.
c***********************************************************************
c

          if (qi(i).gt.pcut) then
           avgmas=qi(i)*rbm(i)/max(.01,conci(i))
           dbari=tp1*sqrt(avgmas)
           if(dbari.gt.dmax)then
            conci(i)=qi(i)*rbm(i)/(dmax/tp1)**2
           endif
          endif

c
c***********************************************************************
c          if mixing ratio of pristine ice exceeds cutting level
c          calculate fall velocity (eq 4.62 Cotton and Anthes (1989)
c          note:  For COAMPS change mean crystal mass (x1) from kg to g 
c***********************************************************************
c
          if (qi(i).le.pcut.or.conci(i).le.0.) then
            fallv(i)=0.0
          else
            pfact=exp(ex1*log(p00/pr(i)))
            x1=1000.*rbm(i)*qi(i)/conci(i)
            if(x1.le.1.e-7)then
             k1=0.515
             k2=304.
             x6=k1*sqrt(x1)
             fallv(i)=-k2*x6*pfact
            elseif(x1.le.1.e-5.and.x1.gt.1.e-7)then
             k1=0.192
             k2=1250.
             x6=k1*sqrt(x1)
             fallv(i)=-k2*x6*pfact
            elseif(x1.gt.1.e-5)then
             pfact=exp(ex2*log(p00/pr(i)))
             k1=.0889
             k2=4.84
             x6=k1*(x1)**(.417)
             fallv(i)=-k2*x6**.25*pfact
            endif
          endif
c
      enddo
c
      return
      end
