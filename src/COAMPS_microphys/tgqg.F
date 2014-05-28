      subroutine tgqg(fallv,qg,sg,pterm4,pcut,abar,bgrp,len)

#include "tvqs.prol"
c
      implicit none
c
c***********************************************************************
c           parameters:
c***********************************************************************
c
      integer len
c
      real fallv  (len)
      real pcut
      real pterm4 (len)
      real qg     (len)
      real sg     (len)
c
c***********************************************************************
c          local variables and dynamic storage:
c***********************************************************************
c
      integer i
c
      real abar
      real bgrp
      real ex1
      real falspd
      real garg
      real pfact
      real tvs1
      real x6
c
c***********************************************************************
c          local constants
c***********************************************************************
c
      ex1=0.4
      garg=4.0+bgrp
      call gamma(garg,tvs1)
      tvs1=abar*tvs1/6.0
c
c***********************************************************************
c          compute terminal velocity of graupel eqa a3 RH 1984
c***********************************************************************
c
      do i=1,len
c
c***********************************************************************
c          if the mixing ratio of graupel does not exeed cutting level
c            then set fall velocity to zero
c***********************************************************************
c
        if (qg(i).le.pcut.or.sg(i).lt.0.) then
          fallv(i)=0.0
c
c***********************************************************************
c          else if mixing ratio of graupel exceeds cutting level
c            calculate fall velocity
c***********************************************************************
c
        else
          falspd=exp(-bgrp*log(sg(i)) )
          fallv(i)=-falspd*pterm4(i)*tvs1
        endif
      enddo
c
c***********************************************************************
c
      return
      end
