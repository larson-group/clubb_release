      subroutine tgqs(fallv,qs,ss,pterm4,pcut,aprpr,bsnow,len)

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
      real qs     (len)
      real ss     (len)
c
c***********************************************************************
c          local variables and dynamic storage:
c***********************************************************************
c
      integer i
c
      real cm2m
      real aprpr
      real bsnow
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
! Michael Falk removed, 12 October 2006, because this line shouldn't be here.
! bsnow is passed in.
!      bsnow=0.11
! eMFc
      ex1=0.4
      garg=4.0+bsnow
      call gamma(garg,tvs1)
      tvs1=aprpr*tvs1/6.0
      cm2m=100.
c
c***********************************************************************
c          compute terminal velocity of snowflakes
c***********************************************************************
c
c***********************************************************************
c          loop over every grid point at each level
c***********************************************************************
c
      do i=1,len
c
c***********************************************************************
c          if the mixing ratio of snowflakes does not exeed cutting level
c            then set fall velocity to zero
c***********************************************************************
c
        if (qs(i).le.pcut) then
          fallv(i)=0.0
c
c***********************************************************************
c          else if mixing ratio of snowflakes exceeds cutting level
c            calculate fall velocity
c***********************************************************************
c
        else
c
c**************************************************************
c Modification by Adam Smith, 17 April 2008
c Updating fall speed equation to match the one used in COAMPS
c The original version may have had an error in the units.
c**************************************************************
c
!          falspd=exp(-bsnow*log(cm2m*sqrt(ss(i))))
          falspd=exp(-bsnow*log(ss(i)))
          fallv(i)=-falspd*pterm4(i)*tvs1
        endif
      enddo
c
c***********************************************************************
c
      return
      end
