      subroutine slope
     1  (sr,srcm,ss,sg,qr,qs,qg,rho,sloper,slopes
     1   ,slopeg,pcut,lice,lgrpl,len)

#include "slope.prol"
c SCCS IDENTIFICATION:  %W% %G%
c
c***********************************************************************
c
      implicit none
c
c***********************************************************************
c
      integer len
c
      logical lice
      logical lgrpl
c
      real pcut
      real qg     (len)
      real qr     (len)
      real qs     (len)
      real rho    (len)
      real slopeg
      real sloper
      real slopes
      real sg     (len)
      real sr     (len)
      real srcm   (len)
      real ss     (len)
c
c***********************************************************************
c local variables
c***********************************************************************
c
      integer i
      real cm2m
c
      cm2m=100.0
c
c***********************************************************************
c          compute slope factor for rain (rs)
c***********************************************************************
c
      do i=1,len
        if (qr(i).le.pcut) then
          sr(i)=-2.0
          srcm(i)=-2.0
        else
          sr(i)=sloper/(rho(i)*qr(i))
          srcm(i)=sqrt(sqrt(sr(i)))
          sr(i)=srcm(i)*cm2m
        endif
      enddo
c
c***********************************************************************
c          compute slope factor for snow (ss)
c***********************************************************************
c
      if (lice) then
c
        do i=1,len
          if (qs(i).le.pcut) then
            ss(i)=-2.0
          else
            ss(i)=slopes/(rho(i)*qs(i))
            ss(i)=sqrt(sqrt(ss(i)))*cm2m
          endif
        enddo
c
        if(lgrpl)then
         do i=1,len
          if (qg(i).le.pcut) then
            sg(i)=-2.0
          else
            sg(i)=slopeg/(rho(i)*qg(i))
            sg(i)=sqrt(sqrt(sg(i)))*cm2m
          endif
         enddo
        endif
c
      else
c
        do i=1,len
          ss(i)=-2.0
          sg(i)=-2.0
        enddo
c
      endif
c
c***********************************************************************
c
      return
      end
