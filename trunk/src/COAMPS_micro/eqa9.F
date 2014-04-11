      subroutine eqa9(pracw,qc,qr,pterm4,sr,sri,rho,nc,icomp,pi,erc
     1               ,rnzero,aa0,aa1,aa2,aa3,gm3,gm4,gm5,gm6
     2               ,pcut,delt,ldrizzle,len)
c
#include "eqa9.prol"
c SCCS IDENTIFICATION:  %W% %G%

      implicit none
c
c***********************************************************************
c           parameters:
c***********************************************************************
c
      integer len
      integer icomp(len)
c
      real aa0
      real aa1
      real aa2
      real aa3
      real delt
      real erc
      real gm3
      real gm4
      real gm5
      real gm6
      real nc     (len)
      real pcut
      real pi
      real pracw  (len)
      real pterm4 (len)
      real qc     (len)
      real qr     (len)
      real rho    (len)
      real rnzero
      real rvc
      real rvcof
      real sr     (len)
      real sri    (len)
c
      logical ldrizzle
c
c***********************************************************************
c          local variables and dynamic storage:
c***********************************************************************
c
      integer i
      integer npts
c
      real a0gm3
      real a1gm4
      real a2gm5
      real a3gm6
      real qlim
      real tp1
      real tp2
      real tp4
      real z1
      real z2
      real z3
      real z4
c
c************************************************************************
c          local constants
c************************************************************************
c
c tp4=1./(den_water*pi*4/3)
c
      tp1=pi*erc*rnzero/4.0
      tp2=1.0/delt
      tp4=2.3873241e-4
      a0gm3=aa0*gm3
      a1gm4=aa1*gm4
      a2gm5=aa2*gm5
      a3gm6=aa3*gm6
c
      npts=0
      do i=1,len
       if (qr(i).ge.pcut.and.qc(i).ge.pcut.and.sr(i)
     1    .gt.0.) then
        npts=npts+1
        icomp(npts)=i
       endif
      enddo
      if(npts.eq.0)return

      if(.not.ldrizzle)then
c
c************************************************************************
c  RH83 accretion
c************************************************************************
c
       do i=1,npts
        sri(i)=1.0/( sr(icomp(i))*sr(icomp(i))*sr(icomp(i)) )
       enddo
c
       do i=1,npts
        z1=a0gm3*sri(i)
        z2=a1gm4*sri(i)/  sr(icomp(i))
        z3=a2gm5*sri(i)/( sr(icomp(i))*sr(icomp(i)) )
        z4=a3gm6*sri(i)*sri(i)
        pracw(icomp(i))=z1+z2+z3+z4
       enddo
c
       do i=1,npts
        qlim=qc(icomp(i))*tp2
        pracw(icomp(i))=pracw(icomp(i))*pterm4(icomp(i))
     1     *qc(icomp(i))*tp1
c
c************************************************************************
c    limit rates
c************************************************************************
c
        pracw(icomp(i))=min(pracw(icomp(i)),qlim)
c
        if (pracw(icomp(i))*delt.le.pcut) pracw(icomp(i))=0.0
       enddo
c
c-->Addition from Jerry Schmidt
      else
c
c************************************************************************
c     Khairoutdinov and Kogan (2000) drizzle parameterization eq(33)
c************************************************************************
c
       do i=1,npts
        qlim=qc(icomp(i))*tp2
        if(nc(icomp(i)).gt.0.)then
         pracw(icomp(i)) = 67. * (qc(icomp(i))*qr(icomp(i)))**1.15
         pracw(icomp(i)) = min(pracw(icomp(i)), qlim )
         if (pracw(icomp(i))*delt.le.pcut) pracw(icomp(i))=0.0
        else
         pracw(icomp(i)) = 0.
        endif
       enddo
c
c<--end of addition
      endif
c
c***********************************************************************
c
      return
      end
