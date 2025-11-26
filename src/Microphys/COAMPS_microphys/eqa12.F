      subroutine eqa12(prevp,qr,qv,qadj,qsatv,rho,tt,sr,pt,apbp
     1                ,nr,pr,w1,w2,w3,w4,w5,w6,w7,iwx1,pi,rnzero,apr
     2                ,pcut,visair,gm3,delt,ldrizzle,len)
c
#include "eqa12.prol"
c SCCS IDENTIFICATION:  %W% %G%

      implicit none
c
c***********************************************************************
c           parameters:
c***********************************************************************
c
      integer len
      integer iwx1(len)
c
      real apbp   (len)
      real apr
      real delt
      real gm3
      real nr     (len)
      real pcut
      real pi
      real pr     (len)
      real prevp  (len)
      real pt     (len)
      real qadj   (len)
      real qr     (len)
      real qsatv  (len)
      real qv     (len)
      real rho    (len)
      real rnzero
      real sr     (len)
      real tt     (len)
      real twodt
      real visair
      real w1     (len)
      real w2     (len)
      real w3     (len)
      real w4     (len)
      real w5     (len)
      real w6     (len)
      real w7     (len)
c
      logical ldrizzle
c
c***********************************************************************
c          local variables and dynamic storage:
c***********************************************************************
c
      integer i
      integer nevap
c
      real c
      real qcond
      real qevap
      real thomo
      real tp1
      real tp2
      real tp3
      real tp4
      real tp5
      real tp6
      real tp7
      real w8     (len)
      real z1
      real z2
      real z3
c
c
c************************************************************************
c          local constants
c************************************************************************
c
      thomo=233.15
      tp1=0.78
      tp2=0.31*gm3
      tp3=apr/visair
      tp4=2.0*pi*rnzero
      tp5=1.0/delt
c
c************************************************************************
c density of water*pi*4/3=4.1887902e3
c************************************************************************
c
      tp7=4.1887902e3
c
c***********************************************************************
c***********************************************************************
c        evaporation of rain droplets
c        count number of points that are below supersaturation
c        and have rain water 
c***********************************************************************
c***********************************************************************
c
      nevap=0
      do i=1,len
        if (qv(i).lt.qsatv(i).and.qr(i).gt.pcut.and.sr(i).gt.0.) then
          nevap=nevap+1
          iwx1(nevap)=i
        endif
      enddo
c
      if(nevap.gt.0)then
c
        do i=1,nevap
          w1(i)=qv(iwx1(i))/qsatv(iwx1(i))-1.0
          w2(i)=sr(iwx1(i))
          w3(i)=pt(iwx1(i))
          w4(i)=qr(iwx1(i))*tp5
          w5(i)=(qv(iwx1(i))-qadj(iwx1(i)))*tp5
          w6(i)=rho(iwx1(i))
          w7(i)=1./(sr(iwx1(i))*sr(iwx1(i)))
          w8(i)=apbp(iwx1(i))
        enddo
c
c***********************************************************************
c     Don't allow evaporation rate to exceed availble rain content
c     (w4) and limit vapor tendency to  (qv-qadj)/twodt per timestep (w5)
c***********************************************************************
c
       if(.not.ldrizzle)then
c
c***********************************************************************
c     RH83 evaporation
c***********************************************************************
c
        do i=1,nevap
            z1=tp1*w7(i)
            z2=tp2*w3(i)*sqrt(tp3*w6(i))*w7(i)/w2(i)
            z3=tp4*w1(i)*(z1+z2)*w8(i)/w6(i)
            prevp(iwx1(i))=MAX(z3,-w4(i),w5(i))
            if(abs(delt*prevp(iwx1(i))).lt.pcut)prevp(iwx1(i))=0.
        enddo
c
c-->Addition from Jerry Schmidt
       else
c
c***********************************************************************
c       Khairoutdinov and Kogan (2000) drizzle parameterization (eq. 22)
c***********************************************************************
c
c       tp6=Fk+Fd  Refer to: A Short Course in Cloud Physics
c                            Rogers and Yau, 1991, page 102
c       c: condensation coef. for drizzle.
c       c=3.*(4./3.*pi*den_water/den_air)**(2./3.)*0.866/(Fk+Fd)
c       1.27=2.042*0.622
c       note: Though the definition of the supersaturation differs here
c             from the original Mechem Code, there is no  minus sign
c             needed in the equation for prevp. This stems from the different
c             treatment of the sign used in the tendency equation for rain
c             used in the original code obtained from OU (subroutine devap)
c check on units of pr
c***********************************************************************
c
        do i=1,nevap
         tp6=.5543903e11/tt(iwx1(i))**2+
     1      1.27e7*tt(iwx1(i))/(qsatv(iwx1(i))*pr(iwx1(i)) )
c
         c = 2.598*(tp7/w6(i))**0.666667/tp6
c
         prevp(iwx1(i)) = c * w1(i)*
     1      ( nr(iwx1(i))*nr(iwx1(i))*qr(iwx1(i)) )**0.3333333
c
         prevp(iwx1(i))=max(w5(i),prevp(iwx1(i)))
         if(abs(delt*prevp(iwx1(i))).lt.pcut)prevp(iwx1(i))=0.
        enddo
c<--end of addition
       endif
c
      endif
c
c***********************************************************************
c***********************************************************************
c        condensation of rain droplets
c        count number of points that are supersaturated
c        and have rain water 
c        Limit condensation to those points greater than thomo
c***********************************************************************
c***********************************************************************
c
      nevap=0
      do i=1,len
        if (qv(i).gt.qsatv(i).and.qr(i).gt.pcut.and.tt(i)
     1     .gt.thomo.and.sr(i).gt.0.) then
          nevap=nevap+1
          iwx1(nevap)=i
        endif
      enddo
c
      if(nevap.gt.0)then
c
        do i=1,nevap
          w1(i)=qv(iwx1(i))/qsatv(iwx1(i))-1.0
          w2(i)=sr(iwx1(i))
          w3(i)=pt(iwx1(i))
          w4(i)=qr(iwx1(i))*tp5
          w5(i)=(qv(iwx1(i))-qadj(iwx1(i)))*tp5
          w6(i)=rho(iwx1(i))
          w7(i)=1./(sr(iwx1(i))*sr(iwx1(i)))
          w8(i)=apbp(iwx1(i))
        enddo
c
c***********************************************************************
c     Limit vapor tendency to  (qv-qadj)/twodt per timestep (w5)
c***********************************************************************
c
       if(.not.ldrizzle)then
c
c***********************************************************************
c       RH83 rain condensation term
c***********************************************************************
c
        do i=1,nevap
            z1=tp1*w7(i)
            z2=tp2*w3(i)*sqrt(tp3*w6(i))*w7(i)/w2(i)
            z3=tp4*w1(i)*(z1+z2)*w8(i)/w6(i)
            prevp(iwx1(i))=min(z3,w5(i))
            if(delt*prevp(iwx1(i)).lt.pcut)prevp(iwx1(i))=0.
        enddo
c
c-->Addition from Jerry Schmidt
       else
c
c***********************************************************************
c       Khairoutdinov and Kogan (2000) drizzle parameterization (eq. 22)
c***********************************************************************
c
c       check on units of pr
c***********************************************************************
c
        do i=1,nevap
         tp6=.5543903e11/tt(iwx1(i))**2+
     1      1.27e7*tt(iwx1(i))/(qsatv(iwx1(i))*pr(iwx1(i)) )
c
         c = 2.598*(tp7/w6(i))**0.666667/tp6
c
         prevp(iwx1(i)) = c * w1(i)*
     1      ( nr(iwx1(i))*nr(iwx1(i))*qr(iwx1(i)) )**0.3333333
c
         prevp(iwx1(i))=min(w5(i),prevp(iwx1(i)))
         if(delt*prevp(iwx1(i)).lt.pcut)prevp(iwx1(i))=0.
        enddo
c<--end of addition
       endif

      endif
c
c***********************************************************************
c
      return
      end
