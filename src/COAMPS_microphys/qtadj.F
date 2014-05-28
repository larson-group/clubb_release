      subroutine qtadj(tt,qv,pr,qsatv,qsati,qadj,tadj,heatlv
     1  ,heatls,cp,rvap,tzero,lice,iflag,pcut,len)
c
      implicit none
c
c***********************************************************************
c           parameters:
c***********************************************************************
c
      integer i
      integer iflag
      integer implflg
      integer len
      integer n
c
      real cp
      real heatls
      real heatlv
      real pcut
      real qadj   (len)
      real qv     (len)
      real qsat
      real qsati  (len)
      real qsatv  (len)
      real pr     (len)
      real rvap
      real tadj   (len)
      real tt     (len)
      real tzero
c
      logical lice
c
c***********************************************************************
c          local variables and dynamic storage:
c***********************************************************************
c
      real a
      real dq 
      real dt 
      real tp1
      real tp2
      real tp3
      real tp4
c
c************************************************************************
c          local constants
c************************************************************************
c
      tp1=heatlv*heatlv/(cp*rvap)
      tp2=heatls*heatls/(cp*rvap)
      tp3=heatls/cp
      tp4=heatlv/cp
      a=7.5*alog(10.)
c
      do i=1,len
       tadj(i)=-999.
       qadj(i)=-999.
      enddo
c
c************************************************************************
c use implicit technique:
c************************************************************************
c
      do i=1,len
c
        if (iflag.eq.2) then
         if (tt(i).lt.tzero) then
          qadj(i)=qsati(i)*(1.+tp2*qv(i)/tt(i)**2.0)/
     1            (1.+tp2*qsati(i)/tt(i)**2.0)
          dq=qv(i)-qadj(i)
          dt=tp3*dq
          qadj(i)=max(pcut,qadj(i))
          tadj(i) =tt(i)+dt
c
         elseif(tt(i).ge.tzero)then
c
          qadj(i)=qsatv(i)*(1.+tp1*qv(i)/tt(i)**2.0)/
     1            (1.+tp1*qsatv(i)/tt(i)**2.0)
          dq=qv(i)-qadj(i)
          dt=tp4*dq
          qadj(i)=max(pcut,qadj(i))
          tadj(i) =tt(i)+dt
c
         endif
c
        elseif(iflag.eq.1)then
c
          qadj(i)=qsatv(i)*(1.+tp1*qv(i)/tt(i)**2.0)/
     1            (1.+tp1*qsatv(i)/tt(i)**2.0)
          dq=qv(i)-qadj(i)
c
          dt=tp4*dq
          qadj(i)=max(pcut,qadj(i))
          tadj(i) =tt(i)+dt
c
        endif
c
       enddo

      return
      end
