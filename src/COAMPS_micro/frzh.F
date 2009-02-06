      subroutine frzh(pchomo,tt,qc,pcut,delt,len)
      integer len
c
      real delt
      real pcut
      real pchomo (len)
      real qc     (len)
      real tt     (len)
c
c***********************************************************************
c          local variables and dynamic storage:
c***********************************************************************
c
      real thomo
      integer i
      real tp1
c
c************************************************************************
c          local constants
c************************************************************************
c
      thomo=233.15
      tp1=1.0/delt
c
      do i=1,len
        if (tt(i).le.thomo.and.qc(i).ge.pcut) pchomo(i)=qc(i)*tp1
      enddo
c
      return
      end
