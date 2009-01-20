       subroutine nrmtqi(qvtemp,qitemp,qgtemp,qstemp,qadj,tttemp
     1            ,pgdep,pdepi,psdep,pint,cond,pcut,twodt
     2            ,lsocp,tzero,lgrpl,ary1d,maxpt1d,kmax,maxvr1d,n1d
     3            ,i1dflg,i1d,j1d,j,icomp,th2t,conci,m,len)
c
c SCCS IDENTIFICATION:  %W% %G%
       implicit none
c
       integer m
       integer i1dflg
       integer icomp(m)
       integer ip
       integer j
       integer jp
       integer kk
       integer kmax
       integer len
       integer maxpt1d
       integer maxvr1d
       integer n1d
       integer np
c
       real ary1d   (maxpt1d,kmax,maxvr1d)
       real i1d     (n1d)
       real j1d     (n1d)
       real lsocp
       real cond    (len)
       real conci   (len)
       real pcut
       real pgdep   (len)
       real pdepi   (len)
       real psdep   (len)
       real qvtemp  (len)
       real pint    (len)
       real qadj    (len)
       real qitemp  (len)
       real qgtemp  (len)
       real qstemp  (len)
       real th2t    (len)
       real tttemp  (len)
       real twodt
       real tzero 
c
       logical lgrpl
c
c***********************************************************************
c      local variables
c***********************************************************************
c
       integer i 
       integer iconsum
       integer k
c
       real dqadj
       real qsum
       real tp3
       real wgt
c
       iconsum=0
c
       do i=1,len
c
         wgt=1.
         dqadj=0.
         if(tttemp(i).lt.tzero)then
c
          qsum= (psdep(i)+pdepi(i)+pint(i)
     1       +pgdep(i) )*twodt
c
          if(abs(qsum).gt.pcut)then
c
           if(qsum.lt.0)then
c
c***********************************************************************
c treat subsaturated conditions and limit ice-to-vapor transition to dqadj
c***********************************************************************
c
            qsum=-qsum
            dqadj=qadj(i)-qvtemp(i)
            if(qsum.gt.dqadj)then
             wgt=dqadj/qsum
            endif
           elseif(qsum.gt.0.)then
c
c***********************************************************************
c treat saturated conditions
c***********************************************************************
c
            dqadj=qvtemp(i)-qadj(i)
c
c***********************************************************************
c note: if it is desired to consume all vapor available during the 
c time step then set iconsum=1
c***********************************************************************
c
            if(iconsum.eq.1 .or. qsum.gt.dqadj)then
             wgt=dqadj/qsum
            else
             wgt=1.0
            endif

           endif
c
c***********************************************************************
c set weight function to normalize rates. The vapor demand cannot exceed
c dqadj over a time step. We may also want all possible vapor consumed so 
c increase each rate by dqadj/qsum if qsum is less than dqadj and total
c consumption is desired (iconsum=1)
c***********************************************************************
c
c***********************************************************************
c   update microphysical tendencies (exclude for the moment those
c   associated with the conversion terms: i.e., praut,pracw,pconv,psaci
c   as they do not impact the temperature or vapor tendencies)
c   Limit updates to the interior regions of the grid
c***********************************************************************
c
c
c be sure pgdep defined for t< 0 only
c
           qvtemp(i)=qvtemp(i)-
     1             wgt*(psdep(i)
     2              +pdepi(i)+pint(i)
     3              +pgdep(i) )*twodt
           qitemp(i)=qitemp(i)
     1             +wgt*(pdepi(i)
     2             +pint(i))*twodt

           qstemp(i)=qstemp(i)
     1             +(wgt*psdep(i))*twodt
c
           qgtemp(i)=qgtemp(i)
     1             +(wgt*pgdep(i))*twodt
c
           tttemp(i)=tttemp(i)+
     1              lsocp*twodt*wgt*(psdep(i)+pdepi(i)
     2             +pint(i)+pgdep(i) ) 

           cond(i)  =cond(i)+
     1             (lsocp*twodt*wgt*(psdep(i)+pdepi(i)
     2             +pint(i)+pgdep(i) ) )
          endif
c
         endif
       enddo
c
c***********************************************************************
c limit values
c***********************************************************************
c
       do i=1,len
         if(tttemp(i).lt.tzero)then
! Michael Falk, 24 August 2006 and 28 August 2006
! put it back the original way, 30 August 2006
          if(qvtemp(i).lt.pcut)qvtemp(i)=pcut
!          if(qvtemp(i).lt.pcut) then
!           qvtemp(i)=pcut
!           print *,'limiting at ',i,' and qvtemp is ',qvtemp(i),
!     .           ' and pcut ',pcut
!          end if
! eMFc
          if(qgtemp(i).lt.pcut)qgtemp(i)=0.0
          if(qitemp(i).lt.pcut)qitemp(i)=0.0
          if(qstemp(i).lt.pcut)qstemp(i)=0.0
         endif
       enddo
c
c***********************************************************************
c
       return
       end

