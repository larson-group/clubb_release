       subroutine nrmtqw(qvtemp,qgtemp,qrtemp,qstemp,qadj,qctemp,tttemp
     1            ,pcond,prevp,pmltge,pmltse,cond,nc,nr,ncn,pcut,twodt
     2            ,lvocp,lgrpl,ary1d,maxpt1d,kmax,maxvr1d,n1d,i1dflg
c-->Chris Golaz modification: added icase
!     3            ,i1d,j1d,j,ldrizzle,th2t,heatlv,cp,rvap,qsatv,len)
     3            ,i1d,j1d,j,ldrizzle,th2t,heatlv,cp,rvap,qsatv,len
     3            ,icase)
c<--end of modification
c
c SCCS IDENTIFICATION:  %W% %G%
       implicit none
c
c-->Chris Golaz modification: added icase
       integer icase
c<--end of modification
       integer i1dflg
       integer j
       integer kk
       integer kmax
       integer len
       integer m
       integer maxpt1d
       integer maxvr1d
       integer n1d
c
       real ary1d(maxpt1d,kmax,maxvr1d)
       real cond    (len)
       real cp
       real heatlv
       real i1d (n1d)
       real j1d (n1d)
       real lvocp
       real nc      (len)
       real ncn     (len)
       real nr      (len)
       real pcond   (len)
       real pcut
       real pmltge  (len)
       real pmltse  (len)
       real prevp   (len)
       real qadj    (len)
       real qctemp  (len)
       real qgtemp  (len)
       real qrtemp  (len)
       real qsatv   (len)
       real qstemp  (len)
       real qvtemp  (len)
       real rvap
       real th2t    (len)
       real tttemp  (len)
       real twodt
c
       logical lgrpl
       logical ldrizzle
c
c***********************************************************************
c      local variables
c***********************************************************************
c
       integer i 
       integer k
c
       real dqadj
       real dqv
       real qcplus
       real qrplus
       real qsum
       real qthresh
       real temp1
       real temp2
       real wgt
       real xx
c
       qthresh=0.1e-3
c
c***********************************************************************
c compute weights used to normalize microphysical rates
c The vapor demand cannot exceed dqadj over a time step.
c***********************************************************************
c
       do i=1,len
        wgt=1.0 
        dqadj=0.  
        qsum= (pcond(i)+prevp(i)
     1      +pmltge(i)+pmltse(i))*twodt
c
        if(abs(qsum).gt.pcut)then
c
         if(qsum.lt.0)then
c
c***********************************************************************
c treat subsaturated conditions
c***********************************************************************
c
          qsum=-qsum
          dqadj=qadj(i)-qvtemp(i)
c
c***********************************************************************
c limit liquid-to-vapor transition to dqadj
c***********************************************************************
c
          if(qsum.gt.dqadj)then
           wgt=dqadj/qsum
          endif
c
         elseif(qsum.gt.0.)then
c
c***********************************************************************
c treat saturated conditions and insure we consume dqadj of vapor
c***********************************************************************
c
          dqadj=qvtemp(i)-qadj(i)
          wgt=dqadj/qsum
         endif
c
c***********************************************************************
c Khairoutdinov and Kogan (2000) scheme needs info on saturation 
c***********************************************************************
c 
         dqv=qadj(i)-qvtemp(i)
         qcplus=max(0.,qctemp(i))
         qrplus=max(0.,qrtemp(i))
c
c***********************************************************************
c update microphysical tendencies and temperature for warm microphysics
c***********************************************************************
c
         qvtemp(i)=qvtemp(i)-
     1             (wgt*pcond(i)
     2             +wgt*prevp(i)
     3             +wgt*pmltse(i)
     4             +wgt*pmltge(i))*twodt
         qctemp(i)=qctemp(i)
     1             +wgt*pcond(i)*twodt
c
         qrtemp(i)=qrtemp(i)
     1             +wgt*prevp(i)*twodt
c
         qstemp(i)=qstemp(i)
     1             +wgt*pmltse(i)*twodt
c
         if(lgrpl)qgtemp(i)=qgtemp(i)
     1             +wgt*pmltge(i)*twodt
c
         tttemp(i)=tttemp(i)
     1             +lvocp*twodt*wgt*(pcond(i)
     2             +prevp(i)
     3             +pmltse(i)
     4             +pmltge(i) )
c
         cond(i)  = cond(i)
     1             +lvocp*twodt*wgt*(pcond(i)
     2             +prevp(i)
     3             +pmltse(i)
     4             +pmltge(i) )
c
         if(ldrizzle)then
c
c***********************************************************************
c     Khairoutdinov and Kogan (2000) drizzle parameterization
c***********************************************************************
c
c     Adjust ncd, ncn according to the empirical formula of
c     Martin et al. (1994) and O'Dowd et al. (1996).
c
c     Martin et al. (1994):
c
c        Nc=-1.15E-3*Ncn**2+0.963*Ncn+5.30      if 36<=Ncn<=280 cm-3
c        Nc=-2.10E-4*Ncn**2+0.568*Ncn-27.9      if 375<=Ncn<=1500 cm-3
c
c     O'Dowd et al. (1996):
c
c        Nc=197.*(1.-exp(-6.13E-3*Ncn))
c

c
c***********************************************************************
c     CCN activation if saturated dqv<0
c***********************************************************************
c
          if(dqv.lt.0.)then

c-->Chris Golaz modification: possible bug fix
!           temp1=max(0.,ncn(i))
           temp1=max(0.,ncn(i)+nc(i))
c<--end of modification
c
c-->Chris Golaz modification: for LES icase=68 and icase=75 use 
c   bulk activation
!           if(temp1.ge.36..and.temp1.le.280.)then
!            temp2=-1.15E-3*temp1*temp1+0.963*temp1+5.30
!           else if(temp1.ge.375.)then
!            temp2=-2.10E-4*temp1*temp1+0.568*temp1-27.9
!           else
!            temp2=197.*(1.-exp(-6.13E-3*temp1))
!           end if

           if ( .not. (icase.eq.68 .or. icase.eq.75) ) then

            if(temp1.ge.36..and.temp1.le.280.)then
             temp2=-1.15E-3*temp1*temp1+0.963*temp1+5.30
            else if(temp1.ge.375.)then
             temp2=-2.10E-4*temp1*temp1+0.568*temp1-27.9
            else
             temp2=197.*(1.-exp(-6.13E-3*temp1))
            end if

           else

             temp2 = temp1

           endif
c<--end of modification
c
c***********************************************************************
c modify CCN based on cloud water threshold (qthresh)
c***********************************************************************
c
           temp2=max(temp2,nc(i))
           temp1=nc(i)
           xx=qcplus/qthresh
           if(xx.lt.1.)temp2=temp2*xx
           if(temp2.gt.temp1)then
            nc(i)=temp2
            ncn(i)=ncn(i)-(temp2-temp1)
           end if
c
          end if
c
c***********************************************************************
c     Reduce CCN if subsaturated
c***********************************************************************
c
          if(dqv.gt.0.and.qcplus.lt.0.05e-3.and.pcond(i).lt.0.)then
           temp1=nc(i)*(-wgt*pcond(i)*twodt)/(qcplus+1.e-20)
           temp2=min(nc(i),temp1)
           ncn(i)=ncn(i)+temp2
           nc(i)=nc(i)-temp2
          end if
c
          if(qctemp(i).le.0.)then
           ncn(i)=ncn(i)+nc(i)
           nc(i)=0.
          endif
         endif
c
c***********************************************************************
c     Reduce Rain concentration if subsaturated (eq. 23)
c***********************************************************************
c
         if(dqv.gt.0..and.qrplus.ge.pcut.and.prevp(i).lt.0.)then
           temp1=nr(i)*(-wgt*prevp(i))*twodt/qrplus
           temp1=min(nr(i),temp1)
           nr(i)=nr(i)-temp1
           ncn(i)=ncn(i)+temp1
         endif

        endif
       enddo
c
c***********************************************************************
c limit values
c***********************************************************************
c
       do i=1,len
! Michael Falk, 24 August 2006 and 28 August 2006
! put it back the original way, 30 August 2006
          if(qvtemp(i).lt.pcut)qvtemp(i)=pcut
!          if(qvtemp(i).lt.pcut) then
!           qvtemp(i)=pcut
!           print *,'limiting in nrmtqw at ',i
!          end if
! eMFc
        if(qgtemp(i).lt.pcut)qgtemp(i)=0.0
        if(qstemp(i).lt.pcut)qstemp(i)=0.0
        if(qrtemp(i).lt.pcut)qrtemp(i)=0.0
        if(qctemp(i).lt.pcut)qctemp(i)=0.0
        if(ldrizzle)then
         if(nc(i).lt.0.)then
c-->Chris Golaz modification: possible bug fix
!          ncn(i)=ncn(i)-nc(i)
          ncn(i)=ncn(i)+nc(i)
c<--end of modification
          nc(i)=0.
         endif

         if(nr(i).lt.0.)then
c-->Chris Golaz modification: possible bug fix
!          ncn(i)=ncn(i)-nr(i)
          ncn(i)=ncn(i)+nr(i)
c<--end of modification
          nr(i)=0.
         endif
        endif
c
        if(qctemp(i).lt.pcut)then
         qctemp(i)=0.0
         if(ldrizzle)then
          ncn(i)=ncn(i)+nc(i)
          nc(i)=0.
         endif
        endif
c
        if(qrtemp(i).lt.pcut)then
         qrtemp(i)=0.0
         if(ldrizzle)then
          ncn(i)=ncn(i)+nr(i)
          nr(i)=0.
         endif
        endif

       enddo
c
      return
      end

