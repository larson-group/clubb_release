       subroutine nrmcol(tttemp,qctemp,qitemp,qgtemp,qrtemp,qstemp
     1            ,pconv,pgaci,pgacs,praci,pracw,praut,psaci,pwacs
     2            ,rho,nc,nr,ncn,pcut,twodt,lice,lgrpl,tt,rvc
     3            ,ary1d,maxpt1d,kmax,maxvr1d,n1d,i1dflg,i1d,j1d,j
     4            ,icomp,ldrizzle,ni,m,len)
c
c SCCS IDENTIFICATION:  %W% %G%
       implicit none
c
       integer m
       integer i1dflg
       integer icomp(m)
       integer j
       integer kmax
       integer len
       integer maxpt1d
       integer maxvr1d
       integer n1d
c
       real ary1d   (maxpt1d,kmax,maxvr1d)
       real i1d     (n1d)
       real j1d     (n1d)
       real nc      (len)
       real ncn     (len)
       real ncplus
       real ni      (len)
       real nr      (len)
       real pcut
       real pconv   (len)
       real pgaci   (len)
       real pgacs   (len)
       real praci   (len)
       real pracw   (len)
       real praut   (len)
       real psaci   (len)
       real pwacs   (len)
       real qcplus
       real qrplus
       real qctemp  (len)
       real qgtemp  (len)
       real qitemp  (len)
       real qrtemp  (len)
       real qstemp  (len)
       real rho     (len)
       real tt      (len)
       real tttemp  (len)
       real rvc     (len)
       real twodt
c
       logical ldrizzle
       logical lice
       logical lgrpl
c
c***********************************************************************
c      local variables
c***********************************************************************
c
       integer i 
       integer k
c
       real abs1
       real abs2
       real abs3
       real abs4
       real abs5
       real abs6
       real delt1
       real delt3
       real qsum
c      real rvc
       real rvcof
       real tem
       real tp2
       real wgt
c
c***********************************************************************
c tem=1./(den_water*pi*4/3)
c***********************************************************************
c
      tem=2.3873241e-4
c
c***********************************************************************
c     update fields due to collection terms
c***********************************************************************
c
      do i=1,len
c
c***********************************************************************
c      determine weights to prevent terms for over collecting a given
c      species
c***********************************************************************
c
        wgt=1.0
        qsum= abs( (praut(i)+pracw(i))*twodt)
c
        if(qsum.lt.pcut)then
         wgt=0.
        elseif(qsum.gt.qctemp(i))then
         wgt= qctemp(i)/qsum
        endif
c
        abs1=wgt
        abs2=wgt
c
        tp2=qrtemp(i)
        qctemp(i)=qctemp(i)-
     1             (abs1*praut(i)+abs2*pracw(i))*twodt
        qrtemp(i)=qrtemp(i)
     1             +(abs1*praut(i)+abs2*pracw(i))*twodt
c
c-->Addition from Jerry Schmidt
c***********************************************************************
c update nc from autoconversion (Kogan and Khairoutdinov 2000) eq(32)
c why 28 microns (as in their original code)? The paper shows 25 microns
c***********************************************************************
c
        if(ldrizzle)then
         rvcof=rho(i)*tem

c-->Chris Golaz comments: 
c   Why not rvc for nc as in Eq (35) ?
c   Is it 0.0025, 0.0028 ?
c<--end of comments

         nc(i)=max(0., nc(i)-rvcof*abs1*praut(i)*twodt
     1                  /0.0025**3)
         nr(i)=nr(i)+rvcof*abs1*twodt*praut(i)/0.0028**3

c
c***********************************************************************
c update nc from accretion (Kogan and Khairoutdinov 2000) eq(35)
c***********************************************************************
c
         qcplus = max(0.0, qctemp(i) )
         qrplus = max(0.0, qrtemp(i) )
         ncplus = max(0.01,nc(i) )
c ---> dschanen made rvc an output for stats purposes
c        rvc=(rvcof*qcplus/ncplus)**0.333333
         rvc(i)=(rvcof*qcplus/ncplus)**0.333333
c <--- end dschanen modification
c***********************************************************************
c  apparently there was a bug here with Mechem's original code autoconv.F
c  as rvc can be zero if qc = 0.
c***********************************************************************
c
         if(rvc(i).ne.0.)then
          nc(i) = max(0.0,nc(i) - rvcof*abs2*twodt*pracw(i)/rvc(i)**3)
         else
          nc(i)=0.0
         endif

c-->Chris Golaz modification: these lines cause problem during initialization
c   if the initial sounding is saturated. They will cause a spurious
c   increase in the cloud droplet number concentration.
!         if (rvc(i).gt.0.0025) then
!           nc(i) = rvcof * qcplus / (0.0025)**3
!         endif
c<--end of modification
        endif
c
c<--end of addition
      enddo
c
      if(lice)then
c
        do i=1,len
c
         if(qrtemp(i).gt..0001.and.lgrpl)then
           delt1=0.0
         else
           delt1=1.0
         endif
c
         if(qctemp(i).gt..0005.and.qstemp(i).gt..0001.and.lgrpl)then
           delt3=0.0
         else
           delt3=1.0
         endif
c
c***********************************************************************
c      determine weights to prevent terms for over collecting a given
c      species
c***********************************************************************
c
         wgt=1.0
         qsum= abs( (pconv(i)+psaci(i)+pgaci(i)
     1            + praci(i) )*twodt)
c
         if(qsum.lt.pcut)then
          wgt=0.
         elseif(qsum.gt.qitemp(i))then
          wgt= qitemp(i)/qsum
         endif
c
         abs1=wgt
         abs2=wgt
         abs3=wgt
         abs4=wgt
c
         wgt=1.0
         qsum=  (pgacs(i)+pwacs(i))*twodt
c
         if(qsum.lt.pcut)then
          wgt=0.
         elseif(qsum.gt.qstemp(i))then
          wgt= qstemp(i)/qsum
         endif
c
         abs5=wgt
         abs6=wgt
c
         qitemp(i)=qitemp(i)
     1             +(-abs1*pconv(i)-abs2*psaci(i)
     2               -abs3*pgaci(i)-abs4*praci(i) )*twodt

         qstemp(i)=qstemp(i)
     1             +(abs1*pconv(i)+abs2*psaci(i)
     2             +abs4*praci(i)*delt1-abs5*pgacs(i)
     3             -abs6*pwacs(i)*(1.0-delt3) )*twodt

         if(lgrpl)qgtemp(i)=qgtemp(i)
     1             +(abs3*pgaci(i)+abs5*pgacs(i)
     2             +abs4*praci(i)*(1.-delt1) 
     3             +abs6*pwacs(i)*(1.-delt3) )*twodt

        enddo
      endif
c
c***********************************************************************
c limit values
c***********************************************************************
c
      do i=1,len
c
       if(qgtemp(i).lt.pcut)qgtemp(i)=0.0
       if(qitemp(i).lt.pcut)qitemp(i)=0.0
       if(qstemp(i).lt.pcut)qstemp(i)=0.0
c
c
c***********************************************************************
c special checks on concentrations
c***********************************************************************
c
       if(ldrizzle)then
        if(nc(i).lt.0.)then
         ncn(i)=ncn(i)-nc(i)
         nc(i)=0.
        endif
        if(nr(i).lt.0.)then
         ncn(i)=ncn(i)-nr(i)

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
