        subroutine adjmlt(delta1,delta2,delta3,qctemp,qgtemp,qitemp
     1             ,qrtemp,qstemp,pchomo,pgacr,pgacrm,pgacw,pgacwm
     2             ,pgmlt,piacr,pracs,prhomo,psacw,psacr,psmlti
     3             ,psmlt,piacw,tttemp,cond,lfocp
     4             ,twodt,pcut,tzero,lgrpl,ary1d,maxpt1d,kmax
     5             ,maxvr1d,n1d,i1dflg,i1d,j1d,j,icomp,m,len)
c
       implicit none
c
       integer m
c
       integer i1dflg
       integer icomp(m)
       integer ip
       integer j
       integer jp
       integer kmax
       integer len
       integer maxpt1d
       integer maxvr1d
       integer n1d
       integer np
c
       real ary1d(maxpt1d,kmax,maxvr1d)
       real cond    (len)     
       real delta1  (len)     
       real delta2  (len)     
       real delta3  (len)     
       real i1d     (n1d)
       real j1d     (n1d)
       real lfocp
       real qctemp  (len)     
       real qgtemp  (len)     
       real qitemp  (len)     
       real qrtemp  (len)     
       real qstemp  (len)     
       real pchomo  (len)     
       real pcut
       real pgacr   (len)     
       real pgacrm  (len)     
       real pgacw   (len)     
       real pgacwm  (len)     
       real pgmlt   (len)     
       real piacr   (len)     
       real piacw   (len)     
       real pracs   (len)     
       real prhomo  (len)     
       real psacr   (len)     
       real psacw   (len)     
       real psmlt   (len)     
       real psmlti  (len)     
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
       integer k
c
       real abs1
       real abs2
       real abs3
       real abs4
       real abs5
       real abs6
       real abs7
       real abs8
       real abs9
       real abs10
       real abs11
       real abs12
       real abs13
       real abs14
       real abs15
       real abs16
       real abs17
       real qsum
       real tst
       real tp3
       real wgt
c
c***********************************************************************
c          adjust temperature and moisture fields over interior of
c          the grid for gain/loss of ice/liquid due to melting
c          and freezing
c***********************************************************************
c
c***********************************************************************
c       set up delta functions (eqa's a30-a32, RH 1984)
c***********************************************************************
c

       do i=1,len
         if(qstemp(i).gt..0001.and.qctemp(i).gt..0005
     1       .and.lgrpl)then
             delta3(i)=0.0
         else
             delta3(i)=1.0
         endif
       enddo
c
       do i=1,len
         if(qstemp(i).gt..0001.and.qrtemp(i).gt..0001
     1      .and.lgrpl)then
             delta2(i)=0.0
         else
             delta2(i)=1.0
         endif
       enddo
c
       do i=1,len
         if(qrtemp(i).gt..0001.and.lgrpl)then
             delta1(i)=0.0
         else
             delta1(i)=1.0
         endif
       enddo
c
c***********************************************************************
c        adjust mixing ratio terms affected by melting/fusion processes
c***********************************************************************
c
       do i=1,len
c
         if(tttemp(i).lt.tzero)then
c
c***********************************************************************
c      determine weights to prevent terms for over collecting a given
c      species
c***********************************************************************
c
          wgt=1.0
          qsum= abs( (-piacw(i)-psacw(i)-pchomo(i)-pgacw(i))
     1        *twodt)
c
          if(qsum.lt.pcut)then
           wgt=0.
          elseif(qsum.gt.qctemp(i))then
           wgt= qctemp(i)/qsum
          endif
c
          abs1=wgt
          abs2=wgt
          abs3=wgt
          abs16=wgt
c
c***********************************************************************
c      determine weights for rain to prevent terms from over 
c      collecting a given species
c***********************************************************************
c
          wgt=1.0
          qsum= abs( (-prhomo(i)-pgacr(i)-piacr(i)
     1            -psacr(i) )*twodt)
c
          if(qsum.lt.pcut)then
           wgt=0.
          elseif(qsum.gt.qrtemp(i))then
           wgt= qrtemp(i)/qsum
          endif
c
          abs4=wgt
          abs5=wgt
          abs6=wgt
          abs7=wgt
c
          qctemp(i)=qctemp(i)
     1                +(-abs1*psacw(i)-abs2*pchomo(i)
     2                  -abs3*pgacw(i) -abs16*piacw(i))*twodt
c
          qstemp(i)=qstemp(i)
     1                 +(abs4*prhomo(i)
     1                 -pracs(i)*(1.-delta2(i))
     1                 +abs6*piacr(i)*delta1(i)
     1                 +abs7*psacr(i)*delta2(i)
     1                 +abs1*psacw(i)*delta3(i) )*twodt
c
          qitemp(i)=qitemp(i)+(abs2*pchomo(i)+abs16*piacw(i))
     1               *twodt
c
          qrtemp(i)=qrtemp(i)+(-abs4*prhomo(i)
     1                -abs5*pgacr(i)-abs6*piacr(i)
     1                -abs7*psacr(i) )*twodt
c
          if(lgrpl)qgtemp(i)=qgtemp(i)+(abs3*pgacw(i)
     1                +abs5*pgacr(i)
     1                +abs6*piacr(i)*(1-delta1(i))
     2                +pracs(i)*(1.-delta2(i))
     4                +abs7*psacr(i)*(1.-delta2(i))
     5                +abs1*psacw(i)*(1-delta3(i)) )*twodt
c
          tttemp(i)=tttemp(i)
     1                +lfocp*twodt*(abs1*psacw(i)+abs2*pchomo(i)
     2                +abs4*prhomo(i)
     3                +abs6*piacr(i)+abs5*pgacr(i)
     4                +abs3*pgacw(i)
     5                +abs7*psacr(i) )
c
          cond(i)=cond(i)
     1                +lfocp*twodt*(abs1*psacw(i)+abs2*pchomo(i)
     2                +abs4*prhomo(i)
     3                +abs6*piacr(i)+abs5*pgacr(i)
     4                +abs3*pgacw(i)
     5                +abs7*psacr(i) )
         else
c
c***********************************************************************
c      determine weights for graupel to prevent terms 
c      for over collecting
c***********************************************************************
c
          wgt=1.0
          qsum= abs( (pgmlt(i)+pgacrm(i)+pgacwm(i))*twodt)
c
          if(qsum.lt.pcut)then
            wgt=0.
          elseif(qsum.gt.qgtemp(i))then
            wgt= qgtemp(i)/qsum
          endif
c
          abs10=wgt
          abs11=wgt
          abs12=wgt
c
c***********************************************************************
c      determine weights for snow to prevent terms 
c      from producing negative snow mixing ratio 
c***********************************************************************
c
          wgt=1.0
          qsum= abs( (psmlt(i)-pracs(i))*twodt)
c
          if(qsum.lt.pcut)then
            wgt=0.
          elseif(qsum.gt.qstemp(i))then
            wgt= qstemp(i)/qsum
          endif
c
          abs13=wgt
          abs14=wgt
c
c***********************************************************************
c      determine weights for cloud water to prevent terms 
c      from producing negative cloud water mixing ratio 
c***********************************************************************
c
          wgt=1.0
          qsum= (pgacw(i)+psacw(i))*twodt

          if(qsum.lt.pcut)then
            wgt=0.
          elseif(qsum.gt.qctemp(i))then
            wgt= qctemp(i)/qsum
          endif
c
          abs15=wgt
          abs17=wgt
c
c***********************************************************************
c note: pgshr = pgacr+pgacw, the budget term for qrtemp is
c               at temperatures above t > 0 is
c               qrtemp=pghsr-pgacr so there is a net gain
c               in rain of pgacw(i) and a net loss in graupel
c               of pgacrm+pgacwm due to the collection of 
c               cloud droplets and rain
c***********************************************************************
c
          qctemp(i)=qctemp(i)+(psmlti(i)-abs17*pgacw(i)
     1               -abs15*psacw(i))*twodt
c
          qrtemp(i)=qrtemp(i)
     1                +(abs15*psacw(i)-abs13*psmlt(i)
     2                -abs11*pgacrm(i)-abs12*pgacwm(i)
     3                +abs14*pracs(i)
     4                -abs10*pgmlt(i)+abs17*pgacw(i) )*twodt
c
          qitemp(i)=qitemp(i)-psmlti(i)*twodt
c
          qgtemp(i)=qgtemp(i)+(abs10*pgmlt(i)
     1                +abs11*pgacrm(i)+abs12*pgacwm(i) )*twodt
c
          qstemp(i)=qstemp(i)+(abs13*psmlt(i)
     1                -abs14*pracs(i) )*twodt
c
          tttemp(i)=tttemp(i)
     1                +lfocp*twodt*(abs13*psmlt(i)-psmlti(i)
     2                +abs10*pgmlt(i)+abs12*pgacwm(i)
     3                +abs11*pgacrm(i)
     4                -abs14*pracs(i) )
c
          cond(i)=cond(i)
     1                +lfocp*twodt*(abs13*psmlt(i)-psmlti(i)
     2                +abs10*pgmlt(i)+abs12*pgacwm(i)
     3                +abs11*pgacrm(i)
     4                -abs14*pracs(i) )
         endif
c
       enddo
c
c***********************************************************************
c assure positive definiteness
c***********************************************************************
c
       do i=1,len
          if(qitemp(i).le.pcut)qitemp(i)=0.
          if(qrtemp(i).le.pcut)qrtemp(i)=0.
          if(qctemp(i).le.pcut)qctemp(i)=0.
          if(qstemp(i).le.pcut)qstemp(i)=0.
          if(qgtemp(i).le.pcut)qgtemp(i)=0.
       enddo

       return
       end
