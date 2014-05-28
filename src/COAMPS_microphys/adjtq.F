      subroutine adjtq
     1  (cond,p3,qc3,qi3,qr3,qg3,qs3,qv3,th3,w3,pr3d
     2  ,qsatv3d,temp3d,qsati3d,th2t3d,wtm,exbm,rbm
     3  ,nc3,nr3,ncn3,ni3,cp,deltf,heatlf,heatls,heatlv
     4  ,pcut,p00,r,rvap,sloper,slopes,slopeg,timea,lice
     5  ,nne,kk,i1d,j1d,ary1d,i1dflg,n1d,maxpt1d,maxvr1d
     6  ,kmax,nrdamp,nipts,nkpts,icomp,kcomp,j
     7  ,xland,aa0,aa1,aa2,aa3,abar,apr,aprpr,bsnow
     8  ,cbeta,cnzero,cimass,cpor,cw,difvap,erc,esi,eic
     9  ,eri,egc,esc,esr,egi,egr,egs,mw,pi,praut1,praut2
     1  ,rholiq,rhosno,rnzero,snzero,gnzero,therco,tice
     2  ,tvr1,tvr2,tvr3,tvr4,tzero,visair,gm3,gm4,gm5,gm6
     3  ,gm7,gm8,gm9,gmbp3,gmbov2,gmbov2g,bgrp,ex1,ex2
     4  ,ex2g,ex3,hlvoka,hkaolf,hlsoka,hlvorv,hlsorv
     5  ,rvochi,lfocp,lvocp,lsocp,ex7,ex7g,ex4,ex4g,ex5
c-->Chris Golaz modification: added icase
!     6  ,ldrizzle,lgrpl,icon,icond,len)
!     6  ,ldrizzle,lgrpl,icon,icond,len,icase)
c-->Adam Smith modification
c    08Aug2006: Added terms for microphysics source budgets
     6  ,ldrizzle,lgrpl,icon,icond,len,icase,snowv,snowslope
     7  ,pcond3,psmlti3,psacw3,pgacw3,piacw3,pchomo3,praut3,pracw3
     8  ,pdepi3,pint3,pgdep3,pconv3,psaci3,pgaci3,praci3,prevp3,psdep3
     9  ,pmltge3,pgmlt3,psmlt3,pgacrm3,pgacwm3,pracs3,pgshr3
! Michael Falk, 3 Oct 2006
!     1  ,pgacr3,psacr3,piacr3,prhomo3,pgacs3,pmltse3,pwacs3)
     1  ,pgacr3,psacr3,piacr3,prhomo3,pgacs3,pmltse3,pwacs3
     2  ,falli,fallg,fallr,falln,falls,rvc,rvr)
! eMFc

c<--end of ajsmith4's modification
c<--end of modification
c                                                                      c
#include "adjtq.prol"
c SCCS IDENTIFICATION:  %W% %G%
c***********************************************************************
c***********************************************************************
c***********************************************************************
c       subroutine adjtq:
c
c          Modified bulk microphysical parameterization developed for 
c          COAMPS3.0.  The details of this routine, together with the 
c          relevant papers upon which it is based, can be found in 
c          Chapter 4 of the COAMPS3.0 user's guide. J.M. Schmidt
c
c***********************************************************************
c***********************************************************************
c***********************************************************************
c
c Michael Falk started changing here, 25 Aug 2006
c      use domdec
#define iminf(nne) 1
#define imaxf(nne) 1
#define jminf(nne) 1
#define jmaxf(nne) 1
c eMFc
      implicit none
c
c***********************************************************************
c           parameters:
c***********************************************************************
c
c-->Chris Golaz modification: added icase
      integer   icase
c<--end of modification
      integer   nipts
      integer   nne
      integer   kk
c
      integer   i1dflg
      integer   icomp(nipts)
      integer   icon
      integer   icond
      integer   kcomp(iminf(nne):imaxf(nne),kk)
      integer   kmax
      integer   len
      integer   maxpt1d
      integer   maxvr1d
      integer   n1d
      integer   nrdamp
      integer   nkpts(iminf(nne):imaxf(nne))
c
      real    aa0
      real    aa1
      real    aa2
      real    aa3
      real    abar
      real    apr
      real    aprpr
      real    ary1d(maxpt1d,kmax-1,maxvr1d)
      real    bgrp
      real    bsnow
      real    cbeta
      real    cimass
      real    cnzero
      real    cpor
      real    cond   (iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk)
      real    cp
      real    cw
      real    deltf
      real    difvap
      real    exbm   (iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk)
      real    egc
      real    egi
      real    egr
      real    egs
      real    eic
      real    erc
      real    eri
      real    esc
      real    esi
      real    esr
      real    ex1
      real    ex2
      real    ex2g
      real    ex3
      real    ex4
      real    ex4g
      real    ex5
      real    ex7
      real    ex7g
      real    gm3
      real    gm4
      real    gm5
      real    gm6
      real    gm7
      real    gm8
      real    gm9
      real    gmbov2
      real    gmbov2g
      real    gmbp3
      real    gnzero
      real    heatlf
      real    heatls
      real    heatlv
      real    hlsorv
      real    hlvoka
      real    hkaolf
      real    hlsoka
      real    hlvorv
      real    i1d (n1d)
      real    j1d (n1d)
      real    lfocp
      real    lsocp
      real    lvocp
      real    mw
      real    nc3     (iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk)
      real    ncn3    (iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk)
      real    ni3     (iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk)
      real    nr3     (iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk)
      real    p00
      real    p3      (iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk)
      real    pcut
      real    pr3d    (iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk+1)
      real    praut1
      real    praut2
      real    qc3     (iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk)
      real    qg3     (iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk)
      real    qi3     (iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk)
      real    qr3     (iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk)
      real    qs3     (iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk)
      real    qsatv3d (iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk+1)
      real    qsati3d (iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk+1)
      real    qv3     (iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk)
      real    r
      real    rbm     (iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk)
      real    rholiq
      real    rhosno
      real    rnzero
      real    rvochi
      real    rvap
      real    sloper
      real    slopes
      real    slopeg
      real    temp3d  (iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk+1)
      real    therco
      real    tice
      real    th2t3d  (iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk+1)
      real    tvr1
      real    tvr2
      real    tvr3
      real    tvr4
      real    tzero
      real    visair
      real    th3     (iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk)
      real    timea
      real    w3      (iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk+1)
      real    wtm     (iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne))
      real    xland   (iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne))
c
c*************************************************************
c Addition by Adam Smith, 08 August 2006
c Adding output variables for snow fallspeed, slope parameter
c Also adding microphysical parameters for qc_mc term
c*************************************************************
c
      real snowv (iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk)
      real snowslope (iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk)

      real pcond3 (iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk)
      real psmlti3(iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk)
      real psacw3 (iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk)
      real pgacw3 (iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk)
      real piacw3 (iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk)
      real pchomo3(iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk)
      real praut3 (iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk)
      real pracw3 (iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk)
      real pdepi3 (iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk)
      real pint3  (iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk)
      real pgdep3 (iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk)
      real pconv3 (iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk)
      real psaci3 (iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk)
      real pgaci3 (iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk)
      real praci3 (iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk)
      real prevp3 (iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk)
      real psdep3 (iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk)
      real pmltge3(iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk)
      real pgmlt3 (iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk)
      real psmlt3 (iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk)
      real pgacrm3(iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk)
      real pgacwm3(iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk)
      real pracs3 (iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk)
      real pgshr3 (iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk)
      real pgacr3 (iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk)
      real psacr3 (iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk)
      real piacr3 (iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk)
      real prhomo3(iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk)
      real pgacs3 (iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk)
      real pmltse3(iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk)
      real pwacs3 (iminf(nne):imaxf(nne),jminf(nne):jmaxf(nne),kk)
c End of ajsmith4's addition

      logical lice
      logical lgrpl
      logical ldrizzle
c
c***********************************************************************
c          local variables and dynamic storage allocation
c***********************************************************************
c
      integer   i
      integer   ic
      integer   iflag
      integer   iindex
      integer   ilen
      integer   kindex
      integer   ip
      integer   j
      integer   jp
      integer   k
      integer   m
      integer   mkk
      integer   n
      integer   ni
      integer   niter
      integer   np
      integer   npts
      integer   iwx1    (len)
c
      real      apbp    (len)
      real      appbpp  (len)
      real      ccnw    (kk)
      real      ccnl    (kk)
      real      cm2m3
      real      delta1  (len)
      real      delta2  (len)
      real      delta3  (len)
      real      dmax
      real      fallr   (len)
      real      falln   (len)
      real      rvc     (len)
      real      rvr     (len)
      real      fallg   (len)
      real      falls   (len)
      real      falli   (len)
      real      nctemp  (len)
      real      nitemp  (len)
      real      nrtemp  (len)
      real      nntemp  (len)
      real      nuci    (len)
      real      pchomo  (len)
      real      pcond   (len)
      real      condrt  (len)
      real      pconv   (len)
      real      pdepi   (len)
      real      pgaci   (len)
      real      pgacr   (len)
      real      pgacrm  (len)
      real      pgacs   (len)
      real      pgacw   (len)
      real      pgacwm  (len)
      real      pgdep   (len)
      real      pgmlt   (len)
      real      pgshr   (len)
      real      piacw   (len)
      real      pi
      real      piacr   (len)
      real      pint    (len)
      real      pmltge  (len)
      real      pmltse  (len)
      real      pr      (len)
      real      praci   (len)
      real      pracs   (len)
      real      pracw   (len)
      real      praut   (len)
      real      prevp   (len)
      real      prhomo  (len)
      real      psaci   (len)
      real      psacr   (len)
      real      psacw   (len)
      real      psdep   (len)
      real      psmlt   (len)
      real      psmlti  (len)
      real      pterm   (len)
      real      pterm4  (len)
      real      pwacs   (len)
      real      qadj    (len)
      real      qadjw   (len)
      real      qctemp  (len)
      real      qdum    (len)
      real      qgtemp  (len)
      real      qitemp  (len)
      real      qrtemp  (len)
      real      qsati   (len)
      real      qsatv   (len)
      real      qsidum  (len)
      real      qvsdum  (len)
      real      qstemp  (len)
      real      qvtemp  (len)
      real      rbar2d  (len)
      real      sat
      real      sfgraup (len)
      real      sfrain  (len)
      real      sfraincm(len)
      real      sfsnow  (len)
      real      snzero
      real      tadj    (len)
      real      tadjw   (len)
      real      tdum    (len)
      real      temp
      real      terma   (len)
      real      termb   (len)
      real      termg   (len)
      real      th      (len)
      real      th2t    (len)
      real      tttemp  (len)
      real      twodt
      real      wgt     (len)
      real      x6      (len)
      real      x7      (len)
      real      x8      (len)
      real      x9      (len)
      real      xa      (len)
      real      xb      (len)
      real      xc      (len)
      real      tmp1    (len)
c
c***********************************************************************
c          local constants
c***********************************************************************
c

!Michael Falk, 29 Aug 2006
!      intent(in) j
! eMFc

      mkk=len
      cm2m3=1.0e6
      niter=1
      dmax=.005
c
c***********************************************************************
c set appropriate time step for hybrid or leapfrog schemes
c***********************************************************************
c
      twodt=deltf
c
      do i=1,len
          apbp  (i)=0.0
          condrt(i)=0.0
          delta1(i)=0.0
          delta2(i)=0.0
          delta3(i)=0.0
          falli (i)=0.0
          falls (i)=0.0
          fallr (i)=0.0
          falln (i)=0.0
          rvc   (i)=0.0
          rvr   (i)=0.0
          nctemp(i)=0.0
          nrtemp(i)=0.0
          nntemp(i)=0.0
          nitemp(i)=0.0
          nuci  (i)=0.0
          pcond (i)=0.0
          pconv (i)=0.0
          pdepi (i)=0.0
          pgaci (i)=0.0
          pgacr (i)=0.0
          pgacrm(i)=0.0
          pgacs (i)=0.0
          pgacw (i)=0.0
          pgacwm(i)=0.0
          pgdep (i)=0.0
          pgmlt (i)=0.0
          pgshr (i)=0.0
          piacr (i)=0.0
          piacw (i)=0.0
          pint  (i)=0.0
          pmltge(i)=0.0
          pmltse(i)=0.0
          pr    (i)=0.0
          prhomo(i)=0.0
          pchomo(i)=0.0
          praut (i)=0.0
          pracw (i)=0.0
          prevp (i)=0.0
          psacw (i)=0.0
          qctemp(i)=0.0
          qgtemp(i)=0.0
          qitemp(i)=0.0
          qrtemp(i)=0.0
          qstemp(i)=0.0
          praci (i)=0.0
          pracs (i)=0.0
          psaci (i)=0.0
          psacr (i)=0.0
          psdep (i)=0.0
          psmlt (i)=0.0
          psmlti(i)=0.0
          pwacs (i)=0.0
          qdum  (i)=0.0
          qvsdum(i)=0.0
          qsidum(i)=0.0
          tdum  (i)=0.0
          qadj  (i)=0.0
          qadjw (i)=0.0
          qsati (i)=0.0
          qsatv (i)=0.0
          rbar2d(i)=0.0
          sfgraup(i)=0.0
          sfrain(i)=0.0
          sfraincm(i)=0.0
          sfsnow(i)=0.0
          tadj  (i)=0.0
          tadjw (i)=0.0
          terma (i)=0.0
          termg (i)=0.0
          th2t  (i)=0.0
          th    (i)=0.0
          tttemp(i)=0.0
      enddo
c
c***********************************************************************
c     compute pressure (assume cp/r=3.5)
c     convert potential temperature to temperature
c     store temperature in temporary working array
c***********************************************************************
c
      ilen=0
      do i=1,nipts
       iindex=icomp(i)
       do k=1,nkpts(iindex)
        kindex=kcomp(iindex,k)
        ilen=ilen+1
        pr(ilen)  = pr3d(iindex,j,kindex)
        th2t(ilen)= th2t3d(iindex,j,kindex)
        th(ilen)  = th3(iindex,j,kindex)
        tttemp(ilen) = temp3d(iindex,j,kindex)
        qsatv(ilen) = qsatv3d(iindex,j,kindex)
        qsati(ilen) = qsati3d(iindex,j,kindex)
        rbar2d(ilen) = pr(ilen)/(r*tttemp(ilen))
       enddo
      enddo
c
c***********************************************************************
c initialize ccn if not using drizzle scheme
c***********************************************************************
c
      if(.not.ldrizzle)then
c
c***********************************************************************
c initialize ccn based on land/water table (#/cm**3)
c***********************************************************************
c
       do k=1,kk
        ccnl(k)=500.*cm2m3
        ccnw(k)=100.*cm2m3
       enddo
c
       ilen=0
       do i=1,nipts
        iindex=icomp(i)
        do k=1,nkpts(iindex)
         kindex=kcomp(iindex,k)
         ilen=ilen+1
         if(xland(iindex,j).gt.0.)then
           nctemp(ilen)  = ccnl(kindex)
         else
           nctemp(ilen)  = ccnw(kindex)
         end if
        enddo
       enddo
      end if
c
c***********************************************************************
c          store mixing ratios in temporary working arrays
c***********************************************************************
c
      ilen=0
      do i=1,nipts
       iindex=icomp(i)
       do k=1,nkpts(iindex)
        kindex=kcomp(iindex,k)
        ilen=ilen+1
        qctemp(ilen)=qc3(iindex,j,kindex)
        qrtemp(ilen)=qr3(iindex,j,kindex)
        qvtemp(ilen)=qv3(iindex,j,kindex)
       enddo
      enddo
c
      if(lice)then
       ilen=0
       do i=1,nipts
        iindex=icomp(i)
        do k=1,nkpts(iindex)
         kindex=kcomp(iindex,k)
         ilen=ilen+1
         qitemp(ilen)=qi3(iindex,j,kindex)
         qstemp(ilen)=qs3(iindex,j,kindex)
         qgtemp(ilen)=qg3(iindex,j,kindex)
         nitemp(ilen)=ni3(iindex,j,kindex)
        enddo
       enddo
      end if
c
      if(ldrizzle)then
       ilen=0
       do i=1,nipts
        iindex=icomp(i)
        do k=1,nkpts(iindex)
         kindex=kcomp(iindex,k)
         ilen=ilen+1
c-->Chris Golaz modification: 
c   nc3, nr3, ncn3 are now in mixing ratio units ( x 1.e-6 #/kg ) 
c   outside microphysics and must be converted to concentration (#/cm3)
c   here.

!         nctemp(ilen)=nc3(iindex,j,kindex)
!         nrtemp(ilen)=nr3(iindex,j,kindex)
!         nntemp(ilen)=ncn3(iindex,j,kindex)

         nctemp(ilen)=nc3(iindex,j,kindex)*
     .    rbm(iindex,j,kindex)
         nrtemp(ilen)=nr3(iindex,j,kindex)*
     .    rbm(iindex,j,kindex)
         nntemp(ilen)=ncn3(iindex,j,kindex)*
     .    rbm(iindex,j,kindex)

c<--end of modification
        enddo
       enddo
      end if
c
c***********************************************************************
c          compute slope factors for rain (sfrain) and for
c          snow (sfsnow) (eqs. 3a and 3b, rh-1983)
c          and graupel eqa 2 RH-1984
c***********************************************************************
c
      call slope(sfrain,sfraincm,sfsnow,sfgraup,qrtemp,qstemp,qgtemp
     1          ,rbar2d,sloper,slopes,slopeg,pcut,lice,lgrpl,mkk)
c
c***********************************************************************
c***********************************************************************
c          compute arrays to be shared in various terms to follow
c          (p00/p)**2
c***********************************************************************
c***********************************************************************
c
      do i=1,mkk
        pterm(i)=exp(ex1*log(p00/pr(i)))
        pterm4(i)=pterm(i)*pterm(i)
      enddo
c
c***********************************************************************
c          1.0/(a'+b'), 1.0/(a''+b''), terma and termb
c***********************************************************************
c
c      jms 04/14/00 bug fix in appbpp stems from error in 
c      eqa a16 of RH 1983
c       appbpp(i)=1.0/(xc(i,1)*hlvoka+appbpp(i,1))
c       replace hlvoka with hlsoka
c
c***********************************************************************
c
      call esatv(tttemp,apbp  ,mkk)
      call esati(tttemp,appbpp,mkk)
c
      do i=1,mkk
        xa(i)=tttemp(i)*rvochi
        apbp  (i)=xa(i)/(apbp  (i))
        appbpp(i)=xa(i)/(appbpp(i))
      enddo
c
      do i=1,mkk
        xa(i)=1.0/tttemp(i)
        xb(i)=(hlvorv*xa(i)-1.0)*xa(i)
        xc(i)=(hlsorv*xa(i)-1.0)*xa(i)
        apbp  (i)=1.0/(xb(i)*hlvoka+apbp  (i))
        appbpp(i)=1.0/(xc(i)*hlsoka+appbpp(i))
      enddo
c
      if(lice)then
c
c***********************************************************************
c compute terms needed in eqa's a18, a22,a25,a26,a27 of RH 1983
c***********************************************************************
c
       do i=1,mkk
        if (sfsnow(i).le.0.0) then
          terma(i)=1.0
          xa(i)=1.0
          termb(i)=1.0
        else
c
c****************************************************************
c Addition by Adam Smith, 27 March 2008
c "terma" is used to calculate PSDEP and other microphysical
c results.  For the Nov.11 and CLEX-9 cases, we use a simplified
c formula for calculating "terma".
c****************************************************************
c
!         terma(i)=0.65/(sfsnow(i)*sfsnow(i))
          if (icase.eq.1000 .or. icase.eq.1002 .or. icase.eq.1003) then
            terma(i)=2.70/(sfsnow(i)*sfsnow(i))
          else
            terma(i)=0.65/(sfsnow(i)*sfsnow(i))
          end if
c
c****************************
c End of ajsmith4's addition
c****************************
c
          temp=log(sfsnow(i))
          xa(i)=temp*ex2
          xa(i)=1.0/exp(xa(i))
          termb(i)=temp*ex3
          termb(i)=1.0/exp(termb(i))
        end if
       enddo
c
       do i=1,mkk
        xa(i)=xa(i)*pterm(i)*ex7
        xb(i)=rbar2d(i)*ex4
       enddo
c
       do i=1,mkk
        xa(i)=xa(i)*sqrt(xb(i))
c
c****************************************************************
c Addition by Adam Smith, 27 March 2008
c "terma" is used to calculate PSDEP and other microphysical 
c results.  For the Nov.11 and CLEX-9 cases, we use a simplified
c formula for calculating "terma".  For other cases, we need to
c modify "terma" using the standard formula.
c****************************************************************
c
        if (icase.ne.1000 .and. icase.ne.1002 .and. icase.ne.1003) then
          terma(i)=terma(i)+xa(i)
        end if
c
c****************************
c End of ajsmith4's addition
c****************************
c
        termb(i)=termb(i)*pterm4(i)*ex5
       enddo
c
c***********************************************************************
c compute terms needed in eqa's a17,a18,a19 of RH 1984
c***********************************************************************
c
       if(lgrpl)then
        do i=1,mkk
         if (sfgraup(i).le.0.0) then
          termg(i)=1.0
          xa(i)=1.0
         else
          termg(i)=0.78/(sfgraup(i)*sfgraup(i))
          temp=log(sfgraup(i))
          xa(i)=temp*ex2g
          xa(i)=1.0/exp(xa(i))
         end if
        enddo
c
        do i=1,mkk
         xa(i)=xa(i)*pterm(i)*ex7g
         xb(i)=rbar2d(i)*ex4g
        enddo
c
        do i=1,mkk
         xa(i)=xa(i)*sqrt(xb(i))
         termg(i)=2.0*pi*gnzero*(termg(i)+xa(i))
        enddo
       end if
c
c***********************************************************************
c         We will compute the melting and fusion terms first. These terms
c         affect the temperature but not the vapor field. Adjust
c         the temperature and saturation values and then compute 
c         the remaing six terms impacting the vapor field
c         Freeze any liquid below the homogeneous freezing temperature
c         (thomo=233 K).Frozen cloud water will be converted to the 
c         pristine ice category. Frozen rain water will be converted 
c         to the snow category
c***********************************************************************
c
       if(lgrpl)then
         call tgqg(fallg,qgtemp,sfgraup,pterm4,pcut,abar,bgrp,mkk)
       end if
c
       call tgqr(fallr,falln,qrtemp,sfraincm,nrtemp,rbar2d,pterm4,pcut,
     1           rvr,tvr1,tvr2,tvr3,tvr4,ldrizzle,mkk)
c
       call tgqs(falls,qstemp,sfsnow,pterm4,pcut,aprpr,bsnow,mkk)
c
       call frzh(pchomo,tttemp,qctemp,pcut,twodt,mkk)
c
c***********************************************************************
c          compute psacw - collection of cloud water by snow
c          (eq. a22, rh-1983). Limit to temperatures greater than thomo.
c***********************************************************************
c
       call eqa22(psacw,qctemp,qstemp,tttemp,termb,esc
     1           ,pcut,twodt,tzero,mkk)
c
c***********************************************************************
c          compute psmlt - melting of snow (eq. a25, rh-1983)
c***********************************************************************
c
       call eqa25(psmlt,qstemp,tttemp,terma,rbar2d,pi,snzero
     1           ,pcut,therco,tzero,heatlf,twodt,mkk)
c
c***********************************************************************
c          compute psmlti - melting of cloud ice (eq. a28, rh-1983)
c***********************************************************************
c
       call eqa28(psmlti,qitemp,tttemp,pcut,twodt,tzero,mkk)
c
c***********************************************************************
c          compute graupel interaction terms
c***********************************************************************
c
c***********************************************************************
c          compute pgacw - collection of cloud water by
c                  graupel (eq. a11 rh-1984) 
c***********************************************************************
c
       if(lgrpl)call eqa11g(pgacw,qgtemp,qctemp,tttemp,sfgraup,xa,iwx1
     1       ,pterm4,tzero,gnzero,pi,egc,twodt,pcut,abar,bgrp,mkk)
c
c***********************************************************************
c       update ice crystal concentration
c***********************************************************************
c
       call conice(nitemp,nuci,tttemp,qvtemp,qsatv,qsati,qctemp,qitemp
     1          ,rbar2d,nctemp,pr,appbpp,icon,pcut,heatlv,rvap,p00
c
c********************************************************************
c Addition by Adam Smith, 27 March 2008
c The CLEX-9 case requires a special setting for ice concentration.
c We add icase to the call statement to allow for the customization.
c********************************************************************
c
!     2          ,visair,cnzero,cbeta,twodt,tice,tzero,cimass,mkk)
     2          ,visair,cnzero,cbeta,twodt,tice,tzero,cimass,mkk,icase)
c
c****************************
c End of ajsmith4's addition
c****************************
c
c***********************************************************************
c       compute pristine ice fall speeds
c***********************************************************************
c
       call tgqi(falli,qitemp,rbar2d,pr,nitemp,pcut,p00,cimass
     1           ,dmax,mkk)
c
c***********************************************************************
c          compute piacw - collection of cloud water by
c                  pristine ice (eq. a27 Reisner 1998) 
c***********************************************************************
c
       call eqa27r(piacw,qctemp,qitemp,tttemp,falli,nitemp
     1      ,rbar2d,pi,eic,pcut,twodt,rholiq,tice,mkk)
c
c***********************************************************************
c          compute piacr - rate of cloud ice collecting 
c                  rain (eq. a7 rh-1984)
c***********************************************************************
c
       call eqa7g(piacr,qitemp,qrtemp,rbar2d,tttemp,sfrain,iwx1,nitemp
     1       ,pterm4,gm6,gm7,gm8,gm9,rnzero,rholiq,pi,eri,aa0,aa1,aa2
     2       ,aa3,twodt,pcut,tzero,mkk)
c
c***********************************************************************
c          compute psacr - rate of snow collecting 
c                  rain (eq. a8 rh-1984)
c***********************************************************************
c
       call eqa8g(psacr,qrtemp,qstemp,rbar2d,tttemp,fallr,falls
     1           ,xa,xb,iwx1,sfrain,sfsnow,pterm4,rnzero,snzero
     2           ,rholiq,pi,esr,twodt,pcut,tzero,mkk)
c
c***********************************************************************
c          compute pracs - rate of rain collecting 
c                  snow (eq. a9 rh-1984)
c***********************************************************************
c
       if(lgrpl)call eqa9g(pracs,qrtemp,qstemp,rbar2d,fallr,falls
     1       ,iwx1,xa,xb,sfsnow,sfrain,pterm4,rnzero,snzero,rhosno
     2       ,pi,esr,twodt,pcut,mkk)
c
c***********************************************************************
c          compute pgacr - collection of rain water by
c                  graupel (eq. a13 rh-1984) 
c***********************************************************************
c
       if(lgrpl)call eqa13g(pgacr,qrtemp,qgtemp,tttemp,rbar2d
     1       ,fallr,fallg,sfrain,sfgraup,iwx1,xa,xb,pterm4,rnzero
     2       ,gnzero,rholiq,pi,egr,twodt,pcut,tzero,mkk)
c
c***********************************************************************
c          compute pgmlt - melting of graupel 
c                  graupel (eq. a18 rh-1984) 
c***********************************************************************
c
       if(lgrpl)call eqa18g(pgmlt,qgtemp,tttemp,iwx1,termg,rbar2d
     1                     ,tzero,twodt,pcut,hkaolf,mkk)
c
c***********************************************************************
c          compute pgacrm - graupel accreting rain
c                  (eq. a21 rh-1984) 
c***********************************************************************
c
       if(lgrpl)call eqa21g(pgacrm,pgacr,tttemp,rbar2d,qgtemp
     1                     ,tzero,twodt,pcut,heatlf,cw,mkk)
c
c***********************************************************************
c          compute pgacwm- enhanced melting dut to the accretion  
c                  of cloud water (eq. a22 rh-1984) 
c***********************************************************************
c
       if(lgrpl)call eqa22g(pgacwm,pgacw,tttemp,rbar2d,qgtemp
     1                     ,tzero,twodt,pcut,heatlf,cw,mkk)
c
c***********************************************************************
c          adjust temperature and moisture fields over interior of 
c          the grid for gain/loss of ice/liquid due to melting
c          and freezing
c***********************************************************************
c
       call adjmlt(delta1,delta2,delta3,qctemp,qgtemp,qitemp,qrtemp
     1         ,qstemp,pchomo,pgacr,pgacrm,pgacw,pgacwm,pgmlt,piacr
     2         ,pracs,prhomo,psacw,psacr,psmlti,psmlt,piacw,tttemp
     3         ,condrt,lfocp,twodt,pcut,tzero,lgrpl,ary1d,maxpt1d
     4         ,kmax,maxvr1d,n1d,i1dflg,i1d,j1d,j,icomp,nipts,mkk)
c
c***********************************************************************
c          adjust supersaturation fields
c***********************************************************************
c
       call qsatvi(tttemp,pr,qsatv,qsati,lice,mkk)
c
      end if
c
c***********************************************************************
c          compute adjusted values of t and q (tadj,qadj) expected from 
c          warm microphysical processes vapor tendencies
c***********************************************************************
c
      iflag=1
      if(niter.gt.1)then
       do i=1,mkk
        tdum(i)=tttemp(i)
        qdum(i)=qvtemp(i)
       enddo
c
       do ni=1,niter
       
        call qsatvi(tdum,pr,qvsdum,qsidum,lice,mkk)
        call qtadj(tdum,qdum,pr,qvsdum,qsidum,qadjw,tadjw
     1   ,heatlv,heatls,cp,rvap,tzero,lice,iflag,pcut,mkk)
c
        do i=1,mkk
         tdum(i)=tadjw(i)
         qdum(i)=qadjw(i)
        enddo
       enddo
c
      else
c
       call qtadj(tttemp,qvtemp,pr,qsatv,qsati,qadjw,tadjw
     1   ,heatlv,heatls,cp,rvap,tzero,lice,iflag,pcut,mkk)
c
      end if
      
c
c***********************************************************************
c                        warm phase microphysics
c***********************************************************************
c
c***********************************************************************
c recompute slope for snow and rain to account for gain/losses 
c           due to melting
c***********************************************************************
c
      call slope(sfrain,sfraincm,sfsnow,sfgraup,qrtemp,qstemp,qgtemp
     1          ,rbar2d,sloper,slopes,slopeg,pcut,lice,lgrpl,mkk)
c
      if(lice)then
       do i=1,mkk
        if (sfsnow(i).le.0.0) then
          terma(i)=1.0
          xa(i)=1.0
          termb(i)=1.0
        else
c
c****************************************************************
c Addition by Adam Smith, 27 March 2008
c "terma" is used to calculate PSDEP and other microphysical 
c results.  For the Nov.11 and CLEX-9 cases, we use a simplified
c formula for calculating "terma".
c****************************************************************
c
!         terma(i)=0.65/(sfsnow(i)*sfsnow(i))
          if (icase.eq.1000 .or. icase.eq.1002 .or. icase.eq.1003) then
            terma(i)=2.70/(sfsnow(i)*sfsnow(i))
          else
            terma(i)=0.65/(sfsnow(i)*sfsnow(i))
          end if
c
c****************************
c End of ajsmith4's addition
c****************************
c
          temp=log(sfsnow(i))
          xa(i)=temp*ex2
          xa(i)=1.0/exp(xa(i))
          termb(i)=temp*ex3
          termb(i)=1.0/exp(termb(i))
        end if
       enddo
c
       do i=1,mkk
        xa(i)=xa(i)*pterm(i)*ex7
        xb(i)=rbar2d(i)*ex4
       enddo
c
       do i=1,mkk
        xa(i)=xa(i)*sqrt(xb(i))
c
c****************************************************************
c Addition by Adam Smith, 27 March 2008
c "terma" is used to calculate PSDEP and other microphysical 
c results.  For the Nov.11 and CLEX-9 cases, we use a simplified
c formula for calculating "terma".  For other cases, we need to
c modify "terma" using the standard formula.
c****************************************************************
c
        if (icase.ne.1000 .and. icase.ne.1002 .and. icase.ne.1003) then
          terma(i)=terma(i)+xa(i)
        end if
c
c****************************
c End of ajsmith4's addition
c****************************
c
        termb(i)=termb(i)*pterm4(i)*ex5
       enddo
c
c***********************************************************************
c compute terms needed in eqa's a17,a18,a19 of RH 1984
c***********************************************************************
c
       if(lgrpl)then
        do i=1,mkk
         if (sfgraup(i).le.0.0) then
          termg(i)=1.0
          xa(i)=1.0
         else
          termg(i)=0.78/(sfgraup(i)*sfgraup(i))
          temp=log(sfgraup(i))
          xa(i)=temp*ex2g
          xa(i)=1.0/exp(xa(i))
         end if
        enddo
c
        do i=1,mkk
         xa(i)=xa(i)*pterm(i)*ex7g
         xb(i)=rbar2d(i)*ex4g
        enddo
c
        do i=1,mkk
         xa(i)=xa(i)*sqrt(xb(i))
         termg(i)=2.0*pi*gnzero*(termg(i)+xa(i))
        enddo
c
       end if
      end if
c
c***********************************************************************
c          compute pcond - condensation and evaporation of
c          cloud water (eq. a6, rh-1983)
c***********************************************************************
c
      call eqa6(pcond,tttemp,qadjw,qvtemp,qctemp,qsatv,qsati,pr
     1         ,x6,x7,iwx1,heatlv,cp,twodt,pcut,j,mkk)
c
c***********************************************************************
c          compute prevp - evaporation of rain drops
c          (eq. a12, rh-1983)
c          budget terms ary1d(np,k,35) ary1d(np,k,45)
c***********************************************************************
c
      call eqa12(prevp,qrtemp,qvtemp,qadjw,qsatv,rbar2d,tttemp,sfrain
     1      ,pterm,apbp,nrtemp,pr,x6,x7,x8,x9,xa,xb,xc,iwx1,pi,rnzero
     2      ,apr,pcut,visair,gm3,twodt,ldrizzle,mkk)
c
      if(lice)then
c
c***********************************************************************
c          compute pmltse - evaporation of melting snow
c          (eq. a27, rh-1983)
c          budget terms ary1d(np,k,37) ary1d(np,k,47)
c***********************************************************************
c
       call eqa27(pmltse,qstemp,qvtemp,qsatv,qadjw,tttemp,terma,apbp
     1            ,rbar2d,snzero,tzero,pcut,twodt,j,mkk)
c
c***********************************************************************
c          compute pmltge - evaporation of melting graupel
c          (eq. a19, rh-1984)
c***********************************************************************
c
c note: pass in qsatv not qsati and qadjw not qadj
c
       if(lgrpl)call eqa19g(pmltge,qgtemp,tttemp,qvtemp,qsatv
     1       ,xa,iwx1,termg,apbp,qadjw,rbar2d,tzero,twodt,pcut,mkk)

      end if
c
c***********************************************************************
c compute weights to used to normalize microphysical rates
c The vapor demand cannot exceed dqadj over a time step. 
c***********************************************************************
c
      call nrmtqw(qvtemp,qgtemp,qrtemp,qstemp,qadjw,qctemp,tttemp
     1            ,pcond,prevp,pmltge,pmltse,condrt,nctemp,nrtemp
     2            ,nntemp,pcut,twodt,lvocp,lgrpl,ary1d,maxpt1d
     3            ,kmax,maxvr1d,n1d,i1dflg,i1d,j1d,j,ldrizzle,th2t
c-->Chris Golaz modification: added icase
!     4            ,heatlv,cp,rvap,qsatv,mkk)
     4            ,heatlv,cp,rvap,qsatv,mkk,icase)
c<--end of modification
c
c***********************************************************************
c          adjust supersaturation fields
c***********************************************************************
c
      call qsatvi(tttemp,pr,qsatv,qsati,lice,mkk)
c
      call slope(sfrain,sfraincm,sfsnow,sfgraup,qrtemp,qstemp,qgtemp
     1          ,rbar2d,sloper,slopes,slopeg,pcut,lice,lgrpl,mkk)
c
      if(lice)then
       do i=1,mkk
        if (sfsnow(i).le.0.0) then
          terma(i)=1.0
          xa(i)=1.0
          termb(i)=1.0
        else
c
c****************************************************************
c Addition by Adam Smith, 27 March 2008
c "terma" is used to calculate PSDEP and other microphysical 
c results.  For the Nov.11 and CLEX-9 cases, we use a simplified
c formula for calculating "terma".
c****************************************************************
c
!         terma(i)=0.65/(sfsnow(i)*sfsnow(i))
          if (icase.eq.1000 .or. icase.eq.1002 .or. icase.eq.1003) then
            terma(i)=2.70/(sfsnow(i)*sfsnow(i))
          else
            terma(i)=0.65/(sfsnow(i)*sfsnow(i))
          end if
c
c****************************
c End of ajsmith4's addition
c****************************
c
          temp=log(sfsnow(i))
          xa(i)=temp*ex2
          xa(i)=1.0/exp(xa(i))
          termb(i)=temp*ex3
          termb(i)=1.0/exp(termb(i))
        end if
       enddo
c
       do i=1,mkk
        xa(i)=xa(i)*pterm(i)*ex7
        xb(i)=rbar2d(i)*ex4
       enddo
c
       do i=1,mkk
        xa(i)=xa(i)*sqrt(xb(i))
c
c****************************************************************
c Addition by Adam Smith, 27 March 2008
c "terma" is used to calculate PSDEP and other microphysical 
c results.  For the Nov.11 and CLEX-9 cases, we use a simplified
c formula for calculating "terma".  For other cases, we need to
c modify "terma" using the standard formula.
c****************************************************************
c
        if (icase.ne.1000 .and. icase.ne.1002 .and. icase.ne.1003) then
          terma(i)=terma(i)+xa(i)
        end if
c
c****************************
c End of ajsmith4's addition
c****************************
c
        termb(i)=termb(i)*pterm4(i)*ex5
       enddo
c
c***********************************************************************
c compute terms needed in eqa's a17,a18,a19 of RH 1984
c***********************************************************************
c
       if(lgrpl)then
        do i=1,mkk
         if (sfgraup(i).le.0.0) then
          termg(i)=1.0
          xa(i)=1.0
         else
          termg(i)=0.78/(sfgraup(i)*sfgraup(i))
          temp=log(sfgraup(i))
          xa(i)=temp*ex2g
          xa(i)=1.0/exp(xa(i))
         end if
        enddo
c
        do i=1,mkk
         xa(i)=xa(i)*pterm(i)*ex7g
         xb(i)=rbar2d(i)*ex4g
        enddo
c
        do i=1,mkk
         xa(i)=xa(i)*sqrt(xb(i))
         termg(i)=2.0*pi*gnzero*(termg(i)+xa(i))
        enddo
       end if
      end if
c
c***********************************************************************
c          Ice-phase Microphysics
c***********************************************************************
c
      if (lice) then
c
c***********************************************************************
c          compute adjusted values of t and q (tadj,qadj) expected from 
c          ice-phase microphysics 
c***********************************************************************
c
       iflag=2
       if(niter.gt.1)then
        do i=1,mkk
         tdum(i)=tttemp(i)
         qdum(i)=qvtemp(i)
        enddo
c
        do ni=1,niter
       
         call qsatvi(tdum,pr,qvsdum,qsidum,lice,mkk)
         call qtadj(tdum,qdum,pr,qvsdum,qsidum,qadj,tadj
     1   ,heatlv,heatls,cp,rvap,tzero,lice,iflag,pcut,mkk)
c
        do i=1,mkk
         tdum(i)=tadj(i)
         qdum(i)=qadj(i)
        enddo
       enddo
      else
       call qtadj(tttemp,qvtemp,pr,qsatv,qsati,qadj,tadj,heatlv
     1   ,heatls,cp,rvap,tzero,lice,iflag,pcut,mkk)
      end if
c
c***********************************************************************
c          compute ice concentration (eq. a13, rh-1983)
c***********************************************************************
c
       call conice(nitemp,nuci,tttemp,qvtemp,qsatv,qsati,qctemp,qitemp
     1          ,rbar2d,nctemp,pr,appbpp,icon,pcut,heatlv,rvap,p00
c
c********************************************************************
c Addition by Adam Smith, 27 March 2008
c The CLEX-9 case requires a special setting for ice concentration.
c We add icase to the call statement to allow for the customization.
c********************************************************************
c
!     2          ,visair,cnzero,cbeta,twodt,tice,tzero,cimass,mkk)
     2          ,visair,cnzero,cbeta,twodt,tice,tzero,cimass,mkk,icase)
c
c****************************
c End of ajsmith4's addition
c****************************
c
c***********************************************************************
c          compute pint - initiation of cloud ice (eq. a15,rh-1983)
c          budget ary1d(np,k,30),ary1d(np,k,40)
c***********************************************************************
c
       call eqa15(pint,tttemp,qvtemp,qsati,qadj,rbar2d,nitemp,tice 
     1            ,pcut,twodt,icon,mkk)
c
c***********************************************************************
c          compute pdepi - depositional growth of cloud ice
c          (eq. a18, rh-1983)
c          budget ary1d(np,k,31),ary1d(np,k,41)
c***********************************************************************
c
       call eqa18(pdepi,qvtemp,qitemp,qsati,qadj,tttemp,rbar2d,nitemp
     1            ,appbpp,tice,twodt,pcut,j,icon,mkk)
c
c***********************************************************************
c          compute psdep - depositional growth of snow
c          (eq. a26, rh-1983)
c          budget ary1d(np,k,32),ary1d(np,k,42)
c***********************************************************************
c
       call eqa26(psdep,qvtemp,qstemp,qsati,qadj,terma,appbpp,tttemp
     1            ,rbar2d,snzero,tzero,pcut,twodt,mkk)
c
c***********************************************************************
c          compute pgdep - depositional growth of graupel
c          (eq. a17, rh-1984)
c***********************************************************************
c
       if(lgrpl)call eqa17g(pgdep,qgtemp,tttemp,qvtemp,qsati,xa
     1       ,iwx1,termg,appbpp,qadj,rbar2d,tzero,twodt,pcut,mkk)
c
c***********************************************************************
c   if total vapor tendency (qsum) exceeds available vapor (qvtemp-qadj)
c   for time step then normalize the demand. Each tendency has already
c   been confined to be no greater than dqadj 
c***********************************************************************
c
       call nrmtqi(qvtemp,qitemp,qgtemp,qstemp,qadj,tttemp,pgdep
     1            ,pdepi,psdep,pint,condrt,pcut,twodt,lsocp
     2            ,tzero,lgrpl,ary1d,maxpt1d,kmax,maxvr1d,n1d,i1dflg
     3            ,i1d,j1d,j,icomp,th2t,nitemp,nipts,mkk)
c
c***********************************************************************
c    Compute conversion terms using updated microphysical field
c    These conversions will not impact temperature or vapor  until
c    the next time step
c***********************************************************************
c
      end if
c
      call slope(sfrain,sfraincm,sfsnow,sfgraup,qrtemp,qstemp,qgtemp
     1          ,rbar2d,sloper,slopes,slopeg,pcut,lice,lgrpl,mkk)
c
c***********************************************************************
c          compute praut - autoconversion of cloud water
c          icond=1 Kessler (eq. a7, rh-1983)
c          icond=2 Manton and Cotton (1977)
c          icond=3 Khairoutdinov and Kogan (2000)
c          icond=4 No autoconversion
c***********************************************************************
c
      call eqa7(praut,qctemp,rbar2d,nctemp,praut1,praut2,pi
     1          ,rholiq,icond,twodt,pcut,mkk)
c
c***********************************************************************
c          compute pracw - collection of cloud water by rain
c          water (eq. a9, rh-1983)
c***********************************************************************
c
      call eqa9(pracw,qctemp,qrtemp,pterm4,sfrain,xa,rbar2d,nctemp
     1         ,iwx1,pi,erc,rnzero,aa0,aa1,aa2,aa3,gm3,gm4,gm5,gm6
     2         ,pcut,twodt,ldrizzle,mkk)
c
      if(lice)then
c
       do i=1,mkk
        if (sfsnow(i).le.0.0) then
          terma(i)=1.0
          xa(i)=1.0
          termb(i)=1.0
        else
c
c****************************************************************
c Addition by Adam Smith, 27 March 2008
c "terma" is used to calculate PSDEP and other microphysical 
c results.  For the Nov.11 and CLEX-9 cases, we use a simplified
c formula for calculating "terma".
c****************************************************************
c
!         terma(i)=0.65/(sfsnow(i)*sfsnow(i))
          if (icase.eq.1000 .or. icase.eq.1002 .or. icase.eq.1003) then
            terma(i)=2.70/(sfsnow(i)*sfsnow(i))
          else
            terma(i)=0.65/(sfsnow(i)*sfsnow(i))
          end if
c
c****************************
c End of ajsmith4's addition
c****************************
c
          temp=log(sfsnow(i))
          xa(i)=temp*ex2
          xa(i)=1.0/exp(xa(i))
          termb(i)=temp*ex3
          termb(i)=1.0/exp(termb(i))
        end if
       enddo
c
       do i=1,mkk
        xa(i)=xa(i)*pterm(i)*ex7
        xb(i)=rbar2d(i)*ex4
       enddo
c
       do i=1,mkk
        xa(i)=xa(i)*sqrt(xb(i))
c
c****************************************************************
c Addition by Adam Smith, 27 March 2008
c "terma" is used to calculate PSDEP and other microphysical 
c results.  For the Nov.11 and CLEX-9 cases, we use a simplified
c formula for calculating "terma".  For other cases, we need to
c modify "terma" using the standard formula.
c****************************************************************
c
        if (icase.ne.1000 .and. icase.ne.1002 .and. icase.ne.1003) then
          terma(i)=terma(i)+xa(i)
        end if
c
c****************************
c End of ajsmith4's addition
c****************************
c
        termb(i)=termb(i)*pterm4(i)*ex5
       enddo
c
c***********************************************************************
c compute terms needed in eqa's a17,a18,a19 of RH 1984
c***********************************************************************
c
       if(lgrpl)then
        do i=1,mkk
         if (sfgraup(i).le.0.0) then
          termg(i)=1.0
          xa(i)=1.0
         else
          termg(i)=0.78/(sfgraup(i)*sfgraup(i))
          temp=log(sfgraup(i))
          xa(i)=temp*ex2g
          xa(i)=1.0/exp(xa(i))
         end if
        enddo
c
        do i=1,mkk
         xa(i)=xa(i)*pterm(i)*ex7g
         xb(i)=rbar2d(i)*ex4g
        enddo
c
        do i=1,mkk
         xa(i)=xa(i)*sqrt(xb(i))
         termg(i)=2.0*pi*gnzero*(termg(i)+xa(i))
        enddo
       end if
c
c***********************************************************************
c
       if(lgrpl)then
        call tgqg(fallg,qgtemp,sfgraup,pterm4,pcut,abar,bgrp,mkk)
c
        call tgqs(falls,qstemp,sfsnow,pterm4,pcut,aprpr,bsnow,mkk)
       end if
c
c***********************************************************************
c          compute pconv - conversion of cloud ice to snow
c          (eq. a19, rh-1983)
c***********************************************************************
c
       call eqa19(pconv,qitemp,rbar2d,nitemp,tttemp,pcut,cimass
     1             ,twodt,mkk)
c
c***********************************************************************
c          compute pwacs - loss of snow due to collisions with cloud
c                  water (eq. a10 rh-1984) (note: might ignore this one)
c***********************************************************************
c
       if(lgrpl)call eqa10g(pwacs,qstemp,qctemp,rbar2d,tttemp
     1       ,sfsnow,xa,nctemp,iwx1,pterm4,tzero,snzero,rhosno
     2       ,pi,esc,twodt,pcut,aprpr,bsnow,mkk)
c
c***********************************************************************
c          compute psaci - collection of cloud ice by snow
c          (eq. a21, rh-1983)
c***********************************************************************
c
       call eqa21(psaci,qitemp,qstemp,sfsnow,termb
     1            ,esi,pcut,twodt,mkk)
c
c***********************************************************************
c          compute praci - collection of cloud ice by rain
c          (eq. a5, rh-1984)
c***********************************************************************
c
       if(lgrpl)call eqa5g(praci,qitemp,qrtemp,tttemp,iwx1,gm3,gm4
     1       ,gm5,gm6,pterm4,rnzero,pi,eri,aa0,aa1,aa2,aa3,sfrain
     2       ,twodt,pcut,tzero,mkk)
c
c***********************************************************************
c          compute pgaci - collection of cloud ice by graupel
c          (eq. a12, rh-1984)
c***********************************************************************
c
       if(lgrpl)call eqa12g(pgaci,qitemp,qgtemp,tttemp,sfgraup,xa,iwx1
     1       ,pterm4,tzero,gnzero,pi,egi,twodt,pcut,abar,bgrp,mkk)
c
c***********************************************************************
c          compute pgacs - collection of snow by graupel
c          (eq. a14, rh-1984)
c***********************************************************************
c
       if(lgrpl)call eqa14g(pgacs,qstemp,qgtemp,tttemp,rbar2d
     1       ,falls,fallg,sfsnow,sfgraup,iwx1,xa,xb,pterm4
     2       ,snzero,gnzero,rhosno,pi,egs,twodt,pcut,tzero,mkk)
c
c***********************************************************************
c          compute pgshr - shedding of accreted water 
c          (eq. a20, rh-1984)
c***********************************************************************
c
       if(lgrpl)call eqa20g(pgshr,pgacr,pgacw,tttemp,rbar2d,tzero
     2                     ,twodt,pcut,mkk)
c
      end if
c
c***********************************************************************
c       compute collection tendencies for all terms not impacting
c       temperature or vapor
c***********************************************************************
c
      call nrmcol(tttemp,qctemp,qitemp,qgtemp,qrtemp,qstemp,pconv,pgaci
     1            ,pgacs,praci,pracw,praut,psaci,pwacs,rbar2d,nctemp
     2            ,nrtemp,nntemp,pcut,twodt,lice,lgrpl,tttemp,rvc,ary1d
     3            ,maxpt1d,kmax,maxvr1d,n1d,i1dflg,i1d,j1d,j,icomp
     4            ,ldrizzle,nitemp,nipts,mkk)
c
c***********************************************************************
c     finalize fields
c***********************************************************************
c
c***********************************************************************
c          convert to potential temperature
c***********************************************************************
c
      do i=1,mkk
        tttemp(i)=tttemp(i)/th2t(i)
      enddo
c
c***********************************************************************
c          add in heating and moistening tendencies 
c***********************************************************************
c
      ilen=0
      do i=1,nipts
       iindex=icomp(i)
       do k=1,nkpts(iindex)
        kindex=kcomp(iindex,k)
        ilen=ilen+1
          th3(iindex,j,kindex)=tttemp(ilen)
          qv3(iindex,j,kindex)=max(pcut,qvtemp(ilen)) 
       enddo
      enddo
c
      ilen=0
      do i=1,nipts
       iindex=icomp(i)
       do k=1,nkpts(iindex)
        kindex=kcomp(iindex,k)
        ilen=ilen+1
! Michael Falk is messing with things
!        qc3(iindex,j,kindex)= qctemp(ilen)
        qc3(1,1,kindex)= qctemp(ilen)
! eMFc

        qr3(iindex,j,kindex)= qrtemp(ilen)
        qs3(iindex,j,kindex)= qstemp(ilen)
        qi3(iindex,j,kindex)= qitemp(ilen)
        if(qi3(iindex,j,kindex).eq.0.)then
          ni3(iindex,j,kindex)=0.
        else
          ni3(iindex,j,kindex)= nitemp(ilen)
        end if
        qg3 (iindex,j,kindex)= qgtemp(ilen)

        cond(iindex,j,kindex)= condrt(ilen)

c Addition by Adam Smith, 08 August 2006
c Adding snow slope factor and fall velocity to conversion
c Also adding variables that affect mc terms in statistical budgets
c Michael Falk started making comments here on 25 Aug 2006
       snowslope(iindex,j,kindex) = sfsnow(ilen)
       snowv(iindex,j,kindex) = falls(ilen)
c        pcond3(iindex,j,kindex)  = pcond(ilen)
c        psmlti3(iindex,j,kindex) = psmlti(ilen)
c        psacw3 (iindex,j,kindex) = psacw (ilen)
c        pgacw3 (iindex,j,kindex) = pgacw (ilen)
c        piacw3 (iindex,j,kindex) = piacw (ilen)
c        pchomo3(iindex,j,kindex) = pchomo(ilen)
c        praut3 (iindex,j,kindex) = praut (ilen)
c        pracw3 (iindex,j,kindex) = pracw (ilen)
c        pdepi3 (iindex,j,kindex) = pdepi (ilen)
c        pint3  (iindex,j,kindex) = pint  (ilen)
c        pgdep3 (iindex,j,kindex) = pgdep (ilen)
c        pconv3 (iindex,j,kindex) = pconv (ilen)
c        psaci3 (iindex,j,kindex) = psaci (ilen)
c        pgaci3 (iindex,j,kindex) = pgaci (ilen)
c        praci3 (iindex,j,kindex) = praci (ilen)
c        prevp3 (iindex,j,kindex) = prevp (ilen)
c        psdep3 (iindex,j,kindex) = psdep (ilen)
c        pmltge3(iindex,j,kindex) = pmltge(ilen)
c        pgmlt3 (iindex,j,kindex) = pgmlt (ilen)
c        psmlt3 (iindex,j,kindex) = psmlt (ilen)
c        pgacrm3(iindex,j,kindex) = pgacrm(ilen)
c        pgacwm3(iindex,j,kindex) = pgacwm(ilen)
c        pracs3 (iindex,j,kindex) = pracs (ilen)
c        pgshr3 (iindex,j,kindex) = pgshr (ilen)
c        pgacr3 (iindex,j,kindex) = pgacr (ilen)
c        psacr3 (iindex,j,kindex) = psacr (ilen)
c        piacr3 (iindex,j,kindex) = piacr (ilen)
c        prhomo3(iindex,j,kindex) = prhomo(ilen)
c        pgacs3 (iindex,j,kindex) = pgacs (ilen)
c        pmltse3(iindex,j,kindex) = pmltse(ilen)
c        pwacs3 (iindex,j,kindex) = pwacs (ilen)
c
cc End of ajsmith4's addition
c eMFc

       enddo
      enddo

      if(ldrizzle)then
       ilen=0
       do i=1,nipts
        iindex=icomp(i)
        do k=1,nkpts(iindex)
         kindex=kcomp(iindex,k)
         ilen=ilen+1
c-->Chris Golaz modification: 
c   Need to convert nc3, nr3, ncn3 from concentration (#/cm3)
c   to mixing ratio units ( x 1.e-6 #/kg ).

!         nc3(iindex,j,kindex) = nctemp(ilen)
!         nr3(iindex,j,kindex) = nrtemp(ilen)
!         ncn3(iindex,j,kindex)= nntemp(ilen)

         nc3(iindex,j,kindex) = nctemp(ilen) /
     .      rbm(iindex,j,kindex)
         nr3(iindex,j,kindex) = nrtemp(ilen) / 
     .      rbm(iindex,j,kindex)
         ncn3(iindex,j,kindex)= nntemp(ilen) / 
     .      rbm(iindex,j,kindex)

c<--end of modification
        enddo
       enddo
      end if
c
c***********************************************************************
c zero out fields in the damping layer
c***********************************************************************
c
      do k=1,nrdamp
       do i=1,nipts
        qc3(icomp(i),j,k) = 0.
        qr3(icomp(i),j,k) = 0.
        qs3(icomp(i),j,k) = 0.
        qi3(icomp(i),j,k) = 0.
        qg3(icomp(i),j,k) = 0.
        nc3(icomp(i),j,k) = 0.
        nr3(icomp(i),j,k) = 0.
        ni3(icomp(i),j,k) = 0.
        qv3(icomp(i),j,k)=max(pcut,qv3(icomp(i),j,k) )
       enddo
      enddo
c
c***********************************************************************
c
      return
      end subroutine adjtq
