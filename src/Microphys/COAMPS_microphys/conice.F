      subroutine conice(ni,nuci,t,qv,qsatv,qsati,qc,qi,rho,nc
     1          ,pr,appbpp,icon,pcut,heatlv,rvap,p00,visair
c
c*******************************************
c Modification by Adam Smith, 10 April 2007
c Adding icase to call statement
c*******************************************
c
c     2          ,cnzero,cbeta,delt,tice,tzero,cimass,len)
     2          ,cnzero,cbeta,delt,tice,tzero,cimass,len,icase)
c
c********************************
c End of ajsmith4's modification
c********************************
c
c
#include "conice.prol"
      implicit none
c
c***********************************************************************
c           parameters:
c***********************************************************************
c
      integer len
      integer icon
      integer icase
      intent(in) icase
c
      real appbpp(len)
      real cbeta
      real cimass
      real cnzero
      real delt
      real heatlv
      real nc    (len)
      real ni    (len)
      real nuci  (len)
      real p00
      real qc    (len)
      real qi    (len)
      real qv    (len)
      real qsatv (len)
      real qsati (len)
      real rho   (len)
      real pcut
      real pr    (len)
      real rvap
      real t     (len)
      real tice
      real tzero
      real visair 
c
c***********************************************************************
c local variables
c***********************************************************************
c
c
c***********************************************************************
c          local variables and dynamic storage:
c***********************************************************************
c
      real aerorad
      real cnn   (len)
      real coef
      real conci (len)
      real connc (len)
      real convd
      real ct    (len)
      real dc    (len)
      real dt
      real f1
      real f2
      real f3
      real ft
      real ftp
      real kapp
      real kappa (len)
      real kbolt
      real knud
      real lambda
      real l2m3
      real mc
      real n12
      real pi
      real qsi
      real rvcof
      real rolv
      real psi
      real ss    (len)
      real t0
      real tcon
      real tdep
      real tmin
      real tmp1  (len)
      real tmp2  (len)
      real tmp3  (len)
c
      integer i
c
      logical lcotton
c
c***************************************************************`
c Addition by Adam Smith, 06 April 2007
c Adding a "fudge" factor to modify Fletcher/Cooper scheme.
c We have found that the M-PACE and CLEX-9: Oct.14 cases produce
c   more accurate ice when we use the fudgefactor.
c Made this a if .. then instead -dschanen 27 Mar 08
c***************************************************************
c
      real fudgefactor
c
c****************************
c End of ajsmith4's addition
c****************************
c
c***********************************************************************
c          end of definitions
c***********************************************************************
c
      lcotton=.false.
      l2m3=1000.
      tcon=tzero-2.0
      tdep=tzero-5.0
      rolv=rvap/heatlv
      pi=4.0*atan(1.0)
      t0=293.15
      lambda=6.6e-8
      aerorad=1.e-7
      kbolt=1.38e-23
c
      if ( icase == 75 ) then
!     if ( .false. ) then
        ! For MPACE A
        fudgefactor = 0.01
      elseif ( icase == 1002) then
        ! For CLEX-9: Oct.14
        fudgefactor = 0.06
      else ! or don't cheat
        fudgefactor = 1.
      end if
c
c***********************************************************************
c aerosol conductivity (kapp) taken from value for clay (Sellers 1965)
c***********************************************************************
c
      kapp=2.92
c
c***********************************************************************
c coef=1./(den_water*pi*4/3)
c***********************************************************************
c
      coef=2.3873241e-4

c
      do i=1,len
        conci(i)=0.
        connc(i)=0.
        tmp1(i)=0.
        tmp2(i)=0.
        tmp3(i)=0.
        nuci(i)=0.
      enddo
c
c***********************************************************************
c Fletcher formulization
c***********************************************************************
c
      if(icon.eq.1)then
c
       do i=1,len
        if (t(i).le.tice) then
          tmin=max(t(i),243.)
          ni(i)=cnzero*exp(cbeta*(tzero-tmin))
        endif
       enddo
c
      elseif(icon.eq.2)then
c
c***********************************************************************
c Diffusiophoresis, Thermophoresis and Brownian motion
c        (Walko, 1995; eqs 52-59)
c***********************************************************************
c
c***********************************************************************
c Deposition/Condensation freezing Meyers et al (1992)
c***********************************************************************
c
c***********************************************************************
c Compute Deposition/Condensation freezing term
c and add to total nucleation only if it exceeds the
c existing crystal concentration
c***********************************************************************
c
       do i=1,len
        qsi=qv(i)/qsati(i)-1.
        if(t(i).le.tdep.and.qsi.gt.0.)then
          convd=min(1.e10,l2m3*exp(-.639+12.96*qsi) )
        else 
          convd=0.
        endif
c
        nuci(i)=max(0.,convd-ni(i))
       enddo
c
c***********************************************************************
c contact freezing terms (Meyers et al (1992);Walko et al 1995)
c***********************************************************************
c
        if (t(i).le.tcon) then
          tmin=max(t(i),243.)
          connc(i)=min(1.e10,l2m3*exp(-2.8+.262*(tzero-tmin) ) )
        endif
c
c***********************************************************************
c compute thermal conductivity and supersaturation with respect to water
c kappa from eq 13-16 Pruppacher and Klett 1978
c***********************************************************************
c
       do i=1,len
        kappa(i)=.04184*(5.69+.017*(t(i)-273.16))
        ss(i)   =qsatv(i)/qv(i)-1.0
       enddo
c
c***********************************************************************
c compute droplet temperature (see Beyers eq. 5.22 and 5.12)
c***********************************************************************
c
       do i=1,len
        if(qc(i).gt.pcut)then
         ct(i)=heatlv*ss(i)*appbpp(i)/kappa(i)+t(i)
        else
         ct(i)=t(i)
        endif
       enddo
c
c***********************************************************************
c compute mean droplet diameter (in microns)
c***********************************************************************
c
       do i=1,len
        if(qc(i).gt.pcut.and.nc(i).gt.0.)then
         dc(i)=.124*1.e6*(qc(i)*rho(i)/nc(i))**(1./3.)
        else
         dc(i)=0.
        endif
       enddo
c
c***********************************************************************
c note: These formulas require mean droplet diameter (m)
c       convert cloud radius from microns to m and convert to diameter
c       2.e-6*dc(i)
c***********************************************************************
c
       do i=1,len
        if (t(i).le.tcon.and.qc(i).gt.pcut) then

         f1=2.0*pi*2.e-6*dc(i)*connc(i)*nc(i)
         f2=kappa(i)*(t(i)-ct(i))/pr(i)
         knud=lambda*t(i)*p00/(t0*pr(i)*aerorad)
         ft=.4*(1.+1.45*knud+.4*knud*exp(-1./knud))
     1     *(kappa(i)+2.5*knud*kapp)
     1     /( (1.+3.*knud)*(2.*kappa(i)+5.0*kapp*knud+kapp) )
         psi=kbolt*(ct(i)-273.16)*(1.+knud)/(6.*pi*aerorad*visair)
c
         tmp1(i)=max(0.,-f1*f2*rolv*t(i)/rho(i))
         tmp2(i)=max(0.,f1*f2*ft/rho(i))
         tmp3(i)=max(0.,f1*psi) 
c
        endif
       enddo
c
       do i=1,len
         if (t(i).le.tcon) then
          nuci(i)=nuci(i)+delt*(tmp1(i)+tmp2(i)+tmp3(i))
         endif
       enddo
      endif
c
       do i=1,len
        ni(i)=max(0.,ni(i)+nuci(i))
       enddo
c
      if(icon.eq.3)then
c
c***********************************************************************
c Hobbs and Rangno (1985) eq 3. converted to /m**3
c diameter is in microns
c***********************************************************************
c
c***********************************************************************
c compute mean droplet diameter (in microns)
c***********************************************************************
c
       do i=1,len
        if(qc(i).gt.pcut.and.nc(i).gt.0.)then
         rvcof=rho(i)*coef
c        dc(i)=.124*1.e6*(qc(i)*rho(i)/nc(i))**(1./3.)
         dc(i)=1.e6*(qc(i)*rvcof/nc(i))**(.33333)
        else
         dc(i)=0.
        endif
       enddo
c
       do i=1,len
        nuci(i)=max(0.,1000.*(dc(i)/18.5)**8.4-ni(i))
       enddo
c
      endif
c
c***********************************************************************
c Cooper formula (1986): convert to #/m**3
c***********************************************************************
c
      if(icon.eq.4)then
        do i=1,len
         qsi=qv(i)/qsati(i)-1.
         if(t(i).lt.tdep.and.qsi.gt.0.)then
          tmin=max(t(i),246.)
          ni(i)=.005*l2m3*exp( 0.304*(tzero-tmin) )
         endif
        enddo
      endif
c
c***********************************************************************
c Use Cooper at warm temps and Fletcher at cold temps
c***********************************************************************
c
      if(icon.eq.5)then
       do i=1,len
        if (t(i).le.tice) then
          tmin=max(t(i),243.)
!
! Michael Falk, 15 Oct 2007
! Use the second line instead of the first for MPACE-A
c
c********************************************************
c Modification by Adam Smith, 28 March 2008
c Re-commenting is no longer necessary for this section,
c because a fudgefactor is applied to all cloud cases.
c********************************************************
c
!         ni(i)=cnzero*exp(cbeta*(tzero-tmin))
          ni(i)=fudgefactor*cnzero*exp(cbeta*(tzero-tmin))
c
c********************************
c End of ajsmith4's modification
c********************************
c
! eMFc
        endif
       enddo
c
       do i=1,len
        qsi=qv(i)/qsati(i)-1.
         if(t(i).lt.tdep.and.qsi.gt.0.)then
          tmin=max(t(i),246.)
! Michael Falk, 15 October 2007
! Use the second and third lines instead of the first for MPACE-A
c
c********************************************************
c Modification by Adam Smith, 28 March 2008
c Re-commenting is no longer necessary for this section,
c because a fudgefactor is applied to all cloud cases.
c********************************************************
c
!         ni(i)=max(ni(i), .005*l2m3*exp(0.304*(tzero-tmin)) )
          ni(i)=max(ni(i), 
     .     fudgefactor*.005*l2m3*exp(0.304*(tzero-tmin)) )
c
c********************************
c End of ajsmith4's modification
c********************************
c
! eMFc
         endif
! Michael Falk, 11 Oct 2007, added this per AJS and VL.
! This fixes the issue where there could be ice particles outside of cloud.
         if (qc(i) .le. pcut) ni(i) = 0.0
! eMFc
       enddo
      endif
c
c***********************************************************************
c  special check: Due to the nozero fallspeed of ice, pristine crystals 
c                 may exist in regions where the diagnosed value is zero.
c                 Prevent this from occurring by using the mean cyrstal
c                 mass to diagnose a new concentration.
c***********************************************************************
c
      do i=1,len
       if(qi(i).gt.pcut.and.ni(i).le.0.)then
        ni(i)=qi(i)*rho(i)/cimass
       endif
      enddo
c
c     do i=1,len
cxx    ni(i)=min(1.e10,max(0.,ni(i)+nuci(i)) )
c     enddo
c
c***********************************************************************
c
      return
      end
