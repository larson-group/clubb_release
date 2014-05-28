      subroutine esat_new(tp,vp,len,l_ice)
c
c SCCS IDENTIFICATION:  %W% %G%
c                                %U% %P% 
c
      implicit none
c
c***********************************************************************
c           parameters:
c***********************************************************************
c
      integer len
c
      real tp     (len)
      real vp     (len)
c      real vd     (len)
c
      logical l_ice
c                                 
c  =================================================================
c  +                                                               +
c  +             CONSTRUCTION OF THE ES TABLE                      +
c  +                                                               +   
c  + THIS TABLE IS CONSTRUCTED FROM ES EQUATIONS FROM THE          +           
c  + SMITHSONIAN TABLES.  THE ES INPUT IS COMPUTED FROM VALUES     +           
c  + (IN ONE-TENTH OF A DEGREE INCREMENTS) OF ES OVER ICE          +
c  + FROM -153C TO 0C AND VALUES OF ES OVER WATER FROM -50C TO 102C+
c  +   NOTE: ALL ES COMPUTATION IS DONE IN MICROBARS AND           +   
c  +          AND THEN CONVERTED TO millibars or Pascals.          +      
c  =================================================================
c                 
c #TRANSITIONED q_cpu (OCT 15, 2001) - NEW GFDL EQUATION FROM GFDL
c
c***********************************************************************
c          local variables and dynamic storage:
c***********************************************************************
c
      integer ntot
      parameter(ntot=2751)
      real vpsat(ntot)      
      real visat(ntot)      
      real varvp(ntot)      
      real varvi(ntot)
c      
      save vpsat,visat
      save varvp,varvi
c                                
      real ESBASW,TBASW,ESBASI,TBASI
      real TEM,AA,B,C,D,E
      real RD
      real ESH2O,WICE,WH2O
      real ESICE,ESBLD
      real VAH2O
      real VAICE
      real VABLD
      real XLH2O
      real XLICE
      real XLBLD
c
      integer  i
      integer ic
      integer ii
      integer  j
      real qcon,rq,vps,var
      real rdf
      real t1,qqq,tqs,temc
      real temp
      real tcutoff
      real varest
      logical lcutoff
      logical lprint
      logical ltest
c
      logical l_xl
      logical l_interp 
      logical l_Pascals 
c
      logical lfirst 
      data lfirst/.true./             
      save lfirst
c
c
      lprint=.false.
c      lprint=.true.
      ltest=.false.
      lcutoff=.false.
      l_interp=.true.
      l_Pascals=.true.
c
      qcon=0.62197
      rq=1.0/(1.0+qcon)
      temp= 1.0/1.622
      rdf=287.04/qcon
c     if lcutoff true then use a
c     temperature cutoff for vapor pressure over water, set at -40C 
c      tcutoff=233.16
c     temperature cutoff for vapor pressure over water, set at -50C 
      tcutoff=223.16
c
c***********************************************************************
c  create lookup table for saturation vapor pressure and it's derivative
c***********************************************************************
c
      if(lfirst)then
c
        ESBASW = 1013246.0
        TBASW =     373.16
        ESBASI =    6107.1
        TBASI =     273.16
        RD    =     287.04
c
        if(lprint) then
          write(6,10)
          write(6,11)
        endif
c
c  COMPUTE ES OVER ICE BETWEEN -153 C AND 0 C AND 
c  COMPUTE ES OVER WATER - SEE SMITHSONIAN METEOROLOGICAL TABLES PAGE 350.
c  BASED ON GOFF-GRATCH FORMULATION, 1946 WITH OLD KELVIN SCALE,
c   TK = TC + 273.16.
c
        do 100 i=1,ntot
c
c         SATURATION VAPOR PRESSURE OVER ICE
          TEM = 120.16+0.1*FLOAT(i-1)
          AA  = -9.09718 *(TBASI/TEM-1.0)
          B   = -3.56654 *ALOG10(TBASI/TEM)
          C   =  0.876793*(1.0-TEM/TBASI)
          E   =  ALOG10(ESBASI)
          ESICE= 10.0**(AA+B+C+E)
c         convert to millibars
c          ESICE= 0.001*ESICE
          if(l_Pascals) then
c           convert to Pascals
            ESICE= 0.1*ESICE
          endif
c
c         SATURATION VAPOR PRESSURE OVER WATER
          AA  = -7.90298*(TBASW/TEM-1.0)
          B   =  5.02808*ALOG10(TBASW/TEM)
          C   = -1.3816E-07*
     &           (10.0**((1.0-TEM/TBASW)*11.344)-1.0)
          D   =  8.1328E-03*
     &           (10.0**((TBASW/TEM-1.0)*(-3.49149))-1.0)
          E   =  ALOG10(ESBASW)
          ESH2O= 10.0**(AA+B+C+D+E)
c         convert to millibars
c         ESH2O= 0.001*ESH2O
          if(l_Pascals) then
c           convert to Pascals
            ESH2O= 0.1*ESH2O
          endif
c
c         DERIVATIVE d(ln es)/ dT  = qcon * Lv(T)/ (Rd*T*T)
c
c         DERIVATIVE OVER ICE
          AA = 5721.9
          B = 3.56654*TEM
          C = -0.0073908*TEM*TEM
          VAICE=(AA+B+C)/(TEM*TEM)
c         latent heat over ice
          XLICE=VAICE*(RD*TEM*TEM)/0.62197
c
c         DERIVATIVE OVER WATER
          AA = 6790.5
          B = -5.02808*TEM
          E = -0.0304*TEM
          C = 4916.826*(10.0**E)*TEM*TEM
          E = -1302.88/TEM
          D = 174207.704*(10.0**E)
          VAH2O=(AA+B+C+D)/(TEM*TEM)
c         latent heat over water
          XLH2O=VAH2O*(RD*TEM*TEM)/0.62197
c
c           print table of values
c
          if(lprint) then
            temc = TEM-TBASI
            if(mod(i,50).eq.31)
     &       write(6,7) i, TEM, temc, ESH2O, ESICE, XLH2O
     &                   , XLICE
          endif
c 
c           save vapor pressures in table
c
c         vapor pressure over liquid water
          vpsat(i)=ESH2O
c         below tcutoff use ESICE
          if(lcutoff.and.TEM.le.tcutoff)vpsat(i)=ESICE
c         vapor pressure over ice
          visat(i)=ESICE
c         above 0C use ESH2O
          if(TEM.ge.273.16)visat(i)=ESH2O
c
c           save variation of vapor pressures
c
c         variation of vapor pressure over liquid water
          varvp(i)=VAH2O
c         below tcutoff use VAICE
          if(lcutoff.and.TEM.le.tcutoff)varvp(i)=VAICE
c         variation of vapor pressure over ice
          varvi(i)=VAICE
c         above 0C use VAH2O
          if(TEM.ge.273.16)varvi(i)=VAH2O
c
c
 100    continue
c
        if(ltest) then
          write(6,10)
          write(6,12)
          do 105 i=1,ntot
            TEM = 120.16+0.1*FLOAT(i-1)
            temc = TEM-TBASI
            if(i.eq.1) then 
               varest = 10.0*alog(vpsat(2)/vpsat(1))
            elseif (i.eq.ntot) then
               varest = 10.0*alog(vpsat(ntot)/vpsat(ntot-1))
            else
               varest = 5.0*alog(vpsat(i+1)/vpsat(i-1))
            endif
            if(mod(i,50).eq.31)
     &       write(6,7) i, TEM, temc, vpsat(i), vpsat(i+1),
     &         vpsat(i-1),varest, varvp(i)
 105      continue
        endif   
c
c***********************************************************************
c       end of lookup table
c***********************************************************************
c
        lfirst=.false.
      endif
c
c
      if(l_interp) then
c
c***********************************************************************
c       linear interpolation using lookup table
c***********************************************************************
c
        do i = 1,len
c
          t1=1.0+10.0*(tp(i)-120.16)
          t1=min(ntot-1.0,t1)
          t1=max(t1,1.0001)
          ic= int(t1)
c
          if(l_ice) then
c
c         saturation vapor pressure over ice
          vp(i)= visat(ic)+(visat(1+ic)-visat(ic))*(t1-float(ic))
c
c         variation of natural log of saturation vapor pressure
c          over ice
c          vd(i)= varvi(ic)+(varvi(1+ic)-varvi(ic))*(t1-float(ic))
c
          else
c
c         saturation vapor pressure over water
          vp(i)= vpsat(ic)+(vpsat(1+ic)-vpsat(ic))*(t1-float(ic))
c
c         variation of natural log of saturation vapor pressure
c          over water
c          vd(i)= varvp(ic)+(varvp(1+ic)-varvp(ic))*(t1-float(ic))
c
c
          endif
c
c         convert variation to latent heat
c
c          if(l_xl) then
c          vd(i)=vd(i)*(rdf*tp(i)*tp(i))
c          endif
c
        enddo
c
      else
c
c***********************************************************************
c         exact calculation
c***********************************************************************
c
        do i = 1, len
c
          TEM =  tp(i)
c
c         SATURATION VAPOR PRESSURE OVER ICE
          AA  = -9.09718 *(TBASI/TEM-1.0)
          B   = -3.56654 *ALOG10(TBASI/TEM)
          C   =  0.876793*(1.0-TEM/TBASI)
          E   =  ALOG10(ESBASI)
          ESICE= 10.0**(AA+B+C+E)
c         convert to millibars
c         ESICE= 0.001*ESICE
          if(l_Pascals) then
c           convert to Pascals
            ESICE= 0.1*ESICE
          endif
c
c         SATURATION VAPOR PRESSURE OVER WATER
          AA  = -7.90298*(TBASW/TEM-1.0)
          B   =  5.02808*ALOG10(TBASW/TEM)
          C   = -1.3816E-07*
     &           (10.0**((1.0-TEM/TBASW)*11.344)-1.0)
          D   =  8.1328E-03*
     &           (10.0**((TBASW/TEM-1.0)*(-3.49149))-1.0)
          E   =  ALOG10(ESBASW)
          ESH2O= 10.0**(AA+B+C+D+E)
c         convert to millibars
c         ESH2O= 0.001*ESH2O
          if(l_Pascals) then
c           convert to Pascals
            ESH2O= 0.1*ESH2O
          endif
c
c         DERIVATIVE d(ln es)/ dT  = qcon * Lv(T)/ (Rd*T*T)
c
c         DERIVATIVE OVER ICE
          AA = 5721.9
          B = 3.56654*TEM
          C = -0.0073908*TEM*TEM
          VAICE=(AA+B+C)/(TEM*TEM)
          XLICE=VAICE*(RD*TEM*TEM)/0.62197
c
c         DERIVATIVE OVER WATER
          AA = 6790.5
          B = -5.02808*TEM
          E = -0.0304*TEM
          C = 4916.826*(10.0**E)*TEM*TEM
          E = -1302.88/TEM
          D = 174207.704*(10.0**E)
          VAH2O=(AA+B+C+D)/(TEM*TEM)
          XLH2O=VAH2O*(RD*TEM*TEM)/0.62197
c
          if(l_ice) then
            vp(i)=ESICE
c            vd(i)=VAICE
c            if(l_xl) vd(i)=XLICE
c           above 0C use ESH2O
            if(TEM.ge.273.16) then
              vp(i)=ESH2O
c              vd(i)=VAH2O
c              if(l_xl) vd(i)=XLH2O
            endif
          else
            vp(i)=ESH2O
c            vd(i)=VAH2O
c            if(l_xl) vd(i)=XLH2O
c           below tcutoff use ESICE
            if(lcutoff.and.TEM.le.tcutoff) then
              vp(i)=ESICE
c              vd(i)=VAICE
c              if(l_xl) vd(i)=XLICE
            endif
          endif
c
        enddo
c
      endif
c
c***********************************************************************
c
c***********************************************************************
c
c
c***********************************************************************
c     format statements
c***********************************************************************
c
    7 format(1x,i5,2(5x,f8.2),4x,4e16.8)
    8 format(1x,i5,2(5x,f8.2),4x,5e16.8)
   10 format(1x,' LEVEL           TEMPERATURE      ',
     :          ' VAPOR PRESSURE (MB)     &          LATENT HEAT',
     :          ' after Goff-Gratch (1946)')
   11 format(1X,'         (DEG C)      (DEG K)  ',
     :          '          Liquid        Ice      Liquid     Ice  ')
   12 format(1X,'         (DEG C)      (DEG K)  ',
     :          '          Liquid        Ice      Blend    Derivative')
      return
      end
