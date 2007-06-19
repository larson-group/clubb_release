!----------------------------------------------------------------------
      subroutine coamps_micro_driver
     .       ( kk, kmax, timea, deltf,
     .         rtm, wmm, p, exner, rhot,
     .         thlm, ricem, rrm, rgraupelm, rsnowm,
     .         nne, nrdamp, j, icon, icond,
     .         lice,
     .         rcm, Ncm, Nrm, Ncnm, Nim,
     .         cond, Vsnow, Vice, Vrr, VNr, Vgraupel,
     .         ritend, rrtend, rgtend, 
     .         rsnowtend, nrmtend,
     .         rttend, thlmtend )

!    Description:
!      Subroutine to compute ice, as it is done in COAMPS, using
!      adjtq.F.
!
!    References:
!      Rutledge and Hobbs, 1984; COAMPS Users Guide.
!----------------------------------------------------------------------
                              
      use constants
      use grid_class, only: zt2zm
#ifdef STATS
      use hoc_stats, only: zt, zm, lstats_samp,
     .    imean_vol_rad_rain,
     .    imean_vol_rad_cloud
#endif

      implicit none

      ! External Calls
      external :: 
     .  gamma, ! From COAMPS, and not the same gamma approx. used in HOC
     .  adjtq  ! COAMPS microphysics subroutine

      real, external ::
     .  rsat,
     .  rsati

      ! Local Constants
      real, parameter ::
     .  aa0 = -0.267,  ! All of these are constants set in COAMPS and used by adjtq.
     .  aa1 = 5150.,
     .  aa2 = -1.0225e6,
     .  aa3 = 7.55e7,
     .  abar = 19.3,
     .  apr = 3000.,
     .  aprpr = 1.139,
     .  bsnow = 0.11,
     .  cbeta = 0.6,
     .  cnzero = 0.01,
     .  cimass = 9.4e-10,
     .  cw = 4218.,
     .  difvap = 2.26e-5,
     .  erc = 1.,
     .  esi = 0.1,
     .  eri = 1.,
     .  egc = 1.,
     .  esc = 1.,
     .  esr = 0.4,
     .  egi = 0.1,
     .  egr = 1.0,
     .  egs = 0.1,
     .  mw = 18.016,
     .  praut1 = 0.001,
     .  praut2 = 0.0004,
     .  rholiq = 1000.,
     .  rhosno = 100.,
     .  rhogrp = 400.,
     .  rnzero = 8.0e6,
     .  snzero = 2.0e7,
     .  gnzero = 4.0e6,
     .  therco = 2.43e-2,
     .  tice = 269.16,
     .  tvr1 = -0.267,
     .  tvr2 = 206.,
     .  tvr3 = -2045.,
     .  tvr4 = 9060.,
     .  tzero = 273.16,
     .  visair = 1.718e-5,
     .  bgrp = 0.37,
     .  ex1 = 0.2,
     .  pcut = 1.0e-10 ! Lower threshold for calculation in COAMPS

      integer, parameter ::
     .  n1d     = 1, ! 1d graphics parameters
     .  i1dflg  = 0, ! 1d graphics parameters
     .  maxpt1d = 0, ! 1d graphics parameters
     .  maxvr1d = 0, ! 1d graphics parameters
     .  ipts    = 1  ! ??

      real, dimension(n1d), parameter ::
     .  i1d = (/0./), ! 1d graphics parameters
     .  j1d = (/0./)  ! 1d graphics parameters

      real, dimension(1,1), parameter ::
     .  xland = 1.0, ! Land/Sea assumption
     .  wtm   = 1.0  ! Weighting array for mass point (never used)

      ! Input Variables
      integer, intent(in) ::
     .  kk,    ! Number of COAMPS m gridpoints in the vertical (gr%nnzp-1)
     .  kmax   ! Maximum array size (kk + ??)

      real, intent(in) ::
     .  timea,        ! Output time             [s]
     .  deltf         ! Timestep, I think       [s]

      integer, intent(in) ::
     .  nne,    ! Has to do with horiz domain; 1 for HOC?
     .  nrdamp, ! Number of upper model levels in sponge layer (8?)
     .  j,      ! j-coordinate (always 1 for HOC)
     .  icon,   ! Ice nucleation type
     .  icond   ! Autoconversion; 1=Kessler, 2=Manton/Cotton, 3=K&K, 4=none

      logical, intent(in) ::
     .  lice           ! Logical: is ice on?

      real, dimension(kk+1), intent(in) ::
     .  rtm,  ! Total water mixing ratio                        [kg/kg]
     .  rcm,  ! Cloud water mixing ratio                        [kg/kg]
     .  wmm,  ! Vertical wind                                   [m/s]
     .  p,    ! Pressure                                        [Pa]
     .  exner,! Mean exner function                             [-]
     .  rhot, ! Mean density                                    [kg/m^3]
     .  thlm  ! Liquid potential temperature                    [K]

      real, dimension(kk+1), intent(in) ::
     .  ricem,     ! Ice water mixing ratio     [kg/kg]
     .  rrm,       ! Rain water mixing ratio    [kg/kg]
     .  rgraupelm, ! Graupel water mixing ratio [kg/kg]
     .  rsnowm,    ! Snow water mixing ratio    [kg/kg]
     .  Nrm        ! Number of rain drops       [count/m^3]

      real, dimension(kk+1), intent(inout) ::
     .  Ncm,       ! Number of cloud droplets   [count/m^3]
     .  Ncnm,      ! Number of cloud nuclei     [count/m^3]
     .  Nim        ! Number of ice crystals     [count/m^3]

      ! Output Variables
      real, dimension(kk+1), intent(out) ::
     .  ritend,    ! d(ri)/dt                   [kg/kg/s]
     .  rrtend,    ! d(rr)/dt                   [kg/kg/s]
     .  rgtend,    ! d(rg)/dt                   [kg/kg/s]
     .  rsnowtend, ! d(rsnow)/dt                [kg/kg/s]
     .  rttend,    ! d(rt)/dt                   [kg/kg/s]
     .  thlmtend,  ! d(thlm)/dt                 [K/s]
     .  nrmtend    ! d(Nrm)/dt                  [count/m^3/s]

      real, dimension(kk+1), intent(out) ::
     .  Vrr,     ! Rain mixing ratio fall speed   [m/s]
     .  VNr,     ! Rain conc. fall speed          [m/s]
     .  Vsnow,   ! Snow fall speed                [m/s]
     .  Vgraupel,! Graupel fall speed             [m/s]
     .  Vice     ! Pristine ice fall speed        [m/s]

      ! Local Variables

      ! Variables on the w grid 
      real, dimension(1,1,kk+1) ::
     .  w3,      ! Vertical wind on the w grid          [m/s]
     .  pr3d,    ! Pressure on w grid                   [Pa]
     .  qsatv3d, ! Saturation mr array?                 [kg/kg]
     .  temp3d,  ! Temperature on w grid                [K]
     .  qsati3d, ! Saturation mr over ice array?        [kg/kg]
     .  th2t3d   ! Mean exner function on w grid        [-]

      ! Variables on the m grid
      real, dimension(1,1,kk) ::
     .  qi3,   ! Pristine ice mixing ratio               [kg/kg]
     .  qr3,   ! Rain water mixing ratio                 [kg/kg]
     .  qg3,   ! Graupel mixing ratio                    [kg/kg]
     .  qs3    ! Snow water mixing ratio                 [kg/kg]

      real ::
     .  gmbov2,
     .  gmbov2g,
     .  gmbp3,
     .  gm3,
     .  gm4,
     .  gm5,
     .  gm6,
     .  gm7,
     .  gm8,
     .  gm9,
     .  ex2,
     .  ex2g,
     .  ex3,
     .  ex7,
     .  ex7g,
     .  ex4,
     .  ex4g,
     .  ex5,
     .  Lf,
     .  hlvoka,
     .  hlsoka,
     .  hlvorv,
     .  hlsorv,
     .  rvochi,
     .  cpor,
     .  lfocp,
     .  lvocp,
     .  lsocp,
     .  hkaolf,
     .  sloper,
     .  slopes,
     .  slopeg,
     .  eic

      real, dimension(kk+1) ::
     .  thm,
     .  rvm

      real, dimension(1,1,kk) ::
     .  qt3,        ! Total water mixing ratio
     .  qv3,        ! Water vapor mixing ratio
     .  qc3,        ! Cloud water mixing ratio
     .  th3,        ! potential temperature
     .  p3,         ! perturbation exner function
     .  nc3,        ! Number of cloud droplets
     .  nr3,        ! Number of rain drops
     .  ncn3,       ! Number of cloud nuclei
     .  ni3,        ! Number of ice crystals
     .  exbm,       ! Mean exner function
     .  rbm,        ! Mean density
     .  snslp,      ! These variables are the individual microphysical terms.  Michael Falk
     .  cond,       ! condensation/evaporation of liquid water
     .  pcond,      ! condensation/evaporation of cloud water
     .  psmlti,     ! melting of cloud ice
     .  psacw,      ! collection of cloud water by snow
     .  pgacw,      ! ???
     .  piacw,      ! ???
     .  pchomo,     ! ???
     .  praut,      ! autoconversion of cloud water
     .  pracw,      ! collection of cloud water by rain water
     .  pdepi,      ! depositional growth of cloud ice
     .  pint,       ! initiation of cloud ice
     .  pgdep,      ! ???
     .  pconv,      ! conversion of cloud ice to snow
     .  psaci,      ! collection of cloud ice by snow
     .  pgaci,      ! ???
     .  praci,      ! ???
     .  prevp,      ! evaporation of raindrops
     .  psdep,      ! depositional growth of snow
     .  pmltge,     ! ???
     .  pgmlt,      ! ???
     .  psmlt,      ! melting of snow
     .  pgacrm,     ! ???
     .  pgacwm,     ! ???
     .  pracs,      ! ???
     .  pgshr,      ! ???
     .  pgacr,      ! ???
     .  psacr,      ! ???
     .  piacr,      ! ???
     .  prhomo,     ! ???
     .  pgacs,      ! ???
     .  pmltse,     ! ???
     .  pwacs       ! ??? end Michael Falk variables

      real, dimension(1,1,kk+1) ::
     .  fallr,  ! Fall speed for rain mixing ratio              [m/s]
     .  falln,  ! Fall speed for rain drop number conc.         [m/s]
     .  falli,  ! Fall speed for pristine ice mixing ratio      [m/s]
     .  snowv,  ! Fall speed for snow mixing ratio              [m/s]
     .  fallg   ! Fall speed for graupel mixing ratio           [m/s]

      real, dimension(1,1,kk) ::
     .  rvc,     ! Cloud droplet radius         [cm]
     .  rvr      ! Rain droplet radius          [cm]

      real, dimension(1,kk,1) ::
     .  ary1d        ! 1d graphics parameters

      integer ::
     .  i,
     .  k, kp1, km1, ! Loop control variables
     .  len,         ! # of saturated points???
     .  icase        ! Which case?

      integer, dimension(1,1) ::
     .  nkpts        

      integer, dimension(kk) ::
     .  icomp          !

      integer, dimension(1,kk) ::
     .  kcomp          !

      logical ::
     .  ldrizzle,      ! is drizzle on?
     .  lgrpl          ! is graupel on?


!----------------------------------------------------------------------

      ! Begin coamps_micro_driver code
      ldrizzle = .false.
      lgrpl    = .true.
      len      = kk
      icase    = 62
      eic      = 1.0

      ! Set up initial fields

      ! Compute quantities for computing tendencies
      rvm = rtm - rcm

      thm(1:kk+1) = thlm(1:kk+1)
     .            + ( Lv /( Cp * exner(1:kk+1) )* rcm(1:kk+1) )

      ! Setup COAMPS w grid variables
      w3(1,1,1:kk+1) = wmm(1:kk+1)

!     do k=1, kk+1, 1
!       pr3d(1,1,k)    = zt2zm( p, k )
!       th2t3d(1,1,k)  = zt2zm( exner, k )
!       temp3d(1,1,k)  = zt2zm( thm, k ) * th2t3d(1,1,k)
!       qsatv3d(1,1,k) = rsat( pr3d(1,1,k), temp3d(1,1,k) )
!       qsati3d(1,1,k) = rsati( pr3d(1,1,k), temp3d(1,1,k) )
!     end do

      ! Setup COAMPS m grid variables that for some inexplicable reason
      ! are put in variable arrays with dimension kk+1
      pr3d(1,1,1:kk)   = p(2:kk+1)
      th2t3d(1,1,1:kk) = exner(2:kk+1)
      temp3d(1,1,1:kk) = thm(2:kk+1) * exner(2:kk+1)

      do k=1, kk, 1
        qsatv3d(1,1,k) = rsat( pr3d(1,1,k), temp3d(1,1,k) )
        qsati3d(1,1,k) = rsati( pr3d(1,1,k), temp3d(1,1,k) )
      end do

      ! Setup COAMPS m grid variables
      qt3(1,1,1:kk)  = rtm(2:kk+1)
      qc3(1,1,1:kk)  = rcm(2:kk+1)
      qr3(1,1,1:kk)  = rrm(2:kk+1)
      qg3(1,1,1:kk)  = rgraupelm(2:kk+1)
      qs3(1,1,1:kk)  = rsnowm(2:kk+1)
      qi3(1,1,1:kk)  = ricem(2:kk+1)
      exbm(1,1,1:kk) = exner(2:kk+1)
      rbm(1,1,1:kk)  = rhot(2:kk+1)
      th3(1,1,1:kk)  = thm(2:kk+1)
      qv3(1,1,1:kk)  = rvm(2:kk+1)

      do k=1,kk
        p3(1,1,k)   = 0.0

        ! Convert from MKS units as needed
        nc3(1,1,k)  = Ncm(k+1) * 1.e-6
        nr3(1,1,k)  = Nrm(k+1) * 1.e-6
        ncn3(1,1,k) = Ncnm(k+1) * 1.e-6

        ni3(1,1,k)  = Nim(k+1)

      end do


      call gamma(3.0,gm3)
      call gamma(4.0,gm4)
      call gamma(5.0,gm5)
      call gamma(6.0,gm6)
      call gamma(7.0,gm7)
      call gamma(8.0,gm8)
      call gamma(9.0,gm9)
      call gamma(bsnow+3.0,gmbp3)
      call gamma(bsnow*0.5 + 2.5,gmbov2)
      call gamma(bgrp*0.5 + 2.5,gmbov2g)

      ex2  = bsnow * 0.5 + 2.5
      ex2g = bgrp * 0.5 + 2.5
      ex3  = bsnow + 3.0
      ex7  = 0.44 * gmbov2
      ex7g = 0.31 * gmbov2g
      ex4  = aprpr/visair
      ex4g = abar/visair
      ex5  = pi*aprpr*snzero*gmbp3/4.0

      Lf     = Ls - Lv   ! Latent heat of fusion
      hlvoka = Lv/therco
      hkaolf = therco/Lf
      hlsoka = Ls/therco
      hlvorv = Lv/Rv
      hlsorv = Ls/Rv
      rvochi = Rv/difvap
      cpor   = Cp / Rd
      lfocp  = Lf/Cp
      lvocp  = Lv/Cp
      lsocp  = Ls/Cp

      ary1d(1,1:kk,1) = 0. ! 1d graphics parameters
      nkpts = kk
      do i=1,ipts
        icomp(i) = 1
        do k=1,kk
          kcomp(1,k) = k
        end do
      end do

      sloper = pi * rholiq * rnzero * 1.0e-8
      slopes = pi * rhosno * snzero * 1.0e-8
      slopeg = pi * rhogrp * gnzero * 1.0e-8

      ! Call the actual COAMPS microphysics scheme
      call adjtq
     1       (cond,p3,qc3,qi3,qr3,qg3,qs3,qv3,th3,w3,pr3d
     2       ,qsatv3d,temp3d,qsati3d,th2t3d,wtm,exbm,rbm
     3       ,nc3,nr3,ncn3,ni3,cp,deltf,Lf,Ls,Lv 
     4       ,pcut,p0,Rd,Rv,sloper,slopes,slopeg,timea,lice 
     5       ,nne,kk,i1d,j1d,ary1d,i1dflg,n1d,maxpt1d,maxvr1d 
     6       ,kmax,nrdamp,ipts,nkpts,icomp,kcomp,j
     7       ,xland,aa0,aa1,aa2,aa3,abar,apr,aprpr,bsnow
     8       ,cbeta,cnzero,cimass,cpor,cw,difvap,erc,esi,eic 
     9       ,eri,egc,esc,esr,egi,egr,egs,mw,pi,praut1,praut2 
     1       ,rholiq,rhosno,rnzero,snzero,gnzero,therco,tice 
     2       ,tvr1,tvr2,tvr3,tvr4,tzero,visair,gm3,gm4,gm5,gm6 
     3       ,gm7,gm8,gm9,gmbp3,gmbov2,gmbov2g,bgrp,ex1,ex2
     4       ,ex2g,ex3,hlvoka,hkaolf,hlsoka,hlvorv,hlsorv 
     5       ,rvochi,lfocp,lvocp,lsocp,ex7,ex7g,ex4,ex4g,ex5 
     6       ,ldrizzle,lgrpl,icon,icond,len,icase,snowv,snslp
     7       ,pcond,psmlti,psacw,pgacw,piacw,pchomo,praut,pracw
     8       ,pdepi,pint,pgdep,pconv,psaci,pgaci,praci,prevp,psdep
     9       ,pmltge,pgmlt,psmlt,pgacrm,pgacwm,pracs,pgshr 
     1       ,pgacr,psacr,piacr,prhomo,pgacs,pmltse,pwacs
     2       ,falli(:,:,2:kk+1),fallg(:,:,2:kk+1)
     3       ,fallr(:,:,2:kk+1),falln(:,:,2:kk+1)
     4       ,rvc,rvr)

      ! Assure positive definiteness in nc3/nr3/ncn3 fields
      do k=1, kk, 1

        if (nr3(1,1,k) < 0.) then
          ncn3(1,1,k) = ncn3(1,1,k) + nr3(1,1,k)
          nr3(1,1,k)  = 0.
        end if

        if (nc3(1,1,k) < 0.) then
          ncn3(1,1,k) = ncn3(1,1,k) + nc3(1,1,k)
          nc3(1,1,k)  = 0.
        end if

        if (ncn3(1,1,k) < 0.) then
          ncn3(1,1,k)  = 0.
        end if
      end do ! k=1..kk


      ! Transfer back to HOC arrays
      do k=1, kk, 1
        ! Convert to MKS as needed
        Ncm(k+1)  = nc3(1,1,k) * 1.e6
        Ncnm(k+1) = ncn3(1,1,k) * 1.e6
        Nim(k+1)  = ni3(1,1,k)
      end do ! k=1..kk+1

      ! Linear extrapolation for the ghost point of fall speeds
      fallr(1,1,1) = .5 * ( fallr(1,1,2) + fallr(1,1,3) )
      falln(1,1,1) = .5 * ( falln(1,1,2) + falln(1,1,3) )
      snowv(1,1,1) = .5 * ( snowv(1,1,2) + snowv(1,1,3) )
      falli(1,1,1) = .5 * ( falli(1,1,2) + falli(1,1,3) )
      fallg(1,1,1) = .5 * ( fallg(1,1,2) + fallg(1,1,3) )

      Vrr      = zt2zm( fallr(1,1,:) )
      VNr      = zt2zm( falln(1,1,:) )
      Vsnow    = zt2zm( snowv(1,1,:) )
      Vice     = zt2zm( falli(1,1,:) )
      Vgraupel = zt2zm( fallg(1,1,:) )


      ! Compute tendencies
      do k=1, kk, 1
        rrtend(k+1)    = ( qr3(1,1,k) - rrm(k+1) ) / deltf
        rgtend(k+1)    = ( qg3(1,1,k) - rgraupelm(k+1) ) / deltf
        ritend(k+1)    = ( qi3(1,1,k) - ricem(k+1) ) / deltf
        nrmtend(k+1)   = ( nr3(1,1,k)*1.e6 - Nrm(k+1) ) / deltf ! Conversion factor
        rsnowtend(k+1) = ( qs3(1,1,k) - rsnowm(k+1) ) / deltf
        rttend(k+1)    = ((qv3(1,1,k) - rvm(k+1)) / deltf)
     .                 + ((qc3(1,1,k) - rcm(k+1)) / deltf)
        thlmtend(k+1) 
     .  = ( ( th3(1,1,k) - (Lv / (Cp * exbm(1,1,k)) * qc3(1,1,k) ) )
     .      - thlm(k+1) ) / deltf
      end do ! k=1..kk

      rrtend(1)    = 0.0
      rgtend(1)    = 0.0
      rgtend(1)    = 0.0
      nrmtend(1)   = 0.0
      rsnowtend(1) = 0.0
      rttend(1)    = 0.0
      thlmtend(1)  = 0.0

#ifdef STATS
      if ( lstats_samp ) then
        if ( imean_vol_rad_rain > 0 ) then
          zt%x(2:,imean_vol_rad_rain) 
     .    = zt%x(2:,imean_vol_rad_rain) + rvr(1,1,:) / 100.0
          zt%n(2:,imean_vol_rad_rain) 
     .    = zt%n(2:,imean_vol_rad_rain) + 1
        end if

        if ( imean_vol_rad_cloud > 0 ) then
          zt%x(2:,imean_vol_rad_cloud) 
     .    = zt%x(2:,imean_vol_rad_cloud) + rvc(1,1,:) / 100.0
          zt%n(2:,imean_vol_rad_cloud) 
     .    = zt%n(2:,imean_vol_rad_cloud) + 1
        end if
      end if
#endif /*STATS*/


      return
      end subroutine coamps_micro_driver
