!----------------------------------------------------------------------
      subroutine coamps_micro_driver
     .       (grnnzp,wtm,timea,time_initial,deltf,
     .        grzt,qt3,w3,w3dw1,exbm,rbm,
     .        thlm,qi3,qr3,qg3,qs3,
     .        nne,kk,kmax,nrdamp,j,icon,icond,
     .        lice,
     .        qc3,nc31d,nr31d,ncn31d,ni31d,
     .        cond,snowv,falli,fallr,fallg,
     .        ritend,rrtend,rgtend,rsnowtend,
     .        rttend,thlmtend)

!    Description:
!      Subroutine to compute ice, as it is done in COAMPS, using
!      adjtq.F.
!
!    References:
!      Rutledge and Hobbs, 1984; COAMPS Users Guide.
!----------------------------------------------------------------------
                              
      use constants
      implicit none

      ! External Calls
      real, external ::
     .  gamma,
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
     .  eic = 1.,
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
     .  pcut = 1.0e-10,
     .  i1d = 0, ! 1d graphics parameters
     .  j1d = 0, ! 1d graphics parameters
     .  n1d = 0  ! 1d graphics parameters

      integer, parameter ::
     .  i1dflg = 0,  ! 1d graphics parameters
     .  maxpt1d = 0, ! 1d graphics parameters
     .  maxvr1d = 0, ! 1d graphics parameters
     .  ipts = 1,    ! ??
     .  xland = 1    ! Land/Sea assumption


      ! Input Variables
      integer, intent(in) ::
     .  grnnzp        ! number of gridpoints in the vertical (this needs to be first)

      real, intent(in) ::
     .  wtm,          ! weighting array for mass point (never used)
     .  timea,        ! Output time
     .  time_initial, ! Initial time
     .  deltf         ! Timestep, I think

      real, dimension(grnnzp), intent(in) ::
     .  grzt,         ! List of heights associated with gridpoints
     .  qt3(grnnzp)   ! cloud water mixing ratio

      real, dimension(1,1,grnnzp), intent(in) ::
     .  w3,           ! vertical wind
     .  w3dw1,        ! aka pr3d    -- pressure array?
     .  exbm,         ! mean exner function
     .  rbm,          ! mean density
     .  thlm,         ! liquid water potential temperature
     .  qi3,          ! pristine ice mixing ratio
     .  qr3,          ! rain water mixing ratio
     .  qg3,          ! graupel mixing ratio
     .  qs3           ! snow water mixing ratio

      integer, intent(in) ::
     .  nne,           ! has to do with horiz domain; 1 for HOC?
     .  kk,            ! equal to grnnzp
     .  kmax,          ! maximum array size (grnnzp + ??)
     .  nrdamp,        ! Number of upper model levels in sponge layer (8?)
     .  j,             ! j-coordinate (always 1 for HOC)
     .  icon,          ! ice nucleation type
     .  icond          ! autoconversion; 1=Kessler, 2=Manton/Cotton, 3=K&K, 4=none

      logical, intent(in) ::
     .  lice           ! Logical: is ice on?


      ! Input-Output Variables
      real, dimension(1,1,grnnzp), intent(inout) ::
     .  qc3        ! cloud water mixing ratio

      real, dimension(grnnzp), intent(inout) ::
     .  nc31d,     ! Number of cloud droplets
     .  nr31d,     ! Number of rain drops
     .  ncn31d,    ! Number of cloud nuclei
     .  ni31d      ! Number of ice crystals


      ! Output Variables
      real, dimension(1,1,grnnzp), intent(out) ::
     .  cond,      ! condensation/evaporation of liquid water
     .  snowv,     ! snow fall speed
     .  falli,     ! pristine ice fall speed
     .  fallr,     ! rain drop fall speed
     .  fallg,     ! graupel fall speed
     .  ritend,    ! d(ri)/dt
     .  rrtend,    ! d(rr)/dt
     .  rgtend,    ! d(rg)/dt
     .  rsnowtend, ! d(rsnow)/dt
     .  rttend,    ! d(rt)/dt
     .  thlmtend   ! d(thlm)/dt


      ! Local Variables

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
     .  slopeg

      real, dimension(1,1,grnnzp) ::
     .  w3dw2,    ! aka qsatv3d -- saturation mr array?
     .  w3dw3,    ! aka temp3d  -- temperature array?
     .  w3dw4,    ! aka qsati3d -- saturation mr over ice array?
     .  w3dw5,    ! aka th2t3d  -- exner function array?
     .  qv3,      ! water vapor mixing ratio
     .  th3,      ! potential temperature
     .  p3,       ! perturbation exner function
     .  nc3,      ! Number of cloud droplets
     .  nr3,      ! Number of rain drops
     .  ncn3,     ! Number of cloud nuclei
     .  ni3,      ! Number of ice crystals
     .  rcbefore,    ! Temporary variable for calculating rc(qc3) tendency
     .  rsnowbefore, ! Temporary variable for calculating rsnow(qs3) tendency
     .  ribefore,    ! Temporary variable for calculating ri(qi3) tendency
     .  rvbefore,    ! Temporary variable for calculating rv(qv3) tendency
     .  rrbefore,    ! Temporary variable for calculating rr(qr3) tendency
     .  rgbefore,    ! Temporary variable for calculating rg(qg3) tendency
     .  rtbefore,    ! Temporary variable for calculating rt(qt3) tendency
     .  qt3_fortend, ! Temporary qt3 tendency
     .  thlmbefore,  ! Temporary variable for calculating thlm tendency
     .  thmbefore,   ! Temporary variable for calculating thm tendency
     .  snslp,       ! These variables are the individual microphysical terms.  Michael Falk
     .  pcond,       ! condensation/evaporation of cloud water
     .  psmlti,      ! melting of cloud ice
     .  psacw,       ! collection of cloud water by snow
     .  pgacw,       ! ???
     .  piacw,       ! ???
     .  pchomo,      ! ???
     .  praut,       ! autoconversion of cloud water
     .  pracw,       ! collection of cloud water by rain water
     .  pdepi,       ! depositional growth of cloud ice
     .  pint,        ! initiation of cloud ice
     .  pgdep,       ! ???
     .  pconv,       ! conversion of cloud ice to snow
     .  psaci,       ! collection of cloud ice by snow
     .  pgaci,       ! ???
     .  praci,       ! ???
     .  prevp,       ! evaporation of raindrops
     .  psdep,       ! depositional growth of snow
     .  pmltge,      ! ???
     .  pgmlt,       ! ???
     .  psmlt,       ! melting of snow
     .  pgacrm,      ! ???
     .  pgacwm,      ! ???
     .  pracs,       ! ???
     .  pgshr,       ! ???
     .  pgacr,       ! ???
     .  psacr,       ! ???
     .  piacr,       ! ???
     .  prhomo,      ! ???
     .  pgacs,       ! ???
     .  pmltse,      ! ???
     .  pwacs        ! ??? end Michael Falk variables

      real, dimension(1,grnnzp,1) ::
     .  ary1d        ! 1d graphics parameters

      integer ::
     .  i,
     .  k,             ! loop control variables
     .  len,           ! # of saturated points???
     .  icase,         ! Which case?
     .  kpts        

      integer, dimension(grnnzp) ::
     .  icomp          !

      integer, dimension(1,grnnzp) ::
     .  kcomp          !

      logical ::
     .  ldrizzle,      ! is drizzle on?
     .  lgrpl          ! is graupel on?


!----------------------------------------------------------------------

      ! Begin coamps_micro_driver code
      ldrizzle = .TRUE.
      lgrpl    = .TRUE.
      len      = grnnzp
      icase    = 62

      ! Set up initial fields
      do k=1,grnnzp
        p3(1,1,k) = 0.
        nc3(1,1,k) = nc31d(k)
        nr3(1,1,k) = nr31d(k)
        ncn3(1,1,k) = ncn31d(k)
        ni3(1,1,k) = ni31d(k)

        qv3(1,1,k) = qt3(k) - qc3(1,1,k) - qi3(1,1,k) - qr3(1,1,k)
     .                  - qg3(1,1,k) - qs3(1,1,k)
        qc3(1,1,k) = qc3(1,1,k) - qi3(1,1,k) - qr3(1,1,k) - qg3(1,1,k)
     .              - qs3(1,1,k)

        if (qc3(1,1,k) < 0.) qc3(1,1,k) = 0.
        th3(1,1,k) = thlm(1,1,k) +
     .               (Lv / (Cp * exbm(1,1,k)) * qc3(1,1,k))
      end do

      if (timea == time_initial) then
        do k=1,grnnzp
          ncn3(1,1,k) = 30. * (1. + exp(-grzt(k)/2000))
        end do
      end if

      gm3 = gamma(3.0,gm3)
      gm4 = gamma(4.0,gm4)
      gm5 = gamma(5.0,gm5)
      gm6 = gamma(6.0,gm6)
      gm7 = gamma(7.0,gm7)
      gm8 = gamma(8.0,gm8)
      gm9 = gamma(9.0,gm9)
      gmbp3 = gamma(bsnow+3.0,gmbp3)
      gmbov2 = gamma(bsnow*0.5 + 2.5,gmbov2)
      gmbov2g = gamma(bgrp*0.5 + 2.5,gmbov2g)

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

      ary1d(1,1:grnnzp,1) = 0. ! 1d graphics parameters
      kpts = grnnzp
      do i=1,ipts
        icomp(i) = 1
        do k=1,kpts
          kcomp(1,k) = k
        end do
      end do

      sloper = pi * rholiq * rnzero * 1.0e-8
      slopes = pi * rhosno * snzero * 1.0e-8
      slopeg = pi * rhogrp * gnzero * 1.0e-8


      ! Saving pre-adjtq values for computation of tendencies
      do k=1,grnnzp
        rcbefore(1,1,k) = qc3(1,1,k)
        rsnowbefore(1,1,k) = qs3(1,1,k)
        ribefore(1,1,k) = qi3(1,1,k)
        rvbefore(1,1,k) = qv3(1,1,k)
        rrbefore(1,1,k) = qr3(1,1,k)
        rgbefore(1,1,k) = qg3(1,1,k)
        rtbefore(1,1,k) = qc3(1,1,k) + qs3(1,1,k) + qi3(1,1,k) +
     .                    qv3(1,1,k) + qr3(1,1,k) + qg3(1,1,k)
        thlmbefore(1,1,k) = thlm(1,1,k)
        thmbefore(1,1,k) = th3(1,1,k)

        w3dw3(1,1,k) = th3(1,1,k) * exbm(1,1,k)
        w3dw2(1,1,k) = rsat(w3dw1(1,1,k), w3dw3(1,1,k))
        w3dw4(1,1,k) = rsati(w3dw1(1,1,k), w3dw3(1,1,k))
        w3dw5(1,1,k) = exbm(1,1,k)
      enddo ! k


      ! Call the actual COAMPS microphysics scheme
      call adjtq
     .       (cond,p3,qc3,qi3,qr3,qg3,qs3,qv3,th3,w3,w3dw1,w3dw2,w3dw3
     .       ,w3dw4,w3dw5,wtm,exbm,rbm,nc3,nr3,ncn3,ni3,cp,deltf,Lf
     .       ,Ls,Lv,pcut,p0,Rd,Rv,sloper,slopes,slopeg
     .       ,timea,lice,nne,kk,i1d,j1d,ary1d,i1dflg,n1d,maxpt1d
     .       ,maxvr1d,kmax,nrdamp,ipts,kpts,icomp,kcomp,j,xland
     .       ,aa0,aa1,aa2,aa3,abar,apr,aprpr,bsnow,cbeta,cnzero,cimass
     .       ,cpor,cw,difvap,erc,esi,eic,eri,egc,esc,esr,egi,egr,egs
     .       ,mw,pi,praut1,praut2,rholiq,rhosno,rnzero,snzero,gnzero
     .       ,therco,tice,tvr1,tvr2,tvr3,tvr4,tzero,visair,gm3,gm4,gm5
     .       ,gm6,gm7,gm8,gm9,gmbp3,gmbov2,gmbov2g,bgrp,ex1,ex2,ex2g
     .       ,ex3,hlvoka,hkaolf,hlsoka,hlvorv,hlsorv,rvochi,lfocp,lvocp
     .       ,lsocp,ex7,ex7g,ex4,ex4g,ex5,ldrizzle,lgrpl,icon,icond,len
     .       ,icase
     .       ,snowv,snslp
     .       ,pcond,psmlti,psacw,pgacw,piacw,pchomo,praut,pracw
     .       ,pdepi,pint,pgdep,pconv,psaci,pgaci,praci,prevp,psdep
     .       ,pmltge,pgmlt,psmlt,pgacrm,pgacwm,pracs,pgshr,pgacr
     .       ,psacr,piacr,prhomo,pgacs,pmltse,pwacs
     .       ,falli,fallg,fallr)


      ! Assure positive definiteness in nc3/nr3/ncn3 fields
      do k=1,grnnzp
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

        nc31d(k) = nc3(1,1,k)
        nr31d(k) = nr3(1,1,k)
        ncn31d(k)= ncn3(1,1,k)
        ni31d(k) = ni3(1,1,k)
      enddo


      ! Compute tendencies
      do k=1,grnnzp
        qt3_fortend(1,1,k) = qr3(1,1,k) + qc3(1,1,k) + qs3(1,1,k) +
     .                     qi3(1,1,k) + qg3(1,1,k) + qv3(1,1,k)
        rrtend(1,1,k)    = (qr3(1,1,k) - rrbefore(1,1,k)) / deltf
        rgtend(1,1,k)    = (qg3(1,1,k) - rgbefore(1,1,k)) / deltf
        ritend(1,1,k)    = (qi3(1,1,k) - ribefore(1,1,k)) / deltf
        rsnowtend(1,1,k) = (qs3(1,1,k) - rsnowbefore(1,1,k))
     .                     / deltf
        rttend(1,1,k)    = ((qv3(1,1,k) - rvbefore(1,1,k)) / deltf)
     .                     + ((qc3(1,1,k) - rcbefore(1,1,k)) / deltf)
        thlmtend(1,1,k) = (th3(1,1,k) -
     .                    (Lv / (Cp * exbm(1,1,k))) * qc3(1,1,k)
     .                    - thmbefore(1,1,k)
     .                    + (Lv / (Cp * exbm(1,1,k))) * rcbefore(1,1,k))
     .                    / deltf
      end do ! k

      return
      end subroutine coamps_micro_driver
