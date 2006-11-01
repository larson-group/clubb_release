!----------------------------------------------------------------------
      subroutine coamps_micro_driver
     .       (grnnzp,p3,w3,w3dw1,wtm,exbm,rbm,timea,
     .        thlm,deltf,nne,kk,kmax,nrdamp,j,icon,icond,
     .        len,icase,lice,ldrizzle,lgrpl,
     .        qc3,qi3,qr3,qg3,qs3,qv3,th3,
     .        nc3,nr3,ncn3,ni3,pcut,
     .        cond,snowv,falli,fallr,fallg,thlmtend,
     .        ritend,rrtend,rgtend,rsnowtend,rttend)

!    Description:
!      Subroutine to compute ice, as it is done in COAMPS, using
!      adjtq.F.
!
!    References:
!      Rutledge and Hobbs, 1984; COAMPS Users Guide.
!----------------------------------------------------------------------
                              
      use constants
      implicit none

!----------------------------------------------------------------------
! INPUT VARIABLES
!----------------------------------------------------------------------

      integer, intent(in) :: grnnzp            ! number of gridpoints in the vertical

      real, intent(in)    :: p3(1,1,grnnzp)    ! perturbation exner function
      real, intent(in)    :: w3(1,1,grnnzp)    ! vertical wind
      real, intent(in)    :: w3dw1(1,1,grnnzp) ! aka pr3d    -- pressure array?
      real, intent(in)    :: wtm               ! weighting array for mass point (never used)
      real, intent(in)    :: exbm(1,1,grnnzp)  ! mean exner function
      real, intent(in)    :: rbm(1,1,grnnzp)   ! mean density
      real, intent(in)    :: timea             ! Output time
      real, intent(in)    :: thlm(1,1,grnnzp)  ! liquid water potential temperature
      real, intent(in)    :: deltf             ! Timestep, I think

      integer, intent(in)     :: nne           ! has to do with horiz domain; 1 for HOC?
      integer, intent(in)     :: kk            ! equal to grnnzp
      integer, intent(in)     :: kmax          ! maximum array size (grnnzp + ??)
      integer, intent(in)     :: nrdamp        ! Number of upper model levels in sponge layer (8?)
      integer, intent(in)     :: j             ! j-coordinate???
      integer, intent(in)     :: icon          ! ice nucleation type
      integer, intent(in)     :: icond         ! autoconversion; 1=Kessler, 2=Manton/Cotton, 3=K&K, 4=none
      integer, intent(in)     :: len           ! # of saturated points???
      integer, intent(in)     :: icase         ! Which case?

      logical, intent(in)     :: lice          ! Logical: is ice on?
      logical, intent(in)     :: ldrizzle      ! is drizzle on?
      logical, intent(in)     :: lgrpl         ! is graupel on?

!----------------------------------------------------------------------
! INPUT/OUTPUT VARIABLES
!----------------------------------------------------------------------
      real, intent(inout) :: qc3(1,1,grnnzp)   ! cloud water mixing ratio
      real, intent(inout) :: qi3(1,1,grnnzp)   ! pristine ice mixing ratio
      real, intent(inout) :: qr3(1,1,grnnzp)   ! rain water mixing ratio
      real, intent(inout) :: qg3(1,1,grnnzp)   ! graupel mixing ratio
      real, intent(inout) :: qs3(1,1,grnnzp)   ! snow water mixing ratio
      real, intent(inout) :: qv3(1,1,grnnzp)   ! water vapor mixing ratio
      real, intent(inout) :: th3(1,1,grnnzp)   ! potential temperature
      real, intent(inout) :: nc3(1,1,grnnzp)   ! Number of cloud droplets
      real, intent(inout) :: nr3(1,1,grnnzp)   ! Number of rain drops
      real, intent(inout) :: ncn3(1,1,grnnzp)  ! Number of cloud nuclei
      real, intent(inout) :: ni3(1,1,grnnzp)   ! Number of ice crystals
      real, intent(inout) :: pcut              ! low threshold for mixing ratio

!----------------------------------------------------------------------
! OUTPUT VARIABLES
!----------------------------------------------------------------------
      real, intent(out)     :: cond(1,1,grnnzp)     ! condensation/evaporation of liquid water

      real, intent(out)     :: snowv(1,1,grnnzp)  ! snow fall speed
      real, intent(out)     :: falli(1,1,grnnzp)  ! pristine ice fall speed
      real, intent(out)     :: fallr(1,1,grnnzp)  !
      real, intent(out)     :: fallg(1,1,grnnzp)  !

      real, intent(out)     :: ritend(1,1,grnnzp)    ! d(ri)/dt
      real, intent(out)     :: rrtend(1,1,grnnzp)    ! d(rr)/dt
      real, intent(out)     :: rgtend(1,1,grnnzp)    ! d(rg)/dt
      real, intent(out)     :: rsnowtend(1,1,grnnzp) ! d(rsnow)/dt
      real, intent(out)     :: rttend(1,1,grnnzp)    ! d(rt)/dt

      real, intent(out)     :: thlmtend(1,1,grnnzp)  ! d(thlm)/dt

!----------------------------------------------------------------------
! EXTERNAL FUNCTIONS
!----------------------------------------------------------------------
      real, external :: gamma
      real, external :: rsat
      real, external :: rsati

!----------------------------------------------------------------------
! INTERNAL VARIABLES
!----------------------------------------------------------------------
      integer i,k ! for loops

      real w3dw2(1,1,grnnzp) ! aka qsatv3d -- saturation mr array?
      real w3dw3(1,1,grnnzp) ! aka temp3d  -- temperature array?
      real w3dw4(1,1,grnnzp) ! aka qsati3d -- saturation mr over ice array?
      real w3dw5(1,1,grnnzp) ! aka th2t3d  -- exner function array?
      real thlmtemp(1,1,grnnzp) ! liquid water potential temperature

!   -------------------------------------------------------------------
!   Constants set at various places with COAMPS
!   -------------------------------------------------------------------
      real, parameter ::       aa0 = -0.267
      real, parameter ::       aa1 = 5150.
      real, parameter ::       aa2 = -1.0225e6
      real, parameter ::       aa3 = 7.55e7
      real, parameter ::       abar = 19.3
      real, parameter ::       apr = 3000.
      real, parameter ::       aprpr = 1.139
      real, parameter ::       bsnow = 0.11
      real, parameter ::       cbeta = 0.6
      real, parameter ::       cnzero = 0.01
      real, parameter ::       cimass = 9.4e-10
      real, parameter ::       cw = 4218.
      real, parameter ::       difvap = 2.26e-5
      real, parameter ::       erc = 1.
      real, parameter ::       esi = 0.1
      real, parameter ::       eic = 1.
      real, parameter ::       eri = 1.
      real, parameter ::       egc = 1.
      real, parameter ::       esc = 1.
      real, parameter ::       esr = 0.4
      real, parameter ::       egi = 0.1
      real, parameter ::       egr = 1.0
      real, parameter ::       egs = 0.1
      real, parameter ::       mw = 18.016
      real, parameter ::       praut1 = 0.001
      real, parameter ::       praut2 = 0.0004
      real, parameter ::       rholiq = 1000.
      real, parameter ::       rhosno = 100.
      real, parameter ::       rhogrp = 400.
      real, parameter ::       rnzero = 8.0e6
      real, parameter ::       snzero = 2.0e7
      real, parameter ::       gnzero = 4.0e6
      real, parameter ::       therco = 2.43e-2
      real, parameter ::       tice = 269.16
      real, parameter ::       tvr1 = -0.267
      real, parameter ::       tvr2 = 206.
      real, parameter ::       tvr3 = -2045.
      real, parameter ::       tvr4 = 9060.
      real, parameter ::       tzero = 273.16
      real, parameter ::       visair = 1.718e-5
      real, parameter ::       bgrp = 0.37
      real, parameter ::       ex1 = 0.2

!   -------------------------------------------------------------------
!   Parameters computed for adjtq.F
!   -------------------------------------------------------------------
      real gmbov2
      real gmbov2g
      real gmbp3
      real gm3, gm4, gm5, gm6, gm7, gm8, gm9
      real ex2,ex2g,ex3
      real ex7
      real ex7g
      real ex4
      real ex4g
      real ex5

      real Lf
      real hlvoka
      real hlsoka
      real hlvorv
      real hlsorv
      real rvochi
      real cpor
      real lfocp
      real lvocp
      real lsocp
      real hkaolf

      real sloper,slopes,slopeg

      real, parameter    :: i1d = 0 ! 1d graphics parameters
      real, parameter    :: j1d = 0 ! 1d graphics parameters
      real, parameter    :: n1d = 0 ! 1d graphics parameters
      integer, parameter :: i1dflg = 0 ! 1d graphics parameters
      integer, parameter :: maxpt1d = 0 ! 1d graphics parameters
      integer, parameter :: maxvr1d = 0 ! 1d graphics parameters
      integer, parameter :: ipts = 1 ! ??
      integer, parameter :: xland = 1 ! We're over land.
      integer kpts
      integer icomp(grnnzp)       ! ???
      integer kcomp(1,grnnzp)  ! ???
      real ary1d(1,grnnzp,1) ! 1d graphics parameters

      real rcbefore(1,1,grnnzp)
      real rsnowbefore(1,1,grnnzp)
      real ribefore(1,1,grnnzp)
      real rvbefore(1,1,grnnzp)
      real rrbefore(1,1,grnnzp)
      real rgbefore(1,1,grnnzp)
      real rtbefore(1,1,grnnzp)
      real qt3(1,1,grnnzp)
      real thlmbefore(1,1,grnnzp)

!   -------------------------------------------------------------------
!   Unused output from adjtq.F
!   -------------------------------------------------------------------
      real snslp(1,1,grnnzp)        ! ajsmith; ARRAY?
      real pcond(1,1,grnnzp)        ! condensation/evaporation of cloud water
      real psmlti(1,1,grnnzp)       ! melting of cloud ice
      real psacw(1,1,grnnzp)        ! collection of cloud water by snow
      real pgacw(1,1,grnnzp)        ! ???
      real piacw(1,1,grnnzp)        ! ???
      real pchomo(1,1,grnnzp)       ! ???
      real praut(1,1,grnnzp)        ! autoconversion of cloud water
      real pracw(1,1,grnnzp)        ! collection of cloud water by rain water
      real pdepi(1,1,grnnzp)        ! depositional growth of cloud ice
      real pint(1,1,grnnzp)         ! initiation of cloud ice
      real pgdep(1,1,grnnzp)        ! ???
      real pconv(1,1,grnnzp)        ! conversion of cloud ice to snow
      real psaci(1,1,grnnzp)        ! collection of cloud ice by snow
      real pgaci(1,1,grnnzp)        ! ???
      real praci(1,1,grnnzp)        ! ???
      real prevp(1,1,grnnzp)        ! evaporation of raindrops
      real psdep(1,1,grnnzp)        ! depositional growth of snow
      real pmltge(1,1,grnnzp)       ! ???
      real pgmlt(1,1,grnnzp)        ! ???
      real psmlt(1,1,grnnzp)        ! melting of snow
      real pgacrm(1,1,grnnzp)       ! ???
      real pgacwm(1,1,grnnzp)       ! ???
      real pracs(1,1,grnnzp)        ! ???
      real pgshr(1,1,grnnzp)        ! ???
      real pgacr(1,1,grnnzp)        ! ???
      real psacr(1,1,grnnzp)        ! ???
      real piacr(1,1,grnnzp)        ! ???
      real prhomo(1,1,grnnzp)       ! ???
      real pgacs(1,1,grnnzp)        ! ???
      real pmltse(1,1,grnnzp)       ! ???
      real pwacs(1,1,grnnzp)        ! ???

!----------------------------------------------------------------------
! COMMENCE EXECUTABLE CODE
!----------------------------------------------------------------------
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

!   -------------------------------------------------------------------
!   Saving pre-adjtq values, for computation of tendencies
!   Computing w3dw2-w3dw5 for use in call to adjtq
!   -------------------------------------------------------------------
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

        w3dw3(1,1,k) = th3(1,1,k) * exbm(1,1,k)
        w3dw2(1,1,k) = rsat(w3dw1(1,1,k), w3dw3(1,1,k))
        w3dw4(1,1,k) = rsati(w3dw1(1,1,k), w3dw3(1,1,k))
        w3dw5(1,1,k) = exbm(1,1,k)
      enddo ! k

!   -------------------------------------------------------------------
!   Actually call adjtq.F
!   -------------------------------------------------------------------

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

!   -------------------------------------------------------------------
!   Assure positive definiteness in nc3/nr3/ncn3 fields
!   -------------------------------------------------------------------
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
      end do ! k

!   -------------------------------------------------------------------
!   Compute tendencies
!   -------------------------------------------------------------------
      do k=1,grnnzp
        qt3(1,1,k)       = qr3(1,1,k) + qc3(1,1,k) + qs3(1,1,k) +
     .                     qi3(1,1,k) + qg3(1,1,k) + qv3(1,1,k)
        rrtend(1,1,k)    = (qr3(1,1,k) - rrbefore(1,1,k)) / deltf
        rgtend(1,1,k)    = (qg3(1,1,k) - rgbefore(1,1,k)) / deltf
        ritend(1,1,k)    = (qi3(1,1,k) - ribefore(1,1,k)) / deltf
        rsnowtend(1,1,k) = (qs3(1,1,k) - rsnowbefore(1,1,k))
     .                     / deltf
        rttend(1,1,k)    = (qt3(1,1,k) - rtbefore(1,1,k)) / deltf
        thlmtend(1,1,k)  = ((th3(1,1,k) -
     .                     (Lv / (Cp * exbm(1,1,k)))
     .                     * qc3(1,1,k)) - thlmbefore(1,1,k)) 
     .                     / deltf
        thlmtemp(1,1,k)  = thlm(1,1,k) + thlmtend(1,1,k)
        if (abs(thlmtend(1,1,k)) < 1e-9) then
          thlmtend(1,1,k) = 0.
        end if
      end do ! k

      return
      end subroutine coamps_micro_driver
