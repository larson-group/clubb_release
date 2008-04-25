!----------------------------------------------------------------------
      subroutine coamps_micro_driver
     .       ( runtype, kk, kmax, timea_in, deltf_in,
     .         rtm, wmm, p, exner, rhot,
     .         thlm, ricem, rrm, rgraupelm, rsnowm,
     .         rcm, Ncm, Nrm, Ncnm, Nim,
     .         cond, Vsnow, Vice, Vrr, VNr, Vgraupel,
     .         ritend, rrtend, rgtend, 
     .         rsnowtend, nrmtend,
     .         rttend, thlmtend )

!    Description:
!      Subroutine to compute ice, as it is done in COAMPS, using
!      adjtq.F.
!      This subroutine assumes the HOC convention that k=1 is at the
!      surface rather than the top of the model domain.  adjtq.F will
!      work whichever way it is called, since it does no advection or
!      sedimentation, but any COAMPS code incorporated into this driver
!      DOES need attention, since COAMPS assumes k=1 to be the top
!      rather than the bottom.
!
!    References:
!      Rutledge and Hobbs, 1984; COAMPS Users Guide.
!----------------------------------------------------------------------
                              
      use constants
      use grid_class, only: zt2zm
      use saturation, only: rsat, rsati
      use stats_prec, only: time_prec

#ifdef STATS
      use stats_hoc, only: zt, zm, lstats_samp,
     .    imean_vol_rad_rain,
     .    imean_vol_rad_cloud
c Addition by Adam Smith, 24 April 2008
c Adding snow particle number concentration and snowslope
     .   ,isnowslope
     .   ,iNsnowm
c End of ajsmith4's addition
#endif

      implicit none

      ! External Calls
      external :: 
     .  gamma, ! From COAMPS, and not the same gamma approx. used in HOC
     .  adjtq  ! COAMPS microphysics subroutine

      ! COAMPS parameters
        integer, parameter :: 
     .  nne = 1,  ! Horizontal domain parameter (always 1 for HOC)
     .  j   = 1,  ! Horizontal grid box (always 1 for HOC)
     .  icon = 5, ! Ice nucleation scheme; 1 = Fletcher, 2 = Meyers, 3 = Hobbs & Rangno, 4 = Cooper, 5 = Cooper/Fletcher (warm/cold)
     .  icond  = 3  ! Autoconversion; 1=Kessler, 2=Manton/Cotton, 3=K&K, 4=none
      

        logical, parameter ::
     .  lice = .true.  ! Whether to produce ice in COAMPS.
                       ! According to Jerry Schmidt of NRL,
                       !   if lice = .true., then we should
                       !   set ldrizzle = .false.
                       !   because collection of drizzle by ice
                       !   is not implemented yet.
   
      ! Local Constants
      real, parameter ::
     .  aa0 = -0.267,  ! All of these are constants set in COAMPS and used by adjtq.
     .  aa1 = 5150.,
     .  aa2 = -1.0225e6,
     .  aa3 = 7.55e7,
     .  abar = 124.1,
     .  apr = 3000.,
     .  aprpr = 2.35,
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
     .  gnzero = 4.0e6,
     .  therco = 2.43e-2,
     .  tice = 269.16,
     .  tvr1 = -0.267,
     .  tvr2 = 206.,
     .  tvr3 = -2045.,
     .  tvr4 = 9060.,
     .  tzero = 273.16,
     .  visair = 1.718e-5,
     .  bgrp = 0.66,
     .  ex1 = 0.2,
     .  pcut = 1.0e-10 ! Lower threshold for calculation in COAMPS

      integer, parameter ::
     .  n1d     = 1, ! 1d graphics parameters
     .  i1dflg  = 0, ! 1d graphics parameters
! Michael Falk, 17 May 2007
     .  maxpt1d = 200, ! 1d graphics parameters
     .  maxvr1d = 200, ! 1d graphics parameters
! eMFc
     .  ipts    = 1  ! Number of COAMPS points in a model height
                     ! that contain liquid or ice (AJS)

      real, dimension(n1d), parameter ::
     .  i1d = (/0./), ! 1d graphics parameters
     .  j1d = (/0./)  ! 1d graphics parameters

      real, dimension(1,1), parameter ::
     .  xland = 0.0, ! Land/Sea assumption
     .  wtm   = 1.0  ! Weighting array for mass point (never used)

      ! Input Variables
      character(len=*), intent(in) :: runtype

      integer, intent(in) ::
     .  kk,    ! Number of COAMPS m gridpoints in the vertical (gr%nnzp-1)
     .  kmax   ! Maximum array size (kk + ??)

      ! Note: Time variables "timea_in" and "deltf_in" need to be passed in 
      !       from CLUBB with precision "time_prec".  I have redefined 
      !       "timea" and "deltf" below to be regular precision variables
      !       that are passed throughout COAMPS microphysics.  Brian; 4/5/2008.
      real(kind=time_prec), intent(in) ::
     .  timea_in,        ! Current model time                   [s]
     .  deltf_in         ! Timestep (i.e. dtmain in CLUBB)      [s]

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
      ! Nrm is now in kg^-1.  Brian.  Sept. 8, 2007.
!     .  Nrm        ! Number of rain drops       [count/m^3]
     .  Nrm        ! Number of rain drops       [count/kg]

      real, dimension(kk+1), intent(inout) ::
      ! Ncm is now in kg^-1.  Brian.  Sept. 8, 2007.
!     .  Ncm,       ! Number of cloud droplets   [count/m^3]
     .  Ncm,       ! Number of cloud droplets   [count/kg]
     .  Ncnm,      ! Number of cloud nuclei     [count/m^3]
     .  Nim        ! Number of ice crystals     [count/m^3]

c Addition by Adam Smith, 24 April 2008
c Adding snow particle number concentration
      real, dimension(kk+1) :: Nsnowm
c End of ajsmith4's addition

      ! Output Variables
      real, dimension(kk+1), intent(out) ::
     .  ritend,    ! d(ri)/dt                   [kg/kg/s]
     .  rrtend,    ! d(rr)/dt                   [kg/kg/s]
     .  rgtend,    ! d(rg)/dt                   [kg/kg/s]
     .  rsnowtend, ! d(rsnow)/dt                [kg/kg/s]
     .  rttend,    ! d(rt)/dt                   [kg/kg/s]
     .  thlmtend,  ! d(thlm)/dt                 [K/s]
      ! Nrm is now in kg^-1.  Brian.  Sept. 8, 2007.
!     .  nrmtend    ! d(Nrm)/dt                  [count/m^3/s]
     .  nrmtend    ! d(Nrm)/dt                  [count/kg/s]

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
     .  th2t3d,  ! Mean exner function on w grid        [-]
! Michael Falk, 13 Jul 2007, added these.  They are flipped versions of these variables;
! that is to say, in the following versions k=1 is the top of the domain and k=kk+1 is the
! sub-ground ghost point.
     .  w3_flip,
     .  pr3d_flip,
     .  qsatv3d_flip,
     .  temp3d_flip,
     .  qsati3d_flip,
     .  th2t3d_flip
! eMFc

      ! Variables on the m grid
      real, dimension(1,1,kk) ::
     .  qi3,   ! Pristine ice mixing ratio               [kg/kg]
     .  qr3,   ! Rain water mixing ratio                 [kg/kg]
     .  qg3,   ! Graupel mixing ratio                    [kg/kg]
     .  qs3,   ! Snow water mixing ratio                 [kg/kg]
! Michael Falk, 13 Jul 2007, added these; they are flipped versions, as above.
     .  qi3_flip,
     .  qr3_flip,
     .  qg3_flip,
     .  qs3_flip
! eMFc

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
     .  eic,
     .  sat

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
     .  snowslope,  ! These variables are the individual microphysical terms.  Michael Falk
                    !  Michael Falk
                    ! (Name modified by Adam Smith, 24 April 2008)
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
     .  pwacs,      ! ???
! Michael Falk, 13 Jul 2007, added these; they are flipped versions of these variables
     .  cond_flip,
     .  p3_flip,
     .  qc3_flip,
     .  qt3_flip,
     .  qv3_flip,
     .  th3_flip,
     .  exbm_flip,
     .  rbm_flip
! eMFc

      real, dimension(1,1,kk+1) ::
     .  fallr,  ! Fall speed for rain mixing ratio              [m/s]
     .  falln,  ! Fall speed for rain drop number conc.         [m/s]
     .  falli,  ! Fall speed for pristine ice mixing ratio      [m/s]
     .  falls,  ! Fall speed for snow mixing ratio              [m/s] (Michael Falk 19 Jul 2007)
     .  snowv,  ! Fall speed for snow mixing ratio              [m/s]
     .  fallg,  ! Fall speed for graupel mixing ratio           [m/s]
! Michael Falk, 13 Jul 2007, added these; they are flipped versions of these variables
     .  fallg_flip,
     .  fallr_flip,
     .  falln_flip,
     .  falli_flip,
     .  falls_flip,
     .  snowv_flip,
! and then added these following variables 18 Jul 2007.
! The in_cloud versions are only assigned values (and passed to and from adjtq) within clouds.
! They are passed to fallr, falln, falli, and falls, which are assigned 0 outside of cloud
! and assigned the in_cloud values within cloud.
! snowv is passed at all gridpoints, 1 to kk+1, regardless of cloud.
     .  fallr_in_cloud,
     .  falln_in_cloud,
     .  falli_in_cloud,
     .  falls_in_cloud,
     .  fallg_in_cloud
! eMFc

      real, dimension(1,1,kk) ::
     .  rvc,     ! Cloud droplet radius         [cm]
     .  rvr      ! Rain droplet radius          [cm]

      real, dimension(1,kk,1) ::
     .  ary1d        ! 1d graphics parameters

      integer ::
     .  i,
     .  k, kp1, km1, ! Loop control variables
     .  len,         ! # of saturated points???
     .  icase,       ! Which case?
     .  nrdamp       ! Number of points in the damping layer (COAMPS)

      real ::
     .  snzero ! ???

      integer, dimension(1) ::
     .  nkpts

      integer, dimension(kk) ::
     .  icomp,         !
     .  kcomp          !

      logical ::
     .  ldrizzle,      ! is drizzle on?
     .  lgrpl          ! is graupel on?


      real :: ! Regular precision to be passed into COAMPS micro. Brian.
     .  timea,        ! Current model time                      [s]
     .  deltf         ! Timestep (i.e. dtmain in CLUBB)         [s]

!----------------------------------------------------------------------

      ! Begin coamps_micro_driver code

      ! Comment by Adam Smith, 25 March 2008
      ! These variables activate rain/drizzle and graupel in the COAMPS
      ! micro scheme.  Set these variables to .TRUE. if you need these
      ! hydrometeors in your simulations.
      ldrizzle = .FALSE.
      lgrpl    = .FALSE.

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Comment by Adam Smith, 25 March 2008
      ! The variable "icase" is used in the COAMPS namelist to identify
      ! the simulation being run.  Normally, each icase value represents
      ! an individual case, in order to allow the user to use specific
      ! model settings and forcings.
      !
      ! Ideally, we would use the character string "runtype" to identify
      ! a simulation, since the variable is already used in hoc.F.
      ! However, with the current structure of CLUBB, it is not realistic
      ! to pass the variable into this driver.  Therefore, we must hard-code
      ! the icase value for now.  Perhaps in the future, it will be possible
      ! to have the "icase" value set automatically based on the cloud case.
      !
      ! Specific icase values used in COAMPS:
      ! icase = 75:   M-PACE: period B
      ! icase = 1000: Nov.11, 1999 altocu
      ! icase = 1001: Jun.25 multilayered altocu
      ! icase = 1002: CLEX-9: Oct.14, 2001 altocu
      ! icase = 1003: CLEX-9: Nov.02, 2001 altocu
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! End of ajsmith4's comment
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! dschanen on 27 Mar 2008 made icase as function of runtype
      select case ( trim( runtype ) )
      case ( "mpace_a", "mpace_b" )
        icase  = 75
        nrdamp = 15

      case ( "nov11_altocu" )
        icase = 1000
        nrdamp = 15

      case ( "jun25_altocu" )
        icase = 1001
        nrdamp = 15

      case ( "clex9_oct14" )
        icase = 1002
        nrdamp = 15

      case ( "clex9_nov02" )
        icase = 1003
        nrdamp = 15

      case default
        icase = -1 ! No idea what to do in this case. -dschanen
        nrdamp = 15

      end select

      ! Something related to ice.  In regular COAMPS it is a parameter but is
      ! then set to another value, which is illegal.
      eic      = 1.0

      ! Brian set regular precision variables "timea" and "deltf" to the values
      ! brought in from CLUBB as precision "time_prec" variables "timea_in" and
      ! "deltf_in".  Brian; 4/5/2008.
      timea = timea_in
      deltf = deltf_in

! Michael Falk changed this, November 2007.
! This is the stock COAMPS value, 2e7:
      if ( trim( runtype ) == "mpace_a" ) then
!     if ( .false. ) then
        ! and this is the value, 2e6, which works for MPACE-A:
        snzero = 2.0e6
      else ! don't cheat
        snzero = 2.0e7
      end if
! eMFc

      ! Set up initial fields

      ! Compute quantities for computing tendencies
      rvm = rtm - rcm

      if ( any( rvm < 0. ) ) then
        write(fstderr,*) "in COAMPS (R) micro driver rvm < 0"
        where ( rvm < 0. ) rvm = 0.
      end if

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

      ! Setup COAMPS m (mass) grid variables that are only used from
      ! gridpoints 1 to kk, but are nonetheless defined from 1 to kk+1
      ! since they are defined as w (momentum) variables elsewhere 
      ! within COAMPS.
      ! Since values 1 to kk are the only ones used, they are the only
      ! ones that are assigned.
      ! When setting up COAMPS variables, we remove the sub-ground ghost
      ! point (at k=1) from the HOC variables.  Since in HOC k=2 is the
      ! first above-ground gridpoint and in COAMPS k=1 is the first
      ! above-ground gridpoint, we assign the variables accordingly.
      ! Comments by Michael Falk, David Schanen, and Vince Larson

      ! The top point is undefined and unreferenced in these '3d' arrays
      pr3d(1,1,1:kk)   = p(2:kk+1) 
      th2t3d(1,1,1:kk) = exner(2:kk+1)
      temp3d(1,1,1:kk) = thm(2:kk+1) * exner(2:kk+1)

      do k=1, kk, 1
        qsatv3d(1,1,k) = rsat( pr3d(1,1,k), temp3d(1,1,k) )
        qsati3d(1,1,k) = rsati( pr3d(1,1,k), temp3d(1,1,k) )
      end do

      ! Setup COAMPS m (mass) grid variables
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
        ! Nrm and Ncm, which are in kg^-1, need to be multiplied by rhot 
        ! to be in units of m^-3, and then converted to cm^-3.
        ! Brian.  Sept. 8, 2007.
!        nc3(1,1,k)  = Ncm(k+1) * 1.e-6
!        nr3(1,1,k)  = Nrm(k+1) * 1.e-6
        nc3(1,1,k)  = ( Ncm(k+1) * rhot(k+1) ) * 1.e-6
        nr3(1,1,k)  = ( Nrm(k+1) * rhot(k+1) ) * 1.e-6
        ncn3(1,1,k) = Ncnm(k+1) * 1.e-6

        ni3(1,1,k)  = Nim(k+1)

      end do

      ! Note that this is much simpler approx. of gamma than the ANL function
      ! used in the rest of the CLUBB code.
      ! It is used only to have a direct comparison with COAMPS-LES -dschanen
      call gamma( 3.0,gm3 )
      call gamma( 4.0,gm4 )
      call gamma( 5.0,gm5 )
      call gamma( 6.0,gm6 )
      call gamma( 7.0,gm7 )
      call gamma( 8.0,gm8 )
      call gamma( 9.0,gm9 )
      call gamma( bsnow+3.0, gmbp3 )
      call gamma( bsnow*0.5 + 2.5, gmbov2 )
      call gamma( bgrp*0.5 + 2.5, gmbov2g )

      ex2  = bsnow * 0.5 + 2.5
      ex2g = bgrp * 0.5 + 2.5
      ex3  = bsnow + 3.0
      ex7  = 0.44 * gmbov2
      ex7g = 0.31 * gmbov2g
      ex4  = aprpr/visair
      ex4g = abar/visair
      ex5  = pi*aprpr*snzero*gmbp3/4.0

!     Lf     = Ls - Lv   ! Latent heat of fusion (occurs in constants)
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
      nkpts = 0
      len = 0

!!! Michael Falk, 6 July 2007
!!! Assigning values to the "flip" variables, in which k=1 is the top of the domain
      
!     cond_flip(1,1,1:kk) = cond(1,1,kk:1:-1) ! this is uninitialized -dschanen
      p3_flip(1,1,1:kk)  = p3(1,1,kk:1:-1)     ! Pressure
      qt3_flip(1,1,1:kk) = qt3(1,1,kk:1:-1)   ! Total water mixing ratio
      qv3_flip(1,1,1:kk) = qv3(1,1,kk:1:-1)   ! Vapor water mixing ratio
      qc3_flip(1,1,1:kk) = qc3(1,1,kk:1:-1)   ! Cloud water mixing ratio
      th3_flip(1,1,1:kk) = th3(1,1,kk:1:-1)   ! Potential temp.
      exbm_flip(1,1,1:kk)= exbm(1,1,kk:1:-1) ! Exner function
      rbm_flip(1,1,1:kk) = rbm(1,1,kk:1:-1)   ! Density

      qi3_flip(1,1,1:kk) = qi3(1,1,kk:1:-1) ! Ice water mixing ratio
      qr3_flip(1,1,1:kk) = qr3(1,1,kk:1:-1) ! Rain water mixing ratio
      qg3_flip(1,1,1:kk) = qg3(1,1,kk:1:-1) ! Graupel water mixing ratio
      qs3_flip(1,1,1:kk) = qs3(1,1,kk:1:-1) ! Snow water mixing ratio

      w3_flip(1,1,1:kk+1) = w3(1,1,kk+1:1:-1) ! not referenced in COAMPS micro

      pr3d_flip(1,1,1:kk)    = pr3d(1,1,kk:1:-1) ! top point is undefined
      qsatv3d_flip(1,1,1:kk) = qsatv3d(1,1,kk:1:-1) ! " " 
      temp3d_flip(1,1,1:kk)  = temp3d(1,1,kk:1:-1) ! " " 
      qsati3d_flip(1,1,1:kk) = qsati3d(1,1,kk:1:-1) ! " " 
      th2t3d_flip(1,1,1:kk)  = th2t3d(1,1,kk:1:-1) ! " " 



! determination of which points are "in_cloud" are which are outside of cloud
      do k=nrdamp,kk
c
c***********************************************************************
c determine if point is saturated as in COAMPS
c***********************************************************************
c
      if ( .not.lice ) then
        sat = qv3_flip(1,1,k)/qsatv3d_flip(1,1,k)-1.0
      else
        if ( temp3d_flip(1,1,k) >= 273.15 ) then
          sat = qv3_flip(1,1,k)/qsatv3d_flip(1,1,k)-1.0
        else
          sat = qv3_flip(1,1,k)/qsati3d_flip(1,1,k)-1.0
        endif
      endif

        if (sat > 0.0 .or.
     .      qc3_flip(1,1,k) .ge. pcut .or.
     .      qr3_flip(1,1,k) .ge. pcut .or.
     .      qs3_flip(1,1,k) .ge. pcut .or.
     .      qi3_flip(1,1,k) .ge. pcut .or.
     .      qg3_flip(1,1,k) .ge. pcut
     .     )
     .      then
          nkpts = nkpts+1
          len   = len+1
          kcomp(nkpts) = k
        end if
      end do

      do i=1,ipts
        icomp(i) = 1
      end do

      sloper = pi * rholiq * rnzero * 1.0e-8
      slopes = pi * rhosno * snzero * 1.0e-8
      slopeg = pi * rhogrp * gnzero * 1.0e-8

! Michael Falk, 17 Jul 2007, is initializing fallspeed arrays
      do k=1,kk+1
        falli(1,1,k) = 0.
        falls(1,1,k) = 0.
        fallg(1,1,k) = 0.
        fallr(1,1,k) = 0.
        falln(1,1,k) = 0.
        snowv(1,1,k) = 0.

        falli_flip(1,1,k) = 0.
        falls_flip(1,1,k) = 0.
        fallg_flip(1,1,k) = 0.
        fallr_flip(1,1,k) = 0.
        falln_flip(1,1,k) = 0.

! Michael Falk added these initializations, 10 Oct 2007
        falli_in_cloud(1,1,k) = 0.
        falls_in_cloud(1,1,k) = 0.
        fallg_in_cloud(1,1,k) = 0.
        fallr_in_cloud(1,1,k) = 0.
        falln_in_cloud(1,1,k) = 0.
! eMFc
      end do
! eMFc

      ! Call the actual COAMPS microphysics scheme
      if ( len > 0 ) call adjtq
     1       (cond_flip(1,1,1:kk),p3_flip(1,1,1:kk),qc3_flip(1,1,1:kk)
     1       ,qi3_flip(1,1,1:kk),qr3_flip(1,1,1:kk),qg3_flip(1,1,1:kk)
     1       ,qs3_flip(1,1,1:kk),qv3_flip(1,1,1:kk),th3_flip(1,1,1:kk)
     1       ,w3_flip(1,1,1:kk+1),pr3d_flip(1,1,1:kk+1)
     1       ,qsatv3d_flip(1,1,1:kk+1),temp3d_flip(1,1,1:kk+1)
     1       ,qsati3d_flip(1,1,1:kk+1),th2t3d_flip(1,1,1:kk+1),wtm
     2       ,exbm_flip(1,1,1:kk),rbm_flip(1,1,1:kk)
!     3       ,nc3,nr3,ncn3,ni3,cp,deltf,Lf,Ls,Lv 
     3       ,nc3(1,1,kk:1:-1),nr3(1,1,kk:1:-1),ncn3(1,1,kk:1:-1)
     3       ,ni3(1,1,kk:1:-1),cp,deltf,Lf,Ls,Lv
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
     6       ,ldrizzle,lgrpl,icon,icond,len,icase
     6       ,snowv_flip(:,:,1:kk)
     6       ,snowslope(1,1,kk:1:-1),pcond(1,1,kk:1:-1)
     7       ,psmlti(1,1,kk:1:-1),psacw(1,1,kk:1:-1),pgacw(1,1,kk:1:-1)
     7       ,piacw(1,1,kk:1:-1),pchomo(1,1,kk:1:-1),praut(1,1,kk:1:-1)
     7       ,pracw(1,1,kk:1:-1),pdepi(1,1,kk:1:-1),pint(1,1,kk:1:-1)
     8       ,pgdep(1,1,kk:1:-1),pconv(1,1,kk:1:-1),psaci(1,1,kk:1:-1)
     8       ,pgaci(1,1,kk:1:-1),praci(1,1,kk:1:-1),prevp(1,1,kk:1:-1)
     9       ,psdep(1,1,kk:1:-1),pmltge(1,1,kk:1:-1),pgmlt(1,1,kk:1:-1)
     9       ,psmlt(1,1,kk:1:-1),pgacrm(1,1,kk:1:-1),pgacwm(1,1,kk:1:-1)
     1       ,pracs(1,1,kk:1:-1),pgshr(1,1,kk:1:-1),pgacr(1,1,kk:1:-1)
     1       ,psacr(1,1,kk:1:-1),piacr(1,1,kk:1:-1),prhomo(1,1,kk:1:-1)
     1       ,pgacs(1,1,kk:1:-1),pmltse(1,1,kk:1:-1),pwacs(1,1,kk:1:-1)
     2       ,falli_in_cloud(:,:,1:nkpts(1))
     2       ,fallg_in_cloud(:,:,1:nkpts(1))
     3       ,fallr_in_cloud(:,:,1:nkpts(1))
     3       ,falln_in_cloud(:,:,1:nkpts(1))
     3       ,falls_in_cloud(:,:,1:nkpts(1))
     4       ,rvc,rvr)

      ! reassigning flipped versions of variables to normal versions
      cond(1,1,1:kk) = cond_flip(1,1,kk:1:-1)
      p3(1,1,1:kk) = p3_flip(1,1,kk:1:-1)
      qt3(1,1,1:kk) = qt3_flip(1,1,kk:1:-1)
      qv3(1,1,1:kk) = qv3_flip(1,1,kk:1:-1)
      qc3(1,1,1:kk) = qc3_flip(1,1,kk:1:-1)
      th3(1,1,1:kk) = th3_flip(1,1,kk:1:-1)
      exbm(1,1,1:kk) = exbm_flip(1,1,kk:1:-1)
      rbm(1,1,1:kk) = rbm_flip(1,1,kk:1:-1)

      qi3(1,1,1:kk) = qi3_flip(1,1,kk:1:-1)
      qr3(1,1,1:kk) = qr3_flip(1,1,kk:1:-1)
      qg3(1,1,1:kk) = qg3_flip(1,1,kk:1:-1)
      qs3(1,1,1:kk) = qs3_flip(1,1,kk:1:-1)

      ! This is unneeded, since these quantities do not change -dschanen
      ! 2 April 2008
!     w3(1,1,1:kk+1) = w3_flip(1,1,kk+1:1:-1)
!     pr3d(1,1,1:kk+1) = pr3d_flip(1,1,kk+1:1:-1)
!     qsatv3d(1,1,1:kk+1) = qsatv3d_flip(1,1,kk+1:1:-1)
!     temp3d(1,1,1:kk+1) = temp3d_flip(1,1,kk+1:1:-1)
!     qsati3d(1,1,1:kk+1) = qsati3d_flip(1,1,kk+1:1:-1)
!     th2t3d(1,1,1:kk+1) = th2t3d_flip(1,1,kk+1:1:-1)

      ! assigning in-cloud fall speeds to _flip arrays, which are zero outside of cloud
      do k=1,nkpts(1)
        falli_flip(:,:,kcomp(k)) = falli_in_cloud(:,:,k)
        falls_flip(:,:,kcomp(k)) = falls_in_cloud(:,:,k)
        fallg_flip(:,:,kcomp(k)) = fallg_in_cloud(:,:,k)
        fallr_flip(:,:,kcomp(k)) = fallr_in_cloud(:,:,k)
        falln_flip(:,:,kcomp(k)) = falln_in_cloud(:,:,k)
      end do

      falli(:,:,kk:1:-1) = falli_flip(:,:,1:kk)
      falls(:,:,kk:1:-1) = falls_flip(:,:,1:kk)
      fallg(:,:,kk:1:-1) = fallg_flip(:,:,1:kk)
      fallr(:,:,kk:1:-1) = fallr_flip(:,:,1:kk)
      falln(:,:,kk:1:-1) = falln_flip(:,:,1:kk)
      snowv(:,:,kk:1:-1) = falls_flip(:,:,1:kk)

      ! Assure positive definiteness in nc3/nr3/ncn3 fields

      ! Should there be a clipping stat for these?
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
        ! nc3, in cm^-3, needs to be converted to m^-3, and then divided by
        ! rhot so that Ncm is in kg^-1.  Brian.  Sept. 8, 2007
!        Ncm(k+1)  = nc3(1,1,k) * 1.e6
        Ncm(k+1)  = ( nc3(1,1,k) * 1.e6 ) / rhot(k+1)
        Ncnm(k+1) = ncn3(1,1,k) * 1.e6
        Nim(k+1)  = ni3(1,1,k)
      end do ! k=1..kk+1

c*******************************************
c Addition by Adam Smith, 24 April 2008
c Adding snow particle number concentration
c Values of snowslope < 1.0 lead to excessive and
c unrealistic Nsnowm outside of the snow region.
c The "if" statement prevents these results.
c*******************************************
      do k = 1, kk, 1
        if (snowslope(1,1,k) < 2.0) then
          snowslope(1,1,k) = 0.
          Nsnowm(k+1) = 0.0
        else
          Nsnowm(k+1) = snzero / snowslope(1,1,k)
        end if
      end do
c****************************
c End of ajsmith4's addition
c****************************

      ! Linear extrapolation for the ghost point of fall speeds
      fallr(1,1,1) = .5 * ( fallr(1,1,2) + fallr(1,1,3) )
      falln(1,1,1) = .5 * ( falln(1,1,2) + falln(1,1,3) )
      snowv(1,1,1) = .5 * ( snowv(1,1,2) + snowv(1,1,3) )
      falli(1,1,1) = .5 * ( falli(1,1,2) + falli(1,1,3) )
      fallg(1,1,1) = .5 * ( fallg(1,1,2) + fallg(1,1,3) )
      falls(1,1,1) = .5 * ( falls(1,1,2) + falls(1,1,3) )

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
        ! Nrm is now in kg^-1, so nr3(1,1,k)*1.e6 needs to be divided by rhot
        ! in order to be in the same units as Nrm.  This will cause nrmtend to
        ! have units of kg^-1 s^-1, which is what we want.  
        ! Brian.  Sept. 8, 2007.
        !nrmtend(k+1)   = ( nr3(1,1,k)*1.e6 - Nrm(k+1) ) / deltf ! Conversion factor
        nrmtend(k+1)   = ( (nr3(1,1,k)*1.e6)/rhot(k+1) - Nrm(k+1) ) 
     .                   / deltf ! Conversion factor
        rsnowtend(k+1) = ( qs3(1,1,k) - rsnowm(k+1) ) / deltf
        rttend(k+1)    = ((qv3(1,1,k) - rvm(k+1)) / deltf)
     .                 + ((qc3(1,1,k) - rcm(k+1)) / deltf)
        thlmtend(k+1) 
     .  = ( ( th3(1,1,k) - (Lv / (Cp * exbm(1,1,k)) * qc3(1,1,k) ) )
     .      - thlm(k+1) ) / deltf
      end do ! k=1..kk

      rrtend(1)    = 0.0
      rgtend(1)    = 0.0
      ritend(1)    = 0.0
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

c Addition by Adam Smith, 24 April 2008
c Adding calculation for snow particle number concentration       
        if (isnowslope > 0) then
          do k = 2,kk,1
            zt%x(k,isnowslope) = zt%x(k,isnowslope) + snowslope(1,1,k)
          end do
        end if

c Addition by Adam Smith, 25 April 2008
c Adding calculation for snow particle number concentration       
        if (iNsnowm > 0) then
          zt%x(2:,iNsnowm) = zt%x(2:,iNsnowm) + Nsnowm(:)
          zt%n(2:,iNsnowm) = zt%n(2:,iNsnowm) + 1
        end if
c End of ajsmith4's additions

      end if
#endif /*STATS*/

      return
      end subroutine coamps_micro_driver
