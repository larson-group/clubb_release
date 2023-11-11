!----------------------------------------------------------------------
! $Id$
module coamps_microphys_driver_module

! This module wraps the adjtq subroutine so that it may be used by
! CLUBB

  implicit none

  public :: coamps_microphys_driver

  private ! Set Default Scope

  contains

#ifdef COAMPS_MICRO
  subroutine coamps_microphys_driver & 
         ( gr, runtype, timea_in, deltf_in, & 
           rtm, wm_zm, p_in_Pa, exner, rho, & 
           thlm, rim, rrm, rgm, rsm, & 
           rcm, Ncm, Nrm, Nim, &
           stats_metadata, &
           stats_zt, &
           Nccnm, cond, &
           Vrs, Vri, Vrr, VNr, Vrg, & 
           ritend, rrtend, rgtend,  & 
           rstend, nrmtend, & 
           ncmtend, nimtend, & 
           rvm_mc, rcm_mc, thlm_mc )

!    Description:
!      Subroutine to compute ice, as it is done in NRL's COAMPS
!      mesoscale model, using adjtq.F.

!      COAMPS was developed by the Naval Research Laboratory,
!      Monterey, California.  COAMPS is a registered trademark of the
!      Naval Research Laboratory.

!      This subroutine assumes the CLUBB convention that k=1 is at the
!      surface rather than the top of the model domain.  adjtq.F will
!      work whichever way it is called, since it does no advection or
!      sedimentation, but any COAMPS code incorporated into this driver
!      DOES need attention, since COAMPS assumes k=1 to be the top
!      rather than the bottom.
!
!    References:
!      Rutledge and Hobbs, 1984; COAMPS Users Guide.
!----------------------------------------------------------------------

    use constants_clubb, only: &
      Cp, Lv, pi, Lf, Ls, Rv, Rd, p0, T_freeze_K, cm3_per_m3, & ! Variable(s)
      fstderr ! Constant(s)
    use saturation, only: sat_mixrat_liq, sat_mixrat_ice ! Procedure(s)
    use clubb_precision, only: time_precision, core_rknd ! Variable(s)
    use grid_class, only: zt2zm ! Procedure(s)


    use grid_class, only: grid ! Type

    use stats_type_utilities, only: stat_update_var_pt ! Procedure(s)

    use stats_variables, only: &
        stats_metadata_type

    use error_code, only: &
        clubb_at_least_debug_level  ! Procedure

    use parameters_microphys, only: l_graupel, l_ice_microphys ! Variable(s)

    use T_in_K_module, only: thlm2T_in_K ! Procedure(s)

    use stats_type, only: stats ! Type

    implicit none

! External Calls
    external ::  & 
      gamma,  & ! From COAMPS, and not the same gamma approx. used in CLUBB
      adjtq     ! COAMPS microphysics subroutine

! COAMPS parameters
    integer, parameter ::  & 
      nne = 1,   & ! Horizontal domain parameter (always 1 for CLUBB)
      j   = 1,   & ! Horizontal grid box (always 1 for CLUBB)
      icon = 5,  & ! Ice nucleation scheme; 1 = Fletcher, 
               !                        2 = Meyers, 
               !                        3 = Hobbs & Rangno, 
               !                        4 = Cooper, 
               !                        5 = Cooper/Fletcher (warm/cold)
      icond  = 3  ! Autoconversion; 1=Kessler, 2=Manton/Cotton, 3=K&K, 4=none



! Local Constants
    real, parameter :: & 
      aa0 = -0.267,   & ! All of these are constants set in COAMPS and used by adjtq.
      aa1 = 5150., & 
      aa2 = -1.0225e6, & 
      aa3 = 7.55e7, & 
      abar = 124.1, & 
      apr = 3000., & 
      aprpr = 2.35, & 
      bsnow = 0.11, & 
      cbeta = 0.6, & 
      cnzero = 0.01, & 
      cimass = 9.4e-10, & 
      cw = 4218., & 
      difvap = 2.26e-5, & 
      erc = 1., & 
      esi = 0.1, & 
      eri = 1., & 
      egc = 1., & 
      esc = 1., & 
      esr = 0.4

    real, parameter :: &
      egi = 0.1, & 
      egr = 1.0, & 
      egs = 0.1, & 
      mw = 18.016, & 
      praut1 = 0.001, & 
      praut2 = 0.0004, & 
      rholiq = 1000., & 
      rhosno = 100., & 
      rhogrp = 400., & 
      rnzero = 8.0e6, & 
      gnzero = 4.0e6, & 
      therco = 2.43e-2, & 
      tice = 269.16, & 
      tvr1 = -0.267, & 
      tvr2 = 206., & 
      tvr3 = -2045., & 
      tvr4 = 9060., & 
      tzero = 273.16

    real, parameter :: & 
      visair = 1.718e-5, & 
      bgrp = 0.66, & 
      ex1 = 0.2, & 
      pcut = 1.0e-10 ! Lower threshold for calculation in COAMPS

    integer, parameter :: & 
      n1d     = 1,  & ! 1d graphics parameters
      i1dflg  = 0,  & ! 1d graphics parameters
! Michael Falk, 17 May 2007
      maxpt1d = 200,  & ! 1d graphics parameters
      maxvr1d = 200,  & ! 1d graphics parameters
! eMFc
      ipts    = 1  ! Number of COAMPS points in a model height
    ! that contain liquid or ice (AJS)

    real, dimension(n1d), parameter :: & 
      i1d = (/0./),  & ! 1d graphics parameters
      j1d = (/0./)  ! 1d graphics parameters

    real, dimension(1,1), parameter :: & 
      xland = 0.0,  & ! Land/Sea assumption
      wtm   = 1.0  ! Weighting array for mass point (never used)

! Input Variables
    type (grid), target, intent(in) :: gr

    character(len=*), intent(in) :: runtype

! Note: Time variables "timea_in" and "deltf_in" need to be passed in
!       from CLUBB with precision "time_precision".  I have redefined
!       "timea" and "deltf" below to be regular precision variables
!       that are passed throughout COAMPS microphysics.  Brian; 4/5/2008.
    real( kind = time_precision ), intent(in) :: & 
      timea_in          ! Current model time                   [s]

    real( kind= core_rknd ), intent(in) :: &
      deltf_in         ! Timestep (i.e. dt_main in CLUBB)      [s]

    real(kind = core_rknd), dimension(gr%nz), intent(in) :: & 
      rtm,     & ! Total water mixing ratio                        [kg/kg]
      rcm,     & ! Cloud water mixing ratio                        [kg/kg]
      wm_zm,   & ! Vertical wind                                   [m/s]
      p_in_Pa, & ! Pressure                                        [Pa]
      exner,   & ! Mean exner function                             [-]
      rho,     & ! Mean density                                    [kg/m^3]
      thlm       ! Liquid potential temperature                    [K]

    real(kind = core_rknd), dimension(gr%nz), intent(in) :: & 
      rim,      & ! Ice water mixing ratio     [kg/kg]
      rrm,     & ! Rain water mixing ratio    [kg/kg]
      rgm,  & ! Graupel water mixing ratio [kg/kg]
      rsm,     & ! Snow water mixing ratio    [kg/kg]
      Nrm,        & ! Number of rain drops       [count/kg]
      Ncm,        & ! Number of cloud droplets   [count/kg]
      Nim           ! Number of ice crystals     [count/kg]

    type (stats_metadata_type), intent(in) :: &
      stats_metadata

    type(stats), target, intent(inout) :: &
      stats_zt 

    real(kind = core_rknd), dimension(gr%nz), intent(inout) :: & 
      Nccnm         ! Number of cloud nuclei     [count/kg]

    real(kind = core_rknd), dimension(1,1,gr%nz-1), intent(inout) :: &
      cond ! condensation/evaporation of liquid water

! Output Variables
    real(kind = core_rknd), dimension(gr%nz), intent(out) :: & 
      ritend,     & ! d(ri)/dt                   [kg/kg/s]
      rrtend,     & ! d(rr)/dt                   [kg/kg/s]
      rgtend,     & ! d(rg)/dt                   [kg/kg/s]
      rstend,  & ! d(rs)/dt                [kg/kg/s]
      rvm_mc,     & ! d(rv)/dt                   [kg/kg/s]
      rcm_mc,     & ! d(rc)/dt                   [kg/kg/s]
      thlm_mc,    & ! d(thlm)/dt                 [K/s]
      nrmtend,    & ! d(Nrm)/dt                  [count/kg/s]
      ncmtend,    & ! d(Ncm)/dt                  [count/kg/s]
      nimtend       ! d(Nim)/dt                  [count/kg/s]

    real(kind = core_rknd), dimension(gr%nz), intent(out) :: & 
      Vrr,       & ! Rain mixing ratio fall speed        [m/s]
      VNr,       & ! Rain conc. fall speed               [m/s]
      Vrs,    & ! Snow mixing ratio fall speed        [m/s]
      Vrg, & ! Graupel fall speed                  [m/s]
      Vri      ! Pristine ice mixing ratio fall speed  [m/s]

! Local Variables

    real, dimension(gr%nz) :: &
      T_in_K     ! Temperature                                     [K]

! Addition by Adam Smith, 24 April 2008
! Adding snow particle number concentration
    real, dimension(gr%nz) :: Nsm ! [count/kg]
! End of ajsmith4's addition

! Variables on the w grid
    real, dimension(1,1,gr%nz) :: & 
      w3,       & ! Vertical wind on the w grid          [m/s]
      pr3d,     & ! Pressure on w grid                   [Pa]
      qsatv3d,  & ! Saturation mr array?                 [kg/kg]
      temp3d,   & ! Temperature on w grid                [K]
      qsati3d,  & ! Saturation mr over ice array?        [kg/kg]
      th2t3d,   & ! Mean exner function on w grid        [-]
! Michael Falk, 13 Jul 2007, added these.  They are flipped versions of these variables;
! that is to say, in the following versions k=1 is the top of the domain and k=kk+1 is the
! sub-ground ghost point.
      w3_flip, & 
      pr3d_flip, & 
      qsatv3d_flip, & 
      temp3d_flip, & 
      qsati3d_flip, & 
      th2t3d_flip
! eMFc

! Variables on the m grid
    real, dimension(1,1,gr%nz-1) :: & 
      qi3,    & ! Pristine ice mixing ratio               [kg/kg]
      qr3,    & ! Rain water mixing ratio                 [kg/kg]
      qg3,    & ! Graupel mixing ratio                    [kg/kg]
      qs3,    & ! Snow water mixing ratio                 [kg/kg]
! Michael Falk, 13 Jul 2007, added these; they are flipped versions, as above.
      qi3_flip, & 
      qr3_flip, & 
      qg3_flip, & 
      qs3_flip
! eMFc

    real :: & 
      gmbov2, & 
      gmbov2g, & 
      gmbp3, & 
      gm3, & 
      gm4, & 
      gm5, & 
      gm6, & 
      gm7, & 
      gm8, & 
      gm9, & 
      ex2, & 
      ex2g, & 
      ex3, & 
      ex7, & 
      ex7g, & 
      ex4, & 
      ex4g, & 
      ex5, & 
      hlvoka

    real :: &
      hlsoka, & 
      hlvorv, & 
      hlsorv, & 
      rvochi, & 
      cpor, & 
      lfocp, & 
      lvocp, & 
      lsocp, & 
      hkaolf, & 
      sloper, & 
      slopes, & 
      slopeg, & 
      eic, & 
      sat

    real, dimension(gr%nz) :: & 
      thm, & 
      rvm

    real, dimension(1,1,gr%nz-1) :: & 
      qt3,         & ! Total water mixing ratio
      qv3,         & ! Water vapor mixing ratio
      qc3,         & ! Cloud water mixing ratio
      th3,         & ! potential temperature
      p3,          & ! perturbation exner function
      nc3,         & ! Number of cloud droplets
      nr3,         & ! Number of rain drops
      ncn3,        & ! Number of cloud nuclei
      ni3,         & ! Number of ice crystals
      exbm,        & ! Mean exner function
      rbm,         & ! Mean density
      snowslope,   & ! These variables are the individual microphysical terms.  Michael Falk
              !  Michael Falk
              ! (Name modified by Adam Smith, 24 April 2008)
      pcond,       & ! condensation/evaporation of cloud water
      psmlti,      & ! melting of cloud ice
      psacw,       & ! collection of cloud water by snow
      pgacw,       & ! ???
      piacw,       & ! ???
      pchomo,      & ! ???
      praut,       & ! autoconversion of cloud water
      pracw,       & ! collection of cloud water by rain water
      pdepi,       & ! depositional growth of cloud ice
      pint,        & ! initiation of cloud ice
      pgdep,       & ! ???
      pconv,       & ! conversion of cloud ice to snow
      psaci,       & ! collection of cloud ice by snow
      pgaci,       & ! ???
      praci,       & ! ???
      prevp,       & ! evaporation of raindrops
      psdep,       & ! depositional growth of snow
      pmltge,      & ! ???
      pgmlt          ! ???

    real, dimension(1,1,gr%nz-1) :: & 
      psmlt,       & ! melting of snow
      pgacrm,      & ! ???
      pgacwm,      & ! ???
      pracs,       & ! ???
      pgshr,       & ! ???
      pgacr,       & ! ???
      psacr,       & ! ???
      piacr,       & ! ???
      prhomo,      & ! ???
      pgacs,       & ! ???
      pmltse,      & ! ???
      pwacs,       & ! ???
! Michael Falk, 13 Jul 2007, added these; they are flipped versions of these variables
      cond_flip, & 
      p3_flip, & 
      qc3_flip, & 
      qt3_flip, & 
      qv3_flip, & 
      th3_flip, & 
      exbm_flip, & 
      rbm_flip
! eMFc

    real, dimension(1,1,gr%nz) :: & 
      fallr,   & ! Fall speed for rain mixing ratio              [m/s]
      falln,   & ! Fall speed for rain drop number conc.         [m/s]
      falli,   & ! Fall speed for pristine ice mixing ratio      [m/s]
      falls,   & ! Fall speed for snow mixing ratio              [m/s] (Michael Falk 19 Jul 2007)
      snowv,   & ! Fall speed for snow mixing ratio              [m/s]
      fallg,   & ! Fall speed for graupel mixing ratio           [m/s]
! Michael Falk, 13 Jul 2007, added these; they are flipped versions of these variables
      fallg_flip, & 
      fallr_flip, & 
      falln_flip, & 
      falli_flip, & 
      falls_flip, & 
      snowv_flip, & 
! and then added these following variables 18 Jul 2007.
! The in_cloud versions are only assigned values (and passed to and from adjtq) within clouds.
! They are passed to fallr, falln, falli, and falls, which are assigned 0 outside of cloud
! and assigned the in_cloud values within cloud.
! snowv is passed at all gridpoints, 1 to kk+1, regardless of cloud.
      fallr_in_cloud, & 
      falln_in_cloud, & 
      falli_in_cloud, & 
      falls_in_cloud, & 
      fallg_in_cloud
! eMFc

    real, dimension(1,1,gr%nz-1) :: & 
      rvc,      & ! Cloud droplet radius         [cm]
      rvr      ! Rain droplet radius          [cm]

    real, dimension(1,gr%nz-1,1) :: & 
      ary1d        ! 1d graphics parameters

    integer :: & 
      i, & 
      k,             & ! Loop control variables
      len,          & ! # of saturated points???
      icase,        & ! Which case?
      nrdamp       ! Number of points in the damping layer (COAMPS)

    real :: & 
      snzero ! ???

    integer, dimension(1) :: & 
      nkpts

    integer, dimension(gr%nz-1) :: & 
      icomp,          & !
      kcomp          !

    logical :: & 
      ldrizzle  ! is drizzle on?


    real ::  & ! Regular precision to be passed into COAMPS microphysics. Brian.
      timea, & ! Current model time                      [s]
      deltf    ! Timestep (i.e. dt_main in CLUBB)        [s]

    integer :: & 
      kk,     & ! Number of COAMPS m gridpoints in the vertical (gr%nz-1)
      kmax  ! Maximum array size (kk + ??)

!----------------------------------------------------------------------

! Begin coamps_microphys_driver code

    kk = gr%nz-1
    kmax = gr%nz

! Comment by Adam Smith, 25 March 2008
! These variables activate rain/drizzle in the COAMPS
! microphysics scheme.  Set these variables to .TRUE. if you need these
! hydrometeors in your simulations.
    ldrizzle = .false.

! Whether to produce ice in COAMPS.
! According to Jerry Schmidt of NRL,
!   if lice = .true., then we should
!   set ldrizzle = .false.
!   because collection of drizzle by ice
!   is not implemented yet.
    if ( l_ice_microphys .and. ldrizzle ) then
      write(fstderr,*) "l_ice_microphys must be false to use ldrizzle"
      error stop
    end if

!------------------------------------------
! Comment by Adam Smith, 25 March 2008
! The variable "icase" is used in the COAMPS namelist to identify
! the simulation being run.  Normally, each icase value represents
! an individual case, in order to allow the user to use specific
! model settings and forcings.
!
! Specific icase values used in COAMPS:
! icase = 75:   M-PACE: period B
! icase = 1000: Nov.11, 1999 altocu
! icase = 1001: Jun.25 multilayered altocu
! icase = 1002: CLEX-9: Oct.14, 2001 altocu
! icase = 1003: CLEX-9: Nov.02, 2001 altocu
!------------------------------------------
! End of ajsmith4's comment
!------------------------------------------

    select case ( trim( runtype ) )
! It appears that using icase = 75 for MPACE B gives spurious results (i.e. no
! ice or snow forms but thlm_mc and rtm_mc are non-zero).
! Therefore, to get results that resemble the COAMPS-LES results, setting 
! icase /= 75 seems to be the best option. 
! Perhaps some bug was introduced in the CLUBB branch of the COAMPS microphysics
! that made icase code different between COAMPS-LES and CLUBB.
! -dschanen 2 June 2011
!   case ( "mpace_a", "mpace_b" )
    case ( "mpace_a" )
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
      icase = -1 ! No special conditions. -dschanen
      nrdamp = 15

    end select

! The eic variable is something related to ice.  In regular COAMPS it is
! declared as a parameter but it is then set to another value, which is illegal.
! -dschanen
    eic      = 1.0

! Brian set regular precision variables "timea" and "deltf" to the values
! brought in from CLUBB as precision "time_precision" variables "timea_in" and
! "deltf_in".  Brian; 4/5/2008.
    timea = real( timea_in )
    deltf = real( deltf_in )

! Michael Falk changed this, November 2007.
! This is the stock COAMPS value, 2e7:
    if ( trim( runtype ) == "mpace_a" ) then
      ! and this is the value, 2e6, which works for MPACE-A:
      snzero = 2.0e6
    else ! don't cheat
      snzero = 2.0e7
    end if
! eMFc

! Set up initial fields

! Compute quantities for computing tendencies
    rvm = real(rtm - rcm)

    if( any( rvm < 0. ) ) then

        if ( clubb_at_least_debug_level( 1 ) ) then
            write(fstderr,*) 'in COAMPS (R) microphys driver: some rvm < 0'
        end if

        where ( rvm < 0. ) rvm = 0.

    end if

      thm(1:kk+1) = real(thlm(1:kk+1) & 
                  + ( Lv /( Cp * exner(1:kk+1) )* rcm(1:kk+1) ))

      ! Determine absolute temperature
      T_in_K = real(thlm2T_in_K( gr%nz, thlm, exner, rcm ))

      ! Setup COAMPS verical velocity / mass grid variables
      w3(1,1,1:kk+1) = real(wm_zm(1:kk+1))


!     do k=1, kk+1, 1
!       pr3d(1,1,k)    = zt2zm( p, k )
!       th2t3d(1,1,k)  = zt2zm( exner, k )
!       temp3d(1,1,k)  = zt2zm( thm, k ) * th2t3d(1,1,k)
!       qsatv3d(1,1,k) = sat_mixrat_liq( pr3d(1,1,k), temp3d(1,1,k) )
!       qsati3d(1,1,k) = sat_mixrat_ice( pr3d(1,1,k), temp3d(1,1,k) )
!     end do

! Setup COAMPS m (mass) grid variables that are only used from
! gridpoints 1 to kk, but are nonetheless defined from 1 to kk+1
! since they are defined as w (momentum) variables elsewhere
! within COAMPS.
! Since values 1 to kk are the only ones used, they are the only
! ones that are assigned.
! When setting up COAMPS variables, we remove the sub-ground ghost
! point (at k=1) from the CLUBB variables.  Since in CLUBB k=2 is the
! first above-ground gridpoint and in COAMPS k=1 is the first
! above-ground gridpoint, we assign the variables accordingly.
! Comments by Michael Falk, David Schanen, and Vince Larson

! The top point is undefined and unreferenced in these '3d' arrays
      pr3d(1,1,1:kk)   = real(p_in_Pa(2:kk+1))
      th2t3d(1,1,1:kk) = real(exner(2:kk+1))
      temp3d(1,1,1:kk) = T_in_K(2:kk+1)

      do k=1, kk, 1
        qsatv3d(1,1,k) = real(sat_mixrat_liq( &
                real(pr3d(1,1,k), kind = core_rknd), &
                real(temp3d(1,1,k), kind = core_rknd) ))
        qsati3d(1,1,k) = real(sat_mixrat_ice( &
                real(pr3d(1,1,k), kind = core_rknd), &
                real(temp3d(1,1,k), kind = core_rknd) ))
      end do

! Setup COAMPS m (mass) grid variables
      qt3(1,1,1:kk)  = real(rtm(2:kk+1))
      qc3(1,1,1:kk)  = real(rcm(2:kk+1))
      qr3(1,1,1:kk)  = real(rrm(2:kk+1))
      qg3(1,1,1:kk)  = real(rgm(2:kk+1))
      qs3(1,1,1:kk)  = real(rsm(2:kk+1))
      qi3(1,1,1:kk)  = real(rim(2:kk+1))
      exbm(1,1,1:kk) = real(exner(2:kk+1))
      rbm(1,1,1:kk)  = real(rho(2:kk+1))
      th3(1,1,1:kk)  = thm(2:kk+1)
      qv3(1,1,1:kk)  = rvm(2:kk+1)

      do k=1,kk
        p3(1,1,k)   = 0.0

        ! Convert from MKS units as needed
        ! Nrm, Ncm, Nccnm are in kg^-1, and need to coverted to
        ! (m^3/cm^3)*kg^-1.  They will be converted to #/cc within adjtq if
        ! ldrizzle is true, or ignored if ldrizzle is false.
        nc3(1,1,k)  = real(Ncm(k+1) / cm3_per_m3)
        nr3(1,1,k)  = real(Nrm(k+1) / cm3_per_m3)
        ncn3(1,1,k) = real(Nccnm(k+1) / cm3_per_m3)

        ! Nim is in #/m^3 within adjtq (See conice.F).
        ni3(1,1,k)  = real(Nim(k+1) * rho(k+1))

      end do

! Note that this is much simpler approx. of gamma than the ANL function
! used in the rest of the CLUBB code.
! It is used only to have a direct comparison with COAMPS-LES -dschanen
      call gamma( 3.0, gm3 )
      call gamma( 4.0, gm4 )
      call gamma( 5.0, gm5 )
      call gamma( 6.0, gm6 )
      call gamma( 7.0, gm7 ) ! Known magic number
      call gamma( 8.0, gm8 ) ! Known magic number
      call gamma( 9.0, gm9 ) ! Known magic number
      call gamma( bsnow+3.0, gmbp3 )
      call gamma( bsnow*0.5 + 2.5, gmbov2 ) ! Known magic number
      call gamma( bgrp*0.5 + 2.5, gmbov2g ) ! Known magic number

      ex2  = bsnow * 0.5 + 2.5 ! Known magic number
      ex2g = bgrp * 0.5 + 2.5 ! Known magic number
      ex3  = bsnow + 3.0
      ex7  = 0.44 * gmbov2 ! Known magic number
      ex7g = 0.31 * gmbov2g ! Known magic number
      ex4  = aprpr/visair
      ex4g = abar/visair
      ex5  = real(pi*aprpr*snzero*gmbp3/4.0)

!     Lf     = Ls - Lv   ! Latent heat of fusion (occurs in constants)
      hlvoka = real(Lv)/therco
      hkaolf = therco/real(Lf)
      hlsoka = real(Ls)/therco
      hlvorv = real(Lv)/real(Rv)
      hlsorv = real(Ls)/real(Rv)
      rvochi = real(Rv)/difvap
      cpor   = real(Cp/Rd)
      lfocp  = real(Lf/Cp)
      lvocp  = real(Lv/Cp)
      lsocp  = real(Ls/Cp)

      ary1d(1,1:kk,1) = 0. ! 1d graphics parameters
      nkpts = 0
      len = 0

! Michael Falk, 6 July 2007
! Assigning values to the "flip" variables, in which k=1 is the top of the domain

      cond_flip(1,1,1:kk) = 0.
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

      w3_flip(1,1,1:kk+1) = w3(1,1,kk+1:1:-1) ! not referenced in COAMPS microphysics

      pr3d_flip(1,1,1:kk)    = pr3d(1,1,kk:1:-1) ! top point is undefined
      qsatv3d_flip(1,1,1:kk) = qsatv3d(1,1,kk:1:-1) ! " "
      temp3d_flip(1,1,1:kk)  = temp3d(1,1,kk:1:-1) ! " "
      qsati3d_flip(1,1,1:kk) = qsati3d(1,1,kk:1:-1) ! " "
      th2t3d_flip(1,1,1:kk)  = th2t3d(1,1,kk:1:-1) ! " "



! determination of which points are "in_cloud" are which are outside of cloud
      do k=nrdamp,kk
!-------------------------------------------------------------------------------
! determine if point is saturated as in COAMPS
!-------------------------------------------------------------------------------
        if ( .not. l_ice_microphys ) then
          sat = qv3_flip(1,1,k)/qsatv3d_flip(1,1,k)-1.0
        else
          if ( temp3d_flip(1,1,k) >= T_freeze_K ) then
            sat = qv3_flip(1,1,k)/qsatv3d_flip(1,1,k)-1.0
          else
            sat = qv3_flip(1,1,k)/qsati3d_flip(1,1,k)-1.0
          endif
        endif

        if (sat > 0.0 .or. & 
            qc3_flip(1,1,k) >= pcut .or. & 
            qr3_flip(1,1,k) >= pcut .or. & 
            qs3_flip(1,1,k) >= pcut .or. & 
            qi3_flip(1,1,k) >= pcut .or. & 
            qg3_flip(1,1,k) >= pcut & 
           ) & 
            then
          nkpts = nkpts+1
          len   = len+1
          kcomp(nkpts) = k
        end if
      end do

      do i=1,ipts
        icomp(i) = 1
      end do

      sloper = real(pi * rholiq * rnzero * 1.0e-8)
      slopes = real(pi * rhosno * snzero * 1.0e-8)
      slopeg = real(pi * rhogrp * gnzero * 1.0e-8)

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
! dschanen added this to avoid uninitialized memory 24 Jul 2008
      snowslope = 0.
      rvr       = 0.
      rvc       = 0.

! Call the actual COAMPS microphysics scheme
      if ( len > 0 ) call adjtq &
             (cond_flip(1,1,1:kk),p3_flip(1,1,1:kk),qc3_flip(1,1,1:kk) &
             ,qi3_flip(1,1,1:kk),qr3_flip(1,1,1:kk),qg3_flip(1,1,1:kk) &
             ,qs3_flip(1,1,1:kk),qv3_flip(1,1,1:kk),th3_flip(1,1,1:kk) &
             ,w3_flip(1,1,1:kk+1),pr3d_flip(1,1,1:kk+1) &
             ,qsatv3d_flip(1,1,1:kk+1),temp3d_flip(1,1,1:kk+1) &
             ,qsati3d_flip(1,1,1:kk+1),th2t3d_flip(1,1,1:kk+1),wtm &
             ,exbm_flip(1,1,1:kk),rbm_flip(1,1,1:kk) &
!     3       ,nc3,nr3,ncn3,ni3,cp,deltf,Lf,Ls,Lv 
             ,nc3(1,1,kk:1:-1),nr3(1,1,kk:1:-1),ncn3(1,1,kk:1:-1) &
             ,ni3(1,1,kk:1:-1),real(Cp),deltf,real(Lf),real(Ls),real(Lv) &
             ,pcut,real(p0),real(Rd),real(Rv),sloper,slopes,slopeg,timea,l_ice_microphys &
             ,nne,kk,i1d,j1d,ary1d,i1dflg,n1d,maxpt1d,maxvr1d &
             ,kmax,nrdamp,ipts,nkpts,icomp,kcomp,j &
             ,xland,aa0,aa1,aa2,aa3,abar,apr,aprpr,bsnow &
             ,cbeta,cnzero,cimass,cpor,cw,difvap,erc,esi,eic &
             ,eri,egc,esc,esr,egi,egr,egs,mw,real(pi),praut1,praut2 &
             ,rholiq,rhosno,rnzero,snzero,gnzero,therco,tice &
             ,tvr1,tvr2,tvr3,tvr4,tzero,visair,gm3,gm4,gm5,gm6 &
             ,gm7,gm8,gm9,gmbp3,gmbov2,gmbov2g,bgrp,ex1,ex2 &
             ,ex2g,ex3,hlvoka,hkaolf,hlsoka,hlvorv,hlsorv &
             ,rvochi,lfocp,lvocp,lsocp,ex7,ex7g,ex4,ex4g,ex5 &
             ,ldrizzle,l_graupel,icon,icond,len,icase &
             ,snowv_flip(:,:,1:kk) &
             ,snowslope(1,1,kk:1:-1),pcond(1,1,kk:1:-1) &
             ,psmlti(1,1,kk:1:-1),psacw(1,1,kk:1:-1),pgacw(1,1,kk:1:-1) &
             ,piacw(1,1,kk:1:-1),pchomo(1,1,kk:1:-1),praut(1,1,kk:1:-1) &
             ,pracw(1,1,kk:1:-1),pdepi(1,1,kk:1:-1),pint(1,1,kk:1:-1) &
             ,pgdep(1,1,kk:1:-1),pconv(1,1,kk:1:-1),psaci(1,1,kk:1:-1) &
             ,pgaci(1,1,kk:1:-1),praci(1,1,kk:1:-1),prevp(1,1,kk:1:-1) &
             ,psdep(1,1,kk:1:-1),pmltge(1,1,kk:1:-1),pgmlt(1,1,kk:1:-1) &
             ,psmlt(1,1,kk:1:-1),pgacrm(1,1,kk:1:-1),pgacwm(1,1,kk:1:-1) &
             ,pracs(1,1,kk:1:-1),pgshr(1,1,kk:1:-1),pgacr(1,1,kk:1:-1) &
             ,psacr(1,1,kk:1:-1),piacr(1,1,kk:1:-1),prhomo(1,1,kk:1:-1) &
             ,pgacs(1,1,kk:1:-1),pmltse(1,1,kk:1:-1),pwacs(1,1,kk:1:-1) &
             ,falli_in_cloud(:,:,1:nkpts(1)) &
             ,fallg_in_cloud(:,:,1:nkpts(1)) &
             ,fallr_in_cloud(:,:,1:nkpts(1)) &
             ,falln_in_cloud(:,:,1:nkpts(1)) &
             ,falls_in_cloud(:,:,1:nkpts(1)) &
             ,rvc,rvr)

! reassigning flipped versions of variables to normal versions
      cond(1,1,1:kk) = real(cond_flip(1,1,kk:1:-1), kind = core_rknd)
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

! Transfer back to CLUBB arrays
      do k=1, kk, 1
        ! Convert to MKS as needed
        ! ncn3 is in (m^3/cm^3)*kg^-1, and needs to be converted to kg^-1.
        Nccnm(k+1) = real(ncn3(1,1,k), kind = core_rknd) * cm3_per_m3
      end do ! k=1..kk+1

!-------------------------------------------------
! Addition by Adam Smith, 24 April 2008
! Adding snow particle number concentration
! Values of snowslope < 1.0 lead to excessive and
! unrealistic Nsm outside of the snow region.
! The "if" statement prevents these results.
!-------------------------------------------------
      do k = 1, kk, 1
        if (snowslope(1,1,k) < 2.0) then
          snowslope(1,1,k) = 0.
          Nsm(k+1) = 0.0
        else
          Nsm(k+1) = snzero / snowslope(1,1,k)
        end if
        ! Convert to #/kg for comparison to Morrison -dschanen 12 Nov 2009
        Nsm(k+1) = Nsm(k+1) / real(rho(k+1))
      end do

!-------------------------------------------------
! End of ajsmith4's addition
!-------------------------------------------------

! Linear extrapolation for the ghost point of fall speeds
      fallr(1,1,1) = .5 * ( fallr(1,1,2) + fallr(1,1,3) )
      falln(1,1,1) = .5 * ( falln(1,1,2) + falln(1,1,3) )
      snowv(1,1,1) = .5 * ( snowv(1,1,2) + snowv(1,1,3) )
      falli(1,1,1) = .5 * ( falli(1,1,2) + falli(1,1,3) )
      fallg(1,1,1) = .5 * ( fallg(1,1,2) + fallg(1,1,3) )
      falls(1,1,1) = .5 * ( falls(1,1,2) + falls(1,1,3) )

      Vrr      = zt2zm( gr, real(fallr(1,1,:), kind = core_rknd) )
      VNr      = zt2zm( gr, real(falln(1,1,:), kind = core_rknd) )
      Vrs    = zt2zm( gr, real(snowv(1,1,:), kind = core_rknd) )
      Vri     = zt2zm( gr, real(falli(1,1,:), kind = core_rknd) )
      Vrg = zt2zm( gr, real(fallg(1,1,:), kind = core_rknd) )

! Compute tendencies
      do k=1, kk, 1
        rrtend(k+1)    = ( real(qr3(1,1,k), kind = core_rknd) - rrm(k+1) ) &
                           / real(deltf, kind = core_rknd)
        rgtend(k+1)    = ( real(qg3(1,1,k), kind = core_rknd) - rgm(k+1) )&
                           / real(deltf, kind = core_rknd)
        ritend(k+1)    = ( real(qi3(1,1,k), kind = core_rknd) - rim(k+1) ) &
                           / real(deltf, kind = core_rknd)
        nrmtend(k+1)   = ( (real(nr3(1,1,k), kind = core_rknd)*cm3_per_m3) - Nrm(k+1) ) &
                           / real(deltf, kind = core_rknd) ! Conversion factor
        rstend(k+1) = ( real(qs3(1,1,k), kind = core_rknd) - rsm(k+1) ) &
                           / real(deltf, kind = core_rknd)
        ! nc3 is in (m^3/cm^3)*kg^-1, and needs to be converted to kg^-1.
        ncmtend(k+1)   = ( ( real( nc3(1,1,k), kind = core_rknd ) * cm3_per_m3 ) &
                             - Ncm(k+1) ) &
                           / real( deltf, kind = core_rknd )
        ! Convert ice number concentration to #/kg
        nimtend(k+1)   = ( ( real( ni3(1,1,k), kind = core_rknd ) / rho(k+1) ) &
                             - Nim(k+1) ) &
                           / real( deltf, kind = core_rknd )
        rvm_mc(k+1)    = real(((qv3(1,1,k) - rvm(k+1)) / deltf), kind = core_rknd)
        rcm_mc(k+1)    = (real(qc3(1,1,k), kind = core_rknd) - rcm(k+1)) &
                           / real(deltf, kind = core_rknd)
        thlm_mc(k+1)  & 
        = real(( ( th3(1,1,k) - (real(Lv) / (real(Cp) * exbm(1,1,k)) * qc3(1,1,k) ) ) & 
            - real(thlm(k+1)) ) / deltf, kind = core_rknd)
      end do ! k=1..kk

      rrtend(1)    = 0.0_core_rknd
      rgtend(1)    = 0.0_core_rknd
      ritend(1)    = 0.0_core_rknd
      nrmtend(1)   = 0.0_core_rknd
      rstend(1) = 0.0_core_rknd
      ncmtend(1)   = 0.0_core_rknd
      nimtend(1)   = 0.0_core_rknd
      rcm_mc(1)    = 0.0_core_rknd
      rvm_mc(1)    = 0.0_core_rknd
      thlm_mc(1)   = 0.0_core_rknd

      if ( stats_metadata%l_stats_samp ) then
        ! Mean volume radius of rain and cloud droplets
        do k=2,kk+1
          call stat_update_var_pt( stats_metadata%im_vol_rad_rain, k, &
               real(rvr(1,1,k-1) / 100.0, kind = core_rknd), stats_zt )
          call stat_update_var_pt( stats_metadata%im_vol_rad_cloud, k, &
               real(rvc(1,1,k-1) / 100.0, kind = core_rknd), stats_zt )
        end do


! Addition by Adam Smith, 24 April 2008
! Adding calculation for snow particle number concentration
        do k = 2,kk,1
          call stat_update_var_pt( stats_metadata%isnowslope, k, real(snowslope(1,1,k), &
            kind = core_rknd), stats_zt)
        end do

! Addition by Adam Smith, 25 April 2008
! Adding calculation for snow particle number concentration
        do k=2, kk+1
          call stat_update_var_pt( stats_metadata%iNsm, k, real(Nsm(k), kind = core_rknd) ,stats_zt )
        end do

      end if



      return
    end subroutine coamps_microphys_driver

#else /* COAMPS_MICRO not defined */
   subroutine coamps_microphys_driver( )

     implicit none

     error stop "Not compiled with COAMPS Microphysics"

   end subroutine coamps_microphys_driver
#endif

  end module coamps_microphys_driver_module
