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
  subroutine coamps_microphys_driver( gr, ngrdcol, & ! In
                                      runtype, timea_in, deltf_in, & ! In
                                      rtm, wm_zm, p_in_Pa, exner, rho, & ! In
                                      thlm, T_in_K, & ! In
                                      rim, rrm, & ! In
                                      rgm, rsm, & ! In
                                      rcm, Ncm, & ! In
                                      Nrm, Nim, & ! In
                                      saturation_formula, & ! In
                                      stats, Nccnm, & ! InOut
                                      Vrs, Vri, & ! Out
                                      Vrr, VNr, & ! Out
                                      Vrg, & ! Out
                                      ritend, rrtend, & ! Out
                                      rgtend, rstend, & ! Out
                                      nrmtend, ncmtend, nimtend, & ! Out
                                      rvm_mc, rcm_mc, thlm_mc ) ! Out

!    Description:
!      Subroutine to compute ice, as it is done in NRL's COAMPS
!      mesoscale model, using adjtq.F. Process all CLUBB columns here while
!      retaining the single-column interface to the imported COAMPS routines.
!      T_in_K is the batch temperature already computed by microphys_driver.

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
    use saturation, only: sat_mixrat_liq_api, sat_mixrat_ice ! Procedure(s)
    use clubb_precision, only: time_precision, core_rknd ! Variable(s)

    use grid_class, only: grid ! Type

    use stats_netcdf, only: &
        stats_type, &
        stats_update

    use error_code, only: &
        clubb_at_least_debug_level_api  ! Procedure

    use parameters_microphys, only: l_graupel, l_ice_microphys ! Variable(s)

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

    !--------------------------- Input Variables ---------------------------
    type(grid), intent(in) :: gr ! CLUBB vertical grids

    integer, intent(in) :: ngrdcol ! Number of grid columns [-]

    character(len=*), intent(in) :: runtype ! Case name

    real( kind = time_precision ), intent(in) :: timea_in ! Current model time [s]
    real( kind = core_rknd ), intent(in) :: deltf_in ! Main timestep [s]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), intent(in) :: &
      rtm ! Total water mixing ratio [kg/kg]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzm), intent(in) :: &
      wm_zm ! Vertical wind on momentum levels [m/s]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), intent(in) :: &
      p_in_Pa, & ! Pressure [Pa]
      exner,   & ! Exner function [-]
      rho,     & ! Density [kg/m^3]
      thlm,    & ! Liquid potential temperature [K]
      T_in_K,  & ! Absolute temperature computed for the full batch [K]
      rim,     & ! Ice water mixing ratio [kg/kg]
      rrm,     & ! Rain water mixing ratio [kg/kg]
      rgm,     & ! Graupel water mixing ratio [kg/kg]
      rsm,     & ! Snow water mixing ratio [kg/kg]
      rcm,     & ! Cloud water mixing ratio [kg/kg]
      Ncm,     & ! Number of cloud droplets [count/kg]
      Nrm,     & ! Number of rain drops [count/kg]
      Nim        ! Number of ice crystals [count/kg]

    integer, intent(in) :: saturation_formula ! Saturation formula selector [-]

    !------------------------ Input/Output Variables ------------------------
    type(stats_type), intent(inout) :: stats ! Statistics buffers and sample controls

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), intent(inout) :: &
      Nccnm ! Number of cloud nuclei [count/kg]

    !--------------------------- Output Variables ---------------------------
    ! All sedimentation velocities are on thermodynamic levels.
    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), intent(out) :: &
      Vrs,     & ! Snow mixing ratio fall speed [m/s]
      Vri,     & ! Pristine ice mixing ratio fall speed [m/s]
      Vrr,     & ! Rain mixing ratio fall speed [m/s]
      VNr,     & ! Rain concentration fall speed [m/s]
      Vrg,     & ! Graupel mixing ratio fall speed [m/s]
      ritend,  & ! Ice water tendency [kg/kg/s]
      rrtend,  & ! Rain water tendency [kg/kg/s]
      rgtend,  & ! Graupel water tendency [kg/kg/s]
      rstend,  & ! Snow water tendency [kg/kg/s]
      nrmtend, & ! Rain number tendency [count/kg/s]
      ncmtend, & ! Cloud droplet number tendency [count/kg/s]
      nimtend, & ! Ice number tendency [count/kg/s]
      rvm_mc,  & ! Vapor water tendency [kg/kg/s]
      rcm_mc,  & ! Cloud water tendency [kg/kg/s]
      thlm_mc    ! Liquid potential temperature tendency [K/s]

    !--------------------------- Local Variables ---------------------------
    integer :: icol, k, kf ! Column, CLUBB level and reversed COAMPS level [-]
    integer :: kk, kmax, icase, nrdamp ! Scheme level counts, case and damping controls [-]
    integer :: nkpts(ngrdcol), len(ngrdcol) ! Saturated point counts in each column [-]
    integer :: kcomp(ngrdcol,gr%nzt) ! Packed saturated level indices in each column [-]
    integer :: icomp(ipts) ! Native horizontal point index, always one [-]
    logical :: ldrizzle ! Enable drizzle in the imported scheme [-]

    real :: timea, deltf ! Scheme time and timestep [s]
    real :: snzero ! Intercept of the snow size distribution [1/m^4]
    real :: sat ! Supersaturation at one point [-]
    real :: gmbov2, gmbov2g, gmbp3, gm3, gm4, gm5, gm6, gm7, gm8, gm9 ! Gamma factors [-]
    real :: ex2, ex2g, ex3, ex7, ex7g, ex4, ex4g, ex5 ! Scheme coefficients [units vary]
    real :: hlvoka, hlsoka, hlvorv, hlsorv, rvochi, cpor ! Thermodynamic coefficients [units vary]
    real :: lfocp, lvocp, lsocp, hkaolf, sloper, slopes, slopeg, eic ! Scheme coefficients [units vary]

    ! Scheme buffers use default real, preserving the original conversion points.
    ! The middle singleton dimension belongs to COAMPS; columns remain contiguous.
    real, dimension(ngrdcol,1,gr%nzm) :: &
      w3_flip, & ! Reversed vertical velocity on momentum levels [m/s]
      pr3d_flip, & ! Reversed pressure on thermodynamic levels; unused top point [Pa]
      qsatv3d_flip, & ! Reversed liquid saturation mixing ratio; unused top point [kg/kg]
      temp3d_flip, & ! Reversed temperature; unused top point [K]
      qsati3d_flip, & ! Reversed ice saturation mixing ratio; unused top point [kg/kg]
      th2t3d_flip ! Reversed Exner function; unused top point [-]

    real, dimension(ngrdcol,1,gr%nzt) :: &
      qi3_flip, & ! Reversed ice water mixing ratio [kg/kg]
      qr3_flip, & ! Reversed rain water mixing ratio [kg/kg]
      qg3_flip, & ! Reversed graupel mixing ratio [kg/kg]
      qs3_flip, & ! Reversed snow water mixing ratio [kg/kg]
      qc3_flip, & ! Reversed cloud water mixing ratio [kg/kg]
      qv3_flip, & ! Reversed vapor mixing ratio [kg/kg]
      th3_flip, & ! Reversed potential temperature [K]
      p3_flip, & ! Reversed perturbation Exner function [-]
      exbm_flip, & ! Reversed mean Exner function [-]
      rbm_flip, & ! Reversed mean density [kg/m^3]
      cond_flip, & ! Reversed condensation diagnostic [kg/kg]
      nc3, & ! Cloud droplet number in native concentration units [(m^3/cm^3)/kg]
      nr3, & ! Rain droplet number in native concentration units [(m^3/cm^3)/kg]
      ncn3, & ! Cloud nuclei in native concentration units [(m^3/cm^3)/kg]
      ni3, & ! Ice crystal number per volume [count/m^3]
      snowslope, & ! Snow distribution slope in CLUBB vertical order [1/m]
      rvc, & ! Cloud droplet mean volume radius in CLUBB vertical order [cm]
      rvr, & ! Rain droplet mean volume radius in CLUBB vertical order [cm]
      falli_flip, & ! Reversed ice fall speed [m/s]
      falls_flip, & ! Reversed snow fall speed [m/s]
      fallg_flip, & ! Reversed graupel fall speed [m/s]
      fallr_flip, & ! Reversed rain fall speed [m/s]
      falln_flip, & ! Reversed rain number fall speed [m/s]
      falli_in_cloud, & ! Packed in-cloud ice fall speed [m/s]
      falls_in_cloud, & ! Packed in-cloud snow fall speed [m/s]
      fallg_in_cloud, & ! Packed in-cloud graupel fall speed [m/s]
      fallr_in_cloud, & ! Packed in-cloud rain fall speed [m/s]
      falln_in_cloud ! Packed in-cloud rain number fall speed [m/s]

    ! These native diagnostics are required by adjtq but are not consumed by CLUBB.
    real, dimension(1,1,gr%nzt) :: &
      pcond, & ! Unused native scheme diagnostic [units vary]
      psmlti, & ! Unused native scheme diagnostic [units vary]
      psacw, & ! Unused native scheme diagnostic [units vary]
      pgacw, & ! Unused native scheme diagnostic [units vary]
      piacw, & ! Unused native scheme diagnostic [units vary]
      pchomo, & ! Unused native scheme diagnostic [units vary]
      praut, & ! Unused native scheme diagnostic [units vary]
      pracw, & ! Unused native scheme diagnostic [units vary]
      pdepi, & ! Unused native scheme diagnostic [units vary]
      pint, & ! Unused native scheme diagnostic [units vary]
      pgdep, & ! Unused native scheme diagnostic [units vary]
      pconv, & ! Unused native scheme diagnostic [units vary]
      psaci, & ! Unused native scheme diagnostic [units vary]
      pgaci, & ! Unused native scheme diagnostic [units vary]
      praci, & ! Unused native scheme diagnostic [units vary]
      prevp, & ! Unused native scheme diagnostic [units vary]
      psdep, & ! Unused native scheme diagnostic [units vary]
      pmltge, & ! Unused native scheme diagnostic [units vary]
      pgmlt, & ! Unused native scheme diagnostic [units vary]
      psmlt, & ! Unused native scheme diagnostic [units vary]
      pgacrm, & ! Unused native scheme diagnostic [units vary]
      pgacwm, & ! Unused native scheme diagnostic [units vary]
      pracs, & ! Unused native scheme diagnostic [units vary]
      pgshr, & ! Unused native scheme diagnostic [units vary]
      pgacr, & ! Unused native scheme diagnostic [units vary]
      psacr, & ! Unused native scheme diagnostic [units vary]
      piacr, & ! Unused native scheme diagnostic [units vary]
      prhomo, & ! Unused native scheme diagnostic [units vary]
      pgacs, & ! Unused native scheme diagnostic [units vary]
      pmltse, & ! Unused native scheme diagnostic [units vary]
      pwacs, & ! Unused native scheme diagnostic [units vary]
      snowv_flip ! Unused native scheme diagnostic [units vary]

    real :: ary1d(1,gr%nzt,1) ! Disabled native graphics buffer [units vary]
    real :: rvm(ngrdcol,gr%nzt) ! Initial water vapor mixing ratio [kg/kg]
    real :: Nsm(ngrdcol,gr%nzt) ! Diagnostic snow number concentration [count/kg]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt) :: &
      mvrr_stats, & ! Rain mean volume radius [m]
      cloud_radius_stats, & ! Cloud droplet mean volume radius [m]
      snowslope_stats, & ! Snow distribution slope [1/m]
      Nsm_stats ! Diagnostic snow number concentration [count/kg]

    !--------------------------- Begin Code ---------------------------
    kk = gr%nzt
    kmax = gr%nzm

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


    sloper = real(pi * rholiq * rnzero * 1.0e-8)
    slopes = real(pi * rhosno * snzero * 1.0e-8)
    slopeg = real(pi * rhogrp * gnzero * 1.0e-8)
    icomp = 1
    nkpts = 0
    len = 0
    snowslope = 0.
    rvr = 0.
    rvc = 0.

    ! Convert the full batch to scheme precision and its top-down vertical order.
    do k = 1, kk
      kf = kk + 1 - k
      do icol = 1, ngrdcol
        rvm(icol,k) = real(rtm(icol,k) - rcm(icol,k))
        th3_flip(icol,1,kf) = real(thlm(icol,k) &
          + ( Lv / ( Cp * exner(icol,k) ) * rcm(icol,k) ))
        pr3d_flip(icol,1,kf) = real(p_in_Pa(icol,k))
        th2t3d_flip(icol,1,kf) = real(exner(icol,k))
        temp3d_flip(icol,1,kf) = real(T_in_K(icol,k))
        qc3_flip(icol,1,kf) = real(rcm(icol,k))
        qr3_flip(icol,1,kf) = real(rrm(icol,k))
        qg3_flip(icol,1,kf) = real(rgm(icol,k))
        qs3_flip(icol,1,kf) = real(rsm(icol,k))
        qi3_flip(icol,1,kf) = real(rim(icol,k))
        exbm_flip(icol,1,kf) = real(exner(icol,k))
        rbm_flip(icol,1,kf) = real(rho(icol,k))
        p3_flip(icol,1,kf) = 0.
        cond_flip(icol,1,kf) = 0.

        ! Concentrations retain their original bottom-up order until the call.
        nc3(icol,1,k) = real(Ncm(icol,k) / cm3_per_m3)
        nr3(icol,1,k) = real(Nrm(icol,k) / cm3_per_m3)
        ncn3(icol,1,k) = real(Nccnm(icol,k) / cm3_per_m3)
        ni3(icol,1,k) = real(Nim(icol,k) * rho(icol,k))
      end do
    end do

    do icol = 1, ngrdcol
      if ( any( rvm(icol,:) < 0. ) .and. clubb_at_least_debug_level_api( 1 ) ) then
        write(fstderr,*) 'in COAMPS (R) microphys driver: some rvm < 0'
      end if
    end do

    do k = 1, kk
      kf = kk + 1 - k
      do icol = 1, ngrdcol
        if ( rvm(icol,k) < 0. ) rvm(icol,k) = 0.
        qv3_flip(icol,1,kf) = rvm(icol,k)
        qsatv3d_flip(icol,1,kf) = real(sat_mixrat_liq_api( &
          real(pr3d_flip(icol,1,kf), kind = core_rknd), &
          real(temp3d_flip(icol,1,kf), kind = core_rknd), saturation_formula ))
        qsati3d_flip(icol,1,kf) = real(sat_mixrat_ice( &
          real(pr3d_flip(icol,1,kf), kind = core_rknd), &
          real(temp3d_flip(icol,1,kf), kind = core_rknd), saturation_formula ))

        falli_flip(icol,1,kf) = 0.
        falls_flip(icol,1,kf) = 0.
        fallg_flip(icol,1,kf) = 0.
        fallr_flip(icol,1,kf) = 0.
        falln_flip(icol,1,kf) = 0.
        falli_in_cloud(icol,1,kf) = 0.
        falls_in_cloud(icol,1,kf) = 0.
        fallg_in_cloud(icol,1,kf) = 0.
        fallr_in_cloud(icol,1,kf) = 0.
        falln_in_cloud(icol,1,kf) = 0.
      end do
    end do

    do k = 1, kk+1
      do icol = 1, ngrdcol
        w3_flip(icol,1,kk+2-k) = real(wm_zm(icol,k))
      end do
    end do

    ! Pack saturated points in the same top-down order within each column.
    do k = nrdamp, kk
      do icol = 1, ngrdcol
        if ( .not. l_ice_microphys ) then
          sat = qv3_flip(icol,1,k)/qsatv3d_flip(icol,1,k)-1.0
        else
          if ( temp3d_flip(icol,1,k) >= T_freeze_K ) then
            sat = qv3_flip(icol,1,k)/qsatv3d_flip(icol,1,k)-1.0
          else
            sat = qv3_flip(icol,1,k)/qsati3d_flip(icol,1,k)-1.0
          end if
        end if

        if ( sat > 0.0 .or. qc3_flip(icol,1,k) >= pcut .or. &
             qr3_flip(icol,1,k) >= pcut .or. qs3_flip(icol,1,k) >= pcut .or. &
             qi3_flip(icol,1,k) >= pcut .or. qg3_flip(icol,1,k) >= pcut ) then
          nkpts(icol) = nkpts(icol) + 1
          len(icol) = len(icol) + 1
          kcomp(icol,nkpts(icol)) = k
        end if
      end do
    end do

    do icol = 1, ngrdcol
      ! The imported ice-collection routine modifies eic; reset it per column.
      eic = 1.0
      ary1d = 0.

      ! Only the imported scalar scheme runs one whole column at a time.
      if ( len(icol) > 0 ) call adjtq &
             (cond_flip(icol,1,1:kk),p3_flip(icol,1,1:kk),qc3_flip(icol,1,1:kk) &
             ,qi3_flip(icol,1,1:kk),qr3_flip(icol,1,1:kk),qg3_flip(icol,1,1:kk) &
             ,qs3_flip(icol,1,1:kk),qv3_flip(icol,1,1:kk),th3_flip(icol,1,1:kk) &
             ,w3_flip(icol,1,1:kk+1),pr3d_flip(icol,1,1:kk+1) &
             ,qsatv3d_flip(icol,1,1:kk+1),temp3d_flip(icol,1,1:kk+1) &
             ,qsati3d_flip(icol,1,1:kk+1),th2t3d_flip(icol,1,1:kk+1),wtm &
             ,exbm_flip(icol,1,1:kk),rbm_flip(icol,1,1:kk) &
             ,nc3(icol,1,kk:1:-1),nr3(icol,1,kk:1:-1),ncn3(icol,1,kk:1:-1) &
             ,ni3(icol,1,kk:1:-1),real(Cp),deltf,real(Lf),real(Ls),real(Lv) &
             ,pcut,real(p0),real(Rd),real(Rv),sloper,slopes,slopeg,timea,l_ice_microphys &
             ,nne,kk,i1d,j1d,ary1d,i1dflg,n1d,maxpt1d,maxvr1d &
             ,kmax,nrdamp,ipts,nkpts(icol:icol),icomp,kcomp(icol,:),j &
             ,xland,aa0,aa1,aa2,aa3,abar,apr,aprpr,bsnow &
             ,cbeta,cnzero,cimass,cpor,cw,difvap,erc,esi,eic &
             ,eri,egc,esc,esr,egi,egr,egs,mw,real(pi),praut1,praut2 &
             ,rholiq,rhosno,rnzero,snzero,gnzero,therco,tice &
             ,tvr1,tvr2,tvr3,tvr4,tzero,visair,gm3,gm4,gm5,gm6 &
             ,gm7,gm8,gm9,gmbp3,gmbov2,gmbov2g,bgrp,ex1,ex2 &
             ,ex2g,ex3,hlvoka,hkaolf,hlsoka,hlvorv,hlsorv &
             ,rvochi,lfocp,lvocp,lsocp,ex7,ex7g,ex4,ex4g,ex5 &
             ,ldrizzle,l_graupel,icon,icond,len(icol),icase &
             ,snowv_flip(:,:,1:kk) &
             ,snowslope(icol,1,kk:1:-1),pcond(1,1,kk:1:-1) &
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
             ,falli_in_cloud(icol:icol,:,1:nkpts(icol)) &
             ,fallg_in_cloud(icol:icol,:,1:nkpts(icol)) &
             ,fallr_in_cloud(icol:icol,:,1:nkpts(icol)) &
             ,falln_in_cloud(icol:icol,:,1:nkpts(icol)) &
             ,falls_in_cloud(icol:icol,:,1:nkpts(icol)) &
             ,rvc(icol:icol,:,:),rvr(icol:icol,:,:))

    end do

    ! Scatter packed fall speeds back to the reversed full grid.
    do k = 1, kk
      do icol = 1, ngrdcol
        if ( k > nkpts(icol) ) cycle
        falli_flip(icol,1,kcomp(icol,k)) = falli_in_cloud(icol,1,k)
        falls_flip(icol,1,kcomp(icol,k)) = falls_in_cloud(icol,1,k)
        fallg_flip(icol,1,kcomp(icol,k)) = fallg_in_cloud(icol,1,k)
        fallr_flip(icol,1,kcomp(icol,k)) = fallr_in_cloud(icol,1,k)
        falln_flip(icol,1,kcomp(icol,k)) = falln_in_cloud(icol,1,k)
      end do
    end do

    ! Keep the original order of cloud/rain/aerosol concentration clipping.
    do k = 1, kk
      do icol = 1, ngrdcol
        if ( nr3(icol,1,k) < 0. ) then
          ncn3(icol,1,k) = ncn3(icol,1,k) + nr3(icol,1,k)
          nr3(icol,1,k) = 0.
        end if
        if ( nc3(icol,1,k) < 0. ) then
          ncn3(icol,1,k) = ncn3(icol,1,k) + nc3(icol,1,k)
          nc3(icol,1,k) = 0.
        end if
        if ( ncn3(icol,1,k) < 0. ) ncn3(icol,1,k) = 0.

        Nccnm(icol,k) = real(ncn3(icol,1,k), kind = core_rknd) * cm3_per_m3
        if ( snowslope(icol,1,k) < 2.0 ) then
          snowslope(icol,1,k) = 0.
          Nsm(icol,k) = 0.0
        else
          Nsm(icol,k) = snzero / snowslope(icol,1,k)
        end if
        Nsm(icol,k) = Nsm(icol,k) / real(rho(icol,k))
      end do
    end do

    ! Compute CLUBB tendencies after converting the same rounded scheme outputs.
    do k = 1, kk
      kf = kk + 1 - k
      do icol = 1, ngrdcol
        Vrr(icol,k) = real(fallr_flip(icol,1,kf), kind = core_rknd)
        VNr(icol,k) = real(falln_flip(icol,1,kf), kind = core_rknd)
        Vrs(icol,k) = real(falls_flip(icol,1,kf), kind = core_rknd)
        Vri(icol,k) = real(falli_flip(icol,1,kf), kind = core_rknd)
        Vrg(icol,k) = real(fallg_flip(icol,1,kf), kind = core_rknd)
        rrtend(icol,k) = ( real(qr3_flip(icol,1,kf), kind = core_rknd) - rrm(icol,k) ) &
          / real(deltf, kind = core_rknd)
        rgtend(icol,k) = ( real(qg3_flip(icol,1,kf), kind = core_rknd) - rgm(icol,k) ) &
          / real(deltf, kind = core_rknd)
        ritend(icol,k) = ( real(qi3_flip(icol,1,kf), kind = core_rknd) - rim(icol,k) ) &
          / real(deltf, kind = core_rknd)
        nrmtend(icol,k) = ( (real(nr3(icol,1,k), kind = core_rknd)*cm3_per_m3) - Nrm(icol,k) ) &
          / real(deltf, kind = core_rknd)
        rstend(icol,k) = ( real(qs3_flip(icol,1,kf), kind = core_rknd) - rsm(icol,k) ) &
          / real(deltf, kind = core_rknd)
        ncmtend(icol,k) = ( (real(nc3(icol,1,k), kind = core_rknd)*cm3_per_m3) - Ncm(icol,k) ) &
          / real(deltf, kind = core_rknd)
        nimtend(icol,k) = ( (real(ni3(icol,1,k), kind = core_rknd)/rho(icol,k)) - Nim(icol,k) ) &
          / real(deltf, kind = core_rknd)
        rvm_mc(icol,k) = real((qv3_flip(icol,1,kf) - rvm(icol,k)) / deltf, kind = core_rknd)
        rcm_mc(icol,k) = ( real(qc3_flip(icol,1,kf), kind = core_rknd) - rcm(icol,k) ) &
          / real(deltf, kind = core_rknd)
        thlm_mc(icol,k) = real(( ( th3_flip(icol,1,kf) &
          - (real(Lv) / (real(Cp)*exbm_flip(icol,1,kf)) * qc3_flip(icol,1,kf)) ) &
          - real(thlm(icol,k)) ) / deltf, kind = core_rknd)
      end do
    end do

    if ( stats%l_sample ) then
      do k = 1, kk
        do icol = 1, ngrdcol
          mvrr_stats(icol,k) = real(rvr(icol,1,k) / 100.0, kind = core_rknd)
          cloud_radius_stats(icol,k) = real(rvc(icol,1,k) / 100.0, kind = core_rknd)
          snowslope_stats(icol,k) = real(snowslope(icol,1,k), kind = core_rknd)
          Nsm_stats(icol,k) = real(Nsm(icol,k), kind = core_rknd)
        end do
      end do
    end if

    if ( stats%l_sample ) then
      call stats_update( "mvrr", mvrr_stats, & ! In
                         stats ) ! InOut
      call stats_update( "m_vol_rad_cloud", cloud_radius_stats, & ! In
                         stats ) ! InOut
      call stats_update( "snowslope", snowslope_stats, & ! In
                         stats ) ! InOut
      call stats_update( "Nsm", Nsm_stats, & ! In
                         stats ) ! InOut
    end if

    return

  end subroutine coamps_microphys_driver

#else /* COAMPS_MICRO not defined */
  subroutine coamps_microphys_driver( gr, ngrdcol, & ! In
                                      runtype, timea_in, deltf_in, & ! In
                                      rtm, wm_zm, p_in_Pa, exner, rho, & ! In
                                      thlm, T_in_K, & ! In
                                      rim, rrm, & ! In
                                      rgm, rsm, & ! In
                                      rcm, Ncm, & ! In
                                      Nrm, Nim, & ! In
                                      saturation_formula, & ! In
                                      stats, Nccnm, & ! InOut
                                      Vrs, Vri, & ! Out
                                      Vrr, VNr, & ! Out
                                      Vrg, & ! Out
                                      ritend, rrtend, & ! Out
                                      rgtend, rstend, & ! Out
                                      nrmtend, ncmtend, nimtend, & ! Out
                                      rvm_mc, rcm_mc, thlm_mc ) ! Out

    ! Description:
    !   Reject a COAMPS request when the imported scheme was not compiled.

    use clubb_precision, only: time_precision, core_rknd
    use grid_class, only: grid
    use stats_netcdf, only: stats_type

    implicit none

    !--------------------------- Input Variables ---------------------------
    type(grid), intent(in) :: gr ! CLUBB vertical grids

    integer, intent(in) :: ngrdcol ! Number of grid columns [-]

    character(len=*), intent(in) :: runtype ! Case name

    real( kind = time_precision ), intent(in) :: timea_in ! Current model time [s]
    real( kind = core_rknd ), intent(in) :: deltf_in ! Main timestep [s]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), intent(in) :: &
      rtm ! Total water mixing ratio [kg/kg]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzm), intent(in) :: &
      wm_zm ! Vertical wind on momentum levels [m/s]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), intent(in) :: &
      p_in_Pa, & ! Pressure [Pa]
      exner,   & ! Exner function [-]
      rho,     & ! Density [kg/m^3]
      thlm,    & ! Liquid potential temperature [K]
      T_in_K,  & ! Absolute temperature computed for the full batch [K]
      rim,     & ! Ice water mixing ratio [kg/kg]
      rrm,     & ! Rain water mixing ratio [kg/kg]
      rgm,     & ! Graupel water mixing ratio [kg/kg]
      rsm,     & ! Snow water mixing ratio [kg/kg]
      rcm,     & ! Cloud water mixing ratio [kg/kg]
      Ncm,     & ! Number of cloud droplets [count/kg]
      Nrm,     & ! Number of rain drops [count/kg]
      Nim        ! Number of ice crystals [count/kg]

    integer, intent(in) :: saturation_formula ! Saturation formula selector [-]

    !------------------------ Input/Output Variables ------------------------
    type(stats_type), intent(inout) :: stats ! Statistics buffers and sample controls

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), intent(inout) :: &
      Nccnm ! Number of cloud nuclei [count/kg]

    !--------------------------- Output Variables ---------------------------
    ! All sedimentation velocities are on thermodynamic levels.
    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), intent(out) :: &
      Vrs,     & ! Snow mixing ratio fall speed [m/s]
      Vri,     & ! Pristine ice mixing ratio fall speed [m/s]
      Vrr,     & ! Rain mixing ratio fall speed [m/s]
      VNr,     & ! Rain concentration fall speed [m/s]
      Vrg,     & ! Graupel mixing ratio fall speed [m/s]
      ritend,  & ! Ice water tendency [kg/kg/s]
      rrtend,  & ! Rain water tendency [kg/kg/s]
      rgtend,  & ! Graupel water tendency [kg/kg/s]
      rstend,  & ! Snow water tendency [kg/kg/s]
      nrmtend, & ! Rain number tendency [count/kg/s]
      ncmtend, & ! Cloud droplet number tendency [count/kg/s]
      nimtend, & ! Ice number tendency [count/kg/s]
      rvm_mc,  & ! Vapor water tendency [kg/kg/s]
      rcm_mc,  & ! Cloud water tendency [kg/kg/s]
      thlm_mc    ! Liquid potential temperature tendency [K/s]

    error stop "Not compiled with COAMPS Microphysics"

  end subroutine coamps_microphys_driver
#endif

  end module coamps_microphys_driver_module
