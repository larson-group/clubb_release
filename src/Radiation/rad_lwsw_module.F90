! $Id$
!-----------------------------------------------------------------------
module rad_lwsw_module

  implicit none

  !public :: rad_lwsw

  public :: sunray_sw

  private ! Default Scope

  contains

!   subroutine rad_lwsw( qc3, rbm, dsigm, & 
!                        coamps_zm, coamps_zt,  & 
!                        Frad, Frad_LW, Frad_SW, & 
!                        radhtk, radht_LW, radht_SW, & 
!                        kk, l_center, xi_abs, F0, F1, kay,  & 
!                        radius, A, gc, Fs0, omega, & 
!                        l_sw_on, l_lw_on)

! ! Description:
! ! For the Larson group altocumulus cases
! ! Written by Vince Larson, Chris Golaz, Adam Smith, Michael Falk, and
! ! others for COAMPS and CLUBB
! !
! ! Subroutine to compute cloud radiation according to GCSS DYCOMS
! ! specifications.
! !
! !
! ! NOTE BY ADAM SMITH, 28 July 2006.
! ! We are attempting to compare CLUBB results with the COAMPS-LES model.
! ! To keep the two codes as similar as possible, this subroutine is a
! ! near-duplicate of the rad_lwsw subroutine  used in COAMPS_LES.  I
! ! have  modified variable declarations as needed to match CLUBB
! ! standards, and sections that do not apply to CLUBB have been removed.
! !
! ! For a full explanation of how this subroutine is implemented in CLUBB,
! ! please see the nov11_altocu_tndcy or jun25_altocu_tndcy subroutines
! ! within gcss.F.
! !
! ! REFERENCES (for rad_lwsw and sunray_sw):
! ! Bluestein, H. B., 1992: Synoptic-Dynamic Meteorology in
! !       Midlatitudes, Volume I.  Oxford Univerity Press, 431 pp.
! ! Duynkerke, P.G. et al., 2005: Observations and numerical simulations
! !       of the diurnal cycle of the EUROCS stratocumulus case.  In
! !       press, special EUROCS issue, Quart. Jour. Roy. Met. Soc.
! ! Salby, M. L., 1996: Fundamentals of Atmospheric Physics.  Academic
! !       Press, 627 pp.
! ! Shettle, E. P., and J. A. Weinman, 1970: The transfer of solar
! !       irradiance through inhomogenous turbid atmospheres evaluated
! !       by Eddington's approximation.  J. Atm. Sci., 27, 1048-1055.
! ! Stevens, B. et al., 2005: Evaluation of large-eddy simulations via
! !       observations of nocturnal marine stratocumulus.  Submitted to
! !       Mon. Wea. Rev.
! ! Wallace, J. M., and P. V. Hobbs, 1977: Atmospheric Science: An
! !       Introductory Survey.  Academic Press, 467 pp.
! !-----------------------------------------------------------------------

!     use constants_clubb, only: Cp ! Variable(s)

!     use interpolation, only: lin_interpolate_two_points ! Procedure(s)

!     use clubb_precision, only: core_rknd ! Variable(s)

!     implicit none

! ! Input Variables
!     integer, intent(in) :: kk ! Number of vertical levels   [-]

!     real( kind = core_rknd ), dimension(kk), intent(in) ::  & 
!       qc3,   & ! Cloud water mixing ratio at time t + dt       [kg/kg]
!       rbm,   & ! Density of reference state at mass levels     [kg/m^3]
!       dsigm    ! Thickness of sigma (mass) levels              [m]

!     real( kind = core_rknd ), dimension(kk+1), intent(in) :: & 
!       coamps_zm,  & ! Altitude of momentum levels w/ COAMPS grid indices      [m]
!       coamps_zt  ! Altitude of thermodynamic levels w/ COAMPS grid indices [m]

!     real( kind = core_rknd ), intent(in) ::  & 
!       xi_abs,  & ! Cosine of the solar zenith angle               [-]
!       F0,      & ! Coefficient for cloud top heating, see Stevens [W/m^2]
!       F1,      & !     "                      "                   [W/m^2]
!       kay,     & ! A "constant" according to Duynkerke eqn.5,
!              ! where his value is 130 m^2/kg.                 [m^2/kg]
!       radius     ! Effective droplet radius.                      [m]

!     real( kind = core_rknd ), intent(in) ::  & 
!       A,     & ! Albedo- sea surface, according to Lenderink             [-]
!       gc,    & ! Asymmetry parameter, "g" in Duynkerke.                  [-] 
!       Fs0,   & ! Incident incoming SW insolation at cloud top in 
!            ! direction of the incoming beam (not the vertical)       [W/m^2]
!       omega ! Single-scattering albedo                                [-]

!     logical, intent(in) ::  & 
!       l_center,  & ! Use centered differences?     [-]
!       l_sw_on,   & ! Is shortwave radiation on?    [-]
!       l_lw_on      ! Is longwave radiation on?     [-]

! ! Output Variables
!     real( kind = core_rknd ), dimension(kk), intent(out) ::  & 
!       radhtk,    & ! Total radiational heating (dT/dt)           [K/s]
!       radht_SW,  & ! Shortwave component of radiational heating  [K/s]
!       radht_LW     ! Longwave component of radiational heating   [K/s]

!     real( kind = core_rknd ), dimension(kk+1), intent(out) ::  & 
!       Frad,     & ! Total radiative flux         [W/m^2]
!       Frad_SW,  & ! Shortwave radiative flux     [W/m^2]
!       Frad_LW     ! Longwave radiative flux      [W/m^2]

! ! Local Variables
!     real( kind = core_rknd ), dimension(kk+1) :: lwp ! Liquid water path from domain top [kg/m^2]
!     real( kind = core_rknd ), dimension(kk+1) :: & 
!       lwp_coamps_zm   ! Liquid water path interpolated
!     ! to COAMPS momentum (or w) levels  [kg/m^2]

!     integer :: k

! !-----------------------------------------------------------------------
! !
! !-> Addition by Adam Smith (09 July 2004)
! !-> NOTE: We want to execute simulations using only COAMPS radiation,
! !-> or only analytic radiation.  If we have lrad = t in the Nov.11
! !-> namelist, we MUST comment out the following call statement and do
! !-> loop.  If lrad = t and the analytic radiation is not commented out,
! !-> we will execute a STOP statement to cancel the simulation.
! !-> End of Adam's comment
! !
! !-----------------------------------------------------------------------

! !      if(lrad) then
! !        print *, "We cannot use both COAMPS and analytic radiation!"
! !        print *, "If you want to use analytic, set lrad = f."
! !        print *, "Otherwise, comment out the analytic radiation."
! !        stop
! !      endif


! !-----------------------------------------------------------------------
! !
! !<-- End of Adam's Addition
! !
! !-----------------------------------------------------------------------

! !-----------------------------------------------------------------------
! !
! !  Computation of liquid water path
! !
! ! Liquid water path is defined by Salby as "the column abundance of
! ! liquid water"-- the mass of water over the length of the column
! ! (eqn.9.44).  This formulation comes from Duynkerke eqn.6.
! ! We compute the path starting at the top of the model domain and
! ! compute a running total as we integrate downward.
! ! Comments by Adam Smith and Michael Falk.
! !
! !-----------------------------------------------------------------------

!     if ( l_lw_on ) then
!       lwp(1) = 0.0_core_rknd
!       do k=2,kk+1
!         lwp(k) = lwp(k-1) & 
!                + rbm(k-1)*qc3(k-1)*dsigm(k-1) !/aoz(i,j)
!       end do


! !-----------------------------------------------------------------------
! !
! !  Computation of radiative fluxes
! !
! ! Frad_LW(k) is the net longwave radiation emitted from cloud base and
! ! cloud top.  This formulation comes from Stevens et al eqn.3.  Their
! ! paper assumes a term related to the height of the inversion layer,
! ! which is not present in this case, so that term is dropped.  The F0
! ! term, for the heat emitted from cloud top, also is used in Duynkerke
! ! eqn.5, where their delta Ft is this case's F0, and their alpha is
! ! this case's kay.
! !
! ! sunray_sw is then called, returning Frad_SW.  The sum of the
! ! shortwave and longwave components then makes the total radiative
! ! flux, Frad.
! !
! ! Comments by Michael Falk, 26 January 2005.
! !
! !-----------------------------------------------------------------------

! !-----------------------------------------------------------------------
! !
! ! Computation of radiative fluxes on staggered grid
! !
! ! Frad (and its components Frad_LW and Frad_SW) should be computed on
! ! w points, not on mass points, which is apparent from its formulation
! ! and from its location in stats_sw instead of stats_sm.  The grid
! ! looks like this:
! !
! !
! ! -----Frad----------------------------------    k = 1  (w level)
! !     /    \            |-dwm
! ! -LWP------radht----------------------------    k = 1  (mass level)
! !     \    /            |-dmw
! ! -----Frad----------------------------------    k = 2  (w level)
! !     /    \
! ! -LWP------radht----------------------------    k = 2  (mass level)
! !     \    /
! ! -----Frad----------------------------------    k = 3  (w level)
! !     /    \
! ! -LWP------radht----------------------------    k = 3  (mass level)
! !
! ! If you consider Frad to take place on mass levels, then computing
! ! LWP is a forward difference and is only first-order accurate, while
! ! if Frad computed in between LWP levels, it is a centered difference
! ! which is second-order accurate.
! !
! ! The coding implementation requires that Frad depend on LWP(k) and
! ! LWP(k-1) since the w level for a given k is at a higher altitude
! ! than the mass level.  radht, back on mass levels, depends on Frad(k)
! ! and Frad(k+1).
! !
! ! ADDITIONAL NOTE: For clarification of terminology, a w level on the
! !                  COAMPS grid is equivalent to a momentum (m) level
! !                  on the CLUBB grid, and a mass level on the COAMPS
! !                  grid is equivalent to a thermodynamic (t) level on
! !                  the CLUBB grid.  Brian Griffin; May 10, 2008.
! !
! ! Additionally, these computations assume that the distance between
! ! mass levels (dsigma) is constant, and that the w levels (spaced by
! ! dsigmw) always fall exactly halfway in between the mass levels.  If
! ! this is not the case, consider dwm to be the distance between a w
! ! level and the mass level below it, and dmw to be the distance
! ! between a mass level and the w level below it.  Then, the
! ! formulation for Frad_LW, for instance, would use a weighted average:
! !
! ! (dwm/(dwm+dmw)) * lwp(k) + (dmw/(dwm+dmw)) * lwp(k-1)
! ! which, for dwm always == dmw, reduces to
! ! (1/2) * (lwp(k)) + (1/2) * (lwp(k-1))
! ! which is identical to the current formulation.
! ! ((lwp(k)+lwp(k-1))/2)
! !
! ! Comments by Michael Falk, 16 February 2005.
! !
! ! ADDITIONAL NOTE: The CLUBB parameterization is now set up to be
! !                  compatible with the use of a stretched
! !                  (or unevenly-spaced) grid, as well as with the use
! !                  of an evenly-spaced grid.  Interpolation functions
! !                  are used to compute any weighted averages, rather
! !                  than using general numbers such as (1/2), which is
! !                  compatible only with an evenly-spaced grid.
! !                  Brian Griffin; May 10, 2008.
! !
! !
! !-----------------------------------------------------------------------

!       ! Interpolate liquid water path (lwp) from COAMPS thermodynamic
!       ! (or mass) levels to COAMPS momentum (or w) levels.
!       do k = 2, kk+1, 1
!         lwp_coamps_zm(k) = lin_interpolate_two_points( coamps_zm(k), coamps_zt(k-1),  & 
!                                     coamps_zt(k), lwp(k-1), lwp(k) )
!       enddo
!       ! The value of liquid water path (lwp) at momentum (or w)
!       ! level 1 (the uppermost level) is defined to be 1/2 of the
!       ! value of liquid water path at thermodynamic level 1.
!       lwp_coamps_zm(1) = lwp(1)/2._core_rknd

!       if ( l_center ) then
!         Frad_LW(1) = F0 * exp( -kay * lwp_coamps_zm(1) ) & 
!                    + F1 * exp( -kay * & 
!                         ( lwp(kk+1) - lwp_coamps_zm(1) ) )

!         do k=2,kk+1
!           Frad_LW(k) = F0 * exp( -kay * lwp_coamps_zm(k) ) & 
!                      + F1 * exp( -kay * & 
!                           ( lwp(kk+1) - lwp_coamps_zm(k) ) )
!         enddo
!       else
!         do k=1,kk+1
!           Frad_LW(k) = F0 * exp( -kay * lwp(k) ) & 
!                      + F1 * exp( -kay * (lwp(kk+1)-lwp(k)) )
!         enddo
!       endif

!     else ! this 'else' means l_lw_on is .FALSE.
!       do k=1,kk+1
!         Frad_LW(k) = 0._core_rknd
!       enddo
!     endif

!     if ( l_sw_on ) then
!       ! TODO: rad_lwsw uses COAMPS top-down ordering; sunray_sw now uses
!       ! CLUBB bottom-up ordering with ngrdcol.  Since rad_lwsw has no
!       ! callers, this call is commented out pending removal of rad_lwsw.
!       !call sunray_sw( qc3, rbm, xi_abs, dsigm, kk, &
!       !                coamps_zm, coamps_zt, &
!       !                radius, A, gc, Fs0, omega, l_center, &
!       !                Frad_SW )
!       Frad_SW(:) = 0._core_rknd
!     else
!       do k=1,kk+1
!         Frad_SW(k) = 0._core_rknd
!       enddo
!     endif

!     do k=1,kk+1
!       Frad(k) = Frad_LW(k) + Frad_SW(k)
!     enddo


! !-----------------------------------------------------------------------
! !
! !  Computation of radiative heating rate
! !
! ! radhtk is a one-dimensional version of the radht(i,k) variable.
! ! Since nov11_rad is being called for each column, we only need to
! ! compute radiation variables in one dimension at a time.  radht(i,k)
! ! is still used in nov11_forcing, and after the call to nov11_rad,
! ! radhtk(k) is converted to radht(i,k).
! !
! ! Radhtk could be derived from Salby 8.72.  More info in Wallace and
! ! Hobbs unit 8.7.  Radhtk would be equal to (q dot) / (Cp), which in
! ! turn is equal to (1 / (rho*Cp)) * (-dF/dz).  The units for this are
! ! K/s.  (units for q dot, diabatic heating rate per mass per time, are
! ! J/kg/s).
! !
! ! The shortwave and longwave contributions to radiative heating are
! ! also computed, so they can be used in statistics.
! !
! ! Comments by Michael Falk, 26 January 2005.
! !
! !-----------------------------------------------------------------------

!     do k=1,kk
!       radhtk(k)   = (-1.0_core_rknd/(Cp*rbm(k))) & 
!                 * (Frad(k)-Frad(k+1))/dsigm(k)
!       radht_SW(k) = (-1.0_core_rknd/(Cp*rbm(k))) & 
!                 * (Frad_SW(k)-Frad_SW(k+1))/dsigm(k)
!       radht_LW(k) = (-1.0_core_rknd/(Cp*rbm(k))) & 
!                 * (Frad_LW(k)-Frad_LW(k+1))/dsigm(k)
!     enddo

!     return
!   end subroutine rad_lwsw

!-----------------------------------------------------------------------
  subroutine sunray_sw( ngrdcol, nzt, rcm, rho, xi_abs, dzt, &
                        zm, zt, &
                        radius, A, gc, Fs0, omega, l_center, &
                        Frad_SW )

! Description:
! for CLEX altocumulus case
! Written by Geert Lenderink for implementation of Shettle and
! Weinman's formulation for radiative flux into the EUROCS
! stratocumulus case.
!
! Adapted by Vince Larson, Chris Golaz, Adam Smith, Michael Falk, and
! others for COAMPS and CLUBB.
!
! Subroutine to compute shortwave radiative flux.
!
!
! The code for sunray_sw was slightly reconstructed in order to more
! closely follow Geert Lenderink's code for the Duynkerke et al.
! EUROCS case.  The original formulation used in that paper comes from
! Shettle and Weinman case.
!
! Tau is now computed inside this routine, as it's not needed in
! nov11_rad.  Tau is also computed for each layer instead of as the
! total optical depth from the top of the domain to the bottom of the
! layer being computed.  This makes tau zero everywhere outside of
! cloud.  F_diff and F_dir are also zero outside of cloud.
!
! Comments by Michael Falk, 26 January 2005.
!
! ADDITION TO COMMENT BY ADAM SMITH, 28 July 2006.
! We are attempting to compare CLUBB results with the COAMPS-LES model.
! To keep the two codes as similar as possible, this subroutine is a
! near-duplicate of the rad_lwsw subroutine  used in COAMPS_LES.  I
! have  modified variable declarations as needed to match CLUBB
! standards, and sections that do not apply to CLUBB have been removed.
!
! For a full explanation of how this subroutine is implemented in CLUBB,
! please see the nov11_altocu_tndcy or jun25_altocu_tndcy subroutines
! within gcss.F.
!
! This version operates natively on CLUBB's bottom-up grid ordering
! (k=1 at surface, k=nzt at top for thermo levels, k=nzm=nzt+1 at
! top for momentum levels) and handles multiple grid columns.
! In CLUBB's grid, thermo level k sits between momentum levels k and
! k+1.  Optical depth (taupath) accumulates from the top of the
! domain downward, i.e. from momentum level nzm down to level 1.

! REFERENCES:
! see subroutine rad_lwsw.
!-----------------------------------------------------------------------

    use constants_clubb, only: &
      rho_lw, & ! Constant(s)
      three_halves

    use interpolation, only: lin_interpolate_two_points ! Procedure(s)
    use clubb_precision, only: core_rknd ! Variable(s)

    implicit none

    ! Input variables

    integer, intent(in) :: ngrdcol, nzt

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) ::  &
      rcm, &
      rho, &
      dzt

    real( kind = core_rknd ), dimension(ngrdcol,nzt+1), intent(in) :: &
      zm  ! Altitude of momentum levels (CLUBB bottom-up)      [m]

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      zt  ! Altitude of thermodynamic levels (CLUBB bottom-up) [m]

    real( kind = core_rknd ), intent(in) ::  &
      xi_abs, &
      radius,  &
      A,  &
      gc,  &
      Fs0,  &
      omega

    logical, intent(in) ::  &
      l_center

    ! Output variables
    real( kind = core_rknd ), dimension(ngrdcol,nzt+1), intent(out) ::  &
      Frad_SW


! Local Variables
    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      tau,     & ! Optical depth of an incremental layer.        [-]
      taude      ! Delta-Eddington transformation of tau.        [-]

    real( kind = core_rknd ), dimension(ngrdcol) :: &
      tauc,    & ! Column total optical depth                    [-]
      taucde,  & ! D-E column total optical depth                [-]
      taupath, & ! Running cumulative D-E optical depth          [-]
      exmu0,   & ! exp(-taucde/xi_abs)                           [-]
      expk,    & ! exp(rk*taucde)                                [-]
      exmk,    & ! 1/expk                                        [-]
      c1, c2     ! Shettle-Weinman coefficients                  [-]

    real( kind = core_rknd ) :: omegade, &
         t1, t2, t3, &
         x1, x2, x3, rk, rk2, xi_abs2, rp, alpha, beta, rtt, &
         xp23p, xm23p, ap23b

    real( kind = core_rknd ) :: F_diff_k, F_dir_k

    integer :: i, k

    real( kind = core_rknd ) :: ff, gcde

    !$acc enter data create( tau, taude, tauc, taucde, taupath, &
    !$acc                    exmu0, expk, exmk, c1, c2 )

!-----------------------------------------------------------------------
!  CONSTANTS/PARAMETERS
!
!  values added by Michael Falk and Adam Smith
!
! ff     : gc^2, denoted "g^2" in Duynkerke.            Unit: none
! gcde   : Delta-Eddington transformation of gc.
!          Notated g-prime in Duynkerke eqn.18.         Unit: None
!
!-----------------------------------------------------------------------

    ff = gc*gc
    gcde = gc/(1.0_core_rknd+gc)

!-----------------------------------------------------------------------
!
!  Computation of tau and omega variables
!
!
! tauc   : column total optical depth                   Unit: none
! taupath: column total Delta-Eddington optical depth.  Unit: none
! omega  : single-scattering albedo                     Unit: none
! omegade: D-E omega-- from Duynkerke eqn.18.           Unit: none
! taucde : D-E tauc -- from Duynkerke eqn.18.           Unit: none
! taude  : D-E tau  -- from Duynkerke eqn.18.           Unit: none
!
! Comments by Michael Falk, 26 January 2005
!
!-----------------------------------------------------------------------

    omegade = (1.0_core_rknd-ff)*omega/(1.0_core_rknd-omega*ff)

    ! Compute per-layer optical depth and column total.
    ! tauc accumulates over k per column, so parallelize over i only.
    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      tauc(i) = 0.0_core_rknd
      do k = 1, nzt
        tau(i,k) = three_halves * rcm(i,k) * rho(i,k) * dzt(i,k)  &
                 / radius / rho_lw
        tauc(i) = tauc(i) + tau(i,k)
      end do
    end do

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nzt
      do i = 1, ngrdcol
        taude(i,k) = (1.0_core_rknd-omega*ff)*tau(i,k)
      end do
    end do


!-----------------------------------------------------------------------
!
!  Computation of constants for radiative transfer equations
!
!  These variables come from Duynkerke eqn.20, which, with slight
!  modifications, matches Shettle and Weinman between eqns.12b and 13.
!  Duynkerke uses slightly different formulations than Shettle and
!  Weinman:
!
! F0(S&W)    = F0(Duynkerke)*pi.
! alpha(S&W) = alpha(Duynkerke)*F0(S&W).
! beta(S&W)  = beta(Duynkerke)*F0(S&W).
! c1(S&W)    = c1(Duynkerke)*F0(S&W).
! c2(S&W)    = c2(Duynkerke)*F0(S&W).
! x3(S&W)    = x3(Duynkerke)*F0(S&W).
! y3(S&W)    = y3(Duynkerke)*F0(S&W).
!
! F0 is divided out of each term in several equations, and then
! reintroduced when the actual radiative flux is computed.  The
! computations here follow Lenderink's original sunray_sw code which
! uses the Duynkerke formulations for these variables.
!
! x1     : term 1 in k equation                         Unit: none
! x2     : term 2 in k equation                         Unit: none
! rk     : k equation                                   Unit: none
! rk2    : k equation squared                           Unit: none
! x3     : term in denominator of alpha and beta        Unit: none
! rp     : p equation                                   Unit: none
! alpha  : alpha equation                               Unit: none
! beta   : beta equation                                Unit: none
!
! The following variables are used by Lenderink to solve for
! Duynkerke's parameters C1 and C2.  They are all dimensionless.
!
! rtt    : 2/3.
! exmu0  : exponential term used in S&W eqn.14-- originally from
!          eqn 1 (also Salby 9.35) in the source function for SW
!          radiation
! expk   : one of the coefficients of C2 on the left side of Shettle
!          and Weinman eqn.14, originally from eqn.12a and eqn.12b
! exmk   : reciprocal of expk, one of the coefficients of C1 on the
!          left side of Shettle and Weinman eqn.14, originally from
!          eqn.12a and eqn.12b
! xp23p  : coefficient of C1 - left side of Shettle and Weinman eqn.13
! xm23p  : coefficient of C2 - left side of Shettle and Weinman eqn.13
! ap23b  : right side of Shettle and Weinman eqn.13
! t1     : the other coefficient of C1 on left side of Shettle and
!          Weinman eqn.14
! t2     : the other coefficient of C2 on left side of Shettle and
!          Weinman eqn.14
! t3     : the coefficient of exmu0 on the right side of Shettle and
!          Weinman eqn.14
!
! Comments by Michael Falk, 26 January 2005
!
!-----------------------------------------------------------------------

    x1 = 1.0_core_rknd-omegade*gcde
    x2 = 1.0_core_rknd-omegade
    rk = sqrt( 3.0_core_rknd*x2*x1 )
    xi_abs2 = xi_abs*xi_abs
    rk2 = rk*rk
    x3 = 4.0_core_rknd*(1.0_core_rknd-rk2*xi_abs2)
    rp = sqrt( 3.0_core_rknd*x2/x1 )
    alpha = 3.0_core_rknd*omegade*xi_abs2*(1.0_core_rknd+gcde*x2)/x3
    beta = 3.0_core_rknd*omegade*xi_abs*(1.0_core_rknd+3.0_core_rknd*gcde*xi_abs2*x2)/x3

    rtt = 2.0_core_rknd/3.0_core_rknd
    xp23p = 1.0_core_rknd+rtt*rp
    xm23p = 1.0_core_rknd-rtt*rp
    ap23b = alpha+rtt*beta

    t1 = 1._core_rknd-A-rtt*(1.0_core_rknd+A)*rp
    t2 = 1._core_rknd-A+rtt*(1.0_core_rknd+A)*rp
    t3 = (1._core_rknd-A)*alpha-rtt*(1._core_rknd+A)*beta+A*xi_abs


!-----------------------------------------------------------------------
!
! Shettle and Weinman 13 and 14, adapted into Duynkerke, give two
! equations and two unknowns which can be solved by linear combination
! (C2) and then substitution (C1).
!
! exmu0, expk, exmk, taucde are per-column since they depend on the
! column total optical depth (tauc).
!
! Comments by Michael Falk, 26 January 2005
!
!-----------------------------------------------------------------------

    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      taucde(i) = (1.0_core_rknd-omega*ff)*tauc(i)
      exmu0(i) = exp( -taucde(i)/xi_abs )
      expk(i) = exp( rk*taucde(i) )
      exmk(i) = 1.0_core_rknd/expk(i)

      c2(i) = (xp23p*t3*exmu0(i)-t1*ap23b*exmk(i)) &
            / (xp23p*t2*expk(i)-xm23p*t1*exmk(i))
      c1(i) = (ap23b-c2(i)*xm23p)/xp23p
    end do


!-----------------------------------------------------------------------
!
!  Computation of diffuse and direct shortwave flux
!
! F_diff is the first term in Duynkerke eqn.19.  The F0 and pi
! constants which are divided out in the Duynkerke formulation are
! reintroduced here.
!
! Duynkerke eqn.19's F_diff term comes from Shettle and Weinman eqn.8,
! where F_diff = F(upward)-F(downward).  Then F_diff = (-4/3)*pi*I1,
! where I1 is solved in Shettle and Weinman eqn.12b.  Capital P in
! Shettle and Weinman eqn.12b should actually be a lowercase p.
!
! F_dir is the second term in Duynkerke eqn.19.
!
! The negative sign for F_diff and F_dir is related to the definition
! of which direction is a positive flux.
!
! Comments by Michael Falk, 26 January 2005
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
! Computation of shortwave fluxes on staggered grid
!
! For a full explanation see the "Computation of radiative fluxes on
! staggered grid" section above.  For Frad_SW to be correctly computed
! on w levels, the non-constant component of F_diff and F_dir, we
! compute taupath on w levels as a centered difference between tau
! values on mass levels.
!
! This uses CLUBB's bottom-up grid ordering (k=1 at surface):
!
! --------taupath-->--Frad_SW----------------    k = nzt+1  (top m level)
!        /                   \
! -taude----------------------               k = nzt    (top t level)
!        \                   /
! --------taupath-->--Frad_SW----------------    k = nzt    (m level)
!        /                   \
! -taude----------------------               k = nzt-1  (t level)
!        \                   /
! --------taupath-->--Frad_SW----------------    k = nzt-1  (m level)
!                        ...
! --------taupath-->--Frad_SW----------------    k = 1      (sfc m level)
!
! At momentum level k, the thermo level above is k and below is k-1.
! taupath accumulates from the top (k=nzt+1) downward to the surface
! (k=1).
!
! Vince Larson changed the F variables to w levels  03 Feb 2005
! Michael Falk changed the loop to start at k=2 and then solved
! separately for k=1 so the array didn't go out of bounds.
!
! This code makes the same assumption as above that dwm=dmw.
!
! Comments by Michael Falk, 16 February 2005.                          c
!                                                                      c
! ADDITIONAL NOTE: The CLUBB parameterization is now set up to be
!                  compatible with the use of a stretched
!                  (or unevenly-spaced) grid, as well as with the use
!                  of an evenly-spaced grid.  Therefore, dwm is not
!                  necessarily equal to dmw.  Interpolation functions
!                  are used to compute any weighted averages, rather
!                  than using general numbers such as (1/2), which is
!                  compatible only with an evenly-spaced grid.
!                  Brian Griffin; May 10, 2008.
!
!-----------------------------------------------------------------------

    ! taupath has a sequential dependency over k per column,
    ! so parallelize over columns with the k sweep sequential inside.
    !$acc parallel loop default(present)
    do i = 1, ngrdcol

      ! Top momentum level (k = nzt+1)
      if ( l_center ) then
        taupath(i) = 0.5_core_rknd*taude(i,nzt)
      else
        taupath(i) = 0._core_rknd
      endif

      F_diff_k = (-4.0_core_rknd/3.0_core_rknd) * Fs0 &
                * ( rp * ( c1(i)*exp( -rk*taupath(i) ) &
                         - c2(i)*exp( rk*taupath(i) ) ) &
                  - beta*exp( -taupath(i)/xi_abs ) )
      F_dir_k = -Fs0*xi_abs*exp( -taupath(i)/xi_abs )
      Frad_SW(i, nzt+1) = F_diff_k + F_dir_k

      ! Interior momentum levels (k = nzt down to 2).
      ! Going from momentum k+1 down to k, we cross thermo level k.
      ! At momentum k, the thermo level above is k and below is k-1.
      do k = nzt, 2, -1
        if ( l_center ) then
          taupath(i) = taupath(i) &
                     + lin_interpolate_two_points( zm(i,k), zt(i,k),  &
                                    zt(i,k-1), taude(i,k), taude(i,k-1) )
        else
          taupath(i) = taupath(i) + taude(i,k-1)
        endif

        F_diff_k = (-4.0_core_rknd/3.0_core_rknd) * Fs0 &
                  * ( rp * ( c1(i)*exp( -rk*taupath(i) ) &
                           - c2(i)*exp( rk*taupath(i) ) ) &
                    - beta*exp( -taupath(i)/xi_abs ) )
        F_dir_k = -Fs0*xi_abs*exp( -taupath(i)/xi_abs )
        Frad_SW(i, k) = F_diff_k + F_dir_k
      end do

      ! Bottom momentum level (k = 1).  The ghost thermo level below the
      ! surface has taude_ghost = taude(:,1), so both centered and
      ! non-centered cases simply add taude(:,1).
      taupath(i) = taupath(i) + taude(i,1)

      F_diff_k = (-4.0_core_rknd/3.0_core_rknd) * Fs0 &
                * ( rp * ( c1(i)*exp( -rk*taupath(i) ) &
                         - c2(i)*exp( rk*taupath(i) ) ) &
                  - beta*exp( -taupath(i)/xi_abs ) )
      F_dir_k = -Fs0*xi_abs*exp( -taupath(i)/xi_abs )
      Frad_SW(i, 1) = F_diff_k + F_dir_k

    end do ! i=1..ngrdcol

    !$acc exit data delete( tau, taude, tauc, taucde, taupath, &
    !$acc                   exmu0, expk, exmk, c1, c2 )

    return
  end subroutine sunray_sw

end module rad_lwsw_module
