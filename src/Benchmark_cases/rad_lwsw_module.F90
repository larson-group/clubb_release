! $Id$
!-----------------------------------------------------------------------
module rad_lwsw_module

  implicit none

  public :: rad_lwsw

  public :: sunray_sw

  private ! Default Scope

  contains

  subroutine rad_lwsw( qc3, rbm, dsigm, & 
                       coamps_zm, coamps_zt,  & 
                       Frad, Frad_LW, Frad_SW, & 
                       radhtk, radht_LW, radht_SW, & 
                       kk, l_center, xi_abs, F0, F1, kay,  & 
                       radius, A, gc, Fs0, omega, & 
                       l_sw_on, l_lw_on)

! Description:
! For the Larson group altocumulus cases
! Written by Vince Larson, Chris Golaz, Adam Smith, Michael Falk, and
! others for COAMPS and CLUBB
!
! Subroutine to compute cloud radiation according to GCSS DYCOMS
! specifications.
!
!
! NOTE BY ADAM SMITH, 28 July 2006.
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
! REFERENCES (for rad_lwsw and sunray_sw):
! Bluestein, H. B., 1992: Synoptic-Dynamic Meteorology in
!       Midlatitudes, Volume I.  Oxford Univerity Press, 431 pp.
! Duynkerke, P.G. et al., 2005: Observations and numerical simulations
!       of the diurnal cycle of the EUROCS stratocumulus case.  In
!       press, special EUROCS issue, Quart. Jour. Roy. Met. Soc.
! Salby, M. L., 1996: Fundamentals of Atmospheric Physics.  Academic
!       Press, 627 pp.
! Shettle, E. P., and J. A. Weinman, 1970: The transfer of solar
!       irradiance through inhomogenous turbid atmospheres evaluated
!       by Eddington's approximation.  J. Atm. Sci., 27, 1048-1055.
! Stevens, B. et al., 2005: Evaluation of large-eddy simulations via
!       observations of nocturnal marine stratocumulus.  Submitted to
!       Mon. Wea. Rev.
! Wallace, J. M., and P. V. Hobbs, 1977: Atmospheric Science: An
!       Introductory Survey.  Academic Press, 467 pp.
!-----------------------------------------------------------------------

    use constants_clubb, only: Cp ! Variable(s)

    use interpolation, only: lin_interpolate_two_points ! Procedure(s)

    use clubb_precision, only: core_rknd ! Variable(s)

    implicit none

! Input Variables
    integer, intent(in) :: kk ! Number of vertical levels   [-]

    real( kind = core_rknd ), dimension(kk), intent(in) ::  & 
      qc3,   & ! Cloud water mixing ratio at time t + dt       [kg/kg]
      rbm,   & ! Density of reference state at mass levels     [kg/m^3]
      dsigm    ! Thickness of sigma (mass) levels              [m]

    real( kind = core_rknd ), dimension(kk+1), intent(in) :: & 
      coamps_zm,  & ! Altitude of momentum levels w/ COAMPS grid indices      [m]
      coamps_zt  ! Altitude of thermodynamic levels w/ COAMPS grid indices [m]

    real( kind = core_rknd ), intent(in) ::  & 
      xi_abs,  & ! Cosine of the solar zenith angle               [-]
      F0,      & ! Coefficient for cloud top heating, see Stevens [W/m^2]
      F1,      & !     "                      "                   [W/m^2]
      kay,     & ! A "constant" according to Duynkerke eqn.5,
             ! where his value is 130 m^2/kg.                 [m^2/kg]
      radius     ! Effective droplet radius.                      [m]

    real( kind = core_rknd ), intent(in) ::  & 
      A,     & ! Albedo- sea surface, according to Lenderink             [-]
      gc,    & ! Asymmetry parameter, "g" in Duynkerke.                  [-] 
      Fs0,   & ! Incident incoming SW insolation at cloud top in 
           ! direction of the incoming beam (not the vertical)       [W/m^2]
      omega ! Single-scattering albedo                                [-]

    logical, intent(in) ::  & 
      l_center,  & ! Use centered differences?     [-]
      l_sw_on,   & ! Is shortwave radiation on?    [-]
      l_lw_on      ! Is longwave radiation on?     [-]

! Output Variables
    real( kind = core_rknd ), dimension(kk), intent(out) ::  & 
      radhtk,    & ! Total radiational heating (dT/dt)           [K/s]
      radht_SW,  & ! Shortwave component of radiational heating  [K/s]
      radht_LW     ! Longwave component of radiational heating   [K/s]

    real( kind = core_rknd ), dimension(kk+1), intent(out) ::  & 
      Frad,     & ! Total radiative flux         [W/m^2]
      Frad_SW,  & ! Shortwave radiative flux     [W/m^2]
      Frad_LW     ! Longwave radiative flux      [W/m^2]

! Local Variables
    real( kind = core_rknd ), dimension(kk+1) :: lwp ! Liquid water path from domain top [kg/m^2]
    real( kind = core_rknd ), dimension(kk+1) :: & 
      lwp_coamps_zm   ! Liquid water path interpolated
    ! to COAMPS momentum (or w) levels  [kg/m^2]

    integer :: k

!-----------------------------------------------------------------------
!
!-> Addition by Adam Smith (09 July 2004)
!-> NOTE: We want to execute simulations using only COAMPS radiation,
!-> or only analytic radiation.  If we have lrad = t in the Nov.11
!-> namelist, we MUST comment out the following call statement and do
!-> loop.  If lrad = t and the analytic radiation is not commented out,
!-> we will execute a STOP statement to cancel the simulation.
!-> End of Adam's comment
!
!-----------------------------------------------------------------------

!      if(lrad) then
!        print *, "We cannot use both COAMPS and analytic radiation!"
!        print *, "If you want to use analytic, set lrad = f."
!        print *, "Otherwise, comment out the analytic radiation."
!        stop
!      endif


!-----------------------------------------------------------------------
!
!<-- End of Adam's Addition
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  Computation of liquid water path
!
! Liquid water path is defined by Salby as "the column abundance of
! liquid water"-- the mass of water over the length of the column
! (eqn.9.44).  This formulation comes from Duynkerke eqn.6.
! We compute the path starting at the top of the model domain and
! compute a running total as we integrate downward.
! Comments by Adam Smith and Michael Falk.
!
!-----------------------------------------------------------------------

    if ( l_lw_on ) then
      lwp(1) = 0.0_core_rknd
      do k=2,kk+1
        lwp(k) = lwp(k-1) & 
               + rbm(k-1)*qc3(k-1)*dsigm(k-1) !/aoz(i,j)
      end do


!-----------------------------------------------------------------------
!
!  Computation of radiative fluxes
!
! Frad_LW(k) is the net longwave radiation emitted from cloud base and
! cloud top.  This formulation comes from Stevens et al eqn.3.  Their
! paper assumes a term related to the height of the inversion layer,
! which is not present in this case, so that term is dropped.  The F0
! term, for the heat emitted from cloud top, also is used in Duynkerke
! eqn.5, where their delta Ft is this case's F0, and their alpha is
! this case's kay.
!
! sunray_sw is then called, returning Frad_SW.  The sum of the
! shortwave and longwave components then makes the total radiative
! flux, Frad.
!
! Comments by Michael Falk, 26 January 2005.
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
! Computation of radiative fluxes on staggered grid
!
! Frad (and its components Frad_LW and Frad_SW) should be computed on
! w points, not on mass points, which is apparent from its formulation
! and from its location in stats_sw instead of stats_sm.  The grid
! looks like this:
!
!
! -----Frad----------------------------------    k = 1  (w level)
!     /    \            |-dwm
! -LWP------radht----------------------------    k = 1  (mass level)
!     \    /            |-dmw
! -----Frad----------------------------------    k = 2  (w level)
!     /    \
! -LWP------radht----------------------------    k = 2  (mass level)
!     \    /
! -----Frad----------------------------------    k = 3  (w level)
!     /    \
! -LWP------radht----------------------------    k = 3  (mass level)
!
! If you consider Frad to take place on mass levels, then computing
! LWP is a forward difference and is only first-order accurate, while
! if Frad computed in between LWP levels, it is a centered difference
! which is second-order accurate.
!
! The coding implementation requires that Frad depend on LWP(k) and
! LWP(k-1) since the w level for a given k is at a higher altitude
! than the mass level.  radht, back on mass levels, depends on Frad(k)
! and Frad(k+1).
!
! ADDITIONAL NOTE: For clarification of terminology, a w level on the
!                  COAMPS grid is equivalent to a momentum (m) level
!                  on the CLUBB grid, and a mass level on the COAMPS
!                  grid is equivalent to a thermodynamic (t) level on
!                  the CLUBB grid.  Brian Griffin; May 10, 2008.
!
! Additionally, these computations assume that the distance between
! mass levels (dsigma) is constant, and that the w levels (spaced by
! dsigmw) always fall exactly halfway in between the mass levels.  If
! this is not the case, consider dwm to be the distance between a w
! level and the mass level below it, and dmw to be the distance
! between a mass level and the w level below it.  Then, the
! formulation for Frad_LW, for instance, would use a weighted average:
!
! (dwm/(dwm+dmw)) * lwp(k) + (dmw/(dwm+dmw)) * lwp(k-1)
! which, for dwm always == dmw, reduces to
! (1/2) * (lwp(k)) + (1/2) * (lwp(k-1))
! which is identical to the current formulation.
! ((lwp(k)+lwp(k-1))/2)
!
! Comments by Michael Falk, 16 February 2005.
!
! ADDITIONAL NOTE: The CLUBB parameterization is now set up to be
!                  compatible with the use of a stretched
!                  (or unevenly-spaced) grid, as well as with the use
!                  of an evenly-spaced grid.  Interpolation functions
!                  are used to compute any weighted averages, rather
!                  than using general numbers such as (1/2), which is
!                  compatible only with an evenly-spaced grid.
!                  Brian Griffin; May 10, 2008.
!
!
!-----------------------------------------------------------------------

      ! Interpolate liquid water path (lwp) from COAMPS thermodynamic
      ! (or mass) levels to COAMPS momentum (or w) levels.
      do k = 2, kk+1, 1
        lwp_coamps_zm(k) = lin_interpolate_two_points( coamps_zm(k), coamps_zt(k-1),  & 
                                    coamps_zt(k), lwp(k-1), lwp(k) )
      enddo
      ! The value of liquid water path (lwp) at momentum (or w)
      ! level 1 (the uppermost level) is defined to be 1/2 of the
      ! value of liquid water path at thermodynamic level 1.
      lwp_coamps_zm(1) = lwp(1)/2._core_rknd

      if ( l_center ) then
        Frad_LW(1) = F0 * exp( -kay * lwp_coamps_zm(1) ) & 
                   + F1 * exp( -kay * & 
                        ( lwp(kk+1) - lwp_coamps_zm(1) ) )

        do k=2,kk+1
          Frad_LW(k) = F0 * exp( -kay * lwp_coamps_zm(k) ) & 
                     + F1 * exp( -kay * & 
                          ( lwp(kk+1) - lwp_coamps_zm(k) ) )
        enddo
      else
        do k=1,kk+1
          Frad_LW(k) = F0 * exp( -kay * lwp(k) ) & 
                     + F1 * exp( -kay * (lwp(kk+1)-lwp(k)) )
        enddo
      endif

    else ! this 'else' means l_lw_on is .FALSE.
      do k=1,kk+1
        Frad_LW(k) = 0._core_rknd
      enddo
    endif

    if ( l_sw_on ) then
      call sunray_sw( qc3, rbm, xi_abs, dsigm, kk, & 
                      coamps_zm, coamps_zt, & 
                      radius, A, gc, Fs0, omega, l_center, & 
                      Frad_SW )
    else
      do k=1,kk+1
        Frad_SW(k) = 0._core_rknd
      enddo
    endif

    do k=1,kk+1
      Frad(k) = Frad_LW(k) + Frad_SW(k)
    enddo


!-----------------------------------------------------------------------
!
!  Computation of radiative heating rate
!
! radhtk is a one-dimensional version of the radht(i,k) variable.
! Since nov11_rad is being called for each column, we only need to
! compute radiation variables in one dimension at a time.  radht(i,k)
! is still used in nov11_forcing, and after the call to nov11_rad,
! radhtk(k) is converted to radht(i,k).
!
! Radhtk could be derived from Salby 8.72.  More info in Wallace and
! Hobbs unit 8.7.  Radhtk would be equal to (q dot) / (Cp), which in
! turn is equal to (1 / (rho*Cp)) * (-dF/dz).  The units for this are
! K/s.  (units for q dot, diabatic heating rate per mass per time, are
! J/kg/s).
!
! The shortwave and longwave contributions to radiative heating are
! also computed, so they can be used in statistics.
!
! Comments by Michael Falk, 26 January 2005.
!
!-----------------------------------------------------------------------

    do k=1,kk
      radhtk(k)   = (-1.0_core_rknd/(Cp*rbm(k))) & 
                * (Frad(k)-Frad(k+1))/dsigm(k)
      radht_SW(k) = (-1.0_core_rknd/(Cp*rbm(k))) & 
                * (Frad_SW(k)-Frad_SW(k+1))/dsigm(k)
      radht_LW(k) = (-1.0_core_rknd/(Cp*rbm(k))) & 
                * (Frad_LW(k)-Frad_LW(k+1))/dsigm(k)
    enddo

    return
  end subroutine rad_lwsw

!-----------------------------------------------------------------------
  subroutine sunray_sw( qc3, rbm, xi_abs, dsigm, kk, & 
                        coamps_zm, coamps_zt, & 
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

    integer, intent(in) :: kk

    real( kind = core_rknd ), dimension(kk), intent(in) ::  & 
      qc3, & 
      rbm, & 
      dsigm

    real( kind = core_rknd ), dimension(kk+1), intent(in) :: & 
      coamps_zm,  & ! Altitude of momentum levels w/ COAMPS grid indices      [m]
      coamps_zt  ! Altitude of thermodynamic levels w/ COAMPS grid indices [m]

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
    real( kind = core_rknd ), dimension(kk+1), intent(out) ::  & 
      Frad_SW


! Local Variables
    real( kind = core_rknd ), dimension(kk+1) :: & 
      tau,     & ! Optical depth of an incremental layer.        [-]
      taude,   & ! Delta-Eddington transformation of tau.        [-]
      F_diff,  & ! Diffuse component of SW radiation             [W/m^2]
      F_dir      ! Diffuse component of LW radiation             [W/m^2]

    real( kind = core_rknd ) :: taupath, tauc, t1, t2, t3, c1, c2, omegade, & 
         x1, x2, x3, rk, rk2, xi_abs2, rp, alpha, beta, rtt, & 
         exmu0, expk, exmk, xp23p, xm23p, ap23b, taucde

    integer :: k

    real( kind = core_rknd ) :: ff, gcde

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


    tauc = 0.0_core_rknd
    do k=1, kk

      tau(k) = three_halves * qc3(k) * rbm(k) * dsigm(k)  & !/ aoz(i,j)
             / radius / rho_lw
      tauc = tauc + tau(k)
    enddo
    tau(kk+1) = tau(kk)

    omegade = (1.0_core_rknd-ff)*omega/(1.0_core_rknd-omega*ff)
    taucde = (1.0_core_rknd-omega*ff)*tauc

    do k=1, kk+1
      taude(k) = (1.0_core_rknd-omega*ff)*tau(k)
    enddo


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
    exmu0 = exp( -taucde/xi_abs )
    expk = exp( rk*taucde )
    exmk = 1.0_core_rknd/expk
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
! Comments by Michael Falk, 26 January 2005
!
!-----------------------------------------------------------------------

    c2 = (xp23p*t3*exmu0-t1*ap23b*exmk) & 
       / (xp23p*t2*expk-xm23p*t1*exmk)
    c1 = (ap23b-c2*xm23p)/xp23p


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
! --------taupath-->--Frad_SW----------------    k = 1  (w level)
!        /                   \        |-dwm
! -taude----------------------radht----------    k = 1  (mass level)
!        \                   /        |-dmw
! --------taupath-->--Frad_SW----------------    k = 2  (w level)
!        /                   \
! -taude----------------------radht----------    k = 2  (mass level)
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

    if ( l_center ) then
      taupath = 0.5_core_rknd*taude(1)
    else
      taupath = 0._core_rknd
    endif

    F_diff(1) = (-4.0_core_rknd/3.0_core_rknd) * Fs0 & 
              * (  & 
                 rp * & 
                     (  & 
                        c1*exp( -rk*taupath ) & 
                      - c2*exp( rk*taupath ) & 
                     ) & 
                 - beta*exp( -taupath/xi_abs )  & 
                )
    F_dir(1) = -Fs0*xi_abs*exp( -taupath/xi_abs )
    Frad_SW(1) = F_diff(1) + F_dir(1)

    do k = 2, kk+1

      if ( l_center ) then
        taupath = taupath  & 
                + lin_interpolate_two_points( coamps_zm(k), coamps_zt(k-1),  & 
                           coamps_zt(k), taude(k-1), taude(k) )
      else
        taupath = taupath + taude(k)
      endif


      F_diff(k) = (-4.0_core_rknd/3.0_core_rknd) * Fs0 & 
                * (  & 
                   rp * & 
                       (  & 
                          c1*exp( -rk*taupath ) & 
                        - c2*exp( rk*taupath ) & 
                       ) & 
                   - beta*exp( -taupath/xi_abs )  & 
                  )
      F_dir(k) = -Fs0*xi_abs*exp( -taupath/xi_abs )
      Frad_SW(k) = F_diff(k) + F_dir(k)

    enddo ! k=2..kk+1

    return
  end subroutine sunray_sw

end module rad_lwsw_module
