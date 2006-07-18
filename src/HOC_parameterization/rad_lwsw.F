cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Subroutine rad_lwsw -- for 11 November 1999 altocumulus case         c
c Written by Vince Larson, Chris Golaz, Adam Smith, Michael Falk, and  c
c others for COAMPS.                                                   c
c                                                                      c
c Subroutine to compute cloud radiation according to GCSS DYCOMS       c
c specifications.                                                      c
c                                                                      c
c SUBROUTINES:                                                         c
c sunray_sw - to calculate diffuse and direct shortwave fluxes         c
c                                                                      c
c ARGUMENTS (IN):                                                      c
c th3   : Potential temperature at time t + dt          Unit: K        c
c qc3   : Cloud water mixing ratio at time t + dt       Unit: kg/kg    c
c rbm   : Density of reference state at mass levels     Unit: kg/m^3   c
c dsigm : Thickness of sigma (mass) levels              Unit: m        c
c aoz   : Derivative of sigma wrt z                     Unit: none     c
c i,j   : Current location on x,y grid                  Unit: none     c
c nne   : ?? (number of nest?)                          Unit: ??       c
c kk    : Number of vertical levels                     Unit: none     c
c lrad  : True/false-- is COAMPS radiation on?		Unit: none     c
c center: True/false-- use centered differences?        Unit: none     c
c                                                                      c
c xi_abs: Cosine of the solar zenith angle              Unit: none     c
c         (.4325 value is specific to Nov. 11 case)                    c
c F0    : Coefficient for cloud top heating, see Stevens               c
c                                                       Unit: W/m^2    c
c F1    : Coefficient for cloud base heating, see Stevens              c
c                                                       Unit: W/m^2    c
c kay   : "a constant" according to Duynkerke eqn.5, where his value   c
c          is 130 m^2/kg.                               Unit: m^2/kg   c
c radius: effective droplet radius.                     Unit: m        c
c A     : Albedo- sea surface, according to Lenderink.  Unit: none     c
c gc    : asymmetry parameter, "g" in Duynkerke.        Unit: none     c
c Fs0   : Incident incoming SW insolation at cloud top in direction    c
c         of the incoming beam (not the vertical).      Unit: W/m^2    c
c omega : single-scattering albedo                      Unit: none     c
c sw_on : Logical variable passed to radiation scheme--                c
c         is shortwave radiation on?                                   c
c lw_on : Logical variable passed to radiation scheme--                c
c         is longwave radiation on?                                    c
c                                                                      
c
c                                                                      c
c ARGUMENTS (OUT):                                                     c
c radhtk: dT/dt caused by radiation.                    Unit: K/s      c
c                                                                      c
c REFERENCES (for nov11_rad and sunray_sw):                            c
c Bluestein, H. B., 1992: Synoptic-Dynamic Meteorology in              c
c       Midlatitudes, Volume I.  Oxford Univerity Press, 431 pp.       c
c Duynkerke, P.G. et al., 2005: Observations and numerical simulations c
c       of the diurnal cycle of the EUROCS stratocumulus case.  In     c
c       press, special EUROCS issue, Quart. Jour. Roy. Met. Soc.       c
c Salby, M. L., 1996: Fundamentals of Atmospheric Physics.  Academic   c
c       Press, 627 pp.                                                 c
c Shettle, E. P., and J. A. Weinman, 1970: The transfer of solar       c
c       irradiance through inhomogenous turbid atmospheres evaluated   c
c       by Eddington's approximation.  J. Atm. Sci., 27, 1048-1055.    c
c Stevens, B. et al., 2005: Evaluation of large-eddy simulations via   c
c       observations of nocturnal marine stratocumulus.  Submitted to  c
c       Mon. Wea. Rev.                                                 c
c Wallace, J. M., and P. V. Hobbs, 1977: Atmospheric Science: An       c
c       Introductory Survey.  Academic Press, 467 pp.                  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine rad_lwsw(th3, qc3, rbm, dsigm, 
     .                  Frad, Frad_LW, Frad_SW,
     .                  radhtk, radht_LW, radht_SW,
     .                  kk, center,
     .                  xi_abs, F0, F1, kay, radius, A, gc, Fs0, omega,
     .                  sw_on, lw_on)

      use constants

c-->Chris Golaz addition for statistics
#ifdef STATS
      use statistics
#endif
c<--end of addition

      implicit none

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c INPUT VARIABLES                                                      c
c                                                                      c
c th3     : liquid potential temperature                Unit: K        c
c qc3     : liquid water mixing ratio                   Unit: kg/kg    c
c rbm     : air density                                 Unit: kg/m^3   c
c dsigm   : spacing between grid heights                Unit: m        c
c                                                                      c
c xi_abs  : cosine of the solar zenith angle            Unit: NONE     c
c F0      : magnitude of cloud top heating                             c
c F1      :                                                            c
c kay     :                                                            c
c radius  :                                                            c
c A       :                                                            c
c gc      :                                                            c
c Fs0     :                                                            c
c omega   :                                                            c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      REAL, DIMENSION(1:kk), INTENT(IN) :: th3
      REAL, DIMENSION(1:kk), INTENT(IN) :: qc3
      REAL, DIMENSION(1:kk), INTENT(IN) :: rbm
      REAL, DIMENSION(1:kk), INTENT(IN) :: dsigm

      REAL, INTENT(IN) :: xi_abs, F0, F1, kay, radius, A, gc, Fs0, omega
      INTEGER, INTENT(IN) :: kk
      LOGICAL, INTENT(IN) :: center, sw_on, lw_on

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c OUTPUT VARIABLES                                                     c
c                                                                      c
c radhtk  : Total radiational heating.                  Unit: K/s      c
c radht_sw: Shortwave component of radiational heating. Unit: K/s      c
c radht_lw: Longwave component of radiational heating.  Unit: K/s      c
c lwp     : Liquid water path from domain top.          Unit: kg/m^2   c
c Frad    : Total radiative flux                        Unit: W/m^2    c
c Frad_SW : Shortwave radiative flux                    Unit: W/m^2    c
c Frad_LW : Longwave radiative flux                     Unit: W/m^2    c
c k       : for loop in the z-direction                 Unit: none     c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      REAL, DIMENSION(1:kk),   INTENT(OUT) :: radhtk
      REAL, DIMENSION(1:kk),   INTENT(OUT) :: radht_sw
      REAL, DIMENSION(1:kk),   INTENT(OUT) :: radht_lw
      REAL, DIMENSION(1:kk+1), INTENT(OUT) :: Frad
      REAL, DIMENSION(1:kk+1), INTENT(OUT) :: Frad_SW ! SW radiative flux
      REAL, DIMENSION(1:kk+1), INTENT(OUT) :: Frad_LW ! LW radiative flux

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c LOCAL VARIABLES                                                      c
c                                                                      c
c lwp     : Liquid water path from domain top.          Unit: kg/m^2   c
c k       : for loop in the z-direction                 Unit: none     c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      real lwp      (kk+1)
      integer k


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c CONSTANTS/PARAMETERS                                                 c
c                                                                      c
c values added by Michael Falk and Adam Smith                          c
c                                                                      c
c rgas    : Gas constant for dry air                    Unit: J/kg/K   c
c cp      : Specific heat constant for dry air          Unit: J/kg/K   c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c      real, parameter :: rgas   = 287.04
c      real, parameter :: cp     = 1005.


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c-> Addition by Adam Smith (09 July 2004)                              c
c-> NOTE: We want to execute simulations using only COAMPS radiation,  c
c-> or only analytic radiation.  If we have lrad = t in the Nov.11     c
c-> namelist, we MUST comment out the following call statement and do  c
c-> loop.  If lrad = t and the analytic radiation is not commented out,c
c-> we will execute a STOP statement to cancel the simulation.         c
c-> End of Adam's comment                                              c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c      if(lrad) then
c        print *, "We cannot use both COAMPS and analytic radiation!"
c        print *, "If you want to use analytic, set lrad = f."
c        print *, "Otherwise, comment out the analytic radiation."
c        stop
c      endif


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c<-- End of Adam's Addition                                            c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c  Computation of liquid water path                                    c
c                                                                      c
c Liquid water path is defined by Salby as "the column abundance of    c
c liquid water"-- the mass of water over the length of the column      c
c (eqn.9.44).  This formulation comes from Duynkerke eqn.6.            c
c We compute the path starting at the top of the model domain and      c
c compute a running total as we integrate downward.                    c
c Comments by Adam Smith and Michael Falk.                             c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if ( lw_on ) then
        lwp(1) = 0.0
        do k=2,kk+1
          lwp(k) = lwp(k-1)
     .           + rbm(k-1)*qc3(k-1)*dsigm(k-1) !/aoz(i,j)
        enddo


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c  Computation of radiative fluxes                                     c
c                                                                      c
c Frad_LW(k) is the net longwave radiation emitted from cloud base and c
c cloud top.  This formulation comes from Stevens et al eqn.3.  Their  c
c paper assumes a term related to the height of the inversion layer,   c
c which is not present in this case, so that term is dropped.  The F0  c
c term, for the heat emitted from cloud top, also is used in Duynkerke c
c eqn.5, where their delta Ft is this case's F0, and their alpha is    c
c this case's kay.                                                     c
c                                                                      c
c sunray_sw is then called, returning Frad_SW.  The sum of the         c
c shortwave and longwave components then makes the total radiative     c
c flux, Frad.                                                          c
c                                                                      c
c Comments by Michael Falk, 26 January 2005.                           c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c Computation of radiative fluxes on staggered grid                    c
c                                                                      c
c Frad (and its components Frad_LW and Frad_SW) should be computed on  c
c w points, not on mass points, which is apparent from its formulation c
c and from its location in stats_sw instead of stats_sm.  The grid     c
c looks like this:                                                     c
c                                                                      c
c                                                                      c
c -----Frad----------------------------------    k = 1  (w level)      c
c     /    \            |-dwm                                          c
c -LWP------radht----------------------------    k = 1  (mass level)   c
c     \    /            |-dmw                                          c
c -----Frad----------------------------------    k = 2  (w level)      c
c     /    \                                                           c
c -LWP------radht----------------------------    k = 2  (mass level)   c
c     \    /                                                           c
c -----Frad----------------------------------    k = 3  (w level)      c
c     /    \                                                           c
c -LWP------radht----------------------------    k = 3  (mass level)   c
c                                                                      c
c If you consider Frad to take place on mass levels, then computing    c
c LWP is a forward difference and is only first-order accurate, while  c
c if Frad computed in between LWP levels, it is a centered difference  c
c which is second-order accurate.                                      c
c                                                                      c
c The coding implementation requires that Frad depend on LWP(k) and    c
c LWP(k-1) since the w level for a given k is at a higher altitude     c
c than the mass level.  radht, back on mass levels, depends on Frad(k) c
c and Frad(k+1).                                                       c
c                                                                      c
c Additionally, these computations assume that the distance between    c
c mass levels (dsigma) is constant, and that the w levels (spaced by   c
c dsigmw) always fall exactly halfway in between the mass levels.  If  c
c this is not the case, consider dwm to be the distance between a w    c
c level and the mass level below it, and dmw to be the distance        c
c between a mass level and the w level below it.  Then, the            c
c formulation for Frad_LW, for instance, would use a weighted average: c
c                                                                      c
c (dwm/(dwm+dmw)) * lwp(k) + (dmw/(dwm+dmw)) * lwp(k-1)                c
c which, for dwm always == dmw, reduces to                             c
c (1/2) * (lwp(k)) + (1/2) * (lwp(k-1))                                c
c which is identical to the current formulation.                       c
c ((lwp(k)+lwp(k-1))/2)                                                c
c                                                                      c
c Comments by Michael Falk, 16 February 2005.                          c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        if (center) then
          Frad_LW(1) = F0 * exp( -kay * (lwp(1)/2))
     .               + F1 * exp( -kay *
     .                    (lwp(kk+1) - (lwp(1)/2) ) )

          do k=2,kk+1
            Frad_LW(k) = F0 * exp( -kay * ((lwp(k)+lwp(k-1))/2))
     .                 + F1 * exp( -kay *
     .                      (lwp(kk+1) - ((lwp(k)+lwp(k-1))/2) ) )
          enddo
        else
          do k=1,kk+1
            Frad_LW(k) = F0 * exp( -kay * lwp(k) )
     .                 + F1 * exp( -kay * (lwp(kk+1)-lwp(k)) )
          enddo
        endif

      else ! this 'else' means lw_on is .FALSE.
        do k=1,kk+1
          Frad_LW(k) = 0.
        enddo
      endif

      if ( sw_on ) then
        call sunray_sw      (qc3,rbm,xi_abs,dsigm,kk,
     .                       Frad_SW,radius,A,gc,Fs0,omega,center)
      else
        do k=1,kk+1
          Frad_SW(k) = 0.
        enddo
      endif

      do k=1,kk+1
        Frad(k) = Frad_LW(k) + Frad_SW(k)
      enddo


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c  Computation of radiative heating rate                               c
c                                                                      c
c radhtk is a one-dimensional version of the radht(i,k) variable.      c
c Since nov11_rad is being called for each column, we only need to     c
c compute radiation variables in one dimension at a time.  radht(i,k)  c
c is still used in nov11_forcing, and after the call to nov11_rad,     c
c radhtk(k) is converted to radht(i,k).                                c
c                                                                      c
c Radhtk could be derived from Salby 8.72.  More info in Wallace and   c
c Hobbs unit 8.7.  Radhtk would be equal to (q dot) / (Cp), which in   c
c turn is equal to (1 / (rho*Cp)) * (-dF/dz).  The units for this are  c
c K/s.  (units for q dot, diabatic heating rate per mass per time, are c
c J/kg/s).                                                             c
c                                                                      c
c The shortwave and longwave contributions to radiative heating are    c
c also computed, so they can be used in statistics.                    c
c                                                                      c
c Comments by Michael Falk, 26 January 2005.                           c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do k=1,kk
        radhtk(k) = (-1./(Cp*rbm(k)))
     .            * (Frad(k)-Frad(k+1))/dsigm(k) !*aoz(i,j)
        radht_sw(k) = (-1./(Cp*rbm(k)))
     .            * (Frad_SW(k)-Frad_SW(k+1))/dsigm(k) !*aoz(i,j)
        radht_lw(k) = (-1./(Cp*rbm(k)))
     .            * (Frad_LW(k)-Frad_LW(k+1))/dsigm(k) !*aoz(i,j)
      enddo

cc-->Chris Golaz addition for statistics
c#ifdef STATS
c      if ( lstats_samp .and. ifrad > 0 .and.
c     .     i.ge.imin_stats .and. i.le.imax_stats .and.
c     .     j.ge.jmin_stats .and. j.le.jmax_stats ) then
c        do k=1,kk+1
c          stats_sw%x(k,ifrad)
c     .    = stats_sw%x(k,ifrad) + Frad(k)
c         stats_sw%n(k,ifrad) = stats_sw%n(k,ifrad) + 1
c        enddo
c      endif
c#endif
cc<--end of addition
c
cc Michael Falk, 28 January 2005, addition for statistics
c#ifdef STATS
c      do k=1,kk
c        if ( lstats_samp .and. ifrad_lw > 0
c     .       .and. i.ge.imin_stats .and. i.le.imax_stats
c     .       .and. j.ge.jmin_stats .and. j.le.jmax_stats ) then
c          stats_sw%x(k,ifrad_lw) 
c     .    = stats_sw%x(k,ifrad_lw) + Frad_LW(k)
c          stats_sw%n(k,ifrad_lw) = stats_sw%n(k,ifrad_lw) + 1
c        endif
c        if ( lstats_samp .and. iradhtk > 0
c     .       .and. i.ge.imin_stats .and. i.le.imax_stats
c     .       .and. j.ge.jmin_stats .and. j.le.jmax_stats ) then
c          stats_sm%x(k,iradhtk) 
c     .    = stats_sm%x(k,iradhtk) + radhtk(k)
c          stats_sm%n(k,iradhtk) = stats_sm%n(k,iradhtk) + 1
c        endif
c        if ( lstats_samp .and. iradht_sw > 0
c     .       .and. i.ge.imin_stats .and. i.le.imax_stats
c     .       .and. j.ge.jmin_stats .and. j.le.jmax_stats ) then
c          stats_sm%x(k,iradht_sw) 
c     .    = stats_sm%x(k,iradht_sw) + radht_sw(k)
c          stats_sm%n(k,iradht_sw) = stats_sm%n(k,iradht_sw) + 1
c        endif
c        if ( lstats_samp .and. iradht_lw > 0
c     .       .and. i.ge.imin_stats .and. i.le.imax_stats
c     .       .and. j.ge.jmin_stats .and. j.le.jmax_stats ) then
c          stats_sm%x(k,iradht_lw) 
c     .    = stats_sm%x(k,iradht_lw) + radht_lw(k)
c          stats_sm%n(k,iradht_lw) = stats_sm%n(k,iradht_lw) + 1
c        endif
c      enddo                                         ! end k loop
c#endif
cc end addition for statistics

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c End of subroutine nov11_rad.                                         c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc










cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Subroutine sunray_sw -- for 11 November 1999 altocumulus case        c
c Written by Geert Lenderink for implementation of Shettle and         c
c Weinman's formulation for radiative flux into the EUROCS             c
c stratocumulus case.                                                  c
c                                                                      c
c Adapted by Vince Larson, Chris Golaz, Adam Smith, Michael Falk, and  c
c others for COAMPS.                                                   c
c                                                                      c
c Subroutine to compute shortwave radiative flux.                      c
c                                                                      c
c The code for sunray_sw was slightly reconstructed in order to more   c
c closely follow Geert Lenderink's code for the Duynkerke et al.       c
c EUROCS case.  The original formulation used in that paper comes from c
c Shettle and Weinman case.                                            c
c                                                                      c
c Tau is now computed inside this routine, as it's not needed in       c
c nov11_rad.  Tau is also computed for each layer instead of as the    c
c total optical depth from the top of the domain to the bottom of the  c
c layer being computed.  This makes tau zero everywhere outside of     c
c cloud.  F_diff and F_dir are also zero outside of cloud.             c
c                                                                      c
c Comments by Michael Falk, 26 January 2005.                           c
c                                                                      c
c SUBROUTINES:                                                         c
c none                                                                 c
c                                                                      c
c ARGUMENTS (IN):                                                      c
c qc3     : Cloud water mixing ratio at time t + dt     Unit: kg/kg    c
c rbm     : Density of reference state at mass levels   Unit: kg/m^3   c
c aoz     : Derivative of sigma wrt z                   Unit: none     c
c           aoz = dsigma / dz.  dsigm/aoz = dz.         Unit: m        c
c xi_abs  : Cosine of the solar zenith angle                           c
c           (.4325 value is specific to Nov. 11 case)   Unit: none     c
c dsigm   : Thickness of sigma (mass) levels            Unit: m        c
c kk      : Number of vertical levels                   Unit: none     c
c nne     : ?? (number of nest?)                        Unit: ??       c
c i,j     : Current location on x,y grid                Unit: none     c
c                                                                      c
c ARGUMENTS (OUT):                                                     c
c Frad_SW : Shortwave radiative flux                    Unit: W/m^2    c
c                                                                      c
c REFERENCES:                                                          c
c see subroutine rad_lwsw.                                             c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine sunray_sw(qc3,rbm,xi_abs,dsigm,kk,Frad_SW,
     .                  radius,A,gc,Fs0,omega,center)

c Michael Falk addition for statistics, 28 January 2005
#ifdef STATS
      use statistics
#endif

      implicit none

      real qc3    (kk)
      real rbm    (kk)

      real xi_abs
      real dsigm  (kk)

      integer kk
      integer nne
      integer i,j,k

      real radius, A, gc, Fs0, omega
      logical center

      real Frad_SW      (kk+1)


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c LOCAL VARIABLES                                                      c
c                                                                      c
c tau    : Optical depth of an incremental layer.      Unit: none      c
c taude  : Delta-Eddington transformation of tau.      Unit: none      c
c F_diff : Diffuse component of SW radiation           Unit: W/m^2     c
c F_dir  : Diffuse component of LW radiation           Unit: W/m^2     c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      real tau          (kk+1)
      real taude        (kk+1)
      real F_diff       (kk+1)
      real F_dir        (kk+1)

      real taupath, tauc, t1, t2, t3, c1, c2, omegade,
     .     x1, x2, x3, rk, rk2, xi_abs2, rp, alpha, beta, rtt,
     .     exmu0, expk, exmk, xp23p, xm23p, ap23b, am23b, taucde,
     .     pi


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  CONSTANTS/PARAMETERS                                                c
c                                                                      c
c  values added by Michael Falk and Adam Smith                         c
c                                                                      c
c rho_l  : the density of liquid water at 0 degrees C.  Unit: kg/m^3.  c
c cp     : Specific heat constant for dry air.          Unit: J/kg/K   c
c ff     : gc^2, notated "g^2" in Duynkerke.            Unit: none     c
c gcde   : Delta-Eddington transformation of gc.                       c
c          Notated g-prime in Duynkerke eqn.18.         Unit: None     c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      real, parameter :: rho_l = 1000.0
      real, parameter :: cp = 1005.0
      real ff, gcde
      ff = gc*gc
      gcde = gc/(1.0+gc)
      pi = 4.0*atan(1.0)


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c  Computation of tau and omega variables                              c
c                                                                      c
c                                                                      c
c tauc   : column total optical depth                   Unit: none     c
c taupath: column total Delta-Eddington optical depth.  Unit: none     c
c omega  : single-scattering albedo                     Unit: none     c
c omegade: D-E omega-- from Duynkerke eqn.18.           Unit: none     c
c taucde : D-E tauc -- from Duynkerke eqn.18.           Unit: none     c
c taude  : D-E tau  -- from Duynkerke eqn.18.           Unit: none     c
c                                                                      c
c Comments by Michael Falk, 26 January 2005                            c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      tauc = 0.0
      do k=1, kk
        tau(k) = 1.5 * qc3(k) * rbm(k) * dsigm(k) !/ aoz(i,j)
     .         / radius / rho_l
        tauc = tauc + tau(k)
      enddo
      tau(kk+1) = tau(kk)

      omegade = (1.0-ff)*omega/(1.0-omega*ff)
      taucde = (1.0-omega*ff)*tauc

      do k=1, kk+1
        taude(k) = (1.0-omega*ff)*tau(k)
      enddo


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c  Computation of constants for radiative transfer equations           c
c                                                                      c
c  These variables come from Duynkerke eqn.20, which, with slight      c
c  modifications, matches Shettle and Weinman between eqns.12b and 13. c
c  Duynkerke uses slightly different formulations than Shettle and     c
c  Weinman:                                                            c
c                                                                      c
c F0(S&W)    = F0(Duynkerke)*pi.                                       c
c alpha(S&W) = alpha(Duynkerke)*F0(S&W).                               c
c beta(S&W)  = beta(Duynkerke)*F0(S&W).                                c
c c1(S&W)    = c1(Duynkerke)*F0(S&W).                                  c
c c2(S&W)    = c2(Duynkerke)*F0(S&W).                                  c
c x3(S&W)    = x3(Duynkerke)*F0(S&W).                                  c
c y3(S&W)    = y3(Duynkerke)*F0(S&W).                                  c
c                                                                      c
c F0 is divided out of each term in several equations, and then        c
c reintroduced when the actual radiative flux is computed.  The        c
c computations here follow Lenderink's original sunray_sw code which   c
c uses the Duynkerke formulations for these variables.                 c
c                                                                      c
c x1     : term 1 in k equation                         Unit: none     c
c x2     : term 2 in k equation                         Unit: none     c
c rk     : k equation                                   Unit: none     c
c rk2    : k equation squared                           Unit: none     c
c x3     : term in denominator of alpha and beta        Unit: none     c
c rp     : p equation                                   Unit: none     c
c alpha  : alpha equation                               Unit: none     c   
c beta   : beta equation                                Unit: none     c
c                                                                      c
c The following variables are used by Lenderink to solve for           c
c Duynkerke's parameters C1 and C2.  They are all dimensionless.       c
c                                                                      c
c rtt    : 2/3.                                                        c
c exmu0  : exponential term used in S&W eqn.14-- originally from       c
c          eqn 1 (also Salby 9.35) in the source function for SW       c
c          radiation                                                   c
c expk   : one of the coefficients of C2 on the left side of Shettle   c
c          and Weinman eqn.14, originally from eqn.12a and eqn.12b     c
c exmk   : reciprocal of expk, one of the coefficients of C1 on the    c
c          left side of Shettle and Weinman eqn.14, originally from    c
c          eqn.12a and eqn.12b                                         c
c xp23p  : coefficient of C1 - left side of Shettle and Weinman eqn.13 c
c xm23p  : coefficient of C2 - left side of Shettle and Weinman eqn.13 c
c ap23b  : right side of Shettle and Weinman eqn.13                    c
c t1     : the other coefficient of C1 on left side of Shettle and     c
c          Weinman eqn.14                                              c
c t2     : the other coefficient of C2 on left side of Shettle and     c
c          Weinman eqn.14                                              c
c t3     : the coefficient of exmu0 on the right side of Shettle and   c
c          Weinman eqn.14                                              c
c                                                                      c
c Comments by Michael Falk, 26 January 2005                            c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      x1 = 1.0-omegade*gcde
      x2 = 1.0-omegade
      rk = sqrt(3.0*x2*x1)
      xi_abs2 = xi_abs*xi_abs
      rk2 = rk*rk
      x3 = 4.0*(1.0-rk2*xi_abs2)
      rp = sqrt(3.0*x2/x1)
      alpha = 3.0*omegade*xi_abs2*(1.0+gcde*x2)/x3
      beta = 3.0*omegade*xi_abs*(1.0+3.0*gcde*xi_abs2*x2)/x3

      rtt = 2.0/3.0
      exmu0 = exp(-taucde/xi_abs)
      expk = exp(rk*taucde)
      exmk = 1.0/expk
      xp23p = 1.0+rtt*rp
      xm23p = 1.0-rtt*rp
      ap23b = alpha+rtt*beta

      t1 = 1.-A-rtt*(1.0+A)*rp
      t2 = 1.-A+rtt*(1.0+A)*rp
      t3 = (1.-A)*alpha-rtt*(1.+A)*beta+A*xi_abs


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c Shettle and Weinman 13 and 14, adapted into Duynkerke, give two      c
c equations and two unknowns which can be solved by linear combination c
c (C2) and then substitution (C1).                                     c
c                                                                      c
c Comments by Michael Falk, 26 January 2005                            c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      c2 = (xp23p*t3*exmu0-t1*ap23b*exmk)
     .   / (xp23p*t2*expk-xm23p*t1*exmk)
      c1 = (ap23b-c2*xm23p)/xp23p


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c  Computation of diffuse and direct shortwave flux                    c
c                                                                      c
c F_diff is the first term in Duynkerke eqn.19.  The F0 and pi         c
c constants which are divided out in the Duynkerke formulation are     c
c reintroduced here.                                                   c
c                                                                      c
c Duynkerke eqn.19's F_diff term comes from Shettle and Weinman eqn.8, c
c where F_diff = F(upward)-F(downward).  Then F_diff = (-4/3)*pi*I1,   c
c where I1 is solved in Shettle and Weinman eqn.12b.  Capital P in     c
c Shettle and Weinman eqn.12b should actually be a lowercase p.        c
c                                                                      c
c F_dir is the second term in Duynkerke eqn.19.                        c
c                                                                      c
c The negative sign for F_diff and F_dir is related to the definition  c
c of which direction is a positive flux.                               c
c                                                                      c
c Comments by Michael Falk, 26 January 2005                            c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c Computation of shortwave fluxes on staggered grid                    c
c                                                                      c
c For a full explanation see the "Computation of radiative fluxes on   c
c staggered grid" section above.  For Frad_SW to be correctly computed c
c on w levels, the non-constant component of F_diff and F_dir, we      c
c compute taupath on w levels as a centered difference between tau     c
c values on mass levels.                                               c
c                                                                      c
c --------taupath-->--Frad_SW----------------    k = 1  (w level)      c
c        /                   \        |-dwm                            c
c -taude----------------------radht----------    k = 1  (mass level)   c
c        \                   /        |-dmw                            c
c --------taupath-->--Frad_SW----------------    k = 2  (w level)      c
c        /                   \                                         c
c -taude----------------------radht----------    k = 2  (mass level)   c
c                                                                      c
c Vince Larson changed the F variables to w levels  03 Feb 2005        c
c Michael Falk changed the loop to start at k=2 and then solved        c
c separately for k=1 so the array didn't go out of bounds.             c
c                                                                      c
c This code makes the same assumption as above that dwm=dmw.           c
c                                                                      c
c Comments by Michael Falk, 16 February 2005.                          c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if (center) then
        taupath = 0.5*taude(1)
      else
        taupath = 0.
      endif

          F_diff(1) = (-4.0/3.0) * Fs0
     .              * ( 
     .                 rp *
     .                     ( 
     .                        c1*exp(-rk*taupath)
     .                      - c2*exp(rk*taupath)
     .                     )
     .                 - beta*exp(-taupath/xi_abs) 
     .                )
          F_dir(1) = -Fs0*xi_abs*exp(-taupath/xi_abs)
          Frad_SW(1) = F_diff(1) + F_dir(1)

        do k = 2, kk+1

          if (center) then
            taupath = taupath + 0.5*(taude(k-1)+taude(k))
          else
            taupath = taupath + taude(k)
          endif


          F_diff(k) = (-4.0/3.0) * Fs0
     .              * ( 
     .                 rp *
     .                     ( 
     .                        c1*exp(-rk*taupath)
     .                      - c2*exp(rk*taupath)
     .                     )
     .                 - beta*exp(-taupath/xi_abs) 
     .                )
          F_dir(k) = -Fs0*xi_abs*exp(-taupath/xi_abs)
          Frad_SW(k) = F_diff(k) + F_dir(k)

cc Michael Falk, 28 January 2005, addition for statistics
c#ifdef STATS
c
c          if ( lstats_samp .and. if_diff > 0
c     .         .and. i.ge.imin_stats .and. i.le.imax_stats
c     .         .and. j.ge.jmin_stats .and. j.le.jmax_stats ) then
c            stats_sm%x(k,if_diff) 
c     .      = stats_sm%x(k,if_diff) + f_diff(k)
c            stats_sm%n(k,if_diff) = stats_sm%n(k,if_diff) + 1
c          endif
c          if ( lstats_samp .and. if_dir > 0
c     .         .and. i.ge.imin_stats .and. i.le.imax_stats
c     .         .and. j.ge.jmin_stats .and. j.le.jmax_stats ) then
c            stats_sm%x(k,if_dir) 
c     .      = stats_sm%x(k,if_dir) + f_dir(k)
c            stats_sm%n(k,if_dir) = stats_sm%n(k,if_dir) + 1
c          endif
c          if ( lstats_samp .and. ifrad_sw > 0
c     .         .and. i.ge.imin_stats .and. i.le.imax_stats
c     .         .and. j.ge.jmin_stats .and. j.le.jmax_stats ) then
c            stats_sw%x(k,ifrad_sw) 
c     .      = stats_sw%x(k,ifrad_sw) + Frad_SW(k)
c            stats_sw%n(k,ifrad_sw) = stats_sw%n(k,ifrad_sw) + 1
c          endif
c#endif
cc end addition for statistics
        enddo                                         ! end k loop

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c End of subroutine sunray_sw.                                         c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
