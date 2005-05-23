! Subroutine to set-up constants

      module constants

!       Physical constants

!        real, parameter :: pi = 3.141592654
!       This is the accuracy used in math.h;  Not sure if it improves anything
        real, parameter :: pi = 3.14159265358979323846
        real, parameter :: Cp = 1004.67   ! Dry air specific heat at constant p
        real, parameter :: Lv = 2.5e6     ! Latent heat of vaporization
        real, parameter :: R  = 287.04    ! Dry air gas constant
        real, parameter :: Rv = 461.5     ! Water vapor gas constant

                                          ! Useful combinations of R and Rv
        real, parameter :: ep  = R/Rv         ! ep = 0.622
        real, parameter :: ep1 = (1.-ep)/ep   ! ep1 = 0.61
        real, parameter :: ep2 = 1./ep        ! ep2 = 1.61

        real, parameter :: kappa = R/Cp   ! kappa

        real, parameter :: g     = 9.81   ! Gravitational acceleration (m/s^2)
        real, parameter :: p0    = 1.e5   ! Reference pressure

        real, parameter :: T0    = 300.     ! Reference temperature
        real, parameter :: alpha = 1. / T0  ! Buoyancy factor

        integer, parameter :: NLM  = 100
        real, parameter    :: vonk = 0.4

        real, parameter :: emin = 1.e-6

!       Model parameters and constraints

        real C1, C2, C4, C5, C6, C7, C8, C10, C11
        real nu1, nu2, nu6, nu8
        real lmin
        real taumin, taumax

!       Vince Larson added a constant to set plume widths for theta_l and rt
!       beta should vary between 0 and 3, with 1.5 the standard value

!        real, parameter :: beta = 1.5
        real beta

!       tolerances below which we consider moments to be zero
!       tolerance on [w, thl, rt]
!       in units of [m/s, K, kg/kg]

        real, parameter :: wtol    = 2.e-2 
        real, parameter :: thltol  = 1.e-2
        real, parameter :: rttol   = 1.e-4
        real, parameter :: qttol   = 1.e-4
        real, parameter :: sstol   = 1.e-8
        real, parameter :: difftol = 0.4

      end module constants

      subroutine setup_constants( deltaz, C1_in, C2_in, C4_in, C5_in, 
     .                            C6_in, C7_in, C8_in, C10_in, C11_in, 
     .                            nu1_in, nu2_in, nu6_in, nu8_in, 
     .                            beta_in, lmin_coef_in, taumin_in, 
     .                            taumax_in ) 

      use constants

      implicit none

      real, intent(in) :: deltaz
      real, intent(in) :: C1_in, C2_in, C4_in, C5_in, C6_in, C7_in, 
     .                    C8_in, C10_in, C11_in
      real, intent(in) :: nu1_in, nu2_in, nu6_in, nu8_in, beta_in
      real, intent(in) :: lmin_coef_in, taumin_in, taumax_in

!----- These now occur in error.in ------------------------------------
!      C1 = 1.7
!      C2 = 0.8 * 1.3
!      C4 = 4.5
!      C5 = 0.
!      C6 = 4.85 * 2.0
!      C7 = 0.8
!      C8 = 0.7 * 0.6 * 6.5
!      C10 = 3. * C2
!      C11 = 0.2

!      nu1 = 20.
!      nu2 = 5.
!      nu6 = 5.
!      nu8 = 20.
!-----------------------------------------------------------------------

      C1 = C1_in
      C2 = C2_in
      C4 = C4_in
      C5 = C5_in
      C6 = C6_in
      C7 = C7_in
      C8 = C8_in
      C10 = C10_in
      C11 = C11_in

      nu1 = nu1_in 
      nu2 = nu2_in 
      nu6 = nu6_in 
      nu8 = nu8_in 
      beta = beta_in
    
      taumin = taumin_in
      taumax = taumax_in
      lmin = lmin_coef_in * deltaz

!     The nu's are chosen for deltaz <= 40 m. Looks like they must 
!     be adjusted for larger grid spacings (Vince Larson)
!
      if ( deltaz > 40. ) then
        nu1 = nu1 * deltaz / 40.
        nu2 = nu2 * deltaz / 40.
        nu6 = nu6 * deltaz / 40.
        nu8 = nu8 * deltaz / 40.
      endif

!     Constraints on mixing length
!     Now in error.in, now 4.0 is the minimum
!      lmin = 20.     ! Minimum value

!     For DYCOMS
!      lmin = 0.5 * deltaz

!     Constraints on dissipation time

!      taumin =   90.
!      taumax = 3600.

      return
      end subroutine setup_constants
