!-----------------------------------------------------------------------
! $Id: constants.F,v 1.23 2006-01-03 19:38:06 dschanen Exp $

!  Contains frequently occuring model constants
!-----------------------------------------------------------------------

      module constants

        implicit none

!       Physical constants
        double precision, parameter 
     .                  :: pi_dp = 3.14159265358979323846D0
        real, parameter :: pi    = 3.141592654
        real, parameter :: Cp    = 1004.67 ! Dry air specific heat at constant p
        real, parameter :: Lv    = 2.5e6   ! Latent heat of vaporization
        real, parameter :: R     = 287.04  ! Dry air gas constant
        real, parameter :: Rv    = 461.5   ! Water vapor gas constant

                                           ! Useful combinations of R and Rv
        real, parameter :: ep    = R/Rv         ! ep = 0.622
        real, parameter :: ep1   = (1.-ep)/ep   ! ep1 = 0.61
        real, parameter :: ep2   = 1./ep        ! ep2 = 1.61

        real, parameter :: kappa = R/Cp   ! kappa

!       Changed g to grav to make it easier to find in the code 5/25/05
!        real, parameter :: grav  = 9.80665 ! Gravitational acceleration (m/s^2)
        real, parameter :: grav  = 9.81    ! Gravitational acceleration (m/s^2)
        real, parameter :: p0    = 1.e5    ! Reference pressure

        real, parameter :: T0    = 300.     ! Reference temperature
        real, parameter :: alpha = 1. / T0  ! Buoyancy factor

        integer, parameter :: NLM    = 100

        real, parameter :: emin      = 1.e-6
        real, parameter :: vonk      = 0.4     ! von Karmans constant
        real, parameter :: rho_lw    = 1000.   ! density of liquid water; kg/m^3

!       BUGSrad Variables
        ! Nov 11 Altocu; cosine of solar zenith angle
        real, parameter :: amu0      = 0.4329
        real, parameter :: sol_const = 1.367e3 ! Solar constant

!       Logical flags
        logical, parameter   :: LH_on    = .false.  ! latin hypercube calc
        logical, parameter   :: local_kk = .false.  ! local drizzle in rain.F

        logical lbugsrad ! BUGSrad radiation scheme

!       Model parameters and constraints

        real C1, C1b, C1c, C2rt, C2thl, C2rtthl, C4, C5
        real C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc
        real C7, C7b, C7c, C8, C8b, C10
        real C11, C11b, C11c, C12
        real nu1, nu2, nu6, nu8, nu_r
        real mu, c_K, gamma_coef
        real lmin
        real taumin, taumax

!       Vince Larson added a constant to set plume widths for theta_l and rt
!       beta should vary between 0 and 3, with 1.5 the standard value

!        real, parameter :: beta = 1.5
        real beta

!       These are referenced together often enough that it made sense to 
!       make a list of them.  Note that lmin_coef is the input parameter, 
!       while the actual lmin model constant is computed from this

        character(len=10), dimension(36), parameter :: constants_list =
     .  (/"C1        ", "C1b       ", "C1c       ", "C2rt      ", 
     .    "C2thl     ", "C2rtthl   ", "C4        ", "C5        ",
     .    "C6rt      ", "C6rtb     ", "C6rtc     ", "C6thl     ",
     .    "C6thlb    ", "C6thlc    ", "C7        ", "C7b       ",
     .    "C7c       ", "C8        ", "C8b       ", "C10       ",
     .    "C11       ", "C11b      ", "C11c      ", "C12       ",
     .    "nu1       ", "nu2       ", "nu6       ", "nu8       ",
     .    "nu_r      ", "beta      ", "gamma_coef", "c_K       ", 
     .    "lmin_coef ", "taumin    ", "taumax    ", "mu        "/)


!       tolerances below which we consider moments to be zero
!       tolerance on [w, thl, rt]
!       in units of [m/s, K, kg/kg]

        real, parameter :: wtol    = 2.e-2 
        real, parameter :: thltol  = 1.e-2
        real, parameter :: rttol   = 1.e-5
        real, parameter :: qttol   = 1.e-4
        real, parameter :: sstol   = 1.e-8
        real, parameter :: difftol = 0.4

        real, parameter :: eps = 1e-10

!       Added for new mixing scheme.  
!       Without this the intrinsic size() would have to be called often

!       integer, parameter :: sclrm_dimension   = 1
        integer, parameter :: sclrm_dimension   = 2

!       Tolerance for new mixing scheme
!       Currently it's elements are equal to rttol and thltol

        real, dimension(sclrm_dimension), 
     .        parameter :: sclrtol = (/thltol, rttol/)

      end module constants

!-----------------------------------------------------------------------
! Subroutine to set-up constants
!-----------------------------------------------------------------------
      subroutine setup_constants( deltaz, C1_in, C1b_in, C1c_in,
     .                            C2rt_in, C2thl_in, C2rtthl_in,
     .                            C4_in, C5_in, 
     .                            C6rt_in, C6rtb_in, C6rtc_in, 
     .                            C6thl_in, C6thlb_in, C6thlc_in,
     .                            C7_in, C7b_in, C7c_in,
     .                            C8_in, C8b_in, C10_in, 
     .                            C11_in, C11b_in, C11c_in,
     .                            C12_in,
     .                            nu1_in, nu2_in, nu6_in, nu8_in,
     .                            nu_r_in, gamma_coef_in, mu_in,
     .                            c_K_in, beta_in, lmin_coef_in,
     .                            taumin_in, taumax_in, lbugsrad_in ) 

      use constants

      implicit none

      real, intent(in) :: deltaz
      real, intent(in) :: C1_in, C1b_in, C1c_in
      real, intent(in) :: C2rt_in, C2thl_in, C2rtthl_in
      real, intent(in) :: C4_in, C5_in
      real, intent(in) :: C6rt_in, C6rtb_in, C6rtc_in
      real, intent(in) :: C6thl_in, C6thlb_in, C6thlc_in
      real, intent(in) :: C7_in, C7b_in, C7c_in
      real, intent(in) :: C8_in, C8b_in, C10_in
      real, intent(in) :: C11_in, C11b_in, C11c_in
      real, intent(in) :: C12_in
      real, intent(in) :: nu1_in, nu2_in, nu6_in, nu8_in, nu_r_in
      real, intent(in) :: beta_in
      real, intent(in) :: lmin_coef_in, taumin_in, taumax_in
      real, intent(in) :: gamma_coef_in, mu_in, c_K_in

      logical, intent(in) :: lbugsrad_in

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
!      C12 = 1.0

!      nu1 = 20.
!      nu2 = 5.
!      nu6 = 5.
!      nu8 = 20.
!-----------------------------------------------------------------------

      C1      = C1_in
      C1b     = C1b_in
      C1c     = C1c_in
      C2rt    = C2rt_in
      C2thl   = C2thl_in
      C2rtthl = C2rtthl_in
      C4      = C4_in
      C5      = C5_in
      C6rt    = C6rt_in
      C6rtb   = C6rtb_in
      C6rtc   = C6rtc_in
      C6thl   = C6thl_in
      C6thlb  = C6thlb_in
      C6thlc  = C6thlc_in
      C7      = C7_in
      C7b     = C7b_in
      C7c     = C7c_in
      C8      = C8_in
      C8b     = C8b_in
      C10     = C10_in
      C11     = C11_in
      C11b    = C11b_in
      C11c    = C11c_in
      C12     = C12_in

      nu1   = nu1_in 
      nu2   = nu2_in 
      nu6   = nu6_in 
      nu8   = nu8_in 
      nu_r  = nu_r_in

      beta  = beta_in
    
      gamma_coef = gamma_coef_in

      taumin = taumin_in
      taumax = taumax_in
      lmin   = lmin_coef_in * deltaz

      mu  = mu_in
      c_K = c_K_in

!     The nu's are chosen for deltaz <= 40 m. Looks like they must 
!     be adjusted for larger grid spacings (Vince Larson)
!
      if ( deltaz > 40. ) then
        nu1 = nu1 * deltaz / 40.
        nu2 = nu2 * deltaz / 40.
        nu6 = nu6 * deltaz / 40.
        nu8 = nu8 * deltaz / 40.
      endif

!     There should be a different formula for determining nu_r for
!     different sized grid spacings.  For DYCOMS2 RF02, nu_r is set
!     to 5.0 for the high-resolution 10 m. grid spacing and to
!     25.0 for the low-resolution 100 m. grid spacing.  The following
!     equation allows for both of those parameters.  Brian.
!
      if ( deltaz > 20. ) then
        nu_r = nu_r * deltaz / 20.
      endif

!     Constraints on mixing length
!     Now in error.in, now 4.0 is the minimum

!      lmin = 20.     ! Minimum value

!     For DYCOMS
!      lmin = 0.5 * deltaz

!     Constraints on dissipation time (now in error.in)

!      taumin =   90.
!      taumax = 3600.
      lbugsrad = lbugsrad_in ! dschanen

      return
      end subroutine setup_constants
