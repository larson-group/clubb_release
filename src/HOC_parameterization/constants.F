!-----------------------------------------------------------------------
! $Id: constants.F,v 1.6 2006-11-02 01:43:41 dschanen Exp $

!  Contains frequently occuring model constants
!-----------------------------------------------------------------------

        module constants

        implicit none

        ! Fortran file unit I/O constants
        integer, parameter :: 
     .  fstderr = 0, fstdin = 5, fstdout = 6

        ! Mathematical Constants
        double precision, parameter :: 
     .  pi_dp = 3.14159265358979323846d0
        real, parameter :: 
     .  pi = 3.141592654 ! The ratio of radii to their circumference

        ! Physical constants
        real, parameter :: 
     .  Cp = 1004.67, ! Dry air specific heat at constant p [J/kg/K]
     .  Lv = 2.5e6,   ! Latent heat of vaporization         [J/kg]
     .  Ls = 2.834e6, ! Latent heat of sublimation          [J/kg]
     .  Rd = 287.04,  ! Dry air gas constant                [J/kg/K]
     .  Rv = 461.5    ! Water vapor gas constant            [J/kg/K]

        ! Useful combinations of R and Rv
        real, parameter :: 
     .  ep  = Rd / Rv,    ! ep = 0.622
     .  ep1 = (1.0-ep)/ep,! ep1 = 0.61
     .  ep2 = 1.0/ep      ! ep2 = 1.61

        real, parameter :: 
     .  kappa = Rd / Cp   ! kappa

        ! Changed g to grav to make it easier to find in the code 5/25/05
!       real, parameter :: grav  = 9.80665 ! Gravitational acceleration [m/s^2]
        real, parameter :: 
     .  grav = 9.81,   ! Gravitational acceleration [m/s^2]
     .  p0   = 1.0e5   ! Reference pressure [P]

        ! Von Karman's constant
        ! Constant of the logarithmic wind profile in the surface layer 
        real, parameter :: 
     .  vonk   = 0.4,    ! Accepted value 0.40 (+/-) 0.01
     .  rho_lw = 1000.0, ! Density of liquid water [kg/m^3]
     .  T0     = 300.0   ! Reference temperature [K]

        real, parameter :: emin = 1.0e-6

        ! BUGSrad Constants
        real, parameter :: 
     .  sol_const = 1.367e3 ! Solar constant


        ! Model parameters and constraints

        real :: 
     .  fcor ! Coriolis parameter

        real ::
     .  C1, C1b, C1c, C2, C2b, C2c, C4, C5,
     .  C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc,
     .  C7, C7b, C7c, C8, C8b, C10,
     .  C11, C11b, C11c, C12,
     .  nu1, nu2, nu6, nu8, nu_r,
     .  mu, c_K, 
     .  gamma_coef, gamma_coefb, gamma_coefc,
     .  lmin,
     .  taumin, taumax

!       Vince Larson added a constant to set plume widths for theta_l and rt
!       beta should vary between 0 and 3, with 1.5 the standard value

!       real, parameter :: beta = 1.5
        real :: beta

!       These are referenced together often enough that it made sense to 
!       make a list of them.  Note that lmin_coef is the input parameter, 
!       while the actual lmin model constant is computed from this

        character(len=11), dimension(38), parameter :: constants_list =
     .  (/"C1         ", "C1b        ", "C1c        ", "C2         ", 
     .    "C2b        ", "C2c        ", "C4         ", "C5         ",
     .    "C6rt       ", "C6rtb      ", "C6rtc      ", "C6thl      ",
     .    "C6thlb     ", "C6thlc     ", "C7         ", "C7b        ",
     .    "C7c        ", "C8         ", "C8b        ", "C10        ",
     .    "C11        ", "C11b       ", "C11c       ", "C12        ",
     .    "nu1        ", "nu2        ", "nu6        ", "nu8        ",
     .    "nu_r       ", "beta       ", "gamma_coef ",
     .    "gamma_coefb", "gamma_coefc", "c_K        ", 
     .    "lmin_coef  ", "taumin     ", "taumax     ", "mu         "/)


        ! Tolerances below which we consider moments to be zero
        real, parameter :: 
     .  wtol    = 2.e-2, ! [m/s]
     .  thltol  = 1.e-2, ! [K]
     .  rttol   = 1.e-5, ! [kg/kg]
     .  qttol   = 1.e-4, ! [?]
     .  sstol   = 1.e-8, ! [kg/kg]
     .  difftol = 0.4    ! [?]

        real, parameter :: 
     .  eps = 1.0e-10 ! Small value to prevent a divide by zero

        ! Added for new mixing scheme.  
        ! Without this the intrinsic size() would have to be called often
        integer, parameter :: sclr_dim = 2

        ! Tolerance for new mixing scheme
        ! Currently it's elements are equal to rttol and thltol
        real, dimension(sclr_dim), parameter :: 
     .  sclrtol = (/thltol, rttol/)

        end module constants

!-----------------------------------------------------------------------
        subroutine setup_constants
     .             ( deltaz, C1_in, C1b_in, C1c_in, 
     .               C2_in, C2b_in, C2c_in, C4_in, C5_in, 
     .               C6rt_in, C6rtb_in, C6rtc_in, 
     .               C6thl_in, C6thlb_in, C6thlc_in,
     .               C7_in, C7b_in, C7c_in,
     .               C8_in, C8b_in, C10_in, 
     .               C11_in, C11b_in, C11c_in, C12_in,
     .               nu1_in, nu2_in, nu6_in, nu8_in, nu_r_in,
     .               gamma_coef_in, gamma_coefb_in, gamma_coefc_in, 
     .               mu_in, c_K_in, beta_in, lmin_coef_in,
     .               taumin_in, taumax_in, fcor_in ) 
!       Description:
!       Subroutine to set-up constants

!       References:
!       None
!-----------------------------------------------------------------------

        use constants

        implicit none

        ! Input Variables
        real, intent(in) :: 
     .  deltaz ! Change per height level [m]

        ! Constant parameters
        real, intent(in) :: 
     .  C1_in, C1b_in, C1c_in,
     .  C2_in, C2b_in, C2c_in, C4_in, C5_in,
     .  C6rt_in, C6rtb_in, C6rtc_in,
     .  C6thl_in, C6thlb_in, C6thlc_in,
     .  C7_in, C7b_in, C7c_in,
     .  C8_in, C8b_in, C10_in,
     .  C11_in, C11b_in, C11c_in, C12_in,
     .  nu1_in, nu2_in, nu6_in, nu8_in, nu_r_in,
     .  gamma_coef_in, gamma_coefb_in, gamma_coefc_in,
     .  mu_in, c_K_in, beta_in, lmin_coef_in, 
     .  taumin_in, taumax_in, fcor_in

       ! These now occur as parameters in a namelist
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

!      nu1 = 20.0
!      nu2 = 5.0
!      nu6 = 5.0
!      nu8 = 20.0

       C1     = C1_in
       C1b    = C1b_in
       C1c    = C1c_in
       C2     = C2_in
       C2b    = C2b_in
       C2c    = C2c_in
       C4     = C4_in
       C5     = C5_in
       C6rt   = C6rt_in
       C6rtb  = C6rtb_in
       C6rtc  = C6rtc_in
       C6thl  = C6thl_in
       C6thlb = C6thlb_in
       C6thlc = C6thlc_in
       C7     = C7_in
       C7b    = C7b_in
       C7c    = C7c_in
       C8     = C8_in
       C8b    = C8b_in
       C10    = C10_in
       C11    = C11_in
       C11b   = C11b_in
       C11c   = C11c_in
       C12    = C12_in

       nu1  = nu1_in 
       nu2  = nu2_in 
       nu6  = nu6_in 
       nu8  = nu8_in 
       nu_r = nu_r_in

       beta = beta_in
    
       gamma_coef  = gamma_coef_in
       gamma_coefb = gamma_coefb_in
       gamma_coefc = gamma_coefc_in

       taumin = taumin_in
       taumax = taumax_in
       lmin   = lmin_coef_in * deltaz

       mu  = mu_in
       c_K = c_K_in

       fcor = fcor_in ! Not currently used 

       ! The nu's are chosen for deltaz <= 40 m. Looks like they must 
       ! be adjusted for larger grid spacings (Vince Larson)
 
       if ( deltaz > 40.0 ) then
         nu1 = nu1 * deltaz / 40.0
         nu2 = nu2 * deltaz / 40.0
         nu6 = nu6 * deltaz / 40.0
         nu8 = nu8 * deltaz / 40.0
       end if

      ! There should be a different formula for determining nu_r for
      ! different sized grid spacings.  For DYCOMS2 RF02, nu_r is set
      ! to 5.0 for the high-resolution 10 m. grid spacing and to
      ! 25.0 for the low-resolution 100 m. grid spacing.  The following
      ! equation allows for both of those parameters.  Brian.
 
      if ( deltaz > 20.0 ) then
        nu_r = nu_r * deltaz / 20.0
      end if

      ! Constraints on mixing length
      ! Now in hoc interface. A value of 4.0 is the current minimum

!     lmin = 20.0 ! Minimum value

      ! For DYCOMS II RF01
!     lmin = 0.5 * deltaz

      ! Constraints on dissipation time

!     taumin =   90.0
!     taumax = 3600.0

      return
      end subroutine setup_constants
