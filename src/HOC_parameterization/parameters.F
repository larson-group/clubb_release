!-----------------------------------------------------------------------
! $Id: parameters.F,v 1.8 2007-08-02 22:51:40 griffinb Exp $
        module parameters 

!       Description:
!       Contains tuneable model parameters.

!       References:
!       None
!-----------------------------------------------------------------------

        use param_index

        implicit none

        ! Default to private
        private

        public :: setup_parameters, read_parameters, read_param_spread

        ! Model parameters and constraints setup in the namelists
        real, public :: 
     .  T0,      ! Reference temperature (usually 300)  [K]
     .  ts_nudge ! Timescale of u/v nudging             [s]

!$omp   threadprivate(T0,ts_nudge)

        ! Model constant parameters 
        real, public ::
     .  C1,         ! Low Skewness in C1 Skewness Function.
     .  C1b,        ! High Skewness in C1 Skewness Function.
     .  C1c,        ! Degree of Slope of C1 Skewness Function.
     .  C2,         ! Low Skewness in C2 Skewness Function.
     .  C2rt,       ! C2 coefficient for the rtp2_dp1 term.
     .  C2thl,      ! C2 coefficient for the thlp2_dp1 term.
     .  C2rtthl,    ! C2 coefficient for the rtpthlp_dp1 term.
     .  C2b,        ! High Skewness in C2 Skewness Function.  
     .  C2c,        ! Degree of Slope of C2 Skewness Function.
     .  C4,         ! Used only when lKhm_aniso is true.
     .  C5,         ! Coefficient in pressure terms in the w'^2 eqn.
     .  C6rt,       ! Low Skewness in C6rt Skewness Function.
     .  C6rtb,      ! High Skewness in C6rt Skewness Function.
     .  C6rtc,      ! Degree of Slope of C6rt Skewness Function.
     .  C6thl,      ! Low Skewness in C6thl Skewness Function.
     .  C6thlb,     ! High Skewness in C6thl Skewness Function.
     .  C6thlc,     ! Degree of Slope of C6thl Skewness Function.
     .  C7,         ! Low Skewness in C7 Skewness Function.
     .  C7b,        ! High Skewness in C7 Skewness Function.
     .  C7c,        ! Degree of Slope of C7 Skewness Function.
     .  C8,         ! Coefficient #1 in C8 Skewness Equation.
     .  C8b,        ! Coefficient #2 in C8 Skewness Equation.  From Golaz.
     .  C10,        ! Currently Not Used in the Model.
     .  C11,        ! Low Skewness in C11 Skewness Function.
     .  C11b,       ! High Skewness in C11 Skewness Function.  From Golaz.
     .  C11c,       ! Degree of Slope of C11 Skewness Function.  From Golaz.
     .  C12,        ! Constant in w'^3 Crank-Nicholson diffusional term.
     .  C13,        ! Not currently used in model.
     .  C14,        ! Constant for u'^2 and v'^2 terms.
     .  nu1,        ! Previously determined eddy viscosity.
     .  nu2,        ! Previously determined eddy viscosity.
     .  nu6,        ! Previously determined eddy viscosity.
     .  nu8,        ! Previously determined eddy viscosity.
     .  nu9,        ! Previously determined eddy viscosity.
     .  nu_r,       ! Background eddy diffusivity for hydrometeors.
     .  gamma_coef, ! Low Skewness in gamma coefficient Skewness Function.
     .  gamma_coefb,! High Skewness in gamma coefficient Skewness Function.
     .  gamma_coefc,! Degree of Slope of gamma coefficient Skewness Function.
     .  c_K,        ! Constant C_mu^(1/4) in Duynkerke & Driedonks 1987.
     .  c_Krrm,     ! Constant for eddy diffusivity equation for hydrometeors.
     .  mu,         ! Fractional entrainment rate per unit altitude.
     .  taumin,     ! Previously determined value.
     .  taumax,     ! Previously determined value.
     .  lmin        ! Minimum value for the length scale.

!$omp   threadprivate(C1, C1b, C1c, C2, C2b, C2c)
!$omp   threadprivate(C2rt, C2thl, C2rtthl, C4, C5, C6rt, C6rtb, C6rtc)
!$omp   threadprivate(C6thl, C6thlb, C6thlc)
!$omp   threadprivate(C7, C7b, C7c, C8, C8b, C10, C11, C11b, C11c, C12)
!$omp   threadprivate(C13, C14)
!$omp   threadprivate(nu1, nu2, nu6, nu8, nu9, nu_r, beta)
!$omp   threadprivate(gamma_coef, gamma_coefb, gamma_coefc)
!$omp   threadprivate(c_K, c_Krrm, taumin, taumax, mu, lmin)

!       Vince Larson added a constant to set plume widths for theta_l and rt
!       beta should vary between 0 and 3, with 1.5 the standard value

        real, public :: beta

!$omp   threadprivate(beta)

        real :: lmin_coef ! Coefficient of lmin

!$omp   threadprivate(lmin_coef)

       ! Since we lack a devious way to do this just once, this namelist 
       ! must be changed as well when a new parameter is added.
        namelist /initvars/ 
     .    C1, C1b, C1c, C2, C2b, C2c, 
     .    C2rt, C2thl, C2rtthl, C4, C5,
     .    C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc,
     .    C7, C7b, C7c, C8, C8b, C10, C11, C11b, C11c,
     .    C12, C13, C14, nu1, nu2, nu6, nu8, nu9, nu_r, beta, 
     .    gamma_coef, gamma_coefb, gamma_coefc,
     .    c_K, c_Krrm, lmin_coef, taumin, taumax, mu

!       These are referenced together often enough that it made sense to 
!       make a list of them.  Note that lmin_coef is the input parameter, 
!       while the actual lmin model constant is computed from this

        character(len=11), dimension(nparams), parameter, public :: 
     .  params_list =
     .  (/"C1         ", "C1b        ", "C1c        ", "C2         ", 
     .    "C2b        ", "C2c        ", "C2rt       ", "C2thl      ",
     .    "C2rtthl    ", "C4         ", "C5         ", "C6rt       ",
     .    "C6rtb      ", "C6rtc      ", "C6thl      ", "C6thlb     ", 
     .    "C6thlc     ", "C7         ", "C7b        ", "C7c        ", 
     .    "C8         ", "C8b        ", "C10        ", "C11        ", 
     .    "C11b       ", "C11c       ", "C12        ", "C13        ",
     .    "nu1        ", "nu2        ", "nu6        ", "nu8        ",
     .    "nu_r       ", "gamma_coef ", "gamma_coefb", "gamma_coefc", 
     .    "mu         ", "c_K        ", "c_Krrm     ", "beta       ",
     .    "lmin_coef  ", "taumin     ", "taumax     ", "C14        ",
     .    "nu9        "/)

        contains

!-----------------------------------------------------------------------
        subroutine setup_parameters
     .             ( deltaz, T0_in, ts_nudge_in, params, isValid ) 
!       Description:
!       Subroutine to setup model parameters

!       References:
!       None
!-----------------------------------------------------------------------

        use constants, only: fstderr

        implicit none

        ! Input Variables
        real, intent(in) :: 
     .  deltaz,     ! Change per height level      [m]
     .  T0_in,      ! Ref. temperature             [K]
     .  ts_nudge_in ! Timescale for u/v nudging    [K]

        real, intent(in), dimension(nparams) :: 
     .  params  ! Tuneable model parameters      [-]

        ! Output Variables
        logical, intent(out) :: isValid ! Error condition

        call unpack_parameters
     .       ( params,
     .         C1, C1b, C1c, C2, C2b, C2c, C2rt, C2thl, C2rtthl,
     .         C4, C5, C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc,
     .         C7, C7b, C7c, C8, C8b, C10,
     .         C11, C11b, C11c, C12, C13, C14,
     .         nu1, nu2, nu6, nu8, nu9, nu_r,
     .         gamma_coef, gamma_coefb, gamma_coefc,
     .         mu, c_K, c_Krrm, beta, lmin_coef,
     .         taumin, taumax )

        ! It was decided after some experimentation, that the best
        ! way to produce grid independent results is to set lmin to be
        ! some fixed value at the surface. -dschanen 21 May 2007
        ! lmin = lmin_coef * deltaz  ! Old
        lmin = lmin_coef * 40.0 ! New fixed value

        T0       = T0_in
        ts_nudge = ts_nudge_in

        ! The nu's are chosen for deltaz <= 40 m. Looks like they must 
        ! be adjusted for larger grid spacings (Vince Larson)
 
        if ( deltaz > 40.0 ) then
          nu1 = nu1 * deltaz / 40.0
          nu2 = nu2 * deltaz / 40.0
          nu6 = nu6 * deltaz / 40.0
          nu8 = nu8 * deltaz / 40.0
          nu9 = nu9 * deltaz / 40.0
        end if

       ! There should be a different formula for determining nu_r for
       ! different sized grid spacings.  For DYCOMS2 RF02, nu_r is set
       ! to 5.0 for the high-resolution 10 m. grid spacing and to
       ! 25.0 for the low-resolution 100 m. grid spacing.  The following
       ! equation allows for both of those parameters.  Brian.
 
       if ( deltaz > 20.0 ) then
         nu_r = nu_r * deltaz / 20.0
       end if

       ! Sanity check
       if ( beta < 0.0 .or. beta > 3.0 ) then
         ! Constraints on beta
         write(fstderr,*) "beta= ", beta
         write(fstderr,*) "beta cannot be < 0 or > 3"
         isValid = .false.
       else if ( lmin < 4.0 ) then
         ! Constraints on mixing length
         write(fstderr,*) "lmin= ", lmin
         write(fstderr,*) "lmin is < 4.0"
         isValid = .false.
       else
         isValid = .true.
       end if

!      write(*,nml=initvars) ! %% debug

       return
       end subroutine setup_parameters

!-----------------------------------------------------------------------
       subroutine read_parameters( iunit, filename, params )

!      Description:
!      Read a namelist containing the model parameters

!      References:
!      None
!-----------------------------------------------------------------------

       use param_index

       implicit none

       ! Input variables
       integer, intent(in) :: iunit

       character(len=*), intent(in) :: filename

       ! Output variables
       real, intent(out), dimension(nparams) :: params

       ! If the filename is empty, assume we're using a `working' set of
       ! parameters that are set statically here (handy for host models).
       if ( filename == "" ) then
         C1          = 2.5
         C1b         = 2.5
         C1c         = 1.0
         C2rt        = 1.0
         C2thl       = 1.0
         C2rtthl     = 2.0
         C2          = 1.3
         C2b         = 1.3
         C2c         = 5.0
         C4          = 5.2
         C5          = 0.3 
         C6rt        = 6.0 
         C6rtb       = 6.0
         C6rtc       = 1.0
         C6thl       = 6.0 
         C6thlb      = 6.0 
         C6thlc      = 1.0 
         C7          = 0.1 
         C7b         = 0.8 
         C7c         = 0.5 
         C8          = 3.0  
         C8b         = 0.005
         C10         = 3.3
         C11         = 0.75
         C11b        = 0.35
         C11c        = 0.5
         C12         = 1.0
         C13         = 0.1
         C14         = 1.0
         nu1         = 20.0
         nu2         = 5.0
         nu6         = 5.0
         nu8         = 20.0
         nu9         = 20.0
         nu_r        = 3.0
         beta        = 1.75
         gamma_coef  = 0.32
         gamma_coefb = 0.32
         gamma_coefc = 5.0
         c_K         = 0.548
         c_Krrm      = 0.075
         taumin      = 90.0
         taumax      = 3600.0
         lmin_coef   = 0.5
         mu          = 6.000E-4

       else
         ! Read the namelist
         open(unit=iunit, file=filename, status='old', action='read')

         read(unit=iunit, nml=initvars)

         close(unit=iunit)

       end if

       ! Put the variables in the output array 
       call pack_parameters
     .      ( C1, C1b, C1c, C2, C2b, C2c, C2rt, C2thl, C2rtthl,
     .        C4, C5, C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc,
     .        C7, C7b, C7c, C8, C8b, C10,
     .        C11, C11b, C11c, C12, C13, C14,
     .        nu1, nu2, nu6, nu8, nu9, nu_r,
     .        gamma_coef, gamma_coefb, gamma_coefc,
     .        mu, c_K, c_Krrm, beta, lmin_coef,
     .        taumin, taumax, params )

       return
       end subroutine read_parameters

!-----------------------------------------------------------------------
       subroutine read_param_spread
     .            ( iunit, filename, nindex, param_spread, ndim )

!      Description:
!      Read a namelist containing the amount to vary model parameters.
!      Used by the downhill simplex / simulated annealing algorithm.

!      References:
!      None
!-----------------------------------------------------------------------

       use param_index

       implicit none

       ! Input variables
       integer, intent(in) :: iunit

       character(len=*), intent(in) :: filename

       ! Output variables

       ! An array of array indices (i.e. which elements of the array `params'
       ! are contained within the simplex and the spread variable)
       integer, intent(out), dimension(nparams) :: nindex  

       real, intent(out), dimension(nparams) :: 
     . param_spread  ! Amount to vary the parameter in the initial simplex
       
       integer, intent(out) :: ndim  ! Dimension of the init simplex

       ! Local variables
       integer :: i

       ! Amount to change each parameter for the initial simplex
       ! This MUST to match the initvars namelist if parameters are added!
        namelist /initspread/ 
     .    C1, C1b, C1c, C2, C2b, C2c, 
     .    C2rt, C2thl, C2rtthl, C4, C5,
     .    C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc,
     .    C7, C7b, C7c, C8, C8b, C10, C11, C11b, C11c,
     .    C12, C13, C14, nu1, nu2, nu6, nu8, nu9, nu_r, beta, 
     .    gamma_coef, gamma_coefb, gamma_coefc,
     .    c_K, c_Krrm, lmin_coef, taumin, taumax, mu

       ! Read the namelist
       open(unit=iunit, file=filename, status='old', action='read')

       read(unit=iunit, nml=initspread)

       close(unit=iunit)

       ! Put the variables in the output array 
       call pack_parameters
     .      ( C1, C1b, C1c, C2, C2b, C2c, C2rt, C2thl, C2rtthl,
     .        C4, C5, C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc,
     .        C7, C7b, C7c, C8, C8b, C10,
     .        C11, C11b, C11c, C12, C13, C14,
     .        nu1, nu2, nu6, nu8, nu9, nu_r,
     .        gamma_coef, gamma_coefb, gamma_coefc,
     .        mu, c_K, c_Krrm, beta, lmin_coef,
     .        taumin, taumax, param_spread )

       ! Initialize to zero
       nindex(1:nparams) = 0
       ndim = 0

       ! Determine how many variables are being changed
       do i = 1, nparams, 1

         if ( param_spread(i) /= 0.0 ) then
           ndim = ndim + 1   ! Increase the total
           nindex(ndim) = i  ! Set the next array index
         end if

       end do

       return
       end subroutine read_param_spread

       end module parameters

!-----------------------------------------------------------------------
       subroutine pack_parameters
     .            ( C1, C1b, C1c, C2, C2b, C2c, C2rt, C2thl, C2rtthl,
     .              C4, C5, C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc,
     .              C7, C7b, C7c, C8, C8b, C10,
     .              C11, C11b, C11c, C12, C13, C14,
     .              nu1, nu2, nu6, nu8, nu9, nu_r,
     .              gamma_coef, gamma_coefb, gamma_coefc,
     .              mu, c_K, c_Krrm, beta, lmin_coef,
     .              taumin, taumax, params )
!       Description:
!       Takes the list of scalar variables and puts them into a 1D
!       vector.  Here for the purposes of keeping the code generalized
!       when new variables are added.

!       References:
!       None
!-----------------------------------------------------------------------

        use param_index

        implicit none

        ! Input variables
        real, intent(in) ::
     .  C1, C1b, C1c, C2, C2b, C2c, C2rt, C2thl, C2rtthl,
     .  C4, C5, C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc,
     .  C7, C7b, C7c, C8, C8b, C10,
     .  C11, C11b, C11c, C12, C13, C14,
     .  nu1, nu2, nu6, nu8, nu9, nu_r,
     .  gamma_coef, gamma_coefb, gamma_coefc,
     .  mu, c_K, c_Krrm, beta, lmin_coef,
     .  taumin, taumax

        ! Output variables
        real, intent(out), dimension(nparams) :: params

        params(iC1)      = C1
        params(iC1b)     = C1b
        params(iC1c)     = C1c
        params(iC2)      = C2
        params(iC2b)     = C2b
        params(iC2c)     = C2c
        params(iC2rt)    = C2rt
        params(iC2thl)   = C2thl
        params(iC2rtthl) = C2rtthl
        params(iC4)      = C4
        params(iC5)      = C5
        params(iC6rt)    = C6rt
        params(iC6rtb)   = C6rtb
        params(iC6rtc)   = C6rtc
        params(iC6thl)   = C6thl
        params(iC6thlb)  = C6thlb
        params(iC6thlc)  = C6thlc
        params(iC7)      = C7
        params(iC7b)     = C7b
        params(iC7c)     = C7c
        params(iC8)      = C8
        params(iC8b)     = C8b
        params(iC10)     = C10
        params(iC11)     = C11
        params(iC11b)    = C11b
        params(iC11c)    = C11c
        params(iC12)     = C12
        params(iC13)     = C13
        params(iC14)     = C14

        params(inu1)  = nu1
        params(inu2)  = nu2
        params(inu6)  = nu6
        params(inu8)  = nu8
        params(inu9)  = nu9
        params(inu_r) = nu_r

        params(igamma_coef)  = gamma_coef
        params(igamma_coefb) = gamma_coefb
        params(igamma_coefc) = gamma_coefc

        params(imu) = mu

        params(ic_K) = c_K 

        params(ic_Krrm) = c_Krrm

        params(ibeta) = beta

        params(ilmin_coef) = lmin_coef

        params(itaumin) = taumin
        params(itaumax) = taumax

        return
        end subroutine pack_parameters

!-----------------------------------------------------------------------
        subroutine unpack_parameters
     .             ( params, 
     .               C1, C1b, C1c, C2, C2b, C2c, C2rt, C2thl, C2rtthl,
     .               C4, C5, C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc,
     .               C7, C7b, C7c, C8, C8b, C10,
     .               C11, C11b, C11c, C12, C13, C14,
     .               nu1, nu2, nu6, nu8, nu9, nu_r,
     .               gamma_coef, gamma_coefb, gamma_coefc,
     .               mu, c_K, c_Krrm, beta, lmin_coef,
     .               taumin, taumax )

!       Description:
!       Takes the 1D vector and returns the list of scalar variables.
!       Here for the purposes of keeping the code generalized
!       when new variables are added.

!       References:
!       None
!-----------------------------------------------------------------------

        use param_index

        implicit none


        ! Input variables
        real, intent(in), dimension(nparams) :: params

        ! Output variables
        real, intent(out) ::
     .  C1, C1b, C1c, C2, C2b, C2c, C2rt, C2thl, C2rtthl,
     .  C4, C5, C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc,
     .  C7, C7b, C7c, C8, C8b, C10,
     .  C11, C11b, C11c, C12, C13, C14,
     .  nu1, nu2, nu6, nu8, nu9, nu_r,
     .  gamma_coef, gamma_coefb, gamma_coefc,
     .  mu, c_K, c_Krrm, beta, lmin_coef,
     .  taumin, taumax

        C1      = params(iC1)
        C1b     = params(iC1b)
        C1c     = params(iC1c)
        C2      = params(iC2)
        C2b     = params(iC2b)
        C2c     = params(iC2c)
        C2rt    = params(iC2rt)
        C2thl   = params(iC2thl)
        C2rtthl = params(iC2rtthl)
        C4      = params(iC4)
        C5      = params(iC5)
        C6rt    = params(iC6rt)
        C6rtb   = params(iC6rtb)
        C6rtc   = params(iC6rtc)
        C6thl   = params(iC6thl)
        C6thlb  = params(iC6thlb)
        C6thlc  = params(iC6thlc)
        C7      = params(iC7)
        C7b     = params(iC7b)
        C7c     = params(iC7c)
        C8      = params(iC8)
        C8b     = params(iC8b)
        C10     = params(iC10)
        C11     = params(iC11)
        C11b    = params(iC11b)
        C11c    = params(iC11c)
        C12     = params(iC12)
        C13     = params(iC13)
        C14     = params(iC14)

        nu1  = params(inu1)
        nu2  = params(inu2)
        nu6  = params(inu6)
        nu8  = params(inu8)
        nu9  = params(inu9)
        nu_r = params(inu_r)

        gamma_coef  = params(igamma_coef)
        gamma_coefb = params(igamma_coefb)
        gamma_coefc = params(igamma_coefc)

        mu = params(imu)

        c_K = params(ic_K)

        c_Krrm = params(ic_Krrm)

        beta = params(ibeta)

        lmin_coef = params(ilmin_coef)

        taumin = params(itaumin)
        taumax = params(itaumax)

        return
        end subroutine unpack_parameters
