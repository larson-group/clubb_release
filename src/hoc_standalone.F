!----------------------------------------------------------------------
! $Id: hoc_standalone.F,v 1.15 2007-02-23 01:25:47 dschanen Exp $

        program hoc_standalone
        use hoc, only: hoc_model
!       Description:
!       This is essentially a minimalist frontend for HOC
!
!----------------------------------------------------------------------

        implicit none

        ! External
        intrinsic :: trim

        ! Constant Variables
        integer, parameter :: funit = 10

        character(len=13), parameter :: fstandalone = "standalone.in"
 
        ! Model constant parameters 
        real ::
     .  C1,         ! Low Skewness in C1 Skewness Function.
     .  C1b,        ! High Skewness in C1 Skewness Function.
     .  C1c,        ! Degree of Slope of C1 Skewness Function.
     .  C2,         ! Low Skewness in C2 Skewness Function.
     .  C2b,        ! High Skewness in C2 Skewness Function.  
     .  C2c,        ! Degree of Slope of C2 Skewness Function.
     .  C2rt,       ! C2 used in rtp2 diffusion term.
     .  C2thl,      ! C2 used in thlp2 diffusion term.
     .  C2rtthl,    ! C2 used in rtpthlp diffusion term.
     .  C4,         ! Currently Not Used in the Model.
     .  C5,         ! Constant from Chris Golaz.
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
     .  C12,        ! Constant from Chris Golaz.
     .  C13,        ! Constant for up2_plus_vp2.
     .  nu1,        ! Previously determined eddy viscosity.
     .  nu2,        ! Previously determined eddy viscosity.
     .  nu6,        ! Previously determined eddy viscosity.
     .  nu8,        ! Previously determined eddy viscosity.
     .  nu_r,       ! Tuned to LES from Surface Value of DYCOMS2 RF02 DS.
     .  beta,       ! Constant from Chris Golaz.
     .  gamma_coef, ! Low Skewness in gamma coefficient Skewness Function.
     .  gamma_coefb,! High Skewness in gamma coefficient Skewness Function.
     .  gamma_coefc,! Degree of Slope of gamma coefficient Skewness Function.
     .  c_K,        ! Constant from Chris Golaz.
     .  mu,         ! Constant from Chris Golaz.
     .  taumin,     ! Previously determined value.
     .  taumax,     ! Previously determined value.
     .  lmin_coef   ! Previously determined value.

        ! Run information
        character(len=50) :: run_file

        logical :: stdout

        ! Internal variables
        logical :: isValid 

        ! Namelists 
        namelist /initvars/ C1, C1b, C1c, C2, C2b, C2c, 
     .    C2rt, C2thl, C2rtthl, C4, C5,
     .    C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc,
     .    C7, C7b, C7c, C8, C8b, C10, C11, C11b, C11c,
     .    C12, C13, nu1, nu2, nu6, nu8, nu_r, beta, 
     .    gamma_coef, gamma_coefb, gamma_coefc,
     .    c_K, lmin_coef, taumin, taumax, mu

        namelist /model/ run_file, stdout

        ! Read in our namelists
        open(unit=funit, file=fstandalone, status='old', action='read')

        read(unit=funit, nml=model)
        read(unit=funit, nml=initvars)

        close(unit=funit)

        ! Initialize status of run 
        isValid = .true.
      
        ! Run the model
        call hoc_model
     .       ( C1, C1b, C1c, C2, C2b, C2c, 
     .         C2rt, C2thl, C2rtthl, C4, C5,
     .         C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc,
     .         C7, C7b, C7c, C8, C8b, C10, 
     .         C11, C11b, C11c, C12, C13,
     .         nu1, nu2, nu6, nu8, nu_r, beta,
     .         gamma_coef, gamma_coefb, gamma_coefc, 
     .         c_K, lmin_coef, taumin, taumax, mu,
     .         trim( run_file ), isValid, stdout )

      if ( .not. isValid ) then 
        stop "Model wasn't valid, check your constants"
      else
        stop "Program exited normally"
      endif 

      end program hoc_standalone
