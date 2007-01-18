!-----------------------------------------------------------------------
! $Id: hoc_inputfields.F,v 1.13 2007-01-18 23:54:32 dschanen Exp $

        program hoc_inputfields

!       Description: 
!       This is essentially a minimalist frontend for the hoc subroutine
!       This version is modified to allow the input of LES, HOC or 
!       some other pre-calculated input data.
!-----------------------------------------------------------------------
        use hoc, only: hoc_model
        use inputfields

        implicit none

        ! Model constants 
        ! Note:  C4 and C10 do not currently play a role in the model
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
        character(len=50) :: run_file ! Text file with the namelists

        logical :: stdout    ! Whether to print iteration number, etc.

        real :: sample_ratio ! Ratio of LES data to model timestep

        logical :: isValid   ! Numerical diagnostic


        ! Namelist definitions
        namelist /model/ run_file, stdout

        namelist /initvars/ C1, C1b, C1c, C2, C2b, C2c, 
     .    C2rt, C2thl, C2rtthl, C4, C5,
     .    C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc,
     .    C7, C7b, C7c, C8, C8b,
     .    C10, C11, C11b, C11c, C12,
     .    nu1, nu2, nu6, nu8, nu_r, beta, 
     .    gamma_coef, gamma_coefb, gamma_coefc,
     .    c_K, lmin_coef, taumin, taumax, mu

        namelist /setfields/ datafile, sample_ratio, input_type,
     .    input_um, input_vm, input_rtm, input_thlm,
     .    input_wp2, input_wprtp, input_wpthlp, 
     .    input_wp3, input_rtp2, input_thlp2, 
     .    input_rtpthlp, input_upwp, input_vpwp

!-----------------------------------------------------------------------

        ! Read in our namelists
        open(unit=10, file='inputfields.in', status='old')

        read(unit=10, nml=model)
        read(unit=10, nml=initvars)
        read(unit=10, nml=setfields)

        close(unit=10)

        ! Initialize the status of the run
        isValid = .true.

        ! Run the model
        call hoc_model( C1, C1b, C1c,
     .                  C2, C2b, C2c,
     .                  C2rt, C2thl, C2rtthl, C4, C5,
     .                  C6rt, C6rtb, C6rtc, 
     .                  C6thl, C6thlb, C6thlc,
     .                  C7, C7b, C7c, 
     .                  C8, C8b, C10, 
     .                  C11, C11b, C11c, C12,
     .                  nu1, nu2, nu6, nu8, nu_r, 
     .                  beta,
     .                  gamma_coef, gamma_coefb, gamma_coefc, 
     .                  c_K, lmin_coef, 
     .                  taumin, taumax, mu,
     .                  run_file, isValid, stdout,
     .                  grads_fields_reader, sample_ratio )

      if ( .not. isValid ) then 
        stop "Model wasn't valid, check your constants and field input"
      else
        stop "Program exited normally"
      endif 

      end program hoc_inputfields
!-----------------------------------------------------------------------
