!----------------------------------------------------------------------
! $Id: hoc_standalone.F,v 1.10 2005-10-29 21:25:07 griffinb Exp $
!  Program HOC_STANDALONE
!
!  This is essentially a minimalist frontend for the hoc subroutine
!
!----------------------------------------------------------------------
      program hoc_standalone
      use hoc

      implicit none

!----- Model constants ------------------------------------------------
! Notes:  C4 and C10 do not currently play a role in the model 3/18/05
!
      real C1, C1b, C1c, C2rt, C2thl, C2rtthl
      real C4, C5, C6rt, C6rtb, C6rtc, C6thl
      real C7, C7b, C7c, C8, C8b, C10, C11, C11b, C11c, C12
      real nu1, nu2, nu6, nu8, nu_r
      real beta, gamma_coef, c_K
      real lmin_coef, taumin, taumax
      real mu

!----- Run information-------------------------------------------------
      character(len=50) run_file

      logical stdout

!----- Internal variables ----------------------------------------------
      logical isValid 

!----- Namelists ------------------------------------------------------
      namelist /initvars/ C1, C1b, C1c, C2rt, C2thl, C2rtthl,
     .                    C4, C5, C6rt, C6rtb, C6rtc, C6thl,
     .                    C7, C7b, C7c, C8, C8b, C10, C11, C11b, C11c,
     .                    C12,
     .                    nu1, nu2, nu6, nu8, nu_r, beta, gamma_coef,
     .                    c_K, lmin_coef, taumin, taumax, mu
      namelist /model/ run_file, stdout

!----- Read in our namelists ------------------------------------------
       open(10, file='standalone.in', status='old')

       read(10, nml=model)
       read(10, nml=initvars)

       close(10)

!----- Initialize status of run ----------------------------------------
      isValid = .true.
      
!----- Run the model----------------------------------------------------
      call hoc_model( C1, C1b, C1c, C2rt, C2thl, C2rtthl,
     .                C4, C5, C6rt, C6rtb, C6rtc, C6thl,
     .                C7, C7b, C7c, C8, C8b, C10, 
     .                C11, C11b, C11c, C12,
     .                nu1, nu2, nu6, nu8, nu_r, beta,
     .                gamma_coef, c_K, lmin_coef, taumin, taumax, mu,
     .                run_file, isValid, stdout )

      if ( .not. isValid ) then 
        stop "Model wasn't valid, check your constants"
      else
        stop "Program exited normally"
      endif 

      end program hoc_standalone
