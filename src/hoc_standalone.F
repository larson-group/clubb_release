!----------------------------------------------------------------------
! $Id: hoc_standalone.F,v 1.7 2005-08-06 13:37:53 griffinb Exp $
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
      real       :: C1, C2rt, C2thl, C2rtthl, C4, C5, C6rt, C6thl
      real       :: C7, C8, C10, C11
      real       :: nu1, nu2, nu6, nu8, nu_r
      real       :: beta
      real       :: gamma_coef, c_K
      real       :: lmin_coef, taumin, taumax
      real       :: mu

!----- Run information-------------------------------------------------
      character(len=50) :: run_file

      logical           :: stdout

!----- Internal variables ----------------------------------------------
      logical isValid 

!----- Namelists ------------------------------------------------------
      namelist /initvars/ C1, C2rt, C2thl, C2rtthl, C4, C5, C6rt, C6thl,
     .                    C7, C8, C10, C11, 
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
      call hoc_model( C1, C2rt, C2thl, C2rtthl, C4, C5, 
     .                C6rt, C6thl, C7, C8, C10, C11, 
     .                nu1, nu2, nu6, nu8, nu_r, beta,
     .                gamma_coef, c_K, lmin_coef, taumin, taumax, mu,
     .                run_file, isValid, stdout )

      if ( .not. isValid ) then 
        stop "Model wasn't valid, check your constants"
      else
        stop "Program exited normally"
      endif 

      end program hoc_standalone
