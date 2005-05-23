!----------------------------------------------------------------------
!	Program HOC_STANDALONE
!
!       This is essentially a minimalist frontend for the hoc subroutine
!
!----------------------------------------------------------------------

      program hoc_standalone

      use hoc

      implicit none

!----- Program variables ----------------------------------------------
      logical isValid 

!----- Model constants ------------------------------------------------
! Notes:  C4 and C10 do not currently play a role in the model 3/18/05
!
      real       :: C1, C2, C4, C5, C6, C7, C8, C10, C11
      real       :: nu1, nu2, nu6, nu8
      real       :: beta
      real       :: gamma_coef, c_K
      real       :: lmin_coef, taumin, taumax
      real       :: mu

!----- Run information-------------------------------------------------
      character(len=40) :: runfile
      character(len=30) :: run_name

      logical           :: stdout

!----- Namelists ------------------------------------------------------
      namelist /initvars/ C1, C2, C4, C5, C6, C7, C8, C10, C11, 
     .                    nu1, nu2, nu6, nu8, beta, gamma_coef,
     .                    c_K, lmin_coef, taumin, taumax, mu
      namelist /model/ runfile, run_name, stdout

!----- Read in our namelists ------------------------------------------
       open(10, file='standalone.in', status='old')

       read(10, nml=model)
       read(10, nml=initvars)

       close(10)

!----- Initialize status of run ----------------------------------------
      isValid = .true.
      
!----- Run the model----------------------------------------------------
      call hoc_sub (C1, C2, C4, C5, C6, C7, C8, C10, C11, 
     .              nu1, nu2, nu6, nu8, beta,
     .              gamma_coef, c_K, lmin_coef, taumin, taumax, mu,
     .              runfile, run_name, isValid, stdout)

      if ( .not. isValid ) then 
        stop "Model wasn't valid, check your constants"
      else
        stop "Program exited normally"
      endif 

      end program hoc_standalone
