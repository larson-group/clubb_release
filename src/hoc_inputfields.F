!-----------------------------------------------------------------------
! $Id: hoc_inputfields.F,v 1.11 2006-04-17 15:03:44 dschanen Exp $

!  Program HOC_INPUTFIELDS
 
!  This is essentially a minimalist frontend for the hoc subroutine
!  This version is modified to allow the input of LES, HOC or 
!  some other Pre-calculated data

! Calls:
! subroutine hoc_model

!-----------------------------------------------------------------------
      program hoc_inputfields
      use hoc
      use inputfields

      implicit none

!----- Model constants -------------------------------------------------
! Notes:  C4 and C10 do not currently play a role in the model 3/18/05
!
      real       :: C1, C1b, C1c, C2, C2b, C2c, C4, C5
      real       :: C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc
      real       :: C7, C7b, C7c, C8, C8b, C10, C11, C11b, C11c
      real       :: C12
      real       :: nu1, nu2, nu6, nu8, nu_r
      real       :: beta
      real       :: gamma_coef, gamma_coefb, gamma_coefc, c_K
      real       :: lmin_coef, taumin, taumax
      real       :: mu

!----- Run information--------------------------------------------------
      character(len=50) :: run_file

      logical           :: stdout
      real              :: sample_ratio

!----- Internal variables ----------------------------------------------
      logical isValid


!----- Namelists -------------------------------------------------------
      namelist /model/ run_file, stdout

      namelist /initvars/ C1, C1b, C1c, C2, C2b, C2c, C4, C5,
     .                    C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc,
     .                    C7, C7b, C7c, C8, C8b,
     .                    C10, C11, C11b, C11c, C12,
     .                    nu1, nu2, nu6, nu8, nu_r, beta, 
     .                    gamma_coef, gamma_coefb, gamma_coefc,
     .                    c_K, lmin_coef, taumin, taumax, mu

      namelist /setfields/ datafile, sample_ratio, input_type,
     .                     input_um, input_vm, input_rtm, input_thlm,
     .                     input_wp2, input_wprtp, input_wpthlp, 
     .                     input_wp3, input_rtp2, input_thlp2, 
     .                     input_rtpthlp, input_upwp, input_vpwp

!----- Read in our namelists -------------------------------------------
      open(10, file='inputfields.in', status='old')

      read(10, nml=model)
      read(10, nml=initvars)
      read(10, nml=setfields)

      close(10)

!----- Initialize status of run ----------------------------------------
      isValid = .true.

!----- Run the model----------------------------------------------------
      call hoc_model( C1, C1b, C1c, C2, C2b, C2c, C4, C5,
     .                C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc,
     .                C7, C7b, C7c, C8, C8b, C10, 
     .                C11, C11b, C11c, C12,
     .                nu1, nu2, nu6, nu8, nu_r, beta,
     .                gamma_coef, gamma_coefb, gamma_coefc, 
     .                c_K, lmin_coef, taumin, taumax, mu,
     .                run_file, isValid, stdout,
     .                grads_fields_reader, sample_ratio )

      if ( .not. isValid ) then 
        stop "Model wasn't valid, check your constants and field input"
      else
        stop "Program exited normally"
      endif 

      end program hoc_inputfields
!-----------------------------------------------------------------------
