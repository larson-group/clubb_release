!-----------------------------------------------------------------------
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
      real       :: C1, C2, C4, C5, C6, C7, C8, C10, C11
      real       :: nu1, nu2, nu6, nu8
      real       :: beta
      real       :: gamma_coef, c_K
      real       :: lmin_coef, taumin, taumax
      real       :: mu

!----- Run information--------------------------------------------------
!      character(len=40) :: datafile
      character(len=40) :: runfile
      character(len=30) :: run_name

      logical           :: stdout
      real              :: sample_ratio

!----- Internal variables ----------------------------------------------
      logical isValid


!----- Namelists -------------------------------------------------------
      namelist /initvars/ C1, C2, C4, C5, C6, C7, C8, C10, C11, 
     .                    nu1, nu2, nu6, nu8, beta, gamma_coef,
     .                    c_K, lmin_coef, taumin, taumax, mu
      namelist /model/ runfile, stdout
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
      call hoc_model( C1, C2, C4, C5, C6, C7, C8, C10, C11, 
     .                nu1, nu2, nu6, nu8, beta,
     .                gamma_coef, c_K, lmin_coef, taumin, taumax, mu,
     .                runfile, isValid, stdout,
     .                grads_fields_reader, sample_ratio )

      if ( .not. isValid ) then 
        stop "Model wasn't valid, check your constants and field input"
      else
        stop "Program exited normally"
      endif 

      end program hoc_inputfields
!-----------------------------------------------------------------------
