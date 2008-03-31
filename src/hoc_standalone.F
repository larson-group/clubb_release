!----------------------------------------------------------------------
! $Id: hoc_standalone.F,v 1.18 2008-03-31 18:10:55 faschinj Exp $

        program hoc_standalone

!       Description:
!       This is essentially a minimalist frontend for HOC.
!
!----------------------------------------------------------------------

        use hoc, only: hoc_model
        
        use error_code
        use param_index, only: nparams
        use parameters, only: read_parameters, params_list

        implicit none

        ! External
        intrinsic :: trim

        ! Constant parameters
        integer, parameter :: iunit = 10

        character(len=13), parameter :: fstandalone = "standalone.in"

        ! Run information
        real, dimension(nparams) ::
     .  params  ! Array of the model constants

        character(len=50) :: run_file ! Text file with the namelists

        logical :: stdout

        ! Internal variables
        integer :: code 

        integer :: i ! Loop iterator

        ! Namelist
        namelist /model/ run_file, stdout

!-----------------------------------------------------------------------

        ! Read in model constant values
        call read_parameters( iunit, fstandalone, params )

        ! Read in model namelist
        open(unit=iunit, file=fstandalone, status='old', action='read')

        read(unit=iunit, nml=model)

        close(unit=iunit)

        ! If standard output (stdout) is selected, print the list of
        ! parameters that are being used to the screen before the run.
        if ( stdout ) then
           write(unit=*,fmt='(4x,A9,5x,11x,A5)') "Parameter", "Value"
           write(unit=*,fmt='(4x,A9,5x,11x,A5)') "---------", "-----"
           do i = 1, nparams, 1
              write(unit=*,fmt='(A18,F27.20)')
     .           params_list(i) // " = ", params(i)
           enddo
        endif

        ! Initialize status of run 
        code = CLUBB_no_error
      
        ! Run the model
        call hoc_model
     .       ( params, trim( run_file ), code, stdout )

      if ( code == CLUBB_singular_matrix .or.
     .     code == CLUBB_var_equals_NaN .or.
     .     code == CLUBB_var_out_of_bounds .or.
     .     code == CLUBB_bad_lapack_arg .or.
     .     code == CLUBB_rtm_level_not_found ) then
      
        stop "Model wasn't valid, check your constants"
      else
        stop "Program exited normally"
      end if 

      end program hoc_standalone
