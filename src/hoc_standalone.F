!----------------------------------------------------------------------
! $Id: hoc_standalone.F,v 1.16 2007-03-14 23:34:47 dschanen Exp $

        program hoc_standalone

!       Description:
!       This is essentially a minimalist frontend for HOC.
!
!----------------------------------------------------------------------

        use hoc, only: hoc_model

        use param_index, only: nparams
        use parameters, only: read_parameters

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
        logical :: isValid 

        ! Namelist
        namelist /model/ run_file, stdout

!-----------------------------------------------------------------------

        ! Read in model constant values
        call read_parameters( iunit, fstandalone, params )

        ! Read in model namelist
        open(unit=iunit, file=fstandalone, status='old', action='read')

        read(unit=iunit, nml=model)

        close(unit=iunit)

        ! Initialize status of run 
        isValid = .true.
      
        ! Run the model
        call hoc_model
     .       ( params, trim( run_file ), isValid, stdout )

      if ( .not. isValid ) then 
        stop "Model wasn't valid, check your constants"
      else
        stop "Program exited normally"
      end if 

      end program hoc_standalone
