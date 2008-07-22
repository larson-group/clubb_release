!-----------------------------------------------------------------------
! $Id: hoc_inputfields.F90,v 1.1 2008-07-22 16:04:13 faschinj Exp $

        program hoc_inputfields

!       Description: 
!       This is essentially a minimalist frontend for the hoc subroutine
!       This version is modified to allow the input of LES, HOC or 
!       some other pre-calculated input data.
!-----------------------------------------------------------------------
        use hoc, only: hoc_model

        use inputfields, only: datafile, input_rtm, input_vm,  & ! Variable(s)
            input_um, input_thlm, input_rtpthlp, input_upwp,  & 
            input_vpwp, input_wp3, input_rtp2, input_thlp2, & 
            input_wp2, input_wprtp, input_wpthlp, input_type

        use inputfields, only: set_filenames ! Procedure(s)

        use param_index, only: nparams ! Variable(s)

        use parameters, only: read_parameters ! Procedure(s)

        use error_code, only: clubb_no_error ! Variable(s)
     
        use error_code, only: fatal_error ! Procedure(s)

        implicit none

        ! Constant parameters
        integer, parameter :: iunit = 10

        character(len=14), parameter :: finputfields = "inputfields.in"

        ! Run information
        real, dimension(nparams) :: & 
        params  ! Array of the model constants

        character(len=50) :: run_file ! Text file with the namelists

        logical :: stdout    ! Whether to print iteration number, etc.

        real :: sample_ratio ! Ratio of LES data to model timestep

        integer :: err_code   ! Numerical diagnostic


        ! Namelist definitions
        namelist /model/ run_file, stdout

        namelist /setfields/ datafile, input_type, & 
          input_um, input_vm, input_rtm, input_thlm, & 
          input_wp2, input_wprtp, input_wpthlp,  & 
          input_wp3, input_rtp2, input_thlp2,  & 
          input_rtpthlp, input_upwp, input_vpwp

!-----------------------------------------------------------------------

        ! Read in model constant values
        call read_parameters( iunit, finputfields, params )

        ! Read in our namelists
        open(unit=iunit, file=finputfields, status='old', action='read')

        read(unit=iunit, nml=model)
        read(unit=iunit, nml=setfields)

        close(unit=iunit)

        ! Initialize the status of the run
        err_code = clubb_no_error

        ! Setup the GrADS file reader
        call set_filenames( )

        ! Run the model
        call hoc_model & 
             ( params, trim( run_file ), err_code, stdout, .true. )

      if ( fatal_error(err_code)) then
       
        stop "Model wasn't valid, check your constants and field input"
      else
        stop "Program exited normally"
      endif 

      end program hoc_inputfields
!-----------------------------------------------------------------------
