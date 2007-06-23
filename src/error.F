!-----------------------------------------------------------------------
! $Id: error.F,v 1.47 2007-06-23 20:20:55 griffinb Exp $

        module error
 
!       Description:

!       subroutine tuner_init: reads in namelists /stats/, /cases/, 
!          /initvars/, & /variance/ from 'error.in' 
!       It then uses them to setup the initial x_array of independent 
!       variables, i.e. the hoc constants, and allocate the runtime arrays
!       for each of the model runs and each of the variables

!       function min_les_hoc_diff:  A driver for the hoc program/module.
!       Calls hoc_model, reads in les & hoc results from GRADS files, and 
!       calculates the average difference between the two over all z-levels.

!       subroutine output_results_stdout : 
!       Prints the results of tuning to the terminal, or
!       to a logfile with redirection(done outside the program in scripting)

!       subroutine output_nml_tuner :  
!       Generates the error.in file using the current constants.

!       subroutine output_nml_standalone :  
!       Generates the standalone.in file using the current constants.
!       Standalone hoc is only configured to run a single model, and so only
!       the first model is used to make the namelist.

!-----------------------------------------------------------------------
        use param_index
        use parameters, only: read_parameters, read_param_spread
      
        implicit none

        integer, public :: ndim ! Size of the simplex

        ! isValid is an important logical used by the min_diff to 
        ! determine whether HOC has become numerically unstable
        logical, public :: 
     .  isValid 
 
        ! inv_count is a modular counter [1-3] used to determine 
        ! which file to output to if stdout_on_invalid is true.
        integer, private :: 
     .  inv_count

!----------------------------------------------------------------------- 
      
        integer, parameter, private :: 
     .  max_run = 12, ! Maximum model runs the tuner can handle at a time
     .  max_variables = 32 ! This number / 2 is maximum variables to tune for

        real :: 
     .  ftol,       ! The precision to tune for
     .  anneal_temp ! Initial temperature for the simulated annealing algorithm

        integer ::
     .  anneal_iter,! Number of annealing iterations to perform
     .  tune_type,  ! Toggle for downhill simplex of simulated annealing
     .  c_total,    ! Total number of simulation cases to tune over
     .  v_total     ! Total number of variables to tune over

        logical ::
     .  results_stdout,
     .  results_file,
     .  stdout_on_invalid

        character(len=10), dimension(:), allocatable :: 
     .  hoc_v, ! Variables in HOC GrADS files
     .  les_v  ! Variables in LES GrADS files

        integer, dimension(:,:), allocatable :: 
     .  time ! Time intervals

        ! Additions for using imposed weights as scaling factors
        logical :: linitialize_sigma

        real, dimension(:,:), allocatable :: 
     .  err_terms,
     .  invsigma2,
     .  min_err_terms,
     .  init_err_terms

        real, dimension(:), allocatable :: weight_case, weight_var
        
        ! End additions for using imposed weights

        integer, dimension(:), allocatable ::
     .  z_i, ! Initial z level for tuning purposes
     .  z_f  ! Final z level for tuning purposes

        character(len=50), dimension(:), allocatable :: 
     .  run_file,       ! Model run files
     .  hoc_stats_file, ! Model GrADS files
     .  les_stats_file  ! Model GrADS files

        ! Various Variables for returning results ----------------------
        integer :: iter ! Total number of iterations amoeba spent 
                        ! calculating optimal values

        real ::
     .  init_err, ! error for the initial constants
     .  min_err   ! the lowest the minimization algorithm could go

        real, dimension(nparams) ::
     .  params  ! Vector of all possible HOC parameters

        integer, dimension(nparams) ::
     .  x_index  ! Index of the params elements that are used in the simplex

        real, allocatable, dimension(:,:) :: 
     .  x_array ! Holds 2D simplex the HOC constant parameters

        real, allocatable, dimension(:) ::
     .  x_spread, ! Amount to vary each respec. constant by
     .  err_array ! cache of differences between the LES and HOC

        real, allocatable, dimension(:) :: 
     .  rand_vect ! A vector of random reals for initializing the x array

        contains

!-----------------------------------------------------------------------
        subroutine tuner_init( read_files )

!       Description:
!       Initializes x_array with constants from error.in
!       Allocates arrays for cases and tuning variables 
!       Initializes grads file names to read in 

!       References:
!       None
!-----------------------------------------------------------------------
        use constants, only: fstdout, fstderr

        implicit none

        ! Constant Variables

        integer, parameter :: 
     .  max_times = 10 ! max number of timesteps to compare

        character(len=8), parameter ::
     .  filename = "error.in"

        ! Input  Variables

        ! Determines whether to read in the namelists and do the 
        ! initial allocation of arrays
        logical, intent(in) :: 
     .  read_files 

        ! Local Variables

        real, dimension(nparams) :: rtmp ! Scratch space

        integer :: i, j ! looping variables
      
!----------------------------------------------------------------------- 

        ! Namelist vars for determining which variable to tune for:

        !  time_nl:   Order pairs of time intervals to analyze
        !  z_i_nl:    initial z-level to begin reading in for tuning
        !  z_f_nl:    final z-level to end reading in for tuning

        integer, dimension(max_run, max_times):: time_nl 
        integer, dimension(max_run)           :: z_i_nl, z_f_nl

        ! Addition to use imposed weights as scaling factors
        real, dimension(max_run)       :: weight_case_nl
        real, dimension(max_variables) :: weight_var_nl

        character(len=50), dimension(max_run) :: 
     .  run_file_nl,
     .  hoc_stats_file_nl,
     .  les_stats_file_nl

        character(len=10), dimension(max_variables) :: t_variables

        ! Namelists read from error.in 
        namelist /stats/ 
     .    ftol, tune_type, anneal_temp, anneal_iter,
     .    results_stdout, results_file, stdout_on_invalid, t_variables,
     .    weight_var_nl

        namelist /cases/ 
     .    les_stats_file_nl, hoc_stats_file_nl,
     .    run_file_nl, z_i_nl, z_f_nl, time_nl, weight_case_nl

        ! Reset iteration counter (set by amoeba)
        iter = 0

        ! Reset invalid run counter
        inv_count = 0

        ! Re-read namelists if requested
        if ( read_files ) then

          ! Initialize all compile time arrays to zero
          time_nl = 0 

          ! Imposed weights as scaling factors
          weight_case_nl = 0.0
          weight_var_nl  = 0.0

          z_i_nl = 0 
          z_f_nl = 0

          ! Initialize variable names to spaces
          t_variables(1:max_variables)  = "          "

          ! Open our namelist input file
          open(unit=10, file=filename, status='old')

          ! Determine which files to read data from based on namelist
          read(unit=10, nml=stats)

          ! Read in the models to be run
          read(unit=10, nml=cases)

          ! Close our input namelist file
          close(unit=10)

          ! Read in initial constant values
          call read_parameters( 10, filename, params )

          ! Allocate the arrays for the tuning variables

          do i = 1, max_variables, 2 ! 1, 3, 5, 7
            if (t_variables(i) == "          ") exit
            v_total = (i + 1) / 2
          end do 

          allocate( hoc_v(v_total), les_v(v_total) )

          ! Allocate the arrays for the run cases
          do i=1, max_run  
            if (z_f_nl(i) == 0 ) exit
            c_total = i
          end do

          allocate(
     .      z_i(c_total), z_f(c_total), time(c_total, max_times), 
     .      run_file(c_total), 
     .      les_stats_file(c_total), hoc_stats_file(c_total) )

          allocate(
     .      err_terms(c_total, v_total),
     .      min_err_terms(c_total, v_total),
     .      init_err_terms(c_total, v_total),
     .      invsigma2(c_total, v_total),
     .      weight_case(c_total), weight_var(v_total) )

          ! Transfer the variable numbers to hoc_v and les_v
          do i=1, v_total
            hoc_v(i) = t_variables(i*2 - 1)
            les_v(i) = t_variables(i*2)
          end do 

          ! Transfer the case information to run-time arrays
          do i = 1, c_total 
            z_i(i)              = z_i_nl(i)
            z_f(i)              = z_f_nl(i)
            les_stats_file(i)   = les_stats_file_nl(i)
            hoc_stats_file(i)   = hoc_stats_file_nl(i)
            run_file(i)         = run_file_nl(i)
            time(i,1:max_times) = time_nl(i,1:max_times)
          end do

         ! Use imposed weights as scaling factors
         weight_case(1:c_total) = weight_case_nl(1:c_total)
         weight_var(1:v_total)  = weight_var_nl(1:v_total)

         ! Setup the simplex

         call read_param_spread( 10, filename, x_index, rtmp, ndim )

         if ( ndim == 0 ) then
           write(fstderr,*) "You must vary at least one parameter"
           stop
         end if

         allocate( rand_vect(ndim), x_array(ndim+1,ndim),
     .             x_spread(ndim), err_array(ndim+1) )

          ! Initialize the HOC parameter spread 
          x_spread(1:ndim)  = rtmp(x_index(1:ndim))

          ! Copy varying parameters into the first row of the simplex
          x_array(1,1:ndim) = params(x_index(1:ndim))
      
          ! Attempt to generate a pseudo-random seed using a file 
          ! generated from /dev/random.  File is an ASCII text file 
          ! and can be edited manually.
          call read_random_seed( "rand_seed.dat" )

        end if  ! read_files
!-----------------------------------------------------------------------

        ! Fill in the remaining values of the array by varying the initial
        ! vector (i.e. the first column of the array) by a small multiple
        do j = 1, ndim

          call random_number( rand_vect(1:ndim) )

          do i = 2, ndim+1, 1
          ! Vince Larson made entries of x_array random  10 Feb 2005
!            x_array(i,j) = x_array(1,j)*
!     .                     (1.0+((real(i)-1.)/real(ndim)*0.5))
            x_array(i,j) = x_array(1,j)*
     .      ( (1.0 - x_spread(j)) + rand_vect(i-1)*x_spread(j)*2 )
          ! End of Vince Larson's change
          end do ! i..ndim+1

        end do ! j..ndim


        ! First call is used to initialize weights

        linitialize_sigma = .true.
        err_array(1) =  min_les_hoc_diff( x_array(1,1:ndim) )
        linitialize_sigma = .false.

        ! Note: min_les_hoc_diff is written to deal with undefined and 
        ! invalid values for variations on the initial vector, but that 
        ! algorithm relies on the initial vector being valid. 

        if ( .not. isValid ) then
          write(fstderr,*) "Initial variable values must be valid."
          stop
        end if

        ! Save initial error

        init_err = err_array(1) 
        init_err_terms = err_terms

        ! Other initialization runs

        ! Initialize the 'y' vector for amoeba 
        ! This is done by calling min_les_hoc_diff with the initial vector 
        do i = 2, ndim+1, 1
          err_array(i) = min_les_hoc_diff( x_array(i,1:ndim) )
        end do 

        write(unit=fstdout,fmt=*) "err_array:"
        write(unit=fstdout,fmt='(6e12.5)') err_array

        return
        end subroutine tuner_init

!-----------------------------------------------------------------------
        real function min_les_hoc_diff( x_vector )

!       Description: 
!       Function that returns the sum of the error between the dependent
!       variable(i.e. the variable we want to match) in each of the models

!       References:
!       _Numerical Recipes in Fortran 77_ P.402-406 (Description)
!       _Numerical Recipes in Fortran 90_ source code (Routine)
!-----------------------------------------------------------------------
 
        use hoc, only: hoc_model
        use grads_common
        use parameters, only: params_list

        implicit none

        ! External
        intrinsic :: achar, modulo, minval, maxval, trim, any, real

        ! Input Variables

        ! The interface declaration in the nr module prevents us from making
        ! this a fixed length array declaration
        real, dimension(:), intent(in) :: 
     .  x_vector ! Tuning vector(ndim dimension) contains
                 ! parameterization constants (C1, C2 etc.)

        ! Local Variables

        real, dimension(nparams) ::
     .  params_local ! Local copy of the HOC parameters fed into hoc_model

        ! These are read after each run from the GrADS control files
        integer :: 
     .  les_nz, ! Extent of the LES domain in the z dimension
     .  hoc_nz  ! Extent of the HOC domain in the z dimension
      
        character(50) :: 
     .  errorfile ! nml filename for invalid runs

        integer :: 
     .  AllocateStatus ! For hoc_zl, les_zl

        real :: 
     .  err_sum  ! scalar sum of all z-levels

        integer :: 
     .  c_terms ! num of terms in err_sum (for normalization)

        ! LES and HOC values over nz z-levels
        real, dimension(:), allocatable :: 
     .  hoc_zl,
     .  hoc2_zl,
     .  les_zl 

        real, dimension(:,:), allocatable :: 
     .  err_sums  ! To save breakdown of cost function

        real :: 
     .  les_minmax ! The largest LES value subtracted from the 
                   ! smallest value of all zlvl's (for normalization)

        logical :: 
     .  error ! Used to det. if reading the variable failed

        logical, dimension(c_total) :: 
     .  run_stat ! isValid over each model case

        integer :: i, j, c_run ! looping variables

!-----------------------------------------------------------------------

        ! Output information every 10 iterations if stdout is enabled;
        ! Amoeba's unusual calling convention makes this happen less 
        ! often than might be expected.
        if ( results_stdout
     .      .and. ( modulo(iter, 10) == 0 )
     .      .and. iter /= 0  ) then

          write(unit=*,fmt='(A12,I10)') "Iteration: ", iter
          write(unit=*,fmt='(A12)') "Parameters: "

          do i = 1, ndim, 1
            j = x_index(i)
            write(unit=*,fmt='(A18,F27.20)') params_list(j)//" = ",
     .        x_vector(i)
          end do

        end if

        allocate( err_sums(c_total, v_total) )        

        ! Initialize
        err_sum  = 0.0
        err_sums = 0.0
        c_terms  = 0
        isValid  = .true.
        run_stat(1:c_total) = .true.

        ! Copy simplex into a vector of all possible HOC parameters
        do i=1, nparams, 1
          ! If the variable isn't in the simplex, leave it as is
          params_local(i) = params(i)
          do j=1, ndim, 1
            if ( i == x_index(j) ) then
              ! Copy variable from x_vector argument
              params_local(i) = x_vector(j)
              exit
            else
              ! Continue searching the list
              cycle
            end if
          end do
        end do

!-----------------------------------------------------------------------

        ! Cycle through all the model cases specified for HOC

        ! OpenMP directives should work as expected now, assuming new
        ! model variables are declared threadprivate -dschanen 31 Jan 2007
!$omp   parallel do default(none), private(c_run),
!$omp.    shared(params_local, run_file, run_stat, c_total)
        do c_run=1, c_total, 1

!         write(unit=*,fmt=*) "Calling HOC with case "//
!    .      trim( run_file(c_run) )

          ! Run the HOC model with parameters as input

          call hoc_model
     .         ( params_local, run_file(c_run), run_stat(c_run), 
     .           .false. )

        end do ! 1..c_run

!$omp   end parallel do

        do c_run = 1, c_total, 1
          if ( .not. run_stat(c_run) ) isValid = .false.
        end do

!-----------------------------------------------------------------------

        ! Now check if hoc has blown up, i.e. if hoc has set a variable to NaN, 
        ! or encountered a failure in the matrix solver routines

        ! If it has, it returns higher value than those previous to
        ! Amoeba (the downhill simplex)
        if ( .not. isValid ) then
          min_les_hoc_diff = 2 * maxval( err_array ) 
     .                     - minval( err_array )

          if ( stdout_on_invalid ) then
            inv_count = modulo( inv_count, 3 ) + 1 ! 1,2,3,1,2,3...
            errorfile = "error_crash_"// achar( inv_count+48 )// ".in"
            call output_nml_tuner( errorfile, x_vector(1:ndim) )
          end if
          return
        end if

!-----------------------------------------------------------------------

        do c_run=1, c_total, 1

          ! Determine how large the GrADS input is
          hoc_nz = grads_zlvl( hoc_stats_file(c_run) )
          les_nz = grads_zlvl( les_stats_file(c_run) )

          ! Allocate the arrays for reading in the GrADS plot data
          allocate( hoc_zl(hoc_nz), hoc2_zl(hoc_nz), 
     .              les_zl(les_nz), stat=AllocateStatus )
          if ( AllocateStatus /= 0 ) then 
            stop "Allocation of arrays in minimization function failed"
          end if

          ! Start with first HOC & LES variables, then loop through and 
          ! calculate the mean squared difference for all the variables
          do i=1, v_total, 1 

            ! Read in LES grads data for one variable, averaged
            ! over specified time intervals
            les_zl = 
     .      grads_average_interval
     .      ( les_stats_file(c_run), les_nz, 
     .        time(c_run,:), les_v(i), 1, error )

            if ( error ) stop "The specified LES variable was invalid"

            ! Read in HOC grads data for one variable, averaged
            ! over specified time intervals
            hoc_zl = 
     .      grads_average_interval
     .      ( hoc_stats_file(c_run), hoc_nz, 
     .        time(c_run,:), hoc_v(i), 1, error )

            if ( error ) stop "The specified HOC variable was invalid"

            ! The same variable, with npower = 2
            hoc2_zl = 
     .      grads_average_interval
     .      ( hoc_stats_file(c_run), hoc_nz,
     .        time(c_run,:), hoc_v(i), 2, error )

            if ( error ) stop "The specified HOC variable was invalid"

!-----------------------------------------------------------------------

            ! Calculate the mean squared difference between the HOC 
            ! and the LES variables

            ! In order to deal with differences in order of magnitude 
            ! between the variables, the err_sum equation has been 
            ! modified to normalize the values with respect to the 
            ! the minimum and maximum in the LES. -Dave Schanen

            les_minmax = maxval( les_zl ) - minval( les_zl )

            if ( les_minmax == 0.0 ) then
              stop "An LES var = 0 over all z-levels"
            end if

            ! Old code
!           err_sum = err_sum 
!    .              + mean_sqr_diff_zt( hoc_nz, les_nz, hoc_zl, 
!    .                                  les_zl, les_minmax )

            ! Chris Golaz modification: mean_sqr_diff_2_zt was designed to try
            ! to limit time noise in tuning simulations.
            ! New code
!           err_sum = err_sum 
!    .              + mean_sqr_diff_2_zt
!    .                ( hoc_nz, les_nz, hoc_zl, hoc2_zl, 
!    .                  les_zl, les_minmax )
            ! End new code
            ! Modification for weighting
            err_sums(c_run,i)
     .      = mean_sqr_diff_2_zt( hoc_nz, les_nz, hoc_zl, 
     .                            hoc2_zl, les_zl, les_minmax )

            c_terms = c_terms + 1

          end do ! i=1..v_total

          deallocate( hoc_zl, hoc2_zl, les_zl )

        end do     ! end of do c_run=1, c_total
!----------------------------------------------------------------------

        ! Return error averaged over all cases, variables, 
        ! and vertical levels
        ! Old Code
!       min_les_hoc_diff = err_sum / real( c_terms )

        !---------------------------------------------------------------
        ! Compute normalization factors to satisfy imposed weights
        !---------------------------------------------------------------
        if ( linitialize_sigma ) then
          do c_run=1,c_total
            do i=1,v_total
              invsigma2(c_run,i) 
     .        = weight_case(c_run)*weight_var(i) / err_sums(c_run,i)
            end do
          end do
        end if

        !---------------------------------------------------------------
        ! Compute normalized error
        !---------------------------------------------------------------
        err_sums = invsigma2 * err_sums
        err_sum  = sum( err_sums )

        !---------------------------------------------------------------
        ! Save total error and error contributions breakdown
        !---------------------------------------------------------------
        err_terms = err_sums
        min_les_hoc_diff = err_sum

        deallocate( err_sums )

        write(*,'(a,f12.5)') "Cost function= ", min_les_hoc_diff

        return
        end function min_les_hoc_diff

!----------------------------------------------------------------------
        subroutine output_results_stdout( )

!       Description:
!       Outputs the results of a tuning run to the terminal

!       References:
!       None
!----------------------------------------------------------------------
        use parameters, only: params_list
        use constants, only: fstdout
      
        implicit none

        ! Local variables

        integer :: i ! Loop iterator

        if ( tune_type == 0 ) print *, "Number of iterations:",  iter
      
        write(unit=*,fmt='(4x,A9,5x,10x,A7,10x,10x,A7)')
     .      "Parameter", "Initial", "Optimal"

        do i = 1, ndim, 1
          write(unit=*,fmt='(A18,2F27.20)') 
     .      params_list(x_index(i))//" = ", 
     .      params(x_index(i)), x_array(1,i)
        end do

        write(fstdout,fmt='(A20)') "Initial cost: "
        write(fstdout,fmt='(F15.6)') init_err
        write(fstdout,fmt='(A20)') "Optimal cost: "
        ! The $$ is here to make it easy to find with grep
        write(fstdout,fmt='(A3,F15.6)') "$$ ", min_err  
     
        write(unit=*,fmt=*) "Approx. percent increase in accuracy:", 
     .    ((init_err - min_err) / init_err*100.0), "%"

        return
        end subroutine output_results_stdout
!-----------------------------------------------------------------------
        subroutine output_nml_tuner( results_f, x_vector )

!       Description:
!       Output namelists to a formatted text file

!       References:
!       None

!       Notes:
!       You can do the same thing with
!       write(unit=<UNIT>,nml=<NAMELIST>), but at the time I wrote this
!       I was more ambitious and didn't like that fact that it came out
!       in all caps on pgf90. -dschanen 28 July 2006
!-----------------------------------------------------------------------
        use parameters, only: params_list

        implicit none

        ! External
        intrinsic :: achar, trim

        ! Parameter Constants
        character(len=1), parameter :: dbqt = '"'

        ! Input Variables
        character(len=*), intent(in) :: 
     .  results_f ! Name of the results file to write to

        real, intent(in), dimension(:) :: 
     .  x_vector ! The current tuning parameters

        ! Local Variables
        real, dimension(nparams) :: params_local

        integer :: i, j    ! loop variable

        character :: i_c   ! loop variable in ASCII

        ! Open a new file
        open(unit=20, file=results_f, 
     .       action="write", access="sequential")
      
        ! Write variables to results file 
        ! All this is based on the previous namelists in error.in, 
        ! except for the constants parameters for HOC

        write(unit=20,fmt=*) "! Parameter file " // results_f
        write(unit=20,fmt=*) "&stats"
        write(unit=20,fmt=*) "ftol = ", ftol
        if (results_stdout) then
          write(unit=20,fmt=*) "results_stdout = " // ".true."
        else
          write(unit=20,fmt=*) "results_stdout = " // ".false."
        end if
        if (results_file) then
          write(unit=20,fmt=*) "results_file = " // ".true."
        else
          write(unit=20,fmt=*) "results_file = " // ".false."
        end if
        if (stdout_on_invalid) then
          write(unit=20,fmt=*) "stdout_on_invalid = " // ".true."
        else
          write(unit=20,fmt=*) "stdout_on_invalid = " // ".false."
        end if
        write(unit=20,fmt=*) "t_variables = "
        do i=1, v_total
          write(unit=20,fmt=*) dbqt, hoc_v(i), dbqt, ",", 
     .      dbqt, les_v(i), dbqt, "," 
        end do
        write(unit=20,fmt=*) "weight_var_nl = ", weight_var
        write(unit=20,fmt=*) "anneal_temp = " , anneal_temp
        write(unit=20,fmt=*) "anneal_iter = " , anneal_iter
        write(unit=20,fmt=*) "tune_type   = " , tune_type
        write(unit=20,fmt=*) "/"

        write(unit=20,fmt=*) "&cases"
        do i = 1, c_total
          i_c = achar( i + 48 )
          write(unit=20,fmt=*) "les_stats_file_nl("// i_c //") = ", 
     .      dbqt, trim( les_stats_file(i) ), dbqt
          write(unit=20,fmt=*) "hoc_stats_file_nl("// i_c //") = ", 
     .      dbqt, trim( hoc_stats_file(i) ), dbqt
          write(unit=20,fmt=*) "run_file_nl("// i_c //") = ", 
     .      dbqt, trim( run_file(i) ), dbqt
          write(unit=20,fmt=*) "z_i_nl("// i_c //")  = " , z_i(i)
          write(unit=20,fmt=*) "z_f_nl("// i_c //")  = " ,  z_f(i)
          write(unit=20,fmt=*) "weight_case_nl("// i_c //") = " , 
     .      weight_case(i)
          write(unit=20,fmt=*) "time_nl("// i_c //",:)  = " , time(i,:)
        end do
        write(unit=20,fmt=*) "/"
      
        write(unit=20,fmt=*) "&initvars"

        ! Copy simplex into a vector of all possible HOC parameters
        do i=1, nparams, 1
          ! If the variable isn't in the simplex, leave it as is
          params_local(i) = params(i)
          do j=1, ndim, 1
            if ( i == x_index(j) ) then
              ! Copy variable from x_vector argument
              params_local(i) = x_vector(j)
              exit
            else
              cycle
            end if
          end do
        end do

        ! Output optimal values and all possible HOC parameters
        do i=1, nparams, 1
          write(unit=20,fmt='(A18,F27.20)')
     .      trim( params_list(i) )//" = ", params_local(i)
        end do
        write(unit=20,fmt=*) "/"

        write(unit=20,fmt=*) "&initspread"
        ! Copy the spread into a vector of all possible HOC parameters
        do i=1, nparams, 1
          ! If the variable isn't being changed, set it to zero
          params_local(i) = 0.0
          do j=1, ndim, 1
            if ( i == x_index(j) ) then
              ! Copy variable from x_vector argument
              params_local(i) = x_spread(j)
              exit
            else
              cycle
            end if
          end do
        end do

        ! Output the amount each variable was changed for the simplex
        do i=1, nparams, 1
          write(unit=20,fmt='(a18,f12.5)')
     .      trim( params_list(i) )//" = ", params_local(i)
        end do
        write(unit=20,fmt=*) "/"

        ! Close new namelist file
        close(unit=20)

        return
        end subroutine output_nml_tuner

!-----------------------------------------------------------------------
        subroutine output_nml_standalone ( results_f, x_vector )

!       Description:
!       Output namelists to a formatted file 

!       References:
!       None

!       Notes:
!       See note for output_nml_tuner, above. dschanen 28 July 2006
!-----------------------------------------------------------------------

        use parameters, only: params_list

        implicit none

        ! External
        intrinsic :: achar, trim

        ! Constant Parameter
        character(len=1), parameter :: dbqt = '"'

        ! Input variables
        character(len=*), intent(in) ::
     .  results_f ! Results file to write to

        real, intent(in), dimension(ndim) ::
     .  x_vector ! the current constants

        ! Local variables
        real, dimension(nparams) :: params_local

        integer   :: i, j  ! loop variables

        ! Open new file
        open(unit=20, file=results_f, 
     .       action="write", access="sequential")
      
        ! Write variables to namelist for standalone hoc. 
        ! All this is based on the previous error.in, except the constants

        write(unit=20,fmt=*) "! Parameter file " // results_f

        write(unit=20,fmt=*) "&model"
        write(unit=20,fmt=*) "run_file = "//dbqt//"hoc.in"//dbqt
        write(unit=20,fmt=*) "stdout   = ", ".true."
        write(unit=20,fmt=*) "/"

        write(unit=20,fmt=*) "&initvars"
        ! Copy simplex into a vector of all possible HOC parameters
        do i=1, nparams, 1
          ! If the variable isn't in the simplex, leave it as is
          params_local(i) = params(i)
          do j=1, ndim, 1
            if ( i == x_index(j) ) then
              ! Copy variable from x_vector argument
              params_local(i) = x_vector(j)
              exit
            else
              cycle
            end if
          end do
        end do

        ! Output optimal values and all possible HOC parameters
        do i=1, nparams, 1
          write(unit=20,fmt='(A18,F27.20)')
     .      trim( params_list(i) )//" = ", params_local(i)
        end do
        write(unit=20,fmt=*) "/"

        ! Close new file 
        close(unit=20)

        return
        end subroutine output_nml_standalone

!-----------------------------------------------------------------------
        real function mean_sqr_diff_zm
     .                ( hoc_nz, les_nz, hoc_zl, les_zl, norm_term )
!       Description
!       Calculate the mean squared difference between two input vectors, 
!       then normalize.

!       References:
!       None

!       Notes:
!       Configured to do interpolation on LES / HOC comparisons on the 
!       momentum grid.
!-----------------------------------------------------------------------
        implicit none

        ! External
        intrinsic sum

        ! Input Variables
        integer, intent(in) :: 
     .  hoc_nz, ! Vertical extent for HOC
     .  les_nz  ! Vertical extent for the LES

        real, intent(in), dimension(hoc_nz) :: 
     .  hoc_zl  ! HOC GrADS variable [units vary]

        real, intent(in), dimension(les_nz) :: 
     .  les_zl  ! The LES GrADS variable [units vary]

        real, intent(in) :: 
     .  norm_term ! normalization term;  
                  !typically maxval(les) - minval(les) 

        ! Local Variables
        real, dimension(hoc_nz) :: 
     .  tmp_zl

!----------------------------------------------------------------------

        select case ( hoc_nz - les_nz )
        case ( 0 ) ! most cases
        ! Due to hoc's lower starting point, we can only use
        ! (total number of z-levels) - 1 (a maximum of 74 for BOMEX).
        ! The code below assumes the LES data are on an evenly spaced grid.
        ! (Need to interpolate hoc to LES' levels.  Right now we just
        !   compare adjacent z levels.  Vince Larson 12 Jan 2005)

          tmp_zl(1:hoc_nz-1) = 
     .    ( hoc_zl(1:hoc_nz) - les_zl(1:(les_nz)) ) / norm_term 

          tmp_zl = tmp_zl**2

          mean_sqr_diff_zm = sum( tmp_zl(1:(hoc_nz-1)), 1 ) 

        case ( 2 )  !  the DYCOMS II RF01 case
          tmp_zl(1:les_nz) = 
     .    ( hoc_zl(1:les_nz) - les_zl(1:les_nz) ) / norm_term

          tmp_zl = tmp_zl**2

          mean_sqr_diff_zm = sum( tmp_zl(1:(les_nz)), 1 )

        case default !
          stop "Not able to handle specified number of HOC z-levels"
        end select

        return
        end function mean_sqr_diff_zm

!-----------------------------------------------------------------------
        real function mean_sqr_diff_zt
     .                ( hoc_nz, les_nz, hoc_zl, les_zl, norm_term )
!       Description
!       Calculate the mean squared difference between two input vectors, 
!       then normalize.

!       References:
!       None

!       Notes:
!       Configured to do interpolation on LES / HOC comparisons on the 
!       thermodynamic grid.  We use a modified version for tuning runs
!       now, see below.
!-----------------------------------------------------------------------
        implicit none

        ! External
        intrinsic sum

        ! Input Variables
        integer, intent(in) :: 
     .  hoc_nz, ! Vertical extent for HOC
     .  les_nz  ! Vertical extent for the LES

        real, intent(in), dimension(hoc_nz) :: 
     .  hoc_zl  ! HOC GrADS variable [units vary]

        real, intent(in), dimension(les_nz) :: 
     .  les_zl  ! The LES GrADS variable [units vary]

        real, intent(in) :: 
     .  norm_term ! normalization term;  
                  !typically maxval(les) - minval(les) 

        ! Local Variables
        real, dimension(hoc_nz) :: 
     .  tmp_zl

!----------------------------------------------------------------------

        select case ( hoc_nz - les_nz )
        case ( 0 ) ! most cases
        ! Due to hoc's lower starting point, we can only use
        ! (total number of z-levels) - 1 (a maximum of 74 for BOMEX).
        ! The code below assumes the LES data are on an evenly spaced grid.
        ! (Need to interpolate hoc to LES' levels.  Right now we just
        !   compare adjacent z levels.  Vince Larson 12 Jan 2005)

          tmp_zl(1:hoc_nz-1) = 
     .    ( hoc_zl(2:hoc_nz) - les_zl(1:(les_nz-1)) ) / norm_term 

          tmp_zl = tmp_zl**2

          mean_sqr_diff_zt = sum( tmp_zl(1:(hoc_nz-1)), 1 ) 

        case ( 2 )  !  the DYCOMS II RF01 case
          tmp_zl(1:les_nz) = 
     .    ( hoc_zl(3:hoc_nz) - les_zl(1:les_nz) ) / norm_term

          tmp_zl = tmp_zl**2

          mean_sqr_diff_zt = sum( tmp_zl(1:(les_nz)), 1 )

        case default !
          write(0,*) "HOC:", hoc_nz, "LES:", les_nz
          stop "Not able to handle specified number of HOC z-levels"
        end select

        return
        end function mean_sqr_diff_zt
!-----------------------------------------------------------------------
        real function mean_sqr_diff_2_zt
     .                ( hoc_nz, les_nz, hoc_zl, 
     .                  hoc2_zl, les_zl, norm_term )
!       Description:
!       Alternate function to compute mean difference between input 
!       fields.
!       It computes:
!       hoc2_zl - 2 * hoc_zl * les_zl + les_zl**2
!       where hoc2_zl = avg( hoc_zl**2 )
!        hoc_zl  = avg( hoc_zl )
!        les_zl  = avg( les_zl )
!       This alternate formulation adds a penalty to the cost function
!       from the time varying noise that might be present in a simulation.
!       It allows the tuner to avoid very noisy simulations, although some
!       noise might still be present.
!  
!       Configured to do interpolation on LES / HOC comparisons on the 
!       thermodynamic grid

!       References:
!       None
!-----------------------------------------------------------------------
      
        implicit none

        ! External
        intrinsic sum

        ! Input Variables
        integer, intent(in) :: 
     .  hoc_nz, ! Vertical extent for HOC
     .  les_nz  ! Vertical extent for the LES

        real, intent(in), dimension(hoc_nz) :: 
     .  hoc_zl, ! HOC GrADS variable [units vary]
     .  hoc2_zl ! HOC GrADS variable [units vary]

        real, intent(in), dimension(les_nz) :: 
     .  les_zl  ! The LES GrADS variable [units vary]

        real, intent(in) :: 
     .  norm_term ! normalization term;  
                  !typically maxval(les) - minval(les) 

        ! Local Variables
        real, dimension(hoc_nz) :: 
     .  tmp_zl

!----------------------------------------------------------------------

        select case ( hoc_nz - les_nz )
        case ( 0 ) ! most cases
        !  Due to hoc's lower starting point, we can only use
        !  (total number of z-levels) - 1 (a maximum of 74 for BOMEX).
        !  The code below assumes the LES data are on an evenly spaced grid.
        !  (Need to interpolate hoc to LES' levels.  Right now we just
        !   compare adjacent z levels.  Vince Larson 12 Jan 2005)

          tmp_zl(1:hoc_nz-1) 
     .    = ( hoc2_zl(2:hoc_nz) 
     .        - 2.0*hoc_zl(2:hoc_nz)*les_zl(1:(les_nz-1))
     .        + les_zl(1:(les_nz-1))*les_zl(1:(les_nz-1))
     .      ) / (norm_term * norm_term)

          mean_sqr_diff_2_zt   = sum( tmp_zl(1:(hoc_nz-1)), 1 ) 

        case ( 2 )  !  the DYCOMS II RF01 case
          tmp_zl(1:les_nz) 
     .    = ( hoc2_zl(3:hoc_nz) 
     .        - 2.0*hoc_zl(3:hoc_nz)*les_zl(1:les_nz)
     .        + les_zl(1:les_nz)*les_zl(1:les_nz)
     .      ) / (norm_term * norm_term)

          mean_sqr_diff_2_zt = sum( tmp_zl(1:(les_nz)), 1 )

        case default !
          stop "Not able to handle specified number of HOC z-levels"
        end select

        return
        end function mean_sqr_diff_2_zt

!-----------------------------------------------------------------------

        subroutine read_random_seed( seed_file )

!       Description:
!       Reads an ASCII flat file passed as an argument

!       References:
!       None
!-----------------------------------------------------------------------
        use constants, only: fstderr

        implicit none

        ! External
        intrinsic :: random_seed

        ! Input
        character(len=*), intent(in) :: 
     .  seed_file ! This should usually contain > 34 integers
 
        ! Local Variables
        integer(kind=4), dimension(:), allocatable :: 
     .  rand_seed  ! Set of 32 bit integers for seeding the generator

        integer :: 
     .  AllocateStatus,
     .  InputStatus,
     .  rand_size ! For the random seed

!-----------------------------------------------------------------------

      call random_seed( size=rand_size )

      allocate(rand_seed( rand_size ), stat = AllocateStatus)

      if ( AllocateStatus /= 0 ) then 
        stop "Allocation of the random seed variable array failed"
      end if

      ! ASCII formatted file, usually generated by int2txt
      open(unit=30, file=seed_file, action='read')

      read(unit=30, fmt=*, iostat=InputStatus) rand_seed(1:rand_size)
      if ( InputStatus /= 0 ) then 
        write(fstderr,*) "Error reading "//seed_file
        stop
      end if

      close(unit=30)
      
      call random_seed( put=rand_seed )

      deallocate( rand_seed )

      return
      end subroutine read_random_seed

!-----------------------------------------------------------------------

      end module error
!-----------------------------------------------------------------------
