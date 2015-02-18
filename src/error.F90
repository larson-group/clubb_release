!-----------------------------------------------------------------------
! $Id$

module error

! Description:

!   subroutine tuner_init: reads in namelists /stats/, /cases/,
!   /initvars/, & /variance/ from 'error.in'
!   It then uses them to setup the initial param_vals_matrix of independent
!   variables, i.e. the CLUBB constants, and allocate the runtime arrays
!   for each of the model runs and each of the variables

!   function min_les_clubb_diff:  A driver for the CLUBB program/module.
!   Calls run_clubb, reads in les & CLUBB results from GRADS files, and
!   calculates the average difference between the two over all z-levels.

!   subroutine write_results :
!   Prints the results of tuning to the terminal or to a file.

!   subroutine output_nml_tuner :
!   Generates the error.in file using the current constants.

!   subroutine output_nml_standalone :
!   Generates the standalone.in file using the current constants.
!   Standalone CLUBB is only configured to run a single model, and so only
!   the first model is used to make the namelist.

!-----------------------------------------------------------------------
  use parameter_indices, only: nparams ! Variable(s)

  use parameters_tunable, only: read_parameters, read_param_spread ! Procedure(s)

  use mt95, only: genrand_real ! Constant

  use clubb_precision, only: &
    core_rknd ! Variable(s)

  implicit none

  ! Constant Parameters
  integer, parameter, private ::  & 
    max_run = 12,  & ! Maximum model runs the tuner can handle at a time
    max_variables = 32 ! This number / 2 is maximum variables to tune for

  integer, parameter, public :: &
    iamoeba = 0, & ! Numerical Recipes downhill simplex
    iamebsa = 1, & ! Numerical Recipes simulated annealing
    iesa    = 2, & ! ESA tuning type
    iflags  = 3    ! Model flags "tuner"

  ! Variables
  integer, public :: &
    ndim ! Number of variables, e.g. rcm, to be tuned. Dimension of the init simplex

  ! 'err_code' is an important integer used by the min_diff to
  ! determine whether CLUBB has become numerically unstable
  integer, public ::  & 
    err_code

  ! inv_count is a modular counter [1-3] used to determine
  ! which file to output to if l_stdout_on_invalid is true.
  integer, private ::  & 
    inv_count = 0

  !-----------------------------------------------------------------------


  real( kind = core_rknd ), public ::  & 
    f_tol = 1e-5_core_rknd, &    ! The precision to tune for
    anneal_temp = 100._core_rknd ! Initial temperature for the simulated annealing algorithm

  integer, public :: & 
    anneal_iter = 0, &    ! Number of annealing iterations to perform
    tune_type = iesa, & ! Toggle for downhill simplex of simulated annealing
    c_total = 0, &      ! Total number of simulation cases to tune over
    v_total = 0         ! Total number of variables to tune over

  logical, public :: & 
    l_results_stdout = .false., & ! Whether to print tuning results to the terminal
    l_results_file = .false.,   & ! Whether to generate a new error.in based on 
                        ! the new tuning constants
    l_stdout_on_invalid = .false. ! Generate a new error.in when the simulation crashes

  logical, parameter, public :: &
    l_save_tuning_run = .true.  ! If true, writes the results of the tuning run to a file

  character(len=50), public :: &
    tuning_filename = '' ! File where results of tuning run are
                         ! written if l_save_tuning_run = .true.

  integer, parameter, public :: &
    file_unit = 15  ! File unit number connected with tuning_filename


  character(len=10), dimension(:), allocatable, private ::  & 
    hoc_v,  & ! Variables in CLUBB GrADS files
    les_v  ! Variables in LES GrADS files

  integer, dimension(:,:), allocatable, private ::  & 
    time ! Time intervals

  ! Additions for using imposed weights as scaling factors
  logical :: l_initialize_sigma = .true.

  real( kind = core_rknd ), dimension(:,:), allocatable, private ::  & 
    err_terms, & 
    invsigma2, & 
    min_err_terms, & 
    init_err_terms

  real( kind = core_rknd ), dimension(:), allocatable, private ::  & 
    weight_case, weight_var

  ! End additions for using imposed weights

  integer, dimension(:), allocatable, private :: & 
    z_i, & ! Initial z level for tuning purposes
    z_f    ! Final z level for tuning purposes

  character(len=100), dimension(:), allocatable, private ::  & 
    run_file,        & ! Model run files
    hoc_stats_file,  & ! Model GrADS files
    les_stats_file  ! Model GrADS files

  ! Various Variables for returning results
  integer, public :: &
    iter = 0 ! Total number of iterations amoeba spent calculating optimal values

  real( kind = core_rknd ), public :: & 
    init_err = -999._core_rknd,  & ! Error for the initial constants
    min_err  = -999._core_rknd     ! The lowest the minimization algorithm could go

  real( kind = core_rknd ), dimension(nparams), private :: & 
    params = -999._core_rknd  ! Vector of all CLUBB tunable parameter values

  integer, dimension(nparams), private :: & 
    params_index = 0  ! Index of the params elements that are used in the simplex

  real( kind = core_rknd ), allocatable, dimension(:,:), public ::  & 
    param_vals_matrix ! Holds 2D simplex the CLUBB constant parameters
    ! The first row contains the initial values of the tunable parameters.
    ! The remaining rows contain random perturbations of those parameter values.

  real( kind = core_rknd ), allocatable, dimension(:), public :: & 
    param_vals_spread,  & ! Amount to vary each respec. constant by
    cost_fnc_vector    ! cache of differences between the LES and CLUBB

  real(kind=genrand_real), allocatable, dimension(:), private ::  & 
    rand_vect ! A vector of random reals for initializing the x array

  logical, public, dimension(:,:), allocatable :: &
    model_flags_array ! Flags we're checking for


  ! Procedures
  public :: tuner_init, min_les_clubb_diff, write_results, & 
    output_nml_standalone, output_nml_tuner

  private ! Default Scope

  contains

  !-----------------------------------------------------------------------
  subroutine tuner_init( l_read_files )

! Description:
!   Initializes param_vals_matrix with constants from error.in
!   Allocates arrays for cases and tuning variables.
!   Initializes grads file names to read in.

! References:
!   None
!-----------------------------------------------------------------------
    use constants_clubb, only: fstdout, fstderr ! Variable(s)

    use error_code, only: fatal_error ! Procedure(s)

    use text_writer, only: write_text ! Subroutine(s)

    use mt95, only: genrand_real1 ! Procedure

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Constant Variables

    integer, parameter ::  & 
      max_times = 50 ! max number of timesteps to compare (arbitrary)

    character(len=8), parameter :: & 
      filename = "error.in"

    ! Input  Variables

    ! Determines whether to read in the namelists and do the
    ! initial allocation of arrays
    logical, intent(in) ::  & 
      l_read_files

    ! Local Variables

    real( kind = core_rknd ), dimension(nparams) :: rtmp ! Scratch space

    integer :: i, j, & ! looping variables
               iunit ! file unit

    !-----------------------------------------------------------------------

    ! Namelist vars for determining which variable to tune for:

    !  time_nl:   Order pairs of time intervals to analyze
    !  z_i_nl:    initial z-level to begin reading in for tuning
    !  z_f_nl:    final z-level to end reading in for tuning

    integer, dimension(max_run, max_times):: time_nl

    integer, dimension(max_run) :: z_i_nl, z_f_nl

    ! Addition to use imposed weights as scaling factors
    real( kind = core_rknd ), dimension(max_run) :: weight_case_nl

    real( kind = core_rknd ), dimension(max_variables) :: weight_var_nl

    character(len=100), dimension(max_run) ::  & 
      run_file_nl, & 
      hoc_stats_file_nl, & 
      les_stats_file_nl

    character(len=10), dimension(max_variables) :: &
      t_variables ! List of variables to be read from the GrADS output

    ! Namelists read from error.in
    namelist /stats/  & 
      f_tol, tune_type, anneal_temp, anneal_iter, & 
      l_results_stdout, l_results_file, l_stdout_on_invalid, &
      t_variables, weight_var_nl

    namelist /cases/  & 
      les_stats_file_nl, hoc_stats_file_nl, & 
      run_file_nl, z_i_nl, z_f_nl, time_nl, weight_case_nl

    ! Reset iteration counter (set by amoeba)
    iter = 0

    ! Reset invalid run counter
    inv_count = 0

    ! Re-read namelists if requested
    if ( l_read_files ) then

      ! Initialize all compile time arrays to zero
      time_nl = 0

      ! Imposed weights as scaling factors
      weight_case_nl = 0.0_core_rknd
      weight_var_nl  = 0.0_core_rknd

      z_i_nl = 0
      z_f_nl = 0

      ! Initialize variable names to spaces
      t_variables(1:max_variables)  = "          "

      iunit = 10

      ! Open our namelist input file
      open(unit=iunit, file=filename, status='old')

      ! Determine which files to read data from based on namelist
      read(unit=iunit, nml=stats)

      ! Read in the models to be run
      read(unit=iunit, nml=cases)

      ! Close our input namelist file
      close(unit=iunit)

      ! Read in initial constant values
      call read_parameters( iunit, filename, params )

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

      allocate( & 
        z_i(c_total), z_f(c_total), time(c_total, max_times),  & 
        run_file(c_total),  & 
        les_stats_file(c_total), hoc_stats_file(c_total) )

      allocate( & 
        err_terms(c_total, v_total), & 
        min_err_terms(c_total, v_total), & 
        init_err_terms(c_total, v_total), & 
        invsigma2(c_total, v_total), & 
        weight_case(c_total), weight_var(v_total) )

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

      call read_param_spread( iunit, filename, params_index,  & 
                              rtmp, ndim )

      if ( ndim == 0 .and. tune_type /= iflags ) then
        write(fstderr,*) "You must vary at least one parameter"
        stop
      end if

      if ( tune_type == iamoeba .or. tune_type == iamebsa ) then
        ! Numerical recipes simulated annealing or downhill simplex
        allocate( rand_vect(ndim), param_vals_matrix(ndim+1,ndim), & 
                  param_vals_spread(ndim), cost_fnc_vector(ndim+1) )

      else ! ESA algorithm or model flags
        allocate( param_vals_matrix(1,ndim), param_vals_spread(ndim), cost_fnc_vector(1) )
      end if
      ! Initialize the CLUBB parameter spread
      param_vals_spread(1:ndim)  = rtmp(params_index(1:ndim))

      ! Copy tunable parameter values into the first row of the simplex
      param_vals_matrix(1,1:ndim) = params(params_index(1:ndim))

      ! Attempt to generate a pseudo-random seed using a file
      ! generated from /dev/random.  File is an ASCII text file
      ! and can be edited manually.
      call read_random_seed( "rand_seed.dat" )

    end if  ! l_read_files
    !-----------------------------------------------------------------------

    if ( tune_type == iamoeba .or. tune_type == iamebsa ) then
      ! Fill in the remaining values of the array by varying the initial
      ! vector (i.e. the first row of the array) by a small random perturbation
      do j = 1, ndim

        call genrand_real1( rand_vect(1:ndim) )

        do i = 2, ndim+1, 1
          ! Vince Larson made entries of param_vals_matrix random  10 Feb 2005
          ! param_vals_matrix(i,j) = param_vals_matrix(1,j)* &
          !     (1.0_core_rknd+((real(i, kind = core_rknd)-1._core_rknd)/ &
          !     real(ndim, kind = core_rknd)*0.5_core_rknd))
          param_vals_matrix(i,j) = param_vals_matrix(1,j)* & 
          ( (1.0_core_rknd - param_vals_spread(j))  & 
           + real( rand_vect(i-1), kind = core_rknd )*param_vals_spread(j)*2._core_rknd )
          ! End of Vince Larson's change
        end do ! i..ndim+1

      end do ! j..ndim
    end if


    ! First call is used to initialize weights

    l_initialize_sigma = .true.
    cost_fnc_vector(1) =  real(min_les_clubb_diff( &
        real(param_vals_matrix(1,1:ndim)) ), kind = core_rknd)
    l_initialize_sigma = .false.

    ! Note: min_les_clubb_diff is written to deal with undefined and
    ! invalid values for variations on the initial vector, but that
    ! algorithm relies on the initial vector being valid.

    if ( fatal_error( err_code ) ) then
      write(fstderr,*) "Initial variable values must be valid."
      stop
    end if

    ! Save initial error

    init_err = cost_fnc_vector(1)
    init_err_terms = err_terms

    ! Other initialization runs

    if ( tune_type == iamoeba .or. tune_type == iamebsa ) then
      ! Initialize the 'y' vector for amoeba
      ! This is done by calling min_les_clubb_diff with the initial vector
      do i = 2, ndim+1, 1
        cost_fnc_vector(i) =  & 
             real(min_les_clubb_diff( real(param_vals_matrix(i,1:ndim))), kind = core_rknd )
      end do
    end if

    ! Save tuning results in file if specified
    if( l_save_tuning_run ) open(unit=file_unit, file=tuning_filename, &
      action="write", position='append')
    call write_text( "cost_fnc_vector:", l_save_tuning_run, file_unit )
    call write_text( '', cost_fnc_vector, l_save_tuning_run, file_unit )
    if( l_save_tuning_run ) close(unit=file_unit)

    return
  end subroutine tuner_init

  !-----------------------------------------------------------------------
  real function min_les_clubb_diff( param_vals_vector_r4 )

! Description:
!   Function that returns the sum of the error between the dependent
!   variable (i.e. the variable we want to match) in each of the models

! References:
!   _Numerical Recipes in Fortran 77_ P.402-406 (Description)
!   _Numerical Recipes in Fortran 90_ source code (Routine)
!-----------------------------------------------------------------------

    use clubb_driver, only: run_clubb ! Procedure(s)

    use stat_file_utils, only: &
      stat_file_num_vertical_levels, & ! Procedure(s)
      stat_file_vertical_levels, &
      stat_file_average_interval

    use parameters_tunable, only: params_list ! Variable(s)

    use error_code, only: clubb_no_error ! Variable(s)

    use error_code, only: fatal_error ! Procedure(s)

    use numerical_check, only: is_nan_2d ! Procedure(s)

    use text_writer, only: write_text ! Subroutine(s)

    use constants_clubb, only: fstderr ! Constant(s)

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: achar, modulo, minval, maxval, trim, any, real

    ! Constant Parameters
    logical, parameter :: l_stdout = .false.

    ! Input Variables

    ! The interface declaration in the nr module prevents us from making
    ! this a fixed length array declaration
    real, dimension(:), intent(in) ::  & 
      param_vals_vector_r4 ! Tuning vector(ndim dimension) contains
                        ! parameterization constants (C1, C2 etc.)

    ! Local Variables

    real( kind = core_rknd ), dimension(size(param_vals_vector_r4)) :: &
      param_vals_vector ! core_rknd version of param_vals_vector_r4

    real( kind = core_rknd ), dimension(nparams) :: & 
      params_local ! Local copy of the CLUBB parameters fed into run_clubb

    ! These are read after each run from the GrADS control files
    integer ::  & 
      clubb_nz   ! Extent of the CLUBB domain in the z dimension

    character(50) ::  & 
      errorfile ! nml filename for invalid runs

    integer ::  & 
      AllocateStatus ! For hoc_zl, les_zl

    real( kind = core_rknd ) ::  & 
      err_sum  ! scalar sum of all z-levels

    integer ::  & 
      c_terms ! num of terms in err_sum (for normalization)

    ! LES and CLUBB values over nz z-levels
    real( kind = core_rknd ), dimension(:), allocatable ::  & 
      clubb_zl, & 
      clubb2_zl, & 
      les_zl, &
      clubb_grid_heights

    real( kind = core_rknd ), dimension(:,:), allocatable ::  & 
      err_sums  ! To save breakdown of cost function

    real( kind = core_rknd ) ::  & 
      les_minmax ! The largest LES value subtracted from the
    ! smallest value of all zlvl's (for normalization)

    logical ::  & 
      l_error ! Used to det. if reading the variable failed

    integer, dimension(c_total) ::  & 
      run_stat ! isValid over each model case

    integer :: i, j, c_run ! looping variables

    !-----------------------------------------------------------------------

    param_vals_vector = real(param_vals_vector_r4, kind = core_rknd)

    l_error = .false.

    ! Output information every 10 iterations if stdout is enabled;
    ! Amoeba's unusual calling convention makes this happen less
    ! often than might be expected.
    if ( l_results_stdout & 
        .and. ( modulo( iter, 10 ) == 0 ) & 
        .and. iter /= 0  ) then

      write(unit=*,fmt='(A12,I10)') "Iteration: ", iter
      write(unit=*,fmt='(A12)') "Parameters: "

      do i = 1, ndim, 1
        j = params_index(i)
        write(unit=*,fmt='(A30,F27.20)') params_list(j)//" = ", & 
          param_vals_vector(i)
      end do

    end if

    ! Do the same as above if the tuning run is being saved to a file
    if( l_save_tuning_run .and. ( modulo( iter, 10 ) == 0 ) .and. iter /= 0 ) then
      open(unit=file_unit, file=tuning_filename, action='write', position='append')

      write(unit=file_unit,fmt='(A12,I10)') "Iteration: ", iter
      write(unit=file_unit,fmt='(A12)') "Parameters: "

      do i = 1, ndim, 1
        j = params_index(i)
        write(unit=file_unit,fmt='(A30,F27.20)') params_list(j)//" = ", & 
          param_vals_vector(i)
      end do      

      close(unit=file_unit)
    end if ! l_save_tuning_run .and. modulo(iter,10) == 0 .and. iter /= 0

    allocate( err_sums(c_total, v_total) )

    ! Initialize
    err_sum  = 0.0_core_rknd
    err_sums = 0.0_core_rknd
    c_terms  = 0
    err_code = clubb_no_error
    run_stat(1:c_total) = clubb_no_error

    ! Copy simplex into a vector of all possible CLUBB parameters
    do i=1, nparams, 1
      ! If the variable isn't in the simplex, leave it as is
      params_local(i) = params(i)
      do j=1, ndim, 1
        if ( i == params_index(j) ) then
          ! Copy variable from param_vals_vector argument
          params_local(i) = param_vals_vector(j)
          exit
        else
          ! Continue searching the list
          cycle
        end if
      end do
    end do

!-----------------------------------------------------------------------

    ! Cycle through all the model cases specified for CLUBB

    ! OpenMP directives should work as expected now, assuming new
    ! model variables are declared threadprivate -dschanen 31 Jan 2007

!$omp parallel do default(none), private(c_run), &
!$omp   shared(params_local, run_file, run_stat, c_total, model_flags_array, iter)
    do c_run=1, c_total, 1

#ifndef _OPENMP 
      ! Write a message about which case we're calling if OpenMP is not enabled
      ! Save tuning results in file if specified
      if( l_save_tuning_run ) open(unit=file_unit, file=tuning_filename, &
        action='write', position='append')
      call write_text( "Calling CLUBB with case "//trim( run_file(c_run) ), &
        l_save_tuning_run, file_unit )
      if( l_save_tuning_run ) close(unit=file_unit)
      
#endif
      ! Run the CLUBB model with parameters as input

      if ( allocated( model_flags_array ) ) then
        call run_clubb &
             ( params_local, run_file(c_run), l_stdout, run_stat(c_run), &
               model_flags_array(iter,:) )
      else
        call run_clubb & 
             ( params_local, run_file(c_run), l_stdout, run_stat(c_run) )
      end if

    end do ! 1..c_run
!$omp end parallel do

    do c_run = 1, c_total, 1
      if ( err_code < run_stat(c_run) ) err_code = run_stat(c_run)
    end do

    !-----------------------------------------------------------------------

    ! Now check if CLUBB has blown up, i.e. if CLUBB has set a variable to NaN,
    ! or encountered a failure in the matrix solver routines

    ! If it has, it returns higher value than those previous to
    ! Amoeba (the downhill simplex)
    if ( fatal_error( err_code ) ) then
      write(fstderr,*) "Warning: the parameter set has caused CLUBB to crash"
      min_les_clubb_diff = real(2._core_rknd * maxval( cost_fnc_vector )  & 
                       - minval( cost_fnc_vector ))

      if ( l_stdout_on_invalid ) then
        inv_count = modulo( inv_count, 3 ) + 1 ! 1,2,3,1,2,3...
        errorfile = "error_crash_"// achar( inv_count+48 ) // ".in" ! Known magic number
        call output_nml_tuner( errorfile,  & 
                               param_vals_vector(1:ndim) )
      end if

      return
    end if

    !----------------------------------------------------------------------- 

    do c_run=1, c_total, 1

      ! Determine how large the GrADS input is
      clubb_nz = stat_file_num_vertical_levels( hoc_v(1), hoc_stats_file(c_run) )

      ! Allocate the arrays for reading in the GrADS plot data
      allocate( clubb_zl(clubb_nz), clubb2_zl(clubb_nz),  & 
                les_zl(clubb_nz), clubb_grid_heights(clubb_nz), stat=AllocateStatus )

      if ( AllocateStatus /= 0 ) then
        stop "Allocation of arrays in minimization function failed"
      end if

      ! Determine the height of GrADS input
      clubb_grid_heights = stat_file_vertical_levels( hoc_v(1), hoc_stats_file(c_run), clubb_nz )

      ! Start with first CLUBB & LES variables, then loop through and
      ! calculate the mean squared difference for all the variables
      do i=1, v_total, 1  

        ! Read in LES grads data for one variable, averaged
        ! over specified time intervals
        les_zl =  & 
        stat_file_average_interval &
        ( les_stats_file(c_run), clubb_nz,  & 
          time(c_run,:), les_v(i), clubb_grid_heights, 1, l_error )

        if ( l_error ) then
          if( l_save_tuning_run ) then
            open(unit=file_unit, file=tuning_filename, &
            action='write', position='append')
            call write_text( "LES variable was invalid or the GrADS file does not exist", &
              l_save_tuning_run, file_unit )
            close(unit=file_unit)
          end if
          write(fstderr,*) "The specified LES variable "//trim( les_v(i) )//" was invalid, "// &
            "or the GrADS file did not exist."
          stop
        end if

        ! Verify that the domain that we're tuning CLUBB over is fully defined in
        ! the LES data.  If not, some points will be NaN
        if ( is_nan_2d( les_zl(z_i(c_run):z_f(c_run)) ) ) then
          write(*,*)
          write(fstderr,*) "The tuning domain exceeds the size of the LES data, "// &
            "or the LES data is NaN"
          write(*,*)
          write(fstderr,*) trim( les_v(i) )//" = ", les_zl(z_i(c_run):z_f(c_run))
          stop "Fatal error"
        end if

        ! Read in CLUBB grads data for one variable, averaged
        ! over specified time intervals
        clubb_zl =  & 
        stat_file_average_interval & 
        ( hoc_stats_file(c_run), clubb_nz,  & 
          time(c_run,:), hoc_v(i), clubb_grid_heights, 1, l_error )

        if ( l_error ) then
          stop "The specified CLUBB variable was invalid"
        end if

        ! The same variable, with npower = 2
        clubb2_zl =  & 
        stat_file_average_interval & 
        ( hoc_stats_file(c_run), clubb_nz, & 
          time(c_run,:), hoc_v(i), clubb_grid_heights, 2, l_error )

        if ( l_error ) then
          stop "The specified CLUBB variable was invalid"
        end if

        !-----------------------------------------------------------------------

        ! Calculate the mean squared difference between the CLUBB
        ! and the LES variables

        ! In order to deal with differences in order of magnitude
        ! between the variables, the err_sum equation has been
        ! modified to normalize the values with respect to the
        ! the minimum and maximum in the LES. -Dave Schanen

        les_minmax = maxval( les_zl(z_i(c_run):z_f(c_run)) ) &
          - minval( les_zl(z_i(c_run):z_f(c_run)) )

        if ( les_minmax == 0.0_core_rknd ) then
          stop "An LES variable was 0 from z_i to z_f."
        end if

        ! Old code
!     err_sum = err_sum &
!      + mean_sqr_diff_zt( clubb_nz, clubb_zl, les_zl, les_minmax )

        ! Chris Golaz modification: mean_sqr_diff_2 was designed to try
        ! to limit time noise in tuning simulations.
        ! New code with Modification for weighting
        err_sums(c_run,i) = mean_sqr_diff_2( clubb_nz, z_i(c_run), z_f(c_run), clubb_zl,  & 
                              clubb2_zl, les_zl, les_minmax )

        c_terms = c_terms + 1

      end do ! i=1..v_total

      ! De-allocate the arrays for reading in the GrADS plot data
      deallocate( clubb_zl, clubb2_zl, les_zl, clubb_grid_heights )

    end do     ! end of do c_run=1, c_total

!----------------------------------------------------------------------

    ! Return error averaged over all cases, variables,
    ! and vertical levels
    ! Old Code
!       min_les_clubb_diff = err_sum / real( c_terms )

    !---------------------------------------------------------------
    ! Compute normalization factors to satisfy imposed weights
    ! This non-dimensionlizes each term in the cost function so 
    ! that the units in the variables(eg. rcm, cloud_frac) don't matter
    ! Invsigma2 is computed once and then reused in all later 
    ! calculations because invsigma2 is declared in the module
    !---------------------------------------------------------------
    if ( l_initialize_sigma ) then
      do c_run=1,c_total
        do i=1,v_total
          invsigma2(c_run,i)  & 
          = weight_case(c_run)*weight_var(i) / err_sums(c_run,i)
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
    min_les_clubb_diff = real(err_sum)

    deallocate( err_sums )

    ! Save tuning results in file if specified
    if( l_save_tuning_run ) open(unit=file_unit, file=tuning_filename, &
      action='write', position='append')
    call write_text( "Cost function= ", real(min_les_clubb_diff, kind = core_rknd), &
      l_save_tuning_run, file_unit, '(a, f12.5)' )
    if( l_save_tuning_run ) close(unit=file_unit)

    return
  end function min_les_clubb_diff

  !----------------------------------------------------------------------
  subroutine write_results( iunit )

! Description:
!   Outputs the results of a tuning run to the terminal if iunit = fstdout
!   or to a file which has already been opened and connected with
!   integer iunit.

! References:
!   None
!----------------------------------------------------------------------

    use parameters_tunable, only: params_list ! Variable(s)

    implicit none

    ! Input variables
    integer, intent(in) :: iunit ! = fstdout for printing to the terminal
                                 ! or = a number connected with a file

    ! Local variables
    integer :: i ! Loop iterator

    if ( tune_type == 0 ) then
      write(unit=iunit,fmt=*) "Number of iterations past initialization:",  iter
    end if ! tune_type == 0

    write(unit=iunit,fmt='(4x,A9,5x,10x,A7,10x,10x,A7)') & 
        "Parameter", "Initial", "Optimal"

    do i = 1, ndim, 1
      write(unit=iunit,fmt='(A30,2F27.20)')  & 
        params_list(params_index(i))//" = ",  & 
        params(params_index(i)), param_vals_matrix(1,i)
    end do

    write(unit=iunit,fmt='(A20)') "Initial cost: "
    write(unit=iunit,fmt='(F15.6)') init_err
    write(unit=iunit,fmt='(A20)') "Optimal cost: "
    ! The $$ is here to make it easy to find with grep
    write(unit=iunit,fmt='(A3,F15.6)') "$$ ", min_err

    write(unit=iunit,fmt=*) "Approx. percent increase in accuracy:",  & 
      ((init_err - min_err) / init_err*100.0_core_rknd), "%" ! Known magic number

    return
  end subroutine write_results
  !-----------------------------------------------------------------------

  subroutine output_nml_tuner( results_f, param_vals_vector )

! Description:
!   Output namelists to a formatted text file

! References:
!   None

! Notes:
!   You can do the same thing with
!   write(unit=<UNIT>,nml=<NAMELIST>), but at the time I wrote this
!   I was more ambitious and didn't like that fact that it came out
!   in all caps on pgf90. -dschanen 28 July 2006
!-----------------------------------------------------------------------

    use parameters_tunable, only: params_list ! Variable(s)

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: achar, trim

    ! Parameter Constants
    character(len=1), parameter :: dbqt = '"'

    ! Input Variables
    character(len=*), intent(in) ::  & 
    results_f ! Name of the results file to write to

    real( kind = core_rknd ), intent(in), dimension(:) ::  & 
    param_vals_vector ! The current tuning parameters

    ! Local Variables
    real( kind = core_rknd ), dimension(nparams) :: params_local

    integer :: i, j, & ! loop variable
               iunit   ! file unit

    character :: i_c   ! loop variable in ASCII
    
    iunit = 20

    ! Open a new file
    open(unit=iunit, file=results_f,  & 
         action="write", access="sequential")

    ! Write variables to results file
    ! All this is based on the previous namelists in error.in,
    ! except for the constants parameters for CLUBB

    write(unit=iunit,fmt=*) "! Parameter file " // results_f
    write(unit=iunit,fmt=*) "&stats"
    write(unit=iunit,fmt=*) "f_tol = ", f_tol
    if (l_results_stdout) then
      write(unit=iunit,fmt=*) "l_results_stdout = " // ".true."
    else
      write(unit=iunit,fmt=*) "l_results_stdout = " // ".false."
    end if
    if (l_results_file) then
      write(unit=iunit,fmt=*) "l_results_file = " // ".true."
    else
      write(unit=iunit,fmt=*) "l_results_file = " // ".false."
    end if
    if (l_stdout_on_invalid) then
      write(unit=iunit,fmt=*) "l_stdout_on_invalid = " // ".true."
    else
      write(unit=iunit,fmt=*) "l_stdout_on_invalid = " // ".false."
    end if
    write(unit=iunit,fmt=*) "t_variables = "
    do i=1, v_total
      write(unit=iunit,fmt=*) dbqt, hoc_v(i), dbqt, ",",  & 
        dbqt, les_v(i), dbqt, ","
    end do
    write(unit=iunit,fmt=*) "weight_var_nl = ", weight_var
    write(unit=iunit,fmt=*) "anneal_temp = " , anneal_temp
    write(unit=iunit,fmt=*) "anneal_iter = " , anneal_iter
    write(unit=iunit,fmt=*) "tune_type   = " , tune_type
    write(unit=iunit,fmt=*) "/"

    write(unit=iunit,fmt=*) "&cases"
    do i = 1, c_total
      i_c = achar( i + 48 )
      write(unit=iunit,fmt=*) "les_stats_file_nl("// i_c //") = ",  & 
        dbqt, trim( les_stats_file(i) ), dbqt
      write(unit=iunit,fmt=*) "hoc_stats_file_nl("// i_c //") = ",  & 
        dbqt, trim( hoc_stats_file(i) ), dbqt
      write(unit=iunit,fmt=*) "run_file_nl("// i_c //") = ",  & 
        dbqt, trim( run_file(i) ), dbqt
      write(unit=iunit,fmt=*) "z_i_nl("// i_c //")  = " , z_i(i)
      write(unit=iunit,fmt=*) "z_f_nl("// i_c //")  = " ,  z_f(i)
      write(unit=iunit,fmt=*) "weight_case_nl("// i_c //") = " ,  & 
        weight_case(i)
      write(unit=iunit,fmt=*) "time_nl("// i_c //",:)  = " , time(i,:)
    end do
    write(unit=iunit,fmt=*) "/"

    write(unit=iunit,fmt=*) "&initvars"

    ! Copy simplex into a vector of all possible CLUBB parameters
    do i=1, nparams, 1
      ! If the variable isn't in the simplex, leave it as is
      params_local(i) = params(i)
      do j=1, ndim, 1
        if ( i == params_index(j) ) then
          ! Copy variable from param_vals_vector argument
          params_local(i) = param_vals_vector(j)
          exit
        else
          cycle
        end if
      end do
    end do

    ! Output optimal values and all possible CLUBB parameters
    do i=1, nparams, 1
      write(unit=iunit,fmt='(A30,F27.20)') & 
        trim( params_list(i) )//" = ", params_local(i)
    end do
    write(unit=iunit,fmt=*) "/"

    write(unit=iunit,fmt=*) "&initspread"
    ! Copy the spread into a vector of all possible CLUBB parameters
    do i=1, nparams, 1
      ! If the variable isn't being changed, set it to zero
      params_local(i) = 0.0_core_rknd
      do j=1, ndim, 1
        if ( i == params_index(j) ) then
          ! Copy variable from param_vals_vector argument
          params_local(i) = param_vals_spread(j)
          exit
        else
          cycle
        end if
      end do
    end do

    ! Output the amount each variable was changed for the simplex
    do i=1, nparams, 1
      write(unit=iunit,fmt='(a30,f12.5)') & 
        trim( params_list(i) )//" = ", params_local(i)
    end do
    write(unit=iunit,fmt=*) "/"

    ! Close new namelist file
    close(unit=iunit)

    return
  end subroutine output_nml_tuner

!-----------------------------------------------------------------------
  subroutine output_nml_standalone ( results_f,  & 
                                     param_vals_vector )

! Description:
!   Output namelists to a formatted file

! References:
!   None

! Notes:
!   See note for output_nml_tuner, above. dschanen 28 July 2006
!-----------------------------------------------------------------------

    use parameters_tunable, only: params_list ! Variable(s)

    implicit none

    ! External
    intrinsic :: achar, trim

    ! Input variables
    character(len=*), intent(in) :: & 
      results_f ! Results file to write to

    real( kind = core_rknd ), intent(in), dimension(ndim) :: & 
      param_vals_vector ! the current constants

    ! Local variables
    real( kind = core_rknd ), dimension(nparams) :: params_local

    integer   :: i, j  ! loop variables

    integer :: iunit ! file unit

    ! ---- Begin Code ---
    write(6,*) "Writing parameters for namelist to "//results_f

    iunit = 20
    ! Open new file
    open(unit=iunit, file=results_f,  & 
         action="write", access="sequential")

    ! Write variables to namelist for standalone CLUBB.
    ! All this is based on the previous error.in, except the constants

    write(unit=iunit,fmt=*) "! Parameter file " // results_f

    write(unit=iunit,fmt=*) "&initvars"
    ! Copy simplex into a vector of all possible CLUBB parameters
    do i=1, nparams, 1
      ! If the variable isn't in the simplex, leave it as is
      params_local(i) = params(i)
      do j=1, ndim, 1
        if ( i == params_index(j) ) then
          ! Copy variable from param_vals_vector argument
          params_local(i) = param_vals_vector(j)
          exit
        else
          cycle
        end if
      end do
    end do

    ! Output optimal values and all possible CLUBB parameters
    do i=1, nparams, 1
      write(unit=iunit,fmt='(A30,F27.20)') & 
        trim( params_list(i) )//" = ", params_local(i)
    end do
    write(unit=iunit,fmt=*) "/"

    ! Close new file
    close(unit=iunit)

    return
  end subroutine output_nml_standalone

!-----------------------------------------------------------------------
  real( kind = core_rknd ) function mean_sqr_diff & 
                ( nz, z_init, z_final, scm_zl, les_zl, norm_term )
! Description
!   Calculate the mean squared difference between two input vectors,
!   then normalize.

! References:
!   None
!-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: sum

    ! Input Variables
    integer, intent(in) ::  & 
      nz,      & ! Vertical extent for CLUBB
      z_init,  & ! Initial point for the purposes of computing the sum
      z_final    ! Final point for the purpose of computign the sum

    real( kind = core_rknd ), intent(in), dimension(nz) ::  & 
      scm_zl, &! CLUBB GrADS variable   [units vary]
      les_zl   ! The LES GrADS variable [units vary]

    real( kind = core_rknd ), intent(in) ::  & 
      norm_term ! normalization term; typically maxval(les) - minval(les)

    real( kind = core_rknd ), dimension(nz) ::  & 
      tmp_zl

    integer :: k

!----------------------------------------------------------------------

    ! ---- Begin Code ----

    tmp_zl = 0.0_core_rknd

    do k = z_init, z_final, 1
      tmp_zl(k) = ( ( scm_zl(k) - les_zl(k) ) / norm_term )**2
    end do

    mean_sqr_diff = sum( tmp_zl )

    return
  end function mean_sqr_diff
!-----------------------------------------------------------------------
  real( kind = core_rknd ) function mean_sqr_diff_2 & 
                ( nz, z_init, z_final, scm_zl,  & 
                  scm2_zl, les_zl, norm_term )
! Description:
!   Alternate function to compute mean difference between input
!   fields.
!   It computes:
!     scm2_zl - 2 * scm_zl * les_zl + les_zl**2
!     where scm2_zl = avg( scm_zl**2 )
!      scm_zl  = avg( scm_zl )
!      les_zl  = avg( les_zl )
!   This alternate formulation adds a penalty to the cost function
!   from the time varying noise that might be present in a simulation.
!   It allows the tuner to avoid very noisy simulations, although some
!   noise might still be present.
!
!   Configured to do interpolation on LES / CLUBB comparisons on the
!   CLUBB grid

! References:
!   None
!-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: sum

    ! Input Variables
    integer, intent(in) ::  & 
      nz,      & ! Vertical extent for CLUBB
      z_init,  & ! Initial point for the purposes of computing the sum
      z_final    ! Final point for the purpose of computign the sum

    real( kind = core_rknd ), intent(in), dimension(nz) ::  & 
      scm_zl,  & ! CLUBB GrADS variable [units vary]
      scm2_zl, & ! CLUBB GrADS variable [units vary]
      les_zl     ! The LES GrADS variable [units vary]

    real( kind = core_rknd ), intent(in) ::  & 
      norm_term ! normalization term. Typically maxval(les) - minval(les)

    ! Local Variables
    real( kind = core_rknd ), dimension(nz) ::  & 
      tmp_zl

    integer :: k

!----------------------------------------------------------------------

    ! ---- Begin Code ----

    tmp_zl = 0.0_core_rknd

    do k = z_init, z_final, 1
      tmp_zl(k) = ( scm2_zl(k) - 2.0_core_rknd * scm_zl(k)*les_zl(k) + les_zl(k)**2 ) &
        / norm_term**2
    end do

    mean_sqr_diff_2 = sum( tmp_zl )

    return
  end function mean_sqr_diff_2

!-----------------------------------------------------------------------

  subroutine read_random_seed( seed_file )

! Description:
!   Reads an ASCII flat file passed as an argument

! References:
!   None
!-----------------------------------------------------------------------
    use constants_clubb, only: fstderr ! Constant

    use mt95, only: genrand_init ! Procedure

    use mt95, only: genrand_intg ! Constant

    implicit none

    ! Parameter constants
    integer, parameter :: rand_size = 34

    ! Input
    character(len=*), intent(in) ::  & 
      seed_file ! This should usually contain >= 34 integers

    ! Local Variables
    integer(kind=genrand_intg), dimension(rand_size) ::  & 
      rand_seed  ! Set of 32 bit integers for seeding the generator

    integer :: & 
      InputStatus, iunit

!-----------------------------------------------------------------------

    iunit = 30

    ! ASCII formatted file, usually generated by int2txt
    open(unit=iunit, file=seed_file, action='read')

    read(unit=iunit, fmt=*, iostat=InputStatus) rand_seed(1:rand_size)
    if ( InputStatus /= 0 ) then
      write(fstderr,*) "Error reading "//seed_file
      stop
    end if

    close(unit=iunit)

    call genrand_init( put=rand_seed )

    return
  end subroutine read_random_seed

!-----------------------------------------------------------------------

end module error
!-----------------------------------------------------------------------
