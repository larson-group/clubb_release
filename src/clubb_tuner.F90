!----------------------------------------------------------------------
! $Id$

program clubb_tuner 

!     Description:
!     ``Tunes'' constants in clubb so that the output matches LES output.
!     Uses amoeba or amebsa to calculate the min of (f_les - f_clubb)^2

!     References:
!     _Numerical Recipes in Fortran 90_ (Chapter 10) 
!     (Amoeba & Amebsa subroutine)
!----------------------------------------------------------------------
  use error, only:  & 
    tuner_init, min_les_clubb_diff,                & ! Subroutines 
    write_results,                                 & ! Subroutine
    output_nml_standalone, output_nml_tuner,       & ! Subroutines
    param_vals_matrix, anneal_temp,                & ! Variables
    l_results_stdout, l_save_tuning_run,           & ! Variables
    l_results_file, tune_type, f_tol, ndim,        & ! Variables
    tuning_filename, file_unit                       ! Variable

  use error, only:  & 
    iamoeba, & ! Constants
    iamebsa, &
    iesa, &
    iflags

  use constants_clubb, only: & 
    fstdout ! Variable

  use text_writer, only: &
    write_text ! Subroutine

  implicit none

  ! External
  external :: enhanced_simann_driver, amoeba_driver, amebsa_driver, &
    logical_flags_driver

  ! Variables
  character(len=10) :: current_time  ! Current time string (no seconds)
  character(len=8)  :: current_date  ! Current date string
  character(len=75) :: results_f     ! Results file

  character(len=3)  :: user_response ! Simple Y/N query

  !----------------- Begin Code -------------------

  ! Determine date and time to be used with file names
  call date_and_time( current_date, current_time )

  ! Create file to save tuning results if specified
  if ( l_save_tuning_run ) then
    ! File where tuning run results are written
    tuning_filename = "../input/tuning_run_results_"// &
      current_date//'_'//current_time(1:4)//".log"
    open(unit=file_unit, file=tuning_filename, action="write")
    write(file_unit,*) "Tuning..."
    close(unit=file_unit)
  end if ! l_save_tuning_run

  ! Read in namelists and define parameters
  call tuner_init( l_read_files=.true. )

  ! Attempt to find the optimal parameter set
  do
    select case ( tune_type )
    case ( iamoeba )
      call amoeba_driver( )

    case ( iamebsa )
      call amebsa_driver( )

    case ( iesa )
      call enhanced_simann_driver( )

    case ( iflags )
      call logical_flags_driver( current_date, current_time )
      stop "Program exited normally"

    case default
       stop "Unknown tuning type"
    end select

    ! Print to stdout if specified
    if ( l_results_stdout ) call write_results( fstdout )

    ! Save results in file if specified
    if ( l_save_tuning_run ) then
      open(unit=file_unit, file=tuning_filename, action="write", position="append")
      call write_results( file_unit )
    end if ! l_save_tuning_run

    ! Query to see if we should exit the loop
    call write_text( "Run Complete.", l_save_tuning_run, file_unit )

    write(unit=fstdout,fmt='(A)', advance='no')  & 
      "Re-run with new parameters?(y/n) "
    read(*,*) user_response
    ! Save tuning results in file if specified
    if( l_save_tuning_run ) then
      write(file_unit,*) "Re-run with new parameters?(y/n) ", user_response
      close(unit=file_unit)
    end if ! l_save_tuning_run

    select case ( trim( user_response ) )
    case ( "Yes", "yes", "y", "Y" )
      continue
    case default
      exit
    end select

    ! Save tuning results in file if specified
    if ( tune_type == iamoeba .or. tune_type == iamebsa ) then

      if( l_save_tuning_run ) open(unit=file_unit, file=tuning_filename, &
        action="write", position="append")

      call write_text( "Current f_tol= ", f_tol, l_save_tuning_run, file_unit )

      write(fstdout,fmt='(A)', advance='no') "Enter new f_tol=   "
      read(*,*) f_tol
      ! Save tuning results in file if specified
      if ( l_save_tuning_run ) write(file_unit,*) "Enter new f_tol=   ", f_tol
      if ( l_save_tuning_run ) close(unit=file_unit)
    end if
    if ( tune_type == iesa .or. tune_type == iamebsa ) then
      write(fstdout,fmt='(A)', advance='no') "New annealing temp =   "
      read(*,*) anneal_temp
    end if

    call tuner_init( l_read_files=.false. )

  end do ! user_response /= 'y', 'Y' or 'yes'

  ! Final namelist file output 

  if ( l_results_file ) then 

    ! Tuner namelist
    ! Save tuning results in file if specified
    if( l_save_tuning_run ) open(unit=file_unit, file=tuning_filename, &
      action="write", position="append")
    call write_text( "Generating new error.in file...", l_save_tuning_run, file_unit )
    if( l_save_tuning_run ) close(unit=file_unit)

    results_f = "../input/error_"//current_date//'_' & 
      //current_time(1:4)//".in" 

    ! Note:
    ! The first column of param_vals_matrix is the optimized result, which
    ! is swapped in by amoeba.
    call output_nml_tuner( results_f, param_vals_matrix(1,1:ndim) )
    ! Save tuning results in file if specified
    if( l_save_tuning_run ) open(unit=file_unit, file=tuning_filename, &
      action="write", position="append")
    call write_text( "New filename is: "//results_f, l_save_tuning_run, file_unit )

    ! Parameters namelist
    call write_text( "Generating new tunable_parameters.in file...", l_save_tuning_run, file_unit )
    if( l_save_tuning_run ) close(unit=file_unit)

    results_f = "../input/tunable_parameters/tunable_parameters_"//current_date//'_' & 
      //current_time(1:4)//".in" 

    call output_nml_standalone( results_f,  & 
                                param_vals_matrix(1,1:ndim) )
    ! Save tuning results in file if specified
    if( l_save_tuning_run ) open(unit=file_unit, file=tuning_filename, &
      action="write", position="append")
    call write_text( "New filename is: "//results_f, l_save_tuning_run, file_unit )
    if( l_save_tuning_run ) close(unit=file_unit)

  end if

  ! Print message if tuning results were saved in a file
  if ( l_save_tuning_run ) then
    print*, "***The tuning results have been saved in the file: "
    print*, "  "//tuning_filename
  end if ! l_save_tuning_run

  ! Exit Program

  stop "Program exited normally"

  end program clubb_tuner

  !------------------------------------------------------------------------
  subroutine amoeba_driver

  !     Description:
  !     Simple interface for the amoeba minimization algorithm

  !     References:
  !     _Numerical Recipes in Fortran 90_.  See full citation above.
  !------------------------------------------------------------------------
#ifdef TUNER
  use nr, only:  & 
      amoeba ! Procedure(s)
  use error, only:  & 
    ! Variable(s)
    ndim,                                & ! Array dimensions
    param_vals_matrix, cost_fnc_vector,  & ! The 'p' matrix and 'y' vector resp.
    f_tol,                                & ! Tolerance of tuning run
    iter,                                & ! Iteration number
    min_les_clubb_diff,                  & ! Cost function
    min_err                                ! Minimum value of the cost function

  use clubb_precision, only: &
    core_rknd ! Variable(s)

  implicit none

  real, dimension(1:ndim+1, 1:ndim) :: &
    param_vals_matrix_r4

  real, dimension(1:ndim+1) :: &
    cost_fnc_vector_r4

  ! ---- Begin Code ----

  param_vals_matrix_r4 = real(param_vals_matrix(1:ndim+1,1:ndim))
  cost_fnc_vector_r4 = real(cost_fnc_vector(1:ndim+1))

  call amoeba( param_vals_matrix_r4,  & 
               cost_fnc_vector_r4,  & 
               real(f_tol), min_les_clubb_diff, iter)

  param_vals_matrix(1:ndim+1,1:ndim) = real(param_vals_matrix_r4, kind = core_rknd)
  cost_fnc_vector(1:ndim+1) = real(cost_fnc_vector_r4, kind = core_rknd)


  ! Note:
  ! Amoeba will make the optimal cost result the first element of
  ! cost_fnc_vector and the optimal parameter set the first column of 
  ! param_vals_matrix, where param_vals_matrix is 'p' in the NR subroutine
  ! and cost_fnc_vector is 'y'.

  min_err = cost_fnc_vector(1)

  return
#else
  stop "Numerical recipes subroutines were disabled at compile time"
#endif
  end subroutine amoeba_driver

  !-----------------------------------------------------------------------

  subroutine amebsa_driver

  ! Description:
  !   Interface for the amoeba simulated annealing minimization algorithm
  !   At the end of the subroutine, the param_vals_matrix's first row gets the
  !   optimal values assigned to it.
  !
  ! References:
  !   None
  !-----------------------------------------------------------------------
#ifdef TUNER
  use nr, only:  & 
      amebsa ! Procedure(s)
  use nrtype, only:  & 
      SP ! Variable(s)
  use error, only: & 
      param_vals_matrix,  & ! Variable(s)
      anneal_temp, & 
      anneal_iter, & 
      iter, & 
      ndim, & 
      cost_fnc_vector, & 
      f_tol, & 
      min_err, & 
      min_les_clubb_diff ! Procedure(s)

  use clubb_precision, only: &
      core_rknd ! Variable(s)

  implicit none

  ! Note:
  ! Most of these are taken from xamebsa, so I only have a vague idea
  ! what the their purpose is and can't pick less obscure names.
  ! dschanen 1 Apr 05
   
  ! Local Variables
  integer ::  & 
    iiter, jiter, & ! Loop variables
    nit ! ???

  real( kind = core_rknd ), dimension(ndim) ::  & 
    pb ! ???
   
  real( kind = core_rknd ) ::  & 
    ybb,   & ! ???
    yb,    & ! ???
    tmptr ! ???

  ! Temporary variables for passing into subroutines that require reals
  real, dimension(1:ndim+1, 1:ndim) :: &
    param_vals_matrix_r4

  real, dimension(1:ndim+1) :: &
    cost_fnc_vector_r4

  real, dimension(1:ndim) :: &
    pb_r4

  real :: &
    yb_r4

  ! ---- Begin Code ----

  ybb   = 1.0e30_core_rknd
  yb    = 1.0e30_core_rknd
  nit   = 0
  iiter = 0
  tmptr = anneal_temp ! anneal_temp taken from /stat/ namelist

  param_vals_matrix_r4 = real(param_vals_matrix(1:ndim+1,1:ndim))
  cost_fnc_vector_r4 = real(cost_fnc_vector(1:ndim+1))

  do jiter = 1, anneal_iter ! anneal_iter taken from /stat/ namelist
    iter  = iiter
    tmptr = tmptr * 0.8_core_rknd

    pb_r4 = real(pb(1:ndim))
    yb_r4 = real(yb)

    call amebsa & 
         ( param_vals_matrix_r4,  & 
           cost_fnc_vector_r4, & 
           pb_r4, yb_r4, real(f_tol), min_les_clubb_diff, iter, real(tmptr) )

    pb(1:ndim) = real(pb_r4, kind = core_rknd)
    yb = real(yb_r4, kind = core_rknd)

    nit = nit + iiter - iter
    if ( yb < ybb ) then
      ybb = yb
    end if
    if ( iter > 0 ) exit
  end do

  param_vals_matrix(1:ndim+1,1:ndim) = real(param_vals_matrix_r4, kind = core_rknd)
  cost_fnc_vector(1:ndim+1) = real(cost_fnc_vector_r4, kind = core_rknd)

  param_vals_matrix(1,1:ndim) = pb(1:ndim)
  min_err = ybb

  return

#else
  stop "Numerical recipes subroutines were disabled at compile time"
#endif
end subroutine amebsa_driver
!----------------------------------------------------------------------
subroutine enhanced_simann_driver

  ! Description:
  !   Wrapper subroutine for the ESA driver

  ! References: 
  !   ``Enhanced Simulated Annealing for Many Globally Minimized Functions
  !   of Many-Continuous Variables'', Siarry, et al. ACMS TOMS Vol. 23,
  !   No. 2, June 1997, pp. 209--228.
  !------------------------------------------------------------------------

  use enhanced_simann, only: esa_driver ! Procedure(s)

  use error, only: & ! Variable(s)
    ndim,               & ! Array dimensions
    param_vals_matrix,  & ! The 'p' matrix
    param_vals_spread,  & ! Used here for the initial value of rostep
    anneal_temp,        & ! Start annealing temperature
    min_err               ! Minimum value of the cost function

  use error, only:  & ! Procedure(s)
    min_les_clubb_diff  ! Cost function

  use clubb_precision, only: &
    core_rknd ! Variable(s)

  implicit none

  ! Local Variables

  real( kind = core_rknd ), dimension(ndim) :: &
    xinit,  & ! Initial values for the tunable parameters
    x0min,  & ! Minimum value for the tunable parameters
    x0max,  & ! Maximum value for the tunable parameters
    rostep, & ! Initial step size
    xopt      ! Final values for the tunable parameters

  real( kind = core_rknd ) :: enopt ! Optimal cost

  ! ---- Begin Code ----

  xinit = param_vals_matrix(1,1:ndim)

  ! Set the minimum for the parameters.  Assume no parameter is < 0 for now
  x0min(1:ndim) = 0._core_rknd 

  ! Set the maximum for the parameters.  Assume parameters will be most 5 
  ! times the current value.
  x0max(1:ndim) = 5._core_rknd * param_vals_matrix(1,1:ndim) 

  rostep(1:ndim) = param_vals_spread(1:ndim)

  call esa_driver( xinit, x0min, x0max, rostep, & ! In
                   anneal_temp, min_les_clubb_diff, & ! In/out, Function
                   xopt, enopt ) ! Out

  param_vals_matrix(1,1:ndim) = xopt(1:ndim)
  min_err = enopt

  return
end subroutine enhanced_simann_driver
!-------------------------------------------------------------------------------
subroutine logical_flags_driver( current_date, current_time )

! Description:
!   While not a true search algorithm in same sense as the simulated annealing
!   or the downhill simplex method is, this driver will try all permutations.
!
! References:
!   None
!-------------------------------------------------------------------------------
  use error, only: &
    param_vals_matrix, & ! Variable(s)
    model_flags_array, &
    iter, &
    l_results_file, &
    l_results_stdout

  use error, only: &
    min_les_clubb_diff ! Procedure(s)

  use constants_clubb, only: &
    fstdout ! Constant(s)

  use quicksort, only: &
    Qsort_flags ! Procedure(s)

  use model_flags, only: &
    get_configurable_model_flags, & ! Procedure(s)
    setup_configurable_model_flags, &
    write_model_flags_to_file

  use clubb_precision, only: &
    core_rknd ! Variable(s)

  implicit none

  ! External
  intrinsic :: btest, selected_int_kind, int

  ! Constant parameters
  integer, parameter :: &
    i8 = selected_int_kind( 15 )

  integer, parameter :: &
    ndim = 11, & ! Temporarily hardwired for a fixed number of flags
    two_ndim = 2**ndim, &
    iunit = 10

  ! Input Variables
  character(len=10), intent(in) :: current_time  ! Current time string (no seconds)
  character(len=8), intent(in)  :: current_date  ! Current date string

  ! Local Variables

  real( kind = core_rknd ), dimension(two_ndim) :: &
    cost_function  ! Values from the cost function

  real( kind = core_rknd ), dimension(ndim) :: &
    cost_func_sum_true,  & ! Averaged cost function when the flag is true
    cost_func_sum_false, & ! Averaged cost function when the flag is false
    cost_func_avg          ! Averaged cost function true - false.

  logical, dimension(ndim) :: model_flags_default

  character(len=128) :: &
    filename_nml, &  ! Namelist file name
    filename_csv     ! Comma seperated values filename

  integer(kind=i8) :: bit_string, bit_iter
  real( kind = core_rknd) :: cost_func_default

  integer :: i, j

  ! ---- Begin Code ----

  ! Determine the current flags
  call get_configurable_model_flags( model_flags_default(1),  &
                                model_flags_default(2),  &
                                model_flags_default(3),  &
                                model_flags_default(4),  &
                                model_flags_default(5),  &
                                model_flags_default(6),  &
                                model_flags_default(7),  &
                                model_flags_default(8),  &
                                model_flags_default(9),  & 
                                model_flags_default(10), &
                                model_flags_default(11) )

  ! This should always be 1.0; it's here as a sanity check
  cost_func_default = real( min_les_clubb_diff( real(param_vals_matrix(1,:)) ), kind = core_rknd )

  allocate( model_flags_array(two_ndim,ndim) )
  bit_string = 0_i8 ! Initialize bits to 00 ... 00
  do i = 1, two_ndim
    do j = 1, ndim
     ! This loop sets 1:n logicals using individual bits, i.e. 0 means
     ! false and 1 means true for the purposes of trying all possibilities
     bit_iter = int( j, i8 )
     model_flags_array(i,j) = btest( bit_string, bit_iter-1_i8 ) 
    end do
    bit_string = bit_string + 1_i8 ! Increment the binary adder
  end do

  ! Compute the cost function with new set of flags.  The model_flags array is
  ! indexed using the iter variable in min_les_clubb_diff to avoid having to
  ! modify the Numerical Recipes code.
  do iter = 1, two_ndim
    cost_function(iter) = real( min_les_clubb_diff( real(param_vals_matrix(1,:)) ), &
         kind = core_rknd )
  end do

  ! Compute a metric of false cost function - true cost function
  cost_func_sum_true = 0.0_core_rknd
  cost_func_sum_false = 0.0_core_rknd
  do i = 1, two_ndim
    do j = 1, ndim
      if ( model_flags_array(i,j) ) then ! Flag is true
        cost_func_sum_true(j) = cost_func_sum_true(j) + cost_function(i)
      else ! Flag is false
        cost_func_sum_false(j) = cost_func_sum_false(j) + cost_function(i)
      end if
    end do
  end do
  cost_func_avg(:) = ( cost_func_sum_false(:) / real( two_ndim / 2, kind = core_rknd ) ) &
                   - ( cost_func_sum_true(:) / real( two_ndim / 2, kind = core_rknd ) )

  ! Sort flags and the cost function in ascending order
  call Qsort_flags( model_flags_array, cost_function )

  if ( l_results_stdout ) then
    ! Output results to the terminal
    write(fstdout,'(A30)') "Default flags:                "
    do j = 1, ndim
      write(fstdout,'(L6,4X)',advance='no') model_flags_default(j)
    end do
    write(fstdout,'(G10.3)') cost_func_default
    write(fstdout,'(A)') "Results from trying all permutations of the flags: "
    do i = 1, ndim
      write(fstdout,'(I6,4X)',advance='no') i
    end do
    write(fstdout,'(A10)') " Cost func" 
    do i = 1, ndim
      write(fstdout,'(G10.3)',advance='no') cost_func_avg(i)
    end do
    write(fstdout,*) ""
    do i = 1, two_ndim
      do j = 1, ndim
        write(fstdout,'(L6,4X)',advance='no') model_flags_array(i,j)
      end do
      write(fstdout,'(G10.3)') cost_function(i)
    end do

    write(fstdout,'(A30)') &
      "Column 1 = upwind_wpxp_ta     ", &
      "Column 2 = upwind_xpyp_ta     ", &
      "Column 3 = upwind_xm_ma       ", &
      "Column 4 = quintic_poly_interp", &
      "Column 5 = vert_avg_closure   ", &
      "Column 6 = single_C2_Skw      ", &
      "Column 7 = standard_term_ta   ", &
      "Column 8 = tke_aniso          ", &
      "Column 9 = use_cloud_cover    "
  end if ! l_results_stdout

  ! Generate CSV file of the results
  filename_csv = "../output/clubb_model_flags_"//current_date//"_"//current_time//".csv"
  open(unit=iunit,file=trim( filename_csv ))
  write(iunit,'(10A20)') "upwind_wpxp_ta, ", "upwind_xpyp_ta, ", "upwind_xm_ma, ", &
    "quintic_poly_interp, ", "vert_avg_closure, ", &
    "single_C2_Skw, ", "standard_term_ta, ", "tke_aniso, ", "use_cloud_cover, ", "Cost func."
  write(iunit,'(A30)') "Default flags:               ,"
  do j = 1, ndim
    write(iunit,'(L20,A2)',advance='no') model_flags_default(j), ", "
  end do
  write(iunit,'(G20.6,A2)') cost_func_default, ", "
  do i = 1, ndim
    write(iunit,'(G20.6,A2)',advance='no') cost_func_avg(i), ", "
  end do
  write(iunit,'(A20)') "Avg false - Avg true ,"
  do i = 1, two_ndim
    do j = 1, ndim
      write(iunit,'(L20,A2)',advance='no') model_flags_array(i,j), ", "
    end do
    write(iunit,'(G20.6,A2)') cost_function(i), ", "
  end do
  close(unit=iunit)
  write(fstdout,*) "Results of tuning model flags written to: ", trim( filename_csv )

  ! Generate namelist file of the optimal result
  if ( l_results_file ) then
   call setup_configurable_model_flags( model_flags_array(1,1), &
                                   model_flags_array(1,2), &
                                   model_flags_array(1,3), &
                                   model_flags_array(1,4), &
                                   model_flags_array(1,5), &
                                   model_flags_array(1,6), &
                                   model_flags_array(1,7), &
                                   model_flags_array(1,8), &
                                   model_flags_array(1,9), & 
                                   model_flags_array(1,10),& 
                                   model_flags_array(1,11)  )

    filename_nml = "../input/tunable_parameters/configurable_model_flags_"//current_date//'_' & 
      //current_time(1:4)//".in"

    call write_model_flags_to_file( iunit, trim( filename_nml ) )

    write(fstdout,*) "New namelist of tuning model flags written to: ", trim( filename_nml )
  end if

  deallocate( model_flags_array )

  return
end subroutine logical_flags_driver
