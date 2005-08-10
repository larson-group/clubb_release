!-----------------------------------------------------------------------
! $Id: error.F,v 1.20 2005-08-10 22:58:52 griffinb Exp $

!       Module ERROR
!
!       Contains:

!       subroutine tuner_init: reads in namelists /stats/, /cases/, 
!          /initvars/, & /variance/ from 'error.in' 
!       It then uses them to setup the initial x_array of independent 
!       variables, i.e. the hoc constants, and allocate the runtime arrays
!       for each of the model runs and each of the variables

!       function funk:  A driver for the hoc program/module.
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
      module error
      
      implicit none

!---------------- Dimensions of the array for amoeba -------------------

      integer, parameter, public  :: ndim = 24
      integer, parameter, public  :: mdim = ndim + 1

      real, dimension(mdim, ndim)  :: x_array
      real, dimension(ndim)        :: x_variance
      real, dimension(mdim)        :: err_array
      real, dimension(ndim)        :: pb           ! only used by amebsa

!----------------------------------------------------------------------- 
!  isValid is an important inout boolean used by funk to determine 
!  whether hoc has blown up 
!  inv_count is a modular counter [1-3] used to determine which file to 
!    output to if stdout_on_invalid is true.
!----------------------------------------------------------------------- 
      logical :: isValid           
 
      integer :: inv_count

!----------------------------------------------------------------------- 
!  max_run: Maximum models the tuner can handle.
!  max_variables: 
!           This number / 2 is maximum variables to tune for
!  hoc_v:   variable names in hoc grads file.  Run-time array
!  les_v:   variable names in les grads file.  Run-time array
!  ftol:    the accuracy to tune for.  Note that a small increase in ftol
!           will result in an exponential increase in computational time.
!  v_total: Total number of variables to tune. Calculated from namelist.
!  c_total: Total number of cases to tune. Calculated from namelist.
!----------------------------------------------------------------------- 
      
      integer, parameter, private       :: max_run       = 12
      integer, parameter, private       :: max_variables = 32

      real                              :: ftol 
      real                              :: anneal_temp

      integer                           :: anneal_iter
      integer                           :: tune_type 
      integer                           :: c_total, v_total

      logical                           :: results_stdout, results_file
      logical                           :: stdout_on_invalid

!----- Hoc constants to tune -----------------------------------------
      real C1, C2rt, C2thl, C2rtthl, C4, C5
      real C6rt, C6thl, C7, C8, C10, C11
      real nu1, nu2, nu6, nu8, nu_r
      real beta
      real gamma_coef, c_K
      real lmin_coef, taumin, taumax
      real mu

!---- Runtime arrays which vary according to cases and variables ------
      character(len=10), dimension(:), allocatable :: hoc_v, les_v

      integer, dimension(:,:), allocatable         :: time
      integer, dimension(:), allocatable           :: z_i, z_f

      character(len=50), dimension(:), allocatable :: run_file
      character(len=50), dimension(:), allocatable :: hoc_stats_file
      character(len=50), dimension(:), allocatable :: les_stats_file

!----- Various Variables for returning results ----------------------
      integer :: iter           ! Total number of iterations amoeba
                                ! spent calculating optimal values

      real    :: init_err       ! error for the initial constants
      real    :: min_err        ! the lowest the minimization algorithm could go

      contains
!-----------------------------------------------------------------------
!     Subroutine TUNER_INIT
!     Parameters: logical read_files
!     Initializes x_array with constants from error.in
!     Allocates arrays for cases and tuning variables 
!     Initializes grads file names to read in 
!-----------------------------------------------------------------------  
      subroutine tuner_init( read_files )

      implicit none

!     Input 
      logical, intent(in)   :: read_files ! whether to read in  nml's and
                                          ! allocate arrays
!     Internal
      integer, parameter    :: max_times = 10 ! max number of timesteps 
                                              ! to compare
      real, dimension(mdim) :: rand_vect      ! vect. of random reals

      integer i, j           ! looping variables
      integer AllocateStatus ! for allocatable arrays

      
!----------------------------------------------------------------------- 
!  Namelist vars for determining which variable to tune for:

!  time_nl:   Order pairs of time intervals to analyze
!  z_i_nl:    initial z-level to begin reading in for tuning
!  z_f_nl:    final z-level to end reading in for tuning

      integer, dimension(max_run, max_times):: time_nl 
      integer, dimension(max_run)           :: z_i_nl, z_f_nl

      character(len=50), dimension(max_run) :: run_file_nl
      character(len=50), dimension(max_run) :: hoc_stats_file_nl
      character(len=50), dimension(max_run) :: les_stats_file_nl

      character(len=10), dimension(max_variables) :: t_variables


!----- Hoc variance for generating array ------------------------------
      real C1_var, C2rt_var, C2thl_var, C2rtthl_var, C4_var, C5_var
      real C6rt_var, C6thl_var, C7_var, C8_var, C10_var, C11_var
      real nu1_var, nu2_var, nu6_var, nu8_var, nu_r_var
      real beta_var
      real gamma_coef_var, c_K_var
      real lmin_coef_var, taumin_var, taumax_var
      real mu_var

!----- Namelists read from error.in ---------------------------------
      namelist /stats/ ftol, tune_type, anneal_temp, anneal_iter,
     .                 results_stdout, results_file, stdout_on_invalid,
     .                 t_variables

      namelist /cases/ les_stats_file_nl, hoc_stats_file_nl, 
     .                 run_file_nl,
     .                 z_i_nl, z_f_nl, time_nl

      namelist /initvars/ C1, C2rt, C2thl, C2rtthl, C4, C5, 
     .                    C6rt, C6thl, C7, C8, C10, C11, 
     .                    nu1, nu2, nu6, nu8, nu_r, beta, gamma_coef,
     .                    c_K, lmin_coef, taumin, taumax, mu

      namelist /variance/ C1_var, C2rt_var, C2thl_var, C2rtthl_var, 
     .                    C4_var, C5_var, C6rt_var, C6thl_var, 
     .                    C7_var, C8_var, C10_var, C11_var, 
     .                    nu1_var, nu2_var, nu6_var, nu8_var, nu_r_var,
     .                    beta_var, gamma_coef_var, c_K_var,
     .                    lmin_coef_var, taumin_var, taumax_var, mu_var

!------------------------- Reset iteration counter ---------------------
      iter      = 0

!------------------------- Reset invalid run counter -------------------
      inv_count = 0

!---------------------- Re-read namelists if requested -----------------
      if ( read_files ) then

!----------------- Initialize all compile time arrays to zero ----------
        time_nl       = 0 
        z_i_nl        = 0 
        z_f_nl        = 0

!----------------- Initialize variable names to spaces ---------------- 
        t_variables(1:max_variables)  = "          "

!----- Open our namelist input file -----------------------------------
        open(10, file='error.in', status='old')

!----- Determine which files to read data from based on namelist ------     
        read(10, nml = stats)

!------------ Read in the models to be run ----------------------------
        read(10, nml = cases)

!--------------- Read in initial constant values ----------------------
        read(10, nml = initvars)

!--------------- Read in initial variance for constants ---------------
        read(10, nml = variance)

!---- Close our input namelist file -----------------------------------
        close(10)

!----- Allocate the arrays for the tuning variables -------------------

        do i=1, max_variables, 2         ! 1, 3, 5, 7
          if (t_variables(i) == "          ") exit
          v_total = (i + 1) / 2
        enddo 

        allocate( hoc_v(v_total), les_v(v_total), stat=AllocateStatus )

        if ( AllocateStatus /= 0 ) then 
          stop "Allocation of the tuning variable arrays failed"
        endif

!----------- Allocate the arrays for the run cases --------------------

        do i=1, max_run  
          if (z_f_nl(i) == 0 ) exit
          c_total = i
        enddo

        allocate( z_i(c_total), z_f(c_total), 
     .            time(c_total, max_times), 
     .            run_file(c_total), les_stats_file(c_total),
     .            hoc_stats_file(c_total),
     .            stat=AllocateStatus )

        if ( AllocateStatus /= 0 ) then 
          stop "Allocation of the case information arrays failed"
        endif

!----- Transfer the variable numbers to hoc_v and les_v ---------------
       do i=1, v_total
         hoc_v(i) = t_variables( i*2 - 1 )
         les_v(i) = t_variables( i*2 )
       enddo 

!------ Transfer the case information to run-time arrays --------------
       do i=1, c_total 
         z_i(i)               = z_i_nl(i)
         z_f(i)               = z_f_nl(i)
         les_stats_file(i)    = les_stats_file_nl(i)
         hoc_stats_file(i)    = hoc_stats_file_nl(i)
         run_file(i)          = run_file_nl(i)
         time(i, 1:max_times) = time_nl(i,1:max_times)
       enddo

!--------------- Initialize the hoc constants -------------------------

        x_array(1,1:ndim) = (/C1, C2rt, C2thl, C2rtthl, C4, C5, 
     .                        C6rt, C6thl, C7, C8, C10, C11, 
     .                        nu1, nu2, nu6, nu8, nu_r, beta,
     .                        gamma_coef, c_K, lmin_coef, taumin,
     .                        taumax, mu/)

!--------------- Initialize the hoc variance -------------------------

        x_variance(1:ndim) = (/C1_var, C2rt_var, C2thl_var, C2rtthl_var,
     .                         C4_var, C5_var, C6rt_var, C6thl_var,
     .                         C7_var, C8_var, C10_var, C11_var, 
     .                         nu1_var, nu2_var, nu6_var, nu8_var, 
     .                         nu_r_var, beta_var, gamma_coef_var, 
     .                         c_K_var, lmin_coef_var, taumin_var, 
     .                         taumax_var, mu_var/)
      
!-----------------------------------------------------------------------
!  Attempt to generate a pseudo-random seed using a file generated from
!  /dev/random.  File is an ASCII text file and can be edited manually.
!-----------------------------------------------------------------------
        call read_random_seed( "rand_seed.dat" )

      endif  ! read_files

!-----------------------------------------------------------------------
!    Fill in the remaining values of the array by varying the initial
!    vector (i.e. the first column of the array) by a small multiple
!-----------------------------------------------------------------------

      do j = 1, ndim
        call random_number( rand_vect(1:mdim) )
        do i = 2, mdim
! Vince Larson made entries of x_array random  10 Feb 2005
!          x_array(i,j) = x_array(1,j)*
!     .                   (1.0+((real(i)-1.)/real(ndim)*0.5))
          x_array(i,j) = x_array(1,j) *( ( 1.0-x_variance(j) ) 
     .                 + rand_vect(i) * x_variance(j)*2 )
! End of Vince Larson's change
        enddo
      enddo

!-----------------------------------------------------------------------
!     Initialize the 'y' vector for amoeba 
!     this is done by calling funk with the initial vector 
!-----------------------------------------------------------------------
      do i = 1, mdim
        err_array(i) = funk(x_array(i,1:ndim))
      enddo 

      write(*,*) "err_array:"
      write(*,'(6e12.5)')  err_array

      return
      end subroutine tuner_init

!-----------------------------------------------------------------------
!  function FUNK
!  Function that returns the sum of the error between the dependent
!  variable (i.e. the variable we want to match) in each of the models
!-----------------------------------------------------------------------
      real function funk (x_vector)
 
      use hoc
      use grads_common

      implicit none

!     Input 
!     real, intent(in), dimension(ndim) :: x_vector ! tuning vector
      real, dimension(:), intent(in) :: x_vector ! tuning vector

!     Internal
      integer                           :: i, j, c_run    ! looping variables
      integer                           :: les_nz, hoc_nz ! read dynamically from GRaDS 
      
      character(50)                     :: errorfile ! nml filename for 
                                                     ! invalid runs

      integer                           :: AllocateStatus ! for hoc_zl, les_zl

      real                              :: err_sum  ! scalar sum of all z-levels
      integer                           :: c_terms  ! num of terms in err_sum

!-------- les and hoc values over nz z-levels --------------------------
      real, dimension(:), allocatable   :: hoc_zl, les_zl

!--------- largest value minus smallest value  of all zlvl -------------
      real                              :: les_minmax

      logical                           :: error ! see if get_var failed

!-----------------------------------------------------------------------
! Output information every 10 iterations if stdout is enabled
! Amoeba's unusual calling convention makes this work less often than
! would be expected.
!-----------------------------------------------------------------------
      if ( results_stdout
     .     .and. modulo(iter, 10) == 0 
     .     .and. iter /= 0  ) then
        write(*,'(A12,I10)') "Iteration: ", iter
        write(*,'(A12)') "Constants: "
        write(*,'(A12,f12.5)') "C1 =", x_vector(1)
        write(*,'(A12,f12.5)') "C2rt =", x_vector(2)
        write(*,'(A12,f12.5)') "C2thl =", x_vector(3)
        write(*,'(A12,f12.5)') "C2rtthl =", x_vector(4)
        write(*,'(A12,f12.5)') "C4 =", x_vector(5)
        write(*,'(A12,f12.5)') "C5 =", x_vector(6)
        write(*,'(A12,f12.5)') "C6rt =", x_vector(7)
        write(*,'(A12,f12.5)') "C6thl =", x_vector(8)
        write(*,'(A12,f12.5)') "C7 =", x_vector(9)
        write(*,'(A12,f12.5)') "C8 =", x_vector(10)
        write(*,'(A12,f12.5)') "C10 =", x_vector(11)
        write(*,'(A12,f12.5)') "C11 =", x_vector(12)
        write(*,'(A12,f12.5)') "nu1 =", x_vector(13)
        write(*,'(A12,f12.5)') "nu2 =", x_vector(14)
        write(*,'(A12,f12.5)') "nu6 =", x_vector(15)
        write(*,'(A12,f12.5)') "nu8 =", x_vector(16)
        write(*,'(A12,f12.5)') "nu_r =", x_vector(17)
        write(*,'(A12,f12.5)') "beta =", x_vector(18)
        write(*,'(A12,f12.5)') "gamma_coef =", x_vector(19)
        write(*,'(A12,f12.5)') "c_K =", x_vector(20)
        write(*,'(A12,f12.5)') "lmin_coef =", x_vector(21)
        write(*,'(A12,f12.5)') "taumin =", x_vector(22)
        write(*,'(A12,f12.5)') "taumax =", x_vector(23)
        write(*,'(A12,f12.5)') "mu =", x_vector(24)
      endif
        
!----Initialize the return value for funk -----------------------------
      err_sum = 0.

!----Initialize the number of terms in the error sum-------------------
      c_terms = 0

!-----------------------------------------------------------------------
!  Cycle through all the model cases specified for HOC

!  If the compiler supports OpenMP, it will try to run the model cases 
!    in parallel.
!-----------------------------------------------------------------------
!$OMP PARALLEL
!$OMP DO PRIVATE(hoc_zl, les_zl, les_minmax, run_file, hoc_stats_file) 
!   . SHARED(c_terms, err_sum)
      do c_run=1, c_total

!-------------- Run the HOC model with x_vector as input ---------------

!      print *, "Calling HOC with case " // run_file ! Debug

        call hoc_model( x_vector(1), x_vector(2), x_vector(3),
     .                  x_vector(4), x_vector(5), x_vector(6), 
     .                  x_vector(7), x_vector(8), x_vector(9),
     .                  x_vector(10), x_vector(11), x_vector(12), 
     .                  x_vector(13), x_vector(14), x_vector(15),
     .                  x_vector(16), x_vector(17), x_vector(18), 
     .                  x_vector(19), x_vector(20), x_vector(21),
     .                  x_vector(22), x_vector(23), x_vector(24),
     .                  run_file(c_run), isValid, .false. )

!-------------- Basic Sanity Check for negative constants --------------
!        if (minval( x_vector ) < 0. ) isValid = .false.

!-----------------------------------------------------------------------
!  Checks if hoc has blown up, i.e. if hoc has set a variable to NaN, 
!  or encountered a failure in subroutines tridag or bandec.  
!  If it has, then return to amoeba with higher value than previously.
!-----------------------------------------------------------------------
        if ( .not. isValid ) then
          funk = 2 * maxval( err_array ) - minval( err_array )
          if ( stdout_on_invalid ) then
            inv_count = modulo(inv_count, 3) + 1          ! 1,2,3,1,2,3...
            errorfile = "error_crash_" // achar( inv_count+48 )// ".in"
            call output_nml_tuner( errorfile, x_vector(1:ndim) )
          endif
          return
        endif

!-----------------------------------------------------------------------
!  Allocate the arrays for reading in the grads plot data
!-----------------------------------------------------------------------
        hoc_nz = grads_zlvl( hoc_stats_file(c_run) )
        les_nz = grads_zlvl( les_stats_file(c_run) )
        allocate( hoc_zl(hoc_nz), les_zl(les_nz), 
     .            stat = AllocateStatus )
        if ( AllocateStatus /= 0 ) then 
          stop "Allocation of z-level arrays in function funk failed"
        endif

!-----------------------------------------------------------------------
!  Start with first hoc & les variable, loop through and calculate the 
!  mean squared difference for all the variables
!-----------------------------------------------------------------------
        do i=1, v_total

!----- Read in LES grads data for one variable, avgd over intervals ----
      
          les_zl = grads_average_interval( les_stats_file(c_run), 
     .                                     les_nz,
     .                                     time(c_run,:), les_v(i),
     .                                     error )
          if ( error ) stop "oops, variable was invalid"

!----- Read in HOC grads data for one variable, avgd over intervals ----
          hoc_zl = grads_average_interval( hoc_stats_file(c_run),
     .                                     hoc_nz,
     .                                     time(c_run,:), hoc_v(i), 
     .                                     error )
          if ( error ) stop "oops, variable was invalid"

!-----------------------------------------------------------------------
!  Calculate the mean squared difference between HOC and LES variables.

!  In order to deal with order of magnitude differences between
!  variables, the err_sum equation has been modified to normalize the
!  the values with respect to the maxval in les_zl.  dschanen 31 Jan 05
!-----------------------------------------------------------------------

          les_minmax = maxval(les_zl) - minval(les_zl)
          if (les_minmax == 0.) stop "An LES var = 0 over all z-levels"

          err_sum = err_sum 
     .            + mean_sqr_diff_zt( hoc_zl, les_zl, les_minmax )
          c_terms = c_terms + 1

        enddo     ! end of do i=1, v_total

        deallocate (hoc_zl, les_zl)

      enddo     ! end of do c_run=1, c_total
!$OMP END PARALLEL

!----------------------------------------------------------------------
! funk returns error avgd over all cases, variables, and vertical levels
!----------------------------------------------------------------------
      funk = err_sum / real( c_terms )

      return

      end function funk

!--------------------Output results to STDOUT --------------------------
      subroutine output_results_stdout( )
      
      implicit none

      if (tune_type == 0) print *, "Number of iterations:",  iter
      
      write(*,'(A12,A12,A12)') "Var", "Initial ", "Optimal"
      write(*,'(A12,f12.5,f12.5)') "C1 =", C1, x_array(1,1)
      write(*,'(A12,f12.5,f12.5)') "C2rt =", C2rt, x_array(1,2)
      write(*,'(A12,f12.5,f12.5)') "C2thl =", C2thl, x_array(1,3)
      write(*,'(A12,f12.5,f12.5)') "C2rtthl =", C2rtthl, x_array(1,4)
      write(*,'(A12,f12.5,f12.5)') "C4 =", C4, x_array(1,5)
      write(*,'(A12,f12.5,f12.5)') "C5 =", C5, x_array(1,6)
      write(*,'(A12,f12.5,f12.5)') "C6rt =", C6rt, x_array(1,7)
      write(*,'(A12,f12.5,f12.5)') "C6thl =", C6thl, x_array(1,8)
      write(*,'(A12,f12.5,f12.5)') "C7 =", C7, x_array(1,9)
      write(*,'(A12,f12.5,f12.5)') "C8 =", C8, x_array(1,10)
      write(*,'(A12,f12.5,f12.5)') "C10 =", C10, x_array(1,11)
      write(*,'(A12,f12.5,f12.5)') "C11 =", C11, x_array(1,12)
      write(*,'(A12,f12.5,f12.5)') "nu1 =", nu1, x_array(1,13)
      write(*,'(A12,f12.5,f12.5)') "nu2 =", nu2, x_array(1,14)
      write(*,'(A12,f12.5,f12.5)') "nu6 =", nu6, x_array(1,15)
      write(*,'(A12,f12.5,f12.5)') "nu8 =", nu8, x_array(1,16)
      write(*,'(A12,f12.5,f12.5)') "nu_r =", nu_r, x_array(1,17)
      write(*,'(A12,f12.5,f12.5)') "beta =", beta, x_array(1,18)
      write(*,'(A12,f12.5,f12.5)') "gamma_coef =", 
     .                             gamma_coef, x_array(1,19)
      write(*,'(A12,f12.5,f12.5)') "c_K =", c_K, x_array(1,20)
      write(*,'(A12,f12.5,f12.5)') "lmin_coef =", 
     .                             lmin_coef, x_array(1,21)
      write(*,'(A12,f12.5,f12.5)') "taumin =", taumin, x_array(1,22)
      write(*,'(A12,f12.5,f12.5)') "taumax =", taumax, x_array(1,23)
      write(*,'(A12,f12.5,f12.5)') "mu =",  mu,   x_array(1,24)
      
      write(*,'(A20,e15.6)') "initial error: ", init_err

      write(*,'(A20,e15.6)') "optimal error: ", min_err
     
      write(*,*) "Approximate percent increase in accuracy:", 
     .          ((init_err - min_err)/init_err * 100.), "%"

      return
      end subroutine output_results_stdout

!----------------- Output namelists to a formatted file ----------------
      subroutine output_nml_tuner( results_f, x_vector )

      implicit none

      character(len=*), intent(in) :: results_f ! Results file to write to

      real, intent(in), dimension(ndim) :: x_vector ! the current constants
    
      integer                   :: i       ! loop variable
      character                 :: i_c     ! loop variable character conversion

!------- Open new file -------------------------------------------------
      open(unit = 20, file=results_f, action="write",
     .     access="sequential")
      
!-----------------------------------------------------------------------
!  Write variables to results file 
!  All this is based on the previous error.in, except the constants
!-----------------------------------------------------------------------

      write(20,*) "! Parameter file " // results_f
      write(20,*) "&stats"
      write(20,*) "ftol = " , ftol
      if (results_stdout) then
          write(20,*) "results_stdout = " // ".true."
      else
          write(20,*) "results_stdout = " // ".false."
      endif
      if (results_file) then
          write(20,*) "results_file = " // ".true."
      else
          write(20,*) "results_file = " // ".false."
      endif
      if (stdout_on_invalid) then
          write(20,*) "stdout_on_invalid = " // ".true."
      else
          write(20,*) "stdout_on_invalid = " // ".false."
      endif
      write(20, *) "t_variables = "
      do i=1, v_total
        write(20, *) achar(34), hoc_v(i), achar(34), ",", 
     .               achar(34), les_v(i), achar(34), "," 
      enddo
      write(20,*) "anneal_temp = " , anneal_temp
      write(20,*) "anneal_iter = " , anneal_iter
      write(20,*) "tune_type   = " , tune_type
      write(20,*) "/"

      write(20,*) "&cases"
      do i = 1, c_total
        i_c = achar( i + 48 )
        write(20,*) "les_stats_file_nl("// i_c //") = ", 
     .              achar(34), trim( les_stats_file(i) ), achar(34)
        write(20,*) "hoc_stats_file_nl("// i_c //") = ", 
     .              achar(34), trim( hoc_stats_file(i) ), achar(34)
        write(20,*) "run_file_nl("// i_c //") = ", 
     .              achar(34), trim( run_file(i) ), achar(34)
        write(20,*) "z_i_nl("// i_c //")  = " , z_i(i)
        write(20,*) "z_f_nl("// i_c //")  = " ,  z_f(i)
        write(20,*) "time_nl("// i_c //",:)  = " , time(i,:)
      enddo
      write(20,*) "/"
      
      write(20,*) "&initvars"
      write(20,*) "C1 =         ", x_vector(1)
      write(20,*) "C2rt =       ", x_vector(2)
      write(20,*) "C2thl =      ", x_vector(3)
      write(20,*) "C2rtthl =    ", x_vector(4)
      write(20,*) "C4 =         ", x_vector(5)
      write(20,*) "C5 =         ", x_vector(6)
      write(20,*) "C6rt =       ", x_vector(7)
      write(20,*) "C6thl =      ", x_vector(8)
      write(20,*) "C7 =         ", x_vector(9)
      write(20,*) "C8 =         ", x_vector(10)
      write(20,*) "C10 =        ", x_vector(11)
      write(20,*) "C11 =        ", x_vector(12)
      write(20,*) "nu1 =        ", x_vector(13)
      write(20,*) "nu2 =        ", x_vector(14)
      write(20,*) "nu6 =        ", x_vector(15)
      write(20,*) "nu8 =        ", x_vector(16)
      write(20,*) "nu_r =       ", x_vector(17)
      write(20,*) "beta =       ", x_vector(18)
      write(20,*) "gamma_coef = ", x_vector(19)
      write(20,*) "c_K =        ", x_vector(20)
      write(20,*) "lmin_coef =  ", x_vector(21)
      write(20,*) "taumin =     ", x_vector(22)
      write(20,*) "taumax =     ", x_vector(23)
      write(20,*) "mu =         ", x_vector(24)
      write(20,*) "/"

      write(20,*) "&variance"
      write(20,*) "C1_var =         ", x_variance(1)
      write(20,*) "C2rt_var =       ", x_variance(2)
      write(20,*) "C2thl_var =      ", x_variance(3)
      write(20,*) "C2rtthl_var =    ", x_variance(4)
      write(20,*) "C4_var =         ", x_variance(5)
      write(20,*) "C5_var =         ", x_variance(6)
      write(20,*) "C6rt_var =       ", x_variance(7)
      write(20,*) "C6thl_var =      ", x_variance(8)
      write(20,*) "C7_var =         ", x_variance(9)
      write(20,*) "C8_var =         ", x_variance(10)
      write(20,*) "C10_var =        ", x_variance(11)
      write(20,*) "C11_var =        ", x_variance(12)
      write(20,*) "nu1_var =        ", x_variance(13)
      write(20,*) "nu2_var =        ", x_variance(14)
      write(20,*) "nu6_var =        ", x_variance(15)
      write(20,*) "nu8_var =        ", x_variance(16)
      write(20,*) "nu_r_var =       ", x_variance(17)
      write(20,*) "beta_var =       ", x_variance(18)
      write(20,*) "gamma_coef_var = ", x_variance(19)
      write(20,*) "c_K_var =        ", x_variance(20)
      write(20,*) "lmin_coef_var =  ", x_variance(21)
      write(20,*) "taumin_var =     ", x_variance(22)
      write(20,*) "taumax_var =     ", x_variance(23)
      write(20,*) "mu_var =         ", x_variance(24)
      write(20,*) "/"

!------- Close new file ------------------------------------------------
      close(20)

      return
      end subroutine output_nml_tuner

!-------------------Output namelists to a formatted file ---------------
      subroutine output_nml_standalone ( results_f, x_vector )

      implicit none

!     Input
      character(len=*), intent(in) :: results_f ! Results file to write to

      real, intent(in), dimension(ndim) :: x_vector ! the current constants

!     Internal    
      integer                   :: i       ! loop variable
      character                 :: i_c     ! loop variable character conversion

!------- Open new file -------------------------------------------------
      open(unit = 20, file=results_f, action="write",
     .     access="sequential")
      
!-----------------------------------------------------------------------
!  Write variables to namelist for standalone hoc. 
!  All this is based on the previous error.in, except the constants
!-----------------------------------------------------------------------

      write(20,*) "! Parameter file " // results_f

      write(20,*) "&model"
      write(20,*) "run_file  = " // achar(34) // "../tune/" //
     .             trim( run_file(1) ) // achar(34)
      write(20,*) "stdout   = ", ".true."
      write(20,*) "/"
      
      write(20,*) "&initvars"
      write(20,*) "C1 =         ", x_vector(1)
      write(20,*) "C2rt =       ", x_vector(2)
      write(20,*) "C2thl =      ", x_vector(3)
      write(20,*) "C2rtthl =    ", x_vector(4)
      write(20,*) "C4 =         ", x_vector(5)
      write(20,*) "C5 =         ", x_vector(6)
      write(20,*) "C6rt =       ", x_vector(7)
      write(20,*) "C6thl =      ", x_vector(8)
      write(20,*) "C7 =         ", x_vector(9)
      write(20,*) "C8 =         ", x_vector(10)
      write(20,*) "C10 =        ", x_vector(11)
      write(20,*) "C11 =        ", x_vector(12)
      write(20,*) "nu1 =        ", x_vector(13)
      write(20,*) "nu2 =        ", x_vector(14)
      write(20,*) "nu6 =        ", x_vector(15)
      write(20,*) "nu8 =        ", x_vector(16)
      write(20,*) "nu_r =       ", x_vector(17)
      write(20,*) "beta =       ", x_vector(18)
      write(20,*) "gamma_coef = ", x_vector(19)
      write(20,*) "c_K =        ", x_vector(20)
      write(20,*) "lmin_coef =  ", x_vector(21)
      write(20,*) "taumin =     ", x_vector(22)
      write(20,*) "taumax =     ", x_vector(23)
      write(20,*) "mu =         ", x_vector(24)
      write(20,*) "/"

!------- Close new file ------------------------------------------------
      close(20)

      return
      end subroutine output_nml_standalone
!-----------------------------------------------------------------------
!  FUNCTION mean_sqr_diff_zt
!  Calculate the mean squared difference between two input vectors, then
!  normalize.
!  Configured to do interpolation on LES / HOC comparisons on the t grid
!-----------------------------------------------------------------------
      real function mean_sqr_diff_zt( hoc_zl, les_zl, norm_term )
      
      implicit none

!     Input
      real, intent(in), dimension(:) :: hoc_zl 
      real, intent(in), dimension(:) :: les_zl 

      real, intent(in) :: norm_term ! normalization term; 
                                    ! typically minval(les) - maxval(les) 
!     Internal
      real, dimension(size( hoc_zl )) :: tmp_zl
      integer nz 

      nz = size( hoc_zl ) 

      select case ( nz )
      case (48, 75, 110, 147, 150, 160) ! most cases
!----------------------------------------------------------------------
!  Due to hoc's lower starting point, we can only use
!  (total number of z-levels) - 1 (a maximum of 74 for BOMEX).
!  The code below assumes the LES data are on an evenly spaced grid.
!  (Need to interpolate hoc to LES' levels.  Right now we just
!   compare adjacent z levels.  Vince Larson 12 Jan 2005)
!----------------------------------------------------------------------
        tmp_zl(1:nz-1) = ( hoc_zl(2:nz) - les_zl(1:(nz-1)) )/norm_term 
        tmp_zl         = tmp_zl**2
        mean_sqr_diff_zt = sum( tmp_zl(1:nz-1), 1 ) 

      case (132)  !  the dycoms II RF01 case
        tmp_zl(1:nz-2) = ( hoc_zl(3:nz) - les_zl(1:nz-2) ) / norm_term 
        tmp_zl = tmp_zl**2
        mean_sqr_diff_zt = sum( tmp_zl(1:nz-2), 1 )

      case (50)   !  the wangara case
        tmp_zl = ( hoc_zl(1:nz) - les_zl(1:nz) ) / norm_term 
        tmp_zl = tmp_zl**2
        mean_sqr_diff_zt = sum( tmp_zl, 1 )

      case default !
        stop "Not able to handle specified number of HOC z-levels"
      end select

      return
      end function mean_sqr_diff_zt
!-----------------------------------------------------------------------
!  FUNCTION mean_sqr_diff_zm
!  Calculate the mean squared difference between two input vectors, then
!  normalize.
!  Configured to do interpolation on LES / HOC comparisons on the m grid
!-----------------------------------------------------------------------
      real function mean_sqr_diff_zm( hoc_zl, les_zl, norm_term )
      
      implicit none

!     Input
      real, intent(in), dimension(:) :: hoc_zl 
      real, intent(in), dimension(:) :: les_zl 

      real, intent(in) :: norm_term ! normalization term; 
                                    ! typically minval(les) - maxval(les) 
!     Internal
      real, dimension(size( hoc_zl )) :: tmp_zl
      integer nz 

      nz = size( hoc_zl ) 

      select case ( nz )
      case (48, 50, 75, 110, 147, 150, 160) ! most cases
        tmp_zl(1:nz) = ( hoc_zl(1:nz) - les_zl(1:nz) )/norm_term 
        tmp_zl       = tmp_zl**2
        mean_sqr_diff_zm  = sum( tmp_zl(1:nz), 1 ) 

      case (132)  !  the dycoms II RF01 case
        tmp_zl(1:nz-2) = ( hoc_zl(1:nz-2) - les_zl(1:nz-2) )/norm_term 
        tmp_zl = tmp_zl**2
        mean_sqr_diff_zm = sum( tmp_zl(1:nz-2), 1 )

      case default !
        stop "Not able to handle specified number of HOC z-levels"
      end select

      return
      end function mean_sqr_diff_zm
!-----------------------------------------------------------------------
!  SUBROUTINE read_random_seed( )
!  Reads a flat file passed as an argument
!-----------------------------------------------------------------------
      subroutine read_random_seed( seed_file )

      implicit none
!     Input
      character(len=*), intent(in) :: seed_file
 
!     Internal
      integer, dimension(:), allocatable :: rand_seed

      integer rand_size      ! for random seed
      integer AllocateStatus
      integer InputStatus
      integer i

      call random_seed( SIZE = rand_size )

      allocate(rand_seed( rand_size ), stat = AllocateStatus)
      if ( AllocateStatus /= 0 ) then 
        stop "Allocation of the random seed variable array failed"
      endif

!     ASCII formatted file, usually generated by int2txt
      open(unit=30, file=seed_file, action='read')
      do i=1, rand_size 
        read(unit=30, fmt=*, iostat=InputStatus) rand_seed(i)
        if ( InputStatus /= 0 ) then 
          print *, "error reading"//seed_file
          stop
        endif
      enddo
      close(30)
      
      call random_seed( PUT = rand_seed )

      deallocate(rand_seed)

      return
      end subroutine read_random_seed

      end module error
!-----------------------------------------------------------------------
