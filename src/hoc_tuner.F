!----------------------------------------------------------------------
! $Id: hoc_tuner.F,v 1.17 2006-12-16 22:17:54 dschanen Exp $

      program hoc_tuner 

!     Description:
!     ``Tunes'' constants in hoc so that the output matches LES output.
!     Uses amoeba or amebsa to calculate the min of (f_les - f_hoc)^2

!     References:
!     _Numerical Recipes in Fortran 90_ (Chapter 10) 
!     (Amoeba & Amebsa subroutine)

!     Subroutines Called :
!      error_init      : Reads in
!                         * GRADS parameters: filenames, z-levels, timesteps, 
!                          variable to tune for
!                         * HOC variables: C1,...C11, nu1,...nu8
!                         * Amoeba tolerance: ftol
!                        Initializes
!                         * The initial dependent variable vector, i.e. the 
!                          total error between the les and hoc models
!                         * The initial independent variable array, formed
!                           from C1,...C11, and nu1,...nu8

!     AMOEBA    : The downhill simplex tuner from _Numerical Recipes_
!     AMEBSA    : The simulated annealing tuner from _Numerical Recipes_
!
!   Functions Called :
!      min_les_hoc_diff  : A parameter of the amoeba subroutine, and 
!                        strictly speaking is directly called from
!                        the hoc_tuner program only initially.  
!                        This is the counterpart to the
!                        minimization function `funk' in amoeba. It returns
!                        the error vector and takes the independent variable
!                        vector as a parameter.

!   Included Common:
!       error.mod: subroutine tuner_init, function min_les_hoc_diff, logical isValid,
!                   integer ndim, integer mdim, real x_array, real err_array 
!                   real init_error, results_stdout, results_file

!       nr.mod 
!         Contain interface's for: 
!         amoeba.o : subroutine amoeba
!         amebsa.o : subroutine amebsa

!       nrtype.mod
!           portability module from _Numerical Recipes In Fortran 90_ 
!----------------------------------------------------------------------
      use error, only: 
     .  tuner_init, min_les_hoc_diff,                 ! Subroutines 
     .  output_results_stdout,                        ! Subroutine
     .  output_nml_standalone, output_nml_tuner,      ! Subroutines
     .  x_array, err_array, init_err, results_stdout, ! Variables
     .  results_file, isValid, tune_type, ftol, ndim, mdim ! Variables

      implicit none

      ! Variables
      character(len=10) :: current_time  ! Current time string (no seconds)
      character(len=8)  :: current_date  ! Current date string
      character(len=50) :: results_f     ! Results file

      character(len=1)  :: user_response ! Simple Y/N query

!-------------------------- Define Constants ---------------------------- 

      call tuner_init(read_files=.true.)

!-------------- Calculate the error for the initial values --------------
!     Notes: min_les_hoc_diff 
!     is written to deal with undefined and invalid values for
!            variations on the initial vector, but that algorithm relies 
!            on the initial vector being valid. 
!------------------------------------------------------------------------

      init_err =  min_les_hoc_diff( x_array(1, 1:ndim) )
      if ( .not. isValid ) then
         print *, "Initial variable values must be valid"
         stop
      end if

!------------------ Attempt to find the optimal values -----------------
      do
        if ( tune_type == 0 ) then 
          call amoeba_driver( )
        else
          call amebsa_driver( )
        end if

!---------------- Print to stdout if specified --------------------------
        if ( results_stdout ) call output_results_stdout

!---------------- Query to see if we should end the loop  ---------------
        write(*,*) "Run Complete."
        write(unit=*,fmt='(A)', advance='no') 
     .        "Re-run with new parameters?(y/n) "
        read *, user_response

        if ( trim(user_response) /= "y" .and.
     .       trim(user_response) /= "Y"   ) then 
          exit 
        end if
       
        write (*,*) "Current ftol= ", ftol 
        write (*, '(A)', advance='no') "Enter new ftol=   "
        read *, ftol

        call tuner_init(read_files=.false.) 

      end do

!     This is done through a script now

!     call output_optimal_results( )

!--------------- Final namelist file output -----------------------------
      if ( results_file ) then 
        ! Tuner namelist
        print *, "Generating new error.in file..."

        call date_and_time( current_date, current_time )

        results_f = "error_" // current_date // '_' 
     .              // current_time(1:4) // ".in" 

        call output_nml_tuner( results_f, x_array(1,1:ndim) )
        print *, "New filename is: ", results_f

        ! Standalone namelist
        print *, "Generating new standalone.in file..."

        results_f = "../standalone/standalone_" // current_date // '_' 
     .              // current_time(1:4) // ".in" 

        call output_nml_standalone( results_f, x_array(1,1:ndim) )
        print *, "New filename is: ", results_f

      end if

!-------------------------- Exit Program  -------------------------------

      stop "Program exited normally"

      end program hoc_tuner

!------------------------------------------------------------------------
      subroutine amoeba_driver

!     Description
!     Simple interface for the amoeba minimization algorithm
!------------------------------------------------------------------------
      use nr, only: amoeba
      use error

      implicit none

      ! External
      intrinsic :: minval

      call amoeba( x_array(1:mdim,1:ndim), err_array(1:mdim), 
     .             ftol, min_les_hoc_diff, iter)

      min_err =  minval( err_array )
        
      return
      end subroutine amoeba_driver

!-----------------------------------------------------------------------

      subroutine amebsa_driver

!     Description:
!     Interface for the amoeba simulated annealing minimization algorithm
!     At the end of the subroutine, the x_array's first row gets the
!     optimal values assigned to it.
!-----------------------------------------------------------------------
      use nr, only: amebsa
      use nrtype, only: SP
      use error

      implicit none

      ! Note:
      ! Most of these are taken from xamebsa, so I only have a vague idea
      ! what the their purpose is and can't pick less obscure names.
      ! dschanen 1 Apr 05
 
      ! Local Variables
      integer :: 
     .  iiter, jiter,! Loop variables
     .  nit ! ???
   
      real :: 
     .  ybb,  ! ???
     .  yb,   ! ???
     .  tmptr ! ???

      ybb   = 1.0e30_SP
      yb    = 1.0e30_SP
      nit   = 0
      iiter = 0
      tmptr = anneal_temp ! anneal_temp taken from /stat/ namelist

      do jiter = 1, anneal_iter ! anneal_iter taken from /stat/ namelist
        iter  = iiter
        tmptr = tmptr * 0.8_SP

        call amebsa
     .       ( x_array(1:mdim,1:ndim), err_array(1:mdim),
     .         pb(1:ndim), yb, ftol, min_les_hoc_diff, iter, tmptr )

        nit = nit + iiter - iter
        if ( yb < ybb ) then
          ybb = yb
        end if
        if ( iter > 0 ) exit
      enddo
      
      x_array(1,1:ndim) = pb(1:ndim)
      min_err = ybb

      end subroutine amebsa_driver
!----------------------------------------------------------------------
