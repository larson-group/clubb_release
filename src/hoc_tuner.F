!----------------------------------------------------------------------
! $Id: hoc_tuner.F,v 1.27 2008-06-20 13:55:35 faschinj Exp $

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

!           Portability module from _Numerical Recipes In Fortran 90_ 
!----------------------------------------------------------------------
      use error, only: 
     .  tuner_init, min_les_hoc_diff,                 ! Subroutines 
     .  output_results_stdout,                        ! Subroutine
     .  output_nml_standalone, output_nml_tuner,      ! Subroutines
     .  param_vals_matrix,           ! Variables
     .  results_stdout,                                ! Variables
     .  results_file, tune_type, ftol, ndim  ! Variables

      use constants, only:
     . fstdout ! Variables

      implicit none

      ! Variables
      character(len=10) :: current_time  ! Current time string (no seconds)
      character(len=8)  :: current_date  ! Current date string
      character(len=50) :: results_f     ! Results file

      character(len=1)  :: user_response ! Simple Y/N query

      ! Read in namelists and define parameters
      call tuner_init( read_files=.true. )

      ! Attempt to find the optimal parameter set
      do
        if ( tune_type == 0 ) then 
          call amoeba_driver( )
        else
          call amebsa_driver( )
        end if

        ! Print to stdout if specified
        if ( results_stdout ) call output_results_stdout( )

        ! Query to see if we should exit the loop
        write(fstdout,*) "Run Complete."
        write(unit=fstdout,fmt='(A)', advance='no') 
     .    "Re-run with new parameters?(y/n) "
        read(*,*) user_response

        if ( trim( user_response ) /= "y" .and.
     .       trim( user_response ) /= "Y"   ) then
          exit 
        end if
       
        write(fstdout,*) "Current ftol= ", ftol 
        write(fstdout,fmt='(A)', advance='no') "Enter new ftol=   "
        read(*,*) ftol

        call tuner_init( read_files=.false. ) 

      end do ! user_response /= 'y', 'Y' or 'yes'

      ! Final namelist file output 

      if ( results_file ) then 

        ! Tuner namelist
        print *, "Generating new error.in file..."

        call date_and_time( current_date, current_time )

        results_f = "error_" // current_date // '_'
     .              // current_time(1:4) // ".in" 

        ! Note:
        ! The first column of param_vals_matrix is the optimized result, which
        ! is swapped in by amoeba.
        call output_nml_tuner( results_f, param_vals_matrix(1,1:ndim) )
        print *, "New filename is: ", results_f

        ! Standalone namelist
        print *, "Generating new standalone.in file..."

        results_f = "../standalone/standalone_" // current_date // '_' 
     .              // current_time(1:4) // ".in" 

        call output_nml_standalone( results_f, 
     .                              param_vals_matrix(1,1:ndim) )
        print *, "New filename is: ", results_f

      end if

      ! Exit Program

      stop "Program exited normally"

      end program hoc_tuner

!------------------------------------------------------------------------
      subroutine amoeba_driver

!     Description:
!     Simple interface for the amoeba minimization algorithm

!     References:
!     _Numerical Recipes in Fortran 90_.  See full citation above.
!------------------------------------------------------------------------
      use nr, only: 
     .    amoeba ! Procedure(s)
      use error, only: 
        ! Variable(s)
     .  ndim,                               ! Array dimensions
     .  param_vals_matrix, cost_fnc_vector, ! The 'p' matrix and 'y' vector resp.
     .  ftol,                               ! Tolerance of tuning run
     .  iter,                               ! Iteration number
     .  min_les_hoc_diff,                   ! Cost function
     .  min_err                             ! Minimum value of the cost function


      implicit none

      ! External
      intrinsic :: minval

      call amoeba( param_vals_matrix(1:ndim+1,1:ndim), 
     .             cost_fnc_vector(1:ndim+1), 
     .             ftol, min_les_hoc_diff, iter)

      ! Note:
      ! Amoeba will make the optimal cost result the first element of
      ! cost_fnc_vector and the optimal parameter set the first column of 
      ! param_vals_matrix, where param_vals_matrix is 'p' in the NR subroutine
      ! and cost_fnc_vector is 'y'.

      min_err = cost_fnc_vector(1)

      return
      end subroutine amoeba_driver

!-----------------------------------------------------------------------

      subroutine amebsa_driver

!     Description:
!     Interface for the amoeba simulated annealing minimization algorithm
!     At the end of the subroutine, the param_vals_matrix's first row gets the
!     optimal values assigned to it.
!-----------------------------------------------------------------------
      use nr, only: 
     .    amebsa ! Procedure(s)
      use nrtype, only: 
     .    SP ! Variable(s)
      use error, only:
     .    param_vals_matrix, ! Variable(s)
     .    anneal_temp,
     .    anneal_iter,
     .    iter,
     .    ndim,
     .    cost_fnc_vector,
     .    ftol,
     .    min_err,
     .    min_les_hoc_diff ! Procedure(s)

      implicit none

      ! Note:
      ! Most of these are taken from xamebsa, so I only have a vague idea
      ! what the their purpose is and can't pick less obscure names.
      ! dschanen 1 Apr 05
 
      ! Local Variables
      integer :: 
     .  iiter, jiter,! Loop variables
     .  nit ! ???

      real, dimension(ndim) :: 
     .  pb ! ???
 
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
     .       ( param_vals_matrix(1:ndim+1,1:ndim), 
     .         cost_fnc_vector(1:ndim+1),
     .         pb(1:ndim), yb, ftol, min_les_hoc_diff, iter, tmptr )

        nit = nit + iiter - iter
        if ( yb < ybb ) then
          ybb = yb
        end if
        if ( iter > 0 ) exit
      enddo
      
      param_vals_matrix(1,1:ndim) = pb(1:ndim)
      min_err = ybb

      end subroutine amebsa_driver
!----------------------------------------------------------------------
