!----------------------------------------------------------------------
!   HOC_TUNER.F
!   ``Tunes" constants in hoc so that the output matches LES output.
!   Use amebsa.f90 to calculate the min of (f_les - f_hoc)^2

!   References:
!   Numerical Recipes in Fortran 90 (Chapter 10) (Amoeba & Amebsa subroutine)

!   Subroutines Called :
!      ERROR_INIT      : Reads in
!                         * GRADS parameters: filenames, z-levels, timesteps, 
!                          variable to tune for
!                         * HOC variables: C1,...C11, nu1,...nu8
!                         * Amoeba tolerance: ftol
!                        Initializes
!                         * The initial dependent variable vector, i.e. the 
!                          total error between the les and hoc models
!                         * The initial independent variable array, formed
!                           from C1,...C11, and nu1,...nu8

!     AMOEBA		: The downhill simplex tuner from Numerical Recipes
!     AMEBSA		: The simulated annealing tuner from Numerical Recipes
!
!   Functions Called :
!      FUNK		:  A parameter of the AMOEBA subroutine, and 
!                        strictly speaking is directly called from
!                        the hoc_tuner program only initially.  
!                        This is the counterpart to the
!                        minimization function FUNK in amoeba. It returns
!                        the error vector and takes the independent variable
!                        vector as a parameter.

!   Included Common:
!       error.mod: subroutine error_init, function funk, logical isValid,
!                   integer ndim, integer mdim, real x_array, real err_array 
!                   real init_error, real init_error, 
!                   results_stdout, results_file

!       nr.mod 
!         Contain interface's for: 
!         amoeba.o : subroutine amoeba
!         amebsa.o : subroutine amebsa

!       nrtype.mod
!           portability module from _Numerical Recipes In Fortran 90_ 
!----------------------------------------------------------------------
      program hoc_tuner 
      use error

      implicit none

      character(4)  :: current_time       ! Current time string (no seconds)
      character(8)  :: current_date       ! Current date string
      character(50) :: results_f          ! Results file

      character(1)    :: user_response   ! simple Y/N query

!-------------------------- Define Constants ---------------------------- 

      call tuner_init(read_files=.true.)

!-------------- Calculate the error for the initial values --------------
!     Notes: funk is written to deal with undefined and invalid values for
!            variations on the initial vector, but that algorithm relies 
!            on the initial vector being valid. 
!------------------------------------------------------------------------

      init_err =  funk( x_array(1, 1:ndim) )
      if ( .not. isValid ) then
         print *, "Initial variable values must be valid"
         stop
      endif

!------------------ Attempt to find the optimal values -----------------
      do
        if ( tune_type == 0 ) then 
          call amoeba_driver
        else
          call amebsa_driver
        endif

!---------------- Print to stdout if specified --------------------------
        if ( results_stdout ) call output_results_stdout

!---------------- Query to see if we should end the loop  ---------------
        write(*,*) "Run Complete."
        write(*,'(A)', advance='no') 
     .        "Re-run with new parameters?(y/n) "
        read *, user_response

        if ( trim(user_response) /= "y" .and.
     .       trim(user_response) /= "Y"   ) then 
          exit 
        then
       
        write (*,*) "Current ftol= ", ftol 
        write (*, '(A)', advance='no') "Enter new ftol=   "
        read *, ftol

        call tuner_init(read_files=.false.) 

      enddo

!--------------- Final namelist file output -----------------------------
      if ( results_file ) then 
!----------------------- Tuner namelist ---------------------------------      
        print *, "Generating new error.in file..."

        call date_and_time (current_date, current_time)

        results_f = "error_" // current_date // '_' 
     .              // current_time // ".in" 

        call output_nml_tuner( results_f, x_array(1,1:ndim) )
        print *, "New filename is: ", results_f

!----------------------- Standalone namelist ----------------------------      
        print *, "Generating new standalone.in file..."

        results_f = "../standalone/standalone_" // current_date // '_' 
     .              // current_time // ".in" 

        call output_nml_standalone( results_f, x_array(1,1:ndim) )
        print *, "New filename is: ", results_f

      endif

!-------------------------- Exit Program  -------------------------------

      stop "program exited normally"

      end program hoc_tuner

!------------------------------------------------------------------------
!   Simple interface for the amoeba minimization algorithm
!------------------------------------------------------------------------
      subroutine amoeba_driver
      use nr, only: amoeba
      use nrtype
      use error

      implicit none
     
      call amoeba( x_array(1:mdim, 1:ndim), err_array(1:mdim), 
     .             ftol, funk, iter)

      min_err = minval( err_array )

      end subroutine amoeba_driver

!------------------------------------------------------------------------
!  Interface for the amoeba simulated annealing minimization algorithm
!  At the end of the subroutine, the x_array's first row gets the optimal
!  values assigned to it.
!------------------------------------------------------------------------

      subroutine amebsa_driver
      use nr, only : amebsa
      use nrtype
      use error

      implicit none

!----------------------------------------------------------------------
!  Internal Variables
!    Most of these are taken from xamebsa, so I only have a vague idea
!    what the their purpose is and can't pick less obscure names.
!    dschanen 1 Apr 05
!----------------------------------------------------------------------

      integer :: jiter ! loop variable
      integer :: iiter ! loop variable
      integer :: nit
   
      real    :: ybb
      real    :: yb
      real    :: tmptr

      ybb   = 1.0e30_sp
      yb    = 1.0e30_sp
      nit   = 0
      iiter = 0
      tmptr = anneal_temp ! anneal_temp taken from /stat/ namelist

      do jiter = 1, anneal_iter ! anneal_iter taken from /stat/ namelist
        iter  = iiter
        tmptr = tmptr * 0.8_sp
        call amebsa( x_array(1:mdim,1:ndim), err_array(1:mdim),
     .               pb(1:ndim), yb, ftol, funk, iter, tmptr )
        nit = nit + iiter - iter
        if (yb < ybb) then
          ybb = yb
        endif
        if (iter > 0) exit
      enddo
      
      x_array(1,1:ndim) = pb(1:ndim)
      min_err = ybb

      end subroutine amebsa_driver
!----------------------------------------------------------------------
