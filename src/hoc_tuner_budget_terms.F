!-----------------------------------------------------------------------
!  $Id: hoc_tuner_budget_terms.F,v 1.7 2005-08-23 21:11:10 dschanen Exp $

! PROGRAM hoc_tuner_budget_terms

! Attempts to tune the constants of the HOC model to best fit the
! individual budget terms of prognostic variables.

! Currently setup for wp2, wp3, wprtp, and wpthlp;
! Can tune buoyancy, pressure, dissipation & advection budgets.
!-----------------------------------------------------------------------

      program hoc_tuner_budget_terms
      use budget_terms
      use error, only: read_random_seed

      implicit none

!     Parameters
      logical, parameter :: c_update = .true.
      logical, parameter :: liter    = .false. ! attempt to optimize for all
                                               ! cases.  Takes some time and
                                               ! may never converge.

!     Internal
      real ftol
      integer i, iterations

      character(10) current_time      ! Current time string
      character(8) current_date       ! Current date string

      character(1) user_response
       

!     Setup the budget terms module
      call setup_budget_terms( "budget.in", ftol)
      
!     Seed the random number generator
      call read_random_seed("../tune/rand_seed.dat") ! from error.mod

      if ( liter ) then
        iterations = 10
      else
        iterations = 1      
      endif

      do i=iterations, 1, -1

        call execute_tuner_runs( )

        if ( i < iterations ) cycle

        write(*,*) "Run Complete."
        write(*,'(A)', advance='no') 
     .       "Re-run with new parameters?(y/n) "
        read *, user_response

        if ( trim(user_response) /= "y" .and.
     .       trim(user_response) /= "Y"   ) then 
          exit 
        endif
       
        write (*,*) "Current ftol= ", ftol 
        write (*, '(A)', advance='no') "Enter new ftol= "
        read *, ftol
      enddo

      call date_and_time (current_date, current_time)
      call generate_namelists( "budget_"//current_date//
     .                         "_"//current_time(1:4)//".in", ftol ) 
      call output_optimal_constants( )

      stop "Program Exited Normally"

      contains
!-----------------------------------------------------------------------
!  SUBROUTINE execute_tuner_runs( )
!-----------------------------------------------------------------------
      subroutine execute_tuner_runs

      implicit none

      if ( lwp2_prdp ) then
        write(*,*) "running wp2(pr+dp)"
        call wpn_run( ftol*i, "2", wp2_min, c_update )
      endif

      if ( lwp3_prdp ) then
        write(*,*) "running wp3(pr+dp)"
        call wpn_run( ftol*i, "3", wp3_min, c_update )
      endif

      if ( lwprtp_prdp ) then
        write(*,*) "running wprtp(pr+dp)"
        call wpxp_run( ftol*i, "rt", wprtp_min, c_update )
      endif

      if ( lwpthlp_prdp ) then
        write(*,*) "running wpthlp(pr+dp)"
        call wpxp_run( ftol*i, "thl", wpthlp_min, c_update )
      endif

      if ( lwp2_bp ) then
        write(*,*) "running wp2_bp"
        call buoy_advect_run( ftol*i, "", wp2_buoy_min, c_update )
      endif

      if ( lwp3_bp ) then
        write(*,*) "running wp3_bp"
        call buoy_advect_run( ftol*i, "", wp3_buoy_min, c_update )
      endif

      if ( lwpthlp_bp ) then
        write(*,*) "running wpthlp_bp"
        call buoy_advect_run( ftol*i, "", wpthlp_buoy_min, c_update )
      endif

      if ( lwprtp_bp ) then
        write(*,*) "running wprtp_bp"
        call buoy_advect_run( ftol*i, "", wprtp_buoy_min, c_update )
      endif

      if ( lwp3_tp ) then
        write(*,*) "running wp3_tp"
        call buoy_advect_run( ftol*i, "", wp3_advect_min, c_update )
      endif

      if ( lwpthlp_tp ) then
        write(*,*) "running wpthlp_tp"
        call buoy_advect_run( ftol*i, "", wpthlp_advect_min, c_update )
      endif

      if ( lwprtp_tp ) then
        write(*,*) "running wprtp_tp"
        call buoy_advect_run( ftol*i, "", wprtp_advect_min, c_update )
      endif

      return
      end subroutine execute_tuner_runs
!-----------------------------------------------------------------------

      end program hoc_tuner_budget_terms
!-----------------------------------------------------------------------
