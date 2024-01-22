!-------------------------------------------------------------------------------
! $Id$
module enhanced_simann
! Description:
!   Implementation of Siarry's Enhanced simulated annealing algorthm in
!   Fortran 90/95.

! References: 
!   ``Enhanced Simulated Annealing for Many Globally Minimized Functions
!   of Many-Continuous Variables'', Siarry, et al. ACMS TOMS Vol. 23,
!   No. 2, June 1997, pp. 209--228. https://dl.acm.org/citation.cfm?id=264043
!
!------------------------------------------------------------------------------

    implicit none

    public :: esa_driver, esa_driver_siarry

    private :: exec_movement, exec_movement_siarry, select_partition, init_random

    private ! Default scope

    
    logical, parameter :: &
        l_esa_debug_statements = .false. ! Verbose debugging output


    contains
!-----------------------------------------------------------------------------
    subroutine esa_driver( xinit, xmin, xmax,                           & ! intent(in)
                           initial_temp, max_final_temp,                & ! intent(inout)
                           xopt, nrgy_opt,                              & ! intent(out)
                           fobj,                                        & ! procedure
                           stp_adjst_shift_in, stp_adjst_factor_in,     & ! optional(in)
                           max_iters_in,                                & ! ^
                           f_tol_in, file_name, file_unit,              & ! ^
                           iter_in                                      ) ! optional(inout)

    ! Description:
    !   Driver subroutine
    ! 
    ! References:
    !   None
    !-----------------------------------------------------------------------------

        use constants_clubb, only: eps

        use clubb_precision, only: &
          core_rknd ! Variable(s)

        implicit none

        ! Interface for fobj
#include "enhanced_simann_fobj.inc"   

        real( kind = core_rknd ), dimension(:), intent(in) :: &
          xmin,     & ! Minimum values for the x vector
          xmax,     & ! Maximum values for the x vector
          xinit       ! Initial argument for fobj

        real( kind = core_rknd ), intent(in) :: &
          initial_temp, &    ! Initial temperature
          max_final_temp    ! Maximum allowed final temperature

        real( kind = core_rknd ), dimension(:), intent(out) :: &
          xopt        ! Array of optimal parameters

        real( kind = core_rknd ), intent(out) :: &
          nrgy_opt    ! Optimal value of the cost function

        ! Input variables
        real( kind = core_rknd ), intent(in), optional :: &
            stp_adjst_shift_in,     & ! If improvement_ratio = this, no step adjust
            stp_adjst_factor_in,    & ! Slope around stp_adjst_shift
            f_tol_in                  ! Threshold value to determine if cost improved

        character(len=50), intent(in), optional :: &
            file_name

        integer, intent(in), optional :: &
            file_unit, &
            max_iters_in        ! Input value for max iterations

        integer, intent(inout), optional :: &
            iter_in

        ! Local variables

        ! Stop conditions parameters
        integer, parameter :: &
          no_improve_max = 3      ! Max num of temperature stages we allow 
                                  ! since the most recent improvement.
                                  ! A stage is a series of iterations at the same temperature.
                                  ! Once stages_w_no_improve reaches no_improve_max, 
                                  ! then the optimization stops.

        ! Parameters for temperature adjustment
        real( kind = core_rknd ), parameter :: &
          use_max  = 0.8_core_rknd,     & ! The maximum chance for a variable to be used
          use_min  = 0.2_core_rknd,     & ! The minimum chance for a variable to be used
          min_cool = 0.7_core_rknd,     & ! We insist that temperature at new stage cools by at 
                                          ! least min_cool times the temperature at the prior stage,
          max_cool = 0.3_core_rknd        ! but no more than max_cool times the prior temperature.

        ! Parameters for step vector adjustment
        real( kind = core_rknd ) :: &
          stp_adjst_shift  = 0.5_core_rknd, & ! If improvement_ratio = this, no step change
          stp_adjst_factor = 1.0_core_rknd, & ! Slope around stp_adjst_shift
          f_tol            = 1.e-4            ! Threshold value to determine if cost improved

        ! Variable names based not on the paper because they were stupid
        real( kind = core_rknd ), dimension(size( xinit )) :: &
          xstart,               & ! Starting point for x 
          xtry,                 & ! Test point for x 
          stp_ini,              & ! Initial step vector
          stp_cur,              & ! Current step vector
          improvement_ratio,    & ! Number of times a variable helped vs tried, improved/attempted
          improved,             & ! Number of times a variable improved results
          attempted               ! Number of times a variable changed

        real( kind = core_rknd ) :: &
          einit,                & ! Initial value of the cost function
          old_nrgy,             & ! Old energy level
          new_nrgy,             & ! New energy level
          delta_nrgy,           & ! Change in energy, i.e. delta fobj( xtry )
          temp                    ! Anneal temperature

        logical, dimension(size( xinit )) :: &
          xpartition              ! List of variables being used this iteration

        integer :: &
          max_iters = 2000,     & ! Max iterations
          iter,                 & ! Total iterations completed
          stages_w_no_improve,  & ! Number of temperature stages since the most recent improvement
          vars,                 & ! Size of the a partion
          n1, n2,               & ! Annealing schedule
          k                       ! Loop variable

        real( kind = core_rknd ) :: &
          min_nrgy,             & ! Minimal fobj value at current temperature stage
          tot_nrgy                ! Sum of successive fobj values at current temperature stage

        real( kind(0.0) ) :: &
          rand                    ! random number from [0,1]

        ! ---- Begin Code ----

        ! initialize random number generator and initial values
        call init_random
        tot_nrgy = 0._core_rknd
        stages_w_no_improve = 0
        vars = size( xinit )
        xopt(:) = xinit(:)
        xstart = xinit     
        temp = initial_temp
        iter = 1
        if ( present(iter_in) ) iter_in = 1

        ! set values if optional arguements if present
        if ( present(f_tol_in) ) f_tol = f_tol_in
        if ( present(stp_adjst_shift_in) ) stp_adjst_shift = stp_adjst_shift_in
        if ( present(stp_adjst_factor_in) ) stp_adjst_factor = stp_adjst_factor_in
        if ( present(max_iters_in) ) max_iters = max_iters_in
        
        ! set annealing schedule, suggested values from the paper
        n1 = 12
        n2 = 50

        ! initialize all energy variables with initial energy value
        einit = real(fobj( real(xinit) ), kind = core_rknd)
        nrgy_opt      = einit     
        old_nrgy      = einit
        min_nrgy      = einit 

        ! set initial step value, 0.35 chosen experimentally
        stp_ini(:) = ( xmax(:) - xmin(:) ) * 0.35_core_rknd
        stp_cur    = stp_ini

        ! intialize improvement and attempts trackers
        improved(:)           = 0._core_rknd  ! Number of improved moves for each var
        attempted(:)          = 0._core_rknd  ! Number of attempted moves for each var
        improvement_ratio(:)  = 1._core_rknd  ! Ratio of improved vs attempted

        ! debugging
        if ( present(file_name) .and. present(file_unit) ) then
            open(unit=file_unit, file=file_name, action='write', position='append')
            write(file_unit,*) "---------------INTITIAL STAGE---------------"
            write(file_unit,*) "Temp = ", temp
            write(file_unit,*) "Step = ", stp_cur
            write(file_unit,*) "Xinit = ", xinit
        end if

        do 
   
            ! randomly choose what variables are used, chance based on how useful variable was
            do k = 1, vars

                call random_number( rand )
                if ( rand < use_min + ( use_max - use_min ) * improvement_ratio(k) ) then
                    xpartition(k) = .true.
                else
                    xpartition(k) = .false.
                end if

            end do

            ! use every variable for the first iteration
            where( abs(attempted) < eps) xpartition = .true.

            ! if no variables were selected to be used, select all
            if ( .not. any( xpartition(:) ) ) xpartition(:) = .true.

            ! randomly change x values a little bit, starting at last accepted values
            ! then call fobj with new x values
            xtry = xstart
            call exec_movement( xpartition, xmax, xmin, stp_cur, fobj, &
                                xtry, new_nrgy )

            ! increment attempt counter for used variables
            where ( xpartition ) attempted = attempted + 1

            ! calculate difference between old and new energy
            delta_nrgy = new_nrgy - old_nrgy                

            ! add new energy to total
            tot_nrgy = tot_nrgy + new_nrgy                                     

            ! increment iterations
            iter = iter + 1
            if ( present(iter_in) ) iter_in = iter

            ! if energy decreased
            if ( delta_nrgy <= 0.0_core_rknd ) then

                ! accept values and increment improvement counter for used variables
                where(xpartition) improved = improved + 1

                ! save energy value and x values
                old_nrgy = new_nrgy
                xstart = xtry

                ! if energy is best for this temperature stage, save energy value only
                if ( new_nrgy < min_nrgy ) then
                  min_nrgy = new_nrgy
                end if

                ! If new energy is the best (smallest) among all temperature stages, 
                ! then save it and the variables
                if ( new_nrgy < nrgy_opt ) then

                    ! If the new best is significantly better than the old best, 
                    ! then reset the num of stages since last improvement to zero.
                    ! This will cause us to iterate more at the same temperature,
                    ! instead of coming closer to reaching the stopping criterion.
                    if ( nrgy_opt - new_nrgy > f_tol * nrgy_opt ) then
                        stages_w_no_improve = 0
                    end if

                    xopt = xtry
                    nrgy_opt = new_nrgy

                end if

            else    
               ! energy increased (answer worsened), but randomly accept it with temp based chance
           
               call random_number( rand )                               
               if ( rand <= exp( -delta_nrgy / temp ) ) then

                    old_nrgy = new_nrgy
                    xstart = xtry

               end if

            end if


            ! If there has not been much improvement, and not many attempts, 
            ! keep going at this step size and temperature.
            ! If there has been either a lot of improvement, or too many attempts, 
            ! move onto the next temperature stage
            if ( sum(improved) <  n1*vars .and. sum(attempted) < n2*vars ) then
                cycle                   
            end if

            ! temp *= decay, where max_cool < decay = min at current temp / average < min_cool
            if ( tot_nrgy > 0.0_core_rknd ) then
                temp = temp * max( min( min_nrgy * sum(attempted)/tot_nrgy, min_cool ), max_cool )
            else 
                temp = temp * min_cool
            end if

            ! improvement ratio = how well a variable did in improving results
            improvement_ratio(:) = improved(:) / attempted(:)

            ! adjust variable based on improvement ratio
            ! if improvement_ratio == stp_adjst_shift, then step size will
            ! remain the same. 
            stp_cur(:) = stp_cur(:) * (  stp_adjst_factor* &
                         ( improvement_ratio(:) - stp_adjst_shift ) + 1 )

            ! too many iterations, so let's stop.
            if ( iter >= max_iters ) then

                ! debugging
                if ( present(file_name) .and. present(file_unit) ) then
                    open(unit=file_unit, file=file_name, action='write', position='append')
                    write(file_unit,*) "Iteration Stop"
                end if

                exit
            end if

            ! too many temperature stages with no improvement, so let's stop.
            if ( temp <= max_final_temp .and. stages_w_no_improve > no_improve_max ) then

                ! debugging
                if ( present(file_name) .and. present(file_unit) ) then
                    open(unit=file_unit, file=file_name, action='write', position='append')
                    write(file_unit,*) "Improvement Stop"
                end if

                exit
            end if

                
            ! debugging 
            if ( present(file_name) .and. present(file_unit) ) then
                open(unit=file_unit, file=file_name, action='write', position='append')
                write(file_unit,*) "---------------NEW TEMP STAGE---------------"
                write(file_unit,*) "Iteration: ", iter
                write(file_unit,*) "Stages with no improvment: ", stages_w_no_improve
                write(file_unit,*) "Improved: ", improvement_ratio
                write(file_unit,*) "Temp = ", temp
                write(file_unit,*) "Step = ", stp_cur
                write(file_unit,*) "Xopt = ", xopt
                write(file_unit,*) "cost = ", nrgy_opt
            end if
                
            ! end of stage, so increment no improvement counter
            stages_w_no_improve = stages_w_no_improve + 1

            ! reset counters for new temp stage
            improved    = 0._core_rknd
            attempted   = 0._core_rknd
            tot_nrgy    = 0._core_rknd
                
            ! set starting x back to best results found
            xstart = xopt
            min_nrgy = nrgy_opt

        end do

        return
    end subroutine esa_driver

    !-----------------------------------------------------------------------------
    subroutine exec_movement( xpartition, xmax, xmin, step, fobj, &
                                       x, cost )
    ! Description:
    !   Execute one movement in domain space

    ! References:
    !   Step 3 in Siarry, et al.
    !-----------------------------------------------------------------------------

        use clubb_precision, only: &
          core_rknd ! Variable(s)

        implicit none

        ! External
        intrinsic :: size

#include "enhanced_simann_fobj.inc"

        ! Input Variables
        logical, dimension(:), intent(in) :: xpartition

        real( kind = core_rknd ), dimension(:), intent(in) :: xmax, xmin, step

        ! Input/Output Variables
        real( kind = core_rknd ), dimension(:), intent(inout) :: x

        ! Output Variables
        real( kind = core_rknd ), intent(out) :: cost

        ! Local variables
        real( kind(0.0) ), dimension(3) :: xrand

        real( kind = core_rknd ) :: xtmp

        integer :: k

        ! --- Begin Code ---

        do k = 1, size( x )
            if ( xpartition(k) ) then ! Augment only values within the partition

                ! Apply a random sign to each xrand value. 
                call random_number( xrand )

                if ( xrand(1) >= 0.5 ) then
                  xtmp = x(k) + real( xrand(2), kind = core_rknd ) * step(k)
                else
                  xtmp = x(k) - real( xrand(3), kind = core_rknd ) * step(k)
                end if

                if ( xtmp > xmax(k) ) then
                    xtmp = xmax(k)
                else if ( xtmp < xmin(k) ) then
                    xtmp = xmin(k)
                end if

                x(k) = xtmp

            end if  ! Points in the partition

        end do ! 1 .. size( k )

        cost = real( fobj( real(x) ), kind = core_rknd )

        return
    end subroutine exec_movement

    subroutine init_random

        use error, only: &
            l_use_prescribed_rand_seed, &   ! constant(s) 
            prescribed_rand_seed

        use constants_clubb, only: fstdout

        implicit none

        integer :: i, n, clock
        integer, dimension(:), allocatable :: seed

        call random_seed(size = n)
        allocate(seed(n))

        if ( l_use_prescribed_rand_seed ) then
            ! Use fixed value as random seed
            seed = (/ (prescribed_rand_seed, i = 1, n) /)
        else
            ! Use clock to randomize seed value
            call system_clock(count=clock)
            seed = clock + 37 * (/ (i - 1, i = 1, n) /)
        end if

        write (fstdout,*) "Enhanced Simann random seed = ", seed

        call random_seed(PUT = seed)

        deallocate(seed)

    end subroutine init_random


  !============================== OLD VERSION ==================================
  !-----------------------------------------------------------------------------
  subroutine esa_driver_siarry( xinit, x0min, x0max, tmpini, fobj, xopt, enopt )

  ! Description:
  !   Driver subroutine
  ! References:
  !   None
  !-----------------------------------------------------------------------------
    use mt95, only: genrand_real1 ! Procedure

    use mt95, only: genrand_real ! Constant

    use error, only: nfobj => iter ! Variable

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: &
      epsilon, & ! Machine dependent epsilon
      log,     & ! Log_e( x )
      size       ! Array size

    ! Define the interface for the minimization function `FOBJ' in the paper.
    ! It uses Fortran 90 assumed shape arrays to be compatable with the existing
    ! cost function used with Numerical Recipes' algorithms
#include "enhanced_simann_fobj.inc"

    ! Parameter constants
    integer, parameter :: &
      nfmax = 5000, & ! Stopping parameter from Siarry
      nvariations = 50 ! Number of points in the problem space used for dgyini

    ! Parameters for temperature adjustment
    real( kind = core_rknd ), parameter :: &
      rmxtmp = 0.9_core_rknd, &
      rmitmp = 0.1_core_rknd

    ! Parameters for step vector adjustment
    real( kind = core_rknd ), parameter :: &
      ratmax = 0.2_core_rknd,  &
      ratmin = 0.05_core_rknd, &
      extstp = 2.0_core_rknd,  &
      shrstp = 0.5_core_rknd

    ! Call the cost function 50 times to determine the correct temperature
    ! based on Siarry, et al.
    logical, parameter :: &
      l_compute_optimal_temp = .false. 

    ! Input variables
    real( kind = core_rknd ), dimension(:), intent(in) :: &
      x0min,  & ! Minimum values for the x vector
      x0max,  & ! Maximum values for the x vector
      xinit     ! Initial argument for fobj
   
    real( kind = core_rknd ), intent(inout) :: &
      tmpini ! Initial temperature

    ! Output variables
    real( kind = core_rknd ), dimension(:), intent(out) :: &
      xopt  ! Optimal point for x

    real( kind = core_rknd ), intent(out) :: &
      enopt   ! Optimal value of the cost function

    ! Local variables

    ! Variable names based on the paper
    real( kind = core_rknd ), dimension(size( xinit )) :: &
      xstart, & ! Starting point for x 
      xtry, &   ! Test point for x 
      stpini, & ! Initial step vector
      stpmst    ! Current step vector

    real( kind = core_rknd ) :: &
      epsrel, epsabs, & ! Stop conditions from paper
      einit, & ! Initial value of the cost function
      oldrgy, & ! Old energy level
      rnewgy, & ! New energy level
      deltae, & ! Change in energy, i.e. delta fobj( xtry )
      probok, & ! User chosen initial acceptance probability
      init_avg, &
      dgyini, &
      tstop,  & ! Final temperarture
      temp,   & ! Anneal temperature
      rftmp

    real( kind = core_rknd ), dimension(nvariations) :: &
      init_moves ! Initial vector of moves

    integer, dimension(size( xinit ),4) :: &
      mokst    ! Initial vector of number of accepted moves (one number for each
               ! variable) at the last 4 temperature stages ( 4 stages added by dschanen )

    integer, dimension(size( xinit )) :: &
      mtotst   ! Vector of numbers of attempted moves at current temp stage

    logical, dimension(size( xinit )) :: &
      spartition ! Parts of R^n space in the current partition

    integer :: &
      np,   &! Size of the a partion
      n1, n2, & ! Annealing schedule
!     inorm,  & ! Normalization
!     nfobj, & ! Number of function evaluations
      mvokst, &
      nmvust, & ! Number of accepted uphill moves at current temperature stage
      nmvst    ! Number of attempted moves at current temp stage

    real( kind = core_rknd ) :: &
      rok,    & ! Step vector adjustment variable
      elowst, & ! Minimal fobj value at current temp stage
      avgyst, & ! Sum of successive fobj values at current temp stage
      sdgyup    ! Sum of accepted uphill fobj variations at current temp stage

    real(kind=genrand_real) :: &
      rand ! random number from [0,1]

    logical :: l_accept_xtry

    integer :: i

    integer :: n, nm1, nm2, nm3, ntmp ! Indices for the last 4 temperature stages

      ! ---- Begin Code ----

      ! Step 1: Initializations 

      ! Attempt to make these machine independent (variable names from paper)
      ! epsrel is 1e2 times that of epsabs
      epsrel = 1.e2_core_rknd * epsilon( xinit ) ! Known magic number
      epsabs = epsilon( xinit )

      ! Siarry's epsilon values in the paper
      !epsrel = 10E-6
      !epsabs = 10E-8

      ! Set temperature stage parameters
      np = size( xinit )  ! np = number of variables in x

      ! Suggested values from the paper
      n1 = 12
      n2 = 100

      ! Set inorm for linear normalization
      ! Normalization is handled wthin min_les_clubb_diff
!     inorm = 1

      probok = 0.5_core_rknd ! Probability of accepting an uphill move

      ! Compute initial cost function
      einit = real(fobj( real(xinit) ), kind = core_rknd)

      ! Step vector
      stpini(:) = ( x0max(:) - x0min(:) ) * 0.25_core_rknd

      spartition(:) = .true. ! Partitioning off for average generation

      if ( l_compute_optimal_temp ) then
        ! Obtain the ``variation average'' of fobj, called dgyini.
        ! Is this the standard deviation?  It looks like it should be. 
        ! -dschanen 8 Dec 2008
        do i = 1, size( init_moves ), 1
          xtry = xinit
          call exec_movement_siarry( spartition, x0max, x0min, stpini, fobj, & ! In
                              xtry, init_moves(i) ) ! In/Out, Out
        end do ! 1 .. size of init_moves

        ! Compute std dev
        init_avg = sum( init_moves ) / real( size( init_moves ), kind = core_rknd )
        dgyini = sqrt( sum( ( init_moves-init_avg )**2 ) &
                       / real( size( init_moves ), kind = core_rknd ) &
                     )
        if ( l_esa_debug_statements ) then
          write(6,*) "Initial moves ="
          write(6,'(8g11.4)') init_moves
        end if

        ! Compute initial temperature.  Value of probok comes from Siarry, et al.
        tmpini = -dgyini / log( probok )
      end if

      if ( l_esa_debug_statements ) then
        print *, "Initial temp = ", tmpini
      end if

      ! Formula from pp 220
      ! Something about this formula seems wrong. The tstop value it gives is
      ! generally much too low, and so the algorithm never exits on condition 2.
      ! -dschanen 21 March 2011
!     tstop = -( epsrel * dgyini + epsabs/log( epsrel * probok + epsabs ) )

      tstop = 0.1_core_rknd ! Use a fixed value

      if ( l_esa_debug_statements ) then
        print *, "Stop temp = ", tstop
      end if

      nfobj = 1

      ! Initial optimal point
      xopt(:) = xinit(:)

      ! Initial optimal fobj value
      enopt   = einit

      xstart = xinit 

      oldrgy = einit

      temp = tmpini

      stpmst = stpini

      mvokst = 0 ! Number of accepted moves at current temperature stage

      mokst(:,:) = 0 ! Vector of numbers of accepted moves at current & 3 prior temp stage
      nmvust = 0     ! Number of accepted uphill moves at current temperature stage
      nmvst = 0      ! Number of attempted moves at current temp stage
      mtotst(:) = 0  ! Vector with numbers of attempted moves at current temp stage
      elowst = einit ! Minimal fobj value at current temp stage
      avgyst = 0._core_rknd    ! Sum of successive fobj values at current temp stage
      sdgyup = 0._core_rknd    ! Sum of accepted uphill fobj variations at current temp stage

      ! Added by dschanen for non-exclusive stop test number 1 pp 217
      n   = 1 ! nth temperature stage
      nm1 = 2 ! n-1th temp stage
      nm2 = 3 ! n-2th "    "
      nm3 = 4 ! n-3th "    "

      do  ! The 4 exit conditions for this loop are handled later in the loop

        l_accept_xtry = .false. ! Initialize

        ! Step 2: Space paritioning
        ! This should be random, uniform, and not over select an element of x
        call select_partition( mtotst, &  ! In
                               spartition ) ! Out

        xtry = xstart
        ! Step 3: Execution of One Movement
        call exec_movement_siarry( spartition, x0max, x0min, stpmst, fobj, & ! In
                                xtry, rnewgy ) ! In/Out, Out

        !print *, "x = ", xtry, "f(x) = ", rnewgy, "xopt", xopt
        !pause

        deltae = rnewgy - oldrgy
        where ( spartition ) mtotst = mtotst + 1
        nmvst = nmvst + 1
        nfobj = nfobj + 1
        avgyst = avgyst + rnewgy

        ! Step 4: Acceptance or Rejection of this movement
        if ( deltae <= 0.0_core_rknd ) then  
            ! Accept xtry
            l_accept_xtry = .true.
            if ( rnewgy < enopt ) then
              xopt = xtry
              enopt = rnewgy
            end if
            if ( rnewgy < elowst ) then
              elowst = rnewgy
            end if
         else
           call genrand_real1( rand )

           ! Accept the number with probability of exp(-deltae / temp)
           if ( real( rand, kind = core_rknd ) <= exp( -deltae/temp ) ) then
             ! Accept xtry
             l_accept_xtry = .true.
             nmvust = nmvust + 1
             sdgyup = sdgyup + deltae
           end if
        end if ! deltae <= 0.0_core_rknd

        if ( l_accept_xtry ) then
          xstart = xtry
          oldrgy = rnewgy
          where ( spartition ) mokst(:,n) = mokst(:,n) + 1
          mvokst = mvokst + 1
        end if

        ! Step 5: Test for the End of the Temperature Stage

        if ( mvokst <  n1*np .and. nmvst < n2*np ) then
          !print *, "cycling"
          cycle 
        else
          ! Cycle indices on mokst (accepted moves)
          ntmp = nm3
          nm3  = nm2
          nm2  = nm1
          nm1  = n
          n    = ntmp

          mokst(:,n) = 0 ! Overwrite the old n-3

          !print *, "new temperature level"
        end if

        ! Step 6: Temperature Adjustment

        avgyst = avgyst / real( nmvst, kind = core_rknd )
        rftmp = max( min( elowst/avgyst, rmxtmp), rmitmp )
        temp = rftmp * temp

        ! Step 7: Step Vector Adjustment
        do i = 1, np, 1
          if ( spartition(i) ) then
            rok = real( mokst(i,n), kind = core_rknd ) / real( mtotst(i), kind = core_rknd )
            if ( rok > ratmax ) then
              stpmst(i) = stpmst(i) * extstp
            else if ( rok < ratmin ) then
              stpmst(i) = stpmst(i) * shrstp
            end if
          end if
        end do ! i = 1 .. np

        ! Step 8: Four non-exclusive stopping tests 

        ! (1) No uphill changes in the last 4 temperature stages
        if ( all( mokst(:,1:4) == 0 ) ) then

          if ( l_esa_debug_statements ) then
            write(6,*) "Exit on condition 1, No uphill changes for the last 4 temperature stages"
            write(6,*) "temp", temp
            write(6,*) "stpmst", stpmst
            write(6,*) "stpini", stpini
          end if

          exit
        end if
           
        ! (2) temp < tstop
        if ( temp < tstop ) then

          if ( l_esa_debug_statements ) then
            write(6,*) "Exit on condition 2, temperature < tstop"
            write(6,*) "temp", temp
            write(6,*) "stpmst", stpmst
            write(6,*) "stpini", stpini
          end if

          exit
        end if

        ! (3) Too close to machine epsilon
        if ( any( stpmst < epsrel * stpini + epsabs ) ) then

          if ( l_esa_debug_statements ) then
            write(6,*) "Exit on condition 3, too close to machine epsilon"
            write(6,*) "temp", temp
            write(6,*) "stpmst", stpmst
            write(6,*) "stpini", stpini
          end if

          exit
        end if

        ! (4) Too many iterations
        if ( nfobj >= nfmax*np ) then

          if ( l_esa_debug_statements ) then
            write(6,*) "Exit on condition 4, too many iterations"
          end if

          exit
        end if

        ! Step 9: Initialization of the new temperature stage
        mvokst = 0
        nmvust = 0
        nmvst  = 0
        elowst = oldrgy
        avgyst = 0._core_rknd
        sdgyup = 0._core_rknd

      end do

      if ( l_esa_debug_statements ) then
        print *, "mtotst = ", mtotst
      end if

    return
  end subroutine esa_driver_siarry

  !-----------------------------------------------------------------------------
  subroutine exec_movement_siarry( spartition, x0max, x0min, step, fobj, &
                            x, cost )
  ! Description:
  !   Execute one movement in domain space

  ! References:
  !   Step 3 in Siarry, et al.
  !-----------------------------------------------------------------------------

    use mt95, only: genrand_real1 ! Procedure

    use mt95, only: genrand_real ! Constant

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: size

#include "enhanced_simann_fobj.inc"

    ! Input Variables
    logical, dimension(:), intent(in) :: spartition

    real( kind = core_rknd ), dimension(:), intent(in) :: x0max, x0min, step

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(:), intent(inout) :: x

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: cost

    ! Local variables
    real(kind=genrand_real), dimension(size( x )) :: xrand, srand

    real( kind = core_rknd ) :: xtmp

    integer :: k

    ! --- Begin Code ---

    ! The random range should be from [0,1]

    call genrand_real1( xrand ) ! Used for the value of x

    ! There is probably a more efficient way to get a random bit in Fortran, but this is easy
    call genrand_real1( srand ) ! Used for the sign of x

    do k = 1, size( x )

      if ( spartition(k) ) then ! Augment only values within the partition

        ! Apply a random sign to each xrand value. 
        if ( srand(k) >= 0.5_genrand_real ) then
          xtmp = x(k) + real( xrand(k), kind = core_rknd ) * step(k)
        else
          xtmp = x(k) - real( xrand(k), kind = core_rknd ) * step(k)
        end if

        ! According to Siarry pp 216 if a number is outside the x range, we
        ! should change the sign of the step accordingly
        if ( xtmp > x0max(k) ) then
          xtmp = x(k) - real( xrand(k), kind = core_rknd ) * step(k)
        else if ( xtmp < x0min(k) ) then
          xtmp = x(k) + real( xrand(k), kind = core_rknd ) * step(k)
        end if

        x(k) = xtmp

      end if  ! Points in the partition

    end do ! 1 .. size( k )

    cost = real(fobj( real(x) ), kind = core_rknd)

    return
  end subroutine exec_movement_siarry

  !-----------------------------------------------------------------------------
  subroutine select_partition( mtotst, spartition )
  ! Description:
  !   Select a partition of p dimension within x
  ! References:
  !  pp. 221-222, Siarry et al. ``SA Parameter Adjustment''
  !-----------------------------------------------------------------------------

    use mt95, only: genrand_int31 ! Procedure

    use mt95, only: genrand_intg ! Constant

    implicit none

    intrinsic :: &
      size, maxval, all, count, mod

    ! Input Variables
    integer, dimension(:), intent(in) :: &
      mtotst ! Vector of which variables were tried in prior iterations

    ! Output Variables
    logical, dimension(:), intent(out) :: &
      spartition  ! Vector of which variables to include

    ! Local Variables
    integer(kind=genrand_intg) :: irand

    integer :: pdim, ndim, elem

    !---- Begin Code ----
    ndim = size( mtotst ) ! size( mtotst ) = size( x )

    ! According to experiments done by Siarry p = n/3 is most efficient for 
    ! complex circuit problems, so for now we'll just use that.
    pdim = max( ndim / 3, 1 ) 

    spartition(:) = .false.

    ! Terminate when we have pdim true entries
    do while ( count( spartition ) < pdim ) 

      call genrand_int31( irand )

      elem = mod( irand, ndim ) + 1  ! Pick a random element

      !print *, elem
      !pause
      ! Attempt to meet Siarry's condition that no element be over-selected
!     if ( mtotst(elem) /= maxval( mtotst ) .or. &
!          all( mtotst == maxval( mtotst ) ) ) then
      spartition(elem) = .true.
!     else
!       print *, mtotst, elem, rand
!     end if

    end do

    return
  end subroutine select_partition

end module enhanced_simann


