!-------------------------------------------------------------------------------
! $Id$
module enhanced_simann
! Description:
!   Implementation of Siarry's Enhanced simulated annealing algorthm in
!   Fortran 90/95.

! References: 
!   ``Enhanced Simulated Annealing for Many Globally Minimized Functions
!   of Many-Continuous Variables'', Siarry, et al. ACMS TOMS Vol. 23,
!   No. 2, June 1997, pp. 209--228.
!------------------------------------------------------------------------------

  implicit none

  public :: esa_driver

  private :: select_partition, exec_movement

  private ! Default scope

  contains

  !-----------------------------------------------------------------------------
  subroutine esa_driver( xinit, x0min, x0max, rostep, tmpini, fobj, xopt, nrgy_opt )

  ! Description:
  !   Driver subroutine
  ! References:
  !   None
  !-----------------------------------------------------------------------------
    use mt95, only: genrand_real1 ! Procedure

    use mt95, only: genrand_real ! Constant

    use error, only: nfobj => iter ! Variable

    use text_writer, only: write_text ! Subroutine(s)

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! External
    !intrinsic :: &
    !  epsilon, & ! Machine dependent epsilon
    !  log,     & ! Log_e( x )
    !  size       ! Array size

    ! Define the interface for the minimization function `FOBJ' in the paper.
    ! It uses Fortran 90 assumed shape arrays to be compatable with the existing
    ! cost function used with Numerical Recipes' algorithms
#include "enhanced_simann_fobj.inc"

    ! Input variables
    real( kind = core_rknd ), dimension(:), intent(in) :: &
      x0min,  & ! Minimum values for the x vector
      x0max,  & ! Maximum values for the x vector
      rostep, & ! Increments for step vector
      xinit     ! Initial argument for fobj
   
    real( kind = core_rknd ), intent(inout) :: &
      tmpini ! Initial temperature

    ! Output variables
    real( kind = core_rknd ), dimension(:), intent(out) :: &
      xopt  ! Optimal point for x

    real( kind = core_rknd ), intent(out) :: &
      nrgy_opt   ! Optimal value of the cost function


    ! Local variables


    ! Parameter constants
    integer, parameter :: &
      nfmax = 5000  ! Stopping parameter from Siarry

    ! Parameters for temperature adjustment
    real( kind = core_rknd ), parameter :: &
      rmxtmp = 0.9_core_rknd, &
      rmitmp = 0.1_core_rknd

    ! Parameters for step vector adjustment
    real( kind = core_rknd ), parameter :: &
       ! If variable is accepted more than this percent of the time, its step will be extended
      extend_chance = 0.8_core_rknd,  &    
      ! If variable is accepted less than this percent of the time, its step will be shortened
      shorten_chance = 0.2_core_rknd, &    
      extstp = 1.1_core_rknd,  &
      shrstp = 0.5_core_rknd


    ! Variable names based not on the paper because they were stupid
    real( kind = core_rknd ), dimension(size( xinit )) :: &
      xstart, & ! Starting point for x 
      xtry, &   ! Test point for x 
      stp_ini, & ! Initial step vector
      stp_cur    ! Current step vector

    real( kind = core_rknd ) :: &
      einit, & ! Initial value of the cost function
      old_nrgy, & ! Old energy level
      new_nrgy, & ! New energy level
      deltae, & ! Change in energy, i.e. delta fobj( xtry )
      tstop,  & ! Final temperarture
      temp      ! Anneal temperature

    integer, dimension(size( xinit )) :: &
      accepted, &
      attempted    ! Vector of numbers of attempted moves at current temp stage

    logical, dimension(size( xinit )) :: &
      spartition ! Parts of R^n space in the current partition

    integer :: &
      np,   &! Size of the a partion
      n1, n2 ! Annealing schedule

    real( kind = core_rknd ) :: &
      percent_accepted,       & ! Step vector adjustment variable
      nrgy_min,  &              ! Minimal fobj value at current temp stage
      nrgy_tot                  ! Sum of successive fobj values at current temp stage

    real(kind=genrand_real) :: &
      rand ! random number from [0,1]

    integer :: i

      ! ---- Begin Code ----

      ! Step 1: Initializations 

      ! Set temperature stage parameters
      np = size( xinit )  ! np = number of variables in x

      ! Suggested values from the paper
      n1 = 12
      n2 = 100

      ! Compute initial cost function
      einit = real(fobj( real(xinit) ), kind = core_rknd)

      ! Step vector
      stp_ini(:) = ( x0max(:) - x0min(:) ) * rostep

      spartition(:) = .true. ! Partitioning off for average generation

      tstop = 0.1e-5_core_rknd ! Use a fixed value

      nfobj = 1

      ! Initial optimal point
      xopt(:) = xinit(:)

      ! Initial optimal fobj value
      nrgy_opt  = einit
      xstart    = xinit 
      old_nrgy  = einit
      temp      = tmpini
      stp_cur   = stp_ini

      accepted(:) = 0               ! Number of accepted moves at current stage for each variable
      attempted(:) = 0              ! Vector with numbers of attempted moves at current temp stage
      nrgy_min = einit              ! Minimal fobj value at current temp stage
      nrgy_tot = 0._core_rknd       ! Sum of successive fobj values at current temp stage

      do  ! The 4 exit conditions for this loop are handled later in the loop


        ! Step 2: Space paritioning
        ! This should be random, uniform, and not over select an element of x
        call select_partition( attempted, &  ! In             ! choose what variables to change 
                               spartition ) ! Out


        ! Step 3: Execution of One Movement
        xtry = xstart

        call exec_movement( spartition, x0max, x0min, stp_cur, fobj, & ! In  
                            xtry, new_nrgy ) ! In/Out, Out

        deltae = new_nrgy - old_nrgy                            ! calc difference in energy
        where ( spartition ) attempted = attempted + 1                                           
        nfobj = nfobj + 1                                   ! increment iteration count
        nrgy_tot = nrgy_tot + new_nrgy



        ! Step 4: Acceptance or Rejection of this movement
        if ( deltae <= 0.0_core_rknd ) then              ! if energy decreased (cost reduced)
    
            xstart = xtry                              ! save the values for next iteration
            old_nrgy = new_nrgy
            where(spartition) accepted = accepted + 1

            if ( new_nrgy < nrgy_opt ) then                  ! if this is the absolute minimum
              xopt = xtry                                ! save the values as final, x and cost
              nrgy_opt = new_nrgy
            end if

            if ( new_nrgy < nrgy_min ) then               ! if it's a local minimum
              nrgy_min = new_nrgy                                     ! save the cost
            end if

        else
           call genrand_real1( rand )                ! if energy increased (cost increased)
            ! Accept the number with probability of exp(-deltae / temp)
           if ( real( rand, kind = core_rknd ) <= exp( -deltae/temp ) ) then

             xstart = xtry                                  ! save the values for next iteration
             old_nrgy = new_nrgy
             where(spartition) accepted = accepted + 1

           end if
        end if ! deltae <= 0.0_core_rknd



        ! Step 5: Test for the End of the Temperature Stage

        ! if accepted < 12 * vars  AND  attempted < 100 * vars
        if ( sum(accepted) <  n1*np .and. sum(attempted) < n2*np ) then  
          cycle                     ! then start loop over (no temp or step adjust)
        end if



        ! Step 6: Temperature Adjustment
        ! max( min(low/average,tempmin), tempmax )       
        temp = temp * max( min( nrgy_min * real( sum(attempted), kind = core_rknd ) &
                    / nrgy_tot, rmxtmp), rmitmp )           




        ! Step 7: Step Vector Adjustment
        do i = 1, np, 1      ! for every variable

            percent_accepted = real( accepted(i), kind = core_rknd )  &
                             / real( attempted(i), kind = core_rknd )                        

            if ( percent_accepted > extend_chance ) then
              stp_cur(i) = stp_cur(i) * extstp

                if( stp_cur(i) > (x0max(i) - x0min(i))/2.0 ) then
                    stp_cur = stp_ini
                end if

            else if ( percent_accepted < shorten_chance ) then
              stp_cur(i) = stp_cur(i) * shrstp
            end if

            

        end do ! i = 1 .. np




        ! Step 8: Four non-exclusive stopping tests 
           
        ! (2) temp < tstop
        if ( temp < tstop ) then
          exit
        end if

        ! (4) Too many iterations
        if ( nfobj >= nfmax*np ) then
          exit
        end if



        ! Step 9: Initialization of the new temperature stage
        accepted = 0
        attempted  = 0
        nrgy_tot = 0._core_rknd

        ! go back to best params when in new temp stage
        xstart = xopt
        nrgy_min = nrgy_opt

      end do

    return
  end subroutine esa_driver

  !-----------------------------------------------------------------------------
  subroutine exec_movement( spartition, x0max, x0min, step, fobj, &
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
          xtmp = x0max(k)
          !xtmp = x(k) - real( xrand(k), kind = core_rknd ) * step(k)
        else if ( xtmp < x0min(k) ) then
          xtmp = x0min(k)
          !xtmp = x(k) + real( xrand(k), kind = core_rknd ) * step(k)
        end if

        x(k) = xtmp

      end if  ! Points in the partition

    end do ! 1 .. size( k )

    cost = real(fobj( real(x) ), kind = core_rknd)

    return
  end subroutine exec_movement

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
