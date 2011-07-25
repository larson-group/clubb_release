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
!-------------------------------------------------------------------------------

  implicit none

  public :: esa_driver

  private :: select_partition, exec_movement

  private ! Default scope

  logical, parameter :: &
    l_esa_debug_statements = .false. ! Verbose debugging output

  contains

  !-----------------------------------------------------------------------------
  subroutine esa_driver( xinit, x0min, x0max, rostep, tmpini, fobj, xopt, enopt )

  ! Description:
  !   Driver subroutine
  ! References:
  !   None
  !-----------------------------------------------------------------------------
    use mt95, only: genrand_real1 ! Procedure

    use mt95, only: genrand_real ! Constant

    use error, only: nfobj => iter ! Variable

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
    real, parameter :: &
      rmxtmp = 0.9, &
      rmitmp = 0.1

    ! Parameters for step vector adjustment
    real, parameter :: &
      ratmax = 0.2,  &
      ratmin = 0.05, &
      extstp = 2.0,  &
      shrstp = 0.5

    ! Call the cost function 50 times to determine the correct temperature
    ! based on Siarry, et al.
    logical, parameter :: &
      l_compute_optimal_temp = .false. 

    ! Input variables
    real, dimension(:), intent(in) :: &
      x0min,  & ! Minimum values for the x vector
      x0max,  & ! Maximum values for the x vector
      rostep, & ! Increments for step vector
      xinit     ! Initial argument for fobj
   
    real, intent(inout) :: &
      tmpini ! Initial temperature

    ! Output variables
    real, dimension(:), intent(out) :: &
      xopt  ! Optimal point for x

    real, intent(out) :: &
      enopt   ! Optimal value of the cost function

    ! Local variables

    ! Variable names based on the paper
    real, dimension(size( xinit )) :: &
      xstart, & ! Starting point for x 
      xtry, &   ! Test point for x 
      stpini, & ! Initial step vector
      stpmst    ! Current step vector

    real :: &
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

    real, dimension(nvariations) :: &
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

    real :: &
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
      epsrel = 1.e2 * epsilon( xinit ) ! Known magic number
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

      probok = 0.5 ! Probability of accepting an uphill move

      ! Compute initial cost function
      einit = fobj( xinit )

      ! Suggested value from the paper
!     rostep(:) = 0.25

      ! Step vector
      stpini(:) = ( x0max(:) - x0min(:) ) * rostep

      spartition(:) = .true. ! Partitioning off for average generation

      if ( l_compute_optimal_temp ) then
        ! Obtain the ``variation average'' of fobj, called dgyini.
        ! Is this the standard deviation?  It looks like it should be. 
        ! -dschanen 8 Dec 2008
        do i = 1, size( init_moves ), 1
          xtry = xinit
          call exec_movement( spartition, x0max, x0min, stpini, fobj, & ! In
                              xtry, init_moves(i) ) ! In/Out, Out
        end do ! 1 .. size of init_moves

        ! Compute std dev
        init_avg = sum( init_moves ) / real( size( init_moves ) )
        dgyini = sqrt( sum( ( init_moves-init_avg )**2 ) &
                       / real( size( init_moves ) ) &
                     )
        if ( l_esa_debug_statements ) then
          write(6,*) "Intial moves ="
          write(6,'(8g10.4)') init_moves
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

      tstop = 0.1 ! Use a fixed value

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
      avgyst = 0.    ! Sum of successive fobj values at current temp stage
      sdgyup = 0.    ! Sum of accepted uphill fobj variations at current temp stage

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

        !spartition(:) = .true. ! uncomment to turn off partitioning

        xtry = xstart
        ! Step 3: Execution of One Movement
        call exec_movement( spartition, x0max, x0min, stpmst, fobj, & ! In
                            xtry, rnewgy ) ! In/Out, Out

        !print *, "x = ", xtry, "f(x) = ", rnewgy, "xopt", xopt
        !pause

        deltae = rnewgy - oldrgy
        where ( spartition ) mtotst = mtotst + 1
        nmvst = nmvst + 1
        nfobj = nfobj + 1
        avgyst = avgyst + rnewgy

        ! Step 4: Acceptance or Rejection of this movement
        if ( deltae <= 0.0 ) then  
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
           if ( rand  <= exp( -deltae/temp ) ) then
             ! Accept xtry
             l_accept_xtry = .true.
             nmvust = nmvust + 1
             sdgyup = sdgyup + deltae
           end if
        end if ! deltae <= 0.0

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

        avgyst = avgyst / nmvst
        rftmp = max( min( elowst/avgyst, rmxtmp), rmitmp )
        temp = rftmp * temp

        ! Step 7: Step Vector Adjustment
        do i = 1, np, 1
          if ( spartition(i) ) then
            rok = real( mokst(i,n) ) / real( mtotst(i) )
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
        avgyst = 0
        sdgyup = 0

      end do

      if ( l_esa_debug_statements ) then
        print *, "mtotst = ", mtotst
      end if

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

    implicit none

    ! External
    intrinsic :: size

#include "enhanced_simann_fobj.inc"

    ! Input Variables
    logical, dimension(:), intent(in) :: spartition

    real, dimension(:), intent(in) :: x0max, x0min, step

    ! Input/Output Variables
    real, dimension(:), intent(inout) :: x

    ! Output Variables
    real, intent(out) :: cost

    ! Local variables
    real(kind=genrand_real), dimension(size( x )) :: xrand, srand

    real :: xtmp

    integer :: k

    ! --- Begin Code ---

    ! The random range should be from [0,1]

    call genrand_real1( xrand ) ! Used for the value of x

    ! There is probably a more efficient way to get a random bit in Fortran, but this is easy
    call genrand_real1( srand ) ! Used for the sign of x

    do k = 1, size( x )

      if ( spartition(k) ) then ! Augment only values within the partition

        ! Apply a random sign to each xrand value. 
        if ( srand(k) >= 0.5 ) then
          xtmp = x(k) + real( xrand(k) ) * step(k)
        else
          xtmp = x(k) - real( xrand(k) ) * step(k)
        end if

        ! According to Siarry pp 216 if a number is outside the x range, we
        ! should change the sign of the step accordingly
        if ( xtmp > x0max(k) ) then
          xtmp = x(k) - real( xrand(k) ) * step(k)
        else if ( xtmp < x0min(k) ) then
          xtmp = x(k) + real( xrand(k) ) * step(k)
        end if

        x(k) = xtmp

      end if  ! Points in the partition

    end do ! 1 .. size( k )

    cost = fobj( x )

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
