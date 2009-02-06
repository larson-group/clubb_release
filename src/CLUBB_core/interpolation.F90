!-------------------------------------------------------------------------------
!$Id$
module interpolation

  implicit none

  private ! Default Scope

  public :: lin_int, binary_search, zlinterp_fnc, & 
    linear_interpolation, factor_interp

  contains

!-------------------------------------------------------------------------------
  pure function lin_int( height_int, height_high, height_low, &
    var_high, var_low )

! Description:
! This function computes a linear interpolation of the value of variable.
! Given two known values of a variable at two height values, the value
! of that variable at a height between those two height levels (rather 
! than a height outside of those two height levels) is computed.
!
! Here is a diagram:
!
!  ################################ Height high, know variable value
!
!
!
!  -------------------------------- Height to be interpolated to; linear interpolation
!
!
!
!
!
!  ################################ Height low, know variable value
!
!
! FORMULA:
!
! variable(@ Height interpolation) =
!
! [ (variable(@ Height high) - variable(@ Height low)) / (Height high - Height low) ]
! * (Height interpolation - Height low)  +  variable(@ Height low)

! Comments from WRF-HOC, Brian Griffin.

! References:
! None
!-------------------------------------------------------------------------------

    implicit none

    ! Input Variables

    real, intent(in) :: &
      height_int,  & ! Height to be interpolated to     [m]
      height_high, & ! Height above the interpolation   [m]
      height_low,  & ! Height below the interpolation   [m]
      var_high,    & ! Variable above the interpolation [units vary]
      var_low        ! Variable below the interpolation [units vary]

    ! Output Variables
    real :: lin_int

    ! Compute linear interpolation

    lin_int = ( ( height_int - height_low )/( height_high - height_low ) ) &
      * ( var_high - var_low ) + var_low

    return
  end function lin_int

!-------------------------------------------------------------------------------------------------
 elemental real function factor_interp( factor, var_high, var_low )
!-------------------------------------------------------------------------------------------------
    implicit none

    real, intent(in) :: &
    factor,   & ! Factor                           [units vary]  
    var_high, & ! Variable above the interpolation [units vary]
    var_low     ! Variable below the interpolation [units vary]

    factor_interp = factor * ( var_high - var_low ) + var_low

    return
  end function factor_interp

!-------------------------------------------------------------------------------
  pure integer function binary_search( n, array, var ) & 
    result( i ) 
!       Description: This subroutine performs a binary search to find
!       the closest value greater than or equal to var in the array.
!        
!-------------------------------------------------------------------------------
    implicit none

    ! Input Variables

    ! Size of the array
    integer, intent(in) :: n

    ! The array being searched
    real, dimension(n), intent(in) :: array

    ! The value being searched for
    real, intent(in) :: var

    ! Local Variables

    ! Has an index been found?
    logical :: l_found

    ! Bounds of the search
    integer :: high
    integer :: low

    ! Initialize local variables

    l_found = .false.

    low = 1

    high = n 

    i = (low + high) / 2

    do while( .not. l_found .and. low <= high )
        
      i = (low + high) / 2

      if ( var > array( i - 1 ) .and. var <= array( i ) )then
        l_found = .true.

      else if (var < array( i ) )then
        high = i - 1
      else if (var > array( i ) )then
        low = i + 1
      end if              

    end do  ! while ( ~l_found & low <= high )
        
    if ( .not. l_found ) i = -1 

    return
  
    end function binary_search
!-------------------------------------------------------------------------------
  function zlinterp_fnc( dim_out, dim_src, grid_out,  & 
                       grid_src, var_src )  & 
  result( var_out )
! Description:
!   Do a linear interpolation in the vertical.  Assumes values that
!   are less than lowest source point are zero and above the highest
!   source point are zero. Also assumes altitude increases linearly.

! References:
!   function LIN_INT from WRF-HOC
!-----------------------------------------------------------------------

    implicit none

    ! Input variables
    integer, intent(in) :: dim_out, dim_src

    real, dimension(dim_src), intent(in) ::  & 
      grid_src,  & ! [m]
      var_src      ! [units vary]

    real, dimension(dim_out), intent(in) :: &
      grid_out ! [m]

    ! Output variable
    real, dimension(dim_out) :: &
      var_out ! [units vary]

    ! Local variables
    integer :: k, kint, km1

!   integer :: tst, kp1

    k = 1

    do kint = 1, dim_out, 1

      ! Set to 0 if we're below the input data's lowest point
      if ( grid_out(kint) < grid_src(1) ) then
        var_out(kint) = 0.0
        cycle
      end if

  ! Increment k until the level is correct
!          do while ( grid_out(kint) > grid_src(k) 
!     .                .and. k < dim_src )
!            k = k + 1
!          end do

  ! Changed so a binary search is used instead of a sequential search
!          tst = binary_search(dim_src, grid_src, grid_out(kint))
      k = binary_search(dim_src, grid_src, grid_out(kint))
  ! Joshua Fasching April 2008

!          print *, "k = ", k
!          print *, "tst = ", tst
!          print *, "dim_src = ", dim_src
!          print *,"------------------------------"
  
    ! If the increment leads to a level above the data, set this
    ! point and all those above it to zero
      !if( k > dim_src ) then
      if ( k == -1 ) then
        var_out(kint:dim_out) = 0.0
        exit
      end if

      km1 = max( 1, k-1 )
      !kp1 = min( k+1, dim_src )

      ! Interpolate
      var_out(kint) = lin_int( grid_out(kint), grid_src(k),  & 
        grid_src(km1), var_src(k), var_src(km1) )
  
!          ( var_src(k) - var_src(km1) ) / &
!          ( grid_src(k) - grid_src(km1) ) &
!            * ( grid_out(kint) - grid_src(km1) ) + var_src(km1) &
!            Changed to use a standard function for interpolation

     !! Note this ends up changing the results slightly because
     !the placement of variables has been changed.
  
!            Joshua Fasching April 2008

    end do ! kint = 1..dim_out

    return
  end function zlinterp_fnc

!-------------------------------------------------------------------------------
  subroutine linear_interpolation & 
             ( nparam, xlist, tlist, xvalue, tvalue )

! Description:
!   Linear interpolation for 25 June 1996 altocumulus case.  

!   For example, to interpolate between two temperatures in space, put
!   your spatial coordinates in x-list and your temperature values in
!   tlist.  The point in question should have its spatial value stored
!   in xvalue, and tvalue will be the temperature at that point.

! Author: Michael Falk for COAMPS.
!-------------------------------------------------------------------------------

    use error_code, only: clubb_debug ! Procedure

    implicit none

    ! Input Variables
    integer, intent(in) :: nparam ! Number of parameters in xlist and tlist

    ! Input/Output Variables
    real, intent(inout), dimension(nparam) ::  & 
      xlist,  & ! List of x-values (independent variable)
      tlist     ! List of t-values (dependent variable)

    real, intent(inout) ::  & 
      xvalue,  & ! x-value at which to interpolate
      tvalue  ! t-value solved by interpolation

    ! Local variables
    integer ::  & 
      i,  & ! Loop control variable for bubble sort- number of the 
            ! lowest yet-unsorted data point.
      j  ! Loop control variable for bubble sort- index of value
         ! currently being tested
    integer ::  & 
      bottombound, & ! Index of the smaller value in the linear interpolation
      topbound,    & ! Index of the larger value in the linear interpolation
      smallest       ! Index of the present smallest value, for bubble sort

    real :: temp ! A temporary variable used for the bubble sort swap

!-------------------------------------------------------------------------------
!
! Bubble Sort algorithm, assuring that the elements are in order so
! that the interpolation is between the two closest points to the
! point in question.
!
!-------------------------------------------------------------------------------

    do i=1,nparam
      smallest = i
      do j=i,nparam
        if ( xlist(j) < xlist(smallest) ) then
          smallest = j
        end if
      end do

      temp = xlist(i)
      xlist(i) = xlist(smallest)
      xlist(smallest) = temp

      temp = tlist(i)
      tlist(i) = tlist(smallest)
      tlist(smallest) = temp
    end do

!-------------------------------------------------------------------------------
!
! If the point in question is larger than the largest x-value or
! smaller than the smallest x-value, crash.
!
!-------------------------------------------------------------------------------

    if ( (xvalue < xlist(1)) .or. (xvalue > xlist(nparam)) ) then
      write(0,'(a)') "linear_interpolation: Value out of range"
      stop
    end if

!-------------------------------------------------------------------------------
!
! Find the correct top and bottom bounds, do the interpolation, return c
! the value.
!
!-------------------------------------------------------------------------------

    topbound = -1
    bottombound = -1

    do i=2,nparam
      if ( (xvalue >= xlist(i-1)) .and. (xvalue <= xlist(i)) ) then
        bottombound = i-1
        topbound    = i
      end if
    end do

    if ( topbound == -1 .or. bottombound == -1 ) then
      call clubb_debug( 1, "Sanity check failed! xlist is not properly sorted" )
      call clubb_debug( 1, "in linear_interpolation.")
    end if

    tvalue =  & 
    lin_int( xvalue, xlist(topbound), xlist(bottombound),  & 
            tlist(topbound), tlist(bottombound) )       

    return
  end subroutine linear_interpolation
 
end module interpolation
