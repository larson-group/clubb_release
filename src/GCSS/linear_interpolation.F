! $Id: linear_interpolation.F,v 1.2 2007-04-24 21:16:44 dschanen Exp $
!-----------------------------------------------------------------------
!      program interp
!      Description:
!      Program interp, to drive subroutine linearinterpolation for testing.
!
!      Only uncomment the following code if you need it to test the
!      interpolation scheme.

!      Author: Michael Falk for COAMPS, 16 March 2005.
!-----------------------------------------------------------------------
!      integer nparam
!      real xlist(nparam),tlist(nparam),xvalue,tvalue
!
!      nparam=5
!
!      xlist(1)=4
!      xlist(2)=3
!      xlist(3)=2
!      xlist(4)=1
!      xlist(5)=2.5
!
!      tlist(1)=400
!      tlist(2)=300
!      tlist(3)=200
!      tlist(4)=100
!      tlist(5)=203
!
!      xvalue = 2.4
!
!      call linearinterpolation(nparam,xlist,tlist,xvalue,tvalue)
!      print *,'Tvalue:', tvalue
!      stop
!      end program interp

!-----------------------------------------------------------------------
      subroutine linear_interpolation
     .           ( nparam, xlist, tlist, xvalue, tvalue )

!     Description:
!     Linear interpolation for 25 June 1996 altocumulus case.  

!     For example, to interpolate between two temperatures in space, put
!     your spatial coordinates in x-list and your temperature values in
!     tlist.  The point in question should have its spatial value stored
!     in xvalue, and tvalue will be the temperature at that point.

!     Author: Michael Falk for COAMPS.
!-----------------------------------------------------------------------

      implicit none

      ! Input Variables
      integer, intent(in) :: nparam ! Number of parameters in xlist and tlist

      ! Input/Output Variables
      real, intent(inout), dimension(nparam) :: 
     .  xlist, ! List of x-values (independent variable)
     .  tlist  ! List of t-values (dependent variable)

      real, intent(inout) :: 
     .  xvalue, ! x-value at which to interpolate
     .  tvalue  ! t-value solved by interpolation

      ! Local variables
      integer :: 
     .  i, ! Loop control variable for bubble sort- number of the 
           ! lowest yet-unsorted data point.
     .  j  ! Loop control variable for bubble sort- index of value
           ! currently being tested
      integer :: 
     .  bottombound, ! Index of the smaller value in the linear interpolation
     .  topbound,    ! Index of the larger value in the linear interpolation
     .  smallest     ! Index of the present smallest value, for bubble sort

      real :: temp ! A temporary variable used for the bubble sort swap


!-----------------------------------------------------------------------
!
! Bubble Sort algorithm, assuring that the elements are in order so
! that the interpolation is between the two closest points to the
! point in question.
!
!-----------------------------------------------------------------------

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

!-----------------------------------------------------------------------
!
! If the point in question is larger than the largest x-value or
! smaller than the smallest x-value, crash.
!
!-----------------------------------------------------------------------

      if ( (xvalue < xlist(1)) .or. (xvalue > xlist(nparam)) ) then
        write(0,'(a)') "linear_interpolation: Value out of range"
        stop
      end if

!-----------------------------------------------------------------------
!
! Find the correct top and bottom bounds, do the interpolation, return c
! the value.
!
!-----------------------------------------------------------------------

      do i=2,nparam
        if ( (xvalue >= xlist(i-1)) .and. (xvalue <= xlist(i)) ) then
          bottombound = i-1
          topbound    = i
        end if
      end do

      tvalue = (xvalue - xlist(bottombound)) /
     .         (xlist(topbound) - xlist(bottombound)) *
     .         (tlist(topbound) - tlist(bottombound)) +
     .         tlist(bottombound)

      return
      end subroutine linear_interpolation
