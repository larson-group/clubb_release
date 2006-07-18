cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Program interp, to drive subroutine linearinterpolation for testing. c
c Written by Michael Falk for COAMPS, 16 March 2005.                   c
c                                                                      c
c Only uncomment the following code if you need it to test the         c
c interpolation scheme.                                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c      program interp
c      integer nparam
c      real xlist(nparam),tlist(nparam),xvalue,tvalue
c
c      nparam=5
c
c      xlist(1)=4
c      xlist(2)=3
c      xlist(3)=2
c      xlist(4)=1
c      xlist(5)=2.5
c
c      tlist(1)=400
c      tlist(2)=300
c      tlist(3)=200
c      tlist(4)=100
c      tlist(5)=203
c
c      xvalue = 2.4
c
c      call linearinterpolation(nparam,xlist,tlist,xvalue,tvalue)
c      print *,'Tvalue:', tvalue
c      stop
c      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Subroutine linearinterpolation.F -- for 25 June 1996 altocumulus     c
c case.  Written by Michael Falk for COAMPS.                           c
c                                                                      c
c ARGUMENTS (IN):                                                      c
c nparam : Number of parameters in xlist and tlist                     c
c xlist  : List of x-values (independent variable)                     c
c tlist  : List of t-values (dependent variable)                       c
c xvalue : x-value at which to interpolate                             c
c                                                                      c
c ARGUMENTS (OUT):                                                     c
c tvalue : t-value solved by interpolation                             c
c                                                                      c
c For example, to interpolate between two temperatures in space, put   c
c your spatial coordinates in x-list and your temperature values in    c
c tlist.  The point in question should have its spatial value stored   c
c in xvalue, and tvalue will be the temperature at that point.         c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine linearinterpolation(nparam,xlist,tlist,xvalue,tvalue)
      integer nparam
      real xlist(nparam)
      real tlist(nparam)
      real xvalue
      real tvalue

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c LOCAL VARIABLES                                                      c
c                                                                      c
c i          : Loop control variable for bubble sort- number of lowest c
c              yet-unsorted data point                                 c
c j          : Loop control variable for bubble sort- index of value   c
c              currently being tested                                  c
c bottombound: Index of the smaller value in the linear interpolation  c
c topbound   : Index of the larger value in the linear interpolation   c
c smallest   : Index of the present smallest value, for bubble sort    c
c temp       : A temporary variable used for the bubble sort swap      c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer i, j, bottombound, topbound, smallest
      real temp

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c Bubble Sort algorithm, assuring that the elements are in order so    c
c that the interpolation is between the two closest points to the      c
c point in question.                                                   c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do i=1,nparam
        smallest = i
        do j=i,nparam
          if (xlist(j) .lt. xlist(smallest)) then
            smallest = j
          endif
        enddo

        temp = xlist(i)
        xlist(i) = xlist(smallest)
        xlist(smallest) = temp

        temp = tlist(i)
        tlist(i) = tlist(smallest)
        tlist(smallest) = temp
      enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c If the point in question is larger than the largest x-value or       c
c smaller than the smallest x-value, crash.                            c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if ((xvalue .lt. xlist(1)) .or. (xvalue .gt. xlist(nparam))) then
        print *,'Value out of range'
        stop
      endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c Find the correct top and bottom bounds, do the interpolation, return c
c the value.                                                           c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do i=2,nparam
        if ((xvalue .ge. xlist(i-1)) .and. (xvalue .le. xlist(i))) then
          bottombound = i-1
          topbound    = i
        endif
      enddo
                                                                                                               
      tvalue = (xvalue - xlist(bottombound)) /
     .         (xlist(topbound) - xlist(bottombound)) *
     .         (tlist(topbound) - tlist(bottombound)) +
     .         tlist(bottombound)
                                                                                                               
      return
      end
