!-------------------------------------------------------------------------------
!$Id$ 
!===============================================================================
module interpolation

  use clubb_precision, only: &
    core_rknd ! Variable(s)

  implicit none

  private ! Default Scope

  public :: lin_interpolate_two_points, binary_search, zlinterp, &
    lin_interpolate_on_grid_api, linear_interp_factor, mono_cubic_interp, &
    lin_interp_between_grids

  contains

!-------------------------------------------------------------------------------
  function lin_interpolate_two_points( height_int, height_high, height_low, &
    var_high, var_low )
  !$acc routine seq

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

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use constants_clubb, only: fstderr ! Constant

    implicit none

    !--------------------------- Input Variables ---------------------------

    real( kind = core_rknd ), intent(in) :: &
      height_int,  & ! Height to be interpolated to     [m]
      height_high, & ! Height above the interpolation   [m]
      height_low,  & ! Height below the interpolation   [m]
      var_high,    & ! Variable above the interpolation [units vary]
      var_low        ! Variable below the interpolation [units vary]

    !--------------------------- Output Variables ---------------------------
    real( kind = core_rknd ) :: lin_interpolate_two_points
    
    ! Check for valid input
#ifndef CLUBB_GPU
    if ( abs(height_low - height_high) < 1.0e-12_core_rknd ) then
      write(fstderr,*) "lin_interpolate_two_points: height_high and height_low cannot be equal."
      error stop
    end if
#endif

    ! Compute linear interpolation

    lin_interpolate_two_points = ( ( height_int - height_low )/( height_high - height_low ) ) &
      * ( var_high - var_low ) + var_low

    return
  end function lin_interpolate_two_points

  !-------------------------------------------------------------------------------------------------
  elemental real( kind = core_rknd ) function linear_interp_factor( factor, var_high, var_low )
  !$acc routine seq
  ! Description:
  !   Determines the coefficient for a linear interpolation
  ! 
  ! References:
  !   None
  !-------------------------------------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    !--------------------------- Input Variables ---------------------------
    real( kind = core_rknd ), intent(in) :: &
      factor,   & ! Factor                           [units vary]  
      var_high, & ! Variable above the interpolation [units vary]
      var_low     ! Variable below the interpolation [units vary]

    linear_interp_factor = factor * ( var_high - var_low ) + var_low

    return
  end function linear_interp_factor
  !-------------------------------------------------------------------------------------------------
  function mono_cubic_interp &
    ( z_in, km1, k00, kp1, kp2, zm1, z00, zp1, zp2, fm1, f00, fp1, fp2 ) result ( f_out )

  !$acc routine seq

  ! Description:
  !   Steffen's monotone cubic interpolation method
  !   Returns monotone cubic interpolated value between x00 and xp1

  ! Original Author:
  !   Takanobu Yamaguchi
  !   tak.yamaguchi@noaa.gov
  !
  !   This version has been modified slightly for CLUBB's coding standards and
  !   adds the 3/2 from eqn 21. -dschanen 26 Oct 2011
  !   We have also added a quintic polynomial option.
  !
  ! References:
  !   M. Steffen, Astron. Astrophys. 239, 443-450 (1990)
  !-------------------------------------------------------------------------------------------------

    use constants_clubb, only: &
        eps ! Constant

    use clubb_precision, only: &
        core_rknd ! Constant
    
    use model_flags, only: &
        l_quintic_poly_interp ! Variable(s)

    implicit none

    ! External
    intrinsic :: sign, abs, min

    !--------------------------- Input Variables ---------------------------
    real( kind = core_rknd ), intent(in) :: & 
      z_in ! The altitude to be interpolated to [m]

    ! k-levels;  their meaning depends on whether we're extrapolating or interpolating
    integer, intent(in) :: &
      km1, k00, kp1, kp2 

    real( kind = core_rknd ), intent(in) :: &
      zm1, z00, zp1, zp2, & ! The altitudes for km1, k00, kp1, kp2      [m]
      fm1, f00, fp1, fp2    ! The field at km1, k00, kp1, and kp2       [units vary]

    !--------------------------- Output Variables ---------------------------
    real( kind = core_rknd ) :: f_out ! The interpolated field
   
    ! Local Variables 
    real( kind = core_rknd ) :: &
      hm1, h00, hp1, &
      sm1, s00, sp1, &
      p00, pp1, &
      dfdx00, dfdxp1, &
      c1, c2, c3, c4, &
      w00, wp1, &
      coef1, coef2, &
      zprime, beta, alpha, zn
   
    ! ---- Begin Code ---- 

    coef1 = 1.0_core_rknd
    coef2 = 1.0_core_rknd

    if ( km1 <= k00 ) then
      hm1 = z00 - zm1
      h00 = zp1 - z00
      hp1 = zp2 - zp1

      if ( km1 == k00 ) then
        s00 = ( fp1 - f00 ) / ( zp1 - z00 )
        sp1 = ( fp2 - fp1 ) / ( zp2 - zp1 )
        dfdx00 = s00
        pp1 = ( s00 * hp1 + sp1 * h00 ) / ( h00 + hp1 )
        dfdxp1 = coef1*( sign( 1.0_core_rknd, s00 ) + sign( 1.0_core_rknd, sp1 ) ) &
          * min( abs( s00 ), abs( sp1 ), coef2*0.5_core_rknd*abs( pp1 ) )

      else if ( kp1 == kp2 ) then
        sm1 = ( f00 - fm1 ) / ( z00 - zm1 )
        s00 = ( fp1 - f00 ) / ( zp1 - z00 )
        p00 = ( sm1 * h00 + s00 * hm1 ) / ( hm1 + h00 )
        dfdx00 = coef1*( sign( 1.0_core_rknd, sm1 ) + sign( 1.0_core_rknd, s00 ) ) &
          * min( abs( sm1 ), abs( s00 ), coef2*0.5_core_rknd*abs( p00 ) )
        dfdxp1 = s00

      else
        sm1 = ( f00 - fm1 ) / ( z00 - zm1 )
        s00 = ( fp1 - f00 ) / ( zp1 - z00 )
        sp1 = ( fp2 - fp1 ) / ( zp2 - zp1 )
        p00 = ( sm1 * h00 + s00 * hm1 ) / ( hm1 + h00 )
        pp1 = ( s00 * hp1 + sp1 * h00 ) / ( h00 + hp1 )
        dfdx00 = coef1*( sign( 1.0_core_rknd, sm1 ) + sign( 1.0_core_rknd, s00 ) ) &
          * min( abs( sm1 ), abs( s00 ), coef2*0.5_core_rknd*abs( p00 ) )
        dfdxp1 = coef1*( sign( 1.0_core_rknd, s00 ) + sign( 1.0_core_rknd, sp1 ) ) &
          * min( abs( s00 ), abs( sp1 ), coef2*0.5_core_rknd*abs( pp1 ) )

      end if

      c1 = ( dfdx00 + dfdxp1 - 2._core_rknd * s00 ) / ( h00 ** 2 )
      c2 = ( 3._core_rknd * s00 - 2._core_rknd * dfdx00 - dfdxp1 ) / h00
      c3 = dfdx00
      c4 = f00

      if ( .not. l_quintic_poly_interp ) then

        ! Old formula
        !f_out = c1 * ( (z_in - z00)**3 ) + c2 * ( (z_in - z00)**2 ) + c3 * (z_in - z00) + c4

        ! Faster nested multiplication
        zprime = z_in - z00
        f_out =  c4 + zprime*( c3 + zprime*( c2 + ( zprime*c1 ) ) )

      else 

       ! Use a quintic polynomial interpolation instead instead of the Steffen formula.
       ! Unlike the formula above, this formula does not guarantee monotonicity.

        beta = 120._core_rknd * ( (fp1-f00) - 0.5_core_rknd * h00 * (dfdx00 + dfdxp1) )

        ! Prevent an underflow by using a linear interpolation
        if ( abs( beta ) < eps ) then 
          f_out = lin_interpolate_two_points( z00, zp1, zm1, &
                           fp1, fm1 )

        else
          alpha = (6._core_rknd/beta) * h00 * (dfdxp1-dfdx00) + 0.5_core_rknd
          zn = (z_in-z00)/h00

          f_out = ( &
                  (( (beta/20._core_rknd)*zn - (beta*(1._core_rknd+alpha) &
                    / 12._core_rknd)) * zn + (beta*alpha/6._core_rknd)) &
                    * zn**2 + dfdx00*h00 &
                  ) * zn + f00
        end if ! beta < eps
      end if ! ~quintic_polynomial

    else
      ! Linear extrapolation
      wp1 = ( z_in - z00 ) / ( zp1 - z00 )
      w00 = 1._core_rknd - wp1
      f_out = wp1 * fp1 + w00 * f00

    end if

    return
  end function mono_cubic_interp

!-------------------------------------------------------------------------------
  integer function binary_search( n, array, var ) & 
    result( i )

    ! Description:
    ! This subroutine performs a binary search to find the closest value greater
    ! than or equal to var in the array.  This function returns the index of the
    ! closest value of array that is greater than or equal to var.  It returns a
    ! value of -1 if var is outside the bounds of array.
    !
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd                  ! Variable(s)

    use constants_clubb, only: &
        fstderr                    ! Variable(s)
    
    implicit none

    ! Input Variables

    ! Size of the array
    integer, intent(in) :: n

    ! The array being searched (must be sorted from least value to greatest
    ! value).
    real( kind = core_rknd ), dimension(n), intent(in) :: array

    ! The value being searched for
    real( kind = core_rknd ), intent(in) :: var

    ! Bounds of the search
    integer :: high
    integer :: low

    ! The initial value of low has been changed from 1 to 2 due to a problem
    ! that was occuring when var was close to the lower bound.
    !
    ! The lowest value in the array (which is sorted by increasing values) is
    ! found at index 1, while the highest value in the array is found at
    ! index n.  Unless the value of var exactly corresponds with one of the
    ! values found in the array, or unless the value of var is found outside of
    ! the array, the value of var will be found between two levels of the array.
    ! In this scenario, the output of function binary_search is the index of the
    ! HIGHER level.  For example, if the value of var is found between array(1)
    ! and array(2), the output of function binary_search will be 2.
    !
    ! Therefore, the lowest index of a HIGHER level in an interpolation is 2.
    ! Thus, the initial value of low has been changed to 2.  This will prevent
    ! the value of variable "i" below from becoming 1.  If the value of "i"
    ! becomes 1, the code below tries to access array(0) (which is array(i-1)
    ! when i = 1) and produces an error.

    low = 2

    high = n

    ! Ensure variable is within range of array and array has more than 1 element
    if ( var < array(1) .or. var > array(n) .or. n < 2 ) then
        i = -1
        return
    end if

    ! Special case to cover var == array(1)
    if ( var >= array(1) .and. var <= array(2) ) then
        i = 2
        return
    end if


    do while( low <= high )

      i = (low + high) / 2

      if ( var > array(i-1) .and. var <= array(i) ) then

        return

      elseif ( var < array(i) ) then

        high = i - 1

      elseif ( var > array(i) ) then

        low = i + 1

      endif

    enddo 

    write(fstderr,*) "Logic error in function binary_search."

    i = -1
    return

  end function binary_search


!-------------------------------------------------------------------------------
  function zlinterp( ngrdcol, dim_out, & ! In
                     dim_src, grid_out, & ! In
                     grid_src, var_src, & ! In
                     nout ) result( var_out ) ! Optional In
  !$acc routine seq
! Description:
!   Interpolate every column independently in the vertical. Optional nout
!   limits the active output levels when merged grids have different sizes.
!   Assumes values that
!   are less than lowest source point are zero and above the highest
!   source point are zero. Also assumes altitude increases linearly.

! References:
!   function LIN_INT from WRF-HOC
!-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    !--------------------------- Input Variables ---------------------------
    integer, intent(in) :: ngrdcol ! Number of grid columns
    integer, intent(in) :: dim_out, dim_src

    real( kind = core_rknd ), dimension(ngrdcol,dim_src), intent(in) ::  &
      grid_src,  & ! [m]
      var_src      ! [units vary]

    real( kind = core_rknd ), dimension(ngrdcol,dim_out), intent(in) :: &
      grid_out ! [m]

    !--------------------------- Output Variables ---------------------------
    real( kind = core_rknd ), dimension(ngrdcol,dim_out) :: &
      var_out ! [units vary]

    !--------------------------- Optional Input Variables ---------------------------
    integer, intent(in), optional :: nout(ngrdcol) ! Active output levels in each column

    !--------------------------- Local Variables ---------------------------
    integer :: i, k, kint, km1 ! Column, source bracket and output level indices [-]
    integer :: last_out(ngrdcol) ! Last active output level in each column [-]

    !--------------------------- Begin Code ---------------------------
    var_out = 0.0_core_rknd
    last_out = dim_out
    if ( present(nout) ) last_out = nout

    do kint = 1, dim_out
      do i = 1, ngrdcol
        if ( kint > last_out(i) ) cycle

        ! Values below the source grid remain zero.
        if ( grid_out(i,kint) < grid_src(i,1) ) cycle

        k = binary_search( dim_src, grid_src(i,:), grid_out(i,kint) )

        ! Preserve the first above-source exit, even for a nonmonotone output grid.
        if ( k == -1 ) then
          last_out(i) = kint - 1
          cycle
        end if

        km1 = max( 1, k-1 )
        var_out(i,kint) = lin_interpolate_two_points( grid_out(i,kint), grid_src(i,k), &
          grid_src(i,km1), var_src(i,k), var_src(i,km1) )

      end do ! i
    end do ! kint

    return
  end function zlinterp


!-------------------------------------------------------------------------------
  subroutine lin_interpolate_on_grid_api &
             ( nparam, xlist, tlist, xvalue, &
               tvalue )

! Description:
!   Linear interpolation for 25 June 1996 altocumulus case.

!   For example, to interpolate between two temperatures in space, put
!   your spatial coordinates in x-list and your temperature values in
!   tlist.  The point in question should have its spatial value stored
!   in xvalue, and tvalue will be the temperature at that point.

! Author: Michael Falk for COAMPS.
!-------------------------------------------------------------------------------

    use constants_clubb, only: fstderr ! Constant

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use error_code, only: &
        clubb_at_least_debug_level_api ! Error indicator

    implicit none

    ! Input Variables
    integer, intent(in) :: nparam ! Number of parameters in xlist and tlist

    ! Input/Output Variables
    real( kind = core_rknd ), intent(inout), dimension(nparam) ::  & 
      xlist,  & ! List of x-values (independent variable)
      tlist     ! List of t-values (dependent variable)

    real( kind = core_rknd ), intent(in) ::  & 
      xvalue  ! x-value at which to interpolate

    real( kind = core_rknd ), intent(inout) ::  & 
      tvalue  ! t-value solved by interpolation

    ! Local variables
    integer ::  & 
      i     ! Loop control variable
 
    integer ::  & 
      bottombound, & ! Index of the smaller value in the linear interpolation
      topbound       ! Index of the larger value in the linear interpolation

!-------------------------------------------------------------------------------
!
! Assure that the elements are in order so that the interpolation is between 
! the two closest points to the point in question.
!
!-------------------------------------------------------------------------------

     if ( clubb_at_least_debug_level_api( 0 ) ) then
       do i=2,nparam
         if ( xlist(i) <= xlist(i-1) ) then
           write(fstderr,*) "xlist must be sorted for lin_interpolate_on_grid_api."
           error stop
         end if
       end do
     end if
!-------------------------------------------------------------------------------
!
! If the point in question is larger than the largest x-value or
! smaller than the smallest x-value, crash.
!
!-------------------------------------------------------------------------------

    if ( clubb_at_least_debug_level_api( 0 ) ) then
        if ( (xvalue < xlist(1)) .or. (xvalue > xlist(nparam)) ) then
          write(fstderr,*) "lin_interpolate_on_grid_api: Value out of range"
          error stop
        end if
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

    if ( clubb_at_least_debug_level_api( 1 ) .and. (topbound == -1 .or. bottombound == -1) ) then
      write(fstderr,*) "Sanity check failed! xlist is not properly sorted"
      write(fstderr,*) "in lin_interpolate_on_grid_api."
    end if

    tvalue = &
    lin_interpolate_two_points( xvalue, xlist(topbound), xlist(bottombound), &
            tlist(topbound), tlist(bottombound) )

    return
  end subroutine lin_interpolate_on_grid_api

  function lin_interp_between_grids( ngrdcol, size_interpolate, size_current, &
                                     interpolate_altitudes, current_altitudes, &
                                     current_values )
    ! Description:
    !   Interpolate all columns to the target grids, using the nearest source
    !   value outside each source grid and preserving the near-level tolerance.

    use constants_clubb, only: fstderr ! Constant

    use error_code, only: &
        clubb_at_least_debug_level_api ! Error indicator

    implicit none
    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
        ngrdcol, & ! Number of grid columns [-]
        size_interpolate, & ! Number of target levels [-]
        size_current ! Number of source levels [-]

    real( kind = core_rknd ), dimension(ngrdcol,size_interpolate), intent(in) :: &
      interpolate_altitudes ! Target heights [m]

    real( kind = core_rknd ), dimension(ngrdcol,size_current), intent(in) :: &
      current_altitudes ! Source heights [m]

    real( kind = core_rknd ), dimension(ngrdcol,size_current), intent(in) :: &
      current_values ! given values on the current_altitudes altitudes that should be interpolated 
                     ! to the interpolate_altitudes altitudes

    !--------------------- Return Variable ---------------------
    real( kind = core_rknd ), dimension(ngrdcol,size_interpolate) :: &
      lin_interp_between_grids ! Field on target heights [units vary]

    !--------------------- Local Variables ---------------------
    
    integer :: i, j, k ! Output level, column and search indices [-]
    
    logical :: calc_done ! Whether a source bracket has been found [-]

    real( kind = core_rknd ) :: tol ! Tolerance for coincident grid levels [m]

    !--------------------- Begin Code ---------------------
    tol = 1.0e-6_core_rknd


    do j = 1, ngrdcol

      !-------------------------------------------------------------------------------
      !
      ! Assure that the elements are in order so that the interpolation is between
      ! the two closest points to the point in question.
      !
      !-------------------------------------------------------------------------------

      if ( clubb_at_least_debug_level_api( 0 ) ) then
        do i=2,size_current
          if ( current_altitudes(j,i) <= current_altitudes(j,i-1) ) then
            write(fstderr,*) "current_altitudes must be sorted for lin_interp_between_grids."
            error stop
          end if
        end do
      end if

      !-------------------------------------------------------------------------------
      !
      ! If the point in question is larger than the largest x-value or
      ! smaller than the smallest x-value, give a warning.
      !
      !-------------------------------------------------------------------------------

      if ( clubb_at_least_debug_level_api( 0 ) ) then
          if ( (minval(interpolate_altitudes(j,:)) < current_altitudes(j,1)) &
                .or. (maxval(interpolate_altitudes(j,:)) > current_altitudes(j,size_current)) ) then
            write(fstderr,*) "lin_interp_between_grids: Value is out of range. Extrapolating..."
          end if
      end if


    end do

    do i = 1, size_interpolate
      do j = 1, ngrdcol

        calc_done = .false.
        ! Use a nearby source level, or interpolate between its neighbors.
        k = 1
        do while ( (.not. calc_done) .and. (k <= size_current) )

          if ( abs(interpolate_altitudes(j,i) - current_altitudes(j,k)) < tol ) then
            lin_interp_between_grids(j,i) = current_values(j,k)
            calc_done = .true.
          else if ( interpolate_altitudes(j,i) < current_altitudes(j,k) ) then

            if ( k > 1 ) then
              lin_interp_between_grids(j,i) = lin_interpolate_two_points( & ! In
                interpolate_altitudes(j,i), current_altitudes(j,k), & ! In
                current_altitudes(j,k-1), current_values(j,k), current_values(j,k-1) ) ! In
            else
              lin_interp_between_grids(j,i) = current_values(j,1)
            end if

            calc_done = .true.

          end if

          k = k + 1

        end do

        if ( (.not. calc_done) .and. (k > size_current) ) then
          lin_interp_between_grids(j,i) = current_values(j,size_current)
          calc_done = .true.
        end if

      end do
    end do

  end function lin_interp_between_grids

end module interpolation
