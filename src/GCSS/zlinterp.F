        module zlinterp
        implicit none

        public :: zlinterp_fnc
        contains
!-----------------------------------------------------------------------
        pure function zlinterp_fnc( dim_out, dim_src, grid_out, 
     .                              grid_src, var_src ) 
     .  result( var_out )
!       Description:
!       Do a linear interpolation in the vertical.  Assumes values that
!       are less than lowest source point are zero and above the highest
!       source point are zero. Also assumes altitude increases linearly.

!       References:
!       function LIN_INT from WRF-HOC
!-----------------------------------------------------------------------

        implicit none

        ! Input variables
        integer, intent(in) :: dim_out, dim_src

        real, dimension(dim_src), intent(in) :: 
     .  grid_src, ! [m]
     .  var_src   ! [generic units]

        real, dimension(dim_out), intent(in) :: grid_out ! [m]

        ! Output variable
        real, dimension(dim_out) :: var_out ! [generic units]

        ! Local variables
        integer :: k, kint, km1, kp1

        k = 1

        do kint = 1, dim_out, 1

          ! Set to 0 if we're below the input data's lowest point
          if ( grid_out(kint) < grid_src(1) ) then
            var_out(kint) = 0.0
            cycle
          end if

          ! Increment k until the level is correct
          do while ( grid_out(kint) > grid_src(k) 
     .                .and. k < dim_src )
            k = k + 1
          end do

          ! If the increment leads to a level above the data, set this
          ! point and all those above it to zero
          if ( k > dim_src ) then
            var_out(kint:dim_out) = 0.0
            exit
          end if

          km1 = max( 1, k-1 )
          !kp1 = min( k+1, dim_src )

          ! Interpolate
          var_out(kint) = 
     .    ( var_src(k) - var_src(km1) ) / 
     .    ( grid_src(k) - grid_src(km1) )
     .      * ( grid_out(kint) - grid_src(km1) ) + var_src(km1)

        end do

        return
        end function zlinterp_fnc

        end module zlinterp
