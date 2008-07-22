!-----------------------------------------------------------------------
! $Id: fill_holes_mult.F90,v 1.1 2008-07-22 16:04:24 faschinj Exp $
        module fill_holes_mult

        implicit none

        public :: fill_holes_multiplicative

        private :: vert_integrate

        private ! Default Scope

        contains

!-----------------------------------------------------------------------
        subroutine fill_holes_multiplicative & 
                  ( threshold, field_grid, field )

!         Description: This subroutine clips values of field that are below
!       threshold (i.e. "fills holes"), but conserves the total integrated 
!       mass of field.  This prevents clipping from acting as a spurious 
!       source.
!         Mass is conserved by reducing the clipped field everywhere by
!       a constant multiplicative coefficient.

!       References:
!          ``Numerical Methods for Wave Equations in Geophysical Fluid
!       Dynamics", Durran (1999), p. 292.  
!-----------------------------------------------------------------------

        use grid_class, only: & 
           gr ! Variable

        implicit none

        real, intent(in) :: & 
        threshold  ! A threshold (e.g. wtol*wtol) below which field must not 
                    ! fall                                       [Units vary]

        character(len=2), intent(in) :: & 
        field_grid ! The grid of the field, either zt or zm

        ! Input/Output variable
        real, dimension(gr%nnzp), intent(inout) ::  & 
        field  ! The field (e.g. wp2) that contains holes 
               !                                    [Units same as threshold]


        ! Local Variables

        real, dimension(gr%nnzp)  ::  & 
        field_clipped  ! The raw field (e.g. wp2) that contains no holes 
                       !                          [Units same as threshold]

        real ::  & 
        original_total_mass,  & ! Vertical integral of field [length * (units of field)]
        clipped_total_mass,   & ! Vertical integral of clipped field [length * (units of field)]
        mass_fraction  ! Coefficient that multiplies clipped field 
                       ! in order to conserve mass.                      []

!-----------------------------------------------------------------------

        ! Compute the field's vertical integral, which must be conserved.
        original_total_mass = vert_integrate( field, field_grid )

        ! Clip small or negative values from field.
        field_clipped = max( threshold, field )

        ! Compute the clipped field's vertical integral.
        ! clipped_total_mass >= original_total_mass
        clipped_total_mass = vert_integrate( field_clipped, field_grid )

        !   Compute coefficient that makes the clipped field have the same
        ! mass as the original field
        mass_fraction = ( original_total_mass - threshold ) / & 
                              ( clipped_total_mass  - threshold )

        ! If original_total_mass < 0, then set field = threshold everywhere.
        ! In this case, the hole filling acts like a net source.
        ! We want 0 <= multiplicative_coef <= 1
        mass_fraction = max( 0.0, mass_fraction )

        ! Output normalized, filled field
        field = mass_fraction * ( field_clipped - threshold )  & 
                     + threshold
        field = max( threshold, field )  ! Just in case there's round-off

        return
        end subroutine fill_holes_multiplicative



!-----------------------------------------------------------------------
        function vert_integrate( field, field_grid ) & 
        result( vert_integral )

!       Description:
!       Compute the vertical integral of a field.  

!       References:
!-----------------------------------------------------------------------

        use grid_class, only: & 
           gr ! Variable

        use error_code, only: & 
           clubb_debug  ! Subroutine

        implicit none

        ! Input variables
        real, dimension(gr%nnzp), intent(in) ::  & 
        field  ! The field (e.g. wp2) that contains holes  [Units vary]

        character(len=2), intent(in) :: & 
        field_grid ! The grid of the field, either zt or zm

        ! Output variable
        real :: & 
        vert_integral  ! Vertical integral of field    [length * (units of field) ]

        ! Local variables
        integer ::  & 
        k   ! Loop index

        ! If field is on the zt grid
        if ( field_grid == "zt" ) then

           vert_integral = 0.0
           do k = 2, gr%nnzp, 1
              vert_integral = vert_integral +  field(k)  / gr%dzt(k)
           end do

        elseif ( field_grid == "zm" )  then

           vert_integral = 0.0
           do k = 1, gr%nnzp-1, 1
              vert_integral = vert_integral +  field(k)  / gr%dzm(k)
           end do

        else

          call clubb_debug( 0,  & 
             "Neither zt nor zm grid is specified in vert_integrate" ) 
          vert_integral = -9.0e20

        end if

        return
        end function vert_integrate
!-----------------------------------------------------------------------

        end module fill_holes_mult
