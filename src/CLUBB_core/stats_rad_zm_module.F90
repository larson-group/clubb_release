!-----------------------------------------------------------------------
! $Id$
!===============================================================================

module stats_rad_zm_module

  implicit none

  private ! Default Scope

  public :: stats_init_rad_zm

! Constant parameters
  integer, parameter, public :: nvarmax_rad_zm = 250 ! Maximum variables allowed

  contains

!-----------------------------------------------------------------------
  subroutine stats_init_rad_zm( vars_rad_zm,                    & ! In
                                l_error,                        & ! In
                                stats_metadata, stats_rad_zm )    ! In/Out

!     Description:
!     Initializes array indices for stats_rad_zm variables
!-----------------------------------------------------------------------

    use constants_clubb, only:  &
        fstderr ! Constant(s)

    use stats_type_utilities, only: & 
        stat_assign ! Procedure

    use stats_type, only: &
        stats ! Type

    use stats_variables, only: &
        stats_metadata_type

    implicit none

    !------------------------ Input Variable ------------------------
    character(len= * ), dimension(nvarmax_rad_zm), intent(in) :: vars_rad_zm

    !--------------------- InOut Variables ---------------------      
    type (stats), intent(inout) :: &
      stats_rad_zm

    type (stats_metadata_type), intent(inout) :: &
      stats_metadata

    logical, intent(inout) :: l_error

    !------------------------ Local Varables ------------------------
    integer :: i, k

    !------------------------ Begin Code ------------------------

    ! Default initialization for array indices for stats_rad_zm

    stats_metadata%iFrad_LW_rad = 0
    stats_metadata%iFrad_SW_rad = 0
    stats_metadata%iFrad_SW_up_rad = 0
    stats_metadata%iFrad_LW_up_rad = 0
    stats_metadata%iFrad_SW_down_rad = 0
    stats_metadata%iFrad_LW_down_rad = 0

    stats_metadata%ifulwcl = 0
    stats_metadata%ifdlwcl = 0
    stats_metadata%ifdswcl = 0
    stats_metadata%ifuswcl = 0

!     Assign pointers for statistics variables stats_rad_zm

    k = 1
    do i=1,stats_rad_zm%num_output_fields

      select case ( trim(vars_rad_zm(i)) )

      case('fulwcl')
        stats_metadata%ifulwcl = k
        call stat_assign( var_index=stats_metadata%ifulwcl, & ! In
             var_name="fulwcl", & ! In
             var_description="Upward clear-sky LW flux [W/m^2]", var_units="W/m^2", & ! In
             l_silhs=.false., & ! In
             grid_kind=stats_rad_zm ) ! In/Out
        k = k + 1

      case( 'fdlwcl' )
        stats_metadata%ifdlwcl = k
        call stat_assign( var_index=stats_metadata%ifdlwcl, & ! In
             var_name="fdlwcl", & ! In
             var_description="Downward clear-sky LW flux [W/m^2]", var_units="W/m^2", & ! In
             l_silhs=.false., & ! In
             grid_kind=stats_rad_zm ) ! In/Out
        k = k + 1

      case( 'fdswcl' )
        stats_metadata%ifdswcl = k
        call stat_assign( var_index=stats_metadata%ifdswcl, & ! In
             var_name="fdswcl", & ! In
             var_description="Downward clear-sky SW flux [W/m^2]", var_units="W/m^2", & ! In
             l_silhs=.false., & ! In
             grid_kind=stats_rad_zm ) ! In/Out
        k = k + 1

      case( 'fuswcl' )
        stats_metadata%ifuswcl = k
        call stat_assign( var_index=stats_metadata%ifuswcl,& ! In
             var_name="fuswcl", & ! In
             var_description="Upward clear-sky SW flux [W/m^2]", var_units="W/m^2", & ! In
             l_silhs=.false., & ! In
             grid_kind=stats_rad_zm ) ! In/Out
        k = k + 1

      case ('Frad_LW_rad')
        stats_metadata%iFrad_LW_rad = k

        call stat_assign( var_index=stats_metadata%iFrad_LW_rad, & ! In
             var_name="Frad_LW_rad", & ! In
             var_description="Net long-wave radiative flux [W/m^2]", var_units="W/m^2", & ! In
             l_silhs=.false., & ! In
             grid_kind=stats_rad_zm ) ! In/Out
        k = k + 1

      case ('Frad_SW_rad')
        stats_metadata%iFrad_SW_rad = k

        call stat_assign( var_index=stats_metadata%iFrad_SW_rad, & ! In
             var_name="Frad_SW_rad", & ! In
             var_description="Net short-wave radiative flux [W/m^2]", var_units="W/m^2", & ! In
             l_silhs=.false., & ! In
             grid_kind=stats_rad_zm ) ! In/Out
        k = k + 1

      case ('Frad_SW_up_rad')
        stats_metadata%iFrad_SW_up_rad = k

        call stat_assign( var_index=stats_metadata%iFrad_SW_up_rad, & ! In
             var_name="Frad_SW_up_rad", & ! In
             var_description="Short-wave upwelling radiative flux [W/m^2]", var_units="W/m^2", &!In
             l_silhs=.false., & ! In
             grid_kind=stats_rad_zm ) ! In/Out
        k = k + 1

      case ('Frad_LW_up_rad')
        stats_metadata%iFrad_LW_up_rad = k

        call stat_assign( var_index=stats_metadata%iFrad_LW_up_rad, & ! In
             var_name="Frad_LW_up_rad", & ! In
             var_description="Long-wave upwelling radiative flux [W/m^2]", var_units="W/m^2", &!In
             l_silhs=.false., & ! In
             grid_kind=stats_rad_zm ) ! In/Out
        k = k + 1

      case ('Frad_SW_down_rad')
        stats_metadata%iFrad_SW_down_rad = k
 
        call stat_assign( var_index=stats_metadata%iFrad_SW_down_rad, & ! In
             var_name="Frad_SW_down_rad", & ! In
             var_description="Short-wave downwelling radiative flux [W/m^2]", & ! In
             var_units="W/m^2", & ! In
             l_silhs=.false., & ! In
             grid_kind=stats_rad_zm ) ! In/Out
        k = k + 1

      case ('Frad_LW_down_rad')
        stats_metadata%iFrad_LW_down_rad = k

        call stat_assign( var_index=stats_metadata%iFrad_LW_down_rad, & ! In
             var_name="Frad_LW_down_rad", & ! In
             var_description="Long-wave downwelling radiative flux [W/m^2]", & ! In
             var_units="W/m^2", & ! In
             l_silhs=.false., & ! In
             grid_kind=stats_rad_zm ) ! In/Out
        k = k + 1

      case default

        write(fstderr,*) 'Error:  unrecognized variable in vars_rad_zm:  ', trim( vars_rad_zm(i) )

        l_error = .true.  ! This will stop the run.


      end select

    end do

    return
  end subroutine stats_init_rad_zm

end module stats_rad_zm_module
