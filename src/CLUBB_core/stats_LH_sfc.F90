!-----------------------------------------------------------------------
! $Id$

module stats_LH_sfc


  implicit none

  private ! Set Default Scope

  public :: stats_init_LH_sfc

  ! Constant parameters
  integer, parameter, public :: nvarmax_LH_sfc = 10  ! Maximum variables allowed

  contains

!-----------------------------------------------------------------------
  subroutine stats_init_LH_sfc( vars_LH_sfc, l_error )

! Description:
!   Initializes array indices for LH_sfc
! References:
!   None
!-----------------------------------------------------------------------

    use constants_clubb, only: &
      fstderr ! Constant(s)

    use stats_variables, only: & 
      LH_sfc ! Variable(s)

    use stats_variables, only: & 
      iLH_morr_snow_rate, & ! Variable(s)
      iLH_vwp, &
      iLH_lwp
      
    use stats_type, only: & 
        stat_assign ! Procedure

    implicit none

    ! External
    intrinsic :: trim

    ! Input Variable
    character(len= * ), dimension(nvarmax_LH_sfc), intent(in) :: vars_LH_sfc

    ! Input / Output Variable        
    logical, intent(inout) :: l_error

    ! Local Varables
    integer :: i, k

    ! ---- Begin Code ----

    ! Default initialization for array indices for sfc is zero (see module
    ! stats_variables)

    ! Assign pointers for statistics variables sfc

    k = 1
    do i = 1, LH_sfc%nn

      select case ( trim( vars_LH_sfc(i) ) )

      case ( 'LH_morr_snow_rate' )
        iLH_morr_snow_rate = k
        call stat_assign( var_index=iLH_morr_snow_rate, var_name="LH_morr_snow_rate", &
             var_description="Snow+Ice+Graupel fallout rate from Morrison scheme [mm/day]", &
             var_units="mm/day", l_silhs=.true., grid_kind=LH_sfc )
        k = k + 1

      case ( 'LH_vwp' )
        iLH_vwp = k
        call stat_assign( var_index=iLH_vwp, var_name="LH_vwp", &
             var_description="Vapor water path [kg/m^2]", var_units="kg/m^2", l_silhs=.true., &
             grid_kind=LH_sfc )
        k = k + 1

      case ( 'LH_lwp' )
        iLH_lwp = k
        call stat_assign( var_index=iLH_lwp, var_name="LH_lwp", &
             var_description="Liquid water path [kg/m^2]", var_units="kg/m^2", l_silhs=.true., &
             grid_kind=LH_sfc )
        k = k + 1

      case default
        write(fstderr,*) 'Error:  unrecognized variable in vars_LH_sfc:  ',  &
              trim( vars_LH_sfc(i) )
        l_error = .true.  ! This will stop the run.

      end select

    end do ! i = 1, LH_sfc%nn

    return
  end subroutine stats_init_LH_sfc

end module stats_LH_sfc

