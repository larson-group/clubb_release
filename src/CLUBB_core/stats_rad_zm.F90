!-----------------------------------------------------------------------
! $Id: stats_rad_zm.F90 4032 2009-08-17 21:45:29Z senkbeil@uwm.edu $

module stats_rad_zm

  implicit none

  private ! Default Scope

  public :: stats_init_rad_zm

! Constant parameters
  integer, parameter, public :: nvarmax_rad_zm = 250 ! Maximum variables allowed

  contains

!-----------------------------------------------------------------------
  subroutine stats_init_rad_zm( vars_rad_zm, l_error )

!     Description:
!     Initializes array indices for rad_zm variables
!-----------------------------------------------------------------------

    use constants, only:  &
        fstderr ! Constant(s)

    use stats_variables, only: & 
        rad_zm, &
        iFrad_LW_rad, & ! Variable(s)
        iFrad_SW_rad, &
        iFrad_SW_up_rad, &
        iFrad_LW_up_rad, &
        iFrad_SW_down_rad, &
        iFrad_LW_down_rad

    use stats_type, only: & 
        stat_assign ! Procedure

!use error_code, only: &
!    clubb_at_least_debug_level ! Function


    implicit none

    ! Input Variable
    character(len= * ), dimension(nvarmax_rad_zm), intent(in) :: vars_rad_zm

    ! Output Variable	
    logical, intent(inout) :: l_error

! Local Varables
    integer :: i, k

! Default initialization for array indices for rad_zm

    iFrad_LW_rad = 0
    iFrad_SW_rad = 0
    iFrad_SW_up_rad = 0
    iFrad_LW_up_rad = 0
    iFrad_SW_down_rad = 0
    iFrad_LW_down_rad = 0

!     Assign pointers for statistics variables rad_zm

    k = 1
    do i=1,rad_zm%nn

      select case ( trim(vars_rad_zm(i)) )
      
      case ('Frad_LW_rad')
        iFrad_LW_rad = k

        call stat_assign( iFrad_LW_rad, "Frad_LW_rad", & 
             "Net long-wave radiative flux [W/m^2]", "W/m^2", rad_zm )
        k = k + 1

      case ('Frad_SW_rad')
        iFrad_SW_rad = k

        call stat_assign( iFrad_SW_rad, "Frad_SW_rad", & 
             "Net short-wave radiative flux [W/m^2]", "W/m^2", rad_zm )
        k = k + 1

      case ('Frad_SW_up_rad')
        iFrad_SW_up_rad = k

        call stat_assign( iFrad_SW_up_rad, "Frad_SW_up_rad", & 
             "Short-wave upwelling radiative flux [W/m^2]", "W/m^2", rad_zm )
        k = k + 1

      case ('Frad_LW_up_rad')
        iFrad_LW_up_rad = k

        call stat_assign( iFrad_LW_up_rad, "Frad_LW_up_rad", & 
             "Long-wave upwelling radiative flux [W/m^2]", "W/m^2", rad_zm )
        k = k + 1

      case ('Frad_SW_down_rad')
        iFrad_SW_down_rad = k

        call stat_assign( iFrad_SW_down_rad, "Frad_SW_down_rad", & 
             "Short-wave downwelling radiative flux [W/m^2]", "W/m^2", rad_zm )
        k = k + 1

      case ('Frad_LW_down_rad')
        iFrad_LW_down_rad = k

        call stat_assign( iFrad_LW_down_rad, "Frad_LW_down_rad", & 
             "Long-wave downwelling radiative flux [W/m^2]", "W/m^2", rad_zm )
        k = k + 1        
        
      case default

        write(fstderr,*) 'Error:  unrecognized variable in vars_rad_zm:  ', trim( vars_rad_zm(i) )

        l_error = .true.  ! This will stop the run.


      end select

    end do

    return
  end subroutine stats_init_rad_zm

end module stats_rad_zm
