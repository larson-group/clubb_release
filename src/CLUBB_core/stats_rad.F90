!-----------------------------------------------------------------------
! $Id: stats_rad.F90 4032 2009-08-17 21:45:29Z senkbeil@uwm.edu $

module stats_rad

  implicit none

  private ! Default Scope

  public :: stats_init_rad

! Constant parameters
  integer, parameter, public :: nvarmax_rad = 250 ! Maximum variables allowed

  contains

!-----------------------------------------------------------------------
  subroutine stats_init_rad( vars_rad, l_error )

!     Description:
!     Initializes array indices for zt
!-----------------------------------------------------------------------

    use constants, only:  &
        fstderr ! Constant(s)

    use stats_variables, only: & 
        rad, &
        iT_in_K_rad, & ! Variable(s)
        ircil_rad, &
        io3l_rad, &
        irsnowm_rad, &
        ircm_in_cloud_rad, &
        icloud_frac_rad, & 
        iradht_rad, &
        iradht_LW_rad, &
        iradht_SW_rad, &
        iFrad_LW_rad, &
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
    character(len= * ), dimension(nvarmax_rad), intent(in) :: vars_rad

    ! Output Variable	
    logical, intent(inout) :: l_error

! Local Varables
    integer :: i, k

! Default initialization for array indices for rad

    iT_in_K_rad = 0
    ircil_rad = 0
    io3l_rad = 0
    irsnowm_rad = 0
    ircm_in_cloud_rad = 0
    icloud_frac_rad = 0
    iradht_rad = 0
    iradht_LW_rad = 0
    iradht_SW_rad = 0
    iFrad_LW_rad = 0
    iFrad_SW_rad = 0
    iFrad_SW_up_rad = 0
    iFrad_LW_up_rad = 0
    iFrad_SW_down_rad = 0
    iFrad_LW_down_rad = 0

!     Assign pointers for statistics variables rad

    k = 1
    do i=1,rad%nn

      select case ( trim(vars_rad(i)) )
      
      case ('T_in_K_rad')
        iT_in_K_rad = k

        call stat_assign( iT_in_K_rad, "T_in_K_rad", & 
             "Temperature [K]", "K", rad )
        k = k + 1
     
      case ('rcil_rad')
        ircil_rad = k

        call stat_assign( ircil_rad, "rcil_rad", & 
             "Ice mixing ratio [kg/kg]", "kg/kg", rad )
        k = k + 1

      case ('o3l_rad')
        io3l_rad = k

        call stat_assign( io3l_rad, "o3l_rad", & 
             "Ozone mixing ratio [kg/kg]", "kg/kg", rad )
        k = k + 1

      case ('rsnowm_rad')
        irsnowm_rad = k

        call stat_assign( irsnowm_rad, "rsnowm_rad", & 
             "Snow water mixing ratio [kg/kg]", "kg/kg", rad )
        k = k + 1

      case ('rcm_in_cloud_rad')
        ircm_in_cloud_rad = k

        call stat_assign( ircm_in_cloud_rad, "rcm_in_cloud_rad", & 
             "rcm in cloud layer [kg/kg]", "kg/kg", rad )
        k = k + 1        

      case ('cloud_frac_rad')
        icloud_frac_rad = k

        call stat_assign( icloud_frac_rad, "cloud_frac_rad", & 
             "Cloud fraction (between 0 and 1) [-]", "count", rad )
        k = k + 1

      case ('radht_rad')
        iradht_rad = k

        call stat_assign( iradht_rad, "radht_rad", & 
             "Total radiative heating rate [K/s]", "K/s", rad )
        k = k + 1

      case ('radht_LW_rad')
        iradht_LW_rad = k

        call stat_assign( iradht_LW_rad, "radht_LW_rad", & 
             "Long-wave radiative heating rate [K/s]", "K/s", rad )
        k = k + 1

      case ('radht_SW_rad')
        iradht_SW_rad = k

        call stat_assign( iradht_SW_rad, "radht_SW_rad", & 
             "Short-wave radiative heating rate [K/s]", "K/s", rad )
        k = k + 1

      case ('Frad_LW_rad')
        iFrad_LW_rad = k

        call stat_assign( iFrad_LW_rad, "Frad_LW_rad", & 
             "Net long-wave radiative flux [W/m^2]", "W/m^2", rad )
        k = k + 1

      case ('Frad_SW_rad')
        iFrad_SW_rad = k

        call stat_assign( iFrad_SW_rad, "Frad_SW_rad", & 
             "Net short-wave radiative flux [W/m^2]", "W/m^2", rad )
        k = k + 1

      case ('Frad_SW_up_rad')
        iFrad_SW_up_rad = k

        call stat_assign( iFrad_SW_up_rad, "Frad_SW_up_rad", & 
             "Short-wave upwelling radiative flux [W/m^2]", "W/m^2", rad )
        k = k + 1

      case ('Frad_LW_up_rad')
        iFrad_LW_up_rad = k

        call stat_assign( iFrad_LW_up_rad, "Frad_LW_up_rad", & 
             "Long-wave upwelling radiative flux [W/m^2]", "W/m^2", rad )
        k = k + 1

      case ('Frad_SW_down_rad')
        iFrad_SW_down_rad = k

        call stat_assign( iFrad_SW_down_rad, "Frad_SW_down_rad", & 
             "Short-wave downwelling radiative flux [W/m^2]", "W/m^2", rad )
        k = k + 1

      case ('Frad_LW_down_rad')
        iFrad_LW_down_rad = k

        call stat_assign( iFrad_LW_down_rad, "Frad_LW_down_rad", & 
             "Long-wave downwelling radiative flux [W/m^2]", "W/m^2", rad )
        k = k + 1        
        
      case default

        write(fstderr,*) 'Error:  unrecognized variable in vars_rad:  ', trim( vars_rad(i) )

        l_error = .true.  ! This will stop the run.


      end select

    end do

    return
  end subroutine stats_init_rad

end module stats_rad
