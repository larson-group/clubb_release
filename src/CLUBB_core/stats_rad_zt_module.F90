!-----------------------------------------------------------------------
! $Id$
!===============================================================================

module stats_rad_zt_module

  implicit none

  private ! Default Scope

  public :: stats_init_rad_zt

  ! Constant parameters
  integer, parameter, public :: nvarmax_rad_zt = 250 ! Maximum variables allowed

  contains

!-----------------------------------------------------------------------
  subroutine stats_init_rad_zt( vars_rad_zt,                    & ! intent(in)
                                l_error,                        & ! intent(inout)
                                stats_metadata, stats_rad_zt )    ! intent(inout)

! Description:
!   Initializes array indices for stats_zt
!
! References:
!   None
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

    !--------------------- Input Variable ---------------------
    character(len= * ), dimension(nvarmax_rad_zt), intent(in) :: &
      vars_rad_zt

    !--------------------- InOut Variables ---------------------      
    type (stats_metadata_type), intent(inout) :: &
      stats_metadata

    type (stats), target, intent(inout) :: &
      stats_rad_zt

    logical, intent(inout) :: l_error

    !--------------------- Local Varables ---------------------
    integer :: i, k

    !--------------------- Begin Code ---------------------

    ! Default initialization for array indices for stats_rad_zt

    stats_metadata%iT_in_K_rad             = 0
    stats_metadata%ircil_rad               = 0
    stats_metadata%io3l_rad                = 0
    stats_metadata%irsm_rad                = 0
    stats_metadata%ircm_in_cloud_rad       = 0
    stats_metadata%icloud_frac_rad         = 0
    stats_metadata%iice_supersat_frac_rad  = 0
    stats_metadata%iradht_rad              = 0
    stats_metadata%iradht_LW_rad           = 0
    stats_metadata%iradht_SW_rad           = 0
    stats_metadata%ip_in_mb_rad            = 0
    stats_metadata%isp_humidity_rad        = 0


    ! Assign pointers for statistics variables stats_rad_zt

    k = 1
    do i=1,stats_rad_zt%num_output_fields

      select case ( trim(vars_rad_zt(i)) )

      case ('T_in_K_rad')
        stats_metadata%iT_in_K_rad = k

        call stat_assign( var_index=stats_metadata%iT_in_K_rad, var_name="T_in_K_rad", & ! intent(in)
             var_description="Temperature [K]", var_units="K", l_silhs=.false., & ! intent(in)
                grid_kind=stats_rad_zt )! intent(inout)
        k = k + 1

      case ('rcil_rad')
        stats_metadata%ircil_rad = k

        call stat_assign( var_index=stats_metadata%ircil_rad, var_name="rcil_rad", & ! intent(in)
             var_description="Ice mixing ratio [kg/kg]", var_units="kg/kg", & ! intent(in)
             l_silhs=.false., & ! intent(in)
             grid_kind=stats_rad_zt ) ! intent(inout)
        k = k + 1

      case ('o3l_rad')
        stats_metadata%io3l_rad = k

        call stat_assign( var_index=stats_metadata%io3l_rad, var_name="o3l_rad", & ! intent(in)
             var_description="Ozone mixing ratio [kg/kg]", var_units="kg/kg", & ! intent(in)
             l_silhs=.false., & ! intent(in)
             grid_kind=stats_rad_zt ) ! intent(inout)
        k = k + 1

      case ('rsm_rad')
        stats_metadata%irsm_rad = k

        call stat_assign( var_index=stats_metadata%irsm_rad, var_name="rsm_rad", & ! intent(in)
             var_description="Snow water mixing ratio [kg/kg]", var_units="kg/kg", & ! intent(in)
             l_silhs=.false., & ! intent(in)
             grid_kind=stats_rad_zt ) ! intent(inout)
        k = k + 1

      case ('rcm_in_cloud_rad')
        stats_metadata%ircm_in_cloud_rad = k

        call stat_assign( var_index=stats_metadata%ircm_in_cloud_rad, var_name="rcm_in_cloud_rad", & ! intent(in)
             var_description="rcm in cloud layer [kg/kg]", var_units="kg/kg", & ! intent(in)
             l_silhs=.false., & ! intent(in)
             grid_kind=stats_rad_zt ) ! intent(inout)
        k = k + 1

      case ('cloud_frac_rad')
        stats_metadata%icloud_frac_rad = k

        call stat_assign( var_index=stats_metadata%icloud_frac_rad, var_name="cloud_frac_rad", & ! intent(in)
             var_description="Cloud fraction (between 0 and 1) [-]", & ! intent(in)
             var_units="count", & ! intent(in)
             l_silhs=.false., & ! intent(in)
             grid_kind=stats_rad_zt ) ! intent(inout)
        k = k + 1
      
      case ('ice_supersat_frac_rad')
        stats_metadata%iice_supersat_frac_rad = k

        call stat_assign( var_index=stats_metadata%iice_supersat_frac_rad, var_name="ice_supersat_frac_rad", & !In
             var_description="Ice cloud fraction (between 0 and 1) [-]", var_units="count", & ! In
             l_silhs=.false., & ! intent(in)
             grid_kind=stats_rad_zt ) ! intent(inout)
        k = k + 1

      case ('radht_rad')
        stats_metadata%iradht_rad = k

        call stat_assign( var_index=stats_metadata%iradht_rad, var_name="radht_rad", & ! intent(in)
             var_description="Total radiative heating rate [K/s]", var_units="K/s", & ! intent(in)
             l_silhs=.false., & ! intent(in)
             grid_kind=stats_rad_zt ) ! intent(inout)
        k = k + 1

      case ('radht_LW_rad')
        stats_metadata%iradht_LW_rad = k

        call stat_assign( var_index=stats_metadata%iradht_LW_rad, var_name="radht_LW_rad", & ! intent(in)
             var_description="Long-wave radiative heating rate [K/s]", var_units="K/s", & ! In
             l_silhs=.false., & ! intent(in)
             grid_kind=stats_rad_zt ) ! intent(inout)
        k = k + 1

      case ('radht_SW_rad')
        stats_metadata%iradht_SW_rad = k

        call stat_assign( var_index=stats_metadata%iradht_SW_rad, var_name="radht_SW_rad", & ! intent(in)
             var_description="Short-wave radiative heating rate [K/s]", var_units="K/s", & ! In
             l_silhs=.false., & ! intent(in)
             grid_kind=stats_rad_zt ) ! intent(inout)
        k = k + 1

      case ('p_in_mb_rad')
        stats_metadata%ip_in_mb_rad = k

        call stat_assign( var_index=stats_metadata%ip_in_mb_rad, var_name="p_in_mb_rad", & ! intent(in)
             var_description="Pressure [hPa]", var_units="hPa", & ! intent(in)
             l_silhs=.false., & ! intent(in)
             grid_kind=stats_rad_zt ) ! intent(inout)
        k = k + 1

      case ('sp_humidity_rad')
        stats_metadata%isp_humidity_rad = k

        call stat_assign( var_index=stats_metadata%isp_humidity_rad, var_name="sp_humidity_rad", & ! intent(in)
             var_description="Specific humidity [kg/kg]", var_units="kg/kg", & ! intent(in)
             l_silhs=.false., & ! intent(in)
             grid_kind=stats_rad_zt ) ! intent(inout)
        k = k + 1

      case default

        write(fstderr,*) 'Error:  unrecognized variable in vars_rad_zt:  ', trim( vars_rad_zt(i) )

        l_error = .true.  ! This will stop the run.


      end select

    end do

    return
  end subroutine stats_init_rad_zt

end module stats_rad_zt_module
