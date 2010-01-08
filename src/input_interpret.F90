!$Id$
module input_interpret

  implicit none  

  public :: &
    read_z_profile, &
    read_theta_profile, &
    read_subs_profile
  private ! Default scope

  contains

  !------------------------------------------------------------------------------
  subroutine read_z_profile( nvar, nsize, retVars, psfc, zm_init, z, p_in_Pa, alt_type)
    !
    !  Description: Searches for the variable specified by either 'z[m]' or
    !  'Press[Pa]' in the collection of retVars. If the subroutine finds the
    !  variable indicated by 'z[m]',
    !  then it returns it. If the subroutine finds 'Press[Pa]' then it converts
    !  it to values of altitude in meters.
    !  If it does not find either or finds both the program using this subroutine
    !  will exit gracefully with a warning message.
    !
    !-------------------------------------------------------------------------------

    use input_reader, only: read_x_profile, & ! Prodedure(s)
      one_dim_read_var ! Type

    use constants, only: kappa, p0, Cp, Lv, zero_threshold, ep2, ep1 ! Variable(s)

    use saturation, only: sat_mixrat_liq, sat_rcm ! Procedure(s)

    use parameters_model, only: T0 ! Variable(s)

    use hydrostatic_mod, only: inverse_hydrostatic ! Procedure(s)

    use input_names, only: &
      z_name, &
      pressure_name, &
      rt_name, &
      temperature_name, &
      thetal_name, &
      theta_name

    implicit none

    ! Input Variable(s)
    integer, intent(in) :: nvar ! Number of elements in retVars

    integer, intent(in) :: nsize

    type(one_dim_read_var), dimension(nvar), intent(in) :: retVars ! Collection
    !                                                                being searched

    real, intent(in) :: &
      psfc, &         ! Pressure at the surface [Pa]
      zm_init         ! Height at zm(1)         [m]

    ! Output Variable(s)

    real, intent(out), dimension(nsize) :: z ! Height sounding profile [m]

    real, intent(out), dimension(nsize) :: p_in_Pa ! Pressure sounding profile [Pa]

    character(len=*), intent(out) :: alt_type ! Indicates where altitudes were
    !                                           gained from

    intrinsic :: max

    ! Local Variables
    real, dimension( nsize ) :: exner, thvm, rcm, theta, rtm

    integer :: nlevels, k

    character(len=40) :: theta_type

    if( count( (/ any(retVars%name == z_name), any(retVars%name == pressure_name) /)) <= 1) then
      if( any(retVars%name == z_name))then
        alt_type = z_name
        z = read_x_profile( nvar, nsize, alt_type, retVars )
        p_in_Pa = -999.9

      elseif( any(retVars%name == pressure_name))then
        alt_type = pressure_name

        p_in_Pa = read_x_profile( nvar, nsize, alt_type, retVars )

        nlevels = size(retVars(1)%values)

        call read_theta_profile(nvar, nsize, retVars, theta_type, theta )

        rtm = read_x_profile(nvar, nsize, rt_name, retVars)

        do k = 1, nlevels
          exner(k) = (p_in_Pa(k)/p0) ** kappa  ! zt
        end do

        if( trim( theta_type ) == temperature_name ) then
          theta = theta / exner
          theta_type = theta_name
        end if

        do k = 1,nlevels
          rcm(k) = &
          max( rtm(k) - sat_mixrat_liq( p_in_Pa(k), theta(k) * exner(k) ), &
            zero_threshold )
        enddo

        ! Compute initial theta-l

        select case ( trim( theta_type ) )
        case ( thetal_name )
          !case ( "dycoms2_rf01", "astex_a209", "nov11_altocu", &
          !      "clex9_nov02", "clex9_oct14", "dycoms2_rf02" )
          ! thlm profile that is initially saturated at points.
          ! thlm profile remains the same as in the input sounding.
          ! use iterative method to find initial rcm.
          do k =1, nlevels, 1
            rcm(k) = sat_rcm( theta(k), rtm(k), p_in_Pa(k), exner(k) )
          end do

        case default ! theta_name
          ! Initial profile is non-saturated thlm or any type of theta.
          theta(1:nlevels) = theta(1:nlevels) &
                           - Lv/(Cp*exner(1:nlevels)) * rcm(1:nlevels)

        end select

        ! Now, compute initial thetav
        do k = 1, nlevels, 1
          thvm(k) = theta(k) + ep1 * T0 * rtm(k)  & 
               + ( Lv/(Cp*exner(k)) - ep2 * T0 ) * rcm(k)
        end do

        call inverse_hydrostatic ( psfc, zm_init, nlevels, thvm, exner, &
                                   z )
      else
        stop "Could not read theta compatable variable"
      endif

    end if

  end subroutine read_z_profile
  !-------------------------------------------------------------------------------------------------
  subroutine read_theta_profile( nvar, nsize, retVars, theta_type, theta )
    !
    !  Description: Searches for the variable specified by either 'thetal[K]' or 'theta[K]' in the
    !  collection of retVars. If the function finds the variable then it returns
    !  it. If it does not the program using this function will exit gracefully
    !  with a warning message.
    !
    !-----------------------------------------------------------------------------------------------
    use input_reader, only: read_x_profile, one_dim_read_var

    use input_names, only: &
      thetal_name, &
      theta_name, &
      temperature_name

    implicit none

    ! Input Variable(s)
    integer, intent(in) :: nvar ! Number of elements in retVars

    integer, intent(in) :: nsize ! Size of value arrays in retVars

    type(one_dim_read_var), dimension(nvar), intent(in) :: retVars ! Collection being
    !                                                                searched through

    character(len=*), intent(out) :: theta_type ! Indicates type read in

    ! Output Variable(s)
    real, dimension(nsize), intent(out) :: theta

    if( count( (/ any(retVars%name == theta_name), &
                  any(retVars%name == thetal_name), &
                  any(retVars%name == temperature_name) /) )<= 1) then
      if( any(retVars%name == theta_name))then
        theta_type = theta_name
      elseif( any(retVars%name == thetal_name))then
        theta_type = thetal_name
      elseif( any(retVars%name == temperature_name))then
        theta_type = temperature_name
      else
        stop "Could not read theta compatable variable"
      endif
      theta = read_x_profile( nvar, nsize, theta_type, retVars )

    end if
  end subroutine read_theta_profile
  !-------------------------------------------------------------------------------------------------
  subroutine read_subs_profile( nvar, nsize, retVars, subs_type, subs )
    !
    !  Description: Searches for the variable specified by either 'w[m\s]' or 'omega[Pa\s]' in the
    !  collection of retVars. If the function finds the variable then it returns
    !  it. If it does not the program using this function will exit gracefully
    !  with a warning message.
    !
    !-----------------------------------------------------------------------------------------------
    use input_reader, only: one_dim_read_var, read_x_profile

    use input_names, only : &
      wm_name, &
      omega_name

    implicit none

    ! Input Variable(s)
    integer, intent(in) :: nvar ! Number of elements in retVars

    integer, intent(in) :: nsize ! Size of the value arrays in retVars

    type(one_dim_read_var), dimension(nvar), intent(in) :: retVars ! Collection being
    !                                                                searched through

    ! Output Variable(s)
    character(len=*), intent(out) :: subs_type ! Indicates type of subsidence measurement

    real, dimension(nsize), intent(out) :: subs ! Subsidence profile [m/s or Pa/s]

    if( count( (/ any(retVars%name == wm_name), any(retVars%name == omega_name) /)) <= 1) then
      if( any(retVars%name == wm_name))then
        subs_type = wm_name
      elseif( any(retVars%name == omega_name))then
        subs_type = omega_name
      else
        stop "Could not read vertical velocity compatable variable"
      endif
      subs = read_x_profile(nvar, nsize, subs_type, retVars)

    end if
  end subroutine read_subs_profile

end module input_interpret
