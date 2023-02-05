!$Id$
module input_interpret

  implicit none  

  public :: &
    read_z_profile, &
    read_theta_profile, &
    read_subs_profile
  private ! Default scope

  contains

  !-----------------------------------------------------------------------------
  subroutine read_z_profile( nvar, nsize, retVars, p_sfc, zm_init, &
                             z, p_in_Pa, alt_type )

    ! Description:
    ! Searches for the variable specified by either 'z[m]' or 'Press[Pa]' in the
    ! collection of retVars.  If the subroutine finds the variable indicated by
    ! 'z[m]', then it returns it.  If the subroutine finds 'Press[Pa]' then it
    ! converts it to values of altitude in meters.  If it does not find either
    ! or finds both the program using this subroutine will exit gracefully with
    ! a warning message.

    !-----------------------------------------------------------------------

    use input_reader, only: &
        read_x_profile, & ! Prodedure(s)
        one_dim_read_var  ! Type

    use constants_clubb, only: &
        kappa, & ! Constant(s)
        p0, &
        Cp, &
        Lv, &
        ep2, &
        one, &
        zero_threshold, &
        fstderr

    use saturation, only: &
        sat_mixrat_liq, & ! Procedure(s)
        rcm_sat_adj

    use hydrostatic_module, only: &
        inverse_hydrostatic ! Procedure(s)

    use calc_pressure, only: &
        calculate_thvm    ! Procedure(s)

    use input_names, only: &
        z_name, & ! Variable(s)
        pressure_name, &
        rt_name, &
        temperature_name, &
        thetal_name, &
        theta_name

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variable(s)
    integer, intent(in) :: nvar ! Number of elements in retVars

    integer, intent(in) :: nsize

    type(one_dim_read_var), dimension(nvar), intent(in) :: &
      retVars    ! Collection being searched

    real( kind = core_rknd ), intent(in) :: &
      p_sfc,   & ! Pressure at the surface    [Pa]
      zm_init    ! Height at zm(1)            [m]

    ! Output Variable(s)
    real( kind = core_rknd ), intent(out), dimension(nsize) :: &
      z ! Height sounding profile [m]

    real( kind = core_rknd ), intent(out), dimension(nsize) :: &
      p_in_Pa ! Pressure sounding profile [Pa]

    character(len=*), intent(out) :: &
      alt_type ! Indicates where altitudes were gained from

    intrinsic :: max

    ! Local Variables
    real( kind = core_rknd ), dimension( 1,nsize ) :: &
      exner,   & ! Exner function at sounding levels                    [-]
      thvm,    & ! Mean theta-v at sounding levels                      [K]
      rcm,     & ! Mean cloud water mixing ratio at sounding levels     [kg/kg]
      theta,   & ! Mean potential temperature at sounding levels        [K]
      rtm,     & ! Mean total water mixing ratio at sounding levels     [kg/kg]
      thlm       ! Mean liquid water potential temp. at sounding levs.  [K]

    integer :: nlevels, k

    character(len=40) :: theta_type


    if( count( (/ any(retVars%name == z_name),  &
                  any(retVars%name == pressure_name) /) ) <= 1) then

      if( any(retVars%name == z_name) ) then

        ! The input sounding is given in terms of altitude.
        alt_type = z_name

        ! Obtain the value of altitude at each sounding level.
        z = read_x_profile( nvar, nsize, alt_type, retVars )

        ! Set the pressure at the sounding levels to the "fill value".
        p_in_Pa = -999.9_core_rknd


      elseif( any(retVars%name == pressure_name) ) then

        ! The input sounding is given in terms of pressure.
        alt_type = pressure_name

        ! Obtain the value of total pressure at each pressure sounding level.
        p_in_Pa = read_x_profile( nvar, nsize, alt_type, retVars )

        nlevels = size(retVars(1)%values)

        ! Obtain the value of theta_type (theta, theta_l, or temperature) at
        ! each pressure sounding level.
        call read_theta_profile(nvar, nsize, retVars, theta_type, theta(1,:) )

        ! Obtain the value of total water mixing ratio at each pressure sounding
        ! level. 
        rtm(1,:) = read_x_profile(nvar, nsize, rt_name, retVars)

        ! Calculate exner from pressure.
        do k = 1, nlevels
          exner(1,k) = ( p_in_Pa(k) / p0 )**kappa
        enddo


        select case ( trim( theta_type ) )

        case ( temperature_name )

           ! The variable "theta" actually contains temperature (in Kelvin) at
           ! this point.

           ! Determine initial cloud water mixing ratio at sounding levels
           ! based on temperature and rtm.  Again, "theta(k)" is actually
           ! "temperature(k)" at this point in the code.
           do k = 1, nlevels
              rcm(1,k) = &
                max( rtm(1,k) - sat_mixrat_liq( p_in_Pa(k), theta(1,k) ), &
                     zero_threshold )
           enddo

           ! Convert temperature to potential temperature, theta.
           theta(1,1:nlevels) = theta(1,1:nlevels) / exner(1,1:nlevels)

           ! Calculate theta_l from theta and cloud water mixing ratio, such
           ! that:  theta_l = theta - [Lv/(Cp*exner)]*rcm.
           thlm(1,1:nlevels) = theta(1,1:nlevels) &
                              - Lv/(Cp*exner(1,1:nlevels)) * rcm(1,1:nlevels)


        case ( theta_name )

           ! This variable "theta" does indeed contain potential temperature.

           ! Determine initial cloud water mixing ratio at sounding levels
           ! based on potential temperature, exner, and rtm.
           do k = 1,nlevels
              rcm(1,k) = &
                max( rtm(1,k) - sat_mixrat_liq( p_in_Pa(k), theta(1,k)*exner(1,k) ), &
                     zero_threshold )
           enddo

           ! Calculate theta_l from theta and cloud water mixing ratio, such
           ! that:  theta_l = theta - [Lv/(Cp*exner)]*rcm.
           thlm(1,1:nlevels) = theta(1,1:nlevels) &
                              - Lv/(Cp*exner(1,1:nlevels)) * rcm(1,1:nlevels)


        case ( thetal_name )

           ! The variable "theta" actually contains liquid water potential
           ! temperature, theta_l, at this point.
           thlm(1,1:nlevels) = theta(1,1:nlevels)

           ! Determine initial cloud water mixing ratio.  If the profile is
           ! unsaturated, then theta = theta_l.  If the profile is saturated at
           ! any level, then cloud water mixing ratio must be determined using
           ! an iterative method involving theta_l, total water mixing ratio,
           ! pressure, and exner.
           do k =1, nlevels, 1
              rcm(1,k) = rcm_sat_adj( thlm(1,k), rtm(1,k), p_in_Pa(k), exner(1,k) )
           enddo

           ! Calculate theta from theta_l and cloud water mixing ratio, such
           ! that:  theta = theta_l + [Lv/(Cp*exner)]*rcm.
           theta(1,1:nlevels) = thlm(1,1:nlevels) &
                              + Lv/(Cp*exner(1,1:nlevels)) * rcm(1,1:nlevels)


        case default

           write(fstderr,*) "Invalid theta_type: ", theta_type
           error stop


        end select


        ! Now, the variable "theta" contains potential temperature, theta.
        ! Compute initial theta_v based on theta, total water mixing ratio, and
        ! cloud water mixing ratio.
        !
        ! theta_v = theta * [ ( 1 + (R_v/R_d)*r_v ) / ( 1 + r_v ) ];
        !
        ! which can be rearranged as:
        !
        ! theta_v = theta * [ 1 + { (R_v/R_d) - 1 } * { r_v / ( 1 + r_v ) } ].
        !
        ! The exact form of the equation for theta_v (including water vapor and
        ! cloud water) is:
        !
        ! theta_v = theta * [ ( 1 + (R_v/R_d)*r_v ) / ( 1 + r_v + r_c ) ];
        !
        ! which can be rearranged as:
        !
        ! theta_v = theta * [ 1 + { (R_v/R_d) - 1 } 
        !                         * { r_v / ( 1 + r_v + r_c ) }
        !                     - { r_c / ( 1 + r_v + r_c ) } ].
        !
        ! This version is written with r_t = r_v + r_c, such that:
        !
        ! theta_v = theta * [ 1 + { (R_v/R_d) - 1 }
        !                         * { ( r_t - r_c ) / ( 1 + r_t ) }
        !                     - { r_c / ( 1 + r_t ) } ].
        !
        ! To use theta_l instead of theta, simply substitute the following for
        ! theta in the above expression:  theta_l + {L_v/(C_p*exner)} * r_c.
        !
        ! The CLUBB code uses a linearized version of the above equation (in
        ! order to calculate thv' terms -- such as w'thv', etc.) for theta_v
        ! throughout the model code, such that:
        !
        ! theta_v = theta_l + { (R_v/R_d) - 1 } * thv_ds * r_t
        !                   + [ {L_v/(C_p*exner)} - (R_v/R_d) * thv_ds ] * r_c;
        !
        ! where thv_ds is used as a reference value to approximate theta_l.
        call calculate_thvm( nlevels, 1, &
             thlm(1,1:nlevels), rtm(1,1:nlevels), rcm(1,1:nlevels), exner(1,1:nlevels), &
             theta(1,1:nlevels) * ( one + ep2 * ( rtm(1,1:nlevels)-rcm(1,1:nlevels) ) )**kappa, &
             thvm(1,1:nlevels) )

        ! Find the altitudes, z, of the pressure sounding levels.
        call inverse_hydrostatic( p_sfc, zm_init, nlevels, thvm, exner, &
                                  z )


      else

         error stop "Could not read sounding vertical-coordinate variable"


      endif ! z_ name or pressure_name

    endif

  end subroutine read_z_profile
  !-----------------------------------------------------------------------------
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

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variable(s)
    integer, intent(in) :: nvar ! Number of elements in retVars

    integer, intent(in) :: nsize ! Size of value arrays in retVars

    type(one_dim_read_var), dimension(nvar), intent(in) :: retVars ! Collection being
    !                                                                searched through

    character(len=*), intent(out) :: theta_type ! Indicates type read in

    ! Output Variable(s)
    real( kind = core_rknd ), dimension(nsize), intent(out) :: theta

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
        error stop "Could not read theta compatable variable"
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
 
    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variable(s)
    integer, intent(in) :: nvar ! Number of elements in retVars

    integer, intent(in) :: nsize ! Size of the value arrays in retVars

    type(one_dim_read_var), dimension(nvar), intent(in) :: &
      retVars ! Collection being searched through

    ! Output Variable(s)
    character(len=*), intent(out) :: &
      subs_type ! Indicates type of subsidence measurement

    real( kind = core_rknd ), dimension(nsize), intent(out) :: &
      subs ! Subsidence profile [m/s or Pa/s]

    if( count( (/ any(retVars%name == wm_name), any(retVars%name == omega_name) /)) <= 1) then
      if( any(retVars%name == wm_name))then
        subs_type = wm_name
      elseif( any(retVars%name == omega_name))then
        subs_type = omega_name
      else
        error stop "Could not read vertical velocity compatable variable"
      endif
      subs = read_x_profile(nvar, nsize, subs_type, retVars)

    end if
  end subroutine read_subs_profile

end module input_interpret
