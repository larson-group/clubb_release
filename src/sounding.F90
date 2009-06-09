! $Id$
module sounding

  implicit none

  public ::  & 
    read_sounding, & 
    read_x_profile, &
    read_profile ! Not currently used in CLUBB

  private :: read_sounding_file, read_sclr_sounding_file, &
    read_edsclr_sounding_file, read_z_profile, &
    read_theta_profile


  ! Constant parameter
  integer, public, parameter :: nmaxsnd = 600
  integer, public, parameter :: sclr_max = 1000

  ! Column identifiers
  character(len=*), public, parameter :: z_name = 'z[m]'
  character(len=*), public, parameter :: pressure_name = 'Press[Pa]'
  character(len=*), public, parameter :: temp_name = 'T[K]'
  character(len=*), public, parameter :: theta_name = 'thm[K]'
  character(len=*), public, parameter :: thetal_name = 'thlm[K]'
  character(len=*), public, parameter :: wm_name = 'w[m\s]'
  character(len=*), public, parameter :: omega_name = 'omega[Pa\s]'

  private ! Default Scope

  contains
  !------------------------------------------------------------------------
  subroutine read_sounding( iunit, runtype, psfc, zm_init,& 
                            thlm, theta_type, rtm, um, vm, ugm, vgm, &
                            alt_type, press, subs_type, wm, &
                            sclrm, edsclrm )

    !       Description:
    !       Subroutine to initialize model variables from a namelist file
    !       References:
    !       None
    !------------------------------------------------------------------------

    use grid_class, only:  & 
        gr ! Variable(s)

    use constants, only:  & 
        fstderr, & ! Constant
        fstdout

    use parameters_model, only: & 
        sclr_dim, &! Variable(s)
        edsclr_dim

    use std_atmosphere_mod, only:  & 
        std_atmosphere ! Procedure(s)

    use interpolation, only:  & 
        lin_int ! Procedure(s)

    use array_index, only: & 
        iisclr_rt, &  ! Variable
        iisclr_thl
    !           ,iisclr_CO2

    use error_code, only: &
      clubb_at_least_debug_level ! Function

    implicit none


    ! Constant parameter
    integer, parameter :: nmaxsnd = 600

    ! Input variables
    integer, intent(in) :: iunit ! File unit to use for namelist

    character(len=*), intent(in) ::  & 
    runtype      ! String for DYCOMS II RF02


    real, intent(in) :: &
      psfc, & ! Pressure at the surface [Pa]
      zm_init ! Height at zm(1)         [m]

    ! Output variables
    real, intent(out), dimension(gr%nnzp) ::  & 
    thlm,  & ! Liquid potential temperature    [K]
    rtm,   & ! Total water mixing ratio        [kg/kg]
    um,    & ! u wind                          [m/s]
    vm,    & ! v wind                          [m/s]
    ugm,   & ! u geostrophic wind              [m/s]
    vgm,   & ! v geostrophic wind              [m/s]
    press, & ! Pressure                        [Pa]
    wm       ! Subsidence                      [m/s or Pa/s]

    character(len=*), intent(out) :: &
      theta_type, &     ! Type of temperature sounding
      alt_type, &       ! Type of independant coordinate
      subs_type         ! Type of subsidence

    ! Optional output variables
    real, intent(out), dimension(gr%nnzp, sclr_dim) ::  & 
      sclrm   ! Passive scalar output      [units vary]

    real, intent(out), dimension(gr%nnzp, edsclr_dim) ::  & 
      edsclrm ! Eddy Passive scalar output [units vary]

    ! Local variables

    ! Input variables from namelist
    integer :: nlevels  ! Levels in the input sounding

    logical :: sounding_exists = .false.
    logical :: sclr_sounding_exists = .false.
    logical :: edsclr_sounding_exists = .false.


    real, dimension(nmaxsnd) :: & 
    z,      & ! Altitude                               [m]
    theta,  & ! Liquid potential temperature sounding  [K]
    rt,     & ! Total water mixing ratio sounding      [kg/kg]
    u,      & ! u wind sounding                        [m/s] 
    v,      & ! v wind sounding                        [m/s]
    ug,     & ! u geostrophic wind sounding            [m/s]
    vg,     & ! v geostrophic wind sounding            [m/s]
    p_in_Pa, &   ! Pressure                               [Pa]
    subs      ! Subsidence                             [m/s or Pa/s]

    real, dimension(nmaxsnd, sclr_max) ::  & 
    sclr, edsclr ! Passive scalar input sounding    [units vary]

    integer :: i, j, k  ! Loop indices

    ! Is this model being extended by 1976 Standard Atmosphere?
    logical :: l_std_atmo

    l_std_atmo = .false.

    theta_type = theta_name ! Default value
    alt_type = z_name ! Default value
    subs_type = wm_name ! Defuault Value

    ! Determine which files exist ahead of time to allow for a graceful exit if
    ! one is missing.
    inquire(file="../input/case_setups/"//trim(runtype)//"_sounding.in", exist=sounding_exists)

    inquire(file="../input/case_setups/"//trim(runtype)//"_sclr_sounding.in", &
      exist=sclr_sounding_exists)

    inquire(file="../input/case_setups/"//trim(runtype)//"_edsclr_sounding.in", &
      exist=edsclr_sounding_exists)

    !---------------------------------------------------------------------------------------------
    ! Status Message
    if( clubb_at_least_debug_level(1) ) then
      print *, "Path to sounding: ", trim(runtype)//'_sounding.in'
      print *, "File exists? ", sounding_exists
      print *, "Path to sclr_sounding: ", trim(runtype)//'_sclr_sounding.in'
      print *, "File exists? ", sclr_sounding_exists
      print *, "Path to sounding: ", trim(runtype)//'_edsclr_sounding.in'
      print *, "File exists? ", edsclr_sounding_exists
    end if
    !----------------------------------------------------------------------------------------------

    if( sounding_exists ) then
      ! Read in SAM-Like <runtype>_sounding.in file
      call read_sounding_file( iunit, runtype, nlevels, psfc, zm_init, & 
                               z, theta, theta_type, rt, u, v, ug, vg, &
                               alt_type, p_in_Pa, subs_type, subs )
    else
      stop 'Cannot open <runtype>_sounding.in file'
      ! sounding namelist is no longer used.
      ! Joshua Fasching April 2009
    end if
    ! Read in a passive scalar sounding, if enabled
    if ( sclr_dim > 0 .or. edsclr_dim > 0 ) then
      ! Initialize to zero
      sclr   = 0.0
      edsclr = 0.0
      ! Read in SAM-Like <runtype>_sclr_sounding.in and
      !                  <runtype>_edsclr_sounding.in
      if( sclr_dim > 0 ) then
        if( sclr_sounding_exists ) then
          call read_sclr_sounding_file( iunit, runtype, sclr )
        else
          stop 'Cannot open <runtype>_sclr_sounding.in file'
        end if
      end if
      if( edsclr_dim > 0 ) then
        if( edsclr_sounding_exists  ) then
          call read_edsclr_sounding_file( iunit, runtype, edsclr )
        else
          stop 'Cannot open <runtype>_edsclr_sounding.in file'
        end if
      end if
    end if

    if ( nlevels > nmaxsnd ) then
      write(fstderr,*) 'Error in sounding: nlevels > nmaxsnd'
      write(fstderr,*) 'nlevels = ',nlevels
      write(fstderr,*) 'nmaxsnd = ',nmaxsnd
      stop 'STOP in read_sounding'
    end if

    ! Error check: if lowest above-model-surface themodynamic grid height
    ! (gr%zt(2)) is lower than the lowest value from the input sounding,
    ! then the linear interpolation scheme will fail.

    if ( gr%zt(2) < z(1) ) then
      write(fstderr,*) "Lowest level of input sounding, z(1), must be",  &
      " below the first above-model-surface thermodynamic level, gr%zt(2)"
      write(fstderr,*) " First sounding level z(1) = ", z(1)
      write(fstderr,*) " First thermodynamic level gr%zt(2) = ", gr%zt(2)
      stop 'STOP in read_sounding'
    endif

    ! First sounding level should be near ground value

    ! dschanen 1 May 2007
    ! We have changed this for Nov. 11 and June 25, both of which
    ! begin above the ground.
    ! if ( abs(z(1)) > 1.e-8 ) then
    if ( .false. ) then
      write(fstderr,*) 'First level of input sounding must be z=0'
      stop 'STOP in read_sounding'
    else
      um(1)   = u(1)
      vm(1)   = v(1)
      ugm(1)  = ug(1)
      vgm(1)  = vg(1)
      thlm(1) = theta(1)
      rtm(1)  = rt(1)
      press(1) = p_in_Pa(1)
      wm(1) = subs(1)
      if ( sclr_dim > 0 ) then
        sclrm(1,1:sclr_dim)   = sclr(1,1:sclr_dim)
      end if
      if ( edsclr_dim > 0 ) then
        edsclrm(1,1:edsclr_dim) = edsclr(1,1:edsclr_dim)
      end if
    end if

    if ( clubb_at_least_debug_level( 1 ) ) then

      write(fstdout,*) "Reading in sounding information"
      !------------Printing Model Inputs-------------------------------
      write(fstdout,*) "z = ", z(1:nlevels)
      write(fstdout,*) "u = ", u(1:nlevels)
      write(fstdout,*) "v = ", v(1:nlevels)
      write(fstdout,*) "ug = ", ug(1:nlevels)
      write(fstdout,*) "vg = ", vg(1:nlevels)
      write(fstdout,*) "theta = ", theta(1:nlevels)
      write(fstdout,*) "rt = ", rt(1:nlevels)
      write(fstdout,*) "p_in_Pa = ", p_in_Pa(1:nlevels)
      write(fstdout,*) "subs = ", subs(1:nlevels)

      do i = 1, sclr_dim, 1
        write(fstdout,'(a5,i2,a2)',advance='no') "sclr(", i,") = "
        write(fstdout,'(8g10.3)') sclr(1:nlevels,i)
      enddo

      do i = 1, edsclr_dim, 1
        write(fstdout,'(a7,i2,a2)',advance='no') "edsclr(", i, ") = "
        write(fstdout,'(8g10.3)') edsclr(1:nlevels,i)
      enddo

    endif ! clubb_at_least_debug_level( 1 )
    !----------------------------------------------------------------------
    ! Use linear interpolation from two nearest prescribed grid points
    ! (one above and one below) to initialize mean quantities in the model
    ! Modified 27 May 2005 -dschanen: eliminated the goto in favor of a do while( )
    do i=2, gr%nnzp
      k=1
      do while ( z(k) < gr%zt(i) .and. .not. l_std_atmo )
        k=k+1
        if ( k > nlevels ) then
!              write(fstderr,*) 'STOP Not enough sounding data to ',
!     .          'initialize grid:'
!              write(fstderr,'(a,f7.1,/a,f7.1)')
!     .          '  highest sounding level', z(nlevels),
!     .          '  should be higher than highest thermodynamic point',
!     .          gr%zt(gr%nnzp)
!              stop 'STOP in sounding'

          l_std_atmo = .true.
          exit
        end if  ! k > nlevels

        ! Regular situation w/ linear int.
        IF ( trim( runtype ) /= "dycoms2_rf02" ) THEN

          um(i)   = lin_int( gr%zt(i), z(k), z(k-1), u(k), u(k-1) )
          vm(i)   = lin_int( gr%zt(i), z(k), z(k-1), v(k), v(k-1) )
          ugm(i)  = lin_int( gr%zt(i), z(k), z(k-1), ug(k), ug(k-1) )
          vgm(i)  = lin_int( gr%zt(i), z(k), z(k-1), vg(k), vg(k-1) )
          thlm(i) = lin_int( gr%zt(i), z(k), z(k-1), theta(k), theta(k-1) )
          rtm(i)  = lin_int( gr%zt(i), z(k), z(k-1), rt(k), rt(k-1) )
          press(i) = lin_int( gr%zt(i), z(k), z(k-1), p_in_Pa(k), p_in_Pa(k-1) )
          wm(i) = lin_int( gr%zt(i), z(k), z(k-1), subs(k), subs(k-1) )

          if ( sclr_dim > 0 ) then
            do j = 1, sclr_dim
              sclrm(i,j) = lin_int( gr%zt(i), z(k), z(k-1),  & 
                                    sclr(k,j), sclr(k-1,j) )
            end do
          end if
          if ( edsclr_dim > 0 ) then
            do j = 1, edsclr_dim
              edsclrm(i,j) = lin_int( gr%zt(i), z(k), z(k-1),  & 
                                      edsclr(k,j), edsclr(k-1,j) )
            end do
          end if

        ELSE  ! DYCOMS II RF02 case

          IF ( gr%zt(i) < 795.0 ) THEN
            um(i)   =  3.0 + (4.3*gr%zt(i))/1000.0
            vm(i)   = -9.0 + (5.6*gr%zt(i))/1000.0
            ugm(i)  = um(i)
            vgm(i)  = vm(i)
            thlm(i) = 288.3
            rtm(i)  = (9.45)/1000.0
            ! Passive Scalars
            ! Change this if they are not equal to theta_l and rt in RF02
            if ( iisclr_thl > 0  ) then
              sclrm(i, iisclr_thl)   = thlm(i)
              edsclrm(i, iisclr_thl) = thlm(i)
            end if
            if ( iisclr_rt > 0 ) then
              sclrm(i, iisclr_rt)    = rtm(i)
              edsclrm(i, iisclr_rt)  = rtm(i)
            end if
            press(i) = lin_int( gr%zt(i), z(k), z(k-1), p_in_Pa(k), p_in_Pa(k-1) )
            wm(i) = lin_int( gr%zt(i), z(k), z(k-1), subs(k), subs(k-1) )
          ELSE
            um(i)   =  3.0 + (4.3*gr%zt(i))/1000.0
            vm(i)   = -9.0 + (5.6*gr%zt(i))/1000.0
            ugm(i)  = um(i)
            vgm(i)  = vm(i)
            thlm(i) = 295.0 + ( (gr%zt(i) - 795.0)**(1.0/3.0) )
            rtm(i)  = (  5.0 - 3.0  & 
            * ( 1.0 - EXP( (795.0 - gr%zt(i))/500.0 ) )  )/1000.0
            ! Passive Scalars
            ! Same as above
            if ( iisclr_thl > 0  ) then
              sclrm(i, iisclr_thl)   = thlm(i)
              edsclrm(i, iisclr_thl) = thlm(i)
            end if
            if ( iisclr_rt > 0 ) then
              sclrm(i, iisclr_rt)    = rtm(i)
              edsclrm(i, iisclr_rt)  = rtm(i)
            end if
            press(i) = lin_int( gr%zt(i), z(k), z(k-1), p_in_Pa(k), p_in_Pa(k-1) )
            wm(i) = lin_int( gr%zt(i), z(k), z(k-1), subs(k), subs(k-1) )
          END IF

        END IF ! runtype

      end do ! do while ( z(k) < gr%zt(i) )

      ! If the grid extends beyond the sounding data, use
      ! Standard Atmosphere
      ! Joshua Fasching April 2009
      if ( l_std_atmo ) then
        call std_atmosphere( gr%zt(i), thlm(i), rtm(i), press(i) )

        um(i) = um(i-1)
        vm(i) = vm(i-1)
        wm(i) = wm(i-1)
        ugm(i) = um(i)
        vgm(i) = vm(i)

      end if

    end do   ! i=2, gr%nnzp

    if ( l_std_atmo .and.  &
         clubb_at_least_debug_level( 1 ) ) then
      write(fstderr,*) "Warning:  1976 Standard Atmosphere "// & 
                       "was used to complete the grid."
    endif

    return
  end subroutine read_sounding

  !-------------------------------------------------------------------------------------------------
  subroutine read_sounding_file( iunit, runtype, nlevels, psfc, zm_init, &
                                 z, theta, theta_type, rt, u, v, ug, vg, &
                                 alt_type, p_in_Pa, subs_type, subs )
    !
    !  Description: This subroutine reads in a <runtype>_sounding.in file and
    !  returns the values contained in that file.
    !
    !-----------------------------------------------------------------------------------------------
    use input_reader, only: read_one_dim_file, fill_blanks_one_dim_vars, &
                            one_dim_read_var, deallocate_one_dim_vars

    implicit none

    ! Input Variable(s)
    integer, intent(in) :: iunit ! I/O unit

    character(len=*), intent(in) :: runtype ! String identifying the model case;
    !                                         e.g. bomex

    real, intent(in) :: &
      psfc, & ! Pressure at the surface [Pa]
      zm_init ! Height at zm(1)         [m]

    ! Output Variable(s)
    integer, intent(out) :: nlevels ! Number of levels from the sounding.in file

    real, intent(out), dimension(nmaxsnd) :: & 
    z,      & ! Altitude                               [m]
    theta,  & ! Liquid potential temperature sounding  [K]
    rt,     & ! Total water mixing ratio sounding      [kg/kg]
    u,      & ! u wind sounding                        [m/s] 
    v,      & ! v wind sounding                        [m/s]
    ug,     & ! u geostrophic wind sounding            [m/s]
    vg,     & ! v geostrophic wind sounding            [m/s]
    p_in_Pa,& ! Pressure sounding                      [Pa]
    subs      ! Subsidence sounding                    [m/s or Pa/s]

    character(len=*), intent(out) :: & 
      theta_type, &     ! Type of temperature sounding
      alt_type, &       ! Type of independant coordinate
      subs_type         ! Type of subsidence

    integer, parameter :: nCol = 8

    type(one_dim_read_var), dimension(nCol) :: retVars

    call read_one_dim_file( iunit, nCol, &
    '../input/case_setups/'//trim(runtype)//'_sounding.in', retVars )

    call fill_blanks_one_dim_vars( nCol, retVars )

    call read_z_profile(nCol, retVars, psfc, zm_init, z, p_in_Pa, alt_type )

    call read_theta_profile(nCol, retVars, theta_type, theta)

    rt = read_x_profile(nCol, nmaxsnd, 'rt[kg\kg]', retVars)

    u = read_x_profile(nCol, nmaxsnd, 'u[m\s]', retVars)

    v = read_x_profile(nCol, nmaxsnd, 'v[m\s]', retVars)

    ug = read_x_profile(nCol, nmaxsnd, 'ug[m\s]', retVars)

    vg = read_x_profile(nCol, nmaxsnd, 'vg[m\s]', retVars)

    call read_subs_profile(nCol, retVars, subs_type, subs)

    nlevels = size(retVars(1)%values)

    call deallocate_one_dim_vars( nCol, retVars )

  end subroutine read_sounding_file

  !-------------------------------------------------------------------------------------------------
  subroutine read_sclr_sounding_file( iunit, runtype, sclr )
    !
    !  Description: This subroutine reads in a <runtype>_sclr_sounding.in file and
    !  returns the values contained in that file.
    !
    !-----------------------------------------------------------------------------------------------
    use input_reader, only: read_one_dim_file, fill_blanks_one_dim_vars, &
                            one_dim_read_var, deallocate_one_dim_vars

    use parameters_model, only: sclr_dim

    use array_index, only: iisclr_rt, iisclr_thl, iisclr_CO2

    implicit none

    ! Input Variable(s)
    integer, intent(in) :: iunit ! I/O unit

    character(len=*), intent(in) :: runtype ! String identifying the model case;
    !                                         e.g. bomex

    ! Output Variable(s)
    real, intent(out), dimension(nmaxsnd,sclr_max) :: & 
    sclr        ! Scalar sounding [?]


    type(one_dim_read_var), dimension(sclr_dim) :: retVars

    integer i

    call read_one_dim_file( iunit, sclr_dim, &
    '../input/case_setups/'//trim(runtype)//'_sclr_sounding.in', retVars )

!    call fill_blanks_one_dim_vars( sclr_dim, retVars )

    do i=1, sclr_dim
      select case(trim(retVars(i)%name))
      case("CO2[ppmv]")
        if( i /= iisclr_CO2 .and. iisclr_CO2 > 0) then
          stop "iisclr_CO2 index does not match column."
        end if
      case('rt[kg\kg]')
        if( i /= iisclr_rt .and. iisclr_rt > 0) then
          stop "iisclr_rt index does not match column."
        end if
      case("thl[K]")
        if( i /= iisclr_thl .and. iisclr_thl > 0) then
          stop "iisclr_thl index does not match column."
        end if
      end select
      sclr(1:size(retVars(i)%values),i) = retVars(i)%values
    end do

    call deallocate_one_dim_vars( sclr_dim, retVars)

  end subroutine read_sclr_sounding_file

  !-------------------------------------------------------------------------------------------------
  subroutine read_edsclr_sounding_file( iunit, runtype, edsclr )
    !
    !  Description: This subroutine reads in a <runtype>_edsclr_sounding.in file and
    !  returns the values contained in that file.
    !
    !-----------------------------------------------------------------------------------------------
    use input_reader, only: read_one_dim_file, fill_blanks_one_dim_vars, &
                            one_dim_read_var, deallocate_one_dim_vars

    use parameters_model, only: edsclr_dim

    use array_index, only: iiedsclr_rt, iiedsclr_thl, iiedsclr_CO2

    implicit none

    ! Input Variable(s)
    integer, intent(in) :: iunit ! I/O unit

    character(len=*), intent(in) :: runtype ! String identifying the model case;
    !                                         e.g. bomex

    ! Output Variable(s)
    real, intent(out), dimension(nmaxsnd,sclr_max) :: & 
    edsclr ! Eddy Scalars [?]

    type(one_dim_read_var), dimension(edsclr_dim) :: retVars

    integer i

    call read_one_dim_file( iunit, edsclr_dim, &
    '../input/case_setups/'//trim(runtype)//'_edsclr_sounding.in', retVars )

 !   call fill_blanks_one_dim_vars( edsclr_dim, retVars )

    do i=1, edsclr_dim
      select case(trim(retVars(i)%name))
      case("CO2[ppmv]")
        if( i /= iiedsclr_CO2 .and. iiedsclr_CO2 > 0) then
          stop "iisclr_CO2 index does not match column."
        end if
      case('rt[kg\kg]')
        if( i /= iiedsclr_rt .and. iiedsclr_rt > 0) then
          stop "iisclr_rt index does not match column."
        end if
      case("thl[K]")
        if( i /= iiedsclr_thl .and. iiedsclr_thl > 0) then
          stop "iisclr_thl index does not match column."
        end if
      end select
      edsclr(1:size(retVars(i)%values),i) = retVars(i)%values
    end do

    call deallocate_one_dim_vars( edsclr_dim, retVars )

  end subroutine read_edsclr_sounding_file


  !-------------------------------------------------------------------------------------------------
  function read_x_profile( nvar, dim_size, target_name, retVars ) result(x)
    !
    !  Description: Searches for the variable specified by target_name in the
    !  collection of retVars. If the function finds the variable then it returns
    !  it. If it does not the program using this function will exit gracefully
    !  with a warning message.
    !
    !-----------------------------------------------------------------------------------------------
    use input_reader, only: one_dim_read_var

    implicit none


    ! Input Variable(s)
    integer, intent(in) :: nvar ! Number of variables in retVars

    integer, intent(in) :: dim_size


    character(len=*), intent(in) :: target_name ! Variable that is being
    !                                             searched for

    type(one_dim_read_var), dimension(nvar), intent(in) :: retVars ! Collection
    !                                                                being searched through

    ! Output Variable(s)
    real, dimension(dim_size) :: x

    ! Local Variables
    integer i

    logical l_found

    l_found = .false.

    i = 1
    do while( i <= nvar .and. .not. l_found)
      if( retVars(i)%name == target_name ) then
        l_found = .true.
        x(1:size(retVars(i)%values)) = retVars(i)%values
      end if
      i=i+1
    end do

    if( .not. l_found ) then
      stop ' Profile could not be found. Check your sounding.in file.'
    end if

  end function read_x_profile

!------------------------------------------------------------------------------
  subroutine read_z_profile(nvar, retVars, psfc, zm_init, z, p_in_Pa, alt_type)
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

    use input_reader, only: one_dim_read_var ! Procedure(s)

    use constants, only: kappa, p0, Cp, Lv, zero_threshold, ep2, ep1 ! Variable(s)

    use saturation, only: sat_mixrat_liq, sat_rcm ! Procedure(s)

    use parameters_model, only: T0 ! Variable(s)

    use hydrostatic_mod, only: inverse_hydrostatic ! Procedure(s)

    implicit none


    ! Input Variable(s)
    integer, intent(in) :: nvar ! Number of elements in retVars

    type(one_dim_read_var), dimension(nvar), intent(in) :: retVars ! Collection
    !                                                                being searched

    real, intent(in) :: &
      psfc, &         ! Pressure at the surface [Pa]
      zm_init         ! Height at zm(1)         [m]

    ! Output Variable(s)

    real, intent(out), dimension(nmaxsnd) :: z ! Height sounding profile [m]

    real, intent(out), dimension(nmaxsnd) :: p_in_Pa ! Pressure sounding profile [Pa]

    character(len=*), intent(out) :: alt_type ! Indicates where altitudes were
    !                                           gained from

    intrinsic :: max

    ! Local Variables
    real, dimension(nmaxsnd) :: exner,thvm, rcm, theta, rtm

    integer :: nlevels, k

    character(len=40) :: theta_type

    if( count( (/ any(retVars%name == z_name), any(retVars%name == pressure_name) /)) <= 1) then
      if( any(retVars%name == z_name))then
        alt_type = z_name
        z = read_x_profile( nvar, nmaxsnd, alt_type, retVars )
        p_in_Pa = -999.9

      elseif( any(retVars%name == pressure_name))then
        alt_type = pressure_name

        p_in_Pa = read_x_profile( nvar, nmaxsnd, alt_type, retVars )

        nlevels = size(retVars(1)%values)

        call read_theta_profile(nvar, retVars, theta_type, theta )

        rtm = read_x_profile(nvar, nmaxsnd, 'rt[kg\kg]', retVars)

        exner(1) = ( psfc/p0 )**kappa

        do k=2, nlevels
          exner(k) = (p_in_Pa(k)/p0) ** kappa  ! zt
        end do

        if( trim( theta_type ) == temp_name ) then
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
          theta = theta - Lv/(Cp*exner) * rcm

        end select

        ! Now, compute initial thetav

        thvm = theta + ep1 * T0 * rtm  & 
                    + ( Lv/(Cp*exner) - ep2 * T0 ) * rcm

        call inverse_hydrostatic ( thvm, zm_init, exner, nlevels, &
                                     z )
      else
        stop "Could not read theta compatable variable"
      endif

    end if

  end subroutine read_z_profile

  !-------------------------------------------------------------------------------------------------
  subroutine read_theta_profile(nvar, retVars, theta_type, theta)
    !
    !  Description: Searches for the variable specified by either 'thetal[K]' or 'theta[K]' in the
    !  collection of retVars. If the function finds the variable then it returns
    !  it. If it does not the program using this function will exit gracefully
    !  with a warning message.
    !
    !-----------------------------------------------------------------------------------------------
    use input_reader, only: one_dim_read_var

    implicit none

    ! Input Variable(s)
    integer, intent(in) :: nvar ! Number of elements in retVars

    type(one_dim_read_var), dimension(nvar), intent(in) :: retVars ! Collection being
    !                                                                searched through

    character(len=*), intent(out) :: theta_type

    ! Output Variable(s)
    real, dimension(nmaxsnd), intent(out) :: theta

    if( count( (/ any(retVars%name == theta_name), &
                  any(retVars%name == thetal_name), &
                  any(retVars%name == temp_name) /) )<= 1) then
      if( any(retVars%name == theta_name))then
        theta_type = theta_name
      elseif( any(retVars%name == thetal_name))then
        theta_type = thetal_name
      elseif( any(retVars%name == temp_name))then
        theta_type = temp_name
      else
        stop "Could not read theta compatable variable"
      endif
      theta = read_x_profile(nvar, nmaxsnd,theta_type, retVars)

    end if
  end subroutine read_theta_profile
  !-------------------------------------------------------------------------------------------------
  subroutine read_subs_profile(nvar, retVars, subs_type, subs )
    !
    !  Description: Searches for the variable specified by either 'w[m\s]' or 'omega[Pa\s]' in the
    !  collection of retVars. If the function finds the variable then it returns
    !  it. If it does not the program using this function will exit gracefully
    !  with a warning message.
    !
    !-----------------------------------------------------------------------------------------------
    use input_reader, only: one_dim_read_var

    implicit none


    ! Input Variable(s)
    integer, intent(in) :: nvar ! Number of elements in retVars

    type(one_dim_read_var), dimension(nvar), intent(in) :: retVars ! Collection being
    !                                                                searched through

    ! Output Variable(s)
    character(len=*), intent(out) :: subs_type ! Indicates type of subsidence measurement

    real, dimension(nmaxsnd), intent(out) :: subs ! Subsidence profile [m/s or Pa/s]

    if( count( (/ any(retVars%name == wm_name), any(retVars%name == omega_name) /)) <= 1) then
      if( any(retVars%name == wm_name))then
        subs_type = wm_name
      elseif( any(retVars%name == omega_name))then
        subs_type = omega_name
      else
        stop "Could not read vertical velocity compatable variable"
      endif
      subs = read_x_profile(nvar, nmaxsnd, subs_type, retVars)

    end if
  end subroutine read_subs_profile

  !------------------------------------------------------------------------
  subroutine read_profile( fname, x )

    !       Description:
    !       Subroutine to initialize one generic model variable from file
    !------------------------------------------------------------------------

    use grid_class, only:  & 
        gr ! Variable(s)

    use interpolation, only:  & 
        lin_int ! Procedure

    use constants, only:  &
        fstderr ! Constant

    implicit none

    ! Constant Parameter

    integer, parameter :: nmaxsnd = 200

    ! Input Variables
    character(len=*), intent(in) :: fname

    ! Output

    real, dimension(gr%nnzp), intent(out) :: x

    ! Local variables

    ! Input variables from namelist

    integer :: nlevels

    real, dimension(nmaxsnd) :: z, var

    namelist /profile/ nlevels, z, var

    ! Loop indices
    integer :: i,k

    ! Read sounding namelist

    open(10, file = trim(fname), status = 'old' )
    read(10, nml = profile )
    close(10)

    if ( nlevels > nmaxsnd ) then
      write(fstderr,*) 'Error in sounding: nlevels > nmaxsnd'
      write(fstderr,*) 'nlevels = ', nlevels
      write(fstderr,*) 'nmaxsnd = ', nmaxsnd
      stop 'STOP in read_profile (sounding.F90)'
    endif

    ! Error check: if lowest above-model-surface themodynamic grid height
    ! (gr%zt(2)) is lower than the lowest value from the input sounding,
    ! then the linear interpolation scheme will fail.

    if ( gr%zt(2) < z(1) ) then
      write(fstderr,*) "Lowest level of input sounding, z(1), must be",  &
      " below the first above-model-surface thermodynamic level, gr%zt(2)"
      write(fstderr,*) " First sounding level z(1) = ", z(1)
      write(fstderr,*) " First thermodynamic level gr%zt(2) = ", gr%zt(2)
      stop 'STOP in read_profile (sounding.F90)'
    endif

    ! Use linear interpolation from two nearest prescribed grid points
    ! (one above and one below) to initialize mean quantities in the model
    ! Modified 27 May 2005 -dschanen: eliminated the goto in favor of a do while( )

    do i = 2, gr%nnzp
      k = 1
      do while ( z(k) < gr%zt(i) )
        k = k + 1
        if ( k > nlevels ) then
          write(fstderr,*) 'STOP Not enough sounding data to ',  & 
                           'initialize grid:'
          write(fstderr,'(a,f7.1,/a,f7.1)') ' Highest sounding level',  &
               z(nlevels),  &
               'should be higher than highest thermodynamic point',  &
               gr%zt(gr%nnzp)
          write(*,*) ' Filename: ', fname
          stop 'STOP in read_profile (sounding.F90)'
        endif
        x(i) = lin_int( gr%zt(i), z(k), z(k-1), var(k), var(k-1) )
      enddo ! while
    enddo ! i=2, gr%nzzp

    return
  end subroutine read_profile

end module sounding
!-----------------------------------------------------------------------
