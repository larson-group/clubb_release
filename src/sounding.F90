! $Id$
module sounding

  implicit none

  public ::  & 
    read_sounding, & 
    read_profile ! Not currently used in CLUBB

  private :: read_sounding_file, read_sclr_sounding_file, &
    read_edsclr_sounding_file


  ! Constant parameter
  integer, public, parameter :: &
    nmaxsnd   = 600, &
    sclr_max  = 1000, &
    n_snd_var = 8

  private ! Default Scope

  contains
  !------------------------------------------------------------------------
  subroutine read_sounding( iunit, runtype, psfc, zm_init,& 
                            thlm, theta_type, rtm, um, vm, ugm, vgm, &
                            alt_type, press, subs_type, wm, &
                            sclrm, edsclrm, sounding_retVars, &
                            sclr_sounding_retVars )

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

    use interpolation, only:  & 
        lin_int ! Procedure(s)

    use array_index, only: & 
        iisclr_rt, &  ! Variable
        iisclr_thl
    !           ,iisclr_CO2

    use error_code, only: &
      clubb_at_least_debug_level ! Function

    use input_names, only: &
      z_name, &
      theta_name, &
      wm_name

    use input_reader, only: &
      one_dim_read_var

    implicit none

    ! External
    intrinsic :: trim, exp

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
      alt_type, &       ! Type of independent coordinate
      subs_type         ! Type of subsidence

    ! Output variables
    real, intent(out), dimension(gr%nnzp, sclr_dim) ::  & 
      sclrm   ! Passive scalar output      [units vary]

    real, intent(out), dimension(gr%nnzp, edsclr_dim) ::  & 
      edsclrm ! Eddy Passive scalar output [units vary]

    type(one_dim_read_var), dimension(n_snd_var), intent(out) :: &
      sounding_retVars ! Sounding Profile

    type(one_dim_read_var), dimension(sclr_dim), intent(out) :: &
      sclr_sounding_retVars ! Sclr Sounding Profile

    ! Local variables

    ! Input variables from namelist
    integer :: nlevels  ! Levels in the input sounding

    logical :: &
      l_sounding_exists        = .false., &
      l_sclr_sounding_exists   = .false., &
      l_edsclr_sounding_exists = .false.


    real, dimension(nmaxsnd) :: & 
    z,      & ! Altitude                               [m]
    theta,  & ! Liquid potential temperature sounding  [K]
    rt,     & ! Total water mixing ratio sounding      [kg/kg]
    u,      & ! u wind sounding                        [m/s] 
    v,      & ! v wind sounding                        [m/s]
    ug,     & ! u geostrophic wind sounding            [m/s]
    vg,     & ! v geostrophic wind sounding            [m/s]
    p_in_Pa, &   ! Pressure                            [Pa]
    subs      ! Subsidence                             [m/s or Pa/s]

    real, dimension(nmaxsnd, sclr_max) ::  & 
    sclr, edsclr ! Passive scalar input sounding    [units vary]

    integer :: i, j, k  ! Loop indices


    theta_type = theta_name ! Default value
    alt_type = z_name ! Default value
    subs_type = wm_name ! Defuault Value

    ! Determine which files exist ahead of time to allow for a graceful exit if
    ! one is missing.
    inquire(file="../input/case_setups/"//trim( runtype )//"_sounding.in", exist=l_sounding_exists)

    inquire(file="../input/case_setups/"//trim( runtype )//"_sclr_sounding.in", &
      exist=l_sclr_sounding_exists)

    inquire(file="../input/case_setups/"//trim( runtype )//"_edsclr_sounding.in", &
      exist=l_edsclr_sounding_exists)

    !---------------------------------------------------------------------------------------------
    ! Status Message
    if( clubb_at_least_debug_level(1) ) then
      print *, "Path to sounding: ", trim( runtype )//'_sounding.in'
      print *, "File exists? ", l_sounding_exists
      print *, "Path to sclr_sounding: ", trim( runtype )//'_sclr_sounding.in'
      print *, "File exists? ", l_sclr_sounding_exists
      print *, "Path to sounding: ", trim( runtype )//'_edsclr_sounding.in'
      print *, "File exists? ", l_edsclr_sounding_exists
    end if
    !----------------------------------------------------------------------------------------------

    if( l_sounding_exists ) then
      ! Read in SAM-Like <runtype>_sounding.in file
      call read_sounding_file( iunit, runtype, nlevels, psfc, zm_init, & 
                               z, theta, theta_type, rt, u, v, ug, vg, &
                               alt_type, p_in_Pa, subs_type, subs, &
                               sounding_retVars )
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
        if( l_sclr_sounding_exists ) then
          call read_sclr_sounding_file( iunit, runtype, sclr, &
          sclr_sounding_retVars )
        else
          stop 'Cannot open <runtype>_sclr_sounding.in file'
        end if
      end if
      if( edsclr_dim > 0 ) then
        if( l_edsclr_sounding_exists  ) then
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
      um   = u(1)
      vm   = v(1)
      ugm  = ug(1)
      vgm  = vg(1)
      thlm = theta(1)
      rtm  = rt(1)
      press = p_in_Pa(1)
      wm = subs(1)
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
      do while ( z(k) < gr%zt(i) )
        k=k+1
        if ( k > nlevels ) then
              write(fstderr,*) 'STOP Not enough sounding data to ',&
               'initialize grid:'
              write(fstderr,'(a,f7.1,/a,f7.1)') &
               '  highest sounding level', z(nlevels),&
               '  should be higher than highest thermodynamic point',&
               gr%zt(gr%nnzp)
              stop 'STOP in sounding'
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

    end do   ! i=2, gr%nnzp

    return
  end subroutine read_sounding

  !-------------------------------------------------------------------------------------------------
  subroutine read_sounding_file( iunit, runtype, nlevels, psfc, zm_init, &
                                 z, theta, theta_type, rt, u, v, ug, vg, &
                                 alt_type, p_in_Pa, subs_type, subs, retVars )
    !
    !  Description: This subroutine reads in a <runtype>_sounding.in file and
    !  returns the values contained in that file.
    !
    !-----------------------------------------------------------------------------------------------
    use input_reader, only: read_one_dim_file, read_x_profile, fill_blanks_one_dim_vars, &
                            one_dim_read_var

    use input_interpret, only: &
     read_z_profile, &
     read_theta_profile, &
     read_subs_profile

    use input_names, only: &
      rt_name, &
      um_name, &
      vm_name, &
      ug_name, &
      vg_name

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

    type(one_dim_read_var), intent(out), dimension(n_snd_var) :: &
      retVars ! Structure containing sounding profile

    character(len=*), intent(out) :: & 
      theta_type, &     ! Type of temperature sounding
      alt_type, &       ! Type of independent coordinate
      subs_type         ! Type of subsidence



    call read_one_dim_file( iunit, n_snd_var, &
    '../input/case_setups/'//trim( runtype )//'_sounding.in', retVars )

    call fill_blanks_one_dim_vars( n_snd_var, retVars )

    call read_z_profile( n_snd_var, nmaxsnd, retVars, psfc, zm_init, z, p_in_Pa, alt_type )

    call read_theta_profile( n_snd_var, nmaxsnd, retVars, theta_type, theta )

    rt = read_x_profile( n_snd_var, nmaxsnd, rt_name, retVars )

    u = read_x_profile( n_snd_var, nmaxsnd, um_name, retVars )

    v = read_x_profile( n_snd_var, nmaxsnd, vm_name, retVars )

    ug = read_x_profile( n_snd_var, nmaxsnd, ug_name, retVars )

    vg = read_x_profile( n_snd_var, nmaxsnd, vg_name, retVars )

    call read_subs_profile( n_snd_var, nmaxsnd, retVars, subs_type, subs )

    nlevels = size( retVars(1)%values )


  end subroutine read_sounding_file

  !-------------------------------------------------------------------------------------------------
  subroutine read_sclr_sounding_file( iunit, runtype, sclr, retVars )
    !
    !  Description: This subroutine reads in a <runtype>_sclr_sounding.in file and
    !  returns the values contained in that file.
    !
    !-----------------------------------------------------------------------------------------------
    use input_reader, only: read_one_dim_file, one_dim_read_var

    use parameters_model, only: sclr_dim

    use array_index, only: iisclr_rt, iisclr_thl, iisclr_CO2

    use input_names, only: &
      CO2_name, &
      rt_name, &
      theta_name, &
      thetal_name, &
      temperature_name

    implicit none

    ! Input Variable(s)
    integer, intent(in) :: iunit ! I/O unit

    character(len=*), intent(in) :: runtype ! String identifying the model case;
    !                                         e.g. bomex

    ! Output Variable(s)
    real, intent(inout), dimension(nmaxsnd,sclr_max) :: & 
      sclr        ! Scalar sounding [?]


    type(one_dim_read_var), dimension(sclr_dim), intent(out) :: &
      retVars ! Structure containing scalar sounding

    integer i

    call read_one_dim_file( iunit, sclr_dim, &
    '../input/case_setups/'//trim( runtype )//'_sclr_sounding.in', retVars )

!    call fill_blanks_one_dim_vars( sclr_dim, retVars )

    do i=1, sclr_dim
      select case ( trim( retVars(i)%name ) )
      case( CO2_name )
        if( i /= iisclr_CO2 .and. iisclr_CO2 > 0) then
          stop "iisclr_CO2 index does not match column."
        end if
      case ( rt_name )
        if( i /= iisclr_rt .and. iisclr_rt > 0) then
          stop "iisclr_rt index does not match column."
        end if
      case ( theta_name, thetal_name, temperature_name )
        if( i /= iisclr_thl .and. iisclr_thl > 0) then
          stop "iisclr_thl index does not match column."
        end if
      end select
      sclr(1:size(retVars(i)%values),i) = retVars(i)%values
    end do


  end subroutine read_sclr_sounding_file

  !-------------------------------------------------------------------------------------------------
  subroutine read_edsclr_sounding_file( iunit, runtype, edsclr )
    !
    !  Description: This subroutine reads in a <runtype>_edsclr_sounding.in file and
    !  returns the values contained in that file.
    !
    !-----------------------------------------------------------------------------------------------
    use input_reader, only: read_one_dim_file, &
                            one_dim_read_var, deallocate_one_dim_vars

    use parameters_model, only: edsclr_dim

    use array_index, only: iiedsclr_rt, iiedsclr_thl, iiedsclr_CO2

    use input_names, only: &
    CO2_name, &
    rt_name, &
    theta_name, &
    thetal_name, &
    temperature_name

    implicit none

    ! External
    intrinsic :: size, trim

    ! Input Variable(s)
    integer, intent(in) :: iunit ! I/O unit

    character(len=*), intent(in) :: runtype ! String identifying the model case;
    !                                         e.g. bomex

    ! Output Variable(s)
    real, intent(inout), dimension(nmaxsnd,sclr_max) :: & 
    edsclr ! Eddy Scalars [?]

    type(one_dim_read_var), dimension(edsclr_dim) :: retVars

    integer i

    call read_one_dim_file( iunit, edsclr_dim, &
    '../input/case_setups/'//trim( runtype )//'_edsclr_sounding.in', retVars )

 !   call fill_blanks_one_dim_vars( edsclr_dim, retVars )

    do i=1, edsclr_dim

      select case ( trim( retVars(i)%name ) )

      case( CO2_name )
        if( i /= iiedsclr_CO2 .and. iiedsclr_CO2 > 0) then
          stop "iisclr_CO2 index does not match column."
        end if
      case( rt_name )
        if( i /= iiedsclr_rt .and. iiedsclr_rt > 0) then
          stop "iisclr_rt index does not match column."
        end if
      case( theta_name, thetal_name, temperature_name )
        if( i /= iiedsclr_thl .and. iiedsclr_thl > 0) then
          stop "iisclr_thl index does not match column."
        end if
      end select
      edsclr(1:size( retVars(i)%values ),i) = retVars(i)%values
    end do

    call deallocate_one_dim_vars( edsclr_dim, retVars )

  end subroutine read_edsclr_sounding_file

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

    ! External
    intrinsic :: trim

    ! Constant Parameter
    integer, parameter :: nmaxsnd = 200

    ! Input Variables
    character(len=*), intent(in) :: fname

    ! Output Variable
    real, dimension(gr%nnzp), intent(out) :: x

    ! Local variables

    ! Input variables from namelist

    integer :: nlevels

    real, dimension(nmaxsnd) :: z, var

    ! Loop indices
    integer :: i,k

    namelist /profile/ nlevels, z, var

    ! ---- Begin Code ----

    ! Read sounding namelist

    open(10, file = trim( fname ), status = 'old' )
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
