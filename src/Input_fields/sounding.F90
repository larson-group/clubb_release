! $Id$
module sounding

  implicit none

  public ::  & 
    read_sounding, & 
    read_profile ! Not currently used in CLUBB

  private :: read_sounding_file, read_sclr_sounding_file, &
    read_edsclr_sounding_file


  integer, private, parameter :: &
    n_snd_var = 8   ! Constant parameter

  integer, public, parameter :: &
    sclr_max  = 1000 ! Maximum number of scalars

  private ! Default Scope

  contains
  !------------------------------------------------------------------------
  subroutine read_sounding( ngrdcol, sclr_dim, edsclr_dim, sclr_idx, &
                            gr, iunit, runtype, p_sfc, zm_init, & 
                            saturation_formula, &
                            l_modify_ic_with_cubic_int, &
                            thlm, theta_type, rtm, um, vm, ugm, vgm, &
                            alt_type, press, subs_type, wm, &
                            rtm_sfc, thlm_sfc, sclrm, edsclrm )

    ! Description:
    !   Subroutine to initialize model variables from a namelist file
    ! References:
    !   None
    !------------------------------------------------------------------------

    use constants_clubb, only:  & 
      fstderr, & ! Constant
      fstdout, &
      g_per_kg

    use interpolation, only:  & 
      lin_interpolate_two_points, & ! Procedure(s)
      mono_cubic_interp, &
      binary_search

    use array_index, only: &
      sclr_idx_type

    use error_code, only: &
      clubb_at_least_debug_level_api  ! Procedure

    use input_names, only: &
      z_name, & ! Variables
      theta_name, &
      wm_name

    use input_reader, only: &
      one_dim_read_var ! Type

    use parameters_radiation, only: &
      l_use_default_std_atmosphere, rad_scheme ! Variable(s)

    use input_reader, only: &
      deallocate_one_dim_vars ! Procedure(s)

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    use grid_class, only: &
        grid

    implicit none

    ! Constant parameter
    integer, parameter :: nmaxsnd = 10000

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
      ngrdcol, &
      sclr_dim, & 
      edsclr_dim

    type (sclr_idx_type), intent(in) :: &
      sclr_idx

    type(grid), intent(in) :: &
      gr

    integer, intent(in) :: iunit ! File unit to use for namelist

    character(len=*), intent(in) ::  & 
      runtype ! Used to determine if this in a DYCOMS II RF02 simulation

    real( kind = core_rknd ), intent(in) :: &
      p_sfc, & ! Pressure at the surface [Pa]
      zm_init ! Height at zm(1)         [m]

    integer, intent(in) :: &
      saturation_formula ! Integer that stores the saturation formula to be used

    ! Flag for interpolating the sounding profile with Steffen's monotone cubic 
    ! method to obtain smoother initial condition profile, which is found to be 
    ! beneficial to achive a better numerical solution convergence. If this flag 
    ! is turned off, the initial conditions will be generated with linear interpolation.
    ! This is done on a case-by-case basis, since using the monotone cubic method
    ! requires a special sounding.in file with many additional sounding levels.
    logical, intent(in) :: &
      l_modify_ic_with_cubic_int

    !--------------------- Output Variables ---------------------
    real( kind = core_rknd ), intent(out), dimension(ngrdcol,gr%nzt) ::  & 
      thlm,  & ! Liquid water potential temperature    [K]
      rtm,   & ! Total water mixing ratio              [kg/kg]
      um,    & ! u wind component                      [m/s]
      vm,    & ! v wind component                      [m/s]
      ugm,   & ! u geostrophic wind component          [m/s]
      vgm,   & ! v geostrophic wind component          [m/s]
      press, & ! Pressure                              [Pa]
      wm       ! Subsidence                            [m/s or Pa/s]

    real( kind = core_rknd ), dimension(ngrdcol), intent(out) ::  &
      rtm_sfc,  & ! Initial surface rtm                [kg/kg]
      thlm_sfc    ! Initial surface thlm               [K]

    character(len=*), intent(out) :: &
      theta_type, &     ! Type of temperature sounding
      alt_type, &       ! Type of independent coordinate
      subs_type         ! Type of subsidence

    real( kind = core_rknd ), intent(out), dimension(ngrdcol, gr%nzt, sclr_dim) ::  & 
      sclrm   ! Passive scalar output      [units vary]

    real( kind = core_rknd ), intent(out), dimension(ngrdcol, gr%nzt, edsclr_dim) ::  & 
      edsclrm ! Eddy Passive scalar output [units vary]

    !--------------------- Local Variables ---------------------

    ! Input variables from namelist
    integer :: nlevels  ! Levels in the input sounding

    logical :: &
      l_sounding_exists, &
      l_sclr_sounding_exists, &
      l_edsclr_sounding_exists

    real( kind = core_rknd ), dimension(nmaxsnd) :: & 
      z,       & ! Altitude                               [m]
      theta,   & ! Liquid potential temperature sounding  [K]
      rt,      & ! Total water mixing ratio sounding      [kg/kg]
      u,       & ! u wind sounding                        [m/s] 
      v,       & ! v wind sounding                        [m/s]
      ug,      & ! u geostrophic wind sounding            [m/s]
      vg,      & ! v geostrophic wind sounding            [m/s]
      p_in_Pa, & ! Pressure                               [Pa]
      subs       ! Vertical velocity sounding             [m/s or Pa/s]

    real( kind = core_rknd ), dimension(:,:), allocatable ::  & 
      sclr_snd, edsclr_snd ! Passive scalar input sounding    [units vary]

    type(one_dim_read_var), dimension(n_snd_var) :: &
      sounding_retVars ! Sounding Profile

    type(one_dim_read_var), dimension(sclr_dim) :: &
      sclr_sounding_retVars ! sclr_snd Sounding Profile

    integer :: i, sclr, edsclr, k, k_above  ! Loop indices

    integer :: km1, kp1, kp2, k00 ! For mono cubic interpolation

    integer :: idx  ! Result of binary search -- sounding level index.

    !--------------------- Begin Code ---------------------

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
    if( clubb_at_least_debug_level_api( 1 ) ) then
      print *, "Path to sounding: ", trim( runtype )//'_sounding.in'
      print *, "File exists? ", l_sounding_exists
      print *, "Path to sclr_sounding: ", trim( runtype )//'_sclr_sounding.in'
      print *, "File exists? ", l_sclr_sounding_exists
      print *, "Path to edsclr_sounding: ", trim( runtype )//'_edsclr_sounding.in'
      print *, "File exists? ", l_edsclr_sounding_exists
    end if
    !----------------------------------------------------------------------------------------------

    if( l_sounding_exists ) then
      ! Read in SAM-Like <runtype>_sounding.in file
      call read_sounding_file( iunit, runtype, nlevels, nmaxsnd, p_sfc, &
                               saturation_formula, &
                               zm_init, z, theta, theta_type, rt, u, v, &
                               ug, vg, alt_type, p_in_Pa, subs_type, subs, &
                               sounding_retVars )
    else
      write(fstderr,*) 'Cannot open ' // trim( runtype ) // '_sounding.in file'
      error stop 'Fatal error in read_sounding'
      ! sounding namelist is no longer used.
      ! Joshua Fasching April 2009
    end if

    ! Read in a passive scalar sounding, if enabled
    if ( sclr_dim > 0 .or. edsclr_dim > 0 ) then

      ! Allocate large arrays
      allocate( sclr_snd(nmaxsnd, sclr_max), edsclr_snd(nmaxsnd, sclr_max) )

      ! Initialize to zero
      sclr_snd   = 0.0_core_rknd
      edsclr_snd = 0.0_core_rknd

      ! Read in SAM-Like <runtype>_sclr_sounding.in and
      !                  <runtype>_edsclr_sounding.in
      if( sclr_dim > 0 ) then
        if( l_sclr_sounding_exists ) then
          call read_sclr_sounding_file( sclr_dim, sclr_idx, iunit, runtype, nmaxsnd, &
                                        sclr_snd, &
                                        sclr_sounding_retVars )
        else
          error stop 'Cannot open <runtype>_sclr_sounding.in file'
        end if

      end if

      if( edsclr_dim > 0 ) then
        if( l_edsclr_sounding_exists  ) then
          call read_edsclr_sounding_file( edsclr_dim, sclr_idx, iunit, runtype, nmaxsnd, &
                                          edsclr_snd )
        else
          error stop 'Cannot open <runtype>_edsclr_sounding.in file'
        end if
      end if

    end if

    if ( nlevels > nmaxsnd ) then
      write(fstderr,*) 'Error in sounding: nlevels > nmaxsnd'
      write(fstderr,*) 'nlevels = ',nlevels
      write(fstderr,*) 'nmaxsnd = ',nmaxsnd
      error stop 'STOP in read_sounding'
    end if

    ! Error check: if lowest above-model-surface themodynamic grid height
    ! (gr%zt(1,1)) is lower than the lowest value from the input sounding,
    ! then the linear interpolation scheme will fail.

    if ( gr%zt(1,1) < z(1) ) then
      write(fstderr,*) "Lowest level of input sounding, z(1), must be",  &
      " below the first above-model-surface thermodynamic level, gr%zt(1,1)"
      write(fstderr,*) " First sounding level z(1) = ", z(1)
      write(fstderr,*) " First thermodynamic level gr%zt(1,1) = ", gr%zt(1,1)
      error stop 'STOP in read_sounding'
    endif

    ! First sounding level should be near ground value

    ! dschanen 1 May 2007
    ! We have changed this for Nov. 11 and June 25, both of which
    ! begin above the ground.
    ! if ( abs(z(1)) > 1.e-8_core_rknd ) then
    if ( .false. ) then
      write(fstderr,*) 'First level of input sounding must be z=0'
      error stop 'STOP in read_sounding'
    else
      um(:,:)     = u(1)
      vm(:,:)     = v(1)
      ugm(:,:)    = ug(1)
      vgm(:,:)    = vg(1)
      thlm(:,:)   = theta(1)
      rtm(:,:)    = rt(1)
      press(:,:)  = p_in_Pa(1)
      wm(:,:)     = subs(1)

      if ( sclr_dim > 0 ) then
        do i = 1, ngrdcol
          sclrm(i,1,1:sclr_dim)   = sclr_snd(1,1:sclr_dim)
        end do
      end if

      if ( edsclr_dim > 0 ) then
        do i = 1, ngrdcol
          edsclrm(i,1,1:edsclr_dim) = edsclr_snd(1,1:edsclr_dim)
        end do
      end if

    end if

    if ( clubb_at_least_debug_level_api( 1 ) ) then

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

      do sclr = 1, sclr_dim, 1
        write(fstdout,'(a5,i2,a2)',advance='no') "sclr_snd(", sclr,") = "
        write(fstdout,'(8g10.3)') sclr_snd(1:nlevels,sclr)
      enddo

      do edsclr = 1, edsclr_dim, 1
        write(fstdout,'(a7,i2,a2)',advance='no') "edsclr_snd(", edsclr, ") = "
        write(fstdout,'(8g10.3)') edsclr_snd(1:nlevels,edsclr)
      enddo

    endif ! clubb_at_least_debug_level_api( 1 )
    !----------------------------------------------------------------------
    ! Use linear interpolation from two nearest prescribed grid points
    ! (one above and one below) to initialize mean quantities in the model
    ! Modified 27 May 2005 -dschanen: eliminated the goto in favor of a do while( )
    do i = 1, ngrdcol
      do k=1, gr%nzt
        k_above=1
        do while ( z(k_above) < gr%zt(i,k) )
          k_above=k_above+1
          if ( k_above > nlevels ) then
                write(fstderr,*) 'STOP Not enough sounding data to ',&
                'initialize grid:'
                write(fstderr,'(a,f7.1,/a,f7.1)') &
                '  highest sounding level', z(nlevels),&
                '  should be higher than highest thermodynamic point',&
                gr%zt(i,gr%nzt)
                error stop 'STOP in sounding'
            exit
          end if  ! k_above > nlevels

          ! situation w/ cubic int. (achieve better numerical solution convergence)
          if (l_modify_ic_with_cubic_int) then 
            !use Steffen's monotone cubic interpolation method to obtain
            !smoothing initial condition profile for convergence test 
            !note: vertical levels in sounding file need to be not too coarse
            if ( k_above == 1 ) then
              km1 = k_above
              k00 = 1
              kp1 = 2
              kp2 = 3
            else if ( k_above == 2 ) then
              km1 = 1
              kp1 = 2
              kp2 = 3
              k00 = 1
            else
              km1 = k_above-2
              kp1 = k_above
              kp2 = k_above+1
              k00 = k_above-1
              !if z(k_above) reaches at the top level in sounding profile,
              !then use the nearest levels for interpolation 
              if ( z(k_above) >= z(nlevels) ) then
                km1 = nlevels-2
                kp1 = nlevels
                kp2 = nlevels
                k00 = nlevels-1
              end if
            end if

            um(i,k)    = mono_cubic_interp( gr%zt(i,k), km1, k00, kp1, kp2, z(km1), z(k00), z(kp1), & 
                                            z(kp2), u(km1), u(k00), u(kp1), u(kp2) )
            vm(i,k)    = mono_cubic_interp( gr%zt(i,k), km1, k00, kp1, kp2, z(km1), z(k00), z(kp1), &
                                            z(kp2), v(km1), v(k00), v(kp1), v(kp2) )
            ugm(i,k)   = mono_cubic_interp( gr%zt(i,k), km1, k00, kp1, kp2, z(km1), z(k00), z(kp1), &
                                            z(kp2), ug(km1), ug(k00), ug(kp1), ug(kp2) )
            vgm(i,k)   = mono_cubic_interp( gr%zt(i,k), km1, k00, kp1, kp2, z(km1), z(k00), z(kp1), &
                                            z(kp2), vg(km1), vg(k00), vg(kp1), vg(kp2) )
            thlm(i,k)  = mono_cubic_interp( gr%zt(i,k), km1, k00, kp1, kp2, z(km1), z(k00), z(kp1), &
                                            z(kp2), theta(km1), theta(k00), theta(kp1), theta(kp2) )
            rtm(i,k)   = mono_cubic_interp( gr%zt(i,k), km1, k00, kp1, kp2, z(km1), z(k00), z(kp1), &
                                            z(kp2), rt(km1), rt(k00), rt(kp1), rt(kp2) )
            press(i,k) = mono_cubic_interp( gr%zt(i,k), km1, k00, kp1, kp2, z(km1), z(k00), z(kp1), &
                                            z(kp2), p_in_Pa(km1), p_in_Pa(k00), p_in_Pa(kp1), &
                                          p_in_Pa(kp2) )
            wm(i,k)    = mono_cubic_interp( gr%zt(i,k), km1, k00, kp1, kp2, z(km1), z(k00), z(kp1), &
                                            z(kp2), subs(km1), subs(k00), subs(kp1), subs(kp2) )

            if ( trim( runtype ) /= "dycoms2_rf02" ) then 
              !initial condition for tracers 
              if ( sclr_dim > 0 ) then
                do sclr = 1, sclr_dim
                  sclrm(i,k,sclr) = mono_cubic_interp( gr%zt(i,k), km1, k00, kp1, kp2, & 
                                                       z(km1), z(k00), z(kp1), z(kp2), & 
                                                       sclr_snd(km1,sclr), sclr_snd(k00,sclr), &
                                                       sclr_snd(kp1,sclr), sclr_snd(kp2,sclr) )
                end do
              end if
              if ( edsclr_dim > 0 ) then
                do edsclr = 1, edsclr_dim
                  edsclrm(i,k,edsclr) = mono_cubic_interp( gr%zt(i,k), km1, k00, kp1, kp2, & 
                                                           z(km1), z(k00), z(kp1), z(kp2), & 
                                                           edsclr_snd(km1,edsclr), edsclr_snd(k00,edsclr), &
                                                           edsclr_snd(kp1,edsclr), edsclr_snd(kp2,edsclr) )
                end do
              end if
            else 
              ! DYCOMS2_RF02 case (use the same treatment as in regular situation w/ linear int.) 
              ugm(i,k)  = um(i,k)
              vgm(i,k)  = vm(i,k)
              ! Passive Scalars
              ! Change this if they are not equal to theta_l and rt in RF02
              if ( sclr_idx%iisclr_thl > 0  ) then
                sclrm(i,k, sclr_idx%iisclr_thl)   = thlm(i,k)
                edsclrm(i,k, sclr_idx%iisclr_thl) = thlm(i,k)
              end if
              if ( sclr_idx%iisclr_rt > 0 ) then
                sclrm(i,k, sclr_idx%iisclr_rt)    = rtm(i,k)
                edsclrm(i,k, sclr_idx%iisclr_rt)  = rtm(i,k)
              end if
            end if  

          else ! default model setup 

            ! Regular situation w/ linear int.
            IF ( trim( runtype ) /= "dycoms2_rf02" ) THEN

              um(i,k)   = lin_interpolate_two_points( gr%zt(i,k), z(k_above), z(k_above-1), u(k_above), u(k_above-1) )
              vm(i,k)   = lin_interpolate_two_points( gr%zt(i,k), z(k_above), z(k_above-1), v(k_above), v(k_above-1) )
              ugm(i,k)  = lin_interpolate_two_points( gr%zt(i,k), z(k_above), z(k_above-1), ug(k_above), ug(k_above-1) )
              vgm(i,k)  = lin_interpolate_two_points( gr%zt(i,k), z(k_above), z(k_above-1), vg(k_above), vg(k_above-1) )
              thlm(i,k) = lin_interpolate_two_points( gr%zt(i,k), z(k_above), z(k_above-1), theta(k_above), theta(k_above-1) )
              rtm(i,k)  = lin_interpolate_two_points( gr%zt(i,k), z(k_above), z(k_above-1), rt(k_above), rt(k_above-1) )
              press(i,k) = lin_interpolate_two_points( gr%zt(i,k), z(k_above), z(k_above-1), p_in_Pa(k_above), &
                                                    p_in_Pa(k_above-1) )
              wm(i,k) = lin_interpolate_two_points( gr%zt(i,k), z(k_above), z(k_above-1), subs(k_above), subs(k_above-1) )

              if ( sclr_dim > 0 ) then
                do sclr = 1, sclr_dim
                  sclrm(i,k,sclr) = lin_interpolate_two_points( gr%zt(i,k), z(k_above), z(k_above-1),  & 
                                        sclr_snd(k_above,sclr), sclr_snd(k_above-1,sclr) )
                end do
              end if
              if ( edsclr_dim > 0 ) then
                do edsclr = 1, edsclr_dim
                  edsclrm(i,k,edsclr) = lin_interpolate_two_points( gr%zt(i,k), z(k_above), z(k_above-1),  & 
                                          edsclr_snd(k_above,edsclr), edsclr_snd(k_above-1,edsclr) )
                end do
              end if

            ELSE  ! DYCOMS II RF02 case
    
              IF ( gr%zt(i,k) < 795.0_core_rknd ) THEN
                ! (Wyant, et al. 2007, eq 1--4)
                um(i,k)   =  3.0_core_rknd + (4.3_core_rknd*gr%zt(i,k))/ &
                    1000.0_core_rknd ! Known magic number
                vm(i,k)   = -9.0_core_rknd + (5.6_core_rknd*gr%zt(i,k))/ &
                    1000.0_core_rknd ! Known magic number
                ugm(i,k)  = um(i,k)
                vgm(i,k)  = vm(i,k)
                thlm(i,k) = 288.3_core_rknd
                rtm(i,k)  = (9.45_core_rknd)/g_per_kg ! Known magic number
                ! Passive Scalars
                ! Change this if they are not equal to theta_l and rt in RF02
                if ( sclr_idx%iisclr_thl > 0  ) then
                  sclrm(i,k, sclr_idx%iisclr_thl)   = thlm(i,k)
                  edsclrm(i,k, sclr_idx%iisclr_thl) = thlm(i,k)
                end if
                if ( sclr_idx%iisclr_rt > 0 ) then
                  sclrm(i,k, sclr_idx%iisclr_rt)    = rtm(i,k)
                  edsclrm(i,k, sclr_idx%iisclr_rt)  = rtm(i,k)
                end if
                press(i,k) = lin_interpolate_two_points( gr%zt(i,k), z(k_above), z(k_above-1), &
                                                      p_in_Pa(k_above), p_in_Pa(k_above-1) )
                wm(i,k) = lin_interpolate_two_points( gr%zt(i,k), z(k_above), z(k_above-1), subs(k_above), subs(k_above-1) )
              ELSE
                ! (Wyant, et al. 2007, eq 1--4)
                um(i,k)   =  3.0_core_rknd + (4.3_core_rknd*gr%zt(i,k))/ &
                              1000.0_core_rknd ! Known magic number
                vm(i,k)   = -9.0_core_rknd + (5.6_core_rknd*gr%zt(i,k))/ &
                              1000.0_core_rknd ! Known magic number
                ugm(i,k)  = um(i,k)
                vgm(i,k)  = vm(i,k)
                thlm(i,k) = 295.0_core_rknd + ( (gr%zt(i,k) - 795.0_core_rknd)** &
                            (1.0_core_rknd/3.0_core_rknd) ) ! Known magic number
                rtm(i,k)  = (  5.0_core_rknd - 3.0_core_rknd  & 
                * ( 1.0_core_rknd - EXP( (795.0_core_rknd - gr%zt(i,k))/ &
                500.0_core_rknd ) )  )/g_per_kg ! Known magic number
                ! Passive Scalars
                ! Same as above
                if ( sclr_idx%iisclr_thl > 0  ) then
                  sclrm(i,k, sclr_idx%iisclr_thl)   = thlm(i,k)
                  edsclrm(i,k, sclr_idx%iisclr_thl) = thlm(i,k)
                end if
                if ( sclr_idx%iisclr_rt > 0 ) then
                  sclrm(i,k, sclr_idx%iisclr_rt)    = rtm(i,k)
                  edsclrm(i,k, sclr_idx%iisclr_rt)  = rtm(i,k)
                end if
                press(i,k) = lin_interpolate_two_points( gr%zt(i,k), z(k_above), z(k_above-1), &
                                                      p_in_Pa(k_above), p_in_Pa(k_above-1) )
                wm(i,k) = lin_interpolate_two_points( gr%zt(i,k), z(k_above), z(k_above-1), subs(k_above), subs(k_above-1) )
              END IF
        
            END IF ! runtype

          end if ! l_modify_ic_with_cubic_int

        end do ! do while ( z(k) < gr%zt(i,k) )

      end do   ! k=1, gr%nzt
    end do


    ! The sounding will be initialized to thermodynamic grid levels successfully
    ! as long as the thermodynamic level 1 (at altitude gr%zt(1,1)) is at or above
    ! the lowest sounding level (z(1)).  However, it is advantageous to know the
    ! initial surface values of a few variables, as long as the sounding extends
    ! to the surface, which is found at momentum level 1 (at altitude gr%zm(1,1)).
    do i = 1, ngrdcol
      if ( gr%zm(i,1) < z(1) ) then

        ! The surface (or model lower boundary) is below the lowest sounding
        ! level.  Initialize the values of rtm_sfc and thlm_sfc to negative
        ! values that will be overwritten later.
        rtm_sfc(i)  = -999.0_core_rknd
        thlm_sfc(i) = -999.0_core_rknd

      else ! gr%zm(1,1) >= z(1)

        ! The surface (or model lower boundary) is above the lowest sounding
        ! level.  Use linear interpolation to find the values of rtm_sfc and
        ! thlm_sfc.

        ! Perform a binary search to find the two sounding levels that the
        ! surface (gr%zm(1,1)) is found between.  The value returned (idx) is the
        ! index of the closest value greater than or equal to gr%zm(1,1).
        idx = binary_search( nlevels, z, gr%zm(i,1) )

        ! The surface is found between sounding levels idx and idx-1.  Find the
        ! value of rtm_sfc.
        rtm_sfc(i) = lin_interpolate_two_points( gr%zm(i,1), z(idx), z(idx-1), &
                                                 rt(idx), rt(idx-1) )

        ! The surface is found between sounding levels idx and idx-1.  Find the
        ! value of thlm_sfc.
        thlm_sfc(i) = lin_interpolate_two_points( gr%zm(i,1), z(idx), z(idx-1), &
                                                  theta(idx), theta(idx-1) )

      end if
    end do

    if ( rad_scheme == "bugsrad" ) then
      ! Prepare extended sounding for radiation
      if ( l_use_default_std_atmosphere ) then

        call load_extended_std_atm( iunit ) ! Intent(in)

      else

        call convert_snd2extended_atm( iunit, runtype, n_snd_var, p_sfc, zm_init, & ! Intent(in)
                                       sounding_retVars, saturation_formula )   ! Intent(in)
      end if
    end if ! rad_scheme == "bugsrad"

    ! Deallocate sounding and scalar sounding profiles.  If this doesn't happen,
    ! then we'll have a memory leak.
    call deallocate_one_dim_vars( n_snd_var, sounding_retVars )
    call deallocate_one_dim_vars( sclr_dim, sclr_sounding_retVars )

    ! Deallocate sclr_snd and edsclr_snd arrays, iff allocated
    if ( allocated(sclr_snd) ) then
      deallocate( sclr_snd )
    end if

    if ( allocated(edsclr_snd) ) then
      deallocate( edsclr_snd )
    end if

    return
  end subroutine read_sounding

  !-----------------------------------------------------------------------------
  subroutine read_sounding_file( iunit, runtype, nlevels, nmaxsnd, p_sfc, &
                                 saturation_formula, &
                                 zm_init, z, theta, theta_type, rt, u, v, &
                                 ug, vg, alt_type, p_in_Pa, subs_type, subs, &
                                 retVars )
    !
    !  Description: This subroutine reads in a <runtype>_sounding.in file and
    !  returns the values contained in that file.

    !-----------------------------------------------------------------------

    use input_reader, only: &
      read_one_dim_file, read_x_profile, fill_blanks_one_dim_vars ! Procedure(s)

    use input_reader, only: &
      one_dim_read_var ! Type

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

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variable(s)
    integer, intent(in) :: iunit ! I/O unit

    character(len=*), intent(in) :: &
      runtype ! String identifying the model case; e.g. bomex

    real( kind = core_rknd ), intent(in) :: &
      p_sfc, & ! Pressure at the surface [Pa]
      zm_init ! Height at zm(1)         [m]

    integer, intent(in) :: &
      nmaxsnd   ! Maximum number of levels allowed in the sounding.in file

    integer, intent(in) :: &
      saturation_formula ! Integer that stores the saturation formula to be used

    ! Output Variable(s)
    integer, intent(out) :: nlevels ! Number of levels from the sounding.in file

    real( kind = core_rknd ), intent(out), dimension(nmaxsnd) :: & 
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

    ! ---- Begin Code ----

    call read_one_dim_file( iunit, n_snd_var, &
    '../input/case_setups/'//trim( runtype )//'_sounding.in', retVars )

    call fill_blanks_one_dim_vars( n_snd_var, retVars )

    call read_z_profile( n_snd_var, nmaxsnd, retVars, p_sfc, zm_init, &
                         saturation_formula, &
                         z, p_in_Pa, alt_type )

    call read_theta_profile( n_snd_var, nmaxsnd, retVars, theta_type, theta )

    rt = read_x_profile( n_snd_var, nmaxsnd, rt_name, retVars )

    u = read_x_profile( n_snd_var, nmaxsnd, um_name, retVars )

    v = read_x_profile( n_snd_var, nmaxsnd, vm_name, retVars )

    ug = read_x_profile( n_snd_var, nmaxsnd, ug_name, retVars )

    vg = read_x_profile( n_snd_var, nmaxsnd, vg_name, retVars )

    call read_subs_profile( n_snd_var, nmaxsnd, retVars, subs_type, subs )

    nlevels = size( retVars(1)%values )

    return
  end subroutine read_sounding_file

  !-------------------------------------------------------------------------------------------------
  subroutine read_sclr_sounding_file( sclr_dim, sclr_idx, iunit, runtype, nmaxsnd, &
                                      sclr_snd, &
                                      retVars )
    !
    ! Description: This subroutine reads in a <runtype>_sclr_sounding.in file and
    !   returns the values contained in that file.

    ! References:
    !   None
    !-----------------------------------------------------------------------------------------------
    use input_reader, only: &
      read_one_dim_file ! Procedure(s)

    use input_reader, only: &
      one_dim_read_var ! Type

    use array_index, only: &
      sclr_idx_type

    use input_names, only: &
      CO2_name, &
      rt_name, &
      theta_name, &
      thetal_name, &
      temperature_name

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
      sclr_dim

    type (sclr_idx_type), intent(in) :: &
      sclr_idx

    integer, intent(in) :: iunit ! I/O unit

    character(len=*), intent(in) :: runtype ! String identifying the model case;
    !                                         e.g. bomex

    integer, intent(in) :: &
      nmaxsnd   ! Maximum number of levels allowed in the sounding.in file

    !--------------------- InOut Variables ---------------------
    real( kind = core_rknd ), intent(inout), dimension(nmaxsnd,sclr_max) :: & 
      sclr_snd        ! Scalar sounding [?]

    !--------------------- Output Variables --------------------
    type(one_dim_read_var), dimension(sclr_dim), intent(out) :: &
      retVars ! Structure containing scalar sounding

    !--------------------- Local Variables --------------------
    integer :: sclr

    call read_one_dim_file( iunit, sclr_dim, &
      '../input/case_setups/'//trim( runtype )//'_sclr_sounding.in', retVars )

!    call fill_blanks_one_dim_vars( sclr_dim, retVars )

    do sclr=1, sclr_dim
      select case ( trim( retVars(sclr)%name ) )
      case( CO2_name )
        if( sclr /= sclr_idx%iisclr_CO2 .and. sclr_idx%iisclr_CO2 > 0) then
          error stop "iisclr_CO2 index does not match column."
        end if
      case ( rt_name )
        if( sclr /= sclr_idx%iisclr_rt .and. sclr_idx%iisclr_rt > 0) then
          error stop "iisclr_rt index does not match column."
        end if
      case ( theta_name, thetal_name, temperature_name )
        if( sclr /= sclr_idx%iisclr_thl .and. sclr_idx%iisclr_thl > 0) then
          error stop "iisclr_thl index does not match column."
        end if
      end select
      sclr_snd(1:size(retVars(sclr)%values),sclr) = retVars(sclr)%values
    end do

    return
  end subroutine read_sclr_sounding_file

  !-------------------------------------------------------------------------------------------------
  subroutine read_edsclr_sounding_file( edsclr_dim, sclr_idx, iunit, runtype, nmaxsnd, &
                                        edsclr_snd )
    !
    !  Description: This subroutine reads in a <runtype>_edsclr_sounding.in file and
    !  returns the values contained in that file.
    !
    !-----------------------------------------------------------------------------------------------
    use input_reader, only: read_one_dim_file, &
                            one_dim_read_var, deallocate_one_dim_vars

    use array_index, only: &
      sclr_idx_type

    use input_names, only: &
      CO2_name, &
      rt_name, &
      theta_name, &
      thetal_name, &
      temperature_name

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    !--------------------- Input Variables ---------------------
    integer, intent(in) :: &
      edsclr_dim

    type (sclr_idx_type), intent(in) :: &
      sclr_idx

    integer, intent(in) :: iunit ! I/O unit

    character(len=*), intent(in) :: &
      runtype ! String identifying the model case; e.g. bomex

    integer, intent(in) :: &
      nmaxsnd   ! Maximum number of levels allowed in the sounding.in file

    !--------------------- Output Variables ---------------------
    real( kind = core_rknd ), intent(inout), dimension(nmaxsnd,sclr_max) :: & 
    edsclr_snd ! Eddy Scalars [?]

    !--------------------- Local Variables ---------------------
    type(one_dim_read_var), dimension(edsclr_dim) :: retVars

    integer :: edsclr 

    !--------------------- Begin Code ---------------------

    call read_one_dim_file( iunit, edsclr_dim, &
      '../input/case_setups/'//trim( runtype )//'_edsclr_sounding.in', retVars )

 !   call fill_blanks_one_dim_vars( edsclr_dim, retVars )

    do edsclr=1, edsclr_dim

      select case ( trim( retVars(edsclr)%name ) )

      case( CO2_name )
        if( edsclr /= sclr_idx%iiedsclr_CO2 .and. sclr_idx%iiedsclr_CO2 > 0) then
          error stop "iisclr_CO2 index does not match column."
        end if
      case( rt_name )
        if( edsclr /= sclr_idx%iiedsclr_rt .and. sclr_idx%iiedsclr_rt > 0) then
          error stop "iisclr_rt index does not match column."
        end if
      case( theta_name, thetal_name, temperature_name )
        if( edsclr /= sclr_idx%iiedsclr_thl .and. sclr_idx%iiedsclr_thl > 0) then
          error stop "iisclr_thl index does not match column."
        end if
      end select
      edsclr_snd(1:size( retVars(edsclr)%values ),edsclr) = retVars(edsclr)%values
    end do

    call deallocate_one_dim_vars( edsclr_dim, retVars )

  end subroutine read_edsclr_sounding_file

  !------------------------------------------------------------------------
  subroutine read_profile( gr, fname, x )

    !       Description:
    !       Subroutine to initialize one generic model variable from file
    !------------------------------------------------------------------------

    use interpolation, only:  & 
        lin_interpolate_two_points ! Procedure

    use constants_clubb, only:  &
        fstderr ! Constant

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    use grid_class, only: grid ! Type

    implicit none

    type(grid), intent(in) :: gr

    ! External
    intrinsic :: trim

    ! Constant Parameter
    integer, parameter :: nmaxsnd = 200

    ! Input Variables
    character(len=*), intent(in) :: fname

    ! Output Variable
    real( kind = core_rknd ), dimension(gr%nzt), intent(out) :: x

    ! Local variables

    ! Input variables from namelist

    integer :: nlevels

    real( kind = core_rknd ), dimension(nmaxsnd) :: z, var

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
      error stop 'STOP in read_profile (sounding.F90)'
    endif

    ! Error check: if lowest above-model-surface themodynamic grid height
    ! (gr%zt(1,1)) is lower than the lowest value from the input sounding,
    ! then the linear interpolation scheme will fail.

    if ( gr%zt(1,1) < z(1) ) then
      write(fstderr,*) "Lowest level of input sounding, z(1), must be",  &
      " below the first above-model-surface thermodynamic level, gr%zt(1,1)"
      write(fstderr,*) " First sounding level z(1) = ", z(1)
      write(fstderr,*) " First thermodynamic level gr%zt(1,1) = ", gr%zt(1,1)
      error stop 'STOP in read_profile (sounding.F90)'
    endif

    ! Use linear interpolation from two nearest prescribed grid points
    ! (one above and one below) to initialize mean quantities in the model
    ! Modified 27 May 2005 -dschanen: eliminated the goto in favor of a do while( )

    do i = 1, gr%nzt
      k = 1
      do while ( z(k) < gr%zt(1,i) )
        k = k + 1
        if ( k > nlevels ) then
          write(fstderr,*) 'STOP Not enough sounding data to ',  & 
                           'initialize grid:'
          write(fstderr,'(a,f7.1,/a,f7.1)') ' Highest sounding level',  &
               z(nlevels),  &
               'should be higher than highest thermodynamic point',  &
               gr%zt(1,gr%nzt)
          write(*,*) ' Filename: ', fname
          error stop 'STOP in read_profile (sounding.F90)'
        endif
        x(i) = lin_interpolate_two_points( gr%zt(1,i), z(k), z(k-1), var(k), var(k-1) )
      enddo ! while
    enddo ! i = 1, gr%nzt

    return
  end subroutine read_profile

  !-------------------------------------------------------------------------------------------------
  subroutine convert_snd2extended_atm( iunit, runtype, n_snd_var, p_sfc, zm_init, &
                                       sounding_profiles, saturation_formula )
    !
    !  Description: This subroutine converts information retrieved from the
    !    sounding files of a case into a format usable for an extended atmosphere.
    !    The extended atmosphere profile is stored in module variables.
    !
    !  References:
    !    none
    !
    !-----------------------------------------------------------------------------------------------
    use input_reader, only: &
      read_x_profile,  & ! Procedure(s)
      read_one_dim_file

    use input_reader, only: &
      one_dim_read_var ! Type

    use input_interpret, only: read_theta_profile, read_z_profile ! Procedure(s)

    use constants_clubb, only: &
      kappa, & ! Constant(s)
      p0, &
      fstderr, &
      pascal_per_mb

    use input_names, only: z_name, temperature_name, ozone_name, rt_name ! Variable(s)

    use clubb_precision, only: &
      core_rknd

    use extended_atmosphere_module, only: &
      extended_atmos_dim, &   ! Size of Extended Atmosphere
      extended_alt, &         ! Altitude, increases with array index    [m]
      extended_T_in_K, &      ! Temperature in degrees Kelvin
      extended_sp_hmdty, &    ! Specific Humidity ( Water Vapor / Density )
      extended_p_in_mb, &     ! Pressure in millibars
      extended_o3l            ! Ozone ( O_3 / Density )

    implicit none

    ! External
    intrinsic :: size, trim

    ! Constant Parameters
    integer, parameter :: &
      n_rad_scalars = 1

    ! Input Variable(s)

    integer, intent(in) :: iunit  ! Fortran file unit

    character(len=*), intent(in) ::  runtype

    integer, intent(in) :: n_snd_var   ! Number of variables from sounding [-]

    real( kind = core_rknd ), intent(in) :: p_sfc  ! Pressure at the surface [Pa]

    real( kind = core_rknd ), intent(in) :: zm_init ! Height at zm(1) [m]

    type(one_dim_read_var), dimension(n_snd_var), intent(in) :: &
      sounding_profiles ! Sounding profile

    integer, intent(in) :: &
      saturation_formula ! Integer that stores the saturation formula to be used

    ! Local Variables

    real( kind = core_rknd ),  dimension(:), allocatable :: &
      extended_p_in_Pa, &
      extended_exner

    type(one_dim_read_var), dimension(n_rad_scalars) :: &
      extended_rad_scalars ! Ozone Sounding profile

    integer :: i, ivar

    character(len=20) :: alt_type, theta_type

    ! -- Begin Code --

    ! Determine the size of the extended atmosphere buffer
    extended_atmos_dim = size( sounding_profiles(1)%values )
   
    ! Allocate variables
    allocate( extended_alt(extended_atmos_dim) )
    allocate( extended_T_in_K(extended_atmos_dim) )
    allocate( extended_sp_hmdty(extended_atmos_dim) )
    allocate( extended_p_in_mb(extended_atmos_dim) )
    allocate( extended_o3l(extended_atmos_dim) )
    allocate( extended_p_in_Pa(extended_atmos_dim) )
    allocate( extended_exner(extended_atmos_dim) )

    ! Either convert to pressure or from pressure

    call read_z_profile( n_snd_var, extended_atmos_dim, sounding_profiles, p_sfc, zm_init, &
                         saturation_formula, &
                         extended_alt , extended_p_in_Pa, alt_type )

    if ( alt_type == z_name ) then
      write(fstderr,*) "Fatal error in convert_snd2extended_atm."
      error stop "Feature not implemented"
    end if

    extended_p_in_mb = extended_p_in_Pa / pascal_per_mb

    ! Convert to temperature from thlm or theta

    call read_theta_profile( n_snd_var, extended_atmos_dim, sounding_profiles, &
      theta_type, extended_T_in_K )

    if( theta_type /= temperature_name ) then
      extended_exner(1) = ( p_sfc/p0 ) ** kappa
      do i = 2, extended_atmos_dim
        extended_exner(i) = (extended_p_in_Pa(i)/p0) ** kappa
      end do
      extended_T_in_K = extended_T_in_K * extended_exner
    end if

    ! Convert rtm to specific humidity

    extended_sp_hmdty = read_x_profile( n_snd_var, extended_atmos_dim, rt_name, &
                                            sounding_profiles )
    extended_sp_hmdty = extended_sp_hmdty / ( extended_sp_hmdty +1._core_rknd )

    ! Read in radiation scalars sounding (currently it only holds O3)
    call read_one_dim_file( iunit, n_rad_scalars, & ! In
      '../input/case_setups/'//trim( runtype )//'_ozone_sounding.in', & ! In
      extended_rad_scalars ) ! Out

    ! Set the array holding the values of o3l (ozone)
    extended_o3l = read_x_profile( 1, extended_atmos_dim, ozone_name, &
                                       extended_rad_scalars )
    ! We would add the setting of new radiation scalar arrays like O3 right
    ! after this. New variable names would have to be added to input_names.

    ! Free Memory
    deallocate( extended_p_in_Pa )
    deallocate( extended_exner )

    do ivar = 1, n_rad_scalars
      deallocate( extended_rad_scalars(ivar)%values )
    end do

    return
  end subroutine convert_snd2extended_atm

  !-------------------------------------------------------------------------------------------------
  subroutine load_extended_std_atm ( iunit )
    !
    !  Description:
    !    Loads in the U.S. Standard atmosphere data from a file.
    !  References:
    !    McClatchey, et al., (1972) _Environmental Research Papers_,
    !    No. 411, p.94
    !
    !-----------------------------------------------------------------------------------------------

    use input_reader, only: &
      read_x_profile, & ! Procedure(s)
      read_one_dim_file, &
      deallocate_one_dim_vars

    use input_reader, only: &
      one_dim_read_var ! Derived type

    use input_names, only: &
      z_name, &
      ozone_name, &
      temperature_name, &
      press_mb_name, &
      sp_humidity_name

    use extended_atmosphere_module, only: &
      extended_atmos_dim, &   ! Size of Extended Atmosphere
      extended_alt, &         ! Altitude, increases with array index    [m]
      extended_T_in_K, &      ! Temperature in degrees Kelvin
      extended_sp_hmdty, &    ! Specific Humidity ( Water Vapor / Density )
      extended_p_in_mb, &     ! Pressure in millibars
      extended_o3l            ! Ozone ( O_3 / Density )
      
    implicit none

    ! Constant Parameters
    integer, parameter :: &
      nCol = 5,  &       ! Number of columns in the text file
      std_atmos_dim = 50 ! Dimension of the standard atmosphere table

    character(len=*), parameter :: &
      atm_input_file = "../input/std_atmosphere/atmosphere.in"

    ! Input Variable(s)
    integer, intent(in) :: iunit ! File I/O unit [-]


    ! Local Variable(s)
    type(one_dim_read_var), dimension(nCol) :: retVars

    ! -- Begin Code --

    extended_atmos_dim = std_atmos_dim

    call read_one_dim_file( iunit, nCol, atm_input_file, retVars )

    ! Allocate and initialize variables for standard atmosphere

    allocate( extended_alt(extended_atmos_dim) )
    allocate( extended_T_in_K(extended_atmos_dim) )
    allocate( extended_sp_hmdty(extended_atmos_dim) )
    allocate( extended_p_in_mb(extended_atmos_dim) )
    allocate( extended_o3l(extended_atmos_dim) )

    extended_alt = read_x_profile( nCol, extended_atmos_dim, z_name, retVars, &
                                 atm_input_file )

    extended_T_in_K = read_x_profile( nCol, extended_atmos_dim, temperature_name, retVars, &
                                    atm_input_file )

    extended_sp_hmdty = read_x_profile( nCol, extended_atmos_dim, sp_humidity_name, retVars, &
                                      atm_input_file )

    extended_p_in_mb = read_x_profile( nCol, extended_atmos_dim, press_mb_name, retVars, &
                                   atm_input_file )

    extended_o3l = read_x_profile( nCol, extended_atmos_dim, ozone_name, retVars, &
                                 atm_input_file )

    ! Deallocate memory
    call deallocate_one_dim_vars( nCol, retVars )

    return
  end subroutine load_extended_std_atm

end module sounding
!-----------------------------------------------------------------------
