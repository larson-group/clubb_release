! $Id$
module sounding

  implicit none

  public ::  & 
    read_sounding, & 
    read_profile ! Not currently used in CLUBB

  private :: read_sounding_file, read_sclr_sounding_file, &
    read_edsclr_sounding_file


  integer, private, parameter :: &
    nmaxsnd   = 600, & ! Constant parameters
    n_snd_var = 8

  integer, public, parameter :: &
    sclr_max  = 1000 ! Maximum number of scalars

  private ! Default Scope

  contains
  !------------------------------------------------------------------------
  subroutine read_sounding( gr, iunit, runtype, p_sfc, zm_init,& 
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

    use parameters_model, only: & 
        sclr_dim, &! Variable(s)
        edsclr_dim

    use interpolation, only:  & 
        lin_interpolate_two_points, & ! Procedure(s)
        mono_cubic_interp, &
        binary_search

    use array_index, only: & 
        iisclr_rt, &  ! Variable
        iisclr_thl
    !           ,iisclr_CO2

    use error_code, only: &
      clubb_at_least_debug_level  ! Procedure

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

    use extended_atmosphere_module, only: &
      convert_snd2extended_atm, & ! Procedure(s)
      load_extended_std_atm

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    use grid_class, only: grid

    implicit none

    type(grid), target, intent(in) :: gr

    ! External
    intrinsic :: trim, exp

    ! Constant parameter
    integer, parameter :: nmaxsnd = 10000

    ! Input variables
    integer, intent(in) :: iunit ! File unit to use for namelist

    character(len=*), intent(in) ::  & 
      runtype ! Used to determine if this in a DYCOMS II RF02 simulation

    real( kind = core_rknd ), intent(in) :: &
      p_sfc, & ! Pressure at the surface [Pa]
      zm_init ! Height at zm(1)         [m]

    ! Output variables
    real( kind = core_rknd ), intent(out), dimension(gr%nz) ::  & 
      thlm,  & ! Liquid water potential temperature    [K]
      rtm,   & ! Total water mixing ratio              [kg/kg]
      um,    & ! u wind component                      [m/s]
      vm,    & ! v wind component                      [m/s]
      ugm,   & ! u geostrophic wind component          [m/s]
      vgm,   & ! v geostrophic wind component          [m/s]
      press, & ! Pressure                              [Pa]
      wm       ! Subsidence                            [m/s or Pa/s]

    real( kind = core_rknd ), intent(out) ::  &
      rtm_sfc,  & ! Initial surface rtm                [kg/kg]
      thlm_sfc    ! Initial surface thlm               [K]

    character(len=*), intent(out) :: &
      theta_type, &     ! Type of temperature sounding
      alt_type, &       ! Type of independent coordinate
      subs_type         ! Type of subsidence

    real( kind = core_rknd ), intent(out), dimension(gr%nz, sclr_dim) ::  & 
      sclrm   ! Passive scalar output      [units vary]

    real( kind = core_rknd ), intent(out), dimension(gr%nz, edsclr_dim) ::  & 
      edsclrm ! Eddy Passive scalar output [units vary]

    ! Local variables

    ! Input variables from namelist
    integer :: nlevels  ! Levels in the input sounding

    ! Flag for interpolating the sounding profile with Steffen's monotone cubic 
    ! method to obtain smoother initial condition profile, which is found to be 
    ! benificial to achive a better numerical solution convergence. If this flage 
    ! is turned off, the initial conditions will be generated with linear interpolation. 
    logical, parameter :: l_modify_ic_with_cubic_int = .true. 

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
      sclr, edsclr ! Passive scalar input sounding    [units vary]

    type(one_dim_read_var), dimension(n_snd_var) :: &
      sounding_retVars ! Sounding Profile

    type(one_dim_read_var), dimension(sclr_dim) :: &
      sclr_sounding_retVars ! Sclr Sounding Profile

    integer :: i, j, k  ! Loop indices

    integer :: km1, kp1, kp2, k00 ! For mono cubic interpolation 

    integer :: idx  ! Result of binary search -- sounding level index.

    ! ---- Begin Code ----

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
    if( clubb_at_least_debug_level( 1 ) ) then
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
      call read_sounding_file( iunit, runtype, nlevels, p_sfc, zm_init, & 
                               z, theta, theta_type, rt, u, v, ug, vg, &
                               alt_type, p_in_Pa, subs_type, subs, &
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
      allocate( sclr(nmaxsnd, sclr_max), edsclr(nmaxsnd, sclr_max) )
      ! Initialize to zero
      sclr   = 0.0_core_rknd
      edsclr = 0.0_core_rknd
      ! Read in SAM-Like <runtype>_sclr_sounding.in and
      !                  <runtype>_edsclr_sounding.in
      if( sclr_dim > 0 ) then
        if( l_sclr_sounding_exists ) then
          call read_sclr_sounding_file( iunit, runtype, sclr, &
                                        sclr_sounding_retVars )
        else
          error stop 'Cannot open <runtype>_sclr_sounding.in file'
        end if
      end if
      if( edsclr_dim > 0 ) then
        if( l_edsclr_sounding_exists  ) then
          call read_edsclr_sounding_file( iunit, runtype, edsclr )
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
    ! (gr%zt(1,2)) is lower than the lowest value from the input sounding,
    ! then the linear interpolation scheme will fail.

    if ( gr%zt(1,2) < z(1) ) then
      write(fstderr,*) "Lowest level of input sounding, z(1), must be",  &
      " below the first above-model-surface thermodynamic level, gr%zt(1,2)"
      write(fstderr,*) " First sounding level z(1) = ", z(1)
      write(fstderr,*) " First thermodynamic level gr%zt(1,2) = ", gr%zt(1,2)
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
    do i=1, gr%nz
      if ( i == 1 .and. gr%zt(1,i) < z(1) ) then
         ! Thermodynamic level 1 is below the model lower boundary.
         ! If it's below the lowest sounding level, just skip setting it,
         ! since it is already initialized above.
         cycle
      endif
      k=1
      do while ( z(k) < gr%zt(1,i) )
        k=k+1
        if ( k > nlevels ) then
              write(fstderr,*) 'STOP Not enough sounding data to ',&
               'initialize grid:'
              write(fstderr,'(a,f7.1,/a,f7.1)') &
               '  highest sounding level', z(nlevels),&
               '  should be higher than highest thermodynamic point',&
               gr%zt(1,gr%nz)
              error stop 'STOP in sounding'
          exit
        end if  ! k > nlevels

        ! situation w/ cubic int. (achieve better numerical solution convergence)
        if (l_modify_ic_with_cubic_int) then 
          !use Steffen's monotone cubic interpolation method to obtain
          !smoothing initial condition profile for convergence test 
          !note: vertical levels in sounding file need to be not too coarse
          if ( k == 1 ) then ! Extrapolation for the ghost point
            km1 = k
            k00 = 1
            kp1 = 2
            kp2 = 3
          else if ( k == 2 ) then
            km1 = 1
            kp1 = 2
            kp2 = 3
            k00 = 1
          else
            km1 = k-2
            kp1 = k
            kp2 = k+1
            k00 = k-1
            !if z(k) reaches at the top level in sounding profile,
            !then use the nearest levels for interpolation 
            if ( z(k) >= z(nlevels) ) then
              km1 = nlevels-2
              kp1 = nlevels
              kp2 = nlevels
              k00 = nlevels-1
            end if
          end if

          um(i)    = mono_cubic_interp( gr%zt(1,i), km1, k00, kp1, kp2, z(km1), z(k00), z(kp1), & 
                                        z(kp2), u(km1), u(k00), u(kp1), u(kp2) )
          vm(i)    = mono_cubic_interp( gr%zt(1,i), km1, k00, kp1, kp2, z(km1), z(k00), z(kp1), &
                                        z(kp2), v(km1), v(k00), v(kp1), v(kp2) )
          ugm(i)   = mono_cubic_interp( gr%zt(1,i), km1, k00, kp1, kp2, z(km1), z(k00), z(kp1), &
                                        z(kp2), ug(km1), ug(k00), ug(kp1), ug(kp2) )
          vgm(i)   = mono_cubic_interp( gr%zt(1,i), km1, k00, kp1, kp2, z(km1), z(k00), z(kp1), &
                                        z(kp2), vg(km1), vg(k00), vg(kp1), vg(kp2) )
          thlm(i)  = mono_cubic_interp( gr%zt(1,i), km1, k00, kp1, kp2, z(km1), z(k00), z(kp1), &
                                        z(kp2), theta(km1), theta(k00), theta(kp1), theta(kp2) )
          rtm(i)   = mono_cubic_interp( gr%zt(1,i), km1, k00, kp1, kp2, z(km1), z(k00), z(kp1), &
                                        z(kp2), rt(km1), rt(k00), rt(kp1), rt(kp2) )
          press(i) = mono_cubic_interp( gr%zt(1,i), km1, k00, kp1, kp2, z(km1), z(k00), z(kp1), &
                                        z(kp2), p_in_Pa(km1), p_in_Pa(k00), p_in_Pa(kp1), p_in_Pa(kp2) )
          wm(i)    = mono_cubic_interp( gr%zt(1,i), km1, k00, kp1, kp2, z(km1), z(k00), z(kp1), &
                                        z(kp2), subs(km1), subs(k00), subs(kp1), subs(kp2) )

          if ( trim( runtype ) /= "dycoms2_rf02" ) then 
            !initial condition for tracers 
            if ( sclr_dim > 0 ) then
              do j = 1, sclr_dim
                sclrm(i,j) = mono_cubic_interp( gr%zt(1,i), km1, k00, kp1, kp2, & 
                                                z(km1), z(k00), z(kp1), z(kp2), & 
                                                sclr(km1,j), sclr(k00,j), &
                                                sclr(kp1,j), sclr(kp2,j) )
              end do
            end if
            if ( edsclr_dim > 0 ) then
              do j = 1, edsclr_dim
                edsclrm(i,j) = mono_cubic_interp( gr%zt(1,i), km1, k00, kp1, kp2, & 
                                                  z(km1), z(k00), z(kp1), z(kp2), & 
                                                  edsclr(km1,j), edsclr(k00,j), &
                                                  edsclr(kp1,j), edsclr(kp2,j) )
              end do
            end if
          else 
            ! DYCOMS2_RF02 case (use the same treatment as in regular situation w/ linear int.) 
            ugm(i)  = um(i)
            vgm(i)  = vm(i)
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
          end if  

        else ! default model setup 

          ! Regular situation w/ linear int.
          IF ( trim( runtype ) /= "dycoms2_rf02" ) THEN

            um(i)   = lin_interpolate_two_points( gr%zt(1,i), z(k), z(k-1), u(k), u(k-1) )
            vm(i)   = lin_interpolate_two_points( gr%zt(1,i), z(k), z(k-1), v(k), v(k-1) )
            ugm(i)  = lin_interpolate_two_points( gr%zt(1,i), z(k), z(k-1), ug(k), ug(k-1) )
            vgm(i)  = lin_interpolate_two_points( gr%zt(1,i), z(k), z(k-1), vg(k), vg(k-1) )
            thlm(i) = lin_interpolate_two_points( gr%zt(1,i), z(k), z(k-1), theta(k), theta(k-1) )
            rtm(i)  = lin_interpolate_two_points( gr%zt(1,i), z(k), z(k-1), rt(k), rt(k-1) )
            press(i) = lin_interpolate_two_points( gr%zt(1,i), z(k), z(k-1), p_in_Pa(k), p_in_Pa(k-1) )
            wm(i) = lin_interpolate_two_points( gr%zt(1,i), z(k), z(k-1), subs(k), subs(k-1) )

            if ( sclr_dim > 0 ) then
              do j = 1, sclr_dim
                sclrm(i,j) = lin_interpolate_two_points( gr%zt(1,i), z(k), z(k-1),  & 
                                      sclr(k,j), sclr(k-1,j) )
              end do
            end if
            if ( edsclr_dim > 0 ) then
              do j = 1, edsclr_dim
                edsclrm(i,j) = lin_interpolate_two_points( gr%zt(1,i), z(k), z(k-1),  & 
                                        edsclr(k,j), edsclr(k-1,j) )
              end do
            end if

          ELSE  ! DYCOMS II RF02 case
  
            IF ( gr%zt(1,i) < 795.0_core_rknd ) THEN
              ! (Wyant, et al. 2007, eq 1--4)
              um(i)   =  3.0_core_rknd + (4.3_core_rknd*gr%zt(1,i))/ &
                  1000.0_core_rknd ! Known magic number
              vm(i)   = -9.0_core_rknd + (5.6_core_rknd*gr%zt(1,i))/ &
                  1000.0_core_rknd ! Known magic number
              ugm(i)  = um(i)
              vgm(i)  = vm(i)
              thlm(i) = 288.3_core_rknd
              rtm(i)  = (9.45_core_rknd)/g_per_kg ! Known magic number
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
              press(i) = lin_interpolate_two_points( gr%zt(1,i), z(k), z(k-1), &
                                                     p_in_Pa(k), p_in_Pa(k-1) )
              wm(i) = lin_interpolate_two_points( gr%zt(1,i), z(k), z(k-1), subs(k), subs(k-1) )
            ELSE
              ! (Wyant, et al. 2007, eq 1--4)
              um(i)   =  3.0_core_rknd + (4.3_core_rknd*gr%zt(1,i))/ &
                            1000.0_core_rknd ! Known magic number
              vm(i)   = -9.0_core_rknd + (5.6_core_rknd*gr%zt(1,i))/ &
                            1000.0_core_rknd ! Known magic number
              ugm(i)  = um(i)
              vgm(i)  = vm(i)
              thlm(i) = 295.0_core_rknd + ( (gr%zt(1,i) - 795.0_core_rknd)** &
                          (1.0_core_rknd/3.0_core_rknd) ) ! Known magic number
              rtm(i)  = (  5.0_core_rknd - 3.0_core_rknd  & 
              * ( 1.0_core_rknd - EXP( (795.0_core_rknd - gr%zt(1,i))/ &
              500.0_core_rknd ) )  )/g_per_kg ! Known magic number
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
              press(i) = lin_interpolate_two_points( gr%zt(1,i), z(k), z(k-1), &
                                                     p_in_Pa(k), p_in_Pa(k-1) )
              wm(i) = lin_interpolate_two_points( gr%zt(1,i), z(k), z(k-1), subs(k), subs(k-1) )
            END IF
       
          END IF ! runtype

        end if ! l_modify_ic_with_cubic_int 

      end do ! do while ( z(k) < gr%zt(1,i) )

    end do   ! i=1, gr%nz


    ! The sounding will be initialized to thermodynamic grid levels successfully
    ! as long as the thermodynamic level 2 (at altitude gr%zt(1,2)) is at or above
    ! the lowest sounding level (z(1)).  However, it is advantageous to know the
    ! initial surface values of a few variables, as long as the sounding extends
    ! to the surface, which is found at momentum level 1 (at altitude gr%zm(1,1)).
    if ( gr%zm(1,1) < z(1) ) then

       ! The surface (or model lower boundary) is below the lowest sounding
       ! level.  Initialize the values of rtm_sfc and thlm_sfc to negative
       ! values that will be overwritten later.
       rtm_sfc  = -999.0_core_rknd
       thlm_sfc = -999.0_core_rknd

    else ! gr%zm(1,1) >= z(1)

       ! The surface (or model lower boundary) is above the lowest sounding
       ! level.  Use linear interpolation to find the values of rtm_sfc and
       ! thlm_sfc.

       ! Perform a binary search to find the two sounding levels that the
       ! surface (gr%zm(1,1)) is found between.  The value returned (idx) is the
       ! index of the closest value greater than or equal to gr%zm(1,1).
       idx = binary_search( nlevels, z, gr%zm(1,1) )

       ! The surface is found between sounding levels idx and idx-1.  Find the
       ! value of rtm_sfc.
       rtm_sfc = lin_interpolate_two_points( gr%zm(1,1), z(idx), z(idx-1), rt(idx), rt(idx-1) )

       ! The surface is found between sounding levels idx and idx-1.  Find the
       ! value of thlm_sfc.
       thlm_sfc &
          = lin_interpolate_two_points( gr%zm(1,1), z(idx), z(idx-1), theta(idx), theta(idx-1) )

    end if

    if ( rad_scheme == "bugsrad" ) then
      ! Prepare extended sounding for radiation
      if ( l_use_default_std_atmosphere ) then

        call load_extended_std_atm( iunit ) ! Intent(in)

      else

        call convert_snd2extended_atm( iunit, runtype, n_snd_var, p_sfc, zm_init, & ! Intent(in)
                                   sounding_retVars )   ! Intent(in)
      end if
    end if ! rad_scheme == "bugsrad"
    ! Deallocate sounding and scalar sounding profiles.  If this doesn't happen,
    ! then we'll have a memory leak.
    call deallocate_one_dim_vars( n_snd_var, sounding_retVars )
    call deallocate_one_dim_vars( sclr_dim, sclr_sounding_retVars )

    ! Deallocate sclr and edsclr arrays, iff allocated
    if ( allocated(sclr) ) then
      deallocate( sclr )
    end if

    if ( allocated(edsclr) ) then
      deallocate( edsclr )
    end if

    return
  end subroutine read_sounding

  !-----------------------------------------------------------------------------
  subroutine read_sounding_file( iunit, runtype, nlevels, p_sfc, zm_init, &
                                 z, theta, theta_type, rt, u, v, ug, vg, &
                                 alt_type, p_in_Pa, subs_type, subs, retVars )
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
  subroutine read_sclr_sounding_file( iunit, runtype, sclr, retVars )
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

    use parameters_model, only: sclr_dim ! Variable

    use array_index, only: iisclr_rt, iisclr_thl, iisclr_CO2

    use input_names, only: &
      CO2_name, &
      rt_name, &
      theta_name, &
      thetal_name, &
      temperature_name

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! Input Variable(s)
    integer, intent(in) :: iunit ! I/O unit

    character(len=*), intent(in) :: runtype ! String identifying the model case;
    !                                         e.g. bomex

    ! Output Variable(s)
    real( kind = core_rknd ), intent(inout), dimension(nmaxsnd,sclr_max) :: & 
      sclr        ! Scalar sounding [?]


    type(one_dim_read_var), dimension(sclr_dim), intent(out) :: &
      retVars ! Structure containing scalar sounding

    integer :: i

    call read_one_dim_file( iunit, sclr_dim, &
      '../input/case_setups/'//trim( runtype )//'_sclr_sounding.in', retVars )

!    call fill_blanks_one_dim_vars( sclr_dim, retVars )

    do i=1, sclr_dim
      select case ( trim( retVars(i)%name ) )
      case( CO2_name )
        if( i /= iisclr_CO2 .and. iisclr_CO2 > 0) then
          error stop "iisclr_CO2 index does not match column."
        end if
      case ( rt_name )
        if( i /= iisclr_rt .and. iisclr_rt > 0) then
          error stop "iisclr_rt index does not match column."
        end if
      case ( theta_name, thetal_name, temperature_name )
        if( i /= iisclr_thl .and. iisclr_thl > 0) then
          error stop "iisclr_thl index does not match column."
        end if
      end select
      sclr(1:size(retVars(i)%values),i) = retVars(i)%values
    end do

    return
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

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: size, trim

    ! Input Variable(s)
    integer, intent(in) :: iunit ! I/O unit

    character(len=*), intent(in) :: &
      runtype ! String identifying the model case; e.g. bomex

    ! Output Variable(s)
    real( kind = core_rknd ), intent(inout), dimension(nmaxsnd,sclr_max) :: & 
    edsclr ! Eddy Scalars [?]

    type(one_dim_read_var), dimension(edsclr_dim) :: retVars

    integer :: i 

    ! ---- Begin Code ----

    call read_one_dim_file( iunit, edsclr_dim, &
      '../input/case_setups/'//trim( runtype )//'_edsclr_sounding.in', retVars )

 !   call fill_blanks_one_dim_vars( edsclr_dim, retVars )

    do i=1, edsclr_dim

      select case ( trim( retVars(i)%name ) )

      case( CO2_name )
        if( i /= iiedsclr_CO2 .and. iiedsclr_CO2 > 0) then
          error stop "iisclr_CO2 index does not match column."
        end if
      case( rt_name )
        if( i /= iiedsclr_rt .and. iiedsclr_rt > 0) then
          error stop "iisclr_rt index does not match column."
        end if
      case( theta_name, thetal_name, temperature_name )
        if( i /= iiedsclr_thl .and. iiedsclr_thl > 0) then
          error stop "iisclr_thl index does not match column."
        end if
      end select
      edsclr(1:size( retVars(i)%values ),i) = retVars(i)%values
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

    type(grid), target, intent(in) :: gr

    ! External
    intrinsic :: trim

    ! Constant Parameter
    integer, parameter :: nmaxsnd = 200

    ! Input Variables
    character(len=*), intent(in) :: fname

    ! Output Variable
    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: x

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
    ! (gr%zt(1,2)) is lower than the lowest value from the input sounding,
    ! then the linear interpolation scheme will fail.

    if ( gr%zt(1,2) < z(1) ) then
      write(fstderr,*) "Lowest level of input sounding, z(1), must be",  &
      " below the first above-model-surface thermodynamic level, gr%zt(1,2)"
      write(fstderr,*) " First sounding level z(1) = ", z(1)
      write(fstderr,*) " First thermodynamic level gr%zt(1,2) = ", gr%zt(1,2)
      error stop 'STOP in read_profile (sounding.F90)'
    endif

    ! Use linear interpolation from two nearest prescribed grid points
    ! (one above and one below) to initialize mean quantities in the model
    ! Modified 27 May 2005 -dschanen: eliminated the goto in favor of a do while( )

    do i = 2, gr%nz
      k = 1
      do while ( z(k) < gr%zt(1,i) )
        k = k + 1
        if ( k > nlevels ) then
          write(fstderr,*) 'STOP Not enough sounding data to ',  & 
                           'initialize grid:'
          write(fstderr,'(a,f7.1,/a,f7.1)') ' Highest sounding level',  &
               z(nlevels),  &
               'should be higher than highest thermodynamic point',  &
               gr%zt(1,gr%nz)
          write(*,*) ' Filename: ', fname
          error stop 'STOP in read_profile (sounding.F90)'
        endif
        x(i) = lin_interpolate_two_points( gr%zt(1,i), z(k), z(k-1), var(k), var(k-1) )
      enddo ! while
    enddo ! i=2, gr%nzzp

    return
  end subroutine read_profile

end module sounding
!-----------------------------------------------------------------------
