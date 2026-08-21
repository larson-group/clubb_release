module readiopdata_module

    use readiopdata_core_module, only: &
        readiopdata_core_dims, readiopdata_core_snd, readiopdata_core_sfc, readiopdata_core_frc

    use netcdf, only: &
        NF90_OPEN, &
        NF90_NOWRITE, &
        NF90_NOERR, &
        NF90_CLOSE, &
        NF90_INQ_VARID, &
        NF90_INQ_DIMID, &
        NF90_GET_VAR, &
        NF90_MAX_NAME, &
        NF90_MAX_VAR_DIMS, &
        NF90_REAL, &
        NF90_FLOAT, &
        NF90_DOUBLE, &
        NF90_INT, &
        NF90_INQUIRE_VARIABLE, &
        NF90_INQUIRE_DIMENSION

    use clubb_precision, only: core_rknd

    implicit none

    public :: readiopdata_dims, readiopdata_snd, readiopdata_sfc, readiopdata_frc

    private

    ! Module level
    logical :: masterproc = .true.
    logical :: scamiop_from_global_cam = .false.
    character(len=120) :: iopfilepath = "../input/case_setups/BOMEX_5day_4scam.nc"

    ! Loading flags
    logical :: l_dims = .false.
    logical :: l_snd = .false.
    logical :: l_sfc = .false.
    logical :: l_frc = .false.

    ! Subroutine variables
    real, allocatable :: dplevs(:)
    integer :: ntime, nlev
    integer, allocatable :: tsec(:)
    logical :: get_add_surface_data, have_omega

    ! Repeatedly accessed variables
    real, allocatable :: dtls(:,:), dqls(:,:)

contains

subroutine readiopdata_dims( lat_val, lon_val )
    implicit none

    ! Outputs
    real( kind=core_rknd ), intent(out) :: lat_val, lon_val

    ! Locals
    real, allocatable :: lat_in(:), lon_in(:)
    integer :: bdate, nlat, nlon

    call readiopdata_core_dims( masterproc, iopfilepath, scamiop_from_global_cam, & ! Intent(in)
                                ntime, nlev, nlat, nlon, tsec,                    & ! Intent(out)
                                bdate, dplevs, lat_in, lon_in )                     ! Intent(out)

    lat_val = real(lat_in(1), core_rknd)
    lon_val = real(lon_in(1), core_rknd)

    l_dims = .true.

    deallocate( lat_in, lon_in )
end subroutine readiopdata_dims

! This needs to return only one timestep for the soundings
subroutine readiopdata_snd( zm_init, saturation_formula, &
                            nlevels, u, v, thlm, rtm, p, psfc, z, ug, vg, subs, subs_type )
    use constants_clubb, only: Rd, Cp, Lv, Ls, kappa, p0, zero_threshold, one, ep2, &
                                g_per_kg, pascal_per_mb
    use hydrostatic_module, only: inverse_hydrostatic
    use calc_pressure, only: calculate_thvm
    use saturation, only: &
        sat_mixrat_liq_api

    implicit none

    ! Needed for pressure -> altitude conversion
    real(kind=core_rknd), intent(in) :: zm_init

    integer, intent(in) :: &
      saturation_formula

    integer, intent(out) :: nlevels
    ! Clubb uses fixed-size static arrays
    real(kind=core_rknd), intent(out) :: psfc, rtm(:), &
        thlm(:), p(:), u(:), v(:), z(:), ug(:), vg(:), subs(:)
    character(len=*), intent(out) :: subs_type

    ! Allocatable arrays for readiopdata_core
    real, allocatable :: usnd(:,:), vsnd(:,:), tsnd(:,:), qsnd(:,:), psnd(:,:), &
        zsnd(:,:), ugls(:,:), vgls(:,:), wgls(:,:)
        
    real(kind=core_rknd), allocatable :: theta(:), thvm(:,:), &
    thlm_temporary(:), rcm(:,:), exner(:)

    integer :: i, k   ! Array index

    logical :: wgls_holds_omega, wgls_holds_w, have_w

    ! Requires that readiopdata_dims has been called
    if ( .not. l_dims ) then
        error stop "readiopdata_dims was not called before readiopdata_snd"
    end if

    call readiopdata_core_snd( masterproc, iopfilepath,                            & ! Intent(in)
                             scamiop_from_global_cam, ntime, nlev,                 & ! Intent(in)
                             dplevs, real(Rd), real(Cp), real(Lv/Cp), real(Ls/Cp), & ! Intent(in)
                             zsnd, usnd, vsnd, tsnd, qsnd, psnd, ugls, vgls, wgls, & ! Intent(out)
                             get_add_surface_data, have_omega, wgls_holds_omega    & ! Intent(out)
#ifdef UWM_MISC
                             , have_w, wgls_holds_w, subs_type                     & ! Intent(out)
#endif
                             )

    if(get_add_surface_data) then
        nlevels = nlev+1
    else
        nlevels = nlev
    endif

    ! Only use the first timestep
    allocate(theta(nlevels), thlm_temporary(nlevels))
    u(1:nlevels) = usnd(1:nlevels,1)
    v(1:nlevels) = vsnd(1:nlevels,1)
    theta(1:nlevels) = tsnd(1:nlevels,1)
    rtm(1:nlevels) = qsnd(1:nlevels,1)
    z(1:nlevels) = zsnd(1:nlevels,1)
    ug(1:nlevels) = ugls(1:nlevels,1)
    vg(1:nlevels) = vgls(1:nlevels,1)
    subs(1:nlevels) = wgls(1:nlevels,1)

    ! Separate psfc from psnd
    psfc = psnd(1,1) * pascal_per_mb

    p(1:nlevels) = psnd(1:nlevels,1) * pascal_per_mb


    ! Calculate zsnd
    ! Look at subroutine read_z_profile
    allocate(exner(nlevels), rcm(nlevels,1), thvm(nlevels, 1))

    ! Calculate exner from pressure.
    do k = 1, nlevels
      exner(k) = ( p(k) / p0 )**kappa
    enddo

    rtm(1:nlevels) = rtm(1:nlevels) / g_per_kg

    ! Return thlm in terms of potential temperature
    thlm(1:nlevels) = theta(1:nlevels)

    ! Determine initial cloud water mixing ratio at sounding levels
    ! based on potential temperature, exner, and rtm.

    do k = 1, 1
        do i = 1, nlevels
            rcm(i, k) &
                = max( rtm(i) &
                        - sat_mixrat_liq_api( p(i), &
                        theta(i) * exner(i), &
                        saturation_formula ), &
                        zero_threshold )
        enddo
    enddo

    ! Calculate theta_l from theta and cloud water mixing ratio, such
    ! that:  theta_l = theta - [Lv/(Cp*exner)]*rcm.
    thlm_temporary(1:nlevels) = theta(1:nlevels) - Lv/(Cp*exner(1:nlevels)) * rcm(1:nlevels, 1)


    call calculate_thvm( nlevels, 1,                                 & ! Intent(in)
                        reshape(thlm_temporary, [nlevels, 1]),       & ! Intent(in)
                        reshape(rtm, [nlevels, 1]),                  & ! Intent(in)
                        rcm,                                         & ! Intent(in)
                        reshape(exner, [nlevels, 1]),                & ! Intent(in)
                        reshape(theta, [nlevels, 1])                 & ! Intent(in)
                        * ( one + ep2 * ( reshape(rtm, [nlevels, 1]) & ! Intent(in)
                                        - rcm ) )**kappa,            & ! Intent(in)
                        thvm )                                         ! Intent(out)
! Find the altitudes, z, of the pressure sounding levels.
    call inverse_hydrostatic( psfc, zm_init, nlevels, thvm, exner, & ! Intent(in)
                              z )                                    ! Intent(out)

    l_snd = .true.

    deallocate( usnd, vsnd, tsnd, qsnd, psnd, zsnd, ugls, vgls, wgls, theta, exner, rcm, thvm )
end subroutine readiopdata_snd

subroutine readiopdata_sfc( tsfc, ubar, shsfc, lhsfc, timesfc )
    implicit none

    real(kind=core_rknd), optional, allocatable, intent(out) :: tsfc(:), ubar(:), &
        shsfc(:), lhsfc(:), timesfc(:)

    real, allocatable :: sstsfc(:), tausfc(:), shsfc_(:), lhsfc_(:)

    ! Requires that readiopdata_dims and readiopdata_snd have been called
    if ( .not. (l_dims .and. l_snd) ) then
        error stop "readiopdata_dims or readiopdata_snd was not called before readiopdata_sfc"
    end if

    call readiopdata_core_sfc( masterproc, iopfilepath,        & ! Intent(in)
                             scamiop_from_global_cam, .false., & ! Intent(in)
                             get_add_surface_data, ntime,      & ! Intent(in)
                             sstsfc, shsfc_, lhsfc_, tausfc )    ! Intent(out)

    if(present(tsfc)) then
        allocate( tsfc(ntime) )
        tsfc(:) = sstsfc(:)
    end if
    if(present(ubar)) then
        allocate( ubar(ntime) )
        ubar(:) = tausfc(:)
    end if
    if(present(shsfc)) then
        allocate( shsfc(ntime) )
        shsfc(:) = shsfc_(:)
    end if
    if(present(lhsfc)) then
        allocate( lhsfc(ntime) )
        lhsfc(:) = lhsfc_(:)
    end if
    if(present(timesfc)) then
        allocate( timesfc(ntime) )
        timesfc(:) = real(tsec(:), kind=core_rknd)
    end if

    l_sfc = .true.

    deallocate(shsfc_, lhsfc_, sstsfc, tausfc)
end subroutine readiopdata_sfc

subroutine readiopdata_frc( rtm, gr, thlm_forcing, rtm_forcing )
    use spec_hum_to_mixing_ratio, only: &
        force_spec_hum_to_mixing_ratio ! Procedure(s)

    use grid_class, only: &
        grid     ! Type

    implicit none

    type (grid), target, intent(in) :: gr

        ! Input Variable
    real( kind = core_rknd ), intent(in), dimension(gr%nzt) :: &
      rtm    ! Total water mixing ratio (thermodynamic levels) 

    real(kind=core_rknd), intent(out) ::  thlm_forcing(:), rtm_forcing(:)

    integer :: nlevels

    integer :: k

        ! Local Variables
    real( kind = core_rknd ), dimension(gr%nzt) :: &
      qtm_forcing  ! Specified total water spec. humidity tendency    [kg/kg/s]


    ! Requires that readiopdata_dims and readiopdata_snd have been called
    if ( .not. (l_dims .and. l_snd) ) then
        error stop "readiopdata_dims or readiopdata_snd was not called before readiopdata_frc"
    end if

    ! Only read from the file once
    if ( .not. l_frc ) then
        call readiopdata_core_frc( masterproc, iopfilepath, get_add_surface_data, & ! Intent(in)
                                 have_omega, dplevs, nlev, ntime,                 & ! Intent(in)
                                 dtls, dqls )                                       ! Intent(out)

        if ( get_add_surface_data ) then
            nlevels = nlev + 1
        else
            nlevels = nlev
        end if
    end if

    thlm_forcing(1:nlevels) = dtls(1:nlevels,1)
    rtm_forcing(1:nlevels) = dqls(1:nlevels,1)

    l_frc = .true.
end subroutine readiopdata_frc


end module readiopdata_module

! Mimic subroutine from sam_clubb
! This subroutine is needed by readiopdata_core
subroutine task_abort()
    error stop 999
end subroutine task_abort