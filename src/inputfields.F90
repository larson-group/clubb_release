!-----------------------------------------------------------------------
! $Id$

! Module inputfields

!  This exists because I wanted to keep the grads_reader code
!  generalized and bypass having to pass the datafile as a parameter,
!  since there may be situations where the fields are be calculated
!  analytically or by calling a different model without reading a datafile.
!  Using a module also saves the trouble of writing an interface definition
!  within the clubb_inputfields code.
!===============================================================================
module inputfields

  implicit none

!----- Run information--------------------------------------------------
  character(len=80), public :: datafile

  character(len=100), public ::  & 
  datafilet, datafilem

  character(3), public      :: input_type

  logical, public :: input_um, input_vm, input_rtm, input_thlm, & 
                     input_wp2, input_wprtp, input_wpthlp,  & 
                     input_wp3, input_rtp2, input_thlp2,  & 
                     input_rtpthlp, input_upwp, input_vpwp, & 
                     input_ug, input_vg, input_rcm,  & 
                     input_wm_zt, input_exner, input_em, & 
                     input_p, input_rho, input_rho_zm, & 
                     input_Lscale, input_Lscale_up, input_Lscale_down, & 
                     input_Kh_zt, input_Kh_zm, input_tau_zm, input_tau_zt, & 
                     input_wpthvp, &
                     input_thl1, input_thl2, input_a, input_s1, input_s2, &
                     input_ss1, input_ss2, input_rc1, input_rc2, &
                     input_thvm, input_rrainm, input_Nrm,  input_Ncm,  & 
                     input_rsnowm, input_ricem, input_rgraupelm, input_Ncnm, input_Nim, & 
                     input_thlm_forcing, input_rtm_forcing, & 
                     input_up2, input_vp2, input_sigma_sqd_w, & 
                     input_cf, input_sigma_sqd_w_zt, &
                     input_veg_T_in_K, input_deep_soil_T_in_K, &
                     input_sfc_soil_T_in_K


  public  :: grads_fields_reader, &
             compute_timestep, &
             set_filenames

  public :: CLUBB_levels_within_LES_domain, &
             LES_grid_to_CLUBB_grid, &
             lin_ext_zm_bottom, &
             lin_ext_zt_bottom

  private ! Default Scope

  contains

!===============================================================================


!-----------------------------------------------------------------------
  subroutine set_filenames( )
! Description: Set the names of the GrADS files to be used.
!   Used by clubb_inputfields and clubb_restart.
!-----------------------------------------------------------------------

    implicit none

    select case ( input_type )
    case ( "les", "rf1" )
      datafilet = trim( datafile )//"_coamps_sm.ctl"
      datafilem = trim( datafile )//"_coamps_sw.ctl"
    case ( "hoc" )
      datafilet = trim( datafile )//"_zt.ctl"
      datafilem = trim( datafile )//"_zm.ctl"
    case default
      write(0,*) "Don't know how to handle input_type = "// & 
        input_type
      stop
    end select

    return
  end subroutine set_filenames
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
  subroutine grads_fields_reader( timestep )
!       Description:
!       Reads in variables for the model from GrADS data

!       Calls:
!       subroutine open_grads_read
!       subroutine get_grads_var
!       subroutine close_grads_read
!-----------------------------------------------------------------------

    use variables_prognostic_module, only: & 
        um,  & ! Variable(s)
        vm, & 
        rtm, & 
        thlm, & 
        wp2, & 
        wp3, & 
        wprtp, & 
        wpthlp, & 
        wpthvp, & 
        rtp2, & 
        thlp2, & 
        rtpthlp, & 
        upwp, & 
        vpwp, & 
        Kh_zt, & 
        p_in_Pa, & 
        exner, & 
        rcm, & 
        wm_zt, & 
        rho, & 
        rho_zm, & 
        thlm_forcing, & 
        rtm_forcing, & 
        cf, & 
        tau_zm, & 
        up2, & 
        vp2, & 
        sigma_sqd_w

    use variables_diagnostic_module, only: & 
        hydromet,  & ! Variable(s)
        tau_zt, & 
        ug, & 
        vg, & 
        Lscale, & 
        Lscale_up, & 
        Lscale_down, & 
        Kh_zm, & 
        thvm, & 
        Ncnm, & 
        sigma_sqd_w_zt, & 
        em

    use variables_prognostic_module, only: & 
        pdf_params ! Variable(s)

    use grid_class, only: & 
        gr,  & ! Variable(s)
        zt2zm ! Procedure(s)

    use constants, only:  &
        rttol,    & ! Variable(s)
        thltol,   &
        wtol_sqd, &
        emin,     &
        fstderr

    use array_index, only:  & 
        iirrainm, iiNrm, iirsnowm, iiricem, iirgraupelm, iiNim, iiNcm

    use stat_file_module, only: & 
        stat_file     ! Type

    use inputfile_class, only: & 
        get_grads_var,  & ! Procedure(s)
        open_grads_read, & 
        close_grads_read

    use interpolation, only: &
        lin_int ! Procedure(s)

    use parameters_microphys, only: &
      micro_scheme ! Variable(s)

    use soil_vegetation, only: deep_soil_T_in_K, sfc_soil_T_in_K, veg_T_in_K

    implicit none

    ! External
    intrinsic :: max, trim, any

    ! Arguments
    integer, intent(in) :: timestep

    ! Local Variables
    logical :: l_read_error, l_fatal_error

    type (stat_file) :: fread_var

    real, dimension(:), allocatable :: LES_tmp1

    real, dimension(gr%nnzp+1) :: tmp1

    integer ::  &
      k_lowest_zt_input, &  ! The lowest CLUBB thermodynamic level that's within the LES domain.
      k_highest_zt_input, & ! The highest CLUBB thermodynamic level that's within the LES domain.
      k_lowest_zm_input, &  ! The lowest CLUBB momentum level that's within the LES domain.
      k_highest_zm_input    ! The highest CLUBB momentum level that's within the LES domain.

    ! Variables used to reconcile CLUBB thermodynamic levels with LES vertical levels.
    integer, dimension(:), allocatable ::  &
      exact_lev_idx_zt, & ! In case of an exact match, index of LES level that is
                        ! exactly even with CLUBB thermodynamic level k.
      lower_lev_idx_zt, & ! In case linear interpolation is needed, index of LES
                        ! level that is immediately below CLUBB thermo. level k.
      upper_lev_idx_zt  ! In case linear interpolation is needed, index of LES
                        ! level that is immediately above CLUBB thermo. level k.

    logical, dimension(:), allocatable ::  &
      l_lin_int_zt  ! Flag that is turned on if linear interpolation is needed.

    ! Variables used to reconcile CLUBB momentum levels with LES vertical levels.
    integer, dimension(:), allocatable ::  &
      exact_lev_idx_zm, & ! In case of an exact match, index of LES level that is
                        ! exactly even with CLUBB momentum level k.
      lower_lev_idx_zm, & ! In case linear interpolation is needed, index of LES
                        ! level that is immediately below CLUBB momentum level k
      upper_lev_idx_zm    ! In case linear interpolation is needed, index of LES
    ! level that is immediately above CLUBB momentum level k

    logical, dimension(:), allocatable ::  &
      l_lin_int_zm  ! Flag that is turned on if linear interpolation is needed.

    integer :: k  ! Array index


    select case( input_type )

    case( "hoc" )

      ! NOTE:  The code is not set up to compensate for grid discrepancies between
      !        the CLUBB GrADS file that is having its variable values passed in
      !        and the CLUBB inputfields or restart run that is using those values
      !        as variable inputs.  Therefore, CLUBB should be set up to match the
      !        number of grid levels and altitude of each grid level found in the
      !        CLUBB GrADS zt and zm files.

      !  Thermo grid - zt file
      call open_grads_read( 15, trim( datafile )//"_zt.ctl",  & 
                            fread_var )


      ! Initialize l_fatal_error for case ( "hoc" )
      l_fatal_error = .false.

      if ( input_um ) then
        call get_grads_var( fread_var, "um", timestep, & 
                      um(1:gr%nnzp), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error
      endif

      if ( input_vm ) then
        call get_grads_var( fread_var, "vm", timestep, & 
                      vm(1:gr%nnzp), l_read_error )
        l_fatal_error = l_fatal_error .or. l_read_error
      endif

      if ( input_rtm ) then
        call get_grads_var( fread_var, "rtm", timestep, & 
                      rtm(1:gr%nnzp),  l_read_error )
        l_fatal_error = l_fatal_error .or. l_read_error
      endif

      if ( input_thlm ) then
        call get_grads_var( fread_var, "thlm",  & 
                      timestep, & 
                      thlm(1:gr%nnzp),  l_read_error )
        l_fatal_error = l_fatal_error .or. l_read_error
      endif

      if ( input_wp3 ) then
        call get_grads_var( fread_var, "wp3", timestep, & 
                      wp3(1:gr%nnzp),  l_read_error )
        l_fatal_error = l_fatal_error .or. l_read_error
      endif
      if ( input_tau_zt ) then
        call get_grads_var( fread_var, "tau_zt", timestep, & 
                      tau_zt(1:gr%nnzp),  l_read_error )
        l_fatal_error = l_fatal_error .or. l_read_error
      endif
      if ( input_rrainm ) then
        call get_grads_var( fread_var, "rrainm", timestep, & 
                      hydromet(1:gr%nnzp,iirrainm),  l_read_error )
        l_fatal_error = l_fatal_error .or. l_read_error
      endif
      if ( input_rsnowm ) then
        call get_grads_var( fread_var, "rsnowm", timestep, & 
                      hydromet(1:gr%nnzp,iirsnowm),  l_read_error )
        l_fatal_error = l_fatal_error .or. l_read_error
      endif
      if ( input_ricem ) then
        call get_grads_var( fread_var, "ricem", timestep, & 
                      hydromet(1:gr%nnzp,iiricem),  l_read_error )
        l_fatal_error = l_fatal_error .or. l_read_error
      endif
      if ( input_rgraupelm ) then
        call get_grads_var( fread_var, "rgraupelm", timestep, & 
                      hydromet(1:gr%nnzp,iirgraupelm),  l_read_error )
        l_fatal_error = l_fatal_error .or. l_read_error
      endif

!--------------------------------------------------------
! Added variables for clubb_restart
      if ( input_p ) then
        call get_grads_var( fread_var, "p_in_Pa", timestep, & 
                      p_in_Pa(1:gr%nnzp),  l_read_error )
        l_fatal_error = l_fatal_error .or. l_read_error
      endif
      if ( input_exner) then
        call get_grads_var( fread_var , "exner", timestep, & 
                      exner(1:gr%nnzp),  l_read_error)
        l_fatal_error = l_fatal_error .or. l_read_error
      endif
      if ( input_ug) then
        call get_grads_var( fread_var , "ug", timestep, & 
                      ug(1:gr%nnzp),  l_read_error)
        l_fatal_error = l_fatal_error .or. l_read_error
      endif
      if ( input_vg) then
        call get_grads_var( fread_var , "vg", timestep, & 
                      vg(1:gr%nnzp),  l_read_error)
        l_fatal_error = l_fatal_error .or. l_read_error
      endif
      if ( input_rcm) then
        call get_grads_var( fread_var , "rcm", timestep, & 
                      rcm(1:gr%nnzp),  l_read_error)
        l_fatal_error = l_fatal_error .or. l_read_error
      endif
      if ( input_wm_zt) then
        call get_grads_var( fread_var , "wm", timestep, & 
                      wm_zt(1:gr%nnzp),  l_read_error)
        l_fatal_error = l_fatal_error .or. l_read_error
      endif
      if ( input_rho) then
        call get_grads_var( fread_var , "rho", timestep, & 
                      rho(1:gr%nnzp), l_read_error)
        l_fatal_error = l_fatal_error .or. l_read_error
      endif
      if ( input_Lscale) then
        call get_grads_var( fread_var , "Lscale", timestep, & 
                      Lscale(1:gr%nnzp), l_read_error)
        l_fatal_error = l_fatal_error .or. l_read_error
      endif
      if ( input_Lscale_up) then
        call get_grads_var( fread_var , "Lscale_up", timestep, & 
                      Lscale_up(1:gr%nnzp), l_read_error)
        l_fatal_error = l_fatal_error .or. l_read_error
      endif
      if ( input_Lscale_down) then
        call get_grads_var( fread_var , "Lscale_down", timestep, & 
                      Lscale_down(1:gr%nnzp), l_read_error)
        l_fatal_error = l_fatal_error .or. l_read_error
      endif
      if ( input_Kh_zt) then
        call get_grads_var( fread_var , "Kh_zt", timestep, & 
                      Kh_zt(1:gr%nnzp), l_read_error)
        l_fatal_error = l_fatal_error .or. l_read_error
      endif
      if ( input_thvm) then
        call get_grads_var( fread_var , "thvm", timestep, & 
                      thvm(1:gr%nnzp), l_read_error)
        l_fatal_error = l_fatal_error .or. l_read_error
      endif
      if ( input_thlm_forcing ) then
        call get_grads_var( fread_var , "thlm_f", timestep, & 
                      thlm_forcing(1:gr%nnzp), l_read_error)
        l_fatal_error = l_fatal_error .or. l_read_error
      endif
      if ( input_rtm_forcing ) then
        call get_grads_var( fread_var , "rtm_f", timestep, & 
                      rtm_forcing(1:gr%nnzp), l_read_error)
        l_fatal_error = l_fatal_error .or. l_read_error
      endif
      if ( input_Ncm) then
        call get_grads_var( fread_var , "Ncm", timestep, & 
                      hydromet(1:gr%nnzp,iiNcm), l_read_error)
        l_fatal_error = l_fatal_error .or. l_read_error
      endif
      if ( input_Ncnm) then
        call get_grads_var( fread_var , "Ncnm", timestep, & 
                      Ncnm(1:gr%nnzp), l_read_error)
        l_fatal_error = l_fatal_error .or. l_read_error
      endif
      if ( input_Nim) then
        call get_grads_var( fread_var , "Nim", timestep, & 
                      hydromet(1:gr%nnzp,iiNim), l_read_error)
        l_fatal_error = l_fatal_error .or. l_read_error
      endif
      if ( input_cf) then
        call get_grads_var( fread_var , "cf", timestep, & 
                      cf(1:gr%nnzp), l_read_error)
        l_fatal_error = l_fatal_error .or. l_read_error
      endif
      if ( input_Nrm ) then
        call get_grads_var( fread_var , "Nrm", timestep, & 
                      hydromet(1:gr%nnzp,iiNrm), l_read_error)
        l_fatal_error = l_fatal_error .or. l_read_error
      endif

      if ( input_sigma_sqd_w_zt ) then
        call get_grads_var( fread_var , "sigma_sqd_w_zt", timestep, & 
                      sigma_sqd_w_zt(1:gr%nnzp), l_read_error)
        l_fatal_error = l_fatal_error .or. l_read_error
      endif

      ! PDF Parameters (needed for K&K microphysics)
      if ( input_thl1 ) then
        call get_grads_var( fread_var , "thl1", timestep, & 
                      pdf_params%thl1(1:gr%nnzp), l_read_error )
        l_fatal_error = l_fatal_error .or. l_read_error
      end if

      if ( input_thl2 ) then
        call get_grads_var( fread_var , "thl2", timestep, & 
                      pdf_params%thl2(1:gr%nnzp), l_read_error )
        l_fatal_error = l_fatal_error .or. l_read_error
      end if

      if ( input_a ) then
        call get_grads_var( fread_var , "a", timestep, & 
                      pdf_params%a(1:gr%nnzp), l_read_error )
        l_fatal_error = l_fatal_error .or. l_read_error
      end if

      if ( input_s1 ) then
        call get_grads_var( fread_var , "s1", timestep, & 
                      pdf_params%s1(1:gr%nnzp), l_read_error )
        l_fatal_error = l_fatal_error .or. l_read_error
      end if

      if ( input_s2 ) then
        call get_grads_var( fread_var , "s2", timestep, & 
                      pdf_params%s2(1:gr%nnzp), l_read_error )
        l_fatal_error = l_fatal_error .or. l_read_error
      end if

      if ( input_ss1 ) then
        call get_grads_var( fread_var , "ss1", timestep, & 
                      pdf_params%ss1(1:gr%nnzp), l_read_error )
        l_fatal_error = l_fatal_error .or. l_read_error
      end if

      if ( input_ss2 ) then
        call get_grads_var( fread_var , "ss2", timestep, & 
                      pdf_params%ss2(1:gr%nnzp), l_read_error )
        l_fatal_error = l_fatal_error .or. l_read_error
      end if

      if ( input_rc1 ) then
        call get_grads_var( fread_var , "rc1", timestep, & 
                      pdf_params%rc1(1:gr%nnzp), l_read_error )
        l_fatal_error = l_fatal_error .or. l_read_error
      end if

      if ( input_rc2 ) then
        call get_grads_var( fread_var , "rc2", timestep, & 
                      pdf_params%rc2(1:gr%nnzp), l_read_error )
        l_fatal_error = l_fatal_error .or. l_read_error
      end if

!--------------------------------------------------------
      call close_grads_read( fread_var )

      ! Read in the zm file
      call open_grads_read( 15, trim(datafile)//"_zm.ctl", & 
                            fread_var )

      if ( input_wp2 ) then
        call get_grads_var( fread_var, "wp2", timestep, & 
                      wp2(1:gr%nnzp),  l_read_error )
        l_fatal_error = l_fatal_error .or. l_read_error
      endif

      if ( input_wprtp ) then
        call get_grads_var( fread_var, "wprtp",  & 
                      timestep, wprtp(1:gr%nnzp), & 
                      l_read_error )
        l_fatal_error = l_fatal_error .or. l_read_error
      endif

      if ( input_wpthlp ) then
        call get_grads_var( fread_var, "wpthlp",  & 
                      timestep,  & 
                      wpthlp(1:gr%nnzp),  & 
                      l_read_error )
        l_fatal_error = l_fatal_error .or. l_read_error
      endif

      if ( input_wpthvp ) then
        call get_grads_var( fread_var, "wpthvp",  & 
                      timestep,  & 
                      wpthvp(1:gr%nnzp),  & 
                      l_read_error )
        l_fatal_error = l_fatal_error .or. l_read_error
      endif

      if ( input_rtp2 ) then
        call get_grads_var( fread_var, "rtp2",  & 
                      timestep, & 
                      rtp2(1:gr%nnzp), l_read_error )
        l_fatal_error = l_fatal_error .or. l_read_error
      endif

      if ( input_thlp2 ) then
        call get_grads_var( fread_var, "thlp2",  & 
                      timestep, & 
                      thlp2(1:gr%nnzp), l_read_error )
        l_fatal_error = l_fatal_error .or. l_read_error
      endif

      if ( input_rtpthlp ) then
        call get_grads_var( fread_var, "rtpthlp",  & 
                      timestep,  & 
                      rtpthlp(1:gr%nnzp), & 
                      l_read_error )
        l_fatal_error = l_fatal_error .or. l_read_error
      endif

      if ( input_upwp ) then
        call get_grads_var( fread_var, "upwp",  & 
                      timestep, & 
                      upwp(1:gr%nnzp), l_read_error )
        l_fatal_error = l_fatal_error .or. l_read_error
      endif

      if ( input_vpwp ) then
        call get_grads_var( fread_var, "vpwp",  & 
                      timestep, & 
                      vpwp(1:gr%nnzp), l_read_error )
        l_fatal_error = l_fatal_error .or. l_read_error
      endif
!-----------------------------------------------------------
      if ( input_em ) then
        call get_grads_var( fread_var, "em", & 
                      timestep, & 
                      em(1:gr%nnzp), l_read_error )
        l_fatal_error = l_fatal_error .or. l_read_error
      endif
      if ( input_rho_zm ) then
        call get_grads_var( fread_var, "rho_zm", & 
                      timestep, & 
                      rho_zm(1:gr%nnzp), l_read_error )
        l_fatal_error = l_fatal_error .or. l_read_error
      endif
      if ( input_Kh_zm ) then
        call get_grads_var( fread_var, "Kh_zm", & 
                      timestep, & 
                      Kh_zm(1:gr%nnzp), l_read_error )
        l_fatal_error = l_fatal_error .or. l_read_error
      endif
      if ( input_tau_zm) then
        call get_grads_var( fread_var, "tau_zm", & 
                      timestep, & 
                      tau_zm(1:gr%nnzp), l_read_error )
        l_fatal_error = l_fatal_error .or. l_read_error
      endif
      if ( input_up2) then
        call get_grads_var( fread_var, "up2", & 
                      timestep, & 
                      up2(1:gr%nnzp), l_read_error )
        l_fatal_error = l_fatal_error .or. l_read_error
      endif
      if ( input_vp2) then
        call get_grads_var( fread_var, "vp2", & 
                      timestep, & 
                      vp2(1:gr%nnzp), l_read_error )
        l_fatal_error = l_fatal_error .or. l_read_error
      endif
      if ( input_sigma_sqd_w ) then
        call get_grads_var( fread_var, "sigma_sqd_w", & 
                      timestep, & 
                      sigma_sqd_w(1:gr%nnzp), l_read_error )
        l_fatal_error = l_fatal_error .or. l_read_error
      endif

!-----------------------------------------------------------
      call close_grads_read( fread_var )

      call open_grads_read( 15, trim( datafile )//"_sfc.ctl",  & 
                            fread_var )

      if ( input_veg_T_in_K ) then
        call get_grads_var( fread_var, "veg_T_in_K", & 
                      timestep, & 
                      tmp1, l_read_error )
        l_fatal_error = l_fatal_error .or. l_read_error
        veg_T_in_K = tmp1(1)
        print *, "Veg T = ", veg_T_in_K
      endif
      if ( input_deep_soil_T_in_K ) then
        call get_grads_var( fread_var, "deep_soil_T_in_", & 
                      timestep, & 
                      tmp1, l_read_error )
        l_fatal_error = l_fatal_error .or. l_read_error
        deep_soil_T_in_K = tmp1(1)
        print *,"Deep soil = ",deep_soil_T_in_K
      endif
      if ( input_sfc_soil_T_in_K ) then
        call get_grads_var( fread_var, "sfc_soil_T_in_K", & 
                      timestep, & 
                      tmp1, l_read_error )
        l_fatal_error = l_fatal_error .or. l_read_error
        sfc_soil_T_in_K = tmp1(1)
        print *,"surface_soil = ", sfc_soil_T_in_K
      endif


      if ( l_fatal_error ) stop "oops, get_grads_var failed in grads_fields_reader"

      call close_grads_read( fread_var )


    case ( "rf1" )   ! special case for COAMPS DYCOMS-II RF01

      ! stats_sm
      call open_grads_read( 15, trim(datafile)//"_coamps_sm.ctl",  & 
                            fread_var )

      ! Temporarily store LES output in variable array LES_tmp1.
      ! Allocate LES_tmp1 based on lowest and highest vertical indices of LES
      ! output.
      allocate( LES_tmp1(fread_var%ia:fread_var%iz) )

      ! Find the lowest and highest indices of CLUBB thermodynamic levels that
      ! fall within the domain of the LES output.
      call CLUBB_levels_within_LES_domain( fread_var, gr%zt,  &
                                           k_lowest_zt_input, k_highest_zt_input )

      allocate( exact_lev_idx_zt(k_lowest_zt_input:k_highest_zt_input) )
      allocate( lower_lev_idx_zt(k_lowest_zt_input:k_highest_zt_input) )
      allocate( upper_lev_idx_zt(k_lowest_zt_input:k_highest_zt_input) )
      allocate( l_lin_int_zt(k_lowest_zt_input:k_highest_zt_input) )

      ! For all CLUBB thermodynamic levels, k, that are within the LES domain,
      ! find either the index of the LES level that exactly matches the altitude
      ! of the CLUBB level, or find the two indices of the LES levels that are on
      ! either side of the CLUBB level.
      do k = k_lowest_zt_input, k_highest_zt_input, 1
        ! CLUBB vertical level k is found at an altitude that is within the
        ! domain of the LES output.
        call LES_grid_to_CLUBB_grid( fread_var, gr%zt, k,  &
                                     exact_lev_idx_zt(k), lower_lev_idx_zt(k),  &
                                     upper_lev_idx_zt(k), l_lin_int_zt(k) )
      enddo

      ! Find the lowest and highest indices of CLUBB momentum levels that fall
      ! within the domain of the LES output.
      call CLUBB_levels_within_LES_domain( fread_var, gr%zm,  &
                                           k_lowest_zm_input, k_highest_zm_input )

      allocate( exact_lev_idx_zm(k_lowest_zm_input:k_highest_zm_input) )
      allocate( lower_lev_idx_zm(k_lowest_zm_input:k_highest_zm_input) )
      allocate( upper_lev_idx_zm(k_lowest_zm_input:k_highest_zm_input) )
      allocate( l_lin_int_zm(k_lowest_zm_input:k_highest_zm_input) )

      ! For all CLUBB momentum levels, k, that are within the LES domain, find
      ! either the index of the LES level that exactly matches the altitude of the
      ! CLUBB level, or find the two indices of the LES levels that are on either
      ! side of the CLUBB level.
      do k = k_lowest_zm_input, k_highest_zm_input, 1
        ! CLUBB vertical level k is found at an altitude that is within the
        ! domain of the LES output.
        call LES_grid_to_CLUBB_grid( fread_var, gr%zm, k,  &
                                     exact_lev_idx_zm(k), lower_lev_idx_zm(k),  &
                                     upper_lev_idx_zm(k), l_lin_int_zm(k) )
      enddo


      ! Initialize l_fatal_error for case ( "rf1" )
      l_fatal_error = .false.

      if ( input_um ) then
        call get_grads_var( fread_var, "um", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of um from the LES GrADS file.
        do k = k_lowest_zt_input, k_highest_zt_input, 1
          if ( l_lin_int_zt(k) ) then
            ! CLUBB thermodynamic level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            um(k) = lin_int( gr%zt(k),  &
                             fread_var%z(upper_lev_idx_zt(k)),  &
                             fread_var%z(lower_lev_idx_zt(k)),  &
                             LES_tmp1(upper_lev_idx_zt(k)),  &
                             LES_tmp1(lower_lev_idx_zt(k)) )
          else
            ! CLUBB thermodynamic level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            um(k) = LES_tmp1(exact_lev_idx_zt(k))
          endif
        enddo
        ! When this is a standard scenario, where CLUBB thermodynamic level 3 is
        ! the first thermodynamic level at or above the lowest LES level, use the
        ! values of um at thermodynamic levels 4 and 3 to find the value at
        ! thermodynamic level 2 through the use of a linear extension.  Then, use
        ! the values of um at thermodynamic levels 3 and 2 to find the value at
        ! thermodynamic level 1 through the use of a linear extension.
        if ( k_lowest_zt_input == 3 ) then
          um(2)  & 
          = lin_ext_zt_bottom( um(4), um(3), & 
                               gr%zt(4), gr%zt(3), gr%zt(2) )
          um(1)  & 
          = lin_ext_zt_bottom( um(3), um(2), & 
                               gr%zt(3), gr%zt(2), gr%zt(1) )
        endif
      endif

      if ( input_vm ) then
        call get_grads_var( fread_var, "vm", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of vm from the LES GrADS file.
        do k = k_lowest_zt_input, k_highest_zt_input, 1
          if ( l_lin_int_zt(k) ) then
            ! CLUBB thermodynamic level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            vm(k) = lin_int( gr%zt(k),  &
                             fread_var%z(upper_lev_idx_zt(k)),  &
                             fread_var%z(lower_lev_idx_zt(k)),  &
                             LES_tmp1(upper_lev_idx_zt(k)),  &
                             LES_tmp1(lower_lev_idx_zt(k)) )
          else
            ! CLUBB thermodynamic level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            vm(k) = LES_tmp1(exact_lev_idx_zt(k))
          endif
        enddo
        ! When this is a standard scenario, where CLUBB thermodynamic level 3 is
        ! the first thermodynamic level at or above the lowest LES level, use the
        ! values of vm at thermodynamic levels 4 and 3 to find the value at
        ! thermodynamic level 2 through the use of a linear extension.  Then, use
        ! the values of vm at thermodynamic levels 3 and 2 to find the value at
        ! thermodynamic level 1 through the use of a linear extension.
        if ( k_lowest_zt_input == 3 ) then
          vm(2)  & 
          = lin_ext_zt_bottom( vm(4), vm(3), & 
                               gr%zt(4), gr%zt(3), gr%zt(2) )
          vm(1)  & 
          = lin_ext_zt_bottom( vm(3), vm(2), & 
                               gr%zt(3), gr%zt(2), gr%zt(1) )
        endif
      endif

      if ( input_rtm) then
        call get_grads_var( fread_var, "qtm", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of rtm from the LES GrADS file.
        do k = k_lowest_zt_input, k_highest_zt_input, 1
          if ( l_lin_int_zt(k) ) then
            ! CLUBB thermodynamic level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            rtm(k) = lin_int( gr%zt(k),  &
                              fread_var%z(upper_lev_idx_zt(k)),  &
                              fread_var%z(lower_lev_idx_zt(k)),  &
                              LES_tmp1(upper_lev_idx_zt(k)),  &
                              LES_tmp1(lower_lev_idx_zt(k)) )
          else
            ! CLUBB thermodynamic level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            rtm(k) = LES_tmp1(exact_lev_idx_zt(k))
          endif
        enddo
        ! When this is a standard scenario, where CLUBB thermodynamic level 3 is
        ! the first thermodynamic level at or above the lowest LES level, set the
        ! values of rtm at thermodynamic levels 2 and 1 to the value at
        ! thermodynamic level 3.
        if ( k_lowest_zt_input == 3 ) then
          rtm(1:2) = rtm(3)
        endif
      endif

      if ( input_thlm) then
        call get_grads_var( fread_var, "thlm", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of thlm from the LES GrADS file.
        do k = k_lowest_zt_input, k_highest_zt_input, 1
          if ( l_lin_int_zt(k) ) then
            ! CLUBB thermodynamic level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            thlm(k) = lin_int( gr%zt(k),  &
                               fread_var%z(upper_lev_idx_zt(k)),  &
                               fread_var%z(lower_lev_idx_zt(k)),  &
                               LES_tmp1(upper_lev_idx_zt(k)),  &
                               LES_tmp1(lower_lev_idx_zt(k)) )
          else
            ! CLUBB thermodynamic level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            thlm(k) = LES_tmp1(exact_lev_idx_zt(k))
          endif
        enddo
        ! When this is a standard scenario, where CLUBB thermodynamic level 3 is
        ! the first thermodynamic level at or above the lowest LES level, set the
        ! values of thlm at thermodynamic levels 2 and 1 to the value at
        ! thermodynamic level 3.
        if ( k_lowest_zt_input == 3 ) then
          thlm(1:2) = thlm(3)
        endif
      endif

      if ( input_wp3) then
        call get_grads_var( fread_var, "wp3", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of wp3 from the LES GrADS file.
        do k = k_lowest_zt_input, k_highest_zt_input, 1
          if ( l_lin_int_zt(k) ) then
            ! CLUBB thermodynamic level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            wp3(k) = lin_int( gr%zt(k),  &
                              fread_var%z(upper_lev_idx_zt(k)),  &
                              fread_var%z(lower_lev_idx_zt(k)),  &
                              LES_tmp1(upper_lev_idx_zt(k)),  &
                              LES_tmp1(lower_lev_idx_zt(k)) )
          else
            ! CLUBB thermodynamic level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            wp3(k) = LES_tmp1(exact_lev_idx_zt(k))
          endif
        enddo
        ! When this is a standard scenario, where CLUBB thermodynamic level 3 is
        ! the first thermodynamic level at or above the lowest LES level, set the
        ! values of wp3 at thermodynamic levels 2 and 1 to 0.
        if ( k_lowest_zt_input == 3 ) then
          wp3(1:2) = 0.
        endif
      endif

      if ( input_wprtp) then
        call get_grads_var( fread_var, "wpqtp", timestep, &
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of wprtp from the LES GrADS file.
        do k = k_lowest_zm_input, k_highest_zm_input, 1
          if ( l_lin_int_zm(k) ) then
            ! CLUBB momentum level k is found at an altitude that is between two
            ! LES levels.  Linear interpolation is required.
            wprtp(k) = lin_int( gr%zm(k),  &
                                fread_var%z(upper_lev_idx_zm(k)),  &
                                fread_var%z(lower_lev_idx_zm(k)),  &
                                LES_tmp1(upper_lev_idx_zm(k)),  &
                                LES_tmp1(lower_lev_idx_zm(k)) )
          else
            ! CLUBB momentum level k is found at an altitude that is an exact
            ! match with an LES level altitude.
            wprtp(k) = LES_tmp1(exact_lev_idx_zm(k))
          endif
        enddo
        ! When this is a standard scenario, where CLUBB momentum level 3 is the
        ! first momentum level above the lowest LES level, use the values of
        ! wprtp at momentum levels 4 and 3 to find the value at momentum level 2
        ! through the use of a linear extension.  Then, use the values of wprtp
        ! at momentum levels 3 and 2 to find the value at momentum level 1 through
        ! the use of a linear extension.  It should be pointed out that the
        ! boundary flux is usually solved in LES or CLUBB via a subroutine
        ! like sfc_var.
        if ( k_lowest_zm_input == 3 ) then
          wprtp(2)  & 
          = lin_ext_zm_bottom( wprtp(4), wprtp(3), & 
                               gr%zm(4), gr%zm(3), gr%zm(2) )
          wprtp(1)  & 
          = lin_ext_zm_bottom( wprtp(3), wprtp(2), & 
                               gr%zm(3), gr%zm(2), gr%zm(1) )
        endif
      endif

      if ( input_wpthlp) then
        call get_grads_var( fread_var, "wpthlp", timestep,  &
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of wpthlp from the LES GrADS file.
        do k = k_lowest_zm_input, k_highest_zm_input, 1
          if ( l_lin_int_zm(k) ) then
            ! CLUBB momentum level k is found at an altitude that is between two
            ! LES levels.  Linear interpolation is required.
            wpthlp(k) = lin_int( gr%zm(k),  &
                                 fread_var%z(upper_lev_idx_zm(k)),  &
                                 fread_var%z(lower_lev_idx_zm(k)),  &
                                 LES_tmp1(upper_lev_idx_zm(k)),  &
                                 LES_tmp1(lower_lev_idx_zm(k)) )
          else
            ! CLUBB momentum level k is found at an altitude that is an exact
            ! match with an LES level altitude.
            wpthlp(k) = LES_tmp1(exact_lev_idx_zm(k))
          endif
        enddo
        ! When this is a standard scenario, where CLUBB momentum level 3 is the
        ! first momentum level above the lowest LES level, use the values of
        ! wpthlp at momentum levels 4 and 3 to find the value at momentum level 2
        ! through the use of a linear extension.  Then, use the values of wpthlp
        ! at momentum levels 3 and 2 to find the value at momentum level 1 through
        ! the use of a linear extension.  It should be pointed out that the
        ! boundary flux is usually solved in LES or CLUBB via a subroutine
        ! like sfc_var.
        if ( k_lowest_zm_input == 3 ) then
          wpthlp(2)  & 
          = lin_ext_zm_bottom( wpthlp(4), wpthlp(3), & 
                               gr%zm(4), gr%zm(3), gr%zm(2) )
          wpthlp(1)  & 
          = lin_ext_zm_bottom( wpthlp(3), wpthlp(2), & 
                               gr%zm(3), gr%zm(2), gr%zm(1) )
        endif
      endif

      if ( input_rtp2) then
        call get_grads_var( fread_var, "qtp2", timestep, &
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of rtp2 from the LES GrADS file.
        do k = k_lowest_zm_input, k_highest_zm_input, 1
          if ( l_lin_int_zm(k) ) then
            ! CLUBB momentum level k is found at an altitude that is between two
            ! LES levels.  Linear interpolation is required.
            rtp2(k) = lin_int( gr%zm(k),  &
                               fread_var%z(upper_lev_idx_zm(k)),  &
                               fread_var%z(lower_lev_idx_zm(k)),  &
                               LES_tmp1(upper_lev_idx_zm(k)),  &
                               LES_tmp1(lower_lev_idx_zm(k)) )
          else
            ! CLUBB momentum level k is found at an altitude that is an exact
            ! match with an LES level altitude.
            rtp2(k) = LES_tmp1(exact_lev_idx_zm(k))
          endif
        enddo
        ! When this is a standard scenario, where CLUBB momentum level 3 is the
        ! first momentum level above the lowest LES level, set the values of rtp2
        ! at momentum levels 1 and 2 to the value at momentum level 3.
        ! Using a linear extension here resulted in negatives.
        if ( k_lowest_zm_input == 3 ) then
          rtp2(1:2) = rtp2(3)
        endif
        if ( any ( rtp2(1:gr%nnzp) < rttol**2 ) ) then
          do k=1, gr%nnzp
            rtp2(k) = max(rtp2(k), rttol**2)
          enddo
        endif
      endif

      if ( input_thlp2 ) then
        call get_grads_var( fread_var, "thlp2", timestep, &
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of thlp2 from the LES GrADS file.
        do k = k_lowest_zm_input, k_highest_zm_input, 1
          if ( l_lin_int_zm(k) ) then
            ! CLUBB momentum level k is found at an altitude that is between two
            ! LES levels.  Linear interpolation is required.
            thlp2(k) = lin_int( gr%zm(k),  &
                                fread_var%z(upper_lev_idx_zm(k)),  &
                                fread_var%z(lower_lev_idx_zm(k)),  &
                                LES_tmp1(upper_lev_idx_zm(k)),  &
                                LES_tmp1(lower_lev_idx_zm(k)) )
          else
            ! CLUBB momentum level k is found at an altitude that is an exact
            ! match with an LES level altitude.
            thlp2(k) = LES_tmp1(exact_lev_idx_zm(k))
          endif
        enddo
        ! When this is a standard scenario, where CLUBB momentum level 3 is the
        ! first momentum level above the lowest LES level, set the values of thlp2
        ! at momentum levels 1 and 2 to the value at momentum level 3.
        ! Using a linear extension here resulted in negatives.
        if ( k_lowest_zm_input == 3 ) then
          thlp2(1:2) = thlp2(3)
        endif
      endif

      if ( input_rtpthlp) then
        call get_grads_var( fread_var, "qtpthlp", timestep, &
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of rtpthlp from the LES GrADS file.
        do k = k_lowest_zm_input, k_highest_zm_input, 1
          if ( l_lin_int_zm(k) ) then
            ! CLUBB momentum level k is found at an altitude that is between two
            ! LES levels.  Linear interpolation is required.
            rtpthlp(k) = lin_int( gr%zm(k),  &
                                  fread_var%z(upper_lev_idx_zm(k)),  &
                                  fread_var%z(lower_lev_idx_zm(k)),  &
                                  LES_tmp1(upper_lev_idx_zm(k)),  &
                                  LES_tmp1(lower_lev_idx_zm(k)) )
          else
            ! CLUBB momentum level k is found at an altitude that is an exact
            ! match with an LES level altitude.
            rtpthlp(k) = LES_tmp1(exact_lev_idx_zm(k))
          endif
        enddo
        ! When this is a standard scenario, where CLUBB momentum level 3 is the
        ! first momentum level above the lowest LES level, use the values of
        ! rtpthlp at momentum levels 4 and 3 to find the value at momentum level 2
        ! through the use of a linear extension.  Then, use the values of rtpthlp
        ! at momentum levels 3 and 2 to find the value at momentum level 1 through
        ! the use of a linear extension.  It should be pointed out that the
        ! boundary flux is usually solved in LES or CLUBB via a subroutine
        ! like sfc_var.
        if ( k_lowest_zm_input == 3 ) then
          rtpthlp(2)  & 
          = lin_ext_zm_bottom( rtpthlp(4), rtpthlp(3), & 
                               gr%zm(4), gr%zm(3), gr%zm(2) )
          rtpthlp(1)  & 
          = lin_ext_zm_bottom( rtpthlp(3), rtpthlp(2), & 
                               gr%zm(3), gr%zm(2), gr%zm(1) )
        endif
      endif

      if ( l_fatal_error ) stop "oops, get_grads_var failed in grads_fields_reader"

      ! Deallocate temporary storage variable LES_tmp1.
      deallocate( LES_tmp1 )

      deallocate( exact_lev_idx_zt )
      deallocate( lower_lev_idx_zt )
      deallocate( upper_lev_idx_zt )
      deallocate( l_lin_int_zt )

      deallocate( exact_lev_idx_zm )
      deallocate( lower_lev_idx_zm )
      deallocate( upper_lev_idx_zm )
      deallocate( l_lin_int_zm )

      call close_grads_read( fread_var )


    case ( "les" )   ! COAMPS LES -- all other cases.

      ! stats_sm
      call open_grads_read( 15, trim(datafile)//"_coamps_sm.ctl",  & 
                            fread_var )
      l_fatal_error = .false.

      ! Temporarily store LES output in variable array LES_tmp1.
      ! Allocate LES_tmp1 based on lowest and highest vertical indices of LES
      ! output.
      allocate( LES_tmp1(fread_var%ia:fread_var%iz) )

      ! Find the lowest and highest indices of CLUBB thermodynamic levels that
      ! fall within the domain of the LES output.
      call CLUBB_levels_within_LES_domain( fread_var, gr%zt,  &
                                           k_lowest_zt_input, k_highest_zt_input )

      allocate( exact_lev_idx_zt(k_lowest_zt_input:k_highest_zt_input) )
      allocate( lower_lev_idx_zt(k_lowest_zt_input:k_highest_zt_input) )
      allocate( upper_lev_idx_zt(k_lowest_zt_input:k_highest_zt_input) )
      allocate( l_lin_int_zt(k_lowest_zt_input:k_highest_zt_input) )

      ! For all CLUBB thermodynamic levels, k, that are within the LES domain,
      ! find either the index of the LES level that exactly matches the altitude
      ! of the CLUBB level, or find the two indices of the LES levels that are on
      ! either side of the CLUBB level.
      do k = k_lowest_zt_input, k_highest_zt_input, 1
        ! CLUBB vertical level k is found at an altitude that is within the
        ! domain of the LES output.
        call LES_grid_to_CLUBB_grid( fread_var, gr%zt, k,  &
                                     exact_lev_idx_zt(k), lower_lev_idx_zt(k),  &
                                     upper_lev_idx_zt(k), l_lin_int_zt(k) )
      enddo

      ! Find the lowest and highest indices of CLUBB momentum levels that fall
      ! within the domain of the LES output.
      call CLUBB_levels_within_LES_domain( fread_var, gr%zm,  &
                                           k_lowest_zm_input, k_highest_zm_input )

      allocate( exact_lev_idx_zm(k_lowest_zm_input:k_highest_zm_input) )
      allocate( lower_lev_idx_zm(k_lowest_zm_input:k_highest_zm_input) )
      allocate( upper_lev_idx_zm(k_lowest_zm_input:k_highest_zm_input) )
      allocate( l_lin_int_zm(k_lowest_zm_input:k_highest_zm_input) )

      ! For all CLUBB momentum levels, k, that are within the LES domain, find
      ! either the index of the LES level that exactly matches the altitude of the
      ! CLUBB level, or find the two indices of the LES levels that are on either
      ! side of the CLUBB level.
      do k = k_lowest_zm_input, k_highest_zm_input, 1
        ! CLUBB vertical level k is found at an altitude that is within the
        ! domain of the LES output.
        call LES_grid_to_CLUBB_grid( fread_var, gr%zm, k,  &
                                     exact_lev_idx_zm(k), lower_lev_idx_zm(k),  &
                                     upper_lev_idx_zm(k), l_lin_int_zm(k) )
      enddo


      ! Initialize l_read_error for case ( "les" )
      l_read_error = .false.

      if ( input_um ) then
        call get_grads_var( fread_var, "um", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of um from the LES GrADS file.
        do k = k_lowest_zt_input, k_highest_zt_input, 1
          if ( l_lin_int_zt(k) ) then
            ! CLUBB thermodynamic level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            um(k) = lin_int( gr%zt(k),  &
                             fread_var%z(upper_lev_idx_zt(k)),  &
                             fread_var%z(lower_lev_idx_zt(k)),  &
                             LES_tmp1(upper_lev_idx_zt(k)),  &
                             LES_tmp1(lower_lev_idx_zt(k)) )
          else
            ! CLUBB thermodynamic level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            um(k) = LES_tmp1(exact_lev_idx_zt(k))
          endif
        enddo
        ! When this is a standard scenario, where CLUBB thermodynamic level 2 is
        ! the first thermodynamic level at or above the lowest LES level, use the
        ! values of um at thermodynamic levels 3 and 2 to find the value at
        ! thermodynamic level 1 through the use of a linear extension.
        if ( k_lowest_zt_input == 2 ) then
          um(1)  & 
          = lin_ext_zt_bottom( um(3), um(2), & 
                               gr%zt(3), gr%zt(2), gr%zt(1) )

        endif
      endif

      if ( input_vm ) then
        call get_grads_var( fread_var, "vm", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of vm from the LES GrADS file.
        do k = k_lowest_zt_input, k_highest_zt_input, 1
          if ( l_lin_int_zt(k) ) then
            ! CLUBB thermodynamic level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            vm(k) = lin_int( gr%zt(k),  &
                             fread_var%z(upper_lev_idx_zt(k)),  &
                             fread_var%z(lower_lev_idx_zt(k)),  &
                             LES_tmp1(upper_lev_idx_zt(k)),  &
                             LES_tmp1(lower_lev_idx_zt(k)) )
          else
            ! CLUBB thermodynamic level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            vm(k) = LES_tmp1(exact_lev_idx_zt(k))
          endif
        enddo
        ! When this is a standard scenario, where CLUBB thermodynamic level 2 is
        ! the first thermodynamic level at or above the lowest LES level, use the
        ! values of vm at thermodynamic levels 3 and 2 to find the value at
        ! thermodynamic level 1 through the use of a linear extension.
        if ( k_lowest_zt_input == 2 ) then
          vm(1)  & 
          = lin_ext_zt_bottom( vm(3), vm(2), & 
                               gr%zt(3), gr%zt(2), gr%zt(1) )
        endif
      endif

      if ( input_rtm ) then
        call get_grads_var( fread_var, "qtm", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of rtm from the LES GrADS file.
        do k = k_lowest_zt_input, k_highest_zt_input, 1
          if ( l_lin_int_zt(k) ) then
            ! CLUBB thermodynamic level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            rtm(k) = lin_int( gr%zt(k),  &
                              fread_var%z(upper_lev_idx_zt(k)),  &
                              fread_var%z(lower_lev_idx_zt(k)),  &
                              LES_tmp1(upper_lev_idx_zt(k)),  &
                              LES_tmp1(lower_lev_idx_zt(k)) )
          else
            ! CLUBB thermodynamic level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            rtm(k) = LES_tmp1(exact_lev_idx_zt(k))
          endif
        enddo
        ! When this is a standard scenario, where CLUBB thermodynamic level 2 is
        ! the first thermodynamic level at or above the lowest LES level, set the
        ! value of rtm at thermodynamic level 1 to the value at thermodynamic
        ! level 2, as it is done in advance_xm_wpxp.
        if ( k_lowest_zt_input == 2 ) then
          rtm(1) = rtm(2)
        endif
      endif

      if ( input_thlm ) then
        call get_grads_var( fread_var, "thlm", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of thlm from the LES GrADS file.
        do k = k_lowest_zt_input, k_highest_zt_input, 1
          if ( l_lin_int_zt(k) ) then
            ! CLUBB thermodynamic level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            thlm(k) = lin_int( gr%zt(k),  &
                               fread_var%z(upper_lev_idx_zt(k)),  &
                               fread_var%z(lower_lev_idx_zt(k)),  &
                               LES_tmp1(upper_lev_idx_zt(k)),  &
                               LES_tmp1(lower_lev_idx_zt(k)) )
          else
            ! CLUBB thermodynamic level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            thlm(k) = LES_tmp1(exact_lev_idx_zt(k))
          endif
        enddo
        ! When this is a standard scenario, where CLUBB thermodynamic level 2 is
        ! the first thermodynamic level at or above the lowest LES level, set the
        ! value of thlm at thermodynamic level 1 to the value at thermodynamic
        ! level 2, as it is done in advance_xm_wpxp.
        if ( k_lowest_zt_input == 2 ) then
          thlm(1) = thlm(2)
        endif
      endif

      ! We obtain wp2 from stats_sw

      if ( input_wp3 ) then
        call get_grads_var( fread_var, "wp3", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of wp3 from the LES GrADS file.
        do k = k_lowest_zt_input, k_highest_zt_input, 1
          if ( l_lin_int_zt(k) ) then
            ! CLUBB thermodynamic level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            wp3(k) = lin_int( gr%zt(k),  &
                              fread_var%z(upper_lev_idx_zt(k)),  &
                              fread_var%z(lower_lev_idx_zt(k)),  &
                              LES_tmp1(upper_lev_idx_zt(k)),  &
                              LES_tmp1(lower_lev_idx_zt(k)) )
          else
            ! CLUBB thermodynamic level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            wp3(k) = LES_tmp1(exact_lev_idx_zt(k))
          endif
        enddo
        ! When this is a standard scenario, where CLUBB thermodynamic level 2 is
        ! the first thermodynamic level at or above the lowest LES level, set the
        ! value of wp3 at thermodynamic level 1 to 0, as it is done in
        ! advance_wp2_wp3.
        if ( k_lowest_zt_input == 2 ) then
          wp3(1) = 0.0
        endif
      endif

      if ( input_wprtp ) then
        call get_grads_var( fread_var, "wpqtp", timestep, &
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of wprtp from the LES GrADS file.
        do k = k_lowest_zm_input, k_highest_zm_input, 1
          if ( l_lin_int_zm(k) ) then
            ! CLUBB momentum level k is found at an altitude that is between two
            ! LES levels.  Linear interpolation is required.
            wprtp(k) = lin_int( gr%zm(k),  &
                                fread_var%z(upper_lev_idx_zm(k)),  &
                                fread_var%z(lower_lev_idx_zm(k)),  &
                                LES_tmp1(upper_lev_idx_zm(k)),  &
                                LES_tmp1(lower_lev_idx_zm(k)) )
          else
            ! CLUBB momentum level k is found at an altitude that is an exact
            ! match with an LES level altitude.
            wprtp(k) = LES_tmp1(exact_lev_idx_zm(k))
          endif
        enddo
        ! When this is a standard scenario, where CLUBB momentum level 2 is the
        ! first momentum level above the lowest LES level, use the values of
        ! wprtp at momentum levels 3 and 2 to find the value at momentum level 1
        ! through the use of a linear extension.  It should be pointed out that
        ! the boundary flux is usually solved in LES or CLUBB via a subroutine
        ! like sfc_var.
        if ( k_lowest_zm_input == 2 ) then
          wprtp(1)  & 
          = lin_ext_zm_bottom( wprtp(3), wprtp(2), & 
                               gr%zm(3), gr%zm(2), gr%zm(1) )
        endif
      endif

      if ( input_wpthlp ) then
        call get_grads_var( fread_var, "wpthlp", timestep,  &
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of wpthlp from the LES GrADS file.
        do k = k_lowest_zm_input, k_highest_zm_input, 1
          if ( l_lin_int_zm(k) ) then
            ! CLUBB momentum level k is found at an altitude that is between two
            ! LES levels.  Linear interpolation is required.
            wpthlp(k) = lin_int( gr%zm(k),  &
                                 fread_var%z(upper_lev_idx_zm(k)),  &
                                 fread_var%z(lower_lev_idx_zm(k)),  &
                                 LES_tmp1(upper_lev_idx_zm(k)),  &
                                 LES_tmp1(lower_lev_idx_zm(k)) )
          else
            ! CLUBB momentum level k is found at an altitude that is an exact
            ! match with an LES level altitude.
            wpthlp(k) = LES_tmp1(exact_lev_idx_zm(k))
          endif
        enddo
        ! When this is a standard scenario, where CLUBB momentum level 2 is the
        ! first momentum level above the lowest LES level, use the values of
        ! wpthlp at momentum levels 3 and 2 to find the value at momentum level 1
        ! through the use of a linear extension.  It should be pointed out that
        ! the boundary flux is usually solved in LES or CLUBB via a subroutine
        ! like sfc_var.
        if ( k_lowest_zm_input == 2 ) then
          wpthlp(1)  & 
          = lin_ext_zm_bottom( wpthlp(3), wpthlp(2), & 
                               gr%zm(3), gr%zm(2), gr%zm(1) )
        endif
      endif

      if ( input_rtp2 ) then
        call get_grads_var( fread_var, "qtp2", timestep, &
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of rtp2 from the LES GrADS file.
        do k = k_lowest_zm_input, k_highest_zm_input, 1
          if ( l_lin_int_zm(k) ) then
            ! CLUBB momentum level k is found at an altitude that is between two
            ! LES levels.  Linear interpolation is required.
            rtp2(k) = lin_int( gr%zm(k),  &
                               fread_var%z(upper_lev_idx_zm(k)),  &
                               fread_var%z(lower_lev_idx_zm(k)),  &
                               LES_tmp1(upper_lev_idx_zm(k)),  &
                               LES_tmp1(lower_lev_idx_zm(k)) )
          else
            ! CLUBB momentum level k is found at an altitude that is an exact
            ! match with an LES level altitude.
            rtp2(k) = LES_tmp1(exact_lev_idx_zm(k))
          endif
        enddo
        ! When this is a standard scenario, where CLUBB momentum level 2 is the
        ! first momentum level above the lowest LES level, set the value of rtp2
        ! at momentum level 1 to the value at momentum level 2.
        ! Using a linear extension here resulted in negatives.
        if ( k_lowest_zm_input == 2 ) then
          rtp2(1) =  rtp2(2)
        endif
        if ( any( rtp2(1:gr%nnzp) < rttol**2 ) ) then
          do k=1, gr%nnzp
            rtp2(k) = max(rtp2(k), rttol**2)
          enddo
        endif
      endif

      if ( input_thlp2 ) then
        call get_grads_var( fread_var, "thlp2", timestep, &
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of thlp2 from the LES GrADS file.
        do k = k_lowest_zm_input, k_highest_zm_input, 1
          if ( l_lin_int_zm(k) ) then
            ! CLUBB momentum level k is found at an altitude that is between two
            ! LES levels.  Linear interpolation is required.
            thlp2(k) = lin_int( gr%zm(k),  &
                                fread_var%z(upper_lev_idx_zm(k)),  &
                                fread_var%z(lower_lev_idx_zm(k)),  &
                                LES_tmp1(upper_lev_idx_zm(k)),  &
                                LES_tmp1(lower_lev_idx_zm(k)) )
          else
            ! CLUBB momentum level k is found at an altitude that is an exact
            ! match with an LES level altitude.
            thlp2(k) = LES_tmp1(exact_lev_idx_zm(k))
          endif
        enddo
        ! When this is a standard scenario, where CLUBB momentum level 2 is the
        ! first momentum level above the lowest LES level, set the value of thlp2
        ! at momentum level 1 to the value at momentum level 2.
        ! Using a linear extension here resulted in negatives.
        if ( k_lowest_zm_input == 2 ) then
          thlp2(1) = thlp2(2)
        endif
        if ( any( thlp2(1:gr%nnzp) < thltol**2 ) ) then
          do k=1, gr%nnzp
            thlp2(k) = max(thlp2(k), thltol**2)
          enddo
        endif
      endif

      if ( input_rtpthlp) then
        call get_grads_var( fread_var, "qtpthlp", timestep, &
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of rtpthlp from the LES GrADS file.
        do k = k_lowest_zm_input, k_highest_zm_input, 1
          if ( l_lin_int_zm(k) ) then
            ! CLUBB momentum level k is found at an altitude that is between two
            ! LES levels.  Linear interpolation is required.
            rtpthlp(k) = lin_int( gr%zm(k),  &
                                  fread_var%z(upper_lev_idx_zm(k)),  &
                                  fread_var%z(lower_lev_idx_zm(k)),  &
                                  LES_tmp1(upper_lev_idx_zm(k)),  &
                                  LES_tmp1(lower_lev_idx_zm(k)) )
          else
            ! CLUBB momentum level k is found at an altitude that is an exact
            ! match with an LES level altitude.
            rtpthlp(k) = LES_tmp1(exact_lev_idx_zm(k))
          endif
        enddo
        ! When this is a standard scenario, where CLUBB momentum level 2 is the
        ! first momentum level above the lowest LES level, use the values of
        ! rtpthlp at momentum levels 3 and 2 to find the value at momentum level 1
        ! through the use of a linear extension.  It should be pointed out that
        ! the boundary flux is usually solved in LES or CLUBB via a subroutine
        ! like sfc_var.
        if ( k_lowest_zm_input == 2 ) then
          rtpthlp(1)  & 
          = lin_ext_zm_bottom( rtpthlp(3), rtpthlp(2), & 
                               gr%zm(3), gr%zm(2), gr%zm(1) )
        endif
      endif

      ! upwp/vpwp in stats_sw
      if ( input_ug ) then
        write(fstderr,*) "The variable ug is not setup for input_type = les"
        l_fatal_error = .true.
      end if

      if ( input_vg ) then
        write(fstderr,*) "The variable vg is not setup for input_type = les"
        l_fatal_error = .true.
      end if

      if ( input_rcm ) then
        call get_grads_var( fread_var, "qcm", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of rcm from the LES GrADS file.
        do k = k_lowest_zt_input, k_highest_zt_input, 1
          if ( l_lin_int_zt(k) ) then
            ! CLUBB thermodynamic level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            rcm(k) = lin_int( gr%zt(k),  &
                              fread_var%z(upper_lev_idx_zt(k)),  &
                              fread_var%z(lower_lev_idx_zt(k)),  &
                              LES_tmp1(upper_lev_idx_zt(k)),  &
                              LES_tmp1(lower_lev_idx_zt(k)) )
          else
            ! CLUBB thermodynamic level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            rcm(k) = LES_tmp1(exact_lev_idx_zt(k))
          endif
        enddo
      endif

      if ( input_wm_zt ) then
        call get_grads_var( fread_var, "wlsm", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of wm_zt from the LES GrADS file.
        do k = k_lowest_zt_input, k_highest_zt_input, 1
          if ( l_lin_int_zt(k) ) then
            ! CLUBB thermodynamic level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            wm_zt(k) = lin_int( gr%zt(k),  &
                              fread_var%z(upper_lev_idx_zt(k)),  &
                              fread_var%z(lower_lev_idx_zt(k)),  &
                              LES_tmp1(upper_lev_idx_zt(k)),  &
                              LES_tmp1(lower_lev_idx_zt(k)) )
          else
            ! CLUBB thermodynamic level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            wm_zt(k) = LES_tmp1(exact_lev_idx_zt(k))
          endif
        enddo
      endif

      if ( input_exner ) then
        call get_grads_var( fread_var, "ex0", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of exner from the LES GrADS file.
        do k = k_lowest_zt_input, k_highest_zt_input, 1
          if ( l_lin_int_zt(k) ) then
            ! CLUBB thermodynamic level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            exner(k) = lin_int( gr%zt(k),  &
                              fread_var%z(upper_lev_idx_zt(k)),  &
                              fread_var%z(lower_lev_idx_zt(k)),  &
                              LES_tmp1(upper_lev_idx_zt(k)),  &
                              LES_tmp1(lower_lev_idx_zt(k)) )
          else
            ! CLUBB thermodynamic level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            exner(k) = LES_tmp1(exact_lev_idx_zt(k))
          endif
        enddo
      endif

      if ( input_em ) then

        ! Read in SGS TKE
        call get_grads_var( fread_var, "em", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of em from the LES GrADS file.
        do k = k_lowest_zm_input, k_highest_zm_input, 1
          if ( l_lin_int_zm(k) ) then
            ! CLUBB momentum level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            em(k) = lin_int( gr%zm(k),  &
                              fread_var%z(upper_lev_idx_zm(k)),  &
                              fread_var%z(lower_lev_idx_zm(k)),  &
                              LES_tmp1(upper_lev_idx_zm(k)),  &
                              LES_tmp1(lower_lev_idx_zm(k)) )
          else
            ! CLUBB momentum level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            em(k) = LES_tmp1(exact_lev_idx_zm(k))
          endif
        enddo

        ! Read in Resolved TKE
        call get_grads_var( fread_var, "tke", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of em from the LES GrADS file.
        do k = k_lowest_zm_input, k_highest_zm_input, 1
          if ( l_lin_int_zm(k) ) then
            ! CLUBB thermodynamic level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            em(k) = em(k) + lin_int( gr%zm(k),  &
                              fread_var%z(upper_lev_idx_zm(k)),  &
                              fread_var%z(lower_lev_idx_zm(k)),  &
                              LES_tmp1(upper_lev_idx_zm(k)),  &
                              LES_tmp1(lower_lev_idx_zm(k)) )
          else
            ! CLUBB thermodynamic level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            em(k) = em(k) + LES_tmp1(exact_lev_idx_zm(k))
          endif
        enddo

        where ( em < emin ) em = emin

      endif ! input_em

      if ( input_p ) then
        call get_grads_var( fread_var, "pm", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of pressure from the LES GrADS file.
        do k = k_lowest_zt_input, k_highest_zt_input, 1
          if ( l_lin_int_zt(k) ) then
            ! CLUBB thermodynamic level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            p_in_Pa(k) = lin_int( gr%zt(k),  &
                              fread_var%z(upper_lev_idx_zt(k)),  &
                              fread_var%z(lower_lev_idx_zt(k)),  &
                              LES_tmp1(upper_lev_idx_zt(k)),  &
                              LES_tmp1(lower_lev_idx_zt(k)) )
          else
            ! CLUBB thermodynamic level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            p_in_Pa(k) = LES_tmp1(exact_lev_idx_zt(k))
          endif
        enddo
      endif

      if ( input_rho ) then
        call get_grads_var( fread_var, "dn0", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of rho from the LES GrADS file.
        do k = k_lowest_zt_input, k_highest_zt_input, 1
          if ( l_lin_int_zt(k) ) then
            ! CLUBB thermodynamic level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            rho(k) = lin_int( gr%zt(k),  &
                              fread_var%z(upper_lev_idx_zt(k)),  &
                              fread_var%z(lower_lev_idx_zt(k)),  &
                              LES_tmp1(upper_lev_idx_zt(k)),  &
                              LES_tmp1(lower_lev_idx_zt(k)) )
          else
            ! CLUBB thermodynamic level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            rho(k) = LES_tmp1(exact_lev_idx_zt(k))
          endif
        enddo
      endif

      ! rho_zm is in stats_sw

      if ( input_Lscale ) then
        write(fstderr,*) "The variable Lscale is not setup for input_type = les"
        l_fatal_error = .true.
      endif

      if ( input_Lscale_up ) then
        write(fstderr,*) "The variable Lscale_up is not setup for input_type = les"
        l_fatal_error = .true.
      end if

      if ( input_Lscale_down ) then
        write(fstderr,*) "The variable Lscale_down is not setup for input_type = les"
        l_fatal_error = .true.
      end if

      if ( input_Kh_zt ) then
        call get_grads_var( fread_var, "kh", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )
        ! LES_tmp1 is the value of mixing length from the LES GrADS file.
        do k = k_lowest_zt_input, k_highest_zt_input, 1
          if ( l_lin_int_zt(k) ) then
            ! CLUBB thermodynamic level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            Kh_zt(k) = lin_int( gr%zt(k),  &
                              fread_var%z(upper_lev_idx_zt(k)),  &
                              fread_var%z(lower_lev_idx_zt(k)),  &
                              LES_tmp1(upper_lev_idx_zt(k)),  &
                              LES_tmp1(lower_lev_idx_zt(k)) )
          else
            ! CLUBB thermodynamic level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            Kh_zt(k) = LES_tmp1(exact_lev_idx_zt(k))
          endif
        enddo
      endif

      if ( input_Kh_zm ) then
        write(fstderr,*) "The variable Kh_zm is not setup for input_type = les"
        l_fatal_error = .true.
      end if

      if ( input_tau_zt ) then
        write(fstderr,*) "The variable tau_zt is not setup for input_type = les"
        l_fatal_error = .true.
      end if

      if ( input_tau_zm ) then
        write(fstderr,*) "The variable tau_zm is not setup for input_type = les"
        l_fatal_error = .true.
      end if

      if ( input_wpthvp ) then
        call get_grads_var( fread_var, "wpthvp", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_fatal_error )
        ! LES_tmp1 is the value of wpthvp from the LES GrADS file.
        do k = k_lowest_zm_input, k_highest_zm_input, 1
          if ( l_lin_int_zm(k) ) then
            ! CLUBB momentum level k is found at an altitude that is between two
            ! LES levels.  Linear interpolation is required.
            wpthvp(k) = lin_int( gr%zm(k),  &
                              fread_var%z(upper_lev_idx_zm(k)),  &
                              fread_var%z(lower_lev_idx_zm(k)),  &
                              LES_tmp1(upper_lev_idx_zm(k)),  &
                              LES_tmp1(lower_lev_idx_zm(k)) )
          else
            ! CLUBB momentum level k is found at an altitude that is an exact
            ! match with an LES level altitude.
            wpthvp(k) = LES_tmp1(exact_lev_idx_zm(k))
          endif
        enddo
      endif

      if ( input_thl1 ) then
        write(fstderr,*) "The variable thl1 is not setup for input_type = les"
        l_fatal_error = .true.
      end if

      if ( input_thl2 ) then
        write(fstderr,*) "The variable thl2 is not setup for input_type = les"
        l_fatal_error = .true.
      end if

      if ( input_a ) then
        write(fstderr,*) "The variable a is not setup for input_type = les"
        l_fatal_error = .true.
      end if

      if ( input_s1 ) then
        write(fstderr,*) "The variable s1 is not setup for input_type = les"
        l_fatal_error = .true.
      end if

      if ( input_s2 ) then
        write(fstderr,*) "The variable s2 is not setup for input_type = les"
        l_fatal_error = .true.
      end if

      if ( input_ss1 ) then
        write(fstderr,*) "The variable ss1 is not setup for input_type = les"
        l_fatal_error = .true.
      end if

      if ( input_ss2 ) then
        write(fstderr,*) "The variable ss2 is not setup for input_type = les"
        l_fatal_error = .true.
      end if

      if ( input_rc1 ) then
        write(fstderr,*) "The variable rc1 is not setup for input_type = les"
        l_fatal_error = .true.
      end if

      if ( input_rc2 ) then
        write(fstderr,*) "The variable rc2 is not setup for input_type = les"
        l_fatal_error = .true.
      end if

      if ( input_thvm ) then
        call get_grads_var( fread_var, "thvm", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of thvm from the LES GrADS file.
        do k = k_lowest_zt_input, k_highest_zt_input, 1
          if ( l_lin_int_zt(k) ) then
            ! CLUBB thermodynamic level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            thvm(k) = lin_int( gr%zt(k),  &
                               fread_var%z(upper_lev_idx_zt(k)),  &
                               fread_var%z(lower_lev_idx_zt(k)),  &
                               LES_tmp1(upper_lev_idx_zt(k)),  &
                               LES_tmp1(lower_lev_idx_zt(k)) )
          else
            ! CLUBB thermodynamic level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            thvm(k) = LES_tmp1(exact_lev_idx_zt(k))
          endif
        enddo
      endif

      if ( input_rrainm ) then
        call get_grads_var( fread_var, "qrm", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of rain water mixing ratio from the LES GrADS file.
        do k = k_lowest_zt_input, k_highest_zt_input, 1

          if ( iirrainm < 1 ) then
            write(fstderr,*) "Rain water mixing ratio cannot be input with"// &
              " micro_scheme = "//micro_scheme
            l_fatal_error = .true.
            exit
          end if

          if ( l_lin_int_zt(k) ) then
            ! CLUBB thermodynamic level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            hydromet(k,iirrainm) = lin_int( gr%zt(k),  &
                               fread_var%z(upper_lev_idx_zt(k)),  &
                               fread_var%z(lower_lev_idx_zt(k)),  &
                               LES_tmp1(upper_lev_idx_zt(k)),  &
                               LES_tmp1(lower_lev_idx_zt(k)) )
          else
            ! CLUBB thermodynamic level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            hydromet(k,iirrainm) = LES_tmp1(exact_lev_idx_zt(k))
          endif
        enddo
      endif

      if ( input_Nrm ) then
        call get_grads_var( fread_var, "nrm", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of rain water mixing ratio from the LES GrADS file.
        do k = k_lowest_zt_input, k_highest_zt_input, 1

          if ( iiNrm < 1 ) then
            write(fstderr,*) "Rain droplet number conc. cannot be input with"// &
              " micro_scheme = "//micro_scheme
            l_fatal_error = .true.
            exit
          end if

          if ( l_lin_int_zt(k) ) then
            ! CLUBB thermodynamic level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            hydromet(k,iiNrm) = lin_int( gr%zt(k),  &
                               fread_var%z(upper_lev_idx_zt(k)),  &
                               fread_var%z(lower_lev_idx_zt(k)),  &
                               LES_tmp1(upper_lev_idx_zt(k)),  &
                               LES_tmp1(lower_lev_idx_zt(k)) )
          else
            ! CLUBB thermodynamic level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            hydromet(k,iiNrm) = LES_tmp1(exact_lev_idx_zt(k))
          endif
        enddo
      endif

      if ( input_Ncm ) then
        call get_grads_var( fread_var, "ncm", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of rain water mixing ratio from the LES GrADS file.
        do k = k_lowest_zt_input, k_highest_zt_input, 1

          if ( iiNcm < 1 ) then
            write(fstderr,*) "Cloud droplet number conc. cannot be input with"// &
              " micro_scheme = "//micro_scheme
            l_fatal_error = .true.
            exit
          end if

          if ( l_lin_int_zt(k) ) then
            ! CLUBB thermodynamic level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            hydromet(k,iiNcm) = lin_int( gr%zt(k),  &
                               fread_var%z(upper_lev_idx_zt(k)),  &
                               fread_var%z(lower_lev_idx_zt(k)),  &
                               LES_tmp1(upper_lev_idx_zt(k)),  &
                               LES_tmp1(lower_lev_idx_zt(k)) )
          else
            ! CLUBB thermodynamic level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            hydromet(k,iiNcm) = LES_tmp1(exact_lev_idx_zt(k))
          endif
        enddo
      endif


      if ( input_rsnowm ) then
        call get_grads_var( fread_var, "qsm", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of rain water mixing ratio from the LES GrADS file.
        do k = k_lowest_zt_input, k_highest_zt_input, 1
          if ( iirsnowm < 1 ) then
            write(fstderr,*) "Snow mixing ratio cannot be input with"// &
              " micro_scheme = "//micro_scheme
            l_fatal_error = .true.
            exit
          end if
          if ( l_lin_int_zt(k) ) then
            ! CLUBB thermodynamic level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            hydromet(k,iirsnowm) = lin_int( gr%zt(k),  &
                               fread_var%z(upper_lev_idx_zt(k)),  &
                               fread_var%z(lower_lev_idx_zt(k)),  &
                               LES_tmp1(upper_lev_idx_zt(k)),  &
                               LES_tmp1(lower_lev_idx_zt(k)) )
          else
            ! CLUBB thermodynamic level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            hydromet(k,iirsnowm) = LES_tmp1(exact_lev_idx_zt(k))
          endif
        enddo
      endif

      if ( input_ricem ) then
        call get_grads_var( fread_var, "qim", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of rain water mixing ratio from the LES GrADS file.
        do k = k_lowest_zt_input, k_highest_zt_input, 1
          if ( iiricem < 1 ) then
            write(fstderr,*) "Ice mixing ratio cannot be input with"// &
              " micro_scheme = "//micro_scheme
            l_fatal_error = .true.
            exit
          end if
          if ( l_lin_int_zt(k) ) then
            ! CLUBB thermodynamic level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            hydromet(k,iiricem) = lin_int( gr%zt(k),  &
                               fread_var%z(upper_lev_idx_zt(k)),  &
                               fread_var%z(lower_lev_idx_zt(k)),  &
                               LES_tmp1(upper_lev_idx_zt(k)),  &
                               LES_tmp1(lower_lev_idx_zt(k)) )
          else
            ! CLUBB thermodynamic level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            hydromet(k,iiricem) = LES_tmp1(exact_lev_idx_zt(k))
          endif
        enddo
      endif

      if ( input_rgraupelm ) then
        call get_grads_var( fread_var, "qgm", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of rain water mixing ratio from the LES GrADS file.
        do k = k_lowest_zt_input, k_highest_zt_input, 1
          if ( iirgraupelm < 1 ) then
            write(fstderr,*) "Graupel mixing ratio cannot be input with"// &
              " micro_scheme = "//micro_scheme
            l_fatal_error = .true.
            exit
          end if
          if ( l_lin_int_zt(k) ) then
            ! CLUBB thermodynamic level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            hydromet(k,iirgraupelm) = lin_int( gr%zt(k),  &
                               fread_var%z(upper_lev_idx_zt(k)),  &
                               fread_var%z(lower_lev_idx_zt(k)),  &
                               LES_tmp1(upper_lev_idx_zt(k)),  &
                               LES_tmp1(lower_lev_idx_zt(k)) )
          else
            ! CLUBB thermodynamic level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            hydromet(k,iirgraupelm) = LES_tmp1(exact_lev_idx_zt(k))
          endif
        enddo
      endif

      if ( input_Ncnm ) then
        call get_grads_var( fread_var, "ncnm", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of rain water mixing ratio from the LES GrADS file.
        do k = k_lowest_zt_input, k_highest_zt_input, 1
          if ( l_lin_int_zt(k) ) then
            ! CLUBB thermodynamic level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            Ncnm(k) = lin_int( gr%zt(k),  &
                               fread_var%z(upper_lev_idx_zt(k)),  &
                               fread_var%z(lower_lev_idx_zt(k)),  &
                               LES_tmp1(upper_lev_idx_zt(k)),  &
                               LES_tmp1(lower_lev_idx_zt(k)) )
          else
            ! CLUBB thermodynamic level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            Ncnm(k) = LES_tmp1(exact_lev_idx_zt(k))
          endif
        enddo
      endif

      if ( input_Nim ) then
        call get_grads_var( fread_var, "nim", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of rain water mixing ratio from the LES GrADS file.
        do k = k_lowest_zt_input, k_highest_zt_input, 1
          if ( iiNim < 1 ) then
            write(fstderr,*) "Ice number conc. cannot be input with"// &
              " micro_scheme = "//micro_scheme
            l_fatal_error = .true.
            exit
          end if
          if ( l_lin_int_zt(k) ) then
            ! CLUBB thermodynamic level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            hydromet(k,iiNim) = lin_int( gr%zt(k),  &
                               fread_var%z(upper_lev_idx_zt(k)),  &
                               fread_var%z(lower_lev_idx_zt(k)),  &
                               LES_tmp1(upper_lev_idx_zt(k)),  &
                               LES_tmp1(lower_lev_idx_zt(k)) )
          else
            ! CLUBB thermodynamic level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            hydromet(k,iiNim) = LES_tmp1(exact_lev_idx_zt(k))
          endif
        enddo
      endif

      if ( input_thlm_forcing ) then
        write(fstderr,*) "The variable thlm_forcing is not setup for input_type = les"
        l_fatal_error = .true.
      end if

      if ( input_rtm_forcing ) then
        write(fstderr,*) "The variable rtm_forcing is not setup for input_type = les"
        l_fatal_error = .true.
      end if

      if ( input_up2 ) then
        call get_grads_var( fread_var, "up2", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of up2 from the LES GrADS file.
        do k = k_lowest_zm_input, k_highest_zm_input, 1
          if ( l_lin_int_zm(k) ) then
            ! CLUBB momentum level k is found at an altitude that is between two
            ! LES levels.  Linear interpolation is required.
            up2(k) = lin_int( gr%zm(k),  &
                              fread_var%z(upper_lev_idx_zm(k)),  &
                              fread_var%z(lower_lev_idx_zm(k)),  &
                              LES_tmp1(upper_lev_idx_zm(k)),  &
                              LES_tmp1(lower_lev_idx_zm(k)) )
          else
            ! CLUBB momentum level k is found at an altitude that is an exact
            ! match with an LES level altitude.
            up2(k) = LES_tmp1(exact_lev_idx_zm(k))
          endif
        enddo

        ! Clip up2 to be no smaller than wtol_sqd
        where ( up2 < wtol_sqd ) up2 = wtol_sqd

      endif

      if ( input_vp2 ) then
        call get_grads_var( fread_var, "vp2", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of vp2 from the LES GrADS file.
        do k = k_lowest_zm_input, k_highest_zm_input, 1
          if ( l_lin_int_zm(k) ) then
            ! CLUBB momentum level k is found at an altitude that is between two
            ! LES levels.  Linear interpolation is required.
            vp2(k) = lin_int( gr%zm(k),  &
                              fread_var%z(upper_lev_idx_zm(k)),  &
                              fread_var%z(lower_lev_idx_zm(k)),  &
                              LES_tmp1(upper_lev_idx_zm(k)),  &
                              LES_tmp1(lower_lev_idx_zm(k)) )
          else
            ! CLUBB momentum level k is found at an altitude that is an exact
            ! match with an LES level altitude.
            vp2(k) = LES_tmp1(exact_lev_idx_zm(k))
          endif
        enddo

        where ( vp2 < wtol_sqd ) vp2 = wtol_sqd

      endif

      if ( input_sigma_sqd_w ) then
        write(fstderr,*) "The variable sigma_sqd_w is not setup for input_type = les"
        l_fatal_error = .true.
      end if

      if ( input_sigma_sqd_w_zt ) then
        write(fstderr,*) "The variable sigma_sqd_w_zt is not setup for input_type = les"
        l_fatal_error = .true.
      end if

      if ( input_cf ) then
        write(fstderr,*) "The variable cf is not setup for input_type = les"
        l_fatal_error = .true.
      end if

      if ( input_veg_T_in_K ) then
        write(fstderr,*) "The variable veg_T_in_K is not setup for input_type = les"
        l_fatal_error = .true.
      end if

      if ( input_deep_soil_T_in_K ) then
        write(fstderr,*) "The variable deep_soil_T_in_K is not setup for input_type = les"
        l_fatal_error = .true.
      end if

      if ( input_sfc_soil_T_in_K ) then
        write(fstderr,*) "The variable sfc_soil_T_in_K is not setup for input_type = les"
        l_fatal_error = .true.
      end if

      if ( l_fatal_error ) stop "oops, get_grads_var failed in grads_fields_reader"

      ! Deallocate temporary storage variable LES_tmp1.
      deallocate( LES_tmp1 )

      deallocate( exact_lev_idx_zt )
      deallocate( lower_lev_idx_zt )
      deallocate( upper_lev_idx_zt )
      deallocate( l_lin_int_zt )

      deallocate( exact_lev_idx_zm )
      deallocate( lower_lev_idx_zm )
      deallocate( upper_lev_idx_zm )
      deallocate( l_lin_int_zm )

      call close_grads_read( fread_var )

    end select


    select case ( input_type )

    case ( "les", "rf1" )

      ! stats_sw
      call open_grads_read( 15, trim(datafile)//"_coamps_sw.ctl",  & 
                            fread_var )

      ! Temporarily store LES output in variable array LES_tmp1.
      ! Allocate LES_tmp1 based on lowest and highest vertical indices of LES
      ! output.
      allocate( LES_tmp1(fread_var%ia:fread_var%iz) )

      ! Find the lowest and highest indices of CLUBB momentum levels that fall
      ! within the domain of the LES output.
      call CLUBB_levels_within_LES_domain( fread_var, gr%zm,  &
                                           k_lowest_zm_input, k_highest_zm_input )

      allocate( exact_lev_idx_zm(k_lowest_zm_input:k_highest_zm_input) )
      allocate( lower_lev_idx_zm(k_lowest_zm_input:k_highest_zm_input) )
      allocate( upper_lev_idx_zm(k_lowest_zm_input:k_highest_zm_input) )
      allocate( l_lin_int_zm(k_lowest_zm_input:k_highest_zm_input) )

      ! For all CLUBB momentum levels, k, that are within the LES domain, find
      ! either the index of the LES level that exactly matches the altitude of the
      ! CLUBB level, or find the two indices of the LES levels that are on either
      ! side of the CLUBB level.
      do k = k_lowest_zm_input, k_highest_zm_input, 1
        ! CLUBB vertical level k is found at an altitude that is within the
        ! domain of the LES output.
        call LES_grid_to_CLUBB_grid( fread_var, gr%zm, k,  &
                                     exact_lev_idx_zm(k), lower_lev_idx_zm(k),  &
                                     upper_lev_idx_zm(k), l_lin_int_zm(k) )
      enddo


      ! Note:  l_fatal_error has already been initialized for both case ( "les" ) and
      !        case ( "rf1" ).

      ! Note:  wpup_sgs and wpvp_sgs must be added to make the u'w' and v'w' terms
      !        as they are in CLUBB.

      if ( input_upwp) then

        call get_grads_var( fread_var, "wpup", timestep, &
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of upwp from the LES GrADS file.
        do k = k_lowest_zm_input, k_highest_zm_input, 1
          if ( l_lin_int_zm(k) ) then
            ! CLUBB momentum level k is found at an altitude that is between two
            ! LES levels.  Linear interpolation is required.
            upwp(k) = lin_int( gr%zm(k),  &
                               fread_var%z(upper_lev_idx_zm(k)),  &
                               fread_var%z(lower_lev_idx_zm(k)),  &
                               LES_tmp1(upper_lev_idx_zm(k)),  &
                               LES_tmp1(lower_lev_idx_zm(k)) )
          else
            ! CLUBB momentum level k is found at an altitude that is an exact
            ! match with an LES level altitude.
            upwp(k) = LES_tmp1(exact_lev_idx_zm(k))
          endif
        enddo

        call get_grads_var( fread_var, "wpup_sgs", timestep, &
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of upwp_sgs from the LES GrADS file.
        do k = k_lowest_zm_input, k_highest_zm_input, 1
          if ( l_lin_int_zm(k) ) then
            ! CLUBB momentum level k is found at an altitude that is between two
            ! LES levels.  Linear interpolation is required.
            upwp(k) = lin_int( gr%zm(k),  &
                               fread_var%z(upper_lev_idx_zm(k)),  &
                               fread_var%z(lower_lev_idx_zm(k)),  &
                               LES_tmp1(upper_lev_idx_zm(k)),  &
                               LES_tmp1(lower_lev_idx_zm(k)) )  &
                      + upwp(k)
          else
            ! CLUBB momentum level k is found at an altitude that is an exact
            ! match with an LES level altitude.
            upwp(k) = LES_tmp1(exact_lev_idx_zm(k)) + upwp(k)
          endif
        enddo

      endif

      if ( input_vpwp) then

        call get_grads_var( fread_var, "wpvp", timestep, &
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of vpwp from the LES GrADS file.
        do k = k_lowest_zm_input, k_highest_zm_input, 1
          if ( l_lin_int_zm(k) ) then
            ! CLUBB momentum level k is found at an altitude that is between two
            ! LES levels.  Linear interpolation is required.
            vpwp(k) = lin_int( gr%zm(k),  &
                               fread_var%z(upper_lev_idx_zm(k)),  &
                               fread_var%z(lower_lev_idx_zm(k)),  &
                               LES_tmp1(upper_lev_idx_zm(k)),  &
                               LES_tmp1(lower_lev_idx_zm(k)) )
          else
            ! CLUBB momentum level k is found at an altitude that is an exact
            ! match with an LES level altitude.
            vpwp(k) = LES_tmp1(exact_lev_idx_zm(k))
          endif
        enddo

        call get_grads_var( fread_var, "wpvp_sgs", timestep, &
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of vpwp_sgs from the LES GrADS file.
        do k = k_lowest_zm_input, k_highest_zm_input, 1
          if ( l_lin_int_zm(k) ) then
            ! CLUBB momentum level k is found at an altitude that is between two
            ! LES levels.  Linear interpolation is required.
            vpwp(k) = lin_int( gr%zm(k),  &
                               fread_var%z(upper_lev_idx_zm(k)),  &
                               fread_var%z(lower_lev_idx_zm(k)),  &
                               LES_tmp1(upper_lev_idx_zm(k)),  &
                               LES_tmp1(lower_lev_idx_zm(k)) )  &
                      + vpwp(k)
          else
            ! CLUBB momentum level k is found at an altitude that is an exact
            ! match with an LES level altitude.
            vpwp(k) = LES_tmp1(exact_lev_idx_zm(k)) + vpwp(k)
          endif
        enddo

      endif

      if ( input_wp2 ) then
        call get_grads_var( fread_var, "wp2", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of wp2 from the LES GrADS file.
        do k = k_lowest_zm_input, k_highest_zm_input, 1
          if ( l_lin_int_zm(k) ) then
            ! CLUBB momentum level k is found at an altitude that is between two
            ! LES levels.  Linear interpolation is required.
            wp2(k) = lin_int( gr%zm(k),  &
                              fread_var%z(upper_lev_idx_zm(k)),  &
                              fread_var%z(lower_lev_idx_zm(k)),  &
                              LES_tmp1(upper_lev_idx_zm(k)),  &
                              LES_tmp1(lower_lev_idx_zm(k)) )
          else
            ! CLUBB momentum level k is found at an altitude that is an exact
            ! match with an LES level altitude.
            wp2(k) = LES_tmp1(exact_lev_idx_zm(k))
          endif
        enddo
        if ( any( wp2(1:gr%nnzp) < wtol_sqd ) ) then
          do k=1, gr%nnzp
            wp2(k) = max( wp2(k), wtol_sqd )
          enddo
        endif
      endif

      if ( input_rho_zm ) then
        call get_grads_var( fread_var, "dn0", timestep, & 
                      LES_tmp1(fread_var%ia:fread_var%iz), l_read_error )

        l_fatal_error = l_fatal_error .or. l_read_error

        ! LES_tmp1 is the value of rho_zm from the LES GrADS file.
        do k = k_lowest_zm_input, k_highest_zm_input, 1
          if ( l_lin_int_zm(k) ) then
            ! CLUBB thermodynamic level k is found at an altitude that is
            ! between two LES levels.  Linear interpolation is required.
            rho_zm(k) = lin_int( gr%zt(k),  &
                              fread_var%z(upper_lev_idx_zm(k)),  &
                              fread_var%z(lower_lev_idx_zm(k)),  &
                              LES_tmp1(upper_lev_idx_zm(k)),  &
                              LES_tmp1(lower_lev_idx_zm(k)) )
          else
            ! CLUBB thermodynamic level k is found at an altitude that is an
            ! exact match with an LES level altitude.
            rho_zm(k) = LES_tmp1(exact_lev_idx_zm(k))
          endif
        enddo
      endif

      if ( l_fatal_error ) stop "get_grads_var failed for stats_sw in grads_fields_reader"

      ! Deallocate temporary storage variable LES_tmp1.
      deallocate( LES_tmp1 )

      deallocate( exact_lev_idx_zm )
      deallocate( lower_lev_idx_zm )
      deallocate( upper_lev_idx_zm )
      deallocate( l_lin_int_zm )

      call close_grads_read( fread_var )


    end select


    return
  end subroutine grads_fields_reader

!===============================================================================
  subroutine CLUBB_levels_within_LES_domain( fread_var, CLUBB_grid,  &
                                             k_lowest_input, k_highest_input )

    ! Description:
    ! This subroutine finds both the lowest and the highest CLUBB grid levels
    ! (for either the thermodynamic grid or the momentum grid) that are with the
    ! domain of the LES grid.

    ! References:
    !   None
    !-----------------------------------------------------------------------


    use stat_file_module, only:  &
        stat_file  ! Variable type

    use constants, only:  &
        fstderr ! Constant

    implicit none

    ! Input Variables.
    type(stat_file), intent(in) ::  &
      fread_var  ! Information about LES run.

    real, dimension(:), intent(in) ::  &
      CLUBB_grid ! Altitude of CLUBB grid levels
                 ! (either thermodynamic or momentum grid levels)  [m]

    ! Output Variables
    integer, intent(out) ::  &
      k_lowest_input,  & ! The lowest CLUBB level that's within the LES domain.
      k_highest_input    ! The highest CLUBB level that's within the LES domain.

    ! Local Variable
    integer :: k, kmax  ! Array index

    ! ---- Begin Code ----

    ! Find the lowest CLUBB level that falls within the LES domain.
    k    = size( CLUBB_grid )
    kmax = size( CLUBB_grid )
    do
      if ( CLUBB_grid(k) < fread_var%z(fread_var%ia) ) then

        if ( k == kmax ) then
          ! The bottom of the LES domain is above the top of the CLUBB
          ! domain.
          write(fstderr,*) "The lowest LES input level is above the top ",  &
                           "of the CLUBB model domain."
          stop "Error in CLUBB_levels_within_LES_domain"
        else
          ! Level k is the first CLUBB level below the LES domain.  Thus, the
          ! lowest CLUBB level within the LES domain has the index k + 1.
          k_lowest_input = k + 1
          exit
        endif

      elseif ( k == 1 ) then

        ! The bottom CLUBB level is within the LES domain.
        k_lowest_input = 1
        exit

      else   ! k > 1 and k <= kmax; level not yet found.

        ! Increment one more CLUBB vertical level down.
        k = k - 1

      endif
    enddo

    ! Find the highest CLUBB level that falls within the LES domain.
    k = 1
    do
      if ( CLUBB_grid(k) > fread_var%z(fread_var%iz) ) then

        if ( k == 1 ) then
          ! The top of the LES domain is below the bottom of the CLUBB
          ! domain.
          write(fstderr,*) "The highest LES input level is below the ",  &
                           "bottom of the CLUBB model domain."
          stop "Error in CLUBB_levels_within_LES_domain"
        else
          ! Level k is the first CLUBB level above the LES domain.  Thus, the
          ! highest CLUBB level within the LES domain has the index k - 1.
          k_highest_input = k - 1
          exit
        endif

      elseif ( k == kmax ) then

        ! The top CLUBB level is within the LES domain.
        k_highest_input = kmax 
        exit

      else   ! k < kmax and k >= 1; level not yet found.

        ! Increment one more CLUBB vertical level up.
        k = k + 1

      endif
    enddo

    return
  end subroutine CLUBB_levels_within_LES_domain

!===============================================================================
  pure subroutine LES_grid_to_CLUBB_grid( fread_var, CLUBB_grid, k,  &
                                     exact_lev_idx, lower_lev_idx,  &
                                     upper_lev_idx, l_lin_int )

    ! Description:
    ! Finds the level on the LES grid that is exactly even with the CLUBB
    ! grid level (either thermodynamic or momentum grid) that is input
    ! (level k).  Else, it finds the two LES levels that sandwich the CLUBB
    ! grid level that is input.

    !-----------------------------------------------------------------------

    use stat_file_module, only:  &
        stat_file  ! Variable type

    implicit none

    ! Input Variables.
    type(stat_file), intent(in) ::  &
      fread_var  ! Information about LES run.

    real, dimension(:), intent(in) ::  &
      CLUBB_grid ! Altitude of CLUBB grid levels
                 ! (either thermodynamic or momentum grid levels)  [m]

    integer, intent(in) ::  &
      k  ! Index of CLUBB vertical level that is being compared to.

    ! Output Variables.
    integer, intent(out) ::  &
      exact_lev_idx, & ! In case of an exact match, index of LES level that is
                       ! exactly even with CLUBB level k.
      lower_lev_idx, & ! In case linear interpolation is needed, index of LES
                       ! level that is immediately below CLUBB level k.
      upper_lev_idx    ! In case linear interpolation is needed, index of LES
                       ! level that is immediately above CLUBB level k.

    logical, intent(out) ::  &
      l_lin_int  ! Flag that is turned on if linear interpolation is needed.

    ! Local Variable.
    integer :: j

    ! Initialize the output quantities.
    exact_lev_idx = 0
    lower_lev_idx = 0
    upper_lev_idx = 0
    l_lin_int     = .false.

    ! Initialize LES vertical grid loop index, j, at the lowest LES grid index,
    ! which is fread_var%ia.
    j = fread_var%ia

    do

      if ( fread_var%z(j) == CLUBB_grid(k) ) then

        ! There is an LES level altitude at LES level j that is an exact
        ! match to the CLUBB level altitude at CLUBB grid level k.
        exact_lev_idx = j
        l_lin_int = .false.

      elseif ( fread_var%z(j) < CLUBB_grid(k) ) then

        ! The LES level altitude at LES level j is lower than the CLUBB level
        ! altitude at CLUBB grid level k.
        lower_lev_idx = j

      else   ! fread_var%z(j) > CLUBB_grid(k)

        ! The LES level altitude at LES level j is higher than the CLUBB level
        ! altitude at CLUBB grid level k.
        upper_lev_idx = j

      endif

      if ( exact_lev_idx > 0 ) exit  ! An exact answer has been found,
      ! exit the loop.

      if ( upper_lev_idx == lower_lev_idx + 1 ) then

        ! CLUBB level k has been found between two successive LES levels.
        ! Linear interpolation is needed.  An answer has been found, exit
        ! the loop.
        l_lin_int = .true.
        exit

      endif

      ! An answer has not been found yet, iterate the j index.
      j = j + 1

    enddo

    return
  end subroutine LES_grid_to_CLUBB_grid

!===============================================================================
  pure function lin_ext_zm_bottom( var_zmp2, var_zmp1,  & 
                                   zmp2, zmp1, zm )  & 
  result( var_zm )

    ! Description:
    ! This function computes the value of a momentum-level variable at a bottom
    ! grid level by using a linear extension of the values of the variable at
    ! the two levels immediately above the level where the result value is
    ! needed.

    !-----------------------------------------------------------------------

    implicit none


    ! Input Variables
    real, intent(in) :: & 
      var_zmp2,    & ! Momentum level variable at level (k+2)  [units vary]
      var_zmp1,    & ! Momentum level variable at level (k+1)  [units vary]
      zmp2,        & ! Altitude at momentum level (k+2)        [m]
      zmp1,        & ! Altitude at momentum level (k+1)        [m]
      zm             ! Altitude at momentum level (k)          [m]

    ! Return Variable
    real :: var_zm   ! Momentum level variable at level (k)    [units vary]

    var_zm = ( ( var_zmp2 - var_zmp1 ) / ( zmp2 - zmp1 ) ) & 
             * ( zm - zmp1 ) + var_zmp1

    return
  end function lin_ext_zm_bottom

!===============================================================================
  pure function lin_ext_zt_bottom( var_ztp2, var_ztp1,  & 
                                   ztp2, ztp1, zt )  & 
  result( var_zt )

    ! Description:
    ! This function computes the value of a thermodynamic-level variable at a
    ! bottom grid level by using a linear extension of the values of the
    ! variable at the two levels immediately above the level where the result
    ! value is needed.

    !-----------------------------------------------------------------------

    implicit none

    ! Input Variables
    real, intent(in) :: & 
      var_ztp2,    & ! Thermodynamic level variable at level (k+2)  [units vary]
      var_ztp1,    & ! Thermodynamic level variable at level (k+1)  [units vary]
      ztp2,        & ! Altitude at thermodynamic level (k+2)        [m]
      ztp1,        & ! Altitude at thermodynamic level (k+1)        [m]
      zt             ! Altitude at thermodynamic level (k)          [m]

    ! Return Variable
    real :: var_zt   ! Thermodynamic level variable at level (k)    [units vary]

    var_zt = ( ( var_ztp2 - var_ztp1 ) / ( ztp2 - ztp1 ) ) & 
             * ( zt - ztp1 ) + var_ztp1

    return
  end function lin_ext_zt_bottom

!===============================================================================
  subroutine compute_timestep( iunit, filename, l_restart, & 
                               time, nearest_timestep )

    ! Description:
    ! Given a time 'time', determines the closest output time in a GrADS file.

    !-----------------------------------------------------------------------

    use stat_file_module, only: & 
        stat_file     ! Type

    use inputfile_class, only: &
        open_grads_read,  & ! Procedure(s)
        close_grads_read

    use constants, only:  & 
        sec_per_min ! Variable(s)

    use stats_precision, only:  & 
        time_precision

    implicit none

    ! Input Variable(s)
    integer, intent(in) :: iunit ! File I/O unit

    character(len=*), intent(in) ::filename ! Path to the file and its name

    logical, intent(in) :: l_restart ! Whether this is a restart run

    real(kind=time_precision), intent(in) ::  & 
      time ! Time near which we want to find GrADS output,
    ! e.g. time_restart     [s]


    ! Output Variable(s)
    integer, intent(out) ::  & 
      nearest_timestep ! Nearest GrADS output time to time [min]

    ! Local Variables
    type (stat_file) :: fread_var

    real(kind=time_precision) :: delta_time   ! In seconds

    call open_grads_read( iunit, trim( filename ), fread_var )

    ! (restart time) - (initial time)
    delta_time = time - (fread_var%time - fread_var%dtwrite)

    !    Joshua Fasching March 2008
!     .        time - fread_var%time

    ! Reporting
    if ( l_restart ) then
      print *, "Initial time of GrADS reference file ", & 
               "[seconds since midnight]: ",  & 
               fread_var%time
      print *, "Model restart time [s]: ", time
      print *, "Elapsed time between ", & 
               "initial time of ref file and restart time [s]: ",  & 
               delta_time
      print *, "GrADS file output time interval [s]: ",  & 
               fread_var%dtwrite

      if ( ( mod( delta_time , fread_var%dtwrite )  > 1e-8 ) .or.  & 
           ( mod( delta_time, fread_var%dtwrite ) < -1e-8 ) ) then
        print*, "Error: Elapsed time is not a multiple ", & 
                "of the reference GrADS output time interval."
        print*, "Elapsed time [s] = ", delta_time
        print*, "GrADS output time interval = ", fread_var%dtwrite
        stop
      end if

      if ( mod( delta_time , sec_per_min ) > 1e-8 & 
            .or. mod( delta_time, sec_per_min ) < -1e-8 ) then
        print*, "Error: Elapsed time is not a multiple ", & 
                "of one minute."
        print*, "Elapsed time [s] = ", delta_time
        stop
      end if

    end if ! l_restart

    ! Determines the closest recorded timestep to the restart
    ! time.
    nearest_timestep = nint( delta_time / sec_per_min )

    if ( l_restart ) then
      print *, "Elapsed time between ", & 
               "initial time of ref file and restart time ", & 
               "rounded to nearest minute: ",  & 
               nearest_timestep

      ! Print the actual record being recalled.
      ! Joshua Fasching March 2008
      print *, "Nearest GrADS output time iteration [ ]: ", & 
               nint( nearest_timestep /  & 
                     (fread_var%dtwrite/sec_per_min) ) - 1
    end if ! l_restart

    call close_grads_read( fread_var )

  end subroutine compute_timestep

!===============================================================================

end module inputfields
