!---------------------------------------------------------------------------
! $Id$
!===============================================================================
module stats_zt_module

  implicit none

  private ! Default Scope

  public :: stats_init_zt

  ! Constant parameters
  integer, parameter, public :: nvarmax_zt = 800 ! Maximum variables allowed

  contains

  !=============================================================================
  subroutine stats_init_zt( vars_zt,                    & ! intent(in)
                            l_error,                    & ! intent(inout)
                            stats_metadata, stats_zt )    ! intent(inout)

    ! Description:
    ! Initializes array indices for stats_zt

    ! Note:
    ! All code that is within subroutine stats_init_zt, including variable
    ! allocation code, is not called if l_stats is false.  This subroutine is
    ! called only when l_stats is true.

    !-----------------------------------------------------------------------

    use constants_clubb, only:  &
        fstderr ! Constant(s)


    use stats_type_utilities, only: &
        stat_assign ! Procedure

    use parameters_model, only: &
        hydromet_dim, & ! Variable(s)
        sclr_dim,     &
        edsclr_dim

    use array_index, only: &
        hydromet_list, &  ! Variable(s)
        l_mix_rat_hm

    use stats_type, only: &
        stats ! Type

    use stats_variables, only: &
        stats_metadata_type

    implicit none

    ! External
    intrinsic :: trim

    ! Local Constants

    !--------------------- Input Variable ---------------------
    character(len= * ), dimension(nvarmax_zt), intent(in) :: &
      vars_zt

    !--------------------- InOut Variables ---------------------      
    type (stats_metadata_type), intent(inout) :: &
      stats_metadata

    type (stats), target, intent(inout) :: &
      stats_zt

    logical, intent(inout) :: l_error
 
    !--------------------- Local Varables ---------------------
    integer :: tot_zt_loops

    integer :: i, j, k

    integer :: hm_idx, hmx_idx, hmy_idx

    character(len=10) :: hm_type, hmx_type, hmy_type

    character(len=50) :: sclr_idx

    !--------------------- Begin Code ---------------------

    ! The default initialization for array indices for stats_zt is zero (see module
    ! stats_variables)

    ! If any of the index arrays are allocated, then we have called this before
    ! to set up stats_metadata, so all we want to do is set stats_zt via stats_assign
    if ( .not. allocated(stats_metadata%ihm_1) ) then

      ! Allocate and initialize hydrometeor statistical variables.
      allocate( stats_metadata%ihm_1(1:hydromet_dim) )
      allocate( stats_metadata%ihm_2(1:hydromet_dim) )
      allocate( stats_metadata%imu_hm_1(1:hydromet_dim) )
      allocate( stats_metadata%imu_hm_2(1:hydromet_dim) )
      allocate( stats_metadata%imu_hm_1_n(1:hydromet_dim) )
      allocate( stats_metadata%imu_hm_2_n(1:hydromet_dim) )
      allocate( stats_metadata%isigma_hm_1(1:hydromet_dim) )
      allocate( stats_metadata%isigma_hm_2(1:hydromet_dim) )
      allocate( stats_metadata%isigma_hm_1_n(1:hydromet_dim) )
      allocate( stats_metadata%isigma_hm_2_n(1:hydromet_dim) )

      allocate( stats_metadata%icorr_w_hm_1(1:hydromet_dim) )
      allocate( stats_metadata%icorr_w_hm_2(1:hydromet_dim) )
      allocate( stats_metadata%icorr_chi_hm_1(1:hydromet_dim) )
      allocate( stats_metadata%icorr_chi_hm_2(1:hydromet_dim) )
      allocate( stats_metadata%icorr_eta_hm_1(1:hydromet_dim) )
      allocate( stats_metadata%icorr_eta_hm_2(1:hydromet_dim) )
      allocate( stats_metadata%icorr_Ncn_hm_1(1:hydromet_dim) )
      allocate( stats_metadata%icorr_Ncn_hm_2(1:hydromet_dim) )
      allocate( stats_metadata%icorr_hmx_hmy_1(1:hydromet_dim,1:hydromet_dim) )
      allocate( stats_metadata%icorr_hmx_hmy_2(1:hydromet_dim,1:hydromet_dim) )

      allocate( stats_metadata%icorr_w_hm_1_n(1:hydromet_dim) )
      allocate( stats_metadata%icorr_w_hm_2_n(1:hydromet_dim) )
      allocate( stats_metadata%icorr_chi_hm_1_n(1:hydromet_dim) )
      allocate( stats_metadata%icorr_chi_hm_2_n(1:hydromet_dim) )
      allocate( stats_metadata%icorr_eta_hm_1_n(1:hydromet_dim) )
      allocate( stats_metadata%icorr_eta_hm_2_n(1:hydromet_dim) )
      allocate( stats_metadata%icorr_Ncn_hm_1_n(1:hydromet_dim) )
      allocate( stats_metadata%icorr_Ncn_hm_2_n(1:hydromet_dim) )
      allocate( stats_metadata%icorr_hmx_hmy_1_n(1:hydromet_dim,1:hydromet_dim) )
      allocate( stats_metadata%icorr_hmx_hmy_2_n(1:hydromet_dim,1:hydromet_dim) )

      allocate( stats_metadata%ihmp2_zt(1:hydromet_dim) )

      allocate( stats_metadata%iwp2hmp(1:hydromet_dim) )

      stats_metadata%ihm_1(:) = 0
      stats_metadata%ihm_2(:) = 0
      stats_metadata%imu_hm_1(:) = 0
      stats_metadata%imu_hm_2(:) = 0
      stats_metadata%imu_hm_1_n(:) = 0
      stats_metadata%imu_hm_2_n(:) = 0
      stats_metadata%isigma_hm_1(:) = 0
      stats_metadata%isigma_hm_2(:) = 0
      stats_metadata%isigma_hm_1_n(:) = 0
      stats_metadata%isigma_hm_2_n(:) = 0

      stats_metadata%icorr_w_hm_1(:) = 0
      stats_metadata%icorr_w_hm_2(:) = 0
      stats_metadata%icorr_chi_hm_1(:) = 0
      stats_metadata%icorr_chi_hm_2(:) = 0
      stats_metadata%icorr_eta_hm_1(:) = 0
      stats_metadata%icorr_eta_hm_2(:) = 0
      stats_metadata%icorr_Ncn_hm_1(:) = 0
      stats_metadata%icorr_Ncn_hm_2(:) = 0
      stats_metadata%icorr_hmx_hmy_1(:,:) = 0
      stats_metadata%icorr_hmx_hmy_2(:,:) = 0

      stats_metadata%icorr_w_hm_1_n(:) = 0
      stats_metadata%icorr_w_hm_2_n(:) = 0
      stats_metadata%icorr_chi_hm_1_n(:) = 0
      stats_metadata%icorr_chi_hm_2_n(:) = 0
      stats_metadata%icorr_eta_hm_1_n(:) = 0
      stats_metadata%icorr_eta_hm_2_n(:) = 0
      stats_metadata%icorr_Ncn_hm_1_n(:) = 0
      stats_metadata%icorr_Ncn_hm_2_n(:) = 0
      stats_metadata%icorr_hmx_hmy_1_n(:,:) = 0
      stats_metadata%icorr_hmx_hmy_2_n(:,:) = 0

      stats_metadata%ihmp2_zt(:) = 0

      stats_metadata%iwp2hmp(:) = 0

      ! Allocate and then zero out passive scalar arrays
      allocate( stats_metadata%isclrm(1:sclr_dim) )
      allocate( stats_metadata%isclrm_f(1:sclr_dim) )

      stats_metadata%isclrm(:)     = 0
      stats_metadata%isclrm_f(:)   = 0

      allocate( stats_metadata%iedsclrm(1:edsclr_dim) )
      allocate( stats_metadata%iedsclrm_f(1:edsclr_dim) )

      stats_metadata%iedsclrm(:)   = 0
      stats_metadata%iedsclrm_f(:) = 0

    end if

    ! Assign pointers for statistics variables stats_zt using stat_assign

    tot_zt_loops = stats_zt%num_output_fields

    if ( any( vars_zt == "hm_i" ) ) then
       ! Correct for number of variables found under "hm_i".
       ! Subtract 2 from the loop size (1st PDF component and 2nd PDF component)
       ! for each hydrometeor.
       tot_zt_loops = tot_zt_loops - 2 * hydromet_dim
       ! Add 1 for "hm_i" to the loop size.
       tot_zt_loops = tot_zt_loops + 1
    endif
    if ( any( vars_zt == "mu_hm_i" ) ) then
       ! Correct for number of variables found under "mu_hm_i".
       ! Subtract 2 from the loop size (1st PDF component and 2nd PDF component)
       ! for each hydrometeor.
       tot_zt_loops = tot_zt_loops - 2 * hydromet_dim
       ! Add 1 for "mu_hm_i" to the loop size.
       tot_zt_loops = tot_zt_loops + 1
    endif
    if ( any( vars_zt == "mu_Ncn_i" ) ) then
       ! Correct for number of variables found under "mu_Ncn_i".
       ! Subtract 2 from the loop size (1st PDF comp. and 2nd PDF comp.).
       tot_zt_loops = tot_zt_loops - 2
       ! Add 1 for "mu_Ncn_i" to the loop size.
       tot_zt_loops = tot_zt_loops + 1
    endif
    if ( any( vars_zt == "mu_hm_i_n" ) ) then
       ! Correct for number of variables found under "mu_hm_i_n".
       ! Subtract 2 from the loop size (1st PDF component and 2nd PDF component)
       ! for each hydrometeor.
       tot_zt_loops = tot_zt_loops - 2 * hydromet_dim
       ! Add 1 for "mu_hm_i_n" to the loop size.
       tot_zt_loops = tot_zt_loops + 1
    endif
    if ( any( vars_zt == "mu_Ncn_i_n" ) ) then
       ! Correct for number of variables found under "mu_Ncn_i_n".
       ! Subtract 2 from the loop size (1st PDF comp. and 2nd PDF comp.).
       tot_zt_loops = tot_zt_loops - 2
       ! Add 1 for "mu_Ncn_i_n" to the loop size.
       tot_zt_loops = tot_zt_loops + 1
    endif
    if ( any( vars_zt == "sigma_hm_i" ) ) then
       ! Correct for number of variables found under "sigma_hm_i".
       ! Subtract 2 from the loop size (1st PDF component and 2nd PDF component)
       ! for each hydrometeor.
       tot_zt_loops = tot_zt_loops - 2 * hydromet_dim
       ! Add 1 for "sigma_hm_i" to the loop size.
       tot_zt_loops = tot_zt_loops + 1
    endif
    if ( any( vars_zt == "sigma_Ncn_i" ) ) then
       ! Correct for number of variables found under "sigma_Ncn_i".
       ! Subtract 2 from the loop size (1st PDF comp. and 2nd PDF comp.).
       tot_zt_loops = tot_zt_loops - 2
       ! Add 1 for "sigma_Ncn_i" to the loop size.
       tot_zt_loops = tot_zt_loops + 1
    endif
    if ( any( vars_zt == "sigma_hm_i_n" ) ) then
       ! Correct for number of variables found under "sigma_hm_i_n".
       ! Subtract 2 from the loop size (1st PDF component and 2nd PDF component)
       ! for each hydrometeor.
       tot_zt_loops = tot_zt_loops - 2 * hydromet_dim
       ! Add 1 for "sigma_hm_i_n" to the loop size.
       tot_zt_loops = tot_zt_loops + 1
    endif
    if ( any( vars_zt == "sigma_Ncn_i_n" ) ) then
       ! Correct for number of variables found under "sigma_Ncn_i_n".
       ! Subtract 2 from the loop size (1st PDF comp. and 2nd PDF comp.).
       tot_zt_loops = tot_zt_loops - 2
       ! Add 1 for "sigma_Ncn_i_n" to the loop size.
       tot_zt_loops = tot_zt_loops + 1
    endif

    if ( any( vars_zt == "corr_w_hm_i" ) ) then
       ! Correct for number of variables found under "corr_whm_i".
       ! Subtract 2 from the loop size (1st PDF component and 2nd PDF component)
       ! for each hydrometeor.
       tot_zt_loops = tot_zt_loops - 2 * hydromet_dim
       ! Add 1 for "corr_whm_i" to the loop size.
       tot_zt_loops = tot_zt_loops + 1
    endif
    if ( any( vars_zt == "corr_w_Ncn_i" ) ) then
       ! Correct for number of variables found under "corr_wNcn_i".
       ! Subtract 2 from the loop size (1st PDF comp. and 2nd PDF comp.).
       tot_zt_loops = tot_zt_loops - 2
       ! Add 1 for "corr_wNcn_i" to the loop size.
       tot_zt_loops = tot_zt_loops + 1
    endif
    if ( any( vars_zt == "corr_chi_hm_i" ) ) then
       ! Correct for number of variables found under "corr_chi_hm_i".
       ! Subtract 2 from the loop size (1st PDF component and 2nd PDF component)
       ! for each hydrometeor.
       tot_zt_loops = tot_zt_loops - 2 * hydromet_dim
       ! Add 1 for "corr_chi_hm_i" to the loop size.
       tot_zt_loops = tot_zt_loops + 1
    endif
    if ( any( vars_zt == "corr_chi_Ncn_i" ) ) then
       ! Correct for number of variables found under "corr_chi_Ncn_i".
       ! Subtract 2 from the loop size (1st PDF comp. and 2nd PDF comp.).
       tot_zt_loops = tot_zt_loops - 2
       ! Add 1 for "corr_chi_Ncn_i" to the loop size.
       tot_zt_loops = tot_zt_loops + 1
    endif
    if ( any( vars_zt == "corr_eta_hm_i" ) ) then
       ! Correct for number of variables found under "corr_eta_hm_i".
       ! Subtract 2 from the loop size (1st PDF component and 2nd PDF component)
       ! for each hydrometeor.
       tot_zt_loops = tot_zt_loops - 2 * hydromet_dim
       ! Add 1 for "corr_eta_hm_i" to the loop size.
       tot_zt_loops = tot_zt_loops + 1
    endif
    if ( any( vars_zt == "corr_eta_Ncn_i" ) ) then
       ! Correct for number of variables found under "corr_eta_Ncn_i".
       ! Subtract 2 from the loop size (1st PDF comp. and 2nd PDF comp.).
       tot_zt_loops = tot_zt_loops - 2
       ! Add 1 for "corr_eta_Ncn_i" to the loop size.
       tot_zt_loops = tot_zt_loops + 1
    endif
    if ( any( vars_zt == "corr_Ncn_hm_i" ) ) then
       ! Correct for number of variables found under "corr_Ncnhm_i".
       ! Subtract 2 from the loop size (1st PDF component and 2nd PDF component)
       ! for each hydrometeor.
       tot_zt_loops = tot_zt_loops - 2 * hydromet_dim
       ! Add 1 for "corr_Ncnhm_i" to the loop size.
       tot_zt_loops = tot_zt_loops + 1
    endif
    if ( any( vars_zt == "corr_hmx_hmy_i" ) ) then
       ! Correct for number of variables found under "corr_hmxhmy_i".
       ! Subtract 2 (1st PDF component and 2nd PDF component) multipled by the
       ! number of correlations of two hydrometeors, which is found by:
       ! (1/2) * hydromet_dim * ( hydromet_dim - 1 ); from the loop size.
       tot_zt_loops = tot_zt_loops - hydromet_dim * ( hydromet_dim - 1 )
       ! Add 1 for "corr_hmxhmy_i" to the loop size.
       tot_zt_loops = tot_zt_loops + 1
    endif

    if ( any( vars_zt == "corr_w_hm_i_n" ) ) then
       ! Correct for number of variables found under "corr_whm_i_n".
       ! Subtract 2 from the loop size (1st PDF component and 2nd PDF component)
       ! for each hydrometeor.
       tot_zt_loops = tot_zt_loops - 2 * hydromet_dim
       ! Add 1 for "corr_whm_i_n" to the loop size.
       tot_zt_loops = tot_zt_loops + 1
    endif
    if ( any( vars_zt == "corr_w_Ncn_i_n" ) ) then
       ! Correct for number of variables found under "corr_wNcn_i_n".
       ! Subtract 2 from the loop size (1st PDF comp. and 2nd PDF comp.).
       tot_zt_loops = tot_zt_loops - 2
       ! Add 1 for "corr_wNcn_i_n" to the loop size.
       tot_zt_loops = tot_zt_loops + 1
    endif
    if ( any( vars_zt == "corr_chi_hm_i_n" ) ) then
       ! Correct for number of variables found under "corr_chi_hm_i_n".
       ! Subtract 2 from the loop size (1st PDF component and 2nd PDF component)
       ! for each hydrometeor.
       tot_zt_loops = tot_zt_loops - 2 * hydromet_dim
       ! Add 1 for "corr_chi_hm_i_n" to the loop size.
       tot_zt_loops = tot_zt_loops + 1
    endif
    if ( any( vars_zt == "corr_chi_Ncn_i_n" ) ) then
       ! Correct for number of variables found under "corr_chi_Ncn_i_n".
       ! Subtract 2 from the loop size (1st PDF comp. and 2nd PDF comp.).
       tot_zt_loops = tot_zt_loops - 2
       ! Add 1 for "corr_chi_Ncn_i_n" to the loop size.
       tot_zt_loops = tot_zt_loops + 1
    endif
    if ( any( vars_zt == "corr_eta_hm_i_n" ) ) then
       ! Correct for number of variables found under "corr_eta_hm_i_n".
       ! Subtract 2 from the loop size (1st PDF component and 2nd PDF component)
       ! for each hydrometeor.
       tot_zt_loops = tot_zt_loops - 2 * hydromet_dim
       ! Add 1 for "corr_eta_hm_i_n" to the loop size.
       tot_zt_loops = tot_zt_loops + 1
    endif
    if ( any( vars_zt == "corr_eta_Ncn_i_n" ) ) then
       ! Correct for number of variables found under "corr_eta_Ncn_i_n".
       ! Subtract 2 from the loop size (1st PDF comp. and 2nd PDF comp.).
       tot_zt_loops = tot_zt_loops - 2
       ! Add 1 for "corr_eta_Ncn_i_n" to the loop size.
       tot_zt_loops = tot_zt_loops + 1
    endif
    if ( any( vars_zt == "corr_Ncn_hm_i_n" ) ) then
       ! Correct for number of variables found under "corr_Ncnhm_i_n".
       ! Subtract 2 from the loop size (1st PDF component and 2nd PDF component)
       ! for each hydrometeor.
       tot_zt_loops = tot_zt_loops - 2 * hydromet_dim
       ! Add 1 for "corr_Ncnhm_i_n" to the loop size.
       tot_zt_loops = tot_zt_loops + 1
    endif
    if ( any( vars_zt == "corr_hmx_hmy_i_n" ) ) then
       ! Correct for number of variables found under "corr_hmxhmy_i_n".
       ! Subtract 2 (1st PDF component and 2nd PDF component) multipled by the
       ! number of normal space correlations of two hydrometeors, which is found
       ! by:  (1/2) * hydromet_dim * ( hydromet_dim - 1 );
       ! from the loop size.
       tot_zt_loops = tot_zt_loops - hydromet_dim * ( hydromet_dim - 1 )
       ! Add 1 for "corr_hmxhmy_i_n" to the loop size.
       tot_zt_loops = tot_zt_loops + 1
    endif

    if ( any( vars_zt == "hmp2_zt" ) ) then
       ! Correct for number of variables found under "hmp2_zt".
       ! Subtract 1 from the loop size for each hydrometeor.
       tot_zt_loops = tot_zt_loops - hydromet_dim
       ! Add 1 for "hmp2_zt" to the loop size.
       tot_zt_loops = tot_zt_loops + 1
    endif

    if ( any( vars_zt == "wp2hmp" ) ) then
       ! Correct for number of variables found under "wp2hmp".
       ! Subtract 1 from the loop size for each hydrometeor.
       tot_zt_loops = tot_zt_loops - hydromet_dim
       ! Add 1 for "wp2hmp" to the loop size.
       tot_zt_loops = tot_zt_loops + 1
    endif

    if ( any( vars_zt == "sclrm" ) ) then
       ! Correct for number of variables found under "sclrm".
       ! Subtract 1 from the loop size for each scalar.
       tot_zt_loops = tot_zt_loops - sclr_dim

       ! Add 1 for "sclrm" to the loop size.
       tot_zt_loops = tot_zt_loops + 1
    endif

    if ( any( vars_zt == "sclrm_f" ) ) then
       ! Correct for number of variables found under "sclrm_f".
       ! Subtract 1 from the loop size for each scalar.
       tot_zt_loops = tot_zt_loops - sclr_dim
       ! Add 1 for "sclrm_f" to the loop size.
       tot_zt_loops = tot_zt_loops + 1
    endif

    if ( any( vars_zt == "edsclrm" ) ) then
       ! Correct for number of variables found under "edsclrm".
       ! Subtract 1 from the loop size for each scalar.
       tot_zt_loops = tot_zt_loops - edsclr_dim
       ! Add 1 for "edsclrm" to the loop size.
       tot_zt_loops = tot_zt_loops + 1
    endif

    if ( any( vars_zt == "edsclrm_f" ) ) then
       ! Correct for number of variables found under "edsclrm_f".
       ! Subtract 1 from the loop size for each scalar.
       tot_zt_loops = tot_zt_loops - edsclr_dim
       ! Add 1 for "edsclrm_f" to the loop size.
       tot_zt_loops = tot_zt_loops + 1
    endif

    k = 1

    do i = 1, tot_zt_loops

      select case ( trim( vars_zt(i) ) )
      case ('thlm')
        stats_metadata%ithlm = k
        call stat_assign( var_index=stats_metadata%ithlm, var_name="thlm", &
             var_description="thlm, Liquid water potential temperature (theta_l)", var_units="K", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('T_in_K')
        stats_metadata%iT_in_K = k
        call stat_assign( var_index=stats_metadata%iT_in_K, var_name="T_in_K", &
             var_description="T_in_K, Absolute temperature", var_units="K", l_silhs=.false., &
             grid_kind=stats_zt )
        k = k + 1

      case ('thvm')
        stats_metadata%ithvm = k
        call stat_assign( var_index=stats_metadata%ithvm, var_name="thvm", &
             var_description="thvm, Virtual potential temperature", &
             var_units="K", l_silhs=.false., &
             grid_kind=stats_zt )
        k = k + 1

      case ('rtm')
        stats_metadata%irtm = k

        call stat_assign( var_index=stats_metadata%irtm, var_name="rtm", &
             var_description="rtm, Total (vapor+liquid) water mixing ratio", &
             var_units="kg/kg", l_silhs=.false., grid_kind=stats_zt )

        k = k + 1

      case ('rcm')
        stats_metadata%ircm = k
        call stat_assign( var_index=stats_metadata%ircm, var_name="rcm", &
             var_description="rcm, Cloud water mixing ratio", var_units="kg/kg", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rfrzm')
        stats_metadata%irfrzm = k
        call stat_assign( var_index=stats_metadata%irfrzm, var_name="rfrzm", &
             var_description="rfrzm, Total ice phase water mixing ratio", var_units="kg/kg", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rvm')
        stats_metadata%irvm = k
        call stat_assign( var_index=stats_metadata%irvm, var_name="rvm", &
             var_description="rvm, Vapor water mixing ratio", var_units="kg/kg", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1
      case ('rel_humidity')
        stats_metadata%irel_humidity = k
        call stat_assign( var_index=stats_metadata%irel_humidity, var_name="rel_humidity", &
             var_description="rel_humidity, Relative humidity w.r.t. liquid (between 0 and 1)", &
             var_units="[-]", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1
      case ('um')
        stats_metadata%ium = k
        call stat_assign( var_index=stats_metadata%ium, var_name="um", &
             var_description="u_bar, Grid-mean eastward (u) wind", &
             var_units="m/s", l_silhs=.false., &
             grid_kind=stats_zt )
        k = k + 1
      case ('vm')
        stats_metadata%ivm = k
        call stat_assign( var_index=stats_metadata%ivm, var_name="vm", &
             var_description="v_bar, Grid-mean northward (v) wind", &
             var_units="m/s", l_silhs=.false., &
             grid_kind=stats_zt )
        k = k + 1
      case ('wm_zt')
        stats_metadata%iwm_zt = k
        call stat_assign( var_index=stats_metadata%iwm_zt, var_name="wm", &
             var_description="w_bar, Grid-mean upward (w) wind", &
             var_units="m/s", l_silhs=.false., &
             grid_kind=stats_zt )
        k = k + 1
      case ('um_ref')
        stats_metadata%ium_ref = k
        call stat_assign( var_index=stats_metadata%ium_ref, var_name="um_ref", &
             var_description="um_ref, Reference u wind", var_units="m/s", l_silhs=.false., &
             grid_kind=stats_zt )
        k = k + 1
      case ('vm_ref')
        stats_metadata%ivm_ref = k
        call stat_assign( var_index=stats_metadata%ivm_ref, var_name="vm_ref", &
             var_description="vm_ref, Reference v wind", var_units="m/s", l_silhs=.false., &
             grid_kind=stats_zt )
        k = k + 1
      case ('ug')
        stats_metadata%iug = k
        call stat_assign( var_index=stats_metadata%iug, var_name="ug", &
             var_description="ug, u geostrophic wind", var_units="m/s", l_silhs=.false., &
             grid_kind=stats_zt )
        k = k + 1
      case ('vg')
        stats_metadata%ivg = k
        call stat_assign( var_index=stats_metadata%ivg, var_name="vg", &
             var_description="vg, v geostrophic wind", var_units="m/s", l_silhs=.false., &
             grid_kind=stats_zt )
        k = k + 1
      case ('cloud_frac')
        stats_metadata%icloud_frac = k
        call stat_assign( var_index=stats_metadata%icloud_frac, var_name="cloud_frac", &
             var_description="cloud_frac, Cloud fraction (between 0 and 1)", var_units="-", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('ice_supersat_frac')
        stats_metadata%iice_supersat_frac = k
        call stat_assign( var_index=stats_metadata%iice_supersat_frac, var_name="ice_supersat_frac", &
             var_description="ice_supersat_frac, Ice cloud fraction (between 0 and 1)", &
             var_units="count", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rcm_in_layer')
        stats_metadata%ircm_in_layer = k
        call stat_assign( var_index=stats_metadata%ircm_in_layer, var_name="rcm_in_layer", &
             var_description="rcm_in_layer, rcm in cloud layer", &
             var_units="kg/kg", l_silhs=.false., &
             grid_kind=stats_zt )
        k = k + 1

      case ('rcm_in_cloud')
        stats_metadata%ircm_in_cloud = k
        call stat_assign( var_index=stats_metadata%ircm_in_cloud, var_name="rcm_in_cloud", &
             var_description="rcm_in_cloud, In-cloud value of rcm (for microphysics)", &
             var_units="kg/kg", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('cloud_cover')
        stats_metadata%icloud_cover = k
        call stat_assign( var_index=stats_metadata%icloud_cover, var_name="cloud_cover", &
             var_description="cloud_cover, Cloud cover (between 0 and 1)", var_units="count", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1
      case ('p_in_Pa')
        stats_metadata%ip_in_Pa = k
        call stat_assign( var_index=stats_metadata%ip_in_Pa, var_name="p_in_Pa", &
             var_description="p_in_Pa, Pressure", &
             var_units="Pa", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1
      case ('exner')
        stats_metadata%iexner = k
        call stat_assign( var_index=stats_metadata%iexner, var_name="exner", &
             var_description="exner, Exner function = (p/p0)**(rd/cp)", var_units="count", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1
      case ('rho_ds_zt')
        stats_metadata%irho_ds_zt = k
        call stat_assign( var_index=stats_metadata%irho_ds_zt, var_name="rho_ds_zt", &
             var_description="rho_ds_zt, Dry static base-state density", var_units="kg m^{-3}", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1
      case ('thv_ds_zt')
        stats_metadata%ithv_ds_zt = k
        call stat_assign( var_index=stats_metadata%ithv_ds_zt, var_name="thv_ds_zt", &
             var_description="thv_ds_zt, Dry base-state theta_v", &
             var_units="K", l_silhs=.false., &
             grid_kind=stats_zt )
        k = k + 1
      case ('Lscale')
        stats_metadata%iLscale = k
        call stat_assign( var_index=stats_metadata%iLscale, var_name="Lscale", &
          var_description="L, Turbulent mixing length", &
          var_units="m", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1
      case ('thlm_forcing')
        stats_metadata%ithlm_forcing = k
        call stat_assign( var_index=stats_metadata%ithlm_forcing, var_name="thlm_forcing", &
             var_description="thlm_forcing, thlm budget: thetal forcing " &
             // "(includes thlm_mc and radht)",&
             var_units="K s^{-1}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1
      case ('thlm_mc')
        stats_metadata%ithlm_mc = k
        call stat_assign( var_index=stats_metadata%ithlm_mc, var_name="thlm_mc", &
             var_description="thlm_mc, Change in thlm due to microphysics (not in budget)", &
             var_units="K s^{-1}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1
      case ('rtm_forcing')
        stats_metadata%irtm_forcing = k
        call stat_assign( var_index=stats_metadata%irtm_forcing, var_name="rtm_forcing", &
             var_description="rtm_forcing, rtm budget: rt forcing (includes rtm_mc)", &
             var_units="kg kg^{-1} s^{-1}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rtm_mc')
        stats_metadata%irtm_mc = k
        call stat_assign( var_index=stats_metadata%irtm_mc, var_name="rtm_mc", &
             var_description="rtm_mc, Change in rt due to microphysics (not in budget)", &
             var_units="kg kg^{-1} s^{-1}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rvm_mc')
        stats_metadata%irvm_mc = k
        call stat_assign( var_index=stats_metadata%irvm_mc, var_name="rvm_mc", &
             var_description="rvm_mc, Time tendency of vapor mixing ratio due to microphysics", &
             var_units="kg/(kg s)", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rcm_mc')
        stats_metadata%ircm_mc = k
        call stat_assign( var_index=stats_metadata%ircm_mc, var_name="rcm_mc", &
             var_description="rcm_mc, Time tendency of liquid water mixing ratio " &
             // "due microphysics",&
             var_units="kg/kg/s", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rcm_sd_mg_morr')
        stats_metadata%ircm_sd_mg_morr = k
        call stat_assign( var_index=stats_metadata%ircm_sd_mg_morr, var_name="rcm_sd_mg_morr", &
             var_description="rcm_sd_mg_morr, rcm sedimentation when using morrision or MG " &
             // "microphysics (not in budget, included in rcm_mc)", &
             var_units="kg kg^{-1} s^{-1}", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('thlm_mfl_min')
        stats_metadata%ithlm_mfl_min = k
        call stat_assign( var_index=stats_metadata%ithlm_mfl_min, var_name="thlm_mfl_min", &
             var_description="thlm_mfl_min, Minimum allowable thlm", var_units="K", &
             l_silhs=.false., &
             grid_kind=stats_zt )
        k = k + 1

      case ('thlm_mfl_max')
        stats_metadata%ithlm_mfl_max = k
        call stat_assign( var_index=stats_metadata%ithlm_mfl_max, var_name="thlm_mfl_max", &
             var_description="thlm_mfl_max, Maximum allowable thlm", var_units="K", &
             l_silhs=.false., &
             grid_kind=stats_zt )
        k = k + 1

      case ('thlm_enter_mfl')
        stats_metadata%ithlm_enter_mfl = k
        call stat_assign( var_index=stats_metadata%ithlm_enter_mfl, var_name="thlm_enter_mfl", &
             var_description="thlm_enter_mfl, Thlm before flux-limiter", var_units="K", &
             l_silhs=.false., &
             grid_kind=stats_zt )
        k = k + 1

      case ('thlm_exit_mfl')
        stats_metadata%ithlm_exit_mfl = k
        call stat_assign( var_index=stats_metadata%ithlm_exit_mfl, var_name="thlm_exit_mfl", &
             var_description="thlm_exit_mfl, Thlm exiting flux-limiter", var_units="K", &
             l_silhs=.false., &
             grid_kind=stats_zt )
        k = k + 1

      case ('thlm_old')
        stats_metadata%ithlm_old = k
        call stat_assign( var_index=stats_metadata%ithlm_old, var_name="thlm_old", &
             var_description="thlm_old, Thlm at previous timestep", var_units="K", &
             l_silhs=.false., &
             grid_kind=stats_zt )
        k = k + 1

      case ('thlm_without_ta')
        stats_metadata%ithlm_without_ta = k
        call stat_assign( var_index=stats_metadata%ithlm_without_ta, var_name="thlm_without_ta", &
             var_description="thlm_without_ta, Thlm without turbulent advection contribution", &
             var_units="K", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rtm_mfl_min')
        stats_metadata%irtm_mfl_min = k
        call stat_assign( var_index=stats_metadata%irtm_mfl_min, var_name="rtm_mfl_min", &
             var_description="rtm_mfl_min, Minimum allowable rtm", var_units="kg/kg", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rtm_mfl_max')
        stats_metadata%irtm_mfl_max = k
        call stat_assign( var_index=stats_metadata%irtm_mfl_max, var_name="rtm_mfl_max", &
             var_description="rtm_mfl_max, Maximum allowable rtm", var_units="kg/kg", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rtm_enter_mfl')
        stats_metadata%irtm_enter_mfl = k
        call stat_assign( var_index=stats_metadata%irtm_enter_mfl, var_name="rtm_enter_mfl", &
             var_description="rtm_enter_mfl, Rtm before flux-limiter", var_units="kg/kg", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rtm_exit_mfl')
        stats_metadata%irtm_exit_mfl = k
        call stat_assign( var_index=stats_metadata%irtm_exit_mfl, var_name="rtm_exit_mfl", &
             var_description="rtm_exit_mfl, Rtm exiting flux-limiter", var_units="kg/kg", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rtm_old')
        stats_metadata%irtm_old = k
        call stat_assign( var_index=stats_metadata%irtm_old, var_name="rtm_old", &
             var_description="rtm_old, Rtm at previous timestep", var_units="kg/kg", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rtm_without_ta')
        stats_metadata%irtm_without_ta = k
        call stat_assign( var_index=stats_metadata%irtm_without_ta, var_name="rtm_without_ta", &
             var_description="rtm_without_ta, Rtm without turbulent advection contribution", &
             var_units="kg/kg", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('wp3')
        stats_metadata%iwp3 = k
        call stat_assign( var_index=stats_metadata%iwp3, var_name="wp3", &
             var_description="w'^3, Third-order moment of vertical air velocity", &
             var_units="m^3/s^3", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('wpup2')
        stats_metadata%iwpup2 = k
        call stat_assign( var_index=stats_metadata%iwpup2, var_name="wpup2", &
             var_description="w'u'^2, Third-order moment from PDF", &
             var_units="m^3/s^3", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('wpvp2')
        stats_metadata%iwpvp2 = k
        call stat_assign( var_index=stats_metadata%iwpvp2, var_name="wpvp2", &
             var_description="w'v'^2, Third-order moment from PDF", &
             var_units="m^3/s^3", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('thlp3')
        stats_metadata%ithlp3 = k
        call stat_assign( var_index=stats_metadata%ithlp3, var_name="thlp3", &
             var_description="thl'^3, Third-order moment of theta_l", var_units="K^3", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rtp3')
        stats_metadata%irtp3 = k
        call stat_assign( var_index=stats_metadata%irtp3, var_name="rtp3", &
             var_description="rt'^3, Third-order moment of total water, rt", var_units="(kg/kg)^3", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('wpthlp2')
        stats_metadata%iwpthlp2 = k
        call stat_assign( var_index=stats_metadata%iwpthlp2, var_name="wpthlp2", &
             var_description="w'theta_l'^2, Vertical turbulent flux of theta_l'^2", &
             var_units="(m K^2)/s", l_silhs=.false., &
             grid_kind=stats_zt )
        k = k + 1

      case ('wp2thlp')
        stats_metadata%iwp2thlp = k
        call stat_assign( var_index=stats_metadata%iwp2thlp, var_name="wp2thlp", &
             var_description="w'^2theta_l', Vertical turbulent flux of w'theta_l'", &
             var_units="(m^2 K)/s^2", l_silhs=.false., &
             grid_kind=stats_zt )
        k = k + 1

      case ('wprtp2')
        stats_metadata%iwprtp2 = k
        call stat_assign( var_index=stats_metadata%iwprtp2, var_name="wprtp2", &
             var_description="w'rt'^2, Vertical turbulent flux of rt'^2", &
             var_units="(m kg)/(s kg)", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('wp2rtp')
        stats_metadata%iwp2rtp = k
        call stat_assign( var_index=stats_metadata%iwp2rtp, var_name="wp2rtp", &
             var_description="w'^2rt', Vertical turbulent flux of w'rt'", &
             var_units="(m^2 kg)/(s^2 kg)", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('Lscale_up')
        stats_metadata%iLscale_up = k
        call stat_assign( var_index=stats_metadata%iLscale_up, var_name="Lscale_up", &
             var_description="Lscale_up, Upward mixing length", var_units="m", l_silhs=.false., &
             grid_kind=stats_zt )
        k = k + 1

      case ('Lscale_down')
        stats_metadata%iLscale_down = k
        call stat_assign( var_index=stats_metadata%iLscale_down, var_name="Lscale_down", &
             var_description="Lscale_down, Downward mixing length", var_units="m", &
             l_silhs=.false.,&
             grid_kind=stats_zt )
        k = k + 1

      case ('Lscale_pert_1')
        stats_metadata%iLscale_pert_1 = k
        call stat_assign( var_index=stats_metadata%iLscale_pert_1, var_name="Lscale_pert_1", &
             var_description="Lscale_pert_1, Mixing length using a perturbed value of rtm/thlm", &
             var_units="m", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('Lscale_pert_2')
        stats_metadata%iLscale_pert_2 = k
        call stat_assign( var_index=stats_metadata%iLscale_pert_2, var_name="Lscale_pert_2", &
             var_description="Lscale_pert_2, Mixing length using a perturbed value of rtm/thlm", &
             var_units="m", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('tau_zt')
        stats_metadata%itau_zt = k
        call stat_assign( var_index=stats_metadata%itau_zt, var_name="tau_zt", &
             var_description="tau_zt, Dissipation time", var_units="s", l_silhs=.false., &
             grid_kind=stats_zt )
        k = k + 1

      case ('invrs_tau_zt')
        stats_metadata%iinvrs_tau_zt = k
        call stat_assign( var_index=stats_metadata%iinvrs_tau_zt, var_name="invrs_tau_zt", &
             var_description="invrs_tau_zt, Inverse of dissipation time", var_units="s^-1", & 
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('Kh_zt')
        stats_metadata%iKh_zt = k
        call stat_assign( var_index=stats_metadata%iKh_zt, var_name="Kh_zt", &
             var_description="Kh_zt, Eddy diffusivity", var_units="m^2/s", l_silhs=.false., &
             grid_kind=stats_zt )
        k = k + 1

      case ('wp2thvp')
        stats_metadata%iwp2thvp = k
        call stat_assign( var_index=stats_metadata%iwp2thvp, var_name="wp2thvp", &
             var_description="w'^2theta_v', Vertical turbulent flux of w'theta_v'", &
             var_units="K m^2/s^2", l_silhs=.false., &
             grid_kind=stats_zt )
        k = k + 1

      case ('wp2rcp')
        stats_metadata%iwp2rcp = k
        call stat_assign( var_index=stats_metadata%iwp2rcp, var_name="wp2rcp", &
             var_description="w'^2rc'", var_units="(m^2 kg)/(s^2 kg)", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1
        
      case ('w_up_in_cloud')
        stats_metadata%iw_up_in_cloud = k
        call stat_assign( var_index=stats_metadata%iw_up_in_cloud, var_name="w_up_in_cloud", &
             var_description="Mean W in saturated updrafts", &
             var_units="m/s", l_silhs=.false., grid_kind=stats_zt )
         k = k + 1

      case ('w_down_in_cloud')
        stats_metadata%iw_down_in_cloud = k
        call stat_assign( var_index=stats_metadata%iw_down_in_cloud, var_name="w_down_in_cloud", &
             var_description="Mean W in saturated downdrafts", &
             var_units="m/s", l_silhs=.false., grid_kind=stats_zt )
         k = k + 1

      case ('cld_updr_frac')
        stats_metadata%icld_updr_frac = k
        call stat_assign( var_index=stats_metadata%icld_updr_frac, var_name="cld_updr_frac", &
             var_description="Cloudy Updraft Fraction", &
             var_units="-", l_silhs=.false., grid_kind=stats_zt )
         k = k + 1

      case ('cld_downdr_frac')
        stats_metadata%icld_downdr_frac = k
        call stat_assign( var_index=stats_metadata%icld_downdr_frac, var_name="cld_downdr_frac", &
             var_description="Cloudy Downdraft Fraction", &
             var_units="-", l_silhs=.false., grid_kind=stats_zt )
         k = k + 1

      case ('wprtpthlp')
        stats_metadata%iwprtpthlp = k
        call stat_assign( var_index=stats_metadata%iwprtpthlp, var_name="wprtpthlp", &
             var_description="w'rt'theta_l', Vertical turbulent flux of rt'theta_l'", &
             var_units="(m kg K)/(s kg)", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rc_coef')
        stats_metadata%irc_coef = k
        call stat_assign( var_index=stats_metadata%irc_coef, var_name="rc_coef", &
             var_description="rc_coef, Coefficient of X'r_c'", &
             var_units="K/(kg/kg)", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('sigma_sqd_w_zt')
        stats_metadata%isigma_sqd_w_zt = k
        call stat_assign( var_index=stats_metadata%isigma_sqd_w_zt, var_name="sigma_sqd_w_zt", &
             var_description="sigma_sqd_w_zt, Nondimensionalized w variance of " &
             // "Gaussian component",&
             var_units="-", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rho')
        stats_metadata%irho = k
        call stat_assign( var_index=stats_metadata%irho, var_name="rho", var_description="rho, Air density", &
             var_units="kg m^{-3}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('Ncm')           ! Brian
        stats_metadata%iNcm = k
        call stat_assign( var_index=stats_metadata%iNcm, var_name="Ncm", &
             var_description="Ncm, Cloud droplet number concentration", var_units="num/kg", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('Nc_in_cloud')
        stats_metadata%iNc_in_cloud = k

        call stat_assign( var_index=stats_metadata%iNc_in_cloud, var_name="Nc_in_cloud", &
             var_description="Nc_in_cloud, In cloud droplet concentration", var_units="num/kg", &
             l_silhs=.false., grid_kind=stats_zt )

        k = k + 1

      case ('Nc_activated')
        stats_metadata%iNc_activated = k

        call stat_assign( var_index=stats_metadata%iNc_activated, var_name="Nc_activated", &
             var_description="Nc_activated, Droplets activated by GFDL activation", &
             var_units="num/kg", l_silhs=.false., grid_kind=stats_zt )

        k = k + 1

      case ('Nccnm')
        stats_metadata%iNccnm = k
        call stat_assign( var_index=stats_metadata%iNccnm, var_name="Nccnm", &
             var_description="Nccnm, Cloud condensation nuclei concentration", &
             var_units="num/kg", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('Nim')           ! Brian
        stats_metadata%iNim = k
        call stat_assign( var_index=stats_metadata%iNim, var_name="Nim", &
             var_description="Nim, Ice crystal number concentration", var_units="num/kg", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('snowslope')     ! Adam Smith, 22 April 2008
        stats_metadata%isnowslope = k
        call stat_assign( var_index=stats_metadata%isnowslope, var_name="snowslope", &
             var_description="snowslope, COAMPS microphysics snow slope parameter", &
             var_units="1/m", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('Nsm')        ! Adam Smith, 22 April 2008
        stats_metadata%iNsm = k
        call stat_assign( var_index=stats_metadata%iNsm, var_name="Nsm", &
             var_description="Nsm, Snow particle number concentration", var_units="num/kg", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('Ngm')
        stats_metadata%iNgm = k
        call stat_assign( var_index=stats_metadata%iNgm, var_name="Ngm", &
             var_description="Ngm, Graupel number concentration", var_units="num/kg", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('sed_rcm')       ! Brian
        stats_metadata%ised_rcm = k
        call stat_assign( var_index=stats_metadata%ised_rcm, var_name="sed_rcm", &
             var_description="sed_rcm, d(rcm)/dt due to cloud sedimentation", &
             var_units="kg / [m^2 s]", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rsat')           ! Brian
        stats_metadata%irsat = k
        call stat_assign( var_index=stats_metadata%irsat, var_name="rsat", &
             var_description="rsat, Saturation mixing ratio over liquid", var_units="kg/kg", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rsati')
        stats_metadata%irsati = k
        call stat_assign( var_index=stats_metadata%irsati, var_name="rsati", &
             var_description="rsati, Saturation mixing ratio over ice", var_units="kg/kg", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rrm')           ! Brian
        stats_metadata%irrm = k
        call stat_assign( var_index=stats_metadata%irrm, var_name="rrm", &
             var_description="rrm, Rain water mixing ratio", var_units="kg/kg", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rsm')
        stats_metadata%irsm = k
        call stat_assign( var_index=stats_metadata%irsm, var_name="rsm", &
             var_description="rsm, Snow water mixing ratio", var_units="kg/kg", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rim')
        stats_metadata%irim = k
        call stat_assign( var_index=stats_metadata%irim, var_name="rim", &
             var_description="rim, Pristine ice water mixing ratio", var_units="kg/kg", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rgm')
        stats_metadata%irgm = k
        call stat_assign( var_index=stats_metadata%irgm, var_name="rgm", &
             var_description="rgm, Graupel water mixing ratio", var_units="kg/kg", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('Nrm')           ! Brian
        stats_metadata%iNrm = k
        call stat_assign( var_index=stats_metadata%iNrm, var_name="Nrm", &
             var_description="Nrm, Rain drop number concentration", var_units="num/kg", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('m_vol_rad_rain')  ! Brian
        stats_metadata%im_vol_rad_rain = k
        call stat_assign( var_index=stats_metadata%im_vol_rad_rain, var_name="mvrr", &
             var_description="mvrr, Rain drop mean volume radius", &
             var_units="m", l_silhs=.false.,&
             grid_kind=stats_zt )
        k = k + 1

      case ('m_vol_rad_cloud')
        stats_metadata%im_vol_rad_cloud = k
        call stat_assign( var_index=stats_metadata%im_vol_rad_cloud, var_name="m_vol_rad_cloud", &
             var_description="m_vol_rad_cloud, Cloud drop mean volume radius", var_units="m", &
             l_silhs=.false., &
             grid_kind=stats_zt )
        k = k + 1

      case ('eff_rad_cloud')
        stats_metadata%ieff_rad_cloud = k
        call stat_assign( var_index=stats_metadata%ieff_rad_cloud, var_name="eff_rad_cloud", &
             var_description="eff_rad_cloud, Cloud drop effective volume radius", &
             var_units="microns", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('eff_rad_ice')
        stats_metadata%ieff_rad_ice = k

        call stat_assign( var_index=stats_metadata%ieff_rad_ice, var_name="eff_rad_ice", &
             var_description="eff_rad_ice, Ice effective volume radius", var_units="microns", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('eff_rad_snow')
        stats_metadata%ieff_rad_snow = k
        call stat_assign( var_index=stats_metadata%ieff_rad_snow, var_name="eff_rad_snow", &
             var_description="eff_rad_snow, Snow effective volume radius", &
             var_units="microns", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('eff_rad_rain')
        stats_metadata%ieff_rad_rain = k
        call stat_assign( var_index=stats_metadata%ieff_rad_rain, var_name="eff_rad_rain", &
             var_description="eff_rad_rain, Rain drop effective volume radius", &
             var_units="microns", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('eff_rad_graupel')
        stats_metadata%ieff_rad_graupel = k
        call stat_assign( var_index=stats_metadata%ieff_rad_graupel, var_name="eff_rad_graupel", &
             var_description="eff_rad_graupel, Graupel effective volume radius", &
             var_units="microns", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('precip_rate_zt')     ! Brian
        stats_metadata%iprecip_rate_zt = k

        call stat_assign( var_index=stats_metadata%iprecip_rate_zt, var_name="precip_rate_zt", &
             var_description="precip_rate_zt, Rain rate", var_units="mm/day", l_silhs=.false., &
             grid_kind=stats_zt )
        k = k + 1

      case ('radht')
        stats_metadata%iradht = k

        call stat_assign( var_index=stats_metadata%iradht, var_name="radht", &
             var_description="radht, Total (sw+lw) radiative heating rate", var_units="K/s", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('radht_LW')
        stats_metadata%iradht_LW = k

        call stat_assign( var_index=stats_metadata%iradht_LW, var_name="radht_LW", &
             var_description="radht_LW, Long-wave radiative heating rate", var_units="K/s", &
             l_silhs=.false., grid_kind=stats_zt )

        k = k + 1

      case ('radht_SW')
        stats_metadata%iradht_SW = k
        call stat_assign( var_index=stats_metadata%iradht_SW, var_name="radht_SW", &
             var_description="radht_SW, Short-wave radiative heating rate", var_units="K/s", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('diam')
        stats_metadata%idiam = k

        call stat_assign( var_index=stats_metadata%idiam, var_name="diam", &
             var_description="diam, Ice crystal diameter", var_units="m", l_silhs=.false., &
             grid_kind=stats_zt )
        k = k + 1

      case ('mass_ice_cryst')
        stats_metadata%imass_ice_cryst = k
        call stat_assign( var_index=stats_metadata%imass_ice_cryst, var_name="mass_ice_cryst", &
             var_description="mass_ice_cryst, Mass of a single ice crystal", var_units="kg", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rcm_icedfs')
        stats_metadata%ircm_icedfs = k
        call stat_assign( var_index=stats_metadata%ircm_icedfs, var_name="rcm_icedfs", &
             var_description="rcm_icedfs, Change in liquid due to ice", var_units="kg/kg/s", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('u_T_cm')
        stats_metadata%iu_T_cm = k
        call stat_assign( var_index=stats_metadata%iu_T_cm, var_name="u_T_cm", &
             var_description="u_T_cm, Ice crystal fallspeed", var_units="cm s^{-1}", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rtm_bt')
        stats_metadata%irtm_bt = k

        call stat_assign( var_index=stats_metadata%irtm_bt, var_name="rtm_bt", &
             var_description="rtm_bt, rtm budget: rtm time tendency", &
             var_units="kg kg^{-1} s^{-1}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rtm_ma')
        stats_metadata%irtm_ma = k

        call stat_assign( var_index=stats_metadata%irtm_ma, var_name="rtm_ma", &
             var_description="rtm_ma, rtm budget: rtm vertical mean advection", &
             var_units="kg kg^{-1} s^{-1}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rtm_ta')
        stats_metadata%irtm_ta = k

        call stat_assign( var_index=stats_metadata%irtm_ta, var_name="rtm_ta", &
             var_description="rtm_ta, rtm budget: rtm turbulent advection", &
             var_units="kg kg^{-1} s^{-1}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rtm_mfl')
        stats_metadata%irtm_mfl = k

        call stat_assign( var_index=stats_metadata%irtm_mfl, var_name="rtm_mfl", &
         var_description="rtm_mfl, rtm budget: rtm correction due to monotonic flux limiter", &
         var_units="kg kg^{-1} s^{-1}", l_silhs=.false., grid_kind=stats_zt)
        k = k + 1

      case ('rtm_tacl')
        stats_metadata%irtm_tacl = k

        call stat_assign( var_index=stats_metadata%irtm_tacl, var_name="rtm_tacl", &
          var_description="rtm_tacl, rtm budget: rtm correction due to ta term" &
          // " (wprtp) clipping", &
          var_units="kg kg^{-1} s^{-1}", l_silhs=.false., grid_kind=stats_zt)

        k = k + 1

      case ('rtm_cl')
        stats_metadata%irtm_cl = k

        call stat_assign( var_index=stats_metadata%irtm_cl, var_name="rtm_cl", &
             var_description="rtm_cl, rtm budget: rtm clipping", &
             var_units="kg kg^{-1} s^{-1}", l_silhs=.false., grid_kind=stats_zt )

        k = k + 1
      case ('rtm_sdmp')
        stats_metadata%irtm_sdmp = k

        call stat_assign( var_index=stats_metadata%irtm_sdmp, var_name="rtm_sdmp", &
             var_description="rtm_sdmp, rtm budget: rtm correction due to sponge damping", &
             var_units="kg kg^{-1} s^{-1}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1


      case ('rtm_pd')
        stats_metadata%irtm_pd = k

        call stat_assign( var_index=stats_metadata%irtm_pd, var_name="rtm_pd", &
             var_description="rtm_pd, rtm budget: rtm positive definite adjustment", &
             var_units="kg kg^{-1} s^{-1}", l_silhs=.false., grid_kind=stats_zt )

        k = k + 1

      case ('thlm_bt')
        stats_metadata%ithlm_bt = k

        call stat_assign( var_index=stats_metadata%ithlm_bt, var_name="thlm_bt", &
             var_description="thlm_bt, thlm budget: thlm time tendency", var_units="K s^{-1}", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('thlm_ma')
        stats_metadata%ithlm_ma = k

        call stat_assign( var_index=stats_metadata%ithlm_ma, var_name="thlm_ma", &
             var_description="thlm_ma, thlm budget: thlm vertical mean advection", &
             var_units="K s^{-1}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('thlm_sdmp')
        stats_metadata%ithlm_sdmp = k

        call stat_assign( var_index=stats_metadata%ithlm_sdmp, var_name="thlm_sdmp", &
             var_description="thlm_sdmp, thlm budget: thlm correction due to sponge damping", &
             var_units="K s^{-1}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1


      case ('thlm_ta')
        stats_metadata%ithlm_ta = k

        call stat_assign( var_index=stats_metadata%ithlm_ta, var_name="thlm_ta", &
             var_description="thlm_ta, thlm budget: thlm turbulent advection", &
             var_units="K s^{-1}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('thlm_mfl')
        stats_metadata%ithlm_mfl = k

        call stat_assign( var_index=stats_metadata%ithlm_mfl, var_name="thlm_mfl", &
             var_description="thlm_mfl, thlm budget: thlm correction due to monotonic " &
             // "flux limiter", &
             var_units="K s^{-1}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('thlm_tacl')
        stats_metadata%ithlm_tacl = k

        call stat_assign( var_index=stats_metadata%ithlm_tacl, var_name="thlm_tacl", &
             var_description="thlm_tacl, thlm budget: thlm correction due to ta term " &
             // "(wpthlp) clipping", &
             var_units="K s^{-1}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('thlm_cl')
        stats_metadata%ithlm_cl = k

        call stat_assign( var_index=stats_metadata%ithlm_cl, var_name="thlm_cl", &
             var_description="thlm_cl, thlm budget: thlm_cl", var_units="K s^{-1}", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('wp3_bt')
        stats_metadata%iwp3_bt = k

        call stat_assign( var_index=stats_metadata%iwp3_bt, var_name="wp3_bt", &
             var_description="wp3_bt, wp3 budget: wp3 time tendency", &
             var_units="m^{3} s^{-4}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('wp3_ma')
        stats_metadata%iwp3_ma = k

        call stat_assign( var_index=stats_metadata%iwp3_ma, var_name="wp3_ma", &
             var_description="wp3_ma, wp3 budget: wp3 vertical mean advection", &
             var_units="m^{3} s^{-4}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('wp3_ta')
        stats_metadata%iwp3_ta = k

        call stat_assign( var_index=stats_metadata%iwp3_ta, var_name="wp3_ta", &
             var_description="wp3_ta, wp3 budget: wp3 turbulent advection", &
             var_units="m^{3} s^{-4}", l_silhs=.false., grid_kind=stats_zt )

        k = k + 1

      case ('wp3_tp')
        stats_metadata%iwp3_tp = k
        call stat_assign( var_index=stats_metadata%iwp3_tp, var_name="wp3_tp", &
             var_description="wp3_tp, wp3 budget: wp3 turbulent transport", &
             var_units="m^{3} s^{-4}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('wp3_ac')
        stats_metadata%iwp3_ac = k
        call stat_assign( var_index=stats_metadata%iwp3_ac, var_name="wp3_ac", &
             var_description="wp3_ac, wp3 budget: wp3 accumulation term", &
             var_units="m^{3} s^{-4}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('wp3_bp1')
        stats_metadata%iwp3_bp1 = k
        call stat_assign( var_index=stats_metadata%iwp3_bp1, var_name="wp3_bp1", &
             var_description="wp3_bp1, wp3 budget: wp3 buoyancy production", &
             var_units="m^{3} s^{-4}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('wp3_pr_tp')
        stats_metadata%iwp3_pr_tp = k
        call stat_assign( var_index=stats_metadata%iwp3_pr_tp, var_name="wp3_pr_tp", &
             var_description= &
               "wp3_pr_tp, wp3 budget: wp3 pressure damping of turbulent production", &
             var_units="m^{3} s^{-4}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('wp3_pr_turb')
        stats_metadata%iwp3_pr_turb = k
        call stat_assign( var_index=stats_metadata%iwp3_pr_turb, var_name="wp3_pr_turb", &
             var_description="wp3_pr_turb, wp3 budget: wp3 pressure-turbulence correlation term", &
             var_units="m^{3} s^{-4}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('wp3_pr_dfsn')
        stats_metadata%iwp3_pr_dfsn = k
        call stat_assign( var_index=stats_metadata%iwp3_pr_dfsn, var_name="wp3_pr_dfsn", &
             var_description="wp3_pr_dfsn, wp3 budget: wp3 pressure diffusion term", &
             var_units="m^{3} s^{-4}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('wp3_pr1')
        stats_metadata%iwp3_pr1 = k
        call stat_assign( var_index=stats_metadata%iwp3_pr1, var_name="wp3_pr1", &
             var_description="wp3_pr1, wp3 budget: wp3 pressure term 1", &
             var_units="m^{3} s^{-4}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('wp3_pr2')
        stats_metadata%iwp3_pr2 = k
        call stat_assign( var_index=stats_metadata%iwp3_pr2, var_name="wp3_pr2", &
             var_description="wp3_pr2, wp3 budget: wp3 pressure term 2", &
             var_units="m^{3} s^{-4}", l_silhs=.false., grid_kind=stats_zt )

        k = k + 1

      case ('wp3_pr3')
        stats_metadata%iwp3_pr3 = k
        call stat_assign( var_index=stats_metadata%iwp3_pr3, var_name="wp3_pr3", &
             var_description="wp3_pr3, wp3 budget: wp3 pressure term 3", &
             var_units="m^{3} s^{-4}", l_silhs=.false., grid_kind=stats_zt )

        k = k + 1

      case ('wp3_dp1')
        stats_metadata%iwp3_dp1 = k
        call stat_assign( var_index=stats_metadata%iwp3_dp1, var_name="wp3_dp1", &
             var_description="wp3_dp1, wp3 budget: wp3 dissipation term 1", &
             var_units="m^{3} s^{-4}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('wp3_sdmp')
        stats_metadata%iwp3_sdmp = k
        call stat_assign( var_index=stats_metadata%iwp3_sdmp, var_name="wp3_sdmp", &
             var_description="wp3_sdmp, wp3 budget: wp3 sponge damping term", &
             var_units="m^{3} s^{-4}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('wp3_cl')
        stats_metadata%iwp3_cl = k
        call stat_assign( var_index=stats_metadata%iwp3_cl, var_name="wp3_cl", &
             var_description="wp3_cl, wp3 budget: wp3 clipping term", &
             var_units="m^{3} s^{-4}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('wp3_splat')
        stats_metadata%iwp3_splat = k
        call stat_assign( var_index=stats_metadata%iwp3_splat, var_name="wp3_splat", &
             var_description="wp3_splat, wp3 budget: wp3 splatting term", &
             var_units="m^3 s^-4", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rtp3_bt')
        stats_metadata%irtp3_bt = k

        call stat_assign( var_index=stats_metadata%irtp3_bt, var_name="rtp3_bt", &
             var_description="rtp3_bt, rtp3 budget: rtp3 time tendency" &
                             // "[kg^{3} kg^{-3} s^{-1}]", &
             var_units="kg^{3} kg^{-3} s^{-1}", l_silhs=.false., &
             grid_kind=stats_zt )
        k = k + 1

      case ('rtp3_tp')
        stats_metadata%irtp3_tp = k

        call stat_assign( var_index=stats_metadata%irtp3_tp, var_name="rtp3_tp", &
             var_description="rtp3_tp, rtp3 budget: rtp3 turbulent production" &
                             // "[kg^{3} kg^{-3} s^{-1}]", &
             var_units="kg^{3} kg^{-3} s^{-1}", l_silhs=.false., &
             grid_kind=stats_zt )
        k = k + 1

      case ('rtp3_ac')
        stats_metadata%irtp3_ac = k

        call stat_assign( var_index=stats_metadata%irtp3_ac, var_name="rtp3_ac", &
             var_description="rtp3_ac, rtp3 budget: rtp3 accumulation" &
                             // "[kg^{3} kg^{-3} s^{-1}]", &
             var_units="kg^{3} kg^{-3} s^{-1}", l_silhs=.false., &
             grid_kind=stats_zt )
        k = k + 1

      case ('rtp3_dp')
        stats_metadata%irtp3_dp = k

        call stat_assign( var_index=stats_metadata%irtp3_dp, var_name="rtp3_dp", &
             var_description="rtp3_dp, rtp3 budget: rtp3 dissipation" &
                             // "[kg^{3} kg^{-3} s^{-1}]", &
             var_units="kg^{3} kg^{-3} s^{-1}", l_silhs=.false., &
             grid_kind=stats_zt )
        k = k + 1

      case ('thlp3_bt')
        stats_metadata%ithlp3_bt = k

        call stat_assign( var_index=stats_metadata%ithlp3_bt, var_name="thlp3_bt", &
             var_description="thlp3_bt, thlp3 budget: thlp3 time tendency" &
                             // "[K^{3} s^{-1}]", &
             var_units="K^{3} s^{-1}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('thlp3_tp')
        stats_metadata%ithlp3_tp = k

        call stat_assign( var_index=stats_metadata%ithlp3_tp, var_name="thlp3_tp", &
             var_description="thlp3_tp, thlp3 budget: thlp3 turbulent production" &
                             // "[K^{3} s^{-1}]", &
             var_units="K^{3} s^{-1}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('thlp3_ac')
        stats_metadata%ithlp3_ac = k

        call stat_assign( var_index=stats_metadata%ithlp3_ac, var_name="thlp3_ac", &
             var_description="thlp3_ac, thlp3 budget: thlp3 accumulation" &
                             // "[K^{3} s^{-1}]", &
             var_units="K^{3} s^{-1}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('thlp3_dp')
        stats_metadata%ithlp3_dp = k

        call stat_assign( var_index=stats_metadata%ithlp3_dp, var_name="thlp3_dp", &
             var_description="thlp3_dp, thlp3 budget: thlp3 dissipation", &
             var_units="K^{3} s^{-1}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rrm_bt')
        stats_metadata%irrm_bt = k
        call stat_assign( var_index=stats_metadata%irrm_bt, var_name="rrm_bt", &
             var_description="rrm_bt, rrm budget: rrm time tendency", &
             var_units="kg kg^{-1} s^{-1}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rrm_ma')
        stats_metadata%irrm_ma = k

        call stat_assign( var_index=stats_metadata%irrm_ma, var_name="rrm_ma", &
             var_description="rrm_ma, rrm budget: rrm vertical mean advection", &
             var_units="kg kg^{-1} s^{-1}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rrm_sd')
        stats_metadata%irrm_sd = k

        call stat_assign( var_index=stats_metadata%irrm_sd, var_name="rrm_sd", &
             var_description="rrm_sd, rrm budget: rrm sedimentation", &
             var_units="kg kg^{-1} s^{-1}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rrm_ts')
        stats_metadata%irrm_ts = k

        call stat_assign( var_index=stats_metadata%irrm_ts, var_name="rrm_ts", &
             var_description="rrm_ts, rrm budget: rrm turbulent sedimentation", &
             var_units="kg kg^{-1} s^{-1}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rrm_sd_morr')
        stats_metadata%irrm_sd_morr = k

        call stat_assign( var_index=stats_metadata%irrm_sd_morr, var_name="rrm_sd_morr", &
             var_description="rrm_sd_morr, rrm sedimentation when using morrision microphysics &
             &(not in budget, included in rrm_mc)", &
             var_units="kg kg^{-1} s^{-1}", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('rrm_ta')
        stats_metadata%irrm_ta = k

        call stat_assign( var_index=stats_metadata%irrm_ta, var_name="rrm_ta", &
             var_description="rrm_ta, rrm budget: rrm turbulent advection", &
             var_units="kg kg^{-1} s^{-1}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rrm_evap')
        stats_metadata%irrm_evap = k

        call stat_assign( var_index=stats_metadata%irrm_evap, var_name="rrm_evap", &
             var_description="rrm_evap, rrm evaporation rate", &
             var_units="kg kg^{-1} s^{-1}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rrm_auto')
        stats_metadata%irrm_auto = k

        call stat_assign( var_index=stats_metadata%irrm_auto, var_name="rrm_auto", &
             var_description="rrm_auto, rrm autoconversion rate", &
             var_units="kg kg^{-1} s^{-1}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rrm_accr')
        stats_metadata%irrm_accr = k
        call stat_assign( var_index=stats_metadata%irrm_accr, var_name="rrm_accr", &
             var_description="rrm_accr, rrm accretion rate", &
             var_units="kg kg^{-1} s^{-1}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rrm_evap_adj')
        stats_metadata%irrm_evap_adj = k

        call stat_assign( var_index=stats_metadata%irrm_evap_adj, var_name="rrm_evap_adj", &
             var_description="rrm_evap_adj, rrm evaporation adjustment due to over-evaporation", &
             var_units="kg kg^{-1} s^{-1}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rrm_src_adj')
        stats_metadata%irrm_src_adj = k

        call stat_assign( var_index=stats_metadata%irrm_src_adj, var_name="rrm_src_adj", &
             var_description="rrm_src_adj, rrm source term adjustment due to over-depletion", &
             var_units="kg kg^{-1} s^{-1}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rrm_mc_nonadj')
        stats_metadata%irrm_mc_nonadj = k

        call stat_assign( var_index=stats_metadata%irrm_mc_nonadj, var_name="rrm_mc_nonadj", &
             var_description="rrm_mc_nonadj, Value of rrm_mc tendency before adjustment", &
             var_units="kg kg^{-1} s^{-1}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rrm_hf')
        stats_metadata%irrm_hf = k
        call stat_assign( var_index=stats_metadata%irrm_hf, var_name="rrm_hf", &
             var_description="rrm_hf, rrm budget: rrm hole-filling term", &
             var_units="kg kg^{-1} s^{-1}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rrm_wvhf')
        stats_metadata%irrm_wvhf = k
        call stat_assign( var_index=stats_metadata%irrm_wvhf, var_name="rrm_wvhf", &
             var_description="rrm_wvhf, rrm budget: rrm water vapor hole-filling term", &
             var_units="kg kg^{-1} s^{-1}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rrm_cl')
        stats_metadata%irrm_cl = k
        call stat_assign( var_index=stats_metadata%irrm_cl, var_name="rrm_cl", &
             var_description="rrm_cl, rrm budget: rrm clipping term", &
             var_units="kg kg^{-1} s^{-1}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rrm_mc')
        stats_metadata%irrm_mc = k

        call stat_assign( var_index=stats_metadata%irrm_mc, var_name="rrm_mc", &
             var_description="rrm_mc, rrm budget: Change in rrm due to microphysics", &
             var_units="kg kg^{-1} s^{-1}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('Nrm_bt')
        stats_metadata%iNrm_bt = k
        call stat_assign( var_index=stats_metadata%iNrm_bt, var_name="Nrm_bt", &
             var_description="Nrm_bt, Nrm budget: Nrm time tendency", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=stats_zt )

        k = k + 1

      case ('Nrm_ma')
        stats_metadata%iNrm_ma = k

        call stat_assign( var_index=stats_metadata%iNrm_ma, var_name="Nrm_ma", &
             var_description="Nrm_ma, Nrm budget: Nrm vertical mean advection", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('Nrm_sd')
        stats_metadata%iNrm_sd = k

        call stat_assign( var_index=stats_metadata%iNrm_sd, var_name="Nrm_sd", &
             var_description="Nrm_sd, Nrm budget: Nrm sedimentation", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=stats_zt )

        k = k + 1

      case ('Nrm_ts')
        stats_metadata%iNrm_ts = k

        call stat_assign( var_index=stats_metadata%iNrm_ts, var_name="Nrm_ts", &
             var_description="Nrm_ts, Nrm budget: Nrm turbulent sedimentation", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('Nrm_ta')
        stats_metadata%iNrm_ta = k
        call stat_assign( var_index=stats_metadata%iNrm_ta, var_name="Nrm_ta", &
             var_description="Nrm_ta, Nrm budget: Nrm turbulent advection", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=stats_zt )

        k = k + 1

      case ('Nrm_evap')
        stats_metadata%iNrm_evap = k

        call stat_assign( var_index=stats_metadata%iNrm_evap, var_name="Nrm_evap", &
             var_description="Nrm_evap, Nrm evaporation rate", var_units="(num/kg)/s", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('Nrm_auto')
        stats_metadata%iNrm_auto = k

        call stat_assign( var_index=stats_metadata%iNrm_auto, var_name="Nrm_auto", &
             var_description="Nrm_auto, Nrm autoconversion rate", var_units="(num/kg)/s", &
             l_silhs=.false., grid_kind=stats_zt )

        k = k + 1

      case ('Nrm_evap_adj')
        stats_metadata%iNrm_evap_adj = k

        call stat_assign( var_index=stats_metadata%iNrm_evap_adj, var_name="Nrm_evap_adj", &
             var_description="Nrm_evap_adj, Nrm evaporation adjustment due to over-evaporation", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('Nrm_src_adj')
        stats_metadata%iNrm_src_adj = k

        call stat_assign( var_index=stats_metadata%iNrm_src_adj, var_name="Nrm_src_adj", &
             var_description="Nrm_src_adj, Nrm source term adjustment due to over-depletion", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('Nrm_cl')
        stats_metadata%iNrm_cl = k
        call stat_assign( var_index=stats_metadata%iNrm_cl, var_name="Nrm_cl", &
             var_description="Nrm_cl, Nrm budget: Nrm clipping term", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('Nrm_mc')
        stats_metadata%iNrm_mc = k
        call stat_assign( var_index=stats_metadata%iNrm_mc, var_name="Nrm_mc", &
             var_description="Nrm_mc, Nrm budget: Change in Nrm due to microphysics " &
             // "(Not in budget)", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rsm_bt')
        stats_metadata%irsm_bt = k
        call stat_assign( var_index=stats_metadata%irsm_bt, var_name="rsm_bt", &
             var_description="rsm_bt, rsm budget: rsm time tendency", &
             var_units="(kg/kg)/s", l_silhs=.false., grid_kind=stats_zt )

        k = k + 1

      case ('rsm_ma')
        stats_metadata%irsm_ma = k

        call stat_assign( var_index=stats_metadata%irsm_ma, var_name="rsm_ma", &
             var_description="rsm_ma, rsm budget: rsm vertical mean advection", &
             var_units="(kg/kg)/s", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rsm_sd')
        stats_metadata%irsm_sd = k
        call stat_assign( var_index=stats_metadata%irsm_sd, var_name="rsm_sd", &
             var_description="rsm_sd, rsm budget: rsm sedimentation", &
             var_units="(kg/kg)/s", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rsm_sd_morr')
        stats_metadata%irsm_sd_morr = k
        call stat_assign( var_index=stats_metadata%irsm_sd_morr, var_name="rsm_sd_morr", &
             var_description="rsm_sd_morr, rsm sedimentation when using morrison microphysics &
             &(Not in budget, included in rsm_mc)", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('rsm_ta')
        stats_metadata%irsm_ta = k

        call stat_assign( var_index=stats_metadata%irsm_ta, var_name="rsm_ta", &
             var_description="rsm_ta, rsm budget: rsm turbulent advection", &
             var_units="(kg/kg)/s", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rsm_mc')
        stats_metadata%irsm_mc = k

        call stat_assign( var_index=stats_metadata%irsm_mc, var_name="rsm_mc", &
             var_description="rsm_mc, rsm budget: Change in rsm due to microphysics", &
             var_units="(kg/kg)/s", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rsm_hf')
        stats_metadata%irsm_hf = k

        call stat_assign( var_index=stats_metadata%irsm_hf, var_name="rsm_hf", &
             var_description="rsm_hf, rsm budget: rsm hole-filling term", &
             var_units="(kg/kg)/s", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rsm_wvhf')
        stats_metadata%irsm_wvhf = k

        call stat_assign( var_index=stats_metadata%irsm_wvhf, var_name="rsm_wvhf", &
             var_description="rsm_wvhf, rsm budget: rsm water vapor hole-filling term", &
             var_units="(kg/kg)/s", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rsm_cl')
        stats_metadata%irsm_cl = k

        call stat_assign( var_index=stats_metadata%irsm_cl, var_name="rsm_cl", &
             var_description="rsm_cl, rsm budget: rsm clipping term", &
             var_units="(kg/kg)/s", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('Nsm_bt')
        stats_metadata%iNsm_bt = k
        call stat_assign( var_index=stats_metadata%iNsm_bt, var_name="Nsm_bt", &
             var_description="Nsm_bt, Nsm budget", var_units="(num/kg)/s", &
             l_silhs=.false., grid_kind=stats_zt )

        k = k + 1

      case ('Nsm_ma')
        stats_metadata%iNsm_ma = k

        call stat_assign( var_index=stats_metadata%iNsm_ma, var_name="Nsm_ma", &
             var_description="Nsm_ma, Nsm budget: Nsm mean advection", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('Nsm_sd')
        stats_metadata%iNsm_sd = k

        call stat_assign( var_index=stats_metadata%iNsm_sd, var_name="Nsm_sd", &
             var_description="Nsm_sd, Nsm budget: Nsm sedimentation", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=stats_zt )

        k = k + 1

      case ('Nsm_ta')
        stats_metadata%iNsm_ta = k
        call stat_assign( var_index=stats_metadata%iNsm_ta, var_name="Nsm_ta", &
             var_description="Nsm_ta, Nsm budget: Nsm turbulent advection", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=stats_zt )

        k = k + 1

      case ('Nsm_mc')
        stats_metadata%iNsm_mc = k
        call stat_assign( var_index=stats_metadata%iNsm_mc, var_name="Nsm_mc", &
             var_description="Nsm_mc, Nsm budget: Nsm microphysics", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=stats_zt )

        k = k + 1

      case ('Nsm_cl')
        stats_metadata%iNsm_cl = k

        call stat_assign( var_index=stats_metadata%iNsm_cl, var_name="Nsm_cl", &
             var_description="Nsm_cl, Nsm budget: Nsm clipping term", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rim_bt')
        stats_metadata%irim_bt = k

        call stat_assign( var_index=stats_metadata%irim_bt, var_name="rim_bt", &
             var_description="rim_bt, rim budget: rim time tendency", &
             var_units="(kg/kg)/s", l_silhs=.false., grid_kind=stats_zt )

        k = k + 1

      case ('rim_ma')
        stats_metadata%irim_ma = k

        call stat_assign( var_index=stats_metadata%irim_ma, var_name="rim_ma", &
             var_description="rim_ma, rim budget: rim vertical mean advection", &
             var_units="(kg/kg)/s", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rim_sd')
        stats_metadata%irim_sd = k

        call stat_assign( var_index=stats_metadata%irim_sd, var_name="rim_sd", &
             var_description="rim_sd, rim budget: rim sedimentation", &
             var_units="(kg/kg)/s", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rim_sd_mg_morr')
        stats_metadata%irim_sd_mg_morr = k

        call stat_assign( var_index=stats_metadata%irim_sd_mg_morr, var_name="rim_sd_mg_morr", &
             var_description="rim_sd_mg_morr, rim sedimentation when using morrison or MG " &
             // "microphysics" &
             // "(not in budget, included in rim_mc)", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('rim_ta')
        stats_metadata%irim_ta = k

        call stat_assign( var_index=stats_metadata%irim_ta, var_name="rim_ta", &
             var_description="rim_ta, rim budget: rim turbulent advection", &
             var_units="(kg/kg)/s", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rim_mc')
        stats_metadata%irim_mc = k

        call stat_assign( var_index=stats_metadata%irim_mc, var_name="rim_mc", &
             var_description="rim_mc, rim budget: Change in rim due to microphysics", &
             var_units="(kg/kg)/s", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rim_hf')
        stats_metadata%irim_hf = k

        call stat_assign( var_index=stats_metadata%irim_hf, var_name="rim_hf", &
             var_description="rim_hf, rim budget: rim hole-filling term", &
             var_units="(kg/kg)/s", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rim_wvhf')
        stats_metadata%irim_wvhf = k

        call stat_assign( var_index=stats_metadata%irim_wvhf, var_name="rim_wvhf", &
             var_description="rim_wvhf, rim budget: rim water vapor hole-filling term", &
             var_units="(kg/kg)/s", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rim_cl')
        stats_metadata%irim_cl = k

        call stat_assign( var_index=stats_metadata%irim_cl, var_name="rim_cl", &
             var_description="rim_cl, rim budget: rim clipping term", &
             var_units="(kg/kg)/s", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rgm_bt')
        stats_metadata%irgm_bt = k

        call stat_assign( var_index=stats_metadata%irgm_bt, var_name="rgm_bt", &
             var_description="rgm_bt, rgm budget: rgm time tendency", &
             var_units="(kg/kg)/s", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rgm_ma')
        stats_metadata%irgm_ma = k

        call stat_assign( var_index=stats_metadata%irgm_ma, var_name="rgm_ma", &
             var_description="rgm_ma, rgm budget: rgm vertical mean advection", &
             var_units="(kg/kg)/s", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rgm_sd')
        stats_metadata%irgm_sd = k

        call stat_assign( var_index=stats_metadata%irgm_sd, var_name="rgm_sd", &
             var_description="rgm_sd, rgm budget: rgm sedimentation", &
             var_units="(kg/kg)/s", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rgm_sd_morr')
        stats_metadata%irgm_sd_morr = k

        call stat_assign( var_index=stats_metadata%irgm_sd_morr, var_name="rgm_sd_morr", &
             var_description="rgm_sd_morr, rgm sedimentation when using morrison microphysics &
             &(not in budget, included in rgm_mc)", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('rgm_ta')
        stats_metadata%irgm_ta = k

        call stat_assign( var_index=stats_metadata%irgm_ta, var_name="rgm_ta", &
             var_description="rgm_ta, rgm budget: rgm turbulent advection", &
             var_units="(kg/kg)/s", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rgm_mc')
        stats_metadata%irgm_mc = k

        call stat_assign( var_index=stats_metadata%irgm_mc, var_name="rgm_mc", &
             var_description="rgm_mc, rgm budget: Change in rgm due to microphysics", &
             var_units="(kg/kg)/s", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rgm_hf')
        stats_metadata%irgm_hf = k

        call stat_assign( var_index=stats_metadata%irgm_hf, var_name="rgm_hf", &
             var_description="rgm_hf, rgm budget: rgm hole-filling term", &
             var_units="(kg/kg)/s", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rgm_wvhf')
        stats_metadata%irgm_wvhf = k

        call stat_assign( var_index=stats_metadata%irgm_wvhf, var_name="rgm_wvhf", &
             var_description="rgm_wvhf, rgm budget: rgm water vapor hole-filling term", &
             var_units="(kg/kg)/s", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rgm_cl')
        stats_metadata%irgm_cl = k

        call stat_assign( var_index=stats_metadata%irgm_cl, var_name="rgm_cl", &
             var_description="rgm_cl, rgm budget: rgm clipping term", &
             var_units="(kg/kg)/s", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('Ngm_bt')
        stats_metadata%iNgm_bt = k
        call stat_assign( var_index=stats_metadata%iNgm_bt, var_name="Ngm_bt", &
             var_description="Ngm_bt, Ngm budget:", var_units="(num/kg)/s", &
             l_silhs=.false., grid_kind=stats_zt )

        k = k + 1

      case ('Ngm_ma')
        stats_metadata%iNgm_ma = k

        call stat_assign( var_index=stats_metadata%iNgm_ma, var_name="Ngm_ma", &
             var_description="Ngm_ma, Ngm budget: Ngm mean advection", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('Ngm_sd')
        stats_metadata%iNgm_sd = k

        call stat_assign( var_index=stats_metadata%iNgm_sd, var_name="Ngm_sd", &
             var_description="Ngm_sd, Ngm budget: Ngm sedimentation", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=stats_zt )

        k = k + 1

      case ('Ngm_ta')
        stats_metadata%iNgm_ta = k
        call stat_assign( var_index=stats_metadata%iNgm_ta, var_name="Ngm_ta", &
             var_description="Ngm_ta, Ngm budget: Ngm turbulent advection", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=stats_zt )

        k = k + 1

      case ('Ngm_mc')
        stats_metadata%iNgm_mc = k

        call stat_assign( var_index=stats_metadata%iNgm_mc, var_name="Ngm_mc", &
             var_description="Ngm_mc, Ngm budget: Ngm microphysics term", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('Ngm_cl')
        stats_metadata%iNgm_cl = k

        call stat_assign( var_index=stats_metadata%iNgm_cl, var_name="Ngm_cl", &
             var_description="Ngm_cl, Ngm budget: Ngm clipping term", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('Nim_bt')
        stats_metadata%iNim_bt = k
        call stat_assign( var_index=stats_metadata%iNim_bt, var_name="Nim_bt", &
             var_description="Nim_bt, Nim budget", var_units="(num/kg)/s", l_silhs=.false., &
             grid_kind=stats_zt )

        k = k + 1

      case ('Nim_ma')
        stats_metadata%iNim_ma = k

        call stat_assign( var_index=stats_metadata%iNim_ma, var_name="Nim_ma", &
             var_description="Nim_ma, Nim budget: Nim mean advection", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('Nim_sd')
        stats_metadata%iNim_sd = k

        call stat_assign( var_index=stats_metadata%iNim_sd, var_name="Nim_sd", &
             var_description="Nim_sd, Nim budget: Nim sedimentation", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=stats_zt )

        k = k + 1

      case ('Nim_ta')
        stats_metadata%iNim_ta = k
        call stat_assign( var_index=stats_metadata%iNim_ta, var_name="Nim_ta", &
             var_description="Nim_ta, Nim budget: Nim turbulent advection", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=stats_zt )

        k = k + 1

      case ('Nim_mc')
        stats_metadata%iNim_mc = k

        call stat_assign( var_index=stats_metadata%iNim_mc, var_name="Nim_mc", &
             var_description="Nim_mc, Nim budget: Nim microphysics term", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('Nim_cl')
        stats_metadata%iNim_cl = k

        call stat_assign( var_index=stats_metadata%iNim_cl, var_name="Nim_cl", &
             var_description="Nim_cl, Nim budget: Nim clipping term", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('Ncm_bt')
        stats_metadata%iNcm_bt = k
        call stat_assign( var_index=stats_metadata%iNcm_bt, var_name="Ncm_bt", &
             var_description="Ncm_bt, Ncm budget: Cloud droplet number concentration budget", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=stats_zt )

        k = k + 1

      case ('Ncm_ma')
        stats_metadata%iNcm_ma = k

        call stat_assign( var_index=stats_metadata%iNcm_ma, var_name="Ncm_ma", &
             var_description="Ncm_ma, Ncm budget: Ncm vertical mean advection", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('Ncm_act')
        stats_metadata%iNcm_act = k

        call stat_assign( var_index=stats_metadata%iNcm_act, var_name="Ncm_act", &
             var_description="Ncm_act, Ncm budget: Change in Ncm due to activation", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=stats_zt )

        k = k + 1

      case ('Ncm_ta')
        stats_metadata%iNcm_ta = k
        call stat_assign( var_index=stats_metadata%iNcm_ta, var_name="Ncm_ta", &
             var_description="Ncm_ta, Ncm budget: Ncm turbulent advection", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=stats_zt )

        k = k + 1

      case ('Ncm_mc')
        stats_metadata%iNcm_mc = k

        call stat_assign( var_index=stats_metadata%iNcm_mc, var_name="Ncm_mc", &
             var_description="Ncm_mc, Ncm budget: Change in Ncm due to microphysics", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('Ncm_cl')
        stats_metadata%iNcm_cl = k

        call stat_assign( var_index=stats_metadata%iNcm_cl, var_name="Ncm_cl", &
             var_description="Ncm_cl, Ncm budget: Ncm clipping term", &
             var_units="(num/kg)/s", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('PSMLT')
        stats_metadata%iPSMLT = k

        call stat_assign( var_index=stats_metadata%iPSMLT, var_name="PSMLT", &
             var_description="PSMLT, Freezing of rain to form snow, +rsm, -rrm", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('EVPMS')
        stats_metadata%iEVPMS = k

        call stat_assign( var_index=stats_metadata%iEVPMS, var_name="EVPMS", &
             var_description="EVPMS, Evaporation of melted snow, +rsm, -rvm", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('PRACS')
        stats_metadata%iPRACS = k

        call stat_assign( var_index=stats_metadata%iPRACS, var_name="PRACS", &
             var_description="PRACS, Collection of rain by snow, +rsm, -rrm", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('EVPMG')
        stats_metadata%iEVPMG = k

        call stat_assign( var_index=stats_metadata%iEVPMG, var_name="EVPMG", &
             var_description="EVPMG, Evaporation of melted graupel, +rgm, -rvm", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('PRACG')
        stats_metadata%iPRACG = k

        call stat_assign( var_index=stats_metadata%iPRACG, var_name="PRACG", &
             var_description="PRACG, Negative of collection of rain by graupel, +rrm, -rgm", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('PGMLT')
        stats_metadata%iPGMLT = k

        call stat_assign( var_index=stats_metadata%iPGMLT, var_name="PGMLT", &
             var_description="PGMLT, Negative of melting of graupel, +rgm, -rrm", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('MNUCCC')
        stats_metadata%iMNUCCC = k

        call stat_assign( var_index=stats_metadata%iMNUCCC, var_name="MNUCCC", &
             var_description="MNUCCC, Contact freezing of cloud droplets, +rim, -rcm", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('PSACWS')
        stats_metadata%iPSACWS = k

        call stat_assign( var_index=stats_metadata%iPSACWS, var_name="PSACWS", &
             var_description="PSACWS, Collection of cloud water by snow, +rsm, -rcm", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('PSACWI')
        stats_metadata%iPSACWI = k

        call stat_assign( var_index=stats_metadata%iPSACWI, var_name="PSACWI", &
             var_description="PSACWI, Collection of cloud water by cloud ice, +rim, -rcm", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('QMULTS')
        stats_metadata%iQMULTS = k

        call stat_assign( var_index=stats_metadata%iQMULTS, var_name="QMULTS", &
             var_description="QMULTS, Splintering from cloud droplets accreted onto snow, " &
             // "+rim, -rcm", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('QMULTG')
        stats_metadata%iQMULTG = k

        call stat_assign( var_index=stats_metadata%iQMULTG, var_name="QMULTG", &
             var_description="QMULTG, Splintering from droplets accreted onto graupel, " &
             // "+rim, -rcm", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('PSACWG')
        stats_metadata%iPSACWG = k

        call stat_assign( var_index=stats_metadata%iPSACWG, var_name="PSACWG", &
             var_description="PSACWG, Collection of cloud water by graupel, +rgm, -rcm", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('PGSACW')
        stats_metadata%iPGSACW = k

        call stat_assign( var_index=stats_metadata%iPGSACW, var_name="PGSACW", &
             var_description="PGSACW, Reclassification of rimed snow as graupel, +rgm, -rcm", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('PRD')
        stats_metadata%iPRD = k

        call stat_assign( var_index=stats_metadata%iPRD, var_name="PRD", &
             var_description="PRD, Depositional growth of cloud ice, +rim, -rvm", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('PRCI')
        stats_metadata%iPRCI = k

        call stat_assign( var_index=stats_metadata%iPRCI, var_name="PRCI", &
             var_description="PRCI, Autoconversion of cloud ice to snow, +rsm, -rim", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('PRAI')
        stats_metadata%iPRAI = k

        call stat_assign( var_index=stats_metadata%iPRAI, var_name="PRAI", &
             var_description="PRAI, Collection of cloud ice by snow, +rsm, -rim", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('QMULTR')
        stats_metadata%iQMULTR = k

        call stat_assign( var_index=stats_metadata%iQMULTR, var_name="QMULTR", &
             var_description="QMULTR, Splintering from rain droplets accreted onto snow, " &
             // "+rim, -rrm", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('QMULTRG')
        stats_metadata%iQMULTRG = k

        call stat_assign( var_index=stats_metadata%iQMULTRG, var_name="QMULTRG", &
             var_description="QMULTRG, Splintering from rain droplets accreted onto graupel, " &
             // "+rim, -rrm", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('MNUCCD')
        stats_metadata%iMNUCCD = k

        call stat_assign( var_index=stats_metadata%iMNUCCD, var_name="MNUCCD", &
             var_description="MNUCCD, Freezing of aerosol, +rim, -rvm", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('PRACI')
        stats_metadata%iPRACI = k

        call stat_assign( var_index=stats_metadata%iPRACI, var_name="PRACI", &
             var_description="PRACI, Collection of cloud ice by rain, +rgm, -rim", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('PRACIS')
        stats_metadata%iPRACIS = k

        call stat_assign( var_index=stats_metadata%iPRACIS, var_name="PRACIS", &
             var_description="PRACIS, Collection of cloud ice by rain, +rsm, -rim", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('EPRD')
        stats_metadata%iEPRD = k

        call stat_assign( var_index=stats_metadata%iEPRD, var_name="EPRD", &
             var_description="EPRD, Negative of sublimation of cloud ice, +rim, -rvm", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('MNUCCR')
        stats_metadata%iMNUCCR = k

        call stat_assign( var_index=stats_metadata%iMNUCCR, var_name="MNUCCR", &
             var_description="MNUCCR, Contact freezing of rain droplets, +rgm, -rrm", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('PIACR')
        stats_metadata%iPIACR = k

        call stat_assign( var_index=stats_metadata%iPIACR, var_name="PIACR", &
             var_description="PIACR, Collection of cloud ice by rain, +rgm, -rrm", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('PIACRS')
        stats_metadata%iPIACRS = k

        call stat_assign( var_index=stats_metadata%iPIACRS, var_name="PIACRS", &
             var_description="PIACRS, Collection of cloud ice by rain, +rsm, -rrm", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('PGRACS')
        stats_metadata%iPGRACS = k

        call stat_assign( var_index=stats_metadata%iPGRACS, var_name="PGRACS", &
             var_description="PGRACS, Collection of rain by snow, +rgm, -rrm", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('PRDS')
        stats_metadata%iPRDS = k

        call stat_assign( var_index=stats_metadata%iPRDS, var_name="PRDS", &
             var_description="PRDS, Depositional growth of snow, +rsm, -rvm", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('EPRDS')
        stats_metadata%iEPRDS = k

        call stat_assign( var_index=stats_metadata%iEPRDS, var_name="EPRDS", &
             var_description="EPRDS, Negative of sublimation of snow, +rsm, -rvm", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('PSACR')
        stats_metadata%iPSACR = k

        call stat_assign( var_index=stats_metadata%iPSACR, var_name="PSACR", &
             var_description="PSACR, Collection of snow by rain, +rgm, -rsm", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('PRDG')
        stats_metadata%iPRDG = k

        call stat_assign( var_index=stats_metadata%iPRDG, var_name="PRDG", &
             var_description="PRDG, Depositional growth of graupel, +rgm, -rvm", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('EPRDG')
        stats_metadata%iEPRDG = k

        call stat_assign( var_index=stats_metadata%iEPRDG, var_name="EPRDG", &
             var_description="EPRDG, Negative of sublimation of graupel, +rgm, -rvm", &
             var_units="(kg/kg)/s", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('NGSTEN')
        stats_metadata%iNGSTEN = k

        call stat_assign( var_index=stats_metadata%iNGSTEN, var_name="NGSTEN", &
             var_description="NGSTEN, Graupel sedimentation tendency", var_units="(#/kg/s)", &
             l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('NRSTEN')
        stats_metadata%iNRSTEN = k

        call stat_assign( var_index=stats_metadata%iNRSTEN, var_name="NRSTEN", &
             var_description="NRSTEN, Rain sedimentation tendency", var_units="(#/kg/s)", &
             l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('NISTEN')
        stats_metadata%iNISTEN = k

        call stat_assign( var_index=stats_metadata%iNISTEN, var_name="NISTEN", &
             var_description="NISTEN, Cloud ice sedimentation tendency", var_units="(#/kg/s)", &
             l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('NSSTEN')
        stats_metadata%iNSSTEN = k

        call stat_assign( var_index=stats_metadata%iNSSTEN, var_name="NSSTEN", &
             var_description="NSSTEN, Snow sedimentation tendency", var_units="(#/kg/s)", &
             l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('NCSTEN')
        stats_metadata%iNCSTEN = k

        call stat_assign( var_index=stats_metadata%iNCSTEN, var_name="NCSTEN", &
             var_description="NCSTEN, Cloud water sedimentation tendency", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('NPRC1')
        stats_metadata%iNPRC1 = k

        call stat_assign( var_index=stats_metadata%iNPRC1, var_name="NPRC1", &
             var_description="NPRC1, Change in Nrm due to autoconversion of droplets, +Nrm", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('NRAGG')
        stats_metadata%iNRAGG = k

        call stat_assign( var_index=stats_metadata%iNRAGG, var_name="NRAGG", &
             var_description="NRAGG, Change in Nrm due to self-collection of raindrops, +Nrm", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('NPRACG')
        stats_metadata%iNPRACG = k

        call stat_assign( var_index=stats_metadata%iNPRACG, var_name="NPRACG", &
             var_description="NPRACG, Collection of rainwater by graupel, -Nrm", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('NSUBR')
        stats_metadata%iNSUBR = k

        call stat_assign( var_index=stats_metadata%iNSUBR, var_name="NSUBR", &
             var_description="NSUBR, Loss of Nrm by evaporation, +Nrm", var_units="(#/kg/s)", &
             l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('NSMLTR')
        stats_metadata%iNSMLTR = k

        call stat_assign( var_index=stats_metadata%iNSMLTR, var_name="NSMLTR", &
             var_description="NSMLTR, Melting of snow to form rain, -Nrm", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('NGMLTR')
        stats_metadata%iNGMLTR = k

        call stat_assign( var_index=stats_metadata%iNGMLTR, var_name="NGMLTR", &
             var_description="NGMLTR, Melting of graupel to form rain, -Nrm", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('NPRACS')
        stats_metadata%iNPRACS = k

        call stat_assign( var_index=stats_metadata%iNPRACS, var_name="NPRACS", &
             var_description="NPRACS, Collection of rainwater by snow, -Nrm", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('NNUCCR')
        stats_metadata%iNNUCCR = k

        call stat_assign( var_index=stats_metadata%iNNUCCR, var_name="NNUCCR", &
             var_description="NNUCCR, Contact freezing of rain, +Ngm, -Nrm", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('NIACR')
        stats_metadata%iNIACR = k

        call stat_assign( var_index=stats_metadata%iNIACR, var_name="NIACR", &
             var_description="NIACR, Collection of cloud ice by rain, +Ngm, -Nrm, -Nim", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('NIACRS')
        stats_metadata%iNIACRS = k

        call stat_assign( var_index=stats_metadata%iNIACRS, var_name="NIACRS", &
             var_description="NIACRS, Collection of cloud ice by rain, +Nsm, -Nrm, -Nim", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('NGRACS')
        stats_metadata%iNGRACS = k

        call stat_assign( var_index=stats_metadata%iNGRACS, var_name="NGRACS", &
             var_description="NGRACS, Collection of rain by snow, +Ngm, -Nrm, -Nsm", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('NSMLTS')
        stats_metadata%iNSMLTS= k

        call stat_assign( var_index=stats_metadata%iNSMLTS, var_name="NSMLTS", &
             var_description="NSMLTS, Melting  of snow, +Nsm", var_units="(#/kg/s)", &
             l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('NSAGG')
        stats_metadata%iNSAGG= k

        call stat_assign( var_index=stats_metadata%iNSAGG, var_name="NSAGG", &
             var_description="NSAGG, Self collection of snow, +Nsm", var_units="(#/kg/s)", &
             l_silhs=.true., grid_kind=stats_zt )

        k = k + 1

      case ('NPRCI')
        stats_metadata%iNPRCI= k

        call stat_assign( var_index=stats_metadata%iNPRCI, var_name="NPRCI", &
             var_description="NPRCI, Autoconversion of cloud ice to snow, -Nim, +Nsm", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('NSCNG')
        stats_metadata%iNSCNG= k

        call stat_assign( var_index=stats_metadata%iNSCNG, var_name="NSCNG", &
             var_description="NSCNG, Conversion of snow to graupel, +Ngm, -Nsm", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('NSUBS')
        stats_metadata%iNSUBS= k

        call stat_assign( var_index=stats_metadata%iNSUBS, var_name="NSUBS", &
             var_description="NSUBS, Loss of snow due to sublimation, +Nsm", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('PRC')
        stats_metadata%iPRC= k

        call stat_assign( var_index=stats_metadata%iPRC, var_name="PRC", &
             var_description="PRC, Autoconversion +rrm -rcm", var_units="(kg/kg/s)", &
             l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('PRA')
        stats_metadata%iPRA= k

        call stat_assign( var_index=stats_metadata%iPRA, var_name="PRA", &
             var_description="PRA, Accretion +rrm -rcm", var_units="(kg/kg/s)", &
             l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('PRE')
        stats_metadata%iPRE= k

        call stat_assign( var_index=stats_metadata%iPRE, var_name="PRE", &
             var_description="PRE, Evaporation of rain -rrm", var_units="(kg/kg/s)", &
             l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('PCC')
        stats_metadata%iPCC= k

        call stat_assign( var_index=stats_metadata%iPCC, var_name="PCC", &
             var_description="PCC, Satuation adjustment -rvm +rcm", var_units="(kg/kg/s)", &
             l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('NNUCCC')
        stats_metadata%iNNUCCC= k

        call stat_assign( var_index=stats_metadata%iNNUCCC, var_name="NNUCCC", &
             var_description="NNUCCC, Contact freezing of drops, -Ncm + Nim", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('NPSACWS')
        stats_metadata%iNPSACWS= k

        call stat_assign( var_index=stats_metadata%iNPSACWS, var_name="NPSACWS", &
             var_description="NPSACWS, Droplet accretion by snow, -Ncm", var_units="(#/kg/s)", &
             l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('NPRA')
        stats_metadata%iNPRA= k

        call stat_assign( var_index=stats_metadata%iNPRA, var_name="NPRA", &
             var_description="NPRA, Droplet accretion by rain, -Ncm", var_units="(#/kg/s)", &
             l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('NPRC')
        stats_metadata%iNPRC= k

        call stat_assign( var_index=stats_metadata%iNPRC, var_name="NPRC", &
             var_description="NPRC, Autoconversion of cloud drops, -Ncm", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('NPSACWI')
        stats_metadata%iNPSACWI= k

        call stat_assign( var_index=stats_metadata%iNPSACWI, var_name="NPSACWI", &
             var_description="NPSACWI, Droplet accretion by cloud ice, -Ncm", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('NPSACWG')
        stats_metadata%iNPSACWG= k

        call stat_assign( var_index=stats_metadata%iNPSACWG, var_name="NPSACWG", &
             var_description="NPSACWG, Collection of cloud droplets by graupel, -Ncm", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('NPRAI')
        stats_metadata%iNPRAI= k

        call stat_assign( var_index=stats_metadata%iNPRAI, var_name="NPRAI", &
             var_description="NPRAI, Accretion of cloud ice by snow, -Nim", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('NMULTS')
        stats_metadata%iNMULTS= k

        call stat_assign( var_index=stats_metadata%iNMULTS, var_name="NMULTS", &
             var_description="NMULTS, Ice multiplication due to riming of cloud droplets " &
             // "by snow, +Nim", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('NMULTG')
        stats_metadata%iNMULTG= k

        call stat_assign( var_index=stats_metadata%iNMULTG, var_name="NMULTG", &
             var_description="NMULTG, Ice multiplication due to accretion of droplets " &
             // "by graupel, +Nim", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('NMULTR')
        stats_metadata%iNMULTR= k

        call stat_assign( var_index=stats_metadata%iNMULTR, var_name="NMULTR", &
             var_description="Ice multiplication due to riming of rain by snow, +Nim", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('NMULTRG')
        stats_metadata%iNMULTRG= k

        call stat_assign( var_index=stats_metadata%iNMULTRG, var_name="NMULTRG", &
             var_description="NMULTR, Ice multiplication due to accretion of rain by " &
             // "graupel, +Nim", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('NNUCCD')
        stats_metadata%iNNUCCD= k

        call stat_assign( var_index=stats_metadata%iNNUCCD, var_name="NNUCCD", &
             var_description="NNUCCD, Primary ice nucleation, freezing of aerosol, +Nim", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('NSUBI')
        stats_metadata%iNSUBI= k

        call stat_assign( var_index=stats_metadata%iNSUBI, var_name="NSUBI", &
             var_description="NSUBI, Loss of ice due to sublimation, -Nim", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('NGMLTG')
        stats_metadata%iNGMLTG= k

        call stat_assign( var_index=stats_metadata%iNGMLTG, var_name="NGMLTG", &
             var_description="NGMLTG, Loss of graupel due to melting, -Ngm", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('NSUBG')
        stats_metadata%iNSUBG= k

        call stat_assign( var_index=stats_metadata%iNSUBG, var_name="NSUBG", &
             var_description="NSUBG, Loss of graupel due to sublimation, -Ngm", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('NACT')
        stats_metadata%iNACT= k

        call stat_assign( var_index=stats_metadata%iNACT, var_name="NACT", &
             var_description="NACT, Cloud drop formation by aerosol activation, +Ncm", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('SIZEFIX_NR')
        stats_metadata%iSIZEFIX_NR= k

        call stat_assign( var_index=stats_metadata%iSIZEFIX_NR, var_name="SIZEFIX_NR", &
             var_description="SIZEFIX_NR, Adjust rain # conc. for large/small drops, +Nrm", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('SIZEFIX_NC')
        stats_metadata%iSIZEFIX_NC= k

        call stat_assign( var_index=stats_metadata%iSIZEFIX_NC, var_name="SIZEFIX_NC", &
             var_description="SIZEFIX_NC, Adjust cloud # conc. for large/small drops, +Ncm", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('SIZEFIX_NI')
        stats_metadata%iSIZEFIX_NI= k

        call stat_assign( var_index=stats_metadata%iSIZEFIX_NI, var_name="SIZEFIX_NI", &
             var_description="SIZEFIX_NI, Adjust ice # conc. for large/small drops, +Nim", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('SIZEFIX_NS')
        stats_metadata%iSIZEFIX_NS= k

        call stat_assign( var_index=stats_metadata%iSIZEFIX_NS, var_name="SIZEFIX_NS", &
             var_description="SIZEFIX_NS, Adjust snow # conc. for large/small drops, +Nsm", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('SIZEFIX_NG')
        stats_metadata%iSIZEFIX_NG= k

        call stat_assign( var_index=stats_metadata%iSIZEFIX_NG, var_name="SIZEFIX_NG", &
             var_description="SIZEFIX_NG, Adjust graupel # conc. for large/small drops,+Ngm",&
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('NEGFIX_NR')
        stats_metadata%iNEGFIX_NR= k

        call stat_assign( var_index=stats_metadata%iNEGFIX_NR, var_name="NEGFIX_NR", &
             var_description="NEGFIX_NR, Removal of negative rain drop number conc., -Nrm", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('NEGFIX_NC')
        stats_metadata%iNEGFIX_NC= k

        call stat_assign( var_index=stats_metadata%iNEGFIX_NC, var_name="NEGFIX_NC", &
             var_description="NEGFIX_NC, Removal of negative cloud drop number conc., -Ncm", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('NEGFIX_NI')
        stats_metadata%iNEGFIX_NI= k

        call stat_assign( var_index=stats_metadata%iNEGFIX_NI, var_name="NEGFIX_NI", &
             var_description="NEGFIX_NI, Removal of negative ice number conc., -Nim", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('NEGFIX_NS')
        stats_metadata%iNEGFIX_NS= k

        call stat_assign( var_index=stats_metadata%iNEGFIX_NS, var_name="NEGFIX_NS", &
             var_description="NEGFIX_NS, Removal of negative snow number conc, -Nsm", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('NEGFIX_NG')
        stats_metadata%iNEGFIX_NG= k

        call stat_assign( var_index=stats_metadata%iNEGFIX_NG, var_name="NEGFIX_NG", &
             var_description="NEGFIX_NG, Removal of negative graupel number conc., -Ngm", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('NIM_MORR_CL')
        stats_metadata%iNIM_MORR_CL= k

        call stat_assign( var_index=stats_metadata%iNIM_MORR_CL, var_name="NIM_MORR_CL", &
             var_description="NIM_MORR_CL, Clipping of large ice concentrations, -Nim", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('QC_INST')
        stats_metadata%iQC_INST= k

        call stat_assign( var_index=stats_metadata%iQC_INST, var_name="QC_INST", &
             var_description="QC_INST, Change in mixing ratio due to instantaneous " &
             // "processes +rcm", &
             var_units="(kg/kg/s)", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('QR_INST')
        stats_metadata%iQR_INST= k

        call stat_assign( var_index=stats_metadata%iQR_INST, var_name="QR_INST", &
             var_description="QR_INST, Change in mixing ratio from instantaneous processes, +rrm",&
             var_units="(kg/kg/s)", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('QI_INST')
        stats_metadata%iQI_INST= k

        call stat_assign( var_index=stats_metadata%iQI_INST, var_name="QI_INST", &
             var_description="QI_INST, Change in mixing ratio from instantaneous processes +rim",&
             var_units="(kg/kg/s)", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('QS_INST')
        stats_metadata%iQS_INST= k

        call stat_assign( var_index=stats_metadata%iQS_INST, var_name="QS_INST", &
             var_description="QS_INST, Change in mixing ratio from instantaneous processes +rsm",&
             var_units="(kg/kg/s)", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('QG_INST')
        stats_metadata%iQG_INST= k

        call stat_assign( var_index=stats_metadata%iQG_INST, var_name="QG_INST", &
             var_description="QG_INST, Change in mixing ratio from instantaneous processes +rgm",&
             var_units="(kg/kg/s)", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('NC_INST')
        stats_metadata%iNC_INST= k

        call stat_assign( var_index=stats_metadata%iNC_INST, var_name="NC_INST", &
             var_description="NC_INST, Change in # conc. from instantaneous processes +Ncm", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('NR_INST')
        stats_metadata%iNR_INST= k

        call stat_assign( var_index=stats_metadata%iNR_INST, var_name="NR_INST", &
             var_description="NR_INST, Change in # conc. from instantaneous processes +Nrm", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('NI_INST')
        stats_metadata%iNI_INST= k

        call stat_assign( var_index=stats_metadata%iNI_INST, var_name="NI_INST", &
             var_description="NI_INST, Change in # conc. from instantaneous processes +Nim", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('NS_INST')
        stats_metadata%iNS_INST= k

        call stat_assign( var_index=stats_metadata%iNS_INST, var_name="NS_INST", &
             var_description="NS_INST, Change in # conc. from instantaneous processes +Nsm", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('NG_INST')
        stats_metadata%iNG_INST= k

        call stat_assign( var_index=stats_metadata%iNG_INST, var_name="NG_INST", &
             var_description="NG_INST, Change in # conc. from instantaneous processes +Ngm", &
             var_units="(#/kg/s)", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1


      case ('T_in_K_mc')
        stats_metadata%iT_in_K_mc= k

        call stat_assign( var_index=stats_metadata%iT_in_K_mc, var_name="T_in_K_mc", &
             var_description="T_in_K_mc, Temperature tendency from Morrison microphysics", &
             var_units="(K/s)", l_silhs=.true., grid_kind=stats_zt )
        k = k + 1

      case ('w_KK_evap_covar_zt')
        stats_metadata%iw_KK_evap_covar_zt = k

        call stat_assign( var_index=stats_metadata%iw_KK_evap_covar_zt, var_name="w_KK_evap_covar_zt", &
             var_description="w_KK_evap_covar_zt, Covariance of w and KK evaporation rate", &
             var_units="m*(kg/kg)/s^2", l_silhs=.false., grid_kind=stats_zt )
       k = k + 1

      case ('rt_KK_evap_covar_zt')
        stats_metadata%irt_KK_evap_covar_zt = k

        call stat_assign( var_index=stats_metadata%irt_KK_evap_covar_zt, var_name="rt_KK_evap_covar_zt", &
             var_description="rt_KK_evap_covar_zt, Covariance of r_t and KK evaporation rate", &
             var_units="(kg/kg)^2/s", l_silhs=.false., grid_kind=stats_zt )
       k = k + 1

      case ('thl_KK_evap_covar_zt')
        stats_metadata%ithl_KK_evap_covar_zt = k

        call stat_assign( var_index=stats_metadata%ithl_KK_evap_covar_zt, var_name="thl_KK_evap_covar_zt", &
             var_description="thl_KK_evap_covar_zt, Covariance of theta_l and KK " &
             // "evaporation rate", &
             var_units="K*(kg/kg)/s", l_silhs=.false., grid_kind=stats_zt )
       k = k + 1

      case ('w_KK_auto_covar_zt')
        stats_metadata%iw_KK_auto_covar_zt = k

        call stat_assign( var_index=stats_metadata%iw_KK_auto_covar_zt, var_name="w_KK_auto_covar_zt", &
             var_description="w_KK_auto_covar_zt, Covariance of w and KK autoconversion rate", &
             var_units="m*(kg/kg)/s^2", l_silhs=.false., grid_kind=stats_zt )
       k = k + 1

      case ('rt_KK_auto_covar_zt')
        stats_metadata%irt_KK_auto_covar_zt = k

        call stat_assign( var_index=stats_metadata%irt_KK_auto_covar_zt, var_name="rt_KK_auto_covar_zt", &
             var_description="rt_KK_auto_covar_zt, Covariance of r_t and KK autoconversion rate",&
             var_units="(kg/kg)^2/s", l_silhs=.false., grid_kind=stats_zt )
       k = k + 1

      case ('thl_KK_auto_covar_zt')
        stats_metadata%ithl_KK_auto_covar_zt = k

        call stat_assign( var_index=stats_metadata%ithl_KK_auto_covar_zt, var_name="thl_KK_auto_covar_zt", &
             var_description="thl_KK_auto_covar_zt, Covariance of theta_l and " &
             // "KK autoconversion rate", &
             var_units="K*(kg/kg)/s", l_silhs=.false., grid_kind=stats_zt )
       k = k + 1

      case ('w_KK_accr_covar_zt')
        stats_metadata%iw_KK_accr_covar_zt = k

        call stat_assign( var_index=stats_metadata%iw_KK_accr_covar_zt, var_name="w_KK_accr_covar_zt", &
             var_description="w_KK_accr_covar_zt, Covariance of w and KK accretion rate", &
             var_units="m*(kg/kg)/s^2", &
             l_silhs=.false., grid_kind=stats_zt )
       k = k + 1

      case ('rt_KK_accr_covar_zt')
        stats_metadata%irt_KK_accr_covar_zt = k

        call stat_assign( var_index=stats_metadata%irt_KK_accr_covar_zt, var_name="rt_KK_accr_covar_zt", &
             var_description="rt_KK_accr_covar_zt, Covariance of r_t and KK accretion rate", &
             var_units="(kg/kg)^2/s", &
             l_silhs=.false., grid_kind=stats_zt )
       k = k + 1

      case ('thl_KK_accr_covar_zt')
        stats_metadata%ithl_KK_accr_covar_zt = k

        call stat_assign( var_index=stats_metadata%ithl_KK_accr_covar_zt, var_name="thl_KK_accr_covar_zt", &
             var_description="thl_KK_accr_covar_zt, Covariance of theta_l and KK accretion rate",&
             var_units="K*(kg/kg)/s", l_silhs=.false., grid_kind=stats_zt )
       k = k + 1

      case ('rr_KK_mvr_covar_zt')
        stats_metadata%irr_KK_mvr_covar_zt = k

        call stat_assign( var_index=stats_metadata%irr_KK_mvr_covar_zt, var_name="rr_KK_mvr_covar_zt", &
             var_description="rr_KK_mvr_covar_zt, Covariance of r_r and KK rain drop mean " &
             // "volume radius", &
             var_units="(kg/kg)m", l_silhs=.false., grid_kind=stats_zt )
       k = k + 1

      case ('Nr_KK_mvr_covar_zt')
        stats_metadata%iNr_KK_mvr_covar_zt = k

        call stat_assign( var_index=stats_metadata%iNr_KK_mvr_covar_zt, var_name="Nr_KK_mvr_covar_zt", &
             var_description="Nr_KK_mvr_covar_zt, Covariance of N_r and KK rain drop mean " &
             // "volume radius", &
             var_units="(num/kg)m", l_silhs=.false., grid_kind=stats_zt )
       k = k + 1

      case ('KK_mvr_variance_zt')
        stats_metadata%iKK_mvr_variance_zt = k

        call stat_assign( var_index=stats_metadata%iKK_mvr_variance_zt, var_name="KK_mvr_variance_zt", &
             var_description="KK_mvr_variance_zt, Variance of KK rain drop mean volume radius", &
             var_units="m^2", l_silhs=.false., grid_kind=stats_zt )
       k = k + 1

      case ('vm_bt')
        stats_metadata%ivm_bt = k

        call stat_assign( var_index=stats_metadata%ivm_bt, var_name="vm_bt", &
             var_description="vm_bt, vm budget: vm time tendency", var_units="m s^{-2}", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('vm_ma')
        stats_metadata%ivm_ma = k
        call stat_assign( var_index=stats_metadata%ivm_ma, var_name="vm_ma", &
             var_description="vm_ma, vm budget: vm vertical mean advection", &
             var_units="m s^{-2}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('vm_gf')
        stats_metadata%ivm_gf = k

        call stat_assign( var_index=stats_metadata%ivm_gf, var_name="vm_gf", &
             var_description="vm_gf, vm budget: vm geostrophic forcing", &
             var_units="m s^{-2}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('vm_cf')
        stats_metadata%ivm_cf = k

        call stat_assign( var_index=stats_metadata%ivm_cf, var_name="vm_cf", &
             var_description="vm_cf, vm budget: vm coriolis forcing", var_units="m s^{-2}", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('vm_ta')
        stats_metadata%ivm_ta = k

        call stat_assign( var_index=stats_metadata%ivm_ta, var_name="vm_ta", &
             var_description="vm_ta, vm budget: vm turbulent transport", &
             var_units="m s^{-2}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('vm_f')
        stats_metadata%ivm_f = k
        call stat_assign( var_index=stats_metadata%ivm_f, var_name="vm_f", &
             var_description="vm_f, vm budget: vm forcing", var_units="m s^{-2}", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('vm_sdmp')
        stats_metadata%ivm_sdmp = k
        call stat_assign( var_index=stats_metadata%ivm_sdmp, var_name="vm_sdmp", &
             var_description="vm_sdmp, vm budget: vm sponge damping", var_units="m s^{-2}", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('vm_ndg')
        stats_metadata%ivm_ndg = k
        call stat_assign( var_index=stats_metadata%ivm_ndg, var_name="vm_ndg", &
             var_description="vm_ndg, vm budget: vm nudging", var_units="m s^{-2}", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('vm_mfl')
        stats_metadata%ivm_mfl = k
        call stat_assign( var_index=stats_metadata%ivm_mfl, var_name="vm_mfl", &
             var_description="vm_mfl, vm budget: vm monotonic flux limiter", &
             var_units="m s^{-2}", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('um_bt')
        stats_metadata%ium_bt = k

        call stat_assign( var_index=stats_metadata%ium_bt, var_name="um_bt", &
             var_description="um_bt, um budget: um time tendency", var_units="m s^{-2}", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('um_ma')
        stats_metadata%ium_ma = k

        call stat_assign( var_index=stats_metadata%ium_ma, var_name="um_ma", &
             var_description="um_ma, um budget: um vertical mean advection", &
             var_units="m s^{-2}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('um_gf')
        stats_metadata%ium_gf = k
        call stat_assign( var_index=stats_metadata%ium_gf, var_name="um_gf", &
             var_description="um_gf, um budget: um geostrophic forcing", &
             var_units="m s^{-2}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('um_cf')
        stats_metadata%ium_cf = k
        call stat_assign( var_index=stats_metadata%ium_cf, var_name="um_cf", &
             var_description="um_cf, um budget: um coriolis forcing", var_units="m s^{-2}", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('um_ta')
        stats_metadata%ium_ta = k
        call stat_assign( var_index=stats_metadata%ium_ta, var_name="um_ta", &
             var_description="um_ta, um budget: um turbulent advection", &
             var_units="m s^{-2}", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('um_f')
        stats_metadata%ium_f = k
        call stat_assign( var_index=stats_metadata%ium_f, var_name="um_f", &
             var_description="um_f, um budget: um forcing", var_units="m s^{-2}", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('um_sdmp')
        stats_metadata%ium_sdmp = k
        call stat_assign( var_index=stats_metadata%ium_sdmp, var_name="um_sdmp", &
             var_description="um_sdmp, um budget: um sponge damping", var_units="m s^{-2}", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('um_ndg')
        stats_metadata%ium_ndg = k
        call stat_assign( var_index=stats_metadata%ium_ndg, var_name="um_ndg", &
             var_description="um_ndg, um budget: um nudging", var_units="m s^{-2}", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('um_mfl')
        stats_metadata%ium_mfl = k
        call stat_assign( var_index=stats_metadata%ium_mfl, var_name="um_mfl", &
             var_description="um_mfl, um budget: um monotonic flux limiter", &
             var_units="m s^{-2}",&
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('mixt_frac')
        stats_metadata%imixt_frac = k
        call stat_assign( var_index=stats_metadata%imixt_frac, var_name="mixt_frac", &
             var_description="mixt_frac, pdf parameter: mixture fraction", var_units="count", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('w_1')
        stats_metadata%iw_1 = k
        call stat_assign( var_index=stats_metadata%iw_1, var_name="w_1", &
             var_description="w_1, pdf parameter: mean w of component 1", var_units="m/s", &
             l_silhs=.false., grid_kind=stats_zt )

        k = k + 1

      case ('w_2')
        stats_metadata%iw_2 = k

        call stat_assign( var_index=stats_metadata%iw_2, var_name="w_2", &
             var_description="w_2, pdf paramete: mean w of component 2", var_units="m/s", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('varnce_w_1')
        stats_metadata%ivarnce_w_1 = k
        call stat_assign( var_index=stats_metadata%ivarnce_w_1, var_name="varnce_w_1", &
             var_description="varnce_w_1, pdf parameter: w variance of component 1", &
             var_units="m^2/s^2", l_silhs=.false., grid_kind=stats_zt )

        k = k + 1

      case ('varnce_w_2')
        stats_metadata%ivarnce_w_2 = k

        call stat_assign( var_index=stats_metadata%ivarnce_w_2, var_name="varnce_w_2", &
             var_description="varnce_w_2, pdf parameter: w variance of component 2", &
             var_units="m^2/s^2", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('thl_1')
        stats_metadata%ithl_1 = k

        call stat_assign( var_index=stats_metadata%ithl_1, var_name="thl_1", &
             var_description="thl_1, pdf parameter: mean thl of component 1", var_units="K", &
             l_silhs=.false., grid_kind=stats_zt )

        k = k + 1

      case ('thl_2')
        stats_metadata%ithl_2 = k

        call stat_assign( var_index=stats_metadata%ithl_2, var_name="thl_2", &
             var_description="thl_2, pdf parameter: mean thl of component 2", var_units="K", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('varnce_thl_1')
        stats_metadata%ivarnce_thl_1 = k

        call stat_assign( var_index=stats_metadata%ivarnce_thl_1, var_name="varnce_thl_1", &
             var_description="varnce_thl_1, pdf parameter: thl variance of component 1", &
             var_units="K^2", &
             l_silhs=.false., grid_kind=stats_zt )

        k = k + 1

      case ('varnce_thl_2')
        stats_metadata%ivarnce_thl_2 = k
        call stat_assign( var_index=stats_metadata%ivarnce_thl_2, var_name="varnce_thl_2", &
             var_description="varnce_thl_2, pdf parameter: thl variance of component 2", &
             var_units="K^2", &
             l_silhs=.false., grid_kind=stats_zt )

        k = k + 1

      case ('rt_1')
        stats_metadata%irt_1 = k
        call stat_assign( var_index=stats_metadata%irt_1, var_name="rt_1", &
             var_description="rt_1, pdf parameter: mean rt of component 1", var_units="kg/kg", &
             l_silhs=.false., grid_kind=stats_zt )

        k = k + 1

      case ('rt_2')
        stats_metadata%irt_2 = k

        call stat_assign( var_index=stats_metadata%irt_2, var_name="rt_2", &
             var_description="rt_2, pdf parameter: mean rt of component 2", var_units="kg/kg", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('varnce_rt_1')
        stats_metadata%ivarnce_rt_1 = k
        call stat_assign( var_index=stats_metadata%ivarnce_rt_1, var_name="varnce_rt_1", &
             var_description="varnce_rt_1, pdf parameter: rt variance of component 1", &
             var_units="(kg^2)/(kg^2)", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('varnce_rt_2')
        stats_metadata%ivarnce_rt_2 = k

        call stat_assign( var_index=stats_metadata%ivarnce_rt_2, var_name="varnce_rt_2", &
             var_description="varnce_rt_2, pdf parameter: rt variance of component 2", &
             var_units="(kg^2)/(kg^2)", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rc_1')
        stats_metadata%irc_1 = k

        call stat_assign( var_index=stats_metadata%irc_1, var_name="rc_1", &
             var_description="rc_1, pdf parameter: mean rc of component 1", var_units="kg/kg", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rc_2')
        stats_metadata%irc_2 = k

        call stat_assign( var_index=stats_metadata%irc_2, var_name="rc_2", &
             var_description="rc_2, pdf parameter: mean rc of component 2", var_units="kg/kg", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rsatl_1')
        stats_metadata%irsatl_1 = k

        call stat_assign( var_index=stats_metadata%irsatl_1, var_name="rsatl_1", &
             var_description="rsatl_1, pdf parameter: sat mix rat based on tl1", &
             var_units="kg/kg", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rsatl_2')
        stats_metadata%irsatl_2 = k

        call stat_assign( var_index=stats_metadata%irsatl_2, var_name="rsatl_2", &
             var_description="rsatl_2, pdf parameter: sat mix rat based on tl2", &
             var_units="kg/kg", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('cloud_frac_1')
        stats_metadata%icloud_frac_1 = k
        call stat_assign( var_index=stats_metadata%icloud_frac_1, var_name="cloud_frac_1", &
             var_description="cloud_frac_1, pdf parameter cloud_frac_1", var_units="-", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('cloud_frac_2')
        stats_metadata%icloud_frac_2 = k

        call stat_assign( var_index=stats_metadata%icloud_frac_2, var_name="cloud_frac_2", &
             var_description="cloud_frac_2, pdf parameter cloud_frac_2", var_units="-", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('chi_1')
        stats_metadata%ichi_1 = k

        call stat_assign( var_index=stats_metadata%ichi_1, var_name="chi_1", &
             var_description="chi_1, pdf parameter: Mellor's s (extended liq) for component 1", &
             var_units="kg/kg", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('chi_2')
        stats_metadata%ichi_2 = k

        call stat_assign( var_index=stats_metadata%ichi_2, var_name="chi_2", &
             var_description="chi_2, pdf parameter: Mellor's s (extended liq) for component 2", &
             var_units="kg/kg", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('stdev_chi_1')
        stats_metadata%istdev_chi_1 = k

        call stat_assign( var_index=stats_metadata%istdev_chi_1, var_name="stdev_chi_1", &
             var_description="stdev_chi_1, pdf parameter: Std dev of chi_1", var_units="kg/kg", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('stdev_chi_2')
        stats_metadata%istdev_chi_2 = k

        call stat_assign( var_index=stats_metadata%istdev_chi_2, var_name="stdev_chi_2", &
             var_description="stdev_chi_2, pdf parameter: Std dev of chi_2", var_units="kg/kg", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('chip2')
        stats_metadata%ichip2 = k
        call stat_assign( var_index=stats_metadata%ichip2, var_name="chip2", &
             var_description="chip2, Variance of chi(s) (overall)", var_units="(kg/kg)^2", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('stdev_eta_1')
        stats_metadata%istdev_eta_1 = k

        call stat_assign( var_index=stats_metadata%istdev_eta_1, var_name="stdev_eta_1", &
             var_description="stdev_eta_1, Standard dev. of eta(t) (1st PDF component)", &
             var_units="kg/kg", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('stdev_eta_2')
        stats_metadata%istdev_eta_2 = k

        call stat_assign( var_index=stats_metadata%istdev_eta_2, var_name="stdev_eta_2", &
             var_description="stdev_eta_2, Standard dev. of eta(t) (2nd PDF component)", &
             var_units="kg/kg", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('covar_chi_eta_1')
        stats_metadata%icovar_chi_eta_1 = k

        call stat_assign( var_index=stats_metadata%icovar_chi_eta_1, var_name="covar_chi_eta_1", &
             var_description="covar_chi_eta_1, Covariance of chi(s) and eta(t) " &
             // "(1st PDF component)", &
             var_units="kg^2/kg^2", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('covar_chi_eta_2')
        stats_metadata%icovar_chi_eta_2 = k

        call stat_assign( var_index=stats_metadata%icovar_chi_eta_2, var_name="covar_chi_eta_2", &
             var_description="covar_chi_eta_2, Covariance of chi(s) and eta(t) " &
             // "(2nd PDF component)", &
             var_units="kg^2/kg^2", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('corr_w_chi_1')
        stats_metadata%icorr_w_chi_1 = k

        call stat_assign( var_index=stats_metadata%icorr_w_chi_1, var_name="corr_w_chi_1", &
                          var_description="corr_w_chi_1, Correlation of w and chi (s)" &
                          // " (1st PDF component)", var_units="-", &
                          l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('corr_w_chi_2')
        stats_metadata%icorr_w_chi_2 = k

        call stat_assign( var_index=stats_metadata%icorr_w_chi_2, var_name="corr_w_chi_2", &
                          var_description="corr_w_chi_2, Correlation of w and chi (s)" &
                          // " (2nd PDF component)", var_units="-", &
                          l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('corr_w_eta_1')
        stats_metadata%icorr_w_eta_1 = k

        call stat_assign( var_index=stats_metadata%icorr_w_eta_1, var_name="corr_w_eta_1", &
                          var_description="corr_w_eta_1, Correlation of w and eta (t)" &
                          // " (1st PDF component)", var_units="-", &
                          l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('corr_w_eta_2')
        stats_metadata%icorr_w_eta_2 = k

        call stat_assign( var_index=stats_metadata%icorr_w_eta_2, var_name="corr_w_eta_2", &
                          var_description="corr_w_eta_2, Correlation of w and eta (t)" &
                          // " (2nd PDF component)", var_units="-", &
                          l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('corr_chi_eta_1')
        stats_metadata%icorr_chi_eta_1 = k

        call stat_assign( var_index=stats_metadata%icorr_chi_eta_1, &
                          var_name="corr_chi_eta_1", &
                          var_description="corr_chi_eta_1, Correlation of chi (s) and" &
                          // " eta (t) (1st PDF component)", &
                          var_units="-", &
                          l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('corr_chi_eta_2')
        stats_metadata%icorr_chi_eta_2 = k

        call stat_assign( var_index=stats_metadata%icorr_chi_eta_2, &
                          var_name="corr_chi_eta_2", &
                          var_description="corr_chi_eta_2, Correlation of chi (s) and" &
                          // " eta (t) (2nd PDF component)", &
                          var_units="-", &
                          l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('corr_w_rt_1')
        stats_metadata%icorr_w_rt_1 = k

        call stat_assign( var_index=stats_metadata%icorr_w_rt_1, var_name="corr_w_rt_1", &
                          var_description="corr_w_rt_1, Correlation of w and rt" &
                          // " (1st PDF component)", var_units="-", &
                          l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('corr_w_rt_2')
        stats_metadata%icorr_w_rt_2 = k

        call stat_assign( var_index=stats_metadata%icorr_w_rt_2, var_name="corr_w_rt_2", &
                          var_description="corr_w_rt_2, Correlation of w and rt" &
                          // " (2nd PDF component)", var_units="-", &
                          l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('corr_w_thl_1')
        stats_metadata%icorr_w_thl_1 = k

        call stat_assign( var_index=stats_metadata%icorr_w_thl_1, var_name="corr_w_thl_1", &
                          var_description="corr_w_thl_1, Correlation of w and thl" &
                          // " (1st PDF component)", var_units="-", &
                          l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('corr_w_thl_2')
        stats_metadata%icorr_w_thl_2 = k

        call stat_assign( var_index=stats_metadata%icorr_w_thl_2, var_name="corr_w_thl_2", &
                          var_description="corr_w_thl_2, Correlation of w and thl" &
                          // " (2nd PDF component)", var_units="-", &
                          l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('corr_rt_thl_1')
        stats_metadata%icorr_rt_thl_1 = k

        call stat_assign( var_index=stats_metadata%icorr_rt_thl_1, var_name="corr_rt_thl_1", &
                          var_description="corr_rt_thl_1, Correlation of rt and thl" &
                          // " (1st PDF component)", var_units="-", &
                          l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('corr_rt_thl_2')
        stats_metadata%icorr_rt_thl_2 = k

        call stat_assign( var_index=stats_metadata%icorr_rt_thl_2, var_name="corr_rt_thl_2", &
                          var_description="corr_rt_thl_2, Correlation of rt and thl" &
                          // " (2nd PDF component)", var_units="-", &
                          l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('crt_1')
        stats_metadata%icrt_1 = k

        call stat_assign( var_index=stats_metadata%icrt_1, var_name="crt_1", &
                          var_description="crt_1, Coefficient on rt in chi/eta" &
                          // " equations (1st PDF comp.)", &
                          var_units="-", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('crt_2')
        stats_metadata%icrt_2 = k

        call stat_assign( var_index=stats_metadata%icrt_2, var_name="crt_2", &
                          var_description="crt_2, Coefficient on rt in chi/eta" &
                          // " equations (2nd PDF comp.)", &
                          var_units="-", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('cthl_1')
        stats_metadata%icthl_1 = k

        call stat_assign( var_index=stats_metadata%icthl_1, var_name="cthl_1", &
                          var_description="cthl_1, Coefficient on theta-l in chi/eta" &
                          // " equations (1st PDF comp.)", &
                          var_units="kg/kg/K", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('cthl_2')
        stats_metadata%icthl_2 = k

        call stat_assign( var_index=stats_metadata%icthl_2, var_name="cthl_2", &
                          var_description="cthl_2, Coefficient on theta-l in chi/eta" &
                          // " equations (2nd PDF comp.)", &
                          var_units="kg/kg/K", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case('F_w')
        stats_metadata%iF_w = k

        call stat_assign( var_index=stats_metadata%iF_w, var_name="F_w", &
                          var_description="F_w, Parameter for the spread of the" &
                          // " PDF component means of w (new PDF)", &
                          var_units="-", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case('F_rt')
        stats_metadata%iF_rt = k

        call stat_assign( var_index=stats_metadata%iF_rt, var_name="F_rt", &
                          var_description="F_rt, Parameter for the spread of the" &
                          // " PDF component means of rt (new PDF)", &
                          var_units="-", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case('F_thl')
        stats_metadata%iF_thl = k

        call stat_assign( var_index=stats_metadata%iF_thl, var_name="F_thl", &
                          var_description="F_thl, Parameter for the spread of the" &
                          // " PDF component means of thl (new PDF)", &
                          var_units="-", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case('min_F_w')
        stats_metadata%imin_F_w = k

        call stat_assign( var_index=stats_metadata%imin_F_w, var_name="min_F_w", &
                          var_description="min_F_w, Minimum allowable value of the" &
                          // " parameter F_w (new PDF)", &
                          var_units="-", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case('max_F_w')
        stats_metadata%imax_F_w = k

        call stat_assign( var_index=stats_metadata%imax_F_w, var_name="max_F_w", &
                          var_description="max_F_w, Maximum allowable value of the" &
                          // " parameter F_w (new PDF)", &
                          var_units="-", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case('min_F_rt')
        stats_metadata%imin_F_rt = k

        call stat_assign( var_index=stats_metadata%imin_F_rt, var_name="min_F_rt", &
                          var_description="min_F_rt, Minimum allowable value of the" &
                          // " parameter F_rt (new PDF)", &
                          var_units="-", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case('max_F_rt')
        stats_metadata%imax_F_rt = k

        call stat_assign( var_index=stats_metadata%imax_F_rt, var_name="max_F_rt", &
                          var_description="max_F_rt, Maximum allowable value of the" &
                          // " parameter F_rt (new PDF)", &
                          var_units="-", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case('min_F_thl')
        stats_metadata%imin_F_thl = k

        call stat_assign( var_index=stats_metadata%imin_F_thl, var_name="min_F_thl", &
                          var_description="min_F_thl, Minimum allowable value of the" &
                          // " parameter F_thl (new PDF)", &
                          var_units="-", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case('max_F_thl')
        stats_metadata%imax_F_thl = k

        call stat_assign( var_index=stats_metadata%imax_F_thl, var_name="max_F_thl", &
                          var_description="max_F_thl, Maximum allowable value of the" &
                          // " parameter F_thl (new PDF)", &
                          var_units="-", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ( 'coef_wprtp2_implicit' )
        stats_metadata%icoef_wprtp2_implicit = k
        call stat_assign( var_index=stats_metadata%icoef_wprtp2_implicit, &
                          var_name="coef_wprtp2_implicit", &
                          var_description="coef_wprtp2_implicit, wprtp2" &
                                          // " = coef_wprtp2_implicit" &
                                          // " * rtp2" &
                                          // " + term_wprtp2_explicit", &
                          var_units="m/s", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ( 'term_wprtp2_explicit' )
        stats_metadata%iterm_wprtp2_explicit = k
        call stat_assign( var_index=stats_metadata%iterm_wprtp2_explicit, &
                          var_name="term_wprtp2_explicit", &
                          var_description="term_wprtp2_explicit, wprtp2" &
                                          // " = coef_wprtp2_implicit" &
                                          // " * rtp2" &
                                          // " + term_wprtp2_explicit", &
                          var_units="m/s kg^2/kg^2", l_silhs=.false., &
                          grid_kind=stats_zt )
        k = k + 1

      case ( 'coef_wpthlp2_implicit' )
        stats_metadata%icoef_wpthlp2_implicit = k
        call stat_assign( var_index=stats_metadata%icoef_wpthlp2_implicit, &
                          var_name="coef_wpthlp2_implicit", &
                          var_description="coef_wpthlp2_implicit, wpthlp2" &
                                          // " = coef_wpthlp2_implicit" &
                                          // " * thlp2" &
                                          // " + term_wpthlp2_explicit", &
                          var_units="m/s", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ( 'term_wpthlp2_explicit' )
        stats_metadata%iterm_wpthlp2_explicit = k
        call stat_assign( var_index=stats_metadata%iterm_wpthlp2_explicit, &
                          var_name="term_wpthlp2_explicit", &
                          var_description="term_wpthlp2_explicit, wpthlp2" &
                                          // " = coef_wpthlp2_implicit" &
                                          // " * thlp2" &
                                          // " + term_wpthlp2_explicit", &
                          var_units="m/s K^2", l_silhs=.false., &
                          grid_kind=stats_zt )
        k = k + 1

      case ( 'coef_wprtpthlp_implicit' )
        stats_metadata%icoef_wprtpthlp_implicit = k
        call stat_assign( var_index=stats_metadata%icoef_wprtpthlp_implicit, &
                          var_name="coef_wprtpthlp_implicit", &
                          var_description="coef_wprtpthlp_implicit, wprtpthlp" &
                                          // " = coef_wprtpthlp_implicit" &
                                          // " * rtpthlp" &
                                          // " + term_wprtpthlp_explicit", &
                          var_units="m/s", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ( 'term_wprtpthlp_explicit' )
        stats_metadata%iterm_wprtpthlp_explicit = k
        call stat_assign( var_index=stats_metadata%iterm_wprtpthlp_explicit, &
                          var_name="term_wprtpthlp_explicit", &
                          var_description="term_wprtpthlp_explicit, wprtpthlp" &
                                          // " = coef_wprtpthlp_implicit" &
                                          // " * rtpthlp" &
                                          // " + term_wprtpthlp_explicit]", &
                          var_units="m/s (kg/kg) K", l_silhs=.false., &
                          grid_kind=stats_zt )
        k = k + 1

      case ( 'coef_wp2rtp_implicit' )
        stats_metadata%icoef_wp2rtp_implicit = k
        call stat_assign( var_index=stats_metadata%icoef_wp2rtp_implicit, &
                          var_name="coef_wp2rtp_implicit", &
                          var_description="coef_wp2rtp_implicit, wp2rtp" &
                                          // " = coef_wp2rtp_implicit" &
                                          // " * wprtp" &
                                          // " + term_wp2rtp_explicit", &
                          var_units="m/s", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ( 'term_wp2rtp_explicit' )
        stats_metadata%iterm_wp2rtp_explicit = k
        call stat_assign( var_index=stats_metadata%iterm_wp2rtp_explicit, &
                          var_name="term_wp2rtp_explicit", &
                          var_description="term_wp2rtp_explicit, wp2rtp" &
                                          // " = coef_wp2rtp_implicit" &
                                          // " * wprtp" &
                                          // " + term_wp2rtp_explicit", &
                          var_units="m^2/s^2 kg/kg", l_silhs=.false., &
                          grid_kind=stats_zt )
        k = k + 1

      case ( 'coef_wp2thlp_implicit' )
        stats_metadata%icoef_wp2thlp_implicit = k
        call stat_assign( var_index=stats_metadata%icoef_wp2thlp_implicit, &
                          var_name="coef_wp2thlp_implicit", &
                          var_description="coef_wp2thlp_implicit, wp2thlp" &
                                          // " = coef_wp2thlp_implicit" &
                                          // " * wpthlp" &
                                          // " + term_wp2thlp_explicit", &
                          var_units="m/s", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ( 'term_wp2thlp_explicit' )
        stats_metadata%iterm_wp2thlp_explicit = k
        call stat_assign( var_index=stats_metadata%iterm_wp2thlp_explicit, &
                          var_name="term_wp2thlp_explicit", &
                          var_description="term_wp2thlp_explicit, wp2thlp" &
                                          // " = coef_wp2thlp_implicit" &
                                          // " * wpthlp" &
                                          // " + term_wp2thlp_explicit", &
                          var_units="m^2/s^2 K", l_silhs=.false., &
                          grid_kind=stats_zt )
        k = k + 1

      case('wp2_zt')
        stats_metadata%iwp2_zt = k

        call stat_assign( var_index=stats_metadata%iwp2_zt, var_name="wp2_zt", &
             var_description="wp2_zt, w'^2 interpolated to thermodynamic levels", &
             var_units="m^2/s^2", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case('thlp2_zt')
        stats_metadata%ithlp2_zt = k

        call stat_assign( var_index=stats_metadata%ithlp2_zt, var_name="thlp2_zt", &
             var_description="thlp2_zt, thl'^2 interpolated to thermodynamic levels", &
             var_units="K^2", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case('wpthlp_zt')
        stats_metadata%iwpthlp_zt = k

        call stat_assign( var_index=stats_metadata%iwpthlp_zt, var_name="wpthlp_zt", &
             var_description="wpthlp_zt, w'thl' interpolated to thermodynamic levels", &
             var_units="(m K)/s", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case('wprtp_zt')
        stats_metadata%iwprtp_zt = k

        call stat_assign( var_index=stats_metadata%iwprtp_zt, var_name="wprtp_zt", &
             var_description="wprtp_zt, w'rt' interpolated to thermodynamic levels", &
             var_units="(m kg)/(s kg)", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case('rtp2_zt')
        stats_metadata%irtp2_zt = k

        call stat_assign( var_index=stats_metadata%irtp2_zt, var_name="rtp2_zt", &
             var_description="rtp2_zt, rt'^2 interpolated to thermodynamic levels", &
             var_units="(kg/kg)^2", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case('rtpthlp_zt')
        stats_metadata%irtpthlp_zt = k

        call stat_assign( var_index=stats_metadata%irtpthlp_zt, var_name="rtpthlp_zt", &
             var_description="rtpthlp_zt, rt'thl' interpolated to thermodynamic levels", &
             var_units="(kg K)/kg", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('up2_zt')
        stats_metadata%iup2_zt = k
        call stat_assign( var_index=stats_metadata%iup2_zt, var_name="up2_zt", &
             var_description="up2_zt, u'^2 interpolated to thermodynamic levels", &
             var_units="m^2/s^2", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('vp2_zt')
        stats_metadata%ivp2_zt = k
        call stat_assign( var_index=stats_metadata%ivp2_zt, var_name="vp2_zt", &
             var_description="vp2_zt, v'^2 interpolated to thermodynamic levels", &
             var_units="m^2/s^2", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('upwp_zt')
        stats_metadata%iupwp_zt = k
        call stat_assign( var_index=stats_metadata%iupwp_zt, var_name="upwp_zt", &
             var_description="upwp_zt, u'w' interpolated to thermodynamic levels", &
             var_units="m^2/s^2", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('vpwp_zt')
        stats_metadata%ivpwp_zt = k
        call stat_assign( var_index=stats_metadata%ivpwp_zt, var_name="vpwp_zt", &
             var_description="vpwp_zt, v'w' interpolated to thermodynamic levels", &
             var_units="m^2/s^2", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('Skw_zt')
        stats_metadata%iSkw_zt = k
        call stat_assign( var_index=stats_metadata%iSkw_zt, var_name="Skw_zt", &
             var_description="Skw_zt, Skewness of w on thermodynamic levels", &
             var_units="-", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('Skthl_zt')
        stats_metadata%iSkthl_zt = k
        call stat_assign( var_index=stats_metadata%iSkthl_zt, var_name="Skthl_zt", &
             var_description="Skthl_zt, Skewness of thl on thermodynamic levels", &
             var_units="-", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('Skrt_zt')
        stats_metadata%iSkrt_zt = k
        call stat_assign( var_index=stats_metadata%iSkrt_zt, var_name="Skrt_zt", &
             var_description="Skrt_zt, Skewness of rt on thermodynamic levels", &
             var_units="-", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rcm_supersat_adj')
        stats_metadata%ircm_supersat_adj = k
        call stat_assign( var_index=stats_metadata%ircm_supersat_adj, var_name="rcm_supersat_adj", &
             var_description="rcm_supersat_adj, rcm adjustment due to spurious supersaturation", &
             var_units="kg/kg", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      ! Hydrometeor overall variances for each hydrometeor type.
      case('hmp2_zt')

         do hm_idx = 1, hydromet_dim, 1

            hm_type = hydromet_list(hm_idx)

            ! The overall variance of the hydrometeor.
            stats_metadata%ihmp2_zt(hm_idx) = k

            if ( l_mix_rat_hm(hm_idx) ) then

               call stat_assign( var_index=stats_metadata%ihmp2_zt(hm_idx), &
                                 var_name=trim( hm_type(1:2) )//"p2_zt", &
                                 var_description="<" &
                                 // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                                 // "'^2> on thermodynamic levels (from " &
                                 // "integration over PDF) [(kg/kg)^2]", &
                                 var_units="(kg/kg)^2", &
                                 l_silhs=.false., grid_kind=stats_zt )

            else ! Concentration

               call stat_assign( var_index=stats_metadata%ihmp2_zt(hm_idx), &
                                 var_name=trim( hm_type(1:2) )//"p2_zt", &
                                 var_description="<" &
                                 // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                                 // "'^2> on thermodynamic levels (from " &
                                 // "integration over PDF) [(num/kg)^2]", &
                                 var_units="(num/kg)^2", &
                                 l_silhs=.false., grid_kind=stats_zt )

            endif ! l_mix_rat_hm(hm_idx)

            k = k + 1

         enddo ! hm_idx = 1, hydromet_dim, 1

      case ('C11_Skw_fnc')
        stats_metadata%iC11_Skw_fnc = k

        call stat_assign( var_index=stats_metadata%iC11_Skw_fnc, var_name="C11_Skw_fnc", &
             var_description="C11_Skw_fnc, C_11 parameter with Sk_w applied", var_units="count", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('chi')
        stats_metadata%ichi = k

        call stat_assign( var_index=stats_metadata%ichi, var_name="chi", &
             var_description="chi, Mellor's s (extended liq)", var_units="kg/kg", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ( 'a3_coef_zt' )
        stats_metadata%ia3_coef_zt = k
        call stat_assign( var_index=stats_metadata%ia3_coef_zt, var_name="a3_coef_zt", &
             var_description="a3_coef_zt, The a3 coefficient interpolated the the zt grid", &
             var_units="count", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ( 'wp3_on_wp2_zt' )
        stats_metadata%iwp3_on_wp2_zt = k
        call stat_assign( var_index=stats_metadata%iwp3_on_wp2_zt, var_name="wp3_on_wp2_zt", &
             var_description="wp3_on_wp2_zt, Smoothed version of wp3 / wp2", var_units="m/s", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      ! Hydrometeor component mean values for each PDF component and hydrometeor
      ! type.
      case ( "hm_i" )

         do hm_idx = 1, hydromet_dim, 1

            hm_type = hydromet_list(hm_idx)

            ! The mean of the hydrometeor in the 1st PDF component.
            stats_metadata%ihm_1(hm_idx) = k

            if ( l_mix_rat_hm(hm_idx) ) then

               call stat_assign( var_index=stats_metadata%ihm_1(hm_idx), &
                                 var_name=trim( hm_type(1:2) )//"_1", &
                                 var_description="Mean of " &
                                 // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                                 // " (1st PDF component) [kg/kg]", &
                                 var_units="kg/kg", &
                                 l_silhs=.false., grid_kind=stats_zt )

            else ! Concentration

               call stat_assign( var_index=stats_metadata%ihm_1(hm_idx), &
                                 var_name=trim( hm_type(1:2) )//"_1", &
                                 var_description="Mean of " &
                                 // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                                 // " (1st PDF component) [num/kg]", &
                                 var_units="num/kg", &
                                 l_silhs=.false., grid_kind=stats_zt )

            endif ! l_mix_rat_hm(hm_idx)

            k = k + 1

            ! The mean of the hydrometeor in the 2nd PDF component.
            stats_metadata%ihm_2(hm_idx) = k

            if ( l_mix_rat_hm(hm_idx) ) then

               call stat_assign( var_index=stats_metadata%ihm_2(hm_idx), &
                                 var_name=trim( hm_type(1:2) )//"_2", &
                                 var_description="Mean of " &
                                 // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                                 // " (2nd PDF component) [kg/kg]", &
                                 var_units="kg/kg", &
                                 l_silhs=.false., grid_kind=stats_zt )

            else ! Concentration

               call stat_assign( var_index=stats_metadata%ihm_2(hm_idx), &
                                 var_name=trim( hm_type(1:2) )//"_2", &
                                 var_description="Mean of " &
                                 // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                                 // " (2nd PDF component) [num/kg]", &
                                 var_units="num/kg", &
                                 l_silhs=.false., grid_kind=stats_zt )

            endif ! l_mix_rat_hm(hm_idx)

            k = k + 1

         enddo ! hm_idx = 1, hydromet_dim, 1

      case ( 'precip_frac' )
        stats_metadata%iprecip_frac = k
        call stat_assign( var_index=stats_metadata%iprecip_frac, var_name="precip_frac", &
             var_description="precip_frac, Precipitation Fraction", var_units="-", &
             l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ( 'precip_frac_1' )
        stats_metadata%iprecip_frac_1 = k
        call stat_assign( var_index=stats_metadata%iprecip_frac_1, var_name="precip_frac_1", &
             var_description="precip_frac_1, Precipitation Fraction (1st PDF component)", &
             var_units="-", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ( 'precip_frac_2' )
        stats_metadata%iprecip_frac_2 = k
        call stat_assign( var_index=stats_metadata%iprecip_frac_2, var_name="precip_frac_2", &
             var_description="precip_frac_2, Precipitation Fraction (2nd PDF component)", &
             var_units="-", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ( 'Ncnm' )
        stats_metadata%iNcnm = k
        call stat_assign( var_index=stats_metadata%iNcnm, var_name="Ncnm", &
             var_description="Ncnm, Cloud nuclei concentration (PDF)", &
             var_units="num/kg", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      ! Hydrometeor component mean values (in-precip) for each PDF component and
      ! hydrometeor type.
      case ( 'mu_hm_i' )

         do hm_idx = 1, hydromet_dim, 1

            hm_type = hydromet_list(hm_idx)

            ! The in-precip mean of the hydrometeor in the 1st PDF component.
            stats_metadata%imu_hm_1(hm_idx) = k

            if ( l_mix_rat_hm(hm_idx) ) then

               call stat_assign( var_index=stats_metadata%imu_hm_1(hm_idx), &
                                 var_name="mu_"//trim( hm_type(1:2) )//"_1", &
                                 var_description="Mean (in-precip) of " &
                                 // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                                 // " (1st PDF component) [kg/kg]", &
                                 var_units="kg/kg", &
                                 l_silhs=.false., grid_kind=stats_zt )

            else ! Concentration

               call stat_assign( var_index=stats_metadata%imu_hm_1(hm_idx), &
                                 var_name="mu_"//trim( hm_type(1:2) )//"_1", &
                                 var_description="Mean (in-precip) of " &
                                 // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                                 // " (1st PDF component) [num/kg]", &
                                 var_units="num/kg", &
                                 l_silhs=.false., grid_kind=stats_zt )

            endif ! l_mix_rat_hm(hm_idx)

            k = k + 1

            ! The in-precip mean of the hydrometeor in the 2nd PDF component.
            stats_metadata%imu_hm_2(hm_idx) = k

            if ( l_mix_rat_hm(hm_idx) ) then

               call stat_assign( var_index=stats_metadata%imu_hm_2(hm_idx), &
                                 var_name="mu_"//trim( hm_type(1:2) )//"_2", &
                                 var_description="Mean (in-precip) of " &
                                 // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                                 // " (2nd PDF component) [kg/kg]", &
                                 var_units="kg/kg", &
                                 l_silhs=.false., grid_kind=stats_zt )

            else ! Concentration

               call stat_assign( var_index=stats_metadata%imu_hm_2(hm_idx), &
                                 var_name="mu_"//trim( hm_type(1:2) )//"_2", &
                                 var_description="Mean (in-precip) of " &
                                 // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                                 // " (2nd PDF component) [num/kg]", &
                                 var_units="num/kg", &
                                 l_silhs=.false., grid_kind=stats_zt )

            endif ! l_mix_rat_hm(hm_idx)

            k = k + 1

         enddo ! hm_idx = 1, hydromet_dim, 1

      case ( 'mu_Ncn_i' )
        stats_metadata%imu_Ncn_1 = k

         call stat_assign( var_index=stats_metadata%imu_Ncn_1, &
                           var_name="mu_Ncn_1", &
                           var_description="mu_Ncn_1, Mean of N_cn (1st PDF component) " &
                           // "[num/kg]", var_units="num/kg", &
                           l_silhs=.false., grid_kind=stats_zt )

         k = k + 1

         stats_metadata%imu_Ncn_2 = k

         call stat_assign( var_index=stats_metadata%imu_Ncn_2, &
                           var_name="mu_Ncn_2", &
                           var_description="mu_Ncn_2, Mean of N_cn (2nd PDF component)", &
                           var_units="num/kg", &
                           l_silhs=.false., grid_kind=stats_zt )

         k = k + 1

      ! Hydrometeor component mean values (in-precip) for ln hm for each PDF
      ! component and hydrometeor type.
      case ( 'mu_hm_i_n' )

         do hm_idx = 1, hydromet_dim, 1

            hm_type = hydromet_list(hm_idx)

            ! The in-precip mean of ln hm in the 1st PDF component.
            stats_metadata%imu_hm_1_n(hm_idx) = k

            if ( l_mix_rat_hm(hm_idx) ) then

               call stat_assign( var_index=stats_metadata%imu_hm_1_n(hm_idx), &
                                 var_name="mu_"//trim( hm_type(1:2) )//"_1_n", &
                                 var_description="Mean (in-precip) of ln " &
                                 // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                                 // " (1st PDF component) [ln(kg/kg)]", &
                                 var_units="ln(kg/kg)", &
                                 l_silhs=.false., grid_kind=stats_zt )

            else ! Concentration

               call stat_assign( var_index=stats_metadata%imu_hm_1_n(hm_idx), &
                                 var_name="mu_"//trim( hm_type(1:2) )//"_1_n", &
                                 var_description="Mean (in-precip) of ln " &
                                 // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                                 // " (1st PDF component) [ln(num/kg)]", &
                                 var_units="ln(num/kg)", &
                                 l_silhs=.false., grid_kind=stats_zt )

            endif ! l_mix_rat_hm(hm_idx)

            k = k + 1

            ! The in-precip mean of ln hm in the 2nd PDF component.
            stats_metadata%imu_hm_2_n(hm_idx) = k

            if ( l_mix_rat_hm(hm_idx) ) then

               call stat_assign( var_index=stats_metadata%imu_hm_2_n(hm_idx), &
                                 var_name="mu_"//trim( hm_type(1:2) )//"_2_n", &
                                 var_description="Mean (in-precip) of ln " &
                                 // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                                 // " (2nd PDF component) [ln(kg/kg)]", &
                                 var_units="ln(kg/kg)", &
                                 l_silhs=.false., grid_kind=stats_zt )

            else ! Concentration

               call stat_assign( var_index=stats_metadata%imu_hm_2_n(hm_idx), &
                                 var_name="mu_"//trim( hm_type(1:2) )//"_2_n", &
                                 var_description="Mean (in-precip) of ln " &
                                 // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                                 // " (2nd PDF component) [ln(num/kg)]", &
                                 var_units="ln(num/kg)", &
                                 l_silhs=.false., grid_kind=stats_zt )

            endif ! l_mix_rat_hm(hm_idx)

            k = k + 1

         enddo ! hm_idx = 1, hydromet_dim, 1

      case ( 'mu_Ncn_i_n' )

         stats_metadata%imu_Ncn_1_n = k

         call stat_assign( var_index=stats_metadata%imu_Ncn_1_n, &
                           var_name="mu_Ncn_1_n", &
                           var_description="mu_Ncn_1_n, Mean of ln N_cn " &
                           // "(1st PDF component)", &
                           var_units="ln(num/kg)", &
                           l_silhs=.false., grid_kind=stats_zt )

         k = k + 1

         stats_metadata%imu_Ncn_2_n = k

         call stat_assign( var_index=stats_metadata%imu_Ncn_2_n, &
                           var_name="mu_Ncn_2_n", &
                           var_description="mu_Ncn_2_n, Mean of ln N_cn " &
                           // "(2nd PDF component)", &
                           var_units="ln(num/kg)", &
                           l_silhs=.false., grid_kind=stats_zt )

         k = k + 1

      ! Hydrometeor component standard deviations (in-precip) for each PDF
      ! component and hydrometeor type.
      case ( 'sigma_hm_i' )

         do hm_idx = 1, hydromet_dim, 1

            hm_type = hydromet_list(hm_idx)

            ! The in-precip standard deviation of the hydrometeor in the 1st PDF
            ! component.
            stats_metadata%isigma_hm_1(hm_idx) = k

            if ( l_mix_rat_hm(hm_idx) ) then

               call stat_assign( var_index=stats_metadata%isigma_hm_1(hm_idx), &
                                 var_name="sigma_" &
                                 // trim( hm_type(1:2) )//"_1", &
                                 var_description="Standard deviation " &
                                 // "(in-precip) of " &
                                 // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                                 // " (1st PDF component) [kg/kg]", &
                                 var_units="kg/kg", &
                                 l_silhs=.false., grid_kind=stats_zt )

            else ! Concentration

               call stat_assign( var_index=stats_metadata%isigma_hm_1(hm_idx), &
                                 var_name="sigma_" &
                                 // trim( hm_type(1:2) )//"_1", &
                                 var_description="Standard deviation " &
                                 // "(in-precip) of " &
                                 // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                                 // " (1st PDF component) [num/kg]", &
                                 var_units="num/kg", &
                                 l_silhs=.false., grid_kind=stats_zt )

            endif ! l_mix_rat_hm(hm_idx)

            k = k + 1

            ! The in-precip standard deviation of the hydrometeor in the 2nd PDF
            ! component.
            stats_metadata%isigma_hm_2(hm_idx) = k

            if ( l_mix_rat_hm(hm_idx) ) then

               call stat_assign( var_index=stats_metadata%isigma_hm_2(hm_idx), &
                                 var_name="sigma_" &
                                 // trim( hm_type(1:2) )//"_2", &
                                 var_description="Standard deviation " &
                                 // "(in-precip) of " &
                                 // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                                 // " (2nd PDF component) [kg/kg]", &
                                 var_units="kg/kg", &
                                 l_silhs=.false., grid_kind=stats_zt )

            else ! Concentration

               call stat_assign( var_index=stats_metadata%isigma_hm_2(hm_idx), &
                                 var_name="sigma_" &
                                 // trim( hm_type(1:2) )//"_2", &
                                 var_description="Standard deviation " &
                                 // "(in-precip) of " &
                                 // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                                 // " (2nd PDF component) [num/kg]", &
                                 var_units="num/kg", &
                                 l_silhs=.false., grid_kind=stats_zt )

            endif ! l_mix_rat_hm(hm_idx)

            k = k + 1

         enddo ! hm_idx = 1, hydromet_dim, 1

      case ( 'sigma_Ncn_i' )

         stats_metadata%isigma_Ncn_1 = k

         call stat_assign( var_index=stats_metadata%isigma_Ncn_1, &
                           var_name="sigma_Ncn_1", &
                           var_description="sigma_Ncn_1, Standard deviation of N_cn " &
                           // "(1st PDF component)", &
                           var_units="num/kg", l_silhs=.false., grid_kind=stats_zt )

         k = k + 1

         stats_metadata%isigma_Ncn_2 = k

         call stat_assign( var_index=stats_metadata%isigma_Ncn_2, &
                           var_name="sigma_Ncn_2", &
                           var_description="sigma_Ncn_2, Standard deviation of N_cn " &
                           // "(2nd PDF component)", &
                           var_units="num/kg", l_silhs=.false., grid_kind=stats_zt )

         k = k + 1

      ! Hydrometeor component standard deviations (in-precip) for ln hm for each
      ! PDF component and hydrometeor type.
      case ( 'sigma_hm_i_n' )

         do hm_idx = 1, hydromet_dim, 1

            hm_type = hydromet_list(hm_idx)

            ! The in-precip standard deviation of ln hm in the 1st PDF
            ! component.
            stats_metadata%isigma_hm_1_n(hm_idx) = k

            call stat_assign( var_index=stats_metadata%isigma_hm_1_n(hm_idx), &
                              var_name="sigma_" &
                              // trim( hm_type(1:2) )//"_1_n", &
                              var_description="Standard deviation " &
                              // "(in-precip) of ln " &
                              // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                              // " (1st PDF component) [-]", &
                              var_units="-", &
                              l_silhs=.false., grid_kind=stats_zt )

            k = k + 1

            ! The in-precip standard deviation of ln hm in the 2nd PDF
            ! component.
            stats_metadata%isigma_hm_2_n(hm_idx) = k

            call stat_assign( var_index=stats_metadata%isigma_hm_2_n(hm_idx), &
                              var_name="sigma_" &
                              // trim( hm_type(1:2) )//"_2_n", &
                              var_description="Standard deviation " &
                              // "(in-precip) of ln " &
                              // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                              // " (2nd PDF component) [-]", &
                              var_units="-", &
                              l_silhs=.false., grid_kind=stats_zt )

            k = k + 1

         enddo ! hm_idx = 1, hydromet_dim, 1

      case ( 'sigma_Ncn_i_n' )
        stats_metadata%isigma_Ncn_1_n = k

         call stat_assign( var_index=stats_metadata%isigma_Ncn_1_n, &
                           var_name="sigma_Ncn_1_n", &
                           var_description="sigma_Ncn_1_n, Standard deviation of ln N_cn " &
                           // "(1st PDF component)", &
                           var_units="-", l_silhs=.false., grid_kind=stats_zt )

         k = k + 1

         stats_metadata%isigma_Ncn_2_n = k

         call stat_assign( var_index=stats_metadata%isigma_Ncn_2_n, &
                           var_name="sigma_Ncn_2_n", &
                           var_description="sigma_Ncn_2_n, Standard deviation of ln N_cn " &
                           // "(2nd PDF component)", &
                           var_units="-", l_silhs=.false., grid_kind=stats_zt )

         k = k + 1

      case ('corr_w_chi_1_ca')
        stats_metadata%icorr_w_chi_1_ca = k

        call stat_assign( var_index=stats_metadata%icorr_w_chi_1_ca, &
                          var_name="corr_w_chi_1_ca", &
                          var_description="corr_w_chi_1_ca, Correlation of w and chi" &
                          // " (1st PDF component) found in the correlation" &
                          // " array", var_units="-", &
                          l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('corr_w_chi_2_ca')
        stats_metadata%icorr_w_chi_2_ca = k

        call stat_assign( var_index=stats_metadata%icorr_w_chi_2_ca, &
                          var_name="corr_w_chi_2_ca", &
                          var_description="corr_w_chi_2_ca, Correlation of w and chi" &
                          // " (2nd PDF component) found in the correlation" &
                          // " array", var_units="-", &
                          l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('corr_w_eta_1_ca')
        stats_metadata%icorr_w_eta_1_ca = k

        call stat_assign( var_index=stats_metadata%icorr_w_eta_1_ca, &
                          var_name="corr_w_eta_1_ca", &
                          var_description="corr_w_eta_1_ca, Correlation of w and eta" &
                          // " (1st PDF component) found in the correlation" &
                          // " array", var_units="-", &
                          l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('corr_w_eta_2_ca')
        stats_metadata%icorr_w_eta_2_ca = k

        call stat_assign( var_index=stats_metadata%icorr_w_eta_2_ca, &
                          var_name="corr_w_eta_2_ca", &
                          var_description="corr_w_eta_2_ca, Correlation of w and eta" &
                          // " (2nd PDF component) found in the correlation" &
                          // " array", var_units="-", &
                          l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      ! Correlation of w and a hydrometeor (in-precip) for each PDF
      ! component and hydrometeor type.
      case ( 'corr_w_hm_i' )

         do hm_idx = 1, hydromet_dim, 1

            hm_type = hydromet_list(hm_idx)

            ! The in-precip correlation of w and the hydrometeor in the
            ! 1st PDF component.
            stats_metadata%icorr_w_hm_1(hm_idx) = k

            call stat_assign( var_index=stats_metadata%icorr_w_hm_1(hm_idx), &
                              var_name="corr_w_"//trim( hm_type(1:2) )//"_1", &
                              var_description="Correlation (in-precip) " &
                              // "of w and " &
                              // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                              // " (1st PDF component) [-]", &
                              var_units="-", l_silhs=.false., grid_kind=stats_zt )

            k = k + 1

            ! The in-precip correlation of w and the hydrometeor in the
            ! 2nd PDF component.
            stats_metadata%icorr_w_hm_2(hm_idx) = k

            call stat_assign( var_index=stats_metadata%icorr_w_hm_2(hm_idx), &
                              var_name="corr_w_"//trim( hm_type(1:2) )//"_2", &
                              var_description="Correlation (in-precip) " &
                              // "of w and " &
                              // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                              // " (2nd PDF component) [-]", &
                              var_units="-", l_silhs=.false., grid_kind=stats_zt )

            k = k + 1

         enddo ! hm_idx = 1, hydromet_dim, 1

      case ( 'corr_w_Ncn_i' )
        stats_metadata%icorr_w_Ncn_1 = k

         call stat_assign( var_index=stats_metadata%icorr_w_Ncn_1, &
                           var_name="corr_w_Ncn_1", &
                           var_description="corr_w_Ncn_1, Correlation of w and N_cn " &
                           // "(1st PDF component)", &
                           var_units="-", l_silhs=.false., grid_kind=stats_zt )

         k = k + 1

         stats_metadata%icorr_w_Ncn_2 = k

         call stat_assign( var_index=stats_metadata%icorr_w_Ncn_2, &
                           var_name="corr_w_Ncn_2", &
                           var_description="corr_w_Ncn_2, Correlation of w and N_cn " &
                           // "(2nd PDF component)", &
                           var_units="-", l_silhs=.false., grid_kind=stats_zt )

         k = k + 1

      case ('corr_chi_eta_1_ca')
        stats_metadata%icorr_chi_eta_1_ca = k

        call stat_assign( var_index=stats_metadata%icorr_chi_eta_1_ca, &
                          var_name="corr_chi_eta_1_ca", &
                          var_description="corr_chi_eta_1_ca, Correlation of chi (s) and" &
                          // " eta (t) (1st PDF component) found in the" &
                          // " correlation array", var_units="-", &
                          l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('corr_chi_eta_2_ca')
        stats_metadata%icorr_chi_eta_2_ca = k

        call stat_assign( var_index=stats_metadata%icorr_chi_eta_2_ca, &
                          var_name="corr_chi_eta_2_ca", &
                          var_description="corr_chi_eta_2_ca, Correlation of chi (s) and" &
                          // " eta (t) (2nd PDF component) found in the" &
                          // " correlation array", var_units="-", &
                          l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      ! Correlation of chi(s) and a hydrometeor (in-precip) for each PDF
      ! component and hydrometeor type.
      case ( 'corr_chi_hm_i' )

         do hm_idx = 1, hydromet_dim, 1

            hm_type = hydromet_list(hm_idx)

            ! The in-precip correlation of chi and the hydrometeor in the
            ! 1st PDF component.
            stats_metadata%icorr_chi_hm_1(hm_idx) = k

            call stat_assign( var_index=stats_metadata%icorr_chi_hm_1(hm_idx), &
                              var_name="corr_chi_"//trim(hm_type(1:2))//"_1", &
                              var_description="Correlation (in-precip) " &
                              // "of chi (s) and " &
                              // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                              // " (1st PDF component) [-]", &
                              var_units="-", l_silhs=.false., grid_kind=stats_zt )

            k = k + 1

            ! The in-precip correlation of chi and the hydrometeor in the
            ! 2nd PDF component.
            stats_metadata%icorr_chi_hm_2(hm_idx) = k

            call stat_assign( var_index=stats_metadata%icorr_chi_hm_2(hm_idx), &
                              var_name="corr_chi_"//trim(hm_type(1:2))//"_2", &
                              var_description="Correlation (in-precip) " &
                              // "of chi (s) and " &
                              // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                              // " (2nd PDF component) [-]", &
                              var_units="-", l_silhs=.false., grid_kind=stats_zt )

            k = k + 1

         enddo ! hm_idx = 1, hydromet_dim, 1

      case ( 'corr_chi_Ncn_i' )
        stats_metadata%icorr_chi_Ncn_1 = k

         call stat_assign( var_index=stats_metadata%icorr_chi_Ncn_1, &
                           var_name="corr_chi_Ncn_1", &
                           var_description="corr_chi_Ncn_1, Correlation of chi and N_cn " &
                           // "(1st PDF component)", &
                           var_units="-", l_silhs=.false., grid_kind=stats_zt )

         k = k + 1

         stats_metadata%icorr_chi_Ncn_2 = k

         call stat_assign( var_index=stats_metadata%icorr_chi_Ncn_2, &
                           var_name="corr_chi_Ncn_2", &
                           var_description="corr_chi_Ncn_2, Correlation of chi and N_cn " &
                           // "(2nd PDF component)", &
                           var_units="-", l_silhs=.false., grid_kind=stats_zt )

         k = k + 1

      ! Correlation of eta(t) and a hydrometeor (in-precip) for each PDF
      ! component and hydrometeor type.
      case ( 'corr_eta_hm_i' )

         do hm_idx = 1, hydromet_dim, 1

            hm_type = hydromet_list(hm_idx)

            ! The in-precip correlation of eta and the hydrometeor in the
            ! 1st PDF component.
            stats_metadata%icorr_eta_hm_1(hm_idx) = k

            call stat_assign( var_index=stats_metadata%icorr_eta_hm_1(hm_idx), &
                              var_name="corr_eta_"//trim(hm_type(1:2))//"_1", &
                              var_description="Correlation (in-precip) " &
                              // "of eta (t) and " &
                              // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                              // " (1st PDF component) [-]", &
                              var_units="-", l_silhs=.false., grid_kind=stats_zt )

            k = k + 1

            ! The in-precip correlation of eta and the hydrometeor in the
            ! 2nd PDF component.
            stats_metadata%icorr_eta_hm_2(hm_idx) = k

            call stat_assign( var_index=stats_metadata%icorr_eta_hm_2(hm_idx), &
                              var_name="corr_eta_"//trim(hm_type(1:2))//"_2", &
                              var_description="Correlation (in-precip) " &
                              // "of eta (t) and " &
                              // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                              // " (2nd PDF component) [-]", &
                              var_units="-", l_silhs=.false., grid_kind=stats_zt )

            k = k + 1

         enddo ! hm_idx = 1, hydromet_dim, 1

      case ( 'corr_eta_Ncn_i' )

         stats_metadata%icorr_eta_Ncn_1 = k

         call stat_assign( var_index=stats_metadata%icorr_eta_Ncn_1, &
                           var_name="corr_eta_Ncn_1", &
                           var_description="corr_eta_Ncn_1, Correlation of eta and N_cn " &
                           // "(1st PDF component)", &
                           var_units="-", l_silhs=.false., grid_kind=stats_zt )

         k = k + 1

         stats_metadata%icorr_eta_Ncn_2 = k

         call stat_assign( var_index=stats_metadata%icorr_eta_Ncn_2, &
                           var_name="corr_eta_Ncn_2", &
                           var_description="corr_eta_Ncn_2, Correlation of eta and N_cn " &
                           // "(2nd PDF component)", &
                           var_units="-", l_silhs=.false., grid_kind=stats_zt )

         k = k + 1

      ! Correlation of Ncn and a hydrometeor (in-precip) for each PDF
      ! component and hydrometeor type.
      case ( 'corr_Ncn_hm_i' )

         do hm_idx = 1, hydromet_dim, 1

            hm_type = hydromet_list(hm_idx)

            ! The in-precip correlation of Ncn and the hydrometeor in the
            ! 1st PDF component.
            stats_metadata%icorr_Ncn_hm_1(hm_idx) = k

            call stat_assign( var_index=stats_metadata%icorr_Ncn_hm_1(hm_idx), &
                              var_name="corr_Ncn_"//trim(hm_type(1:2))//"_1", &
                              var_description="Correlation (in-precip) " &
                              // "of N_cn and " &
                              // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                              // " (1st PDF component) [-]", &
                              var_units="-", l_silhs=.false., grid_kind=stats_zt )

            k = k + 1

            ! The in-precip correlation of Ncn and the hydrometeor in the
            ! 2nd PDF component.
            stats_metadata%icorr_Ncn_hm_2(hm_idx) = k

            call stat_assign( var_index=stats_metadata%icorr_Ncn_hm_2(hm_idx), &
                              var_name="corr_Ncn_"//trim(hm_type(1:2))//"_2", &
                              var_description="Correlation (in-precip) " &
                              // "of N_cn and " &
                              // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                              // " (2nd PDF component) [-]", &
                              var_units="-", l_silhs=.false., grid_kind=stats_zt )

            k = k + 1

         enddo ! hm_idx = 1, hydromet_dim, 1

      ! Correlation (in-precip) of two different hydrometeors (hmx and hmy)
      ! for each PDF component and hydrometeor type.
      case ( 'corr_hmx_hmy_i' )

         do hmx_idx = 1, hydromet_dim, 1

            hmx_type = hydromet_list(hmx_idx)

            do hmy_idx = hmx_idx+1, hydromet_dim, 1

               hmy_type = hydromet_list(hmy_idx)

               ! The in-precip correlation of hmx and hmy in the 1st PDF
               ! component.
               stats_metadata%icorr_hmx_hmy_1(hmy_idx,hmx_idx) = k

               call stat_assign( var_index=stats_metadata%icorr_hmx_hmy_1(hmy_idx,hmx_idx), &
                                 var_name="corr_"//trim( hmx_type(1:2) )//"_" &
                                 // trim( hmy_type(1:2) )//"_1", &
                                 var_description="Correlation (in-precip) " &
                                 // "of " &
                                 // hmx_type(1:1)//"_"//trim( hmx_type(2:2) ) &
                                 // " and " &
                                 // hmy_type(1:1)//"_"//trim( hmy_type(2:2) ) &
                                 // " (1st PDF component) [-]", &
                                 var_units="-", l_silhs=.false., grid_kind=stats_zt )

               k = k + 1

               ! The in-precip correlation of hmx and hmy in the 2nd PDF
               ! component.
               stats_metadata%icorr_hmx_hmy_2(hmy_idx,hmx_idx) = k

               call stat_assign( var_index=stats_metadata%icorr_hmx_hmy_2(hmy_idx,hmx_idx), &
                                 var_name="corr_"//trim( hmx_type(1:2) )//"_" &
                                 // trim( hmy_type(1:2) )//"_2", &
                                 var_description="Correlation (in-precip) " &
                                 // "of " &
                                 // hmx_type(1:1)//"_"//trim( hmx_type(2:2) ) &
                                 // " and " &
                                 // hmy_type(1:1)//"_"//trim( hmy_type(2:2) ) &
                                 // " (2nd PDF component) [-]", &
                                 var_units="-", l_silhs=.false., grid_kind=stats_zt )

               k = k + 1

            enddo ! hmy_idx = hmx_idx+1, hydromet_dim, 1

         enddo ! hmx_idx = 1, hydromet_dim, 1

      ! Correlation (in-precip) of w and ln hm for each PDF component and
      ! hydrometeor type.
      case ( 'corr_w_hm_i_n' )

         do hm_idx = 1, hydromet_dim, 1

            hm_type = hydromet_list(hm_idx)

            ! The in-precip correlation of w and ln hm in the 1st PDF
            ! component.
            stats_metadata%icorr_w_hm_1_n(hm_idx) = k

            call stat_assign( var_index=stats_metadata%icorr_w_hm_1_n(hm_idx), &
                              var_name="corr_w_"//trim(hm_type(1:2))//"_1_n", &
                              var_description="Correlation (in-precip) " &
                              // "of w and ln " &
                              // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                              // " (1st PDF component) [-]", &
                              var_units="-", l_silhs=.false., grid_kind=stats_zt )

            k = k + 1

            ! The in-precip correlation of w and ln hm in the 2nd PDF
            ! component.
            stats_metadata%icorr_w_hm_2_n(hm_idx) = k

            call stat_assign( var_index=stats_metadata%icorr_w_hm_2_n(hm_idx), &
                              var_name="corr_w_"//trim(hm_type(1:2))//"_2_n", &
                              var_description="Correlation (in-precip) " &
                              // "of w and ln " &
                              // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                              // " (2nd PDF component) [-]", &
                              var_units="-", l_silhs=.false., grid_kind=stats_zt )

            k = k + 1

         enddo ! hm_idx = 1, hydromet_dim, 1

      case ( 'corr_w_Ncn_i_n' )

         stats_metadata%icorr_w_Ncn_1_n = k

         call stat_assign( var_index=stats_metadata%icorr_w_Ncn_1_n, &
                           var_name="corr_w_Ncn_1_n", &
                           var_description="corr_w_Ncn_1_n, Correlation of w and " &
                           // "ln N_cn (1st PDF component)", &
                           var_units="-", l_silhs=.false., grid_kind=stats_zt )

         k = k + 1

         stats_metadata%icorr_w_Ncn_2_n = k

         call stat_assign( var_index=stats_metadata%icorr_w_Ncn_2_n, &
                           var_name="corr_w_Ncn_2_n", &
                           var_description="corr_w_Ncn_2_n, Correlation of w and " &
                           // "ln N_cn (2nd PDF component)", &
                           var_units="-", l_silhs=.false., grid_kind=stats_zt )

         k = k + 1

      ! Correlation (in-precip) of chi and ln hm for each PDF component and
      ! hydrometeor type.
      case ( 'corr_chi_hm_i_n' )

         do hm_idx = 1, hydromet_dim, 1

            hm_type = hydromet_list(hm_idx)

            ! The in-precip correlation of chi and ln hm in the 1st PDF
            ! component.
            stats_metadata%icorr_chi_hm_1_n(hm_idx) = k

            call stat_assign( var_index=stats_metadata%icorr_chi_hm_1_n(hm_idx), &
                              var_name="corr_chi_"//trim(hm_type(1:2)) &
                              // "_1_n", &
                              var_description="Correlation (in-precip) " &
                              // "of chi (s) and ln " &
                              // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                              // " (1st PDF component) [-]", &
                              var_units="-", l_silhs=.false., grid_kind=stats_zt )

            k = k + 1

            ! The in-precip correlation of chi(s) and ln hm in the 2nd PDF
            ! component.
            stats_metadata%icorr_chi_hm_2_n(hm_idx) = k

            call stat_assign( var_index=stats_metadata%icorr_chi_hm_2_n(hm_idx), &
                              var_name="corr_chi_"//trim(hm_type(1:2)) &
                              // "_2_n", &
                              var_description="Correlation (in-precip) " &
                              // "of chi (s) and ln " &
                              // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                              // " (2nd PDF component) [-]", &
                              var_units="-", l_silhs=.false., grid_kind=stats_zt )

            k = k + 1

         enddo ! hm_idx = 1, hydromet_dim, 1

      case ( 'corr_chi_Ncn_i_n' )

         stats_metadata%icorr_chi_Ncn_1_n = k

        call stat_assign( var_index=stats_metadata%icorr_chi_Ncn_1_n, &
                           var_name="corr_chi_Ncn_1_n", &
                           var_description="corr_chi_Ncn_1_n, Correlation of chi (s) and " &
                           // "ln N_cn (1st PDF component)", &
                           var_units="-", l_silhs=.false., grid_kind=stats_zt )

         k = k + 1

         stats_metadata%icorr_chi_Ncn_2_n = k

         call stat_assign( var_index=stats_metadata%icorr_chi_Ncn_2_n, &
                           var_name="corr_chi_Ncn_2_n", &
                           var_description="corr_chi_Ncn_2_n, Correlation of chi (s) and " &
                           // "ln N_cn (2nd PDF component)", &
                           var_units="-", l_silhs=.false., grid_kind=stats_zt )

         k = k + 1

      ! Correlation (in-precip) of eta and ln hm for each PDF component and
      ! hydrometeor type.
      case ( 'corr_eta_hm_i_n' )

         do hm_idx = 1, hydromet_dim, 1

            hm_type = hydromet_list(hm_idx)

            ! The in-precip correlation of eta and ln hm in the 1st PDF
            ! component.
            stats_metadata%icorr_eta_hm_1_n(hm_idx) = k

            call stat_assign( var_index=stats_metadata%icorr_eta_hm_1_n(hm_idx), &
                              var_name="corr_eta_"//trim( hm_type(1:2) ) &
                              // "_1_n", &
                              var_description="Correlation (in-precip) " &
                              // "of eta (t) and ln " &
                              // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                              // " (1st PDF component) [-]", &
                              var_units="-", l_silhs=.false., grid_kind=stats_zt )

            k = k + 1

            ! The in-precip correlation of eta and ln hm in the 2nd PDF
            ! component.
            stats_metadata%icorr_eta_hm_2_n(hm_idx) = k

            call stat_assign( var_index=stats_metadata%icorr_eta_hm_2_n(hm_idx), &
                              var_name="corr_eta_"//trim( hm_type(1:2) ) &
                              // "_2_n", &
                              var_description="Correlation (in-precip) " &
                              // "of eta(t) and ln " &
                              // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                              // " (2nd PDF component) [-]", &
                              var_units="-", l_silhs=.false., grid_kind=stats_zt )

            k = k + 1

         enddo ! hm_idx = 1, hydromet_dim, 1

      case ( 'corr_eta_Ncn_i_n' )

         stats_metadata%icorr_eta_Ncn_1_n = k

         call stat_assign( var_index=stats_metadata%icorr_eta_Ncn_1_n, &
                           var_name="corr_eta_Ncn_1_n", &
                           var_description="corr_eta_Ncn_1_n, Correlation of eta (t) and " &
                           // "ln N_cn (1st PDF component)", &
                           var_units="-", l_silhs=.false., grid_kind=stats_zt )

         k = k + 1

         stats_metadata%icorr_eta_Ncn_2_n = k

         call stat_assign( var_index=stats_metadata%icorr_eta_Ncn_2_n, &
                           var_name="corr_eta_Ncn_2_n", &
                           var_description="corr_eta_Ncn_2_n, Correlation of eta (t) and " &
                           // "ln N_cn (2nd PDF component)", &
                           var_units="-", l_silhs=.false., grid_kind=stats_zt )

         k = k + 1

      ! Correlation (in-precip) of ln Ncn and ln hm for each PDF component
      ! and hydrometeor type.
      case ( 'corr_Ncn_hm_i_n' )

         do hm_idx = 1, hydromet_dim, 1

            hm_type = hydromet_list(hm_idx)

            ! The in-precip correlation of ln Ncn and ln hm in the 1st PDF
            ! component.
            stats_metadata%icorr_Ncn_hm_1_n(hm_idx) = k

            call stat_assign( var_index=stats_metadata%icorr_Ncn_hm_1_n(hm_idx), &
                              var_name="corr_Ncn_"//trim(hm_type(1:2)) &
                              // "_1_n", &
                              var_description="Correlation (in-precip) " &
                              // "of ln N_cn and ln " &
                              // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                              // " (1st PDF component) [-]", &
                              var_units="-", l_silhs=.false., grid_kind=stats_zt )

            k = k + 1

            ! The in-precip correlation of ln Ncn and ln hm in the 2nd PDF
            ! component.
            stats_metadata%icorr_Ncn_hm_2_n(hm_idx) = k

            call stat_assign( var_index=stats_metadata%icorr_Ncn_hm_2_n(hm_idx), &
                              var_name="corr_Ncn_"//trim(hm_type(1:2)) &
                              // "_2_n", &
                              var_description="Correlation (in-precip) " &
                              // "of ln N_cn and ln " &
                              // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                              // " (2nd PDF component) [-]", &
                              var_units="-", l_silhs=.false., grid_kind=stats_zt )

            k = k + 1

         enddo ! hm_idx = 1, hydromet_dim, 1

      ! Correlation (in-precip) of ln hmx and ln hmy (hmx and hmy are two
      ! different hydrometeors) for each PDF component and hydrometeor type.
      case ( 'corr_hmx_hmy_i_n' )

         do hmx_idx = 1, hydromet_dim, 1

            hmx_type = hydromet_list(hmx_idx)

            do hmy_idx = hmx_idx+1, hydromet_dim, 1

               hmy_type = hydromet_list(hmy_idx)

               ! The in-precip correlation of ln hmx and ln hmy in the 1st
               ! PDF component.
               stats_metadata%icorr_hmx_hmy_1_n(hmy_idx,hmx_idx) = k

               call stat_assign( var_index=stats_metadata%icorr_hmx_hmy_1_n(hmy_idx,hmx_idx), &
                                 var_name="corr_"//trim( hmx_type(1:2) )//"_" &
                                 // trim( hmy_type(1:2) )//"_1_n", &
                                 var_description="Correlation (in-precip) " &
                                 // "of ln " &
                                 // hmx_type(1:1)//"_"//trim( hmx_type(2:2) ) &
                                 // " and ln " &
                                 // hmy_type(1:1)//"_"//trim( hmy_type(2:2) ) &
                                 // " (1st PDF component) [-]", &
                                 var_units="-", l_silhs=.false., grid_kind=stats_zt )

               k = k + 1

               ! The in-precip correlation of ln hmx and ln hmy in the 2nd
               ! PDF component.
               stats_metadata%icorr_hmx_hmy_2_n(hmy_idx,hmx_idx) = k

               call stat_assign( var_index=stats_metadata%icorr_hmx_hmy_2_n(hmy_idx,hmx_idx), &
                                 var_name="corr_"//trim( hmx_type(1:2) )//"_" &
                                 // trim( hmy_type(1:2) )//"_2_n", &
                                 var_description="Correlation (in-precip) " &
                                 // "of ln " &
                                 // hmx_type(1:1)//"_"//trim( hmx_type(2:2) ) &
                                 // " and ln " &
                                 // hmy_type(1:1)//"_"//trim( hmy_type(2:2) ) &
                                 // " (2nd PDF component) [-]", &
                                 var_units="-", l_silhs=.false., grid_kind=stats_zt )

               k = k + 1

            enddo ! hmy_idx = hmx_idx+1, hydromet_dim, 1

         enddo ! hmx_idx = 1, hydromet_dim, 1

      ! Third-order mixed moment < w'^2 hm' >, where hm is a hydrometeor.
      case ('wp2hmp')

         do hm_idx = 1, hydromet_dim, 1

            hm_type = hydromet_list(hm_idx)

            stats_metadata%iwp2hmp(hm_idx) = k

            if ( l_mix_rat_hm(hm_idx) ) then

               call stat_assign( var_index=stats_metadata%iwp2hmp(hm_idx), &
                                 var_name="wp2"//trim( hm_type(1:2) )//"p", &
                                 var_description="Third-order moment < w'^2 " &
                                 // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                                 // "' > [(m/s)^2 kg/kg]", &
                                 var_units="(m/s)^2 kg/kg", &
                                 l_silhs=.false., grid_kind=stats_zt )

            else ! Concentration

               call stat_assign( var_index=stats_metadata%iwp2hmp(hm_idx), &
                                 var_name="wp2"//trim( hm_type(1:2) )//"p", &
                                 var_description="Third-order moment < w'^2 " &
                                 // hm_type(1:1)//"_"//trim( hm_type(2:2) ) &
                                 // "' > [(m/s)^2 num/kg]", &
                                 var_units="(m/s)^2 num/kg", &
                                 l_silhs=.false., grid_kind=stats_zt )

            endif ! l_mix_rat_hm(hm_idx)

            k = k + 1

         enddo ! hm_idx = 1, hydromet_dim, 1

      case ('cloud_frac_refined')
        stats_metadata%icloud_frac_refined = k
        call stat_assign( var_index=stats_metadata%icloud_frac_refined, var_name="cloud_frac_refined", &
                          var_description="cloud_frac_refined, Cloud fraction computed on " &
                          // "refined grid", &
                          var_units="-", l_silhs=.false., grid_kind=stats_zt )
        k = k + 1

      case ('rcm_refined')
        stats_metadata%ircm_refined = k
        call stat_assign( var_index=stats_metadata%ircm_refined, var_name="rcm_refined", &
                          var_description="rcm_refined, Cloud water mixing ratio computed on " &
                          // "refined grid &
                          &[kg/kg]", var_units="kg/kg", l_silhs=.false., grid_kind=stats_zt)
        k = k + 1

      case ('hl_on_Cp_residual')
        stats_metadata%ihl_on_Cp_residual = k
        call stat_assign( var_index=stats_metadata%ihl_on_Cp_residual, var_name="hl_on_Cp_residual", &
                          var_description="hl_on_Cp_residual, Residual change in HL/Cp from " &
                          // "Morrison microphysics &
                          &not due to sedimentation", &
                          var_units="K", l_silhs=.true., grid_kind=stats_zt)
        k = k + 1

      case ('qto_residual')
        stats_metadata%iqto_residual = k
        call stat_assign( var_index=stats_metadata%iqto_residual, var_name="qto_residual", &
                          var_description="qto_residual, Residual change in total water " &
                          // "from Morrison &
                          &microphysics not due to sedimentation", &
                          var_units="kg/kg", l_silhs=.true., grid_kind=stats_zt)
        k = k + 1

      case ( 'sclrm' )
        do j = 1, sclr_dim, 1
          write(sclr_idx, * ) j
          sclr_idx = adjustl(sclr_idx)
        stats_metadata%isclrm(j) = k
          call stat_assign( var_index=stats_metadata%isclrm(j), var_name="sclr"//trim(sclr_idx)//"m", &
            var_description="passive scalar "//trim(sclr_idx), var_units="unknown", &
            l_silhs=.false., grid_kind=stats_zt )
          k = k + 1
        end do

      case ( 'sclrm_f' )
        do j = 1, sclr_dim, 1
          write(sclr_idx, * ) j
          sclr_idx = adjustl(sclr_idx)
        stats_metadata%isclrm_f(j) = k
          call stat_assign( var_index=stats_metadata%isclrm_f(j), var_name="sclr"//trim(sclr_idx)//"m_f", &
            var_description="passive scalar forcing "//trim(sclr_idx), var_units="unknown", &
            l_silhs=.false., grid_kind=stats_zt )
          k = k + 1
        end do

      case ( 'edsclrm' )
        do j = 1, edsclr_dim, 1
          write(sclr_idx, * ) j
          sclr_idx = adjustl(sclr_idx)
        stats_metadata%iedsclrm(j) = k
          call stat_assign( var_index=stats_metadata%iedsclrm(j), var_name="edsclr"//trim(sclr_idx)//"m", &
            var_description="passive scalar "//trim(sclr_idx), var_units="unknown", &
            l_silhs=.false., grid_kind=stats_zt )
          k = k + 1
        end do

      case ( 'edsclrm_f' )
        do j = 1, edsclr_dim, 1
          write(sclr_idx, * ) j
          sclr_idx = adjustl(sclr_idx)
        stats_metadata%iedsclrm_f(j) = k
          call stat_assign( var_index=stats_metadata%iedsclrm_f(j), var_name="edsclr"//trim(sclr_idx)//"m_f", &
            var_description="passive scalar forcing "//trim(sclr_idx), var_units="unknown", &
            l_silhs=.false., grid_kind=stats_zt )
          k = k + 1
        end do

      case default

        write(fstderr,*) 'Error:  unrecognized variable in vars_zt:  ', trim( vars_zt(i) )
        l_error = .true.  ! This will stop the run.

      end select ! trim( vars_zt )


    end do ! i=1,stats_zt%num_output_fields


    return

  end subroutine stats_init_zt

!===============================================================================

end module stats_zt_module
