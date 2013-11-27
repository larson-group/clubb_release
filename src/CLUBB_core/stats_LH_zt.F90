!-----------------------------------------------------------------------
! $Id$

module stats_LH_zt

  implicit none

  private ! Default Scope

  public :: stats_init_LH_zt

! Constant parameters
  integer, parameter, public :: nvarmax_LH_zt = 100 ! Maximum variables allowed

  contains

!-----------------------------------------------------------------------
  subroutine stats_init_LH_zt( vars_LH_zt, l_error )

! Description:
!   Initializes array indices for zt

! Note:
!   All code that is within subroutine stats_init_zt, including variable
!   allocation code, is not called if l_stats is false.  This subroutine is
!   called only when l_stats is true.

!-----------------------------------------------------------------------

    use constants_clubb, only:  &
      fstderr ! Constant(s)

    use stats_variables, only: & 
      LH_zt    ! Variable

    use stats_variables, only: & 
      iAKm, &  ! Variable(s)
      iLH_AKm, & 
      iAKstd, & 
      iAKstd_cld, & 
      iAKm_rcm, & 
      iAKm_rcc

    use stats_variables, only: &
      iLH_thlm_mc, &  ! Variable(s)
      iLH_rvm_mc, & 
      iLH_rcm_mc, & 
      iLH_Ncm_mc, & 
      iLH_rrainm_mc, & 
      iLH_Nrm_mc, & 
      iLH_rsnowm_mc, & 
      iLH_Nsnowm_mc, & 
      iLH_rgraupelm_mc, & 
      iLH_Ngraupelm_mc, & 
      iLH_ricem_mc, & 
      iLH_Nim_mc, & 
      iLH_Vrr, &
      iLH_VNr, &
      iLH_rcm_avg

    use stats_variables, only: &
      iLH_rrainm, & ! Variable(s)
      iLH_Nrm, &
      iLH_ricem, &
      iLH_Nim, &
      iLH_rsnowm, &
      iLH_Nsnowm, &
      iLH_rgraupelm, &
      iLH_Ngraupelm, &
      iLH_thlm, &
      iLH_rcm, &
      iLH_Ncm, &
      iLH_rvm, &
      iLH_wm, &
      iLH_wp2_zt, &
      iLH_rcp2_zt, &
      iLH_rtp2_zt, &
      iLH_thlp2_zt, &
      iLH_rrainp2_zt, &
      iLH_Nrp2_zt, &
      iLH_Ncp2_zt, &
      iLH_cloud_frac, &
      iLH_s_mellor, &
      iLH_t_mellor, &
      iLH_sp2, &
      iLH_rrainm_auto, &
      iLH_rrainm_accr, &
      iLH_rrainm_evap, &
      iLH_Nrm_auto, &
      iLH_Nrm_cond

    use stats_variables, only: &
      iLH_rrainm_src_adj,  & ! Variable(s)
      iLH_rrainm_cond_adj, &
      iLH_Nrm_src_adj,     &
      iLH_Nrm_cond_adj

    use stats_variables, only: &
      iLH_precip_frac, &
      iLH_mixt_frac

    use stats_type, only: & 
      stat_assign ! Procedure

    implicit none

    ! External
    intrinsic :: trim

    ! Input Variable
    character(len= * ), dimension(nvarmax_LH_zt), intent(in) :: vars_LH_zt

    ! Input / Output Variable        
    logical, intent(inout) :: l_error

    ! Local Varables
    integer :: i, k

    ! ---- Begin Code ----

    ! Default initialization for array indices for LH_zt is zero (see module
    ! stats_variables)

    ! Assign pointers for statistics variables zt

    k = 1
    do i = 1, LH_zt%nn

      select case ( trim( vars_LH_zt(i) ) )
      case ( 'AKm' )           ! Vince Larson 22 May 2005
        iAKm = k
        call stat_assign( var_index=iAKm, var_name="AKm", &
             var_description="Analytic Kessler ac [kg/kg]", var_units="kg/kg", l_silhs=.true., &
             grid_kind=LH_zt )
        k = k + 1

      case ( 'LH_AKm' )       ! Vince Larson 22 May 2005
        iLH_AKm = k

        call stat_assign( var_index=iLH_AKm, var_name="LH_AKm", &
             var_description="LH Kessler estimate  [kg/kg/s]", var_units="kg/kg/s", &
             l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case ( 'AKstd' )
        iAKstd = k

        call stat_assign( var_index=iAKstd, var_name="AKstd", &
             var_description="Exact standard deviation of gba Kessler [kg/kg/s]", &
             var_units="kg/kg/s", l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case ( 'AKstd_cld' )
        iAKstd_cld = k

        call stat_assign( var_index=iAKstd_cld, var_name="AKstd_cld", &
             var_description="Exact w/in cloud std of gba Kessler [kg/kg/s]", &
             var_units="kg/kg/s", l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case ( 'AKm_rcm' )
        iAKm_rcm = k

        call stat_assign( var_index=iAKm_rcm, var_name="AKm_rcm", &
             var_description="Exact local gba auto based on rcm [kg/kg/s]", var_units="kg/kg/s", &
             l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case ( 'AKm_rcc' )
        iAKm_rcc = k

        call stat_assign( var_index=iAKm_rcc, var_name="AKm_rcc", &
             var_description="Exact local gba based on w/in cloud rc [kg/kg/s]", &
             var_units="kg/kg/s", l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case ( 'LH_rvm_mc' )
        iLH_rvm_mc = k

        call stat_assign( var_index=iLH_rvm_mc, var_name="LH_rvm_mc", &
             var_description="Latin hypercube estimate of rvm_mc [kg/kg/s]", var_units="kg/kg/s", &
             l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case ( 'LH_thlm_mc' )
        iLH_thlm_mc = k

        call stat_assign( var_index=iLH_thlm_mc, var_name="LH_thlm_mc", &
             var_description="Latin hypercube estimate of thlm_mc [kg/kg/s]", &
             var_units="kg/kg/s", l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case ( 'LH_rcm_mc' )
        iLH_rcm_mc = k

        call stat_assign( var_index=iLH_rcm_mc, var_name="LH_rcm_mc", &
             var_description="Latin hypercube estimate of rcm_mc [kg/kg/s]", var_units="kg/kg/s", &
             l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case ( 'LH_Ncm_mc' )
        iLH_Ncm_mc = k

        call stat_assign( var_index=iLH_Ncm_mc, var_name="LH_Ncm_mc", &
             var_description="Latin hypercube estimate of Ncm_mc [kg/kg/s]", var_units="kg/kg/s", &
             l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case ( 'LH_rrainm_mc' )
        iLH_rrainm_mc = k

        call stat_assign( var_index=iLH_rrainm_mc, var_name="LH_rrainm_mc", &
             var_description="Latin hypercube estimate of rrainm_mc [kg/kg/s]", &
             var_units="kg/kg/s", l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case ( 'LH_Nrm_mc' )
        iLH_Nrm_mc = k

        call stat_assign( var_index=iLH_Nrm_mc, var_name="LH_Nrm_mc", &
             var_description="Latin hypercube estimate of Nrm_mc [kg/kg/s]", var_units="kg/kg/s", &
             l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case('LH_rsnowm_mc')
        iLH_rsnowm_mc = k

        call stat_assign( var_index=iLH_rsnowm_mc, var_name="LH_rsnowm_mc", &
             var_description="Latin hypercube estimate of rsnowm_mc [kg/kg/s]", &
             var_units="kg/kg/s", l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case ( 'LH_Nsnowm_mc' )
        iLH_Nsnowm_mc = k

        call stat_assign( var_index=iLH_Nsnowm_mc, var_name="LH_Nsnowm_mc", &
             var_description="Latin hypercube estimate of Nsnowm_mc [kg/kg/s]", &
             var_units="kg/kg/s", l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case ( 'LH_rgraupelm_mc' )
        iLH_rgraupelm_mc = k

        call stat_assign( var_index=iLH_rgraupelm_mc, var_name="LH_rgraupelm_mc", &
             var_description="Latin hypercube estimate of rgraupelm_mc [kg/kg/s]", &
             var_units="kg/kg/s", l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case ( 'LH_Ngraupelm_mc' )
        iLH_Ngraupelm_mc = k

        call stat_assign( var_index=iLH_Ngraupelm_mc, var_name="LH_Ngraupelm_mc", &
             var_description="Latin hypercube estimate of Ngraupelm_mc [kg/kg/s]", &
             var_units="kg/kg/s", l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case ( 'LH_ricem_mc' )
        iLH_ricem_mc = k

        call stat_assign( var_index=iLH_ricem_mc, var_name="LH_ricem_mc", &
             var_description="Latin hypercube estimate of ricem_mc [kg/kg/s]", &
             var_units="kg/kg/s", l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case ( 'LH_Nim_mc' )
        iLH_Nim_mc = k

        call stat_assign( var_index=iLH_Nim_mc, var_name="LH_Nim_mc", &
             var_description="Latin hypercube estimate of Nim_mc [kg/kg/s]", var_units="kg/kg/s", &
             l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case ( 'LH_Vrr' )
        iLH_Vrr = k

        call stat_assign( var_index=iLH_Vrr, var_name="LH_Vrr", &
             var_description="Latin hypercube estimate of rrainm sedimentation velocity [m/s]", &
             var_units="m/s", l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case ( 'LH_VNr' )
        iLH_VNr = k

        call stat_assign( var_index=iLH_VNr, var_name="LH_VNr", &
             var_description="Latin hypercube estimate of Nrm sedimentation velocity [m/s]", &
             var_units="m/s", l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case ( 'LH_rcm_avg' )
        iLH_rcm_avg = k

        call stat_assign( var_index=iLH_rcm_avg, var_name="LH_rcm_avg", &
             var_description="Latin hypercube average estimate of rcm [kg/kg]", &
             var_units="kg/kg", l_silhs=.true., grid_kind=LH_zt )

        k = k + 1

      case ( 'LH_rrainm' )
        iLH_rrainm = k

        call stat_assign( var_index=iLH_rrainm, var_name="LH_rrainm", &
             var_description="Latin hypercube estimate of rrainm [kg/kg]", var_units="kg/kg", &
             l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case ( 'LH_Nrm' )
        iLH_Nrm = k

        call stat_assign( var_index=iLH_Nrm, var_name="LH_Nrm", &
             var_description="Latin hypercube estimate of Nrm [count/kg]", var_units="count/kg", &
             l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case ( 'LH_ricem' )
        iLH_ricem = k

        call stat_assign( var_index=iLH_ricem, var_name="LH_ricem", &
             var_description="Latin hypercube estimate of ricem [kg/kg]", var_units="kg/kg", &
             l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case ( 'LH_Nim' )
        iLH_Nim = k

        call stat_assign( var_index=iLH_Nim, var_name="LH_Nim", &
             var_description="Latin hypercube estimate of Nim [count/kg]", var_units="count/kg", &
             l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case ( 'LH_rsnowm' )
        iLH_rsnowm = k

        call stat_assign( var_index=iLH_rsnowm, var_name="LH_rsnowm", &
             var_description="Latin hypercube estimate of rsnowm [kg/kg]", var_units="kg/kg", &
             l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case ( 'LH_Nsnowm' )
        iLH_Nsnowm = k

        call stat_assign( var_index=iLH_Nsnowm, var_name="LH_Nsnowm", &
             var_description="Latin hypercube estimate of Nsnowm [count/kg]", &
             var_units="count/kg", l_silhs=.true., grid_kind=LH_zt )
        k = k + 1


      case ( 'LH_rgraupelm' )
        iLH_rgraupelm = k

        call stat_assign( var_index=iLH_rgraupelm, var_name="LH_rgraupelm", &
             var_description="Latin hypercube estimate of rgraupelm [kg/kg]", var_units="kg/kg", &
             l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case ( 'LH_Ngraupelm' )
        iLH_Ngraupelm = k

        call stat_assign( var_index=iLH_Ngraupelm, var_name="LH_Ngraupelm", &
             var_description="Latin hypercube estimate of Ngraupelm [kg/kg]", var_units="kg/kg", &
             l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case ( 'LH_thlm' )
        iLH_thlm = k

        call stat_assign( var_index=iLH_thlm, var_name="LH_thlm", &
             var_description="Latin hypercube estimate of thlm [K]", var_units="K", &
             l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case ( 'LH_rcm' )
        iLH_rcm = k

        call stat_assign( var_index=iLH_rcm, var_name="LH_rcm", &
             var_description="Latin hypercube estimate of rcm [kg/kg]", var_units="kg/kg", &
             l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case ( 'LH_Ncm' )
        iLH_Ncm = k

        call stat_assign( var_index=iLH_Ncm, var_name="LH_Ncm", &
             var_description="Latin hypercube estimate of Ncm [count/kg]", var_units="count/kg", &
             l_silhs=.true., grid_kind=LH_zt )
        k = k + 1


      case ( 'LH_rvm' )
        iLH_rvm = k

        call stat_assign( var_index=iLH_rvm, var_name="LH_rvm", &
             var_description="Latin hypercube estimate of rvm [kg/kg]", var_units="kg/kg", &
             l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case ( 'LH_wm' )
        iLH_wm = k

        call stat_assign( var_index=iLH_wm, var_name="LH_wm", &
             var_description="Latin hypercube estimate of vertical velocity [m/s]", &
             var_units="m/s", l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case ( 'LH_cloud_frac' )
        iLH_cloud_frac = k

        ! Note: count is the udunits compatible unit
        call stat_assign( var_index=iLH_cloud_frac, var_name="LH_cloud_frac", &
             var_description="Latin hypercube estimate of cloud fraction [count]", &
             var_units="count", l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case ( 'LH_s_mellor' )
        iLH_s_mellor = k
        call stat_assign( var_index=iLH_s_mellor, var_name="LH_s_mellor", &
             var_description="Latin hypercube estimate of Mellor's s (extended liq) [kg/kg]", &
             var_units="kg/kg", l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case ( 'LH_t_mellor' )
        iLH_t_mellor = k
        call stat_assign( var_index=iLH_t_mellor, var_name="LH_t_mellor", &
             var_description="Latin hypercube estimate of Mellor's t [kg/kg]", var_units="kg/kg", &
             l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case ( 'LH_sp2' )
        iLH_sp2 = k
        call stat_assign( var_index=iLH_sp2, var_name="LH_sp2", &
             var_description="Latin hypercube estimate of variance of s [kg/kg]", &
             var_units="kg/kg", l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case ( 'LH_wp2_zt' )
        iLH_wp2_zt = k
        call stat_assign( var_index=iLH_wp2_zt, var_name="LH_wp2_zt", &
             var_description="Variance of the latin hypercube estimate of w [m^2/s^2]", &
             var_units="m^2/s^2", l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case ( 'LH_Ncp2_zt' )
        iLH_Ncp2_zt = k
        call stat_assign( var_index=iLH_Ncp2_zt, var_name="LH_Ncp2_zt", &
             var_description="Variance of the latin hypercube estimate of Nc [count^2/kg^2]", &
             var_units="count^2/kg^2", l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case ( 'LH_Nrp2_zt' )
        iLH_Nrp2_zt = k
        call stat_assign( var_index=iLH_Nrp2_zt, var_name="LH_Nrp2_zt", &
             var_description="Variance of the latin hypercube estimate of Nr [count^2/kg^2]", &
             var_units="count^2/kg^2", l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case ( 'LH_rcp2_zt' )
        iLH_rcp2_zt = k
        call stat_assign( var_index=iLH_rcp2_zt, var_name="LH_rcp2_zt", &
             var_description="Variance of the latin hypercube estimate of rc [kg^2/kg^2]", &
             var_units="kg^2/kg^2", l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case ( 'LH_rtp2_zt' )
        iLH_rtp2_zt = k
        call stat_assign( var_index=iLH_rtp2_zt, var_name="LH_rtp2_zt", &
             var_description="Variance of the latin hypercube estimate of rt [kg^2/kg^2]", &
             var_units="kg^2/kg^2", l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case ( 'LH_thlp2_zt' )
        iLH_thlp2_zt = k
        call stat_assign( var_index=iLH_thlp2_zt, var_name="LH_thlp2_zt", &
             var_description="Variance of the latin hypercube estimate of thl [K^2]", &
             var_units="K^2", l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case ( 'LH_rrainp2_zt' )
        iLH_rrainp2_zt = k
        call stat_assign( var_index=iLH_rrainp2_zt, var_name="LH_rrainp2_zt", &
             var_description="Variance of the latin hypercube estimate of rrain [kg^2/kg^2]", &
             var_units="kg^2/kg^2", l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case ( 'LH_rrainm_auto' )
        iLH_rrainm_auto = k
        call stat_assign( var_index=iLH_rrainm_auto, var_name="LH_rrainm_auto", &
             var_description="Latin hypercube estimate of autoconversion [kg/kg/s]", &
             var_units="kg/kg/s", l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case ( 'LH_rrainm_accr' )
        iLH_rrainm_accr = k
        call stat_assign( var_index=iLH_rrainm_accr, var_name="LH_rrainm_accr", &
             var_description="Latin hypercube estimate of accretion [kg/kg/s]", &
             var_units="kg/kg/s", l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case ( 'LH_rrainm_evap' )
        iLH_rrainm_evap = k
        call stat_assign( var_index=iLH_rrainm_evap, var_name="LH_rrainm_evap", &
             var_description="Latin hypercube estimate of evaporation [kg/kg/s]", &
             var_units="kg/kg/s", l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case ( 'LH_Nrm_auto' )
        iLH_Nrm_auto = k
        call stat_assign( var_index=iLH_Nrm_auto, var_name="LH_Nrm_auto", &
             var_description="Latin hypercube estimate of Nrm autoconversion [num/kg/s]", &
             var_units="num/kg/s", l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case ( 'LH_Nrm_cond' )
        iLH_Nrm_cond = k
        call stat_assign( var_index=iLH_Nrm_cond, var_name="LH_Nrm_cond", &
             var_description="Latin hypercube estimate of Nrm evaporation [num/kg/s]", &
             var_units="num/kg/s", l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case ( 'LH_rrainm_src_adj' )
        iLH_rrainm_src_adj = k
        call stat_assign( var_index=iLH_rrainm_src_adj, var_name="LH_rrainm_src_adj", &
             var_description="Latin hypercube estimate of source adjustment (KK only!) [kg/kg/s]", &
             var_units="kg/kg/s", l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case ( 'LH_rrainm_cond_adj' )
        iLH_rrainm_cond_adj = k
        call stat_assign( var_index=iLH_rrainm_cond_adj, var_name="LH_rrainm_cond_adj", &
             var_description="Latin hypercube estimate of evap adjustment (KK only!) [kg/kg/s]", &
             var_units="kg/kg/s", l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case ( 'LH_Nrm_src_adj' )
        iLH_Nrm_src_adj = k
        call stat_assign( var_index=iLH_Nrm_src_adj, var_name="LH_Nrm_src_adj", &
             var_description="Latin hypercube estimate of Nrm source adjustment (KK only!) &
             &[kg/kg/s]", &
             var_units="kg/kg/s", l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case ( 'LH_Nrm_cond_adj' )
        iLH_Nrm_cond_adj = k
        call stat_assign( var_index=iLH_Nrm_cond_adj, var_name="LH_Nrm_cond_adj", &
             var_description="Latin hypercube estimate of Nrm evap adjustment (KK only!) &
             &[kg/kg/s]", &
             var_units="kg/kg/s", l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case ( 'LH_precip_frac' )
        iLH_precip_frac = k
        call stat_assign( var_index=iLH_precip_frac, var_name="LH_precip_frac", &
             var_description="Latin hypercube estimate of precipitation fraction [-]", &
             var_units="-", l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case ( 'LH_mixt_frac' )
        iLH_mixt_frac = k
        call stat_assign( var_index=iLH_mixt_frac, var_name="LH_mixt_frac", &
             var_description="Latin hypercube estimate of mixture fraction (weight of 1st PDF &
             &component [-]", &
             var_units="-", l_silhs=.true., grid_kind=LH_zt )
        k = k + 1

      case default

        write(fstderr,*) 'Error:  unrecognized variable in vars_LH_zt:  ', trim( vars_LH_zt(i) )

        l_error = .true.  ! This will stop the run.

      end select

    end do ! i = 1, LH_zt%nn

    return
  end subroutine stats_init_LH_zt

end module stats_LH_zt
