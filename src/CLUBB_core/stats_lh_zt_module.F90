!-----------------------------------------------------------------------
! $Id$
!===============================================================================
module stats_lh_zt_module

  implicit none

  private ! Default Scope

  public :: stats_init_lh_zt

! Constant parameters
  integer, parameter, public :: nvarmax_lh_zt = 100 ! Maximum variables allowed

  contains

!-----------------------------------------------------------------------
  subroutine stats_init_lh_zt( vars_lh_zt,                    & ! In
                               l_error,                       & ! In/Out
                               stats_metadata, stats_lh_zt )    ! In/Out

! Description:
!   Initializes array indices for stats_zt

! Note:
!   All code that is within subroutine stats_init_zt, including variable
!   allocation code, is not called if l_stats is false.  This subroutine is
!   called only when l_stats is true.

!-----------------------------------------------------------------------

    use constants_clubb, only:  &
        fstderr ! Constant(s)

    use stats_type_utilities, only: & 
        stat_assign ! Procedure

    use stats_type, only: &
        stats ! Type

    use stats_variables, only: &
        stats_metadata_type

    implicit none

    type (stats), target, intent(inout) :: &
      stats_lh_zt

    ! External
    intrinsic :: trim

    ! Local Constants
    integer, parameter :: &
      silhs_num_importance_categories = 8

    !--------------------- Input Variable ---------------------
    character(len= * ), dimension(nvarmax_lh_zt), intent(in) :: &
      vars_lh_zt

    !--------------------- InOut Variables ---------------------      
    type (stats_metadata_type), intent(inout) :: &
      stats_metadata

    logical, intent(inout) :: l_error

    !--------------------- Local Varables ---------------------
    integer :: i, k, tot_loops, icategory

    character( len = 1 ) :: category_num_as_string

    !--------------------- Begin Code ---------------------

    ! Default initialization for array indices for stats_lh_zt is zero (see module
    ! stats_variables)

    ! If any of the index arrays are allocated, then we have called this before
    ! to set up stats_metadata, so all we want to do is set stats_lh_zt via stats_assign
    if ( .not. allocated(stats_metadata%isilhs_variance_category) ) then
      allocate( stats_metadata%isilhs_variance_category(silhs_num_importance_categories), &
                stats_metadata%ilh_samp_frac_category(silhs_num_importance_categories) )
    end if
    
    stats_metadata%isilhs_variance_category(:) = 0
    stats_metadata%ilh_samp_frac_category(:) = 0

    ! Assign pointers for statistics variables stats_zt

    tot_loops = stats_lh_zt%num_output_fields

    if ( any( vars_lh_zt == "silhs_variance_category" ) ) then
       ! Correct for number of variables found under "silhs_variance_category".
       ! Subtract 1 from the loop size for each SILHS importance category.
       tot_loops = tot_loops - silhs_num_importance_categories
       ! Add 1 for "silhs_variance_category" to the loop size.
       tot_loops = tot_loops + 1
    end if

    if ( any( vars_lh_zt == "lh_samp_frac_category" ) ) then
       ! Correct for number of variables found under "lh_samp_frac_category".
       ! Subtract 1 from the loop size for each SILHS importance category.
       tot_loops = tot_loops - silhs_num_importance_categories
       ! Add 1 for "lh_samp_frac_category" to the loop size.
       tot_loops = tot_loops + 1
    end if

    k = 1
    do i = 1, tot_loops

      select case ( trim( vars_lh_zt(i) ) )
      case ( 'AKm' )           ! Vince Larson 22 May 2005
        stats_metadata%iAKm = k
        call stat_assign( var_index=stats_metadata%iAKm, & ! In
             var_name="AKm", & ! In
             var_description="Analytic Kessler ac [kg/kg]", var_units="kg/kg", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'lh_AKm' )       ! Vince Larson 22 May 2005
        stats_metadata%ilh_AKm = k

        call stat_assign( var_index=stats_metadata%ilh_AKm, & ! In
             var_name="lh_AKm", & ! In
             var_description="LH Kessler estimate  [kg/kg/s]", var_units="kg/kg/s", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'AKstd' )
        stats_metadata%iAKstd = k

        call stat_assign( var_index=stats_metadata%iAKstd, & ! In
             var_name="AKstd", & ! In
             var_description="Exact standard deviation of gba Kessler [kg/kg/s]", & ! In
             var_units="kg/kg/s", l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'AKstd_cld' )
        stats_metadata%iAKstd_cld = k

        call stat_assign( var_index=stats_metadata%iAKstd_cld, & ! In
             var_name="AKstd_cld", & ! In
             var_description="Exact w/in cloud std of gba Kessler [kg/kg/s]", & ! In
             var_units="kg/kg/s", l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'AKm_rcm' )
        stats_metadata%iAKm_rcm = k

        call stat_assign( var_index=stats_metadata%iAKm_rcm, & ! In
             var_name="AKm_rcm", & ! In
             var_description="Exact local gba auto based on rcm [kg/kg/s]", & ! In
             var_units="kg/kg/s", &  ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'AKm_rcc' )
        stats_metadata%iAKm_rcc = k

        call stat_assign( var_index=stats_metadata%iAKm_rcc, & ! In
             var_name="AKm_rcc", & ! In
             var_description="Exact local gba based on w/in cloud rc [kg/kg/s]", & ! In
             var_units="kg/kg/s", l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'lh_rvm_mc' )
        stats_metadata%ilh_rvm_mc = k

        call stat_assign( var_index=stats_metadata%ilh_rvm_mc, & ! In
             var_name="lh_rvm_mc", & ! In
             var_description="Latin hypercube estimate of rvm_mc [kg/kg/s]", & ! In
             var_units="kg/kg/s", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt )  ! In/Out
        k = k + 1

      case ( 'lh_thlm_mc' )
        stats_metadata%ilh_thlm_mc = k

        call stat_assign( var_index=stats_metadata%ilh_thlm_mc, & ! In
             var_name="lh_thlm_mc", & ! In
             var_description="Latin hypercube estimate of thlm_mc [kg/kg/s]", & ! In
             var_units="kg/kg/s", l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'lh_rcm_mc' )
        stats_metadata%ilh_rcm_mc = k

        call stat_assign( var_index=stats_metadata%ilh_rcm_mc, & ! In
             var_name="lh_rcm_mc", & ! In
             var_description="Latin hypercube estimate of rcm_mc [kg/kg/s]", & ! In
             var_units="kg/kg/s", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'lh_Ncm_mc' )
        stats_metadata%ilh_Ncm_mc = k

        call stat_assign( var_index=stats_metadata%ilh_Ncm_mc, & ! In
             var_name="lh_Ncm_mc", & ! In
             var_description="Latin hypercube estimate of Ncm_mc [kg/kg/s]", & ! In
             var_units="kg/kg/s", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'lh_rrm_mc' )
        stats_metadata%ilh_rrm_mc = k

        call stat_assign( var_index=stats_metadata%ilh_rrm_mc, & ! In
             var_name="lh_rrm_mc", & ! In
             var_description="Latin hypercube estimate of rrm_mc [kg/kg/s]", & ! In
             var_units="kg/kg/s", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'lh_Nrm_mc' )
        stats_metadata%ilh_Nrm_mc = k

        call stat_assign( var_index=stats_metadata%ilh_Nrm_mc, & ! In
             var_name="lh_Nrm_mc", & ! In
             var_description="Latin hypercube estimate of Nrm_mc [kg/kg/s]", & ! In
             var_units="kg/kg/s", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case('lh_rsm_mc')
        stats_metadata%ilh_rsm_mc = k

        call stat_assign( var_index=stats_metadata%ilh_rsm_mc, & ! In
             var_name="lh_rsm_mc", & ! In
             var_description="Latin hypercube estimate of rsm_mc [kg/kg/s]", & ! In
             var_units="kg/kg/s", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'lh_Nsm_mc' )
        stats_metadata%ilh_Nsm_mc = k

        call stat_assign( var_index=stats_metadata%ilh_Nsm_mc, & ! In
             var_name="lh_Nsm_mc", & ! In
             var_description="Latin hypercube estimate of Nsm_mc [kg/kg/s]", & ! In
             var_units="kg/kg/s", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'lh_rgm_mc' )
        stats_metadata%ilh_rgm_mc = k

        call stat_assign( var_index=stats_metadata%ilh_rgm_mc, & ! In
             var_name="lh_rgm_mc", & ! In
             var_description="Latin hypercube estimate of rgm_mc [kg/kg/s]", & ! In
             var_units="kg/kg/s", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt )  ! In/Out
        k = k + 1

      case ( 'lh_Ngm_mc' )
        stats_metadata%ilh_Ngm_mc = k

        call stat_assign( var_index=stats_metadata%ilh_Ngm_mc, & ! In
             var_name="lh_Ngm_mc", & ! In
             var_description="Latin hypercube estimate of Ngm_mc [kg/kg/s]", & ! In
             var_units="kg/kg/s", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'lh_rim_mc' )
        stats_metadata%ilh_rim_mc = k

        call stat_assign( var_index=stats_metadata%ilh_rim_mc, & ! In
             var_name="lh_rim_mc", & ! In
             var_description="Latin hypercube estimate of rim_mc [kg/kg/s]", & ! In
             var_units="kg/kg/s", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'lh_Nim_mc' )
        stats_metadata%ilh_Nim_mc = k

        call stat_assign( var_index=stats_metadata%ilh_Nim_mc, & ! In
             var_name="lh_Nim_mc", & ! In
             var_description="Latin hypercube estimate of Nim_mc [kg/kg/s]", & ! In
             var_units="kg/kg/s", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'lh_Vrr' )
        stats_metadata%ilh_Vrr = k

        call stat_assign( var_index=stats_metadata%ilh_Vrr, & ! In
             var_name="lh_Vrr", & ! In
             var_description="Latin hypercube estimate of rrm sedimentation velocity [m/s]", & !In
             var_units="m/s", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'lh_VNr' )
        stats_metadata%ilh_VNr = k

        call stat_assign( var_index=stats_metadata%ilh_VNr, & ! In
             var_name="lh_VNr", & ! In
             var_description="Latin hypercube estimate of Nrm sedimentation velocity [m/s]", & ! In
             var_units="m/s", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'lh_rcm_avg' )
        stats_metadata%ilh_rcm_avg = k

        call stat_assign( var_index=stats_metadata%ilh_rcm_avg, & ! In
             var_name="lh_rcm_avg", & ! In
             var_description="Latin hypercube average estimate of rcm [kg/kg]", & ! In
             var_units="kg/kg", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out

        k = k + 1

      case ( 'lh_rrm' )
        stats_metadata%ilh_rrm = k
 
        call stat_assign( var_index=stats_metadata%ilh_rrm, & ! In
             var_name="lh_rrm", & ! In
             var_description="Latin hypercube estimate of rrm [kg/kg]", & ! In
             var_units="kg/kg", & ! In
             l_silhs=.true., &  ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'lh_Nrm' )
        stats_metadata%ilh_Nrm = k

        call stat_assign( var_index=stats_metadata%ilh_Nrm, & ! In
             var_name="lh_Nrm", & ! In
             var_description="Latin hypercube estimate of Nrm [count/kg]", & ! In
             var_units="count/kg", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'lh_rim' )
        stats_metadata%ilh_rim = k

        call stat_assign( var_index=stats_metadata%ilh_rim, & ! In
             var_name="lh_rim", & ! In
             var_description="Latin hypercube estimate of rim [kg/kg]", & ! In
             var_units="kg/kg", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt )  ! In/Out
        k = k + 1

      case ( 'lh_Nim' )
        stats_metadata%ilh_Nim = k

        call stat_assign( var_index=stats_metadata%ilh_Nim, & ! In
             var_name="lh_Nim", & ! In
             var_description="Latin hypercube estimate of Nim [count/kg]", & ! In
             var_units="count/kg", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'lh_rsm' )
        stats_metadata%ilh_rsm = k

        call stat_assign( var_index=stats_metadata%ilh_rsm, & ! In
             var_name="lh_rsm", & ! In
             var_description="Latin hypercube estimate of rsm [kg/kg]", & ! In
             var_units="kg/kg", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'lh_Nsm' )
        stats_metadata%ilh_Nsm = k

        call stat_assign( var_index=stats_metadata%ilh_Nsm, & ! In
             var_name="lh_Nsm", & ! In
             var_description="Latin hypercube estimate of Nsm [count/kg]", & ! In
             var_units="count/kg", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1


      case ( 'lh_rgm' )
        stats_metadata%ilh_rgm = k

        call stat_assign( var_index=stats_metadata%ilh_rgm, & ! In
             var_name="lh_rgm", & ! In
             var_description="Latin hypercube estimate of rgm [kg/kg]", & ! In
             var_units="kg/kg", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'lh_Ngm' )
        stats_metadata%ilh_Ngm = k

        call stat_assign( var_index=stats_metadata%ilh_Ngm, & ! In
             var_name="lh_Ngm", & ! In
             var_description="Latin hypercube estimate of Ngm [kg/kg]", & ! In
             var_units="kg/kg", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'lh_thlm' )
        stats_metadata%ilh_thlm = k

        call stat_assign( var_index=stats_metadata%ilh_thlm, & ! In
             var_name="lh_thlm", & ! In
             var_description="Latin hypercube estimate of thlm [K]", var_units="K", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'lh_rcm' )
        stats_metadata%ilh_rcm = k

        call stat_assign( var_index=stats_metadata%ilh_rcm, & ! In
             var_name="lh_rcm", & ! In
             var_description="Latin hypercube estimate of rcm [kg/kg]", & ! In
             var_units="kg/kg", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'lh_Ncm' )
        stats_metadata%ilh_Ncm = k

        call stat_assign( var_index=stats_metadata%ilh_Ncm, & ! In
             var_name="lh_Ncm", & ! In
             var_description="Latin hypercube estimate of Ncm [count/kg]", & ! In
             var_units="count/kg", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'lh_Ncnm' )
        stats_metadata%ilh_Ncnm = k

        call stat_assign( var_index=stats_metadata%ilh_Ncnm, & ! In
             var_name="lh_Ncnm", & ! In
             var_description="Latin hypercube estimate of Ncnm [count/kg]", & ! In
             var_units="count/kg", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1


      case ( 'lh_rvm' )
        stats_metadata%ilh_rvm = k

        call stat_assign( var_index=stats_metadata%ilh_rvm, & ! In
             var_name="lh_rvm", & ! In
             var_description="Latin hypercube estimate of rvm [kg/kg]", & ! In
             var_units="kg/kg", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'lh_wm' )
        stats_metadata%ilh_wm = k

        call stat_assign( var_index=stats_metadata%ilh_wm, & ! In
             var_name="lh_wm", & ! In
             var_description="Latin hypercube estimate of vertical velocity [m/s]", & ! In
             var_units="m/s", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1 

      case ( 'lh_cloud_frac' )
        stats_metadata%ilh_cloud_frac = k

        ! Note: count is the udunits compatible unit
        call stat_assign( var_index=stats_metadata%ilh_cloud_frac, & ! In
             var_name="lh_cloud_frac", & ! In
             var_description="Latin hypercube estimate of cloud fraction [count]", & ! In
             var_units="count", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'lh_cloud_frac_unweighted' )
        stats_metadata%ilh_cloud_frac_unweighted = k

        call stat_assign( var_index=stats_metadata%ilh_cloud_frac_unweighted, & ! In
             var_name="lh_cloud_frac_unweighted", var_description="Unweighted fraction of & 
            &silhs sample points that are in cloud [-]", var_units="-", l_silhs=.false., & 
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'lh_chi' )
        stats_metadata%ilh_chi = k
        call stat_assign( var_index=stats_metadata%ilh_chi, & ! In
             var_name="lh_chi", & ! In
             var_description="Latin hypercube estimate of Mellor's s (extended liq) [kg/kg]", & !In
             var_units="kg/kg", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'lh_eta' )
        stats_metadata%ilh_eta = k
        call stat_assign( var_index=stats_metadata%ilh_eta, & ! In
             var_name="lh_eta", & ! In
             var_description="Latin hypercube estimate of Mellor's t [kg/kg]", & ! In
             var_units="kg/kg", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'lh_chip2' )
        stats_metadata%ilh_chip2 = k
        call stat_assign( var_index=stats_metadata%ilh_chip2, & ! In
             var_name="lh_chip2", & ! In
             var_description="Latin hypercube estimate of variance of chi(s) [kg/kg]", &! In
             var_units="kg/kg", & ! In
             l_silhs=.true., &  ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'lh_wp2_zt' )
        stats_metadata%ilh_wp2_zt = k
        call stat_assign( var_index=stats_metadata%ilh_wp2_zt, & ! In
             var_name="lh_wp2_zt", & ! In
             var_description="Variance of the latin hypercube estimate of w [m^2/s^2]",&! In
             var_units="m^2/s^2", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'lh_Ncnp2_zt' )
        stats_metadata%ilh_Ncnp2_zt = k
        call stat_assign( var_index=stats_metadata%ilh_Ncnp2_zt, & ! In
             var_name="lh_Ncnp2_zt", & ! In
             var_description="Variance of the latin hypercube estimate of Ncn [count^2/kg^2]", &!In
             var_units="count^2/kg^2", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'lh_Ncp2_zt' )
        stats_metadata%ilh_Ncp2_zt = k
        call stat_assign( var_index=stats_metadata%ilh_Ncp2_zt, & ! In
             var_name="lh_Ncp2_zt", & ! In
             var_description="Variance of the latin hypercube estimate of Nc [count^2/kg^2]", & !In
             var_units="count^2/kg^2", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'lh_Nrp2_zt' )
        stats_metadata%ilh_Nrp2_zt = k
        call stat_assign( var_index=stats_metadata%ilh_Nrp2_zt, & ! In
             var_name="lh_Nrp2_zt", & ! In
             var_description="Variance of the latin hypercube estimate of Nr [count^2/kg^2]", &!In
             var_units="count^2/kg^2", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'lh_rcp2_zt' )
        stats_metadata%ilh_rcp2_zt = k
        call stat_assign( var_index=stats_metadata%ilh_rcp2_zt, & ! In
             var_name="lh_rcp2_zt", & ! In
             var_description="Variance of the latin hypercube estimate of rc [kg^2/kg^2]", & ! In
             var_units="kg^2/kg^2", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'lh_rtp2_zt' )
        stats_metadata%ilh_rtp2_zt = k
        call stat_assign( var_index=stats_metadata%ilh_rtp2_zt, & ! In
             var_name="lh_rtp2_zt", & ! In
             var_description="Variance of the latin hypercube estimate of rt [kg^2/kg^2]", & ! In
             var_units="kg^2/kg^2", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'lh_thlp2_zt' )
        stats_metadata%ilh_thlp2_zt = k
        call stat_assign( var_index=stats_metadata%ilh_thlp2_zt, & ! In
             var_name="lh_thlp2_zt", & ! In
             var_description="Variance of the latin hypercube estimate of thl [K^2]", & ! In
             var_units="K^2", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'lh_rrp2_zt' )
        stats_metadata%ilh_rrp2_zt = k
        call stat_assign( var_index=stats_metadata%ilh_rrp2_zt, & ! In
             var_name="lh_rrp2_zt", & ! In
             var_description="Variance of the latin hypercube estimate of rr [kg^2/kg^2]", & ! In
             var_units="kg^2/kg^2", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'lh_rrm_auto' )
        stats_metadata%ilh_rrm_auto = k
        call stat_assign( var_index=stats_metadata%ilh_rrm_auto, & ! In
             var_name="lh_rrm_auto", & ! In
             var_description="Latin hypercube estimate of autoconversion [kg/kg/s]", & ! In
             var_units="kg/kg/s", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'lh_rrm_accr' )
        stats_metadata%ilh_rrm_accr = k
        call stat_assign( var_index=stats_metadata%ilh_rrm_accr, & ! In
             var_name="lh_rrm_accr", & ! In
             var_description="Latin hypercube estimate of accretion [kg/kg/s]", & ! In
             var_units="kg/kg/s", l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'lh_rrm_evap' )
        stats_metadata%ilh_rrm_evap = k
        call stat_assign( var_index=stats_metadata%ilh_rrm_evap, & ! In
             var_name="lh_rrm_evap", & ! In
             var_description="Latin hypercube estimate of evaporation [kg/kg/s]", & ! In
             var_units="kg/kg/s", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'lh_Nrm_auto' )
        stats_metadata%ilh_Nrm_auto = k
        call stat_assign( var_index=stats_metadata%ilh_Nrm_auto, & ! In
             var_name="lh_Nrm_auto", & ! In
             var_description="Latin hypercube estimate of Nrm autoconversion [num/kg/s]", & ! In
             var_units="num/kg/s", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'lh_Nrm_evap' )
        stats_metadata%ilh_Nrm_evap = k
        call stat_assign( var_index=stats_metadata%ilh_Nrm_evap, & ! In
             var_name="lh_Nrm_evap", & ! In
             var_description="Latin hypercube estimate of Nrm evaporation [num/kg/s]", & ! In
             var_units="num/kg/s", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'lh_rrm_src_adj' )
        stats_metadata%ilh_rrm_src_adj = k
        call stat_assign( var_index=stats_metadata%ilh_rrm_src_adj, & ! In
             var_name="lh_rrm_src_adj", & ! In
             var_description="Latin hypercube estimate of source adjustment (KK only!) [kg/kg/s]",&
             var_units="kg/kg/s", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'lh_rrm_evap_adj' )
        stats_metadata%ilh_rrm_evap_adj = k
        call stat_assign( var_index=stats_metadata%ilh_rrm_evap_adj, & ! In
             var_name="lh_rrm_evap_adj", & ! In
             var_description="Latin hypercube estimate of evap adjustment (KK only!) [kg/kg/s]", &
             var_units="kg/kg/s", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'lh_Nrm_src_adj' )
        stats_metadata%ilh_Nrm_src_adj = k
        call stat_assign( var_index=stats_metadata%ilh_Nrm_src_adj, & ! In
             var_name="lh_Nrm_src_adj", & ! In
             var_description="Latin hypercube estimate of Nrm source adjustment (KK only!) &
             &[kg/kg/s]", &
             var_units="kg/kg/s", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'lh_Nrm_evap_adj' )
        stats_metadata%ilh_Nrm_evap_adj = k
        call stat_assign( var_index=stats_metadata%ilh_Nrm_evap_adj, & ! In
             var_name="lh_Nrm_evap_adj", & ! In
             var_description="Latin hypercube estimate of Nrm evap adjustment (KK only!) &
             &[kg/kg/s]", &
             var_units="kg/kg/s", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'lh_precip_frac' )
        stats_metadata%ilh_precip_frac = k
        call stat_assign( var_index=stats_metadata%ilh_precip_frac, & ! In
             var_name="lh_precip_frac", & ! In
             var_description="Latin hypercube estimate of precipitation fraction [-]", & ! In
             var_units="-", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'lh_precip_frac_unweighted' )
        stats_metadata%ilh_precip_frac_unweighted = k
        call stat_assign( var_index=stats_metadata%ilh_precip_frac_unweighted, & ! In
             var_name="lh_precip_frac_unweighted", & ! In
             var_description="Unweighted fraction of sample points in precipitation [-]", & ! In
             var_units="-", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'lh_mixt_frac' )
        stats_metadata%ilh_mixt_frac = k
        call stat_assign( var_index=stats_metadata%ilh_mixt_frac, & ! In
             var_name="lh_mixt_frac", & ! In
             var_description="Latin hypercube estimate of mixture fraction (weight of 1st PDF &
             &component [-]", &
             var_units="-", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'lh_mixt_frac_unweighted' )
        stats_metadata%ilh_mixt_frac_unweighted = k
        call stat_assign( var_index=stats_metadata%ilh_mixt_frac_unweighted, & ! In
             var_name="lh_mixt_frac_unweighted", & ! In
             var_description="Unweighted fraction of sample points in first PDF component [-]",&!In
             var_units="-", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'lh_m_vol_rad_rain' )
        stats_metadata%ilh_m_vol_rad_rain = k
        call stat_assign( var_index=stats_metadata%ilh_m_vol_rad_rain, & ! In
             var_name="lh_m_vol_rad_rain", & ! In
             var_description="SILHS est. of rain radius", var_units="m", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'lh_rrm_mc_nonadj' )
        stats_metadata%ilh_rrm_mc_nonadj = k
        call stat_assign( var_index=stats_metadata%ilh_rrm_mc_nonadj, & ! In
             var_name="lh_rrm_mc_nonadj", & ! In
             var_description="SILHS est. of rrm_mc_nonadj [kg/kg/s]", var_units="kg/kg/s", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_zt ) ! In/Out
        k = k + 1

      case ( 'silhs_variance_category' )

        do icategory=1, silhs_num_importance_categories

          stats_metadata%isilhs_variance_category(icategory) = k
          write(category_num_as_string,'(I1)') icategory
          call stat_assign( var_index=stats_metadata%isilhs_variance_category(icategory), & ! In
               var_name="silhs_var_cat_"//category_num_as_string, & ! In
               var_description="Variance of SILHS variable in importance category " // &! In
               category_num_as_string, var_units="various", & ! In
               l_silhs=.false., & ! In
               grid_kind=stats_lh_zt ) ! In/Out
          k = k + 1

        end do

      case ( 'lh_samp_frac_category' )

        do icategory=1, silhs_num_importance_categories

          stats_metadata%ilh_samp_frac_category(icategory) = k
          write(category_num_as_string,'(I1)') icategory
          call stat_assign( var_index=stats_metadata%ilh_samp_frac_category(icategory), & ! In
               var_name="lh_samp_frac_"//category_num_as_string, & ! In
               var_description="Number of samples in importance category " // & ! In
               category_num_as_string // " [-]", var_units="-", l_silhs=.false., & ! In
               grid_kind=stats_lh_zt ) ! In/Out
          k = k + 1

        end do

      case default

        write(fstderr,*) 'Error:  unrecognized variable in vars_lh_zt:  ', trim( vars_lh_zt(i) )

        l_error = .true.  ! This will stop the run.

      end select

    end do ! i = 1, stats_lh_zt%num_output_fields

    return
  end subroutine stats_init_lh_zt

end module stats_lh_zt_module
