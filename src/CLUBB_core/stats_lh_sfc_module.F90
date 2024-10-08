!-----------------------------------------------------------------------
! $Id$
!===============================================================================

module stats_lh_sfc_module

  implicit none

  private ! Set Default Scope

  public :: stats_init_lh_sfc

  ! Constant parameters
  integer, parameter, public :: nvarmax_lh_sfc = 10  ! Maximum variables allowed

  contains

!-----------------------------------------------------------------------
  subroutine stats_init_lh_sfc( vars_lh_sfc,                    & ! In
                                l_error,                        & ! In/Out
                                stats_metadata, stats_lh_sfc )    ! In/Out

! Description:
!   Initializes array indices for stats_lh_sfc
! References:
!   None
!-----------------------------------------------------------------------

    use constants_clubb, only: &
        fstderr ! Constant(s)
      
    use stats_type_utilities, only: & 
        stat_assign ! Procedure

    use stats_type, only: &
        stats ! Type

    use stats_variables, only: &
        stats_metadata_type

    implicit none

    ! External
    intrinsic :: trim

    !--------------------- Input Variable ---------------------
    character(len= * ), dimension(nvarmax_lh_sfc), intent(in) :: vars_lh_sfc

    !--------------------- InOut Variables ---------------------      
    type (stats_metadata_type), intent(inout) :: &
      stats_metadata

    type (stats), intent(inout) :: &
      stats_lh_sfc

    logical, intent(inout) :: l_error

    !--------------------- Local Varables ---------------------
    integer :: i, k

    !--------------------- Begin Code ---------------------

    ! Default initialization for array indices for stats_sfc is zero (see module
    ! stats_variables)

    ! Assign pointers for statistics variables stats_sfc

    k = 1
    do i = 1, stats_lh_sfc%num_output_fields

      select case ( trim( vars_lh_sfc(i) ) )

      case ( 'lh_morr_snow_rate' )
        stats_metadata%ilh_morr_snow_rate = k
        call stat_assign( var_index=stats_metadata%ilh_morr_snow_rate, & ! In
             var_name="lh_morr_snow_rate", & ! In
             var_description="Snow+Ice+Graupel fallout rate from Morrison scheme [mm/day]", & ! In
             var_units="mm/day", l_silhs=.true., & ! In
             grid_kind=stats_lh_sfc ) ! In/Out
        k = k + 1

      case ( 'lh_vwp' )
        stats_metadata%ilh_vwp = k
        call stat_assign( var_index=stats_metadata%ilh_vwp, & ! In
             var_name="lh_vwp", & ! In
             var_description="Vapor water path [kg/m^2]", var_units="kg/m^2", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_sfc ) ! In/Out
        k = k + 1

      case ( 'lh_lwp' )
        stats_metadata%ilh_lwp = k
        call stat_assign( var_index=stats_metadata%ilh_lwp, & ! In
             var_name="lh_lwp", & ! In
             var_description="Liquid water path [kg/m^2]", var_units="kg/m^2", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_sfc ) ! In/Out
        k = k + 1

      case ( 'k_lh_start' )
        stats_metadata%ik_lh_start = k
        call stat_assign( var_index=stats_metadata%ik_lh_start, & ! In
             var_name="k_lh_start", & ! In
             var_description="Index of height level for SILHS sampling preferentially within &
                             &cloud [integer]", var_units="integer", l_silhs=.true., & 
             grid_kind=stats_lh_sfc ) ! In/Out
        k = k + 1

      case ( 'lh_sample_weights_sum' )
        stats_metadata%ilh_sample_weights_sum = k
        call stat_assign( var_index=stats_metadata%ilh_sample_weights_sum, & ! In
             var_name="lh_sample_weights_sum", & ! In
             var_description="Sum of the sample point weights [-]", var_units="-", & ! In
             l_silhs=.true., & ! In
             grid_kind=stats_lh_sfc ) ! In/Out
        k = k + 1
        
      case ( 'lh_sample_weights_avg' )
        stats_metadata%ilh_sample_weights_avg = k
        call stat_assign( var_index=stats_metadata%ilh_sample_weights_avg, & ! In
             var_name="lh_sample_weights_avg", & ! In
             var_description="Average of the sample point weights [-]", &  ! In
             var_units="-", & ! In
             l_silhs=.true., &  ! In
             grid_kind=stats_lh_sfc ) ! In/Out
        k = k + 1

      case default
        write(fstderr,*) 'Error:  unrecognized variable in vars_lh_sfc:  ',  &
              trim( vars_lh_sfc(i) )
        l_error = .true.  ! This will stop the run.

      end select

    end do ! i = 1, stats_lh_sfc%num_output_fields

    return
  end subroutine stats_init_lh_sfc

end module stats_lh_sfc_module

