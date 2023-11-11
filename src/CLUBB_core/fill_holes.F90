!-----------------------------------------------------------------------
! $Id: fill_holes.F90 8738 2018-07-19 19:58:53Z bmg2@uwm.edu $
!===============================================================================
module fill_holes

  implicit none

  public :: fill_holes_driver, &
            fill_holes_vertical, &
            hole_filling_hm_one_lev, &
            fill_holes_hydromet, &
            fill_holes_wv, &
            clip_hydromet_conc_mvr, &
            setup_stats_indices

  private ! Set Default Scope

  contains

  !=============================================================================
  subroutine fill_holes_vertical( nz, ngrdcol, num_draw_pts, threshold, upper_hf_level, &
                                  dz, rho_ds, &
                                  field )

    ! Description:
    !   This subroutine clips values of 'field' that are below 'threshold' as much
    !   as possible (i.e. "fills holes"), but conserves the total integrated mass
    !   of 'field'.  This prevents clipping from acting as a spurious source.
    !
    !   Mass is conserved by reducing the clipped field everywhere by a constant
    !   multiplicative coefficient.
    !
    !   This subroutine does not guarantee that the clipped field will exceed
    !   threshold everywhere; blunt clipping is needed for that.
    !
    !   The lowest level (k=1) should not be included, as the hole-filling scheme
    !   should not alter the set value of 'field' at the surface (for momentum
    !   level variables), or consider the value of 'field' at a level below the
    !   surface (for thermodynamic level variables).  
    !
    !   For momentum level variables only, the hole-filling scheme should not 
    !   alter the set value of 'field' at the upper boundary level (k=nz).
    !   So for momemtum level variables, call with upper_hf_level=nz-1, and
    !   for thermodynamic level variables, call with upper_hf_level=nz.
    !
    ! References:
    !   ``Numerical Methods for Wave Equations in Geophysical Fluid
    !     Dynamics'', Durran (1999), p. 292.
    !-----------------------------------------------------------------------

    use grid_class, only: & 
        grid ! Type

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use constants_clubb, only: &
        eps, &
        one

    implicit none
    
    ! --------------------- Input variables ---------------------
    integer, intent(in) :: &
      nz, &
      ngrdcol
    
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      dz      ! Spacing between thermodynamic grid levels; centered over
              ! momentum grid levels
              ! OR
              ! Spcaing between momentum grid levels; centered over
              ! thermodynamic grid levels
                  
    integer, intent(in) :: & 
      num_draw_pts, & ! The number of points on either side of the hole;
                      ! Mass is drawn from these points to fill the hole.  []
      upper_hf_level  ! Upper grid level of global hole-filling range      []

    real( kind = core_rknd ), intent(in) :: & 
      threshold  ! A threshold (e.g. w_tol*w_tol) below which field must not
                 ! fall                           [Units vary; same as field]

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) ::  & 
      rho_ds       ! Dry, static density on thermodynamic or momentum levels    [kg/m^3]

    ! --------------------- Input/Output variable ---------------------
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(inout) :: & 
      field  ! The field (e.g. wp2) that contains holes [Units same as threshold]
 
    ! --------------------- Local Variables ---------------------
    integer :: & 
      i,             & ! Loop index for column                           []
      k,             & ! Loop index for absolute grid level              []
      j, &
      k_start,     & ! Lower grid level of local hole-filling range    []
      k_end          ! Upper grid level of local hole-filling range    []

    real( kind = core_rknd ), dimension(ngrdcol,nz)  ::  & 
      rho_ds_dz,            & ! rho_ds * dz
      invrs_denom_integral, & ! Inverse of the integral in the denominator (see description)
      field_clipped           ! The raw field (e.g. wp2) that contains no holes
                              !                          [Units same as threshold]

    real( kind = core_rknd ) ::  & 
      field_avg,          & ! Vertical average of field [Units of field]
      field_clipped_avg,  & ! Vertical average of clipped field [Units of field]
      mass_fraction         ! Coefficient that multiplies clipped field
                            ! in order to conserve mass.                      []

    real( kind = core_rknd ), dimension(ngrdcol) ::  & 
      denom_integral_global,  & ! Integral in the denominator for global filling
      numer_integral_global,  & ! Integral in the numerator for global filling
      field_avg_global,       & ! Vertical average of field [Units of field]
      mass_fraction_global      ! Coefficient that multiplies clipped field
                                ! in order to conserve mass.                      []

    logical :: &
      l_field_below_threshold

    ! --------------------- Begin Code --------------------- 

    !$acc enter data create( invrs_denom_integral, field_clipped, denom_integral_global, rho_ds_dz, &
    !$acc                    numer_integral_global, field_avg_global, mass_fraction_global )

    l_field_below_threshold = .false.

    !$acc parallel loop gang vector collapse(2) default(present) &
    !$acc          reduction(.or.:l_field_below_threshold)
    do k = 1, nz
      do i = 1, ngrdcol
        if ( field(i,k) < threshold ) then
          l_field_below_threshold = .true.
        end if
      end do
    end do
    !$acc end parallel loop

    ! If all field values are above the specified threshold, no hole filling is required
    if ( .not. l_field_below_threshold ) then
      !$acc exit data delete( invrs_denom_integral, field_clipped, denom_integral_global, rho_ds_dz, &
      !$acc                   numer_integral_global, field_avg_global, mass_fraction_global )
      return
    end if

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, nz
      do i = 1, ngrdcol
        rho_ds_dz(i,k) = rho_ds(i,k) * dz(i,k)
      end do  
    end do
    !$acc end parallel loop
    

    ! denom_integral does not change throughout the hole filling algorithm
    ! so we can calculate it before hand. This results in unneccesary computations,
    ! but is parallelizable and reduces the cost of the serial k loop
    !$acc parallel loop gang vector collapse(2) default(present)
    do i = 1, ngrdcol
      do k = 2+num_draw_pts, upper_hf_level-num_draw_pts
        k_start = k - num_draw_pts
        k_end   = k + num_draw_pts
        invrs_denom_integral(i,k) = one / sum( rho_ds_dz(i,k_start:k_end) )
      end do  
    end do
    !$acc end parallel loop
    
    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol

      ! Make one pass up the profile, filling holes as much as we can using
      ! nearby mass, ignoring the first level.
      !
      ! This loop can only be done in serial due to the field values required for the next
      ! iteration potentially changing. We could in theory expose more parallelism in cases 
      ! where there are large enough gaps between vertical levels which need hole-filling,
      ! but levels which require hole-filling are often close or consecutive.
      do k = 2+num_draw_pts, upper_hf_level-num_draw_pts

        k_start = k - num_draw_pts
        k_end   = k + num_draw_pts

        if ( any( field(i,k_start:k_end) < threshold ) ) then

          ! Compute the field's vertical average cenetered at k, which we must conserve,
          ! see description of the vertical_avg function in advance_helper_module
          field_avg = sum( rho_ds_dz(i,k_start:k_end) * field(i,k_start:k_end) ) &
                      * invrs_denom_integral(i,k)

          ! Clip small or negative values from field.
          if ( field_avg >= threshold ) then
            ! We know we can fill in holes completely
            field_clipped(i,k_start:k_end) = max( threshold, field(i,k_start:k_end) )
          else
            ! We can only fill in holes partly;
            ! to do so, we remove all mass above threshold.
            field_clipped(i,k_start:k_end) = min( threshold, field(i,k_start:k_end) )
          endif

          ! Compute the clipped field's vertical integral.
          ! clipped_total_mass >= original_total_mass,
          ! see description of the vertical_avg function in advance_helper_module
          field_clipped_avg = sum( rho_ds_dz(i,k_start:k_end) * field_clipped(i,k_start:k_end) ) &
                              * invrs_denom_integral(i,k)

          ! Avoid divide by zero issues by doing nothing if field_clipped_avg ~= threshold
          if ( abs(field_clipped_avg-threshold) > abs(field_clipped_avg+threshold)*eps/2) then
            ! Compute coefficient that makes the clipped field have the same mass as the
            ! original field.  We should always have mass_fraction > 0.
            mass_fraction = ( field_avg - threshold ) &
                            / ( field_clipped_avg - threshold )

            ! Calculate normalized, filled field
            field(i,k_start:k_end) = threshold &
                        + mass_fraction * ( field_clipped(i,k_start:k_end) - threshold )
          endif

        endif

      end do

    end do
    !$acc end parallel loop

    l_field_below_threshold = .false.

    !$acc parallel loop gang vector collapse(2) default(present) reduction(.or.:l_field_below_threshold)
    do k = 1, nz
      do i = 1, ngrdcol
        if ( field(i,k) < threshold ) then
          l_field_below_threshold = .true.
        end if
      end do
    end do
    !$acc end parallel loop

    ! If all field values are above the threshold, no further hole filling is required
    if ( .not. l_field_below_threshold ) then
      !$acc exit data delete( invrs_denom_integral, field_clipped, denom_integral_global, rho_ds_dz, &
      !$acc                   numer_integral_global, field_avg_global, mass_fraction_global )
      return
    end if


    ! Now we fill holes globally to maximize the chance that all holes are filled.
    ! To improve parallelism we assume that global hole filling needs to be done 
    ! for each grid column, perform all calculations required, and only check
    ! if any holes need filling before the final step of updating the field. 

    ! Compute the numerator and denominator integrals
    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      numer_integral_global(i) = 0.0_core_rknd
      denom_integral_global(i) = 0.0_core_rknd
    end do
    !$acc end parallel loop

    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      do k = 2, upper_hf_level
        numer_integral_global(i) = numer_integral_global(i) + rho_ds_dz(i,k) * field(i,k)

        denom_integral_global(i) = denom_integral_global(i) + rho_ds_dz(i,k)
      end do  
    end do
    !$acc end parallel loop

    
    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol

      ! Find the vertical average of field, using the precomputed numerator and denominator,
      ! see description of the vertical_avg function in advance_helper_module
      field_avg_global(i) = numer_integral_global(i) / denom_integral_global(i)

      ! Clip small or negative values from field.
      if ( field_avg_global(i) >= threshold ) then
        ! We know we can fill in holes completely
        field_clipped(i,2:upper_hf_level) = max( threshold, field(i,2:upper_hf_level) )
      else
        ! We can only fill in holes partly;
        ! to do so, we remove all mass above threshold.
        field_clipped(i,2:upper_hf_level) = min( threshold, field(i,2:upper_hf_level) )
      endif

    end do
    !$acc end parallel loop

    ! To compute the clipped field's vertical integral we only need to recompute the numerator
    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol
      numer_integral_global(i) = 0.0_core_rknd
      do k = 2, upper_hf_level
        numer_integral_global(i) = numer_integral_global(i) + rho_ds_dz(i,k) * field_clipped(i,k)
      end do  
    end do
    !$acc end parallel loop

    !$acc parallel loop gang vector default(present)
    do i = 1, ngrdcol

      ! Do not complete calculations or update field values for this 
      ! column if there are no holes that need filling
      if ( any( field(i,2:upper_hf_level) < threshold ) ) then

        ! Compute the clipped field's vertical integral,
        ! see description of the vertical_avg function in advance_helper_module
        field_clipped_avg = numer_integral_global(i) / denom_integral_global(i)

        ! Avoid divide by zero issues by doing nothing if field_clipped_avg ~= threshold
        if ( abs(field_clipped_avg-threshold) > abs(field_clipped_avg+threshold)*eps/2) then

          ! Compute coefficient that makes the clipped field have the same mass as the
          ! original field.  We should always have mass_fraction > 0.
          mass_fraction_global(i) = ( field_avg_global(i) - threshold ) &
                                    / ( field_clipped_avg - threshold )

          ! Calculate normalized, filled field
          field(i,2:upper_hf_level) = threshold + mass_fraction_global(i) &
                                                  * ( field_clipped(i,2:upper_hf_level) &
                                                      - threshold )
        end if
      end if

    end do
    !$acc end parallel loop
    
    !$acc exit data delete( invrs_denom_integral, field_clipped, denom_integral_global, rho_ds_dz, &
    !$acc                   numer_integral_global, field_avg_global, mass_fraction_global )

    return

  end subroutine fill_holes_vertical

  !===============================================================================
  subroutine hole_filling_hm_one_lev( num_hm_fill, hm_one_lev, & ! Intent(in)
                                   hm_one_lev_filled ) ! Intent(out)

  ! Description:
  ! Fills holes between same-phase (i.e. either liquid or frozen) hydrometeors for
  ! one height level.
  !
  ! Warning: Do not input hydrometeors of different phases, e.g. liquid and frozen.
  ! Otherwise heat will not be conserved.
  !
  ! References:
  !
  ! None
  !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one, & ! Variable(s)
        zero, &
        eps

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use error_code, only: &
        clubb_at_least_debug_level  ! Procedure

    implicit none

    ! Input Variables
    integer, intent(in) :: num_hm_fill ! number of hydrometeors involved

    real(kind = core_rknd), dimension(num_hm_fill), intent(in) :: hm_one_lev

    ! Output Variables
    real(kind = core_rknd), dimension(num_hm_fill), intent(out) :: hm_one_lev_filled

    ! Local Variables
    integer :: num_neg_hm ! number of holes

    real(kind = core_rknd) :: &
      total_hole, & ! Size of the hole ( missing mass, less than 0 )
      total_mass    ! Total mass to fill the hole
      ! total mass of water substance = total_mass + total_hole

    integer :: i ! loop iterator

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    ! Initialization
    hm_one_lev_filled = 0._core_rknd
    total_hole = 0._core_rknd
    total_mass = 0._core_rknd
    num_neg_hm = 0

    ! Determine the total size of the hole and the number of neg. hydrometeors
    ! and the total mass of hole filling material
    do i=1, num_hm_fill
!       print *, "hm_one_lev(",i,") = ", hm_one_lev(i)
       if ( hm_one_lev(i) < zero ) then
          total_hole = total_hole + hm_one_lev(i) ! less than zero
          num_neg_hm = num_neg_hm + 1
       else
          total_mass = total_mass + hm_one_lev(i)
       endif

    enddo

!    print *, "total_hole = ", total_hole
!    print *, "total_mass = ", total_mass
!    print *, "num_neg_hm = ", num_neg_hm

    ! There is no water substance at all to fill the hole
    if ( abs(total_mass) < eps ) then

       if ( clubb_at_least_debug_level( 2 ) ) then
          print *, "Warning: One-level hole filling was not successful! total_mass ~= 0"
       endif

       hm_one_lev_filled = hm_one_lev

       return
    endif

    ! Fill the holes and adjust the remaining quantities:
    ! hm_filled(i) = 0, if hm(i) < 0
    ! or
    ! hm_filled(i) = (1 + total_hole/total_mass)*hm(i), if hm(i) > 0
    do i=1, num_hm_fill

       ! if there is not enough material, fill the holes partially with all the material available
       if ( abs(total_hole) > total_mass ) then

          if ( clubb_at_least_debug_level( 2 ) ) then
             print *, "Warning: One-level hole filling was not able to fill holes completely!" // &
                      " The holes were filled partially. |total_hole| > total_mass"
          endif

          hm_one_lev_filled(i) = min(hm_one_lev(i), zero) * ( one + total_mass / total_hole )

       else ! fill holes completely
          hm_one_lev_filled(i) = max(hm_one_lev(i), zero) * ( one + total_hole / total_mass )

       endif

    enddo

    ! Assertion checks (water substance conservation, non-negativity)
    if ( clubb_at_least_debug_level( 2 ) ) then

       if ( abs(sum( hm_one_lev ) - sum(hm_one_lev_filled)) > &
            abs(sum( hm_one_lev ) + sum(hm_one_lev_filled)) * eps/2 ) then
          print *, "Warning: Hole filling was not conservative!"
       endif

       if ( any( hm_one_lev_filled < zero ) ) then
          print *, "Warning: Hole filling failed! A hole could not be filled."
       endif

    endif

    return

  end subroutine hole_filling_hm_one_lev
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  subroutine fill_holes_hydromet( nz, hydromet_dim, hydromet, & ! Intent(in)
                                  hydromet_filled ) ! Intent(out)

  ! Description:
  ! Fills holes between same-phase hydrometeors(i.e. for frozen hydrometeors).
  ! The hole filling conserves water substance between all same-phase (frozen or liquid)
  ! hydrometeors at each height level.
  !
  ! Attention: The hole filling for the liquid phase hydrometeors is not yet implemented
  !
  ! Attention: l_frozen_hm and l_mix_rat_hm need to be set up before this subroutine is called!
  !
  ! References:
  !
  ! None
  !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd

    use array_index, only: &
        l_frozen_hm, & ! Variable(s)
        l_mix_rat_hm

    use constants_clubb, only: &
        zero

    implicit none

    ! Input Variables
    integer, intent(in) :: hydromet_dim, nz

    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(in) :: &
      hydromet

    ! Output Variables
    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(out) :: &
      hydromet_filled

    ! Local Variables
    integer :: i,j ! Loop iterators

    integer :: num_frozen_hm ! Number of frozen hydrometeor mixing ratios

    real( kind = core_rknd ), dimension(:,:), allocatable :: &
      hydromet_frozen,       & ! Frozen hydrometeor mixing ratios
      hydromet_frozen_filled   ! Frozen hydrometeor mixing ratios after hole filling

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    ! Determine the number of frozen hydrometeor mixing ratios
    num_frozen_hm = 0
    do i=1,hydromet_dim
       if ( l_frozen_hm(i) .and. l_mix_rat_hm(i) ) then
          num_frozen_hm = num_frozen_hm + 1
       endif
    enddo

    ! Allocation
    allocate( hydromet_frozen(nz,num_frozen_hm) )
    allocate( hydromet_frozen_filled(nz,num_frozen_hm) )

    ! Determine frozen hydrometeor mixing ratios
    j = 1
    do i = 1,hydromet_dim
       if ( l_frozen_hm(i) .and. l_mix_rat_hm(i) ) then
          hydromet_frozen(:,j) = hydromet(:,i)
          j = j+1
       endif
    enddo

    ! Fill holes for the frozen hydrometeors
    do i=1,nz
       if ( any( hydromet_frozen(i,:) < zero ) ) then
          call hole_filling_hm_one_lev( num_frozen_hm, hydromet_frozen(i,:), & ! Intent(in)
                                     hydromet_frozen_filled(i,:) ) ! Intent(out)
       else
          hydromet_frozen_filled(i,:) = hydromet_frozen(i,:)
       endif
    enddo

    ! Setup the filled hydromet array
    j = 1
    do i=1, hydromet_dim
       if ( l_frozen_hm(i) .and. l_mix_rat_hm(i) ) then
          hydromet_filled(:,i) = hydromet_frozen_filled(:,j)
          j = j+1
       else
          hydromet_filled(:,i) = hydromet(:,i)
       endif
    enddo

    !!! Here we could do the same hole filling for all the liquid phase hydrometeors

    return
  end subroutine fill_holes_hydromet
  !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
  subroutine fill_holes_wv( nz, dt, exner, hydromet_name, & ! Intent(in)
                            rvm_mc, thlm_mc, hydromet )! Intent(inout)

  ! Description:
  ! Fills holes using the cloud water mixing ratio from the current height level.
  !
  ! References:
  !
  ! None
  !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd

    use constants_clubb, only: &
        zero_threshold, &
        Lv, &
        Ls, &
        Cp

    implicit none

    ! Input Variables
    integer, intent(in) :: nz

    real( kind = core_rknd ), intent(in) ::  &
      dt           ! Timestep         [s]

    character(len=10), intent(in) :: hydromet_name

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      exner        ! Exner function                            [-]

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(nz), intent(inout) :: &
      hydromet, &  ! Hydrometeor array                         [units vary]
      rvm_mc, &
      thlm_mc

    ! Local Variables
    integer :: k ! Loop iterator

    real( kind = core_rknd ) :: rvm_clip_tndcy
  !-----------------------------------------------------------------------

    !----- Begin Code -----

    do k = 2, nz, 1

       if ( hydromet(k) < zero_threshold ) then

          ! Set rvm_clip_tndcy to the time tendency applied to vapor and removed
          ! from the hydrometeor.
          rvm_clip_tndcy = hydromet(k) / dt

          ! Adjust the tendency rvm_mc accordingly
          rvm_mc(k) = rvm_mc(k) + rvm_clip_tndcy

          ! Adjust the tendency of thlm_mc according to whether the
          ! effect is an evaporation or sublimation tendency.
          select case ( trim( hydromet_name ) )
          case( "rrm" )
             thlm_mc(k) = thlm_mc(k) - rvm_clip_tndcy * ( Lv / ( Cp*exner(k) ) )
          case( "rim", "rsm", "rgm" )
             thlm_mc(k) = thlm_mc(k) - rvm_clip_tndcy * ( Ls / ( Cp*exner(k) ) )
          case default
             error stop "Fatal error in microphys_driver"
          end select

          ! Set the mixing ratio to 0
          hydromet(k) = zero_threshold

       endif ! hydromet(k,i) < 0

    enddo ! k = 2..gr%nz

    return
  end subroutine fill_holes_wv
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  subroutine fill_holes_driver( gr, nz, dt, hydromet_dim,        & ! Intent(in)
                                l_fill_holes_hm,             & ! Intent(in)
                                rho_ds_zm, rho_ds_zt, exner, & ! Intent(in)
                                stats_metadata,              & ! Intent(in)
                                stats_zt,                    & ! intent(inout)
                                thlm_mc, rvm_mc, hydromet )    ! Intent(inout)

  ! Description:
  ! Fills holes between same-phase hydrometeors(i.e. for frozen hydrometeors).
  ! The hole filling conserves water substance between all same-phase (frozen or liquid)
  ! hydrometeors at each height level.
  !
  ! Attention: The hole filling for the liquid phase hydrometeors is not yet implemented
  !
  ! Attention: l_frozen_hm and l_mix_rat_hm need to be set up before this subroutine is called!
  !
  ! References:
  !
  ! None
  !-----------------------------------------------------------------------

    use grid_class, only: &
        grid ! Type

    use clubb_precision, only: &
        core_rknd   ! Variable(s)

    use constants_clubb, only: &
        zero,            &
        zero_threshold,  &
        Lv,              &
        Ls,              &
        Cp,              &
        fstderr,         &
        num_hf_draw_points

    use array_index, only: &
        hydromet_list, & ! Names of the hydrometeor species
        hydromet_tol

    use array_index, only: &
        l_mix_rat_hm, & ! Variable(s)
        l_frozen_hm

    use index_mapping, only: &
        Nx2rx_hm_idx, & ! Procedure(s)
        mvr_hm_max

    use stats_type_utilities, only: &
        stat_begin_update, & ! Subroutines
        stat_end_update

    use stats_variables, only: &
        stats_metadata_type

    use error_code, only: &
        clubb_at_least_debug_level  ! Procedure

    use stats_type, only: stats ! Type

    implicit none

    !----------------------- Input Variables -----------------------
    type (grid), target, intent(in) :: gr

    integer, intent(in) :: hydromet_dim, nz

    logical, intent(in) :: l_fill_holes_hm

    real( kind = core_rknd ), intent(in) ::  &
      dt           ! Timestep         [s]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      rho_ds_zm, & ! Dry, static density on momentum levels   [kg/m^3]
      rho_ds_zt    ! Dry, static density on thermo. levels    [kg/m^3]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      exner  ! Exner function                                       [-]

    type (stats_metadata_type), intent(in) :: &
      stats_metadata

    !----------------------- Input/Output Variables -----------------------
    type (stats), target, intent(inout) :: &
      stats_zt

    real( kind = core_rknd ), dimension(nz, hydromet_dim), intent(inout) :: &
      hydromet    ! Mean of hydrometeor fields    [units vary]

    real( kind = core_rknd ), dimension(nz), intent(inout) :: &
      rvm_mc,  & ! Microphysics contributions to vapor water           [kg/kg/s]
      thlm_mc    ! Microphysics contributions to liquid potential temp [K/s]

    !----------------------- Local Variables -----------------------
    integer :: i, k ! Loop iterators

    real( kind = core_rknd ), dimension(nz, hydromet_dim) :: &
      hydromet_filled,  & ! Frozen hydrometeor mixing ratios after hole filling
      hydromet_clipped    ! Clipped mean of hydrometeor fields    [units vary]

    character( len = 10 ) :: hydromet_name

    real( kind = core_rknd ) :: &
      max_velocity    ! Maximum sedimentation velocity                     [m/s]

    integer :: ixrm_hf, ixrm_wvhf, ixrm_cl, &
               ixrm_bt, ixrm_mc

    logical :: l_hole_fill = .true.

    !----------------------- Begin Code -----------------------

    ! Start stats output for the _hf variables (changes in the hydromet array
    ! due to fill_holes_hydromet and fill_holes_vertical)
    if ( stats_metadata%l_stats_samp ) then

       do i = 1, hydromet_dim

          ! Set up the stats indices for hydrometeor at index i
          call setup_stats_indices( i, stats_metadata,           & ! Intent(in)
                                    ixrm_bt, ixrm_hf, ixrm_wvhf, & ! Intent(inout)
                                    ixrm_cl, ixrm_mc,            & ! Intent(inout)
                                    max_velocity )                 ! Intent(inout)

          call stat_begin_update( gr%nz, ixrm_hf, hydromet(:,i) / dt, & ! intent(in)
                                  stats_zt ) ! intent(inout)

       enddo ! i = 1, hydromet_dim

    endif ! stats_metadata%l_stats_samp

    ! If we're dealing with negative hydrometeors, we first try to fill the
    ! holes proportionally from other same-phase hydrometeors at each height
    ! level.
    if ( any( hydromet < zero_threshold ) .and. l_fill_holes_hm ) then

       call fill_holes_hydromet( nz, hydromet_dim, hydromet, & ! Intent(in)
                                 hydromet_filled ) ! Intent(out)

       hydromet = hydromet_filled

    endif ! any( hydromet < zero ) .and. l_fill_holes_hm

    hydromet_filled = zero

    do i = 1, hydromet_dim

      ! Set up the stats indices for hydrometeor at index i
      call setup_stats_indices( i, stats_metadata,           & ! Intent(in)
                                ixrm_bt, ixrm_hf, ixrm_wvhf, & ! Intent(inout)
                                ixrm_cl, ixrm_mc,            & ! Intent(inout)
                                max_velocity )                 ! Intent(inout)

      hydromet_name = hydromet_list(i)

      ! Print warning message if any hydrometeor species has a value < 0.
      if ( clubb_at_least_debug_level( 1 ) ) then
         if ( any( hydromet(:,i) < zero_threshold ) ) then

            do k = 1, nz
               if ( hydromet(k,i) < zero_threshold ) then
                  write(fstderr,*) trim( hydromet_name ) //" < ", &
                                   zero_threshold, &
                                   " in fill_holes_driver at k= ", k
               endif ! hydromet(k,i) < 0
            enddo ! k = 1, nz

         endif ! hydromet(:,i) < 0       
      endif ! clubb_at_least_debug_level( 1 )


      ! Store the previous value of the hydrometeor for the effect of the
      ! hole-filling scheme.
!      if ( stats_metadata%l_stats_samp ) then
!         call stat_begin_update( ixrm_hf, hydromet(:,i) &
!                                          / dt, stats_zt )
!      endif

      ! If we're dealing with a mixing ratio and hole filling is enabled,
      ! then we apply the hole filling algorithm
      if ( any( hydromet(:,i) < zero_threshold ) ) then

         if ( hydromet_name(1:1) == "r" .and. l_hole_fill ) then

            !$acc data copyin( gr, gr%dzt, rho_ds_zt ) &
            !$acc        copy( hydromet(:,i) )

            ! Apply the hole filling algorithm
            ! upper_hf_level = nz since we are filling the zt levels
            call fill_holes_vertical( gr%nz, 1, num_hf_draw_points, zero_threshold, gr%nz, & ! In
                                      gr%dzt, rho_ds_zt,                                   & ! In
                                      hydromet(:,i) )                                      ! InOut

            !$acc end data

         endif ! Variable is a mixing ratio and l_hole_fill is true

      endif ! hydromet(:,i) < 0

      ! Enter the new value of the hydrometeor for the effect of the
      ! hole-filling scheme.
      if ( stats_metadata%l_stats_samp ) then
         call stat_end_update( gr%nz, ixrm_hf, hydromet(:,i) / dt, & ! intent(in)
                               stats_zt ) ! intent(inout)
      endif

      ! Store the previous value of the hydrometeor for the effect of the water
      ! vapor hole-filling scheme.
      if ( stats_metadata%l_stats_samp ) then
         call stat_begin_update( gr%nz, ixrm_wvhf, hydromet(:,i) / dt, & ! intent(in)
                                 stats_zt ) ! intent(inout)
      endif

      if ( any( hydromet(:,i) < zero_threshold ) ) then

         if ( hydromet_name(1:1) == "r" .and. l_hole_fill ) then

            ! If the hole filling algorithm failed, then we attempt to fill
            ! the missing mass with water vapor mixing ratio.
            ! We noticed this is needed for ASEX A209, particularly if Latin
            ! hypercube sampling is enabled.  -dschanen 11 Nov 2010
            call fill_holes_wv( nz, dt, exner, hydromet_name, & ! Intent(in)
                                rvm_mc, thlm_mc, hydromet(:,i) )   ! Intent(out)

         endif ! Variable is a mixing ratio and l_hole_fill is true

      endif ! hydromet(:,i) < 0

      ! Enter the new value of the hydrometeor for the effect of the water vapor
      ! hole-filling scheme.
      if ( stats_metadata%l_stats_samp ) then
         call stat_end_update( gr%nz, ixrm_wvhf, hydromet(:,i) / dt, & ! intent(in)
                               stats_zt ) ! intent(inout)
      endif

      ! Clipping for hydrometeor mixing ratios.
      if ( l_mix_rat_hm(i) ) then

         ! Store the previous value of the hydrometeor for the effect of
         ! clipping.
         if ( stats_metadata%l_stats_samp ) then
            call stat_begin_update( gr%nz, ixrm_cl, hydromet(:,i) / dt, & ! intent(in)
                                    stats_zt ) ! intent(inout)
         endif

         if ( any( hydromet(:,i) < zero_threshold ) ) then

            ! Clip any remaining negative values of precipitating hydrometeor
            ! mixing ratios to 0.
            where ( hydromet(:,i) < zero_threshold )
               hydromet(:,i) = zero_threshold
            end where

         endif ! hydromet(:,i) < 0

         ! Eliminate very small values of mean precipitating hydrometeor mixing
         ! ratios by setting them to 0.
         do k = 2, gr%nz, 1

            if ( hydromet(k,i) <= hydromet_tol(i) ) then

               rvm_mc(k) &
               = rvm_mc(k) &
                 + ( hydromet(k,i) / dt )

               if ( .not. l_frozen_hm(i) ) then

                  ! Rain water mixing ratio
   
                  thlm_mc(k) &
                  = thlm_mc(k) &
                    - ( Lv / ( Cp * exner(k) ) ) &
                      * ( hydromet(k,i) / dt )

               else ! Frozen hydrometeor mixing ratio

                  thlm_mc(k) &
                  = thlm_mc(k) &
                    - ( Ls / ( Cp * exner(k) ) ) &
                      * ( hydromet(k,i) / dt )

               endif ! l_frozen_hm(i)

               hydromet(k,i) = zero

            endif ! hydromet(k,i) <= hydromet_tol(i)

         enddo ! k = 2, gr%nz, 1


         ! Enter the new value of the hydrometeor for the effect of clipping.
         if ( stats_metadata%l_stats_samp ) then
            call stat_end_update( gr%nz, ixrm_cl, hydromet(:,i) / dt, & ! intent(in)
                                  stats_zt ) ! intent(inout)
         endif

      endif ! l_mix_rat_hm(i)

    enddo ! i = 1, hydromet_dim, 1

    ! Calculate clipping for hydrometeor concentrations.
    call clip_hydromet_conc_mvr( gr%nz, hydromet_dim, hydromet, & ! Intent(in)
                                 hydromet_clipped )        ! Intent(out)

    ! Clip hydrometeor concentrations and output stats.
    do i = 1, hydromet_dim

       if ( .not. l_mix_rat_hm(i) ) then

          if ( stats_metadata%l_stats_samp ) then

             ! Set up the stats indices for hydrometeor at index i
             call setup_stats_indices( i, stats_metadata,           & ! In
                                       ixrm_bt, ixrm_hf, ixrm_wvhf, & ! In/Out
                                       ixrm_cl, ixrm_mc,            & ! In/Out
                                       max_velocity )                 ! In/Out

             ! Store the previous value of the hydrometeor for the effect of
             ! clipping.
             call stat_begin_update( gr%nz, ixrm_cl, hydromet(:,i) / dt, & ! intent(in)
                                     stats_zt ) ! intent(inout)

          endif ! stats_metadata%l_stats_samp

          ! Apply clipping of hydrometeor concentrations.
          hydromet(:,i) = hydromet_clipped(:,i)

          ! Enter the new value of the hydrometeor for the effect of clipping.
          if ( stats_metadata%l_stats_samp ) then
             call stat_end_update( gr%nz, ixrm_cl, hydromet(:,i) / dt, & ! intent(in)
                                   stats_zt ) ! intent(inout)
          endif

       endif ! .not. l_mix_rat_hm(i)

    enddo ! i = 1, hydromet_dim, 1


    return

  end subroutine fill_holes_driver

  !=============================================================================
  subroutine clip_hydromet_conc_mvr( nz, hydromet_dim, hydromet, & ! Intent(in)
                                     hydromet_clipped )        ! Intent(out)

    ! Description:
    ! Increases the value of a hydrometeor concentration when it is too small,
    ! according to the hydrometeor mixing ratio (which remains unchanged) and
    ! the maximum acceptable drop or particle mean volume radius for that
    ! hydrometeor species.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        grid ! Type

    use constants_clubb, only: &
        pi,          & ! Variable(s)
        four_thirds, &
        one,         &
        zero,        &
        rho_lw,      &
        rho_ice

    use array_index, only: &
        l_mix_rat_hm, & ! Variable(s)
        l_frozen_hm

    use index_mapping, only: &
        Nx2rx_hm_idx, & ! Procedure(s)
        mvr_hm_max

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz
    
    integer, intent(in) :: &
      hydromet_dim    ! Number of hydrometeor fields

    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(in) :: &
      hydromet    ! Mean of hydrometeor fields    [units vary]

    ! Output Variable
    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(out) :: &
      hydromet_clipped    ! Clipped mean of hydrometeor fields    [units vary]

    ! Local Variables
    real( kind = core_rknd ) :: &
      Nxm_min_coef    ! Coefficient for min mean value of a concentration [1/kg]

    integer :: &
      k,   & ! Vertical grid index
      idx    ! Hydrometeor species index


    ! Clipping for hydrometeor concentrations.
    do idx = 1, hydromet_dim

       if ( .not. l_mix_rat_hm(idx) ) then

          ! The field is a hydrometeor concentration.

          if ( .not. l_frozen_hm(idx) ) then

             ! Clipping for mean rain drop concentration, <Nr>.
             !
             ! When mean rain water mixing ratio, <rr>, is found at a grid
             ! level, mean rain drop concentration must be at least a minimum
             ! value so that average rain drop mean volume radius stays within
             ! an upper bound.  Otherwise, mean rain drop concentration is 0.
             ! This can also be applied to values at a sample point, rather than
             ! a grid-box mean.

             ! The minimum mean rain drop concentration is given by:
             !
             ! <Nr> = <rr> / ( (4/3) * pi * rho_lw * mvr_rain_max^3 );
             !
             ! where mvr_rain_max is the maximum acceptable average rain drop
             ! mean volume radius.

             Nxm_min_coef &
             = one / ( four_thirds * pi * rho_lw * mvr_hm_max(idx)**3 )

          else ! l_frozen_hm(idx)

             ! Clipping for mean frozen hydrometeor concentration, <Nx>, where x
             ! stands for any frozen hydrometeor species.
             !
             ! When mean frozen hydrometeor mixing ratio, <rx>, is found at a
             ! grid level, mean frozen hydrometeor concentration must be at
             ! least a minimum value so that average frozen hydrometeor mean
             ! volume radius stays within an upper bound.  Otherwise, mean
             ! frozen hydrometeor concentration is 0.  This can also be applied
             ! to values at a sample point, rather than a grid-box mean.

             ! The minimum mean frozen hydrometeor concentration is given by:
             !
             ! <Nx> = <rx> / ( (4/3) * pi * rho_ice * mvr_x_max^3 );
             !
             ! where mvr_x_max is the maximum acceptable average frozen
             ! hydrometeor mean volume radius for frozen hydrometeor species, x.

             Nxm_min_coef &
             = one / ( four_thirds * pi * rho_ice * mvr_hm_max(idx)**3 )

          endif ! .not. l_frozen_hm(idx)

          ! Loop over vertical levels and increase hydrometeor concentrations
          ! when necessary.
          do k = 2, nz, 1

             if ( hydromet(k,Nx2rx_hm_idx(idx)) > zero ) then

                ! Hydrometeor mixing ratio, <rx>, is found at the grid level.
                hydromet_clipped(k,idx) &
                = max( hydromet(k,idx), &
                       Nxm_min_coef * hydromet(k,Nx2rx_hm_idx(idx)) )

             else ! <rx> = 0

                hydromet_clipped(k,idx) = zero

             endif ! hydromet(k,Nx2rx_hm_idx(idx)) > 0

          enddo ! k = 2, nz, 1

          ! The lowest thermodynamic level is below the model's lower boundary.
          hydromet_clipped(1,idx) = hydromet(1,idx)

       else ! l_mix_rat_hm(idx)

          ! The field is a hydrometeor mixing ratio.
          hydromet_clipped(:,idx) = hydromet(:,idx)

       endif ! .not. l_mix_rat_hm(idx)

    enddo ! idx = 1, hydromet_dim, 1


    return

  end subroutine clip_hydromet_conc_mvr

  !-----------------------------------------------------------------------
  subroutine setup_stats_indices( ihm, stats_metadata,         & ! Intent(in)
                                  ixrm_bt, ixrm_hf, ixrm_wvhf, & ! Intent(inout)
                                  ixrm_cl, ixrm_mc,            & ! Intent(inout)
                                  max_velocity )                 ! Intent(inout)

  ! Description:
  !
  ! Determines the stats output indices depending on the hydrometeor.

  ! Attention: hydromet_list needs to be set up before this routine is called.
  !
  ! Bogus example
  ! References:
  !
  ! None
  !-----------------------------------------------------------------------


    use array_index, only: &
        hydromet_list  ! Names of the hydrometeor species

    use clubb_precision, only: &
        core_rknd

    use constants_clubb, only: &
        zero

    use stats_variables, only: & 
        stats_metadata_type

    implicit none

    ! Input Variables
    integer, intent(in) :: ihm

    type (stats_metadata_type), intent(in) :: &
      stats_metadata

    ! Input/Output Variables
    real( kind = core_rknd ), intent(inout) :: &
      max_velocity ! Maximum sedimentation velocity [m/s]

    integer, intent(inout) :: ixrm_hf, ixrm_wvhf, ixrm_cl, &
                              ixrm_bt, ixrm_mc

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    ! Initializing max_velocity in order to avoid a compiler warning.
    ! Regardless of the case, it will be reset in the 'select case'
    ! statement immediately below.
    max_velocity = zero

    select case ( trim( hydromet_list(ihm) ) )
      case ( "rrm" )
        ixrm_bt   = stats_metadata%irrm_bt
        ixrm_hf   = stats_metadata%irrm_hf
        ixrm_wvhf = stats_metadata%irrm_wvhf
        ixrm_cl   = stats_metadata%irrm_cl
        ixrm_mc   = stats_metadata%irrm_mc

        max_velocity = -9.1_core_rknd ! m/s

      case ( "rim" )
        ixrm_bt   = stats_metadata%irim_bt
        ixrm_hf   = stats_metadata%irim_hf
        ixrm_wvhf = stats_metadata%irim_wvhf
        ixrm_cl   = stats_metadata%irim_cl
        ixrm_mc   = stats_metadata%irim_mc

        max_velocity = -1.2_core_rknd ! m/s

      case ( "rsm" )
        ixrm_bt   = stats_metadata%irsm_bt
        ixrm_hf   = stats_metadata%irsm_hf
        ixrm_wvhf = stats_metadata%irsm_wvhf
        ixrm_cl   = stats_metadata%irsm_cl
        ixrm_mc   = stats_metadata%irsm_mc

        ! Morrison limit
!         max_velocity = -1.2_core_rknd ! m/s
        ! Made up limit.  The literature suggests that it is quite possible
        ! that snow flake might achieve a terminal velocity of 2 m/s, and this
        ! happens in the COAMPS microphysics -dschanen 29 Sept 2009
        max_velocity = -2.0_core_rknd ! m/s

      case ( "rgm" )
        ixrm_bt   = stats_metadata%irgm_bt
        ixrm_hf   = stats_metadata%irgm_hf
        ixrm_wvhf = stats_metadata%irgm_wvhf
        ixrm_cl   = stats_metadata%irgm_cl
        ixrm_mc   = stats_metadata%irgm_mc

        max_velocity = -20._core_rknd ! m/s

      case ( "Nrm" )
        ixrm_bt   = stats_metadata%iNrm_bt
        ixrm_hf   = 0
        ixrm_wvhf = 0
        ixrm_cl   = stats_metadata%iNrm_cl
        ixrm_mc   = stats_metadata%iNrm_mc

        max_velocity = -9.1_core_rknd ! m/s

      case ( "Nim" )
        ixrm_bt   = stats_metadata%iNim_bt
        ixrm_hf   = 0
        ixrm_wvhf = 0
        ixrm_cl   = stats_metadata%iNim_cl
        ixrm_mc   = stats_metadata%iNim_mc

        max_velocity = -1.2_core_rknd ! m/s

      case ( "Nsm" )
        ixrm_bt   = stats_metadata%iNsm_bt
        ixrm_hf   = 0
        ixrm_wvhf = 0
        ixrm_cl   = stats_metadata%iNsm_cl
        ixrm_mc   = stats_metadata%iNsm_mc

        ! Morrison limit
!         max_velocity = -1.2_core_rknd ! m/s
        ! Made up limit.  The literature suggests that it is quite possible
        ! that snow flake might achieve a terminal velocity of 2 m/s, and this
        ! happens in the COAMPS microphysics -dschanen 29 Sept 2009
        max_velocity = -2.0_core_rknd ! m/s

      case ( "Ngm" )
        ixrm_bt   = stats_metadata%iNgm_bt
        ixrm_hf   = 0
        ixrm_wvhf = 0
        ixrm_cl   = stats_metadata%iNgm_cl
        ixrm_mc   = stats_metadata%iNgm_mc

        max_velocity = -20._core_rknd ! m/s

      case ( "Ncm" )
        ixrm_bt   = stats_metadata%iNcm_bt
        ixrm_hf   = 0
        ixrm_wvhf = 0
        ixrm_cl   = stats_metadata%iNcm_cl
        ixrm_mc   = stats_metadata%iNcm_mc

        ! Use the rain water limit, since Morrison has no explicit limit on
        ! cloud water.  Presumably these numbers are never large.
        ! -dschanen 28 Sept 2009
        max_velocity = -9.1_core_rknd ! m/s

      case default
        ixrm_bt   = 0
        ixrm_hf   = 0
        ixrm_wvhf = 0
        ixrm_cl   = 0
        ixrm_mc   = 0

        max_velocity = -9.1_core_rknd ! m/s

    end select


    return

  end subroutine setup_stats_indices
  !-----------------------------------------------------------------------

end module fill_holes
