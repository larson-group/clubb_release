!-------------------------------------------------------------------------
! $Id$
!===============================================================================
module precipitation_fraction

  ! Description:
  ! Sets overall precipitation fraction as well as the precipitation fraction
  ! in each PDF component.

  implicit none

  private

  public :: precip_fraction

  private :: component_precip_frac_weighted, &
             component_precip_frac_specify,  &
             component_precip_frac_ratio

  integer, parameter, public :: &
    precip_frac_calc_type = 1  ! Option used to calculate component precip_frac

  contains

  !=============================================================================
  subroutine precip_fraction( nz, hydromet, cloud_frac, cloud_frac_1, &
                              ice_supersat_frac, ice_supersat_frac_1, &
                              mixt_frac, rho, rc_1, rc_2, l_stats_samp, &
                              precip_frac, precip_frac_1, precip_frac_2 )

    ! Description:
    ! Determines (overall) precipitation fraction over the horizontal domain, as
    ! well as the precipitation fraction within each PDF component, at every
    ! vertical grid level.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one,            & ! Constant(s)
        zero,           &
        cloud_frac_min, &
        fstderr

    use parameters_model, only: &
        hydromet_dim  ! Variable(s)

    use array_index, only: &
        l_mix_rat_hm, & ! Variable(s)
        l_frozen_hm,  &
        hydromet_tol

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz          ! Number of model vertical grid levels

    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(in) :: &
      hydromet    ! Mean of hydrometeor, hm (overall)           [units vary]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      cloud_frac,          & ! Cloud fraction (overall)                      [-]
      cloud_frac_1,        & ! Cloud fraction (1st PDF component)            [-]
      ice_supersat_frac,   & ! Ice supersaturation fraction                  [-]
      ice_supersat_frac_1, & ! Ice supersaturation fraction (1st PDF comp.)  [-]
      mixt_frac,           & ! Mixture fraction                              [-]
      rho,                 & ! Air density                              [kg/m^3]
      rc_1,                & ! Mean cloud wat. mix. rat. (1st PDF comp.) [kg/kg]
      rc_2                   ! Mean cloud wat. mix. rat. (2nd PDF comp.) [kg/kg]

    logical, intent(in) :: &
      l_stats_samp     ! Flag to record statistical output.

    ! Output Variables
    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      precip_frac,   & ! Precipitation fraction (overall)               [-]
      precip_frac_1, & ! Precipitation fraction (1st PDF component)     [-]
      precip_frac_2    ! Precipitation fraction (2nd PDF component)     [-]

    ! Local Variables
    real( kind = core_rknd ) :: &
      precip_frac_tol    ! Minimum precip. frac. when hydromet. are present  [-]

    ! "Maximum allowable" hydrometeor mixing ratio in-precip component mean.
    real( kind = core_rknd ), parameter :: &
      max_hm_ip_comp_mean = 0.0025_core_rknd  ! [kg/kg]

    integer :: &
      k, ivar   ! Loop indices


    ! Initialize the precipitation fraction variables (precip_frac,
    ! precip_frac_1, and precip_frac_2) to 0.
    precip_frac   = zero
    precip_frac_1 = zero
    precip_frac_2 = zero

    ! Set the minimum allowable precipitation fraction when hydrometeors are
    ! found at a grid level.
    precip_frac_tol = max( 0.1_core_rknd * maxval( cloud_frac ), &
                           cloud_frac_min )

    !!! Find overall precipitation fraction.
    do k = nz, 1, -1

       ! The precipitation fraction is the greatest cloud fraction at or above a
       ! vertical level.
       if ( k < nz ) then
          if ( any( l_frozen_hm ) ) then
             ! Ice microphysics included.
             precip_frac(k) = max( precip_frac(k+1), cloud_frac(k), &
                                   ice_supersat_frac(k) )
          else
             ! Warm microphysics.
             precip_frac(k) = max( precip_frac(k+1), cloud_frac(k) )
          endif
       else  ! k = nz
          if ( any( l_frozen_hm ) ) then
             ! Ice microphysics included.
             precip_frac(k) = max( cloud_frac(k), ice_supersat_frac(k) )
          else
             ! Warm microphysics.
             precip_frac(k) = cloud_frac(k)
          endif
       endif

       if ( any( hydromet(k,:) >= hydromet_tol(:) ) &
            .and. precip_frac(k) < precip_frac_tol ) then

          ! In a scenario where we find any hydrometeor at this grid level, but
          ! no cloud at or above this grid level, set precipitation fraction to
          ! a minimum threshold value.
          precip_frac(k) = precip_frac_tol

       elseif ( all( hydromet(k,:) < hydromet_tol(:) ) ) then

          ! The means (overall) of every precipitating hydrometeor are all less
          ! than their respective tolerance amounts.  They are all considered to
          ! have values of 0.  There are not any hydrometeor species found at
          ! this grid level.  There is also no cloud at or above this grid
          ! level, so set precipitation fraction to 0.
          precip_frac(k) = zero

       endif

    enddo ! Overall precipitation fraction loop: k = nz, 1, -1.


    !!! Find precipitation fraction within each PDF component.
    !
    ! The overall precipitation fraction, f_p, is given by the equation:
    !
    ! f_p = a * f_p(1) + ( 1 - a ) * f_p(2);
    !
    ! where "a" is the mixture fraction (weight of PDF component 1), f_p(1) is
    ! the precipitation fraction within PDF component 1, and f_p(2) is the
    ! precipitation fraction within PDF component 2.  Overall precipitation
    ! fraction is found according the method above, and mixture fraction is
    ! already determined, leaving f_p(1) and f_p(2) to be solved for.  The
    ! values for f_p(1) and f_p(2) must satisfy the above equation.
    if ( precip_frac_calc_type == 1 ) then

       ! Calculatate precip_frac_1 and precip_frac_2 based on the greatest
       ! weighted cloud_frac_1 at or above a grid level.
       call component_precip_frac_weighted( nz, hydromet, &
                                            precip_frac, cloud_frac_1, &
                                            ice_supersat_frac_1, mixt_frac, &
                                            precip_frac_tol, &
                                            precip_frac_1, precip_frac_2 )

    elseif ( precip_frac_calc_type == 2 ) then

       ! Specified method.
       call component_precip_frac_specify( nz, hydromet, precip_frac, &
                                           mixt_frac, precip_frac_tol, &
                                           precip_frac_1, precip_frac_2 )

    elseif ( precip_frac_calc_type == 3 ) then

       ! Ratio method:  precip_frac_2 / precip_frac_1 = LWP_2 / LWP_1.
       call component_precip_frac_ratio( nz, hydromet, precip_frac, &
                                         rho, rc_1, rc_2, mixt_frac, &
                                         precip_frac_tol, l_stats_samp, &
                                         precip_frac_1, precip_frac_2 )

    else ! Invalid option selected.

       write(fstderr,*) "Invalid option to calculate precip_frac_1 " &
                        // "and precip_frac_2."
       stop

    endif ! precip_frac_calc_type


    ! Increase Precipiation Fraction under special conditions.
    !
    ! There are scenarios that sometimes occur that require precipitation
    ! fraction to be boosted.  Precipitation fraction is calculated from cloud
    ! fraction and ice supersaturation fraction.  For numerical reasons, CLUBB's
    ! PDF may become entirely subsaturated with respect to liquid and ice,
    ! resulting in both a cloud fraction of 0 and an ice supersaturation
    ! fraction of 0.  When this happens, precipitation fraction drops to 0 when
    ! there aren't any hydrometeors present at that grid level, or to
    ! precip_frac_tol when there is at least one hydrometeor present at that
    ! grid level.  However, sometimes there are large values of hydrometeors
    ! found at that grid level.  When this occurs, the PDF component in-precip
    ! mean of a hydrometeor can become ridiculously large.  This is because the
    ! ith PDF component in-precip mean of a hydrometeor, mu_hm_i,  is given by
    ! the equation:
    !
    ! mu_hm_i = hm_i / precip_frac_i;
    !
    ! where hm_i is the overall ith PDF component mean of the hydrometeor, and
    ! precip_frac_i is the ith PDF component precipitation fraction.  When
    ! precip_frac_i has a value of precip_frac_tol and hm_i is large, mu_hm_i
    ! can be huge.  This can cause enormous microphysical process rates and
    ! result in numerical instability.  It is also very inaccurate.
    !
    ! In order to limit this problem, the ith PDF component precipitation
    ! fraction is increased in order to decrease mu_hm_i.  First, an "upper
    ! limit" is set for mu_hm_i when the hydrometeor is a mixing ratio.  This is
    ! called max_hm_ip_comp_mean.  At every vertical level and for every
    ! hydrometeor mixing ratio, a check is made to try to prevent mu_hm_i from
    ! exceeding the "upper limit".  The check is:
    !
    ! hm_i / precip_frac_i ( which = mu_hm_i )  >  max_hm_ip_comp_mean,
    !
    ! which can be rewritten:
    !
    ! hm_i > precip_frac_i * max_hm_ip_comp_mean.
    !
    ! Since hm_i has not been calculated yet, the assumption for this check is
    ! that all of the hydrometeor is found in one PDF component, which is the
    ! worst-case scenario in violating this limit.  The check becomes:
    !
    ! <hm> / ( mixt_frac * precip_frac_1 )  >  max_hm_ip_comp_mean;
    !    in PDF comp. 1; and
    ! <hm> / ( ( 1 - mixt_frac ) * precip_frac_2 )  >  max_hm_ip_comp_mean;
    !    in PDF comp. 2.
    !
    ! These limits can be rewritten as:
    !
    ! <hm>  >  mixt_frac * precip_frac_1 * max_hm_ip_comp_mean;
    !    in PDF comp. 1; and
    ! <hm>  >  ( 1 - mixt_frac ) * precip_frac_2 * max_hm_ip_comp_mean;
    !    in PDF comp. 2.
    !
    ! When component precipitation fraction is found to be in excess of the
    ! limit, precip_frac_i is increased to:
    !
    ! <hm> / ( mixt_frac * max_hm_ip_comp_mean );
    !    when the limit is exceeded in PDF comp. 1; and
    ! <hm> / ( ( 1 - mixt_frac ) * max_hm_ip_comp_mean );
    !    when the limit is exceeded in PDF comp. 2.
    !
    ! Of course, precip_frac_i is not allowed to exceed 1, so when
    ! <hm> / mixt_frac (or <hm> / ( 1 - mixt_frac )) is already greater than
    ! max_hm_ip_comp_mean, mu_hm_i will also have to be greater than
    ! max_hm_ip_comp_mean.  However, the value of mu_hm_i is still reduced when
    ! compared to what it would have been using precip_frac_tol.  In the event
    ! that multiple hydrometeor mixing ratios violate the check, the code is set
    ! up so that precip_frac_i is increased based on the highest hm_i.
    do k = 1, nz, 1

       do ivar = 1, hydromet_dim, 1

          if ( l_mix_rat_hm(ivar) ) then

             ! The hydrometeor is a mixing ratio.

             if ( hydromet(k,ivar) > mixt_frac(k) * precip_frac_1(k) &
                                     * max_hm_ip_comp_mean ) then

                ! Increase precipitation fraction in the 1st PDF component.
                precip_frac_1(k) &
                = min( hydromet(k,ivar) &
                       / ( mixt_frac(k) * max_hm_ip_comp_mean ), one )

             endif ! <hm>/(mixt_frac*precip_frac_1) > max_hm_ip_comp_mean

             if ( hydromet(k,ivar) > ( one - mixt_frac(k) ) * precip_frac_2(k) &
                                     * max_hm_ip_comp_mean ) then

                ! Increase precipitation fraction in the 2nd PDF component.
                precip_frac_2(k) &
                = min( hydromet(k,ivar) &
                       / ( ( one - mixt_frac(k) ) * max_hm_ip_comp_mean ), one )

             endif ! <hm>/((1-mixt_frac)*precip_frac_2) > max_hm_ip_comp_mean

          endif ! l_mix_rat_hm(ivar)

       enddo ! ivar = 1, hydromet_dim, 1

    enddo ! k = 1, nz, 1

    ! Recalculate overall precipitation fraction for consistency.
    precip_frac = mixt_frac * precip_frac_1 &
                  + ( one - mixt_frac ) * precip_frac_2


    return

  end subroutine precip_fraction

  !=============================================================================
  subroutine component_precip_frac_weighted( nz, hydromet, &
                                             precip_frac, cloud_frac_1, &
                                             ice_supersat_frac_1, mixt_frac, &
                                             precip_frac_tol, &
                                             precip_frac_1, precip_frac_2 )

    ! Description:
    ! Set precipitation fraction in each component of the PDF.  Set
    ! Set precip_frac_1 at a grid level based on the greatest
    ! mixt_frac * cloud_frac_1 at or above the relevant grid level, divided by
    ! mixt_frac at the grid level.  Set precip_frac_2 accordingly.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one,  & ! Constant(s)
        zero

    use parameters_model, only: &
        hydromet_dim  ! Variable(s)

    use array_index, only: &
        l_frozen_hm,  & ! Variable(s)
        hydromet_tol

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz          ! Number of model vertical grid levels

    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(in) :: &
      hydromet    ! Mean of hydrometeor, hm (overall)           [units vary]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      precip_frac,         & ! Precipitation fraction (overall)              [-]
      cloud_frac_1,        & ! Cloud fraction (1st PDF component)            [-]
      ice_supersat_frac_1, & ! Ice supersaturation fraction (1st PDF comp.)  [-]
      mixt_frac              ! Mixture fraction                              [-]

    real( kind = core_rknd ), intent(in) :: &
      precip_frac_tol    ! Minimum precip. frac. when hydromet. are present  [-]

    ! Output Variables
    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      precip_frac_1, & ! Precipitation fraction (1st PDF component)     [-]
      precip_frac_2    ! Precipitation fraction (2nd PDF component)     [-]

    ! Local Variables
    real( kind = core_rknd ), dimension(nz) :: &
      weighted_pfrac_1    ! Product of mixt_frac and precip_frac_1      [-]

    integer :: k  ! Loop index


    !!! Find precipitation fraction within PDF component 1.
    ! The method used to find overall precipitation fraction will also be to
    ! find precipitation fraction within PDF component 1.  In order to do so, it
    ! is assumed (poorly) that PDF component 1 overlaps PDF component 1 at every
    ! vertical level in the vertical profile.
    do k = nz, 1, -1

       ! The weighted precipitation fraction (PDF component 1) is the greatest
       ! value of the product of mixture fraction and cloud fraction (PDF
       ! component 1) at or above a vertical level.
       if ( k < nz ) then
          if ( any( l_frozen_hm ) ) then
             ! Ice microphysics included.
             weighted_pfrac_1(k) = max( weighted_pfrac_1(k+1), &
                                        mixt_frac(k) * cloud_frac_1(k), &
                                        mixt_frac(k) * ice_supersat_frac_1(k) )
          else
             ! Warm microphysics.
             weighted_pfrac_1(k) = max( weighted_pfrac_1(k+1), &
                                        mixt_frac(k) * cloud_frac_1(k) )
          endif
       else  ! k = nz
          if ( any( l_frozen_hm ) ) then
             ! Ice microphysics included.
             weighted_pfrac_1(k) = max( mixt_frac(k) * cloud_frac_1(k), &
                                        mixt_frac(k) * ice_supersat_frac_1(k) )
          else
             ! Warm microphysics.
             weighted_pfrac_1(k) = mixt_frac(k) * cloud_frac_1(k)
          endif
       endif

       precip_frac_1(k) = weighted_pfrac_1(k) / mixt_frac(k)

       ! Special cases for precip_frac_1.
       if ( precip_frac_1(k) > one ) then

          ! Using the above method, it is possible for precip_frac_1 to be
          ! greater than 1.  For example, the mixture fraction at level k+1 is
          ! 0.10 and the cloud_frac_1 at level k+1 is 1, resulting in a
          ! weighted_pfrac_1 of 0.10.  This product is greater than the product
          ! of mixt_frac and cloud_frac_1 at level k.  The mixture fraction at
          ! level k is 0.05, resulting in a precip_frac_1 of 2.  The value of
          ! precip_frac_1 is limited at 1.  The leftover precipitation fraction
          ! (a result of the decreasing weight of PDF component 1 between the
          ! levels) is applied to PDF component 2.
          precip_frac_1(k) = one

       elseif ( precip_frac_1(k) > zero &
                .and. precip_frac_1(k) < precip_frac_tol ) then

          ! In a scenario where we find precipitation in the 1st PDF component
          ! but it is tiny (less than tolerance level), boost 1st PDF component
          ! precipitation fraction to tolerance level.
          precip_frac_1(k) = precip_frac_tol

       endif

    enddo ! Precipitation fraction (1st PDF component) loop: k = nz, 1, -1.


    !!! Find precipitation fraction within PDF component 2.
    ! The equation for precipitation fraction within PDF component 2 is:
    !
    ! f_p(2) = ( f_p - a * f_p(1) ) / ( 1 - a );
    !
    ! given the overall precipitation fraction, f_p (calculated above), the
    ! precipitation fraction within PDF component 1, f_p(1) (calculated above),
    ! and mixture fraction, a.  Any leftover precipitation fraction from
    ! precip_frac_1 will be included in this calculation of precip_frac_2.
    do k = 1, nz, 1

       precip_frac_2(k) &
       = ( precip_frac(k) - mixt_frac(k) * precip_frac_1(k) ) &
         / ( one - mixt_frac(k) )

       ! Special cases for precip_frac_2.
       if ( precip_frac_2(k) > one ) then

          ! Again, it is possible for precip_frac_2 to be greater than 1.  For
          ! example, the mixture fraction at level k+1 is 0.10 and the
          ! cloud_frac_1 at level k+1 is 1, resulting in a weighted_pfrac_1 of
          ! 0.10.  This product is greater than the product of mixt_frac and
          ! cloud_frac_1 at level k.  Additionally, precip_frac (overall) is 1
          ! for level k.  The mixture fraction at level k is 0.5, resulting in
          ! a precip_frac_1 of 0.2.  Using the above equation, precip_frac_2 is
          ! calculated to be 1.8.  The value of precip_frac_2 is limited at 1.
          ! The leftover precipitation fraction (as a result of the increasing
          ! weight of component 1 between the levels) is applied to PDF
          ! component 1.
          precip_frac_2(k) = one

          ! Recalculate the precipitation fraction in PDF component 1.
          precip_frac_1(k) &
          = ( precip_frac(k) - ( one - mixt_frac(k) ) ) / mixt_frac(k)

          ! Double check precip_frac_1
          if ( precip_frac_1(k) > zero &
               .and. precip_frac_1(k) < precip_frac_tol ) then
             precip_frac_1(k) = precip_frac_tol
             precip_frac_2(k) = ( precip_frac(k) &
                                  - mixt_frac(k) * precip_frac_1(k) ) &
                                / ( one - mixt_frac(k) )
          endif

       elseif ( precip_frac_2(k) > zero &
                .and. precip_frac_2(k) < precip_frac_tol ) then

          ! In a scenario where we find precipitation in the 2nd PDF component
          ! but it is tiny (less than tolerance level), boost 2nd PDF component
          ! precipitation fraction to tolerance level.
          precip_frac_2(k) = precip_frac_tol

          ! Recalculate the precipitation fraction in PDF component 1.
          precip_frac_1(k) &
          = ( precip_frac(k) - ( one - mixt_frac(k) ) * precip_frac_2(k) ) &
            / mixt_frac(k)

          ! Double check precip_frac_1
          if ( precip_frac_1(k) > one ) then
             precip_frac_1(k) = one
             precip_frac_2(k) = ( precip_frac(k) - mixt_frac(k) ) &
                                / ( one - mixt_frac(k) )
          endif

       endif

    enddo ! Precipitation fraction (2nd PDF component) loop: k = 1, nz, 1.


    ! When there aren't any hydrometeors found at a grid level, reset the
    ! component precipitation fractions to 0.
    do k = 1, nz, 1
       if ( all( hydromet(k,:) < hydromet_tol(:) ) ) then
          precip_frac_1(k) = zero
          precip_frac_2(k) = zero
       endif
    enddo ! k = 1, nz, 1


    return

  end subroutine component_precip_frac_weighted

  !=============================================================================
  subroutine component_precip_frac_specify( nz, hydromet, precip_frac, &
                                            mixt_frac, precip_frac_tol, &
                                            precip_frac_1, precip_frac_2 )

    ! Description:
    ! Calculates the precipitation fraction in each PDF component.
    !
    ! The equation for precipitation fraction is:
    !
    ! f_p = mixt_frac * f_p(1) + ( 1 - mixt_frac ) * f_p(2);
    !
    ! where f_p is overall precipitation fraction, f_p(1) is precipitation
    ! fraction in the 1st PDF component, f_p(2) is precipitation fraction in the
    ! 2nd PDF component, and mixt_frac is the mixture fraction.  Using this
    ! method, a new specified parameter is introduced, upsilon, where:
    !
    ! upsilon = mixt_frac * f_p(1) / f_p; and where 0 <= upsilon <= 1.
    !
    ! In other words, upsilon is the ratio of mixt_frac * f_p(1) to f_p.  Since
    ! f_p and mixt_frac are calculated previously, and upsilon is specified,
    ! f_p(1) can be calculated by:
    !
    ! f_p(1) = upsilon * f_p / mixt_frac;
    !
    ! and has an upper limit of 1.  The value of f_p(2) can then be calculated
    ! by:
    !
    ! f_p(2) = ( f_p - mixt_frac * f_p(1) ) / ( 1 - mixt_frac );
    !
    ! and also has an upper limit of 1.  When upsilon = 1, all of the
    ! precipitation is found in the 1st PDF component (as long as
    ! f_p <= mixt_frac, otherwise it would cause f_p(1) to be greater than 1).
    ! When upsilon = 0, all of the precipitation is found in the 2nd PDF
    ! component (as long as f_p <= 1 - mixt_frac, otherwise it would cause
    ! f_p(2) to be greater than 1).  When upsilon is between 0 and 1,
    ! precipitation is split between the two PDF components accordingly.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one,  & ! Constant(s)
        zero

    use parameters_tunable, only: &
        upsilon_precip_frac_rat  ! Variable(s)

    use parameters_model, only: &
        hydromet_dim  ! Variable(s)

    use array_index, only: &
        hydromet_tol  ! Variable(s)

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz          ! Number of model vertical grid levels

    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(in) :: &
      hydromet    ! Mean of hydrometeor, hm (overall)           [units vary]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      precip_frac, & ! Precipitation fraction (overall)                      [-]
      mixt_frac      ! Mixture fraction                                      [-]

    real( kind = core_rknd ), intent(in) :: &
      precip_frac_tol    ! Minimum precip. frac. when hydromet. are present  [-]

    ! Output Variables
    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      precip_frac_1, & ! Precipitation fraction (1st PDF component)     [-]
      precip_frac_2    ! Precipitation fraction (2nd PDF component)     [-]

    integer :: k  ! Loop index.


    ! Loop over all vertical levels.
    do k = 1, nz, 1

       if ( any( hydromet(k,:) >= hydromet_tol(:) ) ) then

          ! There are hydrometeors found at this grid level.
          if ( upsilon_precip_frac_rat == one ) then

             if ( precip_frac(k) <= mixt_frac(k) ) then
                ! All the precipitation is found in the 1st PDF component.
                precip_frac_1(k) = precip_frac(k) / mixt_frac(k)
                precip_frac_2(k) = zero
             else ! precip_frac(k) > mixt_frac(k)
                ! Some precipitation is found in the 2nd PDF component.
                precip_frac_1(k) = one
                precip_frac_2(k) = ( precip_frac(k) - mixt_frac(k) ) &
                                   / ( one - mixt_frac(k) )
                if ( precip_frac_2(k) < precip_frac_tol ) then
                   ! Since precipitation is found in the 2nd PDF component, it
                   ! must have a value of at least precip_frac_tol.
                   precip_frac_2(k) = precip_frac_tol
                   ! Recalculate precip_frac_1
                   precip_frac_1(k) &
                   = ( precip_frac(k) &
                       - ( one - mixt_frac(k) ) * precip_frac_2(k) ) &
                     / mixt_frac(k)
                endif ! precip_frac_2(k) < precip_frac_tol
             endif ! precip_frac(k) <= mixt_frac(k)

          elseif ( upsilon_precip_frac_rat == zero ) then

             if ( precip_frac(k) <= ( one - mixt_frac(k) ) ) then
                ! All the precipitation is found in the 2nd PDF component.
                precip_frac_1(k) = zero
                precip_frac_2(k) = precip_frac(k) / ( one - mixt_frac(k) )
             else ! precip_frac(k) > ( 1 - mixt_frac(k) )
                ! Some precipitation is found in the 1st PDF component.
                precip_frac_1(k) = ( precip_frac(k) - ( one - mixt_frac(k) ) ) &
                                   / mixt_frac(k)
                precip_frac_2(k) = one
                if ( precip_frac_1(k) < precip_frac_tol ) then
                   ! Since precipitation is found in the 1st PDF component, it
                   ! must have a value of at least precip_frac_tol.
                   precip_frac_1(k) = precip_frac_tol
                   ! Recalculate precip_frac_2
                   precip_frac_2(k) = ( precip_frac(k) &
                                        - mixt_frac(k) * precip_frac_1(k) ) &
                                      / ( one - mixt_frac(k) )
                endif ! precip_frac_1(k) < precip_frac_tol
             endif ! precip_frac(k) <= ( 1 - mixt_frac(k) )

          else  ! 0 < upsilon_precip_frac_rat < 1

             ! Precipitation is found in both PDF components.  Each component
             ! must have a precipitation fraction that is at least
             ! precip_frac_tol and that does not exceed 1.
             precip_frac_1(k) &
             = upsilon_precip_frac_rat * precip_frac(k) / mixt_frac(k)

             if ( precip_frac_1(k) > one ) then
                precip_frac_1(k) = one
             elseif ( precip_frac_1(k) < precip_frac_tol ) then
                precip_frac_1(k) = precip_frac_tol
             endif

             precip_frac_2(k) = ( precip_frac(k) &
                                  - mixt_frac(k) * precip_frac_1(k) ) &
                                / ( one - mixt_frac(k) )

             if ( precip_frac_2(k) > one ) then
                precip_frac_2(k) = one
                ! Recalculate precip_frac_1
                precip_frac_1(k) = ( precip_frac(k) - ( one - mixt_frac(k) ) ) &
                                   / mixt_frac(k)
             elseif ( precip_frac_2(k) < precip_frac_tol ) then
                precip_frac_2(k) = precip_frac_tol
                ! Recalculate precip_frac_1
                precip_frac_1(k) &
                = ( precip_frac(k) &
                    - ( one - mixt_frac(k) ) * precip_frac_2(k) ) &
                  / mixt_frac(k)
             endif

          endif  ! upsilon_precip_frac_rat


       else ! all( hydromet(k,:) < hydromet_tol(:) )

          ! There aren't any hydrometeors found at the grid level.
          precip_frac_1(k) = zero
          precip_frac_2(k) = zero


       endif ! any( hydromet(k,:) >= hydromet_tol(:) )

    enddo ! k = 1, nz, 1


    return

  end subroutine component_precip_frac_specify

  !=============================================================================
  subroutine component_precip_frac_ratio( nz, hydromet, precip_frac, &
                                          rho, rc_1, rc_2, mixt_frac, &
                                          precip_frac_tol, l_stats_samp, &
                                          precip_frac_1, precip_frac_2 )

    ! Description:
    ! Set precipitation fraction in each component of the PDF.  The
    ! precipitation fractions in each PDF component, f_p(1) and f_p(2), are
    ! related to the overall precipitation fraction, f_p, by:
    !
    ! f_p = a f_p(1) + (1-a) f_p(2);
    !
    ! where "a" is mixture fraction, which is the relative weight of the 1st PDF
    ! component.  This can be rewritten in terms of the ratio f_p(2)/f_p(1):
    !
    ! f_p = f_p(1) ( a + (1-a) f_p(2)/f_p(1) ).
    !
    ! At a grid level that is at least mostly cloudy, the simplest way to handle
    ! the ratio f_p(2)/f_p(1) is to set it equal to the ratio rc_2/rc_1, where
    ! rc_1 is the mean cloud water mixing ratio in PDF component 1 and rc_2 is
    ! the mean cloud water mixing ratio in PDF component 2.  However, a
    ! precipitating hydrometeor sediments, falling from higher altitudes
    ! downwards.  The values of cloud water mixing ratio at a given grid level
    ! are not necessarily indicative of the amount of cloud water at higher
    ! levels.  A precipitating hydrometeor may have been already produced from
    ! cloud water at a higher altitude (vertical level) and fallen downwards to
    ! the given grid level.  Additionally, using grid-level cloud water mixing
    ! ratio especially does not work for a precipitating hydrometeor below cloud
    ! base (near the ground).
    !
    ! However, an alternative to component cloud water mixing ratio is component
    ! liquid water path.  Liquid water path accounts for the cloud water mixing
    ! ratio at the given grid level and at all grid levels higher in altitude.
    !
    ! In a stratocumulus case, the cloud water is spread out over all or almost
    ! all of the horizontal domain over a group of vertical levels.  At a given
    ! vertical level, the component mean cloud water mixing ratios should be
    ! almost equal, although usually slightly larger in the component with the
    ! larger component mean extended liquid water mixing ratio, chi.  Likewise,
    ! the component liquid water paths should be nearly equal, with one
    ! component having a slightly larger liquid water path than the other
    ! component.
    !
    ! In a case of cumulus rising into stratocumulus, the upper portion of the
    ! cloudy domain will be very similar to the stratocumulus case described
    ! above, with similar cloud water mixing ratio and liquid water path
    ! results.  However, below the base of the stratocumulus clouds, where the
    ! cumulus clouds are found, the horizontal domain at each vertical level is
    ! only partially cloudy.  At these levels, any precipitating hydrometeor
    ! that was produced in the stratocumulus clouds above and fallen downwards
    ! is evaporating in the clear-air portions, while not evaporating in the
    ! cloudy portions.  Additionally, new amounts of a hydrometeor are being
    ! produced in the cloudy portions.  The amount of a hydrometeor in the
    ! cloudy portions becomes significantly larger than the amount of a
    ! hydrometeor in the clear portions.  The partially cloudy levels usually
    ! have a PDF where one component is significantly more saturated than the
    ! other component.  By the time the cloud base of the cumulus clouds is
    ! reached, the liquid water path for one PDF component should be
    ! significantly greater than the liquid water path for the other PDF
    ! component.
    !
    ! In a cumulus case, the horizontal domain at each level is usually partly
    ! cloudy.  Throughout the entire vertical domain, at every vertical level,
    ! one component usually is much more saturated than the other component.
    ! The liquid water path for one component is much greater than the liquid
    ! water path in the other component.  Likewise, a precipitating hydrometeor
    ! that is formed in cloud and falls preferentially through cloud will have
    ! large values in a portion of the horizontal domain and very small or 0
    ! values over the rest of the horizontal domain.
    !
    ! In order to estimate the precipitation fraction in each PDF component,
    ! the ratio f_p(2)/f_p(1) is going to be set equal to the ratio LWP_2/LWP_1,
    ! where LWP_1 is the liquid water path in PDF component 1 and LWP_2 is the
    ! liquid water path in PDF component 2.  LWP_1 will be computed by taking
    ! the vertical integral of cloud water (see equation below) through the 1st
    ! PDF component from the given vertical level all the way to the top of the
    ! model.  LWP_2 will be computed in the same manner.   It should be noted
    ! that this method makes the poor assumption that PDF component 1 always
    ! overlaps PDF component 1 between vertical levels, and likewise for PDF
    ! component 2.
    !
    ! Total liquid water path, LWP, is given by the following equation:
    !
    ! LWP(z) = INT(z:z_top) rho_a <r_c> dz';
    !
    ! where z is the altitude of the vertical level for which LWP is desired,
    ! z_top is the altitude at the top of the model domain, and z' is the
    ! dummy variable of integration.  Mean cloud water mixing ratio can be
    ! written as:
    !
    ! <r_c> = a * rc_1 + (1-a) * rc_2.
    !
    ! The equation for liquid water path is rewritten as:
    !
    ! LWP(z) = INT(z:z_top) rho_a ( a rc_1 + (1-a) rc_2 ) dz'; or
    !
    ! LWP(z) = INT(z:z_top) a rho_a rc_1 dz'
    !          + INT(z:z_top) (1-a) rho_a rc_2 dz'.
    !
    ! This can be rewritten as:
    !
    ! LWP(z) = LWP_1(z) + LWP_2(z);
    !
    ! where:
    !
    ! LWP_1(z) = INT(z:z_top) a rho_a rc_1 dz'; and
    ! LWP_2(z) = INT(z:z_top) (1-a) rho_a rc_2 dz'.
    !
    ! The trapezoidal rule will be used to numerically integrate for LWP_1
    ! and LWP_2.
    !
    ! As stated above, precipitation fraction in each PDF component is based on
    ! the liquid water path in each PDF component, such that:
    !
    ! f_p(2)/f_p(1) = LWP_2/LWP_1.
    !
    ! Substituting the ratio LWP_2/LWP_1 for the ratio f_p(2)/f_p(1), the
    ! above equation can be solved for f_p(1):
    !
    ! f_p(1) = f_p / ( a + (1-a) LWP_2/LWP_1 ).
    !
    ! Then, f_p(2) can be solved for according to the equation:
    !
    ! f_p(2) = ( f_p - a f_p(1) ) / (1-a).

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        gr    ! Variable(s)

    use constants_clubb, only: &
        one,      & ! Constant(s)
        one_half, &
        zero

    use parameters_model, only: &
        hydromet_dim  ! Variable(s)

    use array_index, only: &
        hydromet_tol  ! Variable(s)

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    use stats_type_utilities, only: &
        stat_update_var  ! Procedure(s)

    use stats_variables, only : &
        iLWP1, & ! Variable(s)
        iLWP2, &
        stats_zt

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz    ! Number of model vertical grid levels

    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(in) :: &
      hydromet    ! Mean of hydrometeor, hm (overall)           [units vary]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      precip_frac, & ! Precipitation fraction (overall)               [-]
      rho,         & ! Air density                                    [kg/m^3]
      rc_1,        & ! Mean cloud water mixing ratio (1st PDF comp.)  [kg/kg]
      rc_2,        & ! Mean cloud water mixing ratio (2nd PDF comp.)  [kg/kg]
      mixt_frac      ! Mixture fraction                               [-]

    real( kind = core_rknd ), intent(in) :: &
      precip_frac_tol    ! Minimum precip. frac. when hydromet. are present  [-]

    logical, intent(in) :: &
      l_stats_samp     ! Flag to record statistical output.

    ! Output Variables
    real( kind = core_rknd ), dimension(nz), intent(out) :: &
      precip_frac_1, & ! Precipitation fraction (1st PDF component)     [-]
      precip_frac_2    ! Precipitation fraction (2nd PDF component)     [-]

    ! Local Variable
    real( kind = core_rknd ), dimension(nz) :: &
      LWP_1, & ! Liquid water path (1st PDF component) on thermo. levs. [kg/m^2]
      LWP_2    ! Liquid water path (2nd PDF component) on thermo. levs. [kg/m^2]

    integer :: k  ! Array index

    real( kind = core_rknd ), parameter :: &
      LWP_tol = 5.0e-7_core_rknd  ! Tolerance value for component LWP


    !!! Compute component liquid water paths using trapezoidal rule for
    !!! numerical integration.

    ! At the uppermost thermodynamic level (k = nz), use the trapezoidal rule:
    !
    ! 0.5 * (integrand_a + integrand_b) * delta_z,
    !
    ! where integrand_a is the integrand at thermodynamic level k = nz,
    ! integrand_b is the integrand at momentum level k = nz (model upper
    ! boundary), and delta_z = zm(nz) - zt(nz).  At the upper boundary, r_c is
    ! set to 0, and the form of the trapezoidal rule is simply:
    !
    ! 0.5 * integrand_a * delta_z.

    ! Liquid water path in PDF component 1.
    LWP_1(nz) &
    = one_half * mixt_frac(nz) * rho(nz) * rc_1(nz) * ( gr%zm(nz) - gr%zt(nz) )

    ! Liquid water path in PDF component 2.
    LWP_2(nz) &
    = one_half * ( one - mixt_frac(nz) ) * rho(nz) * rc_2(nz) &
      * ( gr%zm(nz) - gr%zt(nz) )

    ! At all other thermodynamic levels, compute liquid water path using the
    ! trapezoidal rule:
    !
    ! 0.5 * (integrand_a + integrand_b) * delta_z,
    !
    ! where integrand_a is the integrand at thermodynamic level k, integrand_b
    ! is the integrand at thermodynamic level k+1, and
    ! delta_z = zt(k+1) - zt(k), or 1/invrs_dzm(k).  The total for the segment
    ! is added to the sum total of all higher vertical segments to compute the
    ! total vertical integral.
    do k = nz-1, 1, -1

       ! Liquid water path in PDF component 1.
       LWP_1(k) &
       = LWP_1(k+1) &
         + one_half * ( mixt_frac(k+1) * rho(k+1) * rc_1(k+1) &
                        + mixt_frac(k) * rho(k) * rc_1(k) ) / gr%invrs_dzm(k)

       ! Liquid water path in PDF component 2.
       LWP_2(k) &
       = LWP_2(k+1) &
         + one_half * ( ( one - mixt_frac(k+1) ) * rho(k+1) * rc_2(k+1) &
                        + ( one - mixt_frac(k) ) * rho(k) * rc_2(k) ) &
           / gr%invrs_dzm(k)

    enddo ! k = nz-1, 1, -1


    !!! Find precip_frac_1 and precip_frac_2 based on the ratio of LWP_2/LWP_1,
    !!! such that:  precip_frac_2/precip_frac_1 = LWP_2/LWP_1.
    do k = 1, nz, 1

       !!! Calculate the component precipitation fractions.
       if ( any( hydromet(k,:) >= hydromet_tol(:) ) ) then

          if ( LWP_1(k) < LWP_tol .and. LWP_2(k) < LWP_tol ) then

             ! Both LWP_1 and LWP_2 are 0 (or an insignificant amount).
             !
             ! Precipitation is found at this level, yet there is no cloud at or
             ! above the current level.  This is usually due to a numerical
             ! artifact.  For example, the hydrometeor is diffused above cloud
             ! top.  Simply set each component precipitation fraction equal to
             ! the overall precipitation fraction.
             precip_frac_1(k) = precip_frac(k)
             precip_frac_2(k) = precip_frac(k)


          elseif ( LWP_1(k) >= LWP_tol .and. LWP_2(k) < LWP_tol ) then

             ! LWP_1 is (significantly) greater than 0, while LWP_2 is 0 (or an
             ! insignificant amount).
             !
             ! Precipitation is found at this level, and all cloud water at or
             ! above this level is found in the 1st PDF component.  All of the
             ! hydrometeors are found in the 1st PDF component.
             precip_frac_1(k) = precip_frac(k) / mixt_frac(k)
             precip_frac_2(k) = zero

             ! Special cases for precip_frac_1.
             if ( precip_frac_1(k) > one ) then

                ! Precipitation fraction cannot be greater than 1.
                precip_frac_1(k) = one

                precip_frac_2(k) &
                = ( precip_frac(k) - mixt_frac(k) * precip_frac_1(k) ) &
                  / ( one - mixt_frac(k) )

                ! Check that precip_frac_2 is at least precip_frac_tol.
                if ( precip_frac_2(k) < precip_frac_tol ) then

                   precip_frac_2(k) = precip_frac_tol

                   precip_frac_1(k) &
                   = ( precip_frac(k) &
                       - ( one - mixt_frac(k) ) * precip_frac_2(k) ) &
                     / mixt_frac(k)

                endif

             endif ! Special cases for precip_frac_1


          elseif ( LWP_2(k) >= LWP_tol .and. LWP_1(k) < LWP_tol ) then

             ! LWP_2 is (significantly) greater than 0, while LWP_1 is 0 (or an
             ! insignificant amount).
             !
             ! Precipitation is found at this level, and all cloud water at or
             ! above this level is found in the 2nd PDF component.  All of the
             ! hydrometeors are found in the 2nd PDF component.
             precip_frac_1(k) = zero
             precip_frac_2(k) = precip_frac(k) / ( one - mixt_frac(k) )

             ! Special cases for precip_frac_2.
             if ( precip_frac_2(k) > one ) then

                ! Precipitation fraction cannot be greater than 1.
                precip_frac_2(k) = one

                precip_frac_1(k) &
                = ( precip_frac(k) &
                    - ( one - mixt_frac(k) ) * precip_frac_2(k) ) &
                  / mixt_frac(k)

                ! Check that precip_frac_1 is at least precip_frac_tol.
                if ( precip_frac_1(k) < precip_frac_tol ) then

                   precip_frac_1(k) = precip_frac_tol

                   precip_frac_2(k) &
                   = ( precip_frac(k) - mixt_frac(k) * precip_frac_1(k) ) &
                     / ( one - mixt_frac(k) )

                endif

             endif ! Special cases for precip_frac_2


          else ! LWP_1(k) >= LWP_tol and LWP_2(k) >= LWP_tol

             ! Both LWP_1 and LWP_2 are (significantly) greater than 0.
             !
             ! Precipitation is found at this level, and there is sufficient
             ! cloud water at or above this level in both PDF components to find
             ! the hydrometeor in both PDF components.  Delegate the
             ! precipitation fraction between the 1st and 2nd PDF components.
             precip_frac_1(k) &
             = precip_frac(k) &
               / ( mixt_frac(k) + ( one - mixt_frac(k) ) * LWP_2(k)/LWP_1(k) )

             ! Special cases for precip_frac_1.
             if ( precip_frac_1(k) > one ) then

                ! Precipitation fraction cannot be greater than 1.
                precip_frac_1(k) = one

             elseif ( precip_frac_1(k) < precip_frac_tol ) then

                ! In a scenario where we find precipitation in the 1st PDF
                ! component but it is tiny (less than tolerance level), boost
                ! 1st PDF component precipitation fraction to tolerance level.
                precip_frac_1(k) = precip_frac_tol

             endif ! Special cases for precip_frac_1

             ! Calculate precip_frac_2
             precip_frac_2(k) &
             = ( precip_frac(k) - mixt_frac(k) * precip_frac_1(k) ) &
               / ( one - mixt_frac(k) )

             ! Special cases for precip_frac_2.
             if ( precip_frac_2(k) > one ) then

                ! Precipitation fraction cannot be greater than 1.
                precip_frac_2(k) = one

                ! Recalculate the precipitation fraction in PDF component 1.
                precip_frac_1(k) &
                = ( precip_frac(k) &
                    - ( one - mixt_frac(k) ) * precip_frac_2(k) ) &
                  / mixt_frac(k)

                ! Double check for errors in PDF component 1.
                if ( precip_frac_1(k) > one ) then
                   precip_frac_1(k) = one
                elseif ( precip_frac_1(k) < precip_frac_tol ) then
                   precip_frac_1(k) = precip_frac_tol
                endif

             elseif ( precip_frac_2(k) < precip_frac_tol ) then

                ! In a scenario where we find precipitation in the 2nd PDF
                ! component but it is tiny (less than tolerance level), boost
                ! 2nd PDF component precipitation fraction to tolerance level.
                precip_frac_2(k) = precip_frac_tol

                ! Recalculate the precipitation fraction in PDF component 1.
                precip_frac_1(k) &
                = ( precip_frac(k) &
                    - ( one - mixt_frac(k) ) * precip_frac_2(k) ) &
                  / mixt_frac(k)

                ! Double check for errors in PDF component 1.
                if ( precip_frac_1(k) > one ) then
                   precip_frac_1(k) = one
                elseif ( precip_frac_1(k) < precip_frac_tol ) then
                   precip_frac_1(k) = precip_frac_tol
                endif

             endif ! Special case for precip_frac_2

          endif ! LWP_1(k) < LWP_tol .and. LWP_2(k) < LWP_tol


       else ! all( hydromet(k,:) < hydromet_tol(:) )

          ! All hydrometeors have a value that is either 0 or below tolerance
          ! value (any postive value is considered to be a numerical artifact).
          ! Simply set each PDF component precipitation fraction equal to 0.
          precip_frac_1(k) = zero
          precip_frac_2(k) = zero

       endif ! any( hydromet(k,:) >= hydromet_tol(:) )

    enddo ! k = 1, nz, 1


    ! Statistics
    if ( l_stats_samp ) then

       if ( iLWP1 > 0 ) then
          ! Liquid water path in PDF component 1.
          call stat_update_var( iLWP1, LWP_1, stats_zt )
       endif

       if ( iLWP2 > 0 ) then
          ! Liquid water path in PDF component 2.
          call stat_update_var( iLWP2, LWP_2, stats_zt )
       endif
       
    endif


    return

  end subroutine component_precip_frac_ratio

!===============================================================================

end module precipitation_fraction
