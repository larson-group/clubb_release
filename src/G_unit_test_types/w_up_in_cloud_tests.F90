module w_up_in_cloud_tests

  implicit none

  private

  public :: w_up_in_cloud_tests_driver

  private :: w_up_in_cloud_setup_tests

contains

  !=============================================================================
  function w_up_in_cloud_tests_driver(gr)

    ! Description:
    !   Returns the driver with exit code 0 if tests were successful or exit 
    !   code 1 if any test failed.
    !   Tests check whether computation of w_up_in_cloud within pdf_closure_module
    !   satisfy some given boundary cases where an exact solution is known.
    ! References:
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd
        
    use grid_class, only: &
        grid ! Type

    implicit none
    
    ! Input var
    type (grid), target, intent(in) :: gr

    ! Output var
    integer :: w_up_in_cloud_tests_driver ! Returns the output of the test

    ! Local vars
    integer :: total_mismatches

    print *, "=================================================="
    print *, " "
    print *, "Performing 'mean vertical velocity in clouds' tests"
    print *, " "
    
    call w_up_in_cloud_setup_tests(gr, total_mismatches)
    
    if ( total_mismatches == 0 ) then
      print *, "Success!"
      w_up_in_cloud_tests_driver = 0 ! Exit Code = 0, Success!
    else  ! total_mismatches > 0
      print *, "The computed result was wrong!"
      w_up_in_cloud_tests_driver = 1 ! Exit Code = 1, Fail
    endif

    print *, " "
    print *, "=================================================="
    print *, " "

return


  end function w_up_in_cloud_tests_driver


  !=============================================================================
  subroutine w_up_in_cloud_setup_tests(gr, total_mismatches)
    
    ! Description:
    !   Tests the subroutine calc_w_up_in_cloud. Test cases for both w_1 and w_2 are:
    !   Case 1: Cloud_frac = 0
    !     if Cloud_frac = 0, then we would divide by 0 in the heaviside function.
    !     Obviously this makes little sense and what we would expect as output
    !     is 0 since if there are no clouds there will be no updraft inside clouds.
    !   
    !   Case 2: sigma = 0
    !     if the variation in a component is zero, then the vertical velocity
    !     is constant in that area. Inside the heaviside function, there are two
    !     terms where we divide by the standard deviation. If it is 0, those
    !     will go to infinity (which is not a problem in exact arithmetic since 
    !     erf(infty) = 1 and exp(-infty) = 0). Thus, for sigma = 0 we expect 
    !     simply max(0, w) to be output from the product of w and the heaviside function
    !
    !   Case 3: w strictly non-Positive WITH small sigma
    !     if w is small enough so that considering the standard deviation, it
    !     will never go above 0, then we expect the average updraft to be 0.
    !
    !  Test cases for testing w_1 and w_2 simultaneously are:
    !  Case 1: C_1=C_2=1, a=0.5, sigma_w_1=sigma_w_2=0
    !     Output is the arithmetic average of w_1 and w_2
    ! References:
    !-----------------------------------------------------------------------
    use grid_class, only: &
        grid ! Type
        
    use constants_clubb, only: &
      one,      &  ! Constant(s)
      two,      &
      three,    &
      zero,     &
      one_half, &
      eps

    use clubb_precision, only: &
      core_rknd
      
    use pdf_closure_module, only: &
      calc_w_up_in_cloud

    implicit none
    
    !----------- Input Variables -----------
    type (grid), target, intent(in) :: gr
    
    !----------- Local Variables -----------
    real( kind = core_rknd ), dimension(gr%nzt) :: result
    
    !----------- Output Variables -----------
    integer, intent(out) :: total_mismatches
    
    real(kind = core_rknd), dimension(gr%nzt) :: &
      mixt_frac,             & ! mixture fraction                             [-]
      cloud_frac_1,          & ! cloud fraction (1st PDF component)           [-]
      cloud_frac_2,          & ! cloud fraction (2nd PDF component)           [-]
      w_1,                   & ! upward velocity (1st PDF component)          [m/s]
      w_2,                   & ! upward velocity (2nd PDF component)          [m/s]
      varnce_w_1,            & ! standard deviation of w (1st PDF component)  [m^2/s^2]
      varnce_w_2,            & ! standard deviation of w (2nd PDF component)  [m^2/s^2]
      result_cmp,            & ! correct result for w_up_in_cloud             [m/s]
      w_down_in_cloud,       & ! Output for mean cloudy downdraft velocity    [m/s]
      cloudy_updraft_frac,   & ! Cloudy updraft fraction                      [-]
      cloudy_downdraft_frac    ! Cloudy downdraft fraction                    [-]

      total_mismatches = 0

      ! For w_1:
      ! Case 1 - No clouds. Positive vertical velocity, with variation
      mixt_frac = one         ! observe only comp. 1
      cloud_frac_1 = zero     ! with no clouds
      cloud_frac_2 = zero
      w_1 = one               ! with updraft
      varnce_w_1 = one        ! and with variation in the vertical velocity
      w_2 = zero
      varnce_w_2 = zero
      result_cmp = zero
      call calc_w_up_in_cloud(gr%nzt, 1, & 
                              mixt_frac, &
                              cloud_frac_1, cloud_frac_2, &
                              w_1, w_2, &
                              varnce_w_1, varnce_w_2, &
                              result, w_down_in_cloud, &
                              cloudy_updraft_frac, cloudy_downdraft_frac)
      total_mismatches = total_mismatches &
                          + COUNT(abs(result - result_cmp) >= eps)
      
      ! Case 2 - No variation. Positive vertical velocity.
      mixt_frac = one
      cloud_frac_1 = one    ! cloud_frac is now 1 for comp. 1
      cloud_frac_2 = zero
      varnce_w_1 = zero     ! and there is no variance in comp. 1
      varnce_w_2 = zero
      w_1 = one
      w_2 = zero
      result_cmp = one
      call calc_w_up_in_cloud(gr%nzt, 1, & 
                              mixt_frac, &
                              cloud_frac_1, cloud_frac_2, &
                              w_1, w_2, &
                              varnce_w_1, varnce_w_2, &
                              result, w_down_in_cloud, &
                              cloudy_updraft_frac, cloudy_downdraft_frac)
      total_mismatches = total_mismatches &
                          + COUNT(abs(result - result_cmp) >= eps)
      
      ! Case 3 - Negative vertical velocity with sufficiently small sigma.
      mixt_frac = one
      cloud_frac_1 = one
      cloud_frac_2 = zero
      varnce_w_1 = zero
      varnce_w_2 = zero
      w_1 = -one            ! now we have purely downdraft
      w_2 = -two        ! since, by convention, w_1 > w_2
      result_cmp = zero
      call calc_w_up_in_cloud(gr%nzt, 1, & 
                              mixt_frac, &
                              cloud_frac_1, cloud_frac_2, &
                              w_1, w_2, &
                              varnce_w_1, varnce_w_2, &
                              result, w_down_in_cloud, &
                              cloudy_updraft_frac, cloudy_downdraft_frac)
      total_mismatches = total_mismatches &
                          + COUNT(abs(result - result_cmp) >= eps)
      
      
      ! For w_2:
      ! Case 1 - No clouds. Positive vertical velocity, with variation
      mixt_frac = zero      ! now we only observe comp. 2
      cloud_frac_1 = zero
      cloud_frac_2 = zero   ! with no clouds
      w_1 = two             ! w_1 has to be greater than w_2 by convention
      w_2 = one             ! and updraft
      varnce_w_1 = zero
      varnce_w_2 = one      ! and variation in vertical velocity
      result_cmp = zero
      call calc_w_up_in_cloud(gr%nzt, 1, & 
                              mixt_frac, &
                              cloud_frac_1, cloud_frac_2, &
                              w_1, w_2, &
                              varnce_w_1, varnce_w_2, &
                              result, w_down_in_cloud, &
                              cloudy_updraft_frac, cloudy_downdraft_frac)
      total_mismatches = total_mismatches &
                          + COUNT(abs(result - result_cmp) >= eps)
      
      ! Case 2 - No variation. Positive vertical velocity.
      mixt_frac = zero
      cloud_frac_1 = zero  
      cloud_frac_2 = one    ! now we have full clouds in comp. 2
      varnce_w_1 = zero
      varnce_w_2 = zero     ! and no variation
      w_1 = two
      w_2 = one
      result_cmp = one
      call calc_w_up_in_cloud(gr%nzt, 1, & 
                              mixt_frac, &
                              cloud_frac_1, cloud_frac_2, &
                              w_1, w_2, &
                              varnce_w_1, varnce_w_2, &
                              result, w_down_in_cloud, &
                              cloudy_updraft_frac, cloudy_downdraft_frac)
      total_mismatches = total_mismatches &
                          + COUNT(abs(result - result_cmp) >= eps)
      
      ! Case 3 - Negative vertical velocity with sufficiently small sigma.
      mixt_frac = zero
      cloud_frac_1 = zero  
      cloud_frac_2 = one
      varnce_w_1 = zero
      varnce_w_2 = zero
      w_1 = zero
      w_2 = -one            ! now we have pure downdraft in comp. 2
      result_cmp = zero
      call calc_w_up_in_cloud(gr%nzt, 1, & 
                              mixt_frac, &
                              cloud_frac_1, cloud_frac_2, &
                              w_1, w_2, &
                              varnce_w_1, varnce_w_2, &
                              result, w_down_in_cloud, &
                              cloudy_updraft_frac, cloudy_downdraft_frac)
      total_mismatches = total_mismatches &
                          + COUNT(abs(result - result_cmp) >= eps)
      
      ! Test w_1 and w_2 together
      ! Equal weights, full clouds, pure updraft, no variation. 
      mixt_frac = one_half  ! weigh both components equally
      cloud_frac_1 = one    ! full clouds in both components
      cloud_frac_2 = one
      w_1 = three           ! updraft in both components
      w_2 = one
      varnce_w_1 = zero     ! no variation in both components
      varnce_w_2 = zero
      result_cmp = two      ! so average updraft is just the arithmetic average
      call calc_w_up_in_cloud(gr%nzt, 1, & 
                              mixt_frac, &
                              cloud_frac_1, cloud_frac_2, &
                              w_1, w_2, &
                              varnce_w_1, varnce_w_2, &
                              result, w_down_in_cloud, &
                              cloudy_updraft_frac, cloudy_downdraft_frac)
      total_mismatches = total_mismatches &
                          + COUNT(abs(result - result_cmp) >= eps)
      

  end subroutine w_up_in_cloud_setup_tests


end module w_up_in_cloud_tests
