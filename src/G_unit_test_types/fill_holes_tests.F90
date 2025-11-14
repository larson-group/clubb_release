!=====================================================================================
! 
! Description:  
!
!   This tests the different hole filling methods we have. 
!
!   The easy_fill_test is as it's named, easy to fill, so all points should be
!   above threshold after the fill.
!
!   The below_thresh_test has a threshold larger than the field average, making
!   it impossible to fill all holes, and instead we just additionally ensure 
!   that no above-threshold levels are driven below.
!   
!   In addition, each test checks:
!       - each method should be mass conservative
!       - no points modified outside of specified range
!       - flipping the fields and running in reverse grid mode is BFB
!=====================================================================================
module fill_holes_tests

    use clubb_precision, only: &
        core_rknd

    use constants_clubb, only: &
        one, &
        zero, &
        two, &
        pi, &
        eps
    
    implicit none

    public :: fill_holes_tests_driver

    private :: easy_fill_test, below_thresh_test

    integer, parameter :: &
        nz = 128, &
        ngrdcol = 16

    integer, parameter :: fill_types(6) = [1, 2, 3, 4, 5, 6]

    character(len=20), parameter :: &
        fill_type_names(6) = [  character(len=20) :: &
                                "global", &
                                "sliding_window", &
                                "widening_windows", &
                                "smart_window", &
                                "smart_window_smooth", &
                                "parallel"]

    logical, parameter :: &                                
        l_print_name = .false.

    contains

    !============================= DRIVER =============================
    function fill_holes_tests_driver()

        implicit none

        integer :: fill_holes_tests_driver

        print *, ""
        print *, "====================== Running fill holes tests ====================== "

        fill_holes_tests_driver = 0

        print *, ""
        print *, "-------------- easy_fill_test --------------"
        if ( easy_fill_test() /= 0 ) then
            print *, "FAIL: easy_fill_test"
            fill_holes_tests_driver = fill_holes_tests_driver + 1
        else 
            print *, "PASS: easy_fill_test"
        end if
        print *, "-----------------------------------------------"

        print *, ""
        print *, "-------------- below_thresh_test --------------"
        if ( below_thresh_test() /= 0 ) then
            print *, "FAIL: below_thresh_test"
            fill_holes_tests_driver = fill_holes_tests_driver + 1
        else 
            print *, "PASS: below_thresh_test"
        end if
        print *, "-----------------------------------------------"

        print *, ""
        print *, "====================== End fill holes tests ====================== "
        print *, ""


    end function fill_holes_tests_driver

    !============================= Easy Fill Test =============================
    ! This is a test where filling should be easy, and we expect all methods to
    ! correct all holes. In theory, our smart_window_smooth or parallel
    ! methods may not completely correct holes when possible, which is due 
    ! to their "smoothing" contrainsts that take precedence, but this is easy
    ! enough to fill that any resonable setting of those soomthing parameters
    ! will not cause the fill to fail.
    function easy_fill_test() 

        use fill_holes, only: &
            fill_holes_vertical_api

        implicit none

        integer :: &
            easy_fill_test

        real( kind = core_rknd ) :: &
            threshold, &        ! Field threshold 
            initial_mass, &     ! The initial mass of the field
            new_mass            ! Mass after the fill

        real( kind = core_rknd ), dimension(ngrdcol,nz) :: & 
            rho_ds, &           ! Dry, static density
            dz, &               ! Spacing between grid levels
            field, &            ! Field to fill
            field_initial, &    ! Save of field
            rho_ds_rev, &       ! Reverse version of rho_ds
            dz_rev, &           ! Reverse version of dz_rev
            field_rev, &        ! Reverse version of field
            field_rev_rev       ! Reverse version of field_rev (should be same as field) 

        integer :: &
            lower_hf_level, &   ! Lower hole-filling level
            upper_hf_level, &   ! Upper hole-filling level
            grid_dir_indx       ! Grid direction (+-1)

        integer :: i, k, fill_type

        !----------------------- Begin Code -----------------------

        ! Initialize return code to 0
        easy_fill_test = 0

        ! Loop over each fill type defined by the "fill_types" vector
        do fill_type = 1, size(fill_types)

            print *, "--- Testing fill type: ", fill_type_names(fill_type), "---"

            ! Set defaults, filling to 2:nz-1
            lower_hf_level = 2
            upper_hf_level = nz-1
            grid_dir_indx = 1

            ! This is small for the field we set, making the fill easy
            threshold = 0.05_core_rknd 

            ! Initialize values, we contructed functions for these that have some fun filling 
            ! properties. There's small hole groups, a large hole group, and a big single hole.
            field_initial   = threshold
            field           = threshold
            dz              = 40.0_core_rknd
            rho_ds          = one
            do k = lower_hf_level, upper_hf_level
                do i = 1, ngrdcol
                    rho_ds(i,k) = 5_core_rknd * exp( - k / 120._core_rknd )
                    field_initial(i,k)  =  10_core_rknd * ( sin( 8.0 * pi * ( k + i ) / nz ) + one ) &
                                                        * exp(- k**2 / 2000.0_core_rknd )
                end do
            end do

            ! Let's make a big hole to force the smart window method to
            ! recalculate the range at least once
            field_initial(:,55) = - 50._core_rknd
            field = field_initial

            ! Ensure that we didn't accidentally set all field values above threshold to start with
            if ( all( field_initial > threshold ) ) then
                print *, "ERROR in easy_fill_test: initial field is holeless"
                print *, "nothing is being tested"
                easy_fill_test = easy_fill_test + 1
            end if

            ! Call the fill
            if ( l_print_name ) print *, "filling holes of: easy_fill_test_"//trim(fill_type_names(fill_type))
            call fill_holes_vertical_api( nz, ngrdcol, threshold, &
                                        lower_hf_level, upper_hf_level, &
                                        dz, rho_ds, grid_dir_indx, &
                                        fill_types(fill_type), &
                                        field )

            ! This case is easily fillable, so no fields should be under threshold
            if ( any( field < threshold ) ) then
                print *, "ERROR in easy_fill_test: there were below threshold values"
                print *, "in the field after the fill. This shouldn't happen, there is"
                print *, "plenty of excess mass to fill the holes with."
                do k = 1, nz
                    do i = 1, ngrdcol
                    if ( field(i,k) < threshold ) then
                        write(*,'(A6,I5,1X,I5,A4,E30.20,A3,E30.20)') "field(", i, k, ") = ", field(i,k), " < ", threshold
                    end if
                    end do
                end do
                easy_fill_test = easy_fill_test + 1
            end if

            ! Calculate the mass, this needs to be conservative, a signficant difference is bad
            initial_mass = sum( field_initial * rho_ds * dz ) 
            new_mass     = sum( field         * rho_ds * dz ) 

            if( two * abs( initial_mass - new_mass ) / ( initial_mass + new_mass ) > eps ) then
                print *, "ERROR in easy_fill_test: method was not conservative"
                print *, "-- initial mass vs after fill_holes: ", initial_mass, new_mass
                easy_fill_test = easy_fill_test + 1
            end if

            ! Ensure that the only the field values in the range [lower_hf_level:upper_hf_level]
            ! were modified
            if ( any( field(:,1:lower_hf_level) /= field(:,1:lower_hf_level) ) &
                 .or. any( field(:,upper_hf_level:nz) /= field(:,upper_hf_level:nz) ) ) then
                print *, "ERROR in easy_fill_test: vertical levels outside the range specified"
                print *, "by [lower_hf_level:upper_hf_level] were modified"
                easy_fill_test = easy_fill_test + 1
            end if


            !---------------------- Reverse grid test ----------------------
            ! Perform the same fill as above, but on flipped arrays in in
            ! reverse grid mode. The result should be BFB, but flipped.

            ! Set new bounds and grid direction to revserse mode
            lower_hf_level = nz-1
            upper_hf_level = 2
            grid_dir_indx = -1

            ! Flip arrays
            do k = 1, nz
                do i = 1, ngrdcol
                    dz_rev(i,nz+1-k)     = dz(i,k)
                    rho_ds_rev(i,nz+1-k) = rho_ds(i,k)
                    field_rev(i,nz+1-k)  = field_initial(i,k)
                end do
            end do

            ! Perform reverse fill
            if ( l_print_name ) print *, "filling holes of: reverse_easy_fill_test_"//trim(fill_type_names(fill_type))
            call fill_holes_vertical_api( nz, ngrdcol, threshold, &
                                          lower_hf_level, upper_hf_level, &
                                          dz_rev, rho_ds_rev, grid_dir_indx, &
                                          fill_types(fill_type), &
                                          field_rev )

            ! This new mass should also be close to original
            new_mass = sum( field_rev * rho_ds_rev * dz_rev ) 

            if( two * abs( initial_mass - new_mass ) / ( initial_mass + new_mass ) > eps ) then
                print *, "ERROR in easy_fill_test REVERSE MODE: method was not conservative"
                print *, "-- initial mass vs after fill_holes: ", initial_mass, new_mass
                easy_fill_test = easy_fill_test + 1
            end if

            ! Flip the flipped array to make for easy comparison with the normal fill
            do k = 1, nz
                do i = 1, ngrdcol
                    field_rev_rev(i,k) = field_rev(i,nz+1-k)
                end do
            end do

            ! The reversed the solution from the reverse fill should be the same as the 
            ! result from the forward fill (field)
            if ( any( field_rev_rev /= field ) ) then
                print *, "ERROR in easy_fill_test REVERSE MODE: grid flipped fill"
                print *, "did not produce the same result"
                easy_fill_test = easy_fill_test + 1
            end if

        end do

        return

    end function easy_fill_test


    !============================= Below Thresh Test =============================
    ! This is a test where filling should be impossible, since the threshold
    ! is set to a value greater than the field average. This means we drop the 
    ! criteria that all holes be filled, and replace it with a check to ensure that
    ! nothing that was above threshold was driven below.
    function below_thresh_test() 

        use fill_holes, only: &
            fill_holes_vertical_api

        implicit none

        integer :: &
            below_thresh_test

        real( kind = core_rknd ) :: &
            threshold, &        ! Field threshold 
            initial_mass, &     ! The initial mass of the field
            new_mass            ! Mass after the fill

        real( kind = core_rknd ), dimension(ngrdcol,nz) :: & 
            rho_ds, &           ! Dry, static density
            dz, &               ! Spacing between grid levels
            field, &            ! Field to fill
            field_initial, &    ! Save of field
            rho_ds_rev, &       ! Reverse version of rho_ds
            dz_rev, &           ! Reverse version of dz_rev
            field_rev, &        ! Reverse version of field
            field_rev_rev       ! Reverse version of field_rev (should be same as field) 

        integer :: &
            lower_hf_level, &   ! Lower hole-filling level
            upper_hf_level, &   ! Upper hole-filling level
            grid_dir_indx       ! Grid direction (+-1)

        integer :: i, k, fill_type

        !----------------------- Begin Code -----------------------

        ! Initialize return code to 0
        below_thresh_test = 0

        ! Loop over each fill type defined by the "fill_types" vector
        do fill_type = 1, size(fill_types)

            print *, "--- Testing fill type: ", fill_type_names(fill_type), "---"

            ! Set defaults, filling to 2:nz-1
            lower_hf_level = 2
            upper_hf_level = nz-1
            grid_dir_indx = 1

            ! This is a large threshold for the field we set. The field average is below
            ! this, so completely filling all holes is impossible
            threshold = 5._core_rknd  

            ! Initialize values, we contructed functions for these that have some fun filling 
            ! properties. There's small hole groups, a large hole group, and a big single hole.
            field_initial   = threshold
            field           = threshold
            dz              = 40.0_core_rknd
            rho_ds          = one
            do k = lower_hf_level, upper_hf_level
                do i = 1, ngrdcol
                    rho_ds(i,k) = 5_core_rknd * exp( - k / 120._core_rknd )
                    field_initial(i,k)  =  10_core_rknd * ( sin( 8.0 * pi * ( k + i ) / nz ) + one ) &
                                                        * exp(- k**2 / 2000.0_core_rknd )
                end do
            end do

            ! Let's make a big hole to force the smart window method to
            ! recalculate the range at least once
            field_initial(:,55) = - 50._core_rknd
            field = field_initial

            ! Ensure that the field average is below the threshold, otherwise this test is
            ! improperly configured
            if ( sum( field_initial * rho_ds * dz ) / sum( rho_ds * dz ) > threshold ) then
                print *, "ERROR in below_thresh_test: initial field average is above thresold"
                print *, "this test is not configured correctly."
                below_thresh_test = below_thresh_test + 1
            end if

            ! Call the fill
            if ( l_print_name ) print *, "filling holes of: below_thresh_test_"//trim(fill_type_names(fill_type))
            call fill_holes_vertical_api( nz, ngrdcol, threshold, &
                                        lower_hf_level, upper_hf_level, &
                                        dz, rho_ds, grid_dir_indx, &
                                        fill_types(fill_type), &
                                        field )

            ! Ensure that no above threshold values were driven below threshold
            if ( any( field_initial >= threshold .and. field < threshold ) ) then
                print *, "ERROR in below_thresh_test: an above threshold value has been"
                print *, "driven below threshold. This is highly frowned upon."
                below_thresh_test = below_thresh_test + 1
            end if


            ! Calculate the mass, this needs to be conservative, a signficant difference is bad
            initial_mass = sum( field_initial * rho_ds * dz ) 
            new_mass     = sum( field         * rho_ds * dz ) 

            if( two * abs( initial_mass - new_mass ) / ( initial_mass + new_mass ) > eps ) then
                print *, "ERROR in below_thresh_test: method was not conservative"
                print *, "-- initial mass vs after fill_holes: ", initial_mass, new_mass
                print *, "-- error", two * abs( initial_mass - new_mass ) / ( initial_mass + new_mass )
                below_thresh_test = below_thresh_test + 1
            end if

            
            ! Ensure that the only the field values in the range [lower_hf_level:upper_hf_level]
            ! were modified
            if ( any( field(:,1:lower_hf_level) /= field(:,1:lower_hf_level) ) &
                 .or. any( field(:,upper_hf_level:nz) /= field(:,upper_hf_level:nz) ) ) then
                print *, "ERROR in below_thresh_test: vertical levels outside the range specified"
                print *, "by [lower_hf_level:upper_hf_level] were modified"
                below_thresh_test = below_thresh_test + 1
            end if


            !---------------------- Reverse grid test ----------------------
            ! Perform the same fill as above, but on flipped arrays in in
            ! reverse grid mode. The result should be BFB, but flipped.

            ! Set new bounds and grid direction to revserse mode
            lower_hf_level = nz-1
            upper_hf_level = 2
            grid_dir_indx = -1

            ! Flip arrays
            do k = 1, nz
                do i = 1, ngrdcol
                    dz_rev               = dz(i,k)
                    rho_ds_rev(i,nz+1-k) = rho_ds(i,k)
                    field_rev(i,nz+1-k)  = field_initial(i,k)
                end do
            end do

            ! Perform reverse fill
            if ( l_print_name ) print *, "filling holes of: reverse_below_thresh_test_"//trim(fill_type_names(fill_type))
            call fill_holes_vertical_api( nz, ngrdcol, threshold, &
                                          lower_hf_level, upper_hf_level, &
                                          dz_rev, rho_ds_rev, grid_dir_indx, &
                                          fill_types(fill_type), &
                                          field_rev )

            ! This new mass should also be close to original
            new_mass = sum( field_rev * rho_ds_rev * dz_rev ) 

            if( two * abs( initial_mass - new_mass ) / ( initial_mass + new_mass ) > eps ) then
                print *, "ERROR in below_thresh_test REVERSE MODE: method was not conservative in grid reverse mode"
                print *, "-- initial mass vs after fill_holes: ", initial_mass, new_mass
                below_thresh_test = below_thresh_test + 1
            end if

            ! Flip the flipped array to make for easy comparison with the normal fill
            do k = 1, nz
                do i = 1, ngrdcol
                    field_rev_rev(i,k) = field_rev(i,nz+1-k)
                end do
            end do

            ! The reversed the solution from the reverse fill should be the same as the 
            ! result from the forward fill (field)
            if ( any( field_rev_rev /= field ) ) then
                print *, "ERROR in below_thresh_test REVERSE MODE: grid flipped fill"
                print *, "did not produce the same result"
                below_thresh_test = below_thresh_test + 1
            end if

        end do

        return

    end function below_thresh_test

end module fill_holes_tests
