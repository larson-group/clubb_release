!---------------------------------------------------------------------------
! $Id$
!===============================================================================
module Nc_Ncn_test

  implicit none

  private ! default scope

  public :: Nc_Ncn_unit_test

  contains

  !=============================================================================
  function Nc_Ncn_unit_test()

    ! Description:
    ! Unit testing for the Nc-Ncn "back-and-forth" code.
    
    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one_hundred, & ! Variable(s)
        sqrt_2,      &
        one,         &
        one_half,    &
        zero,        &
        fstdout

    use Nc_Ncn_eqns, only: &
        Ncnm_to_Nc_in_cloud, & ! Procedure(s)
        Nc_in_cloud_to_Ncnm, &
        Ncnm_to_Ncm, &
        Ncm_to_Ncnm

    use pdf_utilities, only: &
        stdev_L2N,  & ! Procedure(s)
        corr_NL2NN

    use clubb_precision, only: &
        core_rknd  ! Variable(s)

    implicit none

    ! Output Variables
    integer :: Nc_Ncn_unit_test ! Returns the exit code of the test

    ! Local Variables
    integer, parameter :: &
      nz = 1

    real( kind = core_rknd ), dimension(nz) :: &
      mu_chi_1,         & ! Mean of chi(s) (1st PDF component)           [kg/kg]
      mu_chi_2,         & ! Mean of chi(s) (2nd PDF component)           [kg/kg]
      mu_Ncn_1,         & ! Mean of Ncn (1st PDF component)             [num/kg]
      mu_Ncn_2,         & ! Mean of Ncn (2nd PDF component)             [num/kg]
      sigma_chi_1,      & ! Standard deviation of chi(s) (1st PDF comp.) [kg/kg]
      sigma_chi_2,      & ! Standard deviation of chi(s) (2nd PDF comp.) [kg/kg]
      sigma_Ncn_1,      & ! Standard deviation of Ncn (1st PDF comp.)   [num/kg]
      sigma_Ncn_2,      & ! Standard deviation of Ncn (2nd PDF comp.)   [num/kg]
      sigma_Ncn_1_n,    & ! Standard deviation of ln Ncn (1st PDF component) [-]
      sigma_Ncn_2_n,    & ! Standard deviation of ln Ncn (2nd PDF component) [-]
      corr_chi_Ncn_1_n, & ! Correlation of chi(s) and ln Ncn (1st PDF comp.) [-]
      corr_chi_Ncn_2_n, & ! Correlation of chi(s) and ln Ncn (2nd PDF comp.) [-]
      mixt_frac,        & ! Mixture fraction                                 [-]
      cloud_frac_1,     & ! Cloud fraction (1st PDF component)               [-]
      cloud_frac_2        ! Cloud fraction (2nd PDF component)               [-]

    real( kind = core_rknd ), dimension(nz) :: &
      Ncm,         & ! Mean cloud droplet concentration (overall)       [num/kg]
      Nc_in_cloud, & ! Mean cloud droplet concentration (in-cloud)      [num/kg]
      Ncnm           ! Mean simplified cloud nuclei concentration       [num/kg]

    real( kind = core_rknd ) :: &
      const_Ncnp2_on_Ncnm2, & ! Prescribed ratio of <Ncn'^2> to <Ncn>        [-]
      const_corr_chi_Ncn      ! Prescribed correlation of chi(s) and Ncn     [-]

    real( kind = core_rknd ), dimension(nz) :: &
      const_Ncnp2_on_Ncnm2_in    ! Prescribed ratio of <Ncn'^2> to <Ncn>     [-]

    real( kind = core_rknd ), dimension(nz) :: &
      percent_diff    ! Percent difference between Ncnm before and after     [%]

    real( kind = core_rknd ), parameter :: &
      tol = 1.0e-13_core_rknd    ! Tolerance for percent difference          [%]

    integer :: &
      num_mismatches  ! Total number of mismatches

    integer :: iter  ! Loop index

  !-----------------------------------------------------------------------
    !----- Begin Code -----

    write(fstdout,*) ""
    write(fstdout,*) "Nc-Ncn equations unit test"
    write(fstdout,*) "=========================="
    write(fstdout,*) ""

    ! Initialize number of mismatches.
    num_mismatches = 0

    do iter = 1, 10, 1

       write(fstdout,*) "Iteration = ", iter
       write(fstdout,*) ""

       if ( iter == 1 ) then

          mu_chi_1 = -1.0e-5_core_rknd
          mu_chi_2 = -5.0e-5_core_rknd

          sigma_chi_1 = 1.0e-5_core_rknd
          sigma_chi_2 = 1.5e-5_core_rknd

          mu_Ncn_1 = 5.0e+7_core_rknd

          sigma_Ncn_1 = 2.5e+7_core_rknd

          const_corr_chi_Ncn = 0.4_core_rknd

          mixt_frac = 0.3_core_rknd

       elseif ( iter == 2 ) then

          mu_chi_1 = -1.0e-5_core_rknd
          mu_chi_2 = -5.0e-5_core_rknd

          sigma_chi_1 = 1.0e-9_core_rknd
          sigma_chi_2 = 1.5e-9_core_rknd

          mu_Ncn_1 = 5.0e+7_core_rknd

          sigma_Ncn_1 = 2.5e+7_core_rknd

          const_corr_chi_Ncn = 0.4_core_rknd

          mixt_frac = 0.3_core_rknd

       elseif ( iter == 3 ) then

          mu_chi_1 = 1.0e-5_core_rknd
          mu_chi_2 = 5.0e-5_core_rknd

          sigma_chi_1 = 1.0e-5_core_rknd
          sigma_chi_2 = 1.5e-5_core_rknd

          mu_Ncn_1 = 5.0e+7_core_rknd

          sigma_Ncn_1 = 7.5e+7_core_rknd

          const_corr_chi_Ncn = 0.4_core_rknd

          mixt_frac = 0.3_core_rknd

       elseif ( iter == 4 ) then

          mu_chi_1 = 1.0e-8_core_rknd
          mu_chi_2 = -1.0e-8_core_rknd

          sigma_chi_1 = 5.0e-5_core_rknd
          sigma_chi_2 = 5.0e-5_core_rknd

          mu_Ncn_1 = 5.0e+7_core_rknd

          sigma_Ncn_1 = 5.0e+7_core_rknd

          const_corr_chi_Ncn = 0.1_core_rknd

          mixt_frac = 0.5_core_rknd

       elseif ( iter == 5 ) then

          mu_chi_1 = 1.0e-5_core_rknd
          mu_chi_2 = 5.0e-5_core_rknd

          sigma_chi_1 = 1.0e-6_core_rknd
          sigma_chi_2 = 1.5e-6_core_rknd

          mu_Ncn_1 = 5.0e+7_core_rknd

          sigma_Ncn_1 = 2.5e+7_core_rknd

          const_corr_chi_Ncn = 0.5_core_rknd

          mixt_frac = 0.3_core_rknd

       elseif ( iter == 6 ) then

          mu_chi_1 = -1.0e-7_core_rknd
          mu_chi_2 = -5.0e-5_core_rknd

          sigma_chi_1 = 1.0e-5_core_rknd
          sigma_chi_2 = 1.5e-5_core_rknd

          mu_Ncn_1 = 5.0e+7_core_rknd

          sigma_Ncn_1 = 5.0e+7_core_rknd

          const_corr_chi_Ncn = zero

          mixt_frac = 0.2_core_rknd

       elseif ( iter == 7 ) then

          mu_chi_1 = 1.0e-6_core_rknd
          mu_chi_2 = -5.0e-7_core_rknd

          sigma_chi_1 = 1.0e-5_core_rknd
          sigma_chi_2 = 1.5e-5_core_rknd

          mu_Ncn_1 = 5.0e+7_core_rknd

          sigma_Ncn_1 = zero

          const_corr_chi_Ncn = 0.7_core_rknd

          mixt_frac = 0.2_core_rknd

       elseif ( iter == 8 ) then

          mu_chi_1 = 1.0e-6_core_rknd
          mu_chi_2 = -5.0e-7_core_rknd

          sigma_chi_1 = 1.0e-5_core_rknd
          sigma_chi_2 = 1.5e-9_core_rknd

          mu_Ncn_1 = 5.0e+7_core_rknd

          sigma_Ncn_1 = 7.5e+7_core_rknd

          const_corr_chi_Ncn = 0.7_core_rknd

          mixt_frac = 0.2_core_rknd

       elseif ( iter == 9 ) then

          mu_chi_1 = 1.0e-5_core_rknd
          mu_chi_2 = -5.0e-5_core_rknd

          sigma_chi_1 = 1.0e-9_core_rknd
          sigma_chi_2 = 1.5e-5_core_rknd

          mu_Ncn_1 = 5.0e+7_core_rknd

          sigma_Ncn_1 = 7.5e+7_core_rknd

          const_corr_chi_Ncn = 0.7_core_rknd

          mixt_frac = 0.2_core_rknd

       elseif ( iter == 10 ) then

          mu_chi_1 = 1.0e-5_core_rknd
          mu_chi_2 = -5.0e-5_core_rknd

          sigma_chi_1 = 1.0e-5_core_rknd
          sigma_chi_2 = 1.5e-5_core_rknd

          mu_Ncn_1 = 5.0e+7_core_rknd

          sigma_Ncn_1 = 7.5e+7_core_rknd

          const_corr_chi_Ncn = 0.7_core_rknd

          mixt_frac = 0.2_core_rknd

       endif ! iter

       mu_Ncn_2 = mu_Ncn_1

       sigma_Ncn_2 = sigma_Ncn_1

       write(fstdout,*) "mu_chi_1 = ", mu_chi_1
       write(fstdout,*) "mu_chi_2 = ", mu_chi_2
       write(fstdout,*) "sigma_chi_1 = ", sigma_chi_1
       write(fstdout,*) "sigma_chi_2 = ", sigma_chi_2
       write(fstdout,*) "mu_Ncn_1 = mu_Ncn_2 = ", mu_Ncn_1
       write(fstdout,*) "sigma_Ncn_1 = sigma_Ncn_2 = ", sigma_Ncn_1
       write(fstdout,*) "corr_chi_Ncn_1 = corr_chi_Ncn_2 = ", const_corr_chi_Ncn
       write(fstdout,*) "mixture fraction = ", mixt_frac

       const_Ncnp2_on_Ncnm2_in = sigma_Ncn_1**2 / mu_Ncn_1**2

       const_Ncnp2_on_Ncnm2 = const_Ncnp2_on_Ncnm2_in(nz)

       sigma_Ncn_1_n = stdev_L2N( const_Ncnp2_on_Ncnm2_in ) 

       sigma_Ncn_2_n = sigma_Ncn_1_n

       if ( const_Ncnp2_on_Ncnm2 > zero ) then

          corr_chi_Ncn_1_n &
          = corr_NL2NN( const_corr_chi_Ncn, sigma_Ncn_1_n(nz), &
                        const_Ncnp2_on_Ncnm2 ) 

          corr_chi_Ncn_2_n = corr_chi_Ncn_1_n

       endif

       cloud_frac_1 &
       = one_half * erfc( - ( mu_chi_1 / ( sqrt_2 * sigma_chi_1 ) ) )

       cloud_frac_2 &
       = one_half * erfc( - ( mu_chi_2 / ( sqrt_2 * sigma_chi_2 ) ) )

       write(fstdout,*) "cloud fraction (1st PDF component) = ", cloud_frac_1

       write(fstdout,*) "cloud fraction (2nd PDF component) = ", cloud_frac_2

       write(fstdout,*) "cloud fraction = ", &
                        mixt_frac * cloud_frac_1 &
                        + ( one - mixt_frac ) * cloud_frac_2


       ! Calculate Ncm and Nc_in_cloud from PDF parameters (Ncn) and cloud
       ! fraction.
       Nc_in_cloud &
       = Ncnm_to_Nc_in_cloud( nz, mu_chi_1, mu_chi_2, mu_Ncn_1, mu_Ncn_2, &
                              sigma_chi_1, sigma_chi_2, sigma_Ncn_1, &
                              sigma_Ncn_2, sigma_Ncn_1_n, sigma_Ncn_2_n, &
                              corr_chi_Ncn_1_n, corr_chi_Ncn_2_n, mixt_frac, &
                              cloud_frac_1, cloud_frac_2 )

       write(fstdout,*) "Nc_in_cloud (from Ncn) = ", Nc_in_cloud

       Ncm = Ncnm_to_Ncm( nz, mu_chi_1, mu_chi_2, mu_Ncn_1, mu_Ncn_2, &
                          sigma_chi_1, sigma_chi_2, sigma_Ncn_1, &
                          sigma_Ncn_2, sigma_Ncn_1_n, sigma_Ncn_2_n, &
                          corr_chi_Ncn_1_n, corr_chi_Ncn_2_n, mixt_frac )

       write(fstdout,*) "Ncm (from Ncn) = ", Ncm

       ! Calculate Ncnm from Ncm, Nc_in_cloud, PDF parameters, and cloud
       ! fraction.
       Ncnm = Nc_in_cloud_to_Ncnm( nz, mu_chi_1, mu_chi_2, sigma_chi_1, &
                                   sigma_chi_2, mixt_frac, Nc_in_cloud, &
                                   cloud_frac_1, cloud_frac_2, &
                                   const_Ncnp2_on_Ncnm2, const_corr_chi_Ncn )

       write(fstdout,*) "Ncnm (from Nc_in_cloud) = ", Ncnm

       Ncnm = Ncm_to_Ncnm( nz, mu_chi_1, mu_chi_2, sigma_chi_1, sigma_chi_2, &
                           mixt_frac, Ncm, const_Ncnp2_on_Ncnm2, &
                           const_corr_chi_Ncn, Nc_in_cloud )

       write(fstdout,*) "Ncnm (from Ncm) = ", Ncnm

       ! Calculate the percent difference between Ncnm (after the back-and-forth
       ! calculations) and mu_Ncn_1 (which is Ncnm before the back-and-forth
       ! calculations).
       percent_diff = one_hundred * abs( ( Ncnm - mu_Ncn_1 ) / mu_Ncn_1 )

       if ( percent_diff(nz) <= tol ) then

          ! The percentage difference between Ncnm (here, after back-and-forth
          ! calculations) and mu_Ncn_1 (here, Ncnm before back-and-forth
          ! calculations) is very small -- a product of numerical round-off
          ! error.
          write(fstdout,*) ""
          write(fstdout,'(1x,A,I2,A)') "Test iteration ", iter, " is a success!"

       else ! percent_diff > tol

          ! The percentage difference between Ncnm (here, after back-and-forth
          ! calculations) and mu_Ncn_1 (here, Ncnm before back-and-forth
          ! calculations) is beyond an acceptable tolerance -- there may be a
          ! larger problem.
          write(fstdout,*) ""
          write(fstdout,'(1x,A,I2,A)') "Test iteration ", iter, &
                                       " is not successful.  Please check for" &
                                       // " any changes made to the relevant" &
                                       // " portion(s) of the CLUBB model code."

          num_mismatches = num_mismatches + 1

       endif ! percent_diff <= tol

       write(fstdout,*) ""
       write(fstdout,*) "======================================================"
       write(fstdout,*) ""

    enddo ! iter = 1, 10, 1


    if ( num_mismatches == 0 ) then
       write(fstdout,'(1x,A)') "Success!"
       Nc_Ncn_unit_test = 0 ! Exit Code = 0, Success!
    else ! num_mismatches > 0
       write(fstdout,'(1x,A,I2,A)') "There were ", num_mismatches, &
                                    " total mismatches."
       Nc_Ncn_unit_test = 1 ! Exit Code = 1, Fail
    endif ! num_mismatches = 0


    return

  end function Nc_Ncn_unit_test

!===============================================================================

end module Nc_Ncn_test
