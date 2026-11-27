!-----------------------------------------------------------------------
! $Id$
!-----------------------------------------------------------------------
module p3_silhs_rime_test

  ! Description:
  !   Tests how SILHS prepares P3 rime samples and how CLUBB clips mean rime.
  !
  !   construct_p3_rime_samples builds sample rime mass and volume from sampled
  !   total ice at the grid-mean rime fraction and density.
  !   copy_X_nl_into_hydromet_all_pts copies the prepared samples into the
  !   hydrometeor array. clip_hydromet_conc_mvr bounds mean rime mass by total
  !   ice and scales rime volume with it.
  !
  !   The invalid-input assertions in construct_p3_rime_samples (mean rime mass
  !   above total ice, negative sampled ice) end the program with error stop,
  !   so they cannot be tested here.

  implicit none

  private

  public :: p3_silhs_rime_unit_test

  contains

  !=============================================================================
  function p3_silhs_rime_unit_test()

    ! Description:
    !   Run the rime preparation, copying and clipping routines on synthetic
    !   ice fields for two columns, five levels and three samples. The levels
    !   cover rimed, unrimed, ice-free, volume-free and fully rimed ice.
    !   Checks that:
    !     - non-rime samples pass through unchanged;
    !     - sample rime stays within total ice and keeps the grid-mean rime
    !       fraction and density, and equal-weight samples recover the mean;
    !     - zero fallback means (as used by statistics) keep prepared rime;
    !     - with total ice unsampled, rime follows the mean total ice;
    !     - clipping bounds mean rime by total ice, removes rime without ice or
    !       volume, preserves rime density, total ice and rain, and does not
    !       depend on the order of the rime fields;
    !     - without P3 rime indices, samples keep the previous mean fallback.
    !   Returns 0 if all checks pass and 1 otherwise.

    ! References:
    !   None
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use corr_varnce_module, only: &
        hm_metadata_type ! Type

    use latin_hypercube_driver_module, only: &
        construct_p3_rime_samples,       & ! Procedure(s)
        copy_X_nl_into_hydromet_all_pts

    use fill_holes, only: &
        clip_hydromet_conc_mvr ! Procedure(s)

    use constants_clubb, only: &
        fstderr ! Constant(s)

    implicit none

    ! Local Constants
    integer, parameter :: &
      nc = 2, & ! Number of columns [-]
      nz = 5, & ! Number of levels [-]
      ns = 3, & ! Number of samples [-]
      nh = 6, & ! Number of hydrometeors [-]
      np = 10   ! Number of PDF variables [-]

    real( kind = core_rknd ), parameter :: &
      tol = 1.e-12_core_rknd ! Relative tolerance for ratio checks [-]

    ! Output Variable
    integer :: p3_silhs_rime_unit_test ! Returns the exit code of the test

    ! Local Variables
    type(hm_metadata_type) :: &
      meta ! Hydrometeor and PDF index mapping

    real( kind = core_rknd ), dimension(nc,nz,nh) :: &
      means,          & ! Hydrometeor grid means               [units vary]
      clipped,        & ! Clipped grid means                   [units vary]
      swapped_means,  & ! Means with rime fields swapped       [units vary]
      swapped_clipped   ! Clipped means with rime fields swapped [units vary]

    real( kind = core_rknd ), dimension(nc,ns,nz,np) :: &
      points,         & ! PDF samples                          [units vary]
      original_points   ! PDF samples before rime preparation  [units vary]

    real( kind = core_rknd ), dimension(nc,ns,nz,nh) :: &
      samples ! Hydrometeor samples [units vary]

    real( kind = core_rknd ), dimension(nc,ns,nz) :: &
      nuclei ! Cloud nuclei samples [num/kg]

    real( kind = core_rknd ), dimension(nc) :: &
      fraction, & ! Grid-mean rime mass fraction of total ice [-]
      density     ! Grid-mean rime density                    [kg/m^3]

    real( kind = core_rknd ), dimension(ns) :: &
      factors ! Sample total ice as a multiple of the mean [-]

    real( kind = core_rknd ) :: &
      qi,      & ! Sample total ice mass mixing ratio [kg/kg]
      qr,      & ! Sample rime mass mixing ratio      [kg/kg]
      br,      & ! Sample rime volume                 [m^3/kg]
      expected   ! Expected rime mass fraction        [-]

    integer :: i, k, s ! Column, level and sample indices

    !----- Begin Code -----

    p3_silhs_rime_unit_test = 0

    print *, "=================================================="
    print *, ""
    print *, "Performing P3 SILHS rime test"
    print *, ""
    print *, "=================================================="

    meta%iirr = 1; meta%iiNr = 2; meta%iiri = 3; meta%iiNi = 4
    meta%iiqirim = 5; meta%iibirim = 6
    meta%iirs = -1; meta%iirg = -1; meta%iiNs = -1; meta%iiNg = -1
    meta%iiPDF_rr = 5; meta%iiPDF_Nr = 6
    meta%iiPDF_ri = 7; meta%iiPDF_Ni = 8; meta%iiPDF_Ncn = 4
    meta%iiPDF_qirim = 9; meta%iiPDF_birim = 10
    fraction = [0.25_core_rknd, 0.6_core_rknd]
    density = [400._core_rknd, 700._core_rknd]
    factors = [0._core_rknd, 1._core_rknd, 2._core_rknd]

    means = 0._core_rknd
    points = 0._core_rknd
    do k = 1, nz
      do i = 1, nc
        means(i,k,3) = 1.e-4_core_rknd * i * k
        means(i,k,5) = fraction(i) * means(i,k,3)
        means(i,k,6) = means(i,k,5) / density(i)
      enddo
    enddo

    means(:,2,5:6) = 0._core_rknd ! Unrimed ice
    means(:,3,3:6) = 0._core_rknd ! No ice
    means(:,4,6) = 0._core_rknd   ! Rime mass without a defined density
    means(:,5,5) = means(:,5,3)   ! Fully rimed ice
    means(:,5,6) = means(:,5,5) / density

    do s = 1, ns
      points(:,s,:,4) = 123._core_rknd
      points(:,s,:,5) = 2.e-5_core_rknd * s
      points(:,s,:,6) = 50._core_rknd * s
      points(:,s,:,7) = factors(s) * means(:,:,3)
      points(:,s,:,8) = 300._core_rknd * factors(s)**2
      points(:,s,:,9:10) = 999._core_rknd ! Unrelated raw rime samples must not leak through
    enddo

    original_points = points
    call construct_p3_rime_samples( nz, nc, ns, np, nh, meta, means, & ! In
                                    points )                           ! InOut

    call check( any( points(:,:,:,1:8) /= original_points(:,:,:,1:8) ), &
                'Non-rime PDF samples changed' )

    call copy_X_nl_into_hydromet_all_pts( nz, np, ns, nc, points, nh, meta, means, & ! In
                                          samples, nuclei )                          ! Out

    call check( any( samples(:,:,:,1) /= points(:,:,:,5) ), 'Rain mass changed' )
    call check( any( samples(:,:,:,2) /= points(:,:,:,6) ), 'Rain number changed' )
    call check( any( samples(:,:,:,3) /= points(:,:,:,7) ), 'Total ice changed' )
    call check( any( samples(:,:,:,4) /= points(:,:,:,8) ), 'Ice number changed' )
    call check( any( nuclei /= points(:,:,:,4) ), 'Nuclei samples changed' )
    call check( any( samples(:,:,:,5) /= points(:,:,:,9) ), &
                'Prepared rime mass was not copied' )
    call check( any( samples(:,:,:,6) /= points(:,:,:,10) ), &
                'Prepared rime volume was not copied' )

    do s = 1, ns
      do k = 1, nz
        do i = 1, nc
          qi = samples(i,s,k,3)
          qr = samples(i,s,k,5)
          br = samples(i,s,k,6)
          call check( qr < 0._core_rknd .or. qr > qi, 'Rime outside total ice' )
          if ( qi == 0._core_rknd .or. ( k >= 2 .and. k <= 4 ) ) then
            call check( qr /= 0._core_rknd .or. br /= 0._core_rknd, 'Zero-rime case failed' )
          else
            expected = fraction(i)
            if ( k == 5 ) expected = 1._core_rknd
            call check( abs( qr / qi - expected ) > tol, 'Rime fraction changed' )
            call check( abs( qr / br / density(i) - 1._core_rknd ) > tol, &
                        'Rime density changed' )
          endif
        enddo
      enddo
    enddo

    ! These equal-weight ice samples exactly represent the input mean.
    do i = 1, nc
      call check( abs( sum( samples(i,:,1,5) ) / ns / means(i,1,5) - 1._core_rknd ) > tol, &
                  'Mean rime mass not recovered' )
      call check( abs( sum( samples(i,:,1,6) ) / ns / means(i,1,6) - 1._core_rknd ) > tol, &
                  'Mean rime volume not recovered' )
    enddo

    ! Statistics pass zero means for unsampled fields; prepared rime must survive.
    swapped_means = 0._core_rknd
    call copy_X_nl_into_hydromet_all_pts( nz, np, ns, nc, points, nh, meta, swapped_means, & ! In
                                          samples, nuclei )                                  ! Out
    call check( any( samples(:,:,:,5) /= points(:,:,:,9) ), &
                'Zero fallback means changed prepared rime mass' )
    call check( any( samples(:,:,:,6) /= points(:,:,:,10) ), &
                'Zero fallback means changed prepared rime volume' )

    ! With total ice sampling disabled, construct rime from mean total ice.
    meta%iiPDF_ri = -1
    call construct_p3_rime_samples( nz, nc, ns, np, nh, meta, means, & ! In
                                    points )                           ! InOut
    call copy_X_nl_into_hydromet_all_pts( nz, np, ns, nc, points, nh, meta, means, & ! In
                                          samples, nuclei )                          ! Out
    do s = 1, ns
      call check( any( samples(:,s,:,3) /= means(:,:,3) ), 'Mean ice fallback changed' )
      call check( any( abs( samples(:,s,1,5) / means(:,1,5) - 1._core_rknd ) > tol ), &
                  'Mean ice fallback changed rime mass' )
      call check( any( abs( samples(:,s,1,6) / means(:,1,6) - 1._core_rknd ) > tol ), &
                  'Mean ice fallback changed rime volume' )
    enddo
    meta%iiPDF_ri = 7

    ! A sample with excess ice can return a valid end state but leave invalid
    ! bulk rime after averaging increments: (qi,qir)=(10,4), sampled input=(20,8),
    ! sampled end=(10,6) gives the bulk state (0,2), in arbitrary mass units.
    ! Include this state along with over-rimed ice and missing-volume cases.
    allocate( meta%l_mix_rat_hm(nh), meta%l_frozen_hm(nh) )
    meta%l_mix_rat_hm = [.true., .false., .true., .false., .false., .false.]
    meta%l_frozen_hm = [.false., .false., .true., .true., .true., .true.]
    means(:,3,5) = 2.e-5_core_rknd
    means(:,3,6) = means(:,3,5) / density
    means(:,5,5) = 2._core_rknd * means(:,5,3)
    means(:,5,6) = means(:,5,5) / density
    call clip_hydromet_conc_mvr( nc, nz, nh, meta, means, & ! In
                                 clipped )                  ! Out
    call check( any( clipped(:,1,5:6) /= means(:,1,5:6) ), 'Valid bulk rime changed' )
    call check( any( clipped(:,:,3) /= means(:,:,3) ), 'Clipping changed total ice water' )
    call check( any( clipped(:,:,1) /= means(:,:,1) ), 'Clipping changed rain water' )
    call check( any( clipped(:,:,5) > clipped(:,:,3) ), 'Bulk rime exceeds total ice' )
    call check( any( clipped(:,:,5:6) < 0._core_rknd ), 'Negative bulk rime' )
    call check( any( clipped(:,3,5:6) /= 0._core_rknd ), 'Orphan bulk rime retained' )
    call check( any( clipped(:,4,5:6) /= 0._core_rknd ), 'Undefined bulk density retained' )
    do i = 1, nc
      call check( clipped(i,5,5) /= means(i,5,3), 'Bulk rime bound failed' )
      call check( abs( clipped(i,5,5) / clipped(i,5,6) / density(i) - 1._core_rknd ) > tol, &
                  'Bulk correction changed rime density' )
    enddo

    ! The volume field may appear before the mass field in the input array.
    swapped_means = means
    swapped_means(:,:,5) = means(:,:,6)
    swapped_means(:,:,6) = means(:,:,5)
    meta%iiqirim = 6
    meta%iibirim = 5
    call clip_hydromet_conc_mvr( nc, nz, nh, meta, swapped_means, & ! In
                                 swapped_clipped )                  ! Out
    call check( any( swapped_clipped(:,:,6) /= clipped(:,:,5) ), &
                'Rime mass depends on field order' )
    call check( any( swapped_clipped(:,:,5) /= clipped(:,:,6) ), &
                'Rime volume depends on field order' )
    call check( any( swapped_clipped(:,:,1:4) /= clipped(:,:,1:4) ), &
                'Other fields depend on rime order' )

    ! Without P3 rime metadata, retain the existing mean-field fallback.
    meta%iiqirim = -1; meta%iibirim = -1
    original_points = points
    call construct_p3_rime_samples( nz, nc, ns, np, nh, meta, means, & ! In
                                    points )                           ! InOut
    call check( any( points /= original_points ), 'Non-P3 sample preparation changed fields' )

    call copy_X_nl_into_hydromet_all_pts( nz, np, ns, nc, points, nh, meta, means, & ! In
                                          samples, nuclei )                          ! Out
    do s = 1, ns
      call check( any( samples(:,s,:,5:6) /= means(:,:,5:6) ), 'Non-P3 fallback changed' )
    enddo

    if ( p3_silhs_rime_unit_test == 0 ) then
      print *, 'Success!'
    end if

    return

  contains

    !===========================================================================
    subroutine check( l_failed, message )

      ! Description:
      !   Record a failed check and print its message.
      !-------------------------------------------------------------------------

      implicit none

      logical, intent(in) :: &
        l_failed ! The checked condition failed

      character(len=*), intent(in) :: &
        message ! Description of the failed check

      if ( l_failed ) then
        write(fstderr,*) 'P3 SILHS rime test failed: ', message
        p3_silhs_rime_unit_test = 1
      end if

    end subroutine check

  end function p3_silhs_rime_unit_test

end module p3_silhs_rime_test
