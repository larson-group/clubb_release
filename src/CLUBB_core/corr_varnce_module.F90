!-----------------------------------------------------------------------
!$Id$
!-------------------------------------------------------------------------------
module corr_varnce_module

  use clubb_precision, only: &
      core_rknd

  use error_code, only: &
        clubb_at_least_debug_level, &   ! Procedure
        clubb_fatal_error, &            ! Constant
        err_code                        ! Error indicator

  implicit none

  type hmp2_ip_on_hmm2_ip_ratios_type

    ! Prescribed parameters for hydrometeor values of <hm|_ip'^2> / <hm|_ip>^2,
    ! where  <hm|_ip'^2>  is the in-precip. variance of the hydrometeor and
    !        <hm|_ip>     is the in-precip. mean of the hydrometeor
    ! 
    ! These values are dependent on the horizontal grid spacing of the run, and are calculated
    ! using a slope and intercept corresponding to each hydrometer.
    real( kind = core_rknd ) :: &
      rr = 1.0_core_rknd, & ! For rain water mixing ratio [-]
      Nr = 1.0_core_rknd, & ! For rain drop concentration [-]
      ri = 1.0_core_rknd, & ! For ice mixing ratio        [-]
      Ni = 1.0_core_rknd, & ! For ice concentration       [-]
      rs = 1.0_core_rknd, & ! For snow mixing ratio       [-]
      Ns = 1.0_core_rknd, & ! For snow concentration      [-]
      rg = 1.0_core_rknd, & ! For graupel mixing ratio    [-]
      Ng = 1.0_core_rknd    ! For graupel concentration   [-]

  end type hmp2_ip_on_hmm2_ip_ratios_type
!$omp declare mapper (hmp2_ip_on_hmm2_ip_ratios_type::x) map ( &
!$omp  x%ng &
!$omp )

  ! These slopes and intercepts below are used to calculate the hmp2_ip_on_hmm2_ip_ratios_type
  ! values that are defined above. This functionality is described by equations 8, 10, & 11
  ! in ``Parameterization of the Spatial Variability of Rain for Large-Scale Models and 
  ! Remote Sensing`` Lebo, et al, October 2015
  ! https://journals.ametsoc.org/doi/pdf/10.1175/JAMC-D-15-0066.1
  ! see clubb:ticket:830 for more detail
  ! 
  ! hmp2_ip_on_hmm2_ip(iirr) = hmp2_ip_on_hmm2_ip_intrcpt%rr + &
  !                            hmp2_ip_on_hmm2_ip_slope%rr * max( host_dx, host_dy )
  ! 
  ! In Lebo et al. the suggested values were
  !     slope = 2.12e-5 [1/m]
  !     intercept = 0.54 [-]
  ! 
  ! In CLUBB standalone, these parameters can be set based on the value for a
  ! given case in the CASE_model.in file.
  type hmp2_ip_on_hmm2_ip_slope_type

    real( kind = core_rknd ) :: &
      rr = 2.12e-5_core_rknd, & ! For rain water mixing ratio [1/m]
      Nr = 2.12e-5_core_rknd, & ! For rain drop concentration [1/m]
      ri = 2.12e-5_core_rknd, & ! For ice mixing ratio        [1/m]
      Ni = 2.12e-5_core_rknd, & ! For ice concentration       [1/m]
      rs = 2.12e-5_core_rknd, & ! For snow mixing ratio       [1/m]
      Ns = 2.12e-5_core_rknd, & ! For snow concentration      [1/m]
      rg = 2.12e-5_core_rknd, & ! For graupel mixing ratio    [1/m]
      Ng = 2.12e-5_core_rknd    ! For graupel concentration   [1/m]

  end type hmp2_ip_on_hmm2_ip_slope_type
!$omp declare mapper (hmp2_ip_on_hmm2_ip_slope_type::x) map ( &
!$omp  x%ng &
!$omp )

  type hmp2_ip_on_hmm2_ip_intrcpt_type

    real( kind = core_rknd ) :: &
      rr = 0.54_core_rknd, & ! For rain water mixing ratio [-]
      Nr = 0.54_core_rknd, & ! For rain drop concentration [-]
      ri = 0.54_core_rknd, & ! For ice mixing ratio        [-]
      Ni = 0.54_core_rknd, & ! For ice concentration       [-]
      rs = 0.54_core_rknd, & ! For snow mixing ratio       [-]
      Ns = 0.54_core_rknd, & ! For snow concentration      [-]
      rg = 0.54_core_rknd, & ! For graupel mixing ratio    [-]
      Ng = 0.54_core_rknd    ! For graupel concentration   [-]

  end type hmp2_ip_on_hmm2_ip_intrcpt_type
!$omp declare mapper (hmp2_ip_on_hmm2_ip_intrcpt_type::x) map ( &
!$omp  x%ng &
!$omp )



  type hm_metadata_type

    ! Variables
    ! Microphysics mixing ratios
    integer :: &
      iirr,   & ! Hydrometeor array index for rain water mixing ratio, rr
      iirs,   & ! Hydrometeor array index for snow mixing ratio, rs
      iiri,   & ! Hydrometeor array index for ice mixing ratio, ri
      iirg      ! Hydrometeor array index for graupel mixing ratio, rg

    ! Microphysics concentrations
    integer :: &
      iiNr,   & ! Hydrometeor array index for rain drop concentration, Nr
      iiNs,   & ! Hydrometeor array index for snow concentration, Ns
      iiNi,   & ! Hydrometeor array index for ice concentration, Ni
      iiNg      ! Hydrometeor array index for graupel concentration, Ng

    ! Logical fields
    logical, dimension(:), allocatable :: &
      l_frozen_hm, & ! if true, then the hydrometeor is frozen; otherwise liquid
      l_mix_rat_hm   ! if true, then the quantity is a hydrometeor mixing ratio

    character(len=10), dimension(:), allocatable :: & 
      hydromet_list

    real( kind = core_rknd ), dimension(:), allocatable :: &
      hydromet_tol    ! Tolerance values for all hydrometeors    [units vary]

    ! Latin hypercube indices / Correlation array indices
    integer :: &
      iiPDF_chi = -1, &
      iiPDF_eta = -1, &
      iiPDF_w   = -1

    integer :: &
     iiPDF_rr = -1, &
     iiPDF_rs = -1, &
     iiPDF_ri = -1, &
     iiPDF_rg = -1

    integer :: &
     iiPDF_Nr  = -1, &
     iiPDF_Ns  = -1, &
     iiPDF_Ni  = -1, &
     iiPDF_Ng  = -1, &
     iiPDF_Ncn = -1

    real( kind = core_rknd ), dimension(:), allocatable :: &
       hmp2_ip_on_hmm2_ip

    ! Prescribed parameter for <N_cn'^2> / <N_cn>^2.
    ! NOTE: In the case that l_const_Nc_in_cloud is true, Ncn is constant
    !       throughout the entire grid box, so the parameter below should be
    !       ignored.
    real( kind = core_rknd ) :: &
      Ncnp2_on_Ncnm2 = 1.0_core_rknd   ! Prescribed ratio <N_cn'^2> / <N_cn>^2 [-]

  end type hm_metadata_type 
!$omp declare mapper (hm_metadata_type::x) map ( &
!$omp  x%iirg &
!$omp , x%iing &
!$omp , x%l_mix_rat_hm &
!$omp , x%hydromet_list &
!$omp , x%hydromet_tol &
!$omp , x%iipdf_chi &
!$omp , x%iipdf_eta &
!$omp , x%iipdf_w &
!$omp , x%iipdf_rr &
!$omp , x%iipdf_rs &
!$omp , x%iipdf_ri &
!$omp , x%iipdf_rg &
!$omp , x%iipdf_nr &
!$omp , x%iipdf_ns &
!$omp , x%iipdf_ni &
!$omp , x%iipdf_ng &
!$omp , x%iipdf_ncn &
!$omp , x%hmp2_ip_on_hmm2_ip &
!$omp , x%ncnp2_on_ncnm2 &
!$omp )

  ! This "renaming" is used to shorten the matrix declarations below.
  integer, parameter :: c = core_rknd

  integer, parameter, public :: &
    d_var_total = 12 ! Size of the default correlation arrays

  real( kind = core_rknd ), public, dimension(d_var_total,d_var_total) :: &
    corr_array_n_cloud_def = reshape( &
  !  chi   eta    w      Ncn     rr      Nr      ri      Ni      rs      Ns      rg      Ng
  (/1._c,-.6_c, .09_c , .09_c , .788_c, .675_c, .240_c, .222_c, .240_c, .222_c, .240_c, .222_c, & ! chi
    0._c, 1._c, .027_c, .027_c, .114_c, .115_c,-.029_c, .093_c, .022_c, .013_c, 0._c  , 0._c  , & ! eta
    0._c, 0._c, 1._c  , .34_c , .315_c, .270_c, .120_c, .167_c, 0._c  , 0._c  , 0._c  , 0._c  , & ! w
    0._c, 0._c, 0._c  , 1._c  , 0._c  , 0._c  , .464_c, .320_c, .168_c, .232_c, 0._c  , 0._c  , & ! Ncn
    0._c, 0._c, 0._c  , 0._c  , 1._c  , .821_c, 0._c  , 0._c  , .173_c, .164_c, .319_c, .308_c, & ! rr
    0._c, 0._c, 0._c  , 0._c  , 0._c  , 1._c  , .152_c, .143_c, 0._c  , 0._c  , .285_c, .273_c, & ! Nr
    0._c, 0._c, 0._c  , 0._c  , 0._c  , 0._c  , 1._c  , .758_c, .585_c, .571_c, .379_c, .363_c, & ! ri
    0._c, 0._c, 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 1._c  , .571_c, .550_c, .363_c, .345_c, & ! Ni
    0._c, 0._c, 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 1._c  , .758_c, .485_c, .470_c, & ! rs
    0._c, 0._c, 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 1._c  , .470_c, .450_c, & ! Ns
    0._c, 0._c, 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 1._c  , .758_c, & ! rg
    0._c, 0._c, 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 1._c/), & ! Ng
    (/d_var_total,d_var_total/) ), &

    corr_array_n_below_def = reshape( &
  !  chi   eta    w      Ncn     rr      Nr      ri      Ni      rs      Ns      rg      Ng
  (/1._c, .3_c, .09_c , .09_c , .788_c, .675_c, .240_c, .222_c, .240_c, .222_c, .240_c, .222_c, &! chi
    0._c, 1._c, .027_c, .027_c, .114_c, .115_c,-.029_c, .093_c, .022_c, .013_c, 0._c  , 0._c  , &! eta
    0._c, 0._c, 1._c  , .34_c , .315_c, .270_c, .120_c, .167_c, 0._c  , 0._c  , 0._c  , 0._c  , &! w
    0._c, 0._c, 0._c  , 1._c  , 0._c  , 0._c  , .464_c, .320_c, .168_c, .232_c, 0._c  , 0._c  , &! Ncn
    0._c, 0._c, 0._c  , 0._c  , 1._c  , .821_c, 0._c  , 0._c  , .173_c, .164_c, .319_c, .308_c, &! rr
    0._c, 0._c, 0._c  , 0._c  , 0._c  , 1._c  , .152_c, .143_c, 0._c  , 0._c  , .285_c, .273_c, &! Nr
    0._c, 0._c, 0._c  , 0._c  , 0._c  , 0._c  , 1._c  , .758_c, .585_c, .571_c, .379_c, .363_c, &! ri
    0._c, 0._c, 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 1._c  , .571_c, .550_c, .363_c, .345_c, &! Ni
    0._c, 0._c, 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 1._c  , .758_c, .485_c, .470_c, &! rs
    0._c, 0._c, 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 1._c  , .470_c, .450_c, &! Ns
    0._c, 0._c, 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 1._c  , .758_c, &! rg
    0._c, 0._c, 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 0._c  , 1._c/), &! Ng
    (/d_var_total,d_var_total/) )

  private

  public :: hmp2_ip_on_hmm2_ip_ratios_type, &
            hmp2_ip_on_hmm2_ip_slope_type, &
            hmp2_ip_on_hmm2_ip_intrcpt_type, &
            read_correlation_matrix, &
            setup_corr_varnce_array, &
            assert_corr_symmetric, print_corr_matrix, &
            hm_metadata_type, &
            init_pdf_hydromet_arrays


  private :: get_corr_var_index, def_corr_idx


  contains

  !-----------------------------------------------------------------------------
  function def_corr_idx( iiPDF_x, hm_metadata ) result(ii_def_corr)

    ! Description:
    !   Map from a iiPDF index to the corresponding index in the default 
    !   correlation arrays.
    !-----------------------------------------------------------------------------

    implicit none

    ! Constant Parameters

    ! Indices that represent the order in the default corr arrays
    ! (chi (old s), eta (old t), w, Ncn, rr, Nr, ri, Ni, rs, Ns, rg, Ng)
    integer, parameter :: &
    ii_chi = 1, &
    ii_eta = 2, &
    ii_w = 3,   &
    ii_Ncn = 4, &
    ii_rr = 5,  &
    ii_Nr = 6,  &
    ii_ri = 7,  &
    ii_Ni = 8,  &
    ii_rs = 9,  &
    ii_Ns = 10, &
    ii_rg = 11, &
    ii_Ng = 12

    ! Input Variables

    integer, intent(in) :: iiPDF_x

    ! Return Variable

    integer :: ii_def_corr

    ! ---- Local Variables ----
    type (hm_metadata_type) :: &
      hm_metadata

    ! ---- Begin Code ----

    ii_def_corr = -1

      if (iiPDF_x == hm_metadata%iiPDF_chi) then
         ii_def_corr = ii_chi

      elseif (iiPDF_x == hm_metadata%iiPDF_eta) then
        ii_def_corr = ii_eta

      elseif (iiPDF_x == hm_metadata%iiPDF_w) then
        ii_def_corr = ii_w

      elseif (iiPDF_x == hm_metadata%iiPDF_Ncn) then
        ii_def_corr = ii_Ncn

      elseif (iiPDF_x == hm_metadata%iiPDF_rr) then
        ii_def_corr = ii_rr

      elseif (iiPDF_x == hm_metadata%iiPDF_Nr) then
        ii_def_corr = ii_Nr

      elseif (iiPDF_x == hm_metadata%iiPDF_ri) then
        ii_def_corr = ii_ri

      elseif (iiPDF_x == hm_metadata%iiPDF_Ni) then
        ii_def_corr = ii_Ni

      elseif (iiPDF_x == hm_metadata%iiPDF_rs) then
        ii_def_corr = ii_rs

      elseif (iiPDF_x == hm_metadata%iiPDF_Ns) then
        ii_def_corr = ii_Ns

      elseif (iiPDF_x == hm_metadata%iiPDF_rg) then
        ii_def_corr = ii_rg

      elseif (iiPDF_x == hm_metadata%iiPDF_Ng) then
        ii_def_corr = ii_Ng

      endif

  end function def_corr_idx

  !-----------------------------------------------------------------------------
  subroutine set_corr_arrays_to_default( pdf_dim, hm_metadata, &
                                         corr_array_n_cloud, corr_array_n_below ) 

    ! Description:
    !   If there are no corr_array.in files for the current case, default 
    !   correlations are used. 
    !-----------------------------------------------------------------------------
  
    use constants_clubb, only: &
        zero, &
        one

    implicit none

    !------------------------ Input Variables ------------------------
    integer, intent(in) :: &
      pdf_dim

    type (hm_metadata_type), intent(in) :: &
      hm_metadata

    !------------------------ Output Variables ------------------------
    real( kind = core_rknd ), intent(out), dimension(pdf_dim,pdf_dim) :: &
      corr_array_n_cloud, &
      corr_array_n_below

    !------------------------ Local Variables ------------------------
    integer :: i, j, idx_i, idx_j ! Loop iterators

    !------------------------ Begin Code ------------------------

    corr_array_n_cloud = zero
    corr_array_n_below = zero

    do i = 1, pdf_dim
       corr_array_n_cloud(i,i) = one
       corr_array_n_below(i,i) = one
    enddo

    do i = 1, pdf_dim-1
      do j = i+1, pdf_dim

        idx_i = def_corr_idx(i,hm_metadata)
        idx_j = def_corr_idx(j,hm_metadata)

        if ( idx_i > idx_j ) then
          corr_array_n_cloud(j, i) = corr_array_n_cloud_def( idx_i, idx_j )
          corr_array_n_below(j, i) = corr_array_n_below_def( idx_i, idx_j )
        else
          corr_array_n_cloud(j, i) = corr_array_n_cloud_def( idx_j, idx_i )
          corr_array_n_below(j, i) = corr_array_n_below_def( idx_j, idx_i )
        end if
      end do
    end do

    return

  end subroutine set_corr_arrays_to_default


  !-----------------------------------------------------------------------------
  subroutine read_correlation_matrix( iunit, input_file, &
                                      pdf_dim, hm_metadata, &
                                      corr_array_n )

    ! Description:
    !   Reads a correlation variance array from a file and stores it in an array.
    !-----------------------------------------------------------------------------

    use input_reader, only: &
        one_dim_read_var, & ! Variable(s)
        read_one_dim_file, deallocate_one_dim_vars, count_columns ! Procedure(s)

    use matrix_operations, only: set_lower_triangular_matrix ! Procedure(s)

    use constants_clubb, only: fstderr ! Variable(s)

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variable(s)
    integer, intent(in) :: &
      iunit, &    ! File I/O unit
      pdf_dim! number of variables in the array

    type (hm_metadata_type), intent(in) :: &
      hm_metadata

    character(len=*), intent(in) :: input_file ! Path to the file

    ! Input/Output Variable(s)
    real( kind = core_rknd ), dimension(pdf_dim,pdf_dim), intent(inout) :: &
      corr_array_n ! Normal space correlation array

    ! Local Variable(s)

    type(one_dim_read_var), allocatable, dimension(:) :: &
      retVars ! stores the variables read in from the corr_varnce.in file

    integer ::   &
      var_index1,    & ! variable index
      var_index2,    & ! variable index
      nCols,         & ! the number of columns in the file
      i, j         ! Loop index


    !--------------------------- BEGIN CODE -------------------------

    nCols = count_columns( iunit, input_file )

    ! Allocate all arrays based on pdf_dim
    allocate( retVars(1:nCols) )

    ! Initializing to zero means that correlations we don't have are assumed to be 0.
    corr_array_n(:,:) = 0.0_core_rknd

    ! Set main diagonal to 1
    do i=1, pdf_dim
      corr_array_n(i,i) = 1.0_core_rknd
    end do

    ! Read the values from the specified file
    call read_one_dim_file( iunit, nCols, input_file, & ! intent(in)
                            retVars ) ! intent(out)

    if( size( retVars(1)%values ) /= nCols ) then
      write(fstderr, *) "Correlation matrix must have an equal number of rows and cols in file ", &
            input_file
      error stop "Bad data in correlation file."
    end if

    ! Start at 2 because the first index is always just 1.0 in the first row
    ! and the rest of the rows are ignored
    do i=2, nCols
      var_index1 = get_corr_var_index( retVars(i)%name, hm_metadata )
      if( var_index1 > -1 ) then
        do j=1, (i-1)
          var_index2 = get_corr_var_index( retVars(j)%name, hm_metadata )
          if( var_index2 > -1 ) then
            call set_lower_triangular_matrix &
                 ( pdf_dim, var_index1, var_index2, retVars(i)%values(j), & ! intent(in)
                   corr_array_n ) ! intent(inout)
          end if
        end do
      end if
    end do

    call deallocate_one_dim_vars( nCols, & ! intent(in)
                                  retVars ) ! intent(inout)

    return
  end subroutine read_correlation_matrix

  !--------------------------------------------------------------------------
  function get_corr_var_index( var_name, hm_metadata ) result( i )

    ! Definition:
    !   Returns the index for a variable based on its name.
    !--------------------------------------------------------------------------

    implicit none

    character(len=*), intent(in) :: var_name ! The name of the variable

    type (hm_metadata_type), intent(in) :: &
      hm_metadata

    ! Output variable
    integer :: i

    !------------------ BEGIN CODE -----------------------------
    i = -1

    select case( trim(var_name) )

    case( "chi" )
      i = hm_metadata%iiPDF_chi

    case( "eta" )
      i = hm_metadata%iiPDF_eta

    case( "w" )
      i = hm_metadata%iiPDF_w

    case( "Ncn" )
      i = hm_metadata%iiPDF_Ncn

    case( "rr" )
      i = hm_metadata%iiPDF_rr

    case( "Nr" )
      i = hm_metadata%iiPDF_Nr

    case( "ri" )
      i = hm_metadata%iiPDF_ri

    case( "Ni" )
      i = hm_metadata%iiPDF_Ni

    case( "rs" )
      i = hm_metadata%iiPDF_rs

    case( "Ns" )
      i = hm_metadata%iiPDF_Ns
        
    case( "rg" )
      i = hm_metadata%iiPDF_rg

    case( "Ng" )
      i = hm_metadata%iiPDF_Ng

    end select

    return

  end function get_corr_var_index

  !================================================================================================
  ! subroutine init_pdf_hydromet_arrays
  ! 
  ! DESCRIPTION: 
  !     This subroutine intializes the hydromet arrays(iirr, iiNr, etc.) to the values specified by
  !     the input arguements, this determines which hyrometeors are to be used by the microphysics
  !     scheme. It also sets up the corresponding pdf and hydromet arrays, and calculates the 
  !     subgrid variance ratio for each hydrometeor.
  ! 
  ! OPTIONAL FUNCTIONALITY:
  !     The subgrid variance ratio for each hydrometeor is calculated based on the grid spacing 
  !     defined by the host model. The calculation is a linear equation defined by a slope and
  !     intercept, each of which may or may not be passed in to this subroutine. If the slope
  !     and/or intercept are not passed in through the arguement list the default values, which 
  !     are set in the corresponding type definitions, will be used. Otherwise the values
  !     specified by the aruements will be used.
  ! 
  ! NOTES: 
  !     'hmp2_ip_on_hmm2_ip_slope_in' is of type 'hmp2_ip_on_hmm2_ip_slope_type' and
  !     'hmp2_ip_on_hmm2_ip_intrcpt_in' is of type 'hmp2_ip_on_hmm2_ip_intrcpt_in', both of which 
  !     are deinfed in corr_vrnce_module.F90, and made public through this API.
  ! 
  !     If full control over the hydrometeor variance ratios is desired, pass in slopes that are
  !     initialized to 0.0, this causes the ratios to no longer depend on the grid spacing. Then
  !     pass in the intercepts set to the values of the desired ratios.
  ! 
  ! ARGUEMENTS:
  !     host_dx (real) - Horizontal grid spacings
  !     host_dy (real)
  ! 
  !     hydromet_dim (integer) - Number of enabled hydrometeors
  ! 
  !         Each of these is an index value corresponding to a hydrometeor,
  !         used to index the hydrometeor arrays. Each index has to be unqiue
  !         for each different hyrometeor that is enabled. Setting one of these
  !         indices to -1 disables that hydrometeor
  !     iirr (integer) - Index of rain water mixing ratio
  !     iiri (integer) - Index of rain drop concentration
  !     iirs (integer) - Index of ice mixing ratio
  !     iirg (integer) - Index of ice crystal concentration
  !     iiNr (integer) - Index of snow mixing ratio
  !     iiNi (integer) - Index of snow flake concentration
  !     iiNs (integer) - Index of graupel mixing ratio
  !     iiNg (integer) - Index of graupel concentration
  ! 
  !     hmp2_ip_on_hmm2_ip_slope_in (hmp2_ip_on_hmm2_ip_slope_type) - Custom slope values
  !     hmp2_ip_on_hmm2_ip_intrcpt_in (hmp2_ip_on_hmm2_ip_intrcpt_type) - Custom intercept values
  ! 
  !================================================================================================
  subroutine init_pdf_hydromet_arrays( host_dx, host_dy, hydromet_dim,  & ! intent(in)
                                       iirr, iiNr, iiri, iiNi,          & ! intent(in)
                                       iirs, iiNs, iirg, iiNg,          & ! intent(in)
                                       Ncnp2_on_Ncnm2,                  & ! intent(in)
                                       hm_metadata, pdf_dim,         & ! intent(out)
                                       hmp2_ip_on_hmm2_ip_slope_in,     & ! optional(in)
                                       hmp2_ip_on_hmm2_ip_intrcpt_in    ) ! optional(in)
  ! Description:
  ! 
  ! Initialization for the hydromet arrays. How the arrays are initialized is
  ! determined by how the hydromet indicies (iirr, iiNr, etc.) have been set,
  ! which are determined by the microphysics scheme.
  !-------------------------------------------------------------------------------

    use constants_clubb, only: &
        rr_tol, &   ! Constants, tolerances for each hydrometeor
        ri_tol, &
        rs_tol, &
        rg_tol, &
        Nr_tol, &
        Ni_tol, &
        Ns_tol, &
        Ng_tol

    implicit none

    !------------------------ Input Variables ------------------------
    real( kind = core_rknd ), intent(in) :: &
      host_dx, host_dy  ! Horizontal grid spacing, defined by host model [m]

    integer, intent(in) :: &
      hydromet_dim, & ! Total number of hydrometeor species.
      iirr,         & ! Index of rain water mixing ratio
      iiNr,         & ! Index of rain drop concentration
      iiri,         & ! Index of ice mixing ratio
      iiNi,         & ! Index of ice crystal concentration
      iirs,         & ! Index of snow mixing ratio
      iiNs,         & ! Index of snow flake concentration
      iirg,         & ! Index of graupel mixing ratio
      iiNg            ! Index of graupel concentration

    real( kind = core_rknd ), intent(in) :: &
      Ncnp2_on_Ncnm2

    !------------------------ InOut Variables ------------------------
    type (hm_metadata_type), intent(out) :: &
      hm_metadata

    !------------------------ Output Variables ------------------------
    integer, intent(out) :: &
      pdf_dim

    !------------------------ Optional Input Variables ------------------------

    ! Used to overwrite default values of slope and intercept
    type( hmp2_ip_on_hmm2_ip_slope_type ), optional, intent(in) :: &
        hmp2_ip_on_hmm2_ip_slope_in    ! Custom slopes to overwrite defaults [1/m]
      
    type( hmp2_ip_on_hmm2_ip_intrcpt_type ), optional, intent(in) :: &
        hmp2_ip_on_hmm2_ip_intrcpt_in   ! Custom intercepts to overwrite defaults [-]
    
    !------------------------ Local Variables ------------------------    
    integer :: i, pdf_count

    type( hmp2_ip_on_hmm2_ip_slope_type ) :: &
        hmp2_ip_on_hmm2_ip_slope
      
    type( hmp2_ip_on_hmm2_ip_intrcpt_type ) :: &
        hmp2_ip_on_hmm2_ip_intrcpt

    !------------------------ Begin Code ------------------------

    ! Save indices
    hm_metadata%iirr = iirr
    hm_metadata%iiri = iiri
    hm_metadata%iirs = iirs
    hm_metadata%iirg = iirg
    hm_metadata%iiNr = iiNr
    hm_metadata%iiNi = iiNi
    hm_metadata%iiNs = iiNs
    hm_metadata%iiNg = iiNg

    ! Overwrite default Ncnp2_on_Ncnm2
    hm_metadata%Ncnp2_on_Ncnm2 = Ncnp2_on_Ncnm2

    !-----------------------------------------------------------
    ! Calculate the subgrid variances of the hydrometeors
    !-----------------------------------------------------------

    ! If slope and intercept are present in call, then overwrite default values
    if ( present( hmp2_ip_on_hmm2_ip_slope_in ) ) then
        hmp2_ip_on_hmm2_ip_slope = hmp2_ip_on_hmm2_ip_slope_in
    end if

    if ( present( hmp2_ip_on_hmm2_ip_intrcpt_in ) ) then
        hmp2_ip_on_hmm2_ip_intrcpt = hmp2_ip_on_hmm2_ip_intrcpt_in
    end if

    allocate( hm_metadata%hmp2_ip_on_hmm2_ip(hydromet_dim) )
    hm_metadata%hmp2_ip_on_hmm2_ip = 0.0_core_rknd

    if ( iirr > 0 ) then
       hm_metadata%hmp2_ip_on_hmm2_ip(iirr) &
                                = hmp2_ip_on_hmm2_ip_intrcpt%rr  &
                                  + hmp2_ip_on_hmm2_ip_slope%rr * max( host_dx, host_dy )
    endif

    if ( iirs > 0 ) then
      hm_metadata%hmp2_ip_on_hmm2_ip(iirs) &
                                  = hmp2_ip_on_hmm2_ip_intrcpt%rs &
                                  + hmp2_ip_on_hmm2_ip_slope%rs * max( host_dx, host_dy )
    endif

    if ( iiri > 0 ) then
      hm_metadata%hmp2_ip_on_hmm2_ip(iiri) &
                                  = hmp2_ip_on_hmm2_ip_intrcpt%ri &
                                  + hmp2_ip_on_hmm2_ip_slope%ri * max( host_dx, host_dy )
    endif

    if ( iirg > 0 ) then
      hm_metadata%hmp2_ip_on_hmm2_ip(iirg) &
                                  = hmp2_ip_on_hmm2_ip_intrcpt%rg &
                                  + hmp2_ip_on_hmm2_ip_slope%rg * max( host_dx, host_dy )
    endif

    if ( iiNr > 0 ) then
      hm_metadata%hmp2_ip_on_hmm2_ip(iiNr) &
                                  = hmp2_ip_on_hmm2_ip_intrcpt%Nr &
                                  + hmp2_ip_on_hmm2_ip_slope%Nr * max( host_dx, host_dy )
    endif

    if ( iiNs > 0 ) then
      hm_metadata%hmp2_ip_on_hmm2_ip(iiNs) &
                                  = hmp2_ip_on_hmm2_ip_intrcpt%Ns &
                                  + hmp2_ip_on_hmm2_ip_slope%Ns * max( host_dx, host_dy )
    endif

    if ( iiNi > 0 ) then
      hm_metadata%hmp2_ip_on_hmm2_ip(iiNi) &
                                  = hmp2_ip_on_hmm2_ip_intrcpt%Ni &
                                  + hmp2_ip_on_hmm2_ip_slope%Ni * max( host_dx, host_dy )
    endif

    if ( iiNg > 0 ) then
      hm_metadata%hmp2_ip_on_hmm2_ip(iiNg) &
                                  = hmp2_ip_on_hmm2_ip_intrcpt%Ng &
                                  + hmp2_ip_on_hmm2_ip_slope%Ng * max( host_dx, host_dy )
    endif


    !-----------------------------------------------------------
    ! Set up predictive precipitating hydrometeor arrays.
    !-----------------------------------------------------------
    allocate( hm_metadata%hydromet_list(hydromet_dim) )
    allocate( hm_metadata%hydromet_tol(hydromet_dim) )
    allocate( hm_metadata%l_mix_rat_hm(hydromet_dim) )
    allocate( hm_metadata%l_frozen_hm(hydromet_dim) )

    if ( iirr > 0 ) then
       ! The microphysics scheme predicts rain water mixing ratio, rr.
       hm_metadata%hydromet_list(iirr)      = "rrm"
       hm_metadata%l_mix_rat_hm(iirr)       = .true.
       hm_metadata%l_frozen_hm(iirr)        = .false.
       hm_metadata%hydromet_tol(iirr)       = rr_tol
    endif

    if ( iiri > 0 ) then
       ! The microphysics scheme predicts ice mixing ratio, ri.
       hm_metadata%hydromet_list(iiri)      = "rim"
       hm_metadata%l_mix_rat_hm(iiri)       = .true.
       hm_metadata%l_frozen_hm(iiri)        = .true.
       hm_metadata%hydromet_tol(iiri)       = ri_tol
    endif

    if ( iirs > 0 ) then
       ! The microphysics scheme predicts snow mixing ratio, rs.
       hm_metadata%hydromet_list(iirs)      = "rsm"
       hm_metadata%l_mix_rat_hm(iirs)       = .true.
       hm_metadata%l_frozen_hm(iirs)        = .true.
       hm_metadata%hydromet_tol(iirs)       = rs_tol
    endif

    if ( iirg > 0 ) then
       ! The microphysics scheme predicts graupel mixing ratio, rg.
       hm_metadata%hydromet_list(iirg)      = "rgm"
       hm_metadata%l_mix_rat_hm(iirg)       = .true.
       hm_metadata%l_frozen_hm(iirg)        = .true.
       hm_metadata%hydromet_tol(iirg)       = rg_tol
    endif

    if ( iiNr > 0 ) then
       ! The microphysics scheme predicts rain drop concentration, Nr.
       hm_metadata%hydromet_list(iiNr)      = "Nrm"
       hm_metadata%l_frozen_hm(iiNr)        = .false.
       hm_metadata%l_mix_rat_hm(iiNr)       = .false.
       hm_metadata%hydromet_tol(iiNr)       = Nr_tol
    endif

    if ( iiNi > 0 ) then
       ! The microphysics scheme predicts ice concentration, Ni.
       hm_metadata%hydromet_list(iiNi)      = "Nim"
       hm_metadata%l_mix_rat_hm(iiNi)       = .false.
       hm_metadata%l_frozen_hm(iiNi)        = .true.
       hm_metadata%hydromet_tol(iiNi)       = Ni_tol
    endif

    if ( iiNs > 0 ) then
       ! The microphysics scheme predicts snowflake concentration, Ns.
       hm_metadata%hydromet_list(iiNs)      = "Nsm"
       hm_metadata%l_mix_rat_hm(iiNs)       = .false.
       hm_metadata%l_frozen_hm(iiNs)        = .true.
       hm_metadata%hydromet_tol(iiNs)       = Ns_tol
    endif

    if ( iiNg > 0 ) then
       ! The microphysics scheme predicts graupel concentration, Ng.
       hm_metadata%hydromet_list(iiNg)      = "Ngm"
       hm_metadata%l_mix_rat_hm(iiNg)       = .false.
       hm_metadata%l_frozen_hm(iiNg)        = .true.
       hm_metadata%hydromet_tol(iiNg)       = Ng_tol
    endif


    !-----------------------------------------------------------
    ! Set up the PDF indices
    !-----------------------------------------------------------

    hm_metadata%iiPDF_chi = 1 ! Extended liquid water mixing ratio, chi
    hm_metadata%iiPDF_eta = 2 ! 'eta' orthogonal to 'chi'
    hm_metadata%iiPDF_w   = 3 ! vertical velocity
    hm_metadata%iiPDF_Ncn = 4 ! Simplified cloud nuclei concentration or extended Nc.

    pdf_count = hm_metadata%iiPDF_Ncn

    ! Loop over hydrometeors.
    ! Hydrometeor indices in the PDF arrays should be in the same order as
    ! found in the hydrometeor arrays.
    if ( hydromet_dim > 0 ) then

       do i = 1, hydromet_dim, 1

          if ( i == iirr ) then
             pdf_count = pdf_count + 1
             hm_metadata%iiPDF_rr = pdf_count
          endif

          if ( i == iiNr ) then
             pdf_count = pdf_count + 1
             hm_metadata%iiPDF_Nr = pdf_count
          endif

          if ( i == iiri ) then
             pdf_count = pdf_count + 1
             hm_metadata%iiPDF_ri = pdf_count
          endif

          if ( i == iiNi ) then
             pdf_count = pdf_count + 1
             hm_metadata%iiPDF_Ni = pdf_count
          endif

          if ( i == iirs ) then
             pdf_count = pdf_count + 1
             hm_metadata%iiPDF_rs = pdf_count
          endif

          if ( i == iiNs ) then
             pdf_count = pdf_count + 1
             hm_metadata%iiPDF_Ns = pdf_count
          endif

          if ( i == iirg ) then
             pdf_count = pdf_count + 1
             hm_metadata%iiPDF_rg = pdf_count
          endif
        
          if ( i == iiNg ) then
             pdf_count = pdf_count + 1
             hm_metadata%iiPDF_Ng = pdf_count
          endif   

       enddo ! i = 1, hydromet_dim, 1

    endif ! hydromet_dim > 0

    pdf_dim = pdf_count

    return

  end subroutine init_pdf_hydromet_arrays

!===============================================================================
  subroutine setup_corr_varnce_array( input_file_cloud, input_file_below, &
                                      pdf_dim, hm_metadata, iunit, &
                                      l_fix_w_chi_eta_correlations, &
                                      corr_array_n_cloud, corr_array_n_below )

! Description:
!   Setup an array with the x'^2/xm^2 variables on the diagonal and the other
!   elements to be correlations between various variables.

! References:
!   None.
!-------------------------------------------------------------------------------

    use matrix_operations, only: mirror_lower_triangular_matrix ! Procedure

    use constants_clubb, only: &
        fstderr, &  ! Constant(s)
        eps

    implicit none

    ! External
    intrinsic :: max, epsilon, trim


    !------------------------ Input Variables ------------------------
    character(len=*), intent(in) :: &
      input_file_cloud, &    ! Path to the in cloud correlation file
      input_file_below       ! Path to the out of cloud correlation file

    integer, intent(in) :: &
      pdf_dim
      
    type (hm_metadata_type), intent(in) :: &
      hm_metadata

    integer, intent(in) :: &
      iunit ! The file unit

    logical, intent(in) :: &
      l_fix_w_chi_eta_correlations ! Use a fixed correlation for s and t Mellor(chi/eta)

    !------------------------ Output Variables ------------------------
    real( kind = core_rknd ), intent(out), dimension(pdf_dim,pdf_dim) :: &
      corr_array_n_cloud, &
      corr_array_n_below

    !------------------------ Local variables ------------------------
    logical :: l_warning, l_corr_file_1_exist, l_corr_file_2_exist
    integer :: i

    !------------------------ Begin Code ------------------------

    inquire( file = input_file_cloud, exist = l_corr_file_1_exist )
    inquire( file = input_file_below, exist = l_corr_file_2_exist )

    if ( l_corr_file_1_exist .and. l_corr_file_2_exist ) then

       call read_correlation_matrix( iunit, trim( input_file_cloud ), &
                                     pdf_dim, hm_metadata, & ! In
                                     corr_array_n_cloud ) ! Out

       call read_correlation_matrix( iunit, trim( input_file_below ), &
                                     pdf_dim, hm_metadata, & ! In
                                     corr_array_n_below ) ! Out

    else ! Read in default correlation matrices
        
      if ( clubb_at_least_debug_level( 1 ) ) then
        write(fstderr,*) "Warning: "//trim( input_file_cloud )//" was not found! " // &
                         "The default correlation arrays will be used." 
      end if

      call set_corr_arrays_to_default( pdf_dim, hm_metadata, &
                                       corr_array_n_cloud, corr_array_n_below )

    endif
      
    ! Mirror the correlation matrices
    call mirror_lower_triangular_matrix( pdf_dim, & ! intent(in)
                                         corr_array_n_cloud ) ! intent(inout)

    call mirror_lower_triangular_matrix( pdf_dim, & ! intent(in)
                                         corr_array_n_below ) ! intent(inout)
    

    ! Sanity check to avoid confusing non-convergence results.
    if ( clubb_at_least_debug_level( 2 ) ) then

      if ( .not. l_fix_w_chi_eta_correlations .and. hm_metadata%iiPDF_Ncn > 0 ) then
        l_warning = .false.
        do i = 1, pdf_dim
          if ( ( (abs(corr_array_n_cloud(i,hm_metadata%iiPDF_Ncn)) > eps) .or.  &
                 (abs(corr_array_n_below(i,hm_metadata%iiPDF_Ncn))) > eps) .and. &
               i /= hm_metadata%iiPDF_Ncn ) then
            l_warning = .true.
          end if
        end do ! 1..pdf_dim
        if ( l_warning ) then
          write(fstderr,*) "Warning: the specified correlations for chi" &
                           // " (old s) and Ncn are non-zero."
          write(fstderr,*) "The latin hypercube code will not converge to" &
                           // " the analytic solution using these settings."
        end if
       end if ! l_fix_w_chi_eta_correlations .and. iiPDF_Ncn > 0

    end if ! clubb_at_least_debug_level( 2 )


    return

  end subroutine setup_corr_varnce_array

  !-----------------------------------------------------------------------------
  subroutine assert_corr_symmetric( corr_array_n, & ! intent(in)
                                    pdf_dim)        ! intent(in)

    ! Description:
    !   Asserts that corr_matrix(i,j) == corr_matrix(j,i) for all indeces
    !   in the correlation array. If this is not the case, stops the program.
    ! References:
    !   None
    !---------------------------------------------------------------------------

    use constants_clubb, only: fstderr, eps, one ! Constant(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      pdf_dim   ! Number of variables in the correlation array

    real( kind = core_rknd ), dimension(pdf_dim, pdf_dim), &
      intent(in) :: corr_array_n ! Normal space correlation array to be checked

    ! Local Variables

    ! tolerance used for real precision testing
    real( kind = core_rknd ), parameter :: tol = 1.0e-6_core_rknd

    integer :: n_row, n_col ! Indices

    !----- Begin Code -----

    !Do the check
    do n_col = 1, pdf_dim
      do n_row = 1, pdf_dim

        if ( abs(corr_array_n(n_col, n_row) - corr_array_n(n_row, n_col)) > tol ) then
            err_code = clubb_fatal_error
            write(fstderr,*) "in subroutine assert_corr_symmetric: ", &
                             "Correlation array is non symmetric."
            write(fstderr,*) corr_array_n
            return
        end if

        if ( n_col == n_row .and. abs(corr_array_n(n_col, n_row)-one) > eps ) then
            err_code = clubb_fatal_error
            write(fstderr,*) "in subroutine assert_corr_symmetric: ", &
                             "Correlation array is formatted incorrectly."
            write(fstderr,*) corr_array_n
            return
        end if
      end do
    end do

  end subroutine assert_corr_symmetric

  !-----------------------------------------------------------------------------
  subroutine print_corr_matrix( pdf_dim, & ! intent(in)
                                corr_array_n ) ! intent(in)

    ! Description:
    !   Prints the correlation matrix to the console.
    ! References:
    !   None
    !---------------------------------------------------------------------------

    use clubb_precision, only: core_rknd

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      pdf_dim   ! Number of variables in the correlation array

    real( kind = core_rknd ), dimension(pdf_dim, pdf_dim), &
      intent(in) :: corr_array_n ! Normal space correlation array to be printed

    ! Local Variables
    integer :: n, & ! Loop indeces
               m, &
               current_character_index ! keeps track of the position in the string

    character(LEN=72) :: current_line ! The current line to be printed
    character(LEN=10) :: str_array_value

    !----- Begin Code -----

    current_character_index = 0

    do n = 1, pdf_dim
      do m = 1, pdf_dim
        write(str_array_value,'(F5.2)') corr_array_n(m,n)
        current_line = current_line(1:current_character_index)//str_array_value
        current_character_index = current_character_index + 6
      end do
      write(*, *) current_line
      current_line = ""
      current_character_index = 0
    end do

  end subroutine print_corr_matrix
  !-----------------------------------------------------------------------------

end module corr_varnce_module


