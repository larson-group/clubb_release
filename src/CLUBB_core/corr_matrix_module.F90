!-----------------------------------------------------------------------
!$Id$
!-------------------------------------------------------------------------------
module corr_matrix_module

  use clubb_precision, only: &
      core_rknd

  implicit none

  ! Latin hypercube indices / Correlation array indices
  integer, public :: &
    iiPDF_s_mellor = -1, &
    iiPDF_t_mellor = -1, &
    iiPDF_w        = -1
!$omp threadprivate(iiPDF_s_mellor, iiPDF_t_mellor, iiPDF_w)

  integer, public :: &
   iiPDF_rrain    = -1, &
   iiPDF_rsnow    = -1, &
   iiPDF_rice     = -1, &
   iiPDF_rgraupel = -1
!$omp threadprivate(iiPDF_rrain, iiPDF_rsnow, iiPDF_rice, iiPDF_rgraupel)

  integer, public :: &
   iiPDF_Nr       = -1, &
   iiPDF_Nsnow    = -1, &
   iiPDF_Ni       = -1, &
   iiPDF_Ngraupel = -1, &
   iiPDF_Ncn      = -1
!$omp threadprivate(iiPDF_Nr, iiPDF_Nsnow, iiPDF_Ni, iiPDF_Ngraupel, iiPDF_Ncn)

  integer, parameter, public :: &
    d_var_total = 12 ! Size of the default correlation arrays

  integer, public :: &
    d_variables
!$omp threadprivate(d_variables)

  real( kind = core_rknd ), public, dimension(:), allocatable :: &
    sigma2_on_mu2_ip_array_cloud, &
    sigma2_on_mu2_ip_array_below

  real( kind = core_rknd ), public, dimension(:,:), allocatable :: &
    corr_array_cloud, &
    corr_array_below
!$omp threadprivate(sigma2_on_mu2_ip_array_cloud, sigma2_on_mu2_ip_array_below, &
!$omp   corr_array_cloud, corr_array_below)

  real( kind = core_rknd ), public, dimension(:,:), allocatable :: &
      corr_array_cloud_def, &
      corr_array_below_def
!$omp threadprivate( corr_array_cloud_def, corr_array_below_def )


  private

  public :: read_correlation_matrix, setup_pdf_indices, setup_corr_varnce_array, &
            cleanup_corr_matrix_arrays, hm_idx, init_clubb_arrays

  private :: get_corr_var_index, return_pdf_index, def_corr_idx


  contains

  !-----------------------------------------------------------------------------
  subroutine init_default_corr_arrays(  ) 

    ! Description:
    ! Initializes the default correlation arrays.
    !---------------------------------------------------------------------------

    use constants_clubb, only: &
        one,  & ! Constant(s)
        zero

    implicit none

    integer:: indx

    ! This "renaming" is used to shorten the matrix declarations below.
    integer, parameter :: c = core_rknd

    ! ---- Begin Code ----
 
    ! Allocate Arrays.
    allocate( corr_array_cloud_def(d_var_total,d_var_total) )
    allocate( corr_array_below_def(d_var_total,d_var_total) )

    ! Initialize all values to 0.
    corr_array_cloud_def = zero
    corr_array_below_def = zero

    ! Set the correlation of any variable with itself to 1.
    do indx = 1, d_var_total, 1
       corr_array_cloud_def(indx,indx) = one
       corr_array_below_def(indx,indx) = one
    enddo

    ! Set up default correlation arrays.
    ! The default correlation arrays used here are the correlation arrays used
    ! for the ARM 97 case.  Any changes should be made concurrently here and in
    ! ../../input/case_setups/arm_97_corr_array_cloud.in (for "in-cloud") and
    ! in ../../input/case_setups/arm_97_corr_array_cloud.in (for "below-cloud").
    corr_array_cloud_def = reshape( &

(/1._c, -.6_c, .09_c , .09_c , .5_c   , .5_c   , .2_c   , .2_c  , .2_c  , .2_c  , .2_c, .2_c, &! s
  0._c, 1._c , .027_c, .027_c, .0726_c, .0855_c, -.024_c, .084_c, .018_c, .012_c, 0._c, 0._c, &! t
  0._c, 0._c , 1._c  , .34_c , 0.2_c  , 0.2_c  ,  .1_c  , .15_c , 0._c  , 0._c  , 0._c, 0._c, &! w
  0._c, 0._c , 0._c  , 1._c  , 0._c   , 0._c   ,  .39_c , .29_c , .14_c , .21_c , 0._c, 0._c, &! Ncn
  0._c, 0._c , 0._c  , 0._c  , 1._c   , .7_c   ,  0._c  , 0._c  , .1_c  , .1_c  , .2_c, .2_c, &! rr
  0._c, 0._c , 0._c  , 0._c  , 0._c   , 1._c   ,  .1_c  , .1_c  , 0._c  , 0._c  , .2_c, .2_c, &! Nr
  0._c, 0._c , 0._c  , 0._c  , 0._c   , 0._c   ,  1._c  , .7_c  , .5_c  , .5_c  , .3_c, .3_c, &! ri
  0._c, 0._c , 0._c  , 0._c  , 0._c   , 0._c   ,  0._c  , 1._c  , .5_c  , .5_c  , .3_c, .3_c, &! Ni
  0._c, 0._c , 0._c  , 0._c  , 0._c   , 0._c   ,  0._c  , 0._c  , 1._c  , .7_c  , .4_c, .4_c, &! rs
  0._c, 0._c , 0._c  , 0._c  , 0._c   , 0._c   ,  0._c  , 0._c  , 0._c  , 1._c  , .4_c, .4_c, &! Ns
  0._c, 0._c , 0._c  , 0._c  , 0._c   , 0._c   ,  0._c  , 0._c  , 0._c  , 0._c  , 1._c, .7_c, &! rg
  0._c, 0._c , 0._c  , 0._c  , 0._c   , 0._c   ,  0._c  , 0._c  , 0._c  , 0._c  , 0._c, 1._c/),&!Ng

    shape(corr_array_cloud_def) )
!  s     t     w       Ncn     rr       Nr        ri      Ni      rs      Ns      rg    Ng

    corr_array_cloud_def = transpose( corr_array_cloud_def )


    corr_array_below_def = reshape( &

(/1._c, .3_c , .09_c , .09_c , .5_c   , .5_c   , .2_c   , .2_c  , .2_c  , .2_c  , .2_c, .2_c, &! s
  0._c, 1._c , .027_c, .027_c, .0726_c, .0855_c, -.024_c, .084_c, .018_c, .012_c, 0._c, 0._c, &! t
  0._c, 0._c , 1._c  , .34_c , 0.2_c  , 0.2_c  ,  .1_c  , .15_c , 0._c  , 0._c  , 0._c, 0._c, &! w
  0._c, 0._c , 0._c  , 1._c  , 0._c   , 0._c   ,  .39_c , .29_c , .14_c , .21_c , 0._c, 0._c, &! Ncn
  0._c, 0._c , 0._c  , 0._c  , 1._c   , .7_c   ,  0._c  , 0._c  , .1_c  , .1_c  , .2_c, .2_c, &! rr
  0._c, 0._c , 0._c  , 0._c  , 0._c   , 1._c   ,  .1_c  , .1_c  , 0._c  , 0._c  , .2_c, .2_c, &! Nr
  0._c, 0._c , 0._c  , 0._c  , 0._c   , 0._c   ,  1._c  , .7_c  , .5_c  , .5_c  , .3_c, .3_c, &! ri
  0._c, 0._c , 0._c  , 0._c  , 0._c   , 0._c   ,  0._c  , 1._c  , .5_c  , .5_c  , .3_c, .3_c, &! Ni
  0._c, 0._c , 0._c  , 0._c  , 0._c   , 0._c   ,  0._c  , 0._c  , 1._c  , .7_c  , .4_c, .4_c, &! rs
  0._c, 0._c , 0._c  , 0._c  , 0._c   , 0._c   ,  0._c  , 0._c  , 0._c  , 1._c  , .4_c, .4_c, &! Ns
  0._c, 0._c , 0._c  , 0._c  , 0._c   , 0._c   ,  0._c  , 0._c  , 0._c  , 0._c  , 1._c, .7_c, &! rg
  0._c, 0._c , 0._c  , 0._c  , 0._c   , 0._c   ,  0._c  , 0._c  , 0._c  , 0._c  , 0._c, 1._c/),&!Ng

    shape(corr_array_below_def) )
!  s     t     w       Ncn     rr       Nr       ri       Ni      rs      Ns      rg    Ng

    corr_array_below_def = transpose( corr_array_below_def )


    return

  end subroutine init_default_corr_arrays

  !-----------------------------------------------------------------------------
  pure function def_corr_idx( iiPDF_x ) result(ii_def_corr)

    ! Description:
    !   Map from a iiPDF index to the corresponding index in the default 
    !   correlation arrays.
    !-----------------------------------------------------------------------------

    implicit none

    ! Constant Parameters

    ! Indices that represent the order in the default corr arrays
    ! (s, t, w, Ncn, rr, Nr, ri, Ni, rs, Ns, rg, Ng)
    integer, parameter :: &
    ii_s = 1, &
    ii_t = 2, &
    ii_w = 3, &
    ii_Ncn = 4, &
    ii_rr = 5, &
    ii_Nr = 6, &
    ii_ri = 7, &
    ii_Ni = 8, &
    ii_rs = 9, &
    ii_Ns = 10, &
    ii_rg = 11, &
    ii_Ng = 12

    ! Input Variables

    integer, intent(in) :: iiPDF_x

    ! Return Variable

    integer :: ii_def_corr

    ! ---- Begin Code ----

    ii_def_corr = -1

      if (iiPDF_x == iiPDF_s_mellor) then
         ii_def_corr = ii_s

      elseif (iiPDF_x == iiPDF_t_mellor) then
        ii_def_corr = ii_t

      elseif (iiPDF_x == iiPDF_w) then
        ii_def_corr = ii_w

      elseif (iiPDF_x == iiPDF_Ncn) then
        ii_def_corr = ii_Ncn

      elseif (iiPDF_x == iiPDF_rrain) then
        ii_def_corr = ii_rr

      elseif (iiPDF_x == iiPDF_Nr) then
        ii_def_corr = ii_Nr

      elseif (iiPDF_x == iiPDF_rice) then
        ii_def_corr = ii_ri

      elseif (iiPDF_x == iiPDF_Ni) then
        ii_def_corr = ii_Ni

      elseif (iiPDF_x == iiPDF_rsnow) then
        ii_def_corr = ii_rs

      elseif (iiPDF_x == iiPDF_Nsnow) then
        ii_def_corr = ii_Ns

      elseif (iiPDF_x == iiPDF_rgraupel) then
        ii_def_corr = ii_rg

      elseif (iiPDF_x == iiPDF_Ngraupel) then
        ii_def_corr = ii_Ng

      endif
  end function def_corr_idx

  !-----------------------------------------------------------------------------
  subroutine set_corr_arrays_to_default(  ) 

    ! Description:
    !   If there are no corr_array.in files for the current case, default 
    !   correlations are used. 
    !-----------------------------------------------------------------------------
  
    use constants_clubb, only: &
        zero, &
        one

    implicit none

    ! Local Variables
    integer :: i, j ! Loop iterators


    ! ---- Begin Code ----

    corr_array_cloud = zero
    corr_array_below = zero

    do i = 1, d_variables
       corr_array_cloud(i,i) = one
       corr_array_below(i,i) = one
    enddo

    do i = 1, d_variables-1
       do j = i+1, d_variables
          if ( def_corr_idx(i) > def_corr_idx(j) ) then
             corr_array_cloud(j, i) = corr_array_cloud_def(def_corr_idx(j), def_corr_idx(i))
             corr_array_below(j, i) = corr_array_below_def(def_corr_idx(j), def_corr_idx(i))
          else
             corr_array_cloud(j, i) = corr_array_cloud_def(def_corr_idx(i), def_corr_idx(j))
             corr_array_below(j, i) = corr_array_below_def(def_corr_idx(i), def_corr_idx(j))
          endif
       enddo
    enddo

  end subroutine set_corr_arrays_to_default


  !-----------------------------------------------------------------------------
  subroutine read_correlation_matrix( iunit, input_file, d_variables, &
                                      corr_array )

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
      d_variables ! number of variables in the array

    character(len=*), intent(in) :: input_file ! Path to the file

    ! Input/Output Variable(s)
    real( kind = core_rknd ), dimension(d_variables,d_variables), intent(inout) :: &
      corr_array ! Correlation variance array

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

    ! Allocate all arrays based on d_variables
    allocate( retVars(1:nCols) )

    ! Initializing to zero means that correlations we don't have
    ! (e.g. Nc and any variable other than s_mellor ) are assumed to be 0.
    corr_array(:,:) = 0.0_core_rknd

    ! Set main diagonal to 1
    do i=1, d_variables
      corr_array(i,i) = 1.0_core_rknd
    end do

    ! Read the values from the specified file
    call read_one_dim_file( iunit, nCols, input_file, retVars )

    if( size( retVars(1)%values ) /= nCols ) then
      write(fstderr, *) "Correlation matrix must have an equal number of rows and cols in file ", &
            input_file
      stop "Bad data in correlation file."
    end if

    ! Start at 2 because the first index is always just 1.0 in the first row
    ! and the rest of the rows are ignored
    do i=2, nCols
      var_index1 = get_corr_var_index( retVars(i)%name )
      if( var_index1 > -1 ) then
        do j=1, (i-1)
          var_index2 = get_corr_var_index( retVars(j)%name )
          if( var_index2 > -1 ) then
            call set_lower_triangular_matrix &
                 ( d_variables, var_index1, var_index2, retVars(i)%values(j), &
                   corr_array )
          end if
        end do
      end if
    end do

    call deallocate_one_dim_vars( nCols, retVars )

    return
  end subroutine read_correlation_matrix

  !--------------------------------------------------------------------------
  function get_corr_var_index( var_name ) result( i )

    ! Definition:
    !   Returns the index for a variable based on its name.
    !--------------------------------------------------------------------------

    implicit none

    character(len=*), intent(in) :: var_name ! The name of the variable

    ! Output variable
    integer :: i

    !------------------ BEGIN CODE -----------------------------
    i = -1

    select case( trim(var_name) )

    case( "s" )
      i = iiPDF_s_mellor

    case( "t" )
      i = iiPDF_t_mellor

    case( "w" )
      i = iiPDF_w

    case( "Ncn" )
      i = iiPDF_Ncn

    case( "rrain" )
      i = iiPDF_rrain

    case( "Nr" )
      i = iiPDF_Nr

    case( "rice" )
      i = iiPDF_rice

    case( "Ni" )
      i = iiPDF_Ni

    case( "rsnow" )
      i = iiPDF_rsnow

    case( "Nsnow" )
      i = iiPDF_Nsnow
        
    case( "rgraupel" )
      i = iiPDF_rgraupel

    case( "Ngraupel" )
      i = iiPDF_Ngraupel

    end select

    return

  end function get_corr_var_index

  !-----------------------------------------------------------------------
  subroutine setup_pdf_indices( hydromet_dim, iirrainm, iiNrm, &
                                iiricem, iiNim, iirsnowm, iiNsnowm, &
                                iirgraupelm, iiNgraupelm, &
                                l_ice_micro, l_graupel )

  ! Description:
  !
  !   Setup for the iiPDF indices. These indices are used to address s, t, w
  !   and the hydrometeors in the mean/stdev/corr arrays
  !
  ! References:
  !
  !-----------------------------------------------------------------------

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      hydromet_dim    ! Total number of hydrometeor species.

    integer, intent(in) :: &
      iirrainm, & ! Index of rain water mixing ratio
      iiNrm,    & ! Index of rain drop concentration
      iiricem,  & ! Index of ice mixing ratio
      iiNim,    & ! Index of ice crystal concentration
      iirsnowm, & ! Index of snow mixing ratio
      iiNsnowm, & ! Index of snow concentration
      iirgraupelm, & ! Index of graupel mixing ratio
      iiNgraupelm    ! Index of graupel number concentration

    logical, intent(in) :: &
      l_ice_micro, &  ! Whether the microphysics scheme will do ice
      l_graupel       ! True if graupel is used

    ! Local Variables
    integer :: &
      pdf_count, & ! Count number of PDF variables
      i            ! Hydrometeor loop index

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    iiPDF_s_mellor = 1 ! Extended liquid water mixing ratio
    iiPDF_t_mellor = 2 ! 't' orthogonal to 's'
    iiPDF_w        = 3 ! vertical velocity
    iiPDF_Ncn      = 4 ! Cloud nuclei concentration or extended Nc.

    pdf_count = iiPDF_Ncn

    ! Loop over hydrometeors.
    ! Hydrometeor indices in the PDF arrays should be in the same order as
    ! found in the hydrometeor arrays.
    if ( hydromet_dim > 0 ) then

       do i = 1, hydromet_dim, 1

          if ( i == iirrainm ) then
             pdf_count = pdf_count + 1
             iiPDF_rrain = pdf_count
          endif

          if ( i == iiNrm ) then
             pdf_count = pdf_count + 1
             iiPDF_Nr = pdf_count
          endif

          if ( l_ice_micro ) then

             if ( i == iiricem ) then
                pdf_count = pdf_count + 1
                iiPDF_rice = pdf_count
             endif

             if ( i == iiNim ) then
                pdf_count = pdf_count + 1
                iiPDF_Ni = pdf_count
             endif

             if ( i == iirsnowm ) then
                pdf_count = pdf_count + 1
                iiPDF_rsnow = pdf_count
             endif

             if ( i == iiNsnowm ) then
                pdf_count = pdf_count + 1
                iiPDF_Nsnow = pdf_count
             endif

             if ( l_graupel ) then
                if ( i == iirgraupelm ) then
                   pdf_count = pdf_count + 1
                   iiPDF_rgraupel = pdf_count
                endif
        
                if ( i == iiNgraupelm ) then
                   pdf_count = pdf_count + 1
                   iiPDF_Ngraupel = pdf_count
                endif   
             
             else
                
                iiPDF_rgraupel = -1
                iiPDF_Ngraupel = -1
             
             endif

          else

             iiPDF_rice = -1
             iiPDF_Ni = -1
             iiPDF_rsnow = -1
             iiPDF_Nsnow = -1

          endif ! l_ice_micro

       enddo ! i = 1, hydromet_dim, 1

    endif ! hydromet_dim > 0


    d_variables = pdf_count


    return

  end subroutine setup_pdf_indices
  !-----------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  subroutine return_pdf_index( hydromet_index, pdf_count, pdf_index )

  ! Description:
  !   Set the Latin hypercube variable index if the hydrometeor exists
  ! References:
  !   None
  !-------------------------------------------------------------------------

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      hydromet_index

    ! Input/Output Variables
    integer, intent(inout) :: &
      pdf_count

    ! Output Variables
    integer, intent(out) :: &
      pdf_index

    ! ---- Begin Code ----

    if ( hydromet_index > 0 ) then
      pdf_count = pdf_count + 1
      pdf_index = pdf_count
    else
      pdf_index = -1
    end if

    return
  end subroutine return_pdf_index

!===============================================================================
  subroutine setup_corr_varnce_array( input_file_cloud, input_file_below, &
                                      iunit )

! Description:
!   Setup an array with the x'^2/xm^2 variables on the diagonal and the other
!   elements to be correlations between various variables.

! References:
!   None.
!-------------------------------------------------------------------------------

    use parameters_microphys, only: &
      rr_sigma2_on_mu2_ip_cloud, & ! Variables
      Nr_sigma2_on_mu2_ip_cloud, &
      rr_sigma2_on_mu2_ip_below, &
      Nr_sigma2_on_mu2_ip_below, &
      Ncnp2_on_Ncnm2

    use parameters_microphys, only: &
      rs_sigma2_on_mu2_ip_cloud, & ! Variables
      Ns_sigma2_on_mu2_ip_cloud, &
      ri_sigma2_on_mu2_ip_cloud, &
      Ni_sigma2_on_mu2_ip_cloud, &
      rg_sigma2_on_mu2_ip_cloud, &
      Ng_sigma2_on_mu2_ip_cloud, &
      rs_sigma2_on_mu2_ip_below, &
      Ns_sigma2_on_mu2_ip_below, &
      ri_sigma2_on_mu2_ip_below, &
      Ni_sigma2_on_mu2_ip_below, &
      rg_sigma2_on_mu2_ip_below, &
      Ng_sigma2_on_mu2_ip_below

    use parameters_microphys, only: &
      l_fix_s_t_correlations ! Variable(s)

!   use matrix_operations, only: print_lower_triangular_matrix ! Procedure(s)

    use constants_clubb, only: &
      fstderr, &  ! Constant(s)
      zero

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: max, epsilon, trim

    ! Input Variables
    integer, intent(in) :: &
      iunit ! The file unit

    character(len=*), intent(in) :: &
      input_file_cloud, & ! Path to the in cloud correlation file
      input_file_below    ! Path to the out of cloud correlation file

    ! Local variables
    character(len=1) :: response
    logical :: l_warning, corr_file_exist
    integer :: i

    ! ---- Begin Code ----

    allocate( corr_array_cloud(d_variables,d_variables) )
    allocate( corr_array_below(d_variables,d_variables) )

    allocate( sigma2_on_mu2_ip_array_cloud(d_variables) )
    allocate( sigma2_on_mu2_ip_array_below(d_variables) )

    sigma2_on_mu2_ip_array_cloud(:) = zero
    sigma2_on_mu2_ip_array_below(:) = zero

    ! corr_file_exist is true if the *_corr_array_cloud.in file exists
    ! Note: It is assumed that if the *_corr_array_cloud.in file exists
    !       then *_corr_array_below.in also exists
    inquire( file = input_file_cloud, exist = corr_file_exist )

    if ( corr_file_exist ) then

       call read_correlation_matrix( iunit, trim( input_file_cloud ), d_variables, & ! In
                                     corr_array_cloud ) ! Out

       call read_correlation_matrix( iunit, trim( input_file_below ), d_variables, & ! In
                                     corr_array_below ) ! Out

    else ! Read in default correlation matrices

       write(fstderr,*) "Warning: "//trim( input_file_cloud )//" was not found! " // &
                        "The default correlation arrays will be used."

       call init_default_corr_arrays( )

       call set_corr_arrays_to_default( )

    endif

    ! Sanity check to avoid confusing non-convergence results.
    if ( .not. l_fix_s_t_correlations .and. iiPDF_Ncn > 0 ) then
      l_warning = .false.
      do i = 1, d_variables
        if ( ( corr_array_cloud(i,iiPDF_Ncn) /= zero .or.  &
               corr_array_below(i,iiPDF_Ncn) /= zero ) .and. &
             i /= iiPDF_Ncn ) then
          l_warning = .true.
        end if
      end do ! 1..d_variables
      if ( l_warning ) then
        write(fstderr,*) "Warning: the specified correlations for s and Nc are non-zero."
        write(fstderr,*) "The latin hypercube code will not converge to the analytic solution "// &
          "using these settings."
        write(fstderr,'(A)',advance='no') "Continue? "
        read(*,*) response
        if ( response(1:1) /= 'y' .and. response(1:1) /= 'Y' ) then
           stop "Exiting..."
        end if
      end if
    end if ! l_fix_s_t_correlations

    if ( iiPDF_Ncn > 0 ) then
      sigma2_on_mu2_ip_array_cloud(iiPDF_Ncn) = Ncnp2_on_Ncnm2
    end if

    if ( iiPDF_rrain > 0 ) then
      sigma2_on_mu2_ip_array_cloud(iiPDF_rrain) = rr_sigma2_on_mu2_ip_cloud
      if ( iiPDF_Nr > 0 ) then
        sigma2_on_mu2_ip_array_cloud(iiPDF_Nr) = Nr_sigma2_on_mu2_ip_cloud
      end if ! iiPDF_Nr > 0
    end if ! iiPDF_rrain > 0

    if ( iiPDF_rsnow > 0 ) then
      sigma2_on_mu2_ip_array_cloud(iiPDF_rsnow) = rs_sigma2_on_mu2_ip_cloud


      if ( iiPDF_Nsnow > 0 ) then
        sigma2_on_mu2_ip_array_cloud(iiPDF_Nsnow) = Ns_sigma2_on_mu2_ip_cloud


      end if ! iiPDF_Nsnow > 0
    end if ! iiPDF_rsnow > 0

    if ( iiPDF_rice > 0 ) then
      sigma2_on_mu2_ip_array_cloud(iiPDF_rice) = ri_sigma2_on_mu2_ip_cloud


      if ( iiPDF_Ni > 0 ) then
        sigma2_on_mu2_ip_array_cloud(iiPDF_Ni) = Ni_sigma2_on_mu2_ip_cloud

      end if ! iiPDF_Ni > 0
    end if ! iiPDF_rice > 0

    ! Sampling for graupel (disabled)
    if ( iiPDF_rgraupel > 0 ) then
      sigma2_on_mu2_ip_array_cloud(iiPDF_rgraupel) = rg_sigma2_on_mu2_ip_cloud


      if ( iiPDF_Ngraupel > 0 ) then
        sigma2_on_mu2_ip_array_cloud(iiPDF_Ngraupel) = Ng_sigma2_on_mu2_ip_cloud


      end if ! iiPDF_Ngraupel > 0
    end if ! iiPDF_rgraupel > 0

    if ( iiPDF_Ncn > 0 ) then
      sigma2_on_mu2_ip_array_below(iiPDF_Ncn) = Ncnp2_on_Ncnm2
    end if

    if ( iiPDF_rrain > 0 ) then
      sigma2_on_mu2_ip_array_below(iiPDF_rrain) = rr_sigma2_on_mu2_ip_below



      if ( iiPDF_Nr > 0 ) then
        sigma2_on_mu2_ip_array_below(iiPDF_Nr) = Nr_sigma2_on_mu2_ip_below


      end if ! iiPDF_Nr > 0
    end if ! iiPDF_rrain > 0

    if ( iiPDF_rsnow > 0 ) then
      sigma2_on_mu2_ip_array_below(iiPDF_rsnow) = rs_sigma2_on_mu2_ip_below


      if ( iiPDF_Nsnow > 0 ) then
        sigma2_on_mu2_ip_array_below(iiPDF_Nsnow) = Ns_sigma2_on_mu2_ip_below

      end if ! iiPDF_Nsnow > 0
    end if ! iiPDF_rsnow > 0

    if ( iiPDF_rice > 0 ) then
      sigma2_on_mu2_ip_array_below(iiPDF_rice) = ri_sigma2_on_mu2_ip_below


      if ( iiPDF_Ni > 0 ) then
        sigma2_on_mu2_ip_array_below(iiPDF_Ni) =  Ni_sigma2_on_mu2_ip_below
      end if ! iiPDF_Ni > 0

    end if ! iiPDF_rice > 0

    if ( iiPDF_rgraupel > 0 ) then
      sigma2_on_mu2_ip_array_below(iiPDF_rgraupel) = rg_sigma2_on_mu2_ip_below


      if ( iiPDF_Ngraupel > 0 ) then
        sigma2_on_mu2_ip_array_below(iiPDF_Ngraupel) = Ng_sigma2_on_mu2_ip_below


      end if ! iiPDF_Ngraupel > 0
    end if ! iiPDF_rgraupel > 0

    return
  end subroutine setup_corr_varnce_array

  !-----------------------------------------------------------------------------
  subroutine cleanup_corr_matrix_arrays( )

    ! Description:
    !   De-allocate latin hypercube arrays
    ! References:
    !   None
    !---------------------------------------------------------------------------
    implicit none

    ! External
    intrinsic :: allocated

    ! ---- Begin Code ----

    if ( allocated( corr_array_cloud ) ) then
      deallocate( corr_array_cloud )
    end if

    if ( allocated( corr_array_below ) ) then
      deallocate( corr_array_below )
    end if

    if ( allocated( sigma2_on_mu2_ip_array_cloud ) ) then
      deallocate( sigma2_on_mu2_ip_array_cloud )
    end if

    if ( allocated( sigma2_on_mu2_ip_array_below ) ) then
      deallocate( sigma2_on_mu2_ip_array_below )
    end if

    if ( allocated( corr_array_cloud_def ) ) then
      deallocate( corr_array_cloud_def )
    end if

    if ( allocated( corr_array_below_def ) ) then
      deallocate( corr_array_below_def )
    end if

    return
  end subroutine cleanup_corr_matrix_arrays

 !-----------------------------------------------------------------------------
  subroutine init_clubb_arrays( hydromet_dim, iirrainm, iiNrm, iirsnowm, & ! Variables
                                iiricem, iiNsnowm, iiNim, &
                                iirgraupelm, iiNgraupelm, &
                                l_ice_micro, l_graupel, iunit )

    ! Description: This subroutine sets up arrays that are necessary for WRF.
    !   
    !-----------------------------------------------------------------------------

    implicit none

    ! Constant Parameter(s)
    character(len=*), parameter :: &
      lh_file_path_below = "./clubb_corr_array_below.in", &
      lh_file_path_cloud = "./clubb_corr_array_cloud.in"    

    ! Input Variables
    integer, intent(in) :: &
      hydromet_dim, &
      iirrainm, &
      iiNrm, &
      iirsnowm, & 
      iiricem, & 
      iiNsnowm, &
      iiNim, &
      iirgraupelm, &
      iiNgraupelm, &
      iunit

    logical, intent(in) :: l_ice_micro, l_graupel

    ! ---- Begin Code ----

    call setup_pdf_indices( hydromet_dim, iirrainm, iiNrm, &
                            iiricem, iiNim, iirsnowm, iiNsnowm, &
                            iirgraupelm, iiNgraupelm, &
                            l_ice_micro, l_graupel )

    ! Setup the arrays and indices containing the correlations, etc.
    call setup_corr_varnce_array( lh_file_path_cloud, lh_file_path_below, iunit )

    return    

  end subroutine init_clubb_arrays

  !-----------------------------------------------------------------------
  function hm_idx(iiPDF_idx) result(ii_idx)
  ! Description:
  ! Returns the position of a certain hydrometeor within the hydromet arrays
  ! according to its iiPDF index.

  ! References:
  !
  !-----------------------------------------------------------------------

    use array_index, only: &
        iirrainm, &
        iiNrm, &
        iirsnowm, &
        iiNsnowm, &
        iiricem, &
        iiNim, &
        iirgraupelm, &
        iiNgraupelm

      implicit none

    ! Input Variables
    integer, intent(in) :: iiPDF_idx

    ! Return Variable
    integer :: ii_idx

    ! Local Variables

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    ! Get rid of an annoying compiler warning.
    ii_idx = 1

    if ( iiPDF_idx == iiPDF_rrain ) then
       ii_idx = iirrainm
    endif

    if ( iiPDF_idx == iiPDF_Nr ) then
       ii_idx = iiNrm
    endif

    if ( iiPDF_idx == iiPDF_rsnow ) then
       ii_idx = iirsnowm
    endif

    if ( iiPDF_idx == iiPDF_Nsnow ) then
       ii_idx = iiNsnowm
    endif

    if ( iiPDF_idx == iiPDF_rice ) then
       ii_idx = iiricem
    endif

    if ( iiPDF_idx == iiPDF_Ni ) then
       ii_idx = iiNim
    endif

    if ( iiPDF_idx == iiPDF_rgraupel ) then
       ii_idx = iirgraupelm
    endif

    if ( iiPDF_idx == iiPDF_Ngraupel ) then
       ii_idx = iiNgraupelm
    endif

    return

  end function hm_idx
  !-----------------------------------------------------------------------

end module corr_matrix_module
