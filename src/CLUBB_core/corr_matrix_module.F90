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

  integer, public :: &
    d_variables
!$omp threadprivate(d_variables)

  real( kind = core_rknd ), public, dimension(:), allocatable :: &
    xp2_on_xm2_array_cloud, &
    xp2_on_xm2_array_below

  real( kind = core_rknd ), public, dimension(:,:), allocatable :: &
    corr_array_cloud, &
    corr_array_below
!$omp threadprivate(xp2_on_xm2_array_cloud, xp2_on_xm2_array_below, &
!$omp   corr_array_cloud, corr_array_below)

  public :: read_correlation_matrix, setup_pdf_indices, setup_corr_varnce_array, &
            cleanup_corr_matrix_arrays

  private :: get_corr_var_index, return_pdf_index

  private

  contains

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

    end select

    return

  end function get_corr_var_index

  !-----------------------------------------------------------------------
  subroutine setup_pdf_indices( hydromet_dim, iirrainm, iiNrm, &
                                iiricem, iiNim, iirsnowm, iiNsnowm, &
                                l_ice_micro )

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
      iiNsnowm    ! Index of snow concentration

    logical, intent(in) :: &
      l_ice_micro  ! Whether the microphysics scheme will do ice

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

          else

             iiPDF_rice = -1
             iiPDF_Ni = -1
             iiPDF_rsnow = -1
             iiPDF_Nsnow = -1

          endif ! l_ice_micro

       enddo ! i = 1, hydromet_dim, 1

    endif ! hydromet_dim > 0

    ! Disabled until we have values for the correlations of graupel and
    ! other variates in the latin hypercube sampling.
    iiPDF_rgraupel = -1
    iiPDF_Ngraupel = -1

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
                                      input_file_cloud_default, input_file_below_default, &
                                      iunit )

! Description:
!   Setup an array with the x'^2/xm^2 variables on the diagonal and the other
!   elements to be correlations between various variables.

! References:
!   None.
!-------------------------------------------------------------------------------

    use parameters_microphys, only: &
      rrp2_on_rrm2_cloud,   & ! Variables
      Nrp2_on_Nrm2_cloud,   &
      Ncp2_on_Ncm2_cloud => Ncnp2_on_Ncnm2_cloud, &
      rrp2_on_rrm2_below,   &
      Nrp2_on_Nrm2_below,   &
      Ncp2_on_Ncm2_below => Ncnp2_on_Ncnm2_below

    use parameters_microphys, only: &
      rsnowp2_on_rsnowm2_cloud, & ! Variables
      Nsnowp2_on_Nsnowm2_cloud, &
      ricep2_on_ricem2_cloud, &
      Nicep2_on_Nicem2_cloud, &
      rsnowp2_on_rsnowm2_below, &
      Nsnowp2_on_Nsnowm2_below, &
      ricep2_on_ricem2_below, &
      Nicep2_on_Nicem2_below

    use parameters_microphys, only: &
      l_fix_s_t_correlations ! Variable(s)

!   use matrix_operations, only: print_lower_triangular_matrix ! Procedure(s)

    use constants_clubb, only: &
      fstdout, & ! Constant(s)
      fstderr, &
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

    character(len=*), intent(in) :: &
      input_file_cloud_default, &
      input_file_below_default

    ! Local variables
    integer :: iiPDF_Nc

    character(len=1) :: response
    logical :: l_warning, corr_file_exist
    integer :: i

    ! ---- Begin Code ----

    iiPDF_Nc = iiPDF_Ncn

    allocate( corr_array_cloud(d_variables,d_variables) )
    allocate( corr_array_below(d_variables,d_variables) )

    allocate( xp2_on_xm2_array_cloud(d_variables) )
    allocate( xp2_on_xm2_array_below(d_variables) )

    xp2_on_xm2_array_cloud(:) = zero
    xp2_on_xm2_array_below(:) = zero

    ! corr_file_exist is true if the *_corr_array_cloud.in file exists
    ! Note: It is assumed that if the *_corr_array_cloud.in file exists
    !       then *_corr_array_below.in also exists
    inquire( file = input_file_cloud, exist = corr_file_exist )

    if ( corr_file_exist ) then

       call read_correlation_matrix( iunit, trim( input_file_cloud ), d_variables, & ! In
                                     corr_array_cloud ) ! Out

       call read_correlation_matrix( iunit, trim( input_file_below ), d_variables, & ! In
                                     corr_array_below ) ! Out

    else ! Read in default correlation matrices (from arm_97_corr_array_cloud)

       write(fstderr,*) "Warning: "//trim( input_file_cloud )//" was not found! " // &
                        "The default correlation arrays (from arm_97) will be used."

       call read_correlation_matrix( iunit, input_file_cloud_default, d_variables, & ! In
                                     corr_array_cloud ) ! In/Out

       call read_correlation_matrix( iunit, input_file_below_default, d_variables, & ! In
                                     corr_array_below ) ! In/Out

    endif

    ! Sanity check to avoid confusing non-convergence results.
    if ( .not. l_fix_s_t_correlations .and. iiPDF_Nc > 0 ) then
      l_warning = .false.
      do i = 1, d_variables
        if ( ( corr_array_cloud(i,iiPDF_Nc) /= zero .or.  &
               corr_array_below(i,iiPDF_Nc) /= zero ) .and. &
             i /= iiPDF_Nc ) then
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

    if ( iiPDF_Nc > 0 ) then
      xp2_on_xm2_array_cloud(iiPDF_Nc) = Ncp2_on_Ncm2_cloud
    end if

    if ( iiPDF_rrain > 0 ) then
      xp2_on_xm2_array_cloud(iiPDF_rrain) = rrp2_on_rrm2_cloud
      if ( iiPDF_Nr > 0 ) then
        xp2_on_xm2_array_cloud(iiPDF_Nr) = Nrp2_on_Nrm2_cloud
      end if ! iiPDF_Nr > 0
    end if ! iiPDF_rrain > 0

    if ( iiPDF_rsnow > 0 ) then
      xp2_on_xm2_array_cloud(iiPDF_rsnow) = rsnowp2_on_rsnowm2_cloud


      if ( iiPDF_Nsnow > 0 ) then
        xp2_on_xm2_array_cloud(iiPDF_Nsnow) = Nsnowp2_on_Nsnowm2_cloud


      end if ! iiPDF_Nsnow > 0
    end if ! iiPDF_rsnow > 0

    if ( iiPDF_rice > 0 ) then
      xp2_on_xm2_array_cloud(iiPDF_rice) = ricep2_on_ricem2_cloud


      if ( iiPDF_Ni > 0 ) then
        xp2_on_xm2_array_cloud(iiPDF_Ni) = Nicep2_on_Nicem2_cloud

      end if ! iiPDF_Ni > 0
    end if ! iiPDF_rice > 0

    ! Sampling for graupel (disabled)
    if ( iiPDF_rgraupel > 0 ) then
      xp2_on_xm2_array_cloud(iiPDF_rgraupel) = -999._core_rknd


      if ( iiPDF_Ngraupel > 0 ) then
        xp2_on_xm2_array_cloud(iiPDF_Ngraupel) = -999._core_rknd


      end if ! iiPDF_Ngraupel > 0
    end if ! iiPDF_rgraupel > 0

    if ( iiPDF_Nc > 0 ) then
      ! The epsilon is a kluge to prevent a singular matrix in generate_lh_sample
      xp2_on_xm2_array_below(iiPDF_Nc) = &
        max( Ncp2_on_Ncm2_below, epsilon( Ncp2_on_Ncm2_below ) )

    end if

    if ( iiPDF_rrain > 0 ) then
      xp2_on_xm2_array_below(iiPDF_rrain) = rrp2_on_rrm2_below



      if ( iiPDF_Nr > 0 ) then
        xp2_on_xm2_array_below(iiPDF_Nr) = Nrp2_on_Nrm2_below


      end if ! iiPDF_Nr > 0
    end if ! iiPDF_rrain > 0

    if ( iiPDF_rsnow > 0 ) then
      xp2_on_xm2_array_below(iiPDF_rsnow) = rsnowp2_on_rsnowm2_below


      if ( iiPDF_Nsnow > 0 ) then
        xp2_on_xm2_array_below(iiPDF_Nsnow) = Nsnowp2_on_Nsnowm2_below

      end if ! iiPDF_Nsnow > 0
    end if ! iiPDF_rsnow > 0

    if ( iiPDF_rice > 0 ) then
      xp2_on_xm2_array_below(iiPDF_rice) = ricep2_on_ricem2_below


      if ( iiPDF_Ni > 0 ) then
        xp2_on_xm2_array_below(iiPDF_Ni) =  Nicep2_on_Nicem2_below
      end if ! iiPDF_Ni > 0

    end if ! iiPDF_rice > 0

    if ( iiPDF_rgraupel > 0 ) then
      xp2_on_xm2_array_below(iiPDF_rgraupel) = -999._core_rknd


      if ( iiPDF_Ngraupel > 0 ) then
        xp2_on_xm2_array_below(iiPDF_Ngraupel) = -999._core_rknd


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

    if ( allocated( xp2_on_xm2_array_cloud ) ) then
      deallocate( xp2_on_xm2_array_cloud )
    end if

    if ( allocated( xp2_on_xm2_array_below ) ) then
      deallocate( xp2_on_xm2_array_below )
    end if

    return
  end subroutine cleanup_corr_matrix_arrays

end module corr_matrix_module
