!$Id$
!-------------------------------------------------------------------------------
module corr_matrix_module

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

  public :: read_correlation_matrix, setup_pdf_indices

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
  subroutine setup_pdf_indices( iirrainm, iiNrm, iiricem, iiNim, iirsnowm, iiNsnowm, &
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
      iirrainm, & ! Index of rain water mixing ratio
      iiNrm,    & ! Index of rain droplet number conc.
      iiricem,  & ! Index of ice water mixing ratio
      iiNim,    & ! Index of ice crystal number conc.
      iirsnowm, & ! Index snow mixing ratio
      iiNsnowm    ! Index of snow number conc.

    logical, intent(in) :: &
      l_ice_micro  ! Whether the microphysics scheme will do ice

    ! Local Variables
    integer :: i

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    iiPDF_s_mellor = 1 ! Extended rcm
    iiPDF_t_mellor = 2 ! 't' orthogonal to 's'
    iiPDF_w        = 3 ! vertical velocity
    iiPDF_Ncn      = 4 ! Cloud droplet number concentration

    i = iiPDF_Ncn

    call return_pdf_index( iirrainm, i, iiPDF_rrain )
    call return_pdf_index( iiNrm, i, iiPDF_Nr )
    if ( l_ice_micro ) then
      call return_pdf_index( iiricem, i, iiPDF_rice )
      call return_pdf_index( iiNim, i, iiPDF_Ni )
      call return_pdf_index( iirsnowm, i, iiPDF_rsnow )
      call return_pdf_index( iiNsnowm, i, iiPDF_Nsnow )
    else
      iiPDF_rice = -1
      iiPDF_Ni = -1
      iiPDF_rsnow = -1
      iiPDF_Nsnow = -1
    end if
    ! Disabled until we have values for the correlations of graupel and
    ! other variates in the latin hypercube sampling.
    iiPDF_rgraupel = -1
    iiPDF_Ngraupel = -1

    d_variables = i

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

end module corr_matrix_module
