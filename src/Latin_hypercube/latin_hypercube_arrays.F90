! $Id$
module latin_hypercube_arrays

  implicit none

  public :: setup_corr_varnce_array, cleanup_latin_hypercube_arrays
  private :: return_LH_index, read_correlation_matrix
  private

  integer, public :: d_variables
!omp threadprivate(d_variables)

  real, public, dimension(:), allocatable :: &
    xp2_on_xm2_array_cloud, &
    xp2_on_xm2_array_below

  real, public, dimension(:,:), allocatable :: &
    corr_array_cloud, &
    corr_array_below
!$omp threadprivate(xp2_on_xm2_array_cloud, xp2_on_xm2_array_below, &
!$omp   corr_array_cloud, corr_array_below)

  logical, public :: &
    l_fixed_corr_initialized = .false.

!$omp threadprivate(l_fixed_corr_initialized)

  double precision, allocatable, dimension(:,:), target, public :: &
    corr_stw_cloud_Cholesky, & ! Cholesky factorization of the correlation matrix
    corr_stw_below_Cholesky    ! Cholesky factorization of the correlation matrix

!$omp threadprivate(corr_stw_cloud_Cholesky, corr_stw_below_Cholesky)

  double precision, allocatable, dimension(:), public :: &
    corr_stw_cloud_scaling, & ! Scaling factors for the correlation matrix [-]
    corr_stw_below_scaling    ! Scaling factors for the correlation matrix [-]

!$omp threadprivate(corr_stw_cloud_scaling, corr_stw_below_scaling)

  logical, public :: &
    l_corr_stw_cloud_scaling, & ! Whether we're scaling the correlation matrix
    l_corr_stw_below_scaling

!$omp threadprivate(l_corr_stw_cloud_scaling, l_corr_stw_below_scaling)

  integer, allocatable, dimension(:,:,:), public :: & 
    height_time_matrix ! matrix of rand ints

!$omp threadprivate(height_time_matrix)

  ! Latin hypercube indices
  integer, public :: &
    iiLH_s_mellor, iiLH_t_mellor, iiLH_w
!$omp threadprivate(iiLH_s_mellor, iiLH_t_mellor, iiLH_w)

  integer, public :: &
   iiLH_rrain, iiLH_rsnow, iiLH_rice, iiLH_rgraupel
!$omp threadprivate(iiLH_rrain, iiLH_rsnow, iiLH_rice, iiLH_rgraupel)

  integer, public :: &
   iiLH_Nr, iiLH_Nsnow, iiLH_Ni, iiLH_Ngraupel, iiLH_Nc
!$omp threadprivate(iiLH_Nr, iiLH_Nsnow, iiLH_Ni, iiLH_Ngraupel, iiLH_Nc)

  contains
!===============================================================================
  subroutine setup_corr_varnce_array( iiNcm, iirrainm, iiNrm, iiricem, iiNim, iirsnowm, iiNsnowm, &
                                      input_path, iunit, runtype )

! Description:
!   Setup an array with the x'^2/xm^2 variables on the diagonal and the other
!   elements to be correlations between various variables.

! References:
!   None.
!-------------------------------------------------------------------------------

    use parameters_microphys, only: &
      rrp2_on_rrainm2_cloud, Nrp2_on_Nrm2_cloud, Ncp2_on_Ncm2_cloud, & ! Variables
      rrp2_on_rrainm2_below, Nrp2_on_Nrm2_below, Ncp2_on_Ncm2_below

    use parameters_microphys, only: &
      rsnowp2_on_rsnowm2_cloud, & ! Variables
      Nsnowp2_on_Nsnowm2_cloud, & 
      ricep2_on_ricem2_cloud, & 
      Nicep2_on_Nicem2_cloud, &
      rsnowp2_on_rsnowm2_below, & 
      Nsnowp2_on_Nsnowm2_below, & 
      ricep2_on_ricem2_below, & 
      Nicep2_on_Nicem2_below
 
    use matrix_operations, only: set_lower_triangular_matrix_sp ! Procedure(s)

    use constants_clubb, only: fstdout

    implicit none

    ! External
    intrinsic :: max, epsilon

    ! Constant Parameters

    character(len=*), parameter :: &
      cloud_file_name = "_corr_array_cloud.in", & ! File names
      below_file_name = "_corr_array_below.in"
 
    ! Input Variables
    integer, intent(in) :: &
      iiNcm,    & ! Index of cloud droplet number conc.
      iirrainm, & ! Index of rain water mixing ratio
      iiNrm,    & ! Index of rain droplet number conc.
      iiricem,  & ! Index of ice water mixing ratio
      iiNim,    & ! Index of ice crystal number conc.
      iirsnowm, & ! Index snow mixing ratio
      iiNsnowm    ! Index of snow number conc.

    integer, intent(in) :: &
      iunit ! The file unit

    character(len=*), intent(in) :: &
      input_path, & ! Path to the correlation files
      runtype          ! The type of this run.

    ! Local variables

    integer :: i

    character(len=50) :: file_path

    ! ---- Begin Code ----
    iiLH_s_mellor = 1 ! Extended rcm
    iiLH_t_mellor = 2 ! 't' orthogonal to 's'
    iiLH_w        = 3 ! vertical velocity

    i = iiLH_w

    call return_LH_index( iiNcm, i, iiLH_Nc )
    call return_LH_index( iirrainm, i, iiLH_rrain )
    call return_LH_index( iiNrm, i, iiLH_Nr )
    call return_LH_index( iiricem, i, iiLH_rice )
    call return_LH_index( iiNim, i, iiLH_Ni )
    call return_LH_index( iirsnowm, i, iiLH_rsnow )
    call return_LH_index( iiNsnowm, i, iiLH_Nsnow )
    ! Disabled until we have values for the correlations of graupel and
    ! other variates in the latin hypercube sampling.
    iiLH_rgraupel = -1
    iiLH_Ngraupel = -1

    d_variables = i

    allocate( corr_array_cloud(d_variables,d_variables) )
    allocate( corr_array_below(d_variables,d_variables) )

    allocate( xp2_on_xm2_array_cloud(d_variables) )
    allocate( xp2_on_xm2_array_below(d_variables) )

    xp2_on_xm2_array_cloud(:) = 0.0
    xp2_on_xm2_array_below(:) = 0.0

    file_path = input_path//trim( runtype )//cloud_file_name

    call read_correlation_matrix( iunit, file_path, d_variables, &
                                  corr_array_cloud )

    file_path = input_path//trim( runtype )//below_file_name

    call read_correlation_matrix( iunit, file_path, d_variables, &
                                  corr_array_below )

    if ( iiLH_Nc > 0 ) then
      xp2_on_xm2_array_cloud(iiLH_Nc) = Ncp2_on_Ncm2_cloud
    end if

    if ( iiLH_rrain > 0 ) then
      xp2_on_xm2_array_cloud(iiLH_rrain) = rrp2_on_rrainm2_cloud
      if ( iiLH_Nr > 0 ) then
        xp2_on_xm2_array_cloud(iiLH_Nr) = Nrp2_on_Nrm2_cloud
      end if ! iiLH_Nr > 0
    end if ! iiLH_rrain > 0

    if ( iiLH_rsnow > 0 ) then
      xp2_on_xm2_array_cloud(iiLH_rsnow) = rsnowp2_on_rsnowm2_cloud


      if ( iiLH_Nsnow > 0 ) then
        xp2_on_xm2_array_cloud(iiLH_Nsnow) = Nsnowp2_on_Nsnowm2_cloud


      end if ! iiLH_Nsnow > 0
    end if ! iiLH_rsnow > 0

    if ( iiLH_rice > 0 ) then
      xp2_on_xm2_array_cloud(iiLH_rice) = ricep2_on_ricem2_cloud


      if ( iiLH_Ni > 0 ) then
        xp2_on_xm2_array_cloud(iiLH_Ni) = Nicep2_on_Nicem2_cloud

      end if ! iiLH_Ni > 0
    end if ! iiLH_rice > 0

    ! Sampling for graupel (disabled)
    if ( iiLH_rgraupel > 0 ) then
      xp2_on_xm2_array_cloud(iiLH_rgraupel) = -999.


      if ( iiLH_Ngraupel > 0 ) then
        xp2_on_xm2_array_cloud(iiLH_Ngraupel) = -999.


      end if ! iiLH_Ngraupel > 0
    end if ! iiLH_rgraupel > 0

    if ( iiLH_Nc > 0 ) then
      ! The epsilon is a kluge to prevent a singular matrix in generate_lh_sample
      xp2_on_xm2_array_below(iiLH_Nc) = &
        max( Ncp2_on_Ncm2_below, epsilon( Ncp2_on_Ncm2_below ) )

    end if

    if ( iiLH_rrain > 0 ) then
      xp2_on_xm2_array_below(iiLH_rrain) = rrp2_on_rrainm2_below



      if ( iiLH_Nr > 0 ) then
        xp2_on_xm2_array_below(iiLH_Nr) = Nrp2_on_Nrm2_below


      end if ! iiLH_Nr > 0
    end if ! iiLH_rrain > 0

    if ( iiLH_rsnow > 0 ) then
      xp2_on_xm2_array_below(iiLH_rsnow) = rsnowp2_on_rsnowm2_below


      if ( iiLH_Nsnow > 0 ) then
        xp2_on_xm2_array_below(iiLH_Nsnow) = Nsnowp2_on_Nsnowm2_below

      end if ! iiLH_Nsnow > 0
    end if ! iiLH_rsnow > 0

    if ( iiLH_rice > 0 ) then
      xp2_on_xm2_array_below(iiLH_rice) = ricep2_on_ricem2_below


      if ( iiLH_Ni > 0 ) then
        xp2_on_xm2_array_below(iiLH_Ni) =  Nicep2_on_Nicem2_below
      end if ! iiLH_Ni > 0

    end if ! iiLH_rice > 0

    if ( iiLH_rgraupel > 0 ) then
      xp2_on_xm2_array_below(iiLH_rgraupel) = -999.


      if ( iiLH_Ngraupel > 0 ) then
        xp2_on_xm2_array_below(iiLH_Ngraupel) = -999.


      end if ! iiLH_Ngraupel > 0
    end if ! iiLH_rgraupel > 0

    ! Assume rr and Nr are uncorrelated with w for now.
    if ( iiLH_rrain > 0 .and. iiLH_Nr > 0 ) then
      corr_array_cloud(iiLH_rrain,iiLH_w) = 0.0
      corr_array_cloud(iiLH_Nr,iiLH_w) = 0.0
    end if

    ! Fill in values for t_mellor using s_mellor correlations
    do i = 3, d_variables
      corr_array_cloud(i,iiLH_t_mellor) = corr_array_cloud(i,iiLH_s_mellor) &
        * corr_array_cloud(iiLH_t_mellor,iiLH_s_mellor)
      corr_array_below(i,iiLH_t_mellor) = corr_array_below(i,iiLH_s_mellor) &
        * corr_array_below(iiLH_t_mellor,iiLH_s_mellor)
    end do

    return
  end subroutine setup_corr_varnce_array

  !-----------------------------------------------------------------------------
  subroutine cleanup_latin_hypercube_arrays( )

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

    if ( allocated( corr_stw_cloud_Cholesky ) ) then
      deallocate( corr_stw_cloud_Cholesky )
    end if

    if ( allocated( corr_stw_below_Cholesky ) ) then
      deallocate( corr_stw_below_Cholesky )
    end if

    if ( allocated( corr_stw_cloud_scaling ) ) then
      deallocate( corr_stw_cloud_scaling )
    end if

    if ( allocated( corr_stw_below_scaling ) ) then
      deallocate( corr_stw_below_scaling )
    end if

    if ( allocated( height_time_matrix ) ) then
      deallocate( height_time_matrix )
    end if

    return
  end subroutine cleanup_latin_hypercube_arrays

!-------------------------------------------------------------------------------
  subroutine return_LH_index( hydromet_index, LH_count, LH_index )

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
      LH_count

    ! Output Variables
    integer, intent(out) :: &
      LH_index

    ! ---- Begin Code ----

    if ( hydromet_index > 0 ) then
      LH_count = LH_count + 1
      LH_index = LH_count
    else
      LH_index = -1
    end if

    return
  end subroutine return_LH_index

  !-----------------------------------------------------------------------------
  subroutine read_correlation_matrix( iunit, input_file, d_variables, &
                                     corr_array )

  ! Description:
  !   Reads a correlation variance array from a file and stores it in an array.
  !-----------------------------------------------------------------------------

    use input_reader, only: &
      one_dim_read_var, & ! Variable(s)
      read_one_dim_file, deallocate_one_dim_vars, count_columns ! Procedure(s)

    use matrix_operations, only: set_lower_triangular_matrix_sp ! Procedure(s)

    use constants_clubb, only: fstderr ! Variable(s)

    implicit none

    ! Input Variable(s)
    integer, intent(in) :: &
      iunit, &    ! File I/O unit
      d_variables ! number of variables in the array

    character(len=*), intent(in) :: input_file ! Path to the file

    ! Input/Output Variable(s)
    real, dimension(d_variables,d_variables), intent(inout) :: &
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
    corr_array(:,:) = 0.0

    ! Set main diagonal to 1
    do i=1, d_variables
      corr_array(i,i) = 1.0
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
    write(fstderr, *) d_variables
    do i=2, nCols
      var_index1 = get_corr_var_index( retVars(i)%name )
      write(fstderr, *) retVars(i)%name
      if( var_index1 > -1 ) then
        do j=1, (i-1)
          var_index2 = get_corr_var_index( retVars(j)%name )
          if( var_index2 > -1 ) then
            write(fstderr, *) var_index1,", ",var_index2,": ",retVars(i)%values(j)
            call set_lower_triangular_matrix_sp &
                 ( d_variables, var_index1, var_index2, retVars(i)%values(j), &
                   corr_array )
          end if
        end do
      end if
    end do

    call deallocate_one_dim_vars( d_variables, retVars )

    deallocate( retVars )

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

    select case( trim(var_name) )

      case( "s" )
        i = iiLH_s_mellor
 
      case( "t" )
        i = iiLH_t_mellor

      case( "w" )
        i = iiLH_w

      case( "Nc" )
        i = iiLH_Nc

      case( "rrain" )
        i = iiLH_rrain

      case( "Nr" )
        i = iiLH_Nr

      case( "rice" )
        i = iiLH_rice

      case( "Ni" )
        i = iiLH_Ni

      case( "rsnow" )
        i = iiLH_rsnow

      case( "Nsnow" )
        i = iiLH_Nsnow

      case default
        i = -1

    end select

    return

  end function get_corr_var_index

end module latin_hypercube_arrays
