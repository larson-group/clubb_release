!-----------------------------------------------------------------------
! $Id$
!===============================================================================
module parameters_tunable 

  ! Description:
  ! Contains tuneable model parameters_tunable.

  ! References:
  ! None
  !-----------------------------------------------------------------------

  use parameter_indices, only: nparams ! Variable(s)

  implicit none

  ! Default to private
  private

  public :: setup_parameters, read_parameters, read_param_spread, &
            get_parameters

  ! Model parameters and constraints setup in the namelists
  real, public ::  & 
    T0,       & ! Reference temperature (usually 300)  [K]
    ts_nudge    ! Timescale of u/v nudging             [s]

!$omp   threadprivate(T0, ts_nudge)

  integer, public :: & 
    sclr_dim,        & ! Number of passive scalars
    hydromet_dim       ! Number of hydrometeor species

!$omp   threadprivate(sclr_dim, hydromet_dim)

  real, dimension(:), allocatable, public :: & 
    sclrtol ! Threshold(s) on the passive scalars  [units vary]

!$omp   threadprivate(sclrtol)

  ! Model constant parameters 
  real, public :: & 
     C1,          & ! Low Skewness in C1 Skewness Function.
     C1b,         & ! High Skewness in C1 Skewness Function.
     C1c,         & ! Degree of Slope of C1 Skewness Function.
     C2,          & ! Low Skewness in C2 Skewness Function.
     C2rt,        & ! C2 coefficient for the rtp2_dp1 term.
     C2thl,       & ! C2 coefficient for the thlp2_dp1 term.
     C2rtthl,     & ! C2 coefficient for the rtpthlp_dp1 term.
     C2b,         & ! High Skewness in C2 Skewness Function.  
     C2c,         & ! Degree of Slope of C2 Skewness Function.
     C4,          & ! Used only when l_Kh_zm_aniso is true.
     C5,          & ! Coefficient in pressure terms in the w'^2 eqn.
     C6rt,        & ! Low Skewness in C6rt Skewness Function.
     C6rtb,       & ! High Skewness in C6rt Skewness Function.
     C6rtc,       & ! Degree of Slope of C6rt Skewness Function.
     C6thl,       & ! Low Skewness in C6thl Skewness Function.
     C6thlb,      & ! High Skewness in C6thl Skewness Function.
     C6thlc,      & ! Degree of Slope of C6thl Skewness Function.
     C7,          & ! Low Skewness in C7 Skewness Function.
     C7b,         & ! High Skewness in C7 Skewness Function.
     C7c,         & ! Degree of Slope of C7 Skewness Function.
     C8,          & ! Coefficient #1 in C8 Skewness Equation.
     C8b,         & ! Coefficient #2 in C8 Skewness Equation.  From Golaz.
     C10,         & ! Currently Not Used in the Model.
     C11,         & ! Low Skewness in C11 Skewness Function.
     C11b,        & ! High Skewness in C11 Skewness Function.  From Golaz.
     C11c,        & ! Degree of Slope of C11 Skewness Function.  From Golaz.
     C12,         & ! Constant in w'^3 Crank-Nicholson diffusional term.
     C13,         & ! Not currently used in model.
     C14,         & ! Constant for u'^2 and v'^2 terms.
     c_K,         & ! Constant C_mu^(1/4) in Duynkerke & Driedonks 1987.
     c_K1,        & ! Coefficient of Eddy Diffusion for wp2.
     nu1,         & ! Background Coefficient of Eddy Diffusion for wp2.
     c_K2,        & ! Coefficient of Eddy Diffusion for xp2.
     nu2,         & ! Background Coefficient of Eddy Diffusion for xp2.
     c_K6,        & ! Coefficient of Eddy Diffusion for wpthlp and wprtp.
     nu6,         & ! Background Coefficient of Eddy Diffusion for wpxp.
     c_K8,        & ! Coefficient of Eddy Diffusion for wp3.
     nu8,         & ! Background Coefficient of Eddy Diffusion for wp3.
     c_K9,        & ! Coefficient of Eddy Diffusion for up2 and vp2.
     nu9,         & ! Background Coefficient of Eddy Diffusion for up2 and vp2.
     c_Krrainm,   & ! Coefficient of Eddy Diffusion for hydrometeors.
     nu_r,        & ! Background Coefficient of Eddy Diffusion for hydrometeors.
     c_Ksqd,      & ! Constant for scaling effect of value-squared diffusion.
     nu_hd,       & ! Constant coefficient for 4th-order hyper-diffusion.
     gamma_coef,  & ! Low Skewness in gamma coefficient Skewness Function.
     gamma_coefb, & ! High Skewness in gamma coefficient Skewness Function.
     gamma_coefc, & ! Degree of Slope of gamma coefficient Skewness Function.
     mu,          & ! Fractional entrainment rate per unit altitude.
     taumin,      & ! Previously determined value.
     taumax,      & ! Previously determined value.
     lmin        ! Minimum value for the length scale.

!$omp   threadprivate(C1, C1b, C1c, C2, C2b, C2c)
!$omp   threadprivate(C2rt, C2thl, C2rtthl, C4, C5, C6rt, C6rtb, C6rtc)
!$omp   threadprivate(C6thl, C6thlb, C6thlc)
!$omp   threadprivate(C7, C7b, C7c, C8, C8b, C10, C11, C11b, C11c, C12)
!$omp   threadprivate(C13, C14)
!$omp   threadprivate(c_K, c_K1, nu1, c_K2, nu2, c_K6, nu6)
!$omp   threadprivate(c_K8, nu8, c_K9, nu9, c_Krrainm, nu_r, c_Ksqd, nu_hd)
!$omp   threadprivate(gamma_coef, gamma_coefb, gamma_coefc)
!$omp   threadprivate(taumin, taumax, mu, lmin)

  ! Vince Larson added a constant to set plume widths for theta_l and rt
  ! beta should vary between 0 and 3, with 1.5 the standard value

  real, public :: beta

!$omp   threadprivate(beta)

  real :: lmin_coef ! Coefficient of lmin

!$omp   threadprivate(lmin_coef)

  ! Since we lack a devious way to do this just once, this namelist 
  ! must be changed as well when a new parameter is added.
  namelist /initvars/  & 
     C1, C1b, C1c, C2, C2b, C2c,  & 
     C2rt, C2thl, C2rtthl, C4, C5, & 
     C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc, & 
     C7, C7b, C7c, C8, C8b, C10, C11, C11b, C11c, & 
     C12, C13, C14, c_K, c_K1, nu1, c_K2, nu2,  & 
     c_K6, nu6, c_K8, nu8, c_K9, nu9, c_Krrainm, nu_r, c_Ksqd,  & 
     nu_hd, beta, gamma_coef, gamma_coefb, gamma_coefc, & 
     lmin_coef, taumin, taumax, mu

  ! These are referenced together often enough that it made sense to 
  ! make a list of them.  Note that lmin_coef is the input parameter, 
  ! while the actual lmin model constant is computed from this.
  !***************************************************************
  !                    ***** IMPORTANT *****
  ! If you change the order of the parameters in the parameter_indices,
  ! you will need to change the order of this list as well or the
  ! tuner will break!
  !                    ***** IMPORTANT *****
  !***************************************************************
  character(len=11), dimension(nparams), parameter, public ::  & 
  params_list = & 
     (/"C1         ", "C1b        ", "C1c        ", "C2         ", & 
       "C2b        ", "C2c        ", "C2rt       ", "C2thl      ", & 
       "C2rtthl    ", "C4         ", "C5         ", "C6rt       ", & 
       "C6rtb      ", "C6rtc      ", "C6thl      ", "C6thlb     ", & 
       "C6thlc     ", "C7         ", "C7b        ", "C7c        ", & 
       "C8         ", "C8b        ", "C10        ", "C11        ", & 
       "C11b       ", "C11c       ", "C12        ", "C13        ", & 
       "C14        ", "c_K        ", "c_K1       ", "nu1        ", & 
       "c_K2       ", "nu2        ", "c_K6       ", "nu6        ", & 
       "c_K8       ", "nu8        ", "c_K9       ", "nu9        ", & 
       "c_Krrainm  ", "nu_r       ", "c_Ksqd     ", "nu_hd      ", &
       "gamma_coef ", "gamma_coefb", "gamma_coefc", "mu         ", &
       "beta       ", "lmin_coef  ", "taumin     ", "taumax     " /)

contains

  !=============================================================================
  subroutine setup_parameters & 
            ( deltaz, T0_in, ts_nudge_in, hydromet_dim_in, & 
              sclr_dim_in, sclrtol_in, params, nzmax, l_implemented, &
              grid_type, momentum_heights, thermodynamic_heights, &
              err_code ) 

    ! Description:
    ! Subroutine to setup model parameters

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use constants, only:  & 
        fstderr ! Variable(s)

    use error_code, only:  & 
        clubb_var_out_of_bounds,  & ! Variable(s)
        clubb_no_error

    implicit none

    ! External
    intrinsic :: allocated

    ! Input Variables
    real, intent(in) ::  & 
      deltaz,      & ! Change per height level      [m]
      T0_in,       & ! Ref. temperature             [K]
      ts_nudge_in    ! Timescale for u/v nudging    [s]

    integer, intent(in) :: & 
      hydromet_dim_in,  & ! Number of hydrometeor species
      sclr_dim_in         ! Number of passive scalars

    real, intent(in), dimension(sclr_dim_in) :: & 
      sclrtol_in     ! Threshold on passive scalars

    real, intent(in), dimension(nparams) ::  & 
      params  ! Tuneable model parameters      [-]

    ! Grid definition
    integer, intent(in) :: nzmax  ! Vertical grid levels            [#]

    ! Flag to see if CLUBB is running on it's own,
    ! or if it's implemented as part of a host model.
    logical, intent(in) :: l_implemented   ! (T/F)

    ! If CLUBB is running on it's own, this option determines
    ! if it is using:
    ! 1) an evenly-spaced grid,
    ! 2) a stretched (unevenly-spaced) grid entered on the
    !    thermodynamic grid levels (with momentum levels set
    !    halfway between thermodynamic levels), or
    ! 3) a stretched (unevenly-spaced) grid entered on the
    !    momentum grid levels (with thermodynamic levels set
    !    halfway between momentum levels).
    integer, intent(in) :: grid_type

    ! If the CLUBB parameterization is implemented in a host model,
    ! it needs to use the host model's momentum level altitudes
    ! and thermodynamic level altitudes.
    ! If the CLUBB model is running by itself, but is using a
    ! stretched grid entered on thermodynamic levels (grid_type = 2),
    ! it needs to use the thermodynamic level altitudes as input.
    ! If the CLUBB model is running by itself, but is using a
    ! stretched grid entered on momentum levels (grid_type = 3),
    ! it needs to use the momentum level altitudes as input.
    real, intent(in), dimension(nzmax) :: &
      momentum_heights,      & ! Momentum level altitudes (input)      [m]
      thermodynamic_heights    ! Thermodynamic level altitudes (input) [m]

    ! Output Variables
    integer, intent(out) ::  &
      err_code ! Error condition

    ! Local Variable
    real :: avg_deltaz  ! Average grid box height   [m]


    call unpack_parameters( params, & 
                            C1, C1b, C1c, C2, C2b, C2c, C2rt, C2thl, C2rtthl, &
                            C4, C5, C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc, &
                            C7, C7b, C7c, C8, C8b, C10, & 
                            C11, C11b, C11c, C12, C13, C14, & 
                            c_K, c_K1, nu1, c_K2, nu2, c_K6, nu6,  & 
                            c_K8, nu8, c_K9, nu9, c_Krrainm, nu_r, c_Ksqd, & 
                            nu_hd, gamma_coef, gamma_coefb, gamma_coefc, & 
                            mu, beta, lmin_coef, taumin, taumax )


    ! It was decided after some experimentation, that the best
    ! way to produce grid independent results is to set lmin to be
    ! some fixed value at the surface. -dschanen 21 May 2007
    !lmin = lmin_coef * deltaz  ! Old
    lmin = lmin_coef * 40.0 ! New fixed value

    T0       = T0_in
    ts_nudge = ts_nudge_in

    hydromet_dim = hydromet_dim_in
    sclr_dim     = sclr_dim_in

    ! In a tuning run, this array has the potential to be allocated already
    if ( .not. allocated( sclrtol ) ) then
       allocate( sclrtol(1:sclr_dim) )
    else
       deallocate( sclrtol ) 
       allocate( sclrtol(1:sclr_dim) )
    endif

    sclrtol(1:sclr_dim) = sclrtol_in(1:sclr_dim)


    ! ##### Adjust Constant Diffusivity Coefficients Based On Grid Spacing #####

    ! All of the background coefficients of eddy diffusivity, as well as the
    ! constant coefficient for 4th-order hyper-diffusion, must be adjusted based
    ! on the size of the grid spacing.  For a case that uses an evenly-spaced
    ! grid, the adjustment is based on the constant grid spacing deltaz.  For a
    ! case that uses a stretched grid, the adjustment is based on avg_deltaz,
    ! which is the average grid spacing over the vertical domain.

    if ( l_implemented ) then

       ! CLUBB is implemented in a host model.

       ! Find the average deltaz over the grid based on momentum level inputs.

       avg_deltaz  &
          = ( momentum_heights(nzmax) - momentum_heights(1) )  &
            / ( nzmax - 1 )

       ! The nu's are chosen for avg_deltaz <= 40 m. Looks like they must 
       ! be adjusted for larger grid spacings (Vince Larson)

       if ( avg_deltaz > 40.0 ) then
          nu1 = nu1 * avg_deltaz / 40.0
          nu2 = nu2 * avg_deltaz / 40.0
          nu6 = nu6 * avg_deltaz / 40.0
          nu8 = nu8 * avg_deltaz / 40.0
          nu9 = nu9 * avg_deltaz / 40.0
       endif

       ! There should be a different formula for determining nu_r for different
       ! sized grid spacings.  For DYCOMS2 RF02, nu_r is set to 5.0 for the
       ! high-resolution 10 m. grid spacing and to 25.0 for the low-resolution
       ! 100 m. grid spacing.  The following equation allows for both of those
       ! parameters.  Brian.

       if ( avg_deltaz > 20.0 ) then
          nu_r = nu_r * avg_deltaz / 20.0
       endif

       ! The value of nu_hd is based on an average grid box spacing of 40 m.
       ! The value of nu_hd should be adjusted directly proportional to the
       ! average grid box size, whether the average grid box size is less than
       ! 40 m. or greater than 40 m.
       
       nu_hd = nu_hd * ( avg_deltaz / 40.0 )

    else

       ! CLUBB model is running on it's own.

       if ( grid_type == 1 ) then

          ! Evenly-spaced grid.

          ! The nu's are chosen for deltaz <= 40 m. Looks like they must 
          ! be adjusted for larger grid spacings (Vince Larson)

          if ( deltaz > 40.0 ) then
             nu1 = nu1 * deltaz / 40.0
             nu2 = nu2 * deltaz / 40.0
             nu6 = nu6 * deltaz / 40.0
             nu8 = nu8 * deltaz / 40.0
             nu9 = nu9 * deltaz / 40.0
          endif

          ! There should be a different formula for determining nu_r for
          ! different sized grid spacings.  For DYCOMS2 RF02, nu_r is set to 5.0
          ! for the high-resolution 10 m. grid spacing and to 25.0 for the
          ! low-resolution 100 m. grid spacing.  The following equation allows
          ! for both of those parameters.  Brian.

          if ( deltaz > 20.0 ) then
             nu_r = nu_r * deltaz / 20.0
          endif

          ! The value of nu_hd is based on a grid box spacing of 40 m.  The
          ! value of nu_hd should be adjusted directly proportional to the grid
          ! box size, whether the grid box size is less than 40 m. or greater
          ! than 40 m.
       
          nu_hd = nu_hd * ( deltaz / 40.0 )

       elseif ( grid_type == 2 ) then

          ! Stretched (unevenly-spaced) grid:  stretched thermodynamic level
          ! input.

          ! Find the average deltaz over the stretched grid based on
          ! thermodynamic level inputs.

          avg_deltaz  &
             = ( thermodynamic_heights(nzmax) - thermodynamic_heights(1) )  &
               / ( nzmax - 1 )

          ! The nu's are chosen for avg_deltaz <= 40 m. Looks like they must 
          ! be adjusted for larger grid spacings (Vince Larson)

          if ( avg_deltaz > 40.0 ) then
             nu1 = nu1 * avg_deltaz / 40.0
             nu2 = nu2 * avg_deltaz / 40.0
             nu6 = nu6 * avg_deltaz / 40.0
             nu8 = nu8 * avg_deltaz / 40.0
             nu9 = nu9 * avg_deltaz / 40.0
          endif

          ! There should be a different formula for determining nu_r for
          ! different sized grid spacings.  For DYCOMS2 RF02, nu_r is set to 5.0
          ! for the high-resolution 10 m. grid spacing and to 25.0 for the
          ! low-resolution 100 m. grid spacing.  The following equation allows
          ! for both of those parameters.  Brian.

          if ( avg_deltaz > 20.0 ) then
             nu_r = nu_r * avg_deltaz / 20.0
          endif

          ! The value of nu_hd is based on an average grid box spacing of 40 m.
          ! The value of nu_hd should be adjusted directly proportional to the
          ! average grid box size, whether the average grid box size is less
          ! than 40 m. or greater than 40 m.
       
          nu_hd = nu_hd * ( avg_deltaz / 40.0 )

       elseif ( grid_type == 3 ) then

          ! Stretched (unevenly-spaced) grid:  stretched momentum level input.

          ! Find the average deltaz over the stretched grid based on momentum
          ! level inputs.

          avg_deltaz  &
             = ( momentum_heights(nzmax) - momentum_heights(1) )  &
               / ( nzmax - 1 )

          ! The nu's are chosen for avg_deltaz <= 40 m. Looks like they must 
          ! be adjusted for larger grid spacings (Vince Larson)

          if ( avg_deltaz > 40.0 ) then
             nu1 = nu1 * avg_deltaz / 40.0
             nu2 = nu2 * avg_deltaz / 40.0
             nu6 = nu6 * avg_deltaz / 40.0
             nu8 = nu8 * avg_deltaz / 40.0
             nu9 = nu9 * avg_deltaz / 40.0
          endif

          ! There should be a different formula for determining nu_r for
          ! different sized grid spacings.  For DYCOMS2 RF02, nu_r is set to 5.0
          ! for the high-resolution 10 m. grid spacing and to 25.0 for the
          ! low-resolution 100 m. grid spacing.  The following equation allows
          ! for both of those parameters.  Brian.

          if ( avg_deltaz > 20.0 ) then
             nu_r = nu_r * avg_deltaz / 20.0
          endif

          ! The value of nu_hd is based on an average grid box spacing of 40 m.
          ! The value of nu_hd should be adjusted directly proportional to the
          ! average grid box size, whether the average grid box size is less
          ! than 40 m. or greater than 40 m.
       
          nu_hd = nu_hd * ( avg_deltaz / 40.0 )

       endif

    endif  ! l_implemented


    ! Sanity check
    if ( beta < 0.0 .or. beta > 3.0 ) then

       ! Constraints on beta
       write(fstderr,*) "beta= ", beta
       write(fstderr,*) "beta cannot be < 0 or > 3"
       err_code = clubb_var_out_of_bounds    

    elseif ( lmin < 4.0 ) then

       ! Constraints on mixing length
       write(fstderr,*) "lmin= ", lmin
       write(fstderr,*) "lmin is < 4.0"
       err_code = clubb_var_out_of_bounds

    else

       err_code = clubb_no_error

    endif

!    write(*,nml=initvars) ! %% debug


    return

  end subroutine setup_parameters

  !=============================================================================
  subroutine read_parameters( iunit, filename, params )

    ! Description:
    ! Read a namelist containing the model parameters

    ! References:
    ! None
    !-----------------------------------------------------------------------

    implicit none

    ! Input variables
    integer, intent(in) :: iunit

    character(len=*), intent(in) :: filename

    ! Output variables
    real, intent(out), dimension(nparams) :: params

    ! If the filename is empty, assume we're using a `working' set of
    ! parameters that are set statically here (handy for host models).
    if ( filename == "" ) then
       C1          = 2.5
       C1b         = 2.5
       C1c         = 1.0
       C2rt        = 1.0
       C2thl       = 1.0
       C2rtthl     = 2.0
       C2          = 1.3
       C2b         = 1.3
       C2c         = 5.0
       C4          = 5.2
       C5          = 0.3 
       C6rt        = 6.0 
       C6rtb       = 6.0
       C6rtc       = 1.0
       C6thl       = 6.0 
       C6thlb      = 6.0 
       C6thlc      = 1.0 
       C7          = 0.1 
       C7b         = 0.8 
       C7c         = 0.5 
       C8          = 3.0  
       C8b         = 0.005
       C10         = 3.3
       C11         = 0.75
       C11b        = 0.35
       C11c        = 0.5
       C12         = 1.0
       C13         = 0.1
       C14         = 1.0
       c_K         = 0.548
       c_K1        = 0.0
       nu1         = 20.0
       c_K2        = 0.0
       nu2         = 5.0
       c_K6        = 0.0
       nu6         = 5.0
       c_K8        = 0.5
       nu8         = 20.0
       c_K9        = 0.0
       nu9         = 20.0
       c_Krrainm   = 0.075
       nu_r        = 3.0
       c_Ksqd      = 10.0
       nu_hd       = 100000.0
       beta        = 1.75
       gamma_coef  = 0.32
       gamma_coefb = 0.32
       gamma_coefc = 5.0
       taumin      = 90.0
       taumax      = 3600.0
       lmin_coef   = 0.5
       mu          = 6.000E-4

    else

      ! Read the namelist
      open(unit=iunit, file=filename, status='old', action='read')

      read(unit=iunit, nml=initvars)

      close(unit=iunit)

    endif

    ! Put the variables in the output array 
    call pack_parameters( C1, C1b, C1c, C2, C2b, C2c, C2rt, C2thl, C2rtthl, &
                          C4, C5, C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc, &
                          C7, C7b, C7c, C8, C8b, C10, & 
                          C11, C11b, C11c, C12, C13, C14, & 
                          c_K, c_K1, nu1, c_K2, nu2, c_K6, nu6,  & 
                          c_K8, nu8, c_K9, nu9, c_Krrainm, nu_r, c_Ksqd, & 
                          nu_hd, gamma_coef, gamma_coefb, gamma_coefc, & 
                          mu, beta, lmin_coef, taumin, taumax, params )

    return

  end subroutine read_parameters

  !=============================================================================
  subroutine read_param_spread & 
           ( iunit, filename, nindex, param_spread, ndim )

    ! Description:
    ! Read a namelist containing the amount to vary model parameters.
    ! Used by the downhill simplex / simulated annealing algorithm.

    ! References:
    ! None
    !-----------------------------------------------------------------------

    implicit none

    ! Input variables
    integer, intent(in) :: iunit

    character(len=*), intent(in) :: filename

    ! Output variables

    ! An array of array indices (i.e. which elements of the array `params'
    ! are contained within the simplex and the spread variable)
    integer, intent(out), dimension(nparams) :: nindex  

    real, intent(out), dimension(nparams) ::  & 
      param_spread  ! Amount to vary the parameter in the initial simplex

    integer, intent(out) :: ndim  ! Dimension of the init simplex

    ! Local variables
    integer :: i

    ! Amount to change each parameter for the initial simplex
    ! This MUST be changed to match the initvars namelist if parameters are added!
    namelist /initspread/  & 
       C1, C1b, C1c, C2, C2b, C2c,  & 
       C2rt, C2thl, C2rtthl, C4, C5, & 
       C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc, & 
       C7, C7b, C7c, C8, C8b, C10, C11, C11b, C11c, & 
       C12, C13, C14, c_K, c_K1, nu1, c_K2, nu2,  & 
       c_K6, nu6, c_K8, nu8, c_K9, nu9, c_Krrainm, nu_r, c_Ksqd,  & 
       nu_hd, beta, gamma_coef, gamma_coefb, gamma_coefc, & 
       lmin_coef, taumin, taumax, mu

    ! Read the namelist
    open(unit=iunit, file=filename, status='old', action='read')

    read(unit=iunit, nml=initspread)

    close(unit=iunit)

    ! Put the variables in the output array 
    call pack_parameters( C1, C1b, C1c, C2, C2b, C2c, C2rt, C2thl, C2rtthl, &
                          C4, C5, C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc, &
                          C7, C7b, C7c, C8, C8b, C10, & 
                          C11, C11b, C11c, C12, C13, C14, & 
                          c_K, c_K1, nu1, c_K2, nu2, c_K6, nu6,  & 
                          c_K8, nu8, c_K9, nu9, c_Krrainm, nu_r, c_Ksqd, & 
                          nu_hd, gamma_coef, gamma_coefb, gamma_coefc, & 
                          mu, beta, lmin_coef, taumin, taumax, param_spread )

    ! Initialize to zero
    nindex(1:nparams) = 0
    ndim = 0

    ! Determine how many variables are being changed
    do i = 1, nparams, 1

       if ( param_spread(i) /= 0.0 ) then
          ndim = ndim + 1   ! Increase the total
          nindex(ndim) = i  ! Set the next array index
       endif

    enddo

    return

  end subroutine read_param_spread

  !=============================================================================
  subroutine pack_parameters &
             ( C1, C1b, C1c, C2, C2b, C2c, C2rt, C2thl, C2rtthl, &
               C4, C5, C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc, &
               C7, C7b, C7c, C8, C8b, C10, &
               C11, C11b, C11c, C12, C13, C14, &
               c_K, c_K1, nu1, c_K2, nu2, c_K6, nu6,  &
               c_K8, nu8, c_K9, nu9, c_Krrainm, nu_r, c_Ksqd, &
               nu_hd, gamma_coef, gamma_coefb, gamma_coefc, &
               mu, beta, lmin_coef, taumin, taumax, params )

    ! Description:
    ! Takes the list of scalar variables and puts them into a 1D vector.
    ! It is here for the purpose of keeping the code generalized
    ! when new variables are added.

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use parameter_indices, only: & 
        iC1,  & ! Variable(s)
        iC1b, & 
        iC1c, & 
        iC2, & 
        iC2b, & 
        iC2c, & 
        iC2rt, & 
        iC2thl, & 
        iC2rtthl, & 
        iC4, & 
        iC5, & 
        iC6rt, & 
        iC6rtb, & 
        iC6rtc, & 
        iC6thl, & 
        iC6thlb, & 
        iC6thlc, & 
        iC7, & 
        iC7b, & 
        iC7c, & 
        iC8, & 
        iC8b, & 
        iC10, & 
        iC11, & 
        iC11b, & 
        iC11c, & 
        iC12, & 
        iC13, & 
        iC14, & 
        ic_K,  & 
        ic_K1, & 
        inu1, & 
        ic_K2, & 
        inu2, & 
        ic_K6, & 
        inu6, & 
        ic_K8, & 
        inu8, & 
        ic_K9, & 
        inu9, & 
        ic_Krrainm, & 
        inu_r, & 
        ic_Ksqd, &
        inu_hd, & 
        igamma_coef, & 
        igamma_coefb, & 
        igamma_coefc, & 
        imu, & 
        ibeta, & 
        ilmin_coef, & 
        itaumin, & 
        itaumax, & 
        nparams

     implicit none

    ! Input variables
    real, intent(in) :: & 
      C1, C1b, C1c, C2, C2b, C2c, C2rt, C2thl, C2rtthl, & 
      C4, C5, C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc, & 
      C7, C7b, C7c, C8, C8b, C10, & 
      C11, C11b, C11c, C12, C13, C14, & 
      c_K, c_K1, nu1, c_K2, nu2, c_K6, nu6, c_K8, nu8,  & 
      c_K9, nu9, c_Krrainm, nu_r, c_Ksqd, nu_hd, gamma_coef, &
      gamma_coefb, gamma_coefc, mu, beta, lmin_coef, taumin, taumax

    ! Output variables
    real, intent(out), dimension(nparams) :: params

    params(iC1)      = C1
    params(iC1b)     = C1b
    params(iC1c)     = C1c
    params(iC2)      = C2
    params(iC2b)     = C2b
    params(iC2c)     = C2c
    params(iC2rt)    = C2rt
    params(iC2thl)   = C2thl
    params(iC2rtthl) = C2rtthl
    params(iC4)      = C4
    params(iC5)      = C5
    params(iC6rt)    = C6rt
    params(iC6rtb)   = C6rtb
    params(iC6rtc)   = C6rtc
    params(iC6thl)   = C6thl
    params(iC6thlb)  = C6thlb
    params(iC6thlc)  = C6thlc
    params(iC7)      = C7
    params(iC7b)     = C7b
    params(iC7c)     = C7c
    params(iC8)      = C8
    params(iC8b)     = C8b
    params(iC10)     = C10
    params(iC11)     = C11
    params(iC11b)    = C11b
    params(iC11c)    = C11c
    params(iC12)     = C12
    params(iC13)     = C13
    params(iC14)     = C14

    params(ic_K)       = c_K 
    params(ic_K1)      = c_K1
    params(inu1)       = nu1
    params(ic_K2)      = c_K2
    params(inu2)       = nu2
    params(ic_K6)      = c_K6
    params(inu6)       = nu6
    params(ic_K8)      = c_K8
    params(inu8)       = nu8
    params(ic_K9)      = c_K9
    params(inu9)       = nu9
    params(ic_Krrainm) = c_Krrainm
    params(inu_r)      = nu_r
    params(ic_Ksqd)    = c_Ksqd
    params(inu_hd)     = nu_hd

    params(igamma_coef)  = gamma_coef
    params(igamma_coefb) = gamma_coefb
    params(igamma_coefc) = gamma_coefc

    params(imu) = mu

    params(ibeta) = beta

    params(ilmin_coef) = lmin_coef

    params(itaumin) = taumin
    params(itaumax) = taumax

    return
  end subroutine pack_parameters

  !=============================================================================
  subroutine unpack_parameters & 
             ( params, & 
               C1, C1b, C1c, C2, C2b, C2c, C2rt, C2thl, C2rtthl, & 
               C4, C5, C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc, & 
               C7, C7b, C7c, C8, C8b, C10, & 
               C11, C11b, C11c, C12, C13, C14, & 
               c_K, c_K1, nu1, c_K2, nu2, c_K6, nu6, & 
               c_K8, nu8, c_K9, nu9, c_Krrainm, nu_r, c_Ksqd, & 
               nu_hd, gamma_coef, gamma_coefb, gamma_coefc, & 
               mu, beta, lmin_coef, taumin, taumax )

    ! Description:
    ! Takes the 1D vector and returns the list of scalar variables.
    ! Here for the purposes of keeping the code generalized
    ! when new variables are added.

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use parameter_indices, only: & 
        iC1,  & ! Variable(s)
        iC1b, & 
        iC1c, & 
        iC2, & 
        iC2b, & 
        iC2c, & 
        iC2rt, & 
        iC2thl, & 
        iC2rtthl, & 
        iC4, & 
        iC5, & 
        iC6rt, & 
        iC6rtb, & 
        iC6rtc, & 
        iC6thl, & 
        iC6thlb, & 
        iC6thlc, & 
        iC7, & 
        iC7b, & 
        iC7c, & 
        iC8, & 
        iC8b, & 
        iC10, & 
        iC11, & 
        iC11b, & 
        iC11c, & 
        iC12, & 
        iC13, & 
        iC14, & 
        ic_K,  & 
        ic_K1, & 
        inu1, & 
        ic_K2, & 
        inu2, & 
        ic_K6, & 
        inu6, & 
        ic_K8, & 
        inu8, & 
        ic_K9, & 
        inu9, & 
        ic_Krrainm, & 
        inu_r, & 
        ic_Ksqd, & 
        inu_hd, & 
        igamma_coef, & 
        igamma_coefb, & 
        igamma_coefc, & 
        imu, & 
        ibeta, & 
        ilmin_coef, & 
        itaumin, & 
        itaumax, & 
        nparams

    implicit none

    ! Input variables
    real, intent(in), dimension(nparams) :: params

    ! Output variables
    real, intent(out) :: & 
      C1, C1b, C1c, C2, C2b, C2c, C2rt, C2thl, C2rtthl, & 
      C4, C5, C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc, & 
      C7, C7b, C7c, C8, C8b, C10, & 
      C11, C11b, C11c, C12, C13, C14, & 
      c_K, c_K1, nu1, c_K2, nu2, c_K6, nu6, & 
      c_K8, nu8, c_K9, nu9, c_Krrainm, nu_r, c_Ksqd, & 
      nu_hd, gamma_coef, gamma_coefb, gamma_coefc, & 
      mu, beta, lmin_coef, taumin, taumax

    C1      = params(iC1)
    C1b     = params(iC1b)
    C1c     = params(iC1c)
    C2      = params(iC2)
    C2b     = params(iC2b)
    C2c     = params(iC2c)
    C2rt    = params(iC2rt)
    C2thl   = params(iC2thl)
    C2rtthl = params(iC2rtthl)
    C4      = params(iC4)
    C5      = params(iC5)
    C6rt    = params(iC6rt)
    C6rtb   = params(iC6rtb)
    C6rtc   = params(iC6rtc)
    C6thl   = params(iC6thl)
    C6thlb  = params(iC6thlb)
    C6thlc  = params(iC6thlc)
    C7      = params(iC7)
    C7b     = params(iC7b)
    C7c     = params(iC7c)
    C8      = params(iC8)
    C8b     = params(iC8b)
    C10     = params(iC10)
    C11     = params(iC11)
    C11b    = params(iC11b)
    C11c    = params(iC11c)
    C12     = params(iC12)
    C13     = params(iC13)
    C14     = params(iC14)

    c_K       = params(ic_K)
    c_K1      = params(ic_K1)
    nu1       = params(inu1)
    c_K2      = params(ic_K2)
    nu2       = params(inu2)
    c_K6      = params(ic_K6)
    nu6       = params(inu6)
    c_K8      = params(ic_K8)
    nu8       = params(inu8)
    c_K9      = params(ic_K9)
    nu9       = params(inu9)
    c_Krrainm = params(ic_Krrainm)
    nu_r      = params(inu_r)
    c_Ksqd    = params(ic_Ksqd)
    nu_hd     = params(inu_hd)

    gamma_coef  = params(igamma_coef)
    gamma_coefb = params(igamma_coefb)
    gamma_coefc = params(igamma_coefc)

    mu = params(imu)

    beta = params(ibeta)

    lmin_coef = params(ilmin_coef)

    taumin = params(itaumin)
    taumax = params(itaumax)

    return
  end subroutine unpack_parameters

  !=============================================================================
  subroutine get_parameters( params )

    ! Description:
    ! Return an array of all tunable parameters

    ! References:
    ! None
    !-----------------------------------------------------------------------

    implicit none

    ! Input Variables
    real, intent(out), dimension(nparams) :: params

    call pack_parameters( C1, C1b, C1c, C2, C2b, C2c, C2rt, C2thl, C2rtthl, &
                          C4, C5, C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc, &
                          C7, C7b, C7c, C8, C8b, C10, & 
                          C11, C11b, C11c, C12, C13, C14, & 
                          c_K, c_K1, nu1, c_K2, nu2, c_K6, nu6,  & 
                          c_K8, nu8, c_K9, nu9, c_Krrainm, nu_r, c_Ksqd, & 
                          nu_hd, gamma_coef, gamma_coefb, gamma_coefc, & 
                          mu, beta, lmin_coef, taumin, taumax, params )

    return

  end subroutine get_parameters

!===============================================================================

end module parameters_tunable
