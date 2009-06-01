! $Id$
!-------------------------------------------------------------------------------
module latin_hypercube_mod

  implicit none
   
  public :: latin_hypercube_driver

  private ! Default scope

  integer, parameter :: &
    d_variables     = 5,  & ! Number of variables to sample
    n_micro_call    = 12, & ! Number of calls to microphysics per timestep (normally=2)
    sequence_length = 1     ! nt_repeat/n_micro_call; number of timesteps before sequence repeats.

  ! Number of random samples before sequence of repeats (normally=10)
  integer, parameter :: &
    nt_repeat = n_micro_call * sequence_length 

  ! A true/false flag that determines whether
  !     the PDF allows us to construct a sample
  logical :: sample_flag

! integer, dimension(gr%nnzp, nt_repeat, d_variables+1) :: & 
  integer, allocatable, dimension(:,:,:) :: & 
    height_time_matrix ! matrix of rand ints

  contains

!-------------------------------------------------------------------------------
  subroutine latin_hypercube_driver( iter, nnzp, pdf_params, hydromet, cf )
! Description:
!   Do latin hypercube sampling
! References:
!   None
!-------------------------------------------------------------------------------
    use array_index, only: iirrainm ! Variable
    
    use parameters_model, only: hydromet_dim ! Variable

    use permute_height_time_mod, only: & 
      permute_height_time ! Procedure

    use variables_prognostic_module, only: &
      pdf_parameter  ! Type

    use constants, only: fstderr

    implicit none

    ! External
    intrinsic :: allocated, mod

    ! Input Variables
    integer, intent(in) :: &
      iter, & ! Model iteration number
      nnzp    ! Domain dimension

    type(pdf_parameter), intent(in) :: pdf_params

    real, dimension(nnzp,hydromet_dim), intent(in) :: &
      hydromet ! Hydrometeor Species    [units vary]

    real, dimension(nnzp), intent(in) :: cf ! Cloud fraction [%]

    integer :: i_rmd, k

    ! ---- Begin Code ----

    if ( .not. allocated( height_time_matrix ) ) then
      allocate( height_time_matrix(nnzp, nt_repeat, d_variables+1) )
    end if

    ! Latin hypercube sample generation
    ! Generate height_time_matrix, an nnzp x nt_repeat x d_variables array of random integers
    i_rmd = mod( iter-1, sequence_length )

    if ( i_rmd == 0 ) then
      call permute_height_time( nnzp, nt_repeat, d_variables+1, & ! intent(in)
                                height_time_matrix )                 ! intent(out)
    end if
    ! End Latin hypercube sample generation

    ! print*, 'advance_clubb_core: i_rmd=', i_rmd
      !--------------------------------------------------------------
      ! Latin hypercube sampling
      !--------------------------------------------------------------
    if ( iirrainm < 1 ) then
      write(fstderr,*) "Latin hypercube sampling is enabled, but there is"// &
      " no rain water mixing ratio being predicted."
      stop
    end if

    do k = 2, nnzp, 1

      call latin_hypercube_sampling &
           ( k, n_micro_call, d_variables, nt_repeat, i_rmd, &
             pdf_params, hydromet(:,iirrainm), &
             cf, nnzp, sample_flag, height_time_matrix )

    ! print*, 'advance_clubb_core: AKm=', AKm
    ! print*, 'advance_clubb_core: AKm_est=', AKm_est

    end do


    return
  end subroutine latin_hypercube_driver

  !-----------------------------------------------------------------------
  subroutine latin_hypercube_sampling & 
             ( k, n, dvar, nt, i_rmd, & 
               pdf_params, & 
               rrainm, cf, nnzp, l_sflag, height_time_matrix )
    ! Description:
    !   Estimate using Latin Hypercubes.  This is usually disabled by default.
    !   The actual generation of a random matrix is done in a call from the
    !   subroutine advance_clubb_core to permute_height_time()
    ! References:
    !   None
    !-----------------------------------------------------------------------

    use variables_diagnostic_module, only: & 
        AKm_est,  & 
        AKm, & 
        AKstd, & 
        AKstd_cld, & 
        AKm_rcm, & 
        AKm_rcc, & 
        rcm_est

    use variables_prognostic_module, only: & 
      pdf_parameter

    use lh_sampler_mod, only: & 
        lh_sampler ! Procedure

    use micro_calcs_mod, only: & 
        micro_calcs ! Procedure

    implicit none


    ! Input Variables
    integer, intent(in) :: k  ! index
    integer, intent(in) :: n, dvar, i_rmd, nt, nnzp

    type(pdf_parameter), intent(in) :: &
      pdf_params ! PDF parameters [units vary]

    real, dimension(nnzp), intent(in) ::  & 
      rrainm,  & ! Rain water mixing ratio  [kg/kg]
      cf         ! Cloud fraction           [%]

    integer, dimension(nnzp, nt, (dvar+1) ), intent(in) :: & 
      height_time_matrix ! matrix of rand ints

    ! Output Variables
    logical, intent(out) :: l_sflag

    ! Local Variables

    integer :: p_matrix(n, dvar+1)
    ! Sample drawn from uniform distribution
    double precision, dimension(n,(dvar+1)) :: X_u

    ! Sample that is transformed ultimately to normal-lognormal
    double precision, dimension(n,dvar) :: X_nl

    ! Choose which rows of LH sample to feed into closure.
    p_matrix(1:n,1:(dvar+1)) = & 
      height_time_matrix( k,n*i_rmd+1:n*i_rmd+n, 1:(dvar+1) )

    ! print*, 'latin_hypercube_sampling: got past p_matrix'

    ! Generate LH sample, represented by X_u and X_nl, for level k
    call lh_sampler( n, nt, dvar, p_matrix,       & ! intent(in)
                     cf(k), pdf_params, k,        & ! intent(in)
                     rrainm(k),                   & ! intent(in)
                     X_u, X_nl, l_sflag )           ! intent(out)

    !print *, 'latin_hypercube_sampling: got past lh_sampler'

    ! Perform LH and analytic microphysical calculations
    call micro_calcs( n, dvar, X_u, X_nl, l_sflag,                & ! intent(in)
                      pdf_params, k,                              & ! intent(in)
                      AKm_est(k), AKm(k), AKstd(k), AKstd_cld(k), & ! intent(out)
                      AKm_rcm(k), AKm_rcc(k), rcm_est(k) )          ! intent(out)

    !print*, 'k, AKm_est=', k, AKm_est(k)
    !print*, 'k, AKm=', k, AKm(k)

    return
  end subroutine latin_hypercube_sampling
end module latin_hypercube_mod
