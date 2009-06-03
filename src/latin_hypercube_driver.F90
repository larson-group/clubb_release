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

    use lh_sampler_mod, only: & 
      lh_sampler ! Procedure

    use micro_calcs_mod, only: & 
      micro_calcs ! Procedure

    use variables_prognostic_module, only: &
      pdf_parameter  ! Type

    use constants, only: fstderr

    use variables_diagnostic_module, only: & 
      AKm_est,  & 
      AKm, & 
      AKstd, & 
      AKstd_cld, & 
      AKm_rcm, & 
      AKm_rcc, & 
      rcm_est

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

    ! Local variables
    integer :: p_matrix(n_micro_call,d_variables+1)

    ! Sample drawn from uniform distribution
    double precision, dimension(nnzp,n_micro_call,(d_variables+1)) :: X_u

    ! Sample that is transformed ultimately to normal-lognormal
    double precision, dimension(nnzp,n_micro_call,d_variables) :: X_nl

    integer :: i_rmd, k

    ! A true/false flag that determines whether the PDF allows us to construct a sample
    logical, dimension(nnzp) :: l_sample_flag 

    ! ---- Begin Code ----

    if ( .not. allocated( height_time_matrix ) ) then
      allocate( height_time_matrix(nnzp, nt_repeat, d_variables+1) )
    end if

    ! Latin hypercube sample generation
    ! Generate height_time_matrix, an nnzp x nt_repeat x d_variables array of random integers
    i_rmd = mod( iter-1, sequence_length )

    if ( i_rmd == 0 ) then
      call permute_height_time( nnzp, nt_repeat, d_variables+1, & ! intent(in)
                                height_time_matrix )              ! intent(out)
    end if
    ! End Latin hypercube sample generation

    ! print*, 'latin_hypercube_driver: i_rmd=', i_rmd
      !--------------------------------------------------------------
      ! Latin hypercube sampling
      !--------------------------------------------------------------
    if ( iirrainm < 1 ) then
      write(fstderr,*) "Latin hypercube sampling is enabled, but there is"// &
      " no rain water mixing ratio being predicted."
      stop
    end if

    do k = 2, nnzp
      ! Choose which rows of LH sample to feed into closure.
      p_matrix(1:n_micro_call,1:(d_variables+1)) = &
        height_time_matrix(k,n_micro_call*i_rmd+1:n_micro_call*i_rmd+n_micro_call,1:d_variables+1)

      ! print*, 'latin_hypercube_sampling: got past p_matrix'

      ! Generate LH sample, represented by X_u and X_nl, for level k
      call lh_sampler( n_micro_call, nt_repeat, d_variables, p_matrix, & ! intent(in)
                       cf(k), pdf_params, k, &                           ! intent(in)
                       hydromet(k,iirrainm), &                           ! intent(in)
                       X_u(k,:,:), X_nl(k,:,:), l_sample_flag(k) ) ! intent(out)

      ! print *, 'latin_hypercube_sampling: got past lh_sampler'
    end do ! 2..nnzp

      ! Perform LH and analytic microphysical calculations
    call micro_calcs( nnzp, n_micro_call, d_variables, X_u, X_nl, & ! intent(in)
                      l_sample_flag, pdf_params, &                  ! intent(in)
                      AKm_est, AKm, AKstd, AKstd_cld, &             ! intent(out)
                      AKm_rcm, AKm_rcc, rcm_est )                   ! intent(out)

    ! print*, 'latin_hypercube_driver: AKm=', AKm
    ! print*, 'latin_hypercube_driver: AKm_est=', AKm_est

    return
  end subroutine latin_hypercube_driver

end module latin_hypercube_mod
