!-----------------------------------------------------------------------
! $Id: budget_terms.F,v 1.13 2005-09-07 17:17:17 dschanen Exp $     

! MODULE budget_terms

! Description:
! functions and subroutines for tuning for budget terms in the HOC model.

! The general structure for executing a tuning run is as follows:

! subroutine setup_budget_terms ( passed namelist, returns ftol )

! subroutine <varname>_run ( passed function to minimize )
!   |
!   +--- subroutine amoeba ( passed function to minimize )
!                 |
!                 +--- function min( constants to optimize )
!                      ( looped over until ftol is met )
!
! subroutine output_optimal_constants( )  

!-----------------------------------------------------------------------
      module budget_terms

        use inputfields

        implicit none

!-----------------------------------------------------------------------
!       Interfaces

!  <var>_run
!  |
!  +--- amoeba
!       |    
!       +--- wp2_min/wp3_min/wprtp_min...
!            |
!            +--- term_min
!
!  buoy_advect_run
!  |
!  +--- amoeba
!       |    
!       +--- wp2_buoy_min/wp3_buoy_min/wprtp_buoy_min...
!            |
!            +--- gen_betagamma_min
!                 |
!                 +--- term_min
!       or 
!       +--- wp3_advect_min/wprtp_advect_min/wpthlp_advect_min...
!            |
!            +--- gen_betagamma_min
!                 |
!                 +--- term_min
!-----------------------------------------------------------------------

!       Interfaces

        public wpn_run, wpxp_run, buoy_advect_run, xapxbp_run

        public wp2_min, wp3_min, wprtp_min, wpthlp_min

        public wp2_buoy_min, wp3_buoy_min, 
     .         wprtp_buoy_min, wpthlp_buoy_min

        public wp3_advect_min, wprtp_advect_min, wpthlp_advect_min

        public rtp2_advect_min, thlp2_advect_min, rtpthlp_advect_min
        public rtp2_min, thlp2_min, rtpthlp_min

        public output_optimal_constants

        private term_min, gen_betagamma_min
        private calc_les_wpthlp
        private calc_les_rtp2_dp, calc_les_thlp2_dp, calc_les_rtpthlp_dp


!-----------------------------------------------------------------------
!       Internal Run Information

        character(len=50), private :: run_file       ! namelist
        character(len=50), private :: hoc_stats_file ! GRaDS data
        character(len=50), private :: les_stats_file ! GRaDS data

        integer, dimension(10), private :: times_nl

        real, private :: sample_ratio
        real, private :: ftol

        integer, allocatable, dimension(:), private :: times
        integer, private :: z_i
        integer, private :: z_f

!       diff_array is set from hoc_tuner_budget_terms, it is needed 
!       by *_min if isValid == .false.
        real, allocatable, dimension(:), private   :: rand_vect
        real, allocatable, dimension(:,:), private :: const_array ! 'p' matrix in amoeba
        real, allocatable, dimension(:), private   :: diff_array  ! 'y' vector in amoeba

        real, private :: C1, C2rt, C2thl, C2rtthl, C4, C5, C6rt, C6thl
        real, private :: C7, C8, C10, C11
        real, private :: nu1, nu2, nu6, nu8, nu_r
        real, private :: beta, gamma_coef, c_K, lmin_coef
        real, private :: taumin, taumax, mu

        real, private :: C1_var, C2rt_var, C2thl_var, C2rtthl_var
        real, private :: C4_var, C5_var, C6rt_var, C6thl_var
        real, private :: C7_var, C8_var, C10_var, C11_var
        real, private :: nu1_var, nu2_var, nu6_var, nu8_var, nu_r_var
        real, private :: beta_var, gamma_coef_var, c_K_var
        real, private :: lmin_coef_var, taumin_var, taumax_var, mu_var

        logical, public :: lwp2_prdp, lwp3_prdp
        logical, public :: lwprtp_prdp, lwpthlp_prdp
        logical, public :: lrtp2_dp, lthlp2_dp, lrtpthlp_dp

        logical, public :: lwp2_bp, lwp3_bp
        logical, public :: lwprtp_bp, lwpthlp_bp

        logical, public :: lwp2_tp, lwp3_tp
        logical, public :: lwprtp_tp, lwpthlp_tp
        logical, public :: lrtp2_tp, lthlp2_tp, lrtpthlp_tp

!       Internal (misc)

        logical, private :: isValid

!----------------------------------------------------------------------
!       Namelists

        namelist /initvars/ C1, C2rt, C2thl, C2rtthl, C4, C5, 
     .                      C6rt, C6thl, C7, C8, C10, C11, 
     .                      nu1, nu2, nu6, nu8, nu_r, beta, gamma_coef,
     .                      c_K, lmin_coef, taumin, taumax, mu

        namelist /variance/ C1_var, C2rt_var, C2thl_var, C2rtthl_var,
     .                      C4_var, C5_var, C6rt_var, C6thl_var,
     .                      C7_var, C8_var, C10_var, C11_var, 
     .                      nu1_var, nu2_var, nu6_var, nu8_var, 
     .                      nu_r_var, beta_var, gamma_coef_var, 
     .                      c_K_var, lmin_coef_var, taumin_var, 
     .                      taumax_var, mu_var

        namelist /case_tune/ lwp2_prdp, lwp3_prdp, 
     .                       lwprtp_prdp, lwpthlp_prdp, 
     .                       lrtp2_dp, lthlp2_dp, lrtpthlp_dp,
     .                       lwp2_bp, lwp3_bp,
     .                       lwprtp_bp, lwpthlp_bp,
     .                       lwp3_tp,
     .                       lwprtp_tp, lwpthlp_tp,
     .                       lrtp2_tp, lthlp2_tp, lrtpthlp_tp

        namelist /budget_tune/ run_file, hoc_stats_file, les_stats_file,
     .                         times_nl, z_i, z_f, ftol

        namelist /setfields/ datafile, sample_ratio, input_type,
     .                       input_um, input_vm, input_rtm, input_thlm,
     .                       input_wp2, input_wprtp, input_wpthlp, 
     .                       input_wp3, input_rtp2, input_thlp2, 
     .                       input_rtpthlp, input_upwp, input_vpwp

        private initvars, variance, budget_tune, setfields

        public case_tune ! because the elements are public, the nml must also be

        contains
!-----------------------------------------------------------------------
!  SUBROUTINE setup_budget_terms
!  Configures information about the data files and run required for the
!  other functions found in the module.
!-----------------------------------------------------------------------
        subroutine setup_budget_terms( budget_file, ftol_out )

        implicit none

!       Input

        character(len=*), intent(in) :: budget_file

!       Output
        real, intent(out) :: ftol_out

!       Internal

        integer ios

!       Initialize variables
        times_nl = 0

!       Read in the namelists
        open(10, file=trim(budget_file), status='old', iostat=ios)
        if ( ios /= 0 ) stop "budget.in open failed"
        read(10, nml=budget_tune)
        read(10, nml=initvars)
        read(10, nml=variance)
        read(10, nml=setfields)
        read(10, nml=case_tune)
        close(10)

!       Initialize module information
        allocate ( times(maxloc( times_nl, 1 )) )

        times    = times_nl(1:size(times))
        ftol_out = ftol

        isValid = .true.

        return
        end subroutine setup_budget_terms

!-----------------------------------------------------------------------
!  SUBROUTINE wpn_run( ) 
!
!  Minimize the error between HOC and LES budget terms for w'^n
!-----------------------------------------------------------------------
        subroutine wpn_run( ftol, runtype, func_min, c_update )
        use nrtype
        use nr, only: amoeba
        implicit none

!       Input 
        real, intent(in)         :: ftol
        character(*), intent(in) :: runtype

!       Determines whether the module constants are updated
        logical, intent(in)      :: c_update 

!       Internal
        real, dimension(5) :: p_elements 
        real, dimension(5) :: p_variance 

        character(10), dimension(5) :: p_string 

        real opt_error
        integer ndim, mdim
        integer iter, i, j, k

!       Interfaces
        interface
          function func_min(x)
          implicit none
          real, dimension(:), intent(in) :: x
          real func_min
          end function func_min
        end interface

!       Setup constants vector
        select case( runtype )
        case( "2" )
          p_elements = (/C1, C5, beta, gamma_coef, mu/)
          p_string = (/"C1        ", "C5        ",
     .                 "beta      ", "gamma_coef",
     .                 "mu        "/)
          p_variance = (/C1_var, C5_var, beta_var,
     .                   gamma_coef_var, mu_var/)
        case( "3" )
          p_elements = (/C8, C11, beta, gamma_coef, mu/)
          p_string = (/"C8        ", "C11       ",
     .                 "beta      ", "gamma_coef",
     .                 "mu        "/)
          p_variance = (/C8_var, C11_var, beta_var, 
     .                   gamma_coef_var, mu_var/)
        case default
          stop "wpn_run: invalid runtype"
        end select

        ndim = size( p_elements )
        mdim = size( p_elements ) + 1
 
!       Allocate P, Y, and an random factor to generate the starting simplex
        allocate( const_array(mdim, ndim) )
        allocate( diff_array(mdim) )
        allocate( rand_vect(ndim) )

        const_array(1, 1:ndim) = p_elements(1:ndim)

!       Randomize P columns
        do j = 1, ndim
          call random_number( rand_vect(1:ndim) )
          do i = 2, mdim
            const_array(i,j) = const_array(1,j)*
     .                       ( 1.0 - p_variance(j)
     .                         + rand_vect(i-1)*p_variance(j)*2 )
          enddo
        enddo

! %% debug
!        print '(a)', "const_array"
!        print '(5f8.4)', const_array
!        pause
! %% debug

!       Calculate Y-vector
        do i=1, mdim
          diff_array(i) = func_min( const_array(i,:) )
          if ( .not. isValid ) stop "Initial constants must be valid"
        enddo

        write(*,*) "diff_array:"
        do i=1, mdim
          write(*,'(f12.5)') diff_array(i)
        enddo

!       Attempt to find the optimal P-matrix
        call amoeba( const_array, diff_array, ftol, 
     .               func_min, iter )

!       Write output
        write(*,'(a12,i5)') "Iterations: ", iter
        write(*,'(a12,a12,a12)') "Constant", "Initial", "Optimal"
        do i=1, ndim
          write(*,'(a12,2f12.5)') p_string(i)//"= ", 
     .                            p_elements(i), const_array(1,i)
        enddo
!       Update module constants on c_update
        if ( c_update ) then
          select case( runtype )
          case( "2" )
            C1         = const_array(1,1)
            C5         = const_array(1,2)
            beta       = const_array(1,3)
            gamma_coef = const_array(1,4)
            mu         = const_array(1,5)
          case( "3" )
            C8         = const_array(1,1)
            C11        = const_array(1,2)
            beta       = const_array(1,3)
            gamma_coef = const_array(1,4)
            mu         = const_array(1,5)
          end select
        else
!         Generate GRaDS with Optimal constants
          write(*,*) "Running with optimal constants."      
          opt_error = func_min( const_array(1,:) ) 
          write(*,'(a15,f12.5)') "Optimal error=", opt_error

        endif

!       Deallocate & return
        deallocate( const_array, diff_array, rand_vect )

        return
        end subroutine wpn_run

!-----------------------------------------------------------------------
!  SUBROUTINE buoy_advect_run( ) 
!
!  Minimize the error between HOC and LES budget for the 
!  buoyancy or advection terms
!-----------------------------------------------------------------------
        subroutine buoy_advect_run( ftol, runtype, func_min, c_update )
        use nrtype
        use nr, only: amoeba
        implicit none

!       Input 
        real, intent(in)         :: ftol

        ! Isn't actually used, it is just here for consistency
        character(*), intent(in) :: runtype 

!       Determines whether the module constants are updated
        logical, intent(in)      :: c_update 

!       Internal
        real, dimension(2) :: p_elements 
        real, dimension(2) :: p_variance 

        character(10), dimension(2) :: p_string 

        real opt_error
        integer ndim, mdim
        integer iter, i, j, k

!       Interfaces
        interface
          function func_min(x)
          implicit none
          real, dimension(:), intent(in) :: x
          real func_min
          end function func_min
        end interface
        
!       Setup constants vector
        p_elements = (/beta, gamma_coef/)
        p_string   = (/"beta      ", "gamma_coef"/)
        p_variance = (/beta_var, gamma_coef_var/)

        ndim = size( p_elements )
        mdim = size( p_elements ) + 1
 
!       Allocate P, Y, and an random factor to generate the starting simplex
        allocate( const_array(mdim, ndim) )
        allocate( diff_array(mdim) )
        allocate( rand_vect(ndim) )

        const_array(1, 1:ndim) = p_elements(1:ndim)

!       Randomize P columns
        do j = 1, ndim
          call random_number( rand_vect(1:ndim) )
          do i = 2, mdim
            const_array(i,j) = const_array(1,j)*
     .                       ( 1.0 - p_variance(j)
     .                         + rand_vect(i-1)*p_variance(j)*2 )
          enddo
        enddo

!       Calculate Y-vector
        do i=1, mdim
          diff_array(i) = func_min( const_array(i,:) )
          if ( .not. isValid ) stop "Initial constants must be valid"
        enddo

        write(*,*) "diff_array:"
        do i=1, mdim
          write(*,'(f12.5)') diff_array(i)
        enddo

!       Attempt to find the optimal P-matrix
        call amoeba( const_array, diff_array, ftol, 
     .               func_min, iter )

!       Write output
        write(*,'(a12,i5)') "Iterations: ", iter
        write(*,'(a12,a12,a12)') "Constant", "Initial", "Optimal"
        do i=1, ndim
          write(*,'(a12,2f12.5)') p_string(i)//"= ", 
     .                            p_elements(i), const_array(1,i)
        enddo


!       Update module constants on c_update
        if ( c_update ) then
          beta       = const_array(1,1)
          gamma_coef = const_array(1,2)
        else
!         Generate GRaDS with Optimal constants
          write(*,*) "Running with optimal constants."      
          opt_error = func_min( const_array(1,:) ) 
          write(*,'(a15,f12.5)') "Optimal error=", opt_error
        endif

!       Deallocate & return
        deallocate( const_array, diff_array, rand_vect )

        return
        end subroutine buoy_advect_run

!-----------------------------------------------------------------------
!  SUBROUTINE wpxp_run( ) 
!
!  Minimize the error between HOC and LES budget terms for w'rt' or w'th_l
!-----------------------------------------------------------------------
        subroutine wpxp_run( ftol, runtype, func_min, c_update )
        use nrtype
        use nr, only: amoeba
        implicit none

!       Input 
        real, intent(in)         :: ftol
        character(*), intent(in) :: runtype

!       Determines whether the module constants are updated
        logical, intent(in)      :: c_update 

!       Internal
        real, dimension(4) :: p_elements 
        real, dimension(4) :: p_variance 

        character(10), dimension(4) :: p_string 

        real opt_error
        integer ndim, mdim
        integer iter, i, j, k

!       Interfaces
        interface
          function func_min(x)
          implicit none
          real, dimension(:), intent(in) :: x
          real func_min
          end function func_min
        end interface
        
!       Setup constants vector
        select case( runtype )
        case( "rt" )
          p_elements = (/C6rt, C7, beta, gamma_coef/)
          p_string   = (/"C6rt      ", "C7        ",
     .                   "beta      ", "gamma_coef"/)
          p_variance = (/C6rt_var, C7_var, beta_var, gamma_coef_var/)
        case( "thl" )
          p_elements = (/C6thl, C7, beta, gamma_coef/)
          p_string   = (/"C6thl     ", "C7        ",
     .                   "beta      ", "gamma_coef"/)
          p_variance = (/C6thl_var, C7_var, beta_var, gamma_coef_var/)
        case default
          stop "wpxp_run: invalid runtype"
        end select

        ndim = size( p_elements )
        mdim = size( p_elements ) + 1
 
!       Allocate P, Y, and an random factor to generate the starting simplex
        allocate( const_array(mdim, ndim) )
        allocate( diff_array(mdim) )
        allocate( rand_vect(ndim) )

        const_array(1, 1:ndim) = p_elements(1:ndim)

!       Randomize P columns
        do j = 1, ndim
          call random_number( rand_vect(1:ndim) )
          do i = 2, mdim
            const_array(i,j) = const_array(1,j)*
     .                       ( 1.0 - p_variance(j)
     .                         + rand_vect(i-1)*p_variance(j)*2 )
          enddo
        enddo

!       Calculate Y-vector
        do i=1, mdim
          diff_array(i) = func_min( const_array(i,:) )
          if ( .not. isValid ) stop "Initial constants must be valid"
        enddo

        write(*,*) "diff_array:"
        do i=1, mdim
          write(*,'(f12.5)') diff_array(i)
        enddo

!       Attempt to find the optimal P-matrix
        call amoeba( const_array, diff_array, ftol, 
     .               func_min, iter )

!       Write output
        write(*,'(a12,i5)') "Iterations: ", iter
        write(*,'(a12,a12,a12)') "Constant", "Initial", "Optimal"
        do i=1, ndim
          write(*,'(a12,2f12.5)') p_string(i)//"= ", 
     .                            p_elements(i), const_array(1,i)
        enddo


!       Update module constants on c_update
        if ( c_update ) then
          select case( runtype )
          case( "rt" )
            C6rt       = const_array(1,1)
            C7         = const_array(1,2)
            beta       = const_array(1,3)
            gamma_coef = const_array(1,4)
          case( "thl" )
            C6thl      = const_array(1,1)
            C7         = const_array(1,2)
            beta       = const_array(1,3)
            gamma_coef = const_array(1,4)
          end select
        else
!         Generate GRaDS with Optimal constants
          write(*,*) "Running with optimal constants."      
          opt_error = func_min( const_array(1,:) ) 
          write(*,'(a15,f12.5)') "Optimal error=", opt_error
        endif

!       Deallocate & return
        deallocate( const_array, diff_array, rand_vect )

        return
        end subroutine wpxp_run

!-----------------------------------------------------------------------
!  FUNCTION xapxbp_run( )
!  Calculate the total difference between the two GRaDS profiles
!  This particular routine does the variance and covariance budgets,
!  and only for dissipation, since they have no pressure terms.
!-----------------------------------------------------------------------
        subroutine xapxbp_run(ftol, runtype, func_min, c_update)
        use nrtype
        use nr, only: amoeba
        implicit none

!       Input 
        real, intent(in)         :: ftol
        character(*), intent(in) :: runtype

!       Determines whether the module constants are updated
        logical, intent(in)      :: c_update 

!       Internal
        real, dimension(2) :: p_elements 
        real, dimension(2) :: p_variance 

        character(7), dimension(2) :: p_string 

        real opt_error
        integer ndim, mdim
        integer iter, i, j, k

!       Interfaces
        interface
          function func_min(x)
          implicit none
          real, dimension(:), intent(in) :: x
          real func_min
          end function func_min
        end interface
        
!       Setup constants vector
        select case( runtype )
        case("rt")
          p_elements = (/C2rt, mu/)
          p_string   = (/"C2rt   ", "mu     "/)
          p_variance = (/C2rt_var, mu_var/)
        case("thl")
          p_elements = (/C2thl, mu/)
          p_string   = (/"C2thl  ", "mu     "/)
          p_variance = (/C2thl_var, mu_var/)
        case("rtthl")
          p_elements = (/C2rtthl, mu/)
          p_string   = (/"C2rtthl", "mu     "/)
          p_variance = (/C2rtthl_var, mu_var/)
        case default
          stop "xapxbp_run: invalid runtype"
        end select

        ndim = size( p_elements )
        mdim = size( p_elements ) + 1
 
!       Allocate P, Y, and an random factor to generate the starting simplex
        allocate( const_array(mdim, ndim) )
        allocate( diff_array(mdim) )
        allocate( rand_vect(ndim) )

        const_array(1, 1:ndim) = p_elements(1:ndim)

!       Randomize P columns
        do j = 1, ndim
          call random_number( rand_vect(1:ndim) )
          do i = 2, mdim
            const_array(i,j) = const_array(1,j)*
     .                       ( 1.0 - p_variance(j)
     .                         + rand_vect(i-1)*p_variance(j)*2 )
          enddo
        enddo

!       Calculate Y-vector
        do i=1, mdim
          diff_array(i) = func_min( const_array(i,:) )
          if ( .not. isValid ) stop "Initial constants must be valid"
        enddo

        write(*,*) "diff_array:"
        do i=1, mdim
          write(*,'(f12.5)') diff_array(i)
        enddo

!       Attempt to find the optimal P-matrix
        call amoeba( const_array, diff_array, ftol, 
     .               func_min, iter )

!       Write output
        write(*,'(a12,i5)') "Iterations: ", iter
        write(*,'(a12,a12,a12)') "Constant", "Initial", "Optimal"
        do i=1, ndim
          write(*,'(a12,2f12.5)') p_string(i)//"= ", 
     .                            p_elements(i), const_array(1,i)
        enddo


!       Update module constants on c_update
        if ( c_update ) then
          select case( trim(runtype) )
          case("rt")
            C2rt = const_array(1,1)
            mu   = const_array(1,2)
          case("thl")
            C2thl = const_array(1,1)
            mu    = const_array(1,2)
          case("rtthl")
            C2rtthl = const_array(1,1)
            mu      = const_array(1,2)
          end select
        else
!       Generate GRaDS with Optimal constants
          write(*,*) "Running with optimal constants."      
          opt_error = func_min( const_array(1,:) ) 
          write(*,'(a15,f12.5)') "Optimal error=", opt_error
        endif

!       Deallocate & return
        deallocate( const_array, diff_array, rand_vect )

        return
        end subroutine xapxbp_run

!-----------------------------------------------------------------------
!  FUNCTION wp2_min( )
!  Calculate the total difference between the two GRaDS profiles
!-----------------------------------------------------------------------
        real function wp2_min( const_vect )

        use hoc
        use inputfields, only: grads_fields_reader

        implicit none

!       Parameters

        character(7), dimension(6), parameter :: wp2_hoc = 
     .  (/"wp2_pr1", "wp2_pr2", "wp2_pr3", 
     .    "wp2_dp1", "wp2_dp2", "wp2_cl "/) 

        character(6), dimension(2), parameter :: wp2_les = 
     .  (/"wp2_pr", "wp2_dp"/) 

!       Input

        real, dimension(:), intent(in) :: const_vect 


!       The actual subroutine is hoc_model_variable_input, used via overloading

        call hoc_model( const_vect(1), C2rt, C2thl, C2rtthl, C4, 
     .                  const_vect(2), C6rt, C6thl, C7, C8, C10, C11,
     .                  nu1, nu2, nu6, nu8, nu_r, 
     .                  const_vect(3), const_vect(4), ! beta and gamma_coef
     .                  c_K, lmin_coef, taumin, taumax, const_vect(5), ! mu
     .                  run_file, isValid, .false.,
     .                  grads_fields_reader, sample_ratio )

        if ( .not. isValid ) then
          wp2_min = 2 * maxval( diff_array ) - minval( diff_array )
          return
        endif


!       Calculate the error
        wp2_min = term_min( wp2_hoc, wp2_les, "m" )

        print *, "wp2_min=", wp2_min
        return 
        end function wp2_min

!-----------------------------------------------------------------------
!  FUNCTION wp3_min( )
!  Calculate the total difference between the two GRaDS profiles
!-----------------------------------------------------------------------
        real function wp3_min( const_vect )

        use hoc
        use inputfields, only: grads_fields_reader

        implicit none

!       Parameters

        character(7), dimension(4), parameter :: wp3_hoc = 
     .  (/"wp3_pr1", "wp3_pr2", "wp3_dp1", "wp3_cl "/) 

        character(6), dimension(2), parameter :: wp3_les = 
     .  (/"wp3_pr", "wp3_dp"/) 

!       Input

        real, dimension(:), intent(in) :: const_vect 


!       The actual subroutine is hoc_model_variable_input, used via overloading

        call hoc_model( C1, C2rt, C2thl, C2rtthl, C4, C5, 
     .                  C6rt, C6thl, C7, const_vect(1), 
     .                  C10, const_vect(2),
     .                  nu1, nu2, nu6, nu8, nu_r, 
     .                  const_vect(3), const_vect(4),  ! beta and gamma_coef
     .                  c_K, lmin_coef, taumin, taumax, const_vect(5), ! mu
     .                  run_file, isValid, .false.,
     .                  grads_fields_reader, sample_ratio )

        if ( .not. isValid ) then
          wp3_min = 2 * maxval( diff_array ) - minval( diff_array )
          return
        endif


!       Calculate the error
        wp3_min = term_min( wp3_hoc, wp3_les, "t" )

        print *, "wp3_min=", wp3_min

        return 
        end function wp3_min

!-----------------------------------------------------------------------
!  FUNCTION wprtp_min( )
!  Calculate the total difference between the two GRaDS profiles
!-----------------------------------------------------------------------
        real function wprtp_min( const_vect )

        use hoc
        use inputfields, only: grads_fields_reader

        implicit none

!       Parameters

        character(9), dimension(3), parameter :: wprtp_hoc = 
     .  (/"wprtp_pr1", "wprtp_pr2", "wprtp_pr3"/) 

        character(8), dimension(2), parameter :: wprtp_les = 
     .  (/"wpqvp_pr", "wpqcp_pr"/) 

!       Input

        real, dimension(:), intent(in) :: const_vect 


!       The actual subroutine is hoc_model_variable_input, used via overloading

        call hoc_model( C1, C2rt, C2thl, C2rtthl, C4, C5, 
     .                  const_vect(1), C6thl, const_vect(2), 
     .                  C8, C10, C11,
     .                  nu1, nu2, nu6, nu8, nu_r, 
     .                  const_vect(3), const_vect(4), ! beta and gamma_coef
     .                  c_K, lmin_coef, taumin, taumax, mu,
     .                  run_file, isValid, .false.,
     .                  grads_fields_reader, sample_ratio )

        if ( .not. isValid ) then
          wprtp_min = 2 * maxval( diff_array ) - minval( diff_array )
          return
        endif


!       Calculate the error
        wprtp_min = term_min( wprtp_hoc, wprtp_les, "m" )

        print *, "wprtp_min=", wprtp_min

        return 
        end function wprtp_min

!-----------------------------------------------------------------------
!  FUNCTION wpthlp_min( )
!  Calculate the total difference between the two GRaDS profiles
!-----------------------------------------------------------------------
        real function wpthlp_min( const_vect )

        use hoc
        use inputfields, only: grads_fields_reader

        implicit none

!       Parameters

        character(10), dimension(3), parameter :: wpthlp_hoc = 
     .  (/"wpthlp_pr1", "wpthlp_pr2", "wpthlp_pr3"/) 

        character(9), dimension(1), parameter :: wpthlp_les = 
     .  (/"wpthlp_pr"/) 

!       Input

        real, dimension(:), intent(in) :: const_vect 


!       The actual subroutine is hoc_model_variable_input, used via overloading

        call hoc_model( C1, C2rt, C2thl, C2rtthl, C4, C5, 
     .                  C6rt, const_vect(1), const_vect(2), 
     .                  C8, C10, C11,
     .                  nu1, nu2, nu6, nu8, nu_r, 
     .                  const_vect(3), const_vect(4), ! beta and gamma_coef
     .                  c_K, lmin_coef, taumin, taumax, mu,
     .                  run_file, isValid, .false.,
     .                  grads_fields_reader, sample_ratio )

        if ( .not. isValid ) then
          wpthlp_min = 2 * maxval( diff_array ) - minval( diff_array )
          return
        endif


!       Calculate the error
        wpthlp_min = term_min( wpthlp_hoc, wpthlp_les, "m" )

        print *, "wpthlp_min=", wpthlp_min

        return 
        end function wpthlp_min

!-----------------------------------------------------------------------
!  FUNCTION rtp2_min( )
!  Calculate the total difference between the two GRaDS profiles
!  Note that rtp2, thlp2, and rtpthlp have no pressure terms, so this
!  function includes only dissipation.
!-----------------------------------------------------------------------
        real function rtp2_min( const_vect )

        use hoc
        use inputfields, only: grads_fields_reader

        implicit none

!       Parameters

        character(8), dimension(2), parameter :: rtp2_hoc = 
     .  (/"rtp2_dp1", "rtp2_dp2"/)

        character(7), dimension(1), parameter :: rtp2_les = 
     .  (/"qtp2_dp"/) 

!       Input

        real, dimension(:), intent(in) :: const_vect 


!       The actual subroutine is hoc_model_variable_input, used via overloading

        call hoc_model( C1, const_vect(1), C2thl, C2rtthl, C4, ! C2rt
     .                  C5, C6rt, C6thl, C7, C8, C10, C11,
     .                  nu1, nu2, nu6, nu8, nu_r, 
     .                  beta, gamma_coef, c_K, lmin_coef, 
     .                  taumin, taumax, const_vect(2), ! mu
     .                  run_file, isValid, .false.,
     .                  grads_fields_reader, sample_ratio )

        if ( .not. isValid ) then
          rtp2_min = 2 * maxval( diff_array ) - minval( diff_array )
          return
        endif


!       Calculate the error
        rtp2_min = term_min( rtp2_hoc, rtp2_les, "m" )

        print *, "rtp2_min=", rtp2_min

        return 
        end function rtp2_min

!-----------------------------------------------------------------------
!  FUNCTION thlp2_min( )
!  Calculate the total difference between the two GRaDS profiles
!  Note that rtp2, thlp2, and rtpthlp have no pressure terms, so this
!  function includes only dissipation.
!-----------------------------------------------------------------------
        real function thlp2_min( const_vect )

        use hoc
        use inputfields, only: grads_fields_reader

        implicit none

!       Parameters

        character(9), dimension(2), parameter :: thlp2_hoc = 
     .  (/"thlp2_dp1", "thlp2_dp2"/)

        character(8), dimension(1), parameter :: thlp2_les = 
     .  (/"thlp2_dp"/) 

!       Input

        real, dimension(:), intent(in) :: const_vect 


!       The actual subroutine is hoc_model_variable_input, used via overloading

        call hoc_model( C1, C2rt, const_vect(1), C2rtthl, C4, ! C2thl
     .                  C5, C6rt, C6thl, C7, C8, C10, C11,
     .                  nu1, nu2, nu6, nu8, nu_r, 
     .                  beta, gamma_coef, c_K, lmin_coef, 
     .                  taumin, taumax, const_vect(2), ! mu
     .                  run_file, isValid, .false.,
     .                  grads_fields_reader, sample_ratio )

        if ( .not. isValid ) then
          thlp2_min = 2 * maxval( diff_array ) - minval( diff_array )
          return
        endif


!       Calculate the error
        thlp2_min = term_min( thlp2_hoc, thlp2_les, "m" )

        print *, "thlp2_min=", thlp2_min

        return 
        end function thlp2_min

!-----------------------------------------------------------------------
!  FUNCTION rtpthlp_min( )
!  Calculate the total difference between the two GRaDS profiles
!  Note that rtp2, thlp2, and rtpthlp have no pressure terms, so this
!  function includes only dissipation.
!-----------------------------------------------------------------------
        real function rtpthlp_min( const_vect )

        use hoc
        use inputfields, only: grads_fields_reader

        implicit none

!       Parameters

        character(11), dimension(2), parameter :: rtpthlp_hoc = 
     .  (/"rtpthlp_dp1", "rtpthlp_dp2"/)

        character(10), dimension(1), parameter :: rtpthlp_les = 
     .  (/"qtpthlp_dp"/) 

!       Input

        real, dimension(:), intent(in) :: const_vect 


!       The actual subroutine is hoc_model_variable_input, used via overloading

        call hoc_model( C1, C2rt, C2thl, const_vect(1), C4, ! C2rtthl
     .                  C5, C6rt, C6thl, C7, C8, C10, C11,
     .                  nu1, nu2, nu6, nu8, nu_r, 
     .                  beta, gamma_coef, c_K, lmin_coef, 
     .                  taumin, taumax, const_vect(2), ! mu
     .                  run_file, isValid, .false.,
     .                  grads_fields_reader, sample_ratio )

        if ( .not. isValid ) then
          rtpthlp_min = 2 * maxval( diff_array ) - minval( diff_array )
          return
        endif


!       Calculate the error
        rtpthlp_min = term_min( rtpthlp_hoc, rtpthlp_les, "m" )

        print *, "rtpthlp_min=", rtpthlp_min

        return 
        end function rtpthlp_min

!-----------------------------------------------------------------------
!  FUNCTION wp2_buoy_min( )
!  Calculate the mean squared difference between the two GRaDS profiles
!-----------------------------------------------------------------------
        real function wp2_buoy_min( const_vect )

        implicit none

!       Parameters

        character(6), dimension(1), parameter :: wp2_hoc = 
     .  (/"wp2_bp"/) 

        character(6), dimension(1), parameter :: wp2_les = 
     .  (/"wp2_bp"/) 

        character(1), parameter :: grid = "m"

!       Input

        real, dimension(:), intent(in) :: const_vect 

!       Calculate the error
        wp2_buoy_min = gen_betagamma_min( wp2_hoc, wp2_les, 
     .                 const_vect(1), const_vect(2), grid )

        print *, "wp2_buoy_min=", wp2_buoy_min

        return 
        end function wp2_buoy_min

!-----------------------------------------------------------------------
!  FUNCTION wp2_advect_min( )
!  Calculate the mean squared difference between the two GRaDS profiles
!  * Doesn't work -- d( w'3 ) / d ( z ) is computed based on a field
!    that is set from GRaDS data
!-----------------------------------------------------------------------
!        real function wp2_advect_min( const_vect )
!
!        implicit none
!
!       Parameters
!
!        character(6), dimension(3), parameter :: wp2_hoc = 
!     .  (/"wp2_ma", "wp2_ta", "wp2_ac"/) 
!
!        character(6), dimension(1), parameter :: wp2_les = 
!     .  (/"wp2_tp"/) 
!
!        character(1), parameter :: grid = "m"
!
!       Input
!
!        real, dimension(:), intent(in) :: const_vect 
!
!       Calculate the error
!        wp2_advect_min = gen_betagamma_min( wp2_hoc, wp2_les, 
!     .                   const_vect(1), const_vect(2), grid )
!        return 
!        end function wp2_advect_min

!-----------------------------------------------------------------------
!  FUNCTION wp3_buoy_min( )
!  Calculate the mean squared difference between the two GRaDS profiles
!-----------------------------------------------------------------------
        real function wp3_buoy_min( const_vect )

        implicit none

!       Parameters

        character(6), dimension(1), parameter :: wp3_hoc = 
     .  (/"wp3_bp"/) 

        character(6), dimension(1), parameter :: wp3_les = 
     .  (/"wp3_bp"/) 

        character(1), parameter :: grid = "t"

!       Input

        real, dimension(:), intent(in) :: const_vect 

!       Calculate the error
        wp3_buoy_min = gen_betagamma_min( wp3_hoc, wp3_les, 
     .                 const_vect(1), const_vect(2), grid )

        print *, "wp3_buoy_min=", wp3_buoy_min

        return 
        end function wp3_buoy_min

!-----------------------------------------------------------------------
!  FUNCTION wp3_advect_min( )
!  Calculate the mean squared difference between the two GRaDS profiles
!-----------------------------------------------------------------------
        real function wp3_advect_min( const_vect )

        implicit none

!       Parameters

        character(6), dimension(4), parameter :: wp3_hoc = 
     .  (/"wp3_ma", "wp3_ta", "wp3_tp", "wp3_ac"/) 

        character(6), dimension(1), parameter :: wp3_les = 
     .  (/"wp3_tp"/) 

        character(1), parameter :: grid = "t"

!       Input

        real, dimension(:), intent(in) :: const_vect 

!       Calculate the error
        wp3_advect_min = gen_betagamma_min( wp3_hoc, wp3_les, 
     .                   const_vect(1), const_vect(2), grid )

        print *, "wp3_advect_min=", wp3_advect_min

        return 
        end function wp3_advect_min

!-----------------------------------------------------------------------
!  FUNCTION wprtp_buoy_min( )
!  Calculate the mean squared difference between the two GRaDS profiles
!-----------------------------------------------------------------------
        real function wprtp_buoy_min( const_vect )

        implicit none

!       Parameters

        character(8), dimension(1), parameter :: wprtp_hoc = 
     .  (/"wprtp_bp"/) 

        character(8), dimension(2), parameter :: wprtp_les = 
     .  (/"wpqcp_bp", "wpqvp_bp"/) 

        character(1), parameter :: grid = "m"

!       Input

        real, dimension(:), intent(in) :: const_vect 

!       Calculate the error
        wprtp_buoy_min = gen_betagamma_min( wprtp_hoc, wprtp_les, 
     .                   const_vect(1), const_vect(2), grid )

        print *, "wprtp_buoy_min=", wprtp_buoy_min

        return 
        end function wprtp_buoy_min

!-----------------------------------------------------------------------
!  FUNCTION wprtp_advect_min( )
!  Calculate the mean squared difference between the two GRaDS profiles
!-----------------------------------------------------------------------
        real function wprtp_advect_min( const_vect )

        implicit none

!       Parameters

        character(8), dimension(4), parameter :: wprtp_hoc = 
     .  (/"wprtp_ma", "wprtp_ta", "wprtp_tp", "wprtp_ac"/) 

        character(8), dimension(2), parameter :: wprtp_les = 
     .  (/"wpqcp_tp", "wpqvp_tp"/) 

        character(1), parameter :: grid = "m"

!       Input

        real, dimension(:), intent(in) :: const_vect 

!       Calculate the error
        wprtp_advect_min = gen_betagamma_min( wprtp_hoc, wprtp_les, 
     .                     const_vect(1), const_vect(2), grid )

        print *, "wprtp_advect_min=", wprtp_advect_min

        return 
        end function wprtp_advect_min

!-----------------------------------------------------------------------
!  FUNCTION wpthlp_buoy_min( )
!  Calculate the mean squared difference between the two GRaDS profiles
!-----------------------------------------------------------------------
        real function wpthlp_buoy_min( const_vect )

        implicit none

!       Parameters

        character(9), dimension(1), parameter :: wpthlp_hoc = 
     .  (/"wpthlp_bp"/) 

        character(9), dimension(1), parameter :: wpthlp_les = 
     .  (/"wpthlp_bp"/) 

        character(1), parameter :: grid = "m"

!       Input

        real, dimension(:), intent(in) :: const_vect 

!       Calculate the error
        wpthlp_buoy_min = gen_betagamma_min(
     .                    wpthlp_hoc, wpthlp_les, 
     .                    const_vect(1), const_vect(2), grid )

        print *, "wpthlp_buoy_min=", wpthlp_buoy_min

        return 
        end function wpthlp_buoy_min

!-----------------------------------------------------------------------
!  FUNCTION wpthlp_advect_min( )
!  Calculate the mean squared difference between the two GRaDS profiles
!-----------------------------------------------------------------------
        real function wpthlp_advect_min( const_vect )

        implicit none

!       Parameters

        character(9), dimension(4), parameter :: wpthlp_hoc = 
     .  (/"wpthlp_ma", "wpthlp_ta", "wpthlp_tp", "wpthlp_ac"/) 

        character(9), dimension(1), parameter :: wpthlp_les = 
     .  (/"wpthlp_tp"/) 

        character(1), parameter :: grid = "m"

!       Input

        real, dimension(:), intent(in) :: const_vect 

!       Calculate the error
        wpthlp_advect_min = gen_betagamma_min(
     .                      wpthlp_hoc, wpthlp_les, 
     .                      const_vect(1), const_vect(2), grid )

        print *, "wpthlp_advect_min=", wpthlp_advect_min

        return 
        end function wpthlp_advect_min

!-----------------------------------------------------------------------
!  FUNCTION rtp2_advect_min( )
!  Calculate the mean squared difference between the two GRaDS profiles
!-----------------------------------------------------------------------
        real function rtp2_advect_min( const_vect )

        implicit none

!       Parameters

        character(7), dimension(3), parameter :: rtp2_hoc = 
     .  (/"rtp2_ma", "rtp2_ta", "rtp2_tp"/) 

        character(7), dimension(1), parameter :: rtp2_les = 
     .  (/"qtp2_tp"/) 

        character(1), parameter :: grid = "m"

!       Input

        real, dimension(:), intent(in) :: const_vect 

!       Calculate the error

        rtp2_advect_min = gen_betagamma_min(
     .                    rtp2_hoc, rtp2_les, 
     .                    const_vect(1), const_vect(2), grid )

        print *, "rtp2_advect_min=", rtp2_advect_min

        return 
        end function rtp2_advect_min

!-----------------------------------------------------------------------
!  FUNCTION thlp2_advect_min( )
!  Calculate the mean squared difference between the two GRaDS profiles
!-----------------------------------------------------------------------
        real function thlp2_advect_min( const_vect )

        implicit none

!       Parameters

        character(8), dimension(3), parameter :: thlp2_hoc = 
     .  (/"thlp2_ma", "thlp2_ta", "thlp2_tp"/) 

        character(8), dimension(1), parameter :: thlp2_les = 
     .  (/"thlp2_tp"/) 

        character(1), parameter :: grid = "m"

!       Input

        real, dimension(:), intent(in) :: const_vect 

!       Calculate the error

        thlp2_advect_min = gen_betagamma_min(
     .                     thlp2_hoc, thlp2_les, 
     .                     const_vect(1), const_vect(2), grid )

        print *, "thlp2_advect_min=", thlp2_advect_min

        return 
        end function thlp2_advect_min

!-----------------------------------------------------------------------
!  FUNCTION rtpthlp_advect_min( )
!  Calculate the mean squared difference between the two GRaDS profiles
!-----------------------------------------------------------------------
        real function rtpthlp_advect_min( const_vect )

        implicit none

!       Parameters

        character(11), dimension(4), parameter :: rtpthlp_hoc = 
     .  (/"rtpthlp_ma ", "rtpthlp_ta ", "rtpthlp_tp1", "rtpthlp_tp2"/) 

        character(10), dimension(1), parameter :: rtpthlp_les = 
     .  (/"qtpthlp_tp"/) 

        character(1), parameter :: grid = "m"

!       Input

        real, dimension(:), intent(in) :: const_vect 

!       Calculate the error

        rtpthlp_advect_min = gen_betagamma_min(
     .                       rtpthlp_hoc, rtpthlp_les, 
     .                       const_vect(1), const_vect(2), grid )

        print *, "rtpthlp_advect_min=", rtpthlp_advect_min

        return 
        end function rtpthlp_advect_min

!-----------------------------------------------------------------------
!  FUNCTION gen_betagamma_min
!  Calculate the mean squared difference in a term that varies with 
!  respect to beta and gamma
!-----------------------------------------------------------------------
        real function gen_betagamma_min( xvar_hoc, xvar_les, 
     .                                   xbeta, xgamma, grid )

        use hoc
        use inputfields, only: grads_fields_reader

        implicit none

!       Input
        character(len=*), intent(in), dimension(:) :: xvar_hoc
        character(len=*), intent(in), dimension(:) :: xvar_les
        real, intent(in) :: xbeta
        real, intent(in) :: xgamma
        character(len=1), intent(in) :: grid

!       The actual subroutine is hoc_model_variable_input, used via overloading

        call hoc_model( C1, C2rt, C2thl, C2rtthl, C4, C5, 
     .                  C6rt, C6thl, C7, C8, C10, C11,
     .                  nu1, nu2, nu6, nu8, nu_r, 
     .                  xbeta, xgamma,
     .                  c_K, lmin_coef, taumin, taumax, mu,
     .                  run_file, isValid, .false.,
     .                  grads_fields_reader, sample_ratio )

        if ( .not. isValid ) then
          gen_betagamma_min = 2*maxval( diff_array )
     .                           - minval( diff_array )
          return
        endif

        gen_betagamma_min = term_min( xvar_hoc, xvar_les, grid )

        return
        end function gen_betagamma_min

!-----------------------------------------------------------------------
!  FUNCTION term_min
!  Minimize the error between two individual terms in an equation. 
!  Unlike the other procedures this is generalized.
!-----------------------------------------------------------------------
        real function term_min( hoc_variable_names,
     .                          les_variable_names, grid_type )

        use grads_common
        use error, only: mean_sqr_diff_zm, mean_sqr_diff_zt

        implicit none

!       Input

!       These usually are not the same dimension.
!        e.g. wp2_pr1 + wp2_pr2 + wp2_pr3 in HOC (3 elements) 
!        is compared to just wp2_pr in LES.
        character(len=*), dimension(:), intent(in) :: hoc_variable_names
        character(len=*), dimension(:), intent(in) :: les_variable_names

        character(len=1), intent(in) :: grid_type

!       Internal

        integer hoc_nz
        integer les_nz

        integer AllocateStatus
        logical ferror

        real les_minmax
        real, allocatable, dimension(:) :: hoc_zl, hoc_total
        real, allocatable, dimension(:) :: les_zl, les_total
        character(len=50) :: les_stats_file_ext
        character(len=50) :: hoc_stats_file_ext
        integer j, k

!       Allocate the arrays for reading in the grads plot data

        les_stats_file_ext = trim(les_stats_file) // "_sw.ctl"

        if ( grid_type == "m" ) then
          hoc_stats_file_ext = trim(hoc_stats_file) // "_zm.ctl"
        else ! grid_type == "t"
          hoc_stats_file_ext = trim(hoc_stats_file) // "_zt.ctl"
        endif

        hoc_nz = grads_zlvl( hoc_stats_file_ext )
        les_nz = grads_zlvl( les_stats_file_ext )

        allocate( hoc_zl(hoc_nz), hoc_total(hoc_nz), 
     .            les_zl(les_nz), les_total(les_nz), 
     .            stat = AllocateStatus )
        if ( AllocateStatus /= 0 ) then 
          stop "Allocation of arrays in function term_min failed"
        endif

        hoc_total = 0.
        les_total = 0.
        hoc_zl    = 0.
        les_zl    = 0.

!       Read in the values for the specified terms

        do k=1, size( hoc_variable_names )
          hoc_zl = grads_average_interval( hoc_stats_file_ext, hoc_nz, 
     .                                     times, 
     .                                     trim(hoc_variable_names(k)),
     .                                     ferror )
          if (ferror) then
            print *, "grads_average for "//hoc_variable_names(k)
            stop
          endif
          hoc_total = hoc_total + hoc_zl
        enddo

        do k=1, size( les_variable_names )

          select case ( trim( les_variable_names(k) ) )
          case( "wpthlp_pr" )
            les_zl = calc_les_wpthlp( les_nz,
     .                                times, "pr" )
          case( "wpthlp_bp" )
            les_zl = calc_les_wpthlp( les_nz,
     .                                times, "bp" )
          case( "wpthlp_tp" )
            les_zl = calc_les_wpthlp( les_nz,
     .                                times, "tp" )

!         For these cases, the return results has less nz's
          case( "qtp2_tp" )
            les_zl(1:les_nz-1) = calc_les_rtp2_tp( les_nz, times )

          case( "thlp2_tp" )
            les_zl(1:les_nz-1) = calc_les_thlp2_tp( les_nz, times )

          case( "qtpthlp_tp" )
            les_zl(1:les_nz-1) = calc_les_rtpthlp_tp( les_nz, times )

          case( "qtp2_dp" )
            les_zl(1:les_nz-1) = calc_les_rtp2_dp( les_nz, times )

          case( "thlp2_dp" )
            les_zl(1:les_nz-1) = calc_les_thlp2_dp( les_nz, times )

          case( "qtpthlp_dp" )
            les_zl(1:les_nz-1) = calc_les_rtpthlp_dp( les_nz, times )


          case default

            les_zl = grads_average_interval(les_stats_file_ext, les_nz,
     .                                      times, 
     .                                      trim(les_variable_names(k)),
     .                                      ferror )
            if (ferror) then
              print *, "grads_average for "//les_variable_names(k)
              stop
            endif
!          print*, les_variable_names(k), les_zl

          end select

          if (grid_type == "m" ) then
            les_total = les_total + les_zl
          else
            do j=1, les_nz-1
              les_total(j) = les_total(j) 
     .                     + ( les_zl(j) + les_zl(j+1) ) / 2.
            enddo
          endif
        enddo


!       Calculate the difference between the minimum and maximum elements 
!       of LES arrays for normalization purposes

        les_minmax = maxval( les_total(z_i:z_f) ) 
     .             - minval( les_total(z_i:z_f) )
!        les_minmax = maxval( les_total ) - minval( les_total )
        if ( les_minmax == 0. ) stop "les_minmax == 0"

!        print *, "les_minmax", les_minmax
!        print *, "hoc_total", hoc_total
!        print *, "les_total", les_total
!       Return mean squared difference( normalized by les_minmax )

        if ( grid_type == "m" ) then
         term_min = mean_sqr_diff_zm( hoc_total(z_i:z_f), 
     .                                les_total(z_i:z_f), les_minmax )
        else
         term_min = mean_sqr_diff_zt( hoc_total(z_i:z_f), 
     .                                les_total(z_i:z_f), les_minmax )
        endif

        return
        end function term_min

!-----------------------------------------------------------------------
        function calc_les_rtp2_tp( nz, times )

        use grads_common
        use grid_class

        implicit none

!       Input
        integer, intent(in)               :: nz
        integer, dimension(:), intent(in) :: times

!       Return type
        real, dimension(nz) :: calc_les_rtp2_tp

!       Internal

        real, dimension(gr%nnzp) :: qtm
        real, dimension(gr%nnzp) :: wm
        real, dimension(gr%nnzp) :: qtp2
        real, dimension(gr%nnzp) :: wpqcp
        real, dimension(gr%nnzp) :: wpqtp
        real, dimension(gr%nnzp) :: wpqtp2
        real, dimension(gr%nnzp) :: wpqvp

        real, dimension(gr%nnzp) :: qtp2_ma
        real, dimension(gr%nnzp) :: qtp2_ta
        real, dimension(gr%nnzp) :: qtp2_tp

        character(len=50) les_stats_file_sw
        character(len=50) les_stats_file_sm

        logical ferror

!       Read in the needed components
        les_stats_file_sm = trim( les_stats_file )//"_sm.ctl"
        les_stats_file_sw = trim( les_stats_file )//"_sw.ctl"

        qtm(1:gr%nnzp)   = grads_average_interval( les_stats_file_sm,
     .                                             nz, times, "qtm",
     .                                             ferror )

        qtp2(1:gr%nnzp)  = grads_average_interval( les_stats_file_sm,
     .                                             nz, times, "qtp2",
     .                                             ferror )

        wpqvp(1:gr%nnzp) = grads_average_interval( les_stats_file_sw,
     .                                             nz+1, times, 
     .                                             "wpqvp", ferror )

        wpqcp(1:gr%nnzp) = grads_average_interval( les_stats_file_sw,
     .                                             nz+1, times, 
     .                                             "wpqcp", ferror )

        wpqtp(1:gr%nnzp) = grads_average_interval( les_stats_file_sm,
     .                                             nz, times, 
     .                                             "wpqtp", ferror )

        wpqtp2(1:gr%nnzp) = grads_average_interval( les_stats_file_sm,
     .                                             nz, times, 
     .                                             "wpqtp2", ferror )

        wm(1:gr%nnzp)    = grads_average_interval( les_stats_file_sw,
     .                                             nz+1, times, 
     .                                             "wm", ferror )

!       Derive the terms that exist in HOC
        qtp2_tp = -2 * ( wpqvp+wpqcp ) * ddzt( qtm )

        qtp2_ta = -( ddzt( wpqtp2 ) )

        qtp2_ma = -wm * ddzt( qtp2 )

        calc_les_rtp2_tp(1:gr%nnzp) = qtp2_ma + qtp2_ta + qtp2_tp

        return
        end function calc_les_rtp2_tp

!-----------------------------------------------------------------------
        function calc_les_rtpthlp_tp( nz, times )

        use grads_common
        use grid_class

        implicit none

!       Input
        integer, intent(in)               :: nz
        integer, dimension(:), intent(in) :: times

!       Return type
        real, dimension(nz) :: calc_les_rtpthlp_tp

!       Internal

        real, dimension(gr%nnzp) :: qtpthlp
        real, dimension(gr%nnzp) :: wm
        real, dimension(gr%nnzp) :: wpqtpthlp
        real, dimension(gr%nnzp) :: wpqtp
        real, dimension(gr%nnzp) :: thlm
        real, dimension(gr%nnzp) :: qtm
        real, dimension(gr%nnzp) :: wpthlp

        real, dimension(gr%nnzp) :: rtpthlp_ma
        real, dimension(gr%nnzp) :: rtpthlp_ta
        real, dimension(gr%nnzp) :: rtpthlp_tp1
        real, dimension(gr%nnzp) :: rtpthlp_tp2

        character(len=50) les_stats_file_sw
        character(len=50) les_stats_file_sm

        logical ferror

!       Read in the needed components
        les_stats_file_sm = trim( les_stats_file )//"_sm.ctl"
        les_stats_file_sw = trim( les_stats_file )//"_sw.ctl"

        wpthlp(1:gr%nnzp)  = grads_average_interval( les_stats_file_sm,
     .                                               nz, times, 
     .                                               "wpthlp", ferror )

        wpqtp(1:gr%nnzp)   = grads_average_interval( les_stats_file_sm,
     .                                               nz, times, 
     .                                               "wpqtp", ferror )

        wpqtpthlp(1:gr%nnzp) 
     .  = grads_average_interval( les_stats_file_sm,
     .                            nz, times, 
     .                            "wpqtpthlp", ferror )

        thlm(1:gr%nnzp)    = grads_average_interval( les_stats_file_sm,
     .                                               nz, times, 
     .                                               "thlm", ferror )

        qtm(1:gr%nnzp)     = grads_average_interval( les_stats_file_sm,
     .                                               nz, times, 
     .                                               "qtm", ferror )

        qtpthlp(1:gr%nnzp) = grads_average_interval( les_stats_file_sm,
     .                                               nz, times, 
     .                                               "qtpthlp", ferror )

        wm(1:gr%nnzp)      = grads_average_interval( les_stats_file_sw,
     .                                               nz+1, times, 
     .                                               "wm", ferror )
!       Derive the terms that exist in HOC
        rtpthlp_tp1 = -zt2zm( wpqtp ) * ddzt( thlm )

        rtpthlp_tp2 = -zt2zm( wpthlp ) * ddzt( qtm )

        rtpthlp_ta  = -( ddzt( wpqtpthlp ) )

        rtpthlp_ma  = -wm * ddzt( qtpthlp )

        calc_les_rtpthlp_tp(1:gr%nnzp) = rtpthlp_ma + rtpthlp_ta 
     .                                 + rtpthlp_tp1 +  rtpthlp_tp2 

        return
        end function calc_les_rtpthlp_tp

!-----------------------------------------------------------------------
        function calc_les_thlp2_tp( nz, times )

        use grads_common
        use grid_class

        implicit none

!       Input
        integer, intent(in)               :: nz
        integer, dimension(:), intent(in) :: times

!       Return type
        real, dimension(nz) :: calc_les_thlp2_tp

!       Internal

        real, dimension(gr%nnzp) :: thlm
        real, dimension(gr%nnzp) :: wm
        real, dimension(gr%nnzp) :: wpthlp
        real, dimension(gr%nnzp) :: wpthlp2
        real, dimension(gr%nnzp) :: thlp2

        real, dimension(gr%nnzp) :: thlp2_ma
        real, dimension(gr%nnzp) :: thlp2_ta
        real, dimension(gr%nnzp) :: thlp2_tp

        character(len=50) les_stats_file_sw
        character(len=50) les_stats_file_sm

        logical ferror

!       Read in the needed components
        les_stats_file_sm = trim( les_stats_file )//"_sm.ctl"
        les_stats_file_sw = trim( les_stats_file )//"_sw.ctl"

        wpthlp(1:gr%nnzp) = grads_average_interval( les_stats_file_sm,
     .                                              nz, times, 
     .                                              "wpthlp", ferror )

        wpthlp2(1:gr%nnzp) = grads_average_interval( les_stats_file_sm,
     .                                              nz, times, 
     .                                              "wpthlp2", ferror )

        thlm(1:gr%nnzp)   = grads_average_interval( les_stats_file_sm,
     .                                              nz, times, 
     .                                              "thlm", ferror )

        thlp2(1:gr%nnzp)  = grads_average_interval( les_stats_file_sm,
     .                                              nz, times, 
     .                                              "thlp2", ferror )

        wm(1:gr%nnzp)     = grads_average_interval( les_stats_file_sw,
     .                                              nz+1, times, 
     .                                              "wm", ferror )
!       Derive the terms that exist in HOC
        thlp2_tp = -2 * zt2zm( wpthlp ) * ddzt( thlm )

        thlp2_ta = -( ddzt( wpthlp2 ) )

        thlp2_ma = -wm * ddzt( thlp2 )

        calc_les_thlp2_tp(1:gr%nnzp) = thlp2_ma + thlp2_ta + thlp2_tp

        return
        end function calc_les_thlp2_tp

!-----------------------------------------------------------------------
        function calc_les_wpthlp( nz, times, term )

        use grads_common

        implicit none

!       Input
        integer, intent(in)               :: nz
        integer, dimension(:), intent(in) :: times
        character(len=*), intent(in)      :: term

!       Return type
        real, dimension(nz) :: calc_les_wpthlp

!       Internal
        real, dimension(nz) :: wpthlp
        real, dimension(nz) :: wpthp
        real, dimension(nz) :: wpqcp

!        real, dimension(nz-1) :: p_tmp   ! atmos press in Pa (t grid)
        character(len=50) les_stats_file_ext

        logical ferror

        les_stats_file_ext = trim( les_stats_file )//"_sw.ctl"

        wpthp = grads_average_interval( les_stats_file_ext,
     .                                  nz, times, "wpthp_"//term,
     .                                  ferror )

        wpqcp = grads_average_interval( les_stats_file_ext,
     .                                  nz, times, "wpqcp_"//term,
     .                                  ferror )
                
        wpthlp = wpthp 
     .         + wpqcp * theta_on_T( nz, times, les_stats_file_ext )

        calc_les_wpthlp(1:nz) = wpthlp(1:nz)

        return
        end function calc_les_wpthlp

!-----------------------------------------------------------------------
        function calc_les_rtp2_dp( nz, times )

        use grads_common
        use grid_class

        implicit none

!       Input
        integer, intent(in)               :: nz
        integer, dimension(:), intent(in) :: times

!       Return type
        real, dimension(nz) :: calc_les_rtp2_dp

!       Internal
        real, dimension(gr%nnzp) :: qvpqcp_bt
        real, dimension(gr%nnzp) :: qcp2_bt
        real, dimension(gr%nnzp) :: qvp2_bt

        real, dimension(gr%nnzp) :: qtp2_bt

        character(len=50) les_stats_file_sm

        logical ferror

!       Read in the needed components
        les_stats_file_sm = trim( les_stats_file )//"_sm.ctl"

        qvpqcp_bt(1:gr%nnzp) = 
     .    grads_average_interval( les_stats_file_sm,
     .                            nz, times, "qvpqcp_bt", ferror )

        qvp2_bt(1:gr%nnzp) = 
     .    grads_average_interval( les_stats_file_sm,
     .                            nz, times, "qvp2_bt", ferror )

        qcp2_bt(1:gr%nnzp) = 
     .    grads_average_interval( les_stats_file_sm,
     .                            nz, times, "qcp2_bt", ferror )

        qtp2_bt(1:gr%nnzp) = qcp2_bt 
     .                     + 2 * qvpqcp_bt
     .                     + qvp2_bt

        calc_les_rtp2_dp(1:gr%nnzp) = zt2zm( qtp2_bt )
     .                              - calc_les_rtp2_tp( nz, times )

        return
        end function calc_les_rtp2_dp

!-----------------------------------------------------------------------
        function calc_les_thlp2_dp( nz, times )

        use grads_common
        use grid_class

        implicit none

!       Input
        integer, intent(in)               :: nz
        integer, dimension(:), intent(in) :: times

!       Return type
        real, dimension(nz) :: calc_les_thlp2_dp

!       Internal
        real, dimension(gr%nnzp) :: thp2_bt
        real, dimension(gr%nnzp) :: thpqcp_bt
        real, dimension(gr%nnzp) :: qcp2_bt
        real, dimension(gr%nnzp) :: th_on_t

        real, dimension(gr%nnzp) :: thlp2_bt

        character(len=50) les_stats_file_sm
        character(len=50) les_stats_file_sw

        logical ferror

!       Read in the needed components
        les_stats_file_sm = trim( les_stats_file )//"_sm.ctl"
        les_stats_file_sw = trim( les_stats_file )//"_sw.ctl"

        thp2_bt(1:gr%nnzp) = 
     .    grads_average_interval( les_stats_file_sm,
     .                            nz, times, "thp2_bt", ferror )

!        thpqcp_bt(1:gr%nnzp) = 
!     .    grads_average_interval( les_stats_file_sm,
!     .                            nz, times, "thp2_bt", ferror )

        thpqcp_bt(1:gr%nnzp) = 
     .    grads_average_interval( les_stats_file_sm,
     .                            nz, times, "thpqcp_bt", ferror )

        qcp2_bt(1:gr%nnzp) = 
     .    grads_average_interval( les_stats_file_sm,
     .                            nz, times, "qcp2_bt", ferror )

        th_on_T(1:gr%nnzp) = 
     .    zm2zt( theta_on_T( nz, times, les_stats_file_sw ) )

        thlp2_bt(1:gr%nnzp) = thp2_bt 
     .                      + 2 * thpqcp_bt * th_on_t 
     .                      + qcp2_bt * th_on_t**2 

        calc_les_thlp2_dp(1:gr%nnzp) = zt2zm( thlp2_bt )
     .                               - calc_les_thlp2_tp( nz, times )

!        print '(e25.18)', calc_les_thlp2_dp

        return
        end function calc_les_thlp2_dp

!-----------------------------------------------------------------------
        function calc_les_rtpthlp_dp( nz, times )

        use grads_common
        use grid_class

        implicit none

!       Input
        integer, intent(in)               :: nz
        integer, dimension(:), intent(in) :: times

!       Return type
        real, dimension(nz) :: calc_les_rtpthlp_dp

!       Internal
        real, dimension(gr%nnzp) :: qvpqcp_bt
        real, dimension(gr%nnzp) :: qcp2_bt
        real, dimension(gr%nnzp) :: th_on_t
        real, dimension(gr%nnzp) :: thpqvp_bt
        real, dimension(gr%nnzp) :: thpqcp_bt

        real, dimension(gr%nnzp) :: qtpthlp_bt

        character(len=50) les_stats_file_sm
        character(len=50) les_stats_file_sw

        logical ferror

!       Read in the needed components
        les_stats_file_sm = trim( les_stats_file )//"_sm.ctl"
        les_stats_file_sw = trim( les_stats_file )//"_sw.ctl"

        qvpqcp_bt(1:gr%nnzp) = 
     .    grads_average_interval( les_stats_file_sm,
     .                            nz, times, "qvpqcp_bt", ferror )

        thpqcp_bt(1:gr%nnzp) = 
     .    grads_average_interval( les_stats_file_sm,
     .                            nz, times, "thpqcp_bt", ferror )
        thpqvp_bt(1:gr%nnzp) = 
     .    grads_average_interval( les_stats_file_sm,
     .                            nz, times, "thpqvp_bt", ferror )

        qcp2_bt(1:gr%nnzp) = 
     .    grads_average_interval( les_stats_file_sm,
     .                            nz, times, "qcp2_bt", ferror )

        th_on_T(1:gr%nnzp) = 
     .    zm2zt( theta_on_T( nz, times, les_stats_file_sw ) )

        qtpthlp_bt(1:gr%nnzp) = thpqcp_bt 
     .                        + thpqvp_bt
     .                        + th_on_T * qcp2_bt
     .                        + th_on_T * qvpqcp_bt

        calc_les_rtpthlp_dp(1:gr%nnzp) =
     .    zt2zm( qtpthlp_bt ) - calc_les_rtpthlp_tp( nz, times )

!        print '(e25.18)', calc_les_rtpthlp_dp

        return
        end function calc_les_rtpthlp_dp

!-----------------------------------------------------------------------
!  FUNCTION theta_on_T()
!  This equation is used often and so was made its own function.
!  Note that the R exponent used in earlier revisions of this code was
!  incorrect and the exner data in the LES already included it.
!-----------------------------------------------------------------------
        function theta_on_T( nz, times, les_stats_file_ext )
        use grads_common
        use constants, only: Lv, Cp

        implicit none

!       Parameters (included from constants.mod)
!       Lv = latent heat of vaporization at STP
!       Cp = Specific heat of dry air at constant pressure

!       Input
        integer, intent(in)               :: nz
        integer, dimension(:), intent(in) :: times
        character(len=*), intent(in)      :: les_stats_file_ext

!       Return type
        real, dimension(nz) :: theta_on_T

!       Internal
        real, dimension(nz) :: exner
        logical ferror
        
        exner = grads_average_interval( les_stats_file_ext,
     .                                  nz, times, 
     .                                  "ex0", ferror )

        theta_on_T(1:nz) = (Lv/Cp) * (1.0/exner)

        end function theta_on_T

!-----------------------------------------------------------------------
        subroutine output_optimal_constants( )
        use hoc
        use inputfields, only: grads_fields_reader

        implicit none

        real, dimension(24)          :: p_elements
        character(10), dimension(24) :: p_string

        integer i

        p_elements = (/C1, C2rt, C2thl, C2rtthl, C4, C5, 
     .                 C6rt, C6thl, C7, C8, C10, C11,
     .                 nu1, nu2, nu6, nu8, nu_r, 
     .                 beta, gamma_coef, c_K, 
     .                 lmin_coef, taumin, taumax, mu/)

        p_string = (/"C1        ", "C2rt      ", "C2thl     ", 
     .               "C2rtthl   ", "C4        ", "C5        ", 
     .               "C6rt      ", "C6thl     ", "C7        ", 
     .               "C8        ", "C10       ", "C11       ", 
     .               "nu1       ", "nu2       ", "nu6       ", 
     .               "nu8       ", "nu_r      ", "beta      ", 
     .               "gamma_coef", "c_K       ", "lmin_coef ", 
     .               "taumin    ", "taumax    ", "mu        "/)
 
        write(*,*) "Calling HOC with optimal constants"

!       The actual subroutine is hoc_model_variable_input, used via overloading
        call hoc_model( C1, C2rt, C2thl, C2rtthl, C4, C5, 
     .                  C6rt, C6thl, C7, C8, C10, C11,
     .                  nu1, nu2, nu6, nu8, nu_r, 
     .                  beta, gamma_coef,
     .                  c_K, lmin_coef, taumin, taumax, mu,
     .                  run_file, isValid, .false.,
     .                  grads_fields_reader, sample_ratio )

        do i=1, size( p_elements )
          write(*,'(a10,f12.6)') p_string(i), p_elements(i)
        enddo

        return
        end subroutine output_optimal_constants

!-----------------------------------------------------------------------
!  SUBROUTINE generate_namelists
!  Writes the namelists to the parameterized file, usually at the end of
!  a tuning run.
!-----------------------------------------------------------------------
        subroutine generate_namelists( filename, ftol )

        use inputfields

        implicit none

!       Input

        character(len=*), intent(in) :: filename

        real, intent(in)             :: ftol

!       Internal

        integer ios

!       Namelists

        namelist /budget_tune/ run_file, hoc_stats_file, les_stats_file,
     .                         times_nl, z_i, z_f, ftol

        namelist /setfields/ datafile, sample_ratio, input_type,
     .                       input_um, input_vm, input_rtm, input_thlm,
     .                       input_wp2, input_wprtp, input_wpthlp, 
     .                       input_wp3, input_rtp2, input_thlp2, 
     .                       input_rtpthlp, input_upwp, input_vpwp


        open( 20, file=filename, status='unknown', 
     .        iostat=ios, action='write' )
        if ( ios /= 0 ) then
          write(*,*) "Error writing to " // filename
          return
        endif

        write(20,*) "! Namelist generated by hoc_tuner_budget_terms "
        if ( lwp2_prdp ) write(20,*) "! tuned for wp2_pr & wp2_dp"
        if ( lwp3_prdp ) write(20,*) "! tuned for wp3_pr & wp3_dp"
        if ( lwprtp_prdp ) write(20,*) "! tuned for wprtp_pr & wprtp_dp"
        if ( lwpthlp_prdp ) write(20,*) 
     .                      "! tuned for wpthlp_pr & wpthlp_dp"
        if ( lwp2_bp ) write(20,*) "! tuned for wp2_bp "
        if ( lwp3_bp ) write(20,*) "! tuned for wp3_bp "
        if ( lwprtp_bp ) write(20,*) "! tuned for wprtp_bp "
        if ( lwpthlp_bp ) write(20,*) "! tuned for wpthlp_bp "
        if ( lwp2_tp ) write(20,*) "! tuned for wp2_tp "
        if ( lwp3_tp ) write(20,*) "! tuned for wp3_tp "
        if ( lwprtp_tp ) write(20,*) "! tuned for wprtp_tp "
        if ( lwpthlp_tp ) write(20,*) "! tuned for wpthlp_tp "
 

        write(20, nml=budget_tune)
        write(20, nml=initvars)
        write(20, nml=variance)
        write(20, nml=setfields)
        write(20, nml=case_tune)

        close(20)

        return
        end subroutine generate_namelists

      end module budget_terms
