!=======================================================================
!  $Id: hoc_tuner_budget_terms.F,v 1.1 2005-08-04 21:25:54 dschanen Exp $
! PROGRAM hoc_tuner_budget_terms

! Attempts to tune the constants of the HOC model to best fit the
! individual budget terms of prognostic variables.

! Currently only setup for wp2 and wp3
!=======================================================================

      program hoc_tuner_budget_terms
      use budget_terms
      use error, only: read_random_seed
      use inputfields
      use nr, only: amoeba
      use nrtype

      implicit none

!     Internal
      character(len=50) run_file
      character(len=50) hoc_stats_file
      character(len=50) les_stats_file

      real ftol
      real sample_ratio

      real C1, C2rt, C2thl, C2rtthl, C4, C5, C6rt, C6thl, 
     .     C7, C8, C10, C11
      real nu1, nu2, nu6, nu8
      real beta, gamma_coef, c_K, lmin_coef, taumin, taumax, mu

      real C4_var, C5_var, C8_var, C11_var

      real, dimension(ndim)       :: init_const
      real, dimension(mdim)       :: rand_vect
      real, dimension(mdim, ndim) :: const_array
!      real, dimension(mdim)       :: diff_array
      real, dimension(mdim)       :: const_tmp

      integer, dimension(10)      :: times

      real opt_error

      integer j, k  ! loop indices
      integer iter  ! used by amoeba

!     Namelists

      namelist /initvars/ C1, C2rt, C2thl, C2rtthl, C4, C5, 
     .                    C6rt, C6thl, C7, C8, C10, C11, 
     .                    nu1, nu2, nu6, nu8, beta, gamma_coef,
     .                    c_K, lmin_coef, taumin, taumax, mu

      namelist /variance/ C4_var, C5_var, C8_var, C11_var 

      namelist /budget_tune/ run_file, hoc_stats_file, les_stats_file, 
     .                       times, ftol

      namelist /setfields/ datafile, sample_ratio, input_type,
     .                     input_um, input_vm, input_rtm, input_thlm,
     .                     input_wp2, input_wprtp, input_wpthlp, 
     .                     input_wp3, input_rtp2, input_thlp2, 
     .                     input_rtpthlp, input_upwp, input_vpwp

!     Initialize variables
      times = 0

!     Read in the namelists
      open(10, name='budget.in', status='old')
      read(10, nml=budget_tune)
      read(10, nml=initvars)
      read(10, nml=variance)
      read(10, nml=setfields)
      close(10)

      call setup_budget_terms( run_file, sample_ratio, 
     .                         hoc_stats_file, les_stats_file, 
     .                         times(1:maxloc(times, 1)) )
      

!     Create a vectorized version of the constants
      init_const(1:ndim) = (/C1, C2rt, C2thl, C2rtthl, C4, C5, 
     .                       C6rt, C6thl, C7, C8, C10, C11, 
     .                       nu1, nu2, nu6, nu8, beta, gamma_coef, 
     .                       c_K, lmin_coef, taumin, taumax, mu/)

!     Initial constants to all columns of the 'p' matrix
      const_array(1:mdim,:) = spread( init_const, 1, mdim )

!     Generate a random vector
      call read_random_seed("../tune/rand_seed.dat") ! from error.mod
      call random_number( rand_vect(1:mdim) )

!     Vary the constants that influence wp2 and wp3
      const_array(1:mdim, 5)  = const_array(1:mdim, 5) * (1.0 - C4_var
     .                        + rand_vect(1:mdim) * (2.0*C4_var) )

      const_array(1:mdim, 6)  = const_array(1:mdim, 6) * (1.0 - C5_var
     .                        + rand_vect(1:mdim) * (2.0*C5_var) )

      const_array(1:mdim, 10) = const_array(1:mdim, 10) * (1.0 - C8_var
     .                        + rand_vect(1:mdim) * (2.0*C8_var) )

      const_array(1:mdim, 12) = const_array(1:mdim, 12) * (1.0 - C11_var
     .                        + rand_vect(1:mdim) * (2.0*C11_var) )

!     Calculate initial 'y' vector
      do k=1, mdim
        diff_array(k) = budget_min( const_array(k, 1:ndim) )
      enddo

      write(*,*) "diff_array:"
      write(*,'6f12.5') diff_array

      call amoeba( const_array, diff_array, ftol, 
     .             budget_min, iter )

      print '(a12,i5)', "Iterations: ", iter
      print '(a12,a12,a12)', "Constant", "Initial", "Optimal"
      print '(a12,2f12.5)', "C4=", C4, const_array(1,5)
      print '(a12,2f12.5)', "C5=", C5, const_array(1,6)
      print '(a12,2f12.5)', "C8=", C8, const_array(1,10)
      print '(a12,2f12.5)', "C11=", C11, const_array(1,12)

      print *, "Running with optimal constants."      
      opt_error = budget_min( const_array(1,:) ) 
      print '(a15,f12.5)', "Optimal error=", opt_error

      end program hoc_tuner_budget_terms
!-----------------------------------------------------------------------
