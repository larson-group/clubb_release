!-----------------------------------------------------------------------
! $Id: budget_terms.F,v 1.2 2005-08-05 18:03:46 dschanen Exp $     
! MODULE budget_terms
! functions and subroutines for tuning for budget terms in the HOC model.
!-----------------------------------------------------------------------
      module budget_terms

        implicit none

!       Parameters (public)

        integer, parameter :: ndim = 23        ! dimensions of constants array
        integer, parameter :: mdim = ndim + 1  ! dimensions of constants array

!       Run Information

        character(len=50), private :: run_file       ! namelist
        character(len=50), private :: hoc_stats_file ! GRaDS data
        character(len=50), private :: les_stats_file ! GRaDS data

        real, private :: sample_ratio

        integer, allocatable, dimension(:), private :: times

!       diff_array is set from hoc_tuner_budget_terms, it is needed 
!       by budget_min if isValid == .false.

        real, dimension(mdim) :: diff_array ! 'y' vector in amoeba

      
        contains
!-----------------------------------------------------------------------
!  SUBROUTINE setup_budget_terms
!  Configures information about the data files and run required for the
!  other functions found in the module.
!-----------------------------------------------------------------------
        subroutine setup_budget_terms( run_file_in, sample_ratio_in,
     .                                 hoc_stats_file_in, 
     .                                 les_stats_file_in, times_in )

        implicit none

!       Input

        real, intent(in) :: sample_ratio_in

        character(len=*), intent(in) :: run_file_in
        character(len=*), intent(in) :: hoc_stats_file_in
        character(len=*), intent(in) :: les_stats_file_in

        integer, dimension(:), intent(in) :: times_in

!       Initialize module information

        run_file       = run_file_in
        sample_ratio   = sample_ratio_in
        hoc_stats_file = hoc_stats_file_in
        les_stats_file = les_stats_file_in

        allocate ( times(size( times_in )) )
        times = times_in

        return
        end subroutine setup_budget_terms
!-----------------------------------------------------------------------
!  FUNCTION budget_min( )
!  Attempt to minimize the error between HOC and LES of all 
!  parameterized variables.
!-----------------------------------------------------------------------
        real function budget_min( const_vect )

        use hoc
        use inputfields, only: grads_fields_reader

        implicit none

!       Parameters

        character(7), dimension(3), parameter :: wp2_hoc = 
     .  (/"wp2_pr1", "wp2_pr2", "wp2_pr3"/) 

        character(6), dimension(1), parameter :: wp2_les = 
     .  (/"wp2_pr"/) 

        character(7), dimension(2), parameter :: wp3_hoc = 
     .  (/"wp3_pr1", "wp3_pr2"/) 

        character(6), dimension(1), parameter :: wp3_les = 
     .  (/"wp3_pr"/) 

!       Input

        real, dimension(:), intent(in) :: const_vect  ! ndim dimension

!       Internal        

        logical isValid

!       The actual subroutine is hoc_model_variable_input, used via overloading

        call hoc_model( const_vect(1), const_vect(2), const_vect(3),
     .                  const_vect(4), const_vect(5), const_vect(6),
     .                  const_vect(7), const_vect(8), const_vect(9),
     .                  const_vect(10), const_vect(11), const_vect(12),
     .                  const_vect(13), const_vect(14), const_vect(15),
     .                  const_vect(16), const_vect(17), const_vect(18),
     .                  const_vect(19), const_vect(20), const_vect(21),
     .                  const_vect(22), const_vect(23),
     .                  run_file, isValid, .false.,
     .                  grads_fields_reader, sample_ratio )

        if ( .not. isValid ) then
          budget_min = 2 * maxval( diff_array ) - minval( diff_array )
          return
        endif

!       Calculate the error
        budget_min = term_min( wp2_hoc, wp2_les, "m" )
        budget_min = budget_min + term_min( wp3_hoc, wp3_les, "t" )
!        budget_min = budget_min + term_min( rtpthlp_hoc, rtpthlp_les )

!       Divide by the number of calls to term_min()
        budget_min = budget_min / 2.

!        print *,  budget_min   ! %% debug

        return 
        end function budget_min
!-----------------------------------------------------------------------
!  FUNCTION term_min
!  Minimize the error between two individual terms in an equation. 
!  
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

!       Read in the values for the pressure terms

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
          les_zl = grads_average_interval( les_stats_file_ext, les_nz,
     .                                     times, 
     .                                     trim(les_variable_names(k)),
     .                                     ferror )
          if (ferror) then
            print *, "grads_average for "//les_variable_names(k)
            stop
          endif
          if (grid_type == "m" ) then
            les_total = les_total + les_zl
          else
            do j=1, les_nz-1
              les_total(j) = ( les_zl(j) + les_zl(j+1) ) / 2.
            enddo
          endif
        enddo


!       Calculate the difference between the minimum and maximum elements 
!       of LES arrays for normalization purposes

        les_minmax = maxval( les_total ) - minval( les_total )
        if ( les_minmax == 0. ) stop "les_minmax == 0"

!       Return mean squared difference( normalized by les_minmax )
        if ( grid_type == "m" ) then
         term_min = mean_sqr_diff_zm( hoc_total, les_total, les_minmax )
        else
         term_min = mean_sqr_diff_zt( hoc_total, les_total, les_minmax )
        endif

        return
        end function term_min
!-----------------------------------------------------------------------
      end module budget_terms
