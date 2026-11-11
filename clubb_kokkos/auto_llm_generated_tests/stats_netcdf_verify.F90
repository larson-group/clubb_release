! Read closed stats-test files and check native averaging, budgets and batches.
! native_verify_stats_output_c returns zero when both independent contexts agree
! with the exact expected fields and metadata; a nonzero value fails the CTest.
module stats_netcdf_test_verify
  use, intrinsic :: iso_c_binding, only: c_int, c_double
  use netcdf
  implicit none
contains
  integer(c_int) function native_verify_stats_output_c() bind(C) result(result)
    ! Description: Check both closed files; return zero only when both pass.
    result = verify_file('stats_netcdf_a.nc', 4, .true.)
    if ( result == 0 ) result = verify_file('stats_netcdf_b.nc', 2, .false.)
  end function native_verify_stats_output_c

  integer function verify_file(path, ncol, batched) result(result)
    ! Description: Check exact means, budget tendencies and column placement in
    ! one file. Return zero for success, one for an I/O or value mismatch.

    !--------------------------- Input Variables ---------------------------
    character(len=*), intent(in) :: path ! Closed NetCDF file to inspect
    integer, intent(in) :: ncol ! Number of output columns [-]
    logical, intent(in) :: batched ! Whether the file combines two runtime batches

    !--------------------------- Local Variables ---------------------------
    integer :: ncid ! Open NetCDF file identifier [-]
    integer :: varid ! Current NetCDF variable identifier [-]
    integer :: ierr ! NetCDF/read assertion status [-]
    integer :: column, level, sample ! Output column and level indices [-]
    integer :: base ! Column-dependent offset in the expected test data [K]
    integer :: close_error ! NetCDF close status [-]
    real(c_double) :: field(ncol,2,1) ! Mean profile [K]
    real(c_double) :: substep(ncol,1) ! Mean of two substeps [K]
    real(c_double) :: surface(ncol,1) ! Mean surface field [K]
    real(c_double) :: params(ncol,1) ! Saved test parameter [-]
    real(c_double) :: times(1) ! Output time [s]
    real(c_double) :: budget(ncol,2,1) ! Mean closing-minus-opening tendency [K/s]

    real(c_double) :: sample_field(2,2,2,1) ! SILHS field (column,sample,level,time) [-]
    real(c_double) :: sample_zt(2) ! SILHS vertical coordinates [m]
    integer :: sample_var ! Sample output selector [-]
    character(len=20), parameter :: sample_names(4) = &
      [character(len=20) :: 'lh_nl_scalar', 'lh_u_scalar', 'lh_u_component', 'lh_u_weight']

    !----------------------------- Begin Code ------------------------------
    result = 1
    ierr = nf90_open(path,NF90_NOWRITE,ncid)
    if ( ierr /= NF90_NOERR ) return
    ierr = nf90_inq_varid(ncid,'field',varid)
    if ( ierr == NF90_NOERR ) ierr = nf90_get_var(ncid,varid,field)
    if ( ierr == NF90_NOERR ) ierr = nf90_inq_varid(ncid,'surface',varid)
    if ( ierr == NF90_NOERR ) ierr = nf90_get_var(ncid,varid,surface)
    if ( ierr == NF90_NOERR ) ierr = nf90_inq_varid(ncid,'time',varid)
    if ( ierr == NF90_NOERR ) ierr = nf90_get_var(ncid,varid,times)

    if ( batched .and. ierr == NF90_NOERR ) then
      ierr = nf90_inq_varid(ncid,'substep',varid)
      if ( ierr == NF90_NOERR ) ierr = nf90_get_var(ncid,varid,substep)
      if ( ierr == NF90_NOERR ) then
        if ( any(substep /= 4._c_double) ) ierr = -1
      end if
    end if

    if ( batched .and. ierr == NF90_NOERR ) then
      ierr = nf90_inq_varid(ncid,'budget',varid)
      if ( ierr == NF90_NOERR ) ierr = nf90_get_var(ncid,varid,budget)
      if ( ierr == NF90_NOERR ) then
        if ( any(budget /= 2._c_double) ) ierr = -1
      end if
    end if

    if ( batched .and. ierr == NF90_NOERR ) then
      ierr = nf90_inq_varid(ncid,'clubb_params',varid)
      if ( ierr == NF90_NOERR ) ierr = nf90_get_var(ncid,varid,params)
      if ( ierr == NF90_NOERR ) then
        if ( any(params /= 7._c_double) ) ierr = -1
      end if
    end if

    if ( .not. batched .and. ierr == NF90_NOERR ) then
      ierr = nf90_inq_varid(ncid,'lh_zt',varid)
      if ( ierr == NF90_NOERR ) ierr = nf90_get_var(ncid,varid,sample_zt)
      if ( ierr == NF90_NOERR ) then
        if ( any(sample_zt /= [50._c_double,150._c_double]) ) ierr = -1
      end if
      do sample_var = 1, 4
        if ( ierr /= NF90_NOERR ) exit
        ierr = nf90_inq_varid(ncid,trim(sample_names(sample_var)),varid)
        if ( ierr == NF90_NOERR ) ierr = nf90_get_var(ncid,varid,sample_field)
        if ( ierr /= NF90_NOERR ) exit
        do level = 1, 2
          do sample = 1, 2
            do column = 1, 2
              base = 20 + column-1 + 2*(sample-1) + 4*(level-1)
              select case (sample_var)
              case (2)
                base = base + 20
              case (3)
                base = column
              case (4)
                base = base + 40
              end select
              if ( sample_field(column,sample,level,1) /= real(base,c_double) ) ierr = -1
            end do
          end do
        end do
      end do
    end if

    close_error = nf90_close(ncid)
    if ( ierr /= NF90_NOERR .or. close_error /= NF90_NOERR ) return
    if ( times(1) /= 2._c_double ) return

    do column = 1, ncol
      if ( batched ) then
        base = ((column-1)/2)*100 + mod(column-1,2) + 1
      else
        base = 200 + column
      end if
      if ( surface(column,1) /= real(base,c_double) ) return
      do level = 1, 2
        if ( field(column,level,1) /= real(base+10*(level-1),c_double) ) return
      end do
    end do

    result = 0
  end function verify_file
end module stats_netcdf_test_verify
