!----------------------------------------------------------------------
! $Id$

program clubb_standalone_loss

! Description:
!   Standalone front-end for the in-memory CLUBB loss driver.
!-----------------------------------------------------------------------

  use clubb_loss_driver, only: &
    loss_name_len, &
    clubb_get_loss

  use err_info_type_module, only: &
    err_info_type, &
    cleanup_err_info_api

  use clubb_precision, only: &
    core_rknd

  implicit none

  integer, parameter :: &
    success_code = 6

  integer :: &
    arg_len, &
    arg_status

  character(len=256) :: &
    namelist_filename

  real( kind = core_rknd ), allocatable, dimension(:,:,:) :: &
      rmse, &
      correlation, &
      std_ratio, &
      bias_norm
  real( kind = core_rknd ), allocatable, dimension(:) :: &
      case_loss

  character(len=loss_name_len), allocatable, dimension(:) :: &
    clubb_var_names

  type(err_info_type) :: &
    err_info

  character(len=256) :: arg

  namelist_filename = "clubb.in"
  if ( command_argument_count() >= 1 ) then
    call get_command_argument( 1, arg, length = arg_len, status = arg_status )
    if ( arg_status == 0 .and. arg_len > 0 ) then
      namelist_filename = trim( arg(1:arg_len) )
    end if
  end if

  call clubb_get_loss( trim( namelist_filename ), clubb_var_names, rmse, correlation, std_ratio, &
                       bias_norm, case_loss, err_info )

  call print_loss_matrix( clubb_var_names, rmse, correlation, std_ratio, bias_norm, case_loss )

  call cleanup_err_info_api( err_info )
  call exit( success_code )

contains

  subroutine print_loss_matrix( clubb_var_names, rmse, correlation, std_ratio, bias_norm, case_loss )

    ! Description:
  !   Print one row per variable/window with the explicit loss-driver metrics
  !   and best-performing parameter column for each variable.

    implicit none

    character(len=*), dimension(:), intent(in) :: &
      clubb_var_names

    real( kind = core_rknd ), dimension(:,:,:), intent(in) :: &
      rmse, &
      correlation, &
      std_ratio, &
      bias_norm
    real( kind = core_rknd ), dimension(:), intent(in) :: &
      case_loss

    integer :: &
      row_idx, &
      col_idx, &
      window_idx, &
      best_col_idx

    real( kind = core_rknd ) :: &
      best_col_loss

    write( *, '(a)', advance = "no" ) "variable"
    if ( size( rmse, 1 ) > 1 ) write( *, '(1x,a)', advance = "no" ) "window"
    do col_idx = 1, size( rmse, 3 )
      write( *, '(1x,a,i0)', advance = "no" ) "rmse_col", col_idx
    end do
    do col_idx = 1, size( rmse, 3 )
      write( *, '(1x,a,i0)' , advance = "no" ) "corr_col", col_idx
    end do
    do col_idx = 1, size( rmse, 3 )
      write( *, '(1x,a,i0)' , advance = "no" ) "std_ratio_col", col_idx
    end do
    do col_idx = 1, size( rmse, 3 )
      write( *, '(1x,a,i0)' , advance = "no" ) "bias_norm_col", col_idx
    end do
    write( *, '(1x,a)' ) "best_col"

    do window_idx = 1, size( rmse, 1 )
      do row_idx = 1, size( rmse, 2 )
        best_col_idx = 1
        best_col_loss = rmse(window_idx,row_idx,1)
        do col_idx = 2, size( rmse, 3 )
          if ( rmse(window_idx,row_idx,col_idx) < best_col_loss ) then
            best_col_idx = col_idx
            best_col_loss = rmse(window_idx,row_idx,col_idx)
          end if
        end do

        write( *, '(a)', advance = "no" ) trim( clubb_var_names(row_idx) )
        if ( size( rmse, 1 ) > 1 ) write( *, '(1x,a,i0)', advance = "no" ) "window", window_idx
        do col_idx = 1, size( rmse, 3 )
          write( *, '(1x,es16.8)', advance = "no" ) rmse(window_idx,row_idx,col_idx)
        end do
        do col_idx = 1, size( rmse, 3 )
          write( *, '(1x,es16.8)', advance = "no" ) correlation(window_idx,row_idx,col_idx)
        end do
        do col_idx = 1, size( rmse, 3 )
          write( *, '(1x,es16.8)', advance = "no" ) std_ratio(window_idx,row_idx,col_idx)
        end do
        do col_idx = 1, size( rmse, 3 )
          write( *, '(1x,es16.8)', advance = "no" ) bias_norm(window_idx,row_idx,col_idx)
        end do
        write( *, '(1x,a,i0)' ) "col", best_col_idx
      end do
    end do

    write( *, '(a)', advance = "no" ) "case_loss"
    do col_idx = 1, size( case_loss )
      write( *, '(1x,es16.8)', advance = "no" ) case_loss(col_idx)
    end do
    write( *, * )

  end subroutine print_loss_matrix

end program clubb_standalone_loss
