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
    scaled_rmse, &
    correlation, &
    std_ratio, &
    centered_rmse_norm, &
    bias_norm

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

  call clubb_get_loss( trim( namelist_filename ), clubb_var_names, scaled_rmse, correlation, std_ratio, &
                       centered_rmse_norm, bias_norm, err_info )

  call print_loss_matrix( clubb_var_names, scaled_rmse, correlation, std_ratio, centered_rmse_norm, bias_norm )

  call cleanup_err_info_api( err_info )
  call exit( success_code )

contains

  subroutine print_loss_matrix( clubb_var_names, scaled_rmse, correlation, std_ratio, centered_rmse_norm, bias_norm )

    ! Description:
  !   Print one row per variable/window with the explicit loss-driver metrics
  !   and best-performing parameter column for each variable.

    implicit none

    character(len=*), dimension(:), intent(in) :: &
      clubb_var_names

    real( kind = core_rknd ), dimension(:,:,:), intent(in) :: &
      scaled_rmse, &
      correlation, &
      std_ratio, &
      centered_rmse_norm, &
      bias_norm

    integer :: &
      row_idx, &
      col_idx, &
      window_idx, &
      best_col_idx

    real( kind = core_rknd ) :: &
      best_col_loss

    write( *, '(a)', advance = "no" ) "variable"
    if ( size( scaled_rmse, 1 ) > 1 ) write( *, '(1x,a)', advance = "no" ) "window"
    do col_idx = 1, size( scaled_rmse, 3 )
      write( *, '(1x,a,i0)', advance = "no" ) "scaled_rmse_col", col_idx
    end do
    do col_idx = 1, size( scaled_rmse, 3 )
      write( *, '(1x,a,i0)' , advance = "no" ) "corr_col", col_idx
    end do
    do col_idx = 1, size( scaled_rmse, 3 )
      write( *, '(1x,a,i0)' , advance = "no" ) "std_ratio_col", col_idx
    end do
    do col_idx = 1, size( scaled_rmse, 3 )
      write( *, '(1x,a,i0)' , advance = "no" ) "crmse_norm_col", col_idx
    end do
    do col_idx = 1, size( scaled_rmse, 3 )
      write( *, '(1x,a,i0)' , advance = "no" ) "bias_norm_col", col_idx
    end do
    write( *, '(1x,a)' ) "best_col"

    do window_idx = 1, size( scaled_rmse, 1 )
      do row_idx = 1, size( scaled_rmse, 2 )
        best_col_idx = 1
        best_col_loss = scaled_rmse(window_idx,row_idx,1)
        do col_idx = 2, size( scaled_rmse, 3 )
          if ( scaled_rmse(window_idx,row_idx,col_idx) < best_col_loss ) then
            best_col_idx = col_idx
            best_col_loss = scaled_rmse(window_idx,row_idx,col_idx)
          end if
        end do

        write( *, '(a)', advance = "no" ) trim( clubb_var_names(row_idx) )
        if ( size( scaled_rmse, 1 ) > 1 ) write( *, '(1x,a,i0)', advance = "no" ) "window", window_idx
        do col_idx = 1, size( scaled_rmse, 3 )
          write( *, '(1x,es16.8)', advance = "no" ) scaled_rmse(window_idx,row_idx,col_idx)
        end do
        do col_idx = 1, size( scaled_rmse, 3 )
          write( *, '(1x,es16.8)', advance = "no" ) correlation(window_idx,row_idx,col_idx)
        end do
        do col_idx = 1, size( scaled_rmse, 3 )
          write( *, '(1x,es16.8)', advance = "no" ) std_ratio(window_idx,row_idx,col_idx)
        end do
        do col_idx = 1, size( scaled_rmse, 3 )
          write( *, '(1x,es16.8)', advance = "no" ) centered_rmse_norm(window_idx,row_idx,col_idx)
        end do
        do col_idx = 1, size( scaled_rmse, 3 )
          write( *, '(1x,es16.8)', advance = "no" ) bias_norm(window_idx,row_idx,col_idx)
        end do
        write( *, '(1x,a,i0)' ) "col", best_col_idx
      end do
    end do

  end subroutine print_loss_matrix

end program clubb_standalone_loss
