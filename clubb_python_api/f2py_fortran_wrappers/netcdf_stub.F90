! netcdf_stub.F90 — Minimal stub for the netcdf Fortran module.
!
! Provides only the constants and interfaces used by stats_netcdf.F90.
! All functions return NF90_NOERR (success) but do nothing.
! This is safe because stats are initialized with enabled=.false.,
! so no NetCDF I/O routines are actually called at runtime.
!
module netcdf
  implicit none
  private

  integer, parameter, public :: NF90_NOERR     = 0
  integer, parameter, public :: NF90_CLOBBER   = 0
  integer, parameter, public :: NF90_UNLIMITED = 0
  integer, parameter, public :: NF90_DOUBLE    = 6
  integer, parameter, public :: NF90_CHAR      = 2
  integer, parameter, public :: NF90_FILL_DOUBLE = -1  ! placeholder

  public :: nf90_create, nf90_redef, nf90_def_dim, nf90_def_var
  public :: nf90_def_var_fill, nf90_put_att, nf90_enddef
  public :: nf90_put_var, nf90_close, nf90_strerror

  interface nf90_put_var
    module procedure nf90_put_var_r8_1d
    module procedure nf90_put_var_r8_2d
    module procedure nf90_put_var_r8_4d
    module procedure nf90_put_var_char_1d
    module procedure nf90_put_var_scalar_r8
  end interface

  interface nf90_put_att
    module procedure nf90_put_att_text
  end interface

contains

  integer function nf90_create(path, cmode, ncid)
    character(len=*), intent(in)  :: path
    integer,          intent(in)  :: cmode
    integer,          intent(out) :: ncid
    ncid = -1
    nf90_create = NF90_NOERR
  end function

  integer function nf90_redef(ncid)
    integer, intent(in) :: ncid
    nf90_redef = NF90_NOERR
  end function

  integer function nf90_def_dim(ncid, name, len, dimid)
    integer,          intent(in)  :: ncid
    character(len=*), intent(in)  :: name
    integer,          intent(in)  :: len
    integer,          intent(out) :: dimid
    dimid = -1
    nf90_def_dim = NF90_NOERR
  end function

  integer function nf90_def_var(ncid, name, xtype, dimids, varid)
    integer,          intent(in)  :: ncid
    character(len=*), intent(in)  :: name
    integer,          intent(in)  :: xtype
    integer,          intent(in)  :: dimids(:)
    integer,          intent(out) :: varid
    varid = -1
    nf90_def_var = NF90_NOERR
  end function

  integer function nf90_def_var_fill(ncid, varid, no_fill, fill_value)
    integer, intent(in) :: ncid, varid, no_fill
    real(8), intent(in) :: fill_value
    nf90_def_var_fill = NF90_NOERR
  end function

  integer function nf90_put_att_text(ncid, varid, name, values)
    integer,          intent(in) :: ncid, varid
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: values
    nf90_put_att_text = NF90_NOERR
  end function

  integer function nf90_enddef(ncid)
    integer, intent(in) :: ncid
    nf90_enddef = NF90_NOERR
  end function

  integer function nf90_put_var_r8_1d(ncid, varid, values, start, count)
    integer, intent(in) :: ncid, varid
    real(8), intent(in) :: values(:)
    integer, intent(in), optional :: start(:), count(:)
    nf90_put_var_r8_1d = NF90_NOERR
  end function

  integer function nf90_put_var_r8_2d(ncid, varid, values, start, count)
    integer, intent(in) :: ncid, varid
    real(8), intent(in) :: values(:,:)
    integer, intent(in), optional :: start(:), count(:)
    nf90_put_var_r8_2d = NF90_NOERR
  end function

  integer function nf90_put_var_r8_4d(ncid, varid, values, start, count)
    integer, intent(in) :: ncid, varid
    real(8), intent(in) :: values(:,:,:,:)
    integer, intent(in), optional :: start(:), count(:)
    nf90_put_var_r8_4d = NF90_NOERR
  end function

  integer function nf90_put_var_char_1d(ncid, varid, values)
    integer,          intent(in) :: ncid, varid
    character(len=*), intent(in) :: values(:)
    nf90_put_var_char_1d = NF90_NOERR
  end function

  integer function nf90_put_var_scalar_r8(ncid, varid, values, start)
    integer, intent(in) :: ncid, varid
    real(8), intent(in) :: values
    integer, intent(in), optional :: start(:)
    nf90_put_var_scalar_r8 = NF90_NOERR
  end function

  integer function nf90_close(ncid)
    integer, intent(in) :: ncid
    nf90_close = NF90_NOERR
  end function

  function nf90_strerror(ncerr) result(str)
    integer, intent(in)  :: ncerr
    character(len=80)    :: str
    str = "stub: no NetCDF error info"
  end function

end module netcdf
