! Temporary fieldwise conversions for the C++/Fortran call boundary.
! to_native_* makes a local Fortran value from C field descriptors;
! from_native_* returns every field in C-owned buffers. The BIND(C)
! roundtrip procedures exercise both directions without saved state.
module f2k_derived_type_bridge
  use, intrinsic :: iso_c_binding, only: &
    c_ptr, c_null_ptr, c_int, c_int64_t, c_int8_t, c_size_t, c_char, &
    c_f_pointer, c_associated, c_sizeof

  use clubb_precision, only: core_rknd

  use grid_class, only: &
    grid

  use array_index, only: &
    sclr_idx_type

  use model_flags, only: &
    clubb_config_flags_type

  use parameters_tunable, only: &
    nu_vertical_res_dep

  use pdf_parameter_module, only: &
    pdf_parameter, &
    implicit_coefs_terms

  use err_info_type_module, only: &
    err_info_type, &
    f2k_get_err_info_private_api, &
    f2k_set_err_info_private_api

  use corr_varnce_module, only: &
    hm_metadata_type

  use parameters_silhs, only: &
    silhs_config_flags_type

  use hydromet_pdf_parameter_module, only: &
    hydromet_pdf_parameter, &
    precipitation_fractions
  implicit none
  private
  public :: wire_field, allocate_field, empty_field
  public :: to_native_grid, from_native_grid, &
    to_native_sclr_idx_type, from_native_sclr_idx_type, &
    to_native_clubb_config_flags_type, from_native_clubb_config_flags_type, &
    to_native_nu_vertical_res_dep, from_native_nu_vertical_res_dep, &
    to_native_pdf_parameter, from_native_pdf_parameter, &
    to_native_implicit_coefs_terms, from_native_implicit_coefs_terms, &
    to_native_err_info_type, from_native_err_info_type, &
    to_native_hm_metadata_type, from_native_hm_metadata_type, &
    to_native_silhs_config_flags_type, from_native_silhs_config_flags_type, &
    to_native_hydromet_pdf_parameter, from_native_hydromet_pdf_parameter, &
    to_native_precipitation_fractions, from_native_precipitation_fractions
  public :: to_native_hydromet_pdf_parameter_collection, from_native_hydromet_pdf_parameter_collection
  type, bind(C) :: wire_field
    type(c_ptr) :: data
    integer(c_int64_t) :: extent(4)
    integer(c_int) :: rank, allocated, char_len
  end type wire_field
  interface
    function c_malloc(nbytes) bind(C, name='malloc') result(ptr)
      import c_size_t, c_ptr
      integer(c_size_t), value :: nbytes
      type(c_ptr) :: ptr
    end function c_malloc
  end interface
contains

  subroutine allocate_field(rank, extent, nbytes, char_len, field)
    ! Description: Allocate C-owned output bytes for one converted field.

    implicit none

    !--------------------------- Input Variables ---------------------------
    integer, intent(in) :: rank ! Number of array dimensions [-]
    integer(c_int64_t), intent(in) :: extent(4) ! Array extents in native dimension order [-]
    integer(c_size_t), intent(in) :: nbytes ! Output buffer size [bytes]
    integer, intent(in) :: char_len ! Fixed character width; zero for numeric fields [characters]

    !--------------------------- Output Variables ---------------------------
    type(wire_field), intent(out) :: field ! Descriptor of one packed field [-]

    !--------------------------- Begin Code ---------------------------

    field%data = c_malloc(max(1_c_size_t, nbytes))

    if (.not. c_associated(field%data)) error stop 'f2k malloc failed'

    field%rank = rank
    field%extent = extent
    field%allocated = 1
    field%char_len = char_len

  end subroutine allocate_field

  subroutine empty_field(rank, char_len, field)
    ! Description: Mark an unallocated Fortran component in the C descriptor.

    implicit none

    !--------------------------- Input Variables ---------------------------
    integer, intent(in) :: rank ! Number of array dimensions [-]
    integer, intent(in) :: char_len ! Fixed character width; zero for numeric fields [characters]

    !--------------------------- Output Variables ---------------------------
    type(wire_field), intent(out) :: field ! Descriptor of one packed field [-]

    !--------------------------- Begin Code ---------------------------

    field%data = c_null_ptr
    field%extent = 0_c_int64_t
    field%rank = rank
    field%allocated = 0
    field%char_len = char_len

  end subroutine empty_field

  subroutine to_native_grid(input, value)
    ! Description: Copy each C field into a local grid value.

    implicit none

    !--------------------------- Input Variables ---------------------------
    type(wire_field), intent(in) :: input(*) ! Caller-owned input field descriptors [-]

    !--------------------------- Output Variables ---------------------------
    type(grid), intent(out) :: value ! Call-local native argument; units and shape follow the source contract

    !--------------------------- Local Variables ---------------------------
    integer(c_int), pointer :: ip(:) ! Borrowed view of packed integer data [-]
    ! Borrowed view of packed real data; units follow the native field
    real( kind = core_rknd ), pointer :: rp(:)
    integer :: n ! Number of packed elements [-]

    !--------------------------- Begin Code ---------------------------

    ! nzm

    if (input(1)%rank /= 0) error stop 'f2k rank: nzm'

    if (.not. c_associated(input(1)%data)) error stop 'f2k scalar: nzm'

    call c_f_pointer(input(1)%data, ip, [1])
    value%nzm = ip(1)

    ! nzt

    if (input(2)%rank /= 0) error stop 'f2k rank: nzt'

    if (.not. c_associated(input(2)%data)) error stop 'f2k scalar: nzt'

    call c_f_pointer(input(2)%data, ip, [1])
    value%nzt = ip(1)

    ! zm

    if (input(3)%rank /= 2) error stop 'f2k rank: zm'

    if (input(3)%allocated /= 0) then

      allocate(value%zm(int(input(3)%extent(1)), int(input(3)%extent(2))))
      n = int(product(input(3)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(3)%data)) error stop 'f2k data: zm'

        call c_f_pointer(input(3)%data, rp, [n])
        value%zm = reshape(rp, shape(value%zm))

      end if

    end if

    ! zt

    if (input(4)%rank /= 2) error stop 'f2k rank: zt'

    if (input(4)%allocated /= 0) then

      allocate(value%zt(int(input(4)%extent(1)), int(input(4)%extent(2))))
      n = int(product(input(4)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(4)%data)) error stop 'f2k data: zt'

        call c_f_pointer(input(4)%data, rp, [n])
        value%zt = reshape(rp, shape(value%zt))

      end if

    end if

    ! invrs_dzm

    if (input(5)%rank /= 2) error stop 'f2k rank: invrs_dzm'

    if (input(5)%allocated /= 0) then

      allocate(value%invrs_dzm(int(input(5)%extent(1)), int(input(5)%extent(2))))
      n = int(product(input(5)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(5)%data)) error stop 'f2k data: invrs_dzm'

        call c_f_pointer(input(5)%data, rp, [n])
        value%invrs_dzm = reshape(rp, shape(value%invrs_dzm))

      end if

    end if

    ! invrs_dzt

    if (input(6)%rank /= 2) error stop 'f2k rank: invrs_dzt'

    if (input(6)%allocated /= 0) then

      allocate(value%invrs_dzt(int(input(6)%extent(1)), int(input(6)%extent(2))))
      n = int(product(input(6)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(6)%data)) error stop 'f2k data: invrs_dzt'

        call c_f_pointer(input(6)%data, rp, [n])
        value%invrs_dzt = reshape(rp, shape(value%invrs_dzt))

      end if

    end if

    ! dzm

    if (input(7)%rank /= 2) error stop 'f2k rank: dzm'

    if (input(7)%allocated /= 0) then

      allocate(value%dzm(int(input(7)%extent(1)), int(input(7)%extent(2))))
      n = int(product(input(7)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(7)%data)) error stop 'f2k data: dzm'

        call c_f_pointer(input(7)%data, rp, [n])
        value%dzm = reshape(rp, shape(value%dzm))

      end if

    end if

    ! dzt

    if (input(8)%rank /= 2) error stop 'f2k rank: dzt'

    if (input(8)%allocated /= 0) then

      allocate(value%dzt(int(input(8)%extent(1)), int(input(8)%extent(2))))
      n = int(product(input(8)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(8)%data)) error stop 'f2k data: dzt'

        call c_f_pointer(input(8)%data, rp, [n])
        value%dzt = reshape(rp, shape(value%dzt))

      end if

    end if

    ! weights_zm2zt

    if (input(9)%rank /= 3) error stop 'f2k rank: weights_zm2zt'

    if (input(9)%allocated /= 0) then

      allocate(value%weights_zm2zt(int(input(9)%extent(1)), int(input(9)%extent(2)), int(input(9)%extent(3))))
      n = int(product(input(9)%extent(1:3)))

      if (n > 0) then

        if (.not. c_associated(input(9)%data)) error stop 'f2k data: weights_zm2zt'

        call c_f_pointer(input(9)%data, rp, [n])
        value%weights_zm2zt = reshape(rp, shape(value%weights_zm2zt))

      end if

    end if

    ! weights_zt2zm

    if (input(10)%rank /= 3) error stop 'f2k rank: weights_zt2zm'

    if (input(10)%allocated /= 0) then

      allocate(value%weights_zt2zm(int(input(10)%extent(1)), int(input(10)%extent(2)), int(input(10)%extent(3))))
      n = int(product(input(10)%extent(1:3)))

      if (n > 0) then

        if (.not. c_associated(input(10)%data)) error stop 'f2k data: weights_zt2zm'

        call c_f_pointer(input(10)%data, rp, [n])
        value%weights_zt2zm = reshape(rp, shape(value%weights_zt2zm))

      end if

    end if

    ! k_lb_zm

    if (input(11)%rank /= 0) error stop 'f2k rank: k_lb_zm'

    if (.not. c_associated(input(11)%data)) error stop 'f2k scalar: k_lb_zm'

    call c_f_pointer(input(11)%data, ip, [1])
    value%k_lb_zm = ip(1)

    ! k_lb_zt

    if (input(12)%rank /= 0) error stop 'f2k rank: k_lb_zt'

    if (.not. c_associated(input(12)%data)) error stop 'f2k scalar: k_lb_zt'

    call c_f_pointer(input(12)%data, ip, [1])
    value%k_lb_zt = ip(1)

    ! k_ub_zm

    if (input(13)%rank /= 0) error stop 'f2k rank: k_ub_zm'

    if (.not. c_associated(input(13)%data)) error stop 'f2k scalar: k_ub_zm'

    call c_f_pointer(input(13)%data, ip, [1])
    value%k_ub_zm = ip(1)

    ! k_ub_zt

    if (input(14)%rank /= 0) error stop 'f2k rank: k_ub_zt'

    if (.not. c_associated(input(14)%data)) error stop 'f2k scalar: k_ub_zt'

    call c_f_pointer(input(14)%data, ip, [1])
    value%k_ub_zt = ip(1)

    ! grid_dir_indx

    if (input(15)%rank /= 0) error stop 'f2k rank: grid_dir_indx'

    if (.not. c_associated(input(15)%data)) error stop 'f2k scalar: grid_dir_indx'

    call c_f_pointer(input(15)%data, ip, [1])
    value%grid_dir_indx = ip(1)

    ! grid_dir

    if (input(16)%rank /= 0) error stop 'f2k rank: grid_dir'

    if (.not. c_associated(input(16)%data)) error stop 'f2k scalar: grid_dir'

    call c_f_pointer(input(16)%data, rp, [1])
    value%grid_dir = rp(1)

  end subroutine to_native_grid

  subroutine from_native_grid(value, output)
    ! Description: Copy each grid field into C-owned output bytes.

    implicit none

    !--------------------------- Input Variables ---------------------------
    type(grid), intent(in) :: value ! Call-local native argument; units and shape follow the source contract

    !--------------------------- Output Variables ---------------------------
    type(wire_field), intent(out) :: output(*) ! Returned descriptors; C++ releases their C-owned buffers [-]

    !--------------------------- Local Variables ---------------------------
    integer(c_int), pointer :: ip(:) ! Borrowed view of packed integer data [-]
    ! Borrowed view of packed real data; units follow the native field
    real( kind = core_rknd ), pointer :: rp(:)
    integer :: n ! Number of packed elements [-]
    integer(c_int64_t) :: field_shape(4) ! Expected field extents in native dimension order [-]

    !--------------------------- Begin Code ---------------------------

    ! nzm

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(1)) ! Out

    call c_f_pointer(output(1)%data, ip, [1])
    ip(1) = value%nzm

    ! nzt

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(2)) ! Out

    call c_f_pointer(output(2)%data, ip, [1])
    ip(1) = value%nzt

    ! zm

    if (allocated(value%zm)) then

      n = size(value%zm)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%zm, 1, kind=c_int64_t)
      field_shape(2) = size(value%zm, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(3)) ! Out

      call c_f_pointer(output(3)%data, rp, [n])
      rp = reshape(value%zm, [n])

    else

      call empty_field(2, 0, & ! In
                       output(3)) ! Out

    end if

    ! zt

    if (allocated(value%zt)) then

      n = size(value%zt)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%zt, 1, kind=c_int64_t)
      field_shape(2) = size(value%zt, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(4)) ! Out

      call c_f_pointer(output(4)%data, rp, [n])
      rp = reshape(value%zt, [n])

    else

      call empty_field(2, 0, & ! In
                       output(4)) ! Out

    end if

    ! invrs_dzm

    if (allocated(value%invrs_dzm)) then

      n = size(value%invrs_dzm)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%invrs_dzm, 1, kind=c_int64_t)
      field_shape(2) = size(value%invrs_dzm, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(5)) ! Out

      call c_f_pointer(output(5)%data, rp, [n])
      rp = reshape(value%invrs_dzm, [n])

    else

      call empty_field(2, 0, & ! In
                       output(5)) ! Out

    end if

    ! invrs_dzt

    if (allocated(value%invrs_dzt)) then

      n = size(value%invrs_dzt)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%invrs_dzt, 1, kind=c_int64_t)
      field_shape(2) = size(value%invrs_dzt, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(6)) ! Out

      call c_f_pointer(output(6)%data, rp, [n])
      rp = reshape(value%invrs_dzt, [n])

    else

      call empty_field(2, 0, & ! In
                       output(6)) ! Out

    end if

    ! dzm

    if (allocated(value%dzm)) then

      n = size(value%dzm)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%dzm, 1, kind=c_int64_t)
      field_shape(2) = size(value%dzm, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(7)) ! Out

      call c_f_pointer(output(7)%data, rp, [n])
      rp = reshape(value%dzm, [n])

    else

      call empty_field(2, 0, & ! In
                       output(7)) ! Out

    end if

    ! dzt

    if (allocated(value%dzt)) then

      n = size(value%dzt)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%dzt, 1, kind=c_int64_t)
      field_shape(2) = size(value%dzt, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(8)) ! Out

      call c_f_pointer(output(8)%data, rp, [n])
      rp = reshape(value%dzt, [n])

    else

      call empty_field(2, 0, & ! In
                       output(8)) ! Out

    end if

    ! weights_zm2zt

    if (allocated(value%weights_zm2zt)) then

      n = size(value%weights_zm2zt)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%weights_zm2zt, 1, kind=c_int64_t)
      field_shape(2) = size(value%weights_zm2zt, 2, kind=c_int64_t)
      field_shape(3) = size(value%weights_zm2zt, 3, kind=c_int64_t)

      call allocate_field(3, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(9)) ! Out

      call c_f_pointer(output(9)%data, rp, [n])
      rp = reshape(value%weights_zm2zt, [n])

    else

      call empty_field(3, 0, & ! In
                       output(9)) ! Out

    end if

    ! weights_zt2zm

    if (allocated(value%weights_zt2zm)) then

      n = size(value%weights_zt2zm)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%weights_zt2zm, 1, kind=c_int64_t)
      field_shape(2) = size(value%weights_zt2zm, 2, kind=c_int64_t)
      field_shape(3) = size(value%weights_zt2zm, 3, kind=c_int64_t)

      call allocate_field(3, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(10)) ! Out

      call c_f_pointer(output(10)%data, rp, [n])
      rp = reshape(value%weights_zt2zm, [n])

    else

      call empty_field(3, 0, & ! In
                       output(10)) ! Out

    end if

    ! k_lb_zm

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(11)) ! Out

    call c_f_pointer(output(11)%data, ip, [1])
    ip(1) = value%k_lb_zm

    ! k_lb_zt

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(12)) ! Out

    call c_f_pointer(output(12)%data, ip, [1])
    ip(1) = value%k_lb_zt

    ! k_ub_zm

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(13)) ! Out

    call c_f_pointer(output(13)%data, ip, [1])
    ip(1) = value%k_ub_zm

    ! k_ub_zt

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(14)) ! Out

    call c_f_pointer(output(14)%data, ip, [1])
    ip(1) = value%k_ub_zt

    ! grid_dir_indx

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(15)) ! Out

    call c_f_pointer(output(15)%data, ip, [1])
    ip(1) = value%grid_dir_indx

    ! grid_dir

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0.0_core_rknd), 0, & ! In
                        output(16)) ! Out

    call c_f_pointer(output(16)%data, rp, [1])
    rp(1) = value%grid_dir

  end subroutine from_native_grid

  subroutine f2k_roundtrip_grid_c(input, output) bind(C)
    ! Description: Exercise both converters with one local Fortran value.

    implicit none

    !--------------------------- Input Variables ---------------------------
    type(wire_field), intent(in) :: input(*) ! Caller-owned input field descriptors [-]

    !--------------------------- Output Variables ---------------------------
    type(wire_field), intent(out) :: output(*) ! Returned descriptors; C++ releases their C-owned buffers [-]

    !--------------------------- Local Variables ---------------------------
    type(grid) :: value ! Call-local native argument; units and shape follow the source contract

    !--------------------------- Begin Code ---------------------------

    call to_native_grid(input, value)
    call from_native_grid(value, output)

  end subroutine f2k_roundtrip_grid_c

  subroutine to_native_sclr_idx_type(input, value)
    ! Description: Copy each C field into a local sclr_idx_type value.

    implicit none

    !--------------------------- Input Variables ---------------------------
    type(wire_field), intent(in) :: input(*) ! Caller-owned input field descriptors [-]

    !--------------------------- Output Variables ---------------------------
    ! Call-local native argument; units and shape follow the source contract
    type(sclr_idx_type), intent(out) :: value

    !--------------------------- Local Variables ---------------------------
    integer(c_int), pointer :: ip(:) ! Borrowed view of packed integer data [-]

    !--------------------------- Begin Code ---------------------------

    ! iisclr_rt

    if (input(1)%rank /= 0) error stop 'f2k rank: iisclr_rt'

    if (.not. c_associated(input(1)%data)) error stop 'f2k scalar: iisclr_rt'

    call c_f_pointer(input(1)%data, ip, [1])
    value%iisclr_rt = ip(1)

    ! iisclr_thl

    if (input(2)%rank /= 0) error stop 'f2k rank: iisclr_thl'

    if (.not. c_associated(input(2)%data)) error stop 'f2k scalar: iisclr_thl'

    call c_f_pointer(input(2)%data, ip, [1])
    value%iisclr_thl = ip(1)

    ! iisclr_CO2

    if (input(3)%rank /= 0) error stop 'f2k rank: iisclr_CO2'

    if (.not. c_associated(input(3)%data)) error stop 'f2k scalar: iisclr_CO2'

    call c_f_pointer(input(3)%data, ip, [1])
    value%iisclr_CO2 = ip(1)

    ! iiedsclr_rt

    if (input(4)%rank /= 0) error stop 'f2k rank: iiedsclr_rt'

    if (.not. c_associated(input(4)%data)) error stop 'f2k scalar: iiedsclr_rt'

    call c_f_pointer(input(4)%data, ip, [1])
    value%iiedsclr_rt = ip(1)

    ! iiedsclr_thl

    if (input(5)%rank /= 0) error stop 'f2k rank: iiedsclr_thl'

    if (.not. c_associated(input(5)%data)) error stop 'f2k scalar: iiedsclr_thl'

    call c_f_pointer(input(5)%data, ip, [1])
    value%iiedsclr_thl = ip(1)

    ! iiedsclr_CO2

    if (input(6)%rank /= 0) error stop 'f2k rank: iiedsclr_CO2'

    if (.not. c_associated(input(6)%data)) error stop 'f2k scalar: iiedsclr_CO2'

    call c_f_pointer(input(6)%data, ip, [1])
    value%iiedsclr_CO2 = ip(1)

  end subroutine to_native_sclr_idx_type

  subroutine from_native_sclr_idx_type(value, output)
    ! Description: Copy each sclr_idx_type field into C-owned output bytes.

    implicit none

    !--------------------------- Input Variables ---------------------------
    ! Call-local native argument; units and shape follow the source contract
    type(sclr_idx_type), intent(in) :: value

    !--------------------------- Output Variables ---------------------------
    type(wire_field), intent(out) :: output(*) ! Returned descriptors; C++ releases their C-owned buffers [-]

    !--------------------------- Local Variables ---------------------------
    integer(c_int), pointer :: ip(:) ! Borrowed view of packed integer data [-]

    !--------------------------- Begin Code ---------------------------

    ! iisclr_rt

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(1)) ! Out

    call c_f_pointer(output(1)%data, ip, [1])
    ip(1) = value%iisclr_rt

    ! iisclr_thl

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(2)) ! Out

    call c_f_pointer(output(2)%data, ip, [1])
    ip(1) = value%iisclr_thl

    ! iisclr_CO2

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(3)) ! Out

    call c_f_pointer(output(3)%data, ip, [1])
    ip(1) = value%iisclr_CO2

    ! iiedsclr_rt

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(4)) ! Out

    call c_f_pointer(output(4)%data, ip, [1])
    ip(1) = value%iiedsclr_rt

    ! iiedsclr_thl

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(5)) ! Out

    call c_f_pointer(output(5)%data, ip, [1])
    ip(1) = value%iiedsclr_thl

    ! iiedsclr_CO2

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(6)) ! Out

    call c_f_pointer(output(6)%data, ip, [1])
    ip(1) = value%iiedsclr_CO2

  end subroutine from_native_sclr_idx_type

  subroutine f2k_roundtrip_sclr_idx_type_c(input, output) bind(C)
    ! Description: Exercise both converters with one local Fortran value.

    implicit none

    !--------------------------- Input Variables ---------------------------
    type(wire_field), intent(in) :: input(*) ! Caller-owned input field descriptors [-]

    !--------------------------- Output Variables ---------------------------
    type(wire_field), intent(out) :: output(*) ! Returned descriptors; C++ releases their C-owned buffers [-]

    !--------------------------- Local Variables ---------------------------
    type(sclr_idx_type) :: value ! Call-local native argument; units and shape follow the source contract

    !--------------------------- Begin Code ---------------------------

    call to_native_sclr_idx_type(input, value)
    call from_native_sclr_idx_type(value, output)

  end subroutine f2k_roundtrip_sclr_idx_type_c

  subroutine to_native_clubb_config_flags_type(input, value)
    ! Description: Copy each C field into a local clubb_config_flags_type value.

    implicit none

    !--------------------------- Input Variables ---------------------------
    type(wire_field), intent(in) :: input(*) ! Caller-owned input field descriptors [-]

    !--------------------------- Output Variables ---------------------------
    ! Call-local native argument; units and shape follow the source contract
    type(clubb_config_flags_type), intent(out) :: value

    !--------------------------- Local Variables ---------------------------
    integer(c_int), pointer :: ip(:) ! Borrowed view of packed integer data [-]
    integer(c_int8_t), pointer :: bp(:) ! Borrowed view of packed logical bytes [-]

    !--------------------------- Begin Code ---------------------------

    ! iiPDF_type

    if (input(1)%rank /= 0) error stop 'f2k rank: iiPDF_type'

    if (.not. c_associated(input(1)%data)) error stop 'f2k scalar: iiPDF_type'

    call c_f_pointer(input(1)%data, ip, [1])
    value%iiPDF_type = ip(1)

    ! ipdf_call_placement

    if (input(2)%rank /= 0) error stop 'f2k rank: ipdf_call_placement'

    if (.not. c_associated(input(2)%data)) error stop 'f2k scalar: ipdf_call_placement'

    call c_f_pointer(input(2)%data, ip, [1])
    value%ipdf_call_placement = ip(1)

    ! penta_solve_method

    if (input(3)%rank /= 0) error stop 'f2k rank: penta_solve_method'

    if (.not. c_associated(input(3)%data)) error stop 'f2k scalar: penta_solve_method'

    call c_f_pointer(input(3)%data, ip, [1])
    value%penta_solve_method = ip(1)

    ! tridiag_solve_method

    if (input(4)%rank /= 0) error stop 'f2k rank: tridiag_solve_method'

    if (.not. c_associated(input(4)%data)) error stop 'f2k scalar: tridiag_solve_method'

    call c_f_pointer(input(4)%data, ip, [1])
    value%tridiag_solve_method = ip(1)

    ! saturation_formula

    if (input(5)%rank /= 0) error stop 'f2k rank: saturation_formula'

    if (.not. c_associated(input(5)%data)) error stop 'f2k scalar: saturation_formula'

    call c_f_pointer(input(5)%data, ip, [1])
    value%saturation_formula = ip(1)

    ! Lscale_method

    if (input(6)%rank /= 0) error stop 'f2k rank: Lscale_method'

    if (.not. c_associated(input(6)%data)) error stop 'f2k scalar: Lscale_method'

    call c_f_pointer(input(6)%data, ip, [1])
    value%Lscale_method = ip(1)

    ! grid_remap_method

    if (input(7)%rank /= 0) error stop 'f2k rank: grid_remap_method'

    if (.not. c_associated(input(7)%data)) error stop 'f2k scalar: grid_remap_method'

    call c_f_pointer(input(7)%data, ip, [1])
    value%grid_remap_method = ip(1)

    ! grid_adapt_in_time_method

    if (input(8)%rank /= 0) error stop 'f2k rank: grid_adapt_in_time_method'

    if (.not. c_associated(input(8)%data)) error stop 'f2k scalar: grid_adapt_in_time_method'

    call c_f_pointer(input(8)%data, ip, [1])
    value%grid_adapt_in_time_method = ip(1)

    ! fill_holes_type

    if (input(9)%rank /= 0) error stop 'f2k rank: fill_holes_type'

    if (.not. c_associated(input(9)%data)) error stop 'f2k scalar: fill_holes_type'

    call c_f_pointer(input(9)%data, ip, [1])
    value%fill_holes_type = ip(1)

    ! l_use_precip_frac

    if (input(10)%rank /= 0) error stop 'f2k rank: l_use_precip_frac'

    if (.not. c_associated(input(10)%data)) error stop 'f2k scalar: l_use_precip_frac'

    call c_f_pointer(input(10)%data, bp, [1])
    value%l_use_precip_frac = bp(1) /= 0_c_int8_t

    ! l_predict_upwp_vpwp

    if (input(11)%rank /= 0) error stop 'f2k rank: l_predict_upwp_vpwp'

    if (.not. c_associated(input(11)%data)) error stop 'f2k scalar: l_predict_upwp_vpwp'

    call c_f_pointer(input(11)%data, bp, [1])
    value%l_predict_upwp_vpwp = bp(1) /= 0_c_int8_t

    ! l_ho_nontrad_coriolis

    if (input(12)%rank /= 0) error stop 'f2k rank: l_ho_nontrad_coriolis'

    if (.not. c_associated(input(12)%data)) error stop 'f2k scalar: l_ho_nontrad_coriolis'

    call c_f_pointer(input(12)%data, bp, [1])
    value%l_ho_nontrad_coriolis = bp(1) /= 0_c_int8_t

    ! l_ho_trad_coriolis

    if (input(13)%rank /= 0) error stop 'f2k rank: l_ho_trad_coriolis'

    if (.not. c_associated(input(13)%data)) error stop 'f2k scalar: l_ho_trad_coriolis'

    call c_f_pointer(input(13)%data, bp, [1])
    value%l_ho_trad_coriolis = bp(1) /= 0_c_int8_t

    ! l_min_wp2_from_corr_wx

    if (input(14)%rank /= 0) error stop 'f2k rank: l_min_wp2_from_corr_wx'

    if (.not. c_associated(input(14)%data)) error stop 'f2k scalar: l_min_wp2_from_corr_wx'

    call c_f_pointer(input(14)%data, bp, [1])
    value%l_min_wp2_from_corr_wx = bp(1) /= 0_c_int8_t

    ! l_min_xp2_from_corr_wx

    if (input(15)%rank /= 0) error stop 'f2k rank: l_min_xp2_from_corr_wx'

    if (.not. c_associated(input(15)%data)) error stop 'f2k scalar: l_min_xp2_from_corr_wx'

    call c_f_pointer(input(15)%data, bp, [1])
    value%l_min_xp2_from_corr_wx = bp(1) /= 0_c_int8_t

    ! l_C2_cloud_frac

    if (input(16)%rank /= 0) error stop 'f2k rank: l_C2_cloud_frac'

    if (.not. c_associated(input(16)%data)) error stop 'f2k scalar: l_C2_cloud_frac'

    call c_f_pointer(input(16)%data, bp, [1])
    value%l_C2_cloud_frac = bp(1) /= 0_c_int8_t

    ! l_diffuse_rtm_and_thlm

    if (input(17)%rank /= 0) error stop 'f2k rank: l_diffuse_rtm_and_thlm'

    if (.not. c_associated(input(17)%data)) error stop 'f2k scalar: l_diffuse_rtm_and_thlm'

    call c_f_pointer(input(17)%data, bp, [1])
    value%l_diffuse_rtm_and_thlm = bp(1) /= 0_c_int8_t

    ! l_stability_correct_Kh_N2_zm

    if (input(18)%rank /= 0) error stop 'f2k rank: l_stability_correct_Kh_N2_zm'

    if (.not. c_associated(input(18)%data)) error stop 'f2k scalar: l_stability_correct_Kh_N2_zm'

    call c_f_pointer(input(18)%data, bp, [1])
    value%l_stability_correct_Kh_N2_zm = bp(1) /= 0_c_int8_t

    ! l_calc_thlp2_rad

    if (input(19)%rank /= 0) error stop 'f2k rank: l_calc_thlp2_rad'

    if (.not. c_associated(input(19)%data)) error stop 'f2k scalar: l_calc_thlp2_rad'

    call c_f_pointer(input(19)%data, bp, [1])
    value%l_calc_thlp2_rad = bp(1) /= 0_c_int8_t

    ! l_upwind_xpyp_ta

    if (input(20)%rank /= 0) error stop 'f2k rank: l_upwind_xpyp_ta'

    if (.not. c_associated(input(20)%data)) error stop 'f2k scalar: l_upwind_xpyp_ta'

    call c_f_pointer(input(20)%data, bp, [1])
    value%l_upwind_xpyp_ta = bp(1) /= 0_c_int8_t

    ! l_upwind_xm_ma

    if (input(21)%rank /= 0) error stop 'f2k rank: l_upwind_xm_ma'

    if (.not. c_associated(input(21)%data)) error stop 'f2k scalar: l_upwind_xm_ma'

    call c_f_pointer(input(21)%data, bp, [1])
    value%l_upwind_xm_ma = bp(1) /= 0_c_int8_t

    ! l_uv_nudge

    if (input(22)%rank /= 0) error stop 'f2k rank: l_uv_nudge'

    if (.not. c_associated(input(22)%data)) error stop 'f2k scalar: l_uv_nudge'

    call c_f_pointer(input(22)%data, bp, [1])
    value%l_uv_nudge = bp(1) /= 0_c_int8_t

    ! l_rtm_nudge

    if (input(23)%rank /= 0) error stop 'f2k rank: l_rtm_nudge'

    if (.not. c_associated(input(23)%data)) error stop 'f2k scalar: l_rtm_nudge'

    call c_f_pointer(input(23)%data, bp, [1])
    value%l_rtm_nudge = bp(1) /= 0_c_int8_t

    ! l_vert_avg_closure

    if (input(24)%rank /= 0) error stop 'f2k rank: l_vert_avg_closure'

    if (.not. c_associated(input(24)%data)) error stop 'f2k scalar: l_vert_avg_closure'

    call c_f_pointer(input(24)%data, bp, [1])
    value%l_vert_avg_closure = bp(1) /= 0_c_int8_t

    ! l_trapezoidal_rule_zt

    if (input(25)%rank /= 0) error stop 'f2k rank: l_trapezoidal_rule_zt'

    if (.not. c_associated(input(25)%data)) error stop 'f2k scalar: l_trapezoidal_rule_zt'

    call c_f_pointer(input(25)%data, bp, [1])
    value%l_trapezoidal_rule_zt = bp(1) /= 0_c_int8_t

    ! l_trapezoidal_rule_zm

    if (input(26)%rank /= 0) error stop 'f2k rank: l_trapezoidal_rule_zm'

    if (.not. c_associated(input(26)%data)) error stop 'f2k scalar: l_trapezoidal_rule_zm'

    call c_f_pointer(input(26)%data, bp, [1])
    value%l_trapezoidal_rule_zm = bp(1) /= 0_c_int8_t

    ! l_call_pdf_closure_twice

    if (input(27)%rank /= 0) error stop 'f2k rank: l_call_pdf_closure_twice'

    if (.not. c_associated(input(27)%data)) error stop 'f2k scalar: l_call_pdf_closure_twice'

    call c_f_pointer(input(27)%data, bp, [1])
    value%l_call_pdf_closure_twice = bp(1) /= 0_c_int8_t

    ! l_standard_term_ta

    if (input(28)%rank /= 0) error stop 'f2k rank: l_standard_term_ta'

    if (.not. c_associated(input(28)%data)) error stop 'f2k scalar: l_standard_term_ta'

    call c_f_pointer(input(28)%data, bp, [1])
    value%l_standard_term_ta = bp(1) /= 0_c_int8_t

    ! l_partial_upwind_wp3

    if (input(29)%rank /= 0) error stop 'f2k rank: l_partial_upwind_wp3'

    if (.not. c_associated(input(29)%data)) error stop 'f2k scalar: l_partial_upwind_wp3'

    call c_f_pointer(input(29)%data, bp, [1])
    value%l_partial_upwind_wp3 = bp(1) /= 0_c_int8_t

    ! l_godunov_upwind_wpxp_ta

    if (input(30)%rank /= 0) error stop 'f2k rank: l_godunov_upwind_wpxp_ta'

    if (.not. c_associated(input(30)%data)) error stop 'f2k scalar: l_godunov_upwind_wpxp_ta'

    call c_f_pointer(input(30)%data, bp, [1])
    value%l_godunov_upwind_wpxp_ta = bp(1) /= 0_c_int8_t

    ! l_godunov_upwind_xpyp_ta

    if (input(31)%rank /= 0) error stop 'f2k rank: l_godunov_upwind_xpyp_ta'

    if (.not. c_associated(input(31)%data)) error stop 'f2k scalar: l_godunov_upwind_xpyp_ta'

    call c_f_pointer(input(31)%data, bp, [1])
    value%l_godunov_upwind_xpyp_ta = bp(1) /= 0_c_int8_t

    ! l_use_cloud_cover

    if (input(32)%rank /= 0) error stop 'f2k rank: l_use_cloud_cover'

    if (.not. c_associated(input(32)%data)) error stop 'f2k scalar: l_use_cloud_cover'

    call c_f_pointer(input(32)%data, bp, [1])
    value%l_use_cloud_cover = bp(1) /= 0_c_int8_t

    ! l_diagnose_correlations

    if (input(33)%rank /= 0) error stop 'f2k rank: l_diagnose_correlations'

    if (.not. c_associated(input(33)%data)) error stop 'f2k scalar: l_diagnose_correlations'

    call c_f_pointer(input(33)%data, bp, [1])
    value%l_diagnose_correlations = bp(1) /= 0_c_int8_t

    ! l_calc_w_corr

    if (input(34)%rank /= 0) error stop 'f2k rank: l_calc_w_corr'

    if (.not. c_associated(input(34)%data)) error stop 'f2k scalar: l_calc_w_corr'

    call c_f_pointer(input(34)%data, bp, [1])
    value%l_calc_w_corr = bp(1) /= 0_c_int8_t

    ! l_const_Nc_in_cloud

    if (input(35)%rank /= 0) error stop 'f2k rank: l_const_Nc_in_cloud'

    if (.not. c_associated(input(35)%data)) error stop 'f2k scalar: l_const_Nc_in_cloud'

    call c_f_pointer(input(35)%data, bp, [1])
    value%l_const_Nc_in_cloud = bp(1) /= 0_c_int8_t

    ! l_fix_w_chi_eta_correlations

    if (input(36)%rank /= 0) error stop 'f2k rank: l_fix_w_chi_eta_correlations'

    if (.not. c_associated(input(36)%data)) error stop 'f2k scalar: l_fix_w_chi_eta_correlations'

    call c_f_pointer(input(36)%data, bp, [1])
    value%l_fix_w_chi_eta_correlations = bp(1) /= 0_c_int8_t

    ! l_stability_correct_tau_zm

    if (input(37)%rank /= 0) error stop 'f2k rank: l_stability_correct_tau_zm'

    if (.not. c_associated(input(37)%data)) error stop 'f2k scalar: l_stability_correct_tau_zm'

    call c_f_pointer(input(37)%data, bp, [1])
    value%l_stability_correct_tau_zm = bp(1) /= 0_c_int8_t

    ! l_damp_wp2_using_em

    if (input(38)%rank /= 0) error stop 'f2k rank: l_damp_wp2_using_em'

    if (.not. c_associated(input(38)%data)) error stop 'f2k scalar: l_damp_wp2_using_em'

    call c_f_pointer(input(38)%data, bp, [1])
    value%l_damp_wp2_using_em = bp(1) /= 0_c_int8_t

    ! l_do_expldiff_rtm_thlm

    if (input(39)%rank /= 0) error stop 'f2k rank: l_do_expldiff_rtm_thlm'

    if (.not. c_associated(input(39)%data)) error stop 'f2k scalar: l_do_expldiff_rtm_thlm'

    call c_f_pointer(input(39)%data, bp, [1])
    value%l_do_expldiff_rtm_thlm = bp(1) /= 0_c_int8_t

    ! l_Lscale_plume_centered

    if (input(40)%rank /= 0) error stop 'f2k rank: l_Lscale_plume_centered'

    if (.not. c_associated(input(40)%data)) error stop 'f2k scalar: l_Lscale_plume_centered'

    call c_f_pointer(input(40)%data, bp, [1])
    value%l_Lscale_plume_centered = bp(1) /= 0_c_int8_t

    ! l_use_C7_Richardson

    if (input(41)%rank /= 0) error stop 'f2k rank: l_use_C7_Richardson'

    if (.not. c_associated(input(41)%data)) error stop 'f2k scalar: l_use_C7_Richardson'

    call c_f_pointer(input(41)%data, bp, [1])
    value%l_use_C7_Richardson = bp(1) /= 0_c_int8_t

    ! l_use_C11_Richardson

    if (input(42)%rank /= 0) error stop 'f2k rank: l_use_C11_Richardson'

    if (.not. c_associated(input(42)%data)) error stop 'f2k scalar: l_use_C11_Richardson'

    call c_f_pointer(input(42)%data, bp, [1])
    value%l_use_C11_Richardson = bp(1) /= 0_c_int8_t

    ! l_use_shear_Richardson

    if (input(43)%rank /= 0) error stop 'f2k rank: l_use_shear_Richardson'

    if (.not. c_associated(input(43)%data)) error stop 'f2k scalar: l_use_shear_Richardson'

    call c_f_pointer(input(43)%data, bp, [1])
    value%l_use_shear_Richardson = bp(1) /= 0_c_int8_t

    ! l_brunt_vaisala_freq_moist

    if (input(44)%rank /= 0) error stop 'f2k rank: l_brunt_vaisala_freq_moist'

    if (.not. c_associated(input(44)%data)) error stop 'f2k scalar: l_brunt_vaisala_freq_moist'

    call c_f_pointer(input(44)%data, bp, [1])
    value%l_brunt_vaisala_freq_moist = bp(1) /= 0_c_int8_t

    ! l_use_thvm_in_bv_freq

    if (input(45)%rank /= 0) error stop 'f2k rank: l_use_thvm_in_bv_freq'

    if (.not. c_associated(input(45)%data)) error stop 'f2k scalar: l_use_thvm_in_bv_freq'

    call c_f_pointer(input(45)%data, bp, [1])
    value%l_use_thvm_in_bv_freq = bp(1) /= 0_c_int8_t

    ! l_rcm_supersat_adj

    if (input(46)%rank /= 0) error stop 'f2k rank: l_rcm_supersat_adj'

    if (.not. c_associated(input(46)%data)) error stop 'f2k scalar: l_rcm_supersat_adj'

    call c_f_pointer(input(46)%data, bp, [1])
    value%l_rcm_supersat_adj = bp(1) /= 0_c_int8_t

    ! l_damp_wp3_Skw_squared

    if (input(47)%rank /= 0) error stop 'f2k rank: l_damp_wp3_Skw_squared'

    if (.not. c_associated(input(47)%data)) error stop 'f2k scalar: l_damp_wp3_Skw_squared'

    call c_f_pointer(input(47)%data, bp, [1])
    value%l_damp_wp3_Skw_squared = bp(1) /= 0_c_int8_t

    ! l_prescribed_avg_deltaz

    if (input(48)%rank /= 0) error stop 'f2k rank: l_prescribed_avg_deltaz'

    if (.not. c_associated(input(48)%data)) error stop 'f2k scalar: l_prescribed_avg_deltaz'

    call c_f_pointer(input(48)%data, bp, [1])
    value%l_prescribed_avg_deltaz = bp(1) /= 0_c_int8_t

    ! l_lmm_stepping

    if (input(49)%rank /= 0) error stop 'f2k rank: l_lmm_stepping'

    if (.not. c_associated(input(49)%data)) error stop 'f2k scalar: l_lmm_stepping'

    call c_f_pointer(input(49)%data, bp, [1])
    value%l_lmm_stepping = bp(1) /= 0_c_int8_t

    ! l_e3sm_config

    if (input(50)%rank /= 0) error stop 'f2k rank: l_e3sm_config'

    if (.not. c_associated(input(50)%data)) error stop 'f2k scalar: l_e3sm_config'

    call c_f_pointer(input(50)%data, bp, [1])
    value%l_e3sm_config = bp(1) /= 0_c_int8_t

    ! l_vary_convect_depth

    if (input(51)%rank /= 0) error stop 'f2k rank: l_vary_convect_depth'

    if (.not. c_associated(input(51)%data)) error stop 'f2k scalar: l_vary_convect_depth'

    call c_f_pointer(input(51)%data, bp, [1])
    value%l_vary_convect_depth = bp(1) /= 0_c_int8_t

    ! l_use_tke_in_wp3_pr_turb_term

    if (input(52)%rank /= 0) error stop 'f2k rank: l_use_tke_in_wp3_pr_turb_term'

    if (.not. c_associated(input(52)%data)) error stop 'f2k scalar: l_use_tke_in_wp3_pr_turb_term'

    call c_f_pointer(input(52)%data, bp, [1])
    value%l_use_tke_in_wp3_pr_turb_term = bp(1) /= 0_c_int8_t

    ! l_use_tke_in_wp2_wp3_K_dfsn

    if (input(53)%rank /= 0) error stop 'f2k rank: l_use_tke_in_wp2_wp3_K_dfsn'

    if (.not. c_associated(input(53)%data)) error stop 'f2k scalar: l_use_tke_in_wp2_wp3_K_dfsn'

    call c_f_pointer(input(53)%data, bp, [1])
    value%l_use_tke_in_wp2_wp3_K_dfsn = bp(1) /= 0_c_int8_t

    ! l_use_wp3_lim_with_smth_Heaviside

    if (input(54)%rank /= 0) error stop 'f2k rank: l_use_wp3_lim_with_smth_Heaviside'

    if (.not. c_associated(input(54)%data)) error stop 'f2k scalar: l_use_wp3_lim_with_smth_Heaviside'

    call c_f_pointer(input(54)%data, bp, [1])
    value%l_use_wp3_lim_with_smth_Heaviside = bp(1) /= 0_c_int8_t

    ! l_smooth_Heaviside_tau_wpxp

    if (input(55)%rank /= 0) error stop 'f2k rank: l_smooth_Heaviside_tau_wpxp'

    if (.not. c_associated(input(55)%data)) error stop 'f2k scalar: l_smooth_Heaviside_tau_wpxp'

    call c_f_pointer(input(55)%data, bp, [1])
    value%l_smooth_Heaviside_tau_wpxp = bp(1) /= 0_c_int8_t

    ! l_modify_limiters_for_cnvg_test

    if (input(56)%rank /= 0) error stop 'f2k rank: l_modify_limiters_for_cnvg_test'

    if (.not. c_associated(input(56)%data)) error stop 'f2k scalar: l_modify_limiters_for_cnvg_test'

    call c_f_pointer(input(56)%data, bp, [1])
    value%l_modify_limiters_for_cnvg_test = bp(1) /= 0_c_int8_t

    ! l_enable_relaxed_clipping

    if (input(57)%rank /= 0) error stop 'f2k rank: l_enable_relaxed_clipping'

    if (.not. c_associated(input(57)%data)) error stop 'f2k scalar: l_enable_relaxed_clipping'

    call c_f_pointer(input(57)%data, bp, [1])
    value%l_enable_relaxed_clipping = bp(1) /= 0_c_int8_t

    ! l_linearize_pbl_winds

    if (input(58)%rank /= 0) error stop 'f2k rank: l_linearize_pbl_winds'

    if (.not. c_associated(input(58)%data)) error stop 'f2k scalar: l_linearize_pbl_winds'

    call c_f_pointer(input(58)%data, bp, [1])
    value%l_linearize_pbl_winds = bp(1) /= 0_c_int8_t

    ! l_mono_flux_lim_thlm

    if (input(59)%rank /= 0) error stop 'f2k rank: l_mono_flux_lim_thlm'

    if (.not. c_associated(input(59)%data)) error stop 'f2k scalar: l_mono_flux_lim_thlm'

    call c_f_pointer(input(59)%data, bp, [1])
    value%l_mono_flux_lim_thlm = bp(1) /= 0_c_int8_t

    ! l_mono_flux_lim_rtm

    if (input(60)%rank /= 0) error stop 'f2k rank: l_mono_flux_lim_rtm'

    if (.not. c_associated(input(60)%data)) error stop 'f2k scalar: l_mono_flux_lim_rtm'

    call c_f_pointer(input(60)%data, bp, [1])
    value%l_mono_flux_lim_rtm = bp(1) /= 0_c_int8_t

    ! l_mono_flux_lim_um

    if (input(61)%rank /= 0) error stop 'f2k rank: l_mono_flux_lim_um'

    if (.not. c_associated(input(61)%data)) error stop 'f2k scalar: l_mono_flux_lim_um'

    call c_f_pointer(input(61)%data, bp, [1])
    value%l_mono_flux_lim_um = bp(1) /= 0_c_int8_t

    ! l_mono_flux_lim_vm

    if (input(62)%rank /= 0) error stop 'f2k rank: l_mono_flux_lim_vm'

    if (.not. c_associated(input(62)%data)) error stop 'f2k scalar: l_mono_flux_lim_vm'

    call c_f_pointer(input(62)%data, bp, [1])
    value%l_mono_flux_lim_vm = bp(1) /= 0_c_int8_t

    ! l_mono_flux_lim_spikefix

    if (input(63)%rank /= 0) error stop 'f2k rank: l_mono_flux_lim_spikefix'

    if (.not. c_associated(input(63)%data)) error stop 'f2k scalar: l_mono_flux_lim_spikefix'

    call c_f_pointer(input(63)%data, bp, [1])
    value%l_mono_flux_lim_spikefix = bp(1) /= 0_c_int8_t

    ! l_host_applies_sfc_fluxes

    if (input(64)%rank /= 0) error stop 'f2k rank: l_host_applies_sfc_fluxes'

    if (.not. c_associated(input(64)%data)) error stop 'f2k scalar: l_host_applies_sfc_fluxes'

    call c_f_pointer(input(64)%data, bp, [1])
    value%l_host_applies_sfc_fluxes = bp(1) /= 0_c_int8_t

    ! l_wp2_fill_holes_tke

    if (input(65)%rank /= 0) error stop 'f2k rank: l_wp2_fill_holes_tke'

    if (.not. c_associated(input(65)%data)) error stop 'f2k scalar: l_wp2_fill_holes_tke'

    call c_f_pointer(input(65)%data, bp, [1])
    value%l_wp2_fill_holes_tke = bp(1) /= 0_c_int8_t

    ! l_add_dycore_grid

    if (input(66)%rank /= 0) error stop 'f2k rank: l_add_dycore_grid'

    if (.not. c_associated(input(66)%data)) error stop 'f2k scalar: l_add_dycore_grid'

    call c_f_pointer(input(66)%data, bp, [1])
    value%l_add_dycore_grid = bp(1) /= 0_c_int8_t

  end subroutine to_native_clubb_config_flags_type

  subroutine from_native_clubb_config_flags_type(value, output)
    ! Description: Copy each clubb_config_flags_type field into C-owned output bytes.

    implicit none

    !--------------------------- Input Variables ---------------------------
    ! Call-local native argument; units and shape follow the source contract
    type(clubb_config_flags_type), intent(in) :: value

    !--------------------------- Output Variables ---------------------------
    type(wire_field), intent(out) :: output(*) ! Returned descriptors; C++ releases their C-owned buffers [-]

    !--------------------------- Local Variables ---------------------------
    integer(c_int), pointer :: ip(:) ! Borrowed view of packed integer data [-]
    integer(c_int8_t), pointer :: bp(:) ! Borrowed view of packed logical bytes [-]

    !--------------------------- Begin Code ---------------------------

    ! iiPDF_type

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(1)) ! Out

    call c_f_pointer(output(1)%data, ip, [1])
    ip(1) = value%iiPDF_type

    ! ipdf_call_placement

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(2)) ! Out

    call c_f_pointer(output(2)%data, ip, [1])
    ip(1) = value%ipdf_call_placement

    ! penta_solve_method

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(3)) ! Out

    call c_f_pointer(output(3)%data, ip, [1])
    ip(1) = value%penta_solve_method

    ! tridiag_solve_method

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(4)) ! Out

    call c_f_pointer(output(4)%data, ip, [1])
    ip(1) = value%tridiag_solve_method

    ! saturation_formula

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(5)) ! Out

    call c_f_pointer(output(5)%data, ip, [1])
    ip(1) = value%saturation_formula

    ! Lscale_method

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(6)) ! Out

    call c_f_pointer(output(6)%data, ip, [1])
    ip(1) = value%Lscale_method

    ! grid_remap_method

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(7)) ! Out

    call c_f_pointer(output(7)%data, ip, [1])
    ip(1) = value%grid_remap_method

    ! grid_adapt_in_time_method

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(8)) ! Out

    call c_f_pointer(output(8)%data, ip, [1])
    ip(1) = value%grid_adapt_in_time_method

    ! fill_holes_type

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(9)) ! Out

    call c_f_pointer(output(9)%data, ip, [1])
    ip(1) = value%fill_holes_type

    ! l_use_precip_frac

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(10)) ! Out

    call c_f_pointer(output(10)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_use_precip_frac)

    ! l_predict_upwp_vpwp

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(11)) ! Out

    call c_f_pointer(output(11)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_predict_upwp_vpwp)

    ! l_ho_nontrad_coriolis

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(12)) ! Out

    call c_f_pointer(output(12)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_ho_nontrad_coriolis)

    ! l_ho_trad_coriolis

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(13)) ! Out

    call c_f_pointer(output(13)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_ho_trad_coriolis)

    ! l_min_wp2_from_corr_wx

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(14)) ! Out

    call c_f_pointer(output(14)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_min_wp2_from_corr_wx)

    ! l_min_xp2_from_corr_wx

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(15)) ! Out

    call c_f_pointer(output(15)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_min_xp2_from_corr_wx)

    ! l_C2_cloud_frac

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(16)) ! Out

    call c_f_pointer(output(16)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_C2_cloud_frac)

    ! l_diffuse_rtm_and_thlm

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(17)) ! Out

    call c_f_pointer(output(17)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_diffuse_rtm_and_thlm)

    ! l_stability_correct_Kh_N2_zm

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(18)) ! Out

    call c_f_pointer(output(18)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_stability_correct_Kh_N2_zm)

    ! l_calc_thlp2_rad

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(19)) ! Out

    call c_f_pointer(output(19)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_calc_thlp2_rad)

    ! l_upwind_xpyp_ta

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(20)) ! Out

    call c_f_pointer(output(20)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_upwind_xpyp_ta)

    ! l_upwind_xm_ma

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(21)) ! Out

    call c_f_pointer(output(21)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_upwind_xm_ma)

    ! l_uv_nudge

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(22)) ! Out

    call c_f_pointer(output(22)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_uv_nudge)

    ! l_rtm_nudge

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(23)) ! Out

    call c_f_pointer(output(23)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_rtm_nudge)

    ! l_vert_avg_closure

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(24)) ! Out

    call c_f_pointer(output(24)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_vert_avg_closure)

    ! l_trapezoidal_rule_zt

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(25)) ! Out

    call c_f_pointer(output(25)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_trapezoidal_rule_zt)

    ! l_trapezoidal_rule_zm

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(26)) ! Out

    call c_f_pointer(output(26)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_trapezoidal_rule_zm)

    ! l_call_pdf_closure_twice

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(27)) ! Out

    call c_f_pointer(output(27)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_call_pdf_closure_twice)

    ! l_standard_term_ta

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(28)) ! Out

    call c_f_pointer(output(28)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_standard_term_ta)

    ! l_partial_upwind_wp3

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(29)) ! Out

    call c_f_pointer(output(29)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_partial_upwind_wp3)

    ! l_godunov_upwind_wpxp_ta

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(30)) ! Out

    call c_f_pointer(output(30)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_godunov_upwind_wpxp_ta)

    ! l_godunov_upwind_xpyp_ta

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(31)) ! Out

    call c_f_pointer(output(31)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_godunov_upwind_xpyp_ta)

    ! l_use_cloud_cover

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(32)) ! Out

    call c_f_pointer(output(32)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_use_cloud_cover)

    ! l_diagnose_correlations

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(33)) ! Out

    call c_f_pointer(output(33)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_diagnose_correlations)

    ! l_calc_w_corr

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(34)) ! Out

    call c_f_pointer(output(34)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_calc_w_corr)

    ! l_const_Nc_in_cloud

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(35)) ! Out

    call c_f_pointer(output(35)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_const_Nc_in_cloud)

    ! l_fix_w_chi_eta_correlations

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(36)) ! Out

    call c_f_pointer(output(36)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_fix_w_chi_eta_correlations)

    ! l_stability_correct_tau_zm

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(37)) ! Out

    call c_f_pointer(output(37)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_stability_correct_tau_zm)

    ! l_damp_wp2_using_em

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(38)) ! Out

    call c_f_pointer(output(38)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_damp_wp2_using_em)

    ! l_do_expldiff_rtm_thlm

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(39)) ! Out

    call c_f_pointer(output(39)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_do_expldiff_rtm_thlm)

    ! l_Lscale_plume_centered

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(40)) ! Out

    call c_f_pointer(output(40)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_Lscale_plume_centered)

    ! l_use_C7_Richardson

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(41)) ! Out

    call c_f_pointer(output(41)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_use_C7_Richardson)

    ! l_use_C11_Richardson

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(42)) ! Out

    call c_f_pointer(output(42)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_use_C11_Richardson)

    ! l_use_shear_Richardson

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(43)) ! Out

    call c_f_pointer(output(43)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_use_shear_Richardson)

    ! l_brunt_vaisala_freq_moist

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(44)) ! Out

    call c_f_pointer(output(44)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_brunt_vaisala_freq_moist)

    ! l_use_thvm_in_bv_freq

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(45)) ! Out

    call c_f_pointer(output(45)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_use_thvm_in_bv_freq)

    ! l_rcm_supersat_adj

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(46)) ! Out

    call c_f_pointer(output(46)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_rcm_supersat_adj)

    ! l_damp_wp3_Skw_squared

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(47)) ! Out

    call c_f_pointer(output(47)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_damp_wp3_Skw_squared)

    ! l_prescribed_avg_deltaz

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(48)) ! Out

    call c_f_pointer(output(48)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_prescribed_avg_deltaz)

    ! l_lmm_stepping

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(49)) ! Out

    call c_f_pointer(output(49)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_lmm_stepping)

    ! l_e3sm_config

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(50)) ! Out

    call c_f_pointer(output(50)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_e3sm_config)

    ! l_vary_convect_depth

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(51)) ! Out

    call c_f_pointer(output(51)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_vary_convect_depth)

    ! l_use_tke_in_wp3_pr_turb_term

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(52)) ! Out

    call c_f_pointer(output(52)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_use_tke_in_wp3_pr_turb_term)

    ! l_use_tke_in_wp2_wp3_K_dfsn

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(53)) ! Out

    call c_f_pointer(output(53)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_use_tke_in_wp2_wp3_K_dfsn)

    ! l_use_wp3_lim_with_smth_Heaviside

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(54)) ! Out

    call c_f_pointer(output(54)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_use_wp3_lim_with_smth_Heaviside)

    ! l_smooth_Heaviside_tau_wpxp

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(55)) ! Out

    call c_f_pointer(output(55)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_smooth_Heaviside_tau_wpxp)

    ! l_modify_limiters_for_cnvg_test

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(56)) ! Out

    call c_f_pointer(output(56)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_modify_limiters_for_cnvg_test)

    ! l_enable_relaxed_clipping

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(57)) ! Out

    call c_f_pointer(output(57)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_enable_relaxed_clipping)

    ! l_linearize_pbl_winds

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(58)) ! Out

    call c_f_pointer(output(58)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_linearize_pbl_winds)

    ! l_mono_flux_lim_thlm

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(59)) ! Out

    call c_f_pointer(output(59)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_mono_flux_lim_thlm)

    ! l_mono_flux_lim_rtm

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(60)) ! Out

    call c_f_pointer(output(60)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_mono_flux_lim_rtm)

    ! l_mono_flux_lim_um

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(61)) ! Out

    call c_f_pointer(output(61)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_mono_flux_lim_um)

    ! l_mono_flux_lim_vm

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(62)) ! Out

    call c_f_pointer(output(62)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_mono_flux_lim_vm)

    ! l_mono_flux_lim_spikefix

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(63)) ! Out

    call c_f_pointer(output(63)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_mono_flux_lim_spikefix)

    ! l_host_applies_sfc_fluxes

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(64)) ! Out

    call c_f_pointer(output(64)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_host_applies_sfc_fluxes)

    ! l_wp2_fill_holes_tke

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(65)) ! Out

    call c_f_pointer(output(65)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_wp2_fill_holes_tke)

    ! l_add_dycore_grid

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(66)) ! Out

    call c_f_pointer(output(66)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_add_dycore_grid)

  end subroutine from_native_clubb_config_flags_type

  subroutine f2k_roundtrip_clubb_config_flags_type_c(input, output) bind(C)
    ! Description: Exercise both converters with one local Fortran value.

    implicit none

    !--------------------------- Input Variables ---------------------------
    type(wire_field), intent(in) :: input(*) ! Caller-owned input field descriptors [-]

    !--------------------------- Output Variables ---------------------------
    type(wire_field), intent(out) :: output(*) ! Returned descriptors; C++ releases their C-owned buffers [-]

    !--------------------------- Local Variables ---------------------------
    ! Call-local native argument; units and shape follow the source contract
    type(clubb_config_flags_type) :: value

    !--------------------------- Begin Code ---------------------------

    call to_native_clubb_config_flags_type(input, value)
    call from_native_clubb_config_flags_type(value, output)

  end subroutine f2k_roundtrip_clubb_config_flags_type_c

  subroutine to_native_nu_vertical_res_dep(input, value)
    ! Description: Copy each C field into a local nu_vertical_res_dep value.

    implicit none

    !--------------------------- Input Variables ---------------------------
    type(wire_field), intent(in) :: input(*) ! Caller-owned input field descriptors [-]

    !--------------------------- Output Variables ---------------------------
    ! Call-local native argument; units and shape follow the source contract
    type(nu_vertical_res_dep), intent(out) :: value

    !--------------------------- Local Variables ---------------------------
    ! Borrowed view of packed real data; units follow the native field
    real( kind = core_rknd ), pointer :: rp(:)
    integer :: n ! Number of packed elements [-]

    !--------------------------- Begin Code ---------------------------

    ! nu1

    if (input(1)%rank /= 1) error stop 'f2k rank: nu1'

    if (input(1)%allocated /= 0) then

      allocate(value%nu1(int(input(1)%extent(1))))
      n = int(product(input(1)%extent(1:1)))

      if (n > 0) then

        if (.not. c_associated(input(1)%data)) error stop 'f2k data: nu1'

        call c_f_pointer(input(1)%data, rp, [n])
        value%nu1 = reshape(rp, shape(value%nu1))

      end if

    end if

    ! nu2

    if (input(2)%rank /= 1) error stop 'f2k rank: nu2'

    if (input(2)%allocated /= 0) then

      allocate(value%nu2(int(input(2)%extent(1))))
      n = int(product(input(2)%extent(1:1)))

      if (n > 0) then

        if (.not. c_associated(input(2)%data)) error stop 'f2k data: nu2'

        call c_f_pointer(input(2)%data, rp, [n])
        value%nu2 = reshape(rp, shape(value%nu2))

      end if

    end if

    ! nu6

    if (input(3)%rank /= 1) error stop 'f2k rank: nu6'

    if (input(3)%allocated /= 0) then

      allocate(value%nu6(int(input(3)%extent(1))))
      n = int(product(input(3)%extent(1:1)))

      if (n > 0) then

        if (.not. c_associated(input(3)%data)) error stop 'f2k data: nu6'

        call c_f_pointer(input(3)%data, rp, [n])
        value%nu6 = reshape(rp, shape(value%nu6))

      end if

    end if

    ! nu8

    if (input(4)%rank /= 1) error stop 'f2k rank: nu8'

    if (input(4)%allocated /= 0) then

      allocate(value%nu8(int(input(4)%extent(1))))
      n = int(product(input(4)%extent(1:1)))

      if (n > 0) then

        if (.not. c_associated(input(4)%data)) error stop 'f2k data: nu8'

        call c_f_pointer(input(4)%data, rp, [n])
        value%nu8 = reshape(rp, shape(value%nu8))

      end if

    end if

    ! nu9

    if (input(5)%rank /= 1) error stop 'f2k rank: nu9'

    if (input(5)%allocated /= 0) then

      allocate(value%nu9(int(input(5)%extent(1))))
      n = int(product(input(5)%extent(1:1)))

      if (n > 0) then

        if (.not. c_associated(input(5)%data)) error stop 'f2k data: nu9'

        call c_f_pointer(input(5)%data, rp, [n])
        value%nu9 = reshape(rp, shape(value%nu9))

      end if

    end if

    ! nu10

    if (input(6)%rank /= 1) error stop 'f2k rank: nu10'

    if (input(6)%allocated /= 0) then

      allocate(value%nu10(int(input(6)%extent(1))))
      n = int(product(input(6)%extent(1:1)))

      if (n > 0) then

        if (.not. c_associated(input(6)%data)) error stop 'f2k data: nu10'

        call c_f_pointer(input(6)%data, rp, [n])
        value%nu10 = reshape(rp, shape(value%nu10))

      end if

    end if

    ! nu_hm

    if (input(7)%rank /= 1) error stop 'f2k rank: nu_hm'

    if (input(7)%allocated /= 0) then

      allocate(value%nu_hm(int(input(7)%extent(1))))
      n = int(product(input(7)%extent(1:1)))

      if (n > 0) then

        if (.not. c_associated(input(7)%data)) error stop 'f2k data: nu_hm'

        call c_f_pointer(input(7)%data, rp, [n])
        value%nu_hm = reshape(rp, shape(value%nu_hm))

      end if

    end if

  end subroutine to_native_nu_vertical_res_dep

  subroutine from_native_nu_vertical_res_dep(value, output)
    ! Description: Copy each nu_vertical_res_dep field into C-owned output bytes.

    implicit none

    !--------------------------- Input Variables ---------------------------
    ! Call-local native argument; units and shape follow the source contract
    type(nu_vertical_res_dep), intent(in) :: value

    !--------------------------- Output Variables ---------------------------
    type(wire_field), intent(out) :: output(*) ! Returned descriptors; C++ releases their C-owned buffers [-]

    !--------------------------- Local Variables ---------------------------
    ! Borrowed view of packed real data; units follow the native field
    real( kind = core_rknd ), pointer :: rp(:)
    integer :: n ! Number of packed elements [-]
    integer(c_int64_t) :: field_shape(4) ! Expected field extents in native dimension order [-]

    !--------------------------- Begin Code ---------------------------

    ! nu1

    if (allocated(value%nu1)) then

      n = size(value%nu1)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%nu1, 1, kind=c_int64_t)

      call allocate_field(1, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(1)) ! Out

      call c_f_pointer(output(1)%data, rp, [n])
      rp = reshape(value%nu1, [n])

    else

      call empty_field(1, 0, & ! In
                       output(1)) ! Out

    end if

    ! nu2

    if (allocated(value%nu2)) then

      n = size(value%nu2)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%nu2, 1, kind=c_int64_t)

      call allocate_field(1, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(2)) ! Out

      call c_f_pointer(output(2)%data, rp, [n])
      rp = reshape(value%nu2, [n])

    else

      call empty_field(1, 0, & ! In
                       output(2)) ! Out

    end if

    ! nu6

    if (allocated(value%nu6)) then

      n = size(value%nu6)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%nu6, 1, kind=c_int64_t)

      call allocate_field(1, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(3)) ! Out

      call c_f_pointer(output(3)%data, rp, [n])
      rp = reshape(value%nu6, [n])

    else

      call empty_field(1, 0, & ! In
                       output(3)) ! Out

    end if

    ! nu8

    if (allocated(value%nu8)) then

      n = size(value%nu8)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%nu8, 1, kind=c_int64_t)

      call allocate_field(1, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(4)) ! Out

      call c_f_pointer(output(4)%data, rp, [n])
      rp = reshape(value%nu8, [n])

    else

      call empty_field(1, 0, & ! In
                       output(4)) ! Out

    end if

    ! nu9

    if (allocated(value%nu9)) then

      n = size(value%nu9)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%nu9, 1, kind=c_int64_t)

      call allocate_field(1, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(5)) ! Out

      call c_f_pointer(output(5)%data, rp, [n])
      rp = reshape(value%nu9, [n])

    else

      call empty_field(1, 0, & ! In
                       output(5)) ! Out

    end if

    ! nu10

    if (allocated(value%nu10)) then

      n = size(value%nu10)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%nu10, 1, kind=c_int64_t)

      call allocate_field(1, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(6)) ! Out

      call c_f_pointer(output(6)%data, rp, [n])
      rp = reshape(value%nu10, [n])

    else

      call empty_field(1, 0, & ! In
                       output(6)) ! Out

    end if

    ! nu_hm

    if (allocated(value%nu_hm)) then

      n = size(value%nu_hm)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%nu_hm, 1, kind=c_int64_t)

      call allocate_field(1, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(7)) ! Out

      call c_f_pointer(output(7)%data, rp, [n])
      rp = reshape(value%nu_hm, [n])

    else

      call empty_field(1, 0, & ! In
                       output(7)) ! Out

    end if

  end subroutine from_native_nu_vertical_res_dep

  subroutine f2k_roundtrip_nu_vertical_res_dep_c(input, output) bind(C)
    ! Description: Exercise both converters with one local Fortran value.

    implicit none

    !--------------------------- Input Variables ---------------------------
    type(wire_field), intent(in) :: input(*) ! Caller-owned input field descriptors [-]

    !--------------------------- Output Variables ---------------------------
    type(wire_field), intent(out) :: output(*) ! Returned descriptors; C++ releases their C-owned buffers [-]

    !--------------------------- Local Variables ---------------------------
    ! Call-local native argument; units and shape follow the source contract
    type(nu_vertical_res_dep) :: value

    !--------------------------- Begin Code ---------------------------

    call to_native_nu_vertical_res_dep(input, value)
    call from_native_nu_vertical_res_dep(value, output)

  end subroutine f2k_roundtrip_nu_vertical_res_dep_c

  subroutine to_native_pdf_parameter(input, value)
    ! Description: Copy each C field into a local pdf_parameter value.

    implicit none

    !--------------------------- Input Variables ---------------------------
    type(wire_field), intent(in) :: input(*) ! Caller-owned input field descriptors [-]

    !--------------------------- Output Variables ---------------------------
    ! Call-local native argument; units and shape follow the source contract
    type(pdf_parameter), intent(out) :: value

    !--------------------------- Local Variables ---------------------------
    integer(c_int), pointer :: ip(:) ! Borrowed view of packed integer data [-]
    ! Borrowed view of packed real data; units follow the native field
    real( kind = core_rknd ), pointer :: rp(:)
    integer :: n ! Number of packed elements [-]

    !--------------------------- Begin Code ---------------------------

    ! ngrdcol

    if (input(1)%rank /= 0) error stop 'f2k rank: ngrdcol'

    if (.not. c_associated(input(1)%data)) error stop 'f2k scalar: ngrdcol'

    call c_f_pointer(input(1)%data, ip, [1])
    value%ngrdcol = ip(1)

    ! nz

    if (input(2)%rank /= 0) error stop 'f2k rank: nz'

    if (.not. c_associated(input(2)%data)) error stop 'f2k scalar: nz'

    call c_f_pointer(input(2)%data, ip, [1])
    value%nz = ip(1)

    ! w_1

    if (input(3)%rank /= 2) error stop 'f2k rank: w_1'

    if (input(3)%allocated /= 0) then

      allocate(value%w_1(int(input(3)%extent(1)), int(input(3)%extent(2))))
      n = int(product(input(3)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(3)%data)) error stop 'f2k data: w_1'

        call c_f_pointer(input(3)%data, rp, [n])
        value%w_1 = reshape(rp, shape(value%w_1))

      end if

    end if

    ! w_2

    if (input(4)%rank /= 2) error stop 'f2k rank: w_2'

    if (input(4)%allocated /= 0) then

      allocate(value%w_2(int(input(4)%extent(1)), int(input(4)%extent(2))))
      n = int(product(input(4)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(4)%data)) error stop 'f2k data: w_2'

        call c_f_pointer(input(4)%data, rp, [n])
        value%w_2 = reshape(rp, shape(value%w_2))

      end if

    end if

    ! varnce_w_1

    if (input(5)%rank /= 2) error stop 'f2k rank: varnce_w_1'

    if (input(5)%allocated /= 0) then

      allocate(value%varnce_w_1(int(input(5)%extent(1)), int(input(5)%extent(2))))
      n = int(product(input(5)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(5)%data)) error stop 'f2k data: varnce_w_1'

        call c_f_pointer(input(5)%data, rp, [n])
        value%varnce_w_1 = reshape(rp, shape(value%varnce_w_1))

      end if

    end if

    ! varnce_w_2

    if (input(6)%rank /= 2) error stop 'f2k rank: varnce_w_2'

    if (input(6)%allocated /= 0) then

      allocate(value%varnce_w_2(int(input(6)%extent(1)), int(input(6)%extent(2))))
      n = int(product(input(6)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(6)%data)) error stop 'f2k data: varnce_w_2'

        call c_f_pointer(input(6)%data, rp, [n])
        value%varnce_w_2 = reshape(rp, shape(value%varnce_w_2))

      end if

    end if

    ! rt_1

    if (input(7)%rank /= 2) error stop 'f2k rank: rt_1'

    if (input(7)%allocated /= 0) then

      allocate(value%rt_1(int(input(7)%extent(1)), int(input(7)%extent(2))))
      n = int(product(input(7)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(7)%data)) error stop 'f2k data: rt_1'

        call c_f_pointer(input(7)%data, rp, [n])
        value%rt_1 = reshape(rp, shape(value%rt_1))

      end if

    end if

    ! rt_2

    if (input(8)%rank /= 2) error stop 'f2k rank: rt_2'

    if (input(8)%allocated /= 0) then

      allocate(value%rt_2(int(input(8)%extent(1)), int(input(8)%extent(2))))
      n = int(product(input(8)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(8)%data)) error stop 'f2k data: rt_2'

        call c_f_pointer(input(8)%data, rp, [n])
        value%rt_2 = reshape(rp, shape(value%rt_2))

      end if

    end if

    ! varnce_rt_1

    if (input(9)%rank /= 2) error stop 'f2k rank: varnce_rt_1'

    if (input(9)%allocated /= 0) then

      allocate(value%varnce_rt_1(int(input(9)%extent(1)), int(input(9)%extent(2))))
      n = int(product(input(9)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(9)%data)) error stop 'f2k data: varnce_rt_1'

        call c_f_pointer(input(9)%data, rp, [n])
        value%varnce_rt_1 = reshape(rp, shape(value%varnce_rt_1))

      end if

    end if

    ! varnce_rt_2

    if (input(10)%rank /= 2) error stop 'f2k rank: varnce_rt_2'

    if (input(10)%allocated /= 0) then

      allocate(value%varnce_rt_2(int(input(10)%extent(1)), int(input(10)%extent(2))))
      n = int(product(input(10)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(10)%data)) error stop 'f2k data: varnce_rt_2'

        call c_f_pointer(input(10)%data, rp, [n])
        value%varnce_rt_2 = reshape(rp, shape(value%varnce_rt_2))

      end if

    end if

    ! thl_1

    if (input(11)%rank /= 2) error stop 'f2k rank: thl_1'

    if (input(11)%allocated /= 0) then

      allocate(value%thl_1(int(input(11)%extent(1)), int(input(11)%extent(2))))
      n = int(product(input(11)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(11)%data)) error stop 'f2k data: thl_1'

        call c_f_pointer(input(11)%data, rp, [n])
        value%thl_1 = reshape(rp, shape(value%thl_1))

      end if

    end if

    ! thl_2

    if (input(12)%rank /= 2) error stop 'f2k rank: thl_2'

    if (input(12)%allocated /= 0) then

      allocate(value%thl_2(int(input(12)%extent(1)), int(input(12)%extent(2))))
      n = int(product(input(12)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(12)%data)) error stop 'f2k data: thl_2'

        call c_f_pointer(input(12)%data, rp, [n])
        value%thl_2 = reshape(rp, shape(value%thl_2))

      end if

    end if

    ! varnce_thl_1

    if (input(13)%rank /= 2) error stop 'f2k rank: varnce_thl_1'

    if (input(13)%allocated /= 0) then

      allocate(value%varnce_thl_1(int(input(13)%extent(1)), int(input(13)%extent(2))))
      n = int(product(input(13)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(13)%data)) error stop 'f2k data: varnce_thl_1'

        call c_f_pointer(input(13)%data, rp, [n])
        value%varnce_thl_1 = reshape(rp, shape(value%varnce_thl_1))

      end if

    end if

    ! varnce_thl_2

    if (input(14)%rank /= 2) error stop 'f2k rank: varnce_thl_2'

    if (input(14)%allocated /= 0) then

      allocate(value%varnce_thl_2(int(input(14)%extent(1)), int(input(14)%extent(2))))
      n = int(product(input(14)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(14)%data)) error stop 'f2k data: varnce_thl_2'

        call c_f_pointer(input(14)%data, rp, [n])
        value%varnce_thl_2 = reshape(rp, shape(value%varnce_thl_2))

      end if

    end if

    ! corr_w_rt_1

    if (input(15)%rank /= 2) error stop 'f2k rank: corr_w_rt_1'

    if (input(15)%allocated /= 0) then

      allocate(value%corr_w_rt_1(int(input(15)%extent(1)), int(input(15)%extent(2))))
      n = int(product(input(15)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(15)%data)) error stop 'f2k data: corr_w_rt_1'

        call c_f_pointer(input(15)%data, rp, [n])
        value%corr_w_rt_1 = reshape(rp, shape(value%corr_w_rt_1))

      end if

    end if

    ! corr_w_rt_2

    if (input(16)%rank /= 2) error stop 'f2k rank: corr_w_rt_2'

    if (input(16)%allocated /= 0) then

      allocate(value%corr_w_rt_2(int(input(16)%extent(1)), int(input(16)%extent(2))))
      n = int(product(input(16)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(16)%data)) error stop 'f2k data: corr_w_rt_2'

        call c_f_pointer(input(16)%data, rp, [n])
        value%corr_w_rt_2 = reshape(rp, shape(value%corr_w_rt_2))

      end if

    end if

    ! corr_w_thl_1

    if (input(17)%rank /= 2) error stop 'f2k rank: corr_w_thl_1'

    if (input(17)%allocated /= 0) then

      allocate(value%corr_w_thl_1(int(input(17)%extent(1)), int(input(17)%extent(2))))
      n = int(product(input(17)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(17)%data)) error stop 'f2k data: corr_w_thl_1'

        call c_f_pointer(input(17)%data, rp, [n])
        value%corr_w_thl_1 = reshape(rp, shape(value%corr_w_thl_1))

      end if

    end if

    ! corr_w_thl_2

    if (input(18)%rank /= 2) error stop 'f2k rank: corr_w_thl_2'

    if (input(18)%allocated /= 0) then

      allocate(value%corr_w_thl_2(int(input(18)%extent(1)), int(input(18)%extent(2))))
      n = int(product(input(18)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(18)%data)) error stop 'f2k data: corr_w_thl_2'

        call c_f_pointer(input(18)%data, rp, [n])
        value%corr_w_thl_2 = reshape(rp, shape(value%corr_w_thl_2))

      end if

    end if

    ! corr_rt_thl_1

    if (input(19)%rank /= 2) error stop 'f2k rank: corr_rt_thl_1'

    if (input(19)%allocated /= 0) then

      allocate(value%corr_rt_thl_1(int(input(19)%extent(1)), int(input(19)%extent(2))))
      n = int(product(input(19)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(19)%data)) error stop 'f2k data: corr_rt_thl_1'

        call c_f_pointer(input(19)%data, rp, [n])
        value%corr_rt_thl_1 = reshape(rp, shape(value%corr_rt_thl_1))

      end if

    end if

    ! corr_rt_thl_2

    if (input(20)%rank /= 2) error stop 'f2k rank: corr_rt_thl_2'

    if (input(20)%allocated /= 0) then

      allocate(value%corr_rt_thl_2(int(input(20)%extent(1)), int(input(20)%extent(2))))
      n = int(product(input(20)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(20)%data)) error stop 'f2k data: corr_rt_thl_2'

        call c_f_pointer(input(20)%data, rp, [n])
        value%corr_rt_thl_2 = reshape(rp, shape(value%corr_rt_thl_2))

      end if

    end if

    ! alpha_thl

    if (input(21)%rank /= 2) error stop 'f2k rank: alpha_thl'

    if (input(21)%allocated /= 0) then

      allocate(value%alpha_thl(int(input(21)%extent(1)), int(input(21)%extent(2))))
      n = int(product(input(21)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(21)%data)) error stop 'f2k data: alpha_thl'

        call c_f_pointer(input(21)%data, rp, [n])
        value%alpha_thl = reshape(rp, shape(value%alpha_thl))

      end if

    end if

    ! alpha_rt

    if (input(22)%rank /= 2) error stop 'f2k rank: alpha_rt'

    if (input(22)%allocated /= 0) then

      allocate(value%alpha_rt(int(input(22)%extent(1)), int(input(22)%extent(2))))
      n = int(product(input(22)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(22)%data)) error stop 'f2k data: alpha_rt'

        call c_f_pointer(input(22)%data, rp, [n])
        value%alpha_rt = reshape(rp, shape(value%alpha_rt))

      end if

    end if

    ! crt_1

    if (input(23)%rank /= 2) error stop 'f2k rank: crt_1'

    if (input(23)%allocated /= 0) then

      allocate(value%crt_1(int(input(23)%extent(1)), int(input(23)%extent(2))))
      n = int(product(input(23)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(23)%data)) error stop 'f2k data: crt_1'

        call c_f_pointer(input(23)%data, rp, [n])
        value%crt_1 = reshape(rp, shape(value%crt_1))

      end if

    end if

    ! crt_2

    if (input(24)%rank /= 2) error stop 'f2k rank: crt_2'

    if (input(24)%allocated /= 0) then

      allocate(value%crt_2(int(input(24)%extent(1)), int(input(24)%extent(2))))
      n = int(product(input(24)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(24)%data)) error stop 'f2k data: crt_2'

        call c_f_pointer(input(24)%data, rp, [n])
        value%crt_2 = reshape(rp, shape(value%crt_2))

      end if

    end if

    ! cthl_1

    if (input(25)%rank /= 2) error stop 'f2k rank: cthl_1'

    if (input(25)%allocated /= 0) then

      allocate(value%cthl_1(int(input(25)%extent(1)), int(input(25)%extent(2))))
      n = int(product(input(25)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(25)%data)) error stop 'f2k data: cthl_1'

        call c_f_pointer(input(25)%data, rp, [n])
        value%cthl_1 = reshape(rp, shape(value%cthl_1))

      end if

    end if

    ! cthl_2

    if (input(26)%rank /= 2) error stop 'f2k rank: cthl_2'

    if (input(26)%allocated /= 0) then

      allocate(value%cthl_2(int(input(26)%extent(1)), int(input(26)%extent(2))))
      n = int(product(input(26)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(26)%data)) error stop 'f2k data: cthl_2'

        call c_f_pointer(input(26)%data, rp, [n])
        value%cthl_2 = reshape(rp, shape(value%cthl_2))

      end if

    end if

    ! chi_1

    if (input(27)%rank /= 2) error stop 'f2k rank: chi_1'

    if (input(27)%allocated /= 0) then

      allocate(value%chi_1(int(input(27)%extent(1)), int(input(27)%extent(2))))
      n = int(product(input(27)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(27)%data)) error stop 'f2k data: chi_1'

        call c_f_pointer(input(27)%data, rp, [n])
        value%chi_1 = reshape(rp, shape(value%chi_1))

      end if

    end if

    ! chi_2

    if (input(28)%rank /= 2) error stop 'f2k rank: chi_2'

    if (input(28)%allocated /= 0) then

      allocate(value%chi_2(int(input(28)%extent(1)), int(input(28)%extent(2))))
      n = int(product(input(28)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(28)%data)) error stop 'f2k data: chi_2'

        call c_f_pointer(input(28)%data, rp, [n])
        value%chi_2 = reshape(rp, shape(value%chi_2))

      end if

    end if

    ! stdev_chi_1

    if (input(29)%rank /= 2) error stop 'f2k rank: stdev_chi_1'

    if (input(29)%allocated /= 0) then

      allocate(value%stdev_chi_1(int(input(29)%extent(1)), int(input(29)%extent(2))))
      n = int(product(input(29)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(29)%data)) error stop 'f2k data: stdev_chi_1'

        call c_f_pointer(input(29)%data, rp, [n])
        value%stdev_chi_1 = reshape(rp, shape(value%stdev_chi_1))

      end if

    end if

    ! stdev_chi_2

    if (input(30)%rank /= 2) error stop 'f2k rank: stdev_chi_2'

    if (input(30)%allocated /= 0) then

      allocate(value%stdev_chi_2(int(input(30)%extent(1)), int(input(30)%extent(2))))
      n = int(product(input(30)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(30)%data)) error stop 'f2k data: stdev_chi_2'

        call c_f_pointer(input(30)%data, rp, [n])
        value%stdev_chi_2 = reshape(rp, shape(value%stdev_chi_2))

      end if

    end if

    ! stdev_eta_1

    if (input(31)%rank /= 2) error stop 'f2k rank: stdev_eta_1'

    if (input(31)%allocated /= 0) then

      allocate(value%stdev_eta_1(int(input(31)%extent(1)), int(input(31)%extent(2))))
      n = int(product(input(31)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(31)%data)) error stop 'f2k data: stdev_eta_1'

        call c_f_pointer(input(31)%data, rp, [n])
        value%stdev_eta_1 = reshape(rp, shape(value%stdev_eta_1))

      end if

    end if

    ! stdev_eta_2

    if (input(32)%rank /= 2) error stop 'f2k rank: stdev_eta_2'

    if (input(32)%allocated /= 0) then

      allocate(value%stdev_eta_2(int(input(32)%extent(1)), int(input(32)%extent(2))))
      n = int(product(input(32)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(32)%data)) error stop 'f2k data: stdev_eta_2'

        call c_f_pointer(input(32)%data, rp, [n])
        value%stdev_eta_2 = reshape(rp, shape(value%stdev_eta_2))

      end if

    end if

    ! covar_chi_eta_1

    if (input(33)%rank /= 2) error stop 'f2k rank: covar_chi_eta_1'

    if (input(33)%allocated /= 0) then

      allocate(value%covar_chi_eta_1(int(input(33)%extent(1)), int(input(33)%extent(2))))
      n = int(product(input(33)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(33)%data)) error stop 'f2k data: covar_chi_eta_1'

        call c_f_pointer(input(33)%data, rp, [n])
        value%covar_chi_eta_1 = reshape(rp, shape(value%covar_chi_eta_1))

      end if

    end if

    ! covar_chi_eta_2

    if (input(34)%rank /= 2) error stop 'f2k rank: covar_chi_eta_2'

    if (input(34)%allocated /= 0) then

      allocate(value%covar_chi_eta_2(int(input(34)%extent(1)), int(input(34)%extent(2))))
      n = int(product(input(34)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(34)%data)) error stop 'f2k data: covar_chi_eta_2'

        call c_f_pointer(input(34)%data, rp, [n])
        value%covar_chi_eta_2 = reshape(rp, shape(value%covar_chi_eta_2))

      end if

    end if

    ! corr_w_chi_1

    if (input(35)%rank /= 2) error stop 'f2k rank: corr_w_chi_1'

    if (input(35)%allocated /= 0) then

      allocate(value%corr_w_chi_1(int(input(35)%extent(1)), int(input(35)%extent(2))))
      n = int(product(input(35)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(35)%data)) error stop 'f2k data: corr_w_chi_1'

        call c_f_pointer(input(35)%data, rp, [n])
        value%corr_w_chi_1 = reshape(rp, shape(value%corr_w_chi_1))

      end if

    end if

    ! corr_w_chi_2

    if (input(36)%rank /= 2) error stop 'f2k rank: corr_w_chi_2'

    if (input(36)%allocated /= 0) then

      allocate(value%corr_w_chi_2(int(input(36)%extent(1)), int(input(36)%extent(2))))
      n = int(product(input(36)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(36)%data)) error stop 'f2k data: corr_w_chi_2'

        call c_f_pointer(input(36)%data, rp, [n])
        value%corr_w_chi_2 = reshape(rp, shape(value%corr_w_chi_2))

      end if

    end if

    ! corr_w_eta_1

    if (input(37)%rank /= 2) error stop 'f2k rank: corr_w_eta_1'

    if (input(37)%allocated /= 0) then

      allocate(value%corr_w_eta_1(int(input(37)%extent(1)), int(input(37)%extent(2))))
      n = int(product(input(37)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(37)%data)) error stop 'f2k data: corr_w_eta_1'

        call c_f_pointer(input(37)%data, rp, [n])
        value%corr_w_eta_1 = reshape(rp, shape(value%corr_w_eta_1))

      end if

    end if

    ! corr_w_eta_2

    if (input(38)%rank /= 2) error stop 'f2k rank: corr_w_eta_2'

    if (input(38)%allocated /= 0) then

      allocate(value%corr_w_eta_2(int(input(38)%extent(1)), int(input(38)%extent(2))))
      n = int(product(input(38)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(38)%data)) error stop 'f2k data: corr_w_eta_2'

        call c_f_pointer(input(38)%data, rp, [n])
        value%corr_w_eta_2 = reshape(rp, shape(value%corr_w_eta_2))

      end if

    end if

    ! corr_chi_eta_1

    if (input(39)%rank /= 2) error stop 'f2k rank: corr_chi_eta_1'

    if (input(39)%allocated /= 0) then

      allocate(value%corr_chi_eta_1(int(input(39)%extent(1)), int(input(39)%extent(2))))
      n = int(product(input(39)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(39)%data)) error stop 'f2k data: corr_chi_eta_1'

        call c_f_pointer(input(39)%data, rp, [n])
        value%corr_chi_eta_1 = reshape(rp, shape(value%corr_chi_eta_1))

      end if

    end if

    ! corr_chi_eta_2

    if (input(40)%rank /= 2) error stop 'f2k rank: corr_chi_eta_2'

    if (input(40)%allocated /= 0) then

      allocate(value%corr_chi_eta_2(int(input(40)%extent(1)), int(input(40)%extent(2))))
      n = int(product(input(40)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(40)%data)) error stop 'f2k data: corr_chi_eta_2'

        call c_f_pointer(input(40)%data, rp, [n])
        value%corr_chi_eta_2 = reshape(rp, shape(value%corr_chi_eta_2))

      end if

    end if

    ! rsatl_1

    if (input(41)%rank /= 2) error stop 'f2k rank: rsatl_1'

    if (input(41)%allocated /= 0) then

      allocate(value%rsatl_1(int(input(41)%extent(1)), int(input(41)%extent(2))))
      n = int(product(input(41)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(41)%data)) error stop 'f2k data: rsatl_1'

        call c_f_pointer(input(41)%data, rp, [n])
        value%rsatl_1 = reshape(rp, shape(value%rsatl_1))

      end if

    end if

    ! rsatl_2

    if (input(42)%rank /= 2) error stop 'f2k rank: rsatl_2'

    if (input(42)%allocated /= 0) then

      allocate(value%rsatl_2(int(input(42)%extent(1)), int(input(42)%extent(2))))
      n = int(product(input(42)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(42)%data)) error stop 'f2k data: rsatl_2'

        call c_f_pointer(input(42)%data, rp, [n])
        value%rsatl_2 = reshape(rp, shape(value%rsatl_2))

      end if

    end if

    ! rc_1

    if (input(43)%rank /= 2) error stop 'f2k rank: rc_1'

    if (input(43)%allocated /= 0) then

      allocate(value%rc_1(int(input(43)%extent(1)), int(input(43)%extent(2))))
      n = int(product(input(43)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(43)%data)) error stop 'f2k data: rc_1'

        call c_f_pointer(input(43)%data, rp, [n])
        value%rc_1 = reshape(rp, shape(value%rc_1))

      end if

    end if

    ! rc_2

    if (input(44)%rank /= 2) error stop 'f2k rank: rc_2'

    if (input(44)%allocated /= 0) then

      allocate(value%rc_2(int(input(44)%extent(1)), int(input(44)%extent(2))))
      n = int(product(input(44)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(44)%data)) error stop 'f2k data: rc_2'

        call c_f_pointer(input(44)%data, rp, [n])
        value%rc_2 = reshape(rp, shape(value%rc_2))

      end if

    end if

    ! cloud_frac_1

    if (input(45)%rank /= 2) error stop 'f2k rank: cloud_frac_1'

    if (input(45)%allocated /= 0) then

      allocate(value%cloud_frac_1(int(input(45)%extent(1)), int(input(45)%extent(2))))
      n = int(product(input(45)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(45)%data)) error stop 'f2k data: cloud_frac_1'

        call c_f_pointer(input(45)%data, rp, [n])
        value%cloud_frac_1 = reshape(rp, shape(value%cloud_frac_1))

      end if

    end if

    ! cloud_frac_2

    if (input(46)%rank /= 2) error stop 'f2k rank: cloud_frac_2'

    if (input(46)%allocated /= 0) then

      allocate(value%cloud_frac_2(int(input(46)%extent(1)), int(input(46)%extent(2))))
      n = int(product(input(46)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(46)%data)) error stop 'f2k data: cloud_frac_2'

        call c_f_pointer(input(46)%data, rp, [n])
        value%cloud_frac_2 = reshape(rp, shape(value%cloud_frac_2))

      end if

    end if

    ! mixt_frac

    if (input(47)%rank /= 2) error stop 'f2k rank: mixt_frac'

    if (input(47)%allocated /= 0) then

      allocate(value%mixt_frac(int(input(47)%extent(1)), int(input(47)%extent(2))))
      n = int(product(input(47)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(47)%data)) error stop 'f2k data: mixt_frac'

        call c_f_pointer(input(47)%data, rp, [n])
        value%mixt_frac = reshape(rp, shape(value%mixt_frac))

      end if

    end if

    ! ice_supersat_frac_1

    if (input(48)%rank /= 2) error stop 'f2k rank: ice_supersat_frac_1'

    if (input(48)%allocated /= 0) then

      allocate(value%ice_supersat_frac_1(int(input(48)%extent(1)), int(input(48)%extent(2))))
      n = int(product(input(48)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(48)%data)) error stop 'f2k data: ice_supersat_frac_1'

        call c_f_pointer(input(48)%data, rp, [n])
        value%ice_supersat_frac_1 = reshape(rp, shape(value%ice_supersat_frac_1))

      end if

    end if

    ! ice_supersat_frac_2

    if (input(49)%rank /= 2) error stop 'f2k rank: ice_supersat_frac_2'

    if (input(49)%allocated /= 0) then

      allocate(value%ice_supersat_frac_2(int(input(49)%extent(1)), int(input(49)%extent(2))))
      n = int(product(input(49)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(49)%data)) error stop 'f2k data: ice_supersat_frac_2'

        call c_f_pointer(input(49)%data, rp, [n])
        value%ice_supersat_frac_2 = reshape(rp, shape(value%ice_supersat_frac_2))

      end if

    end if

  end subroutine to_native_pdf_parameter

  subroutine from_native_pdf_parameter(value, output)
    ! Description: Copy each pdf_parameter field into C-owned output bytes.

    implicit none

    !--------------------------- Input Variables ---------------------------
    ! Call-local native argument; units and shape follow the source contract
    type(pdf_parameter), intent(in) :: value

    !--------------------------- Output Variables ---------------------------
    type(wire_field), intent(out) :: output(*) ! Returned descriptors; C++ releases their C-owned buffers [-]

    !--------------------------- Local Variables ---------------------------
    integer(c_int), pointer :: ip(:) ! Borrowed view of packed integer data [-]
    ! Borrowed view of packed real data; units follow the native field
    real( kind = core_rknd ), pointer :: rp(:)
    integer :: n ! Number of packed elements [-]
    integer(c_int64_t) :: field_shape(4) ! Expected field extents in native dimension order [-]

    !--------------------------- Begin Code ---------------------------

    ! ngrdcol

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(1)) ! Out

    call c_f_pointer(output(1)%data, ip, [1])
    ip(1) = value%ngrdcol

    ! nz

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(2)) ! Out

    call c_f_pointer(output(2)%data, ip, [1])
    ip(1) = value%nz

    ! w_1

    if (allocated(value%w_1)) then

      n = size(value%w_1)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%w_1, 1, kind=c_int64_t)
      field_shape(2) = size(value%w_1, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(3)) ! Out

      call c_f_pointer(output(3)%data, rp, [n])
      rp = reshape(value%w_1, [n])

    else

      call empty_field(2, 0, & ! In
                       output(3)) ! Out

    end if

    ! w_2

    if (allocated(value%w_2)) then

      n = size(value%w_2)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%w_2, 1, kind=c_int64_t)
      field_shape(2) = size(value%w_2, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(4)) ! Out

      call c_f_pointer(output(4)%data, rp, [n])
      rp = reshape(value%w_2, [n])

    else

      call empty_field(2, 0, & ! In
                       output(4)) ! Out

    end if

    ! varnce_w_1

    if (allocated(value%varnce_w_1)) then

      n = size(value%varnce_w_1)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%varnce_w_1, 1, kind=c_int64_t)
      field_shape(2) = size(value%varnce_w_1, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(5)) ! Out

      call c_f_pointer(output(5)%data, rp, [n])
      rp = reshape(value%varnce_w_1, [n])

    else

      call empty_field(2, 0, & ! In
                       output(5)) ! Out

    end if

    ! varnce_w_2

    if (allocated(value%varnce_w_2)) then

      n = size(value%varnce_w_2)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%varnce_w_2, 1, kind=c_int64_t)
      field_shape(2) = size(value%varnce_w_2, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(6)) ! Out

      call c_f_pointer(output(6)%data, rp, [n])
      rp = reshape(value%varnce_w_2, [n])

    else

      call empty_field(2, 0, & ! In
                       output(6)) ! Out

    end if

    ! rt_1

    if (allocated(value%rt_1)) then

      n = size(value%rt_1)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%rt_1, 1, kind=c_int64_t)
      field_shape(2) = size(value%rt_1, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(7)) ! Out

      call c_f_pointer(output(7)%data, rp, [n])
      rp = reshape(value%rt_1, [n])

    else

      call empty_field(2, 0, & ! In
                       output(7)) ! Out

    end if

    ! rt_2

    if (allocated(value%rt_2)) then

      n = size(value%rt_2)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%rt_2, 1, kind=c_int64_t)
      field_shape(2) = size(value%rt_2, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(8)) ! Out

      call c_f_pointer(output(8)%data, rp, [n])
      rp = reshape(value%rt_2, [n])

    else

      call empty_field(2, 0, & ! In
                       output(8)) ! Out

    end if

    ! varnce_rt_1

    if (allocated(value%varnce_rt_1)) then

      n = size(value%varnce_rt_1)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%varnce_rt_1, 1, kind=c_int64_t)
      field_shape(2) = size(value%varnce_rt_1, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(9)) ! Out

      call c_f_pointer(output(9)%data, rp, [n])
      rp = reshape(value%varnce_rt_1, [n])

    else

      call empty_field(2, 0, & ! In
                       output(9)) ! Out

    end if

    ! varnce_rt_2

    if (allocated(value%varnce_rt_2)) then

      n = size(value%varnce_rt_2)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%varnce_rt_2, 1, kind=c_int64_t)
      field_shape(2) = size(value%varnce_rt_2, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(10)) ! Out

      call c_f_pointer(output(10)%data, rp, [n])
      rp = reshape(value%varnce_rt_2, [n])

    else

      call empty_field(2, 0, & ! In
                       output(10)) ! Out

    end if

    ! thl_1

    if (allocated(value%thl_1)) then

      n = size(value%thl_1)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%thl_1, 1, kind=c_int64_t)
      field_shape(2) = size(value%thl_1, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(11)) ! Out

      call c_f_pointer(output(11)%data, rp, [n])
      rp = reshape(value%thl_1, [n])

    else

      call empty_field(2, 0, & ! In
                       output(11)) ! Out

    end if

    ! thl_2

    if (allocated(value%thl_2)) then

      n = size(value%thl_2)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%thl_2, 1, kind=c_int64_t)
      field_shape(2) = size(value%thl_2, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(12)) ! Out

      call c_f_pointer(output(12)%data, rp, [n])
      rp = reshape(value%thl_2, [n])

    else

      call empty_field(2, 0, & ! In
                       output(12)) ! Out

    end if

    ! varnce_thl_1

    if (allocated(value%varnce_thl_1)) then

      n = size(value%varnce_thl_1)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%varnce_thl_1, 1, kind=c_int64_t)
      field_shape(2) = size(value%varnce_thl_1, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(13)) ! Out

      call c_f_pointer(output(13)%data, rp, [n])
      rp = reshape(value%varnce_thl_1, [n])

    else

      call empty_field(2, 0, & ! In
                       output(13)) ! Out

    end if

    ! varnce_thl_2

    if (allocated(value%varnce_thl_2)) then

      n = size(value%varnce_thl_2)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%varnce_thl_2, 1, kind=c_int64_t)
      field_shape(2) = size(value%varnce_thl_2, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(14)) ! Out

      call c_f_pointer(output(14)%data, rp, [n])
      rp = reshape(value%varnce_thl_2, [n])

    else

      call empty_field(2, 0, & ! In
                       output(14)) ! Out

    end if

    ! corr_w_rt_1

    if (allocated(value%corr_w_rt_1)) then

      n = size(value%corr_w_rt_1)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%corr_w_rt_1, 1, kind=c_int64_t)
      field_shape(2) = size(value%corr_w_rt_1, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(15)) ! Out

      call c_f_pointer(output(15)%data, rp, [n])
      rp = reshape(value%corr_w_rt_1, [n])

    else

      call empty_field(2, 0, & ! In
                       output(15)) ! Out

    end if

    ! corr_w_rt_2

    if (allocated(value%corr_w_rt_2)) then

      n = size(value%corr_w_rt_2)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%corr_w_rt_2, 1, kind=c_int64_t)
      field_shape(2) = size(value%corr_w_rt_2, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(16)) ! Out

      call c_f_pointer(output(16)%data, rp, [n])
      rp = reshape(value%corr_w_rt_2, [n])

    else

      call empty_field(2, 0, & ! In
                       output(16)) ! Out

    end if

    ! corr_w_thl_1

    if (allocated(value%corr_w_thl_1)) then

      n = size(value%corr_w_thl_1)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%corr_w_thl_1, 1, kind=c_int64_t)
      field_shape(2) = size(value%corr_w_thl_1, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(17)) ! Out

      call c_f_pointer(output(17)%data, rp, [n])
      rp = reshape(value%corr_w_thl_1, [n])

    else

      call empty_field(2, 0, & ! In
                       output(17)) ! Out

    end if

    ! corr_w_thl_2

    if (allocated(value%corr_w_thl_2)) then

      n = size(value%corr_w_thl_2)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%corr_w_thl_2, 1, kind=c_int64_t)
      field_shape(2) = size(value%corr_w_thl_2, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(18)) ! Out

      call c_f_pointer(output(18)%data, rp, [n])
      rp = reshape(value%corr_w_thl_2, [n])

    else

      call empty_field(2, 0, & ! In
                       output(18)) ! Out

    end if

    ! corr_rt_thl_1

    if (allocated(value%corr_rt_thl_1)) then

      n = size(value%corr_rt_thl_1)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%corr_rt_thl_1, 1, kind=c_int64_t)
      field_shape(2) = size(value%corr_rt_thl_1, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(19)) ! Out

      call c_f_pointer(output(19)%data, rp, [n])
      rp = reshape(value%corr_rt_thl_1, [n])

    else

      call empty_field(2, 0, & ! In
                       output(19)) ! Out

    end if

    ! corr_rt_thl_2

    if (allocated(value%corr_rt_thl_2)) then

      n = size(value%corr_rt_thl_2)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%corr_rt_thl_2, 1, kind=c_int64_t)
      field_shape(2) = size(value%corr_rt_thl_2, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(20)) ! Out

      call c_f_pointer(output(20)%data, rp, [n])
      rp = reshape(value%corr_rt_thl_2, [n])

    else

      call empty_field(2, 0, & ! In
                       output(20)) ! Out

    end if

    ! alpha_thl

    if (allocated(value%alpha_thl)) then

      n = size(value%alpha_thl)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%alpha_thl, 1, kind=c_int64_t)
      field_shape(2) = size(value%alpha_thl, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(21)) ! Out

      call c_f_pointer(output(21)%data, rp, [n])
      rp = reshape(value%alpha_thl, [n])

    else

      call empty_field(2, 0, & ! In
                       output(21)) ! Out

    end if

    ! alpha_rt

    if (allocated(value%alpha_rt)) then

      n = size(value%alpha_rt)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%alpha_rt, 1, kind=c_int64_t)
      field_shape(2) = size(value%alpha_rt, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(22)) ! Out

      call c_f_pointer(output(22)%data, rp, [n])
      rp = reshape(value%alpha_rt, [n])

    else

      call empty_field(2, 0, & ! In
                       output(22)) ! Out

    end if

    ! crt_1

    if (allocated(value%crt_1)) then

      n = size(value%crt_1)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%crt_1, 1, kind=c_int64_t)
      field_shape(2) = size(value%crt_1, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(23)) ! Out

      call c_f_pointer(output(23)%data, rp, [n])
      rp = reshape(value%crt_1, [n])

    else

      call empty_field(2, 0, & ! In
                       output(23)) ! Out

    end if

    ! crt_2

    if (allocated(value%crt_2)) then

      n = size(value%crt_2)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%crt_2, 1, kind=c_int64_t)
      field_shape(2) = size(value%crt_2, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(24)) ! Out

      call c_f_pointer(output(24)%data, rp, [n])
      rp = reshape(value%crt_2, [n])

    else

      call empty_field(2, 0, & ! In
                       output(24)) ! Out

    end if

    ! cthl_1

    if (allocated(value%cthl_1)) then

      n = size(value%cthl_1)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%cthl_1, 1, kind=c_int64_t)
      field_shape(2) = size(value%cthl_1, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(25)) ! Out

      call c_f_pointer(output(25)%data, rp, [n])
      rp = reshape(value%cthl_1, [n])

    else

      call empty_field(2, 0, & ! In
                       output(25)) ! Out

    end if

    ! cthl_2

    if (allocated(value%cthl_2)) then

      n = size(value%cthl_2)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%cthl_2, 1, kind=c_int64_t)
      field_shape(2) = size(value%cthl_2, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(26)) ! Out

      call c_f_pointer(output(26)%data, rp, [n])
      rp = reshape(value%cthl_2, [n])

    else

      call empty_field(2, 0, & ! In
                       output(26)) ! Out

    end if

    ! chi_1

    if (allocated(value%chi_1)) then

      n = size(value%chi_1)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%chi_1, 1, kind=c_int64_t)
      field_shape(2) = size(value%chi_1, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(27)) ! Out

      call c_f_pointer(output(27)%data, rp, [n])
      rp = reshape(value%chi_1, [n])

    else

      call empty_field(2, 0, & ! In
                       output(27)) ! Out

    end if

    ! chi_2

    if (allocated(value%chi_2)) then

      n = size(value%chi_2)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%chi_2, 1, kind=c_int64_t)
      field_shape(2) = size(value%chi_2, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(28)) ! Out

      call c_f_pointer(output(28)%data, rp, [n])
      rp = reshape(value%chi_2, [n])

    else

      call empty_field(2, 0, & ! In
                       output(28)) ! Out

    end if

    ! stdev_chi_1

    if (allocated(value%stdev_chi_1)) then

      n = size(value%stdev_chi_1)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%stdev_chi_1, 1, kind=c_int64_t)
      field_shape(2) = size(value%stdev_chi_1, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(29)) ! Out

      call c_f_pointer(output(29)%data, rp, [n])
      rp = reshape(value%stdev_chi_1, [n])

    else

      call empty_field(2, 0, & ! In
                       output(29)) ! Out

    end if

    ! stdev_chi_2

    if (allocated(value%stdev_chi_2)) then

      n = size(value%stdev_chi_2)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%stdev_chi_2, 1, kind=c_int64_t)
      field_shape(2) = size(value%stdev_chi_2, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(30)) ! Out

      call c_f_pointer(output(30)%data, rp, [n])
      rp = reshape(value%stdev_chi_2, [n])

    else

      call empty_field(2, 0, & ! In
                       output(30)) ! Out

    end if

    ! stdev_eta_1

    if (allocated(value%stdev_eta_1)) then

      n = size(value%stdev_eta_1)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%stdev_eta_1, 1, kind=c_int64_t)
      field_shape(2) = size(value%stdev_eta_1, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(31)) ! Out

      call c_f_pointer(output(31)%data, rp, [n])
      rp = reshape(value%stdev_eta_1, [n])

    else

      call empty_field(2, 0, & ! In
                       output(31)) ! Out

    end if

    ! stdev_eta_2

    if (allocated(value%stdev_eta_2)) then

      n = size(value%stdev_eta_2)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%stdev_eta_2, 1, kind=c_int64_t)
      field_shape(2) = size(value%stdev_eta_2, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(32)) ! Out

      call c_f_pointer(output(32)%data, rp, [n])
      rp = reshape(value%stdev_eta_2, [n])

    else

      call empty_field(2, 0, & ! In
                       output(32)) ! Out

    end if

    ! covar_chi_eta_1

    if (allocated(value%covar_chi_eta_1)) then

      n = size(value%covar_chi_eta_1)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%covar_chi_eta_1, 1, kind=c_int64_t)
      field_shape(2) = size(value%covar_chi_eta_1, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(33)) ! Out

      call c_f_pointer(output(33)%data, rp, [n])
      rp = reshape(value%covar_chi_eta_1, [n])

    else

      call empty_field(2, 0, & ! In
                       output(33)) ! Out

    end if

    ! covar_chi_eta_2

    if (allocated(value%covar_chi_eta_2)) then

      n = size(value%covar_chi_eta_2)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%covar_chi_eta_2, 1, kind=c_int64_t)
      field_shape(2) = size(value%covar_chi_eta_2, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(34)) ! Out

      call c_f_pointer(output(34)%data, rp, [n])
      rp = reshape(value%covar_chi_eta_2, [n])

    else

      call empty_field(2, 0, & ! In
                       output(34)) ! Out

    end if

    ! corr_w_chi_1

    if (allocated(value%corr_w_chi_1)) then

      n = size(value%corr_w_chi_1)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%corr_w_chi_1, 1, kind=c_int64_t)
      field_shape(2) = size(value%corr_w_chi_1, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(35)) ! Out

      call c_f_pointer(output(35)%data, rp, [n])
      rp = reshape(value%corr_w_chi_1, [n])

    else

      call empty_field(2, 0, & ! In
                       output(35)) ! Out

    end if

    ! corr_w_chi_2

    if (allocated(value%corr_w_chi_2)) then

      n = size(value%corr_w_chi_2)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%corr_w_chi_2, 1, kind=c_int64_t)
      field_shape(2) = size(value%corr_w_chi_2, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(36)) ! Out

      call c_f_pointer(output(36)%data, rp, [n])
      rp = reshape(value%corr_w_chi_2, [n])

    else

      call empty_field(2, 0, & ! In
                       output(36)) ! Out

    end if

    ! corr_w_eta_1

    if (allocated(value%corr_w_eta_1)) then

      n = size(value%corr_w_eta_1)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%corr_w_eta_1, 1, kind=c_int64_t)
      field_shape(2) = size(value%corr_w_eta_1, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(37)) ! Out

      call c_f_pointer(output(37)%data, rp, [n])
      rp = reshape(value%corr_w_eta_1, [n])

    else

      call empty_field(2, 0, & ! In
                       output(37)) ! Out

    end if

    ! corr_w_eta_2

    if (allocated(value%corr_w_eta_2)) then

      n = size(value%corr_w_eta_2)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%corr_w_eta_2, 1, kind=c_int64_t)
      field_shape(2) = size(value%corr_w_eta_2, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(38)) ! Out

      call c_f_pointer(output(38)%data, rp, [n])
      rp = reshape(value%corr_w_eta_2, [n])

    else

      call empty_field(2, 0, & ! In
                       output(38)) ! Out

    end if

    ! corr_chi_eta_1

    if (allocated(value%corr_chi_eta_1)) then

      n = size(value%corr_chi_eta_1)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%corr_chi_eta_1, 1, kind=c_int64_t)
      field_shape(2) = size(value%corr_chi_eta_1, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(39)) ! Out

      call c_f_pointer(output(39)%data, rp, [n])
      rp = reshape(value%corr_chi_eta_1, [n])

    else

      call empty_field(2, 0, & ! In
                       output(39)) ! Out

    end if

    ! corr_chi_eta_2

    if (allocated(value%corr_chi_eta_2)) then

      n = size(value%corr_chi_eta_2)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%corr_chi_eta_2, 1, kind=c_int64_t)
      field_shape(2) = size(value%corr_chi_eta_2, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(40)) ! Out

      call c_f_pointer(output(40)%data, rp, [n])
      rp = reshape(value%corr_chi_eta_2, [n])

    else

      call empty_field(2, 0, & ! In
                       output(40)) ! Out

    end if

    ! rsatl_1

    if (allocated(value%rsatl_1)) then

      n = size(value%rsatl_1)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%rsatl_1, 1, kind=c_int64_t)
      field_shape(2) = size(value%rsatl_1, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(41)) ! Out

      call c_f_pointer(output(41)%data, rp, [n])
      rp = reshape(value%rsatl_1, [n])

    else

      call empty_field(2, 0, & ! In
                       output(41)) ! Out

    end if

    ! rsatl_2

    if (allocated(value%rsatl_2)) then

      n = size(value%rsatl_2)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%rsatl_2, 1, kind=c_int64_t)
      field_shape(2) = size(value%rsatl_2, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(42)) ! Out

      call c_f_pointer(output(42)%data, rp, [n])
      rp = reshape(value%rsatl_2, [n])

    else

      call empty_field(2, 0, & ! In
                       output(42)) ! Out

    end if

    ! rc_1

    if (allocated(value%rc_1)) then

      n = size(value%rc_1)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%rc_1, 1, kind=c_int64_t)
      field_shape(2) = size(value%rc_1, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(43)) ! Out

      call c_f_pointer(output(43)%data, rp, [n])
      rp = reshape(value%rc_1, [n])

    else

      call empty_field(2, 0, & ! In
                       output(43)) ! Out

    end if

    ! rc_2

    if (allocated(value%rc_2)) then

      n = size(value%rc_2)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%rc_2, 1, kind=c_int64_t)
      field_shape(2) = size(value%rc_2, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(44)) ! Out

      call c_f_pointer(output(44)%data, rp, [n])
      rp = reshape(value%rc_2, [n])

    else

      call empty_field(2, 0, & ! In
                       output(44)) ! Out

    end if

    ! cloud_frac_1

    if (allocated(value%cloud_frac_1)) then

      n = size(value%cloud_frac_1)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%cloud_frac_1, 1, kind=c_int64_t)
      field_shape(2) = size(value%cloud_frac_1, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(45)) ! Out

      call c_f_pointer(output(45)%data, rp, [n])
      rp = reshape(value%cloud_frac_1, [n])

    else

      call empty_field(2, 0, & ! In
                       output(45)) ! Out

    end if

    ! cloud_frac_2

    if (allocated(value%cloud_frac_2)) then

      n = size(value%cloud_frac_2)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%cloud_frac_2, 1, kind=c_int64_t)
      field_shape(2) = size(value%cloud_frac_2, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(46)) ! Out

      call c_f_pointer(output(46)%data, rp, [n])
      rp = reshape(value%cloud_frac_2, [n])

    else

      call empty_field(2, 0, & ! In
                       output(46)) ! Out

    end if

    ! mixt_frac

    if (allocated(value%mixt_frac)) then

      n = size(value%mixt_frac)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%mixt_frac, 1, kind=c_int64_t)
      field_shape(2) = size(value%mixt_frac, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(47)) ! Out

      call c_f_pointer(output(47)%data, rp, [n])
      rp = reshape(value%mixt_frac, [n])

    else

      call empty_field(2, 0, & ! In
                       output(47)) ! Out

    end if

    ! ice_supersat_frac_1

    if (allocated(value%ice_supersat_frac_1)) then

      n = size(value%ice_supersat_frac_1)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%ice_supersat_frac_1, 1, kind=c_int64_t)
      field_shape(2) = size(value%ice_supersat_frac_1, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(48)) ! Out

      call c_f_pointer(output(48)%data, rp, [n])
      rp = reshape(value%ice_supersat_frac_1, [n])

    else

      call empty_field(2, 0, & ! In
                       output(48)) ! Out

    end if

    ! ice_supersat_frac_2

    if (allocated(value%ice_supersat_frac_2)) then

      n = size(value%ice_supersat_frac_2)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%ice_supersat_frac_2, 1, kind=c_int64_t)
      field_shape(2) = size(value%ice_supersat_frac_2, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(49)) ! Out

      call c_f_pointer(output(49)%data, rp, [n])
      rp = reshape(value%ice_supersat_frac_2, [n])

    else

      call empty_field(2, 0, & ! In
                       output(49)) ! Out

    end if

  end subroutine from_native_pdf_parameter

  subroutine f2k_roundtrip_pdf_parameter_c(input, output) bind(C)
    ! Description: Exercise both converters with one local Fortran value.

    implicit none

    !--------------------------- Input Variables ---------------------------
    type(wire_field), intent(in) :: input(*) ! Caller-owned input field descriptors [-]

    !--------------------------- Output Variables ---------------------------
    type(wire_field), intent(out) :: output(*) ! Returned descriptors; C++ releases their C-owned buffers [-]

    !--------------------------- Local Variables ---------------------------
    type(pdf_parameter) :: value ! Call-local native argument; units and shape follow the source contract

    !--------------------------- Begin Code ---------------------------

    call to_native_pdf_parameter(input, value)
    call from_native_pdf_parameter(value, output)

  end subroutine f2k_roundtrip_pdf_parameter_c

  subroutine to_native_implicit_coefs_terms(input, value)
    ! Description: Copy each C field into a local implicit_coefs_terms value.

    implicit none

    !--------------------------- Input Variables ---------------------------
    type(wire_field), intent(in) :: input(*) ! Caller-owned input field descriptors [-]

    !--------------------------- Output Variables ---------------------------
    ! Call-local native argument; units and shape follow the source contract
    type(implicit_coefs_terms), intent(out) :: value

    !--------------------------- Local Variables ---------------------------
    integer(c_int), pointer :: ip(:) ! Borrowed view of packed integer data [-]
    ! Borrowed view of packed real data; units follow the native field
    real( kind = core_rknd ), pointer :: rp(:)
    integer :: n ! Number of packed elements [-]

    !--------------------------- Begin Code ---------------------------

    ! ngrdcol

    if (input(1)%rank /= 0) error stop 'f2k rank: ngrdcol'

    if (.not. c_associated(input(1)%data)) error stop 'f2k scalar: ngrdcol'

    call c_f_pointer(input(1)%data, ip, [1])
    value%ngrdcol = ip(1)

    ! nz

    if (input(2)%rank /= 0) error stop 'f2k rank: nz'

    if (.not. c_associated(input(2)%data)) error stop 'f2k scalar: nz'

    call c_f_pointer(input(2)%data, ip, [1])
    value%nz = ip(1)

    ! sclr_dim

    if (input(3)%rank /= 0) error stop 'f2k rank: sclr_dim'

    if (.not. c_associated(input(3)%data)) error stop 'f2k scalar: sclr_dim'

    call c_f_pointer(input(3)%data, ip, [1])
    value%sclr_dim = ip(1)

    ! coef_wp4_implicit

    if (input(4)%rank /= 2) error stop 'f2k rank: coef_wp4_implicit'

    if (input(4)%allocated /= 0) then

      allocate(value%coef_wp4_implicit(int(input(4)%extent(1)), int(input(4)%extent(2))))
      n = int(product(input(4)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(4)%data)) error stop 'f2k data: coef_wp4_implicit'

        call c_f_pointer(input(4)%data, rp, [n])
        value%coef_wp4_implicit = reshape(rp, shape(value%coef_wp4_implicit))

      end if

    end if

    ! coef_wp2rtp_implicit

    if (input(5)%rank /= 2) error stop 'f2k rank: coef_wp2rtp_implicit'

    if (input(5)%allocated /= 0) then

      allocate(value%coef_wp2rtp_implicit(int(input(5)%extent(1)), int(input(5)%extent(2))))
      n = int(product(input(5)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(5)%data)) error stop 'f2k data: coef_wp2rtp_implicit'

        call c_f_pointer(input(5)%data, rp, [n])
        value%coef_wp2rtp_implicit = reshape(rp, shape(value%coef_wp2rtp_implicit))

      end if

    end if

    ! term_wp2rtp_explicit

    if (input(6)%rank /= 2) error stop 'f2k rank: term_wp2rtp_explicit'

    if (input(6)%allocated /= 0) then

      allocate(value%term_wp2rtp_explicit(int(input(6)%extent(1)), int(input(6)%extent(2))))
      n = int(product(input(6)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(6)%data)) error stop 'f2k data: term_wp2rtp_explicit'

        call c_f_pointer(input(6)%data, rp, [n])
        value%term_wp2rtp_explicit = reshape(rp, shape(value%term_wp2rtp_explicit))

      end if

    end if

    ! coef_wp2thlp_implicit

    if (input(7)%rank /= 2) error stop 'f2k rank: coef_wp2thlp_implicit'

    if (input(7)%allocated /= 0) then

      allocate(value%coef_wp2thlp_implicit(int(input(7)%extent(1)), int(input(7)%extent(2))))
      n = int(product(input(7)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(7)%data)) error stop 'f2k data: coef_wp2thlp_implicit'

        call c_f_pointer(input(7)%data, rp, [n])
        value%coef_wp2thlp_implicit = reshape(rp, shape(value%coef_wp2thlp_implicit))

      end if

    end if

    ! term_wp2thlp_explicit

    if (input(8)%rank /= 2) error stop 'f2k rank: term_wp2thlp_explicit'

    if (input(8)%allocated /= 0) then

      allocate(value%term_wp2thlp_explicit(int(input(8)%extent(1)), int(input(8)%extent(2))))
      n = int(product(input(8)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(8)%data)) error stop 'f2k data: term_wp2thlp_explicit'

        call c_f_pointer(input(8)%data, rp, [n])
        value%term_wp2thlp_explicit = reshape(rp, shape(value%term_wp2thlp_explicit))

      end if

    end if

    ! coef_wp2up_implicit

    if (input(9)%rank /= 2) error stop 'f2k rank: coef_wp2up_implicit'

    if (input(9)%allocated /= 0) then

      allocate(value%coef_wp2up_implicit(int(input(9)%extent(1)), int(input(9)%extent(2))))
      n = int(product(input(9)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(9)%data)) error stop 'f2k data: coef_wp2up_implicit'

        call c_f_pointer(input(9)%data, rp, [n])
        value%coef_wp2up_implicit = reshape(rp, shape(value%coef_wp2up_implicit))

      end if

    end if

    ! term_wp2up_explicit

    if (input(10)%rank /= 2) error stop 'f2k rank: term_wp2up_explicit'

    if (input(10)%allocated /= 0) then

      allocate(value%term_wp2up_explicit(int(input(10)%extent(1)), int(input(10)%extent(2))))
      n = int(product(input(10)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(10)%data)) error stop 'f2k data: term_wp2up_explicit'

        call c_f_pointer(input(10)%data, rp, [n])
        value%term_wp2up_explicit = reshape(rp, shape(value%term_wp2up_explicit))

      end if

    end if

    ! coef_wp2vp_implicit

    if (input(11)%rank /= 2) error stop 'f2k rank: coef_wp2vp_implicit'

    if (input(11)%allocated /= 0) then

      allocate(value%coef_wp2vp_implicit(int(input(11)%extent(1)), int(input(11)%extent(2))))
      n = int(product(input(11)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(11)%data)) error stop 'f2k data: coef_wp2vp_implicit'

        call c_f_pointer(input(11)%data, rp, [n])
        value%coef_wp2vp_implicit = reshape(rp, shape(value%coef_wp2vp_implicit))

      end if

    end if

    ! term_wp2vp_explicit

    if (input(12)%rank /= 2) error stop 'f2k rank: term_wp2vp_explicit'

    if (input(12)%allocated /= 0) then

      allocate(value%term_wp2vp_explicit(int(input(12)%extent(1)), int(input(12)%extent(2))))
      n = int(product(input(12)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(12)%data)) error stop 'f2k data: term_wp2vp_explicit'

        call c_f_pointer(input(12)%data, rp, [n])
        value%term_wp2vp_explicit = reshape(rp, shape(value%term_wp2vp_explicit))

      end if

    end if

    ! coef_wprtp2_implicit

    if (input(13)%rank /= 2) error stop 'f2k rank: coef_wprtp2_implicit'

    if (input(13)%allocated /= 0) then

      allocate(value%coef_wprtp2_implicit(int(input(13)%extent(1)), int(input(13)%extent(2))))
      n = int(product(input(13)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(13)%data)) error stop 'f2k data: coef_wprtp2_implicit'

        call c_f_pointer(input(13)%data, rp, [n])
        value%coef_wprtp2_implicit = reshape(rp, shape(value%coef_wprtp2_implicit))

      end if

    end if

    ! term_wprtp2_explicit

    if (input(14)%rank /= 2) error stop 'f2k rank: term_wprtp2_explicit'

    if (input(14)%allocated /= 0) then

      allocate(value%term_wprtp2_explicit(int(input(14)%extent(1)), int(input(14)%extent(2))))
      n = int(product(input(14)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(14)%data)) error stop 'f2k data: term_wprtp2_explicit'

        call c_f_pointer(input(14)%data, rp, [n])
        value%term_wprtp2_explicit = reshape(rp, shape(value%term_wprtp2_explicit))

      end if

    end if

    ! coef_wpthlp2_implicit

    if (input(15)%rank /= 2) error stop 'f2k rank: coef_wpthlp2_implicit'

    if (input(15)%allocated /= 0) then

      allocate(value%coef_wpthlp2_implicit(int(input(15)%extent(1)), int(input(15)%extent(2))))
      n = int(product(input(15)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(15)%data)) error stop 'f2k data: coef_wpthlp2_implicit'

        call c_f_pointer(input(15)%data, rp, [n])
        value%coef_wpthlp2_implicit = reshape(rp, shape(value%coef_wpthlp2_implicit))

      end if

    end if

    ! term_wpthlp2_explicit

    if (input(16)%rank /= 2) error stop 'f2k rank: term_wpthlp2_explicit'

    if (input(16)%allocated /= 0) then

      allocate(value%term_wpthlp2_explicit(int(input(16)%extent(1)), int(input(16)%extent(2))))
      n = int(product(input(16)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(16)%data)) error stop 'f2k data: term_wpthlp2_explicit'

        call c_f_pointer(input(16)%data, rp, [n])
        value%term_wpthlp2_explicit = reshape(rp, shape(value%term_wpthlp2_explicit))

      end if

    end if

    ! coef_wprtpthlp_implicit

    if (input(17)%rank /= 2) error stop 'f2k rank: coef_wprtpthlp_implicit'

    if (input(17)%allocated /= 0) then

      allocate(value%coef_wprtpthlp_implicit(int(input(17)%extent(1)), int(input(17)%extent(2))))
      n = int(product(input(17)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(17)%data)) error stop 'f2k data: coef_wprtpthlp_implicit'

        call c_f_pointer(input(17)%data, rp, [n])
        value%coef_wprtpthlp_implicit = reshape(rp, shape(value%coef_wprtpthlp_implicit))

      end if

    end if

    ! term_wprtpthlp_explicit

    if (input(18)%rank /= 2) error stop 'f2k rank: term_wprtpthlp_explicit'

    if (input(18)%allocated /= 0) then

      allocate(value%term_wprtpthlp_explicit(int(input(18)%extent(1)), int(input(18)%extent(2))))
      n = int(product(input(18)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(18)%data)) error stop 'f2k data: term_wprtpthlp_explicit'

        call c_f_pointer(input(18)%data, rp, [n])
        value%term_wprtpthlp_explicit = reshape(rp, shape(value%term_wprtpthlp_explicit))

      end if

    end if

    ! coef_wpup2_implicit

    if (input(19)%rank /= 2) error stop 'f2k rank: coef_wpup2_implicit'

    if (input(19)%allocated /= 0) then

      allocate(value%coef_wpup2_implicit(int(input(19)%extent(1)), int(input(19)%extent(2))))
      n = int(product(input(19)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(19)%data)) error stop 'f2k data: coef_wpup2_implicit'

        call c_f_pointer(input(19)%data, rp, [n])
        value%coef_wpup2_implicit = reshape(rp, shape(value%coef_wpup2_implicit))

      end if

    end if

    ! term_wpup2_explicit

    if (input(20)%rank /= 2) error stop 'f2k rank: term_wpup2_explicit'

    if (input(20)%allocated /= 0) then

      allocate(value%term_wpup2_explicit(int(input(20)%extent(1)), int(input(20)%extent(2))))
      n = int(product(input(20)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(20)%data)) error stop 'f2k data: term_wpup2_explicit'

        call c_f_pointer(input(20)%data, rp, [n])
        value%term_wpup2_explicit = reshape(rp, shape(value%term_wpup2_explicit))

      end if

    end if

    ! coef_wpvp2_implicit

    if (input(21)%rank /= 2) error stop 'f2k rank: coef_wpvp2_implicit'

    if (input(21)%allocated /= 0) then

      allocate(value%coef_wpvp2_implicit(int(input(21)%extent(1)), int(input(21)%extent(2))))
      n = int(product(input(21)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(21)%data)) error stop 'f2k data: coef_wpvp2_implicit'

        call c_f_pointer(input(21)%data, rp, [n])
        value%coef_wpvp2_implicit = reshape(rp, shape(value%coef_wpvp2_implicit))

      end if

    end if

    ! term_wpvp2_explicit

    if (input(22)%rank /= 2) error stop 'f2k rank: term_wpvp2_explicit'

    if (input(22)%allocated /= 0) then

      allocate(value%term_wpvp2_explicit(int(input(22)%extent(1)), int(input(22)%extent(2))))
      n = int(product(input(22)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(22)%data)) error stop 'f2k data: term_wpvp2_explicit'

        call c_f_pointer(input(22)%data, rp, [n])
        value%term_wpvp2_explicit = reshape(rp, shape(value%term_wpvp2_explicit))

      end if

    end if

    ! coef_wp2sclrp_implicit

    if (input(23)%rank /= 3) error stop 'f2k rank: coef_wp2sclrp_implicit'

    if (input(23)%allocated /= 0) then

      allocate(value%coef_wp2sclrp_implicit(int(input(23)%extent(1)), int(input(23)%extent(2)), int(input(23)%extent(3))))
      n = int(product(input(23)%extent(1:3)))

      if (n > 0) then

        if (.not. c_associated(input(23)%data)) error stop 'f2k data: coef_wp2sclrp_implicit'

        call c_f_pointer(input(23)%data, rp, [n])
        value%coef_wp2sclrp_implicit = reshape(rp, shape(value%coef_wp2sclrp_implicit))

      end if

    end if

    ! term_wp2sclrp_explicit

    if (input(24)%rank /= 3) error stop 'f2k rank: term_wp2sclrp_explicit'

    if (input(24)%allocated /= 0) then

      allocate(value%term_wp2sclrp_explicit(int(input(24)%extent(1)), int(input(24)%extent(2)), int(input(24)%extent(3))))
      n = int(product(input(24)%extent(1:3)))

      if (n > 0) then

        if (.not. c_associated(input(24)%data)) error stop 'f2k data: term_wp2sclrp_explicit'

        call c_f_pointer(input(24)%data, rp, [n])
        value%term_wp2sclrp_explicit = reshape(rp, shape(value%term_wp2sclrp_explicit))

      end if

    end if

    ! coef_wpsclrp2_implicit

    if (input(25)%rank /= 3) error stop 'f2k rank: coef_wpsclrp2_implicit'

    if (input(25)%allocated /= 0) then

      allocate(value%coef_wpsclrp2_implicit(int(input(25)%extent(1)), int(input(25)%extent(2)), int(input(25)%extent(3))))
      n = int(product(input(25)%extent(1:3)))

      if (n > 0) then

        if (.not. c_associated(input(25)%data)) error stop 'f2k data: coef_wpsclrp2_implicit'

        call c_f_pointer(input(25)%data, rp, [n])
        value%coef_wpsclrp2_implicit = reshape(rp, shape(value%coef_wpsclrp2_implicit))

      end if

    end if

    ! term_wpsclrp2_explicit

    if (input(26)%rank /= 3) error stop 'f2k rank: term_wpsclrp2_explicit'

    if (input(26)%allocated /= 0) then

      allocate(value%term_wpsclrp2_explicit(int(input(26)%extent(1)), int(input(26)%extent(2)), int(input(26)%extent(3))))
      n = int(product(input(26)%extent(1:3)))

      if (n > 0) then

        if (.not. c_associated(input(26)%data)) error stop 'f2k data: term_wpsclrp2_explicit'

        call c_f_pointer(input(26)%data, rp, [n])
        value%term_wpsclrp2_explicit = reshape(rp, shape(value%term_wpsclrp2_explicit))

      end if

    end if

    ! coef_wprtpsclrp_implicit

    if (input(27)%rank /= 3) error stop 'f2k rank: coef_wprtpsclrp_implicit'

    if (input(27)%allocated /= 0) then

      allocate(value%coef_wprtpsclrp_implicit(int(input(27)%extent(1)), int(input(27)%extent(2)), int(input(27)%extent(3))))
      n = int(product(input(27)%extent(1:3)))

      if (n > 0) then

        if (.not. c_associated(input(27)%data)) error stop 'f2k data: coef_wprtpsclrp_implicit'

        call c_f_pointer(input(27)%data, rp, [n])
        value%coef_wprtpsclrp_implicit = reshape(rp, shape(value%coef_wprtpsclrp_implicit))

      end if

    end if

    ! term_wprtpsclrp_explicit

    if (input(28)%rank /= 3) error stop 'f2k rank: term_wprtpsclrp_explicit'

    if (input(28)%allocated /= 0) then

      allocate(value%term_wprtpsclrp_explicit(int(input(28)%extent(1)), int(input(28)%extent(2)), int(input(28)%extent(3))))
      n = int(product(input(28)%extent(1:3)))

      if (n > 0) then

        if (.not. c_associated(input(28)%data)) error stop 'f2k data: term_wprtpsclrp_explicit'

        call c_f_pointer(input(28)%data, rp, [n])
        value%term_wprtpsclrp_explicit = reshape(rp, shape(value%term_wprtpsclrp_explicit))

      end if

    end if

    ! coef_wpthlpsclrp_implicit

    if (input(29)%rank /= 3) error stop 'f2k rank: coef_wpthlpsclrp_implicit'

    if (input(29)%allocated /= 0) then

      allocate(value%coef_wpthlpsclrp_implicit(int(input(29)%extent(1)), int(input(29)%extent(2)), int(input(29)%extent(3))))
      n = int(product(input(29)%extent(1:3)))

      if (n > 0) then

        if (.not. c_associated(input(29)%data)) error stop 'f2k data: coef_wpthlpsclrp_implicit'

        call c_f_pointer(input(29)%data, rp, [n])
        value%coef_wpthlpsclrp_implicit = reshape(rp, shape(value%coef_wpthlpsclrp_implicit))

      end if

    end if

    ! term_wpthlpsclrp_explicit

    if (input(30)%rank /= 3) error stop 'f2k rank: term_wpthlpsclrp_explicit'

    if (input(30)%allocated /= 0) then

      allocate(value%term_wpthlpsclrp_explicit(int(input(30)%extent(1)), int(input(30)%extent(2)), int(input(30)%extent(3))))
      n = int(product(input(30)%extent(1:3)))

      if (n > 0) then

        if (.not. c_associated(input(30)%data)) error stop 'f2k data: term_wpthlpsclrp_explicit'

        call c_f_pointer(input(30)%data, rp, [n])
        value%term_wpthlpsclrp_explicit = reshape(rp, shape(value%term_wpthlpsclrp_explicit))

      end if

    end if

  end subroutine to_native_implicit_coefs_terms

  subroutine from_native_implicit_coefs_terms(value, output)
    ! Description: Copy each implicit_coefs_terms field into C-owned output bytes.

    implicit none

    !--------------------------- Input Variables ---------------------------
    ! Call-local native argument; units and shape follow the source contract
    type(implicit_coefs_terms), intent(in) :: value

    !--------------------------- Output Variables ---------------------------
    type(wire_field), intent(out) :: output(*) ! Returned descriptors; C++ releases their C-owned buffers [-]

    !--------------------------- Local Variables ---------------------------
    integer(c_int), pointer :: ip(:) ! Borrowed view of packed integer data [-]
    ! Borrowed view of packed real data; units follow the native field
    real( kind = core_rknd ), pointer :: rp(:)
    integer :: n ! Number of packed elements [-]
    integer(c_int64_t) :: field_shape(4) ! Expected field extents in native dimension order [-]

    !--------------------------- Begin Code ---------------------------

    ! ngrdcol

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(1)) ! Out

    call c_f_pointer(output(1)%data, ip, [1])
    ip(1) = value%ngrdcol

    ! nz

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(2)) ! Out

    call c_f_pointer(output(2)%data, ip, [1])
    ip(1) = value%nz

    ! sclr_dim

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(3)) ! Out

    call c_f_pointer(output(3)%data, ip, [1])
    ip(1) = value%sclr_dim

    ! coef_wp4_implicit

    if (allocated(value%coef_wp4_implicit)) then

      n = size(value%coef_wp4_implicit)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%coef_wp4_implicit, 1, kind=c_int64_t)
      field_shape(2) = size(value%coef_wp4_implicit, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(4)) ! Out

      call c_f_pointer(output(4)%data, rp, [n])
      rp = reshape(value%coef_wp4_implicit, [n])

    else

      call empty_field(2, 0, & ! In
                       output(4)) ! Out

    end if

    ! coef_wp2rtp_implicit

    if (allocated(value%coef_wp2rtp_implicit)) then

      n = size(value%coef_wp2rtp_implicit)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%coef_wp2rtp_implicit, 1, kind=c_int64_t)
      field_shape(2) = size(value%coef_wp2rtp_implicit, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(5)) ! Out

      call c_f_pointer(output(5)%data, rp, [n])
      rp = reshape(value%coef_wp2rtp_implicit, [n])

    else

      call empty_field(2, 0, & ! In
                       output(5)) ! Out

    end if

    ! term_wp2rtp_explicit

    if (allocated(value%term_wp2rtp_explicit)) then

      n = size(value%term_wp2rtp_explicit)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%term_wp2rtp_explicit, 1, kind=c_int64_t)
      field_shape(2) = size(value%term_wp2rtp_explicit, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(6)) ! Out

      call c_f_pointer(output(6)%data, rp, [n])
      rp = reshape(value%term_wp2rtp_explicit, [n])

    else

      call empty_field(2, 0, & ! In
                       output(6)) ! Out

    end if

    ! coef_wp2thlp_implicit

    if (allocated(value%coef_wp2thlp_implicit)) then

      n = size(value%coef_wp2thlp_implicit)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%coef_wp2thlp_implicit, 1, kind=c_int64_t)
      field_shape(2) = size(value%coef_wp2thlp_implicit, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(7)) ! Out

      call c_f_pointer(output(7)%data, rp, [n])
      rp = reshape(value%coef_wp2thlp_implicit, [n])

    else

      call empty_field(2, 0, & ! In
                       output(7)) ! Out

    end if

    ! term_wp2thlp_explicit

    if (allocated(value%term_wp2thlp_explicit)) then

      n = size(value%term_wp2thlp_explicit)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%term_wp2thlp_explicit, 1, kind=c_int64_t)
      field_shape(2) = size(value%term_wp2thlp_explicit, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(8)) ! Out

      call c_f_pointer(output(8)%data, rp, [n])
      rp = reshape(value%term_wp2thlp_explicit, [n])

    else

      call empty_field(2, 0, & ! In
                       output(8)) ! Out

    end if

    ! coef_wp2up_implicit

    if (allocated(value%coef_wp2up_implicit)) then

      n = size(value%coef_wp2up_implicit)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%coef_wp2up_implicit, 1, kind=c_int64_t)
      field_shape(2) = size(value%coef_wp2up_implicit, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(9)) ! Out

      call c_f_pointer(output(9)%data, rp, [n])
      rp = reshape(value%coef_wp2up_implicit, [n])

    else

      call empty_field(2, 0, & ! In
                       output(9)) ! Out

    end if

    ! term_wp2up_explicit

    if (allocated(value%term_wp2up_explicit)) then

      n = size(value%term_wp2up_explicit)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%term_wp2up_explicit, 1, kind=c_int64_t)
      field_shape(2) = size(value%term_wp2up_explicit, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(10)) ! Out

      call c_f_pointer(output(10)%data, rp, [n])
      rp = reshape(value%term_wp2up_explicit, [n])

    else

      call empty_field(2, 0, & ! In
                       output(10)) ! Out

    end if

    ! coef_wp2vp_implicit

    if (allocated(value%coef_wp2vp_implicit)) then

      n = size(value%coef_wp2vp_implicit)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%coef_wp2vp_implicit, 1, kind=c_int64_t)
      field_shape(2) = size(value%coef_wp2vp_implicit, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(11)) ! Out

      call c_f_pointer(output(11)%data, rp, [n])
      rp = reshape(value%coef_wp2vp_implicit, [n])

    else

      call empty_field(2, 0, & ! In
                       output(11)) ! Out

    end if

    ! term_wp2vp_explicit

    if (allocated(value%term_wp2vp_explicit)) then

      n = size(value%term_wp2vp_explicit)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%term_wp2vp_explicit, 1, kind=c_int64_t)
      field_shape(2) = size(value%term_wp2vp_explicit, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(12)) ! Out

      call c_f_pointer(output(12)%data, rp, [n])
      rp = reshape(value%term_wp2vp_explicit, [n])

    else

      call empty_field(2, 0, & ! In
                       output(12)) ! Out

    end if

    ! coef_wprtp2_implicit

    if (allocated(value%coef_wprtp2_implicit)) then

      n = size(value%coef_wprtp2_implicit)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%coef_wprtp2_implicit, 1, kind=c_int64_t)
      field_shape(2) = size(value%coef_wprtp2_implicit, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(13)) ! Out

      call c_f_pointer(output(13)%data, rp, [n])
      rp = reshape(value%coef_wprtp2_implicit, [n])

    else

      call empty_field(2, 0, & ! In
                       output(13)) ! Out

    end if

    ! term_wprtp2_explicit

    if (allocated(value%term_wprtp2_explicit)) then

      n = size(value%term_wprtp2_explicit)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%term_wprtp2_explicit, 1, kind=c_int64_t)
      field_shape(2) = size(value%term_wprtp2_explicit, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(14)) ! Out

      call c_f_pointer(output(14)%data, rp, [n])
      rp = reshape(value%term_wprtp2_explicit, [n])

    else

      call empty_field(2, 0, & ! In
                       output(14)) ! Out

    end if

    ! coef_wpthlp2_implicit

    if (allocated(value%coef_wpthlp2_implicit)) then

      n = size(value%coef_wpthlp2_implicit)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%coef_wpthlp2_implicit, 1, kind=c_int64_t)
      field_shape(2) = size(value%coef_wpthlp2_implicit, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(15)) ! Out

      call c_f_pointer(output(15)%data, rp, [n])
      rp = reshape(value%coef_wpthlp2_implicit, [n])

    else

      call empty_field(2, 0, & ! In
                       output(15)) ! Out

    end if

    ! term_wpthlp2_explicit

    if (allocated(value%term_wpthlp2_explicit)) then

      n = size(value%term_wpthlp2_explicit)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%term_wpthlp2_explicit, 1, kind=c_int64_t)
      field_shape(2) = size(value%term_wpthlp2_explicit, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(16)) ! Out

      call c_f_pointer(output(16)%data, rp, [n])
      rp = reshape(value%term_wpthlp2_explicit, [n])

    else

      call empty_field(2, 0, & ! In
                       output(16)) ! Out

    end if

    ! coef_wprtpthlp_implicit

    if (allocated(value%coef_wprtpthlp_implicit)) then

      n = size(value%coef_wprtpthlp_implicit)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%coef_wprtpthlp_implicit, 1, kind=c_int64_t)
      field_shape(2) = size(value%coef_wprtpthlp_implicit, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(17)) ! Out

      call c_f_pointer(output(17)%data, rp, [n])
      rp = reshape(value%coef_wprtpthlp_implicit, [n])

    else

      call empty_field(2, 0, & ! In
                       output(17)) ! Out

    end if

    ! term_wprtpthlp_explicit

    if (allocated(value%term_wprtpthlp_explicit)) then

      n = size(value%term_wprtpthlp_explicit)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%term_wprtpthlp_explicit, 1, kind=c_int64_t)
      field_shape(2) = size(value%term_wprtpthlp_explicit, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(18)) ! Out

      call c_f_pointer(output(18)%data, rp, [n])
      rp = reshape(value%term_wprtpthlp_explicit, [n])

    else

      call empty_field(2, 0, & ! In
                       output(18)) ! Out

    end if

    ! coef_wpup2_implicit

    if (allocated(value%coef_wpup2_implicit)) then

      n = size(value%coef_wpup2_implicit)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%coef_wpup2_implicit, 1, kind=c_int64_t)
      field_shape(2) = size(value%coef_wpup2_implicit, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(19)) ! Out

      call c_f_pointer(output(19)%data, rp, [n])
      rp = reshape(value%coef_wpup2_implicit, [n])

    else

      call empty_field(2, 0, & ! In
                       output(19)) ! Out

    end if

    ! term_wpup2_explicit

    if (allocated(value%term_wpup2_explicit)) then

      n = size(value%term_wpup2_explicit)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%term_wpup2_explicit, 1, kind=c_int64_t)
      field_shape(2) = size(value%term_wpup2_explicit, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(20)) ! Out

      call c_f_pointer(output(20)%data, rp, [n])
      rp = reshape(value%term_wpup2_explicit, [n])

    else

      call empty_field(2, 0, & ! In
                       output(20)) ! Out

    end if

    ! coef_wpvp2_implicit

    if (allocated(value%coef_wpvp2_implicit)) then

      n = size(value%coef_wpvp2_implicit)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%coef_wpvp2_implicit, 1, kind=c_int64_t)
      field_shape(2) = size(value%coef_wpvp2_implicit, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(21)) ! Out

      call c_f_pointer(output(21)%data, rp, [n])
      rp = reshape(value%coef_wpvp2_implicit, [n])

    else

      call empty_field(2, 0, & ! In
                       output(21)) ! Out

    end if

    ! term_wpvp2_explicit

    if (allocated(value%term_wpvp2_explicit)) then

      n = size(value%term_wpvp2_explicit)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%term_wpvp2_explicit, 1, kind=c_int64_t)
      field_shape(2) = size(value%term_wpvp2_explicit, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(22)) ! Out

      call c_f_pointer(output(22)%data, rp, [n])
      rp = reshape(value%term_wpvp2_explicit, [n])

    else

      call empty_field(2, 0, & ! In
                       output(22)) ! Out

    end if

    ! coef_wp2sclrp_implicit

    if (allocated(value%coef_wp2sclrp_implicit)) then

      n = size(value%coef_wp2sclrp_implicit)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%coef_wp2sclrp_implicit, 1, kind=c_int64_t)
      field_shape(2) = size(value%coef_wp2sclrp_implicit, 2, kind=c_int64_t)
      field_shape(3) = size(value%coef_wp2sclrp_implicit, 3, kind=c_int64_t)

      call allocate_field(3, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(23)) ! Out

      call c_f_pointer(output(23)%data, rp, [n])
      rp = reshape(value%coef_wp2sclrp_implicit, [n])

    else

      call empty_field(3, 0, & ! In
                       output(23)) ! Out

    end if

    ! term_wp2sclrp_explicit

    if (allocated(value%term_wp2sclrp_explicit)) then

      n = size(value%term_wp2sclrp_explicit)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%term_wp2sclrp_explicit, 1, kind=c_int64_t)
      field_shape(2) = size(value%term_wp2sclrp_explicit, 2, kind=c_int64_t)
      field_shape(3) = size(value%term_wp2sclrp_explicit, 3, kind=c_int64_t)

      call allocate_field(3, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(24)) ! Out

      call c_f_pointer(output(24)%data, rp, [n])
      rp = reshape(value%term_wp2sclrp_explicit, [n])

    else

      call empty_field(3, 0, & ! In
                       output(24)) ! Out

    end if

    ! coef_wpsclrp2_implicit

    if (allocated(value%coef_wpsclrp2_implicit)) then

      n = size(value%coef_wpsclrp2_implicit)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%coef_wpsclrp2_implicit, 1, kind=c_int64_t)
      field_shape(2) = size(value%coef_wpsclrp2_implicit, 2, kind=c_int64_t)
      field_shape(3) = size(value%coef_wpsclrp2_implicit, 3, kind=c_int64_t)

      call allocate_field(3, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(25)) ! Out

      call c_f_pointer(output(25)%data, rp, [n])
      rp = reshape(value%coef_wpsclrp2_implicit, [n])

    else

      call empty_field(3, 0, & ! In
                       output(25)) ! Out

    end if

    ! term_wpsclrp2_explicit

    if (allocated(value%term_wpsclrp2_explicit)) then

      n = size(value%term_wpsclrp2_explicit)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%term_wpsclrp2_explicit, 1, kind=c_int64_t)
      field_shape(2) = size(value%term_wpsclrp2_explicit, 2, kind=c_int64_t)
      field_shape(3) = size(value%term_wpsclrp2_explicit, 3, kind=c_int64_t)

      call allocate_field(3, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(26)) ! Out

      call c_f_pointer(output(26)%data, rp, [n])
      rp = reshape(value%term_wpsclrp2_explicit, [n])

    else

      call empty_field(3, 0, & ! In
                       output(26)) ! Out

    end if

    ! coef_wprtpsclrp_implicit

    if (allocated(value%coef_wprtpsclrp_implicit)) then

      n = size(value%coef_wprtpsclrp_implicit)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%coef_wprtpsclrp_implicit, 1, kind=c_int64_t)
      field_shape(2) = size(value%coef_wprtpsclrp_implicit, 2, kind=c_int64_t)
      field_shape(3) = size(value%coef_wprtpsclrp_implicit, 3, kind=c_int64_t)

      call allocate_field(3, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(27)) ! Out

      call c_f_pointer(output(27)%data, rp, [n])
      rp = reshape(value%coef_wprtpsclrp_implicit, [n])

    else

      call empty_field(3, 0, & ! In
                       output(27)) ! Out

    end if

    ! term_wprtpsclrp_explicit

    if (allocated(value%term_wprtpsclrp_explicit)) then

      n = size(value%term_wprtpsclrp_explicit)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%term_wprtpsclrp_explicit, 1, kind=c_int64_t)
      field_shape(2) = size(value%term_wprtpsclrp_explicit, 2, kind=c_int64_t)
      field_shape(3) = size(value%term_wprtpsclrp_explicit, 3, kind=c_int64_t)

      call allocate_field(3, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(28)) ! Out

      call c_f_pointer(output(28)%data, rp, [n])
      rp = reshape(value%term_wprtpsclrp_explicit, [n])

    else

      call empty_field(3, 0, & ! In
                       output(28)) ! Out

    end if

    ! coef_wpthlpsclrp_implicit

    if (allocated(value%coef_wpthlpsclrp_implicit)) then

      n = size(value%coef_wpthlpsclrp_implicit)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%coef_wpthlpsclrp_implicit, 1, kind=c_int64_t)
      field_shape(2) = size(value%coef_wpthlpsclrp_implicit, 2, kind=c_int64_t)
      field_shape(3) = size(value%coef_wpthlpsclrp_implicit, 3, kind=c_int64_t)

      call allocate_field(3, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(29)) ! Out

      call c_f_pointer(output(29)%data, rp, [n])
      rp = reshape(value%coef_wpthlpsclrp_implicit, [n])

    else

      call empty_field(3, 0, & ! In
                       output(29)) ! Out

    end if

    ! term_wpthlpsclrp_explicit

    if (allocated(value%term_wpthlpsclrp_explicit)) then

      n = size(value%term_wpthlpsclrp_explicit)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%term_wpthlpsclrp_explicit, 1, kind=c_int64_t)
      field_shape(2) = size(value%term_wpthlpsclrp_explicit, 2, kind=c_int64_t)
      field_shape(3) = size(value%term_wpthlpsclrp_explicit, 3, kind=c_int64_t)

      call allocate_field(3, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(30)) ! Out

      call c_f_pointer(output(30)%data, rp, [n])
      rp = reshape(value%term_wpthlpsclrp_explicit, [n])

    else

      call empty_field(3, 0, & ! In
                       output(30)) ! Out

    end if

  end subroutine from_native_implicit_coefs_terms

  subroutine f2k_roundtrip_implicit_coefs_terms_c(input, output) bind(C)
    ! Description: Exercise both converters with one local Fortran value.

    implicit none

    !--------------------------- Input Variables ---------------------------
    type(wire_field), intent(in) :: input(*) ! Caller-owned input field descriptors [-]

    !--------------------------- Output Variables ---------------------------
    type(wire_field), intent(out) :: output(*) ! Returned descriptors; C++ releases their C-owned buffers [-]

    !--------------------------- Local Variables ---------------------------
    ! Call-local native argument; units and shape follow the source contract
    type(implicit_coefs_terms) :: value

    !--------------------------- Begin Code ---------------------------

    call to_native_implicit_coefs_terms(input, value)
    call from_native_implicit_coefs_terms(value, output)

  end subroutine f2k_roundtrip_implicit_coefs_terms_c

  subroutine to_native_err_info_type(input, value)
    ! Description: Copy each C field into a local err_info_type value.

    implicit none

    !--------------------------- Input Variables ---------------------------
    type(wire_field), intent(in) :: input(*) ! Caller-owned input field descriptors [-]

    !--------------------------- Output Variables ---------------------------
    ! Call-local native argument; units and shape follow the source contract
    type(err_info_type), intent(out) :: value

    !--------------------------- Local Variables ---------------------------
    integer(c_int), pointer :: ip(:) ! Borrowed view of packed integer data [-]
    ! Borrowed view of packed real data; units follow the native field
    real( kind = core_rknd ), pointer :: rp(:)
    character(c_char), pointer :: cp(:) ! Borrowed view of packed character bytes [-]
    integer :: n ! Number of packed elements [-]
    integer :: j ! Packed array or field index [-]
    integer :: k ! Character or secondary array index [-]
    ! Private error-context column latitudes [degrees north]
    real( kind = core_rknd ), allocatable :: private_host_lat(:)
    ! Private error-context column longitudes [degrees east]
    real( kind = core_rknd ), allocatable :: private_host_lon(:)
    integer :: private_chunk_idx ! Private error-context chunk index [-]
    integer :: private_mpi_rank ! Private error-context MPI rank [-]

    !--------------------------- Begin Code ---------------------------

    ! err_code

    if (input(1)%rank /= 1) error stop 'f2k rank: err_code'

    if (input(1)%allocated /= 0) then

      allocate(value%err_code(int(input(1)%extent(1))))
      n = int(product(input(1)%extent(1:1)))

      if (n > 0) then

        if (.not. c_associated(input(1)%data)) error stop 'f2k data: err_code'

        call c_f_pointer(input(1)%data, ip, [n])
        value%err_code = reshape(ip, shape(value%err_code))

      end if

    end if

    ! host_lat

    if (input(2)%rank /= 1) error stop 'f2k rank: host_lat'

    if (input(2)%allocated /= 0) then

      allocate(private_host_lat(int(input(2)%extent(1))))
      n = int(product(input(2)%extent(1:1)))

      if (n > 0) then

        if (.not. c_associated(input(2)%data)) error stop 'f2k data: host_lat'

        call c_f_pointer(input(2)%data, rp, [n])
        private_host_lat = reshape(rp, shape(private_host_lat))

      end if

    end if

    ! host_lon

    if (input(3)%rank /= 1) error stop 'f2k rank: host_lon'

    if (input(3)%allocated /= 0) then

      allocate(private_host_lon(int(input(3)%extent(1))))
      n = int(product(input(3)%extent(1:1)))

      if (n > 0) then

        if (.not. c_associated(input(3)%data)) error stop 'f2k data: host_lon'

        call c_f_pointer(input(3)%data, rp, [n])
        private_host_lon = reshape(rp, shape(private_host_lon))

      end if

    end if

    ! chunk_idx

    if (input(4)%rank /= 0) error stop 'f2k rank: chunk_idx'

    if (.not. c_associated(input(4)%data)) error stop 'f2k scalar: chunk_idx'

    call c_f_pointer(input(4)%data, ip, [1])
    private_chunk_idx = ip(1)

    ! mpi_rank

    if (input(5)%rank /= 0) error stop 'f2k rank: mpi_rank'

    if (.not. c_associated(input(5)%data)) error stop 'f2k scalar: mpi_rank'

    call c_f_pointer(input(5)%data, ip, [1])
    private_mpi_rank = ip(1)

    ! err_header

    if (input(6)%rank /= 1) error stop 'f2k rank: err_header'

    if (input(6)%allocated /= 0) then

      allocate(value%err_header(int(input(6)%extent(1))))

      if (input(6)%char_len /= 200) error stop 'f2k width: err_header'

      n = int(product(input(6)%extent(1:1)))

      if (n > 0) then

        if (.not. c_associated(input(6)%data)) error stop 'f2k data: err_header'

        call c_f_pointer(input(6)%data, cp, [n * 200])

        do j = 1, n

          do k = 1, 200

            value%err_header(j)(k:k) = cp((j-1)*200+k)

          end do

        end do

      end if

    end if

    ! err_header_global

    if (input(7)%rank /= 0) error stop 'f2k rank: err_header_global'

    if (.not. c_associated(input(7)%data)) error stop 'f2k scalar: err_header_global'

    if (input(7)%char_len /= 214) error stop 'f2k width: err_header_global'

    call c_f_pointer(input(7)%data, cp, [214])

    do k = 1, 214

      value%err_header_global(k:k) = cp(k)

    end do

    call f2k_set_err_info_private_api(private_host_lat, private_host_lon, & ! In
                                      private_chunk_idx, private_mpi_rank, & ! In
                                      value)                                  ! InOut

  end subroutine to_native_err_info_type

  subroutine from_native_err_info_type(value, output)
    ! Description: Copy each err_info_type field into C-owned output bytes.

    implicit none

    !--------------------------- Input Variables ---------------------------
    ! Call-local native argument; units and shape follow the source contract
    type(err_info_type), intent(in) :: value

    !--------------------------- Output Variables ---------------------------
    type(wire_field), intent(out) :: output(*) ! Returned descriptors; C++ releases their C-owned buffers [-]

    !--------------------------- Local Variables ---------------------------
    integer(c_int), pointer :: ip(:) ! Borrowed view of packed integer data [-]
    ! Borrowed view of packed real data; units follow the native field
    real( kind = core_rknd ), pointer :: rp(:)
    character(c_char), pointer :: cp(:) ! Borrowed view of packed character bytes [-]
    integer :: n ! Number of packed elements [-]
    integer :: j ! Packed array or field index [-]
    integer :: k ! Character or secondary array index [-]
    integer(c_int64_t) :: field_shape(4) ! Expected field extents in native dimension order [-]
    ! Private error-context column latitudes [degrees north]
    real( kind = core_rknd ), allocatable :: private_host_lat(:)
    ! Private error-context column longitudes [degrees east]
    real( kind = core_rknd ), allocatable :: private_host_lon(:)
    integer :: private_chunk_idx ! Private error-context chunk index [-]
    integer :: private_mpi_rank ! Private error-context MPI rank [-]

    !--------------------------- Begin Code ---------------------------

    call f2k_get_err_info_private_api(value, &                             ! In
                                      private_host_lat, private_host_lon, & ! Out
                                      private_chunk_idx, private_mpi_rank) ! Out

    ! err_code

    if (allocated(value%err_code)) then

      n = size(value%err_code)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%err_code, 1, kind=c_int64_t)

      call allocate_field(1, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0_c_int), 0, & ! In
                          output(1)) ! Out

      call c_f_pointer(output(1)%data, ip, [n])
      ip = reshape(value%err_code, [n])

    else

      call empty_field(1, 0, & ! In
                       output(1)) ! Out

    end if

    ! host_lat

    if (allocated(private_host_lat)) then

      n = size(private_host_lat)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(private_host_lat, 1, kind=c_int64_t)

      call allocate_field(1, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(2)) ! Out

      call c_f_pointer(output(2)%data, rp, [n])
      rp = reshape(private_host_lat, [n])

    else

      call empty_field(1, 0, & ! In
                       output(2)) ! Out

    end if

    ! host_lon

    if (allocated(private_host_lon)) then

      n = size(private_host_lon)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(private_host_lon, 1, kind=c_int64_t)

      call allocate_field(1, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(3)) ! Out

      call c_f_pointer(output(3)%data, rp, [n])
      rp = reshape(private_host_lon, [n])

    else

      call empty_field(1, 0, & ! In
                       output(3)) ! Out

    end if

    ! chunk_idx

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(4)) ! Out

    call c_f_pointer(output(4)%data, ip, [1])
    ip(1) = private_chunk_idx

    ! mpi_rank

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(5)) ! Out

    call c_f_pointer(output(5)%data, ip, [1])
    ip(1) = private_mpi_rank

    ! err_header

    if (allocated(value%err_header)) then

      n = size(value%err_header)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%err_header, 1, kind=c_int64_t)

      call allocate_field(1, field_shape, & ! In
                          int(n, c_size_t) * 200 * 1_c_size_t, 200, & ! In
                          output(6)) ! Out

      call c_f_pointer(output(6)%data, cp, [n * 200])

      do j = 1, n

        do k = 1, 200

          cp((j-1)*200+k) = value%err_header(j)(k:k)

        end do

      end do

    else

      call empty_field(1, 200, & ! In
                       output(6)) ! Out

    end if

    ! err_header_global

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(214, c_size_t) * 1_c_size_t, 214, & ! In
                        output(7)) ! Out

    call c_f_pointer(output(7)%data, cp, [214])

    do k = 1, 214

      cp(k) = value%err_header_global(k:k)

    end do

  end subroutine from_native_err_info_type

  subroutine f2k_roundtrip_err_info_type_c(input, output) bind(C)
    ! Description: Exercise both converters with one local Fortran value.

    implicit none

    !--------------------------- Input Variables ---------------------------
    type(wire_field), intent(in) :: input(*) ! Caller-owned input field descriptors [-]

    !--------------------------- Output Variables ---------------------------
    type(wire_field), intent(out) :: output(*) ! Returned descriptors; C++ releases their C-owned buffers [-]

    !--------------------------- Local Variables ---------------------------
    type(err_info_type) :: value ! Call-local native argument; units and shape follow the source contract

    !--------------------------- Begin Code ---------------------------

    call to_native_err_info_type(input, value)
    call from_native_err_info_type(value, output)

  end subroutine f2k_roundtrip_err_info_type_c

  subroutine to_native_hm_metadata_type(input, value)
    ! Description: Copy each C field into a local hm_metadata_type value.

    implicit none

    !--------------------------- Input Variables ---------------------------
    type(wire_field), intent(in) :: input(*) ! Caller-owned input field descriptors [-]

    !--------------------------- Output Variables ---------------------------
    ! Call-local native argument; units and shape follow the source contract
    type(hm_metadata_type), intent(out) :: value

    !--------------------------- Local Variables ---------------------------
    integer(c_int), pointer :: ip(:) ! Borrowed view of packed integer data [-]
    integer(c_int8_t), pointer :: bp(:) ! Borrowed view of packed logical bytes [-]
    ! Borrowed view of packed real data; units follow the native field
    real( kind = core_rknd ), pointer :: rp(:)
    character(c_char), pointer :: cp(:) ! Borrowed view of packed character bytes [-]
    integer :: n ! Number of packed elements [-]
    integer :: j ! Packed array or field index [-]
    integer :: k ! Character or secondary array index [-]

    !--------------------------- Begin Code ---------------------------

    ! iirr

    if (input(1)%rank /= 0) error stop 'f2k rank: iirr'

    if (.not. c_associated(input(1)%data)) error stop 'f2k scalar: iirr'

    call c_f_pointer(input(1)%data, ip, [1])
    value%iirr = ip(1)

    ! iirs

    if (input(2)%rank /= 0) error stop 'f2k rank: iirs'

    if (.not. c_associated(input(2)%data)) error stop 'f2k scalar: iirs'

    call c_f_pointer(input(2)%data, ip, [1])
    value%iirs = ip(1)

    ! iiri

    if (input(3)%rank /= 0) error stop 'f2k rank: iiri'

    if (.not. c_associated(input(3)%data)) error stop 'f2k scalar: iiri'

    call c_f_pointer(input(3)%data, ip, [1])
    value%iiri = ip(1)

    ! iirg

    if (input(4)%rank /= 0) error stop 'f2k rank: iirg'

    if (.not. c_associated(input(4)%data)) error stop 'f2k scalar: iirg'

    call c_f_pointer(input(4)%data, ip, [1])
    value%iirg = ip(1)

    ! iiNr

    if (input(5)%rank /= 0) error stop 'f2k rank: iiNr'

    if (.not. c_associated(input(5)%data)) error stop 'f2k scalar: iiNr'

    call c_f_pointer(input(5)%data, ip, [1])
    value%iiNr = ip(1)

    ! iiNs

    if (input(6)%rank /= 0) error stop 'f2k rank: iiNs'

    if (.not. c_associated(input(6)%data)) error stop 'f2k scalar: iiNs'

    call c_f_pointer(input(6)%data, ip, [1])
    value%iiNs = ip(1)

    ! iiNi

    if (input(7)%rank /= 0) error stop 'f2k rank: iiNi'

    if (.not. c_associated(input(7)%data)) error stop 'f2k scalar: iiNi'

    call c_f_pointer(input(7)%data, ip, [1])
    value%iiNi = ip(1)

    ! iiNg

    if (input(8)%rank /= 0) error stop 'f2k rank: iiNg'

    if (.not. c_associated(input(8)%data)) error stop 'f2k scalar: iiNg'

    call c_f_pointer(input(8)%data, ip, [1])
    value%iiNg = ip(1)

    ! l_frozen_hm

    if (input(9)%rank /= 1) error stop 'f2k rank: l_frozen_hm'

    if (input(9)%allocated /= 0) then

      allocate(value%l_frozen_hm(int(input(9)%extent(1))))
      n = int(product(input(9)%extent(1:1)))

      if (n > 0) then

        if (.not. c_associated(input(9)%data)) error stop 'f2k data: l_frozen_hm'

        call c_f_pointer(input(9)%data, bp, [n])
        value%l_frozen_hm = reshape(bp /= 0_c_int8_t, shape(value%l_frozen_hm))

      end if

    end if

    ! l_mix_rat_hm

    if (input(10)%rank /= 1) error stop 'f2k rank: l_mix_rat_hm'

    if (input(10)%allocated /= 0) then

      allocate(value%l_mix_rat_hm(int(input(10)%extent(1))))
      n = int(product(input(10)%extent(1:1)))

      if (n > 0) then

        if (.not. c_associated(input(10)%data)) error stop 'f2k data: l_mix_rat_hm'

        call c_f_pointer(input(10)%data, bp, [n])
        value%l_mix_rat_hm = reshape(bp /= 0_c_int8_t, shape(value%l_mix_rat_hm))

      end if

    end if

    ! hydromet_list

    if (input(11)%rank /= 1) error stop 'f2k rank: hydromet_list'

    if (input(11)%allocated /= 0) then

      allocate(value%hydromet_list(int(input(11)%extent(1))))

      if (input(11)%char_len /= 10) error stop 'f2k width: hydromet_list'

      n = int(product(input(11)%extent(1:1)))

      if (n > 0) then

        if (.not. c_associated(input(11)%data)) error stop 'f2k data: hydromet_list'

        call c_f_pointer(input(11)%data, cp, [n * 10])

        do j = 1, n

          do k = 1, 10

            value%hydromet_list(j)(k:k) = cp((j-1)*10+k)

          end do

        end do

      end if

    end if

    ! hydromet_tol

    if (input(12)%rank /= 1) error stop 'f2k rank: hydromet_tol'

    if (input(12)%allocated /= 0) then

      allocate(value%hydromet_tol(int(input(12)%extent(1))))
      n = int(product(input(12)%extent(1:1)))

      if (n > 0) then

        if (.not. c_associated(input(12)%data)) error stop 'f2k data: hydromet_tol'

        call c_f_pointer(input(12)%data, rp, [n])
        value%hydromet_tol = reshape(rp, shape(value%hydromet_tol))

      end if

    end if

    ! iiPDF_chi

    if (input(13)%rank /= 0) error stop 'f2k rank: iiPDF_chi'

    if (.not. c_associated(input(13)%data)) error stop 'f2k scalar: iiPDF_chi'

    call c_f_pointer(input(13)%data, ip, [1])
    value%iiPDF_chi = ip(1)

    ! iiPDF_eta

    if (input(14)%rank /= 0) error stop 'f2k rank: iiPDF_eta'

    if (.not. c_associated(input(14)%data)) error stop 'f2k scalar: iiPDF_eta'

    call c_f_pointer(input(14)%data, ip, [1])
    value%iiPDF_eta = ip(1)

    ! iiPDF_w

    if (input(15)%rank /= 0) error stop 'f2k rank: iiPDF_w'

    if (.not. c_associated(input(15)%data)) error stop 'f2k scalar: iiPDF_w'

    call c_f_pointer(input(15)%data, ip, [1])
    value%iiPDF_w = ip(1)

    ! iiPDF_rr

    if (input(16)%rank /= 0) error stop 'f2k rank: iiPDF_rr'

    if (.not. c_associated(input(16)%data)) error stop 'f2k scalar: iiPDF_rr'

    call c_f_pointer(input(16)%data, ip, [1])
    value%iiPDF_rr = ip(1)

    ! iiPDF_rs

    if (input(17)%rank /= 0) error stop 'f2k rank: iiPDF_rs'

    if (.not. c_associated(input(17)%data)) error stop 'f2k scalar: iiPDF_rs'

    call c_f_pointer(input(17)%data, ip, [1])
    value%iiPDF_rs = ip(1)

    ! iiPDF_ri

    if (input(18)%rank /= 0) error stop 'f2k rank: iiPDF_ri'

    if (.not. c_associated(input(18)%data)) error stop 'f2k scalar: iiPDF_ri'

    call c_f_pointer(input(18)%data, ip, [1])
    value%iiPDF_ri = ip(1)

    ! iiPDF_rg

    if (input(19)%rank /= 0) error stop 'f2k rank: iiPDF_rg'

    if (.not. c_associated(input(19)%data)) error stop 'f2k scalar: iiPDF_rg'

    call c_f_pointer(input(19)%data, ip, [1])
    value%iiPDF_rg = ip(1)

    ! iiPDF_Nr

    if (input(20)%rank /= 0) error stop 'f2k rank: iiPDF_Nr'

    if (.not. c_associated(input(20)%data)) error stop 'f2k scalar: iiPDF_Nr'

    call c_f_pointer(input(20)%data, ip, [1])
    value%iiPDF_Nr = ip(1)

    ! iiPDF_Ns

    if (input(21)%rank /= 0) error stop 'f2k rank: iiPDF_Ns'

    if (.not. c_associated(input(21)%data)) error stop 'f2k scalar: iiPDF_Ns'

    call c_f_pointer(input(21)%data, ip, [1])
    value%iiPDF_Ns = ip(1)

    ! iiPDF_Ni

    if (input(22)%rank /= 0) error stop 'f2k rank: iiPDF_Ni'

    if (.not. c_associated(input(22)%data)) error stop 'f2k scalar: iiPDF_Ni'

    call c_f_pointer(input(22)%data, ip, [1])
    value%iiPDF_Ni = ip(1)

    ! iiPDF_Ng

    if (input(23)%rank /= 0) error stop 'f2k rank: iiPDF_Ng'

    if (.not. c_associated(input(23)%data)) error stop 'f2k scalar: iiPDF_Ng'

    call c_f_pointer(input(23)%data, ip, [1])
    value%iiPDF_Ng = ip(1)

    ! iiPDF_Ncn

    if (input(24)%rank /= 0) error stop 'f2k rank: iiPDF_Ncn'

    if (.not. c_associated(input(24)%data)) error stop 'f2k scalar: iiPDF_Ncn'

    call c_f_pointer(input(24)%data, ip, [1])
    value%iiPDF_Ncn = ip(1)

    ! hmp2_ip_on_hmm2_ip

    if (input(25)%rank /= 1) error stop 'f2k rank: hmp2_ip_on_hmm2_ip'

    if (input(25)%allocated /= 0) then

      allocate(value%hmp2_ip_on_hmm2_ip(int(input(25)%extent(1))))
      n = int(product(input(25)%extent(1:1)))

      if (n > 0) then

        if (.not. c_associated(input(25)%data)) error stop 'f2k data: hmp2_ip_on_hmm2_ip'

        call c_f_pointer(input(25)%data, rp, [n])
        value%hmp2_ip_on_hmm2_ip = reshape(rp, shape(value%hmp2_ip_on_hmm2_ip))

      end if

    end if

    ! Ncnp2_on_Ncnm2

    if (input(26)%rank /= 0) error stop 'f2k rank: Ncnp2_on_Ncnm2'

    if (.not. c_associated(input(26)%data)) error stop 'f2k scalar: Ncnp2_on_Ncnm2'

    call c_f_pointer(input(26)%data, rp, [1])
    value%Ncnp2_on_Ncnm2 = rp(1)

  end subroutine to_native_hm_metadata_type

  subroutine from_native_hm_metadata_type(value, output)
    ! Description: Copy each hm_metadata_type field into C-owned output bytes.

    implicit none

    !--------------------------- Input Variables ---------------------------
    ! Call-local native argument; units and shape follow the source contract
    type(hm_metadata_type), intent(in) :: value

    !--------------------------- Output Variables ---------------------------
    type(wire_field), intent(out) :: output(*) ! Returned descriptors; C++ releases their C-owned buffers [-]

    !--------------------------- Local Variables ---------------------------
    integer(c_int), pointer :: ip(:) ! Borrowed view of packed integer data [-]
    integer(c_int8_t), pointer :: bp(:) ! Borrowed view of packed logical bytes [-]
    ! Borrowed view of packed real data; units follow the native field
    real( kind = core_rknd ), pointer :: rp(:)
    character(c_char), pointer :: cp(:) ! Borrowed view of packed character bytes [-]
    integer :: n ! Number of packed elements [-]
    integer :: j ! Packed array or field index [-]
    integer :: k ! Character or secondary array index [-]
    integer(c_int64_t) :: field_shape(4) ! Expected field extents in native dimension order [-]

    !--------------------------- Begin Code ---------------------------

    ! iirr

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(1)) ! Out

    call c_f_pointer(output(1)%data, ip, [1])
    ip(1) = value%iirr

    ! iirs

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(2)) ! Out

    call c_f_pointer(output(2)%data, ip, [1])
    ip(1) = value%iirs

    ! iiri

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(3)) ! Out

    call c_f_pointer(output(3)%data, ip, [1])
    ip(1) = value%iiri

    ! iirg

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(4)) ! Out

    call c_f_pointer(output(4)%data, ip, [1])
    ip(1) = value%iirg

    ! iiNr

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(5)) ! Out

    call c_f_pointer(output(5)%data, ip, [1])
    ip(1) = value%iiNr

    ! iiNs

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(6)) ! Out

    call c_f_pointer(output(6)%data, ip, [1])
    ip(1) = value%iiNs

    ! iiNi

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(7)) ! Out

    call c_f_pointer(output(7)%data, ip, [1])
    ip(1) = value%iiNi

    ! iiNg

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(8)) ! Out

    call c_f_pointer(output(8)%data, ip, [1])
    ip(1) = value%iiNg

    ! l_frozen_hm

    if (allocated(value%l_frozen_hm)) then

      n = size(value%l_frozen_hm)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%l_frozen_hm, 1, kind=c_int64_t)

      call allocate_field(1, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0_c_int8_t), 0, & ! In
                          output(9)) ! Out

      call c_f_pointer(output(9)%data, bp, [n])
      bp = merge(1_c_int8_t, 0_c_int8_t, reshape(value%l_frozen_hm, [n]))

    else

      call empty_field(1, 0, & ! In
                       output(9)) ! Out

    end if

    ! l_mix_rat_hm

    if (allocated(value%l_mix_rat_hm)) then

      n = size(value%l_mix_rat_hm)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%l_mix_rat_hm, 1, kind=c_int64_t)

      call allocate_field(1, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0_c_int8_t), 0, & ! In
                          output(10)) ! Out

      call c_f_pointer(output(10)%data, bp, [n])
      bp = merge(1_c_int8_t, 0_c_int8_t, reshape(value%l_mix_rat_hm, [n]))

    else

      call empty_field(1, 0, & ! In
                       output(10)) ! Out

    end if

    ! hydromet_list

    if (allocated(value%hydromet_list)) then

      n = size(value%hydromet_list)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%hydromet_list, 1, kind=c_int64_t)

      call allocate_field(1, field_shape, & ! In
                          int(n, c_size_t) * 10 * 1_c_size_t, 10, & ! In
                          output(11)) ! Out

      call c_f_pointer(output(11)%data, cp, [n * 10])

      do j = 1, n

        do k = 1, 10

          cp((j-1)*10+k) = value%hydromet_list(j)(k:k)

        end do

      end do

    else

      call empty_field(1, 10, & ! In
                       output(11)) ! Out

    end if

    ! hydromet_tol

    if (allocated(value%hydromet_tol)) then

      n = size(value%hydromet_tol)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%hydromet_tol, 1, kind=c_int64_t)

      call allocate_field(1, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(12)) ! Out

      call c_f_pointer(output(12)%data, rp, [n])
      rp = reshape(value%hydromet_tol, [n])

    else

      call empty_field(1, 0, & ! In
                       output(12)) ! Out

    end if

    ! iiPDF_chi

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(13)) ! Out

    call c_f_pointer(output(13)%data, ip, [1])
    ip(1) = value%iiPDF_chi

    ! iiPDF_eta

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(14)) ! Out

    call c_f_pointer(output(14)%data, ip, [1])
    ip(1) = value%iiPDF_eta

    ! iiPDF_w

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(15)) ! Out

    call c_f_pointer(output(15)%data, ip, [1])
    ip(1) = value%iiPDF_w

    ! iiPDF_rr

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(16)) ! Out

    call c_f_pointer(output(16)%data, ip, [1])
    ip(1) = value%iiPDF_rr

    ! iiPDF_rs

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(17)) ! Out

    call c_f_pointer(output(17)%data, ip, [1])
    ip(1) = value%iiPDF_rs

    ! iiPDF_ri

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(18)) ! Out

    call c_f_pointer(output(18)%data, ip, [1])
    ip(1) = value%iiPDF_ri

    ! iiPDF_rg

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(19)) ! Out

    call c_f_pointer(output(19)%data, ip, [1])
    ip(1) = value%iiPDF_rg

    ! iiPDF_Nr

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(20)) ! Out

    call c_f_pointer(output(20)%data, ip, [1])
    ip(1) = value%iiPDF_Nr

    ! iiPDF_Ns

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(21)) ! Out

    call c_f_pointer(output(21)%data, ip, [1])
    ip(1) = value%iiPDF_Ns

    ! iiPDF_Ni

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(22)) ! Out

    call c_f_pointer(output(22)%data, ip, [1])
    ip(1) = value%iiPDF_Ni

    ! iiPDF_Ng

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(23)) ! Out

    call c_f_pointer(output(23)%data, ip, [1])
    ip(1) = value%iiPDF_Ng

    ! iiPDF_Ncn

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(24)) ! Out

    call c_f_pointer(output(24)%data, ip, [1])
    ip(1) = value%iiPDF_Ncn

    ! hmp2_ip_on_hmm2_ip

    if (allocated(value%hmp2_ip_on_hmm2_ip)) then

      n = size(value%hmp2_ip_on_hmm2_ip)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%hmp2_ip_on_hmm2_ip, 1, kind=c_int64_t)

      call allocate_field(1, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(25)) ! Out

      call c_f_pointer(output(25)%data, rp, [n])
      rp = reshape(value%hmp2_ip_on_hmm2_ip, [n])

    else

      call empty_field(1, 0, & ! In
                       output(25)) ! Out

    end if

    ! Ncnp2_on_Ncnm2

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0.0_core_rknd), 0, & ! In
                        output(26)) ! Out

    call c_f_pointer(output(26)%data, rp, [1])
    rp(1) = value%Ncnp2_on_Ncnm2

  end subroutine from_native_hm_metadata_type

  subroutine f2k_roundtrip_hm_metadata_type_c(input, output) bind(C)
    ! Description: Exercise both converters with one local Fortran value.

    implicit none

    !--------------------------- Input Variables ---------------------------
    type(wire_field), intent(in) :: input(*) ! Caller-owned input field descriptors [-]

    !--------------------------- Output Variables ---------------------------
    type(wire_field), intent(out) :: output(*) ! Returned descriptors; C++ releases their C-owned buffers [-]

    !--------------------------- Local Variables ---------------------------
    type(hm_metadata_type) :: value ! Call-local native argument; units and shape follow the source contract

    !--------------------------- Begin Code ---------------------------

    call to_native_hm_metadata_type(input, value)
    call from_native_hm_metadata_type(value, output)

  end subroutine f2k_roundtrip_hm_metadata_type_c

  subroutine to_native_silhs_config_flags_type(input, value)
    ! Description: Copy each C field into a local silhs_config_flags_type value.

    implicit none

    !--------------------------- Input Variables ---------------------------
    type(wire_field), intent(in) :: input(*) ! Caller-owned input field descriptors [-]

    !--------------------------- Output Variables ---------------------------
    ! Call-local native argument; units and shape follow the source contract
    type(silhs_config_flags_type), intent(out) :: value

    !--------------------------- Local Variables ---------------------------
    integer(c_int), pointer :: ip(:) ! Borrowed view of packed integer data [-]
    integer(c_int8_t), pointer :: bp(:) ! Borrowed view of packed logical bytes [-]

    !--------------------------- Begin Code ---------------------------

    ! cluster_allocation_strategy

    if (input(1)%rank /= 0) error stop 'f2k rank: cluster_allocation_strategy'

    if (.not. c_associated(input(1)%data)) error stop 'f2k scalar: cluster_allocation_strategy'

    call c_f_pointer(input(1)%data, ip, [1])
    value%cluster_allocation_strategy = ip(1)

    ! l_lh_importance_sampling

    if (input(2)%rank /= 0) error stop 'f2k rank: l_lh_importance_sampling'

    if (.not. c_associated(input(2)%data)) error stop 'f2k scalar: l_lh_importance_sampling'

    call c_f_pointer(input(2)%data, bp, [1])
    value%l_lh_importance_sampling = bp(1) /= 0_c_int8_t

    ! l_Lscale_vert_avg

    if (input(3)%rank /= 0) error stop 'f2k rank: l_Lscale_vert_avg'

    if (.not. c_associated(input(3)%data)) error stop 'f2k scalar: l_Lscale_vert_avg'

    call c_f_pointer(input(3)%data, bp, [1])
    value%l_Lscale_vert_avg = bp(1) /= 0_c_int8_t

    ! l_lh_straight_mc

    if (input(4)%rank /= 0) error stop 'f2k rank: l_lh_straight_mc'

    if (.not. c_associated(input(4)%data)) error stop 'f2k scalar: l_lh_straight_mc'

    call c_f_pointer(input(4)%data, bp, [1])
    value%l_lh_straight_mc = bp(1) /= 0_c_int8_t

    ! l_lh_clustered_sampling

    if (input(5)%rank /= 0) error stop 'f2k rank: l_lh_clustered_sampling'

    if (.not. c_associated(input(5)%data)) error stop 'f2k scalar: l_lh_clustered_sampling'

    call c_f_pointer(input(5)%data, bp, [1])
    value%l_lh_clustered_sampling = bp(1) /= 0_c_int8_t

    ! l_rcm_in_cloud_k_lh_start

    if (input(6)%rank /= 0) error stop 'f2k rank: l_rcm_in_cloud_k_lh_start'

    if (.not. c_associated(input(6)%data)) error stop 'f2k scalar: l_rcm_in_cloud_k_lh_start'

    call c_f_pointer(input(6)%data, bp, [1])
    value%l_rcm_in_cloud_k_lh_start = bp(1) /= 0_c_int8_t

    ! l_random_k_lh_start

    if (input(7)%rank /= 0) error stop 'f2k rank: l_random_k_lh_start'

    if (.not. c_associated(input(7)%data)) error stop 'f2k scalar: l_random_k_lh_start'

    call c_f_pointer(input(7)%data, bp, [1])
    value%l_random_k_lh_start = bp(1) /= 0_c_int8_t

    ! l_max_overlap_in_cloud

    if (input(8)%rank /= 0) error stop 'f2k rank: l_max_overlap_in_cloud'

    if (.not. c_associated(input(8)%data)) error stop 'f2k scalar: l_max_overlap_in_cloud'

    call c_f_pointer(input(8)%data, bp, [1])
    value%l_max_overlap_in_cloud = bp(1) /= 0_c_int8_t

    ! l_lh_instant_var_covar_src

    if (input(9)%rank /= 0) error stop 'f2k rank: l_lh_instant_var_covar_src'

    if (.not. c_associated(input(9)%data)) error stop 'f2k scalar: l_lh_instant_var_covar_src'

    call c_f_pointer(input(9)%data, bp, [1])
    value%l_lh_instant_var_covar_src = bp(1) /= 0_c_int8_t

    ! l_lh_limit_weights

    if (input(10)%rank /= 0) error stop 'f2k rank: l_lh_limit_weights'

    if (.not. c_associated(input(10)%data)) error stop 'f2k scalar: l_lh_limit_weights'

    call c_f_pointer(input(10)%data, bp, [1])
    value%l_lh_limit_weights = bp(1) /= 0_c_int8_t

    ! l_lh_var_frac

    if (input(11)%rank /= 0) error stop 'f2k rank: l_lh_var_frac'

    if (.not. c_associated(input(11)%data)) error stop 'f2k scalar: l_lh_var_frac'

    call c_f_pointer(input(11)%data, bp, [1])
    value%l_lh_var_frac = bp(1) /= 0_c_int8_t

    ! l_lh_normalize_weights

    if (input(12)%rank /= 0) error stop 'f2k rank: l_lh_normalize_weights'

    if (.not. c_associated(input(12)%data)) error stop 'f2k scalar: l_lh_normalize_weights'

    call c_f_pointer(input(12)%data, bp, [1])
    value%l_lh_normalize_weights = bp(1) /= 0_c_int8_t

    ! l_lh_deterministic_test

    if (input(13)%rank /= 0) error stop 'f2k rank: l_lh_deterministic_test'

    if (.not. c_associated(input(13)%data)) error stop 'f2k scalar: l_lh_deterministic_test'

    call c_f_pointer(input(13)%data, bp, [1])
    value%l_lh_deterministic_test = bp(1) /= 0_c_int8_t

  end subroutine to_native_silhs_config_flags_type

  subroutine from_native_silhs_config_flags_type(value, output)
    ! Description: Copy each silhs_config_flags_type field into C-owned output bytes.

    implicit none

    !--------------------------- Input Variables ---------------------------
    ! Call-local native argument; units and shape follow the source contract
    type(silhs_config_flags_type), intent(in) :: value

    !--------------------------- Output Variables ---------------------------
    type(wire_field), intent(out) :: output(*) ! Returned descriptors; C++ releases their C-owned buffers [-]

    !--------------------------- Local Variables ---------------------------
    integer(c_int), pointer :: ip(:) ! Borrowed view of packed integer data [-]
    integer(c_int8_t), pointer :: bp(:) ! Borrowed view of packed logical bytes [-]

    !--------------------------- Begin Code ---------------------------

    ! cluster_allocation_strategy

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(1)) ! Out

    call c_f_pointer(output(1)%data, ip, [1])
    ip(1) = value%cluster_allocation_strategy

    ! l_lh_importance_sampling

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(2)) ! Out

    call c_f_pointer(output(2)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_lh_importance_sampling)

    ! l_Lscale_vert_avg

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(3)) ! Out

    call c_f_pointer(output(3)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_Lscale_vert_avg)

    ! l_lh_straight_mc

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(4)) ! Out

    call c_f_pointer(output(4)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_lh_straight_mc)

    ! l_lh_clustered_sampling

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(5)) ! Out

    call c_f_pointer(output(5)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_lh_clustered_sampling)

    ! l_rcm_in_cloud_k_lh_start

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(6)) ! Out

    call c_f_pointer(output(6)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_rcm_in_cloud_k_lh_start)

    ! l_random_k_lh_start

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(7)) ! Out

    call c_f_pointer(output(7)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_random_k_lh_start)

    ! l_max_overlap_in_cloud

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(8)) ! Out

    call c_f_pointer(output(8)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_max_overlap_in_cloud)

    ! l_lh_instant_var_covar_src

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(9)) ! Out

    call c_f_pointer(output(9)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_lh_instant_var_covar_src)

    ! l_lh_limit_weights

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(10)) ! Out

    call c_f_pointer(output(10)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_lh_limit_weights)

    ! l_lh_var_frac

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(11)) ! Out

    call c_f_pointer(output(11)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_lh_var_frac)

    ! l_lh_normalize_weights

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(12)) ! Out

    call c_f_pointer(output(12)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_lh_normalize_weights)

    ! l_lh_deterministic_test

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int8_t), 0, & ! In
                        output(13)) ! Out

    call c_f_pointer(output(13)%data, bp, [1])
    bp(1) = merge(1_c_int8_t, 0_c_int8_t, value%l_lh_deterministic_test)

  end subroutine from_native_silhs_config_flags_type

  subroutine f2k_roundtrip_silhs_config_flags_type_c(input, output) bind(C)
    ! Description: Exercise both converters with one local Fortran value.

    implicit none

    !--------------------------- Input Variables ---------------------------
    type(wire_field), intent(in) :: input(*) ! Caller-owned input field descriptors [-]

    !--------------------------- Output Variables ---------------------------
    type(wire_field), intent(out) :: output(*) ! Returned descriptors; C++ releases their C-owned buffers [-]

    !--------------------------- Local Variables ---------------------------
    ! Call-local native argument; units and shape follow the source contract
    type(silhs_config_flags_type) :: value

    !--------------------------- Begin Code ---------------------------

    call to_native_silhs_config_flags_type(input, value)
    call from_native_silhs_config_flags_type(value, output)

  end subroutine f2k_roundtrip_silhs_config_flags_type_c

  subroutine to_native_hydromet_pdf_parameter(input, value)
    ! Description: Copy each C field into a local hydromet_pdf_parameter value.

    implicit none

    !--------------------------- Input Variables ---------------------------
    type(wire_field), intent(in) :: input(*) ! Caller-owned input field descriptors [-]

    !--------------------------- Output Variables ---------------------------
    ! Call-local native argument; units and shape follow the source contract
    type(hydromet_pdf_parameter), intent(out) :: value

    !--------------------------- Local Variables ---------------------------
    ! Borrowed view of packed real data; units follow the native field
    real( kind = core_rknd ), pointer :: rp(:)
    integer :: n ! Number of packed elements [-]
    integer(c_int64_t) :: field_shape(4) ! Expected field extents in native dimension order [-]

    !--------------------------- Begin Code ---------------------------

    ! hm_1

    if (input(1)%rank /= 1) error stop 'f2k rank: hm_1'

    if (input(1)%allocated /= 0) then

      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%hm_1, 1, kind=c_int64_t)

      if (any(input(1)%extent /= field_shape)) &
        error stop 'f2k shape: hm_1'

      n = int(product(input(1)%extent(1:1)))

      if (n > 0) then

        if (.not. c_associated(input(1)%data)) error stop 'f2k data: hm_1'

        call c_f_pointer(input(1)%data, rp, [n])
        value%hm_1 = reshape(rp, shape(value%hm_1))

      end if

    else

      error stop 'f2k missing fixed array: hm_1'

    end if

    ! hm_2

    if (input(2)%rank /= 1) error stop 'f2k rank: hm_2'

    if (input(2)%allocated /= 0) then

      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%hm_2, 1, kind=c_int64_t)

      if (any(input(2)%extent /= field_shape)) &
        error stop 'f2k shape: hm_2'

      n = int(product(input(2)%extent(1:1)))

      if (n > 0) then

        if (.not. c_associated(input(2)%data)) error stop 'f2k data: hm_2'

        call c_f_pointer(input(2)%data, rp, [n])
        value%hm_2 = reshape(rp, shape(value%hm_2))

      end if

    else

      error stop 'f2k missing fixed array: hm_2'

    end if

    ! mu_hm_1

    if (input(3)%rank /= 1) error stop 'f2k rank: mu_hm_1'

    if (input(3)%allocated /= 0) then

      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%mu_hm_1, 1, kind=c_int64_t)

      if (any(input(3)%extent /= field_shape)) &
        error stop 'f2k shape: mu_hm_1'

      n = int(product(input(3)%extent(1:1)))

      if (n > 0) then

        if (.not. c_associated(input(3)%data)) error stop 'f2k data: mu_hm_1'

        call c_f_pointer(input(3)%data, rp, [n])
        value%mu_hm_1 = reshape(rp, shape(value%mu_hm_1))

      end if

    else

      error stop 'f2k missing fixed array: mu_hm_1'

    end if

    ! mu_hm_2

    if (input(4)%rank /= 1) error stop 'f2k rank: mu_hm_2'

    if (input(4)%allocated /= 0) then

      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%mu_hm_2, 1, kind=c_int64_t)

      if (any(input(4)%extent /= field_shape)) &
        error stop 'f2k shape: mu_hm_2'

      n = int(product(input(4)%extent(1:1)))

      if (n > 0) then

        if (.not. c_associated(input(4)%data)) error stop 'f2k data: mu_hm_2'

        call c_f_pointer(input(4)%data, rp, [n])
        value%mu_hm_2 = reshape(rp, shape(value%mu_hm_2))

      end if

    else

      error stop 'f2k missing fixed array: mu_hm_2'

    end if

    ! sigma_hm_1

    if (input(5)%rank /= 1) error stop 'f2k rank: sigma_hm_1'

    if (input(5)%allocated /= 0) then

      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%sigma_hm_1, 1, kind=c_int64_t)

      if (any(input(5)%extent /= field_shape)) &
        error stop 'f2k shape: sigma_hm_1'

      n = int(product(input(5)%extent(1:1)))

      if (n > 0) then

        if (.not. c_associated(input(5)%data)) error stop 'f2k data: sigma_hm_1'

        call c_f_pointer(input(5)%data, rp, [n])
        value%sigma_hm_1 = reshape(rp, shape(value%sigma_hm_1))

      end if

    else

      error stop 'f2k missing fixed array: sigma_hm_1'

    end if

    ! sigma_hm_2

    if (input(6)%rank /= 1) error stop 'f2k rank: sigma_hm_2'

    if (input(6)%allocated /= 0) then

      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%sigma_hm_2, 1, kind=c_int64_t)

      if (any(input(6)%extent /= field_shape)) &
        error stop 'f2k shape: sigma_hm_2'

      n = int(product(input(6)%extent(1:1)))

      if (n > 0) then

        if (.not. c_associated(input(6)%data)) error stop 'f2k data: sigma_hm_2'

        call c_f_pointer(input(6)%data, rp, [n])
        value%sigma_hm_2 = reshape(rp, shape(value%sigma_hm_2))

      end if

    else

      error stop 'f2k missing fixed array: sigma_hm_2'

    end if

    ! corr_w_hm_1

    if (input(7)%rank /= 1) error stop 'f2k rank: corr_w_hm_1'

    if (input(7)%allocated /= 0) then

      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%corr_w_hm_1, 1, kind=c_int64_t)

      if (any(input(7)%extent /= field_shape)) &
        error stop 'f2k shape: corr_w_hm_1'

      n = int(product(input(7)%extent(1:1)))

      if (n > 0) then

        if (.not. c_associated(input(7)%data)) error stop 'f2k data: corr_w_hm_1'

        call c_f_pointer(input(7)%data, rp, [n])
        value%corr_w_hm_1 = reshape(rp, shape(value%corr_w_hm_1))

      end if

    else

      error stop 'f2k missing fixed array: corr_w_hm_1'

    end if

    ! corr_w_hm_2

    if (input(8)%rank /= 1) error stop 'f2k rank: corr_w_hm_2'

    if (input(8)%allocated /= 0) then

      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%corr_w_hm_2, 1, kind=c_int64_t)

      if (any(input(8)%extent /= field_shape)) &
        error stop 'f2k shape: corr_w_hm_2'

      n = int(product(input(8)%extent(1:1)))

      if (n > 0) then

        if (.not. c_associated(input(8)%data)) error stop 'f2k data: corr_w_hm_2'

        call c_f_pointer(input(8)%data, rp, [n])
        value%corr_w_hm_2 = reshape(rp, shape(value%corr_w_hm_2))

      end if

    else

      error stop 'f2k missing fixed array: corr_w_hm_2'

    end if

    ! corr_chi_hm_1

    if (input(9)%rank /= 1) error stop 'f2k rank: corr_chi_hm_1'

    if (input(9)%allocated /= 0) then

      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%corr_chi_hm_1, 1, kind=c_int64_t)

      if (any(input(9)%extent /= field_shape)) &
        error stop 'f2k shape: corr_chi_hm_1'

      n = int(product(input(9)%extent(1:1)))

      if (n > 0) then

        if (.not. c_associated(input(9)%data)) error stop 'f2k data: corr_chi_hm_1'

        call c_f_pointer(input(9)%data, rp, [n])
        value%corr_chi_hm_1 = reshape(rp, shape(value%corr_chi_hm_1))

      end if

    else

      error stop 'f2k missing fixed array: corr_chi_hm_1'

    end if

    ! corr_chi_hm_2

    if (input(10)%rank /= 1) error stop 'f2k rank: corr_chi_hm_2'

    if (input(10)%allocated /= 0) then

      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%corr_chi_hm_2, 1, kind=c_int64_t)

      if (any(input(10)%extent /= field_shape)) &
        error stop 'f2k shape: corr_chi_hm_2'

      n = int(product(input(10)%extent(1:1)))

      if (n > 0) then

        if (.not. c_associated(input(10)%data)) error stop 'f2k data: corr_chi_hm_2'

        call c_f_pointer(input(10)%data, rp, [n])
        value%corr_chi_hm_2 = reshape(rp, shape(value%corr_chi_hm_2))

      end if

    else

      error stop 'f2k missing fixed array: corr_chi_hm_2'

    end if

    ! corr_eta_hm_1

    if (input(11)%rank /= 1) error stop 'f2k rank: corr_eta_hm_1'

    if (input(11)%allocated /= 0) then

      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%corr_eta_hm_1, 1, kind=c_int64_t)

      if (any(input(11)%extent /= field_shape)) &
        error stop 'f2k shape: corr_eta_hm_1'

      n = int(product(input(11)%extent(1:1)))

      if (n > 0) then

        if (.not. c_associated(input(11)%data)) error stop 'f2k data: corr_eta_hm_1'

        call c_f_pointer(input(11)%data, rp, [n])
        value%corr_eta_hm_1 = reshape(rp, shape(value%corr_eta_hm_1))

      end if

    else

      error stop 'f2k missing fixed array: corr_eta_hm_1'

    end if

    ! corr_eta_hm_2

    if (input(12)%rank /= 1) error stop 'f2k rank: corr_eta_hm_2'

    if (input(12)%allocated /= 0) then

      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%corr_eta_hm_2, 1, kind=c_int64_t)

      if (any(input(12)%extent /= field_shape)) &
        error stop 'f2k shape: corr_eta_hm_2'

      n = int(product(input(12)%extent(1:1)))

      if (n > 0) then

        if (.not. c_associated(input(12)%data)) error stop 'f2k data: corr_eta_hm_2'

        call c_f_pointer(input(12)%data, rp, [n])
        value%corr_eta_hm_2 = reshape(rp, shape(value%corr_eta_hm_2))

      end if

    else

      error stop 'f2k missing fixed array: corr_eta_hm_2'

    end if

    ! corr_hmx_hmy_1

    if (input(13)%rank /= 2) error stop 'f2k rank: corr_hmx_hmy_1'

    if (input(13)%allocated /= 0) then

      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%corr_hmx_hmy_1, 1, kind=c_int64_t)
      field_shape(2) = size(value%corr_hmx_hmy_1, 2, kind=c_int64_t)

      if (any(input(13)%extent /= field_shape)) &
        error stop 'f2k shape: corr_hmx_hmy_1'

      n = int(product(input(13)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(13)%data)) error stop 'f2k data: corr_hmx_hmy_1'

        call c_f_pointer(input(13)%data, rp, [n])
        value%corr_hmx_hmy_1 = reshape(rp, shape(value%corr_hmx_hmy_1))

      end if

    else

      error stop 'f2k missing fixed array: corr_hmx_hmy_1'

    end if

    ! corr_hmx_hmy_2

    if (input(14)%rank /= 2) error stop 'f2k rank: corr_hmx_hmy_2'

    if (input(14)%allocated /= 0) then

      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%corr_hmx_hmy_2, 1, kind=c_int64_t)
      field_shape(2) = size(value%corr_hmx_hmy_2, 2, kind=c_int64_t)

      if (any(input(14)%extent /= field_shape)) &
        error stop 'f2k shape: corr_hmx_hmy_2'

      n = int(product(input(14)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(14)%data)) error stop 'f2k data: corr_hmx_hmy_2'

        call c_f_pointer(input(14)%data, rp, [n])
        value%corr_hmx_hmy_2 = reshape(rp, shape(value%corr_hmx_hmy_2))

      end if

    else

      error stop 'f2k missing fixed array: corr_hmx_hmy_2'

    end if

    ! mu_Ncn_1

    if (input(15)%rank /= 0) error stop 'f2k rank: mu_Ncn_1'

    if (.not. c_associated(input(15)%data)) error stop 'f2k scalar: mu_Ncn_1'

    call c_f_pointer(input(15)%data, rp, [1])
    value%mu_Ncn_1 = rp(1)

    ! mu_Ncn_2

    if (input(16)%rank /= 0) error stop 'f2k rank: mu_Ncn_2'

    if (.not. c_associated(input(16)%data)) error stop 'f2k scalar: mu_Ncn_2'

    call c_f_pointer(input(16)%data, rp, [1])
    value%mu_Ncn_2 = rp(1)

    ! sigma_Ncn_1

    if (input(17)%rank /= 0) error stop 'f2k rank: sigma_Ncn_1'

    if (.not. c_associated(input(17)%data)) error stop 'f2k scalar: sigma_Ncn_1'

    call c_f_pointer(input(17)%data, rp, [1])
    value%sigma_Ncn_1 = rp(1)

    ! sigma_Ncn_2

    if (input(18)%rank /= 0) error stop 'f2k rank: sigma_Ncn_2'

    if (.not. c_associated(input(18)%data)) error stop 'f2k scalar: sigma_Ncn_2'

    call c_f_pointer(input(18)%data, rp, [1])
    value%sigma_Ncn_2 = rp(1)

  end subroutine to_native_hydromet_pdf_parameter

  subroutine from_native_hydromet_pdf_parameter(value, output)
    ! Description: Copy each hydromet_pdf_parameter field into C-owned output bytes.

    implicit none

    !--------------------------- Input Variables ---------------------------
    ! Call-local native argument; units and shape follow the source contract
    type(hydromet_pdf_parameter), intent(in) :: value

    !--------------------------- Output Variables ---------------------------
    type(wire_field), intent(out) :: output(*) ! Returned descriptors; C++ releases their C-owned buffers [-]

    !--------------------------- Local Variables ---------------------------
    ! Borrowed view of packed real data; units follow the native field
    real( kind = core_rknd ), pointer :: rp(:)
    integer :: n ! Number of packed elements [-]
    integer(c_int64_t) :: field_shape(4) ! Expected field extents in native dimension order [-]

    !--------------------------- Begin Code ---------------------------

    ! hm_1
    block
      n = size(value%hm_1)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%hm_1, 1, kind=c_int64_t)

      call allocate_field(1, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(1)) ! Out

      call c_f_pointer(output(1)%data, rp, [n])
      rp = reshape(value%hm_1, [n])
    end block

    ! hm_2
    block
      n = size(value%hm_2)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%hm_2, 1, kind=c_int64_t)

      call allocate_field(1, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(2)) ! Out

      call c_f_pointer(output(2)%data, rp, [n])
      rp = reshape(value%hm_2, [n])
    end block

    ! mu_hm_1
    block
      n = size(value%mu_hm_1)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%mu_hm_1, 1, kind=c_int64_t)

      call allocate_field(1, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(3)) ! Out

      call c_f_pointer(output(3)%data, rp, [n])
      rp = reshape(value%mu_hm_1, [n])
    end block

    ! mu_hm_2
    block
      n = size(value%mu_hm_2)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%mu_hm_2, 1, kind=c_int64_t)

      call allocate_field(1, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(4)) ! Out

      call c_f_pointer(output(4)%data, rp, [n])
      rp = reshape(value%mu_hm_2, [n])
    end block

    ! sigma_hm_1
    block
      n = size(value%sigma_hm_1)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%sigma_hm_1, 1, kind=c_int64_t)

      call allocate_field(1, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(5)) ! Out

      call c_f_pointer(output(5)%data, rp, [n])
      rp = reshape(value%sigma_hm_1, [n])
    end block

    ! sigma_hm_2
    block
      n = size(value%sigma_hm_2)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%sigma_hm_2, 1, kind=c_int64_t)

      call allocate_field(1, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(6)) ! Out

      call c_f_pointer(output(6)%data, rp, [n])
      rp = reshape(value%sigma_hm_2, [n])
    end block

    ! corr_w_hm_1
    block
      n = size(value%corr_w_hm_1)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%corr_w_hm_1, 1, kind=c_int64_t)

      call allocate_field(1, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(7)) ! Out

      call c_f_pointer(output(7)%data, rp, [n])
      rp = reshape(value%corr_w_hm_1, [n])
    end block

    ! corr_w_hm_2
    block
      n = size(value%corr_w_hm_2)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%corr_w_hm_2, 1, kind=c_int64_t)

      call allocate_field(1, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(8)) ! Out

      call c_f_pointer(output(8)%data, rp, [n])
      rp = reshape(value%corr_w_hm_2, [n])
    end block

    ! corr_chi_hm_1
    block
      n = size(value%corr_chi_hm_1)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%corr_chi_hm_1, 1, kind=c_int64_t)

      call allocate_field(1, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(9)) ! Out

      call c_f_pointer(output(9)%data, rp, [n])
      rp = reshape(value%corr_chi_hm_1, [n])
    end block

    ! corr_chi_hm_2
    block
      n = size(value%corr_chi_hm_2)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%corr_chi_hm_2, 1, kind=c_int64_t)

      call allocate_field(1, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(10)) ! Out

      call c_f_pointer(output(10)%data, rp, [n])
      rp = reshape(value%corr_chi_hm_2, [n])
    end block

    ! corr_eta_hm_1
    block
      n = size(value%corr_eta_hm_1)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%corr_eta_hm_1, 1, kind=c_int64_t)

      call allocate_field(1, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(11)) ! Out

      call c_f_pointer(output(11)%data, rp, [n])
      rp = reshape(value%corr_eta_hm_1, [n])
    end block

    ! corr_eta_hm_2
    block
      n = size(value%corr_eta_hm_2)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%corr_eta_hm_2, 1, kind=c_int64_t)

      call allocate_field(1, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(12)) ! Out

      call c_f_pointer(output(12)%data, rp, [n])
      rp = reshape(value%corr_eta_hm_2, [n])
    end block

    ! corr_hmx_hmy_1
    block
      n = size(value%corr_hmx_hmy_1)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%corr_hmx_hmy_1, 1, kind=c_int64_t)
      field_shape(2) = size(value%corr_hmx_hmy_1, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(13)) ! Out

      call c_f_pointer(output(13)%data, rp, [n])
      rp = reshape(value%corr_hmx_hmy_1, [n])
    end block

    ! corr_hmx_hmy_2
    block
      n = size(value%corr_hmx_hmy_2)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%corr_hmx_hmy_2, 1, kind=c_int64_t)
      field_shape(2) = size(value%corr_hmx_hmy_2, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(14)) ! Out

      call c_f_pointer(output(14)%data, rp, [n])
      rp = reshape(value%corr_hmx_hmy_2, [n])
    end block

    ! mu_Ncn_1

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0.0_core_rknd), 0, & ! In
                        output(15)) ! Out

    call c_f_pointer(output(15)%data, rp, [1])
    rp(1) = value%mu_Ncn_1

    ! mu_Ncn_2

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0.0_core_rknd), 0, & ! In
                        output(16)) ! Out

    call c_f_pointer(output(16)%data, rp, [1])
    rp(1) = value%mu_Ncn_2

    ! sigma_Ncn_1

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0.0_core_rknd), 0, & ! In
                        output(17)) ! Out

    call c_f_pointer(output(17)%data, rp, [1])
    rp(1) = value%sigma_Ncn_1

    ! sigma_Ncn_2

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0.0_core_rknd), 0, & ! In
                        output(18)) ! Out

    call c_f_pointer(output(18)%data, rp, [1])
    rp(1) = value%sigma_Ncn_2

  end subroutine from_native_hydromet_pdf_parameter

  subroutine f2k_roundtrip_hydromet_pdf_parameter_c(input, output) bind(C)
    ! Description: Exercise both converters with one local Fortran value.

    implicit none

    !--------------------------- Input Variables ---------------------------
    type(wire_field), intent(in) :: input(*) ! Caller-owned input field descriptors [-]

    !--------------------------- Output Variables ---------------------------
    type(wire_field), intent(out) :: output(*) ! Returned descriptors; C++ releases their C-owned buffers [-]

    !--------------------------- Local Variables ---------------------------
    ! Call-local native argument; units and shape follow the source contract
    type(hydromet_pdf_parameter) :: value

    !--------------------------- Begin Code ---------------------------

    call to_native_hydromet_pdf_parameter(input, value)
    call from_native_hydromet_pdf_parameter(value, output)

  end subroutine f2k_roundtrip_hydromet_pdf_parameter_c

  subroutine to_native_precipitation_fractions(input, value)
    ! Description: Copy each C field into a local precipitation_fractions value.

    implicit none

    !--------------------------- Input Variables ---------------------------
    type(wire_field), intent(in) :: input(*) ! Caller-owned input field descriptors [-]

    !--------------------------- Output Variables ---------------------------
    ! Call-local native argument; units and shape follow the source contract
    type(precipitation_fractions), intent(out) :: value

    !--------------------------- Local Variables ---------------------------
    integer(c_int), pointer :: ip(:) ! Borrowed view of packed integer data [-]
    ! Borrowed view of packed real data; units follow the native field
    real( kind = core_rknd ), pointer :: rp(:)
    integer :: n ! Number of packed elements [-]

    !--------------------------- Begin Code ---------------------------

    ! ngrdcol

    if (input(1)%rank /= 0) error stop 'f2k rank: ngrdcol'

    if (.not. c_associated(input(1)%data)) error stop 'f2k scalar: ngrdcol'

    call c_f_pointer(input(1)%data, ip, [1])
    value%ngrdcol = ip(1)

    ! nzt

    if (input(2)%rank /= 0) error stop 'f2k rank: nzt'

    if (.not. c_associated(input(2)%data)) error stop 'f2k scalar: nzt'

    call c_f_pointer(input(2)%data, ip, [1])
    value%nzt = ip(1)

    ! precip_frac

    if (input(3)%rank /= 2) error stop 'f2k rank: precip_frac'

    if (input(3)%allocated /= 0) then

      allocate(value%precip_frac(int(input(3)%extent(1)), int(input(3)%extent(2))))
      n = int(product(input(3)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(3)%data)) error stop 'f2k data: precip_frac'

        call c_f_pointer(input(3)%data, rp, [n])
        value%precip_frac = reshape(rp, shape(value%precip_frac))

      end if

    end if

    ! precip_frac_1

    if (input(4)%rank /= 2) error stop 'f2k rank: precip_frac_1'

    if (input(4)%allocated /= 0) then

      allocate(value%precip_frac_1(int(input(4)%extent(1)), int(input(4)%extent(2))))
      n = int(product(input(4)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(4)%data)) error stop 'f2k data: precip_frac_1'

        call c_f_pointer(input(4)%data, rp, [n])
        value%precip_frac_1 = reshape(rp, shape(value%precip_frac_1))

      end if

    end if

    ! precip_frac_2

    if (input(5)%rank /= 2) error stop 'f2k rank: precip_frac_2'

    if (input(5)%allocated /= 0) then

      allocate(value%precip_frac_2(int(input(5)%extent(1)), int(input(5)%extent(2))))
      n = int(product(input(5)%extent(1:2)))

      if (n > 0) then

        if (.not. c_associated(input(5)%data)) error stop 'f2k data: precip_frac_2'

        call c_f_pointer(input(5)%data, rp, [n])
        value%precip_frac_2 = reshape(rp, shape(value%precip_frac_2))

      end if

    end if

  end subroutine to_native_precipitation_fractions

  subroutine from_native_precipitation_fractions(value, output)
    ! Description: Copy each precipitation_fractions field into C-owned output bytes.

    implicit none

    !--------------------------- Input Variables ---------------------------
    ! Call-local native argument; units and shape follow the source contract
    type(precipitation_fractions), intent(in) :: value

    !--------------------------- Output Variables ---------------------------
    type(wire_field), intent(out) :: output(*) ! Returned descriptors; C++ releases their C-owned buffers [-]

    !--------------------------- Local Variables ---------------------------
    integer(c_int), pointer :: ip(:) ! Borrowed view of packed integer data [-]
    ! Borrowed view of packed real data; units follow the native field
    real( kind = core_rknd ), pointer :: rp(:)
    integer :: n ! Number of packed elements [-]
    integer(c_int64_t) :: field_shape(4) ! Expected field extents in native dimension order [-]

    !--------------------------- Begin Code ---------------------------

    ! ngrdcol

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(1)) ! Out

    call c_f_pointer(output(1)%data, ip, [1])
    ip(1) = value%ngrdcol

    ! nzt

    call allocate_field(0, [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t], & ! In
                        int(1, c_size_t) * c_sizeof(0_c_int), 0, & ! In
                        output(2)) ! Out

    call c_f_pointer(output(2)%data, ip, [1])
    ip(1) = value%nzt

    ! precip_frac

    if (allocated(value%precip_frac)) then

      n = size(value%precip_frac)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%precip_frac, 1, kind=c_int64_t)
      field_shape(2) = size(value%precip_frac, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(3)) ! Out

      call c_f_pointer(output(3)%data, rp, [n])
      rp = reshape(value%precip_frac, [n])

    else

      call empty_field(2, 0, & ! In
                       output(3)) ! Out

    end if

    ! precip_frac_1

    if (allocated(value%precip_frac_1)) then

      n = size(value%precip_frac_1)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%precip_frac_1, 1, kind=c_int64_t)
      field_shape(2) = size(value%precip_frac_1, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(4)) ! Out

      call c_f_pointer(output(4)%data, rp, [n])
      rp = reshape(value%precip_frac_1, [n])

    else

      call empty_field(2, 0, & ! In
                       output(4)) ! Out

    end if

    ! precip_frac_2

    if (allocated(value%precip_frac_2)) then

      n = size(value%precip_frac_2)
      field_shape = [1_c_int64_t, 1_c_int64_t, 1_c_int64_t, 1_c_int64_t]
      field_shape(1) = size(value%precip_frac_2, 1, kind=c_int64_t)
      field_shape(2) = size(value%precip_frac_2, 2, kind=c_int64_t)

      call allocate_field(2, field_shape, & ! In
                          int(n, c_size_t) * 1 * c_sizeof(0.0_core_rknd), 0, & ! In
                          output(5)) ! Out

      call c_f_pointer(output(5)%data, rp, [n])
      rp = reshape(value%precip_frac_2, [n])

    else

      call empty_field(2, 0, & ! In
                       output(5)) ! Out

    end if

  end subroutine from_native_precipitation_fractions

  subroutine f2k_roundtrip_precipitation_fractions_c(input, output) bind(C)
    ! Description: Exercise both converters with one local Fortran value.

    implicit none

    !--------------------------- Input Variables ---------------------------
    type(wire_field), intent(in) :: input(*) ! Caller-owned input field descriptors [-]

    !--------------------------- Output Variables ---------------------------
    type(wire_field), intent(out) :: output(*) ! Returned descriptors; C++ releases their C-owned buffers [-]

    !--------------------------- Local Variables ---------------------------
    ! Call-local native argument; units and shape follow the source contract
    type(precipitation_fractions) :: value

    !--------------------------- Begin Code ---------------------------

    call to_native_precipitation_fractions(input, value)
    call from_native_precipitation_fractions(value, output)

  end subroutine f2k_roundtrip_precipitation_fractions_c

  subroutine to_native_hydromet_pdf_parameter_collection(input, value)
    ! Description: Transfer the driver collection by component without saved state.

    implicit none

    !--------------------------- Input Variables ---------------------------
    type(wire_field), intent(in) :: input(*) ! Caller-owned input field descriptors [-]

    !--------------------------- Output Variables ---------------------------
    ! Call-local native argument; units and shape follow the source contract
    type(hydromet_pdf_parameter), allocatable, intent(out) :: value(:,:)

    !--------------------------- Local Variables ---------------------------
    ! Borrowed view of packed real data; units follow the native field
    real( kind = core_rknd ), pointer :: rp(:)
    integer :: ncol ! Number of columns [-]
    integer :: nzt ! Number of vertical levels [-]
    integer :: cells ! Number of column-level cells [-]
    integer :: n ! Number of packed elements [-]
    integer :: h ! Hydrometeor species index [-]
    integer :: j ! Packed array or field index [-]
    integer :: offset ! Offset of one species component in the packed array [-]
    integer(c_int64_t) :: field_shape(4) ! Expected field extents in native dimension order [-]

    !--------------------------- Begin Code ---------------------------

    ncol = int(input(1)%extent(1))
    nzt = int(input(1)%extent(2))
    allocate(value(ncol,nzt))
    cells = ncol*nzt

    ! hm_1
    field_shape = [int(ncol,c_int64_t), int(nzt,c_int64_t), 8_c_int64_t, 1_c_int64_t]
    n = cells*8

    if (input(1)%rank /= 3 .or. input(1)%allocated == 0) &
      error stop 'f2k collection allocation: hm_1'

    if (any(input(1)%extent /= field_shape)) error stop 'f2k collection shape'

    if (n > 0) then

      call c_f_pointer(input(1)%data,rp,[n])

      do h = 1, 8

        offset = (h-1)*cells
        value(:,:)%hm_1(h) = reshape(rp(offset+1:offset+cells),[ncol,nzt])

      end do

    end if

    ! hm_2
    field_shape = [int(ncol,c_int64_t), int(nzt,c_int64_t), 8_c_int64_t, 1_c_int64_t]
    n = cells*8

    if (input(2)%rank /= 3 .or. input(2)%allocated == 0) &
      error stop 'f2k collection allocation: hm_2'

    if (any(input(2)%extent /= field_shape)) error stop 'f2k collection shape'

    if (n > 0) then

      call c_f_pointer(input(2)%data,rp,[n])

      do h = 1, 8

        offset = (h-1)*cells
        value(:,:)%hm_2(h) = reshape(rp(offset+1:offset+cells),[ncol,nzt])

      end do

    end if

    ! mu_hm_1
    field_shape = [int(ncol,c_int64_t), int(nzt,c_int64_t), 8_c_int64_t, 1_c_int64_t]
    n = cells*8

    if (input(3)%rank /= 3 .or. input(3)%allocated == 0) &
      error stop 'f2k collection allocation: mu_hm_1'

    if (any(input(3)%extent /= field_shape)) error stop 'f2k collection shape'

    if (n > 0) then

      call c_f_pointer(input(3)%data,rp,[n])

      do h = 1, 8

        offset = (h-1)*cells
        value(:,:)%mu_hm_1(h) = reshape(rp(offset+1:offset+cells),[ncol,nzt])

      end do

    end if

    ! mu_hm_2
    field_shape = [int(ncol,c_int64_t), int(nzt,c_int64_t), 8_c_int64_t, 1_c_int64_t]
    n = cells*8

    if (input(4)%rank /= 3 .or. input(4)%allocated == 0) &
      error stop 'f2k collection allocation: mu_hm_2'

    if (any(input(4)%extent /= field_shape)) error stop 'f2k collection shape'

    if (n > 0) then

      call c_f_pointer(input(4)%data,rp,[n])

      do h = 1, 8

        offset = (h-1)*cells
        value(:,:)%mu_hm_2(h) = reshape(rp(offset+1:offset+cells),[ncol,nzt])

      end do

    end if

    ! sigma_hm_1
    field_shape = [int(ncol,c_int64_t), int(nzt,c_int64_t), 8_c_int64_t, 1_c_int64_t]
    n = cells*8

    if (input(5)%rank /= 3 .or. input(5)%allocated == 0) &
      error stop 'f2k collection allocation: sigma_hm_1'

    if (any(input(5)%extent /= field_shape)) error stop 'f2k collection shape'

    if (n > 0) then

      call c_f_pointer(input(5)%data,rp,[n])

      do h = 1, 8

        offset = (h-1)*cells
        value(:,:)%sigma_hm_1(h) = reshape(rp(offset+1:offset+cells),[ncol,nzt])

      end do

    end if

    ! sigma_hm_2
    field_shape = [int(ncol,c_int64_t), int(nzt,c_int64_t), 8_c_int64_t, 1_c_int64_t]
    n = cells*8

    if (input(6)%rank /= 3 .or. input(6)%allocated == 0) &
      error stop 'f2k collection allocation: sigma_hm_2'

    if (any(input(6)%extent /= field_shape)) error stop 'f2k collection shape'

    if (n > 0) then

      call c_f_pointer(input(6)%data,rp,[n])

      do h = 1, 8

        offset = (h-1)*cells
        value(:,:)%sigma_hm_2(h) = reshape(rp(offset+1:offset+cells),[ncol,nzt])

      end do

    end if

    ! corr_w_hm_1
    field_shape = [int(ncol,c_int64_t), int(nzt,c_int64_t), 8_c_int64_t, 1_c_int64_t]
    n = cells*8

    if (input(7)%rank /= 3 .or. input(7)%allocated == 0) &
      error stop 'f2k collection allocation: corr_w_hm_1'

    if (any(input(7)%extent /= field_shape)) error stop 'f2k collection shape'

    if (n > 0) then

      call c_f_pointer(input(7)%data,rp,[n])

      do h = 1, 8

        offset = (h-1)*cells
        value(:,:)%corr_w_hm_1(h) = reshape(rp(offset+1:offset+cells),[ncol,nzt])

      end do

    end if

    ! corr_w_hm_2
    field_shape = [int(ncol,c_int64_t), int(nzt,c_int64_t), 8_c_int64_t, 1_c_int64_t]
    n = cells*8

    if (input(8)%rank /= 3 .or. input(8)%allocated == 0) &
      error stop 'f2k collection allocation: corr_w_hm_2'

    if (any(input(8)%extent /= field_shape)) error stop 'f2k collection shape'

    if (n > 0) then

      call c_f_pointer(input(8)%data,rp,[n])

      do h = 1, 8

        offset = (h-1)*cells
        value(:,:)%corr_w_hm_2(h) = reshape(rp(offset+1:offset+cells),[ncol,nzt])

      end do

    end if

    ! corr_chi_hm_1
    field_shape = [int(ncol,c_int64_t), int(nzt,c_int64_t), 8_c_int64_t, 1_c_int64_t]
    n = cells*8

    if (input(9)%rank /= 3 .or. input(9)%allocated == 0) &
      error stop 'f2k collection allocation: corr_chi_hm_1'

    if (any(input(9)%extent /= field_shape)) error stop 'f2k collection shape'

    if (n > 0) then

      call c_f_pointer(input(9)%data,rp,[n])

      do h = 1, 8

        offset = (h-1)*cells
        value(:,:)%corr_chi_hm_1(h) = reshape(rp(offset+1:offset+cells),[ncol,nzt])

      end do

    end if

    ! corr_chi_hm_2
    field_shape = [int(ncol,c_int64_t), int(nzt,c_int64_t), 8_c_int64_t, 1_c_int64_t]
    n = cells*8

    if (input(10)%rank /= 3 .or. input(10)%allocated == 0) &
      error stop 'f2k collection allocation: corr_chi_hm_2'

    if (any(input(10)%extent /= field_shape)) error stop 'f2k collection shape'

    if (n > 0) then

      call c_f_pointer(input(10)%data,rp,[n])

      do h = 1, 8

        offset = (h-1)*cells
        value(:,:)%corr_chi_hm_2(h) = reshape(rp(offset+1:offset+cells),[ncol,nzt])

      end do

    end if

    ! corr_eta_hm_1
    field_shape = [int(ncol,c_int64_t), int(nzt,c_int64_t), 8_c_int64_t, 1_c_int64_t]
    n = cells*8

    if (input(11)%rank /= 3 .or. input(11)%allocated == 0) &
      error stop 'f2k collection allocation: corr_eta_hm_1'

    if (any(input(11)%extent /= field_shape)) error stop 'f2k collection shape'

    if (n > 0) then

      call c_f_pointer(input(11)%data,rp,[n])

      do h = 1, 8

        offset = (h-1)*cells
        value(:,:)%corr_eta_hm_1(h) = reshape(rp(offset+1:offset+cells),[ncol,nzt])

      end do

    end if

    ! corr_eta_hm_2
    field_shape = [int(ncol,c_int64_t), int(nzt,c_int64_t), 8_c_int64_t, 1_c_int64_t]
    n = cells*8

    if (input(12)%rank /= 3 .or. input(12)%allocated == 0) &
      error stop 'f2k collection allocation: corr_eta_hm_2'

    if (any(input(12)%extent /= field_shape)) error stop 'f2k collection shape'

    if (n > 0) then

      call c_f_pointer(input(12)%data,rp,[n])

      do h = 1, 8

        offset = (h-1)*cells
        value(:,:)%corr_eta_hm_2(h) = reshape(rp(offset+1:offset+cells),[ncol,nzt])

      end do

    end if

    ! corr_hmx_hmy_1
    field_shape = [int(ncol,c_int64_t), int(nzt,c_int64_t), 8_c_int64_t, 8_c_int64_t]
    n = cells*64

    if (input(13)%rank /= 4 .or. input(13)%allocated == 0) &
      error stop 'f2k collection allocation: corr_hmx_hmy_1'

    if (any(input(13)%extent /= field_shape)) error stop 'f2k collection shape'

    if (n > 0) then

      call c_f_pointer(input(13)%data,rp,[n])

      do j = 1, 8

        do h = 1, 8

          offset = ((j-1)*8+h-1)*cells
          value(:,:)%corr_hmx_hmy_1(h,j) = reshape(rp(offset+1:offset+cells),[ncol,nzt])

        end do

      end do

    end if

    ! corr_hmx_hmy_2
    field_shape = [int(ncol,c_int64_t), int(nzt,c_int64_t), 8_c_int64_t, 8_c_int64_t]
    n = cells*64

    if (input(14)%rank /= 4 .or. input(14)%allocated == 0) &
      error stop 'f2k collection allocation: corr_hmx_hmy_2'

    if (any(input(14)%extent /= field_shape)) error stop 'f2k collection shape'

    if (n > 0) then

      call c_f_pointer(input(14)%data,rp,[n])

      do j = 1, 8

        do h = 1, 8

          offset = ((j-1)*8+h-1)*cells
          value(:,:)%corr_hmx_hmy_2(h,j) = reshape(rp(offset+1:offset+cells),[ncol,nzt])

        end do

      end do

    end if

    ! mu_Ncn_1
    field_shape = [int(ncol,c_int64_t), int(nzt,c_int64_t), 1_c_int64_t, 1_c_int64_t]
    n = cells*1

    if (input(15)%rank /= 2 .or. input(15)%allocated == 0) &
      error stop 'f2k collection allocation: mu_Ncn_1'

    if (any(input(15)%extent /= field_shape)) error stop 'f2k collection shape'

    if (n > 0) then

      call c_f_pointer(input(15)%data,rp,[n])
      offset = 0
      value(:,:)%mu_Ncn_1 = reshape(rp(offset+1:offset+cells),[ncol,nzt])

    end if

    ! mu_Ncn_2
    field_shape = [int(ncol,c_int64_t), int(nzt,c_int64_t), 1_c_int64_t, 1_c_int64_t]
    n = cells*1

    if (input(16)%rank /= 2 .or. input(16)%allocated == 0) &
      error stop 'f2k collection allocation: mu_Ncn_2'

    if (any(input(16)%extent /= field_shape)) error stop 'f2k collection shape'

    if (n > 0) then

      call c_f_pointer(input(16)%data,rp,[n])
      offset = 0
      value(:,:)%mu_Ncn_2 = reshape(rp(offset+1:offset+cells),[ncol,nzt])

    end if

    ! sigma_Ncn_1
    field_shape = [int(ncol,c_int64_t), int(nzt,c_int64_t), 1_c_int64_t, 1_c_int64_t]
    n = cells*1

    if (input(17)%rank /= 2 .or. input(17)%allocated == 0) &
      error stop 'f2k collection allocation: sigma_Ncn_1'

    if (any(input(17)%extent /= field_shape)) error stop 'f2k collection shape'

    if (n > 0) then

      call c_f_pointer(input(17)%data,rp,[n])
      offset = 0
      value(:,:)%sigma_Ncn_1 = reshape(rp(offset+1:offset+cells),[ncol,nzt])

    end if

    ! sigma_Ncn_2
    field_shape = [int(ncol,c_int64_t), int(nzt,c_int64_t), 1_c_int64_t, 1_c_int64_t]
    n = cells*1

    if (input(18)%rank /= 2 .or. input(18)%allocated == 0) &
      error stop 'f2k collection allocation: sigma_Ncn_2'

    if (any(input(18)%extent /= field_shape)) error stop 'f2k collection shape'

    if (n > 0) then

      call c_f_pointer(input(18)%data,rp,[n])
      offset = 0
      value(:,:)%sigma_Ncn_2 = reshape(rp(offset+1:offset+cells),[ncol,nzt])

    end if

  end subroutine to_native_hydromet_pdf_parameter_collection

  subroutine from_native_hydromet_pdf_parameter_collection(value, output)
    ! Description: Transfer the driver collection by component without saved state.

    implicit none

    !--------------------------- Input Variables ---------------------------
    ! Call-local native argument; units and shape follow the source contract
    type(hydromet_pdf_parameter), intent(in) :: value(:,:)

    !--------------------------- Output Variables ---------------------------
    type(wire_field), intent(out) :: output(*) ! Returned descriptors; C++ releases their C-owned buffers [-]

    !--------------------------- Local Variables ---------------------------
    ! Borrowed view of packed real data; units follow the native field
    real( kind = core_rknd ), pointer :: rp(:)
    integer :: ncol ! Number of columns [-]
    integer :: nzt ! Number of vertical levels [-]
    integer :: cells ! Number of column-level cells [-]
    integer :: n ! Number of packed elements [-]
    integer :: h ! Hydrometeor species index [-]
    integer :: j ! Packed array or field index [-]
    integer :: offset ! Offset of one species component in the packed array [-]
    integer(c_int64_t) :: field_shape(4) ! Expected field extents in native dimension order [-]

    !--------------------------- Begin Code ---------------------------

    ncol = size(value,1)
    nzt = size(value,2)
    cells = ncol*nzt

    ! hm_1
    field_shape = [int(ncol,c_int64_t), int(nzt,c_int64_t), 8_c_int64_t, 1_c_int64_t]
    n = cells*8

    call allocate_field(3, field_shape, & ! In
                        int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                        output(1)) ! Out

    if (n > 0) then

      call c_f_pointer(output(1)%data,rp,[n])

      do h = 1, 8

        offset = (h-1)*cells
        rp(offset+1:offset+cells) = reshape(value(:,:)%hm_1(h),[cells])

      end do

    end if

    ! hm_2
    field_shape = [int(ncol,c_int64_t), int(nzt,c_int64_t), 8_c_int64_t, 1_c_int64_t]
    n = cells*8

    call allocate_field(3, field_shape, & ! In
                        int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                        output(2)) ! Out

    if (n > 0) then

      call c_f_pointer(output(2)%data,rp,[n])

      do h = 1, 8

        offset = (h-1)*cells
        rp(offset+1:offset+cells) = reshape(value(:,:)%hm_2(h),[cells])

      end do

    end if

    ! mu_hm_1
    field_shape = [int(ncol,c_int64_t), int(nzt,c_int64_t), 8_c_int64_t, 1_c_int64_t]
    n = cells*8

    call allocate_field(3, field_shape, & ! In
                        int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                        output(3)) ! Out

    if (n > 0) then

      call c_f_pointer(output(3)%data,rp,[n])

      do h = 1, 8

        offset = (h-1)*cells
        rp(offset+1:offset+cells) = reshape(value(:,:)%mu_hm_1(h),[cells])

      end do

    end if

    ! mu_hm_2
    field_shape = [int(ncol,c_int64_t), int(nzt,c_int64_t), 8_c_int64_t, 1_c_int64_t]
    n = cells*8

    call allocate_field(3, field_shape, & ! In
                        int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                        output(4)) ! Out

    if (n > 0) then

      call c_f_pointer(output(4)%data,rp,[n])

      do h = 1, 8

        offset = (h-1)*cells
        rp(offset+1:offset+cells) = reshape(value(:,:)%mu_hm_2(h),[cells])

      end do

    end if

    ! sigma_hm_1
    field_shape = [int(ncol,c_int64_t), int(nzt,c_int64_t), 8_c_int64_t, 1_c_int64_t]
    n = cells*8

    call allocate_field(3, field_shape, & ! In
                        int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                        output(5)) ! Out

    if (n > 0) then

      call c_f_pointer(output(5)%data,rp,[n])

      do h = 1, 8

        offset = (h-1)*cells
        rp(offset+1:offset+cells) = reshape(value(:,:)%sigma_hm_1(h),[cells])

      end do

    end if

    ! sigma_hm_2
    field_shape = [int(ncol,c_int64_t), int(nzt,c_int64_t), 8_c_int64_t, 1_c_int64_t]
    n = cells*8

    call allocate_field(3, field_shape, & ! In
                        int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                        output(6)) ! Out

    if (n > 0) then

      call c_f_pointer(output(6)%data,rp,[n])

      do h = 1, 8

        offset = (h-1)*cells
        rp(offset+1:offset+cells) = reshape(value(:,:)%sigma_hm_2(h),[cells])

      end do

    end if

    ! corr_w_hm_1
    field_shape = [int(ncol,c_int64_t), int(nzt,c_int64_t), 8_c_int64_t, 1_c_int64_t]
    n = cells*8

    call allocate_field(3, field_shape, & ! In
                        int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                        output(7)) ! Out

    if (n > 0) then

      call c_f_pointer(output(7)%data,rp,[n])

      do h = 1, 8

        offset = (h-1)*cells
        rp(offset+1:offset+cells) = reshape(value(:,:)%corr_w_hm_1(h),[cells])

      end do

    end if

    ! corr_w_hm_2
    field_shape = [int(ncol,c_int64_t), int(nzt,c_int64_t), 8_c_int64_t, 1_c_int64_t]
    n = cells*8

    call allocate_field(3, field_shape, & ! In
                        int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                        output(8)) ! Out

    if (n > 0) then

      call c_f_pointer(output(8)%data,rp,[n])

      do h = 1, 8

        offset = (h-1)*cells
        rp(offset+1:offset+cells) = reshape(value(:,:)%corr_w_hm_2(h),[cells])

      end do

    end if

    ! corr_chi_hm_1
    field_shape = [int(ncol,c_int64_t), int(nzt,c_int64_t), 8_c_int64_t, 1_c_int64_t]
    n = cells*8

    call allocate_field(3, field_shape, & ! In
                        int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                        output(9)) ! Out

    if (n > 0) then

      call c_f_pointer(output(9)%data,rp,[n])

      do h = 1, 8

        offset = (h-1)*cells
        rp(offset+1:offset+cells) = reshape(value(:,:)%corr_chi_hm_1(h),[cells])

      end do

    end if

    ! corr_chi_hm_2
    field_shape = [int(ncol,c_int64_t), int(nzt,c_int64_t), 8_c_int64_t, 1_c_int64_t]
    n = cells*8

    call allocate_field(3, field_shape, & ! In
                        int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                        output(10)) ! Out

    if (n > 0) then

      call c_f_pointer(output(10)%data,rp,[n])

      do h = 1, 8

        offset = (h-1)*cells
        rp(offset+1:offset+cells) = reshape(value(:,:)%corr_chi_hm_2(h),[cells])

      end do

    end if

    ! corr_eta_hm_1
    field_shape = [int(ncol,c_int64_t), int(nzt,c_int64_t), 8_c_int64_t, 1_c_int64_t]
    n = cells*8

    call allocate_field(3, field_shape, & ! In
                        int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                        output(11)) ! Out

    if (n > 0) then

      call c_f_pointer(output(11)%data,rp,[n])

      do h = 1, 8

        offset = (h-1)*cells
        rp(offset+1:offset+cells) = reshape(value(:,:)%corr_eta_hm_1(h),[cells])

      end do

    end if

    ! corr_eta_hm_2
    field_shape = [int(ncol,c_int64_t), int(nzt,c_int64_t), 8_c_int64_t, 1_c_int64_t]
    n = cells*8

    call allocate_field(3, field_shape, & ! In
                        int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                        output(12)) ! Out

    if (n > 0) then

      call c_f_pointer(output(12)%data,rp,[n])

      do h = 1, 8

        offset = (h-1)*cells
        rp(offset+1:offset+cells) = reshape(value(:,:)%corr_eta_hm_2(h),[cells])

      end do

    end if

    ! corr_hmx_hmy_1
    field_shape = [int(ncol,c_int64_t), int(nzt,c_int64_t), 8_c_int64_t, 8_c_int64_t]
    n = cells*64

    call allocate_field(4, field_shape, & ! In
                        int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                        output(13)) ! Out

    if (n > 0) then

      call c_f_pointer(output(13)%data,rp,[n])

      do j = 1, 8

        do h = 1, 8

          offset = ((j-1)*8+h-1)*cells
          rp(offset+1:offset+cells) = reshape(value(:,:)%corr_hmx_hmy_1(h,j),[cells])

        end do

      end do

    end if

    ! corr_hmx_hmy_2
    field_shape = [int(ncol,c_int64_t), int(nzt,c_int64_t), 8_c_int64_t, 8_c_int64_t]
    n = cells*64

    call allocate_field(4, field_shape, & ! In
                        int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                        output(14)) ! Out

    if (n > 0) then

      call c_f_pointer(output(14)%data,rp,[n])

      do j = 1, 8

        do h = 1, 8

          offset = ((j-1)*8+h-1)*cells
          rp(offset+1:offset+cells) = reshape(value(:,:)%corr_hmx_hmy_2(h,j),[cells])

        end do

      end do

    end if

    ! mu_Ncn_1
    field_shape = [int(ncol,c_int64_t), int(nzt,c_int64_t), 1_c_int64_t, 1_c_int64_t]
    n = cells*1

    call allocate_field(2, field_shape, & ! In
                        int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                        output(15)) ! Out

    if (n > 0) then

      call c_f_pointer(output(15)%data,rp,[n])
      offset = 0
      rp(offset+1:offset+cells) = reshape(value(:,:)%mu_Ncn_1,[cells])

    end if

    ! mu_Ncn_2
    field_shape = [int(ncol,c_int64_t), int(nzt,c_int64_t), 1_c_int64_t, 1_c_int64_t]
    n = cells*1

    call allocate_field(2, field_shape, & ! In
                        int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                        output(16)) ! Out

    if (n > 0) then

      call c_f_pointer(output(16)%data,rp,[n])
      offset = 0
      rp(offset+1:offset+cells) = reshape(value(:,:)%mu_Ncn_2,[cells])

    end if

    ! sigma_Ncn_1
    field_shape = [int(ncol,c_int64_t), int(nzt,c_int64_t), 1_c_int64_t, 1_c_int64_t]
    n = cells*1

    call allocate_field(2, field_shape, & ! In
                        int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                        output(17)) ! Out

    if (n > 0) then

      call c_f_pointer(output(17)%data,rp,[n])
      offset = 0
      rp(offset+1:offset+cells) = reshape(value(:,:)%sigma_Ncn_1,[cells])

    end if

    ! sigma_Ncn_2
    field_shape = [int(ncol,c_int64_t), int(nzt,c_int64_t), 1_c_int64_t, 1_c_int64_t]
    n = cells*1

    call allocate_field(2, field_shape, & ! In
                        int(n,c_size_t)*c_sizeof(0._core_rknd), 0, & ! In
                        output(18)) ! Out

    if (n > 0) then

      call c_f_pointer(output(18)%data,rp,[n])
      offset = 0
      rp(offset+1:offset+cells) = reshape(value(:,:)%sigma_Ncn_2,[cells])

    end if

  end subroutine from_native_hydromet_pdf_parameter_collection

  subroutine f2k_roundtrip_hydromet_pdf_parameter_collection_c(input,output) bind(C)
    ! Description: Check both collection converters with a local native array.

    implicit none

    !--------------------------- Input Variables ---------------------------
    type(wire_field), intent(in) :: input(*) ! Caller-owned input field descriptors [-]

    !--------------------------- Output Variables ---------------------------
    type(wire_field), intent(out) :: output(*) ! Returned descriptors; C++ releases their C-owned buffers [-]

    !--------------------------- Local Variables ---------------------------
    ! Call-local native argument; units and shape follow the source contract
    type(hydromet_pdf_parameter), allocatable :: value(:,:)

    !--------------------------- Begin Code ---------------------------

    call to_native_hydromet_pdf_parameter_collection(input,value)
    call from_native_hydromet_pdf_parameter_collection(value,output)

  end subroutine f2k_roundtrip_hydromet_pdf_parameter_collection_c

end module f2k_derived_type_bridge
