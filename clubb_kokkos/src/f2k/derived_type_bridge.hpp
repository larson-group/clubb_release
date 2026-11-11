#pragma once

#include "CLUBB_core/array_index.hpp"
#include "CLUBB_core/corr_varnce_module.hpp"
#include "CLUBB_core/err_info_type_module.hpp"
#include "CLUBB_core/grid_class.hpp"
#include "CLUBB_core/hydromet_pdf_parameter_module.hpp"
#include "CLUBB_core/model_flags.hpp"
#include "CLUBB_core/parameters_tunable.hpp"
#include "CLUBB_core/pdf_parameter_module.hpp"
#include "SILHS/parameters_silhs.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace f2k {

// C-compatible field descriptor. Array bytes are packed in Fortran index order.
struct WireField {
  void* data{};
  std::int64_t extent[4]{};
  std::int32_t rank{};
  std::int32_t allocated{};
  std::int32_t char_len{};
};

class Transfer {
 public:
  explicit Transfer(std::size_t count, bool output = false);
  ~Transfer();
  Transfer(Transfer&& other) noexcept;
  Transfer(const Transfer&) = delete;
  Transfer& operator=(const Transfer&) = delete;
  Transfer& operator=(Transfer&&) = delete;

  WireField* data() { return fields.data(); }
  const WireField* data() const { return fields.data(); }
  std::vector<WireField> fields;
  std::vector<std::vector<unsigned char>> buffers;

 private:
  bool owns_output_;
};

using Real1D = Kokkos::View<const clubb_kokkos::real_t*, Kokkos::LayoutLeft>;
using Real2D = Kokkos::View<const clubb_kokkos::real_t**, Kokkos::LayoutLeft>;

// Pack standalone array arguments using the same host staging as UDT fields.
WireField pack_real_array(const Real1D& value, Transfer& transfer);
WireField pack_real_array(const Real2D& value, Transfer& transfer);
WireField pack_strings(const std::vector<std::string>& value, Transfer& transfer);
WireField pack_char_array(const Kokkos::View<char**, Kokkos::LayoutLeft>& value,
                          Transfer& transfer);

// The Fortran wrapper around an unported routine calls to_native_* and
// from_native_* with these transfers during that same call.
#define F2K_DECLARE_CONVERTER(Type)                                            \
  Transfer to_fortran(const clubb_kokkos::Type& value);                       \
  Transfer make_output(const clubb_kokkos::Type& value);                      \
  void to_kokkos(clubb_kokkos::Type& value, const Transfer& output);           \
  void roundtrip(clubb_kokkos::Type& value)

F2K_DECLARE_CONVERTER(grid);
F2K_DECLARE_CONVERTER(sclr_idx_type);
F2K_DECLARE_CONVERTER(clubb_config_flags_type);
F2K_DECLARE_CONVERTER(nu_vertical_res_dep);
F2K_DECLARE_CONVERTER(pdf_parameter);
F2K_DECLARE_CONVERTER(implicit_coefs_terms);
F2K_DECLARE_CONVERTER(err_info_type);
F2K_DECLARE_CONVERTER(hm_metadata_type);
F2K_DECLARE_CONVERTER(silhs_config_flags_type);
F2K_DECLARE_CONVERTER(hydromet_pdf_parameter);
F2K_DECLARE_CONVERTER(hydromet_pdf_parameter_collection);
F2K_DECLARE_CONVERTER(precipitation_fractions);

#undef F2K_DECLARE_CONVERTER

}  // namespace f2k
