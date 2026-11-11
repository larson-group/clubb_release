// Caller-owned MT19937 stream translated from mt95.F90.
// Initialization accepts scalar/array seeds or state snapshots. Six draw APIs
// preserve Native signed-integer conversion and Fortran array traversal order;
// genrand_srepr preserves the Native base-62 restart representation.
//
// Native algorithm: Takuji Nishimura and Makoto Matsumoto (1997,2002).
// Fortran conversion/enhancements: Jose Rui Faustino de Sousa (2002,2003).
// This library is free software, without warranty; see src/SILHS/mt95.F90.
#pragma once
#include <Kokkos_Core.hpp>
#include <array>
#include <bit>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace clubb_kokkos {
using genrand_intg = std::int32_t;
using genrand_real = double;

class genrand_state;
struct genrand_srepr {
  std::string repr = std::string(4375, ' ');
  genrand_srepr& operator=(const genrand_state& stt);
};

class mt95;
class genrand_state {
public:
  genrand_state& operator=(const genrand_srepr& rpr);

private:
  friend class mt95;
  friend void genrand_load_state(genrand_state&, const genrand_srepr&);
  friend void genrand_dump_state(genrand_srepr&, const genrand_state&);
  bool ini         = false;
  genrand_intg cnt = 625;
  std::array<std::uint32_t, 624> val{};
};

void genrand_load_state(genrand_state& stt, const genrand_srepr& rpr);
void genrand_dump_state(genrand_srepr& rpr, const genrand_state& stt);

class mt95 {
public:
  // The Native SAVE/threadprivate value belongs to this case-owned object.
  // Copies preserve the entire stream by value; independent cases never alias.
  void genrand_init_api(const genrand_state* put = nullptr, genrand_state* get = nullptr);
  void genrand_init_api(genrand_intg put);
  void genrand_init_api(const genrand_intg* put, int key_length);
  template <class View> void genrand_init_api(View put) {
    // Seeds are control data: stage once, preserving first-index traversal.
    auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace{}, put);
    std::vector<genrand_intg> seed(put.extent(0));
    for (std::size_t i = 0; i < seed.size(); ++i)
      seed[i] = host(i);
    genrand_init_api(seed.data(), static_cast<int>(seed.size()));
  }
  void genrand_int32(genrand_intg& y);
  void genrand_int31(genrand_intg& y);
  void genrand_real1(genrand_real& r);
  void genrand_real2(genrand_real& r);
  void genrand_real3(genrand_real& r);
  void genrand_res53(genrand_real& r);

  // These thin rank adapters mirror Native's rank-1 through rank-7 overloads.
  // Generate the complete ordered stream on the host and transfer the array
  // once; no scalar host/device bridge is introduced for individual draws.
#define CLUBB_MT95_ARRAY_API(NAME)                                                                 \
  template <class View> void NAME(View values) {                                                   \
    static_assert(View::rank >= 1 && View::rank <= 7);                                             \
    using Host =                                                                                   \
        Kokkos::View<typename View::non_const_data_type, Kokkos::LayoutLeft, Kokkos::HostSpace>;   \
    auto host = allocate_rank<Host>(values, std::make_index_sequence<View::rank>{});               \
    NAME##_rank(host);                                                                             \
    using Packed = Kokkos::View<typename View::non_const_data_type, Kokkos::LayoutLeft,            \
                                typename View::memory_space>;                                      \
    auto packed  = allocate_rank<Packed>(values, std::make_index_sequence<View::rank>{});          \
    Kokkos::deep_copy(packed, host);                                                               \
    Kokkos::deep_copy(values, packed);                                                             \
  }
  CLUBB_MT95_ARRAY_API(genrand_int32)
  CLUBB_MT95_ARRAY_API(genrand_int31)
  CLUBB_MT95_ARRAY_API(genrand_real1)
  CLUBB_MT95_ARRAY_API(genrand_real2)
  CLUBB_MT95_ARRAY_API(genrand_real3)
  CLUBB_MT95_ARRAY_API(genrand_res53)
#undef CLUBB_MT95_ARRAY_API

private:
  genrand_state state;
  static std::uint32_t uiadd(std::uint32_t a, std::uint32_t b);
  static std::uint32_t uisub(std::uint32_t a, std::uint32_t b);
  static std::uint32_t uimlt(std::uint32_t a, std::uint32_t b);
  static std::uint32_t uidiv(std::uint32_t a, std::uint32_t b);
  static std::uint32_t uimod(std::uint32_t a, std::uint32_t b);
  void init_by_type(const genrand_state* put, genrand_state* get);
  void init_by_scalar(genrand_intg put);
  void init_by_array(const genrand_intg* put, int key_length);
  void next_state();
  static std::string genrand_encode(std::uint32_t val);
  static std::uint32_t genrand_decode(const std::string& chr);
  friend void genrand_load_state(genrand_state&, const genrand_srepr&);
  friend void genrand_dump_state(genrand_srepr&, const genrand_state&);
  template <class Result, class View, std::size_t... I>
  static Result allocate_rank(const View& view, std::index_sequence<I...>) {
    // Shape adapter for the Native seven supported array ranks.
    return Result(Kokkos::view_alloc(Kokkos::WithoutInitializing, "mt95_draws"), view.extent(I)...);
  }
#define CLUBB_MT95_ORDERED_RANK(NAME)                                                              \
  template <class View> void NAME##_rank(View values) {                                            \
    for (std::size_t i = 0; i < values.size(); ++i)                                                \
      NAME(values.data()[i]);                                                                      \
  }
  // LayoutLeft packed storage has precisely the recursive Native first-index
  // order. The six adapters replace only identical rank-wrapper boilerplate.
  CLUBB_MT95_ORDERED_RANK(genrand_int32)
  CLUBB_MT95_ORDERED_RANK(genrand_int31)
  CLUBB_MT95_ORDERED_RANK(genrand_real1)
  CLUBB_MT95_ORDERED_RANK(genrand_real2)
  CLUBB_MT95_ORDERED_RANK(genrand_real3)
  CLUBB_MT95_ORDERED_RANK(genrand_res53)
#undef CLUBB_MT95_ORDERED_RANK
};
} // namespace clubb_kokkos
