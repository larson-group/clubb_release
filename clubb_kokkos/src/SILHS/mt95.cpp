// MT19937 initialization, restart serialization and Native signed draw formulas.
// See mt95.hpp and the original src/SILHS/mt95.F90 attribution/license.
#include "SILHS/mt95.hpp"
#include <algorithm>

namespace clubb_kokkos {
namespace {
constexpr int n                     = 624;
constexpr int m                     = 397;
constexpr genrand_intg default_seed = 5489;
constexpr double p231               = 2147483648.0;
constexpr double p232               = 4294967296.0;
constexpr double p232_1             = p232 - 1.0;
constexpr double pi232              = 1.0 / p232;
constexpr double pi232_1            = 1.0 / p232_1;
constexpr double pi227              = 1.0 / 134217728.0;
constexpr double pi253              = 1.0 / 9007199254740992.0;
constexpr double p231d232_1         = p231 / p232_1;
constexpr double p231_5d232         = (p231 + 0.5) / p232;
constexpr char alph[]        = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
constexpr std::uint32_t alps = 62;
constexpr std::size_t clen   = 4375;
} // namespace

// Unsigned C++ arithmetic gives the defined modulo-2^32 operations which the
// Native split-word add/subtract/multiply helpers implement without signed overflow. These are
// integer representation adaptations, not floating-point BFB workarounds.
std::uint32_t mt95::uiadd(std::uint32_t a, std::uint32_t b) { return a + b; }
std::uint32_t mt95::uisub(std::uint32_t a, std::uint32_t b) { return a - b; }
std::uint32_t mt95::uimlt(std::uint32_t a, std::uint32_t b) { return a * b; }
std::uint32_t mt95::uidiv(std::uint32_t a, std::uint32_t b) {
  // Retain Native's branch algorithm and signed comparisons of word bits.
  // Only divisor62 is reached by restart encoding; do not silently generalize
  // the private helper to a different unsigned-division contract.
  if (a & 0x80000000U) {
    if (b & 0x80000000U) {
      if (std::bit_cast<genrand_intg>(a) < std::bit_cast<genrand_intg>(b))
        return 0;
      return 1;
    }
    const std::uint32_t dl = ((a >> 1) / b) << 1;
    const std::uint32_t rl = uisub(a, uimlt(b, dl));
    if (std::bit_cast<genrand_intg>(rl) < std::bit_cast<genrand_intg>(b))
      return dl;
    return uiadd(dl, 1);
  }
  if (b & 0x80000000U)
    return 0;
  return a / b;
}
std::uint32_t mt95::uimod(std::uint32_t a, std::uint32_t b) {
  // The corresponding Native remainder uses exactly the same partition.
  if (a & 0x80000000U) {
    if (b & 0x80000000U) {
      if (std::bit_cast<genrand_intg>(a) < std::bit_cast<genrand_intg>(b))
        return a;
      return uisub(a, b);
    }
    const std::uint32_t dl = ((a >> 1) / b) << 1;
    const std::uint32_t rl = uisub(a, uimlt(b, dl));
    if (std::bit_cast<genrand_intg>(rl) < std::bit_cast<genrand_intg>(b))
      return rl;
    return uisub(rl, b);
  }
  if (b & 0x80000000U)
    return a;
  return a % b;
}

void mt95::init_by_type(const genrand_state* put, genrand_state* get) {
  // Native PUT takes precedence over GET, including an uninitialized PUT.
  if (put) {
    if (put->ini)
      state = *put;
  } else if (get) {
    if (!state.ini)
      init_by_scalar(default_seed);
    *get = state;
  } else
    init_by_scalar(default_seed);
}

void mt95::init_by_scalar(genrand_intg put) {
  // Initialize 624 words from a scalar seed (Knuth TAOCP2,3rded,p.106).
  constexpr std::uint32_t mult_a = 1812433253U;
  state.ini                      = true;
  state.val[0]                   = std::bit_cast<std::uint32_t>(put);
  for (int i = 1; i < n; ++i) {
    state.val[i] = state.val[i - 1] ^ (state.val[i - 1] >> 30);
    state.val[i] = uimlt(state.val[i], mult_a);
    state.val[i] = uiadd(state.val[i], std::uint32_t(i));
  }
  state.cnt = n + 1;
}

void mt95::init_by_array(const genrand_intg* put, int key_length) {
  // Nonlinear array seeding follows the original 2002 MT19937 improvement.
  constexpr genrand_intg seed_d  = 19650218;
  constexpr std::uint32_t mult_a = 1664525U;
  constexpr std::uint32_t mult_b = 1566083941U;
  constexpr std::uint32_t msb1_d = 0x80000000U;
  init_by_scalar(seed_d);
  int i = 1, j = 0;
  for (int k = std::max(n, key_length); k >= 1; --k) {
    std::uint32_t tp = state.val[i - 1] ^ (state.val[i - 1] >> 30);
    tp               = uimlt(tp, mult_a);
    state.val[i]     = state.val[i] ^ tp;
    state.val[i] =
        uiadd(state.val[i], uiadd(std::bit_cast<std::uint32_t>(put[j]), std::uint32_t(j)));
    ++i;
    ++j;
    if (i >= n) {
      state.val[0] = state.val[n - 1];
      i            = 1;
    }
    if (j >= key_length)
      j = 0;
  }
  for (int k = n - 1; k >= 1; --k) {
    std::uint32_t tp = state.val[i - 1] ^ (state.val[i - 1] >> 30);
    tp               = uimlt(tp, mult_b);
    state.val[i]     = state.val[i] ^ tp;
    state.val[i]     = uisub(state.val[i], std::uint32_t(i));
    ++i;
    if (i >= n) {
      state.val[0] = state.val[n - 1];
      i            = 1;
    }
  }
  state.val[0] = msb1_d;
}

void mt95::next_state() {
  // Twist in Native order; the second loop intentionally reads updated words.
  constexpr std::uint32_t matrix_a = 0x9908b0dfU;
  if (!state.ini)
    init_by_scalar(default_seed);
  for (int i = 0; i < n - m; ++i) {
    const std::uint32_t mld = (state.val[i + 1] & 0x7fffffffU) | (state.val[i] & 0x80000000U);
    state.val[i]            = state.val[i + m] ^ (mld >> 1);
    if (state.val[i + 1] & 1U)
      state.val[i] ^= matrix_a;
  }
  for (int i = n - m; i < n - 1; ++i) {
    const std::uint32_t mld = (state.val[i + 1] & 0x7fffffffU) | (state.val[i] & 0x80000000U);
    state.val[i]            = state.val[i + m - n] ^ (mld >> 1);
    if (state.val[i + 1] & 1U)
      state.val[i] ^= matrix_a;
  }
  const std::uint32_t mld = (state.val[0] & 0x7fffffffU) | (state.val[n - 1] & 0x80000000U);
  state.val[n - 1]        = state.val[m - 1] ^ (mld >> 1);
  if (state.val[0] & 1U)
    state.val[n - 1] ^= matrix_a;
  state.cnt = 1;
}

std::string mt95::genrand_encode(std::uint32_t val) {
  // Little-endian base-62 digits, as in Native's fixed-width text value.
  std::string chr;
  std::uint32_t d = val;
  do {
    const std::uint32_t m = uimod(d, alps);
    chr += alph[m];
    d = uidiv(d, alps);
  } while (d != 0);
  return chr;
}

std::uint32_t mt95::genrand_decode(const std::string& chr) {
  // Native ignores characters outside its alphabet without advancing exponent.
  std::uint32_t val = 0, e = 1;
  for (char c : chr) {
    const auto p = std::string(alph).find(c);
    if (p != std::string::npos) {
      val = uiadd(val, uimlt(std::uint32_t(p), e));
      e   = uimlt(e, alps);
    }
  }
  return val;
}

void genrand_load_state(genrand_state& stt, const genrand_srepr& rpr) {
  // Decode the source restart record; its counter remains one based.
  std::size_t first = 0;
  int i             = 0;
  for (;;) {
    const auto j = rpr.repr.find('&', first);
    if (j == std::string::npos)
      break;
    stt.val[i] = mt95::genrand_decode(rpr.repr.substr(first, j - first));
    ++i;
    first = j + 1;
  }
  stt.cnt = std::bit_cast<genrand_intg>(mt95::genrand_decode(rpr.repr.substr(first)));
  stt.ini = true;
}

void genrand_dump_state(genrand_srepr& rpr, const genrand_state& stt) {
  // Encode every word plus the counter, padded to Native character length.
  rpr.repr.clear();
  for (int i = 0; i < n; ++i) {
    rpr.repr += mt95::genrand_encode(stt.val[i]);
    rpr.repr += '&';
  }
  rpr.repr += mt95::genrand_encode(std::bit_cast<std::uint32_t>(stt.cnt));
  rpr.repr.resize(clen, ' ');
}

genrand_state& genrand_state::operator=(const genrand_srepr& rpr) {
  // Native generic assignment decodes the public restart representation.
  genrand_load_state(*this, rpr);
  return *this;
}

genrand_srepr& genrand_srepr::operator=(const genrand_state& stt) {
  // Native generic assignment encodes the private generator state.
  genrand_dump_state(*this, stt);
  return *this;
}

void mt95::genrand_int32(genrand_intg& y) {
  // Return the exact signed bit representation of one tempered 32-bit word.
  if (state.cnt > n)
    next_state();
  std::uint32_t word = state.val[state.cnt - 1];
  ++state.cnt;
  word ^= word >> 11;
  word ^= (word << 7) & 0x9d2c5680U;
  word ^= (word << 15) & 0xefc60000U;
  word ^= word >> 18;
  y = std::bit_cast<genrand_intg>(word);
}

void mt95::genrand_int31(genrand_intg& y) {
  // Native ISHFT is a logical right shift, including negative bit patterns.
  genrand_int32(y);
  y = genrand_intg(std::bit_cast<std::uint32_t>(y) >> 1);
}

void mt95::genrand_real1(genrand_real& r) {
  // Signed-word conversion used by Native; this differs from unsigned /2^32.
  genrand_intg a;
  genrand_int32(a);
  r = double(a) * pi232_1 + p231d232_1;
}
void mt95::genrand_real2(genrand_real& r) {
  genrand_intg a;
  genrand_int32(a);
  r = double(a) * pi232 + 0.5;
}
void mt95::genrand_real3(genrand_real& r) {
  genrand_intg a;
  genrand_int32(a);
  r = double(a) * pi232 + p231_5d232;
}
void mt95::genrand_res53(genrand_real& r) {
  // Two draws per real, with the Native 27-/26-bit split and operation order.
  genrand_intg a, b;
  genrand_int32(a);
  genrand_int32(b);
  a = genrand_intg(std::bit_cast<std::uint32_t>(a) >> 5);
  b = genrand_intg(std::bit_cast<std::uint32_t>(b) >> 6);
  r = double(a) * pi227 + double(b) * pi253;
}
void mt95::genrand_init_api(const genrand_state* put, genrand_state* get) {
  init_by_type(put, get);
}
void mt95::genrand_init_api(genrand_intg put) { init_by_scalar(put); }
void mt95::genrand_init_api(const genrand_intg* put, int key_length) {
  init_by_array(put, key_length);
}
} // namespace clubb_kokkos
