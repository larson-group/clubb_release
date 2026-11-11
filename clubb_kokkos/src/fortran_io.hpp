// Shared lexical rules for Fortran-style text and namelist input.
//
// tokenize preserves records and quoted values; typed conversion accepts Native
// decimal/D-exponent spelling without retaining any parser or file state.
// Use character_values for list-directed character READs, and unterminated_tail
// when an unused trailing quote must not invalidate an already completed READ.
// normalized_number canonicalizes real spelling; lower folds identifiers.
// parse_integer checks strict metadata syntax; parse_value converts typed READ
// values and preserves its destination on failure. INTEGER input also accepts
// real spelling, truncating toward zero as in the installed NVHPC reader.
#pragma once

#include <cerrno>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <limits>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

namespace clubb_kokkos::fortran_io {
enum class Kind {
  value,
  comma,
  slash,
  star,
  equal,
  left_paren,
  right_paren,
  colon,
  percent,
  group,
  end_record
};
struct Token {
  Kind kind;
  std::string text;
  bool quoted       = false;
  std::size_t line  = 0;
  std::size_t begin = 0, end = 0;
};

std::vector<Token> tokenize(std::string_view source, bool character_values = false,
                            bool* unterminated_tail = nullptr);
std::string normalized_number(std::string_view value);
std::string lower(std::string_view value);

template <class T> bool parse_integer(std::string_view text, T& value) {
  // Strict integer syntax for metadata such as subscripts/repetition counts.
  static_assert(std::is_integral_v<T> && !std::is_same_v<T, bool>);
  const std::string token(text);
  if (token.empty())
    return false;
  char* end              = nullptr;
  errno                  = 0;
  const long long parsed = std::strtoll(token.c_str(), &end, 10);
  if (end != token.c_str() + token.size() || end == token.c_str() || errno == ERANGE ||
      parsed < std::numeric_limits<T>::lowest() || parsed > std::numeric_limits<T>::max())
    return false;
  value = static_cast<T>(parsed);
  return true;
}

template <class T> bool parse_value(std::string_view text, T& value) {
  // Convert one decoded list-directed token; malformed or overflowing input
  // returns false, retaining the caller's initialized value.
  if constexpr (std::is_same_v<T, std::string>) {
    value = std::string(text);
    return true;
  } else if constexpr (std::is_same_v<T, bool>) {
    std::string token = lower(text);
    if (!token.empty() && token.front() == '.')
      token.erase(token.begin());
    if (token.empty() || (token.front() != 't' && token.front() != 'f'))
      return false;
    value = token.front() == 't';
    return true;
  } else if constexpr (std::is_integral_v<T>) {
    // NVHPC list-directed and namelist INTEGER reads accept real spelling.
    // Preserve exact integer conversion first; real values truncate toward zero.
    if (parse_integer(text, value))
      return true;
    long double real;
    if (!parse_value(text, real) || !std::isfinite(real))
      return false;
    const auto truncated = std::trunc(real);
    if (truncated < std::numeric_limits<T>::lowest() || truncated > std::numeric_limits<T>::max())
      return false;
    value = static_cast<T>(truncated);
    return true;
  } else {
    static_assert(std::is_floating_point_v<T>, "Unsupported Fortran input value");
    const std::string token = normalized_number(text);
    if (token.empty())
      return false;
    char* end = nullptr;
    errno     = 0;
    T parsed;
    if constexpr (std::is_same_v<T, float>)
      parsed = std::strtof(token.c_str(), &end);
    else if constexpr (std::is_same_v<T, double>)
      parsed = std::strtod(token.c_str(), &end);
    else
      parsed = std::strtold(token.c_str(), &end);
    if (end != token.c_str() + token.size() || end == token.c_str() ||
        (errno == ERANGE && !std::isfinite(parsed)))
      return false;
    value = parsed;
    return true;
  }
}
} // namespace clubb_kokkos::fortran_io
