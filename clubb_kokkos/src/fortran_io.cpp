// Lex Fortran input while retaining external-record and quotation boundaries.
#include "fortran_io.hpp"

#include <algorithm>
#include <cctype>
#include <stdexcept>

namespace clubb_kokkos::fortran_io {
std::string lower(std::string_view value) {
  // Fold source identifiers/logical tokens without changing stored text values.
  std::string result(value);
  for (char& c : result)
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  return result;
}

std::string normalized_number(std::string_view value) {
  // Normalize D exponents and Native exponent-letter omission, e.g. 1.5-3.
  // Only decimal numeric syntax and the Native Inf/NaN tokens are accepted.
  std::string result(value);
  const auto begin = result.find_first_not_of(" \t\r\n");
  if (begin == std::string::npos)
    return {};
  result        = result.substr(begin, result.find_last_not_of(" \t\r\n") - begin + 1);
  auto spelling = lower(result);
  if (spelling.front() == '+' || spelling.front() == '-')
    spelling.erase(spelling.begin());
  if (spelling == "inf" || spelling == "infinity" || spelling == "nan")
    return result;
  for (char& c : result)
    if (c == 'd' || c == 'D')
      c = 'e';
  if (result.find_first_of("eE") == std::string::npos) {
    for (std::size_t i = 1; i < result.size(); ++i)
      if ((result[i] == '+' || result[i] == '-') &&
          (std::isdigit(static_cast<unsigned char>(result[i - 1])) || result[i - 1] == '.')) {
        result.insert(i, "e");
        break;
      }
  }
  for (char c : result)
    if (!std::isdigit(static_cast<unsigned char>(c)) && c != '+' && c != '-' && c != '.' &&
        c != 'e' && c != 'E')
      return {};
  return result;
}

std::vector<Token> tokenize(std::string_view source, bool character_values,
                            bool* unterminated_tail) {
  // Return punctuation/record tokens independently of READ or namelist state.
  // Quoted strings decode doubled delimiters and retain embedded punctuation.
  std::vector<Token> result;
  if (unterminated_tail)
    *unterminated_tail = false;
  std::size_t i = 0, line = 1;
  bool repetition_star   = false;
  const auto punctuation = [character_values](char c) {
    switch (c) {
    case ',':
      return Kind::comma;
    case '/':
      return Kind::slash;
    case '*':
      return character_values ? Kind::value : Kind::star;
    case '=':
      return character_values ? Kind::value : Kind::equal;
    case '(':
      return character_values ? Kind::value : Kind::left_paren;
    case ')':
      return character_values ? Kind::value : Kind::right_paren;
    case ':':
      return character_values ? Kind::value : Kind::colon;
    case '%':
      return character_values ? Kind::value : Kind::percent;
    case '&':
    case '$':
      return character_values ? Kind::value : Kind::group;
    default:
      return Kind::value;
    }
  };
  while (i < source.size()) {
    const char c = source[i];
    if (c == '\n') {
      result.push_back({Kind::end_record, {}, false, line++, i, i + 1});
      ++i;
    } else if (std::isspace(static_cast<unsigned char>(c))) {
      ++i;
    } else if (c == '!') {
      while (i < source.size() && source[i] != '\n')
        ++i;
    } else if (c == '\'' || c == '"') {
      const char quote      = c;
      const auto start_line = line;
      const auto start      = i;
      std::string value;
      ++i;
      bool closed = false;
      while (i < source.size()) {
        const char next = source[i++];
        if (next == quote) {
          if (i < source.size() && source[i] == quote) {
            value += quote;
            ++i;
          } else {
            closed = true;
            break;
          }
        } else if (next == '\n') {
          ++line;
        } else if (next != '\r') {
          value += next;
        }
      }
      if (!closed && unterminated_tail) {
        *unterminated_tail = true;
        return result;
      }
      if (!closed)
        throw std::invalid_argument("Unterminated Fortran character input at line " +
                                    std::to_string(start_line));
      result.push_back({Kind::value, std::move(value), true, start_line, start, i});
    } else if (character_values && repetition_star && c == '*') {
      result.push_back({Kind::star, "*", false, line, i, i + 1});
      repetition_star = false;
      ++i;
    } else if (const Kind kind = punctuation(c); kind != Kind::value) {
      result.push_back({kind, std::string(1, c), false, line, i, i + 1});
      ++i;
    } else {
      const auto begin = i;
      while (i < source.size() && !std::isspace(static_cast<unsigned char>(source[i])) &&
             source[i] != '!' && source[i] != '\'' && source[i] != '"' &&
             punctuation(source[i]) == Kind::value) {
        if (character_values && source[i] == '*' && i > begin &&
            std::all_of(source.begin() + begin, source.begin() + i,
                        [](unsigned char digit) { return std::isdigit(digit); })) {
          repetition_star = true;
          break;
        }
        ++i;
      }
      result.push_back(
          {Kind::value, std::string(source.substr(begin, i - begin)), false, line, begin, i});
    }
  }
  return result;
}
} // namespace clubb_kokkos::fortran_io
