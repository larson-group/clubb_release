// Read typed, caller-owned Fortran namelist groups without Native module state.
//
// Group::bind registers initialized scalars or Kokkos Views; read_file applies
// assignments in file order, preserving omitted/null values and fixed character
// widths. Array indices are Native one-based unless explicit bounds are supplied.
// bind_array accepts explicit Native lower bounds and stages strided/device
// Views once. Keep scalar owners and array allocations alive until the read ends.
// read_file returns zero on success, negative for EOF/missing groups and positive
// for malformed input; error provides the corresponding diagnostic. Construct
// and bind a fresh Group for each read so defaults come from the current owner.
//
// Supported syntax includes derived components, repeated/null values, contiguous
// sections and character substrings. Character values require quotes; explicit
// section strides are rejected, matching the installed NVHPC namelist reader.
#pragma once

#include "CLUBB_core/array_types.hpp"
#include "fortran_io.hpp"

#include <algorithm>
#include <functional>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace clubb_kokkos::namelist {
class Group {
public:
  explicit Group(std::string name) : name_(fortran_io::lower(name)) {
    // Identify one case-insensitive group without opening or retaining a file.
  }

  template <class T> void bind(const std::string& name, T&& value, int character_width = 0) {
    // Bind a scalar reference or View while retaining its initialized values.
    // Parsed character scalars are truncated at Native width and right-trimmed,
    // matching the former bridge's unpack_text conversion into std::string.
    using Value = std::decay_t<T>;
    if constexpr (Kokkos::is_view<Value>::value) {
      bind_array(name, value);
    } else {
      static_assert(std::is_lvalue_reference_v<T>, "Namelist scalars need caller-owned storage");
      Binding binding;
      if constexpr (std::is_same_v<Value, std::string>) {
        binding.character       = true;
        binding.character_width = character_width;
        if (character_width > 0)
          assign_text(value, value, character_width, 0, 0);
      }
      binding.assign = [&value, character_width](std::size_t index, const fortran_io::Token& token,
                                                 int first, int last) {
        if (index != 0)
          return false;
        if constexpr (std::is_same_v<Value, std::string>) {
          return token.quoted && assign_text(value, token.text, character_width, first, last);
        } else {
          return first == 0 && !token.quoted && fortran_io::parse_value(token.text, value);
        }
      };
      insert(name, std::move(binding));
    }
  }

  template <class View>
  void bind_array(const std::string& name, View values, std::vector<int> lower_bounds = {}) {
    // Stage each bound array once. Logical indexing preserves Fortran order
    // even when the caller supplies a strided or device-resident View.
    using Element = typename View::non_const_value_type;
    static_assert(View::rank >= 1 && View::rank <= 4, "Unsupported namelist array rank");
    auto host = to_host_space(values);
    Binding binding;
    constexpr bool character      = std::is_same_v<Element, char>;
    constexpr int first_dimension = character ? 1 : 0;
    binding.character             = character;
    if constexpr (character)
      binding.character_width = static_cast<int>(host.extent(0));
    for (int p = first_dimension; p < View::rank; ++p)
      binding.extents.push_back(static_cast<int>(host.extent(p)));
    if (lower_bounds.empty())
      lower_bounds.assign(binding.extents.size(), 1);
    if (lower_bounds.size() != binding.extents.size())
      throw std::invalid_argument("Namelist lower bounds do not match array rank: " + name);
    binding.lower_bounds = std::move(lower_bounds);
    binding.assign       = [host](std::size_t index, const fortran_io::Token& token, int first,
                                  int last) {
      if constexpr (character) {
        if (!token.quoted)
          return false;
        const auto width = host.extent(0);
        std::string value(width, ' ');
        const auto get = [host, index](std::size_t c) -> char& {
          if constexpr (View::rank == 1)
            return host(c);
          else if constexpr (View::rank == 2)
            return host(c, index);
          else if constexpr (View::rank == 3)
            return host(c, index % host.extent(1), index / host.extent(1));
          else
            return host(c, index % host.extent(1), (index / host.extent(1)) % host.extent(2),
                        index / (host.extent(1) * host.extent(2)));
        };
        for (std::size_t c = 0; c < width; ++c)
          value[c] = get(c);
        if (!assign_text(value, token.text, static_cast<int>(width), first, last))
          return false;
        value.resize(width, ' ');
        for (std::size_t c = 0; c < width; ++c)
          get(c) = value[c];
        return true;
      } else {
        if (first != 0 || token.quoted)
          return false;
        const auto assign = [&token](auto& destination) {
          return fortran_io::parse_value(token.text, destination);
        };
        if constexpr (View::rank == 1)
          return assign(host(index));
        else if constexpr (View::rank == 2)
          return assign(host(index % host.extent(0), index / host.extent(0)));
        else if constexpr (View::rank == 3)
          return assign(host(index % host.extent(0), (index / host.extent(0)) % host.extent(1),
                             index / (host.extent(0) * host.extent(1))));
        else
          return assign(host(index % host.extent(0), (index / host.extent(0)) % host.extent(1),
                             (index / (host.extent(0) * host.extent(1))) % host.extent(2),
                             index / (host.extent(0) * host.extent(1) * host.extent(2))));
      }
    };
    binding.finish = [values, host] { from_execution_space(values, host); };
    insert(name, std::move(binding));
  }

  int read_file(const std::string& filename);
  const std::string& error() const { return error_; }

private:
  struct Binding {
    std::vector<int> extents;
    std::vector<int> lower_bounds;
    bool character      = false;
    int character_width = 0;
    std::function<bool(std::size_t, const fortran_io::Token&, int, int)> assign;
    std::function<void()> finish;
    std::size_t size() const;
  };
  struct Subscript {
    bool section   = false;
    bool has_first = false;
    bool has_last  = false;
    int first = 0, last = 0;
  };
  struct Reference {
    std::string name;
    std::vector<std::vector<Subscript>> indices;
    std::size_t next = 0;
  };
  static bool assign_text(std::string& value, const std::string& text, int width, int first,
                          int last);
  void insert(const std::string& name, Binding binding);
  static bool reference(const std::vector<fortran_io::Token>& tokens, std::size_t start,
                        Reference& result);
  Binding* find_binding(const std::string& name, Binding& aggregate);
  bool destinations(const Binding& binding, const Reference& reference,
                    std::vector<std::size_t>& indices, int& first, int& last);
  void finish();
  std::string name_, error_;
  std::map<std::string, Binding> bindings_;
  std::vector<std::string> order_;
};
} // namespace clubb_kokkos::namelist
