// Apply Fortran namelist assignments to registered caller-owned C++ fields.
#include "namelist.hpp"

#include <cctype>
#include <fstream>
#include <iterator>
#include <numeric>

namespace clubb_kokkos::namelist {
using fortran_io::Kind;
using fortran_io::Token;

std::size_t Group::Binding::size() const {
  // A scalar is one value; arrays follow Native first-index-fastest storage.
  std::size_t count = 1;
  for (int extent : extents)
    count *= static_cast<std::size_t>(extent);
  return count;
}

bool Group::assign_text(std::string& value, const std::string& text, int width, int first,
                        int last) {
  // Implement fixed-width character assignment and optional Native substrings.
  std::string result;
  if (first != 0) {
    if (width <= 0 || first < 1 || first > width + 1 || last < 0 || last > width)
      return false;
    if (last < first)
      return true; // A zero-length Native character substring leaves its owner unchanged.
    result = value;
    result.resize(static_cast<std::size_t>(width), ' ');
    std::string piece = text;
    piece.resize(static_cast<std::size_t>(last - first + 1), ' ');
    result.replace(static_cast<std::size_t>(first - 1), piece.size(), piece);
  } else {
    result = text;
    if (width > 0 && result.size() > static_cast<std::size_t>(width))
      result.resize(static_cast<std::size_t>(width));
  }
  const auto end = result.find_last_not_of(' ');
  result.resize(end == std::string::npos ? 0 : end + 1);
  value = std::move(result);
  return true;
}

void Group::insert(const std::string& name, Binding binding) {
  // Keep declaration order for aggregate derived-type assignments.
  const auto key = fortran_io::lower(name);
  if (!bindings_.emplace(key, std::move(binding)).second)
    throw std::invalid_argument("Duplicate namelist field: " + name);
  order_.push_back(key);
}

bool Group::reference(const std::vector<Token>& tokens, std::size_t start, Reference& result) {
  // Recognize name[(indices)][%component...] followed by '=' without consuming
  // value tokens. Contiguous section bounds may be omitted.
  const auto identifier = [](const Token& token) {
    if (token.kind != Kind::value || token.quoted || token.text.empty() ||
        (!std::isalpha(static_cast<unsigned char>(token.text.front())) &&
         token.text.front() != '_'))
      return false;
    for (char c : token.text)
      if (!std::isalnum(static_cast<unsigned char>(c)) && c != '_')
        return false;
    return true;
  };
  const auto skip_records = [&tokens](std::size_t& p) {
    while (p < tokens.size() && tokens[p].kind == Kind::end_record)
      ++p;
  };
  std::size_t p = start;
  skip_records(p);
  if (p == tokens.size() || !identifier(tokens[p]))
    return false;
  result      = Reference{};
  result.name = fortran_io::lower(tokens[p++].text);
  for (;;) {
    skip_records(p);
    if (p < tokens.size() && tokens[p].kind == Kind::left_paren) {
      std::vector<Subscript> dimensions;
      ++p;
      for (;;) {
        skip_records(p);
        Subscript subscript;
        if (p < tokens.size() && tokens[p].kind == Kind::value) {
          if (!fortran_io::parse_integer(tokens[p++].text, subscript.first))
            return false;
          subscript.has_first = true;
        }
        skip_records(p);
        if (p < tokens.size() && tokens[p].kind == Kind::colon) {
          subscript.section = true;
          ++p;
          skip_records(p);
          if (p < tokens.size() && tokens[p].kind == Kind::value) {
            if (!fortran_io::parse_integer(tokens[p++].text, subscript.last))
              return false;
            subscript.has_last = true;
          }
          skip_records(p);
          if (p < tokens.size() && tokens[p].kind == Kind::colon)
            return false; // Explicit section strides are rejected by Native NVHPC namelist input.
        } else if (!subscript.has_first) {
          return false;
        }
        dimensions.push_back(subscript);
        skip_records(p);
        if (p == tokens.size())
          return false;
        if (tokens[p].kind == Kind::right_paren) {
          ++p;
          break;
        }
        if (tokens[p++].kind != Kind::comma)
          return false;
      }
      result.indices.push_back(std::move(dimensions));
    } else if (p < tokens.size() && tokens[p].kind == Kind::percent) {
      ++p;
      skip_records(p);
      if (p == tokens.size() || !identifier(tokens[p]))
        return false;
      result.name += '%' + fortran_io::lower(tokens[p++].text);
    } else {
      break;
    }
  }
  if (p == tokens.size() || tokens[p].kind != Kind::equal)
    return false;
  result.next = p + 1;
  return true;
}

Group::Binding* Group::find_binding(const std::string& name, Binding& aggregate) {
  // Scalar derived objects may be supplied as one list or as %components.
  auto found = bindings_.find(name);
  if (found != bindings_.end())
    return &found->second;
  std::vector<Binding*> members;
  std::size_t total = 0;
  for (const auto& key : order_)
    if (key.compare(0, name.size() + 1, name + '%') == 0) {
      members.push_back(&bindings_.at(key));
      total += members.back()->size();
    }
  if (members.empty())
    return nullptr;
  aggregate.extents      = {static_cast<int>(total)};
  aggregate.lower_bounds = {1};
  aggregate.assign       = [members](std::size_t index, const Token& token, int first, int last) {
    for (const auto* member : members) {
      if (index < member->size())
        return member->assign(index, token, first, last);
      index -= member->size();
    }
    return false;
  };
  return &aggregate;
}

bool Group::destinations(const Binding& binding, const Reference& target,
                         std::vector<std::size_t>& indices, int& first, int& last) {
  // An indexed element continues through remaining array storage. An explicit
  // section instead traverses only its selected elements, in Fortran order.
  std::size_t subscript_group = 0;
  first = last = 0;
  if (!binding.extents.empty() && !target.indices.empty()) {
    const auto& subscripts = target.indices[subscript_group++];
    if (subscripts.size() != binding.extents.size())
      return false;
    bool section = false;
    for (const auto& subscript : subscripts)
      section = section || subscript.section;
    if (!section) {
      std::size_t offset = 0, stride = 1;
      for (std::size_t p = 0; p < subscripts.size(); ++p) {
        const int index = subscripts[p].first - binding.lower_bounds[p];
        if (index < 0 || index >= binding.extents[p])
          return false;
        offset += static_cast<std::size_t>(index) * stride;
        stride *= binding.extents[p];
      }
      for (std::size_t i = offset; i < binding.size(); ++i)
        indices.push_back(i);
    } else {
      indices            = {0};
      std::size_t stride = 1;
      for (std::size_t p = 0; p < subscripts.size(); ++p) {
        const auto& subscript = subscripts[p];
        const int low         = binding.lower_bounds[p];
        const int high        = low + binding.extents[p] - 1;
        const int begin       = subscript.has_first ? subscript.first : low;
        const int end = subscript.section ? (subscript.has_last ? subscript.last : high) : begin;
        if (begin < low || begin > high || end < low || end > high)
          return false;
        std::vector<std::size_t> selected;
        for (int i = begin; i <= end; ++i)
          for (auto previous : indices)
            selected.push_back(previous + static_cast<std::size_t>(i - low) * stride);
        indices = std::move(selected);
        stride *= binding.extents[p];
      }
    }
  } else {
    indices.resize(binding.size());
    std::iota(indices.begin(), indices.end(), std::size_t(0));
  }
  if (subscript_group < target.indices.size()) {
    if (!binding.character || target.indices.size() != subscript_group + 1)
      return false;
    const auto& substring = target.indices[subscript_group];
    if (substring.size() != 1 || !substring[0].section)
      return false;
    first = substring[0].has_first ? substring[0].first : 1;
    last  = substring[0].has_last ? substring[0].last : binding.character_width;
    if (first < 1 || first > binding.character_width + 1 || last < 0 ||
        last > binding.character_width)
      return false;
  }
  return true;
}

void Group::finish() {
  // Restore logical View elements after success or a partially applied read.
  for (auto& field : bindings_)
    if (field.second.finish)
      field.second.finish();
}

int Group::read_file(const std::string& filename) {
  // Read the first matching group, returning zero on success, negative at EOF or
  // a missing group, and positive for malformed input, unknown names or bounds.
  error_.clear();
  // Native OPEN ignores trailing character padding in the file name.
  const auto last = filename.find_last_not_of(' ');
  const auto path = filename.substr(0, last == std::string::npos ? 0 : last + 1);
  std::ifstream file(path, std::ios::binary);
  if (!file) {
    error_ = "Cannot open namelist file: " + filename;
    return 1;
  }
  const std::string source{std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>()};
  std::size_t group_start = std::string::npos;
  for (std::size_t record = 0; record < source.size();) {
    const auto record_end = source.find('\n', record);
    const auto end        = record_end == std::string::npos ? source.size() : record_end;
    std::size_t first     = record;
    while (first < end && std::isspace(static_cast<unsigned char>(source[first])))
      ++first;
    if (first < end && (source[first] == '&' || source[first] == '$')) {
      std::size_t name_end = first + 1;
      while (name_end < end && (std::isalnum(static_cast<unsigned char>(source[name_end])) ||
                                source[name_end] == '_'))
        ++name_end;
      if (fortran_io::lower(std::string_view(source).substr(first + 1, name_end - first - 1)) ==
          name_) {
        group_start = first;
        break;
      }
    }
    record = end == source.size() ? end : end + 1;
  }
  if (group_start == std::string::npos) {
    error_ = "Namelist &" + name_ + " not found in " + filename;
    return -1;
  }
  // Native scans unrelated records without interpreting their quotes/comments.
  // Tokenize only after finding a matching first-nonblank, adjacent &/$name.
  bool unfinished = false;
  auto tokens =
      fortran_io::tokenize(std::string_view(source).substr(group_start), false, &unfinished);
  const auto line_offset =
      static_cast<std::size_t>(std::count(source.begin(), source.begin() + group_start, '\n'));
  for (auto& token : tokens)
    token.line += line_offset;
  std::size_t p = 0;
  p += 2;
  const auto failure = [&](const std::string& message) {
    const auto line = p < tokens.size() ? tokens[p].line : 0;
    error_          = filename + ':' + std::to_string(line) + " (&" + name_ + "): " + message;
    finish();
    return 1;
  };
  const auto ending = [&tokens](std::size_t index) {
    return index < tokens.size() &&
           (tokens[index].kind == Kind::slash ||
            (tokens[index].kind == Kind::group && index + 1 < tokens.size() &&
             tokens[index + 1].kind == Kind::value && !tokens[index + 1].quoted &&
             fortran_io::lower(tokens[index + 1].text).compare(0, 3, "end") == 0));
  };
  for (;;) {
    while (p < tokens.size() &&
           (tokens[p].kind == Kind::end_record || tokens[p].kind == Kind::comma))
      ++p;
    if (ending(p)) {
      finish();
      return 0;
    }
    if (p == tokens.size()) {
      error_ = "Unexpected end of namelist &" + name_ + " in " + filename;
      finish();
      return unfinished ? 1 : -1;
    }
    Reference target;
    if (!reference(tokens, p, target))
      return failure("Expected a namelist assignment");
    Binding aggregate;
    Binding* binding = find_binding(target.name, aggregate);
    if (!binding)
      return failure("Unknown field " + target.name);
    std::vector<std::size_t> indices;
    int first, last;
    if (!destinations(*binding, target, indices, first, last))
      return failure("Invalid indices for " + target.name);
    p                       = target.next;
    std::size_t destination = 0;
    bool expecting_value    = true;
    for (;;) {
      while (p < tokens.size() && tokens[p].kind == Kind::end_record)
        ++p;
      Reference following;
      if (ending(p) || p == tokens.size() || reference(tokens, p, following))
        break;
      if (tokens[p].kind == Kind::comma) {
        ++p;
        if (expecting_value) {
          if (destination >= indices.size())
            return failure("Too many values for " + target.name);
          ++destination;
        }
        expecting_value = true;
        continue;
      }
      std::size_t repeat = 1;
      bool repeated_null = false;
      if (tokens[p].kind == Kind::value && !tokens[p].quoted && !tokens[p].text.empty() &&
          std::all_of(tokens[p].text.begin(), tokens[p].text.end(),
                      [](unsigned char c) { return std::isdigit(c); }) &&
          p + 1 < tokens.size() && tokens[p + 1].kind == Kind::star &&
          tokens[p].end == tokens[p + 1].begin) {
        int count;
        if (tokens[p].quoted || !fortran_io::parse_integer(tokens[p].text, count) || count <= 0)
          return failure("Invalid repetition count for " + target.name);
        repeat                 = static_cast<std::size_t>(count);
        const auto value_begin = tokens[p + 1].end;
        p += 2;
        repeated_null = p == tokens.size() || tokens[p].begin != value_begin;
      }
      if (repeat > indices.size() - std::min(destination, indices.size()))
        return failure("Too many values for " + target.name);
      const bool null = repeated_null || p == tokens.size() || ending(p) ||
                        tokens[p].kind == Kind::comma || reference(tokens, p, following);
      Token value;
      std::size_t value_end = p;
      if (!null) {
        value     = tokens[p];
        value_end = p + 1;
      }
      if (!null && value.kind != Kind::value)
        return failure("Invalid value for " + target.name);
      for (std::size_t i = 0; i < repeat; ++i) {
        if (destination >= indices.size())
          return failure("Too many values for " + target.name);
        if (!null && !binding->assign(indices[destination], value, first, last))
          return failure("Invalid value for " + target.name);
        ++destination;
      }
      if (!null)
        p = value_end;
      expecting_value = false;
    }
  }
}
} // namespace clubb_kokkos::namelist
