// Read numeric case files with their declared number of entries per record.
//
// file_read_1d reads successive full rows and a possible final partial row.
// file_read_2d reads those rows separately for each first-dimension level.
// Outputs retain Native index order and caller-owned storage. Host buffers cross
// to execution-space or strided Views through the shared array owner. Each call
// owns its file; the required Native logical-unit arguments need no registry.
#pragma once

#include "array_types.hpp"
#include <span>
#include <string>
#include <variant>
#include <vector>

namespace clubb_kokkos {

void file_read_1d(int file_unit, const std::string& path_and_filename,                     // In
                  int num_datapts, int entries_per_line,                                   // In
                  Kokkos::View<real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace> variable);  // Out
void file_read_2d(int device, const std::string& file_path,                                // In
                  int file_dimension1, int file_dimension2, int file_per_line,             // In
                  Kokkos::View<real_t**, Kokkos::LayoutLeft, Kokkos::HostSpace> variable); // Out

template <class Variable>
void file_read_1d(int file_unit, const std::string& path_and_filename, // In
                  int num_datapts, int entries_per_line,               // In
                  Variable variable) {                                 // Out
  // Retain caller storage, including values left unchanged by a null item.
  Kokkos::View<real_t*, Kokkos::LayoutLeft, Kokkos::HostSpace> host = to_host_space(variable);
  file_read_1d(file_unit, path_and_filename, num_datapts, entries_per_line, host);
  from_execution_space(variable, host);
}
template <class Variable>
void file_read_2d(int device, const std::string& file_path,                    // In
                  int file_dimension1, int file_dimension2, int file_per_line, // In
                  Variable variable) {                                         // Out
  // The Native routine initializes its complete output before reading records.
  Kokkos::View<real_t**, Kokkos::LayoutLeft, Kokkos::HostSpace> host = to_host_space(variable);
  file_read_2d(device, file_path, file_dimension1, file_dimension2, file_per_line, host);
  from_execution_space(variable, host);
}

namespace file_functions_detail {
// A language I/O adapter shared by the numeric, grid and SAM-like readers.
// It retains physical-record positioning and typed READ destinations; lexical
// rules and numeric conversion remain in the shared fortran_io owner.
struct read_item {
  std::variant<real_t*, int*, std::string*> destination;
  std::size_t character_length = std::string::npos;
  read_item(real_t& value) : destination(&value) {}
  read_item(int& value) : destination(&value) {}
  read_item(std::string& value, std::size_t length = std::string::npos)
      : destination(&value), character_length(length) {}
};
class formatted_file {
public:
  explicit formatted_file(const std::string& filename);
  explicit formatted_file(std::vector<std::string> records);
  int read_record(std::string& record);
  int read(std::span<read_item> items);
  void backspace();

private:
  std::vector<std::string> records_;
  std::size_t position_ = 0;
  bool after_endfile_   = false;
};
} // namespace file_functions_detail
} // namespace clubb_kokkos
