// Reuse case-owned storage for temporary column-first Kokkos Views.
//
// ScratchWorkspace owns stable chunks in the default execution memory space and
// HostSpace. scope opens a nested reservation lifetime; scope.view reserves a
// mutable arithmetic View and initializes it to zero unless WithoutInitializing
// is explicit. Workspace.view reserves an escaping result in the currently
// active caller scope, before the helper opens any local scope. clear releases
// the chunks after all scopes and asynchronous users have finished. The two
// reserved_*_bytes queries report retained chunk capacity, including alignment.
//
// Raw-pointer Views preserve Real1D/2D/3D/4D types and ordinary LayoutLeft strides,
// but do not carry an allocation tracker or individual label. Kernel labels and
// chunk labels remain available to profiling. The Views do not own reservations:
// their uses must finish, or be ordered before reuse on the same default queue,
// before the allocating scope ends. Host scratch used by asynchronous transfers
// must remain reserved until the transfers finish. Scope release never fences.
#pragma once

#include <Kokkos_Core.hpp>
#include <array>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace clubb_kokkos {

class ScratchWorkspace {
private:
  static constexpr std::size_t alignment = 128;
  struct Cursor {
    std::size_t chunk{}, offset{};
  };

  template <class ExecutionSpace> class Pool {
  public:
    using storage_type =
        Kokkos::View<unsigned char*, Kokkos::LayoutLeft, typename ExecutionSpace::memory_space>;

    Pool(std::size_t initial_bytes, std::string_view kind)
        : initial_bytes_(initial_bytes), next_bytes_(initial_bytes), kind_(kind) {
      // The default instance matches the existing default-policy physics kernels.
    }

    unsigned char* reserve(std::size_t bytes, std::string_view label) {
      // Grow by adding chunks; moving View handles never relocates live storage.
      if (bytes == 0)
        return nullptr;
      if (bytes > std::numeric_limits<std::size_t>::max() - (alignment - 1))
        throw std::length_error("Scratch reservation size overflows alignment");
      while (true) {
        if (cursor_.chunk == chunks_.size()) {
          const std::size_t required  = bytes + alignment - 1;
          const std::size_t requested = next_bytes_ < required ? required : next_bytes_;
          const std::size_t tail      = (alignment - requested % alignment) % alignment;
          if (requested > std::numeric_limits<std::size_t>::max() - tail)
            throw std::length_error("Scratch chunk alignment overflows");
          const std::size_t capacity = requested + tail;
          if (capacity > std::numeric_limits<std::size_t>::max() - capacity_bytes_)
            throw std::length_error("Scratch workspace capacity overflows");
          storage_type storage(
              Kokkos::view_alloc(execution_, Kokkos::WithoutInitializing,
                                 "scratch/" + std::string(kind_) + "/" + std::string(label)),
              capacity);
          chunks_.push_back(std::move(storage));
          capacity_bytes_ += capacity;
          next_bytes_ = capacity > std::numeric_limits<std::size_t>::max() / 2
                            ? std::numeric_limits<std::size_t>::max()
                            : 2 * capacity;
        }
        auto& chunk               = chunks_[cursor_.chunk];
        const auto base_remainder = reinterpret_cast<std::uintptr_t>(chunk.data()) % alignment;
        const auto remainder      = (base_remainder + cursor_.offset % alignment) % alignment;
        const auto padding        = remainder == 0 ? 0 : alignment - remainder;
        const auto available      = chunk.extent(0) - cursor_.offset;
        if (padding <= available && bytes <= available - padding) {
          auto* pointer = chunk.data() + cursor_.offset + padding;
          cursor_.offset += padding + bytes;
          return pointer;
        }
        ++cursor_.chunk;
        cursor_.offset = 0;
      }
    }

    Cursor cursor() const {
      // Capture the live prefix before opening a nested routine scope.
      return cursor_;
    }

    void restore(Cursor cursor) noexcept {
      // Release reservations only; queued work and chunks remain untouched.
      cursor_ = cursor;
    }

    void clear() {
      // The caller has completed all users of these chunks at a case boundary.
      chunks_.clear();
      cursor_         = {};
      capacity_bytes_ = 0;
      next_bytes_     = initial_bytes_;
    }

    std::size_t reserved_bytes() const {
      // Return allocated chunk payload, including unused alignment/tail bytes.
      return capacity_bytes_;
    }

    const ExecutionSpace& execution_space() const {
      // Supply the same queue for zero fills as for the default physics policy.
      return execution_;
    }

  private:
    ExecutionSpace execution_{};
    std::vector<storage_type> chunks_;
    Cursor cursor_{};
    std::size_t initial_bytes_{}, next_bytes_{}, capacity_bytes_{};
    std::string_view kind_;
  };

public:
  class Scope {
  public:
    Scope(const Scope&)            = delete;
    Scope& operator=(const Scope&) = delete;
    Scope(Scope&&)                 = delete;
    Scope& operator=(Scope&&)      = delete;

    ~Scope() noexcept {
      // Scopes are a host-side stack; destruction restores both pool cursors.
      if (workspace_.depth_ != depth_)
        std::terminate();
      workspace_.device_.restore(device_);
      workspace_.host_.restore(host_);
      --workspace_.depth_;
    }

    template <class View, class... Dims> View view(std::string_view label, Dims... dims) {
      // Reserve zero-initialized storage for this scope's nonescaping temporaries.
      require_current();
      return workspace_.template allocate_view<View>(true, label, dims...);
    }

    template <class View, class... Dims>
    View view(decltype(Kokkos::WithoutInitializing), std::string_view label, Dims... dims) {
      // Preserve explicit uninitialized allocation when every later read is defined.
      require_current();
      return workspace_.template allocate_view<View>(false, label, dims...);
    }

  private:
    friend class ScratchWorkspace;
    explicit Scope(ScratchWorkspace& workspace)
        : workspace_(workspace), device_(workspace.device_.cursor()),
          host_(workspace.host_.cursor()), depth_(++workspace.depth_) {
      // Bookmark the parent reservations before any local allocation is issued.
    }

    void require_current() const {
      // A parent cannot reserve storage while a child owns the bump cursor.
      if (workspace_.depth_ != depth_)
        throw std::logic_error("Scratch allocation requires the innermost scope");
    }

    ScratchWorkspace& workspace_;
    Cursor device_, host_;
    std::size_t depth_;
  };

  explicit ScratchWorkspace(std::size_t device_chunk_bytes = 1024 * 1024,
                            std::size_t host_chunk_bytes   = 64 * 1024)
      : device_(device_chunk_bytes, "device"), host_(host_chunk_bytes, "host") {
    // Initial capacities are growth hints; chunks are allocated on first demand.
  }
  ScratchWorkspace(const ScratchWorkspace&)            = delete;
  ScratchWorkspace& operator=(const ScratchWorkspace&) = delete;
  ScratchWorkspace(ScratchWorkspace&&)                 = delete;
  ScratchWorkspace& operator=(ScratchWorkspace&&)      = delete;

  ~ScratchWorkspace() {
    // The owning case must outlive every scope and complete outstanding users.
    if (depth_ != 0)
      std::terminate();
  }

  Scope scope() {
    // Open a routine lifetime without allocation, synchronization or zeroing.
    return Scope(*this);
  }

  template <class View, class... Dims> View view(std::string_view label, Dims... dims) {
    // Reserve a returned result in the caller's scope before opening a local one.
    require_scope();
    return allocate_view<View>(true, label, dims...);
  }

  template <class View, class... Dims>
  View view(decltype(Kokkos::WithoutInitializing), std::string_view label, Dims... dims) {
    // Reserve an uninitialized caller-owned result with the same lifetime rule.
    require_scope();
    return allocate_view<View>(false, label, dims...);
  }

  void clear() {
    // Release retained chunks after an explicit case-level completion boundary.
    if (depth_ != 0)
      throw std::logic_error("Cannot clear scratch with active scopes");
    device_.clear();
    host_.clear();
  }

  std::size_t reserved_device_bytes() const {
    // HostSpace default-execution Views use the host pool on a CPU backend.
    return device_.reserved_bytes();
  }

  std::size_t reserved_host_bytes() const {
    // Include all HostSpace temporary storage retained for this case.
    return host_.reserved_bytes();
  }

private:
  template <class Dim> static std::size_t checked_extent(Dim value) {
    // Reject negative extents before converting Native signed dimensions.
    static_assert(std::is_integral_v<Dim>, "Scratch dimensions must be integral");
    if constexpr (std::is_signed_v<Dim>)
      if (value < 0)
        throw std::invalid_argument("Scratch extent is negative");
    if (static_cast<std::uintmax_t>(value) > std::numeric_limits<std::size_t>::max())
      throw std::length_error("Scratch extent is too large");
    return static_cast<std::size_t>(value);
  }

  template <class View, class... Dims>
  View allocate_view(bool initialize, std::string_view label, Dims... dims) {
    // Use the public raw-View sizing contract, including its trailing byte rounding.
    static_assert(View::rank >= 1 && View::rank <= 4);
    static_assert(View::rank == View::rank_dynamic);
    static_assert(sizeof...(Dims) == View::rank);
    static_assert(std::is_same_v<typename View::array_layout, Kokkos::LayoutLeft>);
    static_assert(std::is_same_v<typename View::data_type, typename View::non_const_data_type>);
    using Value = typename View::non_const_value_type;
    static_assert(std::is_arithmetic_v<Value> && alignof(Value) <= alignment);
    const std::array<std::size_t, sizeof...(Dims)> extents{checked_extent(dims)...};
    bool empty = false;
    for (const auto extent : extents)
      empty = empty || extent == 0;
    if (!empty) {
      std::size_t count = 1;
      for (const auto extent : extents) {
        if (count > std::numeric_limits<std::size_t>::max() / extent)
          throw std::length_error("Scratch extent product overflows");
        count *= extent;
      }
      if (count > (std::numeric_limits<std::size_t>::max() - 7) / sizeof(Value))
        throw std::length_error("Scratch View byte size overflows");
    }
    const auto bytes = empty ? 0 : View::required_allocation_size(checked_extent(dims)...);
    auto& pool       = pool_for<View>();
    auto* pointer    = reinterpret_cast<typename View::pointer_type>(pool.reserve(bytes, label));
    View result(pointer, checked_extent(dims)...);
    if (initialize && !empty)
      Kokkos::deep_copy(pool.execution_space(), result, Value{});
    return result;
  }

  template <class View> auto& pool_for() {
    // Keep host and execution storage explicit, without a memory-pointer registry.
    using Memory = typename View::memory_space;
    if constexpr (std::is_same_v<Memory, Kokkos::HostSpace>) {
      return host_;
    } else {
      static_assert(std::is_same_v<Memory, typename Kokkos::DefaultExecutionSpace::memory_space>,
                    "Scratch supports HostSpace and the default execution memory space");
      return device_;
    }
  }

  void require_scope() const {
    // Workspace-only output allocation must have a caller lifetime to reserve into.
    if (depth_ == 0)
      throw std::logic_error("Scratch output requires an active caller scope");
  }

  Pool<Kokkos::DefaultExecutionSpace> device_;
  Pool<Kokkos::DefaultHostExecutionSpace> host_;
  std::size_t depth_{};
};

} // namespace clubb_kokkos
