---
name: clubb-porting-cpp-kokkos
description: Map C++ ownership, indexing and Kokkos kernels, layouts, memory spaces and workspace for a source-faithful portable port.
---

# C++ and Kokkos mapping

Use when C++/Kokkos is source or target with the
[porting strategy](../../SKILL.md). For plain C++, apply ownership/type/index
guidance without imposing Kokkos. These are general considerations, not commands
or support boundaries for an existing model implementation.

## Types, indices and storage

Use native zero-based storage with explicit ranges. Keep counts, physical
levels, option values, sentinels and opaque/serialized IDs distinct. Convert
source storage conventions at shared boundaries instead of keeping one-based
values and subtracting one at nearly every access. Signed sentinel handling
and unsigned extents need particular care.

Choose field ownership before mirroring derived types. C++ object copying,
container copying and `Kokkos::View` copying have different semantics: View
handles alias their allocation. An independent snapshot needs explicit storage
and copying. Avoid unmanaged references outliving their owner, and define
move/copy/reset behavior for case objects containing arrays or handles.

The last index is contiguous in ordinary row-major C++ arrays; a Kokkos View
uses its declared layout. Choose dimension order, LayoutLeft/LayoutRight/Stride
and the thread-varying index together. Match neighboring kernels and bridge
costs rather than blindly copying or reversing source axes. A padded or strided
View is not generally a flat `data()[0:product(extents)]` array: traverse logical
indices or prove density/strides before flat access. Subviews may be noncontiguous.

## Start with simple valid kernels

Map independent work to `RangePolicy` with explicit flattened-index decoding
or `MDRangePolicy`. Flattening can improve some loops, but is not a universal
winner. A column-local recurrence stays sequential inside independent outer
work. Never parallelize a dependency merely to increase launch size.

Compile and run target-device kernels early. Use device-callable functions,
small scalar/View captures and portable numeric traits; avoid host I/O or
host-only calls in device work. Compiler restrictions around templates, extended
lambdas and reference-bound constants can surface only in the device build.
Verify the installed compiler/Kokkos version rather than applying a workaround
or launch policy from another machine blindly.

Prefer a clear initial kernel policy. Introduce teams, vector hierarchy or team
scratch when data reuse/synchronization warrants it, not as an automatic porting
layer. Count barriers, registers, scratch footprint and available work. One-case
latency and large-batch throughput are different objectives.

## Keep data and workspace alive for the right duration

Keep persistent fields in their execution memory across adjacent calls. Stage
only at explicit host/I/O/interop boundaries; wrapping a host pointer in an
unmanaged View does not move it to the device. Never dereference device-only
storage from host code. Count bridge copies separately from native kernel work.

Allocate reusable scratch at a case/context or scoped-workspace owner, with
explicit capacity, shape and lifetime. Avoid per-timestep/routine allocation
when a stable pool suffices. Preserve simultaneous temporary lifetimes, nesting,
alignment, separate snapshots and reset; do not introduce global/TLS caches or
opaque reused variables to save allocations. Team scratch and persistent global
device workspace solve different problems.

Skip initialization only after proving every subsequently read element is
written on all relevant branches. Fence at actual lifetime or host-observation
boundaries, not after every launch. Release owning Views before runtime
finalization; repeated cases can reuse the runtime while resetting case state.

## Keep numerical dependencies visible

Use typed literals and inspect promotions, integer versus real powers,
transcendental/complex libraries, signed zero, FMA and reduction order. `-O0`
does not guarantee cross-compiler or device BFB. Parallel iteration order,
atomic addition and RNG-pool stream assignment can differ from sequential
source behavior. Keep deterministic mapping/order where the contract requires
it, and diagnose differences before proposing a tradeoff.

Apply simple fusion/removal of unnecessary copies when dependencies/options
remain clear. Measure deeper changes with optimized device builds, explicit
measurement synchronization and the intended workload. A device debug build
is unsuitable for performance conclusions. Allocation counts, requested copy
payloads and measured DMA bytes are different quantities.

Read [the performance catalog](references/performance.md) for policy/layout,
transfer or profiling choices and [studies/source examples](references/papers-and-code.md)
when examining tuning designs. Verify APIs against the installed
[Kokkos documentation](https://kokkos.org/kokkos-core-wiki/).
