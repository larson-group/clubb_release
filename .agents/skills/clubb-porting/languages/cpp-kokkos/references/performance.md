# Kokkos GPU mapping and performance catalog

These are design checks and experiments, not performance claims about a model.
Before applying one, verify the installed Kokkos version, compiler/backend,
device architecture and supported algorithm options. Follow the parent language
skill and the port's agreed source-fidelity and numerical contract. Kokkos' [programming guide](https://kokkos.org/kokkos-core-wiki/)
and [CMake integration guide](https://kokkos.org/kokkos-core-wiki/get-started/integrating-kokkos-into-your-cmake-project.html)
are the API references; link `Kokkos::kokkos` and use a compiler configured for
the selected backend. Avoid hard-coding a particular Kokkos release's defaults.

## Expose useful parallel work

| Source pattern | First mapping to investigate | Main check |
| --- | --- | --- |
| Independent column and level calculation | `RangePolicy` over a flattened index or `MDRangePolicy` over column and level | Adjacent lanes access adjacent addresses; launch has enough work. |
| Column-local chain through vertical levels | Independent columns or samples outside, serial level loop inside | Do not parallelize a recurrence as if levels were independent. |
| Several levelwise stages sharing intermediate data | `TeamPolicy` per column/tile, with `TeamThreadRange` or `ThreadVectorRange` inside | Team scratch/reuse benefit exceeds barriers, register use, and lost occupancy. |
| Independent loops with simple nested shape | `MDRangePolicy` or a flattened range | Compare both when a multidimensional reduction is slow. |
| Many tiny sequential calls | Batch independent work; consider fusion only after profiling | Preserve dependencies, BFB-sensitive operation order, and debuggability. |

Kokkos recommends multidimensional ranges for tightly nested independent loops
and team hierarchy when threads need shared intermediate work or synchronization
([hierarchical parallelism](https://kokkos.org/kokkos-core-wiki/ProgrammingGuide/HierarchicalParallelism.html),
[multidimensional parallelism](https://kokkos.org/kokkos-core-wiki/ProgrammingGuide/Multi-Dimensional-Parallelism.html)).
`Kokkos::AUTO` or `team_size_recommended` is a starting point, not a final team
size; scratch and vector length affect launch limits
([TeamPolicy](https://kokkos.org/kokkos-core-wiki/API/core/policies/TeamPolicy.html)).
Treat small single-case latency and large-batch throughput as different
benchmarks. Never increase batch size silently in a numerical comparison.

## Make memory movement visible

- For a Fortran array, `LayoutLeft` retains first-index-contiguous storage;
  `view(column, level)` therefore makes adjacent columns contiguous, matching
  Fortran `(column, level)` in this illustrative mapping. Ordinary C++ row-major
  arrays make the *last* index contiguous; `Kokkos::View` uses its declared
  layout. `view(level, column)` would make columns contiguous with `LayoutRight`.
  With one column, the level stride is the View's leading dimension, which
  may include padding; inspect `stride()` before assuming adjacent levels are
  contiguous. For GPU coalescing, arrange the *thread-varying* index
  contiguously and verify real access patterns. Test
  `LayoutRight`, a transposed/flattened field, or a different iteration order
  when the first extent is tiny or threads vary over another index. Layout and
  iteration order must be chosen together ([View](https://kokkos.org/kokkos-core-wiki/ProgrammingGuide/View.html),
  [MDRange](https://kokkos.org/kokkos-core-wiki/ProgrammingGuide/Multi-Dimensional-Parallelism.html)).

  ```cpp
  // Fortran a(column, level): the column index varies fastest in memory.
  Kokkos::View<double**, Kokkos::LayoutLeft> a("a", ncol, nlev);
  // Inside a compatible Kokkos kernel: a(icol, k).
  ```
- Keep hot loops explicit and cheap: avoid whole-array/slice expressions that
  allocate temporaries or copy data per call. A `Kokkos::subview` is an alias,
  but its stride can be noncontiguous; check the resulting access before using
  it inside a hot kernel
  ([Subviews](https://kokkos.org/kokkos-core-wiki/ProgrammingGuide/Subviews.html)).
- Give hot fields clear device ownership. Use `create_mirror_view_and_copy` and
  `deep_copy` at explicit boundaries; never dereference a device-only View from
  host code. A `View` copy is shallow, not a data copy. Avoid per-level or
  per-routine host copies while a sequence of ported kernels can stay on device
  ([View](https://kokkos.org/kokkos-core-wiki/ProgrammingGuide/View.html),
  [mirror API](https://kokkos.org/kokkos-core-wiki/API/core/view/create_mirror.html)).
- Reuse workspace across timesteps and cases when shape and lifetime allow.
  A normal View allocates and initializes its data; `WithoutInitializing`
  removes initialization only when every element is written before any read.
  Measure persistent footprint, scratch footprint, and temporary traffic
  separately ([View](https://kokkos.org/kokkos-core-wiki/ProgrammingGuide/View.html)).
- Initialize Kokkos once for the standalone/worker process and release owning
  Views before finalization. Repeated worker requests should reuse the runtime
  rather than paying setup for every case; verify case state itself is reset
  ([initialization](https://kokkos.org/kokkos-core-wiki/ProgrammingGuide/Initialization.html)).
- For C++ records containing managed `View` members, do not use an outer
  device `View` without an explicit ownership/lifetime design: a View of objects
  containing Views needs special handling. Consider flat/structure-of-arrays
  hot storage and convert at explicit interop boundaries; confirm fieldwise
  semantics before replacing a record representation
  ([View of Views](https://kokkos.org/kokkos-core-wiki/ProgrammingGuide/View.html#can-i-make-a-view-of-views)).
- Externally owned Fortran buffers can be wrapped as *host* unmanaged Views
  only while the native allocation lives and with its true layout. Wrapping
  does not place host memory on the GPU. Fieldwise converters may copy through
  host buffers; record that cost separately from target-native kernel execution
  ([Kokkos interoperability](https://kokkos.org/kokkos-core-wiki/ProgrammingGuide/Interoperability.html),
  [Fortran interop use case](https://kokkos.org/kokkos-core-wiki/usecases/Kokkos-Fortran-Interoperability.html)).

## Use GPU resources deliberately

- Stage data in team scratch only when multiple threads reuse it. Count loads
  saved, team barriers, scratch bytes, team size, and active teams. Scratch
  and large private arrays can reduce occupancy or spill registers; a larger
  tile can be slower ([hierarchical parallelism](https://kokkos.org/kokkos-core-wiki/ProgrammingGuide/HierarchicalParallelism.html),
  [TERRA-NG study](https://arxiv.org/pdf/2609.21681)).
- Measure whether a kernel is limited by launch latency, global-memory
  bandwidth, cache, arithmetic throughput, branch divergence, synchronization,
  or occupancy before changing it. Check bytes moved and arithmetic intensity
  with a roofline or memory-workload profile, then test the smallest change
  addressing the observed limit
  ([Nsight Compute guide](https://docs.nvidia.com/nsight-compute/ProfilingGuide/),
  [ROCm memory analysis](https://rocm.docs.amd.com/projects/rocprofiler-compute/en/develop/how-to/membw_analysis.html)).
- Fuse adjacent stages only when their dependency order is clear and saved
  launches or intermediate writes matter in the profile. Recheck register
  count, spills, occupancy, and maintainability after fusion. For a column
  physics routine, keeping a serial inner recurrence in one launch may be
  appropriate even while columns run in parallel
  ([SCREAM study](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024MS004314),
  [NICAM preprint](https://arxiv.org/html/2610.41313v1)).
- Source programs often separate short top/bottom boundary loops from the interior.
  Consider combining nearby stages that traverse the same fields, with small
  explicit boundary branches. Keep every supported option and boundary value
  visible. Do not force a large kernel fusion that requires broad inlining or
  extra live temporaries; measure both designs.
- Use `KOKKOS_LAMBDA`/value capture for outer device kernels, and keep device
  code free of host I/O, ordinary host allocation, and host-only function calls.
  Prefer explicit Views and small scalar state over capturing a large owning
  object. Label kernels so a profile can be tied to the source routine
  ([parallel dispatch](https://kokkos.org/kokkos-core-wiki/ProgrammingGuide/ParallelDispatch.html),
  [parallel_for API](https://kokkos.org/kokkos-core-wiki/API/core/parallel-dispatch/parallel_for.html)).


- Compile new template kernels with the GPU compiler while CPU model checks run.
  CUDA extended lambdas require named enclosing template parameters, including
  non-type SFINAE parameters. NVCC also rejects a variable first captured only
  inside `if constexpr`; move rank dispatch outside the launch or give that
  value an unconditional use before the branch, such as `(void)view;`.
  These constraints affect compilation, not arithmetic
  ([CUDA lambda restrictions](https://docs.nvidia.com/cuda/cuda-programming-guide/05-appendices/cpp-language-support.html)).
- Use the installed device-compatible numeric traits, such as
  `Kokkos::finite_max_v<T>`, instead of calling host-only `numeric_limits`
  functions inside kernels. Preserve the actual Fortran argument kind:
  `HUGE(0.0)` need not have the same kind as a core-real array
  ([numeric traits](https://kokkos.org/kokkos-core-wiki/API/core/numerics/numeric-traits.html)).
  Avoid passing global mathematical constants/traits by reference in device
  code. For example, a `Kokkos::complex` constructor takes scalar references;
  copy the value into a local scalar first
  ([known constant-address issues](https://kokkos.org/kokkos-core-wiki/known-issues.html#mathematical-constants-and-numeric-traits)).

## Preserve numerical meaning

- `parallel_for` has no prescribed iteration order. A Kokkos reduction can be
  repeatable on the *same hardware and runtime settings* without matching
  the source's sequential sum or a different backend. For an initial BFB port,
  keep a source-order accumulation when required; separately measure the
  parallel alternative before requesting a numerical tradeoff
  ([Kokkos machine model](https://kokkos.org/kokkos-core-wiki/ProgrammingGuide/Machine-Model.html)).
- Atomic accumulation resolves write conflicts, but concurrent addition order
  is a separate BFB concern. Prefer independent outputs or a fixed-order
  follow-up accumulation when exactness requires it. `ScatterView` trades
  atomics for replicated storage and a later combination; it is an option to
  benchmark, not a guarantee of source order
  ([atomics](https://kokkos.org/kokkos-core-wiki/ProgrammingGuide/Atomic-Operations.html),
  [ScatterView](https://kokkos.org/kokkos-core-wiki/API/containers/ScatterView.html)).
- Stochastic algorithms need a defined mapping from case/column/sample to random
  state. Kokkos documents that acquiring a CUDA RNG pool state uses atomics
  and can assign streams nondeterministically. A fixed pool seed alone does
  not give each sample the same stream. Preserve the agreed source sequence
  for BFB, or design and verify a deterministic indexed stream scheme before
  changing it ([Kokkos RNG](https://kokkos.org/kokkos-core-wiki/API/algorithms/Random-Number.html)).
- Keep numerical build details with the comparison: source and target compilers,
  Kokkos/backend version, architecture, real kind, flags, FMA contraction,
  fast-math/flush-to-zero, transcendental functions, and reduction mapping.
  `-O0` does not enforce cross-language or cross-device BFB; FMA and math flags
  can alter rounding ([NVIDIA floating-point guide](https://docs.nvidia.com/cuda/floating-point/index.html),
  [CUDA math guidance](https://docs.nvidia.com/cuda/cuda-programming-guide/05-appendices/mathematical-functions.html)).

## Profile and report the whole path

1. Label kernels with the source routine/stage name and wrap larger phases in
   a profiling region. Kokkos Tools can measure launch counts, kernel time,
   allocations, and region time without embedding a timer in every kernel
   ([Kokkos Tools](https://kokkos.org/projects/kokkos-tools/),
   [ScopedRegion](https://kokkos.org/kokkos-core-wiki/API/core/profiling/scoped_region.html),
   [simple kernel timer](https://github.com/kokkos/kokkos-tools/wiki/SimpleKernelTimer)).
2. Warm up the same path, repeat runs, and fence the execution space at a
   *measurement boundary*. Kokkos dispatch may be asynchronous; a host timer
   around an unfenced launch records submission time, not completed work.
   Avoid fencing between every intended back-to-back kernel just for timing
   ([fence API](https://kokkos.org/kokkos-core-wiki/API/core/parallel-dispatch/fence.html)).
3. Report case, step count, columns, levels, settings, compiler/backend/device,
   debug or optimized build, repetitions, and whether timing includes startup,
   interop copies, diagnostics, and I/O. Report both total case time and the kernels that
   dominate it. Reuse the project's comparison and reporting
   owners for numerical evidence; do not infer GPU speed from a CPU stand-in.
4. On the target device, inspect a representative hot kernel with Nsight
   Compute or ROCm profiling where available. Compare at least the intended
   GPU backends before making a device-specific tuning rule. If a single-column
   run is launch-bound, retain that result beside batched throughput; a faster
   large batch does not make a one-column run faster.

More aggressive facilities such as [Kokkos Graphs](https://kokkos.org/kokkos-core-wiki/ProgrammingGuide/Graph.html)
may reduce repeated launch overhead, but require a stable launch sequence and
separate correctness/performance evidence. They are later experiments, not a
default porting layer.
