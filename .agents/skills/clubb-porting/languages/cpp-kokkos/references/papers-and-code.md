# Studies and source examples to revisit during tuning

These are case studies, not expected model speedups. Record each workload,
precision, GPU, and baseline before transferring an idea. Keep the agreed source
reference and numerical contract even when a paper uses statistical tests.

## Atmospheric models

**[SCREAM / EAMxx (Donahue et al., 2024)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024MS004314).**
This is the closest atmospheric-physics example. It uses column/team hierarchy,
reusable workspace, reduced memory movement, and CPU vectorization; its
row-normalization example shows a team reduction followed by a team-local
barrier. Process implementations have narrow driver-facing interfaces.
Read the [archived source](https://doi.org/10.5281/zenodo.10724808) and
[performance data](https://github.com/E3SM-Project/perf-data). Its global
3-km grid has far more horizontal work than a one-column target workload, so its
scaling does not predict standalone latency.

**[NICAM Fortran→Kokkos (Ren et al., Oct 2026 preprint)](https://arxiv.org/html/2610.41313v1).**
Listing 1 shows a concrete `LayoutLeft` / `MDRangePolicy` translation with a
contiguous horizontal index, zero-based exclusive bounds, local View handles,
and explicit transfer. Fusion reduced a measured turbulence/surface phase from
58.8 to 28.5 ms; another phase slowed. This is a very recent unreviewed
preprint. Its 100-member ensemble supports statistical consistency, not source-level
BFB; its horizontal grid and OpenACC baseline also differ from the target model.

## Comparative kernel studies

**[ICS 2025 programming-model comparison (Davis et al.)](https://pssg.cs.umd.edu/assets/papers/2025-06-portability-ics.pdf).**
In CloverLeaf, a 2D Kokkos `calc_dt` reduction showed barrier stalls.
Flattening it to 1D improved performance on the tested NVIDIA and AMD systems.
Profile reduction topology before assuming a structured policy wins. Such
regrouping may change bits; the five proxy apps do not establish a general
backend ranking or a target-model speedup.

**[TERRA-NG Kokkos study (Böhm et al., Sep 2026 preprint)](https://arxiv.org/pdf/2609.21681).**
Explicit `LayoutRight` beat GPU-default `LayoutLeft` for radial-first thread
access and a nearly unused first axis. Team scratch staged shared data. A
per-team accumulator cut atomics and raised MI300A throughput 27%, but hurt
MI250X; a private accumulator spilled registers. Browse the staged
[source variants](https://github.com/mantleconvection/TERRA-NG/tree/main/src/terra/fe/wedge/operators/shell/performance_history),
especially [v04 scratch](https://github.com/mantleconvection/TERRA-NG/blob/main/src/terra/fe/wedge/operators/shell/performance_history/epsilon_divdiv_kerngen_v04_shmem_coords.hpp)
and [v09 scatter](https://github.com/mantleconvection/TERRA-NG/blob/main/src/terra/fe/wedge/operators/shell/performance_history/epsilon_divdiv_kerngen_v09_separate_scatter.hpp).
Finite-element assembly differs from column physics. Its atomic strategy
is a performance experiment, not a BFB recipe.

**[neXtSIM-DG GPU study (Jendersie et al., 2025)](https://gmd.copernicus.org/articles/18/3017/2025/gmd-18-3017-2025.html).**
On the tested A100, Kokkos approached CUDA kernel time for large problems but
was about 50% slower on the small problem. The paper excludes host/device
transfer from that comparison. Use its [code snapshot](https://doi.org/10.5281/zenodo.13711171)
to examine launch overhead and problem-size crossover, while measuring the target's
full runtime as well. Sea-ice finite elements are not a full model forecast.

If an optimization conflicts with BFB, keep the original and tuned algorithms
separately measurable until the numerical cause and performance gain are clear.
