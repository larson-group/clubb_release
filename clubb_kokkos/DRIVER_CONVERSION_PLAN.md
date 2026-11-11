# Standalone driver port

Port the supported standalone scientific work and services into C++. The
original driver milestone retained Fortran statistics; the later stats milestone
below replaces that final service boundary.
The authoritative routines are in [the Fortran source](../src/); the existing
[driver](src/clubb_driver.cpp) supplies the call order and owns each case's state.
The supported options remain those agreed for the JAX port, including the
already supported fixed dycore grid. BUGSrad is excluded. The KK/Morrison numerical
engines are included: their C++ results must match, while source resemblance,
formatting cleanup and performance work inside those engine subfolders are waived
for this milestone. Their driver-facing interfaces retain the normal standards.

## Source groups

| Native directory | Top-level files | Work |
| --- | ---: | --- |
| [Input_fields](../src/Input_fields/) | 11 | Sounding and input interpretation, hydrostatic setup, readers and restart data. |
| [Benchmark_cases](../src/Benchmark_cases/) | 33 | Time-dependent input, forcing dispatch, case tendencies and surface fluxes. |
| [Radiation](../src/Radiation/) | 9 | Supported radiation, initialization and extended-atmosphere state; retain the BUGSrad exclusion. |
| [SILHS](../src/SILHS/) | 11 | Sampling, transformations and case-owned random streams. |
| [Microphys](../src/Microphys/) | 15 | Supported initialization, coordination, mathematical helpers and engine interfaces. |

Also port the reachable KK/Morrison engine routines and the scientific
interpretation in [readiopdata_module](../src/readiopdata_module.F90) and its core.
The latter preserves the existing NetCDF input route; those routines contain
thermodynamic calculations and cannot be classified as raw file-I/O services.

The counts describe the initial inventory, not completed ports. Account for
every routine in each group's source outline, including excluded-only code and
already ported code such as the hydrometeor PDF wrapper. Inspect the remaining
driver call tree as well as these directories so services outside them are not
silently missed.

## Integration and acceptance

1. Read each complete Native source and record routines, arguments, branch and
   call order, state, statistics and errors before drafting its port.
2. Keep independent drafts under [drafts](drafts/). Root integration owns the
   live driver, CMake and temporary [f2k](src/f2k/) code. Persistent module state
   becomes case-owned C++ state; scratch and converters remain call-owned.
3. Integrate compatible, source-audited groups. Replace the corresponding
   temporary whole-routine calls and remove obsolete bridge entries. Preserve
   source names, argument groups, comments and directory layout. Keep scientific
   calculations visible in C++; justify retained file-I/O and other service calls
   at their actual boundaries.
4. Use the maintained [Kokkos/Fortran comparison](auto_llm_generated_tests/integration_tests/run_kokkos_vs_fortran_cases.py)
   with short runs and multi-column statistics. Require strict saved-file BFB
   between Serial CPU Kokkos and CPU Fortran. Diagnose a failing group before
   attributing it to an individual file. Record forced arithmetic choices in
   [BFB_ADAPTATIONS.md](BFB_ADAPTATIONS.md).
5. Audit source fidelity and clean up each group separately from its numerical
   check. Finish with the full supported case matrix at native timesteps and
   duration using standard statistics, plus selected option/edge checks.

GPU performance work is deferred. An occasional CUDA build/runtime checkpoint
records the integrated group where saved GPU output changes; GPU numerical-cause
tracing and exact GPU agreement are not acceptance requirements for this phase.
The standalone entry remains the existing `-exe` route.

## Starting evidence

Commit `ac37a3898` completes the supported CLUBB core. Its frozen CPU comparison
passes all 19 cases at native duration/timestep with standard statistics, plus
32 short case/option comparisons. Existing lifecycle/bridge checks and CUDA
build/runtime evidence are recorded in [the core plan](CORE_CONVERSION_PLAN.md).
That evidence establishes the starting point; it does not verify these new
driver ports.


## Integrated source

All supported groups are hooked into the live driver. The source reviews cover
routine inventories, Native leading arguments and call order, supported branches,
comments, stats and fatal returns separately from the numerical checks.

| Source group | Current C++ owner and review |
| --- | --- |
| Input_fields | [Input_fields](src/Input_fields/): all 11 files and 46 routines; [source outline](drafts/driver_input_fields/PACKAGE_OUTLINE.json). |
| Benchmark_cases | [Benchmark_cases](src/Benchmark_cases/): all 33 files and 74 routines, including time-dependent forcing and fixed-dycore remapping; [source outline](drafts/driver_benchmark_cases/SOURCE_OUTLINE.json). |
| Radiation | [Radiation](src/Radiation/): eight non-BUGS files plus common initialization; 22 reviewed contracts, including the explicitly excluded SILHS-radiation entry; [source outline](drafts/driver_radiation_silhs/Radiation/SOURCE_OUTLINE.json). |
| SILHS | [SILHS](src/SILHS/): all 11 files and 116 specifics, including the case-owned MT95 stream; [source outline](drafts/driver_radiation_silhs/SILHS/SOURCE_OUTLINE.json). |
| Microphys | [Microphys](src/Microphys/): supported coordinators/interfaces and the actual KK/Morrison engines; [top-level outline](drafts/driver_microphysics/source_outline.json), [KK closure](drafts/driver_microphysics/KK_engine_closure.json) and [Morrison audit](drafts/driver_microphysics/Morrison_engine_audit.json). Engine calculations currently use host C++ with View staging; GPU performance is deferred. |
| IOP data | [Wrapper](src/readiopdata_module.cpp) and [core](src/readiopdata_core/readiopdata_core_module.cpp): all 14 routines, preserving Native float working precision and data order; [source outline](drafts/driver_input_fields/readiopdata_SOURCE_OUTLINE.json). |

Native module settings now belong to each `ClubbDriver` instance. Appended state
arguments preserve each routine's original leading formals. The obsolete driver
settings snapshot/setter interface is removed. Permanent generated C++ namelist
bindings assign case-owned values, including the six sponge-setting records,
without reconstructing Native model settings or retaining bridge copies.

The only production Fortran boundary is the optional PDF matrix
factorization/diagnostics ABI in `KOKKOS_BFB=ON` mode. Default PDF matrix
operations are custom C++; LAPACK solves and condition-estimate requests return
fatal `err_info`. Statistics, output, debug levels and all supported scientific
routines now execute in C++. The generic procedure manifest and stats bridges
are retired. Native type-conversion helpers and numerical oracles remain
test-only; the local CMake project still configures Native reference targets.

The known Native IOP undefined branches are retained explicitly at their C++
boundaries: repeated forcing reads, add-surface forcing and incomplete input
shapes need separate Native fixes before their values can serve as references.
The short matching IOP case does not establish those branches.

## Driver validation

The Serial CPU port and CPU Fortran match all 19 comparison cases at native
duration and timestep with standard statistics and zero tolerances. Each backend
completes 17,298 main steps, with 38 normal exits and byte-identical NetCDF files
for all 19 pairs. Independent file hashes cover 50,261,711,616 bytes across the
two output trees. The configuration is NVHPC 26.9, real8, Debug/O0 and four C8
parameter columns.

These results come from frozen parallel lanes. The six final KK/SILHS cases use
the rebuilt live source after BFB-031; the earlier non-microphysics cases retain
their full results after a separate source/control compatibility review. The
solar correction is covered by six new full cases. The long GABLS2 simulation
also matches. Later Morrison-only arithmetic changes do not enter those paths.
A separate full Morrison LBA run matches all 360 native 60-second steps with
standard statistics, sampled microphysics disabled and Arctic nucleation off.

Additional strict CPU checks pass:

- All 19 cases for five steps with multi-column statistics, and five-step
  all-stats COBRA and DYCOMS RF02 DO/DS checks.
- The supplied BOMEX IOP dataset and a BOMEX restart from a common Native
  three-step all-stats file, each advanced for five steps.
- ATEX for 190 steps at 30 seconds, with and without the fixed dycore grid,
  reaching the active forcing/remapping branch after 5,400 seconds.
- Cloud sedimentation without another microphysics scheme in DYCOMS RF01,
  plus short Morrison LBA/MPACE-A and stochastic SILHS RICO checks.
- A reset/reinitialization sequence that matches a normal two-step BOMEX run,
  including independent parameter snapshots and a stats-suppressed advance.

All 14 existing CTest checks and 31 focused pytests pass. The small padded-column
regression fails with the previous hole detector and passes with logical indexing;
it also passes on the actual GPU. The four canonical producers regenerate their
outputs without changes. Source inventories, signatures, branches, comments and
formatting were reviewed separately from numerical validation. The engine style
waiver applies only inside the KK/Morrison subdirectories.

The complete CUDA build succeeds. The four-case GPU checkpoint finishes five
30-second steps per case, with normal exits and finite unmasked output. At the
outer-driver integration boundary, GPU files change for COBRA, DYCOMS RF02 DO
and JUN25; BOMEX and all Native reference files stay identical to the preceding
core checkpoint. The subsequent BFB-028–031/stride-fix checkpoint introduces no
additional saved GPU changes in that scope. GPU/CPU differences retain their
strict comparator failure status; no numerical-cause or performance investigation
was performed. Current five-step GPU checks also finish normally with finite unmasked output
for cloudy Morrison MPACE-A, stochastic SILHS RICO, BOMEX IOP input and cloud
sedimentation without another microphysics scheme. The Morrison saved file is
byte-identical to the CPU reference; the other GPU differences remain recorded
diagnostics. Their source hashes and settings are stored separately.

Local commands, source/binary/control hashes, independent audits and preserved
failed/passing results are under `build/validation/remaining_driver_cpu/`,
`build/parallel_cpu_lanes/` and `build/parallel_gpu_lanes/`. In particular,
`FULL19_INDEPENDENT_RESULT_AUDIT.json` records the selected full-case evidence,
and `FROZEN_CASE_SOURCE_COMPATIBILITY_REVIEW.json` explains its compatibility.
These are manual local checks; no Jenkins job or upstream runner selector was
added. Unexercised branches and the agreed exclusions remain explicit in the
source reviews. Debug results do not establish BFB for other compilers or
optimized builds.

## Custom SILHS matrix operations

The default matrix operations now use the JAX-style column-wise Cholesky
algorithm and a Jacobi eigenvalue diagnostic. `KOKKOS_BFB=ON`, local to the
Kokkos CMake project, retains the earlier Fortran/f2k matrix path. The public
interfaces, lower-triangle ownership, scaling policy and ten-attempt retry
coordinator remain unchanged. The BFB coordinator's executable tokens match
the preceding implementation; its four PDF matrix bridge symbols are absent
from default CPU/CUDA archives and present only in the BFB archive.

Both Serial modes pass all 19 five-step comparisons with multi-column stats,
four parameter columns, debug level 2 and zero comparison tolerances. Both also
match a 20-step stochastic LBA_KK_SILHS comparison with standard stats, eight
samples, importance sampling, random start and calculated w correlations. These
short checks do not establish equality for full-duration custom runs. An extra
RICO stress run with those overrides stops in both Native and C++ on the
cloud-normal sample consistency diagnostic; that failed control is retained.

The portable factorization runs successfully in the actual CUDA execution
space in the focused reconstruction/scaling/retry/failure/eigenvalue checks.
The custom-mode GPU RICO_SILHS checkpoint completes five steps per backend,
with finite unmasked saved data and recorded GPU/CPU differences. No numerical
cause or performance investigation was performed, and no claim is made about
which GPU bits changed relative to a preceding GPU build. Host staging remains
in the production factorization coordinator for later performance work.

Build modes and results are recorded under `build/validation/custom_cholesky/`.
The comparison runner records requested and effective `KOKKOS_BFB` values;
configuration checks cover a missing Fortran compiler or either required target.
`compile.py`, upstream configuration and Jenkins definitions remain unchanged.

## Configuration and error output

The two configuration printers and per-case error-header formatting are C++.
The printers retain all 64 CLUBB and 13 SILHS ordered label/field mappings,
including Native's um/vm labels. Caller-owned output streams replace logical
units; obsolete integer-unit overloads and both print bridges are removed.
Printing preserves caller stream flags and width, and keeps existing default
C++ setup-text output. Vendor-specific list-directed whitespace is not emulated.

The error-context lifecycle keeps the Native Aw/Iw/F6.2 records, fixed widths,
line breaks and overflow fields. Only headers change during a formatting update;
coordinates, metadata, error codes and column allocations retain their owners.
Fixed-field conversion uses call-local host storage; its diagnostic-only rounding
choice is recorded in [the BFB log](BFB_ADAPTATIONS.md#diagnostic-header-formatting).
Native's unused module-global process/thread prefix has no executable consumer,
so its dead initialization call is removed without introducing C++ global state.
The four obsolete service selections are removed from the procedure generator.
That configuration/error-output replacement left stats, timers, debug levels
and PDF matrix routines unchanged. The later timer replacement is recorded below.

For this replacement, all 19 five-step Serial/CPU comparisons match saved files
byte-for-byte with four parameter columns, multi-column stats, debug level 2 and
`KOKKOS_BFB=ON`. Full BOMEX also matches all 360 native steps with standard stats.
The focused diagnostic check matches Native global/column headers at field and
rounding edges, preserves error codes and allocation, and checks printer stream
state. The existing NaN error-return check and 43 focused Python checks pass;
all four generators remain idempotent. Source review separately covers routine
order, formal arguments, all printer rows and the unused legacy-prefix exclusion.
Commands, source review and results are under `build/validation/config_error_cpp/`.
The CUDA diagnostic check also passes on the RTX 3080. Short RICO_SILHS and
BOMEX IOP GPU runs finish normally with finite data; all GPU saved files and
Native controls are byte-identical to the preceding input-reader checkpoint.
GPU/CPU comparator differences remain recorded, without numerical-cause tracing.

## C++ input readers

Grid heights, formatted tables, generic profile namelists, driver/microphysics/
radiation namelists and raw NetCDF input now use permanent C++ implementations.
The configuration inventory covers 11 groups and 455 fields, plus the generic
profile's three fields; Native names,
argument order, array lower bounds, kinds and character widths were reviewed
separately from numerical tests. A shared call-local parser preserves defaults,
nulls, repeated values and Native array order. The NetCDF-C service converts
Native dimension and slab conventions at its boundary. The obsolete input f2k
sources and three procedure selections have been removed.

The [reader documentation](README.md#input-readers) records supported syntax,
IOSTAT categories and explicit handling of Native undefined empty-table and
omitted-profile-count inputs. Existing IOP undefined branches listed above are
unchanged; replacing raw I/O does not establish their scientific validity.
Validation for this later input replacement is recorded separately under
`build/validation/input_io_cpp/`; historical full-driver evidence above describes
its original binaries.

The input replacement passes all 19 five-step Serial/CPU comparisons with four
parameter columns, multi-column stats, debug level 2, `KOKKOS_BFB=ON` and zero
tolerances. IOP input, a common Native restart fixture and cloudy MPACE-A each
match for five steps. A full 360-step native-duration/timestep BOMEX run with
standard stats also matches byte-for-byte. Three tiny CTest readers and 43 focused
Python checks pass; all four producers regenerate without output changes.
These checks do not establish a new full-duration 19-case result for the input port.

The CUDA standalone builds, and the three reader checks pass on the RTX 3080.
Five-step RICO_SILHS and BOMEX IOP runs finish normally with finite unmasked GPU
data. Their strict GPU/CPU comparisons still report differences. RICO's GPU
saved file and its Native control are byte-identical to the preceding custom-matrix
checkpoint; no new GPU change is observed in that scope. No GPU numerical-cause
or performance investigation was performed.

## Timers

All 16 routines in [code_timer_module](src/CLUBB_core/code_timer_module.cpp)
now execute in C++, with the Native public API, helper order, nesting/error
checks, reset/disable behavior and output fields. The default backend retains
process CPU time; `USE_GPTL=ON` uses the library's C API directly with the same
options. Private per-thread bookkeeping is the agreed exception to case-owned
model state. No GPU fences or scientific calculations are added.

All 79 live calls retain their names and locations. The five timer entries are
removed from f2k; the later debug port below leaves seven generic stats services. Native
stats have no calls to the timing module or its instrumented parent routines,
so C++ owns the active GPTL lifecycle. The timer replacement did not change PDF
matrices or debug levels.

The unchanged [timing script](../utilities/time_clubb.py) completes CPU and GPTL
sweeps over process counts 1/2 and batch sizes 1/2, with one warmup and two
repetitions per workload. All 36 peer requests succeed; every request records
three core advances and one complete advance region, confirming reset across
persistent worker requests. An additional three-request GPU profile succeeds.

Both Serial modes match all 19 five-step model comparisons byte-for-byte with
multi-column stats, four parameter columns, debug level 2 and `KOKKOS_BFB=ON`.
Full BOMEX also matches all 360 native steps with standard stats in the default
CPU-time build. Timing values themselves are not BFB assertions. The small
timer check covers lifecycle, nesting, growth and output contracts in both modes;
the CUDA toolchain also compiles/runs the host-only GPTL check.

The full CUDA standalone builds. Short RICO_SILHS and BOMEX IOP runs finish
normally with finite data and saved GPU/Native files unchanged from the preceding
configuration/error-output checkpoint. Strict GPU/CPU differences remain
diagnostics; no cause or performance investigation is included. Evidence and
source reviews are under `build/validation/timers_cpp/`. This is not a new
full-duration 19-case result, and no upstream selector or timing-parser change is
included. See [the README](README.md#timers-and-timing-profiles) for commands.

## Debug threshold

Both public accessors in [error_code](src/CLUBB_core/error_code.cpp) now execute
in C++. The setter retains arbitrary signed levels and the check retains
`level <= clubb_debug_level`. Private per-host-thread diagnostic metadata starts
at zero, following the Native threadprivate state. The legacy unused
process/thread header remains excluded. No model fields or arithmetic changed.

All 156 existing debug calls in 42 source/test files use the permanent API with
the same arguments and order. Both generic debug selections and their generated
bindings were removed; the seven remaining selections at that checkpoint were stats calls. The
former stats boundary synchronized its Native debug threshold
when borrowing a mutable handle, after rejecting snapshots. That checkpoint preserved Native output-grid diagnostics. The later statistics
port removes this synchronization and owns its buffers, sampling, budgets and
NetCDF output in C++; PDF matrix routines remain unchanged.

The small [debug check](auto_llm_generated_tests/debug_level.cpp) compares signed
and extreme levels with Native and verifies C++ thread isolation, local setter
ownership, stats refresh and snapshot rejection. It passes in the Serial
CPU-time and GPTL builds. Existing stats-bridge, fatal-error and matrix-wrapper
checks also pass. The matrix fixture now explicitly sets both independent debug
owners to the same level before comparing solvers.

The default Serial CPU-time build matches all 19 five-step cases byte-for-byte
with multi-column stats, four parameter columns, debug level 2 and
`KOKKOS_BFB=ON`. Full BOMEX matches all 360 native steps with standard stats.
These are not new full-duration results for all 19 cases. Source and runtime
evidence are under `build/validation/debug_cpp/`.

The CUDA standalone builds and its host-only debug check passes. Five-step
RICO_SILHS and BOMEX IOP runs finish normally with finite saved data. Their GPU
files and Native controls match the preceding timer checkpoint byte-for-byte;
the strict GPU/CPU differences remain recorded. No GPU numerical-cause or
performance investigation is included.

## Checkpoint before the stats port

The source audit found no substantive unaccepted defect in the integrated
driver, microphysics or build/input/service code. Minor cleanup corrected stale
bridge-ownership comments and continuation indentation, with unchanged
executable tokens. Statistics and its Native diagnostic helpers remain behind
the existing boundary; the accepted option exclusions and GPU differences remain.

The current NVHPC 26.9 real8 Debug/O0 Serial build with `KOKKOS_BFB=ON` passes
all 21 CTests and 55 selected Kokkos/timing-parser Python checks. All 19 model
cases match saved files byte-for-byte through 20 native timesteps with four C8
columns, multi-column stats and debug level 2. A fresh Morrison LBA run also
matches all 360 native timesteps with standard stats, sampled microphysics
disabled and Arctic nucleation off.

The refreshed default `KOKKOS_BFB=OFF` CPU build passes debug/Cholesky checks,
five-step RICO_SILHS and a 20-step stochastic KK/SILHS LBA run. Both model checks
are byte-identical to Native in this scope. The latter uses standard stats,
eight samples, importance sampling, random starting levels and
`l_calc_w_corr=.true.`. This does not promise BFB for every custom-factorization input.
Audit records, effective inputs, commands and hashes are under
`build/validation/pre_stats_checkpoint/`. No new full-duration 19-case matrix
or optimized GPU performance claim is made by this checkpoint.

## C++ statistics

[stats_netcdf](src/CLUBB_core/stats_netcdf.cpp) implements all 40 Native routines,
including the eight generic update/budget bodies in its
[header](src/CLUBB_core/stats_netcdf.hpp). The five Native stats types have
case-owned C++ equivalents. All three
[diagnostic accumulation routines](src/CLUBB_core/stats_clubb_utilities.cpp)
retain their Native signatures, branches and ordered field updates. Source and
cleanup reviews are separate from runtime evidence; the review records live in
[drafts/stats_netcdf](drafts/stats_netcdf/) and
[drafts/stats_clubb_utilities](drafts/stats_clubb_utilities/).

All 61 production caller files switched together. The stats and generic
procedure bridges and their generator are removed. Native type converters and
numerical oracles remain test-only. Production symbol/link checks in both CPU
matrix modes find no Native statistics, scientific libraries or netCDF-Fortran;
only the four optional BFB PDF matrix bindings remain. The local CMake project
still configures Fortran reference targets. Root CMake, `compile.py` and the
upstream runner selectors are unchanged.

Accumulation uses column-first host Views, double statistics/time precision and
Native per-point counts and operation order. The permanent
[NetCDF-C boundary](src/netcdf.cpp) preserves classic-file definitions,
metadata kinds, fill values and write order. Snapshots deep-copy data and cannot
mutate or close the live file. Two execution adaptations are explicit: the
registry-overflow warning is once per context rather than once per process,
and NetCDF calls assume the current sequential host execution contract. No
arithmetic adaptation or tolerance change was needed for this milestone.

The real8 Debug/O0 NVHPC 26.9/Serial CPU gate uses four C8 columns, debug level 2,
`KOKKOS_BFB=ON` and byte-identical saved NetCDF. Both the 19-case five-step
multicolumn-stats matrix and the 19-case twenty-step standard-stats matrix pass.
The final combined standard-stats matrix also passes all 19 cases at their
native durations and timesteps: 17,298 model steps per backend, every saved
NetCDF file byte-identical. Later caller indentation and ownership-comment
cleanup has executable-token identity with the tested binaries.
Additional exact comparisons cover two runtime batches, two passive/three eddy
scalars, fixed-dycore ATEX through the forcing onset, sampled Morrison with
nonzero cloudy/rain/process diagnostics, and driver reinitialization,
suppression, reset, snapshots and successive windows. Surface precipitation in
the short Morrison cases remains zero; the synthetic edge check covers nonzero
rank-one sub-timestep means.

All 21 CPU CTests and 55 focused Kokkos/parser pytests pass. The GPTL build also
passes stats/timer checks and a two-step standard-stats BOMEX comparison; the
unchanged timing utility parses its expected advance/core counts. The stats CTest
checks two contexts, averages, budgets, batches, snapshots, SILHS sample arrays
and exact equal-mass output-grid remapping through a Native NetCDF reader and
the permanent C++ reader. Single-precision model/double-precision stats objects
and public interfaces compile on Serial and CUDA; that is not a single-precision
model/BFB claim. Three canonical generators are idempotent.

The full CUDA build and actual RTX 3080 stats/debug/Cholesky checks pass.
RICO_SILHS and BOMEX IOP complete five multicolumn-stats steps with 45 finite
floating fields. All saved GPU files and Native controls match the pre-stats
checkpoint byte-for-byte. Strict GPU/CPU differences remain accepted; no cause
or performance investigation was made. Default custom-matrix CPU RICO_SILHS
also completes normally, with the accepted SILHS diagnostic differences;
strict CPU acceptance uses BFB mode.

Commands, effective inputs, source/binary hashes, raw comparison evidence,
resolved setup failures and review proofs are local under
`build/validation/stats_cpp/`. The [README](README.md#statistics) gives the
public call sequence and maintained validation commands.


## Case data lifetime

The [driver](src/clubb_driver.cpp) follows the Native case data-region lifetime:
allocate execution fields once, refresh initial conditions on reset, keep those
fields across timesteps and bounded advances, then release them during cleanup.
The [typed field owner](src/driver_data.hpp) binds existing host allocations,
aliases dense CPU Views, and retains execution/packing storage for other layouts.
All 158 runtime field pairs are explicit in the driver; no global registry or
implicit host/device merge is used. The source and phase maps are under
[drafts/driver_data_lifetime](drafts/driver_data_lifetime/).

Existing device-capable forcing, core, PDF, transport and radiation calls share
the execution fields. Six driver loops become column-fast Kokkos loops with the
same per-cell arithmetic; two zeroing loops become scalar deep copies. Native
call/branch order and supported options stay unchanged. Initialization/reset
still run on host fields. Before a later reset, execution changes refresh the
host fields so entries retained by inactive branches keep their values. Case
cleanup and reinitialization fence before releasing device owners. The host
KK/Morrison tendency wrapper, sampled stats and routine-local scratch remain
explicit boundaries/costs for later profiling.

NVHPC 26.9 real8 Debug/O0 Serial with BFB mode passes 19 five-step multicolumn-stats
and 19 twenty-step standard-stats comparisons, plus native-duration BOMEX (360 steps)
with standard stats. All saved files are byte-identical. Additional exact checks
cover two-batch BOMEX/Cobra and scalar fields, reinit/suppression/reset/bounded
advances, active sampled cloudy Morrison (five steps), and fixed-dycore ATEX (190 steps at 30 seconds)
through forcing onset. This milestone does not repeat the full-duration 19-case
stats matrix above. The final comment/alignment rebuild has identical `.text`,
`.rodata` and `.data` sections to the CPU executable used by those comparisons.

On the RTX 3080, four five-step multicolumn-stats cases (RICO_SILHS, BOMEX IOP,
DYCOMS2_RF02_DS and LBA_KK_SILHS) and two-batch BOMEX with standard/scalar stats
match the frozen pre-lifetime GPU outputs byte-for-byte. Native controls also
match their frozen files. A separate two-step GPU driver sequence matches a
normal GPU run. Saved floating fields are finite; no new GPU answer change is
observed in this scope. Existing GPU/CPU differences remain recorded without
cause tracing.

Kokkos Tools counts compare one-step and five-step BOMEX runs with four columns,
native timesteps, `-stats none` and debug level 0. Mean requests per additional
step change as follows:

| Transfer | Payload before | Payload after | Copy calls before → after |
| --- | ---: | ---: | ---: |
| Host → device | 493,376 bytes | 5,600 bytes | 167 → 2 |
| Device → host | 208,048 bytes | 3,504 bytes | 128 → 46 |

CUDA allocation requests fall from 1,280 to 1,103 per extra step; host allocation
requests fall from 293 to 46. These include routine-local scratch and metadata.
They measure requested payload, rather than hardware DMA traffic or elapsed time.
Snapshot allocations inherit source labels, so label prefixes alone do not
prove primary-array ownership. The
[focused owner check](auto_llm_generated_tests/driver_data.cpp) directly verifies
pointer reuse and no redundant allocation/copy, on Serial and CUDA.

Commands, hashes, raw counts, comparisons and independent cleanup proofs are
local under `build/validation/driver_data_lifetime/`. No arithmetic workaround,
tolerance change, fusion, scratch reuse or upstream runner integration was added.
