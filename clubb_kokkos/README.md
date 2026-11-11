# CLUBB Kokkos

This directory holds the C++/Kokkos port of the supported CLUBB standalone
model. The [standalone](src/clubb_standalone.cpp) keeps the namelist and worker
interface of the [Fortran standalone](../src/clubb_standalone.F90). The C++
[driver](src/clubb_driver.cpp) owns the case state and implements the routines
in [clubb_driver.F90](../src/clubb_driver.F90), including initialization,
batch resets, the timestep loop, cleanup, sounding conversion and restart.
The [core advance](src/CLUBB_core/advance_clubb_core_module.cpp) implements both
routines in [advance_clubb_core_module.F90](../src/CLUBB_core/advance_clubb_core_module.F90):
the timestep's solver, PDF, diagnostics and stats sequence, and its surface-flux
helper. Temperature conversion, skewness, PDF width and grid operators also
execute in C++. The [PDF closure](src/CLUBB_core/pdf_closure_module.hpp),
[ADG1 closure](src/CLUBB_core/adg1_adg2_3d_luhar_pdf.hpp),
[PDF utilities](src/CLUBB_core/pdf_utilities.hpp) and
[saturation](src/CLUBB_core/saturation.hpp) own their Kokkos calculations.
The moment solvers, mixing length, parameter and flag initialization, and grid
construction also execute in C++. The supported core source inventory is
promoted; the executable calls the C++ core directly and does not call or link
the Fortran driver.

Supported input/sounding/restart interpretation, benchmark forcing, radiation,
SILHS sampling and microphysics now execute in C++, including the [KK](src/Microphys/KK_microphys/) and
[Morrison](src/Microphys/Morrison_microphys/) numerical engines. The [driver conversion plan](DRIVER_CONVERSION_PLAN.md)
tracks this later port and its acceptance checks. Grid/text input, namelists and
NetCDF input, configuration printing and per-case error headers now execute in
C++. Named timers and debug-level accessors also use C++ implementations. The
[statistics backend](src/CLUBB_core/stats_netcdf.hpp) accumulates on the host and
writes through NetCDF-C. PDF matrix operations use custom C++ algorithms by
default; `KOKKOS_BFB=ON` retains their Fortran/LAPACK bridge for comparisons.
LAPACK solvers and solver condition diagnostics are unsupported. The default
production executable has no Native stats or derived-type bridge dependency.
The local CMake project still configures Fortran comparison targets.

The C++ driver supports the normal standalone call and separate initialization,
reset, bounded advances, parameter and flag overrides, and state getters.
Initialize Kokkos before constructing a `ClubbDriver`. Worker requests run
sequentially in one process; mutable physics settings belong to each C++ case.

The supported feature scope follows the JAX driver: no sponge damping, input
fields mode, generalized-grid testing, adaptive-grid methods, GFDL activation,
BUGSrad or SILHS radiation. The fixed `l_add_dycore_grid` constructor is supported
by [setup_gr_dycore](src/CLUBB_core/grid_adaptation_module.hpp). Microphysics
supports `none`, `khairoutdinov_kogan` and `morrison`; radiation supports `none`,
`simplified`, `simplified_bomex` and `lba`. Restart interpretation and interpolation
and raw NetCDF reads use the C++ input services.
Supported options still require compatible combinations, as in Fortran.
Morrison retains JAX's exclusions for hail, Seifert–Beheng rain, Arctic
nucleation, cloud edge activation and a fixed droplet size distribution shape.
Its supported configuration must disable those selectors.

The [core conversion plan](CORE_CONVERSION_PLAN.md) records the completed core
phase at `ac37a3898`: all 19 native-duration/native-timestep CPU cases with
standard stats produced byte-identical saved outputs. Its CUDA build and
recorded short GPU executions also succeeded, with GPU/CPU bitwise differences
in the recorded scope. Those results do not verify the later driver port.
The recorded integrated driver passes all 19 short CPU comparisons bit-for-bit (BFB),
plus targeted all-stats checks. All 19 standard-statistics cases match at native duration and timestep
in the recorded Serial comparison lanes; their source compatibility is audited
separately. The additional full Morrison LBA comparison also matches. Detailed settings
and source compatibility are recorded in the driver conversion plan. IOP input, restart, Morrison and stochastic SILHS also have
recorded short checks. Source and numerical acceptance are recorded in that plan.
The [BFB adaptation log](BFB_ADAPTATIONS.md) records deliberate arithmetic
choices made to match the compiled Fortran result, with their reasons and evidence.
Those recorded BFB comparisons used the Fortran-backed Cholesky path; select
`KOKKOS_BFB=ON` to retain it. The default custom factorization intentionally
permits rounding differences, particularly in SILHS cases.

## Build and run

Build from the repository root with a Kokkos installation visible to CMake.
Set `Kokkos_ROOT` at configuration time if CMake does not find the installation:

```bash
cmake -S clubb_kokkos -B build/clubb_kokkos -DCMAKE_BUILD_TYPE=Debug
cmake --build build/clubb_kokkos --target clubb_kokkos_standalone
./run_scripts/run_scm.py -exe build/clubb_kokkos/clubb_standalone bomex
```

The separate [CMake project](CMakeLists.txt) builds the Fortran libraries and
the C++ executable from this checkout. It needs a C++20-capable compiler for the
driver and input targets, a compatible Fortran compiler, NetCDF C/Fortran and
`pkg-config` (with `netcdf.pc` on its search path);
it does not use `install/selected`.
Debug is the default configuration for initial correctness comparisons. Use an
explicit optimized configuration for performance measurements later.

SILHS uses a custom, JAX-style Cholesky factorization by default. For strict CPU
comparisons against the current Fortran algorithm, enable the local build option:

```bash
cmake -S clubb_kokkos -B build/clubb_kokkos -DKOKKOS_BFB=ON
cmake --build build/clubb_kokkos --target clubb_kokkos_standalone
```

Use `-DKOKKOS_BFB=OFF` to restore the custom factorization and its Jacobi
eigenvalue diagnostics. Both modes preserve lower-triangle inputs, diagonal
scaling and the Native ten-attempt retry procedure, with reusable temporary host buffers.
The portable factorization kernel also accepts device Views; GPU batching and
transfer tuning remain later work. The option is confined to this CMake project;
it is not exposed by `compile.py`.

BFB mode builds the temporary PDF matrix bindings in
[`src/f2k/lapack_bridge.F90`](src/f2k/lapack_bridge.F90). If requested without
an available Fortran/f2k/LAPACK backend, the
[configuration policy](cmake/KokkosBfb.cmake) warns and uses custom operations,
without defining `KOKKOS_BFB` for C++. Configuration output and
`kokkos_backends.json` report the effective choice separately from the request.
The default standalone does not link Fortran services; BFB mode retains only
these optional matrix bindings. This CMake project still configures Fortran
reference targets for validation.
This option does not enable LAPACK solver methods or condition-number stats,
and does not guarantee GPU/CPU bitwise agreement.

For a short cloudy Morrison run, disable the five unsupported engine selectors
in the MPACE-A case:

```bash
./run_scripts/run_scm.py -exe build/clubb_kokkos/clubb_standalone \
  -max_iters 5 -stats multi_col_stats \
  -override 'l_hail=.false.,l_seifert_beheng=.false.,l_arctic_nucl=.false.,l_cloud_edge_activation=.false.,l_fix_pgam=.false.' \
  mpace_a
```

Use a Serial-only Kokkos installation for the regular CPU BFB checks. This
checkout's local Kokkos 5.2.2 Serial installation is at
`build/kokkos_serial_nvhpc26_debug/install`. In a fresh shell with the machine's
Spack/Lmod setup sourced:

```bash
module load nvhpc/26.9 netcdf-fortran/4.6.4
export NVLOCALRC=/home/guntherhuebler/larson-spack-work/kokkos-gpu/localrc
cmake -S clubb_kokkos -B build/clubb_kokkos_nvhpc26_cpu_debug \
  -DKokkos_ROOT="$PWD/build/kokkos_serial_nvhpc26_debug/install" \
  -DCMAKE_BUILD_TYPE=Debug -DCMAKE_C_COMPILER=nvc \
  -DCMAKE_CXX_COMPILER=nvc++ -DCMAKE_Fortran_COMPILER=nvfortran \
  -DCMAKE_CXX_FLAGS_DEBUG='-O0 -g' -DCMAKE_CXX_STANDARD=20 -DKOKKOS_BFB=ON
cmake --build build/clubb_kokkos_nvhpc26_cpu_debug \
  --target clubb_kokkos_standalone clubb_standalone
```

That installation was built separately from the same Kokkos source with
`Kokkos_ENABLE_SERIAL=ON` and `Kokkos_ENABLE_CUDA=OFF`; it does not modify the
installed GPU package. A CUDA-enabled Kokkos package still selects CUDA as the
default execution space even when it also provides Serial.

For the supplied CUDA installation, load the NVHPC/Kokkos/netCDF modules in
your shell and use a fresh build directory with these compilers:

```bash
module load nvhpc/26.9 kokkos/5.2.2 netcdf-fortran/4.6.4
cmake -S clubb_kokkos -B build/clubb_kokkos_nvhpc26_debug \
  -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_COMPILER=nvcc_wrapper \
  -DCMAKE_Fortran_COMPILER=nvfortran
cmake --build build/clubb_kokkos_nvhpc26_debug
```

The machine's Spack/Lmod setup must already be sourced. Select the GPU with
`CUDA_VISIBLE_DEVICES` when several devices are available; this installation
targets the RTX 3080. The tests and executable need GPU access with this
Kokkos build.

The same runner accepts model options, for example:

```bash
./run_scripts/run_scm.py -exe build/clubb_kokkos/clubb_standalone \
  -max_iters 20 -output_dir output/kokkos_bomex bomex
```

The executable also accepts an aggregate namelist path directly. With
`-worker`, it prints `CLUBB_WORKER_READY`, reads one namelist path per line from
standard input, and prints `CLUBB_WORKER_DONE` after each successful case. The
[shared case runner](../utilities/clubb_case_runner.py) defines this protocol.
Worker requests stay in one process. The C++ driver initializes and cleans
up each case before the next request.

The [driver](src/clubb_driver.hpp) uses column-first `LayoutLeft` Views.
Initialization reads and resets host fields, then copies them to case-owned
execution allocations. Forcing, core, PDF, microphysics transport and radiation
share those allocations across timesteps and bounded advances. Batch resets
refresh the values while retaining the allocation; cleanup releases the case.
Compatible CPU Views share the original host allocation. The
[typed field owner](src/driver_data.hpp) also retains packing storage for strided
host windows, without a global field registry.
The supported core calculations, including all moment advancement, PDF
algorithms, mixing length, diffusion/advection and diagnostics, execute in C++.
The [KK](src/Microphys/KK_microphys/) and
[Morrison](src/Microphys/Morrison_microphys/) engines currently execute in host
C++ with Kokkos View staging; their GPU and performance work is deferred.
Statistics stage sampled fields to host-owned buffers; input readers are C++.
The derived-type mirrors also use Views. Host engine snapshots, sampled stats
transfers and owning helper results remain. Persistent runtime fields remove
the outer per-call staging; this alone does not establish a GPU speedup.
Use optimized builds and measure the complete case before drawing performance
conclusions.

The driver also owns a [scratch workspace](src/scratch_workspace.hpp) for
temporary arrays in the core and physics wrappers. Routines that need temporary
storage receive the workspace explicitly and open a scope. Reservations reuse stable storage;
chunks grow as needed and stay allocated until case cleanup or reinitialization.
Normal reservations retain zero initialization. Use `WithoutInitializing` only
when the routine writes every value it subsequently reads:

```cpp
auto workspace_scope = workspace.scope();
Real2D values = workspace_scope.view<Real2D>("routine::values", ngrdcol, nzt);
Real2D result = workspace_scope.view<Real2D>(
    Kokkos::WithoutInitializing, "routine::result", ngrdcol, nzt);
```

The scope releases reservations without a fence. GPU work must use the same
default execution queue so later reuse follows earlier work. Host buffers must
remain reserved until transfers finish. Returned temporary Views reserve in the
caller's active scope through `workspace.view`; they cannot outlive that scope.
Persistent fields and independent snapshots retain their owning allocations.
Case cleanup fences before clearing the workspace. Pointer-wrapped scratch
Views have no individual allocation label; kernel and backing-chunk labels remain.

Storage coordinates, grid endpoints and species/scalar metadata use zero-based
C++ indexing, with columns contiguous in the first dimension. Namelist inputs,
saved ordinals and the temporary Fortran bridge convert at their boundaries.
Model timesteps, printed ordinals and scientific counts retain their meanings.

The [data-lifetime check](auto_llm_generated_tests/driver_data.cpp) exercises
host/device transfers, CPU aliasing, strided windows, batch rebinding and shape
changes. Kokkos Tools callbacks check that synchronized updates do not copy and
same-shape rebinding does not allocate:

```bash
cmake --build build/clubb_kokkos --target clubb_kokkos_driver_data
ctest --test-dir build/clubb_kokkos -R '^clubb_kokkos_driver_data$' --output-on-failure

# Nested scratch lifetimes, warm reuse and initialization behavior.
cmake --build build/clubb_kokkos --target clubb_kokkos_scratch_workspace
ctest --test-dir build/clubb_kokkos -R '^clubb_kokkos_scratch_workspace$' --output-on-failure
```

The regular CPU acceptance target is BFB agreement with the current Fortran
source. Debug builds are used for correctness checks; `-O0` alone does
not guarantee agreement across compilers or CPU/GPU math libraries.
Use the same precision, settings and stats when comparing executables.
The native standalone currently leaves `upwp_sfc_pert` and `vpwp_sfc_pert`
uninitialized in [its driver](../src/clubb_driver.F90). C++ Views initialize
them to zero, so linearized-wind budget output can differ until the native
inputs are defined. With equal inputs, this core path retains BFB agreement.

## Compare with Fortran

The [Kokkos/Fortran case comparison](auto_llm_generated_tests/integration_tests/run_kokkos_vs_fortran_cases.py)
adapts the [JAX comparison](../human_reviewed_tests/run_jax_vs_fortran_cases.py),
using the same 19 enabled cases and their warm-rain/SILHS variants. Build both
executables from the separate CMake project, then run from the repository root:

```bash
cpu_build=build/clubb_kokkos_nvhpc26_cpu_debug
cmake --build "$cpu_build" --target clubb_kokkos_standalone clubb_standalone
comparison=clubb_kokkos/auto_llm_generated_tests/integration_tests/run_kokkos_vs_fortran_cases.py

# Short development checks with the small multicolumn stats registry.
python3 "$comparison" -build_dir "$cpu_build" -max_iters 5 -stats multi_col_stats

# Milestone checks: native timesteps, at most 600 steps, and standard stats.
python3 "$comparison" -build_dir "$cpu_build" -max_iters 600
```

CPU Fortran versus CPU Kokkos at zero tolerance is the regular gate before the
next live integration. Compare the GPU build occasionally after larger pieces
and record the newly integrated piece where results change, even slightly.
Do not investigate GPU numerical causes unless requested; GPU BFB is not
required before continuing. Keep using `multi_col_stats` for
the usual short cycle and standard/all stats at milestones or for diagnosis.
Cap milestone cases at 600 steps; shorter cases retain their native duration.
Omit `-max_iters` when a longer run is needed for a particular investigation.

Use `-build_dir build/clubb_kokkos_nvhpc26_debug` for a GPU milestone comparison.
The harness selects the chosen directory's C++ executable and its
`fortran/clubb_standalone` reference explicitly; it does not use an installed
executable or prepare a JAX environment. That Fortran reference runs on CPU;
the CUDA build's ported kernels run on GPU. Retain the matching module environment
and the GPU selection described above when choosing the CUDA build.

SILHS sample output inherits the native restriction on parameter batching.
For an all-stats SILHS run, omit `-batch_size` so all parameter columns run
together. Batched BOMEX and scalar cases can use the same all-stats registry.

Select cases and forward shared model settings to both executables:

```bash
python3 "$comparison" -build_dir "$cpu_build" -cases bomex -max_iters 5 -dt_main 30 -stats multi_col_stats
python3 "$comparison" -build_dir "$cpu_build" -cases lba_kk lba_kk_silhs -max_iters 5 -stats multi_col_stats -workers 2
python3 "$comparison" -build_dir "$cpu_build" -cases bomex rico_silhs -stats all
python3 "$comparison" -build_dir "$cpu_build" -cases bomex -max_iters 5 -stats multi_col_stats \
  -flag_config_file input/flag_sets/run_bindiff_w_flags_config_example.json
```

`-cases` accepts comma- or space-separated names, including the named
`lba_kk` and `lba_kk_silhs` variants. Other `run_scm.py` model options, such as
`-batch_size`, `-debug` and repeated `-override`, apply to both runs. User
overrides take precedence over case and flag-set settings. The harness uses
four C8 parameter columns. CPU case concurrency defaults to half the available
logical CPUs; `-workers` sets an explicit number of concurrent case pairs. GPU
builds require one worker. CMake writes `kokkos_backends.json` in the selected
build directory so the runner can make this choice, including CUDA packages
that also enable Serial. Reconfigure an older build before using this runner.
The summary records the configured backends and actual worker count. Native
case durations/timesteps remain the default; JAX-specific roundoff caps and
tolerances are not copied.
An explicit main-timestep override must remain compatible with each case's
radiation timestep; for example, GABLS3 night rejects a 30 s main timestep.
The full list includes multi-day simulations and retains both NetCDF output trees.

Results replace `output/tests/kokkos_driver_test_results/`, with separate
`kokkos_output` and `fortran_output` trees, command/run/diff logs under `logs`,
`case_compare_summary.json`, and a final directory-comparison log. Differences,
failed runs, unequal completed-step counts or a failed final audit return a
nonzero status. With the default zero comparison limits, saved NetCDF files
must be byte-identical; a numeric match after the shared comparator clips tiny
fields cannot pass this BFB gate. For diagnosis, an explicit nonzero
`-bindiff_threshold` or `-bindiff_percent_threshold` selects the shared numeric
criteria and is recorded in the summary. The summary still reports byte identity.
The runner has a documented manual route here; it is not wired into Jenkins.

Focused generator and CLI/report checks use the shared pytest wrapper:

```bash
bash human_reviewed_tests/run_pytests.sh -unit \
  clubb_kokkos/auto_llm_generated_tests/pytests/unit_tests -k kokkos
```

## Derived types and temporary calls

The C++ mirrors under [src/CLUBB_core](src/CLUBB_core/) and
[src/SILHS](src/SILHS/) cover the eleven non-stats derived types used by
[clubb_driver.F90](../src/clubb_driver.F90). Statistics have their own
[case-owned C++ types](src/CLUBB_core/stats_netcdf.hpp). The other types' array
fields are `Kokkos::View`s.
Fortran array order is preserved with `LayoutLeft`; a Fortran array of fixed
length strings becomes a two-dimensional character View whose first extent is
the string width. The test-only [f2k converters](src/f2k/derived_type_bridge.hpp) copy
each field through host buffers for Native conversion and numerical oracles.
A Fortran wrapper creates its derived values locally, calls the native routine, and copies updated fields
back before returning. No derived value is kept in saved Fortran storage.
The driver has a column-by-level collection of `hydromet_pdf_parameter` values.
The [C++ collection](src/CLUBB_core/hydromet_pdf_parameter_module.hpp) stores
each fixed-size component in a separate column-first Kokkos View. It has
allocation and reset operations. The C++ driver owns the active collection;
its converters include the outer column and level dimensions and the fixed
species dimensions, including rank-four correlation fields.

After changing one of the Fortran type declarations, run the
[converter generator](src/f2k/generate_derived_types.py) from the repository
root:

```bash
python3 clubb_kokkos/src/f2k/generate_derived_types.py
python3 clubb_kokkos/src/generate_driver_io.py
python3 clubb_kokkos/src/Microphys/generate_microphys_io.py
```

The generator checks the expected number of fields in each type; update that
inventory when adding or removing a component. Eight marked core/SILHS schema
regions preserve the permanent C++ lifecycle, grid, flag and parameter APIs outside
the generated field declarations. Private fields use the accessor pair in
[err_info_type_module.F90](../src/CLUBB_core/err_info_type_module.F90).
The [C++ radiation interface](src/Radiation/radiation_module.hpp) treats `radht`
as an input/output field because the heating rate persists between radiation timesteps.
`ctest --test-dir build/clubb_kokkos --output-on-failure` runs a fieldwise round
trip through all eleven native types and checks the hydrometeor collection's
shape, column-first layout, native conversion and reset. The
[core transfer check](auto_llm_generated_tests/core_array_transfer.cpp) compares
constants with compiled Fortran and checks execution-space aliasing, strided
column windows, rank-three/four fields and zero-scalar arrays. The
[core error check](auto_llm_generated_tests/core_error_return.cpp) verifies that
a NaN forcing returns fatal status through the C++ API before surface fluxes
or lower solvers change state.
The [radiation cadence check](auto_llm_generated_tests/radiation_bridge.cpp)
uses the C++ case state and API to verify retained heating on a timestep that
does not recompute radiation.
The [PDF kernel check](auto_llm_generated_tests/core_pdf_kernels.cpp) compares
all temperature and skewness interfaces and the one-level PDF width branches
directly with compiled Fortran, including small variances and gamma coefficient
limits. It also checks every standalone constant and tunable-parameter index.
The [LU solver check](auto_llm_generated_tests/core_lu_solvers.cpp) compares
all five tridiagonal/pentadiagonal overloads with compiled Fortran and checks
matrix residuals, boundary rows and input preservation. Their C++ matrices use
`(column, band, row)` so columns remain contiguous; the band is zero-based.
The [matrix wrapper check](auto_llm_generated_tests/core_matrix_wrapper.cpp)
compares all five supported LU solve interfaces and their debug/host/error gates.
LAPACK solve methods and condition-estimate (`rcond`) requests return fatal
`err_info` before any solve runs. This includes stats selections requesting
`*_matrix_condt_num`, such as `-stats all`; use `standard` or `multi_col_stats`
for supported comparisons. The
[LAPACK rejection check](auto_llm_generated_tests/core_lapack_diagnostics.cpp)
checks both retained expert interfaces for the unsupported-feature message,
all-column fatal status, input preservation and cleared outputs. The normal LU
solve uses the selected Kokkos execution space. The four main C++ moment parents and supported C++
microphysics callers use these solver wrappers.
The [Cholesky check](auto_llm_generated_tests/core_cholesky.cpp) runs through
CTest in either build mode. It checks factor reconstruction, lower-triangle
ownership, scaling, retries, invalid pivots and eigenvalue diagnostics, and runs
the portable factorization kernel on the configured execution space. The
[configuration checks](auto_llm_generated_tests/pytests/unit_tests/test_kokkos_bfb_configuration.py)
exercise missing-backend warnings without a Fortran compiler; the comparison
runner records requested and effective BFB mode in its saved report.
The [grid operator check](auto_llm_generated_tests/core_grid_operators.cpp)
covers interpolation, smoothing, gradients and reversal using grids constructed
by the Native reference. Production grid construction, height updates and
interpolation weights, grid-file reading and validation are C++.
The [interpolation check](auto_llm_generated_tests/core_interpolation.cpp)
compares all eight scalar/profile routines with compiled Fortran on the host
and GPU, including signed slopes, source endpoints, tolerance matches, zero
fill, clamping and strided column profiles. It also checks the host input guards.
The [calendar check](auto_llm_generated_tests/core_calendar.cpp)
(`clubb_kokkos_core_calendar`) covers Gregorian century leap rules, Julian
conversions, signed elapsed-time rollovers and the Native upper-day error branch.
The [parameter-validator check](auto_llm_generated_tests/core_parameter_validation.cpp)
(`clubb_kokkos_core_parameter_validation`) compares five-column hard-bound and
C6-equality errors with compiled Native validation, preserving prior error codes
and Native NaN-comparison behavior.

The generic procedure bridge and its generator have been retired. The remaining
[Fortran declaration reader](src/f2k/fortran_procedures.py) supplies the permanent
C++ namelist generators. Generation needs Python 3 and the C preprocessor `cpp`
and uses the standalone preprocessing profile. Defaults, interpretation and
mutable physics state belong to the C++ owners.

Configuration printing uses caller-owned C++ streams in
[model_flags](src/CLUBB_core/model_flags.cpp) and
[parameters_silhs](src/SILHS/parameters_silhs.cpp). The Native labels, field order
and logical values are retained; an output stream replaces the Fortran unit.
The [error-context API](src/CLUBB_core/err_info_type_module.cpp) formats the
214-character global and 200-character column headers in C++, preserving Native
fixed-width padding and overflow markers. Coordinate rounding matches the recorded
NVHPC 26.9 reference, as explained in [the adaptation log](BFB_ADAPTATIONS.md#diagnostic-header-formatting). Headers remain
case-owned; updating them preserves the column allocation and error codes.
The unused Native module-wide process/thread prefix has no C++ copy or registry.
The [debug-level API](src/CLUBB_core/error_code.cpp) retains the Native setter and
threshold comparison, including signed levels. Each host thread starts at level
zero. Physics and statistics checks use the same C++ accessor; no Native
threshold synchronization remains.
The [debug check](auto_llm_generated_tests/debug_level.cpp) compares signed and
extreme thresholds with Native and checks thread isolation and setter ownership:

```bash
cmake --build build/clubb_kokkos --target clubb_kokkos_debug_level
ctest --test-dir build/clubb_kokkos -R '^clubb_kokkos_debug_level$' --output-on-failure
```

The [diagnostic output check](auto_llm_generated_tests/config_error_output.cpp)
compares header edge cases with Native and checks configuration stream state:

```bash
cmake --build build/clubb_kokkos --target clubb_kokkos_config_error_output
ctest --test-dir build/clubb_kokkos -R clubb_kokkos_config_error_output --output-on-failure
```

## Timers and timing profiles

The [timer module](src/CLUBB_core/code_timer_module.hpp) keeps the Fortran
initialize/disable/start/stop/finalize API and named regions. Its default backend
records process CPU time, with nested inclusive/exclusive totals and call counts.
`USE_GPTL=ON` calls GPTL directly through its C API with the Native full-tree
settings. Both backends write the existing `.timing` formats consumed by
[time_clubb.py](../utilities/time_clubb.py); no parser changes are required.

Run the timing script from the repository root with the Python environment used
by the case runner. This short example includes a warmup and two measured runs:

```bash
python3 utilities/time_clubb.py bomex \
  -exe "$PWD/build/clubb_kokkos/clubb_standalone" \
  -max_iters 3 -warmups 1 -repeats 2 -name kokkos_short
```

For a process/batch-size sweep, the same executable remains alive through the
warmup and repetitions for each workload. Each request resets its timer state:

```bash
python3 utilities/time_clubb.py bomex \
  -exe "$PWD/build/clubb_kokkos/clubb_standalone" \
  -process_counts 1,2 -batch_sizes 1,2 \
  -max_iters 3 -warmups 1 -repeats 2 -name kokkos_sweep
```

The script writes compact profiles under `output/timing`, including timer rows,
workload metadata and representative raw timing files. Choose a new name or use
its `-overwrite` option to replace a previous profile. The explicit executable
path and digest identify the tested C++ build; no upstream Kokkos selector is
needed. Use an optimized build for performance measurements; these short Debug
examples check the workflow.

To enable GPTL in an existing configured build:

```bash
cmake -S clubb_kokkos -B build/clubb_kokkos -DUSE_GPTL=ON
cmake --build build/clubb_kokkos --target clubb_kokkos_standalone
```

The existing dependency uses an available installation or builds GPTL 8.1.1.
Set `GPTL_ROOT` when supplying a compiler-compatible installation; C++ needs
`gptl.h` and the C library, while the retained Fortran build needs its module and
bindings. Use `-DUSE_GPTL=OFF` for the process-CPU-time backend. GPTL profiles
report wall time; CPU-time profiles retain their CPU-time basis.

Private per-thread bookkeeping mirrors Native threadprivate state, as an
intentional exception to case-owned model storage. GPTL remains a shared library
service; standalone requests run sequentially within each worker process.
Neither backend adds GPU
fences. Timers around asynchronous regions therefore include the host work and
existing waits; they do not promise isolated device-kernel durations. Model BFB
checks compare saved model fields, not elapsed timing values.

The small [timer check](auto_llm_generated_tests/code_timer.cpp) covers lifecycle,
disabled behavior, nesting errors, registry/stack growth and output headers:

```bash
cmake --build build/clubb_kokkos --target clubb_kokkos_code_timer
ctest --test-dir build/clubb_kokkos -R clubb_kokkos_code_timer --output-on-failure
```

## Input readers

The C++ readers accept the existing standalone input files. The shared
[namelist parser](src/namelist.hpp) applies assignments to initialized scalars
and Kokkos Views. [Driver bindings](src/driver_io.hpp),
[microphysics bindings](src/Microphys/microphys_io.hpp) and
[radiation input](src/Radiation/radiation_io.cpp) preserve Native field names,
array order, character widths and defaults. The first two binding files are
generated from the Native declarations by
[generate_driver_io.py](src/generate_driver_io.py) and
[generate_microphys_io.py](src/Microphys/generate_microphys_io.py).
Run them from the repository root after changing the corresponding Native schema:

```bash
python3 clubb_kokkos/src/generate_driver_io.py
python3 clubb_kokkos/src/Microphys/generate_microphys_io.py
python3 clubb_kokkos/src/Microphys/generate_microphys_io.py
```

Namelists support quoted strings, repeated and null values, indexed assignments,
contiguous array sections, substrings and derived-type components. Missing or null
entries preserve initialized values. As with the installed NVHPC reader, character
values require quotes, explicit section strides are rejected, and real-spelled
integer input is truncated. Reader status uses portable categories: zero for
success, negative for a missing group/end of file, positive for an input error;
vendor-specific IOSTAT numbers are not reproduced. Field and value errors include
file, group and line context. A generic `/profile/` requires a positive `nlevels`; the Native
omitted-count case has undefined storage. Unspecified profile array slots are
initialized to zero rather than copied from indeterminate storage.

[Formatted-file input](src/CLUBB_core/file_functions.hpp) keeps physical-record
READ behavior, including discarded record tails, repeated/null values and the
empty READ in `file_read_2d`. Grid checks retain the Native error paths. The
[table reader](src/Input_fields/input_reader.cpp) retains header widths and
dimension metadata; a two-dimensional table with no data block is rejected
rather than using the Native undefined row count.

[NetCDF input](src/netcdf.hpp) calls NetCDF-C directly. Its logical interface
keeps one-based variable/dimension IDs and fastest-axis-first Native metadata
and slabs; conversion to C conventions happens once inside that service.
Readers retain their existing scientific conversions and View staging. Output
statistics use the C++ backend below. Input buffers, parser
bindings and file handles are owned by their call; no saved bridge copy or
global file-unit registry is used.

The small [namelist](auto_llm_generated_tests/namelist.cpp),
[text/grid](auto_llm_generated_tests/input_text_io.cpp) and
[NetCDF](auto_llm_generated_tests/netcdf_input.cpp) checks cover input edge cases:

```bash
cmake --build build/clubb_kokkos --target clubb_kokkos_namelist clubb_kokkos_input_text_io clubb_kokkos_netcdf_input
ctest --test-dir build/clubb_kokkos -R 'clubb_kokkos_(namelist|input_text_io|netcdf_input)$' --output-on-failure
```

Use the maintained whole-model comparison above for numerical acceptance.
The input replacement matches all 19 five-step CPU cases with multi-column stats,
short IOP/restart/cloudy MPACE-A checks, and all 360 native BOMEX steps with standard
stats in `KOKKOS_BFB=ON` mode. The earlier full 19-case driver evidence predates
this input replacement; it is not a new full-matrix result for these readers.

## Statistics

The [statistics API](src/CLUBB_core/stats_netcdf.hpp) collects named fields,
pointwise means and budgets and writes one classic NetCDF file per case. The
[diagnostic accumulation routines](src/CLUBB_core/stats_clubb_utilities.cpp)
compute the same registered quantities as Fortran. Both run in C++; each case
owns its buffers, lookup caches, grid inputs and file identifiers.

Initialize with `stats_init_api`, bracket model steps with
`stats_begin_timestep_api`/`stats_end_timestep_api`, and close with
`stats_finalize_api`. Use `stats_update` for scalar, surface or profile fields
and the begin/update/finalize budget calls for tendencies. The driver already
performs this sequence. `StatsInitOptions` groups Native optional metadata.
An empty output path retains sampling and accumulation without writing a file.

Persistent buffers use double precision and column-first `LayoutLeft` host
Views. Input Views can reside on the host or device and can be strided; sampled,
registered updates take independent snapshots in the stats context's private
scratch workspace. Buffers become reusable after their update or write completes;
value snapshots returned by `get_stats` retain independent owning storage.
NetCDF-C receives the
Native definition and write order, metadata types, fill values and packed slabs.
Column and level arguments are zero-based; model timestep indices remain
one-based. Sampling/output periods and model times are in seconds, and grid
heights are in meters. The sampling window is `(stats_tstart, stats_tend]`.
Repeated sub-timestep updates use their pointwise mean as one timestep sample.
For rank-one `stats_update`, pass `std::nullopt, true` after the stats context
to request that mean for all columns.

`stats_reset_api` and `start_next_stats_batch_api` reuse the file across column
batches. `get_stats` returns a deep value snapshot that cannot mutate or close
the live file. Optional output-grid remapping and SILHS subcolumn output use
the permanent C++ implementations. NetCDF calls run sequentially on the host;
concurrent host writers are outside the current standalone execution contract.
The registry-overflow warning is once per case-owned context, including resets
and reinitialization; Native warns once per process. Truncation and registry
ordering follow Native.

The [stats CTest](auto_llm_generated_tests/stats_netcdf.cpp) sends execution-space
fields through two contexts. A Native NetCDF reader checks exact averages,
budgets, sub-timestep means, batch columns, metadata and SILHS sample values.
It also checks independent snapshot buffers, read-only file ownership and exact
equal-mass output-grid remapping:

```bash
cmake --build build/clubb_kokkos --target clubb_kokkos_stats_netcdf
ctest --test-dir build/clubb_kokkos -R '^clubb_kokkos_stats_netcdf$' --output-on-failure
```

Use the [whole-model comparison](#compare-with-fortran) for CPU acceptance.
Host accumulation and staging are the initial implementation; GPU performance
work is deferred.
The completed stats milestone matches all 19 native-duration/native-timestep
CPU cases byte-for-byte with standard stats, real8 Debug/O0, four C8 columns
and `KOKKOS_BFB=ON`. The driver plan records additional batch, scalar, remap,
sampled-microphysics and GPU checkpoint coverage.

The [driver sequence check](auto_llm_generated_tests/clubb_driver_sequence.cpp)
reinitializes a case, advances once with stats suppressed, resets, and advances
in successive windows. Run it with a single runtime batch and at least two
timesteps, then compare its stats with a normal standalone run:

```bash
cmake --build build/clubb_kokkos --target clubb_kokkos_driver_sequence
./run_scripts/run_scm.py -exe build/clubb_kokkos/clubb_kokkos_driver_sequence \
  -max_iters 2 -output_dir output/kokkos_driver_sequence bomex
```
