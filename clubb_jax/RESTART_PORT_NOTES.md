# JAX restart port

Standalone restarts retain their initialized full column count; combining a
restart with smaller runtime batches is rejected until the reader selects each
batch's saved column slice. The reset snapshot includes restored fields
and SILHS permutation state. Driver windows preserve the absolute restart
iteration by default. Managed tuning and loss evaluation are separate PRs.

The JAX standalone reads the existing CLUBB NetCDF statistics format during
`init_clubb_case`. It resumes with the original `time_initial`, `time_final`
and absolute iteration number. No compiled Fortran library or new checkpoint
format is required.

```bash
python3 tests/run_restart_test.py bomex -jax -multicol 4
python3 tests/run_restart_test.py rico_silhs -jax
```

The script runs an uninterrupted case and another from the saved interior
record nearest its midpoint, then compares the selected final statistic bit
for bit across every column.
Fortran runs retain the original first-column check; the source restart
reader does not restore distinct column states.
`-var rcm` selects cloud water; `-jax=cpu` and `-jax=gpu` use the normal launcher's
backend selections. Forwarded duration/timestep overrides are included when
determining the restart time. Both runs retain user physics overrides.

## Fortran source and field inventory

- `src/clubb_driver.F90:init_clubb_case`, especially its restart block:
  initialize sounding/reference fields first, restore model state, recompute
  reciprocal dry densities and start at
  `floor((time_restart-time_initial)/dt_main)+1`.
- `src/clubb_driver.F90:restart_clubb`: the port retains all 90 arguments in
  source order, scheme-dependent input flags, eight previous microphysics
  tendencies and four surface fluxes taken from the lowest momentum level.
- `src/Input_fields/input_fields.F90:set_filenames`, the CLUBB branch of
  `stat_fields_reader`, `compute_timestep` and
  `get_clubb_variable_interpolated`: the reader retains 79 arguments and all
  96 source field-read calls, in order. Thermodynamic/momentum fields,
  selected hydrometeors, PDF fields and optional soil temperatures follow
  the existing Fortran selection. Variances retain the source lower bounds.
- The host NetCDF boundary includes the relevant `input_netcdf.F90` and
  `stat_file_utils.F90` reads, unit conversion and linear interpolation.

`advance_clubb_to_end` only changes its starting iteration, remaining step
count and statistics iteration. Forcings, radiation cadence and SILHS seeds
use the resulting absolute iteration without changing their formulas.

## Language adaptations and limits

- NetCDF4 replaces Fortran file handles. Profile arrays include all columns;
  immutable JAX arrays and PDF fields are returned rather than mutated.
  Each saved column is restored separately. The legacy Fortran reader reads
  the first column repeatedly; a single-column input can still be broadcast.
- The reader supports split `_zt.nc`/`_zm.nc`/`_sfc.nc` files and unified
  `_stats.nc` files. It copies matching grids directly, preserving bits;
  remapped grids use the source linear interpolation formula. Valid zeros
  are read raw because CLUBB uses `_FillValue=0`. File/model date offsets and
  one-based record indices are retained.
- Restart time must align with the model timestep and identify a saved
  record. Actual time coordinates select the record, including files with
  delayed output windows; the source instead infers an index from the date
  origin and output interval and additionally requires minute alignment.
- Relative `restart_path_case` values use the repository root; absolute
  prefixes are accepted by JAX. Initialization rejects overwriting its
  input statistics file and closes the output handle if restoration fails.
- SILHS keeps the existing native JAX random generator. The absolute
  iteration restores timestep seeds. Initialization reconstructs a retained
  multi-timestep permutation from the last reshuffle's seed and final column
  key, preserving its prior iteration without replaying model physics.
- Bit-for-bit continuation requires instantaneous statistics containing the
  source-selected restart fields. Use `standard_stats.in` or `all_stats.in`
  with `stats_tsamp=stats_tout=dt_main`; reduced registries may omit required
  fields. Missing fields produce a read error. Averaged/remapped input is an
  initialization state, not an exact checkpoint of an earlier trajectory.
- This ports CLUBB restart input only. General LES `l_input_fields` remains
  unsupported. Existing core/microphysics/radiation support limits still
  apply. Fields absent from the source CLUBB reader are not silently added;
  active frozen-phase Morrison restart parity is not established here.
- The source reader omits passive-scalar and eddy-scalar state and frozen
  number concentrations such as `Nsm` and `Ngm`. Those configurations do not
  have established exact continuation; omitted fields retain initialization
  values. Adding them requires extending the Fortran restart contract first.
- The source reader restores `radht`, but omits radiation fluxes and separate
  shortwave/longwave heating caches. With `dt_rad > dt_main`, resumed output
  can differ until the next radiation update. A six-step DYCOMS RF01 check
  with `dt_main=60`, `dt_rad=180` and a restart after step 3 matches core
  fields, but the next two records lose longwave heating/flux diagnostics.
  Soil/vegetation feedback in this configuration is unvalidated. Use
  `dt_rad=dt_main` for the established bit-for-bit restart tests.

## Validation

The shared `tests/run_restart_test.py` is the real-case restart test for both
backends. It compares the selected final-timestep statistic bit for bit; JAX
checks every column, while Fortran retains its existing first-column check.
Broader all-statistics/all-record comparisons are deferred to a future shared
harness improvement.

The `clubb_restart` Jenkins job adds two JAX variants beside its native tests:

- BOMEX: the full 360-step run at 60 seconds, with four distinct C8 columns.
  Restart from step 180 and compare the final `thlm` in every column.
- RICO SILHS: 74 steps at 60 seconds, two distinct C8 columns and sequence
  length three. Restart from step 37, so step 38 reuses a retained permutation.
  Importance sampling is disabled for this focused sequence check.

Focused provisional pytests cover the host reader and initialization boundary:
raw zeros, units, dimension order, date offsets, saved records, split/unified
files, interpolation, single-column broadcasting, invalid dtype/units, clock
validation and protection against overwriting the reference file. CLI checks
cover backend forwarding, effective duration and comparison of later columns.
These tests do not add another model runner or change the native criteria.

Validation is on CPU with double precision. Hardware GPU execution and exact
continuation for the unsupported/source-omitted configurations above are not
established. Historical broader checks on predecessor revisions are recorded
in the session worklog; they are not additional maintained restart workflows.

On 2026-10-06, the two JAX commands above pass through the shared harness at
unchanged exact comparison criteria. The focused restart suite passes all 52
checks; the shared plain-Python suite passes 304 checks. A separate ten-step
native BOMEX multicolumn run also passes its original first-column criterion.
