# JAX loss-driver port

Canonical source: `src/clubb_loss_driver.F90`. This branch adds standalone
JAX loss evaluation without a Fortran library. Managed tuner integration is
a separate functionality PR.

Source routine outline, in order:

1. `is_finite_core_value`, `set_invalid_field_metric_outputs`
2. `calculate_taylor_metrics`, `stop_with_error`
3. `init_loss_request` (namelist, validation, windows, field descriptors,
   stats bindings and benchmark preparation)
4. `init_clubb_loss`, `get_loss_time_window_count`
5. `clubb_get_loss_for_params` (reset each runtime batch, advance successive
   windows, score native-grid stats, penalize failed columns)
6. `finalize_clubb_loss`, `clubb_get_loss`
7. `prepare_loss_request_for_scoring` (time units, field/grid lookup, height
   bounds, interpolation, open/closed benchmark windows)
8. `calculate_field_loss` (range-scaled squared error and Taylor metrics)

Language adaptations:

- Derived types become dataclasses; Python returns output arguments and raises
  exceptions for fatal configuration errors. Array and field indices are zero based.
- NetCDF and namelist reads, lifecycle and runtime batch slicing stay on the host.
  Profile arithmetic and metrics use JAX. Fixed runtime shapes and module-level
  jitted functions keep parameter values out of compilation cache keys.
- Driver initial conditions are an immutable JAX snapshot instead of repeated
  Fortran assignments. Parameter-dependent derived quantities are recomputed
  on reset, matching `set_case_initial_conditions` in `src/clubb_driver.F90`.
- The loss driver reads the JAX stats banks and sample counts directly, before
  the next window reset; file output is optional. It scores native zt/zm grids.
- A final short batch is padded with a valid candidate and trimmed on return,
  so a smaller candidate count does not create another compiled model shape.
  Variable candidate counts require `stats_output_filename = ""`: NetCDF files
  retain the configured total-column dimension. Full matrices use the source's
  batch numbering so each runtime batch writes its own output-column slice.
- Invalid candidates are replaced by defaults during model advancement and
  receive finite penalties afterward. This host boundary preserves healthy
  neighboring columns while avoiding the source's process-wide `lmin` stop.
- The source derives the global `lmin` and `Skw_max_mag` values from the last
  column in each runtime batch. Varying these parameters can affect neighboring
  columns and make scores depend on batch width; their evaluation is not
  column independent. The JAX port retains this source behavior.
- Microphysics host error checks retain the source's `debug >= 0` guards.
  At `debug = -1`, fatal flags remain per-column outputs for loss penalties,
  allowing healthy neighbors to finish their loss windows.
- NetCDF reads disable automatic masking and scaling, as `nf90_get_var` reads
  raw values. Unit conversions follow the source. Benchmark interpolation calls
  the existing JAX `lin_interpolate_two_points` in source argument order; exact
  levels copy directly. Profile/time reductions vectorize source accumulations.
- JAX scores its active accumulation bank after dividing by sample counts.
  Each loss subwindow must therefore coincide with one complete stats output
  window. The source's printed benchmark timing diagnostic is retained.
- Loss evaluation with `l_restart` is rejected before case initialization.
  The source loss loop begins at timestep 1 even when the driver has a restart;
  choosing valid post-restart loss windows requires a separate lifecycle change.
- Existing standalone unsupported-feature gates apply equally to loss
  evaluation.

Validation on CPU includes source-order/metric/read contracts, the native
Fortran loss oracle, repeated runtime parameter batches and compilation-cache
reuse. The normal-versus-loss consistency script runs both backends with one
full-window output record and checks saved profiles and independently
recomputed metrics at the existing criteria. Shared tuner orchestration is a
separate functionality PR; GPU execution is not validated here.
