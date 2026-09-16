# Flag sets

The JSON files in this directory define **flag sets**: groups of CLUBB model flags to run with.
They are used by tests that check CLUBB under many flag settings (bindiff and JAX-vs-Fortran called by Jenkins).

We can't test CLUBB with every combination of flags because the number of combinations explodes.
Instead in order to reduce the number of combinations, we simultaneously flip the values of a set of flags.  
The sets are chosen so that each flag is tested both true and false across the test runs.
One of the two states may be covered by the `default` run.

- `run_bindiff_w_flags_config_core_flags.json`: core model flags, used by most Jenkins tests.
- `run_bindiff_w_flags_config_host_flags.json`: flag settings used by the generalized vertical grid host-flags test.
- `config_core_flags_wo_never_supported_jax_features.json` : core model flags to test clubb_jax vs fortran clubb without flags  
that jax per design will never support. E.g. LAPACK configurations.

## What is a flag?

A flag is a configuration switch that picks **which code path or scheme** CLUBB uses.
It is either a logical (`l_*`) or an integer method selector (e.g. `saturation_formula`, `iiPDF_type`, `penta_solve_method`).  
A tunable parameter (e.g. `C1`, `C4`) on the other hand is a real-valued coefficient that only **scales a term** in the equations.

Example: `l_diffuse_rtm_and_thlm = .true.` adds an extra diffusion term for `rtm` and `thlm`.
`saturation_formula = 1` or `2` switches to a different saturation vapor pressure formula.

In the namelist, flags and parameters are plain entries in different namelist groups.
Their defaults live in separate files in `input/parameter_and_flag_configs/<config>/`:

| File                          | Namelist groups                                                  |
|-------------------------------|------------------------------------------------------------------|
| `configurable_model_flags.in` | `&configurable_clubb_flags_nl`, `&configurable_silhs_flags_nl`   |
| `tunable_parameters.in`       | `&clubb_params_nl`                                               |
| `silhs_parameters.in`         | `&silhs_params_nl`                                               |

These files are only the source files.
Before a run they are combined into one namelist file (see below).
Fortran reads each group from that one file with `read(unit, nml=<group>)`.

In Fortran, flags end up in `clubb_config_flags_type` (`src/CLUBB_core/model_flags.F90`).
Parameters end up in the `clubb_params` array, which can differ per column and is what the tuner changes.
Overrides from a flag set work the same way for both kinds.


## How flags get to the model

1. `run_scripts/run_clubb_w_varying_flags.py -f <file>.json` (or `tests/run_jax_vs_fortran_cases.py --flag-config-file <file>.json`)
   reads the JSON with `utilities/flag_sets.py`.
   An unmodified `default` run is always added unless `--skip-default-flags` is given, so `default` can't be used as a flag set name.
   Each flag set becomes one `-override` string, e.g. `l_diffuse_rtm_and_thlm=.true.,saturation_formula=2`.
2. For every (flag set, case) pair it calls `run_scripts/run_scm.py -override ... <case>`, writing to `output/<flag_set_name>/`.
3. `run_scm.py` calls `utilities/create_case_namelist.py`, which builds one combined `<case>.in` namelist.
   It pastes together the text of the three config files above, `input/case_setups/<case>_model.in` and the stats file.
   Then it find-and-replaces each override key in that combined text.
   An override key that doesn't match an entry is an error.
4. That namelist is passed to the executable: `clubb_standalone` (Fortran) by default, or the JAX driver with `-jax`.

`tests/run_bindiff_w_flags.py` wraps step 1 for several git refs and bindiffs the per-flag-set outputs.


## Notes

`l_Lscale_plume_centered: true` was removed from `run_bindiff_w_flags_config_core_flags.json`.
It errors because the compile-time constant `l_avg_Lscale` in `src/CLUBB_core/mixing_length.F90` is set to false.
For more context see https://github.com/larson-group/clubb/issues/1354#issuecomment-5640458314
