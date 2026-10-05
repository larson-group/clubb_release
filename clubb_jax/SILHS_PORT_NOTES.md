# SILHS JAX port

Canonical sources are `src/SILHS` and the sampled-microphysics routines in
`src/Microphys`. Native JAX random draws and permutations intentionally replace
MT95 and CLUBB random-generator machinery. Numerical sampling transformations,
sequence-length permutation reuse, overlap, clipping, weighting and feedback
follow the Fortran. Shape metadata/configuration are static; permutation arrays
and random keys are explicit JAX values. Output/inout arguments become returns.
Initialization owns configuration and initial sampling storage. The timestep
adapter only passes evolving state and the current iteration.
Initialization also calls the SILHS host output API after opening statistics.
As in Fortran, enabled SILHS defines `lh_zt` and `lh_sample_number` even when
the statistics list contains no SILHS fields and optional sample output is off.

Native random-stream matching with Fortran is deferred. The deterministic
comparison mode below supplies common sampling inputs; native randomness is
checked with invariant/statistical tests.
SILHS radiation remains outside the microphysics port.

## Fortran routine outline

### src/SILHS/est_kessler_microphys_module.F90

- Line 18: `subroutine est_kessler_microphys_api`
- Line 267: `subroutine calc_estimate( num_samples, mixt_frac,`

### src/SILHS/generate_uniform_sample_module.F90

- Line 20: `function rand_uniform_real( )`
- Line 68: `subroutine generate_uniform_lh_sample( iter, num_samples, sequence_length, n_vars,`
- Line 160: `function choose_permuted_random( nt_repeat, p_matrix_element )`
- Line 201: `subroutine permute_height_time( nt_repeat, n_vars, one_height_time_matrix )`
- Line 240: `subroutine rand_permute( n, pvect )`

### src/SILHS/latin_hypercube_arrays.F90

- Line 20: `subroutine cleanup_latin_hypercube_arrays( )`

### src/SILHS/latin_hypercube_driver_module.F90

- Line 27: `subroutine generate_silhs_sample(`
- Line 489: `subroutine generate_random_pool( nzt, ngrdcol, pdf_dim, num_samples, d_uniform_extra,`
- Line 632: `subroutine generate_all_uniform_samples(`
- Line 938: `subroutine compute_k_lh_start( gr, nzt, ngrdcol, rcm_pdf, pdf_params,`
- Line 1142: `subroutine clip_transform_silhs_output( nzt, ngrdcol, num_samples,           & ! In`
- Line 1330: `subroutine assert_consistent_cloud_frac( chi_1, chi_2,`
- Line 1394: `subroutine assert_consistent_cf_component( mu_chi_i, sigma_chi_i, cloud_frac_i,`
- Line 1528: `subroutine assert_correct_cloud_normal( num_samples,`
- Line 1637: `subroutine latin_hypercube_2D_output_api( nzt, zt, pdf_dim, num_samples, hm_metadata,`
- Line 1785: `subroutine compute_arb_overlap( nzt, ngrdcol, num_samples, pdf_dim, d_uniform_extra,`
- Line 1913: `subroutine stats_accumulate_lh_api(`
- Line 2214: `subroutine stats_accumulate_uniform_lh( nzt, num_samples, ngrdcol, l_in_precip_all_levs,`
- Line 2410: `subroutine copy_X_nl_into_hydromet_all_pts(`

### src/SILHS/lh_microphys_var_covar_module.F90

- Line 15: `subroutine lh_microphys_var_covar_driver_api(`

### src/SILHS/math_utilities.F90

- Line 18: `pure function compute_sample_mean( n_levels, n_samples, ngrdcol,`
- Line 68: `pure function compute_sample_variance( n_levels, n_samples, ngrdcol,`
- Line 119: `pure function compute_sample_covariance( n_levels, n_samples, ngrdcol,`
- Line 175: `function rand_integer_in_range(low, high)`

### src/SILHS/parameters_silhs.F90

- Line 95: `subroutine set_default_silhs_config_flags_api( cluster_allocation_strategy,`
- Line 171: `subroutine initialize_silhs_config_flags_type_api( cluster_allocation_strategy,`
- Line 245: `subroutine print_silhs_config_flags_api( iunit, silhs_config_flags )`

### src/SILHS/silhs_api_module.F90

- Line 148: `subroutine generate_silhs_sample_api(`
- Line 312: `subroutine clip_transform_silhs_output_api(`

### src/SILHS/silhs_importance_sample_module.F90

- Line 31: `subroutine importance_sampling_driver`
- Line 224: `function define_importance_categories( )`
- Line 283: `function compute_category_real_probs( importance_categories,`
- Line 376: `function compute_category_sample_weights( category_real_probs, category_prescribed_probs )`
- Line 429: `subroutine limit_category_weights( category_real_probs, category_prescribed_probs )`
- Line 555: `function pick_sample_categories( num_samples, category_prescribed_probs,`
- Line 685: `subroutine scale_sample_to_category( category, cloud_frac_1, cloud_frac_2,`
- Line 788: `function two_cluster_cp_nocp( importance_categories, category_real_probs,`
- Line 884: `function eight_cluster_allocation( importance_categories, category_real_probs,`
- Line 990: `function four_cluster_no_precip( importance_categories, category_real_probs,`
- Line 1115: `function compute_clust_category_probs`
- Line 1178: `function clust_cat_probs_frm_var_fracs`
- Line 1280: `function clust_cat_probs_frm_presc_prb`
- Line 1433: `function cloud_importance_sampling( importance_categories, category_real_probs,`
- Line 1522: `subroutine importance_sampling_assertions`
- Line 1699: `subroutine cloud_weighted_sampling_driver`
- Line 1861: `function generate_strat_uniform_variate( num_samples )`
- Line 1919: `subroutine choose_X_u_scaled`
- Line 2061: `function determine_sample_categories( num_samples, pdf_dim, hm_metadata,`

### src/SILHS/transform_to_pdf_module.F90

- Line 15: `subroutine transform_uniform_samples_to_pdf(`
- Line 195: `subroutine cdfnorminv( pdf_dim, nzt, ngrdcol, num_samples, X_u_all_levs,`
- Line 286: `function ltqnorm( p_core_rknd )`
- Line 490: `subroutine multiply_Cholesky( nzt, ngrdcol, num_samples, pdf_dim, std_normal,`
- Line 581: `subroutine chi_eta_2_rtthl( nzt, ngrdcol, num_samples,`

### src/Microphys/lh_microphys_driver_module.F90

- Line 15: `subroutine lh_microphys_driver(`

### src/Microphys/estimate_scm_microphys_module.F90

- Line 15: `subroutine est_silhs_tndcy(`
- Line 371: `subroutine adjust_KK_src_means(`

### src/Microphys/silhs_category_variance_module.F90

- Line 15: `subroutine silhs_category_variance_driver(`
- Line 89: `subroutine silhs_sample_category_variance(`

### src/Microphys/pdf_hydromet_microphys_wrapper.F90

- Line 13: `subroutine pdf_hydromet_microphys_prep`

## Configuration and ownership

Use the existing namelist controls: `lh_microphys_type="interactive"` feeds
sampled tendencies back to CLUBB; `"non-interactive"` collects sampled
microphysics diagnostics while retaining the grid-mean physics; `"disabled"`
uses the ordinary path. Supported schemes are `khairoutdinov_kogan` and
`morrison`. Interactive KK requires `l_local_kk=.true.`. For example:

```bash
./run_scripts/run_scm.py -jax -max_iters 3 -stats none rico_silhs
```

`lh_num_samples`, `lh_seed`, the importance-sampling flags, cluster strategy,
weight normalization, overlap and variance/covariance flags are read during
initialization. Counts must be positive. Importance sampling requires
`lh_sequence_length=1`, as in Fortran. With importance sampling disabled,
longer sequences preserve the full permutation between reshuffles and reuse
its first `lh_num_samples` rows each timestep, exactly as the source does.

`sampling_state` is a case-owned pytree holding the permutation and prior
iteration. Initialization allocates it; the PDF wrapper returns its updated
value each step. There is no mutable JAX generator or module-level permutation.
Per-timestep native keys use the source `lh_seed * itime` seed convention;
separate draw/permutation/column keys replace the original random stream.
Restarts using multi-timestep sequences must retain this sampling state.

### Deterministic comparison inputs

`configurable_silhs_flags_nl.l_lh_deterministic_test=.true.` enables the same
test inputs in Fortran and JAX. `generate_uniform_lh_sample` constructs ordered
permutations and uses the hard-coded cycle `(0.125, 0.625, 0.375, 0.875)` as
within-stratum offsets, retaining the existing sequence storage/reuse. Its
zero-based cycle index is `(sample + variate + iteration - 1) % 4`.
`generate_random_pool` uses the same cycle with zero-based index
`(column + sample + level + variate) % 4`. These binary-exact fractions avoid
endpoints; explicit indices avoid dependence on array storage order.
The low-level random helpers remain unchanged. This mode defaults to false and
requires `l_lh_importance_sampling=.false.` and `l_random_k_lh_start=.false.`;
their independent draws are not covered by these two sampling-level branches.

The inverse-normal transform also preserves Fortran's coefficient precision:
its default-real literals round to float32 before promotion to the input dtype.
Using double-precision literals directly introduced small sample differences
that amplified during the LBA comparisons. A live source-routine oracle checks
both central and tail regions, in eager and compiled JAX execution.

The default comparison suite uses 60-second timesteps, interactive sampled
microphysics, four C8 columns and standard statistics:

| Case | Steps | Samples | Microphysics | Percentage tolerance |
| --- | ---: | ---: | --- | ---: |
| `rico_silhs` | 360 | 8 | Native local KK | 1e-7% |
| `lba_kk_silhs` | 275 | 64 | Local warm-rain KK variant of LBA | 1e-7% |
| `lba_silhs` | 240 | 64 | Native Morrison, ice/graupel enabled | 1e-3% |

All cases retain the existing absolute tolerance of 1e-7; Morrison retains
the existing float32 percentage policy. With the repeating cycle and matching
inverse-normal coefficient precision, fresh RICO and LBA KK comparisons pass
all 360 and 275 timesteps against both Debug/O0 and Release/O2 native builds.
Morrison LBA passes all 240 steps against Release/O2, the normal native build
configuration. Against Debug/O0 it first fails at saved prefix 78, with 17
fields failing the full comparison. The native Debug and Release runs also
exceed the same tolerance in 17 fields when compared directly without JAX;
this case remains sensitive to compiler arithmetic. No tolerance or timestep
limit was changed to accommodate the debug build.
These are curated comparison durations, not full native-duration LBA runs.
Native random-stream parity remains untested.

Run these cases with
`./tests/run_jax_vs_fortran_cases.py -cases rico_silhs lba_kk_silhs lba_silhs -workers 1`.
Fortran reports a tiny negative incoming `rtm` in RICO at debug level 2,
but explicitly clears that input error and continues normally. RICO's native
300-second timestep completes 60 steps but its `Nrp2` comparison exceeds
tolerance; the default test uses the passing 60-second timestep.

Among the other native interactive SILHS cases, LBA's Morrison ice/graupel
and prescribed radiation paths are supported. `arm_97`, `mc3e`, and `twp_ice`
require unported BUGSrad radiation; `mpace_b_silhs` requires unported Arctic
nucleation. Those native cases remain excluded pending these ports.

## Necessary JAX adaptations

- Array, category and starting-level indices become zero based. Mixture
  component labels remain 1 and 2. Statistics retain source one-based
  starting-level values.
- Input/output and output arguments become functional returns. Native random
  routines take explicit keys. `rand_uniform_real` also accepts a draw shape.
- Source `DO` loops with evolving state use `lax.scan`/`fori_loop`; independent
  sample/column algebra is batched. Source routine order and input argument
  names/order are audited against the Fortran files.
- Source fatal returns use `lax.cond`, preserving the boundaries before mixed
  moments and between sampling and clipping. Traced sampling assertions return
  error status through `ErrInfo` instead of executing `ERROR STOP`. Weight
  limiting reports impossible transfers and the source debug-2 assertion.
  Kessler diagnostics validate fractions and empty samples, accumulate in
  source component/sample order, and return failure through both microphysics
  interfaces. Host stops retain the source debug gates; negative-debug tracing
  preserves returned failure status.
- Sample-weighted microphysics diagnostic updates use the shared
  `JaxStats.average_subtimesteps` operation to count one model timestep while
  preserving previously accumulated statistics.
- Host initialization appends the ordered microphysics/SILHS settings and
  approximate physical-space correlation matrices to `<prefix>_setup.txt` at
  debug level 1 or higher. The source fixed-variance-ratio condition gates the
  correlation report. Python streams/scalar formatting replace Fortran units
  and `write_text`; source report names and ordering are retained.
- Timers/OpenACC regions have no JAX numerical counterpart. The source's
  compile-time-disabled 2D sample output stays inactive; its initialization
  routine uses the existing host writer. SILHS radiation remains gated.

No Fortran/F2PY calls are used by these kernels. MT95 is excluded from the port.

## Validation and remaining scope

The focused checks cover native LHS stratification and reproducibility,
permutation reuse, category allocations/weights, overlap against a scalar
reference, inverse-normal transforms against SciPy, Cholesky multiplication,
variance/covariance formulas, category RMS, statistics averaging and fatal
boundaries. Driver checks exercise three real timesteps for interactive and
non-interactive KK (`rico_silhs`) and Morrison (`lba`), two columns with a
three-timestep sequence, and the complete statistics registry. Non-interactive
runs are checked against disabled runs for unchanged prognostic trajectories.

Native random streams intentionally differ from Fortran. Deterministic bindiff
comparisons use the common repeating inputs and durations documented above;
they do not validate native random-stream identity. Short driver checks also
exercise both straight-MC overlap branches. Whole-driver JIT tracing succeeds
at debug=-1 with statistics disabled. Reverse-mode tracing succeeds with the
README's tau-based mixing-length setting; default parcel loops still prevent
reverse-mode differentiation. This is tracing evidence, not validation of
long-run gradients, every microphysics configuration or GPU execution.
Existing Morrison restrictions (including
arctic nucleation and predicted-Nc aerosol activation) still apply; the stock
`mpace_b_silhs` namelist encounters those pre-existing gates. Its unsupported
features are not silently disabled.

One source behavior is preserved deliberately: non-interactive local KK
replaces the mean tendencies but retains sampled variance/covariance tendencies
when `l_var_covar_src=.true.`. The no-feedback KK test uses the case's disabled
variance-source setting. Non-interactive Morrison explicitly clears those
sampled variance sources, following its Fortran dispatch branch.

Re-run the SILHS checks from the repository root:

```bash
.venv-jax/bin/python -m pytest -q \
  clubb_jax/tests/test_silhs_sampling.py \
  clubb_jax/tests/test_silhs_diagnostics.py \
  clubb_jax/tests/test_silhs_fortran_oracle.py \
  clubb_jax/tests/test_silhs_driver.py \
  clubb_jax/tests/test_silhs_port_structure.py
```

The ordinary comparison harness remains applicable to disabled-SILHS runs:
`./tests/run_jax_vs_fortran_cases.py -cases bomex rico lba -max_iters 5 -workers 1`.
Keep its existing
standard-statistics selection and case-specific thresholds when reproducing
this short regression check.

Final audit checks (2026-10-04): 60 sampling, diagnostic, live Fortran-oracle and
host-error tests pass; 74 real driver, source-contract, surface, interface and
diagnostic tests pass; 14 host/dispatch/report checks pass. These counts overlap.
All 20 ordinary comparison cases also pass five timesteps with four C8 columns,
standard statistics and their existing tolerances. Reports are in the ignored
`output/silhs_audit_fixes/` directory. GPU execution, active frozen-phase Morrison
and long-run gradients remain unvalidated.
