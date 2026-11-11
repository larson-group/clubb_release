# Deliberate BFB adaptations

This records C++ arithmetic choices made specifically to reproduce the compiled
Fortran result. Review these choices before changing compilers, optimization
flags or math backends. They preserve the current reference's rounding; they
are not claims that one spelling is scientifically more accurate or faster.

The initial evidence configuration is NVHPC 26.9, Debug `-O0`, double precision,
CPU Fortran versus Kokkos 5.2.2/CUDA 13.3 on an RTX 3080 (`sm_86`). The regular
acceptance gate is now CPU Fortran versus a separate Serial Kokkos build at
zero tolerance. GPU comparisons are occasional milestone checks; document
the newly integrated piece where results change, even slightly. Do not
investigate GPU numerical causes unless requested. GPU BFB is not required
before the next port. Existing CUDA adaptations remain reviewable here;
their CPU necessity has not yet been established. Fortran `-O0` still contracts some
expressions. Explicit `Kokkos::fma(a,b,c)` selects a single rounded multiply-add
rather than leaving contraction and the selected product to the C++ compiler.
No comparison threshold or fast-math flag was changed for these adaptations.

Keep the entry identifiers stable. For each new adaptation, record the routine
and expression, the ordinary translation it replaces, the observed cause,
and before/after evidence. Distinguish a proposed fix from a verified one.
Ordinary index/layout conversion and direct translations of integer powers
need no entry unless a particular evaluation order was selected to fix bits.

## Optional SILHS factorization compatibility

`KOKKOS_BFB` is a local [CMake option](CMakeLists.txt), default `OFF`. The default
[Cholesky factorization](src/CLUBB_core/matrix_operations.hpp) uses the JAX-style
column algorithm and preserves Native scaling, lower-triangle ownership and ten
total attempts. Its arithmetic order deliberately differs from bundled LAPACK;
SILHS rounding changes are accepted for this mode. Failure-only eigenvalues use
cyclic Jacobi rotations. No global comparison tolerances are changed.

Configure with `-DKOKKOS_BFB=ON` for the existing Fortran/f2k scaling,
factorization and eigenvalue path used by strict CPU comparisons. This selection
preserves that library implementation rather than forcing the custom algorithm
to imitate its rounding. If the backend is unavailable, CMake warns and reports
custom operations as the effective choice. The comparison runner records both
requested and effective modes. GPU BFB is still not promised; solver LAPACK
selectors and condition-number diagnostics remain unsupported in both modes.

## Recorded choices

| ID | C++ location | Deliberate choice and reason | Evidence and verification |
| --- | --- | --- | --- |
| BFB-001 | [Skx_func](src/CLUBB_core/Skx_module.cpp) | Evaluate native `sqrt(xp2+tol)**(-3)` as `1/(s*(s*s))`, then multiply by `xp3`, rather than a real-exponent `pow(s,-3)`. The latter differed by one bit. | NVHPC integer-power runtime disassembly showed square, multiply, reciprocal. The corrected direct native/GPU PDF-kernel check passes; the pre-fix executable was retained locally in `/tmp/clubb_core_pow_baseline`. |
| BFB-002 | [linear_interpolated_azm_2D and linear_interpolated_azt_2D](src/CLUBB_core/grid_class.cpp) | Explicitly round the above-level weighted product and fuse the below-level product with the sum. CUDA's implicit contraction can select the other product. Includes extrapolated momentum boundaries. | Native assembly and failing direct comparisons identified the contraction. Both grid directions and model comparisons pass exactly after the change. Local evidence: `/tmp/clubb_core_lower_validation/grid_native_azt_assembly.txt` and `/tmp/clubb_grid_fma_baseline`. |
| BFB-003 | [linear_interpolated_azm_2D optional floor](src/CLUBB_core/grid_class.cpp) | Use `value > minimum ? value : minimum` instead of a generic maximum that keeps the first argument on ties. This retains the reference's second argument, including opposite signed zeros. | A signed-zero direct native comparison failed before the change and passes afterward. Both grid directions and floor cases are covered by the existing grid check. Local pre-fix log: `/tmp/clubb_core_lower_validation/grid_operators_before_zero_tie.log`. |
| BFB-004 | [mono_cubic_interp](src/CLUBB_core/interpolation.hpp) | Explicit FMA for weighted slope means (`p00`, `pp1`), cubic coefficients (`c1`, `c2`), nested cubic evaluation and linear extrapolation. Retain the precise fused product and nested order rather than relying on the compiler to contract the written sum. | Broader cubic grid stencils exposed differences. Native assembly established each contraction. Existing scalar interpolation and linear/cubic grid checks pass exactly. Local evidence: `/tmp/clubb_core_lower_validation/interpolation_native_cubic_assembly.txt`. |
| BFB-005 | [mono_cubic_interp quintic branch](src/CLUBB_core/interpolation.hpp) | Explicit FMA in `beta`, `alpha` and nested quintic evaluation, preserving the native grouping and intermediate `zn*zn`. This branch is compile-time disabled in the normal standalone build. | An isolated build of the actual native quintic-enabled body and corresponding C++ branch passes exact comparisons. This is separate branch evidence, not coverage from ordinary model runs. Local assembly/results: `/tmp/clubb_core_lower_validation/interpolation_native_quintic_assembly.txt` and `interpolation_quintic/`. |
| BFB-006 | [pdf_closure_driver: Skw_velocity](src/CLUBB_core/pdf_closure_module.cpp) | Spell the native reciprocal-times-quotient as one division: `wp3_zm / (wp2_max*(1-sigma_sqd_w))`. A literal pair of divisions rounds differently. | Native assembly shows a multiplied denominator and one division. Before: only `Skw_velocity` differed in BOMEX/COBRA, with `chi` also differing in DYCOMS2_RF02_DO. After this and BFB-007: all three five-step, 30 s, all-stats comparisons are byte-identical. Local assembly: `/tmp/clubb_pdf_parent_native_assembly.txt`. |
| BFB-007 | [pdf_closure_driver: sampled chi; compute_cloud_cover: chi_mean](src/CLUBB_core/pdf_closure_module.cpp) | Round `mixt_frac*chi_1`, then evaluate `fma(1-mixt_frac,chi_2,first_product)`. The ordinary weighted sum can fuse the other component. | Native parent/cloud-cover assembly establishes this order. The targeted three-case all-stats rerun described in BFB-006 is byte-identical. Local evidence: `/tmp/clubb_cloud_cover_native_assembly.txt` and `/tmp/clubb_core_lower_validation/pdf_parent_after_rounding_fix/`. |
| BFB-008 | [calc_trapezoid_zt and calc_trapezoid_zm](src/CLUBB_core/pdf_closure_module.cpp) | Explicitly round the first (upper) trapezoid term, then fuse the second (lower) term's final multiplication by inverse grid spacing. Other products retain native order. | Short multicolumn BOMEX/COBRA runs with vertical averaging differed; two-call-only and cloud-cover-only paths matched. Native assembly shows the lower term in the final `vfmadd213sd`. An initial reading fused the opposite term and did not resolve the difference; independent review corrected the operand mapping. After correction, both cases with default, two-call-only and averaging/cloud-cover settings are byte-identical (six five-step, 30 s pairs). Local assembly: `/tmp/clubb_trapezoid_{zt,zm}_native_assembly.txt`; results: `/tmp/clubb_core_lower_validation/pdf_parent_flags_final/`. |
| BFB-009 | [pdf_closure: rcm/rcp2 floors and Luhar diagnostic denominators](src/CLUBB_core/pdf_closure_module.cpp) | Explicit ternaries select the second argument on equality for native `max(zero_threshold,rcm)`, `max(zero_threshold,rcp2)`, `max(wm,eps)` and `max(rtm,eps)`. This extends the signed-zero convention established in BFB-003. | The native tie behavior is verified by the earlier grid check; this is a source-wide consistency choice, not a newly observed PDF mismatch. Integrated CPU closure checks pass all19 valid short cases and full360-step BOMEX standard stats; the Luhar branch is not exercised by that matrix. |
| BFB-010 | [transform_pdf_chi_eta_component: corr_rt_thl_term](src/CLUBB_core/pdf_closure_module.cpp) | On CUDA, use `__dmul_rn` for the correlation term's final multiplication so it rounds before the following chi/eta variance sums. Other backends retain the source product. This local CUDA intrinsic avoids a volatile memory spill. | All sampled transform inputs match in the first-step DYCOMS comparison, but six `stdev_chi_1` values differ. Native assembly rounds the product separately; CUDA SASS fuses it across the assignment. Recomputing a failing cell reproduces the native bits with the rounded product and C++ bits with the fused expression. Corrected targeted GPU transform diagnostics are byte-identical. Integrated CPU checks pass19 valid short cases and full360-step BOMEX standard stats. Local evidence: `/tmp/clubb_pdf_closure_native_assembly/transform_pdf_chi_eta_component.txt` and `/tmp/clubb_pdf_transform_cuda_2.sass`. |
| BFB-011 | [pdf_closure: cloud_frac and rcm mixtures](src/CLUBB_core/pdf_closure_module.cpp) | Round the first component product, then explicitly fuse the second component product with it. This is the same contraction pattern as BFB-007, now in the physical cloud/water means. | First-step DYCOMS has 39 differing `rcm` cells and two differing `cloud_frac` cells whose sampled component inputs are identical. Recomputing those cells with each FMA choice reproduces the Native/C++ bits respectively. Native assembly at source1016/1018 confirms the second-product contraction. Corrected targeted GPU cloud/water diagnostics match; CPU checks pass19 valid short cases and full360-step BOMEX standard stats. Local evidence: `/tmp/clubb_pdf_closure_native_assembly/pdf_closure.txt` and `/tmp/clubb_core_lower_validation/pdf_closure_first/mixture_diagnostics/`. |
| BFB-012 | [pdf_closure: six water-moment mixtures](src/CLUBB_core/pdf_closure_module.cpp) | For `wprcp`, `wp2rcp`, `rtprcp`, `thlprcp`, `uprcp` and `vprcp`, round component1's weighted product and explicitly fuse component2 with it. | Native assembly at source1126/1129/1132/1135/1138/1141 has the same second-product FMA for all six. First-step DYCOMS `wp2rcp` recomputation reproduces the Native/C++ bits with the second/first-product FMA respectively. Targeted GPU `wprcp`/`wp2rcp` now match; some other covariance diagnostics still differ. CPU checks pass19 valid short cases and full360-step BOMEX standard stats. Local evidence: `/tmp/clubb_pdf_closure_native_assembly/pdf_closure.txt` and `/tmp/clubb_pdf_moment_diagnostics.log`. |
| BFB-013 | [calc_wp2xp2_pdf: each component's final variance bracket](src/CLUBB_core/pdf_closure_module.cpp) | In `dx**2 + (1+2*corr**2)*var_x`, round the variance product, then explicitly fuse `dx*dx` with that product. Ordinary C++ contracts the variance product into the rounded square instead. | CPU Fortran versus CPU Kokkos full BOMEX differs only in `wp2up2`/`wp2vp2`, starting at the 120 s sample. Independent assembly review identifies opposite contractions in both components: Native `3a1b0–3a1b4` / `3a5d0–3a5d4`, C++ `112d9c–112da0` / `11328b–11328f`. After correction, three-step standard-stats BOMEX and the full360-step/native60s standard-stats case are byte-identical. This is a demonstrated CPU/CPU arithmetic adaptation; the final full 19-case native-duration/standard-stats CPU matrix also passes. All-stats evidence remains scoped to the recorded targeted checks. Local assembly: `/tmp/clubb_pdf_closure_native_assembly/calc_wp2xp2_pdf.txt`, `/tmp/clubb_wp2xp2_cpp_assembly.txt`; pre-fix report: `/tmp/clubb_core_lower_validation/pdf_closure_cpu/full_bomex_before_fourth_moment_fix/`. |
| BFB-014 | [ADG1_ADG2_responder_params: width_factor_1](src/CLUBB_core/adg1_adg2_3d_luhar_pdf.cpp) | Round `two*mixt_frac*(1-two_thirds*beta)`, then fuse `two_thirds*beta` with that term. Ordinary C++ fuses the second term's last product with a rounded first term. | Before correction, short CPU BOMEX/COBRA match but DYCOMS/JUN25 differ. First-step PDF diagnostics show matching mixture fractions and water variances, but differing temperature-component variances. Native assembly `11356–11367` rounds the second term and fuses the first; C++ `35512–3554d` chooses the other product. After correction, all four initial five-step CPU pairs match, both previously failing first-step diagnostic selections match, and all19 native-timestep five-step multi_col_stats CPU pairs are byte-identical. Local assembly: `/tmp/clubb_adg1_assembly/`; pre-fix summaries: `/tmp/clubb_core_lower_validation/adg1_cpu/before_fix/`. |
| BFB-015 | [calc_sfc_varnce: surface temperature/moisture/scalar power bases](src/CLUBB_core/sfc_varnce_module.cpp) | Materialize each `flux/uf` quotient before squaring in both convective-depth branches and passive-scalar variance. Ordinary duplicated C++ quotients were compiled as `flux^2/uf^2`; Native `**2` squares rounded quotients. | NVHPC26.9 CPU assembly proves the different divisions/multiplications at `-O0`. The first sampled COBRA temperature/moisture and DYCOMS moisture surface budgets differ before the change; other surface velocity budgets and saved flux inputs match. After materializing the quotients, all four5step30s multicolumn CPU cases are byte-identical14.4s. The same quotient boundary also fixes passive-scalar variance: COBRA/DYCOMS5step30s all-stats checks at debug_level2 are byte-identical2/2(10.2s), including the actual Native surface checker and C++ length/input checks. The combined full 19-case native-duration/standard-stats CPU matrix also passes. Alternate-depth runtime evidence remains outside this claim. Evidence: `/tmp/clubb_sfc_varnce_{native,cpp}_assembly.txt` and `/tmp/clubb_core_lower_validation/sfc_varnce_cpu/`. The nondefault branch uses the same source power-base rule but requires its own runtime check. |
| BFB-016 | [xp2_xpyp_lhs: main diagonal](src/CLUBB_core/advance_xp2_xpyp_module.cpp) | Use `fma(gamma,lhs_ta,lhs_dp1)+(lhs_diff+lhs_ma)`. Ordinary C++ contracts the product into the diffusion/advection pair before adding dissipation. | First changed shared LHS has identical input arrays. Actual-input replay reproduces both compiled sides at all1376 interior points, including47 scalar and36 wind differences. Clean CPU4/4 short cases and19/19 native-timestep five-step MC cases are byte-identical after correction; restored Native sources and final binaries are trace-free. Evidence: `build/parallel_cpu_lanes/mixing/variance/BFB-016.json`. The combined full 19-case native-duration/standard-stats CPU matrix also passes; untested options remain outside this claim. |
| BFB-017 | [xm_wpxp_lhs: momentum main diagonal](src/CLUBB_core/advance_xm_wpxp_module.cpp) | Round the mean-advection/diffusion pair, then add it to `fma(gamma, ta+pr1, ac)`. Ordinary C++ contracts the product with `(ma+diff)+ac`; Native contracts it with `ac` before adding `ma+diff`. | Actual first BOMEX call has identical five operands and C6/C7/Skw/parameters. Both associations reproduce their respective compiled sides at all 344 interior cells; 56 diagonal cells differ. Native assembly confirms FMA then the separately rounded pair. After the explicit association, clean Serial CPU checks are byte-identical 4/4 five-step30s cases and 19/19 five-step native-timestep multicolumn cases. Source restored exactly and both executables checked free of temporary tracing. Evidence: `build/parallel_cpu_lanes/xp3/lane_artifacts/xm/BFB017.json`, `operand_replay.json` and `native_assembly.txt`. No Gaussian, solver or tolerance change. The combined full 19-case native-duration/standard-stats CPU matrix also passes; untested options remain outside this claim. |
| BFB-018 | [calc_comp_mu_sigma_hm: both emergency-variance denominators](src/CLUBB_core/setup_clubb_pdf_params.cpp) | Use `fma((a*fp1)*(1+zeta), mu1*mu1, ((1-a)*fp2)*(mu2*mu2))`. Ordinary C++ contracts the right product into the rounded left product; Native contracts the left product into the rounded right product. Preserve the common numerator and both surrounding branches. | Actual identical-input replay reproduces all four changed cells, including both original variance ratios and physical sigmas. Original Native/C++ assembly independently confirms opposite contractions. Corrected clean isolated CPU checks pass 4/4 short and 19/19 native-timestep MC pairs, cloudy DYCOMS all-stats and four RICO/LBA/SILHS all-stats pairs byte-identically. Both restored models were force-recompiled and verified free of tracing; combined live-tree 19-case MC/debug2 and four-case all-stats/debug2 gates also pass. Evidence: `build/parallel_cpu_lanes/mixing/metadata_pdf_batch/BFB-018.json` and `diagnosis/emergency_ratio_same_input_replay.json`. |
| BFB-019 | [calculate_mixture_fraction](src/CLUBB_core/new_hybrid_pdf.hpp) | Fuse `Skw*Skw` with the neighboring six-F numerator term and thirty-six-F radicand term before adding the earlier F terms. In the radicand's earlier pair use `fma(1-F_w,12*F_w^2,4*F_w^3)`. Ordinary C++ contracts Skw with the whole prefix, and shortening that prefix can choose the opposite F product. | Native303–315 and CPU assembly establish all three groups. First BOMEX actual-input replay reproduces both original model sides348/348(56 mixture differences); COBRA's first two calls reproduce both sides800/800 and expose two more changed mixtures requiring the inner F pair. Final candidate agrees with Native800/800. Clean combined model evidence is recorded in `build/parallel_cpu_lanes/xp3/PDF_ALGORITHM_LANE_RESULT.json`; detailed replay/assembly is in `lane_artifacts/pdf_algorithms/diagnosis/`. |
| BFB-020 | [calculate_coef_wp4_implicit](src/CLUBB_core/new_hybrid_pdf.hpp) | Retain Native's separately rounded first-component pair, nested FMA groups for the following three terms, and final sixth fraction: `(last_fraction+component_pair_2)+component_pair_1`. Ordinary C++ nests the first pair into the later contraction chain. | Native815–820 and CPU assembly identify the grouping. Actual first BOMEX input replay reproduces both compiled model sides348/348, with134 changed wp4 coefficients; candidate agrees with Native348/348. The same final clean batch gates cover the correction. Evidence: `diagnosis/wp4_same_input_replay.json`, Native/C++ wp4 disassembly and the lane result. |
| BFB-021 | [calc_responder_driver: both positive-factor bounds](src/CLUBB_core/new_hybrid_pdf_main.cpp) | Materialize `(1+mixt_frac)/sqrt(mixt_frac*(1-mixt_frac))` before multiplying the cubed correlation and dividing by `F_w^three_halves`. Ordinary C++ merges the quotient into one denominator. Retain ordinary `Kokkos::pow` and all source clipping branches. | Extra COBRA hybrid/predictor=false all-stats exposes10 changed Skv clips among400 cells, followed by component-v variance and moment changes. Actual Native/C++ private-object replay reproduces their respective five model outputs400/400; the rounded-quotient candidate matches all seven Native outputs400/400. A power rewrite alone was insufficient and was not applied. Evidence: `diagnosis/bounds_actual_objects_replay.json` and both private-driver disassemblies. Clean combined gates and tested options are in the lane result. |
| BFB-022 | [pdf_hydromet_microphys_prep: SILHS seed product](src/Microphys/pdf_hydromet_microphys_wrapper.hpp) | Multiply the seed and iteration as `uint32_t`, then `bit_cast<int32_t>` to retain the Native signed 32-bit seed bits. Ordinary signed C++ multiplication can overflow with undefined behavior on long runs. This is an integer language boundary, with no floating-point workaround. | Native declares `lh_seed_custom` as `genrand_intg` and assigns `lh_seed*itime`. Integrated stochastic RICO passes five-step all-stats raw-file equality with importance sampling and random starting levels enabled. This checks that stream scope; it does not establish all overflow inputs or other Fortran compilers. |
| BFB-023 | [bivar_LL_mean](src/Microphys/KK_microphys/PDF_integrals_means.hpp) | Group the variance/cross terms with explicit FMAs, contract the mean pair separately, then add before `exp`. Ordinary C++ contracts the entire written sum in another order. | NVHPC26.9 CPU/O0 actual-input replay matches both Native mixture components and the weighted radius; original C++ replay reproduces its different bits. Native assembly confirms the grouping. Combined023–026 DO/DS five-step MC files match2/2; the final six-case native-duration/dt standard gate also matches6/6. [Replay record](drafts/driver_microphysics/BFB-023.json). |
| BFB-024 | [quadrivar_NNLL_covar](src/Microphys/KK_microphys/PDF_integrals_covar.hpp) | Preserve separate variance/cross and mean groups in the first exponent, plus Native fused correlation complements/difference. The second exponent and all other factors are unchanged. | Actual-input replay reproduces both original implementations and matches Native2/2 components after the change; assembly identifies each contraction. Same combined MC gate; the separate early evaporation-covariance difference was subsequently corrected by BFB-028, and the live DO/DS five-step all-stats files match. [Replay record](drafts/driver_microphysics/BFB-024.json). |
| BFB-025 | [bivar_LL_covar_partial](src/Microphys/KK_microphys/KK_upscaled_turbulent_sed.hpp) | Form `alpha²+2*alpha` with Native's FMA and keep variance/cross terms separate from the contracted mean pair before `exp`. | All35 differing actual inputs reproduce both original sides and match Native after the change, including31 swapped number-concentration calls. Native assembly confirms the groups; same combined MC gate. [Replay record](drafts/driver_microphysics/BFB-025.json). |
| BFB-026 | [trivar_NLL_mean_const_x1](src/Microphys/KK_microphys/PDF_integrals_means.hpp) | In the constant-chi branch, separately contract variance/cross and mean groups, then add before `exp`; retain `pow`, thresholds and weights. | Eight actual calls reproduce their Native/C++ originals and all eight match Native after grouping. Native assembly confirms the separate contractions; same combined MC gate. [Replay record](drafts/driver_microphysics/BFB-026.json). |
| BFB-027 | [cos_solar_zen: solar declination](src/Radiation/cos_solar_zen_module.cpp) | Round the five sine/higher-harmonic terms separately from `c0+c1*cos(t)`, then add the groups. Ordinary C++ chains all seven terms; parentheses alone select the Native NVHPC contractions. | Actual first changed solar input reproduces both implementations: grouped zenith cosine matches Native exactly. Reconstructed baseline reproduces the fixed-SST prefix; grouped candidate passes all540 standard-statistics steps byte-for-byte against both Native references. Only fixed-SST and JUN25 call this function in the default matrix. [Evidence](drafts/driver_radiation_silhs/Radiation/BFB-027.json). |
| BFB-028 | [quadrivar_NNLL_covar constant branches](src/Microphys/KK_microphys/PDF_integrals_covar.hpp) | Group variance/cross terms separately from the mean pair in the `const_x2`, `const_x1` and `const_x1x2` exponentials. Preserve prefactors and final subtractions. | Native assembly confirms all three groups; eight differing actual constant-x2 calls reproduce both originals and match Native after the change. Combined with the independent strided hole-filling correction, DO/DS/RICO20-step standard files match3/3. [Evidence](drafts/driver_microphysics/BFB-028.json). |
| BFB-029 | [est_kessler_microphys_api: AKstd_cld](src/SILHS/est_kessler_microphys_module.cpp) | Materialize the rounded `AKm/cloud_frac` quotient before squaring. Repeated View expressions otherwise compile as a quotient of squares. Preserve the existing positive-cloud guard and all other formulas. | Actual View-access replay reproduces original and Native/candidate bits; Native/C++ kernel assembly identifies the change. Scalar-only extraction did not reproduce it and is not evidence. RICO_SILHS15-step standard output is byte-identical against both Native references after the fix; initial/trace controls reproduce all635 frozen variables. [Evidence](drafts/driver_radiation_silhs/SILHS/BFB-029.json). |
| BFB-030 | [Morrison rain-slope N0RR](src/Microphys/Morrison_microphys/module_mp_graupel.hpp) | At six supported rain-slope sites, compute the fourth power with two `float` squarings before the original product/division. `Kokkos::pow(float,int)` returns double and promoted the rest of the expression; Native integer exponentiation retains float. | All196 actual engine input/context records match; the first changed NR3D output is one float ULP. Exact operands reproduce both old and Native N0RR/NR values. The clean candidate matches all641 saved variables through49 steps; the combined030/032/033 full360 files match. Four similar sites belong to the explicitly excluded Seifert–Beheng branch and remain unexercised. [Evidence](drafts/driver_microphysics/BFB-030.json). |
| BFB-031 | [Parabolic expaelem: phi(2)](src/Microphys/KK_microphys/Parabolic.hpp) | Contract the low-power and high-power polynomial groups separately, then add. Preserve all coefficients, normalization, recurrences, gamma and surrounding PDF factors. | Native assembly and actual-input probes isolate the first difference to phi(2), then sphip/U/Dv. Candidate matches four changed evaporation means, two covariances and the independent late-accretion parent/six children, including U/V derivatives. RICO100-step standard file and all six final full-duration KK/SILHS files are byte-identical. [Evidence](drafts/driver_microphysics/BFB-031.json). |
| BFB-032 | [Morrison hl_before/hl_after](src/Microphys/morrison_microphys_module.hpp) | Round `grav*z` before fusing `Cp*T` into each energy sum. Ordinary C++ rounds `Cp*T` and fuses the other product. Leave the residual formula and all remaining terms unchanged. | Actual View replay reproduces all508 original C++ points and matches all508 Native points after the change. Native object constants/instructions identify the order. Morrison100-step standard output is byte-identical; the combined030/032/033 full360 run also matches. [Evidence](drafts/driver_microphysics/BFB-032.json). |
| BFB-033 | [Morrison warm rain conservation limiter](src/Microphys/Morrison_microphys/module_mp_graupel.hpp) | Balance the seven-term numerator in Native order: `(PRC+(PRACG+PRA)) + ((QR3D/DT)+PRACS) -PSMLT-PGMLT`. Preserve the denominator and other limiters. | All660 engine input/context prefixes match. Exact warm-point operands reproduce the old and Native final PRE bits; Native object/DWARF establishes the groups. Combined030/032/033 files match through165 and the full360 native steps, including the rebuilt live executable. [Evidence](drafts/driver_microphysics/BFB-033.json). |

## Documented CPU/GPU math-library differences

An isolated probe using actual first-step JUN25 inputs found these results in
the recorded configuration. The user chose CPU BFB as the regular gate, with
GPU differences recorded at milestones. These results do not change the
model comparison's zero default limits or excuse unrelated numerical changes.

| Function and exact input | NVHPC CPU result | CUDA result | 110-digit calculation rounded to binary64 |
| --- | --- | --- | --- |
| `erf(-0x1.0470794f5ed90p-2)` | `-0x1.1fa8abdc20552p-2` | `-0x1.1fa8abdc20553p-2` | `-0x1.1fa8abdc20552p-2` |
| `exp(-0x1.16ecb6c757b1ep+0)` | `0x1.5870daa09d50dp-2` | `0x1.5870daa09d50ep-2` | `0x1.5870daa09d50ep-2` |

CPU symbols are `__mth_i_derf` and `__fd_exp_1_avx2` from the reference
executable's `libnvcpumath.so`; CUDA uses the ordinary Kokkos scalar functions.
The probe source/build/output is under `/tmp/clubb_pdf_math_probe/`. The
high-precision check uses Python Decimal with 110 digits, the convergent erf
power series and Decimal's exp; it is diagnostic evidence, not a substitute
for the maintained model test. Keep the ordinary GPU math implementation;
no CPU emulation or additional math implementation is required for these
documented device/library differences.

### ADG1 and PDF-utility GPU checkpoint

With the same NVHPC/Kokkos/CUDA configuration, five steps at30 s and
`multi_col_stats` are byte-identical for BOMEX and COBRA. DYCOMS2_RF02_DO and
JUN25_ALTOCU complete but differ in five saved fields each. Their largest
cloud-fraction absolute difference is `1.1102230246251565e-16`; remaining fields
and exact difference counts are preserved in
`/tmp/clubb_core_lower_validation/pdf_utilities_gpu/field_differences.json`.
The corresponding CPU checks pass all19 short native-timestep cases and four
five-step all-stats cases at zero tolerance.

The isolated erf/exp discrepancies above establish actual library differences;
this checkpoint does not attribute every saved-field difference to those two
functions. No further per-field GPU tracing is planned under the accepted CPU gate. No GPU tolerance, special math implementation or CPU execution
fallback was introduced to make these comparisons pass.

### PDF moments, saturation and cloudy-updraft checkpoint

The following integration group added four public PDF moment routines, all ten
saturation routines and the cloudy-updraft helper. The next GPU milestone
completed all four short multi-column cases, with differences now present in
BOMEX and COBRA as well as DYCOMS2_RF02_DO and JUN25_ALTOCU. The first saved
differences occur at step4 in BOMEX and step2 in the other three cases.
The group is recorded together; no finer introduction boundary was tested.
Reports are under `/tmp/clubb_core_lower_validation/pdf_supported_final/`.
CPU checks of this group pass all19 short cases and full360-step BOMEX exactly.

Earlier diagnostic evidence is retained in
`/tmp/clubb_pdf_gpu_milestone_diagnosis.json`; no GPU arithmetic adaptation was
added for this checkpoint. Following the user's updated direction, subsequent
GPU checks only record the newly integrated piece and changed saved results.
No further GPU numerical investigation is required to advance the CPU port.


### Moment parents, lower physics and metadata preparation checkpoint

The next group integrates the helper/pressure/mixing work, moment-advancement
parents and lower physics, numerical checks, type lifecycle, metadata/PDF
preparation and calendar. All four representative GPU case pairs complete five
30-second steps with four C8 columns and multi-column stats. All four saved GPU
files differ from the previous PDF/saturation checkpoint; the four Native
reference files remain byte-identical before and after the group. Both GPU
checkpoints differ from CPU Fortran in all four cases. The group is recorded
together; no finer introduction boundary or numerical cause was investigated.

Reproducible commands, frozen binary/source hashes, saved outputs and the
file comparison are under `build/validation/gpu_physics_metadata_checkpoint/`
and `build/checkpoints/physics_metadata_gpu/`. The combined CPU gate passes all
19 short native-timestep cases and four cloudy/scalar/SILHS all-stats cases.

### Final PDF and initialization checkpoint

The remaining PDF algorithms, grid/parameter/flag initialization, fixed dycore
constructor and final call hookups compile in the complete CUDA executable.
The same four default-ADG1 cases complete five 30-second steps with four C8
columns and multi-column stats. All four saved GPU files and all four Native
reference files are byte-identical to the preceding physics/metadata checkpoint.
This group introduces no additional saved-output changes in that scope.
GPU versus CPU Fortran still differs in all four cases; no numerical causes
were investigated. PDF7 has active CPU comparison evidence; this GPU checkpoint
does not exercise it or the input-fields-gated PDF selectors.

The frozen binaries, runtime/input hashes, commands and file comparison are
under `build/parallel_gpu_lanes/final_core_checkpoint/`, particularly
`BINARY_FREEZE.json` and `GROUP_OUTPUT_CHANGES.json`. Remapping's retained Native
outer-package callers do not exercise its C++ bodies in this checkpoint.

### Active Hybrid7 GPU checkpoint

A separate COBRA check exercises default ADG1 and active Hybrid7 before and
after the PDF/initialization integration. Each of the four model pairs completes
five 30-second steps with four C8 columns and multi-column stats. Generated
namelists confirm selectors 1 and 7 on both sides. The two stages use identical
Native reference binaries and inputs.

The Hybrid7 GPU saved file changes with this group; the default ADG1 GPU file
and both Native reference files remain byte-identical. This records the change
at the active hybrid integration boundary without investigating numerical values
or causes. Both stages retain the maintained comparator's exit 1 for their two
GPU/CPU differences. This is GPU runtime evidence, not GPU/CPU BFB or performance
evidence; input-fields-gated selectors remain unexercised. All unmasked numeric
values in the eight saved files are finite.

Frozen binaries, source/input hashes, actual selectors, return codes and file
comparisons are under `build/parallel_gpu_lanes/final_hybrid_checkpoint/`,
including `BINARY_FREEZE.json`, `RUN_RESULTS.json` and `GROUP_OUTPUT_CHANGES.json`.

## Core-phase CPU integration evidence (`ac37a3898`)

The final Serial CPU Kokkos and CPU Fortran comparison passes all 19 cases at
native duration/timestep with standard stats and zero comparison limits. Both
models complete the independently expected 17,298 main steps; all 19 saved
NetCDF file pairs are byte-identical, with normal exits and unchanged frozen
binaries. The configuration is NVHPC 26.9, real8 and Debug flags. Each case's
configured PDF path is covered; active Hybrid7 evidence remains scoped to its
separate short CPU/GPU checks. Source-only and untested branches remain explicit
in the [core inventory](CORE_CONVERSION_PLAN.md).

The combined short checks, targeted all-stats/predictor variants, fixed dycore
cases and reset/reinitialization sequence also pass. Results, expected counts
and frozen source/compiler/binary records are under
`build/validation/final_supported_core_cpu/` and
`build/checkpoints/final_supported_core_cpu/`.

## Standalone driver CPU evidence

The later complete driver passes all 19 native-duration/timestep standard-stats
cases in compatible frozen Serial CPU lanes: 17,298 steps per backend, 38 normal
exits and byte-identical saved files. The rebuilt live Morrison LBA executable
also matches all 360 native steps after BFB-030/032/033. Short all-stats and
lifecycle/IOP/restart/fixed-dycore checks supplement the full matrix. The
[driver plan](DRIVER_CONVERSION_PLAN.md) records settings, scope, source reviews,
retained Native services and GPU checkpoints. Earlier failed controls and each
arithmetic adaptation remain preserved; no comparison tolerance was increased.

## Diagnostic header formatting

The host-only `f6_2_field` conversion in
[err_info_type_module.cpp](src/CLUBB_core/err_info_type_module.cpp) deliberately
matches Native F6.2 output rather than ordinary C++ `std::fixed` formatting.
On NVHPC 26.9, the exact tie `0.125` prints `0.13`; the ordinary C++ stream prints
`0.12`. Conversely, binary `0.015` prints `0.01`, so multiplying in core precision
and rounding the resulting apparent tie would incorrectly print `0.02`.

The conversion scales the exact binary significand with integer arithmetic and
rounds hundredths away from zero at exact ties. It preserves signed zero and uses Native nonfinite spelling and
all-star overflow. This affects diagnostic character records only; model fields,
coordinate normalization and scientific arithmetic are unchanged. Actual Native
controls are recorded under `build/validation/config_error_cpp/`; the focused
public-API check also compares complete global/column records, including padding,
metadata and field overflow. This choice is independent of `KOKKOS_BFB`.

## Review notes

The source expressions remain in [Fortran](../src/CLUBB_core/). Source comments
mark the explicit arithmetic boundaries. Existing tests establish equality only
for their recorded cases and branches; they do not establish BFB for every
compiler or optimized build. Local `/tmp` evidence paths are investigation
artifacts, not prerequisites for building or running the port.

Update this record when an adaptation is added, revised or removed. If a later
compiler makes it unnecessary, retain the reason and replacement evidence here
so the decision can still be reviewed.
