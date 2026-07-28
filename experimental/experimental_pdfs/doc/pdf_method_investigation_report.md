# CLUBB PDF-method investigation

**Status:** 15 July 2026. This report summarizes the completed non-neural PDF
investigation. Lower normalized root-mean-square error (NRMSE) is better; a
positive joint-density skill is better than a moment-matched single Gaussian.
Results from different evidence tiers are kept separate because they answer
different questions.

## Executive conclusion

The main error is not that two Gaussians are mathematically unable to describe
a tilted cloud plume. It is that the CLUBB moments supplied to the closure do
not uniquely determine **whether a rare cloudy tail exists, how much probability
it contains, or how far its vertical-velocity centroid lies from the bulk**.
Exact preservation of the grid mean, full `(w, rt, thl)` covariance, `wp3`, and
positive-semidefinite (PSD) component covariances leaves too little freedom in
the tested implementable 2G diagnoses to place the rare tail reliably.

The strongest current architecture is therefore a small, explicitly diagnosed
third tail component plus an exact residual ADG1 pair. The best current
*hypothesis* is the five-coefficient calibrated prescribed-skew-tail ADG3,
enabled only when a moment-Gaussian cloud-probability guard says a tail is
plausible. In a frozen ten-case profile-statistics sensitivity it reduces the
four direct cloud-transport NRMSE from `0.42360` for safe ADG1 to `0.39667`
(`6.36%`) and essentially ties the eight-diagnostic score (`0.40975` versus
`0.41105`). This is promising but **not yet a validated closure**: the guard
threshold was selected after inspecting a declared sensitivity set, the method
is not casewise dominant, and raw-density fidelity remains inconsistent.
The first larger regime-holdout attempt did not produce a score: every fitted
ADG1 candidate had at least one audited-invalid state, so the deliberately
strict 100%-coverage selector stopped at the first family before reading the
ARM outer fold. This is a realizability/coverage blocker, not evidence for or
against the tail family.

Strict Bergmann 3G is the safer low-risk shape baseline. It is comparatively
easy to map into CLUBB and usually preserves density geometry better, but its
mean-centered third component does not isolate the detached rare-cloud tail as
effectively. No tested 2G method was simultaneously exact, low-parameter,
portable, and competitive on tail transport.

## The physical and mathematical problem

The original LES-advance diagnosis localized the `wp2` budget error to
`wpthvp`, then to `wprcp`: LES-consistent `wpthlp` and `wprtp` did not repair
the buoyancy-flux discrepancy. Direct SAM `WPRCP` and its independent cloudy-
transport field agree to within converter precision, ruling out a unit or
averaging convention error.

The repaired CLUBB New PDF sharpened this diagnosis. Over the ARM 65–75 ks,
0.9–2.3 km layer it predicted the unthresholded `w-chi` covariance almost
exactly (median model/LES `wpchi=1.036`, profile correlation `0.9997`) but
converted only `0.284` of it into `wprcp`, versus `0.658` in SAM. Its `wprcp`
RMSE (`1.877e-5`) was better than ADG1 (`2.463e-5`) but still deficient; the
tested New Hybrid was worse (`2.022e-5`). The modeled cloud-bearing component
was a broad, partly clear population rather than a small saturated updraft.
This established that the missing object is conditional cloud-tail geometry,
not the total covariance itself. New PDF and New Hybrid were subsequently
excluded from the focused Bake-off by design.

The motivating ARM plane (SAM minute 529, 2180 m) contains only 231 cloudy
cells out of 9,216, and just 19 cells (`0.206%` of the domain) carry half of the
net upward `w'rc'`. The raw SAM cloud-water-weighted vertical centroid is
`+2.396 m/s` from the grid mean. A diagnosis using only the supplied mean and
covariance places it at `+0.492 m/s`: a factor of `4.87` too small. The scalar
displacement is only a factor of `1.40` too small, the direction is wrong by
`24 degrees`, and the full covariance-metric mismatch is `3.22 RMS`.

This is an **anisotropic placement error**, not a missing universal distance
multiplier. Scaling the complete saturation-tail vector enough to fix `w`
would over-displace its scalar coordinate.

For a Gaussian mixture, total `w-x` covariance is the sum of within-component
tilt and covariance created by separated component centers. ADG1 places most
of that covariance into center separation because its within-component
`w-scalar` covariances are zero or tightly constrained. A general full-
covariance 2G can tilt, but after imposing the three means, six independent
covariances, and `wp3`, a nominal 19-parameter 2G has nine remaining shape
degrees of freedom that CLUBB does not prognose. A diagnosis must replace
those degrees of freedom with assumptions.

The raw ARM analyzer quantifies the resulting realizability bottleneck. The
raw tail centroid lies inside the interval allowed by `wp2/wp3` center
separation on `63.0%` of diagnosable planes, but inside the interval that also
satisfies the complete covariance and component-PSD constraints on only
`31.7%`. This is not a proof that every 2G is impossible—the cloud-water-
weighted centroid also includes internal tilt—but it explains why forcing a
2G tail center repeatedly exhausted its Schur/PSD capacity.

## What was tested

All operational candidates use only CLUBB-available pressure,
thermodynamics, means, covariance, `wp3`, and, where stated, `rtp3` and
`thlp3`. Raw cloud water and target transports are used only for fitting or
verification. Accepted mixtures reproduce the supplied mean, complete
covariance, and `wp3` and retain PSD component covariances.

The additive registry now contains 28 independently selectable families; new
experiments no longer overwrite prior ones. The main lines of inquiry were:

- **Baselines and capacity bounds:** ADG1, strict Bergmann 3G, and fully
  flexible per-plane 2G/3G oracles. The flexible oracles establish what a
  Gaussian count can represent, not what CLUBB can diagnose.
- **Tilt allocation:** tilted Bergmann, covariance-skew, tilt-boost, and
  tilt-allocation families redistributed `w-(rt,thl)` covariance while
  retaining exact grid moments. PSD limits usually capped the useful transfer.
- **Cloud-aware 2G placement:** saturation-partitioned, cloud-lifecycle,
  transport-background, skew-center, centroid-tilt, and asymmetric-centroid
  ADG2 variants tried to separate a correlated transport population from a
  quiet background. They improved individual planes but did not yield a
  portable global diagnosis.
- **Residual-tail 3G placement:** cloud-anchored, skew-rescue,
  Skw-vertical, anisotropic, centroid-boundary, prescribed-skew-tail, and its
  calibrated variant dedicate one small component to the rare tail and close
  the residual moments with ADG1. This was the most productive family.
- **Causal memory:** previous-below saturation information improved tail
  occurrence detection, but did not robustly locate its vertical centroid.
  A one-minute previous-below signal reduced repeated-split activation Brier
  score (`0.0357` to `0.0260`), while tail-w RMSE was essentially unchanged
  (`1.155` to `1.173`). Fitted memory gains were localized by height and often
  neutral, so memory is useful as a possible gate, not yet as a universal
  displacement closure.
- **Non-Gaussian saturation hinge:** a mode constructed to affect nonlinear
  cloud transport without changing the supplied mean or covariance was
  realizable analytically, but `wp3` could not determine its amplitude. The
  bounded form still underperformed the scalar-skew reference and was not
  promoted into the Bake-off.

## Evidence, kept in tiers

### 1. ARM raw-3D development and same-case holdout

The raw-3D ARM run is the hypothesis-development source. Complete 15-minute
time blocks were assigned to train/validation/test, so individual test results
are temporally held out but the family formulas were still invented while
studying ARM.

The three-control prescribed-skew-tail ADG3 scored `0.26346` test aggregate
NRMSE on the primary split. Across three alternate complete-block splits it
scored `0.22484`, `0.31418`, and `0.56022` (mean `0.36641`) and beat the
centroid-boundary and Skw-vertical alternatives on every split. Its density
skill remained negative on all three (`-0.1057`, `-0.2741`, `-0.1104`).

Allowing two more controls produced calibrated prescribed-skew-tail ADG3. It
improved the primary ARM test score only from `0.26346` to `0.25539` and kept
negative density skill (`-0.06457`). That modest same-case gain did not alone
justify the extra controls; the independent tests below are more informative.

Flexible-oracle studies show capacity but also overfitting danger. A matched
per-level 3G oracle improves the four fitted diagnostics over flexible 2G by
`23.5%`, yet is worse on every withheld cloud diagnostic and drives component
correlations near one. On the exceptional ARM reference plane, a freely fitted
tail 3G can nearly eliminate transport error, but that result uses per-plane
target information and is not an operational diagnosis.

Authoritative artifacts:

- [`output_pdf_bakeoff/global_prescribed_skew_tail_comparison_06f93c9ab073ca15.json`](output_pdf_bakeoff/global_prescribed_skew_tail_comparison_06f93c9ab073ca15.json)
- [`output_pdf_bakeoff/repeated_holdout/repeated_holdout_dc8c7c1f69241cc5.json`](output_pdf_bakeoff/repeated_holdout/repeated_holdout_dc8c7c1f69241cc5.json)
- [`output_pdf_bakeoff/sam_w_tail_skew_study_multifeature.md`](output_pdf_bakeoff/sam_w_tail_skew_study_multifeature.md)

### 2. Frozen ten-case profile-statistics comparison

ARM-selected coefficient vectors were frozen and evaluated without refitting
on 600 target-blind planes from each of ARM, BOMEX, DYCOMS RF01, four DYCOMS
RF02 variants, LBA, RICO, and cloud-free Wangara. Profile files provide the
four direct transports and secondary cloud summaries, but not raw joint-PDF
density. Non-PSD input states were excluded rather than projected.

| Frozen method | Direct-four NRMSE | All-eight NRMSE | Cloud-active direct-four |
|---|---:|---:|---:|
| Safe ADG1 defaults | `0.42360` | `0.41105` | `0.49477` |
| Bergmann 3G, ARM vector | `0.56813` | `0.56624` | `0.55862` |
| Prescribed-tail ADG3, ARM vector | `0.48902` | `0.49655` | `0.47109` |
| Calibrated-tail ADG3, ARM vector | `0.45830` | `0.47124` | `0.43562` |
| Calibrated tail + fixed `p>=1e-4` guard | **`0.39667`** | **`0.40975`** | `0.46717` |

The calibrated family is substantially more portable than the original
prescription and helps cloud-active direct transports, but ungated activation
causes large dry-regime errors. The `1e-4` guard removes Wangara activation and
gives the best equal-case direct-four score; it improves cloud-active
direct-four by `5.58%` but worsens cloud-active all-eight by `0.33%`. RICO and
some stratocumulus variants remain counterexamples. Because `1e-4` was the
best member after inspecting a predeclared threshold sweep, this is a
candidate to validate, not a fair final generalization estimate.

Authoritative report:
[`output_pdf_bakeoff/crosscase/frozen_arm_profile_stats.md`](output_pdf_bakeoff/crosscase/frozen_arm_profile_stats.md).

### 3. Sparse independent raw SAM 3-D geometry

Twelve unfamiliar snapshots—three DYCOMS RF02, three RICO, and six LBA—were
scored with frozen ARM vectors. These data test actual joint density and tail
geometry, but are too sparse for an independent time-block fit.

Calibrated versus original prescribed-tail ADG3 improved both transport NRMSE
and density skill in all three regimes. Against the stronger baselines, however,
the result is mixed: it improves on ADG1 and Bergmann in DYCOMS (`0.4077`) but
loses there to Skw-vertical ADG3 (`0.3044`); it is best in RICO (`0.5249`);
and it loses to Bergmann in LBA (`0.3169` versus `0.2220`). Its density skill
is negative in every case (`-0.1718`, `-0.2632`, `-0.2807`), whereas Bergmann
is safer and positive in RICO (`+0.2625`). This confirms that transport fitting
and joint-shape fidelity are distinct objectives.

Authoritative artifact:
[`output_pdf_bakeoff/sparse_raw_frozen_crosscase.json`](output_pdf_bakeoff/sparse_raw_frozen_crosscase.json).

### 4. FIRE and the unresolved “FITT” name

No literal FITT case, dataset, paper, or history entry was found in the
repository, benchmark archive, `~/Code`, or `~/Downloads`. The most plausible
local interpretation is the GCSS **FIRE** stratocumulus case. Its COAMPS
profile file contains the complete operational moment inputs and the four
direct cloud-water targets. All 2,400 planes in the registered 61–120 minute,
up-to-1000 m comparison window were finite and moment-PSD. It does not contain
raw 3-D samples or conditional in-cloud fields, so it can support a direct-four
profile audit but not density or all-eight scoring. No completed FIRE closure
score is claimed here.

The alternative published “FiTT” fixed-tropopause-temperature dataset was
also inspected; its DAM profiles contain none of the moments or cloud-water
transport targets needed for this study. If “FITT” means another source, its
exact path or expanded name is still required.

### 5. Regime-group holdout status

[`experimental/pdf_bakeoff/utilities/sam_pdf_group_holdout.py`](experimental/pdf_bakeoff/utilities/sam_pdf_group_holdout.py)
implements the required next protocol: related RF02 variants form one outer
group; coefficient vectors come from a target-independent bounded design;
each outer fold selects coefficients and the discrete guard using only the
remaining groups; the held-out group is read once for reporting.

A very small protocol smoke completed, proving the loading, grouping,
candidate-bank, and reporting path, but its resolution is non-reportable. The
larger fast run evaluated its candidate bank on 1,760 planes, then stopped at
the first outer selection because every fitted-ADG1 candidate had at least one
audited-invalid prediction in the non-ARM training groups. The selector
correctly requires zero invalid states; the held-out ARM targets were therefore
never used and no result artifact was written. **No group-holdout performance
claim is made.** The immediate blocker is a predeclared, exact safe-fallback
and coverage policy, not permission to weaken the moment or PSD audits.

## What calibrated-skew-tail ADG3 still needs

The five-control calibrated family is a research hypothesis, not yet an
operational closure. Its third Gaussian represents a rare scalar-skew tail;
the remaining mean, full covariance, and `wp3` are closed by a residual ADG1
pair. The calibration exposes residual `gamma`/`beta`, tail capacity, the
covariance-regression displacement scale, and the normalized-`wp3` vertical
gain. Before promotion it still needs:

1. **A fair final comparison with centroid-tilt ADG2.** Both families must use
   the same frozen case set, regime-group outer folds, targets, fallback
   contract, and raw-density checks. The current evidence does not establish
   that the extra Gaussian is preferable to the strongest 2G diagnosis.
2. **Universal exact/PSD coverage.** Every eligible PSD moment state must
   return either a valid 3G mixture or an exact, PSD-safe ADG1 result. Fallback
   frequency and location must be reported separately from fit error; hard
   planes must not be silently removed.
3. **A genuinely held-out activation guard.** The current promising guard
   activates the tail only for unsaturated but tail-capable moment states, but
   its probability threshold was examined retrospectively. It must be chosen
   inside each training fold, remain inactive in Wangara, and avoid abrupt
   timestep-to-timestep switching, potentially through a smooth transition or
   carefully tested hysteresis.
4. **A physical tail-occupancy diagnosis.** Tail weight is presently
   `tail_capacity` times the largest weight that leaves a realizable residual
   ADG1 pair. That is primarily a mathematical capacity, not a direct plume
   occupancy estimate. Scalar skewness, Gaussian cloud probability, and the
   previous-below occurrence signal should be tested as target-blind occupancy
   predictors without using raw cloud truth at runtime.
5. **An explicit scalar-third-moment contract.** `rtp3` and `thlp3` diagnose
   the tail direction, but the completed mixture is not constrained to
   reproduce those two moments. The discrepancy must either be made small by
   construction or accepted and documented as part of the closure.
6. **Independent raw-3D and prognostic validation.** More complete SAM regimes
   must test joint-density geometry and conditional cloud statistics. Then
   unforced ARM, BOMEX, DYCOMS, RICO, Wangara, and FIRE SCM runs must test
   stability, budget closure, guard continuity, feedback on the prognosed
   moments, and downstream buoyancy and microphysics—not only offline LES
   reconstruction.
7. **Parameter-identifiability and reduction tests.** ARM ablations indicate
   that `reference_scale` supplies most of the calibrated improvement while
   `vertical_skw_gain` remains close to its prescribed value. Parameters that
   do not retain tight held-out boxes or demonstrable independent value should
   be fixed or removed.

Fortran implementation should be staged. A closure-only prototype is
moderate work: add an experimental `iiPDF` option, port the moment diagnosis,
reuse ADG1 for the residual pair, sum the existing analytic Gaussian cloud
integrals over three components, and retain an exact ADG1 fallback. A complete
first-class implementation is substantially larger because CLUBB's
`pdf_parameter` storage, accelerator data movement, restart serialization,
statistics, SILHS category sampling, and parts of microphysics are hard-coded
for two components. The displaced, correlated tail cannot generally be folded
into a single mixture-weight correction as readily as Bergmann's mean-centered
third component. Full three-component infrastructure should therefore follow,
not precede, a successful closure-only prognostic prototype.

### Coefficient-reduction check

A final projection ablation gives one clean simplification. Replacing the
fitted `vertical_skw_gain=2.5385` by **2.5**, while leaving the other four
controls unchanged, moves ARM validation NRMSE from 0.36439 to 0.36608 and
test NRMSE from 0.25539 to 0.25576. Test density skill is effectively
unchanged (-0.06457 to -0.06470). In a reduced 100-plane-per-case frozen
ten-case check with the structural gate but no probability floor, direct-four
NRMSE changes from 0.46146 to 0.46014 and all-eight NRMSE from 0.48941 to
0.48834. These are projection results rather than a new full fit, but they
agree with the independent tail analysis and justify fixing the gain at 2.5
in the next same-protocol fit.

The other controls are not yet safely removable. Projecting
`reference_scale` from 0.6124 to 0.87 raises ARM test NRMSE to 0.26049 and
worsens density skill to -0.07395; the independently fitted three-control
parent reaches 0.26346. Fixing residual ADG1 to `gamma=0.25`, `beta=1` can
slightly improve the eight diagnosed scalar/cloud terms, but its ARM density
skill degrades to about -0.22 to -0.26. That is evidence of a diagnostic-score
versus joint-geometry tradeoff, not evidence that gamma and beta are
irrelevant.

Therefore the recommended reduced continuous vector is
**`gamma`, `beta`, `tail_capacity`, and `reference_scale`**, with
`vertical_skw_gain=2.5`. Tail occurrence should use the defensible structural
condition already in the analyzer—activate only for `alpha_chi < 0`—and an
exact PSD-safe fallback. Do not promote an empirical hard Gaussian-cloud-
probability floor yet: `1e-4` mainly repairs dry-case averages and worsens the
cloud-active cohort, while the low-resolution group smoke selected either no
floor or `1e-8`. The floor remains an inner-fold validation question.

## Recommendation

### New causal prototype: upwind-state ADG2

The bake-off now includes **Upwind-state ADG2**, a five-coefficient causal
predictor-corrector experiment.  Its first three coefficients are the existing
centroid-tilt ADG2 controls.  `upwind_gain` interpolates the previous labeled
transport Gaussian toward the physically upwind neighboring level using
`clip(|w_transport| dt/dz, 0, 1)`.  `history_strength` blends that advected
prior with the current local centroid-tilt diagnosis.  Gaussian sufficient
statistics are blended, not the plotted ellipse parameters.

The corrector then uses the current supplied mean, complete covariance, and
`wp3` to rebuild both components on the exact two-Gaussian realizability
manifold.  If the predicted center lies outside that manifold, the historical
correction is reduced toward the guaranteed-realizable local parent; a final
fallback is exactly centroid-tilt ADG2.  The offline context contains only the
previous minute's `(w, rt, thl)` means, covariance, marginal third moments,
pressure, and linearized-saturation coefficients at the same/lower/upper
levels.  It contains no raw cloud water or fitted cloud target.

This is deliberately a first operationally plausible test, not yet a new
prognostic CLUBB closure.  The bake-off reconstructs each predecessor's
transport component from predecessor moments.  A true CLUBB implementation
should retain the labeled component state across timesteps, advect it, define
initialization/restart and component-label continuity, and test whether the
moment projection repeatedly erases the useful memory.  Held-out improvement
over centroid-tilt ADG2 is required before adding that state to Fortran.

### Upward moist-flux packet follow-up

**Flux-packet ADG2** tests a more physical version of that causal idea without
assuming that the lower-level transport population is already represented by
one labeled Gaussian.  For both predecessor Gaussians it analytically
moment-matches the distribution weighted by `(w-w_cutoff)+`.  A component
contributes only when that positive-velocity part carries a positive total-
water anomaly.  Consequently, a near-mean Gaussian with strong internal
`w-rt` tilt can initiate the packet near the surface, while a separated wet
Gaussian naturally dominates aloft.  The incoming lower-level packet is
combined with the retained previous same-level transport component, blended
with the current local centroid diagnosis, and projected back onto the exact
current mean, full covariance, `wp3`, and component-PSD manifold.  The two new
controls are packet strength and a normalized vertical-velocity cutoff;
packet strength zero is exactly Centroid-tilt ADG2.

On the unchanged 58-block ARM study, validation selected that exact parent:
`packet_strength=0` and `w_cutoff_sigma=0`.  The frozen test aggregate NRMSE
was `0.48434` with density skill `-0.49734`; validation objective was
`0.40624`.  Height-local diagnostic fits did select strong packet influence at
220--380 m and weak influence at a few higher levels, but one global causal
coefficient did not generalize.  This is useful negative evidence: allowing
both predecessor components to seed an analytically correct upward moist
packet fixes the low-level initiation objection, yet the exact current-moment
projection and a global strength still leave no held-out aggregate advantage.
The family remains in the dashboard for geometry and localized-event study,
but it is not a Fortran candidate on this result.

### Causal-packet ADG3 follow-up

**Causal-packet ADG3** separates the causal predictor from the current-level
skew constraint more cleanly.  Gaussian 1 is the retained same-level plus
incoming lower-level moist-flux packet diagnosed above.  Its center,
covariance, and orientation are held fixed.  Only its occupancy may be reduced
to keep the residual realizable; a fresh ADG1 pair then reproduces the
remaining current mean, full covariance, and `wp3` exactly.  This allows a
small positive moist updraft to coexist with negative aggregate `wp3`, which
an exact two-Gaussian corrector often accomplishes only by reversing the wet
component's vertical center.

The construction succeeds on the motivating 498-minute / 2140 m ARM plane:
despite negative grid-box `wp3`, it retains a packet near `w=+2.1 m/s` and
reduces its feasible occupancy to about 2%.  All supplied moments remain exact
and every component remains PSD.  This confirms that the third component fixes
the specific **constraint collision** we intended to isolate.

The unchanged fast 58-block ARM study nevertheless gives a mixed global
result.  Validation selects a nonzero packet (`packet_strength=0.199`,
`w_cutoff_sigma=1.209`), so the causal signal is not simply discarded.  Frozen
test aggregate NRMSE is `0.56678`, with `wprcp=0.47045`, `wp2rcp=0.40611`,
`rtprcp=0.68107`, and `thlprcp=0.71776`; density skill is `-0.12696`.  This is
better than ADG1/Bergmann on the same split and matches some difficult onset
planes well, but it trails cloud-anchored ADG3 (`0.4279`) and upwind-state ADG2
(`0.4552`).  At the rare high-cloud plane it preserves the correct sign of the
updraft but overshoots its transport magnitude.  Height-local fits favor
nonzero packet strength mainly from about 1180--2300 m, while several lower
levels return to the ADG1 parent.

The conclusion is narrower than “causality failed”: predecessor state can
identify a geometrically meaningful packet, and a third residual component
prevents current `wp3` from destroying it.  The remaining failure is packet
**amplitude/occupancy calibration**—one global strength and cutoff cannot
distinguish a useful arriving plume from a high-level packet that must be much
more strongly attenuated.  Keep this family as the clean causal-geometry
diagnostic in the dashboard, but do not promote it to Fortran on the present
held-out result.

A plane-by-plane held-out audit localizes that statement further.  The worst
errors are concentrated in the 496-minute upper cloud (1820--2300 m) and the
436-minute rising plume (1500--1980 m).  At 496 min / 1980 m, the reconstructed
"retained packet" has weight 0.098 but cloud probability only `0.00019`; its
center is `(w,rt)=(+0.41,0.00879)`, whereas the raw cloud-water-weighted SAM
tail is near `(+1.68,0.01479)`.  Exact residual closure then puts nearly all
cloud in a `w=-1.81 m/s` residual Gaussian and predicts negative `wprcp`.  The
cause is label loss: the offline predictor re-diagnoses the previous
centroid-tilt component and retains it unconditionally even when it is a large
dry population, while the genuinely moist incoming flux is much smaller.

The 436-minute failure is the opposite.  At 1820 m the packet is genuinely
cloudy, but the `1.21 sigma` velocity cutoff conditions on only its fastest
part: predicted `w=+3.92 m/s` versus the raw tail's `+1.84 m/s`.  Global packet
strength then leaves only 0.19% occupancy versus 2.34% raw cloud fraction.  It
therefore predicts a cloud that is too fast and too sparse.  Across the test
set, realizability is almost irrelevant (no fallback and only 1.6% occupancy
limiting); the nominal predictor weight reaches its 0.49 cap on about half the
planes.  The next causal prototype should retain the **actual labeled packet
state**, use saturation/coherent-transport information to reject a dry
retained component, and separate donor-occurrence gating from the packet mean
and covariance rather than using the same hard velocity truncation for both.

The positive counterexample is equally important.  At 418 min / 1300 m,
Causal-packet ADG3 cleanly separates a fully cloudy rising packet
(`w=+1.50 m/s`), an intermediate cloudy shoulder (`w=+0.75 m/s`), and the
weakly cloudy background (`w=-0.25 m/s`).  This is the first tested family to
recover that three-part raw-SAM geometry without forcing the tail/background
pair into nearly singular 2G ellipses.  Its amplitude is still imperfect—the
packet is too fast and too lightly weighted, so `wp2rcp` remains low—but the
example confirms that preserving a causal third component supplies useful
representational capacity rather than merely improving a scalar score.

### Coherent-packet ADG3 follow-up

**Coherent-packet ADG3** tests the two changes implied by that audit while
keeping Causal-packet ADG3 available as a separate comparison.  Its exact
neutral parent is Cloud-anchored ADG3.  A reconstructed predecessor component
may replace the local saturation tail only when it carries positive upward
total-water transport and has either analytic cloud support or positive
internal `w-rt` correlation above a fitted coherence threshold.  The
positive-velocity integral controls only how much probability crosses the
level boundary; the donor component's unconditioned center, covariance, and
tilt define packet geometry.  A residual ADG1 pair still closes the current
mean, full covariance, `wp3`, and component PSD exactly.

The unchanged 58-block ARM fit selected a nonzero packet strength of `0.677`
and coherence threshold `0.291`, along with `gamma=0.576`, `beta=0.734`, local
tail displacement `1.567`, local capacity `0.490`, `r_chi=0.383`, and centroid
skew gain `0.933`.  Frozen test aggregate NRMSE is `0.32054`, substantially
below Cloud-anchored ADG3 (`0.4279`), Upwind-state ADG2 (`0.4552`), and the
original Causal-packet ADG3 (`0.5668`).  Direct test NRMSEs are `0.351` for
`wprcp`, `0.368` for `wp2rcp`, `0.251` for `rtprcp`, and `0.252` for
`thlprcp`; density skill is `-0.0948`.  Validation independently improves to
`0.3393`, so this is not a test-only accident or a collapse to the parent.

The mechanism matters more than the headline score.  The coherent packet is
active on only 21 of 192 untouched test planes (10.9%); otherwise the method
returns its local parent.  It rejects all three reconstructed contributors at
the known 496 min / 1980 m dry-packet failure and therefore avoids the causal
model's negative-`wprcp` catastrophe.  At 436 min / 1820 m it keeps the cloudy
packet but lowers the predecessor-derived vertical center from the original
causal diagnosis near `3.92` to `1.26 m/s` under fitted controls, closer to the
raw plume.  At 418 min / 1300 m it retains the useful three-population
geometry and raises predicted `wp2rcp` from about `3.17e-5` to `4.96e-5`
against SAM's `8.00e-5`.

This version is promising but not finished.  Its principal remaining failure
is the opposite edge of the occurrence/geometry separation: a predecessor
component can pass the coherence gate because its *positive-velocity part*
carries moist transport while its unconditioned Gaussian center is still at
negative `w`.  Preserving that center literally produces a downward packet at
some upper-cloud planes (for example 498 min / 2140 m and 631 min / 2300 m).
The next controlled extension should leave scalar center, covariance, and
tilt unconditioned, but apply a one-sided vertical-center rescue only when a
qualified donor center lies below the crossing velocity.  That would
interpolate between the two established limits: full truncation is too fast
and sparse, while fully unconditioned geometry can point a real crossing
packet downward.  It should remain a separate family until held-out blocks
show that the extra rule improves active-event skill without weakening the
new gate.

#### Dense time-height packet audit

A dedicated Notes laboratory now evaluates the frozen global vector at all
870 saved ARM minutes and all 64 native heights from 60--2580 m, rather than
only the 928 planes in the fitting protocol.  Its linked altitude profiles,
time-height maps, raw-SAM contours, and donor ledger separate occurrence,
placement, geometry, occupancy, and residual closure.  On the dense subset
inside untouched test-time blocks, 1195 of 1478 active packets improve the
direct score.  The gate is therefore useful but less uniformly precise than
the sparse 19-of-21 audit suggested.

The new instrumentation identifies three distinct remaining failures:

1. **Donor-mass dilution is the clearest correctable geometry error.** Of 45
   held-out active planes with all three reconstructed contributors, 39 are
   harmful and their median packet `rt`-center error is 5.55 g/kg.  The packet
   is currently assembled with crossing probability as its geometry weight.
   A broad, nearly dry lower background can therefore overwhelm a rare wet
   donor even when the wet donor carries more positive moisture flux.  At
   497 min / 2140 m this produces `rt=7.58 g/kg` instead of the raw tail's
   `14.40 g/kg`; weighting the same accepted donor geometry by positive moist
   flux moves the diagnostic center more than 5 g/kg toward the wet tail.
   The next target-blind variant should use moist flux for geometry and retain
   crossing mass only for occupancy.
2. **Moment re-diagnosis loses packet identity in coherent cloud columns.** At
   496 minutes the raw cloudy-tail center remains near `w=1.8 m/s,
   rt=15.1 g/kg` at 1820 m and decays smoothly aloft, but reconstructed
   predecessor components lose both cloud support and internal `w-rt` tilt.
   The gate falls back continuously from roughly 1740--2380 m, then abruptly
   reactivates above it.  Similar gaps occur around 441--449 minutes.  Lowering
   the global gate would also admit dry backgrounds; this evidence instead
   favors retaining an explicitly labeled packet state across calls.
3. **Positive crossing does not guarantee a positive Gaussian center.** Among
   170 downward-centered held-out activations, 90 are harmful.  This remains a
   separate one-sided vertical-center problem.  Only 169 of 1478 activations
   use less than 95% of requested occupancy, so realizability clipping is
   secondary to donor selection and geometry.

This changes the next experiment order.  First test flux-weighted geometry on
the already accepted donor set; it adds no new source variable and directly
targets the dominant three-contributor failure.  Next preserve labeled packet
state only to bridge diagnosed continuity gaps.  Test a one-sided `w`-center
rescue separately.  Do not weaken the occurrence gate or add another global
capacity coefficient until those three mechanisms have been isolated.

### Focused centroid-conditioning and occurrence-prior tests

Two controlled follow-ups separate local two-Gaussian geometry from causal
tail occurrence.

First, **Conditioned-centroid ADG2** uses exactly the same three fitted controls
and analytic centroid target as Centroid-tilt ADG2, but selects the exact
feasible center/tilt allocation with the largest normalized component Schur
log product.  It removes the most obvious nearly singular ellipses: on the ARM
reference plane the weakest normalized Schur margin rises from numerical zero
to `0.756`.  On the unchanged held-out split, however, aggregate NRMSE worsens
from `0.46325` to `0.48908`.  Density skill improves from `-0.5191` to
`-0.4238`, `wprcp` improves from `0.453` to `0.444`, and in-cloud-w NRMSE from
`0.914` to `0.756`, but `wp2rcp` worsens from `0.354` to `0.525`.  The skinny
shape was therefore partly a coefficient-free selection artifact, but maximum
conditioning alone is not the transport closure.

Second, an isolated **upwind occurrence prior** freezes the already fitted
Centroid-tilt coefficients for the current and predecessor diagnoses and
advects only the labeled transport-component weight.  Validation selects a
moment-only positive transport gate, full upwind interpolation, and history
strength `0.5`.  It changes only 14 of 192 untouched test planes and improves
aggregate NRMSE from `0.46325` to `0.45564` (`1.64%`).  On nine arriving-plume
planes it wins seven and loses one, improving all four direct transport
diagnostics.  Suppressing seven departing tails improves cloud fraction but
worsens `wprcp` and `wp2rcp`, so overall `wprcp` NRMSE rises slightly.  This is
mechanism evidence for a one-sided arrival prior, not permission to tune that
new branch on the already inspected split; an arrival-only rule needs a fresh
block assignment.

A subsequent **saturation-front phase experiment** tested a stricter version
of that idea after diagnosing the difficult 294-minute ARM onset layer.  It
used the independently derived causal probability
`logit(P)=4.66+0.90 alpha_current+2.62(alpha_previous_below-alpha_current)` to
reduce tail capacity, vertical displacement, or both in Centroid-tilt ADG2 and
Calibrated-skew-tail ADG3.  Two fresh complete-block assignments refit both
parents before validation selected among 12 predeclared responses.  In both
assignments, **the literal parent won validation for both families**.  The
modulation also lost within the descriptive moving-front cohort itself.  A
small vertical-displacement benefit appeared only in already mature
`P>=0.95` states, which is not enough to justify a new gate.  This rules out
using the present saturation-front probability as a broad multiplier on tail
weight or center displacement; it does not rule out transporting a persistent
labeled component state.

That remaining possibility was then tested directly in a standalone
**persistent-component experiment**.  At three-minute cadence, the accepted
tail component's weight, mean anomaly, and covariance were retained through
each 15-minute block, optionally blended with the vertically upwind prior, and
projected back onto the current exact mean/covariance/`wp3` manifold.  The
candidate bank included full-state and center-only memory down to strength
`0.02`; state was reset at every split boundary.  On two independent block
assignments, Centroid-tilt ADG2 selected either the parent or a barely better
validation candidate that worsened test NRMSE by `0.98%`.  Calibrated-tail
ADG3 selected weak temporal-only memory, but test NRMSE changed by `+0.14%`
and `+2.49%`.  Stronger memory and vertical upwind blending were consistently
worse.  Persistence often reduced center jitter or improved density, but it
degraded the diagnosed transport terms.  Therefore this direct state-retention
form is not promoted either: smooth component identity is not the same as
predictive conditional plume information.

Together these tests sharpen the modeling division: current moments should
diagnose local placement, while previous-below state may diagnose the arrival
or persistence of a rare transport population.  That history should not be
collapsed into a generic scalar multiplier on the current diagnosis.  None of
the 2G follow-ups closes the large ARM diagnostic gap to calibrated-skew-tail
ADG3 (`0.25539`).

Do not replace CLUBB's operational PDF yet. Preserve three finalists:

1. **Safe ADG1** as the operational and dry-regime baseline.
2. **Strict Bergmann 3G** as the conservative implementable shape baseline.
3. **Guarded calibrated prescribed-skew-tail ADG3** as the leading transport
   hypothesis.

For the next implementation-oriented bake-off, use the reduced four-control
ADG3 vector (`gamma`, `beta`, `tail_capacity`, `reference_scale`) with
`vertical_skw_gain=2.5` as the primary candidate.  Retain conditioned-centroid
ADG2 as a geometry/conditioning control and test a predeclared one-sided
arrival prior on a fresh split.  Do not combine the new prior and ADG3 until
each has independently survived that test.

The immediate next step is to define and test a target-blind coverage contract:
every coefficient family must return an exact, PSD safe fallback on every
eligible PSD input, and the report must expose fallback rate separately from
fit quality. Do not discard difficult planes or loosen the exactness audit.
Then rerun the leakage-free regime-group holdout over the profile cases, with
FIRE included as direct-four-only evidence and RF02 variants kept in one group.
The guard threshold and all continuous coefficients must be selected inside
each outer fold. Acceptance should require:

- lower held-out direct-four error in most physical regime groups, not only a
  lower equal-case average;
- no false tail in Wangara;
- no material degradation of all-eight cloud diagnostics where available;
- frozen-vector raw-density/tail checks on DYCOMS, RICO, and LBA;
- exact supplied moments and PSD components on every accepted state; and
- reasonably tight outer-fold coefficients, especially tail capacity,
  reference scale, and vertical-skew gain.

If that test fails, the defensible conclusion is that the currently prognosed
local moments cannot robustly identify rare-tail occupancy and displacement.
The next model change should then be an additional prognostic/causal plume
state—not another globally tuned Gaussian-shape coefficient. The previous-
below saturation signal is the most promising occurrence predictor, but the
evidence says it should gate a tail rather than determine the tail's vertical
center.

## Reproducibility and implementation map

- Family definitions, exact-moment constructors, PSD guards, and raw tail
  diagnostics: [`experimental/pdf_bakeoff/utilities/sam_pdf_bakeoff.py`](experimental/pdf_bakeoff/utilities/sam_pdf_bakeoff.py)
- ARM complete-block fitting and height-local parameter boxes:
  [`experimental/pdf_bakeoff/utilities/sam_pdf_global_bakeoff.py`](experimental/pdf_bakeoff/utilities/sam_pdf_global_bakeoff.py)
- Frozen profile cross-case evaluator:
  [`experimental/pdf_bakeoff/utilities/sam_pdf_crosscase.py`](experimental/pdf_bakeoff/utilities/sam_pdf_crosscase.py)
- Sparse raw-3D evaluator:
  [`experimental/pdf_bakeoff/utilities/sam_pdf_sparse_raw_evaluator.py`](experimental/pdf_bakeoff/utilities/sam_pdf_sparse_raw_evaluator.py)
- Group-holdout protocol:
  [`experimental/pdf_bakeoff/utilities/sam_pdf_group_holdout.py`](experimental/pdf_bakeoff/utilities/sam_pdf_group_holdout.py)
- Standalone coefficient-reduction analyzer:
  [`experimental/pdf_bakeoff/utilities/analyze_calibrated_tail_reduction.py`](experimental/pdf_bakeoff/utilities/analyze_calibrated_tail_reduction.py)
- Standalone matched-parent upwind-occurrence experiment:
  [`experimental/pdf_bakeoff/utilities/sam_pdf_upwind_occurrence_experiment.py`](experimental/pdf_bakeoff/utilities/sam_pdf_upwind_occurrence_experiment.py)
- Standalone saturation-front phase experiment:
  [`experimental/pdf_bakeoff/utilities/sam_pdf_phase_front_experiment.py`](experimental/pdf_bakeoff/utilities/sam_pdf_phase_front_experiment.py)
- Standalone persistent-component experiment:
  [`experimental/pdf_bakeoff/utilities/sam_pdf_persistent_component_experiment.py`](experimental/pdf_bakeoff/utilities/sam_pdf_persistent_component_experiment.py)
- Full chronological evidence and verification log:
  [`worklog.md`](worklog.md)
