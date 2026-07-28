# Design ideas: vertically coherent PDF diagnosis

## Purpose

This note sketches small, robust ways to make a diagnostic PDF geometry more
continuous in height.  The target is a sequence of local conditional PDFs,

$$
p(w,r_t,\theta_l \mid z_k),
$$

not a literal four-variate PDF with height as a fluctuating random variable.

The motivation is physical.  A transport population can rise, entrain, cool,
lose moisture, detrain, and sometimes reverse vertical velocity over several
model levels.  Diagnosing each level independently can make a component jump,
rotate abruptly, or exchange identity with the other component even when the
underlying transport structure is vertically coherent.

The constraint is equally important: CLUBB's predicted moments are local in
height.  Any vertical coupling must remain a weak diagnostic preference.  It
must not overwrite the local grid mean, covariance, or a genuine cloud-base,
inversion, detrainment, or regime transition.

**Relevant source:**
[`src/CLUBB_core/mixing_length.F90`](../src/CLUBB_core/mixing_length.F90)

---

## 1. What should be coupled

Do not directly smooth component means or component covariance entries after
the PDF has been diagnosed.  That would generally break the exact mixture
mean/covariance identities.

Instead, couple a small set of **pre-reconstruction geometry controls**.  For
a two-Gaussian PDF, a useful local state is

$$
\boldsymbol{q}_k =
\begin{bmatrix}
d_{w,k} \\
d_{r_t,k} \\
d_{\theta_l,k} \\
\operatorname{logit}(a_k)
\end{bmatrix},
$$

where:

- $(d_w,d_{r_t},d_{\theta_l})$ is the standardized separation direction and
  magnitude;
- $a$ is the rare-component weight.

The locally diagnosed PDF supplies
$\boldsymbol{q}_{k,\mathrm{local}}$.  A vertically coherent procedure
returns a lightly adjusted $\boldsymbol{q}_k$.  The ordinary local PDF
reconstruction then uses that adjusted geometry and re-enforces the supplied
mean, full covariance, component-width constraints, and PSD caps at each
level.

This makes the coupling a hypothesis about *how to allocate* an already known
local moment state, rather than a modification of the moment state itself.

For trivariate transport PDF, a prudent first version should couple only the three separation
coordinates.  Keeping the component weight local avoids spreading a weak plume
population into a layer where the local skewness does not support one.

---

## 2. What `Lscale` can and cannot tell us

### 2.1 Current default meaning

The usual CLUBB setting has `l_diag_Lscale_from_tau = .true.`.  In that path,
the code diagnoses a local dissipation time scale from background, surface,
shear, and static-stability contributions.  It then defines

$$
L_k = \tau_k \sqrt{e_k},
$$

where $e_k$ is TKE.  The time scale is capped by the allowed mixing length;
in a host model that upper limit is one quarter of the smaller horizontal grid
spacing.  In standalone CLUBB the cap is effectively very large.

Thus, the default `Lscale` is a physically useful estimate of the size of
energy-containing turbulent motion.  It responds to mixing, shear, stability,
surface forcing, and TKE.  It is a sensible input for the *amount* of vertical
coherence to prefer.

It is not, however:

- a prognosed vertical covariance of PDF-component properties;
- a guarantee that the same plume persists unchanged over $L_k$;
- a reason to couple across a sharp inversion or a layer with no local
  transport evidence.

The direct, parcel-based mixing-length option can additionally diagnose
`Lscale_up` and `Lscale_down`.  The default tau-based path deliberately sets
those directional outputs to zero.  Therefore an initial implementation
should use the ordinary scalar `Lscale`; it should not require directional
lengths that are unavailable in the usual configuration.

### 2.2 Use a bounded coherence length

Define a conservative coherence length at an interface between levels,

$$
\ell_{k+1/2}
= \min\left[
c_L\,\frac{L_k+L_{k+1}}{2},
\ell_{\max}
\right].
$$

Here $c_L$ is a small tunable fraction and $\ell_{\max}$ is an explicit
safety cap.  The cap matters: a large turbulent length should make nearby
levels more mutually informative, not create a rigid vertically uniform PDF
through an entire boundary layer.

On nonuniform grids, all couplings should use the physical interface distance
$\Delta z_{k+1/2}$ rather than a number of model levels.

---

## 3. Design A: local three-level consistency correction

### Idea

This is the simplest and safest experiment.  At level $k$, form a prior from
only the immediate neighbors:

$$
\boldsymbol{q}_{k,\mathrm{prior}}
=
\frac{
w^-_k\boldsymbol{q}_{k-1,\mathrm{local}}
+ w^+_k\boldsymbol{q}_{k+1,\mathrm{local}}
}{w^-_k+w^+_k},
$$

and blend it with the local diagnosis,

$$
\boldsymbol{q}_k
= (1-\alpha_k)\boldsymbol{q}_{k,\mathrm{local}}
+ \alpha_k\boldsymbol{q}_{k,\mathrm{prior}}.
$$

The weights decay with distance relative to `Lscale`, for example,

$$
w^+_k = g_{k+1/2}
\exp\left[-\left(\frac{\Delta z_{k+1/2}}{\ell_{k+1/2}}\right)^2\right].
$$

The bounded blend factor $0\le\alpha_k\le\alpha_{\max}$ should remain
modest.  A first test might use only $\alpha_{\max}\le 0.2$.

### Gate and barriers

The interface gate $g_{k+1/2}$ should be zero unless both levels support a
meaningful transport geometry.  It can combine:

- a minimum TKE or $w$-variance condition;
- a smooth magnitude threshold for the local separation/skewness signal;
- a covariance-coherence condition, such as consistent positive
  $w'r_t'$ for a rising-moist plume branch;
- a transition barrier based on a strong static-stability jump, a cloud-top
  transition, or a large local change in the diagnosed rare-component weight.

The exact gate should be intentionally conservative.  It is better for the
first experiment to leave many levels unchanged than to impose a false plume
continuation.

### Computational character

This method has a strict one-level half-bandwidth.  It needs no column solve,
no new prognostic variable, and no storage beyond neighboring local diagnoses.
It is the best first diagnostic test because it isolates whether even a weak
continuity preference improves component placement.

### Limitations

It is symmetric and therefore does not distinguish upward from downward
propagation.  It also cannot resolve a component-label swap by itself if both
neighbors are ambiguous.  Those limitations motivate Design B.

---

## 4. Design B: gated block-tridiagonal MAP diagnosis

### Idea

Treat the local PDF diagnosis as an observation and seek the vertically most
coherent geometry that remains close to it.  For one column, minimize

$$
J =
\frac{1}{2}\sum_k
(\boldsymbol{q}_k-\boldsymbol{q}_{k,\mathrm{local}})^\mathsf{T}
\mathbf{W}_k
(\boldsymbol{q}_k-\boldsymbol{q}_{k,\mathrm{local}})

 + \frac{1}{2}\sum_k \gamma_{k+1/2}
 \left\|
 \mathbf{S}_{k+1/2}
 (\boldsymbol{q}_{k+1}-\boldsymbol{q}_k)
 \right\|^2.
$$

$\mathbf{W}_k$ is a local-confidence matrix.  It should be large when local
moments strongly constrain the geometry and small only where the local problem
is genuinely ambiguous.  $\mathbf{S}$ nondimensionalizes the coupled
coordinates so moisture units cannot dominate the solve.

The interface coefficient can be

$$
\gamma_{k+1/2}
= \gamma_0\,g_{k+1/2}
\exp\left[-\left(\frac{\Delta z_{k+1/2}}
{\ell_{k+1/2}}\right)^2\right].
$$

The normal equations are block tridiagonal.  For the three-coordinate
separation vector, every block is only 3 by 3; including the component weight
makes it 4 by 4.  The solve is $O(n_z)$ per column and uses only neighboring
off-diagonal blocks.

### Why a tridiagonal solve can still be appropriate

The matrix has only one block off-diagonal, even though the solution has a
decaying influence over a connected segment.  This is desirable: it produces
smooth, physically scaled coherence without constructing a wide dense matrix
or introducing an arbitrary many-level stencil.

To keep the *effective* influence local, do all three of the following:

1. cap $\ell$ as described above;
2. set $g_{k+1/2}=0$ at physical barriers, splitting the column into
   independent short segments;
3. retain a substantial local-confidence floor in $\mathbf{W}_k$.

With these choices, the procedure is a weak regularizer, not a vertical
averaging operation.

### Component identity

Before forming each interface penalty, compare the two possible component
pairings across the interface.  Choose the pairing with the smaller
standardized mean/covariance distance, but only if both pairings are locally
realizable.  This lets a physical moist transport component retain an identity
when the numerical labels G1/G2 would otherwise swap.

For trivariate transport PDF, label continuity should be a *diagnostic association* only.  The
final local component weights and covariances must still be reconstructed from
the adjusted local geometry and checked for PSD.

### Limitation

This is a genuine cross-level diagnostic solve.  It is fast, but it adds a new
vertical dependency and should not be the first production trivariate transport PDF option.
It should first be evaluated offline against LES component-placement metrics
and cloud/transport diagnostics.

---

## 5. Optional later extension: directed one-step transport prior

If trustworthy directional reach information becomes available, use it only as
an additional prior.  For a rising transport component, a level can receive a
one-sided prediction from below; for a descending component, from above.  The
prior should decay over a direction-specific reach length and be gated by the
component's diagnosed vertical velocity.

This has an appealing plume interpretation, but it should be deferred:

- default `Lscale_up` and `Lscale_down` are not available in the standard
  tau-based mixing-length path;
- a component can slow, detrain, or reverse, so vertical velocity does not
  identify an unbroken trajectory;
- it introduces more assumptions than Designs A or B.

The parcel trajectory ideas already present in the PDF-9 and direct
mixing-length infrastructure are useful inspiration, but should not silently
become a trivariate transport PDF dependence.

---

## 6. Recommended staged experiment

1. **Offline diagnosis only.**  Save the local trivariate transport PDF geometry over height,
   apply Design A after the fact, reconstruct each level, and compare with
   matched SAM planes.  Do not affect CLUBB tendencies or the running model.
2. **Measure the right quantities.**  Score component-center trajectories,
   component-label swaps, covariance orientation, cloud fraction, and
   $w'r_c'$ alongside ordinary moment conservation.  A smoother picture alone
   is not evidence of improvement.
3. **Use `Lscale` only as a bounded scale.**  Sweep $c_L$ and
   $\ell_{\max}$ across ARM, BOMEX, and a stratocumulus case.  Verify that
   barriers prevent cross-inversion and cross-cloud-layer contamination.
4. **Only then test Design B.**  Keep the solve block tridiagonal, start by
   coupling the three separation coordinates only, and require a local
   confidence floor plus PSD fallback at every level.
5. **Decide whether a prognostic extension is warranted.**  If the gains come
   only from smoothing ambiguous local geometry, a diagnostic solve is likely
   enough.  If useful coherence has persistent time memory, then consider a
   separate prognostic component-geometry state; that is a much larger model
   change and should not be inferred merely from this experiment.

---

## 7. Bottom-line recommendation

Start with Design A: a gated, one-neighbor, `Lscale`-bounded correction to the
trivariate transport PDF separation geometry, followed by the existing local moment/PSD
reconstruction.  It has the smallest stencil, no matrix solve, and gives a
clean test of whether vertical coherence contains useful information.

If that helps without smearing physical transitions, Design B is the natural
next step.  Its block-tridiagonal system is still small and fast, while being
substantially more principled than applying repeated local smoothing.
