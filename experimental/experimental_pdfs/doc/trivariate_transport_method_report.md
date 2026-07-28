# trivariate transport PDF: A trivariate transport-oriented two-Gaussian closure for CLUBB

## Technical review report

**Status:** implemented research closure (`iiPDF_type=9`), not yet a
validated replacement for ADG1.  This document describes the algorithm that
is currently in the code, its intended atmospheric interpretation, and the
scientific questions that should be resolved before a paper or operational
adoption.

**Primary implementation:**
[`src/CLUBB_core/trivariate_transport_pdf.F90`](../src/CLUBB_core/trivariate_transport_pdf.F90)

**Scope:** the thermodynamic portion of CLUBB's two-Gaussian PDF for

$$
\boldsymbol{x} = (w, r_t, \theta_l)^\mathsf{T},
$$

where $w$ is vertical velocity, $r_t$ total water mixing ratio, and
$\theta_l$ liquid-water potential temperature.  trivariate transport PDF diagnoses two
trivariate Gaussian components from the predicted mean, covariance, and three
marginal skewnesses.  It then uses CLUBB's existing analytic Gaussian cloud
integrals.  It does **not** introduce numerical quadrature or a nonlinear
optimization.  trivariate transport PDF now also supplies an analytic scalar-third-moment flux
to a standard one-dimensional prognostic advance; it does not require a new
global or nonlinear solver.

---

## 1. Executive summary

### 1.1 Scientific problem

Many cloudy boundary-layer and lower-free-tropospheric LES samples contain a
transport population that is distinctive primarily in moisture rather than in
its mean vertical velocity.  Such a population can be far from the environment
in $r_t$, broad in $w$, and strongly tilted in the $w$-$r_t$ plane.
It may therefore contribute substantially to moisture and cloud transport even
when its mean $w$ is not extremely separated from the grid mean.  The
desired geometry is a slanted, moisture-displaced transport lobe adjacent to a
broad environmental population.

ADG1 is intentionally simple and robust, but its construction is effectively
anchored in the vertical-velocity distribution.  Its component
$w$-scalar covariances are zero; grid $w$-scalar covariance is represented
principally through a separation of component means.  This makes it natural to
produce vertically stacked components and makes the target moisture-oriented,
internally tilted lobe difficult to represent.

### 1.2 Central idea

trivariate transport PDF retains two Gaussians, but diagnoses their centers jointly as a single
trivariate separation vector.  It then assigns residual covariance and
width contrast to the components while enforcing positive definiteness of both
component covariance matrices.  The rare component, called G1, is diagnosed
from the magnitude of moisture skewness and can carry a positive internal
$w$-$r_t$ covariance.  Thus, a strongly moisture-separated G1 with modest
mean-$w$ displacement is structurally allowed.

### 1.3 Main properties of the implemented method

| Property | trivariate transport PDF treatment |
| --- | --- |
| Mixture family | Two trivariate normal components in $(w,r_t,\theta_l)$. |
| Inputs | Means, full 3-by-3 covariance, and $w$, $r_t$, and $\theta_l$ skewnesses already available to CLUBB. |
| Weight | A bounded, moisture-skewness-based rare-component weight. |
| Center geometry | One covariance-metric-normalized trivariate direction, using all three skewnesses. |
| Coherent-plume structure | A gated, signed transport-plume blend for fresh warm/moist and mature moist/cool plumes, with a well-mixed background preference. |
| $w$ separation | Explicitly attenuated by a parameter, so scalar separation need not imply large mean-$w$ separation. |
| G1 internal tilt | Positive residual $w$-$r_t$ covariance can be allocated to G1. |
| Moment conservation | The grid mean and full 3-by-3 covariance are reconstructed to roundoff for every non-fallback state. |
| Third moments | Each marginal third moment is requested through an analytic width contrast, but can be softened by realizability limits. |
| Realizability | Fixed-size Cholesky checks and closed-form 3-by-3 eigenvalue caps; no iteration. |
| Failure behavior | Deterministic collapsed two-copy Gaussian fallback. |

The closure is therefore best viewed as a **constrained, moment-informed
geometric hypothesis**.  It is not an exact inversion of all available third
moments, nor does it claim that diagonal third moments uniquely identify a
physical plume population.

---

## 2. Atmospheric motivation and design requirements

### 2.1 Target transport geometry

The motivating LES geometry has four features:

1. A relatively low-probability transport population is displaced strongly in
   total water.
2. Its mean $w$ can remain close to the grid mean relative to its moisture
   displacement.
3. It has large internal $w$ variance, consistent with active rising
   parcels that slow, detrain, and can return downward.
4. Within that population, moist anomalies and upward motion are positively
   correlated.  The relevant structure is therefore a tilted ellipse, not
   merely two circular or axis-aligned blobs.

This distinction matters for cloud and moisture transport.  A closure can
reconstruct a grid $\overline{w'r_t'}$ either by separating component centers
or by retaining covariance within a component.  These choices are not
equivalent after nonlinear saturation and cloud-water calculations.  In
particular, a component's covariance with the saturation variable affects its
contribution to $\overline{w'r_c'}$, not merely the location of its center.

### 2.2 Why a two-Gaussian closure is retained

The aim is not to reproduce every cloud population in LES.  CLUBB needs a
pointwise closure that is fast, robust over many regimes, and at least as safe
to use as its standard PDF machinery.  A two-component Gaussian mixture
offers several practical advantages:

- it keeps analytic cloud integrals already used by CLUBB;
- all required operations are local in height and column;
- mean and covariance constraints have simple exact identities;
- realizability can be checked with small dense matrices;
- a conservative limiting state is obvious: two identical copies of one
  Gaussian.

The price is non-uniqueness.  A two-Gaussian distribution cannot independently
fit an arbitrary set of third and higher moments.  trivariate transport PDF makes that
non-uniqueness explicit by specifying how the limited degrees of freedom are
allocated.

### 2.3 Design requirements

The implemented design was selected to meet the following requirements.

| Requirement | Consequence in trivariate transport PDF |
| --- | --- |
| Moisture-oriented transport geometry | Moisture skewness controls the rare-component mass; the center direction gives $r_t$ full standing. |
| Modest mean-$w$ displacement remains possible | A bounded `w_direction_scale` attenuates the $w$ part of the direction. |
| Broad and tilted transport component | Width contrasts and positive residual $w$-$r_t$ covariance can be placed in G1. |
| No nonlinear solve | Algebraic center, width, and covariance formulas plus direct PSD caps. |
| Exact low-order second-moment consistency | Means and the entire trivariate covariance are constructed to reconstruct the supplied state. |
| Robust invalid-state behavior | Cholesky checks and a local collapsed fallback. |
| Existing CLUBB infrastructure | Reuse `pdf_closure_module` and its analytic Gaussian cloud integrations. |

---

## 3. Notation and available moments

All derivations below are pointwise in grid column and height.  Define the
grid mean and centered state as

$$
\overline{\boldsymbol{x}}
= (\overline{w},\overline{r_t},\overline{\theta_l})^\mathsf{T},
\qquad
\boldsymbol{x}' = \boldsymbol{x}-\overline{\boldsymbol{x}}.
$$

The supplied grid covariance is

$$
\mathbf{C} = \overline{\boldsymbol{x}'\boldsymbol{x}'^\mathsf{T}}
= \begin{bmatrix}
\overline{w'^2} & \overline{w'r_t'} & \overline{w'\theta_l'}\\
\overline{w'r_t'} & \overline{r_t'^2} & \overline{r_t'\theta_l'}\\
\overline{w'\theta_l'} & \overline{r_t'\theta_l'} & \overline{\theta_l'^2}
\end{bmatrix}.
$$

trivariate transport PDF also receives the three existing marginal skewnesses

$$
S_j = \frac{\overline{x_j'^3}}{\sigma_j^3},
\qquad
\sigma_j = \sqrt{C_{jj}},
\qquad j\in\{w,r_t,\theta_l\}.
$$

The method intentionally does **not** require mixed third moments such as
$\overline{w'^2r_t'}$ or $\overline{w'r_t'^2}$.  That is a practical
choice, not a statement that those moments are unimportant.  They would be
especially informative for diagnosing whether a moist population has enhanced
$w$ variance and internal tilt; see Section 13.

### 3.1 Standardized coordinates

Let

$$
\mathbf{D}=\mathrm{diag}(\sigma_w,\sigma_{r_t},\sigma_{\theta_l}),
\qquad
\mathbf{R}=\mathbf{D}^{-1}\mathbf{C}\mathbf{D}^{-1}.
$$

The algorithm works in the standardized correlation coordinates
$\mathbf{R}$.  Physical-unit component covariances are restored at the end
with $\mathbf{D}$.  This prevents a direction in the mixed units
$(\mathrm{m\ s^{-1}},\mathrm{kg\ kg^{-1}},\mathrm{K})$ from depending on
arbitrary unit scales.

The Fortran implementation symmetrizes $\mathbf{R}$, sets its diagonal to
one, and requires a positive-definite Cholesky factorization.  It does not
silently repair a non-positive-definite grid covariance.  Such a state takes
the fallback path described in Section 9.

---

## 4. The two-Gaussian family and exact second-moment identities

Let $a$ be the G1 mixture weight and $b=1-a$ the G2 weight.  The
component means use one standardized separation vector $\boldsymbol{d}$:

$$
\boldsymbol{\mu}_1 = \overline{\boldsymbol{x}} + b\mathbf{D}\boldsymbol{d},
\qquad
\boldsymbol{\mu}_2 = \overline{\boldsymbol{x}} - a\mathbf{D}\boldsymbol{d}.
$$

Thus, the difference of the component means is

$$
\boldsymbol{\mu}_1-\boldsymbol{\mu}_2 = \mathbf{D}\boldsymbol{d}.
$$

This construction exactly conserves the grid mean:

$$
a\boldsymbol{\mu}_1+b\boldsymbol{\mu}_2
=\overline{\boldsymbol{x}}.
$$

Let $\boldsymbol{\Sigma}_1$ and $\boldsymbol{\Sigma}_2$ denote the
standardized component covariance matrices.  The exact standardized mixture
covariance identity is

$$
\mathbf{R}
=a\boldsymbol{\Sigma}_1+b\boldsymbol{\Sigma}_2
+ab\boldsymbol{d}\boldsymbol{d}^\mathsf{T}.
\qquad (1)
$$

The final term is the rank-one between-component covariance.  It is the
mechanism that makes a single trivariate displacement vector scientifically
important: center separation affects all six covariance entries coherently.

Define the residual covariance after the centers are placed,

$$
\mathbf{R}_{\mathrm{res}}
=\mathbf{R}-ab\boldsymbol{d}\boldsymbol{d}^\mathsf{T}.
\qquad (2)
$$

Any two component covariances satisfying

$$
a\boldsymbol{\Sigma}_1+b\boldsymbol{\Sigma}_2
=\mathbf{R}_{\mathrm{res}}
\qquad (3)
$$

reconstruct the supplied covariance exactly.  trivariate transport PDF parameterizes a
mean-preserving contrast about this residual, which makes both conservation
and realizability transparent.

---

## 5. Diagnosing the transport weight and center separation

### 5.1 Bounded skewness preprocessing

To avoid allowing rare extreme samples to control the closure, the raw
skewnesses are first capped and softened:

$$
S_j^{\mathrm{cap}}=\min\bigl(12,\max(-12,S_j)\bigr),
\qquad
\widetilde{S}_j=\tanh\left(\frac{S_j^{\mathrm{cap}}}{1.25}\right).
\qquad (4)
$$

The cap prevents a pathological skewness from producing arbitrarily large
width requests, while the hyperbolic tangent makes the response bounded and
continuous.  These are deliberately pragmatic regularizations; they should
be tuned and tested rather than interpreted as physical constants.

### 5.2 Moisture-primary rare-component weight

G1 is assigned a moisture-tail weight from the **magnitude** of total-water
skewness:

$$
q_a=\tanh\left(
g_m\frac{|S_{r_t}^{\mathrm{cap}}|}{1.25}
\right),
\qquad
a=\mathrm{clip}\left[\frac{1-q_a}{2},\ 0.035,\ 0.5\right],
\qquad (5)
$$

where $g_m$ is `trivar_moisture_tail_gain`.  Stronger moisture skewness
therefore produces a rarer G1, but the weight cannot become smaller than
0.035 or exceed one half.

The sign of moisture skewness does not change which component is called G1;
it changes the direction of its moisture displacement.  A positive
$S_{r_t}$ makes G1 the moist tail in the usual convention.  A negative
$S_{r_t}$ makes the same rare G1 a dry tail.  The label “transport G1” is
therefore a model role, not an unconditional declaration that G1 is cloudy,
rising, or moist.

### 5.3 Coherent transport-plume structure blend

The moisture-primary construction is retained as the baseline.  A coherent
transport plume, however, is more generally a positive-$w$ tail with positive
$w$-$r_t$ association.  Its $w$-$\theta_l$ association evolves: it can be
positive in a fresh warm/moist plume and negative after the plume rises into
warmer environmental air.  That signed information is already present in the
supplied **second** moments.  Define the normalized correlations

$$
\rho_{wr_t}=R_{wr_t},
\qquad
\rho_{w\theta_l}=R_{w\theta_l}.
\qquad (6a)
$$

trivariate transport PDF uses a gate that identifies the common transport signature without
requiring a fixed $w$-$\theta_l$ sign:

$$
G=\min\left[1,
\frac{
\max(\widetilde S_w,0)
\max(\rho_{wr_t},0)
}{r_g^2}\right],
\qquad r_g=0.25,
\qquad \eta=f_{\mathrm{struct}}G.
\qquad (6b)
$$

Here $f_{\mathrm{struct}}=$ `trivar_plume_structure_strength`, bounded to
$[0,1]$.  The gate is local rather than a prescribed altitude switch.  It is
nonzero for a positive-$w$ skewed, moisture-transporting tail whether
$\rho_{w\theta_l}$ is positive or negative.  Thus it represents the same
plume family through its fresh and mature thermodynamic stages.

For an active gate, define a signed plume direction and a
$w$-skewness-derived small-component weight:

$$
\boldsymbol q_{\mathrm{plume}}=
\begin{bmatrix}
1\\ \rho_{wr_t}\\ \rho_{w\theta_l}
\end{bmatrix},
\qquad
a_w=\mathrm{clip}\left[
\frac12\left(1-
\frac{|S_w^{\mathrm{cap}}|}{\sqrt{4+(S_w^{\mathrm{cap}})^2}}
\right),\ 0.035,\ 0.5\right].
\qquad (6c)
$$

If $a_{r_t}$ denotes Eq. (5)'s moisture-tail weight, the actual weight and
pre-normalization direction become

$$
a=(1-\eta)a_{r_t}+\eta a_w,
\qquad
\boldsymbol q_0=(1-\eta)
\begin{bmatrix}
g_w\widetilde S_w\\ \widetilde S_{r_t}\\ \widetilde S_{\theta_l}
\end{bmatrix}
+\eta\boldsymbol q_{\mathrm{plume}}.
\qquad (6d)
$$

The $w$-weight expression is the two-point mixture relation for a specified
$w$ skewness.  It is used only as a bounded transport-plume proxy; the
subsequent width contrast still determines how much of each marginal third
moment can be realized.  This construction introduces neither a mixed
third-moment input nor an iterative fit.

The center-separation signal is blended consistently:

$$
\chi=(1-\eta)|\widetilde S_{r_t}|
+\eta\max\left(|\widetilde S_{r_t}|,|\widetilde S_w|\right).
\qquad (6e)
$$

### 5.4 A joint, covariance-aware direction

The raw standardized direction is

$$
\boldsymbol{q}_0=
\begin{bmatrix}
g_w\widetilde{S}_w\\
\widetilde{S}_{r_t}\\
\widetilde{S}_{\theta_l}
\end{bmatrix},
\qquad (6)
$$

where $g_w$ is `trivar_w_direction_scale`.  This is the key moisture-oriented
choice: $r_t$ and $\theta_l$ enter directly, while the $w$ contribution
is deliberately attenuated.  With the default $g_w=0.20$, a scalar tail can
be far separated without demanding an equally large mean-$w$ offset.

When $\boldsymbol{q}_0$ is non-negligible, trivariate transport PDF normalizes it in the
Mahalanobis metric of the grid correlation matrix:

$$
\boldsymbol{q}=
\frac{\boldsymbol{q}_0}
{\sqrt{\boldsymbol{q}_0^\mathsf{T}\mathbf{R}^{-1}\boldsymbol{q}_0}}.
\qquad (7)
$$

Consequently $\boldsymbol{q}^\mathsf{T}\mathbf{R}^{-1}\boldsymbol{q}=1$.
The direction respects the already predicted correlation structure: a
displacement along a strongly correlated thermodynamic direction consumes less
independent covariance budget than the same Euclidean displacement across it.
The code computes this using a 3-by-3 Cholesky solve.

If the raw direction is effectively zero, trivariate transport PDF sets
$\boldsymbol{q}=\boldsymbol{0}$ and the component centers coincide.  This
is the symmetric limit; marginal skewness cannot then be reconstructed through
the two-Gaussian width formula because it requires a nonzero mean separation.

### 5.5 Covariance-budgeted separation

The standardized component separation is

$$
\boldsymbol{d}
=\frac{c\,\chi}{\sqrt{ab}}\boldsymbol{q},
\qquad (8)
$$

where $c$ is `trivar_center_budget`, bounded to
$[0.02,0.97]$.  Substitution into Eq. (1) gives

$$
ab\boldsymbol{d}^\mathsf{T}\mathbf{R}^{-1}\boldsymbol{d}
=c^2\chi^2<1.
\qquad (9)
$$

Thus, the between-component covariance consumes a controlled fraction of the
covariance-metric room inside $\mathbf{R}$.  The moisture signal, or the
gated coherent-plume signal in the structure branch, makes the actual
fraction smaller than $c^2$ near symmetric states.  This is why a strong
center displacement disappears smoothly instead of being imposed solely by a
parameter.

Equation (9) guarantees that the residual in Eq. (2) is positive semidefinite
in exact arithmetic.  The code nevertheless performs a Cholesky check, both
for numerical safety and to reject ill-conditioned inputs.

---

## 6. Component widths, internal covariance, and marginal skewness

### 6.1 The diagonal third-moment identity

For any one standardized variable $x_j$, let
$V_{1j}$ and $V_{2j}$ be its two component variances, and let $d_j$
be the corresponding center separation.  A two-Gaussian mixture obeys

$$
S_j
=ab\left[(b-a)d_j^3+3d_j(V_{1j}-V_{2j})\right].
\qquad (10)
$$

For nonzero $d_j$, the requested variance difference is therefore

$$
\Delta V_j\equiv V_{1j}-V_{2j}
=\frac{
S_j^{\mathrm{cap}}/(ab)-(b-a)d_j^3
}{3d_j}.
\qquad (11)
$$

This is an exact algebraic identity, not an iterative solve.  It provides the
important extra degree of freedom that ADG1's equal or tightly constrained
component widths do not provide: G1 can be broader in $w$ even when
$d_w$ is modest.

A zero $d_j$ cannot carry a nonzero third moment through a width difference
alone.  A mixture with identical centers is symmetric in that marginal even if
its component widths differ.  trivariate transport PDF consequently makes no width request for
that variable when $|d_j|$ is below a small numerical threshold.

### 6.2 Diagonal contrast construction

Let the base covariance of each component be temporarily
$\mathbf{R}_{\mathrm{res}}$.  Define a contrast matrix $\mathbf{A}$
whose diagonal entries are

$$
A_{jj}=b\Delta V_j.
\qquad (12)
$$

and construct a contrast path

$$
\boldsymbol{\Sigma}_1(t)=\boldsymbol{\Sigma}_{1,0}+t\mathbf{A},
\qquad
\boldsymbol{\Sigma}_2(t)=\boldsymbol{\Sigma}_{2,0}
-\frac{a}{b}t\mathbf{A}.
\qquad (13)
$$

For every $t$,

$$
a\boldsymbol{\Sigma}_1(t)+b\boldsymbol{\Sigma}_2(t)
=a\boldsymbol{\Sigma}_{1,0}+b\boldsymbol{\Sigma}_{2,0}.
$$

The diagonal part of Eq. (13) asks for the marginal third moments without
altering the conserved mixture covariance.  trivariate transport PDF first clips each requested
$\Delta V_j$ to retain a small positive marginal-variance margin.  It then
applies a full-matrix PSD cap (Section 8).  Hence the diagonal third moments
are exactly reproduced only when the requested contrast survives both caps.

This distinction is essential in model evaluation:

- **Mean and covariance error** from the closure construction should be at
  numerical roundoff for a non-fallback point.
- **Third-moment error** can be structural and should be diagnosed as an
  active realizability constraint, not automatically treated as a coding
  error.

### 6.3 Positive residual $w$-$r_t$ covariance in G1

For a positive residual $R_{\mathrm{res},wr_t}>0$, trivariate transport PDF asks G1 to
capture a fraction $f_{wr_t}$ of the residual grid covariance through the
off-diagonal contrast

$$
A_{wr_t}=A_{r_tw}
=f_{wr_t}\frac{b}{a}R_{\mathrm{res},wr_t},
\qquad 0\le f_{wr_t}\le1.
\qquad (14)
$$

Using Eq. (13), if this request is fully realizable and $f_{wr_t}=1$, G2
has zero internal residual $w$-$r_t$ covariance and G1 carries all of it:

$$
a\Sigma_{1,wr_t}=R_{\mathrm{res},wr_t},
\qquad
\Sigma_{2,wr_t}=0.
$$

This is the mathematical expression of the central transport hypothesis:

> After allowing a moisture-oriented separation of centers, positive residual
> $w$-$r_t$ covariance belongs preferentially inside the transport
> component.

The user-facing control is `trivar_g1_wrt_capture`; its default is one.
The final PSD cap can reduce the requested allocation.

### 6.4 Negative grid $w$-$r_t$ covariance

The intended plume interpretation is specifically a positive G1 internal
tilt.  When the grid correlation is negative, trivariate transport PDF follows a separate,
deliberate convention: it attempts to force the **G1 residual**
$w$-$r_t$ covariance to zero, leaving the negative residual covariance in
G2.  In standardized form, before width and positive-tilt contrasts,

$$
\Sigma_{1,wr_t}=0,
\qquad
\Sigma_{2,wr_t}=
\frac{R_{\mathrm{res},wr_t}}{b}.
\qquad (15)
$$

The requested reallocation itself is subjected to the same PSD-safe scaling.
Therefore a rare incompatible state may retain a small G1 negative tilt, but
the algorithm aims for zero rather than allocating a negative internal tilt to
the transport component.

### 6.5 Coherent-plume thermodynamic covariance allocation

When the structure gate in Eq. (6b) is active, the requested G1 capture of a
positive residual $w$-$r_t$ covariance is at least $\eta$; its effective
capture fraction is $\max(f_{wr_t},\eta)$.  The same branch uses the residual
component-2 covariance before the final contrast, denoted
$\boldsymbol\Sigma_{2,0}$, to request

$$
A_{w\theta_l}=\eta\frac{b}{a}\Sigma_{2,0,w\theta_l},
\qquad
A_{r_t\theta_l}=\eta\frac{b}{a}
\max\left(\Sigma_{2,0,r_t\theta_l},0\right).
\qquad (15a)
$$

By Eq. (13), a fully realizable request gives

$$
\Sigma_{2,w\theta_l}=(1-\eta)\Sigma_{2,0,w\theta_l}.
\qquad (15b)
$$

Thus G2's $w$-$\theta_l$ tilt moves toward zero regardless of sign.  Positive
G2 $r_t$-$\theta_l$ tilt also moves toward zero, while already-negative G2
$r_t$-$\theta_l$ covariance is retained.  The weighted opposite contrasts in
G1 preserve the supplied covariance exactly.  This encodes a simple division
of roles: G1 carries the coherent signed plume structure; G2 can retain a
nonpositive thermodynamic background relation but is not asked to represent an
internally warm/moist transport plume.

These requests and the diagonal width contrasts enter the same fixed-size PSD
cap.  The preference is therefore algebraic rather than hard: incompatible
moments reduce the entire contrast continuously, with no iterative fit.

### 6.6 Mixed-third moments: evidence for the allocation, not yet a replacement

The coherent-plume branch above is a deliberately low-order prior.  It uses
the signs and geometry already available from the predicted second moments to
prefer a structured G1 and a comparatively well-mixed G2.  It does **not**
know whether energetic motion and moisture departures actually occur in the
same parcels.  The resolved mixed third moments
`WP2RTP = <w'^2 r_t'>` and `WPRTPTHLP = <w' r_t' theta_l'>` supply precisely
that missing conditional information in SAM.

For fixed weights, centers, and component variances, either quantity is
linear in the three internal off-diagonal covariance entries.  For example, a
mean-preserving G1/G2 contrast in `C_wr` changes `WP2RTP` in proportion to
the difference between the two component w offsets.  However, the low-level
ARM V exposes that the already diagnosed width/tilt contrast often exhausts
the available PSD room.  The current Notes laboratory therefore uses the
mixed values first as *center-direction compatibility information*: it takes
one bounded local solve that rotates the r_t and theta_l components of the
center direction at fixed covariance-metric length, then rebuilds the normal
trivariate transport PDF width and covariance allocation.  `WPRTP2` is intentionally held out
as an independent check.

This is a useful diagnostic experiment but not a production replacement for
the plume-structure prior:

- It holds the center-separation metric length fixed, but rotates its
  thermodynamic direction.  A weak or counterproductive result says that the
  error cannot be repaired by center orientation alone.
- The raw mixed moments are not yet reliable prognostic CLUBB inputs.  Using
  them directly in a free-running closure would make an LES comparison
  circular.
- The finite-difference compatibility step is intentionally one-shot and
  bounded.  It compares a fixed small set of fractions of that one step and
  retains only a fraction that lowers the combined two-target residual across
  the ordinary reconstruction's PSD/cap branches.  It is a laboratory
  diagnostic of the missing conditional information, not a hidden iterative
  fit or a new production solve.

If the predictor experiment shows consistent benefit across cases, the more
principled successor is a bounded analytic direction predictor, using mixed
moment proxies to regularize the existing coherent-plume direction before the
ordinary trivariate transport PDF reconstruction.  A secondary joint covariance allocation can
then be tested only where the rebuilt component widths leave PSD room.  Either
path needs externally supplied or independently predicted mixed moments before
it replaces `trivar_plume_structure_strength` in CLUBB.

---

## 7. Recovering physical component statistics

After the standardized means and covariances are diagnosed, the physical
component state is

$$
\boldsymbol{\mu}_1
=\overline{\boldsymbol{x}}+b\mathbf{D}\boldsymbol{d},
\qquad
\boldsymbol{\mu}_2
=\overline{\boldsymbol{x}}-a\mathbf{D}\boldsymbol{d},
\qquad (16)
$$

$$
\mathbf{C}_1=\mathbf{D}\boldsymbol{\Sigma}_1\mathbf{D},
\qquad
\mathbf{C}_2=\mathbf{D}\boldsymbol{\Sigma}_2\mathbf{D}.
\qquad (17)
$$

The reported component correlations are calculated from the entries of
$\mathbf{C}_1$ and $\mathbf{C}_2$.  Importantly, these G1 and G2
correlations are not overwritten by CLUBB's historical equal-correlation
reconstruction.  That preservation is necessary for trivariate transport PDF's internal tilt
to have any downstream effect.

For a successful diagnosis, substitution of Eqs. (16) and (17) into the
mixture identities gives

$$
a\boldsymbol{\mu}_1+b\boldsymbol{\mu}_2
=\overline{\boldsymbol{x}},
$$

$$
a\left[\mathbf{C}_1+
(\boldsymbol{\mu}_1-\overline{\boldsymbol{x}})
(\boldsymbol{\mu}_1-\overline{\boldsymbol{x}})^\mathsf{T}\right]
+b\left[\mathbf{C}_2+
(\boldsymbol{\mu}_2-\overline{\boldsymbol{x}})
(\boldsymbol{\mu}_2-\overline{\boldsymbol{x}})^\mathsf{T}\right]
=\mathbf{C}.
\qquad (18)
$$

Equation (18) is one of trivariate transport PDF's strongest implementation invariants and is
checked by the dedicated unit test.

---

## 8. Realizability: direct positive-definite caps without iteration

### 8.1 Why a full-matrix constraint is required

Positive marginal variances are insufficient.  A three-variable covariance
matrix may have all positive diagonal entries but still be physically
impossible because its correlations are mutually inconsistent.  trivariate transport PDF must
therefore keep **both** component covariance matrices positive definite.

The challenge is that the diagonal width contrast and G1 $w$-$r_t$ tilt
are requested together.  A scalar correlation clip would not preserve their
joint geometry or the mixture covariance.  The implementation instead scales
the entire requested contrast by a single analytically determined factor.

### 8.2 Whitening-based contrast cap

Consider a component covariance path

$$
\mathbf{C}(t)=\mathbf{C}_0+t\mathbf{A},
$$

where $\mathbf{C}_0$ is positive definite.  Let
$\mathbf{C}_0=\mathbf{L}\mathbf{L}^\mathsf{T}$ be its Cholesky factorization.
Then

$$
\mathbf{C}(t)\succ0
\quad\Longleftrightarrow\quad
\mathbf{I}+t\mathbf{L}^{-1}\mathbf{A}\mathbf{L}^{-\mathsf{T}}\succ0.
\qquad (19)
$$

The eigenvalues of the whitened contrast determine the largest admissible
$t$.  If its most negative eigenvalue is $\lambda_-<0$, the boundary is

$$
t< -\frac{1}{\lambda_-}.
\qquad (20)
$$

trivariate transport PDF evaluates this constraint for **both** paths in Eq. (13), chooses the
most restrictive value, and retains a small interior margin.  The 3-by-3
symmetric eigenvalues are evaluated with a closed-form trigonometric formula;
no iterative eigensolver or nonlinear root finder is used.

There are two applications of this cap:

1. the negative-covariance allocation that seeks a zero G1
   $w$-$r_t$ tilt; and
2. the combined diagonal width plus positive G1-tilt contrast.

The second cap preserves the relative allocation among the requested width
differences and positive tilt, but it may reduce all of them together.  This
is a conscious preference for smooth, conservative behavior over an exact
third-moment fit at the edge of the realizable set.

### 8.3 Bounds and numerical margins

Several layered safeguards keep the closure in the interior of the feasible
set.

| Safeguard | Purpose |
| --- | --- |
| Positive input-variance floor | Avoids division by a nearly zero standard deviation. |
| Weight floor $a\ge0.035$ | Prevents singular rare-component factors. |
| Center budget $c\le0.97$ | Keeps the residual covariance away from its rank-deficient boundary. |
| Diagonal width bounds | Retain a small positive marginal variance before the full-matrix cap. |
| Cholesky tests | Reject non-positive-definite grid, residual, base, or final covariance matrices. |
| PSD interior margin | Avoids component covariances that are mathematically valid but numerically singular. |
| Collapsed fallback | Ensures a deterministic result even if the diagnosis cannot proceed. |

This hierarchy should be judged scientifically as well as numerically.  A cap
that activates frequently is not merely a stability success; it says the
diagnosed moment combination and the chosen two-Gaussian allocation are in
tension.  Future diagnostic statistics should record cap activation and
strength by regime.

---

## 9. Algorithm and fallback state

### 9.1 Pointwise algorithm

For each model column and level, trivariate transport PDF performs the following sequence.

1. Form the grid mean $\overline{\boldsymbol{x}}$, variance vector, full
   covariance $\mathbf{C}$, and skewness vector.
2. Require positive marginal variances; standardize to $\mathbf{R}$ and
   require a Cholesky factorization.
3. Soft-cap the three skewnesses using Eq. (4).
4. Diagnose the G1 weight from moisture skewness using Eq. (5).
5. Form and covariance-normalize the trivariate direction using Eqs. (6)--(7).
6. Form the budgeted center separation using Eq. (8) and check the residual
   correlation in Eq. (2).
7. Construct the negative-tilt base allocation when the grid
   $w$-$r_t$ correlation is negative; PSD-cap it if needed.
8. Diagnose the three requested diagonal width differences from Eq. (11),
   apply marginal bounds, and add the requested positive G1 $w$-$r_t$
   capture when applicable.
9. Apply the whitening/eigenvalue cap to the complete contrast, producing two
   positive-definite standardized component covariances.
10. Transform centers and covariances to physical units and report their
    variances and correlations to the existing PDF closure.

The computational work is fixed size: a few 3-by-3 matrix products, Cholesky
factorizations, triangular solves, and closed-form eigenvalue evaluations.
There is no vertical coupling, root solve, iterative projection, or sampling
step.

### 9.2 Collapsed fallback

If any required state is invalid or a safety check fails, trivariate transport PDF returns

$$
a=b=\tfrac12,
\qquad
\boldsymbol{\mu}_1=\boldsymbol{\mu}_2=\overline{\boldsymbol{x}},
\qquad
\mathbf{C}_1=\mathbf{C}_2=\mathbf{C}_{\mathrm{safe}}.
\qquad (21)
$$

For a valid supplied covariance,
$\mathbf{C}_{\mathrm{safe}}=\mathbf{C}$.  If the supplied covariance is
itself not positive definite, the fallback uses its positive diagonal
variances and zero cross-covariances.  This fallback conserves the mean and,
when the input covariance was valid, the covariance, but it deliberately
removes all skewness and component distinction.

The fallback is preferable to a local numerical failure because CLUBB can then
continue safely.  It must nevertheless be monitored in LES and production
experiments: frequent fallback would invalidate the scientific interpretation
of the PDF in that regime.

---

## 10. Relationship to ADG1

trivariate transport PDF should not be described as “ADG1 with one extra correlation.”  It makes
three coupled structural changes.

| Aspect | ADG1 tendency | trivariate transport PDF tendency |
| --- | --- | --- |
| Primary asymmetry | $w$-anchored through $w'^3$. | Moisture-primary rare-component mass, with all three skewnesses setting one direction. |
| Center diagnosis | $w$ separation first; scalar response follows. | One trivariate separation vector diagnosed jointly. |
| $w$-scalar covariance | Within-component values are zero. | Positive residual $w$-$r_t$ covariance may reside in G1. |
| Width allocation | Simplicity-oriented component-width structure. | Independent algebraic width contrasts requested in all three marginals, then jointly capped. |
| Thermodynamic correlations | Historical component-correlation reconstruction. | Full G1/G2 thermodynamic covariance matrices diagnosed and retained. |
| Intended geometry | Distinct rising/sinking or vertically separated components. | A moisture-separated, internally tilted transport lobe with modest mean-$w$ separation. |

Both methods remain two-Gaussian closures and therefore share basic limits:
they cannot represent arbitrary multimodality, strong non-Gaussian tails, or
all independent third and fourth moments.  trivariate transport PDF simply places its limited
freedom in a different part of the geometry.

---

## 11. Connection to CLUBB cloud and transport diagnostics

### 11.1 Existing analytic integration is retained

trivariate transport PDF is a new component-state diagnosis inside
`pdf_closure_module`; it is not a new cloud-integration method.  The existing
CLUBB machinery transforms each component from $(w,r_t,\theta_l)$ to the
linearized thermodynamic saturation variables and evaluates Gaussian cloud
fraction, cloud water, and moments analytically.

The preservation of distinct component correlations is crucial here.  For a
component saturation variable that is a local linear combination of
thermodynamic fluctuations,

$$
\chi'_i=c_{r_t,i}r_{t,i}'+c_{\theta_l,i}\theta_{l,i}',
$$

its within-component velocity covariance is

$$
\mathrm{cov}_i(w,\chi)
=c_{r_t,i}\mathrm{cov}_i(w,r_t)
+c_{\theta_l,i}\mathrm{cov}_i(w,\theta_l).
\qquad (22)
$$

The coefficient signs and values are set by CLUBB's saturation linearization.
Equation (22) shows why a positive internal $w$-$r_t$ tilt does not by
itself guarantee a larger cloud-water flux: the $w$-$\theta_l$ term can
reinforce or offset it.

### 11.2 Cloud-water flux includes the within-component covariance

For PDF types with nonzero component $w$-$\chi$ correlation, CLUBB's
component contribution to cloud-water flux contains both a between-mean and a
within-component term of the schematic form

$$
\left\langle w'r_c'\right\rangle_i
= (\mu_{w,i}-\overline w)(\mu_{r_c,i}-\overline{r_c})
+\mathrm{cov}_i(w,\chi)\,P_i(\chi>0).
\qquad (23)
$$

This is an important distinction from an axis-aligned ADG-like interpretation.
trivariate transport PDF's internal tilt has a direct route into diagnosed $w'r_c'$; it is not
discarded after the component contours are drawn.  Conversely, underprediction
of $w'r_c'$ cannot be diagnosed solely as missing tilt.  It can also result
from saturation linearization, cloud probability, a PSD-limited covariance,
or cancellation in Eq. (22).

### 11.3 Higher-order transport moments

trivariate transport PDF requests direct PDF integration of the $w$-scalar and
scalar-scalar higher moments used by the closure, including
$\overline{w'r_t'^2}$, $\overline{w'\theta_l'^2}$, and
$\overline{w'r_t'\theta_l'}$.  This ensures the PDF is available before
the relevant advancement terms even when the usual call placement would
otherwise defer it.  The change is infrastructural, but scientifically it
means the diagnosed tilted state participates in the moment tendencies rather
than being only an output diagnostic.

### 11.4 Variables outside the trivariate diagnosis

trivariate transport PDF currently diagnoses only the $(w,r_t,\theta_l)$ state.  Horizontal
velocity and passive-scalar component means and variances are held at their
grid values, with existing responder-correlation machinery used where
available.  It does not diagnose the special implicit horizontal-wind
responder coefficients required by the optional
`l_predict_upwp_vpwp = .true.` path.  CLUBB's configuration check therefore
requires that flag to be false for trivariate transport PDF.

This is not a small bookkeeping detail.  It defines the present domain of
validity: trivariate transport PDF is a trivariate thermodynamic transport closure, not yet a
complete replacement for every multivariate responder treatment in CLUBB.

---

## 12. Semi-implicit turbulent-advection responders

### 12.1 Status, purpose, and central assumption

trivariate transport PDF now supplies a deliberately limited set of local responders to CLUBB's
existing semi-implicit banded advancement code.  The first implementation
converts the turbulent advection of $w'^3$, $r_t'^2$, $\theta_l'^2$, and
$r_t'\theta_l'$; it leaves the $w'r_t'$ and $w'\theta_l'$ equations on their
established explicit PDF10 path.  This boundary is intentional: the converted
moments have nonsingular coefficient-plus-residual forms based on direct PDF
moments, whereas a robust responder for $\overline{w'^2x'}$ needs a separate
tangent diagnosis.

The central time-discretization assumption is modest but important:

> During one timestep, freeze the diagnosed trivariate transport PDF mixture geometry,
> realizability caps, and responder coefficients at time level $n$; advance
> the selected prognostic low-order moment at time level $n+1$.

Thus ``implicit'' here means *linearly semi-implicit*.  It does not mean that
the full PDF10 diagnosis, including its clipped weights, covariance caps, and
component Cholesky checks, is reevaluated inside a nonlinear solve.  The PDF
is diagnosed before the solve, supplies coefficients, and is rediagnosed on
the next model step.  This is the same practical distinction that permits
ADG1-like closures to take stable timesteps without an iterative PDF fit.

For a generic scalar $x$, the relevant turbulent-advection fluxes are

$$
\begin{aligned}
\frac{D\overline{w'^3}}{Dt}\bigg|_{\rm TA}
  &= -\frac{1}{\rho}\frac{\partial\rho\overline{w'^4}}{\partial z},\\
\frac{D\overline{w'x'}}{Dt}\bigg|_{\rm TA}
  &= -\frac{1}{\rho}\frac{\partial\rho\overline{w'^2x'}}{\partial z},\\
\frac{D\overline{x'^2}}{Dt}\bigg|_{\rm TA}
  &= -\frac{1}{\rho}\frac{\partial\rho\overline{w'x'^2}}{\partial z}.
\end{aligned}
\qquad (24)
$$

For the $(r_t,\theta_l)$ system, there is additionally

$$
\frac{D\overline{r_t'\theta_l'}}{Dt}\bigg|_{\rm TA}
=-\frac{1}{\rho}\frac{\partial
\rho\overline{w'r_t'\theta_l'}}{\partial z}.
\qquad (25)
$$

The code already has coefficient-plus-explicit-residual interfaces for every
flux in Eqs. (24)--(25).  The scientific task is therefore to decide how the
PDF10 geometry responds when one supplied low-order moment changes, not to
introduce a new vertical solver.

### 12.2 Why the original ADG1 formulas cannot simply be copied

Larson and Golaz (2005) obtain their compact formulas from a deliberately
restricted family: a common within-component $w$ variance fraction
$\widetilde\sigma_w^2$, coupled scalar-width assumptions, and a diagnostic
scalar-skewness ansatz.  For example, their fourth-moment relation is

$$
\overline{w'^4}
=a_3\overline{w'^2}^{,2}
+a_1\frac{\overline{w'^3}^{,2}}{\overline{w'^2}},
\qquad
a_1=\frac{1}{1-\widetilde\sigma_w^2},
$$

$$
a_3=3\widetilde\sigma_w^4
+6(1-\widetilde\sigma_w^2)\widetilde\sigma_w^2
+(1-\widetilde\sigma_w^2)^2.
\qquad (26)
$$

Their nonsingular scalar-variance transport formula is

$$
\overline{w'x'^2}
=\frac{\overline{w'^3}}
{(1-\widetilde\sigma_w^2)\overline{w'^2}}
\left[
\frac{\beta}{3}\overline{x'^2}
+\frac{1-\beta}{3(1-\widetilde\sigma_w^2)}
\frac{\overline{w'x'}^{,2}}{\overline{w'^2}}
\right].
\qquad (27)
$$

The first term inside the brackets is implicit in $\overline{x'^2}$;
the second is an explicit residual.  This is why the original formulation can
avoid explicit turbulent advection without predicting $\overline{x'^3}$.

trivariate transport PDF intentionally breaks several of those simplifying assumptions.  It
allows unequal component $w$ widths, a moisture-primary center direction,
component-specific $w$-$r_t$ and $w$-$\theta_l$ tilt, and a PSD-capped
allocation of $r_t$-$\theta_l$ covariance.  A single
$\widetilde\sigma_w^2$ or one universal $\beta$ therefore cannot in general
describe its actual higher moments.  Directly inserting Eq. (26) or (27)
would make the implicit transport correspond to a *different PDF* from the
one used for cloud and microphysical integration.

That mismatch is avoidable: PDF10 can obtain its responders from its own
diagnosed Gaussian components.

### 12.3 PDF10 fourth-moment responder for $w'^3$ transport

Let $\delta_{w,i}=\mu_{w,i}-\overline w$ and
$V_{w,i}=C_{i,ww}$.  The exact fourth moment of the diagnosed PDF10 mixture
is

$$
M_4 \equiv \overline{w'^4}
=\sum_{i=1}^{2} p_i
\left(\delta_{w,i}^4+6\delta_{w,i}^2V_{w,i}+3V_{w,i}^2\right),
\qquad (p_1,p_2)=(a,b).
\qquad (28)
$$

Define the bounded local shape factor

$$
\Gamma_4^n=
\max\!\left(
\frac{M_4^n}{(\overline{w'^2}^{,n})^2},0
\right)
\quad\hbox{for }\overline{w'^2}^{,n}>0,
\qquad
\Gamma_4^n=0\quad\hbox{otherwise}.
\qquad (29)
$$

The implemented PDF10 semi-implicit approximation is

$$
\overline{w'^4}^{,n+1}
\approx
\Gamma_4^n\,\overline{w'^2}^{,n}\,
\overline{w'^2}^{,n+1}.
\qquad (30)
$$

Equation (30) exactly recovers the diagnosed $M_4^n$ at the coefficient
state, preserves a positive fourth-moment coefficient, and uses the existing
new-PDF-style two-diagonal band contribution.  It is a one-step Picard
linearization of the actual PDF10 fourth moment, not an assertion that
PDF10 obeys Eq. (26).

There is a more accurate but less attractive option.  A tangent linearization
of the diagnosed mapping can retain the $w'^3$ response:

$$
M_4^{n+1}\approx M_4^n
+\left.\frac{\partial M_4}{\partial\overline{w'^2}}\right|_n
\Delta\overline{w'^2}
+\left.\frac{\partial M_4}{\partial\overline{w'^3}}\right|_n
\Delta\overline{w'^3}.
\qquad (31)
$$

This produces the familiar coupled five-diagonal $w'^2$--$w'^3$ system used
by ADG1.  It remains a linear solve when the derivatives are frozen, but the
PDF10 mapping contains floors and PSD caps and is only piecewise smooth.
Equation (30) is therefore the appropriate first implementation; Eq. (31)
should be considered only if timestep-sensitivity tests show that the frozen
shape factor is inadequate.

### 12.4 Scalar-flux and scalar-variance responders

For each $x\in\{r_t,\theta_l\}$, write the desired PDF10 responder forms as

$$
\overline{w'^2x'}
=A_{w^2x}\,\overline{w'x'}+B_{w^2x},
\qquad (32)
$$

$$
\overline{w'x'^2}
=A_{wx^2}\,\overline{x'^2}+B_{wx^2}.
\qquad (33)
$$

Equation (32) is the PDF10 analogue of Larson and Golaz's simple
$\overline{w'^2x'}\propto\overline{w'x'}$ result.  Equation (33) is the
analogue of Eq. (27).  In both cases the coefficients and residual are
evaluated from the **same** two component means and covariance matrices used
by PDF10's cloud calculation.

The implemented Eq. (33) responder uses the frozen-normalized-geometry
secant.  Holding mixture weight, standardized center direction, width
contrast, and component correlations fixed, a common rescaling of the scalar
fluctuation makes $\overline{w'x'^2}$ homogeneous in
$\overline{x'^2}$.  Hence, away from a variance floor,

$$
A_{wx^2}^{\rm raw}
=\frac{\overline{w'x'^2}^{,n}}
{\overline{x'^2}^{,n}},
\qquad
B_{wx^2}=\overline{w'x'^2}^{,n}
-A_{wx^2}\overline{x'^2}^{,n}.
\qquad (34)
$$

Here $A_{wx^2}$ is a bounded version of $A_{wx^2}^{\rm raw}$ and the second
identity defines the residual after that bound.  The code caps the magnitude
at $4\sqrt{\overline{w'^2}}$ and uses zero coefficient plus the full direct
moment as residual for nonpositive scalar variance.  Consequently the
diagnosed flux is reproduced at time $n$, extreme responder velocities can be
clipped, and a near-zero variance leaves the flux explicit rather than
singular.  This plays the same numerical role as Eq. (27), but does not
require a guessed scalar-skewness relation or a division by
$\overline{w'x'}$.

Equation (32) remains explicit in the first implementation.  Its preferred
future formulation is the local directional derivative
of the PDF10 moment map at frozen geometry,

$$
A_{w^2x}=
\left.\frac{\partial\overline{w'^2x'}}
{\partial\overline{w'x'}}\right|_{\rm frozen\ PDF10},
\qquad
B_{w^2x}=
\overline{w'^2x'}-A_{w^2x}\overline{w'x'}.
\qquad (35)
$$

The derivative can be derived analytically from the two-Gaussian moment
formula or evaluated by a small bounded local perturbation.  It must hold
the chosen center/covariance allocation rule fixed; otherwise the derivative
would accidentally treat a change in grid flux as a demand to rediagnose an
entirely different plume.  A regularized secant is an acceptable first
fallback when this derivative is poorly conditioned.

### 12.5 The mixed $r_t$--$\theta_l$ covariance transport is essential

Removing explicit transport from both scalar variances while retaining an
explicit mixed covariance would leave the thermodynamic covariance system
inconsistent.  PDF10 therefore also needs

$$
\overline{w'r_t'\theta_l'}
=A_{wr_t\theta_l}\,\overline{r_t'\theta_l'}
+B_{wr_t\theta_l}.
\qquad (36)
$$

The implemented mixed-covariance responder uses the same bounded secant
prescription rather than the derivative,

$$
A_{wr_t\theta_l}^{\rm raw}
=\frac{\overline{w'r_t'\theta_l'}^{,n}}
{\overline{r_t'\theta_l'}^{,n}},
\qquad
B_{wr_t\theta_l}=
\overline{w'r_t'\theta_l'}^{,n}
-A_{wr_t\theta_l}\overline{r_t'\theta_l'}^{,n}.
\qquad (37)
$$

The coefficient is limited by the same $4\sqrt{\overline{w'^2}}$ speed cap.
When $|\overline{r_t'\theta_l'}|$ is less than $10^{-6}$ times the product of
the scalar standard deviations, the code uses zero coefficient and retains
the full direct moment as the residual.  This avoids a sign-unstable division
near a vanishing mixed covariance while preserving the frozen-state flux.

A future tangent prescription would instead be

$$
A_{wr_t\theta_l}
=\left.\frac{\partial\overline{w'r_t'\theta_l'}}
{\partial\overline{r_t'\theta_l'}}\right|_{\rm frozen\ PDF10},
\qquad
B_{wr_t\theta_l}=
\overline{w'r_t'\theta_l}
-A_{wr_t\theta_l}\overline{r_t'\theta_l'}.
\qquad (38)
$$

This is where PDF10's trivariate formulation matters most.  The responder
must preserve the diagnosed division between center separation, G1 internal
tilt, and G2's well-mixed-background preference.  A scalar copy of the ADG1
coefficient would erase that distinction.  The existing CLUBB interface
already accepts precisely the $A_{wr_t\theta_l}$ and $B_{wr_t\theta_l}$
pair in Eq. (36).

### 12.6 Consequences, safeguards, and implementation boundary

The implemented responders make four assumptions that should be evaluated
rather than hidden:

| Assumption | Benefit | Consequence to test |
| --- | --- | --- |
| Freeze PDF10 geometry over one step | Retains a linear, local band solve. | Rapidly evolving plume regimes may lag by one step. |
| Bound responder velocities and retain the remainder in $B$ | Avoids singular/anti-diffusive coefficients. | A strongly capped state is partly explicit and should be diagnosed. |
| Use the actual PDF10 mixture moments | Transport, cloud integration, and plotted geometry describe one PDF. | Coefficients inherit piecewise behavior at PSD-cap transitions. |
| Use trivariate transport PDF only for the missing $\overline{w'x'^3}$ flux | Retains the existing scalar-third-moment state while adding conservative vertical transport. | The scalar-third-moment forcing covariance remains unclosed. |

No new global matrix type, quadrature, or nonlinear root solve is needed.
The existing vertical band solvers remain the authority.  PDF10's only new
responsibility is to return the local coefficients and explicit residuals,
with finite checks and a safe zero-coefficient fallback.  Exact agreement
between the explicit and semi-implicit tendencies is expected only at the
frozen coefficient state; stability, convergence with timestep, and retained
PDF realizability are the relevant validation criteria.

### 12.7 trivariate transport PDF scalar-third-moment transport advance

trivariate transport PDF also uses its two-Gaussian geometry to supply the previously omitted
fourth-order turbulent fluxes

$$
F_x = \overline{w'x'^3}, \qquad x\in\{r_t,\theta_l\}.
$$

For one Gaussian component, define $d_w=\mu_w-\overline{w}$,
$d_x=\mu_x-\overline{x}$, scalar variance $V_x$, and covariance $C_{wx}$.
Gaussian moment identities give

$$
F_{x,\mathrm{comp}}
= d_w\left(d_x^3+3d_xV_x\right)
+3C_{wx}\left(d_x^2+V_x\right).
\qquad (39)
$$

The mixture-weighted sum is exact for the frozen diagnosed PDF.  It is
evaluated from `pdf_params` on thermodynamic levels, where the ordinary
single-call PDF closure always stores the authoritative component state.  The
flux is then interpolated to momentum levels before CLUBB applies its
density-weighted divergence explicitly,

$$
T_{\rm ta}=-\rho_d^{-1}\frac{\partial(\rho_d F_x)}{\partial z},
\qquad (40)
$$

alongside the existing accumulation and turbulent-production terms.  Mean
advection, diffusion, and tau dissipation are implicit:

$$
\left[\frac{1}{\Delta t}+\frac{1}{\tau}
+\mathcal{A}_{\rm ma}+\mathcal{A}_{\rm diff}\right]x_3^{n+1}
=\frac{x_3^n}{\Delta t}+T_{\rm ta}^n+T_{\rm tp}^n+T_{\rm ac}^n.
\qquad (41)
$$

Here $x_3=\overline{x'^3}$ and
$K_{xp3}=\mathtt{c_K3}\,K_h$ with background diffusivity `nu3`.  Equation
(41) is one scalar tridiagonal system for $r_t'^3$ and one for
$\theta_l'^3$, using the same level ordering and solver wrapper as the other
thermodynamic-level CLUBB equations.  The top $x_3$ value remains fixed at
zero.  The diagnosed flux is set to zero at both boundaries because no surface
or top boundary closure for $\overline{w'x'^3}$ has been introduced.

This grid choice is intentional.  `pdf_params_zm` is an optional
momentum-level container that is populated only when the two-call PDF option
is active.  Reading it in the ordinary one-call configuration would turn the
new flux into a zero compatibility artifact rather than the diagnosed PDF
transport.  The thermodynamic-level diagnosis followed by the standard
interpolation makes the scalar-third-moment advance independent of that
optional execution mode.

For runtime verification, trivariate transport PDF writes the complete discrete tendency set to
the standard statistics stream: `rtp3_bt/ta/tp/ac/ma/df/dp/rs` and the
corresponding `thlp3_*` fields.  `rs` is the stored discrete time tendency
minus the named terms and should be roundoff-sized unless a future unbudgeted
adjustment is explicitly added.  In an anelastic SCM with $\overline{w}=0$,
`*_ma` is expected to be zero; that is a physical property of the configured
mean flow, not a missing budget calculation.

This is deliberately a first conservative stage.  The unresolved forcing term
$3\overline{x'^2F_x'}$ (e.g., microphysical or radiative covariance) is still
omitted, rather than approximated inconsistently.  The new tunable diffusion
pair has defaults `c_K3 = 0.025` and `nu3 = 1 m2 s-1`, matching the existing
scalar-variance scale until data support a different choice.

---

## 13. What the diagonal third moments do and do not identify

### 13.1 What they provide

The three marginal skewnesses provide three useful signals:

- $S_{r_t}$ identifies the strength of a moisture-tail asymmetry and sets
  the G1 mass.
- $S_{\theta_l}$ rotates the same transport direction thermodynamically.
- $S_w$ contributes a limited mean-$w$ direction and requests a $w$
  width contrast.

Together with the full covariance, these signals are enough to create a
plausible, realizable trivariate transport geometry without a nonlinear solve.

### 13.2 What they cannot determine

They do not uniquely determine:

- how much grid $w$-$r_t$ covariance should arise from center separation
  versus internal G1 tilt;
- which component should carry each marginal width contrast;
- the mixed third moments that distinguish a broad, moist, high-variance
  plume from another distribution with the same three marginal skewnesses;
- whether a rare moist tail is a cloud transport plume, an environmental
  remnant, or a different subgrid population.

trivariate transport PDF resolves these ambiguities by explicit modeling choices: moisture sets
the rare mass, a joint skewness direction sets the centers, and positive
residual $w$-$r_t$ covariance preferentially enters G1.  These choices
are hypotheses to be evaluated against LES, not consequences forced by the
available moments.

### 13.3 Most informative future inputs

If a later closure needs to diagnose rather than prescribe the allocation of
internal covariance and width, the natural candidates are mixed third moments
such as

$$
\overline{w'^2r_t'},\qquad
\overline{w'r_t'^2},\qquad
\overline{w'^2\theta_l'},\qquad
\overline{w'\theta_l'^2}.
$$

Of these, $\overline{w'^2r_t'}$ is particularly relevant to the motivating
case: it directly indicates whether large moisture anomalies preferentially
occupy broad $w$ distributions.  Adding such moments should be justified by
demonstrable LES information gain and by a simple, stable diagnostic relation;
it should not be added merely because it makes the PDF more flexible.

---

## 14. Tunable parameters and their interpretation

The five geometry controls are regular CLUBB tunable parameters (indices
103--107).  trivariate transport PDF's scalar-third-moment transport additionally uses the
diffusion pair at indices 113--114.  Their defaults and admissible bounds are
listed below.

| Parameter | Default | Allowed range | Physical/modeling role |
| --- | ---: | ---: | --- |
| `trivar_moisture_tail_gain` | 1.00 | 0--3 | Strength with which $|S_{r_t}|$ makes G1 rare. |
| `trivar_center_budget` | 0.72 | 0.02--0.97 | Maximum covariance-metric budget available to center separation. |
| `trivar_w_direction_scale` | 0.20 | 0--1.2 | Relative contribution of $S_w$ to the trivariate separation direction. |
| `trivar_g1_wrt_capture` | 1.00 | 0--1 | Requested fraction of positive residual $w$-$r_t$ covariance placed in G1. |
| `trivar_plume_structure_strength` | 0.00 compiled default; 0.50 default configuration | 0--1 | Strength of the gated, signed coherent-plume placement and G2-background covariance preference. |
| `c_K3` | 0.025 | nonnegative | Multiplier on $K_h$ in the trivariate transport PDF $r_t'^3$/$\theta_l'^3$ diffusion operator. |
| `nu3` | 1 m2 s-1 | nonnegative | Resolution-adjusted background diffusion in that same operator. |

These parameters should not be tuned as independent curve-fitting constants
without examining the diagnosed geometry.  Their principal interactions are:

- Larger moisture-tail gain makes G1 rarer, which amplifies the center offset
  implied by a fixed covariance budget.
- Larger center budget allocates more covariance to separation and less to
  within-component covariance.
- Larger $w$-direction scale turns the transport displacement more vertical;
  it can improve $w$ skewness fit but work against the target
  moisture-first geometry.
- Larger G1 capture transfers positive residual $w$-$r_t$ covariance from
  G2 to G1, subject to the joint PSD cap.
- Larger plume-structure strength has no effect unless the gate in Eq. (6b)
  detects positive-$w$, positive-$w$-$r_t$ transport.  When active, it moves
  G1 along the signed $w$-$\theta_l$ plume direction and asks G2 to retain a
  well-mixed $w$-$\theta_l$ structure and nonpositive $r_t$-$\theta_l$ tilt,
  subject to the same covariance-budget and PSD checks.

A sensible evaluation workflow is to diagnose distributions first, inspect
where caps activate, and only then tune parameters against statistically
independent LES cases.  A low scalar loss alone could conceal an unphysical
allocation of covariance between centers and component interiors.

---

## 15. Limiting behavior and interpretive checks

The following limits are useful for reviewing the method and for debugging
profiles.

| Situation | Expected trivariate transport PDF behavior |
| --- | --- |
| All three skewnesses near zero | $a\approx0.5$, the center direction vanishes, and the components approach the same centered Gaussian. |
| Large positive $S_{r_t}$ | G1 becomes rarer and its center moves toward positive $r_t$, with $w$ displacement limited by `w_direction_scale`. |
| Large negative $S_{r_t}$ | G1 becomes rarer and moves toward negative $r_t$; it is a dry-tail G1, not a moist plume. |
| Positive residual $w$-$r_t$ covariance | G1 can become internally positively tilted; at full capture G2 is requested to have no residual internal $w$-$r_t$ covariance. |
| Negative grid $w$-$r_t$ covariance | G1 is requested to be untilted; G2 carries the negative residual covariance. |
| Large third-moment request near a covariance boundary | The width/tilt contrast is reduced by PSD caps; mean and covariance remain exact, marginal skewness may not. |
| Invalid/ill-conditioned covariance | The closure returns the collapsed two-copy Gaussian fallback. |

These checks should be part of every LES diagnostic notebook or Notes-tab
comparison.  A plotted ellipse alone is insufficient: its component weight,
center contribution, internal covariance, and cap state must be shown
together.

---

## 16. Relationship to LES-advance and free-running tests

trivariate transport PDF depends directly on $S_{r_t}$ and $S_{\theta_l}$, unlike ADG1,
which is much more directly controlled by $w$ skewness.  This creates an
important distinction between two experiment types.

### 16.1 LES-advance experiments

When the optional `xp3` LES-advance group is enabled and the source contains
`RTP3` and `THLP3`, the scalar third moments can be constrained directly from
LES.  Such experiments are useful for isolating the geometry and cloud
diagnostics of trivariate transport PDF because the new scalar-skewness inputs are not themselves
being tested as a closure.

Without `xp3`, trivariate transport PDF remains defined: CLUBB supplies the existing diagnosed
scalar third moments.  The experiment is then less isolated, because a trivariate transport PDF
placement or flux error can arise from the scalar-third-moment diagnosis as
well as from the two-Gaussian geometry.

### 16.2 Free-running experiments

For free-running simulations, $r_t'^3$ and $\theta_l'^3$ are part of the
closure chain rather than merely diagnostic inputs.  For trivariate transport PDF,
`l_advance_xp3` now uses the conservative transport advance in Eq. (41),
including implicit mean advection, diffusion, and tau damping.  It is still
not a complete independently validated third-moment model because forcing
covariances are absent and the diagnosed turbulent flux is frozen explicitly
over each timestep.

The appropriate staged question is therefore:

1. Does trivariate transport PDF improve placement and cloud transport when the scalar third
   moments are LES constrained?
2. Does the existing scalar-third-moment path provide sufficient quality for
   the same improvement to survive free-running tests?
3. Only if not, is a more predictive scalar-third-moment treatment warranted?

This separation prevents PDF geometry from being blamed for an upstream
moment-closure error, or vice versa.

---

## 17. Current evidence and validation status

### 17.1 What has been verified

The present implementation has the following software and algebraic evidence:

- A dedicated Fortran unit test reconstructs the prescribed trivariate mean
  and full covariance, verifies both returned component covariance matrices
  are positive semidefinite/definite within tolerance, and checks unequal
  nonzero $w$-$r_t$ tilts for a positive-covariance test case.
- The unit test also checks the negative-covariance convention: G1 is not
  given an internal $w$-$r_t$ tilt when the grid covariance is negative.
- A Notes-tab Python laboratory applies the same closure family to raw SAM
  moment samples and overlays the component placement on $w$-$r_t$,
  $w$-$\chi$, and $w$-$r_c$ distributions.
- A short BOMEX SCM smoke run completed with
  `iiPDF_type=9,l_predict_upwp_vpwp=false`.

These checks establish basic algebraic safety and integration.  They do **not**
establish superior climate, cloud, or transport performance.

### 17.2 Claims that remain unproven

The following should be treated as hypotheses until demonstrated over an LES
suite:

- trivariate transport PDF places the relevant transport population better than ADG1 in the
  intended regimes.
- Its component covariance allocation improves cloud fraction, liquid water,
  and $w'r_c'$, not merely contour appearance.
- The four tunable parameters have stable values across cloud regimes and do
  not compensate for case-specific errors in scalar third moments.
- PSD caps and fallbacks are infrequent enough that the nominal closure,
  rather than its safety machinery, controls the answer.
- The free-running scalar-third-moment closure is accurate enough for the
  PDF's increased reliance on $S_{r_t}$ and $S_{\theta_l}$.

---

## 18. Recommended scientific evaluation plan

### Stage A — Diagnose the inputs before judging the PDF

For each LES time-height sample, compare CLUBB and LES for:

- all six entries of $\mathbf{C}$;
- $S_w$, $S_{r_t}$, and $S_{\theta_l}$;
- the derived weight $a$, separation vector $\boldsymbol{d}$, and
  between-component covariance fraction;
- the positive G1 $w$-$r_t$ capture and both PSD scales;
- fallback occurrence.

This stage answers whether the closure is receiving a plausible moment state.

### Stage B — Evaluate component geometry against raw LES/SAM samples

For representative regimes, overlay component contours and centers on raw
joint distributions in at least:

$$
(w,r_t),\qquad (w,\chi),\qquad (w,r_c).
$$

Assess separately:

- center placement;
- orientation and width of the transport component;
- weight of the transport component;
- whether apparent improvement survives the nonlinear cloud transform.

Conditioning the LES analysis on cloud, signed transport, or moisture-tail
samples can be informative, but the conditioning definition must be fixed
before scoring to avoid selecting a favorable geometry after the fact.

### Stage C — Evaluate closure-relevant diagnostics

Score at least cloud fraction, $r_c$, $w'r_c'$, $w'r_t'^2$,
$w'\theta_l'^2$, and $w'r_t'\theta_l'$, together with the basic moments.
The score should be stratified by cloud regime and height rather than relying
only on one column-integrated error metric.  trivariate transport PDF's reason for existing is
not a generic improvement in every statistic; it is a targeted improvement in
moist tilted transport structure and its cloud consequences.

### Stage D — Separate LES-constrained from free-running conclusions

Run the same samples in two modes:

1. moments constrained from LES where available, including `xp3`; and
2. fully free-running CLUBB moments.

The difference quantifies how much of the remaining error belongs to the PDF
versus its scalar-third-moment inputs.

### Stage E — Tune only after diagnosing cap statistics

Tune the five trivariate transport PDF parameters with training/test case separation.  Include
penalties or reporting for frequent PSD caps and fallback use.  A parameter
set that lowers a bulk error by repeatedly driving the closure against its
realizability cap should be viewed skeptically.

---

## 19. Key review questions before a full paper

The following questions provide a concise checklist for scientific review.

1. Is moisture skewness the right primary signal for a rare transport
   component across the target regimes, or should the weight depend jointly on
   a trivariate skewness measure?
2. Is the fixed attenuation of the $w$ direction physically justified, and
   can one parameter remain stable from boundary layer to elevated cloud?
3. Is the “positive residual $w$-$r_t$ covariance belongs in G1” rule
   supported by conditional LES statistics?
4. Are the three marginal third moments sufficient for the intended geometry,
   or does a mixed third moment measurably improve the center/width allocation
   enough to justify its prediction cost?
5. How often do center, negative-tilt, and final contrast caps activate in
   realistic CLUBB and LES states?
6. Does improved placement in $(w,r_t)$ lead to improved $w'r_c'$, or do
   saturation-variable covariance and cloud probability remain the dominant
   limitation?
7. Can the current scalar-third-moment treatment support free-running trivariate transport PDF,
   or is its error larger than the benefit of the new geometry?
8. What additional responder treatment is required before trivariate transport PDF can support
   `l_predict_upwp_vpwp = .true.`?

---

## 20. Bottom line

trivariate transport PDF is a deliberately simple but materially different two-Gaussian
closure.  It replaces a primarily $w$-anchored view of subgrid asymmetry
with a trivariate, moisture-oriented transport geometry.  Its main innovation
is not merely an added correlation parameter: it jointly diagnoses a
covariance-budgeted center vector, unequal component widths, and a protected
positive internal G1 $w$-$r_t$ tilt while preserving the full supplied
second-moment state.

The method is attractive because its mathematical complexity remains bounded:
it uses local algebra and fixed 3-by-3 linear algebra rather than a nonlinear
fit.  Its scientific risk is equally clear: diagonal third moments do not
uniquely identify the desired transport population, so the covariance
allocation rules are hypotheses.  The next step is therefore not to add
complexity reflexively, but to test those hypotheses against conditional LES
statistics, cloud diagnostics, cap frequencies, and free-running moment
quality.
