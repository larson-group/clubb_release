# Choosing a tuning algorithm

The algorithms propose parameter sets for the [tuner](../README.md). Every method
uses the same case losses and ranking policy. Their differences are how they
choose points, whether they need derivatives, and how much work they request.
The [JSON catalog](algorithms.json) supplies their options to the CLI and Dash;
use `python run_scripts/run_tuner_job.py -list_algorithms` from the repository
root to see current defaults, types, categories and requirements. A `sampler`
chooses points without adapting to their losses; an `optimizer` uses losses or
gradients to guide subsequent proposals.

| Method / CLI name | Category | How points are chosen | Required model output | Main controls |
| --- | --- | --- | --- | --- |
| [Random](random.py), `random` | Sampler | Independent uniform samples across the bounds. | Scalar loss; Fortran or JAX. | `max_samples` (200). |
| [Latin hypercube](latin_hypercube.py), `latin_hypercube` | Sampler | A fixed stratified ensemble, covering each interval of every individual parameter range once. | Scalar loss; Fortran or JAX. | `num_points` (200). |
| [Resolve](resolve.py), `resolve` | Sampler | Every combination of evenly spaced coordinate values, including bounds. | Scalar loss; Fortran or JAX. | Optional `spacing`, default half the smallest nonzero parameter range width. |
| [SimAnn](simann.py), `simann` | Optimizer | Independent annealing chains; sometimes accept worse points, then cool and adapt steps. | Scalar loss; Fortran or JAX. | `max_iters` (200), temperature and step adaptation (`stp_adjst_shift` 0.7); optional `chain_count`. |
| [Adam (SPSA)](adam_spsa.py), `adam` | Optimizer | Estimate a direction with paired loss probes, then make bounded Adam updates. | Scalar loss; Fortran or JAX. | `max_updates` (30), `learning_rate`, `perturbation`, `spsa_pairs`. |
| [Adam (grad)](adam_gradient.py), `adam_gradient` | Optimizer | Apply bounded Adam updates using derivatives of the combined objective. | Physical-parameter loss gradient; JAX required for CLUBB. | `max_updates` (30), `learning_rate` (0.01), `beta1` (0.8), `beta2`, `epsilon`. |

## Coverage versus adaptive search

Random is a simple reference for comparing other searches. Latin hypercube is
useful when the goal is a space-covering ensemble of a chosen size. Its sampling
is stratified along each parameter individually; it does not visit every cell
of a multidimensional grid. Neither sampler changes proposals after seeing a
loss. The tuner still ranks their evaluated points and retains the best rows.
There is no separate adaptive Latin hypercube optimizer in this catalog.

Resolve gives deterministic grid coverage and can help inspect a small number
of parameters. Its evaluation count multiplies across dimensions, making a fine
grid expensive quickly. Its single spacing applies in physical units, which
matters when parameters have very different scales. Omit spacing (leave the Dash
field blank) to use half the smallest nonzero range width. Fixed parameters do
not constrain spacing and an entirely fixed search evaluates one point.

SimAnn can accept an uphill move, which helps explore beyond the current local
minimum. Temperature settings depend on the scale of the objective. Cooling,
finite budgets and the chosen starts can still leave it in a local minimum.

Both Adam methods search in coordinates scaled to each parameter's `[0,1]`
range and project updates back into those bounds. Gradient Adam uses the
executed model's local derivative; SPSA estimates a direction by changing all
active coordinates together. A useful derivative can guide local refinement,
but neither method guarantees decreasing loss at every step or a global minimum.
Nonsmooth branches, clipped probes and the learning rate can affect progress.
SPSA also has direction-sampling noise. Whether AD is faster than SPSA depends
on model differentiation cost, parameter count and the requested evaluation
budget; these implementations make no universal speed or convergence claim.

## Budgets, chains and batching

A **chain** is an independent search trajectory with its own current point and
adaptive state. Several chains expose work the tuner can evaluate concurrently.
They share the selected experiment and ranking policy, not optimizer moments.

- Random and LHS evaluate their requested number of points. Untuned parameters
  keep their defaults; linked parameters share one sampled coordinate.
- Resolve evaluates the complete Cartesian grid. Fixed coordinates contribute
  one grid value.
- SimAnn defaults to `workers × batch_size` chains, with an explicit chain
  count available. `max_iters` includes the initial evaluation; early temperature
  stopping can use fewer evaluations.
- Gradient Adam always uses `workers × batch_size` chains. Each chain scores
  its initial point and every updated center, so 30 updates cost 31 evaluations.
- SPSA Adam needs `2 × spsa_pairs` columns per chain's probe block. Batch size
  must be divisible by that number. Chains per batch are the quotient; concurrent
  batches are `ceil(workers / number_of_cases)`. A chain scores initial/final
  centers plus every paired probe. With two pairs and 30 updates it uses 122
  evaluations. Intermediate updated centers are not separately scored.

The first Adam chain starts at clipped defaults. Additional Adam chains use
Latin hypercube starts; SimAnn uses Latin hypercube starts for all chains.
The standalone LHS ensemble contains only stratified points and adds no default
point. The scheduler evaluates model baselines separately from candidate counts.
Both Adam methods reduce an all-fixed search to one evaluation.

Workers limit concurrent candidate-case evaluations; batch size is the number
of parameter columns in a model evaluation. Warm case processes can outnumber
active workers. Increasing either control can increase memory consumption,
especially with JAX gradients. Use a small worker/batch configuration when
trying a new experiment. The shared `-seed` controls reproducible starts and
proposals; it is not a separate setting for each algorithm.

## Assessing the methods

The generated [objective benchmarks](../auto_llm_generated_tests/pytests/unit_tests/test_algorithm_objectives.py)
run Goldstein–Price, five-dimensional Rastrigin, Himmelblau, Eggholder and
Schaffer from the [Fortran G_unit suite](../../src/G_unit_test_types/tuner_tests.F90)
through every registered method using default options. All methods must produce
finite bounded proposals, unique sample IDs and finite losses, make progress,
and finish within their declared budgets. For methods categorized as `optimizer`,
at least one trial must also find a loss strictly below that function's predefined
starting-point loss. Samplers need not improve. These checks use the same
proposal/feedback loop and catalog metadata, with no method-specific rules.
Global-minimum hit percentages are reported for comparison and do not decide
pass/fail. Each function uses the actual Fortran driver's global-hit rule.

From the repository root, run the quick default (ten trials per pair):

```bash
bash human_reviewed_tests/run_pytests.sh -unit -k algorithm_objectives -q --tb=short
```

For more useful success-rate estimates, increase the number of trials; narrow
`-k` to a method or function when exploring:

```bash
bash human_reviewed_tests/run_pytests.sh -unit -k algorithm_objectives -q --tb=short \
  --tuner-trials=100 --tuner-seed=7
bash human_reviewed_tests/run_pytests.sh -unit \
  -k 'test_algorithms_search_fortran_objectives and adam_gradient' -q --tb=short \
  --tuner-trials=100 --tuner-seed=7
bash human_reviewed_tests/run_pytests.sh -unit \
  -k 'test_algorithms_search_fortran_objectives and eggholder' -q --tb=short \
  --tuner-trials=20 --tuner-seed=100
```

The pytest plugin options use `--` (pytest reserves lowercase single-dash
options). `--tuner-trials` defaults to 10; `--tuner-seed` defaults to 7, with
consecutive seeds from that starting value. Ten trials resolve success rates
in steps of 10%, so larger trial counts are useful for manual assessment.
The terminal table appears with ordinary capture and `-q`; `-s` is unnecessary.
It shows each method's category, global hits/runs, global-hit and improvement
percentages, starting/best losses, median gap to the known minimum, average
candidate count and PASS/FAIL.
The exact coordinate/loss hit conditions are printed below the table. JUnit
properties retain per-trial points, losses, seeds and hit/improvement decisions.
For example, `Rastrigin / Adam (grad): 1/3` means one of three complete searches
found a point satisfying Rastrigin's hit rule. Each row covers one function and
one method; this fraction counts neither functions nor optimizer updates.

The harness supplies empty algorithm options, batch size 8, one worker and one
synthetic case. Gradients are calculated and supplied only when the catalog
requires them. Budgets come from each method's defaults and differ between
methods; increasing trials adds complete searches without enlarging each
search's budget.
These default settings are not calibrated to the Fortran ESA benchmark or its
10,000 trials. Failures to improve remain visible in the metrics; the harness
does not adjust algorithm settings or starting points to make them pass.
A method need not return the same coordinates as another: Himmelblau alone has
several equally good global minimizers, and its hit rule checks loss.

The suite also verifies copied objective values, analytic gradient references
and strict/inclusive hit-rule boundaries independently. Separate focused checks
cover Adam updates, linked coordinates and continuation. Real CLUBB loss and
gradient transport is checked by the tuner's Jenkins integration stage.

## Implementation and extension

Equations, state transitions, stopping rules and numerical limitations are
documented in each algorithm's `.py` file linked above. The
[shared base](base.py) defines the model-independent `fill`/`tell` contract;
[gradient helpers](../gradients.py) explain multi-case aggregation and the
chain rule for linked coordinates. See [adding an algorithm](../README.md#adding-an-algorithm)
for registration and feedback details. A new method needs a Python file and a
catalog entry; it does not need algorithm-specific Dash callbacks.
