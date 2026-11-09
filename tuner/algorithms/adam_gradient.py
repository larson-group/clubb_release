"""Minimize a scalar objective with independent bounded Adam chains and supplied gradients.

AdamGradientStrategy proposes parameter sets through fill and accepts their
losses and derivatives through tell. For CLUBB, the scheduler combines BOMEX,
ATEX or other case derivatives into total_gradient before delivering feedback.
The strategy itself is objective-independent: an analytic test function can
supply the same total_loss/total_gradient fields without running CLUBB or JAX.
The catalog requires JAX for production model derivatives; no JAX code runs here.

Coordinates and equations
-------------------------
Logical coordinate x_j lies in [0,1], with physical p_j = l_j + s_j*x_j and
range width s_j = u_j-l_j. For a linked range, every physical target shares p_j.
range_gradient applies the chain rule: g_j = s_j * sum_k dF/dp_k, summing the
physical targets k belonging to logical range j. Fixed ranges have g_j=0.
A chain performs ordinary bias-corrected Adam at update t=1,...,T:

    m_t = beta1*m_(t-1) + (1-beta1)*g_t
    v_t = beta2*v_(t-1) + (1-beta2)*g_t**2
    mhat_t = m_t/(1-beta1**t), vhat_t = v_t/(1-beta2**t)
    x_t = clip(x_(t-1) - learning_rate*mhat_t/(sqrt(vhat_t)+epsilon), 0, 1).

Products, squares and divisions are coordinatewise. m_0=v_0=0. Defaults are
beta1=0.8, beta2=0.999 and epsilon=1e-8. The first moment therefore emphasizes
roughly the last five updates, which limits the influence of older directions
in short searches. Learning rate 0.01 is a step in the
normalized range, not 0.01 in every physical parameter's units. Projection
clips the new point; it does not reset accumulated moments at the boundary.
Nonfinite active derivatives, loss values or moments stop the method.
A finite derivative follows the executed numerical branch; it need not reveal
changes hidden by branch equalities or identify the global minimum.

Starts, feedback and completion
-------------------------------
The first chain starts at clipped configured defaults. Other chains use a
seeded Latin hypercube over active coordinates; the extra default point is
outside that stratified subset. prepare_options sets chain_count to workers
* batch_size. Direct construction can instead supply a chain count.
Each chain has at most one outstanding center. Sample IDs map feedback back
to the correct chain, so partial or reordered arrivals cannot mix Adam moments.
A chain must score its initial point and every updated point: T updates cost
T+1 candidates. All-fixed ranges produce only one candidate overall. Scoring
the final center is required before is_exhausted becomes true.

Usage
-----
Construct with names, defaults, ranges, seed and catalog options. fill queues
centers; tell consumes total_loss, total_gradient and gradient_parameters in
physical-name order. is_exhausted, estimated_sample_count and pending_batch_limit
report completion, candidate budget and a useful scheduling limit. Static
prepare_options/sample_count expose derived scheduling and budget to frontends.
_Chain stores x, moments, update index and pending sample ID for current-format
checkpoints. The scheduler separately retains the best scored point; an Adam
step can temporarily increase loss and never guarantees global convergence.

Reference: Kingma and Ba, Adam: A Method for Stochastic Optimization (2014),
arXiv:1412.6980. Parameter-bound projection is applied after each Adam update.
"""

from collections import deque
from dataclasses import dataclass
import math

import numpy as np

from tuner.gradients import range_gradient
from tuner.algorithms.base import BaseTuningStrategy, latin_hypercube_points


@dataclass
class _Chain:
    """Keep one chain's normalized center, moments and outstanding evaluation."""

    x: np.ndarray
    m: np.ndarray
    v: np.ndarray
    updates: int = 0
    pending: int | None = None
    done: bool = False


class AdamGradientStrategy(BaseTuningStrategy):
    """Propose bounded centers and update them using the final score's gradient."""

    def __init__(self, *, param_names, default_params_row, parameter_ranges,
                 max_updates=30, learning_rate=0.01, chain_count=1,
                 beta1=0.8, beta2=0.999, epsilon=1.e-8, seed=None):
        """Initialize zero moments and independent centers in normalized logical coordinates.

        Only positive-width ranges participate in Latin hypercube starts. Defaults
        come from each range's first physical target; base construction preserves
        links when proposals expand back into full parameter rows.
        """
        super().__init__(param_names=param_names, default_params_row=default_params_row,
                         parameter_ranges=parameter_ranges)
        self.max_updates = int(max_updates)
        self.learning_rate = float(learning_rate)
        self.beta1, self.beta2, self.epsilon = float(beta1), float(beta2), float(epsilon)
        self.minimum = np.asarray([spec["min"] for spec in self.parameter_ranges], dtype=float)
        self.span = np.asarray([spec["max"] - spec["min"] for spec in self.parameter_ranges], dtype=float)
        self.active = self.span > 0.0
        default = np.asarray([
            self.default_params_row[self.param_index[spec["targets"][0]]]
            for spec in self.parameter_ranges
        ])
        start = np.clip((default - self.minimum) / np.where(self.active, self.span, 1.0), 0.0, 1.0)
        start[~self.active] = 0.0
        starts = [start]
        if self.active.any() and int(chain_count) > 1:
            samples = latin_hypercube_points(np.random.default_rng(seed), int(chain_count) - 1,
                                             int(np.count_nonzero(self.active)))
            for sample in samples:
                point = np.zeros_like(start)
                point[self.active] = sample
                starts.append(point)
        self.chains = [_Chain(point, np.zeros_like(point), np.zeros_like(point)) for point in starts]
        self.pending_context = {}

    def fill(self, pending_samples: deque, capacity: int) -> None:
        """Queue one unscored center per ready chain until the pending queue reaches capacity.

        Record sample ID -> chain index before emission. A chain with an outstanding
        center cannot advance until tell returns its loss and gradient.
        """
        for index, chain in enumerate(self.chains):
            if len(pending_samples) >= capacity:
                break
            if chain.done or chain.pending is not None:
                continue
            physical = self.minimum + self.span * chain.x
            sample = self._sample_from_selected_params(dict(zip(
                (spec["name"] for spec in self.parameter_ranges), physical)))
            chain.pending = sample["sample_id"]
            self.pending_context[chain.pending] = index
            pending_samples.append(sample)

    def tell(self, completed_samples: list[dict]) -> None:
        """Apply one Adam update for each recognized completed center.

        Read the physical-name order from feedback, sum linked partials and scale by
        range widths before updating moments. After max_updates, the next completed
        center marks the chain done instead of creating an unscored final move.
        Unknown sample IDs are ignored; accepted IDs are removed from pending context.
        """
        for completed in completed_samples:
            sample_id = int(completed["sample_id"])
            index = self.pending_context.get(sample_id)
            if index is None:
                continue
            chain = self.chains[index]
            if not math.isfinite(float(completed["total_loss"])):
                raise ValueError("Adam requires a finite total_loss")
            gradient = range_gradient(completed["total_gradient"], completed["gradient_parameters"],
                                      self.parameter_ranges)
            del self.pending_context[sample_id]
            chain.pending = None
            if chain.updates >= self.max_updates or not self.active.any():
                chain.done = True
                continue
            chain.updates += 1
            chain.m = self.beta1 * chain.m + (1.0 - self.beta1) * gradient
            chain.v = self.beta2 * chain.v + (1.0 - self.beta2) * gradient**2
            if not np.isfinite(chain.m).all() or not np.isfinite(chain.v).all():
                raise ValueError("Adam gradient moments overflowed")
            m_hat = chain.m / (1.0 - self.beta1**chain.updates)
            v_hat = chain.v / (1.0 - self.beta2**chain.updates)
            chain.x = np.clip(chain.x - self.learning_rate * m_hat / (np.sqrt(v_hat) + self.epsilon), 0.0, 1.0)
            chain.x[~self.active] = 0.0

    def is_exhausted(self) -> bool:
        """Return true only after every chain has scored its final center.
        """
        return all(chain.done for chain in self.chains)

    def estimated_sample_count(self) -> int:
        """Count initial plus updated centers for each active chain, or one all-fixed center.
        """
        return len(self.chains) * (self.max_updates + 1 if self.active.any() else 1)

    def pending_batch_limit(self, batch_size: int) -> int:
        """Return enough batches to hold one outstanding center for every independent chain.
        """
        return max(1, math.ceil(len(self.chains) / batch_size))

    @staticmethod
    def prepare_options(options, request):
        """Set chain_count to workers times batch_size; saved counts do not override capacity.
        """
        options["chain_count"] = max(1, request["max_workers"] * request["batch_size"])

    @staticmethod
    def sample_count(options, ranges):
        """Report chain_count*(max_updates+1), or one candidate when every range is fixed.
        """
        if ranges and all(spec["min"] == spec["max"] for spec in ranges):
            return 1
        return options["chain_count"] * (options["max_updates"] + 1)
