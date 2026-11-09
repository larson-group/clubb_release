"""Combine cloud-case loss gradients into one direction for bounded parameter tuning.

A candidate can have separate loss derivatives from cases such as BOMEX and
ATEX. The optimizer needs derivatives of the same single score used to rank
that candidate. combine_case_gradients applies the chain rule to the selected
cross-case aggregation policy; range_gradient then converts the physical
parameter derivatives into the optimizer's linked search coordinates.
Neither function runs or differentiates a model.

Case aggregation
----------------
Let L_c(p) be a case loss and F(p)=A(L_1(p),...,L_C(p)) the final score.
The aggregator supplies a_c=dA/dL_c in case order. For P physical parameters,
gradients has shape (C,P), coefficients shape (C,), and the returned vector is

    dF/dp = sum_c a_c * dL_c/dp.

For example, coefficients [0.25,0.75] and BOMEX/ATEX partials [2,4] for C8
produce dF/dC8=3.5. Excluded cases have zero coefficient; even a NaN derivative
from such a case is omitted. Active nonfinite derivatives or coefficients fail
clearly. Coefficients must come from the actual ranking policy, not a separate
averaging rule. Sorting/ranking boundaries can be nonsmooth; this helper uses
the coefficients of the branch chosen by the aggregator without smoothing it.

Linked and normalized coordinates
---------------------------------
For logical range j with physical targets S_j, tying their values to q_j gives

    dF/dq_j = sum_(k in S_j) dF/dp_k.

When normalized=True, q_j=lower_j+(upper_j-lower_j)*x_j and x_j lies in [0,1],
so dF/dx_j=(upper_j-lower_j)*dF/dq_j. For linked C6rt=C6thl with partials 2
and 3 and range [0,4], the physical linked derivative is 5 and the normalized
derivative is 20. A fixed normalized range returns zero. normalized=False
returns the physical linked derivative and still requires its target partials.

Usage
-----
combine_case_gradients(coefficients, gradients) returns a physical vector.
range_gradient(gradient, gradient_parameters, parameter_ranges, normalized=True)
returns one value per logical range, in range order. Physical names may arrive
in any explicit order but must be distinct and complete for the used ranges.
Both functions check shapes and reject unusable active derivatives.
"""

import numpy as np


def combine_case_gradients(coefficients, gradients):
    """Combine one physical gradient per case using the final score's case coefficients.

    Input shapes are (C,) and (C,P); output is (P,). Zero-weight rows are omitted
    before checking derivative finiteness, so excluded-case NaNs cannot poison F.
    """
    coefficients = np.asarray(coefficients, dtype=float)
    gradients = np.asarray(gradients, dtype=float)
    if gradients.ndim != 2 or coefficients.shape != (gradients.shape[0],):
        raise ValueError("Case gradient shape does not match aggregation coefficients")
    active = coefficients != 0.0
    if not np.isfinite(coefficients).all() or not np.isfinite(gradients[active]).all():
        raise ValueError("Cannot use nonfinite gradients from an active case")
    result = np.sum(coefficients[active, None] * gradients[active], axis=0)
    if not np.isfinite(result).all():
        raise ValueError("Combined loss gradient is nonfinite")
    return result


def range_gradient(gradient, gradient_parameters, parameter_ranges, *, normalized=True):
    """Sum linked physical partials and optionally convert them to [0,1] range coordinates.

    Match gradient entries through gradient_parameters, then emit one value per
    parameter_ranges entry. Normalized fixed ranges contribute zero without
    requiring their unused derivatives. Other missing/nonfinite targets fail.
    """
    values = np.asarray(gradient, dtype=float)
    names = list(gradient_parameters)
    if values.shape != (len(names),) or len(names) != len(set(names)):
        raise ValueError("Loss gradient must match distinct physical parameter names")
    by_name = dict(zip(names, values))
    result = []
    for spec in parameter_ranges:
        width = float(spec["max"]) - float(spec["min"])
        if normalized and width == 0.0:
            result.append(0.0)
            continue
        targets = spec.get("targets", [spec["name"]])
        missing = set(targets) - set(by_name)
        if missing:
            raise ValueError("Missing loss gradient target(s): " + ", ".join(sorted(missing)))
        partials = np.asarray([by_name[target] for target in targets])
        if not np.isfinite(partials).all():
            raise ValueError(f"Nonfinite loss gradient for {spec['name']}")
        result.append(float(np.sum(partials)) * (width if normalized else 1.0))
    result = np.asarray(result)
    if not np.isfinite(result).all():
        raise ValueError("Tuning-coordinate loss gradient is nonfinite")
    return result
