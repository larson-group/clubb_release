#!/usr/bin/env python3
"""lock in the tunable-parameter name list + index mirror (parameters_tunable.py)."""


from clubb_jax.src.CLUBB_core.parameters_tunable import get_param_names
from clubb_jax.src.CLUBB_core import parameter_indices

# Critical zero-based index constants → the PARAM_NAMES entry they must point at.
_INDEX_CONSTANTS = {
    "ibeta": "beta", "imu": "mu", "iSkw_denom_coef": "Skw_denom_coef",
    "ilambda0_stability_coef": "lambda0_stability_coef", "iC1": "C1", "iC2rt": "C2rt",
    "iC8": "C8", "iC11": "C11", "iC14": "C14", "igamma_coef": "gamma_coef",
}


def test_index_constant_self_consistency():
    """Each i<name> index constant maps PARAM_NAMES[i] == <name> (no oracle needed)."""
    names = list(get_param_names())
    bad = []
    for const, expected in _INDEX_CONSTANTS.items():
        idx = getattr(parameter_indices, const, None)
        if idx is None:
            bad.append(f"{const}: not defined in parameter_indices"); continue
        got = names[idx] if 0 <= idx < len(names) else f"<out of range {idx}>"
        if got != expected:
            bad.append(f"{const}={idx} -> PARAM_NAMES[{idx}]={got!r}, expected {expected!r}")
    assert not bad, "param-index constant(s) drifted from PARAM_NAMES:\n  " + "\n  ".join(bad)
    print(f"  {len(_INDEX_CONSTANTS)} index constants (ibeta/imu/iSkw_denom_coef/…) all map correctly  PASS")
