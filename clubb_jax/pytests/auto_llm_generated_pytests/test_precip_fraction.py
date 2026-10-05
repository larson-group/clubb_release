"""Focused synthetic precipitation fraction assertion contracts; real Rico checks live in tests/."""
import numpy as np
import jax

jax.config.update("jax_enable_x64", True)


from clubb_jax.pytests.microphysics_test_inputs import precip_fraction_from_fields
from clubb_jax.src.CLUBB_core.precipitation_fraction import precip_frac_assert_check


def test_precip_frac_assert_check():
    """precip_frac_assert_check accepts precip_fraction's own (self-consistent, in-range) output and rejects
    corrupted inputs — and serves as a cross-check that the JAX precip_fraction satisfies the Fortran assertions."""
    rng = np.random.default_rng(0)
    ng, nzt, hd = 1, 20, 2
    hydromet = np.abs(rng.normal(size=(ng, nzt, hd))) * 1e-4
    cf = rng.uniform(0, 1, (ng, nzt)); mf = rng.uniform(0.3, 0.7, (ng, nzt)); isf = np.zeros((ng, nzt))
    hmtol = np.array([1e-10, 1.9e-7])
    pf, pf1, pf2, pftol = precip_fraction_from_fields(
        hydromet, cf, cf, cf, isf, isf, isf, mf,
        np.array([True, False]), np.array([False, False]), hmtol, 1.0)
    args = (hydromet[0], hmtol, mf[0], pf[0], pf1[0], pf2[0], float(pftol[0]))
    assert precip_frac_assert_check(*args) is True, "valid precip_fraction output rejected"
    bad = np.array(pf[0]); bad[5] = 1.5            # precip_frac > 1
    assert precip_frac_assert_check(hydromet[0], hmtol, mf[0], bad, pf1[0], pf2[0], float(pftol[0])) is False
    bad2 = np.array(pf1[0]); bad2[3] += 0.5        # breaks the mixt_frac-weighted consistency
    assert precip_frac_assert_check(hydromet[0], hmtol, mf[0], pf[0], bad2, pf2[0], float(pftol[0])) is False
    print("  precip_frac_assert_check: valid PASS / precip_frac>1 + inconsistent FAIL  PASS")
