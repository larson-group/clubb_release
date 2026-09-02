"""Focused scheduling checks for ``Radiation/radiation_module.py``."""

from clubb_jax.src.CLUBB_core.jax_stats import JaxStats
from clubb_jax.src.advance_clubb_to_end import _advance_radiation
from clubb_jax.src.clubb_case_initalization import init_clubb_case


def test_sampled_radiation_stats_do_not_depend_on_radiation_cadence():
    """Fortran updates radiation stats on every sampled timestep."""
    state = init_clubb_case("input/case_setups/bomex_model.in")
    state["_jax_stats"] = JaxStats.empty(
        l_sample=True,
        names=("radht",),
        grids=("zt",),
        ncol=state["ngrdcol"],
        max_nlev=state["nzt"],
        grid_nlev=(state["nzt"], state["nzt"], 1, state["nzt"], 1, state["nzt"], state["nzt"]),
    )

    _advance_radiation(state, 0.0, l_rad_itime=False)

    assert sum(int(bank.sum()) for bank in state["_jax_stats"].nsamples) == state["nzt"]
