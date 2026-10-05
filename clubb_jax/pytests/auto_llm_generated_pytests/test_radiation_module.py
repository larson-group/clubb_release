"""Focused scheduling checks for ``Radiation/radiation_module.py``."""

from dataclasses import replace

import jax
import jax.numpy as jnp
import numpy as np
import pytest

from clubb_jax.src import clubb_case_initalization
from clubb_jax.src.CLUBB_core.jax_stats import JaxStats
from clubb_jax.src.Radiation.simple_rad_module import simple_rad_lba
from clubb_jax.src.Radiation.soil_vegetation import advance_soil_veg
from clubb_jax.src.advance_clubb_to_end import _advance_radiation


RADIATION_FIELDS = (
    "radht", "Frad", "Frad_SW_up", "Frad_LW_up", "Frad_SW_down", "Frad_LW_down",
    "radht_SW", "radht_LW", "Frad_SW", "Frad_LW",
)
SOIL_FIELDS = ("deep_soil_T_in_K", "sfc_soil_T_in_K", "veg_T_in_K")


@pytest.fixture(scope="module", params=["arm", "bomex", "lba", "dycoms2_rf01"])
def radiation_state(request):
    read_namelist = clubb_case_initalization.read_namelist
    # Exercise radiation independently of unsupported microphysics/SILHS and
    # avoid opening case output files. Keep the real case's radiation settings.
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(
            clubb_case_initalization, "read_namelist",
            lambda path: dict(
                read_namelist(path), l_stats=False,
                microphys_scheme="none", lh_microphys_type="disabled",
            ),
        )
        state = clubb_case_initalization.init_clubb_case(
            f"input/case_setups/{request.param}_model.in"
        )
    state["_jax_stats"] = JaxStats.empty(
        l_sample=True,
        names=("radht", "Frad"),
        grids=("zt", "zm"),
        ncol=state["ngrdcol"],
        max_nlev=state["nzm"],
        grid_nlev=(state["nzt"], state["nzm"], 1, state["nzt"], 1, state["nzt"], state["nzt"]),
    )
    return state


@pytest.mark.parametrize("radiation_state", ["lba"], indirect=True)
def test_radiation_schedule_can_toggle_between_calls(radiation_state):
    """Static schedule variants still consume the current time and array values."""
    state = dict(radiation_state)
    previous = np.asarray(state["radht"]).copy()
    for step, advance in enumerate((True, False, True, False, True), start=1):
        time_current = state["time_initial"] + step * 600.0
        _advance_radiation(state, time_current, advance)
        if advance:
            expected = simple_rad_lba(
                state["gr"], state["ngrdcol"], time_current,
                state["time_initial"], state["radiation_parameters"],
            )
            np.testing.assert_allclose(state["radht"], expected, rtol=2e-6, atol=1e-12)
            assert np.any(np.asarray(state["radht"]) != previous)
        else:
            np.testing.assert_array_equal(state["radht"], previous)
        bank, slot = state["_jax_stats"].name_to_slot["radht"]
        np.testing.assert_array_equal(state["_jax_stats"].nsamples[bank][slot], step)
        previous = np.asarray(state["radht"]).copy()


@pytest.mark.parametrize("l_rad_itime", [False, True])
@pytest.mark.parametrize("l_sample", [False, True])
def test_radiation_schedule_preserves_outputs_and_stats(radiation_state, l_rad_itime, l_sample):
    """Compiled scheduling preserves retained values and samples every requested step."""
    state = dict(radiation_state)
    state["_jax_stats"] = state["_jax_stats"].begin_timestep(
        l_sample=l_sample, reset_accumulators=True,
    )
    # Nonzero, distinct values catch accidental clearing on a retained step.
    for index, name in enumerate(RADIATION_FIELDS, start=1):
        state[name] = jnp.full_like(state[name], index)
    before = dict(state)
    expected = dict(state)
    time_current = state["time_initial"] + 300.0
    with jax.disable_jit():
        _advance_radiation(expected, time_current, l_rad_itime)
    _advance_radiation(state, time_current, l_rad_itime)

    for name in (*RADIATION_FIELDS, *SOIL_FIELDS, "_jax_stats", "err_info"):
        for actual, reference in zip(
            jax.tree_util.tree_leaves(state[name]), jax.tree_util.tree_leaves(expected[name]),
        ):
            # Compiled float32 arithmetic can differ slightly from eager
            # evaluation (including fused operations on CPU).
            tolerance = 5e-6 if actual.dtype == jnp.float32 else 1e-12
            np.testing.assert_allclose(actual, reference, rtol=tolerance, atol=1e-12, equal_nan=False)
    if not l_rad_itime:
        for name in RADIATION_FIELDS:
            np.testing.assert_array_equal(state[name], before[name])
    for name, levels in (("radht", state["nzt"]), ("Frad", state["nzm"])):
        stats = state["_jax_stats"]
        bank, slot = stats.name_to_slot[name]
        np.testing.assert_array_equal(stats.nsamples[bank][slot], int(l_sample))
        np.testing.assert_array_equal(
            stats.buffers[bank][slot], state[name] if l_sample else np.zeros((state["ngrdcol"], levels)),
        )


@pytest.mark.parametrize("radiation_state", ["arm"], indirect=True)
@pytest.mark.parametrize("l_rad_itime", [False, True])
def test_soil_updates_use_incoming_fluxes_on_every_step(radiation_state, l_rad_itime):
    """Soil advances before radiation, including steps which retain radiation."""
    state = dict(radiation_state)
    state["radiation_parameters"] = replace(state["radiation_parameters"], l_soil_veg=True)
    for name, value in zip(SOIL_FIELDS, (288.58, 295.0, 300.0)):
        state[name] = jnp.full_like(state[name], value)
    for name, value in (("Frad_SW_up", 40.0), ("Frad_SW_down", 200.0), ("Frad_LW_down", 320.0)):
        state[name] = jnp.full_like(state[name], value)
    before = {name: state[name] for name in SOIL_FIELDS}
    _, *expected = advance_soil_veg(
        state["ngrdcol"], state["dt_main"], state["rho_zm"][:, 0],
        state["Frad_SW_up"][:, 0], state["Frad_SW_down"][:, 0], state["Frad_LW_down"][:, 0],
        state["wpthlp_sfc"], state["wprtp_sfc"], state["p_sfc"], state["_jax_stats"],
        *(state[name] for name in SOIL_FIELDS),
    )
    _advance_radiation(state, state["time_initial"], l_rad_itime)
    for name, value in zip(SOIL_FIELDS, expected):
        np.testing.assert_allclose(state[name], value, rtol=1e-7, atol=1e-12)
        assert np.any(np.asarray(state[name]) != np.asarray(before[name]))
