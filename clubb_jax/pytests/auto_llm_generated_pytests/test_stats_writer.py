"""Focused contract tests for the F2PY-free stats implementation."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import jax
import jax.numpy as jnp
import pytest
from netCDF4 import Dataset

from clubb_jax.src.CLUBB_core.jax_stats import JaxStats
from clubb_jax.src.CLUBB_core.stats_netcdf import StatsWriter, _parse_registry


def _registry(tmp_path: Path) -> Path:
    path = tmp_path / "stats.in"
    path.write_text(
        """&clubb_stats_nl
  entry(1) = "profile | zt | K | Profile value"
  entry(2) = "moment | zm | m2/s2 | Moment value"
  entry(3) = "surface | sfc | K | Surface value"
  entry(4) = "sclrm | zt | kg/kg | Scalar mean"
  entry(5) = "wpedsclrp | zm | m/s | ED scalar flux"
  entry(6) = "hm_i | zt | kg/kg | Hydrometeor component"
  entry(7) = "silhs_variance_category | lh_zt | 1 | SILHS category"
/
""",
        encoding="ascii",
    )
    return path


def _writer(tmp_path: Path, **overrides) -> StatsWriter:
    args = dict(
        registry_path=str(_registry(tmp_path)),
        output_path="",
        nzt=3,
        nzm=4,
        ngrdcol=2,
        zt=np.array([10.0, 20.0, 30.0]),
        zm=np.array([0.0, 10.0, 20.0, 30.0]),
        stats_tsamp=60.0,
        stats_tout=120.0,
        dt_main=60.0,
        day=1,
        month=1,
        year=2000,
        time_initial=0.0,
    )
    args.update(overrides)
    return StatsWriter(**args)


def _jax_stats(writer: StatsWriter) -> JaxStats:
    return JaxStats.from_layout(writer.get_jax_layout(), ncol=writer.ngrdcol).begin_timestep(
        l_sample=writer.l_sample,
        reset_accumulators=writer.l_reset,
    )


def test_registry_expansion_matches_source_families(tmp_path):
    registry = _parse_registry(
        str(_registry(tmp_path)),
        sclr_dim=2,
        edsclr_dim=2,
        hydromet_list=("rr", "N_r"),
    )

    assert list(registry) == [
        "profile",
        "moment",
        "surface",
        "sclr1m",
        "sclr2m",
        "wpedsclr1p",
        "wpedsclr2p",
        "rr_1",
        "rr_2",
        "N__1",
        "N__2",
        *[f"silhs_var_cat_{i}" for i in range(1, 9)],
    ]


@pytest.mark.parametrize(
    "entry, message",
    (
        ('entry(1) = "profile | zt"', "expected name | grid | units | long_name"),
        ('entry(1) = "profile | bad_grid | K | Profile"', "Invalid stats grid"),
        ('entry(1) = "profile | zt | K | Profile', "Malformed stats registry entry"),
        ('entry(1) = " | zt | K | Profile"', "name and grid must be non-empty"),
    ),
)
def test_registry_errors_are_not_silently_skipped(tmp_path, entry, message):
    registry_path = tmp_path / "invalid_stats.in"
    registry_path.write_text(f"&clubb_stats_nl\n  {entry}\n/\n", encoding="ascii")

    with pytest.raises(ValueError, match=message):
        _parse_registry(str(registry_path))


def test_single_quoted_registry_entry_is_supported(tmp_path):
    registry_path = tmp_path / "single_quote_stats.in"
    registry_path.write_text(
        "&clubb_stats_nl\n  entry(1) = 'profile | zt | K | Profile value'\n/\n",
        encoding="ascii",
    )

    assert _parse_registry(str(registry_path)) == {
        "profile": ("zt", "K", "Profile value")
    }


def test_sampling_uses_open_closed_absolute_time_window(tmp_path):
    writer = _writer(tmp_path, stats_tstart=60.0, stats_tend=180.0)

    assert writer.begin_timestep(0) == (False, False)  # end time = tstart
    assert writer.begin_timestep(1) == (True, False)
    assert writer.begin_timestep(2) == (True, True)
    assert writer.begin_timestep(3) == (False, False)  # after tend


def test_rank_and_selector_dispatch_updates_pointwise_counts(tmp_path):
    writer = _writer(tmp_path)
    assert writer.begin_timestep(0) == (True, False)

    stats = _jax_stats(writer)
    stats = stats.update("profile", 7.0, icol=1, level=2)
    stats = stats.update("profile", np.array([1.0, 2.0, 3.0]), icol=0)
    stats = stats.update("surface", np.array([4.0, 5.0]))
    stats = stats.update("moment", np.full((2, 4), 6.0))
    buffers, nsamples, _, _ = stats.to_host_banks()
    profile_bank, profile_slot = stats.name_to_slot["profile"]
    surface_bank, surface_slot = stats.name_to_slot["surface"]
    moment_bank, moment_slot = stats.name_to_slot["moment"]

    np.testing.assert_array_equal(
        buffers[profile_bank][profile_slot],
        [[1.0, 2.0, 3.0], [0.0, 0.0, 7.0]],
    )
    np.testing.assert_array_equal(
        nsamples[profile_bank][profile_slot],
        [[1, 1, 1], [0, 0, 1]],
    )
    np.testing.assert_array_equal(buffers[surface_bank][surface_slot, :, 0], [4.0, 5.0])
    np.testing.assert_array_equal(nsamples[surface_bank][surface_slot, :, 0], [1, 1])
    np.testing.assert_array_equal(buffers[moment_bank][moment_slot], 6.0)
    np.testing.assert_array_equal(nsamples[moment_bank][moment_slot], 1)


def test_budget_sequence_and_l_count_sample(tmp_path):
    writer = _writer(tmp_path)
    writer.begin_timestep(0)

    stats = _jax_stats(writer)
    stats = stats.begin_budget("profile", np.array([10.0, 20.0, 30.0]), icol=0)
    stats = stats.update_budget("profile", np.array([1.0, 2.0, 3.0]), icol=0)
    stats = stats.finalize_budget(
        "profile",
        np.array([14.0, 25.0, 36.0]),
        icol=0,
        l_count_sample=False,
    )

    buffers, nsamples, _, l_in_budget = stats.to_host_banks()
    bank, slot = stats.name_to_slot["profile"]
    var_id = stats.name_to_id["profile"]
    np.testing.assert_array_equal(buffers[bank][slot, 0], [5.0, 7.0, 9.0])
    np.testing.assert_array_equal(nsamples[bank][slot, 0], 0)
    assert not l_in_budget[var_id]
    assert bool(np.asarray(stats.l_budget)[var_id])


def test_jitted_banks_preserve_column_and_level_selectors(tmp_path):
    writer = _writer(tmp_path)
    writer.begin_timestep(0)
    stats = _jax_stats(writer)

    @jax.jit
    def record(current):
        current = current.update("profile", jnp.asarray(9.0), icol=1, level=2)
        current = current.update("surface", jnp.asarray([3.0, 4.0]))
        return current

    result = record(stats)
    buffers, nsamples, _, _ = result.to_host_banks()
    profile_bank, profile_slot = result.name_to_slot["profile"]
    surface_bank, surface_slot = result.name_to_slot["surface"]

    assert buffers[profile_bank][profile_slot, 1, 2] == 9.0
    assert nsamples[profile_bank][profile_slot, 1, 2] == 1
    np.testing.assert_array_equal(buffers[surface_bank][surface_slot, :, 0], [3.0, 4.0])


def test_jax_layout_reuses_cached_slot_maps(tmp_path):
    writer = _writer(tmp_path)
    first = JaxStats.from_layout(writer.get_jax_layout(), ncol=2)
    second = JaxStats.from_layout(writer.get_jax_layout(), ncol=2)

    assert first.name_to_id is second.name_to_id
    assert first.name_to_slot is second.name_to_slot


def test_rank_one_column_selector_is_bounds_checked(tmp_path):
    writer = _writer(tmp_path)
    writer.begin_timestep(0)
    stats = _jax_stats(writer)

    with pytest.raises(ValueError, match="column selector is out of bounds"):
        stats.update("profile", np.ones(3), icol=2)


def test_direct_banks_accumulate_across_sampled_and_unsampled_steps(tmp_path):
    output = tmp_path / "window.nc"
    writer = _writer(
        tmp_path,
        output_path=str(output),
        stats_tsamp=120.0,
        stats_tout=240.0,
    )
    stats = None
    for itime, expected_sample in enumerate((False, True, False, True)):
        l_sample, l_last_sample = writer.begin_timestep(itime)
        assert l_sample is expected_sample
        if stats is None:
            stats = JaxStats.from_layout(writer.get_jax_layout(), ncol=2)
        stats = stats.begin_timestep(
            l_sample=l_sample,
            reset_accumulators=writer.l_reset,
        )
        if l_sample:
            value = 2.0 if itime == 1 else 4.0
            stats = stats.update("profile", jnp.full((2, 3), value))
        if l_last_sample:
            stats = writer.end_timestep((itime + 1) * 60.0, jax_stats=stats)
    writer.finalize()

    with Dataset(output) as dataset:
        np.testing.assert_array_equal(dataset.variables["profile"][0], 3.0)


def test_netcdf_record_uses_end_step_time(tmp_path):
    output = tmp_path / "stats.nc"
    writer = _writer(tmp_path, output_path=str(output), stats_tout=60.0)
    assert writer.begin_timestep(0) == (True, True)
    stats = _jax_stats(writer).update("profile", np.full((2, 3), 12.0))
    writer.end_timestep(60.0, stats)
    writer.finalize()

    with Dataset(output) as dataset:
        np.testing.assert_array_equal(dataset.variables["time"][:], [60.0])
        np.testing.assert_array_equal(dataset.variables["time_bnds"][:], [[0.0, 60.0]])
        np.testing.assert_array_equal(dataset.variables["profile"][0, :, :], 12.0)
        assert dataset.variables["profile"].units == "K"
        assert dataset.variables["time"].bounds == "time_bnds"
        assert dataset.variables["time"].units == "seconds since 2000-01-01 00:00:00.0"
        assert len(dataset.dimensions["sfc"]) == 1
        assert dataset.stats_tsamp == 60.0
        assert dataset.stats_tout == 60.0


def test_oversampling_is_reported_at_the_host_output_boundary(tmp_path):
    writer = _writer(tmp_path, stats_tout=60.0)
    writer.begin_timestep(0)
    stats = _jax_stats(writer)
    stats = stats.update("profile", np.ones((2, 3)))
    stats = stats.update("profile", np.ones((2, 3)))

    with pytest.warns(RuntimeWarning, match="stats oversampling warning for profile"):
        writer.end_timestep(60.0, stats)


def test_config_metadata_data_and_source_grid_contract(tmp_path):
    writer = _writer(tmp_path)

    assert writer.get_stats_config() == (1, 2, 11, 1, 2, 2, 0, 0, 0)
    assert writer.get_stats_var_meta(0) == (
        "profile",
        "zt",
        "K",
        "Profile value",
        1,
        3,
    )
    zt, zm = writer.get_source_grid()
    assert zt.shape == (2, 3)
    assert zm.shape == (2, 4)
    np.testing.assert_array_equal(zt[0], [10.0, 20.0, 30.0])

    writer.finalize()
    assert writer.get_stats_config()[0] == 0


def test_netcdf_metadata_matches_fortran_schema(tmp_path):
    output = tmp_path / "metadata.nc"
    writer = _writer(
        tmp_path,
        output_path=str(output),
        time_initial=78000.0,
        stats_tstart=78000.0,
        clubb_params_vals=np.ones((2, 2)),
        param_names=["one", "parameter_name_longer_than_28"],
    )
    writer.finalize()

    with Dataset(output) as dataset:
        assert dataset.variables["time"].units == "seconds since 2000-01-01 00:00:00.0"
        assert len(dataset.dimensions["param_strlen"]) == 32
        assert len(dataset.dimensions["sfc"]) == 1
        assert dataset.variables["param"].units == "index"


def test_radiation_registry_entries_require_radiation_grids(tmp_path):
    registry = tmp_path / "radiation_stats.in"
    registry.write_text(
        '&clubb_stats_nl\n entry(1) = "heating | rad_zt | K/s | Heating"\n/\n',
        encoding="ascii",
    )
    writer = _writer(tmp_path, registry_path=str(registry))
    assert "heating" not in writer.registry

    with pytest.raises(ValueError, match="supplied together"):
        _writer(tmp_path, registry_path=str(registry), rad_zt=np.arange(3.0))

    writer = _writer(
        tmp_path,
        registry_path=str(registry),
        rad_zt=np.arange(5.0),
        rad_zm=np.arange(6.0),
    )
    assert writer.get_stats_var_meta(0)[1:] == (
        "rad_zt",
        "K/s",
        "Heating",
        3,
        5,
    )


def test_multiple_runtime_batches_fill_one_output_file(tmp_path):
    output = tmp_path / "batched.nc"
    writer = _writer(
        tmp_path,
        output_path=str(output),
        stats_tout=60.0,
        ncol_total=4,
    )

    writer.begin_timestep(0)
    stats = _jax_stats(writer).update("profile", np.array([[1.0] * 3, [2.0] * 3]))
    stats = writer.end_timestep(60.0, stats)
    writer.reset()
    stats = stats.reset()
    writer.start_next_batch()
    writer.begin_timestep(0)
    stats = stats.begin_timestep(l_sample=writer.l_sample, reset_accumulators=writer.l_reset)
    stats = stats.update("profile", np.array([[3.0] * 3, [4.0] * 3]))
    writer.end_timestep(60.0, stats)
    writer.finalize()

    with Dataset(output) as dataset:
        np.testing.assert_array_equal(
            dataset.variables["profile"][0],
            np.array([[1.0, 2.0, 3.0, 4.0]] * 3),
        )


def test_output_remapping_preserves_values_on_identical_grids(tmp_path):
    output = tmp_path / "remapped.nc"
    zt = np.array([5.0, 15.0, 25.0])
    zm = np.array([0.0, 10.0, 20.0, 30.0])
    writer = _writer(
        tmp_path,
        output_path=str(output),
        zt=zt,
        zm=zm,
        output_zt=zt,
        output_zm=zm,
        grid_remap_method=1,
        stats_tout=60.0,
    )
    rho_levels = np.tile(np.linspace(0.0, 40.0, 9), (2, 1))
    rho_vals = 1.2 * np.exp(-rho_levels / 8000.0)
    writer.update_grid(zt, zm, rho_vals, rho_levels, np.full(2, 1.0e5))

    profile = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    moment = np.array([[7.0, 8.0, 9.0, 10.0], [11.0, 12.0, 13.0, 14.0]])
    writer.begin_timestep(0)
    stats = _jax_stats(writer)
    stats = stats.update("profile", profile)
    stats = stats.update("moment", moment)
    writer.end_timestep(60.0, stats)
    writer.finalize()

    with Dataset(output) as dataset:
        np.testing.assert_allclose(dataset.variables["profile"][0], profile.T, atol=1.0e-12)
        np.testing.assert_allclose(dataset.variables["moment"][0], moment.T, atol=1.0e-12)


def test_silhs_sample_fields_share_the_stats_record(tmp_path):
    output = tmp_path / "silhs.nc"
    writer = _writer(tmp_path, output_path=str(output), stats_tout=60.0)
    writer.stats_lh_samples_init(
        num_samples=2,
        nzt=3,
        nl_var_names=("chi",),
        u_var_names=("u1", "mixture", "weight"),
        zt_vals=np.array([10.0, 20.0, 30.0]),
    )
    writer.begin_timestep(0)
    lognormal = np.arange(12.0).reshape(2, 2, 3, 1)
    uniform = np.arange(12.0, 24.0).reshape(2, 2, 3, 1)
    mixture = np.arange(12).reshape(2, 2, 3)
    weights = np.full((2, 2, 3), 0.25)
    writer.stats_lh_samples_write_lognormal(lognormal)
    writer.stats_lh_samples_write_uniform(uniform, mixture, weights)
    writer.end_timestep(60.0, _jax_stats(writer))
    writer.finalize()

    with Dataset(output) as dataset:
        assert dataset.variables["lh_nl_chi"].dimensions == (
            "time",
            "lh_zt",
            "lh_sample_number",
            "col",
        )
        np.testing.assert_array_equal(
            dataset.variables["lh_nl_chi"][0],
            lognormal[:, :, :, 0].transpose(2, 1, 0),
        )
        np.testing.assert_array_equal(
            dataset.variables["lh_u_mixture"][0],
            mixture.transpose(2, 1, 0),
        )
        np.testing.assert_array_equal(dataset.variables["lh_u_weight"][0], 0.25)
