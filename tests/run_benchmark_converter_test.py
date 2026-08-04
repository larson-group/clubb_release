#!/usr/bin/env python3
"""Smoke test for the standalone benchmark normalizer.

This test builds a tiny SAM-like NetCDF with time/z/y/x variables, runs the
converter, and checks that aliases and simple formulas are written correctly.
It does not require external benchmark files.
"""

from __future__ import annotations

from pathlib import Path
import sys
import tempfile

import netCDF4
import numpy as np


TESTS_DIR = Path(__file__).resolve().parent
CLUBB_ROOT = TESTS_DIR.parent
if str(CLUBB_ROOT) not in sys.path:
    sys.path.insert(0, str(CLUBB_ROOT))

from utilities.benchmark_converter import convert_benchmark_file  # noqa: E402
from utilities.les_chi_moments import derive_chi_moments  # noqa: E402


def _write_coord(ds: netCDF4.Dataset, name: str, values: np.ndarray, units: str) -> None:
    var = ds.createVariable(name, "f8", (name,))
    var[:] = values
    var.units = units


def _write_var(ds: netCDF4.Dataset, name: str, data: np.ndarray, units: str = "") -> None:
    var = ds.createVariable(name, "f8", ("time", "z", "y", "x"))
    var[:] = data
    if units:
        var.units = units


def build_sam_fixture(path: Path) -> dict[str, np.ndarray]:
    time = np.array([0.0, 1.0], dtype=float)
    z = np.array([10.0, 30.0, 60.0], dtype=float)
    y = np.arange(2, dtype=float)
    x = np.arange(2, dtype=float)
    shape = (len(time), len(z), len(y), len(x))

    base = np.arange(np.prod(shape), dtype=float).reshape(shape)
    u2 = 1.0 + base
    v2 = 2.0 + base
    wp2 = 10.0 + base
    w2 = 1000.0 + base
    qc = 2.0 + base
    qt2 = 3.0e6 + base
    rtp2_sgs = np.full(shape, 0.25)
    tl2 = 7.0 + base
    thlp2_sgs = np.full(shape, 0.5)
    up2_sgs = np.full(shape, 0.1)
    vp2_sgs = np.full(shape, 0.2)
    wp2_sgs = np.full(shape, 0.3)
    wp3 = 6.0 + base
    rho = 1.2 - 0.001 * base
    rcm = 1.0e-4 + 1.0e-7 * base
    rcp2 = 2.0e-8 + 1.0e-10 * base
    wprcp = 1.0e-5 + 1.0e-7 * base
    wp2rcp = 2.0e-5 + 1.0e-7 * base
    rtprcp = 3.0e-8 + 1.0e-10 * base
    thlprcp = -4.0e-5 - 1.0e-7 * base
    wpthvp = 0.01 + 0.001 * base
    wp2_bp = 1.0e-4 + 1.0e-6 * base
    upwp_bp = 2.0e-4 + 1.0e-6 * base
    vpwp_bp = 3.0e-4 + 1.0e-6 * base
    wprtp_bp = 4.0e-6 + 1.0e-8 * base
    wpthlp_bp = 5.0e-4 + 1.0e-6 * base
    shear = 6.0e-4 + 1.0e-6 * base
    pressure_mb = 950.0 - 0.1 * base
    w_in_cloud = 0.25 + 0.001 * base
    rtm_in_cloud_gkg = 12.0 + 0.01 * base
    rcm_in_cloud_gkg = 0.3 + 0.001 * base
    rtp3 = 1.0e-9 + 1.0e-12 * base
    thlp3 = -0.2 - 0.001 * base
    w_in_cloud[0, 0, 0, 0] = -9999.0
    rtm_in_cloud_gkg[0, 0, 0, 0] = -9999.0
    rcm_in_cloud_gkg[0, 0, 0, 0] = -9999.0

    with netCDF4.Dataset(path, "w") as ds:
        ds.createDimension("time", len(time))
        ds.createDimension("z", len(z))
        ds.createDimension("y", len(y))
        ds.createDimension("x", len(x))
        _write_coord(ds, "time", time, "hours since 2000-01-01 00:00:00")
        _write_coord(ds, "z", z, "m")
        _write_coord(ds, "y", y, "index")
        _write_coord(ds, "x", x, "index")
        _write_var(ds, "U2", u2)
        _write_var(ds, "V2", v2)
        _write_var(ds, "WP2", wp2, "m2/s2")
        _write_var(ds, "W2", w2, "m2/s2")
        _write_var(ds, "UP2_SGS", up2_sgs)
        _write_var(ds, "VP2_SGS", vp2_sgs)
        _write_var(ds, "WP2_SGS", wp2_sgs)
        _write_var(ds, "CLD", 0.1 + 0.01 * base)
        _write_var(ds, "QC", qc, "g/kg")
        _write_var(ds, "RCM", rcm, "kg/kg")
        _write_var(ds, "RCP2", rcp2, "kg2/kg2")
        _write_var(ds, "QT2", qt2)
        _write_var(ds, "RTP2_SGS", rtp2_sgs)
        _write_var(ds, "TL2", tl2)
        _write_var(ds, "THLP2_SGS", thlp2_sgs)
        _write_var(ds, "WPRTP", 4.0 + base)
        _write_var(ds, "WPTHLP", 5.0 + base)
        _write_var(ds, "WP3", wp3)
        _write_var(ds, "RHO", rho, "kg/m3")
        _write_var(ds, "WPRCP", wprcp, "m/s kg/kg")
        _write_var(ds, "WP2RCP", wp2rcp, "m2/s2 kg/kg")
        _write_var(ds, "RTPRCP", rtprcp, "kg2/kg2")
        _write_var(ds, "THLPRCP", thlprcp, "K kg/kg")
        _write_var(ds, "WPTHVP", wpthvp, "m/s K")
        _write_var(ds, "TVFLUX", np.full(shape, 999.0), "W/m2")
        _write_var(ds, "W2BUOY", wp2_bp, "m2/s3")
        _write_var(ds, "WUBUOY", upwp_bp, "m2/s3")
        _write_var(ds, "WVBUOY", vpwp_bp, "m2/s3")
        _write_var(ds, "QWBUOY", wprtp_bp, "m/s2")
        _write_var(ds, "THLWBUOY", wpthlp_bp, "m K/s2")
        _write_var(ds, "SHEAR", shear, "m2/s3")
        _write_var(ds, "PRES", pressure_mb, "mb")
        _write_var(ds, "WCLD", w_in_cloud, "m/s")
        _write_var(ds, "QTCLD", rtm_in_cloud_gkg, "g/kg")
        _write_var(ds, "QCCLD", rcm_in_cloud_gkg, "g/kg")
        _write_var(ds, "QCWCLD", wprcp * 1000.0, "g/kg m/s")
        _write_var(ds, "RTP3", rtp3, "kg3/kg3")
        _write_var(ds, "THLP3", thlp3, "K3")

    u2_total = np.mean(u2 + up2_sgs, axis=(2, 3))
    v2_total = np.mean(v2 + vp2_sgs, axis=(2, 3))
    w2_total = np.mean(w2 + wp2_sgs, axis=(2, 3))
    wp2_profile = np.mean(wp2, axis=(2, 3))
    wp3_profile = np.mean(wp3, axis=(2, 3))
    return {
        "em": 0.5 * (u2_total + v2_total + w2_total),
        "W2": w2_total,
        "Skw_zt": wp3_profile / (wp2_profile + 1.6e-3) ** 1.5,
        "wp2": wp2_profile,
        "rcm": np.mean(rcm, axis=(2, 3)),
        "rcp2": np.mean(rcp2, axis=(2, 3)),
        "rho": np.mean(rho, axis=(2, 3)),
        "wprcp": np.mean(wprcp, axis=(2, 3)),
        "wp2rcp": np.mean(wp2rcp, axis=(2, 3)),
        "rtprcp": np.mean(rtprcp, axis=(2, 3)),
        "thlprcp": np.mean(thlprcp, axis=(2, 3)),
        "wpthvp": np.mean(wpthvp, axis=(2, 3)),
        "wp2_bp": np.mean(wp2_bp, axis=(2, 3)),
        "upwp_bp": np.mean(upwp_bp, axis=(2, 3)),
        "vpwp_bp": np.mean(vpwp_bp, axis=(2, 3)),
        "wprtp_bp": np.mean(wprtp_bp, axis=(2, 3)),
        "wpthlp_bp": np.mean(wpthlp_bp, axis=(2, 3)),
        "shear": np.mean(shear, axis=(2, 3)),
        "p_in_Pa": np.mean(pressure_mb, axis=(2, 3)) * 100.0,
        "rtp2": np.mean(qt2, axis=(2, 3)) / 1.0e6 + np.mean(rtp2_sgs, axis=(2, 3)),
        "thlp2": np.mean(tl2, axis=(2, 3)) + np.mean(thlp2_sgs, axis=(2, 3)),
        "w_in_cloud": np.nanmean(np.where(w_in_cloud <= -9000.0, np.nan, w_in_cloud), axis=(2, 3)),
        "rtm_in_cloud": np.nanmean(
            np.where(rtm_in_cloud_gkg <= -9000.0, np.nan, rtm_in_cloud_gkg), axis=(2, 3)
        ) * 0.001,
        "rcm_in_cloud": np.nanmean(
            np.where(rcm_in_cloud_gkg <= -9000.0, np.nan, rcm_in_cloud_gkg), axis=(2, 3)
        ) * 0.001,
        "rtp3": np.mean(rtp3, axis=(2, 3)),
        "thlp3": np.mean(thlp3, axis=(2, 3)),
        "wrc_cloud": np.mean(wprcp, axis=(2, 3)),
    }


def assert_close(label: str, actual: np.ndarray, expected: np.ndarray) -> None:
    actual = np.squeeze(np.asarray(actual, dtype=float))
    if not np.allclose(actual, expected, rtol=0.0, atol=1.0e-12):
        raise AssertionError(f"{label} mismatch\nactual={actual}\nexpected={expected}")


def check_wprcp_flux_fallback(tmp_path: Path) -> None:
    src = tmp_path / "sam_wprcp_fallback.nc"
    out = tmp_path / "sam_wprcp_fallback_normalized.nc"
    qcflux = np.full((2, 2, 1, 1), 25.104)
    rho = np.full((2, 2, 1, 1), 1.0)
    with netCDF4.Dataset(src, "w") as ds:
        ds.createDimension("time", 2)
        ds.createDimension("z", 2)
        ds.createDimension("y", 1)
        ds.createDimension("x", 1)
        _write_coord(ds, "time", np.array([0.0, 1.0]), "seconds")
        _write_coord(ds, "z", np.array([10.0, 20.0]), "m")
        _write_var(ds, "QCFLUX", qcflux, "W/m2")
        _write_var(ds, "RHO", rho, "kg/m3")

    status = convert_benchmark_file(src, out, source_type="sam", fields=["wprcp"])
    if status["wprcp"] != "written":
        raise AssertionError(f"wprcp fallback was not written: {status['wprcp']}")
    with netCDF4.Dataset(out) as ds:
        assert_close("wprcp fallback", ds.variables["wprcp"][:], np.full((2, 2), 1.0e-5))
        if ds.variables["wprcp"].benchmark_formula != "QCFLUX/(RHO*2.5104e6)":
            raise AssertionError("wprcp fallback formula metadata is incorrect")


def check_thlp2_exact_precedence(tmp_path: Path) -> None:
    src = tmp_path / "sam_thlp2_precedence.nc"
    out = tmp_path / "sam_thlp2_precedence_normalized.nc"
    shape = (2, 2, 1, 1)
    exact_thlp2 = np.full(shape, 0.125)
    with netCDF4.Dataset(src, "w") as ds:
        ds.createDimension("time", 2)
        ds.createDimension("z", 2)
        ds.createDimension("y", 1)
        ds.createDimension("x", 1)
        _write_coord(ds, "time", np.array([0.0, 1.0]), "seconds")
        _write_coord(ds, "z", np.array([10.0, 20.0]), "m")
        _write_var(ds, "THLP2", exact_thlp2, "K2")
        _write_var(ds, "TL2", np.full(shape, 7.0), "K2")
        _write_var(ds, "THLP2_SGS", np.full(shape, 0.5), "K2")

    status = convert_benchmark_file(src, out, source_type="sam", fields=["thlp2_zt"])
    if status["thlp2_zt"] != "written":
        raise AssertionError(f"exact thlp2 was not written: {status['thlp2_zt']}")
    with netCDF4.Dataset(out) as ds:
        result = ds.variables["thlp2_zt"]
        assert_close("exact thlp2 precedence", result[:], np.full((2, 2), 0.125))
        if result.benchmark_source_variables != "THLP2":
            raise AssertionError("thlp2_zt should prefer exact THLP2 over TL2-based reconstruction")


def check_chi_moment_fields(tmp_path: Path) -> None:
    src = tmp_path / "sam_chi_moments.nc"
    out = tmp_path / "sam_chi_moments_normalized.nc"
    shape = (2, 2, 1, 1)
    values = {
        "QT": 11.0,
        "THETAL": 298.0,
        "PRES": 900.0,
        "QT2": 2.5,
        "TL2": 0.25,
        "RTPTHLP": -1.0e-4,
        "WPRTP": 2.0e-4,
        "WPTHLP": -3.0e-2,
        "WP2": 0.5,
    }
    with netCDF4.Dataset(src, "w") as ds:
        ds.createDimension("time", 2)
        ds.createDimension("z", 2)
        ds.createDimension("y", 1)
        ds.createDimension("x", 1)
        _write_coord(ds, "time", np.array([0.0, 1.0]), "seconds")
        _write_coord(ds, "z", np.array([10.0, 20.0]), "m")
        for name, value in values.items():
            _write_var(ds, name, np.full(shape, value))

    status = convert_benchmark_file(
        src,
        out,
        source_type="sam",
        fields=["chi", "chip2", "wpchi"],
    )
    if any(status[name] != "written" for name in ("chi", "chip2", "wpchi")):
        raise AssertionError(f"derived chi fields were not written: {status}")

    expected = derive_chi_moments(
        mean_rt=values["QT"] / 1000.0,
        mean_thl=values["THETAL"],
        pressure_pa=values["PRES"] * 100.0,
        var_rt=values["QT2"] / 1.0e6,
        var_thl=values["TL2"],
        covar_rt_thl=values["RTPTHLP"],
        covar_w_rt=values["WPRTP"],
        covar_w_thl=values["WPTHLP"],
        var_w=values["WP2"],
    )
    with netCDF4.Dataset(out) as ds:
        assert_close("chi", ds.variables["chi"][:], np.full((2, 2), expected["mean_chi"]))
        assert_close("chip2", ds.variables["chip2"][:], np.full((2, 2), expected["var_chi"]))
        assert_close("wpchi", ds.variables["wpchi"][:], np.full((2, 2), expected["covar_w_chi"]))
        if ds.variables["wpchi"].benchmark_formula != "crt*wprtp-cthl*wpthlp":
            raise AssertionError("wpchi formula metadata is incorrect")


def run_test() -> None:
    with tempfile.TemporaryDirectory(prefix="clubb_benchmark_converter_") as tmp:
        tmp_path = Path(tmp)
        src = tmp_path / "sam_fixture.nc"
        out = tmp_path / "sam_normalized.nc"
        expected = build_sam_fixture(src)

        status = convert_benchmark_file(
            src,
            out,
            source_type="sam",
        )

        for field in (
            "wp2", "wp2_zt", "W2", "Skw_zt", "em", "rcm", "rcp2", "rho",
            "rtp2_zt", "thlp2_zt", "cloud_frac", "wp3", "wprcp", "wp2rcp",
            "rtprcp", "thlprcp", "wpthvp", "wp2_bp", "upwp_bp", "vpwp_bp",
            "wprtp_bp", "wpthlp_bp", "shear", "p_in_Pa",
            "w_in_cloud", "rtm_in_cloud", "rcm_in_cloud", "rtp3", "thlp3",
            "wrc_cloud",
        ):
            if status.get(field) != "written":
                raise AssertionError(f"{field} was not written: {status.get(field)}")

        with netCDF4.Dataset(out) as ds:
            if ds.data_model != "NETCDF3_64BIT_OFFSET":
                raise AssertionError(
                    f"normalized benchmark should be NETCDF3_64BIT_OFFSET, got {ds.data_model}"
                )
            if ds.variables["wp2"].dimensions != ("time", "z", "y", "x"):
                raise AssertionError("normalized variables should use time,z,y,x dimensions")
            if getattr(ds.variables["cloud_frac"], "units", None) != "1":
                raise AssertionError("normalized variables should always include a units attribute")
            assert_close("time", ds.variables["time"][:], np.array([0.0, 3600.0]))
            assert_close("z", ds.variables["z"][:], np.array([10.0, 30.0, 60.0]))
            assert_close("wp2", ds.variables["wp2"][:], expected["wp2"])
            assert_close("wp2_zt", ds.variables["wp2_zt"][:], expected["wp2"])
            assert_close("W2", ds.variables["W2"][:], expected["W2"])
            assert_close("Skw_zt", ds.variables["Skw_zt"][:], expected["Skw_zt"])
            assert_close("em", ds.variables["em"][:], expected["em"])
            assert_close("rcm", ds.variables["rcm"][:], expected["rcm"])
            assert_close("rcp2", ds.variables["rcp2"][:], expected["rcp2"])
            assert_close("rho", ds.variables["rho"][:], expected["rho"])
            assert_close("wprcp", ds.variables["wprcp"][:], expected["wprcp"])
            assert_close("wp2rcp", ds.variables["wp2rcp"][:], expected["wp2rcp"])
            assert_close("rtprcp", ds.variables["rtprcp"][:], expected["rtprcp"])
            assert_close("thlprcp", ds.variables["thlprcp"][:], expected["thlprcp"])
            assert_close("wpthvp", ds.variables["wpthvp"][:], expected["wpthvp"])
            assert_close("wp2_bp", ds.variables["wp2_bp"][:], expected["wp2_bp"])
            assert_close("upwp_bp", ds.variables["upwp_bp"][:], expected["upwp_bp"])
            assert_close("vpwp_bp", ds.variables["vpwp_bp"][:], expected["vpwp_bp"])
            assert_close("wprtp_bp", ds.variables["wprtp_bp"][:], expected["wprtp_bp"])
            assert_close("wpthlp_bp", ds.variables["wpthlp_bp"][:], expected["wpthlp_bp"])
            assert_close("shear", ds.variables["shear"][:], expected["shear"])
            assert_close("p_in_Pa", ds.variables["p_in_Pa"][:], expected["p_in_Pa"])
            if ds.variables["p_in_Pa"].units != "Pa":
                raise AssertionError("p_in_Pa should be normalized to Pa")
            assert_close("rtp2_zt", ds.variables["rtp2_zt"][:], expected["rtp2"])
            assert_close("thlp2_zt", ds.variables["thlp2_zt"][:], expected["thlp2"])
            assert_close("w_in_cloud", ds.variables["w_in_cloud"][:], expected["w_in_cloud"])
            assert_close("rtm_in_cloud", ds.variables["rtm_in_cloud"][:], expected["rtm_in_cloud"])
            assert_close("rcm_in_cloud", ds.variables["rcm_in_cloud"][:], expected["rcm_in_cloud"])
            assert_close("rtp3", ds.variables["rtp3"][:], expected["rtp3"])
            assert_close("thlp3", ds.variables["thlp3"][:], expected["thlp3"])
            assert_close("wrc_cloud", ds.variables["wrc_cloud"][:], expected["wrc_cloud"])

            if ds.variables["wp2"].benchmark_source_variables != "WP2":
                raise AssertionError("wp2 should prefer SAM WP2 over W2")
            if ds.variables["wp2_zt"].benchmark_source_variables != "WP2":
                raise AssertionError("wp2_zt should use the same SAM WP2 source")
            if "QT2" not in ds.variables["rtp2_zt"].benchmark_source_variables:
                raise AssertionError("rtp2_zt should record QT2 as a source variable")
            if ds.variables["wprcp"].benchmark_source_variables != "WPRCP":
                raise AssertionError("wprcp should prefer SAM's direct WPRCP covariance")
            if ds.variables["wpthvp"].benchmark_source_variables != "WPTHVP":
                raise AssertionError("wpthvp should prefer direct WPTHVP over TVFLUX fallback")

        check_wprcp_flux_fallback(tmp_path)
        check_thlp2_exact_precedence(tmp_path)
        check_chi_moment_fields(tmp_path)

        print(f"PASS benchmark converter smoke test: {out}")


if __name__ == "__main__":
    run_test()
