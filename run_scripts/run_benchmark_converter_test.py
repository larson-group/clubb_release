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


RUN_SCRIPTS = Path(__file__).resolve().parent
CLUBB_ROOT = RUN_SCRIPTS.parent
if str(CLUBB_ROOT) not in sys.path:
    sys.path.insert(0, str(CLUBB_ROOT))

from utilities.benchmark_converter import convert_benchmark_file  # noqa: E402


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
        _write_var(ds, "QT2", qt2)
        _write_var(ds, "RTP2_SGS", rtp2_sgs)
        _write_var(ds, "TL2", tl2)
        _write_var(ds, "THLP2_SGS", thlp2_sgs)
        _write_var(ds, "WPRTP", 4.0 + base)
        _write_var(ds, "WPTHLP", 5.0 + base)
        _write_var(ds, "WP3", wp3)

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
        "rcm": np.mean(qc, axis=(2, 3)) * 0.001,
        "rtp2": np.mean(qt2, axis=(2, 3)) / 1.0e6 + np.mean(rtp2_sgs, axis=(2, 3)),
        "thlp2": np.mean(tl2, axis=(2, 3)) + np.mean(thlp2_sgs, axis=(2, 3)),
    }


def assert_close(label: str, actual: np.ndarray, expected: np.ndarray) -> None:
    actual = np.squeeze(np.asarray(actual, dtype=float))
    if not np.allclose(actual, expected, rtol=0.0, atol=1.0e-12):
        raise AssertionError(f"{label} mismatch\nactual={actual}\nexpected={expected}")


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

        for field in ("wp2", "wp2_zt", "W2", "Skw_zt", "em", "rcm", "rtp2_zt", "thlp2_zt", "cloud_frac", "wp3"):
            if status.get(field) != "written":
                raise AssertionError(f"{field} was not written: {status.get(field)}")

        with netCDF4.Dataset(out) as ds:
            if ds.variables["wp2"].dimensions != ("time", "z", "y", "x"):
                raise AssertionError("normalized variables should use time,z,y,x dimensions")
            if getattr(ds.variables["cloud_frac"], "units", None) != "1":
                raise AssertionError("normalized variables should always include a units attribute")
            assert_close("time", ds.variables["time"][:], np.array([0.0, 1.0]))
            assert_close("z", ds.variables["z"][:], np.array([10.0, 30.0, 60.0]))
            assert_close("wp2", ds.variables["wp2"][:], expected["wp2"])
            assert_close("wp2_zt", ds.variables["wp2_zt"][:], expected["wp2"])
            assert_close("W2", ds.variables["W2"][:], expected["W2"])
            assert_close("Skw_zt", ds.variables["Skw_zt"][:], expected["Skw_zt"])
            assert_close("em", ds.variables["em"][:], expected["em"])
            assert_close("rcm", ds.variables["rcm"][:], expected["rcm"])
            assert_close("rtp2_zt", ds.variables["rtp2_zt"][:], expected["rtp2"])
            assert_close("thlp2_zt", ds.variables["thlp2_zt"][:], expected["thlp2"])

            if ds.variables["wp2"].benchmark_source_variables != "WP2":
                raise AssertionError("wp2 should prefer SAM WP2 over W2")
            if ds.variables["wp2_zt"].benchmark_source_variables != "WP2":
                raise AssertionError("wp2_zt should use the same SAM WP2 source")
            if "QT2" not in ds.variables["rtp2_zt"].benchmark_source_variables:
                raise AssertionError("rtp2_zt should record QT2 as a source variable")

        print(f"PASS benchmark converter smoke test: {out}")


if __name__ == "__main__":
    run_test()
