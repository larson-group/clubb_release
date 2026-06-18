"""Normalize SAM/COAMPS benchmark files to CLUBB-facing field names.

This module is intentionally standalone for now.  It creates a small NetCDF
whose variables are named after CLUBB stats fields, while normalizing the LES
time coordinate to CLUBB-style absolute seconds from midnight on the case date.
Consumers can then do their own time averaging, vertical interpolation, and
plotting/loss calculations without knowing SAM/COAMPS field names.
"""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Callable, Iterable

import numpy as np
from netCDF4 import Dataset


TIME_DIM_NAMES = {"t", "time"}
Z_DIM_NAMES = {"z", "zm", "zt", "altitude", "height", "lev"}
SOURCE_TYPES = {"sam", "coamps"}
CONVERTER_VERSION = "0.3"
_TIME_UNITS_RE = re.compile(r"^\s*([A-Za-z]+)\s+since\s+(.+?)\s*$")


class MissingFieldError(KeyError):
    """Raised when a field cannot be resolved from the input benchmark file."""


def _time_units_factor(unit_name: str) -> float | None:
    unit = str(unit_name or "").strip().lower()
    if unit in {"s", "sec", "secs", "second", "seconds"}:
        return 1.0
    if unit in {"min", "mins", "minute", "minutes"}:
        return 60.0
    if unit in {"h", "hr", "hrs", "hour", "hours"}:
        return 3600.0
    return None


def _parse_time_origin(origin_text: str) -> datetime | None:
    text = str(origin_text or "").strip().replace("T", " ")
    if text.endswith("Z"):
        text = text[:-1].strip()
    for fmt in (
        "%Y-%m-%d %H:%M:%S.%f",
        "%Y-%m-%d %H:%M:%S",
        "%Y-%m-%d %H:%M",
        "%Y-%m-%d",
    ):
        try:
            return datetime.strptime(text, fmt)
        except ValueError:
            continue
    return None


def _normalize_time_to_midnight_seconds(values: np.ndarray, units: str) -> dict[str, object]:
    arr = np.asarray(values, dtype=float)
    raw_units = str(units or "")
    match = _TIME_UNITS_RE.match(raw_units)
    if match:
        factor = _time_units_factor(match.group(1))
        origin = _parse_time_origin(match.group(2))
        if factor is not None and origin is not None:
            midnight = origin.replace(hour=0, minute=0, second=0, microsecond=0)
            origin_offset = (origin - midnight).total_seconds()
            return {
                "values": arr * factor + origin_offset,
                "units": f"seconds since {midnight:%Y-%m-%d %H:%M:%S}.0",
                "origin_offset_seconds": float(origin_offset),
            }

    factor = _time_units_factor(raw_units)
    if factor is not None:
        return {
            "values": arr * factor,
            "units": "seconds",
            "origin_offset_seconds": None,
        }
    return {
        "values": arr,
        "units": raw_units,
        "origin_offset_seconds": None,
    }


@dataclass(frozen=True)
class FieldData:
    """One normalized time-height benchmark field."""

    data: np.ndarray
    time: np.ndarray
    z: np.ndarray
    units: str = ""
    long_name: str = ""
    source_variables: tuple[str, ...] = ()
    formula: str = "raw"
    source_detail: str = ""


class BenchmarkContext:
    """Thin NetCDF reader with helpers for raw profile fields and formulas."""

    def __init__(self, dataset: Dataset, *, source_type: str, input_path: str):
        source_type = str(source_type).lower()
        if source_type not in SOURCE_TYPES:
            raise ValueError(f"Unsupported source type: {source_type}")
        self.dataset = dataset
        self.source_type = source_type
        self.input_path = str(input_path)
        self.time_dim = self._find_dim(TIME_DIM_NAMES)
        self.z_dim = self._find_dim(Z_DIM_NAMES)
        if self.time_dim is None:
            raise ValueError("Could not find a time dimension in benchmark file")
        if self.z_dim is None:
            raise ValueError("Could not find a vertical dimension in benchmark file")
        self.native_time = self._read_coord(self.time_dim)
        self.z = self._read_coord(self.z_dim)
        self.native_time_units = self._coord_units(self.time_dim)
        normalized_time = _normalize_time_to_midnight_seconds(self.native_time, self.native_time_units)
        self.time = normalized_time["values"]
        self.time_units = normalized_time["units"]
        self.time_origin_offset_seconds = normalized_time["origin_offset_seconds"]
        self.z_units = self._coord_units(self.z_dim)

    def _find_dim(self, candidates: set[str]) -> str | None:
        lowered = {name.lower(): name for name in self.dataset.dimensions}
        for candidate in candidates:
            if candidate in lowered:
                return lowered[candidate]
        return None

    def _read_coord(self, dim_name: str) -> np.ndarray:
        if dim_name in self.dataset.variables:
            return np.asarray(self.dataset.variables[dim_name][:], dtype=float)
        return np.arange(len(self.dataset.dimensions[dim_name]), dtype=float)

    def _coord_units(self, dim_name: str) -> str:
        if dim_name in self.dataset.variables:
            return str(getattr(self.dataset.variables[dim_name], "units", "") or "")
        return ""

    def has_raw(self, name: str) -> bool:
        return name in self.dataset.variables

    def read_raw(self, name: str, *, scale: float = 1.0, units: str | None = None) -> FieldData:
        """Read one raw LES variable as a time-height profile.

        Non-time, non-vertical dimensions are horizontally averaged.  This lets
        the converter accept both profile stats and SAM-style time/z/y/x fields.
        """

        if name not in self.dataset.variables:
            raise MissingFieldError(name)

        var = self.dataset.variables[name]
        dims = tuple(var.dimensions)
        if self.time_dim not in dims or self.z_dim not in dims:
            raise MissingFieldError(f"{name} does not have time and z dimensions")

        arr = np.ma.filled(var[:], np.nan).astype(float, copy=False)
        time_axis = dims.index(self.time_dim)
        z_axis = dims.index(self.z_dim)
        arr = np.moveaxis(arr, (time_axis, z_axis), (0, 1))
        if arr.ndim > 2:
            arr = np.nanmean(arr, axis=tuple(range(2, arr.ndim)))
        arr = np.asarray(arr, dtype=float) * float(scale)

        return FieldData(
            data=arr,
            time=np.asarray(self.time, dtype=float),
            z=np.asarray(self.z, dtype=float),
            units=(str(getattr(var, "units", "") or "") if units is None else units),
            long_name=str(getattr(var, "long_name", "") or name),
            source_variables=(name,),
            formula="raw",
            source_detail=name,
        )

    def first_raw(self, candidates: Iterable[str | tuple[str, float] | tuple[str, float, str]]) -> FieldData:
        errors = []
        for item in candidates:
            if isinstance(item, tuple):
                name = item[0]
                scale = item[1]
                units = item[2] if len(item) > 2 else None
            else:
                name, scale, units = item, 1.0, None
            try:
                return self.read_raw(name, scale=scale, units=units)
            except MissingFieldError as exc:
                errors.append(str(exc))
        raise MissingFieldError(", ".join(errors))

    def optional_raw_like(self, name: str, base: FieldData, *, scale: float = 1.0) -> np.ndarray:
        try:
            optional = self.read_raw(name, scale=scale)
        except MissingFieldError:
            return np.zeros_like(base.data)
        _require_same_grid(base, optional)
        return optional.data


Resolver = Callable[[BenchmarkContext], FieldData]


def _require_same_grid(*fields: FieldData) -> None:
    first = fields[0]
    for field in fields[1:]:
        if first.data.shape != field.data.shape:
            raise MissingFieldError("field shapes differ")
        if not np.allclose(first.time, field.time):
            raise MissingFieldError("time coordinates differ")
        if not np.allclose(first.z, field.z):
            raise MissingFieldError("z coordinates differ")


def _derived(
    base: FieldData,
    data: np.ndarray,
    *,
    formula: str,
    source_variables: Iterable[str],
    units: str | None = None,
) -> FieldData:
    return FieldData(
        data=np.asarray(data, dtype=float),
        time=base.time,
        z=base.z,
        units=base.units if units is None else units,
        long_name=base.long_name,
        source_variables=tuple(source_variables),
        formula=formula,
        source_detail=formula,
    )


def _source_vars(*fields: FieldData, extra: Iterable[str] = ()) -> tuple[str, ...]:
    names = []
    for field in fields:
        names.extend(field.source_variables)
    names.extend(extra)
    return tuple(dict.fromkeys(names))


def _sam_add_optional(
    ctx: BenchmarkContext,
    base_name: str,
    optional_name: str,
    *,
    scale: float = 1.0,
    units: str | None = None,
) -> FieldData:
    base = ctx.read_raw(base_name, scale=scale)
    optional = ctx.optional_raw_like(optional_name, base)
    used = list(base.source_variables)
    if ctx.has_raw(optional_name):
        used.append(optional_name)
    return _derived(
        base,
        base.data + optional,
        formula=f"{base_name}*{scale:g}+optional({optional_name})",
        source_variables=used,
        units=units,
    )


def _sam_wpthlp_flux(ctx: BenchmarkContext) -> FieldData:
    tlflux = ctx.read_raw("TLFLUX")
    rho = ctx.read_raw("RHO")
    _require_same_grid(tlflux, rho)
    data = tlflux.data / (rho.data * 1004.0)
    if ctx.has_raw("WPTHLP_SGS"):
        data = data + ctx.optional_raw_like("WPTHLP_SGS", tlflux)
        extra = ("WPTHLP_SGS",)
    else:
        extra = ()
    return _derived(
        tlflux,
        data,
        formula="TLFLUX/(RHO*1004)+optional(WPTHLP_SGS)",
        source_variables=_source_vars(tlflux, rho, extra=extra),
    )


def _sam_wprtp_flux(ctx: BenchmarkContext) -> FieldData:
    qtflux = ctx.read_raw("QTFLUX")
    rho = ctx.read_raw("RHO")
    _require_same_grid(qtflux, rho)
    data = qtflux.data / (rho.data * 2.5104e6)
    if ctx.has_raw("WPRTP_SGS"):
        data = data + ctx.optional_raw_like("WPRTP_SGS", qtflux)
        extra = ("WPRTP_SGS",)
    else:
        extra = ()
    return _derived(
        qtflux,
        data,
        formula="QTFLUX/(RHO*2.5104e6)+optional(WPRTP_SGS)",
        source_variables=_source_vars(qtflux, rho, extra=extra),
    )


def _sam_rtm(ctx: BenchmarkContext) -> FieldData:
    if ctx.has_raw("RTM"):
        return ctx.read_raw("RTM")
    qt = ctx.read_raw("QT")
    if ctx.has_raw("QI"):
        qi = ctx.read_raw("QI")
        _require_same_grid(qt, qi)
        data = (qt.data - qi.data) / 1000.0
        data = np.where(np.isnan(qi.data), qt.data / 1000.0, data)
        return _derived(
            qt,
            data,
            formula="(QT-QI)/1000",
            source_variables=_source_vars(qt, qi),
            units="kg/kg",
        )
    return _derived(
        qt,
        qt.data / 1000.0,
        formula="QT/1000",
        source_variables=qt.source_variables,
        units="kg/kg",
    )


def _sam_rtp2(ctx: BenchmarkContext) -> FieldData:
    if ctx.has_raw("QT2"):
        qt2 = ctx.read_raw("QT2", scale=1.0e-6)
        data = qt2.data
        used = list(qt2.source_variables)
        if ctx.has_raw("RTP2_SGS"):
            data = data + ctx.optional_raw_like("RTP2_SGS", qt2)
            used.append("RTP2_SGS")
        return _derived(
            qt2,
            data,
            formula="QT2/1e6+optional(RTP2_SGS)",
            source_variables=used,
            units="kg2/kg2",
        )
    return ctx.read_raw("RTP2")


def _sam_thlp2(ctx: BenchmarkContext) -> FieldData:
    if ctx.has_raw("TL2"):
        return _sam_add_optional(ctx, "TL2", "THLP2_SGS")
    return ctx.read_raw("THLP2")


def _sam_wp2_equivalent(ctx: BenchmarkContext) -> FieldData:
    # Prefer SAM's central-moment WP2 when present; fall back to legacy W2.
    return ctx.first_raw(("WP2", "W2"))


def _sam_wp3_equivalent(ctx: BenchmarkContext) -> FieldData:
    return ctx.first_raw(("WP3", "W3", "wp3"))


def _sam_wprtp(ctx: BenchmarkContext) -> FieldData:
    if ctx.has_raw("WPRTP"):
        return ctx.read_raw("WPRTP")
    return _sam_wprtp_flux(ctx)


def _sam_wpthlp(ctx: BenchmarkContext) -> FieldData:
    if ctx.has_raw("WPTHLP"):
        return ctx.read_raw("WPTHLP")
    return _sam_wpthlp_flux(ctx)


def _sam_thlm(ctx: BenchmarkContext) -> FieldData:
    if ctx.has_raw("THETAL"):
        thetal = ctx.read_raw("THETAL")
        if ctx.has_raw("THETA") and ctx.has_raw("TABS") and ctx.has_raw("QI"):
            theta = ctx.read_raw("THETA")
            tabs = ctx.read_raw("TABS")
            qi = ctx.read_raw("QI")
            _require_same_grid(thetal, theta, tabs, qi)
            data = thetal.data + 2500.4 * (theta.data / tabs.data) * (qi.data / 1000.0)
            return _derived(
                thetal,
                data,
                formula="THETAL+2500.4*(THETA/TABS)*(QI/1000)",
                source_variables=_source_vars(thetal, theta, tabs, qi),
                units="K",
            )
        return thetal
    return ctx.read_raw("THETA")


def _sam_bv_freq_sqd(ctx: BenchmarkContext) -> FieldData:
    thetav = ctx.read_raw("THETAV")
    if thetav.z.size < 2:
        raise MissingFieldError("THETAV needs at least two z levels")
    grad = np.gradient(thetav.data, thetav.z, axis=1)
    data = 9.81 / thetav.data * grad
    return _derived(
        thetav,
        data,
        formula="9.81/THETAV*d(THETAV)/dz",
        source_variables=thetav.source_variables,
        units="1/s2",
    )


def _raw_by_source(
    ctx: BenchmarkContext,
    *,
    sam: Iterable[str | tuple[str, float] | tuple[str, float, str]] = (),
    coamps: Iterable[str | tuple[str, float] | tuple[str, float, str]] = (),
) -> FieldData:
    if ctx.source_type == "sam":
        return ctx.first_raw(sam)
    return ctx.first_raw(coamps)


def _sam_number_over_rho(ctx: BenchmarkContext, candidates: Iterable[str], *, formula_name: str) -> FieldData:
    number = ctx.first_raw(candidates)
    rho = ctx.read_raw("RHO")
    _require_same_grid(number, rho)
    return _derived(
        number,
        number.data * 1.0e6 / rho.data,
        formula=f"{formula_name}*1e6/RHO",
        source_variables=_source_vars(number, rho),
        units="#/kg",
    )


def _sam_nc_in_cloud(ctx: BenchmarkContext) -> FieldData:
    return _sam_number_over_rho(ctx, ("NC",), formula_name="NC")


def _sam_ncm(ctx: BenchmarkContext) -> FieldData:
    rho = ctx.read_raw("RHO")
    if ctx.has_raw("NC") and ctx.has_raw("CLD"):
        nc = ctx.read_raw("NC")
        cld = ctx.read_raw("CLD")
        _require_same_grid(rho, nc, cld)
        data = cld.data * nc.data * 1.0e6 / rho.data
        if not np.all(np.isclose(data, 0.0, equal_nan=True)):
            return _derived(
                rho,
                data,
                formula="CLD*NC*1e6/RHO",
                source_variables=_source_vars(rho, nc, cld),
                units="#/kg",
            )
    gcssnc = ctx.read_raw("GCSSNC")
    _require_same_grid(rho, gcssnc)
    return _derived(
        rho,
        gcssnc.data * 1.0e6 / rho.data,
        formula="GCSSNC*1e6/RHO",
        source_variables=_source_vars(rho, gcssnc),
        units="#/kg",
    )


def _sam_skw(ctx: BenchmarkContext) -> FieldData:
    numerator = ctx.first_raw(("WP3", "W3", "wp3"))
    denominator = ctx.first_raw(("WP2", "W2", "wp2"))
    _require_same_grid(numerator, denominator)
    data = numerator.data / (denominator.data + 1.6e-3) ** 1.5
    return _derived(
        numerator,
        data,
        formula="WP3/(WP2+1.6e-3)^1.5",
        source_variables=_source_vars(numerator, denominator),
    )


def _sam_skrt(ctx: BenchmarkContext) -> FieldData:
    numerator = ctx.first_raw(("RTP3", "qtp3", "rtp3"))
    denominator = ctx.first_raw(("RTP2", "qtp2", "rtp2", "rlp2"))
    _require_same_grid(numerator, denominator)
    data = numerator.data / (denominator.data + 4.0e-16) ** 1.5
    return _derived(
        numerator,
        data,
        formula="RTP3/(RTP2+4e-16)^1.5",
        source_variables=_source_vars(numerator, denominator),
    )


def _sam_skthl(ctx: BenchmarkContext) -> FieldData:
    numerator = ctx.first_raw(("THLP3", "thlp3"))
    denominator = ctx.first_raw(("THLP2", "thlp2"))
    _require_same_grid(numerator, denominator)
    data = numerator.data / (denominator.data + 4.0e-4) ** 1.5
    return _derived(
        numerator,
        data,
        formula="THLP3/(THLP2+4e-4)^1.5",
        source_variables=_source_vars(numerator, denominator),
    )


def _sam_wpthvp(ctx: BenchmarkContext) -> FieldData:
    if ctx.has_raw("TVFLUX") and ctx.has_raw("RHO"):
        tvflux = ctx.read_raw("TVFLUX")
        rho = ctx.read_raw("RHO")
        _require_same_grid(tvflux, rho)
        return _derived(
            tvflux,
            tvflux.data / (rho.data * 1004.0),
            formula="TVFLUX/(RHO*1004)",
            source_variables=_source_vars(tvflux, rho),
        )
    return ctx.read_raw("WPTHVP")


def _coamps_sum_or_raw(ctx: BenchmarkContext, first_name: str, second_name: str, fallback_name: str) -> FieldData:
    if ctx.has_raw(first_name) and ctx.has_raw(second_name):
        first = ctx.read_raw(first_name)
        second = ctx.read_raw(second_name)
        _require_same_grid(first, second)
        return _derived(
            first,
            first.data + second.data,
            formula=f"{first_name}+{second_name}",
            source_variables=_source_vars(first, second),
        )
    return ctx.read_raw(fallback_name)


def get_bv_freq_sqd(ctx: BenchmarkContext) -> FieldData:
    if ctx.source_type == "sam":
        if ctx.has_raw("BV_FREQ_SQD"):
            return ctx.read_raw("BV_FREQ_SQD")
        return _sam_bv_freq_sqd(ctx)
    raise MissingFieldError("bv_freq_sqd has no COAMPS resolver yet")


def get_bv_freq_sqd_smth(ctx: BenchmarkContext) -> FieldData:
    return _raw_by_source(ctx, sam=("BV_FREQ_SQD_SMTH",), coamps=("bv_freq_sqd_smth",))


def get_c1_skw_fnc(ctx: BenchmarkContext) -> FieldData:
    return _raw_by_source(ctx, sam=("C1_SKW_FNC",), coamps=("C1_Skw_fnc",))


def get_c6_term(ctx: BenchmarkContext) -> FieldData:
    return _raw_by_source(ctx, sam=("C6_TERM",), coamps=("C6_term",))


def get_c6rt_skw_fnc(ctx: BenchmarkContext) -> FieldData:
    return _raw_by_source(ctx, sam=("C6RT_SKW_FNC",), coamps=("C6rt_Skw_fnc",))


def get_c6thl_skw_fnc(ctx: BenchmarkContext) -> FieldData:
    return _raw_by_source(ctx, sam=("C6THL_SKW_FNC",), coamps=("C6thl_Skw_fnc",))


def get_c7_skw_fnc(ctx: BenchmarkContext) -> FieldData:
    return _raw_by_source(ctx, sam=("C7_SKW_FNC",), coamps=("C7_Skw_fnc",))


def get_cloud_frac(ctx: BenchmarkContext) -> FieldData:
    if ctx.source_type == "sam":
        return ctx.first_raw(("CLD_FRAC_CLUBB", "CLD"))
    return ctx.read_raw("cf")


def get_em(ctx: BenchmarkContext) -> FieldData:
    if ctx.source_type == "sam":
        u2 = _sam_add_optional(ctx, "U2", "UP2_SGS")
        v2 = _sam_add_optional(ctx, "V2", "VP2_SGS")
        w2 = _sam_add_optional(ctx, "W2", "WP2_SGS")
        _require_same_grid(u2, v2, w2)
        return _derived(
            u2,
            0.5 * (u2.data + v2.data + w2.data),
            formula="0.5*((U2+optional(UP2_SGS))+(V2+optional(VP2_SGS))+(W2+optional(WP2_SGS)))",
            source_variables=_source_vars(u2, v2, w2),
        )
    return ctx.read_raw("em")


def get_ice_supersat_frac(ctx: BenchmarkContext) -> FieldData:
    return ctx.read_raw("ice_supersat_frac")


def get_invrs_tau_no_n2_zm(ctx: BenchmarkContext) -> FieldData:
    return _raw_by_source(ctx, sam=("INVRS_TAU_NO_N2_ZM",), coamps=("invrs_tau_no_N2_zm",))


def get_invrs_tau_wp2_zm(ctx: BenchmarkContext) -> FieldData:
    return _raw_by_source(ctx, sam=("INVRS_TAU_WP2_ZM",), coamps=("invrs_tau_wp2_zm",))


def get_invrs_tau_wp3_zm(ctx: BenchmarkContext) -> FieldData:
    return _raw_by_source(ctx, sam=("INVRS_TAU_WP3_ZM",), coamps=("invrs_tau_wp3_zm",))


def get_invrs_tau_wpxp_zm(ctx: BenchmarkContext) -> FieldData:
    return _raw_by_source(ctx, sam=("INVRS_TAU_WPXP_ZM",), coamps=("invrs_tau_wpxp_zm",))


def get_invrs_tau_xp2_zm(ctx: BenchmarkContext) -> FieldData:
    return _raw_by_source(ctx, sam=("INVRS_TAU_XP2_ZM",), coamps=("invrs_tau_xp2_zm",))


def get_invrs_tau_zm(ctx: BenchmarkContext) -> FieldData:
    return _raw_by_source(ctx, sam=("INVRS_TAU_ZM",), coamps=("invrs_tau_zm",))


def get_kh_zm(ctx: BenchmarkContext) -> FieldData:
    return _raw_by_source(ctx, sam=("KH_ZM",), coamps=("Kh_zm",))


def get_nc_in_cloud(ctx: BenchmarkContext) -> FieldData:
    if ctx.source_type == "sam":
        return _sam_nc_in_cloud(ctx)
    return ctx.read_raw("Nc_in_cloud")


def get_ncm(ctx: BenchmarkContext) -> FieldData:
    if ctx.source_type == "sam":
        return _sam_ncm(ctx)
    return ctx.first_raw(("Ncm", "ncm"))


def get_ngm(ctx: BenchmarkContext) -> FieldData:
    if ctx.source_type == "sam":
        return ctx.read_raw("NG", scale=1.0e6, units="#/kg")
    return ctx.read_raw("Ngm")


def get_nim(ctx: BenchmarkContext) -> FieldData:
    if ctx.source_type == "sam":
        return _sam_number_over_rho(ctx, ("NI",), formula_name="NI")
    return ctx.first_raw(("Nim", "nim"))


def get_nrm(ctx: BenchmarkContext) -> FieldData:
    if ctx.source_type == "sam":
        return _sam_number_over_rho(ctx, ("NR", "CONP"), formula_name="NR_or_CONP")
    return ctx.read_raw("Nrm")


def get_nsm(ctx: BenchmarkContext) -> FieldData:
    if ctx.source_type == "sam":
        return _sam_number_over_rho(ctx, ("NS",), formula_name="NS")
    return ctx.read_raw("Nsm")


def get_precip_frac(ctx: BenchmarkContext) -> FieldData:
    return ctx.read_raw("precip_frac")


def get_radht(ctx: BenchmarkContext) -> FieldData:
    if ctx.source_type == "sam":
        return ctx.read_raw("RADQR", scale=1.15741e-5, units="K/s")
    return ctx.read_raw("radht")


def get_rcm(ctx: BenchmarkContext) -> FieldData:
    if ctx.source_type == "sam":
        return ctx.first_raw((("QCL", 0.001, "kg/kg"), ("QC", 0.001, "kg/kg")))
    return ctx.first_raw(("qcm", "rcm"))


def get_rcp2(ctx: BenchmarkContext) -> FieldData:
    return _raw_by_source(ctx, sam=("QC2",), coamps=("qcp2", "rcp2", "rlp2"))


def get_rgm(ctx: BenchmarkContext) -> FieldData:
    return _raw_by_source(ctx, sam=(("QG", 0.001, "kg/kg"),), coamps=("rgm", "qgm"))


def get_ri_zm(ctx: BenchmarkContext) -> FieldData:
    return _raw_by_source(ctx, sam=("RI_ZM",), coamps=("Ri_zm",))


def get_rim(ctx: BenchmarkContext) -> FieldData:
    return _raw_by_source(ctx, sam=(("QI", 0.001, "kg/kg"), ("QCI", 0.001, "kg/kg")), coamps=("rim", "qim"))


def get_rrm(ctx: BenchmarkContext) -> FieldData:
    return _raw_by_source(ctx, sam=(("QPL", 0.001, "kg/kg"),), coamps=("rrm",))


def get_rrm_accr(ctx: BenchmarkContext) -> FieldData:
    return _raw_by_source(ctx, sam=("ACCRM",), coamps=("rrm_accr",))


def get_rrm_auto(ctx: BenchmarkContext) -> FieldData:
    return _raw_by_source(ctx, sam=("AUTOM",), coamps=("rrm_auto",))


def get_rrm_evap(ctx: BenchmarkContext) -> FieldData:
    return _raw_by_source(ctx, sam=("EVAPM",), coamps=("rrm_evap",))


def get_rsm(ctx: BenchmarkContext) -> FieldData:
    return _raw_by_source(ctx, sam=(("QS", 0.001, "kg/kg"), ("QPI", 0.001, "kg/kg")), coamps=("rsm", "qsm"))


def get_rtm(ctx: BenchmarkContext) -> FieldData:
    if ctx.source_type == "sam":
        return _sam_rtm(ctx)
    return ctx.first_raw(("qtm", "rtm"))


def get_rtp2(ctx: BenchmarkContext) -> FieldData:
    if ctx.source_type == "sam":
        return _sam_rtp2(ctx)
    return ctx.first_raw(("qtp2", "rtp2"))


def get_rtp2_zt(ctx: BenchmarkContext) -> FieldData:
    return get_rtp2(ctx)


def get_rtp3(ctx: BenchmarkContext) -> FieldData:
    return _raw_by_source(ctx, sam=("RTP3",), coamps=("qtp3", "rtp3"))


def get_rtpthlp(ctx: BenchmarkContext) -> FieldData:
    return _raw_by_source(ctx, sam=("RTPTHLP_SGS", "RTPTHLP", "TQ"), coamps=("qtpthlp", "rtpthlp"))


def get_rtpthvp(ctx: BenchmarkContext) -> FieldData:
    return _raw_by_source(ctx, sam=("RTPTHVP",), coamps=("qtpthvp", "rtpthvp"))


def get_sclr1m(ctx: BenchmarkContext) -> FieldData:
    if ctx.source_type == "sam":
        raise MissingFieldError("sclr1m has no SAM resolver")
    return ctx.read_raw("sclr1m")


def get_sclr1p2(ctx: BenchmarkContext) -> FieldData:
    if ctx.source_type == "sam":
        raise MissingFieldError("sclr1p2 has no SAM resolver")
    return ctx.read_raw("sclr1p2")


def get_sclr2m(ctx: BenchmarkContext) -> FieldData:
    if ctx.source_type == "sam":
        raise MissingFieldError("sclr2m has no SAM resolver")
    return ctx.read_raw("sclr2m")


def get_sclr2p2(ctx: BenchmarkContext) -> FieldData:
    if ctx.source_type == "sam":
        raise MissingFieldError("sclr2p2 has no SAM resolver")
    return ctx.read_raw("sclr2p2")


def get_skrt_zt(ctx: BenchmarkContext) -> FieldData:
    return _sam_skrt(ctx)


def get_skthl_zt(ctx: BenchmarkContext) -> FieldData:
    return _sam_skthl(ctx)


def get_skw_zt(ctx: BenchmarkContext) -> FieldData:
    return _sam_skw(ctx)


def get_thlp2(ctx: BenchmarkContext) -> FieldData:
    if ctx.source_type == "sam":
        return _sam_thlp2(ctx)
    return ctx.read_raw("thlp2")


def get_thlp2_zt(ctx: BenchmarkContext) -> FieldData:
    return get_thlp2(ctx)


def get_thlp3(ctx: BenchmarkContext) -> FieldData:
    return _raw_by_source(ctx, sam=("THLP3",), coamps=("thlp3",))


def get_thlpthvp(ctx: BenchmarkContext) -> FieldData:
    return _raw_by_source(ctx, sam=("thlpthvp", "THLPTHVP"), coamps=("thlpthvp",))


def get_wp2(ctx: BenchmarkContext) -> FieldData:
    if ctx.source_type == "sam":
        return _sam_wp2_equivalent(ctx)
    return ctx.first_raw(("wp2", "W2"))


def get_wp2_zt(ctx: BenchmarkContext) -> FieldData:
    return get_wp2(ctx)


def get_w2(ctx: BenchmarkContext) -> FieldData:
    if ctx.source_type == "sam":
        return _sam_add_optional(ctx, "W2", "WP2_SGS")
    return ctx.first_raw(("wp2", "W2"))


def get_w_up_in_cloud(ctx: BenchmarkContext) -> FieldData:
    if ctx.source_type == "sam":
        return ctx.read_raw("WSUP")
    raise MissingFieldError("w_up_in_cloud has no COAMPS resolver")


def get_wlsm(ctx: BenchmarkContext) -> FieldData:
    return _raw_by_source(ctx, sam=("WOBS", "WM"), coamps=("wlsm", "wm"))


def get_wm(ctx: BenchmarkContext) -> FieldData:
    return get_wlsm(ctx)


def get_wp3(ctx: BenchmarkContext) -> FieldData:
    if ctx.source_type == "sam":
        return _sam_wp3_equivalent(ctx)
    return ctx.first_raw(("wp3", "W3", "WP3"))


def get_upwp(ctx: BenchmarkContext) -> FieldData:
    if ctx.source_type == "sam":
        return _sam_add_optional(ctx, "UW", "UPWP_SGS")
    return _coamps_sum_or_raw(ctx, "wpup", "wpup_sgs", "upwp")


def get_vpwp(ctx: BenchmarkContext) -> FieldData:
    if ctx.source_type == "sam":
        return _sam_add_optional(ctx, "VW", "VPWP_SGS")
    return _coamps_sum_or_raw(ctx, "wpvp", "wpvp_sgs", "vpwp")


def get_wp2rtp(ctx: BenchmarkContext) -> FieldData:
    return _raw_by_source(ctx, sam=("WP2RTP",), coamps=("wp2qtp", "wp2rtp"))


def get_wp2thlp(ctx: BenchmarkContext) -> FieldData:
    return _raw_by_source(ctx, sam=("WP2THLP",), coamps=("wp2thlp",))


def get_wp2thvp(ctx: BenchmarkContext) -> FieldData:
    return _raw_by_source(ctx, sam=("WP2THVP",), coamps=("wp2thvp",))


def get_wp2up2(ctx: BenchmarkContext) -> FieldData:
    return _raw_by_source(ctx, sam=("WP2UP2",), coamps=("wp2up2",))


def get_wp2vp2(ctx: BenchmarkContext) -> FieldData:
    return _raw_by_source(ctx, sam=("WP2VP2",), coamps=("wp2vp2",))


def get_wp4(ctx: BenchmarkContext) -> FieldData:
    return _raw_by_source(ctx, sam=("WP4",), coamps=("wp4",))


def get_wpnrp(ctx: BenchmarkContext) -> FieldData:
    return _raw_by_source(ctx, sam=("WPNRP",), coamps=("wpNrp",))


def get_wprrp(ctx: BenchmarkContext) -> FieldData:
    return _raw_by_source(ctx, sam=("WPRRP",), coamps=("wprrp",))


def get_wprtp(ctx: BenchmarkContext) -> FieldData:
    if ctx.source_type == "sam":
        return _sam_wprtp(ctx)
    return ctx.first_raw(("wpqtp", "wprtp"))


def get_wprtp_zt(ctx: BenchmarkContext) -> FieldData:
    return get_wprtp(ctx)


def get_wprtp2(ctx: BenchmarkContext) -> FieldData:
    return _raw_by_source(ctx, sam=("WPRTP2",), coamps=("wpqtp2", "wprtp2"))


def get_wprtpthlp(ctx: BenchmarkContext) -> FieldData:
    return _raw_by_source(ctx, sam=("WPRTPTHLP",), coamps=("wpqtpthlp", "wprtpthlp"))


def get_wpthlp(ctx: BenchmarkContext) -> FieldData:
    if ctx.source_type == "sam":
        return _sam_wpthlp(ctx)
    return ctx.read_raw("wpthlp")


def get_wpthlp_zt(ctx: BenchmarkContext) -> FieldData:
    return get_wpthlp(ctx)


def get_wpthlp2(ctx: BenchmarkContext) -> FieldData:
    return _raw_by_source(ctx, sam=("WPTHLP2",), coamps=("wpthlp2",))


def get_wpthvp(ctx: BenchmarkContext) -> FieldData:
    if ctx.source_type == "sam":
        return _sam_wpthvp(ctx)
    return ctx.read_raw("wpthvp")


def get_wpup2(ctx: BenchmarkContext) -> FieldData:
    if ctx.source_type == "sam":
        raise MissingFieldError("wpup2 has no SAM resolver")
    return ctx.read_raw("wpup2")


def get_wpvp2(ctx: BenchmarkContext) -> FieldData:
    if ctx.source_type == "sam":
        raise MissingFieldError("wpvp2 has no SAM resolver")
    return ctx.read_raw("wpvp2")


def get_up2(ctx: BenchmarkContext) -> FieldData:
    if ctx.source_type == "sam":
        return _sam_add_optional(ctx, "U2", "UP2_SGS")
    return ctx.read_raw("up2")


def get_vp2(ctx: BenchmarkContext) -> FieldData:
    if ctx.source_type == "sam":
        return _sam_add_optional(ctx, "V2", "VP2_SGS")
    return ctx.read_raw("vp2")


def get_thlm(ctx: BenchmarkContext) -> FieldData:
    if ctx.source_type == "sam":
        return _sam_thlm(ctx)
    return ctx.read_raw("thlm")


def get_um(ctx: BenchmarkContext) -> FieldData:
    if ctx.source_type == "sam":
        return ctx.read_raw("U")
    return ctx.read_raw("um")


def get_vm(ctx: BenchmarkContext) -> FieldData:
    if ctx.source_type == "sam":
        return ctx.read_raw("V")
    return ctx.read_raw("vm")


FIELD_RESOLVERS: dict[str, Resolver] = {
    "BV_FREQ_SQD": get_bv_freq_sqd,
    "C1_Skw_fnc": get_c1_skw_fnc,
    "C6_term": get_c6_term,
    "C6rt_Skw_fnc": get_c6rt_skw_fnc,
    "C6thl_Skw_fnc": get_c6thl_skw_fnc,
    "C7_Skw_fnc": get_c7_skw_fnc,
    "Kh_zm": get_kh_zm,
    "Nc_in_cloud": get_nc_in_cloud,
    "Ncm": get_ncm,
    "Ngm": get_ngm,
    "Nim": get_nim,
    "Nrm": get_nrm,
    "Nsm": get_nsm,
    "Ri_zm": get_ri_zm,
    "Skrt_zt": get_skrt_zt,
    "Skthl_zt": get_skthl_zt,
    "Skw_zt": get_skw_zt,
    "W2": get_w2,
    "bv_freq_sqd": get_bv_freq_sqd,
    "bv_freq_sqd_smth": get_bv_freq_sqd_smth,
    "cloud_frac": get_cloud_frac,
    "em": get_em,
    "ice_supersat_frac": get_ice_supersat_frac,
    "invrs_tau_no_N2_zm": get_invrs_tau_no_n2_zm,
    "invrs_tau_wp2_zm": get_invrs_tau_wp2_zm,
    "invrs_tau_wp3_zm": get_invrs_tau_wp3_zm,
    "invrs_tau_wpxp_zm": get_invrs_tau_wpxp_zm,
    "invrs_tau_xp2_zm": get_invrs_tau_xp2_zm,
    "invrs_tau_zm": get_invrs_tau_zm,
    "ncm": get_ncm,
    "precip_frac": get_precip_frac,
    "radht": get_radht,
    "rcm": get_rcm,
    "rcp2": get_rcp2,
    "rgm": get_rgm,
    "rim": get_rim,
    "rrm": get_rrm,
    "rrm_accr": get_rrm_accr,
    "rrm_auto": get_rrm_auto,
    "rrm_evap": get_rrm_evap,
    "rsm": get_rsm,
    "rtm": get_rtm,
    "rtp2": get_rtp2,
    "rtp2_zt": get_rtp2_zt,
    "rtp3": get_rtp3,
    "rtpthlp": get_rtpthlp,
    "rtpthvp": get_rtpthvp,
    "sclr1m": get_sclr1m,
    "sclr1p2": get_sclr1p2,
    "sclr2m": get_sclr2m,
    "sclr2p2": get_sclr2p2,
    "thlp3": get_thlp3,
    "thlpthvp": get_thlpthvp,
    "thlm": get_thlm,
    "thlp2": get_thlp2,
    "thlp2_zt": get_thlp2_zt,
    "um": get_um,
    "up2": get_up2,
    "upwp": get_upwp,
    "vm": get_vm,
    "vp2": get_vp2,
    "vpwp": get_vpwp,
    "w_up_in_cloud": get_w_up_in_cloud,
    "wlsm": get_wlsm,
    "wm": get_wm,
    "wp2": get_wp2,
    "wp2_zt": get_wp2_zt,
    "wp2rtp": get_wp2rtp,
    "wp2thlp": get_wp2thlp,
    "wp2thvp": get_wp2thvp,
    "wp2up2": get_wp2up2,
    "wp2vp2": get_wp2vp2,
    "wp3": get_wp3,
    "wp4": get_wp4,
    "wpNrp": get_wpnrp,
    "wprrp": get_wprrp,
    "wprtp": get_wprtp,
    "wprtp_zt": get_wprtp_zt,
    "wprtp2": get_wprtp2,
    "wprtpthlp": get_wprtpthlp,
    "wpthlp": get_wpthlp,
    "wpthlp_zt": get_wpthlp_zt,
    "wpthlp2": get_wpthlp2,
    "wpthvp": get_wpthvp,
    "wpup2": get_wpup2,
    "wpvp2": get_wpvp2,
}


def supported_fields() -> tuple[str, ...]:
    return tuple(sorted(FIELD_RESOLVERS))


def detect_source_type(path: str | Path) -> str:
    with Dataset(path, "r") as ds:
        names = set(ds.variables)
        sam_score = len(names & {"CLD", "CLD_FRAC_CLUBB", "QC", "QCL", "WP2", "W2", "QT2", "TL2"})
        coamps_score = len(names & {"cf", "qcm", "rcm", "wp2", "qtp2", "thlp2", "up2", "vp2"})
    if sam_score > coamps_score:
        return "sam"
    if coamps_score > sam_score:
        return "coamps"
    raise ValueError("Could not auto-detect LES source type; pass --source sam or --source coamps")


def _create_coord(out_ds: Dataset, name: str, values: np.ndarray, units: str) -> None:
    var = out_ds.createVariable(name, "f8", (name,))
    var[:] = values
    if units:
        var.units = units


def _write_field(out_ds: Dataset, field_name: str, field: FieldData) -> None:
    var = out_ds.createVariable(field_name, "f8", ("time", "z", "y", "x"), zlib=True)
    var[:, :, 0, 0] = field.data
    var.clubb_field_name = field_name
    var.benchmark_formula = field.formula
    var.benchmark_source_variables = ",".join(field.source_variables)
    var.benchmark_source_detail = field.source_detail
    var.units = field.units or "1"
    var.long_name = field.long_name or field_name


def convert_benchmark_file(
    input_path: str | Path,
    output_path: str | Path,
    *,
    source_type: str = "auto",
    fields: Iterable[str] | None = None,
    clobber: bool = True,
) -> dict[str, str]:
    """Create a normalized benchmark NetCDF and return field resolution status."""

    input_path = Path(input_path)
    output_path = Path(output_path)
    if source_type == "auto":
        source_type = detect_source_type(input_path)
    source_type = source_type.lower()
    if source_type not in SOURCE_TYPES:
        raise ValueError(f"Unsupported source type: {source_type}")

    requested_fields = list(fields or supported_fields())
    unknown = sorted(set(requested_fields) - set(FIELD_RESOLVERS))
    if unknown:
        raise ValueError("Unsupported normalized field(s): " + ", ".join(unknown))

    mode = "w" if clobber else "x"
    status: dict[str, str] = {}

    with Dataset(input_path, "r") as in_ds, Dataset(output_path, mode) as out_ds:
        ctx = BenchmarkContext(in_ds, source_type=source_type, input_path=str(input_path))
        out_ds.createDimension("time", len(ctx.time))
        out_ds.createDimension("z", len(ctx.z))
        out_ds.createDimension("y", 1)
        out_ds.createDimension("x", 1)
        _create_coord(out_ds, "time", ctx.time, ctx.time_units)
        _create_coord(out_ds, "z", ctx.z, ctx.z_units)
        _create_coord(out_ds, "y", np.array([1.0]), "index")
        _create_coord(out_ds, "x", np.array([1.0]), "index")

        out_ds.normalized_benchmark_converter = "utilities.benchmark_converter"
        out_ds.normalized_benchmark_converter_version = CONVERTER_VERSION
        out_ds.normalized_benchmark_source_file = str(input_path)
        out_ds.normalized_benchmark_source_type = source_type
        out_ds.normalized_benchmark_source_time_units = ctx.native_time_units
        if ctx.time_origin_offset_seconds is not None:
            out_ds.normalized_benchmark_source_time_origin_offset_seconds = ctx.time_origin_offset_seconds
        out_ds.normalized_benchmark_status_note = "Variables that could not be resolved are omitted."

        for field_name in requested_fields:
            resolver = FIELD_RESOLVERS[field_name]
            try:
                field = resolver(ctx)
                _require_same_grid(
                    FieldData(np.empty((len(ctx.time), len(ctx.z))), ctx.time, ctx.z),
                    field,
                )
            except MissingFieldError as exc:
                status[field_name] = f"missing: {exc}"
                continue
            _write_field(out_ds, field_name, field)
            status[field_name] = "written"

        out_ds.normalized_benchmark_written_fields = ",".join(
            name for name, result in status.items() if result == "written"
        )
        out_ds.normalized_benchmark_missing_fields = ",".join(
            f"{name}={result}" for name, result in status.items() if result != "written"
        )

    return status


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Convert a SAM/COAMPS benchmark file to CLUBB-facing field names.")
    parser.add_argument("input", nargs="?", help="Input SAM or COAMPS benchmark NetCDF file.")
    parser.add_argument("output", nargs="?", help="Output normalized benchmark NetCDF file.")
    parser.add_argument("--source", choices=("auto", "sam", "coamps"), default="auto", help="Input LES source type.")
    parser.add_argument("--fields", nargs="+", help="Optional CLUBB field names to convert. Defaults to all supported fields.")
    parser.add_argument("--no-clobber", action="store_true", help="Fail if the output file already exists.")
    parser.add_argument("--list-fields", action="store_true", help="Print supported normalized field names and exit.")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if args.list_fields:
        for field_name in supported_fields():
            print(field_name)
        return 0
    if not args.input or not args.output:
        raise SystemExit("input and output are required unless --list-fields is used")
    status = convert_benchmark_file(
        args.input,
        args.output,
        source_type=args.source,
        fields=args.fields,
        clobber=not args.no_clobber,
    )
    for field_name in sorted(status):
        print(f"{field_name}: {status[field_name]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
