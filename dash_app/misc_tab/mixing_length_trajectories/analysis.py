"""Python replica of ``compute_mixing_length`` with parcel-path diagnostics.

The calculation below intentionally follows the Fortran routine's control
flow.  Besides the directional and final length profiles, it retains the
remaining parcel energy at every traversed model level for visualization.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import netCDF4 as nc
import numpy as np

from dash_app.shared.output_discovery import discover_stats_files

CP = 1004.67
LV = 2.5e6
RD = 287.04
RV = 461.5
EP = RD / RV
EP1 = (1.0 - EP) / EP
EP2 = 1.0 / EP
GRAV = 9.81
LV2_COEF = EP * LV**2 / (RD * CP)
ZLMIN = 0.1
LSCALE_SURFACE_LAYER_DEPTH = 500.0
LSCALE_MAX_STANDALONE = 1.0e5
LMIN_DELTAZ = 40.0

REQUIRED_VARIABLES = (
    "zt",
    "zm",
    "time",
    "thlm_old",
    "rtm_old",
    "thlm",
    "rtm",
    "thvm",
    "p_in_Pa",
    "exner",
    "thv_ds_zt",
    "em",
    "Lscale",
    "Lscale_up",
    "Lscale_down",
)


@dataclass(frozen=True)
class ParcelTrajectory:
    """Remaining-energy samples for one parcel launch."""

    launch_index: int
    altitude: np.ndarray
    energy: np.ndarray
    buoyancy: np.ndarray
    parcel_thl: np.ndarray
    parcel_rt: np.ndarray
    energy_endpoint: float
    raw_length: float


@dataclass(frozen=True)
class MixingLengthResult:
    """Replica output and the paths used to obtain it."""

    z: np.ndarray
    tke: np.ndarray
    raw_up: np.ndarray
    raw_down: np.ndarray
    enveloped_up: np.ndarray
    enveloped_down: np.ndarray
    lscale_up: np.ndarray
    lscale_down: np.ndarray
    lscale: np.ndarray
    upward_paths: tuple[ParcelTrajectory, ...]
    downward_paths: tuple[ParcelTrajectory, ...]


@dataclass(frozen=True)
class DatasetRecord:
    """One NetCDF record in the state consumed by ``compute_mixing_length``."""

    path: Path
    record: int
    column: int
    time_seconds: float
    z: np.ndarray
    zm: np.ndarray
    thlm: np.ndarray
    rtm: np.ndarray
    mean_thlm: np.ndarray
    mean_rtm: np.ndarray
    thvm: np.ndarray
    pressure: np.ndarray
    exner: np.ndarray
    thv_ds: np.ndarray
    em: np.ndarray
    fortran_lscale: np.ndarray
    fortran_up: np.ndarray
    fortran_down: np.ndarray
    mu: float
    lmin: float


def saturation_mixing_ratio_flatau(pressure, temperature):
    """Return CLUBB's liquid saturation mixing ratio for formula 3."""
    temperature_c = np.maximum(-85.0, np.asarray(temperature) - 273.15)
    temperature_c2 = temperature_c**2
    saturation_pressure = (
        -3.21582393e-14
        * (temperature_c - 646.5835252598777)
        * (temperature_c + 90.72381630364440)
        * (
            temperature_c2
            + 111.0976961559954 * temperature_c
            + 6459.629194243118
        )
        * (
            temperature_c2
            + 152.3131930092453 * temperature_c
            + 6499.774954705265
        )
        * (
            temperature_c2
            + 174.4279584934021 * temperature_c
            + 7721.679732114084
        )
    )
    return EP * saturation_pressure / (
        np.asarray(pressure) - saturation_pressure
    )


def _parcel_buoyancy(
    level,
    thl_par,
    rt_par,
    thvm,
    pressure,
    exner,
    thv_ds,
):
    """Compute the moist parcel buoyancy at one target level."""
    liquid_temperature = thl_par * exner[level]
    rsatl_par = saturation_mixing_ratio_flatau(
        pressure[level], liquid_temperature
    )
    liquid_temperature_sqd = liquid_temperature**2
    saturation_excess = (
        (rt_par - rsatl_par)
        * liquid_temperature_sqd
        / (liquid_temperature_sqd + LV2_COEF * rsatl_par)
    )
    rc_par = max(float(saturation_excess), 0.0)
    lv_coef = LV / (exner[level] * CP) - EP2 * thv_ds[level]
    thv_par = (
        thl_par
        + EP1 * thv_ds[level] * rt_par
        + lv_coef * rc_par
    )
    return GRAV * (thv_par - thvm[level]) / thvm[level]


def _entrain_one_layer(parcel, environment_lower, environment_upper, dz, mu):
    """Apply the exact layer-linear entrainment recurrence used by Fortran."""
    exp_mu_dz = np.exp(-mu * dz)
    entrain_coef = (1.0 - exp_mu_dz) / (mu * dz)
    return (
        environment_upper
        - environment_lower * exp_mu_dz
        - (environment_upper - environment_lower) * entrain_coef
        + parcel * exp_mu_dz
    )


def _energy_root(energy, buoyancy_start, buoyancy_end, dz, work_sign):
    """Locate K=0 in one layer with the Fortran linear-buoyancy assumption."""
    linear = work_sign * buoyancy_start
    quadratic = (
        0.5
        * work_sign
        * (buoyancy_end - buoyancy_start)
        / dz
    )
    if abs(quadratic) <= np.finfo(float).eps * max(abs(linear), 1.0):
        distance = -energy / linear
    else:
        discriminant = max(linear * linear - 4.0 * quadratic * energy, 0.0)
        root_1 = (-linear + np.sqrt(discriminant)) / (2.0 * quadratic)
        root_2 = (-linear - np.sqrt(discriminant)) / (2.0 * quadratic)
        candidates = [
            root
            for root in (root_1, root_2)
            if np.isfinite(root) and -1.0e-9 <= root <= dz + 1.0e-9
        ]
        distance = min(candidates) if candidates else root_1
    return float(np.clip(distance, 0.0, dz))


def energy_on_thermodynamic_grid(z, zm, em):
    """Use CLUBB's linear zm-to-zt interpolation weights."""
    z = np.asarray(z)
    zm = np.asarray(zm)
    em = np.asarray(em)
    if zm.size != z.size + 1:
        raise ValueError("Expected the momentum grid to have one extra level.")
    fraction_above = (z - zm[:-1]) / (zm[1:] - zm[:-1])
    return (1.0 - fraction_above) * em[:-1] + fraction_above * em[1:]


def compute_mixing_length(
    z,
    zm,
    thvm,
    thlm,
    rtm,
    em,
    pressure,
    exner,
    thv_ds,
    mu,
    lmin,
    *,
    lscale_max=LSCALE_MAX_STANDALONE,
):
    """Replicate ascending-grid ``compute_mixing_length`` and retain K paths.

    This is the standalone, Flatau-saturation path represented by CLUBB SCM
    statistics files.  The internal 0.1 m directional floor, endpoint
    envelopes, surface ``lmin`` floor, top boundary copy, and final cap are
    applied in the same order as the Fortran routine.
    """
    z = np.asarray(z, dtype=float)
    zm = np.asarray(zm, dtype=float)
    thvm = np.asarray(thvm, dtype=float)
    thlm = np.asarray(thlm, dtype=float)
    rtm = np.asarray(rtm, dtype=float)
    em = np.asarray(em, dtype=float)
    pressure = np.asarray(pressure, dtype=float)
    exner = np.asarray(exner, dtype=float)
    thv_ds = np.asarray(thv_ds, dtype=float)
    nzt = z.size
    if nzt < 3 or not np.all(np.diff(z) > 0.0):
        raise ValueError("The trajectory replica requires an ascending zt grid.")
    if abs(mu) < 1.0e-10:
        raise ValueError("Entrainment rate mu cannot be zero.")

    tke_i = energy_on_thermodynamic_grid(z, zm, em)
    raw_up = np.full(nzt, ZLMIN)
    raw_down = np.full(nzt, ZLMIN)
    upward_paths = []
    downward_paths = []

    # ---------------- Upwards length-scale calculation ----------------
    # As in Fortran, launches stop two levels below the upper boundary and
    # parcels can traverse only through the level immediately below it.
    for k in range(0, nzt - 2):
        energy = float(tke_i[k])
        path_altitude = [z[k]]
        path_energy = [energy]
        path_buoyancy = [0.0]
        parcel_thl = float(thlm[k])
        parcel_rt = float(rtm[k])
        path_thl = [parcel_thl]
        path_rt = [parcel_rt]
        buoyancy_previous = 0.0
        energy_endpoint = float(z[k])

        for j in range(k + 1, nzt - 1):
            dz = float(z[j] - z[j - 1])
            parcel_thl_previous = parcel_thl
            parcel_rt_previous = parcel_rt
            parcel_thl = _entrain_one_layer(
                parcel_thl, thlm[j - 1], thlm[j], dz, mu
            )
            parcel_rt = _entrain_one_layer(
                parcel_rt, rtm[j - 1], rtm[j], dz, mu
            )
            buoyancy = _parcel_buoyancy(
                j,
                parcel_thl,
                parcel_rt,
                thvm,
                pressure,
                exner,
                thv_ds,
            )
            cape_increment = 0.5 * (buoyancy_previous + buoyancy) * dz
            if energy + cape_increment <= 0.0:
                partial = _energy_root(
                    energy, buoyancy_previous, buoyancy, dz, 1.0
                )
                energy_endpoint = float(z[j - 1] + partial)
                path_altitude.append(energy_endpoint)
                path_energy.append(0.0)
                path_buoyancy.append(
                    buoyancy_previous
                    + (buoyancy - buoyancy_previous) * partial / dz
                )
                fraction = partial / dz
                path_thl.append(
                    _entrain_one_layer(
                        parcel_thl_previous,
                        thlm[j - 1],
                        thlm[j - 1] + fraction * (thlm[j] - thlm[j - 1]),
                        partial,
                        mu,
                    )
                )
                path_rt.append(
                    _entrain_one_layer(
                        parcel_rt_previous,
                        rtm[j - 1],
                        rtm[j - 1] + fraction * (rtm[j] - rtm[j - 1]),
                        partial,
                        mu,
                    )
                )
                raw_up[k] += energy_endpoint - z[k]
                break
            energy += cape_increment
            energy_endpoint = float(z[j])
            path_altitude.append(z[j])
            path_energy.append(energy)
            path_buoyancy.append(buoyancy)
            path_thl.append(parcel_thl)
            path_rt.append(parcel_rt)
            buoyancy_previous = buoyancy
        else:
            raw_up[k] += energy_endpoint - z[k]

        upward_paths.append(
            ParcelTrajectory(
                launch_index=k,
                altitude=np.asarray(path_altitude),
                energy=np.asarray(path_energy),
                buoyancy=np.asarray(path_buoyancy),
                parcel_thl=np.asarray(path_thl),
                parcel_rt=np.asarray(path_rt),
                energy_endpoint=energy_endpoint,
                raw_length=raw_up[k],
            )
        )

    # ---------------- Downwards length-scale calculation ----------------
    # The top boundary is a valid downward launch; the bottom boundary is not.
    for k in range(nzt - 1, 0, -1):
        energy = float(tke_i[k])
        path_altitude = [z[k]]
        path_energy = [energy]
        path_buoyancy = [0.0]
        parcel_thl = float(thlm[k])
        parcel_rt = float(rtm[k])
        path_thl = [parcel_thl]
        path_rt = [parcel_rt]
        buoyancy_previous = 0.0
        energy_endpoint = float(z[k])

        for j in range(k - 1, -1, -1):
            dz = float(z[j + 1] - z[j])
            parcel_thl_previous = parcel_thl
            parcel_rt_previous = parcel_rt
            parcel_thl = _entrain_one_layer(
                parcel_thl, thlm[j + 1], thlm[j], dz, mu
            )
            parcel_rt = _entrain_one_layer(
                parcel_rt, rtm[j + 1], rtm[j], dz, mu
            )
            buoyancy = _parcel_buoyancy(
                j,
                parcel_thl,
                parcel_rt,
                thvm,
                pressure,
                exner,
                thv_ds,
            )
            cape_increment = 0.5 * (buoyancy_previous + buoyancy) * dz
            if energy - cape_increment <= 0.0:
                partial = _energy_root(
                    energy, buoyancy_previous, buoyancy, dz, -1.0
                )
                energy_endpoint = float(z[j + 1] - partial)
                path_altitude.append(energy_endpoint)
                path_energy.append(0.0)
                path_buoyancy.append(
                    buoyancy_previous
                    + (buoyancy - buoyancy_previous) * partial / dz
                )
                fraction = partial / dz
                path_thl.append(
                    _entrain_one_layer(
                        parcel_thl_previous,
                        thlm[j + 1],
                        thlm[j + 1] + fraction * (thlm[j] - thlm[j + 1]),
                        partial,
                        mu,
                    )
                )
                path_rt.append(
                    _entrain_one_layer(
                        parcel_rt_previous,
                        rtm[j + 1],
                        rtm[j + 1] + fraction * (rtm[j] - rtm[j + 1]),
                        partial,
                        mu,
                    )
                )
                raw_down[k] += z[k] - energy_endpoint
                break
            energy -= cape_increment
            energy_endpoint = float(z[j])
            path_altitude.append(z[j])
            path_energy.append(energy)
            path_buoyancy.append(buoyancy)
            path_thl.append(parcel_thl)
            path_rt.append(parcel_rt)
            buoyancy_previous = buoyancy
        else:
            raw_down[k] += z[k] - energy_endpoint

        downward_paths.append(
            ParcelTrajectory(
                launch_index=k,
                altitude=np.asarray(path_altitude),
                energy=np.asarray(path_energy),
                buoyancy=np.asarray(path_buoyancy),
                parcel_thl=np.asarray(path_thl),
                parcel_rt=np.asarray(path_rt),
                energy_endpoint=energy_endpoint,
                raw_length=raw_down[k],
            )
        )

    # The nonlocal endpoint envelopes are the same cumulative scans that are
    # interleaved with the launch loops in Fortran.  Keeping a separate copy
    # lets the UI show what each parcel did before a neighboring launch wins.
    enveloped_up = raw_up.copy()
    lscale_up_max_alt = 0.0
    for k in range(0, nzt - 2):
        if z[k] + enveloped_up[k] < lscale_up_max_alt:
            enveloped_up[k] = lscale_up_max_alt - z[k]
        else:
            lscale_up_max_alt = z[k] + enveloped_up[k]

    enveloped_down = raw_down.copy()
    lscale_down_min_alt = float(z[-1])
    for k in range(nzt - 1, 0, -1):
        if z[k] - enveloped_down[k] > lscale_down_min_alt:
            enveloped_down[k] = z[k] - lscale_down_min_alt
        else:
            lscale_down_min_alt = z[k] - enveloped_down[k]

    # ---------------- Final length-scale calculation ----------------
    lminh = (
        np.maximum(0.0, LSCALE_SURFACE_LAYER_DEPTH - z)
        * lmin
        / LSCALE_SURFACE_LAYER_DEPTH
    )
    lscale_up = np.maximum(lminh, enveloped_up)
    lscale_down = np.maximum(lminh, enveloped_down)
    lscale = np.sqrt(lscale_up * lscale_down)
    lscale[-1] = lscale[-2]
    lscale = np.minimum(lscale, lscale_max)

    return MixingLengthResult(
        z=z,
        tke=tke_i,
        raw_up=raw_up,
        raw_down=raw_down,
        enveloped_up=enveloped_up,
        enveloped_down=enveloped_down,
        lscale_up=lscale_up,
        lscale_down=lscale_down,
        lscale=lscale,
        upward_paths=tuple(upward_paths),
        downward_paths=tuple(downward_paths),
    )


def _profile(variable, record, column, vertical_dimension):
    """Read one profile without assuming a particular dimension order."""
    index = []
    for dimension in variable.dimensions:
        if dimension == "time":
            index.append(record)
        elif dimension in {"col", "column", "columns", "ngrdcol"}:
            index.append(column)
        elif dimension == vertical_dimension:
            index.append(slice(None))
        else:
            index.append(0)
    return np.asarray(variable[tuple(index)], dtype=float).reshape(-1)


def _parameter_value(dataset, name, column, fallback):
    if not {"param_name", "clubb_params"} <= set(dataset.variables):
        return float(fallback)
    names = [
        item.strip()
        for item in nc.chartostring(dataset["param_name"][:]).tolist()
    ]
    if name not in names:
        return float(fallback)
    values = np.asarray(dataset["clubb_params"][names.index(name), :])
    return float(values[min(column, values.size - 1)])


def inspect_dataset(path):
    """Return selectable record/column metadata after validating a file."""
    path = Path(path).expanduser().resolve()
    with nc.Dataset(path) as dataset:
        missing = [
            name for name in REQUIRED_VARIABLES if name not in dataset.variables
        ]
        if missing:
            raise ValueError(
                "Missing required CLUBB fields: " + ", ".join(missing)
            )
        times = np.asarray(dataset["time"][:], dtype=float)
        if times.size == 0:
            raise ValueError(
                f"{path} contains no time records yet; wait for the writer or choose another file."
            )
        column_count = len(dataset.dimensions.get("col", (0,))) or 1
    return {
        "path": path,
        "times": times,
        "record_count": times.size,
        "column_count": column_count,
    }


def load_dataset_record(path, record, column=0):
    """Load one record plus tunable parameters embedded in the NetCDF file."""
    metadata = inspect_dataset(path)
    record = int(np.clip(record, 0, metadata["record_count"] - 1))
    column = int(np.clip(column, 0, metadata["column_count"] - 1))
    with nc.Dataset(metadata["path"]) as dataset:
        # The selected output may still be written by CLUBB. Re-check inside
        # the load handle so a truncate/recreate between inspection and this
        # read becomes a useful validation error instead of netCDF4 IndexError.
        available_records = int(dataset["time"].shape[0])
        if available_records == 0:
            raise ValueError(
                f"{metadata['path']} contains no time records yet; "
                "wait for the writer or choose another file."
            )
        record = int(np.clip(record, 0, available_records - 1))
        mu = _parameter_value(dataset, "mu", column, 1.0e-3)
        lmin_coef = _parameter_value(dataset, "lmin_coef", column, 0.5)
        return DatasetRecord(
            path=metadata["path"],
            record=record,
            column=column,
            time_seconds=float(dataset["time"][record]),
            z=np.asarray(dataset["zt"][:], dtype=float),
            zm=np.asarray(dataset["zm"][:], dtype=float),
            thlm=_profile(dataset["thlm_old"], record, column, "zt"),
            rtm=_profile(dataset["rtm_old"], record, column, "zt"),
            mean_thlm=_profile(dataset["thlm"], record, column, "zt"),
            mean_rtm=_profile(dataset["rtm"], record, column, "zt"),
            thvm=_profile(dataset["thvm"], record, column, "zt"),
            pressure=_profile(dataset["p_in_Pa"], record, column, "zt"),
            exner=_profile(dataset["exner"], record, column, "zt"),
            thv_ds=_profile(dataset["thv_ds_zt"], record, column, "zt"),
            em=_profile(dataset["em"], record, column, "zm"),
            fortran_lscale=_profile(
                dataset["Lscale"], record, column, "zt"
            ),
            fortran_up=_profile(
                dataset["Lscale_up"], record, column, "zt"
            ),
            fortran_down=_profile(
                dataset["Lscale_down"], record, column, "zt"
            ),
            mu=mu,
            lmin=lmin_coef * LMIN_DELTAZ,
        )


def compute_record(record):
    """Run the replica for a loaded NetCDF record."""
    return compute_mixing_length(
        record.z,
        record.zm,
        record.thvm,
        record.thlm,
        record.rtm,
        record.em,
        record.pressure,
        record.exner,
        record.thv_ds,
        record.mu,
        record.lmin,
    )


def profile_metrics(calculated, reference):
    """Return concise absolute-error metrics for one profile."""
    calculated = np.asarray(calculated, dtype=float)
    reference = np.asarray(reference, dtype=float)
    difference = calculated - reference
    finite = np.isfinite(calculated) & np.isfinite(reference)
    calculated_for_correlation = calculated[finite]
    reference_for_correlation = reference[finite]
    if (
        calculated_for_correlation.size < 2
        or reference_for_correlation.size < 2
        or float(np.std(calculated_for_correlation)) == 0.0
        or float(np.std(reference_for_correlation)) == 0.0
    ):
        correlation = float("nan")
    else:
        correlation = float(
            np.corrcoef(calculated_for_correlation, reference_for_correlation)[0, 1]
        )
    return {
        "rmse": float(np.sqrt(np.mean(difference**2))),
        "mae": float(np.mean(np.abs(difference))),
        "max_abs": float(np.max(np.abs(difference))),
        "correlation": correlation,
    }


def discover_netcdf_file_records(repo_root):
    """Return trajectory-compatible stats files with shared unique labels."""
    root = Path(repo_root).resolve()
    compatible = []
    # Keep the established ``output*`` root policy, but share traversal and
    # stable per-file identity with other output-oriented views.
    roots = sorted(path for path in root.glob("output*") if path.is_dir())
    for record in discover_stats_files(roots, exclude_dir_names=()):
        path = Path(record["path"])
        try:
            inspect_dataset(path)
        except (OSError, ValueError):
            continue
        compatible.append({**record, "path": str(path.resolve())})
    return compatible


def discover_netcdf_files(repo_root):
    """Compatibility wrapper returning only trajectory-compatible paths."""
    return [Path(record["path"]) for record in discover_netcdf_file_records(repo_root)]
