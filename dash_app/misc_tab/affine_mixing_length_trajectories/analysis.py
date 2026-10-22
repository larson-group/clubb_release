"""Recompute parcel travel and coherent mixing lengths from a saved column.

This NumPy diagnostic supplies the Dash Affine Coherence Mixing Length page.
It follows parcels through a fixed environmental profile, computes the greatest
adverse buoyancy-work barrier along each path, and uses a launch-energy survival
weight to diagnose travel distance. An affine layer update multiplies the
incoming parcel property by a factor and adds an environmental contribution.

Usage:
  - ``affine_parcel_matrix`` evaluates one property for every launch/target pair.
  - ``compute_affine_parcel_fields`` constructs ``AffineParcelFields`` from
    profiles; ``compute_record_fields`` takes a loaded DatasetRecord instead.
  - ``compute_required_energy_barriers`` returns path barriers for those fields.
  - ``inverse_square_coherence_fraction`` evaluates the survival weight directly.
  - ``compute_inverse_square_mixing_length`` returns a
    ``CoherenceMixingLengthResult`` using one of ``OVERWRITE_MODES``.
  File discovery, loading, and profile metrics are reused from the neighboring
  ``mixing_length_trajectories.analysis`` module.

Understanding:
  Grids must ascend and contain at least three thermodynamic levels. Inputs
  describe one column/record: heights [m], temperatures [K], water [kg/kg],
  pressure [Pa], momentum-level TKE [m^2/s^2], and entrainment rate [1/m].
  Dense matrices have axes (launch, target), with upward paths above the
  diagonal in plotted altitude coordinates. They require quadratic storage.
  Directional launch properties are separate; displayed property diagonals
  retain environmental means, and the buoyancy diagonal is zero. Energy
  continues beyond deterministic stops for display; barriers sample endpoints.

  The reducer includes the first layer with mean coherence <= 0.01, then masks
  the tail. It excludes the lid launch from pooled averages. Raw lengths have
  a 0.1 m numerical minimum; final lengths also apply the surface floor,
  directional geometric mean, cap, and copied combined-length lid value.
  No model timestep is run and no dataset is modified. Saved statistics may
  be time averages, so recomputing from them is not an exact timestep replay
  or a Fortran/JAX regression comparison.
  The shared plotting saturation helper uses a factored Flatau polynomial
  without the production parcel helper's low-pressure cap. It therefore differs
  when pressure minus saturation vapor pressure is below 1 Pa.

References:
  ``doc/coherence_mixing_length.md`` and
  ``src/CLUBB_core/mixing_length.F90`` define the method and active bounds.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from dash_app.misc_tab.mixing_length_trajectories.analysis import (
    CP,
    EP1,
    EP2,
    GRAV,
    LSCALE_MAX_STANDALONE,
    LSCALE_SURFACE_LAYER_DEPTH,
    LV,
    LV2_COEF,
    ZLMIN,
    DatasetRecord,
    _energy_root,
    discover_netcdf_file_records,
    discover_netcdf_files,
    energy_on_thermodynamic_grid,
    inspect_dataset,
    load_dataset_record,
    profile_metrics,
    saturation_mixing_ratio_flatau,
)

COHERENCE_CUTOFF = 0.01
NO_OVERWRITE = "no_overwrite"
LEGACY_MAX_OVERWRITE = "legacy_max"
COHERENCE_AVERAGED_OVERWRITE = "coherence_averaged"
OVERWRITE_MODES = (
    NO_OVERWRITE,
    LEGACY_MAX_OVERWRITE,
    COHERENCE_AVERAGED_OVERWRITE,
)


@dataclass(frozen=True)
class AffineParcelFields:
    """Store one ascending column and its dense parcel diagnostics.

    Matrices use (launch, target) axes: thl [K], rt [kg/kg], buoyancy [m/s^2],
    and energy [m^2/s^2]. Directional launch buoyancies are upward-positive.
    Stop-altitude vectors [m] describe deterministic exhaustion along these
    same paths; untracked launch levels are NaN.
    """

    z: np.ndarray
    tke: np.ndarray
    mu: float
    parcel_thl_matrix: np.ndarray
    parcel_rt_matrix: np.ndarray
    buoyancy_matrix: np.ndarray
    launch_buoyancy_up: np.ndarray
    launch_buoyancy_down: np.ndarray
    energy_matrix: np.ndarray
    legacy_upward_stops: np.ndarray
    legacy_downward_stops: np.ndarray


@dataclass(frozen=True)
class CoherenceMixingLengthResult:
    """Store one survival-weight reduction of the parcel diagnostics.

    Barrier, endpoint coherence, and layer-integral matrices use (launch,target)
    axes and units [m^2/s^2], [-], and [m]. Raw directional lengths already
    include the 0.1 m numerical minimum; final lengths also include the selected
    aggregation and surface floor. All length vectors and z are in m.
    """

    z: np.ndarray
    tke: np.ndarray
    overwrite_mode: str
    tke_multiplier: float
    up_weight: float
    barrier_matrix: np.ndarray
    coherence_matrix: np.ndarray
    layer_integrals: np.ndarray
    raw_up: np.ndarray
    raw_down: np.ndarray
    raw_lscale: np.ndarray
    lscale_up: np.ndarray
    lscale_down: np.ndarray
    lscale: np.ndarray


def _layer_affine_constant(environment_start, environment_end, dz, mu):
    """Return retention and environmental contribution for one entraining layer.

    Environment endpoints use the same units as the parcel property; positive
    dz [m] and mu [1/m] define the layer-linear analytic entrainment update.
    """
    retention = np.exp(-mu * dz)
    response = -np.expm1(-mu * dz) / (mu * dz)
    constant = (
        environment_end
        - retention * environment_start
        - response * (environment_end - environment_start)
    )
    return retention, constant


def _affine_prefix(environment, z, mu, direction):
    """Return the zero-boundary parcel solution along an up/down traversal.

    The prefix and environmental vector have the same property units; z is an
    ascending height vector [m]. Invalid direction names raise ValueError.
    """
    prefix = np.zeros_like(environment)
    if direction == "up":
        targets = range(1, z.size)
        offset = -1
    elif direction == "down":
        targets = range(z.size - 2, -1, -1)
        offset = 1
    else:
        raise ValueError(f"Unknown affine traversal direction: {direction}")
    for target in targets:
        source = target + offset
        retention, constant = _layer_affine_constant(
            environment[source],
            environment[target],
            abs(z[target] - z[source]),
            mu,
        )
        prefix[target] = retention * prefix[source] + constant
    return prefix


def affine_parcel_matrix(
    environment,
    z,
    mu,
    upward_launch=None,
    downward_launch=None,
):
    """Return one parcel-property matrix with axes (launch, target).

    Environment, upward_launch, and downward_launch are equal-length profiles
    in the same units; omitted launches use the environment. z ascends in m
    and mu is nonzero [1/m]. Analytic affine updates solve entrainment through
    a layer-linear environment. The matrix diagonal retains the environment.
    """

    environment = np.asarray(environment, dtype=float)
    z = np.asarray(z, dtype=float)
    if environment.shape != z.shape:
        raise ValueError("Environment and z must have the same shape.")
    if abs(mu) < 1.0e-10:
        raise ValueError("Entrainment rate mu cannot be zero.")
    upward_launch = np.asarray(
        environment if upward_launch is None else upward_launch, dtype=float
    )
    downward_launch = np.asarray(
        environment if downward_launch is None else downward_launch, dtype=float
    )
    if upward_launch.shape != z.shape or downward_launch.shape != z.shape:
        raise ValueError("Directional launch profiles must match z.")

    prefix_up = _affine_prefix(environment, z, mu, "up")
    prefix_down = _affine_prefix(environment, z, mu, "down")
    retention = np.exp(-mu * np.abs(z[:, None] - z[None, :]))
    upward = prefix_up[None, :] + retention * (
        upward_launch - prefix_up
    )[:, None]
    downward = prefix_down[None, :] + retention * (
        downward_launch - prefix_down
    )[:, None]
    launch, target = np.indices((z.size, z.size))
    matrix = np.where(target > launch, upward, downward)
    np.fill_diagonal(matrix, environment)
    return matrix


def _buoyancy_values(thl_parcel, rt_parcel, thvm, pressure, exner, thv_ds):
    """Return upward-positive parcel buoyancy [m/s^2] for broadcastable profiles.

    Use liquid-water potential temperature [K], total water [kg/kg], and target
    environmental thvm/thv_ds [K], pressure [Pa], and dimensionless Exner factor.
    Liquid water follows the local Flatau saturation and saturation-excess formula.
    """
    liquid_temperature = thl_parcel * exner
    saturation_rt = saturation_mixing_ratio_flatau(pressure, liquid_temperature)
    liquid_temperature_sqd = liquid_temperature**2
    liquid_water = np.maximum(
        (rt_parcel - saturation_rt)
        * liquid_temperature_sqd
        / (liquid_temperature_sqd + LV2_COEF * saturation_rt),
        0.0,
    )
    latent_coefficient = LV / (exner * CP) - EP2 * thv_ds
    thv_parcel = (
        thl_parcel
        + EP1 * thv_ds * rt_parcel
        + latent_coefficient * liquid_water
    )
    return GRAV * (thv_parcel - thvm) / thvm


def _energy_matrix(z, tke, buoyancy, launch_up, launch_down):
    """Return TKE plus trapezoidal work [m^2/s^2] for every launch/target pair.

    z ascends in m; all buoyancies are upward-positive [m/s^2]. Energy starts
    from launch tke, adds upward work, and subtracts downward work. Paths
    continue past zero energy so the display can show the adverse-work tail.
    """
    energy = np.empty_like(buoyancy)
    np.fill_diagonal(energy, tke)
    for launch in range(z.size):
        parcel_energy = float(tke[launch])
        previous_buoyancy = float(launch_up[launch])
        for target in range(launch + 1, z.size):
            dz = float(z[target] - z[target - 1])
            target_buoyancy = float(buoyancy[launch, target])
            parcel_energy += 0.5 * (
                previous_buoyancy + target_buoyancy
            ) * dz
            energy[launch, target] = parcel_energy
            previous_buoyancy = target_buoyancy

        parcel_energy = float(tke[launch])
        previous_buoyancy = float(launch_down[launch])
        for target in range(launch - 1, -1, -1):
            dz = float(z[target + 1] - z[target])
            target_buoyancy = float(buoyancy[launch, target])
            parcel_energy -= 0.5 * (
                previous_buoyancy + target_buoyancy
            ) * dz
            energy[launch, target] = parcel_energy
            previous_buoyancy = target_buoyancy
    return energy


def _legacy_stop_altitudes(
    z,
    tke,
    buoyancy,
    launch_up,
    launch_down,
):
    """Return deterministic first-zero-energy heights [m] along the supplied paths.

    Use ascending z, launch TKE, and upward-positive endpoint/launch buoyancies.
    Resolve the first impassable layer with the shared energy root. NaN marks
    untracked launches. These are stopping references for the supplied launch
    states, not a rerun of all legacy Fortran launch and overwrite behavior.
    """

    upward = np.full(z.size, np.nan)
    downward = np.full(z.size, np.nan)
    for launch in range(0, z.size - 2):
        parcel_energy = float(tke[launch])
        previous_buoyancy = float(launch_up[launch])
        endpoint = float(z[launch])
        for target in range(launch + 1, z.size - 1):
            dz = float(z[target] - z[target - 1])
            target_buoyancy = float(buoyancy[launch, target])
            energy_change = 0.5 * (
                previous_buoyancy + target_buoyancy
            ) * dz
            if parcel_energy + energy_change <= 0.0:
                endpoint = float(
                    z[target - 1]
                    + _energy_root(
                        parcel_energy,
                        previous_buoyancy,
                        target_buoyancy,
                        dz,
                        1.0,
                    )
                )
                break
            parcel_energy += energy_change
            endpoint = float(z[target])
            previous_buoyancy = target_buoyancy
        upward[launch] = endpoint

    for launch in range(z.size - 1, 0, -1):
        parcel_energy = float(tke[launch])
        previous_buoyancy = float(launch_down[launch])
        endpoint = float(z[launch])
        for target in range(launch - 1, -1, -1):
            dz = float(z[target + 1] - z[target])
            target_buoyancy = float(buoyancy[launch, target])
            energy_change = 0.5 * (
                previous_buoyancy + target_buoyancy
            ) * dz
            if parcel_energy - energy_change <= 0.0:
                endpoint = float(
                    z[target + 1]
                    - _energy_root(
                        parcel_energy,
                        previous_buoyancy,
                        target_buoyancy,
                        dz,
                        -1.0,
                    )
                )
                break
            parcel_energy -= energy_change
            endpoint = float(z[target])
            previous_buoyancy = target_buoyancy
        downward[launch] = endpoint
    return upward, downward


def compute_affine_parcel_fields(
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
    thl_launch_up=None,
    rt_launch_up=None,
    thl_launch_down=None,
    rt_launch_down=None,
):
    """Return parcel properties, buoyancy, energy, and stopping references.

    Inputs are one ascending column: z/zm [m], thvm/thlm/thv_ds [K], rtm [kg/kg],
    em on zm [m^2/s^2], pressure [Pa], dimensionless exner, and mu [1/m].
    Optional thl/rt launch vectors select different upward/downward parcels;
    otherwise use environmental means. Return AffineParcelFields with dense
    (launch,target) matrices; a short/nonascending grid or near-zero mu raises.
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
    if z.size < 3 or not np.all(np.diff(z) > 0.0):
        raise ValueError("The affine lab requires an ascending zt grid.")
    if abs(mu) < 1.0e-10:
        raise ValueError("Entrainment rate mu cannot be zero.")

    thl_launch_up = np.asarray(
        thlm if thl_launch_up is None else thl_launch_up, dtype=float
    )
    rt_launch_up = np.asarray(
        rtm if rt_launch_up is None else rt_launch_up, dtype=float
    )
    thl_launch_down = np.asarray(
        thlm if thl_launch_down is None else thl_launch_down, dtype=float
    )
    rt_launch_down = np.asarray(
        rtm if rt_launch_down is None else rt_launch_down, dtype=float
    )
    parcel_thl = affine_parcel_matrix(
        thlm, z, mu, thl_launch_up, thl_launch_down
    )
    parcel_rt = affine_parcel_matrix(
        rtm, z, mu, rt_launch_up, rt_launch_down
    )
    buoyancy = _buoyancy_values(
        parcel_thl,
        parcel_rt,
        thvm[None, :],
        pressure[None, :],
        exner[None, :],
        thv_ds[None, :],
    )
    np.fill_diagonal(buoyancy, 0.0)
    launch_buoyancy_up = _buoyancy_values(
        thl_launch_up, rt_launch_up, thvm, pressure, exner, thv_ds
    )
    launch_buoyancy_down = _buoyancy_values(
        thl_launch_down, rt_launch_down, thvm, pressure, exner, thv_ds
    )
    tke = energy_on_thermodynamic_grid(z, zm, em)
    upward_stops, downward_stops = _legacy_stop_altitudes(
        z,
        tke,
        buoyancy,
        launch_buoyancy_up,
        launch_buoyancy_down,
    )
    return AffineParcelFields(
        z=z,
        tke=tke,
        mu=float(mu),
        parcel_thl_matrix=parcel_thl,
        parcel_rt_matrix=parcel_rt,
        buoyancy_matrix=buoyancy,
        launch_buoyancy_up=launch_buoyancy_up,
        launch_buoyancy_down=launch_buoyancy_down,
        energy_matrix=_energy_matrix(
            z,
            tke,
            buoyancy,
            launch_buoyancy_up,
            launch_buoyancy_down,
        ),
        legacy_upward_stops=upward_stops,
        legacy_downward_stops=downward_stops,
    )


def compute_record_fields(record: DatasetRecord, pdf_launch_frac=None):
    """Return parcel fields from one loaded DatasetRecord.

    Blend environmental thl/rt toward PDF component 1 upward and component 2
    downward. The optional fraction overrides the record's saved launch fraction;
    it does not change the record. Mixture fraction is not an additional weight.
    """

    launch_fraction = float(
        record.pdf_launch_frac if pdf_launch_frac is None else pdf_launch_frac
    )
    return compute_affine_parcel_fields(
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
        record.thlm + launch_fraction * (record.thl_1 - record.thlm),
        record.rtm + launch_fraction * (record.rt_1 - record.rtm),
        record.thlm + launch_fraction * (record.thl_2 - record.thlm),
        record.rtm + launch_fraction * (record.rt_2 - record.rtm),
    )


def compute_required_energy_barriers(fields):
    """Return the greatest adverse-work deficit at every launch/target [m^2/s^2].

    Accumulate trapezoidal work from directional launch buoyancies in fields.
    Running barriers start at zero, never decrease along a path, and sample
    level endpoints only; they do not search for within-layer energy extrema.
    """

    z = fields.z
    buoyancy = fields.buoyancy_matrix
    barriers = np.zeros_like(buoyancy)
    for launch in range(z.size):
        work = 0.0
        max_barrier = 0.0
        previous_buoyancy = float(fields.launch_buoyancy_up[launch])
        for target in range(launch + 1, z.size):
            dz = float(z[target] - z[target - 1])
            target_buoyancy = float(buoyancy[launch, target])
            work += 0.5 * (previous_buoyancy + target_buoyancy) * dz
            max_barrier = max(max_barrier, -work)
            barriers[launch, target] = max_barrier
            previous_buoyancy = target_buoyancy

        work = 0.0
        max_barrier = 0.0
        previous_buoyancy = -float(fields.launch_buoyancy_down[launch])
        for target in range(launch - 1, -1, -1):
            dz = float(z[target + 1] - z[target])
            target_buoyancy = -float(buoyancy[launch, target])
            work += 0.5 * (previous_buoyancy + target_buoyancy) * dz
            max_barrier = max(max_barrier, -work)
            barriers[launch, target] = max_barrier
            previous_buoyancy = target_buoyancy
    return barriers


def inverse_square_coherence_fraction(barrier, energy_scale):
    """Return the dimensionless survival weight (E/(E+H))^2.

    barrier H and energy_scale E [m^2/s^2] broadcast together. At zero barrier
    return one even for E=0; a positive barrier and nonpositive E return zero.
    The path reducer separately disables all travel for zero-energy launches.
    """

    barrier, energy_scale = np.broadcast_arrays(
        np.asarray(barrier, dtype=float),
        np.asarray(energy_scale, dtype=float),
    )
    coherence = np.ones_like(barrier)
    challenged = barrier > 0.0
    positive = challenged & (energy_scale > 0.0)
    coherence[challenged] = 0.0
    coherence[positive] = (
        energy_scale[positive]
        / (energy_scale[positive] + barrier[positive])
    ) ** 2
    return coherence


def _full_layer_integrals(fields, barriers, energy_scale):
    """Return untruncated layer travel integrals [m] and layer depths [m].

    For each launch/target in barriers, integrate survival assuming the barrier
    is linear between endpoints. The diagonal integral is zero. Zero-barrier
    layers get unit mean coherence, including zero-energy cases later masked.
    """
    z = fields.z
    launch, target = np.indices(barriers.shape)
    upward = target > launch
    downward = target < launch
    dz = np.ones_like(barriers)
    barrier_start = np.zeros_like(barriers)
    up_launch, up_target = np.where(upward)
    up_previous = up_target - 1
    dz[up_launch, up_target] = z[up_target] - z[up_previous]
    barrier_start[up_launch, up_target] = barriers[up_launch, up_previous]
    down_launch, down_target = np.where(downward)
    down_previous = down_target + 1
    dz[down_launch, down_target] = z[down_previous] - z[down_target]
    barrier_start[down_launch, down_target] = barriers[
        down_launch, down_previous
    ]
    scale = energy_scale[:, None]
    denominator = (scale + barrier_start) * (scale + barriers)
    average_coherence = np.zeros_like(barriers)
    average_coherence[(barrier_start == 0.0) & (barriers == 0.0)] = 1.0
    np.divide(
        scale**2,
        denominator,
        out=average_coherence,
        where=(scale > 0.0) & (denominator > 0.0),
    )
    return np.where(upward | downward, dz * average_coherence, 0.0), dz


def _truncate_coherence_paths(fields, barriers, energy_scale):
    """Return endpoint weights and layer integrals with the 1% tail cutoff.

    fields supplies the grid, barriers the full endpoint deficits [m^2/s^2],
    and energy_scale one value per launch [m^2/s^2]. Zero-energy launches have
    no travel; the local diagonal retains unit bookkeeping weight.
    """

    full_integrals, dz = _full_layer_integrals(fields, barriers, energy_scale)
    full_coherence = inverse_square_coherence_fraction(
        barriers, energy_scale[:, None]
    )
    return _truncate_paths(
        fields,
        full_coherence,
        full_integrals,
        dz,
        energy_scale > 0.0,
    )


def _truncate_paths(fields, full_coherence, full_integrals, dz, active_launch):
    """Retain the first layer at mean coherence <= 0.01, then discard its tail.

    Read full endpoint weights, layer integrals [m], depths [m], and active
    launch flags. Return truncated (launch,target) weights/integrals with unit
    diagonal weight; upward paths end before the lid, downward paths at ground.
    """

    coherence = np.zeros_like(full_coherence)
    integrals = np.zeros_like(full_integrals)
    np.fill_diagonal(coherence, 1.0)
    for launch in range(fields.z.size):
        if not active_launch[launch]:
            continue
        for target in range(launch + 1, fields.z.size - 1):
            coherence[launch, target] = full_coherence[launch, target]
            integrals[launch, target] = full_integrals[launch, target]
            if integrals[launch, target] / dz[launch, target] <= COHERENCE_CUTOFF:
                break
        for target in range(launch - 1, -1, -1):
            coherence[launch, target] = full_coherence[launch, target]
            integrals[launch, target] = full_integrals[launch, target]
            if integrals[launch, target] / dz[launch, target] <= COHERENCE_CUTOFF:
                break
    return coherence, integrals


def _raw_directional_lengths(layer_integrals):
    """Sum retained layer integrals by direction; return launch lengths [m] floored at 0.1 m."""
    nzt = layer_integrals.shape[0]
    launch, target = np.indices((nzt, nzt))
    upward = np.sum(
        np.where((target > launch) & (target < nzt - 1), layer_integrals, 0.0),
        axis=1,
    )
    downward = np.sum(
        np.where(target < launch, layer_integrals, 0.0), axis=1
    )
    return np.maximum(ZLMIN, upward), np.maximum(ZLMIN, downward)


def _legacy_endpoint_overwrite(z, raw_up, raw_down):
    """Return directional lengths [m] extended by neighboring maximum/minimum endpoints."""
    upward = raw_up.copy()
    maximum_endpoint = 0.0
    for launch in range(0, z.size - 2):
        endpoint = z[launch] + upward[launch]
        if endpoint < maximum_endpoint:
            upward[launch] = maximum_endpoint - z[launch]
        else:
            maximum_endpoint = endpoint
    downward = raw_down.copy()
    minimum_endpoint = float(z[-1])
    for launch in range(z.size - 1, 0, -1):
        endpoint = z[launch] - downward[launch]
        if endpoint > minimum_endpoint:
            downward[launch] = z[launch] - minimum_endpoint
        else:
            minimum_endpoint = endpoint
    return upward, downward


def _coherence_weighted_averages(z, coherence, raw_up, raw_down):
    """Return arrival-weighted remaining upward/downward distances [m].

    Use retained endpoint coherence to reconstruct layer integrals and remaining
    reach at each target. z ascends in m; raw lengths initialize boundary values.
    Only active launches below the lid contribute to pooled averages, which
    receive the 0.1 m numerical minimum.
    """

    nzt = z.size
    pooled_up = raw_up.copy()
    coherence_numerator = np.zeros(nzt)
    coherence_denominator = np.zeros(nzt)
    for launch in range(0, nzt - 1):
        remaining_length = 0.0
        for target in range(nzt - 2, launch - 1, -1):
            target_coherence = coherence[launch, target]
            coherence_numerator[target] += remaining_length
            coherence_denominator[target] += target_coherence
            if target != launch:
                previous = target - 1
                remaining_length += np.sqrt(
                    coherence[launch, previous] * coherence[launch, target]
                ) * (z[target] - z[previous])
    pooled_up[:-1] = np.maximum(
        ZLMIN,
        coherence_numerator[:-1] / coherence_denominator[:-1],
    )
    pooled_down = raw_down.copy()
    coherence_numerator.fill(0.0)
    coherence_denominator.fill(0.0)
    for launch in range(0, nzt - 1):
        remaining_length = 0.0
        for target in range(0, launch + 1):
            target_coherence = coherence[launch, target]
            coherence_numerator[target] += remaining_length
            coherence_denominator[target] += target_coherence
            if target != launch:
                previous = target + 1
                remaining_length += np.sqrt(
                    coherence[launch, target] * coherence[launch, previous]
                ) * (z[previous] - z[target])
    pooled_down[:-1] = np.maximum(
        ZLMIN,
        coherence_numerator[:-1] / coherence_denominator[:-1],
    )
    return pooled_up, pooled_down


def _finalize(z, upward, downward, lmin, up_weight, lscale_max, *, floor):
    """Return bounded directional/combined length profiles [m] without mutating inputs.

    Optionally apply the surface floor from lmin before the geometric mean
    with upward exponent up_weight. Cap only the combined length at lscale_max
    and copy its highest active value to the lid.
    """
    upward = upward.copy()
    downward = downward.copy()
    if floor:
        surface_floor = (
            np.maximum(0.0, LSCALE_SURFACE_LAYER_DEPTH - z)
            * lmin
            / LSCALE_SURFACE_LAYER_DEPTH
        )
        upward = np.maximum(surface_floor, upward)
        downward = np.maximum(surface_floor, downward)
    lscale = upward**up_weight * downward ** (1.0 - up_weight)
    lscale = np.minimum(lscale, lscale_max)
    lscale[-1] = lscale[-2]
    return upward, downward, lscale


def _reduce_coherence_mixing_length(
    fields,
    lmin,
    barriers,
    coherence,
    layer_integrals,
    multiplier,
    up_weight,
    lscale_max,
    overwrite_mode,
):
    """Return one result after the selected aggregation and directional bounds.

    fields/barriers supply the common parcel solution; coherence and
    layer_integrals supply retained weights/travel. Apply numerical floors,
    the overwrite_mode, surface floor lmin [m], up_weight, and lscale_max [m].
    Unknown aggregation names raise ValueError.
    """

    raw_up, raw_down = _raw_directional_lengths(layer_integrals)
    _, _, raw_lscale = _finalize(
        fields.z,
        raw_up,
        raw_down,
        lmin,
        up_weight,
        lscale_max,
        floor=False,
    )
    if overwrite_mode == NO_OVERWRITE:
        upward, downward = raw_up, raw_down
    elif overwrite_mode == LEGACY_MAX_OVERWRITE:
        upward, downward = _legacy_endpoint_overwrite(
            fields.z, raw_up, raw_down
        )
    elif overwrite_mode == COHERENCE_AVERAGED_OVERWRITE:
        upward, downward = _coherence_weighted_averages(
            fields.z, coherence, raw_up, raw_down
        )
    else:
        raise ValueError(f"Unknown overwrite method: {overwrite_mode}")
    lscale_up, lscale_down, lscale = _finalize(
        fields.z,
        upward,
        downward,
        lmin,
        up_weight,
        lscale_max,
        floor=True,
    )
    return CoherenceMixingLengthResult(
        z=fields.z,
        tke=fields.tke,
        overwrite_mode=overwrite_mode,
        tke_multiplier=multiplier,
        up_weight=up_weight,
        barrier_matrix=barriers,
        coherence_matrix=coherence,
        layer_integrals=layer_integrals,
        raw_up=raw_up,
        raw_down=raw_down,
        raw_lscale=raw_lscale,
        lscale_up=lscale_up,
        lscale_down=lscale_down,
        lscale=lscale,
    )


def compute_inverse_square_mixing_length(
    fields,
    lmin,
    tke_multiplier=1.0,
    *,
    barriers=None,
    up_weight=0.5,
    lscale_max=LSCALE_MAX_STANDALONE,
    overwrite_mode=COHERENCE_AVERAGED_OVERWRITE,
):
    """Return coherent mixing lengths for one fixed parcel solution.

    fields supplies the ascending grid and TKE; optional barriers reuse a prior
    endpoint-barrier calculation. A finite nonnegative tke_multiplier sets mean
    launch energy/TKE; up_weight must lie in [0,1]. lmin and lscale_max are in m.
    The default overwrite_mode pools remaining reach by arrival coherence.
    Return CoherenceMixingLengthResult, preserving all supplied fields/arrays.
    """

    multiplier = float(tke_multiplier)
    up_weight = float(up_weight)
    if not np.isfinite(multiplier) or multiplier < 0.0:
        raise ValueError("TKE multiplier must be finite and nonnegative.")
    if not 0.0 <= up_weight <= 1.0:
        raise ValueError("Upward mixing-length weight must be in [0, 1].")
    barriers = (
        compute_required_energy_barriers(fields)
        if barriers is None
        else np.asarray(barriers, dtype=float)
    )
    energy_scale = multiplier * np.maximum(fields.tke, 0.0)
    coherence, layer_integrals = _truncate_coherence_paths(
        fields, barriers, energy_scale
    )
    return _reduce_coherence_mixing_length(
        fields,
        lmin,
        barriers,
        coherence,
        layer_integrals,
        multiplier,
        up_weight,
        lscale_max,
        overwrite_mode,
    )


__all__ = [
    "AffineParcelFields",
    "COHERENCE_AVERAGED_OVERWRITE",
    "CoherenceMixingLengthResult",
    "LEGACY_MAX_OVERWRITE",
    "NO_OVERWRITE",
    "OVERWRITE_MODES",
    "affine_parcel_matrix",
    "compute_affine_parcel_fields",
    "compute_inverse_square_mixing_length",
    "compute_record_fields",
    "compute_required_energy_barriers",
    "discover_netcdf_file_records",
    "discover_netcdf_files",
    "inspect_dataset",
    "inverse_square_coherence_fraction",
    "load_dataset_record",
    "profile_metrics",
]
