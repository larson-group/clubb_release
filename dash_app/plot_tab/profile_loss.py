import numpy as np

from utilities.loss_metrics import calculate_scaled_rmse


def _normalize_term(var_name, benchmark_values):
    norm = float(np.max(benchmark_values) - np.min(benchmark_values))
    if var_name == "cloud_frac":
        return max(0.01, norm), None
    if var_name == "rcm":
        return max(1.0e-6, norm), None
    if abs(norm) < np.finfo(float).eps:
        return None, f"{var_name} benchmark profile is flat over the selected height range"
    return norm, None


def _prepare_interp_profile(z_values, profile):
    z_arr = np.asarray(z_values, dtype=float)
    profile_arr = np.asarray(profile, dtype=float)
    finite = np.isfinite(z_arr) & np.isfinite(profile_arr)
    if not np.any(finite):
        return None, None
    z_arr = z_arr[finite]
    profile_arr = profile_arr[finite]
    order = np.argsort(z_arr)
    z_arr = z_arr[order]
    profile_arr = profile_arr[order]
    unique_z, unique_indices = np.unique(z_arr, return_index=True)
    return unique_z, profile_arr[unique_indices]


def compute_profile_loss(var_name, clubb_profile, clubb_z, benchmark_profile, benchmark_z, height_range):
    clubb_z_arr = np.asarray(clubb_z, dtype=float)
    clubb_profile_arr = np.asarray(clubb_profile, dtype=float)
    finite_clubb = np.isfinite(clubb_z_arr) & np.isfinite(clubb_profile_arr)
    if not np.any(finite_clubb):
        return None, "No finite CLUBB profile data are available"

    clubb_z_arr = clubb_z_arr[finite_clubb]
    clubb_profile_arr = clubb_profile_arr[finite_clubb]

    benchmark_z_arr, benchmark_profile_arr = _prepare_interp_profile(benchmark_z, benchmark_profile)
    if benchmark_z_arr is None or benchmark_z_arr.size == 0:
        return None, "No finite benchmark profile data are available"

    low = float(min(height_range))
    high = float(max(height_range))
    within_height = (clubb_z_arr >= low) & (clubb_z_arr <= high)
    within_benchmark = (clubb_z_arr >= benchmark_z_arr[0]) & (clubb_z_arr <= benchmark_z_arr[-1])
    selected = within_height & within_benchmark
    if not np.any(selected):
        return None, "No overlapping plotted levels fall inside the selected height range"

    target_z = clubb_z_arr[selected]
    clubb_values = clubb_profile_arr[selected]
    benchmark_values = np.interp(target_z, benchmark_z_arr, benchmark_profile_arr)

    norm_term, norm_error = _normalize_term(var_name, benchmark_values)
    if norm_error is not None:
        return None, norm_error

    loss = calculate_scaled_rmse(clubb_values, benchmark_values, norm=norm_term)
    return loss, None
