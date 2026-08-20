"""Worker process loop for multi-case CLUBB tuning."""

from __future__ import annotations

import os
from pathlib import Path
import multiprocessing
import sys
import traceback

import numpy as np

from tuner.paths import REPO_ROOT, RUN_SCRIPTS
from tuner.taylor_metrics import INVALID_SCALED_RMSE_PENALTY, LOSS_METRIC_NAMES

if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from utilities.create_case_namelist import create_loss_case_namelist  # noqa: E402


def _actionable_error_message(exc: BaseException) -> str:
    """Return an error string with the common module workaround when relevant."""
    message = f"{exc}\n{traceback.format_exc(limit=10)}"
    lower = message.lower()
    if (
        "f2py_clubb_get_loss_for_params" in lower
        and "expected 3" in lower
    ):
        message += (
            "\nThe loaded CLUBB Python API does not expose Taylor-metric loss outputs. "
            "Rebuild the CLUBB Python API so f2py_clubb_get_loss_for_params returns "
            "scaled_rmse, correlation, std_ratio, centered_rmse_norm, and bias_norm."
        )
    elif "libnvomp" in lower or "netcdf" in lower or "no module named" in lower or "importerror" in lower:
        message += (
            "\nCLUBB Python import/init failed. Check that the compiler and "
            "NetCDF modules used for the Python build are available in this shell."
        )
    return message


def _baseline_worker_main(conn, init_payload: dict) -> None:
    """Evaluate one-column baselines in an isolated CLUBB process."""
    case_name = init_payload.get("case_name")
    selected_fields = list(init_payload.get("selected_fields", []))
    num_time_windows = int(init_payload.get("num_time_windows", 1))
    case_defaults = dict(init_payload.get("case_defaults") or {})
    config = init_payload.get("config") or "default"
    override = init_payload.get("override") or None
    worker_dir = Path(init_payload.get("worker_dir")).resolve()
    original_cwd = os.getcwd()

    try:
        from clubb_python import clubb_api

        baselines = {}
        baseline_specs = [("clubb_default", "default", None)]
        if override:
            baseline_specs.append(("override_defaults", config, override))
        for baseline_name, baseline_config, baseline_override in baseline_specs:
            try:
                baseline_path, _, _ = create_loss_case_namelist(
                    case_name,
                    worker_dir / "baselines" / baseline_name,
                    selected_fields,
                    num_time_windows=num_time_windows,
                    config=baseline_config,
                    case_defaults=case_defaults or None,
                    batch_size=1,
                    duplicate_params_for_batch=True,
                    disable_stats_storage=True,
                    override=baseline_override,
                )
                os.chdir(RUN_SCRIPTS)
                baseline_fields, baseline_params = clubb_api.init_clubb_loss(
                    str(baseline_path),
                    return_default_params=True,
                )
                try:
                    baseline_metrics = clubb_api.clubb_get_loss_for_params(
                        np.asarray(baseline_params, dtype=np.float64)
                    )
                    scaled_rmse = np.asarray(baseline_metrics[0], dtype=float)
                    baseline_available = not (
                        np.any(~np.isfinite(scaled_rmse))
                        or np.any(scaled_rmse >= INVALID_SCALED_RMSE_PENALTY)
                    )
                finally:
                    clubb_api.finalize_clubb_loss()
            except Exception as exc:
                raise RuntimeError(
                    f"{case_name} {baseline_name} baseline evaluation failed: {exc}"
                ) from exc
            baselines[baseline_name] = {
                "available": baseline_available,
                "unavailable_reason": (
                    "CLUBB rejected this baseline parameter set"
                    if not baseline_available
                    else ""
                ),
                "field_names": list(baseline_fields),
                "loss_metric_names": list(LOSS_METRIC_NAMES),
                "loss_metrics_by_metric_window_field_and_column": {
                    metric_name: np.asarray(values, dtype=float).tolist()
                    for metric_name, values in zip(LOSS_METRIC_NAMES, baseline_metrics)
                },
            }
        conn.send({"baselines": baselines})
    except Exception as exc:
        conn.send({"error_message": _actionable_error_message(exc)})
    finally:
        os.chdir(original_cwd)
        conn.close()


def _evaluate_baselines_in_subprocess(init_payload: dict, control_conn) -> dict | None:
    """Keep one-column baseline allocations out of the candidate session."""
    context = multiprocessing.get_context("spawn")
    parent_conn, child_conn = context.Pipe()
    process = context.Process(target=_baseline_worker_main, args=(child_conn, init_payload))
    process.start()
    child_conn.close()
    message = None
    while message is None:
        if parent_conn.poll(0.1):
            try:
                message = parent_conn.recv()
            except EOFError as exc:
                process.join()
                raise RuntimeError(
                    f"Baseline worker for {init_payload.get('case_name')} exited with code {process.exitcode}"
                ) from exc
            break
        if control_conn.poll() and (control_conn.recv() or {}).get("type") == "stop":
            process.terminate()
            process.join()
            parent_conn.close()
            return None
        if not process.is_alive():
            process.join()
            parent_conn.close()
            raise RuntimeError(
                f"Baseline worker for {init_payload.get('case_name')} exited with code {process.exitcode}"
            )
    parent_conn.close()
    process.join()
    if process.exitcode != 0:
        raise RuntimeError(
            f"Baseline worker for {init_payload.get('case_name')} exited with code {process.exitcode}"
        )
    if message.get("error_message"):
        raise RuntimeError(message["error_message"])
    return dict(message.get("baselines") or {})


def worker_main(conn, init_payload: dict) -> None:
    """Own one initialized CLUBB loss session and evaluate parameter batches."""
    case_name = init_payload.get("case_name")
    selected_fields = list(init_payload.get("selected_fields", []))
    batch_size = int(init_payload.get("batch_size", 1))
    num_time_windows = int(init_payload.get("num_time_windows", 1))
    case_defaults = dict(init_payload.get("case_defaults") or {})
    config = init_payload.get("config") or "default"
    override = init_payload.get("override") or None
    evaluate_baselines = bool(init_payload.get("evaluate_baselines"))
    worker_dir = Path(init_payload.get("worker_dir")).resolve()
    original_cwd = os.getcwd()
    initialized = False

    try:
        baselines = _evaluate_baselines_in_subprocess(init_payload, conn) if evaluate_baselines else {}
        if baselines is None:
            return
        from clubb_python import clubb_api

        worker_dir.mkdir(parents=True, exist_ok=True)
        aggregate_path, _, _ = create_loss_case_namelist(
            case_name,
            worker_dir,
            selected_fields,
            num_time_windows=num_time_windows,
            config=config,
            case_defaults=case_defaults or None,
            batch_size=batch_size,
            duplicate_params_for_batch=True,
            disable_stats_storage=True,
            override=override,
        )

        os.chdir(RUN_SCRIPTS)
        field_names, default_params = clubb_api.init_clubb_loss(
            str(aggregate_path),
            return_default_params=True,
        )
        initialized = True
        param_names = list(clubb_api.get_param_names())
        nparams = len(clubb_api.get_param_names())
        hard_parameter_bounds = list(clubb_api.get_parameter_hard_bounds(nparams))

        conn.send(
            {
                "type": "initialized",
                "case_name": case_name,
                "worker_dir": str(worker_dir),
                "field_names": list(field_names),
                "param_names": param_names,
                # The worker owns the compiled F2PY module that will score
                # candidates.  Return its Fortran-owned hard envelope so the
                # scheduler cannot propose values that this exact build will
                # reject during model setup.
                "hard_parameter_bounds": hard_parameter_bounds,
                "default_params": np.asarray(default_params, dtype=float).tolist(),
                "baselines": baselines,
            }
        )

        while True:
            message = conn.recv()
            message_type = message.get("type")

            if message_type == "stop":
                break

            if message_type != "evaluate_batch":
                raise RuntimeError(f"Unknown worker message type: {message_type}")

            params_batch = np.asarray(message["params_batch"], dtype=np.float64)
            (
                scaled_rmse,
                correlation,
                std_ratio,
                centered_rmse_norm,
                bias_norm,
            ) = clubb_api.clubb_get_loss_for_params(params_batch)
            metric_arrays = {
                "scaled_rmse": scaled_rmse,
                "correlation": correlation,
                "std_ratio": std_ratio,
                "centered_rmse_norm": centered_rmse_norm,
                "bias_norm": bias_norm,
            }
            conn.send(
                {
                    "type": "result",
                    "case_name": case_name,
                    "batch_id": int(message["batch_id"]),
                    "field_names": list(field_names),
                    "loss_metric_names": list(LOSS_METRIC_NAMES),
                    "loss_metrics_by_metric_window_field_and_column": {
                        metric_name: np.asarray(metric_arrays[metric_name], dtype=float).tolist()
                        for metric_name in LOSS_METRIC_NAMES
                    },
                }
            )

    except EOFError:
        pass
    except Exception as exc:
        try:
            conn.send(
                {
                    "type": "error",
                    "case_name": case_name,
                    "error_message": _actionable_error_message(exc),
                }
            )
        except Exception:
            pass
    finally:
        try:
            os.chdir(original_cwd)
        except Exception:
            pass
        if initialized:
            try:
                from clubb_python import clubb_api

                clubb_api.finalize_clubb_loss()
            except Exception:
                pass
        try:
            conn.close()
        except Exception:
            pass
