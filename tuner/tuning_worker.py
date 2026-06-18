"""Worker process loop for multi-case CLUBB tuning."""

from __future__ import annotations

import os
from pathlib import Path
import sys
import traceback

import numpy as np

from tuner.paths import RUN_SCRIPTS
from tuner.taylor_metrics import LOSS_METRIC_NAMES

if str(RUN_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(RUN_SCRIPTS))

from create_case_namelist import create_loss_case_namelist  # noqa: E402


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


def worker_main(conn, init_payload: dict) -> None:
    """Own one initialized CLUBB loss session and evaluate parameter batches."""
    case_name = init_payload.get("case_name")
    selected_fields = list(init_payload.get("selected_fields", []))
    batch_size = int(init_payload.get("batch_size", 1))
    num_time_windows = int(init_payload.get("num_time_windows", 1))
    case_defaults = dict(init_payload.get("case_defaults") or {})
    worker_dir = Path(init_payload.get("worker_dir")).resolve()
    original_cwd = os.getcwd()
    initialized = False

    try:
        from clubb_python import clubb_api

        worker_dir.mkdir(parents=True, exist_ok=True)
        aggregate_path, _, _ = create_loss_case_namelist(
            case_name,
            worker_dir,
            selected_fields,
            num_time_windows=num_time_windows,
            case_defaults=case_defaults or None,
            batch_size=batch_size,
            duplicate_params_for_batch=True,
            disable_stats_storage=True,
        )

        os.chdir(RUN_SCRIPTS)
        field_names, default_params = clubb_api.init_clubb_loss(
            str(aggregate_path),
            return_default_params=True,
        )
        initialized = True
        param_names = list(clubb_api.get_param_names())

        conn.send(
            {
                "type": "initialized",
                "case_name": case_name,
                "worker_dir": str(worker_dir),
                "field_names": list(field_names),
                "param_names": param_names,
                "default_params": np.asarray(default_params, dtype=float).tolist(),
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
