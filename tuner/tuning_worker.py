"""Worker process loop for multi-case CLUBB tuning."""

from __future__ import annotations

import os
from pathlib import Path
import multiprocessing
import sys
import traceback

import numpy as np

from tuner.paths import REPO_ROOT, RUN_SCRIPTS
from tuner.taylor_metrics import INVALID_RMSE_PENALTY, LOSS_METRIC_NAMES
from utilities.clubb_loss_runner import clubb_loss_runner, validate_gradient_backend

if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))



def _actionable_error_message(exc: BaseException) -> str:
    """Return an error string with the common module workaround when relevant."""
    message = f"{exc}\n{traceback.format_exc(limit=10)}"
    lower = message.lower()
    if "libnvomp" in lower or "netcdf" in lower or "no module named" in lower or "importerror" in lower:
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
        baselines = {}
        baseline_specs = [('clubb_default', 'default', None)]
        if override:
            baseline_specs.append(('override_defaults', config, override))
        for baseline_name, baseline_config, baseline_override in baseline_specs:
            os.chdir(RUN_SCRIPTS)
            with clubb_loss_runner(
                init_payload.get('backend', 'fortran'), case_name, selected_fields,
                worker_dir / 'baselines' / baseline_name, batch_size=1,
                case_defaults=dict(case_defaults, num_time_windows=num_time_windows),
                config=baseline_config, override=baseline_override,
                loss_weights=init_payload.get("loss_weights"), field_weights=init_payload.get("field_weights"),
            ) as runner:
                result = runner.evaluate(runner.default_params)
                rmse = result['metrics']['rmse']
                available = bool(np.isfinite(rmse).all()
                                 and np.all(rmse < INVALID_RMSE_PENALTY))
                baselines[baseline_name] = {
                    'available': available,
                    'unavailable_reason': '' if available else 'CLUBB rejected this baseline parameter set',
                    'case_loss': result['loss'].tolist(),
                    'field_loss': result['field_loss'].tolist(),
                    'field_names': list(runner.field_names),
                    'loss_metric_names': list(LOSS_METRIC_NAMES),
                    'loss_metrics_by_metric_window_field_and_column': {
                        name: values.tolist() for name, values in result['metrics'].items()},
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
        backend = init_payload.get("backend", "fortran")
        validate_gradient_backend(backend, init_payload.get('grad', False))
        os.chdir(RUN_SCRIPTS)
        runner = clubb_loss_runner(
            backend, case_name, selected_fields, worker_dir, batch_size=batch_size,
            case_defaults=dict(case_defaults, num_time_windows=num_time_windows),
            config=config, override=override,
            loss_weights=init_payload.get("loss_weights"), field_weights=init_payload.get("field_weights")).init()
        initialized = True
        field_names, default_params = runner.field_names, runner.default_params
        param_names, hard_parameter_bounds = runner.param_names, runner.bounds

        conn.send(
            {
                "type": "initialized",
                "case_name": case_name,
                "worker_dir": str(worker_dir),
                "field_names": list(field_names),
                "param_names": param_names,
                # The worker owns the runtime that will score candidates.
                # Return its source-matched hard envelope so the
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
            result = runner.evaluate(params_batch, grad=init_payload.get('grad', False),
                                     gradient_parameters=init_payload.get('gradient_parameters'),
                                     gradient_mode=init_payload.get('gradient_mode', 'auto'))
            metric_arrays = result['metrics']
            conn.send(
                {
                    "case_loss": result["loss"].tolist(),
                    "field_loss": result["field_loss"].tolist(),
                    "gradient_target": result.get("gradient_target"),
                    "gradient_parameters": result.get("gradient_parameters", []),
                    "gradient_mode": result.get('gradient_mode'),
                    "gradient": {name: values.tolist() for name, values in result.get("gradient", {}).items()},
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
                runner.close()
            except Exception:
                pass
        try:
            conn.close()
        except Exception:
            pass
