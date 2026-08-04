"""Chunked numeric sample-history output for tuner analysis workflows."""

from __future__ import annotations

from pathlib import Path

import numpy as np


SCHEMA_VERSION = "sample_history_observation_v1"
DEFAULT_CHUNK_SIZE = 1000
SAMPLE_HISTORY_DIR = "sample_history"
CHUNK_PREFIX = "sample_history_"


class SampleHistoryWriter:
    """Write all completed tuner samples as compact NumPy chunks.

    The stored data is intentionally tuner-policy agnostic.  One row represents
    one evaluated parameter set.  Loss diagnostics are indexed by an
    observation axis, where each observation is a specific case/field/window.
    """

    def __init__(
        self,
        *,
        job_dir: Path,
        param_names: list[str],
        case_names: list[str],
        field_names: list[str],
        metric_names: list[str],
        case_configs: list[dict],
        case_window_counts: dict[str, int],
        chunk_size: int = DEFAULT_CHUNK_SIZE,
    ):
        self.job_dir = Path(job_dir)
        self.history_dir = self.job_dir / SAMPLE_HISTORY_DIR
        self.param_names = [str(name) for name in param_names]
        self.case_names = [str(name) for name in case_names]
        self.field_names = [str(name) for name in field_names]
        self.metric_names = [str(name) for name in metric_names]
        self.chunk_size = int(chunk_size)
        if self.chunk_size < 1:
            raise ValueError("sample history chunk_size must be >= 1")

        self.case_config_by_name = {
            str(config.get("case_name")): dict(config)
            for config in case_configs
        }
        self.case_window_counts = {
            str(case_name): int(count)
            for case_name, count in case_window_counts.items()
        }
        self.observations = self._build_observations()
        self.metadata = self._build_metadata()
        self.buffer: list[dict] = []
        # A stopped Tune revision may later continue in the same immutable
        # execution directory.  Continue chunk numbering rather than
        # overwriting the earlier raw observations.
        existing_chunks = sample_history_paths(self.job_dir)
        self.chunk_index = len(existing_chunks) + 1

    def append(self, entries: list[dict]) -> None:
        """Buffer completed sample entries and flush full chunks."""
        for entry in entries:
            self.buffer.append(entry)
            if len(self.buffer) >= self.chunk_size:
                self._flush_entries(self.buffer[: self.chunk_size])
                del self.buffer[: self.chunk_size]

    def close(self) -> None:
        """Flush any buffered samples."""
        if self.buffer:
            self._flush_entries(self.buffer)
            self.buffer = []

    def _build_observations(self) -> list[dict]:
        observations = []
        for case_idx, case_name in enumerate(self.case_names):
            window_count = int(self.case_window_counts.get(case_name, 1))
            config = self.case_config_by_name.get(case_name, {})
            time_range = list(config.get("time_average_range", []))
            if len(time_range) == 2 and window_count > 0:
                start = float(time_range[0])
                end = float(time_range[1])
                window_width = (end - start) / float(window_count)
            else:
                start = np.nan
                window_width = np.nan

            for field_idx, field_name in enumerate(self.field_names):
                for window_idx in range(window_count):
                    if np.isfinite(start) and np.isfinite(window_width):
                        window_start = start + window_idx * window_width
                        window_end = window_start + window_width
                    else:
                        window_start = np.nan
                        window_end = np.nan
                    observations.append(
                        {
                            "case_index": case_idx,
                            "field_index": field_idx,
                            "window_index": window_idx,
                            "case_name": case_name,
                            "field_name": field_name,
                            "window_start_seconds": window_start,
                            "window_end_seconds": window_end,
                        }
                    )
        return observations

    def _build_metadata(self) -> dict[str, np.ndarray]:
        case_window_counts = [
            int(self.case_window_counts.get(case_name, 1))
            for case_name in self.case_names
        ]
        return {
            "schema_version": np.asarray([SCHEMA_VERSION], dtype="U64"),
            "param_names": np.asarray(self.param_names, dtype="U96"),
            "case_names": np.asarray(self.case_names, dtype="U96"),
            "field_names": np.asarray(self.field_names, dtype="U96"),
            "metric_names": np.asarray(self.metric_names, dtype="U64"),
            "case_window_counts": np.asarray(case_window_counts, dtype=np.int32),
            "obs_case": np.asarray(
                [obs["case_index"] for obs in self.observations],
                dtype=np.int32,
            ),
            "obs_field": np.asarray(
                [obs["field_index"] for obs in self.observations],
                dtype=np.int32,
            ),
            "obs_window": np.asarray(
                [obs["window_index"] for obs in self.observations],
                dtype=np.int32,
            ),
            "obs_window_start_seconds": np.asarray(
                [obs["window_start_seconds"] for obs in self.observations],
                dtype=np.float64,
            ),
            "obs_window_end_seconds": np.asarray(
                [obs["window_end_seconds"] for obs in self.observations],
                dtype=np.float64,
            ),
        }

    def _flush_entries(self, entries: list[dict]) -> None:
        arrays = self._entries_to_arrays(entries)
        arrays.update(self.metadata)

        self.history_dir.mkdir(parents=True, exist_ok=True)
        path = self.history_dir / f"{CHUNK_PREFIX}{self.chunk_index:06d}.npz"
        tmp_path = path.with_suffix(path.suffix + ".tmp")
        with tmp_path.open("wb") as handle:
            np.savez_compressed(handle, **arrays)
        tmp_path.replace(path)
        self.chunk_index += 1

    def _entries_to_arrays(self, entries: list[dict]) -> dict[str, np.ndarray]:
        sample_count = len(entries)
        all_params = np.empty((sample_count, len(self.param_names)), dtype=np.float64)
        sample_id = np.empty(sample_count, dtype=np.int64)
        batch_id = np.empty(sample_count, dtype=np.int64)
        loss_metrics = np.empty(
            (sample_count, len(self.observations), len(self.metric_names)),
            dtype=np.float64,
        )

        for sample_idx, entry in enumerate(entries):
            all_param_values = entry.get("all_params", {})
            all_params[sample_idx, :] = [
                float(all_param_values[name])
                for name in self.param_names
            ]
            sample_id[sample_idx] = int(entry["sample_id"])
            batch_id[sample_idx] = int(entry["batch_id"])

            for obs_idx, obs in enumerate(self.observations):
                subwindows = (
                    entry["field_metrics"][obs["case_name"]][obs["field_name"]]["subwindows"]
                )
                subwindow = subwindows[int(obs["window_index"])]
                loss_metrics[sample_idx, obs_idx, :] = [
                    float(subwindow[metric_name])
                    for metric_name in self.metric_names
                ]

        return {
            "all_params": all_params,
            "sample_id": sample_id,
            "batch_id": batch_id,
            "loss_metrics": loss_metrics,
        }


def sample_history_paths(job_dir: Path) -> list[Path]:
    """Return sorted sample-history chunk paths in a job directory.

    New tuner runs store chunks under sample_history/.  The job-root lookup is
    kept as a compatibility fallback for runs created before that directory
    existed.
    """
    root = Path(job_dir)
    pattern = f"{CHUNK_PREFIX}[0-9][0-9][0-9][0-9][0-9][0-9].npz"
    chunks = sorted((root / SAMPLE_HISTORY_DIR).glob(pattern))
    if chunks:
        return chunks
    return sorted(root.glob(pattern))
