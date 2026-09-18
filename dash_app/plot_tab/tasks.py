"""Broker-owned, bounded queues for transient native Plot requests.

Only registered Plot callbacks can run here. Workers are warm spawned processes,
so NetCDF stays out of HTTP threads and unchanged extraction caches survive.
Each browser/card has one newest request; superseded queued work is discarded.
These replaceable UI requests deliberately do not create scientific JobStore or
artifact records. A broker restart returns "missing" and the browser resubmits.
"""
from collections import OrderedDict
from concurrent.futures import ProcessPoolExecutor
import atexit
import json
import multiprocessing
import threading
import time

from dash._callback_context import context_value
from dash._utils import AttributeDict, to_json

ALLOWED_TASKS = frozenset({"catalog", "case", "params", "profile", "budget", "custom",
                           "pdf_contour", "subcolumn", "timeheight", "timeseries"})
_HANDLERS = {}


def execute(request):
    """Execute a fixed native callback with its original Dash context."""
    if not _HANDLERS:
        from .callbacks_case import register_case_callbacks
        from .callbacks_params import register_param_callbacks
        from .plot_types.registry import register_plot_callbacks

        class Registry:
            worker_registry = True
            def callback(self, *a, **kw):
                return lambda fn: fn
            def clientside_callback(self, *a, **kw):
                pass

        registry = Registry()
        register_case_callbacks(registry)
        register_param_callbacks(registry)
        register_plot_callbacks(registry)
    task = request["task"]
    if task not in ALLOWED_TASKS or task not in _HANDLERS:
        raise ValueError("Unknown Plot task")
    # Render state belongs to the browser's accepted figure, never another job.
    if task not in {"case", "catalog", "params"}:
        from .plot_types.registry import PLOT_TYPES
        plot = PLOT_TYPES[task]
        clear = getattr(plot, "clear_render_state", None)
        if clear:
            clear()
        if request.get("can_patch") and hasattr(plot, "_mark_full_render"):
            plot._mark_full_render(request["plot_index"])
    token = context_value.set(AttributeDict(triggered_inputs=request.get("triggered") or []))
    try:
        return to_json(_HANDLERS[task](*request["args"]))
    finally:
        context_value.reset(token)


class PlotTasks:
    """Two independent lanes: metadata and figures; one running job per lane."""
    def __init__(self, executor_factory=None, max_scopes=128, max_bytes=64 * 1024 * 1024, ttl=120):
        self.lock = threading.RLock()
        self.records = OrderedDict()
        self.pending = OrderedDict()
        self.running = {}
        self.pools = {}
        self.closed = False
        self.max_scopes, self.max_bytes, self.ttl = max_scopes, max_bytes, ttl
        self.executor_factory = executor_factory or (lambda: ProcessPoolExecutor(
            max_workers=1, mp_context=multiprocessing.get_context("spawn")))

    @staticmethod
    def lane(request):
        return "metadata" if request["task"] in {"catalog", "case", "params"} else "figures"

    def _prune(self):
        now = time.monotonic()
        for scope, record in list(self.records.items()):
            if now - record["touched"] > self.ttl:
                self.records.pop(scope)
                self.pending.pop(scope, None)
        size = sum(r.get("result_bytes", 0) for r in self.records.values())
        while len(self.records) > self.max_scopes or size > self.max_bytes:
            scope, record = self.records.popitem(last=False)
            self.pending.pop(scope, None)
            size -= record.get("result_bytes", 0)

    def submit(self, request):
        scope, revision = request.get("scope"), request.get("revision")
        if (not isinstance(scope, str) or not 1 <= len(scope) <= 180
                or not isinstance(revision, int) or revision < 1
                or request.get("task") not in ALLOWED_TASKS
                or not isinstance(request.get("args"), list)):
            raise ValueError("Invalid Plot request")
        with self.lock:
            if self.closed:
                raise RuntimeError("Plot workers are stopping")
            self._prune()
            page = scope.split(":", 1)[0]
            selection = int(request.get("selection", 0))
            if request["task"] != "catalog":
                related = [(key, value) for key, value in self.records.items()
                           if key.split(":", 1)[0] == page and value.get("task") != "catalog"]
                if any(value.get("selection", 0) > selection for _, value in related):
                    return {"state": "superseded", "revision": revision}
                for key, value in related:
                    if value.get("selection", 0) < selection:
                        self.records.pop(key)
                        self.pending.pop(key, None)
            previous = self.records.get(scope)
            if previous and previous["revision"] >= revision:
                return {"state": previous["state"], "revision": previous["revision"]}
            self.records[scope] = {"revision": revision, "task": request["task"], "selection": selection,
                                   "state": "queued", "touched": time.monotonic()}
            self.records.move_to_end(scope)
            self.pending[scope] = request
            self._prune()
            self._schedule()
            return {"state": "queued", "revision": revision}

    def _schedule(self):
        if self.closed:
            return
        for scope, request in list(self.pending.items()):
            if scope not in self.pending:
                continue
            lane = self.lane(request)
            if lane in self.running:
                continue
            self.pending.pop(scope)
            if lane not in self.pools:
                self.pools[lane] = self.executor_factory()
            try:
                future = self.pools[lane].submit(execute, request)
            except RuntimeError as exc:
                self.records[scope].update(state="error", error=str(exc))
                # A crashed process pool is replaceable on the next request.
                self.pools.pop(lane).shutdown(wait=False, cancel_futures=True)
                continue
            self.running[lane] = future
            record = self.records[scope]
            record["state"] = "running"
            future.add_done_callback(lambda f, s=scope, r=record, l=lane: self._done(s, r, l, f))

    def _done(self, scope, record, lane, future):
        with self.lock:
            self.running.pop(lane, None)
            if self.records.get(scope) is record:
                try:
                    result = future.result()
                    result_bytes = len(result.encode("utf-8"))
                    if result_bytes > self.max_bytes:
                        raise ValueError("Plot result exceeds the result-cache limit; reduce the selection")
                    record.update(state="ready", result=result, result_bytes=result_bytes)
                except Exception as exc:
                    record.update(state="error", error=str(exc))
                record["touched"] = time.monotonic()
            self._prune()
            self._schedule()

    def poll(self, scope, revision):
        with self.lock:
            self._prune()
            record = self.records.get(scope)
            if record is None:
                return {"state": "missing", "revision": revision}
            if record["revision"] != revision:
                return {"state": "superseded", "revision": record["revision"]}
            record["touched"] = time.monotonic()
            self.records.move_to_end(scope)
            response = {k: v for k, v in record.items() if k not in {"touched", "result", "result_bytes"}}
            if "result" in record:
                response["result"] = json.loads(record["result"])
            return response

    def poll_many(self, requests):
        """One native UI round trip for all cards currently awaiting results."""
        if (not isinstance(requests, list) or len(requests) > self.max_scopes
                or any(not isinstance(item, dict) for item in requests)):
            raise ValueError("Invalid Plot polling batch")
        with self.lock:
            return {"results": [self.poll(str(item.get("scope") or ""), int(item.get("revision") or 0))
                                for item in requests]}

    def close(self):
        with self.lock:
            self.closed = True
            self.pending.clear()
            pools = list(self.pools.values())
            self.pools.clear()
        for pool in pools:
            # Python 3.14 supplies immediate worker termination. Older supported
            # versions finish only the bounded in-flight tasks during shutdown.
            terminate = getattr(pool, "terminate_workers", None)
            if terminate:
                terminate()
            else:
                pool.shutdown(wait=True, cancel_futures=True)


_MANAGER = None
_MANAGER_LOCK = threading.Lock()


def manager():
    global _MANAGER
    with _MANAGER_LOCK:
        if _MANAGER is None:
            _MANAGER = PlotTasks()
            atexit.register(_MANAGER.close)
        return _MANAGER


def close_workers():
    if _MANAGER is not None:
        _MANAGER.close()
