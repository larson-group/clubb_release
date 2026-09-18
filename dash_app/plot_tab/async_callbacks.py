"""Register native Plot callbacks as fixed broker tasks, with browser-side delivery."""
import json

from dash import Output, State
from flask import jsonify, request

from .tasks import ALLOWED_TASKS, _HANDLERS


def task_callback(app, name, *dependencies, **options):
    if name not in ALLOWED_TASKS:
        raise ValueError(f"Unregistered Plot task: {name}")

    def register(function):
        _HANDLERS[name] = function
        if getattr(app, "worker_registry", False):
            return function
        outputs = [item for item in dependencies if isinstance(item, Output)]
        inputs = [item for item in dependencies if not isinstance(item, Output)]
        specs = [{"id": out.component_id, "property": out.component_property} for out in outputs]
        # MATCH output ids must be resolved using the mounted card's actual id.
        graph = next((out for out in outputs if out.component_property == "figure"), None)
        if graph:
            inputs.append(State(graph.component_id, "id"))
            specs = [{"id": {**s["id"], "index": None}, "property": s["property"]} for s in specs]
        script = """function() {
            var args = Array.from(arguments), graph = GRAPH ? args.pop() : null;
            return window.clubbPlotTasks.run(NAME, args, OUTPUTS, graph,
                                             dash_clientside.callback_context);
        }""".replace("GRAPH", json.dumps(bool(graph))).replace("NAME", json.dumps(name)).replace("OUTPUTS", json.dumps(specs))
        app.clientside_callback(script, *outputs, *inputs,
                                prevent_initial_call=options.get("prevent_initial_call", False))
        return function
    return register


def register_task_routes(app):
    """Same-origin native UI proxy; broker credentials never enter the browser."""
    from dash_app.shared.broker_client import perform_action

    @app.server.post("/plots/tasks")
    def plot_task():
        if request.headers.get("Origin", request.host_url.rstrip("/")) != request.host_url.rstrip("/"):
            return jsonify(error="Cross-origin Plot requests are not accepted"), 403
        if request.content_length and request.content_length > 8 * 1024 * 1024:
            return jsonify(error="Plot request is too large"), 413
        data = request.get_json(silent=True) or {}
        if not isinstance(data, dict):
            return jsonify(error="Expected a Plot request object"), 400
        action = data.pop("operation", None)
        if action not in {"submit", "poll"}:
            return jsonify(error="Unknown Plot operation"), 400
        try:
            return jsonify(perform_action("plot_task_" + action, data, timeout_seconds=3,
                                          ensure_running=False))
        except (RuntimeError, ValueError, TypeError) as exc:
            return jsonify(error=str(exc)), 503
