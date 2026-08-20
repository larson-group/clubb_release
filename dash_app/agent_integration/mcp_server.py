#!/usr/bin/env python3
"""MCP adapter for the agent-neutral CLUBB Dash local gateway."""

from __future__ import annotations

import asyncio
import argparse
from contextlib import redirect_stdout
import hmac
import json
import os
import sys
import threading
import time
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from dash_app.shared import actions
from dash_app.agent_integration.client import connect, perform_action
from dash_app.agent_integration.broker_endpoint import (
    cleanup_stable_record,
    process_started_at as stable_process_started_at,
    stable_endpoint_record_is_live,
)
from dash_app.agent_integration.endpoint import (
    ENDPOINT_PATH,
    cleanup_instance_record,
    endpoint_record_is_live,
    load_endpoint_config,
    process_started_at,
)
from dash_app.services import CompileRequest, LeaderboardRerunRequest, ProfileArtifactRequest, ScmRunBatchRequest, ScmRunRequest, TuneRequest


def _domain_result(call, *args):
    """Keep expected service failures structured and free of traceback text."""
    try:
        return call(*args)
    except ValueError as exc:
        message = str(exc)[:500]
        code = "REQUEST_ID_CONFLICT" if message.startswith("REQUEST_ID_CONFLICT:") else "INVALID_REQUEST"
        return {"error": {"code": code, "message": message.removeprefix("REQUEST_ID_CONFLICT: ")}}
    except OSError:
        return {"error": {"code": "SERVICE_IO_ERROR", "message": "local service data is unavailable"}}
    except Exception:
        return {"error": {"code": "SERVICE_FAILURE", "message": "local CLUBB service operation failed"}}


def _gateway_result(call, *args, **kwargs):
    """Keep a transient Dash/broker outage from aborting the MCP tool call.

    FastMCP wraps an exception raised by a tool as ``ToolError``.  The custom
    stdio transport cannot distinguish that execution failure from argument
    validation, and some MCP hosts close the tool transport after receiving
    it.  Browser-handoff tools therefore return the same bounded structured
    error style as domain tools while leaving the MCP session usable.
    """

    try:
        return call(*args, **kwargs)
    except RuntimeError as exc:
        return {
            "error": {
                "code": "DASHBOARD_CONNECTION_FAILED",
                "message": str(exc)[:500],
            }
        }
    except (OSError, ValueError):
        return {
            "error": {
                "code": "DASHBOARD_REQUEST_FAILED",
                "message": "the local dashboard request could not be completed",
            }
        }


def _broker_domain_result(
    operation: str,
    payload: dict[str, Any],
    *,
    connection: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Execute one typed mutation in the detached broker process.

    The stdio adapter is owned by the agent host and may be restarted or
    closed at any time.  Watchers for Compile, SCM, Tune, and artifact jobs
    must therefore be created by the durable broker, not by this process.
    """

    try:
        kwargs: dict[str, Any] = {
            "internal": True,
            "timeout_seconds": 120.0,
        }
        if connection is not None:
            kwargs["connection"] = connection
        return perform_action(f"domain_{operation}", payload, **kwargs)
    except RuntimeError as exc:
        return {
            "error": {
                "code": "BROKER_OPERATION_FAILED",
                "message": str(exc)[:500],
            }
        }


def create_server(
    *,
    broker_connection: dict[str, Any] | None = None,
    instance_info: dict[str, Any] | None = None,
):
    try:
        from mcp.server.fastmcp import FastMCP
        from mcp.types import ToolAnnotations
    except ModuleNotFoundError as exc:
        raise RuntimeError("Install dash_app/requirements.txt before starting the CLUBB Dash MCP adapter.") from exc

    server_kwargs: dict[str, Any] = {}
    if broker_connection is not None:
        server_kwargs.update({"json_response": True, "stateless_http": True, "host": "127.0.0.1"})
    server = FastMCP("CLUBB Dash", **server_kwargs)
    read_only = ToolAnnotations(readOnlyHint=True, destructiveHint=False, idempotentHint=True, openWorldHint=False)
    # Domain submissions below are idempotent only because their strict
    # request_id is persisted by JobStore.  Browser handoffs remain transient
    # requests and do not create an agent session.
    mutate = ToolAnnotations(readOnlyHint=False, destructiveHint=False, idempotentHint=False, openWorldHint=False)
    idempotent_mutate = ToolAnnotations(readOnlyHint=False, destructiveHint=False, idempotentHint=True, openWorldHint=False)
    cancel = ToolAnnotations(readOnlyHint=False, destructiveHint=True, idempotentHint=True, openWorldHint=False)

    def gateway_connect() -> dict[str, Any]:
        result = connect(connection=broker_connection)
        if instance_info:
            result["dashboard_instance"] = {
                "instance_id": instance_info.get("instance_id"),
                "dashboard_pid": instance_info.get("dashboard", {}).get("pid"),
                "broker_pid": instance_info.get("broker", {}).get("pid"),
            }
        return result

    def gateway_action(action: str, payload: dict[str, Any]) -> dict[str, Any]:
        return perform_action(action, payload, connection=broker_connection)

    def broker_domain(operation: str, payload: dict[str, Any]) -> dict[str, Any]:
        request_payload = dict(payload)
        if operation in {"submit_scm_run", "submit_scm_batch"}:
            request_payload["submission_origin"] = "mcp"
        return _broker_domain_result(operation, request_payload, connection=broker_connection)

    @server.tool(annotations=read_only)
    def connect_to_dashboard() -> dict[str, Any]:
        """Verify one short-lived authenticated connection to the local broker."""
        return _gateway_result(gateway_connect)

    @server.tool(annotations=read_only)
    def inspect_dashboard(tab: str | None = None) -> dict[str, Any]:
        """List agent-ready dashboard tabs and their safe, typed semantic operations."""
        return _gateway_result(
            gateway_action,
            "inspect_dashboard",
            {"tab": tab} if tab else {},
        )

    @server.tool(annotations=mutate)
    def invoke_dashboard(tab: str, operation: str, arguments: dict[str, Any] | None = None) -> dict[str, Any]:
        """Invoke one operation returned by inspect_dashboard; no shell or browser automation is accepted."""
        return _gateway_result(
            gateway_action,
            "invoke_dashboard",
            {"tab": tab, "operation": operation, "arguments": arguments or {}},
        )

    @server.tool(annotations=mutate)
    def add_budget_plot(
        case: str,
        budget_group: str = "wp2",
        run_id: str | None = None,
        output_dir: str | None = None,
        time_start_seconds: float | None = None,
        average_minutes: float | None = None,
        window_preset: str | None = None,
    ) -> dict[str, Any]:
        """Add a validated budget plot to the active Plot tab; WP2 is the default group."""
        arguments = {
            "case": case,
            "budget_group": budget_group,
            "run_id": run_id,
            "output_dir": output_dir,
            "time_start_seconds": time_start_seconds,
            "average_minutes": average_minutes,
            "window_preset": window_preset,
        }
        return _gateway_result(
            gateway_action,
            "invoke_dashboard",
            {"tab": "plots", "operation": "add_budget", "arguments": arguments},
        )

    @server.tool(annotations=read_only)
    def list_plots() -> dict[str, Any]:
        """List active Plot cards and their stable IDs from the current dashboard."""
        return _gateway_result(
            gateway_action,
            "invoke_dashboard",
            {"tab": "plots", "operation": "list", "arguments": {}},
        )

    @server.tool(annotations=mutate)
    def remove_plot(plot_id: int) -> dict[str, Any]:
        """Remove one active Plot card by the stable ID returned by list_plots."""
        return _gateway_result(
            gateway_action,
            "invoke_dashboard",
            {"tab": "plots", "operation": "remove", "arguments": {"plot_id": plot_id}},
        )

    @server.tool(annotations=read_only)
    def get_server_info() -> dict[str, Any]:
        """Read the local-only CLUBB service capabilities and compatibility policy."""
        return _domain_result(actions.get_server_info)

    @server.tool(annotations=read_only)
    def list_cases() -> dict[str, Any]:
        """List supported SCM cases and checked-in stats choices."""
        return _domain_result(actions.list_cases)

    @server.tool(annotations=read_only)
    def list_tunable_parameters(config: str = "default") -> dict[str, Any]:
        """Read current Tune parameter names and suggested ranges from one selected config."""
        return _domain_result(actions.list_tunable_parameters, config)

    @server.tool(annotations=read_only)
    def resolve_clubb_settings(
        config: str = "default",
        flags: dict[str, Any] | None = None,
        parameters: dict[str, Any] | None = None,
        auto_correct: bool = False,
    ) -> dict[str, Any]:
        """Preview flag conflicts, linked parameters, and active controls without a model run."""
        return _domain_result(
            actions.resolve_clubb_settings_request,
            config,
            flags,
            parameters,
            auto_correct,
        )

    @server.tool(annotations=read_only)
    def list_tune_workspaces() -> dict[str, Any]:
        """List saved Tune workspaces, immutable revisions, status, timestamps, and disk use."""
        return _domain_result(actions.list_tuning_workspaces)

    @server.tool(annotations=read_only)
    def get_tune_workspace(workspace_id: str, revision_id: str) -> dict[str, Any]:
        """Read one saved Tune revision's immutable request, status, results, and paths."""
        return _domain_result(actions.load_tuning_workspace, workspace_id, revision_id)

    @server.tool(annotations=mutate)
    def create_tune_revision(workspace_id: str, revision_id: str) -> dict[str, Any]:
        """Clone a saved revision into a new editable, unstarted revN draft."""
        return broker_domain(
            "create_tune_revision",
            {"workspace_id": workspace_id, "revision_id": revision_id},
        )

    @server.tool(annotations=mutate)
    def restart_tune_revision(workspace_id: str, revision_id: str) -> dict[str, Any]:
        """Create and start a fresh exact <revision>_restartN attempt from a saved revision."""
        return broker_domain(
            "restart_tune_revision",
            {"workspace_id": workspace_id, "revision_id": revision_id},
        )

    @server.tool(annotations=cancel)
    def reset_tune_revision(workspace_id: str, revision_id: str) -> dict[str, Any]:
        """Delete an inactive revision's generated data and restore it as an editable draft."""
        return broker_domain(
            "reset_tune_revision",
            {"workspace_id": workspace_id, "revision_id": revision_id},
        )

    @server.tool(annotations=mutate)
    def start_tune_draft_revision(workspace_id: str, revision_id: str) -> dict[str, Any]:
        """Start an unchanged, unstarted saved draft revision through the validated Tune service."""
        return broker_domain(
            "start_tune_draft_revision",
            {"workspace_id": workspace_id, "revision_id": revision_id},
        )

    @server.tool(annotations=mutate)
    def resume_tune_revision(workspace_id: str, revision_id: str) -> dict[str, Any]:
        """Continue a gracefully stopped Tune revision from its durable checkpoint."""
        return broker_domain(
            "resume_tune_revision",
            {"workspace_id": workspace_id, "revision_id": revision_id},
        )

    @server.tool(annotations=mutate)
    def rename_tune_workspace(workspace_id: str, display_name: str) -> dict[str, Any]:
        """Rename a Tune workspace display label without moving its files."""
        return broker_domain(
            "rename_tune_workspace",
            {"workspace_id": workspace_id, "display_name": display_name},
        )

    @server.tool(annotations=cancel)
    def delete_tune_workspace(workspace_id: str) -> dict[str, Any]:
        """Delete an inactive Tune workspace lineage; active revisions are rejected."""
        return broker_domain(
            "delete_tune_workspace",
            {"workspace_id": workspace_id},
        )

    @server.tool(annotations=idempotent_mutate)
    def submit_compile(request: CompileRequest) -> dict[str, Any]:
        """Submit an idempotent CLUBB compile.  Reuse request_id on retry."""
        return broker_domain("submit_compile", {"request": request.model_dump(mode="json")})

    @server.tool(annotations=idempotent_mutate)
    def submit_scm_run(request: ScmRunRequest) -> dict[str, Any]:
        """Submit one reproducible SCM case through the shared batch service."""
        return broker_domain("submit_scm_run", {"request": request.model_dump(mode="json")})

    @server.tool(annotations=idempotent_mutate)
    def submit_scm_batch(request: ScmRunBatchRequest) -> dict[str, Any]:
        """Queue multiple SCM cases into one flat, discoverable output group."""
        return broker_domain("submit_scm_batch", {"request": request.model_dump(mode="json")})

    @server.tool(annotations=idempotent_mutate)
    def submit_tune(request: TuneRequest) -> dict[str, Any]:
        """Submit one bounded Tune job using the existing scientific validator."""
        return broker_domain("submit_tune", {"request": request.model_dump(mode="json")})

    @server.tool(annotations=idempotent_mutate)
    def submit_leaderboard_rerun(request: LeaderboardRerunRequest) -> dict[str, Any]:
        """Idempotently rerun current Tune leaderboard rows as a cancellable job."""
        return broker_domain(
            "submit_leaderboard_rerun",
            {"request": request.model_dump(mode="json")},
        )

    @server.tool(annotations=idempotent_mutate)
    def create_profile_artifact(request: ProfileArtifactRequest) -> dict[str, Any]:
        """Export native profile PNGs from one immutable SCM run selection."""
        return broker_domain(
            "create_profile_artifact",
            {"request": request.model_dump(mode="json")},
        )

    @server.tool(annotations=read_only)
    def get_job(job_id: str) -> dict[str, Any]:
        """Read a durable job record and its compact provenance."""
        return _domain_result(actions.get_job, job_id)

    @server.tool(annotations=read_only)
    def get_run_manifest(run_id: str) -> dict[str, Any]:
        """Read the immutable manifest for a submitted SCM run."""
        return _domain_result(actions.get_run_manifest, run_id)

    @server.tool(annotations=read_only)
    def get_artifact(artifact_id: str) -> dict[str, Any]:
        """Read one artifact manifest; binary payloads remain artifact resources."""
        return _domain_result(actions.get_artifact, artifact_id)

    @server.tool(annotations=read_only)
    def read_job_log(job_id: str, cursor: int = 0, max_bytes: int = 8192) -> dict[str, Any]:
        """Read a bounded incremental job log chunk."""
        return _domain_result(actions.read_job_log, job_id, cursor, max_bytes)

    @server.tool(annotations=cancel)
    def cancel_job(job_id: str) -> dict[str, Any]:
        """Request safe cancellation for a broker-owned compile, SCM, or Tune job."""
        return broker_domain("cancel_job", {"job_id": job_id})

    @server.resource("clubb-artifact://{artifact_id}/manifest.json", mime_type="application/json")
    def artifact_manifest_resource(artifact_id: str) -> str:
        """Expose compact manifests as MCP resources, never NetCDF payloads."""
        return __import__("json").dumps(actions.get_artifact(artifact_id)["manifest"], indent=2, sort_keys=True)

    @server.resource("clubb-artifact://{artifact_id}/{filename}")
    def artifact_file_resource(artifact_id: str, filename: str) -> bytes:
        """Read a whitelisted PNG/text artifact with a bounded transfer size."""
        return actions.read_artifact_file(artifact_id, filename)

    return server


def _content_json(content) -> list[dict[str, Any]]:
    return [item.model_dump(by_alias=True, exclude_none=True) for item in content]


def _rpc_result(request_id: Any, result: dict[str, Any]) -> None:
    sys.stdout.write(json.dumps({"jsonrpc": "2.0", "id": request_id, "result": result}, separators=(",", ":")) + "\n")
    sys.stdout.flush()


def _rpc_error(request_id: Any, code: int, message: str, *, data: dict[str, Any] | None = None) -> None:
    error: dict[str, Any] = {"code": code, "message": message}
    if data:
        error["data"] = data
    sys.stdout.write(json.dumps({"jsonrpc": "2.0", "id": request_id, "error": error}, separators=(",", ":")) + "\n")
    sys.stdout.flush()


def _serve_stdio(server) -> None:
    """Serve the standard line-delimited MCP protocol without async stdio.

    MCP 1.28's bundled async stdio transport hangs on CPython 3.14 in this
    environment (including for a one-tool upstream FastMCP probe).  This small
    synchronous transport implements standard discovery, tools, and resources
    while retaining FastMCP for schema generation and tool invocation.
    """
    from mcp.shared.version import SUPPORTED_PROTOCOL_VERSIONS
    from mcp.server.fastmcp.exceptions import ToolError
    from pydantic import ValidationError

    initialized = False
    for raw_line in sys.stdin:
        try:
            request = json.loads(raw_line)
        except json.JSONDecodeError:
            _rpc_error(None, -32700, "parse error")
            continue
        method = str(request.get("method") or "")
        request_id = request.get("id")
        params = dict(request.get("params") or {})
        if request_id is None:
            # A compliant client sends notifications/initialized after the
            # initialize response.  No action is required because this small
            # local server has no deferred capability setup.
            continue
        try:
            if method == "initialize":
                requested_version = str(params.get("protocolVersion") or "")
                if requested_version not in SUPPORTED_PROTOCOL_VERSIONS:
                    _rpc_error(
                        request_id,
                        -32602,
                        "unsupported MCP protocol version",
                        data={"code": "UNSUPPORTED_PROTOCOL", "supported": list(SUPPORTED_PROTOCOL_VERSIONS)},
                    )
                    continue
                initialized = True
                _rpc_result(request_id, {"protocolVersion": requested_version, "capabilities": {"tools": {}, "resources": {"subscribe": False, "listChanged": False}}, "serverInfo": {"name": "CLUBB Dash", "version": "3"}})
            elif not initialized:
                _rpc_error(request_id, -32002, "MCP session is not initialized", data={"code": "NOT_INITIALIZED"})
            elif method == "tools/list":
                tools = asyncio.run(server.list_tools())
                _rpc_result(request_id, {"tools": [tool.model_dump(by_alias=True, exclude_none=True) for tool in tools]})
            elif method == "tools/call":
                # Domain launchers predate MCP and may print a concise command
                # or subprocess notice.  Stdout is the JSON-RPC transport, so
                # even one such line makes an MCP host report ``Transport
                # closed``.  Preserve diagnostics on stderr and reserve
                # stdout exclusively for protocol frames.
                with redirect_stdout(sys.stderr):
                    content, structured = asyncio.run(
                        server.call_tool(
                            str(params.get("name") or ""),
                            dict(params.get("arguments") or {}),
                        )
                    )
                result: dict[str, Any] = {"content": _content_json(content)}
                if structured is not None:
                    result["structuredContent"] = structured
                _rpc_result(request_id, result)
            elif method == "resources/list":
                _rpc_result(request_id, {"resources": []})
            elif method == "resources/templates/list":
                templates = asyncio.run(server.list_resource_templates())
                _rpc_result(request_id, {"resourceTemplates": [template.model_dump(by_alias=True, exclude_none=True) for template in templates]})
            elif method == "resources/read":
                uri = str(params.get("uri") or "")
                prefix = "clubb-artifact://"
                if not uri.startswith(prefix):
                    raise ValueError("unknown resource URI")
                artifact_id, filename = uri[len(prefix):].split("/", 1)
                data = actions.read_artifact_file(artifact_id, filename)
                if filename.endswith(".png"):
                    import base64

                    contents = [{"uri": uri, "mimeType": "image/png", "blob": base64.b64encode(data).decode("ascii")}]
                else:
                    contents = [{"uri": uri, "mimeType": "application/json" if filename.endswith(".json") else "text/plain", "text": data.decode("utf-8", errors="replace")}]
                _rpc_result(request_id, {"contents": contents})
            elif method == "ping":
                _rpc_result(request_id, {})
            else:
                _rpc_error(request_id, -32601, "method not found")
        except ToolError:
            # FastMCP wraps its Pydantic argument-validation exception in a
            # ToolError before it reaches this transport.  Do not surface its
            # implementation-specific traceback/text to another MCP client.
            _rpc_error(
                request_id,
                -32602,
                "invalid tool arguments",
                data={"code": "INVALID_REQUEST", "message": "tool arguments do not match the published schema"},
            )
        except ValidationError as exc:
            details = [
                {
                    "path": ".".join(str(part) for part in item.get("loc") or ()),
                    "message": str(item.get("msg") or "invalid value"),
                    "type": str(item.get("type") or "validation_error"),
                }
                for item in exc.errors()[:16]
            ]
            _rpc_error(
                request_id,
                -32602,
                "invalid tool arguments",
                data={"code": "INVALID_REQUEST", "details": details},
            )
        except ValueError as exc:
            _rpc_error(
                request_id,
                -32602,
                "invalid request",
                data={"code": "INVALID_REQUEST", "message": str(exc)[:500]},
            )
        except Exception:
            _rpc_error(request_id, -32603, "internal local MCP service error")


class _AuthenticatedEndpoint:
    """Loopback bearer-token and owner-liveness boundary for MCP HTTP."""

    def __init__(self, app, config: dict[str, Any], record: dict[str, Any]):
        self.app = app
        self.config = config
        self.record = record

    @staticmethod
    async def _respond(send, status: int, payload: dict[str, Any]) -> None:
        body = json.dumps(payload, separators=(",", ":")).encode("utf-8")
        await send(
            {
                "type": "http.response.start",
                "status": int(status),
                "headers": [
                    [b"content-type", b"application/json"],
                    [b"content-length", str(len(body)).encode("ascii")],
                ],
            }
        )
        await send({"type": "http.response.body", "body": body})

    def _authorized(self, scope: dict[str, Any]) -> bool:
        headers = {
            key.decode("latin-1").lower(): value.decode("latin-1")
            for key, value in scope.get("headers", [])
        }
        supplied = str(headers.get("authorization") or "")
        expected = f"Bearer {self.config.get('endpoint_token')}"
        return bool(self.config.get("endpoint_token")) and hmac.compare_digest(supplied, expected)

    def _owner_is_live(self) -> bool:
        if self.record.get("owner"):
            return stable_endpoint_record_is_live(self.record)
        return endpoint_record_is_live(self.record)

    async def __call__(self, scope, receive, send):
        if scope.get("type") != "http":
            await self.app(scope, receive, send)
            return
        if not self._authorized(scope):
            await self._respond(
                send,
                401,
                {"error": {"code": "AUTH_REQUIRED", "message": "dashboard MCP bearer token required"}},
            )
            return
        if self.record.get("endpoint_token") != self.config.get("endpoint_token") or not self._owner_is_live():
            if self.record.get("owner"):
                cleanup_stable_record(expected_endpoint_pid=os.getpid())
            else:
                cleanup_instance_record(str(self.record["instance_id"]), expected_endpoint_pid=os.getpid())
            await self._respond(
                send,
                503,
                {
                    "error": {
                    "code": "BROKER_ENDPOINT_UNAVAILABLE",
                    "message": "the owning broker endpoint is no longer live",
                    }
                },
            )
            return
        if scope.get("path") == "/health":
            await self._respond(
                send,
                200,
                {
                    "status": "ok",
                    "instance_id": self.record["instance_id"],
                    "dashboard_pid": (self.record.get("dashboard") or {}).get("pid"),
                    "broker_pid": (self.record.get("broker") or {}).get("pid"),
                    "owner_pid": (self.record.get("owner") or {}).get("pid"),
                },
            )
            return
        await self.app(scope, receive, send)


def _serve_http(config_path: Path) -> None:
    config = load_endpoint_config(config_path)
    endpoint_started_at = stable_process_started_at(os.getpid()) if config.get("owner") else process_started_at(os.getpid())
    initial_record: dict[str, Any] = {
        "version": config["version"],
        "instance_id": config["instance_id"],
        "endpoint_pid": os.getpid(),
        "endpoint_started_at": endpoint_started_at,
        "broker": config["broker"],
    }
    if config.get("owner"):
        initial_record["owner"] = config["owner"]
        if not stable_endpoint_record_is_live(initial_record, check_broker=True):
            raise RuntimeError("the owning broker is not live")
    else:
        initial_record["dashboard"] = config["dashboard"]
        if not endpoint_record_is_live(initial_record, check_broker=True):
            raise RuntimeError("the owning dashboard or broker is not live")

    import uvicorn

    instance_info = {
        "instance_id": config["instance_id"],
        "dashboard": config.get("dashboard") or {},
        "broker": config["broker"],
    }
    server = create_server(broker_connection=dict(config["broker"]), instance_info=instance_info)
    mcp_app = server.streamable_http_app()
    port = int(config["endpoint_port"])
    record = {
        **initial_record,
        "endpoint_url": f"http://127.0.0.1:{port}{ENDPOINT_PATH}",
        "endpoint_token": config["endpoint_token"],
        "started_at": time.time(),
    }
    registry_path = Path(str(config["registry_path"]))
    registry_path.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
    from dash_app.shared.runtime import atomic_write_json

    atomic_write_json(registry_path, record)
    wrapped_app = _AuthenticatedEndpoint(mcp_app, config, record)

    def monitor_owner() -> None:
        while True:
            live = stable_endpoint_record_is_live(record) if record.get("owner") else endpoint_record_is_live(record)
            if not live:
                if record.get("owner"):
                    cleanup_stable_record(expected_endpoint_pid=os.getpid())
                else:
                    cleanup_instance_record(str(record["instance_id"]), expected_endpoint_pid=os.getpid())
                os._exit(0)
            time.sleep(1.0)

    threading.Thread(target=monitor_owner, name="dashboard-mcp-owner-monitor", daemon=True).start()
    try:
        uvicorn.run(wrapped_app, host="127.0.0.1", port=port, log_level="warning", access_log=False)
    finally:
        if record.get("owner"):
            cleanup_stable_record(expected_endpoint_pid=os.getpid())
        else:
            cleanup_instance_record(str(record["instance_id"]), expected_endpoint_pid=os.getpid())


def main() -> None:
    parser = argparse.ArgumentParser(description="CLUBB Dash MCP adapter")
    parser.add_argument("--http-config", type=Path, help="run one broker-owned authenticated Streamable HTTP endpoint")
    args = parser.parse_args()
    if args.http_config is not None:
        _serve_http(args.http_config)
        return
    _serve_stdio(create_server())


if __name__ == "__main__":
    main()
