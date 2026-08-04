"""Flask route for static report bundle assets, registered once at app startup."""

from __future__ import annotations

from flask import abort, send_from_directory

from .catalog import REPORTS_ROOT, report_by_id, resolve_report_asset


def register_static_report_routes(app) -> None:
    """Serve validated static report files without exposing the docs tree."""
    endpoint = "serve_static_report_asset"
    if endpoint in app.server.view_functions:
        return

    @app.server.route("/static-reports/<report_id>/<path:asset_path>", endpoint=endpoint)
    def serve_static_report_asset(report_id: str, asset_path: str):
        # A bundle reaches the browser only after it has been atomically added
        # to index.json; leave staging directories and orphaned artifacts dark.
        if report_by_id(report_id) is None:
            abort(404)
        asset = resolve_report_asset(report_id, asset_path)
        if asset is None:
            abort(404)
        return send_from_directory(REPORTS_ROOT / report_id, asset_path, max_age=3600)
