"""Serve completed PyPlotGen galleries from the fixed export root."""

from __future__ import annotations

import re

from flask import abort, send_from_directory

from .pyplotgen_runtime import PYPLOTGEN_OUTPUT_ROOT


_FOLDER_RE = re.compile(r"^[0-9]{8}_[0-9]{6}_[a-f0-9]{8}$")


def register_pyplotgen_routes(app) -> None:
    endpoint = "serve_clubb_pyplotgen_asset"
    if endpoint in app.server.view_functions:
        return

    @app.server.get("/_clubb-pyplotgen/<folder>/<path:asset_path>", endpoint=endpoint)
    def serve_pyplotgen_asset(folder: str, asset_path: str):
        if not _FOLDER_RE.fullmatch(folder):
            abort(404)
        export_dir = PYPLOTGEN_OUTPUT_ROOT / folder
        if not (export_dir / "index.html").is_file():
            abort(404)
        return send_from_directory(export_dir, asset_path, max_age=3600)
