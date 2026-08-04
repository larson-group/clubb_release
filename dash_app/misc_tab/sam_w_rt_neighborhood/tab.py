"""Instant 5-by-5 browser over a pre-rendered raw-SAM w-rt image atlas."""

from __future__ import annotations

from functools import lru_cache
import json
from pathlib import Path

from dash import ALL, Input, Output, State, ctx, dcc, html, no_update
from flask import send_from_directory
import numpy as np

from dash_app.shared.reporting import empty_state, report_header

from utilities.sam_3d_reference import DEFAULT_SAM_RUN

from .generation import job_status, start_atlas_generation


TIME_SLIDER_ID = "notes-sam-wrt-time"
HEIGHT_SLIDER_ID = "notes-sam-wrt-height"
TIME_LABEL_ID = "notes-sam-wrt-time-label"
HEIGHT_LABEL_ID = "notes-sam-wrt-height-label"
STATUS_ID = "notes-sam-wrt-status"
MANIFEST_STORE_ID = "notes-sam-wrt-atlas-manifest"
SPRITE_TYPE = "notes-sam-wrt-sprite"
SPRITE_LABEL_TYPE = "notes-sam-wrt-sprite-label"
GENERATE_ID = "notes-sam-wrt-generate"
GENERATION_STATUS_ID = "notes-sam-wrt-generation-status"
GENERATION_INTERVAL_ID = "notes-sam-wrt-generation-interval"
GENERATION_RELOAD_ID = "notes-sam-wrt-generation-reload"

DEFAULT_TIME_MINUTES = 294.0
DEFAULT_HEIGHT_M = 820.0
GRID_RADIUS = 2
ATLAS_DIR = Path(DEFAULT_SAM_RUN).expanduser() / "PDF_TILES" / "w_rt_signed"
MANIFEST_PATH = ATLAS_DIR / "manifest.json"
ATLAS_URL = "/sam-wrt-atlas"


def _nearest_index(values, target: float) -> int:
    return int(np.argmin(np.abs(np.asarray(values, dtype=float) - float(target))))


def _valid_center_bounds(size: int) -> tuple[int, int]:
    if size < 2 * GRID_RADIUS + 1:
        raise ValueError("The SAM atlas needs at least five time and height records.")
    return GRID_RADIUS, size - GRID_RADIUS - 1


def _slider_marks(values, low: int, high: int, formatter) -> dict[int, str]:
    indices = np.unique(np.rint(np.linspace(low, high, 7)).astype(int))
    return {int(index): formatter(values[int(index)]) for index in indices}


def _validated_manifest(payload: dict) -> dict:
    """Validate the small contract shared by the generator and browser."""

    required = {
        "version",
        "projection",
        "generated_utc",
        "time_seconds",
        "height_m",
        "tile_pixels",
        "time_block",
        "height_block",
        "atlas_filename",
        "x_range_rt_g_per_kg",
        "y_range_w_m_per_s",
    }
    missing = sorted(required.difference(payload))
    if missing:
        raise ValueError(f"SAM atlas manifest is missing: {', '.join(missing)}")
    if payload["projection"] != "w_rt":
        raise ValueError("The Misc explorer requires a w_rt atlas.")
    time_count = len(payload["time_seconds"])
    height_count = len(payload["height_m"])
    _valid_center_bounds(time_count)
    _valid_center_bounds(height_count)
    if len(payload["x_range_rt_g_per_kg"]) != height_count:
        raise ValueError("SAM atlas r_t ranges do not match its height coordinate.")
    if len(payload["y_range_w_m_per_s"]) != height_count:
        raise ValueError("SAM atlas w ranges do not match its height coordinate.")
    if int(payload["time_block"]) < 5 or int(payload["height_block"]) < 5:
        raise ValueError("SAM sprite blocks must contain at least five panels per axis.")
    return {
        **payload,
        "asset_base_url": ATLAS_URL,
    }


def _validate_sprite_files(manifest: dict) -> None:
    """Ensure a copied atlas is complete before exposing the browser controls."""

    missing = []
    for time_start in range(0, len(manifest["time_seconds"]), int(manifest["time_block"])):
        for height_start in range(0, len(manifest["height_m"]), int(manifest["height_block"])):
            filename = f"atlas_t{time_start:04d}_z{height_start:03d}.png"
            if not (ATLAS_DIR / filename).is_file():
                missing.append(filename)
    if missing:
        preview = ", ".join(missing[:3])
        suffix = "…" if len(missing) > 3 else ""
        raise FileNotFoundError(
            f"SAM atlas is incomplete at {ATLAS_DIR}: missing {preview}{suffix}. "
            "Click Generate / refresh atlas."
        )


@lru_cache(maxsize=1)
def _load_manifest() -> dict:
    if not MANIFEST_PATH.exists():
        raise FileNotFoundError(
            f"No pre-rendered atlas exists at {MANIFEST_PATH}. "
            "Run: python -m dash_app.misc_tab.sam_w_rt_neighborhood.atlas"
        )
    manifest = _validated_manifest(json.loads(MANIFEST_PATH.read_text(encoding="utf-8")))
    _validate_sprite_files(manifest)
    return manifest


def _controls(manifest: dict):
    times = manifest["time_seconds"]
    heights = manifest["height_m"]
    time_low, time_high = _valid_center_bounds(len(times))
    height_low, height_high = _valid_center_bounds(len(heights))
    time_value = int(
        np.clip(
            _nearest_index(times, DEFAULT_TIME_MINUTES * 60.0),
            time_low,
            time_high,
        )
    )
    height_value = int(
        np.clip(_nearest_index(heights, DEFAULT_HEIGHT_M), height_low, height_high)
    )
    return html.Div(
        [
            html.Div(
                [
                    html.Div(
                        [html.Span("Center time"), html.Strong(id=TIME_LABEL_ID)],
                        className="sam-wrt-control-heading",
                    ),
                    dcc.Slider(
                        id=TIME_SLIDER_ID,
                        min=time_low,
                        max=time_high,
                        step=1,
                        value=time_value,
                        marks=_slider_marks(
                            times,
                            time_low,
                            time_high,
                            lambda value: f"{float(value) / 60.0:.0f}m",
                        ),
                        updatemode="drag",
                        persistence=True,
                        persistence_type="local",
                    ),
                ],
                className="sam-wrt-control",
            ),
            html.Div(
                [
                    html.Div(
                        [html.Span("Center altitude"), html.Strong(id=HEIGHT_LABEL_ID)],
                        className="sam-wrt-control-heading",
                    ),
                    dcc.Slider(
                        id=HEIGHT_SLIDER_ID,
                        min=height_low,
                        max=height_high,
                        step=1,
                        value=height_value,
                        marks=_slider_marks(
                            heights,
                            height_low,
                            height_high,
                            lambda value: f"{float(value):.0f}m",
                        ),
                        updatemode="drag",
                        persistence=True,
                        persistence_type="local",
                    ),
                ],
                className="sam-wrt-control",
            ),
        ],
        className="sam-wrt-controls",
    )


def _sprite_grid():
    cells = []
    for index in range(25):
        classes = "sam-wrt-sprite-cell"
        if index == 12:
            classes += " sam-wrt-sprite-center"
        cells.append(
            html.Div(
                html.Div(
                    id={"type": SPRITE_LABEL_TYPE, "index": index},
                    className="sam-wrt-sprite-label",
                ),
                id={"type": SPRITE_TYPE, "index": index},
                className=classes,
            )
        )
    return html.Div(cells, className="sam-wrt-sprite-grid")


def build_layout():
    """Build a fixed sprite grid; only the background crop changes while browsing."""

    try:
        manifest = _load_manifest()
        body = [
            _controls(manifest),
            html.Div(id=STATUS_ID, className="sam-wrt-status"),
            dcc.Store(id=MANIFEST_STORE_ID, data=manifest, storage_type="memory"),
            _sprite_grid(),
        ]
    except (FileNotFoundError, json.JSONDecodeError, OSError, TypeError, ValueError) as error:
        body = [
            html.Div(str(error), id=STATUS_ID, className="sam-wrt-status-error"),
            empty_state(
                "Generate the SAM image atlas once",
                "Run `python -m dash_app.misc_tab.sam_w_rt_neighborhood.atlas`; the report then "
                "browses static tiles without reopening NetCDF files or redrawing plots.",
            ),
        ]

    return html.Article(
        [
            dcc.Location(id=GENERATION_RELOAD_ID, refresh=True),
            dcc.Interval(
                id=GENERATION_INTERVAL_ID,
                interval=1_000,
                n_intervals=0,
                disabled=True,
            ),
            report_header(
                "SAM w–rₜ time-height neighborhood",
                "Twenty-five consecutive raw SAM planes from a pre-rendered image atlas. "
                "Time increases left to right, altitude increases bottom to top, and the "
                "cyan frame marks the slider-selected center.",
                eyebrow="Raw 3-D SAM explorer",
                badges=("ARM", "5×5 neighborhood", "Static atlas"),
            ),
            html.Div(
                [
                    html.Button(
                        "Generate / refresh atlas",
                        id=GENERATE_ID,
                        n_clicks=0,
                        className="coherent-action-button",
                    ),
                    html.Div(
                        "The atlas is generated from ARM raw-SAM snapshots. "
                        "Use this when it is absent or you want fresh tiles.",
                        id=GENERATION_STATUS_ID,
                        className="coherent-control-note",
                    ),
                ],
                className="sam-wrt-generation",
            ),
            *body,
            html.Div(
                [
                    html.Strong("Color: "),
                    "gold is populated air, blue is positive local w′rᶜ′ transport, red is "
                    "negative local w′rᶜ′ transport, and pale bins contain both. Each height "
                    "keeps the same axes and color scale through every saved SAM minute.",
                ],
                className="sam-wrt-caption",
            ),
        ],
        className="notes-report sam-wrt-report",
    )


def register_callbacks(app):
    """Serve atlas PNGs and crop them entirely in the browser."""

    @app.callback(
        Output(GENERATION_STATUS_ID, "children"),
        Output(GENERATION_INTERVAL_ID, "disabled"),
        Output(GENERATION_RELOAD_ID, "href"),
        Input(GENERATE_ID, "n_clicks"),
        Input(GENERATION_INTERVAL_ID, "n_intervals"),
    )
    def generate_or_poll(_clicks, _ticks):
        if ctx.triggered_id == GENERATE_ID:
            return start_atlas_generation(), False, no_update
        running, succeeded, message = job_status("sam_wrt_atlas")
        if running:
            return message, False, no_update
        if succeeded:
            _load_manifest.cache_clear()
            return message + " Reloading the report…", True, "/"
        return message or no_update, True, no_update

    endpoint = "notes_sam_wrt_atlas"
    if endpoint not in app.server.view_functions:
        app.server.add_url_rule(
            f"{ATLAS_URL}/<path:filename>",
            endpoint,
            lambda filename: send_from_directory(
                str(ATLAS_DIR), filename, max_age=86400
            ),
        )

    app.clientside_callback(
        """
        function(timeIndex, heightIndex, manifest) {
            if (!manifest || timeIndex === null || heightIndex === null) {
                return [
                    window.dash_clientside.no_update,
                    window.dash_clientside.no_update,
                    window.dash_clientside.no_update,
                    window.dash_clientside.no_update,
                    window.dash_clientside.no_update,
                    "Atlas unavailable"
                ];
            }
            const radius = 2;
            const times = manifest.time_seconds;
            const heights = manifest.height_m;
            const timeBlock = Number(manifest.time_block);
            const heightBlock = Number(manifest.height_block);
            const pad = (value, width) => String(value).padStart(width, "0");
            const version = encodeURIComponent(manifest.generated_utc || manifest.version);
            const styles = [];
            const labels = [];
            const titles = [];
            const urls = new Set();

            timeIndex = Math.max(radius, Math.min(times.length - radius - 1, Number(timeIndex)));
            heightIndex = Math.max(radius, Math.min(heights.length - radius - 1, Number(heightIndex)));
            for (let row = 0; row < 5; row += 1) {
                const zIndex = heightIndex + radius - row;
                for (let column = 0; column < 5; column += 1) {
                    const tIndex = timeIndex + column - radius;
                    const timeStart = Math.floor(tIndex / timeBlock) * timeBlock;
                    const heightStart = Math.floor(zIndex / heightBlock) * heightBlock;
                    const filename = `atlas_t${pad(timeStart, 4)}_z${pad(heightStart, 3)}.png`;
                    const url = `${manifest.asset_base_url}/${filename}?v=${version}`;
                    urls.add(url);
                    const localTime = tIndex - timeStart;
                    const localHeight = zIndex - heightStart;
                    styles.push({
                        backgroundImage: `url("${url}")`,
                        backgroundSize: `${timeBlock * 100}% ${heightBlock * 100}%`,
                        backgroundPosition: `${100 * localTime / (timeBlock - 1)}% ${100 * localHeight / (heightBlock - 1)}%`,
                        backgroundRepeat: "no-repeat"
                    });
                    const minutes = Number(times[tIndex]) / 60.0;
                    labels.push(`${minutes.toFixed(0)} min · ${Number(heights[zIndex]).toFixed(0)} m`);
                    const xr = manifest.x_range_rt_g_per_kg[zIndex];
                    const yr = manifest.y_range_w_m_per_s[zIndex];
                    titles.push(
                        `SAM ${minutes.toFixed(0)} min, ${Number(heights[zIndex]).toFixed(0)} m\n` +
                        `r_t ${Number(xr[0]).toPrecision(4)} to ${Number(xr[1]).toPrecision(4)} g/kg\n` +
                        `w ${Number(yr[0]).toPrecision(4)} to ${Number(yr[1]).toPrecision(4)} m/s`
                    );
                }
            }
            const result = [
                styles,
                labels,
                titles,
                `${Number(times[timeIndex]) / 60.0} min`,
                `${Number(heights[heightIndex]).toFixed(0)} m`,
                ""
            ];
            // Apply all 25 crops together after their few shared atlas PNGs are cached.
            const loads = Array.from(urls).map((url) => new Promise((resolve) => {
                const image = new Image();
                image.onload = resolve;
                image.onerror = resolve;
                image.src = url;
            }));
            return Promise.all(loads).then(() => result);
        }
        """,
        Output({"type": SPRITE_TYPE, "index": ALL}, "style"),
        Output({"type": SPRITE_LABEL_TYPE, "index": ALL}, "children"),
        Output({"type": SPRITE_TYPE, "index": ALL}, "title"),
        Output(TIME_LABEL_ID, "children"),
        Output(HEIGHT_LABEL_ID, "children"),
        Output(STATUS_ID, "children"),
        Input(TIME_SLIDER_ID, "value"),
        Input(HEIGHT_SLIDER_ID, "value"),
        State(MANIFEST_STORE_ID, "data"),
    )
