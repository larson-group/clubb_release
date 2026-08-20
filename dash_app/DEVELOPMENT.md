# Dash App Development Notes

## Design Intent

The dash app should be an interface for doing command-line work quickly. Features should remain implemented and accessible through the CLI first; the dash app should make those commands easier to run, easier to parameterize, and easier to inspect in different views.

Avoid making the dash app the only place a workflow exists. If a workflow matters, make sure there is a script, command, or CLI path that can run it without the app.

## Dropdown Styling

Some Dash dropdowns can render a small input box around the first letter typed into the search field. The compile tab's extra modules dropdown has hit this because it is a searchable `dcc.Dropdown`.

Fix this the same way the other tabs handle dropdowns: scope CSS to the tab/container and clear the native input styling on `.Select-input > input`.

```css
#app-root .some-container .dash-dropdown .Select-input > input {
  caret-color: transparent;
  background: transparent !important;
  border: none !important;
  padding: 0 !important;
  margin: 0 !important;
  box-shadow: none !important;
  outline: none !important;
}

#app-root .some-container .dash-dropdown .Select-input > input:focus {
  outline: none !important;
  box-shadow: none !important;
  border: none !important;
  background: transparent !important;
}
```

Also add the matching single-value/placeholder padding rule when needed so the text aligns with neighboring controls.

New shared dropdowns should use `dash_app.shared.components.styled_dropdown`.
It applies the common `clubb-dropdown` repair while leaving each tab's ids,
options, callback contract, and layout local.  Do not replace a working
tab-specific dropdown merely for consistency; use the wrapper when adding a
control or removing a demonstrated duplicate styling fix.

## Service boundaries

Put reusable UI-neutral validation, models, time-window selection, provenance,
artifact, and job lifecycle code in `dash_app/shared/` or `dash_app/services/`.
Tabs own their controls, figures, and callbacks, but call the service for work
that an MCP client could also request.  Extract existing behavior before moving
callers so Dash remains the visual regression baseline.  Shared helpers should
be small and named by a real repeated contract; `shared/` is not a catch-all
directory.

Agent-owned execution is local-only. Public MCP requests use strict Pydantic
models, typed overrides, checked-in stats choices, explicit physical time
windows, and request-id idempotency. A typed MCP SCM batch writes scientific
output by default to `output/mcp_runs/<batch-id>/<case>_stats.nc`; the broker
owns that batch directory, while the JobStore is the canonical batch/child
record. MCP `out_dir` may select another path, but it must resolve below the
repository `output/` directory. Native Run multi-select uses the same service
but honors the Run-tab's repository-resolved output directory. Admission is
atomic and scoped by canonical `(case, output directory)`, so one case may run
in different directories concurrently but may not overlap in the same one.
`output/agent_artifacts/` holds private manifests and temporary evidence only:
protect a live bundle from cleanup, but do not use it as a durable result,
report asset, or raw-experiment archive.
Do not add new free-form shell, path, or NetCDF-returning endpoints. The
generic browser handoff actions remain compatibility-only and must not become a
second job runner.

### Local MCP endpoint lifecycle

Each checkout has one broker-owned Streamable HTTP MCP endpoint. The broker
persists its loopback URL, random bearer credential, and endpoint identity under
the private runtime directory, so an ordinary Dash restart does not require
reconfiguring the Codex chat. The dashboard utility drawer shows the details
once the broker migration is ready. Dash registers its PID/start time and
heartbeat with the broker; browser-targeted operations return
`DASHBOARD_UNAVAILABLE` when no current browser is registered. Compile, Run,
Tune, and artifact jobs remain broker-owned and continue without a browser.
Under `./launch_dashboard.sh`, a foreground manager supervises both processes:
it restarts an unhealthy Dash every 10 seconds for at most 5 minutes, while the
broker enforces a 30-second manager lease so an abruptly lost launcher cannot
leave the dashboard or scientific workers running indefinitely.

The older `mcp-instances/` rotating records are migration-only. A broker
startup terminates only the exact recorded legacy endpoint process (matching
PID and process start time), removes its private record, and starts the stable
endpoint without changing its persisted URL or bearer credential. Multiple live
Dash processes are rejected rather than selected unpredictably.

## Controls

Do not add extra `+` and `-` buttons for entry boxes unless they are truly useful for the workflow. Numeric inputs and text boxes should stay simple when direct editing is sufficient.

## Notecards

In this app, a **notecard** means the reusable modal-style overlay used for
context that should not permanently crowd a tab: plot `?` explanations,
tutorial background notes, detailed warnings, and logs such as a failed source
standards check. It is not an ordinary card embedded in the page layout.

Build these overlays with `dash_app.shared.notecard.notecard`; use
`information_body` when the content follows the usual concise overview followed
by detailed sections and nuances. The component supports `small`, `medium`,
`large`, and `full` sizes. Keep the open/close callback beside the feature that
owns the content, while keeping the presentation and accessibility behavior in
the shared component. New code should use “notecard” in function and component
names.

## Tutorial Equation Guide

The interactive CLUBB equation reference lives in
`dash_app/tutorial_tab/clubb_equations_demo/` and is embedded by
`dash_app/tutorial_tab/clubb_equations.py`. Keep authoritative equation and
term metadata in `content.py`, reusable renderers in `components.py`, and Dash
callbacks/layout assembly in `app.py`. Each clickable term needs a unique
occurrence ID because the same budget label can appear in several equations.

Equation colors follow the official CLUBB-SILHS document rather than a local
decorative palette: red is host/microphysics/radiation, blue is returned by
CLUBB to the host, purple is CLUBB-internal, green is diagnosed from the PDF,
and brown is closed by a classical non-PDF parameterization. Physical families
such as transport, buoyancy, pressure, and dissipation are secondary inspector
badges so they do not overwrite the official color meaning.

The cards are a curated quick-reference view of `doc/CLUBBeqns.tex`, not a
second independent equation source. When a branch changes an equation, update
the guide metadata, its displayed source routine, and the focused test in
`dash_app/pytests/test_equations.py` together.

## ADG1 Explorer Teaching Model

The Tutorial ADG1 explorer uses the reusable normalized mixture mathematics in
`dash_app/shared/two_gaussian_model.py`. It is a teaching visualization rather
than a runtime PDF diagnostic.

## Misc subtabs

The Misc tab is a lightweight left-rail catalog for living investigations and
focused diagnostics. Each visible Misc page owns one direct subdirectory under
`dash_app/misc_tab/` (for example, `sam_w_rt_neighborhood/`). Its package registers a `SubtabSpec`, and the Misc shell
discovers those packages automatically for the persistent vertical directory.
Keep each subtab's data loading, callbacks, generation helpers, and focused
focused tests in `dash_app/pytests/`.

Use the small layout builders in `dash_app/shared/reporting.py` for report
headers, fact grids, sections, and intentional empty states. They are shared
presentation components, not sources of scientific truth. Expensive or
scientifically meaningful analysis must continue to have a reproducible CLI or
saved artifact; a Misc subtab should visualize and explain that result rather
than silently recompute a different dashboard-only result.

The SAM w–r_t neighborhood follows this pattern through its local
`atlas.py` generator: it writes a reproducible pre-rendered atlas while the
browser only changes image crops. Keep that generator as an importable module
and a tracked background subprocess rather than implementing unlogged,
dashboard-only generation.

## Static report bundles

The top-level Reports tab is intentionally separate from Misc.  Misc contains
interactive Python/Dash tools; Reports displays immutable HTML bundles under
`doc/reports/`.  Adding a report must not create a Python module, touch
`dash_app/assets/`, or otherwise alter the dashboard source tree.  Publish
only a complete bundle and its atomically updated `doc/reports/index.json`
catalog so a running Dash process can discover it through its lightweight poll
without a development reload.  See `doc/reports/README.md` and use
`dash_app.reports_tab.publisher.ReportBuilder` for the standard staging and
publication flow.

The supported production Misc subtab is the SAM w–rₜ time-height
neighborhood. It serves a pre-rendered sprite atlas generated by
`dash_app.misc_tab.sam_w_rt_neighborhood.atlas`; its sliders change browser-side image
crops rather than reopening NetCDF files or redrawing 25 Plotly figures. Keep
the atlas contract and generator in production because this subtab must work
without the experimental tree.

The SAM subtab's Generate / refresh control launches its generator as a
tracked background subprocess and exposes its log tail in the tab. Keep this
command in the repository rather than implementing unlogged dashboard-only
data generation.

Prototype work belongs outside `dash_app/misc_tab/`; the production app must
neither import nor register it. This keeps the dashboard usable when
prototype-only directories are absent from a checkout.
