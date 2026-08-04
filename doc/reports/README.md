# Static CLUBB reports

This directory holds static investigation bundles shown in Dash’s **Reports** tab. A report keeps its conclusion, figures, code excerpts, supporting data, and provenance together and can be opened directly as HTML.

Do not create a Python report module or edit `dash_app/` to add a report. The dashboard polls `index.json`; a new indexed bundle appears within a few seconds without a Dash reload.

## Bundle contract

```text
doc/reports/
  index.json                         # ordered, published-report catalog
  <report-id>/
    manifest.json                    # same report metadata, stored with evidence
    report.html                      # standalone static page
    figures/                         # copied PNG/SVG/PDF-derived figures
    snippets/                        # copied code excerpts or commands
    data/                            # compact supporting CSV/JSON/etc.
```

`index.json` is the only discovery mechanism. Its `reports` list contains the report ID, title, summary, timestamp, tags, static entrypoint, and optional source revision. A report is visible only when both the index entry and entrypoint exist.

Report IDs use lowercase letters, digits, and hyphens. A report can be deliberately updated in place without restarting Dash: use the same ID and `replace=True`. The publisher swaps the complete staged bundle into that directory, updates its one catalog entry, and the agent's `open_report` handoff forces the iframe to reload. Use a new ID only for a genuinely separate investigation.

## Temporary staging is not report evidence

`output/agent_artifacts/` is a disposable, ignored staging area for live
agent/MCP work. It may contain transient PNGs, compact manifests, logs, and
isolated SCM output while a job is running, but it is not a durable source for
a report and may be cleared after jobs finish. Never reference it from
`report.html` or leave a report dependent on it.

During an investigation, an agent may render or analyze material there. Before
publishing, copy every figure, compact data file, code snippet, and exact
reproduction script worth preserving into this report bundle. A report reader
must be able to open the bundle after `output/agent_artifacts/` has been
deleted. Use a named directory under `output/` instead if a large raw model
run, rather than a static report asset, needs to remain available.

## Rendered evidence is part of the report

If a report uses generated figures, commit the final PNG/SVG files under that
report's `figures/` directory and reference them with relative paths from
`report.html`.  The dashboard only serves these static artifacts: it never
runs a plotting script when a reader opens a report.  This keeps the evidence
viewable after source data move or a Python environment changes.

### Raster-size policy

Published PNG/JPEG evidence is capped at **1200 pixels in either dimension**.
`ReportBuilder.figure(...)` enforces this automatically while staging a new or
replacement bundle, preserving aspect ratio and using high-quality Lanczos
resampling. SVG remains vector evidence and is not rasterized. Keep plotting
scripts at a sensible DPI as well—the cap protects the published bundle, not
the temporary image they create.

To bring legacy report assets into the same policy, run the checked-in
maintenance command from repository root:

```bash
PYTHONPATH=. .venv-dash/bin/python doc/reports/resize_raster_assets.py
```

It changes only report PNG/JPEG files whose largest dimension exceeds 1200.

Keep the documented reproduction script beside the report (or preserve the
exact command and inputs in `snippets/`), and commit its compact output
metadata under `data/`.  A complete generated-figure report therefore contains
all three pieces: the embedded rendered image for readers, the script for a
deliberate rerender, and the metadata needed to identify the source plane and
normalization.  Re-run the script only when intentionally revising the report,
then inspect and commit the updated assets along with the HTML.

## Preferred authoring path

Use the small builder after completing the analysis. It stages the entire bundle under `doc/reports/.staging/`, then atomically moves it into place and rewrites the catalog only once the HTML and metadata are complete.

```python
from dash_app.reports_tab.publisher import ReportBuilder

report = ReportBuilder(
    "advance-xp3-stability-2026-07-22",
    "RTP3 and THLP3 stability with advance_xp3",
    summary="Tests whether third-moment behavior remains bounded and informative.",
    tags=("advance_xp3", "rtp3", "thlp3", "ARM"),
    source_revision="<git commit or working-tree note>",
)
report.heading("Answer at a glance")
report.callout("Conclusion", "State the conclusion before the supporting detail.", tone="success")
report.metrics([
    ("Cases", "ARM + BOMEX", "same configuration used for the figures"),
    ("Runs", "8", "explicitly listed in the appendix"),
])
report.figure("output/my_diagnostic.png", caption="An informative caption.")
report.equation(r"\overline{w' r_t'^2}", caption="Example MathText equation rendered to SVG.")
report.code_file("src/some_routine.F90", language="fortran", caption="Relevant closure branch")
report.publish()
```

For an update to the same investigation, retain the report ID:

```python
report = ReportBuilder(
    "advance-xp3-stability",
    "RTP3 and THLP3 stability with advance_xp3",
    summary="Revised results with a corrected time alignment.",
    replace=True,
)
```

The builder uses only HTML, copied artifacts, and JSON. Its equation helper renders MathText to SVG so reports remain static and do not depend on JavaScript or external CDNs.

## After publishing: open the report in Dash

Creating the bundle is not the final agent step.  After `report.publish()`
returns, open the exact updated report in the active dashboard so the user
can inspect it immediately:

```text
invoke_dashboard(
  tab="reports",
  operation="open_report",
  arguments={"report_id": "advance-xp3-stability-2026-07-22"},
)
```

The Reports tab receives this semantic request, selects the report, and loads
its static `report.html`.  Do this only after publication succeeds: a staging
directory is deliberately not discoverable or served.  If no dashboard agent
connection is active, still publish the report and state its report ID and path
in the normal agent reply; Dash will add it from `index.json` within a few
seconds once opened.

## Recommended report shape

1. Answer at a glance.
2. Question and configuration.
3. Evidence and figures.
4. Method, equations, and selected code.
5. Physical interpretation and implementation implications.
6. Limitations and next tests.
7. Reproducibility appendix: runs, flags, input files, data paths, and source revision.

The checked-in [`hello-world`](hello-world/report.html) bundle is a compact visual example.
