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
models, typed overrides, validated stats categories or local namelist choices, explicit physical time
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

## Startup and page loading

Hidden `dcc.Tab` children still participate in Dash's callback graph. Use
`lazy_tabs.LazyTabs` to register page builders without constructing their bodies
at startup. Register every callback and route before serving the app; only
layout construction is deferred. A browser-side gate requests each page once,
including a persisted initial selection. The loader applies workspace
persistence to new controls and keeps visited bodies in the layout, preserving
form/runtime state across navigation. Loaded/request stores are memory-only.
An explicit dashboard handoff also loads its target, even with `preserve_tab`,
so requests can update an unvisited workspace without changing the visible tab.

Compile, Run, Profile and Tune load independently. The shared runner selector
lives outside tab bodies and owns its catalog, runtime preferences and metadata.
Consumers use its published stores; optional badge inputs tolerate unmounted
tabs. Named-config saves publish the shared catalog rather than writing another
tab's rendered controls. Tutorial and
Misc have separate gates for each lesson/tool, also requiring the parent tab to
be active. Component JavaScript still loads eagerly: lazy bundles in the current
Dash version caused slider-tooltip callbacks to access unmounted React nodes
during rapid navigation. Defer page bodies, but retain `eager_loading=True`
until that browser regression is resolved.

First visits can still take time to initialize a selected diagnostic. This is
separate from broker-owned background work; do not enable Flask request threads
to work around NetCDF loading. When adding a cross-tab callback, check that its
inputs, state, and outputs exist together, including on a fresh browser page.

## Avoiding repeated work

Profile publishes a compact result signature when its manifest/CSV files or job
state change. Its graphs depend on that signal and comparison controls, not on
every log update or interval tick. File signatures include modification time,
change time, size, and inode; timestamps cross the browser boundary as strings
because JavaScript numbers cannot preserve nanoseconds. Timing CSVs derive the
summary and process tables together. Keep file-change detection working for
same-size overwrites and atomic replacements, not just appended rows.

Documents and Tune workspace activity polling pause outside their tabs and refresh
on reentry. Tune job supervision/leases continue independently. Activity state
still uses its lock and atomic writes, but unchanged transactions do not rewrite
the file. Trajectory discovery caches only compatibility booleans in a bounded
LRU keyed by file version; it retains no NetCDF handles/arrays and retries I/O
errors. Resize handlers attach window listeners only during a drag, so tab
remounts cannot accumulate handlers retaining old layouts.

See [PERFORMANCE_AUDIT.md](PERFORMANCE_AUDIT.md) for measurements and remaining
opportunities from the September 2026 audit.

## Plot data and workspace state

`plots-case-data` is a memory-only browser descriptor. Large time arrays and
variable catalogs live in the process-shared `plot_tab/case_cache.py` cache;
callbacks needing those fields call `resolve_case_data` before using them.
UI-neutral services still accept and return full metadata. Persist only intent
in `plots-case-selection`, not the derived snapshot or its cache key. Keep the
legacy workspace migration in `assets/37_workspace_persistence.js` compatible
with exported workspaces.

Read time/column selections before materializing NetCDF fields. Cache reduced
profiles and images through the byte-bounded extraction cache, keyed by file
signature and selection. Do not reintroduce full-cube preloads or request
threads for HDF5. Profile, budget, and subcolumn playback can update trace values
with Dash patches while preserving manual zoom.

`plot_tab/output_chooser.py` owns the output chooser's two-sided layout, stores,
discovery, refresh, extra paths and confirmed deletion. Its browser renderer and
selection transitions live in `assets/34_plot_chooser.js`, with styles in
`assets/34_output_chooser.css`. The component publishes the ordered list in
`plots-output-dirs.data`; case loading and PyPlotGen consume that list, without
reading its catalog or manipulating its internal display. The chooser does not
call the Plot task runtime; the case callback invalidates old plot requests when
it receives a changed selection. Saved selections retain the same persistent
store id, including temporarily unavailable folders.

The chooser reserves a fixed compact height, with a five-entry Available list
and wrapping Active bubbles. Its inner surface is absolutely positioned in both
states; expansion grows that surface downward over adjacent content without
changing the chooser or parent's layout height. Refresh/Add extra are top-bar
siblings of the expansion button. Add extra expands the surface and reveals its
path field; collapse hides that field. Expansion is memory-only component state
and does not change selection or trigger case loading. Expanded lists scroll
after their viewport cap, with guarded delete controls available there. Browser
rendering remains the sole list renderer. Catalog scans run in the existing metadata worker on
mount/refresh; selection changes inspect only unknown paths once a catalog
exists. Initial hydration with an empty catalog performs full discovery so a
restored selection cannot hide other available folders. Internal catalog
refreshes never prune the selected list or trigger scientific case reloads.

Directory selection and chooser rendering run in the browser. Register expensive
native Plot callbacks with `plot_tab/async_callbacks.py:task_callback`; its fixed
task names execute in broker-owned spawned processes through `plot_tab/tasks.py`.
There are separate single-worker lanes for metadata, figures and LES acquisition,
so downloads do not delay figures or case discovery. Warm workers retain extraction caches. The
same-origin `/plots/tasks` proxy accepts only submit/poll for these registered
tasks; it does not expose broker credentials or general function execution.

Browser callbacks return immediately, keep existing figures visibly marked, and
batch polling across pending cards without repeating the broker status probe.
Unchanged control values do not resubmit figures. Keep cards mounted on ordinary
case switches; grid callbacks need card ids, not serialized figure children.
Page/card revisions and selection generations
discard superseded results; queued work is replaced or dropped before execution.
The broker retains at most 128 request scopes and 64 MiB of serialized results,
expiring inactive records after two minutes. A missing record is resubmitted,
including after a broker restart. These transient UI tasks do not create durable
scientific job/artifact records. Broker shutdown closes their worker pools.
Playback patches require a successful accepted figure with matching inputs;
worker render flags alone cannot establish a valid patch base across browsers.

SAM/COAMPS overlay metadata retains the canonical collection paths from the
case definitions, including missing files. The `benchmark_sam` and
`benchmark_coamps` tasks call `benchmark_overlay.acquire_benchmark_source`, which
uses the shared `confirm_les_data` helper in the LES lane. The browser starts a
task for a selected source or a custom expression referencing LES, and shows its
loading/error state beside the source button. These tasks share the existing
revision checks, bounded queues and batched polling; no separate polling loop
or job store is needed.

`plots-enabled-benchmarks` retains selection intent. `plots-ready-benchmarks`
carries current-case readiness and the prepared subset of that selection;
figure callbacks and time-height sizing use this subset. Custom expressions
also accept their prepared source without enabling an overlay on other plots.
`benchmark_overlay._source_paths` uses the shared read-only `find_les_data`, so
rendering never downloads. Missing files omit that source until preparation
completes. Readiness changes rerender plots automatically, with case/path checks
preventing stale results from enabling another case. Failed preparation remains
visible at the button and off/on retries it. A running download may complete
after deselection; its result cannot re-enable the overlay, and its file remains
reusable. Pyplotgen's own data lookup remains unchanged.

The PyPlotGen header category uses `pyplotgen_runtime.pyplotgen_help()` to run
the actual CLI with `-h`. The `pyplotgen_help` metadata task loads it without
blocking tab mounting and survives case-selection revisions, like the output
catalog. CLI text is cached per process until the entry-point signature changes.
Its verbatim stdout supplies the help notecard; option declaration lines supply
names, value hints and wrapped descriptions for the searchable multi-selector.
No second PyPlotGen option inventory is maintained. Dashboard-managed input,
destination and help options are excluded from selection.

Selected names and value fields pass through the existing `launch_pyplotgen_request`
broker owner. `pyplotgen_arguments` checks names against the help-derived catalog,
requires nonoptional values, and splits quoted value text into literal argv tokens.
It rejects flags in value fields so those cannot override managed destinations.
The PyPlotGen parser retains responsibility for types, choices and option
compatibility. Generation still uses its existing process, progress, cancellation
and fixed output-serving lifecycle; help discovery does not generate plots.

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
`dash_app/auto_llm_generated_tests/pytests/unit_tests/test_equations.py` together.

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
focused tests in `dash_app/auto_llm_generated_tests/pytests/unit_tests/`.

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

## Document browser

The Documents tab lists `.pdf`, `.md` and `.txt` files recursively under `doc/`, including
`README.md` but excluding hidden paths. The list groups PDFs, Markdown and plain
text in that order, with alphabetical paths inside each group. Its search matches filenames and relative paths.
Other formats, including notebooks and HTML, are ignored. Discovery needs no
index, manifest or application change when adding a document.

Complete an investigation first, then save its final PDF, Markdown or plain-text document
under `doc/`. Keep raw runs and other generated research artifacts in a named
`output/` directory. The viewer reads files; it never runs research code.
PDFs use the browser viewer and Markdown uses Dash's renderer with HTML disabled.
Plain text uses a wrapping preformatted element with its original spacing and
line breaks; both text formats are read as UTF-8.
The document wrapper records its repository-relative path. The delegated link
handler resolves rendered Markdown anchors against that path, including reference
links, without parsing or rewriting Markdown source. `/repository-files/` serves
visible `doc/` files and Git-tracked repository files, rejecting outside paths and
hidden files. Source and other UTF-8 text is served as `text/plain` with sniffing
disabled; binary text targets are rejected. Source code is never executed.
The catalog refreshes while this tab is visible. The semantic browser handoff is
`reports/open_document` with a `path` relative to `doc/`.

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

## Runtime selection and JAX job settings

[`compile_tab/runner_selector.py`](compile_tab/runner_selector.py) owns the shared
chooser. Its menu opens and closes in the browser. Build discovery and JAX
inspection use the existing bounded native task transport in a separate runner
worker lane at page initialization. The root owns `compile-discovery` and
`compile-build-statuses`; Compile reads and updates those same stores. Mounting
Compile does not scan again. A completed compile refreshes the catalog through
Compile when mounted, or through the shared worker otherwise. Deletion updates
it directly. Refresh runners explicitly rescans changes made outside Dash.
Opening the menu, changing implementation/device or changing allocation
preferences submits no metadata task; there is no expiry timer.

The metadata callback builds both Fortran and JAX panels and their badge views.
The selector asset handles native click/change events, paints those cached
panels directly and publishes the existing preference stores after the first
frame. Pending values are flushed before any other control's click, so Launch
reads the current target. Metadata updates merge pending preferences instead of
restoring an older selection. Selection buttons do not feed Dash callbacks:
waiting for even clientside callbacks was visible with Run's large form loaded.
For remaining pattern-button callbacks, `clicked_trigger_id` resolves each
positive event's own ID; a zero-click mount at the front of a batch is not the
user's selection.

The existing handoff heartbeat checks managed environment marker timestamps.
A changed stamp requests only installed-package readiness through
[`refresh_environment`](../clubb_jax/runtime_info.py), preserving discovered
hardware. CUDA tile compatibility is derived from the cached inventory using
the backend's existing compatibility check. The selected install alias remains
authoritative for launches. A newer Compile catalog invalidates any older
in-flight discovery result.

The runtime chooser uses the read-only report from
[`clubb_jax/run_jax.py`](../clubb_jax/run_jax.py). Hardware discovery does not
require an installed JAX environment. A compatible device whose environment
needs setup remains selectable.

GPU choices are stored by UUID. Until a specific GPU is chosen, the inherited
selection retains the server's `CUDA_VISIBLE_DEVICES` setting (or CUDA's
default); an unmapped active default has its own tile. Explicit selections
override visibility for that job and all its cases or profiling workers. They
do not distribute work across GPUs. CPU runs ignore a saved GPU selection.

The preallocation checkbox is enabled only for a compatible CUDA selection,
including one that needs environment setup. It is disabled during checks and
for CPU, Metal, or unavailable selections. Both checked and unchecked values
explicitly override the server's `XLA_PYTHON_CLIENT_PREALLOCATE` setting. Freeze
the effective GPU and preallocation settings into each submitted Run/Profile/Tune
request and include them in command previews. The launcher arguments carry both enabled and disabled preallocation choices.

Compiled implementations follow `install/selected`, falling back to
`install/latest`. Explicit runner `-exe` or `-install_dir` options override
that default. Tune converts the chooser selection into the canonical tuner
`backend`, `jax_options`, `jax_gpu` and `jax_xla_prealloc` request fields. The
broker, typed MCP service and native Start callback preserve the same fields.
Loading a saved revision restores its runtime; Continue and leaderboard reruns
use the immutable request rather than the current chooser. The broker owns
result-run processes and terminal polling as well as tuning workers.

Tune algorithm buttons and typed option panels are generated from
[`tuner/algorithms/algorithms.json`](../tuner/algorithms/algorithms.json), read by
[`tuner/algorithm_registry.py`](../tuner/algorithm_registry.py). Algorithm classes
live under `tuner/algorithms/`; the Tune callbacks only map wildcard input
values to canonical `strategy.options`. Saved requests hydrate those same
controls. The Algorithm question-mark help reads
[`tuner/algorithms/README.md`](../tuner/algorithms/README.md) on each open and renders it in the shared notecard; repository links point to
the checkout commit on GitHub. Algorithm explanations stay in that README.
The broker and typed service forward `strategy_options` through the
registry, so a new algorithm file and JSON entry do not require Dash-specific
validation or option callbacks. Algorithm requirements enable mandatory JAX
gradients before launch; the shared runtime chooser still selects the backend.
The native form derives gradient computation from that requirement and only
enables its mode selector for gradient algorithms. Its preview validates the
actual shared runtime selection: incompatible algorithms are red, their tooltip
and summary explain the requirement, and Start is disabled. Running and readonly
revisions validate their saved request's backend instead of the shared chooser.
Numeric options use compact
manual entries beside their labels. Shared notecards explain gradient modes,
RMSE weights and across-case aggregation.

Taylor diagrams, parameter spread and parameter landscapes are permanent.
Only correlations, sensitivity and interactions have optional controls. The
`tune-diagnostic-plots` store combines the fixed primary set with those choices;
individual Tune controls are not restored from component persistence.

Plots refresh automatically through the native task transport's separate Tune
lane. `read_tuning_status` adds a compact history-publication signature from
chunk file metadata, without loading arrays. Landscape/control tasks follow
that signature; Taylor/spread tasks follow retained leaderboard data. Timer-only
status changes do not submit work. An in-flight data refresh finishes instead
of being continually superseded; its latest queued data update runs immediately
afterward, including when final-job polling has already stopped. Explicit
control changes, revision changes and Refresh now still request a fresh view.
The worker computes only enabled plots and retains at most 64 MiB of decompressed
history for one job. Hidden tabs do no plot work; browser-local uptime remains
independent, and status polling pauses while another tab is visible.

Sample-history publication remains unchanged: the scheduler publishes the first completed
set immediately, then pending rows after 500 samples or 20 seconds, whichever
comes first, and flushes them on stopping.

`clubb_jax/run_jax.py` owns the public selection interface: `runtime_selection`
validates saved settings, `runtime_arguments` encodes one opaque selection value
for `-options` or `-jax`, and `runtime_configuration` supplies launch provenance without
setup. Discovery uses the same arguments with `-info=json`; its report contains
available devices, the resolved selection and runtime readiness. Dash and the
tuner forward settings rather than constructing GPU environment variables.
The launcher delegates device visibility and memory settings to its backends
before inspection or importing JAX. Rerun parameter/physics override precedence
is shared through `input/scripts/create_case_namelist.py`.

## Profile results

The Profile tab is a browser interface to `utilities/time_clubb.py`. Its top
benchmark panel configures the case, process/per-process-batch-size sweep,
repetitions, executable, configuration, overrides, and additional
`run_scm.py` arguments. While visible, the tab checks active results once per
second and refreshes figures when timing files change; warmups remain hidden.
Logs and job status update independently, without rebuilding unchanged plots.
Browser stores retain compact progress and timer choices rather than the
growing raw timing table. Stored profiles can be overlaid, compared with a baseline, or
viewed as process distributions and exclusive-cost decompositions. The right
rail has a profile-selection section above a separate set of shared comparison
controls; plot-specific options remain beside the plot they affect. The
profile chooser shows the three newest unselected results by default, expands
to the full library, and displays active comparisons as removable pills. The
benchmark-label field indicates when its normalized profile name already
exists. Starting that benchmark asks for confirmation, then replaces the
existing profile in place instead of creating a timestamped version; any older
same-label/same-case versions are removed from the active comparison selection.

The selected directory is a collection of compact, directly commit-able
profile folders. Each benchmark creates `<profile-name>/` containing
`README.md`, profile-wide `profile.json` provenance, one workload row per
process-count/batch-size point in `batches.csv`, raw timer observations in
`timings.csv`, and one representative input/setup/log/native-timing set under
`logs/<batch-id>/`. Child processes otherwise run in temporary directories,
which are deleted after each workload is aggregated. Warmups are retained with
`phase=warmup` but excluded from the default plots. Dash derives statistical
summaries in memory instead of storing duplicate summary files. Failed warmups
and repetitions show a child log excerpt in the console and retain its complete
runner/model diagnostics in a `*_failed.log` file, even if another process
succeeded. **Export
selected** downloads complete profiles as a ZIP; **Import** accepts those ZIPs
on another machine or checkout. Provenance includes the effective vertical
level count, observed model steps, source revision, executable checksum, host,
timer backend, and time basis so Dash can flag potentially incomparable runs.
