# CLUBB Dash App

The Dash app is a browser interface for common CLUBB workflows. It provides a
Run tab for launching CLUBB cases, a Profile tab for configurable process-based
timing sweeps, and a Plots tab for inspecting CLUBB NetCDF output. The app is
intended to be run from an existing CLUBB checkout.

The Run, Profile, and Tune action areas each show the effective default CLUBB
build next to their launch buttons. This read-only badge follows
`install/selected` (or the same `install/latest` fallback used by
`run_scm.py`) and updates when the Compile tab selects another build. Hover it
for the resolved install and CMake paths, Fortran compiler, build type,
precision, accelerator, OpenMP, and GPTL details. Explicit runner `-exe` or
`-install_dir` options still override the displayed default.

The Profile tab is a browser interface to `utilities/time_clubb.py`. Its top
benchmark panel configures the case, process/per-process-batch-size sweep,
repetitions, executable, configuration, overrides, and additional
`run_scm.py` arguments. A direct one-second polling path reads the active
summary and process rows and renders figures server-side, so results appear
after each measured repetition while the broker-owned job is running; warmups
remain hidden. Browser stores retain only compact timer/process choices rather
than the growing raw timing table. The running row counter and all four figures
are returned by the same callback response, so visible progress cannot advance
independently of the plots. Stored profiles can be overlaid, compared with a baseline, or
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
summaries in memory instead of storing duplicate summary files. **Export
selected** downloads complete profiles as a ZIP; **Import** accepts those ZIPs
on another machine or checkout. Provenance includes the effective vertical
level count, observed model steps, source revision, executable checksum, host,
timer backend, and time basis so Dash can flag potentially incomparable runs.

## Install

The top-level launcher can create the local virtualenv, install dependencies,
and start the foreground dashboard manager:

```bash
./launch_dashboard.sh
```

The manager starts the runtime broker and Dash as child processes. Arguments
are passed through to `dash_app/app.py`, for example:

```bash
./launch_dashboard.sh --port 23404 -debug
```

If Dash crashes or stops reporting its broker heartbeat, the manager retries it
every 10 seconds for up to 5 minutes. A successful restart is selected without
opening another browser tab. If Dash does not recover within that window, the
manager reports the last failure, gracefully stops broker-owned work, stops the
broker, and exits nonzero. `SIGINT`, `SIGTERM`, and terminal hangup use the same
ordered shutdown.

The broker also watches a private manager heartbeat. If the manager is killed
without a chance to clean up, a replacement launcher can adopt the broker for
30 seconds. After that grace period, the broker stops the orphaned Dash process
group and active Compile/Run/Tune work, then exits.

Dash serializes ordinary callbacks by default to protect NetCDF/HDF5 access.
Explicitly expensive callbacks use isolated background worker processes. Use
`--threaded` only for short diagnostics on a stack known to be thread-safe.

For manual setup, run this from the repository root:

```bash
python3 -m pip install -r dash_app/requirements.txt
```

Run the Dash test suite with the same environment:

```bash
tests/run_pytests.sh -dash
```

Compile CLUBB before using the run tab:

```bash
./compile.py
```

Plotting existing NetCDF output does not require a fresh compile.

Dash builds the Tune controls from a checked-in, Fortran-validated bound table;
it does not need the Python/F2PY interface merely to start. Actual Tune jobs
use CLUBB's in-memory F2PY loss driver, so compile with `-python` before
running a tuning workflow:

```bash
./compile.py -python
```

When Dash is launched with `./launch_dashboard.sh`, build the extension with
the same virtualenv Python that runs Dash. This avoids loading an F2PY module
compiled against NumPy 1.x into a Tune worker's NumPy 2 environment:

```bash
.venv-dash/bin/python compile.py -python
```

Use the corresponding `bin/python` path if `CLUBB_DASH_VENV` names a different
virtualenv.

## Run

From the repository root:

```bash
python3 dash_app/app.py
```
or
```bash
python3 dash_app/app.py &
```

By default the app opens in a browser at port `23404`, or the next available
port. Starting it again while this checkout's dashboard is already running
reopens the registered dashboard instead of starting a second process. Use
`python3 dash_app/app.py --help` for host, port, debug, and threading options.

## JULY_2017 statistics vs. 3-D recreation viewer

For a deliberately small, standalone comparison of horizontally averaged
fields only, run:

```bash
.venv-dash/bin/python dash_app/july_2017_les_comparison.py
```

It overlays the original JULY_2017 SAM profile statistic with the matching
resolved horizontal average recalculated from each 3-D recreation snapshot.
The recreated curve does not include any native SAM SGS contribution.

## Local agent integration

### Runtime boundary

The durable local broker is dashboard runtime infrastructure, not an agent
implementation. Its canonical modules live under `dash_app/shared/`:

- `broker.py`, `gateway.py`, and `broker_protocol.py` own the detached local
  service and its compatibility contract.
- `activity.py` owns durable handoff events and Compile/Run/Tune job state.
- `actions.py` owns the application services and semantic dashboard operations
  used by both native Dash callbacks and external adapters.
- `broker_client.py` is the internal Dash client; it has no agent session or
  model/adapter dependency.

The `dash_app/agent_integration/` package contains the stdio MCP adapter, the
dashboard-owned Streamable HTTP endpoint, a small transient generic client,
and hidden browser-handoff polling. There is no persistent agent session, chat
drawer, agent presence list, or bridge process.
Runtime and application services live under `dash_app/shared/` and are imported
directly by the dashboard and MCP adapter.

Start the launcher first. Its manager starts (or reuses) a small localhost-only
broker sidecar, prints the path to its private connection record, and opens the
dashboard.
Agents use the MCP adapter for one operation at a time; closing the adapter
ends that transient connection. The durable broker continues to own
Compile/Run/Tune/artifact workers and recovery across dashboard or adapter
restarts.

### Add the running dashboard to a Codex chat

The manager-owned broker starts one loopback-only Streamable HTTP MCP endpoint
for the checkout. Open the dashboard's bottom-left utilities menu
and copy the values under **MCP connection** into the chat's manual MCP-server
setup:

1. Choose **Streamable HTTP**.
2. Enter the displayed **Server URL** (it ends in `/mcp`).
3. Supply the displayed **Bearer token** when prompted for authentication.

The displayed instance ID identifies the selected checkout broker.
The endpoint is authenticated with a random per-instance bearer token, checks
the owning broker PID and start time, and routes every typed request to the
broker connection captured at startup. It is valid only while the checkout's
broker is alive; Dash itself may restart around it. Normal
shutdown removes its private record; crash recovery reconciles records by
checking endpoint, dashboard, manager, and broker liveness. The
endpoint is stateless HTTP and does not create a persistent agent chat/session.
No Codex configuration file is edited automatically, and another checkout has
a different URL, token, and instance ID.

The durable broker is intentionally separate from this endpoint. Removing an
MCP endpoint does not stop active broker-owned jobs; a later dashboard or
transient adapter can recover them through the normal typed service boundary.
If it cannot start, the broker's local diagnostic log is
the `broker` entry in the private connection record's `log_paths` object.
The same object lists the Dash application log and MCP endpoint log. When
started through `./launch_dashboard.sh`, Dash stdout and stderr remain live in
the launching terminal and are also captured in the private rotating `app`
log. Logs are mode-0600, capped at 5 MiB with three backups, and contain no
bearer credentials or request payloads.
For an intentional broker-code update (not a normal Dash restart), stop the
idle broker and then relaunch Dash:

```bash
python -m dash_app.shared.broker stop
```
The `connect_to_dashboard` MCP tool is a short-lived authenticated status check;
it does not register an agent or create a session.

Another local adapter can use the small generic client directly:

```bash
python -m dash_app.agent_integration.client connect
```

The generic client also exposes the browser-handoff boundary for adapters that
do not use MCP:

```bash
python -m dash_app.agent_integration.client action inspect_dashboard \
  --payload-json '{"tab":"plots"}'
```

Scientific execution is deliberately not available through that generic
action command. Use the typed MCP tools for Compile, SCM, Tune, artifacts, and
cancellation.

The preferred browser-handoff interface is two generic actions:
`inspect_dashboard` returns every top-level tab's typed operation
manifest plus lightweight live choices, and `invoke_dashboard` invokes one
declared operation with an arguments object. These are browser-handoff
compatibility tools: they can make work visible, but are not the authority for
scientific execution or job state. The generic invoker
therefore accepts only navigation/view handoff; use typed domain tools for
compile, SCM, Tune, artifacts, logs, and cancellation.

The documented MCP interface is purpose-specific: `get_server_info`,
`list_cases`, `submit_compile`, `submit_scm_run`, `submit_scm_batch`, `submit_tune`,
`submit_leaderboard_rerun`, `create_profile_artifact`, `get_job`, `get_run_manifest`, `get_artifact`,
bounded `read_job_log`, and `cancel_job`. Mutating requests require a stable
`request_id`; an identical retry returns the original job, while a changed
request with the same ID is rejected. Typed MCP SCM runs write their scientific
output to the controlled, plot-discoverable
`output/mcp_runs/<batch-id>/` location by default. Each case writes directly beneath
that directory as `<case>_stats.nc`; the broker JobStore records the durable
group and child job/run statuses. Uniquely scoped MCP artifact manifests remain
available through the artifact API, but no manifest is written into the public
scientific-output directory. `submit_scm_run` remains the
backward-compatible one-case wrapper over the same batch service. The MCP
client may optionally supply `out_dir`; it is resolved below the repository's
`output/` directory, for example `dash_default` becomes
`output/dash_default`.

`submit_scm_batch` accepts `{request_id, cases, stats_file, config, overrides,
run_options, max_workers, out_dir}`. `cases` must be a nonempty list of unique checked-in
case names; the other settings are common to every child. The returned parent
`job_id`/`batch_id` can be passed to `get_job`, while each child retains its own
`job_id` and `run_id` for immutable case-level provenance.
Private immutable manifests and temporary execution evidence still live under
ignored, owner-private `output/agent_artifacts/`. This is **ephemeral staging**,
not an experiment or report archive: active bundles are protected from broker
cleanup, completed bundles have bounded retention, and the root may be cleared
between jobs. Copy evidence that must survive into `doc/reports/<report-id>/`
(or a named `output/` directory for a raw run).
Plot artifacts are selected by `run_id`, exact requested/actual time windows,
and coordinate metadata, so a later same-case run cannot silently change them.

`submit_compile` accepts the validated request shape
`{request_id, debug, python_bindings, fresh, gptl}`. All fields except
`request_id` default to `false` (with `debug` defaulting to `true`); setting
`gptl: true` adds the same `-gptl` compile option exposed by the native
Compile-tab checklist. The typed MCP request and the Dash button both pass
this field through the shared broker-owned compile launcher.

The broker is loopback-only. Its connection and job state is held in a private
runtime directory, and every request requires the private bearer token. The
connection record also carries a content fingerprint for the non-test Python
runtime. On a new dashboard start, an idle broker or MCP endpoint with a
different fingerprint is replaced automatically; active jobs are preserved and
the dashboard warns until a later safe restart. This detects uncommitted and
untracked runtime edits as well as committed changes.
For an intentional runtime refresh, use `./launch_dashboard.sh --restart-runtime`
after closing the current dashboard. The command refuses to replace a live
dashboard or interrupt active/queued work; retry after those conditions clear.
gateway accepts only the two safe browser-handoff actions from external MCP
clients; typed domain mutations cross the internal broker boundary. It never
accepts arbitrary shell commands. This private Flask broker is not a remote API:
any future remote deployment must use an authenticated MCP transport plus
explicit authorization.

The older convenience actions (`compile_clubb`, `run_scm`, `plot_profiles`,
`save_profile_png`, `open_dashboard`, `open_note`, `launch_tuning`,
`inspect_tuning`, `run_tuning_loss`, `stop_tuning`, `stop_compile`, and
`stop_run`) remain internal Dash wrappers while callbacks migrate. They are not
accepted as external broker actions; adapters use the typed MCP service.
The typed equivalents carry explicit configuration and time-window choices
rather than inheriting the visible browser state.  They
return stable job/run IDs, so a running request can be inspected or cancelled
without relying on a current tab selection.

The Plot-tab `set_view` operation also accepts an optional
`benchmark_sources` array, for example `{"benchmark_sources":["sam"]}`.
This selects exactly the available SAM overlay for the requested case and is
particularly useful with `run_id`, which keeps the profile output tied to one
immutable SCM run.  Omit the field to preserve the normal/default Plot-tab
selection; an unavailable source is rejected.  The native SAM/COAMPS toggle
buttons and this MCP operation share the same UI-neutral source validator, so
neither path simulates clicks or bypasses case availability checks.

The Plot budget family has a separate typed operation so its controls do not
share a giant all-purpose schema. Use `plots.add_budget` through
`invoke_dashboard`, for example:

```json
{"case":"arm","budget_group":"wp2"}
```

The direct MCP adapter exposes the same request as `add_budget_plot` with
`case`, optional `budget_group` (default `wp2`), optional immutable `run_id`,
or an optional validated `output_dir` below the repository `output/` root,
and the common `time_start_seconds`, `average_minutes`, or `window_preset`
controls. Prefer `run_id` for immutable provenance; when both selectors are
provided they must resolve to the same directory. With neither selector, the
legacy top-level `output/` default is retained. The group is validated against
the selected output before a handoff is published; the native Add budget plot button and this MCP
operation both use the shared Plot state-transition service. WP2 renders its
registered `wp2_*` budget terms through the existing budget plugin. Budget
plots remain single-column and do not introduce a benchmark-overlay control.
Use the typed `plots.list` operation (or direct MCP `list_plots`) to inspect
the currently mounted cards. It returns stable card IDs plus the plot family
and selection. Remove one with `plots.remove`, for example
`{"plot_id":4}`, or direct MCP `remove_plot(plot_id=4)`. IDs are validated
against the current dashboard-owned Plot state before the handoff and again by
the native removal callback; removed IDs are not reused during that dashboard
session. There is intentionally no generic Plot update operation yet.
All typed mutations cross the private internal broker boundary before launch.
The manager-owned broker—not the agent-owned stdio adapter—therefore owns process
watchers and terminal job updates. Closing or replacing an MCP adapter does not
orphan active work; another adapter can reconnect, query the same `job_id`, and
cancel it.

The ordinary **Start tuning** button uses that same broker path; it is not only
for agent-launched jobs. On startup the broker also scans `output/tuner/` (and
the legacy `output_tuner/` location) for
the newest still-running file-backed Tune worker, so a worker that outlived an
older dashboard reload can be reattached and have its keepalive renewed.

For terminal/IDE MCP hosts that are configured from the repository root, the
stdio adapter remains available as a static fallback:

```bash
codex mcp add clubb-dash -- .venv-dash/bin/python dash_app/agent_integration/mcp_server.py
```

### ADG1 two-Gaussian explorer

The Tutorial tab includes the active ADG1 two-Gaussian explorer. It visualizes
the normalized ADG1 diagnosis and a direct-control trivariate comparison using
the same grid moments. It is a teaching visualization, not a replacement for
the full Fortran PDF diagnosis.

## Basic Workflows

- **Run tab:** choose benchmark cases and settings, launch CLUBB, and watch the
  run output in the browser.
- **Plots tab:** load one or more CLUBB output directories and make profile,
  time-height, time-series, budget, and subcolumn plots from the NetCDF files.
- **Tune tab:** configure and monitor tuner runs when the branch and local build
  support the tuner workflow. This requires a CLUBB build compiled with
  `./compile.py -python`.
- **Tutorial tab:** follow short interactive explanations of CLUBB concepts.
  The welcome page suggests a path through the lessons; a vertical page rail
  opens each lesson. The CLUBB Equations page provides a clickable quick
  reference for the core prognostic budgets, PDF transport closures, and
  cloud/buoyancy diagnostics. Its colors preserve the official equation
  document's ownership convention, while a stable inspector explains each
  term's physical role, source, closure path, and implementation relevance.
  The next lesson opens the ADG1 two-Gaussian explorer, where shared moments
  and the supplied moments show how ADG1 component placement and covariance
  allocation change the PDF geometry.
- **Reports tab:** browse immutable, static investigation bundles from
  `doc/reports/`. Each bundle carries its own HTML, figures, excerpts, data,
  and provenance. Dash polls the published JSON catalog, so an agent can add a
  completed report without editing dashboard source or restarting the app.
- **Misc tab:** browse living investigations and focused diagnostics from a
  persistent left-side vertical directory, including the SAM w–rₜ time-height
  neighborhood and the Mixing Length Trajectories explorer. The neighborhood
  browses a pre-rendered 5×5 atlas generated by
  `python -m dash_app.misc_tab.sam_w_rt_neighborhood.atlas`; the trajectory
  explorer reconstructs upward and downward parcel-energy paths from a
  compatible CLUBB statistics file and compares them with its stored
  `Lscale_up`, `Lscale_down`, and `Lscale` profiles.

The equation-guide content is curated in
`dash_app/tutorial_tab/clubb_equations_demo/` from
`doc/CLUBBeqns.tex` and the corresponding current routines in
`src/CLUBB_core/`. It deliberately emphasizes the continuous closed equations;
an on-page caveat distinguishes them from implicit discretization, host
coupling, surface forcing, clipping, limiters, and other enabled adjustments.

## LES Benchmark Overlays

The plots tab can overlay LES benchmark data for cases that define SAM or
COAMPS benchmark files in:

```text
postprocessing/pyplotgen/config/Case_definitions.py
```

Currently those case definitions use the repository benchmark symlink:

```text
input/les_and_clubb_benchmark_runs/
```

`input/les_and_clubb_benchmark_runs` defaults to
`/home/pub/les_and_clubb_benchmark_runs`. The app only shows LES overlay
options when the configured files exist on the local machine. If the archive is
mounted somewhere else, retarget that symlink locally.

## Shared UI components

Reusable themed overlay notecards live in `dash_app/shared/notecard.py`, with
their common styles in `dash_app/assets/05_shared_modal.css`. `notecard`
supports small, medium, large, and full-window panels with arbitrary Dash
content. Plot-card help and tutorial explanations use the same component as the
Compile tab source-check log; plot-family help text is centralized in
`dash_app/plot_tab/plot_types/help_content.py`.
