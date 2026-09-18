# CLUBB Dash App

The Dash app is a browser interface for common CLUBB workflows. It provides a
Run tab for launching CLUBB cases, a Profile tab for configurable process-based
timing sweeps, and a Plots tab for inspecting CLUBB NetCDF output. The app is
intended to be run from an existing CLUBB checkout.

## Quick start

From the repository root:

```bash
./launch_dashboard.sh
```

The launcher installs the dashboard dependencies and opens the app in your
browser, normally at port `23404`. Keep the launching terminal open while using
the app. Running the command again reopens this checkout's existing dashboard.

To run a case:

1. Open **Run** and click the runtime badge beside the launch button to choose
   Fortran, Python, or JAX.
2. Prepare the selected implementation: Fortran needs a build from the
   **Compile** tab or `./compile.py`; Python needs a
   [Python build](#python-runs-and-tuning). JAX sets up its own environment on
   first use and needs no Fortran build.
3. Choose a case and output directory, launch the run, and watch its log.
4. Open **Plots** to view the output. You can also plot existing NetCDF files
   without compiling or running anything.

### Choosing CPU or GPU for JAX

In the runtime chooser, select JAX and then the CPU or a compatible GPU tile.
NVIDIA tiles show the physical GPU index, model, and memory. The chooser
remembers your preference in this browser. The first run may install the
selected runtime; a missing or incompatible GPU produces an error.

A selected GPU is used for all cases in that job; selecting it does not spread
the cases across multiple GPUs. Leave **Preallocate GPU memory** off for a
shared NVIDIA GPU. Enable it when you want up-front memory reservation; it is
unavailable for CPU and Metal. The chooser's **?** button explains the options.
See the [JAX guide](../clubb_jax/README.md#gpu-running) for
hardware requirements and command-line equivalents.

Fortran and Python use the selected compiled build, which you can change or
rebuild through the chooser. Tune uses its own Fortran/Python worker regardless
of the JAX selection.

## Basic Workflows

- **Run tab:** choose benchmark cases and settings, launch CLUBB, and watch the
  run output in the browser.
- **Profile tab:** measure runtime across process counts and batch sizes, then
  compare saved profiles. See [profiling](#profiling) for details.
- **Plots tab:** load one or more CLUBB output directories and make profile,
  time-height, time-series, budget, and subcolumn plots from the NetCDF files.
- **Tune tab:** configure and monitor tuning runs. This requires a
  [Python build](#python-runs-and-tuning).
- **Tutorial tab:** explore CLUBB concepts through interactive lessons,
  including a guide to the model equations and the ADG1 two-Gaussian explorer.
- **Reports tab:** browse saved investigation reports from `doc/reports/`,
  including their figures, data, and provenance.
- **Misc tab:** open focused diagnostics such as the SAM w–rₜ neighborhood
  viewer and Mixing Length Trajectories explorer. Setup and implementation
  notes are in [DEVELOPMENT.md](./DEVELOPMENT.md#misc-subtabs).

## Advanced usage

### Python runs and tuning

Python runs and Tune jobs require CLUBB's Python/F2PY interface. After the
launcher has prepared the Dash environment, build it with that same Python:

```bash
.venv-dash/bin/python compile.py -python
```

Use the corresponding `bin/python` path if `CLUBB_DASH_VENV` names a different
virtual environment. Using the same environment avoids NumPy compatibility
problems when loading the compiled interface. The interface is not needed just
to open the dashboard or configure Tune controls.

### Launch options and manual setup

Pass application options through the launcher, for example:

```bash
./launch_dashboard.sh --port 23404 -debug
```

For manual setup in your chosen Python environment:

```bash
python3 -m pip install -r dash_app/requirements.txt
python3 dash_app/app.py
```

Use `python3 dash_app/app.py --help` for host, port, debug, and threading options.
Dash serializes ordinary callbacks to protect NetCDF/HDF5 access; use
`--threaded` only for diagnostics on a stack known to be thread-safe.
`-debug` enables developer tools and reloading without changing that threading
default. Leave it off when comparing interaction performance.

Pages load on their first visit and retain their state when you switch tabs.
Tutorial lessons and Misc tools initialize individually; Compile, Run, Profile,
and Tune initialize together because they share the build selector. Unvisited
pages do not scan outputs or generate figures during startup. Saved tab and
control selections still restore on refresh.

The launcher supervises Dash and attempts recovery after a crash. Lifecycle
and broker details are in the
[development notes](./DEVELOPMENT.md#local-mcp-endpoint-lifecycle).

### Profiling

Use **Profile** to choose a case, process counts, batch sizes, and repetitions.
Results appear after each measured repetition; warmups are excluded from the
default plots. Saved profiles can be overlaid, compared with a baseline, or
viewed as process distributions and exclusive-cost decompositions.
JAX does not yet emit the native timer files required by Profile.

Reusing a profile name asks for confirmation before replacing the existing
profile. **Export selected** downloads complete profiles as a ZIP; **Import**
loads those ZIPs on another machine or checkout. Profiles include run and build
metadata so Dash can flag potentially incomparable results. Storage and update
details are documented in the
[development notes](./DEVELOPMENT.md#profile-results).

### Plot loading

The output chooser discovers folders when Plots first loads, when opened, or
when **Refresh outputs** is clicked. Adding and removing folders updates the
selection immediately; metadata and figures load in background workers.
Existing figures stay visible with an updating label until replacements are
ready. Rapid changes apply only the latest selection. The selected case, plot
cards, and view controls survive a browser refresh; older saved workspaces
migrate automatically.

### JULY_2017 statistics vs. 3-D recreation viewer

For a deliberately small, standalone comparison of horizontally averaged
fields only, run:

```bash
.venv-dash/bin/python dash_app/july_2017_les_comparison.py
```

It overlays the original JULY_2017 SAM profile statistic with the matching
resolved horizontal average recalculated from each 3-D recreation snapshot.
The recreated curve does not include any native SAM SGS contribution.

### Local agent integration

#### Runtime boundary

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

#### Add the running dashboard to a Codex chat

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

### LES Benchmark Overlays

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

### Shared UI components

Reusable themed overlay notecards live in `dash_app/shared/notecard.py`, with
their common styles in `dash_app/assets/05_shared_modal.css`. `notecard`
supports small, medium, large, and full-window panels with arbitrary Dash
content. Plot-card help and tutorial explanations use the same component as the
Compile tab source-check log; plot-family help text is centralized in
`dash_app/plot_tab/plot_types/help_content.py`.

### Development and tests

Run the Dash test suite with the dashboard environment:

```bash
tests/run_pytests.sh -dash
```

See [DEVELOPMENT.md](./DEVELOPMENT.md) for UI conventions, service boundaries,
and how runtime selections are recorded in jobs.
