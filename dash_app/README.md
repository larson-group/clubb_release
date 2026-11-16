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
   Fortran or JAX.
2. Prepare the selected implementation: Fortran needs a build from the
   **Compile** tab or `./compile.py`. JAX sets up its own environment on
   first use and needs no Fortran build.
3. Choose a case and output directory, launch the run, and watch its log.
4. Open **Plots** to view the output. You can also plot existing NetCDF files
   without compiling or running anything.

### Creating and editing cases

The Run tab's **Create Case** button in the **Custom cases** row opens a model-settings editor.
Choose **Defaults**, a catalog case or a saved case in **Start from**, then
enter a case name. Edit the name/value rows under each namelist heading;
booleans use a **false | default | true** slider when a catalog default exists,
and **false | unset | true** otherwise. Orange markers show changes from the starting point.
Blank cells restore catalog defaults or omit fields without a catalog default;
enter `""` for an explicit empty string and JSON notation for arrays.
**Default** restores the catalog's explicit value.
**Unset** omits an optional boolean that has no catalog default and uses the model
default; **false** and **true** explicitly set it. Click or drag between the three
positions, or use the arrow keys while the slider is focused.
The optional **Namelist preview → Update preview** displays the generated input.
The editor/converter does not check model validity or whether input files exist;
model validation happens downstream when you run.
**Close** keeps the current draft ready for reopening. Selecting a different
starting point replaces it; reloading the page clears unsaved edits.

**Save case** writes overrides relative to the catalog defaults in
`input/<name>.case.json`, an ignored local file, and selects it for the next run.
Saved copies inherit later changes to shared defaults. They appear in the **Custom cases** row below the
catalog buttons. To replace a catalog or saved case, use its existing name
and confirm the before/after values in the replacement popup. Catalog replacement
updates that case's overrides; other cases and shared defaults are preserved.

The saved name identifies the run and output; **runtype** identifies its physical
inputs and case-specific behavior. A BOMEX copy can have a new name while keeping
`runtype = bomex`. Defaults can be saved with incomplete or invalid model settings.
For a runnable case, supply its physical runtype, inputs and timing. LES
comparison defaults are separate from diagnostic `setfields.datafile` inputs;
copies retain their source comparison defaults while the physical runtype stays
the same. Copies with a changed runtype require comparison configuration before
they can be tuned against a benchmark.

Normal runs resolve cases solely from saved custom JSON or the
catalog. Dash edits and saves JSON; native namelists are generated for each run. Case settings are resolved when a queued
run starts, so edits made while waiting affect that run. The generated namelist
is retained with its output. Reload other already-open workflow tabs after
creating cases to refresh their case dropdowns.

### Choosing statistics

The Run tab offers **All**, **Standard**, and **Core** from
[`input/stats.jsonc`](../input/stats.jsonc). **Custom…** opens its
checkbox tree, with all branches initially collapsed. Expand a category to
choose individual variables. Each compact row shows the variable, units, grid
(thermo, momentum or surface), and long name; long rows scroll horizontally.
Selecting a parent
includes every variable below it, and mixed parents show a dash. The tall,
scrollable editor previews the selected registry definitions.
Search refines the tree as you type, matching names, descriptions, units, grids
and category paths (case insensitive; multiple words must all match). Matching
branches open automatically. **Clear search** restores the previous expanded
branches. Filtering preserves selections; category checkboxes still select the
whole category, including hidden variables. Each popup opening starts with an
empty search.
Draft checkbox changes, previews and popup reopening run locally in the browser;
opening and closing branches uses native disclosures. Apply and Save use the
shared Python stats helpers.
**Clear all**, then **Use selection**, disables statistics.

**Saved lists & presets** selects a named JSON preset or saved snapshot. In the custom editor,
**Save stats config** writes a new variable-name list to `input/<name>.stats.json`
and adds it to that list. These local snapshots are ignored by Git.
Names may contain letters, numbers, underscores and
hyphens; existing files are never overwritten. Snapshots retain their exact
variables even if the catalog changes, and can also be passed to `-stats` on
the command line. Reload the dashboard after editing the JSON hierarchy or
changing a saved/named preset outside the app.

MCP `stats_file` accepts category expressions, individual `var:NAME` selectors
(including qualified paths and comma/plus unions), `none`, or a relative
saved/named preset name. MCP and the Run tab default to the **Standard** category.
Presets select only their own fields; combining a preset and category selects
their union, with each variable emitted once. Shared `-config` settings do not
change the selected stats.

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

Fortran uses the selected compiled build, which you can change or
rebuild through the chooser. New Tune jobs use the selected Fortran or JAX loss
backend. Each saved revision keeps its runtime choice; Continue and result
reruns use that saved choice even if the chooser changes later.

## Basic Workflows

- **Run tab:** choose benchmark cases and settings, launch CLUBB, and watch the
  run output in the browser.
- **Profile tab:** measure runtime across process counts and batch sizes, then
  compare saved profiles. See [profiling](#profiling) for details.
- **Plots tab:** load one or more CLUBB output directories and make profile,
  time-height, time-series, budget, and subcolumn plots from the NetCDF files.
- **Tune tab:** configure and monitor tuning runs. Choose JAX CPU/GPU in the
  runtime chooser, or use a [Python build](#python-runs-and-tuning) for the
  Fortran loss backend. Taylor diagrams, parameter spread and the parameter
  landscape are always shown and update automatically as new results arrive.
  Correlations, sensitivity and interactions are optional under **Additional
  plots**. Plot calculations run in the background, separately from live status,
  uptime and best loss.
- **Tutorial tab:** explore CLUBB concepts through interactive lessons,
  including a guide to the model equations and the ADG1 two-Gaussian explorer.
- **Documents tab:** browse PDF, Markdown and plain-text (`.txt`) files under [`doc/`](../doc/),
  including `README.md` files. Search by filename or relative path.
- **Misc tab:** open focused diagnostics such as the SAM w–rₜ neighborhood
  viewer, Mixing Length Trajectories explorer, and
  [Affine Coherence Mixing Length](#examine-coherent-mixing-lengths) diagnostic.
  Setup and implementation
  notes are in [DEVELOPMENT.md](./DEVELOPMENT.md#misc-subtabs).

### Examine coherent mixing lengths

Open **Misc → Mixing length → Affine Coherence Mixing Length** to see how
parcel launch properties and energy barriers affect the mixing distance.
The [diagnostic](misc_tab/affine_mixing_length_trajectories/) recomputes lengths
from one saved column and compares them with saved CLUBB profiles. The
[method description](../doc/coherence_mixing_length.md) explains coherence as
the fraction of assumed launch energies able to reach a level.

To create suitable input, first build Fortran as described in the quick start,
then run from the repository root:

```bash
./run_scripts/run_scm.py -stats all arm
```

Launch Dash, select that run's **Output directory** and **Case**, then choose
the grid column and record. **Add file by path → Load** opens a compatible
NetCDF file outside the discovered output folders. **Refresh** rescans files.
The loader requires the thermodynamic/PDF fields and saved directional lengths
listed in [`REQUIRED_VARIABLES`](misc_tab/mixing_length_trajectories/analysis.py);
a reduced statistics list may omit them.

The plots show parcel energy/buoyancy, the greatest energy deficit encountered
along each path, survival weights, directional/final lengths, and changes
between nearby records. In each matrix, launch altitude is horizontal and
target altitude vertical. A weight is statistical survival, not retained
parcel material. The energy display continues past deterministic stopping
points so adverse work remains visible.

The four parcel controls initialize from saved parameters when available.
The **Overwrite method** control changes how the recomputed trajectories
are combined: keep local reaches, extend maximum endpoints, or average by
arrival coherence. This is separate from the model's `Lscale_method` flag.
Changing **Maximum calculated altitude** creates a shorter calculation with
a new upper boundary; it does not just zoom the figure.

Controls update the diagnostic, leaving the saved curves and model files as
recorded. Saved statistics can be time averages, and this reader estimates the
surface floor using a fixed 40 m grid scale. The comparison therefore helps
explore behavior; it does not establish numerical parity with a model timestep.

## Advanced usage

### Python runs and tuning

Python SCM runs and Fortran-backed Tune jobs require CLUBB's Python/F2PY
interface. JAX tuning uses its managed environment and needs no Fortran build.
Build the Fortran interface with:

```bash
./compile.py -python
```

The build prepares and uses the shared Python environment automatically, so
the compiled interface matches Dash's NumPy version. The interface is not
needed just to open the dashboard or configure Tune controls.

### Launch options

Pass application options through the launcher, for example:

```bash
./launch_dashboard.sh --port 23404 -debug
```

The dashboard launcher installs Dash packages into that environment with uv on
first launch. Use `./launch_dashboard.sh --help` for host, port, debug, and
threading options.
Dash serializes ordinary callbacks to protect NetCDF/HDF5 access; use
`--threaded` only for diagnostics on a stack known to be thread-safe.
`-debug` enables developer tools and reloading without changing that threading
default. Leave it off when comparing interaction performance.

Pages load on their first visit and retain their state when you switch tabs.
Tutorial lessons and Misc tools initialize individually; Compile, Run, Profile,
and Tune initialize independently. Their shared runner chooser stays available
across tabs, opens immediately, and checks builds and JAX readiness in the
background once per page session. Compile and the chooser share one build
catalog; compilation and deletion update it. Opening the chooser or changing
CPU/GPU preferences reuses that state. Selection paints from cached panels
without waiting for dashboard callbacks. **Refresh runners** rescans changes made
outside Dash. Runtime setup updates readiness without rediscovering devices.
Unvisited pages do not scan outputs or generate figures during startup. Saved
tab and control selections still restore on refresh.

The launcher supervises Dash and attempts recovery after a crash. Lifecycle
and broker details are in the
[development notes](./DEVELOPMENT.md#local-mcp-endpoint-lifecycle).

### Profiling

Use **Profile** to choose a case, process counts, batch sizes, and repetitions.
The batch-size sweep defaults to `1,2,3,4,5,6,7,8,12,16,32,64` columns per
process for CPU targets and `256,512,1024,2048,4096` for GPU targets. Switching
between CPU and GPU resets that field; switching within the same device type
preserves your edits. JAX uses its selected CPU/GPU target. Compiled builds
currently use their OpenACC/OpenMP GPU build setting as a heuristic, even if
that build is run on a CPU. An explicit executable override does not change
this chooser-based default. Process counts are a separate sweep.

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

### Documents

Save finished reports as PDF, Markdown or plain text (`.txt`) under [`doc/`](../doc/), for example
[`sam-heatmap-color-signals.pdf`](../doc/sam-heatmap-color-signals.pdf). The
Documents tab discovers them automatically; search matches filenames and paths
relative to `doc/`. Files are grouped as PDFs, Markdown, then Text, with paths
sorted alphabetically within each group. Markdown includes `README.md` files.
Hidden paths and other formats are excluded. Plain text preserves spacing and line breaks without interpreting
Markdown or HTML. Text files use UTF-8. No index, manifest or publisher is required.

Finish the research before writing the document. Embed figures in a PDF, or keep
Markdown images beside the document. Raw runs, logs and other generated research
artifacts belong in a named `output/` directory.

Relative links in Markdown open their target file in a new browser tab, resolved
from the Markdown file's folder. Linked source files such as `.F90` and `.py`
open as plain text; PDFs and images use the browser viewer. Files outside `doc/`
must be tracked by Git. External links keep their original destinations.

To open a document through the semantic browser interface, use
`reports/open_document` with `{"path": "sam-heatmap-color-signals.pdf"}`.

### Plot loading

Enable **SAM LES** or **COAMPS LES** to compare against that case's configured
benchmark data. Enabling a source starts the shared
[LES finder/downloader](../input/scripts/confirm_les_data.py) in the background,
even if no plots are open yet. The button shows a spinner and **Loading…** while
its files are being found or downloaded. Model plots and controls remain usable;
the overlay appears automatically when ready. A failed download shows an error
beside the button; toggle the source off and on to retry. Turning it off or
switching cases prevents a late result from adding an unwanted overlay. An
already running download may finish and remain cached for later use.

Shared files and local downloads in `input/les_and_clubb_benchmark_runs` are reused.
Opening a case with LES sources disabled or listing outputs does not download
data. Custom expressions referencing `les.<field>` also prepare their LES source
in the background, then update automatically.


Time-height comparisons stack one panel per selected output directory, followed
by any enabled LES benchmark panels. Panels share axes and a color scale and
use the selected averaging interval. In comparison mode, choose a variable
with time and height dimensions in every selected output.

The output chooser keeps **Available** and **Active** folders together inside
one box, separated by a vertical divider. Click an available folder to move it
into the active side as a bubble; remove a bubble to return the folder to the
available side. Active bubbles flow onto additional lines. The compact view has
a fixed height and shows five available entries. Click the **Output folders ▾**
heading to expand both sides downward over the content below, without resizing
the header or plot area; click again to collapse. The expanded view includes all
available entries and scrolls for longer lists.

The chooser discovers folders when it first loads or when **Refresh outputs**
is clicked. **Refresh outputs** and **Add extra folder** sit in the chooser's
top bar. Adding an extra folder opens a path field in the expanded overlay and
accepts paths outside the discovered list.
Moving folders updates the selection immediately; metadata and figures load in background workers.
Existing figures stay visible with an updating label until replacements are
ready. Rapid changes apply only the latest selection. The selected case, plot
cards, and view controls survive a browser refresh; older saved workspaces
migrate automatically.

### PyPlotGen export

The **PyPlotGen** section in the top bar creates a static plot gallery from the
active output folders. Search the options selector, press **Enter** to add an
option, and repeat to select more. For example, search for `pdf` to select
`-pdf`, or add `-cases` and enter `bomex arm` in its value field to limit the
gallery to those cases. Options that take values show an additional field;
quote paths containing spaces. `-les` is selected initially and can be removed.

Click **Generate PyPlotGen** to start; the same button cancels an active export.
The completed gallery opens in a new tab and is stored under `output/pyplotgen/`.
Click **(?)** to see what PyPlotGen does and its complete command-line help.
The option list and help come from
[`pyplotgen.py -h`](../pyplotgen/pyplotgen.py), including requirements and option
restrictions. Dash supplies `-clubb` and `-output_dir` from the active folders
and export destination, so those and the help flags are excluded from selection.

### JULY_2017 statistics vs. 3-D recreation viewer

For a deliberately small comparison of horizontally averaged fields only,
use the JULY_2017 viewer in Dash.

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
run_options, max_workers, out_dir}`. `cases` must be a nonempty list of unique
checked-in or locally saved case names; the other settings are common to every
child. The returned parent
`job_id`/`batch_id` can be passed to `get_job`, while each child retains its own
`job_id` and `run_id` for recorded case-level provenance.
Private immutable manifests and temporary execution evidence still live under
ignored, owner-private `output/agent_artifacts/`. This is **ephemeral staging**,
not an experiment or report archive: active bundles are protected from broker
cleanup, completed bundles have bounded retention, and the root may be cleared
between jobs. Save final reports as PDF, Markdown or plain text under `doc/`. Keep raw runs and other
generated research artifacts in a named `output/` directory.
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
for agent-launched jobs. On startup the broker also scans `output/tuner/` for
the newest still-running file-backed Tune worker, so a worker that outlived an
older dashboard reload can be reattached and have its keepalive renewed.

For terminal/IDE MCP hosts that are configured from the repository root, the
stdio adapter remains available as a static fallback:

```bash
codex mcp add clubb-dash -- python3 dash_app/agent_integration/mcp_server.py
```

### LES Benchmark Overlays

The plots tab can overlay LES benchmark data for cases that define SAM or
COAMPS benchmark files in:

```text
pyplotgen/config/Case_definitions.py
```

The app offers the sources configured for each case, even before their files
are available locally. Enable **SAM** or **COAMPS** to include that source in
the plots. A background worker uses the shared
[LES finder/downloader](../input/scripts/confirm_les_data.py) to reuse existing
files or download missing ones; an input symlink is optional. See
[plot loading](#plot-loading) for progress and retry behavior, and
[acquiring LES benchmark data](../input/README.md#acquiring-les-benchmark-data)
for search locations, download placement and standalone use.

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
bash human_reviewed_tests/run_pytests.sh -dash
```

This includes the bounded MCP stdio protocol regression; it does not run SCM
cases or launch the full dashboard workflow.

For the manual parcel-trajectory comparison, supply fresh Fortran ARM
`all` output from 625 timesteps, including saved record 612:

```bash
python3 dash_app/auto_llm_generated_tests/integration_tests/run_dash_mixing_length_test.py \
  -stats_file output/tests/dash_arm/arm_stats.nc
```

The check compares the mixing-length profiles, parcel paths and figure data.
It does not generate the reference run, and missing input or a failed numerical
assertion fails the check. It is separate from the `clubb_dash` pytest job.

See [DEVELOPMENT.md](./DEVELOPMENT.md) for UI conventions, service boundaries,
and how runtime selections are recorded in jobs.
