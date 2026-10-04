# CLUBB Dash-first workflow

Use this guidance when a request can be completed, demonstrated, or reviewed
through the local CLUBB Dash app.  Treat Dash as the preferred shared surface:
the user should be able to see the selected tab, inputs, job state, console
output, and plots while the work is happening.  Do not replace a suitable Dash
interaction with an invisible one-off shell workflow merely because it is
familiar.

This is workflow guidance, not permission to start processes, launch large
jobs, change model inputs, or overwrite user work beyond what the user asked
for.

## 1. Decide whether Dash is the right surface

Use Dash first for requests involving:

- compilation, SCM runs, output/profile plots, PDF contours, and saved plot
  images;
- Tune configuration, launch, status, stopping, or inspecting parameter and
  loss settings;
- opening Tutorial, Misc, Reports, or another dashboard page;
- an investigation whose useful evidence can be displayed in the app;
- a request to create or inspect a static report.

For a report, **always read `doc/reports/README.md` before authoring it**.
Follow its bundle contract and publication procedure.  A completed report uses
embedded, committed rendered figures plus the reproduction script and compact
provenance data.  After publishing, open the exact report through the active
dashboard using its `reports/open_report` semantic operation.

Do not force Dash on a purely conceptual question, a narrow source-code edit
with no useful interactive view, or a task the dashboard does not expose. A request to run a named Jenkins job
belongs to the [Jenkins workflow](jenkins_workflow.md). Say
briefly why the normal local workflow is necessary in those cases.

## 2. Discover before connecting or starting

First make a read-only, checkout-aware discovery pass.  Look for a live Dash
browser/client and its matching local agent broker, not merely a stale process
or connection file.  This checkout records a recent browser heartbeat under
`/tmp/clubb_dash_client_<checkout-hash>.heartbeat`; the dashboard broker is
durable across normal Dash reloads, so a live broker alone is **not** proof
that a browser view is open.

Use the available local Dash connection capability (normally
`connect_to_dashboard`) to verify the matching broker and then
`get_server_info` / `list_cases` to learn the documented local service.  Use
purpose-specific MCP operations (`submit_compile`, `submit_scm_run`, `submit_scm_batch`,
`submit_tune`, `create_profile_artifact`, status/log/artifact reads, and
`cancel_job`) for scientific work.  Keep the same `request_id` when retrying a
mutation.  `inspect_dashboard` / `invoke_dashboard` are deprecated
compatibility-only browser handoff tools: use them to make a result visible,
not as the authority for execution, job state, or a scientific selection.
Never invent component IDs, browser clicks, raw shell payloads, or operations
outside the declared service.

Treat an MCP **"Transport closed"** error as an adapter-launch failure, *not*
as evidence that Dash is absent.  In particular, a Dash debug/source reload
does not normally stop the manager-owned broker.  Do not launch a second dashboard
in response to that error.  Refresh or restart the local agent's MCP adapter,
then retry discovery; if the problem persists, report that the adapter itself
needs attention separately from the dashboard/broker state.

Be explicit about ambiguity:

- **Exactly one live dashboard/browser for this checkout:** connect to it and
  use it.
- **More than one plausible live dashboard or checkout:** report the choices
  (at least checkout/path, port or PID if known, and whether a browser is
  active) and ask the user which one to use.  Do not connect or start a
  replacement by guesswork.
- **No confirmed open dashboard:** tell the user that none was found and ask
  permission before launching one.  Do not silently start your own Dash app.
  State the proposed command/configuration and why Dash would help.

Only after the user declines that launch permission may you fall back to a
non-Dash workflow for a request that otherwise belongs in Dash.  Make the
fallback visible in your reply and preserve useful artifacts so they can be
opened later.  If the user permits a launch, use the repository launcher
(`./launch_dashboard.sh`) rather than recreating its environment or starting
an ad-hoc server.

## 3. Connect through the current broker

Read the current `dash_app/README.md` and `dash_app/DEVELOPMENT.md` connection
sections when connection/lifecycle details matter. There is no persistent agent
session, chat drawer, presence list or bridge process. Do not recreate or wait
for those retired facilities. Report connection/work status in the current
conversation.

The broker owns a stable loopback Streamable HTTP MCP endpoint for this checkout.
Use the available purpose-specific MCP adapter, or the endpoint shown in the
dashboard's **MCP connection** utility. Never copy its bearer credential into a
worklog, report or commit. The agent connection is transient; the broker owns
Compile/Run/Tune/artifact workers and recovery across adapter/Dash reloads.
Ending an agent turn is not proof that broker-owned work stopped, and a live
broker is not proof that a browser rendered an artifact. Inspect durable job
state through the typed service and the matching browser for visible results.

## 4. Operate semantically and make work visible

Before a new class of action, inspect the relevant typed service schema and
use its structured request.  Then use the compatibility handoff only when it
will make the selected native tab, controls, log, or artifact visible.

- For compile/run/plot/tune work, submit one structured job, retain its
  `job_id` / `run_id`, and then open the appropriate native view.  Do not
  separately run a hidden duplicate command just to obtain output.
- Treat agent-visible configuration as the source of truth.  Use structured
  arguments for cases, fields, overrides, time windows, and parameters; never
  smuggle arbitrary shell text through the dashboard boundary.
- Before `submit_tune`, call `list_tunable_parameters(config)` and use the
  returned names and suggested bounds.  Do not reuse a parameter name from
  chat history, a saved result, or a different configuration: the selected
  namelist is authoritative.
- For a running Tune request, inspect the displayed request/state before
  acting.  Stop it before asking Dash to change its immutable settings; then
  create a replacement request through normal controls.
- After work completes, open the requested plot, tab, or report and reply with
  the outcome plus any important limitation.  The dashboard UI action is part
  of completion, not an optional afterthought.

For multi-case SCM work, use `submit_scm_batch` and the shared Run service.
Follow the README's current output contract: one flat, public, plot-discoverable
`output/mcp_runs/<batch-id>/` group by default, with `<case>_stats.nc` directly
inside it. Keep private manifests in their existing store. Do not invent a
second output hierarchy or put long-lived scientific results in ephemeral
agent staging. Distinguish the CLI `-output_dir` convention from the current
structured MCP schema; use the schema's actual field names.

A successful handoff/request is not proof of a visible plot. Inspect the active
browser's selected tab, figure and relevant settings before saying it was
shown; rapid updates may replace an earlier requested view. If browser checks
are unavailable, report scientific job/artifact completion separately from
unconfirmed presentation. Avoid launching another job to repair a view. Resolve conflicting status from
the current job/request ID, durable state and fresh logs before calling work
hung or complete; an old record or quiet browser is insufficient evidence.

For Dash code changes, keep model/job behavior in the existing UI-neutral
services and shared utilities documented in `dash_app/DEVELOPMENT.md`. Browser
callbacks and MCP adapters should call the same owner, not grow separate logic.
Store user intent and clear stale derived/override state when its selection
changes. Preserve current cancellation/admission semantics in the broker/JobStore
rather than adding a UI-only runner or another public manifest. When changing
profiling/polling, check the actual live graphs as well as retained compact
timing data and failure logs; status updates alone do not establish a live view.

## 5. Reports: publish atomically, then show them

For static reports, `doc/reports/README.md` is authoritative.  In particular:

1. Complete the investigation and render final PNG/SVG evidence into the
   report bundle; do not make a report whose figures are generated on viewer
   load.
2. Use `ReportBuilder` to stage and publish a complete bundle atomically.
   Update an existing investigation using the same report ID with
   `replace=True`; choose a new ID only for a genuinely separate report.
3. Keep the reproduction script and compact metadata beside the embedded
   evidence.
4. After a successful publish, use the active dashboard's `reports` tab
   `open_report` operation for the exact report ID.  This forces the browser
   iframe to reload even if that report was already selected.

`output/agent_artifacts/` is only ephemeral staging for active agent/MCP work.
It can be cleared after jobs finish, so never leave a report linked to or
dependent on it. Copy final figures, compact data, snippets, and reproduction
scripts into `doc/reports/<report-id>/` before publication. Use a named
`output/` directory, not agent staging, for a raw run that must be retained.

If no dashboard is active and the user declined permission to launch one,
still publish a valid report when requested, but state its ID and path and
that it was not opened in Dash.

## 6. Communication and safety checks

Keep the user oriented without narrating every low-level request.  Say which
dashboard was selected, what is about to run, and where the result is visible.
For a potentially expensive compile, long SCM run, or tuning job, confirm that
the requested scope matches the visible settings before starting it.  Do not
claim that a plot or report was displayed until its handoff succeeded and the
matching browser view was confirmed.

When Dash cannot satisfy a request, inspect its declared capabilities once,
state the gap, and then use the least-surprising local alternative only under
the discovery/fallback rule above.  Prefer extending the declarative Dash tab
operation surface when maintaining the application; do not accumulate a
parallel collection of agent-specific shell shortcuts.
