/* Native Plot work: immediate callbacks, bounded asynchronous requests, latest result wins. */
(function () {
    "use strict";
    const page = (globalThis.crypto && crypto.randomUUID) ? crypto.randomUUID() : Math.random().toString(36).slice(2);
    const scopes = new Map();
    const pending = new Map();
    const polls = [];
    let polling = false;
    let selectionRevision = 0;
    let awaitingCase = false;
    const noUpdate = () => window.dash_clientside.no_update;
    const componentId = id => typeof id === "string" ? id : JSON.stringify(id, Object.keys(id).sort());
    const set = (id, property, value) => window.dash_clientside.set_props(id, {[property]: value});
    const isNoUpdate = value => value && value._dash_no_update !== undefined;
    function benchmarkRequest(caseData, source) {
        // Identify source readiness by case name and configured file paths.
        return {name: caseData?.name ?? null, source, paths: caseData?.benchmarks?.[source] || {}};
    }
    function benchmarkNeeded(caseData, selected, plotState, source) {
        // Prepare selected overlays and sources referenced by custom expressions.
        const available = caseData?.benchmarks?.available_sources || [];
        if (!available.includes(source)) return false;
        if ((selected || []).includes(source)) return true;
        // Custom expressions choose SAM, or the sole other configured source.
        const customSource = available.includes("sam") ? "sam" : available.length === 1 ? available[0] : null;
        return source === customSource && Object.values(plotState || {}).some(
            state => state.plot_type === "custom" && /\bles\s*\./.test(state.expression || ""));
    }
    function sameBenchmarkRequest(actual, expected) {
        // HTTP serialization may reorder fields and file labels without changing the source.
        if (!actual || actual.name !== expected.name || actual.source !== expected.source) return false;
        const paths = actual.paths || {}, expectedPaths = expected.paths || {};
        return Object.keys(paths).length === Object.keys(expectedPaths).length &&
            Object.entries(expectedPaths).every(([label, path]) => paths[label] === path);
    }
    function benchmarkReadiness(caseData, selected, sam, coamps, previous) {
        // Publish matching prepared sources; unchanged readiness cannot restart plots.
        const ready = [sam, coamps].filter(state => state?.state === "ready" &&
            sameBenchmarkRequest(state.request, benchmarkRequest(caseData, state.request?.source)))
            .map(state => state.request.source);
        const value = {name: caseData?.name ?? null, benchmarks: caseData?.benchmarks || null,
            ready, sources: (selected || []).filter(source => ready.includes(source))};
        return JSON.stringify(value) === JSON.stringify(previous) ? noUpdate() : value;
    }
    const taskMessages = new Map();
    function status(task, graph, text, error = false) {
        // Loading and errors share an overlay; neither changes the figure's size.
        if (graph) window.dash_clientside.set_props({type: "plots-load-message", index: graph.index}, {
            children: text,
            className: "plots-load-message" + (text ? (error ? " plots-load-message--error" : " plots-load-message--busy") : ""),
            role: error ? "alert" : "status",
            "aria-live": error ? "assertive" : "polite",
            "aria-label": text,
        });
        else if (task === "case") set("plots-load-status", "children", text);
        else if (task === "catalog") set("plots-catalog-status", "children", text);
        else if (task === "pyplotgen_help") set("plots-pyplotgen-help-status", "children", text);
        else if (task.startsWith("runner_")) {
            taskMessages.set(task, text ? (error ? text : "Checking runners…") : "");
            set("runner-loading", "children", [...taskMessages].filter(([key]) => key.startsWith("runner_")).map(([, value]) => value).filter(Boolean).join(" "));
        }
        else if (task.startsWith("tune_")) {
            taskMessages.set(task, text ? (error ? text : "Preparing plots…") : "");
            set("tune-postprocessing-status", "children", [...taskMessages].filter(([key]) => key.startsWith("tune_")).map(([, value]) => value).filter(Boolean).join(" "));
        }
        else if (task === "params") {
            set("plots-param-status", "children", text);
            set("plots-param-panel", "className", text ? "plots-param-panel--loading" : "");
        }
    }
    async function exchange(payload) {
        const response = await fetch("/plots/tasks", {
            method: "POST", headers: {"Content-Type": "application/json"}, body: JSON.stringify(payload)
        });
        const data = await response.json();
        if (!response.ok || data.error && !data.state) throw new Error(data.error || "Could not load Plot data");
        return data;
    }
    function poll(request) {
        return new Promise((resolve, reject) => {
            polls.push({request, resolve, reject});
            schedulePoll();
        });
    }
    function schedulePoll() {
        if (polling || !polls.length) return;
        polling = true;
        setTimeout(async () => {
            const batch = polls.splice(0, 128);
            try {
                const response = await exchange({operation: "poll", requests: batch.map(item => item.request)});
                if (!Array.isArray(response.results) || response.results.length !== batch.length) {
                    throw new Error("Unexpected Plot response; restart the dashboard runtime.");
                }
                batch.forEach((item, index) => item.resolve(response.results[index]));
            } catch (error) {
                batch.forEach(item => item.reject(error));
            } finally {
                polling = false;
                schedulePoll();
            }
        }, 120);
    }
    function patchFigure(base, patch) {
        // Native playback patches currently assign complete trace arrays/ranges.
        const result = structuredClone(base);
        for (const op of patch.operations || []) {
            if (op.operation !== "Assign") throw new Error("Unsupported plot patch: " + op.operation);
            let target = result;
            const path = op.location;
            for (const key of path.slice(0, -1)) {
                if (target[key] === undefined) target[key] = {};
                target = target[key];
            }
            target[path[path.length - 1]] = op.params.value;
        }
        return result;
    }
    async function launch(entry) {
        const current = () => {
            if (entry.graph && typeof document !== "undefined" && !document.getElementById(componentId(entry.graph))) {
                if (scopes.get(entry.scope) === entry) scopes.delete(entry.scope);
                return false;
            }
            return !entry.cancelled && scopes.get(entry.scope) === entry &&
                (["catalog", "pyplotgen_help", "runner_info", "runner_catalog", "tune_controls", "tune_landscape", "tune_diagnostics"].includes(entry.task) || entry.selection === selectionRevision);
        };
        try {
            const accepted = await exchange({operation: "submit", ...entry.request});
            if (accepted.state === "superseded") return;
            while (current()) {
                const result = await poll({scope: entry.scope, revision: entry.revision});
                if (!current()) return;
                if (result.state === "missing") {
                    await exchange({operation: "submit", ...entry.request});
                    continue;
                }
                if (result.state === "superseded") return;
                if (result.state === "error") throw new Error(result.error || "Loading failed");
                if (result.state !== "ready") continue;
                const values = entry.outputs.length === 1 ? [result.result] : result.result;
                const updates = [];
                values.forEach((value, index) => {
                    if (isNoUpdate(value)) return;
                    const out = entry.outputs[index];
                    if (out.property === "figure") {
                        if (value && value.__dash_patch_update) value = patchFigure(entry.figure, value);
                        entry.figure = value;
                    }
                    updates.push([out, value]);
                });
                if (!current()) return;
                if (entry.task === "case") {
                    const selected = updates.find(([out]) => out.id === "plots-case-data");
                    if (selected && !selected[1]?.preserve_plot_view) {
                        // Cards stay mounted; reset the previous case's zoom
                        // state before the new figure callbacks read it.
                        for (const old of scopes.values()) {
                            if (old.graph) set(old.graph, "relayoutData", null);
                        }
                    }
                    awaitingCase = false;
                }
                updates.forEach(([out, value]) => set(out.id, out.property, value));
                status(entry.task, entry.graph, "");
                entry.done = true;
                entry.success = true;
                return;
            }
        } catch (error) {
            if (current()) {
                if (entry.benchmark) set("plots-benchmark-" + entry.benchmark.source, "data",
                    {request: entry.benchmark, state: "error", error: error.message});
                status(entry.task, entry.graph, "Could not load: " + error.message, true);
                if (entry.task === "case") {
                    status("params", null, "Could not update selection.");
                    for (const old of scopes.values()) {
                        if (old.graph) status(old.task, old.graph, "Could not update selection.", true);
                    }
                }
                entry.done = true;
            }
        } finally {
            entry.done = true;
            const latest = entry.tunePendingUpdate;
            if (latest && current()) run(latest.task, latest.args, latest.specs, latest.graph, latest.context);
        }
    }
    function tuneDataKey(task, args) {
        const status = args[task === "tune_diagnostics" ? args.length - 2 : 0] || {};
        if (task === "tune_diagnostics") return JSON.stringify([status.job_dir, args[0] || [], args[1] || {}]);
        return JSON.stringify([status.job_dir,
            status.sample_history_signature ?? status.samples_evaluated ?? 0,
            status.loss_weights, status.aggregation_mode, status.case_window_counts, status.selected_fields]);
    }
    function run(task, args, specs, graph, context) {
        const unchanged = () => specs.length === 1 ? noUpdate() : specs.map(noUpdate);
        if (task.startsWith("runner_") && !args[0]) return unchanged();
        if (task.startsWith("tune_")) {
            const triggered = context.triggered_id;
            const inputValues = (context.inputs_list || []).filter(item => !Array.isArray(item));
            const enabled = inputValues.find(item => item.id === "tune-diagnostic-plots")?.value || [];
            const tab = inputValues.find(item => item.id === "dashboard-tabs")?.value;
            const families = task === "tune_diagnostics" ? ["taylor", "parameters"] : ["landscape", "correlation", "sensitivity", "interaction"];
            if (tab !== "tune" || !enabled.some(key => families.includes(key))) {
                const old = scopes.get(page + ":" + task);
                if (old) { clearTimeout(pending.get(old.scope)); old.revision++; old.cancelled = true; old.done = true; }
                status(task, null, "");
                return unchanged();
            }
            const previous = scopes.get(page + ":" + task);
            const statusIndex = task === "tune_diagnostics" ? args.length - 2 : 0;
            const jobDir = args[statusIndex]?.job_dir || null;
            const dataTriggers = ["tune-status", "tune-best-results", "tune-best-results-by-case"];
            const dataRefresh = dataTriggers.includes(triggered) && (context.triggered || []).every(item =>
                dataTriggers.includes(item.prop_id?.slice(0, item.prop_id.lastIndexOf("."))));
            if (dataRefresh && previous && previous.jobDir === jobDir) {
                const key = tuneDataKey(task, args);
                if (!previous.done) {
                    if (key !== previous.tuneDataKey) previous.tunePendingUpdate = {
                        task, args, specs, graph,
                        context: {triggered_id: context.triggered_id,
                            triggered: (context.triggered || []).map(item => ({...item})),
                            inputs_list: context.inputs_list}
                    };
                    return unchanged();
                }
                if (previous.success && previous.tuneDataKey === key) return unchanged();
            }
        }
        const benchmark = task.startsWith("benchmark_") ? benchmarkRequest(args[0], task.slice(10)) : null;
        if (task === "case") {
            // A directory change can coincide with zero-click button mounts.
            // Preserve the real selection event when Dash batches those inputs.
            const directory = (context.triggered || []).find(item => item.prop_id === "plots-output-dirs.data");
            if (directory) context = {...context, triggered_id: "plots-output-dirs",
                triggered: [directory, ...context.triggered.filter(item => item !== directory)]};
            // Mounting/rebuilding case buttons emits zero-click events.
            if (context.triggered_id && typeof context.triggered_id === "object" &&
                !(Number((context.triggered || [])[0]?.value) > 0)) return unchanged();
            if (context.triggered_id === "dashboard-request" &&
                (!args[2] || args[2].tab !== "plots" || !["set_view", "add_budget"].includes(args[2].operation))) return unchanged();
            invalidateSelection();
        }
        if ((graph || benchmark) && awaitingCase) {
            status(task, graph, "Updating selection…");
            return unchanged();
        }
        // Release completed figures for cards that have been removed.
        if (typeof document !== "undefined") {
            for (const [key, old] of scopes) {
                if (old.graph && !document.getElementById(componentId(old.graph))) scopes.delete(key);
            }
        }
        const scope = page + ":" + task + (graph ? ":" + graph.index : "");
        const previous = scopes.get(scope);
        const needed = benchmark && benchmarkNeeded(args[0], args[1], args[2], benchmark.source);
        const inputKey = benchmark ? JSON.stringify([benchmark, needed]) : JSON.stringify(context.inputs_list || args);
        // Parameter controls can emit their unchanged initial values after a
        // case loads. They must not restart a pending or completed figure.
        if ((graph || benchmark) && previous && previous.selection === selectionRevision &&
            previous.inputKey === inputKey && (!previous.done || previous.success)) return unchanged();
        const revision = previous ? previous.revision + 1 : 1;
        const outputs = specs.map(out => ({...out, id: graph ? {...out.id, index: graph.index} : out.id}));
        const triggered = (context.triggered || []).map(item => ({prop_id: item.prop_id, value: item.value}));
        const key = JSON.stringify((context.inputs_list || []).filter(input => input.id !== "plots-global-time-point"));
        const entry = {scope, revision, task, graph, outputs, benchmark,
            tuneDataKey: task.startsWith("tune_") ? tuneDataKey(task, args) : null,
            jobDir: task.startsWith("tune_") ? (args[task === "tune_diagnostics" ? args.length - 2 : 0]?.job_dir || null) : null,
            selection: selectionRevision, key, inputKey,
            figure: previous && previous.figure, done: false, success: false};
        entry.request = {task, scope, revision, args, triggered, selection: selectionRevision,
            plot_index: graph ? graph.index : null,
            can_patch: Boolean(graph && entry.figure && previous && previous.done && previous.success && previous.key === key)};
        scopes.set(scope, entry);
        if (benchmark) {
            clearTimeout(pending.get(scope));
            pending.delete(scope);
            if (!needed) {
                entry.done = entry.success = true;
                return {request: benchmark, state: "disabled"};
            }
            set("plots-benchmark-" + benchmark.source, "data", {request: benchmark, state: "loading"});
        }
        status(task, graph, "Updating…");
        if (task === "case" && typeof context.triggered_id === "object" && context.triggered_id?.name) {
            const name = context.triggered_id.name;
            status(task, null, "Loading " + name + "…");
            if (typeof document !== "undefined") {
                document.getElementById("plots-case-button-container")?.querySelectorAll("button").forEach(button => {
                    button.style.borderColor = JSON.parse(button.id).name === name ? "#f59e0b" : "transparent";
                });
            }
        }
        // Coalesce rapid clicks/slider changes before reaching the broker.
        clearTimeout(pending.get(scope));
        pending.set(scope, setTimeout(() => {pending.delete(scope); launch(entry);}, graph ? 45 : 0));
        return specs.length === 1 ? noUpdate() : specs.map(noUpdate);
    }
    function cancel(task) {
        const entry = scopes.get(page + ":" + task);
        if (!entry) return;
        clearTimeout(pending.get(entry.scope));
        pending.delete(entry.scope);
        entry.revision++;
        entry.cancelled = entry.done = true;
        status(task, null, "");
    }
    function invalidateSelection() {
        selectionRevision++;
        awaitingCase = true;
        status("params", null, "Updating selection…");
        for (const entry of scopes.values()) {
            if (entry.graph) status(entry.task, entry.graph, "Selection changed; updating…");
        }
    }
    window.clubbPlotTasks = {run, cancel, invalidateSelection, patchFigure, scopes, benchmarkReadiness};
})();
