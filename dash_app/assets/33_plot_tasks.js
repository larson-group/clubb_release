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
    function status(task, graph, text) {
        if (graph) set({type: "plots-load-message", index: graph.index}, "children", text);
        else if (task === "case") set("plots-load-status", "children", text);
        else if (task === "catalog") set("plots-catalog-status", "children", text);
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
            return scopes.get(entry.scope) === entry &&
                (entry.task === "catalog" || entry.selection === selectionRevision);
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
                status(entry.task, entry.graph, "Could not load: " + error.message);
                if (entry.task === "case") {
                    status("params", null, "Could not update selection.");
                    for (const old of scopes.values()) {
                        if (old.graph) status(old.task, old.graph, "Could not update selection.");
                    }
                }
                entry.done = true;
            }
        }
    }
    function run(task, args, specs, graph, context) {
        const unchanged = () => specs.length === 1 ? noUpdate() : specs.map(noUpdate);
        if (task === "catalog" && context.triggered_id === "plots-output-menu-expanded" && !args[0]) return unchanged();
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
        if (graph && awaitingCase) {
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
        const inputKey = JSON.stringify(context.inputs_list || args);
        // Parameter controls can emit their unchanged initial values after a
        // case loads. They must not restart a pending or completed figure.
        if (graph && previous && previous.selection === selectionRevision &&
            previous.inputKey === inputKey && (!previous.done || previous.success)) return unchanged();
        const revision = previous ? previous.revision + 1 : 1;
        const outputs = specs.map(out => ({...out, id: graph ? {...out.id, index: graph.index} : out.id}));
        const triggered = (context.triggered || []).map(item => ({prop_id: item.prop_id, value: item.value}));
        const key = JSON.stringify((context.inputs_list || []).filter(input => input.id !== "plots-global-time-point"));
        const entry = {scope, revision, task, graph, outputs, selection: selectionRevision, key, inputKey,
            figure: previous && previous.figure, done: false, success: false};
        entry.request = {task, scope, revision, args, triggered, selection: selectionRevision,
            plot_index: graph ? graph.index : null,
            can_patch: Boolean(graph && entry.figure && previous && previous.done && previous.success && previous.key === key)};
        scopes.set(scope, entry);
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
    function invalidateSelection() {
        selectionRevision++;
        awaitingCase = true;
        status("params", null, "Updating selection…");
        for (const entry of scopes.values()) {
            if (entry.graph) status(entry.task, entry.graph, "Selection changed; updating…");
        }
    }
    window.clubbPlotTasks = {run, invalidateSelection, patchFigure, scopes};
})();
