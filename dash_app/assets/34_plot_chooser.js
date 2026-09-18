/* Selection and menu feedback stays in the browser while workers load data. */
(function () {
    "use strict";
    const no = () => window.dash_clientside.no_update;
    const el = (type, props, children) => ({namespace: "dash_html_components", type, props: {...props, children}});
    const ctx = () => window.dash_clientside.callback_context;
    const clicked = () => ctx().triggered && Number(ctx().triggered[0]?.value) > 0;
    function normalize(path, paths) {
        path = String(path || "").trim();
        if (path === "~" || path.startsWith("~/")) path = paths.home + path.slice(1);
        if (!path.startsWith("/")) path = paths.repo + "/" + path;
        const parts = [];
        path.split("/").forEach(part => {
            if (part === "..") parts.pop();
            else if (part && part !== ".") parts.push(part);
        });
        return "/" + parts.join("/");
    }
    function label(record, path) {
        const value = String(record.label || record.relative_path || "").replace(/^\/+|\/+$/g, "").replace(/^output\//, "");
        return value || path.split("/").filter(Boolean).pop() || path;
    }
    function summary(record) {
        return record.error || "Cases: " + ((record.case_names || []).join(", ") || "none currently available");
    }
    const buttonText = (record, path) => label(record, path) + " · " + (record.case_count || (record.case_names || []).length);
    function deleteButton(path, confirmation, paths) {
        const allowed = path.startsWith(paths.output + "/");
        const confirming = allowed && confirmation && confirmation.path === path;
        return el("Button", {
            id: {type: "plots-delete-output-dir", path}, n_clicks: 0, n_clicks_timestamp: -1,
            disabled: !allowed,
            title: allowed ? (confirming ? "Click again within 3 seconds to permanently delete this output" : "Permanently delete this output directory") : "Only subdirectories inside output/ can be deleted here",
            className: "plots-output-dir-delete" + (confirming ? " plots-output-dir-delete--confirm" : "")
        }, confirming ? "Confirm" : "×");
    }
    function render(records, selected, expanded, confirmation, paths) {
        selected = selected || []; records = records || [];
        const byPath = new Map(records.map(record => [record.path, record]));
        let available = records.filter(record => record.available !== false && !selected.includes(record.path));
        if (!expanded) available = available.slice(0, 3);
        available = available.map(record => {
            const path = record.path;
            const button = el("Button", {id: {type: "plots-add-output-dir", path}, n_clicks: 0, n_clicks_timestamp: -1,
                title: summary(record), className: "plots-output-available-button"}, buttonText(record, path));
            return expanded ? el("Div", {className: "plots-output-menu-row"}, [button, deleteButton(path, confirmation, paths)]) : button;
        });
        if (!available.length) available = [el("Div", {className: "plots-output-menu-empty"}, "All outputs selected")];
        const active = selected.map(path => {
            const record = byPath.get(path);
            const unavailable = record && record.available === false;
            const children = [el("Span", {className: "plots-output-active-label"}, record ? buttonText(record, path) : path)];
            if (!record || unavailable) children.push(el("Span", {className: "plots-output-active-status"}, record ? "Unavailable" : "Loading…"));
            if (expanded) children.push(deleteButton(path, confirmation, paths));
            children.push(el("Button", {id: {type: "plots-remove-output-dir", path}, n_clicks: 0, n_clicks_timestamp: -1,
                title: "Remove from the comparison", className: "plots-output-dir-remove"}, "−"));
            return el("Div", {title: record ? summary(record) : "Loading directory…",
                className: "plots-output-active-item" + (unavailable ? " plots-output-active-item--unavailable" : "")}, children);
        });
        if (!active.length) active.push(el("Small", {className: "plots-mutable-output-warning"}, "No output folders selected."));
        return [available, active, "plots-output-menu" + (expanded ? " plots-output-menu--expanded" : "")];
    }
    function choose(_add, _extra, _remove, current, extra, catalog, paths) {
        if (!clicked()) return [no(), no(), no(), no()];
        let selected = (current || []).slice();
        const trigger = ctx().triggered_id;
        if (trigger === "plots-add-extra-dir") {
            if (!String(extra || "").trim()) return [no(), no(), "Enter a directory path before adding it.", no()];
            const path = normalize(extra, paths);
            if (selected.includes(path)) return [no(), no(), "That folder is already selected.", no()];
            selected.push(path);
        } else if (trigger && trigger.type === "plots-add-output-dir") {
            if (selected.includes(trigger.path) || !(catalog || []).some(record => record.path === trigger.path)) return [no(), no(), no(), no()];
            selected.push(trigger.path);
        } else if (trigger && trigger.type === "plots-remove-output-dir") {
            selected = selected.filter(path => path !== trigger.path);
        } else return [no(), no(), no(), no()];
        window.clubbPlotTasks.invalidateSelection();
        return [selected, trigger === "plots-add-extra-dir" ? "" : no(), "", no()];
    }
    window.dash_clientside = Object.assign({}, window.dash_clientside, {plotsChooser: {
        toggle: function (_timestamp, expanded) {return clicked() ? !expanded : no();},
        showExtra: function (clicks, current) {return clicks ? {...current, display: "flex", marginTop: "8px", alignItems: "center", gap: "8px"} : no();},
        render, choose
    }});
})();
