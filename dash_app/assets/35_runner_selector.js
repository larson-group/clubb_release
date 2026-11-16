(function () {
  "use strict";
  const ids = ["run-selected-build-badge", "tune-selected-build-badge", "profile-selected-build-badge"];
  let currentAnchor = null;
  let view = {views: null, anchor: null, implementation: "fortran", profile: "cpu", gpu: "", prealloc: false};
  let pending = {};
  let published = {implementation: "fortran", profile: "cpu", gpu: "", prealloc: false};
  const selectionStores = {gpu: "compile-run-jax-gpu", profile: "compile-run-jax-profile",
    implementation: "compile-run-implementation", prealloc: "compile-run-jax-xla-prealloc"};
  let commitFrame = null;
  function commitSelection() {
    Object.entries(selectionStores).forEach(([key, id]) => {
      if (Object.hasOwn(pending, key) && pending[key] !== published[key]) {
        published[key] = pending[key];
        dash_clientside.set_props(id, {data: pending[key]});
      }
    });
  }
  function choose(changes) {
    Object.assign(pending, changes);
    Object.assign(view, changes);
    renderView();
    afterPaint();
  }
  function afterPaint() {
    if (commitFrame !== null) return;
    commitFrame = requestAnimationFrame(() => {
      commitFrame = requestAnimationFrame(() => {commitFrame = null; commitSelection();});
    });
  }
  function show(anchor) {
    currentAnchor = anchor;
    view.anchor = anchor;
    const api = window.dash_clientside;
    api.set_props("compile-build-selector-popover", {className: "compile-build-selector-popover" + (anchor ? "" : " compile-build-selector-popover-closed")});
    if (!anchor) return;
    const width = Math.min(520, window.innerWidth - 16);
    const below = window.innerHeight - anchor.bottom - 8;
    const above = anchor.top - 8;
    const style = {width: width + "px", left: Math.max(8, Math.min(anchor.left, window.innerWidth - width - 8)) + "px"};
    if (below < 220 && above > below) {
      style.bottom = window.innerHeight - anchor.top + 6 + "px";
      style.maxHeight = Math.max(120, above - 6) + "px";
    } else {
      style.top = anchor.bottom + 6 + "px";
      style.maxHeight = Math.max(120, below - 6) + "px";
    }
    api.set_props("runner-surface", {style: style});
  }
  function nextRequest(previous, changes) {
    return Object.assign({}, previous, changes, {revision: ((previous || {}).revision || 0) + 1});
  }
  function markChoice(type, selected) {
    document.querySelectorAll('[id*="' + type + '"]').forEach(button => {
      let id;
      try { id = JSON.parse(button.id); } catch (_) { return; }
      if (id.type !== type) return;
      const active = selected(id.index);
      button.classList.toggle("compile-run-implementation-choice-selected", active);
      button.setAttribute("aria-pressed", String(active));
    });
  }
  function renderView() {
    const {views, anchor, implementation, profile, gpu, prealloc} = view;
    const jax = implementation === "jax";
    const fortranPanel = document.getElementById("runner-fortran-options");
    const jaxPanel = document.getElementById("runner-jax-options");
    if (fortranPanel) fortranPanel.style.display = jax ? "none" : "";
    if (jaxPanel) jaxPanel.style.display = jax ? "" : "none";
    markChoice("compile-run-implementation-choice", value => value === implementation);
    const selected = (((views || {}).jax || {})[profile === "cpu" ? "cpu" : (gpu || "gpu")]);
    markChoice("compile-run-jax-profile-choice", value => profile === "cpu" ? value === "cpu" : value === (gpu || "gpu"));
    const checkbox = document.querySelector('[id*="compile-jax-prealloc-choice"]');
    if (checkbox) {
      const input = checkbox.querySelector("input");
      const disabled = !selected || !selected.can_preallocate || profile !== "gpu";
      if (input) { input.disabled = disabled; input.checked = Boolean(prealloc); }
    }
    const note = document.getElementById("runner-context-note");
    if (note) note.textContent = anchor && anchor.trigger_id === "profile-selected-build-badge" && jax ?
        "JAX does not yet write the timer files required by Profile." :
        anchor && anchor.trigger_id === "tune-selected-build-badge" ? "Tune uses the selected Fortran or JAX loss backend." : "";
    window.clubbRunnerSelector.paint(jax ? selected : (views || {}).fortran, anchor);
  }
  window.clubbRunnerSelector = {
    chooseImplementation: function (implementation) {
      if (["fortran", "jax"].includes(implementation)) choose({implementation});
    },
    chooseCompute: function (device) {
      if (device === "cpu") choose({profile: "cpu"});
      else if (device === "gpu") choose({profile: "gpu", gpu: ""});
      else if (((view.views || {}).jax || {})[device]) choose({profile: "gpu", gpu: device});
    },
    chooseAllocation: function (enabled) {
      const selected = ((view.views || {}).jax || {})[view.gpu || "gpu"];
      if (view.profile === "gpu" && selected && selected.can_preallocate) choose({prealloc: Boolean(enabled)});
    },
    synchronize: function (views, anchor, implementation, profile, gpu, prealloc) {
      const incoming = {implementation, profile, gpu, prealloc};
      Object.keys(incoming).forEach(key => {
        if (pending[key] === incoming[key]) delete pending[key];
        if (!Object.hasOwn(pending, key)) published[key] = incoming[key];
      });
      view = Object.assign({views, anchor}, incoming, pending);
      renderView();
      return dash_clientside.no_update;
    },
    runtimeRequest: function (refresh, environment, previous) {
      refresh = Number(refresh) || 0;
      const revision = Number((environment || {}).revision) || 0;
      if (!previous || refresh !== previous.refresh) {
        return nextRequest(previous, {kind: "inventory", refresh: refresh, environment: revision});
      }
      if (revision !== previous.environment) {
        return nextRequest(previous, {kind: "environment", environment: revision});
      }
      return window.dash_clientside.no_update;
    },
    catalogRequest: function (refresh, jobs, mounted, previous) {
      refresh = Number(refresh) || 0;
      const record = (jobs || {}).compile || {};
      const job = record.job || {};
      const completion = JSON.stringify([record.state || "", job.pid, job.start_time]);
      if (!previous || refresh !== previous.refresh) {
        return nextRequest(previous, {refresh: refresh, completion: completion});
      }
      if (!mounted && ["finished", "failed", "error", "cancelled", "stopped"].includes(record.state) &&
          completion !== previous.completion) {
        return nextRequest(previous, {completion: completion});
      }
      return window.dash_clientside.no_update;
    },
    compileState: function (jobs, previous) {
      const record = (jobs || {}).compile || {};
      const next = {compile: {state: record.state || "", job: record.job || null}};
      return JSON.stringify(next) === JSON.stringify(previous) ? window.dash_clientside.no_update : next;
    },
    position: function () {
      const ctx = window.dash_clientside.callback_context || {};
      const event = (ctx.triggered || []).find(item => Number(item.value) > 0);
      if (!event || !event.prop_id) return window.dash_clientside.no_update;
      const raw = event.prop_id.slice(0, event.prop_id.lastIndexOf("."));
      const id = raw.startsWith("{") ? JSON.parse(raw) : raw;
      if (id === "compile-build-selector-dismiss" || typeof id === "object") { show(null); return null; }
      const button = document.getElementById(id);
      if (!button) return window.dash_clientside.no_update;
      if (currentAnchor && currentAnchor.trigger_id === id) { show(null); return null; }
      const rect = button.getBoundingClientRect();
      const anchor = {trigger_id: id, left: rect.left, top: rect.top, bottom: rect.bottom,
        viewport_width: window.innerWidth, viewport_height: window.innerHeight};
      show(anchor);
      return anchor;
    },
    paint: function (state, anchor) {
      ids.forEach(function (id) {
        const button = document.getElementById(id);
        if (!button) return;
        button.setAttribute("aria-expanded", String(Boolean(anchor && anchor.trigger_id === id)));
        if (!state) return;
        button.className = state.className + (anchor && anchor.trigger_id === id ? " selected-build-badge-open" : "");
        button.title = state.title;
        const label = button.querySelector(".selected-build-badge-label");
        const name = button.querySelector(".selected-build-badge-name");
        if (label && label.textContent !== state.implementation.toUpperCase()) label.textContent = state.implementation.toUpperCase();
        if (name && name.textContent !== state.name) name.textContent = state.name;
      });
      return dash_clientside.no_update;
    }
  };
  document.addEventListener("click", function (event) {
    const button = event.target.closest && event.target.closest("button");
    let id;
    try { id = button && JSON.parse(button.id); } catch (_) {}
    if (id && ["compile-run-implementation-choice", "compile-run-jax-profile-choice"].includes(id.type)) {
      if (button.disabled) return;
      event.preventDefault();
      event.stopImmediatePropagation();
      if (id.type === "compile-run-implementation-choice") window.clubbRunnerSelector.chooseImplementation(id.index);
      else window.clubbRunnerSelector.chooseCompute(id.index);
    } else {
      // Publish before another control, including Launch, can read the stores.
      commitSelection();
    }
  }, true);
  document.addEventListener("change", function (event) {
    const control = event.target.closest && event.target.closest('[id*="compile-jax-prealloc-choice"]');
    if (!control || event.target.disabled) return;
    event.stopImmediatePropagation();
    window.clubbRunnerSelector.chooseAllocation(event.target.checked);
  }, true);
}());
