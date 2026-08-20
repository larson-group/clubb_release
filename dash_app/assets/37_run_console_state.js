(function () {
  "use strict";

  const CONTAINER_ID = "run-console-container";
  const TELEMETRY_URL = "/_clubb-run-telemetry";
  const ACTIVE_POLL_MS = 500;
  const IDLE_POLL_MS = 5000;
  const HIDDEN_CHECK_MS = 1000;
  const MAX_LINES = 5000;
  const NEAR_BOTTOM_PX = 32;
  const ACTIVE_STATES = new Set(["queued", "running", "stopping"]);
  const STATE_LABELS = {
    queued: "queued",
    running: "running",
    stopping: "stopping",
    finished: "finished",
    error: "failed",
    cancelled: "cancelled",
  };

  const runs = new Map();
  const panels = new Map();
  const buffers = new Map();
  const progressByRun = new Map();
  const openRuns = new Set();
  const transientErrors = new Map();
  let selectedCases = new Set();
  let knownRevision = null;
  let requestActive = false;
  let pollTimer = null;
  let submitting = false;
  let targetOutputDirectory = null;
  let renderSequence = 0;

  function element(tag, className, text) {
    const node = document.createElement(tag);
    if (className) node.className = className;
    if (text !== undefined) node.textContent = text;
    return node;
  }

  function runTabVisible() {
    const layout = document.getElementById("run-tab-layout");
    return Boolean(layout && layout.offsetParent !== null);
  }

  function formatDuration(seconds) {
    const total = Math.max(0, Math.round(Number(seconds) || 0));
    const hours = Math.floor(total / 3600);
    const minutes = Math.floor((total % 3600) / 60);
    const remainder = total % 60;
    if (hours) return `${hours}h ${minutes}m ${remainder}s`;
    if (minutes) return `${minutes}m ${remainder}s`;
    return `${remainder}s`;
  }

  function setText(node, value) {
    if (node && node.textContent !== value) node.textContent = value;
  }

  function newBuffer(node) {
    return {
      node: node,
      cursor: null,
      decoder: new TextDecoder("utf-8"),
      terminalEof: false,
      newlineCount: 0,
    };
  }

  function createPanel(run) {
    const panel = element("details", "run-console-panel");
    panel.dataset.runId = run.run_id;
    const summary = element("summary", "run-console-summary");
    const left = element("div", "run-summary-left");
    left.appendChild(element("span", "run-console-case-name", run.case));
    left.appendChild(element("span", "run-lifecycle-status", "queued"));
    left.appendChild(element("span", "run-muted-text run-console-origin", ""));
    left.appendChild(element("span", "run-muted-text run-console-runtime", ""));
    const progressBlock = element("div", "run-progress-block");
    progressBlock.appendChild(element("span", "run-muted-text run-console-eta", "eta --"));
    const track = element("div", "run-progress-track");
    track.appendChild(element("div", "run-progress-fill"));
    progressBlock.appendChild(track);
    summary.appendChild(left);
    summary.appendChild(progressBlock);
    panel.appendChild(summary);

    const commandRow = element("div", "run-command-row");
    commandRow.hidden = true;
    commandRow.appendChild(element("div", "run-command-text", ""));
    const copy = element("button", "run-copy-command-button", "Copy run command");
    copy.type = "button";
    commandRow.appendChild(copy);
    panel.appendChild(commandRow);
    const error = element("div", "run-console-error", "");
    error.hidden = true;
    panel.appendChild(error);
    const consoleNode = element("pre", "run-console", "");
    panel.appendChild(consoleNode);

    panel.addEventListener("toggle", function () {
      if (panel.open) {
        openRuns.add(run.run_id);
        if (!buffers.has(run.run_id)) buffers.set(run.run_id, newBuffer(consoleNode));
        schedulePoll(0, true);
      } else {
        openRuns.delete(run.run_id);
        buffers.delete(run.run_id);
        consoleNode.replaceChildren();
      }
    });
    copy.addEventListener("click", function () {
      const command = panel.querySelector(".run-command-text");
      if (command && command.textContent) navigator.clipboard.writeText(command.textContent);
    });
    return panel;
  }

  function updateStatus(node, state) {
    if (!node) return;
    setText(node, STATE_LABELS[state] || state || "unknown");
    node.style.backgroundImage = state === "running"
      ? "linear-gradient(90deg, #dc2626, #16a34a)"
      : "";
    node.style.webkitBackgroundClip = state === "running" ? "text" : "";
    node.style.webkitTextFillColor = state === "running" ? "transparent" : "";
    node.style.color = state === "running" ? "" : ({
      queued: "#d97706",
      stopping: "#ea580c",
      finished: "#16a34a",
      error: "#dc2626",
      cancelled: "#dc2626",
    }[state] || "#2563eb");
  }

  function updatePanel(panel, run) {
    panel.dataset.runState = run.state || "";
    panel.dataset.runStartTime = String(run.started_at || "");
    panel.dataset.runFinishedAt = String(run.finished_at || "");
    panel.dataset.logAvailable = run.log_available ? "true" : "false";
    setText(panel.querySelector(".run-console-case-name"), run.case || "unknown case");
    updateStatus(panel.querySelector(".run-lifecycle-status"), run.state);
    setText(panel.querySelector(".run-console-origin"), run.origin ? `via ${run.origin}` : "");
    const command = String(run.command || "");
    const commandRow = panel.querySelector(".run-command-row");
    if (commandRow) commandRow.hidden = !command;
    setText(panel.querySelector(".run-command-text"), command);
    const message = String(run.message || "");
    const error = panel.querySelector(".run-console-error");
    if (error) error.hidden = !message;
    setText(error, message);
  }

  function renderRuns() {
    const container = document.getElementById(CONTAINER_ID);
    if (!container) return;
    const displayed = new Map(runs);
    transientErrors.forEach(function (value, key) { displayed.set(key, value); });

    Array.from(panels.keys()).forEach(function (runId) {
      if (displayed.has(runId)) return;
      const panel = panels.get(runId);
      if (panel) panel.remove();
      panels.delete(runId);
      buffers.delete(runId);
      progressByRun.delete(runId);
      openRuns.delete(runId);
    });

    const orderedRuns = Array.from(displayed.entries()).sort(function (left, right) {
      function priority(run) {
        if (run.state === "running" || run.state === "stopping") return 0;
        return run.state === "queued" ? 1 : 2;
      }
      return priority(left[1]) - priority(right[1])
        || Number(right[1].updated_at || 0) - Number(left[1].updated_at || 0);
    });
    orderedRuns.forEach(function (entry, index) {
      const runId = entry[0];
      const run = entry[1];
      let panel = panels.get(runId);
      if (!panel) {
        panel = createPanel(run);
        panels.set(runId, panel);
        container.appendChild(panel);
      }
      updatePanel(panel, run);
      panel.style.order = String(index);
    });
    const empty = document.getElementById("run-console-empty");
    if (empty) empty.style.display = displayed.size ? "none" : "block";
    updateCaseButtons();
    updateDisplays();
    updateControls();
  }

  function matchingCaseRun(caseName) {
    if (!targetOutputDirectory) return null;
    let best = null;
    runs.forEach(function (run) {
      if (run.case !== caseName || run.output_directory !== targetOutputDirectory) return;
      const active = ACTIVE_STATES.has(run.state);
      const bestActive = best && ACTIVE_STATES.has(best.state);
      if (!best || (active && !bestActive)
          || (active === bestActive && Number(run.updated_at || 0) > Number(best.updated_at || 0))) {
        best = run;
      }
    });
    return best;
  }

  function updateCaseButtons() {
    document.querySelectorAll(".run-case-button[data-case-name]").forEach(function (button) {
      const caseName = button.dataset.caseName || "";
      const selected = selectedCases.has(caseName);
      const run = matchingCaseRun(caseName);
      const state = run && run.state;
      let color = "#2563eb";
      if (state === "running") color = "linear-gradient(135deg, #dc2626 0%, #16a34a 100%)";
      else if (state === "queued") color = "#d97706";
      else if (state === "stopping") color = "#ea580c";
      else if (state === "finished") color = "#16a34a";
      else if (state === "error" || state === "cancelled") color = "#dc2626";
      button.style.background = color;
      button.style.backgroundColor = color.includes("gradient") ? "#dc2626" : color;
      button.style.border = `${selected ? 3 : 1}px solid ${selected ? "#f59e0b" : "transparent"}`;
      button.dataset.runState = state || "";
      button.title = run
        ? `${caseName}: ${STATE_LABELS[state] || state} in ${targetOutputDirectory}`
        : caseName;
    });
  }

  function updateControls() {
    const active = Array.from(runs.values()).some(function (run) {
      return ACTIVE_STATES.has(run.state);
    });
    const eligible = Array.from(selectedCases).some(function (caseName) {
      const run = matchingCaseRun(caseName);
      return !run || !ACTIVE_STATES.has(run.state);
    });
    const runButton = document.getElementById("run-button");
    const cancelButton = document.getElementById("run-cancel");
    if (runButton) {
      runButton.disabled = submitting || selectedCases.size === 0
        || !targetOutputDirectory || !eligible;
      runButton.textContent = submitting
        ? "Submitting…"
        : selectedCases.size > 0 && !targetOutputDirectory
        ? "Invalid output"
        : selectedCases.size > 0 && !eligible
        ? "No eligible cases"
        : "Run selected";
    }
    if (cancelButton) cancelButton.disabled = !active;
  }

  function observeProgress(runId, iteration, total) {
    if (!(total > 0)) return;
    const fraction = Math.max(0, Math.min(1, iteration / total));
    const now = performance.now() / 1000;
    let progress = progressByRun.get(runId);
    if (!progress || fraction < progress.fraction) {
      progress = {
        fraction: fraction,
        observedAt: now,
        rate: null,
        increasingSamples: 0,
        reset: Boolean(progress),
      };
    } else if (fraction > progress.fraction && now > progress.observedAt) {
      const instantaneous = (fraction - progress.fraction) / (now - progress.observedAt);
      progress.rate = progress.rate === null
        ? instantaneous
        : 0.3 * instantaneous + 0.7 * progress.rate;
      progress.increasingSamples += 1;
      progress.fraction = fraction;
      progress.observedAt = now;
    }
    progressByRun.set(runId, progress);
  }

  function updateDisplays() {
    const now = Date.now() / 1000;
    panels.forEach(function (panel, runId) {
      const run = runs.get(runId) || transientErrors.get(runId) || {};
      const active = ACTIVE_STATES.has(run.state);
      const start = Number(run.started_at);
      const finish = Number(run.finished_at);
      const progress = progressByRun.get(runId);
      setText(panel.querySelector(".run-console-runtime"), start > 0
        ? `${active ? "elapsed" : "runtime"}: ${formatDuration((active ? now : (finish > 0 ? finish : start)) - start)}`
        : "");
      const fill = panel.querySelector(".run-progress-fill");
      if (fill) {
        const fraction = run.state === "finished" ? 1 : (progress && progress.fraction) || 0;
        fill.style.width = `${Math.round(100 * fraction)}%`;
      }
      const eta = panel.querySelector(".run-console-eta");
      if (eta) {
        const elapsed = start > 0 ? now - start : 0;
        const canEstimate = active && progress && !progress.reset && progress.rate > 0
          && progress.increasingSamples >= 2 && elapsed >= 2;
        setText(eta, canEstimate
          ? `eta ${formatDuration((1 - progress.fraction) / progress.rate)}`
          : "eta --");
        eta.style.display = active ? "" : "none";
      }
    });
  }

  function decodeBase64(value) {
    const binary = atob(value || "");
    const bytes = new Uint8Array(binary.length);
    for (let index = 0; index < binary.length; index += 1) bytes[index] = binary.charCodeAt(index);
    return bytes;
  }

  function appendLog(runId, payload) {
    const buffer = buffers.get(runId);
    const panel = panels.get(runId);
    if (!buffer || !panel || !panel.open) return;
    const bytes = decodeBase64(payload.data);
    const chunk = buffer.decoder.decode(bytes, { stream: !payload.eof });
    const node = buffer.node;
    if (chunk) {
      const nearBottom = node.scrollHeight - node.scrollTop - node.clientHeight <= NEAR_BOTTOM_PX;
      let textNode = node.lastChild;
      if (!textNode || textNode.nodeType !== Node.TEXT_NODE) {
        textNode = document.createTextNode("");
        node.appendChild(textNode);
      }
      textNode.appendData(chunk);
      buffer.newlineCount += (chunk.match(/\n/g) || []).length;
      let excess = buffer.newlineCount - MAX_LINES;
      if (excess > 0) {
        let cut = 0;
        while (excess > 0) {
          cut = textNode.data.indexOf("\n", cut) + 1;
          excess -= 1;
        }
        textNode.deleteData(0, cut);
        buffer.newlineCount = MAX_LINES;
      }
      if (nearBottom) node.scrollTop = node.scrollHeight;
    }
    buffer.cursor = Number(payload.next_cursor) || 0;
    buffer.terminalEof = Boolean(payload.eof) && !ACTIVE_STATES.has(payload.state);
  }

  function logCursors() {
    const result = {};
    openRuns.forEach(function (runId) {
      const run = runs.get(runId);
      const buffer = buffers.get(runId);
      if (run && run.log_available && buffer && !buffer.terminalEof) {
        result[runId] = buffer.cursor;
      }
    });
    return result;
  }

  function schedulePoll(delay, replace) {
    if (pollTimer !== null) {
      if (!replace) return;
      window.clearTimeout(pollTimer);
    }
    pollTimer = window.setTimeout(function () {
      pollTimer = null;
      pollTelemetry();
    }, Math.max(0, delay));
  }

  async function pollTelemetry() {
    if (requestActive) return;
    if (!runTabVisible()) {
      schedulePoll(HIDDEN_CHECK_MS, false);
      return;
    }
    requestActive = true;
    let payload = null;
    try {
      const response = await fetch(TELEMETRY_URL, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        cache: "no-store",
        body: JSON.stringify({
          known_revision: knownRevision,
          log_cursors: logCursors(),
        }),
      });
      if (response.ok) payload = await response.json();
    } catch (_error) {
      // A later poll retries transient Dash reloads and broker replacements.
    } finally {
      requestActive = false;
    }
    if (payload) {
      knownRevision = Number(payload.revision) || 0;
      if (Array.isArray(payload.runs)) {
        runs.clear();
        payload.runs.forEach(function (run) {
          if (run && run.run_id) runs.set(String(run.run_id), run);
        });
      }
      Object.entries(payload.progress || {}).forEach(function (item) {
        observeProgress(item[0], Number(item[1].iteration), Number(item[1].total));
      });
      Object.entries(payload.logs || {}).forEach(function (item) {
        appendLog(item[0], item[1]);
      });
      renderRuns();
      if ((payload.deferred_logs || []).length) {
        schedulePoll(0, false);
        return;
      }
    }
    const active = Array.from(runs.values()).some(function (run) {
      return ACTIVE_STATES.has(run.state);
    });
    schedulePoll(active ? ACTIVE_POLL_MS : IDLE_POLL_MS, false);
  }

  window.dashboardRunConsoleState = window.dashboardRunConsoleState || {};
  window.dash_clientside = window.dash_clientside || {};
  window.dash_clientside.dashboardRunConsoleState = window.dashboardRunConsoleState;
  window.dashboardRunConsoleState.synchronizeInputs = function (selection, actionResult, outputDirectory) {
    selectedCases = new Set(selection || []);
    targetOutputDirectory = outputDirectory ? String(outputDirectory) : null;
    updateCaseButtons();
    if (actionResult) {
      submitting = false;
      const action = String(actionResult.action || "");
      if (action === "clear") {
        transientErrors.clear();
        buffers.clear();
        progressByRun.clear();
        openRuns.clear();
        knownRevision = null;
      } else if (action === "run") {
        transientErrors.clear();
      } else if (action === "error") {
        const at = Number(actionResult.at) || Date.now() / 1000;
        (actionResult.cases || ["submission"]).forEach(function (caseName, index) {
          const runId = `submission-error-${at}-${index}`;
          transientErrors.set(runId, {
            run_id: runId,
            case: String(caseName),
            state: "error",
            origin: "dash",
            message: String(actionResult.message || "Run submission failed"),
            updated_at: at,
          });
        });
      }
      renderRuns();
      schedulePoll(0, true);
    } else {
      updateControls();
    }
    renderSequence += 1;
    return renderSequence;
  };

  document.addEventListener("click", function (event) {
    const button = event.target && event.target.closest
      ? event.target.closest("#run-button")
      : null;
    if (button && !button.disabled) {
      submitting = true;
      updateControls();
    }
  }, true);

  function init() {
    updateControls();
    updateCaseButtons();
    schedulePoll(0, true);
    window.setInterval(updateDisplays, 1000);
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init, { once: true });
  } else {
    init();
  }
})();
