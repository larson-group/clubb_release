(function () {
  "use strict";

  const PANEL_SELECTOR = ".run-console-panel[data-case-name]";
  const LOG_SELECTOR = ".run-console[id]";
  const CONTAINER_ID = "run-console-container";
  const STORE_ID = "run-open-cases";
  const NEAR_BOTTOM_PX = 32;
  let openCases = new Set();
  let logScrollById = new Map();
  let lastContainerScrollTop = null;
  let lastWindowScrollY = null;

  function storeNode() {
    return document.getElementById(STORE_ID);
  }

  function containerNode() {
    return document.getElementById(CONTAINER_ID);
  }

  function writeStore() {
    const store = storeNode();
    if (!store || !window.dash_clientside || !window.dash_clientside.set_props) return;
    window.dash_clientside.set_props(STORE_ID, { data: Array.from(openCases).sort() });
  }

  function rememberContainerScroll() {
    const container = containerNode();
    if (container) {
      lastContainerScrollTop = container.scrollTop;
    }
  }

  function rememberWindowScroll() {
    lastWindowScrollY = window.scrollY;
  }

  function restoreScroll() {
    const container = containerNode();
    if (container && lastContainerScrollTop !== null) {
      const maxScroll = Math.max(0, container.scrollHeight - container.clientHeight);
      container.scrollTop = Math.min(lastContainerScrollTop, maxScroll);
    }
    if (lastWindowScrollY !== null) {
      window.scrollTo(window.scrollX, lastWindowScrollY);
    }
    restoreLogScroll();
  }

  function logIsNearBottom(logNode) {
    return logNode.scrollHeight - logNode.scrollTop - logNode.clientHeight <= NEAR_BOTTOM_PX;
  }

  function rememberLogScroll(logNode) {
    if (!logNode || !logNode.id) return;
    logScrollById.set(logNode.id, {
      scrollTop: logNode.scrollTop,
      nearBottom: logIsNearBottom(logNode),
    });
  }

  function restoreLogScroll() {
    document.querySelectorAll(LOG_SELECTOR).forEach(function (logNode) {
      const saved = logScrollById.get(logNode.id);
      if (!saved || saved.nearBottom) return;
      const maxScroll = Math.max(0, logNode.scrollHeight - logNode.clientHeight);
      logNode.scrollTop = Math.min(saved.scrollTop, maxScroll);
    });
  }

  function syncPanelsFromStore() {
    let changed = false;
    document.querySelectorAll(PANEL_SELECTOR).forEach(function (panel) {
      const caseName = panel.getAttribute("data-case-name");
      if (!caseName) return;
      const isRunning = panel.getAttribute("data-case-running") === "true";
      if (isRunning && openCases.delete(caseName)) {
        changed = true;
      }
      if (!isRunning && openCases.has(caseName) && !panel.open) {
        panel.open = true;
      }
      if (!isRunning && !openCases.has(caseName) && panel.open) {
        panel.open = false;
        changed = true;
      }
    });
    if (changed) {
      writeStore();
    }
  }

  document.addEventListener(
    "toggle",
    function (event) {
      const panel = event.target;
      if (!panel || !panel.matches || !panel.matches(PANEL_SELECTOR)) return;
      const caseName = panel.getAttribute("data-case-name");
      if (!caseName) return;
      const isRunning = panel.getAttribute("data-case-running") === "true";
      if (panel.open && !isRunning) {
        openCases.add(caseName);
      } else {
        openCases.delete(caseName);
      }
      writeStore();
    },
    true
  );

  const observer = new MutationObserver(function () {
    syncPanelsFromStore();
    window.requestAnimationFrame(restoreScroll);
  });

  function init() {
    const container = containerNode();
    if (container) {
      container.addEventListener("scroll", rememberContainerScroll, { passive: true });
    }
    document.addEventListener(
      "scroll",
      function (event) {
        const target = event.target;
        if (target && target.matches && target.matches(LOG_SELECTOR)) {
          rememberLogScroll(target);
        }
      },
      true
    );
    window.addEventListener("scroll", rememberWindowScroll, { passive: true });
    observer.observe(document.body, { childList: true, subtree: true });
    syncPanelsFromStore();
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init, { once: true });
  } else {
    init();
  }
})();
