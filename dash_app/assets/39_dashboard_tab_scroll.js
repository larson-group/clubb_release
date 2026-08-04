/* Preserve each top-level dashboard tab's reading position across tab remounts. */
(function () {
  "use strict";

  const TAB_SELECTOR = "#dashboard-tabs .tab--selected";
  const PANE_IDS = ["run-right-pane", "tune-right-pane", "plots-right-pane"];
  const scrollByTab = new Map();
  let activeTab = null;
  let restoreScheduled = false;

  function selectedTabKey() {
    const tab = document.querySelector(TAB_SELECTOR);
    if (!tab) return null;
    return tab.getAttribute("value") || tab.getAttribute("data-value") || tab.textContent.trim();
  }

  function stateFor(tabKey) {
    if (!scrollByTab.has(tabKey)) {
      scrollByTab.set(tabKey, { windowY: 0, panes: {} });
    }
    return scrollByTab.get(tabKey);
  }

  function remember(tabKey) {
    if (!tabKey) return;
    const state = stateFor(tabKey);
    state.windowY = window.scrollY;
    PANE_IDS.forEach(function (paneId) {
      const pane = document.getElementById(paneId);
      if (pane) state.panes[paneId] = pane.scrollTop;
    });
  }

  function restore(tabKey) {
    const state = scrollByTab.get(tabKey);
    if (!state) return;
    window.scrollTo(window.scrollX, state.windowY || 0);
    Object.keys(state.panes).forEach(function (paneId) {
      const pane = document.getElementById(paneId);
      if (!pane) return;
      const maximum = Math.max(0, pane.scrollHeight - pane.clientHeight);
      pane.scrollTop = Math.min(maximum, state.panes[paneId]);
    });
  }

  function scheduleRestore(tabKey) {
    if (restoreScheduled) return;
    restoreScheduled = true;
    window.requestAnimationFrame(function () {
      restoreScheduled = false;
      restore(tabKey);
      // Dash may mount a selected tab's children in a follow-up React pass.
      window.setTimeout(function () { restore(tabKey); }, 0);
    });
  }

  function syncActiveTab() {
    const nextTab = selectedTabKey();
    if (!nextTab) return;
    if (activeTab && activeTab !== nextTab && !scrollByTab.has(activeTab)) {
      // Pointer/click capture normally saved the old tab first.  If navigation
      // was programmatic before the tab was ever scrolled, capture its initial
      // position without overwriting a previously saved reading position.
      remember(activeTab);
    }
    if (activeTab !== nextTab) {
      activeTab = nextTab;
      scheduleRestore(nextTab);
    }
  }

  function rememberBeforeTabChange(event) {
    if (event.target && event.target.closest && event.target.closest("#dashboard-tabs .tab")) {
      remember(activeTab || selectedTabKey());
    }
  }

  document.addEventListener("pointerdown", rememberBeforeTabChange, true);
  document.addEventListener("click", rememberBeforeTabChange, true);
  window.addEventListener("scroll", function () { remember(activeTab || selectedTabKey()); }, { passive: true });
  document.addEventListener("scroll", function (event) {
    const target = event.target;
    if (!target || !target.id || PANE_IDS.indexOf(target.id) === -1) return;
    const tabKey = activeTab || selectedTabKey();
    if (!tabKey) return;
    stateFor(tabKey).panes[target.id] = target.scrollTop;
  }, true);

  let mutationScheduled = false;
  function scheduleTabSync() {
    if (mutationScheduled) return;
    mutationScheduled = true;
    window.requestAnimationFrame(function () {
      mutationScheduled = false;
      syncActiveTab();
    });
  }

  function observeTabs() {
    const tabs = document.getElementById("dashboard-tabs");
    if (!tabs) {
      window.setTimeout(observeTabs, 50);
      return;
    }
    new MutationObserver(scheduleTabSync).observe(tabs, {
      attributes: true,
      subtree: true,
      attributeFilter: ["class"],
    });
    syncActiveTab();
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", observeTabs);
  } else {
    observeTabs();
  }
}());
