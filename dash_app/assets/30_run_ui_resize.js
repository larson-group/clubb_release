(function () {
  "use strict";

  function getLayout() {
    return document.getElementById("run-tab-layout");
  }

  function getDivider() {
    return document.getElementById("run-pane-divider");
  }

  function getCurrentRightWidth(layout) {
    if (!layout) return 360;
    var tpl = layout.style.gridTemplateColumns || "";
    var match = tpl.match(/(\d+)px\s*$/);
    if (match) {
      return Math.max(0, Number(match[1]) || 360);
    }
    return 360;
  }

  function saveWidth(width) {
    try {
      window.localStorage.setItem("run-ui-right-width", String(Math.round(width)));
    } catch (e) {}
  }

  function loadWidth() {
    try {
      var raw = window.localStorage.getItem("run-ui-right-width");
      if (!raw) return null;
      var num = Number(raw);
      if (!Number.isFinite(num)) return null;
      return num;
    } catch (e) {
      return null;
    }
  }

  function maxRightWidth(layout) {
    var layoutW = layout ? layout.getBoundingClientRect().width : window.innerWidth;
    return Math.max(300, Math.min(Math.round(window.innerWidth * 0.8), Math.round(layoutW - 360)));
  }

  function clampRightWidth(layout, width) {
    var minW = 300;
    var maxW = maxRightWidth(layout);
    var rightW = Number(width) || minW;
    if (rightW < minW) rightW = minW;
    if (rightW > maxW) rightW = maxW;
    return rightW;
  }

  function applyRightWidth(layout, width) {
    layout.style.gridTemplateColumns = "minmax(0,1fr) 8px " + Math.round(width) + "px";
  }

  function bindDivider() {
    var divider = getDivider();
    var layout = getLayout();
    if (!divider || !layout || divider.dataset.bound === "1") return;
    divider.dataset.bound = "1";

    var dragging = false;

    function startDrag(ev) {
      if (ev.button !== 0) return;
      dragging = true;
      document.body.classList.add("run-ui-dragging");
      ev.preventDefault();
    }

    function moveDrag(ev) {
      if (!dragging) return;
      var rect = layout.getBoundingClientRect();
      var rightW = clampRightWidth(layout, rect.right - ev.clientX);
      applyRightWidth(layout, rightW);
    }

    function endDrag(ev) {
      if (!dragging) return;
      dragging = false;
      document.body.classList.remove("run-ui-dragging");
      var rect = layout.getBoundingClientRect();
      var rightW = clampRightWidth(layout, rect.right - ev.clientX);
      applyRightWidth(layout, rightW);
      saveWidth(rightW);
    }

    divider.addEventListener("pointerdown", startDrag);
    window.addEventListener("pointermove", moveDrag);
    window.addEventListener("pointerup", endDrag);
  }

  function hydrateWidth() {
    var layout = getLayout();
    if (!layout || layout.dataset.widthHydrated === "1") return;
    layout.dataset.widthHydrated = "1";
    var saved = loadWidth();
    if (saved === null) return;
    applyRightWidth(layout, clampRightWidth(layout, saved));
  }

  function tick() {
    bindDivider();
    hydrateWidth();
  }

  document.addEventListener("DOMContentLoaded", tick);
  window.setInterval(tick, 500);
})();
