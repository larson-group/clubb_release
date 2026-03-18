(function () {
  "use strict";

  function getLayout() {
    return document.getElementById("plots-tab-layout");
  }

  function getDivider() {
    return document.getElementById("plots-pane-divider");
  }

  function saveWidth(width) {
    try {
      window.localStorage.setItem("plots-ui-right-width", String(Math.round(width)));
    } catch (e) {}
  }

  function notifyResize() {
    try {
      window.dispatchEvent(new Event("resize"));
    } catch (e) {}
  }

  function loadWidth() {
    try {
      var raw = window.localStorage.getItem("plots-ui-right-width");
      if (!raw) return null;
      var num = Number(raw);
      if (!Number.isFinite(num)) return null;
      return num;
    } catch (e) {
      return null;
    }
  }

  function bindDivider() {
    var divider = getDivider();
    var layout = getLayout();
    if (!divider || !layout || divider.dataset.bound === "1") return;
    divider.dataset.bound = "1";

    var dragging = false;
    var minW = 260;
    var maxW = Math.max(420, Math.round(window.innerWidth * 0.8));

    function startDrag(ev) {
      if (ev.button !== 0) return;
      maxW = Math.max(420, Math.round(window.innerWidth * 0.8));
      dragging = true;
      document.body.classList.add("plots-ui-dragging");
      ev.preventDefault();
    }

    function moveDrag(ev) {
      if (!dragging) return;
      var rect = layout.getBoundingClientRect();
      var rightW = rect.right - ev.clientX;
      if (rightW < minW) rightW = minW;
      if (rightW > maxW) rightW = maxW;
      layout.style.gridTemplateColumns = "minmax(0,1fr) 8px " + Math.round(rightW) + "px";
    }

    function endDrag(ev) {
      if (!dragging) return;
      dragging = false;
      document.body.classList.remove("plots-ui-dragging");
      var rect = layout.getBoundingClientRect();
      var rightW = rect.right - ev.clientX;
      if (rightW < minW) rightW = minW;
      if (rightW > maxW) rightW = maxW;
      layout.style.gridTemplateColumns = "minmax(0,1fr) 8px " + Math.round(rightW) + "px";
      saveWidth(rightW);
      notifyResize();
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
    var minW = 260;
    var maxW = Math.max(420, Math.round(window.innerWidth * 0.8));
    var rightW = saved;
    if (rightW < minW) rightW = minW;
    if (rightW > maxW) rightW = maxW;
    layout.style.gridTemplateColumns = "minmax(0,1fr) 8px " + Math.round(rightW) + "px";
    window.setTimeout(notifyResize, 0);
  }

  function tick() {
    bindDivider();
    hydrateWidth();
  }

  document.addEventListener("DOMContentLoaded", tick);
  window.setInterval(tick, 500);
})();
