(function () {
  "use strict";

  const GEAR_FRAME_INTERVAL_MS = 160;
  const GEAR_SIZE = 32;
  let gearTimer = null;
  let gearAngle = 0;
  let gearLink = null;
  let idleFavicon = null;

  function stopGear() {
    if (gearTimer !== null) {
      window.clearInterval(gearTimer);
      gearTimer = null;
    }
  }

  function normalizedCount(value) {
    const count = Number(value);
    return Number.isFinite(count) && count > 0 ? Math.floor(count) : 0;
  }

  function faviconLink() {
    if (gearLink && gearLink.isConnected) return gearLink;

    const existingLinks = Array.from(document.querySelectorAll('link[rel~="icon"]'));
    gearLink = existingLinks.shift() || document.createElement("link");
    idleFavicon = {
      href: gearLink.getAttribute("href"),
      rel: gearLink.getAttribute("rel"),
      type: gearLink.getAttribute("type"),
      sizes: gearLink.getAttribute("sizes"),
    };
    existingLinks.forEach(function (link) {
      link.remove();
    });
    gearLink.id = "dashboard-job-favicon";
    if (!gearLink.isConnected) {
      document.head.appendChild(gearLink);
    }
    return gearLink;
  }

  function setOptionalAttribute(element, name, value) {
    if (value === null || value === undefined) {
      element.removeAttribute(name);
    } else {
      element.setAttribute(name, value);
    }
  }

  function restoreIdleFavicon() {
    const link = faviconLink();
    setOptionalAttribute(link, "href", idleFavicon.href);
    setOptionalAttribute(link, "rel", idleFavicon.rel || "icon");
    setOptionalAttribute(link, "type", idleFavicon.type);
    setOptionalAttribute(link, "sizes", idleFavicon.sizes);
  }

  function showGearFrame(dataUrl) {
    const link = faviconLink();
    link.rel = "icon";
    link.type = "image/png";
    link.setAttribute("sizes", "32x32");
    link.href = dataUrl;
  }

  function drawGear(angle, active) {
    const canvas = document.createElement("canvas");
    canvas.width = GEAR_SIZE;
    canvas.height = GEAR_SIZE;
    const context = canvas.getContext("2d");
    if (!context) return;

    const fill = active ? "#22c55e" : "#64748b";
    const edge = active ? "#15803d" : "#334155";
    context.translate(GEAR_SIZE / 2, GEAR_SIZE / 2);
    context.rotate(angle);
    context.shadowColor = active ? "rgba(34, 197, 94, 0.55)" : "rgba(15, 23, 42, 0.35)";
    context.shadowBlur = active ? 2.5 : 1.5;

    context.fillStyle = fill;
    context.strokeStyle = edge;
    context.lineWidth = 1.25;
    for (let tooth = 0; tooth < 10; tooth += 1) {
      context.save();
      context.rotate((tooth * Math.PI * 2) / 10);
      context.beginPath();
      context.rect(-2.7, -15, 5.4, 7);
      context.fill();
      context.stroke();
      context.restore();
    }

    context.beginPath();
    context.arc(0, 0, 10.5, 0, Math.PI * 2);
    context.fill();
    context.stroke();
    context.shadowBlur = 0;

    context.globalCompositeOperation = "destination-out";
    context.beginPath();
    context.arc(0, 0, 4.2, 0, Math.PI * 2);
    context.fill();
    context.globalCompositeOperation = "source-over";
    context.beginPath();
    context.arc(0, 0, 4.2, 0, Math.PI * 2);
    context.strokeStyle = active ? "#86efac" : "#94a3b8";
    context.lineWidth = 1.3;
    context.stroke();

    showGearFrame(canvas.toDataURL("image/png"));
  }

  function showIdleFavicon() {
    stopGear();
    gearAngle = 0;
    restoreIdleFavicon();
  }

  function startActiveGear() {
    stopGear();
    gearAngle = 0;
    drawGear(gearAngle, true);
    gearTimer = window.setInterval(function () {
      gearAngle = (gearAngle + Math.PI / 12) % (Math.PI * 2);
      drawGear(gearAngle, true);
    }, GEAR_FRAME_INTERVAL_MS);
  }

  window.dash_clientside = Object.assign({}, window.dash_clientside, {
    dashboardTitle: {
      sync: function (activeCount, baseTitle) {
        const title = String(baseTitle || "CLUBB Dash");
        const count = normalizedCount(activeCount);

        if (!count) {
          document.title = title;
          showIdleFavicon();
          return { active_count: 0 };
        }

        document.title = title + " - " + count;
        startActiveGear();
        return { active_count: count };
      },
    },
  });

  window.addEventListener("beforeunload", stopGear);
})();
