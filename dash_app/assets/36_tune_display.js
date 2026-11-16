(function () {
  "use strict";
  let anchor = null;
  const plots = {
    taylor: "tune-taylor-diagram", parameters: "tune-parameter-box-plot",
    landscape: "tune-landscape-plot", correlation: "tune-parameter-correlation-plot",
    sensitivity: "tune-field-sensitivity-plot", interaction: "tune-field-interaction-plot"
  };
  function paintClock() {
    const node = document.getElementById("tune-runtime-uptime");
    if (!node || !anchor) return;
    const elapsed = Math.max(0, Math.floor(anchor.elapsed + (anchor.running ? (performance.now() - anchor.at) / 1000 : 0)));
    const hours = Math.floor(elapsed / 3600), minutes = Math.floor(elapsed / 60) % 60, seconds = elapsed % 60;
    node.textContent = hours ? hours + "h " + minutes + "m " + seconds + "s" : minutes ? minutes + "m " + seconds + "s" : seconds + "s";
  }
  window.clubbTuneDisplay = {
    clock: function (status) {
      status = status || {};
      const job = status.job_dir || "";
      const elapsed = Number(status.elapsed_seconds) || 0;
      const running = ["running", "initializing", "stopping"].includes(status.state);
      if (!anchor || anchor.job !== job || anchor.observed !== elapsed || anchor.running !== running) {
        anchor = {job: job, observed: elapsed, elapsed: elapsed, running: running, at: performance.now()};
      }
      paintClock();
      return window.dash_clientside.no_update;
    },
    visibility: function (enabled) {
      enabled = enabled || [];
      const set = window.dash_clientside.set_props;
      set("tune-taylor-container", {style: {display: enabled.some(key => ["taylor", "parameters"].includes(key)) ? "" : "none"}});
      set("tune-landscape-container", {style: {display: enabled.some(key => ["landscape", "correlation", "sensitivity", "interaction"].includes(key)) ? "" : "none"}});
      Object.entries(plots).forEach(function (pair) {
        set(pair[1], {style: {width: "100%", minWidth: 0, height: "430px", display: enabled.includes(pair[0]) ? "block" : "none"}});
        const node = document.getElementById(pair[1]);
        const card = node && node.closest(".tune-diagnostics-taylor, .tune-diagnostics-params");
        if (card) card.style.display = enabled.includes(pair[0]) ? "" : "none";
      });
      return window.dash_clientside.no_update;
    }
  };
  window.setInterval(paintClock, 500);
}());
