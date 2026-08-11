(function () {
  "use strict";

  var loadedGeneration = window.__CLUBB_DASH_GENERATION__ || null;
  var reloadStarted = false;

  function checkGeneration() {
    window.fetch("/_clubb-dashboard-generation", { cache: "no-store" })
      .then(function (response) {
        if (!response.ok) {
          throw new Error("Dashboard generation check failed");
        }
        return response.json();
      })
      .then(function (payload) {
        var current = payload && payload.generation;
        if (!current) {
          return;
        }
        if (loadedGeneration === null) {
          loadedGeneration = current;
          return;
        }
        if (current !== loadedGeneration && !reloadStarted) {
          reloadStarted = true;
          window.location.reload();
        }
      })
      .catch(function () {
        // A short outage is expected while the manager replaces Dash. The
        // next poll will compare the new process with the page we loaded.
      });
  }

  checkGeneration();
  window.setInterval(checkGeneration, 2500);
}());
