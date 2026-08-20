(function () {
  "use strict";

  window.dash_clientside = window.dash_clientside || {};
  window.dash_clientside.plotsPyplotgen = {
    reserveWindow: function (clicks) {
      if (!clicks) {
        return window.dash_clientside.no_update;
      }
      const name = "clubb_pyplotgen_" + Date.now();
      return {name: name, available: false, click: clicks};
    },

    openCompleted: function (jobs, action, popup, openedRun) {
      const brokerJob = (jobs && jobs.pyplotgen) || {};
      const actionJob = (action && action.job) || {};
      const job = Number(actionJob.started_at || 0) > Number(brokerJob.started_at || 0) ? actionJob :
        (Object.keys(brokerJob).length ? brokerJob : actionJob);
      const runId = String(job.run_id || "");
      const actionRunId = String(actionJob.run_id || "");
      if (job.state !== "finished" || !job.html_url || !runId || runId === openedRun) {
        return window.dash_clientside.no_update;
      }
      if (!action || action.kind !== "started" || actionRunId !== runId) {
        return window.dash_clientside.no_update;
      }
      const target = popup && popup.available ? popup.name : "_blank";
      window.open(job.html_url, target);
      return runId;
    }
  };
}());
