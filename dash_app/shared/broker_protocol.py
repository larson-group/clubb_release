"""Small dashboard-runtime broker compatibility contract."""

# Bump this integer when the Dash source requires a new internal broker action
# or changes a broker response contract.  ``ensure_broker`` replaces an idle
# lower-version sidecar automatically and preserves an active one.
# Version 13 moves all typed mutating MCP operations into the detached broker.
# A reload must replace an idle older sidecar so watcher threads cannot remain
# owned by a short-lived agent adapter process.
# Version 14 moves the Streamable HTTP endpoint under the durable broker.
# Version 15 adds the typed SCM batch submission action used by the Run tab.
# Version 16 lets native Run batches honor the selected output directory.
# Version 17 lets public MCP SCM requests select a directory below output/.
# Version 18 lets MCP child runs use any directory below output/.
# Version 19 fixes artifact retention ordering after active bundles are released.
# Version 20 adds manager-owned graceful shutdown for every broker job family.
# Version 21 adds the broker-owned Profile timing job used by the Profile tab.
# Version 22 adds atomic SCM-batch cancellation and browser-visible queued children.
# Version 23 makes SCM state compact, session-scoped, and directly log-streamed.
# Version 24 adds broker-owned cancellation of every Run-tab SCM job.
# Version 25 adds broker-owned PyPlotGen exports from the Plot tab.
# Version 26 launches PyPlotGen unbuffered for live progress reporting.
# Version 27 adds immediate broker-owned PyPlotGen cancellation.
# Version 28 propagates unbuffered stdout into PyPlotGen pool workers.
# Version 29 scopes SCM admission to canonical case/output pairs and returns
# accepted/skipped batch cases without public output-directory manifests.
# Version 30 serializes artifact retention and keeps SCM batch queues alive
# after an individual child setup failure.
# Runtime source fingerprints are metadata/lifecycle policy, not an API
# contract, so they do not require a protocol bump.
# Version 31 adds bounded native Plot tasks in broker-owned worker processes.
# Version 32 batches native Plot result polling across cards.
# Version 33 adds a dedicated lane for the short gradient response plot.
# Version 34 reports live gradient-response stages and elapsed times.
# Version 35 adds the case and averaging window to gradient-response requests.
# Version 36 evaluates gradient-response samples together in 16 columns.
# Version 37 adds the independent multi-field gradient map and response lane.
# Version 38 replaces the old gradient response task with map-anchored sweeps.
# Version 39 carries entered parameter ranges through gradient maps and sweeps.
# Version 40 adds JAX loss-gradient requests and persistent explorer loss sessions.
# Version 41 moves explorer JAX process ownership into Dash runtime executors.
# Version 42 adds all-sample gradient overlays using the retained 16-column worker.
# Version 43 carries automatic/forward/reverse loss-gradient requests.
# Version 45 carries fixed Taylor objective weights and per-case loss gradients.
BROKER_PROTOCOL_VERSION = 45
