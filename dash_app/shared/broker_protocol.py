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
# Runtime source fingerprints are metadata/lifecycle policy, not an API
# contract, so they do not require a protocol bump.
BROKER_PROTOCOL_VERSION = 20
