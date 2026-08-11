(function () {
  "use strict";

  function noUpdate() {
    return window.dash_clientside.no_update;
  }

  function ownedKey(key, meta) {
    var exact = (meta.store_ids || []).concat(meta.extra_keys || []);
    exact.push(meta.schema_key);
    return exact.indexOf(key) >= 0 || key.indexOf(String(meta.token || "")) >= 0;
  }

  function ownedKeys(meta) {
    var keys = [];
    for (var index = 0; index < window.localStorage.length; index += 1) {
      var key = window.localStorage.key(index);
      if (key && ownedKey(key, meta || {})) {
        keys.push(key);
      }
    }
    return keys;
  }

  function clearRunControlPersistence() {
    // Dash stores automatically persisted component values under implementation
    // keys that include the configured persistence token and serialized ID.
    // The exact key prefix is Dash-version dependent, so identify only the
    // dynamic Run controls by their stable pattern-ID fragments.  Do not clear
    // Run Stores such as selected cases or selected config: those are updated
    // by their owning Dash callbacks and should persist normally.
    var markers = [
      'run-param',
      'run-linked-param',
      'run-flag',
      'run-hr-param',
      'run-hr-min',
      'run-hr-max',
      'run-hr-npoints'
    ];
    var keys = [];
    for (var index = 0; index < window.localStorage.length; index += 1) {
      var key = window.localStorage.key(index);
      if (key && markers.some(function (marker) { return key.indexOf(marker) >= 0; })) {
        keys.push(key);
      }
    }
    keys.forEach(function (key) { window.localStorage.removeItem(key); });
  }

  function clearRunTransientPersistence() {
    // These Stores are normally kept in local storage by Dash.  Clear is an
    // explicit request to discard their UI state, so remove the durable copy
    // as well as letting the Run callback clear the live values.  This avoids
    // a stale selected-case list being restored during the next page load.
    ["run-selected-cases", "run-open-cases"].forEach(function (key) {
      window.localStorage.removeItem(key);
    });
    var consoleState = window.dashboardRunConsoleState;
    if (consoleState && typeof consoleState.clearOpenCases === "function") {
      consoleState.clearOpenCases();
    }
  }

  function decodeUpload(contents) {
    var encoded = String(contents || "").split(",", 2)[1];
    if (!encoded) {
      throw new Error("The selected file has no JSON payload.");
    }
    var binary = window.atob(encoded);
    var bytes = new Uint8Array(binary.length);
    for (var index = 0; index < binary.length; index += 1) {
      bytes[index] = binary.charCodeAt(index);
    }
    return new TextDecoder("utf-8").decode(bytes);
  }

  window.dash_clientside = Object.assign({}, window.dash_clientside, {
    dashboardWorkspace: {
      initializeWorkspace: function (meta) {
        if (!meta || !meta.schema_key) {
          return noUpdate();
        }
        var expected = String(meta.schema_version);
        var current = window.localStorage.getItem(meta.schema_key);
        if (current !== null && current !== expected) {
          ownedKeys(meta).forEach(function (key) {
            window.localStorage.removeItem(key);
          });
          window.localStorage.setItem(meta.schema_key, expected);
          window.setTimeout(function () { window.location.reload(); }, 100);
          return "Discarded an incompatible saved workspace.";
        }
        window.localStorage.setItem(meta.schema_key, expected);
        return "Autosave active (schema " + expected + ").";
      },

      exportWorkspace: function (clicks, meta) {
        if (!clicks) {
          return noUpdate();
        }
        var entries = {};
        ownedKeys(meta || {}).forEach(function (key) {
          entries[key] = window.localStorage.getItem(key);
        });
        var payload = {
          schema_version: meta.schema_version,
          token: meta.token,
          repo_name: meta.repo_name,
          exported_at: new Date().toISOString(),
          local_storage: entries
        };
        var blob = new Blob([JSON.stringify(payload, null, 2) + "\n"], {
          type: "application/json"
        });
        var link = document.createElement("a");
        link.href = URL.createObjectURL(blob);
        link.download = String(meta.repo_name || "clubb") + "_dashboard_workspace.json";
        document.body.appendChild(link);
        link.click();
        document.body.removeChild(link);
        URL.revokeObjectURL(link.href);
        return "Exported " + Object.keys(entries).length + " saved entries.";
      },

      importWorkspace: function (contents, filename, meta) {
        if (!contents) {
          return noUpdate();
        }
        try {
          var payload = JSON.parse(decodeUpload(contents));
          if (Number(payload.schema_version) !== Number(meta.schema_version)) {
            throw new Error(
              "Workspace schema " + payload.schema_version +
              " is incompatible with schema " + meta.schema_version + "."
            );
          }
          ownedKeys(meta || {}).forEach(function (key) {
            window.localStorage.removeItem(key);
          });
          var entries = payload.local_storage || {};
          var restored = 0;
          Object.keys(entries).forEach(function (key) {
            if (ownedKey(key, meta || {}) && typeof entries[key] === "string") {
              window.localStorage.setItem(key, entries[key]);
              restored += 1;
            }
          });
          window.setTimeout(function () { window.location.reload(); }, 100);
          return "Restored " + restored + " entries from " + (filename || "workspace") + ".";
        } catch (error) {
          return "Import failed: " + error.message;
        }
      },

      resetWorkspace: function (clicks, meta) {
        if (!clicks) {
          return noUpdate();
        }
        var keys = ownedKeys(meta || {});
        keys.forEach(function (key) { window.localStorage.removeItem(key); });
        window.setTimeout(function () { window.location.reload(); }, 100);
        return "Cleared " + keys.length + " saved entries.";
      },

      selectSavedRunConfig: function (selection) {
        if (!selection || !selection.name) {
          return noUpdate();
        }
        clearRunControlPersistence();
        return selection;
      },

      resetRunConfigControls: function (_clicks, buttonIds) {
        var context = window.dash_clientside.callback_context || {};
        var triggered = context.triggered_id;
        if (!triggered || triggered.type !== 'run-config-button' || !triggered.name) {
          return noUpdate();
        }
        var triggeredIndex = (buttonIds || []).findIndex(function (buttonId) {
          return buttonId &&
            buttonId.type === triggered.type &&
            buttonId.name === triggered.name;
        });
        if (triggeredIndex < 0 || Number((_clicks || [])[triggeredIndex] || 0) < 1) {
          // Adding/rebuilding config buttons invokes an ALL-pattern callback
          // with zero clicks. It is layout hydration, not user selection.
          return noUpdate();
        }
        clearRunControlPersistence();
        return { name: String(triggered.name), nonce: Date.now() };
      },

      clearRunTransientState: function (clicks) {
        if (!clicks) {
          return noUpdate();
        }
        clearRunTransientPersistence();
        return { nonce: Date.now() };
      }
    }
  });
}());
