(function () {
  var state = new WeakMap();
  var stickById = new Map();
  var rafPending = false;
  var NEAR_BOTTOM_PX = 32;

  function isNearBottom(el) {
    return el.scrollHeight - el.scrollTop - el.clientHeight <= NEAR_BOTTOM_PX;
  }

  function ensureState(el) {
    var entry = state.get(el);
    if (entry) {
      return entry;
    }
    var id = el.id || null;
    var stick = true;
    if (id && stickById.has(id)) {
      stick = stickById.get(id);
    }
    entry = { stick: stick, id: id };
    state.set(el, entry);
    el.addEventListener(
      "scroll",
      function (evt) {
        var target = evt && evt.target ? evt.target : el;
        entry.stick = isNearBottom(target);
        if (entry.id) {
          stickById.set(entry.id, entry.stick);
        }
      },
      { passive: true }
    );
    return entry;
  }

  function maybeScroll(el) {
    if (!el) {
      return;
    }
    var entry = ensureState(el);
    if (entry.stick || isNearBottom(el)) {
      el.scrollTop = el.scrollHeight;
      entry.stick = true;
      if (entry.id) {
        stickById.set(entry.id, true);
      }
    }
  }

  function tick() {
    rafPending = false;
    var nodes = document.querySelectorAll(".run-console-active");
    for (var i = 0; i < nodes.length; i += 1) {
      var el = nodes[i];
      if (!el) {
        continue;
      }
      maybeScroll(el);
    }
  }

  function scheduleTick() {
    if (rafPending) {
      return;
    }
    rafPending = true;
    window.requestAnimationFrame(tick);
  }

  function observe() {
    var container =
      document.getElementById("run-console-container") || document.body;
    if (!container) {
      setTimeout(observe, 500);
      return;
    }
    var observer = new MutationObserver(function () {
      scheduleTick();
    });
    observer.observe(container, {
      childList: true,
      subtree: true,
      characterData: true,
    });
    scheduleTick();
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", observe);
  } else {
    observe();
  }
})();
