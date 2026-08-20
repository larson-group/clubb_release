(function () {
  "use strict";

  const CONSOLE_ID = "compile-console";
  const NEAR_BOTTOM_PX = 32;
  let consoleNode = null;
  let observer = null;
  let stick = true;
  let framePending = false;

  function isNearBottom(node) {
    return node.scrollHeight - node.scrollTop - node.clientHeight <= NEAR_BOTTOM_PX;
  }

  function scrollIfNeeded() {
    framePending = false;
    if (consoleNode && (stick || isNearBottom(consoleNode))) {
      consoleNode.scrollTop = consoleNode.scrollHeight;
      stick = true;
    }
  }

  function scheduleScroll() {
    if (framePending) return;
    framePending = true;
    window.requestAnimationFrame(scrollIfNeeded);
  }

  function attach() {
    const nextNode = document.getElementById(CONSOLE_ID);
    if (nextNode === consoleNode) return;
    if (observer) observer.disconnect();
    consoleNode = nextNode;
    if (!consoleNode) return;
    stick = true;
    consoleNode.addEventListener("scroll", function () {
      stick = isNearBottom(consoleNode);
    }, { passive: true });
    observer = new MutationObserver(scheduleScroll);
    observer.observe(consoleNode, {
      childList: true,
      subtree: true,
      characterData: true,
    });
    scheduleScroll();
  }

  function init() {
    attach();
    window.setInterval(attach, 1000);
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init, { once: true });
  } else {
    init();
  }
})();
