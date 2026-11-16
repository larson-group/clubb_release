(function () {
  "use strict";

  // Use Markdown's rendered anchors, including reference links and nested labels.
  // Resolve before the browser follows the link so normal modifier/middle clicks
  // and the context menu keep their native new-tab behavior.
  function resolveDocumentLink(event) {
    const anchor = event.target.closest && event.target.closest(".reports-markdown a[href]");
    const documentView = anchor && anchor.closest("[data-repository-path]");
    if (!documentView) return;
    const href = anchor.getAttribute("href");
    if (!href || href.startsWith("//") || /^[a-z][a-z0-9+.-]*:/i.test(href)) return;
    if (["/repository-files/", "/documents/", "/document-images/"].some(prefix => href.startsWith(prefix))) return;
    const base = new URL("/repository-files/" +
      documentView.dataset.repositoryPath.split("/").map(encodeURIComponent).join("/"), window.location.origin);
    const target = new URL(href.startsWith("/") ? "/repository-files" + href : href, base);
    if (!target.pathname.startsWith("/repository-files/")) return;
    anchor.setAttribute("href", target.href);
    anchor.setAttribute("rel", "noopener noreferrer");
  }

  ["pointerdown", "focusin", "click", "auxclick", "contextmenu"].forEach(function (name) {
    document.addEventListener(name, resolveDocumentLink, true);
  });
})();
