(function () {
  "use strict";

  function collapseOutputMenu() {
    if (!window.dash_clientside || !window.dash_clientside.set_props) return;
    window.dash_clientside.set_props("plots-output-menu-expanded", { data: false });
  }

  document.addEventListener("click", function (event) {
    const menu = document.getElementById("plots-output-menu");
    if (!menu || !menu.classList.contains("plots-output-menu--expanded")) return;

    const target = event.target;
    if (target && target.closest && target.closest(".plots-output-available-button")) {
      collapseOutputMenu();
      return;
    }
    if (!menu.contains(target)) collapseOutputMenu();
  });
})();
