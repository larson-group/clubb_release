document.addEventListener("keydown", function (event) {
    if (event.key !== "Escape") {
        return;
    }
    const popover = document.getElementById("compile-build-selector-popover");
    const dismiss = document.getElementById("compile-build-selector-dismiss");
    if (popover && dismiss && !popover.classList.contains("compile-build-selector-popover-closed")) {
        dismiss.click();
    }
});
