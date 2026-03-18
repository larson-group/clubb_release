(function() {
  const GRID_SELECTOR = '.run-param-list';
  const ROW_SELECTOR = '.run-param-container';

  function colorFor(themeDark, stripeA) {
    if (themeDark) {
      return stripeA ? 'rgba(148, 163, 184, 0.08)' : 'rgba(148, 163, 184, 0.03)';
    }
    return stripeA ? 'rgba(148, 163, 184, 0.14)' : 'rgba(148, 163, 184, 0.06)';
  }

  function detectColumnCount(rows) {
    if (!rows.length) {
      return 1;
    }
    const firstTop = rows[0].offsetTop;
    for (let idx = 1; idx < rows.length; idx += 1) {
      if (rows[idx].offsetTop !== firstTop) {
        return idx;
      }
    }
    return rows.length;
  }

  function updateGrid(grid) {
    const rows = Array.from(grid.querySelectorAll(ROW_SELECTOR));
    if (!rows.length) {
      return;
    }
    const appRoot = document.getElementById('app-root');
    const themeDark = !!(appRoot && appRoot.classList.contains('theme-dark'));
    const ncols = Math.max(detectColumnCount(rows), 1);
    rows.forEach((row, idx) => {
      const stripeA = Math.floor(idx / ncols) % 2 === 0;
      row.style.backgroundColor = colorFor(themeDark, stripeA);
      row.style.borderRadius = '4px';
    });
  }

  function updateAll() {
    document.querySelectorAll(GRID_SELECTOR).forEach(updateGrid);
  }

  let raf = null;
  function scheduleUpdate() {
    if (raf !== null) {
      return;
    }
    raf = window.requestAnimationFrame(() => {
      raf = null;
      updateAll();
    });
  }

  const observer = new MutationObserver(scheduleUpdate);

  function init() {
    observer.observe(document.body, {
      childList: true,
      subtree: true,
      attributes: true,
      attributeFilter: ['class', 'style'],
    });
    window.addEventListener('resize', scheduleUpdate);
    scheduleUpdate();
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init, { once: true });
  } else {
    init();
  }
})();
