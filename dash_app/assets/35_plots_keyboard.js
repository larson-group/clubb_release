(function() {
  document.addEventListener('keydown', function(e) {
    if (e.key !== ' ') return;
    var tag = (document.activeElement || {}).tagName || '';
    if (tag === 'INPUT' || tag === 'TEXTAREA' || tag === 'SELECT') return;
    var btn = document.getElementById('plots-playback-toggle');
    if (btn && !btn.disabled) {
      e.preventDefault();
      btn.click();
    }
  });
})();
