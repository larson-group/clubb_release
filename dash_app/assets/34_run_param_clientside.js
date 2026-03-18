(function() {
  function normalize(rawValue) {
    const text = rawValue === null || rawValue === undefined ? '' : String(rawValue).trim();
    if (!text) {
      return text;
    }

    const lower = text.toLowerCase();
    const boolTrue = new Set(['.true.', 'true', 't']);
    const boolFalse = new Set(['.false.', 'false', 'f']);
    if (boolTrue.has(lower) || boolFalse.has(lower)) {
      return text;
    }
    if (lower.includes('e')) {
      return text;
    }
    if (!text.includes('.')) {
      return text;
    }

    let trimmed = text.replace(/0+$/, '');
    if (trimmed.endsWith('.')) {
      trimmed += '0';
    }
    return trimmed;
  }

  const clientside = window.dash_clientside = window.dash_clientside || {};
  clientside.runTab = clientside.runTab || {};

  clientside.runTab.syncMulticolDisabled = function(multicolParamValues, tunableIds, currentDisabled) {
    var claimed = new Set();
    (multicolParamValues || []).forEach(function(v) {
      if (v !== null && v !== undefined) {
        var s = String(v).trim();
        if (s) claimed.add(s);
      }
    });
    return (tunableIds || []).map(function(id, i) {
      var shouldDisable = claimed.has((id || {}).name);
      if (shouldDisable === !!(currentDisabled || [])[i]) {
        return window.dash_clientside.no_update;
      }
      return shouldDisable;
    });
  };

  clientside.runTab.syncParamRowClass = function(value, disabled, inputId, defaultsByKey) {
    const key = inputId && inputId.file && inputId.name
      ? `${inputId.file}:${inputId.name}`
      : '';
    const defaults = defaultsByKey || {};
    const currentValue = normalize(value);
    const defaultValue = normalize(defaults[key]);
    const classes = ['run-param-container'];

    if (disabled) {
      classes.push('run-param-row--changed-disabled');
    } else if (currentValue !== defaultValue) {
      classes.push('run-param-row--changed');
    } else {
      classes.push('run-param-row--default');
    }

    return classes.join(' ');
  };
})();
