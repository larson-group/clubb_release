(function() {

  const clientside = window.dash_clientside = window.dash_clientside || {};
  clientside.runTab = clientside.runTab || {};

  function classForState(state, baseClass) {
    const classes = [baseClass];
    switch ((state || {}).state) {
      case 'invalid': classes.push('run-param-row--invalid'); break;
      case 'inactive': classes.push('run-param-row--inactive-disabled'); break;
      case 'linked-changed': classes.push('run-param-row--changed'); break;
      case 'linked': classes.push('run-param-row--default'); break;
      case 'changed': classes.push('run-param-row--changed'); break;
      default: classes.push('run-param-row--default'); break;
    }
    return classes.join(' ');
  }

  // The shared Python settings evaluator owns every semantic decision.  This
  // browser callback only maps its already-resolved state to CSS classes.
  clientside.runTab.syncAllParamRowClasses = function(resolution, multicolParams, inputIds) {
    const states = (resolution || {}).control_states || {};
    const claimed = new Set((multicolParams || []).map(String).map(value => value.trim()).filter(Boolean));
    const classes = [];
    const titles = [];
    (inputIds || []).forEach(function(inputId) {
      const key = inputId && inputId.file && inputId.name ? `${inputId.file}:${inputId.name}` : '';
      const state = states[key] || {};
      const reserved = inputId && inputId.file === 'tunable' && claimed.has(inputId.name);
      if (reserved && state.state !== 'invalid') {
        classes.push('run-param-container run-param-row--multicol-disabled');
        titles.push('Varied by the Multicol grid above. Remove it there to edit this value directly.');
      } else {
        classes.push(classForState(state, 'run-param-container'));
        titles.push((state.state === 'inactive' || state.state === 'invalid') ? String(state.reason || 'Disabled by the current CLUBB settings.') : '');
      }
    });
    return [classes, titles];
  };

  clientside.runTab.syncAllFlagRowClasses = function(resolution, flagIds) {
    const states = (resolution || {}).control_states || {};
    return (flagIds || []).map(function(flagId) {
      const name = flagId && flagId.name ? flagId.name : '';
      return classForState(states[`flag:${name}`], 'run-param-container run-flag-container');
    });
  };

  clientside.runTab.syncLinkedParamRowClasses = function(resolution, containerIds) {
    const states = (resolution || {}).linked_group_states || {};
    return (containerIds || []).map(function(containerId) {
      const group = containerId && containerId.group ? containerId.group : '';
      return classForState(states[group], 'run-linked-param-container');
    });
  };

  // Apply disabled states synchronously.  A server round trip here can finish
  // after a config switch and reconcile newly rendered inputs with the old
  // component tree, restoring the prior config's values.
  clientside.runTab.syncTunableDisabled = function(multicolParams, resolution, inputIds) {
    const states = (resolution || {}).parameter_states || {};
    const claimed = new Set((multicolParams || []).map(String).map(value => value.trim()).filter(Boolean));
    return (inputIds || []).map(function(inputId) {
      const name = inputId && inputId.name ? String(inputId.name) : '';
      const state = states[name] || {};
      return claimed.has(name) || String(state.state || 'active') !== 'active';
    });
  };

})();
