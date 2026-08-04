(function() {

  const clientside = window.dash_clientside = window.dash_clientside || {};
  clientside.runTab = clientside.runTab || {};

  function classForState(state, baseClass) {
    const classes = [baseClass];
    switch ((state || {}).state) {
      case 'invalid': classes.push('run-param-row--invalid'); break;
      case 'inactive': classes.push('run-param-row--changed-disabled'); break;
      case 'linked-changed': classes.push('run-param-row--changed'); break;
      case 'linked': classes.push('run-param-row--default'); break;
      case 'changed': classes.push('run-param-row--changed'); break;
      default: classes.push('run-param-row--default'); break;
    }
    return classes.join(' ');
  }

  // The shared Python settings evaluator owns every semantic decision.  This
  // browser callback only maps its already-resolved state to CSS classes.
  clientside.runTab.syncAllParamRowClasses = function(resolution, inputIds) {
    const states = (resolution || {}).control_states || {};
    return (inputIds || []).map(function(inputId) {
      const key = inputId && inputId.file && inputId.name ? `${inputId.file}:${inputId.name}` : '';
      return classForState(states[key], 'run-param-container');
    });
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

})();
