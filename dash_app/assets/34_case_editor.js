/* Simple browser-local cells; Python owns validation and saved cases. */
window.dash_clientside = Object.assign({}, window.dash_clientside, {
  caseEditor: {
    showEditor: function() {
      return window.dash_clientside.callback_context.triggered_id === 'run-case-cancel'
        ? 'shared-notecard-overlay run-case-modal--hidden' : 'shared-notecard-overlay run-case-modal';
    },
    editDraft: function(values, baseline, name, ids, schema, previous) {
      const noUpdate = window.dash_clientside.no_update;
      if (!baseline || !schema || !ids.length || values.length !== ids.length)
        return [noUpdate, noUpdate, noUpdate, noUpdate, noUpdate];
      const trigger = window.dash_clientside.callback_context.triggered_id;
      const get = (object, path) => path.split('.').reduce((value, key) => value && value[key], object);
      const equal = (a, b) => JSON.stringify(a) === JSON.stringify(b);
      const sameSource = previous && previous.source === baseline.key && previous.sha256 === baseline.sha256;
      const initialize = trigger === 'run-case-baseline' || (previous && !sameSource);
      const nextValues = ids.map(() => noUpdate);
      values = values.slice();
      if (initialize) ids.forEach((id, index) => {
        const value = get(baseline.namelists, id.path), spec = schema[id.path];
        const formatted = spec.kind === 'logical' ? (value === undefined ? 1 : value ? 2 : 0)
          : value === undefined ? '' : Array.isArray(value) || value === '' ? JSON.stringify(value) : String(value);
        if (!equal(values[index], formatted)) nextValues[index] = formatted;
        values[index] = formatted;
      });
      if (trigger === 'run-case-name' && name && name !== baseline.key) {
        const index = ids.findIndex(id => id.path === 'stats_setting.fname_prefix');
        if (index >= 0 && values[index] !== name) {
          values[index] = name;
          nextValues[index] = name;
        }
      }
      const draft = JSON.parse(JSON.stringify(baseline.namelists));
      const errors = [], changed = [], groups = {}, logicalStates = {};
      const write = (path, value, present) => {
        const parts = path.split('.');
        let object = draft;
        const parents = [];
        parts.slice(0, -1).forEach(key => {
          if (!object[key]) object[key] = {};
          parents.push([object, key]);
          object = object[key];
        });
        if (present) object[parts.at(-1)] = value;
        else {
          delete object[parts.at(-1)];
          for (let i = parents.length - 1; i >= 0; --i) {
            const [parent, key] = parents[i];
            if (Object.keys(parent[key]).length || (i === 0 && Object.hasOwn(baseline.namelists, key))) break;
            delete parent[key];
          }
        }
      };
      ids.forEach((id, index) => {
        const spec = schema[id.path], before = get(baseline.namelists, id.path);
        let value = values[index];
        let present = value !== null && value !== '';
        try {
          if (spec.kind === 'logical') {
            const state = Number(value);
            if (![0, 1, 2].includes(state)) throw new Error('Choose false, unset or true');
            logicalStates[id.path] = ['false', 'unset', 'true'][state];
            present = state !== 1;
            value = state === 2;
          }
          if (present) {
            if (spec.kind === 'array') {
              value = JSON.parse(value);
              const numeric = item => item === null || (typeof item === 'number' && Number.isFinite(item))
                || (Array.isArray(item) && item.every(numeric));
              if (value === null || !numeric(value)) throw new Error('Use a number or JSON array of numbers');
            } else if (spec.kind === 'number') {
              if (!Number.isFinite(Number(value))) throw new Error('Enter a number');
              value = Number(value);
            } else if (spec.kind === 'text' && value === '""') value = '';
          }
          write(id.path, value, present);
        } catch (error) { errors.push(`${id.path}: ${error.message}`); }
        if (present !== (before !== undefined) || (present && !equal(value, before))) {
          changed.push(id.path);
          const group = id.path.split('.')[0];
          groups[group] = (groups[group] || 0) + 1;
        }
      });
      const sameChanges = sameSource && equal(changed, previous.changed);
      const changedPaths = new Set(changed), beforePaths = new Set(previous ? previous.changed : []);
      const dirtyRows = [...new Set([...changedPaths, ...beforePaths])]
        .filter(path => changedPaths.has(path) !== beforePaths.has(path));
      const sliders = initialize || !previous ? ids.filter(id => schema[id.path].kind === 'logical').map(id => id.path)
        : trigger && schema[trigger.path] && schema[trigger.path].kind === 'logical' ? [trigger.path] : [];
      // No DOM scan or paint for ordinary typing in an already changed cell.
      if (dirtyRows.length || sliders.length) window.requestAnimationFrame(() => {
        dirtyRows.forEach(path => {
          const row = document.getElementById(`run-case-row-${path}`);
          if (row) row.classList.toggle('run-case-field--changed', changedPaths.has(path));
        });
        sliders.forEach(path => {
          const row = document.getElementById(`run-case-row-${path}`);
          const slider = row && row.querySelector('.run-case-bool input');
          if (slider) {
            const state = logicalStates[path];
            if (slider.parentElement.dataset.state !== state) slider.parentElement.dataset.state = state;
            if (slider.getAttribute('aria-valuetext') !== state) slider.setAttribute('aria-valuetext', state);
          }
        });
        new Set(dirtyRows.map(path => path.split('.')[0])).forEach(group => {
          const label = document.getElementById(`run-case-group-changes-${group}`);
          if (label) label.textContent = groups[group] ? `${groups[group]} changed` : '';
        });
      });
      const data = {namelists: draft, errors, name, changed, source: baseline.key, sha256: baseline.sha256};
      const canSave = (draftName, paths, problems) => !!draftName && !problems.length
        && (draftName !== baseline.key || paths.length > 0);
      const saveDisabled = !canSave(name, changed, errors);
      return [sameSource && equal(data, previous) ? noUpdate : data,
        sameChanges ? noUpdate : `${changed.length} ${changed.length === 1 ? 'entry' : 'entries'} changed`,
        sameSource && equal(errors, previous.errors) ? noUpdate : errors.join('\n'),
        sameSource && saveDisabled === !canSave(previous.name, previous.changed, previous.errors) ? noUpdate : saveDisabled,
        nextValues.some(value => value !== noUpdate) ? nextValues : noUpdate];
    }
  }
});
