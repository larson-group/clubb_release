/* Draft editing stays in the browser; Python resolves, applies and saves lists. */
window.dash_clientside = Object.assign({}, window.dash_clientside, {
  statsSelection: {
    filterTree: function (query, clear, open) {
      const trigger = window.dash_clientside.callback_context.triggered_id;
      const reset = trigger === 'run-stats-search-clear' || trigger === 'run-stats-custom-open';
      const terms = (reset ? '' : query || '').toLowerCase().trim().split(/\s+/).filter(Boolean);
      const tree = document.getElementById('run-stats-tree');
      if (!tree) return [window.dash_clientside.no_update, ''];
      const nodes = tree.querySelectorAll('.run-stats-node');
      const searching = terms.length > 0;
      const wasSearching = tree.dataset.searching === 'true';
      if (searching && !wasSearching) {
        nodes.forEach(node => { node.dataset.searchOpen = String(node.open); });
      }
      const visiblePaths = new Set();
      const matches = new Set();
      // Hide existing labels instead of rebuilding controls or changing their values.
      tree.querySelectorAll('.run-stats-variable-label').forEach(label => {
        const matched = terms.every(term => label.dataset.search.includes(term));
        const row = label.closest('label');
        if (row.hidden === matched) row.hidden = !matched;
        if (searching && matched) {
          matches.add(label.dataset.name);
          const parts = label.dataset.path.split('/');
          parts.forEach((_, index) => visiblePaths.add(parts.slice(0, index + 1).join('/')));
        }
      });
      nodes.forEach(node => {
        const visible = !searching || visiblePaths.has(node.dataset.path);
        if (node.hidden === visible) node.hidden = !visible;
        if (searching && visible) node.open = true;
        else if (!searching && wasSearching) {
          node.open = node.dataset.searchOpen === 'true';
          delete node.dataset.searchOpen;
        }
      });
      tree.dataset.searching = String(searching);
      const message = !searching ? '' : matches.size ?
        `${matches.size} matching variable${matches.size === 1 ? '' : 's'}` :
        'No variables match. Try a different search.';
      return [reset ? '' : window.dash_clientside.no_update, message];
    },
    editSelection: function (open, cancel, apply, clear, values, variableValues,
                             ids, variableIds, draft, catalog, seed) {
      const noUpdate = window.dash_clientside.no_update;
      const unchanged = () => [noUpdate, noUpdate, ids.map(() => noUpdate),
                               variableIds.map(() => noUpdate), noUpdate];
      const trigger = window.dash_clientside.callback_context.triggered_id;
      if (trigger === 'run-stats-custom-cancel' || trigger === 'run-stats-apply') {
        const result = unchanged();
        result[0] = 'shared-notecard-overlay run-stats-modal--hidden';
        return result;
      }
      if (!catalog || !catalog.all) return unchanged();
      const same = (a, b) => a.length === b.length && a.every((value, index) => value === b[index]);
      let selected = new Set(draft || []);
      if (trigger === 'run-stats-custom-open') {
        selected = new Set((seed || {}).names || []);
      } else if (trigger === 'run-stats-clear') {
        selected.clear();
      } else if (trigger && typeof trigger === 'object') {
        const members = catalog[trigger.path] || [];
        if (trigger.type === 'run-stats-variables') {
          const index = variableIds.findIndex(id => id.path === trigger.path);
          members.forEach(name => selected.delete(name));
          (variableValues[index] || []).forEach(name => selected.add(name));
        } else {
          const index = ids.findIndex(id => id.path === trigger.path);
          const checked = (values[index] || []).length > 0;
          members.forEach(name => checked ? selected.add(name) : selected.delete(name));
        }
      } else {
        return unchanged();
      }
      const ordered = catalog.all.filter(name => selected.has(name));
      if (trigger !== 'run-stats-custom-open' && same(ordered, draft || [])) return unchanged();

      // Skip untouched controls, including their hidden descendants.
      const categoryChecks = ids.map((id, index) => {
        const members = catalog[id.path];
        const next = members.length && members.every(name => selected.has(name)) ? ['selected'] : [];
        return same(next, values[index] || []) ? noUpdate : next;
      });
      const variableChecks = variableIds.map((id, index) => {
        const next = catalog[id.path].filter(name => selected.has(name));
        return same(next, variableValues[index] || []) ? noUpdate : next;
      });
      const div = (children, props) => ({namespace: 'dash_html_components', type: 'Div',
                                        props: Object.assign({children: children}, props)});
      const chips = variableIds.flatMap(id => {
        const members = catalog[id.path];
        const count = members.filter(name => selected.has(name)).length;
        if (!count) return [];
        const name = id.path.split('/').pop().replaceAll('_', ' ');
        const label = name[0].toUpperCase() + name.slice(1) +
                      (count === members.length ? '' : ` · ${count}/${members.length}`);
        return [div(label, {title: id.path, className: 'run-stats-chip'})];
      });
      const preview = [div(`${ordered.length} stats selected`, {className: 'run-stats-preview-count'}),
                       ordered.length ? div(chips, {className: 'run-stats-chips'}) :
                                        div('No statistics will be written.')];
      if (seed && seed.missing && seed.missing.length) {
        preview.push(div('This list contains variables outside the catalog. Applying the tree omits: ' +
                         seed.missing.join(', '), {className: 'run-stats-feedback-warning'}));
      }
      return ['shared-notecard-overlay run-stats-modal', ordered, categoryChecks, variableChecks, preview];
    },
    syncMixedCheckboxes: function (names, catalog) {
      const selected = new Set(names || []);
      // Wait for React to apply checkbox values; mixed states need no tree rerender.
      window.requestAnimationFrame(function () {
        document.querySelectorAll('.run-stats-node > summary .run-stats-check input[type="checkbox"]').forEach(function (input) {
          const path = JSON.parse(input.closest('.run-stats-node').id).path;
          const members = (catalog || {})[path] || [];
          const count = members.filter(name => selected.has(name)).length;
          const mixed = count > 0 && count < members.length;
          input.indeterminate = mixed;
          if (mixed) input.setAttribute('aria-checked', 'mixed');
          else input.removeAttribute('aria-checked');
        });
      });
      return window.dash_clientside.no_update;
    }
  }
});
