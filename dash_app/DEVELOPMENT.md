# Dash App Development Notes

## Design Intent

The dash app should be an interface for doing command-line work quickly. Features should remain implemented and accessible through the CLI first; the dash app should make those commands easier to run, easier to parameterize, and easier to inspect in different views.

Avoid making the dash app the only place a workflow exists. If a workflow matters, make sure there is a script, command, or CLI path that can run it without the app.

## Dropdown Styling

Some Dash dropdowns can render a small input box around the first letter typed into the search field. The compile tab's extra modules dropdown has hit this because it is a searchable `dcc.Dropdown`.

Fix this the same way the other tabs handle dropdowns: scope CSS to the tab/container and clear the native input styling on `.Select-input > input`.

```css
#app-root .some-container .dash-dropdown .Select-input > input {
  caret-color: transparent;
  background: transparent !important;
  border: none !important;
  padding: 0 !important;
  margin: 0 !important;
  box-shadow: none !important;
  outline: none !important;
}

#app-root .some-container .dash-dropdown .Select-input > input:focus {
  outline: none !important;
  box-shadow: none !important;
  border: none !important;
  background: transparent !important;
}
```

Also add the matching single-value/placeholder padding rule when needed so the text aligns with neighboring controls.

## Controls

Do not add extra `+` and `-` buttons for entry boxes unless they are truly useful for the workflow. Numeric inputs and text boxes should stay simple when direct editing is sufficient.
