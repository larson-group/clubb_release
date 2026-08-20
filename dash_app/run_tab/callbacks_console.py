"""Connect Run launch inputs to the browser-owned broker telemetry view."""

from dash import ClientsideFunction, Input, Output


def register_console_callbacks(app):
    """Pass launch selection and action acknowledgements to browser telemetry."""
    app.clientside_callback(
        ClientsideFunction(
            namespace="dashboardRunConsoleState",
            function_name="synchronizeInputs",
        ),
        Output("run-ui-render-signal", "data"),
        Input("run-selected-cases", "data"),
        Input("run-action-result", "data"),
        Input("run-resolved-output-dir", "data"),
    )
