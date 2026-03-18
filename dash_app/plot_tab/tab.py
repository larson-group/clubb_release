from dash import dcc

from .callbacks_case import register_case_callbacks
from .callbacks_controls import register_control_callbacks
from .callbacks_grid import register_grid_callbacks
from .callbacks_params import register_param_callbacks
from .layout import build_layout
from .plot_types.registry import register_plot_callbacks
from .state import DEFAULT_OUTPUT_DIR, build_initial_plots_state


def build_tab(app):

    """Build the plots tab and register its callback groups."""

    initial_state = build_initial_plots_state([DEFAULT_OUTPUT_DIR])

    # Wire case and directory selection first because the rest of the tab depends on case_data.
    register_case_callbacks(app)

    # Register time, height, and playback controls that drive the visible plots.
    register_control_callbacks(app)

    # Register parameter and column-selection callbacks used by the right-hand control pane.
    register_param_callbacks(app)

    # Register plot-grid state, add/remove behavior, and dropdown-to-store synchronization.
    register_grid_callbacks(app)

    # Finally register the plot-family figure callbacks that render individual cards.
    register_plot_callbacks(app)

    return dcc.Tab(label="Plots", value="plots", children=build_layout(initial_state))
 
