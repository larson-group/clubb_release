from .budget_plot import PLOT as BUDGET_PLOT
from .profile_plot import PLOT as PROFILE_PLOT
from .subcolumn_plot import PLOT as SUBCOLUMN_PLOT
from .timeheight_plot import PLOT as TIMEHEIGHT_PLOT
from .timeseries_plot import PLOT as TIMESERIES_PLOT

PLOT_TYPES = {
    BUDGET_PLOT.plot_type_id: BUDGET_PLOT,
    PROFILE_PLOT.plot_type_id: PROFILE_PLOT,
    TIMESERIES_PLOT.plot_type_id: TIMESERIES_PLOT,
    TIMEHEIGHT_PLOT.plot_type_id: TIMEHEIGHT_PLOT,
    SUBCOLUMN_PLOT.plot_type_id: SUBCOLUMN_PLOT,
}


def register_plot_callbacks(app):
    for plot in PLOT_TYPES.values():
        plot.register_callbacks(app)
