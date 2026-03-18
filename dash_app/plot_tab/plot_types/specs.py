from dataclasses import dataclass


@dataclass(frozen=True)
class PlotFamilySpec:
    plot_type_id: str
    add_button_id: str
    dropdown_type: str
    time_dependent: bool
    compare_allowed: bool
    default_in_single_mode: bool


PLOT_FAMILY_SPECS = {
    "budget": PlotFamilySpec(
        plot_type_id="budget",
        add_button_id="plots-add-budget",
        dropdown_type="budget-var",
        time_dependent=True,
        compare_allowed=False,
        default_in_single_mode=False,
    ),
    "profile": PlotFamilySpec(
        plot_type_id="profile",
        add_button_id="plots-add-profile",
        dropdown_type="profile-var",
        time_dependent=True,
        compare_allowed=True,
        default_in_single_mode=True,
    ),
    "timeseries": PlotFamilySpec(
        plot_type_id="timeseries",
        add_button_id="plots-add-timeseries",
        dropdown_type="timeseries-var",
        time_dependent=False,
        compare_allowed=True,
        default_in_single_mode=False,
    ),
    "timeheight": PlotFamilySpec(
        plot_type_id="timeheight",
        add_button_id="plots-add-timeheight",
        dropdown_type="timeheight-var",
        time_dependent=False,
        compare_allowed=False,
        default_in_single_mode=False,
    ),
    "subcolumn": PlotFamilySpec(
        plot_type_id="subcolumn",
        add_button_id="plots-add-subcolumn",
        dropdown_type="subcolumn-var",
        time_dependent=True,
        compare_allowed=False,
        default_in_single_mode=False,
    ),
}
